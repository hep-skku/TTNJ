#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/ConsumesCollector.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/Common/interface/Handle.h"

#include "TTNJ/EventSelection/interface/LeptonConsumer.h"
#include "TTNJ/EventSelection/interface/JetConsumer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include <algorithm>

using namespace std;

class TopDileptonProducer : public edm::EDFilter
{
public:
  TopDileptonProducer(const edm::ParameterSet& pset);
  ~TopDileptonProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup&) override final;

private:
  template<typename TColl>
  bool isOverlapToAny(const reco::Candidate& cand, const TColl& coll, const double dR) const;

private:
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<pat::METCollection> metToken_;
  const LeptonConsumer<pat::Muon> muonC_;
  const LeptonConsumer<pat::Electron> electronC_;
  const JetConsumer jetC_;

  const int cutStepToAccept_;

};

TopDileptonProducer::TopDileptonProducer(const edm::ParameterSet& pset):
  vertexToken_(consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"))),
  metToken_(consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("met"))),
  muonC_(pset.getParameter<edm::ParameterSet>("muon"), consumesCollector()),
  electronC_(pset.getParameter<edm::ParameterSet>("electron"), consumesCollector()),
  jetC_(pset.getParameter<edm::ParameterSet>("jet"), consumesCollector()),
  cutStepToAccept_(pset.getParameter<int>("cutStepToAccept"))
{
  produces<int>("mode");
  produces<int>("nVertex");
  produces<int>("step");

  produces<std::vector<reco::LeafCandidate> >("leptons");
  produces<pat::JetCollection>("jets");
  produces<edm::RefVector<pat::JetCollection> >("bjets");
//  produces<pat::METCollection>("mets");
}

bool TopDileptonProducer::filter(edm::Event& event, const edm::EventSetup&)
{
  // Define collections for output
  int mode = 0, passedCutStep = 0;

  const auto muons = muonC_.load(event);
  const auto electrons = electronC_.load(event);
  const auto jets = jetC_.load(event);

  // Count good vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  const int nVertex = vertexHandle->size();
  if ( nVertex == 0 ) return false;
  //const reco::Vertex& pv = vertexHandle->at(0);

  edm::Handle<pat::METCollection> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& met = metHandle->at(0);

  typedef std::vector<reco::LeafCandidate> LeafCands;
  std::auto_ptr<LeafCands> out_leptons(new LeafCands);
  std::auto_ptr<pat::JetCollection> out_jets(new pat::JetCollection);
  std::vector<int> bjetIndex;

  return true;

  do
  {
    // Step 1   Dilepton pair choice
    // After the full ID and Isolation requirements for leptons,
    // one choose the pair of opposite charge maximizing pT_lepton1+pT_lepton2 
    // (or equivalently, the 2 leptons with highest pT)
    // Warning, important remove mass < 20 GeV/c2 at this step
    // This step includes all event-level selections as in the grand table of selections: Electron, Muon, Dilepton, Trigger, Event Cleaning sections.
    for ( auto muPtr : *muons ) out_leptons->push_back(reco::LeafCandidate(*muPtr));
    for ( auto elPtr : *electrons ) out_leptons->push_back(reco::LeafCandidate(*elPtr));
    std::sort(out_leptons->begin(), out_leptons->end(), GreaterByPt<reco::Candidate>());

    // At least two leptons in acceptance and id/isolation cut
    if ( out_leptons->size() < 2 ) break;
    const auto& lepton1 = out_leptons->at(0);
    const auto& lepton2 = out_leptons->at(1);

    // Split channels
    if ( lepton1.isMuon() and lepton2.isMuon() ) mode = 1;
    else if ( lepton1.isElectron() and lepton2.isElectron() ) mode = 2;
    else mode = 3;
    
    // Choose event with opposite signed pair
    if ( lepton1.charge() == lepton2.charge() )
    {
      mode *= -1; // Flip mode number in case of SS case
      break;
    }

    // Minimum dilepton mass
    const auto dileptonP4 = lepton1.p4() + lepton2.p4();
    const double dileptonMass = dileptonP4.mass();
    if ( dileptonMass < 20 ) break;

    ++passedCutStep; // Now this passes cut step 1.

    // Step 2   Z mass veto   Dilepton mass not in 76-106 GeV range (Zmass+-15) for ee/mm
    if ( mode != 3 and (dileptonMass >= 76 and dileptonMass <= 106) ) break;

    ++passedCutStep;

    for ( auto& jet : *jets )
    {
      if ( isOverlapToAny(*jet, *muons, 0.3) ) continue;
      if ( isOverlapToAny(*jet, *electrons, 0.3) ) continue;
      out_jets->push_back(*jet);
      if ( jet->bDiscriminator("CombinedSecondaryVertexTagsV0") > 0.244 ) bjetIndex.push_back(out_jets->size()-1);
    }

    // Step 3   Minimal jet multiplicity  NJets>=2
    if ( out_jets->size() < 2 ) break;

    ++passedCutStep;

    // Step 4   MET cuts  MET > 40 GeV for ee/mm
    if ( mode != 3 and met.pt() <= 40 ) break;

    ++passedCutStep;

    // Step 5   BTagging cuts   DiscriminatorValue > 0.244 for CSV Tagger with loose working point 
    if ( bjetIndex.size() < 2 ) break;

    ++passedCutStep;

  } while ( false );

  event.put(std::auto_ptr<int>(new int(mode)), "mode");
  event.put(std::auto_ptr<int>(new int(nVertex)), "nVertex");
  event.put(std::auto_ptr<int>(new int(passedCutStep)), "step");

  event.put(out_leptons, "leptons");
  auto out_jetHandle = event.put(out_jets, "jets");
  std::auto_ptr<edm::RefVector<pat::JetCollection> > out_bjets(new edm::RefVector<pat::JetCollection>);
  for ( auto i : bjetIndex ) out_bjets->push_back(pat::JetRef(out_jetHandle, i));
  event.put(out_jets, "bjets");

  return passedCutStep >= cutStepToAccept_;
}

template<typename TColl>
bool TopDileptonProducer::isOverlapToAny(const reco::Candidate& cand, const TColl& coll, const double dR) const
{
  const double dR2 = dR*dR;
  const double eta = cand.eta(), phi = cand.phi();
  for ( auto x : coll )
  {
    if ( deltaR2(eta, phi, x->eta(), x->phi()) < dR2 ) return true;
  }
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopDileptonProducer);
