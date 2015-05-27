#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

using namespace std;

class TopDileptonObjectProducer : public edm::EDFilter
{
public:
  TopDileptonObjectProducer(const edm::ParameterSet& pset);
  ~TopDileptonObjectProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup&) override final;

private:
  bool isGoodVertex(const reco::Vertex& vtx) const;
  bool isGoodMuon(const pat::Muon& m, const reco::Vertex& vtx) const;
  bool isGoodElectron(const pat::Electron& e, const reco::Vertex& vtx) const;
  bool isGoodJet(const pat::Jet& j) const;

  template<typename TColl>
  bool isOverlapToAny(const reco::Candidate& cand, const TColl& coll, const double dR) const;

private:
  typedef std::vector<reco::Vertex> TVertices;
  typedef std::vector<pat::Muon> TMuons;
  typedef std::vector<pat::Electron> TElectrons;
  typedef std::vector<pat::Jet> TJets;
  typedef std::vector<pat::MET> TMETs;

  edm::EDGetTokenT<TVertices> vertexToken_;
  edm::EDGetTokenT<TMuons> muonToken_;
  edm::EDGetTokenT<TElectrons> electronToken_;
  edm::EDGetTokenT<TJets> jetToken_;
  edm::EDGetTokenT<TMETs> metToken_;

  const int cutStepToAccept_;

};

TopDileptonObjectProducer::TopDileptonObjectProducer(const edm::ParameterSet& pset):
  vertexToken_(consumes<TVertices>(pset.getParameter<edm::InputTag>("vertex"))),
  muonToken_(consumes<TMuons>(pset.getParameter<edm::InputTag>("muon"))),
  electronToken_(consumes<TElectrons>(pset.getParameter<edm::InputTag>("electron"))),
  jetToken_(consumes<TJets>(pset.getParameter<edm::InputTag>("jet"))),
  cutStepToAccept_(pset.getParameter<int>("cutStepToAccept"))
{
  produces<int>("mode");
  produces<int>("nVertex");
  produces<int>("step");
//  produces<TMuons>("muons");
//  produces<TElectrons>("electrons");
//  produces<TJets>("jets");
//  produces<TMETs>("mets");
}

bool TopDileptonObjectProducer::filter(edm::Event& event, const edm::EventSetup&)
{
  // Define collections for output
  int mode = 0, passedCutStep = 0;

  // Get objects from event
  edm::Handle<TVertices> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);

  edm::Handle<TMuons> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<TElectrons> electronHandle;
  event.getByToken(electronToken_, electronHandle);

  edm::Handle<TJets> jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<TMETs> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& met = metHandle->at(0);

  // Count good vertices
  int nVertex = 0;
  const reco::Vertex* pv = 0;
  for ( const auto& vtx : *vertexHandle )
  {
    if ( !isGoodVertex(vtx) ) continue;
    if ( !pv ) pv = &vtx;
    ++nVertex;
  }
  if ( !pv ) return false;

  // Objects to be used in plotting
  const reco::Candidate* lepton1 = 0, * lepton2 = 0;
  std::vector<const pat::Jet*> goodJets, bJets;
  do
  {
    // Select leptons - requre full ID and Isolation cuts
    // We are using pair of opposite charged maximizing pt1+pt2, so no need to store 3rd leptons
    for ( const auto& mu : *muonHandle )
    {
      const double pt = mu.pt();
      if ( pt < 20 or std::abs(mu.eta()) > 2.4 ) continue;
      if ( !isGoodMuon(mu, *pv) ) continue;
      //if ( !mu.isIso() ) continue; // FIXME: Add isolation cut here

      if ( mu.charge() > 0 and (!lepton1 or lepton1->pt() < pt) ) lepton1 = &mu;
      else if ( !lepton2 or lepton2->pt() < pt ) lepton2 = &mu;
    }
    for ( const auto& el : *electronHandle )
    {
      const double pt = el.pt();
      if ( pt < 20 or std::abs(el.eta()) > 2.4 ) continue;
      if ( !isGoodElectron(el, *pv) ) continue;
      //if ( el.isIso() ) continue; // FIXME: Add isolation cut here

      if ( el.charge() > 0 and (!lepton1 or lepton1->pt() < pt) ) lepton1 = &el;
      else if ( !lepton2 or lepton2->pt() < pt ) lepton2 = &el;
    }

    // Collect good jets
    for ( const auto& jet : *jetHandle )
    {
      if ( jet.pt() < 30 or std::abs(jet.eta()) > 2.5 ) continue;
      if ( !isGoodJet(jet) ) continue;
      if ( isOverlapToAny(jet, *muonHandle, 0.3) ) continue;
      if ( isOverlapToAny(jet, *electronHandle, 0.3) ) continue;

      goodJets.push_back(&jet);
      //if ( jet.bDiscriminator("CombinedSecondaryVertexTagsV0") > 0.244 ) bJets.push_back(&jet);
    }

    // Now we have all ingredients. Go the cut step determination
    if ( !lepton1 or !lepton2 ) break;
    if      ( lepton1->isMuon() and lepton2->isMuon() ) mode = 1;
    else if ( lepton1->isElectron() and lepton2->isElectron() ) mode = 2;
    else mode = 3;

    // Step 1   Dilepton pair choice
    // After the full ID and Isolation requirements for leptons,
    // one choose the pair of opposite charge maximizing pT_lepton1+pT_lepton2 
    // (or equivalently, the 2 leptons with highest pT)
    // Warning, important remove mass < 20 GeV/c2 at this step
    // This step includes all event-level selections as in the grand table of selections: Electron, Muon, Dilepton, Trigger, Event Cleaning sections.
    const auto& dileptonP4 = lepton1->p4() + lepton2->p4();
    const double dileptonMass = dileptonP4.mass();
    if ( dileptonMass < 20 ) break;

    ++passedCutStep;

    // Step 2   Z mass veto   Dilepton mass not in 76-106 GeV range (Zmass+-15) for ee/mm
    if ( mode != 3 and (dileptonMass >= 76 and dileptonMass <= 106) ) break;

    ++passedCutStep;

    // Step 3   Minimal jet multiplicity  NJets>=2
    if ( goodJets.size() < 2 ) break;

    ++passedCutStep;

    // Step 4   MET cuts  MET > 40 GeV for ee/mm
    if ( mode != 3 and met.pt() <= 40 ) break;

    ++passedCutStep;

    // Step 5   BTagging cuts   DiscriminatorValue > 0.244 for CSV Tagger with loose working point 
    if ( bJets.size() < 2 ) break;

    ++passedCutStep;

  } while ( false );

  event.put(std::auto_ptr<int>(new int(mode)), "mode");
  event.put(std::auto_ptr<int>(new int(nVertex)), "nVertex");
  event.put(std::auto_ptr<int>(new int(passedCutStep)), "step");

  return passedCutStep >= cutStepToAccept_;
}

template<typename TColl>
bool TopDileptonObjectProducer::isOverlapToAny(const reco::Candidate& cand, const TColl& coll, const double dR) const
{
  const double dR2 = dR*dR;
  const double eta = cand.eta(), phi = cand.phi();
  for ( auto x : coll )
  {
    if ( deltaR2(eta, phi, x.eta(), x.phi()) < dR2 ) return true;
  }
  return false;
}

bool TopDileptonObjectProducer::isGoodVertex(const reco::Vertex& vtx) const
{
  return true;
}

bool TopDileptonObjectProducer::isGoodMuon(const pat::Muon& m, const reco::Vertex& vtx) const
{
  return true;
}

bool TopDileptonObjectProducer::isGoodElectron(const pat::Electron& e, const reco::Vertex& vtx) const
{
  return true;
}

bool TopDileptonObjectProducer::isGoodJet(const pat::Jet& j) const
{
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopDileptonObjectProducer);
