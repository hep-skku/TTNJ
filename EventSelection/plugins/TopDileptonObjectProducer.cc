#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
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
  template<typename TColl>
  bool isOverlapToAny(const reco::Candidate& cand, const TColl& coll, const double dR) const;

private:
  typedef std::vector<reco::Vertex> TVertices;
  typedef std::vector<pat::Muon> TMuons;
  typedef std::vector<pat::Electron> TElectrons;
  typedef std::vector<pat::Jet> TJets;
  typedef std::vector<pat::MET> TMETs;
  typedef edm::ValueMap<bool> VMapB;
  typedef edm::ValueMap<float> VMapF;

  template<typename TLepton>
  struct LeptonConsumer
  {
    LeptonConsumer(const edm::ParameterSet& pset, edm::ConsumesCollector && iC)
    {
      srcToken_ = iC.consumes<std::vector<TLepton> >(pset.getParameter<edm::InputTag>("src"));
      idToken_ = iC.consumes<VMapB>(pset.getParameter<edm::InputTag>("id"));
      chIsoToken_ = iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("chIso"));
      nhIsoToken_ = iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("nhIso"));
      phIsoToken_ = iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("phIso"));
      puIsoToken_ = iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("puIso"));
    };

    void load(const edm::Event& event)
    {
      event.getByToken(srcToken_, srcHandle_);
      event.getByToken(idToken_, idHandle_);
      event.getByToken(chIsoToken_, chIsoHandle_);
      event.getByToken(nhIsoToken_, nhIsoHandle_);
      event.getByToken(phIsoToken_, phIsoHandle_);
      event.getByToken(puIsoToken_, puIsoHandle_);
    };

    typedef std::vector<TLepton> TLeptons;

    edm::EDGetTokenT<TLeptons> srcToken_;
    edm::EDGetTokenT<VMapB> idToken_;
    edm::EDGetTokenT<VMapF> chIsoToken_;
    edm::EDGetTokenT<VMapF> nhIsoToken_;
    edm::EDGetTokenT<VMapF> phIsoToken_;
    edm::EDGetTokenT<VMapF> puIsoToken_;

    edm::Handle<TLeptons> srcHandle_;
    edm::Handle<VMapB> idHandle_;
    edm::Handle<VMapF> chIsoHandle_;
    edm::Handle<VMapF> nhIsoHandle_;
    edm::Handle<VMapF> phIsoHandle_;
    edm::Handle<VMapF> puIsoHandle_;
  };

  const edm::EDGetTokenT<TVertices> vertexToken_;
  const edm::EDGetTokenT<TJets> jetToken_;
  const edm::EDGetTokenT<TMETs> metToken_;
  LeptonConsumer<pat::Muon> muonC_;
  LeptonConsumer<pat::Electron> electronC_;

  const int cutStepToAccept_;

};

TopDileptonObjectProducer::TopDileptonObjectProducer(const edm::ParameterSet& pset):
  vertexToken_(consumes<TVertices>(pset.getParameter<edm::InputTag>("vertex"))),
  jetToken_(consumes<TJets>(pset.getParameter<edm::InputTag>("jet"))),
  metToken_(consumes<TMETs>(pset.getParameter<edm::InputTag>("met"))),
  muonC_(pset.getParameter<edm::ParameterSet>("muon"), consumesCollector()),
  electronC_(pset.getParameter<edm::ParameterSet>("electron"), consumesCollector()),
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

  edm::Handle<TJets> jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<TMETs> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& met = metHandle->at(0);

  muonC_.load(event);
  electronC_.load(event);

  // Count good vertices
  const int nVertex = vertexHandle->size();
  if ( nVertex == 0 ) return false;
  //const reco::Vertex& pv = vertexHandle->at(0);

  // Objects to be used in plotting
  const reco::Candidate* lepton1 = 0, * lepton2 = 0;
  std::vector<const pat::Jet*> goodJets, bJets;
  do
  {
    // Select leptons - requre full ID and Isolation cuts
    // We are using pair of opposite charged maximizing pt1+pt2, so no need to store 3rd leptons
    int nGoodLepton = 0;
    for ( size_t i=0, n=muonC_.srcHandle_->size(); i<n; ++i )
    {
      const auto mu = pat::MuonRef(muonC_.srcHandle_, i);
      const double pt = mu->pt();

      const bool id = (*muonC_.idHandle_)[mu];
      const float chIso = (*muonC_.chIsoHandle_)[mu];
      const float nhIso = (*muonC_.nhIsoHandle_)[mu];
      const float phIso = (*muonC_.phIsoHandle_)[mu];
      const float puIso = (*muonC_.puIsoHandle_)[mu];
      const float relIso = (chIso+max(0., nhIso+phIso-0.5*puIso))/pt;

      if ( pt < 20 or std::abs(mu->eta()) > 2.4 ) continue;
      if ( !id || relIso > 0.2 ) continue;
      ++nGoodLepton;

      if ( mu->charge() > 0 and (!lepton1 or lepton1->pt() < pt) ) lepton1 = &*mu;
      else if ( !lepton2 or lepton2->pt() < pt ) lepton2 = &*mu;
    }
    for ( size_t i=0, n=electronC_.srcHandle_->size(); i<n; ++i )
    {
      const auto el = pat::ElectronRef(electronC_.srcHandle_, i);

      const double pt = el->pt();
      const bool id = (*electronC_.idHandle_)[el];
      const float chIso = (*electronC_.chIsoHandle_)[el];
      const float nhIso = (*electronC_.nhIsoHandle_)[el];
      const float phIso = (*electronC_.phIsoHandle_)[el];
      const float puIso = (*electronC_.puIsoHandle_)[el];
      const float relIso = (chIso+max(0., nhIso+phIso-0.5*puIso))/pt;

      if ( pt < 20 or std::abs(el->eta()) > 2.4 ) continue;
      if ( !id || relIso > 0.2 ) continue;
      ++nGoodLepton;

      if ( el->charge() > 0 and (!lepton1 or lepton1->pt() < pt) ) lepton1 = &*el;
      else if ( !lepton2 or lepton2->pt() < pt ) lepton2 = &*el;
    }

    // Collect good jets
    for ( const auto& jet : *jetHandle )
    {
      if ( jet.pt() < 30 or std::abs(jet.eta()) > 2.5 ) continue;
      if ( isOverlapToAny(jet, *muonC_.srcHandle_, 0.3) ) continue;
      if ( isOverlapToAny(jet, *electronC_.srcHandle_, 0.3) ) continue;

      goodJets.push_back(&jet);
      //if ( jet.bDiscriminator("CombinedSecondaryVertexTagsV0") > 0.244 ) bJets.push_back(&jet);
    }

    // Now we have all ingredients. Go the cut step determination
    if ( !lepton1 or !lepton2 ) break;
    //if ( nGoodLepton != 2 ) break; // Require exactly 2 leptons, veto events with 3rd lepton
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

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopDileptonObjectProducer);
