#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include <algorithm>

using namespace std;

class TopGenInfoProducer : public edm::stream::EDProducer<>
{
public:
  TopGenInfoProducer(const edm::ParameterSet& pset);
  ~TopGenInfoProducer() {};

  void produce(edm::Event& event, const edm::EventSetup&) override final;

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> partonToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> pseudoToken_;

};

TopGenInfoProducer::TopGenInfoProducer(const edm::ParameterSet& pset):
  partonToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("parton"))),
  pseudoToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("pseudo")))
{
  produces<int>("partonChannel");
  produces<int>("partonAccept");
  produces<int>("pseudoChannel");
}

void TopGenInfoProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<reco::GenParticleCollection> partonHandle;
  event.getByToken(partonToken_, partonHandle);

  edm::Handle<reco::GenParticleCollection> pseudoHandle;
  event.getByToken(pseudoToken_, pseudoHandle);

  // parton T->BW, W->PQ
  //const reco::GenParticle* partonT1 = 0, * partonT2 = 0;
  const reco::GenParticle* partonB1 = 0, * partonB2 = 0;
  //const reco::GenParticle* partonW1 = 0, * partonW2 = 0;
  const reco::GenParticle* partonP1 = 0, * partonP2 = 0;
  const reco::GenParticle* partonQ1 = 0, * partonQ2 = 0;

  int partonNLepton = 0;
  int partonNtau = 0;
  for ( auto& x : *partonHandle )
  {
    const reco::GenParticle* p = &x;
    const int pdgId = p->pdgId();
    if ( abs(pdgId) == 6 ) continue;

    if ( abs(p->mother()->pdgId()) == 6 )
    {
      if      ( pdgId > 0 and pdgId <  6 ) partonB1 = p;
      else if ( pdgId < 0 and pdgId > -6 ) partonB2 = p;
    }
    else
    {
      if ( abs(pdgId) == 11 or abs(pdgId) == 13 ) ++partonNLepton;
      else if ( abs(pdgId) == 15 ) ++partonNtau;

      int motherId = p->mother()->pdgId();
      if ( abs(motherId) == 15 ) motherId = p->mother()->mother()->pdgId();

      if      ( motherId ==  24 ) (partonP1 ? partonQ1 : partonP1 ) = p;
      else if ( motherId == -24 ) (partonP2 ? partonQ2 : partonP2 ) = p;
    }
  }
  bool partonAccept = true;
  int partonChannel = 0; // CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON
  if ( !partonB1 or !partonB2 or !partonP1 or !partonQ1 or !partonP2 or !partonQ2 )
  {
    partonAccept = false;
  }
  else
  {
    partonChannel = 1 + partonNLepton;
    if ( partonNtau != 0 ) partonChannel += 3; // Shift ch. if tau is in the decay chain
    
    if ( abs(partonP1->pdgId()) > abs(partonQ1->pdgId()) ) std::swap(partonP1, partonQ1);
    if ( abs(partonP2->pdgId()) > abs(partonQ2->pdgId()) ) std::swap(partonP2, partonQ2);

    if ( partonB1->pt() < 30 or std::abs(partonB1->eta()) > 2.4 ) partonAccept = false;
    if ( partonB2->pt() < 30 or std::abs(partonB2->eta()) > 2.4 ) partonAccept = false;

    if ( partonChannel == 1 )
    {
      if ( partonP1->pt() < 30 or std::abs(partonP1->eta()) > 2.4 ) partonAccept = false;
      if ( partonQ1->pt() < 30 or std::abs(partonQ1->eta()) > 2.4 ) partonAccept = false;
      if ( partonP2->pt() < 30 or std::abs(partonP2->eta()) > 2.4 ) partonAccept = false;
      if ( partonQ2->pt() < 30 or std::abs(partonQ2->eta()) > 2.4 ) partonAccept = false;
    }
    else if ( partonChannel == 2 )
    {
      if ( partonP1->pt() < 30 or std::abs(partonP1->eta()) > 2.4 ) partonAccept = false;
      if ( partonQ1->pt() < 30 or std::abs(partonQ1->eta()) > 2.4 ) partonAccept = false;
      if ( partonP2->pt() < 30 or std::abs(partonP2->eta()) > 2.4 ) partonAccept = false;
      if ( partonQ2->pt() < 30 or std::abs(partonQ2->eta()) > 2.4 ) partonAccept = false;
    }
    else if ( partonChannel == 3 )
    {
      if ( partonP1->pt() < 20 or std::abs(partonP1->eta()) > 2.4 ) partonAccept = false;
      if ( partonP2->pt() < 20 or std::abs(partonP2->eta()) > 2.4 ) partonAccept = false;
    }
  }

  // pseudo T->BW, W->PQ
  //const reco::GenParticle* pseudoT1 = 0, * pseudoT2 = 0;
  //const reco::GenParticle* pseudoB1 = 0, * pseudoB2 = 0;
  //const reco::GenParticle* pseudoW1 = 0, * pseudoW2 = 0;
  //const reco::GenParticle* pseudoP1 = 0, * pseudoP2 = 0;
  //const reco::GenParticle* pseudoQ1 = 0, * pseudoQ2 = 0;

  int pseudoNLepton = 0;
  for ( auto& x : *pseudoHandle )
  {
    const reco::GenParticle* p = &x;
    const int pdgId = p->pdgId();
    if ( abs(pdgId) == 11 or abs(pdgId) == 13 ) ++pseudoNLepton;

    //const int motherId = p->mother()->pdgId();
    //if      ( motherId ==  24 ) (pseudoP1 ? pseudoQ1 : pseudoP1 ) = p;
    //else if ( motherId == -24 ) (pseudoP2 ? pseudoQ2 : pseudoP2 ) = p;
  }
  const int pseudoChannel = pseudoNLepton+1; // CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON

  event.put(std::auto_ptr<int>(new int(partonChannel)), "partonChannel");
  event.put(std::auto_ptr<int>(new int(partonAccept)), "partonAccept");
  event.put(std::auto_ptr<int>(new int(pseudoChannel)), "pseudoChannel");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopGenInfoProducer);
