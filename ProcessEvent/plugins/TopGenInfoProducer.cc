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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TH2F.h>
#include <TMath.h>
#include <Math/Boost.h>
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

  TH2F* hFH_pt_pt_, * hFH_y_y_, * hFH_ptCM_ptCM_, * hFH_dphi_dphi_;
  TH2F* hFH_pt1_pt1_, * hFH_pt2_pt2_;
  TH2F* hFH_ttpt_ttpt_, * hFH_tty_tty_, * hFH_ttm_ttm_;

  TH2F* hSL_pt_pt_, * hSL_y_y_, * hSL_ptCM_ptCM_, * hSL_dphi_dphi_;
  TH2F* hSL_pt1_pt1_, * hSL_pt2_pt2_;
  TH2F* hSL_ttpt_ttpt_, * hSL_tty_tty_, * hSL_ttm_ttm_;

  TH2F* hFL_pt_pt_, * hFL_y_y_, * hFL_ptCM_ptCM_, * hFL_dphi_dphi_;
  TH2F* hFL_pt1_pt1_, * hFL_pt2_pt2_;
  TH2F* hFL_ttpt_ttpt_, * hFL_tty_tty_, * hFL_ttm_ttm_;
};

TopGenInfoProducer::TopGenInfoProducer(const edm::ParameterSet& pset):
  partonToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("parton"))),
  pseudoToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("pseudo")))
{
  produces<int>("partonChannel");
  produces<int>("partonAccept");
  produces<int>("pseudoChannel");

  edm::Service<TFileService> fs;
  auto dFH = fs->mkdir("FullHadron");

  hFH_pt_pt_ = dFH.make<TH2F>("hpt_pt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_y_y_ = dFH.make<TH2F>("hy_y", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFH_ptCM_ptCM_ = dFH.make<TH2F>("hptCM_ptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_dphi_dphi_ = dFH.make<TH2F>("hdphi_dphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hFH_pt1_pt1_ = dFH.make<TH2F>("hpt1_pt1", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_pt2_pt2_ = dFH.make<TH2F>("hpt2_pt2", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hFH_ttpt_ttpt_ = dFH.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_tty_tty_ = dFH.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hFH_ttm_ttm_ = dFH.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);

  auto dSL = fs->mkdir("SemiLepton");

  hSL_pt_pt_ = dSL.make<TH2F>("hpt_pt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_y_y_ = dSL.make<TH2F>("hy_y", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hSL_ptCM_ptCM_ = dSL.make<TH2F>("hptCM_ptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_dphi_dphi_ = dSL.make<TH2F>("hdphi_dphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hSL_pt1_pt1_ = dSL.make<TH2F>("hpt1_pt1", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_pt2_pt2_ = dSL.make<TH2F>("hpt2_pt2", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hSL_ttpt_ttpt_ = dSL.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_tty_tty_ = dSL.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hSL_ttm_ttm_ = dSL.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);

  auto dFL = fs->mkdir("FullLepton");

  hFL_pt_pt_ = dFL.make<TH2F>("hpt_pt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_y_y_ = dFL.make<TH2F>("hy_y", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFL_ptCM_ptCM_ = dFL.make<TH2F>("hptCM_ptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_dphi_dphi_ = dFL.make<TH2F>("hdphi_dphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hFL_pt1_pt1_ = dFL.make<TH2F>("hpt1_pt1", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_pt2_pt2_ = dFL.make<TH2F>("hpt2_pt2", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hFL_ttpt_ttpt_ = dFL.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_tty_tty_ = dFL.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hFL_ttm_ttm_ = dFL.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);
}

void TopGenInfoProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<reco::GenParticleCollection> partonHandle;
  event.getByToken(partonToken_, partonHandle);

  edm::Handle<reco::GenParticleCollection> pseudoHandle;
  event.getByToken(pseudoToken_, pseudoHandle);

  // parton T->BW, W->PQ
  const reco::GenParticle* partonT1 = 0, * partonT2 = 0;
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
    if      ( pdgId ==  6 ) { partonT1 = p; continue; }
    else if ( pdgId == -6 ) { partonT2 = p; continue; }

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
  const reco::GenParticle* pseudoT1 = 0, * pseudoT2 = 0;
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
    else if ( abs(pdgId) == 6 )
    {
      if ( !pseudoT1 ) pseudoT1 = p;
      else if ( !pseudoT2 ) pseudoT2 = p;
    }

    //const int motherId = p->mother()->pdgId();
    //if      ( motherId ==  24 ) (pseudoP1 ? pseudoQ1 : pseudoP1 ) = p;
    //else if ( motherId == -24 ) (pseudoP2 ? pseudoQ2 : pseudoP2 ) = p;
  }
  const int pseudoChannel = pseudoNLepton+1; // CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON
  if ( pseudoT1 and pseudoT2 and pseudoT1->pt() < pseudoT2->pt() ) std::swap(pseudoT1, pseudoT2);

  if ( partonAccept and pseudoT1 and pseudoT2 )
  {
    const reco::GenParticle* partonT1Tmp = partonT1, * partonT2Tmp = partonT2;
    if ( partonT1Tmp->pt() < partonT2Tmp->pt() ) std::swap(partonT1Tmp, partonT2Tmp);
    const auto partonTT = partonT1->p4()+partonT2->p4();
    const auto pseudoTT = pseudoT1->p4()+pseudoT2->p4();
    const double partonPtCM = (ROOT::Math::Boost(partonTT.BoostToCM())(partonT1->p4())).pt();
    const double pseudoPtCM = (ROOT::Math::Boost(pseudoTT.BoostToCM())(pseudoT1->p4())).pt();
    const double partonDphi = deltaPhi(partonT1->phi(), partonT2->phi());
    const double pseudoDphi = deltaPhi(pseudoT1->phi(), pseudoT2->phi());

    if ( partonChannel == 1 )
    {
      hFH_pt_pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFH_pt_pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hFH_y_y_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hFH_y_y_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hFH_ptCM_ptCM_->Fill(partonPtCM, pseudoPtCM);
      hFH_dphi_dphi_->Fill(partonDphi, pseudoDphi);

      hFH_pt1_pt1_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFH_pt2_pt2_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hFH_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hFH_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hFH_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());
    }
    else if ( partonChannel == 2 )
    {
      hSL_pt_pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hSL_pt_pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hSL_y_y_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hSL_y_y_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hSL_ptCM_ptCM_->Fill(partonPtCM, pseudoPtCM);
      hSL_dphi_dphi_->Fill(partonDphi, pseudoDphi);

      hSL_pt1_pt1_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hSL_pt2_pt2_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hSL_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hSL_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hSL_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());
    }
    else if ( partonChannel == 3 )
    {
      hFL_pt_pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFL_pt_pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hFL_y_y_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hFL_y_y_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hFL_ptCM_ptCM_->Fill(partonPtCM, pseudoPtCM);
      hFL_dphi_dphi_->Fill(partonDphi, pseudoDphi);

      hFL_pt1_pt1_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFL_pt2_pt2_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hFL_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hFL_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hFL_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());
    }
  }

  event.put(std::auto_ptr<int>(new int(partonChannel)), "partonChannel");
  event.put(std::auto_ptr<int>(new int(partonAccept)), "partonAccept");
  event.put(std::auto_ptr<int>(new int(pseudoChannel)), "pseudoChannel");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopGenInfoProducer);
