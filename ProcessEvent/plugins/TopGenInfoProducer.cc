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

  TH2F* hFH_tpt_tpt_, * hFH_ty_ty_, * hFH_tptCM_tptCM_, * hFH_ttdphi_ttdphi_;
  TH2F* hFH_t1pt_t1pt_, * hFH_t2pt_t2pt_;
  TH2F* hFH_ttpt_ttpt_, * hFH_tty_tty_, * hFH_ttm_ttm_;

  TH2F* hSL_tpt_tpt_, * hSL_ty_ty_, * hSL_tptCM_tptCM_, * hSL_ttdphi_ttdphi_;
  TH2F* hSL_t1pt_t1pt_, * hSL_t2pt_t2pt_;
  TH2F* hSL_ttpt_ttpt_, * hSL_tty_tty_, * hSL_ttm_ttm_;

  TH2F* hFL_tpt_tpt_, * hFL_ty_ty_, * hFL_tptCM_tptCM_, * hFL_ttdphi_ttdphi_;
  TH2F* hFL_t1pt_t1pt_, * hFL_t2pt_t2pt_;
  TH2F* hFL_ttpt_ttpt_, * hFL_tty_tty_, * hFL_ttm_ttm_;

  TH2F* hSL_lpt_lpt_, * hSL_leta_leta_;
  TH2F* hSL_bpt_bpt_, * hSL_beta_beta_, * hSL_bbpt_bbpt_, * hSL_bbm_bbm_;

  TH2F* hFL_lpt_lpt_, * hFL_leta_leta_, * hFL_llpt_llpt_, * hFL_llm_llm_;
  TH2F* hFL_bpt_bpt_, * hFL_beta_beta_, * hFL_bbpt_bbpt_, * hFL_bbm_bbm_;

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

  hFH_tpt_tpt_ = dFH.make<TH2F>("htpt_tpt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_ty_ty_ = dFH.make<TH2F>("hty_ty", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFH_tptCM_tptCM_ = dFH.make<TH2F>("htptCM_tptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_ttdphi_ttdphi_ = dFH.make<TH2F>("httdphi_ttdphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hFH_t1pt_t1pt_ = dFH.make<TH2F>("ht1pt_t1pt", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_t2pt_t2pt_ = dFH.make<TH2F>("ht2pt_t2pt", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hFH_ttpt_ttpt_ = dFH.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFH_tty_tty_ = dFH.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hFH_ttm_ttm_ = dFH.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);

  auto dSL = fs->mkdir("SemiLepton");

  hSL_tpt_tpt_ = dSL.make<TH2F>("htpt_tpt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_ty_ty_ = dSL.make<TH2F>("hty_ty", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hSL_tptCM_tptCM_ = dSL.make<TH2F>("htptCM_tptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_ttdphi_ttdphi_ = dSL.make<TH2F>("httdphi_ttdphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hSL_t1pt_t1pt_ = dSL.make<TH2F>("ht1pt_t1pt", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_t2pt_t2pt_ = dSL.make<TH2F>("ht2pt_t2pt", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hSL_ttpt_ttpt_ = dSL.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_tty_tty_ = dSL.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hSL_ttm_ttm_ = dSL.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);

  auto dFL = fs->mkdir("FullLepton");

  hFL_tpt_tpt_ = dFL.make<TH2F>("htpt_tpt", "pt vs pt;Parton top p_{T} (GeV);Pseudo top p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_ty_ty_ = dFL.make<TH2F>("hty_ty", "y vs y;Parton top y;Pseudo top y", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFL_tptCM_tptCM_ = dFL.make<TH2F>("htptCM_tptCM", "pt at CM vs pt at CM;Parton top p^{*}_{T} (GeV);Pseudo top p^{*}_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_ttdphi_ttdphi_ = dFL.make<TH2F>("httdphi_ttdphi", "#delta#phi vs #delta#phi;Parton top #delta#phi;Pseudo top #delta#phi", 100, 0, TMath::Pi(), 100, 0, TMath::Pi());

  hFL_t1pt_t1pt_ = dFL.make<TH2F>("ht1pt_t1pt", "pt1 vs pt1;Parton top p_{T}^{1st lead} (GeV);Pseudo top p_{T}^{1st lead} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_t2pt_t2pt_ = dFL.make<TH2F>("ht2pt_t2pt", "pt2 vs pt2;Parton top p_{T}^{2nd lead} (GeV);Pseudo top p_{T}^{2nd lead} (GeV)", 500, 0, 500, 500, 0, 500);

  hFL_ttpt_ttpt_ = dFL.make<TH2F>("httpt_ttpt", "ttbar pt vs ttbar pt;Parton t#bar{t} p_{T} (GeV);Pseudo t#bar{t} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_tty_tty_ = dFL.make<TH2F>("htty_tty", "ttbar y vs ttbar y;Parton t#bar{t} y;Pseudo t#bar{t} y", 100, -2.5, 2.5, 100, 2.5, -2.5);
  hFL_ttm_ttm_ = dFL.make<TH2F>("httm_ttm", "ttbar mass vs ttbar mass;Parton t#bar{t} mass (GeV);Pseudo t#bar{t} mass (GeV)", 1600, 0, 1600, 1600, 0, 1600);

  hSL_lpt_lpt_ = dSL.make<TH2F>("hSL_lpt_lpt", "lep pt vs pt;Parton lepton p_{T} (GeV);Pseudo lepton p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_leta_leta_ = dSL.make<TH2F>("hSL_leta_leta", "lep eta vs eta;Parton lepton #eta;Pseudo lepton #eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hSL_bpt_bpt_ = dSL.make<TH2F>("hSL_bpt_bpt", "b pt vs pt;Parton b p_{T} (GeV);Pseudo bjet p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_beta_beta_ = dSL.make<TH2F>("hSL_beta_beta", "lep eta vs eta;Parton b #eta;Pseudo bjet #eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hSL_bbpt_bbpt_ = dSL.make<TH2F>("hSL_bbpt_bbpt", "bb pt vs pt;Parton b#bar{b} p_{T} (GeV);Pseudo b#bar{b} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hSL_bbm_bbm_ = dSL.make<TH2F>("hSL_bbm_bbm", "bb m vs m;Parton b#bar{b} mass (GeV);Pseudo b#bar{b} mass (GeV)", 500, 0, 500, 500, 0, 500);

  hFL_lpt_lpt_ = dFL.make<TH2F>("hFL_lpt_lpt", "lep pt vs pt;Parton lepton p_{T} (GeV);Pseudo lepton p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_leta_leta_ = dFL.make<TH2F>("hFL_leta_leta", "lep eta vs eta;Parton lepton #eta;Pseudo lepton #eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFL_llpt_llpt_ = dFL.make<TH2F>("hFL_llpt_llpt", "ll pt vs pt;Parton ll p_{T} (GeV);Pseudo ll p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_llm_llm_ = dFL.make<TH2F>("hFL_llm_llm", "ll m vs m;Parton ll mass (GeV);Pseudo ll mass (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_bpt_bpt_ = dFL.make<TH2F>("hFL_bpt_bpt", "b pt vs pt;Parton b p_{T} (GeV);Pseudo bjet p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_beta_beta_ = dFL.make<TH2F>("hFL_beta_beta", "lep eta vs eta;Parton b #eta;Pseudo bjet #eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
  hFL_bbpt_bbpt_ = dFL.make<TH2F>("hFL_bbpt_bbpt", "bb pt vs pt;Parton b#bar{b} p_{T} (GeV);Pseudo b#bar{b} p_{T} (GeV)", 500, 0, 500, 500, 0, 500);
  hFL_bbm_bbm_ = dFL.make<TH2F>("hFL_bbm_bbm", "bb m vs m;Parton b#bar{b} mass (GeV);Pseudo b#bar{b} mass (GeV)", 500, 0, 500, 500, 0, 500);
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
      if ( abs(partonQ1->pdgId()) < 6 and (partonQ1->pt() < 30 or std::abs(partonQ1->eta()) > 2.4) ) partonAccept = false;
      if ( partonP2->pt() < 30 or std::abs(partonP2->eta()) > 2.4 ) partonAccept = false;
      if ( abs(partonQ2->pdgId()) < 6 and (partonQ2->pt() < 30 or std::abs(partonQ2->eta()) > 2.4) ) partonAccept = false;
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
  const reco::GenParticle* pseudoB1 = 0, * pseudoB2 = 0;
  //const reco::GenParticle* pseudoW1 = 0, * pseudoW2 = 0;
  const reco::GenParticle* pseudoP1 = 0, * pseudoP2 = 0;
  const reco::GenParticle* pseudoQ1 = 0, * pseudoQ2 = 0;

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
      continue;
    }

    const int motherId = p->mother()->pdgId();
    if      ( motherId ==  24 ) (pseudoP1 ? pseudoQ1 : pseudoP1 ) = p;
    else if ( motherId == -24 ) (pseudoP2 ? pseudoQ2 : pseudoP2 ) = p;
    else if ( motherId ==  6 and abs(pdgId) == 5 ) pseudoB1 = p;
    else if ( motherId == -6 and abs(pdgId) == 5 ) pseudoB2 = p;
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

    const auto partonBB = partonB1->p4()+partonB2->p4();
    const auto pseudoBB = pseudoB1->p4()+pseudoB2->p4();

    if ( partonChannel == 1 )
    {
      hFH_tpt_tpt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFH_tpt_tpt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hFH_ty_ty_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hFH_ty_ty_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hFH_tptCM_tptCM_->Fill(partonPtCM, pseudoPtCM);
      hFH_ttdphi_ttdphi_->Fill(partonDphi, pseudoDphi);

      hFH_t1pt_t1pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFH_t2pt_t2pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hFH_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hFH_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hFH_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());
    }
    else if ( partonChannel == 2 )
    {
      hSL_tpt_tpt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hSL_tpt_tpt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hSL_ty_ty_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hSL_ty_ty_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hSL_tptCM_tptCM_->Fill(partonPtCM, pseudoPtCM);
      hSL_ttdphi_ttdphi_->Fill(partonDphi, pseudoDphi);

      hSL_t1pt_t1pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hSL_t2pt_t2pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hSL_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hSL_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hSL_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());

      if ( abs(partonP1->pdgId()) == 11 or abs(partonP1->pdgId()) == 13 )
      {
        hSL_lpt_lpt_->Fill(partonP1->pt(), pseudoP1->pt());
        hSL_leta_leta_->Fill(partonP1->eta(), pseudoP1->eta());
      }
      else if ( abs(partonP2->pdgId()) == 11 or abs(partonP2->pdgId()) == 13 )
      {
        hSL_lpt_lpt_->Fill(partonP2->pt(), pseudoP2->pt());
        hSL_leta_leta_->Fill(partonP2->eta(), pseudoP2->eta());
      }

      hSL_bpt_bpt_->Fill(partonB1->pt(), pseudoB1->pt());
      hSL_bpt_bpt_->Fill(partonB2->pt(), pseudoB2->pt());
      hSL_beta_beta_->Fill(partonB1->eta(), pseudoB1->eta());
      hSL_beta_beta_->Fill(partonB2->eta(), pseudoB2->eta());
      hSL_bbpt_bbpt_->Fill(partonBB.pt(), pseudoBB.pt());
      hSL_bbm_bbm_->Fill(partonBB.mass(), pseudoBB.mass());
    }
    else if ( partonChannel == 3 )
    {
      hFL_tpt_tpt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFL_tpt_tpt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());
      hFL_ty_ty_->Fill(partonT1Tmp->p4().Rapidity(), pseudoT1->p4().Rapidity());
      hFL_ty_ty_->Fill(partonT2Tmp->p4().Rapidity(), pseudoT2->p4().Rapidity());
      hFL_tptCM_tptCM_->Fill(partonPtCM, pseudoPtCM);
      hFL_ttdphi_ttdphi_->Fill(partonDphi, pseudoDphi);

      hFL_t1pt_t1pt_->Fill(partonT1Tmp->pt(), pseudoT1->pt());
      hFL_t2pt_t2pt_->Fill(partonT2Tmp->pt(), pseudoT2->pt());

      hFL_ttpt_ttpt_->Fill(partonTT.pt(), pseudoTT.pt());
      hFL_tty_tty_->Fill(partonTT.Rapidity(), pseudoTT.Rapidity());
      hFL_ttm_ttm_->Fill(partonTT.mass(), pseudoTT.mass());

      hFL_lpt_lpt_->Fill(partonP1->pt(), pseudoP1->pt());
      hFL_lpt_lpt_->Fill(partonP2->pt(), pseudoP2->pt());
      hFL_lpt_lpt_->Fill(partonP1->eta(), pseudoP1->eta());
      hFL_lpt_lpt_->Fill(partonP2->eta(), pseudoP2->eta());
      const auto partonll = partonP1->p4()+partonP2->p4();
      const auto pseudoll = pseudoP1->p4()+pseudoP2->p4();
      hFL_llpt_llpt_->Fill(partonll.pt(), pseudoll.pt());
      hFL_llm_llm_->Fill(partonll.mass(), pseudoll.mass());

      hFL_bpt_bpt_->Fill(partonB1->pt(), pseudoB1->pt());
      hFL_bpt_bpt_->Fill(partonB2->pt(), pseudoB2->pt());
      hFL_beta_beta_->Fill(partonB1->eta(), pseudoB1->eta());
      hFL_beta_beta_->Fill(partonB2->eta(), pseudoB2->eta());
      hFL_bbpt_bbpt_->Fill(partonBB.pt(), pseudoBB.pt());
      hFL_bbm_bbm_->Fill(partonBB.mass(), pseudoBB.mass());
    }
  }

  event.put(std::auto_ptr<int>(new int(partonChannel)), "partonChannel");
  event.put(std::auto_ptr<int>(new int(partonAccept)), "partonAccept");
  event.put(std::auto_ptr<int>(new int(pseudoChannel)), "pseudoChannel");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopGenInfoProducer);
