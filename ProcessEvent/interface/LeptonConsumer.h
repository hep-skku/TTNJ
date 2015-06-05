#ifndef TTNJ_ProcessEvent_LeptonConsumer_H
#define TTNJ_ProcessEvent_LeptonConsumer_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <vector>
#include <memory>

template<typename TLepton>
class LeptonConsumer
{
public:
  typedef std::auto_ptr<edm::PtrVector<reco::Candidate> > ACandPtr;
  typedef edm::ValueMap<bool> VMapB;
  typedef edm::ValueMap<float> VMapF;

  LeptonConsumer(const edm::ParameterSet& pset, edm::ConsumesCollector && iC):
    srcToken_(iC.consumes<std::vector<TLepton> >(pset.getParameter<edm::InputTag>("src"))),
    idToken_(iC.consumes<VMapB>(pset.getParameter<edm::InputTag>("id"))),
    chIsoToken_(iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("chIso"))),
    nhIsoToken_(iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("nhIso"))),
    phIsoToken_(iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("phIso"))),
    puIsoToken_(iC.consumes<VMapF>(pset.getParameter<edm::InputTag>("puIso"))),
    minPt_(pset.getUntrackedParameter<double>("minPt", 20)),
    maxEta_(pset.getUntrackedParameter<double>("maxEta", 2.4)),
    maxRiso_(pset.getUntrackedParameter<double>("maxRiso", 0.2)) {};

  ACandPtr load(const edm::Event& event) const
  {
    ACandPtr out(new edm::PtrVector<reco::Candidate>);

    edm::Handle<std::vector<TLepton> > srcHandle;
    edm::Handle<VMapB> idHandle;
    edm::Handle<VMapF> chIsoHandle, nhIsoHandle, phIsoHandle, puIsoHandle;

    event.getByToken(srcToken_, srcHandle);
    event.getByToken(idToken_, idHandle);
    event.getByToken(chIsoToken_, chIsoHandle);
    event.getByToken(nhIsoToken_, nhIsoHandle);
    event.getByToken(phIsoToken_, phIsoHandle);
    event.getByToken(puIsoToken_, puIsoHandle);

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      const edm::Ref<std::vector<TLepton> > r(srcHandle, i);
      const double pt = r->pt();

      const bool id = (*idHandle)[r];
      const float chIso = (*chIsoHandle)[r];
      const float nhIso = (*nhIsoHandle)[r];
      const float phIso = (*phIsoHandle)[r];
      const float puIso = (*puIsoHandle)[r];
      const float sumIso = chIso+std::max(0., nhIso+phIso-0.5*puIso);

      if ( pt < minPt_ or std::abs(r->eta()) > maxEta_ ) continue;
      if ( !id || sumIso > maxRiso_*pt ) continue;

      out->push_back(edm::Ptr<TLepton>(srcHandle, i));
    }

    return out;
  }

private:
  const edm::EDGetTokenT<std::vector<TLepton> > srcToken_;
  const edm::EDGetTokenT<VMapB> idToken_;
  const edm::EDGetTokenT<VMapF> chIsoToken_;
  const edm::EDGetTokenT<VMapF> nhIsoToken_;
  const edm::EDGetTokenT<VMapF> phIsoToken_;
  const edm::EDGetTokenT<VMapF> puIsoToken_;

  const double minPt_, maxEta_, maxRiso_;

};

#endif

