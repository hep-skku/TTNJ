#ifndef TTNJ_ProcessEvent_JetConsumer_H
#define TTNJ_ProcessEvent_JetConsumer_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>
#include <memory>

class JetConsumer
{
public:
  typedef std::auto_ptr<edm::PtrVector<pat::Jet> > ACandPtr;
  typedef edm::ValueMap<bool> VMapB;
  typedef edm::ValueMap<float> VMapF;

  JetConsumer(const edm::ParameterSet& pset, edm::ConsumesCollector && iC):
    srcToken_(iC.consumes<std::vector<pat::Jet> >(pset.getParameter<edm::InputTag>("src"))),
    minPt_(pset.getUntrackedParameter<double>("minPt", 20)),
    maxEta_(pset.getUntrackedParameter<double>("maxEta", 2.4)) {};

  ACandPtr load(const edm::Event& event) const
  {
    ACandPtr out(new edm::PtrVector<pat::Jet>);

    edm::Handle<std::vector<pat::Jet> > srcHandle;
    event.getByToken(srcToken_, srcHandle);

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      edm::Ref<std::vector<pat::Jet> > r(srcHandle, i);
      const double pt = r->pt();

      if ( pt < minPt_ or std::abs(r->eta()) > maxEta_ ) continue;

      out->push_back(edm::Ptr<pat::Jet>(srcHandle, i));
    }

    return out;
  }

private:
  const edm::EDGetTokenT<std::vector<pat::Jet> > srcToken_;
  const double minPt_, maxEta_;

};

#endif

