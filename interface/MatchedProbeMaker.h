#ifndef MuonAnalysis_TagAndProbe_MatchedProbeMaker_H
#define MuonAnalysis_TagAndProbe_MatchedProbeMaker_H

// system include files
#include <memory>
#include <vector>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//
// class decleration
//

template< typename T >
class MatchedProbeMaker : public edm::EDProducer 
{
   public:
      typedef std::vector< T > collection;

      explicit MatchedProbeMaker(const edm::ParameterSet& iConfig);

      ~MatchedProbeMaker();
      
   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      edm::InputTag m_candidateSource;
      edm::InputTag m_referenceSource;
      edm::InputTag m_resMatchMapSource;

      bool matched_;    

};

template< typename T >
MatchedProbeMaker<T>::MatchedProbeMaker(const edm::ParameterSet& iConfig) :
   m_candidateSource(iConfig.getUntrackedParameter<edm::InputTag>("CandidateSource")),
   m_referenceSource(iConfig.getUntrackedParameter<edm::InputTag>("ReferenceSource")),
   m_resMatchMapSource(iConfig.getUntrackedParameter<edm::InputTag>("ResMatchMapSource")),
   matched_(iConfig.getUntrackedParameter< bool >("Matched",true))
{
   //register your products
   produces< edm::RefVector< collection > >();
}


template< typename T >
MatchedProbeMaker<T>::~MatchedProbeMaker(){}


template< typename T >
void MatchedProbeMaker<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr< edm::RefVector< collection > > outputCollection_matched( new edm::RefVector< collection > );
  std::auto_ptr< edm::RefVector< collection > > outputCollection_unmatched(new edm::RefVector< collection > );
  
  // Get the candidates from the event
  edm::Handle< edm::RefVector< collection > > Cands;
  iEvent.getByLabel(m_candidateSource,Cands);

  // Get the resolution matching map from the event
  edm::Handle<reco::CandViewMatchMap> ResMatchMap;
  iEvent.getByLabel(m_resMatchMapSource,ResMatchMap);
  
  //std::cout<< "Cands size " << Cands->size() << std::endl;
  //std::cout<< "Map size " << ResMatchMap->size() << std::endl;

  // Loop over the candidates looking for a match
  for (unsigned i=0; i<Cands->size(); i++)
  {
     const edm::Ref< collection > CandRef = (*Cands)[i];

     reco::CandidateBaseRef candBaseRef( CandRef );

     // Loop over match map
     reco::CandViewMatchMap::const_iterator f = ResMatchMap->find( candBaseRef );
     if( f!=ResMatchMap->end() )
     {
	//std::cout << "Match found" << std::endl;
	outputCollection_matched->push_back(CandRef);      
     }
     else 
     {
	//std::cout<< "Match not found" << std::endl;
	outputCollection_unmatched->push_back(CandRef);      
     }
  }  

  //iEvent.put(outputCollection_matched,"matched");
  //iEvent.put(outputCollection_unmatched,"unmatched");

  if( matched_ ) iEvent.put( outputCollection_matched );
  else           iEvent.put( outputCollection_unmatched );

}

// ------------ method called once each job just before starting event loop  ------------
template< typename T >
void MatchedProbeMaker<T>::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T >
void MatchedProbeMaker<T>::endJob() 
{
}

#endif


