// -*- C++ -*-
//
// Package:    OniaPhotonKinematicFit
// Class:      OniaPhotonKinematicFit
// 
/**

 Description: Kinematic Fit for Onia + Photon

 Implementation:
     Original work from Stefano Argiro, and the Torino group
     Adapter for MINIAOD by Alberto Sanchez-Hernandez
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h> 
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class OniaPhotonKinematicFit : public edm::EDProducer {
  public:
    explicit OniaPhotonKinematicFit(const edm::ParameterSet&);
    ~OniaPhotonKinematicFit() override {};
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;
	void endJob() override;

	double meson_nS_mass_;
	bool isDebug_;
	std::string product_name_;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label;
	
	int candidates;
	
	template<typename T>
	struct GreaterByVProb {
	  typedef T first_argument_type;
	  typedef T second_argument_type;
	  bool operator()( const T & t1, const T & t2 ) const {
			return t1.userFloat("vProb") > t2.userFloat("vProb");
	  }
	};
};

OniaPhotonKinematicFit::OniaPhotonKinematicFit(const edm::ParameterSet& iConfig) {
  chi_Label     = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("chi_cand"));
  meson_nS_mass_ = iConfig.getParameter<double>("meson_nS_mass");
  product_name_ = iConfig.getParameter<std::string>("product_name");
  produces<pat::CompositeCandidateCollection>(product_name_);
  isDebug_ = iConfig.getParameter<bool>("is_Debug");
  candidates = 0;
}

// ------------ method called to produce the data  ------------
void OniaPhotonKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // Grab parameters
  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByToken(chi_Label, chiCandHandle);
  
  //Kinematic refit collection
  std::unique_ptr<pat::CompositeCandidateCollection> chicCompCandRefitColl(new pat::CompositeCandidateCollection);
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> TrkBuilder; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TrkBuilder); 
  
  int indexConversion=-1;
  
  //Empty dummy pat
  reco::CompositeCandidate recoempty(0, math::XYZTLorentzVector(0,0,0,0), math::XYZPoint(0,0,0), -99);
  pat::CompositeCandidate patempty(recoempty);
  
  //Internal counter
  int validchi = 0;

  for (pat::CompositeCandidateCollection::const_iterator chiCand=chiCandHandle->begin();chiCand!=chiCandHandle->end();++chiCand) { 
    indexConversion++;
	if(isDebug_) std::cout<<">>>>> Looping over chicand <<<<<"<<std::endl;
    reco::TrackRef JpsiTk[2]={
      ( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon2") ) )->innerTrack()
    };
      
    reco::TrackCollection convTracks;
    const reco::Track tk0=*(dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("photon"))->userData<reco::Track>("track0"));
    const reco::Track tk1=*(dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("photon"))->userData<reco::Track>("track1"));
    convTracks.push_back(tk0);
    convTracks.push_back(tk1);

    std::vector<reco::TransientTrack> MuMuTT;
    MuMuTT.push_back((*TrkBuilder).build(&JpsiTk[0]));
    MuMuTT.push_back((*TrkBuilder).build(&JpsiTk[1]));
      
    std::vector<reco::TransientTrack> EETT;
    EETT.push_back((*TrkBuilder).build(convTracks[0]));
    EETT.push_back((*TrkBuilder).build(convTracks[1]));
      
    const ParticleMass zero_mass(0);
    float zero_sigma = 1E-6;
      
    const ParticleMass eleMass(0.000511);
    float eleSigma = 1E-6;
      
    KinematicParticleFactoryFromTransientTrack pFactory;
    std::vector<RefCountedKinematicParticle> PhotonParticles;
    PhotonParticles.push_back(pFactory.particle(EETT[0],eleMass,float(0),float(0),eleSigma));
    PhotonParticles.push_back(pFactory.particle(EETT[1],eleMass,float(0),float(0),eleSigma));
    
    KinematicParticleVertexFitter fitter;
    RefCountedKinematicTree photonVertexFitTree;
    photonVertexFitTree = fitter.fit(PhotonParticles);
    
	//try to fit again
    if (!photonVertexFitTree->isValid()) { 
      edm::ParameterSet pSet;
      pSet.addParameter<double>("maxDistance", 3);
      pSet.addParameter<int>("maxNbrOfIterations", 10000);
      KinematicParticleVertexFitter fitter2(pSet);
      photonVertexFitTree = fitter2.fit(PhotonParticles);
    }  
	//If the process is stopped because kinematic fit / cut failed, push_back an empty CompositeCandidate into chicCompCandRefitColl, thus, each chiCand would have its own refit1P (even though it might be empty) to access in the tuplizer
    if (!photonVertexFitTree->isValid()) {
		if(isDebug_) std::cout<<"Photon vertex fit is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}

	// now apply Photon mass constraint
	KinematicParticleFitter csFitterPhoton;
	KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);

	// add mass constraint to the photon fit to do a constrained fit:
	photonVertexFitTree->movePointerToTheTop();
	photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
	if (!photonVertexFitTree->isValid()) {
		if(isDebug_) std::cout<<"Constrained photon vertex fit is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}

	const ParticleMass muonMass(0.1056584);
	float muonSigma = muonMass*1E-6;
		  
	photonVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
		  
	std::vector<RefCountedKinematicParticle> allChiDaughters;
	allChiDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
	allChiDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
	allChiDaughters.push_back(fitPhoton);
		  
	KinematicConstrainedVertexFitter constVertexFitter;
		  
	MultiTrackKinematicConstraint *meson_nS_mtc = new TwoTrackMassKinematicConstraint(meson_nS_mass_);
	RefCountedKinematicTree ChiTree = constVertexFitter.fit(allChiDaughters,meson_nS_mtc);

	if (ChiTree->isEmpty()) {
		if(isDebug_) std::cout<<"Constrained chi vertex fit is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}
		  
	ChiTree->movePointerToTheTop();
	RefCountedKinematicParticle fitChi = ChiTree->currentParticle();
	RefCountedKinematicVertex ChiDecayVertex = ChiTree->currentDecayVertex();

	if (!fitChi->currentState().isValid()) {
		if(isDebug_) std::cout<<"Fitted chi is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}

	//Get chi      
	float ChiM_fit  = fitChi->currentState().mass();
	float ChiPx_fit = fitChi->currentState().kinematicParameters().momentum().x();
	float ChiPy_fit = fitChi->currentState().kinematicParameters().momentum().y();
	float ChiPz_fit = fitChi->currentState().kinematicParameters().momentum().z();
	float ChiVtxX_fit = ChiDecayVertex->position().x();
	float ChiVtxY_fit = ChiDecayVertex->position().y();
	float ChiVtxZ_fit = ChiDecayVertex->position().z();
	float ChiVtxP_fit = ChiSquaredProbability((double)(ChiDecayVertex->chiSquared()),
											   (double)(ChiDecayVertex->degreesOfFreedom()));
		  
	reco::CompositeCandidate recoChi(0, math::XYZTLorentzVector(ChiPx_fit, ChiPy_fit, ChiPz_fit,
										  sqrt(ChiM_fit*ChiM_fit + ChiPx_fit*ChiPx_fit + ChiPy_fit*ChiPy_fit +
										  ChiPz_fit*ChiPz_fit)), math::XYZPoint(ChiVtxX_fit,
										  ChiVtxY_fit, ChiVtxZ_fit), 50551);
		  
	pat::CompositeCandidate patChi(recoChi);
	patChi.addUserFloat("vProb",ChiVtxP_fit);
	patChi.addUserInt("Index",indexConversion);  // this also holds the index of the current chiCand

	//get first muon
	bool child = ChiTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle fitMu1 = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Muon 1 is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}

	float mu1M_fit  = fitMu1->currentState().mass();
	float mu1Q_fit  = fitMu1->currentState().particleCharge();
	float mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
	float mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
	float mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoMu1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit, 
										 sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit + 
										 mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 13);
	pat::CompositeCandidate patMu1(recoMu1);
		  
	//get second muon
	child = ChiTree->movePointerToTheNextChild();
	RefCountedKinematicParticle fitMu2 = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Muon 2 daughter is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}

	float mu2M_fit  = fitMu2->currentState().mass();
	float mu2Q_fit  = fitMu2->currentState().particleCharge();
	float mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
	float mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
	float mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoMu2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
										 sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit + 
										 mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 13);
	pat::CompositeCandidate patMu2(recoMu2);
				  
	//Define Onia from two muons
	pat::CompositeCandidate meson_1S;
	meson_1S.addDaughter(patMu1,"muon1");
	meson_1S.addDaughter(patMu2,"muon2");	
	meson_1S.setP4(patMu1.p4()+patMu2.p4());
	  
	//get photon
	child = ChiTree->movePointerToTheNextChild();
	RefCountedKinematicParticle fitGamma = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Photon is invalid, continue"<<std::endl;
		chicCompCandRefitColl->push_back(patempty);
		continue;
	}
		  
	float gammaM_fit  = fitGamma->currentState().mass();
	float gammaPx_fit = fitGamma->currentState().kinematicParameters().momentum().x();
	float gammaPy_fit = fitGamma->currentState().kinematicParameters().momentum().y();
	float gammaPz_fit = fitGamma->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoGamma(0, math::XYZTLorentzVector(gammaPx_fit, gammaPy_fit, gammaPz_fit, 
										   sqrt(gammaM_fit*gammaM_fit + gammaPx_fit*gammaPx_fit + gammaPy_fit*gammaPy_fit +
										   gammaPz_fit*gammaPz_fit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 22);
	pat::CompositeCandidate patGamma(recoGamma);

	patChi.addDaughter(meson_1S,"dimuon");
	patChi.addDaughter(patGamma,"photon");

	chicCompCandRefitColl->push_back(patChi);	
	validchi++;
	candidates++;
  }  
  // End kinematic fit
  
  if(chicCompCandRefitColl->size() > 0 && isDebug_){
	std::cout<<"----- We have: "<<indexConversion+1<<" of ChiCands being input into Chi fitter -----"<<std::endl;
	std::cout<<"----- We have: "<<validchi<<" of patChis fitted by Chi fitter -----"<<std::endl;
	std::cout<<"----- Length of chicCompCandRefitColl is: "<<chicCompCandRefitColl->size()<<" -----"<<std::endl;
  }

  // ...ash not sorted, since we will use the best un-refitted candidate
  // now sort by vProb
  //OniaPhotonKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  //std::sort(chicCompCandRefitColl->begin(),chicCompCandRefitColl->end(), vPComparator);

  iEvent.put(std::move(chicCompCandRefitColl),product_name_); 
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaPhotonKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void OniaPhotonKinematicFit::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "chi Fit report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " fitted chi candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaPhotonKinematicFit);
