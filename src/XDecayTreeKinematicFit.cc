// -*- C++ -*-
//
// Package:    XDecayTreeKinematicFit
// Class:      XDecayTreeKinematicFit
// 
/**

 Description: Kinematic Fit for Onia + Photon + pion(charged)

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
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h> 
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"    
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"        
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
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "TLorentzVector.h"
#include <utility>
#include <string>
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class XDecayTreeKinematicFit : public edm::EDProducer {
  public:
    explicit XDecayTreeKinematicFit(const edm::ParameterSet&);
    ~XDecayTreeKinematicFit() override {};
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	bool IsTheSameMu(const reco::Candidate& tk1, const pat::Muon* tk2);
	bool IsTheSameEle(const reco::Candidate& tk1, const reco::Track& tk2);

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;
	void endJob() override;

	edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> meson_nS_Label;
	edm::EDGetTokenT<std::vector<pat::PackedCandidate>> track_Label; //only work if it is std::vector<pat::PackedCandidate>...
	edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
	edm::EDGetTokenT<reco::VertexCollection> thePVs_;

    double meson_nS_mass_;
	std::string chi_product_name_;
	std::string XCand_product_name_;
    std::string X_product_name_;
	bool isDebug_;
	bool resolvePVAmbiguity_;
	
	double deltaMass_;
	double dzMax_;
	double deltaR_pi_;
	
	int chi_fitted_candidates;
	int X_candidates;
	int X_fitted_candidates;
	
	int fit_fail;
	int daughter_fail;
    

    template<typename T>
    struct GreaterByVProb {
      typedef T first_argument_type;
      typedef T second_argument_type;
      bool operator()( const T & t1, const T & t2 ) const {
            return t1.userFloat("vProb") > t2.userFloat("vProb");
      }
    };
};

XDecayTreeKinematicFit::XDecayTreeKinematicFit(const edm::ParameterSet& iConfig) {
  chi_Label       = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("chi_cand"));
  meson_nS_Label  = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("meson_nS_cand"));
  track_Label     = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter< edm::InputTag>("track"));
  thebeamspot_    = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
  thePVs_         = consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertexTag")),
  
  meson_nS_mass_      = iConfig.getParameter<double>("meson_nS_mass");
  chi_product_name_   = iConfig.getParameter<std::string>("chi_product_name");
  deltaMass_          = iConfig.getParameter<double>("deltaMass"),
  dzMax_              = iConfig.getParameter<double>("dzmax"),
  deltaR_pi_          = iConfig.getParameter<double>("deltaR_pi"),
  XCand_product_name_ = iConfig.getParameter<std::string>("XCand_product_name");
  X_product_name_     = iConfig.getParameter<std::string>("X_product_name");
  produces<pat::CompositeCandidateCollection>(chi_product_name_);
  produces<pat::CompositeCandidateCollection>(XCand_product_name_);
  produces<pat::CompositeCandidateCollection>(X_product_name_);
  isDebug_ = iConfig.getParameter<bool>("is_Debug");
  resolvePVAmbiguity_ = iConfig.getParameter<bool>("resolvePVAmbiguity");
  chi_fitted_candidates = 0;
  X_candidates = 0;
  X_fitted_candidates = 0;
  fit_fail = 0;
  daughter_fail = 0;
}

// ------------ method called to produce the data  ------------
void XDecayTreeKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // Grab parameters
  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByToken(chi_Label, chiCandHandle);
  
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();
  
  reco::Vertex theBeamSpotV;
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());
  
  reco::Vertex thePrimaryV;
  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  // if(priVtxs->begin() != priVtxs->end())
    // thePrimaryV = reco::Vertex(*(priVtxs->begin()));
  // else
    // thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  
  //Kinematic refit collection
  std::unique_ptr<pat::CompositeCandidateCollection> chicCompCandRefitColl(new pat::CompositeCandidateCollection);
  std::unique_ptr<pat::CompositeCandidateCollection> XCandCompCandRefitColl(new pat::CompositeCandidateCollection);
  std::unique_ptr<pat::CompositeCandidateCollection> XCompCandRefitColl(new pat::CompositeCandidateCollection);
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> TTrkBuilder; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTrkBuilder); 
  
  int indexChiCand=-1;
  
  //Internal counter
  int validchi = 0;
  int preX = 0;
  int validX = 0;

  for (pat::CompositeCandidateCollection::const_iterator chiCand=chiCandHandle->begin();chiCand!=chiCandHandle->end();++chiCand) { 
    indexChiCand++;
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
    MuMuTT.push_back((*TTrkBuilder).build(&JpsiTk[0]));
    MuMuTT.push_back((*TTrkBuilder).build(&JpsiTk[1]));
      
    std::vector<reco::TransientTrack> EETT;
    EETT.push_back((*TTrkBuilder).build(convTracks[0]));
    EETT.push_back((*TTrkBuilder).build(convTracks[1]));
      
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
      
    //If the process is stopped because kinematic fit / cut failed, push_back an empty CompositeCandidate into XCompCandRefitColl, thus, each chiCand would have its own refit1P (even though it might be empty) to access in the tuplizer
    if (!photonVertexFitTree->isValid()) {
		if(isDebug_) std::cout<<">>>>> Photon vertex fit is invalid, continue <<<<<"<<std::endl;
		continue;
	}
	  
	// now apply Photon mass constraint
	KinematicParticleFitter csFitterPhoton;
	KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);

	// add mass constraint to the photon fit to do a constrained fit:
	photonVertexFitTree->movePointerToTheTop();
	photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
	if (!photonVertexFitTree->isValid()) {
		if(isDebug_) std::cout<<">>>>> Constrained photon vertex fit is invalid, continue <<<<<"<<std::endl;
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
		if(isDebug_) std::cout<<">>>>> Constrained chi vertex fit is invalid, continue <<<<<"<<std::endl;
		continue;
	}
		  
	ChiTree->movePointerToTheTop();
	RefCountedKinematicParticle fitChi = ChiTree->currentParticle();
	RefCountedKinematicVertex ChiDecayVertex = ChiTree->currentDecayVertex();

	if (!fitChi->currentState().isValid()) {
		if(isDebug_) std::cout<<">>>>> Fitted chi is invalid, continue <<<<<"<<std::endl;
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
		  
	reco::CompositeCandidate recoChi_3trkfit(0, math::XYZTLorentzVector(ChiPx_fit, ChiPy_fit, ChiPz_fit,
										  sqrt(ChiM_fit*ChiM_fit + ChiPx_fit*ChiPx_fit + ChiPy_fit*ChiPy_fit +
										  ChiPz_fit*ChiPz_fit)), math::XYZPoint(ChiVtxX_fit,
										  ChiVtxY_fit, ChiVtxZ_fit), 50551);
		  
	pat::CompositeCandidate patChi_3trkfit(recoChi_3trkfit);
	patChi_3trkfit.addUserFloat("vProb",ChiVtxP_fit);
	patChi_3trkfit.addUserInt("Index",indexChiCand);  // this also holds the index of the current chiCand
	TLorentzVector rf_chi_p4;
	rf_chi_p4.SetXYZM(ChiPx_fit,ChiPx_fit,ChiPx_fit,ChiM_fit);																		   
	//get first muon
	bool child = ChiTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle fitMu1_3trkfit = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Muon 1 is invalid in chi 3-track fit, continue"<<std::endl;
		continue;
	}

	float mu1M_3trkfit  = fitMu1_3trkfit->currentState().mass();
	float mu1Q_3trkfit  = fitMu1_3trkfit->currentState().particleCharge();
	float mu1Px_3trkfit = fitMu1_3trkfit->currentState().kinematicParameters().momentum().x();
	float mu1Py_3trkfit = fitMu1_3trkfit->currentState().kinematicParameters().momentum().y();
	float mu1Pz_3trkfit = fitMu1_3trkfit->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoMu1_3trkfit(mu1Q_3trkfit, math::XYZTLorentzVector(mu1Px_3trkfit, mu1Py_3trkfit, mu1Pz_3trkfit, 
										 sqrt(mu1M_3trkfit*mu1M_3trkfit + mu1Px_3trkfit*mu1Px_3trkfit + mu1Py_3trkfit*mu1Py_3trkfit + 
										 mu1Pz_3trkfit*mu1Pz_3trkfit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 13);
	pat::CompositeCandidate patMu1_3trkfit(recoMu1_3trkfit);
		  
	//get second muon
	child = ChiTree->movePointerToTheNextChild();
	RefCountedKinematicParticle fitMu2_3trkfit = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Muon 2 is invalid in chi 3-track fit, continue"<<std::endl;
		continue;
	}

	float mu2M_3trkfit  = fitMu2_3trkfit->currentState().mass();
	float mu2Q_3trkfit  = fitMu2_3trkfit->currentState().particleCharge();
	float mu2Px_3trkfit = fitMu2_3trkfit->currentState().kinematicParameters().momentum().x();
	float mu2Py_3trkfit = fitMu2_3trkfit->currentState().kinematicParameters().momentum().y();
	float mu2Pz_3trkfit = fitMu2_3trkfit->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoMu2_3trkfit(mu2Q_3trkfit, math::XYZTLorentzVector(mu2Px_3trkfit, mu2Py_3trkfit, mu2Pz_3trkfit, 
										 sqrt(mu2M_3trkfit*mu2M_3trkfit + mu2Px_3trkfit*mu2Px_3trkfit + mu2Py_3trkfit*mu2Py_3trkfit + 
										 mu2Pz_3trkfit*mu2Pz_3trkfit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 13);
	pat::CompositeCandidate patMu2_3trkfit(recoMu2_3trkfit);
				  
	//Define Onia from two muons
	pat::CompositeCandidate meson_1S_3trkfit;
	meson_1S_3trkfit.addDaughter(patMu1_3trkfit,"muon1");
	meson_1S_3trkfit.addDaughter(patMu2_3trkfit,"muon2");	
	meson_1S_3trkfit.setP4(patMu1_3trkfit.p4() + patMu2_3trkfit.p4());
	  
	//get photon
	child = ChiTree->movePointerToTheNextChild();
	RefCountedKinematicParticle fitGamma_3trkfit = ChiTree->currentParticle();
	if (!child) {
		if(isDebug_) std::cout<<"Photon is invalid in chi 3-track fit, continue"<<std::endl;
		continue;
	}
		  
	float gammaM_3trkfit  = fitGamma_3trkfit->currentState().mass();
	float gammaPx_3trkfit = fitGamma_3trkfit->currentState().kinematicParameters().momentum().x();
	float gammaPy_3trkfit = fitGamma_3trkfit->currentState().kinematicParameters().momentum().y();
	float gammaPz_3trkfit = fitGamma_3trkfit->currentState().kinematicParameters().momentum().z();
	reco::CompositeCandidate recoGamma_3trkfit(0, math::XYZTLorentzVector(gammaPx_3trkfit, gammaPy_3trkfit, gammaPz_3trkfit, 
										   sqrt(gammaM_3trkfit*gammaM_3trkfit + gammaPx_3trkfit*gammaPx_3trkfit + gammaPy_3trkfit*gammaPy_3trkfit +
										   gammaPz_3trkfit*gammaPz_3trkfit)), math::XYZPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit), 22);
	pat::CompositeCandidate patGamma_3trkfit(recoGamma_3trkfit);
	
	//locate the PV closed to the fitted chi, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOnia2MuMuPAT
	// const reco::Vertex *dimu_vertex = dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("dimuon"))->userData<reco::Vertex>("commonVertex");
	// reco::Candidate::LorentzVector dimu_p4 = dynamic_cast<const pat::CompositeCandidate*>(chiCand->daughter("dimuon"))->p4();
	unsigned int chi_PV_index = 0; //default PV is the first one in PV collection -> highest pt sum
	float minDz = 999999;
	if(resolvePVAmbiguity_){
		TwoTrackMinimumDistance ttmd;
		bool status = ttmd.calculate(GlobalTrajectoryParameters(GlobalPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit),GlobalVector(ChiPx_fit,ChiPy_fit,ChiPz_fit),TrackCharge(0),&(*magneticField)), GlobalTrajectoryParameters(GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),GlobalVector(bs.dxdz(), bs.dydz(), 1),TrackCharge(0),&(*magneticField)));
		float extrapZ = -9E20;
		if(status) 
			extrapZ=ttmd.points().first.z();
		for(unsigned int ipv=0; ipv<priVtxs->size(); ipv++){
			float deltaZ = fabs(extrapZ - priVtxs->at(ipv).position().z()) ;
			if (deltaZ < minDz){
				minDz = deltaZ;    
				chi_PV_index = ipv;
			}
		}
	}
	const reco::Vertex chi_PV = priVtxs->at(chi_PV_index); //the PV which closest to J/psi in dz becomes the selected PV

	patChi_3trkfit.addDaughter(meson_1S_3trkfit,"dimuon");
	patChi_3trkfit.addDaughter(patGamma_3trkfit,"photon");
	patChi_3trkfit.addUserInt("Index",indexChiCand);  // this also holds the index of the current chiCand
	patChi_3trkfit.addUserInt("PV_index",chi_PV_index);
	patChi_3trkfit.addUserFloat("dz_PVdimu",minDz);

	chicCompCandRefitColl->push_back(patChi_3trkfit);	
	validchi++;
	chi_fitted_candidates++;

	edm::Handle<std::vector<pat::PackedCandidate>> trackHandle;
	iEvent.getByToken(track_Label, trackHandle);
	
	int trk_index = -1;
	for (std::vector<pat::PackedCandidate>::const_iterator iTrack = trackHandle->begin(); iTrack != trackHandle->end(); ++iTrack){ //pair the fitted chi_c with a track
		trk_index++;
		//quality cuts track, charged
		if (iTrack->charge() == 0) continue;
		if (fabs(iTrack->pdgId()) != 211) continue;
		if (iTrack->pt() < 0.6) continue;
		if (!(iTrack->trackHighPurity())) continue;
		if (iTrack->numberOfPixelHits() < 1) continue;
		if (iTrack->numberOfHits() < 5) continue;
		TLorentzVector X_prefit_p4;
		TLorentzVector itrk_p4;
		itrk_p4.SetPtEtaPhiE(iTrack->pt(),iTrack->eta(),iTrack->phi(),iTrack->energy());
		X_prefit_p4 = rf_chi_p4 + itrk_p4;
		if(X_prefit_p4.M() < 3.872 - deltaMass_ || X_prefit_p4.M() > 3.872 + deltaMass_) continue;
		if(X_prefit_p4.DeltaR(itrk_p4) > deltaR_pi_) continue;
		
		//Now let's checks if the track we are selecting is the muon or conversion electron we already selected
		const pat::Muon *mu0 = dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon1"));
		const pat::Muon *mu1 = dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon2"));
		if (IsTheSameMu(*iTrack,mu0) || IsTheSameMu(*iTrack,mu1) ) continue;
		if (IsTheSameEle(*iTrack,tk0) || IsTheSameEle(*iTrack,tk1) ) continue;
		
		//make sure track comes from the same PV as the dimuon does
		if(fabs(iTrack->vertexRef()->x()-chi_PV.x()) > 0.002 || fabs(iTrack->vertexRef()->y()-chi_PV.y()) > 0.002 || fabs(iTrack->vertexRef()->z()-chi_PV.z()) > 0.002)
			continue;
		if(iTrack->fromPV() < 2) //track association to the vertexRef(), used for isolation calculations, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
			continue;
		
		//DCA
		// TrajectoryStateClosestToPoint chiTS = t_tks[0].impactPointTSCP();
		// TrajectoryStateClosestToPoint pionTS = t_tks[1].impactPointTSCP();
		// float dca = 1E20;
		// if (mu1TS.isValid() && mu2TS.isValid()) {
			// ClosestApproachInRPhi cApp;
			// cApp.calculate(mu1TS.theState(), mu2TS.theState());
		// if (cApp.status() ) dca = cApp.distance();
		// }
		// myCand.addUserFloat("DCA", dca );
		
		//DCA
		int pion_charge = iTrack->charge();
		FreeTrajectoryState chiFTS(GlobalPoint(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit),GlobalVector(ChiPx_fit,ChiPy_fit,ChiPz_fit),TrackCharge(0),&(*magneticField));
		FreeTrajectoryState trkFTS(GlobalPoint(iTrack->vertexRef()->x(), iTrack->vertexRef()->y(), iTrack->vertexRef()->z()),GlobalVector(iTrack->px(),iTrack->py(),iTrack->pz()),TrackCharge(pion_charge),&(*magneticField));
		float dca = 1E20;
		ClosestApproachInRPhi cApp; //closest approach in the transverse plane
		cApp.calculate(chiFTS, trkFTS);
		if(cApp.status())
			dca = cApp.distance();
		
		//dz
		double dz = iTrack->dz(reco::Candidate::Point(ChiVtxX_fit, ChiVtxY_fit, ChiVtxZ_fit)); //IP using the fitted chi vertex as a reference
		
		pat::CompositeCandidate XCand;
		XCand.addDaughter(*chiCand,"chi");
		XCand.addDaughter(*iTrack,"pionpm");
		XCand.addUserInt("Index",indexChiCand);  // this also holds the index of the current chiCand
		XCand.addUserInt("PV_association",iTrack->fromPV()); 
		XCand.addUserFloat("dz_chitrk",dz);
		XCand.addUserFloat("DCA_chitrk", dca);
		XCand.addUserFloat("HcalFrac",iTrack->hcalFraction());
		XCand.addUserInt("isIsoChgHad", iTrack->isIsolatedChargedHadron());
		reco::Candidate::LorentzVector vX;
		vX.SetPxPyPzE(X_prefit_p4.Px(),X_prefit_p4.Py(),X_prefit_p4.Pz(),X_prefit_p4.E());
		XCand.setP4(vX);
		
		XCandCompCandRefitColl->push_back(XCand);	
		preX++;
		X_candidates++;
  
		reco::TransientTrack pionpmTT = (*TTrkBuilder).build(iTrack->pseudoTrack());
		const ParticleMass pionpmMass(0.1395704); //Assuming this track is a charged pion n: 0.1349766 pm: 0.1395704
		float pionpmSigma = pionpmMass*1E-6;

		//Currently we perform 4-track fit with di-muon J/psi mass constraints , maybe can do 2-track fit with chi_c mass constraints?
		std::vector<RefCountedKinematicParticle> allXDaughters;
		allXDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
		allXDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
		allXDaughters.push_back(fitPhoton);
		allXDaughters.push_back(pFactory.particle (pionpmTT, pionpmMass, float(0), float(0), pionpmSigma));

		RefCountedKinematicTree XTree = constVertexFitter.fit(allXDaughters,meson_nS_mtc);

		if (XTree->isEmpty()) {
			if(isDebug_) std::cout<<">>>>> Constrained X vertex fit is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			fit_fail++;
			continue;
		}

		XTree->movePointerToTheTop();
		RefCountedKinematicParticle fittedX = XTree->currentParticle();
		RefCountedKinematicVertex XDecayVertex = XTree->currentDecayVertex();

		if (fittedX->currentState().isValid()) {
			if(isDebug_) std::cout<<">>>>> Fitted X is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			fit_fail++;
			continue;
		}

		float XM_fit  = fittedX->currentState().mass();
		float XPx_fit = fittedX->currentState().kinematicParameters().momentum().x();
		float XPy_fit = fittedX->currentState().kinematicParameters().momentum().y();
		float XPz_fit = fittedX->currentState().kinematicParameters().momentum().z();
		float XVtxX_fit = XDecayVertex->position().x();
		float XVtxY_fit = XDecayVertex->position().y();
		float XVtxZ_fit = XDecayVertex->position().z();
		float XVtxP_fit = ChiSquaredProbability((double)(XDecayVertex->chiSquared()),
												   (double)(XDecayVertex->degreesOfFreedom()));
			  
		reco::CompositeCandidate recoX(0, math::XYZTLorentzVector(XPx_fit, XPy_fit, XPz_fit,
											  sqrt(XM_fit*XM_fit + XPx_fit*XPx_fit + XPy_fit*XPy_fit +
											  XPz_fit*XPz_fit)), math::XYZPoint(XVtxX_fit,
											  XVtxY_fit, XVtxZ_fit), 88888);
			  
		pat::CompositeCandidate patX(recoX);
		patX.addUserFloat("vProb",XVtxP_fit);
		patX.addUserInt("Index",indexChiCand);  // this also holds the index of the current chiCand

		//get first muon
		bool child = XTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle fitMu1 = XTree->currentParticle();
		if (!child) {
			if(isDebug_) std::cout<<">>>>> Muon 1 is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			daughter_fail++;
			continue;
		}

		float mu1M_fit  = fitMu1->currentState().mass();
		float mu1Q_fit  = fitMu1->currentState().particleCharge();
		float mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
		float mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
		float mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
		reco::CompositeCandidate recoMu1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit, 
											 sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit + 
											 mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);
		pat::CompositeCandidate patMu1(recoMu1);
			  
		//get second muon
		child = XTree->movePointerToTheNextChild();
		RefCountedKinematicParticle fitMu2 = XTree->currentParticle();
		if (!child) {
			if(isDebug_) std::cout<<">>>>> Muon 2 is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			daughter_fail++;
			continue;
		}

		float mu2M_fit  = fitMu2->currentState().mass();
		float mu2Q_fit  = fitMu2->currentState().particleCharge();
		float mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
		float mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
		float mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
		reco::CompositeCandidate recoMu2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
											 sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit + 
											 mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 13);
		pat::CompositeCandidate patMu2(recoMu2);
		
		//get photon
		child = XTree->movePointerToTheNextChild();
		RefCountedKinematicParticle fitGamma = XTree->currentParticle();
		if (!child) {
			if(isDebug_) std::cout<<">>>>> Photon is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			daughter_fail++;
			continue;
		}
			  
		float gammaM_fit  = fitGamma->currentState().mass();
		float gammaPx_fit = fitGamma->currentState().kinematicParameters().momentum().x();
		float gammaPy_fit = fitGamma->currentState().kinematicParameters().momentum().y();
		float gammaPz_fit = fitGamma->currentState().kinematicParameters().momentum().z();
		reco::CompositeCandidate recoGamma(0, math::XYZTLorentzVector(gammaPx_fit, gammaPy_fit, gammaPz_fit, 
											   sqrt(gammaM_fit*gammaM_fit + gammaPx_fit*gammaPx_fit + gammaPy_fit*gammaPy_fit +
											   gammaPz_fit*gammaPz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 22);
		pat::CompositeCandidate patGamma(recoGamma);
		
		//get pion
		child = XTree->movePointerToTheNextChild();
		RefCountedKinematicParticle fitpion = XTree->currentParticle();
		if (!child) {
			if(isDebug_) std::cout<<">>>>> Pion is invalid with #"<<trk_index<<" track, continue <<<<<"<<std::endl;
			daughter_fail++;
			continue;
		}

		float pionM_fit  = fitpion->currentState().mass();
		float pionQ_fit  = fitpion->currentState().particleCharge();
		float pionPx_fit = fitpion->currentState().kinematicParameters().momentum().x();
		float pionPy_fit = fitpion->currentState().kinematicParameters().momentum().y();
		float pionPz_fit = fitpion->currentState().kinematicParameters().momentum().z();
		reco::CompositeCandidate recopion(pionQ_fit, math::XYZTLorentzVector(pionPx_fit, pionPy_fit, pionPz_fit, 
											 sqrt(pionM_fit*pionM_fit + pionPx_fit*pionPx_fit + pionPy_fit*pionPy_fit + 
											 pionPz_fit*pionPz_fit)), math::XYZPoint(XVtxX_fit, XVtxY_fit, XVtxZ_fit), 211);
		pat::CompositeCandidate patpion(recopion);
					  
		//Define 1S Onia from two muons
		pat::CompositeCandidate meson_1S;
		meson_1S.addDaughter(patMu1,"muon1");
		meson_1S.addDaughter(patMu2,"muon2");	
		meson_1S.setP4(patMu1.p4() + patMu2.p4());
		
		//Define 1P Onia from two muons and photon
		pat::CompositeCandidate patChi ;
		patChi.addDaughter(meson_1S,"dimuon");
		patChi.addDaughter(patGamma,"photon");
		patChi.setP4(meson_1S.p4() + patGamma.p4());
		  
		//Define X from chi and pion
		patX.addDaughter(patChi,"chi");
		patX.addDaughter(patpion,"pion_pm");
		patX.setP4(patChi.p4() + patpion.p4());
		
		XCompCandRefitColl->push_back(patX);
		validX++;
		std::cout<<">>>>> We got 1 charged X candidate! <<<<<"<<std::endl;
		X_fitted_candidates++;
	}
  }  
    
  // End kinematic fit
  if(XCompCandRefitColl->size() > 0 && isDebug_){
	std::cout<<">>>>> We have: "<<indexChiCand+1<<" of ChiCands <<<<<"<<std::endl;
	std::cout<<">>>>> We have: "<<validchi<<" of fitted chis <<<<<"<<std::endl;
	std::cout<<">>>>> We have: "<<validX<<" of patXs <<<<<"<<std::endl;
	std::cout<<">>>>> Length of XCompCandRefitColl is: "<<XCompCandRefitColl->size()<<" <<<<<"<<std::endl;
	std::cout<<">>>>> END OF KINEMATIC FITS <<<<<"<<std::endl;
  }

  // ...ash not sorted, since we will use the best un-refitted candidate
  // now sort by vProb
  //XDecayTreeKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  //std::sort(XCompCandRefitColl->begin(),XCompCandRefitColl->end(), vPComparator);

  iEvent.put(std::move(chicCompCandRefitColl),chi_product_name_); 
  iEvent.put(std::move(XCandCompCandRefitColl),XCand_product_name_);
  iEvent.put(std::move(XCompCandRefitColl),X_product_name_); 
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void XDecayTreeKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool XDecayTreeKinematicFit::IsTheSameMu(const reco::Candidate& tk1, const pat::Muon* tk2){
  double DeltaEta = fabs(tk1.eta() - tk2->eta());
  double DeltaP   = fabs(tk1.p() - tk2->p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) 
	return true;
  else
	return false;
}

bool XDecayTreeKinematicFit::IsTheSameEle(const reco::Candidate& tk1, const reco::Track& tk2){
  double DeltaEta = fabs(tk1.eta()-tk2.eta());
  double DeltaP   = fabs(tk1.p()-tk2.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) 
	return true;
  else
	return false;
}

void XDecayTreeKinematicFit::endJob(){
  std::cout << "-------------------------------" << std::endl;
  std::cout << "X DecayTree Kinematic Fit report:" << std::endl;
  std::cout << "-------------------------------" << std::endl;
  std::cout << "Found " << chi_fitted_candidates << " fitted chi candidates from 3-track fit" << std::endl;
  std::cout << "Found " << X_candidates << " X candidates" << std::endl;
  std::cout << "Found " << X_fitted_candidates << " fitted X candidates from 4-track fit" << std::endl;
  std::cout << fit_fail << " X candidates failed during fit" << std::endl;
  std::cout << daughter_fail << " X candidates failed during daughter finding" << std::endl;
  std::cout << "-------------------------------" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(XDecayTreeKinematicFit);
