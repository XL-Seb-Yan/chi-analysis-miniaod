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

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "TLorentzVector.h"
#include <utility>
#include <string>

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
	bool IsTheSameMu(const reco::Candidate& tk1, const reco::TrackRef& tk2);
	bool IsTheSameEle(const reco::Candidate& tk1, const reco::TrackRef& tk2);
	void endJob();

  private:
    void produce(edm::Event&, const edm::EventSetup&) override;

	edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_Label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> meson_nS_Label;
	edm::EDGetTokenT<pat::PackedCandidate> track_Label;
    double meson_nS_mass_;
    std::string product_name_;
	
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

XDecayTreeKinematicFit::XDecayTreeKinematicFit(const edm::ParameterSet& iConfig) {
  chi_Label       = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("chi_cand"));
  meson_nS_Label  = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("Onia2MuMuFiltered"));
  track_Label     = consumes<pat::PackedCandidate>(iConfig.getParameter< edm::InputTag>("track")),
  meson_nS_mass_  = iConfig.getParameter<double>("meson_nS_mass");
  product_name_   = iConfig.getParameter<std::string>("product_name");
  produces<pat::CompositeCandidateCollection>(product_name_);
  candidates = 0;
}

// ------------ method called to produce the data  ------------
void XDecayTreeKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // Grab parameters
  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByToken(chi_Label, chiCandHandle);
  
  //Kinematic refit collection
  std::unique_ptr<pat::CompositeCandidateCollection> XCompCandRefitColl(new pat::CompositeCandidateCollection);
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> TTrkBuilder; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTrkBuilder); 
  
  int indexChiCand=-1;

  for (pat::CompositeCandidateCollection::const_iterator chiCand=chiCandHandle->begin();chiCand!=chiCandHandle->end();++chiCand) { 
    indexChiCand++;
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
      
    if (photonVertexFitTree->isValid()) {
	  
      // now apply Photon mass constraint
      KinematicParticleFitter csFitterPhoton;
      KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);

      // add mass constraint to the photon fit to do a constrained fit:
      photonVertexFitTree->movePointerToTheTop();
      photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
	  
      if (photonVertexFitTree->isValid()) {

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

		if (!ChiTree->isEmpty()) {
			  
		  ChiTree->movePointerToTheTop();
		  RefCountedKinematicParticle fitChi = ChiTree->currentParticle();
		  RefCountedKinematicVertex ChiDecayVertex = ChiTree->currentDecayVertex();
		  
		  if (fitChi->currentState().isValid()) { //make sure the fitted chi is valid 
		  
			  edm::Handle<pat::PackedCandidate> trackHandle;
			  iEvent.getByToken(track_Label, trackHandle);
			  
			  for (pat::PackedCandidate::const_iterator iTrack = trackHandle->begin(); iTrack != trackHandle->end(); ++iTrack){ //pair the fitted chi_c with a track
				//quality cuts track
				if (iTrack->charge() == 0) continue;
				if (fabs(iTrack->pdgId()) != 211) continue;
				if (iTrack->pt() < 0.6) continue;
				if (!(iTrack->trackHighPurity())) continue;
				if (iTrack->numberOfPixelHits() < 1) continue;
				if (iTrack->numberOfHits() < 5) continue;
				
				//Now let's checks if the track we are selecting is the muon or conversion electron we already selected
				reco::TrackRef mu0 = dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon1"))->innerTrack();
				reco::TrackRef mu1 = dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon2"))->innerTrack();
				if (IsTheSameMu(*iTrack,mu0) || IsTheSameMu(*iTrack,mu1) ) continue;
				if (IsTheSameEle(*iTrack,tk0) || IsTheSameEle(*iTrack,tk1) ) continue;
			
				const ParticleMass pionpmMass(0.1395704); //Assuming this track is a charged pion
				float pionpmSigma = pionpmMass*1E-6;
				
				// Do a simple check of chi+pion invariant mass before kinematic fit
				TLorentzVector chi_p4, pionpm_p4;
				chi_p4.SetXYZM(fitChi->currentState().kinematicParameters().momentum().x(),
								fitChi->currentState().kinematicParameters().momentum().y(),
								fitChi->currentState().kinematicParameters().momentum().z(),
								fitChi->currentState().mass());

				// TLorentzVector mu1_p4, mu2_p4, photon_p4, pionpm_p4;
				// mu1_p4.SetXYZM(chiCand->daughter("dimuon")->daughter("muon1")->px(),
								// chiCand->daughter("dimuon")->daughter("muon1")->py(),
								// chiCand->daughter("dimuon")->daughter("muon1")->pz(),
								// muonMass);
								
				// mu2_p4.SetXYZM(chiCand->daughter("dimuon")->daughter("muon2")->px(),
								// chiCand->daughter("dimuon")->daughter("muon2")->py(),
								// chiCand->daughter("dimuon")->daughter("muon2")->pz(),
								// muonMass);
								
				// photon_p4.SetXYZM(fitPhoton->currentState().kinematicParameters().momentum().x(),
									// fitPhoton->currentState().kinematicParameters().momentum().y(),
									// fitPhoton->currentState().kinematicParameters().momentum().z(),
									// zero_mass);
									
				pionpm_p4.SetXYZM(iTrack->px(),iTrack->py(),iTrack->pz(),pionpmMass);
				
				TLorentzVector X_prefit_p4;
				X_prefit_p4 = chi_p4 + pionpm_p4;
				if(X_prefit_p4.M() < 3.872-0.5 || X_prefit_p4.M() > 3.872+0.5) continue;
				reco::TransientTrack pionpmTT = (*TTrkBuilder).build(iTrack->pseudoTrack());
			  
			    //Currently we perform 4-track fit with di-muon J/psi mass constraints , maybe can do 2-track fit with chi_c mass constraints?
				std::vector<RefCountedKinematicParticle> allXDaughters;
				allXDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
				allXDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
				allXDaughters.push_back(fitPhoton);
				allXDaughters.push_back(pFactory.particle (pionpmTT, pionpmMass, float(0), float(0), pionpmSigma));
			  
				RefCountedKinematicTree XTree = constVertexFitter.fit(allXDaughters,meson_nS_mtc);
			
				if (!XTree->isEmpty()) {
			  
					XTree->movePointerToTheTop();
					RefCountedKinematicParticle fittedX = XTree->currentParticle();
					RefCountedKinematicVertex XDecayVertex = XTree->currentDecayVertex();
			
					if (fittedX->currentState().isValid()) { //make sure the fitted X is valid 
		   
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
							std::cout<<"Missing the first daughter (muon1) of this fitted X, exiting!"<<std::endl;
							break;
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
							std::cout<<"Missing the second daughter (muon2) of this fitted X, exiting!"<<std::endl;
							break;
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
						child = ChiTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitGamma = ChiTree->currentParticle();
						if (!child) {
							std::cout<<"Missing the third daughter (photon) of this fitted X, exiting!"<<std::endl;
							break;
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
						child = XTree->movePointerToTheFirstChild();
						RefCountedKinematicParticle fitpion = XTree->currentParticle();
						if (!child) {
							std::cout<<"Missing the fourth daughter (pion) of this fitted X, exiting!"<<std::endl;
							break;
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
						candidates++;
					}
				}
			  }
		  }
		}
      }
    }
	else{
		reco::CompositeCandidate recoempty(0, math::XYZTLorentzVector(0,0,0,0), math::XYZPoint(0,0,0), -99);
		pat::CompositeCandidate patempty(recoempty);
		XCompCandRefitColl->push_back(patempty);
	}
	//if photon vertex is not valid, nothing is write ti the XCompCandRefitColl, thus, the length of refit1P is not equal to chiCand
  }  
  // End kinematic fit
  
  if(XCompCandRefitColl->size() > 0){
	std::cout<<"======= We have: "<<indexChiCand+1<<" of ChiCands ======="<<std::endl;
	std::cout<<"======= We have: "<<XCompCandRefitColl->size()<<" of patXs ======="<<std::endl;
  }

  // ...ash not sorted, since we will use the best un-refitted candidate
  // now sort by vProb
  //XDecayTreeKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  //std::sort(XCompCandRefitColl->begin(),XCompCandRefitColl->end(), vPComparator);

  iEvent.put(std::move(XCompCandRefitColl),product_name_); 
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void XDecayTreeKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool XDecayTreeKinematicFit::IsTheSameMu(const reco::Candidate& tk1, const reco::TrackRef& tk2){
  double DeltaEta = fabs(tk1.eta()-tk2.eta());
  double DeltaP   = fabs(tk1.p()-tk2.p());
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
  std::cout << "###########################" << std::endl;
  std::cout << "X DecayTree Kinematic Fit report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " fitted X candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(XDecayTreeKinematicFit);
