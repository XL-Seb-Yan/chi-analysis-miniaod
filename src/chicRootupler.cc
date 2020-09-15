#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/PatCandidates/interface/UserData.h> 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

class chicRootupler:public edm::EDAnalyzer {
      public:
		explicit chicRootupler(const edm::ParameterSet &);
		~chicRootupler() override {};
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
		void analyze(const edm::Event &, const edm::EventSetup &) override;

		std::string file_name;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> meson_nS_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit1P_;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> X_;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> refitX_;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;
		bool isMC_;
		edm::EDGetTokenT<reco::GenParticleCollection>       genCands_;
		
		TTree *chi_tree;
		
		//For evt info
		UInt_t    run;
        ULong64_t event;
        UInt_t    lumiblock;
		UInt_t    numPrimaryVertices;
		UInt_t    trigger;
		Int_t     chi_cand_N;
		Int_t     chi_fit_N;
		Int_t     X_cand_N;
		Int_t     X_fit_N;
		
		//For Chi_cand(before kinematic fit)
		std::vector<Double_t> chi_pt_arr;
		std::vector<Double_t> chi_eta_arr;
		std::vector<Double_t> chi_phi_arr;
		std::vector<Double_t> chi_mass_arr;
		std::vector<Double_t> dimuon_pt_arr;
		std::vector<Double_t> dimuon_eta_arr;
		std::vector<Double_t> dimuon_phi_arr;
		std::vector<Double_t> dimuon_mass_arr;
		std::vector<Double_t> muonP_pt_arr;
		std::vector<Double_t> muonP_eta_arr;
		std::vector<Double_t> muonP_phi_arr;
		std::vector<Double_t> muonP_mass_arr;
		std::vector<Double_t> muonM_pt_arr;
		std::vector<Double_t> muonM_eta_arr;
		std::vector<Double_t> muonM_phi_arr;
		std::vector<Double_t> muonM_mass_arr;
		std::vector<Double_t> photon_pt_arr;
		std::vector<Double_t> photon_eta_arr;
		std::vector<Double_t> photon_phi_arr;
		std::vector<Double_t> photon_mass_arr;
		
		std::vector<Double_t> ele_lowerPt_pt_arr;
		std::vector<Double_t> ele_higherPt_pt_arr;
		std::vector<Double_t> ctpv_arr;
		std::vector<Double_t> ctpv_error_arr;
		std::vector<Double_t> conv_vertex_arr;
		std::vector<Double_t> chi_dz_arr;
		std::vector<UInt_t>   photon_flags_arr;
		std::vector<Double_t> invm1P_arr;
		//Double_t meson_1P_nsigma; //This can be access later on with dimuon_p4.Rapidity() and dimuon_p4.M()

		//For fitted Chi
		std::vector<Double_t> rf1P_chi_pt_arr;
		std::vector<Double_t> rf1P_chi_eta_arr;
		std::vector<Double_t> rf1P_chi_phi_arr;
		std::vector<Double_t> rf1P_chi_mass_arr;
        std::vector<Double_t> probFit1P_arr;
		std::vector<Int_t> rf1P_chi_index_arr;
		
		//For X Cand, since info for chi is already stored, only need to store the track info and X cand info 
		std::vector<Double_t> X_pt_arr;
		std::vector<Double_t> X_eta_arr;
		std::vector<Double_t> X_phi_arr;
		std::vector<Double_t> X_mass_arr;
		std::vector<Double_t> trk_pt_arr;
		std::vector<Double_t> trk_eta_arr;
		std::vector<Double_t> trk_phi_arr;
		std::vector<Double_t> trk_e_arr;
		std::vector<Int_t> X_index_arr;
		
		//For fitted X
		std::vector<Double_t> rfX_pt_arr;
		std::vector<Double_t> rfX_eta_arr;
		std::vector<Double_t> rfX_phi_arr;
		std::vector<Double_t> rfX_mass_arr;
        std::vector<Double_t> probFitX_arr;
		std::vector<Int_t> rfX_index_arr;
		
		//MC
		TLorentzVector gen_X_p4;
		TLorentzVector gen_chi_p4;
		TLorentzVector gen_meson_nS_p4;
        TLorentzVector gen_dimuon_p4;
		TLorentzVector gen_photon_p4;
		TLorentzVector gen_muonP_p4;
		TLorentzVector gen_muonM_p4;
		TLorentzVector gen_pionPM_p4;

        // TTree *meson_nS_tree;
        // TLorentzVector mumu_p4, muP_p4, muM_p4;
        // UInt_t mumu_rank;

};

static const double pipm_mass = 0.1395704;
static const double meson_1SMass = 3.0969;

chicRootupler::chicRootupler(const edm::ParameterSet & iConfig): 
    chi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
    meson_nS_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("meson_nS_cand"))),
    refit1P_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1P"))),
	X_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("X_cand"))),
	refitX_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refitX"))),
    primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
    triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))), 
    isMC_(iConfig.getParameter < bool > ("isMC")),
	genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles")))
{
    edm::Service < TFileService > fs;
    chi_tree = fs->make < TTree > ("chiTree", "Tree of chic");

    chi_tree->Branch("run",                &run);
    chi_tree->Branch("event",              &event);
    chi_tree->Branch("lumiblock",          &lumiblock);
	chi_tree->Branch("numPrimaryVertices", &numPrimaryVertices);
    chi_tree->Branch("trigger",            &trigger);
	chi_tree->Branch("chi_N",              &chi_cand_N);
	chi_tree->Branch("rfchi_N",            &chi_fit_N);
	chi_tree->Branch("X_N",                &X_cand_N);
	chi_tree->Branch("rfX_N",              &X_fit_N);
	
	chi_tree->Branch("chi_pt",             &chi_pt_arr);
	chi_tree->Branch("chi_eta",            &chi_eta_arr);
	chi_tree->Branch("chi_phi",            &chi_phi_arr);
	chi_tree->Branch("chi_mass",           &chi_mass_arr);
	chi_tree->Branch("dimuon_pt",          &dimuon_pt_arr);
	chi_tree->Branch("dimuon_eta",         &dimuon_eta_arr);
	chi_tree->Branch("dimuon_phi",         &dimuon_phi_arr);
	chi_tree->Branch("dimuon_mass",        &dimuon_mass_arr);
	chi_tree->Branch("muonP_pt",           &muonP_pt_arr);
	chi_tree->Branch("muonP_eta",          &muonP_eta_arr);
	chi_tree->Branch("muonP_phi",          &muonP_phi_arr);
	chi_tree->Branch("muonP_mass",         &muonP_mass_arr);
	chi_tree->Branch("muonM_pt",           &muonM_pt_arr);
	chi_tree->Branch("muonM_eta",          &muonM_eta_arr);
	chi_tree->Branch("muonM_phi",          &muonM_phi_arr);
	chi_tree->Branch("muonM_mass",         &muonM_mass_arr);
	chi_tree->Branch("photon_pt",          &photon_pt_arr);
	chi_tree->Branch("photon_eta",         &photon_eta_arr);
	chi_tree->Branch("photon_phi",         &photon_phi_arr);
	chi_tree->Branch("photon_mass",        &photon_mass_arr);
	
	chi_tree->Branch("ele_lowerPt_pt",     &ele_lowerPt_pt_arr);
	chi_tree->Branch("ele_higherPt_pt",    &ele_higherPt_pt_arr);
	chi_tree->Branch("ctpv",               &ctpv_arr);
	chi_tree->Branch("ctpv_error",         &ctpv_error_arr);
	chi_tree->Branch("conv_vertex",        &conv_vertex_arr);
	chi_tree->Branch("chi_dz",             &chi_dz_arr);
	chi_tree->Branch("photon_flags",       &photon_flags_arr);
	chi_tree->Branch("invm1P",             &invm1P_arr);
	
	chi_tree->Branch("rf1P_chi_pt",        &rf1P_chi_pt_arr);
	chi_tree->Branch("rf1P_chi_eta",       &rf1P_chi_eta_arr);
	chi_tree->Branch("rf1P_chi_phi",       &rf1P_chi_phi_arr);
	chi_tree->Branch("rf1P_chi_mass",      &rf1P_chi_mass_arr);
    chi_tree->Branch("probFit1P",          &probFit1P_arr);
	chi_tree->Branch("rf1P_chi_index",     &rf1P_chi_index_arr);
	
	chi_tree->Branch("trk_pt",             &trk_pt_arr);
	chi_tree->Branch("trk_eta",            &trk_eta_arr);
	chi_tree->Branch("trk_phi",            &trk_phi_arr);
	chi_tree->Branch("trk_e",              &trk_e_arr);
	chi_tree->Branch("X_pt",               &X_pt_arr);
	chi_tree->Branch("X_eta",              &X_eta_arr);
	chi_tree->Branch("X_phi",              &X_phi_arr);
	chi_tree->Branch("X_mass",             &X_mass_arr);
	chi_tree->Branch("X_index",            &X_index_arr);
	
	chi_tree->Branch("rfX_pt",             &rfX_pt_arr);
	chi_tree->Branch("rfX_eta",            &rfX_eta_arr);
	chi_tree->Branch("rfX_phi",            &rfX_phi_arr);
	chi_tree->Branch("rfX_mass",           &rfX_mass_arr);
    chi_tree->Branch("probFitX",           &probFitX_arr);
	chi_tree->Branch("rfX_index",          &rfX_index_arr);

    if (isMC_) {
	   chi_tree->Branch("gen_X_p4",        &gen_X_p4);
       chi_tree->Branch("gen_chi_p4",      &gen_chi_p4);
	   chi_tree->Branch("gen_pionPM_p4",   &gen_pionPM_p4);
       chi_tree->Branch("gen_meson_nS_p4", &gen_meson_nS_p4);
       chi_tree->Branch("gen_photon_p4",   &gen_photon_p4);
       chi_tree->Branch("gen_muonP_p4",    &gen_muonP_p4);
       chi_tree->Branch("gen_muonM_p4",    &gen_muonM_p4);
    }

    // meson_nS_tree = fs->make<TTree>("JpsiTree","Tree of Jpsi");
    // meson_nS_tree->Branch("mumu_p4",  "TLorentzVector", &mumu_p4);
    // meson_nS_tree->Branch("muP_p4",   "TLorentzVector", &muP_p4);
    // meson_nS_tree->Branch("muM_p4",   "TLorentzVector", &muM_p4);
    // meson_nS_tree->Branch("trigger",  &trigger,         "trigger/i");
    // meson_nS_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
    // meson_nS_tree->Branch("mumu_rank",&mumu_rank,       "mumu_rank/i"); 
}

//Check recursively if any ancestor of particle is the given one
bool chicRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void chicRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle < pat::CompositeCandidateCollection >chi_cand_handle;
  iEvent.getByToken(chi_, chi_cand_handle);

  edm::Handle < pat::CompositeCandidateCollection >meson_nS_hand;
  iEvent.getByToken(meson_nS_, meson_nS_hand);

  edm::Handle < pat::CompositeCandidateCollection >refit1P_handle;
  iEvent.getByToken(refit1P_, refit1P_handle);
  
  edm::Handle < pat::CompositeCandidateCollection >X_handle;
  iEvent.getByToken(X_, X_handle);
  
  edm::Handle < pat::CompositeCandidateCollection >refitX_handle;
  iEvent.getByToken(refitX_, refitX_handle);

  edm::Handle < reco::VertexCollection  >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_, triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_,pruned);
  
  if (isMC_ && pruned.isValid()) {
	gen_X_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_chi_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_pionPM_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_meson_nS_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_photon_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_muonP_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	gen_muonM_p4.SetPtEtaPhiM(-999,-999,-999,-999);
	for (size_t i=0; i<pruned->size(); i++) {
		const reco::Candidate *genX = &(*pruned)[i];
		if(abs(genX->pdgId()) == 888888 && genX->status() == 2){ //X
			gen_X_p4.SetPtEtaPhiM(genX->pt(), genX->eta(), genX->phi(), genX->mass());
			for (size_t j=0; j<genX->numberOfDaughters(); j++) {
				const reco::Candidate *dau_l1 = genX->daughter(j);
				if ((dau_l1->pdgId() == 10441 || dau_l1->pdgId() == 20443 || dau_l1->pdgId() == 445) && dau_l1->status() == 2) { //chi
					gen_chi_p4.SetPtEtaPhiM(dau_l1->pt(), dau_l1->eta(), dau_l1->phi(), dau_l1->mass());
					for (size_t k=0; k<dau_l1->numberOfDaughters(); k++) {
						const reco::Candidate *dau_l2 = dau_l1->daughter(k);
						if ((dau_l2->pdgId() == 443 && dau_l2->status() == 2)) { //Jpsi
							gen_meson_nS_p4.SetPtEtaPhiM(dau_l2->pt(), dau_l2->eta(), dau_l2->phi(), dau_l2->mass());
							for (size_t l=0; l<dau_l2->numberOfDaughters(); l++) {
								const reco::Candidate *dau_l3 = dau_l2->daughter(l);
								if ((dau_l3->pdgId() == 13 && dau_l3->status() == 1)) { //muon
									gen_muonM_p4.SetPtEtaPhiM(dau_l3->pt(), dau_l3->eta(), dau_l3->phi(), dau_l3->mass());
								}
								else if((dau_l3->pdgId() == -13 && dau_l3->status() == 1)){
									gen_muonP_p4.SetPtEtaPhiM(dau_l3->pt(), dau_l3->eta(), dau_l3->phi(), dau_l3->mass());
								}
							}
						}
						else if(dau_l2->pdgId() == 22 && dau_l2->status() == 1){
							gen_photon_p4.SetPtEtaPhiM(dau_l2->pt(), dau_l2->eta(), dau_l2->phi(), dau_l2->mass());
						}
					}
				}
				else if(abs(dau_l1->pdgId()) == 211 && (dau_l1->status() == 2 || dau_l1->status() == 1)){
					gen_pionPM_p4.SetPtEtaPhiM(dau_l1->pt(), dau_l1->eta(), dau_l1->phi(), dau_l1->mass());
				}
			}
		}
		if(gen_X_p4.Pt() > 0 && gen_chi_p4.Pt() > 0 && gen_pionPM_p4.Pt() > 0 && gen_meson_nS_p4.Pt() > 0 && gen_photon_p4.Pt() > 0 && gen_muonP_p4.Pt() > 0 && gen_muonM_p4.Pt() > 0)
			break;
	}
	
	if(gen_X_p4.Pt() < -990 || gen_chi_p4.Pt() < -990 || gen_pionPM_p4.Pt() < -990 || gen_meson_nS_p4.Pt() < -990 || gen_photon_p4.Pt() < -990 || gen_muonP_p4.Pt() < -990 || gen_muonM_p4.Pt() < -990)
		std::cout << "Rootupler does not found the given decay " << run << "," << event << std::endl;
  }
	

   //grab Trigger informations
   // save it in variable trigger, 26 triggers are accessed and stored with binary bits
   // (pass 11)(pass 8)(pass 7)(pass 5)
   // es. 11 = pass 5, 7 and 11
   // es. 4 = pass only 8

   trigger = 0;
   if (triggerResults_handle.isValid()) {
	const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	unsigned int NTRIGGERS = 26;
	std::string TriggersToTest[NTRIGGERS] = {"HLT_Dimuon0_Jpsi3p5_Muon2",
										  "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
										  "HLT_Dimuon0_Jpsi_L1_NoOS",
										  "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
										  "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
										  "HLT_Dimuon0_Jpsi_NoVertexing",
										  "HLT_Dimuon0_Jpsi",
										  "HLT_Dimuon0_LowMass_L1_0er1p5R",
										  "HLT_Dimuon0_LowMass_L1_0er1p5",
										  "HLT_Dimuon0_LowMass_L1_4R",
										  "HLT_Dimuon0_LowMass_L1_4",
										  "HLT_Dimuon0_LowMass",
										  "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
										  "HLT_Dimuon25_Jpsi_noCorrL1",
										  "HLT_Dimuon25_Jpsi",
										  "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
										  "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
										  "HLT_DoubleMu4_3_Jpsi",
										  "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
										  "HLT_DoubleMu4_JpsiTrk_Displaced",
										  "HLT_DoubleMu4_Jpsi_Displaced",
										  "HLT_DoubleMu4_Jpsi_NoVertexing",
										  "HLT_Mu7p5_L2Mu2_Jpsi",
										  "HLT_Mu7p5_Track2_Jpsi",
										  "HLT_Mu7p5_Track3p5_Jpsi",
										  "HLT_Mu7p5_Track7_Jpsi"
										  };

	for (unsigned int i = 0; i < NTRIGGERS; i++) {
		for (int version = 0; version < 29; version++) {
			std::stringstream ss;
			ss << TriggersToTest[i] << "_v" << version;
			unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
				if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
				   trigger += (1<<i);
				   break;
				}
			}
		}
   }
	else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
	
	// Clear the vector for each event
	chi_pt_arr.clear();
	chi_eta_arr.clear();
	chi_phi_arr.clear();
	chi_mass_arr.clear();
	dimuon_pt_arr.clear();
	dimuon_eta_arr.clear();
	dimuon_phi_arr.clear();
	dimuon_mass_arr.clear();
	muonP_pt_arr.clear();
	muonP_eta_arr.clear();
	muonP_phi_arr.clear();
	muonP_mass_arr.clear();
	muonM_pt_arr.clear();
	muonM_eta_arr.clear();
	muonM_phi_arr.clear();
	muonM_mass_arr.clear();
	photon_pt_arr.clear();
	photon_eta_arr.clear();
	photon_phi_arr.clear();
	photon_mass_arr.clear();
	ele_higherPt_pt_arr.clear();
	ele_lowerPt_pt_arr.clear();
	ctpv_arr.clear();
	ctpv_error_arr.clear();
	conv_vertex_arr.clear();
	chi_dz_arr.clear();
	photon_flags_arr.clear();
	invm1P_arr.clear();
	rf1P_chi_pt_arr.clear();
	rf1P_chi_eta_arr.clear();
	rf1P_chi_phi_arr.clear();
	rf1P_chi_mass_arr.clear();
	probFit1P_arr.clear();
	rf1P_chi_index_arr.clear();
	trk_pt_arr.clear();
	trk_eta_arr.clear();
	trk_phi_arr.clear();
	trk_e_arr.clear();
	X_pt_arr.clear();
	X_eta_arr.clear();
	X_phi_arr.clear();
	X_mass_arr.clear();
	X_index_arr.clear();
	rfX_pt_arr.clear();
	rfX_eta_arr.clear();
	rfX_phi_arr.clear();
	rfX_mass_arr.clear();
	probFitX_arr.clear();
	rfX_index_arr.clear();
	chi_cand_N = 0;
	chi_fit_N = 0;
	X_cand_N = 0;
	X_fit_N = 0;
	bool bestCandidateOnly_ = false;
    // Accessing X_cand information (made of the prefit chi_cand and the selected track (pion)) when using chi_cand as a reference (i.e., a chi_cand corresponds to a X_cand) in a single event
    if (chi_cand_handle.isValid() && !chi_cand_handle->empty()) {

       unsigned int csize = chi_cand_handle->size();
       if (bestCandidateOnly_) csize = 1;

       for (unsigned int i = 0; i < csize; i++) {
         pat::CompositeCandidate chi_cand = chi_cand_handle->at(i);
		 
		 chi_pt_arr.push_back(chi_cand.pt());
		 chi_eta_arr.push_back(chi_cand.eta());
		 chi_phi_arr.push_back(chi_cand.phi());
		 chi_mass_arr.push_back(chi_cand.mass());
		 dimuon_pt_arr.push_back(chi_cand.daughter("dimuon")->pt());
		 dimuon_eta_arr.push_back(chi_cand.daughter("dimuon")->eta());
		 dimuon_phi_arr.push_back(chi_cand.daughter("dimuon")->phi());
		 dimuon_mass_arr.push_back(chi_cand.daughter("dimuon")->mass());
		 photon_pt_arr.push_back(chi_cand.daughter("photon")->pt());
		 photon_eta_arr.push_back(chi_cand.daughter("photon")->eta());
		 photon_phi_arr.push_back(chi_cand.daughter("photon")->phi());
		 photon_mass_arr.push_back(chi_cand.daughter("photon")->mass());
		 
		 if (chi_cand.daughter("dimuon")->daughter("muon1")->charge() < 0) {
			muonP_pt_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->pt());
			muonP_eta_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->eta());
			muonP_phi_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->phi());
			muonP_mass_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->mass());
			muonM_pt_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->pt());
			muonM_eta_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->eta());
			muonM_phi_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->phi());
			muonM_mass_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->mass());
         }

		muonP_pt_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->pt());
		muonP_eta_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->eta());
		muonP_phi_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->phi());
		muonP_mass_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon1")->mass());
		muonM_pt_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->pt());
		muonM_eta_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->eta());
		muonM_phi_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->phi());
		muonM_mass_arr.push_back(chi_cand.daughter("dimuon")->daughter("muon2")->mass());
		 
		Double_t ele_higherPt_pt = -99;
		Double_t ele_lowerPt_pt = -99;

		Double_t ele1_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track0"))->pt();
		Double_t ele2_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track1"))->pt();

		if (ele1_pt > ele2_pt) {
			ele_higherPt_pt = ele1_pt;
			ele_lowerPt_pt = ele2_pt;
		} else {
			ele_higherPt_pt = ele2_pt;
			ele_lowerPt_pt = ele1_pt;
		}

		Double_t ctpv = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlPV"); //pseudo proper decay length
		Double_t ctpv_error = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlErrPV");
		Double_t conv_vertex = chi_cand.daughter("photon")->vertex().rho();
		Double_t dz = chi_cand.userFloat("dz");
		UInt_t photon_flags = (UInt_t) dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userInt("flags");

		TLorentzVector chi_p4;
		TLorentzVector dimuon_p4;
		chi_p4.SetPtEtaPhiM(chi_cand.pt(),chi_cand.eta(),chi_cand.phi(),chi_cand.mass());
		dimuon_p4.SetPtEtaPhiM(chi_cand.daughter("dimuon")->pt(),chi_cand.daughter("dimuon")->eta(),chi_cand.daughter("dimuon")->phi(),chi_cand.daughter("dimuon")->mass());
		double QValue = chi_p4.M() - dimuon_p4.M();
		Double_t invm1P = QValue + meson_1SMass;

		ele_higherPt_pt_arr.push_back(ele_higherPt_pt);
		ele_lowerPt_pt_arr.push_back(ele_lowerPt_pt);
		ctpv_arr.push_back(ctpv);
		ctpv_error_arr.push_back(ctpv_error);
		conv_vertex_arr.push_back(conv_vertex);
		chi_dz_arr.push_back(dz);
		photon_flags_arr.push_back(photon_flags);
		invm1P_arr.push_back(invm1P);
		
		chi_cand_N = invm1P_arr.size();

		if(refit1P_handle.isValid() && !refit1P_handle->empty()){
			for(unsigned int refit_i = 0; refit_i < refit1P_handle->size(); refit_i++){
				pat::CompositeCandidate chifitted = refit1P_handle->at(refit_i);
				if(chifitted.userInt("Index") == int(i)){ //matching the corresponding fitted chi, if any
					rf1P_chi_pt_arr.push_back(chifitted.pt());
					rf1P_chi_eta_arr.push_back(chifitted.eta());
					rf1P_chi_phi_arr.push_back(chifitted.phi());
					rf1P_chi_mass_arr.push_back(chifitted.mass());
					probFit1P_arr.push_back(chifitted.userFloat("vProb"));
					rf1P_chi_index_arr.push_back(i);
				}
			}
		}
		chi_fit_N = rf1P_chi_mass_arr.size();
		// else{ //if current chi_cand did not pass the kinematic fit   //seems there is no need to push a dummy entry if current candidate failed the fit, as the index can be retrieved from the userInt
			// rf1P_chi_pt_arr.push_back(0);
			// rf1P_chi_eta_arr.push_back(0);
			// rf1P_chi_phi_arr.push_back(0);
			// rf1P_chi_mass_arr.push_back(0);
			// probFit1P_arr.push_back(0);
			// rf1P_chi_index_arr.push_back(i)
		// }
		
		if(X_handle.isValid() && !X_handle->empty()){
			for(unsigned int prefit_i = 0; prefit_i < X_handle->size(); prefit_i++){
				pat::CompositeCandidate X_cand = X_handle->at(prefit_i);
				if(X_cand.userInt("Index") == int(i)){ //matching the corresponding X cand, if any
					X_pt_arr.push_back(X_cand.pt());
					X_eta_arr.push_back(X_cand.eta());
					X_phi_arr.push_back(X_cand.phi());
					X_mass_arr.push_back(X_cand.mass());
					trk_pt_arr.push_back(X_cand.daughter("pionpm")->pt());
					trk_eta_arr.push_back(X_cand.daughter("pionpm")->eta());
					trk_phi_arr.push_back(X_cand.daughter("pionpm")->phi());
					trk_e_arr.push_back(X_cand.daughter("pionpm")->energy());
					X_index_arr.push_back(i);
				}
			}
		}
		X_cand_N = X_mass_arr.size();
		// else{
			// trk_pt_arr.push_back(0);
			// trk_eta_arr.push_back(0);
			// trk_phi_arr.push_back(0);
			// trk_e_arr.push_back(0);
			// X_pt_arr.push_back(0);
			// X_eta_arr.push_back(0);
			// X_phi_arr.push_back(0);
			// X_mass_arr.push_back(0);
			// X_index_arr.push_back(i);
		// }

		if(refitX_handle.isValid() && !refitX_handle->empty()){
			for(unsigned int refit_i = 0; refit_i < refitX_handle->size(); refit_i++){
				pat::CompositeCandidate Xfitted = refitX_handle->at(refit_i);
				if(Xfitted.userInt("Index") == int(i)){ //matching the corresponding fitted X, if any
					rfX_pt_arr.push_back(Xfitted.pt());
					rfX_eta_arr.push_back(Xfitted.eta());
					rfX_phi_arr.push_back(Xfitted.phi());
					rfX_mass_arr.push_back(Xfitted.mass());
					probFitX_arr.push_back(Xfitted.userFloat("vProb"));
					rfX_index_arr.push_back(i);
				}
			}
		}
		X_fit_N = rfX_mass_arr.size();
		// else{
			// rfX_pt_arr.push_back(0);
			// rfX_eta_arr.push_back(0);
			// rfX_phi_arr.push_back(0);
			// rfX_mass_arr.push_back(0);
			// probFitX_arr.push_back(0);
			// rfX_index_arr.push_back(i);
		// }

      }	// for i on chi_cand_handle
    } 
	chi_tree->Fill();
    
    // mumu_rank = 0;
    // if (meson_nS_hand.isValid() && !meson_nS_hand->empty()) {
      // for (unsigned int i=0; i< meson_nS_hand->size(); i++) {
        // pat::CompositeCandidate meson_nS_ = meson_nS_hand->at(i);
        // mumu_p4.SetPtEtaPhiM(meson_nS_.pt(), meson_nS_.eta(), meson_nS_.phi(), meson_nS_.mass());

        // reco::Candidate::LorentzVector vP = meson_nS_.daughter("muon1")->p4();
        // reco::Candidate::LorentzVector vM = meson_nS_.daughter("muon2")->p4();
        // if (meson_nS_.daughter("muon1")->charge() < 0) {
           // vP = meson_nS_.daughter("muon2")->p4();
           // vM = meson_nS_.daughter("muon1")->p4();
        // }

        // muP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
        // muM_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
        // meson_nS_tree->Fill();
        // mumu_rank++;
      // }
    // } 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void chicRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(chicRootupler);
