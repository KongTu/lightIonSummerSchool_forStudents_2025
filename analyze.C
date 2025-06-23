/*
To start, download the pseudo-data:
https://drive.google.com/file/d/1ofgtVJOQZuD-fsEqRfIrLZhoGxKIj9dp/view?usp=sharing

Goal of this task:

- process all the data that has a mixture of coherent and incoherent VM production. 
- remove all the incoherent production as much as possible to look at coherent production
- reconstruct the t distribution and perform imagining. 
- the final questions are: 
	what ion are we looking at? 
	what is the radius?
  (the answers have to be based on reconstructed "pseudo-data")

Remarks:
-1 the MC particles are provided. Make use of that. The events that are generated 
	are Q2>1 GeV2

-2 "particles" branch has all charged particles in the main detector and B0;
   EEMC cluster is backward EMCal for scattered electron
   RP, OMD, ZDC are three far forward detectors with hit or cluster information.

-3 daughter of VM decay is kaon+ and kaon-.

-4 the best resolution of scattered electron is by using:
  	1. the EEMC energy cluster (e' energy), 
  	2. the eta, phi from tracks ("partciles")
  	3. and the known electron mass.

*/
#include "analyze.h"

void analyze(){

	//input file
	TFile* file = new TFile("tree_all.root");

	//read tree
	auto tree = (TTree*) file->Get("miniTree");
	Event* event = nullptr;
	tree->SetBranchAddress("event", &event);
	
	//number of events
	int nEvents=tree->GetEntries();
	
	//beam momentum with 25 mrad crossing angle.
	TLorentzVector eIn(0,0,-18,18);
	TLorentzVector aIn(-3.438,0,137.457,137.503);

	TFile* output = new TFile("results.root","RECREATE");
	//example hist
	TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",100,0,20);
	// TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,20);
	
	// TH1D* h_EoverP_REC = new TH1D("h_EoverP_REC",";E/p",200,0,2);
	// TH1D* h_Epz_REC = new TH1D("h_Epz_REC",";E - p_{z} (GeV)",200,0,50);

	// TH1D* h_massMC = new TH1D("h_massMC",";mass",100,0.9,1.2);
	// TH1D* h_massREC = new TH1D("h_massREC",";mass",100,0.9,1.2);

	TH1D* h_tMC = new TH1D("h_tMC",";t_{MC} (GeV^{2})",100,0,0.3);
	TH1D* h_tREC = new TH1D("h_tREC",";t_{REC} (GeV^{2})",100,0,0.3);

	//event loop
	for(int ievt=0;ievt<nEvents;ievt++){
		tree->GetEntry(ievt);
		
		// don't forget to remove incoherent in the MC, if you want to look at MC.

		//mc particle loop
	
			TLorentzVector escat(0,0,0,0);
			TLorentzVector Kaon(0,0,0,0);
			vector<TLorentzVector> kdaug;

		bool incoherent_MC=false; // if particles are found in forward detectors, use "incoherent_MC" to veto. 

		for(int i=0;i<event->mcp.size();i++){
			const MCp& p = event->mcp[i];
			int status=p.status;
			float mass=p.mass;
				//the tree for mc particles has px,py,pz,mass, and status.
			if(status != 1) continue; // status == 1 means stable final-state particles.
			
			//do your MC analysis here.//
				// find the scattered electron and the two kaon daughters
				// check the invariant mass, and calculate t.
			// hint:
				// in exclusive reactions, only e', K+ and K- are in the main detector. 
				// the rest are the breakup product from ions, mostly in the forward detectors.
				
		}


   	/*phase space*/	
		TLorentzVector qbeam=eIn-escat;
		double Q2=-(qbeam).Mag2();  
		double pq=aIn.Dot(qbeam);
		double y= pq/aIn.Dot(eIn);

	    //MC level phase space cut; 
		if(Q2<2.||Q2>10.) continue;
		if(y<0.01||y>0.85) continue;

		h_Q2_e->Fill(Q2);
		/*end phase space*/

			//try to get t_MC by using the MC e' and VM. with 
			// a chosen t reconstructed method. provided in the analyze.h

		double t_MC=999.;
			if(!incoherent_MC) h_tMC->Fill(t_MC); //only save coherent.



		/********separation of MC (above) and REC (below) analysis*********/

		bool incoherent_REC = false;
		//define the condition for incoherent at the REC level.
			//hint: use forward detector information

						// example macro for accessing FF detectors.

						// //reconstructed ZDC loop
						// for(int i=0;i<event->clusters_zdc.size();i++){
						// 	const Cluster_ZDC& clus = event->clusters_zdc[i];
						// }

						// //reconstructed RP hit loop
						// for(int i=0;i<event->hit_rp.size();i++){
						// 	const Hit_RP& hit = event->hit_rp[i];
						// }

						// //reconstructed OMD hit loop
						// for(int i=0;i<event->hit_omd.size();i++){
						// 	const Hit_OMD& hit = event->hit_omd[i];
						// }

		if(incoherent_REC) continue;

		//reco events:
		TLorentzVector escat_REC(0,0,0,0);
		TLorentzVector K1_REC(0,0,0,0);
		TLorentzVector K2_REC(0,0,0,0);

		//reconstructed EEMC loop (why do you need this? look at Remark #4)
		for(int i=0;i<event->clusters_eemc.size();i++){
			const Cluster_EEMC& clus = event->clusters_eemc[i];
		
			//hint: find the leading energy cluster for scattered electron's energy
				//   how do we know its an electron? the technique we usually use (without a PID detector) is called E/p
				//   where the energy divide by momentum should be peaked at 1 due to the small mass of electron. 
				//   momentum p comes from tracking
		}

		double maxpz=0.;
		TLorentzVector elec_trk(0,0,0,0); //this is track only scattered electron for E/p later.
		TLorentzVector hfs_e(0,0,0,0); // this has all charged particles's four vectors
		TLorentzVector part(0,0,0,0);

		//reconstructed charged particle loop
		for(int i=0;i<event->particles.size();i++){
			const Particle& p = event->particles[i];
			TVector3 trk(p.px,p.py,p.pz);
			if(fabs(p.pz)>maxpz && trk.Eta()<-1.5) {
				maxpz=p.pz;
				
				// get the best scattered electron //

					//escat_REC.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
					//elec_trk.SetVectM(trk,MASS_ELECTRON);
			}

			// daughter of vm.
			if(fabs(trk.Eta())<1.5){
				//get the daughter Kaon vector; limit our acceptance in +- 1.5 to be easier.
			
				//K1_REC.SetPxPyPzE();
				//K2_REC.SetPxPyPzE();

			}
			//hfs + scattered electron
			if(fabs(trk.Eta())<5) {
				part.SetVectM(trk,MASS_PION);
				hfs_e += part;
			}
		}

		//selections
		//don't forget removing incoherent.

		//E over p
		double EoverP=escat_REC.E()/elec_trk.P();

		//E - pz (use the better scat elec four vector)
		// this is another trick to help select the DIS scattered electron 
		// one should work out the math, where E-pz should peak at 2xE(incoming electron), 
		// which is 18x2 = 36 GeV.
		double EpzREC= (hfs_e-elec_trk+escat_REC).E() - (hfs_e-elec_trk+escat_REC).Pz();

		TLorentzVector qbeamREC=eIn-escat_REC;
		double Q2REC=-(qbeamREC).Mag2();  
		double pqREC=aIn.Dot(qbeamREC);
		double yREC= pqREC/aIn.Dot(eIn);

		//Event selection:
		if( EpzREC<32||EpzREC>38 ) continue;
		if( EoverP<0.8||EoverP>1.16 ) continue;		

		//REC level phase space cut
		if(Q2REC<2.||Q2REC>10.) continue;
		if(yREC<0.01||yREC>0.85) continue;
		
		double mass_REC=(K1_REC+K2_REC).M();
			//hint: check mass if the same as MC.
		
		double t_REC = 999.; //find out what t_REC is.

		// if( )// if you want to cut on the selection
		// {
			double t_L=999.;// one can use method L or any other method. method L is in analyze.h.
			h_tREC->Fill(t_L);
		// }

	}

	output->Write();
	output->Close();

	

}