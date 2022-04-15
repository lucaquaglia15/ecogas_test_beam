#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>

using namespace std;

const char ext[10] =".root";

TF1 *gamma_flat_X = new TF1("gamma_flat_X","pol0(0)",500,3000); //Staright line to estimate the muon background in the general time profile (X plane)
TF1 *gamma_flat_Y = new TF1("gamma_flat_Y","pol0(0)",500,3000); //Staright line to estimate the muon background in the general time profile (Y plane)
TF1 *muon_gamma_X = new TF1("muon_gamma_X","pol0(0)+gaus(1)",0,5000); //Gaussian + constant (X plane)
TF1 *muon_gamma_Y = new TF1("muon_gamma_Y","pol0(0)+gaus(1)",0,5000); //Gaussian + constant (Y plane)
TF1 *muon_gamma[2];
TF1 *gamma_flat[2];
TF1 *fit_eff; 
// fit_eff = new TF1("fit_eff","[0]/(1+TMath::Exp(-[1]*(x-[2])))",8800,10600); //Sigmoid for eff curve fit

void twodeff(const int scan) {

	//Histograms cosmetics
	gStyle->SetHistFillColor(kBlue);
	gStyle->SetHistFillStyle(3357);
	gStyle->SetHistLineColor(kBlue);

	//Show fit results on the plots
	gStyle->SetOptFit(1111);

	//Strip mapping for X and Y plane 
	ifstream mappingX;
	ifstream mappingY;
	mappingX.open("mappingX.txt");
	mappingY.open("mappingY.txt");

	//X mapping 
	int chX, strpX, muteX;
	vector <int> tdc_chX, rpc_strpX;

	while (mappingX >> muteX >> chX >> strpX) {

		if (muteX == 1) {
			tdc_chX.push_back(chX);
			rpc_strpX.push_back(strpX);
		}

		else continue;
	}

	//Y mapping
	int chY, strpY, muteY;
	vector <int> tdc_chY, rpc_strpY;

	while (mappingY >> muteY >> chY >> strpY) {

		if (muteY == 1) {
			tdc_chY.push_back(chY);
			rpc_strpY.push_back(strpY);
		}

		else continue;
	}

	//First and last TDC channels for X and Y plane
	int instripX = tdc_chX.front();
	int finstripX = tdc_chX.back();
	int instripY = tdc_chY.front();
	int finstripY = tdc_chY.back();

	//X strip mapping
	map <int, int> strp_mapX;
	for (unsigned int i = 0; i < tdc_chX.size();i++) {
		strp_mapX.insert(pair<int, int>(tdc_chX.at(i),rpc_strpX.at(i)));
	}

	//Y strip mapping
	map <int, int> strp_mapY;
	for (unsigned int i = 0; i < tdc_chY.size();i++) {
		strp_mapY.insert(pair<int, int>(tdc_chY.at(i),rpc_strpY.at(i)));
	}	//Enter in the scan folder 

	ifstream ftimeX;
	ftimeX.open(("Timestamps_X_scan"+to_string(scan)+".txt").c_str());
	double start_muonX, end_muonX, meanX, muon_windowX;
	vector <double> vstartX, vendX, vmeanX, vwindowX;

	while (ftimeX >> start_muonX >> meanX >> end_muonX >> muon_windowX) {
		vstartX.push_back(start_muonX);
		vmeanX.push_back(meanX);
		vendX.push_back(end_muonX);
		vwindowX.push_back(muon_windowX);
	}

	ifstream ftimeY;
	ftimeY.open(("Timestamps_Y_scan"+to_string(scan)+".txt").c_str());
	double start_muonY, end_muonY, meanY, muon_windowY;
	vector <double> vstartY, vendY, vmeanY, vwindowY;

	while (ftimeY >> start_muonY >> meanY >> end_muonY >> muon_windowY) {
		vstartY.push_back(start_muonY);
		vmeanY.push_back(meanY);
		vendY.push_back(end_muonY);
		vwindowY.push_back(muon_windowY);
	}

	string fol = "/media/luca/Elements/TB_October_2021/STD/Scan_00" + to_string(scan);

	//Create .root file with all the output data
	TFile *fout = new TFile(("outfile_ALICE_2D_"+to_string(scan)+".root").c_str(),"RECREATE");

	fout->cd();
	TDirectory *cdtof[12]; //12 directories in the root file, one for each HV point

	gSystem->cd(fol.c_str()); //Enter folder

	//////////////////////
	//					//
	//	CANVAS BOOKING	//
	//					//
	//////////////////////

	TCanvas *c_time_profile_ALICEX = new TCanvas(); //Time profile X
	TCanvas *c_time_profile_ALICEY = new TCanvas(); //Time profile Y
	TCanvas *c_strp_profile_ALICEX = new TCanvas(); //Strip profile of all hits X
	TCanvas *c_strp_profile_ALICEY = new TCanvas(); //Strip profile of all hits Y

	/////////////////////
	//				   //
	//HISTOGRAM BOOKING//
	//				   //
	/////////////////////

	//Define histograms in arrays of 12 elements, one per HV point
	TH1F *time_profile_ALICEX[12];
	TH1F *time_profile_ALICEY[12];
	TH1F *strp_profile_ALICEX[12];
	TH1F *strp_profile_ALICEY[12];

	//HV, current, error on HV, error on current, current density and error on current density + env parameters
	vector<double> veff, imon, eveff, eimon, currdens, ecurrdens,vmon, evmon, veff_eco; 
	vector<double> vT, vp;
	double p0 = 990., T0 = 293., T, p, tens_mon, tens_eff; //mbar, K

	//Muon efficiency and error, muon efficiency with gamma correction and error, fake efficiency and error
	vector<double> eff, e_eff;

	vector<float> muon_timesX, muon_hitsX, muon_timesY, muon_hitsY;

	//Count how many .root files are in the folder
	const char* entry; 
	const char* filename;
    int file_count = 0;
    TString str;
      
    char* dir = gSystem->ExpandPathName(fol.c_str());
    void* dirp = gSystem->OpenDirectory(dir);
    
    while((entry=gSystem->GetDirEntry(dirp))) { //if the file ends with ".root" -> increase the file count by 1
    	str = entry;
      	if(str.EndsWith(ext)) {
      		file_count++;
    	}
    }

    cout << "Number of HV points: " << file_count/3 << endl;

	//Read DAQ root files: cycle on all 12 HV points --- MAIN part of the code
	for (int point = 0; point < file_count/3; point++) {

		string cartella = "HV" + to_string(point+1); //Create folders for different HV points in the out root file
		cdtof[point] = fout->mkdir(cartella.c_str()); 

		//Create histograms + set axes title
		time_profile_ALICEX[point] = new TH1F(("a"+cartella).c_str(),("Time profile ALICE X"+cartella).c_str(),5000,0,5000);
		time_profile_ALICEX[point]->GetXaxis()->SetTitle("Time[ns]");
		time_profile_ALICEX[point]->GetYaxis()->SetTitle("Hits"); 
		time_profile_ALICEY[point] = new TH1F(("b"+cartella).c_str(),("Time profile ALICE Y"+cartella).c_str(),5000,0,5000);
		time_profile_ALICEY[point]->GetXaxis()->SetTitle("Time[ns]");
		time_profile_ALICEY[point]->GetYaxis()->SetTitle("Hits"); 
		strp_profile_ALICEX[point] = new TH1F(("c"+cartella).c_str(),("Strip profile ALICE X"+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICEX[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICEX[point]->GetYaxis()->SetTitle("Hits");
		strp_profile_ALICEY[point] = new TH1F(("d"+cartella).c_str(),("Strip profile ALICE Y"+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICEY[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICEY[point]->GetYaxis()->SetTitle("Hits");

		//Open the different DAQ .root files
		string hvpoint = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_DAQ.root"; //For data
		string caen = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_CAEN.root"; //For HV and I
		string dip = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_DIP.root"; //For temperature and pressure

		//Define the file to be opened
		TFile *f = new TFile(hvpoint.c_str(),"READ");

		//Open the tree in the file
		TTree *t = (TTree*)f->Get("RAWData");

		//Define the elements in the tree
		int nev, eventry, quality, trig_type, trigger = 0, count_Y = 0, count_X = 0, count_both = 0; //trig_type = 1 -> during spill, trigger_type = 0 -> outside spill
		bool hit_muonX = 0, hit_muonY = 0, hit_gamma = 0; //To count if the strips have been fired in an event
		vector <int> *channel = 0; //TDC channels 
		vector <float> *timestamp = 0; //TDC timestamps

		//Set branch addresses
		t->SetBranchAddress("EventNumber",&nev);
		t->SetBranchAddress("number_of_hits", &eventry);
		t->SetBranchAddress("Quality_flag",&quality);
		t->SetBranchAddress("TDC_channel", &channel);
		t->SetBranchAddress("TDC_TimeStamp",&timestamp);
		t->SetBranchAddress("TriggerTag",&trig_type);

		int nentries = t->GetEntries(); //Number of entries in the TTree (=number of triggers for given HV point)

		//Cycle on all the tree entries to take the data from the tree
		for (int i = 0; i < nentries; i++) { 
			t->GetEntry(i);
			count_Y = 0;
			count_X = 0;
			//cout << "event number " << i+1 << endl; //Debug
			if (trig_type == 1) trigger++; //count the number of muon triggers

			for (unsigned int k = 0; k < channel->size(); k++) { //Go through the vector that contains all the channels of the event

				//All hits profile
				if (channel->at(k) <= finstripY && channel->at(k) >= instripY && trig_type == 1 && timestamp->at(k) > 100) { //Y plane
					strp_profile_ALICEY[point]->Fill(strp_mapY.lower_bound(channel->at(k))->second); //Strip profile of all hits Y
					time_profile_ALICEY[point]->Fill(timestamp->at(k));//Time profile of the chamber Y
					count_Y++;
					//cout << "strips in the event (Y): " << strp_mapY.lower_bound(channel->at(k))->second << " "; //Debug
				}

				else if (channel->at(k) <= finstripX && channel->at(k) >= instripX && trig_type == 1 && timestamp->at(k) > 100) { //X plane
					strp_profile_ALICEX[point]->Fill(strp_mapX.lower_bound(channel->at(k))->second); //Strip profile of all hits Y
					time_profile_ALICEX[point]->Fill(timestamp->at(k));//Time profile of the chamber Y
					count_X++;
					//cout << "strips in the event (X): " << strp_mapX.lower_bound(channel->at(k))->second << " "; //Debug
				}
			} //end cycle on channel vector			
			/*
			//Debug printouts
			cout << "size timestamp: " << timestamp->size() << ", size timestamp: " << timestamp->size() << endl;
			cout << "Event number: " << nev << ". Channels fired: " << ", quality flag: " << quality << endl;
			for (unsigned int j = 0; j < channel->size(); j++) {
				cout << channel->at(j) << " at: " << timestamp->at(j) << "-";
			}
			cout << endl << endl << endl;
			*/
		} //End of cycle on all tree entries

		//Draw the tim profile with more complicated gaussian + plo0 fit to extract the signal window
		//Y-plane
		/*muon_gamma[0] = new TF1("muon_gammaY","pol0(0)+gaus(1)",0,5000);
		gamma_flat[0] = new TF1("gamma_flat_Y","pol0(0)",500,3000);
		c_time_profile_ALICEY->cd();
		//muon_gamma_Y->SetParLimits(2,4500.,4600.);
		muon_gamma_Y->SetParLimits(2,4460.,4540.);
		muon_gamma_Y->SetParLimits(3,0.,40.);
		time_profile_ALICEY[point]->Fit("gamma_flat_Y","RM+"); //Horizontal line fit to gamma background
		double avg_gamma_Y = gamma_flat_Y->GetParameter(0); //Get the value of the fit
		muon_gamma_Y->SetParameter(0,avg_gamma_Y); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		//muon_gamma_Y->FixParameter(0,avg_gamma_Y); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
		time_profile_ALICEY[point]->Fit("muon_gamma_Y","RM+"); //Fit with gaussian + constant
		time_profile_ALICEY[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICEY[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICEY[point]->Draw("HISTO");
		//gamma_flat_Y->Draw("SAME");
		muon_gamma_Y->Draw("SAME");

		//muon_gamma[0]->SetParLimits(2,4500.,4600.);
		//muon_gamma[0]->SetParLimits(1,0.,4540.);
		muon_gamma[0]->SetParLimits(2,4460.,4540.);
		muon_gamma[0]->SetParLimits(3,0.,40.);
		time_profile_ALICEY[point]->Fit(gamma_flat[0],"RM+"); //Horizontal line fit to gamma background
		double avg_gamma_Y = gamma_flat[0]->GetParameter(0); //Get the value of the fit
		muon_gamma[0]->SetParameter(0,avg_gamma_Y); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		//muon_gamma[0]->FixParameter(0,avg_gamma_Y); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
		time_profile_ALICEY[point]->Fit(muon_gamma[0],"RM+"); //Fit with gaussian + constant
		time_profile_ALICEY[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICEY[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICEY[point]->Draw("HISTO");
		//gamma_flat[0]->Draw("SAME");
		muon_gamma[0]->Draw("SAME");

		//X-plane
		muon_gamma[1] = new TF1("muon_gammaX","pol0(0)+gaus(1)",0,5000);
		gamma_flat[1] = new TF1("gamma_flat_X","pol0(0)",500,3000);
		c_time_profile_ALICEX->cd();
		//muon_gamma_X->SetParLimits(2,4500.,4600.);
		muon_gamma_X->SetParLimits(2,4460.,4540.);
		muon_gamma_X->SetParLimits(3,0.,40.);
		time_profile_ALICEX[point]->Fit("gamma_flat_X","RM+"); //Horizontal line fit to gamma background
		double avg_gamma_X = gamma_flat_X->GetParameter(0); //Get the value of the fit
		muon_gamma_X->SetParameter(0,avg_gamma_X); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		//muon_gamma_X->FixParameter(0,avg_gamma_X); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
		time_profile_ALICEX[point]->Fit("muon_gamma_X","RM+"); //Fit with gaussian + constant
		time_profile_ALICEX[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICEX[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICEX[point]->Draw("HISTO");
		//gamma_flat_X->Draw("SAME");
		muon_gamma_X->Draw("SAME");

		//muon_gamma[1]->SetParLimits(2,4500.,4600.);
		muon_gamma[1]->SetParLimits(1,0.,4540.);
		muon_gamma[1]->SetParLimits(2,4460.,4540.);
		muon_gamma[1]->SetParLimits(3,0.,40.);
		time_profile_ALICEX[point]->Fit(gamma_flat[1],"RM+"); //Horizontal line fit to gamma background
		double avg_gamma_X = gamma_flat[1]->GetParameter(0); //Get the value of the fit
		muon_gamma[1]->SetParameter(2,avg_gamma_X); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		//muon_gamma_X->FixParameter(2,avg_gamma_X); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
		time_profile_ALICEX[point]->Fit(muon_gamma[1],"RM+"); //Fit with gaussian + constant
		time_profile_ALICEX[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICEX[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICEX[point]->Draw("HISTO");
		//gamma_flat[1]->Draw("SAME");
		muon_gamma[1]->Draw("SAME");*/

		//Start and end time of the muon window
		//X-plane
		/*double start_muonX = muon_gamma_X->GetParameter(2) - 3*muon_gamma_X->GetParameter(3); 
		double end_muonX = muon_gamma_X->GetParameter(2) + 3*muon_gamma_X->GetParameter(3);
		double meanX = muon_gamma_X->GetParameter(2);
		//Y-plane
		double start_muonY = muon_gamma_Y->GetParameter(2) - 3*muon_gamma_Y->GetParameter(3); 
		double end_muonY = muon_gamma_Y->GetParameter(2) + 3*muon_gamma_Y->GetParameter(3);
		double meanY = muon_gamma_Y->GetParameter(2);*/

		//X-plane
		/*double start_muonX = muon_gamma[1]->GetParameter(2) - 3*muon_gamma[1]->GetParameter(3); 
		double end_muonX = muon_gamma[1]->GetParameter(2) + 3*muon_gamma[1]->GetParameter(3);
		double meanX = muon_gamma[1]->GetParameter(2);
		//Y-plane
		double start_muonY = muon_gamma[0]->GetParameter(2) - 3*muon_gamma[0]->GetParameter(3); 
		double end_muonY = muon_gamma[0]->GetParameter(2) + 3*muon_gamma[0]->GetParameter(3);
		double meanY = muon_gamma[0]->GetParameter(2);

		double muon_windowX = (end_muonX - start_muonX)*1e-9; //Muon window in s (X-plane)
		double muon_windowY = (end_muonY - start_muonY)*1e-9; //Muon window in s (Y-plane)
		cout << endl << "muon window X: " << start_muonX << " ns < " << meanX << " ns < " << end_muonX << " ns" << endl; //X-plane
		cout << endl << "muon window Y: " << start_muonY << " ns < " << meanY << " ns < " << end_muonY << " ns" << endl; //Y-plane*/

		//Draw strip profile, all hits
		c_strp_profile_ALICEX->cd();
		strp_profile_ALICEX[point]->Draw("HISTO");
		c_strp_profile_ALICEY->cd();
		strp_profile_ALICEY[point]->Draw("HISTO");

		for (int i = 0; i < nentries; i++) { 
			muon_timesY.clear();
			muon_hitsY.clear();
			muon_timesX.clear();
			muon_hitsX.clear();

			t->GetEntry(i);

			for (unsigned int k = 0; k < channel->size(); k++) { //Go through the vector that contains all the channels of the event

				//Muon hits profile Y plane
				if ((channel->at(k) <= finstripY && channel->at(k) >= instripY && trig_type == 1) && (timestamp->at(k) <= vendY.at(point) && timestamp->at(k) >= vstartY.at(point))) { //Y plane
					//strp_profile_ALICE_muonsY[point]->Fill(strp_mapY.lower_bound(channel->at(k))->second); //Strip profile of muons
					hit_muonY = 1; //The chamber saw something, set the counter to true, it wil be reset at the end of the event
					//muon_counterY++;
					//cout << "K: " << k << ", time:" << timestamp->at(k) << ", channel: " << channel->at(k) << " "; //Debug for cluster size
					muon_timesY.push_back(timestamp->at(k));
					muon_hitsY.push_back(channel->at(k));
				}
				
				//Muon hits profile X plane
				if ((channel->at(k) <= finstripX && channel->at(k) >= instripX && trig_type == 1) && (timestamp->at(k) <= vendX.at(point) && timestamp->at(k) >= vstartX.at(point))) { //X plane
					//strp_profile_ALICE_muonsX[point]->Fill(strp_mapX.lower_bound(channel->at(k))->second); //Strip profile of muons
					hit_muonX = 1; //The chamber saw something, set the counter to true, it wil be reset at the end of the event
					//muon_counterX++;
					//cout << "K: " << k << ", time:" << timestamp->at(k) << ", channel: " << channel->at(k) << " "; //Debug for cluster size
					muon_timesX.push_back(timestamp->at(k));
					muon_hitsX.push_back(channel->at(k));	
				}
			}

			if (hit_muonY == 1) count_Y++; //Increase the number of events seen by the RPC (Y plane)
			if (hit_muonX == 1) count_X++; //Increase the number of events seen by the RPC (Y plane)
			if (hit_muonX == 1 && hit_muonY == 1) count_both++;
			hit_muonY = 0; //reset the bool for muon events counting
			hit_muonX = 0;
		} //End of cycle on tree entries

		//2D muon efficiency without gamma contribution
		double raw2Deff = (count_both/(double)trigger);
		double eraw2Deff = TMath::Sqrt((raw2Deff*(1-raw2Deff))/trigger);

		//Push back in vectors for plots
		eff.push_back(raw2Deff*100);
		e_eff.push_back(eraw2Deff*100);

		////////////////////
		//				  //
		//	EFFECTIVE HV  //
		//				  //
		////////////////////

		//Get current and HV values
		TFile *f1 = new TFile(caen.c_str(),"READ");
		TH1F *HVeff = (TH1F*)f1->Get("HVeff_ALICE-2-0");
		TH1F *HVmon = (TH1F*)f1->Get("HVmon_ALICE-2-0");
		TH1F *Imon = (TH1F*)f1->Get("Imon_ALICE-2-0");

		tens_mon = HVmon->GetMean();
		veff.push_back(HVeff->GetMean());
		vmon.push_back(tens_mon);

		if (HVeff->GetStdDev() > 1) eveff.push_back(HVeff->GetStdDev());
		else eveff.push_back(1.);
		if (HVmon->GetStdDev() > 1) evmon.push_back(HVmon->GetStdDev());
		else evmon.push_back(1.);
		imon.push_back(Imon->GetMean());
		if (Imon->GetStdDev() > 0.1) eimon.push_back(Imon->GetStdDev());
		else eimon.push_back(0.1);

		TFile *f2 = new TFile(dip.c_str(),"READ");
		TH1F *T_in = (TH1F*)f2->Get("TIN");
		TH1F *p_in = (TH1F*)f2->Get("P");

		T = T_in->GetMean() + 273.15;
		p = p_in->GetMean();

		vT.push_back(T);
		vp.push_back(p);

		tens_eff = tens_mon*(T/T0)*(p0/p); //Effective HV calculated with ecogas definition

		veff_eco.push_back(tens_eff);

		fout->cd(); //open .root file
		cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
		c_time_profile_ALICEX->Write(("Time_profile_ALICE_X_"+cartella).c_str()); //Save time profile X
		c_time_profile_ALICEY->Write(("Time_profile_ALICE_Y_"+cartella).c_str()); //Save time profile Y
		c_strp_profile_ALICEX->Write(("Strip_profile_ALICE_all_hits_X_"+cartella).c_str()); //Save strip profile (all hits X)
		c_strp_profile_ALICEY->Write(("Strip_profile_ALICE_all_hits_Y_"+cartella).c_str()); //Save strip profile (all hits Y)
	} //End of cycle on HV points

	//Set parameters of the functions used for the fit of eff curve
	fit_eff = new TF1("fit_eff","[0]/(1+TMath::Exp(-[1]*(x-[2])))",veff_eco.front(),veff_eco.back()); //Sigmoid for eff curve fit
	fit_eff->SetParLimits(0,45.,110.);
	fit_eff->SetParLimits(2,9000.,10000.);

	//Efficiency with gamma background contribution
	TCanvas *ceff = new TCanvas();
	ceff->cd();

	TGraphErrors *Raw2DEff = new TGraphErrors(veff_eco.size(),&veff_eco[0],&eff[0],&eveff[0],&e_eff[0]);
	Raw2DEff->SetTitle(("2D eff curve for scan " + to_string(scan)).c_str());
	Raw2DEff->GetXaxis()->SetTitle("HV_{eff} [V]");
	Raw2DEff->GetYaxis()->SetTitle("#epsilon");
	Raw2DEff->GetYaxis()->SetRangeUser(veff_eco.at(0),veff_eco.at(veff_eco.size()-1));
	Raw2DEff->GetYaxis()->SetRangeUser(0.,100.);
	Raw2DEff->SetMarkerStyle(8);
	Raw2DEff->SetMarkerSize(1);
	Raw2DEff->Fit(fit_eff,"RM+");
	Raw2DEff->Draw("AP");
	fit_eff->Draw("SAME");

	//Fit parameters
	double effmax = fit_eff->GetParameter(0);
	double lambda = fit_eff->GetParameter(1);
	double HV50 = fit_eff->GetParameter(2);

	//double WP = (TMath::Log(19)/lambda)+HV50+150; //CMS definition
	double WP = fit_eff->GetX(0.95*effmax,veff_eco.front(),veff_eco.back()); //Definition for ecogas paper = 95% of maximum efficiency
	double effwp = fit_eff->Eval(WP);

	TPaveText *teff = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	teff->AddText(Form("Eff max: %3.1lf",effmax));
	teff->AddText(Form("#lambda: %3.4lf",lambda));
	teff->AddText(Form("HV 50%% (V): %3.1lf",HV50));
	teff->AddText(Form("WP (V): %3.1lf",WP)); 
	teff->AddText(Form("Eff WP: %3.1lf",effwp));
	teff->Draw("SAME");

	fout->cd();
	ceff->Write(("Eff_curve_2D_" + to_string(scan)).c_str());

	fout->Close();

}