#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
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
#include "Fit/FitResult.h"
#include "TBox.h"
#include "TLine.h"

using namespace std;

const char ext[10] =".root";

TF1 *muon_sig = new TF1("muon_sig","gaus(0)",4450,4650); //Gaussian fit to estimate the muon window in the RPC time profile with variable x-range in order to test the std devition values while changing the range
TF1 *gamma_flat = new TF1("gamma_flat","pol0(0)",500,3000); //Staright line to estimate the muon background in the general time profile
TF1 *muon_gamma = new TF1("muon_gamma","pol0(0)+gaus(1)",0,5000); //Gaussian + constant (value of the constant determined with the linear fit of the background)
TF1 *gamma_flat_temp = new TF1("gamma_flat_temp","pol0(0)",500,3000); //Staright line to estimate the muon background in the general time profile
TF1 *muon_gamma_temp = new TF1("muon_gamma_temp","pol0(0)+gaus(1)",0,5000); //Gaussian + constant (value of the constant determined with the linear fit of the background)

TF1 *muon_strip = new TF1("muon_strip","gaus(0)",0,16); //Gaussian fit for beam profile in strip profile
TF1 *fit_eff = new TF1("fit_eff","[0]/(1+TMath::Exp(-[1]*(x-[2])))",8800,10600); //Sigmoid for eff curve fit with gamma background
TF1 *fit_eff_gamma = new TF1("fit_eff_gamma","[0]/(1+TMath::Exp(-[1]*(x-[2])))",8800,10600); //Sigmoid for eff curve fit without gamma background

void time_profile_analysis(const int scan, const double clus_time) {
	
	//string plane = "X";
	string plane = "Y";

	//Histograms cosmetics
	gStyle->SetHistFillColor(kBlue);
	gStyle->SetHistFillStyle(3357);
	gStyle->SetHistLineColor(kBlue);

	//Show fit results on the plots
	gStyle->SetOptFit(1111);

	//Strip mapping
	ifstream mapping;
	if (plane == "X") mapping.open("mappingX.txt");
	else if (plane == "Y") mapping.open("mappingY.txt");
	
	int ch, strp, mute;
	vector <int> tdc_ch, rpc_strp;

	while (mapping >> mute >> ch >> strp) {
		if (mute == 1) {
			tdc_ch.push_back(ch);
			rpc_strp.push_back(strp);
		}
		else continue;
	}

	int instrip = tdc_ch.front();
	int finstrip = tdc_ch.back();

	//cout << instrip << "\t" << finstrip << endl; //Debug

	map <int, int> strp_map;
	for (unsigned int i = 0; i < tdc_ch.size();i++) {
		strp_map.insert(pair<int, int>(tdc_ch.at(i),rpc_strp.at(i)));
	}

	//RPC parameters to calculate rates
	double strip_area = 150; //cm2, area of a strip
	double meas_time = 5000; //ns, time of a single measurement
	double area = 2500; //Area of the ALICE RPC

	ifstream ftimesX;
	ifstream ftimesY;
	if (plane == "X") {
		ftimesX.open(("Timestamps_X_scan"+to_string(scan)+".txt").c_str());
	}
	else if (plane == "Y") {
		ftimesY.open(("Timestamps_Y_scan"+to_string(scan)+".txt").c_str());	
	} 

	vector<double> mu_window_start, mu_window_end, mu_window_mean, mu_window; 
	double muS, muM, muE, muW; 

	if (plane == "X") {
		while (ftimesX >> muS >> muM >> muE >> muW) {
			mu_window_start.push_back(muS);
			mu_window_mean.push_back(muM);
			mu_window_end.push_back(muE);
			mu_window.push_back(muW);
		}
	}
	else if (plane == "Y") {
		while (ftimesY >> muS >> muM >> muE >> muW) {
			mu_window_start.push_back(muS);
			mu_window_mean.push_back(muM);
			mu_window_end.push_back(muE);
			mu_window.push_back(muW);
		}
	} 

	for (unsigned int i = 0; i < mu_window_end.size(); i++) cout << mu_window_start.at(i) << endl;

	vector<double> strips{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

	//Enter in the scan folder 
	string fol = "/media/luca/Elements/TB_October_2021/STD/Scan_00" + to_string(scan);

	//Create .root file with all the output data
	TFile *fout = new TFile(("Time_profile_analysis"+plane+to_string(scan)+".root").c_str(),"RECREATE");

	fout->cd();
	TDirectory *cdtof[12]; //12 directories in the root file, one for each HV point
	TDirectory *cstrip[16]; //16 subdirectories for strips

	gSystem->cd(fol.c_str()); //Enter folder

 	//////////////////////
	//					//
	//	CANVAS BOOKING	//
	//					//
	//////////////////////

	TCanvas *c_time_profile_ALICE = new TCanvas(); //Time profile
	TCanvas *c_time_profile_analysis_ALICE = new TCanvas(); //Time profile
	TCanvas *c_strp_profile_ALICE = new TCanvas(); //Strip profile of all hits
	TCanvas *c_strp_profile_ALICE_muons = new TCanvas(); //Strip profile muons
	TCanvas *c_strp_profile_ALICE_muons_no8_no16 = new TCanvas(); //Strip profile muons if cs is not 16 or 8
	TCanvas *c_strp_profile_ALICE_gamma = new TCanvas(); //Strip profile gamma
	TCanvas *c_strp_profile_ALICE_gamma_rate = new TCanvas(); //Strip profile gamma rate
	TCanvas *c_cluster_size_muon_ALICE = new TCanvas(); //Muon cluster size
	TCanvas *c_cluster_mult_muon_ALICE = new TCanvas(); //Muon cluster multiplicity
	TCanvas *c_cluster_size_gamma_ALICE = new TCanvas(); //Gamma cluster size
	TCanvas *c_cluster_mult_gamma_ALICE = new TCanvas(); //Gamma cluster multiplicity
	TCanvas *c_time_profile_analysis_strip[12][16];
	TCanvas *c_time_profile_strips[16]; //Time profile for each strip
	TCanvas *c_time_profile_ALICE_tmp = new TCanvas();
	TCanvas *c_time_spread = new TCanvas();
	TCanvas *c_first_and_second = new TCanvas();

	/////////////////////
	//				   //
	//HISTOGRAM BOOKING//
	//				   //
	/////////////////////

	//Define histograms in arrays of 12 elements, one per HV point
	TH1F *time_profile_ALICE[12];
	TH1F *time_profile_analysis_ALICE[12];
	TH1F *strp_profile_ALICE[12];
	TH1F *strp_profile_ALICE_gamma[12];
	TH1F *strp_profile_ALICE_gamma_rate[12]; 
	TH1F *strp_profile_ALICE_muons[12];
	TH1F *strp_profile_ALICE_muons_no8_no16[12]; //Muon strip profile if cs is not 16 or 8 
	TH1F *muon_cluster_size_ALICE[12];
	TH1F *muon_cluster_mult_ALICE[12];
	TH1F *gamma_cluster_size_ALICE[12];
	TH1F *gamma_cluster_mult_ALICE[12];
	TH1F *time_profile_strips[12][16];
	TH1I *h_first_and_second[12];
	TH1F *strp_profile_ALICE_tmp; //Temporary histo if the fit fails
	TH1F *time_profile_ALICE_tmp; //Temporary histo if the fit fails
	TH1F *mult_ALICE = new TH1F("h7","Multiplicity ALICE Y",33,-0.5,32.5);
	TH2F *event_display = new TH2F("h8","Event display ALICE Y",600,0,600,16,16.5,32.5);
	//TH1F *event_time_profile[12][16];

	/////////////////////
	//				   //
	//BOXES 4 WINDOW   //
	//				   //
	/////////////////////
	TBox *mu_window_total[12];
	TLine *line_start[12];
	TLine *line_mean[12];
	TLine *line_end[12];

	/////////////////////
	//				   //
	//		GRAPHS     //
	//				   //
	/////////////////////
	TGraphErrors *strip_distribution[12];

	////////////////////
	//				  //
	//USEFUL VARIABLES//
	//				  //
	////////////////////

	//Gamma window definition
	double gamma_start = 0, gamma_end = 4350; //Start and finish of the gamma window
	double gamma_window = (gamma_end-gamma_start)*1e-9;

	//HV, current, error on HV, error on current, current density and error on current density + env parameters
	vector<double> veff, imon, eveff, eimon, currdens, ecurrdens,vmon, evmon, veff_eco; 
	vector<double> vT, vp;
	double p0 = 990., T0 = 293., T, p, tens_mon, tens_eff; //mbar, K

	//To estimate fake efficiency: tot time = total measuring time, gamma_rate = in Hz, gamma rate outside mu window
	//muon_window_gamma = expected # of gammas in the mu window according to the rate
	//p_gamma = prob to have detected 0 gammas (Poisson statistics)
	//Gamma_eff = efficiency for gamma detection (fake efficiency) with error e_Gamma_eff
	double tot_time = 0, gamma_rate = 0, muon_window_gamma = 0, p_gamma = 0, Gamma_eff = 0, e_Gamma_eff = 0;

	//Muon efficiency and error, muon efficiency with gamma correction and error, fake efficiency and error
	vector<double> eff, e_eff, eff_gamma, e_eff_gamma, gammaefficiency, e_gammaefficiency;

	//Clustering: 
	//clus_time = clusterization time, muon_times/hits = time and strip of muon hits (inside mu window)
	//gamma_times/hits = time and strip of gamma hits (outside spill) 
	//muon_time_hit, gamma_time_hit = strip ordered pairs of strip-hit for muons and gammas
	//clusters = to save all the clusters in an event (each element is the number of adjacent strips)
	//coppia_strip/time = auxiliary vectors used to calculare cs
	//avg_muon/gamma.... vectors to save the avg cluster size/multiplicity in a scan (one element per HV value)
	//double clus_time = 15.;
	vector<float> muon_times, muon_hits, gamma_times, gamma_hits;
	vector<pair<float,float>> muon_time_hit, gamma_time_hit;
	vector <int> muon_clusters, gamma_clusters; 
	vector <int> coppia_strip;
	vector <double> coppia_time;
	vector <double> avg_muon_clus, e_avg_muon_clus, avg_gamma_clus, e_avg_gamma_clus;
	vector <double> avg_muon_clus_mul, e_avg_muon_clus_mul, avg_gamma_clus_mul, e_avg_gamma_clus_mul;

	//To calculate gamma rate: 
	//gamma_counter-muon_counter = to count how many gammas and muons in an HV point
	//rate_bin... = to calculate average rate with error in an HV point
	vector<double> avg_gamma_rate, e_avg_gamma_rate; 
	int gamma_counter = 0, muon_counter = 0; 
	double rate_bin = 0, rate_sum = 0, e_avg_rate = 0, partial = 0, rate_error = 0, rate_dev = 0;

	//////////////////////////
	//						//
	//		ANALYSIS		//
	//						//
	//////////////////////////	

	//ofstream width; //Output file for debug, in order to save the values of std deviation for the time profile fit

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

    	int secPeakCounter = 0;
    	int totalCounts = 0;

		rate_sum = 0; //Reset rate sum for average rate measurements

		string cartella = "HV" + to_string(point+1); //Create folders for different HV points in the out root file
		cdtof[point] = fout->mkdir(cartella.c_str());  
		cdtof[point]->cd();
		for (int i = 0; i < 16; i++) cstrip[i] = cdtof[point]->mkdir(("Strip " + to_string(i+1)).c_str());

		//Create histograms + set axes title
		time_profile_ALICE[point] = new TH1F(("a"+cartella).c_str(),("Time profile ALICE "+plane+cartella).c_str(),5000,0,5000);
		time_profile_ALICE[point]->GetXaxis()->SetTitle("Time[ns]");
		time_profile_ALICE[point]->GetYaxis()->SetTitle("Hits");
		strp_profile_ALICE[point] = new TH1F(("b"+cartella).c_str(),("Strip profile ALICE "+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE[point]->GetYaxis()->SetTitle("Hits");
		strp_profile_ALICE_muons[point] = new TH1F(("c"+cartella).c_str(),("Strip profile ALICE muons "+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_muons[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_muons[point]->GetYaxis()->SetTitle("Muon hits");
		time_profile_ALICE_tmp = new TH1F(("k"+cartella).c_str(),("Time profile ALICE (temporary)"+plane+cartella).c_str(),5000,0,5000);
		time_profile_ALICE_tmp->GetXaxis()->SetTitle("Time[ns]");
		time_profile_ALICE_tmp->GetYaxis()->SetTitle("Hits");
		strp_profile_ALICE_tmp = new TH1F(("l"+cartella).c_str(),("Strip profile ALICE (temporary)"+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_tmp->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_tmp->GetYaxis()->SetTitle("Hits");
		time_profile_analysis_ALICE[point] = new TH1F(("m"+cartella).c_str(),("Time profile ALICE (analysis)"+plane+cartella).c_str(),5000,0,5000);
		time_profile_analysis_ALICE[point]->GetXaxis()->SetTitle("Time[ns]");
		time_profile_analysis_ALICE[point]->GetYaxis()->SetTitle("Hits");
		h_first_and_second[point] = new TH1I(("n"+cartella).c_str(),("Hits first peak and second peak"+plane+cartella).c_str(),2,0,2);
		h_first_and_second[point]->GetXaxis()->SetTitle("0 = no, 1 = yes");
		h_first_and_second[point]->GetYaxis()->SetTitle("Counts");
		
		//Divide the canvas for the time profile for each strip
		c_time_profile_strips[point] = new TCanvas();
		c_time_profile_strips[point]->Divide(4,4);
		//Set titles of the axes for all the histograms of the strip time profile
		for (int i = 0; i < 16; i++) {
			c_time_profile_analysis_strip[point][i] = new TCanvas();
			time_profile_strips[point][i] = new TH1F(("HV"+to_string(point+1)+"strip"+to_string(i+1)).c_str(),("Time profile strip "+to_string(i+1)).c_str(),5000,0,5000);
			time_profile_strips[point][i]->GetXaxis()->SetTitle("Time [ns]");
			time_profile_strips[point][i]->GetYaxis()->SetTitle("Counts");
			//event_time_profile[point][i] = new TH1F(("HV"+to_string(point+1)+"strip"+to_string(i+1)).c_str(),("Time profile strip "+to_string(i+1)).c_str(),5000,0,5000);
		}
		
		//Open the different DAQ .root files
		string hvpoint = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_DAQ.root"; //For data
		string caen = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_CAEN.root"; //For HV and I
		string dip = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_DIP.root"; //For temperature and pressure
		//cout << "DAQ file:" << hvpoint << " CAEN file: " << caen << endl; //Debug

		//Define the file to be opened
		TFile *f = new TFile(hvpoint.c_str(),"READ");

		//Open the tree in the file
		TTree *t = (TTree*)f->Get("RAWData");

		//Define the elements in the tree
		int nev, eventry, quality, trig_type, trigger = 0, count_Y = 0, count_Y_gamma = 0, gamma_trigger = 0; //trig_type = 1 -> during spill, trigger_type = 0 -> outside spill
		bool hit_muon = 0, hit_gamma = 0; //To count if the strips have been fired in an event
		vector <int> *channel = 0; //TDC channels 
		vector <float> *timestamp = 0; //TDC timestamps
		vector <int> alice_hit;
		vector <double> alice_time;

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
			//cout << "event number " << i+1 << endl; //Debug
			if (trig_type == 1) trigger++; //count the number of muon triggers
			if (trig_type == 0) gamma_trigger++;

			for (unsigned int k = 0; k < channel->size(); k++) { //Go through the vector that contains all the channels of the event

				//All hits profile
				if (channel->at(k) <= finstrip && channel->at(k) >= instrip && trig_type == 1 && timestamp->at(k) > 100) { //ALICE chamber Y -> All the data is analyzed in the spill
					strp_profile_ALICE[point]->Fill(strp_map.lower_bound(channel->at(k))->second); //Strip profile of all hits
					time_profile_ALICE[point]->Fill(timestamp->at(k));//Time profile of the chamber
					//alice_hit.push_back(strp_map.lower_bound(channel->at(k))->second);
					alice_hit.push_back(channel->at(k));
					alice_time.push_back(timestamp->at(k));
					count_Y++;
					//cout << "strips in the event: " << strp_map.lower_bound(channel->at(k))->second << " "; //Debug
					//cout << channel->at(k) << "-" << timestamp->at(k) << "\t"; //Debug

					//Fill time profile strip per strip
					/*switch (channel->at(k)) {
						case 4064:
						time_profile_strips[point][0]->Fill(timestamp->at(k));
						break;
						case 4065:
						time_profile_strips[point][1]->Fill(timestamp->at(k));
						break;
						case 4066:
						time_profile_strips[point][2]->Fill(timestamp->at(k));
						break;
						case 4067:
						time_profile_strips[point][3]->Fill(timestamp->at(k));
						break;
						case 4068:
						time_profile_strips[point][4]->Fill(timestamp->at(k));
						break;
						case 4069:
						time_profile_strips[point][5]->Fill(timestamp->at(k));
						break;
						case 4070:
						time_profile_strips[point][6]->Fill(timestamp->at(k));
						break;
						case 4071:
						time_profile_strips[point][7]->Fill(timestamp->at(k));
						break;
						case 4072:
						time_profile_strips[point][8]->Fill(timestamp->at(k));
						break;
						case 4073:
						time_profile_strips[point][9]->Fill(timestamp->at(k));
						break;
						case 4074:
						time_profile_strips[point][10]->Fill(timestamp->at(k));
						break;
						case 4075:
						time_profile_strips[point][11]->Fill(timestamp->at(k));
						break;
						case 4076:
						time_profile_strips[point][12]->Fill(timestamp->at(k));
						break;
						case 4077:
						time_profile_strips[point][13]->Fill(timestamp->at(k));
						break;
						case 4078:
						time_profile_strips[point][14]->Fill(timestamp->at(k));
						break;
						case 4079:
						time_profile_strips[point][15]->Fill(timestamp->at(k));
						break;
					}*/		
				}
			} //end cycle on channel vector	

			/*int first_and_second = 0;

			for (unsigned int z = 0; z < alice_hit.size(); z++) {
				if (alice_time.at(z) > mu_window_start[point] && alice_time.at(z) < mu_window_end[point]) { //Hit in the muon window
					for (unsigned int zz = 0; zz < alice_hit.size(); zz++) {
						if (alice_time.at(zz) > 4545 && alice_time.at(zz) < 4600) {
							first_and_second = 1;
							break;
						}
					}
				}
			} //End for on alice_hit

			h_first_and_second[point]->Fill(first_and_second);*/
			int first_peak_count = 0;
			int second_peak_count = 0;


			for (unsigned int z = 0; z < alice_time.size(); z++) {
				if (alice_time.at(z) > 4480 && alice_time.at(z) < 4540) { //Hit in the muon window
					for (unsigned int zz = 0; zz < alice_time.size(); zz++) {
						time_profile_analysis_ALICE[point]->Fill(alice_time.at(zz));
					}
				}
				break;
			} //End for on alice_hit

			//Analysis of the secondary peaks in the time profile
			/*for (unsigned int z = 0; z < alice_hit.size(); z++) {
				if (alice_time.at(z) > 4545 && alice_time.at(z) < 4600) { //2nd peak -> if there is at least one hit in this time, draw all hits in the time profile
				//if (alice_time.at(z) > 4600 && alice_time.at(z) < 4640) { //3rd peak -> if there is at least one hit in this time, draw all hits in the time profile
					//time_profile_analysis_ALICE[point]->FillN(alice_time.size(),&alice_time[0],NULL,1);
					cout << "For event #: " << nev << " there is a hit in the 2nd peak!" << endl;
					secPeakCounter++;
					int numev = alice_time.size();
					totalCounts += alice_time.size();
					//cout << "Size of the timestamp vector: " << numev << endl;

					for (unsigned int zz = 0; zz < alice_time.size(); zz++) {
						//if (alice_time.at(zz) > mu_window_start[point] && alice_time.at(zz) < mu_window_end[point]) {
							time_profile_analysis_ALICE[point]->Fill(alice_time.at(zz));
							if (alice_hit.at(zz) == tdc_ch.at(0)) time_profile_strips[point][0]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(1)) time_profile_strips[point][1]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(2)) time_profile_strips[point][2]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(3)) time_profile_strips[point][3]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(4)) time_profile_strips[point][4]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(5)) time_profile_strips[point][5]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(6)) time_profile_strips[point][6]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(7)) time_profile_strips[point][7]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(8)) time_profile_strips[point][8]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(9)) time_profile_strips[point][9]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(10)) time_profile_strips[point][10]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(11)) time_profile_strips[point][11]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(12)) time_profile_strips[point][12]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(13)) time_profile_strips[point][13]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(14)) time_profile_strips[point][14]->Fill(alice_time.at(zz));
							else if (alice_hit.at(zz) == tdc_ch.at(15)) time_profile_strips[point][15]->Fill(alice_time.at(zz));
						//}
					} //End of for on zz
					break;
				} //End if ALICE time is in the 2nd peak
				
				if (alice_time.at(z) > mu_window_start[point] && alice_time.at(z) < mu_window_end[point]) { //Hit in the muon window

					if (alice_hit.at(z) == tdc_ch.at(0)) time_profile_strips[point][0]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(1)) time_profile_strips[point][1]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(2)) time_profile_strips[point][2]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(3)) time_profile_strips[point][3]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(4)) time_profile_strips[point][4]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(5)) time_profile_strips[point][5]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(6)) time_profile_strips[point][6]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(7)) time_profile_strips[point][7]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(8)) time_profile_strips[point][8]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(9)) time_profile_strips[point][9]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(10)) time_profile_strips[point][10]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(11)) time_profile_strips[point][11]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(12)) time_profile_strips[point][12]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(13)) time_profile_strips[point][13]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(14)) time_profile_strips[point][14]->Fill(alice_time.at(z));
					else if (alice_hit.at(z) == tdc_ch.at(15)) time_profile_strips[point][15]->Fill(alice_time.at(z));
					for (unsigned int zz = 0; zz < alice_time.size(); zz++) {
						time_profile_analysis_ALICE[point]->Fill(alice_time.at(zz));

						if (alice_hit.at(zz) == 4064) time_profile_strips[point][0]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4065) time_profile_strips[point][1]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4066) time_profile_strips[point][2]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4067) time_profile_strips[point][3]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4068) time_profile_strips[point][4]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4069) time_profile_strips[point][5]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4070) time_profile_strips[point][6]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4071) time_profile_strips[point][7]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4072) time_profile_strips[point][8]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4073) time_profile_strips[point][9]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4074) time_profile_strips[point][10]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4075) time_profile_strips[point][11]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4076) time_profile_strips[point][12]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4077) time_profile_strips[point][13]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4078) time_profile_strips[point][14]->Fill(alice_time.at(zz));
						else if (alice_hit.at(zz) == 4079) time_profile_strips[point][15]->Fill(alice_time.at(zz));
					} //End of for on zz
					//break;
				}
			}*/ //End for on alice_hit
			alice_hit.clear();
			alice_time.clear();	

			fout->cd(); //open .root file
			cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
			cstrip[0]->cd();
			//event_time_profile[point][0]->SetOption("COLZ");
			//event_time_profile[point][0]->Write(("Muon event #: " + to_string(count_Y) + " HV point: " + to_string(point)).c_str()); //Save time profile
			//event_time_profile[point][0]->Reset();	
		} //End of cycle on all tree entries

		cout << endl << endl << "Times in 2nd peak: " << secPeakCounter << endl;
		cout << endl << endl << "Total hits: " << totalCounts << endl;

		c_first_and_second->cd();
		h_first_and_second[point]->Draw("HISTO");

		//Put a shadowed box over the mu window
		mu_window_total[point] = new TBox(mu_window_start[point],0.,mu_window_end[point],time_profile_ALICE[point]->GetMaximum());
		mu_window_total[point]->SetFillColor(kOrange);
		mu_window_total[point]->SetFillStyle(3144);

		//Horizontal lines to identify the mean and the start/end
		line_start[point] = new TLine(strips.front(),mu_window_start.at(point),strips.back(),mu_window_start.at(point));
		line_start[point]->SetLineColor(kGreen);
		line_mean[point] = new TLine(strips.front(),mu_window_mean.at(point),strips.back(),mu_window_mean.at(point));
		line_mean[point]->SetLineColor(kRed);
		line_end[point] = new TLine(strips.front(),mu_window_end.at(point),strips.back(),mu_window_end.at(point));
		line_end[point]->SetLineColor(kGreen);

		//Vectors to store the mean value and error of mu window for each strip
		vector<double> mu_mean_time, mu_mean_time_error;

		for (unsigned int j = 0; j < strips.size(); j++) {
			mu_mean_time.push_back(time_profile_strips[point][j]->GetMean()); //Mean 
			//mu_mean_time_error.push_back(time_profile_strips[point][j]->GetStdDev()); //Std dev = old error
			//mu_mean_time_error.push_back(time_profile_strips[point][j]->GetStdDev()/time_profile_strips[point][j]->GetEntries());//Mean error = new error
			mu_mean_time_error.push_back(time_profile_strips[point][j]->GetMeanError());//Mean error = new error
			//cout << "Mio error (diviso n): " << time_profile_strips[point][j]->GetStdDev()/TMath::Sqrt(time_profile_strips[point][j]->GetEntries()) << endl;
			//cout << "Mio error (diviso n-1): " << time_profile_strips[point][j]->GetStdDev()/TMath::Sqrt(time_profile_strips[point][j]->GetEntries() - 1) << endl;
			//cout << "Histo error: " << time_profile_strips[point][j]->GetMeanError() << endl;
		}

		strip_distribution[point] = new TGraphErrors(strips.size(),&strips[0],&mu_mean_time[0],NULL,&mu_mean_time_error[0]);
		strip_distribution[point]->SetTitle(("Mean mu window HV " + to_string(point+1)).c_str());
		strip_distribution[point]->GetXaxis()->SetTitle("Strip");
		strip_distribution[point]->GetYaxis()->SetTitle("Mean mu window time [ns]");
		strip_distribution[point]->GetYaxis()->SetRangeUser(mu_window_start.at(point)-5,mu_window_end.at(point)+5);
		strip_distribution[point]->SetMarkerStyle(8);
		strip_distribution[point]->SetMarkerSize(1);

		c_time_spread->cd();
		strip_distribution[point]->Draw("AP");
		line_start[point]->Draw("SAME");
		line_mean[point]->Draw("SAME");
		line_end[point]->Draw("SAME");

		c_time_profile_analysis_ALICE->cd();
		time_profile_analysis_ALICE[point]->Draw("HISTO");
		//mu_window_total[point]->Draw("SAME");

		//Draw the time profile of each strip of the RPC
		for (int i = 0; i < 16; i++) {
			c_time_profile_strips[point]->cd(i+1);
			time_profile_strips[point][i]->Draw("HISTO");

			c_time_profile_analysis_strip[point][i]->cd();
			time_profile_strips[point][i]->Draw("HISTO");
		}

		double start_muon, end_muon, mean;
		
		//Draw the tim profile with more complicated gaussian + plo0 fit to extract the signal window
		c_time_profile_ALICE->cd();
		muon_gamma->SetParLimits(2,4485.,4535.);
		muon_gamma->SetParLimits(1,0.,time_profile_ALICE[point]->GetBinContent(time_profile_ALICE[point]->GetMaximumBin()));
		muon_gamma->SetParLimits(3,0.,40.); //This works
		time_profile_ALICE[point]->Fit("gamma_flat","RM+"); //Horizontal line fit to gamma background
		double avg_gamma = gamma_flat->GetParameter(0); //Get the value of the fit
		muon_gamma->SetParameter(0,avg_gamma); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		time_profile_ALICE[point]->Fit("muon_gamma","RM+"); //Fit with gaussian + constant
		string minuitstatus = string(gMinuit->fCstatu);
		cout << endl << endl << "Status: " << minuitstatus << endl << endl;
		int r = 0;
		if(minuitstatus.compare("CONVERGED ") != 0 && minuitstatus.compare("OK        ") != 0) r = -1;
		time_profile_ALICE[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICE[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICE[point]->Draw("HISTO");
		//gamma_flat->Draw("SAME");
		muon_gamma->Draw("SAME");
		//mu_window_total[point]->Draw("SAME");

		if (r == 0) {
			start_muon = muon_gamma->GetParameter(2) - 3*muon_gamma->GetParameter(3); 
			end_muon = muon_gamma->GetParameter(2) + 3*muon_gamma->GetParameter(3);
			mean = muon_gamma->GetParameter(2);
		}

		if (r != 0) { //Fit fails, open the last HV point file, do the analysis there and use this as the muon window
			c_time_profile_ALICE_tmp->cd();
			cout << "Fit of hvpoint: " << point+1 << " failed, analyzing the last point" << endl;
			cout << "Last point: " << file_count/3 << endl;

			vector <int> *channel_last = 0; //TDC channels 
			vector <float> *timestamp_last = 0; //TDC timestamps
			int trig_type_last;

			//Define the file to be opened
			TFile *flast = new TFile(("Scan00" +to_string(scan) + "_HV" + to_string(file_count/3) + "_DAQ.root").c_str(),"READ");
			cout << "Scan00" +to_string(scan) + "_HV" + to_string(file_count/3) + "_DAQ.root" << endl;

			//Open the tree in the file
			TTree *tlast = (TTree*)flast->Get("RAWData");
			tlast->SetBranchAddress("TDC_channel", &channel_last);
			tlast->SetBranchAddress("TDC_TimeStamp",&timestamp_last);
			tlast->SetBranchAddress("TriggerTag",&trig_type_last);

			for (int i = 0; i < nentries; i++) { 
			tlast->GetEntry(i);

				for (unsigned int k = 0; k < channel_last->size(); k++) { //Go through the vector that contains all the channels of the event

					//All hits profile
					if (channel_last->at(k) <= finstrip && channel_last->at(k) >= instrip && trig_type_last == 1 && timestamp_last->at(k) > 100) { //ALICE chamber Y -> All the data is analyzed in the spill
						strp_profile_ALICE_tmp->Fill(strp_map.lower_bound(channel_last->at(k))->second); //Strip profile of all hits
						time_profile_ALICE_tmp->Fill(timestamp_last->at(k));//Time profile of the chamber
					}
				}
			}

			muon_gamma_temp->SetParLimits(1,0.,time_profile_ALICE_tmp->GetBinContent(time_profile_ALICE_tmp->GetMaximumBin()));
			muon_gamma_temp->SetParLimits(2,4485.,4535.); //This works
			muon_gamma->SetParLimits(3,0.,40.); //This works
			time_profile_ALICE_tmp->Fit("gamma_flat_temp","RM+"); //Horizontal line fit to gamma background
			double avg_gamma_temp = gamma_flat_temp->GetParameter(0); //Get the value of the fit
			muon_gamma_temp->SetParameter(0,avg_gamma_temp); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
			//muon_gamma->FixParameter(0,avg_gamma); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
			time_profile_ALICE_tmp->Fit("muon_gamma_temp","RM+"); //Fit with gaussian + constant

			//Start and end time of the muon window
			start_muon = muon_gamma_temp->GetParameter(2) - 3*muon_gamma_temp->GetParameter(3); 
			end_muon = muon_gamma_temp->GetParameter(2) + 3*muon_gamma_temp->GetParameter(3);
			mean = muon_gamma_temp->GetParameter(2);

		}  //End if on the fit failure

		double muon_window = (end_muon - start_muon)*1e-9; //Muon window in s
		cout << endl << "muon window: " << start_muon << " ns < " << mean << " ns < " << end_muon << " ns" << endl;

		//Draw strip profile, all hits
		c_strp_profile_ALICE->cd();
		strp_profile_ALICE[point]->Draw("HISTO");


		/////////////////
		//	PROFILES   //
		/////////////////	

		//Strip profile, muons
		c_strp_profile_ALICE_muons->cd();
		muon_strip->SetParLimits(1,10.,13.);
		strp_profile_ALICE_muons[point]->Fit("muon_strip","RM+");
		strp_profile_ALICE_muons[point]->Draw("HISTO");
		muon_strip->Draw("SAME");

		fout->cd(); //open .root file
		cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
		c_time_profile_ALICE->Write(("Time_profile_ALICE_"+plane+"_"+cartella).c_str()); //Save time profile
		c_time_profile_analysis_ALICE->Write(("Time_profile_ALICE_analysis_"+plane+"_"+cartella).c_str());
		c_time_profile_strips[point]->Write(("Strip_time_profile_ALICE_"+plane+"_"+cartella).c_str());
		c_time_spread->Write(("Time spread in mu window HV " + to_string(point+1)).c_str());
		c_first_and_second->Write(("Hits_first_and_second_peak"+plane+cartella).c_str());
		for (int i = 0; i < 16; i++) {
			c_time_profile_analysis_strip[point][i]->Write(("Time_profile_ALICE_"+plane+"_strip_"+to_string(i+1)+"_"+cartella).c_str());
		}
		//c_strp_profile_ALICE->Write(("Strip_profile_ALICE_all_hits_"+plane+"_"+cartella).c_str()); //Save strip profile (all hits)
		//c_strp_profile_ALICE_muons->Write(("Strip_profile_ALICE_muons_"+plane+"_"+cartella).c_str());
		//c_strp_profile_ALICE_gamma->Write(("Strip_profile_ALICE_gammas_"+plane+"_"+cartella).c_str());
		//c_strp_profile_ALICE_gamma_rate->Write(("Strip_profile_ALICE_gammas_rate_"+plane+"_"+cartella).c_str());
		//c_cluster_size_muon_ALICE->Write(("Muon_cluster_Size_ALICE_"+plane+"_"+cartella).c_str());
		//c_cluster_mult_muon_ALICE->Write(("Muon_cluster_Mult_ALICE_"+plane+"_"+cartella).c_str());
		//c_cluster_size_gamma_ALICE->Write(("Gamma_cluster_Size_ALICE_"+plane+"_"+cartella).c_str());
		//c_cluster_mult_gamma_ALICE->Write(("Gamma_cluster_Mult_ALICE_"+plane+"_"+cartella).c_str());

	} //End of cycle on the HV points in a scan
}