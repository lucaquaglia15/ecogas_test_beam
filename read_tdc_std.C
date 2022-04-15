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

void read_tdc_std(const int scan, const double clus_time) {
	
	//string plane = "X";
	string plane = "Y";

	//Histograms cosmetics
	gStyle->SetHistFillColor(kBlue);
	gStyle->SetHistFillStyle(3357);
	gStyle->SetHistLineColor(kBlue);

	//Show fit results on the plots
	gStyle->SetOptFit(1111);
	gStyle->SetDrawOption("COLZ");

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

	instrip = 4072;
	finstrip = 4079;

	//cout << instrip << "\t" << finstrip << endl; //Debug

	map <int, int> strp_map;
	for (unsigned int i = 0; i < tdc_ch.size();i++) {
		strp_map.insert(pair<int, int>(tdc_ch.at(i),rpc_strp.at(i)));
	}

	//RPC parameters to calculate rates
	double strip_area = 150; //cm2, area of a strip
	double meas_time = 5000; //ns, time of a single measurement
	double area = 2500; //Area of the ALICE RPC

	ofstream ftimesX;
	ofstream ftimesY;
	if (plane == "X") {
		remove(("Timestamps_X_scan"+to_string(scan)+".txt").c_str());
		ftimesX.open(("Timestamps_X_scan"+to_string(scan)+".txt").c_str(),ios::app);
	}
	else if (plane == "Y") {
		remove(("Timestamps_Y_scan"+to_string(scan)+".txt").c_str());
		ftimesY.open(("Timestamps_Y_scan"+to_string(scan)+".txt").c_str(),ios::app);	
	} 

	//Enter in the scan folder 
	string fol = "/media/luca/Elements/TB_October_2021/STD/Scan_00" + to_string(scan);

	//Create .root file with all the output data
	TFile *fout = new TFile(("outfile_ALICE"+plane+to_string(scan)+".root").c_str(),"RECREATE");

	fout->cd();
	TDirectory *cdtof[12]; //12 directories in the root file, one for each HV point

	gSystem->cd(fol.c_str()); //Enter folder

 	//////////////////////
	//					//
	//	CANVAS BOOKING	//
	//					//
	//////////////////////

	TCanvas *c_time_profile_ALICE = new TCanvas(); //Time profile
	TCanvas *c_strp_profile_ALICE = new TCanvas(); //Strip profile of all hits
	TCanvas *c_strp_profile_ALICE_muons = new TCanvas(); //Strip profile muons
	TCanvas *c_strp_profile_ALICE_muons_no8_no16 = new TCanvas(); //Strip profile muons if cs is not 16 or 8
	TCanvas *c_strp_profile_ALICE_gamma = new TCanvas(); //Strip profile gamma
	TCanvas *c_strp_profile_ALICE_gamma_rate = new TCanvas(); //Strip profile gamma rate
	TCanvas *c_cluster_size_muon_ALICE = new TCanvas(); //Muon cluster size
	TCanvas *c_cluster_mult_muon_ALICE = new TCanvas(); //Muon cluster multiplicity
	TCanvas *c_cluster_size_gamma_ALICE = new TCanvas(); //Gamma cluster size
	TCanvas *c_cluster_mult_gamma_ALICE = new TCanvas(); //Gamma cluster multiplicity
	TCanvas *c_time_profile_strips[16]; //Time profile for each strip
	TCanvas *c_time_profile_ALICE_tmp = new TCanvas();
	TCanvas *c_event_display = new TCanvas(); //Canvas for event display
	TCanvas *c_deltaTime = new TCanvas();

	/////////////////////
	//				   //
	//HISTOGRAM BOOKING//
	//				   //
	/////////////////////

	//Define histograms in arrays of 12 elements, one per HV point
	TH1F *time_profile_ALICE[12];
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
	TH1F *strp_profile_ALICE_tmp; //Temporary histo if the fit fails
	TH1F *time_profile_ALICE_tmp; //Temporary histo if the fit fails
	TH1F *mult_ALICE = new TH1F("h7","Multiplicity ALICE Y",33,-0.5,32.5);
	TH2F *event_display = new TH2F("h8","Event display ALICE Y",300,4400,4700,16,4064.5,4079.5);
	TH1F *deltaTime[12];
	TPaveText *tdelta[12];
	TLine *after15[12];

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
	vector<double> muon_times, muon_hits, gamma_times, gamma_hits;
	vector<pair<double,double>> muon_time_hit, gamma_time_hit;
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

		rate_sum = 0; //Reset rate sum for average rate measurements

		string cartella = "HV" + to_string(point+1); //Create folders for different HV points in the out root file
		cdtof[point] = fout->mkdir(cartella.c_str());   

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
		strp_profile_ALICE_muons_no8_no16[point] = new TH1F(("d"+cartella).c_str(),("Strip profile ALICE muons (no cs 8 or 16) "+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_muons_no8_no16[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_muons_no8_no16[point]->GetYaxis()->SetTitle("Muon hits");
		strp_profile_ALICE_gamma[point] = new TH1F(("e"+cartella).c_str(),("Strip profile ALICE gamma "+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_gamma[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_gamma[point]->GetYaxis()->SetTitle("Gamma hits");
		strp_profile_ALICE_gamma_rate[point] = new TH1F(("f"+cartella).c_str(),("Strip profile ALICE gamma rate "+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_gamma_rate[point]->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_gamma_rate[point]->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
		muon_cluster_size_ALICE[point] = new TH1F(("g"+cartella).c_str(),("Muon cluster size ALICE "+plane+cartella).c_str(),16,0,16);
		muon_cluster_size_ALICE[point]->GetXaxis()->SetTitle("Muon cluster size [strips]");
		muon_cluster_size_ALICE[point]->GetYaxis()->SetTitle("Counts");
		muon_cluster_mult_ALICE[point] = new TH1F(("h"+cartella).c_str(),("Muon cluster mult ALICE "+plane+cartella).c_str(),16,0,16);
		muon_cluster_mult_ALICE[point]->GetXaxis()->SetTitle("Muon cluster multiplicity [strips]");
		muon_cluster_mult_ALICE[point]->GetYaxis()->SetTitle("Counts");
		gamma_cluster_size_ALICE[point] = new TH1F(("i"+cartella).c_str(),("Gamma cluster size ALICE "+plane+cartella).c_str(),16,0,16);
		gamma_cluster_size_ALICE[point]->GetXaxis()->SetTitle("Gamma cluster size [strips]");
		gamma_cluster_size_ALICE[point]->GetYaxis()->SetTitle("Counts");
		gamma_cluster_mult_ALICE[point] = new TH1F(("j"+cartella).c_str(),("Gamma cluster mult ALICE "+plane+cartella).c_str(),16,0,16);
		gamma_cluster_mult_ALICE[point]->GetXaxis()->SetTitle("Gamma cluster multiplicity [strips]");
		gamma_cluster_mult_ALICE[point]->GetYaxis()->SetTitle("Counts");
		time_profile_ALICE_tmp = new TH1F(("k"+cartella).c_str(),("Time profile ALICE (temporary)"+plane+cartella).c_str(),5000,0,5000);
		time_profile_ALICE_tmp->GetXaxis()->SetTitle("Time[ns]");
		time_profile_ALICE_tmp->GetYaxis()->SetTitle("Hits");
		strp_profile_ALICE_tmp = new TH1F(("l"+cartella).c_str(),("Strip profile ALICE (temporary)"+plane+cartella).c_str(),16,0.5,16.5);
		strp_profile_ALICE_tmp->GetXaxis()->SetTitle("Strip");
		strp_profile_ALICE_tmp->GetYaxis()->SetTitle("Hits");
		deltaTime[point] = new TH1F(("m"+cartella).c_str(),("Maximum delta time in cluster"+plane+cartella).c_str(),500,-0.5,50.5);
		deltaTime[point]->GetXaxis()->SetTitle("#DELTA T[ns]");
		deltaTime[point]->GetYaxis()->SetTitle("Hits");
		tdelta[point] = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");

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
		int muon_ev, gamma_ev;
		bool hit_muon = 0, hit_gamma = 0; //To count if the strips have been fired in an event
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
			//cout << "event number " << i+1 << endl; //Debug
			if (trig_type == 1) trigger++; //count the number of triggers in spill
			if (trig_type == 0) gamma_trigger++; //count the number of triggers out of spill

			for (unsigned int k = 0; k < channel->size(); k++) { //Go through the vector that contains all the channels of the event

				//All hits profile
				if (channel->at(k) <= finstrip && channel->at(k) >= instrip && trig_type == 1 && timestamp->at(k) > 100) { //ALICE chamber Y -> All the data is analyzed in the spill
						strp_profile_ALICE[point]->Fill(strp_map.lower_bound(channel->at(k))->second); //Strip profile of all hits
						time_profile_ALICE[point]->Fill(timestamp->at(k));//Time profile of the chamber
						count_Y++;
						//cout << "strips in the event: " << strp_map.lower_bound(channel->at(k))->second << " "; //Debug
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

		/*
		//Draw time profile with simple gaussian fit to extract the signal window
		c_time_profile_ALICE->cd();
		muon_sig->SetParLimits(1,4500.,4600.);
		muon_sig->SetParLimits(2,0.,30.);
		time_profile_ALICE[point]->Fit("muon_sig","RM+"); //Fit with pure gaussian 
		time_profile_ALICE[point]->GetXaxis()->SetTitle("Time [ns]");
		time_profile_ALICE[point]->GetYaxis()->SetTitle("Events");
		time_profile_ALICE[point]->Draw("HISTO");
		muon_sig->Draw("SAME");

		//Write the width of the pure gaussian to an output txt file
		width.open("std_dev.txt",std::ios_base::app);
		width << to_string(point+1) << "\t" << muon_sig->GetParameter(2) << "\n";
		width.close();

		//Start and end time of the muon window
		double start_muon = muon_sig->GetParameter(1) - 3*muon_sig->GetParameter(2); 
		double end_muon = muon_sig->GetParameter(1) + 3*muon_sig->GetParameter(2);
		double mean = muon_sig->GetParameter(1);
		*/

		double start_muon, end_muon, mean;
		
		//Draw the time profile with more complicated gaussian + plo0 fit to extract the signal window
		c_time_profile_ALICE->cd();
		muon_gamma->SetParLimits(1,0.,time_profile_ALICE[point]->GetBinContent(time_profile_ALICE[point]->GetMaximumBin()));
		muon_gamma->SetParLimits(2,4485.,4535.);
		muon_gamma->SetParLimits(3,0.,40.); //This works
		time_profile_ALICE[point]->Fit("gamma_flat","RM+"); //Horizontal line fit to gamma background
		double avg_gamma = gamma_flat->GetParameter(0); //Get the value of the fit
		muon_gamma->SetParameter(0,avg_gamma); //Set the constant of the gaussian to the value obtained with the pol0 fit and leave it as a free parameter in the fit
		//muon_gamma->FixParameter(0,avg_gamma); //Fix the constant of the gaussian to the value obtained with the pol0 fit and not change it in the new fit
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

		if (r == 0) {
			start_muon = muon_gamma->GetParameter(2) - 3*muon_gamma->GetParameter(3); 
			end_muon = muon_gamma->GetParameter(2) + 3*muon_gamma->GetParameter(3);
			mean = muon_gamma->GetParameter(2);
		}

		/*
		//Write the width of the gaussian (plus the pol0) to an output txt file
		width.open("std_dev.txt",std::ios_base::app);
		width << to_string(point+1) << "\t" << muon_gamma->GetParameter(3) << "\n";
		width.close();
		*/

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

		//Start and end time of the muon window
		//double start_muon = muon_gamma->GetParameter(2) - 3*muon_gamma->GetParameter(3); 
		//double end_muon = muon_gamma->GetParameter(2) + 3*muon_gamma->GetParameter(3);
		//double mean = muon_gamma->GetParameter(2);

		double muon_window = (end_muon - start_muon)*1e-9; //Muon window in s
		cout << endl << "muon window: " << start_muon << " ns < " << mean << " ns < " << end_muon << " ns" << endl;

		if (plane == "X") ftimesX << start_muon << "\t" << mean << "\t" << end_muon << "\t" << muon_window << "\n";
		else if (plane == "Y") ftimesY << start_muon << "\t" << mean << "\t" << end_muon << "\t" << muon_window << "\n";

		//Draw strip profile, all hits
		c_strp_profile_ALICE->cd();
		strp_profile_ALICE[point]->Draw("HISTO");

		///////////////////
		//Muon-gamma hits//
		///////////////////

		for (int i = 0; i < nentries; i++) { 
			muon_times.clear();
			muon_hits.clear();
			muon_time_hit.clear();
			gamma_times.clear();
			gamma_hits.clear();
			gamma_time_hit.clear();

			t->GetEntry(i); //Analuze the event

			for (unsigned int k = 0; k < channel->size(); k++) { //Go through the vector that contains all the channels of the event

				//Muon hits profile
				if ((channel->at(k) <= finstrip && channel->at(k) >= instrip && trig_type == 1) && (timestamp->at(k) <= end_muon && timestamp->at(k) >= start_muon)) { //ALICE chamber -> muons
						strp_profile_ALICE_muons[point]->Fill(strp_map.lower_bound(channel->at(k))->second); //Strip profile of muons
						//VERY IMPORTANT to uncomment the following line if we want to use the normal eff definition. Use the following if we want the hit to be in the beam "strips"
						hit_muon = 1; //The chamber saw something, set the counter to true, it wil be reset at the end of the event
						muon_counter++;
						//cout << "K: " << k << ", time:" << timestamp->at(k) << ", channel: " << channel->at(k) << " "; //Debug for cluster size
						muon_times.push_back(timestamp->at(k));
						muon_hits.push_back(channel->at(k));
						muon_ev = nev;
				}
				
				//Gamma hits profile (in spill but before mu window)
				if ((channel->at(k) <= finstrip && channel->at(k) >= instrip && trig_type == 1) && (timestamp->at(k) <= gamma_end && timestamp->at(k) >= gamma_start) && timestamp->at(k) > 100) { //ALICE chamber Y -> gamma profile
					//strp_profile_ALICE_gamma[point]->Fill(strp_map.lower_bound(channel->at(k))->second); //Strip profile of gammas
					gamma_counter++;
					hit_gamma = 1;
					//gamma_times.push_back(timestamp->at(k));
					//gamma_hits.push_back(channel->at(k));
				}
				
				//Gamma hits profile (outside spill)
				if ((channel->at(k) <= finstrip && channel->at(k) >= instrip && trig_type == 0) && timestamp->at(k) > 100) { //ALICE chamber Y -> gamma profile
					strp_profile_ALICE_gamma[point]->Fill(strp_map.lower_bound(channel->at(k))->second); //Strip profile of gammas
					//gamma_counter++;
					gamma_times.push_back(timestamp->at(k));
					gamma_hits.push_back(channel->at(k));
					gamma_ev = nev;
				}
			}

			if (hit_muon == 1) count_Y++; //Increase the number of events seen by the RPC
			if (hit_gamma == 1) count_Y_gamma++;

			//Event display
			if (hit_muon != 0) {
				event_display->SetTitle(("Muon event #: " + to_string(count_Y) + " HV point: " + to_string(point)).c_str());
				event_display->FillN(muon_times.size(),&muon_times[0],&muon_hits[0],NULL,1);
				fout->cd(); //open .root file
				cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
				event_display->SetOption("COLZ");
				//event_display->Write(("Muon event #: " + to_string(count_Y) + " HV point: " + to_string(point)).c_str()); //Save time profile
				event_display->Reset();
			}

			hit_muon = 0; //reset the bool for muon events counting
			hit_gamma = 0; //reset the bool for gamma events counting

			//Cluster Size = mean number of strips next to each other that fired in an event
			//Multiplicity = total numer of strips fired in an event
			//Cluster Multiplicity = number of clusters in an event

			//////////////////////
			//	Muon Clustering	//
			//////////////////////

			int muon_pass = 0, muon_tmpcluster = 1; //integer to keep track of the pass along the hits vector
			unsigned int muon_counter_tmp = 1;
			//if (hits.size() == 0) {
			//	muon_cluster_size_ALICE[point]->Fill(0);
			//}
			
			if (muon_hits.size() == 1) {
			//else if (muon_hits.size() == 1) {
				muon_cluster_size_ALICE[point]->Fill(1);
				muon_cluster_mult_ALICE[point]->Fill(1);
				//deltaTime[point]->Fill(0);
			}

			else if (muon_hits.size() > 1) {
				//cout << "Size hit: " << muon_hits.size() << endl;
				muon_counter_tmp = 1;
				for (unsigned int hh = 0; hh < muon_times.size(); hh++) {
					muon_time_hit.push_back(make_pair(muon_hits.at(hh),muon_times.at(hh))); //hit - time
					//cout << "I'm here" << endl; //Debug
					//cout << endl; //Debug to be uncommented
				}

				//cout << "Time_hit size: " << muon_time_hit.size() << endl;

				/*for (unsigned int hh = 0; hh < muon_time_hit.size(); hh++) { //Normal printout
					cout << muon_time_hit.at(hh).first << "\t" << muon_time_hit.at(hh).second << "\t";
				}
				cout << endl;*/
				
				sort(muon_time_hit.begin(),muon_time_hit.end());

				/*for (unsigned int hh = 0; hh < muon_time_hit.size(); hh++) { //Ordered printout
					cout << muon_time_hit.at(hh).first << "\t" << muon_time_hit.at(hh).second << "\t";
				}
				cout << endl;*/

				do {
					//cout << "Size hit_time decreasing: " << muon_time_hit.size() << endl;
					double muon_first_strip = 0, muon_first_time = 0, muon_strip = 0, muon_time = 0;
					if (muon_pass == 0) { //Define the first element of the vector, compare all the others to that
						//cout << "Muon pass: " << muon_pass << endl;
						muon_first_strip = muon_time_hit.at(0).first;
						muon_first_time = muon_time_hit.at(0).second;
					}

					if (muon_pass == 0 && muon_time_hit.size() == 1) { //To take care of the last element of the vetor if it's left alone
						muon_clusters.push_back(1);
						muon_time_hit.erase(muon_time_hit.begin());
						muon_counter_tmp = 1;
						continue;
					}

					//cout << "Assigning the other strips!" << endl;
					
					if (muon_counter_tmp < muon_time_hit.size()) {
						muon_strip = muon_time_hit.at(muon_counter_tmp).first;
						muon_time = muon_time_hit.at(muon_counter_tmp).second;						
					}

					if (muon_pass == 0 && (TMath::Abs(muon_strip-muon_first_strip) == 1 || TMath::Abs(muon_strip-muon_first_strip) == 0) && TMath::Abs(muon_time-muon_first_time) <= clus_time && muon_time_hit.size() == 2) {
						muon_clusters.push_back(2);
						deltaTime[point]->Fill(TMath::Abs(muon_time-muon_first_time));
						break;
					}

					if (muon_pass == 0 && (TMath::Abs(muon_strip-muon_first_strip) == 1 || TMath::Abs(muon_strip-muon_first_strip) == 0) && TMath::Abs(muon_time-muon_first_time) <= clus_time) { //first muon_pass on the array, compare first element to all other vector elements
						//cout << "Adjacent pair!" << endl;
						coppia_strip.push_back(muon_first_strip);
						coppia_strip.push_back(muon_strip);
						coppia_time.push_back(muon_first_time);
						coppia_time.push_back(muon_time);
						muon_pass++;
						muon_tmpcluster++;
						muon_time_hit.erase(muon_time_hit.begin()+muon_counter_tmp);
						muon_time_hit.erase(muon_time_hit.begin());
						muon_counter_tmp = 1;
						continue;
					}

					else if (muon_pass > 0) { //There is a couple created
						for (unsigned int tmp = 0; tmp < coppia_strip.size(); tmp++) {
							double tmpstrip = coppia_strip.at(tmp);
							double tmptime = coppia_time.at(tmp);

							for (unsigned int tmp2 = 0; tmp2 < muon_time_hit.size(); tmp2++) {
								double tmpstrip2 = muon_time_hit.at(tmp2).first;
								double tmptime2 = muon_time_hit.at(tmp2).second;

								if((TMath::Abs(tmpstrip2-tmpstrip) == 1 || TMath::Abs(tmpstrip2-tmpstrip) == 0) && TMath::Abs(tmptime2-tmptime) <= clus_time) {
									coppia_strip.push_back(tmpstrip2);
									coppia_time.push_back(tmptime2);
									muon_time_hit.erase(muon_time_hit.begin()+tmp2);
									muon_tmpcluster++;
									muon_counter_tmp = 1;
									continue;
								}
							}
						}
						muon_clusters.push_back(muon_tmpcluster); //Uncomment
						//cout << " " << muon_tmpcluster << " "; //Debug
						muon_tmpcluster = 1; //Uncomment
						//cout << "Size: " << coppia_strip.size() << endl; //Debug
						double deltaMax = *max_element(coppia_time.begin(), coppia_time.end()) - *min_element(coppia_time.begin(), coppia_time.end());
						//cout << "Max Element = " << *max_element(coppia_time.begin(), coppia_time.end()); //Debug
						deltaTime[point]->Fill(deltaMax);
						/*if (deltaMax <= 15) {
							muon_clusters.push_back(muon_tmpcluster);
						}*/
						//muon_tmpcluster = 1; //Uncomment if total clustering time can be over 15 ns
						coppia_strip.clear();
						coppia_time.clear();
						muon_counter_tmp = 1;
						muon_pass = 0;
						continue;
					}

					else if ((muon_pass == 0 && TMath::Abs(muon_strip-muon_first_strip) > 1) || (muon_pass == 0 && TMath::Abs(muon_time-muon_first_time) > clus_time)) {
						//cout << "Not adjacent pair on the first go!" << endl;
						muon_counter_tmp++;
						if (muon_counter_tmp == muon_time_hit.size()) {
							muon_time_hit.erase(muon_time_hit.begin());
							muon_clusters.push_back(1);
							muon_counter_tmp = 1;
							continue;
						}
						continue;
					}
				} while (muon_time_hit.size() > 0);
				//for (unsigned int i = 0; i < muon_clusters.size(); i++) cout << muon_clusters.at(i) << " ";
				//cout << endl;
				//cout << accumulate(muon_clusters.begin(), muon_clusters.end(), 0.0) / muon_clusters.size() << endl;
				muon_cluster_size_ALICE[point]->Fill(accumulate(muon_clusters.begin(), muon_clusters.end(), 0.0) / muon_clusters.size(),muon_clusters.size());
				muon_cluster_mult_ALICE[point]->Fill(muon_clusters.size());
				//if (muon_clusters.size() == 1) cout << "C.M. = 1" << endl;
				//double cs = accumulate(muon_clusters.begin(), muon_clusters.end(), 0.0) / muon_clusters.size();
				//if (cs != 8 || cs != 16) {
					//cout << "1= from 8 or 16" << endl;
				//	strp_profile_ALICE_muons_no8_no16[point]->FillN(muon_hits.size(),muon_hits[0],NULL);
				//}
				muon_clusters.clear();
				//cout << endl;
			}

			//////////////////////
			// Gamma Clustering //
			//////////////////////

			int gamma_pass = 0, gammma_tmpcluster = 1; //integer to keep track of the pass along the hits vector
			unsigned int gamma_counter_tmp = 1;
			coppia_strip.clear();
			coppia_time.clear();
			//if (gamma_hits.size() == 0) {
			//	muon_cluster_size_ALICE[point]->Fill(0);
			//}
			
			if (gamma_hits.size() == 1) {
			//else if (muon_hits.size() == 1) {
				gamma_cluster_size_ALICE[point]->Fill(1);
				gamma_cluster_mult_ALICE[point]->Fill(1);
			}

			else if (gamma_hits.size() > 1) {
				//cout << "Gamma size hit: " << gamma_hits.size() << endl;
				gamma_counter_tmp = 1;
				for (unsigned int hh = 0; hh < gamma_times.size(); hh++) {
					gamma_time_hit.push_back(make_pair(gamma_hits.at(hh),gamma_times.at(hh))); //hit - time
					//cout << "I'm here" << endl; //Debug
					//cout << endl; //Debug to be uncommented
				}

				//cout << "Time_hit size: " << muon_time_hit.size() << endl;
				
				sort(gamma_time_hit.begin(),gamma_time_hit.end());
			
				//for (unsigned int hh = 0; hh < gamma_time_hit.size(); hh++) {
				//	cout << gamma_time_hit.at(hh).first << "\t" << gamma_time_hit.at(hh).second << "\t";
				//}
				//cout << endl;

				do {
					//cout << "Size hit_time decreasing: " << gamma_time_hit.size() << endl;
					double first_strip = 0, first_time = 0, strip = 0, time = 0;
					if (gamma_pass == 0) { //Define the first element of the vector, compare all the others to that
						//cout << "Gamma pass: " << gamma_pass << endl;
						first_strip = gamma_time_hit.at(0).first;
						first_time = gamma_time_hit.at(0).second;
					}

					if (gamma_pass == 0 && gamma_time_hit.size() == 1) {
						gamma_clusters.push_back(1);
						gamma_time_hit.erase(gamma_time_hit.begin());
						gamma_counter_tmp = 1;
						continue;
					}

					//cout << "Assigning the other strips!" << endl;
					
					if (gamma_counter_tmp < gamma_time_hit.size()) {
						strip = gamma_time_hit.at(gamma_counter_tmp).first;
						time = gamma_time_hit.at(gamma_counter_tmp).second;						
					}

					if (gamma_pass == 0 && (TMath::Abs(strip-first_strip) == 1 || TMath::Abs(strip-first_strip) == 0) && TMath::Abs(time-first_time) <= clus_time && gamma_time_hit.size() == 2) {
						gamma_clusters.push_back(2);
						break;
					}

					if (gamma_pass == 0 && (TMath::Abs(strip-first_strip) == 1 || TMath::Abs(strip-first_strip) == 0) && TMath::Abs(time-first_time) <= clus_time) { //first pass on the array, compare first element to all other vector elements
						//cout << "Adjacent pair!" << endl;
						coppia_strip.push_back(first_strip);
						coppia_strip.push_back(strip);
						coppia_time.push_back(first_time);
						coppia_time.push_back(time);
						gamma_pass++;
						gammma_tmpcluster++;
						gamma_time_hit.erase(gamma_time_hit.begin()+gamma_counter_tmp);
						gamma_time_hit.erase(gamma_time_hit.begin());
						gamma_counter_tmp = 1;
						continue;
					}

					else if (gamma_pass > 0) { //There is a couple created
						for (unsigned int tmp = 0; tmp < coppia_strip.size(); tmp++) {
							double tmpstrip = coppia_strip.at(tmp);
							double tmptime = coppia_time.at(tmp);

							for (unsigned int tmp2 = 0; tmp2 < gamma_time_hit.size(); tmp2++) {
								double tmpstrip2 = gamma_time_hit.at(tmp2).first;
								double tmptime2 = gamma_time_hit.at(tmp2).second;

								if((TMath::Abs(tmpstrip2-tmpstrip) == 1 || TMath::Abs(tmpstrip2-tmpstrip) == 0) && TMath::Abs(tmptime2-tmptime) <= clus_time) {
									coppia_strip.push_back(tmpstrip2);
									coppia_time.push_back(tmptime2);
									gamma_time_hit.erase(gamma_time_hit.begin()+tmp2);
									gammma_tmpcluster++;
									gamma_counter_tmp = 1;
									continue;
								}
							}
						}
						gamma_clusters.push_back(gammma_tmpcluster);
						//cout << " " << gammma_tmpcluster << " ";
						gammma_tmpcluster = 1;
						coppia_strip.clear();
						coppia_time.clear();
						gamma_counter_tmp = 1;
						gamma_pass = 0;
						continue;
					}

					else if ((gamma_pass == 0 && TMath::Abs(strip-first_strip) > 1) || (gamma_pass == 0 && TMath::Abs(time-first_time) > clus_time)) {
						//cout << "Not adjacent pair on the first go!" << endl;
						gamma_counter_tmp++;
						if (gamma_counter_tmp == gamma_time_hit.size()) {
							gamma_time_hit.erase(gamma_time_hit.begin());
							gamma_clusters.push_back(1);
							gamma_counter_tmp = 1;
							continue;
						}
						continue;
					}
				} while (gamma_time_hit.size() > 0);
				//for (unsigned int i = 0; i < gamma_clusters.size(); i++) cout << gamma_clusters.at(i) << " ";
				//cout << endl;
				//cout << accumulate(gamma_clusters.begin(), gamma_clusters.end(), 0.0) / gamma_clusters.size() << endl;
				gamma_cluster_size_ALICE[point]->Fill(accumulate(gamma_clusters.begin(), gamma_clusters.end(), 0.0) / gamma_clusters.size(), gamma_clusters.size());
				gamma_cluster_mult_ALICE[point]->Fill(gamma_clusters.size());
				gamma_clusters.clear();
				//cout << endl;
			}
		} //End of cycle on the tree entries

		//Push back values for each HV point
		avg_muon_clus.push_back(muon_cluster_size_ALICE[point]->GetMean());
		e_avg_muon_clus.push_back(muon_cluster_size_ALICE[point]->GetMeanError());
		avg_muon_clus_mul.push_back(muon_cluster_mult_ALICE[point]->GetMean()); 
		e_avg_muon_clus_mul.push_back(muon_cluster_mult_ALICE[point]->GetMeanError());
		avg_gamma_clus.push_back(gamma_cluster_size_ALICE[point]->GetMean());
		e_avg_gamma_clus.push_back(gamma_cluster_size_ALICE[point]->GetMeanError());
		avg_gamma_clus_mul.push_back(gamma_cluster_mult_ALICE[point]->GetMean()); 
		e_avg_gamma_clus_mul.push_back(gamma_cluster_mult_ALICE[point]->GetMeanError());

		//Gamma rates for gamma/fake efficiency calculations
		tot_time = 5000*1e-9*trigger; //in s
		gamma_rate = gamma_counter/(gamma_window*trigger); //rate outside muon window
		muon_window_gamma = gamma_rate*muon_window*trigger; //number of gammas inside muon window
		p_gamma = TMath::Exp(-muon_window_gamma); 
		cout << point+1 << ") Gamma counter: " << gamma_counter << ", triggers: " << trigger << ", total time: " << tot_time << " s, gamma rate: " << gamma_rate << " Hz, avg gammas in mu window: " << muon_window_gamma << " p_gamma: " << p_gamma << endl; //Debug
		//Fake efficiency estimation, rescaled to the ratio of the muon window to the gamma window
		Gamma_eff = (count_Y_gamma/(double)trigger)*(muon_window/gamma_window);
		//Error on fake efficiency
		e_Gamma_eff = TMath::Sqrt(Gamma_eff*(1-Gamma_eff)/(double)trigger);
		gamma_counter = 0;
		muon_counter = 0;

		/////////////////
		//	PROFILES   //
		/////////////////	

		//Strip profile, muons
		c_strp_profile_ALICE_muons->cd();
		muon_strip->SetParLimits(0,0.,strp_profile_ALICE_muons[point]->GetMaximum());
		muon_strip->SetParLimits(1,strp_profile_ALICE_muons[point]->GetMean()-1,strp_profile_ALICE_muons[point]->GetMean()+1);
		//muon_strip->SetParLimits(2,0.,4.);
		strp_profile_ALICE_muons[point]->Fit("muon_strip","RM+");
		strp_profile_ALICE_muons[point]->Draw("HISTO");
		muon_strip->Draw("SAME");

		double mean_strip = muon_strip->GetParameter(1);
		double sigma_strip = muon_strip->GetParameter(2);
		//Muon efficiency without gamma contribution
		double raweff = (count_Y/(double)trigger);
		//Error on muon efficiency without gamma contribution
		double eraweff = TMath::Sqrt((raweff*(1-raweff))/trigger);
		//Efficiency corrected for gamma background inside muon window
		double gamma_eff = (raweff-Gamma_eff)/(1-Gamma_eff); 
		//Error on efficiency corrected for gammas
		double term_one = (1/(1-Gamma_eff))*(1/(1-Gamma_eff))*e_Gamma_eff*e_Gamma_eff;
		double term_two = ((raweff-1)/((1-Gamma_eff)*(1-Gamma_eff)))*((raweff-1)/((1-Gamma_eff)*(1-Gamma_eff)))*eraweff*eraweff;
		double er_gamma_eff = TMath::Sqrt(term_one + term_two);

		//Push back in vectors for plots
		eff.push_back(raweff*100);
		e_eff.push_back(eraweff*100);
		eff_gamma.push_back(gamma_eff*100);
		e_eff_gamma.push_back(er_gamma_eff*100);

		gammaefficiency.push_back(Gamma_eff*100); //Fake efficiency (gamma detection efficiency)
		e_gammaefficiency.push_back(e_Gamma_eff*100); //Error on fake efficiency (gamma detection efficiency)

		//Print efficiency
		cout << "HV point: " << point+1 << " RPC counts: " << count_Y << " Trigger num: " << trigger << " Raw efficiency: " << raweff << " Efficiency without gamma: " << gamma_eff << endl;
		cout << "Gamma eff: " << Gamma_eff << " +- " << e_Gamma_eff << endl;

		//Strip profile, gammas
		c_strp_profile_ALICE_gamma->cd();
		strp_profile_ALICE_gamma[point]->Draw("HISTO");

		//Gamma rate (in the muon window) + average rate calculations
		for (int strip = 1; strip <= 16; strip++) {
			double gamma_hits = strp_profile_ALICE_gamma[point]->GetBinContent(strip);
			//During a spill, utside of mu window 
			//strp_profile_ALICE_Y_gamma_rate[point]->SetBinContent(strip,gamma_hits/(strip_area*(gamma_window)*trigger)); //gamma_window*trigger = total time of the measurements
			//Outside spills, no mu
			strp_profile_ALICE_gamma_rate[point]->SetBinContent(strip,gamma_hits/(strip_area*(4900*1e-9)*gamma_trigger)); //meas_time*trigger = total time of the measurements
			rate_bin = strp_profile_ALICE_gamma_rate[point]->GetBinContent(strip);
			if (strip >= 9 && strip <= 16) rate_sum += rate_bin; 
			//cout << "strip: " << strip << ", rate: " << rate_bin << ", rate sum:" << rate_sum << endl; 
		}

		//calculate the mean + error
		double avg_rate = 0;

		//avg_rate = rate_sum/16;
		avg_rate = rate_sum/8;
		//for (int i = 1; i <= 16; i++) { //All strips
		for (int i = 9; i <= 16; i++) { //Only strips in second board
			rate_bin = strp_profile_ALICE_gamma_rate[point]->GetBinContent(i);
			partial += (rate_bin-avg_rate)*(rate_bin-avg_rate);
		}
		//All strips
		//rate_dev = TMath::Sqrt(partial/15);
		//e_avg_rate = rate_dev/TMath::Sqrt(16);

		//Only strips in second board
		rate_dev = TMath::Sqrt(partial/7);
		e_avg_rate = rate_dev/TMath::Sqrt(8);

		partial = 0; //Reset variable to calculate std deviation and then error on the mean
		avg_gamma_rate.push_back(avg_rate);
		e_avg_gamma_rate.push_back(e_avg_rate);

		//cout << "Rate: (" << avg_rate << " +- " << e_avg_rate << ") Hz" << endl;

		TPaveText *trate = new TPaveText(0.111051,0.913306,0.225894,0.960685,"NDC");
		trate->AddText(Form("Gama rate: %3.1lf  HZ/cm^{2}",avg_rate));

		//Strip profile, gamma rate
		c_strp_profile_ALICE_gamma_rate->cd();
		strp_profile_ALICE_gamma_rate[point]->Draw("HISTO");
		trate->Draw("SAME");

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

		printf("T [K]: %3f, p [mbar]: %3f, HV mon [V]: %3f, HV eff [V]: %3f \n",T,p,tens_mon,tens_eff); //Debug printout

		//Muon cluster size
		c_cluster_size_muon_ALICE->cd();
		muon_cluster_size_ALICE[point]->Scale(1/muon_cluster_size_ALICE[point]->Integral());
		muon_cluster_size_ALICE[point]->Draw("HISTO");

		//Muon cluster multiplicity
		c_cluster_mult_muon_ALICE->cd();
		muon_cluster_mult_ALICE[point]->Scale(1/muon_cluster_mult_ALICE[point]->Integral());
		muon_cluster_mult_ALICE[point]->Draw("HISTO");

		//Gamma cluster size
		c_cluster_size_gamma_ALICE->cd();
		gamma_cluster_size_ALICE[point]->Scale(1/gamma_cluster_size_ALICE[point]->Integral());
		gamma_cluster_size_ALICE[point]->Draw("HISTO");

		//Gamma cluster multiplicity
		c_cluster_mult_gamma_ALICE->cd();
		gamma_cluster_mult_ALICE[point]->Scale(1/gamma_cluster_mult_ALICE[point]->Integral());
		gamma_cluster_mult_ALICE[point]->Draw("HISTO");

		//Muon time delta in cluster
		c_deltaTime->cd();
		deltaTime[point]->Draw("HISTO");
		after15[point] = new TLine(clus_time,0,clus_time,deltaTime[point]->GetMaximum());
		after15[point]->SetLineColor(kRed);
		int countsAfter15 = deltaTime[point]->Integral(148,500);
		int totalCounts = deltaTime[point]->Integral();
		double percentageOut = countsAfter15/(double)totalCounts;
		cout << endl << "Counts after 15 ns: " << countsAfter15 << endl;
		cout << endl << "Counts total: " << totalCounts << endl;
		cout << endl << "% 15 ns: " << percentageOut*100 << endl;
		tdelta[point]->AddText(Form("Counts after 15 ns: %3.1d",countsAfter15));
		tdelta[point]->AddText(Form("Total counts: %3.1d",totalCounts));
		tdelta[point]->AddText(Form("After 15 ns: %3.1lf %%",percentageOut*100));
		tdelta[point]->Draw("SAME");
		after15[point]->Draw("SAME");

		fout->cd(); //open .root file
		cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
		c_time_profile_ALICE->Write(("Time_profile_ALICE_"+plane+"_"+cartella).c_str()); //Save time profile
		c_strp_profile_ALICE->Write(("Strip_profile_ALICE_all_hits_"+plane+"_"+cartella).c_str()); //Save strip profile (all hits)
		c_strp_profile_ALICE_muons->Write(("Strip_profile_ALICE_muons_"+plane+"_"+cartella).c_str());
		c_strp_profile_ALICE_gamma->Write(("Strip_profile_ALICE_gammas_"+plane+"_"+cartella).c_str());
		c_strp_profile_ALICE_gamma_rate->Write(("Strip_profile_ALICE_gammas_rate_"+plane+"_"+cartella).c_str());
		c_cluster_size_muon_ALICE->Write(("Muon_cluster_Size_ALICE_"+plane+"_"+cartella).c_str());
		c_cluster_mult_muon_ALICE->Write(("Muon_cluster_Mult_ALICE_"+plane+"_"+cartella).c_str());
		c_cluster_size_gamma_ALICE->Write(("Gamma_cluster_Size_ALICE_"+plane+"_"+cartella).c_str());
		c_cluster_mult_gamma_ALICE->Write(("Gamma_cluster_Mult_ALICE_"+plane+"_"+cartella).c_str());
		c_deltaTime->Write(("Delta_time_cluster_ALICE_"+plane+"_"+cartella).c_str());
		//c_time_profile_strips[point]->Write(("Strip_time_profile_ALICE_"+plane+"_"+cartella).c_str());
	} //End of cycle on the HV points in a scan

	if (plane == "X") ftimesX.close(); //close file where the timestamps are saved
	else if (plane == "Y") ftimesY.close(); //close file where the timestamps are saved

	//Set parameters of the functions used for the fit of eff curve
	fit_eff->SetParLimits(0,50.,110.);
	fit_eff->SetParLimits(2,9000.,10000.);
	fit_eff_gamma->SetParLimits(0,50.,110.);
	fit_eff_gamma->SetParLimits(2,9000.,10000.);

	//Raw efficiency (with gamma background contribution)
	TCanvas *ceff = new TCanvas();
	ceff->cd();

	TGraphErrors *RawEff = new TGraphErrors(veff_eco.size(),&veff_eco[0],&eff[0],&eveff[0],&e_eff[0]);
	RawEff->SetTitle(("Eff curve for scan " + plane + " " + to_string(scan)).c_str());
	RawEff->GetXaxis()->SetTitle("HV_{eff} [V]");
	RawEff->GetYaxis()->SetTitle("#epsilon");
	RawEff->GetYaxis()->SetRangeUser(veff_eco.at(0),veff_eco.at(veff_eco.size()-1));
	RawEff->GetYaxis()->SetRangeUser(0.,100.);
	RawEff->SetMarkerStyle(8);
	RawEff->SetMarkerSize(1);
	RawEff->Fit(fit_eff,"RM+");
	RawEff->Draw("AP");
	fit_eff->Draw("SAME");

	//Fit parameters
	double effmax = fit_eff->GetParameter(0);
	double lambda = fit_eff->GetParameter(1);
	double HV50 = fit_eff->GetParameter(2);

	//double WP = (TMath::Log(19)/lambda)+HV50+150; //CMS definition
	double WP = fit_eff->GetX(0.95*effmax,veff_eco.front(),veff_eco.back()); //Definition for ecogas paper = 95% of maximum efficiency
	double effwp = fit_eff->Eval(WP);

	cout << "WP: " << WP << "\t eff WP: " << effwp << endl;

	TPaveText *teff = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	teff->AddText(Form("Eff max: %3.1lf",effmax));
	teff->AddText(Form("#lambda: %3.4lf",lambda));
	teff->AddText(Form("HV 50%% (V): %3.1lf",HV50));
	teff->AddText(Form("WP (V): %3.1lf",WP)); 
	teff->AddText(Form("Eff WP: %3.1lf",effwp));
	teff->Draw("SAME");

	//Efficiency without gamma background contribution
	TCanvas *ceffgamma = new TCanvas();
	ceffgamma->cd();

	TGraphErrors *GammaEff = new TGraphErrors(veff_eco.size(),&veff_eco[0],&eff_gamma[0],&eveff[0],&e_eff[0]);
	GammaEff->SetTitle(("Eff curve for scan, gamma corrected " + plane + " " + to_string(scan)).c_str());
	GammaEff->GetXaxis()->SetTitle("HV_{eff} [V]");
	GammaEff->GetYaxis()->SetTitle("#epsilon");
	GammaEff->GetYaxis()->SetRangeUser(veff_eco.at(0),veff_eco.at(veff_eco.size()-1));
	GammaEff->GetYaxis()->SetRangeUser(0.,100.);
	GammaEff->SetMarkerStyle(8);
	GammaEff->SetMarkerSize(1);
	GammaEff->Fit(fit_eff_gamma,"RM+");
	GammaEff->Draw("AP");
	fit_eff_gamma->Draw("SAME");

	//Fit parameters
	double effmax_gamma = fit_eff_gamma->GetParameter(0);
	double lambda_gamma = fit_eff_gamma->GetParameter(1);
	double HV50_gamma = fit_eff_gamma->GetParameter(2);

	//double WP_gamma = (TMath::Log(19)/lambda_gamma)+HV50_gamma+150; //CMS definition
	double WP_gamma = fit_eff_gamma->GetX(0.95*effmax_gamma,veff_eco.front(),veff_eco.back()); //Definition for ecogas paper = 95% of maximum efficiency
	double effwp_gamma = fit_eff_gamma->Eval(WP_gamma);

	TPaveText *teffgamma = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	teffgamma->AddText(Form("Eff max: %3.1lf",effmax_gamma));
	teffgamma->AddText(Form("#lambda: %3.4lf",lambda_gamma));
	teffgamma->AddText(Form("HV 50%% (V): %3.1lf",HV50_gamma));
	teffgamma->AddText(Form("WP (V): %3.1lf",WP_gamma)); 
	teffgamma->AddText(Form("Eff WP: %3.1lf",effwp_gamma));
	teffgamma->Draw("SAME");

	fout->cd();
	ceff->Write(("Eff_curve_" + plane + " " + to_string(scan)).c_str());
	ceffgamma->Write(("Eff_curve_gamma_corrected_" + plane + " " + to_string(scan)).c_str());

	//Fake efficiency (Gamma detection efficiency)
	TCanvas *cfake = new TCanvas();
	cfake->cd();
	TGraphErrors *gfake = new TGraphErrors(veff_eco.size(),&veff_eco[0],&gammaefficiency[0],&eveff[0],&e_gammaefficiency[0]);
	gfake->SetTitle(("Fake eff curve for scan " + plane + " " + to_string(scan)).c_str());
	gfake->GetXaxis()->SetTitle("HV_{eff} [V]");
	gfake->GetYaxis()->SetTitle("Fake/Gamma efficiency [%]");
	gfake->GetYaxis()->SetRangeUser(veff_eco.at(0),veff_eco.at(veff_eco.size()-1));
	gfake->GetYaxis()->SetRangeUser(0.,5.);
	gfake->SetMarkerStyle(8);
	gfake->SetMarkerSize(1);
	gfake->Draw("AP");

	fout->cd();
	cfake->Write(("Gamma_fake_eff_curve_" + plane + " " + to_string(scan)).c_str());

	//I(V) curve 
	TCanvas *civ = new TCanvas();
	civ->cd();
	TGraphErrors *iv = new TGraphErrors(veff_eco.size(),&veff_eco[0],&imon[0],&eveff[0],&eimon[0]);
	iv->SetTitle(("I(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	iv->GetXaxis()->SetTitle("HV_{eff} [V]");
	iv->GetYaxis()->SetTitle("I_{mon} [#muA]");
	iv->GetYaxis()->SetRangeUser(0,imon.back()+0.05*imon.back());
	iv->SetMarkerStyle(8);
	iv->SetMarkerSize(1);
	iv->SetMarkerColor(kRed);
	iv->Draw("AP");

	double i_wp = iv->Eval(WP);
	TPaveText *tiv = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tiv->AddText(Form("I (WP) [#muA]: %3.1lf",i_wp));
	tiv->Draw("SAME");

	fout->cd();
	civ->Write(("I(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Calculate current density
	for (unsigned int i = 0; i < imon.size(); i++) {
		currdens.push_back(imon.at(i)/area);
		ecurrdens.push_back(eimon.at(i)/area);
	}

	//Current density(V) curve 
	TCanvas *civdens = new TCanvas();
	civdens->cd();
	TGraphErrors *ivdens = new TGraphErrors(veff_eco.size(),&veff_eco[0],&currdens[0],&eveff[0],&ecurrdens[0]);
	ivdens->SetTitle(("Current density(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	ivdens->GetXaxis()->SetTitle("HV_{eff} [V]");
	ivdens->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	ivdens->GetYaxis()->SetRangeUser(0,currdens.back()+0.05*currdens.back());
	ivdens->SetMarkerStyle(8);
	ivdens->SetMarkerSize(1);
	ivdens->SetMarkerColor(kRed);
	ivdens->Draw("AP");

	double i_dens_wp = ivdens->Eval(WP);
	TPaveText *tivdens = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tivdens->AddText(Form("Curr dens (WP) [#muA/cm^{2}]: %3.1lf",i_dens_wp));
	tivdens->Draw("SAME");

	fout->cd();
	civdens->Write(("Current_density(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Rate(HV) curve
	TCanvas *crate = new TCanvas();
	crate->cd();
	TGraphErrors *ratev = new TGraphErrors(veff_eco.size(),&veff_eco[0],&avg_gamma_rate[0],&eveff[0],&e_avg_gamma_rate[0]);
	ratev->SetTitle(("Rate(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	ratev->GetXaxis()->SetTitle("HV_{eff} [V]");
	ratev->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	ratev->GetYaxis()->SetRangeUser(0,avg_gamma_rate.back()+0.05*avg_gamma_rate.back());
	ratev->SetMarkerStyle(8);
	ratev->SetMarkerSize(1);
	ratev->SetMarkerColor(kRed);
	ratev->Draw("AP");

	double rate_wp = ratev->Eval(WP);
	TPaveText *trate = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	trate->AddText(Form("Rate (WP) [Hz/cm^{2}]: %3.1lf",rate_wp));
	trate->Draw("SAME");

	fout->cd();
	crate->Write(("Rate(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Clutster rate (HV) curve
	vector<double> vclusrate, evclusrate;
	for (unsigned int i = 0; i < avg_gamma_rate.size(); i++) {
		vclusrate.push_back(avg_gamma_rate.at(i)/avg_gamma_clus.at(i));
		evclusrate.push_back((1/avg_gamma_clus.at(i))*TMath::Sqrt(e_avg_gamma_rate.at(i)*e_avg_gamma_rate.at(i)+((avg_gamma_rate.at(i)*avg_gamma_rate.at(i))/(avg_gamma_clus.at(i)*avg_gamma_clus.at(i)))*e_avg_gamma_clus.at(i)*e_avg_gamma_clus.at(i)));
	}

	TCanvas *cclusrate = new TCanvas();
	cclusrate->cd();
	TGraphErrors *clusratev = new TGraphErrors(veff_eco.size(),&veff_eco[0],&vclusrate[0],&eveff[0],&evclusrate[0]);
	clusratev->SetTitle(("Cluster Rate(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	clusratev->GetXaxis()->SetTitle("HV_{eff} [V]");
	clusratev->GetYaxis()->SetTitle("Gamma cluster rate [Hz/cm^{2}]");
	clusratev->GetYaxis()->SetRangeUser(0,vclusrate.back()+0.05*vclusrate.back());
	clusratev->SetMarkerStyle(8);
	clusratev->SetMarkerSize(1);
	clusratev->SetMarkerColor(kRed);
	clusratev->Draw("AP");

	double cluster_rate_wp = clusratev->Eval(WP);
	TPaveText *tcrate = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tcrate->AddText(Form("Cluster rate (WP) [Hz/cm^{2}]: %3.1lf",cluster_rate_wp));
	tcrate->Draw("SAME");

	fout->cd();
	cclusrate->Write(("Cluster_rate(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Muon cluster size (HV) curve
	TCanvas *cmuonclustersize = new TCanvas();
	cmuonclustersize->cd();
	TGraphErrors *muonclutsersizev = new TGraphErrors(veff_eco.size(),&veff_eco[0],&avg_muon_clus[0],&eveff[0],&e_avg_muon_clus[0]);
	muonclutsersizev->SetTitle(("Muon_cluster_size(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	muonclutsersizev->GetXaxis()->SetTitle("HV_{eff} [V]");
	muonclutsersizev->GetYaxis()->SetTitle("Muon cluster size [strips]");
	muonclutsersizev->GetYaxis()->SetRangeUser(0,avg_muon_clus.back()+0.05*avg_muon_clus.back());
	muonclutsersizev->SetMarkerStyle(8);
	muonclutsersizev->SetMarkerSize(1);
	muonclutsersizev->SetMarkerColor(kBlue);
	muonclutsersizev->Draw("AP");

	double muon_cluster_size_wp = muonclutsersizev->Eval(WP);
	TPaveText *tmuonclustersize = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tmuonclustersize->AddText(Form("Muon c.s. (WP) [strips]: %3.1lf",muon_cluster_size_wp));
	tmuonclustersize->Draw("SAME");

	fout->cd();
	cmuonclustersize->Write(("Muon_cluster_size(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Muon cluster multiplicity (HV) curve
	TCanvas *cmuonclustermult = new TCanvas();
	cmuonclustermult->cd();
	TGraphErrors *muonclutsermultv = new TGraphErrors(veff_eco.size(),&veff_eco[0],&avg_muon_clus_mul[0],&eveff[0],&e_avg_muon_clus_mul[0]);
	muonclutsermultv->SetTitle(("Muon_cluster_mult(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	muonclutsermultv->GetXaxis()->SetTitle("HV_{eff} [V]");
	muonclutsermultv->GetYaxis()->SetTitle("Muon cluster multiplicity");
	muonclutsermultv->GetYaxis()->SetRangeUser(0,avg_muon_clus_mul.back()+0.05*avg_muon_clus_mul.back());
	muonclutsermultv->SetMarkerStyle(8);
	muonclutsermultv->SetMarkerSize(1);
	muonclutsermultv->SetMarkerColor(kBlue);
	muonclutsermultv->Draw("AP");

	double muon_cluster_mult_wp = muonclutsermultv->Eval(WP);
	TPaveText *tmuonclustermult = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tmuonclustermult->AddText(Form("Muon cmp (WP): %3.1lf",muon_cluster_mult_wp));
	tmuonclustermult->Draw("SAME");

	fout->cd();
	cmuonclustermult->Write(("Muon_cluster_mult(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Gamma cluster size (HV) curve
	TCanvas *cgammaclustersize = new TCanvas();
	cgammaclustersize->cd();
	TGraphErrors *gammaclutsersizev = new TGraphErrors(veff_eco.size(),&veff_eco[0],&avg_gamma_clus[0],&eveff[0],&e_avg_gamma_clus[0]);
	gammaclutsersizev->SetTitle(("Gamma_cluster_size(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	gammaclutsersizev->GetXaxis()->SetTitle("HV_{eff} [V]");
	gammaclutsersizev->GetYaxis()->SetTitle("Gamma cluster size [strips]");
	gammaclutsersizev->GetYaxis()->SetRangeUser(0,avg_gamma_clus.back()+0.05*avg_gamma_clus.back());
	gammaclutsersizev->SetMarkerStyle(8);
	gammaclutsersizev->SetMarkerSize(1);
	gammaclutsersizev->SetMarkerColor(kRed);
	gammaclutsersizev->Draw("AP");

	double gamma_cluster_size_wp = gammaclutsersizev->Eval(WP);
	TPaveText *tgammaclustersize = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tgammaclustersize->AddText(Form("Gamma c.s. (WP) [strips]: %3.1lf",gamma_cluster_size_wp));
	tgammaclustersize->Draw("SAME");

	fout->cd();
	cgammaclustersize->Write(("Gamma_cluster_size(V)_curve_" + plane + " " + to_string(scan)).c_str());

	//Gamma cluster multiplicity (HV) curve
	TCanvas *cgammaclustermult = new TCanvas();
	cgammaclustermult->cd();
	TGraphErrors *gammaclutsermultv = new TGraphErrors(veff_eco.size(),&veff_eco[0],&avg_gamma_clus_mul[0],&eveff[0],&e_avg_gamma_clus_mul[0]);
	gammaclutsermultv->SetTitle(("Gamma_cluster_mult(HV_{eff}) " + plane + " " + to_string(scan)).c_str());
	gammaclutsermultv->GetXaxis()->SetTitle("HV_{eff} [V]");
	gammaclutsermultv->GetYaxis()->SetTitle("Gamma cluster multiplicity");
	gammaclutsermultv->GetYaxis()->SetRangeUser(0,avg_gamma_clus_mul.back()+0.05*avg_gamma_clus_mul.back());
	gammaclutsermultv->SetMarkerStyle(8);
	gammaclutsermultv->SetMarkerSize(1);
	gammaclutsermultv->SetMarkerColor(kRed);
	gammaclutsermultv->Draw("AP");

	double gamma_cluster_mult_wp = gammaclutsermultv->Eval(WP);
	TPaveText *tgammaclustermult = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
	tgammaclustermult->AddText(Form("Gamma cmp (WP): %3.1lf",gamma_cluster_mult_wp));
	tgammaclustermult->Draw("SAME");

	fout->cd();
	cgammaclustermult->Write(("Gamma_cluster_mult(V)_curve_" + plane + " " + to_string(scan)).c_str());

	fout->Close();

	gSystem->cd("/home/luca/Desktop/PhD/TB_october_2021/runs/std/Currents");

	//Write HV and I to a txt file
	ofstream currfile;
	currfile.open(("I(HV)_ALICE_" + plane + to_string(scan) + ".txt").c_str());
	currfile << "HV eff [V]" << "\t" << "I [muA]" << "\t" << "Efficiency [%]" << "Eff error [%]" << "\t" << "Cluster rate [Hz/cm^{2}]" << "\t" << "Err Clus rate" << endl;
	for (unsigned int icurr = 0; icurr < veff_eco.size(); icurr++) {
		currfile << veff_eco.at(icurr) << "\t" << imon.at(icurr) << "\t" << eff_gamma.at(icurr) << "\t" << e_eff.at(icurr) << "\t" << vclusrate.at(icurr) << "\t" << evclusrate.at(icurr) << endl;
	}

	currfile.close();

}