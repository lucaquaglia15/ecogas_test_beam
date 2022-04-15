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
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

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

void plot_window(const int scanStd, const int scanEco2, const int scanEco3) {

	//STD gas mixture
	gSystem->cd("/home/luca/Desktop/PhD/TB_october_2021/runs/std");
	
	//X Timestamp file
	ifstream timestampXstd;
	timestampXstd.open(("Timestamps_X_scan"+to_string(scanStd)+".txt").c_str());
	double startXstd, meanXstd, endXstd, windowXstd;
	vector <double> vstartXstd, vmeanXstd, vendXstd, vwindowXstd;

	while (timestampXstd >> startXstd >> meanXstd >> endXstd >> windowXstd) {
		vstartXstd.push_back(startXstd);
		vmeanXstd.push_back(meanXstd);
		vendXstd.push_back(endXstd);
		vwindowXstd.push_back(windowXstd);
	}
	
	//Y Timestamp file
	ifstream timestampYstd;
	timestampYstd.open(("Timestamps_Y_scan"+to_string(scanStd)+".txt").c_str());
	double startYstd, meanYstd, endYstd, windowYstd;
	vector <double> vstartYstd, vmeanYstd, vendYstd, vwindowYstd;

	while (timestampYstd >> startYstd >> meanYstd >> endYstd >> windowYstd) {
		vstartYstd.push_back(startYstd);
		vmeanYstd.push_back(meanYstd);
		vendYstd.push_back(endYstd);
		vwindowYstd.push_back(windowYstd);
	}

	//ECO2 gas mixture
	gSystem->cd("/home/luca/Desktop/PhD/TB_october_2021/runs/eco2");
	
	//X Timestamp file
	ifstream timestampXeco2;
	timestampXeco2.open(("Timestamps_X_scan"+to_string(scanEco2)+".txt").c_str());
	double startXeco2, meanXeco2, endXeco2, windowXeco2;
	vector <double> vstartXeco2, vmeanXeco2, vendXeco2, vwindowXeco2;

	while (timestampXeco2 >> startXeco2 >> meanXeco2 >> endXeco2 >> windowXeco2) {
		vstartXeco2.push_back(startXeco2);
		vmeanXeco2.push_back(meanXeco2);
		vendXeco2.push_back(endXeco2);
		vwindowXeco2.push_back(windowXeco2);
	}
	
	//Y Timestamp file
	ifstream timestampYeco2;
	timestampYeco2.open(("Timestamps_Y_scan"+to_string(scanEco2)+".txt").c_str());
	double startYeco2, meanYeco2, endYeco2, windowYeco2;
	vector <double> vstartYeco2, vmeanYeco2, vendYeco2, vwindowYeco2;

	while (timestampYeco2 >> startYeco2 >> meanYeco2 >> endYeco2 >> windowYeco2) {
		vstartYeco2.push_back(startYeco2);
		vmeanYeco2.push_back(meanYeco2);
		vendYeco2.push_back(endYeco2);
		vwindowYeco2.push_back(windowYeco2);
	}

	//ECO3 gas mixture
	gSystem->cd("/home/luca/Desktop/PhD/TB_october_2021/runs/eco3");
	
	//X Timestamp file
	ifstream timestampXeco3;
	timestampXeco3.open(("Timestamps_X_scan"+to_string(scanEco3)+".txt").c_str());
	double startXeco3, meanXeco3, endXeco3, windowXeco3;
	vector <double> vstartXeco3, vmeanXeco3, vendXeco3, vwindowXeco3;

	while (timestampXeco3 >> startXeco3 >> meanXeco3 >> endXeco3 >> windowXeco3) {
		vstartXeco3.push_back(startXeco3);
		vmeanXeco3.push_back(meanXeco3);
		vendXeco3.push_back(endXeco3);
		vwindowXeco3.push_back(windowXeco3);
	}
	
	//Y Timestamp file
	ifstream timestampYeco3;
	timestampYeco3.open(("Timestamps_Y_scan"+to_string(scanEco3)+".txt").c_str());
	double startYeco3, meanYeco3, endYeco3, windowYeco3;
	vector <double> vstartYeco3, vmeanYeco3, vendYeco3, vwindowYeco3;

	while (timestampYeco3 >> startYeco3 >> meanYeco3 >> endYeco3 >> windowYeco3) {
		vstartYeco3.push_back(startYeco3);
		vmeanYeco3.push_back(meanYeco3);
		vendYeco3.push_back(endYeco3);
		vwindowYeco3.push_back(windowYeco3);
	}

	TMultiGraph *m = new TMultiGraph();

	vector<double> stdnum = {1,2};
	vector<double> std_mean = {vmeanXstd.at(6), vmeanYstd.at(6)};
	vector<double> std_low = {vmeanXstd.at(6)-vstartXstd.at(6), vmeanXstd.at(6)-vstartYstd.at(6)};
	vector<double> std_high = {vendXstd.at(6)-vmeanXstd.at(6), vendYstd.at(6)-vmeanXstd.at(6)};
	TGraphAsymmErrors *std = new TGraphAsymmErrors(stdnum.size(),&stdnum[0],&std_mean[0],NULL,NULL,&std_low[0],&std_high[0]);
	std->SetMarkerStyle(8);
	std->SetMarkerSize(1);
	std->SetMarkerColor(kBlack);
	std->SetLineColor(kBlack);
	m->Add(std);

	vector<double> eco2num = {3,4};
	vector<double> eco2_mean = {vmeanXeco2.at(6), vmeanYeco2.at(6)};
	vector<double> eco2_low = {vmeanXeco2.at(6)-vstartXeco2.at(6), vmeanXeco2.at(6)-vstartYeco2.at(6)};
	vector<double> eco2_high = {vendXeco2.at(6)-vmeanXeco2.at(6), vendYeco2.at(6)-vmeanXeco2.at(6)};
	TGraphAsymmErrors *eco2 = new TGraphAsymmErrors(eco2num.size(),&eco2num[0],&eco2_mean[0],NULL,NULL,&eco2_low[0],&eco2_high[0]);
	eco2->SetMarkerStyle(8);
	eco2->SetMarkerSize(1);
	eco2->SetMarkerColor(kRed);
	eco2->SetLineColor(kRed);
	m->Add(eco2);

	vector<double> eco3num = {5,6};
	vector<double> eco3_mean = {vmeanXeco3.at(6), vmeanYeco3.at(6)};
	vector<double> eco3_low = {vmeanXeco3.at(6)-vstartXeco3.at(6), vmeanXeco3.at(6)-vstartYeco3.at(6)};
	vector<double> eco3_high = {vendXeco3.at(6)-vmeanXeco3.at(6), vendYeco3.at(6)-vmeanXeco3.at(6)};
	TGraphAsymmErrors *eco3 = new TGraphAsymmErrors(eco3num.size(),&eco3num[0],&eco3_mean[0],NULL,NULL,&eco3_low[0],&eco3_high[0]);
	eco3->SetMarkerStyle(8);
	eco3->SetMarkerSize(1);
	eco3->SetMarkerColor(kGreen);
	eco3->SetLineColor(kGreen);
	m->Add(eco3);

	for (unsigned int i = 0; i < eco2num.size(); i++) {
		cout << "std: " << std_low.at(i) << "\t" << std_mean.at(i) << "\t" << std_high.at(i) << endl;
		cout << "eco2: " << eco2_low.at(i) << "\t" << eco2_mean.at(i) << "\t" << eco2_high.at(i) << endl;
		cout << "eco3: " << eco3_low.at(i) << "\t" << eco3_mean.at(i) << "\t" << eco3_high.at(i) << endl;
	}


	new TCanvas();
	m->SetTitle("Muon window for different gas mixtures");
	m->GetXaxis()->SetTitle("Mixture");
	m->GetYaxis()->SetTitle("Muon window [ns]");
	//m->GetYaxis()->SetRangeUser(0,7);
	//m->GetYaxis()->SetRangeUser(4000.,5000.);
	m->Draw("AP");

}