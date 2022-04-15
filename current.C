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
#include "TSpectrum.h"
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

bool verbose = true; //If true comments will be printed, otherwise no comments are printed

//TF1 *peakFun = new TF1("peakFun","gaus(0)+gaus(3)",0,100);
//TF1 *peakFun = new TF1("peakFun","gaus(0)",0,100);

//TF1 *peakFun[2];
// = new TF1("peakFun","gaus(0)",0,100);

void current(int scan) {

	double area = 2500; //cm^2 = RPC area

	string fol = "/media/luca/Elements/TB_October_2021/STD/Scan_00" + to_string(scan);

	//Create .root file with all the output data
	TFile *fout = new TFile(("current_ALICE"+to_string(scan)+".root").c_str(),"RECREATE");

	fout->cd();
	TDirectory *cdtof[12]; //12 directories in the root file, one for each HV point

	gSystem->cd(fol.c_str()); //Enter folder

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

    if (verbose) cout << "Number of HV points: " << file_count/3 << endl;

    double sxCurr = 0, dxCurr = 0;
    double avgCurr = 0, eAvgCurr = 0;
    double hv = 0, eHv = 0;
    vector<double> curr, eCurr;
    vector<double> voltage, eVoltage;
    int binSx = 0, binDx = 0;

	//Read DAQ root files: cycle on all 12 HV points --- MAIN part of the code
	for (int point = 0; point < file_count/3; point++) {

		if (verbose) cout << "HV point: " << point+1 << endl;
		
		string cartella = "HV" + to_string(point+1); //Create folders for different HV points in the out root file
		cdtof[point] = fout->mkdir(cartella.c_str());  

		string caen = "Scan00" +to_string(scan) + "_HV" + to_string(point+1) + "_CAEN.root"; //For HV and I

		//Get current and HV values
		TFile *f1 = new TFile(caen.c_str(),"READ");
		TH1F *HVeff = (TH1F*)f1->Get("HVeff_ALICE-2-0");
		TH1F *Imon = (TH1F*)f1->Get("Imon_ALICE-2-0");
		Imon->Rebin(2); //Rebin to find the peaks better -> this works
		//Imon->Rebin(8); //Rebin to find the peaks better
		
		//TH1F *ImonNoZoom = (TH1F*)f1->Get("Imon_ALICE-2-0");

		TCanvas *c = new TCanvas();

		hv = HVeff->GetMean();
		eHv = HVeff->GetMeanError();

		voltage.push_back(hv);
		eVoltage.push_back(eHv);

		int nonEmptyCounter = 0; //Check how many bins are not empty in the current histogram

		for (int i = 0; i < Imon->GetNbinsX(); i++) {
			if (Imon->GetBinContent(i) != 0) nonEmptyCounter++;
		}

		if (verbose) cout << "Non empty: " << nonEmptyCounter << endl;


		if (nonEmptyCounter == 1 || nonEmptyCounter == 2) { //Only one bin in the histogram, no need to fit just take the mean and its errors for the current

			avgCurr = Imon->GetMean();
   			eAvgCurr = 0.1; //if only one value is measured -> the error is the sensibility of the HV module

   			if (verbose) cout << "Avg curr: " << avgCurr << "+-" << eAvgCurr << "\t HV: " << hv << "+-" << eHv << endl; 

   			curr.push_back(avgCurr);
   			eCurr.push_back(eAvgCurr);

   			c->cd();
   			Imon->Draw("HISTO");
		}

		else if (nonEmptyCounter > 2) {

			double maxCounts = Imon->GetMaximum();  //Find the maximum number of counts to use as an estimate for the gaussian constant
			int maxBin = Imon->GetMaximumBin(); //Find bin with the most counts 
			double maxValue = Imon->GetBinCenter(maxBin); //Find x value with the most counts in order to use it as an estimate for the gaussian mean

			if (verbose) {
				cout << endl << maxCounts << endl;
				cout << endl << maxValue << endl;
			}
			
			TF1 *peakFun = new TF1("peakFun","gaus(0)",maxValue-4,maxValue+4);
			peakFun->SetNpx(1000);
			peakFun->SetParLimits(0,maxCounts-5,maxCounts+5);
			peakFun->SetParLimits(1,maxValue,maxValue+0.01);
			//peakFun->SetParameter(1,maxValue+0.1);

			Imon->Fit("peakFun","0M+");

			string minuitstatus = string(gMinuit->fCstatu);
			cout << endl << endl << "Status: " << minuitstatus << endl << endl;
			//int r = 0;
			if(minuitstatus.compare("CONVERGED ") != 0 && minuitstatus.compare("OK        ") != 0) {
				cout << endl << endl << "Fit failed!! Trying again!" << endl << endl;
				peakFun->SetParLimits(0,maxCounts-5,maxCounts+5);
				peakFun->SetParLimits(1,maxValue-0.01,maxValue+0.05);
				Imon->Fit("peakFun","M+");
			}
   		
   			sxCurr = peakFun->GetParameter(1) - 5*peakFun->GetParameter(2);
   			dxCurr = peakFun->GetParameter(1) + 5*peakFun->GetParameter(2);
   			if (verbose) cout << "Current range: " << sxCurr << "----" << dxCurr << endl;

   			binSx = Imon->FindBin(sxCurr);
   			binDx = Imon->FindBin(dxCurr);
   			if (verbose) cout << "Bin range: " << binSx << "----" << binDx << endl;

   			//Imon->GetXaxis()->SetRange(binSx,binDx); //This works
  			Imon->GetXaxis()->SetRangeUser(sxCurr,dxCurr);

   			avgCurr = Imon->GetMean();
   			eAvgCurr = Imon->GetMeanError();
   			if (eAvgCurr > 0.1) eCurr.push_back(eAvgCurr);
   			else if (eAvgCurr <= 0.1) eCurr.push_back(0.1);

   			if (verbose) cout << "Avg curr: " << avgCurr << "+-" << eAvgCurr << "\t HV: " << hv << "+-" << eHv << endl; 

   			curr.push_back(avgCurr);
   			//eCurr.push_back(eAvgCurr);

   			c->cd();
   			Imon->Draw("HISTO");
   			peakFun->Draw("SAME");

		}

		/*TSpectrum *currSpect = new TSpectrum(1);
		Int_t nfound = currSpect->Search(Imon,1,"",0.5);*/
		
		//TSpectrum *currSpect = new TSpectrum(2);
		//Int_t nfound = currSpect->Search(Imon,0.9,"",0.01); //0.5 - 0.15 This works, maybe also 0.96 - 0.01. Now for eco2 using 1-0.01
		//Int_t nfound = currSpect->Search(Imon,0.9,"",0.01); //Trials

		//printf("Found %d candidate peaks to fit\n",nfound);

		//First peak has to be higher than second and has also to be higher in counts

		//Double_t *xpeaks, *ypeaks;
   		//xpeaks = currSpect->GetPositionX();
   		//ypeaks = currSpect->GetPositionY();

   		/*for (int i = 0; i < nfound; i++) {
   			cout << "X position peak " << i+1 << ":" << xpeaks[i] << "\t" << ypeaks[i] << endl;
   		}*/

   		/*peakFun->SetParLimits(0,ypeaks[0]-5,ypeaks[0]+5);
   		peakFun->SetParLimits(1,xpeaks[0]-0.01,xpeaks[0]+0.01);
   		
   		peakFun->SetParLimits(3,ypeaks[1]-5,ypeaks[1]+5);
   		peakFun->SetParLimits(4,xpeaks[1]-0.01,xpeaks[1]+0.01);

   		Imon->Fit("peakFun","M+");
   		
   		sxCurr = peakFun->GetParameter(1) - 5*peakFun->GetParameter(2);
   		dxCurr = peakFun->GetParameter(1) + 5*peakFun->GetParameter(2);
   		cout << "Current range: " << sxCurr << "----" << dxCurr << endl;

   		binSx = Imon->FindBin(sxCurr);
   		binDx = Imon->FindBin(dxCurr);
   		cout << "Bin range: " << binSx << "----" << binDx << endl;

   		Imon->GetXaxis()->SetRange(binSx,binDx);

   		Imon->Fit("peakFun","RM+");

   		sxCurr = peakFun->GetParameter(1) - 5*peakFun->GetParameter(2);
   		dxCurr = peakFun->GetParameter(1) + 5*peakFun->GetParameter(2);
   		cout << "Current range: " << sxCurr << "----" << dxCurr << endl;

   		binSx = Imon->FindBin(sxCurr);
   		binDx = Imon->FindBin(dxCurr);
   		cout << "Bin range: " << binSx << "----" << binDx << endl;

   		Imon->GetXaxis()->SetRange(binSx,binDx);

   		avgCurr = Imon->GetMean();
   		eAvgCurr = Imon->GetMeanError();

   		cout << "Avg curr: " << avgCurr << "+-" << eAvgCurr << endl; 

   		curr.push_back(avgCurr);
   		eCurr.push_back(eAvgCurr);

   		TCanvas *c = new TCanvas();
   		c->cd();
   		Imon->Draw("HISTO");
   		peakFun->Draw("SAME");*/
   		
		fout->cd(); //open .root file
		cdtof[point]->cd();  //Enter the folder corresponding to the current HV point
		c->Write(("Current_ALICE_"+cartella).c_str());
		//Imon->Write(("Current_ALICE_"+cartella).c_str()); //Save time profile
	}

	TCanvas *cCurr = new TCanvas();
	cCurr->cd();
	TGraphErrors *currHv = new TGraphErrors(voltage.size(),&voltage[0],&curr[0],&eVoltage[0],&eCurr[0]);
	currHv->SetMarkerStyle(8);
	currHv->SetMarkerSize(0.8);
	currHv->SetMarkerColor(kRed);
	//currHv->GetXaxis()->SetRangeUser(0,,voltage.back()+0.05*voltage.back());
	currHv->GetYaxis()->SetRangeUser(0,curr.back()+0.05*curr.back());
	currHv->GetXaxis()->SetTitle("HV_{eff} [V]");
	currHv->GetYaxis()->SetTitle("I_{mon} [#muA]");
	currHv->Draw("AP");

	vector<double> currDens, eCurrDens;

	for (int i = 0; i < currHv->GetN(); i++) {
		cout << i << endl;
		currDens.push_back(curr.at(i)/area);
		eCurrDens.push_back(eCurr.at(i)/area);
	}

	TCanvas *cCurrDens = new TCanvas();
	cCurrDens->cd();
	TGraphErrors *currDensHv = new TGraphErrors(voltage.size(),&voltage[0],&currDens[0],&eVoltage[0],&eCurrDens[0]);
	//currDensHv->SetName("Graph");
	currDensHv->SetMarkerStyle(8);
	currDensHv->SetMarkerSize(0.8);
	currDensHv->SetMarkerColor(kRed);
	//currDensHv->GetXaxis()->SetRangeUser(0,,voltage.back()+0.05*voltage.back());
	//currDensHv->GetYaxis()->SetRangeUser(0,curr.back()+0.05*curr.back());
	currDensHv->GetXaxis()->SetTitle("HV_{eff} [V]");
	currDensHv->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	currDensHv->Draw("AP");

	fout->cd();
	cCurr->Write(("I(HV)_curve " + to_string(scan)).c_str());
	cCurrDens->Write(("Current_density(V)_curve " + to_string(scan)).c_str());


	fout->Close();

}