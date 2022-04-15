#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include <TMultiGraph.h>
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TAttText.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>

using namespace std;
using namespace TMath;

void plotter() {

	TFile *fout = new TFile("ALICE_results.root","RECREATE");
	fout->cd();
	TDirectory *cdtof[3];

	TPaveText *ttitle = new TPaveText(0.0406284,0.902218,0.428494,0.927419,"NDC");
	ttitle->SetFillColorAlpha(0, 1);
	ttitle->AddText("ECOgas@GIF++ (ALICE, ATLAS, CMS, EP-DT, SHiP/LHCb)");

	TLegend *eff_comp = new TLegend();
	TLine *l = new TLine(0,0.6,0,2.2);
	l->SetLineColor(kBlack);

	//vector<double> hvMinusFifty;
	double hvFifty, wpOff;
	vector<double> iwp_std, iwp_eco2, iwp_eco3, idens_std, idens_eco2, idens_eco3, ratewp_std, ratewp_eco2, ratewp_eco3, effwp_std, effwp_eco2, effwp_eco3;
	vector<double> effwp_std_200, effwp_eco2_200, effwp_eco3_200;
	vector<double> e_iwp_std, e_iwp_eco2, e_iwp_eco3, e_idens_std, e_idens_eco2, e_idens_eco3, vrate, vmuon_clus_size; 
	vector<double> muon_clus_size_wp_std, muon_clus_size_wp_eco2, muon_clus_size_wp_eco3, e_muon_clus_size_wp_std, e_muon_clus_size_wp_eco2, e_muon_clus_size_wp_eco3;
	vector<double> muon_clus_size_wp_std_200, muon_clus_size_wp_eco2_200, muon_clus_size_wp_eco3_200;
	vector<double> e_ratewp_std, e_ratewp_eco2, e_ratewp_eco3, e_effwp_std, e_effwp_eco2, e_effwp_eco3;
	vector<double> e_effwp_std_200, e_effwp_eco2_200, e_effwp_eco3_200;
	vector<double> dose {0,510,2070};
	vector<double> vWp, vEffWp, vFilterTot;
	double iwp, idenswp, ratewp, effwp, effwp_200, knee, lambda, eff_max, hv50, wp, muon_clus_sizewp, muon_clus_sizewp_200, e_muon_clus_sizewp;
	double e_iwp, e_idenswp, e_ratewp, e_effwp;

	double area = 2500; //cm2
	//double area = 0.25; //m2 
	double strip_area = 150; //Area of the ALICE RPC, in cm2

	TMultiGraph *mcurrent_std = new TMultiGraph();
	TMultiGraph *mcurrent_eco2 = new TMultiGraph();
	TMultiGraph *mcurrent_eco3 = new TMultiGraph();

	TMultiGraph *mcurrentdens_std = new TMultiGraph();
	TMultiGraph *mcurrentdens_eco2 = new TMultiGraph();
	TMultiGraph *mcurrentdens_eco3 = new TMultiGraph();

	TMultiGraph *mrate_std = new TMultiGraph();
	TMultiGraph *mrate_eco2 = new TMultiGraph();
	TMultiGraph *mrate_eco3 = new TMultiGraph();

	TMultiGraph *meff_std = new TMultiGraph();
	TMultiGraph *meff_curr_std = new TMultiGraph();
	TMultiGraph *meff_eco2 = new TMultiGraph();
	TMultiGraph *meff_curr_eco2 = new TMultiGraph();
	TMultiGraph *meff_eco3 = new TMultiGraph();
	TMultiGraph *meff_curr_eco3 = new TMultiGraph();

	TMultiGraph *mmuon_cs_std = new TMultiGraph();
	TMultiGraph *mmuon_cs_eco2 = new TMultiGraph();
	TMultiGraph *mmuon_cs_eco3 = new TMultiGraph();

	TMultiGraph *meff_wp = new TMultiGraph();
	TMultiGraph *meff_wp_abs = new TMultiGraph();
	TMultiGraph *mi_wp = new TMultiGraph();
	TMultiGraph *mi_wp_abs = new TMultiGraph();
	TMultiGraph *mi_dens_wp = new TMultiGraph();
	TMultiGraph *mi_dens_wp_abs = new TMultiGraph();
	TMultiGraph *mi_rate_wp = new TMultiGraph();
	TMultiGraph *mmuon_cs_wp = new TMultiGraph();
	TMultiGraph *mmuon_cs_wp_abs = new TMultiGraph();

	TMultiGraph *mEffMixComp = new TMultiGraph();

	TLegend *labs_std = new TLegend();
	TLegend *labs_eco2 = new TLegend();
	TLegend *labs_eco3 = new TLegend();

	TLegend *lrate_std = new TLegend();
	TLegend *lrate_eco2 = new TLegend();
	TLegend *lrate_eco3 = new TLegend();
	TLegend *lMix = new TLegend();

	int bin_min = 0, bin_max = 0;
	double min, max;
	double err_min, err_max, errpoint;

	//string plane = "X";
	string plane1 = "Y";
	string plane2 = "X";

	string fol = " ";
	string mix = " ";

	vector<double> vfilterup_std, vfilterdown_std, vfilterup_eco2, vfilterdown_eco2, vfilterup_eco3, vfilterdown_eco3;

	TF1 *fiteff1[3];
	TPaveText *tEff[3];
 
	for (int j = 0; j < 3; j++) {

		if (j == 0) {
			mix = "std";
			fol = "/home/luca/cernbox/PhD/TB_october_2021/runs/std";
			fout->cd();
			string cartella = "STD";
			cdtof[0] = fout->mkdir(cartella.c_str());  
		}

		else if (j == 1) {
			mix = "eco2";
			fol = "/home/luca/cernbox/PhD/TB_october_2021/runs/eco2";
			fout->cd();
			string cartella = "ECO2";
			cdtof[1] = fout->mkdir(cartella.c_str());
		}

		else if (j == 2) {
			mix = "eco3";
			fol = "/home/luca/cernbox/PhD/TB_october_2021/runs/eco3";
			fout->cd();
			string cartella = "ECO3";
			cdtof[2] = fout->mkdir(cartella.c_str());
		}

		cout << fol << endl;

		gSystem->cd(fol.c_str());
		remove("ALICE_data.txt");
		ofstream aliceData;
		aliceData.open("ALICE_data.txt"); //output txt file with results from all three mixtures
		aliceData << mix << endl;
		aliceData << "HV [V]" << "\t" << "I [uA]" << "\t" << "I_dens [uA/cm2]" << "\t" << "Eff [%]" << "\t" << "c.s. [strips]" << "\t" << "rate [Hz/cm2]" << endl;

		ifstream hrun;
		hrun.open("runs_filters.txt");

		int run, points, accept; 
		vector<int> vrun, vpoints; 

		double filterup, filterdown;
		vector <double> vfilterup, vfilterdown;

		while(hrun >> accept >> run >> filterup >> filterdown >> points) {
			if (accept == 1) {
				vrun.push_back(run);
				vfilterup.push_back(filterup);
				vfilterdown.push_back(filterdown);
				vpoints.push_back(points);
			}

			if (j == 0 && accept == 1) {
				vfilterup_std.push_back(filterup);
				vfilterdown_std.push_back(filterdown);
			}

			else if (j == 1 && accept == 1) {
				vfilterup_eco2.push_back(filterup);
				vfilterdown_eco2.push_back(filterdown);
			}
			
			else if (j == 2 && accept == 1) {
				vfilterup_eco3.push_back(filterup);
				vfilterdown_eco3.push_back(filterdown);	
			}
		}

		for (unsigned int i = 0; i < vrun.size(); i++) {

			aliceData << endl << vrun.at(i) << endl;

			//Open the various .root files with the plots relative to ALICE
			//2D efficiency (HV)
			cout << "Opening file: outfile_ALICE_2D_"+to_string(vrun.at(i))+".root" << endl;
			TFile *f = new TFile(("outfile_ALICE_2D_"+to_string(vrun.at(i))+".root").c_str(),"READ");

			f->cd();
			TCanvas *ceff = (TCanvas*)f->Get(("Eff_curve_2D_"+to_string(vrun.at(i))).c_str());
			TGraphErrors *effv = (TGraphErrors*)ceff->GetListOfPrimitives()->FindObject("Graph");
			TF1 *fiteff = (TF1*)ceff->GetListOfPrimitives()->FindObject("fit_eff");

			if (vfilterup.at(i) == 0 && (vrun.at(i) <= 4818 && vrun.at(i) >= 4783)) labs_std->AddEntry(effv,"Source OFF - STD","p");
			else if (vfilterup.at(i) == 0 && (vrun.at(i) <= 4992 && vrun.at(i) >= 4946)) labs_eco2->AddEntry(effv,"Source OFF - ECO2","p");
			else if (vfilterup.at(i) == 0 && (vrun.at(i) <= 4897 && vrun.at(i) >= 4868)) labs_eco3->AddEntry(effv,"Source OFF - ECO3","p");

			else if ((vfilterup.at(i) != 0 && (vrun.at(i) <= 4818 && vrun.at(i) >= 4783))) labs_std->AddEntry(effv,Form("ABS: %3.1lf - STD mix",vfilterup.at(i)),"p");
			else if ((vfilterup.at(i) != 0 && (vrun.at(i) <= 4992 && vrun.at(i) >= 4946))) labs_eco2->AddEntry(effv,("ABS " + to_string((int)vfilterup.at(i)) + " - ECO2 mix").c_str(),"p");
			else if ((vfilterup.at(i) != 0 && (vrun.at(i) <= 4897 && vrun.at(i) >= 4868))) labs_eco3->AddEntry(effv,("ABS " + to_string((int)vfilterup.at(i)) + " - ECO3 mix").c_str(),"p");

			effv->GetXaxis()->SetTitle("HV_{eff} [V]");
			effv->GetYaxis()->SetTitle("#varepsilon [%]");
			effv->SetMarkerColor(i+1);
			effv->SetMarkerStyle(i+35);
			effv->SetMarkerSize(2.5);
			fiteff->SetLineColor(i+1);
			effv->Fit(fiteff);

			TGraphErrors *effv1 = (TGraphErrors*)effv->Clone();

			//compute WP & eff @ WP
			Double_t *x = effv->GetX();
			int numPoints = effv->GetN();
			Double_t *yEff = effv->GetY();

			cout << fiteff->GetParameter(0) << endl;
			cout << fiteff->GetParameter(1) << endl;
			cout << fiteff->GetParameter(2) << endl;

			eff_max = fiteff->GetParameter(0);
			lambda = fiteff->GetParameter(1);
			hv50 = fiteff->GetParameter(2);
			//wp = (TMath::Log(19)/lambda)+hv50+150; //CMS definition
			wp = fiteff->GetX(0.95*eff_max,x[0],x[numPoints-1]); //Ecogas paper definition, 95% of the maximum efficiency

			//effwp = fiteff->Eval(wp);	

			vWp.push_back(wp);
			vEffWp.push_back(effwp);
			vFilterTot.push_back(vfilterup.at(i));
			
			cout << "WP " << wp << endl;

			double hvMinusFifty[numPoints-1];

			hvFifty = fiteff->GetX(50,x[0],x[numPoints-1]); //HV where efficiency is 50% (not of the maximum, 50% absolute value)
			/*for (int zz = 0; zz < numPoints-1; zz++) {
				//hvMinusFifty.push_back(x[zz] - wp);
				hvMinusFifty[zz] = x[zz] - wp;
			}*/
			//cout << "WP " << x[sizeof(x)/sizeof(x[0])] << endl;

			if (vrun.at(i) == 4811) {
				//effv1->Fit(fiteff1);
				effv1->SetMarkerColor(kBlack);
				//fiteff->SetLineColor(kBlack);
				fiteff1[0] = fiteff;
				mEffMixComp->Add(effv1);
				lMix->AddEntry(effv1,"STD mix","p");
				wpOff = wp;
			}

			if (vrun.at(i) == 4960) {
				//effv1->Fit(fiteff1);
				effv1->SetMarkerColor(kRed);
				//fiteff->SetLineColor(kRed);
				fiteff1[1] = fiteff;
				mEffMixComp->Add(effv1);
				lMix->AddEntry(effv1,"ECO2 mix","p");
				wpOff = wp;
			}

			if (vrun.at(i) == 4897) {
				//effv1->Fit(fiteff1);
				effv1->SetMarkerColor(kGreen);
				//fiteff->SetLineColor(kGreen);
				fiteff1[2] = fiteff;
				mEffMixComp->Add(effv1);
				lMix->AddEntry(effv1,"ECO3 mix","p");
				wpOff = wp;
			}

			for (int zz = 0; zz < numPoints-1; zz++) {
				//hvMinusFifty.push_back(x[zz] - wp);
				//hvMinusFifty[zz] = x[zz] - wp;
				hvMinusFifty[zz] = x[zz] - wpOff;
			}

			cout << endl << endl << "WP OFF: " << wpOff << endl << endl;

			effwp = fiteff->Eval(wpOff);
			effwp_200 = fiteff->Eval(wpOff+200);
			
			min = x[0];
			max = x[0];
			for (int k = 0; k < effv->GetN(); k++) {
				//cout << x[k] << endl;
				if (x[k] <= wp) {
					min = x[k];
					bin_min = k;
					continue;
				}
				else if (x[k] >= wp) {
					max = x[k];
					bin_max = k;
					cout << "min: " << min << " avg " << wp << ", max: " << max << endl;
					cout << "bin min: " << bin_min << ", bin max: " << bin_max << endl;
					break;
				}
			}	

			err_min = effv->GetErrorY(bin_min);
			err_max = effv->GetErrorY(bin_max);
			e_effwp = (err_min + err_max)/2;

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) {
				effwp_std.push_back(effwp);
				effwp_std_200.push_back(effwp_200);
				e_effwp_std.push_back(e_effwp);
				meff_std->Add(effv);
				meff_curr_std->Add(effv);
			}
			
			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) {
				effwp_eco2.push_back(effwp);
				effwp_eco2_200.push_back(effwp_200);
				e_effwp_eco2.push_back(e_effwp);
				meff_eco2->Add(effv);
				meff_curr_eco2->Add(effv);
			}
			
			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) {
				effwp_eco3.push_back(effwp);
				effwp_eco3_200.push_back(effwp_200);
				e_effwp_eco3.push_back(e_effwp);
				meff_eco3->Add(effv);
				meff_curr_eco3->Add(effv);
			}

			//I(HV) .root file
			cout << "Opening file: current_ALICE"+to_string(vrun.at(i))+".root" << endl;
			cout << "I(HV)_curve " + to_string(vrun.at(i)) << endl;
			TFile *f1 = new TFile(("current_ALICE"+to_string(vrun.at(i))+".root").c_str(),"READ");

			//Current
			TCanvas *civ = (TCanvas*)f1->Get(("I(HV)_curve " + to_string(vrun.at(i))).c_str()); //Get the canvas with the I(V) curve
			TGraphErrors *iv = (TGraphErrors*)civ->GetListOfPrimitives()->FindObject("Graph"); 
			Double_t *yCurr = iv->GetY();
			iv->SetMarkerColor(i+1);
			iv->SetMarkerStyle(i+35);
			iv->SetMarkerSize(2);
			
			iwp = iv->Eval(wpOff);
			cout << endl << iwp << endl;
			err_min = iv->GetErrorY(bin_min);
			err_max = iv->GetErrorY(bin_max);
			e_iwp = (err_min + err_max)/2;

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) {
				iwp_std.push_back(iwp);
				e_iwp_std.push_back(e_iwp);
				mcurrent_std->Add(iv);
				meff_curr_std->Add(iv);
			}
			
			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) {
				iwp_eco2.push_back(iwp);
				e_iwp_eco2.push_back(e_iwp);
				mcurrent_eco2->Add(iv);
				meff_curr_eco2->Add(iv);
			}
			
			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) {
				iwp_eco3.push_back(iwp);
				e_iwp_eco3.push_back(e_iwp);
				mcurrent_eco3->Add(iv);
				meff_curr_eco3->Add(iv);
			}

			//Current density
			TCanvas *civdens = (TCanvas*)f1->Get(("Current_density(V)_curve "+to_string(vrun.at(i))).c_str());
			TGraphErrors *ivdens = (TGraphErrors*)civdens->GetListOfPrimitives()->FindObject("Graph"); 
			Double_t *yCrrDens = ivdens->GetY();
			ivdens->SetMarkerColor(i+1);
			ivdens->SetMarkerStyle(i+35);
			ivdens->SetMarkerSize(2);
			
			idenswp = ivdens->Eval(wp);
			err_min = ivdens->GetErrorY(bin_min);
			err_max = ivdens->GetErrorY(bin_max);
			e_idenswp = (err_min + err_max)/2;

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) {
				idens_std.push_back(idenswp);
				e_idens_std.push_back(e_idenswp);
				mcurrentdens_std->Add(ivdens);
			}

			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) {
				idens_eco2.push_back(idenswp);
				e_idens_eco2.push_back(e_idenswp);
				mcurrentdens_eco2->Add(ivdens);
			}

			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) {
				idens_eco3.push_back(idenswp);
				e_idens_eco3.push_back(e_idenswp);
				mcurrentdens_eco3->Add(ivdens);
			}

			cout << "Opening file: outfile_ALICE"+plane1+to_string(vrun.at(i))+".root" << endl;
			TFile *f2 = new TFile(("outfile_ALICE"+plane1+to_string(vrun.at(i))+".root").c_str(),"READ");

			//Rate as a function of HV
			//TCanvas *crate = (TCanvas*)f2->Get(("Rate(V)_curve_" + plane1 + " " +to_string(vrun.at(i))).c_str()); //Get the canvas with the rate(V) curve for different abs 
			TCanvas *crate = (TCanvas*)f2->Get(("Cluster_rate(V)_curve_" + plane1 + " " +to_string(vrun.at(i))).c_str()); //Get the canvas with the rate(V) curve for different abs 
			TGraphErrors *ratev = (TGraphErrors*)crate->GetListOfPrimitives()->FindObject("Graph");
			Double_t *yRate = ratev->GetY();
			ratev->SetMarkerColor(i+1);
			ratev->SetMarkerStyle(i+35);

			ratewp = ratev->Eval(wp);
			err_min = ratev->GetErrorY(bin_min);
			err_max = ratev->GetErrorY(bin_max);
			e_ratewp = (err_min + err_max)/2;

			vrate.push_back(ratewp);

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) {
				ratewp_std.push_back(ratewp);
				e_ratewp_std.push_back(e_ratewp);
				mrate_std->Add(ratev);
			}

			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) {
				ratewp_eco2.push_back(ratewp);
				e_ratewp_eco2.push_back(e_ratewp);
				mrate_eco2->Add(ratev);
			}

			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) {
				ratewp_eco3.push_back(ratewp);
				e_ratewp_eco3.push_back(e_ratewp);	
				mrate_eco3->Add(ratev);		
			}

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) lrate_std->AddEntry(effv,Form("%3.0lf Hz/cm^{2} STD mix",ratewp),"p");
			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) lrate_eco2->AddEntry(effv,Form("%3.0lf Hz/cm^{2} ECO2 mix",ratewp),"p");
			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) lrate_eco3->AddEntry(effv,Form("%3.0lf Hz/cm^{2} ECO3 mix",ratewp),"p");


			//Muon cluster size as a function of HV
			cout << "Opening file: outfile_ALICE"+plane2+to_string(vrun.at(i))+".root" << endl;
			TFile *f3 = new TFile(("outfile_ALICE"+plane2+to_string(vrun.at(i))+".root").c_str(),"READ");

			TCanvas *cclus_size_muon = (TCanvas*)f3->Get(("Muon_cluster_size(V)_curve_" + plane2 + " " + to_string(vrun.at(i))).c_str()); //Get the canvas with the cluster size(V) curve for different abs 
			TGraphErrors *muon_clus_sizev = (TGraphErrors*)cclus_size_muon->GetListOfPrimitives()->FindObject("Graph");
			Double_t *yClus = muon_clus_sizev->GetY();
			Double_t *yClusError = muon_clus_sizev->GetEY();
			muon_clus_sizev->SetMarkerColor(i+1);
			muon_clus_sizev->SetMarkerStyle(i+35);

			//muon_clus_sizewp = muon_clus_sizev->Eval(wp);

			muon_clus_sizewp = muon_clus_sizev->Eval(wpOff);
			muon_clus_sizewp_200 = muon_clus_sizev->Eval(wpOff+200);

			err_min = muon_clus_sizev->GetErrorY(bin_min);
			err_max = muon_clus_sizev->GetErrorY(bin_max);
			e_muon_clus_sizewp = (err_min + err_max)/2;

			vmuon_clus_size.push_back(muon_clus_sizewp);

			TGraphErrors *muon_clus_sizev_minus50 = new TGraphErrors(numPoints-1,hvMinusFifty,yClus,NULL,yClusError);
			muon_clus_sizev_minus50->GetXaxis()->SetTitle("HV-wp");
			muon_clus_sizev_minus50->GetYaxis()->SetTitle("Muon C.S. [strips]");
			muon_clus_sizev_minus50->SetMarkerColor(i+1);
			muon_clus_sizev_minus50->SetMarkerStyle(i+35);
			muon_clus_sizev_minus50->SetMarkerSize(2);

			if (vrun.at(i) <= 4818 && vrun.at(i) >= 4783) {
				muon_clus_size_wp_std.push_back(muon_clus_sizewp);
				muon_clus_size_wp_std_200.push_back(muon_clus_sizewp_200);
				e_muon_clus_size_wp_std.push_back(e_muon_clus_sizewp);
				//mmuon_cs_std->Add(muon_clus_sizev);
				mmuon_cs_std->Add(muon_clus_sizev_minus50);
			}

			else if (vrun.at(i) <= 4992 && vrun.at(i) >= 4946) {
				muon_clus_size_wp_eco2.push_back(muon_clus_sizewp);
				muon_clus_size_wp_eco2_200.push_back(muon_clus_sizewp_200);
				e_muon_clus_size_wp_eco2.push_back(e_muon_clus_sizewp);
				//mmuon_cs_eco2->Add(muon_clus_sizev);
				mmuon_cs_eco2->Add(muon_clus_sizev_minus50);
			}

			else if (vrun.at(i) <= 4897 && vrun.at(i) >= 4868) {
				muon_clus_size_wp_eco3.push_back(muon_clus_sizewp);
				muon_clus_size_wp_eco3_200.push_back(muon_clus_sizewp_200);
				e_muon_clus_size_wp_eco3.push_back(e_muon_clus_sizewp);	
				//mmuon_cs_eco3->Add(muon_clus_sizev);	
				mmuon_cs_eco3->Add(muon_clus_sizev_minus50);	
			}

			for (int kk = 0; kk < numPoints; kk++) {
				aliceData << x[kk] << "\t" << yCurr[kk] << "\t" << yCrrDens[kk] << "\t" << yEff[kk] << "\t" << yClus[kk] << "\t" << yRate[kk] << endl;
			}
		}
		aliceData.close();

		tEff[j] = new TPaveText(0.65,0.2,0.85,0.4,"brNDC");
		for (int i = 0; i < vrun.size(); i++) {
			tEff[j]->AddText(Form("Abs: %3.1lf, WP: %3.1lf, eff WP: %3.1lf",vFilterTot.at(i),vWp.at(i),vEffWp.at(i)));
		}
		vWp.clear();
		vEffWp.clear();
		vFilterTot.clear();
	}

	///////////////////////
	//  				 //
	//  Legend with ABS  //  
	//					 //
	///////////////////////
	
	//Efficiency (V) curve
	TCanvas *cefficiency_std = new TCanvas();
	cefficiency_std->cd();
	meff_std->SetTitle("STD gas mixture");
	meff_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_std->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_std->GetYaxis()->SetRangeUser(0,100);
	meff_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[0]->Draw("SAME");

	TCanvas *cefficiency_eco2 = new TCanvas();
	cefficiency_eco2->cd();
	meff_eco2->SetTitle("ECO2 gas mixture");
	meff_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_eco2->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_eco2->GetYaxis()->SetRangeUser(0,100);
	meff_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[1]->Draw("SAME");

	TCanvas *cefficiency_eco3 = new TCanvas();
	cefficiency_eco3->cd();
	meff_eco3->SetTitle("ECO3 gas mixture");
	meff_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_eco3->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_eco3->GetYaxis()->SetRangeUser(0,100);
	meff_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[2]->Draw("SAME");

	//Efficiency and current (V) curve
	/*TCanvas *cefficiency_curr_std = new TCanvas();
	cefficiency_curr_std->cd();
	meff_curr_std->SetTitle("STD gas mixture");
	meff_curr_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_curr_std->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_curr_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[0]->Draw("SAME");

	TCanvas *cefficiency_curr_eco2 = new TCanvas();
	cefficiency_curr_eco2->cd();
	meff_curr_eco2->SetTitle("ECO2 gas mixture");
	meff_curr_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_curr_eco2->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_curr_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[1]->Draw("SAME");

	TCanvas *cefficiency_curr_eco3 = new TCanvas();
	cefficiency_curr_eco3->cd();
	meff_curr_eco3->SetTitle("ECO3 gas mixture");
	meff_curr_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	meff_curr_eco3->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_curr_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");
	tEff[2]->Draw("SAME");*/

	//I(V) curve
	TCanvas *ccurr_std = new TCanvas();
	ccurr_std->cd();
	mcurrent_std->SetTitle("STD gas mixture");
	mcurrent_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrent_std->GetYaxis()->SetTitle("I [#muA]");
	mcurrent_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurr_eco2 = new TCanvas();
	ccurr_eco2->cd();
	mcurrent_eco2->SetTitle("ECO2 gas mixture");
	mcurrent_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrent_eco2->GetYaxis()->SetTitle("I [#muA]");
	mcurrent_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurr_eco3 = new TCanvas();
	ccurr_eco3->cd();
	mcurrent_eco3->SetTitle("ECO3 gas mixture");
	mcurrent_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrent_eco3->GetYaxis()->SetTitle("I [#muA]");
	mcurrent_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Current density (V) curve
	TCanvas *ccurrdens_std = new TCanvas();
	ccurrdens_std->cd();
	mcurrentdens_std->SetTitle("STD gas mixture");
	mcurrentdens_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrentdens_std->GetYaxis()->SetTitle("Current density [#muA/m^{2}]");
	mcurrentdens_std->GetYaxis()->SetRangeUser(0,0.07);
	mcurrentdens_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *c = new TCanvas();

	TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
 	TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
 	p2->SetFillStyle(4000); // will be transparent
  
  	p1->Draw();
  	p1->cd();
  	meff_std->Draw("AP");
  	gPad->Update();
  
  	Double_t xmin = p1->GetUxmin();
 	Double_t xmax = p1->GetUxmax();
  	Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
  	Double_t ymin = mcurrentdens_std->GetHistogram()->GetMinimum();
  	Double_t ymax = mcurrentdens_std->GetHistogram()->GetMaximum();
  	Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
  	p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
  	p2->Draw();
  	p2->cd();
  	mcurrentdens_std->Draw("P");
  	gPad->Update();
  
  	TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
  	axis->SetTitleFont(42);
  	axis->SetLabelFont(42);
  	axis->SetTitle("Current Density [#muA/cm^{2}]");
  	axis->Draw();
  	gPad->Update();
	
	TCanvas *ccurrdens_eco2 = new TCanvas();
	ccurrdens_eco2->cd();
	mcurrentdens_eco2->SetTitle("ECO2 gas mixture");
	mcurrentdens_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrentdens_eco2->GetYaxis()->SetTitle("Current density [#muA/m^{2}]");
	mcurrentdens_eco2->GetYaxis()->SetRangeUser(0,0.07);
	mcurrentdens_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *c1 = new TCanvas();

	TPad *p3 = new TPad("p3", "", 0, 0, 1, 1);
 	TPad *p4 = new TPad("p4", "", 0, 0, 1, 1);
 	p4->SetFillStyle(4000); // will be transparent
  
  	p3->Draw();
  	p3->cd();
  	meff_eco2->Draw("AP");
  	gPad->Update();
  
  	Double_t xmin1 = p3->GetUxmin();
 	Double_t xmax1 = p3->GetUxmax();
  	Double_t dx1 = (xmax1 - xmin1) / 0.8; // 10 percent margins left and right
  	Double_t ymin1 = mcurrentdens_eco2->GetHistogram()->GetMinimum();
  	Double_t ymax1 = mcurrentdens_eco2->GetHistogram()->GetMaximum();
  	Double_t dy1 = (ymax1 - ymin1) / 0.8; // 10 percent margins top and bottom
  	p4->Range(xmin1-0.1*dx1, ymin1-0.1*dy1, xmax1+0.1*dx1, ymax1+0.1*dy1);
  	p4->Draw();
  	p4->cd();
  	mcurrentdens_eco2->Draw("P");
  	gPad->Update();
  
  	TGaxis *axis1 = new TGaxis(xmax1, ymin1, xmax1, ymax1, ymin1, ymax1, 510, "+L");
  	axis1->SetTitleFont(42);
  	axis1->SetLabelFont(42);
  	axis1->SetTitle("Current Density [#muA/cm^{2}]");
  	axis1->Draw();
  	gPad->Update();

	TCanvas *ccurrdens_eco3 = new TCanvas();
	ccurrdens_eco3->cd();
	mcurrentdens_eco3->SetTitle("ECO3 gas mixture");
	mcurrentdens_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	mcurrentdens_eco3->GetYaxis()->SetTitle("Current density [#muA/m^{2}]");
	mcurrentdens_eco3->GetYaxis()->SetRangeUser(0,0.07);
	mcurrentdens_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *c2 = new TCanvas();

	TPad *p5 = new TPad("p5", "", 0, 0, 1, 1);
 	TPad *p6 = new TPad("p6", "", 0, 0, 1, 1);
 	p6->SetFillStyle(4000); // will be transparent
  
  	p5->Draw();
  	p5->cd();
  	meff_eco3->Draw("AP");
  	gPad->Update();
  
  	Double_t xmin2 = p5->GetUxmin();
 	Double_t xmax2 = p5->GetUxmax();
  	Double_t dx2 = (xmax2 - xmin2) / 0.8; // 10 percent margins left and right
  	Double_t ymin2 = mcurrentdens_eco3->GetHistogram()->GetMinimum();
  	Double_t ymax2 = mcurrentdens_eco3->GetHistogram()->GetMaximum();
  	Double_t dy2 = (ymax2 - ymin2) / 0.8; // 10 percent margins top and bottom
  	p6->Range(xmin2-0.1*dx2, ymin2-0.1*dy2, xmax2+0.1*dx2, ymax2+0.1*dy2);
  	p6->Draw();
  	p6->cd();
  	mcurrentdens_eco3->Draw("P");
  	gPad->Update();
  
  	TGaxis *axis2 = new TGaxis(xmax2, ymin2, xmax2, ymax2, ymin2, ymax2, 510, "+L");
  	axis2->SetTitleFont(42);
  	axis2->SetLabelFont(42);
  	axis2->SetTitle("Current Density [#muA/cm^{2}]");
  	axis2->Draw();
  	gPad->Update();

	//Gamma rate (V) curve
	TCanvas *cgamma_std = new TCanvas();
	cgamma_std->cd();
	mrate_std->SetTitle("STD gas mixture");
	mrate_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	mrate_std->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	mrate_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cgamma_eco2 = new TCanvas();
	cgamma_eco2->cd();
	mrate_eco2->SetTitle("ECO2 gas mixture");
	mrate_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	mrate_eco2->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	mrate_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cgamma_eco3 = new TCanvas();
	cgamma_eco3->cd();
	mrate_eco3->SetTitle("ECO3 gas mixture");
	mrate_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	mrate_eco3->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	mrate_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Muon cs (V) curve
	TCanvas *cmuon_cs_std = new TCanvas();
	cmuon_cs_std->cd();
	mmuon_cs_std->SetTitle("STD gas mixture");
	//mmuon_cs_std->GetXaxis()->SetTitle("HV_{eff} [V]");
	mmuon_cs_std->GetXaxis()->SetTitle("(HV_{eff} - knee) [V]");
	mmuon_cs_std->GetYaxis()->SetTitle("Muon cluster size [strips]");
	mmuon_cs_std->GetYaxis()->SetRangeUser(0.6,2.2);
	mmuon_cs_std->Draw("AP");
	labs_std->Draw("SAME");
	ttitle->Draw("SAME");
	l->Draw("SAME");

	TCanvas *cmuon_cs_eco2 = new TCanvas();
	cmuon_cs_eco2->cd();
	mmuon_cs_eco2->SetTitle("ECO2 gas mixture");
	//mmuon_cs_eco2->GetXaxis()->SetTitle("HV_{eff} [V]");
	mmuon_cs_eco2->GetXaxis()->SetTitle("(HV_{eff} - knee) [V]");
	mmuon_cs_eco2->GetYaxis()->SetTitle("Muon cluster size [strips]");
	mmuon_cs_eco2->GetYaxis()->SetRangeUser(0.6,2.2);
	mmuon_cs_eco2->Draw("AP");
	labs_eco2->Draw("SAME");
	ttitle->Draw("SAME");
	l->Draw("SAME");

	TCanvas *cmuon_cs_eco3 = new TCanvas();
	cmuon_cs_eco3->cd();
	mmuon_cs_eco3->SetTitle("ECO3 gas mixture");
	//mmuon_cs_eco3->GetXaxis()->SetTitle("HV_{eff} [V]");
	mmuon_cs_eco3->GetXaxis()->SetTitle("(HV_{eff} - knee) [V]");
	mmuon_cs_eco3->GetYaxis()->SetTitle("Muon cluster size [strips]");
	mmuon_cs_eco3->GetYaxis()->SetRangeUser(0.6,2.2);
	mmuon_cs_eco3->Draw("AP");
	labs_eco3->Draw("SAME");
	ttitle->Draw("SAME");
	l->Draw("SAME");

	/////////////////////////
	//   				   //
	//   Legend with rate  //  
	//					   //
	/////////////////////////

	//I(V) curve
	TCanvas *ccurr_rate_std = new TCanvas();
	ccurr_rate_std->cd();
	mcurrent_std->Draw("AP");
	lrate_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurr_rate_eco2 = new TCanvas();
	ccurr_rate_eco2->cd();
	mcurrent_eco2->Draw("AP");
	lrate_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurr_rate_eco3 = new TCanvas();
	ccurr_rate_eco3->cd();
	mcurrent_eco3->Draw("AP");
	lrate_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Current density (V) curve
	TCanvas *ccurrdens_rate_std = new TCanvas();
	ccurrdens_rate_std->cd();
	mcurrentdens_std->Draw("AP");
	lrate_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurrdens_rate_eco2 = new TCanvas();
	ccurrdens_rate_eco2->cd();
	mcurrentdens_eco2->Draw("AP");
	lrate_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *ccurrdens_rate_eco3 = new TCanvas();
	ccurrdens_rate_eco3->cd();
	mcurrentdens_eco3->Draw("AP");
	lrate_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Efficiency (V) curve
	TCanvas *cefficiency_rate_std = new TCanvas();
	cefficiency_rate_std->cd();
	meff_std->Draw("AP");
	lrate_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cefficiency_rate_eco2 = new TCanvas();
	cefficiency_rate_eco2->cd();
	meff_eco2->Draw("AP");
	lrate_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cefficiency_rate_eco3 = new TCanvas();
	cefficiency_rate_eco3->cd();
	meff_eco3->Draw("AP");
	lrate_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Muon cs (V) curve
	TCanvas *cmuon_cs_rate_std = new TCanvas();
	cmuon_cs_rate_std->cd();
	mmuon_cs_std->Draw("AP");
	lrate_std->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cmuon_cs_rate_eco2 = new TCanvas();
	cmuon_cs_rate_eco2->cd();
	mmuon_cs_eco2->Draw("AP");
	lrate_eco2->Draw("SAME");
	ttitle->Draw("SAME");

	TCanvas *cmuon_cs_rate_eco3 = new TCanvas();
	cmuon_cs_rate_eco3->cd();
	mmuon_cs_eco3->Draw("AP");
	lrate_eco3->Draw("SAME");
	ttitle->Draw("SAME");

	//Efficiency at working point as a function of rate (eco+std) wp source OFF
	
	//TGraphErrors *g_eff_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&effwp_std[0],NULL,&e_effwp_std[0]);
	TGraphErrors *g_eff_wp_std = new TGraphErrors(ratewp_std.size(),&dose[0],&effwp_std[0],NULL,&e_effwp_std[0]);
	g_eff_wp_std->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_std->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_std->SetMarkerStyle(8);
	g_eff_wp_std->SetMarkerSize(2);
	g_eff_wp_std->SetMarkerColor(kBlack);
	meff_wp->Add(g_eff_wp_std);
	eff_comp->AddEntry(g_eff_wp_std,"STD mix knee","p");

	//TGraphErrors *g_eff_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&effwp_eco2[0],NULL,&e_effwp_eco2[0]);
	TGraphErrors *g_eff_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&effwp_eco2[0],NULL,&e_effwp_eco2[0]);
	g_eff_wp_eco2->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_eco2->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco2->SetMarkerStyle(8);
	g_eff_wp_eco2->SetMarkerSize(2);
	g_eff_wp_eco2->SetMarkerColor(kRed);
	meff_wp->Add(g_eff_wp_eco2);
	eff_comp->AddEntry(g_eff_wp_eco2,"ECO2 mix knee","p");
	
	//TGraphErrors *g_eff_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&effwp_eco3[0],NULL,&e_effwp_eco3[0]);
	TGraphErrors *g_eff_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&effwp_eco3[0],NULL,&e_effwp_eco3[0]);
	g_eff_wp_eco3->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_eco3->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco3->SetMarkerStyle(8);
	g_eff_wp_eco3->SetMarkerSize(2);
	g_eff_wp_eco3->SetMarkerColor(kBlue);
	meff_wp->Add(g_eff_wp_eco3);
	eff_comp->AddEntry(g_eff_wp_eco3,"ECO3 mix knee","p");

	//Efficiency at working point as a function of rate (eco+std) wp source OFF + 200 V

	//TGraphErrors *g_eff_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&effwp_std[0],NULL,&e_effwp_std[0]);
	TGraphErrors *g_eff_wp_std_200 = new TGraphErrors(ratewp_std.size(),&dose[0],&effwp_std_200[0],NULL,&e_effwp_std[0]);
	g_eff_wp_std_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_std_200->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_std_200->SetMarkerStyle(26);
	g_eff_wp_std_200->SetMarkerSize(2);
	g_eff_wp_std_200->SetMarkerColor(kBlack);
	meff_wp->Add(g_eff_wp_std_200);
	eff_comp->AddEntry(g_eff_wp_std_200,"STD mix knee + 200 V","p");

	//TGraphErrors *g_eff_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&effwp_eco2[0],NULL,&e_effwp_eco2[0]);
	TGraphErrors *g_eff_wp_eco2_200 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&effwp_eco2_200[0],NULL,&e_effwp_eco2[0]);
	g_eff_wp_eco2_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_eco2_200->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco2_200->SetMarkerStyle(26);
	g_eff_wp_eco2_200->SetMarkerSize(2);
	g_eff_wp_eco2_200->SetMarkerColor(kRed);
	meff_wp->Add(g_eff_wp_eco2_200);
	eff_comp->AddEntry(g_eff_wp_eco2_200,"ECO2 mix knee + 200 V","p");
	
	//TGraphErrors *g_eff_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&effwp_eco3[0],NULL,&e_effwp_eco3[0]);
	TGraphErrors *g_eff_wp_eco3_200 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&effwp_eco3_200[0],NULL,&e_effwp_eco3[0]);
	g_eff_wp_eco3_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_eff_wp_eco3_200->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco3_200->SetMarkerStyle(26);
	g_eff_wp_eco3_200->SetMarkerSize(2);
	g_eff_wp_eco3_200->SetMarkerColor(kBlue);
	meff_wp->Add(g_eff_wp_eco3_200);
	eff_comp->AddEntry(g_eff_wp_eco3_200,"ECO3 mix knee + 200 V","p");

	TCanvas *c_eff_wp = new TCanvas();
	c_eff_wp->cd();
	//meff_wp->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	meff_wp->GetXaxis()->SetTitle("Dose [#muS/h]");
	meff_wp->GetYaxis()->SetTitle("#varepsilon at fixed voltage [%]");
	meff_wp->GetYaxis()->SetRangeUser(0,100);
	meff_wp->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");
	
	//Efficiency at working point as a function of abs (eco+std)
	/*TGraphErrors *g_eff_wp_std_abs = new TGraphErrors(vfilterup_std.size(),&vfilterup_std[0],&effwp_std[0],NULL,&e_effwp_std[0]);
	g_eff_wp_std_abs->GetXaxis()->SetTitle("ABS");
	g_eff_wp_std_abs->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_std_abs->SetMarkerStyle(8);	
	g_eff_wp_std_abs->SetMarkerSize(2);
	g_eff_wp_std_abs->SetMarkerColor(kBlack);
	meff_wp_abs->Add(g_eff_wp_std_abs);

	TGraphErrors *g_eff_wp_eco2_abs = new TGraphErrors(vfilterup_eco2.size(),&vfilterup_eco2[0],&effwp_eco2[0],NULL,&e_effwp_eco2[0]);
	g_eff_wp_eco2_abs->GetXaxis()->SetTitle("ABS");
	g_eff_wp_eco2_abs->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco2_abs->SetMarkerStyle(8);
	g_eff_wp_eco2_abs->SetMarkerSize(2);
	g_eff_wp_eco2_abs->SetMarkerColor(kRed);
	meff_wp_abs->Add(g_eff_wp_eco2_abs);

	TGraphErrors *g_eff_wp_eco3_abs = new TGraphErrors(vfilterup_eco3.size(),&vfilterup_eco3[0],&effwp_eco3[0],NULL,&e_effwp_eco3[0]);
	g_eff_wp_eco3_abs->GetXaxis()->SetTitle("ABS");
	g_eff_wp_eco3_abs->GetYaxis()->SetTitle("#varepsilon [%]");
	g_eff_wp_eco3_abs->SetMarkerStyle(8);
	g_eff_wp_eco3_abs->SetMarkerSize(2);
	g_eff_wp_eco3_abs->SetMarkerColor(kBlue);
	meff_wp_abs->Add(g_eff_wp_eco3_abs);

	TCanvas *c_eff_wp_abs = new TCanvas();
	c_eff_wp_abs->cd();
	meff_wp_abs->GetXaxis()->SetTitle("ABS");
	meff_wp_abs->GetYaxis()->SetTitle("#varepsilon [%]");
	meff_wp_abs->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");*/

	//Current at working point as a function of rate
	TGraphErrors *g_i_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&iwp_std[0],NULL,&e_iwp_std[0]);
	//TGraphErrors *g_i_wp_std = new TGraphErrors(ratewp_std.size(),&dose[0],&iwp_std[0],NULL,&e_iwp_std[0]);
	g_i_wp_std->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_i_wp_std->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_std->SetMarkerStyle(8);
	g_i_wp_std->SetMarkerSize(2);
	g_i_wp_std->SetMarkerColor(kBlack);
	mi_wp->Add(g_i_wp_std);

	TGraphErrors *g_i_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&iwp_eco2[0],NULL,&e_iwp_eco2[0]);
	//TGraphErrors *g_i_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&iwp_eco2[0],NULL,&e_iwp_eco2[0]);
	g_i_wp_eco2->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_i_wp_eco2->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_eco2->SetMarkerStyle(8);
	g_i_wp_eco2->SetMarkerSize(2);
	g_i_wp_eco2->SetMarkerColor(kRed);
	mi_wp->Add(g_i_wp_eco2);

	TGraphErrors *g_i_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&iwp_eco3[0],NULL,&e_iwp_eco3[0]);
	//TGraphErrors *g_i_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&iwp_eco3[0],NULL,&e_iwp_eco3[0]);
	g_i_wp_eco3->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_i_wp_eco3->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_eco3->SetMarkerStyle(8);
	g_i_wp_eco3->SetMarkerSize(2);
	g_i_wp_eco3->SetMarkerColor(kBlue);
	mi_wp->Add(g_i_wp_eco3);

	TCanvas *c_i_wp = new TCanvas();
	c_i_wp->cd();
	//mi_wp->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	mi_wp->GetXaxis()->SetTitle("Dose [#muS/h]");
	mi_wp->GetYaxis()->SetTitle("I [#muA]");
	mi_wp->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");

	//Current at working point as a function of ABS
	/*TGraphErrors *g_i_wp_std_abs = new TGraphErrors(vfilterup_std.size(),&vfilterup_std[0],&iwp_std[0],NULL,&e_iwp_std[0]);
	g_i_wp_std_abs->GetXaxis()->SetTitle("ABS");
	g_i_wp_std_abs->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_std_abs->SetMarkerStyle(8);
	g_i_wp_std_abs->SetMarkerSize(2);
	g_i_wp_std_abs->SetMarkerColor(kBlack);
	mi_wp_abs->Add(g_i_wp_std_abs);

	TGraphErrors *g_i_wp_eco2_abs = new TGraphErrors(vfilterup_eco2.size(),&vfilterup_eco2[0],&iwp_eco2[0],NULL,&e_iwp_eco2[0]);
	g_i_wp_eco2_abs->GetXaxis()->SetTitle("ABS");
	g_i_wp_eco2_abs->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_eco2_abs->SetMarkerStyle(8);
	g_i_wp_eco2_abs->SetMarkerSize(2);
	g_i_wp_eco2_abs->SetMarkerColor(kRed);
	mi_wp_abs->Add(g_i_wp_eco2_abs);

	TGraphErrors *g_i_wp_eco3_abs = new TGraphErrors(vfilterup_eco3.size(),&vfilterup_eco3[0],&iwp_eco3[0],NULL,&e_iwp_eco3[0]);
	g_i_wp_eco3_abs->GetXaxis()->SetTitle("ABS");
	g_i_wp_eco3_abs->GetYaxis()->SetTitle("I [#muA]");
	g_i_wp_eco3_abs->SetMarkerStyle(8);
	g_i_wp_eco3_abs->SetMarkerSize(2);
	g_i_wp_eco3_abs->SetMarkerColor(kBlue);
	mi_wp_abs->Add(g_i_wp_eco3_abs);

	TCanvas *c_i_wp_abs = new TCanvas();
	c_i_wp_abs->cd();
	mi_wp_abs->GetXaxis()->SetTitle("ABS");
	mi_wp_abs->GetYaxis()->SetTitle("I [#muA]");
	mi_wp_abs->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");

	//Current density at working point as a function of rate
	//TGraphErrors *g_idens_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&idens_std[0],NULL,&e_idens_std[0]);
	TGraphErrors *g_idens_wp_std = new TGraphErrors(ratewp_std.size(),&dose[0],&idens_std[0],NULL,&e_idens_std[0]);
	g_idens_wp_std->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_idens_wp_std->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_std->SetMarkerStyle(8);
	g_idens_wp_std->SetMarkerSize(2);
	g_idens_wp_std->SetMarkerColor(kBlack);
	mi_dens_wp->Add(g_idens_wp_std);

	//TGraphErrors *g_idens_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&idens_eco2[0],NULL,&e_idens_eco2[0]);
	TGraphErrors *g_idens_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&idens_eco2[0],NULL,&e_idens_eco2[0]);
	g_idens_wp_eco2->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_idens_wp_eco2->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_eco2->SetMarkerStyle(8);
	g_idens_wp_eco2->SetMarkerSize(2);
	g_idens_wp_eco2->SetMarkerColor(kRed);
	mi_dens_wp->Add(g_idens_wp_eco2);

	//TGraphErrors *g_idens_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&idens_eco3[0],NULL,&e_idens_eco3[0]);
	TGraphErrors *g_idens_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&idens_eco3[0],NULL,&e_idens_eco3[0]);
	g_idens_wp_eco3->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_idens_wp_eco3->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_eco3->SetMarkerStyle(8);
	g_idens_wp_eco3->SetMarkerSize(2);
	g_idens_wp_eco3->SetMarkerColor(kBlue);
	mi_dens_wp->Add(g_idens_wp_eco3);

	TCanvas *c_idens_wp = new TCanvas();
	c_idens_wp->cd();
	//mi_dens_wp->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	mi_dens_wp->GetXaxis()->SetTitle("Dose [#muS/h]");
	mi_dens_wp->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	mi_dens_wp->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");

	//Current density at working point as a function of ABS
	TGraphErrors *g_idens_wp_std_abs = new TGraphErrors(vfilterup_std.size(),&vfilterup_std[0],&idens_std[0],NULL,&e_idens_std[0]);
	g_idens_wp_std_abs->GetXaxis()->SetTitle("ABS");
	g_idens_wp_std_abs->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_std_abs->SetMarkerStyle(8);
	g_idens_wp_std_abs->SetMarkerSize(2);
	g_idens_wp_std_abs->SetMarkerColor(kBlack);
	mi_dens_wp_abs->Add(g_idens_wp_std_abs);

	TGraphErrors *g_idens_wp_eco2_abs = new TGraphErrors(vfilterup_eco2.size(),&vfilterup_eco2[0],&idens_eco2[0],NULL,&e_idens_eco2[0]);
	g_idens_wp_eco2_abs->GetXaxis()->SetTitle("ABS");
	g_idens_wp_eco2_abs->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_eco2_abs->SetMarkerStyle(8);
	g_idens_wp_eco2_abs->SetMarkerSize(2);
	g_idens_wp_eco2_abs->SetMarkerColor(kRed);
	mi_dens_wp_abs->Add(g_idens_wp_eco2_abs);

	TGraphErrors *g_idens_wp_eco3_abs = new TGraphErrors(vfilterup_eco3.size(),&vfilterup_eco3[0],&idens_eco3[0],NULL,&e_idens_eco3[0]);
	g_idens_wp_eco3_abs->GetXaxis()->SetTitle("ABS");
	g_idens_wp_eco3_abs->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	g_idens_wp_eco3_abs->SetMarkerStyle(8);
	g_idens_wp_eco3_abs->SetMarkerSize(2);
	g_idens_wp_eco3_abs->SetMarkerColor(kBlue);
	mi_dens_wp_abs->Add(g_idens_wp_eco3_abs);

	TCanvas *c_idens_wp_abs = new TCanvas();
	c_idens_wp_abs->cd();
	mi_dens_wp_abs->GetXaxis()->SetTitle("ABS");
	mi_dens_wp_abs->GetYaxis()->SetTitle("Current density [#muA/cm^{2}]");
	mi_dens_wp_abs->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");*/

	//Muon cluster size at working point as a function of rate
	//TGraphErrors *g_muon_cs_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&muon_clus_size_wp_std[0],NULL,&e_muon_clus_size_wp_std[0]);
	TGraphErrors *g_muon_cs_wp_std = new TGraphErrors(ratewp_std.size(),&dose[0],&muon_clus_size_wp_std[0],NULL,&e_muon_clus_size_wp_std[0]);
	g_muon_cs_wp_std->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_std->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_std->SetMarkerStyle(8);
	g_muon_cs_wp_std->SetMarkerSize(2);
	g_muon_cs_wp_std->SetMarkerColor(kBlack);
	mmuon_cs_wp->Add(g_muon_cs_wp_std);

	//TGraphErrors *g_muon_cs_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&muon_clus_size_wp_eco2[0],NULL,&e_muon_clus_size_wp_eco2[0]);
	TGraphErrors *g_muon_cs_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&muon_clus_size_wp_eco2[0],NULL,&e_muon_clus_size_wp_eco2[0]);
	g_muon_cs_wp_eco2->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_eco2->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco2->SetMarkerStyle(8);
	g_muon_cs_wp_eco2->SetMarkerSize(2);
	g_muon_cs_wp_eco2->SetMarkerColor(kRed);
	mmuon_cs_wp->Add(g_muon_cs_wp_eco2);

	//TGraphErrors *g_muon_cs_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&muon_clus_size_wp_eco3[0],NULL,&e_muon_clus_size_wp_eco3[0]);
	TGraphErrors *g_muon_cs_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&muon_clus_size_wp_eco3[0],NULL,&e_muon_clus_size_wp_eco3[0]);
	g_muon_cs_wp_eco3->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_eco3->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco3->SetMarkerStyle(8);
	g_muon_cs_wp_eco3->SetMarkerSize(2);
	g_muon_cs_wp_eco3->SetMarkerColor(kBlue);
	mmuon_cs_wp->Add(g_muon_cs_wp_eco3);

	//TGraphErrors *g_muon_cs_wp_std = new TGraphErrors(ratewp_std.size(),&ratewp_std[0],&muon_clus_size_wp_std[0],NULL,&e_muon_clus_size_wp_std[0]);
	TGraphErrors *g_muon_cs_wp_std_200 = new TGraphErrors(ratewp_std.size(),&dose[0],&muon_clus_size_wp_std_200[0],NULL,&e_muon_clus_size_wp_std[0]);
	g_muon_cs_wp_std_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_std_200->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_std_200->SetMarkerStyle(26);
	g_muon_cs_wp_std_200->SetMarkerSize(2);
	g_muon_cs_wp_std_200->SetMarkerColor(kBlack);
	mmuon_cs_wp->Add(g_muon_cs_wp_std_200);

	//TGraphErrors *g_muon_cs_wp_eco2 = new TGraphErrors(ratewp_eco2.size(),&ratewp_eco2[0],&muon_clus_size_wp_eco2[0],NULL,&e_muon_clus_size_wp_eco2[0]);
	TGraphErrors *g_muon_cs_wp_eco2_200 = new TGraphErrors(ratewp_eco2.size(),&dose[0],&muon_clus_size_wp_eco2_200[0],NULL,&e_muon_clus_size_wp_eco2[0]);
	g_muon_cs_wp_eco2_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_eco2_200->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco2_200->SetMarkerStyle(26);
	g_muon_cs_wp_eco2_200->SetMarkerSize(2);
	g_muon_cs_wp_eco2_200->SetMarkerColor(kRed);
	mmuon_cs_wp->Add(g_muon_cs_wp_eco2_200);

	//TGraphErrors *g_muon_cs_wp_eco3 = new TGraphErrors(ratewp_eco3.size(),&ratewp_eco3[0],&muon_clus_size_wp_eco3[0],NULL,&e_muon_clus_size_wp_eco3[0]);
	TGraphErrors *g_muon_cs_wp_eco3_200 = new TGraphErrors(ratewp_eco3.size(),&dose[0],&muon_clus_size_wp_eco3_200[0],NULL,&e_muon_clus_size_wp_eco3[0]);
	g_muon_cs_wp_eco3_200->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	g_muon_cs_wp_eco3_200->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco3_200->SetMarkerStyle(26);
	g_muon_cs_wp_eco3_200->SetMarkerSize(2);
	g_muon_cs_wp_eco3_200->SetMarkerColor(kBlue);
	mmuon_cs_wp->Add(g_muon_cs_wp_eco3_200);

	TCanvas *c_muon_cs_wp = new TCanvas();
	c_muon_cs_wp->cd();
	//mmuon_cs_wp->GetXaxis()->SetTitle("Rate [Hz/cm^{2}]");
	mmuon_cs_wp->GetXaxis()->SetTitle("Dose [#muS/h]");
	mmuon_cs_wp->GetYaxis()->SetTitle("Muon cluster size at fixed HV [strips]");
	mmuon_cs_wp->GetYaxis()->SetRangeUser(0.8,1.8);
	mmuon_cs_wp->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");

	//Muon cluster size at working point as a function of ABS
	/*TGraphErrors *g_muon_cs_wp_std_abs = new TGraphErrors(vfilterup_std.size(),&vfilterup_std[0],&muon_clus_size_wp_std[0],NULL,&e_muon_clus_size_wp_std[0]);
	g_muon_cs_wp_std_abs->GetXaxis()->SetTitle("ABS");
	g_muon_cs_wp_std_abs->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_std_abs->SetMarkerStyle(8);
	g_muon_cs_wp_std_abs->SetMarkerSize(2);
	g_muon_cs_wp_std_abs->SetMarkerColor(kBlack);
	mmuon_cs_wp_abs->Add(g_muon_cs_wp_std_abs);

	TGraphErrors *g_muon_cs_wp_eco2_abs = new TGraphErrors(vfilterup_eco2.size(),&vfilterup_eco2[0],&muon_clus_size_wp_eco2[0],NULL,&e_muon_clus_size_wp_eco2[0]);
	g_muon_cs_wp_eco2_abs->GetXaxis()->SetTitle("ABS");
	g_muon_cs_wp_eco2_abs->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco2_abs->SetMarkerStyle(8);
	g_muon_cs_wp_eco2_abs->SetMarkerSize(2);
	g_muon_cs_wp_eco2_abs->SetMarkerColor(kRed);
	mmuon_cs_wp_abs->Add(g_muon_cs_wp_eco2_abs);

	TGraphErrors *g_muon_cs_wp_eco3_abs = new TGraphErrors(vfilterup_eco3.size(),&vfilterup_eco3[0],&muon_clus_size_wp_eco3[0],NULL,&e_muon_clus_size_wp_eco3[0]);
	g_muon_cs_wp_eco3_abs->GetXaxis()->SetTitle("ABS");
	g_muon_cs_wp_eco3_abs->GetYaxis()->SetTitle("Muon cluster size [strips]");
	g_muon_cs_wp_eco3_abs->SetMarkerStyle(8);
	g_muon_cs_wp_eco3_abs->SetMarkerSize(2);
	g_muon_cs_wp_eco3_abs->SetMarkerColor(kBlue);
	mmuon_cs_wp_abs->Add(g_muon_cs_wp_eco3_abs);

	TCanvas *c_muon_cs_wp_abs = new TCanvas();
	c_muon_cs_wp_abs->cd();
	mmuon_cs_wp_abs->GetXaxis()->SetTitle("ABS");
	mmuon_cs_wp_abs->GetYaxis()->SetTitle("Muon cluster size [strips]");
	mmuon_cs_wp_abs->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");

	//Rate at WP as a function of ABS
	TGraphErrors *g_rate_wp_std = new TGraphErrors(vfilterup_std.size(),&vfilterup_std[0],&ratewp_std[0],NULL,&e_ratewp_std[0]);
	g_rate_wp_std->GetXaxis()->SetTitle("ABS");
	g_rate_wp_std->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	g_rate_wp_std->SetMarkerStyle(8);
	g_rate_wp_std->SetMarkerSize(2);
	g_rate_wp_std->SetMarkerColor(kBlack);
	mi_rate_wp->Add(g_rate_wp_std);

	TGraphErrors *g_rate_wp_eco2 = new TGraphErrors(vfilterup_eco2.size(),&vfilterup_eco2[0],&ratewp_eco2[0],NULL,&e_ratewp_eco2[0]);
	g_rate_wp_eco2->GetXaxis()->SetTitle("ABS");
	g_rate_wp_eco2->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	g_rate_wp_eco2->SetMarkerStyle(8);
	g_rate_wp_eco2->SetMarkerSize(2);
	g_rate_wp_eco2->SetMarkerColor(kRed);
	mi_rate_wp->Add(g_rate_wp_eco2);

	TGraphErrors *g_rate_wp_eco3 = new TGraphErrors(vfilterup_eco3.size(),&vfilterup_eco3[0],&ratewp_eco3[0],NULL,&e_ratewp_eco3[0]);
	g_rate_wp_eco3->GetXaxis()->SetTitle("ABS");
	g_rate_wp_eco3->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	g_rate_wp_eco3->SetMarkerStyle(8);
	g_rate_wp_eco3->SetMarkerSize(2);
	g_rate_wp_eco3->SetMarkerColor(kBlue);
	mi_rate_wp->Add(g_rate_wp_eco3);

	TCanvas *c_rate_wp = new TCanvas();
	c_rate_wp->cd();
	mi_rate_wp->GetXaxis()->SetTitle("ABS");
	mi_rate_wp->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
	mi_rate_wp->Draw("AP");
	eff_comp->Draw("SAME");
	ttitle->Draw("SAME");
	*/	

	//Eff curve comparison between different mixtures
	TCanvas *cEffCurveComp = new TCanvas();
	cEffCurveComp->cd();
	mEffMixComp->SetTitle("Eff - curves, source OFF");
	mEffMixComp->GetXaxis()->SetTitle("HV_{eff} [V]");
	mEffMixComp->GetYaxis()->SetTitle("#varepsilon [%]");
	mEffMixComp->Draw("AP");
	for (int i = 0; i < 3; i++) {
		fiteff1[i]->SetLineColor(i+1);
		fiteff1[i]->Draw("SAME");
	}
	lMix->Draw("SAME");
	ttitle->Draw("SAME");

	fout->cd(); //open .root file
	cdtof[0]->cd(); 
	cefficiency_std->Write("Eff(V)_abs");
	ccurr_std->Write("I(V)_abs");
	ccurrdens_std->Write("I_dens(V)_abs");
	cgamma_std->Write("Rate(V)_abs");
	cmuon_cs_std->Write("Muon_CS(V)_abs");
	cefficiency_rate_std->Write("Eff(V)_rate");
	ccurr_rate_std->Write("I(V)_rate");
	ccurrdens_std->Write("I_dens(V)_rate");
	cmuon_cs_rate_std->Write("Muon_CS(V)_rate");
	c->Write("Eff(v)_Idens(V)_curves_std");

	//fout->cd();
	cdtof[1]->cd();  //Enter the folder corresponding to the current HV point
	cefficiency_eco2->Write("Eff(V)_abs");
	ccurr_eco2->Write("I(V)_abs");
	ccurrdens_eco2->Write("I_dens(V)_abs");
	cgamma_eco2->Write("Rate(V)_abs");
	cmuon_cs_eco2->Write("Muon_CS(V)_abs");
	cefficiency_rate_eco2->Write("Eff(V)_rate");
	ccurr_rate_eco2->Write("I(V)_rate");
	ccurrdens_eco2->Write("I_dens(V)_rate");
	cmuon_cs_rate_eco2->Write("Muon_CS(V)_rate");
	c1->Write("Eff(v)_Idens(V)_curves_eco2");

	//fout->cd();
	cdtof[2]->cd();  //Enter the folder corresponding to the current HV point
	cefficiency_eco3->Write("Eff(V)_abs");
	ccurr_eco3->Write("I(V)_abs");
	ccurrdens_eco3->Write("I_dens(V)_abs");
	cgamma_eco3->Write("Rate(V)_abs");
	cmuon_cs_eco3->Write("Muon_CS(V)_abs");
	cefficiency_rate_eco3->Write("Eff(V)_rate");
	ccurr_rate_eco3->Write("I(V)_rate");
	ccurrdens_eco3->Write("I_dens(V)_rate");
	cmuon_cs_rate_eco3->Write("Muon_CS(V)_rate");
	c2->Write("Eff(v)_Idens(V)_curves_eco3");

	fout->cd();
	cEffCurveComp->Write("Eff_curve(V)_std_eco2_eco3");
	c_eff_wp->Write("Eff(WP)_std_eco2_eco3_dose");
	c_muon_cs_wp->Write("Muon_cs(WP)_std_eco2_eco3_dose");
	c_i_wp->Write("I_wp_rate_std_eco2_eco3");

	fout->Close();
}