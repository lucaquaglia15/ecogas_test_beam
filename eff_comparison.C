#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TMath.h"
#include "TDatime.h"
#include "TMarker.h"
#include <string>
#include <vector>

TF1 *fit_eff = new TF1("fit_eff","[0]/(1+TMath::Exp(-[1]*(x-[2])))",8000,10600);

using namespace std;
using namespace TMath;

void eff_comparison() {

	string fol = "/home/luca/Desktop/PhD/TB_october_2021/runs/std/ALICEX"; 
	gSystem->cd(fol.c_str());

	fit_eff->SetParLimits(1,0,0.1);
	fit_eff->SetParLimits(2,8000,10600);
	double eff_max, eff_fifty, eff_95, tens_eff_95, tens_fifty, lambda, wp, eff_wp;

	vector <double> efficiency_wp, efficiency_wp_rates;
	//vector <double> rates_wp{13,150,120,63};
	//vector <double> i_wp{1.78,6.13,10.12,15.82};
	//vector <double> cluster_rate{9.22,43.33,85.11,109.49};
	
	int filters_eff[15] = {0,1,2,3,4,6,10,22,33,46,69,100,220,460,1000};
	//double filters_eff[15] = {(int)0,(int)1,(int)2.2,3.3,4.6,6.9,(int)10,(int)22,(int)33,(int)46,(int)69,(int)100,(int)220,(int)460,(int)1000};

	TMultiGraph *m = new TMultiGraph();
	m->GetXaxis()->SetTitle("HV_{eff}");
	m->GetYaxis()->SetTitle("Efficiency");
	m->GetXaxis()->SetLimits(8600,10600);

	TGraphErrors *eff_filtri[15];

	TLegend *l = new TLegend();
	TLegend *l_rate = new TLegend();

	for (int i = 0; i < 15; i++) {
		string filtro = to_string(filters_eff[i]) + ".txt";

		cout << filtro << endl;

		ifstream hfile(filtro.c_str());

		vector <double> HV, eff, e_eff, cs, e_cs, rndm1, rndm2;
		double high_volt, efficiency, e_efficiency, cluster_size, e_cluster_size, random1, random2;

		while (hfile >> high_volt >> efficiency >> e_efficiency >> cluster_size >> e_cluster_size >> random1 >> random2) {
			HV.push_back(high_volt);
			eff.push_back(efficiency);
			e_eff.push_back(e_efficiency);
			cs.push_back(cluster_size);
			e_cs.push_back(e_cluster_size);
			rndm1.push_back(random1);
			rndm2.push_back(random2);
		}

	 cout << "Filtro: " +to_string(filters_eff[i]) << endl;
	 for (unsigned int i = 0; i < HV.size(); i++) {
	 	cout << "HV: " << HV[i] << " eff: " << eff[i] << endl;
	 }	

	 eff_filtri[i] = new TGraphErrors(HV.size(), &HV[0], &eff[0], NULL, &e_eff[0]);
	 eff_filtri[i]->GetXaxis()->SetTitle("HV_{eff}");
	 eff_filtri[i]->GetYaxis()->SetTitle("Efficiency");
	 eff_filtri[i]->GetXaxis()->SetLimits(8600,10600);
	 eff_filtri[i]->SetTitle(filtro.c_str());
	 eff_filtri[i]->SetMarkerStyle(8);
	 eff_filtri[i]->SetMarkerSize(1);
	 eff_filtri[i]->SetMarkerColor(27+i+1); //To get more colors
	 fit_eff->SetLineColor(27+i+1);
	 eff_filtri[i]->Fit(fit_eff,"RM+");
	 m->Add(eff_filtri[i]);

	 if (filters_eff[i] == 0) l->AddEntry(eff_filtri[i],"Source OFF","p");
	 else if (filters_eff[i] == 2) l->AddEntry(eff_filtri[i],"2.2","p");
	 else if (filters_eff[i] == 3) l->AddEntry(eff_filtri[i],"3.3","p");
	 else if (filters_eff[i] == 4) l->AddEntry(eff_filtri[i],"4.6","p");
	 else if (filters_eff[i] == 6) l->AddEntry(eff_filtri[i],"6.9","p");
	 else l->AddEntry(eff_filtri[i],to_string(filters_eff[i]).c_str(),"p"); 

	 new TCanvas();
	 eff_filtri[i]->Draw("AP");
	 //eff_filtri[i]->Draw("SAME");
	 l->Draw("SAME");

	 //Get sigmoid values
	 //eff_max = 1/fit_eff->GetParameter(0);
	 eff_max = fit_eff->GetParameter(0);
	 lambda = fit_eff->GetParameter(1);
	 tens_fifty = fit_eff->GetParameter(2);
	 eff_fifty = eff_max/2;
	 eff_95 = 0.95*eff_max;
	 tens_eff_95 = fit_eff->GetX(eff_95, 8000, 10600);

	 cout << "Eff max: " << eff_max << ", eff 50%: " << eff_fifty << ", eff 95%: " << eff_95 << ", voltage for 95% eff: " << tens_eff_95 << ", lambda: " << lambda << endl;

	 wp = (Log(19)/lambda) + tens_fifty +150;

	 cout << "Working point: " << wp << endl;
	 eff_wp = fit_eff->Eval(wp);
	 efficiency_wp.push_back(eff_wp);
	 /*if (filters_eff[i] == 0 || filters_eff[i] == 69 || filters_eff[i] == 100 || filters_eff[i] == 220) {
	 	if (filters_eff[i] == 0) l_rate->AddEntry(eff_filtri[i],"Source OFF","p");
	    else l_rate->AddEntry(eff_filtri[i],to_string(filters_eff[i]).c_str(),"p");
	 	efficiency_wp_rates.push_back(eff_wp);
	 }*/

	 HV.clear();
	 eff.clear();
	 e_eff.clear();
	 cs.clear();
	 e_cs.clear();
	 rndm1.clear();
 	 rndm2.clear();
	}

	new TCanvas();
	m->Draw("AP");
	l->Draw("SAME");

	/*
	new TCanvas();
	TGraphErrors *eff_rate = new TGraphErrors(rates_wp.size(),&rates_wp[0],&efficiency_wp_rates[0],NULL,NULL);
	eff_rate->GetXaxis()->SetTitle("Rate [HZ/cm^{2}]");
	eff_rate->GetYaxis()->SetTitle("Efficiency");
	eff_rate->SetMarkerStyle(8);
	eff_rate->SetMarkerSize(2);
	eff_rate->GetYaxis()->SetRange(0.,1.);
	eff_rate->Draw("AP");

	TMarker *mark;
   	for (int j=0; j<4; j++) {
      eff_rate->GetPoint(j,rates_wp[j],efficiency_wp_rates[j]);
      mark = new TMarker(rates_wp[j],efficiency_wp_rates[j],8);
      mark->SetMarkerSize(2);
      mark->SetMarkerColor(j+1);
      mark->Draw("SAME");
    }
	l_rate->Draw("SAME");

	new TCanvas();
	TGraphErrors *i_rate = new TGraphErrors(i_wp.size(),&cluster_rate[0],&i_wp[0],NULL,NULL);
	i_rate->GetXaxis()->SetTitle("Cluster rate [HZ/cm^{2}]");
	i_rate->GetYaxis()->SetTitle("Current at WP [#muA]");
	i_rate->SetMarkerStyle(8);
	i_rate->SetMarkerSize(2);
	i_rate->Draw("AP");
	*/
}