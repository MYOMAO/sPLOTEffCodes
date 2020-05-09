#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "TCanvas.h"



void DrawWeight(){

	
	
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);





	TFile * fin = new TFile("DataWeight/sPLOTEffWeight0_90.root");
	
	fin->cd();




	TTree * t = (TTree * ) fin->Get("EffWeightTreesPLOT");

	TH1D * EffSplotWeighthis = new TH1D("EffSplotWeighthis","",50,-0.5,1.5);
	EffSplotWeighthis->GetXaxis()->SetTitle("sPlot Weight");
	EffSplotWeighthis->GetYaxis()->SetTitle("B^{0}_{s} Candidates Counts");
	EffSplotWeighthis->GetYaxis()->SetTitleOffset(1.4);
	EffSplotWeighthis->GetXaxis()->CenterTitle();
	EffSplotWeighthis->GetYaxis()->CenterTitle();


	t->Project("EffSplotWeighthis","EffSplotWeight");



	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	EffSplotWeighthis->Draw();
	c->SaveAs("sPLOTEffDistribution.png");


}
