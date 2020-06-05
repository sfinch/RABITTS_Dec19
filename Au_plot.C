
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

using std::cout;
using std::cerr;
using std::endl;

void Au_plot(int run_num=0){

    //Variables
    const int num_det = 2;

    //double range[2] = {0, 60};
    double range[2] = {1, 59};
    

    //in file
    TFile *fHist = new TFile(Form("data_hist/RABBIT_%i.root", run_num));

    TH1F *hEnIrr   = (TH1F*)fHist->Get("hEnIrr");
    TH1F *hEnCount = (TH1F*)fHist->Get("hEnCount");
    TH1F *hEnCoinc = (TH1F*)fHist->Get("hEnCoinc");

    TH1F *hCycle = (TH1F*)fHist->Get("hCycle");
    TH1F *hROI = (TH1F*)fHist->Get("hROI");
    TH1F *hBkgdL = (TH1F*)fHist->Get("hBkgdL");
    TH1F *hBkgdR = (TH1F*)fHist->Get("hBkgdR");
    TH1F *hNet = (TH1F*)fHist->Get("hNet");

    //cuts
    double time_irr[2] = {0, 20};
    double time_count[2] = {21, 81};
    double ROI[2] = {275, 285};
    double BkgdL[2] = {270, 275};
    double BkgdR[2] = {285, 290};

    //TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1])", range[0], range[1]);
    TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1]) + [2]", range[0], range[1]);
    fDecay->SetParameter(0, 5000);
    fDecay->SetParameter(1, 7);
    fDecay->SetParameter(2, 0);

    TCanvas *cSpec = new TCanvas("cSpec", "cSpec", 800, 500);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetLogy();
    hEnCount->SetTitle("");
    hEnIrr->SetTitle("");
    hEnCount->SetLineColor(2);
    hEnCount->Rebin(10);
    hEnIrr->Rebin(10);
    hEnIrr->GetXaxis()->SetRange(100,600);
    hEnCount->GetXaxis()->SetRange(100,600);
    hEnIrr->SetMinimum(100);
    hEnIrr->SetMaximum(40000);
    hEnIrr->GetXaxis()->SetTitle("E_{#gamma} [keV]");
    hEnIrr->GetYaxis()->SetTitle("Counts/keV");
    hEnIrr->GetXaxis()->SetTitleSize(0.06);
    hEnIrr->GetYaxis()->SetTitleSize(0.06);
    hEnIrr->SetLineWidth(2);
    hEnCount->SetLineWidth(2);

    hEnIrr->Draw();
    hEnCount->Draw("same");

    TCanvas *cDecay= new TCanvas("cDecay", "cDecay", 800, 500);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);

    hNet->SetTitle("");
    hNet->Rebin(10);
    hNet->GetXaxis()->SetRangeUser(range[0], range[1]);
    hNet->SetMarkerColor(2);
    hNet->SetMarkerStyle(20);
    hNet->SetMarkerSize(1.2);
    hNet->GetXaxis()->SetTitle("Time since last cycle start [s]");
    hNet->GetYaxis()->SetTitle("Counts/s");
    hNet->GetXaxis()->SetTitleSize(0.06);
    hNet->GetYaxis()->SetTitleSize(0.06);
    hNet->Draw("e");
    hNet->Fit("fDecay");
    
    //new TCanvas();
    //hROI->Rebin(10);
    //hROI->GetXaxis()->SetRange(range[0], range[1]);
    //hROI->Draw("e");
    //hROI->Fit("fDecay");
    
    //new TCanvas();
    //hEnCoinc->SetTitle("Given 279 keV in another crystal");
    //hEnCoinc->Rebin(10);
    //hEnCoinc->Draw();
    
    
    
}
