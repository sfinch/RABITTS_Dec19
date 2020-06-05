
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/processed.h"

using std::cout;
using std::cerr;
using std::endl;

void Zr_plot(int run_num=3233){

    //Variables
    //cuts
    double range[2] = {4.2, 11.0};

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

    TH1F *hRes = (TH1F*)hNet->Clone();

    TCanvas *c1 = new TCanvas("c1", "Zr Energy", 800, 500);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetLogy();
    hEnCount->Rebin(2*5);
    hEnIrr->Rebin(2*5);
    hEnCount->GetXaxis()->SetRangeUser(2100, 2700);
    hEnCount->SetLineColor(2);
    hEnCount->SetLineWidth(2);
    hEnIrr->SetLineWidth(2);
    hEnCount->SetTitle("");
    hEnCount->GetXaxis()->SetTitle("E_{#gamma} [keV]");
    hEnCount->GetYaxis()->SetTitle("Counts / 0.5 keV");
    hEnCount->GetXaxis()->SetTitleSize(0.06);
    hEnCount->GetYaxis()->SetTitleSize(0.06);
    hEnCount->GetXaxis()->SetLabelSize(0.05);
    hEnCount->GetYaxis()->SetLabelSize(0.05);

    hEnCount->Draw();
    hEnIrr->Draw("same");

    TLegend *leg = new TLegend(0.8, 0.8, 0.95, 0.95);
    leg->SetLineColor(0);
    leg->AddEntry(hEnIrr, "Irradiation, beam on", "l");
    leg->AddEntry(hEnCount, "Counting, beam off", "l");
    leg->SetMargin(0.15);
    //leg->Draw("same");

    TCanvas *c2 = new TCanvas("c2", "Zr decay", 600, 800);
    c2->Divide(1,2);

    TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1])", range[0], range[1]);
    //TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1]) + [2]", range[0], range[1]);
    //TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1]) + [2]*2^(-1*x/[3])", range[0], range[1]);
    fDecay->SetParameter(0, 1000);
    fDecay->SetParameter(1, 0.8);
    fDecay->SetParameter(2, 0.);
    fDecay->SetParameter(3, 1.);
    
    c2->cd(1);
    //hROI->Rebin(10);
    hROI->GetXaxis()->SetRangeUser(range[0], range[1]);
    hROI->Draw("e");
    hROI->SetTitle("");
    hROI->GetXaxis()->SetTitle("Cycle time [s]");
    hROI->GetYaxis()->SetTitle("Gross counts / 0.1 s");
    hROI->GetXaxis()->SetTitleSize(0.06);
    hROI->GetYaxis()->SetTitleSize(0.06);
    hROI->GetXaxis()->SetLabelSize(0.05);
    hROI->GetYaxis()->SetLabelSize(0.05);
    
    c2->cd(2);
    //hNet->Rebin(10);
    hNet->GetXaxis()->SetRangeUser(range[0], range[1]);
    hNet->Draw("e");
    hNet->Fit("fDecay");
    hNet->SetTitle("");
    hNet->GetXaxis()->SetTitle("Time since last cycle start [s]");
    hNet->GetYaxis()->SetTitle("Counts / 0.1 s");
    hNet->GetXaxis()->SetTitleSize(0.06);
    hNet->GetYaxis()->SetTitleSize(0.06);
    hNet->GetXaxis()->SetLabelSize(0.05);
    hNet->GetYaxis()->SetLabelSize(0.05);

    TPaveText *l = new TPaveText(3.5, 900, 7.5, 1300);
    l->Draw("same");
    l->SetTextAlign(11);
    l->AddText("Measured T_{1/2} = 787 #pm 8 ms");
    l->AddText("NNDC T_{1/2} = 809 #pm 2 ms");
    l->SetFillColor(0);
    l->Draw("same");

    TCanvas *c3 = new TCanvas("c3", "Residuals", 600, 800);
    c3->Divide(1,2);
    c3->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    hNet->GetXaxis()->SetRangeUser(range[0], range[1]);
    hNet->Draw("e");
    hNet->Fit("fDecay");
    hNet->SetTitle("");
    hNet->GetXaxis()->SetTitle("Cycle time [s]");
    hNet->GetYaxis()->SetTitle("Net counts / 0.1 s");
    hNet->GetXaxis()->SetTitleSize(0.06);
    hNet->GetYaxis()->SetTitleSize(0.06);
    hNet->GetXaxis()->SetLabelSize(0.05);
    hNet->GetYaxis()->SetLabelSize(0.05);
    l->Draw("same");

    c3->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    double res, val;
    for (int i=0; i<hRes->GetXaxis()->GetNbins(); i++){
        val =fDecay->Eval(hRes->GetBinCenter(i));
        res = hRes->GetBinContent(i) - val;
        //cout << hRes->GetBinCenter(i) << " " << val << " " << res << endl;
        hRes->SetBinContent(i, res);
        hRes->SetBinError(i, hNet->GetBinError(i));
    } 

    hRes->GetXaxis()->SetRangeUser(range[0], range[1]);
    hRes->Draw("e");
    hRes->SetTitle("");
    hRes->GetXaxis()->SetTitle("Cycle time [s]");
    hRes->GetYaxis()->SetTitle("Fit residuals / 0.1 s");
    hRes->GetXaxis()->SetTitleSize(0.06);
    hRes->GetYaxis()->SetTitleSize(0.06);
    hRes->GetXaxis()->SetLabelSize(0.05);
    hRes->GetYaxis()->SetLabelSize(0.05);

    TCanvas *c4 = new TCanvas("c4", "Net", 800, 500);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    hNet->GetXaxis()->SetRangeUser(range[0], range[1]);
    hNet->Draw("e");
    hNet->Fit("fDecay");
    hNet->SetTitle("");
    hNet->GetXaxis()->SetTitle("Time since last cycle start [s]");
    hNet->GetYaxis()->SetTitle("Net counts / 0.1 s");
    hNet->GetXaxis()->SetTitleSize(0.06);
    hNet->GetYaxis()->SetTitleSize(0.06);
    hNet->GetXaxis()->SetLabelSize(0.05);
    hNet->GetYaxis()->SetLabelSize(0.05);

}
