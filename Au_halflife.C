#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

using std::cout;
using std::cerr;
using std::endl;

void Au_halflife(int run_num=0){

    //Variables
    const int num_det = 2;

    double range[2] = {1, 59};
    double startrange[2] = {0, 15};
    double endrange[2] = {30, 60};
    TH1F *hHalflife = new TH1F("hHalflife", "hHalflife", 2*150, 0, 15);
    double t12;
    double mean = 0;
    double counter = 0;
    
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

    //TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1])", range[0], range[1]);
    TF1 *fDecay = new TF1("fDecay", "[0]*2^(-1*x/[1]) + [2]", range[0], range[1]);
    fDecay->SetParameter(0, 5000);
    fDecay->SetParameter(1, 7);
    fDecay->SetParameter(2, 0);

    TCanvas *cDecay= new TCanvas("cDecay", "cDecay", 800, 500);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);

    hNet->SetTitle("");
    hNet->Rebin(10);
    hNet->SetMarkerColor(2);
    hNet->SetMarkerStyle(20);
    hNet->SetMarkerSize(1.2);
    hNet->GetXaxis()->SetTitle("Time since last cycle start [s]");
    hNet->GetYaxis()->SetTitle("Counts/s");
    hNet->GetXaxis()->SetTitleSize(0.06);
    hNet->GetYaxis()->SetTitleSize(0.06);

    for (int sr=startrange[0]; sr<startrange[1]; sr++){
        for (int er=endrange[0]; er<endrange[1]; er++){
            hNet->GetXaxis()->SetRangeUser(sr, er);
            fDecay->SetRange(sr, er);
            hNet->Fit("fDecay");
            t12 = fDecay->GetParameter(1);
            hHalflife->Fill(t12);
            mean+=t12;
            counter++;
        }
    }
    hNet->Draw("e");
    mean = mean/counter;

    new TCanvas();
    hHalflife->Rebin(1);
    hHalflife->Draw();
    hHalflife->Fit("gaus");
    cout << mean << endl;
    
    
}
