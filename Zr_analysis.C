
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "include/processed.h"

using std::cout;
using std::cerr;
using std::endl;

void Zr_analysis(){

    //Variables
    int run_num = 33;
    const int num_det = 2;

    TH1F *hEnIrr   = new TH1F("hEnIrr", "hEnIrr", 30000, 10, 3010);
    TH1F *hEnCount = new TH1F("hEnCount", "hEnCount", 30000, 10, 3010);
    TH1F *hEnCoinc = new TH1F("hEnCoinc", "hEnCoinc", 30000, 10, 3010);

    TH1F *hCycle = new TH1F("hCycle", "hCycle", 1000, -10, 90);
    TH1F *hROI = new TH1F("hROI", "hROI", 1000, -10, 90);
    TH1F *hBkgdL = new TH1F("hBkgdL", "hBkgdL", 1000, -10, 90);
    TH1F *hBkgdR = new TH1F("hBkgdR", "hBkgdR", 1000, -10, 90);
    TH1F *hNet = new TH1F("hNet", "hNet", 1000, -10, 90);

    //cuts
    double time_irr[2] = {0.0, 3.0};
    double time_count[2] = {4.0, 12.0};

    double ROI[2] = {2308, 2322};
    double BkgdL[2] = {2300, 2307};
    double BkgdR[2] = {2322, 2329};

    //double ROI[2] = {2179, 2191};
    //double BkgdL[2] = {2169, 2175};
    //double BkgdR[2] = {2194, 2200};
    
    //double ROI[2] = {655, 665};
    //double BkgdL[2] = {645, 650};
    //double BkgdR[2] = {670, 675};

    //in file
    processed rabbit(run_num);

    //out file
    TFile *fHist = new TFile(Form("data_hist/RABBIT_%i.root", run_num), "RECREATE");

    //loop over data
    Long64_t nentries = rabbit.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        nb = rabbit.GetEntry(jentry);   nbytes += nb;
        if (jentry%100000==0){
            cout << '\r' << "Processing event " << jentry;
        }
        hCycle->Fill(rabbit.cycle_time);

        //irradiation time
        if ((rabbit.cycle_time>time_irr[0]) && (rabbit.cycle_time<time_irr[1])){
            for (int det=0; det<num_det; det++){
                hEnIrr->Fill(rabbit.En[det]);            
            }
        }
        //counting time
        else if ((rabbit.cycle_time>time_count[0]) && (rabbit.cycle_time<time_count[1])){
            for (int det=0; det<num_det; det++){
                hEnCount->Fill(rabbit.En[det]);            

                //look for coincidences
                //if ((rabbit.En[det]>ROI[0]) && (rabbit.En[det]<ROI[1])){
                //    for (int coinc=0; coinc<num_det; coinc++){
                //        if (coinc!=det){
                //            hEnCoinc->Fill(rabbit.En[coinc]);
                //        }
                //    }
                //}//end coinc
            }
        }

        ///ROI cut
        for (int det=0; det<num_det; det++){
            if ((rabbit.En[det]>ROI[0]) && (rabbit.En[det]<ROI[1])){
                hROI->Fill(rabbit.cycle_time);
            }
            else if ((rabbit.En[det]>BkgdL[0]) && (rabbit.En[det]<BkgdL[1])){
                hBkgdL->Fill(rabbit.cycle_time);
            }
            else if ((rabbit.En[det]>BkgdR[0]) && (rabbit.En[det]<BkgdR[1])){
                hBkgdR->Fill(rabbit.cycle_time);
            }
        }

    }

    hNet->Add(hROI);
    hNet->Add(hBkgdL, -1);
    hNet->Add(hBkgdR, -1);

    //write histos to file
    fHist->cd();

    hCycle->Write();
    hEnIrr->Write();
    hEnCount->Write();
    hEnCoinc->Write();
    hROI->Write();
    hBkgdL->Write();
    hBkgdR->Write();
    hNet->Write();

    fHist->Write();
    fHist->Close();
    
}
