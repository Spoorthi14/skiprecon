#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include <iostream>
#include <stdio.h>

double sigma;

TH1F *spec;
TH1F *reduced;
TH1F *model;

void loaddata(TString infile){

    TFile *f = TFile::Open(infile,"READ");
    spec = (TH1F*) f->Get("spec_hist");
    spec->SetDirectory(0);
    double minVal = TMath::Max(-10.5,spec->GetBinLowEdge(1));
    double bw = spec->GetBinWidth(1);
    double maxVal = 50.5; //Hardcoded - careful
    spec->GetXaxis()->SetRangeUser(minVal,maxVal);

    // We make a reduced spec. so as not to fit 20K ADU spectrum
    reduced = new TH1F("reduced","reduced", int((maxVal-minVal)/bw)+1,minVal,maxVal);

    int firstbin = spec->FindBin(minVal);
    int lastbin = spec->FindBin(maxVal);

    int cnt = 1;
    for (int i=firstbin; i<=lastbin; i++){
        reduced->SetBinContent(cnt,spec->GetBinContent(i));
        cnt++;
    }

    TFitResultPtr ft = reduced->Fit("gaus","SQ0NR","",-4,4);
    sigma = ft->Parameter(2);

    model = (TH1F*) reduced->Clone();
    model->SetDirectory(0);

}

void fillmodel(double A, double l, double o, double cb, bool toplot){

    model->Reset();

    //Standard Histogram Parameters
    double start = model->GetBinLowEdge(1);
    double bw    = model->GetBinWidth(1);
    double nbins = model->GetNbinsX();

    for(int k=1; k<=nbins; k++){
        double v = 0;
        double he = start + k*bw;
        double le = he - bw;
             
        int j=0;
        while(TMath::Poisson(j,l)>1.E-9){
            double uv = j*cb + o;
            v += A*TMath::Poisson(j,l)*0.5*(TMath::Erf((he-uv)/sigma/TMath::Sqrt2()) - TMath::Erf((le-uv)/sigma/TMath::Sqrt2()) );
            j++;
        }
    
        model->SetBinContent(k,v);
    }

    if (toplot){
        model->SetLineColor(2);
        model->SetLineWidth(3);
        model->SetLineStyle(2);

        TCanvas *cv = new TCanvas("cv","cv",700,500);
        cv->cd();
        reduced->Draw("H");
        model->Draw("H same");
        cv->SetLogy();
    }

}

double LogLikelihood(const double *xx){

    double A = xx[0];
    double l = xx[1];
    double o = xx[2];
    double cb = xx[3];
    
    //Standard Histogram Parameters
    double start = reduced->GetBinLowEdge(1);
    double bw    = reduced->GetBinWidth(1);
    double nbins = reduced->GetNbinsX();

    double ll = 0;
    
    fillmodel(A, l, o, cb, false);

    for(int k=1; k<=nbins; k++){
    
        double bc = reduced->GetBinContent(k);
        double v = model->GetBinContent(k);
    
        if(v==0 && bc==0) continue;
        if(v==0 && bc>0) ll += 1e9;
        else
            ll += bc*TMath::Log(v) - v - TMath::LnGamma(bc+1);
    }
    return -ll;
}

 // Takes in file processed by skipAnalyzer
 // run e.g. darkFit("1631_spectrum_2000.root", 22.5, 5290) *calib in ADU/e-*
void darkFit(TString infile, double sec_exposure){

    loaddata(infile);

    // Minimizer + techincal junk
    ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    ROOT::Math::Functor fn(&LogLikelihood,4);
    mini->SetMaxFunctionCalls(100000); 
    mini->SetTolerance(1E-7);
    mini->SetPrintLevel(0);
    mini->SetErrorDef(0.5); 
    mini->SetFunction(fn);
    
    mini->SetFixedVariable(0 , "A",     reduced->Integral());  
    mini->SetLimitedVariable(1,	"l",	0.3, 1e-6, 0., 10.); // in e-/pix/img
    mini->SetLimitedVariable(2, "o",    0, 1e-4, -10., 10.); // offset ADU
    mini->SetLimitedVariable(3, "cb", 12, 1e-3, 2, 50.); // ADU/e-
    
    mini->Minimize();
    const double *xscan = mini->X();
    
    double Lambda = xscan[1]*(86400./sec_exposure); 
    cout << "DC rate: " << Lambda << " e-/pix/day" << endl;
    cout << "Offset: " << xscan[2] << " ADU" << endl;
    cout << "Calib: " << xscan[3] << " ADU/e-" << endl;

    // Plot
    fillmodel(xscan[0], xscan[1], xscan[2], xscan[3], true);


}
