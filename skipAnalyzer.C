#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TParameter.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include <iostream>
#include <stdio.h>

// Takes in infile processed by skipRead.C
//e.g. skipAnalyzer("1631.root", 5, 965, 40)
void skipAnalyzer(TString infile, int iskipstart, int iskipend, double sguess){

    double maxVal = 25000.5;  // Hardcoded approximation of image range im ADU
    double minVal = -0.5;
    int nbins = 5*(maxVal-minVal+1); // X bins per ADU

    double zeroadu; // Location of 0 pix peak
    int zeroind;
    
    int imStart = 1;
    int imEnd = 1000;
    int ovscStart = 1040;
    int ovscEnd = 1100;
    int nX = imEnd-imStart+1;
    int nLast;

    // Get TTree & finfo (CCDDrone should save ovsc. info at some point..)
    TFile *f = TFile::Open(infile,"READ");
    TParameter<double>* ncP = NULL;
    TParameter<double>* nrP = NULL;
    TParameter<double>* ndP = NULL;
    TTree *info = (TTree*) f->Get("finfo");
    info->SetBranchAddress("NAXIS1",&ncP);
    info->SetBranchAddress("NAXIS2",&nrP);
    info->SetBranchAddress("NDCMS",&ndP);
    info->GetEntry(0);
    int ncolumns = (int) ncP->GetVal();
    int nrows = (int) nrP->GetVal();
    int ndcms = (int) ndP->GetVal();

    nLast = ncolumns/ndcms; // Tottal number of columns
    nX = imEnd-imStart+1; // Num. columns in image under consideration

    // Get Image "cube of size (columns, rows, skips)"    
    TH3F* cube = (TH3F*) f->Get("image_cube");
    cube->SetDirectory(0); // break link with file/parent

    vector<double> skipsV;  
    vector<double> sigmaV; // These contain the 1/sqrt(N) graph
    vector<double> idealV;
    
    // Compute average image & further analyze %100 images
    TH2F* image_avg = new TH2F("image_avg","image_avg",nX,0,nX,nrows,0,nrows); // to store final average
    TH2F* running_img; // to store %100 etc. images
    TH1F* running_hist = new TH1F("image_hist","image_hist", nbins, minVal, maxVal); // to store histograms as we proceed

    int cnt = 1; // where we are in the loop
    for (int k = iskipstart; k<=iskipend; k++){
        for (int i = imStart; i<=imEnd; i++){
            for (int j = 4; j<=nrows-1; j++){
                image_avg->SetBinContent(i,j,image_avg->GetBinContent(i,j)+cube->GetBinContent(i,j,k));
            }
        }
        if (cnt == 1 || cnt == 10 || cnt%100 == 0){
            cout << "Processing Skip: " << cnt << endl;
            running_hist->Reset();
            skipsV.push_back(cnt); 
            running_img = (TH2F*) image_avg->Clone();
            running_img->SetDirectory(0); // Breaks pointer link between image_avg and running_img
            running_img->Scale(1./cnt); // This is now the average up to that skip
            for (int ii = imStart; ii<=imEnd; ii++){
                for (int jj = 1; jj<=nrows; jj++){
                    running_hist->Fill(running_img->GetBinContent(ii,jj)); // 1D spectrum populated
                }
            }
            // Find bin with highest content and *assume* it corresponds to 0 peak
            running_hist->GetXaxis()->SetRangeUser(2,maxVal); // Don't include saturation 0 ADU line for next step
            zeroind = running_hist->GetMaximumBin();
            zeroadu = running_hist->GetBinCenter(zeroind); // ADU value
         
            // Fit the largest peak we found
            double fitR = (cnt==1) ? sguess:(idealV[0]/TMath::Sqrt(cnt)); //+- sigma to fit around zeroadu. If first, use guess, otherwise use ideal
            //S: store, Q: quiet, 0N: don't draw, R: use specified range
            TFitResultPtr ft = running_hist->Fit("gaus","SQ0NR","",zeroadu-fitR, zeroadu+2.*fitR); //Fit gaussian, weighting + more than -
            double sigmaFit = ft->Parameter(2);
            sigmaV.push_back(sigmaFit);
            (cnt==1) ? idealV.push_back(sigmaFit):idealV.push_back(idealV[0]/TMath::Sqrt(cnt));
        }
        cnt+=1;
    }
    image_avg->Scale(1./(cnt-1)); // contains total average

    // Populate 1D
    TH1F* avg_hist = new TH1F("avg_hist","avg_hist", nbins, minVal, maxVal);
    for (int i = imStart; i<=imEnd; i++){
        for (int j = 4; j<=nrows-1; j++){
            avg_hist->Fill(image_avg->GetBinContent(i,j));
        }
    }
    // Flip it to 
    avg_hist->GetXaxis()->SetRangeUser(2,maxVal); // Don't include saturation 0ADu
    zeroind = avg_hist->GetMaximumBin();
    zeroadu = avg_hist->GetBinCenter(zeroind); // Repeat same process as above loop
    TH1F* spec_hist = new TH1F("spec_hist","Single e^{-} resolution spectrum",nbins, zeroadu-maxVal, zeroadu-minVal); // Final spectrum
    spec_hist->GetXaxis()->SetTitle("Zeroed Pixel Value [ADU]");
    spec_hist->GetYaxis()->SetTitle("Count [/0.2 ADU]");
    cout << "Zero ADU: " << zeroadu << endl;
    for (int i = 1; i<=nbins; i++){
        spec_hist->SetBinContent(i,avg_hist->GetBinContent(nbins-i+1));
    }
    spec_hist->GetXaxis()->SetRangeUser(-50.,100.); // To zoom in on spectrum when plotting (check negative side to make sure we're not hiding peaks)

    // Create 1/sqrt(N) graph
    TGraph *sqrtN = new TGraph(skipsV.size(),&skipsV[0],&sigmaV[0]);
    TGraph *ideal = new TGraph(skipsV.size(),&skipsV[0],&idealV[0]);
    sqrtN->GetXaxis()->SetTitle("Skip Number");
    sqrtN->GetYaxis()->SetTitle("#sigma_{ADU}");
    sqrtN->SetTitle(Form("Sigma vs. Skips [%i NDCM]",ndcms));

    // Plots
    gStyle->SetOptStat(0);
    TCanvas *cv = new TCanvas("cv","cv",700,500);
    cv->cd();
    spec_hist->SetLineWidth(2);
    spec_hist->SetLineColor(4);
    spec_hist->Draw("H");
    cv->SetLogy();

    TCanvas *cv2 = new TCanvas("cv2","cv2",700,500);
    cv2->cd();
    sqrtN->SetMarkerStyle(20);
    sqrtN->SetLineColor(2);
    sqrtN->SetLineWidth(2);
    ideal->SetLineColor(1);
    ideal->SetLineStyle(2);
    ideal->SetLineWidth(2);
    sqrtN->Draw("ALP");
    ideal->Draw("L same");
    sqrtN->GetXaxis()->SetLimits(1,5000);
    sqrtN->GetHistogram()->SetMaximum(200);   // along          
    sqrtN->GetHistogram()->SetMinimum(0.1);
    cv2->SetLogy();
    cv2->SetLogx();

    // Output file name to save avg. spectrum (and anything else)
    TString ofname = infile.ReplaceAll("/",""); 
    ofname = ofname.ReplaceAll(".root",""); // If you send in folder names, expect a very long name
    TFile *outfile = new TFile(ofname+Form("_spectrum_%i.root",ndcms),"RECREATE");
    outfile->cd();
    spec_hist->Write();
    image_avg->Write();
    outfile->Close();

}
