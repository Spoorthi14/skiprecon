#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TParameter.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TFITS.h"
#include "TVector.h"
#include "TMath.h"
#include <iostream>
#include <stdio.h>

// Takes in a fits file and spits out a root file 
// run in Root interpreter e.g. skipRead("Imagefile1631.root", "1631", 0)
void skipRead(TString fits_fname, TString out_name, Int_t data_extension){
    
    //*********************
    Int_t offset = 2; // HARDCODED UNTIL SKIPPER FIX
    //*********************

    //file to output data
    out_name += ".root";
    
    TFile* f = new TFile(out_name, "RECREATE");
    
    //object from .fits
	TFITSHDU *hdu = new TFITSHDU(fits_fname, data_extension);

    //Now read the list with the header information
    ifstream hinfo("header.txt");
    
    TTree* finfo = new TTree("finfo","finfo");
    TList* vars = new TList();
    
    //First add the most general information
    FILE* sysout;
    char ct[80];
    int lastchar;
    sysout = popen("pwd","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString pwd(ct);
    pclose(sysout);
    sysout = popen("hostname","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString hostname(ct);
    pclose(sysout);
    sysout = popen("svnversion","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString svnversion(ct);
    pclose(sysout);
    
    vars->Add(new TNamed("HOSTNAME",hostname.Data()));
    vars->Add(new TNamed("PWD",pwd.Data()));
    vars->Add(new TNamed("SVNVERSION",svnversion.Data()));
    
    //Now the information from the .fits header
    TString vname, vtype, vval;
    while(hinfo>>vname && hinfo>>vtype){
        
        vval = hdu->GetKeywordValue(vname);
        
        if(vval.Length()==0)
            std::cout << "Warning: " << vname << " does not exist in .fits header!" << std::endl;
        
        if(vtype=="D" || vtype=="I")
            vars->Add(new TParameter<double>(vname,vval.Atof()));
        
        else if(vtype=="S")
            vars->Add(new TNamed(vname,vval));
    }
    hinfo.close();
    
    //Header info not from .fits file
    vars->AddFirst(new TNamed(TString("ROOTFPATH"), out_name));
    vars->AddFirst(new TNamed(TString("FITSFPATH"), fits_fname));
    
    finfo->Branch("EXTID", &data_extension, "EXTID/I");
    finfo->Branch(vars);
    finfo->Fill();
    
    Int_t nbins_x_skip = (Int_t) ((TParameter<double>*) vars->FindObject("NAXIS1"))->GetVal();
    Int_t nbins_y = (Int_t) ((TParameter<double>*) vars->FindObject("NAXIS2"))->GetVal();
    Int_t ndcms = (Int_t) ((TParameter<double>*) vars->FindObject("NDCMS"))->GetVal();
    Int_t nbins_x_actual = nbins_x_skip/ndcms;
    
    //histogram to keep raw image
    TH2F* image_raw = new TH2F("image_raw", "image_raw", nbins_x_skip, 0, nbins_x_skip, nbins_y, 0, nbins_y);
    TH3F* image_cube = new TH3F("image_cube","image_cube", nbins_x_actual,0,nbins_x_actual, nbins_y, 0, nbins_y, ndcms, 0, ndcms);
    
    //fill histogram
    Int_t actual = 0;
    for(Int_t i=1; i<=nbins_x_skip; i++){
                
        TVectorD* currcol = hdu->GetArrayColumn(i-1);
        
            Int_t slice = (i-offset)%ndcms;
            if (slice==1) actual+=1;

            for(Int_t j=1; j<=nbins_y; j++){
         
                Int_t cont = (Int_t) (*currcol)[j-1];
                image_raw->SetBinContent(i,j,cont);
                if (i>offset){
                    image_cube->SetBinContent(actual,j,slice,cont);
                    
                }
            }
        
        currcol->Delete();
    }
    
    //delete the hdu object to clear memory
    hdu->Delete();
    f->cd();

    //write the image
    finfo->Write();
    image_raw->Write();
    image_cube->Write();
    
    f->Close();

}
