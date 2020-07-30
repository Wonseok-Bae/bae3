#include "TMath.h"

#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TString.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>

using namespace std;
double sum=0;
double random_value=0;
void getBinContent_75bin()
{
       
      TFile * inputfile = new TFile("./straight2_outputFigs_75bin.root");

      TH1D* NEOS = (TH1D*)inputfile->Get("outPrediction[3]");
      NEOS->Draw();
      NEOS->Fit("pol9");
      TF1* func = NEOS->GetFunction("pol9"); 
      
      for(int i=0; i<75; i++){
	 cout<<+i*0.1<<"~"<<(i+1)*0.1<<"\t"<<func->Eval(+0.1/2+i*0.1)<<endl; 
	 sum = sum + func->Eval(+0.1/2+i*0.1);
	 cout<<" sum : "<<sum<<endl;
      }

      TFile* outputFile = new TFile("Escalefraction_75bin.root","RECREATE");
      TH1D * outPrediction = new TH1D("outPrediction[3]","title",75,0.5,8);
      TH1D * outPrediction_shifted = new TH1D("outPrediction_shifted[3]","title",75,0.5,8);	
     
      for(int i=1; i<=36251000; i++){//xxx
	 gRandom->SetSeed(i*1234);
         //cout<< func->GetRandom() <<endl;
	 random_value = func->GetRandom();
         outPrediction->Fill(random_value,1);
	 outPrediction_shifted->Fill(1.048*random_value,1);
      }
            
      TH1D* NEOS1 = (TH1D*)outputFile->Get("outPrediction[3]");
      TH1D* NEOS2 = (TH1D*)outputFile->Get("outPrediction_shifted[3]");
         
      TH1D * outPrediction_subtraction = new TH1D("Escalefraction","Escalefraction",75,0.5,8);

      for(Int_t i=0;i<65;i++)
      {
              //cout<<" i: "<<i<<" value:" << NEOS1->GetBinContent(i) - NEOS2->GetBinContent(i) <<endl;
              outPrediction_subtraction->SetBinContent(i+1,  (NEOS1->GetBinContent(i+1) - NEOS2->GetBinContent(i+1))/NEOS1->GetBinContent(i+1)  );
      }

      for(Int_t i=65;i<75;i++)
      {
             outPrediction_subtraction->SetBinContent(i+1,0);
      }
      
      outPrediction->Write();
      outPrediction_shifted->Write();  
      outPrediction_subtraction->Write();
      outputFile->Close();
}

