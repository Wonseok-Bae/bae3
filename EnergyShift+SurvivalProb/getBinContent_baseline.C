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
void getBinContent_baseline()
{
      TFile * inputfile1 = new TFile("./straight_outputFigs2_23.7.root");
      TH1D* NEOS1 = (TH1D*)inputfile1->Get(Form("scaledOscPred_%s_%s", gApplication->Argv(4), gApplication->Argv(5) ));
      /*
      for(int i=0; i<75; i++){
	 cout<<"i+1 : "<<i+1<<" NEOS1->GetBinContent(i) : "<<NEOS1->GetBinContent(i+1)<<endl;
      }
      */
      TFile * inputfile2 = new TFile("./straight_outputFigs2_26.7.root");
      TH1D* NEOS2 = (TH1D*)inputfile2->Get(Form("scaledOscPred_%s_%s", gApplication->Argv(4), gApplication->Argv(5) ));
      /*
      for(int i=0; i<75; i++){
         cout<<"i+1 : "<<i+1<<" NEOS2->GetBinContent(i) : "<<NEOS2->GetBinContent(i+1)<<endl;
      }
      */
      
      TFile* outputFile = new TFile("Baselinefraction.root","UPDATE");     
      TH1D * outPrediction_subtraction = new TH1D("Baselinefraction","",75,0.5,8);

      for(Int_t i=0;i<75;i++)
      {
              //cout<<" i: "<<i<<" value: " << NEOS1->GetBinContent(i+1) << " - " << NEOS2->GetBinContent(i+1) <<endl;
              outPrediction_subtraction->SetBinContent(i+1,  (NEOS1->GetBinContent(i+1) - NEOS2->GetBinContent(i+1))/NEOS1->GetBinContent(i+1)  );
      }
      /*
      for(Int_t i=65;i<75;i++)
      {
             outPrediction_subtraction->SetBinContent(i+1,0);
      }
      */
      //outPrediction_subtraction->Write();
      outPrediction_subtraction->Write(Form("Baselinefraction_%s_%s",  gApplication->Argv(4), gApplication->Argv(5) ));
}

