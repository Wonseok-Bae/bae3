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

void Survival_Spectrum()
{
      double delta_31=0;
      double dm=0.9;
      double ds=6.;
      double L=20;
      double E[300000];
      double survival_prob=0;
      double number_of_neutrino[300000];
      double total_number_of_neutrino=0;
      double weight[300000];

      TFile * inputfile = new TFile("./outputFigs.root");

      TH1D* NEOS = (TH1D*)inputfile->Get("outPrediction[3]");
      NEOS->Draw();
      NEOS->Fit("pol9");
      TF1* func = NEOS->GetFunction("pol9"); 

      TFile* outputFile = new TFile("Survival_Spectrum.root","RECREATE");
     
      //====================================================================================================================================================================================

      //for 300000bin 
      TH1D * oscillated_spectrum_300000 = new TH1D("oscillated_spectrum_300000","title",300000,0.5,8);
      TH1D * unoscillated_spectrum_300000 = new TH1D("unoscillated_spectrum_300000","title",300000,0.5,8);
      TH1D * prob_300000 = new TH1D("prob_300000","title",300000,0.5,8);

      for(int i=1; i<300001; i++){
         //cout<<"i :"<<i<<" "<<0.5-0.025+i*0.025<<"~"<<0.5+i*0.025<<"\t"<<func->Eval(0.5-(0.025/2)+i*0.025)<<endl;
         E[i]=0.5-(0.000025/2)+i*0.000025;
         number_of_neutrino[i]=func->Eval(0.5-(0.000025/2)+i*0.000025);
         //cout<<"i :"<<i<<"E[i] : "<<0.5-(0.025/2)+i*0.025<<" number_of_neutrino[i] : "<<func->Eval(0.5-(0.025/2)+i*0.025)<<endl;       
         //cout<<i<<" "<<0.5-(0.025/2)+i*0.025<<" "<<func->Eval(0.5-(0.025/2)+i*0.025)<<endl;
      }

      TH1D* NEOS300000 = (TH1D*)outputFile->Get("oscillated_spectrum_300000");
      TH1D* NEOS300000_unosc = (TH1D*)outputFile->Get("unoscillated_spectrum_300000");
      TH1D* NEOS300000_prob = (TH1D*)outputFile->Get("prob_300000");

      for(Int_t i=1;i<300001;i++)
      {
              delta_31 = dm * (L/E[i])*1.27;
              survival_prob = 1 - ds *  TMath::Power(TMath::Sin(delta_31),2);
              
	      //cout<<"i: "<<i<<" E[i] : "<<E[i]<<" L : "<<L<<" dm : "<<dm<<" ds : "<<ds<< " survival_prob : "<<survival_prob<<" value:" << survival_prob*number_of_neutrino[i] <<endl;
              NEOS300000->SetBinContent(i, survival_prob*number_of_neutrino[i]);
	      NEOS300000_unosc->SetBinContent(i, number_of_neutrino[i] );
	      NEOS300000_prob->SetBinContent(i, survival_prob );



      }

      //=====================================================================================================================================================================================

      //for 300bin 
      TH1D * oscillated_spectrum_300 = new TH1D("oscillated_spectrum_300","title",300,0.5,8);
      TH1D * prob_300 = new TH1D("prob_300","title",250,1.75,8);

      int j=0;
      for(int i=1; i<251; i++){
	 //cout<<"i :"<<i<<" "<<0.5-0.025+i*0.025<<"~"<<0.5+i*0.025<<"\t"<<func->Eval(0.5-(0.025/2)+i*0.025)<<endl;
	 E[i]=1.75-(0.025/2)+i*0.025;
	 number_of_neutrino[i]=func->Eval(1.75-(0.025/2)+i*0.025);
	 total_number_of_neutrino = total_number_of_neutrino + number_of_neutrino[i];
         cout<<"i : "<<i<<", j : "<<j<<", E[i] : "<<E[i]<<", number_of_neutrino[i] : "<<number_of_neutrino[i]<<", total_number_of_neutrino : "<<total_number_of_neutrino<<endl;
         j++;
	 if(j==10){
		 
		 for(int k=9; k>=0; k--){
		 weight[i-k] = number_of_neutrino[i-k]*10./total_number_of_neutrino;
		 //cout<<"i-k : "<<i-k<<", weight[i-k] : "<<weight[i-k]<<endl;
		 }
		 j=0;
		 total_number_of_neutrino=0;
	 }	 
       }


       for(int i=1; i<251; i++){
       cout<<"i : "<<i<<", E[i] : "<<E[i]<<", weight[i] : "<<weight[i]<<endl;
       }


      //Maybe Trash 
      /*  
      for(int i=1; i<301; i++){

         double weight = number_of_neutrino[i]*10/total_number_of_neutrino;
	 cout<<"i : "<<i<<" weight : "<<weight<<endl;

      }
      */

       
      TH1D* NEOS300 = (TH1D*)outputFile->Get("oscillated_spectrum_300");
      TH1D* NEOS300_prob = (TH1D*)outputFile->Get("prob_300");

      for(Int_t i=1;i<251;i++)
      {
	      delta_31 = dm * (L/E[i])*1.27;
	      survival_prob = 1 - ds *  TMath::Power(TMath::Sin(delta_31),2);
              	      
	      //cout<<"i: "<<i<<" E[i] : "<<E[i]<<" L : "<<L<<" dm : "<<dm<<" ds : "<<ds<< " survival_prob : "<<survival_prob<<" value:" << survival_prob*number_of_neutrino[i] <<endl;
	      NEOS300->SetBinContent(i, survival_prob * number_of_neutrino[i] );
	      NEOS300_prob->SetBinContent(i, survival_prob );
     
      }
      
      //======================================================================================================================================================================================

      //for 30bin
      TH1D * oscillated_spectrum_30 = new TH1D("oscillated_spectrum_30","title",30,0.5,8);
      TH1D * unoscillated_spectrum_30 = new TH1D("unoscillated_spectrum_30","title",30,0.5,8);
      TH1D * prob_30 = new TH1D("prob_30","title",30,0.5,8);	

      for(int i=1; i<31; i++){
         //cout<<"i :"<<i<<" "<<0.5-0.25+i*0.25<<"~"<<0.5+i*0.25<<"\t"<<func->Eval(0.5-(0.25/2)+i*0.25)<<endl;
         E[i]=0.5-(0.25/2)+i*0.25;
         number_of_neutrino[i]=func->Eval(0.5-(0.25/2)+i*0.25);
         //cout<<"i :"<<i<<"E[i] : "<<0.5-(0.25/2)+i*0.25<<" number_of_neutrino[i] : "<<func->Eval(0.5-(0.25/2)+i*0.25)<<endl;
      }


      TH1D* NEOS30 = (TH1D*)outputFile->Get("oscillated_spectrum_30");
      TH1D* NEOS30_unosc = (TH1D*)outputFile->Get("unoscillated_spectrum_30");
      TH1D* NEOS30_prob = (TH1D*)outputFile->Get("prob_30");

      for(Int_t i=1;i<31;i++)
      {
              delta_31 = dm * (L/E[i])*1.27;
              survival_prob = 1 - ds *  TMath::Power(TMath::Sin(delta_31),2);
              	      
              //cout<<"i: "<<i<<" E[i] : "<<E[i]<<" L : "<<L<<" dm : "<<dm<<" ds : "<<ds<< " survival_prob : "<<survival_prob<<" value:" << survival_prob*number_of_neutrino[i] <<endl;
	      //cout<<E[i]<<" "<<survival_prob<<" "<<survival_prob*number_of_neutrino[i]<<" "<<number_of_neutrino[i]<<endl;
	      NEOS30->SetBinContent(i, survival_prob * number_of_neutrino[i] );
	      NEOS30_unosc->SetBinContent(i, number_of_neutrino[i] );
	      NEOS30_prob->SetBinContent(i, survival_prob );
      }

      oscillated_spectrum_30->Write();
      unoscillated_spectrum_30->Write();
      prob_30->Write();
      
      oscillated_spectrum_300->Write();
      prob_300->Write();

      oscillated_spectrum_300000->Write();
      unoscillated_spectrum_300000->Write();
      prob_300000->Write();
      
      outputFile->Close();
}

