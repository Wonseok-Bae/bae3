//May,20,2020 : Scaling for Fitted plot
//May,23,2020 : Input Order ds,dm,experiment / Output Order dm,ds,experiment

#include "simple_t2k_wonseok.hh"
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

 int main(int argc, char**argv){


 RooFitResult* res;
 Sterile * rep = new Sterile ("_rep");
   char formula[10];

 std::cout<<"start to run "<<std::endl;

 rep->setAtmBaseline(0);
 rep->setDensity(0);
 TH1D* vecInput1 = new TH1D("","",11,0,11);
 TH1D* vecInput2 = new TH1D("","",11,0,11);
 //vecInput1->Print();

 vecInput1 ->SetBinContent(1, TMath::ASin(TMath::Sqrt(0.85))/2.);
 vecInput1 ->SetBinContent(2, 0); //TMath::ASin(TMath::Sqrt(0.95))/2.);
 vecInput1 ->SetBinContent(3, 0); //0.1);
 vecInput1 ->SetBinContent(4,-1.5);
 vecInput1 ->SetBinContent(5,0.000075);
 vecInput1 ->SetBinContent(6, 0 ); //0.00238);
 vecInput1 ->SetBinContent(7, 0 ); //0.00244);
 vecInput1 ->SetBinContent(8,1);
 vecInput1 ->SetBinContent(9,1);
 vecInput1 ->SetBinContent(10,1);
 //vecInput1 ->SetBinContent(11,1);//wrong
 vecInput1 ->SetBinContent(11,0);//right

 vecInput2 ->SetBinContent(1,TMath::ASin(TMath::Sqrt(0.85+0.021))/2.-TMath::ASin(TMath::Sqrt(0.85))/2.);
 vecInput2 ->SetBinContent(2, 100); //2);
 vecInput2 ->SetBinContent(3,100);
 vecInput2 ->SetBinContent(4,2);
 vecInput2 ->SetBinContent(5,0.0000018);
 vecInput2 ->SetBinContent(6, 100); //0.0001);
 vecInput2 ->SetBinContent(7, 100); //0.00009);
 //vecInput2 ->SetBinContent(8,0.009);
 //vecInput2 ->SetBinContent(9,0.03);
 vecInput2 ->SetBinContent(8, 100);
 vecInput2 ->SetBinContent(9, 100);
 vecInput2 ->SetBinContent(10,0.04); 
 vecInput2 ->SetBinContent(11,1);//right
 // vecInput2 ->SetBinContent(11,0.027);//wrong

 std::cout<<"'ve set some inputs "<<std::endl;

 // binned from 0.5-9 MeV with bin width of 0.25 MeV
 /* 
 int nBins = 30;
 TH1D* binHist = new TH1D("","",nBins+1,0,nBins+1);
 for(Int_t i=0;i<nBins+1;i++){
   binHist->SetBinContent(i+1, 0.5 + 0.25*i);
 }
 */
 
 int nBins = 75;
 TH1D* binHist = new TH1D("","",nBins+1,0,nBins+1);
 for(Int_t i=0;i<nBins+1;i++){
   binHist->SetBinContent(i+1, 0.5 + 0.1*i);
 }

 //TString fileLocation = "/gpfs/projects/McGrewGroup/gyang/REACTOR-related/";
 TString fileLocation = "./";
 rep->setFileLocation(fileLocation);

 TH1D* fissionHist = new TH1D("","",16,0,16);
 double fissionList[16]={0.49,0.09,0.35,0.07,    0.501,0.073,0.353,0.072,    0.679,0.072,0.215,0.035,      0.679,0.072,0.215,0.035};
 for(Int_t i=0;i<16;i++){
   fissionHist->SetBinContent(i+1, fissionList[i]);
 }

 rep->SetMatrixNameDC(fileLocation+"data/ReactorMatrixND_shape_rho0_NU.root");
 rep->SetMatrixNameDYB(fileLocation+"data/NEOSCorrelationMatrix.root");
 rep->SetMatrixNameRENO(fileLocation+"data/NEOSCorrelationMatrix.root");
 rep->SetMatrixNameNEOS(fileLocation+"data/NEOSCorrelationMatrix.root");
 rep->SetMatrixNamePROS(fileLocation+"data/NEOSCorrelationMatrix.root");

 rep->SetBinning(binHist);
 rep->SetFissionFraction(fissionHist);

 std::cout<<" setting the conversion matrix . "<<std::endl;
 std::vector<TMatrixD*> cfMatrix = rep->ConversionMatrix(fileLocation+"data/FinalFitIBDTree_DCIV_ND.root","FinalFitIBDTree");
 //rep->ConversionMatrix("data/FinalFitIBDTree_DCIV_ND.root","FinalFitIBDTree");
 //TH2D* fHist = rep->GetConversionMatrix();

 //std::cout<<"'ve set some more inputs, conversiont matrix summation is : "<<fHist->Integral()<<std::endl;

 std::vector<TString> mList;
 TString r235 = (fileLocation+"data/huber235.txt");
 TString r238 = (fileLocation+"data/mueller238.txt");
 TString r239 = (fileLocation+"data/huber239.txt");
 TString r241 = (fileLocation+"data/huber241.txt");

 mList.push_back(r235);
 mList.push_back(r238); 
 mList.push_back(r239);
 mList.push_back(r241);

 rep->SetModelList(mList);
 //rep->preparePrediction(rep->getPullList(), false);
 rep->ifEShiftHist(false);

 std::vector<TH1D*> tempPredList = rep->preparePrediction(rep->getPullList(), false);
 rep->prepareData(rep->preparePrediction(rep->getPullList(), false));

 rep->setPull(vecInput1); 
 rep->setPullUnc(vecInput2);
 rep->addSK(true);

 rep->setNBins(34);
 rep->setTime(6);

 std::cout<<"ended up with setting us basic stuff "<<std::endl;

 rep->FillEv(rep->getPullList());

 std::cout<<"'ve set some more & more inputs "<<std::endl;

 RooArgList list("list");
 list.add(*rep);
 //std::cout<<"formula[0] : "<< formula[0] <<" formula[1] : "<< formula[1] <<" formula[2] : "<< formula[2] << " formula[3] : "<<formula[3] <<" formula[4] : "<< formula[4] << " formula[5] : "<< formula[5] << " formula[6] : "<< formula[6] << " formula[7] : "<< formula[7] << " formula[8] : "<< formula[8] << " formula[9] : "<<formula[9] <<" formula[10] : "<< formula[10] << std::endl;
 sprintf(formula,"%s","@0");
 std::cout<<"sprintf(formula, s, @0); "<< formula <<std::endl;
 RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

 // ******************************** Important setup here *************************************
 // *******************************************************************************************
 rep->fitSingleExp(Form("%s",argv[3]));
 rep->ifEqualIso(false);
 rep->setSysts(true);

 rep->ifEShiftHist(true);
 rep->setEShiftHist("Escalefraction.root");
 // *******************************************************************************************
 // ******************************************************************************************* 

 rep->setBaselineDYB(560.0);
 rep->setBaselineDC(400.0);
 rep->setBaselineNEOS(20.0);
 rep->setBaselineRENO(410.6);
 rep->setBaselinePROS(7.9);

 //rep->getParVar(3)->setVal(2.5);
 rep->getParVar(2)->setConstant(true);
 rep->getParVar(6)->setConstant(true);

 rep->getParVar(0)->setConstant(true);
 rep->getParVar(1)->setConstant(true);
 rep->getParVar(3)->setConstant(true);
 rep->getParVar(4)->setConstant(true);
 rep->getParVar(5)->setConstant(true);
 rep->getParVar(7)->setConstant(true);
 rep->getParVar(8)->setConstant(true);
 rep->getParVar(9)->setConstant(true);
 //rep->getParVar(10)->setConstant(true);
 if (!rep->getIfEShiftHist())
   rep->getParVar(10)->setConstant(true);
 else
   rep->getParVar(10)->setConstant(false);

/*
 rep->getParVar(12)->setConstant(true);  // 0.5 - 0.75
 rep->getParVar(13)->setConstant(true);  // 0.75 - 1
 rep->getParVar(14)->setConstant(true);  // 1 - 1.25
 rep->getParVar(15)->setConstant(true);  // 1.25 - 1.5
 rep->getParVar(16)->setConstant(true);  // 1.5 - 1.75
 rep->getParVar(17)->setConstant(true);  // 1.75 - 2
 rep->getParVar(18)->setConstant(true);  // 2 - 2.25
 rep->getParVar(19)->setConstant(true);  // 2.25 - 2.5
 rep->getParVar(20)->setConstant(true);  // 2.5 - 2.75
 rep->getParVar(21)->setConstant(true);  // 2.75 - 3
 rep->getParVar(22)->setConstant(true);  // 3 - 3.25
 rep->getParVar(23)->setConstant(true);  // 3.25 - 3.5
 rep->getParVar(24)->setConstant(true);  // 3.5 - 3.75
 rep->getParVar(25)->setConstant(true);  // 3.75 - 4
 rep->getParVar(26)->setConstant(true);  // 4 - 4.25
 rep->getParVar(27)->setConstant(true);  // 4.25 - 4.5
 rep->getParVar(28)->setConstant(true);  // 4.5 - 4.75
 rep->getParVar(29)->setConstant(true);  // 4.75 - 5
 rep->getParVar(30)->setConstant(true);  // 5 - 5.25
 rep->getParVar(31)->setConstant(true);  // 5.25 - 5.5
 rep->getParVar(32)->setConstant(true);  // 5.5 - 5.75
 rep->getParVar(33)->setConstant(true);  // 5.75 - 6
 rep->getParVar(34)->setConstant(true);  // 6 - 6.25
 rep->getParVar(35)->setConstant(true);  // 6.25 - 6.5
 rep->getParVar(36)->setConstant(true);  // 6.5 - 6.75
 rep->getParVar(37)->setConstant(true);  // 6.75 - 7
 rep->getParVar(38)->setConstant(true);  // 7 - 7.25
 rep->getParVar(39)->setConstant(true);  // 7.25 - 7.5
 rep->getParVar(40)->setConstant(true);  // 7.5 - 7.75
 rep->getParVar(41)->setConstant(true);  // 7.75 - 8
 rep->getParVar(42)->setConstant(true);  // 8 - 8.25
 rep->getParVar(43)->setConstant(false);  // 8.25 - 8.5
 rep->getParVar(44)->setConstant(true);  // 8.5 - 8.75
 rep->getParVar(45)->setConstant(true);  // 8.75 - 9
*/
 rep->getParVar(11)->setConstant(true);
 rep->getParVar(12)->setConstant(true);  // 0.5 - 0.75
 rep->getParVar(13)->setConstant(true);  // 0.75 - 1
 rep->getParVar(14)->setConstant(true);  // 1 - 1.25
 rep->getParVar(15)->setConstant(true);  // 1.25 - 1.5
 rep->getParVar(16)->setConstant(true);  // 1.5 - 1.75
 rep->getParVar(17)->setConstant(true);  // 1.75 - 2
 rep->getParVar(18)->setConstant(true);  // 2 - 2.25
 rep->getParVar(19)->setConstant(true);  // 2.25 - 2.5
 rep->getParVar(20)->setConstant(true);  // 2.5 - 2.75
 rep->getParVar(21)->setConstant(true);  // 2.75 - 3
 rep->getParVar(22)->setConstant(true);  // 3 - 3.25
 rep->getParVar(23)->setConstant(true);  // 3.25 - 3.5
 rep->getParVar(24)->setConstant(true);  // 3.5 - 3.75
 rep->getParVar(25)->setConstant(true);  // 3.75 - 4
 rep->getParVar(26)->setConstant(true);  // 4 - 4.25
 rep->getParVar(27)->setConstant(true);  // 4.25 - 4.5
 rep->getParVar(28)->setConstant(true);  // 4.5 - 4.75
 rep->getParVar(29)->setConstant(true);  // 4.75 - 5
 rep->getParVar(30)->setConstant(true);  // 5 - 5.25
 rep->getParVar(31)->setConstant(true);  // 5.25 - 5.5
 rep->getParVar(32)->setConstant(true);  // 5.5 - 5.75
 rep->getParVar(33)->setConstant(true);  // 5.75 - 6
 rep->getParVar(34)->setConstant(true);  // 6 - 6.25
 rep->getParVar(35)->setConstant(true);  // 6.25 - 6.5
 rep->getParVar(36)->setConstant(true);  // 6.5 - 6.75
 rep->getParVar(37)->setConstant(true);  // 6.75 - 7
 rep->getParVar(38)->setConstant(true);  // 7 - 7.25
 rep->getParVar(39)->setConstant(true);  // 7.25 - 7.5
 rep->getParVar(40)->setConstant(true);  // 7.5 - 7.75
 rep->getParVar(41)->setConstant(true);  // 7.75 - 8
 rep->getParVar(42)->setConstant(true);  // 8 - 8.25
 rep->getParVar(43)->setConstant(true);  // 8.25 - 8.5
 rep->getParVar(44)->setConstant(true);  // 8.5 - 8.75
 rep->getParVar(45)->setConstant(true);  // 8.75 - 9

/*
 rep->getParVar(12)->setConstant(true);  // 0.5 - 0.75
 rep->getParVar(13)->setConstant(true);  // 0.75 - 1
 rep->getParVar(14)->setConstant(true);  // 1 - 1.25
 rep->getParVar(15)->setConstant(true);  // 1.25 - 1.5
 rep->getParVar(16)->setConstant(true);  // 1.5 - 1.75
 rep->getParVar(17)->setConstant(true);  // 1.75 - 2
 rep->getParVar(18)->setConstant(true);  // 2 - 2.25
 rep->getParVar(19)->setConstant(true);  // 2.25 - 2.5
 rep->getParVar(20)->setConstant(true);  // 2.5 - 2.75
 rep->getParVar(20)->setConstant(false);  // 2.5 - 2.75
 rep->getParVar(21)->setConstant(true);  // 2.75 - 3
 rep->getParVar(22)->setConstant(true);  // 3 - 3.25
 rep->getParVar(23)->setConstant(true);  // 3.25 - 3.5
 rep->getParVar(24)->setConstant(true);  // 3.5 - 3.75
 rep->getParVar(25)->setConstant(true);  // 3.75 - 4
 rep->getParVar(26)->setConstant(true);  // 4 - 4.25
 rep->getParVar(27)->setConstant(true);  // 4.25 - 4.5
 rep->getParVar(28)->setConstant(true);  // 4.5 - 4.75
 rep->getParVar(29)->setConstant(true);  // 4.75 - 5
 rep->getParVar(30)->setConstant(true);  // 5 - 5.25
 rep->getParVar(31)->setConstant(true);  // 5.25 - 5.5
 rep->getParVar(32)->setConstant(true);  // 5.5 - 5.75
 rep->getParVar(33)->setConstant(true);  // 5.75 - 6
 rep->getParVar(34)->setConstant(true);  // 6 - 6.25
 rep->getParVar(35)->setConstant(true);  // 6.25 - 6.5
 rep->getParVar(36)->setConstant(true);  // 6.5 - 6.75
 rep->getParVar(37)->setConstant(true);  // 6.75 - 7
 rep->getParVar(38)->setConstant(true);  // 7 - 7.25
 rep->getParVar(39)->setConstant(true);  // 7.25 - 7.5
 rep->getParVar(40)->setConstant(true);  // 7.5 - 7.75
 rep->getParVar(41)->setConstant(true);  // 7.75 - 8
 rep->getParVar(42)->setConstant(true);  // 8 - 8.25
 rep->getParVar(43)->setConstant(true);  // 8.25 - 8.5
 rep->getParVar(44)->setConstant(true);  // 8.5 - 8.75
 rep->getParVar(45)->setConstant(true);  // 8.75 - 9
*/

 // energy scales for four exp.
 rep->getParVar(45)->setConstant(true);  
 rep->getParVar(46)->setConstant(true);  
 rep->getParVar(47)->setConstant(true);  
 rep->getParVar(48)->setConstant(true);  
 rep->getParVar(49)->setConstant(true);  
 rep->getParVar(50)->setConstant(true);  
 rep->getParVar(51)->setConstant(true);  
 rep->getParVar(52)->setConstant(true);
 rep->getParVar(53)->setConstant(true);
 rep->getParVar(54)->setConstant(true);
 rep->getParVar(55)->setConstant(true);  // 1.25 - 1.5
 rep->getParVar(56)->setConstant(true);  // 1.5 - 1.75
 rep->getParVar(57)->setConstant(true);  // 1.75 - 2
 rep->getParVar(58)->setConstant(true);  // 2 - 2.25
 rep->getParVar(59)->setConstant(true);  // 2.25 - 2.5
 rep->getParVar(60)->setConstant(true);  // 2.5 - 2.75
 rep->getParVar(61)->setConstant(true);  // 2.75 - 3
 rep->getParVar(62)->setConstant(true);  // 3 - 3.25
 rep->getParVar(63)->setConstant(true);  // 3.25 - 3.5
 rep->getParVar(64)->setConstant(true);  // 3.5 - 3.75
 rep->getParVar(65)->setConstant(true);  // 3.75 - 4
 rep->getParVar(66)->setConstant(true);  // 4 - 4.25
 rep->getParVar(67)->setConstant(true);  // 4.25 - 4.5
 rep->getParVar(68)->setConstant(true);  // 4.5 - 4.75
 rep->getParVar(69)->setConstant(true);  // 4.75 - 5
 rep->getParVar(70)->setConstant(true);  // 5 - 5.25
 rep->getParVar(71)->setConstant(true);  // 5.25 - 5.5
 rep->getParVar(72)->setConstant(true);  // 5.5 - 5.75
 rep->getParVar(73)->setConstant(true);  // 5.75 - 6
 rep->getParVar(74)->setConstant(true);  // 6 - 6.25
 rep->getParVar(75)->setConstant(true);  // 6.25 - 6.5
 rep->getParVar(76)->setConstant(true);  // 6.5 - 6.75
 rep->getParVar(77)->setConstant(true);  // 6.75 - 7
 rep->getParVar(78)->setConstant(true);  // 7 - 7.25
 rep->getParVar(79)->setConstant(true);  // 7.25 - 7.5
 rep->getParVar(80)->setConstant(true);  // 7.5 - 7.75
 rep->getParVar(81)->setConstant(true);  // 7.75 - 8
 rep->getParVar(82)->setConstant(true);  // 8 - 8.25
 rep->getParVar(83)->setConstant(true);  // 8.25 - 8.5
 rep->getParVar(84)->setConstant(true);  // 8.5 - 8.75
 rep->getParVar(85)->setConstant(true);  // 1.25 - 1.5
 rep->getParVar(86)->setConstant(true);  // 1.5 - 1.75
 rep->getParVar(87)->setConstant(true);  // 1.75 - 2
 rep->getParVar(88)->setConstant(true);  // 2 - 2.25
 rep->getParVar(89)->setConstant(true);  // 2.25 - 2.5
 rep->getParVar(90)->setConstant(true);  // 2.5 - 2.75
 rep->getParVar(91)->setConstant(true);  // 2.75 - 3
 rep->getParVar(92)->setConstant(true);  // 3 - 3.25
 rep->getParVar(93)->setConstant(true);  // 3.25 - 3.5
 rep->getParVar(94)->setConstant(true);  // 3.5 - 3.75
 rep->getParVar(95)->setConstant(true);  // 3.75 - 4
 rep->getParVar(96)->setConstant(true);  // 4 - 4.25
 //rep->getParVar(97)->setConstant(true);  // 4.25 - 4.5
 
 //    for(Int_t nn= 0; nn< npts ; ++nn) {

 //rep->setNBins(binSetup );
 //rep->setTime(atof(argv[5]));

 std::cout<<"------------  Getting current spectra "<<std::endl;
 std::vector<TH1D*> outPrediction = tempPredList; // rep->GetCurrentPrediction();
 std::vector<TH1D*> outData = rep->GetCurrentData(outPrediction);
 std::cout<<"------------  Have got current spectra "<<std::endl;

 // save the un-oscillated standard prediction and data spectra
 TFile* outputFile = new TFile("straight_outputFigs.root","RECREATE");
 for(Int_t i=0;i<outPrediction.size();i++)
 {
   outPrediction[i]->Write(Form("outPrediction[%d]",i));
   if(i<5) outData[i]->Write(Form("outData[%d]",i));
 }
 cfMatrix[0]->Write("cfMatrix[0]");
 outputFile -> Close();
 //exit(1);
 std::cout<<"first thing saved "<<std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // variations of the predicted spectra with different oscillation parameters
 if(false){
	 
   TFile* outputFileOsc = new TFile("outputFigVar1.root","RECREATE"); 
   float dmVar = 0.001 ;
   float stVar = 0.1 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar); 
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     cout<<"-----------dmVar-------------- "<<dmVar<<endl;
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
   }
   outputFileOsc -> Close();

   TFile* outputFileOsc2 = new TFile("outputFigVar2.root","RECREATE");
   dmVar = 0.1 ;
   stVar = 0.1 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar); 
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar, stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i, dmVar,stVar));
   }
   outputFileOsc2 -> Close();

   TFile* outputFileOsc3 = new TFile("outputFigVar3.root","RECREATE");
   dmVar = 1 ;
   stVar = 0.1 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar); 
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar ));
   }
   outputFileOsc3 -> Close();

   TFile* outputFileOsc4 = new TFile("outputFigVar4.root","RECREATE");
   dmVar = 0.001 ;
   stVar = 0.01 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar); 
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar ));
   }
   outputFileOsc4 -> Close();

   TFile* outputFileOsc5 = new TFile("outputFigVar5.root","RECREATE");
   dmVar = 0.01 ;
   stVar = 0.01 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar);
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
   }
   outputFileOsc5 -> Close();

   TFile* outputFileOsc6 = new TFile("outputFigVar6.root","RECREATE");
   dmVar = 1 ;
   stVar = 0.01 ;
   rep->getParVar(2)->setVal(stVar);
   rep->getParVar(6)->setVal(dmVar); 
   outPrediction = rep->GetCurrentPrediction();
   //outData = rep->GetCurrentData(outPrediction);
   for(Int_t i=0;i<outPrediction.size();i++)
   {
     outPrediction[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
     //outData[i]->Write(Form("outOsc[%d]_dm%f_st%f",i,dmVar,stVar));
   }
   outputFileOsc6 -> Close();
 }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 ofstream outText;
 outText.open(Form("./straight_scan2Dnew_%s_fissionFree_ISO%d_SYST%d_%s_%s.txt",argv[3], rep->GetEqualIso(), rep->GetSysts(), argv[1], argv[2]) );
 //outText.open(Form(fileLocation+"scan2Dnew_%s_fissionFixBump_ISO%d_SYST%d_%d_%d.txt",argv[3], rep->GetEqualIso(), rep->GetSysts(), atoi(argv[1]), atoi(argv[2]) ));
 //outText.open(Form("./scan2Dnew_%s_fissionFree_ISO%d_SYST%d_%s_%s_%s.txt",argv[4], rep->GetEqualIso(), rep->GetSysts(), argv[1], argv[2], argv[3]) );

 double iDM = atof(argv[2]);//May,23,2020//dm
 double iST = atof(argv[1]);//May,23,2020//ds
// double iSTfrom = atof(argv[2]);
// double iSTto = atof(argv[3]);
  
 //while( iSTfrom < iSTto )
 //{
 //for(Int_t iDM=6;iDM<21;iDM++)
 //{
 //    for(Int_t iST=0;iST<21;iST++)
 //    {
        // means that s2t14 0.001 - 1 and dm2 0.01 - 10
        rep->getParVar(2)->setVal(iST);//ds
        rep->getParVar(6)->setVal(iDM);//dm
	//rep->getParVar(2)->setVal(iST*0.02);
        //rep->getParVar(6)->setVal(iDM*0.4);
        //rep->getParVar(2)->setVal(iST);
	//rep->getParVar(6)->setVal(iDM);
	rep->getParVar(2)->setConstant(true);
 	rep->getParVar(6)->setConstant(true);
	rep->getParVar(7)->setConstant(true);
	rep->getParVar(8)->setConstant(true);

	//rep->prepareData(rep->preparePrediction(rep->getPullList(), false));
        //double scaling4 = rep->getScaling4();//June,3,2020 before the fitting
        //std::cout<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& scaling4 in app : "<< scaling4 << std::endl;
 	RooMinuit m(*fcn);
 	m.setStrategy(2);
 	Double_t callsEDM[2] = {10500., 1.e-6};
 	Int_t irf = 0;

 	gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);
 	m.migrad();
 	//m.hesse();
 	//m.minos(); 
 	res = m.save();
 	double bestFit = res->minNll();
        std::cout<<"fit status code is : "<< res->status()<<std::endl;	
	//status = 0    : OK
  	//status = 1    : Covariance was mad  epos defined
   	//status = 2    : Hesse is invalid
   	//status = 3    : Edm is above max
   	//status = 4    : Reached call limit
   	//status = 5    : Any other failure
	std::cout<<"quiality code of covariance matrix is : "<< res->covQual()<<std::endl;
	//status = -1 :  not available (inversion failed or Hesse failed)
        //status =  0 : available but not positive defined
        //status =  1 : covariance only approximate
        //status =  2 : full matrix but forced pos def
        //status =  3 : full accurate matrix
        double scaling4 = rep->getScaling4();//June,3,2020 after the fitting
	outText<<iDM<<" "<<iST<<" "<<bestFit<<" "<<res->status()<<" "<<res->covQual()<<" "<<scaling4<<std::endl;
//     }
// }
	//iSTfrom++;
 //}

 //////////////////////	
 //inline Int_t status() const 
 //inline Int_t covQual() const 
 //inline Int_t numInvalidNLL() const 
 //inline Double_t edm() const 
 //inline Double_t minNll() const 
 /////////////////////

 outPrediction = rep->GetCurrentPrediction();
 //std::cout<<"outPrediction.size() : "<<outPrediction.size()<<" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;

 //rep->GetWonseokCurrentPrediction();//May,20,2020

 //outData = rep->GetCurrentData(outPrediction);


 //std::cout<<"list of reactor pulls : "<<std::endl;
 //for(Int_t i=18;i<41;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}

 TFile* outputFile2 = new TFile("straight_outputFigs2.root","RECREATE");
 for(Int_t i=0;i<outPrediction.size();i++)
 //for(Int_t i=0;i<75;i++)
 {
   outPrediction[i]->Write(Form("outFit[%d]",i));	 
 }

 //outputFile2->Close();//May,20,2020

 std::cout<<"size of output prediction list "<<outPrediction.size()<<std::endl;
 std::cout<<"size of output data list "<<outData.size()<<std::endl;
 //std::cout<<atoi(argv[1])<<" "<<atoi(argv[2])<<" "<<bestFit<<std::endl;

 //std::cout<<"result list: "<<std::endl;
 std::cout<<"chi2: "<<bestFit <<std::endl;

 double bb =  rep->getParVar(2)->getAsymErrorLo();
 double dd =  rep->getParVar(2)->getAsymErrorHi();
 
 TVectorD* outVec = rep->getTestVec();//May,21,2020
 //std::cout<<(*outVec)[0]<<" "<<(*outVec)[29]<<std::endl;//May,21,2020
 
 
 TH1D* scaledOscPred = new TH1D("scaledOscPred","", 75 ,0.5, 8);
 for(Int_t i=0;i< scaledOscPred->GetNbinsX();i++)
 {
 	scaledOscPred->SetBinContent(i+1,(*outVec)[i]);
 }
 scaledOscPred->Write("scaledOscPred");
 outputFile2->Close();//May,20,2020
 
 //cout<<"errors are "<<bb<<" "<<dd<<endl;

 //for(Int_t i=0;i<11;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}
 //int aa = (atof(argv[5])+0.000001)*100;
 //cout<<atoi(argv[2])<<" "<<aa<<" "<<dd<<endl; 
 //  }
 }
