#include <iostream>
#include <TH1D.h>
#include <TVirtualFitter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include "TLatex.h"

void add6cent(){

  using namespace std;
  
  const char *pairs[4] = { "KpPro", "KmAPro", "KpAPro", "KmPro" };
  const char *pairs2[4] = { "K^{+}p", "K^{-}#bar{p}", "K^{+}#bar{p}", "K^{-}p" };




  TH1D *hpart[4];


 // TFile *f1=new TFile("corrs.alice.5TeV.1904.plus.root");




  TFile *f2=new TFile(Form("plus.%s.m6.root",pairs[0]));
  hpart[0] = (TH1D*)f2->Get(Form("CfnReYlm00cylm%stpcM6",pairs[0]));

  TFile *f3=new TFile(Form("plus.%s.m6.root",pairs[1]));
  hpart[1] = (TH1D*)f3->Get(Form("CfnReYlm00cylm%stpcM6",pairs[1]));

  TFile *f4=new TFile(Form("plus.%s.m6.root",pairs[2]));
  hpart[2] = (TH1D*)f4->Get(Form("CfnReYlm00cylm%stpcM6",pairs[2]));

  TFile *f5=new TFile(Form("plus.%s.m6.root",pairs[3]));
  hpart[3] = (TH1D*)f5->Get(Form("CfnReYlm00cylm%stpcM6",pairs[3]));
  //hpart[i]->SetName("low3");
  TFile *f = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root","UPDATE");
  hpart[0]->Write();  
  hpart[1]->Write();  
  hpart[2]->Write();  
  hpart[3]->Write();  


}


