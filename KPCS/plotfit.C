void preparepad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.003);
  gPad->SetRightMargin(0.005);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.11);
  
}

void preparehist(TH1D *hist)
{
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.9);
  hist->SetLineColor(kRed);
  hist->SetMarkerColor(kRed);

  hist->GetXaxis()->SetRange(1,75);
  
  hist->SetTitle(";k* (GeV/#it{c});C(k*)");

  hist->GetXaxis()->SetNdivisions(610);
  
  hist->GetXaxis()->SetLabelSize(0.051);
  hist->GetYaxis()->SetLabelSize(0.051);

  hist->GetXaxis()->SetTitleSize(0.051);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleSize(0.051);

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}

void plotfit()
{
  gStyle->SetOptStat(1000000000);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"X");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetFrameLineWidth(1);

  
 //  TH1D *cf = (TH1D *) gDirectory->Get("CfnReYlm00KmPcent1NonIdYlms;1");
 //  if (cf == 0) {
 //    cf = (TH1D *) gDirectory->Get("CfnReYlm00KpApcent1NonIdYlms;1");
 //  }
 //  if (cf == 0) {
 //    cf = (TH1D *) gDirectory->Get("CfnReYlm00KpApcent2NonIdYlms;1");
 //  }
 //  if (cf == 0) {
 //    cf = (TH1D *) gDirectory->Get("CfnReYlm00KmPcent2NonIdYlms;1");
 //  }
 //   if (cf == 0) 
 //     {
	
 // cf = (TH1D *) gDirectory->Get("CfnReYlm00cylmKpAProtpcM0");
 //     }

  TH1D *cf1 = (TH1D *) gDirectory->Get("datarfpmbg");
  TH1D *cf2 = (TH1D *) gDirectory->Get("datarfmpbg");
  TH1D *cf3 = (TH1D *) gDirectory->Get("dataffpmbg");
  TH1D *cf4 = (TH1D *) gDirectory->Get("dataffmpbg");

  preparehist(cf1);
  preparehist(cf2);
  preparehist(cf3);
  preparehist(cf4);
  
  TGraph *gr = (TGraph *) gDirectory->Get("grappfit;1");

  TCanvas *canfit = new TCanvas("canfit","canfit",1200,950);
  preparepad();
  
  canfit->Divide(2,2,0.0001,0.0001);

  canfit->cd(1); preparepad();

  cf1->Draw();

  gr->SetLineColor(kBlue);
  gr->SetLineStyle(2);
  gr->Draw("CP");


  canfit->cd(2); preparepad();

  cf2->Draw();

  gr->SetLineColor(kBlue);
  gr->SetLineStyle(2);
  gr->Draw("CP");

  canfit->cd(3); preparepad();
  
  cf3->Draw();

  gr->SetLineColor(kBlue);
  gr->SetLineStyle(2);
  gr->Draw("CP");


  canfit->cd(4); preparepad();

  
  cf4->Draw();

  gr->SetLineColor(kBlue);
  gr->SetLineStyle(2);
  gr->Draw("CP");
}
