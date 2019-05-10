void preparepad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.003);
  gPad->SetRightMargin(0.005);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.11);

  gPad->SetGridx();
  gPad->SetGridy();
}

void preparehist(TH1D *hist, int htype=0, int hvar=0)
{
  int cols[4] = { kGreen+2, kRed, kViolet+1, kOrange+7 };
  int mars[4] = { 20, 21, 25, 24 };
  
  if (htype == 0)  {
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.9);
  hist->SetLineColor(kRed);
  hist->SetMarkerColor(kRed);

  hist->GetXaxis()->SetRange(1,200);
  
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
  if (htype == 1)  {
    hist->SetMarkerStyle(mars[hvar]);
    hist->SetMarkerSize(0.9);
    hist->SetLineColor(cols[hvar]);
    hist->SetMarkerColor(cols[hvar]);

    hist->GetXaxis()->SetRange(1,75);
    
    hist->SetTitle(";k* (GeV/#it{c});C(k*)");
    
    hist->GetXaxis()->SetNdivisions(310);
    hist->GetYaxis()->SetNdivisions(310);
    
    hist->GetXaxis()->SetLabelSize(0.051);
    hist->GetYaxis()->SetLabelSize(0.051);
    
    hist->GetXaxis()->SetTitleSize(0.051);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleSize(0.051);
    
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
  }
  
}

void plotpubl()
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

  // Centrality 0
  
  TH1D *cf1;
  TH1D *cf2;
  TH1D *cf3;
  TH1D *cf4;

  TH1D *cf1ss;
  TH1D *cf2ss;
  TH1D *cf3ss;
  TH1D *cf4ss;

  TGraphErrors *grvalerr = (TGraphErrors *) gDirectory->Get("fitresult");
  
  TLatex lat;
  lat.SetTextSize(0.06);
  lat.SetTextFont(42);

  double leglow = 0.6;
  double legtop = 0.95;
  
  for (int icent=0; icent<6; icent++) {
  
    cf1 = (TH1D *) gDirectory->Get(Form("datarfpmbgC%i", icent));
    cf2 = (TH1D *) gDirectory->Get(Form("datarfmpbgC%i", icent));
    cf3 = (TH1D *) gDirectory->Get(Form("dataffpmbgC%i", icent));
    cf4 = (TH1D *) gDirectory->Get(Form("dataffmpbgC%i", icent));

    if (cf1) {
      
      preparehist(cf1,1,0);
      preparehist(cf2,1,1);
      preparehist(cf3,1,2);
      preparehist(cf4,1,3);

      Double_t hulmax = 1.27;
      Double_t hulmin = 0.94;
      
      cf1->SetMaximum(hulmax);
      cf2->SetMaximum(hulmax);
      cf3->SetMaximum(hulmax);
      cf4->SetMaximum(hulmax);

      cf1->SetMinimum(hulmin);
      cf2->SetMinimum(hulmin);
      cf3->SetMinimum(hulmin);
      cf4->SetMinimum(hulmin);
      
      TGraph *grC0 = (TGraph *) gDirectory->Get(Form("graulfitC%i;1", icent));
      grC0->SetLineColor(kBlue);
      grC0->SetLineWidth(2);
      grC0->SetLineStyle(2);

      TCanvas *canfitC0 = new TCanvas(Form("canfitC%i", icent) ,Form("canfitC%i", icent) ,1200,750);
      preparepad();
      
      canfitC0->Divide(2,2,0.0001,0.0001);

      canfitC0->cd(1); preparepad();
      cf1->Draw();
      grC0->Draw("CP");
      
      canfitC0->cd(2); preparepad();
      cf2->Draw();
      grC0->Draw("CP");
      
      canfitC0->cd(3); preparepad();
      cf3->Draw();
      grC0->Draw("CP");
      
      canfitC0->cd(4); preparepad();
      cf4->Draw();
      grC0->Draw("CP");

      canfitC0->cd(1);

      TLegend *leg0 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg0->SetFillColor(0);
      leg0->SetFillStyle(0);
      leg0->SetBorderSize(0);

      leg0->AddEntry(cf1, "p K^{-} field rf", "");
      leg0->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[21+4*icent], grvalerr->GetEY()[21+4*icent]), "");
      leg0->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg0->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg0->Draw();

      canfitC0->cd(2);

      TLegend *leg1 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg1->SetFillColor(0);
      leg1->SetBorderSize(0);
      leg1->SetFillStyle(0);
      
      leg1->AddEntry(cf1, "#bar{p} K^{+} field rf", "");
      leg1->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[22+4*icent], grvalerr->GetEY()[22+4*icent]), "");
      leg1->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg1->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg1->Draw();

      canfitC0->cd(3);

      TLegend *leg2 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      leg2->SetFillStyle(0);

      leg2->AddEntry(cf1, "p K^{-} field ff", "");
      leg2->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[23+4*icent], grvalerr->GetEY()[23+4*icent]), "");
      leg2->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg2->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg2->Draw();

      canfitC0->cd(4);

      TLegend *leg3 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg3->SetFillColor(0);
      leg3->SetBorderSize(0);
      leg3->SetFillStyle(0);
      
      leg3->AddEntry(cf1, "#bar{p} K^{+} field ff", "");
      leg3->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[24+4*icent], grvalerr->GetEY()[24+4*icent]), "");
      leg3->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg3->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg3->Draw();

      canfitC0->cd(1);
      
      TLegend *legf0v = new TLegend(0.2, 0.3, 0.7, leglow-0.02);
      legf0v->SetFillColor(0);
      legf0v->SetBorderSize(0);
      legf0v->SetFillStyle(0);
      
      legf0v->AddEntry(cf1, Form("#Rgothic f_{0} = %.3f #pm %.3f", grvalerr->GetY()[1], grvalerr->GetEY()[1]), "");
      legf0v->AddEntry(cf1, Form("#Jgothic f_{0} = %.3f #pm %.3f", grvalerr->GetY()[2], grvalerr->GetEY()[2]), "");

      legf0v->Draw();

      canfitC0->cd(2);
      
      TLegend *legrad = new TLegend(0.2, 0.3, 0.8, leglow-0.02);
      legrad->SetFillColor(0);
      legrad->SetBorderSize(0);
      legrad->SetFillStyle(0);
      
      legrad->AddEntry(cf1, Form("Radius = %.3f #pm %.3f", grvalerr->GetY()[5+icent*3], grvalerr->GetEY()[5+icent*3]), "");
      //      legrad->AddEntry(cf1, Form("#Jgothic f_{0} = %.3f #pm %.3f", grvalerr->GetY()[2], grvalerr->GetEY()[2]), "");

      legrad->Draw();
      // lat.DrawLatex(0.1, 1.15, "p K^{-}");
      // lat.DrawLatex(0.1, 1.08, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[21], grvalerr->GetEY()[21]));
      
      canfitC0->SaveAs(Form("canfitC%i.root",icent));
      canfitC0->SaveAs(Form("canfitC%i.eps",icent));
      canfitC0->SaveAs(Form("canfitC%i.png",icent));
    }

    cf1ss = (TH1D *) gDirectory->Get(Form("datarfppbgC%i", icent));
    cf2ss = (TH1D *) gDirectory->Get(Form("datarfmmbgC%i", icent));
    cf3ss = (TH1D *) gDirectory->Get(Form("dataffppbgC%i", icent));
    cf4ss = (TH1D *) gDirectory->Get(Form("dataffmmbgC%i", icent));

    if (cf1ss) {
      
      preparehist(cf1ss,1,0);
      preparehist(cf2ss,1,1);
      preparehist(cf3ss,1,2);
      preparehist(cf4ss,1,3);
    
      Double_t hssmax = 1.03;
      Double_t hssmin = 0.77;
      
      cf1ss->SetMaximum(hssmax);
      cf2ss->SetMaximum(hssmax);
      cf3ss->SetMaximum(hssmax);
      cf4ss->SetMaximum(hssmax);

      cf1ss->SetMinimum(hssmin);
      cf2ss->SetMinimum(hssmin);
      cf3ss->SetMinimum(hssmin);
      cf4ss->SetMinimum(hssmin);

      TGraph *grC0ss = (TGraph *) gDirectory->Get(Form("grassfitC%i;1", icent));
      grC0ss->SetLineColor(kBlue);
      grC0ss->SetLineStyle(2);
      grC0ss->SetLineWidth(2);

      TCanvas *canfitC0ss = new TCanvas(Form("canfitC%iss", icent) ,Form("canfitC%iss", icent) ,1200,750);
      preparepad();
      
      canfitC0ss->Divide(2,2,0.0001,0.0001);

      canfitC0ss->cd(1); preparepad();
      cf1ss->Draw();
      grC0ss->Draw("CP");
      
      canfitC0ss->cd(2); preparepad();
      cf2ss->Draw();
      grC0ss->Draw("CP");
      
      canfitC0ss->cd(3); preparepad();
      cf3ss->Draw();
      grC0ss->Draw("CP");
      
      canfitC0ss->cd(4); preparepad();
      cf4ss->Draw();
      grC0ss->Draw("CP");

      canfitC0ss->cd(1);
      
      TLegend *leg0 = new TLegend(0.4, 0.3, 0.95, 0.75);
      leg0->SetFillColor(0);
      leg0->SetBorderSize(0);

      leg0->AddEntry(cf1, "K^{+} p field rf", "");
      leg0->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[45+4*icent], grvalerr->GetEY()[45+4*icent]), "");
      leg0->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg0->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg0->Draw();

      canfitC0ss->cd(2);

      TLegend *leg1 = new TLegend(0.4, 0.3, 0.95, 0.75);
      leg1->SetFillColor(0);
      leg1->SetBorderSize(0);

      leg1->AddEntry(cf1, "K^{-} #bar{p} field rf", "");
      leg1->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[46+4*icent], grvalerr->GetEY()[46+4*icent]), "");
      leg1->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg1->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg1->Draw();

      canfitC0ss->cd(3);

      TLegend *leg2 = new TLegend(0.4, 0.3, 0.95, 0.75);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);

      leg2->AddEntry(cf1, "K^{+} p field ff", "");
      leg2->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[47+4*icent], grvalerr->GetEY()[47+4*icent]), "");
      leg2->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg2->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg2->Draw();

      canfitC0ss->cd(4);

      TLegend *leg3 = new TLegend(0.4, 0.3, 0.95, 0.75);
      leg3->SetFillColor(0);
      leg3->SetBorderSize(0);

      leg3->AddEntry(cf1, " K^{-} #bar{p} field ff", "");
      leg3->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[48+4*icent], grvalerr->GetEY()[48+4*icent]), "");
      leg3->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg3->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg3->Draw();

      canfitC0ss->cd(1);

      TLegend *legrad = new TLegend(0.2, 0.1, 0.8, 0.28);
      legrad->SetFillColor(0);
      legrad->SetBorderSize(0);

      legrad->AddEntry(cf1, Form("BGScale = %.3f #pm %.3f", grvalerr->GetY()[69+icent], grvalerr->GetEY()[69+icent]), "");

      legrad->Draw();
      
      

      canfitC0ss->SaveAs(Form("canfitC%iss.root",icent));
      canfitC0ss->SaveAs(Form("canfitC%iss.eps",icent));
      canfitC0ss->SaveAs(Form("canfitC%iss.png",icent));
    }
  }
  
  // TGraph *grC3 = (TGraph *) gDirectory->Get("grappfitC3;1");
  // grC3->SetLineColor(kBlue);
  // grC3->SetLineStyle(2);

  
  // TH1D *cf1 = (TH1D *) gDirectory->Get("datarfpmbgC3");
  // TH1D *cf2 = (TH1D *) gDirectory->Get("datarfmpbgC3");
  // TH1D *cf3 = (TH1D *) gDirectory->Get("dataffpmbgC3");
  // TH1D *cf4 = (TH1D *) gDirectory->Get("dataffmpbgC3");

  // preparehist(cf1);
  // preparehist(cf2);
  // preparehist(cf3);
  // preparehist(cf4);

  // TCanvas *canfitC3 = new TCanvas("canfitC3","canfitC3",1200,950);
  // preparepad();
  
  // canfitC3->Divide(2,2,0.0001,0.0001);

  // canfitC3->cd(1); preparepad();
  // cf1->Draw();
  // grC3->Draw("CP");

  // canfitC3->cd(2); preparepad();
  // cf2->Draw();
  // grC3->Draw("CP");

  // canfitC3->cd(3); preparepad();
  // cf3->Draw();
  // grC3->Draw("CP");

  // canfitC3->cd(4); preparepad();
  // cf4->Draw();
  // grC3->Draw("CP");

  Double_t xym[6] = { 0.025, 0.075, 0.15, 0.25, 0.35, 0.45 };
  Double_t ykr[6];
  Double_t ykre[6];

  for (int icent=0; icent<6; icent++) {
    ykr[icent] = grvalerr->GetY()[5+icent*3];
    ykre[icent] = grvalerr->GetEY()[5+icent*3];
  }
  
  TGraphErrors *grrads = new TGraphErrors(6, xym, ykr, NULL, ykre);
  
  TCanvas *canradfit = new TCanvas(Form("canradfit") ,Form("canradfit") ,600,380);
  preparepad();

  TH1D *hramka = new TH1D ("hramka",";centrality;R (fm)",100,0.001,0.61);
  preparehist(hramka);
  hramka->SetMaximum(6.1);
  hramka->SetMinimum(1.1);

  
  
  canradfit->cd();
  hramka->Draw();

  grrads->SetMarkerSize(0.8);
  grrads->SetMarkerColor(kRed);
  grrads->SetMarkerStyle(24);

  grrads->Draw("P");

  canradfit->SaveAs("canradpubl.eps");
  canradfit->SaveAs("canradpubl.png");
  canradfit->SaveAs("canradpubl.root");
      
  

}
