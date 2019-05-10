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

void preparehist(TH1D *hist)
{
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

void plotfitcentall()
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
      
      preparehist(cf1);
      preparehist(cf2);
      preparehist(cf3);
      preparehist(cf4);
    
      TGraph *grC0 = (TGraph *) gDirectory->Get(Form("graulfitC%i;1", icent));
      grC0->SetLineColor(kBlue);
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

      leg1->AddEntry(cf1, "#bar{p} K^{+} field rf", "");
      leg1->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[22+4*icent], grvalerr->GetEY()[22+4*icent]), "");
      leg1->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg1->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg1->Draw();

      canfitC0->cd(3);

      TLegend *leg2 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);

      leg2->AddEntry(cf1, "p K^{-} field ff", "");
      leg2->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[23+4*icent], grvalerr->GetEY()[23+4*icent]), "");
      leg2->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg2->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg2->Draw();

      canfitC0->cd(4);

      TLegend *leg3 = new TLegend(0.4, leglow, 0.95, 0.95);
      leg3->SetFillColor(0);
      leg3->SetBorderSize(0);

      leg3->AddEntry(cf1, "#bar{p} K^{+} field ff", "");
      leg3->AddEntry(cf1, Form("Norm: %.5f #pm %.5f", grvalerr->GetY()[24+4*icent], grvalerr->GetEY()[24+4*icent]), "");
      leg3->AddEntry(cf1, Form("Lmbd: %.5f #pm %.5f", grvalerr->GetY()[3+3*icent], grvalerr->GetEY()[3+3*icent]), "");
      leg3->AddEntry(cf1, Form("Wdth: %.5f #pm %.5f", grvalerr->GetY()[4+3*icent], grvalerr->GetEY()[4+3*icent]), "");
      
      leg3->Draw();

      canfitC0->cd(1);
      
      TLegend *legf0v = new TLegend(0.2, 0.3, 0.7, leglow-0.02);
      legf0v->SetFillColor(0);
      legf0v->SetBorderSize(0);

      legf0v->AddEntry(cf1, Form("#Rgothic f_{0} = %.3f #pm %.3f", grvalerr->GetY()[1], grvalerr->GetEY()[1]), "");
      legf0v->AddEntry(cf1, Form("#Jgothic f_{0} = %.3f #pm %.3f", grvalerr->GetY()[2], grvalerr->GetEY()[2]), "");

      legf0v->Draw();

      canfitC0->cd(2);
      
      TLegend *legrad = new TLegend(0.2, 0.3, 0.8, leglow-0.02);
      legrad->SetFillColor(0);
      legrad->SetBorderSize(0);

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
      
      preparehist(cf1ss);
      preparehist(cf2ss);
      preparehist(cf3ss);
      preparehist(cf4ss);
    
      TGraph *grC0ss = (TGraph *) gDirectory->Get(Form("grassfitC%i;1", icent));
      grC0ss->SetLineColor(kBlue);
      grC0ss->SetLineStyle(2);

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


}
