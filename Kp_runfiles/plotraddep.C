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

void preparegraph(TGraph *hist, int iter)
{
  const int cols[8] = { kGreen+2, kRed, kBlue, kCyan-3, kViolet+1, kOrange, kMagenta, kTeal+1 };
  
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.9);
  hist->SetLineColor(cols[iter]);
  hist->SetMarkerColor(cols[iter]);

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

void plotraddep()
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

  TFile *infiles[8];

  infiles[0] = new TFile("kpcoulstrongwgtlhc.rad2.5.re-0.700.im0.800.root");
  infiles[1] = new TFile("kpcoulstrongwgtlhc.rad3.0.re-0.700.im0.800.root");
  infiles[2] = new TFile("kpcoulstrongwgtlhc.rad3.5.re-0.700.im0.800.root");
  infiles[3] = new TFile("kpcoulstrongwgtlhc.rad4.0.re-0.700.im0.800.root");
  infiles[4] = new TFile("kpcoulstrongwgtlhc.rad4.5.re-0.700.im0.800.root");
  infiles[5] = new TFile("kpcoulstrongwgtlhc.rad5.0.re-0.700.im0.800.root");
  infiles[6] = new TFile("kpcoulstrongwgtlhc.rad5.5.re-0.700.im0.800.root");
  infiles[7] = new TFile("kpcoulstrongwgtlhc.rad6.0.re-0.700.im0.800.root");

  TGraph *grs[8];

  for (int iter=0; iter<8; iter++) {
    grs[iter] = (TGraph *) infiles[iter]->Get("grkaonproton");

    preparegraph(grs[iter], iter);
  }

  TCanvas *canraddep = new TCanvas("canraddep" ,"canraddep" ,1200,750);
  preparepad();

  TH1D *ramka = new TH1D("ramka",";k^{*} (GeV/c);C(k^{*})",100,0.0,0.4);
  preparehist(ramka);
  ramka->SetMaximum(1.59);
  ramka->SetMinimum(0.72);

  canraddep->cd();
  ramka->Draw();

  for (int iter=0; iter<8; iter++) {
    grs[iter]->Draw("C");
  }

  TLegend *leg = new TLegend(0.2,0.45, 0.8, 0.95);
  leg->SetBorderSize(0);
  leg->AddEntry(grs[0], "K^{-}p Re f_{0} = -0.7 fm, Im f_{0} = 0.8 fm","");
  leg->AddEntry(grs[0], "r_{0} = 2.5  fm", "l");
  leg->AddEntry(grs[1], "r_{0} = 3.0  fm", "l");
  leg->AddEntry(grs[2], "r_{0} = 3.5  fm", "l");
  leg->AddEntry(grs[3], "r_{0} = 4.0  fm", "l");
  leg->AddEntry(grs[4], "r_{0} = 4.5  fm", "l");
  leg->AddEntry(grs[5], "r_{0} = 5.0  fm", "l");
  leg->AddEntry(grs[6], "r_{0} = 5.5  fm", "l");
  leg->AddEntry(grs[7], "r_{0} = 6.0  fm", "l");
  
  leg->Draw();

  canraddep->SaveAs("canraddep.root");
  canraddep->SaveAs("canraddep.eps");
  canraddep->SaveAs("canraddep.png");
  
}
