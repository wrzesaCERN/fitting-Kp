void preparepad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.003);
  gPad->SetRightMargin(0.11);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.11);
  
}

void preparehist(TH1D *hist)
{
  //  hist->GetXaxis()->SetRange(xbinmin, xbinmax);
  //  hist->GetXaxis()->SetTitle("sign(k*_{Out}) * k*");
  //  hist->GetYaxis()->SetTitle("C(k*)");
  //  hist->SetTitle("");
  // hist->SetMarkerStyle(20);
  // hist->SetMarkerSize(1.0);

  hist->GetXaxis()->SetLabelSize(0.051);
  hist->GetYaxis()->SetLabelSize(0.051);

  hist->GetXaxis()->SetTitleSize(0.051);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleSize(0.051);

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}

void plotgrs()
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

  TFile *infiles[6];

  TGraph *grs[6];

  Double_t rads[6] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

  Double_t af0re = -0.7;
  Double_t af0im = 0.5;

  Int_t cols[6] = { kBlue, kRed, kBlack, kGreen+2, kOrange+7, kViolet+1 };
  
  for (int ifile=0; ifile<6; ifile++)
    {
      infiles[ifile] = new TFile(Form("kpcoulstrong.rad%.1f.re%.3f.im%.3f.root",rads[ifile],af0re, af0im));
      grs[ifile] = (TGraph *) infiles[ifile]->Get("grkaonproton");
      grs[ifile]->SetName(Form("grrad%i",ifile));

      grs[ifile]->SetLineColor(cols[ifile]);
      grs[ifile]->SetMarkerColor(cols[ifile]);
    }

  TCanvas *cankpth = new TCanvas("cankpth","cankpth",1000,700);
  cankpth->cd(); preparepad();
  
  TH1D *hramka = new TH1D ("hramka",";k* (GeV/#it{c});C(k*)",100,0.0,0.29);
  preparehist(hramka);

  hramka->SetMaximum(1.025);
  hramka->SetMinimum(0.901);
  hramka->Draw();

  for (int ifile=0; ifile<6; ifile++)
    grs[ifile]->Draw("CP");
}
