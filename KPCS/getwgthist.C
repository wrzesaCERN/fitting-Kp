void getwgthist()
{
  //  TFile *infile = new TFile("Therm2EventsPiKProofOutput.b3.R-1.0.I0.45.flt.root");
  TFile *infile = new TFile("Therm2EventsPiKProofOutput.b3.R-1.1.I0.45.flt.root");
  TH2D *hcos = infile->Get("hKStarRStar");

  Double_t rsum;
  
  for (int ix=1; ix<=hcos->GetNbinsX(); ix++) {
    rsum = 0;
    for (int iy=1; iy<=hcos->GetNbinsY(); iy++) {
      rsum += hcos->GetBinContent(ix, iy);
    }
    for (int iy=1; iy<=hcos->GetNbinsY(); iy++) {
      hcos->SetBinContent(ix, iy, hcos->GetBinContent(ix, iy)/rsum);
    }
  }

  TFile *outfile = new TFile("hcoswghts.root","RECREATE");
  outfile->cd();
  hcos->Write();
}
