int getfuns()
{
  TFile *infile = new TFile("CollectivelyKP.dir.root");

  TH1D *nums[4];
  TH1D *dens[4];

  TH1D *corr[4];

  const char* pairs[4] = { "KpProt", "KmProt", "KpAProt", "KmAProt" };

  Double_t scale;
  
  for (int iter=0; iter<4; iter++) {
    nums[iter] = (TH1D *) infile->Get(Form("NumReYlm00cylm%spcM0",pairs[iter]));
    dens[iter] = (TH1D *) infile->Get(Form("DenReYlm00cylm%spcM0",pairs[iter]));

    corr[iter] = new TH1D (*nums[iter]);
    corr[iter]->SetName(Form("CfnReYlm00cylm%spcM0",pairs[iter]));
    corr[iter]->SetTitle(Form("CfnReYlm00cylm%spcM0;k* (GeV/c);C(k*)",pairs[iter]));

    scale = dens[iter]->Integral(40,200)/nums[iter]->Integral(40,200);

    cout << "Scale is " << scale << endl;
    
    corr[iter]->Divide(nums[iter], dens[iter],1.0, 1.0);
    corr[iter]->Scale(scale);
  }

  TFile *outfile = new TFile("corrs.root","RECREATE");
  outfile->cd();
  for (int iter=0; iter<4; iter++)
    corr[iter]->Write();
  
  
  return 0;
}
