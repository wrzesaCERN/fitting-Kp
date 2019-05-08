int getfuns()
{
  TFile *infile1 = new TFile("AnalysisResults_tpconly_c0.plus.dir.root");
  TFile *infile2 = new TFile("AnalysisResults_tpconly_c1.plus.dir.root");
  TFile *infile3 = new TFile("AnalysisResults_tpconly_c3.plus.dir.root");
  TFile *infile4 = new TFile("AnalysisResults_tpconly_c5.plus.dir.root");

  TH1D *nums[4*7];
  TH1D *dens[4*7];

  TH1D *corr[4*7];

  const char* pairs[4] = { "KpProt", "KmProt", "KpAProt", "KmAProt" };
  
  Double_t scale;

  TFile *outfile = new TFile("corrs.root","RECREATE");
  outfile->cd();

  for (int icent=0; icent<6; icent++) {
    for (int iter=0; iter<4; iter++) {
      if (icent == 0) {
	cout << "Reading " << Form("NumReYlm00cylm%stpcM%i",pairs[iter],icent) << endl;
	nums[icent*4+iter] = (TH1D *) infile1->Get(Form("NumReYlm00cylm%spcM%i",pairs[iter],icent));
	dens[icent*4+iter] = (TH1D *) infile1->Get(Form("DenReYlm00cylm%spcM%i",pairs[iter],icent));
      }
      if ((icent == 1) || (icent == 2)) {
	cout << "Reading " << Form("NumReYlm00cylm%stpcM%i",pairs[iter],icent) << endl;
	nums[icent*4+iter] = (TH1D *) infile2->Get(Form("NumReYlm00cylm%spcM%i",pairs[iter],icent));
	dens[icent*4+iter] = (TH1D *) infile2->Get(Form("DenReYlm00cylm%spcM%i",pairs[iter],icent));
      }
      if ((icent == 3) || (icent == 4)) {
	nums[icent*4+iter] = (TH1D *) infile3->Get(Form("NumReYlm00cylm%spcM%i",pairs[iter],icent));
	dens[icent*4+iter] = (TH1D *) infile3->Get(Form("DenReYlm00cylm%spcM%i",pairs[iter],icent));
      }
      if ((icent == 5) || (icent == 6)) {
	nums[icent*4+iter] = (TH1D *) infile4->Get(Form("NumReYlm00cylm%spcM%i",pairs[iter],icent));
	dens[icent*4+iter] = (TH1D *) infile4->Get(Form("DenReYlm00cylm%spcM%i",pairs[iter],icent));
      }

      // cout << "Got num den " << icent << " " << iter << endl;
      // cout << nums[icent*4+iter] << endl;
      // cout << dens[icent*4+iter] << endl;
      
      // cout << nums[icent*4+iter]->GetEntries() << endl;
      
      // cout << dens[icent*4+iter]->GetEntries() << endl;

      cout << nums[icent*4+iter]->GetBinContent(3)/dens[icent*4+iter]->GetBinContent(3) << endl;
      cout << nums[icent*4+iter]->GetBinContent(7)/dens[icent*4+iter]->GetBinContent(7) << endl;
      
      corr[icent*4+iter] = new TH1D (*nums[icent*4+iter]);
      corr[icent*4+iter]->SetName(Form("CfnReYlm00cylm%spcM%i",pairs[iter], icent));
      corr[icent*4+iter]->SetTitle(Form("CfnReYlm00cylm%spcM%i;k* (GeV/c);C(k*)",pairs[iter],icent));
      
      scale = dens[icent*4+iter]->Integral(40,200)/nums[icent*4+iter]->Integral(40,200);
      
      cout << "Scale is " << scale << endl;
      
      corr[icent*4+iter]->Divide(nums[icent*4+iter], dens[icent*4+iter],1.0, 1.0);
      cout << corr[icent*4+iter]->GetBinContent(3) << " " << corr[icent*4+iter]->GetBinContent(7) << endl;
      corr[icent*4+iter]->Scale(scale);
      cout << corr[icent*4+iter]->GetBinContent(3) << " " << corr[icent*4+iter]->GetBinContent(7) << endl;

      
    }
  }

  for (int iter=0; iter<24; iter++) {
    corr[iter]->Write();
    cout << "Writing " << endl;
    cout << corr[iter]->GetBinContent(3) << " " << corr[iter]->GetBinContent(7) << endl;

    cout << "Written " << endl;

  }

  outfile->Close();
  
  return 0;
}
