#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TGraph.h>

#include <TGraph2D.h>

//#define RMAX 9
#define RMAX 9
#define IMAX 8

using namespace std;

Double_t rvals[RMAX] = { -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.001, 0.1, 0.2 };
//Double_t rvals[RMAX] = { -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, 0.001, 0.5 };
Double_t ivals[IMAX] = { 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };

TGraph *cfs00sbgs[IMAX*RMAX];
TGraph *cfs00sbgsmr[IMAX*RMAX];
TH1D *cfs00ppbg;

TH1D *cfs00ppbgmr;

TH1D *cdata;
TH1D *cdatapp;

TGraph *smearmr(TGraph *cfppdv)
{
  double width = 0.003;
  //  double purity = 0.9;
  double kvt, kvm, wgt, wgtsum, valsum;

  TGraph *cfppmr = new TGraph(*cfppdv);
  cfppmr->SetName(Form("%s_mr",cfppdv->GetName()));
  
  double bwdt = cfppdv->GetX()[1] - cfppdv->GetX()[0];

  for (int ib=0; ib<cfppdv->GetN(); ib++) {
    wgtsum = 0.0;
    valsum = 0.0;
    kvm = cfppdv->GetX()[ib];
    for (int ix = -17; ix<18; ix++) {
      kvt = cfppmr->GetX()[TMath::Abs(ib+ix)+1];
      wgt = TMath::Exp(-(ix*bwdt)*(ix*bwdt)/(2*width*width))*kvt*kvt;
      if ((ib+ix)>cfppdv->GetN()) wgt = 0;
      valsum += cfppdv->GetY()[TMath::Abs(ib+ix)+1]*wgt;
      wgtsum += wgt;
    }

    cfppmr->GetY()[ib] = valsum/wgtsum;
  }

  return cfppmr;
}


void getlinext(double rev, double imv, double *vals)
{
  int rbin=0, ibin=0;

  if (rev<rvals[0]) rbin = 0;
  else if (rev>rvals[RMAX-1]) rbin=RMAX-2;
  else {
    for (rbin=0; rbin<RMAX; rbin++)
      if (rev < rvals[rbin]) {rbin--; break;}
    //      else cout << " Next r" << endl;

  }

  if (imv<ivals[0]) ibin = 0;
  else if (imv>ivals[IMAX-1]) ibin=IMAX-2;
  else {
    for (ibin=0; ibin<IMAX; ibin++)
      if (imv < ivals[ibin]) {ibin--; break;}
    //      else cout << " Next i" << endl;
      
  }

  // cout << "r i " << rbin << " " << ibin << endl;
  
  // cout << " " << rev << " " << rvals[rbin] << "  - "  << rvals[rbin+1] << endl;
  // cout << " " << imv << " " << ivals[ibin] << "  - "  << ivals[ibin+1] << endl;
  
  // TH1D *cf00 = cfs00sbgs[rbin*IMAX + ibin];
  // TH1D *cf01 = cfs00sbgs[rbin*IMAX + ibin+1];
  // TH1D *cf10 = cfs00sbgs[(rbin+1)*IMAX + ibin];
  // TH1D *cf11 = cfs00sbgs[(rbin+1)*IMAX + ibin+1];
  
  TGraph *cf00 = cfs00sbgsmr[rbin*IMAX + ibin];
  TGraph *cf01 = cfs00sbgsmr[rbin*IMAX + ibin+1];
  TGraph *cf10 = cfs00sbgsmr[(rbin+1)*IMAX + ibin];
  TGraph *cf11 = cfs00sbgsmr[(rbin+1)*IMAX + ibin+1];
  
  double wgt00, wgt01, wgt10, wgt11;
  double intv;

  wgt11 = (rev - rvals[rbin])  *(imv  - ivals[ibin])  /((rvals[rbin+1] - rvals[rbin])*(ivals[ibin+1]-ivals[ibin]));
  wgt10 = (rev - rvals[rbin])  *(-imv  + ivals[ibin+1])/((rvals[rbin+1] - rvals[rbin])*(ivals[ibin+1]-ivals[ibin]));
  wgt01 = (-rev + rvals[rbin+1])*(imv  - ivals[ibin])  /((rvals[rbin+1] - rvals[rbin])*(ivals[ibin+1]-ivals[ibin]));
  wgt00 = (-rev + rvals[rbin+1])*(-imv + ivals[ibin+1])/((rvals[rbin+1] - rvals[rbin])*(ivals[ibin+1]-ivals[ibin]));

  for (int ib=0; ib<cf00->GetN(); ib++) 
    {
      intv = (wgt00 * cf00->GetY()[ib] +
	      wgt01 * cf01->GetY()[ib] +
	      wgt10 * cf10->GetY()[ib] +
	      wgt11 * cf11->GetY()[ib]);
      vals[ib-1] = intv;

//       cout << "  wgt " << wgt00 << " " << wgt01 << " " << wgt10 << " " << wgt11 << endl;
//       cout << "int " << cf00->GetBinContent(ib) << " " << cf01->GetBinContent(ib) << " " << cf10->GetBinContent(ib) << " " << cf11->GetBinContent(ib) << "     " << intv << endl;
    }
}

double getchi2(TH1D *data, TH1D *model, double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;

  int maxbin = 65;
  int minbin = 3;

  for (int ib=minbin; ib<=maxbin; ib++)
    {
      if (ib == 47) continue;
      if (ib == 48) continue;
      if (ib == 49) continue;
      
      val = data->GetBinContent(ib);
      ksv = data->GetXaxis()->GetBinCenter(ib);
      bgv = norm + lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
      val += -bgv + 1.0;
      val = (val-1.0)/purity + 1.0;

      thv = model->GetBinContent(ib);
      err = data->GetBinError(ib);

      chi = (val-thv)*(val-thv)/(err*err);
      
      //      cout << "val thv chi " << val << "   " << thv << "   " << chi << endl;

      sumchi += chi;
    }

  cout << "full chi2 " << sumchi << "     ----       " << norm << " "<< lambda << " " << radbg << " " << purity << endl;

  return sumchi;
}

double getchi2(TH1D *data, double *model, double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;

  int maxbin = 65;
  int minbin = 3;

  int lobin, hibin;

  //  cout << "Fitting for n p l r " << norm << " " << purity << " " << lambda << " "<< radbg << endl;
  
  Double_t lokval, hikval;
  
  for (int ib=minbin; ib<=maxbin; ib++)
    {
      if (ib == 47) continue;
      if (ib == 48) continue;
      if (ib == 49) continue;

      val = data->GetBinContent(ib);
      ksv = data->GetXaxis()->GetBinCenter(ib);
      
      //      cout << "Finding bin " << endl;
      lobin = (int) floor((ksv-0.001)/0.002);
      hibin = lobin+1;

      lokval = lobin*0.002+0.001;
      hikval = hibin*0.002+0.001;

      thv = (model[lobin] * (hikval-ksv)/0.002 +
	     model[hibin] * (ksv-lokval)/0.002);
      
      // cout << "Interpolating between " << lobin << " " << lokval << " " << model[lobin] << endl;
      // cout << "                  and " << hibin << " " << hikval << " " << model[hibin] << endl;
      // cout << "                   is " << thv << endl;
      // cout << "                  for " << ksv << endl;
      
      bgv = norm + lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
      val += -bgv + 1.0;
      val = (val-1.0)/purity + 1.0;

      err = data->GetBinError(ib);

      chi = (val-thv)*(val-thv)/(err*err);
      
      //      cout << "val thv chi " << val << "   " << thv << "   " << chi << "   ----    " << norm << " " << purity << " " << lambda << " " << radbg <<  endl;

      sumchi += chi;
    }

  cout << "full chi2 " << sumchi << "     ----       " << norm << " "<< lambda << " " << radbg << " " << purity << endl;
  //  cout << "full chi2 " << sumchi << endl;

  return sumchi;
}

double getchi2pp(double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;

  int maxbin = 20;
  int minbin = 3;

  for (int ib=minbin; ib<=maxbin; ib++)
    {
      val = cdatapp->GetBinContent(ib);
      ksv = cdatapp->GetXaxis()->GetBinCenter(ib);
      bgv = norm + 1.4*lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
      val += -bgv + 1.0;
      val = (val-1.0)/purity + 1.0;

      thv = cfs00ppbgmr->GetBinContent(ib);
      err = cdatapp->GetBinError(ib);

      chi = (val-thv)*(val-thv)/(err*err);
      
      sumchi += chi;
    }

  return sumchi;

}

void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t Nv  = par[0];
  Double_t Lv  = par[1];
  Double_t Rv  = par[2];
  Double_t Pv  = par[3];
  
  Double_t Rev = par[4];
  Double_t Imv = par[5];

  Double_t chi2 = 0.0;
  Double_t tv, qv;
    
  Double_t vals[500];

  cout << "Funcf called " << Nv << " " << Lv << " " << Rv << " " << Pv << " " << Rev << " " << Imv << endl;
  
  getlinext(Rev, Imv, vals);

  chi2 = getchi2(cdata, vals, Nv, Pv, Lv, Rv);

  //  cout << " chi2 now " << chi2 << endl;

  //  chi2 += getchi2pp(Nv, Pv, Lv, Rv);

  //  cout << " chi2 now " << chi2 << endl;

  //  cout << "Funcf called " << Nv << " " << Lv << " " << Rv << " " << Pv << " " << Rev << " " << Imv << "    " << chi2 << endl;

  f = chi2;
}


int main(int argc, char **argv)
{
  // Declare functions for interpolation
  
  TFile *inthfiles[IMAX*RMAX];
  TGraph *cfs00s[IMAX*RMAX];

  // Get baseline file
  
  TFile *inbasefile = new TFile("shmout.recalc.kp.base.root");
  TH1D *cbase = (TH1D *) inbasefile->Get("CfnReYlm00CYlmBase");

  Double_t val;

  // Load functions for interpolation
  // And divide by baseline
  
  for (int ir=0; ir<RMAX; ir++)
    for (int ii=0; ii<IMAX; ii++) {
      inthfiles[ir*IMAX + ii] = new TFile (Form("kpcoulstrong.rad4.0.re%.3f.im%.3f.root",rvals[ir], ivals[ii]));
      //      inthfiles[ir*IMAX + ii] = new TFile (Form("shmout.recalc.kp.R%i.I%i.root",ir+1+(ir>1)*1, ii+1));

      cfs00s[ir*IMAX + ii] = (TGraph *) inthfiles[ir*IMAX + ii]->Get("grkaonproton");

      if (cfs00s[ir*IMAX + ii]) {
	cfs00sbgs[ir*IMAX + ii] = new TGraph(*cfs00s[ir*IMAX + ii]);
	//	cfs00sbgs[ir*IMAX + ii]->Divide(cbase);
	cfs00sbgsmr[ir*IMAX + ii] = smearmr(cfs00sbgs[ir*IMAX + ii]);
      }
      else {
	cfs00sbgs[ir*IMAX + ii] = 0;
	cfs00sbgsmr[ir*IMAX + ii] = 0;
      }
    }

  // Get data for plus-plus combination
  
  TFile *infilepp = new TFile("shmout.recalc.kp.pp.root");
  TH1D *cfs00pp = (TH1D *) infilepp->Get("CfnReYlm00CYlm");
  cfs00ppbg = new TH1D(*cfs00pp);
  cfs00ppbg->Divide(cbase);

  // Smear by momentum resolution
  
  //  cfs00ppbgmr = smearmr(cfs00ppbg);
  cfs00ppbgmr = new TH1D(*cfs00ppbg);
  
  // Get data for plus-minus combination
  
  // TFile *indatafile = new TFile("shmout.recalc.c1.mp.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KmPcent1NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.c1.pm.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KpApcent1NonIdYlms");
  TFile *indatafile = new TFile("shmout.recalc.c2.mp.root");
  cdata = (TH1D *) indatafile->Get("CfnReYlm00KmPcent2NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.c2.pm.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KpApcent2NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.mpos.m1.kmpp.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKmProtpcM1");
  // TFile *indatafile = new TFile("shmout.recalc.ff.cut1.KmPp.m0.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKmProtpcM0");
  // TFile *indatafile = new TFile("shmout.recalc.ff.cut1.KpPm.m0.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKpAProtpcM0");
  // TFile *indatafile = new TFile("shmout.recalc.mpos.m1.kppm.root");
  //cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKpAProtpcM1");

  Double_t xmax = cdata->GetXaxis()->GetXmax();
  
  TF1 *funn = new TF1("funn","[0]");
  funn->SetParameter(0,1.000);
  cdata->Fit(funn, "", "", 0.75*xmax, xmax);
  
  TF1 *fung = new TF1("fung","[0] + [1]*exp(-x*x*[2]*[2])");
  fung->SetParameters(funn->GetParameter(0),0.002,2.0);
  fung->FixParameter(0,funn->GetParameter(0));
  //  fung->FixParameter(2,2.0);
  cdata->Fit(fung, "","",0.25, xmax);

  // Get data for plus-plus combination
  
  // TFile *indatppfile = new TFile("shmout.recalc.mpos.m1.kppp.root");
  // cdatapp = (TH1D *) indatppfile->Get("CfnReYlm00cylmKpProtpcM1");
  TFile *indatppfile = new TFile("shmout.recalc.c1.pp.root");
  cdatapp = (TH1D *) indatppfile->Get("CfnReYlm00KpPcent1NonIdYlms");

  // Create Chi2 map
  
  TGraph2D *grchi = new TGraph2D();
  grchi->SetName("grchi");
  grchi->SetTitle("#chi^{2} Map;#Rgothicf_{0} (fm);#Jgothicf_{0} (fm);#chi^{2}");

  double chival;
  int pcount = 0;

  double minchi = 1e9;

  int minir, minii;

  double purval = 0.80;
  double nrmval = fung->GetParameter(0);
  double lbdval = fung->GetParameter(1)*0.3;
  double wdtval = fung->GetParameter(2);

  // Create a fit chi2 map
  
  for (int ir=0; ir<RMAX; ir++) 
    for (int ii=0; ii<IMAX; ii++) {
      if (cfs00sbgsmr[ir*IMAX + ii]) {
	//	chival = getchi2(cdata, cfs00sbgsmr[ir*IMAX + ii], nrmval, purval, lbdval, wdtval);
	chival = getchi2(cdata, cfs00sbgsmr[ir*IMAX + ii]->GetY(), nrmval, purval, lbdval, wdtval);
	grchi->SetPoint(pcount++, rvals[ir], ivals[ii], chival);

	if (chival < minchi) 
	  { minchi = chival;
	    minir = ir;
	    minii = ii;
	    cout << "Now minimum is " << minchi << endl;
	  }
      }
    }

  // Create the best-fir graph
  
  TH1D *cfit = new TH1D(*cdata);
  cfit->SetName("databg");

  double ksv, bgv;
  for (int ib=1; ib<=cdata->GetNbinsX(); ib++)
    {
      val = cdata->GetBinContent(ib);
      ksv = cdata->GetXaxis()->GetBinCenter(ib);
      bgv = nrmval + lbdval*TMath::Exp(-ksv*ksv*wdtval*wdtval);
      
      val += -bgv + 1.0;
      val = (val-1.0)/purval + 1.0;

      cfit->SetBinContent(ib, val);
    }

  TFile *outfile = new TFile("outchi.root","RECREATE");
  outfile->cd();
  grchi->Write();

  cdata->Write();

  cfs00sbgs[minir*IMAX + minii]->SetName("chiminbg");
  cfs00sbgs[minir*IMAX + minii]->Write();

  cfs00sbgsmr[minir*IMAX + minii]->SetName("chiminbgmr");
  cfs00sbgsmr[minir*IMAX + minii]->Write();

  cfs00s[minir*IMAX + minii]->SetName("chimin");
  cfs00s[minir*IMAX + minii]->Write();

  cfit->Write();
  
  
  cout << "Minimum chi " << minchi << endl;

  cout << "Re Im " << rvals[minir] << " " << ivals[minii] << endl;

  if (minir == RMAX-1) minir--;
  if (minii == IMAX-1) minii--;
  
  double vals[500];
  getlinext(-0.32, 0.55, vals);

  //  exit(0);

  cout << "Now fitting " << endl;

  cout << "Setting values " << nrmval << " " << lbdval << " " << wdtval << " " << purval << endl;
  TVirtualFitter *fitter=TVirtualFitter::Fitter(0, 9);
  fitter->SetFCN(myfuncf);
  fitter->SetParameter(0, "Norm"    ,nrmval        ,0.001, 0.9,   1.1);
  fitter->SetParameter(1, "Lambda"  ,lbdval        ,0.001, 0.0,   0.1);
  fitter->SetParameter(2, "Width"   ,wdtval        ,0.001, 1.0,   10.0);
  fitter->SetParameter(3, "Purity"  ,purval        ,0.0001,0.2,   1.0);
  fitter->SetParameter(4, "Ref0"    ,rvals[minir]  ,0.0001,-5.0, 5.0);
  fitter->SetParameter(5, "Imf0"    ,ivals[minii]  ,0.0001,-5.0, 5.0);
  //  fitter->SetParameter(5, "Imf0"    ,0.25  ,0.001,0.0, 5.0);
  //  fitter->SetParameter(6, "norm"     ,1.0       ,0.0001,0.0,   2.0);


  
  Double_t arglist[100];
  arglist[0] = 1;
  fitter->ExecuteCommand("CALL FCN", arglist, 1);
  fitter->FixParameter(0);
  //  fitter->FixParameter(1);
  fitter->FixParameter(2);
  //  fitter->FixParameter(3);
  //  fitter->FixParameter(4);
  //  fitter->FixParameter(5);
//   fitter->FixParameter(3);
//   fitter->FixParameter(11);
//   fitter->FixParameter(14);

//   fitter->FixParameter(3);
//   fitter->FixParameter(4);
//   fitter->FixParameter(15);
  
//    fitter->FixParameter(7);
//    fitter->FixParameter(8);
//   fitter->FixParameter(9);
//   fitter->FixParameter(10);
//   fitter->FixParameter(12);
//   fitter->FixParameter(13);
//   fitter->FixParameter(14);

  // fitter->SetParameter(5,0.25);
  // fitter->FixParameter(5);

  //fitter->FixParameter(15);
  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT", arglist, 1);
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
  
//   f1 ->SetParameters(fitter->GetParameter(0), fitter->GetParameter(1), fitter->GetParameter(2), 1.0, n0);
//   f20->SetParameters(Ro, Rs, Rl, 1.0, n0*n1);
//   f22->SetParameters(Ro, Rs, Rl, 1.0, n0*n2);

  cout << "Fitting done " << " " << fitter->GetParameter(4) << " " << fitter->GetParameter(5) << endl;


  double fitnorm = fitter->GetParameter(0);
  double fitlabd = fitter->GetParameter(1);
  double fitwdth = fitter->GetParameter(2);
  double fitpuri = fitter->GetParameter(3);

  double fitref0 = fitter->GetParameter(4);
  double fitimf0 = fitter->GetParameter(5);

  cout << "Final paramerers:" << endl;
  cout << "Normalization   : " << fitnorm << endl;
  cout << "BG lambda       : " << fitlabd << endl;
  cout << "BG width        : " << fitwdth << endl;
  cout << "Purity          : " << fitpuri << endl;
  cout << "Re f0           : " << fitref0 << endl;
  cout << "Im f0           : " << fitimf0 << endl;
  
  outfile->cd();

  getlinext(fitref0, fitimf0, vals);

  double xy[200];
  for (int ib=1; ib<=cdata->GetNbinsX(); ib++)
    xy[ib-1] = cdata->GetXaxis()->GetBinCenter(ib);

  TGraph *grfit = new TGraph(500, cfs00s[0]->GetX(), vals);
  grfit->SetTitle("grappfit");
  grfit->SetName("grappfit");

  grfit->Write();
  cfs00ppbg->SetName("cfs00ppbg");
  cfs00ppbg->Write();
  cfs00ppbgmr->SetName(Form("%s_mr",cfs00ppbg->GetName()));
  cfs00ppbgmr->Write();

  TH1D *cfitpmlin = new TH1D(*cdata);
  cfitpmlin->SetName("datapmbg");

  for (int ib=1; ib<=cdata->GetNbinsX(); ib++)
    {
      val = cdata->GetBinContent(ib);
      ksv = cdata->GetXaxis()->GetBinCenter(ib);
      bgv = fitnorm + fitlabd*TMath::Exp(-ksv*ksv*fitwdth*fitwdth);
      
      val += -bgv + 1.0;
      val = (val-1.0)/fitpuri + 1.0;

      //      cout << "val " << val << "    ----     " << nrmval << " " << purval << " " << lbdval << " " << wdtval << endl; 
      
      cfitpmlin->SetBinContent(ib, val);
    }

  cfitpmlin->Write();

  TH1D *cfitpplin = new TH1D(*cdatapp);
  cfitpplin->SetName("datappbg");

  for (int ib=1; ib<=cdatapp->GetNbinsX(); ib++)
    {
      val = cdatapp->GetBinContent(ib);
      ksv = cdatapp->GetXaxis()->GetBinCenter(ib);
      bgv = nrmval + lbdval*TMath::Exp(-ksv*ksv*wdtval*wdtval);
      
      val += -bgv + 1.0;
      val = (val-1.0)/purval + 1.0;

      cfitpplin->SetBinContent(ib, val);
    }

  cfitpplin->Write();

  
}
 
