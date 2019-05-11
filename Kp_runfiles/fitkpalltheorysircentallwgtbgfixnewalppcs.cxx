#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>

#define RMAX 9
#define IMAX 10
#define SMAX 8

using namespace std;

//Re_scatlength - real part of scattering length
//Im_scatlength - imaginary part of scattering length
//source - source size
Double_t Re_scatlength[RMAX] = { -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.001, 0.1, 0.2 }; 
Double_t Im_scatlength[IMAX] = { 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 }; 
Double_t source[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0 }; 

//docent - centralities which we take into acount (1-yes, 0-no)
Int_t docent[7] = { 1, 1, 1, 1, 1, 1, 1};

//th_oppo - calculated theoretical correlation functions (for opposite sign pairs)
//th_same - calculated theoretical correlation functions (for same sign pairs)
//th_oppo_smr, th_same_smr - calculated theoretical correlation functions smeared by MR (assumed by Gaussian distribution)
TGraph *th_oppo[IMAX*RMAX*SMAX];
TGraph *th_oppo_smr[IMAX*RMAX*SMAX];

TGraph *th_same[RMAX]; 
TGraph *th_same_smr[RMAX]; 

//below data
//data: all pairs; both fields: negative - mf, positive - pf; particle charge : plus -p and minus - m (i.e. pm -> kaon minus and antideuteron)
TH1D *cdatapfpm[7];
TH1D *cdatapfmp[7];
TH1D *cdatapfpp[7];
TH1D *cdatapfmm[7];

TH1D *cdatamfpm[7];
TH1D *cdatamfmp[7];
TH1D *cdatamfpp[7];
TH1D *cdatamfmm[7];

//Gaussian smearing function MR
TGraph *smearmr(TGraph *cfppdv)
{
  double width = 0.0041;
  double kvt, kvm, wgt, wgtsum, valsum;

  TGraph *cfppmr = new TGraph(*cfppdv);
  cfppmr->SetName(Form("%s_mr",cfppdv->GetName()));
  
  double bwdt = cfppdv->GetX()[1] - cfppdv->GetX()[0];

  for (int ib=0; ib<cfppdv->GetN(); ib++) {
    wgtsum = 0.0;
    valsum = 0.0;
    kvm = cfppdv->GetX()[ib];
    for (int ix = -17; ix<18; ix++) {
      if ((TMath::Abs(ib+ix)+1) > (cfppmr->GetN()-1)) continue;
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

//Get linear extrapolation. Linear interpretation of the linear correlation function between points on the grid 3D (each point is a 1D function in k*). In other words, we multiply every point of the theoretical function by linear coefficients and we get functions at any point in space (R, Im, Re).
void getlinexth_oppo(double rev, double imv, double siz, double *vals) 
{

  //coefficient corresponding to the elements of the Re_scatlength, Im_scatlength and source arrays which are right for "rev", "iev" and "siz" (**)
  int rbin=0, ibin=0, sbin=0;

  if (rev<Re_scatlength[0]) rbin = 0;
  else if (rev>Re_scatlength[RMAX-1]) rbin=RMAX-2;
  else 
    for (rbin=0; rbin<RMAX; rbin++)
      if (rev < Re_scatlength[rbin]) {rbin--; break;}

  if (imv<Im_scatlength[0]) ibin = 0;
  else if (imv>Im_scatlength[IMAX-1]) ibin=IMAX-2;
  else 
    for (ibin=0; ibin<IMAX; ibin++)
      if (imv < Im_scatlength[ibin]) {ibin--; break;}

  if (siz<source[0]) sbin = 0;
  else if (siz>source[SMAX-1]) sbin=SMAX-2;
  else 
    for (sbin=0; sbin<SMAX; sbin++)
      if (siz < source[sbin]) {sbin--; break;}
  

//theoretical function for selected rbin, ibin and sbin and their neighbours (in other words, we select the theoretical functions closest to our values)
  TGraph *cf000 = th_oppo_smr[sbin*IMAX*RMAX + rbin*IMAX + ibin];
  TGraph *cf001 = th_oppo_smr[sbin*IMAX*RMAX + rbin*IMAX + ibin+1];
  TGraph *cf010 = th_oppo_smr[sbin*IMAX*RMAX + (rbin+1)*IMAX + ibin];
  TGraph *cf011 = th_oppo_smr[sbin*IMAX*RMAX + (rbin+1)*IMAX + ibin+1];
  TGraph *cf100 = th_oppo_smr[(sbin+1)*IMAX*RMAX + rbin*IMAX + ibin];
  TGraph *cf101 = th_oppo_smr[(sbin+1)*IMAX*RMAX + rbin*IMAX + ibin+1];
  TGraph *cf110 = th_oppo_smr[(sbin+1)*IMAX*RMAX + (rbin+1)*IMAX + ibin];
  TGraph *cf111 = th_oppo_smr[(sbin+1)*IMAX*RMAX + (rbin+1)*IMAX + ibin+1];

  double wgt000, wgt001, wgt010, wgt011, wgt100, wgt101, wgt110, wgt111;//weights to (*), the product of the differences between "siz rev and imv" und the values appropriate to bins which we selected above (**)
  double intv;//intv -- the sum of coefficients multiplied by values of the theoretical function at a given point of Y and in the end divided by denval, 
  double denval; //denval -- denominator value, the product of the differences between neighbour values of the Re_scatlength, Im_scatlength and source arrays for bins which we selected above (sbin,ibin,rbin).
  denval = ((source[sbin+1] - source[sbin])*(Re_scatlength[rbin+1] - Re_scatlength[rbin])*(Im_scatlength[ibin+1]-Im_scatlength[ibin]));
  
  wgt111 = (siz  - source[sbin])  *(rev  - Re_scatlength[rbin])  *(imv  - Im_scatlength[ibin])  ;
  wgt110 = (siz  - source[sbin])  *(rev  - Re_scatlength[rbin])  *(-imv + Im_scatlength[ibin+1]);
  wgt101 = (siz  - source[sbin])  *(-rev + Re_scatlength[rbin+1])*(imv  - Im_scatlength[ibin])  ;
  wgt100 = (siz  - source[sbin])  *(-rev + Re_scatlength[rbin+1])*(-imv + Im_scatlength[ibin+1]);
  wgt011 = (-siz + source[sbin+1])*(rev  - Re_scatlength[rbin])  *(imv  - Im_scatlength[ibin])  ;
  wgt010 = (-siz + source[sbin+1])*(rev  - Re_scatlength[rbin])  *(-imv + Im_scatlength[ibin+1]);
  wgt001 = (-siz + source[sbin+1])*(-rev + Re_scatlength[rbin+1])*(imv  - Im_scatlength[ibin])  ;
  wgt000 = (-siz + source[sbin+1])*(-rev + Re_scatlength[rbin+1])*(-imv + Im_scatlength[ibin+1]);

//try to get the function at any point of space
  for (int ib=0; ib<cf000->GetN(); ib++) 
    {
      intv = (wgt000 * cf000->GetY()[ib] +
	      wgt001 * cf001->GetY()[ib] +
	      wgt010 * cf010->GetY()[ib] +
	      wgt011 * cf011->GetY()[ib] +
	      wgt100 * cf100->GetY()[ib] +
	      wgt101 * cf101->GetY()[ib] +
	      wgt110 * cf110->GetY()[ib] +
	      wgt111 * cf111->GetY()[ib])/denval;
      vals[ib] = intv;
    }
}

//Similar to the above but this function is to same sign pairs and this is a 1D interpolation. We assume known interaction so only R parameter space is used 
void getlinexth_same(double siz, double *vals)
{
  int sbin=0;

  if (siz<source[0]) sbin = 0;
  else if (siz>source[SMAX-1]) sbin=SMAX-2;
  else {
    for (sbin=0; sbin<SMAX; sbin++)
      if (siz < source[sbin]) {sbin--; break;}
  }

  TGraph *cf0 = th_same_smr[sbin];
  TGraph *cf1 = th_same_smr[sbin+1];
  
  double wgt0, wgt1;
  double intv;

  double denval = (source[sbin+1] - source[sbin]);
  wgt1 = (siz  - source[sbin]);
  wgt0 = (-siz + source[sbin+1]);
  
  for (int ib=0; ib<cf0->GetN(); ib++) 
    {
      intv = (wgt0 * cf0->GetY()[ib] +
	      wgt1 * cf1->GetY()[ib])/denval;
      vals[ib] = intv;
    }
}

//chi2 between experimental and theoretical function, taking into account background (Gaussian) and normalization
//lambda - Gausian background magnitude
//radbg - Gaussian background width, variance
double getchi2(TH1D *data, double *model, double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;//thv- theoretical value, bgv- background value, ksv -- k* value

  //fit range
  int maxbin = 35;
  int minbin = 1;

  for (int ib=minbin; ib<=maxbin; ib++)
    {

      if (ib == 46) continue;
      if (ib == 47) continue;
      if (ib == 48) continue;
      if (ib == 49) continue;
      if (ib == 50) continue;

      val = data->GetBinContent(ib);
      ksv = data->GetXaxis()->GetBinCenter(ib);
      thv = model[ib-1];
      bgv = norm + lambda*TMath::Exp(-ksv*ksv*radbg*radbg);     
      val += -bgv + 1.0;
      val = (val-1.0)/purity + 1.0;
      err = data->GetBinError(ib);
      chi = (val-thv)*(val-thv)/(err*err);
      sumchi += chi;
    }

  cout << "full chi2 " << sumchi << "     ----       " << norm << " "<< lambda << " " << radbg << " " << purity << endl;

  return sumchi;
}


//fit function
void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t Pv  = par[0];//purity
  Double_t Rev = par[1];
  Double_t Imv = par[2];

  Double_t Lvs[7];//bg height
  Double_t Rvs[7];//bg width
  Double_t Svs[7];//size

  Double_t Nvs[56];//normalization, different for each exp. function

  for (int iter=0; iter<7; iter++) {
    Lvs[iter] = par[3+iter*3+0];
    Rvs[iter] = par[3+iter*3+1];
    Svs[iter] = par[3+iter*3+2];
  }

  for (int iter=0; iter<56; iter++) {
    Nvs[iter] = par[24+iter];
  }

  Double_t bgscales[7];
  
  for (int iter=0; iter<7; iter++) {
    bgscales[iter] = par[80+iter];
  }
  
  Double_t chi2 = 0.0;
  Double_t tv, qv;   
  Double_t vals[500];//interpolated theory

  //loop over centralities
  for (int icent=0; icent<7; icent++) {
    if (docent[icent]) {

      //theoretical interpolation for a given fit point
      //background is assumed the same for each pair
      getlinexth_oppo(Rev, Imv, Svs[icent], vals);

      chi2 += getchi2(cdatamfmp[icent], vals, Nvs[icent*4+0], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatamfpm[icent], vals, Nvs[icent*4+1], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatapfmp[icent], vals, Nvs[icent*4+2], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatapfpm[icent], vals, Nvs[icent*4+3], Pv, Lvs[icent], Rvs[icent]);

      //the only difference is additional scaling of the background magnitude bgscales
      getlinexth_same(Svs[icent], vals);
      
      chi2 += getchi2(cdatamfpp[icent], vals, Nvs[28+icent*4+0], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatamfmm[icent], vals, Nvs[28+icent*4+1], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatapfpp[icent], vals, Nvs[28+icent*4+2], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatapfmm[icent], vals, Nvs[28+icent*4+3], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
    }
  }

  f = chi2; //chi2 for all 56 functions
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////



int main(int argc, char **argv)
{
  // Declare functions for interpolation
  
  TFile *th_oppo_files[IMAX*RMAX*SMAX];
  TGraph *cf_th_oppo[IMAX*RMAX*SMAX];

  TFile *th_same_files[RMAX];
  TGraph *cf_th_same[RMAX];
  
  // Get baseline file
  

  Double_t val;

  // Load functions for interpolation
  // And divide by baseline
  
  for (int is=0; is<SMAX; is++) {
    for (int ir=0; ir<RMAX; ir++) {
      for (int ii=0; ii<IMAX; ii++) {
	
	th_oppo_files[is*RMAX*IMAX + ir*IMAX + ii] = new TFile (Form("./theoretical_cf/kpcoulstrongwgtlhc.rad%.1f.re%.3f.im%.3f.root",source[is], Re_scatlength[ir], Im_scatlength[ii]));

	cf_th_oppo[is*RMAX*IMAX + ir*IMAX + ii] = (TGraph *) th_oppo_files[is*RMAX*IMAX + ir*IMAX + ii]->Get("grkaonproton");
	
	if (cf_th_oppo[is*RMAX*IMAX + ir*IMAX + ii]) {
	  th_oppo[is*RMAX*IMAX + ir*IMAX + ii] = new TGraph(*cf_th_oppo[is*RMAX*IMAX + ir*IMAX + ii]);
	  th_oppo_smr[is*RMAX*IMAX + ir*IMAX + ii] = smearmr(th_oppo[is*RMAX*IMAX + ir*IMAX + ii]);
	}
	else {
	  th_oppo[is*RMAX*IMAX + ir*IMAX + ii] = 0;
	  th_oppo_smr[is*RMAX*IMAX + ir*IMAX + ii] = 0;
	}
      }
    }

    th_same_files[is] = new TFile (Form("./theoretical_cf/kpsscoulstrongwgtlhc.rad%.1f.re-0.360.im0.000.root",source[is]));    
    cf_th_same[is] = (TGraph *) th_same_files[is]->Get("grkaonproton");

    
    if (cf_th_same[is]) {
      th_same[is] = new TGraph(*cf_th_same[is]);
      th_same_smr[is] = smearmr(th_same[is]);
    }
    else {
      th_same[is] = 0;
      th_same_smr[is] = 0;
    }
    cout << "Got same sign graphs " << cf_th_same[is] << " " << th_same_smr[is] << endl;

  }

  // New Data from Wioleta from 24.04.2019 TPC-only tracks
  
  TFile *indatafilemfmp = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilemfpm = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilepfmp = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  TFile *indatafilepfpm = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  
  // Centrality 0
  cdatamfmp[0] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatamfpm[0] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdatapfmp[0] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatapfpm[0] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 1
  cdatamfmp[1] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM1");
  cdatamfpm[1] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM1");
  cdatapfmp[1] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM1");
  cdatapfpm[1] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM1");

  // Centrality 2
  cdatamfmp[2] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM2");
  cdatamfpm[2] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM2");
  cdatapfmp[2] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM2");
  cdatapfpm[2] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM2");

  // Centrality 3
  cdatamfmp[3] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM3");
  cdatamfpm[3] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM3");
  cdatapfmp[3] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM3");
  cdatapfpm[3] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM3");

  // Centrality 4
  cdatamfmp[4] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM4");
  cdatamfpm[4] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM4");
  cdatapfmp[4] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM4");
  cdatapfpm[4] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM4");

  // Centrality 5
  cdatamfmp[5] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM5");
  cdatamfpm[5] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM5");
  cdatapfmp[5] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM5");
  cdatapfpm[5] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM5");

  // Centrality 6
  cdatamfmp[6] = (TH1D *) indatafilemfmp->Get("CfnReYlm00cylmKmProtpcM6");
  cdatamfpm[6] = (TH1D *) indatafilemfpm->Get("CfnReYlm00cylmKpAProtpcM6");
  cdatapfmp[6] = (TH1D *) indatafilepfmp->Get("CfnReYlm00cylmKmProtpcM6");
  cdatapfpm[6] = (TH1D *) indatafilepfpm->Get("CfnReYlm00cylmKpAProtpcM6");

  // Data on same-sign correlations
  
  TFile *indatafilemfpp = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilemfmm = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilepfpp = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  TFile *indatafilepfmm = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  
  // Centrality 0
  cdatamfpp[0] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatamfmm[0] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdatapfpp[0] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatapfmm[0] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 1
  cdatamfpp[1] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM1");
  cdatamfmm[1] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM1");
  cdatapfpp[1] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM1");
  cdatapfmm[1] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM1");

  // Centrality 2
  cdatamfpp[2] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM2");
  cdatamfmm[2] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM2");
  cdatapfpp[2] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM2");
  cdatapfmm[2] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM2");

  // Centrality 3
  cdatamfpp[3] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM3");
  cdatamfmm[3] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM3");
  cdatapfpp[3] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM3");
  cdatapfmm[3] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM3");

  // Centrality 4
  cdatamfpp[4] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM4");
  cdatamfmm[4] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM4");
  cdatapfpp[4] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM4");
  cdatapfmm[4] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM4");

  // Centrality 5
  cdatamfpp[5] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM5");
  cdatamfmm[5] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM5");
  cdatapfpp[5] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM5");
  cdatapfmm[5] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM5");

  // Centrality 6
  cdatamfpp[6] = (TH1D *) indatafilemfpp->Get("CfnReYlm00cylmKpProtpcM6");
  cdatamfmm[6] = (TH1D *) indatafilemfmm->Get("CfnReYlm00cylmKmAProtpcM6");
  cdatapfpp[6] = (TH1D *) indatafilepfpp->Get("CfnReYlm00cylmKpProtpcM6");
  cdatapfmm[6] = (TH1D *) indatafilepfmm->Get("CfnReYlm00cylmKmAProtpcM6");

  Double_t xmax = cdatamfmp[0]->GetXaxis()->GetXmax();

  double purval = 0.75;

//"normalizacja + lamda*exp(-x*x*width*width)"
  double normval_oppo[7];
  double lambdaval_oppo[7];
  double widthval_oppo[7];

  double normval_oppo_fix[28];  

  TF1 *fnorm = new TF1("fnorm","[0]");//function used for normalization fitting, high k*
  TF1 *fgauss = new TF1("fgauss","[0] + [1]*exp(-x*x*[2]*[2])"); //function used for Gaussian fitting, medium k*


  //1st step - fitting of background only to fix normalization and background parameters
  //opposite sign fitting


  for (int icent=0; icent<7; icent++) {
    normval_oppo[icent]=0;
    lambdaval_oppo[icent]=0;
    widthval_oppo[icent]=0;
    for (int ipair=0; ipair<4; ipair++) {
       if (docent[icent]) {

      cout << "Fitting opposite sign background " << icent << " " << ipair << endl;
      
      fnorm->SetParameter(0,1.000);
      cdatamfmp[icent]->Fit(fnorm, "", "", 0.85*xmax, xmax); //k* fit region for normalization - very high k* interval
      
      fgauss->SetParameters(fnorm->GetParameter(0),0.002,2.0);
      fgauss->FixParameter(0,fnorm->GetParameter(0));
      //  fgauss->FixParameter(2,2.0);
      cdatamfmp[icent]->Fit(fgauss, "","",0.35, xmax); //k* fit of background slope, medium k* region
      normval_oppo_fix[icent*4+ipair] = fgauss->GetParameter(0);
      cout << "tu" << icent << " " << ipair << endl;
      //saving background fits
      normval_oppo[icent] += fgauss->GetParameter(0)/4;
      lambdaval_oppo[icent] += fgauss->GetParameter(1)/4;
      widthval_oppo[icent] += fgauss->GetParameter(2)/4;
      cout << "sie cos dzieje " << icent << " " << ipair << endl;
      }
    }
  }

  
  double normval_same[7];
  double lambdaval_same[7];
  double widthval_same[7];
  
  double normval_same_fix[28];

  //same sign background fitting
  for (int icent=0; icent<7; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
       if (docent[icent]) {

     	  cout << "Fitting same sign background " << icent << " " << ipair << endl;
    
    	  fnorm->SetParameter(0,1.000);
     	  cdatamfpp[icent]->Fit(fnorm, "", "", 0.85*xmax, xmax);
      
    	  fgauss->SetParameters(fnorm->GetParameter(0),0.002,2.0);
     	  fgauss->FixParameter(0,fnorm->GetParameter(0));
      	  //fgauss->FixParameter(2,2.0);
          cdatamfpp[icent]->Fit(fgauss, "","",0.3, xmax);
          normval_same_fix[icent*4+ipair] = fnorm->GetParameter(0);

          //saving background fits
          normval_same[icent] += fgauss->GetParameter(0)/4;
          lambdaval_same[icent] += fgauss->GetParameter(1)/4;
          widthval_same[icent] += fgauss->GetParameter(2)/4;
       }
    }
  }
 
  
  TGraph2D *grchi = new TGraph2D();
  grchi->SetName("grchi");
  grchi->SetTitle("#chi^{2} Map;#Rgothicf_{0} (fm);#Jgothicf_{0} (fm);#chi^{2}");

  double chival;
  int pcount = 0;

  double minchi = 1e9;

  int miniChi_Re, miniChi_Im;
  int miniChi_Re_sum=0, miniChi_Im_sum=0, centcount=0;
  int miniChi_Re_s[7], miniChi_Im_s[7], miniChi_size[7];

  for (int icent=0; icent<7; icent++) {
    miniChi_Re_s[icent] = 1;
    miniChi_Im_s[icent] = 1;
    miniChi_size[icent] = 1;
  }
  
  // Create a fit chi2 map

  //2nd step, we calculate chi2 between precomputed theoretical function (not interpolated but MR smeared) and experimental ones, background and normalization assumed from the previous fit
  //this is to calculate starting values for final fitting
  for (int icent=0; icent<7; icent++) {
    if (docent[icent]) {
      minchi = 1e9;
      for (int is=0; is<SMAX; is++) {
	for (int ir=0; ir<RMAX; ir++) {
	  for (int ii=0; ii<IMAX; ii++) {
	    if (th_oppo_smr[is*RMAX*IMAX + ir*IMAX + ii]) {
	      chival  = getchi2(cdatamfmp[icent], th_oppo_smr[is*RMAX*IMAX + ir*IMAX + ii]->GetY(), normval_oppo[icent], purval, lambdaval_oppo[icent], widthval_oppo[icent]);
	      // grchi->SetPoint(pcount++, Re_scatlength[ir], Im_scatlength[ii], chival);	
	      if (chival < minchi) 
		{ minchi = chival;
		  miniChi_Re_s[icent] = ir;//the bin number of Re scat.lenght with minimum chival for a given size of the source and centrality
		  miniChi_Im_s[icent] = ii;//the bin number of Im scat.lenght with minimum chival for a given size of the source and centrality
		  miniChi_size[icent] = is;//the bin number of source size with minimum chival for a given centrality
		  cout << "Now minimum is " << minchi << "    s r i " << is << " "<< ir << " "<< ii << endl;
		}
	    }
	  }
	}
	miniChi_Re_sum += miniChi_Re_s[icent];
	miniChi_Im_sum += miniChi_Im_s[icent];
	centcount++;
      }
    }
  }

  miniChi_Re= miniChi_Re_sum/centcount;
  miniChi_Im = miniChi_Im_sum/centcount;
  

  TFile *outfile = new TFile("outchicentallwgt.root","RECREATE");
  outfile->cd();  
  
  cout << "Minimum chi " << minchi << endl;
  cout << "Re Im Size " << Re_scatlength[miniChi_Re] << " " << Im_scatlength[miniChi_Im] << " ";

  for (int icent=0; icent<7; icent++)
    if (docent[icent])
      cout << "Centrality " << icent << "  " << source[miniChi_size[icent]] << endl;

  //if the point on the edge of 3D lattice, we move it inside
  if (miniChi_Re== RMAX-1) miniChi_Re--;
  if (miniChi_Im == IMAX-1) miniChi_Im--;
  
  for (int icent=0; icent<7; icent++)
    if (miniChi_size[icent] == SMAX-1) miniChi_size[icent]--;

  double vals[500];

  cout << "Now fitting " << endl;
  
  // Values from Therminator2 fits
  // lbdval = 0.00555;
  // wdtval = 1.92;
  
  cout << "Setting values " << normval_oppo[0] << " " << lambdaval_oppo[0] << " " << widthval_oppo[0] << " " << purval << endl;
  
  TVirtualFitter *fitter=TVirtualFitter::Fitter(0,87);
  fitter->SetFCN(myfuncf);
  fitter->SetParameter(0,  "Purity"  ,purval        ,0.0001,0.2,   0.85);
  fitter->SetParameter(1,  "Ref0"    ,Re_scatlength[miniChi_Re]  ,0.0001,-5.0, 5.0);
  fitter->SetParameter(2,  "Imf0"    ,Im_scatlength[miniChi_Im]  ,0.0001,0.00001, 5.0);
  // fitter->SetParameter(2,  "Imf0"    ,0.2 ,0.0001,0.00001, 5.0);
  // fitter->FixParameter(2);
  
  
  for (int icent=0; icent<7; icent++) {
    fitter->SetParameter(3+icent*3+0,  Form("Lambda%i", icent)  ,TMath::Abs(lambdaval_same[icent])   ,0.001, 0.0,   0.1);
    fitter->SetParameter(3+icent*3+1,  Form("Width%i", icent)   ,TMath::Abs(widthval_same[icent])   ,0.001, 1.0,   10.0);
    fitter->SetParameter(3+icent*3+2,  Form("Radius%i", icent)  ,source[miniChi_size[icent]]  ,0.01  ,0.9, 8.0);

    if (!docent[icent]) {
      fitter->FixParameter(3+icent*3+0);
      fitter->FixParameter(3+icent*3+1);
      fitter->FixParameter(3+icent*3+2);
    }

    //fitter->FixParameter(3+icent*3+0);
    fitter->FixParameter(3+icent*3+1);
  }

  //normalizations, set from background fit
  for (int icent=0; icent<7; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //      fitter->SetParameter(21+icent*4+ipair,  Form("Norm%i%i", icent, ipair)    ,normval_oppo[icent]    ,0.001, 0.9,   1.1);
      fitter->SetParameter(24+icent*4+ipair,  Form("Norm%i%i", icent, ipair)    ,normval_oppo_fix[icent*4+ipair]    ,0.001, 0.9,   1.1);  
      if (!docent[icent])
         fitter->FixParameter(24+icent*4+ipair);
      fitter->FixParameter(24+icent*4+ipair);
    }
  }
  
  for (int icent=0; icent<7; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //      fitter->SetParameter(45+icent*4+ipair,  Form("Normss%i%i", icent, ipair)    ,normval_same[icent]    ,0.001, 0.9,   1.1);
      fitter->SetParameter(52+icent*4+ipair,  Form("Normss%i%i", icent, ipair)    ,normval_same_fix[icent*4+ipair]    ,0.001, 0.9,   1.1);      
      if (!docent[icent])
	fitter->FixParameter(52+icent*4+ipair);
      fitter->FixParameter(52+icent*4+ipair);
    }
  }

  //fixing of background width and height
  //  fitter->SetParameter(69, "BGsstoul", lambdaval_same[0]/lambdaval_oppo[0], 0.001, 0.1, 2.0);
  for (int icent=0; icent<7; icent++) 
    fitter->SetParameter(80+icent, Form("BGsstoul%i", icent), 0.6, 0.001, 0.1, 5.0);
  
  // fitter->FixParameter(69+icent);
  
  // fitter->SetParameter(0,  "Norm"    ,normval_oppo[0]    ,0.001, 0.9,   1.1);
  // fitter->SetParameter(7,  "nrmsc1"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(8,  "nrmsc2"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(9,  "nrmsc3"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(10, "Norm3"   ,normval_oppo[3]    ,0.001, 0.9,   1.1);
  // fitter->SetParameter(11, "Lambda3" ,lambdaval_oppo[3]    ,0.001, 0.0,   0.1);
  // fitter->SetParameter(12, "Width3"  ,widthval_oppo[3]    ,0.001, 1.0,   10.0);
  // fitter->SetParameter(13, "Radius3" ,source[minsz3] ,0.01  ,0.9, 8.0);
  // fitter->SetParameter(14, "nrmsc4"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(15, "nrmsc5"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(16, "nrmsc6"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(5, "Imf0"    ,0.25  ,0.001,0.0, 5.0);
  // fitter->SetParameter(6, "norm"     ,1.0       ,0.0001,0.0,   2.0);


  
  Double_t arglist[100];
  arglist[0] = 1;
  fitter->ExecuteCommand("CALL FCN", arglist, 1);
  fitter->FixParameter(0);

  //  fitter->FixParameter(1);
  //  fitter->FixParameter(2);
  //  fitter->FixParameter(3);
  //  fitter->FixParameter(4);
  //  fitter->FixParameter(5);
  //  fitter->FixParameter(6);
  //  fitter->FixParameter(3);
  //  fitter->FixParameter(11);
  //  fitter->FixParameter(14);
  //  fitter->FixParameter(3);
  //  fitter->FixParameter(4);
  //  fitter->FixParameter(15);
  //  fitter->FixParameter(7);
  //  fitter->FixParameter(8);
  //  fitter->FixParameter(9);
  //  fitter->FixParameter(10);
  //  fitter->FixParameter(12);
  //  fitter->FixParameter(13);
  //  fitter->FixParameter(14);
  //  fitter->SetParameter(5,0.25);
  //  fitter->FixParameter(5);
  //  fitter->FixParameter(15);

  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT", arglist, 1);
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
  
//   f1 ->SetParameters(fitter->GetParameter(0), fitter->GetParameter(1), fitter->GetParameter(2), 1.0, n0);
//   f20->SetParameters(Ro, Rs, Rl, 1.0, n0*n1);
//   f22->SetParameters(Ro, Rs, Rl, 1.0, n0*n2);

  cout << "Fitting done " << " " << fitter->GetParameter(4) << " " << fitter->GetParameter(5) << "   " << fitter->GetParameter(6) << endl;

  // double fitnorm = fitter->GetParameter(0);
  // double fitlabd = fitter->GetParameter(1);
  // double fitwdth = fitter->GetParameter(2);

  double fit_puri = fitter->GetParameter(0); 
  double fit_Ref0 = fitter->GetParameter(1);
  double fit_Imf0 = fitter->GetParameter(2);

  // double fitrad0 = fitter->GetParameter(6);
  // double fitnrm1 = fitter->GetParameter(7);
  // double fitnrm2 = fitter->GetParameter(8);
  // double fitnrm3 = fitter->GetParameter(9);
  // double fitnorm3 = fitter->GetParameter(10);
  // double fitlabd3 = fitter->GetParameter(11);
  // double fitwdth3 = fitter->GetParameter(12);
  // double fitrad3 = fitter->GetParameter(13);
  // double fitnrm4 = fitter->GetParameter(14);
  // double fitnrm5 = fitter->GetParameter(15);
  // double fitnrm6 = fitter->GetParameter(16);

  //  double fitnorms[6];
  double fit_lambda[7];
  double fit_width[7];
  double fit_radius[7];

  double fit_normScales_oppo[4*7];
  double fit_normScales_same[4*7];

  for (int icent=0; icent<7; icent++) {
    fit_lambda[icent] = fitter->GetParameter(3+icent*3+0);
    fit_width[icent] = fitter->GetParameter(3+icent*3+1);
    fit_radius[icent] = fitter->GetParameter(3+icent*3+2);

    for (int ipair=0; ipair<4; ipair++) {
      fit_normScales_oppo[icent*4+ipair] = fitter->GetParameter(24+icent*4+ipair);
      fit_normScales_same[icent*4+ipair] = fitter->GetParameter(52+icent*4+ipair);
    }
  }

  double fit_baScales[7];

  for (int icent=0; icent<7; icent++) {
    fit_baScales[icent] = fitter->GetParameter(80+icent);
  }
  
  double fiterr_puri = fitter->GetParError(0); 
  double fiterr_Ref0 = fitter->GetParError(1);
  double fiterr_Imf0 = fitter->GetParError(2);

  // double fitrad0 = fitter->GetParError(6);
  // double fitnrm1 = fitter->GetParError(7);
  // double fitnrm2 = fitter->GetParError(8);
  // double fitnrm3 = fitter->GetParError(9);
  // double fitnorm3 = fitter->GetParError(10);
  // double fitlabd3 = fitter->GetParError(11);
  // double fitwdth3 = fitter->GetParError(12);
  // double fitrad3 = fitter->GetParError(13);
  // double fitnrm4 = fitter->GetParError(14);
  // double fitnrm5 = fitter->GetParError(15);
  // double fitnrm6 = fitter->GetParError(16);

  // double fitnorms[6];
  double fiterr_lambda[7];
  double fiterr_wdth[7];
  double fiterr_radius[7];

  double fiterr_normScales_oppo[4*7];//normalization errors for opposite sign
  double fiterr_normScales_same[4*7];

  for (int icent=0; icent<7; icent++) {
    fiterr_lambda[icent] = fitter->GetParError(3+icent*3+0);
    fiterr_wdth[icent] = fitter->GetParError(3+icent*3+1);
    fiterr_radius[icent] = fitter->GetParError(3+icent*3+2);

    for (int ipair=0; ipair<4; ipair++) {
      fiterr_normScales_oppo[icent*4+ipair] = fitter->GetParError(24+icent*4+ipair);
      fiterr_normScales_same[icent*4+ipair] = fitter->GetParError(52+icent*4+ipair);
    }
  }

  double fiterr_bgScales[7];

  for (int icent=0; icent<7; icent++) {
    fiterr_bgScales[icent] = fitter->GetParError(80+icent);
  }
    
  cout << "Final paramerers:" << endl;
  cout << "Purity          : " << fit_puri << " +/- " << fiterr_puri <<  endl;
  cout << "Re f0           : " << fit_Ref0 << " +/- " << fiterr_Ref0 <<  endl;
  cout << "Im f0           : " << fit_Imf0 << " +/- " << fiterr_Imf0 <<  endl;

  for (int icent=0; icent<7; icent++) {
    cout << "Radius C" << icent << "     :" << fit_radius[icent]  << " +/- " << fiterr_radius[icent] << endl;
  }

  for (int icent=0; icent<7; icent++) {
    cout << "Backgroung parameters centrality " << icent << endl;
    cout << "alfa width " << fit_lambda[icent] << " " << fit_width[icent] << endl;
    cout << "Normalizations :"
	 << fit_normScales_oppo[icent*4+0] << " "
	 << fit_normScales_oppo[icent*4+1] << " "
	 << fit_normScales_oppo[icent*4+2] << " "
	 << fit_normScales_oppo[icent*4+3] << " "
	 << endl;
    cout << "Normalizations same-sign:"
	 << fit_normScales_same[icent*4+0] << " "
	 << fit_normScales_same[icent*4+1] << " "
	 << fit_normScales_same[icent*4+2] << " "
	 << fit_normScales_same[icent*4+3] << " "
	 << endl;

  }

  for (int icent=0; icent<7; icent++) {
    cout << "Background ss to ul scale " << icent << " : " << fit_baScales[icent] << " +/- " << fiterr_bgScales[icent] << endl;
  }
  
  outfile->cd();

  TGraph *grcfits_oppo[7];
  TGraph *grcfits_same[7];
  
  for (int icent=0; icent<7; icent++) {
    if (docent[icent]) {
      
      getlinexth_oppo(fit_Ref0, fit_Imf0, fit_radius[icent], vals);    
      cout << "First bins " << vals[0] << " " << vals[1] << " " << vals[2] << endl;
      grcfits_oppo[icent] = new TGraph(200, th_oppo_smr[0]->GetX(), vals);
      cout << "Have graph " << grcfits_oppo[icent] << " with " << grcfits_oppo[icent]->GetN() << " bins " << endl;
      grcfits_oppo[icent]->SetTitle(";k^ (GeV/#it{c});C(k*)");
      grcfits_oppo[icent]->SetName(Form("graulfitC%i", icent));
      grcfits_oppo[icent]->Write();

      cout << "Doing linear extrapolation " << endl;
      
      getlinexth_same(fit_radius[icent], vals);
      cout << "First bins " << vals[0] << " " << vals[1] << " " << vals[2] << endl;    
      grcfits_same[icent] = new TGraph(200, th_same_smr[0]->GetX(), vals);
      cout << "Have graph " << grcfits_same[icent] << " with " << grcfits_same[icent]->GetN() << " bins " << endl;
      grcfits_same[icent]->SetTitle(";k^ (GeV/#it{c});C(k*)");
      grcfits_same[icent]->SetName(Form("grassfitC%i", icent));
      grcfits_same[icent]->Write();

    }
  }
  

  TH1D *cfitmfpmlin[7];
  TH1D *cfitmfmplin[7];
  TH1D *cfitpfpmlin[7];
  TH1D *cfitpfmplin[7];

  TH1D *cfitmfpplin[7];
  TH1D *cfitmfmmlin[7];
  TH1D *cfitpfpplin[7];
  TH1D *cfitpfmmlin[7];

  double ksv, bgv;//k* and background values 
  
  for (int icent=0; icent<7; icent++) {
    if (docent[icent]) {
    
      cfitmfpmlin[icent] = new TH1D(*cdatamfpm[icent]);
      cfitmfpmlin[icent]->SetName(Form("datamfpmbgC%i",icent));
      cfitmfmplin[icent] = new TH1D(*cdatamfmp[icent]);
      cfitmfmplin[icent]->SetName(Form("datamfmpbgC%i",icent));
      cfitpfpmlin[icent] = new TH1D(*cdatapfpm[icent]);
      cfitpfpmlin[icent]->SetName(Form("datapfpmbgC%i",icent));
      cfitpfmplin[icent] = new TH1D(*cdatapfmp[icent]);
      cfitpfmplin[icent]->SetName(Form("datapfmpbgC%i",icent));
      
      for (int ib=1; ib<=cdatamfmp[0]->GetNbinsX(); ib++) {
	ksv = cdatamfmp[icent]->GetXaxis()->GetBinCenter(ib);
	bgv = fit_normScales_oppo[icent*4+0] + fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
      
	val = cdatamfmp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitmfmplin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_oppo[icent*4+1] + fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatamfpm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitmfpmlin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_oppo[icent*4+2] + fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatapfmp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitpfmplin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_oppo[icent*4+3] + fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatapfpm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitpfpmlin[icent]->SetBinContent(ib, val);
      }
    
      cfitmfmplin[icent]->Write();
      cfitmfpmlin[icent]->Write();
      cfitpfmplin[icent]->Write();
      cfitpfpmlin[icent]->Write();

      cfitmfpplin[icent] = new TH1D(*cdatamfpp[icent]);
      cfitmfpplin[icent]->SetName(Form("datamfppbgC%i",icent));
      cfitmfmmlin[icent] = new TH1D(*cdatamfmm[icent]);
      cfitmfmmlin[icent]->SetName(Form("datamfmmbgC%i",icent));
      cfitpfpplin[icent] = new TH1D(*cdatapfpp[icent]);
      cfitpfpplin[icent]->SetName(Form("datapfppbgC%i",icent));
      cfitpfmmlin[icent] = new TH1D(*cdatapfmm[icent]);
      cfitpfmmlin[icent]->SetName(Form("datapfmmbgC%i",icent));
      
      for (int ib=1; ib<=cdatamfmm[0]->GetNbinsX(); ib++) {
	ksv = cdatamfpp[icent]->GetXaxis()->GetBinCenter(ib);
	bgv = fit_normScales_same[icent*4+0] + fit_baScales[icent]*fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
      
	val = cdatamfpp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitmfpplin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_same[icent*4+1] + fit_baScales[icent]*fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatamfmm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitmfmmlin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_same[icent*4+2] + fit_baScales[icent]*fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatapfpp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitpfpplin[icent]->SetBinContent(ib, val);
	
	bgv = fit_normScales_same[icent*4+3] + fit_baScales[icent]*fit_lambda[icent]*TMath::Exp(-ksv*ksv*fit_width[icent]*fit_width[icent]);
	
	val = cdatapfmm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fit_puri + 1.0;
	
	cfitpfmmlin[icent]->SetBinContent(ib, val);
      }
    
      cfitmfmmlin[icent]->Write();
      cfitmfpplin[icent]->Write();
      cfitpfmmlin[icent]->Write();
      cfitpfpplin[icent]->Write();

    }
  }

  Double_t fvx[87];
  Double_t fvy[87];
  Double_t fve[87];

  for (int ipar=0; ipar<88; ipar++) {
    fvx[ipar] = ipar;
    fvy[ipar] = fitter->GetParameter(ipar);
    fve[ipar] = fitter->GetParError(ipar);
  }

  TGraphErrors *fitresult = new TGraphErrors(87, fvx, fvy, NULL, fve);
  fitresult->SetName("fitresult");
  fitresult->Write();
  
  
}
 
