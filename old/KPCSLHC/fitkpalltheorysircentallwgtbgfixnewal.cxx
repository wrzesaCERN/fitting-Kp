#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TGraph2D.h>

//#define RMAX 8
#define RMAX 9
#define IMAX 10
//#define SMAX 10
#define SMAX 8

using namespace std;

Double_t rvals[RMAX] = { -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.001, 0.1, 0.2 };
//Double_t rvals[RMAX] = { -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.2 };
//Double_t rvals[RMAX] = { -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, 0.001, 0.5 };
Double_t ivals[IMAX] = { 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
//Double_t svals[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0 };
Double_t svals[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0 };
//Double_t svals[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5 };

Int_t docent[6] = { 1, 0, 0, 0, 0, 0 };

TGraph *cfs00sbgs[IMAX*RMAX*SMAX];
TGraph *cfs00sbgsmr[IMAX*RMAX*SMAX];
TH1D *cfs00ppbg;

TGraph *cfs00sssbgs[RMAX];
TGraph *cfs00sssbgsmr[RMAX];

TH1D *cfs00ppbgmr;

//TH1D *cdata;

TH1D *cdataffpm[6];
TH1D *cdataffmp[6];
TH1D *cdatarfpm[6];
TH1D *cdatarfmp[6];

TH1D *cdataffpp[6];
TH1D *cdataffmm[6];
TH1D *cdatarfpp[6];
TH1D *cdatarfmm[6];

TH1D *cdatapp;

TGraph *smearmr(TGraph *cfppdv)
{
  double width = 0.008;
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


void getlinext(double rev, double imv, double siz, double *vals)
{
  int rbin=0, ibin=0, sbin=0;

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

  if (siz<svals[0]) sbin = 0;
  else if (siz>svals[SMAX-1]) sbin=SMAX-2;
  else {
    for (sbin=0; sbin<SMAX; sbin++)
      if (siz < svals[sbin]) {sbin--; break;}
    //      else cout << " Next i" << endl;
      
  }

  // cout << "r i " << rbin << " " << ibin << endl;
  
  // cout << " " << rev << " " << rvals[rbin] << "  - "  << rvals[rbin+1] << endl;
  // cout << " " << imv << " " << ivals[ibin] << "  - "  << ivals[ibin+1] << endl;
  
  // TH1D *cf00 = cfs00sbgs[rbin*IMAX + ibin];
  // TH1D *cf01 = cfs00sbgs[rbin*IMAX + ibin+1];
  // TH1D *cf10 = cfs00sbgs[(rbin+1)*IMAX + ibin];
  // TH1D *cf11 = cfs00sbgs[(rbin+1)*IMAX + ibin+1];
  
  TGraph *cf000 = cfs00sbgsmr[sbin*IMAX*RMAX + rbin*IMAX + ibin];
  TGraph *cf001 = cfs00sbgsmr[sbin*IMAX*RMAX + rbin*IMAX + ibin+1];
  TGraph *cf010 = cfs00sbgsmr[sbin*IMAX*RMAX + (rbin+1)*IMAX + ibin];
  TGraph *cf011 = cfs00sbgsmr[sbin*IMAX*RMAX + (rbin+1)*IMAX + ibin+1];
  TGraph *cf100 = cfs00sbgsmr[(sbin+1)*IMAX*RMAX + rbin*IMAX + ibin];
  TGraph *cf101 = cfs00sbgsmr[(sbin+1)*IMAX*RMAX + rbin*IMAX + ibin+1];
  TGraph *cf110 = cfs00sbgsmr[(sbin+1)*IMAX*RMAX + (rbin+1)*IMAX + ibin];
  TGraph *cf111 = cfs00sbgsmr[(sbin+1)*IMAX*RMAX + (rbin+1)*IMAX + ibin+1];
  
  double wgt000, wgt001, wgt010, wgt011, wgt100, wgt101, wgt110, wgt111;
  double intv;

  double denval = ((svals[sbin+1] - svals[sbin])*(rvals[rbin+1] - rvals[rbin])*(ivals[ibin+1]-ivals[ibin]));
  
  wgt111 = (siz  - svals[sbin])  *(rev  - rvals[rbin])  *(imv  - ivals[ibin])  ;
  wgt110 = (siz  - svals[sbin])  *(rev  - rvals[rbin])  *(-imv + ivals[ibin+1]);
  wgt101 = (siz  - svals[sbin])  *(-rev + rvals[rbin+1])*(imv  - ivals[ibin])  ;
  wgt100 = (siz  - svals[sbin])  *(-rev + rvals[rbin+1])*(-imv + ivals[ibin+1]);
  wgt011 = (-siz + svals[sbin+1])*(rev  - rvals[rbin])  *(imv  - ivals[ibin])  ;
  wgt010 = (-siz + svals[sbin+1])*(rev  - rvals[rbin])  *(-imv + ivals[ibin+1]);
  wgt001 = (-siz + svals[sbin+1])*(-rev + rvals[rbin+1])*(imv  - ivals[ibin])  ;
  wgt000 = (-siz + svals[sbin+1])*(-rev + rvals[rbin+1])*(-imv + ivals[ibin+1]);

  // cout << "  wgt "
  //      << wgt000/denval << " " << wgt001/denval << " "
  //      << wgt010/denval << " " << wgt011/denval << " "
  //      << wgt100/denval << " " << wgt101/denval << " "
  //      << wgt110/denval << " " << wgt111/denval << endl;
  
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

      // cout << "int "
      // 	   << cf000->GetY()[ib] << " " << cf001->GetY()[ib] << " "
      // 	   << cf010->GetY()[ib] << " " << cf011->GetY()[ib] << " " 
      // 	   << cf100->GetY()[ib] << " " << cf101->GetY()[ib] << " "
      // 	   << cf110->GetY()[ib] << " " << cf111->GetY()[ib] << " " 
      // 	   << "     " << intv << endl;
      
      // cout << " --- " << sbin << " " << rbin << " " << ibin << endl;

    }
}

void getlinextss(double siz, double *vals)
{
  int sbin=0;

  if (siz<svals[0]) sbin = 0;
  else if (siz>svals[SMAX-1]) sbin=SMAX-2;
  else {
    for (sbin=0; sbin<SMAX; sbin++)
      if (siz < svals[sbin]) {sbin--; break;}
  }

  TGraph *cf0 = cfs00sssbgsmr[sbin];
  TGraph *cf1 = cfs00sssbgsmr[sbin+1];
  
  double wgt0, wgt1;
  double intv;

  double denval = (svals[sbin+1] - svals[sbin]);
  
  wgt1 = (siz  - svals[sbin]);
  wgt0 = (-siz + svals[sbin+1]);
  
  for (int ib=0; ib<cf0->GetN(); ib++) 
    {
      intv = (wgt0 * cf0->GetY()[ib] +
	      wgt1 * cf1->GetY()[ib])/denval;
      vals[ib] = intv;
    }
}

double getchi2(TH1D *data, double *model, double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;

  int maxbin = 35;
  int minbin = 3;

  //  int lobin, hibin;

  //  cout << "Fitting for n p l r " << norm << " " << purity << " " << lambda << " "<< radbg << endl;
  
  //  Double_t lokval, hikval;

  for (int ib=minbin; ib<=maxbin; ib++)
    {

      if (ib == 46) continue;
      if (ib == 47) continue;
      if (ib == 48) continue;
      if (ib == 49) continue;
      if (ib == 50) continue;

      val = data->GetBinContent(ib);
      ksv = data->GetXaxis()->GetBinCenter(ib);
      
      //      cout << "Finding bin " << endl;
      // lobin = (int) floor((ksv-0.001)/0.002);
      // hibin = lobin+1;

      // lokval = lobin*0.002+0.001;
      // hikval = hibin*0.002+0.001;

      // thv = (model[lobin] * (hikval-ksv)/0.002 +
      // 	     model[hibin] * (ksv-lokval)/0.002);

      thv = model[ib-1];
      
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

// double getchi2ss(TH1D *data, double *model, double norm, double purity, double lambda, double radbg)
// {
//   double sumchi = 0;
//   double val, err, chi, thv, bgv, ksv;

//   int maxbin = 55;
//   int minbin = 3;

//   int lobin, hibin;

//   Double_t lokval, hikval;

//   for (int ib=minbin; ib<=maxbin; ib++)
//     {
//       if (ib == 47) continue;
//       if (ib == 48) continue;
//       if (ib == 49) continue;

//       val = data->GetBinContent(ib);
//       ksv = data->GetXaxis()->GetBinCenter(ib);
      
//       lobin = (int) floor((ksv-0.001)/0.002);
//       hibin = lobin+1;

//       lokval = lobin*0.002+0.001;
//       hikval = hibin*0.002+0.001;

//       thv = (model[lobin] * (hikval-ksv)/0.002 +
// 	     model[hibin] * (ksv-lokval)/0.002);
            
//       bgv = norm + lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
//       val += -bgv + 1.0;
//       val = (val-1.0)/purity + 1.0;

//       err = data->GetBinError(ib);

//       chi = (val-thv)*(val-thv)/(err*err);
      
//       sumchi += chi;
//     }

//   cout << "full chi2 " << sumchi << "     ----       " << norm << " "<< lambda << " " << radbg << " " << purity << endl;
//   //  cout << "full chi2 " << sumchi << endl;

//   return sumchi;
// }

// double getchi2pp(double norm, double purity, double lambda, double radbg)
// {
//   double sumchi = 0;
//   double val, err, chi, thv, bgv, ksv;

//   int maxbin = 20;
//   int minbin = 3;

//   for (int ib=minbin; ib<=maxbin; ib++)
//     {
//       val = cdatapp->GetBinContent(ib);
//       ksv = cdatapp->GetXaxis()->GetBinCenter(ib);
//       bgv = norm + 1.4*lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
//       val += -bgv + 1.0;
//       val = (val-1.0)/purity + 1.0;

//       thv = cfs00ppbgmr->GetBinContent(ib);
//       err = cdatapp->GetBinError(ib);

//       chi = (val-thv)*(val-thv)/(err*err);
      
//       sumchi += chi;
//     }

//   return sumchi;

// }

void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t Pv  = par[0];
  Double_t Rev = par[1];
  Double_t Imv = par[2];

  Double_t Lvs[6];
  Double_t Rvs[6];
  Double_t Svs[6];

  Double_t Nvs[48];

  for (int iter=0; iter<6; iter++) {
    Lvs[iter] = par[3+iter*3+0];
    Rvs[iter] = par[3+iter*3+1];
    Svs[iter] = par[3+iter*3+2];
  }

  for (int iter=0; iter<48; iter++) {
    Nvs[iter] = par[21+iter];
  }

  Double_t bgscales[6];
  
  for (int iter=0; iter<6; iter++) {
    bgscales[iter] = par[69+iter];
  }
  //    Double_t bgscale = par[69];
  
  Double_t chi2 = 0.0;
  Double_t tv, qv;
    
  Double_t vals[500];

  //  cout << "Funcf called " << Nv << " " << Lv << " " << Rv << " " << Pv << " " << Rev << " " << Imv << " " << Szv << endl;

  for (int icent=0; icent<6; icent++) {
    if (docent[icent]) {
  
      getlinext(Rev, Imv, Svs[icent], vals);
      
      chi2 += getchi2(cdatarfmp[icent], vals, Nvs[icent*4+0], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatarfpm[icent], vals, Nvs[icent*4+1], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdataffmp[icent], vals, Nvs[icent*4+2], Pv, Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdataffpm[icent], vals, Nvs[icent*4+3], Pv, Lvs[icent], Rvs[icent]);

      getlinextss(Svs[icent], vals);

      cout << "Now chi2 for ss " << endl;
      
      chi2 += getchi2(cdatarfpp[icent], vals, Nvs[24+icent*4+0], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdatarfmm[icent], vals, Nvs[24+icent*4+1], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdataffpp[icent], vals, Nvs[24+icent*4+2], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
      chi2 += getchi2(cdataffmm[icent], vals, Nvs[24+icent*4+3], Pv, bgscales[icent]*Lvs[icent], Rvs[icent]);
    }
  }
      // chi2 += getchi2(cdatarfpm[0], vals, Nv*Ns1, Pv, Lv, Rv);
      // chi2 += getchi2(cdataffmp[0], vals, Nv*Ns2, Pv, Lv, Rv);
      // chi2 += getchi2(cdataffpm[0], vals, Nv*Ns3, Pv, Lv, Rv);

  // getlinext(Rev, Imv, Sz3, vals);

  // chi2 += getchi2(cdatarfmp[3], vals, Nv3,     Pv, Lv3, Rv3);
  // chi2 += getchi2(cdatarfpm[3], vals, Nv3*Ns4, Pv, Lv3, Rv3);
  // chi2 += getchi2(cdataffmp[3], vals, Nv3*Ns5, Pv, Lv3, Rv3);
  // chi2 += getchi2(cdataffpm[3], vals, Nv3*Ns6, Pv, Lv3, Rv3);

  //  cout << " chi2 now " << chi2 << endl;

  //  chi2 += getchi2pp(Nv, Pv, Lv, Rv);

  //  cout << " chi2 now " << chi2 << endl;

  //  cout << "Funcf called " << Nv << " " << Lv << " " << Rv << " " << Pv << " " << Rev << " " << Imv << " " << Szv << "    " << chi2 << endl;

  f = chi2;
}


int main(int argc, char **argv)
{
  // Declare functions for interpolation
  
  TFile *inthfiles[IMAX*RMAX*SMAX];
  TGraph *cfs00s[IMAX*RMAX*SMAX];

  TFile *inthssfiles[RMAX];
  TGraph *cfs00sss[RMAX];
  
  // Get baseline file
  
  // TFile *inbasefile = new TFile("shmout.recalc.kp.base.root");
  // TH1D *cbase = (TH1D *) inbasefile->Get("CfnReYlm00CYlmBase");

  Double_t val;

  // Load functions for interpolation
  // And divide by baseline
  
  for (int is=0; is<SMAX; is++) {
    for (int ir=0; ir<RMAX; ir++) {
      for (int ii=0; ii<IMAX; ii++) {
	
	inthfiles[is*RMAX*IMAX + ir*IMAX + ii] = new TFile (Form("kpcoulstrongwgtlhc.rad%.1f.re%.3f.im%.3f.root",svals[is], rvals[ir], ivals[ii]));
	//      inthfiles[is*RMAX*IMAX + ir*IMAX + ii] = new TFile (Form("shmout.recalc.kp.R%i.I%i.root",ir+1+(ir>1)*1, ii+1));
	
	cfs00s[is*RMAX*IMAX + ir*IMAX + ii] = (TGraph *) inthfiles[is*RMAX*IMAX + ir*IMAX + ii]->Get("grkaonproton");
	
	if (cfs00s[is*RMAX*IMAX + ir*IMAX + ii]) {
	  cfs00sbgs  [is*RMAX*IMAX + ir*IMAX + ii] = new TGraph(*cfs00s[is*RMAX*IMAX + ir*IMAX + ii]);
	  //	cfs00sbgs[is*RMAX*IMAX + ir*IMAX + ii]->Divide(cbase);
	  cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii] = smearmr(cfs00sbgs[is*RMAX*IMAX + ir*IMAX + ii]);
	}
	else {
	  cfs00sbgs[is*RMAX*IMAX + ir*IMAX + ii] = 0;
	  cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii] = 0;
	}
      }
    }
    inthssfiles[is] = new TFile (Form("kpcoulwgtlhc.rad%.1f.re0.001.im0.001.root",svals[is]));
    
    cfs00sss[is] = (TGraph *) inthssfiles[is]->Get("grkaonproton");

    
    if (cfs00sss[is]) {
      cfs00sssbgs  [is] = new TGraph(*cfs00sss[is]);
      cfs00sssbgsmr[is] = smearmr(cfs00sssbgs[is]);
    }
    else {
      cfs00sssbgs[is] = 0;
      cfs00sssbgsmr[is] = 0;
    }
    cout << "Got ss graphs " << cfs00sss[is] << " " << cfs00sssbgsmr[is] << endl;

  }

  // Get data for plus-plus combination
  
  // TFile *infilepp = new TFile("shmout.recalc.kp.pp.root");
  // TH1D *cfs00pp = (TH1D *) infilepp->Get("CfnReYlm00CYlm");
  // cfs00ppbg = new TH1D(*cfs00pp);
  // cfs00ppbg->Divide(cbase);

  // Smear by momentum resolution
  
  //  cfs00ppbgmr = smearmr(cfs00ppbg);
  // cfs00ppbgmr = new TH1D(*cfs00ppbg);
  
  // Get data for plus-minus combination
  
  // TFile *indatafile = new TFile("shmout.recalc.c1.mp.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KmPcent1NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.c1.pm.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KpApcent1NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.c2.mp.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KmPcent2NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.c2.pm.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00KpApcent2NonIdYlms");
  // TFile *indatafile = new TFile("shmout.recalc.mpos.m1.kmpp.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKmProtpcM1");
  // TFile *indatafile = new TFile("shmout.recalc.ff.cut1.KmPp.m0.root");
  // cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKmProtpcM0");
//  TFile *indatafile = new TFile("shmout.recalc.ff.cut1.KpPm.m0.root");
//  cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKpAProtpcM0");
   // TFile *indatafile = new TFile("shmout.recalc.mpos.m1.kppm.root");
  //cdata = (TH1D *) indatafile->Get("CfnReYlm00cylmKpAProtpcM1");

  // Recent ALICE data
  // TFile *indatafilerfmp = new TFile("shmout.recalc.mmin.m0.kmpp.5percentMFcut.root");
  // cdatarfmp = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  // TFile *indatafilerfpm = new TFile("shmout.recalc.mmin.m0.kppm.5percentMFcut.root");
  // cdatarfpm = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  // TFile *indatafileffmp = new TFile("shmout.recalc.mpos.m0.kmpp.5percentMFcut.root");
  // cdataffmp = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  // TFile *indatafileffpm = new TFile("shmout.recalc.mpos.m0.kppm.5percentMFcut.root");
  // cdataffpm = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // // Data from Wioletta

  // TFile *indatafilerfmp;
  // TFile *indatafilerfpm;
  // TFile *indatafileffmp;
  // TFile *indatafileffpm;
  
  // // Centrality 0
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c0.mp.root");
  // cdatarfmp[0] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c0.pm.root");
  // cdatarfpm[0] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c0.mp.root");
  // cdataffmp[0] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c0.pm.root");
  // cdataffpm[0] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // // Centrality 1
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c1.mp.root");
  // cdatarfmp[1] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM1");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c1.pm.root");
  // cdatarfpm[1] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM1");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c1.mp.root");
  // cdataffmp[1] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM1");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c1.pm.root");
  // cdataffpm[1] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM1");

  // // Centrality 2
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c2.mp.root");
  // cdatarfmp[2] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM2");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c2.pm.root");
  // cdatarfpm[2] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM2");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c2.mp.root");
  // cdataffmp[2] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM2");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c2.pm.root");
  // cdataffpm[2] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM2");

  // // Centrality 3
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c3.mp.root");
  // cdatarfmp[3] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM3");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c3.pm.root");
  // cdatarfpm[3] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM3");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c3.mp.root");
  // cdataffmp[3] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM3");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c3.pm.root");
  // cdataffpm[3] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM3");

  // // Centrality 4
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c4.mp.root");
  // cdatarfmp[4] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM4");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c4.pm.root");
  // cdatarfpm[4] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM4");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c4.mp.root");
  // cdataffmp[4] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM4");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c4.pm.root");
  // cdataffpm[4] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM4");

  // // Centrality 5
  // indatafilerfmp = new TFile("shmout.recalc.mmin.c5.mp.root");
  // cdatarfmp[5] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM5");
  // indatafilerfpm = new TFile("shmout.recalc.mmin.c5.pm.root");
  // cdatarfpm[5] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM5");
  // indatafileffmp = new TFile("shmout.recalc.mplu.c5.mp.root");
  // cdataffmp[5] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM5");
  // indatafileffpm = new TFile("shmout.recalc.mplu.c5.pm.root");
  // cdataffpm[5] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM5");

  // // Data on same-sign correlations
  
  // TFile *indatafilerfpp;
  // TFile *indatafilerfmm;
  // TFile *indatafileffpp;
  // TFile *indatafileffmm;
  
  // // Centrality 0
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c0.pp.root");
  // cdatarfpp[0] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c0.mm.root");
  // cdatarfmm[0] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c0.pp.root");
  // cdataffpp[0] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c0.mm.root");
  // cdataffmm[0] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // // Centrality 1
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c1.pp.root");
  // cdatarfpp[1] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM1");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c1.mm.root");
  // cdatarfmm[1] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM1");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c1.pp.root");
  // cdataffpp[1] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM1");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c1.mm.root");
  // cdataffmm[1] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM1");

  // // Centrality 2
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c2.pp.root");
  // cdatarfpp[2] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM2");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c2.mm.root");
  // cdatarfmm[2] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM2");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c2.pp.root");
  // cdataffpp[2] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM2");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c2.mm.root");
  // cdataffmm[2] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM2");

  // // Centrality 3
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c3.pp.root");
  // cdatarfpp[3] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM3");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c3.mm.root");
  // cdatarfmm[3] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM3");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c3.pp.root");
  // cdataffpp[3] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM3");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c3.mm.root");
  // cdataffmm[3] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM3");

  // // Centrality 4
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c4.pp.root");
  // cdatarfpp[4] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM4");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c4.mm.root");
  // cdatarfmm[4] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM4");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c4.pp.root");
  // cdataffpp[4] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM4");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c4.mm.root");
  // cdataffmm[4] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM4");

  // // Centrality 5
  // indatafilerfpp = new TFile("shmout.recalc.mmin.c5.pp.root");
  // cdatarfpp[5] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM5");
  // indatafilerfmm = new TFile("shmout.recalc.mmin.c5.mm.root");
  // cdatarfmm[5] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM5");
  // indatafileffpp = new TFile("shmout.recalc.mplu.c5.pp.root");
  // cdataffpp[5] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM5");
  // indatafileffmm = new TFile("shmout.recalc.mplu.c5.mm.root");
  // cdataffmm[5] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM5");

  // New Data from Wioletta from 10.04.2019

  // TFile *indatafilerfmp = new TFile("corrs.alice.5TeV.1904.minus.root");
  // TFile *indatafilerfpm = new TFile("corrs.alice.5TeV.1904.minus.root");
  // TFile *indatafileffmp = new TFile("corrs.alice.5TeV.1904.plus.root");
  // TFile *indatafileffpm = new TFile("corrs.alice.5TeV.1904.plus.root");

  // New Data from Wioletta from 24.04.2019 TPC-only tracks
  
  TFile *indatafilerfmp = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilerfpm = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafileffmp = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  TFile *indatafileffpm = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  
  // Centrality 0
  cdatarfmp[0] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[0] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[0] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[0] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 1
  cdatarfmp[1] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[1] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[1] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[1] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 2
  cdatarfmp[2] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[2] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[2] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[2] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 3
  cdatarfmp[3] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[3] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[3] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[3] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 4
  cdatarfmp[4] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[4] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[4] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[4] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 5
  cdatarfmp[5] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdatarfpm[5] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  cdataffmp[5] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  cdataffpm[5] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Data on same-sign correlations
  
  TFile *indatafilerfpp = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafilerfmm = new TFile("corrs.alice.5TeV.1904.tpconly.minus.root");
  TFile *indatafileffpp = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  TFile *indatafileffmm = new TFile("corrs.alice.5TeV.1904.tpconly.plus.root");
  
  // Centrality 0
  cdatarfpp[0] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[0] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[0] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[0] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 1
  cdatarfpp[1] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[1] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[1] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[1] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 2
  cdatarfpp[2] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[2] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[2] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[2] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 3
  cdatarfpp[3] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[3] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[3] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[3] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 4
  cdatarfpp[4] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[4] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[4] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[4] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  // Centrality 5
  cdatarfpp[5] = (TH1D *) indatafilerfpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdatarfmm[5] = (TH1D *) indatafilerfmm->Get("CfnReYlm00cylmKmAProtpcM0");
  cdataffpp[5] = (TH1D *) indatafileffpp->Get("CfnReYlm00cylmKpProtpcM0");
  cdataffmm[5] = (TH1D *) indatafileffmm->Get("CfnReYlm00cylmKmAProtpcM0");

  //  cdata = cdatarfmp[0];
    
  // shmout.recalc.mmin.m0.kmpm.5percentMFcut.root
  //   shmout.recalc.mmin.m0.kppp.5percentMFcut.root
  //   shmout.recalc.mpos.m0.kmpm.5percentMFcut.root
  //   shmout.recalc.mpos.m0.kppp.5percentMFcut.root
  
  
  Double_t xmax = cdatarfmp[0]->GetXaxis()->GetXmax();

  double purval = 0.85;

  double nrmvals[6];
  double lbdvals[6];
  double wdtvals[6];
  
  TF1 *funn = new TF1("funn","[0]");
  TF1 *fung = new TF1("fung","[0] + [1]*exp(-x*x*[2]*[2])");

  double normvalsfix[24];
  
  for (int icent=0; icent<6; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //    if (docent[icent]) {

      cout << "Fitting ul back " << icent << " " << ipair << endl;
      
      funn->SetParameter(0,1.000);
      cdatarfmp[icent]->Fit(funn, "", "", 0.85*xmax, xmax);
      
      fung->SetParameters(funn->GetParameter(0),0.002,2.0);
      fung->FixParameter(0,funn->GetParameter(0));
      //  fung->FixParameter(2,2.0);
      cdatarfmp[icent]->Fit(fung, "","",0.25, xmax);

      normvalsfix[icent*4+ipair] = fung->GetParameter(0);
      
      nrmvals[icent] = fung->GetParameter(0);
      lbdvals[icent] = fung->GetParameter(1);
      wdtvals[icent] = fung->GetParameter(2);
      //    }
    }
  }

  // for (int icent=0; icent<6; icent++) {
  //   cdatarfmp[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsfix[icent*4+0] = funn->GetParameter(0);
  //   cdatarfpm[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsfix[icent*4+1] = funn->GetParameter(0);
  //   cdataffmp[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsfix[icent*4+2] = funn->GetParameter(0);
  //   cdataffpm[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsfix[icent*4+3] = funn->GetParameter(0);
  // }
  
  double nrmvalsss[6];
  double lbdvalsss[6];
  double wdtvalsss[6];
  
  double normvalsssfix[24];
  
  for (int icent=0; icent<6; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //    if (docent[icent]) {

      //      cout << "Fitting ss " << icent << endl;
      cout << "Fitting ss back " << icent << " " << ipair << endl;
    
      funn->SetParameter(0,1.000);
      cdatarfpp[icent]->Fit(funn, "", "", 0.85*xmax, xmax);
      
      fung->SetParameters(funn->GetParameter(0),0.002,2.0);
      fung->FixParameter(0,funn->GetParameter(0));
      //  fung->FixParameter(2,2.0);
      cdatarfpp[icent]->Fit(fung, "","",0.3, xmax);
      
      normvalsssfix[icent*4+ipair] = funn->GetParameter(0);
      nrmvalsss[icent] = fung->GetParameter(0);
      lbdvalsss[icent] = fung->GetParameter(1);
      wdtvalsss[icent] = fung->GetParameter(2);
      //    }
    }
  }
  
  // for (int icent=0; icent<6; icent++) {
  //   cdatarfpp[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsssfix[icent*4+0] = funn->GetParameter(0);
  //   cdatarfmm[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsssfix[icent*4+1] = funn->GetParameter(0);
  //   cdataffpp[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsssfix[icent*4+2] = funn->GetParameter(0);
  //   cdataffmm[icent]->Fit(funn,"","",0.85*xmax,xmax);
  //   normvalsssfix[icent*4+3] = funn->GetParameter(0);
  // }
  

  
  // funn->SetParameter(0,1.000);
  // cdatarfmp[3]->Fit(funn, "", "", 0.75*xmax, xmax);
  
  // fung->SetParameters(funn->GetParameter(0),0.002,2.0);
  // fung->FixParameter(0,funn->GetParameter(0));
  // //  fung->FixParameter(2,2.0);
  // cdatarfmp[3]->Fit(fung, "","",0.25, xmax);
  
  // nrmvals[3] = fung->GetParameter(0);
  // lbdvals[3] = fung->GetParameter(1);
  // wdtvals[3] = fung->GetParameter(2);

  // Get data for plus-plus combination
  
  // TFile *indatppfile = new TFile("shmout.recalc.mpos.m1.kppp.root");
  // cdatapp = (TH1D *) indatppfile->Get("CfnReYlm00cylmKpProtpcM1");
  // TFile *indatppfile = new TFile("shmout.recalc.c1.pp.root");
  // cdatapp = (TH1D *) indatppfile->Get("CfnReYlm00KpPcent1NonIdYlms");

  // Create Chi2 map
  
  TGraph2D *grchi = new TGraph2D();
  grchi->SetName("grchi");
  grchi->SetTitle("#chi^{2} Map;#Rgothicf_{0} (fm);#Jgothicf_{0} (fm);#chi^{2}");

  double chival;
  int pcount = 0;

  double minchi = 1e9;

  int minir, minii;
  int minirsum=0, miniisum=0, centcount=0;
  int minirs[6], miniis[6], minszs[6];

  for (int icent=0; icent<6; icent++) {
    minirs[icent] = 1;
    miniis[icent] = 1;
    minszs[icent] = 1;
  }
  
  // Create a fit chi2 map

  for (int icent=0; icent<6; icent++) {
    if (docent[icent]) {
      minchi = 1e9;
      for (int is=0; is<SMAX; is++) {
	for (int ir=0; ir<RMAX; ir++) {
	  for (int ii=0; ii<IMAX; ii++) {
	    if (cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii]) {
	      chival  = getchi2(cdatarfmp[icent], cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii]->GetY(), nrmvals[icent], purval, lbdvals[icent], wdtvals[icent]);
	  
	      //	  grchi->SetPoint(pcount++, rvals[ir], ivals[ii], chival);
	      
	      if (chival < minchi) 
		{ minchi = chival;
		  minirs[icent] = ir;
		  miniis[icent] = ii;
		  minszs[icent] = is;
		  cout << "Now minimum is " << minchi << "    s r i " << is << " "<< ir << " "<< ii << endl;
		}
	    }
	  }
	}
	minirsum += minirs[icent];
	miniisum += miniis[icent];
	centcount++;
      }
    }
  }

  minir = minirsum/centcount;
  minii = miniisum/centcount;
  

  TFile *outfile = new TFile("outchicentallwgt.root","RECREATE");
  outfile->cd();  
  
  cout << "Minimum chi " << minchi << endl;

  cout << "Re Im Size " << rvals[minir] << " " << ivals[minii] << " ";

  for (int icent=0; icent<6; icent++)
    if (docent[icent])
      cout << "Centrality " << icent << "  " << svals[minszs[icent]] << endl;

  if (minir == RMAX-1) minir--;
  if (minii == IMAX-1) minii--;
  
  for (int icent=0; icent<6; icent++)
    if (minszs[icent] == SMAX-1) minszs[icent]--;

  double vals[500];
  //  getlinext(-0.32, 0.55, 3.5, vals);

  // getlinext(-0.7, 0.3, 3.5, vals);
  // TGraph *grr1 = new TGraph(500, cfs00s[0]->GetX(), vals);
  // grr1->SetName("grr1");
  
  // getlinext(-0.5, 0.3, 4.0, vals);
  // TGraph *grr2 = new TGraph(500, cfs00s[0]->GetX(), vals);
  // grr2->SetName("grr2");
  
  // getlinext(-0.6, 0.3, 3.70, vals);
  // TGraph *grr3 = new TGraph(500, cfs00s[0]->GetX(), vals);
  // grr3->SetName("grr3");

  //  exit(0);

  cout << "Now fitting " << endl;
  
  // Values from Therminator2 fits
  // lbdval = 0.00555;
  // wdtval = 1.92;
  
  cout << "Setting values " << nrmvals[0] << " " << lbdvals[0] << " " << wdtvals[0] << " " << purval << endl;
  
  TVirtualFitter *fitter=TVirtualFitter::Fitter(0, 75);
  fitter->SetFCN(myfuncf);
  fitter->SetParameter(0,  "Purity"  ,purval        ,0.0001,0.2,   0.85);
  fitter->SetParameter(1,  "Ref0"    ,rvals[minir]  ,0.0001,-5.0, 5.0);
  fitter->SetParameter(2,  "Imf0"    ,ivals[minii]  ,0.0001,0.00001, 5.0);
  // fitter->SetParameter(2,  "Imf0"    ,0.2 ,0.0001,0.00001, 5.0);
  // fitter->FixParameter(2);
  
  
  for (int icent=0; icent<6; icent++) {
    // fitter->SetParameter(3+icent*3+0,  Form("Lambda%i", icent)  ,(TMath::Abs(lbdvals[icent])+TMath::Abs(lbdvalsss[icent]))/2.0    ,0.001, 0.0,   0.1);
    // fitter->SetParameter(3+icent*3+1,  Form("Width%i", icent)   ,(TMath::Abs(wdtvals[icent])+TMath::Abs(wdtvalsss[icent]))/2.0    ,0.001, 1.0,   10.0);
    fitter->SetParameter(3+icent*3+0,  Form("Lambda%i", icent)  ,TMath::Abs(lbdvalsss[icent])   ,0.001, 0.0,   0.1);
    fitter->SetParameter(3+icent*3+1,  Form("Width%i", icent)   ,TMath::Abs(wdtvalsss[icent])   ,0.001, 1.0,   10.0);
    fitter->SetParameter(3+icent*3+2,  Form("Radius%i", icent)  ,svals[minszs[icent]]  ,0.01  ,0.9, 8.0);

    if (!docent[icent]) {
      fitter->FixParameter(3+icent*3+0);
      fitter->FixParameter(3+icent*3+1);
      fitter->FixParameter(3+icent*3+2);
    }

    fitter->FixParameter(3+icent*3+0);
    fitter->FixParameter(3+icent*3+1);
  }

  for (int icent=0; icent<6; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //      fitter->SetParameter(21+icent*4+ipair,  Form("Norm%i%i", icent, ipair)    ,nrmvals[icent]    ,0.001, 0.9,   1.1);
      fitter->SetParameter(21+icent*4+ipair,  Form("Norm%i%i", icent, ipair)    ,normvalsfix[icent*4+ipair]    ,0.001, 0.9,   1.1);      if (!docent[icent])
	fitter->FixParameter(21+icent*4+ipair);

      fitter->FixParameter(21+icent*4+ipair);
    }
  }
  
  for (int icent=0; icent<6; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      //      fitter->SetParameter(45+icent*4+ipair,  Form("Normss%i%i", icent, ipair)    ,nrmvalsss[icent]    ,0.001, 0.9,   1.1);
      fitter->SetParameter(45+icent*4+ipair,  Form("Normss%i%i", icent, ipair)    ,normvalsssfix[icent*4+ipair]    ,0.001, 0.9,   1.1);      
      if (!docent[icent])
	fitter->FixParameter(45+icent*4+ipair);
      fitter->FixParameter(45+icent*4+ipair);
    }
  }

  //  fitter->SetParameter(69, "BGsstoul", lbdvalsss[0]/lbdvals[0], 0.001, 0.1, 2.0);
  for (int icent=0; icent<6; icent++) {
    fitter->SetParameter(69+icent, Form("BGsstoul%i", icent), 0.8, 0.001, 0.1, 5.0);

    //    fitter->FixParameter(69+icent);
  }
  
  
  // fitter->SetParameter(0,  "Norm"    ,nrmvals[0]    ,0.001, 0.9,   1.1);
  // fitter->SetParameter(7,  "nrmsc1"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(8,  "nrmsc2"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(9,  "nrmsc3"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(10, "Norm3"   ,nrmvals[3]    ,0.001, 0.9,   1.1);
  // fitter->SetParameter(11, "Lambda3" ,lbdvals[3]    ,0.001, 0.0,   0.1);
  // fitter->SetParameter(12, "Width3"  ,wdtvals[3]    ,0.001, 1.0,   10.0);
  // fitter->SetParameter(13, "Radius3" ,svals[minsz3] ,0.01  ,0.9, 8.0);
  // fitter->SetParameter(14, "nrmsc4"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(15, "nrmsc5"  ,1.0           ,0.01  ,0.9, 1.1);
  // fitter->SetParameter(16, "nrmsc6"  ,1.0           ,0.01  ,0.9, 1.1);
  //  fitter->SetParameter(5, "Imf0"    ,0.25  ,0.001,0.0, 5.0);
  //  fitter->SetParameter(6, "norm"     ,1.0       ,0.0001,0.0,   2.0);


  
  Double_t arglist[100];
  arglist[0] = 1;
  fitter->ExecuteCommand("CALL FCN", arglist, 1);
  fitter->FixParameter(0);
  //  fitter->FixParameter(1);
  //fitter->FixParameter(2);
  //fitter->FixParameter(3);
  //  fitter->FixParameter(4);
  //  fitter->FixParameter(5);
  //  fitter->FixParameter(6);
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

  cout << "Fitting done " << " " << fitter->GetParameter(4) << " " << fitter->GetParameter(5) << "   " << fitter->GetParameter(6) << endl;

  // double fitnorm = fitter->GetParameter(0);
  // double fitlabd = fitter->GetParameter(1);
  // double fitwdth = fitter->GetParameter(2);
  double fitpuri = fitter->GetParameter(0); 

  double fitref0 = fitter->GetParameter(1);
  double fitimf0 = fitter->GetParameter(2);
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
  double fitlabds[6];
  double fitwdths[6];
  double fitradis[6];

  double fitnrmscales[4*6];
  double fitnrmscalesss[4*6];

  for (int icent=0; icent<6; icent++) {
    fitlabds[icent] = fitter->GetParameter(3+icent*3+0);
    fitwdths[icent] = fitter->GetParameter(3+icent*3+1);
    fitradis[icent] = fitter->GetParameter(3+icent*3+2);

    for (int ipair=0; ipair<4; ipair++) {
      fitnrmscales[icent*4+ipair] = fitter->GetParameter(21+icent*4+ipair);
      fitnrmscalesss[icent*4+ipair] = fitter->GetParameter(45+icent*4+ipair);
    }
  }

  double fitbgscales[6];

  for (int icent=0; icent<6; icent++) {
    fitbgscales[icent] = fitter->GetParameter(69+icent);
  }
  
  double fitepuri = fitter->GetParError(0); 

  double fiteref0 = fitter->GetParError(1);
  double fiteimf0 = fitter->GetParError(2);
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

  //  double fitnorms[6];
  double fitelabds[6];
  double fitewdths[6];
  double fiteradis[6];

  double fitenrmscales[4*6];
  double fitenrmscalesss[4*6];

  for (int icent=0; icent<6; icent++) {
    fitelabds[icent] = fitter->GetParError(3+icent*3+0);
    fitewdths[icent] = fitter->GetParError(3+icent*3+1);
    fiteradis[icent] = fitter->GetParError(3+icent*3+2);

    for (int ipair=0; ipair<4; ipair++) {
      fitenrmscales[icent*4+ipair] = fitter->GetParError(21+icent*4+ipair);
      fitenrmscalesss[icent*4+ipair] = fitter->GetParError(45+icent*4+ipair);
    }
  }

  double fitebgscales[6];

  for (int icent=0; icent<6; icent++) {
    fitebgscales[icent] = fitter->GetParError(69+icent);
  }
    
  cout << "Final paramerers:" << endl;
  cout << "Purity          : " << fitpuri << " +/- " << fitepuri <<  endl;
  cout << "Re f0           : " << fitref0 << " +/- " << fiteref0 <<  endl;
  cout << "Im f0           : " << fitimf0 << " +/- " << fiteimf0 <<  endl;

  for (int icent=0; icent<6; icent++) {
    cout << "Radius C" << icent << "     :" << fitradis[icent]  << " +/- " << fiteradis[icent] << endl;
  }

  for (int icent=0; icent<6; icent++) {
    cout << "Backgroung parameters centrality " << icent << endl;
    cout << "alfa width " << fitlabds[icent] << " " << fitwdths[icent] << endl;
    cout << "Normalizations :"
	 << fitnrmscales[icent*4+0] << " "
	 << fitnrmscales[icent*4+1] << " "
	 << fitnrmscales[icent*4+2] << " "
	 << fitnrmscales[icent*4+3] << " "
	 << endl;
    cout << "Normalizations same-sign:"
	 << fitnrmscalesss[icent*4+0] << " "
	 << fitnrmscalesss[icent*4+1] << " "
	 << fitnrmscalesss[icent*4+2] << " "
	 << fitnrmscalesss[icent*4+3] << " "
	 << endl;

  }

  for (int icent=0; icent<6; icent++) {
    cout << "Background ss to ul scale " << icent << " : " << fitbgscales[icent] << " +/- " << fitebgscales[icent] << endl;
  }
  
  outfile->cd();

  TGraph *grcfits[6];
  TGraph *grcfitsss[6];
  
  for (int icent=0; icent<6; icent++) {
    if (docent[icent]) {
      
      getlinext(fitref0, fitimf0, fitradis[icent], vals);
      
      cout << "First bins " << vals[0] << " " << vals[1] << " " << vals[2] << endl;

      grcfits[icent] = new TGraph(200, cfs00sbgsmr[0]->GetX(), vals);
      cout << "Have graph " << grcfits[icent] << " with " << grcfits[icent]->GetN() << " bins " << endl;
      grcfits[icent]->SetTitle(";k^ (GeV/#it{c});C(k*)");
      grcfits[icent]->SetName(Form("graulfitC%i", icent));
      
      grcfits[icent]->Write();

      cout << "Doing linear extrapolation " << endl;
      
      getlinextss(fitradis[icent], vals);

      cout << "First bins " << vals[0] << " " << vals[1] << " " << vals[2] << endl;
      
      grcfitsss[icent] = new TGraph(200, cfs00sssbgsmr[0]->GetX(), vals);
      cout << "Have graph " << grcfitsss[icent] << " with " << grcfitsss[icent]->GetN() << " bins " << endl;

      grcfitsss[icent]->SetTitle(";k^ (GeV/#it{c});C(k*)");
      grcfitsss[icent]->SetName(Form("grassfitC%i", icent));
      
      grcfitsss[icent]->Write();

    }
  }
  
  // getlinext(fitref0, fitimf0, fitrad3, vals);
  
  // TGraph *grfitc3 = new TGraph(500, cfs00s[0]->GetX(), vals);
  // grfitc3->SetTitle(";k^ (GeV/#it{c});C(k*)");
  // grfitc3->SetName("grappfitC3");

  // grfitc3->Write();


  // cfs00ppbg->SetName("cfs00ppbg");
  // cfs00ppbg->Write();
  // cfs00ppbgmr->SetName(Form("%s_mr",cfs00ppbg->GetName()));
  // cfs00ppbgmr->Write();

  TH1D *cfitrfpmlin[6];
  TH1D *cfitrfmplin[6];
  TH1D *cfitffpmlin[6];
  TH1D *cfitffmplin[6];

  TH1D *cfitrfpplin[6];
  TH1D *cfitrfmmlin[6];
  TH1D *cfitffpplin[6];
  TH1D *cfitffmmlin[6];

  double ksv, bgv;
  
  for (int icent=0; icent<6; icent++) {
    if (docent[icent]) {
    
      cfitrfpmlin[icent] = new TH1D(*cdatarfpm[icent]);
      cfitrfpmlin[icent]->SetName(Form("datarfpmbgC%i",icent));
      cfitrfmplin[icent] = new TH1D(*cdatarfmp[icent]);
      cfitrfmplin[icent]->SetName(Form("datarfmpbgC%i",icent));
      cfitffpmlin[icent] = new TH1D(*cdataffpm[icent]);
      cfitffpmlin[icent]->SetName(Form("dataffpmbgC%i",icent));
      cfitffmplin[icent] = new TH1D(*cdataffmp[icent]);
      cfitffmplin[icent]->SetName(Form("dataffmpbgC%i",icent));
      
      for (int ib=1; ib<=cdatarfmp[0]->GetNbinsX(); ib++) {
	ksv = cdatarfmp[icent]->GetXaxis()->GetBinCenter(ib);
	bgv = fitnrmscales[icent*4+0] + fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
      
	val = cdatarfmp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitrfmplin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscales[icent*4+1] + fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdatarfpm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitrfpmlin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscales[icent*4+2] + fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdataffmp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitffmplin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscales[icent*4+3] + fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdataffpm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitffpmlin[icent]->SetBinContent(ib, val);
      }
    
      cfitrfmplin[icent]->Write();
      cfitrfpmlin[icent]->Write();
      cfitffmplin[icent]->Write();
      cfitffpmlin[icent]->Write();

      cfitrfpplin[icent] = new TH1D(*cdatarfpp[icent]);
      cfitrfpplin[icent]->SetName(Form("datarfppbgC%i",icent));
      cfitrfmmlin[icent] = new TH1D(*cdatarfmm[icent]);
      cfitrfmmlin[icent]->SetName(Form("datarfmmbgC%i",icent));
      cfitffpplin[icent] = new TH1D(*cdataffpp[icent]);
      cfitffpplin[icent]->SetName(Form("dataffppbgC%i",icent));
      cfitffmmlin[icent] = new TH1D(*cdataffmm[icent]);
      cfitffmmlin[icent]->SetName(Form("dataffmmbgC%i",icent));
      
      for (int ib=1; ib<=cdatarfmm[0]->GetNbinsX(); ib++) {
	ksv = cdatarfpp[icent]->GetXaxis()->GetBinCenter(ib);
	bgv = fitnrmscalesss[icent*4+0] + fitbgscales[icent]*fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
      
	val = cdatarfpp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitrfpplin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscalesss[icent*4+1] + fitbgscales[icent]*fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdatarfmm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitrfmmlin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscalesss[icent*4+2] + fitbgscales[icent]*fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdataffpp[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitffpplin[icent]->SetBinContent(ib, val);
	
	bgv = fitnrmscalesss[icent*4+3] + fitbgscales[icent]*fitlabds[icent]*TMath::Exp(-ksv*ksv*fitwdths[icent]*fitwdths[icent]);
	
	val = cdataffmm[icent]->GetBinContent(ib);
	val += -bgv + 1.0;
	val = (val-1.0)/fitpuri + 1.0;
	
	cfitffmmlin[icent]->SetBinContent(ib, val);
      }
    
      cfitrfmmlin[icent]->Write();
      cfitrfpplin[icent]->Write();
      cfitffmmlin[icent]->Write();
      cfitffpplin[icent]->Write();

    }
  }

  Double_t fvx[75];
  Double_t fvy[75];
  Double_t fve[75];

  for (int ipar=0; ipar<76; ipar++) {
    fvx[ipar] = ipar;
    fvy[ipar] = fitter->GetParameter(ipar);
    fve[ipar] = fitter->GetParError(ipar);
  }

  TGraphErrors *fitresult = new TGraphErrors(75, fvx, fvy, NULL, fve);
  fitresult->SetName("fitresult");
  fitresult->Write();
  
  // TH1D *cfitpplin = new TH1D(*cdatapp);
  // cfitpplin->SetName("datappbg");

  // for (int ib=1; ib<=cdatapp->GetNbinsX(); ib++)
  //   {
  //     val = cdatapp->GetBinContent(ib);
  //     ksv = cdatapp->GetXaxis()->GetBinCenter(ib);
  //     bgv = nrmvals[0] + lbdvals[0]*TMath::Exp(-ksv*ksv*wdtvals[0]*wdtvals[0]);
      
  //     val += -bgv + 1.0;
  //     val = (val-1.0)/purval + 1.0;

  //     cfitpplin->SetBinContent(ib, val);
  //   }

  // cfitpplin->Write();

  // grr1->Write();
  // grr2->Write();
  // grr3->Write();
  
}
 
