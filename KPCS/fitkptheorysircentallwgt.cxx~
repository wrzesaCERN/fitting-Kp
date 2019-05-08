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
#define SMAX 8

using namespace std;

Double_t rvals[RMAX] = { -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.001, 0.1, 0.2 };
//Double_t rvals[RMAX] = { -1.0, -0.7, -0.5, -0.4, -0.3, -0.25, 0.001, 0.5 };
Double_t ivals[IMAX] = { 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
Double_t svals[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0 };
//Double_t svals[SMAX] = { 2.5, 3.0, 3.5, 4.0, 4.5 };

Int_t docent[6] = { 1, 1, 1, 1, 1, 1 };

TGraph *cfs00sbgs[IMAX*RMAX*SMAX];
TGraph *cfs00sbgsmr[IMAX*RMAX*SMAX];
TH1D *cfs00ppbg;

TH1D *cfs00ppbgmr;

//TH1D *cdata;

TH1D *cdataffpm[6];
TH1D *cdataffmp[6];
TH1D *cdatarfpm[6];
TH1D *cdatarfmp[6];

TH1D *cdatapp;

TGraph *smearmr(TGraph *cfppdv)
{
  double width = 0.006;
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

//  double getchi2(TH1D *data, TH1D *model, double norm, double purity, double lambda, double radbg)
// {
//   double sumchi = 0;
//   double val, err, chi, thv, bgv, ksv;

//   int maxbin = 65;
//   int minbin = 3;

//   for (int ib=minbin; ib<=maxbin; ib++)
//     {
//       if (ib == 47) continue;
//       if (ib == 48) continue;
//       if (ib == 49) continue;
      
//       val = data->GetBinContent(ib);
//       ksv = data->GetXaxis()->GetBinCenter(ib);
//       bgv = norm + lambda*TMath::Exp(-ksv*ksv*radbg*radbg);
      
//       val += -bgv + 1.0;
//       val = (val-1.0)/purity + 1.0;

//       thv = model->GetBinContent(ib);
//       err = data->GetBinError(ib);

//       chi = (val-thv)*(val-thv)/(err*err);
      
//       //      cout << "val thv chi " << val << "   " << thv << "   " << chi << endl;

//       sumchi += chi;
//     }

//   cout << "full chi2 " << sumchi << "     ----       " << norm << " "<< lambda << " " << radbg << " " << purity << " " << endl;

//   return sumchi;
// }

double getchi2(TH1D *data, double *model, double norm, double purity, double lambda, double radbg)
{
  double sumchi = 0;
  double val, err, chi, thv, bgv, ksv;

  int maxbin = 85;
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
  Double_t Pv  = par[0];
  Double_t Rev = par[1];
  Double_t Imv = par[2];

  Double_t Lvs[6];
  Double_t Rvs[6];
  Double_t Svs[6];

  Double_t Nvs[24];

  for (int iter=0; iter<6; iter++) {
    Lvs[iter] = par[3+iter*3+0];
    Rvs[iter] = par[3+iter*3+1];
    Svs[iter] = par[3+iter*3+2];
  }

  for (int iter=0; iter<24; iter++) {
    Nvs[iter] = par[21+iter];
  }
  
  // Double_t Nv  = par[0];
  // Double_t Lv  = par[1];
  // Double_t Rv  = par[2];
  
  // Double_t Szv = par[6]; 

  // Double_t Ns1 = par[7];
  // Double_t Ns2 = par[8];
  // Double_t Ns3 = par[9];

  // Double_t Nv3 = par[10];
  // Double_t Lv3 = par[11];
  // Double_t Rv3 = par[12];

  // Double_t Sz3 = par[13];

  // Double_t Ns4 = par[14];
  // Double_t Ns5 = par[15];
  // Double_t Ns6 = par[16];
  
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

  // Get baseline file
  
  TFile *inbasefile = new TFile("shmout.recalc.kp.base.root");
  TH1D *cbase = (TH1D *) inbasefile->Get("CfnReYlm00CYlmBase");

  Double_t val;

  // Load functions for interpolation
  // And divide by baseline
  
  for (int is=0; is<SMAX; is++) {
    for (int ir=0; ir<RMAX; ir++) {
      for (int ii=0; ii<IMAX; ii++) {
	
	inthfiles[is*RMAX*IMAX + ir*IMAX + ii] = new TFile (Form("kpcoulstrong.rad%.1f.re%.3f.im%.3f.root",svals[is], rvals[ir], ivals[ii]));
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

  // Data from Wioletta

  TFile *indatafilerfmp;
  TFile *indatafilerfpm;
  TFile *indatafileffmp;
  TFile *indatafileffpm;
  
  // Centrality 0
  indatafilerfmp = new TFile("shmout.recalc.mmin.c0.mp.root");
  cdatarfmp[0] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM0");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c0.pm.root");
  cdatarfpm[0] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM0");
  indatafileffmp = new TFile("shmout.recalc.mplu.c0.mp.root");
  cdataffmp[0] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM0");
  indatafileffpm = new TFile("shmout.recalc.mplu.c0.pm.root");
  cdataffpm[0] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM0");

  // Centrality 1
  indatafilerfmp = new TFile("shmout.recalc.mmin.c1.mp.root");
  cdatarfmp[1] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM1");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c1.pm.root");
  cdatarfpm[1] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM1");
  indatafileffmp = new TFile("shmout.recalc.mplu.c1.mp.root");
  cdataffmp[1] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM1");
  indatafileffpm = new TFile("shmout.recalc.mplu.c1.pm.root");
  cdataffpm[1] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM1");

  // Centrality 2
  indatafilerfmp = new TFile("shmout.recalc.mmin.c2.mp.root");
  cdatarfmp[2] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM2");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c2.pm.root");
  cdatarfpm[2] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM2");
  indatafileffmp = new TFile("shmout.recalc.mplu.c2.mp.root");
  cdataffmp[2] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM2");
  indatafileffpm = new TFile("shmout.recalc.mplu.c2.pm.root");
  cdataffpm[2] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM2");

  // Centrality 3
  indatafilerfmp = new TFile("shmout.recalc.mmin.c3.mp.root");
  cdatarfmp[3] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM3");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c3.pm.root");
  cdatarfpm[3] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM3");
  indatafileffmp = new TFile("shmout.recalc.mplu.c3.mp.root");
  cdataffmp[3] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM3");
  indatafileffpm = new TFile("shmout.recalc.mplu.c3.pm.root");
  cdataffpm[3] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM3");

  // Centrality 4
  indatafilerfmp = new TFile("shmout.recalc.mmin.c4.mp.root");
  cdatarfmp[4] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM4");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c4.pm.root");
  cdatarfpm[4] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM4");
  indatafileffmp = new TFile("shmout.recalc.mplu.c4.mp.root");
  cdataffmp[4] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM4");
  indatafileffpm = new TFile("shmout.recalc.mplu.c4.pm.root");
  cdataffpm[4] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM4");

  // Centrality 5
  indatafilerfmp = new TFile("shmout.recalc.mmin.c5.mp.root");
  cdatarfmp[5] = (TH1D *) indatafilerfmp->Get("CfnReYlm00cylmKmProtpcM5");
  indatafilerfpm = new TFile("shmout.recalc.mmin.c5.pm.root");
  cdatarfpm[5] = (TH1D *) indatafilerfpm->Get("CfnReYlm00cylmKpAProtpcM5");
  indatafileffmp = new TFile("shmout.recalc.mplu.c5.mp.root");
  cdataffmp[5] = (TH1D *) indatafileffmp->Get("CfnReYlm00cylmKmProtpcM5");
  indatafileffpm = new TFile("shmout.recalc.mplu.c5.pm.root");
  cdataffpm[5] = (TH1D *) indatafileffpm->Get("CfnReYlm00cylmKpAProtpcM5");

  
  
  //  cdata = cdatarfmp[0];
    
  // shmout.recalc.mmin.m0.kmpm.5percentMFcut.root
  //   shmout.recalc.mmin.m0.kppp.5percentMFcut.root
  //   shmout.recalc.mpos.m0.kmpm.5percentMFcut.root
  //   shmout.recalc.mpos.m0.kppp.5percentMFcut.root
  
  
  Double_t xmax = cdatarfmp[0]->GetXaxis()->GetXmax();

  double purval = 0.75;

  double nrmvals[6];
  double lbdvals[6];
  double wdtvals[6];
  
  TF1 *funn = new TF1("funn","[0]");
  TF1 *fung = new TF1("fung","[0] + [1]*exp(-x*x*[2]*[2])");

  for (int icent=0; icent<6; icent++) {
    //    if (docent[icent]) {
      
      funn->SetParameter(0,1.000);
      cdatarfmp[icent]->Fit(funn, "", "", 0.75*xmax, xmax);
      
      fung->SetParameters(funn->GetParameter(0),0.002,2.0);
      fung->FixParameter(0,funn->GetParameter(0));
      //  fung->FixParameter(2,2.0);
      cdatarfmp[icent]->Fit(fung, "","",0.25, xmax);
      
      nrmvals[icent] = fung->GetParameter(0);
      lbdvals[icent] = fung->GetParameter(1);
      wdtvals[icent] = fung->GetParameter(2);
      //    }
  }
  
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
  
  // int minir3, minii3, minsz3;

  // for (int is=0; is<SMAX; is++) {
  //   for (int ir=0; ir<RMAX; ir++) {
  //     for (int ii=0; ii<IMAX; ii++) {
  // 	if (cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii]) {
  // 	  //	chival = getchi2(cdata, cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii], nrmval, purval, lbdval, wdtval);
  // 	  //	  chival  = getchi2(cdatarfmp[0], cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii]->GetY(), nrmvals[0], purval, lbdvals[0], wdtvals[0]);
  // 	  chival = getchi2(cdatarfmp[3], cfs00sbgsmr[is*RMAX*IMAX + ir*IMAX + ii]->GetY(), nrmvals[3], purval, lbdvals[3], wdtvals[3]);
  // 	  grchi->SetPoint(pcount++, rvals[ir], ivals[ii], chival);
	  
  // 	  if (chival < minchi) 
  // 	    { minchi = chival;
  // 	      minir3 = ir;
  // 	      minii3 = ii;
  // 	      minsz3 = is;
  // 	      cout << "Now minimum is " << minchi << "    s r i " << is << " "<< ir << " "<< ii << endl;
  // 	    }
  // 	}
  //     }
  //   }
  // }
  
  // Create the best-fit graph
  
  // TH1D *cfit = new TH1D(*cdatarfmp[0]);
  // cfit->SetName("databg");

  // double ksv, bgv;
  // for (int ib=1; ib<=cdatarfmp[0]->GetNbinsX(); ib++)
  //   {
  //     val = cdatarfmp[0]->GetBinContent(ib);
  //     ksv = cdatarfmp[0]->GetXaxis()->GetBinCenter(ib);
  //     bgv = nrmvals[0] + lbdvals[0]*TMath::Exp(-ksv*ksv*wdtvals[0]*wdtvals[0]);
      
  //     val += -bgv + 1.0;
  //     val = (val-1.0)/purval + 1.0;

  //     cfit->SetBinContent(ib, val);
  //   }

  TFile *outfile = new TFile("outchicentall.root","RECREATE");
  outfile->cd();
  //  grchi->Write();

  //  cdata->Write();

  // cfs00sbgs[minsz*RMAX*IMAX + minir*IMAX + minii]->SetName("chiminbg");
  // cfs00sbgs[minsz*RMAX*IMAX + minir*IMAX + minii]->Write();

  // cfs00sbgsmr[minsz*RMAX*IMAX + minir*IMAX + minii]->SetName("chiminbgmr");
  // cfs00sbgsmr[minsz*RMAX*IMAX + minir*IMAX + minii]->Write();

  // cfs00s[minsz*RMAX*IMAX + minir*IMAX + minii]->SetName("chimin");
  // cfs00s[minsz*RMAX*IMAX + minir*IMAX + minii]->Write();

  //  cfit->Write();
  
  
  cout << "Minimum chi " << minchi << endl;

  cout << "Re Im Size " << rvals[minir] << " " << ivals[minii] << " ";

  for (int icent=0; icent<6; icent++)
    if (docent[icent])
      cout << "Centrality " << icent << "  " << svals[minszs[icent]] << endl;

  if (minir == RMAX-1) minir--;
  if (minii == IMAX-1) minii--;
  
  // if (minir3 == RMAX-1) minir3--;
  // if (minii3 == IMAX-1) minii3--;
  // if (minsz3 == SMAX-1) minsz3--;

  for (int icent=0; icent<6; icent++)
    if (minszs[icent] == SMAX-1) minszs[icent]--;

  
  double vals[500];
  //  getlinext(-0.32, 0.55, 3.5, vals);

  getlinext(-0.7, 0.3, 3.5, vals);
  TGraph *grr1 = new TGraph(500, cfs00s[0]->GetX(), vals);
  grr1->SetName("grr1");
  
  getlinext(-0.5, 0.3, 4.0, vals);
  TGraph *grr2 = new TGraph(500, cfs00s[0]->GetX(), vals);
  grr2->SetName("grr2");
  
  getlinext(-0.6, 0.3, 3.70, vals);
  TGraph *grr3 = new TGraph(500, cfs00s[0]->GetX(), vals);
  grr3->SetName("grr3");

  //  exit(0);

  cout << "Now fitting " << endl;
  
  // Values from Therminator2 fits
  // lbdval = 0.00555;
  // wdtval = 1.92;
  
  cout << "Setting values " << nrmvals[0] << " " << lbdvals[0] << " " << wdtvals[0] << " " << purval << endl;
  
  TVirtualFitter *fitter=TVirtualFitter::Fitter(0, 45);
  fitter->SetFCN(myfuncf);
  fitter->SetParameter(0,  "Purity"  ,purval        ,0.0001,0.2,   0.85);
  fitter->SetParameter(1,  "Ref0"    ,rvals[minir]  ,0.0001,-5.0, 5.0);
  fitter->SetParameter(2,  "Imf0"    ,ivals[minii]  ,0.0001,0.00001, 5.0);
  // fitter->SetParameter(2,  "Imf0"    ,0.2 ,0.0001,0.00001, 5.0);
  // fitter->FixParameter(2);
  
  
  for (int icent=0; icent<6; icent++) {
    fitter->SetParameter(3+icent*3+0,  Form("Lambda%i", icent)  ,lbdvals[icent]    ,0.001, 0.0,   0.1);
    fitter->SetParameter(3+icent*3+1,  Form("Width%i", icent)   ,wdtvals[icent]    ,0.001, 1.0,   10.0);
    fitter->SetParameter(3+icent*3+2,  Form("Radius%i", icent)  ,svals[minszs[icent]]  ,0.01  ,0.9, 8.0);

    if (!docent[icent]) {
      fitter->FixParameter(3+icent*3+0);
      fitter->FixParameter(3+icent*3+1);
      fitter->FixParameter(3+icent*3+2);
    }

    fitter->FixParameter(3+icent*3+1);
  }

  for (int icent=0; icent<6; icent++) {
    for (int ipair=0; ipair<4; ipair++) {
      fitter->SetParameter(21+icent*4+ipair,  Form("Norm%i%i", icent, ipair)    ,nrmvals[icent]    ,0.001, 0.9,   1.1);
      if (!docent[icent])
	fitter->FixParameter(21+icent*4+ipair);
    }
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

  for (int icent=0; icent<6; icent++) {
    fitlabds[icent] = fitter->GetParameter(3+icent*3+0);
    fitwdths[icent] = fitter->GetParameter(3+icent*3+1);
    fitradis[icent] = fitter->GetParameter(3+icent*3+2);

    for (int ipair=0; ipair<4; ipair++) {
      fitnrmscales[icent*4+ipair] = fitter->GetParameter(21+icent*4+ipair);
    }
  }


  
  // fitnorms[0] = fitnorm;
  // fitnorms[3] = fitnorm3;
  // fitlabds[0] = fitlabd;
  // fitlabds[3] = fitlabd3;
  // fitwdths[0] = fitwdth;
  // fitwdths[3] = fitwdth3;
  // fitradis[0] = fitrad0;
  // fitradis[3] = fitrad3;

  // fitnrmscales[0*4+0] = fitnorm;
  // fitnrmscales[0*4+1] = fitnorm*fitnrm1;
  // fitnrmscales[0*4+2] = fitnorm*fitnrm2;
  // fitnrmscales[0*4+3] = fitnorm*fitnrm3;
  
  // fitnrmscales[3*4+0] = fitnorm3;
  // fitnrmscales[3*4+1] = fitnorm3*fitnrm4;
  // fitnrmscales[3*4+2] = fitnorm3*fitnrm5;
  // fitnrmscales[3*4+3] = fitnorm3*fitnrm6;
  
  cout << "Final paramerers:" << endl;
  cout << "Purity          : " << fitpuri << endl;
  cout << "Re f0           : " << fitref0 << endl;
  cout << "Im f0           : " << fitimf0 << endl;

  for (int icent=0; icent<6; icent++) {
    cout << "Radius C" << icent << "     :" << fitradis[icent] << endl;
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
  }

  
  
  // cout << "Normalization   : " << fitnorm << endl;
  // cout << "BG lambda       : " << fitlabd << endl;
  // cout << "BG width        : " << fitwdth << endl;
  // cout << "Purity          : " << fitpuri << endl;
  // cout << "Re f0           : " << fitref0 << endl;
  // cout << "Im f0           : " << fitimf0 << endl;
  // cout << "Radius          : " << fitrad0 << endl;
  // cout << "Norm scale 1    : " << fitnrm1 << endl;
  // cout << "Norm scale 2    : " << fitnrm2 << endl;
  // cout << "Norm scale 3    : " << fitnrm3 << endl;
  
  // cout << "Normalization 3 : " << fitnorm3 << endl;
  // cout << "BG lambda 3     : " << fitlabd3 << endl;
  // cout << "BG width 3      : " << fitwdth3 << endl;
  // cout << "Radius 3        : " << fitrad3 << endl;
  // cout << "Norm scale 1    : " << fitnrm4 << endl;
  // cout << "Norm scale 2    : " << fitnrm5 << endl;
  // cout << "Norm scale 3    : " << fitnrm6 << endl;
  
  outfile->cd();

  TGraph *grcfits[6];
  
  for (int icent=0; icent<6; icent++) {
    if (docent[icent]) {
      
      getlinext(fitref0, fitimf0, fitradis[icent], vals);
      
      grcfits[icent] = new TGraph(500, cfs00s[0]->GetX(), vals);
      grcfits[icent]->SetTitle(";k^ (GeV/#it{c});C(k*)");
      grcfits[icent]->SetName(Form("grappfitC%i", icent));
      
      grcfits[icent]->Write();
    }
  }
  
  // getlinext(fitref0, fitimf0, fitrad3, vals);
  
  // TGraph *grfitc3 = new TGraph(500, cfs00s[0]->GetX(), vals);
  // grfitc3->SetTitle(";k^ (GeV/#it{c});C(k*)");
  // grfitc3->SetName("grappfitC3");

  // grfitc3->Write();


  cfs00ppbg->SetName("cfs00ppbg");
  cfs00ppbg->Write();
  cfs00ppbgmr->SetName(Form("%s_mr",cfs00ppbg->GetName()));
  cfs00ppbgmr->Write();

  TH1D *cfitrfpmlin[6];
  TH1D *cfitrfmplin[6];
  TH1D *cfitffpmlin[6];
  TH1D *cfitffmplin[6];

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
    }
  }
  
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

  grr1->Write();
  grr2->Write();
  grr3->Write();
  
}
 
