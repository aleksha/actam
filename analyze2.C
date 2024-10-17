#include <fstream>
#include <iostream>
#include <iomanip>

#include "TH1F.h"
#include "TH1S.h"
#include "TH2S.h"
#include "TF1.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TLatex.h"

double Digi[125];

double R1=10;
double R2=30;
double R3=50;

TH1F *h1;
TH1F *h2;
TH1F *h3;

void fillFADC(double E_step, double t_anod, double ll){

  double tt=2.;

  if(ll<R1*R1){
    for(int iii = 0 ; iii<125; iii++  ){
      h1->Fill( t_anod + tt , E_step*Digi[iii] );
      tt = tt + 4 ;
    }
  }

  if(ll>=R1*R1 && ll<R2*R2){
    for(int iii = 0 ; iii<125; iii++  ){
      h2->Fill( t_anod + tt, E_step*Digi[iii] );
      tt = tt + 4 ;
      }
  }

  if(ll>=R2*R2 && ll<R3*R3){
    for(int iii = 0 ; iii<125; iii++  ){
      h3->Fill( t_anod + tt, E_step*Digi[iii] );
      tt = tt + 4 ;
    }
  }

}

void findTPCtracks(){

  TRandom *rE = new TRandom();

//  double tau = 500.; // 2.0 MHz
//  double tau = 10000.; // 0.1 MHz
  double tau = 1000.; // 1.0 MHz

  int process_ev = 10;
  int EVENT = 0;

  double startTPC = 40000. ;

  h1 = new TH1F("h1"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h2 = new TH1F("h2"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h3 = new TH1F("h3"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );

  TH1F* h1p;
  TH1F* h2p;
  TH1F* h3p;


  TH1F* hE1 = new TH1F("hE1"," ; energy, a.u.;entries",  100, 0.1, 10.1 );
  TH1F* hE2 = new TH1F("hE2"," ; energy, a.u.;entries",  100, 0.1, 10.1 );
  TH1F* hE3 = new TH1F("hE3"," ; energy, a.u.;entries",  100, 0.1, 10.1 );

//  h1->SetMinimum(0);
//  h2->SetMinimum(0);
//  h3->SetMinimum(0);

  int    ev, tr ;
  long int code;
  double ed;
  double x ,y ,z ,t ;
  double xi,yi,zi,ti;
  double xf,yf,zf,tf;
  double t_anod, tt, ll;



  double W1 = 2.0*0.001; // mm/ns
  double W2 = 4.0*0.001; // mm/ns
  double z_anod = -50.;

//==============================================================================
// create digitization
//==============================================================================
  TF1 *f_e1 = new TF1("f_e1"," 0.02207760*exp(-3.64674*0.673286*x)"                        ,0,5);
  TF1 *f_e2 = new TF1("f_e2","-0.02525950*exp(-3.35196*0.673286*x)*cos(1.74266*0.673286*x)",0,5);
  TF1 *f_e3 = new TF1("f_e3"," 0.00318191*exp(-2.32467*0.673286*x)*cos(3.57102*0.673286*x)",0,5);
  TF1 *f_e4 = new TF1("f_e4"," 0.01342450*exp(-3.35196*0.673286*x)*sin(1.74266*0.673286*x)",0,5);
  TF1 *f_e5 = new TF1("f_e5","-0.00564406*exp(-2.32467*0.673286*x)*sin(3.57102*0.673286*x)",0,5);
  TF1 *ft   = new TF1("ft"  ,"636.252*(f_e1 + f_e2 + f_e3 + f_e4 + f_e5)"                  ,0,5);

  tt = 0.02;
  double DS = 0 ;
  for(int ii=0; ii<125; ii++){
    Digi[ii] = 0.04*ft->Eval(tt);
    DS = DS + Digi[ii];
    tt = tt + 0.04;
  }



//==============================================================================
// loop on trigger event
//==============================================================================

  std::ifstream frTPC( "./out.data"          , std::ios::in ); // signal from recoil
  std::ifstream fBEAM( "./out.data.beam"     , std::ios::in ); // pile-up from the beam
  std::ifstream fNOIS( "./noise_events.data" , std::ios::in ); // electronic noise

  int n_ev = 0;
  float E_step = 0.001;
  int n_steps;

  float beam_offset = -100000.;
  int beam_ev=0;

  ev = 0;

  int ev_b;
  float xib,yib,zib,tib;
  float xfb,yfb,zfb,tfb,ed_b;
  int tr_b,code_b;
  float noise;

  fBEAM >> ev_b >> tr_b >> code_b >> ed_b >> xib >> yib >> zib >> tib  >> xfb >> yfb >> zfb >> tfb;

  while(n_ev<process_ev){

    frTPC >> ev >> tr >> code >> ed >> xi >> yi >> zi >> ti  >> xf >> yf >> zf >> tf;

    if(ev!=EVENT){
      if( !(n_ev%1) ) std::cout << n_ev << " event processed\n" ;
      n_ev++;


      while(beam_offset<100000){

        while(beam_ev==ev_b){

          n_steps = int(ed_b/E_step)+1;
          for(int step=0; step < n_steps; step++){

            x = xib + (xfb-xib)*(0.5+step)/n_steps;
            y = yib + (yfb-yib)*(0.5+step)/n_steps;
            z = zib + (zfb-zib)*(0.5+step)/n_steps;
            t = tib + (tfb-tib)*(0.5+step)/n_steps;

            t_anod = 0.1*( beam_offset + startTPC + (z-z_anod) / W1 + 3./W2 );
            ll = x*x+y*y;

            fillFADC(E_step, t_anod, ll);

          }

          fBEAM >> ev_b >> tr_b >> code_b >> ed_b >> xib >> yib >> zib >> tib  >> xfb >> yfb >> zfb >> tfb;
        }

        beam_offset += rE->Exp( tau );
        beam_ev = ev_b;

      }

      beam_offset = -100000.;



      for(int ch=0;ch<2692;ch++){
        fNOIS >> noise;
        h1->SetBinContent( ch+1, 27.5*0.001*0.001*noise+h1->GetBinContent(ch+1) );
      }

      for(int ch=0;ch<2692;ch++){
        fNOIS >> noise;
        h2->SetBinContent( ch+1, 27.5*0.001*0.001*noise+h2->GetBinContent(ch+1) );
      }

      for(int ch=0;ch<2692;ch++){
        fNOIS >> noise;
        h3->SetBinContent( ch+1, 27.5*0.001*0.001*noise+h3->GetBinContent(ch+1) );
      }



      hE1->Fill( h1->Integral() );
      hE2->Fill( h2->Integral() );
      hE3->Fill( h3->Integral() );

//TH1F *hnew = (TH1F*)h->Clone("hnew");
      h1p = (TH1F*)h1->Clone("h1p");
      h2p = (TH1F*)h2->Clone("h2p");
      h3p = (TH1F*)h3->Clone("h3p");

      h1->Reset();
      h2->Reset();
      h3->Reset();
      EVENT=ev;
    }

    n_steps = int(ed/E_step)+1;
    for(int step=0; step < n_steps; step++){

      x = xi + (xf-xi)*(0.5+step)/n_steps;
      y = yi + (yf-yi)*(0.5+step)/n_steps;
      z = zi + (zf-zi)*(0.5+step)/n_steps;
      t = ti + (tf-ti)*(0.5+step)/n_steps;

      t_anod = 0.1*( startTPC + (z-z_anod) / W1 + 3./W2 );
      ll = x*x+y*y;

      fillFADC(E_step, t_anod, ll);

    }

  }  frTPC.close();


  double hMAX=0;
  if(h1p->GetMaximum()>hMAX) hMAX = h1p->GetMaximum();
  if(h2p->GetMaximum()>hMAX) hMAX = h2p->GetMaximum();
  if(h3p->GetMaximum()>hMAX) hMAX = h3p->GetMaximum();

  hMAX = hMAX*1.05;
  h1p->SetMaximum(hMAX);
  h2p->SetMaximum(hMAX);
  h3p->SetMaximum(hMAX);

  double hMIX=0;
  if(h1p->GetMinimum()<hMIX) hMIX = h1p->GetMinimum();
  if(h2p->GetMinimum()<hMIX) hMIX = h2p->GetMinimum();
  if(h3p->GetMinimum()<hMIX) hMIX = h3p->GetMinimum();

  hMIX = hMIX*1.05;
  h1p->SetMinimum(hMIX);
  h2p->SetMinimum(hMIX);
  h3p->SetMinimum(hMIX);



  TCanvas* canv = new TCanvas("canv","canv",900,600);
  canv->Divide(1,3);
  canv->cd(1);
  h1p->Draw("hist");
  canv->cd(2);
  h2p->Draw("hist");
  canv->cd(3);
  h3p->Draw("hist");

  canv->Print("c.png");
//  canv->Close();

//  hTime->Draw();


}



void analyze2(){


  findTPCtracks();

//  gSystem->Exit(0);

}
