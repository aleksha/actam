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


void findTPCtracks(){

  int process_ev = 999;
  int EVENT = 0;

  double startTPC = 0. ;

  TH1F* h1 = new TH1F("h1"," ;time, 10*ns; energy, a.u.", 2550, 0., 4.*2550. );
  TH1F* h2 = new TH1F("h2"," ;time, 10*ns; energy, a.u.", 2550, 0., 4.*2550. );
  TH1F* h3 = new TH1F("h3"," ;time, 10*ns; energy, a.u.", 2550, 0., 4.*2550. );

  TH1F* hE1 = new TH1F("hE1"," ; energy, a.u.;entries",  50, 0., 5 );
  TH1F* hE2 = new TH1F("hE2"," ; energy, a.u.;entries",  50, 0., 5 );
  TH1F* hE3 = new TH1F("hE3"," ; energy, a.u.;entries",  50, 0., 5 );


  h1->SetMinimum(0);
  h2->SetMinimum(0);
  h3->SetMinimum(0);

  int    ev, tr ;
  long int code;
  double ed;
  double x,y,z;
  double t_anod, tt, ll;

  double Digi[125];


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

  std::ifstream frTPC("./out.data"      , std::ios::in);

  int n_ev = 0;

  ev = 0;
  while(n_ev<process_ev){

    frTPC >> ev >> tr >> code >> ed >> x >> y >> z ;


    if(ev!=EVENT){
      n_ev++;

      hE1->Fill( h1->Integral() );
      hE2->Fill( h2->Integral() );
      hE3->Fill( h3->Integral() );

      h1->Reset();
      h2->Reset();
      h3->Reset();
      EVENT=ev;
    }

    t_anod = 0.1*( (z-z_anod) / W1 + 3./W2 );
    tt = 2;
    ll = x*x+y*y;

    if(ll<10.*10.){
      for(int iii = 0 ; iii<125; iii++  ){
        h1->Fill( t_anod + tt, ed*Digi[iii] );
        tt = tt + 4 ;
      }
    }

    if(ll>=10.*10. && ll<30.*30.){
      for(int iii = 0 ; iii<125; iii++  ){
        h2->Fill( t_anod + tt, ed*Digi[iii] );
        tt = tt + 4 ;
      }
    }

    if(ll>=30.*30.){
      for(int iii = 0 ; iii<125; iii++  ){
        h3->Fill( t_anod + tt, ed*Digi[iii] );
        tt = tt + 4 ;
      }
    }


  }  frTPC.close();

  TCanvas* canv = new TCanvas("canv","canv",900,600);
  canv->Divide(1,3);
  canv->cd(1);
  hE1->Draw("hist");
  canv->cd(2);
  hE2->Draw("hist");
  canv->cd(3);
  hE3->Draw("hist");

  canv->Print("c.png");
//  canv->Close();

//  hTime->Draw();


}



void analyze(){


  findTPCtracks();

//  gSystem->Exit(0);

}
