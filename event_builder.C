#include <fstream>
#include <iostream>
#include <iomanip>

#include "TBranch.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2S.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

TString FILENAME = "mc_fadc.root";

//  double tau = 500.; // 2.0 MHz
//  double tau = 10000.; // 0.1 MHz
double TAU = 1000.; // 1.0 MHz
int PROCESS_EV = 100;
int EVERY      = 10;

double R1=10; // mm
double R2=30; // mm
double R3=50; // mm

double W1 = 2.0*0.001; // mm/ns
double W2 = 4.0*0.001; // mm/ns
double z_anod = -50.;

double Digi[125];

TTree *tree;
TTree *mc_true;
TFile *file;

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

  int EVENT = 0;

  double startTPC = 40000. ;

  h1 = new TH1F("h1"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h2 = new TH1F("h2"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h3 = new TH1F("h3"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );

  tree = new TTree("tree", "fadc_tree");
  TBranch *br1 = tree->Branch("h1", "TH1F", &h1, 640000, 0);
  TBranch *br2 = tree->Branch("h2", "TH1F", &h2, 640000, 0);
  TBranch *br3 = tree->Branch("h3", "TH1F", &h3, 640000, 0);

  int    ev, tr ;
  long int code;
  double ed;
  double x ,y ,z ,t ;
  double xi,yi,zi,ti;
  double xf,yf,zf,tf;
  double t_anod, tt, ll;




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

  std::ifstream frGEN( "./gen.data"          , std::ios::in ); // generator data
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

  while(n_ev<PROCESS_EV){

    frTPC >> ev >> tr >> code >> ed >> xi >> yi >> zi >> ti  >> xf >> yf >> zf >> tf;

    if(ev!=EVENT){
      if( !(n_ev%EVERY) ) std::cout << n_ev << " event processed\n" ;
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

        beam_offset += rE->Exp( TAU );
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

      tree->Fill();

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


  int mc_ev, mc_id;
  float mc_TR, mc_theta, mc_vx, mc_vy,mc_vz;

  mc_true = new TTree("mc_true", "mc_true");
  TBranch *br_ev    = mc_true->Branch("ev", &mc_ev, "ev/I");
  TBranch *br_id    = mc_true->Branch("id", &mc_id, "id/I");
  TBranch *br_TR    = mc_true->Branch("T", &mc_TR, "T/F");
  TBranch *br_Theta = mc_true->Branch("theta", &mc_theta, "theta/F");
  TBranch *br_vx    = mc_true->Branch("vx", &mc_vx, "vx/F");
  TBranch *br_vy    = mc_true->Branch("vy", &mc_vy, "vy/F");
  TBranch *br_vz    = mc_true->Branch("vz", &mc_vz, "vz/F");

  while( frGEN >> mc_ev >> mc_id >> mc_TR >> mc_theta >> mc_vx >> mc_vy >> mc_vz){
    mc_true->Fill();
  } frGEN.close();


  file = new TFile(FILENAME,"RECREATE");
  tree->Write();
  mc_true->Write();
  file->Write();
  file->Close();

}



void event_builder(){

  findTPCtracks();

  gSystem->Exit(0);

}
