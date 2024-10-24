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

bool PRINT_RECO=false;

//  double tau = 500.; // 2.0 MHz
//  double tau = 10000.; // 0.1 MHz
double TAU = 1000.; // 1.0 MHz
int PROCESS_EV = 500;
int EVERY      = 10;

double R1=10; // mm
double R2=30; // mm
double R3=50; // mm

double W1 = 2.0*0.001; // mm/ns
double W2 = 4.0*0.001; // mm/ns
double z_anod = -50.;

double Digi[125];

TTree *tree;
TTree *reco;
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

struct cv
{
    UInt_t t; // time channel
    double v; // value
};


struct signal
{
    UInt_t status;
    UInt_t start_bin ; // begin of signal
    UInt_t peak_bin  ; // begin of maximum
    UInt_t end_bin   ; // end of signal
    UInt_t length_bin;
    double energy    ;
    double peak      ;
    double start     ;
    double end       ;
    double base      ;
};

  struct signal s1;
  struct signal s2;
  struct signal s3;
  struct signal s ;


double avg( TH1F *h ){
  double sum = 0;
  for (int t = 1; t < h->GetNbinsX()+1; t++) sum += h->GetBinContent(t);
  return ( sum / h->GetNbinsX() );
}


const int    WINDOW    = 10    ;
const double THRESHOLD =  0.1  ;

UInt_t check_fired( TH1F *h ){

  int TRACELENGTH = h->GetNbinsX();
  double av = avg(h);
  cv pmin,pmax;
  pmin.t = h->GetMinimumBin() ;
  pmax.t = h->GetMaximumBin() ;
  pmin.v = h->GetMinimum()    ;
  pmax.v = h->GetMaximum()    ;

  UInt_t isf = 0;

  double running_sum=0;
  for(int t=0; t<WINDOW; t++)
    running_sum += ( h->GetBinContent(t+1)*(1./WINDOW) );

  if( av-pmin.v <= pmax.v-av ){
    for(int t=0; t < TRACELENGTH - WINDOW ; t++ ){
      if( running_sum > av*(1.+THRESHOLD) ) isf = 1 ;
      running_sum -= ( h->GetBinContent( t + 1              )  / WINDOW );
      running_sum += ( h->GetBinContent( t + 1 + WINDOW + 1 )  / WINDOW );
    }
  return isf;
  }

  for(int t=0; t < TRACELENGTH - WINDOW ; t++ ){
    if( running_sum < av*(1.-THRESHOLD) ) isf = 2 ;
    running_sum -= ( h->GetBinContent( t + 1              )  / WINDOW );
    running_sum += ( h->GetBinContent( t + 1 + WINDOW + 1 )  / WINDOW );
  }
  return isf;

}

const int    MARGIN_LEFT  = 150 ;
const int    MARGIN_RIGHT = 150 ;

double eval_base( TH1F *h ){

  int TRACELENGTH = h->GetNbinsX();
  UInt_t status = check_fired( h );

  double base = 0;
  int counted_bins = 0;
  Int_t pos;

  if( status == 0 ) return avg(h);
  if( status == 1 ) pos = h->GetMaximumBin();
  else              pos = h->GetMinimumBin();

  int bin_low  = pos - 1 - MARGIN_LEFT  ;
  int bin_high = pos - 1 + MARGIN_RIGHT ;

  for( int t=0; t<TRACELENGTH; t++ ){
    if( (t < bin_low) || (t > bin_high) ){
      base += h->GetBinContent(t+1);
      counted_bins++;
    }
  }

 return base/counted_bins;

}



bool get_s( TH1F *h ){

  s.start_bin=0; s.end_bin=0; s.peak=0; s.start=0; s.end=0; s.energy=0;
  s.status = check_fired( h );

  s.base = eval_base( h );
  double base = s.base;

  if( s.status!=1 ) return false;

  Int_t pos = h->GetMaximumBin(); s.peak_bin=pos; s.peak = h->GetMaximum();

  s.energy = s.peak - base;
  //left
  int t=pos-1;
  double height = s.peak - base;
  cv p80, p25; p80.t=0; p25.t=0; p80.v=s.peak; p25.v = base ;
  while( h->GetBinContent(t)>base ){
     s.energy += ( h->GetBinContent(t)-base );
     if( ( h->GetBinContent(t)-base ) / height < 0.80 && p80.t==0 ){
         p80.t = t ; p80.v = h->GetBinContent(t);
     }
     if( (h->GetBinContent(t)-base ) / height < 0.25 && p25.t==0 ){
         p25.t = t ; p25.v = h->GetBinContent(t);
     }
     t--;
  } s.start_bin = t+1;

  s.start = p25.t - double(p80.t-p25.t)*double(p25.v-UInt_t( base) )/double(p80.v-p25.v);

  //right
  t=pos+1;
  while( h->GetBinContent(t)>base ){
     s.energy += ( h->GetBinContent(t)-base );
     t++;
  } s.end_bin = t-1; s.end = s.end_bin;

  s.length_bin = s.end_bin - s.start_bin ;

  return true;
}


bool get_sr( TH1F *h , int&  start_bin, int& end_bin, int& peak_bin, double& peak, double& start, double& end, double& energy, UInt_t& status, double& base){

  start_bin=0; end_bin=0; peak_bin = 0;
  peak=0; start=0; end=0; energy=0;
  status = check_fired( h ); 
  base = eval_base( h );

  if( status!=1 ) return false;

  Int_t pos = h->GetMaximumBin(); peak_bin=pos; peak = h->GetMaximum();

  energy = peak - base;
  //left
  int t=pos-1;
  double height = peak - base;
  cv p80, p25; p80.t=0; p25.t=0; p80.v=peak; p25.v = base ;
  while( h->GetBinContent(t)>base && t>0){
     energy += ( h->GetBinContent(t)-base );
     if( ( h->GetBinContent(t)-base ) / height < 0.80 && p80.t==0 ){
         p80.t = t ; p80.v = h->GetBinContent(t);
     }
     if( (h->GetBinContent(t)-base ) / height < 0.25 && p25.t==0 ){
         p25.t = t ; p25.v = h->GetBinContent(t);
     }
     t--;
  } start_bin = t+1;

  start = p25.t - double(p80.t-p25.t)*double(p25.v- base )/double(p80.v-p25.v);

  //right
  t=pos+1;
  while( h->GetBinContent(t)>base && t<2693){
     energy += ( h->GetBinContent(t)-base );
     t++;
  } end_bin = t-1; end = end_bin;


  return true;
}




void findTPCtracks(){


  TRandom *rE = new TRandom();

  int EVENT = 0;

  double startTPC = 40000. ;

  h1 = new TH1F("h1"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h2 = new TH1F("h2"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );
  h3 = new TH1F("h3"," ;time, 10*ns; energy, a.u.", 2692, 0., 4.*2692. );

  tree = new TTree("tree", "fadc_tree");
  TBranch *br1 = tree->Branch("h1", "TH1F", &h1, 1280000, 0);
  TBranch *br2 = tree->Branch("h2", "TH1F", &h2, 1280000, 0);
  TBranch *br3 = tree->Branch("h3", "TH1F", &h3, 1280000, 0);

  double ss1,ss2,ss3;
  int    ev, tr ;
  long int code;
  double ed;
  double x ,y ,z ,t ;
  double xi,yi,zi,ti;
  double xf,yf,zf,tf;
  double t_anod, tt, ll;


  int start_bin[3], end_bin[3], peak_bin[3];
  double peak[3], start[3], end[3], energy[3], base[3], s_length[3];
  double energy_tot;
  UInt_t status[3];

  reco = new TTree("reco", "reco_tree");
  TBranch *bR_s1     = reco->Branch("status1", &status[0], "status1/I");
  TBranch *bR_s2     = reco->Branch("status2", &status[1], "status2/I");
  TBranch *bR_s3     = reco->Branch("status3", &status[2], "status3/I");
  TBranch *bR_E1     = reco->Branch("E1", &energy[0], "E1/D");
  TBranch *bR_E2     = reco->Branch("E2", &energy[1], "E2/D");
  TBranch *bR_E3     = reco->Branch("E3", &energy[2], "E3/D");
  TBranch *bR_ER     = reco->Branch("ER", &energy_tot, "ER/D");
  TBranch *bR_start1 = reco->Branch("start1", &start[0], "start1/D");
  TBranch *bR_start2 = reco->Branch("start2", &start[1], "start2/D");
  TBranch *bR_start3 = reco->Branch("start3", &start[2], "start3/D");
  TBranch *bR_stop1  = reco->Branch("stop1", &end[0],  "stop1/D");
  TBranch *bR_stop2  = reco->Branch("stop2", &end[1],  "stop2/D");
  TBranch *bR_stop3  = reco->Branch("stop3", &end[2],  "stop3/D");
  TBranch *bR_peak1  = reco->Branch("peak1", &peak[0], "peak1/D");
  TBranch *bR_peak2  = reco->Branch("peak2", &peak[1], "peak2/D");
  TBranch *bR_peak3  = reco->Branch("peak3", &peak[2], "peak3/D");
  TBranch *bR_base1  = reco->Branch("base1", &base[0], "base1/D");
  TBranch *bR_base2  = reco->Branch("base2", &base[1], "base2/D");
  TBranch *bR_base3  = reco->Branch("base3", &base[2], "base3/D");



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

      ss1 = get_sr(h1,start_bin[0], end_bin[0], peak_bin[0], peak[0], start[0], end[0], energy[0], status[0], base[0] );
      ss2 = get_sr(h2,start_bin[1], end_bin[1], peak_bin[1], peak[1], start[1], end[1], energy[1], status[1], base[1] );
      ss3 = get_sr(h3,start_bin[2], end_bin[2], peak_bin[2], peak[2], start[2], end[2], energy[2], status[2], base[2] );

      energy_tot = energy[0]+energy[1]+energy[2];
      reco->Fill();
      if(PRINT_RECO){
        std::cout << start_bin[0] << "\t"  << end_bin[0] << "\t" << peak_bin[0]  << "\t" << start[0] << "\t" << end[0] << " "<< peak[0] << "\t"<<energy[0] << "\t"<< status[0] << "\t"<< base[0] << "\n";
        std::cout << start_bin[1] << "\t"  << end_bin[1] << "\t" << peak_bin[1]  << "\t" << start[1] << "\t" << end[1] << " "<< peak[1] << "\t"<<energy[1] << "\t"<< status[1] << "\t"<< base[1] << "\n";
        std::cout << start_bin[2] << "\t"  << end_bin[2] << "\t" << peak_bin[2]  << "\t" << start[2] << "\t" << end[2] << " "<< peak[2] << "\t"<<energy[2] << "\t"<< status[2] << "\t"<< base[2] << "\n\n";
      }

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
  reco->Write();
  mc_true->Write();
  file->Write();
  file->Close();

}



void event_builder(){

  findTPCtracks();

  gSystem->Exit(0);

}
