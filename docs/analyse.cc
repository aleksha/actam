struct cv
{
    UInt_t t; // time channel
    UInt_t v; // value
};

struct signal
{
    UInt_t status;
    UInt_t start_bin ; // begin of signal
    UInt_t end_bin   ; // end of signal
    UInt_t length_bin;
    double energy    ;
    double peak      ;
    double start     ;
    double end       ;
    double noise[ NOISESAMPLES ];
    double nback[ NOISEBACK    ];
};


double av( UInt_t h[][TRACELENGTH] , Int_t ch){
  UInt_t sum = 0;
  for (int t = 0; t < TRACELENGTH; t++) sum += h[ch][t];
  return ( double(sum) / TRACELENGTH );
}

cv min( UInt_t h[][TRACELENGTH] , Int_t ch){
  cv p; p.t = 0; p.v=20000;
  for (int t = 0; t < TRACELENGTH; t++){
    if( h[ch][t] < p.v){ p.v=h[ch][t]; p.t=t;}
  }
  return p;
}

cv max( UInt_t h[][TRACELENGTH] , Int_t ch){
  cv p; p.t = 0; p.v=0;
  for (int t = 0; t < TRACELENGTH; t++){
    if( h[ch][t] > p.v){ p.v=h[ch][t]; p.t=t;}
  }
  return p;
}

UInt_t check_fired( UInt_t h[][TRACELENGTH] , UInt_t ch, double av, cv pmin, cv pmax){

  UInt_t isf = 0;

  double running_sum=0;
  for(int t=0; t<WINDOW; t++)
    running_sum += ( h[ch][t]*(1./WINDOW) );

  if( av-pmin.v <= pmax.v-av ){
    for(int t=0; t < TRACELENGTH - WINDOW ; t++ ){
      if( running_sum > av*(1.+THRESHOLD) ) isf = 1 ;
      running_sum -= ( h[ch][ t ]               / WINDOW );
      running_sum += ( h[ch][ t + WINDOW + 1 ]  / WINDOW );
    }
  return isf;
  }

  for(int t=0; t < TRACELENGTH - WINDOW ; t++ ){
    if( running_sum < av*(1.-THRESHOLD) ) isf = 2 ;
    running_sum -= ( h[ch][ t ]               / WINDOW );
    running_sum += ( h[ch][ t + WINDOW + 1 ]  / WINDOW );
  }
  return isf;

}

double eval_base(  UInt_t h[][TRACELENGTH] , UInt_t ch, double av, cv pmin, cv pmax, UInt_t status, bool start_only = false ){
  double base = 0;
  int counted_bins = 0;
  Int_t pos;

  if( !start_only ){
    if( status == 0 ) return av;
    if( status == 1 ) pos=pmax.t;
    else              pos=pmin.t;
    int bin_low  = pos - MARGIN_LEFT  ;
    int bin_high = pos + MARGIN_RIGHT ;

    for( int t=0; t<TRACELENGTH; t++ ){
      if( (t < bin_low) || (t > bin_high) ){
        base += h[ch][t];
        counted_bins++;
      }
    }
  return base/counted_bins;
  }

  for(Int_t t=0; t<(NOISESAMPLES*NOISE_LENGTH); t++) base += h[ch][t];

  return base/(NOISESAMPLES*NOISE_LENGTH);
  
}

struct signal get_s(  UInt_t h[][TRACELENGTH] , UInt_t ch, double base, cv pmin, cv pmax, UInt_t status ){

  struct signal s; s.start_bin=0; s.end_bin=0; s.peak=0; s.start=0; s.end=0; s.energy=0;
  s.status = status;

  if( status==0 ) return s;

  Int_t pos;
  if( status == 1 ) pos=pmax.t;
  else              pos=pmin.t;

  s.peak   = pos;
  s.energy = h[ch][pos]-base;
  //left
  int t=pos-1;
  double height = h[ch][t] - base;
  cv p80, p25; p80.t=0; p25.t=0; p80.v=h[ch][t]; p25.v = UInt_t( base );
  while( (status==1 && h[ch][t]>base) || (status==2 && h[ch][t]<base) ){
     s.energy += ( h[ch][t]-base );
     if( double(h[ch][t]-base ) / height < 0.80 && p80.t==0 ){
         p80.t = t ; p80.v = h[ch][t];
     }
     if( double(h[ch][t]-base ) / height < 0.25 && p25.t==0 ){
         p25.t = t ; p25.v = h[ch][t];
     }
     t--;
  } s.start_bin = t+1; 
  if(status==1)
    s.start = p25.t - double(p80.t-p25.t)*double(p25.v-UInt_t( base) )/double(p80.v-p25.v);

  //right
  t=pos+1;
  while( (status==1 && h[ch][t]>base) || (status==2 && h[ch][t]<base) ){
     s.energy += ( h[ch][t]-base );
     t++;
  } s.end_bin = t-1; s.end = s.end_bin;

  s.length_bin = s.end_bin - s.start_bin ;

  for(Int_t ns=0; ns<NOISESAMPLES; ns++){
    s.noise[ns]=0;
    for(pos=ns*NOISE_LENGTH;pos<(ns+1)*NOISE_LENGTH;pos++)
      s.noise[ns] += ( h[ch][pos]-base );
  }

  for(Int_t ns=0; ns<NOISEBACK; ns++){
    s.nback[ns]=0;
    for(pos=TRACELENGTH-ns*NOISE_LENGTH;pos>( TRACELENGTH - (ns+1)*NOISE_LENGTH);pos--)
      s.nback[ns] += ( h[ch][pos]-base );
  }

  return s;
}

struct signal get_n(  UInt_t h[][TRACELENGTH] , UInt_t ch, double base, cv pmin, cv pmax, UInt_t status ){

  struct signal s; s.start_bin=0; s.end_bin=0; s.peak=0; s.start=0; s.end=0; s.energy=0;
  s.status = status;

  Int_t pos;  pos=1900;  s.peak   = pos;

  for(Int_t ns=0; ns<NOISESAMPLES; ns++){
    s.noise[ns]=0;
    for(pos=ns*NOISE_LENGTH;pos<(ns+1)*NOISE_LENGTH;pos++)
      s.noise[ns] += ( h[ch][pos]-base );
  }

  for(Int_t ns=0; ns<NOISEBACK; ns++){
    s.nback[ns]=0;
    for(pos=TRACELENGTH-ns*NOISE_LENGTH;pos>( TRACELENGTH - (ns+1)*NOISE_LENGTH);pos--)
      s.nback[ns] += ( h[ch][pos]-base );
  }

  return s;
}



int sasha2status( int s ){
  if(s<16)         return s-1;
  if(s>15 && s<31) return s  ;
  if(s>30 && s<46) return s+1;
  if(s>45 && s<61) return s+2;
  return 63;
}

int status2sasha( int s ){
  if(s<15)         return s+1;
//  if(s==15) return 0;
  if(s>15 && s<31) return s  ;
//  if(s==31) return 0;
  if(s>31 && s<47) return s-1;
//  if(s==45) return 0;
  if(s>47 && s<63) return s-2;
  return 0;
}

int get_good_fired(){
  int n_fired = 0;
  for(int ch=0;ch<NOFCH;ch++){
    h[ch]->SetMaximum( MaxInEvent ) ;
     if(StatusArray[ch]==1 && GoodChannel[ch]){
       h[ch]->SetLineColor( n_fired+1 );
       CrownColor[ status2sasha(ch) ] = n_fired+1;
       //cout << ch << "\t" << status2sasha(ch) << "\t" << CrownColor[ status2sasha(ch) ] << "\n";
       n_fired++;
     }
  }
  return n_fired;
}

/*
int status2sasha( int s ){
  if(s<16)         return s+1;
  if(s>15 && s<31) return s  ;
  if(s>30 && s<45) return s-1;
  if(s>44 && s<61) return s-2;
  return 63;
}
*/

