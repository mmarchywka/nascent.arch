#ifndef MJM_TIMING_H__
#define MJM_TIMING_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <time.h>
/***

Mike Marchywka May 2011

*/

#include <sys/time.h>
#include <sys/resource.h>
#include "mjm_globals.h"

template < int MAX_SIZE=10 > class  mjm_timing 
{
private:
typedef mjm_timing<MAX_SIZE> Myt;
typedef std::runtime_error ErTy;
typedef char ChTy;
typedef unsigned int FlagTy;
typedef unsigned int IdxTy;
public:

enum {MAXDTA= MAX_SIZE }; 
static  double getcputime(void)        
             { struct timeval tim;        
               struct rusage ru;        
               getrusage(RUSAGE_SELF, &ru);        
               tim=ru.ru_utime;        
               double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
               tim=ru.ru_stime;        
               t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
               return t; }    

//http://rabbit.eng.miami.edu/info/functions/time.html#getrusage
// perf meteres.
static void mems(OsTy & os)
{
rusage ru;
getrusage(RUSAGE_SELF, &ru);
// also can get diusk IO etc in other fields 
os<<MM_MARK<<" cpu="<<getcputime()<<" res="<<ru.ru_maxrss<<" code="<<ru.ru_ixrss<<" globs="
<<ru.ru_idrss <<" stack="<<ru.ru_isrss<<CRLF;
os<<MM_MARK<<" faults="<<ru.ru_minflt<<" "<<ru.ru_majflt<<" "<<ru.ru_nswap<<CRLF;



}
// alt only has microsecond resolution 


//rt
typedef timeval ClockTy;
//typedef int64 ClockerTy;
typedef long  ClockerTy;

ClockTy m_clock_origin;
ClockTy get_clock()
{
ClockTy t;
::gettimeofday(&t, NULL);
return t;

}
// this returns microseconds, 
// this seems to have returned a negative value at
// recent run??? jjjj
ClockerTy  delta_time(const ClockTy & t, const ClockTy & o)
{
typedef ClockerTy Tx;
const Tx s1=1000;
// torn here, want fp to make mult but div will be
// exact? Fixed point? 
// const Tx s2=1000;
Tx sec= t.tv_sec-o.tv_sec;
Tx musec=t.tv_usec-o.tv_usec;
// doh, seconds mult by 1k musec divice
// this seems to create millisecond resolution, bad 
//return (ClockerTy)(sec*s1+musec/s2);
return (ClockerTy)(sec*s1*s1+musec);

}
ClockerTy elapsed_rt()
{ return delta_time(get_clock(),m_clock_origin); } 

#define USE_SYS_TIME
#ifdef USE_SYS_TIME

typedef double TimeTy;
typedef timespec TimerTy;


public:

TimerTy now() 
{
TimerTy tm;
// I want a wall clock but limited to 1 microsecond re apparently 
// clock_gettime(CLOCK_MONOTONIC,&tm);
clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tm);
return tm; 
}
TimeTy elapsed() {
return timex(now(),time_origin); 

} 
private:
public:
// no operator for time spec?
// consider private etc to get inline
TimeTy timex(const TimerTy & t1, const TimerTy &t2)
{
__time_t ts=t1.tv_sec-t2.tv_sec;
long int tn=t1.tv_nsec-t2.tv_nsec;
TimeTy tt=((TimeTy)ts)+((TimeTy)tn)*1e-9;
// osx()<<MM_MARK<<" "<< (tt)<<" "<<" "<<CRLF; 
return tt;
//return (TimeTy)((t1-t2))/CLOCKS_PER_SEC;
}
void init_times(const ChTy ** lbl )
{

for ( IdxTy i=0; i<MAXDTA; ++i) {
act[i]=0; dta[i]=time_origin; dtai[i]=time_origin; dtaf[i]=time_origin;
dta_lbl[i]=NULL; } 
IdxTy i=0;
while ( lbl[i]!=0)
{

dta_lbl[i]=lbl[i];
++i;
if ( i>=MAXDTA) return; 
}

dta_lbl[i]=NULL;

}
void start_time(IdxTy i) { dtai[i]=now(); }
void mark_time(IdxTy i) 
{ dtaf[i]=now();
//act[i]+=((double)(dtaf[i]-dtai[i]))/(double)CLOCKS_PER_SEC;
act[i]+=timex(dtaf[i],dtai[i]);

}
void dump_times(OsTy & os)
{
for ( IdxTy i=0; i<MAXDTA; ++i)
{ 
if ( dta_lbl[i]!=NULL) { os<<dta_lbl[i]<<"="<<act[i]<<" "; } 
}
os<<CRLF;
}

#else

typedef double TimeTy;
typedef clock_t TimerTy;

public:
TimerTy now() { return ::clock(); }
TimeTy elapsed() {
return timex(::clock(),time_origin); 

} 
private:
// consider private etc to get inline
TimeTy timex(const clock_t & t1, const clock_t &t2)
{
osx()<<MM_MARK<<" "<< (t1-t2)<<" "<<CLOCKS_PER_SEC<<CRLF; 
return ((TimeTy)(t1-t2))/(TimeTy)CLOCKS_PER_SEC;
//return (TimeTy)((t1-t2))/CLOCKS_PER_SEC;
}
void init_times()
{

for ( IdxTy i=0; i<MAXDTA; ++i) {
act[i]=0; dta[i]=time_origin; dtai[i]=time_origin; dtaf[i]=time_origin;
dta_lbl[i]=NULL; } 
//dta_lbl[0]="READ";
//dta_lbl[1]="PROC";
//dta_lbl[2]="ACCU";
//dta_lbl[3]="ERR";

}
void start_time(IdxTy i) { dtai[i]=clock(); }
void mark_time(IdxTy i) 
{ dtaf[i]=clock();
//act[i]+=((double)(dtaf[i]-dtai[i]))/(double)CLOCKS_PER_SEC;
act[i]+=timex(dtaf[i],dtai[i]);

}
void dump_times(OsTy & os)
{
for ( IdxTy i=0; i<MAXDTA; ++i)
{ 
if ( dta_lbl[i]!=NULL) { os<<dta_lbl[i]<<"="<<act[i]<<" "; } 
}
os<<CRLF;
}
 
#endif

// enum {MAXDTA=10, DTAREAD=0, DTAPROC=1, DTAACC=2,DTAERR=3};
TimerTy dta[MAXDTA],dtai[MAXDTA], dtaf[MAXDTA],time_origin;
TimeTy act[MAXDTA];
//const ChTy ** dta_lbl;
const ChTy * dta_lbl[MAXDTA];



mjm_timing (): m_clock_origin(get_clock()),  
time_origin(now())// , dta_lbl(0)
{
//setOs(&std::cout,0);
/*
setOs(&std::cerr,1);
setOs(&std::cout,2);
setOs(&std::cerr,3);
init_times();
*/

}




~mjm_timing()
{
}

}; // mjm_timing


#endif

