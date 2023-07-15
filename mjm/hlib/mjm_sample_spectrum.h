#ifndef MJM_SAMPLE_SPECTRUM_H__
#define MJM_SAMPLE_SPECTRUM_H__


#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
//#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>
#include <memory.h>
//#include <cblas.h>
#include "mjm_conflicts.h" 
// not sure this works with linpack?
#include <complex>



/*
extern "C"
{       void zheevr_(char* JOBZ, char* RANGE, char* UPLO, int * N, CoTy* A, int * LDA,
                double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
                double* W, CoTy* Z, int* LDZ, int* ISUPPZ, CoTy* WORK, int* LWORK,
                double* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO);

void cgesv_(int * N, int * NRHS, CoTy * A, int* LDA, int * IPIV, CoTy * B,
int * LDB, int* INFO);


// FORTRAN adds _ after all the function names 
// and all variables are called by reference
double ddot_( const int *N, const double *a, const int *inca, const double *b, const int *incb );

double ddot( int N, double *a, int inca, double *b, int incb ){
  return ddot_( &N, a, &inca, b, &incb );
};
};
*/

/*
int main( int argc, char** argv ){
  // you can define the arrays in one of two ways
  // on the heap
  double *a = (double*) malloc( 3 * sizeof(double) );
  a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
  // on the stack
  double b[3] = { 4.0, 5.0, 6.0 };

  double dot_product = ddot( 3, a, 1, b, 1 );
  printf(" The dot product is: %f \n",dot_product );

  return 0;
};

*/



class mjm_sample_spectrum
{
typedef mjm_sample_spectrum Myt;

typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;

IdxTy m_sig_length;
D m_carrier_omega, m_mod_omega, m_dt,m_a;
D * m_signal;
void Assign(const Myt & that)
{
m_sig_length=that.m_sig_length;
m_carrier_omega= that.m_carrier_omega;
m_mod_omega=that.m_mod_omega;
m_dt=that.m_dt;
m_a= that.m_a;
delete [] m_signal;
m_signal= new D[m_sig_length];
::memcpy(m_signal, that.m_signal,m_sig_length*sizeof(D));

}
public:
template <class Os> void dump_signal(Os & os, const StrTy & nm, const IdxTy flags)
{
const bool name=((flags&1)==0);
const bool time=((flags&2)==0);
const bool idx=((flags&4)==0);
for ( IdxTy i=0; i<m_sig_length; ++i)
{
if ( name) os<<nm<<" ";
if (time) os<<i*m_dt<<" ";
if (idx) os<<i<<" ";
os<<m_signal[i];

os<<CRLF;
}

}
template <class Os> void dump(Os & os)
{
const D et=m_sig_length*m_dt;
const D pi2=::atan(1)*8;
const D carrier_cycles=m_carrier_omega*et/pi2;
const D mod_cycles=m_mod_omega*et/pi2;

os<<MM_MARK<<" m_sig_length="<<m_sig_length;
os<<" m_carrier_omega="<<m_carrier_omega;
os<<" m_mod_omega="<<m_mod_omega;
os<<" m_dt="<<m_dt;
os<<" m_a="<<m_a;
os<<CRLF;
os<<MM_MARK<<" carrier cycles = "<<carrier_cycles<<" mo cycles = "<<
mod_cycles<<" t = "<<et<<CRLF;

}
public:
//mjm_sample_spectrum():m_signal(0) {  } 
mjm_sample_spectrum(): m_sig_length(0),m_carrier_omega(0), m_mod_omega(0),
m_dt(0), m_a(0),
m_signal(new D[m_sig_length]) {dump(std::cout);  } 

mjm_sample_spectrum(const IdxTy len, const D co, const D mo,
const D dt, const D ma):
m_sig_length(len),m_carrier_omega(co), m_mod_omega(mo),
m_dt(dt), m_a(ma),
m_signal(new D[m_sig_length]) {dump(std::cout);  } 

mjm_sample_spectrum ( const Myt & that)
{
Assign(that);

}
// need to make a copy ctor etc...
~mjm_sample_spectrum() { delete [] m_signal; } 
D mxx(const D t)
{ return 1+m_a*cos(m_mod_omega*t);} 
D signal(const D t)
{
const D mx=mxx(t);
return mx*::cos(m_carrier_omega*mx*t);
}
void make_signal()
{
for ( IdxTy i=0; i<m_sig_length; ++i)
{
// rather slow to make double from int... 
const D t=m_dt*i; // could just add but that would be cum roundoff error
m_signal[i]=signal(t);

}

}

void rectify()
{
D tot=0;
for ( IdxTy i=0; i<m_sig_length; ++i)
{
// rather slow to make double from int... 
//m_signal[i]=::fabs(m_signal[i]);
// power law
m_signal[i]=m_signal[i]*m_signal[i];
tot+=m_signal[i];
}
tot/=m_sig_length;

for ( IdxTy i=0; i<m_sig_length; ++i) m_signal[i]=m_signal[i]-tot;

}



// probably faster to batch these, do many in each traversal
void freq(const D test_omega)
{
D in_phase=0;
D q_phase=0;
for ( IdxTy i=0; i<m_sig_length; ++i)
{
const D t=m_dt*i; // could just add but that would be cum roundoff error
in_phase+=::cos(test_omega*t)*m_signal[i];
q_phase+=::sin(test_omega*t)*m_signal[i];

}
const D m=::sqrt(in_phase*in_phase+q_phase*q_phase);
std::cout<<"DATA"<<" omega= "<<test_omega<<" r.i= "<<in_phase<<" "<<q_phase
<<" m= "<<m<<" p= "<<::atan2(in_phase,q_phase)*45/::atan(1)<<CRLF;

}

}; // sample spectrum


#endif

