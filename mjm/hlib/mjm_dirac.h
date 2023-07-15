#ifndef MJM_TCL_DIRAC_H__
#define MJM_TCL_DIRAC_H__


#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
#include <mjm_38spatial.h>
#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>

#include "mjm_conflicts.h" // #include <cblas.h>
// not sure this works with linpack?
#include <complex>

#include "mjm_sample_spectrum.h"
#include "mjm_io.h"
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



class mjm_dirac

{

typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;


// this is the point holder
typedef std::vector<D> V;
// this is for reference staacks of points. 
typedef std::vector<V> StackTy;
typedef std::vector<StrTy> Plist;

//typedef mjm_molecule MolTy;
//typedef std::map<StrTy,MolTy> MolMap;
typedef mjm_sample_spectrum SigTy;
typedef std::map<StrTy,SigTy*> SigMap;

static SigMap & sigs() 
{
static SigMap m_mols;
return m_mols;
}
static StrTy & current() 
{static StrTy  m_current;
return m_current; }

static SigTy & mol() { static SigTy m_mol; return m_mol; } 





public:

static void  get(Plist & p, const int n, Tcl_Obj * const objv[])
{
p.clear();
for (IdxTy i=0; i<(IdxTy)n; ++i)
{
p.push_back(StrTy(Tcl_GetString(objv[i])));

}

}
typedef ClientData CdTy;
typedef Tcl_Interp InTy;
typedef Tcl_Obj * const ArTy;
typedef int I;
enum { ORIGIN=1 }; 
#undef TCL_PRERAMBLE
#undef TCL_PRERAMBLE_SHORT
#define TCL_PRERAMBLE debug(clientData,interp,objc,objv); Plist p; get(p,objc,objv);  V v; xyz(v,p,ORIGIN,std::min(3+ORIGIN,objc));
 // #define TCL_PRERAMBLE_SHORT debug(cd,in,objc,objv); Plist p; get(p,objc,objv); C & s=state(); V v; xyz(v,p,ORIGIN,std::min(3+ORIGIN,objc));
#define TCL_PRERAMBLE_SHORT debug(cd,in,objc,objv); Plist p; get(p,objc,objv);  V v; xyz(v,p,ORIGIN,objc);

static int debug(CdTy clientData, InTy * interp, I objc,  ArTy  objv[])
//static int debug(ClientData clientData, Tcl_Interp * interp, 
//int objc,  Tcl_Obj * const  objv[])
{
osx()<<MM_MARK<<" data is "<<objc<<CRLF;
for ( IdxTy i=0; i<(IdxTy)objc; ++i) osx()<<MM_MARK<<" "<<Tcl_GetString(objv[i])<<CRLF;
 Tcl_SetObjResult(interp, Tcl_NewStringObj("999", -1));
return 0;
}

static int cursor(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
// return coordinate in p(0) from current location. 
// Tcl_SetObjResult(in, Tcl_NewStringObj(mjm::toString(s.pc(IdxTy(v[0]))).c_str(), -1));
return 0; 
}
/*
// return an error for the estimate of E given well pair params
static D dwell(const D & l, const D & C, const D & v, const D & E)
{

//osx()<<MM_MARK<<" v="<<v<<" E="<<E<<CRLF;
const D  k1=::sqrt(E);
const D  k2=::sqrt(v-E);
const D  A =k1*C;
const D  B =k2*(l*.5-C);
//osx()<<MM_MARK<<" A="<<A<<" B="<<B<<CRLF;

D c1=::cos(A)/::sin(A); c1=c1*c1;
const D x1=::exp(B); 
const D x2=1.0/x1;
D c2=(x1+x2)/(x1-x2);  c2=c2*c2;
//osx()<<MM_MARK<<" c1="<<c1<<" c2="<<c2<<CRLF;
 
return E*(1+c1*c2)-v;
}

static int well_bond(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
if ( v.size()<5)
osx()<<MM_MARK<<" "<<dwell(v[0],v[1],v[2],v[3])<<CRLF;
else
{
D emin=0;
D emax=v[3];
const D del=1e-9;
const IdxTy imax=100;
D offs= 1e-3;
IdxTy iter=0;
D guess=(emax+emin)/2;
D res=0;
while ( (emax-emin>del)&&(iter<imax))
{
res=dwell(v[0],v[1],v[2],guess);
//osx()<<MM_MARK<<"i="<<iter<<" guess="<<guess<<" res="<<res<<CRLF;
// this needs to check nan etc. 
if ( isnan(res)) {emin+=offs; emax-=offs;} else
if ( res<0) emin=guess;
else if ( res>0) emax=guess;
else if ( res==0) {emin=guess; emax=guess; } 
guess=(emax+emin)/2;
++iter;
} 
osx()<<MM_MARK<<" i="<<iter<<" guess="<<guess<<" res="<<res<<" "<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<CRLF;


}
////s.normal(v);
return 0; 
}



static int blasdot(CdTy clientData, InTy * interp, I objc,  ArTy  objv[])
{
TCL_PRERAMBLE
//s.basis(v);
  double *a = (double*) malloc( 3 * sizeof(double) );
  a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
  // on the stack
  double b[3] = { 4.0, 5.0, 6.0 };

  double dot_product = ddot( 3, a, 1, b, 1 );
  printf(" The dot product is: %f \n",dot_product );

  return 0;
}

static int dir(CdTy clientData, InTy * interp, I objc,  ArTy  objv[])
{
TCL_PRERAMBLE
//s.basis(v);
return 0; 
}
static int len(CdTy clientData, InTy * interp, I objc,  ArTy  objv[])
{
TCL_PRERAMBLE
//s.length(v[0]);
return 0; 
}
static int phi(CdTy clientData, InTy * interp, I objc,  ArTy  objv[])
{
TCL_PRERAMBLE
//s.angle(v[0]);
return 0; 
}
*/

static int set(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
current()=p[ORIGIN];
//s.digits()=IdxTy(v[0]);
return 0; 
}

static int load(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
//SigTy & m=sigs()[current()];
std::cout<<MM_MARK<<" loading "<<current()<<" from "<<p[0]<<CRLF;
//mjm_molecule_io mio;
//mio.load(m,p[ORIGIN]);
//cmap()[p[ORIGIN]]= s; 
return 0; 
}


static int store(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
SigTy * m=sigs()[current()];
std::cout<<MM_MARK<<" storgin "<<current()<<" to "<<p[ORIGIN]<< " flags "
<<v[1]<<CRLF;
OsTy * osp = mjm_io::get_out(p[ORIGIN].c_str(),std::cout);
m->dump_signal(*osp,current(),IdxTy(v[1]));
osp->flush();
delete osp;
//mjm_molecule_io mio;
//mio.store(m,std::cout);
//cmap()[p[ORIGIN]]= s; 
return 0; 
}

static int fft(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
typedef mjm_sample_spectrum MSD;
const IdxTy ofset=0;
MSD * msd=sigs()[current()];
//if ( flag==1) msd->rectify();
D oc=v[0+ofset];
for (IdxTy i=0; i< IdxTy(v[2+ofset]); ++i)
{
msd->freq(oc);
oc+=v[1+ofset];
}


return 0;

return 0;
}

static int make(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
typedef mjm_sample_spectrum MSD;
const IdxTy ofset=1;
std::cout<<MM_MARK<< 
" m_sig_length(len),m_carrier_omega(co), m_mod_omega(mo), m_dt(dt), m_a(ma)"
<< " del omega n " <<CRLF;

std::cout<<MM_MARK<<" "<<IdxTy(v[0+ofset])<<" "<<v[1+ofset]<<" "<<v[2+ofset]<<" "<<v[3+ofset]
<<" "<<v[4+ofset] <<CRLF;
MSD * msd = new MSD(IdxTy(v[0+ofset]),v[1+ofset],v[2+ofset],v[3+ofset],v[4+ofset] );
msd->make_signal();
current()=p[ORIGIN];
sigs()[current()]=msd;
return 0;
}
static int dirac_spectrum(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
do_signal(v,0);
return 0;
}

static int dirac_rectify(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
do_signal(v,1);
return 0;
}



/*
mjm_sample_spectrum(const IdxTy len, const D co, const D mo,
const D dt, const D ma):
m_sig_length(len),m_carrier_omega(co), m_mod_omega(mo), m_dt(dt), m_a(ma),
m_signal(new D[m_sig_length]) {  } 

~mjm_sample_spectrum() { delete [] m_signal; } 
D mxx(const D t)
{ return 1+m_a*cos(m_mod_omega*t);} 
D signal(const D t)
*/
template <class Ty > static int  do_signal( Ty & v, const IdxTy flag)
{
typedef mjm_sample_spectrum MSD;
std::cout<<MM_MARK<< 
" m_sig_length(len),m_carrier_omega(co), m_mod_omega(mo), m_dt(dt), m_a(ma)"
<< " del omega n " <<CRLF;

std::cout<<MM_MARK<<" "<<IdxTy(v[0])<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]
<<" "<<v[4] <<CRLF;
MSD * msd = new MSD(IdxTy(v[0]),v[1],v[2],v[3],v[4] );
msd->make_signal();
if ( flag==1) msd->rectify();
D oc=v[5];
for (IdxTy i=0; i< IdxTy(v[7]); ++i)
{
msd->freq(oc);
oc+=v[6];
}


return 0;

}


/*I
static int analyses(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
MolTy & m=mols()[current()];
OsTy & os=std::cout;
typedef mjm_38spatial::analyze_param Tp;
Tp results=Tp(v,p);
MolTy::iterator  mi=m.itor();
std::cout<<MM_MARK<<" mol size is "<<m.size()<<CRLF;
//std::cout<<MM_MARK<<" first atom is  "<<std::vector<double>(*mi).size()<<CRLF;

mjm_38spatial::moment_matrix(results, mi);
results.dump(os);
std::cout<<MM_MARK<<" exiting mol size is "<<m.size()<<CRLF;
// list of string and numerical parameters
//mjm_38spatial::misc_stats(results, m,v,p);
//mjm_38spatial::misc_stats(results, m);
return 0; 
}



static int distances(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
MolTy & m=mols()[current()];
OsTy & os=std::cout;
const StrTy sep=" ";
std::cout<<MM_MARK<<" distancing "<<current()<<CRLF; // <<" from "<<p[0]<<CRLF;
D* edm =mjm_38spatial::exhaustive_distance_matrix(m);
D* d=edm;
std::cout<<MM_MARK<<" distancing  dun calculating now prtin "<<current()<<CRLF; // <<" from "<<p[0]<<CRLF;
const IdxTy sz=m.size();
// this should return it... 
for ( IdxTy i=0; i<sz; ++i)
{
const StrTy ni=StrTy(m[i]);
for ( IdxTy j=0; j<sz; ++j)
{
const StrTy nj=StrTy(m[j]);
if ( (*d)!=0)
{
os<<"distance "<<i<<" "<<j<<" "<< ni<<" "<<nj<<" "<<::sqrt(*(d++))<<sep;
os<<CRLF;
} else ++d;
} // j
} // i 
delete [] edm;
//mjm_molecule_io mio; mio.load(m,p[ORIGIN]);
//cmap()[p[ORIGIN]]= s; 
return 0; 
}



static int jump(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
//s=cmap()[p[ORIGIN]]; 
return 0; 
}

static int csolve(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
//s=cmap()[p[ORIGIN]]; 
// Ax=b 
int d=4;
CoTy*  a=new CoTy[16];
CoTy* b= new CoTy[4];
for ( IdxTy i=0; i<16; ++i) a[i]= CoTy(i,i);
for ( IdxTy i=0; i<4; ++i) b[i]= CoTy(i,i);
std::cout<<MM_MARK << " calling "<<" "<<CRLF;
int* ipiv=new int[4];
//int ldb=0;
int info=999;
// variables, equants, A, N ( LDA), IPIV
cgesv_(&d,&d,a,&d,ipiv,b,&d,&info );
std::cout<<MM_MARK << " info="<<info<<CRLF;
for ( IdxTy i=0; i<4; ++i) std::cout<<MM_MARK<<" b[i]="<<b[i]<<" "; 
for ( IdxTy i=0; i<4; ++i) std::cout<<MM_MARK<<" piv="<<ipiv[i]<<" "; 
std::cout<<CRLF;

return 0; 
}



// using last bond dir as normal, add a group of atoms slightly
// out of that plane.
static int cluster(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
// this requires a name and number as well as theta and phi
// oritnetiaion rlative to last "n"
//s=cmap()[p[ORIGIN]];
// atom name, angle from plane , angel from last normal, N
//s.cluster(p[ORIGIN], v[1], v[2], v[3],IdxTy(v[4]));

return 0; 
}
*/

/*
static int clear(CdTy cd, InTy * in, I objc, ArTy objv[])
{
TCL_PRERAMBLE_SHORT
const IdxTy n=p.size();
for ( IdxTy i=0; i<n ; ++i)
{
const StrTy & param=p[i];
//if ( param=="map") cmap().clear(); 


} // for
return 0; 
}

*/



static int reference(ClientData clientData, Tcl_Interp * interp, 
int objc,  Tcl_Obj * const  objv[])
{
debug(clientData,interp,objc,objv);
return 0; 
}





static int set( Tcl_Interp * tcl_interp)
{
//Tcl_CreateObjCommand(tcl_interp, "well_bond", well_bond,0,0);
//Tcl_CreateObjCommand(tcl_interp, "blasdot", blasdot,0,0);
//Tcl_CreateObjCommand(tcl_interp, "csolve", csolve,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_set", set,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_store", store,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_load", load,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_make", make,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_fft", fft,0,0);
//Tcl_CreateObjCommand(tcl_interp, "distancesmol", distances,0,0);
//Tcl_CreateObjCommand(tcl_interp, "analyze", analyses,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_spectrum_new", dirac_spectrum,0,0);
Tcl_CreateObjCommand(tcl_interp, "dirac_rectify_new", dirac_rectify,0,0);
//Tcl_CreateObjCommand(tcl_interp, "mjm_test", mjm_test,0,0);
/*
Tcl_CreateObjCommand(tcl_interp, "atom", add_atom,0,0);
Tcl_CreateObjCommand(tcl_interp, "ref", reference,0,0);
Tcl_CreateObjCommand(tcl_interp, "normal", normal,0,0);
Tcl_CreateObjCommand(tcl_interp, "phi", phi,0,0);
Tcl_CreateObjCommand(tcl_interp, "len", len,0,0);
Tcl_CreateObjCommand(tcl_interp, "dir", dir,0,0);
Tcl_CreateObjCommand(tcl_interp, "digits", digits,0,0);
Tcl_CreateObjCommand(tcl_interp, "mark", mark,0,0);
Tcl_CreateObjCommand(tcl_interp, "jump", jump,0,0);
Tcl_CreateObjCommand(tcl_interp, "cluster", cluster,0,0);
*/

return 0; 
}
static D atof(const StrTy & f) { return ::atof(f.c_str());} 
static void xyz(V & v, Plist & p, const IdxTy start, const IdxTy f )
{
//const IdxTy f=start+3;
for( IdxTy i=start; i<f; ++i)
{
v.push_back(atof(p[i]));
}
}

} ; //mjm_misc
#endif

