#ifndef PROB_CURRENT_H__
#define PROB_CURRENT_H__
/***

Just reduce this to brute force P current, 
negf_xue_2001

First-Principles Based Matrix-Green’s Function Approach to Molecular Electronic
Devices: General Formalism
Yongqiang Xue
∗
Department of Chemistry and Materials Research Center, Northwestern University, Evanston, IL 60208 and School of
Electrical and Computer Engineering, Purdue University, West Lafayette, IN 47907
Supriyo Datta
School of Electrical and Computer Engineering, Purdue University, West Lafayette, IN 47907
Mark A. Ratner
Department of Chemistry and Materials Research Center, Northwestern University, Evanston, IL 60208
(February 1, 2008)
Transport in molecular electronic devices is different from
i

uthor to whom correspondence should be addressed. Present address: Department of Chemistry and Materials Research
Center, Northwestern University, Evanston, IL 60208. Electronic mail: ayxue@chem.nwu.edu




i

since there is no “lesser” self-energy operator associated with Vxc . And we can express the correlation function in
terms of the distribution in each electrodes:
G< (E) = i[GR (E)ΓL (E)GA (E)]f (E − μL ) + i[GR (E)ΓR (E)GA (E)]f (E − μR )
(3.17)
where the products within the brackets are matrix products. Every physical observable of interest can be computed
from the matrix correlation function G< . In particular, the current density is:
ij
J(r) =
′
dEG< (E) lim (∇ − ∇)φi (r)φ∗ (r′ )
j
ij
′
dEJ(r; E) = 1/2
r →r
ij
(3.18)
The terminal current can be calculated numerically by integrating the current density J(r) over the boundary
surface between the molecule and the electrodes or any cross sectional area in the “extended molecule” due to the
current continuity in the “extended molecule” region. 

marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp$ googlelink -history | grep 0112136
Wed Oct 31 10:30:31 EDT 2012 /home/marchywka/MyDocs/reprints http://arxiv.org/pdf/cond-mat/0112136 PARAMS rc=0 Saving to: `0112136' -> 0112136.pdf
marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp$ svn copy jdft_geo.h negf_xue_2001.h
A         negf_xue_2001.h


jdft_geo.h



jdft_tool : from mjm_perf for manipulating the binary output and
other stuff 
mjm_perf.cpp: launch a bunch of threads and dump timing stats.
derived from, used copies of headers etc:
logs_db_lut.cpp : integrate the various log and othere tasks inot one place

duphus_agg.cpp: derioved from 
x-ui-util.cpp : Phluant x-windows manipulator for qa and other
image needs. 

Mike Marchywka May 2011
May 30, 2011 : ran fine overnight but mem diag is zed,
start adding "Break" after load LOL, input escapes, temp files
on too big etc


*/

#include "mjm_globals.h"
#include "mjm_strings.h"
#include "mjm_io.h"
#include "mjm_timing.h"
#include "mjm_config.h"

#include <pthread.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <memory.h>

#include "mjm_interchanges.h"
#include "spatial_util.h"
//typedef unsigned int IdxTy;
class  prob_current 
{
private:
typedef prob_current  Myt;

typedef mjm_generic_traits Tr;
typedef computational_traits Tc;

typedef Tc::compute_type  DoTy;

typedef Tr::ChTy ChTy; 
typedef Tr::ErTy ErTy;
// in traits? 
typedef unsigned int FlagTy;
typedef unsigned int IdxTy;
typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;

typedef  sample_blob BlobTy;
typedef  simple_complex WeightTy;
typedef  simple_complex CoTy;
typedef std::vector<BlobTy> Blobs;
typedef std::vector<WeightTy> Weights;

void usage(); 




// have prinln array 
// OsTy * m_ss[4];
IoTy m_io;
 IdxTy m_id; 

// raw wavefunctions
Blobs m_inputs;
// scaled or derived etc.
Blobs m_blobs,m_star,m_grad_dot,m_grad_star;
Weights m_weights;
std::vector<IdxTy> m_dims;

public:
class SetSelector
{
typedef unsigned int IdxTy;
IdxTy m_start, m_n, m_actual_n;
bool m_default;
public:
SetSelector(): m_start(0),m_n(0),m_actual_n(0), m_default(true)
{}
void select_range(const IdxTy s, const IdxTy n)
{
m_start=s;
m_n=n;
m_default=false;
}
bool take(const IdxTy i ) const 
{
if ( m_default) return true;
if ( (i>=m_start)&&(i<=(m_start+m_n))) return true;
return false;
}


};
typedef SetSelector Set_Selector;
prob_current() //:m_data(0),m_sz(0)
{
// memset anyone LOL? 

setOs(&std::cout,0);
setOs(&std::cerr,1);
setOs(&std::cout,2);
setOs(&std::cerr,3);

}


~prob_current()
{
//std::cout<<MM_MARK<<" dtor for duphus called at " // <<elapsed() <<CRLF; 
//delete [] m_data;
}
void add_energies(IsTy * is, const IdxTy flags, const Set_Selector & ss )
{
double e;
IdxTy i=0;
IdxTy in=0;
while ( !is->eof()&&is->good())
{
(*is)>>e;
bool add_this_one=ss.take(in);
++in;
if ( !(!is->eof()&&is->good())) break; 
if ( add_this_one)
{
std::cout<<MM_MARK<<" adding energfy "<<i<<" energ "<<e<<" "<<i<<CRLF;
m_inputs[i].energy(e);

++i;
} // add
}


}

// add things from stream to the grid
void add ( IsTy * is)
{
add(is, Set_Selector());

}

bool add_raw_wfn ( const StrTy & nm, const IdxTy in, const Set_Selector & s,
const IdxTy n1, const IdxTy n2, const IdxTy n3, const bool is_real=false)
{
if ( m_dims.size()==0)
{m_dims.push_back(n1);
m_dims.push_back(n2);
m_dims.push_back(n3);
}

bool take_this_one=s.take(in);
if ( take_this_one)
{

std::cout<<MM_MARK<<" adding "<<nm<<" idx "<<in<<CRLF;
BlobTy s;
m_inputs.push_back(s);
BlobTy & x=m_inputs.back();
if ( !is_real) x.load(nm,0); 
else x.load(nm,1);
x.add_dim(n1);
 x.add_dim(n2);
 x.add_dim(n3);

}
return take_this_one;

}


void add ( IsTy * is, const Set_Selector & s)
{
StrTy nm;
double r,i;
IdxTy in =0;
// dimensions of grid followed by names and wieghts
IdxTy n1,n2,n3;
(*is)>>n1>>n2>>n3;
m_dims.push_back(n1);
m_dims.push_back(n2);
m_dims.push_back(n3);
bool fixing=false;
bool is_real=false;
IdxTy f1,f2,f3;
while ( !is->eof()&&is->good())
{
(*is)>>nm;

if ( !(!is->eof()&&is->good())) break; 
if ( nm=="#" ) continue;
if ( nm=="!REMAP" ) { fixing=true; (*is)>>f1>>f2>>f3; continue; }
if ( nm=="!REAL" ) { is_real=true; continue; }


(*is)>>r>>i;
bool take_this_one=s.take(in);
if ( take_this_one)
{

std::cout<<MM_MARK<<" adding "<<nm<<" weights "<<r<<" "<<i<<CRLF;
BlobTy s;
m_inputs.push_back(s);
BlobTy & x=m_inputs.back();
if ( !is_real) x.load(nm,0); 
else x.load(nm,1);
x.add_dim(n1);
 x.add_dim(n2);
 x.add_dim(n3);
if ( fixing)
{
x.remap(f1,f2,f3);


}
/*
m_inputs.back().load(nm);
m_inputs.back().add_dim(n1);
m_inputs.back().add_dim(n2);
m_inputs.back().add_dim(n3);
*/

m_weights.push_back(WeightTy(r,i));
} // take this one.
++in;
}


std::cout<<MM_MARK<<" done with adding  "<<m_inputs.size()<<CRLF;
}
void load_selected_weights(IsTy * is, const Set_Selector & s)
{
IdxTy in =0;

while ( !is->eof()&&is->good())
{
DoTy i=0;
DoTy r;
(*is)>>r;
if ( s.take(in)) 
{
std::cout<<MM_MARK<<" adding weigth  "<<in<<" "<<r<<CRLF;
m_weights.push_back(WeightTy(r,i));
}

++in;
} 

}



void reset(){}
void setup()
{


}
// selection made during loading.. 
void setup_noithing(const Set_Selector & s )
{
const IdxTy sz=m_inputs.size();
const IdxTy szw=m_weights.size();
std::cout<<MM_MARK<<" setup_noithing  "<< szw<< " wieghts and inp "<<sz<<CRLF;
reset();
for ( IdxTy i=0; i<sz; ++i)
{
// normalize to 8 iirc not 1
BlobTy & xx=m_inputs[i]; // *CoTy(0,0);
m_blobs.push_back(xx);
m_blobs.back().weight(m_weights[i]);
//m_inputs.pop_front();
m_inputs[i]=BlobTy();
//m_star.push_back(m_blobs[i].star());
//m_star.back().weight(m_weights[i]);
}

}

 void setup_simple( )
{
const IdxTy sz=m_inputs.size();
std::cout<<MM_MARK<<" setup_simple   no derivates"<<sz<<CRLF;
reset();
for ( IdxTy i=0; i<sz; ++i)
{
BlobTy xx=m_inputs[i].normalize(); // *CoTy(0,0);
if (m_blobs.size()!=0) m_blobs.push_back(m_blobs.back()+=xx*m_weights[i]);
 else m_blobs.push_back(xx*m_weights[i]);
m_star.push_back(m_blobs[i].star());
//m_grad_dot.push_back(m_blobs[i].grad_dot(dir));
//m_grad_star.push_back(m_star[i].grad_dot(dir));
} // for
}
template <class Ty> void setup(const Ty & dir )
{

//Blobs m_inputs;
// scaled or derived etc.
//Blobs m_blobs,m_star,m_grad_dot_d,m_grad_star;
//Weights m_weights;
const IdxTy sz=m_inputs.size();
std::cout<<MM_MARK<<" setup  "<<sz<<CRLF;
// should clear the old crap first really
reset();
for ( IdxTy i=0; i<sz; ++i)
{
//std::cout<<MM_MARK<<" add entry "<<i<<CRLF;

//BlobTy xx=m_inputs[i].normalize().fft_quads(); // *CoTy(0,0);
// revert fftquad for 2013-01-03
BlobTy xx=m_inputs[i].normalize(); // *CoTy(0,0);

//BlobTy xx=m_inputs[i]*CoTy(0,0);
//m_blobs.push_back(m_inputs[i].fft_quads()*m_weights[i]);
if (m_blobs.size()!=0) m_blobs.push_back(m_blobs.back()+=xx*m_weights[i]); else 
m_blobs.push_back(xx*m_weights[i]);
//m_blobs.push_back(m_inputs[i]*m_weights[i]);
//std::cout<<MM_MARK<<" add entry  start "<<i<<CRLF;
m_star.push_back(m_blobs[i].star());
//std::cout<<MM_MARK<<" add entry grad "<<i<<CRLF;
m_grad_dot.push_back(m_blobs[i].grad_dot(dir));
//std::cout<<MM_MARK<<" add entry grad dot"<<i<<CRLF;
m_grad_star.push_back(m_star[i].grad_dot(dir));

//m_blobs[i]=m_blobs[i].grad_dot<Ty,BlobTy::smear>(dir);
//m_star[i]=m_star[i].grad_dot<Ty,BlobTy::smear>(dir);

}





} 


typedef CoTy  (*FuTy )(const CoTy &, const CoTy &, const CoTy &, const CoTy &); 


static CoTy  current(const CoTy & a , const CoTy & b , const CoTy & c , const CoTy & d) { return a*b-c*d;}
static CoTy  mag(const CoTy & a , const CoTy & b , const CoTy & c , const CoTy & d) { return b*c;}
//static CoTy  mo(const CoTy & a , const CoTy & b , const CoTy & c , const CoTy & d) { return d*c;}
static CoTy  mo(const CoTy & a , const CoTy & b , const CoTy & c , const CoTy & d) { return d;}

class Pi
{
public:
IdxTy dim;
IdxTy vali;
IdxTy thickness;


}; 



void integrate(const IdxTy dim, const IdxTy vali, const IdxTy t=1 ) { ix<current>(dim,vali,t); }
void n(const IdxTy dim, const IdxTy vali, const IdxTy t=1 ) { ix<mag>(dim,vali, t); }
void p(const IdxTy dim, const IdxTy vali, const IdxTy t=1) { ix<mo>(dim,vali,t); }
void integrate(const Pi & p ) { ix<current>(p); }
void n(const Pi & p ) { ix<mag>(p); }
void p(const Pi & p ) { ix<mo>(p); }

// integrate in a plane in dim =x,y, or z at value val.
template <FuTy fu= current> void ix(const Pi & p)  {}
template <FuTy fu= current> void ix(const IdxTy dim, const IdxTy vali, const IdxTy thickness=1)
{
CoTy tot=CoTy(0,0);
const IdxTy sz=m_blobs.size();
const IdxTy first=sz-1;
std::cout<<MM_MARK<<" input size in integrate is "<<sz<<" thick="<<thickness<<CRLF;
// for each wavefunction weighted by occupancy, 
typedef double D; 
D * acc=NULL; 
CoTy * dacc=NULL;
IdxTy ax=0, ay=0;
const IdxTy valf=vali+thickness;
for (IdxTy val=vali; val<valf; ++val)
{

for ( IdxTy i=first; i<sz; ++i)
{
// should return a complex
typedef double * X;
X blob=m_blobs[i].d();
X star=m_star[i].d();
X g=m_grad_dot[i].d();
X gstar=m_grad_star[i].d();
// these are comlex, doh. 
// these need iterators, 
threeDIterator<double> * itor= NULL;
// this needs to be yanked and the itor reset... 
switch ( dim)
{
case 0: {ax=m_dims[1]; ay=m_dims[2];  itor= new X3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 
case 1: {ax=m_dims[0]; ay=m_dims[2]; itor= new Y3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 
case 2: {ax=m_dims[0]; ay=m_dims[1]; itor= new Z3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 

}
if ( acc==NULL) {
IdxTy asz=ax*ay;
// this seems to call default ctor LOL
 dacc = new CoTy[asz]; acc = new D[asz]; 
::memset(acc,0,sizeof(double)*asz);
for (IdxTy  jj=0; jj<asz; ++jj) dacc[jj]=CoTy(0,0); 
 }

CoTy totwf=CoTy(0,0);
//for ( IdxTy j=j0; j<jend; j+=jd)
int cnt=0;
while (!(*itor).done())
{
IdxTy k=**itor; // j<<1; /// arrgggg ... 
CoTy a=CoTy(gstar+k);
CoTy  b=CoTy(blob+k);
CoTy  c=CoTy(star+k);
CoTy d=CoTy(g+k);
//std::cout<<MM_MARK<<" "<<" "<<a.r()<<" "<<b.r()<<" "<<c.r()<<" "<<d.r()<<CRLF;
//std::cout<<MM_MARK<<" "<<" "<<a.i()<<" "<<b.i()<<" "<<c.i()<<" "<<d.i()<<CRLF;
//CoTy xtot=a*b-c*d;
CoTy xtot=fu(a,b,c,d); /// a*b-c*d;
//xtot=b*c;
tot+=xtot;
totwf+=xtot;
std::cout<<MM_MARK<<" MApdetail  "<<(*itor).x()<<" "<<(*itor).y()<<" "<<(*itor).z()<<" "<<xtot.r()<<" "<<xtot.i()<<" "<<sqrt(xtot.abs())<<CRLF;
std::cout<<MM_MARK<<" ITEMS  "<<i<<" "<<(*itor).x()<<" "<<(*itor).y()<<" "<<(*itor).z()<<" "<<a.abs()<<" "<<b.abs()<<" "<<(c.abs())<<" "<<d.abs()<<CRLF;
// if ( i==5) if ( (xtot.r()!=0) ||(xtot.i()!=0)) std::cout<<MM_MARK<<" JUNK  "<<(*itor).x()<<" "<<(*itor).y()<<" "<<(*itor).z()<<" "<<xtot.r()<<" "<<xtot.i()<<" "<<sqrt(xtot.abs())<<CRLF;
// complex
//double dx=gstar[i]*blob[i]-star[i]*g[i];
acc[cnt]+=xtot.i();
dacc[cnt]+=xtot;
++cnt;
if ( !(itor->is_ok()))
std::cout<<MM_MARK<<" itor is bad "<<(*itor).val()<<" "<<(*itor).ptr()<<CRLF;
++*itor;
} // j
delete itor;
std::cout<<MM_MARK<<" "<<cnt<<" a="<<totwf.abs()<<" r="<<tot.r()<<" i="<<tot.i()<<" tot="<<tot.abs()<<CRLF;

}
}
IdxTy idx=0;
for (IdxTy i=0; i<ay; ++i)
{for (IdxTy j=0; j<ax; ++j)
{ ++idx;
std::cout<<MM_MARK<<" MAP  "<<j<<" "<<i<<" "<<vali<<" "<<0<<" "<<acc[idx]<<" "<<acc[idx]<<CRLF;
std::cout<<MM_MARK<<" MAPC  "<<j<<" "<<i<<" "<<vali<<" "<<0<<" "<<dacc[idx].r()<<" "
<<dacc[idx].i()<<" "<<::sqrt(dacc[idx].abs())<<CRLF;
}}


delete [] acc;
delete [] dacc;

}

// iterate and sum or do something with each result
class isum 
{
typedef isum Myt;
typedef double D;
typedef simple_complex  CoTy;

threeDIterator<double> * m_itor; // = NULL;
CoTy * m_acc; D * m_dacc;
D  m_px,m_py,m_pz;
IdxTy m_ax, m_ay,m_az,m_cnt; // two d space of ieration
IdxTy m_x0, m_y0,m_z0;
public:
virtual ~isum()
{
delete [] m_acc; delete [] m_dacc;
}
template <class Ty> void init(const IdxTy val, const IdxTy dim, Ty & _dims)
{
switch ( dim)
{
case 0: {m_ax=_dims[1]; m_ay=_dims[2];  m_az=_dims[0]; m_itor= new X3Iterator(val,_dims[0],_dims[1],_dims[2],2,NULL); break  ;} 
case 1: {m_ax=_dims[0]; m_ay=_dims[2]; m_az=_dims[1]; m_itor= new Y3Iterator(val,_dims[0],_dims[1],_dims[2],2,NULL); break  ;} 
case 2: {m_ax=_dims[0]; m_ay=_dims[1]; m_az=_dims[2]; m_itor= new Z3Iterator(val,_dims[0],_dims[1],_dims[2],2,NULL); break  ;} 

}
IdxTy asz=m_ax*m_ay;
m_x0=m_ax>>1;
m_y0=m_ay>>1;
m_z0=m_az>>1;

// this seems to call default ctor LOL
m_acc = new CoTy[asz]; 
m_dacc = new D[asz]; 
::memset(m_dacc,0,sizeof(double)*asz);
for (IdxTy  jj=0; jj<asz; ++jj) m_acc[jj]=CoTy(0,0); 


}
template <class Ty> isum(const IdxTy val, const IdxTy dim, Ty & _dims, const Myt& that )
: m_itor(NULL), m_px(that.m_px), m_py(that.m_py), m_pz(that.m_pz)
{ init(val,dim,_dims); }
template <class Ty> isum(const IdxTy val, const IdxTy dim, Ty & _dims)
: m_itor(NULL), m_px(0), m_py(0), m_pz(0)
{ init(val,dim,_dims); }





const IdxTy & operator*() { return **m_itor; }
bool done() { return m_itor->done(); }
// wtf
virtual Myt & operator++()
{
++(*m_itor);  ++m_cnt;
return *this;
}
void add(const CoTy & b)
{
// here we want to update px and py as abs value squared time index. 
int  kx=(*m_itor).x()-m_x0;
int  ky=(*m_itor).y()-m_y0;
int kz=(*m_itor).z()-m_z0;
// b is the fft coef of this momemtum eigenvector. 
// the mag squared is amt, 
const D m=b.abs();
//std::cout<<MM_MARK<<" PCOMP "<<(*m_itor).x()<<" "<<m_x0<<" "<<kx<<" "<<m<<CRLF;
std::cout<<MM_MARK<<" PCOMP "<<kx<<" "<<ky<<" "<<kz<<" "<<m<<CRLF;
//std::cout<<MM_MARK<<" PCOMP "<<m_x0<<" "<<m_y0<<" "<<m_z0<<" "<<m<<CRLF;
m_px+=m*kx;
m_py+=m*ky;
m_pz+=m*kz;

//m_acc[m_cnt]+=b;
//m_dacc[m_cnt]+=b.i();

}
void dump()
{
OsTy & x = std::cout;
x<<MM_MARK<<
" px=";du(x,m_px)<<" py= ";du(x,m_py)<<" pz=";;du(x,m_pz)<<CRLF;

}
template <class Ty> Ty &  du(Ty & os, const D  & xx)
{
os<<" "<<xx;
return os;
}

template <class Ty> Ty &  du(Ty & os, const CoTy & x)
{
os<<" "<<x.r()<<" "<<x.i();
return os;
}

} ; // isum,

void dipole_etc(OsTy & osd=std::cout)
{
const IdxTy sz=m_blobs.size();
const IdxTy first=0;
const IdxTy xs=m_dims[0];
const IdxTy ys=m_dims[1];
const IdxTy zs=m_dims[2];
const IdxTy gsz=xs*ys*zs;
OsTy & os = std::cout;
os<<MM_MARK<<" find dipoles  "<<sz<<" eleemnts="<<gsz<<CRLF;
typedef double D;
for ( IdxTy i=first; i<sz; ++i)
{
//typedef double * X;
typedef CoTy  * X;
X  blobi=(X)m_blobs[i].d();
//D pr=0;
//D pi=0;
CoTy tx=CoTy(0,0);
CoTy ty=CoTy(0,0);
CoTy tz=CoTy(0,0);
for ( IdxTy j=i+1; j<sz; ++j)
{
IdxTy loc=0;

//X  blobj=(X)m_star[j].d();
X  blobj=(X)m_blobs[j].d();

for ( IdxTy x=0; x<xs; ++x)
{
for ( IdxTy y=0; y<ys; ++y)
{
for ( IdxTy z=0; z<zs; ++z)
{
const CoTy p=blobi[loc].conj(blobj[loc]);
// you think compiler can remove the zero mult LOL?
tx+=p*CoTy(x,0);
ty+=p*CoTy(y,0);
tz+=p*CoTy(z,0);
++loc;
} // z

}  // y


}  // x
D de=m_blobs[j].energy()-m_blobs[i].energy();
// these are sqraured
D txa = tx.abs();
D tya = ty.abs();
D tza = tz.abs();
D tt=::sqrt(txa+tya+tza);
osd<<MM_MARK<<" "<<i<<" "<<j<<" "<<de*27.2<<" "<<m_blobs[i].weight().mag()<<" "<<m_blobs[j].weight().mag()<<" "<<tt <<" "<<tx.r()<<" "<<tx.i()<<" "<<tx.mag();
os<<" "<<ty.r()<<" "<<ty.i()<<" "<<ty.mag();
os<<" "<<tz.r()<<" "<<tz.i()<<" "<<tz.mag()<<CRLF;


}


}
}

void ftrs(const IdxTy dim, const IdxTy vali, const IdxTy t=1 ) 
{ 
const IdxTy sz=m_blobs.size();
const IdxTy first=sz-1;
const IdxTy thickness=t;
std::cout<<MM_MARK<<" input size in integrate is "<<sz<<" thick="<<thickness<<CRLF;
// for each wavefunction weighted by occupancy, 
//typedef double D;
const IdxTy _vali=0;
//const IdxTy valf=vali+thickness;
// this is not right 
const IdxTy _valf=m_dims[2];
isum * old_itor=NULL;
//for (IdxTy val=vali; val<valf; ++val)
for (IdxTy val=_vali; val<_valf; ++val)
{
// this must only use 1
for ( IdxTy i=first; i<sz; ++i)
{
//if ( val==_vali) m_blobs[i].fft(); //.fft_quads();
if ( val==_vali) m_blobs[i]=m_blobs[i].fft().fft_quads();
// should return a complex
typedef double * X;
X blob=m_blobs[i].d();

isum * itor = NULL;
if ( old_itor!=NULL) itor= new  isum( val, dim, m_dims,*old_itor); else 
itor= new  isum( val, dim, m_dims);
 delete  old_itor; 
old_itor=itor;  


//int cnt=0;
while (!(*itor).done())
{
IdxTy k=**itor; // j<<1; /// arrgggg ... 
CoTy  b=CoTy(blob+k);
// 
(*itor).add(b);

++*itor;
}
(*itor).dump();
//delete itor;

}
}

delete old_itor;
}


///////////////////////////////////////////////////////////////////////

void planes()
{
//CoTy tot=CoTy(0,0);
const IdxTy sz=m_blobs.size();
const IdxTy first=sz-1;
const IdxTy thickness=1;
std::cout<<MM_MARK<<" input size in planes() is "<<sz<<" thick="<<thickness<<CRLF;
// for each wavefunction weighted by occupancy, 
typedef double D; 
for ( IdxTy dim=0; dim<3; ++dim)
{ 
IdxTy min_location=0;
D min_value=1e99;
const IdxTy planes=m_dims[dim];
//IdxTy ax=0, ay=0;
for (IdxTy val=0; val<planes; ++val)
{

for ( IdxTy i=first; i<sz; ++i)
{
// should return a complex
typedef double * X;
X blob=m_blobs[i].d();
X star=m_star[i].d();
//X g=m_grad_dot[i].d();
//X gstar=m_grad_star[i].d();
// these are comlex, doh. 
// these need iterators, 
threeDIterator<double> * itor= NULL;
// this needs to be yanked and the itor reset... 
switch ( dim)
{
case 0: {//ax=m_dims[1]; ay=m_dims[2];  
itor= new X3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 
case 1: {//ax=m_dims[0]; ay=m_dims[2]; 
itor= new Y3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 
case 2: {//ax=m_dims[0]; ay=m_dims[1]; 
itor= new Z3Iterator(val,m_dims[0],m_dims[1],m_dims[2],2,NULL); break  ;} 

}

CoTy totwf=CoTy(0,0);
D tot=0;
//for ( IdxTy j=j0; j<jend; j+=jd)
//int cnt=0;
// iterates over 1 plane
while (!(*itor).done())
{
IdxTy k=**itor; // j<<1; /// arrgggg ... 
//CoTy a=CoTy(gstar+k);
CoTy  b=CoTy(blob+k);
CoTy  c=CoTy(star+k);
CoTy xtot=b*c; //fu(a,b,c,d); /// a*b-c*d;
//xtot=b*c;
tot+=xtot.abs();
totwf+=xtot;
++*itor;
}
if ( tot<min_value)  { min_location=val; min_value=tot; } 
std::cout<<MM_MARK<<" dim= "<<dim<<" plane= "<<val<<" total= "<<tot<<CRLF;
delete itor;
//std::cout<<MM_MARK<<" "<<cnt<<" a="<<totwf.abs()<<" r="<<tot.r()<<" i="<<tot.i()<<" tot="<<tot.abs()<<CRLF;

} //i ( thickness)  
} // plane 
std::cout<<MM_MARK<<" min for "<<dim<<" plane= "<<min_location<<" value= "<<min_value<<CRLF;
} // dim 


} // method
/*
// iterate and sum or do something with each result
class isum 
{
typedef isum Myt;
typedef double D;
typedef simple_complex  CoTy;

threeDIterator<double> * m_itor; // = NULL;
CoTy * m_acc; D * m_dacc;

*/





//////////////////////////////////////////////////////////////////////////


OsTy & os(const IdxTy i ) { return m_io.os(i); } 
void setOs(OsTy * os,const IdxTy i ) { m_io.setOs(os,i); } 

//void config(IsTy & is) { m_config.load(is); }

}; // logs_db_lut


//////////////////////////////////////////////////////////////
// impl

#endif


