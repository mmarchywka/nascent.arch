#ifndef MJM_BURMANN_UTIL_H__
#define MJM_BURMANN_UTIL_H__

#include "mjm_globals.h"


 #define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
 #include "mjm_bigint_mmpr.h"
//typedef int mjm_bigint_mmpr ;


//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
#include "mjm_integrals.h"
#include "mjm_closed_families.h"

#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"

/*
2017-11-1
 2408  g++ -DTEST_BURMANN3__ -Wall -Wno-unused-function  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -gdwarf-3 -x c++ mjm_burmann.h 
2017-10-15
 2408  g++ -DTEST_BURMANN2__ -Wall -Wno-unused-function  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -gdwarf-3 -x c++ mjm_burmann.h 
 2409  ./a.out 10 0 1 0 > xxx
 2410  vi xxx
 2411  ../delbackup.tex 


Based largely on this ,
[1]
http://www.mathematica-journal.com/2014/11/on-burmanns-theorem-and-its-application-to-problems-of-linear-and-nonlinear-heat-transfer-and-diffusion/#Eq:21

 The Mathematica Journal
Volume 16
H. M. Schöpf, P. H. Supancic
On Bürmann’s Theorem and Its Application to Problems of Linear and Nonlinear Heat Transfer and Diffusion	
Expanding a Function in Powers of Its Derivative


*/

class mjm_burmann
{

typedef mjm_burmann Myt;

public:
class Tr{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
//typedef LocTy location_type;
//typedef mjm_rational RatTy;
}; // Tr

typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<IdxTy> MyBlockInt;
typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
//typedef mjm_generic_itor<2> LoopItor2;
//typedef mjm_generic_itor<3> LoopItor3;
//typedef mjm_generic_itor<4> LoopItor4;

typedef mjm_fixed_sum_itor GenItor;
mjm_burmann () {}

~mjm_burmann () {}


////////////////////////////////////////////////////////////////////////////

// [1] http://www.mathematica-journal.com/2014/11/on-burmanns-theorem-and-its-application-to-problems-of-linear-and-nonlinear-heat-transfer-and-diffusion/#more-39602/ https://en.wikipedia.org/wiki/Error_function#cite_note-8
// http://mathworld.wolfram.com/BuermannsTheorem.html

template <class Tvec, class Tf>
 void burmann_ref1_eq21(Tvec & dest, Tf & func, const D & z0, const IdxTy mu, const IdxTy nu, const IdxTy kmax, const IdxTy flags=0)
{
m_cm.mark("eq21");
for (IdxTy k=0; k<kmax; ++k)
{
D bc=0;
m_cm.mark("genitor");
GenItor gi(k,mu,nu);
m_cm.cumdel("genitortotal", "genitor");
const IdxTy sz=gi.k(); 
while (gi.ok())
{
m_cm.mark("multi");
	D fac=gi.multiplicity();
m_cm.cumdel("multitotal", "multi");
	for (IdxTy s=1; s<=sz; ++s)
	{	
		const D rat=func.ratio(s,0,z0);
		const D f1=(rat/gi.factorial(s));
		// captured uin multiplicity 
		const D f= gi.power(f1,gi[s]); // /gi.factorial(gi[s]);
		//MM_MSG(MMPR(rat)<<MMPR4(fsz0,f0z0,f1,f)<<MMPR3(gi[s],s,gi.factorial(s)))
		fac*=f;
	}
	const IdxTy sgn=(k-gi.p0())&1;
	if (sgn!=0) bc-=fac; else bc+=fac;	
//	MM_MSG(MMPR(fac)<<MMPR4(bc,sgn,k,gi.p0()))
	gi.inc();
} // gi.ok()
dest.push_back(bc);
} // k 
m_cm.cumdel("eq21total", "eq21");
} // burmann_ref1_eq21
///////////////////////////////////////////////////////////

template <class Tvec, class Tf,class Tpf,class Tpf2>
 void burmann_ref1_eq14(Tvec & dest, Tf & func, Tpf & phifunc, Tpf2 & phistarfunc, const D & z0, const IdxTy nu, const IdxTy nmax, const IdxTy flags=0)
{
m_cm.mark("eq14");
// TODO just used as a source for power and factorial lol WTF 
m_cm.mark("genitor");
GenItor gi(0,0,0);
m_cm.cumdel("genitortotal", "genitor");
D c1=::pow(gi.factorial(nu+1)/abs(phifunc(nu+1,z0)),D(1)/D(nu+1));
//MM_MSG(MMPR4(c1,phifunc(nu+1,z0),z0,(1.0/(nu+1))))
//const D z0=0;
for (IdxTy n=1; n<nmax; ++n)
{
D bc=0;
const IdxTy rmax=n; // gi.size();
Tvec R;
const D r_mu=n;
const D r_nu=nu+1;
//MM_MSG(MMPR4( z0,  r_mu,  r_nu,  n))
burmann_ref1_eq21(R,  phistarfunc,  z0,  r_mu,  r_nu,  nmax+1, flags);
for(IdxTy r=0; r<rmax; ++r)
{
	D den=n*gi.factorial(r); // yes this can be done incrementally wtf 
	const D frz0=func(r+1,z0);
	bc+=frz0*R[n-r-1]/den;
	//MM_MSG(MMPR4(r,n,frz0,R[n-r-1])<<MMPR(den))
} // r
const D c=gi.power(c1,n);
const D bcc=bc*c;
//MM_MSG(" erf coef "<<MMPR4(n,c,bc,bcc))
dest.push_back(bcc);
} //n  
m_cm.cumdel("eq14total", "eq14");
} // burmann_ref1_e14q

// hard coded erf at z0=0
const bool odd(const IdxTy & x) { return ((x&1)!=0); } 
const bool odd(const int  & x) { return ((x&1)!=0); } 
const bool odd2(const IdxTy & x) { return ((x&2)!=0); } 
const bool odd2(const int  & x) { return ((x&2)!=0); } 


template <class Tvec>
 void burmann_ref1_eq14_erf(Tvec & dest,  const IdxTy nmax, const IdxTy flags=0)
{
m_cm.mark("eq14erf");
// 2m= n+1 from eqn 14, the odd terms are zero 
for (IdxTy m=1; m<nmax; ++m)
{
D bcc=0;
// this is n-1 or n-1=2m-2, 
const int  n=2*m-1;
const int  rlim=n-1; // 2*m-2;
for(int  r=0; r<=rlim; ++r)
{
// first get the f(r+1)/f! term, never mind see fac at end 
const D pref=1;
const int k=n-r-1; // 2*m-r-2;
const int  mu=n; // 2*m-1;
const int  nu=2;
D part=0;
GenItor gi(k,mu,nu);
// eqn 14
 int ord=r;
ord=ord/2;
D fac=0;
if (r==0) fac=1.0/D(n);
else if (!odd(r))  fac=2.0*(odd(ord)?-1:1)*gi.factorial(2*ord-1)/gi.factorial(ord-1)/(D(n))/gi.factorial(r);
//if (odd(r+1)) fac=0;

// this is all eqn 21i, mu is n, v+1=2
while (gi.ok())
{
const int  p0=gi.p0(); // gi[0];
const int qlim=k-p0;
D parity=odd(qlim)?-1:1;
const D qr=D(mu)/D(nu);
D qterm=1;
for (int  q=0; q<qlim; ++q)
{
qterm*=(q+qr);
} // q
const int  slim=k;
bool sodd=true;
D sterm=1;
for (int  s=1; s<=slim; ++s)
{
sodd=odd(s);
if (sodd) { if (gi[s]!=0) sterm=0; else sterm=sterm/gi.factorial(gi[s]); }
else{
 const D fsz0=(odd2(s)?-1.0:1.0)*gi.factorial(s)/gi.factorial((s>>1)+1);
 sterm*=gi.power(fsz0/gi.factorial(s),gi[s])/gi.factorial(gi[s]);
}

sodd=!sodd;
} // s
part+=parity*qterm*sterm;
gi.inc();
} // gi 
bcc+=fac*part*pref;
} // r

//MM_MSG(" erf coef "<<MMPR4(n,c,bc,bcc))
dest.push_back(bcc);
} //n  
m_cm.cumdel("eq14erftotal", "eq14erf");
} // burmann_ref1_e14q
//////////////////////////////////////////////////////

//CounterMap m_cm;

////////////////////////////////////////////////////////////////////////////
//}; // mjm_burmann
///////////////////////////////////////////////////////////////
/*

based on [1] eqn 14 and 21
*/

template <class Tvec, class Tfunc>
 void burmann_1(Tvec & dest,  Tfunc & f, const IdxTy nmax, const IdxTy flags=0)
{
typedef typename Tfunc::calc_type Rv;
typedef typename Tfunc::int_type It; // int
typedef typename Tfunc::index_type Itx; // int
typedef std::vector<IdxTy> ExcludeList;
typedef std::vector<Rv> ResultCache;
//m_cm.mark("eq14erf"); // put into f 
m_cm.mark(f.perf_name()); // put into f 
f.enter();

const Itx  realnmax=f.n_from_m(nmax); // 2*m-1;
// this may NOT always be better for fast calculations
// move to "f" 
ExcludeList el;
const bool tp=f.trim_partitions();
if (tp)
{
for( Itx i=0; i<=realnmax; ++i) {
		if ( Rv(0)==f.phistar(i)) {el.push_back(1) ;} 
		else { el.push_back(0); }  }
}
// 2m= n+1 from eqn 14, the odd terms are zero 
for (Itx m=1; m<Itx(nmax); ++m)
{
Rv bcc=0;
// this is n-1 or n-1=2m-2, 
const Itx  n=f.n_from_m(m); // 2*m-1;

const It  mu=n; 
const It  nu=f.nu()+1;  
const Rv qr=Rv(mu)/Rv(nu); // moved out of loop
Rv qcterm;
qcterm=1;
ResultCache qcache;
qcache.push_back(qcterm);
//for (Itx  q=0; q<qlim; ++q) { qterm=qterm*(Rv(q)+qr); } // q
while (It(qcache.size())<=n) 
{Itx q=qcache.size()-1; qcterm=qcache[q]*(Rv(q)+qr); qcache.push_back(qcterm);  }

const Itx  rlim=n-1; 
for(Itx  r=0; r<=rlim; ++r)
{
const Itx  k=n-r-1;  

Rv part=0;
// should come from "f" 
GenItor gi=tp?GenItor(IdxTy(k),IdxTy(mu),IdxTy(nu),el)
 				:GenItor(IdxTy(k),IdxTy(mu),IdxTy(nu));
// eqn 14
Rv fac=f.f(r+1)/(Rv(n)*f.factorial(r));
// this is all eqn 21, mu is n, v+1=2
while (gi.ok())
{
const It  p0=gi.p0();  
const Itx qlim=k-Itx(p0);
Rv parity=odd(qlim)?-1:1;
//const Rv qr=Rv(mu)/Rv(nu); // move out of loop
// right now the rational class lacks some ctor's and operators etc. 
// This should call a factory in or "ctor" in f to get options right for
// things like reducing and policies etc
//Rv qterm;
//qterm=1;
//const Rv & qterm=qcache[qlim];
// these can be cached 
//for (Itx  q=0; q<qlim; ++q) { qterm=qterm*(Rv(q)+qr); } // q
const Itx  slim=k;
//Rv sterm;
 Rv  sterm=qcache[qlim];
//sterm=1;
for (Itx  s=1; s<=slim; ++s)
{
 const Rv fsz0=f.phistar(s) ;
 sterm=sterm*f.power_t(fsz0/f.factorial(s),gi[s])/f.factorial(gi[s]);
} // s
//part=part+parity*qterm*sterm;
part=part+parity*sterm;
gi.inc();
} // gi 
bcc=bcc+fac*part;
} // r

//MM_MSG(" erf coef "<<MMPR4(n,c,bc,bcc))
dest.push_back(bcc);
} //n  
f.exit();
//m_cm.mark(f.perf_name()); // put into f 
//f.enter();
//m_cm.cumdel("eq14erftotal", "eq14erf");
m_cm.cumdel(f.perf_total(), f.perf_name());
} // burmann_1
//////////////////////////////////////////////////////
/*
based on [1] eqn 1,5, and 6
AFAICT this does not work 
*/
template <class Tvec, class Tfunc>
 void burmann_2(Tvec & dest,  Tfunc & f, const IdxTy nmax, const IdxTy flags=0)
{
typedef typename Tfunc::calc_type Rv;
//typedef typename Tfunc::int_type It; // int
typedef typename Tfunc::index_type Itx; // int
//typedef std::vector<IdxTy> ExcludeList;
//typedef std::vector<Rv> ResultCache;
//m_cm.mark("eq14erf"); // put into f 
m_cm.mark(f.perf_name()); // put into f 
f.enter();
Rv bcc=f.f(1)/f.phi(1);
MM_MSG(MMPR(bcc.to_string()))

dest.push_back(f.evaluate(bcc));
const Itx  realnmax=f.n_from_m(nmax); // 2*m-1;
for (Itx m=1; m<Itx(nmax); ++m)
{
const Itx  n=f.n_from_m(m); // 2*m-1;
// this is just eqn 5 but the symbolic math needs to be determined by f
MM_MSG(MMPR(bcc.to_string()))
MM_MSG(MMPR((f.der(bcc)).to_string()))
MM_MSG(MMPR2(n,(f.phi(1)).to_string()))
MM_MSG(MMPR((f.phi(1)*n).to_string()))
bcc=f.der(bcc)/(f.phi(1)*n);
//MM_MSG(MMPR3(f.der(bcc).to_string(),( f.phi(1)*n).to_string() , bcc.to_string()))
MM_MSG(MMPR(bcc.to_string()))

//MM_MSG(" erf coef "<<MMPR4(n,c,bc,bcc))
dest.push_back(f.evaluate(bcc));
} //n  
f.exit();
m_cm.cumdel(f.perf_total(), f.perf_name());
} // burmann_2

/*
find polynomial ratios in delta for erf from fd equation
*/

template <class Ty> void erf_fd_poly(Ty & dest, const IdxTy n)
{
typedef mjm_rational RatTy;
//typedef RatTy::IntTy IntTy;
typedef RatTy Fty;
//typedef double Fty;
typedef mjm_closed_families::exp_poly_ratio<Fty> Nr;
//typedef std::vector<Fty> Cache;
//typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<Fty> Po;
//typedef std::vector<Nr> Dstack;
Fty one=Fty(1);
one.set_reduce();
Fty mone=Fty(-1);
mone.set_reduce();

Fty zzz=Fty(0);
zzz.set_reduce();

Nr minusone=Nr();
Po p=Po();
// these are polynomials in detlta, the minus one value is -delta
//p.push_back(Fty(0));
p.push_back(zzz);
//p.push_back(Fty(-1));
p.push_back(mone);
minusone.poly(p);
Nr zed=Nr();
p[1]=zzz; // 0;
zed.poly(p);
//Nr y=zed;
for(IdxTy i=0; i<n; ++i)
{
dest.push_back(zed);
// so the next term is x= (i+1)*delta
Po p=Po();
Po pn=Po();
//p.push_back(int(1));
p.push_back(one);
p.push_back(Fty(i));
pn.push_back(int(2));
Nr pf=Nr();
pf.poly(pn,p);
//MM_MSG(MMPR4(p[0],p[1],pn[0],pf.to_string()))
Nr next=minusone+ pf*(zed-minusone);
//MM_MSG(MMPR4(next.to_string(),pf.to_string(),zed.to_string(),minusone.to_string()))
//MM_MSG(MMPR(next.to_string()))
minusone=zed;
zed=next;
}
}
///////////////////////////////////////////////////

template <class Ty,class Trat > void erf_fd_rat(Ty & dest, const IdxTy n, const Trat & delta)
{
typedef mjm_rational RatTy;
typedef RatTy Fty;
//typedef shape_function_integrals SFi;
//typedef SFi::simple_polynomial<Fty> Po;
//typedef std::vector<Nr> Dstack;
Fty two=Fty(2);
two.set_reduce();
Fty one=Fty(1);
one.set_reduce();
Fty mone=Fty(-1);
mone.set_reduce();

//Fty zzz=Fty(0);
//zzz.set_reduce();
Fty  minusone=-delta;
minusone.set_reduce();
Fty zed=Fty(0);
zed.set_reduce();
for(IdxTy i=0; i<n; ++i)
{
dest.push_back(zed);
// so the next term is x= (i+1)*delta
Fty pf=two/(one+delta*i);
pf.set_reduce();
//MM_MSG(MMPR4(p[0],p[1],pn[0],pf.to_string()))
Fty next=minusone+ pf*(zed-minusone);
//MM_MSG(MMPR4(next.to_string(),pf.to_string(),zed.to_string(),minusone.to_string()))
//MM_MSG(MMPR(next.to_string()))
//MM_MSG(MMPR4(i,p.to_string(),v,double(v))<<MMPR((erf(z2*i)*0.8862269))) 
MM_MSG(MMPR4(i,next,double(next),(double(next)/.8862269))<<MMPR(erf(double(delta*(i+1)))))
minusone=zed;
zed=next;
}
}



//////////////////////////////////////////////////////
class erf_demo
{
typedef erf_demo Myt;
typedef mjm_burmann Tr;
typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
// TODO FIXME this is a kluge for factorial and power right now really dumb
//typedef Tr::GenItor GenItor;

typedef D Fty;
typedef std::vector<Fty> Cache;

public:
typedef Fty calc_type;
typedef int index_type;
typedef int int_type;
erf_demo(): m_z0(0) // ,m_gi(0,0,0)
{

}
// actually twice as fast when doing same work lol 
//int n_from_m(const IdxTy m ) { return 2*m-1;}
int n_from_m(const IdxTy m ) { return m;}
// this is defininitely slower 
bool trim_partitions() { return false; } 

int nu( ) { return 1;}

// these are always cached so ref is not a temp 
// these are deriviative evealuated at z0
const Fty & f(const IdxTy n  ) { Makef(n); return m_f[n]; }
const Fty & phi(const IdxTy n  ) { Makephi(n); return m_phi[n]; }
const Fty & phistar(const IdxTy n  ) { Makephistar(n); return m_phi_star[n]; }
const Fty & factorial(const IdxTy n) { return Fact(n); } 
// for perf monitor 
void enter() {}
void exit() {}

//m_cm.cumdel("eq14erftotal", "eq14erf");
StrTy perf_name() const { return "erf_demo"; }
StrTy perf_total() const { return "erf_demototal"; }

template <class Tv> Tv power_t( const Tv & b,const IdxTy & e) {
Tv v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; }

 return v; }




private:
void Makef(const IdxTy n)
{
while (m_f.size()<=n)
{
const IdxTy _sz=m_f.size();
// this is not quite right but does work with r+ 1
if (_sz==0) {m_f.push_back(1); continue; } 
if (_sz==1) {m_f.push_back(1); continue; } 
const IdxTy sz=_sz-1;
const bool odd=((sz&1)!=0);
if (odd) {m_f.push_back(0); continue; }
const IdxTy ord=sz>>1;
const bool oddo=((ord&1)!=0);
//const Fty  fac=2.0*(oddo?-1.0:1.0)*m_gi.factorial(2*ord-1)/m_gi.factorial(ord-1);
const Fty  fac=2.0*(oddo?-1.0:1.0)*Fact(2*ord-1)/Fact(ord-1);
m_f.push_back(fac);
}

} //Makef
void Makephi(const IdxTy n)
{

}
void Makephistar(const IdxTy n)
{
while (m_phi_star.size()<=n)
{
const IdxTy _sz=m_phi_star.size();
const IdxTy sz=_sz-0; // kluge for the code dopy lol
const bool odd=((sz&1)!=0);
if (odd) { m_phi_star.push_back(0); continue; } 
const bool odd2=((sz&2)!=0);
 //const D fsz0= (odd2?-1.0:1.0)*m_gi.factorial(sz)/m_gi.factorial((sz>>1)+1);
 const D fsz0= (odd2?-1.0:1.0)*Fact(sz)/Fact((sz>>1)+1);
m_phi_star.push_back(fsz0);
}

} //Makephistar
const Fty &  Fact(const IdxTy n)
{
while (m_fact.size()<=n)
{
const IdxTy sz=m_fact.size();
if (sz==0) { m_fact.push_back(1); continue; } 
m_fact.push_back(m_fact[sz-1]*sz);
}
return m_fact[n];
}

Cache m_f,m_phi,m_phi_star,m_fact;
Fty m_z0;
//GenItor m_gi;
}; // erf_demo
////////////////////////////////////////////////////////////////

class erf_demo_alg
{
typedef erf_demo_alg Myt;
typedef mjm_burmann Tr;
typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
// TODO FIXME this is a kluge for factorial and power right now really dumb
//typedef Tr::GenItor GenItor;

typedef D Fty;
typedef std::vector<Fty> Cache;

typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
typedef SFi::norm_rat<D> Nr;
typedef std::vector<Nr> Dstack;
typedef std::vector<D> Vstack;

typedef std::vector<Po> Fstack;
public:
typedef Fty calc_type;
typedef int index_type;
typedef int int_type;
erf_demo_alg(): m_z0(0) { }
erf_demo_alg(const D & z0): m_z0(z0) { }
// actually twice as fast when doing same work lol 
//int n_from_m(const IdxTy m ) { return 2*m-1;}
int n_from_m(const IdxTy m ) { return m;}
bool trim_partitions() { return !false; } 
int nu( ) {if (m_z0==0)  return 1; return 0;}

// these are always cached so ref is not a temp 
// these are deriviative evealuated at z0
const Fty & f(const IdxTy n  ) { Makef(n); return m_f[n]; }
const Fty & phi(const IdxTy n  ) { Makephi(n); return m_phi[n]; }
const Fty & phistar(const IdxTy n  ) { Makephistar(n); return m_phi_star[n]; }
const Fty & factorial(const IdxTy n) { return Fact(n); } 
// for perf monitor 
void enter() {}
void exit() {}

//m_cm.cumdel("eq14erftotal", "eq14erf");
StrTy perf_name() const { return "erf_demoalg"; }
StrTy perf_total() const { return "erf_demoalgtotal"; }

template <class Tv> Tv power_t( const Tv & b,const IdxTy & e) {
Tv v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; }

 return v; }




private:
void Makef(const IdxTy n)
{
while (m_f.size()<=n)
{
const IdxTy _sz=m_f.size();
// this is not quite right but does work with r+ 1
// the first one is a dummy, the next is the first gaussiant poly 
//if (m_der.size()==0){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 
//if (m_der.size()==1){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 

if (_sz==0) {Po crap; crap.push_back(1); m_f_der.push_back(crap); m_f.push_back(::sqrt(M_PI)*.5*erf(m_z0)); continue; } 
if (_sz==1) {Po crap; crap.push_back(1); m_f_der.push_back(crap); m_f.push_back(exp(-m_z0*m_z0)); continue; } 
//if (_sz==1) {m_f.push_back(1); continue; } 
/*
const IdxTy sz=_sz-1;
const bool odd=((sz&1)!=0);
if (odd) {m_f.push_back(0); continue; }
const IdxTy ord=sz>>1;
const bool oddo=((ord&1)!=0);
//const Fty  fac=2.0*(oddo?-1.0:1.0)*m_gi.factorial(2*ord-1)/m_gi.factorial(ord-1);
const Fty  fac=2.0*(oddo?-1.0:1.0)*Fact(2*ord-1)/Fact(ord-1);
*/
//while (m_der.size()<=n) {
Po p;
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
mi.next_gauss_der_polynomial(p, m_f_der[m_f_der.size()-1]);
m_f_der.push_back(p);
m_f.push_back(mi.evaluate_polynomial(p,m_z0));
}
// this is the polynomial
//m_f_der.push_back();
// this is teh value at z0
//m_f.push_back(fac);

} //Makef



void Makephi(const IdxTy n)
{

}
void Makephistar(const IdxTy n)
{


/*
while (m_phi_star.size()<=n)
{
const IdxTy _sz=m_phi_star.size();
const IdxTy sz=_sz-0; // kluge for the code dopy lol
const bool odd=((sz&1)!=0);
if (odd) { m_phi_star.push_back(0); continue; } 
const bool odd2=((sz&2)!=0);
 //const D fsz0= (odd2?-1.0:1.0)*m_gi.factorial(sz)/m_gi.factorial((sz>>1)+1);
 const D fsz0= (odd2?-1.0:1.0)*Fact(sz)/Fact((sz>>1)+1);
*/
///////////////////////////////////////////////
	if (m_der.size()==0)
	{ 
		Nr crap; 
		if (m_z0==0) crap.norm_minus_1_over_xn(2);  
		else  crap.norm_minus_1_over_xn(1,m_z0);  

		m_der.push_back(crap);
		m_phi_star.push_back(1);
//		m_ratio_cache.push_back(1);
 	} 
	while (m_der.size()<=n) 
	{
		Nr p; Po n2,d2;
		//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
 		m_der[m_der.size()-1].der(p,n2,d2);
		MM_MSG(" adding der phistar "<<MMPR2(m_der.size(),p.to_string()))
		m_der.push_back(p);
///////////////////////////////////////////
		//Nr &  p=m_der[n+m_offset];
		// l'hopital?  mais no moise numerical est le merde wtf??? 
		const D fsz0=(p/m_der[0]).evaluate(m_z0); // evaluate_at_zero();
//		const D fsz0=p.evaluate(m_z0)/m_phi_star[0];
//MM_MSG(MMPR4(x,n,z,p.to_string()))
		m_phi_star.push_back(fsz0);
	}

} //Makephistar
const Fty &  Fact(const IdxTy n)
{
while (m_fact.size()<=n)
{
const IdxTy sz=m_fact.size();
if (sz==0) { m_fact.push_back(1); continue; } 
m_fact.push_back(m_fact[sz-1]*sz);
}
return m_fact[n];
}

Cache m_f,m_phi,m_phi_star,m_fact;
Fty m_z0;

Dstack m_der;
Fstack m_f_der;
//Vstack m_ratio_cache;
Mi mi;

//GenItor m_gi;
}; // erf_demo_alg



//////////////////////////////////////////////

class erf_demo_alg_rat
{
typedef erf_demo_alg_rat Myt;
typedef mjm_burmann Tr;
typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;

typedef mjm_rational RatTy;
typedef RatTy::IntTy IntTy;
//typedef D Fty;
//typedef std::vector<Fty> Cache;

typedef RatTy Fty;
typedef std::vector<Fty> Cache;

typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<Fty> Po;
typedef SFi::norm_rat<Fty> Nr;
typedef std::vector<Nr> Dstack;
typedef std::vector<D> Vstack;

typedef std::vector<Po> Fstack;
public:
typedef Fty calc_type;
typedef int index_type;
typedef mjm_bigint_mmpr int_type;
//typedef int int_type;
erf_demo_alg_rat(): m_z0(0),m_exp_z0z0(1) {init(); }
erf_demo_alg_rat(const D & z0): m_z0(z0),m_exp_z0z0(D(exp(-double(m_z0*m_z0)))) { init();}
erf_demo_alg_rat(const int & z0n,const int & z0d): m_z0(z0n,z0d),m_exp_z0z0(D(exp(-double(m_z0*m_z0)))) { init();}

void init()
{
m_z0.simplify();
m_exp_z0z0.simplify();

}
//erf_demo_alg_rat(const double & z0): m_z0(z0),m_exp_z0z0(D(exp(-double(m_z0*m_z0)))) { }
// actually twice as fast when doing same work lol 
//int n_from_m(const IdxTy m ) { return 2*m-1;}
int_type  n_from_m(const IdxTy m ) { return IntTy(m);}
bool trim_partitions() { return !false; } 
int_type nu( ) {if (m_z0==0)  return IntTy(1); return IntTy(0);}

// these are always cached so ref is not a temp 
// these are deriviative evealuated at z0
const Fty & f(const IdxTy n  ) { Makef(n); return m_f[n]; }
const Fty & phi(const IdxTy n  ) { Makephi(n); return m_phi[n]; }
const Fty & phistar(const IdxTy n  ) { Makephistar(n); return m_phi_star[n]; }
const Fty & factorial(const IdxTy n) { return Fact(n); } 
// for perf monitor 
void enter() {}
void exit() {}

//m_cm.cumdel("eq14erftotal", "eq14erf");
StrTy perf_name() const { return "erf_demoalgrat"; }
StrTy perf_total() const { return "erf_demoalgrattotal"; }

template <class Tv> Tv power_t( const Tv & b,const IdxTy & e) {
Tv v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; }
v.simplify();

 return v; }




private:
void Makef(const IdxTy n)
{
while (m_f.size()<=n)
{
const IdxTy _sz=m_f.size();
// this is not quite right but does work with r+ 1
// the first one is a dummy, the next is the first gaussiant poly 
//if (m_der.size()==0){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 
//if (m_der.size()==1){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 

if (_sz==0) {Po crap; crap.push_back(1); m_f_der.push_back(crap); m_f.push_back(::sqrt(M_PI)*.5*erf(double(m_z0))); continue; } 
if (_sz==1) {Po crap; crap.push_back(1); m_f_der.push_back(crap); m_f.push_back(m_exp_z0z0); continue; } 
//if (_sz==1) {m_f.push_back(1); continue; } 
/*
const IdxTy sz=_sz-1;
const bool odd=((sz&1)!=0);
if (odd) {m_f.push_back(0); continue; }
const IdxTy ord=sz>>1;
const bool oddo=((ord&1)!=0);
//const Fty  fac=2.0*(oddo?-1.0:1.0)*m_gi.factorial(2*ord-1)/m_gi.factorial(ord-1);
const Fty  fac=2.0*(oddo?-1.0:1.0)*Fact(2*ord-1)/Fact(ord-1);
*/
//while (m_der.size()<=n) {
Po p;
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
mi.next_gauss_der_polynomial(p, m_f_der[m_f_der.size()-1]);
m_f_der.push_back(p);
m_f.push_back(mi.evaluate_polynomial(p,m_z0));
}
// this is the polynomial
//m_f_der.push_back();
// this is teh value at z0
//m_f.push_back(fac);

} //Makef



void Makephi(const IdxTy n)
{

}
void Makephistar(const IdxTy n)
{


/*
while (m_phi_star.size()<=n)
{
const IdxTy _sz=m_phi_star.size();
const IdxTy sz=_sz-0; // kluge for the code dopy lol
const bool odd=((sz&1)!=0);
if (odd) { m_phi_star.push_back(0); continue; } 
const bool odd2=((sz&2)!=0);
 //const D fsz0= (odd2?-1.0:1.0)*m_gi.factorial(sz)/m_gi.factorial((sz>>1)+1);
 const D fsz0= (odd2?-1.0:1.0)*Fact(sz)/Fact((sz>>1)+1);
*/
///////////////////////////////////////////////
	if (m_der.size()==0)
	{ 
		Nr crap; 
		if (m_z0==0) crap.norm_minus_1_over_xn(2);  
		else  crap.norm_minus_1_over_xn_exact(1,m_exp_z0z0,-m_z0);  
//		crap.simplify();
		m_der.push_back(crap);
		m_phi_star.push_back(1);
//		m_ratio_cache.push_back(1);
 	} 
	while (m_der.size()<=n) 
	{
		Nr p; Po n2,d2;
		//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
 		m_der[m_der.size()-1].der(p,n2,d2);
		//p.simplify();
		MM_MSG(" adding der phistar "<<MMPR2(m_der.size(),p.to_string()))
		m_der.push_back(p);
///////////////////////////////////////////
		//Nr &  p=m_der[n+m_offset];
		// l'hopital?  mais no moise numerical est le merde wtf??? 
		const Fty  fsz0=(p/m_der[0]).evaluate(m_z0,m_exp_z0z0); // evaluate_at_zero();
//		const D fsz0=p.evaluate(m_z0)/m_phi_star[0];
//MM_MSG(MMPR4(x,n,z,p.to_string()))
		m_phi_star.push_back(fsz0);
	}

} //Makephistar
const Fty &  Fact(const IdxTy n)
{
while (m_fact.size()<=n)
{
const IdxTy sz=m_fact.size();
if (sz==0) { m_fact.push_back(1); continue; } 
m_fact.push_back(m_fact[sz-1]*sz);
}
return m_fact[n];
}

Cache m_f,m_phi,m_phi_star,m_fact;
Fty m_z0,m_exp_z0z0;

Dstack m_der;
Fstack m_f_der;
//Vstack m_ratio_cache;
Mi mi;

//GenItor m_gi;
}; // erf_demo_alg



//////////////////////////////////////////////




//////////////////////////////////////////////




////////////////////////////////////////////////////////////////
// for burnmann_2
// this goes out of the domain of (Pexp+q)/(Rexp+s)
// needs general polynomial of (exp(-x*x)^n accepting rationals
// and preserviing roots lol 
class erf_demo2_alg
{
typedef erf_demo2_alg Myt;
typedef mjm_burmann Tr;
typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;

//typedef D Fty;
//typedef std::vector<Fty> Cache;

typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
// these are ratio of (P(x)*exp(-x*x)+Q(x)/(R(x)*exp()+Q(x)
// note this is not closed under some operations as exp(-2x*x) terms result soeneed care 
//typedef SFi::norm_rat<D> Nr;
typedef mjm_closed_families::exp_poly_ratio<D> Nr;
typedef Nr Fty;
typedef std::vector<Fty> Cache;



typedef std::vector<Nr> Dstack;
typedef std::vector<D> Vstack;

typedef std::vector<Po> Fstack;
public:
typedef Fty calc_type;
typedef D eval_type;
typedef int index_type;
typedef int int_type;
erf_demo2_alg(): m_z0(0) { }
erf_demo2_alg(const D & z0): m_z0(z0) { }
// actually twice as fast when doing same work lol 
//int n_from_m(const IdxTy m ) { return 2*m-1;}
int n_from_m(const IdxTy m ) { return m;}
bool trim_partitions() { return !false; } 
int nu( ) {if (m_z0==0)  return 1; return 0;}

// these are always cached so ref is not a temp 
// these are deriviative evealuated at z0
// user may simplofy aarrrrggghhh 
//const 
Fty & f(const IdxTy n  ) { Makef(n); return m_f[n]; }
//const 
Fty & phi(const IdxTy n  ) { Makephi(n); return m_phi[n]; }
//const Fty & phistar(const IdxTy n  ) { Makephistar(n); return m_phi_star[n]; }
//const Fty & factorial(const IdxTy n) { return Fact(n); } 
Fty  der( Fty & x) { return Der(x); } 
const eval_type evaluate(Fty & x) { return x.evaluate(m_z0); } 
// for perf monitor 
void enter() {}
void exit() {}

//m_cm.cumdel("eq14erftotal", "eq14erf");
StrTy perf_name() const { return "erf_demo2alg"; }
StrTy perf_total() const { return "erf_demo2algtotal"; }


private:
//const Fty  Der( Fty & x) { Fty p; Po n,d; x.der(p,n,d); return (p); } 
const Fty  Der( Fty & x) { Fty p;  x.der(p); return (p); } 
void Makef(const IdxTy n)
{
while (m_f.size()<=n)
{
const IdxTy _sz=m_f.size();

//if (_sz==0) {Fty crap; crap.push_back(1); m_f_der.push_back(crap); m_f.push_back(::sqrt(M_PI)*.5*erf(m_z0)); continue; } 
if (_sz==0) {  m_f.push_back(Fty(::sqrt(M_PI)*.5*erf(m_z0))); continue; } 
if (_sz==1) {Fty crap; crap.normal();  m_f.push_back(crap); continue; } 
//Po p;
Fty p;
//Po n,d;
//m_f[m_f.size()-1].der(p,n,d);
m_f[m_f.size()-1].der(p);
MM_MSG(MMPR2(m_f.size(),p.to_string()))
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
//mi.next_gauss_der_polynomial(p, m_f_der[m_f_der.size()-1]);
m_f.push_back(p);
//m_f.push_back(mi.evaluate_polynomial(p,m_z0));
}
} //Makef



void Makephi(const IdxTy n)
{

while (m_phi.size()<=n)
{
const IdxTy _sz=m_phi.size();
MM_MSG(MMPR2(m_phi.size(),f(_sz+1).to_string()))
if (_sz==0)
{
Fty x;
x=f(1);
//void add_single_term(const IdxTy  expn, const IdxTy p, const IdxTy flag, const value_type & c)
//{ PolyPoly & x=((flag&1)==0)?m_n:m_d; x.add_single_term(expn,p,flag,c);

x.add_single_term(0,0,0,-x.evaluate(m_z0));
MM_MSG(MMPR2(m_phi.size(),x.to_string()))

m_phi.push_back(x);
continue;
}
m_phi.push_back(f(_sz+1));
}

}


void Makephistar(const IdxTy n)
{

/*
///////////////////////////////////////////////
	if (m_der.size()==0)
	{ 
		Nr crap; 
		if (m_z0==0) crap.norm_minus_1_over_xn(2);  
		else  crap.norm_minus_1_over_xn(1,m_z0);  

		m_der.push_back(crap);
		m_phi_star.push_back(1);
//		m_ratio_cache.push_back(1);
 	} 
	while (m_der.size()<=n) 
	{
		Nr p; Po n2,d2;
		//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
 		m_der[m_der.size()-1].der(p,n2,d2);
		m_der.push_back(p);
///////////////////////////////////////////
		//Nr &  p=m_der[n+m_offset];
		// l'hopital? 
		const D fsz0=(p/m_der[0]).evaluate(m_z0); // evaluate_at_zero();
//		const D fsz0=p.evaluate(m_z0)/m_phi_star[0];
//MM_MSG(MMPR4(x,n,z,p.to_string()))
		m_phi_star.push_back(fsz0);
	}
*/

} //Makephistar
/*
const Fty &  Fact(const IdxTy n)
{
while (m_fact.size()<=n)
{
const IdxTy sz=m_fact.size();
if (sz==0) { m_fact.push_back(1); continue; } 
m_fact.push_back(m_fact[sz-1]*sz);
}
return m_fact[n];
}
*/
///Cache m_f,m_phi,m_phi_star,m_fact;
Cache m_f,m_phi;
eval_type m_z0;

//Dstack m_der;
//Fstack m_f_der;
Mi mi;

}; // erf_demo2_alg






/////////////////////////////////////////////////////////////////////

class erf_demo_rat
{
typedef erf_demo Myt;
typedef mjm_burmann Tr;
typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;

typedef mjm_rational RatTy;
typedef RatTy::IntTy IntTy;
//typedef std::vector<RatTy> Bk;

// TODO FIXME this is a kluge for factorial and power right now really dumb
//typedef Tr::GenItor GenItor;

typedef RatTy Fty;
typedef std::vector<Fty> Cache;

public:
typedef Fty calc_type;
typedef int index_type;
typedef mjm_bigint_mmpr int_type;
erf_demo_rat(): m_z0(D(0)) // ,m_gi(0,0,0)
{

}
StrTy perf_name() const { return "erf_demorat"; }
StrTy perf_total() const { return "erf_demorattotal"; }

// actually twice as fast when doing same work lol 
//int n_from_m(const IdxTy m ) { return 2*m-1;}
// needs template parameters etc 
index_type n_from_m(const IdxTy m ) { return m;}
//index_type n_from_m(const IdxTy m ) { return 2*m-1;}
int_type nu( ) { return IntTy(1);}
bool trim_partitions() { return !false; } 

// these are always cached so ref is not a temp 
// these are deriviative evealuated at z0
// for flat out speed this lazy thing kills the float version
// need to do this on enter or before 
const Fty & f(const IdxTy n  ) { Makef(n); return m_f[n]; }
const Fty & phi(const IdxTy n  ) { Makephi(n); return m_phi[n]; }
const Fty & phistar(const IdxTy n  ) { Makephistar(n); return m_phi_star[n]; }
const Fty & factorial(const IdxTy n) { return Fact(n); } 
// for perf monitor 
void enter() {}
void exit() {}
template <class Tv> Tv power_t( const Tv & b,const IdxTy & e) {
Tv v=1;
IdxTy j=e;
while (j!=0) { v=v*b; --j; }
v.simplify();
 return v; }


private:
void Makef(const IdxTy n)
{
while (m_f.size()<=n)
{
const IdxTy _sz=m_f.size();
// this is not quite right but does work with r+ 1
if (_sz==0) {m_f.push_back(IntTy(1)); continue; } 
if (_sz==1) {m_f.push_back(IntTy(1)); continue; } 
const IdxTy sz=_sz-1;
const bool odd=((sz&1)!=0);
if (odd) {m_f.push_back(0); continue; }
const IdxTy ord=sz>>1;
const bool oddo=((ord&1)!=0);
//const Fty  fac=2.0*(oddo?-1.0:1.0)*m_gi.factorial(2*ord-1)/m_gi.factorial(ord-1);
//const Fty  fac=(oddo?-2.0:2.0)*Fact(2*ord-1)/Fact(ord-1);
 Fty  fac=Fact(2*ord-1)/Fact(ord-1)*(oddo?-2:2);
fac.set_reduce();
m_f.push_back(fac);
}

} //Makef
void Makephi(const IdxTy n)
{

}
void Makephistar(const IdxTy n)
{
while (m_phi_star.size()<=n)
{
const IdxTy _sz=m_phi_star.size();
const IdxTy sz=_sz-0; // kluge for the code dopy lol
const bool odd=((sz&1)!=0);
if (odd) { m_phi_star.push_back(IntTy(0)); continue; } 
const bool odd2=((sz&2)!=0);
 //const D fsz0= (odd2?-1.0:1.0)*m_gi.factorial(sz)/m_gi.factorial((sz>>1)+1);
 //const D fsz0= (odd2?-1.0:1.0)*Fact(sz)/Fact((sz>>1)+1);
 Fty fsz0= Fact(sz)/Fact((sz>>1)+1)*(odd2?-1:1);
fsz0.set_reduce();
m_phi_star.push_back(fsz0);
}

} //Makephistar
const Fty &  Fact(const IdxTy n)
{
while (m_fact.size()<=n)
{
const IdxTy sz=m_fact.size();
if (sz==0) { m_fact.push_back(IntTy(1)); m_fact[0].set_reduce(); continue; } 
m_fact.push_back(m_fact[sz-1]*sz);
}
return m_fact[n];
}

Cache m_f,m_phi,m_phi_star,m_fact;
Fty m_z0;
//GenItor m_gi;
}; // erf_demo_rat





////////////////////////////////////////////////////////////////////
template <class Ts> void dump_counter(const Ts & lbl )
{
    m_cm.dump(lbl , std::cout);
    m_cm.dump(lbl , std::cerr);
    MM_MSG(lbl<<CRLF<<m_cm.time_table("solver_times"))
    MM_ERR(lbl<<CRLF<<m_cm.time_table("solver_times"))
}

CounterMap m_cm;

////////////////////////////////////////////////////////////////////////////
}; // mjm_burmann




//////////////////////////////////////////////////////////



#ifdef  TEST_BURMANN2__ 
class func_class
{
//typedef mjm_fixed_sum_itor GenItor;
typedef func_class Myt;
typedef mjm_integrals Mi;
typedef unsigned int IdxTy;
typedef double D ;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
typedef std::vector<Po> Dstack;
public: 
func_class(): m_one(false) {}
D operator()(const IdxTy n, const D & z) //const
{
if (m_one) return 1.0; 
make_der(n);
const Po &  p=m_der[n];;
const D x=mi.evaluate_polynomial(p,z);
return x; 
}
D ratio(const IdxTy n, const IdxTy d, const D & z)
{ Myt & x=(*this); return x(n,z)/x(0,z); }

void make_der(const IdxTy n)
{
while (m_der.size()<=n) 
{
Po p;
mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
m_der.push_back(p);
}

}
void set_example1()
{
m_one=false;
m_poly= Po(4);
m_poly[0]=(1);
m_poly[1]=(7);
m_poly[2]=(8);
m_poly[3]=(11);
m_der.clear();
m_der.push_back(m_poly);
}

void set_unity()
{
m_one=true;
}
void set_1n()
{
m_poly.clear();
for(IdxTy i=0; i<10; ++i) m_poly.push_back(1);

m_der.clear();
m_der.push_back(m_poly);
}


Po m_poly;
Dstack m_der;
Mi mi;
bool m_one;
}; // func_class
////////////////////////////////////////////////////////
class func_erf_class
{
//typedef mjm_fixed_sum_itor GenItor;
typedef func_erf_class Myt;
typedef mjm_integrals Mi;
typedef unsigned int IdxTy;
typedef double D ;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
typedef std::vector<Po> Dstack;
public: 
func_erf_class(): m_one(false) {}
D operator()(const IdxTy n, const D& z) //const
{
//if (m_one) return 1.0; 
static const D fa=sqrt(M_PI)*.5;
if (n==0) return fa*erf(z);
make_der(n);
const Po &  p=m_der[n];
const D x=mi.evaluate_polynomial(p,z)*exp(-z*z);
return x; 
}
D ratio(const IdxTy n, const IdxTy d, const D & z)
{ Myt & x=(*this); return x(n,z)/x(0,z); }

void make_der(const IdxTy n)
{
// the first one is a dummy, the next is the first gaussiant poly 
if (m_der.size()==0){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 
if (m_der.size()==1){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 
while (m_der.size()<=n) 
{
Po p;
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
mi.next_gauss_der_polynomial(p, m_der[m_der.size()-1]);
m_der.push_back(p);
}

}
void set_example1() { }

void set_unity() {  }
void set_1n() { }


Po m_poly;
Dstack m_der;
Mi mi;
bool m_one;
}; // 

class func_erfphistar_class
{
//typedef mjm_fixed_sum_itor GenItor;
typedef func_erfphistar_class Myt;
typedef mjm_integrals Mi;
typedef unsigned int IdxTy;
typedef double D ;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
typedef SFi::norm_rat<D> Nr;
//void norm_minus_1_over_xn(const IdxTy n)
typedef std::vector<Nr> Dstack;
typedef std::vector<D> Vstack;
public: 
func_erfphistar_class(): m_one(false),m_offset(0) {}
const IdxTy & offset() const { return m_offset;}
const IdxTy & offset(const IdxTy off)  { m_offset=off; return m_offset;}

D operator()(const IdxTy n, const D& z) //const
{
//if (m_one) return 1.0; 
//static const D fa=sqrt(M_PI)*.5;
//if (n==0) return fa*erf(z);
make_der(n+m_offset);
Nr &  p=m_der[n+m_offset];
const D x=p.evaluate(z);
//MM_MSG(MMPR4(x,n,z,p.to_string()))
return x; 
}
D ratio(const IdxTy n, const IdxTy d, const D & z)
{//  Myt & x=(*this); 
make_der(n+m_offset);
make_der(d+m_offset);
if (z==0) if (d==0) if (m_offset==0) return m_ratio_cache[n];
Nr &  pn=m_der[n+m_offset];
Nr &  pd=m_der[d+m_offset];
const D x=(z==0)?((pn/pd).evaluate_at_zero()):((pn/pd).evaluate(z));
//MM_MSG(MMPR4(x,z,n,d)<<MMPR3(pn.to_string(),pd.to_string(),(pn/pd).to_string()))

return x; 

}

void make_der(const IdxTy n)
{
//void norm_minus_1_over_xn(const IdxTy n)
if (m_der.size()==0)
	{ Nr crap; crap.norm_minus_1_over_xn(2);  m_der.push_back(crap);
m_ratio_cache.push_back(1);


 } 
while (m_der.size()<=n) 
{
Nr p;
Po n2,d2;
bool dummy_thing=true; 
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
if (!dummy_thing) m_der[m_der.size()-1].der(p,n2,d2);
else MM_ONCE(" warning only the z=0 thing is evaluated for speed",)
m_der.push_back(p);
// TODO FIXME this ratio ignored m_offset
// STILL slow
//const D x=(p/m_der[0]).evaluate_at_zero();
D x=0;
const IdxTy neff=m_ratio_cache.size();
bool odd=((neff&1)!=0);
if (!odd)
{
 // (2n+1)!/(n+1)!
const IdxTy n2=neff>>1;
const IdxTy ii=n2+1;
bool odd2=((ii&1)!=0);
if (!odd2) x=-1; else x=1;
const IdxTy jend=2*neff-1;
for(IdxTy j=1; j<=neff; ++j) x=x*D(j);
for(IdxTy j=1; j<=ii; ++j) x=x/D(j);

}
m_ratio_cache.push_back(x);
//MM_MSG(MMPR3(m_der.size(),p.to_string(),x))
}

}
void set_example1() { }

void set_unity() {  }
void set_1n() { }
void dump()
{
const IdxTy sz=m_der.size();
for(IdxTy i=0; i<sz; ++i)
{
MM_MSG(MMPR(i)<<MMPR(m_der[i].to_string()))

} // i 


}

Po m_poly;
Dstack m_der;
Vstack m_ratio_cache;
Mi mi;
bool m_one;
IdxTy m_offset;
}; // 
/////////////////////////////////////////////////////////
class func_erfphi_class
{
//typedef mjm_fixed_sum_itor GenItor;
typedef func_erfphi_class Myt;
typedef mjm_integrals Mi;
typedef unsigned int IdxTy;
typedef double D ;
typedef shape_function_integrals SFi;
typedef SFi::simple_polynomial<D> Po;
typedef std::vector<Po> Dstack;
public: 
func_erfphi_class(): m_one(false),m_offset(0) {}
const IdxTy & offset() const { return m_offset;}
const IdxTy & offset(const IdxTy off)  { m_offset=off; return m_offset;}

D operator()(const IdxTy n, const D z) //const
{
//if (m_one) return 1.0; 
//static const D fa=sqrt(M_PI)*.5;
//if (n==0) return fa*erf(z);
make_der(n+m_offset);
const Po &  p=m_der[n+m_offset];
const D x=mi.evaluate_polynomial(p,z)*exp(-z*z);
//MM_MSG(MMPR4(x,n,z,mi.print_polynomial(p)))
return x; 
}
void make_der(const IdxTy n)
{
if (m_der.size()==0){ Po crap; crap.push_back(1);  m_der.push_back(crap); } 
while (m_der.size()<=n) 
{
Po p;
//mi.differentiate_polynomial(p,m_der[m_der.size()-1]);
mi.next_gauss_der_polynomial(p, m_der[m_der.size()-1]);
m_der.push_back(p);
}

}
void set_example1() { }

void set_unity() {  }
void set_1n() { }


Po m_poly;
Dstack m_der;
Mi mi;
bool m_one;
IdxTy m_offset;
}; // 
/////////////////////////////////////////////////////////

void compare_erf(std::vector<double> & coef)
{
mjm_integrals mi;
std::vector<double>  deste,e1,destes;
e1.push_back(1);
e1.push_back(-1);
mi.composite_polynomial(deste,coef,e1);
mi.multiply_polynomial(destes,deste,2.0/sqrt(M_PI));

double x=0;
unsigned int sz=coef.size();
for (int i=0; i<100; ++i)
{ 
	double expx=(exp(-x*x));
	double theta=sqrt(1-expx);
	double v=theta;
	double ve=theta; // 1; // expx;
	double sum=0;
	double sume=0;
	for ( unsigned int j=0; j<sz; ++j)
	{
		sum+=v*coef[j];
		v=v*theta;
	} // j 
	for ( unsigned int j=0; j<destes.size(); ++j)
	{
		sume+=ve*destes[j];
		ve=ve*expx;
	}

	sum=sum*2.0/sqrt(M_PI);
	double exact=erf(x);
	double del=exact-sum;
	double dele=exact-sume;
	double delee=sum-sume;
	MM_MSG(MMPR4(x,exact,sum, del)<<MMPR3(sume,dele,delee))
	x=x+.1;
}

}




//////////////////////////////////////////////////////////
int main(int argc,char **args)
{

typedef mjm_integrals Myt;
typedef mjm_burmann Myb;
Myt mi;
Myb mib;
typedef unsigned int IdxTy;
typedef double D ;
typedef func_class Func;
typedef func_erf_class FuncErf;
typedef func_erfphi_class FuncErfPhi;
typedef func_erfphistar_class FuncErfPhiStar;
const bool main1=false;
const bool main2=!false;


//const IdxTy nu=1;
const IdxTy flags=0;
typedef std::vector<D> Dest;
Dest dest,desthc,destalg;
Func func;
FuncErf funcerf;
FuncErfPhi funcerfphi;
FuncErfPhiStar funcerfphistar;
if (main1)
{MM_MSG(" args are v,mu,nu,kmax");
int  pos=1;
const D v=(argc>pos)?atof(args[pos]):0; ++pos;
const IdxTy  mu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  nu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  kmax=(argc>pos)?atoi(args[pos]):0; ++pos;
const D  z0=(argc>pos)?atof(args[pos]):0; ++pos;



func.set_example1();
//func.set_unity();
//Burman _eq21
const D m=D(mu)/D(nu);
const D fzed=pow(func(0,z0),m);
MM_MSG(MMPR4(mu,nu,kmax,flags)<<MMPR2(v,z0))
//Myt::burmann_ref1_eq21(dest,  func,z0,  mu,  nu,  kmax,  flags);
mi.burmann_ref1_eq21(dest,  func,z0,  mu,  nu,  kmax,  flags);
for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG(MMPR2(i,dest[i])) }
const D dz=.2;
for (D z=-1; z<1; z+=dz)
{
const D f=pow(func(0,z),m);
const D fberm=fzed/mi.evaluate_polynomial(dest,z-z0);
MM_MSG(MMPR4(z,(z-z0),f,fberm)<<MMPR((f-fberm)))
}
return 0;
} // main1

if (main2)
{MM_MSG(" erf bermann args are nmax, z0,nu,offset");
int  pos=1;
//const D v=(argc>pos)?atof(args[pos]):0; ++pos;
//const IdxTy  mu=(argc>pos)?atoi(args[pos]):0; ++pos;
//const IdxTy  nu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  nmax=(argc>pos)?atoi(args[pos]):0; ++pos;
const D  z0=(argc>pos)?atof(args[pos]):0; ++pos;
const IdxTy nu =(argc>pos)?atoi(args[pos]):1; ++pos;
const IdxTy off =(argc>pos)?atoi(args[pos]):0; ++pos;
if (nu==0)
{
funcerfphistar(10,0);
funcerfphistar.dump();
return 0;
}
MM_MSG(MMPR4( z0, nu,  nmax, flags)<<MMPR(off));
funcerfphi.offset(off);
//Myt::burmann_ref1_eq14(dest, funcerf, funcerfphi, funcerfphistar, z0, nu,  nmax, flags);
//mib.burmann_ref1_eq14(dest, funcerf, funcerfphi, funcerfphistar, z0, nu,  nmax, flags);
// now too slow to use 
if (false) mib.burmann_ref1_eq14(dest, funcerf, funcerfphi, funcerfphistar, z0, nu,  nmax, flags);
//mib.burmann_ref1_eq14_erf( desthc, nmax);
auto erff=Myb::erf_demo();
auto erffrat=Myb::erf_demo_rat();
typedef std::vector<Myb::erf_demo_rat::calc_type> RatVec;
RatVec destrat;
auto erfalg=Myb::erf_demo_alg(.001);
 //mib.burmann_1( desthc,  erff,  nmax,flags);
 mib.burmann_1( dest,  erff,  nmax,flags);
 mib.burmann_1( destrat,  erffrat,  nmax,flags);
 mib.burmann_1( destalg,  erfalg,  nmax,flags);
compare_erf(dest);
while (desthc.size()<dest.size()) desthc.push_back(0); 
while (destalg.size()<dest.size()) destalg.push_back(0); 
while (dest.size()<destrat.size()) dest.push_back(0); 
for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG(" erf "<<MMPR4(i,dest[i],destrat[i],destalg[i])) }
//for (IdxTy i=0; i<destrat.size(); ++i) { MM_MSG("bermann erf "<<MMPR2(i,destrat[i])) }
MM_MSG(MMPR(mi.print_polynomial(dest)))
Dest deste,e1,destes;
e1.push_back(1);
e1.push_back(-1);
mi.composite_polynomial(deste,dest,e1);
mi.multiply_polynomial(destes,deste,2.0/sqrt(M_PI));
MM_MSG(MMPR(mi.print_polynomial(deste)))
MM_MSG(MMPR(mi.print_polynomial(destes)))
mi.dump_counter("final");
mib.dump_counter("final");
return 0;
} // main2
} // main 
#endif //   TEST_BURMANN2__ 
#ifdef TEST_BURMANN3__
int main(int argc,char **args)
{

typedef mjm_integrals Myt;
typedef mjm_burmann Myb;
Myt mi;
Myb mib;
typedef unsigned int IdxTy;
typedef double D ;

const IdxTy flags=0;
typedef std::vector<D> Dest;
{MM_MSG(" args are nmax,z0");
int  pos=1;
//const D v=(argc>pos)?atof(args[pos]):0; ++pos;
//const IdxTy  mu=(argc>pos)?atoi(args[pos]):0; ++pos;
//const IdxTy  nu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  nmax=(argc>pos)?atoi(args[pos]):0; ++pos;
const D  z0=(argc>pos)?atof(args[pos]):0; ++pos;
const D  z0n=(argc>pos)?atof(args[pos]):0; ++pos;
const D  z0d=(argc>pos)?atof(args[pos]):0; ++pos;
MM_MSG(MMPR4( nmax,z0,z0n,z0d));

auto erff=Myb::erf_demo();
auto erffrat=Myb::erf_demo_rat();
auto erffalgrat=(z0d==0)?Myb::erf_demo_alg_rat(z0)
:Myb::erf_demo_alg_rat(z0n,z0d);
Dest dest,desthc,destalg,dest2alg;
typedef std::vector<Myb::erf_demo_rat::calc_type> RatVec;
RatVec destrat,destalgrat;
auto erfalg=Myb::erf_demo_alg(z0);
auto erf2alg=Myb::erf_demo2_alg(z0);
 //mib.burmann_1( desthc,  erff,  nmax,flags);
 mib.burmann_1( dest,  erff,  nmax,flags);
 mib.burmann_1( destrat,  erffrat,  nmax,flags);
// mib.burmann_1( destalgrat,  erffalgrat,  nmax,flags);
 mib.burmann_1( destalg,  erfalg,  nmax,flags);
 mib.burmann_2( dest2alg,  erf2alg,  nmax,flags);
//compare_erf(dest);
while (dest.size()<destrat.size()) dest.push_back(0); 
while (dest.size()<destalgrat.size()) dest.push_back(0); 
while (dest.size()>destalgrat.size()) destalgrat.push_back(0); 
while (desthc.size()<dest.size()) desthc.push_back(0); 
while (destalg.size()<dest.size()) destalg.push_back(0); 
while (dest2alg.size()<dest.size()) dest2alg.push_back(0); 

typedef std::stringstream Ss;
for (IdxTy i=0; i<dest.size(); ++i) 
{ 
Ss ss;
ss<<" erf "<<MMPR2(i,z0) ;
ss<<MMPR(dest[i]);
ss<<MMPR(destrat[i]);
ss<<MMPR(destalg[i]);
ss<<MMPR(destalgrat[i]);
ss<<MMPR(double(destalgrat[i])) ;
ss<<MMPR(dest2alg[i]);
MM_MSG(ss.str())
}

//for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG(" erf "<<MMPR2(i,z0)<<MMPR3(dest[i],destrat[i],destalg[i])<<MMPR(dest2alg[i])) }
//for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG(" erf "<<MMPR2(i,z0)<<MMPR3(dest[i],destrat[i],destalg[i])) }
//for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG(" erf "<<MMPR2(i,z0)<<MMPR4(dest[i],destrat[i],destalg[i],destalgrat[i])<<MMPR(double(destalgrat[i]))) }


//for (IdxTy i=0; i<destrat.size(); ++i) { MM_MSG("bermann erf "<<MMPR2(i,destrat[i])) }
if (false)
{
MM_MSG(MMPR(mi.print_polynomial(dest)))
Dest deste,e1,destes;
e1.push_back(1);
e1.push_back(-1);
mi.composite_polynomial(deste,dest,e1);
mi.multiply_polynomial(destes,deste,2.0/sqrt(M_PI));
MM_MSG(MMPR(mi.print_polynomial(deste)))
MM_MSG(MMPR(mi.print_polynomial(destes)))
}
//mi.dump_counter("final");
mib.dump_counter("final");
return 0;
} // main2
} // main 
#endif //   TEST_BURMANN3__ 


#ifdef TEST_BURMANN4__

int main(int argc,char **args)
{

//typedef mjm_integrals Myt;
typedef mjm_burmann Myb;
//Myt mi;
Myb mib;
typedef unsigned int IdxTy;
typedef double D ;

const IdxTy flags=0;
//typedef std::vector<D> Dest;
MM_MSG(" args are nmax,z0n,z0d");
int  pos=1;
const IdxTy  nmax=(argc>pos)?atoi(args[pos]):0; ++pos;
const D  z0n=(argc>pos)?atof(args[pos]):0; ++pos;
const D  z0d=(argc>pos)?atof(args[pos]):0; ++pos;

typedef mjm_rational RatTy;
//typedef RatTy::IntTy IntTy;
typedef RatTy Fty;
//typedef D Fty;
typedef mjm_closed_families::exp_poly_ratio<Fty> Nr;
//typedef std::vector<Fty> Cache;
//typedef mjm_integrals Mi;
//typedef shape_function_integrals SFi;
//typedef SFi::simple_polynomial<Fty> Po;
typedef std::vector<Nr> Dstack;
Fty z0= Fty(int(z0n),int(z0d));
z0.set_reduce();
double z2=double(z0);
Fty ze= Fty(int(0),int(1));
ze.set_reduce();
MM_MSG(MMPR4(nmax,z0n,z0d,z0))
if (false) { 
Dstack dest=Dstack();
Nr minusone=Nr();
mib.erf_fd_poly(dest,  nmax);
for (IdxTy i=0; i<dest.size(); ++i)
{
Nr & p=dest[i];
Fty v=p.evaluate(z0,ze);
MM_MSG(MMPR4(i,p.to_string(),v,double(v))<<MMPR((erf(z2*i)*0.8862269))) 

}
} // alt
// void erf_fd_rat(Ty & dest, const IdxTy n, const Trat & delta)
typedef std::vector<Fty> Dt;
Dt dest=Dt();
mib.erf_fd_rat(dest, nmax, z0);
return 0;
}

#endif //   TEST_BURMANN4__ 

#endif // guard 

