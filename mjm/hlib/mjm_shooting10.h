#ifndef MJM_SHOOTING5_H__
#define MJM_SHOOTING5_H__

#include "mjm_globals.h"
#include "mjm_shooting_logic.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_diffuse_iterators.h"
#include "mjm_defect_levels.h"

#include <algorithm>
#include <signal.h>
#define USING_PETSC
#ifdef USING_PETSC
// deviate from shooting to matrix 
#include "mjm_petsc_util.h"
#include "mjm_petsc_fd_base.h"
#endif
/*

FIXME kluges in drift bc code 
TODO the trap gr is not FEM ized and neither is velocity

TODO verify the gr integral makes sense

TODO check the drift code for very small and zero velocity. 
Overlap diffusion is a big issue for real accuracy. 

TODO positive terminal still has blip on right

On quick impulse test with low vt this looks right - the
diffusion is messing up boundary conditions. 
The 100 point version is not working but it looks like the 1k point
boundaries are almost good except for a blip on the positive right
after the last run. 


*/

class debug_levels
{
typedef debug_levels Myt;
enum {BAD=~0};
public:
debug_levels() : m_level(0), m_enable(false),m_mom(0) {}
debug_levels(const bool x ) : m_level(0), m_enable(x),m_mom(0) {}
// this works I guess 
debug_levels( const Myt & x): m_level(x.m_level-1), m_enable(x.m_enable), m_mom(&x) {}
operator bool() const { return m_enable; }
void gdb() const { raise(SIGINT); } 
// this needs a way to keep track of levels and only enable
// soo deep 
int m_level;
bool m_enable;
const Myt * m_mom;
}; //debug_levels


/*
mjm_shooting10.h : trying to understand FEM problems with drift diffusion. 
Starting with FD, try to develop Ficks Law from microscopic minimal system.
See what happens, relate to macroscopic properies and see if results work.
The microscopic Fick derivation leads to Lancozs integral derivaties
so the approach is called Fick-Lancozs. The local physics should lead
to better local equations and the integration naturally leads to FEM- like
things. 

Determine if retarded time is needed for stability.  Finite propogation velocity
seems to be an important thing but retarded effects may be too. 

The DD operators do not appear to attribute the results of motion calculations
leading to bumps at the boundaries and probably more subtle problems. 
Introduce shape functions and attribute the pices to the nodes. 

*/

/*

backup from liquid to pnj for quick validation
fudd 
p(Na)  | n(Nd) 

g++ -DTEST_SHOOTING__ -Wall  -std=gnu++11 -gdwarf-3 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_shooting10.h

g++ -DTEST_SHOOTING__ -Wall  -s=gnu++11 -gdwarf-3 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -Wno-unused-function -Wno-non-template-friend -x c++ mjm_shooting10.h

*/
#ifdef USING_PETSC

////////////////////////////////////////////////////////
#endif // USING_PETSC


///////////////////////////////////////////////////////////////
class quadratic_fick_lanc_element
{
/*
 This is a consistent means of interpolating between FD points and allowing
integrals over surrounding areas for diffusion calculations. Note that
"drift" is often used interchangably with diffusion lol. This is really two
kinds of propogators from t0 to t1 but just moving stuff around. 

There is a recurring problem with density and total carriers. The values are
densities but integrals measure total count into and out of a region. Getting
1/h right has not been checked. 

*/


class Tr
{
public:

typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
typedef std::stringstream Ss;


}; // Tr

// this needs a StrTy def somewhere? Where wtf? 
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;

public:


template <class Ti> D isselflin(const Ti & di) const
{ return islin(di.x0(),di.x2(),di.x0(),di.x1(),di.x2(), di.n0(),di.n1(),di.n2()); }


// quadratic fit  returns integrated amount in [a,b] that are clipped into [x0,x2].
D is(const D & _a, const D & _b, const D & x0, const D & x1, const D & x2, const D & n0, const D & n1, const D & n2) const
{
D res=0;
D a=_a; if (a<x0) a=x0; if (a>x2) a=x2;
D b=_b; if (b>x2) b=x2; if (b<x0) b=x0;
// the thing has clipped oor 
if (a>=b) return res;
// fit to Ax^2 + Bx + C = N
// C=N0, A(x1-x0)^2+ B(x1-x0)+C=N1
// C=N0, A(x2-x0)^2+ B(x2-x0)+C=N2
//  A(x2-x0)^2=- B(x2-x0)-C+N2, 
//  A=(-B(x2-x0)-C+N2)/(x2-x0)^2, 
//  B(x1-x0)=-A(x1-x0)^2-C+N1
//  B(x1-x0)=(B(x2-x0)+C-N2 )(x1-x0)^2/(x2-x0)^2 -C+N1
//  B(x1-x0-F)=(+C-N2 )(x1-x0)^2/(x2-x0)^2 -C+N1, F=(x2-x0)^(-1)(x1-x0)^2;
const D C=n0;
const D x10=x1-x0;
const D x20=x2-x0;
const D  F=(x1-x0)*(x1-x0)/(x2-x0);
const D B=(-C+n1+(C-n2)*F/(x2-x0))/(x1-x0-F);
const D A=(-C+n1-B*(x1-x0))/(x1-x0)/(x1-x0);
const D b0=b-x0; 
const D a0=a-x0;
res+=1.0/3.0*A*b0*b0*b0+.5*B*b0*b0+C*b0;
res-=1.0/3.0*A*a0*a0*a0+.5*B*a0*a0+C*a0;

return res;
}



/////////////////////////////////////////////
// linear fit and integrate
D islin(const D & _a, const D & _b, const D & x0, const D & x1, const D & x2, const D & n0, const D & n1, const D & n2) const
{
D res=0;
D a=_a; if (a<x0) a=x0; if (a>x2) a=x2;
D b=_b; if (b>x2) b=x2; if (b<x0) b=x0;
// the thing has clipped oor 
if (a>=b) return res;
D A,B,C;
// left side first 
fitlin(A, B, C, x0, x1, x2, n0, n1, n2, false); // const
const D xleft=x1-x0;
res=(C+.5*B*xleft)*xleft;
fitlin(A, B, C, x0, x1, x2, n0, n1, n2, !false); // const
const D xright1=x1-x0;
const D xright2=x2-x0;
const D xmid=(xright2+xright1)*.5;
// these are fit from x0,
res+=(C+B*xmid)*(xright2-xright1);
//MM_MSG(" fit "<<res<<" n0="<<n0<<" n1="<<n1<<" n2="<<n2)

return res;
}




////////////////////////////////////////////////


// the quadratic fit needs to be to ln(n) NOT n but this is too confusing for now.
// so doe a linear fit even though it adds to pwl confusion the boundaries coincide with sf bounds

template <class Di> 
void fitlin(D & A,  D & B, D & C,const Di & di, const bool right) const
{
const D x0=di.x0();
const D x1=di.x1();
const D x2=di.x2();
const D n0=di.n0();
const D n1=di.n1();
const D n2=di.n2();
fitlin(A,B,C,x0,x1,x2,n0,n1,n2,right);

}
void fitlin(D & A,  D & B, D & C,const D & x0, const D & x1, const D & x2, const D & n0, const D & n1, const D & n2, const bool right) const
{
// right means use x2 and x1 for coefs bu relative to x0,
// n=b*(x-x0) + c 
A=0;
if (right)
{
// n1=b*(x1-x0) +c
// n2=b*(x2-x0) +c
// n2-n1=b*(x2-x1) 
B=(n2-n1)/(x2-x1);
C=n1-B*(x1-x0);
return; 
}
B=(n1-n0)/(x1-x0);
C=n0;

}
void fitlinbase(D & A,  D & B, D & C,const D & x0, const D & x1, const D & x2, const D & n0, const D & n1, const D & n2, const bool right,const D & xbase) const
{
// right means use x2 and x1 for coefs bu relative to x0,
A=0;
if (right)
{
B=(n2-n1)/(x2-x1);
// this normally makes C into n2 
C=n1-B*(x1-xbase);
return; 
}
B=(n1-n0)/(x1-x0);
C=n1-B*(x1-xbase); // for xbase x0 should be same 
//C=n0;

}



template <class Di> 
void fitq(D & A,  D & B, D & C,const Di & di) const
{
C=di.n0();
const D x0=di.x0();
const D x1=di.x1();
const D x2=di.x2();
const D n0=di.n0();
const D n1=di.n1();
const D n2=di.n2();

const D x10=x1-x0;
const D x20=x2-x0;
const D F=(x1-x0)*(x1-x0)/(x2-x0);
B=(-C+n1+(C-n2)*F/(x2-x0))/(x1-x0-F);
A=(-C+n1-B*(x1-x0))/(x1-x0)/(x1-x0);
// at least check
const D n2check=A*(x2-x0)*(x2-x0) + B*(x2-x0) +C;
if (n2check!=n2) { MM_ONCE(" some diff "<<n2check<<" vs "<<n2<<" "<<(n2-n2check),) } 

}

void fitq(D & A,  D & B, D & C,const D & x0, const D & x1, const D & x2, const D & n0, const D & n1, const D & n2) const
{
C=n0;
const D x10=x1-x0;
const D x20=x2-x0;
const D F=(x1-x0)*(x1-x0)/(x2-x0);
B=(-C+n1+(C-n2)*F/(x2-x0))/(x1-x0-F);
A=(-C+n1-B*(x1-x0))/(x1-x0)/(x1-x0);
}

/*
// integrate the g(v) function times a quadratic fit to the carrier density
// note this is not in retarded time yet.n(x)g(v)dx 
*/



/*
TODO move this comment block somewhere wtf 
Calculate the number of carriers attributable to x1 that move to xp1
due to drift amount vt.  Note that linear shape functions are included
in both source and destination so this is not just integrated
number of carriers that move during time t.  See notes.

First, vt>0 implies xp>x
Quadratic fit to carrier concentration is currently used with
linear shape functions centered on x1. Work relative to x0 and
define z=xp-vt-x0 and zp=xp-xp0 and then k=vt+x0-xp0 and zp=z+k
Also note that xp =vt+x or z=x-x0
There are 2x2 breakpoints around x1 and xp1 or
z<x1-x0 and zp<xp1-xp0 or z<xp1-xp0-k 
Integration limits are on the source, between x0 and x2 but
clipped by desination lmits xp0 and xp2 


*/


// integral over source and dest as source ntegral as variable endpoints.
// cxmax could be xmax times some coef such as -1, cabs may be redundant
// but is what multuiplies x-xpZZ
// xp shape function is s1*xp + s0 and x is evaluated at a*xp+b
/*
This is the working core code for the FL diffusion with linear verlocity
profile. It is a double integral over source x and destination x' or xp
evaluated at a single point xpj with the inner x integral limit of
a*x'+b where  is either 1 or zero and b contains an origin shift term
and some breakpoint such as xmax. 
The source shape function is x1*x+c0 and the desination is s1*x' + s0.
The n(x) is ax^2+Bx+C , caller needs to make sure the origins are right

*/

D piece_sf_var(const StrTy &lbl, const D & xpj, const D & a, const D & b, const D & c0, const D & c1,const D & k, const D & cx, const D & cxp
, const D & s1, const D & s0
, const D & A, const D & B, const D & C, const bool dbg=true) const
{  
if (dbg) {MM_MSG(lbl<<" xpj="<<xpj<<" a="<<a<<" b="<<b<<" c0="<<c0<<" c1="<<c1<<" k="<<k<<" cx="<<cx<<" cxp="<<cxp<<" s1="<<s1<<" s0="<<s0<<" A="<<A<<" B="<<B<<" C="<<C) }
if ((a==0)&&(b==0)) return 0; 
const D xpjab=(a*xpj+b);
if ((xpjab==0)) return 0; 
D res=0;
const D P4=A*c1*cx;
const D P3=B*c1*cx+A*c1*k+A*c0*cx;
const D Q3=A*c1*cxp;
const D P2=A*c0*k+B*c1*k+B*c0*cx+C*c1*cx;
const D Q2=A*c0*cxp+c1*B*cxp;
const D P1=B*c0*k + C*c0*cx + C*c1*k ;
const D Q1=B*c0*cxp+C*c1*cxp;
const D P0=C*c0*k ;
const D Q0=cxp*C*c0;
const D c13=1.0/3.0; const D c16=1.0/6.0; const D c17=1.0/7.0; 
const D c15=1.0/5.0; const D c14=1.0/4.0;  const D c12=1.0/2.0; 
if (dbg)
{
MM_MSG(" P4="<<P4<<" P3="<<P3<<" Q3="<<Q3<<" P2="<<P2<<" Q2="<<Q2<<" P1="<<P1<<" Q1="<<Q1<<" P0="<<P0<<" !0="<<Q0 )
}
// now evaluate at x =xa*ap+b and work in variable "x" if a!=0,
// the limit of the x integration is a constant so this is a simple polynomial in xp
if (a==0)
{
const D R=.2*P4*b*b*b*b*b + .25*P3*b*b*b*b+ 1.0/3.0*P2*b*b*b + .5*P1*b*b+P0*b;
const D S=.25*b*b*b*b*Q3+ 1.0/3.0*Q2*b*b*b + .5*Q1*b*b+ Q0*b;
// the shape function multiples the linear x' puting a cubic term in result
// this looks right 
const D xpp2=S*s1;
const D xpp1=R*s1+S*s0;
const D xpp0=R*s0;
const D t3nx=c13*xpp2;
const D t2nx=.5*xpp1;
const D t1nx= xpp0;
const D sumnx=((t3nx*xpj+t2nx)*xpj+t1nx)*xpj; //  c13*xpp2*xpj*xpj*xpj+.5*xpp1*xpj*xpj+ xpp0*xpj;
if (dbg) {

const D t3=c13*xpp2*xpj*xpj*xpj;
const D t2=.5*xpp1*xpj*xpj;
const D t1= xpp0*xpj;
const D sum=t3+t2+t1; //  c13*xpp2*xpj*xpj*xpj+.5*xpp1*xpj*xpj+ xpp0*xpj;
MM_MSG("a==0 piece "<<lbl<<" t3="<<t3<<" t2="<<t2<<" t1="<<t1<<" sum="<<sum) 
MM_MSG("a==0 piece "<<lbl<<" t3nx="<<t3nx<<" t2nx="<<t2<<" t1nx="<<t1<<" delsum="<<(sum-sumnx)) }
return sumnx; // c13*xpp2*xpj*xpj*xpj+.5*xpp1*xpj*xpj+ xpp0*xpj;
}
// if a!=0, change variables to x=axp+b and then dx=a dxp
// so, the s1*xp+ s0 shape function becomes s1*x/a-s1*b/a+s0 
//const D xpjab=(xpj-b)/a; // (a*xpj+b);
//const D xpjab=(a*xpj+b);
// this is the second inteagrl times the x' from first integral AND shape function
// note that theere is a 1/a AND 1/xmax^2 being ignored 
// P*cs1* x 
const D cs1=s1/a; const D cs0=s0-b*s1/a;

const D cq0=-b*(cs0)/a;
const D cq1=(cs0/a-b/a*cs1);
const D cq2=cs1/a;

// this is the alt impl which should be similar to the original one 
const D F7=cs1*c15*c17*P4               +cq2*c14*c17*Q3;
const D F6=cs1*c14*c16*P3+cs0*c15*c16*P4+cq2*c13*c16*Q2+cq1*c14*c16*Q3;
const D F5=cs1*c13*c15*P2+cs0*c14*c15*P3+cq2*c12*c15*Q1+cq1*c13*c15*Q2+cq0*c14*c15*Q3;
const D F4=cs1*c12*c14*P1+cs0*c13*c14*P2+cq2*c14*Q0    +cq1*c12*c14*Q1+cq0*c13*c14*Q2;
const D F3=cs1*c13*P0    +cs0*c12*c13*P1               +cq1*c13*Q0    +cq0*c12*c13*Q1;
const D F2=               cs0*c12*P0                                  +cq0*c12*Q0;
const D w=xpjab;
const D sumnx=(((((F7*w+F6)*w+F5)*w+F4)*w+F3)*w+F2)*w*w/a;

if (dbg) {
const D xpjab2=xpjab*xpjab;
const D xpjab3=xpjab*xpjab2;
const D xpjab4=xpjab2*xpjab2;
const D xpjab5=xpjab4*xpjab;
const D xpjab6=xpjab3*xpjab3;
const D xpjab7=xpjab3*xpjab4;

const D I1= cs1*(c15*c17*P4*xpjab7+c14*c16*P3*xpjab6+c13*c15*P2*xpjab5+c12*c14*P1*xpjab4+c13*P0*xpjab3);
const D I2= cs0*(c15*c16*P4*xpjab6+c14*c15*P3*xpjab5+c13*c14*P2*xpjab4+c12*c13*P1*xpjab3+c12*P0*xpjab2);
//const D cq0=-b*(cs0)/a;
//const D cq1=(cs0/a-b/a*cs1);
//const D cq2=cs1/a;
const D I3=cq2*(c14*c17*Q3*xpjab7+c13*c16*Q2*xpjab6+c12*c15*Q1*xpjab5+c14*Q0*xpjab4);
const D I4=cq1*(c14*c16*Q3*xpjab6+c13*c15*Q2*xpjab5+c12*c14*Q1*xpjab4+c13*Q0*xpjab3);
const D I5=cq0*(c14*c15*Q3*xpjab5+c13*c14*Q2*xpjab4+c12*c13*Q1*xpjab3+c12*Q0*xpjab2);
const D I6=0;
// the factor of x,ax^2 is left off yet 
res=(I1+I2+I3+I4+I5+I6)/a;

MM_MSG("a!=0 piece "<<lbl<<" I1="<<I1 <<" I2="<<I2 <<" I3="<<I3 <<" I4="<<I4
<<" I5="<<I5 <<" I6="<<I6<<" res="<<res<<" vs "<<sumnx<<" diff "<<(res-sumnx)) 
}
return sumnx; // res;
}
/*
template <class Td>
D diff_sf(const Td & di, const D & xmax, const D & vt, const bool dbg=true )
{
return diff_sf(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,dbg);
}
*/

template <class Td>
D diff_sf_loop(const Td & di, const D & xmax, const D & vt, const bool no_self, const bool dbg=true )
{
MM_ONCE( " use slower but more accurate sf_bc instead of loop for now  ",)
//return diff_sf_loop(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,no_self,dbg);
return diff_sf_bc(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,no_self,dbg);
}
// compute self to self integral regardless of what des the iterator points to 
template <class Td>
D diff_sf_self_loop(const Td & di, const D & xmax, const D & vt, const bool dbg=true )
{
MM_ONCE( " self use slower but more accurate sf_bc instead of loop for now  ",)
//return diff_sf_loop(di.x0(),di.x1(),di.x2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,false,dbg);
return diff_sf_bc(di.x0(),di.x1(),di.x2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,false,dbg);
}


// find sf weighted self amount of charge
// this is not right either, as it is self charge attirutable to self
// just want one sf.
template <class Td>
D self_sf_loop(const Td & di , const bool dbg=true )
{
//return self_sf_loop(di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),dbg);
drift_flags df(false,true,dbg);
//return drift_sf_bc(di.x0(),di.x1(),di.x2(),di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),0,false,dbg);
return drift_sf_bc(di.x0(),di.x1(),di.x2(),di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),0,df);
}



void make_sf(D & s1, D & s0, const D & x0, const D & x1, const D & x2, const IdxTy i) const 
{
// these are not consistent.  not suew about s0
if (i==1) { s0=(x2-1.0*x0)/(x2-x1);  s1=-1.0/(x2-x1);  } // xp>xp1 right
else if( i==0) {  s1=1.0/(x1-x0);  s0=0;   } // xp<xp1 left 
else {s0=0; s1=0; }
//MM_MSG(" make_sf i="<<i<<" x0="<<x0<<" "<<x1<<" "<<x2<<" s1="<<s1<<" s0="<<s0)
return; 
}
void make_sf(D & s1, D & s0, const D & x0, const D & x1, const D & x2, const IdxTy i, const D & xbase) const 
{
// these are not consistent.  not suew about s0
//if (i==1) { s0=(x2-1.0*xbase)/(x2-x0);  s1=-1.0/(x2-x1);  } // xp>xp1 right
if (i==1) { s0=(x2-1.0*xbase)/(x2-x1);  s1=-1.0/(x2-x1);  } // xp>xp1 right
else if( i==0) {  s1=1.0/(x1-x0);  s0=(xbase-x0)/(x1-x0);   } // xp<xp1 left 
else {s0=0; s1=0; }
//MM_MSG(" make_sf i="<<i<<" x0="<<x0<<" "<<x1<<" "<<x2<<" s1="<<s1<<" s0="<<s0)
return; 
}




void make_cx(D & cx, D & cxp, const D & x, const D & xp) const
{
if (x<xp) cx=1; else cx=-1;
cxp=-cx;
}
// the limits is consant, a=0, if the source completely fillx the range
// this is the limit on x expressed as x=a*xp+b
// since this is the UPPER limit, does that effect the equal condition?
void make_ul(D & a,D & b,const D & x2,const D &xalt, const D & xm) const
{
if (x2<xalt) {a=0; b=x2; }
else { a=1; b=xm; }
}
void make_ll(D & a,D & b,const D & x0,const D &xalt, const D & xm) const 
{
if (x0>=xalt) {a=0; b=x0; }
else { a=1; b=-xm; }
}

// segment the doble integral over source and destination weighted by
// "shape functions" to attribute density to middle node x1 and xp1.
// Hopefully this conserves things better... 
// x is the source element and integration variable, xp or x prime x' 
// is the desination area
// may be useful for testing 
#if 0 

// this should be the distance of closest approach 
// for disjoint, either |    |       |     |
//                     xp0  xp2      x0    x2
// so this is x0-xp2 and the source disto is non-zero from x0-xmax to x2+xmax 
// but the source integral changes if x2-xmax>xp0
// for disjoint, either |    |       |     |
//                     x0    x2     xp0   xp2
// xp0-x2 and the source distro is the same with int limit hitting xp2>x0+xmax

#endif

/*
D self_sf_loop(const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2,   const bool dbg=false) const
{
D res=0;
D A,B,C,c1,c0;
const bool lin_fit=true;
if (!lin_fit) fitq(A,B,C,x0,x1,x2,n0,n1,n2);
//if (lin_fit) fitlin(A,B,C,x0,x1,x2,n0,n1,n2,src_right);
bool src_right=true;
bool src_left=!true;
make_sf(c1,c0, x0, x1, x2,(src_left)?1:0);
if (lin_fit) fitlin(A,B,C,x0,x1,x2,n0,n1,n2,src_left);
// this is simply the integral of shape function times the charge polynomial
// (c1x+c0)(Ax2+Bx+c)
D v=x0-x0;
D I1= .25*c1*A*v*v*v*v + 1.0/3.0*(c1*B+c0*A)*v*v*v+.5*(B*c0+c1*C)*v*v+C*c0*v;
v=x1-x0;
D I2= .25*c1*A*v*v*v*v + 1.0/3.0*(c1*B+c0*A)*v*v*v+.5*(B*c0+c1*C)*v*v+C*c0*v;
// the right side should be ok the sf is based on x-x0
if (lin_fit) fitlin(A,B,C,x0,x1,x2,n0,n1,n2,src_right);
make_sf(c1,c0, x0, x1, x2,(src_right)?1:0);
v=x1-x0;
D I3= .25*c1*A*v*v*v*v + 1.0/3.0*(c1*B+c0*A)*v*v*v+.5*(B*c0+c1*C)*v*v+C*c0*v;
v=x2-x0;
D I4= .25*c1*A*v*v*v*v + 1.0/3.0*(c1*B+c0*A)*v*v*v+.5*(B*c0+c1*C)*v*v+C*c0*v;
res=I2-I1+I4-I3;

// right now the 2x is a kluge fudge factor but there is probably a reason for it...
res=2.0*res/(x2-x0); 
if (dbg) { MM_MSG(" self_sr_loop res="<<res<<" x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" A="<<A<<" B="<<B<<" C="<<C) } 

return res;

} //self_sf_loop
*/

/*
Evaluate the double x and x' integral at the 4 corners and sum
into res, return the segment reslt for debuggin
*/

template <class Ts> 
D corners(const Ts & lbl, D & res, const D & cs,const D & ulp, const D & llp, const D & au, const D & bu, const D &al, const D & bl, const D & c0, const D & c1, const D & k, const D & cx, const D & cxp, const D & s1, const D & s0, const D & A, const D & B, const D & C, const bool dbg) const
{
 D p1=piece_sf_var("p1", ulp, au, bu, c0, c1, k, cx, cxp, s1, s0, A, B, C,dbg);  
 D p2=piece_sf_var("p2", llp, au, bu, c0, c1, k, cx, cxp, s1, s0, A, B, C,dbg);  
 D p3=piece_sf_var("p3", ulp, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C,dbg);  
 D p4=piece_sf_var("p4", llp, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C,dbg);  
const D segr=cs*(p1-p2-p3+p4);
res+=segr; 
if (dbg) {MM_MSG(lbl<<" llp="<<llp<<" ulp="<<ulp<<" p1="<<p1<<" p2="<<p2<<" p3="<<p3<<" p4="<<p4<<" segr="<<segr<<" res="<<res) } 
if (dbg) if (segr<(-1e-3)) raise(SIGINT);
return segr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// debug diff_sf_loop for the case of boundary linked to infinte pixel beyond end.

//rangeck("Case2s1b",ulp, ulpx, au, bu, al, bl, k,xmax, xmins,xmaxs,xpmin,xpmax dbg);

/*
Verify that the ranges for the segments seem to make sense

*/
void rangeck(const StrTy & lbl,const D & ulp, const D & llp
	, const D & au, const D & bu, const D & al, const D & bl,  const D & k, const D & cx, const D & cxp
	,const D & xbase, const D & xpbase, const D & xmax, const D & xmins, const D & xmaxs, const D & xpmin, const D & xpmax
	, const bool dbg) const
{
if (!dbg)
{
MM_ONCE(" range checking disabled as dbg is false", )
return; 
}
const D xuu=au*ulp+bu; //+xbase; // this must be xbase relative 
const D xul=al*ulp+bl;//+xbase; // this must be xbase relative 
const D xlu=au*llp+bu;//+xbase; // this must be xbase relative 
const D xll=al*llp+bl; //+xbase; // this must be xbase relative 

const D kiuu=k+cx*xuu+cxp*ulp;
const D kiul=k+cx*xul+cxp*ulp;
const D kilu=k+cx*xlu+cxp*llp;
const D kill=k+cx*xll+cxp*llp;

const D xuub=xuu+xbase;
const D xulb=xul+xbase;
const D xlub=xlu+xbase;
const D xllb=xll+xbase;

bool bad=false; 
bool realbad=false; 

// tolerance here is a huge problem 
if (kiuu<0) {bad=true; } 
if (kiul<0) {bad=true; } 
if (kilu<0) {bad=true; } 
if (kill<0) {bad=true; } 
const D kerr=-1e-4; // -1e-6*(xmaxs-xmins);
if (kiuu<kerr) {realbad=true; } 
if (kiul<kerr) {realbad=true; } 
if (kilu<kerr) {realbad=true; } 
if (kill<kerr) {realbad=true; } 


if (realbad&&bad)
{
MM_MSG("rangeck "<<lbl<<" "<<MMPR(kiuu)<<MMPR(kiul)<<MMPR(kilu)<<MMPR(kill))
MM_MSG("rangeck "<<lbl<<" "<<MMPR(xbase)<<MMPR(xpbase)<<MMPR(ulp)<<MMPR(llp)<<MMPR(au)<<MMPR(bu)<<MMPR(al)<<MMPR(bl))
// we need to supply a cm pointer.. 
//MM_INC_MSG(m_cm,"rangeck" ) 
// rangeck(const StrTy & lbl,const D & ulp, const D & llp
//	, const D & au, const D & bu, const D & al, const D & bl,  const D & k, constg D & cx, const D & cxp
//	,const D & xbase, const D & xpbase, const D & xmax, const D & xmins, const D & xmaxs, const D & xpmin, const D & xpmax
//	, const bool dbg)

if (dbg) if (realbad) raise(SIGINT); 
}



} //rangeck
// try to adapt diff_sf_bc to drift_sf_bc
// this appears equal to the earlier drift_sf but more compact
// although the other code may be faster and even easier to understand
/*
The drift integral turned out to be a lot simpler than 
i thought and most of this is a waste. 

*/
class geo_set
{
public:
geo_set(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2, const D & _xmax)
//: m_xp0(xp0),m_xp1(xp1), m_xp2(xp2), m_x0(x0),m_x1(x1), m_x2(x2), xmax(_xmax) {}
: m_xp0(xp0),m_xp1(xp1), m_xp2(xp2), m_x0(x0),m_x1(x1), m_x2(x2), m_vt(_xmax) {}

void set_sd(const bool src_right, const bool dest_right)
{
if (!src_right)  { xmins=m_x0; xmaxs=m_x1; xbase=m_x0;} else { xmins=m_x1; xmaxs=m_x2; xbase=m_x2; }
if (!dest_right)  { xpmin=m_xp0; xpmax=m_xp1;} else { xpmin=m_xp1; xpmax=m_xp2; }
//const D dxbase=xpbase-xbase; // add to x' to convert to x
// how make loop consts const?? 
xpbase=m_xp0;dxbase=xpbase-xbase; // add to x' to convert to x
// strict equality should catch self as the limits should jsut be 
// bitcopied
self=(xpmin==xmins)&&(xpmax==xmaxs);

if (xpmax<=xmins) {m_case=0; return; } // case 0, x' is to the left of x the src region 
if ((xpmin<xmins)&&(xpmax>xmins)&&(xpmax<=xmaxs)){ m_case=1; return; }  // case 1, x' is left of x and into part of it  
if ((xpmin<=xmins)&&(xpmax>xmaxs)) { m_case=2; return; } // case 2, x' spans x the src region   
if ((xpmin>=xmins)&&(xpmax<=xmaxs)){ m_case=3; return; }  // case 3a, x' is contained in the src region   
if ((xpmin<xmaxs)&&(xpmin>=xmins)&&(xpmax>xmaxs)) {m_case=4; return; } // case 4, x' overlaps and extends to the right of src region   
if (xpmin>=xmaxs) {m_case=5; return; }  // case 5, x' is to the right  of x the src region 
MM_ERR(" fall thgourh ") 

}
StrTy to_string() const
{
Ss s;
s<<MMPR(m_xp0);
s<<MMPR(m_xp1);
s<<MMPR(m_xp2);
s<<MMPR(m_x0);
s<<MMPR(m_x1);
s<<MMPR(m_x2);
s<<MMPR(m_vt);
s<<MMPR(xmins);
s<<MMPR(xmins);
s<<MMPR(xmaxs);
s<<MMPR(xpmin);
s<<MMPR(xpmax);
s<<MMPR(xbase);
s<<MMPR(xpbase);
s<<MMPR(dxbase);
s<<MMPR(m_case);
s<<MMPR(self);

return s.str();

}
//D m_xp0,m_xp1, m_xp2, m_x0,m_x1, m_x2,xmax;
D m_xp0,m_xp1, m_xp2, m_x0,m_x1, m_x2,m_vt;
D xmins,xmaxs,xpmin,xpmax,xbase,xpbase,dxbase;
IdxTy m_case;
bool self;

}; // geo_set

class seg_set
{
public:

bool start( const geo_set & g, const bool no_self) 
{
// integrate over the variable x' which indexes the s shape function
// from xpmin to xpmax but limited by the offset vt and xmaxs and xmins
if (g.self&&no_self) return true; // skip 
 llp=g.xpmin-g.xpbase;  // lower limit of x' integration relaive to base
 ulp=g.xpmax-g.xpbase; //upper x' limit 
D srcmax=g.xmaxs-g.xpbase+g.m_vt;
D srcmin=g.xmins-g.xpbase+g.m_vt;
if (llp<srcmin) llp=srcmin;
if (ulp>srcmax) ulp=srcmax;
return (ulp<=llp);

}


// with the variable origin chosen for numerical issues, get 
// the right shape functions for that origin with possibly the
// c src function being forced to uniform. 
template <class Tt > void make_sf(const Tt & that, const geo_set & g,const bool src_right,const bool dest_right, const bool uniform_src)
{
// gee, make a member of g? lol 
that.make_sf(s1,s0, g.m_xp0, g.m_xp1, g.m_xp2,(dest_right)?1:0);
if (uniform_src) { c1=0; c0=1; } 
else that.make_sf(c1,c0, g.m_x0, g.m_x1, g.m_x2,(src_right)?1:0,g.xbase);

}

StrTy to_string() const
{
Ss s;
s<<MMPR(llp);
s<<MMPR(ulp);
s<<MMPR(s1);
s<<MMPR(s0);
s<<MMPR(c1);
s<<MMPR(c0);
return s.str();
}

D ulp,llp;
D au; D al;
D bu; D bl;
// c0, const term of src shape function
D c0; D c1; // x term of src shape function
D cx; // x coef in velocity distron
D cxp; // x coef in velocity distron
D k; // xmax+c1*x0+cxp*xp0 ; // const part of verlocity distro 
// cpeffs of destination shape function
D s0; D s1; 
 D cs; // normalization const for integral

} ; // set_set


// 
class drift_flags
{
public:
drift_flags(): no_self(false), uniform_src(true), dbg(false)
,only_left_src(false), only_right_src(false) {}
drift_flags(const bool ns, const bool us, const bool db):
 no_self(ns), uniform_src(us), dbg(db)
,only_left_src(false), only_right_src(false) {}


bool no_self;
bool uniform_src;
bool dbg; 
bool only_left_src;
bool only_right_src;
}; // drift_flags

/*
Main code for drift integral. Move from source x into x' as 
translated by amount vt. This is basically a convolution 
*/

D drift_sf_bc(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	//, const D & n0, const D & n1, const D & n2,  const D & vt, const bool no_self, const bool dbg) const
	, const D & n0, const D & n1, const D & n2,  const D & vt, 
const drift_flags & flags) const
{
D res=0;
const bool db=flags.dbg;
const bool dbg_self=(false&&(xp0==x0)); // for testing uniform_src self term 
typedef geo_set Geo;
typedef seg_set Seg;
// this needs a different vt for left and right sides 
Geo geo(xp0,xp1, xp2,x0,x1, x2,vt);
Seg sego;
const bool lin_fit=true;
D A,B,C;
if (!lin_fit) fitq(A,B,C,x0,x1,x2,n0,n1,n2);

for (IdxTy seg=0; seg<4; ++seg)
{
const bool  src_right=((seg&1)!=0);
if (flags.only_left_src&&src_right) continue;
if (flags.only_right_src&&!src_right) continue;

const bool  dest_right=((seg&2)!=0);
// xbase is zero for "x" and is  nearest boundary x0 or x2
geo.set_sd(src_right,dest_right);
// create ulp and llp check for overlap 
if (flags.no_self) { MM_ONCE(" no_self set so ignore self diffusion",) } 
if (sego.start(geo,flags.no_self)) continue;
// make part of geo? 
if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,geo.xbase);
sego.make_sf(*this, geo,src_right,dest_right, flags.uniform_src);
if (dbg_self)
{
MM_MSG(" attempting to get non-zero drift for nigh v self ????")
MM_MSG(" geo "<<geo.to_string()<<" sego "<<sego.to_string() ) 
}

// convert the s (x') shape function to coords x'- vt
// and then find limits of integration and do it analytically.

// this should just be a simple integral now, and A==0 lol.
// but convert it to x' offset by vt
const D k=geo.xbase-geo.xpbase+geo.m_vt;
const D s1p=sego.s1;
const D s0p=sego.s0+sego.s1*k;
const D P4=A*s1p*sego.c1;
const D P3=A*s1p*sego.c0+A*s0p*sego.c1+B*s1p*sego.c1;
// doh typo, caused center to zero out lol 
//const D P2=A*s0p*sego.c0+B*s1p*sego.c0+B*s0p*sego.c1+C*s1p*sego.c0;
const D P2=A*s0p*sego.c0+B*s1p*sego.c0+B*s0p*sego.c1+C*s1p*sego.c1;
const D P1=B*s0p*sego.c0+C*s1p*sego.c0+C*s0p*sego.c1;
const D P0=C*s0p*sego.c0;
if (flags.dbg)
{
MM_MSG(MMPR(k)<<MMPR(vt)<<MMPR(geo.xbase)<<MMPR(geo.xpbase)<<MMPR(s1p)<<MMPR(s0p))
MM_MSG(MMPR(geo.xmins)<<MMPR(geo.xmaxs)<<MMPR(geo.xpmin)<<MMPR(geo.xpmax))
MM_MSG(MMPR(sego.s1)<<MMPR(sego.s0)<<MMPR(sego.c1)<<MMPR(sego.c0))
MM_MSG(MMPR(P4)<<MMPR(P3)<<MMPR(P2)<<MMPR(P1)<<MMPR(P0)<<MMPR(A)<<MMPR(B)<<MMPR(C))
}
// the whole integrand is now in terms of x+vt NOT x' 
// this should be zero at xbase, negative values are bad... 
D xeff=sego.llp+geo.xpbase-geo.xbase-vt;
static const D c13=1.0/3.0;
const D p1=((((.2*P4*xeff+.25*P3)*xeff+c13*P2)*xeff+.5*P1)*xeff+P0)*xeff;
if (flags.dbg) { MM_MSG(MMPR(xeff)<<MMPR(sego.llp)<<MMPR(p1)) }
 xeff=sego.ulp+geo.xpbase-geo.xbase-vt;
const D p2=((((.2*P4*xeff+.25*P3)*xeff+c13*P2)*xeff+.5*P1)*xeff+P0)*xeff;
if (flags.dbg) { MM_MSG(MMPR(xeff)<<MMPR(sego.ulp)<<MMPR(p2)) }
const D segr=p2-p1; // missing anything? 
if (segr>0) res+=segr;
else
{

//MM_MSG(" negative segr "<<MMPR(segr))
/*
MM_MSG(MMPR(k)<<MMPR(vt)<<MMPR(geo.xbase)<<MMPR(geo.xpbase)<<MMPR(s1p)<<MMPR(s0p))
MM_MSG(MMPR(geo.xmins)<<MMPR(geo.xmaxs)<<MMPR(geo.xpmin)<<MMPR(geo.xpmax))
MM_MSG(MMPR(sego.s1)<<MMPR(sego.s0)<<MMPR(sego.c1)<<MMPR(sego.c0))
MM_MSG(MMPR(P4)<<MMPR(P3)<<MMPR(P2)<<MMPR(P1)<<MMPR(P0)<<MMPR(A)<<MMPR(B)<<MMPR(C))

{ MM_MSG(MMPR(xeff)<<MMPR(sego.llp)<<MMPR(p1)) }
{ MM_MSG(MMPR(xeff)<<MMPR(sego.ulp)<<MMPR(p2)) }
*/
}

} // seg

return res;
} // drift_sf_bc

/*
Main FL diffusion calculation from element x to xp with carriers n
and displacement profile linear with length xmax. 

*/

D diff_sf_bc(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2,  const D & xmax, const bool no_self, const bool dbg) const
{
D res=0;
const bool db=dbg;
if (xmax<=0)
{
MM_ERR(" xmax is unusual, make sure parameters are right ")
{MM_ERR("x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" xp0="<<xp0<<" xp1="<<xp1<<" xp2="<<xp2
<<" n0="<<n0<<" n1="<<n1<<" n2="<<n2<<" xmax="<<xmax<<" no_self="<<no_self<<" dbg="<<dbg) }

}
const bool lin_fit=true;
D A,B,C;
if (!lin_fit) fitq(A,B,C,x0,x1,x2,n0,n1,n2);
// limits of inner x integral, a polynomial in x=a*x'+b referienced to respecitve
// zero bases xbase and xp0. 
D au=0; D al=0;
D bu=0; D bl=0;
// c0, const term of src shape function
D c0=  0; D c1=  0; // x term of src shape function
D cx= 1; // x coef in velocity distron
D cxp= -cx; // x coef in velocity distron
D k=0; // xmax+c1*x0+cxp*xp0 ; // const part of verlocity distro 
// cpeffs of destination shape function
D s0=0; D s1=0; 
const D cs=1; // normalization const for integral
// each element has a left and right part for n and shape functions
for (IdxTy seg=0; seg<4; ++seg)
{
const bool  src_right=((seg&1)!=0);
const bool  dest_right=((seg&2)!=0);
// xbase is zero for "x" and is  nearest boundary x0 or x2
D xmins=0,xmaxs=0,xpmin=0,xpmax=0,xbase=0,xpbase=0;
xpbase=xp0;
if (!src_right)  { xmins=x0; xmaxs=x1; xbase=x0;} else { xmins=x1; xmaxs=x2; xbase=x2; }
if (!dest_right)  { xpmin=xp0; xpmax=xp1;} else { xpmin=xp1; xpmax=xp2; }
const D dxbase=xpbase-xbase; // add to x' to convert to x

 D llp=xpmin-xpbase;  // lower limit of x' integration relaive to base
 D ulp=xpmax-xpbase; //upper x' limit 
if (db) { MM_MSG("******  seg "<<seg<<" "<<MMPR(src_right)<<MMPR(dest_right)) }
// these can not extend more than xmax away from the source, offset by bases 
const D maxxp=xmaxs+xmax-xpbase;
const D minxp=xmins-xmax-xpbase;
if (ulp> (maxxp) ){ulp=maxxp;}
if (ulp< (minxp) ){continue;}
if (llp< (minxp) ){llp=minxp;}
if (llp> (maxxp) ){continue;}
if (ulp<=llp) continue;

if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,xbase);
// the shape functions are const over the segment
make_sf(s1,s0, xp0, xp1, xp2,(dest_right)?1:0);
make_sf(c1,c0, x0, x1, x2,(src_right)?1:0,xbase);
// siz mutually exclusive cases, each with some subregions possible 



if (xpmax<=xmins) // case 0, x' is to the left of x the src region 
{
if (db) { MM_MSG(" CASE 0 "<<MMPR(xpmax)<<MMPR(xmins)) }
// x>x', x-x'<xmax
cxp=1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
// depending on xmax, there could be three regions
const D bpmax=xmaxs-xmax-xpbase;
const D bpmin=xmins-xmax-xpbase;
D llpx=llp;
D ulpx=ulp;
// for x'>bpmax, the upper x limit is fixed at xmaxs 
if (bpmax<ulpx) ulpx=bpmax;
if (bpmin>llpx) llpx=bpmin;
// the integrand is only nonzero for x'>bpmin
// should hve been checked earlier 
if (ulp<bpmin ) continue;  // can not return 0 doh 
// the maximum x is x'+xmax 
// the lowest x is the lower range of the src or else everthing is oor and zero 
al=0; bl=xmins-xbase;
// a middle range is only in range for part of src 
if (ulpx<llpx) ulpx=llpx;
if (ulpx>llpx)
{
	au=1; bu=xmax+dxbase;
	corners("Case0s1", res, cs,ulpx, llpx, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
}
// the entier src region contributes so the upper limiton x is xmaxs
if (ulpx<ulp)
{
au=0; bu=xmaxs-xbase;
	corners("Case0s2", res, cs,ulp, ulpx, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
}
continue; // return res;
} // case  0

if ((xpmin<xmins)&&(xpmax>xmins)&&(xpmax<=xmaxs)) // case 1, x' is left of x and into part of it  
{
if (db) { MM_MSG(" CASE 1a x' left and olap x  "<<MMPR(xpmax)<<MMPR(xmins)) }
// there is always some segment for x>x'
// x>x', x-x'<xmax
cxp=1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
// xp region, further limited by xmax away from x region boundaries. 
D llpx=llp;
D ulpx=ulp;
D bpxl=xmins-xpbase;
D bpxu=xmaxs-xmax-xpbase;
D clim=llpx;
while (clim<ulpx)
{
D uulim=ulpx;

if (bpxl>clim) if (bpxl<uulim) uulim=bpxl;
if (bpxu>clim) if (bpxu<uulim) uulim=bpxu;

const D ulim=uulim;
const D llim=clim;
if (ulim<=clim)
{
MM_MSG(" loop stuck at "<<MMPR(ulim)<<MMPR(llim)<<MMPR(ulp)<<MMPR(llp)<<MMPR(bpxl)<<MMPR(bpxu))
break;
}
if (ulim>clim)
{
const bool ubpxu=((ulim>=bpxu)&&(llim>=bpxu));
const bool lbpxu=(ulim<=bpxu)&&(llim<=bpxu);
const bool ubpxl=((ulim>=bpxl)&&(llim>=bpxl));
const bool lbpxl=(ulim<=bpxl)&&(llim<=bpxl);
if (ubpxl) { al=1; bl=dxbase; }// lower x limit is lower range of src
else if (lbpxl) { al=0; bl=xmins-xbase; }//but x no lower than this 
else { MM_MSG(" bad limits "); raise(SIGINT); } 

if (ubpxu){ au=0; bu=xmaxs-xbase; }//x no larger than this 
else if (lbpxu) { au=1; bu=xmax+dxbase; }//x no more than xmax greter than x' 
else { MM_MSG(" bad limits "); raise(SIGINT); } 

StrTy lbl="Case1axx";
lbl=lbl+StrTy(ubpxu?"1":"0");
 lbl=lbl+StrTy(lbpxu?"1":"0");
 lbl=lbl+StrTy(ubpxl?"1":"0");
 lbl=lbl+StrTy(lbpxl?"1":"0");
corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim, llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);

clim=ulim;
} // if  
} // while 
/////////////////////////////////////





// x'>x
cxp=-1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 

 llpx=llp;
 ulpx=ulp;
// this relies on the above bp to prevent 
if (llpx<bpxl) llpx=bpxl; 
// lower x lim is xp-xmax or xmins
 bpxl=xmins+xmax-xpbase;
// x upper lim is x' until x'>xmax
 bpxu=xmaxs-xpbase;
 clim=llpx;
while (clim<ulpx)
{
D uulim=ulpx;
if (bpxl>clim) if (bpxl<uulim) uulim=bpxl;
if (bpxu>clim) if (bpxu<uulim) uulim=bpxu;
const D ulim=uulim;
const D llim=clim;
if (ulim<=clim) {
MM_MSG(" loop stuck at "<<MMPR(ulim)<<MMPR(llim)<<MMPR(ulp)<<MMPR(llp)<<MMPR(bpxl)<<MMPR(bpxu))
break;
}
if (ulim>clim)
{
const bool ubpxu=((ulim>=bpxu)&&(llim>=bpxu));
const bool lbpxu=(ulim<=bpxu)&&(llim<=bpxu);
const bool ubpxl=((ulim>=bpxl)&&(llim>=bpxl));
const bool lbpxl=(ulim<=bpxl)&&(llim<=bpxl);
if (lbpxu) { au=1; bu=dxbase; }// x no bigger than xp  
else if (ubpxu) { au=0; bu=xmaxs-xbase; }// x no bigger than xp  
else { MM_MSG(" bad limits "); raise(SIGINT); } 

if (ubpxl) { al=1; bl=-xmax+dxbase; }// x no bigger than xp  
else if (lbpxl) { al=0; bl=xmins-xbase; }// x no bigger than xp  
else { MM_MSG(" bad limits "); raise(SIGINT); } 

StrTy lbl="Case1ayx";
lbl=lbl+StrTy(ubpxu?"1":"0");
 lbl=lbl+StrTy(lbpxu?"1":"0");
 lbl=lbl+StrTy(ubpxl?"1":"0");
 lbl=lbl+StrTy(lbpxl?"1":"0");
corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim, llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);

clim=ulim;
} // if  
} // while 


//////////////////////

continue;
} // case1a


if ((xpmin<=xmins)&&(xpmax>xmaxs)) // case 2, x' spans x the src region   
{
if (db) { MM_MSG(" CASE 2a x' contains x  "<<MMPR(xpbase)<<MMPR(xpmin)<<MMPR(xpmax)<<MMPR(xbase)<<MMPR(xmins)<<MMPR(xmaxs)) }
// ulp is either xpmax-xpbase or xmaxs+xmax-xbase
// llp is either xpmin-xpbase or xmin-xmax-xbase
D ulp=xpmax-xpbase;
D llp=xpmin-xpbase;
D ulpx=xmaxs+xmax-xbase;
if (ulpx<ulp) ulp=ulpx;
D llpx=xmins-xmax-xbase;
if (llpx>llp) llp=ulpx;

// x>x', starting from the left
cxp=1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
// x>x'
al=1; bl=dxbase;
// the lower limit on xp is whatever llp we had.
// the upper limit is xmax
// the upper limit on x is either xmaxs or xp+xmax
// 
// for x>x', the lower limit on x' must by xmins
D ulpxu=xmaxs-xpbase;
// llpx=llp; // xmins-xpbase;
 ulpx= xmaxs-xmax-xpbase;
if (ulpx>ulpxu) ulpx=ulpxu;
if (ulpx<llp) ulpx=llp;

if (ulpx>llp)
{ 
// for xp small enough, the upper limit on x ix xp+xmaax
au=1; bu=xmax+dxbase;
const StrTy lbl="Case2as1a";
if (dbg) { MM_MSG(lbl<<" "<<MMPR(xbase)<<MMPR(xpbase)) } 
corners(lbl, res, cs,ulpx, llp, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulpx, llp, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}
if (ulpx<ulpxu)
{ 
// close to the end, ul is just xmaxs
au=0; bu=xmaxs-xbase;
const StrTy lbl="Case2as1a";
if (dbg) { MM_MSG(lbl<<" "<<MMPR(xbase)<<MMPR(xpbase)) } 
corners(lbl, res, cs,ulpxu, ulpx, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulpxu, ulpx, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}



// x<x'
cxp=-1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
au=1; bu=dxbase;
// as above just reverse everything...
 llp=xmins-xpbase;
ulp=xmaxs+xmax-xpbase;
if (ulp>(xpmax-xpbase)) ulp=xpmax-xpbase;

// integration is now from llp to ulp.
// lower x limit is either xmins( xp<xmins+xmax)
// or xp-xmax (xp>xmins+xmax
D bpxl=xmins+xmax-xpbase;
if (bpxl>ulp) bpxl=ulp;
// the mxximum x is xp until getitng to xmaxs
D bpxu=xmaxs-xpbase;
// depending on xmax vs (xmaxs-xmins), orders will vary

D cllp=llp;
while (cllp<ulp)
{
// jj
D uulim=ulp;
if (bpxl>cllp) if (bpxl<uulim) uulim=bpxl;
if (bpxu>cllp) if (bpxu<uulim) uulim=bpxu;

const D ulim=uulim;
const D llim=cllp;
const bool ubpxu=((ulim>=bpxu)&&(llim>=bpxu));
const bool lbpxu=(ulim<=bpxu)&&(llim<=bpxu);
const bool ubpxl=((ulim>=bpxl)&&(llim>=bpxl));
const bool lbpxl=(ulim<=bpxl)&&(llim<=bpxl);
if (ulim<llim)
{
MM_MSG(" loop stuck at "<<MMPR(ulim)<<MMPR(llim)<<MMPR(ulp)<<MMPR(llp)<<MMPR(bpxl)<<MMPR(bpxu))
break; 
}
if (ulim>=llim)
{
StrTy lbl=StrTy("Case2asy1v");
lbl=lbl+StrTy(ubpxu?"1":"0");
 lbl=lbl+StrTy(lbpxu?"1":"0");
 lbl=lbl+StrTy(ubpxl?"1":"0");
 lbl=lbl+StrTy(lbpxl?"1":"0");
// x==xp, xp<xmaxs
if (lbpxu) {  au=1; bu=dxbase; } 
// x must be no greater than xmaxs
else if (ubpxu ) { au=0; bu=xmaxs-xbase;} // xp>xmaxs
else { MM_MSG(" bad limits "); raise(SIGINT); } 
// x lower limit is xp-xmax which must be more than xmins
if (ubpxl) { al=1; bl=-xmax+dxbase;} // xp>xmin+xmax
else if (lbpxl) {al=0; bl=xmins-xbase;}
else { MM_MSG(" bad limits "); raise(SIGINT); } 
corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim,llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
cllp=ulim;
}
}

continue;
} // case2


if ((xpmin>=xmins)&&(xpmax<=xmaxs)) // case 3a, x' is contained in the src region   
{

if (db) { MM_MSG(" CASE 3a x' contained in x  "<<MMPR(xpbase)<<MMPR(xpmin)<<MMPR(xpmax)<<MMPR(xbase)<<MMPR(xmins)<<MMPR(xmaxs)) }
// x' > x
cxp=-1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
ulp=xpmax-xpbase;
llp=xpmin-xpbase;
const D bp1=xmins+xmax-xpbase; // value of xp at bl switch 

D clim = llp;
while (clim<ulp)
{
D uulim=ulp;
if (bp1>clim) if (bp1<uulim) uulim=bp1;
const D ulim=uulim;
const D llim=clim;
const bool uxl=(ulim<=bp1)&&(llim<=bp1);
const bool lxl=(ulim>=bp1)&&(llim>=bp1);
au=1; bu=dxbase; // x less tham x'
//au=0; bu=xmaxs-xbase; // this limit can not be raeched  
if (lxl) {al=1; bl=-xmax+dxbase;} // x closer than xmax to x' 
else if (uxl) { al=0; bl=xmins-xbase;} // x no smaller than min  
else 
{
MM_MSG(" bad limits "<<MMPR(ulim)<<MMPR(clim)<<MMPR(bp1))
raise(SIGINT); 
}
if (ulim>llim)
{
StrTy lbl="Case3as1x";
lbl=lbl+StrTy(uxl?"1":"0");
lbl=lbl+StrTy(lxl?"1":"0");

corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim,llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}
else
{
MM_MSG(" bad limits "<<MMPR(ulim)<<MMPR(clim)<<MMPR(bp1))
}

clim=ulim;
}

// x> x'
cxp=1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
ulp=xpmax-xpbase;
llp=xpmin-xpbase;

const D bp2=xmaxs-xmax-xpbase; // value of xp at bu switch 

 clim = llp;
while (clim<ulp)
{
D uulim=ulp;
if (bp2>clim) if (bp2<uulim) uulim=bp2;
const D ulim=uulim;
const D llim=clim;
const bool uxl=(ulim<=bp2)&&(llim<=bp2);
const bool lxl=(ulim>=bp2)&&(llim>=bp2);
al=1; bl=dxbase; // x > x'
//au=0; bu=xmaxs-xbase; // this limit can not be raeched  
if (uxl) {au=1; bu=xmax+dxbase;} // x closer than xmax to x' 
else if (lxl) { au=0; bu=xmaxs-xbase;} // x no smaller than min  
else 
{
MM_MSG(" bad limits "<<MMPR(ulim)<<MMPR(clim)<<MMPR(bp1))
raise(SIGINT); 
}
if (ulim>llim)
{
StrTy lbl="Case3as2x";
lbl=lbl+StrTy(uxl?"1":"0");
lbl=lbl+StrTy(lxl?"1":"0");

corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim,llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}
else
{
MM_MSG(" bad limits "<<MMPR(ulim)<<MMPR(clim)<<MMPR(bp1))
}

clim=ulim;
}


continue; // return res;
} // case3a


if ((xpmin<xmaxs)&&(xpmin>=xmins)&&(xpmax>xmaxs)) // case 4, x' overlaps and extends to the right of src region   
{
if (db) { MM_MSG(" CASE 4a "<<MMPR(xpmax)<<MMPR(xmins)) }
// x> x'
cxp=1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
llp=xpmin-xpbase;
ulp=xpmax-xpbase;
const D ulp2=xmaxs+xmax-xpbase;
const D ulp3=xmaxs-xpbase;
// upper limit no higher than src region max plus xmax
if ( ulp>ulp2) ulp=ulp2;
D llpx=llp;
D ulpx=ulp;
// for x>x' no point integrating beyond xmaxs
if (ulpx>ulp3) ulpx=ulp3;
//the x src var is at lest x' and not less than xmins ( should coincide with xpmins )  
D bpxl=xmins-xpbase;
// the upper x limit is either xmaxs or xp+xmax
D bpxu=xmaxs-xmax-xpbase;
D clim=llpx;
while (clim<ulpx)
{
D uulim=ulpx;
// this should never trigger as clim is bpxl
if (bpxl>clim) if (bpxl<uulim) uulim=bpxl;
	// when closer to xmaxs than xmax, the upper x limit changes 
if (bpxu>clim) if (bpxu<uulim) uulim=bpxu;

const D ulim=uulim;
const D llim=clim;
if (ulim<=clim)
{
MM_MSG(" loop stuck at "<<MMPR(ulim)<<MMPR(llim)<<MMPR(ulp)<<MMPR(llp)<<MMPR(bpxl)<<MMPR(bpxu))
break;
}
if (ulim>clim)
{
const bool ubpxu=((ulim>=bpxu)&&(llim>=bpxu));
const bool lbpxu=(ulim<=bpxu)&&(llim<=bpxu);
const bool ubpxl=((ulim>=bpxl)&&(llim>=bpxl));
const bool lbpxl=(ulim<=bpxl)&&(llim<=bpxl);
//  this should remain fixed 
if (ubpxl) { al=1; bl=dxbase; }// lower x limit is lower range of src
else if (lbpxl) { al=0; bl=xmins-xbase; }//but x no lower than this 
else { MM_MSG(" bad limits "); raise(SIGINT); } 
// closer to xmaxs than xmax, fix the limit to xmaxs
if (ubpxu){ au=0; bu=xmaxs-xbase; }//x no larger than this 
else if (lbpxu) { au=1; bu=xmax+dxbase; }//x no more than xmax greter than x' 
else { MM_MSG(" bad limits "); raise(SIGINT); } 

StrTy lbl="Case4axx";
lbl=lbl+StrTy(ubpxu?"1":"0");
 lbl=lbl+StrTy(lbpxu?"1":"0");
 lbl=lbl+StrTy(ubpxl?"1":"0");
 lbl=lbl+StrTy(lbpxl?"1":"0");
corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim, llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);

clim=ulim;
} // if  
} // while 
/////////////////////////////////////




///////////////////////////


// x<x' x' extends to the right of x which is now a relevant region 
cxp=-1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 


 llpx=llp;
 ulpx=ulp; // this now extends as far as useful not just to x'
 bpxl=xmaxs-xpbase;
 bpxu=xmins+xmax-xpbase;
 clim=llpx;
while (clim<ulpx)
{
D uulim=ulpx;

if (bpxl>clim) if (bpxl<uulim) uulim=bpxl;
if (bpxu>clim) if (bpxu<uulim) uulim=bpxu;

const D ulim=uulim;
const D llim=clim;
if (ulim<=clim)
{
MM_MSG(" loop stuck at "<<MMPR(ulim)<<MMPR(llim)<<MMPR(ulp)<<MMPR(llp)<<MMPR(bpxl)<<MMPR(bpxu))
break;
}
if (ulim>clim)
{
const bool ubpxu=((ulim>=bpxu)&&(llim>=bpxu));
const bool lbpxu=(ulim<=bpxu)&&(llim<=bpxu);
const bool ubpxl=((ulim>=bpxl)&&(llim>=bpxl));
const bool lbpxl=(ulim<=bpxl)&&(llim<=bpxl);
// x can not be greater than xp
if (lbpxl) { au=1; bu=dxbase; }//
// and x can not be greater than xmaxs 
else if (ubpxl) { au=0; bu=xmaxs-xbase; }// 
else { MM_MSG(" bad limits "); raise(SIGINT); } 

if (lbpxu){ al=0; bl=xmins-xbase; }//x no larger than this 
else if (ubpxu) { al=1; bl=-xmax+dxbase; }//x no more than xmax greter than x' 
else { MM_MSG(" bad limits "); raise(SIGINT); } 

StrTy lbl="Case4axx";
lbl=lbl+StrTy(ubpxu?"1":"0");
 lbl=lbl+StrTy(lbpxu?"1":"0");
 lbl=lbl+StrTy(ubpxl?"1":"0");
 lbl=lbl+StrTy(lbpxl?"1":"0");
corners(lbl, res, cs,ulim, llim, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck(lbl,ulim, llim, au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);

clim=ulim;
} // if  
} // while 
/////////////////////////////////////



continue;
} // case4 a


if (xpmin>=xmaxs) // case 5, x' is to the right  of x the src region 
{
if (db) { MM_MSG(" CASE 5 "<<MMPR(xpbase)<<MMPR(xpmin)<<MMPR(xpmax)<<MMPR(xbase)<<MMPR(xmins)<<MMPR(xmaxs)) }
// x<x', x'-x<xmax
cxp=-1; cx=-cxp;
k=xmax+cx*xbase+cxp*xpbase ; // const part of verlocity distro 
// depending on xmax, there could be three regions
llp=xpmin-xpbase;
ulp=xpmax-xpbase;
//D llpx=llp;
//D ulpx=ulp;
// for xp>bpmax, the integrand is zero 
const D bpmax=xmaxs+xmax-xpbase;
// for xp>bpmin, the lower x limit is now xp-xmax instead of xmins
const D bpmin=xmins-xpbase;
if (ulp>bpmax) ulp=bpmax;
if (llp<bpmin) llp=bpmin;
//ulpx=ulp;
D bpxp=xmins+xmax-xpbase;
if (bpxp>ulp) bpxp=ulp;
if (bpxp<llp) bpxp=llp;
// the natural limits are the ends of the destination region clipped by 
// limits on source and xmax
// for x'>bpmax, the upper x limit is fixed at xmaxs 
// the integrand is only nonzero for x'<bpmax
// should hve been checked earlier 
if (llp>ulp ) return 0; 

// integration in xp extends from  xpmin to min(xpmax, xmaxs+xmax) and the
// lower limit on x is  xp-xmax or xmins
// these are the natural limits, 
//au=1; bu=xmax+dxbase;
au=0; bu=xmaxs-xbase;
// in all cases, xp>x
if (bpxp>llp) {
al=0; bl=xmins-xbase;
if (dbg) MM_MSG(" Case5s0     xmins<x<xp+xmax  note x'>x")
corners("Case5s0", res, cs,bpxp, llp, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck("Case5s0",bpxp,llp,au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}

if (bpxp<ulp)
{
if (db) MM_MSG(" Case5s1     xp-xmax<x<xmaxs  note x'>x")
al=1; bl=-xmax+dxbase;
corners("Case5s1", res, cs,ulp, bpxp, au, bu, al, bl, c0, c1, k, cx, cxp, s1, s0, A, B, C, dbg);
rangeck("Case5s1",ulp,bpxp,au, bu, al, bl, k,cx,cxp,xbase,xpbase,xmax, xmins,xmaxs,xpmin,xpmax, dbg);
}
continue; // return res;
} // case 5 
MM_ERR(" danger will robinson fall through "<<MMPR(xmins)<<MMPR(xmaxs)<<MMPR(xpmin)<<MMPR(xpmax))
MM_MSG(" danger will robinson fall through "<<MMPR(xmins)<<MMPR(xmaxs)<<MMPR(xpmin)<<MMPR(xpmax))

} // seg 

return res;
} // diff_sf_bc


//////////////////////////////////////////////////////////////////////////
#if 0

D ax=(x1>xp1)?x0:x1;
// this should be the distance of closest approach 
// for disjoint, either |    |       |     |
//                     xp0  xp2      x0    x2
// so this is x0-xp2 or
// for disjoint, either |    |       |     |
//                     x0    x2     xp0   xp2
// xp0-x2

#endif 

///////////////////////////////////////////

private:


}; // quadratic_fick_lanc_element
///////////////////////////////////////////////////////////////////////////////////////////



class gerrymander
{

typedef gerrymander Myt;
class Tr
{
public:

typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;

typedef std::stringstream Ss;

}; // Tr

// this needs a StrTy def somewhere? Where wtf? 
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;
typedef Tr::Ss Ss;

//typedef diffuse_iterator Di;
typedef diffuse_iterator_sf Di;
typedef quadratic_fick_lanc_element QFLE;
//typedef quadratic_fick_lanc_element_sf QFLESF;
//typedef quadratic_fick_lanc_element_sf QFLE;

typedef std::vector<D> Pts;
public:
gerrymander() {}

//using fixed(,ix) for positions and values sol(,iy) generate points based on
// subdividing between changes in sol(,iy)
enum { AB_CURVE=1, SLOPE=2,OOL=4};
void max_delta_refine
	(MyBlock & newpoints, const MyBlock & fixed, const MyBlock & sol, const IdxTy ix, const IdxTy iy, const D & thr, const IdxTy flags=0)
{
	const IdxTy ndims=sol.n_dims();
	if (ndims<1) {MM_ERR(" resizing from "<<ndims<<" is bad ")}
	const IdxTy oldpoints = sol.size(0);
	if (oldpoints<3) {MM_ERR(" resizing from oldpoints="<<oldpoints<<" is bad ")}
	QFLE qfle;
	Pts pts;
	pts.push_back(fixed(0,ix));
	Di di( oldpoints, ix, iy,  sol, fixed);
	while (di.ok())
	{
		D A,B,C;
		qfle.fitq(A,B,C,di);
		const D dx=di.x2()-di.x0();
		const D dn=di.n2()-di.n0();
		const D dxm=di.x1()-di.x0();
		const D dnm=di.n1()-di.n0();
		bool refine=false;
		// for now probably want to calculate each criterion and print as binary string for interactive mode 
		if (0!=(flags&&(AB_CURVE))) { refine|=((A*A*dx*dx)>(B*B));}
		if (0!=(flags&&(SLOPE))) { refine|=((B*B*dx*dx)>thr);}
		if (0!=(flags&&(OOL))) { const D d1=dnm/dxm; const D d2=dn/dx; const D ds=d1-d2; refine|=((ds*ds)>thr);}

		if (refine)
		{
			MM_MSG("refining "<<di.to_string()<<" based on "<<A<<" "<<B<<" "<<C<<" "<<dx<<" thr "<<thr)
		 	AddNPoints(pts,di.x0(),di.x1(),1); // pts.push_back(.5*(di.x2()+di.x1()));
		//	pts.push_back(.5*(di.x0()+di.x1()));
		}
		// adding points between x0 and x1
		pts.push_back(di.x1());
		// adding points between x1 and 2x
	if (refine)		 AddNPoints(pts,di.x1(),di.x2(),1); // pts.push_back(.5*(di.x2()+di.x1()));
		di.inc();
	} // while 
	pts.push_back(fixed(oldpoints-1,ix));

	newpoints=pts;
}

private:

void AddNPoints( Pts & pts, const D & xi, const D & xf, const IdxTy n)
{
const D del=(xf-xi)/D(n+1);
for (IdxTy i=1; i<=n; ++i) pts.push_back(xi+del*D(i)); 

}


}; // gerrymander





////////////////////////////////////////////////////////////////



class ficklanc // petsc_log_fd :public petsc_fd_base //  public Petnl::jr_base
{
/*




*/


class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef ficklanc Myt;
//typedef petsc_fd_base Super;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;


	typedef  quadratic_fick_lanc_element QFLE;
typedef QFLE::drift_flags drift_flags;
typedef mjm_defect_levels Defects;
//typedef quadratic_fick_lanc_element_sf QFLESF;
//	typedef  quadratic_fick_lanc_element_sf QFLE;
//	typedef diffuse_iterator Di;
	typedef diffuse_iterator_sf Di;
public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
ficklanc():m_points(3),m_save_sol(false) {Init();}
// make this m_points now, wtf
ficklanc(const IdxTy s): m_points(s), m_save_sol(false) {Init();}
ficklanc(int argc,char **args) : m_save_sol(false)
{
m_size=1;
m_points=2000; // this really should not have a default.. 
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
m_tree.config("-tree",i,argc,args);
m_flp.config("-params",i,argc,args);
configi(m_points,"-points",i,argc,args);
m_flp.config_set("-set-param",  i,  argc, args);
m_tree.config_set("-set-branch",  i,  argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
Init();
}
////////////////////////////////////////////////////////
// command block

// this should be in the parameters map, nothing special here... 
 void configi(IdxTy & dest, const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
//const StrTy nm=StrTy(args[i]);
dest=::atoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}



// the hardcoded bersions had no benefit except maybe for speed but it was just
// extra clutter to maintain and the map is low overhead. 
// the thing needs a consistent means between config file, cmd line, and command interpretter
typedef BranchesVar BraTy;
typedef FickLancParamsVar FlpTy;

// moe this crap somehwere for others lol 
const IdxTy arbor(const StrTy & label, const bool b0, const bool b1=false, const bool b2=false, 
const bool b3=false, const bool b4=false, const bool b5=false, const bool b6=false, const bool b7=false, 
const bool b8=false, const bool b9=false, const bool b10=false, const bool b11=false ) const
{
IdxTy val=0;
IdxTy mask=1;
if (b0) val|=mask; mask<<=1;
if (b1) val|=mask; mask<<=1;
if (b2) val|=mask; mask<<=1;
if (b3) val|=mask; mask<<=1;
if (b4) val|=mask; mask<<=1;
if (b5) val|=mask; mask<<=1;
if (b6) val|=mask; mask<<=1;
if (b7) val|=mask; mask<<=1;
MM_ONCE(" made branach value "<<val<<" for "<<label, ); 
return val;
}


void command_mode()
{
//LineIterator 
CommandInterpretter li(&std::cin);
command_mode(li);

}
void command_mode(const StrTy & cmd)
{

CommandInterpretter li;
li.set(cmd,1);
command_mode(li);
}


void command_mode(CommandInterpretter & li)
{
StrTy local_label="fick";
while (li.nextok())
{
const IdxTy sz=li.size();
MM_ERR(" processing "<<li.dump())
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
if (cmd=="legend") { dump_legend(std::cout,local_label,1);  continue; } 
if (cmd=="dump-human") { dump_state_etc(std::cout,local_label,0);  continue; } 
if (cmd=="init") { Init(); /*  m_cm.clear(); */  continue; } 
if (cmd=="save-solution") { m_save_sol=true;  continue; } 
if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="test") { test(li);  continue; } 
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="quit") { clean_up(); return; } 

MM_ERR(" unrecignized command "<<li.dump())

}




} //command_mode
// when exiting the interpretter
void clean_up()
{


} 
// tests, document what this crap should do lol

void test_flow_symmetry( CommandInterpretter & li,const bool dbg=false,const D & fac=1.5 )
{
const StrTy id="test_flow_symm";
QFLE qfle;
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D A,B,C;
qfle.fitq(A,B,C, x0, x1, x2, n0, n1, n2) ;
//ngvdx(di.a(),di.b(),di.x0(),di.x1(),di.x2(), di.n0(),di.n1(),di.n2(),xmax,vt);
D a=-1; D b=0; D vt=0; D xmax=1;
for (IdxTy j=0; j<10; ++j)
{
for (D i=.1; i<10; i=fac*i)
{
xmax=i;
// removed in version 9
D f1=0; // qfle.ngvdx(a,b, x0, x1, x2, n0, n1, n2,xmax,vt,dbg);
D a2=1.0-b;
D b2=1.0-a;
D f2=0; // qfle.ngvdx(a2,b2, x0, x1, x2, n0, n1, n2,xmax,vt,dbg);
MM_MSG(id<<" a="<<a<<" b="<<b<<" a2="<<a2<<" b2="<<b2<<" x0="<<x0<<" x2="<<x2<<" xmax="<<xmax<<" f1="<<f1<<" f2="<<f2);
}
a=a-.2;
b=b-.2;
}
}

void test_diffsf( CommandInterpretter & li )
{
const StrTy id="test_driftfs";
const IdxTy sz=li.size();
//D step =d(li,  sz, 2);
//int  steps  =n(li,  sz, 3);
//if ( sz>4 ) foo  =d(li,  sz, 4);
// fudding wrong fudder 
//QFLESF qflefs;
QFLE qflefs;
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D xp0=0; D xp1=.5; D xp2=1; 
D del=.1;
D xmax=.5;
IdxTy niter=5;
D del0=1;

 { del=.1; del0=1; xmax=.5; niter=5; }
 { del=.1; del0=1; xmax=2.0; niter=25; }
 { del=.09333; del0=1; xmax=3.14159; niter=80; }
 { del=.1; del0=1; xmax=1; niter=20; }


 xp0+=del0; xp1+=del0;  xp2+=del0; 
//del=.2; xmax=2;
const bool no_self=false;
for (IdxTy i=0; i<niter;  ++i)
{
MM_MSG(" .......   new del ..................")
// diff_sf(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,dbg);
//D nf=qflefs.diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
// removed in version 9
D nf=0; // qflefs.diff_sf_loop(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax,no_self, !false); 
MM_MSG(" .......  nr integral  ..................")
//D nf=qflefs.disjoint_diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
//D nr=qflefs.diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
D nr=0; // qflefs.diff_sf_loop(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, no_self,!false); 
//D nr=qflefs.disjoint_diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
const D dxi=x0-xp2;
const D dxj=xp0-x2;
MM_MSG(" .......  nr and nf results  ..................")

MM_MSG(" x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" xp0="<<xp0<<" xp1="<<xp1<<" xp2="<<xp2<<" xmax="<< xmax) 
if ((nr+nf)!=0) 
	{ MM_MSG("dxi= "<<dxi<<" dxj="<<dxj<<" xmax="<<xmax<<" nr="<<nr<<" nf="<<nf<<" diff="<<(nf-nr)<<" "<<((nf-nr)/(nf+nr)*.5))
}

 xp0+=del; xp1+=del;  xp2+=del; 

} //i 

}
////////////////////////////////////////////////////////

void test_diffsf_olap( CommandInterpretter & li )
{
const StrTy id="test_driftfs_olap";
const IdxTy sz=li.size();
//D step =d(li,  sz, 2);
//int  steps  =n(li,  sz, 3);
//if ( sz>4 ) foo  =d(li,  sz, 4);
// fudding wrong fudder 
//QFLESF qflefs;
QFLE qflefs;
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D xp0=0; D xp1=.5; D xp2=1; 
D del=.1;
D xmax=.5;
D dxmax=.1;
IdxTy niter=5;
D del0=0;

 { del=.1; del0=1; xmax=.5; niter=5; }
 { del=.1; del0=1; xmax=2.0; niter=25; }
 { del=.09333; del0=1; xmax=3.14159; niter=80; }
 { del=0;  del0=0; xmax=.1; niter=30; }
 { del=0;  del0=0; xmax=1e-5; dxmax=1e-3; niter=30; }
const bool no_self=false;

 xp0+=del0; xp1+=del0;  xp2+=del0; 
//del=.2; xmax=2;
for (IdxTy i=0; i<niter;  ++i)
{
MM_MSG(" .......   new del ..................")
// diff_sf(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,dbg);
//D nf=qflefs.diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
// version 9 removed
D nf=0; // qflefs.diff_sf_loop(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, no_self,false); 
MM_MSG(" .......  nr integral  ..................")
//D nf=qflefs.disjoint_diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
//D nr=qflefs.diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
D nr=0; // qflefs.diff_sf_loop(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, no_self,false); 
//D nr=qflefs.disjoint_diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
const D dxi=x0-xp2;
const D dxj=xp0-x2;
MM_MSG(" .......  nr and nf results  ..................")

MM_MSG(" x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" xp0="<<xp0<<" xp1="<<xp1<<" xp2="<<xp2<<" xmax="<< xmax) 
if ((nr+nf)!=0) 
	{ MM_MSG("dxi= "<<dxi<<" dxj="<<dxj<<" xmax="<<xmax<<" nr="<<nr<<" nf="<<nf<<" diff="<<(nf-nr)<<" "<<((nf-nr)/(nf+nr)*.5))
}

 xp0+=del; xp1+=del;  xp2+=del; 
xmax+=dxmax;

} //i 

}



///////////////////////////////////////////////////////
// these are now used to make pseudo elements for the bc crap 
class pixelt { public: D x0,x1,x2; 
	void set(const D & x, const D & h) {
if (h>0) {x0=x;x1=x0+h; x2=x1+h;}
else {x2=x;x1=x2+h; x0=x1+h;}
}
	void setmid(const D & x, const D & h) {
 {x0=x-h;x1=x; x2=x+h;}
}


	void set(const D & x0, const D & x1, const D & x2) 
	{this->x0=x0; this->x1=x1; this->x2=x2; } 

void move(const D & dx) { x0+=dx; x1+=dx; x2+=dx; } 
template <class Os> Os & 
	dump(Os & os) const { os<<" x0="<<x0<<" x1="<<x1<<" x2="<<x2; return os; }
	StrTy dump() const { std::stringstream s; s<<" x0="<<x0<<" x1="<<x1<<" x2="<<x2; return s.str(); }
	StrTy dc(const pixelt & y) const 
	{ std::stringstream s; s<<(x1-y.x1); return s.str(); }
	const D dx(const pixelt & y) const { const D dx20=x2-y.x0; const D dx02=x0-y.x2;
return -666; }  
template <class Tq > const D diff_sf_bc(pixelt & xp, Tq & qflefs, const D & n0, const D & n1, const D & n2, 
const D & xmax, const bool no_selfs, const bool dbg)
{
return qflefs.diff_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
}

template <class Tq > const D drift_sf_bc(pixelt & xp, Tq & qflefs, const D & n0, const D & n1, const D & n2, 
//const D & vt, const bool no_selfs, const bool dbg)
const D & vt, const drift_flags & df)
{
//return qflefs.drift_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, no_selfs, dbg); // /xmax; 
return qflefs.drift_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, df); // /xmax; 
}

template <class Tq > const D drift_sf(pixelt & xp, Tq & qflefs, const D & n0, const D & n1, const D & n2, 
const D & vt,  const drift_flags & df, const bool dbg)
{
//return qflefs.drift_sf(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, dbg); // /xmax; 
//return qflefs.drift_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, false, dbg); // /xmax; 
//return qflefs.drift_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, drift_flags()); // /xmax; 
return qflefs.drift_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, vt, df); // /xmax; 
}

template <class Tq > const D diff_sf(pixelt & xp, Tq & qflefs, const D & n0, const D & n1, const D & n2, 
const D & xmax, const bool no_selfs, const bool dbg)
{ // removed version 9 
//return qflefs.diff_sf_loop(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
return qflefs.diff_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
}


};
/*
Check an array of fixed sized pixels and determine if charge will be conserved
in real life. Right now it is clsoe except for intermeidiate, around .5 or so xmax.



*/

void test_gr_integrator( CommandInterpretter & li )
{
//bool dbg=false;
const IdxTy ncol=0;
const IdxTy pcol=1;
const IdxTy points=100; // m_points;
const D fudge=m_flp.gr_xsection();
const D tau=m_flp.gr_tau();
const D coef=fudge/tau;

//if (li.bad()!=li.find("debug")) dbg=true;

MyBlock rsol,ndiff;
rsol.resize(points,2);
for (IdxTy i=0; i<points; i+=2)
{
rsol(i,pcol)=1e14*i*i/points;
rsol(i,ncol)=1e8;
rsol(i+1,ncol)=1e14*i*i/points;
rsol(i+1,pcol)=1e8;

}
for (IdxTy j=0; j<20; ++j)
{
ndiff=rsol;
D dt=1e-10;
D t=dt*j;
 gr_integrating(rsol,t,ncol,pcol,points);
//	if (integrate_gr) gr_integrating(m_state,t,in,ip);
for (IdxTy i=0; i<points; ++i) 
{
D ni2=rsol(i,ncol)*rsol(i,pcol)/m_ni2;
D ni2zed=ndiff(i,ncol)*ndiff(i,pcol)/m_ni2;

MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(t)<<MMPR(ni2zed)<<MMPR(ni2)<<MMPR(rsol(i,ncol))<<MMPR(rsol(i,pcol))<<MMPR(ndiff(i,ncol))<<MMPR(ndiff(i,pcol))<<MMPR(coef))

} // i
rsol=ndiff; 
} // j 
} // test_gr_integrator



void test_diff_impulse( CommandInterpretter & li )
{

bool dbg=false;
//bool dbgbc=false;
if (li.bad()!=li.find("debug")) dbg=true;
//if (li.bad()!=li.find("debugbc")) dbgbc=true;

MyBlock rsol,ndiff;
const IdxTy ncol=0;
const IdxTy vcol=1;
const IdxTy points=m_points;
rsol.resize(points,2);
for (IdxTy i=0; i<points; ++i)
{
rsol(i,ncol)=0;
if (i==(points/2)) rsol(i,ncol)=1;
//if ((1+i)==(points/2)) rsol(i,ncol)=1;
//if ((i-1)==(points/2)) rsol(i,ncol)=1;
rsol(i,vcol)=m_flp.vmax()/2;

}
const bool left_src=true;
const bool right_src=true;

if (left_src) rsol(0,ncol)=1;
if (right_src) rsol(points-1,ncol)=1;

D v=m_flp.vmax();
D dx=m_fixed(points/2,m_x)-m_fixed(points/2 -1, m_x);
// the mu is needed for DIFFUSION to create einstein factor but NOT for drift with v column 
//D dt=.1*dx/v/m_mun; // 1e-14;
// this created a flow of about 5 pixels per itor
//D dt=10.0*dx/v; // 1e-14;
D dt=.05*dx/v; // 1e-14;
const bool uniform_src=false;
const bool skip_self=false;
	const drift_flags df(skip_self,uniform_src,dbg);
	debug_levels dbl(dbg);
for (IdxTy j=0; j<1500; ++j)
{
D t=dt*(j+1);
//void do_sf_drift(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, const int iv)
//void fick_lanc_flow_sf(MyBlock & ndiff, const MyBlock & sol, const D & t, const int in, const D vfac=1.0)
ndiff.resize(points);
//drift_flags(): no_self(false), uniform_src(true), dbg(false) {}
	//const drift_flags df(false,true,false);
do_sf_drift( ndiff, rsol, dt, ncol, vcol,df, dbl);
//fick_lanc_flow_sf( ndiff, rsol,  dt, ncol,m_mun,dbg);
for (IdxTy i=0; i<points; ++i) rsol(i,ncol)=ndiff(i);
D total=0;
for (IdxTy i=0; i<points; ++i) 
{
total+=rsol(i,ncol);
D ival=rsol(i,ncol);
const D tar=0;
if ((ival-tar)*(ival-tar)>1e-58) { MM_MSG(" t= "<<t<<" i= "<<i<<" x= "<<m_fixed(i,m_x)<<" n= "<<rsol(i,ncol)<<" v= "<<rsol(i,vcol)<<" j= "<<j // <<MMPR(uniform_src) 
)  } 
} 
MM_MSG(MMPR(total))
} // j
//	fick_lanc_flow_sf( ndiff_fl, rsol,  t, m_n,m_mun);
// 	do_sf_drift(ndiff, rsol t, m_n, m_v_n);

}


void test_diffsf_bc( CommandInterpretter & li )
{

const StrTy id="test_driftfs_bc";
if (li.help())
{
// pass this back to li and then give up 
MM_ERR(" user invoked witha help command pass helpback to li")
return;

}

bool dbg=false;
bool dbgbc=false;
if (li.bad()!=li.find("debug")) dbg=true;
if (li.bad()!=li.find("debugbc")) dbgbc=true;

//const IdxTy scene=0; // check a few xmax for 4 bc sitatuiotions 
//const IdxTy scene=1; // make plot of one bc for xmax vaiation 
const IdxTy scene=2; // matrix of offsets and xmax
//const IdxTy scene=3; // DRIFT matrix of offsets and xmax

QFLE qflefs;
D n0=1e13; D n1=1e13; D n2=1e13;
D nleft=1e13;
D nright=1e13;
const IdxTy npix=10;
pixelt p[npix],x,z;
D n[npix];
//D del=.1;
// h is the disgtance between x_i and x_(i+1) so each pixel size is 2h
const D h=.5;
const D hbc=5;
x.set(0,h);
z.set(0,h);
// { del=0;  del0=.5; 
//{xmax=1e-4; dxmax=1e-4; niter=30; }
IdxTy i=0;
D xi=h*((i>>1)+1);
// x is the reference pixel 
x.set(0,h);
// right border overlaps 
p[0].set(xi,hbc);
// left border overlaps 
p[1].set(xi,-hbc);

//p[i+1].set(-xi,h);
MM_MSG(" reference pixel "<<"x"<<" "<<x.dump()<<" "<<x.dc(x))
MM_MSG(" right pixel "<<0<<" "<<p[0].dump()<<" "<<p[0].dc(x))
MM_MSG(" left pixel "<<1<<" "<<p[1].dump()<<" "<<p[1].dc(x))
MM_MSG(id<<MMPR(n0)<<MMPR(n1)<<MMPR(n2)<<MMPR(nright)<<MMPR(nleft))
//MM_MSG(" pixel "<<(i+1)<<" "<<p[i+1].dump()<<" "<<p[i+1].dc(x))
const bool no_self=false;
const bool no_selfs=!false;
// xp0+=del0; xp1+=del0;  xp2+=del0; 
pixelt & lbc=p[1];
pixelt & rbc=p[0];
D xmax=1e-4; D dxmax=.2; IdxTy niter=10; 
{ xmax=1e-4; dxmax=.2; niter=7;  }
if (0==scene) { // alt
for (IdxTy i=0; i<niter; ++i)
{
// these total numbers need to be turned into source density
const D bcfo=1.0/h/xmax/xmax;
const D bcf=1.0/h/xmax/xmax;
//D diff_sf_bc(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
///	, const D & n0, const D & n1, const D & n2,  const D & xmax, const bool no_self, const bool dbg) const
if (dbg ) {MM_MSG("----------------"<<i<<"     out_right")}
const D out_right=x.diff_sf(rbc,qflefs,n0,n1,n2,xmax, no_self, dbg);
//if (dbgbc) raise(SIGINT);
// probably working now 
const D out_right_bc=bcfo*x.diff_sf_bc(rbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
if (dbg ) {MM_MSG("----------------"<<i<<"     out_left") }
const D out_left=x.diff_sf(lbc,qflefs,n0,n1,n2,xmax, no_self, dbg);
const D out_left_bc=bcfo*x.diff_sf_bc(lbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
if (dbg ) {MM_MSG("----------------"<<i<<"     in_right") }
const D in_right=rbc.diff_sf(x,qflefs,nright,nright,nright,xmax, no_self, dbg);
const D in_right_bc=bcf*rbc.diff_sf_bc(x,qflefs,nright,nright,nright,xmax, no_self, dbgbc);
if (dbg ) {MM_MSG("----------------"<<i<<"     in_left") }
const D in_left=lbc.diff_sf(x,qflefs,nleft,nleft,nleft,xmax, no_self, dbg);
const D in_left_bc=bcf*lbc.diff_sf_bc(x,qflefs,nleft,nleft,nleft,xmax, no_self, dbgbc);
//MM_MSG(MMPR(xmax)<<MMPR(out_right)<<MMPR(in_right)<<MMPR(out_left)<<MMPR(in_left))
MM_MSG(MMPR(xmax)<<MMPR(bcf)<<MMPR(bcfo)<<MMPR(out_right_bc)<<MMPR(in_right_bc)<<MMPR(out_left_bc)<<MMPR(in_left_bc))
xmax+=dxmax;
}
} // alt
if (1==scene)
{
{ xmax=1e-4; dxmax=.01; niter=1000;  }
for (IdxTy i=0; i<niter; ++i)
{
const D bcfo=1.0/h/xmax/xmax;
const D bcf=1.0/h/xmax/xmax;
//const D out_right_bc=bcfo*x.diff_sf_bc(rbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
//const D in_right_bc=bcf*rbc.diff_sf_bc(x,qflefs,nright,nright,nright,xmax, no_self, dbgbc);
//const D out_left_bc=bcfo*x.diff_sf_bc(lbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
const D in_left_bc=bcf*lbc.diff_sf_bc(x,qflefs,nleft,nleft,nleft,xmax, no_self, dbgbc);
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(out_right_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(in_right_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(out_left_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(in_left_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
xmax+=dxmax;
}
}

if ((2==scene)||(3==scene))
{
drift_flags df;
D dxx=.03;
IdxTy miter=200;
z.move(-dxx*miter/2);
for (IdxTy j=0; j<miter; ++j)
{
{ xmax=1e-5; dxmax=.012*5; niter=300;  }
for (IdxTy i=0; i<niter; ++i)
{
const D bcfo=1.0/h/xmax/xmax;
const D bcf=1.0/h/xmax/xmax;
//const D out_right_bc=bcfo*x.diff_sf_bc(rbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
//const D in_right_bc=bcf*rbc.diff_sf_bc(x,qflefs,nright,nright,nright,xmax, no_self, dbgbc);
//const D out_left_bc=bcfo*x.diff_sf_bc(lbc,qflefs,n0,n1,n2,xmax, no_self, dbgbc);
const D in_left_bc=(scene==2)?(bcf*z.diff_sf_bc(x,qflefs,nleft,nleft,nleft,xmax, no_self, dbgbc)):0;
//const D drifter=(scene==3)?(z.drift_sf_bc(x,qflefs,nleft,nleft,nleft,xmax, no_self, dbgbc)):0;
const D drifter=(scene==3)?(z.drift_sf_bc(x,qflefs,nleft,nleft,nleft,xmax, df)):0;
const D drifterold=(scene==3)?(z.drift_sf(x,qflefs,nleft,nleft,nleft,xmax, df, dbgbc)):0;
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(out_right_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(in_right_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
//MM_MSG(MMPR(xmax)<<MMPR(bcfo)<<MMPR(out_left_bc)<<" @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ")
if (scene==2) { MM_MSG(MMPR(xmax)<<MMPR(z.x0)<<MMPR(in_left_bc)<<" scene2 ")}
if (scene==3) { MM_MSG(MMPR(xmax)<<MMPR(z.x0)<<MMPR(drifter)<<MMPR(drifterold)<<" scene3 ") }
xmax+=dxmax;
} //i 
z.move(dxx);
} // j
} // scene=2



}

void test_diffsf_conserve( CommandInterpretter & li )
{
const StrTy id="test_driftfs_conserve";
const IdxTy sz=li.size();
//D step =d(li,  sz, 2);
//int  steps  =n(li,  sz, 3);
//if ( sz>4 ) foo  =d(li,  sz, 4);
// fudding wrong fudder 
//QFLESF qflefs;
QFLE qflefs;
//D x0=0; D x1=.5; D x2=1;  
//D n0=1e6; D n1=1e13; D n2=1e13;
D n0=1e13; D n1=1e13; D n2=1e13;
//D xp0=0; D xp1=.5; D xp2=1; 
//D yp0=0; D yp1=.5; D yp2=1; 
const IdxTy npix=6;
pixelt p[npix],x;
D n[npix];
//D del=.1;
D xmax=.5;
D dxmax=1;
IdxTy niter=5;
//D del0=0;

const D h=.5;
x.set(0,h);
// { del=0;  del0=.5; 
{xmax=1e-4; dxmax=1e-4; niter=30; }
{xmax=1e-2; dxmax=1e-1; niter=30; }

// not symmetric but conservative over range 
{xmax=1e-5; dxmax=1e-4; niter=30; }
{xmax=1e-5; dxmax=.08; niter=30; }

MM_MSG(" carriers n0"<<n0<<" n1="<<n1<<" n2="<<n2)
for (IdxTy i=0; i<npix; i+=2)
{
const D xi=h*((i>>1)+1);
p[i].set(xi,h);
p[i+1].set(-xi,h);
MM_MSG(" pixel "<<i<<" "<<p[i].dump()<<" "<<p[i].dc(x))
MM_MSG(" pixel "<<(i+1)<<" "<<p[i+1].dump()<<" "<<p[i+1].dc(x))

}

const bool no_self=false;
const bool no_selfs=!false;
// xp0+=del0; xp1+=del0;  xp2+=del0; 
//del0=-.5;
// yp0+=del0; yp1+=del0;  yp2+=del0; 
//del=.2; xmax=2;
// obsolete in version 9 
const D neff=0; // qflefs.self_sf_loop( x.x0, x.x1, x.x2, n0, n1, n2,false); // const
for (IdxTy j=0; j<niter;  ++j)
{
D ntotal=0;
const D nself=0; // qflefs.diff_sf_loop(x.x0,x.x1,x.x2, x.x0, x.x1, x.x2, n0, n1, n2, xmax, no_selfs, false); // /xmax; 
MM_MSG(" .......   new del ..................")
for (IdxTy i=0; i<npix;  ++i)
{
// diff_sf(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,dbg);
//D nf=qflefs.diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
 MM_MSG(" ------------ i= "<<i<<" --------------------")
n[i]=0; // qflefs.diff_sf_loop(p[i].x0, p[i].x1, p[i].x2, x.x0, x.x1, x.x2, n0, n1, n2, xmax,no_self, (i<4)); 
ntotal+=n[i];
}
//MM_MSG(" .......  nr integral  ..................")
//D nf=qflefs.disjoint_diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
//D nr=qflefs.diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
//D nr=qflefs.diff_sf_loop(yp0, yp1, yp2, x0, x1, x2, n0, n1, n2, xmax, false); 
//D nr=qflefs.disjoint_diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
//const D dxi=x0-xp2;
//const D dxj=xp0-x2;
//MM_MSG(" .......  nr and nf results  ..................")
std::stringstream ss;
for (IdxTy i=0; i<npix;  ++i) ss<<" n["<<i<<"]="<<n[i];
const D t1=nself+ntotal;
const D f1=t1/neff;
const D f2=ntotal/nself;

// 0.00291 neff 1e+13 nself 1.93714e+10 ntotal 6.66187e+12 f1 0.668124 f2 343.902 1.00075 1.00219
// marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$
// echo test diffsfconserve | ./a.out 2> xxx  | grep h2076 | sed -e 's/=/ /g' | mjm zed 2 
// | awk '{print $0" "($11*$1)" "($9*1.5)}' | sort -g -k 12


MM_MSG(" xmax="<<xmax<<" neff="<<neff<<" nself="<<nself<<" ntotal="<<ntotal<<" diff="<<(neff-nself-ntotal)<<ss.str()); 
MM_MSG( MMPR(xmax)<<MMPR(neff)<<MMPR(nself)<<MMPR(ntotal)<<MMPR(f1)<<MMPR(f2))
//MM_MSG(" x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" xp0="<<xp0<<" xp1="<<xp1<<" xp2="<<xp2<<" xmax="<< xmax) 
//if ((nr+nf)!=0) 
//	{ MM_MSG("dxi= "<<dxi<<" dxj="<<dxj<<" xmax="<<xmax<<" nr="<<nr<<" nf="<<nf<<" diff="<<(nf-nr)<<" "<<((nf-nr)/(nf+nr)*.5))
xmax+=dxmax;
}

// xp0+=del; xp1+=del;  xp2+=del; 

} //i 





///////////////////////////////////////////////////////


void test_diffsf_mat( CommandInterpretter & li )
{
const StrTy id="test_driftfs";
const IdxTy sz=li.size();
//D step =d(li,  sz, 2);
//int  steps  =n(li,  sz, 3);
//if ( sz>4 ) foo  =d(li,  sz, 4);
// fudding wrong fudder 
//QFLESF qflefs;
QFLE qflefs;
D xmax=.5;
for (IdxTy xi=0; xi<50; ++xi)
{
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D xp0=0; D xp1=.5; D xp2=1; 
D del=.1;
IdxTy niter=5;
D del0=1;

 { del=.1; del0=1; xmax=.5; niter=5; }
 { del=.1; del0=1; xmax=2.0; niter=25; }
 { del=.09333; del0=1; xmax=3.14159; niter=80; }
 { del=.09333; del0=1; xmax=xi*.1+.1; niter=80; }

const bool no_self=false;
 xp0+=del0; xp1+=del0;  xp2+=del0; 
//del=.2; xmax=2;
for (IdxTy i=0; i<niter;  ++i)
{
MM_MSG(" .......   new del ..................")
// diff_sf(di.xp0(),di.xp1(),di.xp2(), di.x0(),di.x1(),di.x2(), di.n0(), di.n1(), di.n2(),xmax,dbg);
//D nf=qflefs.diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
D nf=0; // qflefs.diff_sf_loop(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, no_self, !false); 
MM_MSG(" .......  nr integral  ..................")
//D nf=qflefs.disjoint_diff_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, xmax, !false); 
//D nr=qflefs.diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
D nr=0; // qflefs.diff_sf_loop(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, no_self, !false); 
//D nr=qflefs.disjoint_diff_sf(x0, x1, x2, xp0, xp1, xp2, n0, n1, n2, xmax, !false); 
const D dxi=x0-xp2;
const D dxj=xp0-x2;
MM_MSG(" .......  nr and nf results  ..................")

MM_MSG(" x0="<<x0<<" x1="<<x1<<" x2="<<x2<<" xp0="<<xp0<<" xp1="<<xp1<<" xp2="<<xp2<<" xmax="<< xmax) 
if ((nr+nf)!=0) 
	{ MM_MSG("dxi= "<<dxi<<" dxj="<<dxj<<" xmax="<<xmax<<" nr="<<nr<<" nf="<<nf<<" diff="<<(nf-nr)<<" "<<((nf-nr)/(nf+nr)*.5))
}

 xp0+=del; xp1+=del;  xp2+=del; 

} //i 
} // xi 

} // test_diffsf_mat



void test_driftsf( CommandInterpretter & li )
{
const StrTy id="test_driftfs";
const IdxTy sz=li.size();
// fudding wrong fudder 
//QFLESF qflefs;
QFLE qflefs;
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D xp0=0; D xp1=.5; D xp2=1; 
D del=.1;
D vt=.5;

for (IdxTy i=0; i<20; ++i)
{
// version 9 
D n=0; // qflefs.drift_sf(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, vt, !false); 
MM_MSG(" i="<<i<<" xp0="<<xp0<<" x0="<<x0<<" vt="<<vt<<" n="<<n)
 xp0+=del; xp1+=del;  xp2+=del; 
}


 xp0=0; xp1=.5;  xp2=1; 
del=-.2;
vt=-.5;
const drift_flags df(false,false,!false);
for (IdxTy i=0; i<10; ++i)
{ // changed to bc in version 9 
//D n=qflefs.drift_sf_bc(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, vt,false, !false); 
D n=qflefs.drift_sf_bc(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, vt,df); 
MM_MSG(" i="<<i<<" xp0="<<xp0<<" x0="<<x0<<" vt="<<vt<<" n="<<n)
 xp0+=del; xp1+=del;  xp2+=del; 
}

 xp0=3; xp1=3.5;  xp2=4; 
del=.2;
vt=5;
for (IdxTy i=0; i<25; ++i)
{
//D n=qflefs.drift_sf_bc(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, vt, false, !false); 
D n=qflefs.drift_sf_bc(xp0, xp1, xp2, x0, x1, x2, n0, n1, n2, vt, df); 
MM_MSG(" i="<<i<<" xp0="<<xp0<<" x0="<<x0<<" vt="<<vt<<" n="<<n)
 xp0+=del; xp1+=del;  xp2+=del; 
}




} // test_driftsf


void test_flow( CommandInterpretter & li )
{
const StrTy id="test_flow";
const IdxTy sz=li.size();
QFLE qfle;
D x0=0; D x1=.5; D x2=1;  D n0=1; D n1=1; D n2=1;
D A,B,C;
qfle.fitq(A,B,C, x0, x1, x2, n0, n1, n2) ;
//ngvdx(di.a(),di.b(),di.x0(),di.x1(),di.x2(), di.n0(),di.n1(),di.n2(),xmax,vt);
D a=-1;
D b=0;
D vt=0;
D xmax=1;
if (sz<4)
{
// version 9 
D f=0; // qfle.ngvdx(a,b, x0, x1, x2, n0, n1, n2,xmax,vt,true) ;
MM_MSG(id<<" a="<<a<<" b="<<b<<" x0="<<x0<<" x2="<<x2<<" xmax="<<xmax<<" flow="<<f);
a=1; b=2;
// version 9 
 f=0; // qfle.ngvdx(a,b, x0, x1, x2, n0, n1, n2,xmax,vt,true) ;
MM_MSG(id<<" a="<<a<<" b="<<b<<" x0="<<x0<<" x2="<<x2<<" xmax="<<xmax<<" flow="<<f);
MM_MSG(" ish assumed to br 1")
return; 
}

a =li.d(2); // d(li,  sz, 2);
b =li.d(3); // d(li,  sz, 3);
if ( sz>4 ) xmax  =li.d(4); // d(li,  sz, 4);
 D f=0; // qfle.ngvdx(a,b, x0, x1, x2, n0, n1, n2,xmax,vt,true) ;
MM_MSG(id<<" a="<<a<<" b="<<b<<" x0="<<x0<<" x2="<<x2<<" xmax="<<xmax<<" flow="<<f);
MM_MSG(" ish assumed to br 1")

} // test_flow 

// the command object needs a shift function or something
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
if (cmd=="flow") { test_flow(li); li.pop(); return ; } 
if (cmd=="symm") { test_flow_symmetry(li,true); li.pop(); return ; } 
if (cmd=="symmplot") { test_flow_symmetry(li,!false,1.01); li.pop(); return ; } 
if (cmd=="driftsf") { test_driftsf(li); li.pop(); return ; } 
if (cmd=="diffsf") { test_diffsf(li); li.pop(); return ; } 
if (cmd=="diffsfmat") { test_diffsf_mat(li); li.pop(); return ; } 
if (cmd=="diffsfolap") { test_diffsf_olap(li); li.pop(); return ; } 

//void test_diffsf_conserve( CommandInterpretter & li )
if (cmd=="diffsfconserve") { test_diffsf_conserve(li); li.pop(); return ; } 
if (cmd=="diffsfbc") { test_diffsf_bc(li); li.pop(); return ; } 
if (cmd=="diffimpulse") { test_diff_impulse(li); li.pop(); return ; } 
if (cmd=="gr") { test_gr_integrator(li); li.pop(); return ; } 
//void test_gr_integrator( CommandInterpretter & li )
MM_ERR(" unrecignized TEST command "<<li.dump())

li.pop();
} // test


void steps(const IdxTy n)
{

IdxTy iter=0;
for (iter=0; iter<n; ++iter)
{
	m_cm.inc("iterations");	
	step(iter);
	per_step_outputs(iter);	

} // iter


}


void refine(const D & thr )
{
gerrymander gm;
MM_ERR(" attempting to refine with thr="<<thr)
//(MyBlock & newpoints, const MyBlock & fixed, const MyBlock & sol, const IdxTy ix, const IdxTy iy, const D & thr)
MyBlock np(3); // this should also take an std::vector 
gm.max_delta_refine(np,m_fixed,m_state,m_x,m_u,thr); 
MM_ERR(" attempt to reinit with pointsize = "<<np.size()) 
Init(0,np);
}


void config_banner()
{
MM_INC_MSG(m_cm,"test test" ) 
MM_MSG(" configuration "<<m_flp.to_string()<<" points "<<m_points)
MM_MSG(" logic "<<m_tree.to_string())
}

void solve()
{
const IdxTy miter= m_flp.max_iterations(); // 10000;
//#define MM_INC_MSG(cm,msg) cm.inc_location( __FILE__,__LINE__,msg);
config_banner();

IdxTy iter=0;
for (iter=0; iter<miter; ++iter)
{
	m_cm.inc("iterations");	
	step(iter );
	per_step_outputs(iter);	

} // iter

final_outputs(iter);



}
void final_outputs(const IdxTy iter )
{
	if (m_tree.dump_state_atend()) 
	{ 
		dump_state_etc(std::cout,"fick",0);
// 		const StrTy ival=(m_cm).format_entry("iterations");
//		const StrTy lbl="fick"+ival;
//		std::cout<<m_state.to_string_cols(lbl,2,8," "); 
	} 


}

void per_step_outputs(const IdxTy iter )
{
	const FlpTy & flp= m_flp;
	const IdxTy freq=flp.output_freq(); // 100;
	const IdxTy offset=flp.output_offset(); // freq-1;
	
	const bool everyfew= ((iter%freq)==offset);
 	const StrTy ival=(m_cm).format_entry("iterations");
	// the termina point had been forced but is not now
	// actually it still is but not by the shooter or solver
	if (m_tree.show_iter_status()) 
	{	MM_STATUS(" iteration "<<ival<<" "<<m_state(m_points-2+0,m_u)<<status_string()) }
	if (m_tree.always_dump_state())  dump_state_etc();
 	else if (m_tree.dump_state_everyfew()&&everyfew)  dump_state_etc();

	if (m_tree.dump_counter_everyfew()) if (everyfew)
	{
		m_cm.dump("solve_step"  , std::cout);
		MM_MSG("solve_times"<<CRLF<<m_cm.time_table("solver_times"))
	}
	if (m_tree.cerr_counter_everyfew()) if (everyfew)
	{
		m_cm.dump("solve_step"  , std::cerr);
	}

}
StrTy status_string() const
{
Ss ss;
D nimax=0;
D nimin=1e20*m_ni2;
	for (IdxTy i=0; i<m_points; ++i)
		{
			D np=m_state(i,m_n)*m_state(i,m_p)/m_ni2;
			if (np>nimax) nimax=np;
			if (np<nimin) nimin=np;
		}
ss<<MMPR(nimax);
ss<<MMPR(nimin);
ss<<"    "; // guard, although the status macro should erase the line doh
return ss.str();

}
//void dump_legend() { dump_state_etc(std::cout,"fick",1); } 
template <class Os> void dump_legend(Os & os, const StrTy & label="legend",const IdxTy flags=0)
{

// the extraction code must include the iter info to R
//  and it needs to be called at the tight time or the iter count is fudded. 
	const IdxTy ival=(m_cm).get_count("iterations");
	const StrTy lbl=" Rlegend ";
	// must keep label
	os<<lbl<<" label "<<label<<" "<<ival<<CRLF;

}


void dump_to_file(const StrTy & fn , const StrTy & label, const IdxTy flags=0)
{
std::ofstream fout(fn);
MM_ERR(" dumping to "<<fn);
// dump_state_etc(fout,label,flags);
dump_state_R(fout,label,flags);

}

void dump_state_etc() { dump_state_etc(std::cout,"fick",1); } 

template <class Os> void dump_state_R(Os & os, const StrTy & label="",const IdxTy flags=0)
{
	const bool column_names=true;
	const IdxTy ival=(m_cm).get_count("iterations");
		MyBlock qoi;
		compute_qoi(qoi);

	if (column_names) { os<<"iters"<<""<< " "<< m_state.col_names()<<" "<<qoi.col_names()<<CRLF; }
	os<<m_state.to_string_cols_R(ival,12,12,qoi," ",false);   

}


template <class Os> void dump_state_etc(Os & os, const StrTy & label="fick",const IdxTy flags=0)
{
	//const StrTy lbl="fick"+ival;
	const StrTy ival=(m_cm).format_entry("iterations");
	const StrTy lbl=label+ival;
	const bool human_readable=(0==(flags&1));
	const bool column_names=(1==(flags&2));
	switch (flags)
	{
		case 0: 
		case 1: 
		default:
		MyBlock qoi;
		compute_qoi(qoi);
		if (column_names) { os<<"cols"<<lbl<< " "<< m_state.col_names()<<" "<<qoi.col_names()<<CRLF; }
		if (human_readable)  { os<<m_state.to_string_cols(lbl,2,8,qoi," ");  }
		// machine readable more precisions 
		else { os<<m_state.to_string_cols(lbl,12,12,qoi," ",false);  } 

	} // switch 
}

// move this down into real code 
void compute_qoi( MyBlock & qoi)
{
	const IdxTy nqoi=2;
	qoi.resize(m_points,nqoi);
	qoi.col_name("current",0);
	qoi.col_name("x",1);

	const IdxTy bc=1;
	if (m_points<bc) return; 
	const MyBlock & s=m_state;
	const IdxTy df=1; const IdxTy db=1;
	const IdxTy sz=m_points-bc;
	for (IdxTy i=bc; i<sz; ++i)
	{
		const D ih=1.0/(m_fixed(i+df,m_x)-m_fixed(i-db,m_x));
		const D gradU=ih*(s(i+df,m_u)-s(i-db,m_u));
		const D gradP=ih*(s(i+df,m_p)-s(i-db,m_p));
		const D gradN=ih*(s(i+df,m_n)-s(i-db,m_n));
		// this needs to use the calculated velocitiies 
		//const D phi_n=m_mun*s(i,m_n)*(gradU) -m_Vt*m_mun*gradN;
		const D phi_n=s(i,m_v_n)*s(i,m_n) -m_Vt*m_mun*gradN;
		//const D phi_p=-m_mup*s(i,m_p)*(gradU)-m_Vt*m_mup*gradP;
		const D phi_p=s(i,m_v_p)*s(i,m_p)-m_Vt*m_mup*gradP;
		// this is charge current, the phis are particle current
		const D jc=phi_n-phi_p;
		qoi(i,0)=jc;
		qoi(i,1)=m_fixed(i,m_x); // no easy way to include location lol 
	} // i 
	// kluge
	qoi(0,0)=qoi(1,0);
	qoi(sz,0)=qoi(sz-1,0);
	qoi(0,1)=m_fixed(0,m_x); // no easy way to include location lol 
	qoi(sz,1)=m_fixed(sz,m_x); // no easy way to include location lol 



} // compute qoi 

void step(const IdxTy iter ) // const step_param & sp)
{
const FlpTy & flp= m_flp;
const D t=flp.fl_t(); //1e-11; // sp.dtau;
// this should not be used now 
//const D dxmax=flp.fl_dmax(); // m_h*35; // sp.dxmax;
const bool shoot_voltage=m_tree.shoot_voltage(); // !false; // probably faster and could be improved
const bool shoot_itor_voltage=m_tree.shoot_itor_voltage(); // !false; // probably faster and could be improved
const bool integrate_gr=m_tree.integrate_gr(); // !false; // probably faster and could be improved
const bool integrate_sf_gr=m_tree.integrate_sf_gr(); // !false; // probably faster and could be improved
const bool enable_gr=m_tree.enable_gr(); // !false; // probably faster and could be improved
const bool enable_fl_flow=m_tree.enable_fl_flow(); // !false; // probably faster and could be improved
const bool enable_subpixel_fl_flow=m_tree.enable_subpixel_fl_flow(); // !false; // probably faster and could be improved
const bool enable_sf_fl_flow=m_tree.enable_sf_fl_flow(); // !false; // probably faster and could be improved
const bool enable_drifting=m_tree.enable_drifting(); // !false; // probably faster and could be improved
const bool enable_sub_drift=m_tree.enable_sub_drift(); // !false; // probably faster and could be improved
const bool enable_sub_drift_itor=m_tree.enable_sub_drift_itor(); // !false; // probably faster and could be improved
const bool enable_sf_drift=m_tree.enable_sf_drift(); // !false; // probably faster and could be improved
const bool def_drift=!enable_sub_drift&&!enable_sub_drift_itor&&!enable_sf_drift;
const bool sanity_check=m_tree.check_for_negs();
const bool defer_dd_updates=m_tree.defer_dd_updates();
const bool output_after_drift=m_tree.output_after_drift();
const bool enable_trapping=m_tree.enable_trapping();
//const IdxTy pmax=m_points-0;
//for (IdxTy point=0; point<pmax; ++point)
{
	m_cm.mark("step");
	// examine field to come up with velocity distributions
	// as drift shifts this at each point
	m_cm.mark("find_pot");
	// if the shotting fixed the initial e field it would be ok 
	if (shoot_itor_voltage)	drifts_itor(m_state);
	else if (shoot_voltage)	{MM_ONCE(" old drift effecitvely disabled",) } // drifts(m_state);
	else	drifts_petsc(m_state);
	m_cm.cumdel("pottotal","find_pot");
	// move all the carriers at this point to their respective desinations 
	// eventuallly tally those going into boundaries
	MyBlock ndiff_fl,pdiff_fl;
	if( enable_fl_flow) {	
	m_cm.mark("fick_lanc");

	enforce_bc(m_state); // this is now needed as bc carriers diffuse back on above flow  
	//fick_lanc_flow( ndiff, m_state,  t,  dxmax, m_n, m_v_n);
	//fick_lanc_flow( pdiff, m_state,  t,  dxmax, m_p, m_v_p);

	if(! enable_subpixel_fl_flow&&!enable_sf_fl_flow) {	
	fick_lanc_flow_no_drift( ndiff_fl, m_state,  t, m_n);
	fick_lanc_flow_no_drift( pdiff_fl, m_state,  t, m_p);
}
	if( enable_subpixel_fl_flow) {	
	MM_ONCE(" fick_lanc_flow_subpixel left in v 8",)
//	fick_lanc_flow_subpixel( ndiff_fl, m_state,  t, m_n,m_mun);
//	fick_lanc_flow_subpixel( pdiff_fl, m_state,  t, m_p,m_mup);
}

	if( enable_sf_fl_flow) {	
	fick_lanc_flow_sf( ndiff_fl, m_state,  t, m_n,m_mun);
	fick_lanc_flow_sf( pdiff_fl, m_state,  t, m_p,m_mup);
}

	// need to assign columns in m_state from ndiff
	if (!defer_dd_updates)
	{		
		m_state.assign_column(pdiff_fl,m_p);
		m_state.assign_column(ndiff_fl,m_n);
	}
	m_cm.cumdel("fltotal","fick_lanc");

	if (sanity_check){  check_carriers("flow")	; MM_MSG(" did FL sanity check ") } 
	enforce_bc(m_state); // this is now needed as bc carriers diffuse back on above flow  
	}

if (enable_drifting) { 
	MyBlock ndiff,pdiff;
	const drift_flags df(false,!true,false);
	if (!df.uniform_src) { MM_ONCE( "uniform source DISABLED",) }
	const debug_levels dbl(false);
	m_cm.mark("drifts");
	const IdxTy val=arbor("drifting mode",def_drift,enable_sub_drift,enable_sub_drift_itor);
//	do_drift(ndiff, m_state, t, dxmax, m_n, m_v_n);
//	do_drift(pdiff, m_state, t, dxmax, m_p, m_v_p);
	if (def_drift) do_drift_one(ndiff, m_state, t,  m_n, m_v_n);
// 	else if (enable_sub_drift) do_drift_subpixel(ndiff, m_state, t, m_n, m_v_n);
// 	else if (enable_sub_drift_itor) do_drift_subpixel_itor(ndiff, m_state, t, m_n, m_v_n);
 	else if (enable_sf_drift) do_sf_drift(ndiff, m_state, t, m_n, m_v_n,df,dbl);

	if (def_drift) do_drift_one(pdiff, m_state, t,  m_p, m_v_p);
//	else if (enable_sub_drift) do_drift_subpixel(pdiff, m_state, t, m_p, m_v_p);
//	else if (enable_sub_drift_itor) do_drift_subpixel_itor(pdiff, m_state, t, m_p, m_v_p);
	else if (enable_sf_drift) do_sf_drift(pdiff, m_state, t, m_p, m_v_p,df,dbl);

	if (defer_dd_updates)
{
	m_state.assign_column((pdiff+pdiff_fl)*.5,m_p);
	m_state.assign_column((ndiff+ndiff_fl)*.5,m_n);
}
else
{
if (false) MM_MSG(" pdiff(499)="<<pdiff(499)<<" ndiff(499)="<<ndiff(499))
	m_state.assign_column(pdiff,m_p);
	m_state.assign_column(ndiff,m_n);
}
	m_cm.cumdel("drifttotal","drifts");
}


	if (sanity_check ) check_carriers("drift")	;

	enforce_bc(m_state); // ???? 
	// apply GR to help approach equillibrium 
	if (output_after_drift) per_step_outputs(iter);

	if (enable_gr) { 
	m_cm.mark("find_gr");
	if (integrate_gr&&!integrate_sf_gr) 
		gr_integrating(m_state,t,m_n,m_p,m_points);
	else 
		if (integrate_sf_gr)  gr_integrating_sf(m_state, t, m_n,m_p, m_points);
	else { MM_ONCE(" NOTE gr effecitvely disabled", ) }// gr(m_state,t);
	m_cm.cumdel("grtotal","find_gr");
	}

	if (sanity_check) check_carriers("gr")	;

	enforce_bc(m_state); // ????? 

if (enable_trapping)
{
	m_cm.mark("find_traps");
//void step( const MyBlock & fixed, const MyBlock & sol, const D & dti, const IdxTy & points, const IdxTy & in, const IdxTy &ip, const IdxTy & iu, const IdxTy & ix)

	m_defects.step(m_fixed,m_state,t,m_points,m_n,m_p,m_u,m_x);
	m_defects.net_charge(m_state,m_trap,m_points);


	m_cm.cumdel("traptotal","find_traps");
} // enable_trapping


	m_cm.cumdel("steptotal","step");

} // point


} // step

// this is probably not really needed or makes confusing results 
void enforce_bc(MyBlock & sol)
{
	const int iu=m_u; const int in=m_n; const int ip=m_p;
	const IdxTy sz=m_points-1;
// the solver is the only one to change these so it should be ok 
//	sol(0,iu)= m_v0; 
//	sol(sz,iu)= m_vd;
 	sol(0,in)= m_n_p; 
	sol(0,ip)=m_p_p;  
	sol(sz,in)= m_n_n; 
	sol(sz,ip)=m_p_n; 
//	sol(point,iz)=1;


}

// needs nancheck
IdxTy count_negs(MyBlock & state, const IdxTy col)
{
	//typedef mjm_numb Ck;
	IdxTy k=0; 
	const IdxTy sz=m_points;
	for (IdxTy i=0; i<sz; ++i) {
		 if ( state(i,col)<0) {++k;
		MM_MSG(" negative carrier "<<i<<" col "<<col<<" value "<<state(i,col))

} } // for


return k;

}

// make this a member of MyBlock and do something with FILE and LINE lol 
IdxTy count_nans(MyBlock & state, const IdxTy col)
{
	typedef mjm_numb Ck;
	IdxTy k=0; 
	const IdxTy sz=m_points;
	for (IdxTy i=0; i<sz; ++i) 
	{ 
		if (Ck::denan(state(i,col) ,__FILE__,__LINE__," fudd "))
		{ ++k; 
			MM_MSG(" nan found  for carrier  "<<i<<" "<<col<<" "<<state(i,col)<<" "<<"")
		}
	}

return k;
}
bool check_carriers(const StrTy & label  )
{
const IdxTy nneg=count_negs(m_state,m_n);
const IdxTy npos=count_negs(m_state,m_p);
const IdxTy nnan=count_nans(m_state,m_n);
const IdxTy npan=count_nans(m_state,m_p);

const bool bad=((nneg+npos)>0);
const bool badder=((nnan+npan)>0);
if (bad||badder)
{ 
MM_MSG(label<<" neagative carrier counts n="<<nneg<<" "<<npos)
MM_ERR(label<<" neagative carrier counts n="<<nneg<<" "<<npos)
m_cm.inc(label+"nneg",nneg);
m_cm.inc(label+"npos",npos);
MM_MSG(label<<" nan carrier counts n="<<nnan<<" "<<npan)
MM_ERR(label<<" nan carrier counts n="<<nnan<<" "<<npan)
m_cm.inc(label+"nneg",nnan);
m_cm.inc(label+"npos",npan);

}
return bad||badder;
}

/////////////////////////////////////////////////////////
// start of real code that does things 
////////////////////////////////////////////////////////

// probably just integrate this, it should not cross equ or collect 200 dollars lol 
// can be done as a function of position however physics becomes
// a question- maybe it smears out anyway and nothing is added 

void gr_integrating_sf(MyBlock & sol, const D & dt, const IdxTy in, const IdxTy ip, const IdxTy points)
{

	const FlpTy & flp= m_flp;
	const D fudge=flp.gr_xsection(); // 1e-3;
	const bool debug_gr=true;
	const D tau=flp.gr_tau();
	const D phitau=fudge/tau;
	Di di(points,m_x,in,~0,sol,m_fixed);
	Di dip(points,m_x,ip,~0,sol,m_fixed);
	MyBlock diff(points);
	const D f1=.5;
//	const drift_flags df(false,true,false);
	while (di.ok())
	{
		const IdxTy i=di.i();
		const D ish=1.0/di.sh();
//		const D neff=2.0*f1*qfle.isselflin(di)*ish;
		const D guard=di.sh();
		const D x0=di.x0(); 
		const D x1=di.x1(); 
		const D x2=di.x2();
		const D n0=di.n0(); 
		const D n1=di.n1(); 
		const D n2=di.n2();
		const D p0=dip.n0(); 
		const D p1=dip.n1(); 
		const D p2=dip.n2();
		D dnp=0;
		// numerically integrate with or without sf 
		const IdxTy npts=23;
		// for uniform spacing WTF
		const D dh=1.0/npts;
		for (IdxTy j=0; j<npts; ++j)
		{
			const D f=1.0*j/(npts-1);
			const D x=(x2-x0)*f+x0;
			D n=0;
			D p=0;
			if ( x>x1) 
			{ D fr=(x-x1)/(x2-x1); n=n1*(1.0-fr)+n2*fr; 
			 p=p1*(1.0-fr)+p2*fr; 
				}
			else { D fr=(x-x0)/(x1-x0); n=n0*(1.0-fr)+n1*fr; 
			 p=p0*(1.0-fr)+p1*fr; 
				}
			D dn=gr_b2b_point(n,p,m_ni2,phitau,dt);
			if ( true)
			{
				if (x>x1) dn=dn*(x2-x)/(x2-x1);
				else dn=dn*(x-x0)/(x1-x0);
			}
			dnp+=dn;
		}
		dnp=dnp*dh;
		if (-dnp>n1) dnp=-.99999*n1;
		if (-dnp>p1) dnp=-.99999*p1;
		diff(i)=dnp;
		// this needs to update a temp and then copy back like all the others
		di.inc(); 
		dip.inc(); 

	} // di.ok(0
	const D ni=::sqrt(m_ni2);
		for (IdxTy j=1; j<(points-1); ++j) 
		{ sol(j,in)+=diff(j); sol(j,ip)+=diff(j); 
	//		if (sol(j,
	} 
}

// this is a copy of the code in the point wise gr method
// which can be removed once this seems to work and
// oor handling is fixed 
D gr_b2b_point(const D & n , const D & p, const D & ni2
	, const D & phitau, const D & t)
{
D dn=0;
const D expmax=300;
const D expmin=-300;
IdxTy oflow=0;
IdxTy uflow=0;
	const bool pgtn=(p>n);
	const D np=n*p/ni2;
	const D gr=phitau*(ni2-n*p);
	const D p0=pgtn?(p-n):(n-p);
	const D npinf=pgtn?n:p;
	const D disc=.5*::sqrt(p0*p0+4*ni2);
	const D a=-.5*p0-disc; // <0 
	const D b=-.5*p0+disc; // >0
	const D z=1.0/(b-a); const D w=-z;
	const D na=npinf-a;
	const D nb=npinf-b;
	if ((na*na)<(nb*nb))
	{ 
		const D A=(npinf-a)/(npinf-b);
		D tw=-t/w*phitau;
		const bool ovf=(tw>expmax);
		if (ovf){ ++oflow;  tw=expmax;}
		if (tw<expmin){ ++uflow; }
		D ex=(tw<expmin)?0:exp(tw);
		D gri=(a-b*A*ex)/(1-A*ex);
		if (ovf) { gri=b; }  // added 2017-04-11
		dn=gri-npinf; // gr*t;
		return dn;
	}
	if ((na*na)>(nb*nb))
	{
		const D Ainv=(npinf-b)/(npinf-a);
		D tw=-t/w*phitau;
		const bool ovf=(tw>expmax);
		if (ovf){ ++oflow;  tw=expmax;}
		if (tw<expmin){ ++uflow; }
		const D ex=(tw<expmin)?0:exp(tw);
		D gri=(a*Ainv-b*ex)/(Ainv-ex);
		if (ovf) { gri=b; } 
		dn=gri-npinf; // gr*t;
		return dn;
	} 
	const D ni=::sqrt(m_ni2); // save this ? lol 
	D dp=0;	
	if (p!=ni)
	{
		const D c=.5*::log((ni+p)/fabs(ni-p));
		D e=ni*phitau*t+c;
		const bool ovf=(e>expmax);
		if (ovf) { return ni; } 
		if (e>expmax){ ++oflow;  e=expmax;}
		if (e<expmin){ ++uflow; }
		const D A=::exp(e);
		D dp=ni*(A-1)/(A+1);
		return dp;
	}

return dn;
}

void gr_integrating(MyBlock & state, const D & t, const IdxTy in, const IdxTy ip, const IdxTy points)
{
	const FlpTy & flp= m_flp;
	const D fudge=flp.gr_xsection(); // 1e-3;
	const bool debug_gr=true;
	const D tau=flp.gr_tau();
	const D expmax=300;
	const D expmin=-300;
	//const D fudge=1e-13;
	// in places sz is points MINUS 1.... 
	const IdxTy sz=points; // m_points;
	IdxTy oflow=0;
	IdxTy uflow=0;

	for (IdxTy i=0; i<sz; ++i)
	{
	//	D n=state(i,m_n);
		D n=state(i,in);
		//D p=state(i,m_p);
		D p=state(i,ip);
		const bool pgtn=(p>n);
		const D np=n*p/m_ni2;
		// this is both dn and dp 
		const D gr=fudge*(m_ni2-n*p)/tau;
		// write p=p0+n and integrate by partial fractions with 
		// roots -p0/2  +/- .5*sqrt(po*po+4ni*ni)
		// note that po is NOT positivie definite
		// this only seems to work for p>n
		const D p0=pgtn?(p-n):(n-p);
		const D npinf=pgtn?n:p;
		const D disc=.5*::sqrt(p0*p0+4*m_ni2);
		const D a=-.5*p0-disc; // <0 
		const D b=-.5*p0+disc; // >0
		const D z=1.0/(b-a); const D w=-z;
		// w/(n-a) + z/(n-b) integrates to w*ln(n-a) + z*ln(n-b)a=t+C
		// exp(t+C)=((n-a)^w)*(n-b)^z
		// doh exp(t+C)=((n-a)/(n-b))^w
		// exp(t/w)=(n-a)/(n-b), (n-b)exp(t/w)=(n-a), n(1-Aexp(t/w))=a-b*Aexp(t/w) 
		// n(0)=(a-b*A)/(1-A), n(0)(1-A)=a-b*A, n(0)-a=A(n(0)-b), A=(n(0)-a)/(n(0)-b)
		//  if pn==ni^2, this is a constant and ex falls out but A can be zero or infinity.
		// pick the better A but putting smaller one on top 
		const D na=npinf-a;
		const D nb=npinf-b;
		if ((na*na)<(nb*nb))
		{ 
			const D A=(npinf-a)/(npinf-b);
			D tw=-t/w*fudge/tau;
			const bool ovf=(tw>expmax);
			if (ovf){ ++oflow;  tw=expmax;}
			if (tw<expmin){ ++uflow; }
			D ex=(tw<expmin)?0:exp(tw);
			D gri=(a-b*A*ex)/(1-A*ex);
			if (ovf) { gri=b; }  // added 2017-04-11
			D dn=gri-npinf; // gr*t;
	//	MM_ERR(" final dn is "<<dn)
			//state(i,m_n)+=dn;
			state(i,in)+=dn;
			//state(i,m_p)+=dn;
			state(i,ip)+=dn;

			if (debug_gr)
			{	
				D npfinal=state(i,in)*state(i,ip)/m_ni2;
				bool bad_gr=false;
				bad_gr|= ((np>1)&&(npfinal>=np));
				bad_gr|= ((np<1)&&(npfinal<=np));
				bad_gr|= ((np<1)&&(npfinal>=1));
				bad_gr|= ((np>1)&&(npfinal<=1));
				bad_gr|= ((np==1)&&(npfinal!=1));
				if (bad_gr)
				{
					MM_INC_MSG(m_cm,"grfail" ) 
					MM_MSG(" gr_err bigger "<<MMPR(np)<<MMPR(npfinal)<<MMPR(n)<<MMPR(p)<<MMPR(dn))
				}
			}
		
		}
		// strict equality can occur if p==n==ni2. 
		else if ((na*na)>(nb*nb))
		{
			const D Ainv=(npinf-b)/(npinf-a);
			D tw=-t/w*fudge/tau;
			const bool ovf=(tw>expmax);
			if (ovf){ ++oflow;  tw=expmax;}
			if (tw<expmin){ ++uflow; }
			const D ex=(tw<expmin)?0:exp(tw);
			D gri=(a*Ainv-b*ex)/(Ainv-ex);
			if (ovf) { gri=b; } 
			D dn=gri-npinf; // gr*t;
	//	MM_ERR(" final dn is "<<dn)
			//state(i,m_n)+=dn;
			state(i,in)+=dn;
			//state(i,m_p)+=dn;
			state(i,ip)+=dn;

		} else
		{ // |a| == |b|  and if disc is zero, 
		// but it still integrates ok
			const D ni=::sqrt(m_ni2); // save this ? lol 
			//t*(fudge/tau)*ni2+c=.5*ln((ni+p)/(|ni-p|)) = tfn+c
			// exp(tfn)=sqrt(ni+p)/sqrt(|ni-p|)
			// exp(2*tfn)*(ni-p)=(ni+p); A=exp(2*tfn)
			// ni(A-1)=p(A+1); p=ni(A-1)/(A+1) 
			// A=Cexp(ni*fudge/tau*t) , C=.5*ln((ni+p0)/(|ni-p0|))
			D dp=0;	
			if (p!=ni)
			{
			const D c=.5*::log((ni+p)/fabs(ni-p));
			D e=ni*fudge/tau*t+c;
			if (e>expmax){ ++oflow;  e=expmax;}
			if (e<expmin){ ++uflow; }
			const D A=::exp(e);
			D dp=ni*(A-1)/(A+1);
			//state(i,m_n)+=dp;
			state(i,in)+=dp;
			//state(i,m_p)+=dp;
			state(i,ip)+=dp;

			if (debug_gr)
			{	
				D npfinal=state(i,in)*state(i,ip)/m_ni2;
				bool bad_gr=false;
				bad_gr|= ((np>1)&&(npfinal>=np));
				bad_gr|= ((np<1)&&(npfinal<=np));
				bad_gr|= ((np<1)&&(npfinal>=1));
				bad_gr|= ((np>1)&&(npfinal<=1));
				bad_gr|= ((np==1)&&(npfinal!=1));
				if (bad_gr)
				{
					MM_MSG(" gr_err bigger "<<MMPR(np)<<MMPR(npfinal)<<MMPR(n)<<MMPR(p)<<MMPR(dp))
					MM_INC_MSG(m_cm,"grfail" ) 
				}
			}
		



			}
			MM_MSG(" equality a "<<a<<" b "<<b<<" p0 "<<p0<<" n "<<n<<" p "<<p<<" dp="<<dp)
			MM_ERR(" equality a "<<a<<" b "<<b<<" p0 "<<p0<<" n "<<n<<" p "<<p<<" dp="<<dp)
			MM_INC_MSG(m_cm,"grint_equality" ) 
		} 



	}
if (( uflow+oflow)!=0)
{
 	MM_MSG(" gr_integrate oflow "<<oflow<<" uflow "<<uflow)
	// this is slow, only when needed 
	m_cm.inc("grint_uflow",uflow);
	m_cm.inc("grint_oflow",oflow);

}


}

#if 0 
// use the diffuse_iterator and interpolated charges
// in progress, not working not used. 
 void shoot_voltage_q(MyBlock & state, const IdxTy sz, const bool fwd)
{
	//diffuse_iterator(const IdxTy points, const IdxTy _x, const IdxTy _in, const MyBlock & sol, const MyBlock & f)
	Di din(m_points,m_x,m_n,m_state,m_fixed);
	Di dip(m_points,m_x,m_p,m_state,m_fixed);

	QFLE qfle;
	while (din.ok())
	{
		D A,B,C;
		qfle.fitq(A,B,C,din);
 		D dvn=qfle.ngvdxi2(din.x2(),A,B,C) -qfle.ngvdxi2(din.x0(),A,B,C) ;
		qfle.fitq(A,B,C,dip);
 		D dvp=qfle.ngvdxi2(dip.x2(),A,B,C) -qfle.ngvdxi2(dip.x0(),A,B,C) ;
		// need to do something with fixed charge...ubpixel

		din.inc();
		dip.inc();
	}

}
#endif 

 void shoot_voltage(MyBlock & state, const IdxTy sz, const bool fwd)
{

	int st=1;
	int end=sz;
	int del=1;
	if (!fwd) { end=1; st=sz-1; del=-1; } 
	//for (IdxTy i=1; i<sz;  ++i)
	for (int  i=st; i<end;  i+=del)
	{
		// back to shooting from initial conditions, need lancozs here too lol
		const D x=m_fixed(i,m_x);
		const D xrev=m_fixed(i-1,m_x);
		const D xfwd=m_fixed(i+1,m_x);
		const D hf=xfwd-x;
		const D hr=x-xrev;
		const D hs=hr+hf;	
		const D rho=state(i,m_p)-state(i,m_n)+rho_fixed(i);
		if (fwd)
		{
			const D rhod=(-m_qe*rho*hr*hf*.5*hs+state(i,m_u)*hs-hf*state(i-1,m_u))/hr;
			// this would be easier with a linear solver but we have boundary and charge distro, 
			// only the e fields can be fixed.  
			state(i+1,m_u)=rhod;
		}
		else
		{
			const D rhod=(-m_qe*rho*hr*hf*.5*hs+state(i,m_u)*hs-hr*state(i+1,m_u))/hf;
			// this would be easier with a linear solver but we have boundary and charge distro, 
			// only the e fields can be fixed.  
			state(i-1,m_u)=rhod;
		}

	} // i 

} // shoot_voltage




void drifts_itor(MyBlock & state)
{
	const FlpTy & flp= m_flp;
	const bool do_rev_too=m_tree.do_reverse_shoot();	
	shooting_iterator si,sib;
	si.params(m_vd,flp.shooting_vi_min(), flp.shooting_vi_max()
		,flp.shooting_vd_tol_max(),flp.shooting_vd_tol_min()
		,flp.shooting_maxiter() );
	si.cm(&m_cm);
	const IdxTy sz=m_points-1;
	bool done=false;
	while (!done)
	{
		state(0,m_u)=m_v0; // should be safe to pull out of loop 
		state(1,m_u)=si.guess();
		shoot_voltage(state,sz,true);		
		done=si.result(state(sz,m_u));
	} // while
	if (do_rev_too)
	{
	sib.params(m_v0,m_vd+flp.shooting_vi_min(), m_vd+flp.shooting_vi_max()
		,flp.shooting_vd_tol_max(),flp.shooting_vd_tol_min()
		,flp.shooting_maxiter() );
	MyBlock temp(m_points);
	state.copy_column(temp,m_u);
	done=false;
	while (!done)
	{
		state(sz,m_u)=m_vd; // should be safe to pull out of loop 
		state(sz-1,m_u)=sib.guess();
		shoot_voltage(state,sz,!true);		
		done=sib.result(state(0,m_u));
	} // while
	state.axby_column(temp,m_u,.5,.5);
	} 	
//	m_cm.inc("shots",si.iterations());
//	m_cm.inc("backshots",sib.iterations());
	update_drift(state);

} // drifts_itor

// This seems to create noitcable prolems in the integration direction, need the 
// linear solver . Actually doing this in both directions and changing e field 
/// probably works better and faster. 
// it looks good on the left, better than the solver, but noisy on right
// need bidirectional etc. 

// this really is no better and probably a lot slower. 
void drifts_petsc(MyBlock & state)
{
	const FlpTy & flp= m_flp;
    const IdxTy iu=m_u;
	typedef mjm_petsc_solve Solver;
	MySparse lhs;	
	MyBlock rhs;
	rhs.resize(m_points);
	// setup rhs and lhs
    const IdxTy sz=m_points-1;
    for (IdxTy i=1; i<sz;  ++i)
{
	        const D x=m_fixed(i,m_x);
        const D xrev=m_fixed(i-1,m_x);
        const D xfwd=m_fixed(i+1,m_x);
        const D hf=xfwd-x;
        const D hr=x-xrev;
        const D hs=hr+hf;
		const D h2=2.0/(hs);
        const D rho=state(i,m_p)-state(i,m_n)+rho_fixed(i);
//        const D rhod=(-m_qe*rho*hr*hf*.5*hs+state(i,m_u)*hs-hf*state(i-1,m_u))/hr;
        // this would be easier with a linear solver but we have boundary and charge distro, 
        // only the e fields can be fixed.  
		lhs(i,i-1)=h2/hf; lhs(i,i)=-hs*h2/(hr*hf); lhs(i,i+1)=h2/hr;
		rhs(i)=-rho*m_qe;
	} // for 
	lhs(0,0)=1; rhs(0)=m_v0;
	lhs(sz,sz)=1; rhs(sz)=m_vd;

	Solver solver;
	const D rtol=flp.drift_solver_rtol(); // PETSC_DEFAULT;
	const D abstol=flp.drift_solver_abstol(); // PETSC_DEFAULT;
	const D dtol=flp.drift_solver_dtol(); // PETSC_DEFAULT;
//	const D abstol=PETSC_DEFAULT;
//	const D dtol=1e40;
	const D maxits=flp.drift_solver_maxits(); // 5000;

	solver.set(lhs,rhs,rtol,abstol,dtol,maxits);
	int reason=solver.solve();
	 solver.get(rhs);
	//if (reson>=0)
	m_state.assign_column(rhs,m_u);

	update_drift(state);
}

// split up 
void update_drift(MyBlock & state)
{
    // this is bad with zero points, fwiw
//    state(0,m_u)=m_v0; state(sz,m_u)=m_vd; state(1,m_u)=m_v0;
    const D vsat=m_flp.vsat(); // 1e5;
    const IdxTy iu=m_u;
	const IdxTy sz=m_points-1;
	const bool do_2_vel=true;
	if (do_2_vel) { MM_ONCE(" doing split velocity ",) } 
	IdxTy iinit=1;
	if (do_2_vel ) iinit=0;
    for (IdxTy i=iinit; i<sz;  ++i)
    {
        // back to shooting from initial conditions, need lancozs here too lol
	        const D x=m_fixed(i,m_x);
        const D xfwd=m_fixed(i+1,m_x);
        const D hf=xfwd-x;
// in the do_2_vel case, this is fwd velocity
	if (do_2_vel)
{
        const D e=(state(i+1,iu)-state(i,iu))/hf;
        D vx=m_mun*e;
        if (vx>vsat) vx=vsat;
        if (vx<-vsat) vx=-vsat;
        state(i,m_v_n)=vx;//m_mun*e;
        vx=-m_mup*e;
        if (vx>vsat) vx=vsat;
        if (vx<-vsat) vx=-vsat;
        state(i,m_v_p)=vx;

}
else {
        const D xrev=m_fixed(i-1,m_x);
        const D hr=x-xrev;
        const D hs=hr+hf;
        const D e=(state(i+1,iu)-state(i-1,iu))/hs;
        D vx=m_mun*e;
        if (vx>vsat) vx=vsat;
        if (vx<-vsat) vx=-vsat;
        state(i,m_v_n)=vx;//m_mun*e;
        vx=-m_mup*e;
        if (vx>vsat) vx=vsat;
        if (vx<-vsat) vx=-vsat;
        state(i,m_v_p)=vx;
		}
		// for the heck of it track rho too
		state(i,m_q)=rho_fixed(i)+state(i,m_p)-state(i,m_n);
    } // i 
}



// addd traps as variable fixed amounts lol. 
D rho_fixed(const IdxTy point)
{
return m_fixed(point,m_nd)-m_fixed(point, m_na)+m_state(point,m_trap);

}

// normally this could be local but due to drifts ti can be unpredictable 
// diffuse the carriers at in using drift velocities in iv to shift startic
// thermal velocity distro 
// this is the only thing that was interesting lol 
//void old_fick_lanc_flow(MyBlock & ndiff, const MyBlock & sol, const D & t, const D & dxmax, const int in, const int iv)

//////////////////////////////////////////////////////////////
/*
 x=which(df$V1==2)
> df$V7[x[1:10]]
 [1] 4500000 2250000 4500000 4500000 4500000 4500000 4500000 4500000 4500000
[10] 4500000
> x=which(df$V1==3)
> df$V7[x[1:10]]
 [1] 4.50e+06 1.50e+10 1.50e+10 1.50e+10 1.05e+07 4.30e+06 4.14e+06 4.18e+06
 [9] 4.22e+06 4.26e+06



*/


#define DEBUG_CHECK_FL	const D check=n*gdvself+nfwd+nrev; const bool check_lo=(check<n*0*(1.0-1e-2)); const bool check_hi=(check>n*(1.0+1e-10)); if ( check_lo||check_hi) { if (check_hi) { MM_INC_MSG(m_cm,"flconservehi" )  } if (check_lo) { MM_INC_MSG(m_cm,"flconservelo" )  } MM_ERR(" fl_flow err "<<i<<" n "<<n<<" check "<< check <<" gdvself "<<gdvself<<" fwd "<<nfwd<<" nrev "<<nrev)  }









///////////////////////////////////////////////////////////////////////

 // sf_bc
// these thigns fudding overlap just like the other fudding ones.... doh!
/*
 This is the BC handler for the most complete impl so far- integrated over sf.
There is no field and the BC just needs to find the isotropic flow into and 
out of the infinite void. 

// USING

*/
//return qflefs.diff_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
template <class Q> D diff_vec_bc(const IdxTy br, Q & qfle, const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2,  const D & xmax, const bool no_self, const bool dbg) const
{
if (br==1) { MM_MSG(" removed in v 9") return 0; }  //  return qfle.diff_sf_loop(xp0,xp1,xp2,x0,x1, x2,n0,n1,n2 , xmax,no_self,!true);
return qfle.diff_sf_bc(xp0,xp1,xp2,x0,x1, x2,n0,n1,n2 , xmax,no_self,!true);

}

class fl_branches
{
public:
template <class Tc> fl_branches(const Tc & tree):
	use_perf(tree.track_fl_delta()),
	 keep_the_lost(true), add_from_beyond(true), debug_1(false ), dbg(false ),
	 dbg2(false ), dbgg(false), dbgxo(false ), no_self(!true) {}
template <class Tc> fl_branches(const Tc & tree, const bool dbx):
	use_perf(tree.track_fl_delta()),
	 keep_the_lost(true), add_from_beyond(true), debug_1(false ), dbg(dbx ),
	 dbg2(false ), dbgg(false), dbgxo(false ), no_self(!true) {}


	bool use_perf; // =m_tree.track_fl_delta();
	bool keep_the_lost; // =true;
	bool add_from_beyond; // =true;
	bool debug_1; // =false; 
	bool dbg; // =false; 
	bool dbg2; // =false; 
	bool dbgg; // =false; 
	bool dbgxo; // =false; 
	bool no_self; // =!true; 
}; //fl_branches



class pixel_block
{
public:
template <class Tq> void in_and_out_left(D & out, D & in, pixelt & x,const D & xmax, const Di & di , const fl_branches & flb, const Tq & qfle
,const D & nbc, const D & nbc1, const D & ish)
{
	D selff=(true)?.5:1;
	const D x0=di.x0();
	const D x1=di.x1();
	const D x2=di.x2();
	const D offset=x0-di.left();
	const D dumdim=x1-x0; // need to think about this now
	const D dx=xmax-offset;
	const D faccin=selff*ish/xmax/xmax;
	const D faccout=selff*ish/xmax/xmax;
	const IdxTy pixmax=IdxTy((dx)/(x1-x0)+2) ;
	pixelt bc[pixmax];
	for (IdxTy k=0; k<pixmax; ++k)
		bc[k].setmid(di.left()-dumdim*k,dumdim);
	for (IdxTy k=0; k<pixmax; ++k){

	 if (k==0) in+=faccin*bc[k].diff_sf(x,qfle,nbc,nbc,nbc1,xmax,flb.no_self,flb.dbg);
	 else in+=faccin*bc[k].diff_sf(x,qfle,nbc,nbc,nbc,xmax,flb.no_self,flb.dbg);
	 out+=faccout*x.diff_sf(bc[k],qfle,di.n0(),di.n1(),di.n2(),xmax,flb.no_self,flb.dbg);
//		net_left=in_from_left-out_to_left;
	}

}  // in_and_out_left

template <class Tq> void in_and_out_right(D & out, D & in, pixelt & x,const D & xmax, const Di & di , const fl_branches & flb, const Tq & qfle
,const D & nbc, const D & nbc1, const D & ish)
{
	D selff=(true)?.5:1;
	const D x0=di.x0();
	const D x1=di.x1();
	const D x2=di.x2();
	//const D offset=x0-di.left();
	const D offset=-x2+di.right();
	const D dumdim=x1-x0; // need to think about this now
	const D dx=xmax-offset;
	const D faccin=selff*ish/xmax/xmax;
	const D faccout=selff*ish/xmax/xmax;
	const IdxTy pixmax=IdxTy((dx)/(x1-x0)+2) ;
	pixelt bc[pixmax];
	for (IdxTy k=0; k<pixmax; ++k)
		bc[k].setmid(di.right()+dumdim*k,dumdim);
	for (IdxTy k=0; k<pixmax; ++k){

	 if (k==0) in+=faccin*bc[k].diff_sf(x,qfle,nbc1,nbc,nbc,xmax,flb.no_self,flb.dbg);
	 else in+=faccin*bc[k].diff_sf(x,qfle,nbc,nbc,nbc,xmax,flb.no_self,flb.dbg);
	 out+=faccout*x.diff_sf(bc[k],qfle,di.n0(),di.n1(),di.n2(),xmax,flb.no_self,flb.dbg);
//		net_left=in_from_left-out_to_left;
	}

}  // in_and_out_left




}; // pix_block

D fi_bc_loop(MyBlock & ndiff, Di & di, const MyBlock & sol,  const IdxTy in, const D & xmax, const D & vt,const D & ish,const D & n , const fl_branches & flb )
{
	const IdxTy i=di.i();
	const D x0=di.x0();
	const D x1=di.x1();
	const D x2=di.x2();
	const bool hits_left=di.hits_left(xmax+(x2-x0));
	const bool hits_right=di.hits_right(xmax+(x2-x0));
	if (!hits_left&&!hits_right) return 0 ; 
	D selff=(true)?.5:1;
	QFLE qfle;
	pixelt lbc, rbc,x;
	x.set(x0,x1,x2);	
	D dumdim=xmax;
	D mindumdim=2.0*(x2-x0);
	IdxTy pixmax=IdxTy((xmax)/(x2-x0)+2) ;
	if (dumdim<mindumdim) dumdim=mindumdim; 
	// these need to be on the middle to overlap the rel pixels doh
	dumdim=x1-x0; // need to think about this now
	MM_ONCE(" need to thing about the diff bc doh",)
	// unfonrunately we now need many of these as with drift
	// the size should not matter to the physics, there is 
	// something simply wrong here. 
	lbc.setmid(di.left(),dumdim);
//	lbc.set(di.left(),-dumdim);
	//rbc.set(di.right(),dumdim);
	rbc.setmid(di.right(),dumdim);
	const D n0=di.n0();
	const D n1=di.n1();
	const D n2=di.n2();
	const D faccin=selff*ish/xmax/xmax;
	const D faccout=selff*ish/xmax/xmax;
	D net_left=0;
	D net_right=0;
	// amount from lbc INTO dest pixel x. 
	// do NOT do this for every pixel lol 
	if (hits_left)
	{
		const D nbc=sol(0,in);
		const D nbc1=sol(1,in);
		//const D in_from_left=faccin*lbc.diff_sf(x,qfle,nbc,nbc,nbc,xmax,flb.no_self,flb.dbg);
	//	const D in_from_left=faccin*lbc.diff_sf(x,qfle,nbc,nbc,nbc1,xmax,flb.no_self,flb.dbg);
//		const D out_to_left=faccout*x.diff_sf(lbc,qfle,n0,n1,n2,xmax,flb.no_self,flb.dbg);
	pixel_block bc;
	D in_from_left=0, out_to_left=0;
	bc.in_and_out_left(out_to_left, in_from_left, x, xmax,  di,  flb
	,qfle,nbc,nbc1, ish);
		net_left=in_from_left-out_to_left;
		ndiff(i)+=net_left; // in_from_left-out_to_left;		 	
		MM_MSG(MMPR(i)<<MMPR(net_left)<<MMPR(in_from_left)<<MMPR(out_to_left)<<MMPR(ndiff(i))<<MMPR(nbc)<<MMPR(n0)<<MMPR(n1)<<MMPR(n2)) 
	} // hits_left
	if (hits_right)
	{
		const D nbc=sol(m_points-1,in);
		const D nbc2=sol(m_points-2,in);
		pixel_block bc;
	D in_from_right=0, out_to_right=0;
	bc.in_and_out_right(out_to_right, in_from_right, x, xmax,  di,  flb
	,qfle,nbc,nbc2, ish);
		//const D in_from_right=faccin*rbc.diff_sf(x,qfle,nbc,nbc,nbc,xmax,flb.no_self,flb.dbg);
//		const D in_from_right=faccin*rbc.diff_sf(x,qfle,nbc2,nbc,nbc,xmax,flb.no_self,flb.dbg);
//		const D out_to_right=faccout*x.diff_sf(rbc,qfle,n0,n1,n2,xmax,flb.no_self,flb.dbg);
		net_right=in_from_right-out_to_right;		 	
		ndiff(i)+=net_right; // in_from_right-out_to_right;		 	
		MM_MSG(MMPR(i)<<MMPR(net_right)<<MMPR(in_from_right)<<MMPR(out_to_right)<<MMPR(ndiff(i))<<MMPR(nbc)<<MMPR(n0)<<MMPR(n1)<<MMPR(n2)) 

	} // hits_left
//template <class Tq > const D diff_sf(pixelt & xp, Tq & qflefs, const D & n0, const D & n1, const D & n2, 
//const D & xmax, const bool no_selfs, const bool dbg)
//return qflefs.diff_sf_loop(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
//return qflefs.diff_sf_bc(xp.x0,xp.x1,xp.x2, x0, x1, x2, n0, n1, n2, xmax, no_selfs, dbg); // /xmax; 
//}

//return qfle.diff_sf_bc(xp0,xp1,xp2,x0,x1, x2,n0,n1,n2 , xmax,no_self,!true);


return net_left+net_right;
} // fi_bc_loop



//////////////////////////////////////////////////////////////////////


/*
Raison d'etre. Right now the problem is boundary condiations probably due to
inconsistent flow using integrated element but subtracting from the center only.

This naturally leads to shape functions lol 

*/


/*
void fick_lanc_flow_subpixel(MyBlock & ndiff, const MyBlock & sol, const D & t, const int in, const D vfac=1.0)
{
*/


/////////////////////////////////////////////////////////////////////////////////
// USING
/*
More or less complete impl of FL flow with "shape functions."
Not sure what to do with overlap yet and small xmax.  Right now
this seems to more or less work except at boundary, mostly the left
side with the pnj test case. Fix BC handler and try to get quantitative
result, add traps etc.  

*/



void fick_lanc_flow_sf(MyBlock & ndiff, const MyBlock & sol, const D & t, const int in, 
	const D vfac=1.0, const bool ddbg=false)
{
	MySparse tofrom;
	const bool track_to_from=!false;
	CounterMap perf;
	QFLE qfle;
	fl_branches flb(m_tree);
	const StrTy pname="perf";
	const IdxTy sz=m_points-1;
	ndiff.resize(m_points);
	const FlpTy  & flp=m_flp;
	const D vmax=flp.vmax()  *vfac; // 1e8;
	MM_INC_MSG(m_cm," enable_einstein");
	const D xmax=vmax*t;
	const D vt=0;
	Di di(m_points,m_x,in,sol,m_fixed);
	while (di.ok())
	{
		// this goes invaldi at bad times in itor 
		const int i=di.i();
		D nfwd=0;
		D nrev=0;
		const D ish=1.0/di.sh();
		// this appears to be too low by a factor of two 
		const D n=2.0*qfle.self_sf_loop(di,!true)*ish;
		// the other user of this has a factor of .5 wtf? 
		const D neff=qfle.isselflin(di)*ish;
		D bc_gain=0;	
		if (flb.add_from_beyond)  
			bc_gain=fi_bc_loop( ndiff,di,sol,in, xmax, vt, ish,n,flb );
		if (n==0) { di.inc(); continue; }
		const D facc=.5*ish/xmax/xmax;
		while (0==di.state())
		{
			const int j=di.j();
			const D gdv= facc*qfle.diff_sf_loop( di, xmax, vt,flb.no_self,false); 
			// this makes the itor fudding invalid 
			const bool bump=di.result(gdv);
			if (j!=0) { ndiff(j)+=gdv;	nrev+=gdv;    
			if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=gdv;
			if (track_to_from) if ((i>(m_points-3))||(j>(m_points-3)))tofrom(i,j)+=gdv;
			}
			if (flb.use_perf) perf.add(pname,j-i,gdv/neff);
			// fudding iterator is fudded ... 
			if (bump) break;
		} // 0==state
		while (1==di.state())
		{
			const int j=di.j();
			if (j==sz) { di.result(0); break ; } // wtf ??? 
			const D gdv=   facc*qfle.diff_sf_loop( di, xmax, vt,flb.no_self,false); 
			const bool bump=di.result(gdv);
		//	if ( j!=sz) 
			{ nfwd+=gdv; ndiff(j)+=gdv;	}		
			if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=gdv;
			if (track_to_from) if ((i>(m_points-3))||(j>(m_points-3)))tofrom(i,j)+=gdv;
			if (flb.use_perf) perf.add(pname,j-i,gdv/neff);
			if (bump) break; 
		} // 1==state
		const D gdvself= facc*qfle.diff_sf_self_loop(di,  xmax, vt, !true );
		const D fmov= ((nfwd+nrev)/neff);
		if ((fmov>=1.0)||(nfwd<0)||(nrev<0)||(gdvself<0)) 
		{ 
		 	{ MM_INC_MSG(m_cm,"fracmoved" )  }
			MM_MSG(" fraction moved i="<<i<<" frac="<<fmov<<" fwd "<<nfwd<<" rev "<<nrev<<" n "<<n<<" self "<<gdvself)
		}
// wow this is fatal 
//		if (flb.keep_the_lost) ndiff(i)+=n-nfwd-nrev;
		if (flb.keep_the_lost){
 //ndiff(i)+=sol(i,in)-nfwd-nrev+bc_gain;
 const D selfn=sol(i,in)-nfwd-nrev+bc_gain;
 ndiff(i)+=selfn; // sol(i,in)-nfwd-nrev+bc_gain;
			if (track_to_from) if ((i<5)||(i<5))tofrom(i,i)+=selfn;

}
		else ndiff(i)+=gdvself;
		if (flb.use_perf) perf.add(pname,0,(gdvself/neff));
		if (bc_gain!=0)
			{
			MM_MSG("bc diff "<<i<<MMPR(bc_gain)<<MMPR(ndiff(i))<<MMPR(nfwd)<<MMPR(nrev)<<MMPR(n)<<MMPR(neff)<<MMPR(sol(i,in)))
			}
			// the itor is invalid 	
	} // while
	// these should be ok now 
	ndiff(0)=sol(0,in);
	ndiff(sz)=sol(sz,in);
	for(IdxTy fudd=0; fudd<sz; ++fudd) 
	{
		const D v=ndiff(fudd); 	
		if (v<0) { { MM_MSG(" negative carrier found at "<<fudd<<" n="<<ndiff(fudd) ) }  
			ndiff(fudd)=::sqrt(m_ni2); }  } 

	if (flb.use_perf) { MM_MSG(" fick_lanc_perf "<<perf.array_to_string(pname)) }
	if (track_to_from)
	{
		tofrom.dump(std::cout,7,"drift_to_from");
	}

} // fick_lanc_flow_sf

/////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
/*
void fick_lanc_flow_subpixel_old(MyBlock & ndiff, const MyBlock & sol, const D & t, const int in)

*/

void fick_lanc_flow_no_drift(MyBlock & ndiff, const MyBlock & sol, const D & t, const int in)
{
	CounterMap perf;
	MySparse tofrom;
	const bool track_to_from=!false;
	const bool use_perf=m_tree.track_fl_delta();
	// more physical to make a fake bulk off both ends returning particles
	// but this should be ok for now. The ends are becoming depleted lol 
	const bool keep_the_lost=true;
	const bool replace_the_lost=!true;
	const bool debug_1=false; 
	const StrTy pname="perf";

	ndiff.resize(m_points);
	const FlpTy  & flp=m_flp;
	const D vmax=flp.vmax(); // 1e8;
//	const D norm=flp.norm(); // .5/vmax;
	const D xmax=vmax*t;
	const IdxTy sz=m_points-1;
	for (IdxTy i=1; i<sz;  ++i)
//	for (IdxTy i=0; i<m_points;  ++i)
	{
		const D n=sol(i,in);
		const D x=m_fixed(i,m_x);
		D nfwd=0;
		D nrev=0;
		int j=i-1;
		while (j>0)
		{
			const D dx=x-m_fixed(j,m_x);
			// this is an artificial loss, not at boundary, 
			// if (false) if (dx>dxmax) { break;}	
			const D gdv=int_gdv(j,i,xmax,0);
			if (gdv<=0) { break; } // we could break except that vdrift may vary 
			const D frac=n*gdv;
//			MM_MSG(" fudd "<<frac<<" n "<<n<<" gdv "<<gdv<<" dh "<<dh<<" dxd "<<dxd<<" xmax "<<xmax<<" vmax "<< vmax<<" vdrift "<<vdrift)
				if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=frac;
			nrev+=frac; ndiff(j)+=frac;		 	
			if (use_perf) perf.add(pname,i-j,gdv);
		if (j<10) 
		{ // these seem to sum ok 
			if (debug_1)
				{ MM_ERR(" j ffrac "<<j<<" from "<<i<<"ndiff(j) "<<ndiff(j)<<" gdv "<<gdv<<" ngdv "<<(frac))}
		}
			--j;
		} // j 
		j=i+1;
		while (j<sz)
		{
			const D dx=m_fixed(j,m_x)-x;
			// if (false) if (dx>dxmax) 	{ break;	}	
			const D gdv=int_gdv(j,i,xmax,0);
			if (gdv<=0) { break; } // we could break except that vdrift may vary 
			const D frac=n*gdv;
				if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=frac;
			nfwd+=frac; ndiff(j)+=frac;		 	
			if (use_perf) perf.add(pname,i-j,gdv);
		if (j<10) 
		{ // these seem to sum ok 
			if (debug_1)
			{ MM_ERR(" j ffrac "<<j<<" from "<<i<<"ndiff(j) "<<ndiff(j)<<" gdv "<<gdv<<" ngdv "<<(frac)) }
		}
			++j;
		} // j 
		const D fmov= ((nfwd+nrev)/n);
		if (fmov>=1.0) { MM_MSG(" fraction moved "<<fmov<<" "<<nfwd<<" "<<nrev)}
		const bool i_ok=(i>0)&&(i<sz);	
		const D gdvself=i_ok?int_gdv(i,i,xmax,0):0;
		if (use_perf) perf.add(pname,0,gdvself);
		if (keep_the_lost) ndiff(i)+=n-nfwd-nrev;
		else ndiff(i)+=n*gdvself;
		if (i<10) 
		{ // these seem to sum ok 
			if (debug_1)
			{ MM_ERR(" i frac "<<i<<" gdvself "<<gdvself<<" fwd "<<(nfwd/n)<<" rev "<<(nrev/n)) } 
		}
		const D check=n*gdvself+nfwd+nrev;
// this trips out at 1 and 998 boundaried 
//		const bool check_lo=(check<n*(1.0-1e-2));
		const bool check_lo=(check<n*0*(1.0-1e-2));
		const bool check_hi=(check>n*(1.0+1e-10));
		//if ((check>n*(1.0+1e-10)) ||(check<0))
		if ( check_lo||check_hi)
		{
			if (check_hi) { MM_INC_MSG(m_cm,"flconservehi" )  }
			if (check_lo) { MM_INC_MSG(m_cm,"flconservelo" )  }
			MM_ERR(" fl_flow err "<<i<<" n "<<n<<" check "<< check <<" gdvself "<<gdvself<<" fwd "<<nfwd<<" nrev "<<nrev)  
		}
		if (replace_the_lost) // this does not get those lost due to dxmax cutoff
		{	// this pixel is nearer than xmax to a bc. It obtains n over a collecting
			// region .5*(h_forward+h_back) from the continuum beyond the bc with
			// bc carrier density. For now just ignore density and n*h lol
//		const D x0=m_fixed(0,m_x);
//		const D xd=m_fixed(sz,m_x);
		const D n0=sol(0,in);
		const D nsz=sol(sz,in);
		IdxTy kluge=0; 
		while (true)
		{
		// integrate over fake points ASSUMING const h doh 
			if (kluge>=i) break;	
			const D gdvlost=int_gdv(i-kluge,sz,xmax,0);
			if (gdvlost<=0) break;
			ndiff(i)+=nsz*gdvlost;
			++kluge;
		}
		kluge=0;
		while (true)
		{
		// integrate over fake points ASSUMING const h doh 
		// this looks like it may be working but the integral needs to be over
		// source and destination or there is some other error. 
			if ((kluge+i)>=sz) break;	
			const D gdvlost=int_gdv(i+kluge,0,xmax,0);
			if (gdvlost<=0) break;
			const D old=ndiff(i);
			ndiff(i)+=n0*gdvlost;
			if (debug_1)
			{ MM_ERR(" adding lost "<<i<<" n0 "<<n0<<" kluge "<<kluge<<" gdv "<<gdvlost<<" ndiff "<<ndiff(i)<<" was "<<old <<" pct "<<(ndiff(i)/n0))
			}
			++kluge;
		}


		} // replace_the_lost
		// this is not right, the ones that moved too far need to be eliminated 
		//ndiff(i)-=(nfwd+nrev);
	} // i 
	if (use_perf) { MM_MSG(" fick_lanc_perf "<<perf.array_to_string(pname)) }

} // fick_lanc_flow






/////////////////////////////////////////////////////////////
D limp( const D & x) { if (x>1.0) return 1.0; return x; }
D limn( const D & x) { if (x<-1.0) return -1.0; return x; }

// the flow was now expanded to terminal pixels. So +/- need to be
// checked. 
D int_gdv(const IdxTy point, const IdxTy from, const D & xmax, const D & xoffset)
{
	const D d=m_fixed(point,m_x);
	const D c=m_fixed(from,m_x);
	const D b=.5*(m_fixed(point+1,m_x)+d);
	const D a=.5*(d+m_fixed(point-1,m_x));
	const D peak=c+xoffset;
	const D xnb=(b-peak)/xmax;
	const D xna=(a-peak)/xmax;
	// round off error lol 
	if (xnb>=1) if (xna>=1) return 0;
	if (xnb<=-1) if (xna<=-1) return 0;

	D dxd=0;
	if (b>peak) { D xn=limp(xnb); dxd+=xn-.5*xn*xn; }
	else { D xn=limn(xnb); dxd+=xn+.5*xn*xn; } 
	if (a>peak) { D xn=limp(xna); dxd-=xn-.5*xn*xn; }
	else { D xn=limn(xna); dxd-=xn+.5*xn*xn; } 
	if ((a-peak)*(b-peak)<0)
	{
		// two integrals, both evaulated at zero do not seem to matter. 

	}
	D gdv=dxd; //	const D gdv=(xmax-dxd)/(xmax*xmax);
	if ((gdv<0)||(gdv>1))
	{
		MM_ERR(" gdv "<<gdv<<" point "<<point<<" from "<<from<<" xmax "<<xmax<<" xoffset "<<xoffset<<" b "<<b<<" a "<<a<<" peak "<< peak )
		MM_ERR(" gdv2 "<<gdv<<" point "<<point<<" from "<<from<<" dxb "<<(b-peak)<<" dxa "<<(peak-a)<<" dpeak "<< (peak-c)<<" dxpoint "<<(d-c) )
	}
//	if (use_perf) perf.add(pname,i+j,gdv);
	return gdv;
}

 


// find the elemen less than or equal to target 
IdxTy floor(const MyBlock & loc, const IdxTy max, const IdxTy col, const D target, const IdxTy start)
{
IdxTy i=start;
D v=loc(i,col);
while (v>target)
{
if (i==0) return ~0;
--i;
v=loc(i,col);
if (v<=target) return i;
}
while (v<target)
{
// max is normally size-1
if (i==max) return max+1;
++i;
v=loc(i,col);
if (v==target) return i;
if (v>target) return i-1;
}

return i;
}

// move the carrier ballistically in the field, should not go more than a pixel or two 
/*
void do_drift(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, const int iv)

*/

/*
// this needs to save all the i dependent crap 
D  drift_subpixel(const IdxTy to, const IdxTy from,  const MyBlock & sol, const D & t,  const int in, const int iv, const IdxTy points)

return remove ;
}
*/

	 		//D dnji= drift_sf(di, sol,  t,  in,  iv);
D  drift_sf(Di & di,  const MyBlock & sol, const D & t,  const int in, const int iv
	, const drift_flags & df)
{
QFLE qfle;
	const int from=di.i();
	const int to=di.j();
	// the iterator should have this since it has the x coords etc. 
	const D v0=di.v0(); // .5*(sol(from,iv)+sol(from-1,iv));
	const D v1=di.v1(); // .5*(sol(from+1,iv)+sol(from,iv));

	const D bv=.5*(v0+v1); // yes, this is mf(from, iv) with uniform h lol 
	const D delta_x=t*bv;
const bool do_2_vel=true;
if (do_2_vel) { MM_ONCE(" doing split velocity ",) } 
  	//removing=is(a,  b, xleft, x, xright, nleft, nzed, nright)/(xright-xleft);
//	D remove=qfle.(xip0,xip1,x0,xm,x1,.5*(nn1+n0),n0,.5*(n0+n1))/(x1-x0);


//	D remove=qfle.drift_sf(di.xp0(), di.xp1(), di.xp2(), di.x0(), di.x1(), di.x2()
//		, di.n0(), di.n1(), di.n2(), delta_x, false); 
//	MM_ONCE(" should be using the new simplyer drift thing may be all fudded up",)
//const drift_flags df(false,true,false);
if (!do_2_vel)
{
	D removenew=2.0*qfle.drift_sf_bc(di.xp0(), di.xp1(), di.xp2(), di.x0(), di.x1(), di.x2()
		//, di.n0(), di.n1(), di.n2(), delta_x, false,false); 
		, di.n0(), di.n1(), di.n2(), delta_x, df) ; // ZZ, dbg); 
return removenew;
}
const D xleft=sol(from-1,iv)*t;
drift_flags df2=df;
df2.only_left_src=true;
	D removenew_left=2.0*qfle.drift_sf_bc(di.xp0(), di.xp1(), di.xp2(), di.x0(), di.x1(), di.x2()
		//, di.n0(), di.n1(), di.n2(), delta_x, false,false); 
		, di.n0(), di.n1(), di.n2(), xleft, df2) ; // ZZ, dbg); 

df2.only_right_src=true;
df2.only_left_src=!true;
const D xright=sol(from,iv)*t;
	D removenew_right=2.0*qfle.drift_sf_bc(di.xp0(), di.xp1(), di.xp2(), di.x0(), di.x1(), di.x2()
		//, di.n0(), di.n1(), di.n2(), delta_x, false,false); 
		, di.n0(), di.n1(), di.n2(), xright, df2) ; // ZZ, dbg); 



//	D removenewbw=2.0*qfle.drift_sf_bc(di.xp0(), di.xp1(), di.xp2(), di.x0(), di.x1(), di.x2()
//		, di.n0(), di.n1(), di.n2(), -delta_x, false,false); 


// quick spot check shows similarly and the plots of test output overlap.
// note that the bw one varies a lot however, may be useful to know  
//	MM_MSG(" drift bc compare "<<MMPR(remove)<<MMPR(removenew)<<MMPR(removenewbw))
//return remove ;
return removenew_left+removenew_right ;
}





//	diffuse_iterator(const IdxTy points, const IdxTy _x, const IdxTy _in, const MyBlock & sol, const MyBlock & f)
//	: m_state(~0),m_points(points), m_x(_x),m_in(_in),m_sol(sol), m_f(f) {Init();}

D  drift_subpixel( Di & di, const MyBlock & sol, const D & t,  const int in, const int iv )
{
	QFLE qfle;
	const int from=di.i();
	const int to=di.j();
	// the iterator should have this since it has the x coords etc. 
	const D v0=di.v0(); // .5*(sol(from,iv)+sol(from-1,iv));
	const D v1=di.v1(); // .5*(sol(from+1,iv)+sol(from,iv));

	const D av=0;
	const D bv=.5*(v0+v1); // yes, this is mf(from, iv) with uniform h lol 
//	const D delta_x=t*bv;
	D xip0=(di.a()-t*bv)/(1.0+t*av);
	D xip1=(di.b()-t*bv)/(1.0+t*av);
//	MM_MSG(" itor "<<di.state_string()<<" "<<di.const_string()<<" xip0="<<xip0<<" xip1="<<xip1<<" tbv="<<(t*bv))
	D remove=qfle.is(xip0,xip1,di.x0(),di.x1(),di.x2(),di.n0(),di.n1(),di.n2())/di.sh();
	return remove;

}

// USING 
void do_sf_drift(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, 
const int iv, const drift_flags & df, const debug_levels &  dbg=false)
{
	// this was really meant for diffusion not drift
	fl_branches flb(m_tree,dbg);
	CounterMap perf;
	MySparse tofrom;
	QFLE qfle;
	const StrTy pname="perf";
	const bool track_to_from=!false;
	const IdxTy sz=m_points-1;
	ndiff.resize(m_points);
	Di di(m_points,m_x,in,iv,sol,m_fixed);
	const D f1=.5;
//	const drift_flags df(false,true,false);
	while (di.ok())
	{
		const IdxTy i=di.i();
		const D ish=1.0/di.sh();
		/// the test has zdero leftovers and this was not important
		// will need df wft
		const D neff=2.0*f1*qfle.isselflin(di)*ish;
		const D guard=di.sh();
		const D dxmax=t*di.v1(); // +guard;
		const D xmax=di.x2()+dxmax +guard;
		const D xmin=di.x0()+dxmax -guard;
		const D n0=di.n0(); 
		const D n1=di.n1(); 
		const D n2=di.n2();
//		if (neff!=0) {MM_MSG(MMPR(i)<<MMPR(n0)<<MMPR(n1)<<MMPR(n2)<<MMPR(neff)) } 		
		D gone=0;
// must be done before itor increments. 
		if (flb.add_from_beyond) 
			gone=sf_drift_bc( ndiff, di,  sol, dxmax,  xmin,  xmax, sz, ish, in,flb, df ); 
	//	if (track_to_from) 
		if (gone!=0) {MM_MSG(" bc losses "<<MMPR(i)<<MMPR(gone))}
		if (gone<0) gone=0; 
		// neff -> n1 is zero but not the other way 
		if (neff==0) { di.inc(); continue; } 
		// not sure what this fickes up
//		if (n1==0) { di.inc(); continue; } 
		// testhing would not show an ish factor doh
		D facc=f1;
		if (!df.uniform_src) facc=facc*2.0;
		while (0==di.state())
		{
			const int j=di.j();
			const D jx=m_fixed(j,m_x); // these limits are really risky 
	 		D dnji= facc*drift_sf(di, sol,  t,  in,  iv,df)*ish;
			if (dnji<0) dnji=0; // should just be numerical noise 
			if ((j+1)==i)
			{
	 			D dnjized= facc*drift_sf(di, sol,  0,  in,  iv,df)*ish;
				//MM_MSG( " FUDDING EXECUTE "<<MMPR(i)<<MMPR(j)<<MMPR(dnji)<<MMPR(dnjized))
				dnji-=dnjized;
				if (dnji<0) dnji=0;
			}
			//if (flb.use_perf) perf.add(pname,int(j)-int(i),dnji/neff);
			if (j!=i)
			{
				//if (flb.use_perf) perf.add(pname,int(j)-int(i),dnji/n1);
				if (flb.use_perf) perf.add(pname,int(j)-int(i),dnji/neff);
				if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=dnji;
				if (track_to_from) if ((i>(m_points-4))||(j>(m_points-4)))tofrom(i,j)+=dnji;
				gone+=dnji; ndiff(j)+=dnji;
			}
			//if ((j>1996)||(i>1996)) {MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			//if (dnji!=0) {MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			//{MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			const bool bump=di.result(jx-xmin);
			if (bump) break; 
		}
		while (1==di.state())
		{
			const int j=di.j();
			const D jx=m_fixed(j,m_x);
	 		D dnji= facc*drift_sf(di, sol,  t,  in,  iv,df)*ish;
			if (dnji<0) dnji=0; // should just be numerical noise 
			if ((j)==(i+1))
			{
	 			D dnjized= facc*drift_sf(di, sol,  0,  in,  iv,df)*ish;
				//MM_MSG( " FUDDING EXECUTE "<<MMPR(i)<<MMPR(j)<<MMPR(dnji)<<MMPR(dnjized))
				dnji-=dnjized;
				if (dnji<0) dnji=0;
			}
			// we could break here but the iterator again... fudd 
			if (j!=i) 
			{
				//if (flb.use_perf) perf.add(pname,int(j)-int(i),dnji/n1);
				if (flb.use_perf) perf.add(pname,int(j)-int(i),dnji/neff);
				//if (track_to_from) tofrom(i,j)+=dnji;
				if (track_to_from) if ((i<5)||(j<5))tofrom(i,j)+=dnji;
				if (track_to_from) if ((i>(m_points-4))||(j>(m_points-4)))tofrom(i,j)+=dnji;
				gone+=dnji; ndiff(j)+=dnji;
			}
			//if ((j>1996)||(i>1996)) {MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			//if (dnji!=0) {MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			//{MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(dnji)) }
			const bool bump=di.result(xmax-jx);
			if (bump) break; 
		}
//		const D twelth=(2.0*n1-n0-n2)/12.0;
		//ndiff(i)+=neff+twelth-gone; // neff-gone;
		MM_ONCE(" neff has worked well except for small drifts ",)
		const D leftover=n1-gone;
		//const D leftover=neff-gone;
		if (leftover>0){ ndiff(i)+=leftover; // n1+0*twelth-gone; // neff-gone;
				//if (track_to_from) tofrom(i,i)+=leftover;
				if (track_to_from) if ((i<5)||(i<5))tofrom(i,i)+=leftover;
//			MM_MSG(" leftover "<<i<<MMPR(neff)<<MMPR(gone))
		 //if (flb.use_perf) perf.add(pname,0,leftover/n1); } 
		 if (flb.use_perf) perf.add(pname,0,leftover/neff); } 

	} // ok 
	ndiff(0)=sol(0,in);
	ndiff(sz)=sol(sz,in);
	for(IdxTy fudd=0; fudd<sz; ++fudd) 
	{
	
		if (ndiff(fudd)<0) 
			{ { MM_MSG(" negative carrier found at "<<fudd<<" n="<<ndiff(fudd) ) }  
			ndiff(fudd)=0*::sqrt(m_ni2); }  
	} 

	if (flb.use_perf) { MM_MSG(" sf drift perf "<<perf.array_to_string(pname)) }
	if (track_to_from)
	{
		tofrom.dump(std::cout,7,"drift_to_from");

	}

} // do_sf_drift



/*
//D ndiff=0;
void do_drift_subpixel_itor(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, const int iv)
*/
class bc_pixels
{

public:
bc_pixels() {}
bool set(pixelt & x, const D & ref0, const D & ref1, const D & vt, const bool right)
{
bool nearby=false;
m_missing=0;
m_right=right;
D x0=0; D x1=0; D x2=0;
// for now all elements are same size but will later need to fix 
// note that these OVERLAP 
D pitch=1*(ref1-ref0); 
D absvt=fabs(vt);
if ((!right)&&(x.x0<((ref1+absvt)))) nearby=true;
if ((right)&&(x.x2>((ref0-absvt)))) nearby=true;
if (!nearby) return false;
// need to make sure that these are right to capture
// orphan carriers based on vt
if (right)
{ 
D xmax=absvt+x.x2;
D off=0;
// these need to be quantized to pick n later 
// and set m_missing... 
D bound=ref1+2*pitch;
if (xmax>bound) off=xmax-bound;
if ( off>(2*pitch)) {m_missing = 2; } 
else if ( off>(pitch)) { m_missing = 2; off=2*pitch; } 
else if (off>0) { off=pitch; m_missing=1; } 
//else off=0; 
x0=ref0+off; x1=ref1+off; x2=ref1+(ref1-ref0)+off; }
else { 
D xmin=x.x0-absvt;
D off=0;
D bound=ref0-2*pitch;
if (xmin<bound) off=bound-xmin;
if ( off>(2*pitch)) {m_missing = 2; } 
else if ( off>(pitch)) { m_missing = 2; off=2*pitch; } 
else if (off>0) { off=pitch; m_missing=1; } 

x2=ref1-off; x1=ref0-off; x0=ref0-(ref1-ref0)-off; }

if (!right) pitch=-pitch;
for (IdxTy i=0; i<4; ++i)
{
p[i].set(x0,x1,x2);
// on the LEFT these DECREASE, the first element may touch grid 
x0+=pitch; x1+=pitch; x2+=pitch;
}

return nearby;
}
// this also needs the n values  - for loss these are on x
template <class QF, class DF > D loss(pixelt & x, const QF & qfle,
	const D & n0, const D & n1, const D & n2, const D & dmef, const DF & df  )
{
D sum=0;
for (IdxTy i=0; i<4; ++i)
{
	//D removing=facr*x.drift_sf_bc(p[i],qfle,n0,n1,n2, dmef, df); 
	const D removing=x.drift_sf_bc(p[i],qfle,n0,n1,n2, dmef, df); 
	const D removingzed=x.drift_sf_bc(p[i],qfle,n0,n1,n2, 0, df); 
	//MM_MSG(" removing p["<<i<<"]="<<removing)
	sum+=removing-removingzed;
}
return sum;
}

		//D adding=(faca*lbc.drift_sf_bc(x,qfle,nl,nl,nl, dmef, df)); 
template <class QF, class DF > D gain(pixelt & x, const QF & qfle,
	const D & nref0, const D & nref1, const D & nbc, const D & dmef, const DF & df  )
{
D sum=0;
	D n0=nref0;
	D n1=nref1;
MM_ONCE("subtracting zero field from drift ",)
// on the right, the first two values for the first pixel
// could be lattice values 
if (m_right)
{ 
	if (m_missing==1)  { n0=n1; n1=nbc; }
	if (m_missing==2)  { n0=nbc; n1=nbc; }
	const D adding0=p[0].drift_sf_bc(x,qfle,n0,n1,nbc, dmef, df); 
	const D adding1=p[1].drift_sf_bc(x,qfle,n1,nbc,nbc, dmef, df); 
	const D adding2=p[2].drift_sf_bc(x,qfle,nbc,nbc,nbc, dmef, df); 
	const D adding3=p[3].drift_sf_bc(x,qfle,nbc,nbc,nbc, dmef, df); 

	const D adding0z=p[0].drift_sf_bc(x,qfle,n0,n1,nbc, 0, df); 
	const D adding1z=p[1].drift_sf_bc(x,qfle,n1,nbc,nbc, 0, df); 
	const D adding2z=p[2].drift_sf_bc(x,qfle,nbc,nbc,nbc, 0, df); 
	const D adding3z=p[3].drift_sf_bc(x,qfle,nbc,nbc,nbc, 0, df); 

	sum+=adding0+adding1+adding2+adding3;
	const D sumzed=adding0z+adding1z+adding2z+adding3z;
	sum-=sumzed;
}
if (!m_right)
{
	if (m_missing==1)  { n0=nbc; n1=n0; }
	if (m_missing==2)  { n0=nbc; n1=nbc; }
	const D adding0=p[0].drift_sf_bc(x,qfle,nbc,n0,n1, dmef, df); 
	const D adding1=p[1].drift_sf_bc(x,qfle,nbc,nbc,n0, dmef, df); 
	const D adding2=p[2].drift_sf_bc(x,qfle,nbc,nbc,nbc, dmef, df); 
	const D adding3=p[3].drift_sf_bc(x,qfle,nbc,nbc,nbc, dmef, df); 

	const D adding0z=p[0].drift_sf_bc(x,qfle,nbc,n0,n1, 0, df); 
	const D adding1z=p[1].drift_sf_bc(x,qfle,nbc,nbc,n0, 0, df); 
	const D adding2z=p[2].drift_sf_bc(x,qfle,nbc,nbc,nbc, 0, df); 
	const D adding3z=p[3].drift_sf_bc(x,qfle,nbc,nbc,nbc, 0, df); 

//	MM_MSG(MMPR(adding0)<<MMPR(adding1)<<MMPR(adding2)<<MMPR(adding3))
	sum+=adding0+adding1+adding2+adding3;
	const D sumzed=adding0z+adding1z+adding2z+adding3z;
	sum-=sumzed;
}

return sum;
}

pixelt p[4];
bool m_right, m_mearby;
IdxTy m_missing;

}; // bc_pixels
D sf_drift_bc(MyBlock & ndiff, Di & di, const MyBlock & sol, const D & dxmax , const D & xmin, const D & xmax
	, const IdxTy sz, const D & ish, const IdxTy in
	, const fl_branches & flb, const drift_flags & df )
{
//	MyBlock & sol=m_state;
	const bool going_left=(dxmax<0);
	const bool going_right=(dxmax>0);
	const D field_decay_factor=2*.5; 
	MM_ONCE(" field decay factor = "<<field_decay_factor, )
	const D n0=di.n0();
	const D n1=di.n1();
	const D n2=di.n2();
	const D x0=di.x0();
	const D x1=di.x1();
	const D x2=di.x2();
//	pixelt rbc,lbc,rbc2,lbc2,lbc3,x;
	pixelt x; 
	bc_pixels rbc,lbc;
	x.set(x0,x1,x2);
	const D dmef=field_decay_factor*dxmax; 
	const IdxTy rrp0=m_points-2;
	const IdxTy rrp1=m_points-1;
	const bool near_right=rbc.set(x,m_fixed(rrp0,m_x),m_fixed(rrp1,m_x),dmef,true);
	const bool near_left=lbc.set(x,m_fixed(0,m_x),m_fixed(1,m_x),dmef,!true);
	const IdxTy i=di.i();
	D nr=sol(sz,in);
	D nl=sol(0,in);
	QFLE qfle;
	D gfac=1.0;
	if (!df.uniform_src) gfac=2.0;
	const D facr=gfac/di.sh();
	const D faca=gfac/di.sh();
	if (df.no_self!=flb.no_self) 
	{ MM_ONCE(" mismatch "<<MMPR(df.no_self)<<MMPR(flb.no_self),) } 
	if (df.dbg!=flb.dbgg) 
	{ MM_ONCE(" mismatch "<<MMPR(df.dbg)<<MMPR(flb.dbgg),) }
	// const drift_flags df(flb.no_self,true,flb.dbgg);
	if (near_right&&going_right)
	{  // this will not work if the decay factor is not used above 
//		MM_ONCE(" kluge foacter of .5 on right bc drift",)
//		D removing=.5*facr*rbc.loss(x,qfle,n0,n1,n2, dmef, df); 
		D removing=facr*rbc.loss(x,qfle,n0,n1,n2, dmef, df); 
		//MM_MSG(" rotal removed at i="<<i<<" is "<<removing<<" d "<<((rbc.p[0].x0-x.x2)*ish)<<" ndiff(i) was "<<ndiff(i))
		// all the other stuff is ADDITIVE, with remainders left.
		// this is a logic problem. need to just 
		// add back whatever was originally there. 
//		ndiff(i)=ndiff(i)-removing;
		// in theory we could be near BOTH right and left but not a concern now
		return removing;
	}
	if (near_left&&going_left)
	{
		//D dmef=field_decay_factor*dxmax; 
//		MM_ONCE(" another kluge for drift bc fac 2 ",)
		//D removing=2.0*facr*lbc.loss(x,qfle,n0,n1,n2, dmef, df); 
		D removing=facr*lbc.loss(x,qfle,n0,n1,n2, dmef, df); 
		//ndiff(i)-=removing;
// make it a leftover
//		ndiff(i)=ndiff(i)-removing;
		return removing;
	}
	if (near_left&&going_right)
	{
//		const D sum=(nearer_left)?(adding+adding2):adding3; 
		D sum=facr*lbc.gain(x,qfle,sol(0,in),sol(1,in),nl, dmef, df); 
		ndiff(i)+=sum;
		return -sum; 	
	}
	if (near_right&&going_left)
	{
		D adding =facr*rbc.gain(x,qfle,sol(rrp0,in),sol(rrp1,in),nr, dmef, df); 
		ndiff(i)+=adding;
		return -adding; 	
	}

return 0; 

}

/*
void itor_drift_bc(MyBlock & ndiff, Di & di, const D & dxmax , const D & xmin, const D & xmax, const IdxTy sz, const D & ish, const IdxTy in )

*/
/*

void do_drift_subpixel(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, const int iv)
*/


///////////////////////////////////////

void do_drift_one(MyBlock & ndiff, const MyBlock & sol, const D & t,  const int in, const int iv)
{
	CounterMap perf;
	const bool use_perf=m_tree.track_drift_delta();
	const bool keep_the_lost=!true;
	const bool add_from_beyond=!true;
	const bool add_from_beyond_e=true;
	const StrTy pname="perf";
	//typedef mjm_numb Ck;
	ndiff.resize(m_points);
	const D n0=sol(0,in);
	const int sz=m_points-1;
	const D nd=sol(sz,in);
	for (int  i=1; i<sz;  ++i)
	{
		const D n=sol(i,in);
		const D x=m_fixed(i,m_x);
		const D vdrift= sol(i,iv);
		const D dx=vdrift*t;
		const D target=dx+x;
			// while originally suspeicious, this is probably safe but al the checks
			// were to find a typoe now fixed. 
			const int idx=floor(m_fixed, sz, m_x, target, i);
			int idx1=idx;
			int idx2=idx+1;
			const bool ok=(idx2<m_points)&&(idx!=~0);
			if (ok) {
				//MM_MSG(" moving by "<<i<<" to "<<idx1<<" "<<(idx1-i)<<" dx "<<(m_fixed(idx1,m_x)-x))
				const D delx=m_fixed(idx1,m_x);
				const D delx2=m_fixed(idx2,m_x);
				const D f=(target-delx)/(delx2-delx);
				const D foo=ndiff(idx1);
				ndiff(idx1)+=n*(1.0-f);
				ndiff(idx2)+=n*f;
				if (add_from_beyond_e)
				{
					int di=idx2-i;
					int di1=idx1-i;
					int die=(i-idx1);
					int die2=(i-idx2);
					// in essence this is idx-2*i or "did they move fatther than i is close to edge"
					if (di>=i) { ndiff(i)+=f*n0; }
					if (di1>=i) { ndiff(i)+=(1.0-f)*n0; }
					// this looks wrong now for minority carriers on right 
					if ((die2)>=(sz-i)) { ndiff(i)+=f*nd; }
					if ((die)>=(sz-i)) { ndiff(i)+=(1.0-f)*nd; }

				}
				if (use_perf) perf.add(pname,idx1-i,1.0-f);
				if (use_perf) perf.add(pname,idx2-i,f);
			}
			else 
			{ // these oterhwise deplete the contacts of carriers 
				MM_MSG(" not ok move "<<i<<" dx "<<dx<<" "<<idx)
				if (keep_the_lost) ndiff(i)+=n;
				// the order here is reversed , this loop main part finds the source of all carriers
				// endings up a given pixel. It misses those off stage. 
				if (add_from_beyond)
				{


				} 
			}
	// it was never set 
	//	ndiff(i)-=n;
	} // i 
		D vdrift= sol(1,iv);
	if (add_from_beyond&&(vdrift>0))
	{ // this is where the stupid density and h get fed up
		//if (vdrift<=0) break;
		const D n=sol(0,in);
		IdxTy i=1;
		D tdrift=(m_fixed(i,m_x)-m_fixed(i-1,m_x))/vdrift;
//		if (vdrift
// 		 while transit time from bc is less then time step, add more carriers to each
// 	near-edge pixel
		while (tdrift<t)
		{
				const D target=m_fixed(i,m_x);
//				const D delx=m_fixed(idx1,m_x);
//				const D delx2=m_fixed(idx2,m_x);
				const D f=1.0; // (target-delx)/(delx2-delx);
			ndiff(i)+=n*(1.0-f);
			if (i>=sz) break; 
//			ndiff(i+1)+=n*(f);
			++i;
			vdrift= sol(i,iv);
			if (vdrift<=0) break;
			tdrift+=(m_fixed(i,m_x)-m_fixed(i-1,m_x))/vdrift;
		}

	}
		 vdrift= sol(sz,iv);
	if (add_from_beyond&&(vdrift<0))
	{ // this is where the stupid density and h get fed up
		//if (vdrift<=0) break;
		const D n=sol(sz,in);
		IdxTy i=sz;
		D tdrift=(m_fixed(i-1,m_x)-m_fixed(i,m_x))/vdrift;
//		if (vdrift
// 		 while transit time from bc is less then time step, add more carriers to each
// 	near-edge pixel
		while (tdrift<t)
		{
				const D target=m_fixed(i,m_x);
//				const D delx=m_fixed(idx1,m_x);
//				const D delx2=m_fixed(idx2,m_x);
				const D f=1.0; // (target-delx)/(delx2-delx);
			ndiff(i)+=n*(1.0-f);
			if (i==0) break; 
//			ndiff(i+1)+=n*(f);
			--i;
			vdrift= sol(i,iv);
			if (vdrift>=0) break;
			tdrift+=(m_fixed(i-1,m_x)-m_fixed(i,m_x))/vdrift;
		}
	}




	ndiff(0)=sol(0,in);
	ndiff(sz)=sol(sz,in);
	if (use_perf) { MM_MSG(" drift perf "<<perf.array_to_string(pname)) }
}



///////////////////////////////////////
 void initial_guess(   MyBlock & sol) {
//MM_ERR(" right initial_guess callled doh ")
//sol.resize(size());
sol.resize(m_points,m_vars);
/*
// available if needed 
typedef mjm_interpolants Miu;
//typedef mjm_exp_interp Mi;
*/
const IdxTy iu=m_u;
const IdxTy ip=m_p;
const IdxTy in=m_n;
const IdxTy iz=m_z;
const IdxTy pmax=m_points-0;
for (IdxTy point=0; point<pmax; ++point)
{
	const D fac=1.0*point/(pmax-1);
	if (point<m_mid) { sol(point,iu)= m_v0; } else { sol(point,iu)= m_vd; }
	if (point<m_mid) { sol(point,in)= m_n_p; sol(point,ip)=m_p_p; } 
	else { sol(point,in)= m_n_n; sol(point,ip)=m_p_n; }
	sol(point,iz)=1;
}
	if (m_tree.enable_trapping())
	m_defects.net_charge(m_state,m_trap,m_points);
}

void SetCols()
{
m_fixed.col_name("x",m_x);
m_fixed.col_name("Na",m_na);
m_fixed.col_name("Nd",m_nd);


m_state.col_name("voltage",m_u);
m_state.col_name("z",m_z);
m_state.col_name("h",m_p);
m_state.col_name("n",m_n);
m_state.col_name("v_h",m_v_p);
m_state.col_name("v_e",m_v_n);
m_state.col_name("rho",m_q);

}
void InitBase()
{
 	m_u=0; m_z=1;m_p=2; m_n=3;m_v_p=4;m_v_n=5,m_q=6,m_trap=7;
	m_x=0; m_na=1; m_nd=2;
    //m_best_l2=-1;
    //m_last_l2=-1;
    m_vars=8;
    //m_points=m_size/m_vars;
    m_size=m_points*m_vars;
    m_iter=0;

} // init

D neutral(const D & N, bool plus)
{
const D a=.5*N;
const D b=.5*::sqrt(N*N+4.0*m_ni2);
if (plus) return a+b ;

return a-b;

}
void InitConsts()
{

	m_fixed.resize(m_points,3);
	m_h=m_flp.total_length()/D(m_points-1); // 1e-6; // good value
	m_x_junction=m_flp.total_length()/D(2.0); // 1e-6; // good value
	m_Vt=m_flp.Vt(); // .0259;
	m_qe=m_flp.qe(); // 1.6e-19/8.854e-14;
	m_ni2=m_flp.ni2(); // 1.5e10*1.5e10;
	m_mup=m_flp.mup(); m_mun=m_flp.mun(); m_muz=m_flp.muz();
	const IdxTy midd=0+(m_points>>1);
	// this needs to be eliminated 
	m_mid=midd;
	m_v0=m_flp.v0(); // 0;
	m_vd=m_flp.vd(); // .5-.1*.5;
	//m_vd=.7;
}

void InitGeometry(const MyBlock & x)
{
	if ( x.size()<3)
	{
		for ( IdxTy i=0; i<m_points; ++i) { m_fixed(i,m_x)=m_h*D(i);	} 
	}
	else
	{
		if (x.size()!=m_points)
		{
			MM_INC_MSG(m_cm,"point_mismatch");
			MM_ERR(" points counts do not agree, wtf "<<m_points<<" vs "<<x.size()) 
		}
		for ( IdxTy i=0; i<m_points; ++i) { m_fixed(i,m_x)=x(i);	} 

	}	

}

void InitFixeds()
{
	const D Nd=m_flp.Nd(); // f*1e15; // + 
	const D Na=m_flp.Na(); // f*1e14; // - 
	// doh, these need to differ from Nx for neutrality. 
 	m_n_n= Fwd(neutral(Nd,true )); // Fwd(1e15*f);
 	m_p_n=Fwd(m_ni2/Rev(m_n_n));
	m_p_p=Fwd(neutral(Na,true)); // Fwd(1e14*f);
	m_n_p=Fwd(m_ni2/Rev(m_p_p)); //1e14;
	//const D midx=
	IdxTy ilast=0;
	for ( IdxTy i=0; i<m_points; ++i)
	{
		if (m_fixed(i,m_x)>m_x_junction) {  m_fixed(i,m_nd)=Rev(Nd); m_fixed(i,m_na)=0;	} 
		else {ilast =i; m_mid=i+1;  m_fixed(i,m_na)=Rev(Na); m_fixed(i,m_nd)=0;	} 
	}
	if (m_tree.enable_trapping())
	{
		m_defects.set(1,m_points);
		m_defects.add(1e13);
m_defects.set_to_equ( m_n_n, m_p_n, m_ni2, 0,  Nd,m_points-1);
m_defects.set_to_equ( m_n_p, m_p_p, m_ni2, Na,  0,0);
MM_MSG("traps bc "<< MMPR(m_n_n)<<MMPR( m_p_n)<<MMPR(m_p_p)<<MMPR(m_n_p)<<MMPR(m_x_junction))
//void set_to_equ(const IdxTy first, const IdxTy end, const D & ni2, const D & Na, const D & Nd)
m_defects.set_to_equ(0, ilast+1, m_ni2, Na, 0);
m_defects.set_to_equ(ilast+1, m_points, m_ni2, 0,  Nd);
MM_MSG(MMPR(ilast));
// does nothing as it is reset when initial_guess is done 
	m_defects.net_charge(m_state,m_trap,m_points);

	}
}
void InitSol()
{
	if (!m_save_sol) initial_guess(m_state);
	else
	{
			if (!m_state.form_is(m_points,m_vars))
			{
				MM_ERR(" resizing m_state from"<<m_state.form().to_string())	
				m_state.resize(m_points,m_vars);
				// there needs to be a projection operation som where. 
				// as well as saved location info 
			}

	}
} // InitSol
void InitVerify()
{
// this needs to be part of banner and have more for dev crap 
MM_MSG( MMPR(m_n_n)<<MMPR( m_p_n)<<MMPR(m_p_p)<<MMPR(m_n_p)<<MMPR(m_x_junction))
// points monotonic, no nans or negatives etc. 
// perhaps even hmax/hmin diagnose and critique

}
void Init(const IdxTy points=0, const MyBlock & x= MyBlock(1) ) {

	if (points!=0) m_points=points;
	if (x.size()>3) m_points=x.size();
	InitBase();
	// the solution has not yet been resized 
	//SetCols(); // guess resizes trashing cols possible. 
	InitConsts();
	InitGeometry(x);
	InitFixeds();
	InitSol();
	SetCols(); // guess resizes trashing cols possible. 
	InitVerify();
}


D Fwd(const D & x) const { return (x); } 
D Rev(const D & x) const { return (x); } 


BraTy m_tree;
FlpTy m_flp; // right now this is being SET in Init not being used to set m_h
IdxTy m_size;
MyBlock m_state;
MyBlock m_fixed;
//MyBlock m_best_sol;
//MyBlock m_current;
//D m_best_l2,m_last_l2;
IdxTy m_iter;
IdxTy m_vars, m_points;
CounterMap m_cm;

// misc stuff
D m_Vt,m_h,m_qe;
//D m_ni2,m_tau,m_mup, m_muz, m_mun;
D m_ni2,m_mup, m_muz, m_mun;
D m_n_n, m_p_n,m_p_p,m_n_p,m_v0,m_vd,m_x_junction;
IdxTy m_mid;
//MyBlock m_Na, m_Nd;
int m_u,m_z,m_p,m_n,m_v_n,m_v_p,m_q,m_trap;
int m_na,m_nd,m_x;
bool m_save_sol;

Defects m_defects;
}; // ficklanc




/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////

#ifdef  TEST_FICK__


static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
int main(int argc,char **args)
{
int  m_ierr = PetscInitialize(&argc,&args,(char*)0,help);if (m_ierr) return m_ierr;
typedef ficklanc Myt;

//Myt x(2*2000/6);
//Myt x(2000);
Myt x(argc,args);

//x.solve();
x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

