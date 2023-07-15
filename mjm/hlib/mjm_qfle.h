#ifndef MJM_QFLE_H__
#define MJM_QFLE_H__

#include "mjm_globals.h"
//#include "mjm_shooting_logic.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
// improve resolution of mb integral in place of commplementary normal distro
#include "mjm_numeric_tricks.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_sparse_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_diffuse_iterators.h"
//#include "mjm_defect_levels.h"

#include <algorithm>
#include <signal.h>
#include <sstream>
//#define USING_PETSC
//#ifdef USING_PETSC
// deviate from shooting to matrix 
//#include "mjm_petsc_util.h"
//#include "mjm_petsc_fd_base.h"
//#endif
/*
Copied from mjm_shooting10.h

*/
// TODO FIXME mbflow is noisy for beta small, certainly 1e-8 is awful. 
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
template <class Tt > void make_sf(const Tt & that, const geo_set & g,const bool src_right,const bool dest_right, const bool uniform_src, const bool uniform_dest=false)
{
// gee, make a member of g? lol 
if (uniform_dest) { s1=0; s0=1; } 
else that.make_sf(s1,s0, g.m_xp0, g.m_xp1, g.m_xp2,(dest_right)?1:0);
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
,only_left_src(false), only_right_src(false)
,only_left_dest(false),only_right_dest(false),uniform_dest(false)
,do_taylor(false),taylor_flag(0)
 {}
drift_flags(const bool ns, const bool us, const bool db):
 no_self(ns), uniform_src(us), dbg(db)
,only_left_src(false), only_right_src(false)
,only_left_dest(false),only_right_dest(false),uniform_dest(false)
,do_taylor(false),taylor_flag(0)
 {}
// alt, the source is an element either left or right of source node
void xleft() { only_left_src=true; }
void xright() { only_right_src=true; }
// the destination is just a sink between xp0 and xp1
void xpleft() { only_left_dest=true; }
// as above but use xp1 to xp2
void xpright() { only_right_dest=true; }
bool do_taylor_if_small() const { return ((taylor_flag&1) !=0); }
void do_taylor_if_small(const bool x )  { bitrs(taylor_flag,1,x); }
void bitrs(IdxTy & f, const IdxTy bits, const bool rs)
{
if (rs) f|=bits;
else f&=~bits;
}
void lr() // from left src to right dest
{
only_left_src=true; only_right_src=!true;
only_left_dest=!true; only_right_dest=true;
}
void ll() // from left src to left dest
{
only_left_src=true; only_right_src=!true;
only_left_dest=true; only_right_dest=!true;
}
void rl() // from right src to left dest
{
only_left_src=!true; only_right_src=true;
only_left_dest=true; only_right_dest=!true;
}
void rr() // from right src to left dest
{
only_left_src=!true; only_right_src=true;
only_left_dest=!true; only_right_dest=true;
}





// but the original point was to get both shape functions for the sink,
// or just make the dest shape function uniform... 
bool no_self;
bool uniform_src;
bool dbg; 
bool only_left_src;
bool only_right_src;
bool only_left_dest;
bool only_right_dest;
bool uniform_dest;
bool do_taylor;
IdxTy taylor_flag;
}; // drift_flags


template <class Td>
D mb_sf(const Td & di , const D & beta, const D & v0, const D & tau
,const bool only_node, const drift_flags & df, const bool dbg=true )
{
//drift_flags df(false,true,dbg);
if (only_node) return mb_sf(di.xp0(),di.xp1(),di.xp2(),di.x0(),di.x1(),di.x2()
, 0, di.n1(), 0,beta,v0,tau,df);
return mb_sf(di.xp0(),di.xp1(),di.xp2(),di.x0(),di.x1(),di.x2()
, di.n0(), di.n1(), di.n2(),beta,v0,tau,df);
}



#if 0

// maxwell boltzmann shape function etc as with others from x to xp
D mb_sf_old(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2
	,  const D & beta , const D & v0, const D & tau
	,const drift_flags & flags) const
{
const bool db=flags.dbg;
const bool dbg_self=(false&&(xp0==x0)); // for testing uniform_src self term 
D vt=3e10; // This still has relevance but not reall ?
D res=0;
typedef geo_set Geo;
typedef seg_set Seg;
// this needs a different vt for left and right sides 
Geo geo(xp0,xp1, xp2,x0,x1, x2,vt);
Seg sego;
D A,B,C;
// this is the pdf normalization, NOT a carrier concentration
//const D N=::sqrt(beta/M_PI)/tau;
// this is an integral over V, not x 
const D N=::sqrt(beta/M_PI);
const bool lin_fit=true;
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
// these ALWAS have some overlap now 
//if (sego.start(geo,flags.no_self)) {} // continue;
if (geo.self&&flags.no_self) return true; // skip 
// make part of geo? 
if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,geo.xbase);
sego.make_sf(*this, geo,src_right,dest_right, flags.uniform_src);
if (dbg_self)
{
MM_MSG(" attempting to get non-zero drift for nigh v self ????")
MM_MSG(" geo "<<geo.to_string()<<" sego "<<sego.to_string() ) 
}
// these all need origin shift
const D c1=sego.c1;
const D c0=sego.c0;

const D s1=0;
const D s0=0;

const D vt0=v0*tau;
const D alpha=beta/tau/tau;
// all the origins are differnt lol 
const D Q3=A*c1; // zero for now
const D Q2=A*c0+B*c1;
const D Q1=B*c0+C*c1;
const D Q0=C*c0;
// x origins need to be consistent
const D xf=geo.xmaxs;
const D xpf=geo.xpmax;
const D xi=geo.xmins;
const D xpi=geo.xpmin;
const D alpha2=2.0*alpha;
const D D1=Q0+Q3/alpha2;
const D D2=Q1/alpha2;
const D D3=Q2/alpha2;

//const D yf=xf-x'+v0t

const D kf=xf-vt0;
const D ki=xi-vt0;
const D zff=xpf+kf; // there are really 4 of these doh 
const D zfi=xpf+ki; // there are really 4 of these doh 
const D zif=xpi+kf; // there are really 4 of these doh 
const D zii=xpi+ki; // there are really 4 of these doh 
// g(v)=sqrt(beta/pi)*exp(-beta*v*v)
// beta=m/(2kT)
// exp(-beta*v0^2)s(x')c(x)n(x)exp(-beta*P(x'))exp(-beta*y^2)
// y=x-x'+V0*tau so limits change but dx=dy
// apparently erf(x) is the integral from x=0 to x and prefactor is 2/sqrt(pi)


//http://www.cplusplus.com/reference/cmath/erf/
// http://nvlpubs.nist.gov/nistpubs/jres/73B/jresv73Bn1p1_A1b.pdf
 // now just do the intergral of the normal v distro P(x)exp(-beta*(x-a)^2)
// and then the x' integral 
const D sqalpha=sqrt(alpha);
static const D ps=sqrt(M_PI);
const D E2=D2*(exp(-alpha*kf*kf)-exp(-alpha*ki*ki));
// wtf???
const D I2a=0; // E2*s1*((-exp(-alpha*zf*zf)+exp(-alpha*zi*zi) -exp(-alpha*zf*zf)+exp(-alpha*zi*zi) ))/2.0/alpha);
const D s0t=s0+0; // translation adds up 
const D I2b=E2*s0t*.5*sqrt(M_PI/alpha)*(
erf(sqalpha*zff)-erf(sqalpha*zif) +erf(sqalpha*zii)-erf(sqalpha*zfi) );

const D I3=D3*(xf-xi)*(.5*s1*xpf*xpf-.5*s1*xpi*xpi+s0*(xpf-xpi));
const D i1fac=s1*D1*.5*ps/sqalpha;
const D i1fac0=s0*D1*.5*ps/sqalpha;
const D I1a=i1fac*(Izerf(zff,sqalpha)-Izerf(zif,sqalpha)-Izerf(zfi,sqalpha)+Izerf(zii,sqalpha));
const D I1b=i1fac0*(Ierf(zff,sqalpha)-Ierf(zif,sqalpha)-Ierf(zfi,sqalpha)+Ierf(zii,sqalpha));
const D I1c=0;
const D I1d=0; 

const D segr=I2a+I2b+I3+I1a+I1b+I1c+I1d;
res+=N*segr;
MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(I2a)<<MMPR(I2b)<<MMPR(I3)<<MMPR(I1a)<<MMPR(I1b)<<MMPR(I1c)<<MMPR(I1d))
// the x integral with limtis generated erf's an expoentials
} // seg

return res;

}
#endif

// wht wolfrom integrals fixed some things, still negative
// numbers with larger beta and noise but back drift
// not obvious now. 
// integral of z*erf(az) from see nist ref   4.1.4
// http://nvlpubs.nist.gov/nistpubs/jres/73B/jresv73Bn1p1_A1b.pdf
D Izerf(const D & z, const D & a) const
{
static D ps=::sqrt(M_PI);
return (.5*z*z-.25/a/a)*erf(a*z)+.5*z/a/ps*exp(-a*a*z*z);
}
// integral of z*erf(az) from see nist ref   4.1.1
D Ierf(const D & z, const D & a) const
{
// http://www.wolframalpha.com/input/?i=integrate+erf%28a*x%29
static D ps=::sqrt(M_PI);
return z*erf(a*z)+1.0/a/ps*exp(-a*a*z*z);
}
// 4.1.8 in nist ref 
D Izzerf(const D & z, const D & a) const
{
static D ps=::sqrt(M_PI);
static D third=1.0/3.0;
//return .5*third/(a*a)*Ierf(z,a)+(z*z-1.0/(a*a))*z*third*erf(a*z)+1.0/(a*ps*3.0)*z*z*exp(-a*a*z*z);
return third*( (1.0+a*a*z*z)*exp(-z*z*a*a)/(ps*a*a*a)+z*z*z*erf(a*z));
}

// http://www.wolframalpha.com/input/?i=integrate+x^3*erf%28a*x%29
D Izzzerf(const D & z, const D & a) const
{
static D ps=::sqrt(M_PI);
const D z3=z*z*z;
const D z4=z3*z;
const D a3=a*a*a;
//return .5*third/(a*a)*Ierf(z,a)+(z*z-1.0/(a*a))*z*third*erf(a*z)+1.0/(a*ps*3.0)*z*z*exp(-a*a*z*z);
return (exp(-a*a*z*z)*(4.0*a3*z3+6.0*a*z)+ps*erf(a*z)*(4.0*z4*a3*a-3.0))/(16*ps*a3*a);
}


// as above, but etf(az)*z^(n+2)
D Iznerf(const D & z, const D & a, const IdxTy _n ) const
{
MM_ONCE(" probably wrong ",)
if (_n==2) return Izzerf(z,a);
if (_n==1) return Izerf(z,a);
if (_n==0) return Ierf(z,a);
const D n=_n-2;
static D ps=::sqrt(M_PI);
static D n3=1.0/(n+3.0);
const D zn1=::pow(z,n+1);
return .5*(n+2.0)*(n+1.0)*n3/(a*a)*Iznerf(z,a,_n-2)+(z*z-.5*(n+2.0)/(a*a))*zn1*n3*erf(a*z)+n3/(a*ps)*z*zn1*exp(-a*a*z*z);
}

// new analysis Ki is integral of z^i erf(az) 
//D K3(const D & z, const D & a) const { return Iznerf(z,a,3); } 
D K3(const D & z, const D & a) const { return Izzzerf(z,a); } 
D K2(const D & z, const D & a) const { return Izzerf(z,a); } 
D K1(const D & z, const D & a) const { return Izerf(z,a); } 
D K0(const D & z, const D & a) const { return Ierf(z,a); } 

//const D j0k0=(K0(yff,aeff)-K0(yfi,aeff))*c0i0f 
D K2K1K0(const D& z, const D& a, const D& c2, const D & c1, const D& c0) const
{
//
static D ps=::sqrt(M_PI);
static D third=1.0/3.0;
const D az=a*z;
const D az2=az*az;
const D cerf=erf(az);
const D cexp=exp(-az2)/(ps*a);
const D k0=c0*( z*cerf+cexp);
const D k1=c1*((.5*z*z-.25/a/a)*cerf+.5*z*cexp);
const D k2=c2*third*((1.0+az2)*cexp/(a*a)+z*z*z*cerf);
return k0+k1+k2;
}

D K2K1K0(const D& yf, const D& yi,const D& a, const D& c2, const D & c1, const D& c0) const
{
//const D ref= K2K1K0(yf,a,c2,c1,c0)- K2K1K0(yi,a,c2,c1,c0);
static D ps=::sqrt(M_PI);
static D third=1.0/3.0;
const D azi=a*yi;
const D azi2=azi*azi;
const D azf=a*yf;
const D azf2=azf*azf;

const D cerfi=erf(azi);
const D cerff=erf(azf);
const D cexpi=exp(-azi2)/(ps*a);
const D cexpf=exp(-azf2)/(ps*a);

const D k0=c0*( yf*cerff+cexpf-cexpi-yi*cerfi);
const D k1=c1*(.5*yf*yf*cerff-.5*yi*yi*cerfi- .25/a/a*(cerff-cerfi)+.5*(yf*cexpf-yi*cexpi));
const D k2=c2*third*(((1.0+azf2)*cexpf-(1.0+azi2)*cexpi) /(a*a)
		+yf*yf*yf*cerff-yi*yi*yi*cerfi);



return k0+k1+k2;
//return ref;

}


// TODO need to check a and a^2 and sqrt(a) FUDD 
// integrals of exp(-a*z^2)*z^n 
// http://www.wolframalpha.com/input/?i=integrate+x^2+exp%28-x*x%29
D Ke2(const D & z, const D & a) const 
// G and R disagree? 
//{const static D ps=sqrt(M_PI);   return -.5/(a*a)*(z-.5*ps/(a)*erf((a)*z)); } 
{const static D ps=sqrt(M_PI);   return -.5/(a*a)*(z*exp(-a*a*z*z)-.5*ps/(a)*erf((a)*z)); } 
D Ke1(const D & z, const D & a) const 
{   return -.5/(a*a)*exp(-a*a*z*z); } 
D Ke0(const D & z, const D & a) const 
{const static D ps=sqrt(M_PI);   return (.5*ps/(a)*erf((a)*z)); } 

D Ke2Ke1Ke0(const D& z, const D& a, const D& c2, const D & c1, const D& c0) const
{
static D ps=::sqrt(M_PI);
const D az=a*z;
const D az2=az*az;
const D cexp=exp(-az2);
const D cerf=.5*ps*erf(az);
//  -.5/(a*a)*(z*exp(-a*a*z*z)-.5*ps/(a)*erf((a)*z));   // Ke2
const D k2=c2*(-.5/(a*a)*(z*cexp-cerf/a));   // Ke2
// -.5/(a*a)*exp(-a*a*z*z); // Ke1 
const D k1=c1*( -.5/(a*a)*cexp); // Ke1 
//{(.5*ps/(a)*erf((a)*z)); // Ke0 
const D k0=c0*((cerf/a)); // Ke0 
return k2+k1+k0;


}
D Ke2Ke1Ke0(const D& yf, const D& yi,const D& a, const D& c2, const D & c1, const D& c0) const
{
//const D ref= Ke2Ke1Ke0(yf,a,c2,c1,c0)- Ke2Ke1Ke0(yi,a,c2,c1,c0);

static D ps=::sqrt(M_PI);
const D azi=a*yi;
const D azf=a*yf;
const D azf2=azf*azf;
const D azi2=azi*azi;
const D cexpf=exp(-azf2);
const D cexpi=exp(-azi2);
const D cerff=.5*ps*erf(azf);
const D cerfi=.5*ps*erf(azi);
//  -.5/(a*a)*(z*exp(-a*a*z*z)-.5*ps/(a)*erf((a)*z));   // Ke2
const D k2=c2*(-.5/(a*a)*(yf*cexpf-(cerff-cerfi)/a-yi*cexpi));   // Ke2
// -.5/(a*a)*exp(-a*a*z*z); // Ke1 
const D k1=c1*( -.5/(a*a)*(cexpf-cexpi)); // Ke1 
//{(.5*ps/(a)*erf((a)*z)); // Ke0 
const D k0=c0*(((cerff-cerfi)/a)); // Ke0 
return k2+k1+k0;
//return ref;
}
// doh, here "a" and "y" are at equal powers whereas
// the taylor expansion uses beta*y^2 fudd 
D K2K1K0Ke2Ke1Ke0(const D& yf, const D& yi,const D& a
, const D & c3, const D& c2, const D & c1, const D& c0
, const D& ce2, const D & ce1, const D& ce0) const
{

static D ps=::sqrt(M_PI);
static D third=1.0/3.0;
const D azi=a*yi;
const D azi2=azi*azi;
const D azf=a*yf;
const D azf2=azf*azf;

const D cerfi=erf(azi);
const D cerff=erf(azf);
const D cexpi=exp(-azi2); // /(ps*a);
const D cexpf=exp(-azf2); // /(ps*a);
const D expd=1.0/(ps*a);


//const D k2=c2*third*(((1.0+azf2)*cexpf-(1.0+azi2)*cexpi) /(a*a)
//		+yf*yf*yf*cerff-yi*yi*yi*cerfi);

//const D k0=c0*( yf*cerff+cexpf-cexpi-yi*cerfi);
const D k0=c0*( yf*cerff+cexpf*expd-cexpi*expd-yi*cerfi);
//const D k1=c1*(.5*yf*yf*cerff-.5*yi*yi*cerfi- .25/a/a*(cerff-cerfi)+.5*(yf*cexpf-yi*cexpi));
const D k1=c1*(.5*yf*yf*cerff-.5*yi*yi*cerfi- .25/a/a*(cerff-cerfi)+.5*expd*(yf*cexpf-yi*cexpi));
const D k2=c2*third*(((1.0+azf2)*cexpf*expd-(1.0+azi2)*cexpi*expd) /(a*a)
		+yf*yf*yf*cerff-yi*yi*yi*cerfi);

// k3
//return (exp(-a*a*z*z)*(4.0*a3*z3+6.0*a*z)+ps*erf(a*z)*(4.0*z4*a3*a-3.0))/(16*ps*a3*a);
D k3=0;
//if (c3!=0)
//{
const D  a3=a*a*a;
const D  a4=a3*a;
const D z3f=yf*yf*yf;
const D z3i=yi*yi*yi;
const D z4f=z3f*yf;
const D z4i=z3i*yi;
// http://www.wolframalpha.com/input/?i=integrate+x^3+erf%28a*x%29
//k3= c3*(cexpf*(4.0*a3*z3f+6.0*azf)+ps*cerff*(4.0*z4f*a4-3.0))/(16*ps*a4);
//k3= k3-c3*(cexpi*(4.0*a3*z3i+6.0*azi)+ps*cerfi*(4.0*z4i*a4-3.0))/(16*ps*a4);
k3= c3*(4.0*a3*(z3f*cexpf-z3i*cexpi)+6.0*(cexpf*azf-cexpi*azi)
	+ps*a4*4.0*(cerff*z4f-cerfi*z4i)-3.0*ps*(cerff-cerfi))/(16*ps*a4);

//k3=c3/(16.0*ps*a4)*(6.0*a*(cexpf*yf-cexpi*yi)+4.0*a3*(cexpf*z3f-cexpi*z3i) +ps*4.0*a4*(cerff*z4f-cerfi*z4i)-3.0*ps*(cerff-cerfi));

// c3*exp*(4*a3*z3+6*az)/(16*ps*a4)
// c3*erf*(ps*a4*z4*4-3*ps)/(16*ps*a4)
//}

//static D ps=::sqrt(M_PI);
//const D azi=a*yi;
//const D azf=a*yf;
//const D azf2=azf*azf;
//const D azi2=azi*azi;

// moved up and change facto r
//const D cexpf=exp(-azf2);
//const D cexpi=exp(-azi2);
//const D cerff=erf(azf); // *ps*.5;
//const D cerfi=erf(azi); // *ps*.5;
const D erfn=.5*ps;
//  -.5/(a*a)*(z*exp(-a*a*z*z)-.5*ps/(a)*erf((a)*z));   // Ke2
//const D ke2=ce2*(-.5/(a*a)*(yf*cexpf-(cerff*yf-cerfi*yi)*erfn/a-yi*cexpi));   // Ke2
// http://www.wolframalpha.com/input/?i=integrate+x^2+exp%28-x*x*a*a%29
const D ke2=ce2*(-.5/(a*a)*(yf*cexpf-(cerff-cerfi)*erfn/a-yi*cexpi));   // Ke2
// -.5/(a*a)*exp(-a*a*z*z); // Ke1 
const D ke1=ce1*( -.5/(a*a)*(cexpf-cexpi)); // Ke1 
//{(.5*ps/(a)*erf((a)*z)); // Ke0 
const D ke0=ce0*(((cerff-cerfi)*erfn/a)); // Ke0 

const D sh1=c2*third*expd;
//cexpf*(c0*expd+.5*c1*expd*yf+c2*third*(1.0+azf2)*expd/(a*a)-.5/(a*a)*(ce2*yf+ce1))
//cexpf*(c0*expd+.5*c1*expd*yf+sh1*yf*yf+sh1/(a*a)-.5/(a*a)*(ce2*yf+ce1))
// these look similar now but now sure about noise floor 
// c3*exp*(4*a3*z3+6*az)/(16*ps*a4)
// c3*erf*(ps*a4*z4*4-3*ps)/(16*ps*a4)
const D c3s=c3/(16*ps*a4);

const D cgf=(4.0*a3*z3f*c3s+c3s*6.0*azf+ c0*expd+.5*c1*expd*yf+0*sh1*yf*yf-.5/(a*a)*(ce2*yf+ce1-2.0*sh1*(1.0+azf2))); 
const D cgfx3=( c0*expd+.5*c1*expd*yf+0*sh1*yf*yf-.5/(a*a)*(ce2*yf+ce1-2.0*sh1*(1.0+azf2))); 

const D altxf=cexpf*cgf;
const D cgi=(4.0*a3*z3i*c3s+c3s*6.0*azi+c0*expd+.5*c1*expd*yi+0*sh1*yi*yi-.5/(a*a)*(ce2*yi+ce1-2.0*sh1*(1.0+azi2))); 
const D cgix3=(c0*expd+.5*c1*expd*yi+0*sh1*yi*yi-.5/(a*a)*(ce2*yi+ce1-2.0*sh1*(1.0+azi2))); 

const D altxi=cexpi*cgi;
IdxTy nnorm=0;
if (a<1e-3) nnorm=10;
D n3fac=.5*(azf+azi);
if (n3fac==0) n3fac=1;
//D trickedc3s=6.0*c3s*mjm_numeric_tricks::sqdiffe( -azf2,-azi2,azf/n3fac,azi/n3fac,20)*n3fac;
D trickedc3s=6.0*c3s*mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,azf/n3fac,azi/n3fac,nnorm)*n3fac;
D n3fac2=.5*(z3f+z3i);
 trickedc3s+=4.0*c3s*a3*mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,z3f/n3fac2,z3i/n3fac2,nnorm)*n3fac2;


D nx3fac=.5*(cgfx3+cgix3);
if (nx3fac==0) nx3fac=1;
//D trickedxc3s=mjm_numeric_tricks::sqdiffe( -azf2,-azi2,cgfx3/nx3fac,cgix3/nx3fac,20)*nx3fac;
D trickedxc3s=mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,cgfx3/nx3fac,cgix3/nx3fac,nnorm)*nx3fac;
D talt=trickedxc3s+trickedc3s;
D nfac=.5*(cgf+cgi);
if (nfac==0) nfac=1.0;
IdxTy niter=0;
//if ( cexpf>.9999) if (cexpi>.9999) niter=5; 
if ( a<1e-3) niter=10; 
// not mattering 
//D tricked=mjm_numeric_tricks::sqdiffe( -azf2,-azi2,cgf/nfac,cgi/nfac,niter)*nfac;
D tricked=mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,cgf/nfac,cgi/nfac,niter)*nfac;

//D sqdiffetaylor(const D & a, const D & b, const D & ca, const D & cb, const IdxTy n)


//const D altyf= cerff*(c3s*ps*a4*z4f*4.0-c3s*3.0*ps+c0*yf+.5*c1*yf*yf-.25*c1/(a*a)+c2*third*yf*yf*yf+ce2*.5/(a*a*a)*erfn+ce0*erfn/a) ;
const D cegf= (c3s*ps*a4*z4f*4.0-c3s*3.0*ps+c0*yf+.5*c1*yf*yf-.25*c1/(a*a)+c2*third*yf*yf*yf+ce2*.5/(a*a*a)*erfn+ce0*erfn/a) ;
const D cegfx3= (c3s*ps*a4*z4f*4.0+.5*c1*yf*yf) ;
const D altyf= cerff*cegf ;

//const D altyi= cerfi*(c3s*ps*a4*z4i*4.0-c3s*3.0*ps+c0*yi+.5*c1*yi*yi-.25*c1/(a*a)+c2*third*yi*yi*yi+ce2*.5/(a*a*a)*erfn+ce0*erfn/a) ;
const D cegi=(c3s*ps*a4*z4i*4.0-c3s*3.0*ps+c0*yi+.5*c1*yi*yi-.25*c1/(a*a)+c2*third*yi*yi*yi+ce2*.5/(a*a*a)*erfn+ce0*erfn/a) ;
const D cegix3=(c3s*ps*a4*z4i*4.0+.5*c1*yi*yi) ;

const D altyi= cerfi*cegi;
const D nerffac=.5*(cegix3+cegfx3);
IdxTy nitererf=0;
//if ( a<1e-3) if( azf2<1) if (azi2<1)  nitererf=20; 
if ( a<1e-3)  nitererf=10; 
//D trickederf=mjm_numeric_tricks::sqdifferftaylor( azf,azi,cegf/nerffac,cegi/nerffac,nitererf)*nerffac;
D trickederf=0;
trickederf+=mjm_numeric_tricks::sqdifferftaylor( azf,azi,cegfx3/nerffac,cegix3/nerffac,nitererf)*nerffac ;
const D f31=.5*((yi+yf));
trickederf+= mjm_numeric_tricks::sqdifferftaylor( azf,azi,yf/f31,yi/f31,nitererf)*(c0*f31) ;
const D f33=.5*((yi*yi*yi+yf*yf*yf));
trickederf+= mjm_numeric_tricks::sqdifferftaylor( azf,azi,yf*yf*yf/f33,yi*yi*yi/f33,nitererf)*(c2*f33*third) ;
trickederf+= mjm_numeric_tricks::sqdifferftaylor( azf,azi,1.0,1.0,nitererf)*(-.25*c1/(a*a)-c3s*3.0*ps) ;
trickederf+= mjm_numeric_tricks::sqdifferftaylor( azf,azi,1.0,1.0,nitererf)*(ce0*erfn/a) ;
trickederf+= mjm_numeric_tricks::sqdifferftaylor( azf,azi,1.0,1.0,nitererf)*(.5*ce2*erfn/(a*a*a)) ;

//+ mjm_numeric_tricks::sqdifferftaylor( azf,azi,1.0,1.0,nitererf)*(-c3s*3.0*ps);
// more important
//const D c3s=c3/(16*ps*a4);
const D altyifd= 4.0*c3s*ps*a4*(cerff*z4f-cerfi*z4i)
+c0*(cerff*yf-cerfi*yi)
+.5*c1*(cerff*yf*yf-cerfi*yi*yi)
-.25*c1/(a*a)*(cerff-cerfi)
+c2*third*(cerff*yf*yf*yf-cerfi*yi*yi*yi)
+(ce2*.5/(a*a*a)*erfn -c3s*3.0*ps)*(cerff-cerfi)
+ce0*erfn/a*(cerff-cerfi) ;


const D da= altxf-altxi+altyf-altyi;
const D daold= k3+k2+k1+k0+ke2+ke1+ke0;
if (false)
{
MM_MSG(MMPR(altxf)<<MMPR(altxi)<<MMPR(altyf)<<MMPR(altyi)<<MMPR(da)<<MMPR(daold)<<MMPR(altxf-altxi)<<MMPR(altxf-altxi-talt)<<MMPR(altyf-altyi)<<MMPR(altxf+altyf)<<MMPR(altyf-altyi-altyifd))
MM_MSG(MMPR(cexpf)<<MMPR(cexpi)<<MMPR(cerff)<<MMPR(cerfi)<<MMPR(a*yf)<<MMPR(a*yi))
}

if (true) return altxf-altxi+altyf-altyi;
if (!true) return tricked+altyf-altyi;
if (!true) return tricked+altyifd;
if (!true) return tricked+trickederf;
// not working 
if (!true) return talt+trickederf;
//if (true) return talt+altyifd;

//if (false) 
return k3+k2+k1+k0+ke2+ke1+ke0;

}



///////////////////////////////////////////////////////////////////

#if 0 
D mb_sf_old(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2
	,  const D & beta , const D & v0, const D & tau
	,const drift_flags & flags) const
{
const bool db=flags.dbg;
const bool dbg_self=(false&&(xp0==x0)); // for testing uniform_src self term 
D vt=3e10; // This still has relevance but not reall ?
D res=0;
typedef geo_set Geo;
typedef seg_set Seg;
// this needs a different vt for left and right sides 
Geo geo(xp0,xp1, xp2,x0,x1, x2,vt);
Seg sego;
D A,B,C;
// this is the pdf normalization, NOT a carrier concentration
//const D N=::sqrt(beta/M_PI)/tau;
// this is an integral over V, not x 
// TODO FIXME the tau does not belong here 
const D N=::sqrt(beta/M_PI)/tau;
static const  D ps=::sqrt(M_PI);
const bool lin_fit=true;
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
// these ALWAS have some overlap now 
//if (sego.start(geo,flags.no_self)) {} // continue;
if (geo.self&&flags.no_self) return true; // skip 
// make part of geo? 
if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,geo.xbase);
sego.make_sf(*this, geo,src_right,dest_right, flags.uniform_src);
if (dbg_self)
{
MM_MSG(" attempting to get non-zero drift for nigh v self ????")
MM_MSG(" geo "<<geo.to_string()<<" sego "<<sego.to_string() ) 
}
// these all need origin shift
const D xbase=geo.xbase;
const D xpbase=geo.xpbase;
const D c1=sego.c1;
const D c0=sego.c0;

const D s1=sego.s1;
// for nos make this relative to simpler origin 
// this would make it absolate, removing xpbase
//const D s0=sego.s0-xpbase*s1;

const D s0=sego.s0-xpbase*s1+xbase*s1;
//const D s0=sego.s0+xpbase*s1-xbase*s1;

// FIXME there is a fudding sign mistake in here some where 
const D vt0=-v0*tau;
const D alpha=beta/tau/tau;
const D aeff=::sqrt(alpha);
// do any xbase translate here... 
const D a=B;
// the erf and exp only depend on differences anyway
// keep this xbase. 
//const D b=C-a*xbase;
const D b=C; // -a*xbase;
const D c=c1;
const D d=c0;

//const D xf=geo.xmaxs;
//const D xi=geo.xmins;
//const D xpf=geo.xpmax;
//const D xpi=geo.xpmin;
// keep these numbers smalll 
const D xf=geo.xmaxs-xbase;
const D xi=geo.xmins-xbase;
MM_MSG(" carrieres nf="<<(a*xf+b)<<" ni="<<(a*xi+b))
// for now these work best with same base 
const D xpf=geo.xpmax-xbase;
const D xpi=geo.xpmin-xbase;


// copied from old code, check these 
const D kf=xf-vt0;
const D ki=xi-vt0;

//const D zff=-xpf+kf; // there are really 4 of these doh 
//const D zfi=-xpf+ki; // there are really 4 of these doh 
//const D zif=-xpi+kf; // there are really 4 of these doh 
//const D zii=-xpi+ki; // there are really 4 of these doh 
/*

// this is based on a new analysis  
// the alpha square root leel may be fudded up 

const D z1fj2=s1*a*c;
const D z1ij2=s1*a*c;
const D z0fj2=a*c*(s0-s1*kf);
const D z0ij2=a*c*(s0-s1*ki);

//{ 
const D J2ff1=z1fj2*Ke1(zff,aeff);
const D J2ff0=z0fj2*Ke0(zff,aeff);
const D J2ii1=z1ij2*Ke1(zii,aeff);
const D J2ii0=z0ij2*Ke0(zii,aeff);

const D J2if1=z1ij2*Ke1(zfi,aeff);
const D J2if0=z0ij2*Ke0(zfi,aeff);
const D J2fi1=z1fj2*Ke1(zif,aeff);
const D J2fi0=z0fj2*Ke0(zif,aeff);

//}
const D J2= J2ff1+J2ff0+J2ii1+J2ii0
			-(J2if1+J2if0+J2fi1+J2fi0); 
//MM_MSG(MMPR(J2)<<MMPR(J2ff1)<<MMPR(J2ff0)<<MMPR(J2ii1)<<MMPR(J2ii0)<<MMPR(J2if1)<<MMPR(J2if0)<<MMPR(J2fi1)<<MMPR(J2fi0))
//const D J2=z1fj2*Ke1(zff,aeff)+ z0fj2*Ke0(zff,aeff)
//		+z1ij2*Ke1(zii,aeff)+ z0ij2*Ke0(zii,aeff)
//	-(	+z1ij2*Ke1(zfi,aeff)+ z0ij2*Ke0(zfi,aeff)
//		+z1fj2*Ke1(zif,aeff)+ z0fj2*Ke0(zif,aeff));


const D R1=a*d+b*c-2.0*a*c*vt0;
const D z2fj1=s1*a*c/aeff;
const D z2ij1=s1*a*c/aeff;
const D z1fj1=-.5/aeff*(s1*R1+s1*4.0*a*c*kf+2.0*s0*a*c);
const D z1ij1=-.5/aeff*(s1*R1+s1*4.0*a*c*ki+2.0*s0*a*c);
const D z0fj1=.5/aeff*(s1*kf*R1+s0*R1+2.0*a*c*s1*kf*kf+2.0*a*c*kf*s0);
const D z0ij1=.5/aeff*(s1*ki*R1+s0*R1+2.0*a*c*s1*ki*ki+2.0*a*c*ki*s0);

const D J1=z2fj1*Ke2(zff,aeff)+z1fj1*Ke1(zff,aeff)+z0fj1*Ke0(zff,aeff)
+z2ij1*Ke2(zii,aeff)+z1ij1*Ke1(zii,aeff)+z0ij1*Ke0(zii,aeff)
-(+z2ij1*Ke2(zfi,aeff)+z1ij1*Ke1(zfi,aeff)+z0ij1*Ke0(zfi,aeff)
+z2fj1*Ke2(zif,aeff)+z1fj1*Ke1(zif,aeff)+z0fj1*Ke0(zif,aeff));
*/




//const D w2=a*c;
//const D w1=-2.0*vt0+a*c+a*d+b*c;
//const D w0=b*d+a*c*vt0*vt0-vt0*(a*d+b*c);
//const D z2f=.5*ps/aeff*w2;
//const D z2i=z2f;
//const D z1f=.5*ps/aeff*(-w2-2.0*w2*kf);
//const D z1i=.5*ps/aeff*(-w2-2.0*w2*ki);
//const D z0f=w0+w1*kf+w2*kf*kf;
//const D z0i=w0+w1*ki+w2*ki*ki;


// when speed matters, note that some integrals are duplicated
// now... 
//const D K2ff=K2(zff,aeff);
//const D K2ii=K2(zii,aeff);

// this is actually NEGATIVE and needs to be multiplied by s0
// this needs to be var changed with y=k-x' and so it depends 
// on k and x'... doh 
//const D s0ff=s0+kf*s1;
//const D s0fi=s0+ki*s1;
//const D s0if=s0+kf*s1;
//const D s0ii=s0+ki*s1;
const D syi=s0+ki*s1;
const D syf=s0+kf*s1;

// use c1=c, c0=d
const D cy3f=s1*a*c;
const D cy3i=s1*a*c;
const D cy2f= -s1*a*(d+c*xf)-s1*c*(b+a*xf)-a*c*syf;
const D cy2i= -s1*a*(d+c*xi)-s1*c*(b+a*xi)-a*c*syi;
const D cy1f=s1*(b+a*xf)*(d+c*xf)+a*syf*(d+c*xf)+c*syf*(b+a*xf);
// TODO verify fixing the syi did not wreck anythign lol 
const D cy1i=s1*(b+a*xi)*(d+c*xi)+a*syi*(d+c*xi)+c*syi*(b+a*xi);
const D cy0f=-syf*(b+a*xf)*(d+c*xf);
const D cy0i=-syi*(b+a*xi)*(d+c*xi);
const D jf0=ps*.5/aeff;

const D yff=xf-xpf-vt0;
const D yii=xi-xpi-vt0;
const D yif=xi-xpf-vt0;
const D yfi=xf-xpi-vt0;

// this appears to work 
const D j0k3=(K3(yff,aeff)-K3(yfi,aeff))*cy3f + (K3(yii,aeff)-K3(yif,aeff))*cy3i;
const D j0k2=(K2(yff,aeff)-K2(yfi,aeff))*cy2f + (K2(yii,aeff)-K2(yif,aeff))*cy2i;
const D j0k1=(K1(yff,aeff)-K1(yfi,aeff))*cy1f + (K1(yii,aeff)-K1(yif,aeff))*cy1i;
const D j0k0=(K0(yff,aeff)-K0(yfi,aeff))*cy0f + (K0(yii,aeff)-K0(yif,aeff))*cy0i;
const D J0alt=(j0k3+j0k2+j0k1+j0k0)*jf0;
MM_MSG(MMPR(J0alt)<<MMPR(j0k0)<<MMPR(j0k1)<<MMPR(j0k2)<<MMPR(j0k3)<<MMPR(jf0))

const D jf2=.5*a*c/(aeff*aeff);
const D f2=.5*ps/aeff;
const D cy3f2=-s1/3.0;
const D cy3i2=-s1/3.0;
const D cy2f2=syf*.5;
const D cy2i2=syi*.5;
const D cy1f2=-f2*s1;
const D cy1i2=-f2*s1;
const D cy0f2=f2*syf;
const D cy0i2=f2*syi;

D J2alt = 
( yff*yff*yff-yfi*yfi*yfi)*cy3f2 +(yii*yii*yii-yif*yif*yif)*cy3i2;
J2alt += 
( yff*yff-yfi*yfi)*cy2f2 +(yii*yii-yif*yif)*cy2i2;
J2alt += 
(K1(yff,aeff)-K1(yfi,aeff))*cy1f2 + (K1(yii,aeff)-K1(yif,aeff))*cy1i2;
J2alt +=
(K0(yff,aeff)-K0(yfi,aeff))*cy0f2 + (K0(yii,aeff)-K0(yif,aeff))*cy0i2;
J2alt=J2alt*jf2;



const D jf1=.5/(aeff*aeff);
//const D ca1=a+a*c;
const D ca1=c*a+a*c;
const D cjf=c*b+a*d+ca1*xf;
const D cji=c*b+a*d+ca1*xi;
const D cy2f1=s1*ca1;
const D cy2i1=s1*ca1;
const D cy1f1=-s1*cjf-syf*ca1;
const D cy1i1=-s1*cji-syi*ca1;
const D cy0f1=cjf*syf;
const D cy0i1=cji*syi;

//D J1alt = 
const D j1k2=(Ke2(yff,aeff)-Ke2(yfi,aeff))*cy2f1 + (Ke2(yii,aeff)-Ke2(yif,aeff))*cy2i1;
//J1alt += 
const D j1k1=(Ke1(yff,aeff)-Ke1(yfi,aeff))*cy1f1 + (Ke1(yii,aeff)-Ke1(yif,aeff))*cy1i1;
//J1alt +=
const D j1k0=(Ke0(yff,aeff)-Ke0(yfi,aeff))*cy0f1 + (Ke0(yii,aeff)-Ke0(yif,aeff))*cy0i1;
//J1alt=J1alt*jf1;
const D J1alt=(j1k2+j1k1+j1k0)*jf1;
MM_MSG(MMPR(J1alt)<<MMPR(j1k2)<<MMPR(j1k1)<<MMPR(j1k0)<<MMPR(jf1)) 

/*
const D J00=s0ff*(K2(zff,aeff)*z2f+K1(zff,aeff)*z1f+K0(zff,aeff)*z0f)
	+s0ii*(K2(zii,aeff)*z2i+K1(zii,aeff)*z1i+K0(zii,aeff)*z0i)
	-(s0if*(+K2(zfi,aeff)*z2i+K1(zfi,aeff)*z1i+K0(zfi,aeff)*z0i)
	+s0fi*(K2(zif,aeff)*z2f+K1(zif,aeff)*z1f+K0(zif,aeff)*z0f));
const D J01=K3(zff,aeff)*z2f+K2(zff,aeff)*z1f+K1(zff,aeff)*z0f
	+K3(zii,aeff)*z2i+K2(zii,aeff)*z1i+K1(zii,aeff)*z0i
	-(+K3(zfi,aeff)*z2i+K2(zfi,aeff)*z1i+K1(zfi,aeff)*z0i
	+K3(zif,aeff)*z2f+K2(zif,aeff)*z1f+K1(zif,aeff)*z0f);
// the x' in s needs to be shifted into integration variable y...
// that has not been done yet 
*/

//const D J0= J0alt; // -J00 + s1*J01;


const D segr=J0alt+J1alt+J2alt; // I2a+I2b+I3+I1a+I1b+I1c+I1d;
res+=N*segr;
//MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(I2a)<<MMPR(I2b)<<MMPR(I3)<<MMPR(I1a)<<MMPR(I1b)<<MMPR(I1c)<<MMPR(I1d))
//MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0alt)<<MMPR(J0)<<MMPR(J00)<<MMPR(J01)<<MMPR(s0)<<MMPR(s1))

if (!false)
{
MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0alt)<<MMPR(J1alt)<<MMPR(J2alt)<<MMPR(s0)<<MMPR(s1))
}
//MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0)<<MMPR(J1)<<MMPR(J2)<<MMPR(0)<<MMPR(0)<<MMPR(0)<<MMPR(0))
// the x integral with limtis generated erf's an expoentials
} // seg

return res;

}
#endif

// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ old 
/////////////////////////////////////////////////////////////////
// VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV new 

/*
2017-04-29 The integrals K0-3 and Ke0-3 rely on various sources
that did not all seem to work. 
A copy of Gradstein and Rhyhik , nist ref link , and wolfram site
The nist paper has nice formulas but not sure I got them to work,

http://nvlpubs.nist.gov/nistpubs/jres/73B/jresv73Bn1p1_A1b.pdf
Right now the results of the online integrator appear to 
work great,

http://www.wolframalpha.com/input/?i=integrate+x^2+exp%28-x*x%29

but need to find numerical noise as beta -> zero. 

the integrand is P(x)exp(-beta*(x-y-v*t)^2)
*/

D mb_sf_old(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2
	,  const D & beta , const D & v0, const D & tau
	,const drift_flags & flags) const
{
const bool db=flags.dbg;
const bool dbg_self=(false&&(xp0==x0)); // for testing uniform_src self term 
const bool dbdump=false;
D vt=3e10; // This still has relevance but not reall ?
D res=0;
typedef geo_set Geo;
typedef seg_set Seg;
// this needs a different vt for left and right sides 
Geo geo(xp0,xp1, xp2,x0,x1, x2,vt);
Seg sego;
D A,B,C;
// this is the pdf normalization, NOT a carrier concentration
//const D N=::sqrt(beta/M_PI)/tau;
// this is an integral over V, not x 
// TODO FIXME the tau does not belong here 
//const D N=::sqrt(beta/M_PI)/tau;
// the caller was wrong
//const D N=::sqrt(beta*beta/M_PI);
const D N=::sqrt(beta/M_PI);
static const  D ps=::sqrt(M_PI);
const bool lin_fit=true;
for (IdxTy seg=0; seg<4; ++seg)
{
// right now this is kind of a bastard between
// nodes and elements. The elements sink but the
// nodes are the sources ... I guess a seg filter woud work  
const bool  src_right=((seg&1)!=0);
if (flags.only_left_src&&src_right) continue;
if (flags.only_right_src&&!src_right) continue;

const bool  dest_right=((seg&2)!=0);
if (flags.only_left_dest&&dest_right) continue;
if (flags.only_right_dest&&!dest_right) continue;




// xbase is zero for "x" and is  nearest boundary x0 or x2
geo.set_sd(src_right,dest_right);
// create ulp and llp check for overlap 
if (flags.no_self) { MM_ONCE(" no_self set so ignore self diffusion",) } 
// these ALWAS have some overlap now 
// WTF this converted a bool to a D... 
//if (sego.start(geo,flags.no_self)) {} // continue;
//if (geo.self&&flags.no_self) return true; // skip 
if (geo.self&&flags.no_self) continue; // skip 
// make part of geo? 
if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,geo.xbase);
sego.make_sf(*this, geo,src_right,dest_right, flags.uniform_src,flags.uniform_dest);
if (dbg_self)
{
MM_MSG(" attempting to get non-zero drift for nigh v self ????")
MM_MSG(" geo "<<geo.to_string()<<" sego "<<sego.to_string() ) 
}
// these all need origin shift
const D xbase=geo.xbase;
const D xpbase=geo.xpbase;
const D c1=sego.c1;
const D c0=sego.c0;

const D s1=sego.s1;

const D s0=sego.s0-xpbase*s1+xbase*s1;
//const D s0=sego.s0+xpbase*s1-xbase*s1;

// FIXME there is a fudding sign mistake in here some where 
const D vt0=-v0*tau;
// this should have been fixed, done by caller 
//const D alpha=beta/tau/tau;
const D alpha=beta*beta;
//const D aeff=::sqrt(alpha);
const D aeff=::sqrt(beta);
// do any xbase translate here... 
const D a=B;
// the erf and exp only depend on differences anyway
// keep this xbase. 
//const D b=C-a*xbase;
const D b=C; // -a*xbase;
const D c=c1;
const D d=c0;

// keep these numbers smalll 
const D xf=geo.xmaxs-xbase;
const D xi=geo.xmins-xbase;
//MM_MSG(" carrieres nf="<<(a*xf+b)<<" ni="<<(a*xi+b))
// for now these work best with same base 
const D xpf=geo.xpmax-xbase;
const D xpi=geo.xpmin-xbase;
// copied from old code, check these 

const D yff=xf-xpf-vt0;
const D yii=xi-xpi-vt0;
const D yif=xi-xpf-vt0;
const D yfi=xf-xpi-vt0;

const D syf=-s0-s1*(xf-vt0);
const D syi=-s0-s1*(xi-vt0);

// fix later lol 
const D ac=a*c; 
const D ca=a*c; 

const D cy1f=(-2.0*xf*ac-d*a-b*c);
const D cy1i=(-2.0*xi*ac-d*a-b*c);
const D cy0f=(d*b+ac*xf*xf+d*a*xf+b*c*xf);
const D cy0i=(d*b+ac*xi*xi+d*a*xi+b*c*xi);

const D jf0=ps*.5/aeff;
const D c0i3f=s1*ac;
const D c0i3i=s1*ac;
const D c0i2f=ac*syf+s1*cy1f;
const D c0i2i=ac*syi+s1*cy1i;
const D c0i1f=s1*(cy0f)+syf*(cy1f);
const D c0i1i=s1*(cy0i)+syi*(cy1i);
const D c0i0f=syf*cy0f;
const D c0i0i=syi*cy0i;

const bool do_old_int=false;
const bool do_alt_int=!false; 
D J0alt=0;
if (do_old_int)
{
// this appears to work 
const D j0k3=(ac==0)?0:((K3(yff,aeff)-K3(yfi,aeff))*c0i3f 
			+ (K3(yii,aeff)-K3(yif,aeff))*c0i3i);
//D K2K1K0( y, a,  c2,  c1,  c0) const
// K2K1K0(yff,aeff,c0i2f,c0i1f,c0i0f);
// K2K1K0(yfi,aeff,c0i2f,c0i1f,c0i0f);
// K2K1K0(yii,aeff,c0i2i,c0i1i,c0i0i);
// K2K1K0(yif,aeff,c0i2i,c0i1i,c0i0i);
// K2K1K0(yff,yfi,aeff,c0i2f,c0i1f,c0i0f);
// K2K1K0(yii,yif,aeff,c0i2i,c0i1i,c0i0i);
const bool j0new=true;
const D j0newi= (j0new)?(K2K1K0(yff,yfi,aeff,c0i2f,c0i1f,c0i0f)+ K2K1K0(yii,yif,aeff,c0i2i,c0i1i,c0i0i)):0;

const D j0k2=(j0new)?0:((K2(yff,aeff)-K2(yfi,aeff))*c0i2f 
			+ (K2(yii,aeff)-K2(yif,aeff))*c0i2i);
const D j0k1=(j0new)?0:((K1(yff,aeff)-K1(yfi,aeff))*c0i1f 
			+ (K1(yii,aeff)-K1(yif,aeff))*c0i1i);
const D j0k0=(j0new)?0:((K0(yff,aeff)-K0(yfi,aeff))*c0i0f 
			+ (K0(yii,aeff)-K0(yif,aeff))*c0i0i);


//const D 
J0alt=(j0k3+j0k2+j0k1+j0k0+j0newi)*jf0;
//MM_MSG(MMPR(c0i3f)<<MMPR(c0i3i)<<MMPR(c0i2f)<<MMPR(c0i2i)) // <<MMPR(j0k3)<<MMPR(jf0))
if (dbdump) {
MM_MSG(MMPR(c0i1f)<<MMPR(c0i1i)<<MMPR(c0i0f)<<MMPR(c0i0i)<<MMPR(j0newi))
MM_MSG(MMPR(K1(yff,aeff))<<MMPR(K1(yfi,aeff))<<MMPR(K0(yff,aeff))<<MMPR(K0(yfi,aeff)))
MM_MSG(MMPR(J0alt)<<MMPR(j0k0)<<MMPR(j0k1)<<MMPR(j0k2)<<MMPR(j0k3)<<MMPR(jf0))
}

} // do_old_int


const D jf1=-.5/(aeff*aeff);
const D c1i2f=s1*(-2.0*ca);
const D c1i2i=s1*(-2.0*ca);
const D c1i1f=s1*(-cy1f)+syf*(-2.0*ca);
const D c1i1i=s1*(-cy1i)+syi*(-2.0*ca);
const D c1i0f=syf*(-cy1f);
const D c1i0i=syi*(-cy1i);

const bool j1new=true;

D J1alt=0;
if (do_old_int)
{
if (j1new) { 
const D j1newi= Ke2Ke1Ke0(yff,yfi,aeff,c1i2f,c1i1f,c1i0f)+ Ke2Ke1Ke0(yii,yif,aeff,c1i2i,c1i1i,c1i0i);

J1alt=j1newi*jf1;
if (dbdump) {
MM_MSG(MMPR(J1alt)<<MMPR(j1newi)<<MMPR(jf1)) 
}

}
else {
const D j1k2=(ca==0)?0:((Ke2(yff,aeff)-Ke2(yfi,aeff))*c1i2f 
			+ (Ke2(yii,aeff)-Ke2(yif,aeff))*c1i2i);
const D j1k1=(Ke1(yff,aeff)-Ke1(yfi,aeff))*c1i1f 
			+ (Ke1(yii,aeff)-Ke1(yif,aeff))*c1i1i;
const D j1k0=(Ke0(yff,aeff)-Ke0(yfi,aeff))*c1i0f 
			+ (Ke0(yii,aeff)-Ke0(yif,aeff))*c1i0i;
 J1alt=(j1k2+j1k1+j1k0)*jf1;
if (dbdump) {
MM_MSG(MMPR(J1alt)<<MMPR(j1k2)<<MMPR(j1k1)<<MMPR(j1k0)<<MMPR(jf1)) 
}
}
} // do_old_int



// old code, only needed  for non uniform src shape function  
D J2alt=0;
const D jf2=1; // .5*a*c/(aeff*aeff);
const D a3=(1.0/(aeff*aeff*aeff));
const D a2=(1.0/(aeff*aeff));
const D c2i3f=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i3i=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i2f=-syf*ac*.5*a2; // /aeff/aeff;
const D c2i2i=-syi*ac*.5*a2; // /aeff/aeff;
const D c2i1f=.25*ps*a3*s1*ac;
const D c2i1i=.25*ps*a3*s1*ac;
const D c2i0f=.25*ps*a3*syf*ac;
const D c2i0i=.25*ps*a3*syi*ac;

if ( do_old_int) if (ac!=0) {
const D j2k3= (Ke2(yff,aeff)-Ke2(yfi,aeff))*c2i3f 
			+ (Ke2(yii,aeff)-Ke2(yif,aeff))*c2i3i;
const D j2k2= (Ke1(yff,aeff)-Ke1(yfi,aeff))*c2i2f 
			+ (Ke1(yii,aeff)-Ke1(yif,aeff))*c2i2i;
const D j2k1= (K1(yff,aeff)-K1(yfi,aeff))*c2i1f 
			+ (K1(yii,aeff)-K1(yif,aeff))*c2i1i;
const D j2k0= (K0(yff,aeff)-K0(yfi,aeff))*c2i0f 
			+ (K0(yii,aeff)-K0(yif,aeff))*c2i0i;

J2alt=(j2k0+j2k1+j2k2+j2k3)*jf2;
if (dbdump) {
MM_MSG(MMPR(J2alt)<<MMPR(j2k3)<<MMPR(j2k2)<<MMPR(j2k1)<<MMPR(j2k0)<<MMPR(jf1)) 
}
}

if (do_alt_int)
{
//const D j0k3ai=(ac==0)?0:((K3(yff,aeff)-K3(yfi,aeff))*c0i3f 
//			+ (K3(yii,aeff)-K3(yif,aeff))*c0i3i);
D  alt1=K2K1K0Ke2Ke1Ke0(yff,  yfi, aeff ,  c0i3f*jf0,c0i2f*jf0,  c0i1f*jf0+c2i1f*jf2,  c0i0f*jf0+c2i0f*jf2 ,
	  c1i2f*jf1+c2i3f*jf2,  c1i1f*jf1+c2i2f*jf2,  c1i0f*jf1);
D  alt2= K2K1K0Ke2Ke1Ke0(yii,  yif, aeff ,c0i3i*jf0,  c0i2i*jf0,  c0i1i*jf0+c2i1i*jf2,  c0i0i*jf0+c2i0i*jf2 ,
	  c1i2i*jf1+c2i3i*jf2,  c1i1i*jf1+c2i2i*jf2,  c1i0i*jf1);
J0alt=0;J1alt=alt1+alt2 ; // +j0k3ai*jf0;
J2alt=0;
if (dbdump) {
MM_MSG(MMPR(J1alt)<<MMPR(alt1)<<MMPR(alt2))
}
}
const D segr=J0alt+J1alt+J2alt; 
const D ab=beta;
const D area=ab*(yff*yff+yii*yii-yif*yif-yfi*yfi);
const D fom1=ab*yff*yff;
const D fom2=ab*yii*yii;
//MM_MSG(" fom "<<MMPR(alpha)<<MMPR(area)<<MMPR(yff)<<MMPR(yii)<<MMPR(fom1)<<MMPR(fom2))
const bool do_full_taylor=((flags.do_taylor)
		||(flags.do_taylor_if_small()&&(fom1<1)&&fom2<1)); // true;
if (do_full_taylor)
{
if (false) MM_MSG(" fom "<<MMPR(alpha)<<MMPR(area)<<MMPR(yff)<<MMPR(yii)<<MMPR(fom1)<<MMPR(fom2))
 //trickedc3s+=4.0*c3s*a3*mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,z3f/n3fac2,z3i/n3fac2,nnorm)*n3fac2;
// poly coefs with origin tranlsated to x'+vt0 (sign?)
const D nterms=10;
// the normalization is different here...
const D jft=ps*.5/aeff;
//const D beta=alpha; // aeff;
// these need the polynomial origin to match the exp arg 
//const D c0x=c0;
//const D s0x=s0;
//const D bx=b;
// except the s transaltion is limite dependent 
//const D cxzed3=-s1*c1*a;
//const D cxzed2=s0x*c1*a+s1*c0x*a+s1*c1*bx;
//const D cxzed1=s1*c0x*bx+s0x*c1*bx+s0x*c0x*a;
//const D cxzed0=s0x*c0x*bx;
//do the math lol 
//const D yif=xi-xpf-vt0;
const D p2=c1*a;
const D p1f=c1*b+a*c0+2.0*xf*c1*a;//  FUDD split up the y part of p1... 
const D p1i=c1*b+a*c0+2.0*xi*c1*a;
const D p1y=-2.0*c1*a;
const D cab=c1*b+c0*a;
const D p0f=c0*b+cab*xf+c1*a*xf*xf;
const D p0i=c0*b+cab*xi+c1*a*xi*xi;
const D p0yf=-cab-2.0*xf*c1*a;
const D p0yi=-cab-2.0*xi*c1*a;
const D p0y2=c1*a;
// this really depends on x^2 somehwat independe tof beta the way the
// cdoe is 
IdxTy niter=100;
if (beta<.1) niter=10;
bool failed=false;
 D segrft=jft*mjm_numeric_tricks::mb_taylor
	//(beta,yff,yii, yfi, yif,s1,syf,syi,p2,p1,p0,10);
	(beta,yff,yii,yif,yfi,s1,syf,syi,p2,p1f,p1i,p1y,p0f,p0i,p0yf,p0yi,p0y2,niter, failed);
//static D mb_taylor(const D & beta,const D & yff,const D & yii,
//    const D & yfi,const D & yif,const D & s1, const D & s0f,const D & s0i
//,const D & p2, const D & p1, const D & p0,  const IdxTy nterms)

//D ii3= cxzed3*mjm_numeric_tricks::mb_taylor(beta,yff,yii,yfi,yif,3,nterms);
//D ii2= cxzed2*mjm_numeric_tricks::mb_taylor(beta,yff,yii,yfi,yif,2,nterms);
//D ii1= cxzed1*mjm_numeric_tricks::mb_taylor(beta,yff,yii,yfi,yif,1,nterms);
//D ii0= cxzed0*mjm_numeric_tricks::mb_taylor(beta,yff,yii,yfi,yif,0,nterms);
// D segrft=jft*(ii3+ii2+ii1+ii0);
if (segrft<0) { MM_MSG(" danger will robinson segrft < 0 "<<MMPR(segrft)) }
//if ( beta<1e-1) 
const D scaled=2.0*segrft/M_PI*beta; // *N*N; // the 2.0 is a kluge factor something wrong 
if (!failed) res+=scaled;
else { res+=N*segr;}
//if (!false) MM_MSG(MMPR(beta)<<MMPR(segrft)<<MMPR(scaled)<<MMPR(segr)<<MMPR(segr*N))
const D fige=(segr*N-scaled)/(segr*N+scaled); // 
if ((fige>.1)||(fige<-.1)) MM_MSG("tayloer mismatch "<<MMPR(fige)<<MMPR(beta)<<MMPR(segrft)<<MMPR(scaled)<<MMPR(segr)<<MMPR(segr*N))
//else res+=N*segr;
}

else { res+=N*segr;}

if (dbdump)
{
MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0alt)<<MMPR(J1alt)<<MMPR(J2alt)<<MMPR(s0)<<MMPR(s1))
}
//MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0)<<MMPR(J1)<<MMPR(J2)<<MMPR(0)<<MMPR(0)<<MMPR(0)<<MMPR(0))
// the x integral with limtis generated erf's an expoentials
} // seg

return res;

}
///////////////////////////////////////////////////////////////////////
 // new one

D mb_sf(const D & xp0, const D & xp1, const D & xp2, const D & x0, const D & x1, const D & x2
	, const D & n0, const D & n1, const D & n2
	,  const D & beta , const D & v0, const D & tau
	,const drift_flags & flags) const
{
const bool db=flags.dbg;
const bool dbg_self=(false&&(xp0==x0)); // for testing uniform_src self term 
const bool dbdump=false;
D vt=3e10; // This still has relevance but not reall ?
D res=0;
typedef geo_set Geo;
typedef seg_set Seg;
// this needs a different vt for left and right sides 
Geo geo(xp0,xp1, xp2,x0,x1, x2,vt);
Seg sego;
D A,B,C;
// this is the pdf normalization, NOT a carrier concentration
//const D N=::sqrt(beta/M_PI)/tau;
// this is an integral over V, not x 
// TODO FIXME the tau does not belong here 
//const D N=::sqrt(beta/M_PI)/tau;
// the caller was wrong
//const D N=::sqrt(beta*beta/M_PI);
const D N=::sqrt(beta/M_PI);
static const  D ps=::sqrt(M_PI);
const bool lin_fit=true;
for (IdxTy seg=0; seg<4; ++seg)
{
// right now this is kind of a bastard between
// nodes and elements. The elements sink but the
// nodes are the sources ... I guess a seg filter woud work  
const bool  src_right=((seg&1)!=0);
if (flags.only_left_src&&src_right) continue;
if (flags.only_right_src&&!src_right) continue;

const bool  dest_right=((seg&2)!=0);
if (flags.only_left_dest&&dest_right) continue;
if (flags.only_right_dest&&!dest_right) continue;

// xbase is zero for "x" and is  nearest boundary x0 or x2
geo.set_sd(src_right,dest_right);
// create ulp and llp check for overlap 
if (flags.no_self) { MM_ONCE(" no_self set so ignore self diffusion",) } 
// these ALWAS have some overlap now 
// WTF this converted a bool to a D... 
//if (sego.start(geo,flags.no_self)) {} // continue;
//if (geo.self&&flags.no_self) return true; // skip 
if (geo.self&&flags.no_self) continue; // skip 
// make part of geo? 
if (lin_fit) fitlinbase(A,B,C,x0,x1,x2,n0,n1,n2,src_right,geo.xbase);
sego.make_sf(*this, geo,src_right,dest_right, flags.uniform_src,flags.uniform_dest);
if (dbg_self)
{
MM_MSG(" attempting to get non-zero drift for nigh v self ????")
MM_MSG(" geo "<<geo.to_string()<<" sego "<<sego.to_string() ) 
}
// these all need origin shift
const D xbase=geo.xbase;
const D xpbase=geo.xpbase;
const D c1=sego.c1;
const D c0=sego.c0;

const D s1=sego.s1;

const D s0=sego.s0-xpbase*s1+xbase*s1;
//const D s0=sego.s0+xpbase*s1-xbase*s1;

// FIXME there is a fudding sign mistake in here some where 
const D vt0=-v0*tau;
// this should have been fixed, done by caller 
//const D alpha=beta/tau/tau;
const D alpha=beta*beta;
//const D aeff=::sqrt(alpha);
const D aeff=::sqrt(beta);
// do any xbase translate here... 
const D a=B;
// the erf and exp only depend on differences anyway
// keep this xbase. 
//const D b=C-a*xbase;
const D b=C; // -a*xbase;
const D c=c1;
const D d=c0;

// keep these numbers smalll 
const D xf=geo.xmaxs-xbase;
const D xi=geo.xmins-xbase;
//MM_MSG(" carrieres nf="<<(a*xf+b)<<" ni="<<(a*xi+b))
// for now these work best with same base 
const D xpf=geo.xpmax-xbase;
const D xpi=geo.xpmin-xbase;
// copied from old code, check these 

const D yff=xf-xpf-vt0;
const D yii=xi-xpi-vt0;
const D yif=xi-xpf-vt0;
const D yfi=xf-xpi-vt0;

const D syf=-s0-s1*(xf-vt0);
const D syi=-s0-s1*(xi-vt0);

// fix later lol 
const D ac=a*c; 
const D ca=a*c; 

const D cy1f=(-2.0*xf*ac-d*a-b*c);
const D cy1i=(-2.0*xi*ac-d*a-b*c);
const D cy0f=(d*b+ac*xf*xf+d*a*xf+b*c*xf);
const D cy0i=(d*b+ac*xi*xi+d*a*xi+b*c*xi);

const D jf0=ps*.5/aeff;
const D c0i3f=s1*ac;
const D c0i3i=s1*ac;
const D c0i2f=ac*syf+s1*cy1f;
const D c0i2i=ac*syi+s1*cy1i;
const D c0i1f=s1*(cy0f)+syf*(cy1f);
const D c0i1i=s1*(cy0i)+syi*(cy1i);
const D c0i0f=syf*cy0f;
const D c0i0i=syi*cy0i;

//const bool do_old_int=false;
//const bool do_alt_int=!false; 
D J0alt=0;

const D jf1=-.5/(aeff*aeff);
const D c1i2f=s1*(-2.0*ca);
const D c1i2i=s1*(-2.0*ca);
const D c1i1f=s1*(-cy1f)+syf*(-2.0*ca);
const D c1i1i=s1*(-cy1i)+syi*(-2.0*ca);
const D c1i0f=syf*(-cy1f);
const D c1i0i=syi*(-cy1i);

const bool j1new=true;

D J1alt=0;


// old code, only needed  for non uniform src shape function  
D J2alt=0;
const D jf2=1; // .5*a*c/(aeff*aeff);
const D a3=(1.0/(aeff*aeff*aeff));
const D a2=(1.0/(aeff*aeff));
const D c2i3f=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i3i=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i2f=-syf*ac*.5*a2; // /aeff/aeff;
const D c2i2i=-syi*ac*.5*a2; // /aeff/aeff;
const D c2i1f=.25*ps*a3*s1*ac;
const D c2i1i=.25*ps*a3*s1*ac;
const D c2i0f=.25*ps*a3*syf*ac;
const D c2i0i=.25*ps*a3*syi*ac;

//if (do_alt_int)
{
//const D j0k3ai=(ac==0)?0:((K3(yff,aeff)-K3(yfi,aeff))*c0i3f 
//			+ (K3(yii,aeff)-K3(yif,aeff))*c0i3i);
D  alt1=K2K1K0Ke2Ke1Ke0(yff,  yfi, aeff ,  c0i3f*jf0,c0i2f*jf0,  c0i1f*jf0+c2i1f*jf2,  c0i0f*jf0+c2i0f*jf2 ,
	  c1i2f*jf1+c2i3f*jf2,  c1i1f*jf1+c2i2f*jf2,  c1i0f*jf1);
D  alt2= K2K1K0Ke2Ke1Ke0(yii,  yif, aeff ,c0i3i*jf0,  c0i2i*jf0,  c0i1i*jf0+c2i1i*jf2,  c0i0i*jf0+c2i0i*jf2 ,
	  c1i2i*jf1+c2i3i*jf2,  c1i1i*jf1+c2i2i*jf2,  c1i0i*jf1);
J0alt=0;J1alt=alt1+alt2 ; // +j0k3ai*jf0;
J2alt=0;
if (dbdump) {
MM_MSG(MMPR(J1alt)<<MMPR(alt1)<<MMPR(alt2))
}
}
const D segr=J0alt+J1alt+J2alt; 
const D ab=beta;
const D area=ab*(yff*yff+yii*yii-yif*yif-yfi*yfi);
const D fom1=ab*yff*yff;
const D fom2=ab*yii*yii;
//MM_MSG(" fom "<<MMPR(alpha)<<MMPR(area)<<MMPR(yff)<<MMPR(yii)<<MMPR(fom1)<<MMPR(fom2))
const bool do_full_taylor=((flags.do_taylor)
		||(flags.do_taylor_if_small()&&(fom1<1)&&fom2<1)); // true;
if (do_full_taylor)
{
if (false) MM_MSG(" fom "<<MMPR(alpha)<<MMPR(area)<<MMPR(yff)<<MMPR(yii)<<MMPR(fom1)<<MMPR(fom2))
 //trickedc3s+=4.0*c3s*a3*mjm_numeric_tricks::sqdiffetaylor( -azf2,-azi2,z3f/n3fac2,z3i/n3fac2,nnorm)*n3fac2;
// poly coefs with origin tranlsated to x'+vt0 (sign?)
//const D nterms=10;
// the normalization is different here...
const D jft=ps*.5/aeff;
//const D beta=alpha; // aeff;
//do the math lol 
//const D yif=xi-xpf-vt0;
const D p2=c1*a;
const D p1f=c1*b+a*c0+2.0*xf*c1*a;//  FUDD split up the y part of p1... 
const D p1i=c1*b+a*c0+2.0*xi*c1*a;
const D p1y=-2.0*c1*a;
const D cab=c1*b+c0*a;
const D p0f=c0*b+cab*xf+c1*a*xf*xf;
const D p0i=c0*b+cab*xi+c1*a*xi*xi;
const D p0yf=-cab-2.0*xf*c1*a;
const D p0yi=-cab-2.0*xi*c1*a;
const D p0y2=c1*a;
// this really depends on x^2 somehwat independe tof beta the way the
// cdoe is 
IdxTy niter=1;
if (beta<.1) niter=1;
bool failed=false;
 D segrft=jft*mjm_numeric_tricks::mb_taylor
	//(beta,yff,yii, yfi, yif,s1,syf,syi,p2,p1,p0,10);
	(beta,yff,yii,yif,yfi,s1,syf,syi,p2,p1f,p1i,p1y,p0f,p0i,p0yf,p0yi,p0y2,niter, failed);
//static D mb_taylor(const D & beta,const D & yff,const D & yii,
// D segrft=jft*(ii3+ii2+ii1+ii0);
if (segrft<0) { MM_MSG(" danger will robinson segrft < 0 "<<MMPR(segrft)) }
//if ( beta<1e-1) 
const D scaled=2.0*segrft/M_PI*beta; // *N*N; // the 2.0 is a kluge factor something wrong 
D picking=0;
D ascaled=N*segr;
if (ascaled>0) if ( scaled>ascaled) { MM_MSG(" assuming failure..."<<MMPR2(ascaled,scaled) ) failed=true; }
if (scaled<0) if ( ascaled>0) { MM_MSG(" assuming failure due to sign ..."<<MMPR2(ascaled,scaled) ) failed=true; }

if (!failed) picking=scaled;
else { picking=ascaled;}
res+=picking;
//if (!false) MM_MSG(MMPR(beta)<<MMPR(segrft)<<MMPR(scaled)<<MMPR(segr)<<MMPR(segr*N))
const D fige=(segr*N-scaled)/(segr*N+scaled); // 
if ((fige>.1)||(fige<-.1)) MM_MSG("tayloer mismatch "<<MMPR(fige)<<MMPR(beta)<<MMPR(segrft)<<MMPR(scaled)<<MMPR(picking)<<MMPR(segr)<<MMPR(segr*N))
//else res+=N*segr;
}

else { res+=N*segr;}

if (dbdump)
{
MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0alt)<<MMPR(J1alt)<<MMPR(J2alt)<<MMPR(s0)<<MMPR(s1))
}
//MM_MSG(MMPR(seg)<<MMPR(segr)<<MMPR(N)<<MMPR(J0)<<MMPR(J1)<<MMPR(J2)<<MMPR(0)<<MMPR(0)<<MMPR(0)<<MMPR(0))
// the x integral with limtis generated erf's an expoentials
} // seg

return res;

}
///////////////////////////////////////////////////////////////////////

























/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
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


#endif

