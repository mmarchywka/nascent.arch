#ifndef MJM_INTEGRALS_UTIL_H__
#define MJM_INTEGRALS_UTIL_H__

#include "mjm_globals.h"
//#include "mjm_geometry.h"
//#include "mjm_libmesh.h"
//#include "mjm_shapes.h"
//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"

/*
2017-10-15
 2408  g++ -DTEST_BURMANN__ -Wall -Wno-unused-function  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -gdwarf-3 -x c++ mjm_integrals.h 
 2409  ./a.out 10 0 1 0 > xxx
 2410  vi xxx
 2411  ../delbackup.tex 
 2412  svn commit -m " it seems to work matches all terms in (32) in http://www.mathematica-journal.com/2014/11/on-burmanns-theorem-and-its-application-to-problems-of-linear-and-nonlinear-heat-transfer-and-diffusion/#Eq:21"
 2413  cp mjm_integrals.h mjm_integrals_matches_eqn32in1.h 




g++ -DTEST_BURMANN__ -Wall  -Wno-unused-function  -std=gnu++11 -gdwarf-3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -x c++ mjm_integrals.h 
Interpolate exponentials on rhs, 
 g++ -DTEST_EXP_SF__ -Wall  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -x c++ mjm_integrals.h 

when I added the catalog feature to the block matrix it required the right libs, 
 g++ -DTEST_INTEGRAL_SF__ -Wall  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -x c++ mjm_integrals.h 


 g++ -DTEST_INTEGRAL__ -Wall  -std=gnu++11 -I.. -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -x c++ mjm_petsc_util.h 
 g++ -DTEST_INTEGRAL__ -Wall  -std=gnu++11 -I.. -x c++ mjm_integrals.h 
 g++ -DTEST_INTEGRAL_AS2__ -Wall  -std=gnu++11 -I.. -x c++ mjm_integrals.h 

 g++ -DTEST_PETSC_CLASS_UTIL_MAIN__ -Wall  -std=gnu++11 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -x c++ mjm_petsc_util.h 


*/

class mjm_integrals
{

typedef mjm_integrals Myt;

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
typedef mjm_generic_itor<2> LoopItor2;
typedef mjm_generic_itor<3> LoopItor3;
typedef mjm_generic_itor<4> LoopItor4;

mjm_integrals () {}

~mjm_integrals () {}


class exp_sf
{
// Integrals of the form P(y)exp(Q(y)) that come from integrating
// expoenentials of FEM quantifies interpolated with shape functions.
// Speciically, do the X integral first ( x/a ) and then Y (y/b).
// Q = cY+d  and P includes the shpe function Y component but
// also the results of the X integration in terms of K(y)= aY+b.
// This would be easy except when K=0,  easier to just use power expansion lol
// Integration by parts works when p(n)=0 or ||c||>1. Otherwsie
// integrate by parts in other directions for ||c||<1 which hopefully works ...
typedef exp_sf Myt;
typedef mjm_integrals Super;
typedef Super Tr;
typedef Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef mjm_generic_itor<2> LoopItor2;

public:


// try to evaulate approximately by parts
// see notes for meanings but there are variables a,b,c,d that
// are derived from the values . a and b also used for element
//lengths but here it is normlaized to 1x1.
// do the x integral first then y 
template <class Tv>
D int_parts(const Tv & values, const IdxTy si, const IdxTy nx, const IdxTy ny ) 
{
D I=0;
const D b=values[1]-values[0];
const D a=values[2]-values[3]-b;
const D c=values[3]-values[0];
const D d=values[0];
// depending on the ||c||, different schemes can work. There is also an issue
//  The c==0 singularity should be ok if the ||c||<1 thing works. 

return I;
}

// integrate by brute force trapezoid, \int_{0}^{1} dX dY S_iexp( \sigma S_j u_j)
// or the expenential of the interpolated value times one of the shape functions.
//  for now just use indexes instead of polynomials to get sf...
template <class Tv>
D int_trapezoid(const Tv & values, const IdxTy si, const IdxTy nx, const IdxTy ny ) 
{
D I=0;
const D a=1;
const D b=1;
const D dx=a/D(nx);
const D dy=b/D(ny);
const D cx=dx*.5;
const D cy=dy*.5;
LoopItor2 itor(nx,ny);
if (!itor) { MM_ERR(" itor should have at least one itor wtf "<<nx<<" "<<ny) } 
while (itor)
{
const D x=itor[0]*dx+cx;
const D y=itor[1]*dy+cy;
const D v=interpolate(values,x,y);
// exp likely slow lol 
I+=sf(si,x,y)*exp(v); 
++itor;
} // itor
return I*dx*dy; 
}
template <class Tv> 
D interpolate(const Tv & values , const D& x, const D& y) 
{
D v=0;
const IdxTy sz=values.size();
if (sz!=4) { MM_ERR(" value size should be 4 but is "<<sz) } 
for (IdxTy i=0; i<sz; ++i)
{ v+=sf(i,x,y)*values[i]; } // i 
return v;
}

D sf(const IdxTy i, const D& x, const D& y)
{
switch (i)
{
case 0:{ return (1.0-x)*(1.0-y); }
case 1:{ return (x)*(1.0-y); }
case 2:{ return (x)*(y); }
case 3:{ return (1.0-x)*(y); }

} // i
MM_ERR(" sf fall through should not happen out of range "<<i) 
return 0;

} // sf


}; // exp_sf
class as2
{
// Integral of form r^alpha*Sin(beta theta)
// Integrate over a rectangle from (0,0) to (a,b)
typedef as2 Myt;
typedef mjm_integrals Super;
typedef Super Tr;
typedef Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef mjm_generic_itor<2> LoopItor2;

public:

as2() : m_alpha(1),m_beta(1),m_xc(0),m_yc(0),m_rot(0) {}
as2(const D & alpha, const D & beta) 
	: m_alpha(alpha),m_beta(beta),m_xc(0),m_yc(0),m_rot(0){}
as2(const D & alpha, const D & beta, const D & xc, const D & yc, const D & rot) 
	: m_alpha(alpha),m_beta(beta),m_xc(xc),m_yc(yc),m_rot(rot){}


D integrand(const D & rho, const D & alpha, const D & theta, const D & beta)
{ return std::pow(rho,alpha)*std::sin(beta*theta); }
// the chain rule produces cos theta without beta
D integrand_dx(const D & rho, const D & alpha, const D & theta, const D & beta)
{ return alpha*std::pow(rho,alpha-1)*(std::sin(beta*theta)*std::cos(theta)
-std::cos(beta*theta)*std::sin(theta)) ; }
D integrand_dy(const D & rho, const D & alpha, const D & theta, const D & beta)
{ return alpha*std::pow(rho,alpha-1)*(std::sin(beta*theta)*std::sin(theta)
+std::cos(beta*theta)*std::cos(theta)) ; }


D integrand(const D & rho,const D & theta)
{ return integrand(rho, m_alpha, theta, m_beta); } 

D integrand_dx(const D & rho,const D & theta)
{ return integrand_dx(rho, m_alpha, theta, m_beta); } 


D integrand_dy(const D & rho,const D & theta)
{ return integrand_dy(rho, m_alpha, theta, m_beta); } 



// the surface is along x=0  so theta is 0 to pi/2
D int_trapezoid(const D & a, const D & b, const IdxTy nx, const IdxTy ny)
{
D I =0;
IdxTy n=0; 
const D dx=a/D(nx);
const D dy=a/D(ny);
const D cx=dx*.5;
const D cy=dy*.5;
LoopItor2 itor(nx,ny);
while (itor)
{
const D x=itor[0]*dx+cx;
const D y=itor[1]*dy+cy;
const D rho=sqrt(x*x+y*y);
const D theta=atan2(y,x);
I+=integrand(rho,theta); ++n;
++itor;
}
I=I*a*b/D(n);
return I;
}

// this is supposed to be the integral of the model equation for potential
// times a shape function. The laplacian is zero and all we need is the grad
// which is multipled by grad carrier concentration. 
D int_trapezoid(const D & a, const D & b, const IdxTy nx, const IdxTy ny
	, const IdxTy mod, const IdxTy shape, const IdxTy dshape)
{

D I =0;
IdxTy n=0; 
const D dx=a/D(nx);
const D dy=b/D(ny);
const D cx=dx*.5-m_xc;
const D cy=dy*.5-m_yc;
LoopItor2 itor(nx,ny);
const D adb=1.0/(a*b);
// use the bump feature and only redo sonme numbers after carry out
// or better use the other integration scheme lol 
while (itor)
{
const D x=itor[0]*dx+cx;
const D y=itor[1]*dy+cy;
const D rho=sqrt(x*x+y*y);
const D theta=atan2(y,x)-m_rot;
D imodel=1; // integrand(rho,theta); 
D sf=1;
D dsf=1;

switch (mod)
{
case 0: {  imodel=integrand(rho,theta);  break; }
case 1: {  imodel=integrand_dx(rho,theta);  break; }
case 2: {  imodel=integrand_dy(rho,theta);  break; }


default:  
	{ MM_ERR(" bad shape specification "<<mod) }

} // modl 
switch (shape)
{
case 0:{ sf=(a-x)*(b-y);  break; }
case 1:{ sf=(x)*(b-y);  break; }
case 2:{ sf=(x)*(y);  break; }
case 3:{ sf=(a-x)*(y);  break; }

default:  
	{ MM_ERR(" bad shape specification "<<shape) }
} // shape

switch (dshape)
{
case 0:{dsf=(-1.0)*(b-y);  break; }
case 1:{ dsf=(1.0)*(b-y);  break; }
case 2:{ dsf=(1.0)*(y);  break; }
case 3:{ dsf=(-1.0)*(y);  break; }

case 4:{ dsf=(a-x)*(-1.0);  break; }
case 5:{ dsf=(x)*(-1.0);  break; }
case 6:{ dsf=(x)*(1.0);  break; }
case 7:{ dsf=(a-x)*(1.0);  break; }

default:  
	{ MM_ERR(" bad shape specification "<<dshape) }
} // shape


I+=imodel*sf*dsf*adb;
// der shape
++n;
++itor;
}
I=I*a*b/D(n);
return I;



} // int_Trap


D int_trapezoid(const D & a, const D & b, const IdxTy nx, const IdxTy ny, const IdxTy shape)
{
D I =0;
IdxTy n=0; 
const D dx=a/D(nx);
const D dy=b/D(ny);
const D cx=dx*.5-m_xc;
const D cy=dy*.5-m_yc;
LoopItor2 itor(nx,ny);
const D adb=1.0/(a*b);
// use the bump feature and only redo sonme numbers after carry out
// or better use the other integration scheme lol 
while (itor)
{
const D x=itor[0]*dx+cx;
const D y=itor[1]*dy+cy;
const D rho=sqrt(x*x+y*y);
const D theta=atan2(y,x)-m_rot;
switch (shape)
{
case 0:{ I+=integrand(rho,theta)*(a-x)*(b-y)*adb;  break; }
case 1:{ I+=integrand(rho,theta)*(x)*(b-y)*adb;  break; }
case 2:{ I+=integrand(rho,theta)*(x)*(y)*adb;  break; }
case 3:{ I+=integrand(rho,theta)*(a-x)*(y)*adb;  break; }

// der shape
case 4:{ I+=integrand_dx(rho,theta)*(a-x)*(b-y)*adb;  break; }
case 5:{ I+=integrand_dx(rho,theta)*(x)*(b-y)*adb;  break; }
case 6:{ I+=integrand_dx(rho,theta)*(x)*(y)*adb;  break; }
case 7:{ I+=integrand_dx(rho,theta)*(a-x)*(y)*adb;  break; }

case 8:{ I+=integrand_dy(rho,theta)*(a-x)*(b-y)*adb;  break; }
case 9:{ I+=integrand_dy(rho,theta)*(x)*(b-y)*adb;  break; }
case 10:{ I+=integrand_dy(rho,theta)*(x)*(y)*adb;  break; }
case 11:{ I+=integrand_dy(rho,theta)*(a-x)*(y)*adb;  break; }

// der der
case 12:{ I+=integrand_dx(rho,theta)*(-1.0)*(b-y)*adb;  break; }
case 13:{ I+=integrand_dx(rho,theta)*(1.0)*(b-y)*adb;  break; }
case 14:{ I+=integrand_dx(rho,theta)*(1.0)*(y)*adb;  break; }
case 15:{ I+=integrand_dx(rho,theta)*(-1.0)*(y)*adb;  break; }

case 16:{ I+=integrand_dy(rho,theta)*(a-x)*(-1.0)*adb;  break; }
case 17:{ I+=integrand_dy(rho,theta)*(x)*(-1.0)*adb;  break; }
case 18:{ I+=integrand_dy(rho,theta)*(x)*(1.0)*adb;  break; }
case 19:{ I+=integrand_dy(rho,theta)*(a-x)*(1.0)*adb;  break; }

// der der sf



default:  
	{ MM_ERR(" bad shape specification "<<shape) }
}; // switch 

++n;
++itor;
}
I=I*a*b/D(n);
return I;


}


// integrate by parts and ignore remainder integral
// this is only one term and creates a complex value
// evaluated at arccos z=a/p for the lower BC.
// evaluate from rhomin=a, rhomax=sqrt(a*a+b*b)
// note that this is not an integration interval 
// quite but shoudl be as above
// it turns out the sin case is same as cos except w -> iw

// evaluate w^exp_top/(w*w+1)^exp_den note that ||w|| ==1
void  w_over_w2( D & re,  D & im, const D & _w_re, const D & _w_im,
const D & exp_top, const D & exp_den, const bool sine)
{
D w_re=0;
D w_im=0;
// sin is just multiplied by -i 
// also note that denominitor pole moves
if (sine)
{
w_re=-_w_im;
w_im=_w_re;
}
else { w_re=_w_re; w_im=_w_im; }

const D dex_re=w_re*w_re-w_im*w_im+(sine?(-1):1);
const D dex_im=2.0*w_re*w_im;
const D dex_r=std::pow(dex_re*dex_re+dex_im*dex_im,exp_den*.5);
const D theta_d=exp_den*std::atan2(dex_im,dex_re);
const D nex_r=std::pow(w_re*w_re+w_im*w_im,exp_top*.5);
const D theta_n=exp_top*std::atan2(w_im,w_re);

const D r=nex_r/dex_r;
const D theta=theta_n-theta_d;
re=r*std::cos(theta);
im=r*std::sin(theta);

}
void  cint_parts_cos( D & re, D & im, const D & alpha, const D & beta,
const D & e, const D & zmin, const IdxTy iter,
const bool add_min, const bool sine)
{
re=0;
im=0;
const D zmax=zmin;
 D w_re_max=zmax;
const D z2=zmax*zmax;
const bool rev=(z2>1);
 D w_im_max=std::sqrt(rev?(z2-1):(1-z2));
 D w_re_min=zmin;
 D w_im_min=std::sqrt(1-zmin*zmin);

D exp_top=e ;
D exp_den=alpha+2;
D fac=-.5/(exp_den);
for (IdxTy i=0; i<iter; ++i)
{
D rterm=0;
D iterm=0;
w_over_w2(rterm,iterm,w_re_max,w_im_max,exp_top,exp_den,sine);
re+=fac*rterm;
im+=fac*iterm;
if (false) if (zmax!=zmin)
{
w_over_w2(rterm,iterm,w_re_min,w_im_min,exp_top,exp_den,sine);
if (add_min) { re+=fac*rterm; im+=fac*iterm; } 
else { re-=fac*rterm; im-=fac*iterm; }
}

fac=fac*exp_top;
exp_den=exp_den-1;
exp_top=exp_top-2;
fac=fac*.5/(exp_den);

}


} 

////////////////////////////////////////////////

void  cint_parts_cos_alt( D & re, D & im, const D & alpha, const D & beta,
const D & e, const D & zmin, const IdxTy iter,
const bool add_min, const bool sine)
{
re=0;
im=0;
const D zmax=zmin;
 D w_re_max=zmax;
const D z2=zmax*zmax;
const bool rev=(z2>1);
 D w_im_max=std::sqrt(rev?(z2-1):(1-z2));
 D w_re_min=zmin;
 D w_im_min=std::sqrt(1-zmin*zmin);

D exp_top=e+1 ;
D exp_den=alpha+2;
//D fac=-.5/(exp_den);
D fac=1.0/(exp_top);
for (IdxTy i=0; i<iter; ++i)
{
D rterm=0;
D iterm=0;
w_over_w2(rterm,iterm,w_re_max,w_im_max,exp_top,exp_den,sine);
re+=fac*rterm;
im+=fac*iterm;
if (false) if (zmax!=zmin)
{
w_over_w2(rterm,iterm,w_re_min,w_im_min,exp_top,exp_den,sine);
if (add_min) { re+=fac*rterm; im+=fac*iterm; } 
else { re-=fac*rterm; im-=fac*iterm; }
}

fac=fac*exp_den;
exp_den=exp_den+1;
exp_top=exp_top+2;
fac=fac*2.0/(exp_den);

}


} 

////////////////////////////////////////////////

void  cint_parts_sin_alt( D & re, D & im, const D & alpha, const D & beta,
const D & e, const D & zmin, const IdxTy iter,
const bool add_min, const bool sine)
{
re=0;
im=0;
const D zmax=zmin;
 D w_re_max=zmax;
const D z2=zmax*zmax;
const bool rev=(z2>1);
 D w_im_max=std::sqrt(rev?(z2-1):(1-z2));
 D w_re_min=zmin;
 D w_im_min=std::sqrt(1-zmin*zmin);

D exp_top=e+1 ;
D exp_den=alpha+2;
//D fac=-.5/(exp_den);
D fac=1.0/(exp_top);
for (IdxTy i=0; i<iter; ++i)
{
D rterm=0;
D iterm=0;
w_over_w2(rterm,iterm,w_re_max,w_im_max,exp_top,exp_den,sine);
re+=fac*rterm;
im+=fac*iterm;
if (false) if (zmax!=zmin)
{
w_over_w2(rterm,iterm,w_re_min,w_im_min,exp_top,exp_den,sine);
if (add_min) { re+=fac*rterm; im+=fac*iterm; } 
else { re-=fac*rterm; im-=fac*iterm; }
}

fac=fac*exp_den;
exp_den=exp_den+1;
exp_top=exp_top+2;
fac=fac*2.0/(exp_den);

}

}



////////////////////////////////////////////////


void  cint_parts_sin( D & re, D & im, const D & alpha, const D & beta,
const D & e, const D & zmin, const IdxTy iter,
const bool add_min, const bool sine)
{
re=0;
im=0;
const D zmax=zmin;
 D w_re_max=zmax;
const D z2=zmax*zmax;
const bool rev=(z2>1);
 D w_im_max=std::sqrt(rev?(z2-1):(1-z2));
 D w_re_min=zmin;
 D w_im_min=std::sqrt(1-zmin*zmin);

D exp_top=e ;
D exp_den=alpha+2;
D fac=-.5/(exp_den);
for (IdxTy i=0; i<iter; ++i)
{
D rterm=0;
D iterm=0;
w_over_w2(rterm,iterm,w_re_max,w_im_max,exp_top,exp_den,sine);
re+=fac*rterm;
im+=fac*iterm;
if (false) if (zmax!=zmin)
{
w_over_w2(rterm,iterm,w_re_min,w_im_min,exp_top,exp_den,sine);
if (add_min) { re+=fac*rterm; im+=fac*iterm; } 
else { re-=fac*rterm; im-=fac*iterm; }
}

fac=fac*exp_top;
exp_den=exp_den-1;
exp_top=exp_top-2;
fac=fac*.5/(exp_den);

}


} 



// the integral has been re arranged and can be integrated 
// with the lower limit arccos(a/rho) for b>a and rho>a
// Uppse limit for arcsin(b/rho) similar with some modifications
D int_parts(const D & a, const D & b, const IdxTy nx, const IdxTy ny)
{
	D re=0, im=0,dre=0,dim=0;	
	const D alpha=m_alpha;
	const D beta=m_beta;
	const D ba=m_beta+alpha;
	const D mba=-m_beta+alpha;
	//const D dmin=(a<b)?a:b;
	//const D dmax=(a<b)?b:a;
	const D zmin=1; // a; 
	//const D zmax=a/sqrt(a*a+b*b); 
	const IdxTy iter=nx;
	const bool add_min=true;
	const IdxTy n_exp=6;
	bool sine=false;
	const D rmax=sqrt(a*a+b*b);
	const D theta_diag=atan2(b,a);
	const D pref=1.0/beta/(alpha+2.0);
	const D zmaxa=a/rmax; // a; 
//	const D zmaxb=b/rmax; // a; 
	// for small r, the angle limits are 0 to pi/2.
	re+=pref*std::pow(b,alpha+2)*(cos(beta*M_PI*.5));
	re-=pref*std::pow(b,alpha+2)*(cos(beta*theta_diag));
	MM_MSG(" complex integral "<<re<<" "<<im)
	re+=pref*std::pow(a,alpha+2)*(cos(beta*theta_diag));
	re-=pref*std::pow(a,alpha+2)*(cos(0));
	MM_MSG(" complex integral "<<re<<" "<<im)


	D exps[n_exp],ac[n_exp],as[n_exp];
//		exps[0]=ba+3.0; exps[1]=mba+3.0; exps[2]=ba+1.0; exps[3]=mba+1.0;
		exps[0]=ba+4.0; exps[1]=ba+3.0; exps[2]=ba+2.0;
		exps[3]=mba+4.0; exps[4]=mba+3.0; exps[5]=mba+2.0;
		ac[0]=1; ac[1]=2; ac[2]=-1; ac[3]=1; ac[4]=2; ac[5]=-1;
		as[0]=1; as[1]=-2; as[2]=1; as[3]=1; as[4]=-2; as[5]=1;

	const D cf=.5/beta*std::pow(2.0*a,alpha+2);

	const D sf=.5/beta*std::pow(2.0*b,alpha+2);
	const D sftheta=M_PI*.5*(alpha+2);
	const D sfr=sf*cos(sftheta);
	const D sfi=sf*sin(sftheta);

	// the exponents soone need to be rationals
	for (IdxTy i=0; i<n_exp; ++i)
	{
	cint_parts_cos_alt(dre,dim,m_alpha,m_beta,exps[i],zmin,iter,add_min,sine);
	re-=cf*ac[i]*dre; im-=cf*ac[i]*dim;
	MM_MSG(" dre "<<dre<<" "<<dim)
	MM_MSG(" complex integral "<<re<<" "<<im)
	cint_parts_cos(dre,dim,m_alpha,m_beta,exps[i],zmaxa,iter,add_min,sine);
	re+=cf*ac[i]*dre; im+=cf*ac[i]*dim;
	MM_MSG(" dre "<<dre<<" "<<dim)
	MM_MSG(" complex integral "<<re<<" "<<im)
	
}
	sine=true;
	for (IdxTy i=0; i<n_exp; ++i)
	{
	cint_parts_sin_alt(dre,dim,m_alpha,m_beta,exps[i],zmin,iter,add_min,sine);
	re+=sfr*as[i]*dre; im+=sfi*as[i]*dim;
	MM_MSG(" dre "<<dre<<" "<<dim)
	//cint_parts_cos(dre,dim,m_alpha,m_beta,exps[i],zmaxb,iter,add_min,sine);
	//re-=sfr*as[i]*dre; im-=sfi*as[i]*dim;
	//MM_MSG(" dre "<<dre<<" "<<dim)

	}


/*
	cint_parts_cos(dre,dim,m_alpha,m_beta,ba+3.0,rhomin,rhomax,iter,add_min,sine);
	re+=dre; im+=dim;
	cint_parts_cos(dre,dim,m_alpha,m_beta,mba+3.0,rhomin,rhomax,iter,add_min,sime);
	re+=dre; im+=dim;
	cint_parts_cos(dre,dim,m_alpha,m_beta,ba+1.0,rhomin,rhomax,iter,add_min,sinr);
	re+=dre; im+=dim;
	cint_parts_cos(dre,dim,m_alpha,m_beta,mba+1.0,rhomin,rhomax,iter,add_min);
	re+=dre; im+=dim;
*/
	MM_MSG(" complex integral "<<re<<" "<<im)
	return re;
}

D m_alpha,m_beta;
D m_xc,m_yc;
D m_rot;
} ; // as2

class angle_scale
{
typedef angle_scale Myt;
typedef mjm_integrals Super;
typedef Super Tr;
typedef Tr::D D;
typedef Tr::IdxTy IdxTy;
public:

angle_scale() : m_alpha(0),m_a(1), m_z(1.0/(1.0+m_alpha)) {}
angle_scale(const D & alpha, const D & a) : m_alpha(alpha),m_a(a), m_z(1.0/(1.0+m_alpha)) {}

D integrand(const D & rho, const D & alpha, const D & a)
{ return std::pow(rho,alpha)*sqrt(rho*rho-a*a); }

D integrand(const D & rho)
{ return integrand(rho, m_alpha, m_a); } 

D integral(const D & rho, const D & alpha, const D & a)
//{ const D z=1.0/(1.0+alpha); const D v=std::pow(rho,z)-a*a; return .25*std::pow(v*v*v*v+a*a,2.0+2.0*alpha); }
{ const D z=1.0/(1.0+alpha); const D v=std::pow(rho,z)-a*a; return .25*std::pow(v*v*v*v+a*a,8.0+8.0*alpha); }

D integral(const D & rho)
{ return integral(rho, m_alpha, m_a); } 


// this should integrate r^(gamma)sqrt(r*r-b*b) from rhomin to rhomax
// for r>b, integrate once by parts so we hage a derivative 
D integrate_rect_bc(const D & gamma,const D & a,const D & rhomax,const D & rhomin)
{
const D third=1.0/3.0;
const D r2max=rhomax*rhomax;
const D r2min=rhomin*rhomin;
const D a1=gamma-1.0;
const D a2=a*a;
const D iparts=third*std::pow(rhomax,a1)*std::pow(r2max-a2,1.5)
	-third*std::pow(rhomin,a1)*std::pow(r2min-a2,1.5);

if (rhomax==rhomin) return 0;
const D tol=1e-7;
D I=0;
D drmin=(rhomax-rhomin)*.0001;
D drmax=(rhomax-rhomin)*.001;
D r=rhomin;
D last=0;
IdxTy steps=0;
while (r<=rhomax)
{
const D r2=r*r;
const D di=third*a1*std::pow(r,gamma-2.0)*std::pow(r2-a2,1.5);
const D d2i=a1*std::pow(r,gamma-3.0)*sqrt((r2-a2))*(third*(r2-a2)*(gamma-2.0)+r*r);
D dr=0;
if (d2i==0) dr=drmax;
else {
D dx=tol*di/d2i;
if (dx<0) dx=-dx;
if (dx>drmax) dr=drmax;
else if (dx<drmin) dr=drmin;
else dr=dx;
} // dr
I-=di*(last+dr)*.5;
last=dr;
r+=dr;
//MM_MSG( steps<<" "<<r<<" "<<dr<<" "<<I<<" "<<di<<" "<<d2i)
++steps;
}
//MM_MSG(" steps = "<<steps<<" i "<<I<<" parts: "<<iparts<<" frac "<<((iparts)/(I+iparts)))
I=I+iparts;

return I;
}
// for integral exponent gamma the thing could terminate or be differnet, need
// to check that. 1
D integrate_rect_bc_po(const D & gamma,const D & a,const D & rhomax,const D & rhomin,const IdxTy iter)
{
const D third=1.0/3.0;
const D r2max=rhomax*rhomax;
const D r2min=rhomin*rhomin;
const D r2maxi=1.0/r2max;
const D r2mini=1.0/r2min;
// now usnused 
//const D a1=gamma-1.0;
const D a2=a*a;
// this ie no really needed, just sqrt at min and max as well as gamma power 
// after that all recursion of integral powers. 
//const D iparts=third*std::pow(rhomax,a1)*std::pow(r2max-a2,1.5)
//	-third*std::pow(rhomin,a1)*std::pow(r2min-a2,1.5)

// -a1*third*.2*std::pow(rhomax,gamma-3)*std::pow(r2max-a2,2.5)
// +a1*third*.2*std::pow(rhomin,gamma-3)*std::pow(r2min-a2,2.5)

//;

//D last=-a1*third*.2*(gamma-3);
//D last=-a1*third *.2; // (gamma-3);
D last=third ; // (gamma-3);

//IdxTy j=7; // 7
IdxTy j=3; // 7
D maxn=std::pow(r2max-a2,j*.5); // 2.5
D minn=std::pow(r2min-a2,j*.5);
D maxd=std::pow(rhomax,gamma-j+2); // 3
D mind=std::pow(rhomin,gamma-j+2);

D maxf=maxd*maxn;
D minf=mind*minn;
// this only works for j=3 otherwise there is a sum etc. 
const D iparts= last*maxf-last*minf; 

D I=iparts;

j=j+2;

// these are the factors for each iteration

D maxrecur=(r2max-a2)*r2maxi;
D minrecur=(r2min-a2)*r2mini;
// the reverse sum does show some dufferences but nothing relevant yet 
//D reverse_sum[iter];
// the maxf version appears to work or match trapeziod ad absurdium
for (IdxTy i=0; i<iter; ++i)
{
const D fac=1.0/j;
const D next=-last*fac* (gamma-j+4);
// recursion here will be very very fast ... 
// this seems to work, now try recursion
//const D inext= next*std::pow(rhomax,gamma-j+2)*std::pow(r2max-a2,.5*j)
//-next*std::pow(rhomin,gamma-j+2)*std::pow(r2min-a2,.5*j);
//maxn=maxn*(r2max-a2);
//minn=minn*(r2min-a2);
//maxd=maxd*r2maxi;
//mind=mind*r2mini;

maxf=maxf*maxrecur;
minf=minf*minrecur;
// this works up to a point 
//const D inext= next*maxd*maxn -next*mind*minn;
const D inext= next*(maxf -minf);
//reverse_sum[i]=inext;
I+=inext;
//last=next*(gamma-j+2);
last=next;
j=j+2;
}
//D J=0;
//for (IdxTy i=iter; i!=0; --i) J+=reverse_sum[i-1];
//J+=iparts;
//MM_MSG(" summation "<<I<<" "<<J<<" "<<(I-J))

return I;


}







// 2 integration by parts, should though change vars to x=b/r... 
D integrate_rect_bc_2(const D & gamma,const D & a,const D & rhomax,const D & rhomin)
{
const D third=1.0/3.0;
const D r2max=rhomax*rhomax;
const D r2min=rhomin*rhomin;
const D a1=gamma-1.0;
const D a2=a*a;
const D iparts=third*std::pow(rhomax,a1)*std::pow(r2max-a2,1.5)
	-third*std::pow(rhomin,a1)*std::pow(r2min-a2,1.5)
-a1*third*.2*std::pow(rhomax,gamma-3)*std::pow(r2max-a2,2.5)
+a1*third*.2*std::pow(rhomin,gamma-3)*std::pow(r2min-a2,2.5);


if (rhomax==rhomin) return 0;
const D tol=1e-4;
D I=0;
D drmin=(rhomax-rhomin)*.00001;
D drmax=(rhomax-rhomin)*.01;
D r=rhomin;
D last=0;
IdxTy steps=0;
while (r<=rhomax)
{
const D r2=r*r;
const D di=third*a1*(gamma-3.0)*.2*std::pow(r,gamma-4.0)*std::pow(r2-a2,2.5);
const D d2i=third*a1*(gamma-3.0)*.2*((gamma-4)*std::pow(r,gamma-5.0)*std::pow(r2-a2,2.5)
		+5*std::pow(r,gamma-3.0)*std::pow(r2-a2,1.5) );

//const D d2i=a1*third*.2*(gamma-4.0)*std::pow(r,gamma-5.0)*std::pow((r2-a2),1.5)*(1.5*(r2-a2)*(gamma-2.0)+r*r);
D dr=0;
if (d2i==0) dr=drmax;
else {
D dx=tol*di/d2i;
//MM_MSG(" dx "<<dx<<" di "<<di<<" d2i "<<d2i<<" "<<drmax<<" "<<drmin)
if (dx<0) dx=-dx;
if (dx>drmax) dr=drmax;
else if (dx<drmin) dr=drmin;
else dr=dx;
} // dr
I+=di*(last+dr)*.5;
last=dr;
r+=dr;
//MM_MSG( steps<<" "<<r<<" "<<dr<<" "<<I<<" "<<di<<" "<<d2i)
++steps;
}
//MM_MSG(" steps = "<<steps<<" pf= "<<(iparts/I)<<" "<<iparts<<" i "<<I)
//MM_MSG(" v2  steps = "<<steps<<" i "<<I<<" parts: "<<iparts<<" frac "<<((iparts)/(I+iparts)))
I=I+iparts;

return I;
}





private:
D m_alpha,m_a, m_z;

}; // angle_scale
///////////////////////////////////////////////////////

// this is really a special case of compositing, 
template<class Td > 
	void reverse_polynomial(Td & dest, const Td & src)
{
const IdxTy sz=src.size();
dest= Td(sz);
for (IdxTy j=0; j<sz; ++j)
{ dest[j]=src[sz-1-j];
} 
} // reverse_polynomial



//////////////////////////////////////////////////////

// use binomial formula and rearrnage 
// this is really a special case of compositing, 
template<class Td,class Tv > 
	void shift_polynomial(Td & dest, const Td & src, const Tv & a)
{
const IdxTy sz=src.size();
dest= Td(sz);
Tv aij[sz];
aij[0]=1;
for (IdxTy j=1; j<sz; ++j) aij[j]=a*aij[j-1];
for (IdxTy j=0; j<sz; ++j)
{
auto & dj=dest[j];
IdxTy fac=1; // TODO FIXME need to do the math here... 
for (IdxTy i=j; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ 
	dj+=src[i]*aij[i-j]*Tv(fac);
	fac=fac*(i+1)/(i-j+1); //   this may be right now 

} 
} 

} // shift_polynomial

// return g(f(x))
template<class Td > 
	void composite_polynomial(Td & dest, const Td & g, const Td & f)
{
//typedef typename Td::value_type Tv;
const IdxTy szf=f.size();
const IdxTy szg=g.size();
if (szf==0) return;
if (szg==0) return;
Td fn=f;
const IdxTy sz=1+(szf-1)*(szg-1);
dest= Td(sz);
if (sz<1) return;
dest[0]=g[0];
//Tv aij[sz];
//aij[0]=1;
//for (IdxTy j=1; j<sz; ++j) aij[j]=a*aij[j-1];
for (IdxTy j=1; j<szg; ++j)
{
// this is really an accumulating polynomical add 
const IdxTy szfn=fn.size();
for (IdxTy i=0; i<szfn; ++i)
{
dest[i]+=g[j]*fn[i];
}
Td next;
multiply_polynomials(next,f,fn);
fn=next;
}

} // composite_polynomial

// integrate x^(n)exp(-b*x^2) by parts and return the polynomale
// factor for exp(-bx^2) and the coeffieinct of the erf term 
//template<class Td >  typename Td::value_type x_n_exp_integral
template<class Td,class Tvalt  >  Tvalt x_n_exp_integral
	//(Td & dest, const IdxTy  & n, const typename Td::value_type & b, const bool
	(Td & dest, const IdxTy  & n, const Tvalt  & b, const bool
init_dest=true, const Tvalt & scale=1)
{
//typedef typename Td::value_type Tval;
typedef Tvalt Tval;
Tval coef_erf=0;
// returning zero size shoule be ok 
if (init_dest) dest= Td(n+0*1);
else
{
if (dest.size()<(n+0*1)) { MM_ERR(" wrong size for exp series "<<MMPR2(n,dest.size())) } 
}
Tval fac=scale; // 1;
IdxTy ni=n;
const Tval f1=.5/b;
const Tval f2=sqrt(M_PI/b)*.5;
while ( ni>1)
{
dest[ni-1]+=fac*(-f1); // (-.5/b);
fac=fac*(ni-1)*f1; // /2.0/b;
ni-=2;
}
if (ni==0) coef_erf=fac*f2; // sqrt(M_PI)/2.0/sqrt(b);
//if (ni==1) dest[0]=-f1; // .5*fac/(b);
if (ni==1) dest[0]+=-fac*f1; // .5*fac/(b);

return coef_erf;
}




// integrate x^(n)erf(sqrt(b)*x) by parts and return the polynomale
// factor for exp(-bx^2) and the coeffieinct of the erf term 
// not that one integraion by parts converts this into an exp integral. 
// the result is the polynomial time exp(-bx^2) and erf()(x^(n+1) + c ) 
// so just return an erf poly too  even if wasteful 
template<class Td, class Tvalt >  void x_n_erf_integral
	//(Td & derf, Td & dexp, const IdxTy  & n, const typename Td::value_type & b
	(Td & derf, Td & dexp, const IdxTy  & n, const Tvalt & b
, const bool init_dest=true, const Tvalt & scale=1 )
{
//typedef typename Td::value_type Tval;
//typedef typename Tvalt  Tval;
typedef  Tvalt  Tval;
const  Tval expf=-1.0/(n+1)*2*sqrt(b/M_PI);
if (init_dest) { derf = Td(n+2); dexp=Td(n+1); } 
else 
{ // TODO use asserts here 
if (derf.size()<(n+2)) { MM_ERR(" wrong size for erf series "<<MMPR3(n,derf.size(),dexp.size())) } 
if (dexp.size()<(n+1)) { MM_ERR(" wrong size for exp series "<<MMPR3(n,derf.size(),dexp.size())) } 

}
derf[n+1]+=1.0/(n+1)*scale;
derf[0]+= x_n_exp_integral(dexp, n+1,  b,init_dest,expf*scale);
// this fails now is dexp had not been inited so use scale . 
//const IdxTy sze=dexp.size();
//for (IdxTy i=0; i<sze; ++i) dexp[i]*=expf;

}


//////////////////////////////////////////////////////////////////////////
// versions of above using rational erf or I call rerf lacking the sqrt prefactor
// remove 2sqrt(b)/sqrt(pi)

template<class Td,class Tvalt  >  Tvalt x_n_rexp_integral
	(Td & dest, const IdxTy  & n, const Tvalt  & b, const bool
init_dest=true, const Tvalt & scale=1)
{
typedef Tvalt Tval;
Tval coef_rerf=0;
if (init_dest) dest= Td(n);
else { if (dest.size()<(n)) { MM_ERR(" wrong size for exp series "<<MMPR2(n,dest.size())) } }
Tval fac=scale; // 1;
IdxTy ni=n;
const Tval f1=.5/b;
// (-.5/b); // /2.0/b;
while ( ni>1) { dest[ni-1]+=fac*(-f1); fac=fac*(ni-1)*f1; ni-=2; } 
if (ni==0) coef_rerf=fac; // *f2; // sqrt(M_PI)/2.0/sqrt(b);
if (ni==1) dest[0]+=-fac*f1; // .5*fac/(b);
return coef_rerf;
}

template<class Td, class Tvalt >  void x_n_rerf_integral
	(Td & derf, Td & dexp, const IdxTy  & n, const Tvalt & b
, const bool init_dest=true, const Tvalt & scale=1 )
{
typedef  Tvalt  Tval;
const  Tval expf=-1.0/(n+1); // *2*sqrt(b/M_PI);
if (init_dest) { derf = Td(n+2); dexp=Td(n+1); } 
else 
{ // TODO use asserts here 
if (derf.size()<(n+2)) { MM_ERR(" wrong size for erf series "<<MMPR3(n,derf.size(),dexp.size())) } 
if (dexp.size()<(n+1)) { MM_ERR(" wrong size for exp series "<<MMPR3(n,derf.size(),dexp.size())) } 

}
derf[n+1]+=1.0/(n+1)*scale;
derf[0]+= x_n_rexp_integral(dexp, n+1,  b,init_dest,expf*scale);

}













// dense vector with v[i] being coefficient on x^i
template<class Td>  void multiply_polynomials(Td & dest, const Td & v1, const Td & v2)
{
const IdxTy sz1=v1.size();
const IdxTy sz2=v2.size();
const IdxTy sz=sz1+sz2;
if (sz==0) { dest=Td(); return ; } 
dest= Td(sz-1);
for (IdxTy i=0; i<sz1; ++i)
for (IdxTy j=0; j<sz2; ++j)
// I don't have a += operator for rationals lol 
{ dest[i+j]=dest[i+j]+v1[i]*v2[j]; } 
} // multiply_polynomials

template<class Td>  void power_polynomial(Td & dest, const Td & v1, const IdxTy n )
{
Td q,r;
q.push_back(1);
for(IdxTy i=0; i<n; ++i)
{
multiply_polynomials(r,q,v1);
q=r;
}
dest=q;
trim_polynomial(dest);
}


template<class Td>  void multiply_polynomials_nz(Td & dest, const Td & v1, const Td & v2)
{
const IdxTy sz1=v1.size();
const IdxTy sz2=v2.size();
const IdxTy sz=sz1+sz2;
if (sz==0) { dest=Td(); return ; } 
dest= Td(sz-1);
for (IdxTy i=0; i<sz1; ++i)
{
const auto & v1i=v1[i];
if (v1i==0) continue; 
for (IdxTy j=0; j<sz2; ++j)
// I don't have a += operator for rationals lol 
{ if (v2[j]!=0) dest[i+j]=dest[i+j]+v1[i]*v2[j]; } 
}
} // multiply_polynomials




template<class Td>  void multiply_polynomial(Td & dest, const Td & v1, const typename Td::value_type & x)
{
const IdxTy sz=v1.size();
dest= Td(sz);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ dest[i]=v1[i]*x; } 
} // multiply_polynomials


template<class Td>  void multiply_polynomial_nz(Td & dest, const Td & v1, const typename Td::value_type & x)
{
const IdxTy sz=v1.size();
dest= Td(sz);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ if (v1[i]!=0) dest[i]=v1[i]*x; } 
} // multiply_polynomials




template<class Td>  void add_polynomials(Td & dest, const Td & v1, const Td & v2)
{
const IdxTy sz1=v1.size();
const IdxTy sz2=v2.size();
const IdxTy sz=(sz1>sz2)?sz1:sz2;
dest= Td(sz);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ 
if(i<sz1) dest[i]+=v1[i];//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
if(i<sz2) dest[i]+=v2[i];//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
} 


} // add_polynomials
template<class Td>  void subtract_polynomials(Td & dest, const Td & v1, const Td & v2)
{
const IdxTy sz1=v1.size();
const IdxTy sz2=v2.size();
const IdxTy sz=(sz1>sz2)?sz1:sz2;
dest= Td(sz);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ 
if(i<sz1) dest[i]+=v1[i];//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
if(i<sz2) dest[i]-=v2[i];//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
} 


} // add_polynomials



template<class Td>  void accumulate_polynomial(Td & dest, const Td & v1)
{
typedef typename Td::value_type Tvv;
const IdxTy sz1=v1.size();
//const IdxTy sz2=dest.size();
//const bool too_small=(sz2<sz1);
const IdxTy sz=sz1; // (too_small)?sz1:sz2;
while (dest.size()<sz1) { dest.push_back(Tvv()); } 
//dest= Td(sz);
//if (too_small) dest.resize(sz1);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ dest[i]+=v1[i]; } 
} // accumulate_polynomial




template<class Td,class Tvv >  void add_polynomials(Td & dest, const Td & v1, const Td & v2, const Tvv & scale1, const Tvv & scale2=Tvv(1))
{
const IdxTy sz1=v1.size();
const IdxTy sz2=v2.size();
const IdxTy sz=(sz1>sz2)?sz1:sz2;
dest= Td(sz);
for (IdxTy i=0; i<sz; ++i)
// I don't have a += operator for rationals lol 
{ 
if(i<sz1) dest[i]+=v1[i]*scale1;//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
if(i<sz2) dest[i]+=v2[i]*scale2;//	dest[i+j]=dest[i+j]+v1[i]*v2[j]; 
} 


} // add_polynomials



template<class Tx, class Ty, class Tn > bool  add_term( Tx & terms, Ty & ss, const Tn & t)
{
	if (t==0) return false;
 	if (terms!=0) { if ( t>0) ss<<"+"; } 
	if (t==(-1)) { ss<<"-"; }
	else if (t==1) {}
	else ss<<t; 
	++terms; 
	return true;  
//	return false; 
}
template<class Td>  std::string print_polynomial( const Td & x2, const IdxTy flags=0, const StrTy & var="x")
{
std::stringstream ss;
const int x2s=x2.size();
if (x2s==0) return ss.str();
int terms=0;
for (int j=x2s-1; j>0; --j)
//	  if ( x2[j]!=0) 
		//{ if (terms!=0) ss<<"+"; ss<<x2[j]; ++terms; { ss<<"x^"<<j;} }
		{ if ( add_term(terms,ss,x2[j]))   
			{ if ((j==0) ||((x2[j]*x2[j])==1)) {ss<<"1"; } if (j==1) {ss<<var; } else {ss<<var<<"^"<<j;}} }

//if (x2s>0) {if (x2[0]!=0) {  if (terms!=0) ss<<"+"; ss<<x2[0];} }
if (x2s>0) {add_term(terms,ss,x2[0]); }
if  (((x2[0]*x2[0])==1)) {ss<<"1"; } 
if (terms==0) if  (((x2[0]*x2[0])==0)) {ss<<"0"; } 
return ss.str();
}

// this should work with rational coeffs too
template<class Td>  void integrate_polynomial(Td & dest, const Td & v1)
{
const IdxTy sz1=v1.size();
dest= Td(sz1+1);
dest[0]=0;
for (IdxTy i=0; i<sz1; ++i)
{
// I guess I should define a rational divide by int operator
dest[i+1]=v1[i]/(1+i);

} // sz1

}

template<class Td,class Tvalt>  void integrate_polynomial_e(Td & derf, Td& dexp
		, const Td & v1,const Tvalt & b, const IdxTy factor, const bool init_dest=true )
{
const IdxTy sz1=v1.size();
//if (init_dest) { dexp= Td(sz1+1); derf=Td(sz1+2); } 
if (init_dest) { dexp= Td(sz1); derf=Td(sz1+1); } 
//dest[0]=0;
for (IdxTy i=0; i<sz1; ++i)
{
// I guess I should define a rational divide by int operator
//dest[i+1]=v1[i]/(1+i);
// Tvalt x_n_exp_integral (dest, n, b,  init_dest,scale);
if (factor==1) derf[0]+= x_n_exp_integral (dexp, i, b,  false,v1[i]);
// x_n_erf_integral ( derf, dexp,  n,  b , init_dest, scale );
else  x_n_erf_integral ( derf, dexp,  i,  b , false, v1[i] );


} // sz1

}


template<class Td>  void differentiate_polynomial(Td & dest, const Td & v1)
{
const IdxTy sz1=v1.size();
if (sz1<1) return; 
dest= Td(sz1-1);
for (IdxTy i=0; i<(sz1-1); ++i)
{
// I guess I should define a rational divide by int operator
dest[i]=v1[i+1]*(1+i);
} // sz1
}

template<class Td>  void differentiate_polynomial(Td & dest, const Td & v1,const IdxTy n )
{
const IdxTy sz1=v1.size();
if (sz1<n) return; 
dest= Td(sz1-n);
for (IdxTy i=0; i<(sz1-n); ++i)
{
// I guess I should define a rational divide by int operator
//dest[i]=v1[i+n]*(1+i);
typename Td::value_type c=1; // +i;
for(IdxTy j=i+1; j<=(i+n); ++j) c=c*j;
dest[i]=v1[i+n]*c;
} // sz1
}




// remove higher order zero 
template<class Td>  void trim_polynomial(Td & dest)
{
Td d;
polynomial_trim(d,dest);
dest=d;
//if (d!=dest) 
//		MM_MSG(MMPR2(print_polynomial(dest),print_polynomial(d)))
}

template<class Td>  void polynomial_trim(Td & dest, const Td & v1)
{
dest.clear();
//typedef typename Td::value_type Tp;
const IdxTy sz=v1.size();
if (sz==0) return;
int i=sz-1;
for (; i>=0; --i) if (v1[i]!=0)  break;
if (i<0) { dest=v1; return; } 
for(int j=0; j<=i; ++j) dest.push_back(v1[j]);
}
// remove factors of "x" common to all polynomials, nothing more lol 
template<class Td>  void polynomial_list_reduce(Td & dest, const Td & v1)
{
dest.clear();
typedef typename Td::value_type Tp;
const IdxTy sz=v1.size();
if (sz==0) return;
IdxTy highest=v1[0].size();
for (IdxTy i=0; i<sz; ++i)
{
const Tp & p=v1[i];
const IdxTy szp=p.size();
for(IdxTy j=0; j<szp; ++j) 
{
if (p[j]!=0)  { if (j<highest) highest=j; break; }
} // j 
} //i 
for (IdxTy i=0; i<sz; ++i)
{
dest.push_back(Tp());
Tp & d=dest[i];
const Tp & p=v1[i];
const IdxTy szp=p.size();
for(IdxTy j=highest; j<szp; ++j) 
{
d.push_back(p[j]);
//if (p[j]!=0)  { if (j<highest) highest=j; break; }
} // j 
} //i 



}

template<class Td>  void gauss_taylor_polynomial(Td & dest, const IdxTy n )
{
typedef typename Td::value_type  Tval;
Tval fac=1;
for (IdxTy i=0; i<n; ++i)
{
//const bool odd=((i&1)!=0);
dest.push_back(1.0/fac);
dest.push_back(0);
fac=-fac*Tval(i+1);
}

}
// this apparently gets term n from a polynomial times a gaussian expansion 
template<class Td>  typename Td::value_type  gauss_term_polynomial( const Td & p,const IdxTy n )
{
typedef typename Td::value_type  Tval;
Tval fac=1;
Tval sum=0;
IdxTy tmax=(n+1)>>1;
const IdxTy sz=p.size();
for (IdxTy i=0; i<tmax; ++i)
{
int j=n-2*i;
if (j<0) return sum;
if (j<int(sz)) sum+=p[j]/fac;
//const bool odd=((i&1)!=0);
//dest.push_back(1.0/fac);
//dest.push_back(0);
fac=-fac*Tval(i+1);
}
return sum;
}  // gauss_term_polynomial


// the polynomial prefactor obtained by differentiating p(x)exp(-x^2)
// for sequetial derivatives of exp(-x^2)
template<class Td>  void next_gauss_der_polynomial(Td & dest, const Td & v1, const typename Td::value_type & scale=1.0)
{
typedef typename Td::value_type Tvv;
const IdxTy sz1=v1.size();
if (sz1<1) return; 
const IdxTy dsize=sz1+1;
dest= Td(dsize);
if (sz1>1) dest[0]=v1[1];
for (IdxTy i=1; i<(sz1-1); ++i)
{
//dest[i]=v1[i+1]*(1+i)-2.0*v1[i-1];
dest[i]=v1[i+1]*(1+i)-v1[i-1]-v1[i-1];

} // sz1
//if (sz1>1) dest[sz1-1]=-2.0*v1[sz1-2];
if (sz1>1) dest[sz1-1]=-v1[sz1-2]-v1[sz1-2];
//dest[sz1]= -2.0*v1[sz1-1];
dest[sz1]= -v1[sz1-1]-v1[sz1-1];
if (scale!=Tvv(1))
{
for (IdxTy i=0; i<(dsize); ++i) dest[i]=scale*dest[i];


}
}
template<class Td>  void next_gauss_der_beta_polynomial(Td & dest, const Td & v1, const typename Td::value_type & beta, const typename Td::value_type & scale=1.0)
{
typedef typename Td::value_type Tvv;
const IdxTy sz1=v1.size();
if (sz1<1) return; 
const IdxTy dsize=sz1+1;
dest= Td(dsize);
if (sz1>1) dest[0]=v1[1];
for (IdxTy i=1; i<(sz1-1); ++i)
{
//dest[i]=v1[i+1]*(1+i)-2.0*v1[i-1];
const Tvv bder=(-v1[i-1]-v1[i-1])*beta;
//dest[i]=v1[i+1]*(1+i)-v1[i-1]-v1[i-1];
dest[i]=v1[i+1]*(1+i)+bder; // v1[i-1]-v1[i-1];

} // sz1
//if (sz1>1) dest[sz1-1]=-2.0*v1[sz1-2];
if (sz1>1) dest[sz1-1]=(-v1[sz1-2]-v1[sz1-2])*beta;
//dest[sz1]= -2.0*v1[sz1-1];
dest[sz1]= (-v1[sz1-1]-v1[sz1-1])*beta;
if (scale!=Tvv(1))
{
for (IdxTy i=0; i<(dsize); ++i) dest[i]=scale*dest[i];


}
}






template <class Td,class Tv> Tv evaluate_polynomial(const Td & p,const Tv & x)
{
Tv v;
v=0;
const IdxTy sz=p.size();
if (sz==0) return v; 
Tv xx;
xx=1;
// this overflows oftern
if (false){
for (IdxTy i=0; i<sz; ++i) { v=v+p[i]*xx; xx=xx*x; } }
//MM_ERR(" fudding poly x "<<x<<" xx "<<xx<<" p "<<p[i]<<" v "<<v)
// I think this is the Horner method... 
for (IdxTy i=(sz-1); i>0; --i) { v=(v+p[i])*x; } 
v=v+p[0];
return v;
} 
template <class Td,class Tv> Tv evaluate_polynomial_nz(const Td & p,const Tv & x)
{
Tv v;
v=0;
const IdxTy sz=p.size();
if (sz==0) return v; 
Tv xx;
xx=1;
// this overflows oftern
if (false){
for (IdxTy i=0; i<sz; ++i) { v=v+p[i]*xx; xx=xx*x; } }
//MM_ERR(" fudding poly x "<<x<<" xx "<<xx<<" p "<<p[i]<<" v "<<v)
for (IdxTy i=(sz-1); i>0; --i) { if (p[i]==0) v=v*x; else v=(v+p[i])*x; } 
v=v+p[0];
return v;
} 


/////////////////////////////////////////////////////////////////////////////


// difference of two erf's and two normals
// erf(xf)-erf(xi) and exp(xf2)-exp(xi^2)
int doen_old( D & derf, D & dnorm, const D & xf, const D & xi, const IdxTy n, const bool debug=false)
{
//typedef unsigned int IdxTy;
typedef double D ;
//typedef std::string StrTy;
//typedef std::stringstream Ss;
//typedef mjm_integrals Myt;
typedef std::vector<D>  P;
P x1,x2;
const D alpha=.1;
//Myt mi;
const D delx=(xf-xi)*alpha;
const D mean=.5*(xf+xi)*alpha;
Myt &  mi= *this;
x1.push_back(1);
D dxn=delx*.5;
D fac=1.0;
D sum=0;
const D tol=1e-14;
const bool dump=debug; // (n==0);
const bool auto_stop=(n==0);
const int neff=(n==0)?1000:n;
for (int i=0; i<neff; ++i)
{
D coef=mi.evaluate_polynomial(x1,mean);
D term=coef*dxn/fac;
if ((i&1)==0) sum+=term;
fac=fac*(i+2);
dxn=dxn*delx*.5;
if (dump) {
	MM_MSG(mi.print_polynomial(x1));
	MM_MSG(" terms "<<MMPR4(i,term,coef,dxn/fac)<<MMPR(sum))
	}
if (i>3) if (auto_stop)
{
const D t=(term<0)?(-term):term;
const D v=(sum<0)?(-sum):sum;
if (t<(v*tol)) break;
}
//template<class Td>  void next_gauss_der_polynomial(Td & dest, const Td & v1)
mi.next_gauss_der_polynomial(x2,x1);
x1=x2;
}
static const D ps=sqrt(M_PI);
derf=2.0*sum*exp(-mean*mean)*2.0/ps/alpha;
//derf=sum*erf(mean);
return 0;
}
// erf(xf)-erf(xi) and exp(xf2)-exp(xi^2)
////////////////////////////////////////////////////////////////////////////
// this seems to return derf ok, need to do the norm now 
int doen_split( D & derf, D & dnorm, const D & xf, const D & xi, const IdxTy n, const bool debug=false)
{
//D xxx=0;
D xd=xf-xi;
D res=0;
D resexp=0;
D res2=0;
D res2exp=0;
if (xd>1)
{
	D reslut=0;
	D reslutexp=0;
	int imax=xf;
	while ((imax+.5)>xf ) --imax;
	int imin=xi;
	while ((imin-.5)<xi ) ++imin;
	if (imax>=imin){ // break; // get at least one 
	const D xtabu=D(imax)+.5;
	const D xtabl=D(imin)-.5;
	for (int i=imin; i<=imax; ++i)
	{
		 //const D lv=doen_lut(i);
		D lv=0, lvexp=0;
		 //const IdxTy lutrc=
			doen_lut(lv,lvexp,i);
	//	MM_MSG(MMPR3(i,lv,reslut))
		reslut+=lv;
		reslutexp+=lvexp;
	}	
	doen(res,resexp,xf,xtabu,n,false);
	int rc=doen(res2,res2exp,xtabl,xi,n,false);
	derf= reslut+res+res2;
	dnorm= reslutexp+resexp+res2exp;
	return rc; 
}}
	int rc=doen(res2,res2exp,xf,xi,n,debug);

	derf=res+res2;
	dnorm=resexp+res2exp;
	return rc; 
}
// this should actually return a fail for oor
//D doen_lut(const int n)
int  doen_lut(D & erfi, D & expi,const int n)
{
static bool init=false;
static const int nmin=-100;
static const int nmax=100;
static const int ntot=nmax-nmin+1;

//static D lut[ntot];
static D lut_erf[ntot];
static D lut_exp[ntot];
//D xxx;
if (!init)
{
	// these are the differences at n+/- .5
	for (int i=nmin; i<=nmax; ++i) 
	{
		const D ivalf=D(i)+.5;
		const D ivali=D(i)-.5;
		//doen(lut[i-nmin],xxx,ivalf,ivali,0,false);
		doen(lut_erf[i-nmin],lut_exp[i-nmin],ivalf,ivali,0,false);
	}
	MM_ERR(" should only get here once ")

}
init=true;
IdxTy idx=n-nmin;
//if (n>=ntot) return lut[ntot-1];
if (n>=ntot)  idx=ntot-1;
//if (n<nmin) return lut[0];
if (n<nmin) idx=0;
//return lut[n-nmin];
erfi=lut_erf[idx];
expi=lut_exp[idx];
return 0;

}
/*
Evaluate a difference between two error functions and normals at xf and xi
derf=erf(xf)-erf(xi) and dnorm=exp(-xf*xf)-exp(-xi*xi)
See test code at bottom in main .
This is based on expanding an integral that evaluates to erf or normal
exaluated at xf and xi. The integral is expressed as a taylor series
at midpoint, (xf+xi)/2, The even terms cancel, (dx)^2-(-dx)^2 . As erf()
is the integral of exp^2 up to a scale factor, the terms
alternate between these two. The polynomials are predictable and only
a finite number likely needed it may be better to tabulate them  


*/

int doen( D & derf, D & dnorm, const D & xf, const D & xi, const IdxTy n, const bool debug=false)
{
//typedef unsigned int IdxTy;
typedef double D ;
//typedef std::string StrTy;
//typedef std::stringstream Ss;
//typedef mjm_integrals Myt;
typedef std::vector<D>  P;
P x1,x2;
//Myt mi;
const D delx=(xf-xi);
const D mean=.5*(xf+xi);
// this delta x is the full range NOT x2-xmean
if (false) { MM_MSG(" est nmax "<<MMPR3(delx,mean,(delx*mean*1.0))) } 
if (xf<=xi) { derf=0; dnorm=0; return 1; } 
Myt &  mi= *this;
D dxn=delx*.5;
D fac=1.0; // *dxn;
x1.push_back(dxn);
D sum=0;
D sumexp=0;
const D tol=1e-14;
const bool dump=debug; // (n==0);
const bool auto_stop=(n==0);
const int neff=(n==0)?1000:n;
for (int i=0; i<neff; ++i)
{
// include these in the poly
D coef=mi.evaluate_polynomial(x1,mean);
D term=coef; // *dxn/fac;
if ((i&1)==0) sum+=term; // these are the erf terms
else  // these are the exp terms 
// evntuall move this const division to the bottom 
//{ sumexp+=term*D(i+1)/dxn; } 
{ sumexp+=term*D(i+1)/dxn; } 
//fac=fac*(i+2);
//dxn=dxn*delx*.5;
if (dump) {
	MM_MSG(mi.print_polynomial(x1));
	MM_MSG(" terms "<<MMPR4(i,term,coef,dxn/fac)<<MMPR(sum))
	}
if (i>3) if (auto_stop)
{
// this only does one of the results should at least guess which one likely worse 
const D t=(term<0)?(-term):term;
const D v=(sum<0)?(-sum):sum;
if (t<(v*tol)) break;
}
//template<class Td>  void next_gauss_der_polynomial(Td & dest, const Td & v1)
mi.next_gauss_der_polynomial(x2,x1,dxn/(i+2.0));
x1=x2;
}
static const D ps=sqrt(M_PI);
const D expm=2.0*exp(-mean*mean);
//derf=2.0*sum*exp(-mean*mean)*2.0/ps;
derf=sum*expm*2.0/ps;
//dnorm=2.0*sumexp*exp(-mean*mean);
dnorm=sumexp*expm;
//derf=sum*erf(mean);
return 0;
}
////////////////////////////////////////////////////////////////////////////

// [1] http://www.mathematica-journal.com/2014/11/on-burmanns-theorem-and-its-application-to-problems-of-linear-and-nonlinear-heat-transfer-and-diffusion/#more-39602/ https://en.wikipedia.org/wiki/Error_function#cite_note-8
// http://mathworld.wolfram.com/BuermannsTheorem.html

template <class Tvec, class Tf>
 void burmann_ref1_eq21(Tvec & dest, Tf & func, const D & z0, const IdxTy mu, const IdxTy nu, const IdxTy kmax, const IdxTy flags=0)
{
m_cm.mark("eq21");
//const D z0=0;
for (IdxTy k=0; k<kmax; ++k)
{
D bc=0;
//typedef mjm_generic_itor<1> GenItor;
typedef mjm_fixed_sum_itor GenItor;
m_cm.mark("genitor");
GenItor gi(k,mu,nu);
m_cm.cumdel("genitortotal", "genitor");
const IdxTy sz=gi.k(); // gi.size();
while (gi.ok())
{
m_cm.mark("multi");
	D fac=gi.multiplicity();
m_cm.cumdel("multitotal", "multi");
	for (IdxTy s=1; s<=sz; ++s)
	{	
		//const D fsz0=func(s,z0);
	//	const D f0z0=func(0,z0);
//m_cm.inc("ratio");
//m_cm.mark("rat");
		const D rat=func.ratio(s,0,z0);
// this takes about 10 percent of the time lol
// at n-18, 11:20 reduced to 10:20 min:sec
//m_cm.cumdel("rattotal", "rat");
		//const D f1=(func(s,z0)/func(0,z0)/gi.factorial(s));
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
typedef mjm_fixed_sum_itor GenItor;
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
typedef mjm_fixed_sum_itor GenItor;
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




/////////////////////////////////////////////////////////

static D erf_burmann(const D & x, const IdxTy flags=0)
{
D res=0;
// for large x, this still has the diff of numbers close to 1 wtf 
// however it may be easier to fix 
const D theta=sqrt(1.0-exp(-x*x));
// [1] eqn 32 and prior equations:
// f(z)=sqrt(pi)/2*erf(z)=\int_0^z exp(-s*s)ds -> phi(z)=f'(z)=exp(-z*z)
// phi' = f'' = -2z  exp(-z*z)
// phi^(2n)(0)=2*(-1)^n\frac{(2n-1)!}{(n-1)!} ; \phi{2n+1) (0)=0 
// erf(x)=2/sqrt(pi)*(thata - (1/12)theta^3 - (7/480)theta^5
// - (5/896)theta^7 -(787/276480)theta^9... 
D phi_top=1;
D phi_bottom=1;
D sign=1;
D thetan=theta;
IdxTy n=0;
while (true)
{

//const D phi=2.0*sign*phi_top/phi_bottom;
// now back to (14) in [1]

++n;
phi_top=phi_top*(2*n-1);
phi_bottom=phi_bottom*(n-1);
sign=-sign;
thetan=thetan*theta;
}
res=res*2.0/sqrt(M_PI);
return res;
}
public:
template <class Ts> void dump_counter(const Ts & lbl )
{
    m_cm.dump(lbl , std::cout);
    m_cm.dump(lbl , std::cerr);
    MM_MSG(lbl<<CRLF<<m_cm.time_table("solver_times"))
    MM_ERR(lbl<<CRLF<<m_cm.time_table("solver_times"))
}


CounterMap m_cm;

////////////////////////////////////////////////////////////////////////////
}; // mjm_integrals


class shape_function_integrals
{
typedef shape_function_integrals Myt;
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
typedef mjm_generic_itor<2> LoopItor2;
typedef mjm_generic_itor<4> LoopItor4;
typedef mjm_generic_itor<3> LoopItor3;
// sf for node i in dir =0 is x and =1 is y
template <class Tp> void make_shape(Tp & shape, const IdxTy i, const IdxTy dir)
{
const IdxTy sf= (i<<1)+dir;
shape=Tp(2);
// in reality there are two types here that could be assigned based
// on bit patterns 
switch (sf)
{
// node 0 x comp
case 0: { shape[0]=1; shape[1]=-1; return ; }
// y comp
case 1: { shape[0]=1; shape[1]=-1; return ; }
// node 1
case 2: { shape[0]=0; shape[1]=1; return ; }
// y comp
case 3: { shape[0]=1; shape[1]=-1; return ; }
// node 2 x comp
case 4: { shape[0]=0; shape[1]=1; return ; }
// y comp
case 5: { shape[0]=0; shape[1]=1; return ; }
// node3 
case 6: { shape[0]=1; shape[1]=-1; return ; }
// y comp
case 7: { shape[0]=0; shape[1]=1; return ; }






} // switch
MM_ERR(" invalid shape request i= "<<i<<" dir "<<dir)
}


void laplacian()
{
typedef mjm_integrals Mi;
Mi mi;
// pairs of derivatives after integration by parts
//typedef double Tf;
typedef mjm_rational Tf;
typedef std::vector<Tf> Tpoly;
Tpoly destx,desty,x(2),y(2),xb,yb,x2i,y2i;
Tpoly x2,y2,x2b,y2b;
// loopitor? 
for (IdxTy i=0; i<4; ++i)
{
Tpoly res(4);
for (IdxTy j=0; j<4; ++j)
{
make_shape(x,i,0); make_shape(y,i,1);
make_shape(xb,j,0); make_shape(yb,j,1);

mi.differentiate_polynomial(x2i, x);
mi.differentiate_polynomial(x2b, xb);
mi.multiply_polynomials(x2,x2i,x2b);
mi.multiply_polynomials(y2,y,yb);
mi.integrate_polynomial(destx, x2);
mi.integrate_polynomial(desty, y2);
// the temopate needs to be fixed.. 
const Tf xx=mi.evaluate_polynomial(destx,Tf(1.0));
const Tf yy=mi.evaluate_polynomial(desty,Tf(1.0));
res[j]=xx*yy;
//MM_ERR(" results "<<i<<" "<<j<<" "<<xx<<" "<<yy<<" "<<(xx*yy))

} // j
MM_ERR(" lap "<<i<<" "<<res[0]<<" "<<res[1]<<" "<<res[2]<<" "<<res[3])
} // i 

} // laplacian

D grad_avg(const D & x1, const D & x2, const D & n1, const D & n2)
{
typedef mjm_integrals Mi;
Mi mi;
typedef D Tf;
typedef std::vector<Tf> Tpoly;
Tpoly s(2),sc(2),n(1),x,xc;
// sf and complementary
make_shape(s,0,0);
make_shape(sc,1,0);
n[0]=n1; 
mi.multiply_polynomials(x,s,n);
n[0]=n2;
mi.multiply_polynomials(xc,sc,n);
Tpoly xt;
mi.add_polynomials(xt,x,xc);
n[0]=x1;
mi.multiply_polynomials(x,s,n);
n[0]=x2;
mi.multiply_polynomials(xc,sc,n);
Tpoly nt;
mi.add_polynomials(nt,x,xc);
Tpoly ii;
mi.multiply_polynomials(ii,nt,xt);
Tpoly dest;
mi.integrate_polynomial(dest, ii);
const Tf yy=mi.evaluate_polynomial(dest,Tf(1.0));
return yy;
} // grad_avg


// integrals of current with grads as 3-0 and 2-1 with signs
// invertable with sign of d. Shape function integrals are normally
// 1/3 for sames and 1/6 for cross for quad elelemtn. 
// there is a balance between n3 and n0 equal to that of n2 and n1
// which specifies the position perpendicular to direction of integration
// which by default is the middle but coudl be moved 
MyBlock  j_mat( const D & s1, const D& s2, const D& s12, const D & d)
{
MyBlock mat(4,4);
const D f=.5/d;
const D S1=s1*f;
const D S2=s2*f;
const D S12=s12*f;
mat(0,0)=-S1; mat(0,1)=-S12; mat(0,2)=S12; mat(0,3)=S1;
mat(1,0)=-S12; mat(1,1)=-S2; mat(1,2)=S2; mat(1,3)=S12;
mat(2,0)=-S12; mat(2,1)=-S2; mat(2,2)=S2; mat(2,3)=S12;
mat(3,0)=-S1; mat(3,1)=-s12; mat(3,2)=S12; mat(3,3)=S1;

return mat;

}




void grad()
{
typedef mjm_integrals Mi;
Mi mi;
// pairs of derivatives after integration by parts
//typedef double Tf;
typedef mjm_rational Tf;
typedef std::vector<Tf> Tpoly;
Tpoly destx,desty,x(2),y(2),xb,yb,x2i,y2i;
Tpoly x2,y2,x2b,y2b;
// loopitor? 
for (IdxTy i=0; i<4; ++i)
{
Tpoly res(4);
for (IdxTy j=0; j<4; ++j)
{
make_shape(x,i,0); make_shape(y,i,1);
make_shape(xb,j,0); make_shape(yb,j,1);

//mi.differentiate_polynomial(x2i, x);
x2i= x;
mi.differentiate_polynomial(x2b, xb);
mi.multiply_polynomials(x2,x2i,x2b);
mi.multiply_polynomials(y2,y,yb);
mi.integrate_polynomial(destx, x2);
mi.integrate_polynomial(desty, y2);
// the temopate needs to be fixed.. 
const Tf xx=mi.evaluate_polynomial(destx,Tf(1.0));
const Tf yy=mi.evaluate_polynomial(desty,Tf(1.0));
res[j]=xx*yy;
//MM_ERR(" results "<<i<<" "<<j<<" "<<xx<<" "<<yy<<" "<<(xx*yy))

} // j
MM_ERR(" grad "<<i<<" "<<res[0]<<" "<<res[1]<<" "<<res[2]<<" "<<res[3])
} // i 

} // grad



#define SP_ITOR_I for (IdxTy i=0; i<x.m_n; ++i) {{ 
#define SP_ITOR_IJ for (IdxTy i=0; i<x.m_n; ++i) { for (IdxTy j=0; j<x.m_n; ++j) {
#define SP_ITOR_IJP for (IdxTy ip=0; ip<y.m_n; ++ip) { for (IdxTy jp=0; jp<y.m_n; ++jp) {
#define SP_ITOR_END }}

template <class Tv> class simple_polynomial : public std::vector<Tv>
{
typedef std::vector<Tv> Super;
typedef simple_polynomial<Tv> Myt;
public:
typedef typename Super::value_type value_type;
simple_polynomial(const IdxTy sz) : Super(sz) {}
simple_polynomial() : Super() {}
// users of this will not be able to accept std::vector AFAIK 
Myt & operator=(const value_type & v) 
{ Super::clear(); Super::push_back(v); return *this; }
bool operator==(const Myt & x)
{
const IdxTy sz=Super::size();
const IdxTy xsz=x.size();
IdxTy maxx=(sz>xsz)?sz:xsz;
while (maxx>0)
{
--maxx;
if (maxx>=sz) {if (x[maxx]!=0) { return false;  } else continue;}
if (maxx>=xsz){  if ((*this)[maxx]!=0) return false; else continue; } 
if ((*this)[maxx]!=x[maxx]) return false; 
}
return true; 
}

StrTy to_string(const IdxTy flags=0,const StrTy & var=StrTy("x")) const
{
mjm_integrals mi;
return mi.print_polynomial(*this,flags,var);

}

}; // simple_polynomial 

// a polynomial of x and y up to max degreee n-1. 
// use block matrix for later extension to multinomial 
// TODO FIXME this is almost always going to be sparse wtf 
template <class Tv> class dense_2_polynomial : public mjm_block_matrix<Tv>
{
typedef dense_2_polynomial<Tv> Myt;
typedef mjm_block_matrix<Tv> Super;
typedef Tv value_type;
public:
// this is the sie of the matrix, degree 0 to n-1
dense_2_polynomial(const IdxTy n): Super(n,n),m_n(n) {}
dense_2_polynomial(): Super(0),m_n(0) {}
dense_2_polynomial(const Super & m ): Super(m),m_n(m.size(0)) {}
Myt&  operator=(const Myt & that)
{
Myt&  x= *this; // (that.m_n);
x.m_n=that.m_n;
x.resize(x.m_n,x.m_n);
SP_ITOR_IJ
x(i,j)=that(i,j);
SP_ITOR_END
return x;
}

const IdxTy get_size() const { return m_n; } 
Myt&  set_size(const IdxTy n )
{
Myt&  x= *this; // (that.m_n);
x.m_n=n;
x.resize(x.m_n,x.m_n);
SP_ITOR_IJ
x(i,j)=0; // that(i,j);
SP_ITOR_END
return x;
}

//dense_2_polynomial(const IdxTy n): Super(n,n),m_n(n) {}
template <class Tpoly, class Tval > dense_2_polynomial(const Tpoly & px, const Tval  & cxp, const  Tval & cy , const Tval & c)
{
(*this)=  from_convolution( px, cxp, cy , c);
}
 dense_2_polynomial(const Tv  & cx, const  Tv & cy , const Tv & c)
:Super(3,3),m_n(3)
{
(*this)(1,0)=cx;
(*this)(0,1)=cy;
(*this)(0,0)=c;
}


template <class Tpoly > dense_2_polynomial(const Tpoly & px, const Tpoly & py)
{
(*this)=  from_outer_product( px,py);
}

template<class Td> static 
	void get_vector( Td & dest, const Myt & x, const IdxTy k, const IdxTy dir)
{
dest=Td(x.m_n);
SP_ITOR_I
if (dir==0) dest[i]=x(i,k);
else if (dir==1) dest[i]=x(k,i);
SP_ITOR_END

}
template<class Td> static 
	void set_vector( Myt & x, const Td & dest,  const IdxTy k, const IdxTy dir)
{
const IdxTy sz=dest.size();
//dest=Td(x.m_n);
SP_ITOR_I
if (i<sz)  {  
if (dir==0) x(i,k)=dest[i];
else if (dir==1) x(k,i)=dest[i];
} else 
{
if (dir==0) x(i,k)=0;
else if (dir==1) x(k,i)=0;

}
SP_ITOR_END

}
template<class Td> static 
	void add_vector( Myt & x, const Td & dest,  const IdxTy k, const IdxTy dir, const value_type & a)
{
//dest=Td(x.m_n);
const IdxTy sz=dest.size();
SP_ITOR_I
if (i<sz)  {  
if (dir==0) x(i,k)+=dest[i]*a;
else if (dir==1) x(k,i)+=dest[i]*a;
}
SP_ITOR_END

}



// this only works for same size, AFAICT but should
// check the block_matrix operator 
static void add( Myt & d, const Myt & x, const Myt & y)
{
const IdxTy szx=x.m_n;
const IdxTy szy=y.m_n;

const bool xbigger=(szx>szy);
const IdxTy szmax=(xbigger)?szx:szy;
const IdxTy szmin=(!xbigger)?szx:szy;
d= Myt(szmax);
//d=x+y;
if (xbigger)
{
SP_ITOR_IJ
d(i,j)=x(i,j);
if (i<szmin) if (j<szmin) d(i,j)+=y(i,j);
SP_ITOR_END
}
else
{
SP_ITOR_IJP
d(ip,jp)=y(ip,jp); // +y(ip,jp);
if (ip<szmin) if (jp<szmin) d(ip,jp)+=x(ip,jp);
SP_ITOR_END
}
}

 void accumulate(  const Myt & y)
{
Myt & x =(*this);
const IdxTy szx=x.m_n;
const IdxTy szy=y.m_n;

const bool xbigger=(szx>szy);
const IdxTy szmax=(xbigger)?szx:szy;
//d=x+y;
if (!xbigger)
{
SP_ITOR_IJ
x(i,j)+=y(i,j);
SP_ITOR_END
}
else
{
SP_ITOR_IJP
x(ip,jp)+=y(ip,jp);
SP_ITOR_END
}
}
template <class Tp> 
static void multiply( Myt & d, const Myt & x, const Tp & y, const bool dir, const value_type & scale=1)
{
const IdxTy szy=y.size();
if ((x.m_n==0) || (szy==0)) return; 
const IdxTy szd=1+(x.m_n-1)+(szy-1);
d= Myt(szd);
const bool dirx=dir;
//SP_ITOR_IJP
//for (IdxTy ip=0; ip<m_n; ++ip) for (IdxTy jp=0; jp<m_n; ++jp)
SP_ITOR_IJ
// TODO FIXME only do this if x(i,j) not zero 
for (IdxTy ip=0; ip<szy; ++ip)
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dirx) d(i+ip,j)+=x(i,j)*y[ip]*scale;
else d(i,j+ip)+=x(i,j)*y[ip]*scale;

//SP_ITOR_END
SP_ITOR_END

}




static void multiply( Myt & d, const Myt & x, const Myt & y, const value_type & scale=1)
{
if ((x.m_n==0) || (y.m_n==0)) return; 
const IdxTy szd=1+(x.m_n-1)+(y.m_n-1);
d= Myt(szd);
SP_ITOR_IJP
//for (IdxTy ip=0; ip<m_n; ++ip) for (IdxTy jp=0; jp<m_n; ++jp)
// TODO FIXME only do this if x(ip,jp) not zero 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
d(i+ip,j+jp)+=x(i,j)*y(ip,jp)*scale;

SP_ITOR_END
SP_ITOR_END

}




// considering x as a polynomial in variables x and y, replace the
// powers of "y" in x with the definition in y as a polynomial in 2 variables
// on or which is the x variable in x lol. 
template <class Txx> static Txx cpo(const Txx & v, const IdxTy n)
{
Txx r=1;
IdxTy i=n;
while (i>0){  r=v*r; --i; }  
return r;
}
static void composite( Myt & d, const Myt & x, const Myt & y, const IdxTy flags=1)
{
if ((x.m_n==0) || (y.m_n==0)) return; 
// may need to trim later but the biggest dimension is 
// x(m_n-1,m_n-1) and y(m_n-1,m_n-1)
const IdxTy szd=1+(x.m_n-1)*(y.m_n-1); // 1+(x.m_n-1)+(y.m_n-1);
const IdxTy imax=x.m_n;
const IdxTy jmax=x.m_n;
d= Myt(szd);
Myt yn=Myt(1);
yn(0,0)=1;
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
// the first index of "x" is for the variable being replcaed by y
const IdxTy xyz=flags;
switch (xyz)
{
case 0://{ SP_ITOR_IJP SP_ITOR_IJ d(i*ip,j+jp*i)+=x(i,j)*cpo(y(ip,jp),i); SP_ITOR_END SP_ITOR_END break; } 
{
//for (IdxTy j=0; j<jmax; ++j) {d(0,j)=x(0,j); }
for(IdxTy i=0; i<imax; ++i)
{
const IdxTy ipmax=yn.m_n;
const IdxTy jpmax=yn.m_n;

for (IdxTy j=0; j<jmax; ++j)
{
// TODO FIXME only for yn!=0 or just make these all SPARSE 
for (IdxTy ip=0; ip<ipmax; ++ip)
for (IdxTy jp=0; jp<jpmax; ++jp)

// for each row just do multiply times the thing to n power
// not tested, the dest power is just the yn power and NOT the sum as it is replcaced
d(ip,j+jp)+=x(i,j)*yn(ip,jp);

} // j 
//if ( i!=0) 
{ Myt ynext; multiply(ynext,yn,y); yn=ynext;}

} // i 
break; 
}
case 1://{SP_ITOR_IJP SP_ITOR_IJ d(i+j*ip,j*jp)+=x(i,j)*cpo(y(ip,jp),j);SP_ITOR_END SP_ITOR_END break;} 
{
//for (IdxTy j=0; j<jmax; ++j) {d(j,0)=x(j,0); }
for(IdxTy i=0; i<imax; ++i)
{
const IdxTy ipmax=yn.m_n;
const IdxTy jpmax=yn.m_n;

for (IdxTy j=0; j<jmax; ++j)
{
for (IdxTy ip=0; ip<ipmax; ++ip)
for (IdxTy jp=0; jp<jpmax; ++jp)
// this may work, the dest is NOT the sum of sources but is replcaed by the y power
// for each row just do multiply times the thing to n power
d(j+jp,ip)+=x(j,i)*yn(jp,ip);

} // j 
//if ( i!=0) 
{ Myt ynext; multiply(ynext,yn,y); yn=ynext;}

} // i 
break; 
}

default : MM_ERR(" bad flag in composite "<<MMPR(xyz))
}; // switch 

}

static void integrate( Myt & d, const Myt & x, const IdxTy dir)
{
// this is stupid to have SQUARE matrix when you do this crap 
const IdxTy szd=1+(x.m_n);
d= Myt(szd);
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dir==0) d(i+1,j)=x(i,j)/(i+1);
else if (dir==1) d(i,j+1)=x(i,j)/(j+1);
SP_ITOR_END
}
static void differentiate( Myt & y, const Myt & x, const IdxTy dir)
{
// this is stupid to have SQUARE matrix when you do this crap 
if (x.m_n==0) { y=Myt(0); return; } 
const IdxTy szd=(x.m_n)-0;
const IdxTy szmax=(x.m_n)-1;
y= Myt(szd);
SP_ITOR_IJP
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dir==0) {if (ip<szmax) { y(ip,jp)=(ip+1.0)*x(ip+1,jp); }}
else {if (dir==1) {if (jp<szmax) y(ip,jp)=(jp+1.0)*x(ip,jp+1);}}
SP_ITOR_END
}




// this probably does not work right, use _new
static void integrate_e( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
MM_ONCE(" obsolete make work with _e_new", )
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(sz);

// this is stupid to have SQUARE matrix when you do this crap 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
	const IdxTy loc=(dir==0)?i:j;
	const IdxTy ord=(dir==0)?j:i;
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	// note that these could be cached but pretty easy to compute 
	//D res=mi.x_n_exp_integral(p, n, b);
	mi.x_n_erf_integral (perf, p,  ord, b,false);
	//derf.add_vector(derf,perf,loc,(dir==1)?1:0,x(i,j));
	derf.add_vector(derf,perf,ord,(dir==1)?1:0,x(i,j));
	//derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,dir,x(i,j));
//	dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	dexp.add_vector(dexp,p,ord,(dir==1)?1:0,x(i,j));
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	for (IdxTy ii=0; ii<p.size(); ++ii) p[ii]=0;
	D res=mi.x_n_exp_integral(p, ord, b,false);
	perf[0]=res*x(i,j);
//	mi.x_n_erf_integral (perf, p,  n, b);
//	derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	// dir==0 means add up and down  vary fiest index  i order is i  
	dexp.add_vector(dexp,p,ord,(dir==1)?1:0,x(i,j));
	//dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	if (dir==1) derf(i,0)+=perf[0];
	else derf(0,j)+=perf[0];

} // dir==1 
SP_ITOR_END
}
///////////////////////////////////////////////////////////

// this appears to work in the one case I tested
static void integrate_e_new( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(erfsz);
const IdxTy dirv=dir;
const bool rev=(dirv==0);
for (IdxTy i=0; i<sz; ++i)
{
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	mi.x_n_erf_integral (perf, p, i, b,false);
	MM_MSG(MMPR4(perf.size(),p.size(),derf.m_n,dexp.m_n))
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		derf.add_vector(derf,perf,j,dirv,v);
		dexp.add_vector(dexp,p,j,dirv,v);
	} // j 
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	D res=mi.x_n_exp_integral(p,i, b,false);
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		perf[0]=res*v;
		dexp.add_vector(dexp,p,j,dirv,v);
		if (dir==1) derf(j,0)+=perf[0];
		else derf(0,j)+=perf[0];
	}

} // dir==1 

} // i

}
/////////////////////////////////////////////////////////////////////
// user the rationzlied d "rerf" function which does integration ignoring the
// conversion of 2sqrt(b)/sqrt(pi) but user needs to add that back before calling erf 
static void integrate_e_rnew( Myt & dexp,Myt & derf, const Myt & x, const value_type & b
, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
const IdxTy sz=x.m_n; // this means that the max power is (m_n-1) 
// for erf the power can go up one in erf but not in exp as int by parts
// puts (n+1)  pwer in integrand going back to n after integration 
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(erfsz);
const IdxTy dirv=dir;
const bool rev=(dirv==0);
for (IdxTy i=0; i<sz; ++i)
{
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	mi.x_n_rerf_integral (perf, p, i, b,false);
	// TODO this can not be trimmed right away as it could be accumulating although
	// should fix it 
	if (false) MM_MSG(MMPR4(perf.size(),p.size(),derf.m_n,dexp.m_n))
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		derf.add_vector(derf,perf,j,dirv,v);
		dexp.add_vector(dexp,p,j,dirv,v);
	} // j 
	}

	else if (factor==1)
{ // integrating the exp only produces a single rerf term 
//	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	D res=mi.x_n_rexp_integral(p,i, b,false);
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		const D perf0=res*v; // 	perf[0]=res*v;
		dexp.add_vector(dexp,p,j,dirv,v);
		if (dir==1) derf(j,0)+=perf0; // perf[0];
		else derf(0,j)+=perf0; // perf[0];
	}

} // dir==1 

} // i

}
/////////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////

static void integrate_e_old( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
MM_ONCE(" obsolete make work with _e_new", )
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(sz);

// this is stupid to have SQUARE matrix when you do this crap 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
	const IdxTy loc=(dir==0)?i:j;
	const IdxTy ord=(dir==0)?j:i;
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	// note that these could be cached but pretty easy to compute 
	//D res=mi.x_n_exp_integral(p, n, b);
	mi.x_n_erf_integral (perf, p,  ord, b,false);
	derf.add_vector(derf,perf,loc,(dir==0)?1:0,x(i,j));
	//derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,dir,x(i,j));
	dexp.add_vector(dexp,p,loc,dir,x(i,j));
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	for (IdxTy ii=0; ii<p.size(); ++ii) p[ii]=0;
	D res=mi.x_n_exp_integral(p, ord, b,false);
	perf[0]=res*x(i,j);
//	mi.x_n_erf_integral (perf, p,  n, b);
//	derf.add_vector(derf,perf,loc,dir,x(i,j));
	dexp.add_vector(dexp,p,loc,dir,x(i,j));
	if (dir==0) derf(i,0)+=perf[0];
	else derf(0,j)+=perf[0];

} // dir==1 
SP_ITOR_END
}






// evaluate the 2 nomial at value v for one of the variable and 
// return a polynomial in dest
template<class Td> static void evaluate( Td & dest, const Myt & x, const value_type & v, const IdxTy dir)
{
dest=Td(x.m_n);
value_type vp[x.m_n];
vp[0]=1;
SP_ITOR_IJ
if (dir==0) dest[i]+=x(i,j)*vp[i];
else if (dir==1) dest[j]+=x(i,j)*vp[j];
SP_ITOR_END

}

value_type evaluate(  const value_type & a, const value_type & b ) const
{
const Myt & x=(*this);
const IdxTy sz=x.m_n;
if (sz==0) return 0; 
value_type vx[sz],vy[sz];
vx[0]=1; vy[0]=1;
for (IdxTy i=1; i<sz; ++i) { vx[i]=vx[i-1]*a; vy[i]=vy[i-1]*b; }
// TODO  consider making the power tables exportable 
value_type res=0;
SP_ITOR_IJ
res+=x(i,j)*vx[i]*vy[j];
SP_ITOR_END
return res;


}
// these are actually backwards but ok 
value_type sum_squares(  const value_type & x0, const value_type & x1
,  const value_type & xp0, const value_type & xp1,const value_type & d
) const 
{
value_type res=0;
value_type v=evaluate(x0-d,xp0-d); res+=v*v;
v=evaluate(x1-d,xp0-d); res+=v*v;
v=evaluate(x1-d,xp1-d); res+=v*v;
v=evaluate(x0-d,xp1-d); res+=v*v;
return res;

}
value_type evaluate_four(  const value_type & x0, const value_type & x1
,  const value_type & xp0, const value_type & xp1,
  const value_type & e00, const value_type & e11,
  const value_type & e01, const value_type & e10
, IdxTy * pflags=0
) const 
{
//value_type res=0;
//value_type v=evaluate(x0-d,xp0-d); res+=v*v;
//v=evaluate(x1-d,xp0-d); res+=v*v;
//v=evaluate(x1-d,xp1-d); res+=v*v;
//v=evaluate(x0-d,xp1-d); res+=v*v;
const Myt & x=(*this);
const IdxTy sz=x.m_n;
if (sz==0) return 0; 
bool no_x=true;
value_type vx[sz],vy[sz],vx2[sz],vy2[sz];
vx[0]=1; vy[0]=1;
vx2[0]=1; vy2[0]=1;
for (IdxTy i=1; i<sz; ++i) { vx[i]=vx[i-1]*x0; vy[i]=vy[i-1]*xp0; }
for (IdxTy i=1; i<sz; ++i) { vx2[i]=vx2[i-1]*x1; vy2[i]=vy2[i-1]*xp1; }
// TODO  consider making the power tables exportable 
value_type res=0;
SP_ITOR_IJ
if (x(i,j)==0) continue;
if (i!=0) if (j!=0)  no_x=false;
res+=x(i,j)*(vx[i]*vy[j]*e00+vx2[i]*vy2[j]*e11-vx[i]*vy2[j]*e10-vx2[i]*vy[j]*e01);
SP_ITOR_END
if (pflags!=0) *pflags=(no_x)?1:0;
return res;

}



// evaluate with a polynomial in the other variable, 
// the resulting polynomial dest is a higher order
template<class Td> static void evaluate( Td & dest, const Myt & x, const Td & vpoly, const IdxTy dir)
{
const IdxTy szx=x.m_n;
const IdxTy szp=vpoly.size();
if (szx==0) return;
if (szp==0) return;
const IdxTy sz=1+(szx-1)*(szp-1);
dest= Td(sz);
if (sz<1) return;
mjm_integrals mi;
const IdxTy vpoly_size=vpoly.size();
dest=Td(sz);
//dest[0]=x(0,0);
//value_type vp[x.m_n];
//vp[0]=1;
// composite
for (IdxTy i=0; i<szx; ++i) { 
Td fn; // =vpoly;
fn.push_back(1);
for (IdxTy k=0; k<fn.size(); ++k) {
for (IdxTy j=0; j<szx; ++j) {
if (dir==0) dest[i+k]+=x(i,j)*fn[i];
else if (dir==1) dest[i+k]+=x(i,j)*fn[j];
}
}
Td next;
mi.multiply_polynomials(next,vpoly,fn);
fn=next;

}
}

template <class Tpoly > static Myt from_outer_product(const Tpoly & px, const Tpoly & py)
{
const IdxTy szx=px.size();
const IdxTy szy=py.size();
const IdxTy sz=(szx>szy)?szx:szy;
//Tval cp[szx]; cp[0]=1;
Myt mt(sz);
Myt & x=mt;
SP_ITOR_IJ
if ((i>=szx) ||(j>=szy)) continue;
mt(i,j)=px[i]*py[j];
SP_ITOR_END
return mt;
}

// make a _2_poly from a polynomail p(x) and composite with 
// x= cxp*x' + cy*y+ c 
template <class Tpoly, class Tval > static Myt from_convolution(const Tpoly & px, const Tval  & cx, const  Tval & cy , const Tval & c)
{
const IdxTy szx=px.size();
//Tval cp[szx]; cp[0]=1;
Tval cxp[szx], cyp[szx],cp[szx];
Myt mt(szx);
cxp[0]=1;  cyp[0]=1; cp[0]=1;
for (IdxTy i=1; i<szx; ++i) 
	{cxp[i]=cxp[i-1]*cx; cyp[i]=cyp[i-1]*cy; cp[i]=cp[i-1]*c; }
Myt & x=mt;
SP_ITOR_IJ
for (IdxTy ix=(i+j); ix<szx; ++ix) 
{
mt(i,j)+=px[ix]*cxp[i]*cyp[j]*cp[ix-i-j];

} // ix
SP_ITOR_END

return mt;
}




static void coef_table( Super & table, const IdxTy szx, const value_type & cx, 
	const value_type & cy, const value_type & c)
{
Super m(szx,szx,szx);
//typedef value_type Tval;
trinomial_coef(m,szx);
//Tval cx[szx], cyp[szx],cp[szx];
// any of these could be zero 
//cxp[0]=1;  cyp[0]=1; cp[0]=1;
//for (IdxTy i=1; i<szx; ++i) 
//	{cxp[i]=cxp[i-1]*cx; cyp[i]=cyp[i-1]*cy; cp[i]=cp[i-1]*c; }


table=m;
}
static IdxTy bffac(IdxTy n) { IdxTy i=1; while (n>0) { i=i*n; --n; }return i; }
static void trinomial_coefs(Super & t, const IdxTy sz)
{
const IdxTy szi=sz;
const IdxTy szj=sz;
const IdxTy szk=sz;
value_type fac=1;
for (IdxTy i=0; i<szi; ++i) { 
//fac=1;
for (IdxTy j=0; i<szj; ++j) { 
fac=1;
for (IdxTy k=0; i<szk; ++k) { 
t(i,j,k)=fac;
IdxTy bf=bffac(i)/(bffac(j)*bffac(k)*bffac(i-j-k));
if (t(i,j,k)!=bf) MM_ERR(MMPR4(t(i,j,k),bf,i,j)<<MMPR(k))
fac=fac*(i-j-k)/(k+1);
}
}
}


}

template<class Tx, class Ty, class Tn > bool  add_term( Tx & terms, Ty & ss, const Tn & t, const bool force=!true) const 
{
	if ((t!=0)||(force)) { if (terms!=0) { if ( t>=0) ss<<"+"; } ss<<t; ++terms; return true;  }
	return false; 
}
//StrTy to_string( const IdxTy flags=0)
StrTy to_string(const IdxTy flags=0,const StrTy & var=StrTy("x"),
	const StrTy & var2=StrTy("y")) const
{
std::stringstream ss;
//const int x2s=x2.size();
const Myt & x=*this;
int terms=0;
const IdxTy sz=x.m_n;
SP_ITOR_IJ
const IdxTy ii=sz-i-1;
const IdxTy jj=sz-j-1;
const value_type & t=x(ii,jj);
 if ( add_term(terms,ss,x(ii,jj)))   
{
//if (jj!=0) {if (jj==1) {ss<<"x";} else   {ss<<"x^"<<jj;} } 
if (ii!=0) {ss<<var; if (ii!=1) {ss<<"^"<<ii;}  } 
if (jj!=0) {ss<<var2; if (jj!=1) {ss<<"^"<<jj;}  } 
//if (ii!=0) {if (ii==1) {ss<<"y";} else { ss<<"y^"<<ii;} } 

} 
//if (x2s>0) {if (x2[0]!=0) {  if (terms!=0) ss<<"+"; ss<<x2[0];} }
//if (x2s>0) {add_term(terms,ss,x2[0]); }
SP_ITOR_END
return ss.str();
}

private:
IdxTy m_n;

}; // dense_2_polynomial

// rational function (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x))
template <class Tv=double> class norm_rat
{

class Tr{
public:
//typedef double D;
typedef Tv  D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
//typedef LocTy location_type;
//typedef mjm_rational RatTy;
}; // Tr

typedef typename Tr::D D;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::StrTy StrTy;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<IdxTy> MyBlockInt;
typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
typedef mjm_generic_itor<2> LoopItor2;
typedef mjm_generic_itor<3> LoopItor3;
typedef mjm_generic_itor<4> LoopItor4;

typedef norm_rat<Tv>  Myt;
typedef Tv value_type;
typedef simple_polynomial<Tv> Tpoly;
typedef mjm_integrals Mi;
public:
norm_rat() { Init(); Default(); }
norm_rat(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S) 
:m_P(P),m_Q(Q),m_R(R),m_S(S)
{ Init();  }
norm_rat(const value_type & x) {  from_const(x); } 

void from_const(const value_type & x)
{
Mono(m_P,0,0); Mono(m_Q,0,x); Mono(m_R,0,0); Mono(m_S,0,1);
}
void normal()
{
Mono(m_P,0,1); Mono(m_Q,0,0); Mono(m_R,0,0); Mono(m_S,0,1);
}


void norm_minus_1_over_xn(const IdxTy n)
{
Mono(m_P,0,1); Mono(m_Q,0,-1); Mono(m_R,0,0); Mono(m_S,n,1);

}
void norm_minus_1_over_xn(const IdxTy n, const D & z0)
{
Mono(m_P,0,1); Mono(m_Q,0,-D(exp(-double(z0*z0)))); Mono(m_R,0,0); Mono(m_S,n,1); m_S[0]=-z0;
}
void norm_minus_1_over_xn_exact(const IdxTy n, const D & z0,const D & zc)
{
Mono(m_P,0,1); Mono(m_Q,0,z0); Mono(m_R,0,0); Mono(m_S,n,1); m_S[0]=zc;
}


bool is_zero(const Tpoly & x)
{
const IdxTy sz=x.size();
for (IdxTy i=0; i<sz; ++i) if (x[i]!=0) return false;
return true;

}
Myt operator*(const value_type & that) //const
{
Myt x;
Tpoly x1;
x1.push_back(that);
mi.multiply_polynomials(x.m_Q,m_Q,x1); // zero E numerator
mi.multiply_polynomials(x.m_P,m_P,x1); // zero E numerator
//MM_MSG(MMPR3(mi.print_polynomial(m_Q),mi.print_polynomial(x1),mi.print_polynomial(x.m_Q)))
//MM_MSG(MMPR3(mi.print_polynomial(m_P),mi.print_polynomial(x1),mi.print_polynomial(x.m_P)))

x.m_R=m_R;
x.m_S=m_S;
mi.trim_polynomial(x.m_P);
mi.trim_polynomial(x.m_Q);
mi.trim_polynomial(x.m_R);
mi.trim_polynomial(x.m_S);
return x; 
}



Myt operator/(const Myt & that) //const
{
Myt x;
// (kkkkkkkkk
Tpoly x1,x2;
mi.multiply_polynomials(x.m_Q,m_Q,that.m_S); // zero E numerator
mi.multiply_polynomials(x1,m_Q,that.m_R); //  E numerator
mi.multiply_polynomials(x2,m_P,that.m_S); //  E numerator
mi.add_polynomials(x.m_P,x1,x2);
// need to check E^2 terms are zed...
mi.multiply_polynomials(x.m_S,that.m_Q,m_S); // zero E numerator
mi.multiply_polynomials(x1,that.m_Q,m_R); //  E numerator
mi.multiply_polynomials(x2,that.m_P,m_S); //  E numerator
mi.add_polynomials(x.m_R,x1,x2);
mi.multiply_polynomials(x1,that.m_P,m_R); //  E^2 denr
mi.multiply_polynomials(x1,that.m_R,m_P); //  E^2 numerator
if ( is_zero(x.m_S) && is_zero(x.m_Q))
{
x.m_S= x.m_R; x.m_Q=x.m_P;
x.m_R=x1; x.m_P=x2;
}

// no idea why these were getting into the hundreds or more wtf TODO FIXME 
mi.trim_polynomial(x.m_P);
mi.trim_polynomial(x.m_Q);
mi.trim_polynomial(x.m_R);
mi.trim_polynomial(x.m_S);
//MM_MSG(MMPR4(x.m_P.size(),x.m_Q.size(),x.m_R.size(),x.m_S.size())<<MMPR(x.to_string()))

return x;
}

// derivative and put the higher powers into e2 for top and bottom 
/*
TODO FIXME FUDD see if this should work cancelling exp terms lol 
mjm_burmann.h350  bcc.to_string()=(1)*exp(-x*x)+0
----------------------------------------------------------------
(-2x)*exp(-x*x)+

mjm_burmann.h352  f.der(bcc).to_string()=()*exp(-x*x)+
----------------------------------------------------------------
()*exp(-x*x)+

*/

 void der( Myt & y, Tpoly & e2n, Tpoly & e2d)
{
// keep within this constraint, punt the higher order polys
// (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x)) ->
// (N'D- D'N)/D^2  from quotient rule.
// N' = P'E-2xPE+Q', D'=R'E-2xRE+S'
Tpoly n1,d1,n2,d2,t1,t2,t3,t4;
mi.next_gauss_der_polynomial(n1,m_P);
mi.next_gauss_der_polynomial(d1,m_R);
mi.differentiate_polynomial(n2,m_Q);
mi.differentiate_polynomial(d2,m_S);
mi.multiply_polynomials(y.m_S,m_S,m_S);
mi.multiply_polynomials(t1,m_R,m_S);
mi.add_polynomials(y.m_R,t1,t1,1.0,1.0);
mi.multiply_polynomials(e2d,m_R,m_R);
// N' is n1*E+n2 and D' is d1*E+d2
// N = PE+Q, D=RE+S
// N'D= E(R*n2+S*n1)+n2*S+EER*n1 
// ND'= EE*(d1*P) + E*(d1*Q+d2*P) +Q*d2  
// N'D- ND'
mi.multiply_polynomials(t1,d1,m_P);
mi.multiply_polynomials(t2,n1,m_R);
mi.add_polynomials(e2n,t2,t1,1.0,-1.0);
mi.multiply_polynomials(t1,d1,m_Q);
mi.multiply_polynomials(t2,d2,m_P);
mi.add_polynomials(t3,t1,t2); //  ND' E term 
mi.multiply_polynomials(t1,n1,m_S);
mi.multiply_polynomials(t2,n2,m_R);
mi.add_polynomials(t4,t1,t2); //  N'D E term 
mi.add_polynomials(y.m_P,t3,t4,-1.0,1.0); //  N'D E term 

mi.multiply_polynomials(t1,d2,m_Q);
//mi.multiply_polynomials(y.m_Q,d2,m_Q);
mi.multiply_polynomials(t2,n2,m_S);
mi.add_polynomials(y.m_Q,t1,t2,-1.0,1.0); //  N'D E term 
if ( is_zero(y.m_S) && is_zero(y.m_Q))
{
y.m_S= y.m_R; y.m_Q=y.m_P;
y.m_R=e2d; y.m_P=e2n;
e2d.clear();
e2n.clear();
}
if ( !is_zero(e2n) ||!is_zero(e2d))
{
MM_MSG(" der genearated e 2 terms "<<MMPR2(e2n.to_string(),e2d.to_string()))
}
//template<class Td>  void polynomial_list_reduce(Td & dest, const Td & v1)
std::vector<Tpoly> dest,v1;
v1.push_back(y.m_P);
v1.push_back(y.m_Q);
v1.push_back(y.m_R);
v1.push_back(y.m_S);
mi.polynomial_list_reduce( dest,  v1);
y.m_P=dest[0];
y.m_Q=dest[1];
y.m_R=dest[2];
y.m_S=dest[3];

}
// this is going to get evaluated at zero's of the denominator,
// however l'hopital can save it 
// (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x)) ->
Tv eval(const Tpoly & p, const Tpoly & q, const Tv & x, const Tv & ex)
{
return mi.evaluate_polynomial(p,x)*ex+
			mi.evaluate_polynomial(q,x);
}
Tv evaluate(const Tv & x) //const
{
const Tv ex=exp(-x*x);
return evaluate(x,ex);
}

Tv evaluate(const Tv & x, const Tv & ex) //const
{
const Tv n=eval(m_P,m_Q,x,ex); 
const Tv d=eval(m_R,m_S,x,ex); 
// ok this is not right lol as small non zero prevent l'hopital 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( d!=0) return (n/d);
Tpoly P,Q,R,S;
Tpoly P2,Q2,R2,S2;
P=m_P; Q=m_Q; R=m_R; S=m_S;
P2=P; Q2=Q; R2=R; S2=S;

while (true)
{
mi.next_gauss_der_polynomial(P2,P);
mi.differentiate_polynomial(Q2,Q);
mi.next_gauss_der_polynomial(R2,R);
mi.differentiate_polynomial(S2,S);
const Tv n=eval(P2,Q2,x,ex); 
const Tv d=eval(R2,S2,x,ex); 
//	MM_MSG(MMPR(print(P,Q,R,S)))	
//	MM_MSG(MMPR4(print(P2,Q2,R2,S2),x,n,d))	
if (d!=0)  return (n/d);
// TODO decide if you wan to do this FIXME 
//if ((n!=0)||(d!=0)) return (n/d);
// return the nan apropos for Tv 
// not sure cycles are possible with the exp poly 
if ((P2==P) &&(Q2==Q) && (R2==R)&&(S2==S)) return (n/d);
P=P2; Q=Q2; R=R2; S=S2;
} // while 
// not reachable 
return 0; 

} // evaluate

Tv evaluate_at_zero() //const
{
// FIXME TODO this will bomb if the sizes are zero... 
const IdxTy szp=m_P.size();
const IdxTy szq=m_Q.size();
const IdxTy szr=m_R.size();
const IdxTy szs=m_S.size();
const Tv p0=(szp>0)?m_P[0]:0;
const Tv q0=(szq>0)?m_Q[0]:0;
const Tv r0=(szr>0)?m_R[0]:0;
const Tv s0=(szs>0)?m_S[0]:0;
const Tv n=p0+q0;  
const Tv d=r0+s0; // eval(m_R,m_S,x,ex); 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( d!=0) return (n/d);
//Tpoly expn,expd,temp,neff,deff;
// TODO FIXME need to find real limits ... 
// int nterms=int(szq)-int(szp);
//	if (nterms<10) nterms=10;
// int dterms=int(szs)-int(szr);
//	if (dterms<10) dterms=10;
//mi.gauss_taylor_polynomial(expn, nterms);
//mi.gauss_taylor_polynomial(expd, dterms);
//MM_MSG(MMPR(mi.print_polynomial(expn)))
//MM_MSG(MMPR(mi.print_polynomial(expd)))
//mi.multiply_polynomials(temp,expn,m_P);
//mi.add_polynomials(neff,temp,m_Q);
//mi.multiply_polynomials(temp,expd,m_R);
//mi.add_polynomials(deff,temp,m_S);
// TODO FIXME KLUGE fix the limit lol 
const IdxTy szdeff=szp+szq+10 ; // neff.size();
const IdxTy szneff=szr+szs+10  ; // deff.size();

IdxTy i=0;
while (true)
{
Tv n=0;
Tv d=0; 
if (i>=szneff) if (i>=szdeff)
{
MM_ERR(" no limit found "<<MMPR4(szp,szq,szr,szs)<<MMPR4(szneff,szdeff,mi.print_polynomial(m_P),mi.print_polynomial(m_R)))
 return (n/d);
}
d=mi.gauss_term_polynomial( m_R,i );
const Tv si=(szs>i)?m_S[i]:0;
d+=si;
// this is not really right as 1/0 is a problem however that should
// not occur but the user needs to be careful.
// Danger Will Robinson !!!!!! FIXED 
// and it is no faster anyway 
//if (d!=0)
{
n=mi.gauss_term_polynomial( m_P,i );
const Tv qi=(szq>i)?m_Q[i]:0;
n+=qi;
}

//if (i<szneff) n=neff[i];
//if (i<szdeff) d=deff[i];
if (d!=0) return (n/d);
++i;
} // true 
MM_ERR(" should not get here ")
return 0; // unreadable
} // evaluate_at_zero

StrTy print(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S)
{
SsTy ss;
const StrTy exps="exp(-x*x)";
ss<<"("<<mi.print_polynomial(P)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(Q)<<CRLF;
ss<<"----------------------------------------------------------------"<<CRLF;
ss<<"("<<mi.print_polynomial(R)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(S)<<CRLF;

return ss.str();
}

StrTy to_string()  { return print(m_P,m_Q,m_R,m_S); }

private:

const Tpoly &  Init()
{
Mono(m_default,0,1);
return m_default;
}
void Default() { 
Default(m_P); Default(m_Q); Default(m_R); Default(m_S);
 }
void Mono(Tpoly & p, const IdxTy n, const Tv & v)
{
p = Tpoly(n+1);
p[n]=v;
}

void Default(Tpoly & x) { x=m_default; }
Tpoly m_default;
// rational function (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x))
Tpoly m_P,m_Q,m_R,m_S;
Mi mi;

}; // norm_rat


// 2d shape function 
class sf2d
{
typedef sf2d Myt;
public:
typedef mjm_rational Tf;
//typedef std::vector<Tf> Tpoly;
typedef simple_polynomial<Tf> Tpoly;

template <class Ti> 
void make_sf(Ti & mi, const IdxTy i)
{ mi.make_shape(x,i,0); mi.make_shape(y,i,1); }

template <class Ti> void differentiate_x(Ti & mi )
{
Tpoly d;
mi.differentiate_polynomial(d, x);
x=d;
}

template <class Ti> void differentiate_y(Ti & mi )
{
Tpoly d;
mi.differentiate_polynomial(d, y);
y=d;
}

template <class Ti> void integrate(Ti & mi )
{
Tpoly d;
mi.integrate_polynomial(d, x);
x=d;
mi.integrate_polynomial(d, y);
y=d;
}

template <class Ti> void integrate_x(Ti & mi )
{
Tpoly d;
mi.integrate_polynomial(d, x);
x=d;
}

template <class Ti> void integrate_y(Ti & mi )
{
Tpoly d;
mi.integrate_polynomial(d, y);
y=d;
}


template <class Ti> Myt times(Ti & mi,const Myt & that )
{
Myt x;
mi.multiply_polynomials(x.x,(*this).x,that.x);
mi.multiply_polynomials(x.y,(*this).y,that.y);
return x; 
}

// this works as long as something has been integrated lol 
template <class Ti>  Tf evaluate(Ti & mi )
{
const Tf xx=mi.evaluate_polynomial(x,Tf(1.0));
const Tf yy=mi.evaluate_polynomial(y,Tf(1.0));
return xx*yy;
}
template <class Ti>  Tf evaluate(Ti & mi, const Tf & _x, const Tf & _y )
{
const Tf xx=mi.evaluate_polynomial(x,Tf(_x));
const Tf yy=mi.evaluate_polynomial(y,Tf(_y));
return xx*yy;
}




Tpoly x,y;

};// sf2d

void lap2()
{
typedef mjm_integrals Mi;
Mi mi;
LoopItor2 itor(4,4);
sf2d si,sj,sp;
std::vector<sf2d::Tf> ans(4);
while (itor)
{
si.make_sf(*this,itor[0]);
sj.make_sf(*this,itor[1]);
si.differentiate_x(mi);
sj.differentiate_x(mi);
sp=si.times(mi,sj);
sp.integrate(mi);
ans[itor[1]]=sp.evaluate(mi);
if ( itor[1]==3)  {
	MM_ERR(" lap2 "<<itor.string()<<" "<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<" "<<ans[3]) }

++itor;
} // itor

} // lap2



void lap2y()
{
typedef mjm_integrals Mi;
Mi mi;
LoopItor2 itor(4,4);
sf2d si,sj,sp;
std::vector<sf2d::Tf> ans(4);
while (itor)
{
si.make_sf(*this,itor[0]);
sj.make_sf(*this,itor[1]);
si.differentiate_y(mi);
sj.differentiate_y(mi);
sp=si.times(mi,sj);
sp.integrate(mi);
ans[itor[1]]=sp.evaluate(mi);
if ( itor[1]==3)  {
	MM_ERR(" lap2y "<<itor.string()<<" "<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<" "<<ans[3]) }

++itor;
} // itor

} // lap2y

//////////////////////////////////////////////////////////


void gradnewy()
{
typedef mjm_integrals Mi;
mjm_block_matrix<sf2d::Tf> mat(4,4);
Mi mi;
LoopItor2 itor(4,4);
sf2d si,sj,sp;
std::vector<sf2d::Tf> ans(4);
while (itor)
{
si.make_sf(*this,itor[0]);
sj.make_sf(*this,itor[1]);
si.differentiate_y(mi);
//sj.differentiate_y(mi);
sp=si.times(mi,sj);
sp.integrate(mi);
ans[itor[1]]=sp.evaluate(mi);
mat(itor.cursor())=ans[itor[1]];
if ( itor[1]==3)  {
	MM_ERR(" lap2y "<<itor.string()<<" "<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<" "<<ans[3]) }

if ( itor[1]==3)  {
std::cout<<generate(itor,ans,2)<<CRLF;
}
++itor;
} // itor
mat.catalog();
} // gradnew2y




//////////////////////////////////////////////////////////




void grad2(const IdxTy dir )
{

typedef mjm_integrals Mi;
Mi mi;
// pairs of derivatives after integration by parts

LoopItor3 itor(4,4,4);
sf2d si,sj,sk,sp;
std::vector<sf2d::Tf> ans(4);
mjm_block_matrix<sf2d::Tf> mat(4,4,4);
while (itor)
{
si.make_sf(*this,itor[2]);
sj.make_sf(*this,itor[1]);
sk.make_sf(*this,itor[0]);
mat(itor.cursor()).of_float();
if (dir==0){ si.differentiate_x(mi); sj.differentiate_x(mi); }

if (dir==1){ si.differentiate_y(mi); sj.differentiate_y(mi); }

sp=si.times(mi,sj).times(mi,sk);
sp.integrate(mi);
ans[itor[2]]=sp.evaluate(mi);
mat(itor.cursor())=ans[itor[2]];
// not work for ints lol
ans[itor[2]].of_float();
if ( itor[2]==3)  {
std::cout<<generate(itor,ans,2)<<CRLF;
//	MM_ERR(" grad2 "<<itor.string()<<" "<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<" "<<ans[3]) }
}
++itor;
} // itor

mat.catalog();

} // grad2

// periphery terms after integration by parts of laplacian 
void parts_peri()
{

typedef mjm_integrals Mi;
Mi mi;
// pairs of derivatives after integration by parts

LoopItor3 itor(4,4,4);
sf2d si,sj,sp; // ,sk,sp;
std::vector<sf2d::Tf> ans(4);
// side, eqn and grad terms 
mjm_block_matrix<sf2d::Tf> mat(4,4,4);
while (itor)
{
const IdxTy side=itor[0];
si.make_sf(*this,itor[1]); // equation 
sj.make_sf(*this,itor[2]); // variable 
//sk.make_sf(*this,itor[0]);
// the one differentiated is the variable not the equation formed by multiplying by sf
const IdxTy dir=((side==0)||(side==2))?1:0;
if (dir==0){  sj.differentiate_x(mi); }
if (dir==1){  sj.differentiate_y(mi); }

sp=si.times(mi,sj);
// this needs to be evaluated at xy0 or 1...

if (dir==1) { sp.integrate_x(mi); }
if (dir==0) { sp.integrate_y(mi); }
double x=0;
double y=0;
double sign=1;
switch (side)
{
case 0: { x=1; y=0; sign=-1; break; }
case 1: { x=1; y=1; break; }
case 2: { x=1; y=1; break; }
case 3: { x=0; y=1; sign=-1; break; }

} // side
ans[itor[2]]=sp.evaluate(mi,x,y)*sign;


mat(itor.cursor())=ans[itor[2]];
mat(itor.cursor()).of_float();
// not work for ints lol
ans[itor[2]].of_float();
if ( itor[2]==3)  {
std::cout<<generate(itor,ans,2)<<CRLF;
//	MM_ERR(" grad2 "<<itor.string()<<" "<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<" "<<ans[3]) }
}
++itor;
} // itor
// FIXME wtf this does not generate decimal 
mat.catalog();

} // peri_part 



/////////////////////////////////////////////

template <class Tf, class Tv>
StrTy generate(const Tf & inicies, const Tv & vals, const IdxTy szz=0)
{
const IdxTy sz=(szz==0)?inicies.size():szz;
const IdxTy szv=vals.size();
std::stringstream ss;
for (IdxTy j=0; j<szv; ++j)
{
ss<<"mat(";
for (IdxTy i=0; i<sz; ++i) { ss<<inicies[i];  ss<<","; }
ss<<j<<")="<<vals[j]<<"; ";
} // j
return ss.str();
} // generate
/////////////////////////////////////////////////////////////




}; //shape_function_integrals
#ifdef  TEST_SHIFT__

int main(int argc,char **args)
{
typedef unsigned int IdxTy;
typedef double D ;
//typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_integrals Myt;
typedef std::vector<D>  P;
Myt mi;
D a=atof(args[1]);
//const bool dump_debug=false;
P x1,x2,x3,ux,x4;
ux.push_back(a);
ux.push_back(1);
x1.push_back(1.23);
x1.push_back(2.45);
x1.push_back(3.99);
x1.push_back(-4);
x1.push_back(5);
x1.push_back(-2);
x1.push_back(-.5);
mi.shift_polynomial(x2,x1,a);
mi.shift_polynomial(x3,x2,-a);
mi.composite_polynomial(x4,x1,ux);
const D x=1.234;
const D x1e= mi.evaluate_polynomial(x1, x);
const D x2e= mi.evaluate_polynomial(x2, x-a);
const D x3e= mi.evaluate_polynomial(x3, x);
const D x4e= mi.evaluate_polynomial(x4, x-a);
Ss ss;
const IdxTy sz=x2.size();
MM_MSG(MMPR(x4e)<<MMPR4(x3e,x2e,x1e,x))
for (IdxTy i=0; i<sz; ++i)
	{ ss<<MMPR2(x4[i],x3[i])<<MMPR4(x2[i],x1[i],a,i)<<CRLF; }
MM_MSG(":"<<CRLF<<ss.str())
return 0;
}

#endif

#ifdef  TEST_IBP__
//template<class Td >  typename Td::value_type x_n_exp_integral
//	(Td & dest, const IdxTy  & n, const typename Td::value_type & b)

int main(int argc,char **args)
{
typedef unsigned int IdxTy;
typedef double D ;
//typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_integrals Myt;
typedef std::vector<D>  P;
Myt mi;
D b=atof(args[1]);
IdxTy n=atoi(args[2]);
//const bool dump_debug=false;
P dest,desterfexp,desterferf;
D res=mi.x_n_exp_integral(dest, n, b);
mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
IdxTy sz=dest.size();
IdxTy szerf=desterferf.size();
IdxTy szexp=desterfexp.size();
MM_MSG(MMPR4(b,n,res,sz))
Ss ss;


//for (IdxTy i=0; i<sz; ++i) { ss<<MMPR2(dest[i],i)<<CRLF; }
for (IdxTy i=0; i<szerf; ++i) { 
const D ex=(i<szexp)?desterfexp[i]:0;
ss<<MMPR3(desterferf[i],ex,i)<<CRLF; }


MM_MSG(":"<<CRLF<<ss.str())
return 0;
}

#endif





#ifdef  TEST_GAUSS_DER__

int main(int argc,char **args)
{
typedef unsigned int IdxTy;
typedef double D ;
//typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_integrals Myt;
typedef std::vector<D>  P;
Myt mi;
const bool dump_debug=false;
// this should be somewhere else lol 
const bool make_table=false;
const bool one_point=false;
if ( make_table)
{
P x1,x2;
x1.push_back(1);
for (IdxTy i=0; i<100; ++i)
{
//template<class Td>  void next_gauss_der_polynomial(Td & dest, const Td & v1)
mi.next_gauss_der_polynomial(x2,x1);
Ss ss;
const IdxTy x2s=x2.size();
for (int j=x2s-1; j>0; --j)  if ( x2[j]!=0) { ss<<x2[j]; { ss<<"x^"<<j<<"+ ";} } 
if (x2s>0) {if (x2[0]!=0) {  ss<<x2[0];} }
MM_MSG(i<<": "<<ss.str())
x1=x2;
} // for
} // make_Table
if (one_point)
{
D dnorm=0; D derf=0;
const D xf=atof(args[1]);
const D xi=atof(args[2]);
const int  npts=atoi(args[3]);

//int doen( D & derf, D & dnorm, const D & xf, const D & xi, const IdxTy n)
for (int i=0; i<=npts; ++i)
{
const IdxTy nn=i;
//int  test=mi.doen(derf,dnorm, xf,xi,nn);
//mi.doen(derf,dnorm, xf,xi,nn,true);
mi.doen_split(derf,dnorm, xf,xi,nn,dump_debug);
//doen_split(  derf,  dnorm,  xf,  xi,  n,  debug=false)
D erff= erf(xf)-erf(xi);
D normf= exp(-xf*xf)-exp(-xi*xi);
D normden=normf+dnorm;
D releexp=(normden==0)?(normf-dnorm):((normf-dnorm)/(normden)*2);
D erfden=erff+derf;
D rele=(erfden==0)?(erff-derf):((erff-derf)/(erfden)*2);
MM_MSG(" diffs "<<MMPR(i)<<MMPR4(erff,derf,xf,xi)<<MMPR2(erff-derf,rele)<<MMPR4(normf,dnorm,normf-dnorm,releexp))

} // i 
} // one_point

if( true )
{
IdxTy i=0; 
for (D xi=-10; xi<10; xi+=.1) { 
for (D xf=xi+.1; xf<10; xf+=.1) { 

D dnorm=0; D derf=0;
mi.doen_split(derf,dnorm, xf,xi,0,dump_debug);
//doen_split(  derf,  dnorm,  xf,  xi,  n,  debug=false)
D erff= erf(xf)-erf(xi);
D normf= exp(-xf*xf)-exp(-xi*xi);
D normden=normf+dnorm;
D releexp=(normden==0)?(normf-dnorm):((normf-dnorm)/(normden)*2);
D erfden=erff+derf;
D rele=(erfden==0)?(erff-derf):((erff-derf)/(erfden)*2);
MM_MSG(" diffs "<<MMPR(i)<<MMPR4(erff,derf,xf,xi)<<MMPR2(erff-derf,rele)<<MMPR4(normf,dnorm,normf-dnorm,releexp))
++i;
} // xf 
} // xi 


} // true



return 0;
}

#endif // TEST_GAUSS_DER__

#ifdef  TEST_INTEGRAL_SF__

int main(int argc,char **args)
{
typedef  shape_function_integrals Ty;

Ty x;
x.laplacian();
x.lap2();
x.lap2y();
//x.grad();
//x.grad2(0);
if (false) x.grad2(1);
x.gradnewy();
//MM_MSG(" grad_avg "<<x.grad_avg(0,1,0,1))
MM_MSG("parts_peri")
if (false) x.parts_peri();
return 0; 
}
#endif // TEST_INTEGRAL_SF__

#ifdef  TEST_EXP_SF__

int main(int argc,char **args)
{
typedef mjm_integrals Myt;
typedef  Myt::exp_sf Ty;
typedef mjm_generic_itor<2> LoopItor2;
typedef unsigned int IdxTy;
typedef double D ;

Ty x;
LoopItor2 itor(5,5);
std::vector<D> values(4);
values[0]=1;
values[1]=2;
values[2]=3;
values[3]=4;
D avg=0;
D eavg=0;
for (IdxTy i=0; i<4; ++i){   avg+=.25*values[i]; eavg+=.25*exp(values[i]); } 
while(itor)
{
const IdxTy nx=itor[0]*30+5;
const IdxTy ny=itor[1]*30+5;
D answer=0;
for (IdxTy i=0; i<4; ++i) 
 answer+=x.int_trapezoid(values,i,nx,ny);
MM_MSG(" exp_sf_test "<<nx<<" "<<ny<<" "<<answer<<" avg "<<exp(avg)<<" "<<eavg<<" geomean "<<sqrt(eavg*exp(avg)))

++itor;
} // itor 

return 0; 
}
#endif // TEST_EXP_SF__




#ifdef  TEST_INTEGRAL_AS2__

int main(int argc,char **args)
{

typedef mjm_integrals Myt;
typedef unsigned int IdxTy;
typedef double D ;
//typedef mjm_block_matrix<double> MyBlock;
const D power=-1.0/3.0,beta=2.0/3.0; // 1.0/3.0;
Myt::as2  x(power,beta);
const D a=.01;
const D b=.01;


for (D rho=.1; rho<10; rho+=.1)
{
const D beff=b*rho;
// just moving the stupid singularity off the grid point
// would help a lot lol... 
D trap10=x.int_trapezoid(a,  beff, 10,10);
D trap100=x.int_trapezoid(a,  beff, 100,100);
D trap400=x.int_trapezoid(a,  beff, 400,400);
D pa1=x.int_parts(a,  beff, 1,1);
D pa10=x.int_parts(a,  beff, 10,10);
D pa20=x.int_parts(a,  beff, 300,300);
MM_MSG( " rho "<<rho<<" "<<trap10<<" "<<trap100<<" t4 "<<trap400<<" pa1 "<<pa1<<" "<<pa10<<" "<<pa20)

} // rho

return 0;
}

#endif

 
#ifdef  TEST_INTEGRAL__ 

int main(int argc,char **args)
{

typedef mjm_integrals Myt;
typedef unsigned int IdxTy;
typedef double D ;
//typedef mjm_block_matrix<double> MyBlock;
const D power=.9; // 1.0/3.0;
Myt::angle_scale  x(power,1);

for (D rho=1; rho<10; rho+=1)
{
// this is wrong 
//const D anl=0; // x.integral(rho)-x.integral(1);

//Real integrate_rect_bc(const Real & gamma,const Real & a,const Real & rhomax,const Real & rhomin)
const D n2=x.integrate_rect_bc(power,1.0,rho,1);
const D n3=x.integrate_rect_bc_2(power,1.0,rho,1);
const D npo=x.integrate_rect_bc_po(power,1.0,rho,1,250);

D num=0;
for (D i=0; i<=1000000; ++i) { const D dr=(rho-1)/1000000.0; const D ri=1+dr*i; num +=dr*x.integrand(ri); } 

MM_MSG(" rho "<<rho<<" anl "<<(n2)<<" "<<n3<<" po "<<npo<<" num " << num);

//const D npoias=x.integrate_rect_bc_po(power,1.0,rho,1,100000);
//for (IdxTy iter=2; iter<10000; iter=iter*2)
{
//const D npoi=x.integrate_rect_bc_po(power,1.0,rho,1,iter);

//MM_MSG( rho<<" "<<iter<<" "<<npoi<<" "<<(npoi-npoias))
}


}


return 0; 
} // main

#endif // test class

#ifdef  TEST_BURMANN__ 
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
Myt mi;
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
Dest dest,desthc;
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
mi.burmann_ref1_eq14(dest, funcerf, funcerfphi, funcerfphistar, z0, nu,  nmax, flags);
mi.burmann_ref1_eq14_erf( desthc, nmax);
compare_erf(dest);
while (desthc.size()<dest.size()) desthc.push_back(0); 
for (IdxTy i=0; i<dest.size(); ++i) { MM_MSG("bermann erf "<<MMPR3(i,dest[i],desthc[i])) }
MM_MSG(MMPR(mi.print_polynomial(dest)))
Dest deste,e1,destes;
e1.push_back(1);
e1.push_back(-1);
mi.composite_polynomial(deste,dest,e1);
mi.multiply_polynomial(destes,deste,2.0/sqrt(M_PI));
MM_MSG(MMPR(mi.print_polynomial(deste)))
MM_MSG(MMPR(mi.print_polynomial(destes)))
mi.dump_counter("final");
return 0;
} // main2
} // main 
#endif


#endif // guard 

