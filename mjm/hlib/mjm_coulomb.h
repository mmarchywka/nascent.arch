#ifndef MJM_COULOMB_H__
#define MJM_COULOMB_H__

 
#include <stdlib.h>
#include <gsl_sf.h>
#include <mjm_globals.h>
#include <mjm_templates.h>
#include <stdio.h>
#include <iostream>


/* contour integration to test some coulomb sto integral appraoches

*/
class mjm_coulomb{
typedef double D;
public:

static D integrand_1d(const D z, const D a, const D xprime, const D b)
{
D b2=b*b; 
const D a2=a*a;
const D zp=z-xprime;
return 1.0/((z*z+a2)*::sqrt(zp*zp+b2));

}

static D integrate_1d( const D a, const D xprime, const D b)
{
D i=0;
D z1=0;
D z2=0;
const D steps=.01;
const D scale=(a*a<b*b)?( ::fabs(a)) : (::fabs(b));
const D maxstep=scale*steps;
// shoulc pickup the tqails with expansion. 
// poles at +/-iA and xprime+/-ib
// of course this need not recompute squares on each call wtf.
 D v1=integrand_1d(z1,  a, xprime,  b);
 D v2=integrand_1d(z2,  a, xprime,  b);

i+=(v1+v2)*.5*(z2-z1);

return i;
}


static D contour_1d( const D a, const D xprime, const D b)
{
D i=0;
return i;
}
// atom a is at origin, b is on x axis at zb
#define LLOP(i,llim,lim,di) for ( D i=llim; i<lim; i+=di) { 
static D brute_force( const D alphaa, const D alphab, const D xb)
{
D i=0;
D extremes=alphaa;
if( (alphab)>extremes) extremes=alphab;
if( (xb)>extremes) extremes=xb;
extremes*=5;
const D mmin=-extremes;
const D mmax=extremes;
D dx=extremes*1e-1;
D dv=dx*dx*dx*dx*dx*dx;
LLOP(x,mmin,mmax,dx)
LLOP(y,mmin,mmax,dx)
LLOP(z,mmin,mmax,dx)
const D r=::sqrt(x*x+y*y+z*z);
const D dir=::exp(-2.0*alphaa*r)/(r*r);
LLOP(xp,mmin,mmax,dx)
const D xp2b=(xp-xb)*(xp-xb);
const D xp2=(xp)*(xp);
LLOP(yp,mmin,mmax,dx)
const D yp2=yp*yp;
LLOP(zp,mmin,mmax,dx)
const D rp=::sqrt(xp2+yp2+zp*zp);
const D rpb=::sqrt(xp2b+yp2+zp*zp);
D di=::exp(-2.0*alphaa*r);
i+=di*dv;

}}}}}}
return i;
}



}; // mjm_coulomb


#endif //guard


