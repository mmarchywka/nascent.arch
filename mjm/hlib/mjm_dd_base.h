#ifndef MJM_DD_BASE_H__
#define MJM_DD_BASE_H__
 
#include "mjm_text_data.h"
#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
//#include <mjm_templates.h>
// needed for file io crap
#include <mjm_csv_ops.h>
// finally added rule hits
//#include <mjm_sequence_hits.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

class dd_base_typedefs
{
public:
typedef unsigned int IdxTy;
typedef  int IntTy;
typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::stringstream SsTy;
};

namespace mjm_dd_base 
{
typedef dd_base_typedefs Tr;

typedef Tr::IdxTy IdxTy;
typedef  Tr::IntTy IntTy;
typedef Tr::ChTy ChTy;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;
enum { BAD=~0U}; 

typedef double Point;
typedef double Value;
// rather than a bunch of classes, use if/else for now
class Grid : public std::vector<Point>
{

public:
Grid() {Init();}
Grid(const Point & i, const Point & e, const IdxTy n) 
{Set(i,e,n); }

Grid(const Point & i) { Set(i,i,1); } 

void set(const Point & i, const Point & e, const IdxTy n) 
{Set(i,e,n); }

const bool & uniform() const { return m_uniform; }
const Value & min() const { return m_min; }
const Value & max() const { return m_max; }
const Value & min_dx() const { return m_min_dx; }
const Value & max_dx() const { return m_max_dx; }

void check() const { Check(); } 
private:

Value m_min, m_max, m_min_dx, m_max_dx;
bool m_uniform;
void Init()
{
m_min=0;
m_max=0;
m_min_dx=0;
m_max_dx=0;
m_uniform=true;
}
void Set(const Point & ii, const Point & e, const IdxTy n) 
{
clear();
MM_TRACE
if ( n==0) { Init(); return; } 
const Point del=(e-ii);
const double nf=n;
//MM_MSG(" n = "<<n<<" "<<nf)
for (IdxTy i=0; i<n; ++i) { push_back(ii+(i*del)/(nf));  }
m_min=(*this)[0];
m_max=(*this)[n-1];
// this is not quite righ, but ok for now
m_min_dx=del/nf; // ignores round of etc
m_max_dx=m_min_dx; // ignores round of etc
m_uniform=true;
}
void Check() const
{
//OsTy & os=std::cout();
const IdxTy sz=size();
Value dminx=1e30;
Value minx=0;
Value maxx=0;
Value last=0;
for (IdxTy i=0; i<sz; ++i) 
{
const Value & v=(*this)[i];
if ( v>maxx) maxx=v;
if ( v<minx) minx=v;
if ( last!=0) if ( (v-last)<dminx) dminx=v-last;
last=v;
}
MM_MSG( " grid check " << sz<<" "<<minx<<" "<<maxx<<" "<<dminx)

}

}; // Grid


class Values : public std::vector<Value>
{
typedef Values Myt;
public:
Values(): m_grid(0) {}
Values(const Grid & g): m_grid(&g) {resize(g.size());}

Myt resample(const Grid & gd) const { return Resample(gd); } 
Myt resample2x(const Grid & gd) const { return Resample2x(gd); } 
void clip(const Value & v, const IdxTy f) { Clip(v,f); } 
void set(const Value & v) { Set(v); } 
void add(const Value & v) { Add(v); } 
void div(const Value & v) { Div(v); } 
void mul(const Value & v) { Mul(v); } 
Value min() const  { return Min();  } 
Value max() const  { return Max();  } 
Value total() const  { return Total(); } 
Value sum_squares() const  { return Total2(); } 
//Value sum() const  { return Sum();  } 
const Grid & grid() const { return *m_grid;} 

// make operators and check grid at some point
Myt plus(const Values & that ) const { return Plus(that); } 
Myt plus(const Value & that ) const { return Plus(that); } 
Myt operator+(const Values & that) const { return Plus(that); } 
Myt operator+(const Value & that) const { return Plus(that); } 
Myt minus(const Values & that ) const { return Minus(that); } 
Myt operator-(const Values & that) const { return Minus(that); } 
Myt operator-(const Value & that) const { return Minus(that); } 
Myt times(const Values & that ) const { return Times(that); } 
Myt times(const Value & that ) const { return Times(that); } 
Myt operator*(const Values & that ) const { return Times(that); } 
Myt operator*(const Grid & that ) const { return Times(that); } 
Myt operator/(const Values & that ) const { return Div(that); } 
Myt operator*(const Value & that ) const { return Times(that); } 
Myt ln( ) const { return Ln(); } 
Myt ln(const Value lim ,const Value limv)  const { return Ln(lim,limv); }
void mix(const Values & that, const Value & fold, const IdxTy first, const IdxTy last)
{
Mix(that,fold,first,last);
}
void axby(Values & d, const Values & that,const Value & a, const Value & b ) 
{ return Axby(d,that,a,b); } 
Myt axby(const Values & that,const Value & a, const Value & b ) 
{ return Axby(that,a,b); } 
void boundaries(const Value & i, const Value & f)
{ Boundaries(i,f); }
void boundaries(const Value & i)
{ Boundaries(i); }


private:
const Grid * m_grid;
const 
void CheckGrid(const Values & that ) const
{ if ( m_grid!=that.m_grid) {MM_MSG(" grids not same object fwiw ") } }

void Axby(Values & d, const Values & that,const Value & a, const Value & b )  const
{
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) d[i]=a*(*this)[i]+b*that[i]; } 


Myt Axby(const Values & that,const Value & a, const Value & b )  const
{ Myt d(*m_grid); CheckGrid(that);
Axby(d,that,a,b);
return d;
}

Myt Ln( )  const
{ Myt d(*m_grid); 
const IdxTy sz=size(); 
for ( IdxTy i=0; i<sz; ++i) 
{
const Value x=(*this)[i];
d[i]=::log(x);  
}

return d;
}
Myt Ln(const Value lim ,const Value limv)  const
{ Myt d(*m_grid); 
const IdxTy sz=size(); 
for ( IdxTy i=0; i<sz; ++i) 
{
const Value x=(*this)[i];
if ( x>lim)d[i]=::log(x);  
else d[i]=limv;
}

return d;
}







void Mix(const Values & that, const Value & fold, const IdxTy first, const IdxTy last)
{
Myt & d=(*this);
const Myt& r=(*this);
const Value fnew=1.0-fold;
//const IdxTy sz=size(); 
for ( IdxTy i=first; i<last; ++i) d[i]=fold*r[i]+fnew*that[i];  



}
// apply linear xform to fit beginning and end
void Boundaries(const Value & ii, const Value & f)
{
Myt & x=*this;
const IdxTy sz=size(); 
const Value a=x[0]-ii;
const Value b=(f-x[sz-1]+a)/(sz-1);
for ( IdxTy i=0; i<sz; ++i) x[i]=x[i]-a+b*i;  

}

void Boundaries(const Value & ii)
{
Myt & x=*this;
const IdxTy sz=size(); 
// bombs if size is zero...
const Value a=x[0]-ii;
for ( IdxTy i=0; i<sz; ++i) x[i]=x[i]-a;  

}





Myt Plus(const Values & that ) const
{ Myt d(*m_grid); CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) {d[i]=(*this)[i]+that[i]; } 
return d;
} 
Myt Plus(const Value &  that ) const
{ Myt d(*m_grid); // CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]+that; } 
return d;
} 


Myt Minus(const Values & that ) const
{ Myt d(*m_grid); CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]-that[i]; } 
return d;
} 
Myt Minus(const Value & that ) const
{ Myt d(*m_grid); // CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]-that; } 
return d;
} 

Myt Times(const Values & that ) const
{ Myt d(*m_grid); CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]*that[i]; } 
return d;
} 

Myt Times(const Grid & that ) const
{ Myt d(*m_grid); //CheckGrid(that);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]*that[i]; } 
return d;
} 



Myt Times(const Value & that ) const
{ Myt d(*m_grid); 
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) { d[i]=(*this)[i]*that; } 
return d;
} 

Value Total() const  { 
Value t=0;
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) t+=(*this)[i];
return t;
 } 
Value Total2() const  { 
Value t=0;
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) t+=(*this)[i]*(*this)[i];
return t;
 } 



void Add(const Value & v)  { 
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) (*this)[i]+=v; } 

void Set(const Value & v) { 
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) (*this)[i]=v; } 

void Clip(const Value & v, const IdxTy f) { 
const IdxTy sz=size(); 
switch ( f) 
{

case 0:
{  for ( IdxTy i=0; i<sz; ++i) { Value & x= (*this)[i]; if ( x<v) x=v; }
break; 
}

case 1:
{  for ( IdxTy i=0; i<sz; ++i) { Value & x= (*this)[i]; if ( x>v) x=v; }
break; 
}




} // switch f

} 

void Mul(const Value & v) { 
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) (*this)[i]*=v; } 


Myt Div(const Values & v) const { 
Myt d(*m_grid);
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) d[i]=(*this)[i]/v[i]; 
return d;
} 

void Div(const Value & v) { 
const IdxTy sz=size(); for ( IdxTy i=0; i<sz; ++i) (*this)[i]/=v; } 

Value Min() const { 
const IdxTy sz=size(); 
if ( sz==0) return Value();
Value x=(*this)[0];
for ( IdxTy i=0; i<sz; ++i) if ((*this)[i]<x) x=(*this)[i];
return x;  
} 

Value Max() const { 
const IdxTy sz=size(); 
if ( sz==0) return Value();
Value x=(*this)[0];
for ( IdxTy i=0; i<sz; ++i) if ((*this)[i]>x) x=(*this)[i];
return x;  
} 


// assume or verify simple 2x dilation
Myt Resample2x(const Grid & gd) const 
{
Myt d(gd);
const Values & s=*this;
if ( m_grid==0) MM_MSG(" Resample needs a grid will bomb now ")
const Grid & gs=*m_grid;
if ( !(gs.uniform()&&gd.uniform())) { MM_MSG(" resample for uniforms only") } 
const IdxTy szd=gd.size();
const IdxTy szs=gs.size();
if ( (2*szs)!=szd){  MM_MSG(" resample 2x size "<<szs<<" vs "<<szd) } 

IdxTy j=0;
for ( IdxTy i=0; i<szs; ++i) { d[j]=s[i]; d[j+1]=s[i]; j+=2;} 

return d;
}


Myt Resample(const Grid & gd) const 
{ 
Myt d(gd);
const Values & s=*this;
if ( m_grid==0) MM_MSG(" Resample needs a grid will bomb now ")
const Grid & gs=*m_grid;
if ( !(gs.uniform()&&gd.uniform())) MM_MSG(" resample for uniforms only")
//const Value  sdx=gs.min_dx();
//const Value  ddx=gd.min_dx();
const IdxTy szd=gd.size();
const IdxTy szs=gs.size();
const IdxTy szs1=szs-1;
 IdxTy j=0;
//const Value&  dv=gd[0];
for ( IdxTy i=0; i<szd; ++i)
{
while ( j<szs1)
{
const Value d1=gd[i]-gs[j];
const Value d2=gd[i]-gs[j+1];
// the both sources are greater, the lower oneis closer
if ( (d1<0)&&(d2<0)) { break ;}
// both sources less, the ++ is nearer but ma not be nearest
if ( (d1>0)&&(d2>0)) {++j; continue;}
// straddles, return closest
if ( (d1>0)&&(d2<0)) {if ( (d1+d2)>0)++j; break;}
if ( d1==0) break;
if ( d2==0) { ++j; break; }
MM_MSG(" fell through resample lol "<<i<<" "<<j<<" "<<gd[i]<<" "<<gd[j])
++j;
}


if ( j>=szs) j=szs1;
d[i]=s[j];
//MM_MSG( " resamplig "<<i<<" to "<< j)
} // i 



return d; 
} 


}; // values
class GeneralParameters
{

public:
virtual Value q() const { return 1.6e-19; } 
virtual Value mu() const { return 10; } 
virtual Value vt() const { return .0259; } 
virtual Value time_step() const { return 1e-10; } 
//virtual Value time_step() const { return 1e-7; } 
// this is now silicon lol 
virtual Value epsilon() const { return 12*8.854e-14; } 
virtual Value epsilon_vacu() const { return 8.854e-14; } 
virtual Value aRichardson() const { return 1.2017e2; } 
virtual Value dwFac() const { return q()/epsilon()/4.0/M_PI; } 



}; // parameters
typedef GeneralParameters GP;
typedef Value value_type;
typedef Values value_array;

// right now these are only reverse-complement hits
// well if there was stuff in there,,,
class dummy
{

public:
// drift diffusion equaiton
static void J(Values & j, const Grid & g, const Values & n, const Values & gradn, const Values & e, const GP & gp)
{
const Value k=gp.q()*gp.mu();
const Value vt=gp.vt();
const IdxTy sz=g.size();
for  ( IdxTy i=0; i<sz; ++i) { j[i]=k*(n[i]*e[i]+0*vt*gradn[i]); }

}
static Value thermionic(const Value & Vt, const Value & sf, const Value & e,
const GP & gp)
{
const Value aRich=gp.aRichardson();
const Value f=gp.dwFac();
 Value dw=(f*e);
const Value T=Vt/.0259*300.0;
dw=(dw<0)?(-::sqrt(-dw)): ::sqrt(dw);
return aRich*.5*T*T*::exp(-(sf-dw)/Vt); 
}
static void step_rho(Values & rhoi, const Grid & g, const Values & gradji, const Values & Gi, const Values & Ri, const GP & gp)
{
const Value dt=gp.time_step();
// thse are off by afctor of q
const IdxTy szg=Gi.size();
const IdxTy szr=Ri.size();
const IdxTy sz=g.size();
if ( (szg==sz)&&(szr==sz))
{for  ( IdxTy i=0; i<sz; ++i) { rhoi[i]-=(gradji[i]+Gi[i]-Ri[i])*dt; }
}
else{ for  ( IdxTy i=0; i<sz; ++i) { rhoi[i]-=(gradji[i])*dt; } }

}
// filter
static void F(Values & d, const Grid & g, const Values & n)
{
const IdxTy sz=g.size();
if (sz==0) return; 
const Value ap1=.01;
const Value an1=ap1;
const Value a0=1.0-ap1-an1;
const IdxTy sze=sz-1;
for  ( IdxTy i=1; i<sze; ++i) { d[i]=(ap1*n[i+1]+an1*n[i-1]+a0*n[i]);  }
d[0]=n[0];
d[sze]=n[sze];

}

// finite diff Derivative for now 
static void D(Values & d, const Grid & g, const Values & n)
{

const IdxTy sz=g.size();
if (sz==0) return; 
const IdxTy sze=sz-1;
for  ( IdxTy i=1; i<sze; ++i) { d[i]=(n[i+1]-n[i-1])/(g[i+1]-g[i-1]);  }
if ( sz>1) d[0]=(n[1]-n[0])/(g[1]-g[0]);
if ( sz>1) d[sze]=(n[sze]-n[sze-1])/(g[sze]-g[sze-1]);
else d[0]=0;

}

// trapezoid integrate
// this was off a factor of 2, may have botched prior usages.. 
static void Iold(Values & d, const Grid & g, const Values & n)
{
const IdxTy sz=g.size();
if (sz==0) return; 
const IdxTy sze=sz-1;
d[0]=0;  
if ( sze>1) d[0]=.25*((n[1]+n[0])*(g[1]-g[0]));  
if ( sze>1) d[sze]=-.25*((n[sze-1]+n[sze])*(g[sze-1]-g[sze]));  
for  ( IdxTy i=1; i<sze; ++i) 
{ 
// this part is wrong, it shifts the limit by 1/2 a thingy
d[i]=.25*((n[i+1]+n[i])*(g[i+1]-g[i]))+d[i-1];  
d[i]-=.25*((n[i-1]+n[i])*(g[i-1]-g[i]));  
}

if ( sze>1) d[sze]+=d[sze-1];
} // I 

static void I(Values & d, const Grid & g, const Values & n)
{
const IdxTy sz=g.size();
if (sz==0) return; 
const IdxTy sze=sz-0;
d[0]=0;  
//if ( sze>1) d[sze]=-.25*((n[sze-1]+n[sze])*(g[sze-1]-g[sze]));  
for  ( IdxTy i=1; i<sze; ++i) 
{ 
// this part is wrong, it shifts the limit by 1/2 a thingy
d[i]=-.5*((n[i-1]+n[i])*(g[i-1]-g[i]))+d[i-1];  
}

//if ( sze>1) d[sze]+=d[sze-1];
}


static Value I( const Grid & g, const Values & n)
{

const IdxTy sz=g.size();
if (sz==0) return 0 ; 
Value d=0;
const IdxTy sze=sz-1;
if ( sze>1) d=.5*((n[1]+n[0])*(g[1]-g[0]));  
if ( sze>1) d+=-.5*((n[sze-1]+n[sze])*(g[sze-1]-g[sze]));  
for  ( IdxTy i=1; i<sze; ++i) 
{ 
d+=.5*((n[i+1]+n[i])*(g[i+1]-g[i]));  
d-=.5*((n[i-1]+n[i])*(g[i-1]-g[i]));  
}
//if ( sze>1) d[sze]+=d[sze-1];
return d;
} // I 




static Value image(const Values & rho,const Value & p1, const Value & p2)
{
Value x=0;
const Value dx=p2-p1;
const IdxTy sz=rho.size();
const Grid & g=rho.grid();
const Value cgrid=g[1]-g[0]; // fix this for all grids, not fixed larger
for ( IdxTy i=0; i<sz; ++i)
{
x+=rho[i]*(p2-g[i])/dx*cgrid;

}

return x;
}

/////////////////////////////////////////////////////
//step(rho,too_left,too_right,e,gp,dt);
static void step(Values & rho, Value & too_left, Value & too_right,
 const Values & e, const GP & gp, const Value dt)
{
const Grid & g=rho.grid();
if ( !g.uniform()) MM_MSG(" needs uniform grid")
Values rho2(g);
//  move charge based on veloctiy distro, need a model for this
const IdxTy sz=g.size();
const Value dnx=dt/g.min_dx();
//const Value left=g[0];
//const Value right=g[sz-1];
too_left=0;
too_right=0;
for ( IdxTy i=0; i<sz; ++i)
{
rho2[i]+=rho[i];
// negative field moves eletrons to right
// rho is negative, negative field moves negative charge right
Grid v_grid(-gp.mu()*e[i]);
Values distro(v_grid);
distro[0]=1;
const IdxTy vsz=v_grid.size();
for ( IdxTy j=0; j<vsz; ++j)
{
// the grid points are velocities samples, the velocityes are fraction
Value  dn=dnx*v_grid[j];
Value drho=distro[j]*rho[i];
IntTy dni=0;
//MM_MSG(" stepping "<<i<<" "<<dn<<" "<<drho)
if ( (dn)<(-1.0*i-.5))  {
MM_MSG(" LEFT SHOT FUDD "<<dn<<" "<<i<<" ASSSSSSSS "<<(-i-.5))

too_left+=drho; rho2[i]-=drho; } 
else if ( dn>(sz+.5-i)) {
MM_MSG(" RIGHT SHOT FUDD "<<dn<<" "<<i)

 too_right+=drho; rho2[i]-=drho; }
else
{
// just make if positive kind of dumb.. 
dni=IdxTy(dn+sz+.5)-sz;
Value ty=dn-dni; // try to split it up a bi t
if (( dni<0)||(dni>=IntTy(sz)))  MM_MSG(" bad dn "<<dni<<" for "<<dn<<" and "<<sz)
if ( ty>0){ 
//MM_MSG(" ASS FUCVK "<<ty<<" "<<i<<" "<<dn<<" "<<dni)
if (( i+dni+1)<sz) { rho2[i+dni]+=(1.0-ty)*drho; rho2[i+dni+1]+=ty*drho;}
else rho2[i+dni]+=drho;

 }
else{ 
// ASS FUDDING SIGNED FUDDIG INT AGAIN FUDD 
//MM_MSG(" SHOT FUCVK "<<ty<<" "<<i<<" "<<dn<<" "<<dni)
// slightly too left
 rho2[i+dni]+=(1.0+ty)*drho;
if ( IntTy(i+dni-1)>=0) {  rho2[i+dni-1]-=ty*drho;} 
else too_left-=ty*drho;
//else rho2[i+dni]+=drho;

 }
// rho2[i+dn]+=drho; 
rho2[i]-=drho; 
} // else 
} // j 
} // i 
rho=rho2;
} 

////////////////////////////////////////////////////



class scl_state
{
typedef scl_state Myt;

public:
scl_state() {}
scl_state(const Grid & gg, const GP & gpp)
: g(gg),gp(gpp),j(g),n(g),rho(g),e(g),gradj(g),gradn(g),temp(g)
{

}
IdxTy npts;
 Value  dist;
Value Eext;

Grid g;
GP gp;
Values j, n,rho;
Values e,gradj,gradn;
Values temp;

StrTy dump1() const
{ SsTy ss;
ss<<" j0="<<j[0]<<" j "<<j.min()<<" "<<j.max()<<" e "<<e.min()
<<" "<<e.max()<<" n "<<n.min()<<" "<<n.max()<<" ntot="<<I(g,n);
return ss.str();
}
// bad idea for long things... 
StrTy dump2(const IdxTy i ) const
{ SsTy ss;
const IdxTy sz=g.size();
OsTy & os=ss;
for ( IdxTy k=0; k<sz; ++k) os<<"state "<<i<<" "<<k<<" "<<g[k]<<" "<<j[k]<<" "<<gradj[k]<<" "<<n[k]<<" "<<gradn[k]<<" "<<e[k]<<CRLF; 


return ss.str();
}


Myt resample ( const Grid & gd) // , const Myt & s)
{
const Myt & s=*this; 
Myt d(gd,s.gp);
// not for the grid, for the run
d.npts=s.npts;
d.dist=s.dist;
d.Eext=s.Eext;

d.gp=s.gp;
d.j=s.j.resample(gd);
d.n=s.n.resample(gd);
d.rho=s.rho.resample(gd);
d.e=s.e.resample(gd);
d.gradj=s.gradj.resample(gd);
d.gradn=s.gradn.resample(gd);
//d.temp=s.temp.resample(gd);

return d; 
}



} ; // scl_state

// a thermionic cathode that injects essentially fixed
// current into a mobility controlled medium create
// IV curves with varied Eapplied and length.
// In low bias, most of the current is reflected
// as E*gradE groes through zero as e==sqrt(x). 
static void test_scl( scl_state & ss, const IdxTy iter)
{
// makeing shorthand to reference would make other code shorter
Value t=0;
IdxTy i=0;
//const bool debug=false;
//n[10]=.01;
while (true)
{
Value Ezed=-image(ss.rho,ss.g[0],ss.g[ss.g.size()-1])/ss.gp.epsilon();
I(ss.e,ss.g,ss.rho); ss.e.div(ss.gp.epsilon());
ss.e.add(ss.Eext+Ezed);
J(ss.j,ss.g,ss.n,ss.gradn,ss.e,ss.gp);
//if ( i==0) j[0]=1e-6;
const Value icathode=thermionic(.0259*4 ,2,ss.e[0],ss.gp);
ss.j[0]-=icathode;
D(ss.gradj,ss.g,ss.j);
step_rho(ss.rho,ss.g,ss.gradj,Values(),Values(),ss.gp);
ss.rho.clip(0,1);
ss.n=ss.rho;
ss.n.div(-ss.gp.q());
D(ss.gradn,ss.g,ss.n);
if (( i%10000) ==0 ){
 MM_MSG(" scl "<<i<<" icath="<<icathode<< ss.dump1()<<" Ezed="<<Ezed)
}

const bool dump_state=!true; 
if ( dump_state) { std::cout<<ss.dump2(i); } 

++i;
if ( iter!=0) if ( i>iter){ ss.j[0]+=icathode;  return; } 
t+=ss.gp.time_step();

} // trye

} // scl_state 

// this looks like it matches the drift only simulation somewhat
// but not exactly. This is based on uniform e*grad_e product tokeep
// current uniform giving a sqrt(x) field that changes direction
// of current when ti goes through zero near the cathode and
// accounts for reflection ofexcess. . This is somewhat
// sensitive to cathode model but details of contact are
// vague too. The "turning point" is probably the first point
// in the mesh but not exactly due to derivative at edges. It
// may be interesting to make this model more realistic.
// In any case, adding diffusion will mess this up but helpful
// to understand if numerical code is close.
static void test_scl3( const IdxTy npts, const Value & dist, const Value
& Eext)
{

const Value icath=.58;
GP gp;
const Value step=::sqrt(::sqrt(::sqrt(10.0)));
//const Value cx=1.0/(gp.q()*gp.mu()*gp.epsilon());
const Value A=2.0/3.0*::sqrt(2.0/(gp.mu()*gp.epsilon()));
//const Value ic=icath*cx;
const Value l=1.0;
const Value dx=1e-3;
const Value  dx32=::sqrt(dx)*dx;
const Value  dxl32=::sqrt(l-dx)*(l-dx);
//Value jmin=0;
//Value jmax=.38;
//Value j=(jmax+jmin)*.5;;
Value jdevice=1e-12;
while ( true)
{
const Value   jdelta=icath-jdevice;
const Value rhs= A*(::sqrt(jdelta)*dx32-::sqrt(jdevice)*dxl32);
MM_MSG(" ivpoint icath="<<icath<<" jdev="<<jdevice<<" V="<<rhs<<" Eapp="<<(rhs/l))
jdevice*=step;
}



} //3 

static void test_scl2( const IdxTy npts, const Value & dist, const Value
& Eext)
{

GP gp;
Grid g(0,dist,npts);
scl_state ss(g,gp);

ss.dist=dist;
ss.npts=npts;
ss.Eext=Eext;


const bool a_once=false;
// this looks like it restarts from a fudd 
if ( a_once)
{
test_scl(ss,1000000);
Grid g2(0,dist,10*npts);
scl_state s2=ss.resample(g2);
MM_MSG(" resample begins "<<s2.dump1())
test_scl(s2,00);
}
else
{
test_scl(ss,2000000);
Grid g2(0,dist,2*npts);
scl_state s2=ss.resample(g2);
test_scl(s2,100000);
Grid g3(0,dist,4*npts);
scl_state s3=s2.resample(g3);
test_scl(s3,100000);
Grid g4(0,dist,8*npts);
scl_state s4=s3.resample(g4);
test_scl(s3,10000);
Grid g5(0,dist,16*npts);
scl_state s5=s4.resample(g5);
test_scl(s4,0);




}




} //scl_2


// mobility 10, 300k etc;  should be cgs through out
static void test_scl( const IdxTy npts, const Value & dist, const Value
& Eext_, const Value tol )
{
Grid g(0,dist,npts);
g.check();
GP gp;
Values j(g), n(g),rho(g);
Values e(g),gradj(g),gradn(g);
Values temp(g);
Value Eext=Eext_;
while ( true)
{
Value t=0;
IdxTy i=0;
const bool debug=false;
//n[10]=.01;
while (true)
{
// calculate the image charges on teh cathode to get
// the boundary conditiosn right
Value Ezed=-image(rho,g[0],g[g.size()-1])/gp.epsilon();
// n is really rho...
I(e,g,rho); e.div(gp.epsilon());
e.add(Eext+Ezed);
J(j,g,n,gradn,e,gp);
//if ( i==0) j[0]=1e-6;
const Value icathode=thermionic(.0259*4 ,2,e[0],gp);
j[0]-=icathode;
if ( debug) MM_MSG(" j[0]="<<j[0]<<" field at cathode "<<e[0]<<" ezed="<<Ezed)
//if ( i==0) j.set(j[0]);
//j[1]=j[0];
//j[g.size()-1]=j[0];
D(gradj,g,j);
//F(temp,g,gradj); gradj=temp;
if ( debug) MM_MSG(" gradj "<<gradj.min()<<" "<<gradj.max())

// ok, rho is NOT n lol. 
const bool step_continuity = !false; 
if (step_continuity)
{
// this is needed to pickup the incoming current and
// convert to accumulated charge. 
step_rho(rho,g,gradj,Values(),Values(),gp);
}
else
{
// this should besame as above but more exact once
// incoming is dealt with. This conserves charge, except for
// outgoing. 
const Value dt=gp.time_step();
// injected needs to be handled better but insteady state ok I guess
rho[0]-=icathode*dt; // arrgh, so this just comes in and stops....
Value too_left=0; Value too_right=0;
// these too_xxx are terminal currents. 
step(rho,too_left,too_right,e,gp,dt);

}


//rho[0]=0; rho[g.size()-1]=0;
//rho[0]+=j[0]*gp.time_step()/(g[1]-g[0]);

rho.clip(0,1);
//F(temp,g,rho);rho=temp;

n=rho;
n.div(-gp.q());
D(gradn,g,n);


if ( debug) MM_MSG(" gradn "<<gradn.min()<<" "<<gradn.max())


if (( i%10000) ==0 ){
 MM_MSG(" scl "<<i<<" icath="<<icathode<<" j0="<<j[0]<<" j "<<j.min()<<" "<<j.max()<<" e "<<e.min()
<<" "<<e.max()<<" n "<<n.min()<<" "<<n.max()<<" ntot="<<I(g,n)<<" Ezed="<<Ezed)
//const bool dump_state=!true; 
const Value jmin=j.min();
const Value jmax=j.max();
 bool jclose=(((::fabs(jmax-jmin))/(::fabs(jmax)+::fabs(jmin)))<tol);
//MM_MSG("jclose is "<<jclose)
//const bool dump_state=(i>1000000)&&(jclose); 
const bool dump_state=(i>100000)&&(jclose||(tol==0))&&false; 
if ( dump_state)
{
MM_MSG("dumping")
const IdxTy sz=g.size();
OsTy & os=std::cout;
for ( IdxTy k=0; k<sz; ++k) {os<<"state "<<i<<" "<<k<<" "<<g[k]<<" "<<j[k]<<" "<<gradj[k]<<" "<<n[k]<<" "<<gradn[k]<<" "<<e[k]<<CRLF; 
}
MM_MSG(" converged ext="<<Eext<<" jmin="<<j.min()<<" jmax="<<j.max())
if ( tol!=0) break; // return;
} // dump_state
} // dump


++i;
t+=gp.time_step();
} // true
// change Eext if needed
if ( 1==1 ) break;
Eext=Eext/::sqrt(2.0);

} // true 2....


} // test+Sc;


// in areas of G/R. keep J fixed but allow Jn and Jp to vary,

// for updating n, change sign of kt
// e is ngative the field
static void fixed_j_total_x(Values & nnew, const Values & v,const Values & n,const Values & p,const Values & rhodiveps, const Values & e , const Value & ni2, const Value & G, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q)
{
const Value dx2=dx*dx;
//const Value dx2i=1.0/dx2;
const Value kT2=kT*2.0;
//const Value c1=dx2/q/mu;
const Value c1=-G*dx2/mu;
const Value ni=::sqrt(ni2);
for ( IdxTy i=2; i<(sz-1); ++i)
{
//const Value poc=kT2+dx2*rhodiveps[i];
// move the npg n part to rhs. 
const Value poc=-c1*p[i]-kT2+(2.0*v[i]-v[i-1]-v[i+1]);
const Value gi=c1*(-ni2);
//const Value gi=c1*(n[i]*p[i]-ni2);
// make e as grad v note sign change 
//pnew[i]=1.0/poc*(c1*gi-(-.5*dx*e[i]-kT)*(p[i+1]-p[i-1]));
const Value dv=v[i+1]-v[i-1];
const Value dn=n[i+1]-n[i-1];
const Value sn=n[i+1]+n[i-1];
if ( false)
{

nnew[i]=1.0/poc*(gi+.25*(dn)*(dv)-kT*(sn));
}
else{

const Value dvp=v[i+1]-v[i];
const Value dvn=v[i]-v[i-1];
const Value dnp=n[i+1]-n[i];
const Value dnn=n[i]-n[i-1];
if ( false)
// this seems to explode, 
nnew[i]=1.0/poc*(gi+.5*(dnp)*(dvp)+.5*dnn*dvn-kT*(sn));
else
{
// this seems similar to the original, 
// AFAICT this prevents nan but still not right 
const Value srhkluge=1.0; // (p[i]+n[i]+2*ni)/ni;
const Value poc2=-c1*p[i]/srhkluge-kT2+(v[i]-.5*v[i-1]-.5*v[i+1]);
nnew[i]=1.0/poc2*(gi/srhkluge-.5*(v[i]*n[i+1]+n[i-1]*v[i]-n[i+1]*v[i+1]-n[i-1]*v[i-1])-kT*(sn));
const bool del_rho_del_t=!true;
if ( del_rho_del_t)
{
const Value n_=n[i];
const Value v_=v[i];
const Value delndelE=.5*((n[i+1]-n_)*(v_-v[i+1])+(n_-n[i-1])*(v[i-1]-v_));
const Value ndelE=-n_*(v[i+1]+v[i-1]-2*v_);
const Value lapn=kT*(n[i+1]+n[i-1]-2*n_);
const Value delrho=mu*(delndelE+ndelE+lapn)/dx2-.01*G*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni)*ni;
if ( kT<0) 
nnew[i]=n[i]-delrho*1e-18;
else nnew[i]=n[i]+delrho*1e-18;
//else nnew[i]=n[i]-delrho*1e-19;
//else nnew[i]-=delrho*1e-14;

} // del rho

// tack in srh lol 
//nnew[i]=1.0/poc2*(gi/(p[i]+n[i]+2*ni)-.5*(v[i]*n[i+1]+n[i-1]*v[i]-n[i+1]*v[i+1]-n[i-1]*v[i-1])-kT*(sn));
// discard higher order terms. 
// worse
//nnew[i]=1.0/poc2*(gi-.5*(v[i]*n[i+1]+n[i-1]*v[i])-kT*(sn));


} // solve for n again... 

} 


if ( nnew[i]<0)
{
MM_MSG(" clipping charge "<<i<<" "<<nnew[i])
 nnew[i]=0;
}


} //i 
} 

static void fixed_j_total(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p,const Values & rhodiveps, const Values & e , const Value & ni2, const Value & G, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q)
{
// this should use mu p and then mu n 
fixed_j_total_x(nnew,v,n,p,rhodiveps,e,ni2,G,mu,dx,sz,-kT,q);
fixed_j_total_x(pnew,v,p,n,rhodiveps,e,ni2,-G,mu,dx,sz,kT,q);

}



// ================ 55555555555555555555555555555555 ==================
// j is the reduced current, J div q and mu. 
static Value nn( const Value & n0, const Value & dv,  const Value & kT, const Value & ni2, const Value & j)
{

if ( j==0) return n0*::exp(dv/kT);
Value nf=n0;
//const Value lhs= ::pow(n0,kT)*::pow(n0*n0+ni2,-j)*::exp(dv);

//const Value rhs=::pow(nf,kT)*::pow(nf*nf+ni2,-j);
MM_MSG(" did not eiomplement this yet as the behavrir seems wrong")
return nf;

}
static void sum_and_diff5(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p, const Value & ni2, const Value & G_, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q, const Value & J,const Value Neff, const Values & Nteff)
{
//const Value ni=::sqrt(ni2);
const Value j=J/q/mu;
//const Value kti=1.0/kT;
//const Value G2=G_/mu;
//const Grid & g= n.grid();
static IdxTy iter=0;
//const Value dxx2=dx*dx;
//Values xnew(g),ynew(g),G(g);
//for ( IdxTy i=0; i<sz; ++i) { G[i]=G2*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni);}  


for ( IdxTy i=0; i<sz; ++i) {  
nnew[i]=nn(n[0],v[i]-v[0],kT,ni2,j);
pnew[i]=nn(p[sz-1],v[i]-v[sz-1],-kT,ni2,j);
if ( pnew[i]>nnew[i]) 
{
nnew[i]=ni2/pnew[i]; // nn(ni2/n[0],v[i]-v[0],-kT,ni2,j);
}
else{
pnew[i]=ni2/nnew[i]; // nn(ni2/n[0],v[i]-v[0],-kT,ni2,j);
}

MM_MSG(iter<<" "<<i<<" nnew= "<<nnew[i]<<" pnew= "<<pnew[i]<<" dv= "<<(v[i]-v[0])<<" n= "<<n[i]<<" p= "<<p[i])
} //i 
 ++iter;

pnew[0]=p[0];
pnew[1]=p[1];
pnew[sz-1]=p[sz-1];
pnew[sz-2]=p[sz-2];

nnew[0]=n[0];
nnew[1]=n[1];
nnew[sz-1]=n[sz-1];
nnew[sz-2]=n[sz-2];


} 

// =============== 44444444444444444444444444444444 ===================
static void sum_and_diff4(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p, const Value & ni2, const Value & G_, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q, const Value & J,const Value Neff)
{
const Value ni=::sqrt(ni2);
const Value j=J/q/mu;
const Value kti=1.0/kT;
const Value G2=G_/mu;
const Grid & g= n.grid();
//const Value dxx2=dx*dx;
Values xnew(g),ynew(g),G(g);
//Values & x=(n);
//Values & y=(p);
// absorb mu here as n and p aer same 
for ( IdxTy i=0; i<sz; ++i) { G[i]=G2*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni);} //i 
//for ( IdxTy i=0; i<sz; ++i) { G[i]=G2*(p[i]*n[i]-ni2)/1e13;} //i 
Values e(g),n1(g),n2(g);
D(e,g,v);
D(n1,g,n);
D(n2,g,n1);
Values en=e*n;
D(n1,g,en);
Values dn=(G-n1+n2*kti);
D(n1,g,p);
D(n2,g,n1);
en=e*p;
D(n1,g,en);
Values dp=(G+n1+n2*kti);
Value dt=mu*1e-19;
if ( true)
{
D(n1,g,n-p);
Values jold=e*(n+p)+n1*kT;
// not clipped for less than zed yet
D(n1,g,dn-dp);
Values jinc=e*(dn+dp)+n1*kT;
// rms
Value s1=0; Value s2=0; Value s3=0;
// remove edge cap
for ( IdxTy i=2; i<(sz-2); ++i) 
{
s1+=jinc[i]*jinc[i];
s2+=jinc[i]*jold[i];
s3+=jinc[i]*j;

}
if ( s1!=0) dt=-.5*(s2-s3)/s1;
//const Value jto=jinc.total();
//c*jinc.total()+jold.total()=j*sz;
//if ( jto!=0) dt=(1.0/jto)*(j*sz-jold.total());


}


Value hi=0; Value hf=0; Value hclip=0;
Value ei=0; Value ef=0; Value eclip=0;
for ( IdxTy i=0; i<sz; ++i) 
{
// dj is actually current div charge making flux but also div mu 
// e is negative field 
//Value dj=-e[i]*(pnew[i]+nnew[i])+kT*gradY-j; 
//Value dj=-e[i]*(p[i]+n[i]+dp[i]+dn[i])+kT*gradY-j; 
// e[i]*(dp[i]+dn[i])*a=-e[i]*(p[i]+n[i])+kT*gradY-j;a
Value dtt=dt;
/*
if ( false) { if ( e[i]!=0){ if ( (i>0) &&(i<(sz-1)))
{
Value gradY=(-p[i+1]+n[i+1]-n[i-1]+p[i-1])/2.0/dx; 
Value da=(-e[i]*(p[i]+n[i])+kT*gradY-j)/(e[i]*dtt*(dp[i]+dn[i])); 
if ( da>100) da=100;
if ( da<.01) da=.01;
if ( false) { 
Value djzed=-e[i]*(p[i]+n[i])+kT*gradY-j; 
Value djzedx=-e[i]*(p[i]+n[i]+dtt*(dp[i]+dn[i]))+kT*gradY-j; 
Value djzedo=-e[i]*(p[i]+n[i]+dtt*da*(dp[i]+dn[i]))+kT*gradY-j; 
dtt=dtt*da;
 { MM_MSG(i<<" da= "<<da<<" jz= "<<djzed<<" jx= "<<(djzedx-djzed)<<" jo= " <<(djzedo-djzed) ) } 

} else 
dtt=dtt*da;

}
}}
*/


pnew[i]=p[i]+dp[i]*dtt;
nnew[i]=n[i]+dn[i]*dtt;

hi+=p[i]; ei+=n[i];
if ( pnew[i]<0) { hclip+=pnew[i]; pnew[i]=0; if (n[i]>0) pnew[i]=ni2/n[i];}
if ( nnew[i]<0 ) { eclip+=nnew[i]; nnew[i]=0; if ( p[i]>0) nnew[i]=ni2/p[i];}
hf+=pnew[i]; ef+=nnew[i];
} //i 
//MM_MSG(" q ex bounds hd= "<<(hf-hi)<<" ed= "<<(ef-ei)<<" ce= "<<eclip<<" hc= "<<hclip)

pnew[0]=p[0];
pnew[1]=p[1];
pnew[sz-1]=p[sz-1];
pnew[sz-2]=p[sz-2];

nnew[0]=n[0];
nnew[1]=n[1];
nnew[sz-1]=n[sz-1];
nnew[sz-2]=n[sz-2];


} // 4

// =============== 33333333333333333333333 ===================
static void sum_and_diff3(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p, const Value & ni2, const Value & G_, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q, const Value & J,const Value Neff)
{
const Value ni=::sqrt(ni2);
const Value j=.5*J/q;
const Value kti=1.0/kT;
const Value G2=G_;
const Grid & g= n.grid();
//const Value dxx2=dx*dx;
Values xnew(g),ynew(g),G(g);
Values x=(n)*mu;
Values y=(p)*mu;
for ( IdxTy i=0; i<sz; ++i) { G[i]=-G2*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni);} //i 
Values e(g),e2(g) ;
// note this is NEGATIVE e field 
D(e,g,v);
D(e2,g,e);
Values gi(g);
I(gi,g,G);
Values tt(g);
tt=((gi+e*x+j)*kti);
//for ( IdxTy i=0; i<sz; ++i) { if ( tt[i]<0) tt[i]=0;} //i 
I(xnew,g,tt);
Values tt2(g);
tt2=((gi+e*y+(-j))*kti);
//for ( IdxTy i=0; i<sz; ++i) { if ( tt[i]<0) tt[i]=0;} //i 
I(ynew,g,tt2);
// xnew is now actually "x" but without boundary condictions 
// note that ndeff was multiplied by mu by caller lol
//xnew.boundaries(n[0]*mu,n[sz-1]*mu);
//ynew.boundaries(p[0]*mu,p[sz-1]*mu);

xnew.boundaries(n[0]*mu);
//ynew.boundaries(-ynew[sz-1]+Neff);
ynew.boundaries(-ynew[sz-1]+p[sz-1]*mu);
//xnew[0]=Neff;
//ynew[sz-1]=Neff;


//ynew.boundaries(Neff,-Neff); // there is only one integral however lol 


xnew.axby(nnew,ynew,1.0/mu,0);
xnew.axby(pnew,ynew,0,1.0/mu);
static IdxTy iter=0; 
for ( IdxTy i=0; i<sz; ++i) 
{ 
if ( pnew[i]<0) pnew[i]=0; if ( nnew[i]<0) nnew[i]=0; 
// this can make pn<ni2 and make a huge mess although
// once p/n is extreme, noise in x and y may be a pronlem 
//if ( pnew[i]>1e21) pnew[i]=1e21; if ( nnew[i]>1e21) nnew[i]=1e21; 

 MM_MSG(" iter= "<<iter<< " i= "<<i<<" nn= "<<nnew[i]<<" pn= "<<pnew[i]<<" x= "<<x[i]<<" y= "<<y[i]<<" g= "<<gi[i]<<" G= "<<G[i]<<" p*n= "<<(p[i]*n[i]/ni2)<<" n= "<<n[i]<<" p= "<<p[i]<<" e= "<<e[i])

} //i 
nnew[0]=n[0];
pnew[0]=p[0];
nnew[sz-1]=n[sz-1];
pnew[sz-1]=p[sz-1];
++iter;


} // sum_and_dif3f



// ===============================   2222222222222222222 =======================
static void sum_and_diff2(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p, const Value & ni2, const Value & G_, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q, const Value & J,const Value Neff)
{
const Value ni=::sqrt(ni2);
const Value j=J/q;
const Value kti=1.0/kT;
const Value G2=2.0*G_;
const Grid & g= n.grid();
//const Value dxx2=dx*dx;
Values xnew(g),ynew(g),G(g);
Values x=(n+p)*mu;
Values y=(n-p)*mu;
for ( IdxTy i=0; i<sz; ++i) { G[i]=-G2*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni);} //i 
Values e(g),e2(g) ;
// note this is NEGATIVE e field 
D(e,g,v);
D(e2,g,e);
Values gi(g);
I(gi,g,G);
I(xnew,g,(gi+e*y)*kti);
// xnew is now actually "x" but without boundary condictions 
// note that ndeff was multiplied by mu by caller lol
//xnew.boundaries(Neff,Neff);
xnew.boundaries(Neff);
xnew[0]=Neff;
xnew[sz-1]=Neff;
for ( IdxTy i=0; i<sz; ++i) { if ( xnew[i]<0) xnew[i]=0;} //i 
I(ynew,g,((e*x+j)*kti));
ynew.boundaries(Neff); // this also should have a second boundary"?
//Value base=0;
for ( IdxTy i=0; i<sz; ++i) {
// should probably rest base, the eqn is delta 

 if ( ynew[i]<(-xnew[i])) {
//const Value dy=ynew[i]+xnew[i];
//for ( IdxTy j=i+1; j<sz; ++j) ynew[j]=ynew[j]-dy;
ynew[i]=-xnew[i]; 

}
 else if ( ynew[i]>(xnew[i])){ 
//const Value dy=ynew[i]-xnew[i];
//for ( IdxTy j=i+1; j<sz; ++j) ynew[j]=ynew[j]-dy;

ynew[i]=xnew[i]; 

}

} //i 


//ynew.boundaries(Neff,-Neff); // there is only one integral however lol 


xnew.axby(nnew,ynew,.5/mu,.5/mu);
xnew.axby(pnew,ynew,.5/mu,-.5/mu);

for ( IdxTy i=0; i<sz; ++i) 
{ 
if ( pnew[i]<0) pnew[i]=0; if ( nnew[i]<0) nnew[i]=0; 
// this can make pn<ni2 and make a huge mess although
// once p/n is extreme, noise in x and y may be a pronlem 
//if ( pnew[i]>1e21) pnew[i]=1e21; if ( nnew[i]>1e21) nnew[i]=1e21; 

 MM_MSG("i= "<<i<<" pn= "<<pnew[i]<<" nn= "<<nnew[i]<<" xn= "<<xnew[i]<<" yn= "<<ynew[i]<<" x= "<<x[i]<<" y= "<<y[i]<<" g= "<<gi[i]<<" G= "<<G[i]<<" p*n= "<<(p[i]*n[i]/ni2)<<" n= "<<n[i]<<" p= "<<p[i]<<" e= "<<e[i])

} //i 
nnew[0]=n[0];
pnew[0]=p[0];
nnew[sz-1]=n[sz-1];
pnew[sz-1]=p[sz-1];


} // sum_and_diff2



static void sum_and_diff(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p, const Value & ni2, const Value & G_, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q, const Value & J,const Value Neff)
{
Values x=(n+p)*mu;
Values y=(n-p)*mu;
const Value ni=::sqrt(ni2);
const Value j=J/q;
const Value kti=1.0/kT;
const Value G2=2.0*G_;
const Grid & g= n.grid();
//const Value dxx2=dx*dx;
Values xnew(g),ynew(g),G(g);
// const IdxTy sz=g.size();
for ( IdxTy i=0; i<sz; ++i) { G[i]=G2*(p[i]*n[i]-ni2)/(p[i]+n[i]+2*ni);} //i 
Values e(g),e2(g),dx2(g) ,temp(g),temp2(g);
// note this is NEGATIVE e field 
D(e,g,v);
D(e2,g,e);
D(temp,g,y);
I(temp2,g,e*y);

//Values d2x=e2.times(x).plus(G2).plus(e.mul(kti).times(e.times(x).add(j))).mul(kti);

// thisone was sort of working 
//Values d2x=(((e2*y)+G)+(e*temp))*(kti);
Values d2x=((G))*(kti);
Values t1(g),t2(g);
I(t1,g,d2x);
I(t2,g,t1);
t2=t2+temp2;
Value x0=t2[1];
Value xf=t2[sz-2];
// now set boundaries...
Value a=x0-Neff;
Value b=(Neff-xf+a)/(sz-3.0);

Values dy=((e*x)+(j))*(kti);
I(temp,g,dy);
Value ay=temp[1]-Neff;
//Value by=(-3e16-temp[sz-1]+ay)/(sz-3.0);


Value base=y[0];
for ( IdxTy i=1; i<(sz-1); ++i) {


// finite diff integration 
// xnew[i]=.5*(x[i-1]+x[i+1]-d2x[i]*dxx2);
xnew[i]=t2[i]-a+b*(i-1);

if ( xnew[i]<0) xnew[i]=0;
// ynew[i]=.5*(dy[i+1]+dy[i-1])*dx+base; 
ynew[i]=temp[i]-ay; // -base+y[0]; // +by*(i-1);
if ( ynew[i]>xnew[i]) { ynew[i]=xnew[i]; base=xnew[i];}
else if ( ynew[i]<(-xnew[i])) { ynew[i]=-xnew[i]; base=-xnew[i];}
else base=ynew[i];


MM_MSG("i= "<<i<<" x= "<<x[i]<<" d2x= "<<d2x[i]<<" xnew= "<<xnew[i]<<" ynew= "<<ynew[i]<<" base= "<<base)

} //i 
ynew[0]=y[0];
ynew[sz-1]=y[sz-1];
xnew[0]=x[0];
xnew[sz-1]=x[sz-1];

// do not assign? 
xnew.axby(nnew,ynew,.5/mu,.5/mu);
xnew.axby(pnew,ynew,.5/mu,-.5/mu);

for ( IdxTy i=0; i<sz; ++i) 
{ if ( pnew[i]<0) pnew[i]=0; if ( nnew[i]<0) nnew[i]=0; 

MM_MSG("i= "<<i<<" pnew= "<<pnew[i]<<" nnew= "<<nnew[i])

} //i 
}




static void srh_total_x(Values & nnew, const Values & v,const Values & n,const Values & p,const Values & rhodiveps, const Values & e , const Value & ni2, const Value & G, const Value & Gsrh,const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q)
{

const Value dx2=dx*dx;
const Value kT2=kT*2.0;
const Value c1=-G*dx2/mu;
const Value ni=::sqrt(ni2);
for ( IdxTy i=2; i<(sz-1); ++i)
{

//const Value poc=-c1*p[i]-kT2+(2.0*v[i]-v[i-1]-v[i+1]);
//const Value gi=c1*(-ni2);
//const Value gi=c1*(n[i]*p[i]-ni2);
// make e as grad v note sign change 
//pnew[i]=1.0/poc*(c1*gi-(-.5*dx*e[i]-kT)*(p[i+1]-p[i-1]));

//const Value dv=v[i+1]-v[i-1];
//const Value dn=n[i+1]-n[i-1];
const Value sn=n[i+1]+n[i-1];

// this is probably not consistent order, as product of derivatives
// should saparate but had no impact on the "fixed" results.
//const Value poc2=-c1*p[i]-kT2+(v[i]-.5*v[i-1]-.5*v[i+1]);
//nnew[i]=1.0/poc2*(gi-.5*(v[i]*n[i+1]+n[i-1]*v[i]-n[i+1]*v[i+1]-n[i-1]*v[i-1])-kT*(sn));
//const Value A=.25*dv*dn-kT*sn+G*dx2/mu*ni2;
const Value A=G*dx2/mu/ni2-.5*(v[i]*n[i+1]+n[i-1]*v[i]-n[i+1]*v[i+1]-n[i-1]*v[i-1])-kT*(sn);

//const Value B= -c1*p[i]-kT2+(2.0*v[i]-v[i-1]-v[i+1]);
const Value B= -c1*p[i]-kT2+(v[i]-.5*v[i-1]-.5*v[i+1]);

const Value a=B;
const Value b=B*p[i]+B*2*ni-A+p[i]*dx2*Gsrh/mu;
const Value c=-ni2*dx2/mu*Gsrh-A*(p[i]+2*ni);
//nnew[i]=1.0/poc*(gi+.25*(dn)*(dv)-kT*(sn));
if ( b<0) nnew[i]=.5/a*(-b+::sqrt(b*b-4.0*a*c));
else nnew[i]=.5/a*(-b-::sqrt(b*b-4.0*a*c));




if ( nnew[i]<0) nnew[i]=0;


} // for

} // srh 


static void srh_total(Values & nnew,  Values & pnew, const Values & v,const Values & n,const Values & p,const Values & rhodiveps, const Values & e , const Value & ni2, const Value & G, const Value & Gsrh, const Value & mu,const Value & dx, const IdxTy sz, const Value kT, const Value & q)
{
// this should use mu p and then mu n 
srh_total_x(nnew,v,n,p,rhodiveps,e,ni2,G,Gsrh,mu,dx,sz,-kT,q);
srh_total_x(pnew,v,p,n,rhodiveps,e,ni2,-G,-Gsrh,mu,dx,sz,kT,q);

}




// update charge assumign fixed currnet zed
static void fixed_j(const Values & v,Values & n,Values & p,Values & nnew, Values & pnew,const GeneralParameters & gp, const IdxTy sz, const Value kT)
{

for ( IdxTy i=2; i<(sz-1); ++i)
{
// eh this makes zero difference in the result lol 
//Value dvt=::exp(gp.q()*(v[i]-v[i-1])/kT);
Value dvt=::exp((v[i]-v[i-1])/kT);
Value dvtf=::exp((v[i+1]-v[i])/kT);
//MM_MSG(" ASS FUDD "<<dvt<<" "<<v[i]<<" "<<v[i-1]<<" "<<(gp.q()*(v[i]-v[i-1])/kT))
if ( dvt>2) dvt=2;
if ( dvt<.5) dvt=.5;
Value n1=n[i-1]*dvt;
Value n2=n[i+1]/dvtf;
nnew[i]=::sqrt(n1*n2);
// this forgets about the fermi level and doping etc, 
// the junction creeps right lol 
//nnew[i]=::sqrt(n[i-1]*n[i-1])*dvt;

n1=p[i-1]/dvt;
n2=p[i+1]*dvtf;
pnew[i]=::sqrt(n1*n2);
//pnew[i]=::sqrt(p[i-1]*p[i-1])/dvt;


} // i
} // fixed j  

// efn and efp need to be 
static void equill(const Values & v,Values & n,Values & p,Values & nnew, Values & pnew,const GeneralParameters & gp, const IdxTy sz, const Value & kT,const IdxTy junction, const Value & Na, const Value & Nd,const Value & Nc, const Value & Nv, const Value & efn, const Value & efp, const Value & Eg  )
{

/// this should just be two loops, right an ledt
//for ( IdxTy i=2; i<(sz-1); ++i)
for ( IdxTy i=2; i<(junction); ++i)
{
// these are both wrt the Ec level, so efn needs to be fixed 
const Value dv=v[i]-v[0];
nnew[i]=Nc*::exp((-Eg+efn+dv)/kT);
pnew[i]=Nv*::exp((-efn-dv)/kT);

} // left
for ( IdxTy i=junction; i<(sz-1); ++i)
{
const Value dv=-v[i]+v[sz-1];
nnew[i]=Nc*::exp((-Eg+efp-dv)/kT);
pnew[i]=Nv*::exp((-efp+dv)/kT);



} // right  
} // equill  








static Value fd(const Value & arg)
{ return 1.0/(1+::exp(arg)); }

static Value rhox(const Value & Eft, const Value & Na, const Value & Eat, const Value & Nd, const Value & Edt, const Value & Nc, const Value & Nv,  const Value & Ect, const Value & Evt)
{

const Value Na_=Na* fd(Eat-Eft);
const Value Nd_=Nd*fd(Eft-Edt);
const Value n=Nc*::exp(Eft-Ect);
const Value p=Nv*::exp(Evt-Eft);
return Na_-Nd_+n-p;

}

//fermi_level( Na, Ea,  Nd,  Ed,  Nc,  Nv,  Vt, Ec,Ev)

static Value fermi_level(const Value Na, const Value Ea, const Value Nd, const Value Ed, const Value Nc, const Value Nv, const Value Vt, const Value Ec, const Value Ev)
{

const Value Eat=Ea/Vt;
const Value Edt=Ed/Vt;
const Value Ect=Ec/Vt;
const Value Evt=Ev/Vt;

const Value Ef_tol=1e-8;

const Value rhov=rhox(Evt,Na,Eat,Nd,Edt,Nc,Nv,Ect,Evt);
const Value rhoc=rhox(Ect,Na,Eat,Nd,Edt,Nc,Nv,Ect,Evt);
const bool updown=(rhov>0)&&(rhoc<0);
const bool downup=(rhov<0)&&(rhoc>0);
if (( updown&&!downup)||(downup&&!updown))  {}
else MM_MSG(" boundary for fermi confusing, rhov= "<<rhov<<" rhoc= "<<rhoc);

Value Efmin=Ev;
Value Efmax=Ec;
Value Ef=(Efmax+Efmin)/2.0;
while ( true)
{
const Value Eft=Ef/Vt;
const Value rho=rhox(Eft,Na,Eat,Nd,Edt,Nc,Nv,Ect,Evt);
if ( rho>0) { if (updown) Efmin=Ef; else Efmax=Ef; } 
else { if (updown) Efmax=Ef; else Efmin=Ef; } 
Ef=(Efmax+Efmin)/2.0;
if ( Efmax-Efmin<Ef_tol)
{
MM_MSG(" returning Ef= "<<Ef<<" with rho= "<<rho<<" and rhov= "<< rhov<<" rhoc= "<<rhoc) 
 break; 
}

} // true

return Ef;
}





// try to get the debye tails for an abrupt junction
// using depletion approx as starting point
static void test_abrupt( const IdxTy npts, const Value & dist, const Value
& Eext_, const Value tol ,const IdxTy iter_count=4000)
{
const bool dump_each=false;
// put the junction in the center, go back and forth 
Grid g(0,dist,npts);
// these are now changed due to resample
IdxTy junction=npts>>1;
IdxTy sz=npts;
g.check();
GP gp;
//Values v(g),vx(g), n(g),p(g),rho(g),n0g(g),p0g(g),nnew(g),pnew(g);
Values v(g),vx(g), n(g),p(g),nnew(g),pnew(g),Nteff(g);
//Values e(g),gradj(g),gradn(g);
//Values temp(g);
//Value Eext=Eext_;
// put Na or p side on the right 
//const Value Na=1e16;
// forgot dielectic const, lol
const Value Na=3e16;
//const Value Nd=3e15;
const Value Nd=Na;
const Value kT=gp.vt();
const Value eps=gp.epsilon();
const Value Nc=2.8e19;
const Value Nv=1.83e19;
const Value ni=1.4e10;
const Value ni2=ni*ni;
const Value Ev=0;
const Value egap=1.1;
const Value Ec=Ev+egap;
// the efn is measured from confuction band
// wtf, jkust make them both from Ec and figure it out later lol
//const Value efn=kT*::log(Nd/Nc);
// relative to Ec...
//const Value efn=Ev-Ec+fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efn=fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efp=fermi_level( Na, Ev+.001,  0,0 ,  Nc,  Nv,  kT , Ec,Ev);
// from valuence band 
//const Value efp=kT*::log(Na/Nv);
//const Value vlim=egap;
//const Value Ev=-egap;
//const Value builtin=::fabs(efn+efp+egap);
// the efn is now refernced to Ec just like the efp so just subtact 
// this is positive, as efn>efp 
const Value builtin=(efn-efp); // it looks forward biases lol 

//const Value phi_target=0; // applied voltage measured

Value phi=0; // applied voltage measured
// should assert uniform grid
// but now we resample 
Value dx=g.min_dx();
Value dxf=dx*dx/eps*gp.q();
// try to guess the depltion region
const Value xp=::sqrt(2.0*builtin*Nd/(Na*(Na+Nd))*eps/gp.q());
const Value xn=::sqrt(2.0*builtin*Na/(Nd*(Na+Nd))*eps/gp.q());
Value Ezed= .5*gp.q()*Na*xp/eps;
 Ezed+= .5*gp.q()*Nd*xn/eps;
//Value Ezed_max=Ezed*10;
//Value Ezed_min=-Ezed_max;
MM_MSG("builtin= "<<builtin<<" efp= "<<efp<<" efn= "<<efn<<" xp= "<<xp<<" xn= "<<xn)
//const Value phitol=::fabs(builtin)*1e-7;
IdxTy cnt=0;
const Value Nd2=Nd/2;
const Value Na2=Na/2;
const Value Ndeff=Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Npdeff=-Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Naeff=Na2+::sqrt(Na2*Na2+ni2);
const Value Nnaeff=-Na2+::sqrt(Na2*Na2+ni2);
for ( IdxTy i=0; i<(sz); ++i)
{
const bool left=(i<junction);
// this is negaive of "q"
Value Nt=(left)?(-Nd):(Na);
Nteff[i]=Nt;
}

// this integrates out to about the builtin voltage but it
// seems to not sum the charge exactly to zero so the end drifts
// back but this goes away with finer grid, 
Value qini=0;
bool was_depleted=false;
for ( IdxTy i=0; i<sz; ++i)
{
const bool left=(i<junction);
const Value dn=g[junction]-xn-g[i];
const Value dp=g[junction]+xp-g[i];
const bool depleted= ( left)?(g[i]>(g[junction]-xn)):(g[i]<(g[junction]+xp));
if ( !depleted){
if ( left) { n[i]=Ndeff; p[i]=Npdeff;v[i]=0; qini+=n[i]-p[i]-Nd;   }
else { p[i]=Naeff; n[i]=Nnaeff; v[i]=-builtin; qini+=n[i]-p[i]+Na;}
}

if ( depleted) { n[i]=ni; p[i]=ni; 
if ( left  ) {v[i]=-.5*gp.q()*Nd/eps*dn*dn; qini-=Nd; } 
else { v[i]=-builtin+.5*gp.q()*Na/eps*dp*dp; qini+=Na; } 
} // depleted
// dump the charge here for now 
if ( !depleted&&was_depleted)
{
// simplyusing these values kept charge at 1e12 instead of creeping
// to 1e16 in a few iters. 
const Value Na2x=(Na+3e-13*qini)/2.0;
const Value Naeffx=Na2x+::sqrt(Na2x*Na2x+ni2);
const Value Nnaeffx=-Na2x+::sqrt(Na2x*Na2x+ni2);
p[i-1]=Naeffx; n[i-1]=Nnaeffx;
//p[i]+=qini;
if (false) if ( p[i]!=0)
{
Value qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;
qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;
qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;

}


//n[i]-=0*qini/2;
} 
was_depleted=depleted;
} //i init
// if this did not achieve neutrality, dump more charge somewhere. 
const bool some_resamples=!true;
IdxTy rescnt=0;
Grid * gnew=NULL;
while ( true) // phi-rho loop 
{
const bool do_resample=some_resamples&&(rescnt<6)&&((cnt%1000)==999);
if ( do_resample)
{
const IdxTy sznew=g.size()*2;
MM_MSG(" begin resample sz="<<sznew<<" "<<g.size())
//Grid g(0,dist,npts);
Grid * gold= gnew;
gnew= new  Grid(0,dist,sznew);
MM_TRACE
Grid & gn= * gnew;
// we can not give ctor address of a temp lol 
// except now the next time through, wtf... 
Values xx(gn);
xx=v.resample2x(gn); v=xx;
xx=vx.resample2x(gn); vx=xx;
xx=p.resample2x(gn); p=xx;
xx=pnew.resample2x(gn); pnew=xx;
xx=n.resample2x(gn); n=xx;
xx=nnew.resample2x(gn); nnew=xx;
xx=Nteff.resample2x(gn); Nteff=xx;
junction=junction*2;
sz=gn.size();
g=gn;
dx=g.min_dx();
dxf=dx*dx/eps*gp.q();
// try to guess the depltion region
++rescnt;
delete gold;
MM_MSG(" done resampl "<<rescnt)
} // resamples

v[0]=0;
v[1]=0;
Value qtot=0;
// this does not gie good asymptotes because total charge is not zero.
// that should occur automatically. 
for ( IdxTy i=0; i<sz; ++i) vx[i]=v[i];

const bool udate_finite_diff=false;
if ( udate_finite_diff)
{
for ( IdxTy i=2; i<(sz-1); ++i)
//for ( IdxTy i=1; i<(sz-2); ++i)
{
// this code is notusednow anyway... 
const bool left=(i<junction);
//Value Nt=(left)?(-Nd):(Na);
Value Nt=Nteff[i];
const Value nq=Nt+n[i]-p[i];
vx[i]=.5*(v[i-1]+v[i+1]-dxf*(nq));
//vx[i+1]=2*v[i]-v[i-1]+dxf*(nq);

qtot+=nq;
if ( dump_each) MM_MSG(cnt<<" "<<left<<" depl "<<i<<" "<<g[i]<<" v= "<<v[i]<<" n= "<<n[i]<<" p= "<<p[i]<<" n+p= "<<(n[i]+p[i])<<" "<<(n[i]*p[i]/ni2)<<" qtot= "<<qtot<<" nq= "<<nq)
// update n and p for zero J constraint 
}
} // update+finite
else
{
const bool force_voltage=!true;
const bool force_voltage_bound=!true;
const bool force_neutral=!true;
Values ee(g),rho(g),qn(g);
Value qnet=0;
for ( IdxTy i=0; i<(sz); ++i) {
//const bool left=(i<junction);
//Value Nt=(left)?(-Nd):(Na);
const Value Nt=Nteff[i]; 
const Value nq=Nt+n[i]-p[i];
qn[i]=nq;
// this starts out ok and drifts off by Nx times dx
qnet+=nq;
rho[i]=nq*gp.q()/eps;
}
if ( force_neutral)
{
Value qnet2=0;
const Value nq=-qnet/sz;
//MM_MSG("ASSFUDD "<<nq<<" qnet= "<<qnet)
for ( IdxTy i=0; i<(sz); ++i) {
if ( nq<0) p[i]-=nq; else n[i]+=nq;
rho[i]+=nq*gp.q()/eps;
qnet2+=rho[i];
}
//MM_MSG(" qnet was "<<qnet<<" now is "<<((eps/gp.q())*qnet2))
} // force_neutral

I(ee,g,rho);
I(vx,g,ee);
// boundary conditions should create field of zer at both ends..
// if not there is a net charge
Value a=vx[0];
Value b=(vx[sz-1]-vx[0]+builtin)/(sz-1.0);
if ( (cnt&255) == 0 ) 
{MM_MSG(cnt<<" qnet= "<<qnet<<" eerror= "<<b<<" "<<vx[sz-1])}

if (force_voltage) for ( IdxTy i=0; i<(sz); ++i) { vx[i]-=a+b*i; } 
if (force_voltage_bound) { vx[sz-1]=-builtin; vx[sz-2]=vx[sz-1];} 


//vx[0]=v[0]; vx[1]=v[1];
//vx[0]=v[0]; vx[1]=v[1];
// this is slower than the computations, not surprising
if (false) for ( IdxTy i=0; i<(sz); ++i) {
MM_MSG(" "<<cnt<<" dumpin i= "<<i<<" ee= "<<ee[i]<<"  rho= "<<rho[i]<<" vx= "<<vx[i]<<" qnet= "<<qn[i])
}


}
const bool do_backwards_too=!true;
if ( do_backwards_too)
{
for ( IdxTy i=(sz-2); i>1; --i)
{
const bool left=(i<junction);
Value Nt=(left)?(-Nd):(Na);
const Value nq=Nt+n[i]-p[i];
vx[i]=.5*vx[i]+.5*(.5*(v[i-1]+v[i+1]-dxf*(nq)));

} // i 

} // backwards

const Value fold=.1;
const Value fnew=1.0-fold;
for ( IdxTy i=2; i<sz; ++i) v[i]=fold*v[i]+fnew*vx[i];

const bool fixed_current=!true; 
const bool fixed_current_too=!true; 
const bool use_equ=false; // !fixed_current; 
const bool fixed_x=false; // !fixed_current; 
const bool use_srh=false; // !fixed_current; 
const bool use_sum_and_diff=!false; // !fixed_current; 
if ( use_sum_and_diff){

const Value G=1e14;
const Value mu=gp.mu();
const Value q=gp.q();

 //sum_and_diff(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
// sum_and_diff2(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
 //sum_and_diff3(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
 // sum_and_diff4(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0,Ndeff*mu);
 sum_and_diff5(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0,Ndeff*mu, Nteff);
}

if ( fixed_x)
{
//const Value G=400000000;
// negative values explode, positive sign is right AFAICT
// 1e2 with the finit diff, 1e13 with the integrator lol 
//const Value G=1e2;
//const Value G=1e13;
//const Value G=1e15;
const Value G=1e12;
//const Value G=00;
const Value mu=gp.mu();
const Value q=gp.q();
Values rhodiveps(g),e(g);
const bool need_garbage=false;
if ( need_garbage)
{
D(e,g,v); // sign is wrong here. 
for ( IdxTy i=0; i<sz; ++i) {
const bool left=(i<junction);
Value Nt=(left)?(-Nd):(Na);
rhodiveps[i]=-q*(Nt+n[i]-p[i])/gp.epsilon();
}
} // garbage
// q negative"? 
fixed_j_total(nnew,pnew,v,n,p,rhodiveps,e,ni2,G,mu,dx,sz,kT,q);
}

if ( use_srh)
{
const Value G=1e17+0*1e2;
const Value Gsrh=1e18;
const Value mu=gp.mu();
const Value q=gp.q();
Values junk;

srh_total(nnew,pnew,v,n,p,junk,junk,ni2,G,Gsrh,mu,dx,sz,kT,q);
} // srh 

// this causes junction creep 
if ( fixed_current) fixed_j(v,n,p,nnew,pnew,gp,sz,kT);
// efn and efp referenced to Ec and Ev, 
if ( use_equ ) equill(v,n,p,nnew,pnew,gp,sz,kT,junction,Na,Nd,Nc,Nv,efn,efp,egap);
//Value d1=0;
//Value d2=0;
//for ( IdxTy i=0; i<(sz); ++i) d1+=(n[i]-nnew[i])*(n[i]-nnew[i]);

for ( IdxTy i=2; i<(sz-1); ++i)
{ n[i]=(fold*n[i]+fnew*nnew[i]); p[i]=(fold*p[i]+fnew*pnew[i]); } 


if ( fixed_current_too) { fixed_j(v,n,p,nnew,pnew,gp,sz,kT);
for ( IdxTy i=2; i<(sz-1); ++i)
{ n[i]=(.75*n[i]+.25*nnew[i]); 
p[i]=(.75*p[i]+.25*pnew[i]);
//if ( n[i]>0) p[i]=ni2/n[i];
 } 
} 

++cnt;
phi=v[0]-v[sz-1];
//const Value dphi=phi-builtin;
const bool converged=true; // (::fabs(dphi)<phitol);
if ( (cnt%31) ==0 ) MM_MSG("cnt "<<cnt<<" phi "<<phi<<" v0= "<<v[0]<<" vd= "<<v[sz-1]<<" builtin= "<<builtin)

if (converged&&(cnt>iter_count))
{ 
//static void D(Values & d, const Grid & g, const Values & n)
Values e(g), grn(g), grp(g);
D(e,g,v); D(grn,g,n); D(grp,g,p);
const Value f1=gp.q()*gp.mu();
for ( IdxTy i=0; i<sz; ++i){
const Value Jn=f1*(-n[i]*e[i]+kT*grn[i]);
const Value Jp=f1*(-p[i]*e[i]-kT*grp[i]);
MM_MSG("cnt= "<<cnt<<" i= " <<i<<" "<<g[i]<< " phi= "<<v[i]<<" n= "<<n[i]<<" p= "<<p[i]<<" Jp= "<<(Jp)<<" Jn= "<<(Jn)<<" Jt= "<<(Jn+Jp)<<" ndrift= "<<(-f1*n[i]*e[i])<<" "<<(f1*kT*grn[i])) 

}
break;
}
} // phi-rho loop 
delete gnew;
//const bool debug=false;
#ifdef OLD_JUNK
while (false)
{
// set boundary condition, assume minmum points although
// shouldbe asseted earlier.
// this now uses junction and junction-1 to
// set a field.
// this drifts upward, jj
v[junction]=0; // .5*(v[0]+v[sz-1]);;
v[junction-1]=v[junction]-Ezed*dx;
 

for ( IdxTy i=junction+1; i<sz; ++i)
{
// this "v" has to be measured from the bulk 
Value vbend=-(v[i-1]); // -v[sz-1]);
// instead, use FD stats?
if ( vbend>vlim) vbend=vlim;
if ( vbend<(-vlim)) vbend=-vlim;

const Value phin=(-efp-egap-vbend)/kT;
const Value phip=(efp+vbend)/kT;
const Value n0=Nc*::exp(phin);
 Value p0=Nv*::exp(phip);
if ( p0>Na) p0=Na; 
MM_MSG(g[i]<<"  vbend= "<<vbend<<" phin= "<<phin<<" phip= "<<phip<<" n0= "<<n0<<" p0= "<<p0<<" np/ni^2= "<<(n0*p0/1.4e10/1.4e10)<<" n+p="<<(n0+p0))
n0g[i]=n0; p0g[i]=p0;
//MM_MSG(" p size "<<i<<" no "<<n0<<" p0 "<<p0<<" vi1 "<<v[i-1])
// put Na or p side on the right 
v[i]=2*v[i-1]-v[i-2]+dxf*(Na+n0-p0);

} // integrate right
// sz is now the reference level, 
for ( IdxTy i=junction+1; i<sz; ++i) v[i]=v[i]-v[sz-1];
for ( IntTy  i=junction; i>=0; --i)
{
Value vbend=-(v[i+1]);
// instead, use FD stats?
if ( vbend>vlim) vbend=vlim;
if ( vbend<(-vlim)) vbend=-vlim;

const Value phin=(efn-vbend)/kT;
const Value phip=(-efn-egap+vbend)/kT;
MM_MSG(" - vbend= "<<vbend<<" phin= "<<phin<<" phip= "<<phip)
 Value n0=Nc*::exp(phin);
if ( n0>Nd) n0=Nd;
 Value p0=Nv*::exp(phip);
n0g[i]=n0; p0g[i]=p0;
v[i]=2*v[i+1]-v[i+2]+dxf*(-Nd+n0-p0);
} // integrateleft

phi=v[0]-v[sz-1];
// this needs to get the sign right, 
const Value dphi=phi-builtin;
const bool converged=(::fabs(dphi)<phitol);
if ( phi>builtin) Ezed_min=Ezed;
//if ( phi<builtin) Ezed_min=Ezed;
else  Ezed_max=Ezed;
Ezed=(Ezed_min+Ezed_max)/2;

MM_MSG("cnt "<<cnt<<" phi "<<phi<<" v0= "<<v[0]<<" vd= "<<v[sz-1]<<" builtin= "<<builtin<<" Ezed= "<<Ezed<<" "<<Ezed_min<<" "<<Ezed_max )
if (converged||(cnt>100))
{ 
for ( IdxTy i=0; i<sz; ++i){MM_MSG("cnt= "<<cnt<<" i= " <<i<<" "<<g[i]<< " phi= "<<v[i]<<" n0= "<<n0g[i]<<" p0= "<<p0g[i]) }
break;
}
++cnt; 
} // true
#endif // old_junk
} // test_abrupt

}; // dummy ?

}; // ns  mjm_dd_base


#endif


