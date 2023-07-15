#ifndef MJM_DEFECT_LEVELS_H__
#define MJM_DEFECT_LEVELS_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
// finding the equillibrium levels
#include "mjm_generic_iterators.h"
// this also needs things suitable for pseudo FEM

#include <algorithm>

/*
introduced in  mjm_shooting10.h :  terrible problems with b2b only
and want to see if these dynamics are any better  

*/



class mjm_defect_level 
{
/*
 Properties of a vriable occupany fixed thing. 
*/

protected:
class Tr
{
public:

typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
typedef std::string StrTy;
typedef std::stringstream Ss;

}; // Tr

// this needs a StrTy def somewhere? Where wtf? 
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;
typedef Tr::Ss Ss;


public:

 mjm_defect_level(const StrTy & name
	,const D &  phi_e, const D & phi_h, const D & e_e
	,const D &  e_h, const D &  rho_empty,const D & rho_full  )
: m_name(name),m_phi_e(phi_e),m_phi_h(phi_h),m_e_e(e_e)
, m_e_h(e_h), m_rho_empty(rho_empty), m_rho_full(rho_full)
{}
 mjm_defect_level(  )
: m_name(StrTy("default")),m_phi_e(1e-6),m_phi_h(1e-6),m_e_e(1e7)
, m_e_h(1e7), m_rho_empty(0), m_rho_full(-1)
{}
 mjm_defect_level(const D & speed  )
: m_name(StrTy("default"))
,m_phi_e(1e-6*speed),m_phi_h(1e-6*speed),m_e_e(1e7*speed)
, m_e_h(1e7*speed), m_rho_empty(0), m_rho_full(-1)
{}



const StrTy name() const { return m_name; } 
const StrTy name(const StrTy & name)  {m_name=name;  return m_name; } 
const D & phi_e() const { return m_phi_e; } 
const D & phi_h() const { return m_phi_h; } 
const D & e_e() const { return m_e_e; } 
const D & e_h() const { return m_e_h; } 

// can add field dependence  multipliciity etc later
const D step( const D & nfull, const D & N, const D & dt, const D & n, const D & p) const 
{
const D nempty=N-nfull; 
// dhdt=m_phi_e*n*nempty-m_phi_h*nfull*p-m_e_e*nfull+m_e_h*nempty;
//const D  dhdt=m_phi_e*n*nempty-m_phi_h*nfull*p-m_e_e*nfull+m_e_h*nempty;
//const D  tau=m_phi_e*n*(N-nfull)-m_phi_h*nfull*p-m_e_e*nfull+m_e_h*(N-nfull);
//const D  tau=N*(m_phi_e*n+m_e_h) -nfull*(m_phi_e*n+m_phi_h*p+m_e_e+m_e_h);
const D  C=N*(m_phi_e*n+m_e_h);
const D a= (m_phi_e*n+m_phi_h*p+m_e_e+m_e_h);

const D k= C -nfull*a;
const D nfinal=(C-k*exp(-a*dt))/a;
return nfinal; 

// this is great but can violate boundaries, integrate
// exact for const params, 
//return nfull+dt*dhdt;

}

const D rho( const D & nfull, const D & N) const 
{
D x=m_rho_full*nfull+(N-nfull)*m_rho_empty;
return x;
}

// steady state occupancy 
const D nzed(  const D & N, const D & n, const D & p) const 
{
const D  C=N*(m_phi_e*n+m_e_h);
const D a= (m_phi_e*n+m_phi_h*p+m_e_e+m_e_h);
// a could be zero... 
return C/a;

}


StrTy to_string() const
{
Ss ss;
ss<<MMPR(m_name);
ss<<MMPR(m_phi_e)<<MMPR(m_phi_h)<<MMPR(m_e_e)<<MMPR( m_e_h);
ss<<MMPR(m_rho_empty)<<MMPR( m_rho_full); 

return ss.str();
}

private:
StrTy m_name;
D m_phi_e, m_phi_h, m_e_e, m_e_h;
D m_rho_empty, m_rho_full; 

}; //mjm_defect_level 
//////////////////////////////////////////////////////////////


class mjm_defect_levels 
{
/*
  Collection of defect levels. 
*/
typedef mjm_defect_levels Myt;
typedef mjm_defect_level Level;
typedef std::vector<Level> LevelInfo;
protected:
class Tr
{
public:

typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
typedef std::string StrTy;
typedef std::stringstream Ss;

}; // Tr

// this needs a StrTy def somewhere? Where wtf? 
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;
typedef Tr::Ss Ss;


public:
typedef Level defect_level;

mjm_defect_levels() {}

D net_charge( const IdxTy i) const
{
D n=0;
for (IdxTy j=0; j<m_sz; ++j)
{
n+=m_levels[j].rho(m_info(i,2*j),m_info(i,1+(2*j)));

}
return n;
}

void net_charge( MyBlock & dest, const IdxTy col,const IdxTy sz) const
{
for(IdxTy i=0; i<sz; ++i) dest(i,col)=net_charge(i);

}
void net_charge( MyBlock & dest, const IdxTy col,const IdxTy sz,const IdxTy bn) const
{
for(IdxTy i=bn; i<(sz-bn); ++i) dest(i,col)=net_charge(i);

}


void net_charge( MyBlock & dest, const IdxTy sz) const
{
for(IdxTy i=0; i<sz; ++i) dest(i)=net_charge(i);

}
// needs to find n,p, and maybe x or v 
void step( const MyBlock & fixed, const MyBlock & sol
, const D & dti, const IdxTy & points
, const IdxTy & in, const IdxTy &ip, const IdxTy & iu, const IdxTy & ix
, const IdxTy boundaries=1)
{
// for now in place is ok as there is no field etc effect
// just local n and p 
const IdxTy sz=m_sz;
const IdxTy pts=points-boundaries;
for (IdxTy i=boundaries; i<pts; ++i)
{
for (IdxTy j=0; j<sz; ++j)
{

D & nfull=m_info(i,2*j);
const D & N=m_info(i,1+(2*j));
const D & n=sol(i,in);
const D & p=sol(i,ip);
const D nnew= m_levels[j].step( nfull,  N,  dti,  n, p);
nfull=nnew;
}

} // i 
}

template <class Di> 
void step( Di & din, Di & dip, const MyBlock & fixed, const MyBlock & sol
, const D & dti 
,  const IdxTy & iu, const IdxTy & ix)
{
// for now in place is ok as there is no field etc effect
// just local n and p 
const IdxTy sz=m_sz;
while (din.ok())
{
const IdxTy i=din.i();
for (IdxTy j=0; j<sz; ++j)
{

D & nfull=m_info(i,2*j);
D & nfullnext=m_info(i+1,2*j);
// this depends on the shape functins, needs to be elsewhere 
// I think the Di provides a matrix for this 
const D & N0=m_info(i-1,1+(2*j));
const D & N1=m_info(i,1+(2*j));
const D & N2=m_info(i+1,1+(2*j));
const D & n0=din.n0(); // sol(i,in);
const D & n1=din.n1(); // sol(i,in);
const D & n2=din.n2(); // sol(i,in);
const D & p0=dip.n0(); // sol(i,in);
const D & p1=dip.n1(); // sol(i,in);
const D & p2=dip.n2(); // sol(i,in);
const D & x0=din.x0();
const D & x1=din.x1();
const D & x2=din.x2();

D nnew=0;
const IdxTy npts=11;
const D ish=1.0/(x2-x0);
const D h=(x2-x0)/(npts-1);
for( IdxTy k=0; k<npts; ++k)
{
D x=x0+h*k;
const bool left=(x<x1);
const D xbase=(left)?x0:x1;
const D xend=(left)?x1:x2;
//const D slope=(left)?(1.0/(x1-x0)):(-1.0/(x2-x1));
D eta=(x-xbase)/(xend-xbase);
D xi=1.0-eta;
// do brute force integral for now as this allows the
// gr kindetics to be unknown 
const D  n=left?(n0*eta+n1*xi):(n1*eta+n2*xi); // sol(i,ip);
const D  p=left?(p0*eta+p1*xi):(p1*eta+p2*xi); // sol(i,ip);
const D  N=left?(N0*eta+N1*xi):(N1*eta+N2*xi); // sol(i,ip);
const D w=left?((x-x0)/(x1-x0)):((x2-x)/(x2-x1));
// eh, TODO FIXME nfull needs to be interpolated too 
const D nnewdx= w*m_levels[j].step( nfull,  N,  dti,  n, p);
nnew+=2.0*nnewdx*h*ish;
//MM_MSG(MMPR(i)<<MMPR(k)<<MMPR(N)<<MMPR(nfull)<<MMPR(nnew)<<MMPR(n)<<MMPR(p)<<MMPR(nfullnext) )
} // j


nfull=nnew;
}

din.inc();
dip.inc();
} // din.ok() 
}


/////////////////////////////////////////////////////
void set(const IdxTy nlevels, const IdxTy points)
{
m_points=points;
m_sz=nlevels;
m_info.resize(m_points,2*m_sz);
m_levels.clear();
MM_MSG(" traps sizes set to "<<nlevels<<" and points= "<<points)
}

void add(const D & N)
{
Level defect;
add(N,defect);
}
void add(const D & N, const Level & defect)
{
m_levels.push_back(defect);
const IdxTy iN=m_levels.size()*2-1;
const IdxTy in=iN-1;
MM_MSG(" adding defect amount "<<N<<" at idx "<<iN<<" and occupancy at "<<in)
MM_MSG(" seeting pops to half "<<m_points<<MMPR(m_info.size()))
for (IdxTy i=0; i<m_points; ++i) 
{
m_info(i,iN)=N;
m_info(i,in)=.5*N;

}
}
// project one point as a linear combination from old defect thing
void project_from( const Myt & old, const IdxTy point,const IdxTy  i1, const IdxTy i2, const D & f1, const D & f2,const bool do_N)
{
for (IdxTy i=0; i<m_sz; ++i) 
{
const IdxTy in=2*i;
const IdxTy iN=in+1;
const D N=old.m_info(i1,iN)*f1+f2*old.m_info(i2,iN);
const D newN=m_info(point,iN);
if (do_N) m_info(point,iN)=N;
const D nfull=old.m_info(i1,in)*f1+f2*old.m_info(i2,in);
// assume that the occupancy rates are intepolated 
if (do_N) m_info(point,in)=nfull;
else m_info(point,in)=nfull*newN/N;
}
}

void project_from( const Myt & old, const IdxTy point,
const IdxTy  i1, const IdxTy i2, const IdxTy i3, const IdxTy i4,
const D & f1, const D & f2, const D & f3, const D & f4, const bool do_N)
{
for (IdxTy i=0; i<m_sz; ++i) 
{
const IdxTy in=2*i;
const IdxTy iN=in+1;
	// the init code should have set this up more accuately right? 
const D N=old.m_info(i1,iN)*f1+f2*old.m_info(i2,iN)
			+old.m_info(i3,iN)*f3+f4*old.m_info(i4,iN);
const D newN=m_info(point,iN);
if (do_N) m_info(point,iN)=N;
const D nfull=old.m_info(i1,in)*f1+f2*old.m_info(i2,in)
			+old.m_info(i3,in)*f3+f4*old.m_info(i4,in);
// assume that the occupancy rates are intepolated 
if (do_N) m_info(point,in)=nfull;
else m_info(point,in)=nfull*newN/N;
}
}


void set_to_equ(const D & n, const D & p, const Level & lev
, const IdxTy in, const IdxTy iN,const IdxTy first, const IdxTy end) 
{
//const D nzed(  const D & N, const D & n, const D & p) const 
for (IdxTy i=first; i<end; ++i) 
{
const D N=m_info(i,iN);
const D nfull=lev.nzed(N,n,p);
if (i==first) { MM_MSG(" set to equ "<<MMPR(i)<<MMPR(n)<<MMPR(p)<<lev.to_string()<<MMPR(nfull)<<MMPR(N)<<MMPR((nfull/N))<<MMPR(iN)) } 
m_info(i,in)=nfull;

}

}
void set_to_equ(const D & n, const D & p
, const IdxTy first, const IdxTy end) 
{
//const D nzed(  const D & N, const D & n, const D & p) const 
for (IdxTy i=0; i<m_sz; ++i) 
{
set_to_equ(n,p,m_levels[i],2*i, 2*i+1,first,end); 
}


}

void set_to_equ(const IdxTy first, const IdxTy end, const D & ni2, const D & Na, const D & Nd)
{
D n,p;
set_to_equ(n,p,ni2,Na,Nd,first);
set_to_equ(n,p,first,end);

}
// given ni2 and fixed dopants set up everything from equ
// and return the n and p 
void set_to_equ(D & n, D & p, const D & ni2, const D & Na, const D & Nd, const IdxTy point, const bool dbg=!false)
{
// the 
shooting_iterator si;
D nmin=1;
D nmax=ni2;
//if (Nd!=0) if (Na==0) { nmin=Nd*1e-10; nmax=Nd*10; }
//if (Na!=0) if (Nd==0) { nmin=ni2/Na*1e-10; nmax=ni2/Na*10; }
bool done=false;
const D tol=.1;
//bool params( const D & tgt,const D & inmin, const D & inmax, const D & tmax, const D & tmin, const IdxTy maxiter)
// this can take 100 iterations because t needs to divide the log
// not the linear thing... argh jjjjjjjjjj
si.params(0,nmin,nmax,tol,-tol,300);
//si.strategy(2);
si.strategy(0);
D rho=0;
D nfinal=0;
MM_MSG(" setting equ levels="<<m_sz)
while (!done)
{
n=si.guess();
p=ni2/n;
rho=p-n+Nd-Na;
for(IdxTy i=0; i<m_sz; ++i) {
const Level level= m_levels[i];
const D N=m_info(point, 2*i+1);
nfinal=level.nzed(   N,  n, p); 
//rho+=level.rho(level.nzed(   N,  n, p),N ); 
rho+=level.rho(nfinal,N ); 
//MM_MSG(MMPR(rho)<<MMPR(n)<<" si: "<<si.to_string())
}

//MM_MSG("  rho "<< MMPR(rho)<<MMPR(n)<<" si: "<<si.to_string())
// this assumes sort of that goes up with rho
done=si.result(-rho);
}

if (dbg) MM_MSG("  final rho "<< MMPR(rho)<<MMPR(n)<<MMPR(p)<<MMPR(Na)<<MMPR(Nd)<<" si: "<<si.to_string())
} //set_to_equ

StrTy to_string() const
{
Ss ss;
ss<<MMPR(m_points)<<MMPR(m_sz);
for (IdxTy i=0; i<m_sz; ++i) ss<<CRLF<<m_levels[i].to_string();
//ss<<CRLF;
return ss.str();

}
private:
MyBlock m_info;
IdxTy m_points;
LevelInfo m_levels;
IdxTy m_sz;
}; //mjm_defect_levels


#endif

