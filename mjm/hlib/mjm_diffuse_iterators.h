#ifndef MJM_DIFFUSE_ITERATORS_H__
#define MJM_DIFFUSE_ITERATORS_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"

#include <algorithm>
#define USING_PETSC
#ifdef USING_PETSC
// deviate from shooting to matrix 
#include "mjm_petsc_util.h"
#include "mjm_petsc_fd_base.h"
#endif
/*
copied from mjm_shooting7.h : trying to understand FEM problems with drift diffusion. 
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



class diffuse_iterator
{
/*
This is now a junk bin for iterating over FD points ( or a grid if you prefer )
 pretending to be fininite elemetns.
First it was just a way to make interpolation easier and define nodes consisitently
across diffusion mechanisms. It then added optional velocty info andginally
a interpolation scheme for user provided matrix and location. 

As most diffusion or drift is implemented as from-to pixels it is very wasteful
but this itor is predicated on that iterating first over i and then plus/minus
j with limits set by results. 

Some grid constant like left and right are here too. As it is an iterator, they
are handy to have.  

When iterating with a drift velocity component, j==i is included. 
For diffusion the remainder is just added back although all this can be more consistent.

*/
/*
 This is used as an element iterator for FL and drift carrier moves. However, it is not
yet used in potential and GR. The interpolation logic however is external. 
There is an optional velocity thing but no way yet to find fixed charges. 

It only deals with one carrier type and defines geometry of element and boundary
proximity. Creates pseudo element for boundary. 

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
	diffuse_iterator(const IdxTy points, const IdxTy _x, const IdxTy _in, const MyBlock & sol, const MyBlock & f)
	: m_state(~0),m_points(points), m_x(_x),m_in(_in),m_iv(~0),m_sol(sol), m_f(f) {Init();}
	diffuse_iterator(const IdxTy points, const IdxTy _x, const IdxTy _in, const IdxTy _iv,const MyBlock & sol, const MyBlock & f)
	: m_state(~0),m_points(points), m_x(_x),m_in(_in),m_iv(_iv),m_sol(sol), m_f(f) {Init();}


	bool ok() const { return m_ok; }
	bool result(const D & gdv) 
	{
		const bool bump=(gdv<=0);	
		if (bump) { ++m_state; BumpInto(); return bump; }    
		if ( Bump()) { ++m_state; BumpInto(); return true; } 
		return bump; 
	}
	// these are left,mid,right to be passed to inerpolator or fit 
	const D & x0() const { return x[0]; } 
	const D & x1() const { return x[1]; } 
	const D & x2() const { return x[2]; } 
	const D & n0() const { return n[0]; } 
	const D & n1() const { return n[1]; } 
	const D & n2() const { return n[2]; } 
	// optional drift velocty info
	const D & v0() const { return v[0]; } 
	const D & v1() const { return v[1]; } 
	const D & v2() const { return v[2]; } 
	// from 
	const int & i() const { return m_i; } 
	// to 
	const int & j() const { return m_j; } 
	// iterator state 
	const IdxTy & state() const { return m_state; } 
	// destination or j limits 
	const D & a() const { return m_a; } 
	const D & b() const { return m_b; } 
	const D dh() const { return m_b-m_a; }  // destination h 
	const D sh() const { return x[2]-x[0]; }  // src or i h value 
	// this needs a general interpolation function 
//	const D ivalue(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_i );  }
//	const D jvalue(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_j );  }

	const D ivalue0(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_i,0 );  }
	const D ivalue1(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_i,1 );  }
	const D ivalue2(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_i,2 );  }

	const D jvalue0(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_j,0 );  }
	const D jvalue1(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_j,1 );  }
	const D jvalue2(const MyBlock & data, const IdxTy col) const { return Value(data, col, m_j,2 );  }

	// these are the terminal locations of the grid over which these things are defined. 
	const D & left() const { return m_x_0; } 
	const D & right() const { return m_x_d; } 
	// this is not quite right as these sf's overlap.. doh
	const bool hits_left(const D & xmax ) const { return (x0()-xmax)<=left(); } 
	const bool hits_right(const D & xmax ) const { return (x2()+xmax)>=right(); } 

	const IdxTy iv() const { return m_iv; }
	StrTy to_string() const
	{ return state_string()+const_string(); } 
	StrTy state_string() const
	{
		Ss ss; 
		ss<<" state="<<state();
		ss<<" i="<<i();
		ss<<" j="<<j();
		ss<<" a="<<a();
		ss<<" b="<<b();
		return ss.str();
	}
	StrTy const_string() const
	{
		Ss ss; 
		ss<<" x0="<<x0();
		ss<<" x1="<<x1();
		ss<<" x2="<<x2();
		ss<<" n0="<<n0();
		ss<<" n1="<<n1();
		ss<<" n2="<<n2();
		if (m_iv!=(~0))
		{
		ss<<" v0="<<v0();
		ss<<" v1="<<v1();
		ss<<" v2="<<v2();
		}
		ss<<" dh="<<dh();
		ss<<" sh="<<sh();
		ss<<" left="<<left();
		ss<<" right="<<right();
		return ss.str();
	}

	void inc() { IncI(); } 

protected:


bool  Bump()
{
	switch (m_state)
	{
		case 0: {  if (m_j<=1) return true; --m_j; BumpJ(m_j,m_sol,m_f); return false; } 
		case 1: { ++m_j; if (m_j>=m_sz) return true; BumpJ(m_j,m_sol,m_f); return false; } 
	} ; // switch 
return false; 
}
void IncI()
{
 	++m_i; 
	BumpI(m_i,m_sol,m_f); 
}
void BumpInto()
{
		switch (m_state)
		{
			case 1: { m_j=m_i+m_j_offset; BumpJ(m_j,m_sol,m_f); break; } // now incrementing j 
			case 2: { ++m_i; BumpI(m_i,m_sol,m_f); break; } // incrementn i  
		} ; // switch 
}

virtual void Init()
{
m_sz=m_points-1;
m_ok=(m_points>3);
if (!ok()) return; 
m_i=1;
m_j=1; // really not needed 
m_j_offset=1;
if (m_iv!=(~0)) m_j_offset=0; 
m_x_0=.5*(m_f(0,m_x)+m_f(1,m_x));
m_x_d=.5*(m_f(m_sz,m_x)+m_f(m_sz-1,m_x));
BumpI(m_i,m_sol,m_f);
}

// this is still not right, need left/mid/right etc. 
virtual 	D Value(const MyBlock & data, const IdxTy col, const IdxTy i ) const
{ 
MM_ERR(" do notcall this crap ")

return .5*data(i,col)+.25*(data(i+1,col)+data(i-1,col)); }
// better
virtual	D Value(const MyBlock & data, const IdxTy col, const IdxTy i, const IdxTy node ) const
{ 
MM_ERR(" do notcall this crap ")
switch (node)
{
case 0: { return .5*(data(i,col)+data(i-1,col)); } 
case 1: { return .5*(data(i,col)+data(i-0,col)); } 
case 2: { return .5*(data(i,col)+data(i+1,col)); } 

}

return .5*data(i,col)+.25*(data(i+1,col)+data(i-1,col)); 
}
virtual void BumpJ(const int j, const MyBlock & sol, const MyBlock & f )
{
	const D tz=f(j,m_x);
	const D tm=f(j-1,m_x);
	const D tp=f(j+1,m_x);

	m_a=.5*(tz+tm);
	m_b=.5*(tz+tp);

}


virtual void BumpI(const int  i, const MyBlock & sol, const MyBlock & f )
{
	if (m_i>=m_sz) { m_state=~0; m_ok=false; return; } 
	// need to check the bumpj later
	m_j=i-1;
	const D nm=sol(i-1,m_in);
	const D nz=sol(i,m_in);
	const D np=sol(i+1,m_in);
	const D xz=f(i,m_x);
	const D xm=f(i-1,m_x);
	const D xp=f(i+1,m_x);

	x[0]=.5*(xm+xz); n[0]=.5*(nm+nz);	
	x[1]=.5*(xz+xz); n[1]=.5*(nz+nz);	
	x[2]=.5*(xp+xz); n[2]=.5*(np+nz);	
	if (m_iv!=(~0))
	{
		v[0]=.5*(sol(i-1,m_iv)+sol(i,m_iv));
		v[1]=.5*(sol(i,m_iv)+sol(i,m_iv));
		v[2]=.5*(sol(i+1,m_iv)+sol(i,m_iv));
	}
	m_state=0; //  dec j
	if (m_j<1) { ++m_state; BumpInto(); }
	else BumpJ(m_j,sol,f);
}



bool m_ok;
int m_i,m_j,m_j_offset;
IdxTy m_sz,m_state;
//D  x0,x1,x2, n0,n1,n2;
D  x[3], n[3],v[3],m_a,m_b,m_x_0,m_x_d;
const IdxTy m_points;
const IdxTy m_x,m_in,m_iv;
const MyBlock & m_sol;
const MyBlock & m_f;

}; // diffuse_iterator
//////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////



class diffuse_iterator_sf : public diffuse_iterator
{
/*
A somewhat generalized version of the parent it should allow a shape
function to make attributions more accurate and fix among other things
boundary problems. 

*/

typedef diffuse_iterator_sf Myt;
typedef diffuse_iterator Super;
typedef Super::Tr Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;
typedef Tr::Ss Ss;

public:
// this gave nan etc without Init call I guess vft is not setup yet when Super() is called  
	diffuse_iterator_sf(const IdxTy points, const IdxTy _x, const IdxTy _in, const MyBlock & sol, const MyBlock & f)
	: Super(points,_x,_in,sol,f) {Init(); }
	//: m_state(~0),m_points(points), m_x(_x),m_in(_in),m_iv(~0),m_sol(sol), m_f(f) {Init();}
	diffuse_iterator_sf(const IdxTy points, const IdxTy _x, const IdxTy _in, const IdxTy _iv,const MyBlock & sol, const MyBlock & f)
	: Super(points,_x,_in,_iv,sol,f) {Init();  }
//m_state(~0),m_points(points), m_x(_x),m_in(_in),m_iv(_iv),m_sol(sol), m_f(f) {Init();}

	const D & xp0() const { return xp[0]; } 
	const D & xp1() const { return xp[1]; } 
	const D & xp2() const { return xp[2]; } 

protected:

// instead of all these virtuals, really just need mixing coefficients.... 
void InitMat()
{
m_mat.resize(3,3);
// quicj regression test worked... 
const bool old_thing=!true;
if (old_thing)
{
// negative is hald and half
// "elements" or integration ranges are distinct uniform up to midpoints. 
m_mat(0,0)=.5; m_mat(0,1)=.5;
m_mat(1,1)=1;
m_mat(2,2)=.5; m_mat(2,1)=.5;

}
else 
{
// just take adjacent elements as boundary and use shape functions
// for attributions in integrals. 
m_mat(0,0)=1; m_mat(1,1)=1; m_mat(2,2)=1;
}
}

// convert grid points into sort of nodes that bound sort of elements lol
void xform(D * p0, D * p1, D * p2, const D & v0, const D & v1, const D & v2)
{
if (p0!=0){ *p0=m_mat(0,0)*v0+m_mat(0,1)*v1+m_mat(0,2)*v2;}
if (p1!=0){ *p1=m_mat(1,0)*v0+m_mat(1,1)*v1+m_mat(1,2)*v2;}
if (p2!=0){ *p2=m_mat(2,0)*v0+m_mat(2,1)*v1+m_mat(2,2)*v2;}

}
//void xform(&m_a,0,&m_b,f,j,m_x);
void xform(D * p0, D * p1, D * p2, const MyBlock & f, const IdxTy j, const IdxTy col )
{ xform(p0, p1, p2, f(j-1,col),f(j,col),f(j+1,col)); }

virtual void Init()
{
InitMat();
m_sz=m_points-1;
m_ok=(m_points>3);
if (!ok()) return; 
m_i=1;
m_j=1; // really not needed 
m_j_offset=1;
if (m_iv!=(~0)) m_j_offset=0; 
xform(&m_x_0,0,0,m_f,1,m_x); 
//m_x_0=m_f(0,m_x); // .5*(m_f(0,m_x)+m_f(1,m_x));
//m_x_d=m_f(m_sz-1,m_x); // .5*(m_f(m_sz,m_x)+m_f(m_sz-1,m_x));
xform(0,0,&m_x_d,m_f,m_sz-1,m_x); 
BumpI(m_i,m_sol,m_f);
}

// this is still not right, need left/mid/right etc. 
virtual 	D Value(const MyBlock & data, const IdxTy col, const IdxTy i ) const
{ 
MM_ERR(" do notcall this crap ")
return .5*data(i,col)+.25*(data(i+1,col)+data(i-1,col)); }
// better
virtual	D Value(const MyBlock & data, const IdxTy col, const IdxTy i, const IdxTy node ) const
{ 
MM_ERR(" do notcall this crap ")
switch (node)
{
case 0: { return .5*(data(i,col)+data(i-1,col)); } 
case 1: { return .5*(data(i,col)+data(i-0,col)); } 
case 2: { return .5*(data(i,col)+data(i+1,col)); } 

}

return .5*data(i,col)+.25*(data(i+1,col)+data(i-1,col)); 
}

virtual void BumpJ(const int j, const MyBlock & sol, const MyBlock & f )
{
//	const D tz=f(j,m_x);
//	const D tm=f(j-1,m_x);
//	const D tp=f(j+1,m_x);

//	m_a=tm; // .5*(tz+tm);
//	m_b=tz; // .5*(tz+tp);
	xform(&m_a,0,&m_b,f,j,m_x);
	xform(&xp[0],&xp[1],&xp[2],f,j,m_x);
}
virtual void BumpI(const int  i, const MyBlock & sol, const MyBlock & f )
{
	if (m_i>=m_sz) { m_state=~0; m_ok=false; return; } 
	// need to check the bumpj later
	m_j=i-1;
//	const D nm=sol(i-1,m_in);
//	const D nz=sol(i,m_in);
//	const D np=sol(i+1,m_in);
//	const D xz=f(i,m_x);
//	const D xm=f(i-1,m_x);
//	const D xp=f(i+1,m_x);

	xform(&x[0],&x[1],&x[2],f,i,m_x);
	xform(&n[0],&n[1],&n[2],sol,i,m_in);
//	x[0]=xm; n[0]=nm; // .5*(xm+xz); n[0]=.5*(nm+nz);	
//	x[1]=xz; n[1]=nz; // .5*(xz+xz); n[1]=.5*(nz+nz);	
//	x[2]=xp; n[2]=np; // .5*(xp+xz); n[2]=.5*(np+nz);	
	if (m_iv!=(~0))
	{
		xform(&v[0],&v[1],&v[2],sol,i,m_iv);
//		v[0]=sol(i-1,m_iv); // .5*(sol(i-1,m_iv)+sol(i,m_iv));
//		v[1]=sol(i,m_iv); // .5*(sol(i,m_iv)+sol(i,m_iv));
//		v[2]=sol(i+1,m_iv); // .5*(sol(i+1,m_iv)+sol(i,m_iv));
	}
	m_state=0; //  dec j
	if (m_j<1) { ++m_state; BumpInto(); }
	else BumpJ(m_j,sol,f);
}

MyBlock m_mat;
D xp[3];

}; // diffuse_iterator_sf


#endif

