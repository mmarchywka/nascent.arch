#ifndef MJM_FLOW_KERNELS_H__
#define MJM_FLOW_KERNELS_H__

#include "mjm_globals.h"
//#include "mjm_qfle.h"
//#include "mjm_logic_base.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
//#include "mjm_interpolation.h"
#include "mjm_instruments.h"
//#include "mjm_diffuse_iterators.h"
//#include "mjm_defect_levels.h"
#include "mjm_gauss.h"

//#include "mjm_fem_indexing.h"
//#include "mjm_convergence_tracking.h"
//#include "mjm_fick_equations.h"

#include <algorithm>
#include <map>
#include <signal.h>
#include <cstdlib>
//#define USING_PETSC
//#ifdef USING_PETSC
// deviate from shooting to matrix 
//#include "mjm_petsc_util.h"
//#include "mjm_petsc_fd_base.h"
//#endif

// libmesh crap 
//#pragma GCC diagnostic ignored "-Wunused-parameter"
//#include "mjm_libmesh.h"
// fick 
//#pragma GCC diagnostic ignored "-Wunused-parameter"



/*

2017-08-13 


*/







template <class Tparams,class Tbr > class gauss_kernel2
{
/*
Calculate the flow from "point" i or src to j or x' xp or destination.
The nood itor finds the neighors of i and makes one pass first 
before making final calculations allowing for better conservation etc. 

*/

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 

//typedef petsc_fd_base Super;

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::MySparse MySparse;

typedef Tparams FlpTy;
typedef Tbr BraTy;


class state { public: D gni; };
typedef int NodeTy;
typedef state StateTy;
typedef std::map<NodeTy,StateTy> NodeMap;
typedef mjm_gauss QuadTy;
// evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a,b,  beta)
public:
gauss_kernel2( const FlpTy & flp, const D & _dt,const D & mu
	, const D & dx, const D & dy, const BraTy & tree, CounterMap * pcm)
{
Init(flp,_dt,mu,dx,dy,pcm);
Init(tree);
}
gauss_kernel2( const FlpTy & flp, const D & _dt,const D & mu, const D & dx, const D & dy, const BraTy & tree)
{
//mp_cm=0; a=dx; b=dy;
Init(flp,_dt,mu,dx,dy,0);
Init(tree);
}

gauss_kernel2( const FlpTy & flp, const D & _dt,const D & mu, const D & dx, const D & dy)
{
//mp_cm=0; a=dx; b=dy;
Init(flp,_dt,mu,dx,dy,0);
}
~gauss_kernel2()
{
if (rejects!=0) MM_MSG("bad kernel clinks " << MMPR(rejects));
}
void set_bc( const D & minus, const D & plus) { nminus=minus; nplus=plus; } 

void set(CounterMap * pcm ) { mp_cm=pcm; } 
// a new point from the out node iterator. This may be a real mesh or a boundary pixel
// carrier velocity and number are needed but may need to expand this later 
//void point(const NodeTy _id, const D & _n0,const D & _vx,const D & _vy)
void point(const NodeTy _id, const D & _n0,const D & _vx,const D & _vy, const D & x, const D & y)
{
	m_map.clear();
	id=_id; n0=_n0; vx=_vx; vy=_vy;
px=x; py=y;
	tvx=vx*dt; tvy=vy*dt; vxt=::sqrt(tvx*tvx+tvy*tvy);
	gnorm=0;
	if (debug_result) { MM_MSG(MMPR4(id,n0,vx,vy)<<MMPR3(tvx,tvy,vxt)) }
}
/* right now we do not need the velocity at the dest node but will later and
in the case of bulk off grid may be useful 
const D v= m_quad.evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a,b,  beta)

*/

void clink( const D & v, const NodeTy _dnode,const D & _r2 , const D&  _dx, const D& _dy, const D & _ni
, const D & vxf, const D & vyf, const D bcf)
{
	//{ MM_MSG("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm)) }
	++rejects; MM_INCP_MSG(mp_cm, "kernel_clinksc1" )
 	{ MM_ONCE("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm),) }  

}
bool do_olap(const NodeTy _dnode)
{
 return noverlap;
// return noverlap&&(_dnode!=id);
}
void calculate( const NodeTy _dnode,const D & _r2 , const D& _dx, const D& _dy, const D & _ni)
{
	if (_dnode==id) return;  // remainder term is left behind 
	const D v= 
		m_quad.evaluate_node_pair(  -_dx, - _dy, 0, 0, tvx, tvy, a, b, beta,pattern,noverlap); 
	if (debug_calculate) 
		{ MM_MSG(MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm)) }
	if (v<0) { clink(v,_dnode,_r2,_dx,_dy,_ni,tvx,tvy,0); }
	else { if (_dnode!=id) gnorm+=v; m_map[_dnode].gni+=n0* v; }
}
// for off mesh nodes where the fields are measured at on-mesh node or element 
// fix code to always call this but only use the v and bcf if the node is off grid
// or attempt to find integated time of flight lol 
void calculate( const NodeTy _dnode,const D & _r2 , const D&  _dx, const D& _dy, const D & _ni
, const D & vxf, const D & vyf, const D bcf)
{
	if (_dnode==id) return;  // remainder term is left behind 
 	D tvxf=bcf*vxf*dt;
	D tvyf=bcf*vyf*dt;
// check convention on dx and dy doh 
	const D v= 
		m_quad.evaluate_node_pair(  -_dx,  -_dy, 0, 0, tvxf, tvyf, a, b, beta,pattern,noverlap); 
	if (debug_calculate) 
		{ MM_MSG(MMPR4(_dnode,id,_dx,_dy)<<MMPR4(_ni,vxf,vyf,bcf)<<MMPR4(tvxf,tvyf,v,gnorm)) }
	if (v<0) { clink(v,_dnode,_r2,_dx,_dy,_ni,tvxf,tvyf,bcf); }
	else { if (_dnode!=id) gnorm+=v; m_map[_dnode].gni+=n0* v; }

}
void compile()
{
// as opposed to the first version, just take whatever the difference is for
// the source remainder value 
	const D dnorm=jnorm-gnorm;
	if (dnorm>0) { m_map[id].gni=n0* dnorm; }
	else { m_map[id].gni=0; }
	jvsgnorm(); // check the calculated total and what we found on mesh traverse
}
bool jvsgnorm()
{
	const bool normx=(jnorm<gnorm);
	if (!false) if (id>=0) if (normx)
	{
		MM_INCP_MSG(mp_cm, "jnorm_lt_gnorm" )
		if (debug_norm) 
		{ MM_MSG(" the real normal is less than point total "<<MMPR3(id,jnorm,gnorm)) } 
		MM_ONCE(" the real normal is less than point total "<<MMPR3(id,jnorm,gnorm), )
	}
	MM_ONCE(" compare norma "<<MMPR3(gnorm,jnorm,id),)
return normx;

}
// if not calcualable, return n0 to change nothing
// thes user needs to copy this before next cycle 
//const D  result(const NodeTy  _dnode, const D & _r2,const int _delx, const int _dely,const D & _ni)
const int result(D & gni, const NodeTy  _dnode)
{
	int rc=2; // if nothing is done 
	if (!false) if (id<0)
	{
		// this is still larger than the original norm making a problem 
		const D bnorm=a*b*a*b*.25;
		MM_ONCE(" setting boundary node to infinite source "<<MMPR2(jnorm,bnorm),)
		jnorm=bnorm;
	}
	const D xnorm=jnorm;
	auto ii=m_map.find(_dnode);
	if (ii!=m_map.end()) { gni=(*ii).second.gni/xnorm; rc= 0; }
	else rc= 1;
	if (debug_result) { MM_MSG(MMPR4(_dnode,id,gni,gnorm)<<MMPR(rc)) }
	return rc;
}
// the original design goal was the calculate and compile  calls could be slow 
// and results fast but no reason for that. This is kind of an odd one. 
const int bulk ( D & gni)
{
int rc=2;
if (id<0)
{
	MM_INCP_MSG(mp_cm, "bulk_to_bulk" )
	MM_MSG(" dnger will robinson bulk to bulk request "<<id)
	MM_ERR(" dnger will robinson bulk to bulk request "<<id)
}
	const int dir=1;
		//  integrate the flow from a fixed distributed source such as bulk
		// boundary edge 
		// right now we do not track the point location 
	const D vplus
		=nplus*m_quad.evaluate_line( px, py,  ymax,  2 ,  tvx,  tvy, a,  b,  beta, pattern);
// evaluate_line( xpc, ypc,  xc,  dir ,  vxt,  vyt, a,  b,  beta, pattern)

	//	nplus*m_quad.fixed_line(  px,  py, dir, ymax, tvxf, tvyf, a, b, beta,pattern,noverlap); 
	const D vminus=nminus*m_quad.evaluate_line( px, py,  ymin,  0 ,  tvx,  tvy, a,  b,  beta, pattern);
		//nminus*m_quad.fixed_line(  px,  py, dir, ymin, tvxf, tvyf, a, b, beta,pattern,noverlap); 
	gni=vplus+vminus; rc=0;
return rc; 
}

void Init( const FlpTy & flp, const D & _dt,const D & _mu,const D & dx, const D & dy,CounterMap * pcm )
{
	debug_init=!true;
	debug_result=!true;
	debug_calculate=!true;
	debug_norm=!true;
	mp_cm=pcm;
	a=dx;
	b=dy;
//	static const D ps=::sqrt(M_PI);
 	bereal_cgs=1e-4*9.1e-31/4.11e-21;
	MM_ONCE(" had beb fudinging bereal_cgs from "<<bereal_cgs<<" to 1e8",)
	mu=_mu;
 	vmax=flp.vmax();
	eta=flp.fl_cutoff(); // 10e-2;
 	dt=_dt; // this is just the time of flow. ilonger more diffuse
	//const D den=mu*mu*vmax*vmax*dt*dt;
	//const D den=mu*mu*vmax*dt;
	//const D den=mu*mu*vmax*dt*dt/1e-12;
	//const D den=mu*vmax*dt*dt/1e-12;
	const D den=300.0*vmax*dt*dt/1e-12;
    //beta= bereal_cgs*mu/den; // 1e-4;   // m/(2kT) lol 
	// move to flp once the kluges are fixed 
    beta= bereal_cgs*1000.0/den; // 1e-4;   // m/(2kT) lol 
 	radius=::sqrt(-log(eta)/beta);
	// these should be cached as member vars 
	// doh - a and b are members for how 
	const IdxTy mesh_nx=flp.mesh_nx(); // 
	const IdxTy mesh_ny=flp.mesh_ny(); // 
	const IdxTy layers_y=flp.bulk_layers_y();
//	const D mesh_xmin=flp.mesh_xmin(); //
//	const D mesh_ymin=flp.mesh_ymin(); // 
ymin=flp.mesh_ymin(); // 
//	const D mesh_xmax=flp.mesh_xmax(); // 
//	const D mesh_ymax=flp.mesh_ymax(); // 
ymax=flp.mesh_ymax(); // 
//    const D a=(mesh_xmax-mesh_xmin)/(mesh_nx-1);
//    const D b=(mesh_ymax-mesh_ymin)/(mesh_ny-1);
// this must be much larger than the actual pattern rounded up 
// in case decay is faster thn a pixel size 
// to maintain consistent boundary behaviour for fixed
// niter*dt or total time elapsed, the norm needs
// to reflect loss through baths at terminals but not at
// other mesh sides. Norming over whole mess is probably the
// best thing to do. 
// this also depends on POSITION  now 
	const D ainf=a*mesh_nx; // radius*10+a*10.0;
	const D binf=b*(mesh_ny+2*layers_y); // radius*10+b*10.0;
	if ((radius>ainf)|(radius>binf)){ MM_INCP_MSG(mp_cm, "kernel_rad" )  
	if (mp_cm!=0) if (((*mp_cm).get_count("kernel_rad")%100)==1)  
			MM_MSG( "kernel_rad  ( every hundered ) "<<MMPR3(radius,ainf,binf) ) 
	} 
// so far so good, this appears to work. 
// in somes cases however it can be .001 or so too small 
 	jnorm=m_quad.normalize(a,  b,  ainf,  binf ,  beta, !true);
	pattern=flp.quad_pattern();
	noverlap=true ;// !flp.enable_overlap(); // !true;
	rejects=0;
	MM_ONCE(" make sure other init called to get overlap right ",)
//	MM_ONCE(MMPR4(radius,beta,bereal_cgs,mu)<<MMPR3(vmax,dt,m_integrate),)
if (debug_init) 	MM_MSG(MMPR4(radius,beta,bereal_cgs,mu)<<MMPR3(vmax,dt,noverlap))
}
void Init(const BraTy & tree)
{
	noverlap=!tree.enable_overlap(); // !true;
	MM_ONCE(" make sure right init has been called lol  "<<MMPR(noverlap),)
}

StrTy to_string()
{
Ss ss;
ss<<MMPR4(bereal_cgs,mu,vmax,dt);
ss<<MMPR4(beta,eta,radius,pattern);
ss<<MMPR(noverlap);
ss<<MMPR4(id,n0,vx,vy)<<MMPR4(tvx,tvy,vxt,gnorm);
ss<<MMPR(rejects);
return ss.str();
}
D bereal_cgs; // from physics at target temperature 
D mu,vmax,dt,beta,eta,radius; // vmax is scale, eta cutoff probability
NodeMap m_map;
NodeTy id;
D n0; D vx;  D vy;
D nplus,nminus,px,py; // values of infinte bulk at ymax and ymin
D ymin,ymax;
D tvx,tvy,vxt;
D gnorm,jnorm;
QuadTy m_quad; // this is a quadrature altoighm , NOT a quadrant for integration over 
D a,b ; // lattice parameters in x and y  for integrator 
IdxTy pattern;
bool noverlap;
bool debug_init,debug_result,debug_calculate,debug_norm;
IdxTy rejects;
CounterMap * mp_cm;
}; // gauss_kernel


//////////////////////////////////////////////////////////////////////////////
// 

//////////////////////////////////////////////////////////////
#if 0 
class gauss_kernel
{
/*
Calculate the flow from "point" i or src to j or x' xp or destination.
The nood itor finds the neighors of i and makes one pass first 
before making final calculations allowing for better conservation etc. 

*/
class state { public: D gni; };
typedef int NodeTy;
typedef state StateTy;
typedef std::map<NodeTy,StateTy> NodeMap;
typedef mjm_gauss QuadTy;
// evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a,b,  beta)
public:
gauss_kernel( const FlpTy & flp, const D & _dt,const D & mu, const D & dx, const D & dy, const BraTy & tree)
{
mp_cm=0;
a=dx;
b=dy;
Init(flp,_dt,mu);
Init(tree);
}
gauss_kernel( const FlpTy & flp, const D & _dt,const D & mu, const D & dx, const D & dy)
{
mp_cm=0;
a=dx;
b=dy;
Init(flp,_dt,mu);
}
~gauss_kernel()
{
if (rejects!=0) MM_MSG("bad kernel clinks " << MMPR(rejects));
}
void set(CounterMap * pcm ) { mp_cm=pcm; } 


// a new point from the out node iterator. This may be a real mesh or a boundary pixel
// carrier velocity and number are needed but may need to expand this later 
void point(const NodeTy _id, const D & _n0,const D & _vx,const D & _vy)
{
m_map.clear();
id=_id;
n0=_n0;
vx=_vx;
vy=_vy;
tvx=vx*dt;
tvy=vy*dt;
vxt=::sqrt(tvx*tvx+tvy*tvy);
gnorm=0;
if (debug_result) { MM_MSG(MMPR4(id,n0,vx,vy)<<MMPR3(tvx,tvy,vxt)) }
}
/* right now we do not need the velocity at the dest node but will later and
in the case of bulk off grid may be useful 
const D v= m_quad.evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a,b,  beta)

*/

void calculate( const NodeTy _dnode,const D & _r2 , const D& _dx, const D& _dy, const D & _ni)
{
	once_per=true;
if (!m_integrate)
{
	const D q2=(-beta*(_r2+vxt*vxt));
	const D vexp=exp(q2-2.0*beta*(tvx*_dx+tvy*_dy));
	gnorm+=vexp;
	// this map lut is going to be stuid slow 
	m_map[_dnode].gni+=(n0*vexp); 
	return;
		}
// check convention on dx and dy doh 
	//const D v= m_quad.evaluate_node_pair(  -_dx, - _dy, 0, 0, tvx, tvy, a, b, beta,pattern,noverlap&&(id!=_dnode)); 
	const D v= m_quad.evaluate_node_pair(  -_dx, - _dy, 0, 0, tvx, tvy, a, b, beta,pattern,noverlap); 
	if (debug_calculate) { MM_MSG(MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm)) }
// this needs to count these ... 	
	//if (v<0) { MM_MSG("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm)) }
	if (v<0) {++rejects; 
MM_INCP_MSG(mp_cm, "kernel_clinksc1" )

 { MM_ONCE("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm),) }  }
	else { gnorm+=v; m_map[_dnode].gni+=n0* v; }
}
// for off mesh nodes where the fields are measured at on-mesh node or element 
// fix code to always call this but only use the v and bcf if the node is off grid
// or attempt to find integated time of flight lol 
void calculate( const NodeTy _dnode,const D & _r2 , const D&  _dx, const D& _dy, const D & _ni
, const D & vxf, const D & vyf, const D bcf)
{
	once_per=true;
 D tvxf=bcf*vxf*dt;
D tvyf=bcf*vyf*dt;
if (!m_integrate) {
D vxtf2=(tvxf*tvxf+tvyf*tvyf);
	const D q2=(-beta*(_r2+vxtf2));
	const D vexp=exp(q2-2.0*beta*(tvxf*_dx+tvyf*_dy));
	gnorm+=vexp;
	// this map lut is going to be stuid slow 
	m_map[_dnode].gni+=(n0*vexp); 
	return;
	}
// check convention on dx and dy doh 
//	const D v= m_quad.evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a,b,  beta)
	//const D v= m_quad.evaluate_node_pair(  -_dx,  -_dy, 0, 0, tvxf, tvyf, a, b, beta,pattern,noverlap&&(id!=_dnode)); 
	const D v= m_quad.evaluate_node_pair(  -_dx,  -_dy, 0, 0, tvxf, tvyf, a, b, beta,pattern,noverlap); 
	//if (v<0) { MM_MSG("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm)) }
	if (debug_calculate) { MM_MSG(MMPR4(_dnode,id,_dx,_dy)<<MMPR4(_ni,vxf,vyf,bcf)<<MMPR4(tvxf,tvyf,v,gnorm)) }
	if (v<0) {
MM_INCP_MSG(mp_cm, "kernel_clinksc2" )

++rejects;  { MM_ONCE("v<0"<< MMPR4(_dnode,id,_dx,_dy)<<MMPR(_ni)<<MMPR4(tvx,tvy,v,gnorm),) }  }
	else { gnorm+=v; m_map[_dnode].gni+=n0* v; }

//	gnorm+=v;
//	m_map[_dnode].gni+=n0* v;


}
void compile()
{
// take all the values and normalize and impement conservation criteria 
// use a fixed normalization and then add back anything missing to the "from" node
// to preserve ccarrier and avoid edge effects
// gnorm is what we have, jnorm is what we should have. 
// Add anything missing back to node of origin
const D dnorm=jnorm-gnorm;
if (dnorm>0)
{
	 m_map[id].gni+=n0* dnorm; 

}
}
// if not calcualable, return n0 to change nothing
// thes user needs to copy this before next cycle 
//const D  result(const NodeTy  _dnode, const D & _r2,const int _delx, const int _dely,const D & _ni)
#endif
#if 0
const int result(D & gni, const NodeTy  _dnode)
{
int rc=2; // if nothing is done 
//D gni=n0; // ni?
if (jnorm<gnorm)
{
MM_INCP_MSG(mp_cm, "jnorm_lt_gnorm" )
MM_MSG(" the real normal is less than point total "<<MMPR2(jnorm,gnorm))
}
MM_ONCE(" compare norma "<<MMPR3(fnorm,gnorm,jnorm),)
if (id<0)
{
const D bnorm=a*b;
MM_ONCE(" boundary norms are ignored per latest notes "<<MMPR3(jnorm,gnorm,bnorm),)
jnorm=bnorm;
gnorm=bnorm;
}
//if ( once_per ) { MM_MSG(" compare norma "<<MMPR3(fnorm,gnorm,jnorm)) once_per=false; } 
if (gnorm!=0) 
{
	const D xnorm=real_normal?jnorm:gnorm;
	auto ii=m_map.find(_dnode);
	if (ii!=m_map.end()) { gni=(*ii).second.gni/xnorm; rc= 0; }
	else rc= 1;
}
if (debug_result) { MM_MSG(MMPR4(_dnode,id,gni,gnorm)<<MMPR2(fnorm,rc)) }
return rc;
}

void Init( const FlpTy & flp, const D & _dt,const D & _mu)
{
	debug_init=!true;
	debug_result=!true;
	debug_calculate=!true;
	// this has a negative bias due to cutoff radius 
	real_normal=true;
//	static const D ps=::sqrt(M_PI);
 	bereal_cgs=1e-4*9.1e-31/4.11e-21;
	MM_ONCE(" had beb fudinging bereal_cgs from "<<bereal_cgs<<" to 1e8",)
 	//bereal_cgs=1e8;
	mu=_mu;
 	vmax=flp.vmax();
 	dt=_dt; // this is just the time of flow. ilonger more diffuse
	//const D den=2000*vmax*vmax*dt*dt;
	const D den=mu*mu*vmax*vmax*dt*dt;
	// this is ficking backwards lol 
    //beta= bereal_cgs*mu/den; // 1e-4;   // m/(2kT) lol 
    beta= bereal_cgs*1000.0*1000.0/den; // 1e-4;   // m/(2kT) lol 
	eta=10e-2;
 	radius=::sqrt(-log(eta)/beta);
	// kluge for cutoff, probably neds to be to power or something but so discrete 
	// likelyt to make a mess; This should be sqrt(pi/b) SQUARED for 2D but
	// still off by power of 2? 
	// doh, this ignores shape functions which can be a big deal	
	fnorm=M_PI/beta; // beta/M_PI; // /(1.0-eta);
	fnorm=400.0*1.0/(1.0*beta*beta); // beta/M_PI; // /(1.0-eta);
	const D infs=radius*100; // 10 worked here, 1 is way too small
	const D area=D(flp.mesh_nx()*flp.mesh_ny());
	// considering shape functions it is amazing this is close  lol. 
//	static const D fudg= ::sqrt(M_PI);
//	jnorm= fudg*m_quad.evaluate_node_pair(0,0,0,0, 0, 0, infs, infs, beta,0,false)/area; 
#endif
#if 0
const IdxTy mesh_nx=flp.mesh_nx(); // 
const IdxTy mesh_ny=flp.mesh_ny(); // 
const D mesh_xmin=flp.mesh_xmin(); //
const D mesh_ymin=flp.mesh_ymin(); // 
const D mesh_xmax=flp.mesh_xmax(); // 
const D mesh_ymax=flp.mesh_ymax(); // 
    const D a=(mesh_xmax-mesh_xmin)/(mesh_nx-1);
    const D b=(mesh_ymax-mesh_ymin)/(mesh_ny-1);
// this must be much larger than the actual pattern rounded up 
// in case decay is faster thn a pixel size 
	const D ainf=radius*10+a*10.0;
	const D binf=radius*10+b*10.0;
// so far so good, this appears to work. 
 	jnorm=m_quad.normalize(a,  b,  ainf,  binf ,  beta, !true);

	once_per=true;
	m_integrate=true;
	pattern=flp.quad_pattern();
	noverlap=true ;// !flp.enable_overlap(); // !true;
	rejects=0;
	MM_ONCE(" make sure other init called ",)
//	MM_ONCE(MMPR4(radius,beta,bereal_cgs,mu)<<MMPR3(vmax,dt,m_integrate),)
if (debug_init) 	MM_MSG(MMPR4(radius,beta,bereal_cgs,mu)<<MMPR4(vmax,dt,m_integrate,noverlap)<<MMPR(real_normal))
}
void Init(const BraTy & tree)
{
	noverlap=true; // !tree.enable_overlap(); // !true;
	noverlap=!tree.enable_overlap(); // !true;
	MM_ONCE(" make sure right init has been called lol  ",)
}

StrTy to_string()
{
Ss ss;
ss<<MMPR4(bereal_cgs,mu,vmax,dt);
ss<<MMPR4(beta,eta,radius,pattern);
ss<<MMPR2(m_integrate,noverlap);
ss<<MMPR4(id,n0,vx,vy)<<MMPR4(tvx,tvy,vxt,gnorm);
ss<<MMPR(rejects);
return ss.str();
}
D bereal_cgs; // from physics at target temperature 
D mu,vmax,dt,beta,eta,radius; // vmax is scale, eta cutoff probability
NodeMap m_map;
NodeTy id;
D n0; D vx;  D vy;
D tvx,tvy,vxt;
D gnorm,fnorm,jnorm;
QuadTy m_quad; // this is a quadrature altoighm , NOT a quadrant for integration over 
D a,b ; // lattice parameters in x and y  for integrator 
bool m_integrate;
IdxTy pattern;
bool noverlap, real_normal;
bool debug_init,debug_result,debug_calculate,once_per;
IdxTy rejects;
CounterMap * mp_cm;
}; // gauss_kernel
#endif


//////////////////////////////////////////////////////////////////////////////
// 
/*
Destination functors. These right now either peform the flows on an update
matrix to be sent back to m_state or matrix eleemnts to find the
steady state solution self consistently with u. 

*/
#endif

