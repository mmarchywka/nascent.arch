#ifndef MJM_FICK_LANCZOS_2D_H__
#define MJM_FICK_LANCZOS_2D_H__

#include "mjm_fick_equations.h"
#include "mjm_globals.h"
//#include "mjm_qfle.h"
#include "mjm_logic_base.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_diffuse_iterators.h"
#include "mjm_defect_levels.h"


#include <algorithm>
#include <map>
#include <signal.h>
#include <cstdlib>
#define USING_PETSC
#ifdef USING_PETSC
// deviate from shooting to matrix 
#include "mjm_petsc_util.h"
#include "mjm_petsc_fd_base.h"
#endif

// libmesh crap 
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mjm_libmesh.h"
// fudd 
#pragma GCC diagnostic ignored "-Wunused-parameter"




/*
2017-06-22 copy from mjm_shooting11.h and remove all the 1D crap debate about adding libmesh

2017-06-25 skeleton flmb runs miserably . Node location calculation
is wrong, see notbook for element numbering. Preint itor info to get some 
idea wtf
 put real code here void update_drift(MyBlock & state)

2017-06-17 make velocity approximation by recgangle or element
and take an incorrect average between src and dest. Fixed some 
things but did make low point solutions oscillate and left features
in depletion region but currents may be better at high point
density.

2017-05-25
the uniform refine urefine script seems to work but there is
ringing on the left. Suspect either shooting init error or
v(x) problem. 
TODO

The non-uniform refinement seems to produce problems at
trnasition regions that vary with the uniform_src or 
node_only setting. Hopefully this points to a code error
somewhere that fixes all. 

The BC ringing was likely due to retaining V0 well into the neutral
regions. This should be zero and the shooter needs to reflect
zero field in the bulk although maybe a contradiction with
current continuity the return current is purely diffusion. 

1) left ringing - persists after right shoot but petsc explodes
2) refine traps  - coded no idea if it works


2017-04-29 The MB maxwell-boltzmann integrals seem to work
but getting the flow and pnj not yet done. there is terrible
numerical noise as T-> ininity ( beta goes to zero as things
spread out ). 
2017-04-16 100 pt example pn abrupt looks pretty good but
there are no oscillations on the left wall of depletion edge. 
Currents not equal but the caclulation is coarse too

TODO check the drift code for very small and zero velocity. 
Overlap diffusion is a big issue for real accuracy. 
Make the code  subtracting zero field "drift" optional. 


*/

class fl_params : public mjm_logic_base 
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
//FickLancParamsVar( const StrTy & nm) : m_map(nm) {}
fl_params( const StrTy & nm) : Super(nm) {}
fl_params() : Super() {}
// these need to be GENERATED for consistency 
//IdxTy output_freq() const { return m_map.get_uint("output_freq",100); } // // 100;
//D total_length() const { return m_map.get_double("total_length",.25e-6*2000); }
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
D time_step() const { return m_map.get_double("time_step",1e-13); }
IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",100); } // // 100;
IdxTy mesh_ny() const { return m_map.get_uint("mesh_ny",100); } // // 100;
// cgs mostly 
D mesh_xmin() const { return m_map.get_double("mesh_xmin",0); }
D mesh_ymin() const { return m_map.get_double("mesh_ymin",0); }
D mesh_xmax() const { return m_map.get_double("mesh_xmax",1e-4); }
D mesh_ymax() const { return m_map.get_double("mesh_ymax",1e-4); }
IdxTy max_iterations() const { return m_map.get_uint("max_iterations",10000); }
IdxTy output_freq() const { return m_map.get_uint("output_freq",100); } // // 100;
IdxTy output_offset() const { return m_map.get_uint("output_offset",output_freq()-1); } // freq-1;
StrTy output_label() const { return m_map.get_string("output_label","fick"); }

StrTy dump_file() 
	const { return m_map.get_string("dump_file","fl_default_dump.txt"); }
StrTy dump_label()
	const { return m_map.get_string("dump_label","fl_default_label"); }
IdxTy dump_flags()
	const { return m_map.get_uint("dump_flags",0); }


D gr_xsection() const { return m_map.get_double("gr_xsection",1e-12); } // 1e-3;
D gr_tau() const { return m_map.get_double("gr_tau",1e-8); }


D vsat() const { return m_map.get_double("vsat",1e7); } // 1e7 for Si

D Vt() const { return m_map.get_double("Vt",.0259); }
D qe() const { return m_map.get_double("qe",1.6e-19/8.854e-14); }
D ni2() const { return m_map.get_double("ni2",1.5e10*1.5e10); }
D mun() const { return m_map.get_double("mun",1000); }
D mup() const { return m_map.get_double("mun",100); }
D muz() const { return m_map.get_double("mun",10); }
D f() const { return m_map.get_double("f",.5); }
D Nd() const { return m_map.get_double("Nd",f()*1e15); }
D Na() const { return m_map.get_double("Na",f()*1e14); }
D v0() const { return m_map.get_double("v0",0); }
D vd() const { return m_map.get_double("vd",.45); }
D trap_speed() const { return m_map.get_double("trap_speed",1.0); }
D trap_N() const { return m_map.get_double("trap_N",1e13); }
D vmax() const { return m_map.get_double("vmax",1e0); }




StrTy to_string(const StrTy & sep=" ") const
{
// should be geneated code, do not edit  wtf ?
Ss ss;

ss<<"time_step="<<time_step()<<sep;
ss<<"mesh_nx="<<mesh_nx()<<sep;
ss<<"mesh_ny="<<mesh_ny()<<sep;
ss<<"mesh_xmin="<<mesh_xmin()<<sep;
ss<<"mesh_ymin="<<mesh_ymin()<<sep;
ss<<"mesh_xmax="<<mesh_xmax()<<sep;
ss<<"mesh_ymax="<<mesh_ymax()<<sep;
ss<<"max_iterations="<<max_iterations()<<sep;
ss<<"output_freq="<<output_freq()<<sep;
ss<<"output_offset="<<output_offset()<<sep;
ss<<"output_label="<<output_label()<<sep;
// maybe a macro would work? 
ss<<"dump_file="<<dump_file()<<sep;
ss<<"dump_label="<<dump_label()<<sep;
ss<<"dump_flags="<<dump_flags()<<sep;

ss<<"gr_xsection="<<gr_xsection()<<sep;
ss<<"gr_tau="<<gr_tau()<<sep;
ss<<"vast="<<vsat()<<sep;
ss<<"Vt="<<Vt()<<sep;
ss<<"qe="<<qe()<<sep;
ss<<"ni2="<<ni2()<<sep;
ss<<"mun="<<mun()<<sep;
ss<<"mup="<<mup()<<sep;
ss<<"muz="<<muz()<<sep;
ss<<"f="<<f()<<sep;
ss<<"Nd="<<Nd()<<sep;
ss<<"Na="<<Na()<<sep;
ss<<"v0="<<v0()<<sep;
ss<<"vd="<<vd()<<sep;
ss<<"trap_speed="<<trap_speed()<<sep;
ss<<"trap_N="<<trap_N()<<sep;
ss<<"vmax="<<vmax()<<sep;
return ss.str();
}

}; // fl_params

class fl_branches : public mjm_logic_base 
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
//FickLancParamsVar( const StrTy & nm) : m_map(nm) {}
fl_branches( const StrTy & nm) : Super(nm) {}
fl_branches() : Super() {}
// these need to be GENERATED for consistency 
//bool show_iter_status() const { return m_map.get_bool("show_iter_status",true); }
//bool always_dump_state() const { return m_map.get_bool("always_dump_state",true); }
bool show_iter_status() const { return m_map.get_bool("show_iter_status",true); }
bool always_dump_state() const { return m_map.get_bool("always_dump_state",!true); }
bool always_dump_to_file() const { return m_map.get_bool("always_dump_to_file",!true); }
bool dump_state_everyfew() const { return m_map.get_bool("dump_state_everyfew",true); }
bool dump_state_atend() const { return m_map.get_bool("dump_state_atend",true); }
bool dump_counter_everyfew() const { return m_map.get_bool("dump_counter_everyfew",true); }
bool cerr_counter_everyfew() const { return m_map.get_bool("cerr_counter_everyfew",true); }
bool enable_trapping() const { return m_map.get_bool("enable_trapping",!true); }
bool enable_gr() const { return m_map.get_bool("enable_gr",true); }
bool enable_flmb() const { return m_map.get_bool("enable_flmb",true); }
bool enable_dump4() const { return m_map.get_bool("enable_dump4",true); }

StrTy to_string(const StrTy & sep=" ") const
{
// should be geneated code, do not edit  wtf ?
Ss ss;
ss<<"show_iter_status="<<show_iter_status()<<sep;
ss<<"always_dump_state="<<always_dump_state()<<sep;
ss<<"always_dump_to_file"<<always_dump_to_file()<<sep;
ss<<"dump_state_everyfew="<<dump_state_everyfew()<<sep;
ss<<"dump_state_atend="<<dump_state_atend()<<sep;
ss<<"dump_counter_everyfew="<<dump_counter_everyfew()<<sep;
ss<<"cerr_counter_everyfew="<<cerr_counter_everyfew()<<sep;
ss<<"enable_trapping="<<enable_trapping()<<sep;
ss<<"enable_gr="<<enable_gr()<<sep;
ss<<"enable_flmb="<<enable_flmb()<<sep;
ss<<"enable_dump4="<<enable_dump4()<<sep;

//ss<<"vmax="<<vmax()<<sep;
return ss.str();
}

}; // fl_params
//////////////////////////////////////////////////////////////////

#if 0 
#endif

////////////////////////////////////////////////////////////////////////
//
//  Actually return the right nodes for neighbors using a LUT
// between nodes numbered on a grid and actual mesh node number


class mesh_indexing_new 
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<int> MyIntBlock;
typedef mjm_block_matrix<IdxTy> MyuIntBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 
typedef mesh_indexing_new Myt;
public:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MyIntBlock  MyIntBlock;
typedef Tr::MyuIntBlock  MyuIntBlock;
typedef Tr::MySparse MySparse;
typedef int intt;

typedef fl_params FlpTy;
typedef libMesh::SerialMesh MyMesh;
typedef libMesh::Node Node;
//typedef std::vector<Node> node_list;
typedef std::vector<IdxTy> node_list;
typedef libMesh::MeshBase::node_iterator NodeItor;
typedef std::map<libMesh::dof_id_type,Real> DiriMap;
// index is the geometric id, x+nx*y , with value of the libmesh node id 
// this is now SIGNED as the negative values are peripheral nodes. 
typedef MyIntBlock NodeIndex;
// right now this is a tic tack intero not round 
class round_itor
{
public:

typedef unsigned int IdxTy;
typedef double D;
round_itor(): m_ok(false) {}
// nx and ny already are padded with 2*padding 
round_itor(const D & distance, const IdxTy _nx, const IdxTy _ny, const D & _dx,
const D & _dy, const IdxTy padding)
: d(distance),nx(_nx),ny(_ny),dx(_dx),dy(_dy), layers(padding) {Init(); } 
// these are in mesh corrdinates
void center(const int cx, const int cy) { cnx=cx; cny=cy; } 
bool ok() { return m_ok; } 
void inc()
{
if (!ok()) return;
++posx;
if ( posx>nmaxx) 
{ 
	posx=nminx; 
	if (posy==nminy) { m_ok=false; }
	else --posy; // unsigned fudd again 
}
}
// posx is relative to total grid, NOT just real grod
int relx() const { return int(posx)-int(cnxl); } 
int rely() const { return int(posy)-int(cnyl); } 
D delx() const { return dx*relx(); } 
D dely() const { return dy*rely(); } 
int meshx() const { return int(posx)-layers; }
int meshy() const { return int(posy)-layers; }
D r2() const { return relx()*relx()*dx*dx+rely()*rely()*dy*dy; } 
// these are all medh parameters and this will return the right geo node
// for on mesh nodes. However oor makes this wrong 
//int node() const { return (posx-m_layers)+nx*(posy-m_layers); } 
// this is a geo node starting with zero at padding  lower left.
//int node() const { return (posx)+nx*(posy); } 

void start() {
// zero on the mesh is layers,layers on the composite grid 
cnxl=cnx+layers;
cnyl=cny+layers;
// these are already padded with 2layers
nxl=nx; //+layers;
nyl=ny;//+layers;
// this is signed??? no
// these limits need to extend by layers on each edge
intt xx=intt(d/dx+.5)+1;

//if (xx>cnxl) nminx=-layers; else nminx=cnx-xx;
if (xx>cnxl) nminx=0; else nminx=cnxl-xx;
nmaxx=cnxl+xx;
if (nmaxx>= nxl ) nmaxx=nxl-1;

intt yy=intt(d/dy+.5)+1;
//if (yy>cnyl) nminy=-layers; else nminy=cny-yy;
if (yy>cnyl) nminy=0; else nminy=cnyl-yy;
nmaxy=cnyl+yy;
if (nmaxy>= nyl ) nmaxy=nyl-1;

posx=nminx;
posy=nmaxy;
// it is possible that there are no elements in the range, hence !true
// center could be invalid, handle run off in boundary code
m_ok=(nminx<=nmaxx)&&(nminy<=nmaxy); // true; // minx<=maxx, miny<=maxy
}
void Init() { m_ok=false; } 
StrTy to_string() const
{
Ss ss;
ss<<MMPR4(nx,ny,cnx,cny);
ss<<MMPR2(cnxl,cnyl);
ss<<MMPR4(nminx,nminy,nmaxx,nmaxy);
ss<<MMPR3(dx,dy,d);
ss<<MMPR3(posx,posy,m_ok);

return ss.str();

}
int nx,ny,nminx,nminy,nmaxx,nmaxy;
int cnx,cny;
int cnxl,cnyl,nxl,nyl;
D dx,dy,d;
int posx,posy;
bool m_ok;
int layers;
}; // round_itor

typedef round_itor NayItor;
class peripheral_nodes
{

public:


}; // peripheral_nodes

mesh_indexing_new(MyMesh & m, const FlpTy & p, const fl_branches & b) 
	: m_mesh(&m), m_flp(p),m_branches(b) {}
mesh_indexing_new() : m_mesh(0), m_flp(),m_branches() {}
// these are now SIGNED , try to get non sign extend 
int  bad() const { return ~ int((~0U)>>1); } 

// neighbors
// these take a mesh node id, translate to geo location, and then back
// to a node id 
int north( const IdxTy node ) const { return neighbormesh(node,0,1); }
int south( const IdxTy node ) const { return neighbormesh(node,0,-1); }
int east( const IdxTy node ) const { return neighbormesh(node,1,0); }
int west( const IdxTy node ) const { return neighbormesh(node,-1,0); }

bool off_mesh_x(const int n)
{
if (n>=m_n_x) return true;
if (n<0) return true;
return false;
}
bool off_mesh_y(const int n)
{
if (n>=m_n_y) return true;
if (n<0) return true;
return false;
}

int neighbor( const int  nodeid, const int dx, const int dy ) const
{
	const int node= nodeid;
	int x0=((node%m_n_x));
	int x=(x0+dx);
	if (x>=m_n_xl) return bad();
	if (x< -m_layers) return bad();
	int y0=((node/m_n_x));
	int y=(y0+dy);
	if (y>=m_n_yl) return bad();
	if (y<-m_layers) return bad();
	//return x+m_n_x*y;
	return index(x,y);
}
// this does not quote work on first rows due to element forming 
// FALSE turnes a node id into a geo node 
int neighbormesh( const int  nodeid, const int dx, const int dy ) const
{ 	const bool  fr=false; 
	int node=lut(nodeid,fr);
	int x0=((node%m_n_x));
	int x=(x0+dx);
	if (x>=m_n_x) return bad();
	if (x< 0) return bad();
	int y0=((node/m_n_x));
	int y=(y0+dy);
	if (y>=m_n_y) return bad();
	if (y<0) return bad();
	//MM_MSG(MMPR4(x0,x,y0,y)<<MMPR4(m_n_x,m_n_y,nodeid,node))
	//return x+m_n_x*y;
	// this needs to return an on-grid index
	// geo location id becomes a libmesh node id 
	return lut(x+ m_n_x*y,!fr); // index(x,y);
}

// return the approximate node number offset by dx and dy from node node
IdxTy node_from( const IdxTy node, const D & dx, const D & dy ) const
{
	const bool xneg=(dx<0);
	const bool yneg=(dy<0);

	D dnx=((xneg)?(-dx):dx)/m_dx;
	D dny=((yneg)?(-dy):dy)/m_dy;
	intt dnxi=intt(dnx+.5);
	intt dnyi=intt(dny+.5);
	return node+(xneg?(-dnxi):dnxi)+m_n_x*((yneg)?(-dnyi):dnyi);	
}
typedef round_itor hood_iterator;
// this is an on grid geo node center relative to grid ex padding.
// when true, bc means that the node is off grid and encoded on larger mesh 
hood_iterator hood_itor(const D & r, const int node, const D & offx, const D & offy, const bool bc=false)
{
const int sx=(offx>0)?1:(-1);
const int sy=(offy>0)?1:(-1);
// these assume node is on grid 
int ll=m_n_x;
if (bc ) ll+=2*m_layers;
int cx=(node%ll)+ int(sx*offx/m_dx+.5)*sx;
int cy=(node/ll)+ int(sy*offy/m_dy+.5)*sy;
// now put back on grid if needed
if (bc ) { cx-=m_layers; cy-=m_layers;}
// the iterator however should use full sizes
//round_itor ri(r, m_n_x, m_n_y, m_dx,m_dy,m_layers);
round_itor ri(r, m_n_x+2*m_layers, m_n_y+2*m_layers, m_dx,m_dy,m_layers);
ri.center( cx,cy);
ri.start();
return ri;

}

void make_index()
{



} // mahke_index

// setup variables for a regular mesh 
void load_regulars() { load_regulars(m_flp); } 

void load_regulars(const FlpTy & flp)
{
	m_n_x=flp.mesh_nx(); // 
	m_n_y=flp.mesh_ny(); // 
	const D mesh_xmin=flp.mesh_xmin(); //
	const D mesh_ymin=flp.mesh_ymin(); // 
	const D mesh_xmax=flp.mesh_xmax(); // 
	const D mesh_ymax=flp.mesh_ymax(); // 
	m_minx=mesh_xmin;
	m_miny=mesh_ymin;
	m_dx=(mesh_xmax-mesh_xmin)/(m_n_x-1);
	m_dy=(mesh_ymax-mesh_ymin)/(m_n_y-1);
} // load_regulars 
// index nodes is a hige reduncant thing vs the m_fixed
// FALSE takes x as libmesh nodeid and creates a geo location x+mesh_n_x_*y
// right now this only works on grid, oor return a negatice for now 
int lut(const int x, const bool fwd) const 
{
const IdxTy dir=fwd?0:1;
if (x<0) return x;
// these it no const accessor for this thing... 
int y= m_node_idx(x,dir);
//MM_MSG(MMPR3(x,y,dir))
//if (y==0) {  MM_MSG(" maps to zero "<<MMPR3(x,dir,y))}
return y;
}
void dump_lut(const StrTy & label)
{
std::ostream & os=std::cout;
const IdxTy sz=m_node_idx.size(0);
for(IdxTy i=0; i<sz; ++i)
{
MM_MSG(MMPR4(label,i,m_node_idx(i,0),m_node_idx(i,1)))

}

}


void index_nodes()
{
//	m_node_idx.resize(nx*ny-m_index_base,3);
// either use element iterator and assumed traversal and node order
// or guess the location based on x and y and dx and dy.... or make
// grid based on rational numbers lol.
// it is very slow but just do everything and then fix the mesh entries
IdxTy serial=0;
IdxTy ms=m_node_idx.dim(0); 
//for (IdxTy i=-m_layers; i<m_n_yl; ++i)
for (IdxTy i=0; i<m_n_y; ++i)
{
	const int dy=m_n_x*i;
	//for (IdxTy j=-m_layers; j<m_n_xl; ++j)
	for (IdxTy j=0; j<m_n_x; ++j)
	{
		const int luti=j+dy; // index(j,i);
		//m_node_idx(luti-m_index_base,0)=serial;	
		if ((luti<0)||(luti>ms)) { MM_MSG(" out of range "<<MMPR4(luti,ms,i,j))}
		if ((serial<0)||(serial>ms)) { MM_MSG(" out of range "<<MMPR4(serial,ms,i,j))}
		m_node_idx(luti,0)=serial;	
		m_node_idx(serial,1)=luti;	
		++serial;
	} // j
} // i

	const NodeItor ie=(*m_mesh).nodes_end();
	for (NodeItor ii=(*m_mesh).nodes_begin(); ii!=ie; ++ii)
	{
		const Node & np=*(*ii);
		const IdxTy id=np.id();
		const IdxTy nx=IdxTy((np(0)-m_minx)/m_dx+.5);
		const IdxTy ny=IdxTy((np(1)-m_miny)/m_dy+.5);
		//m_node_idx(index(nx,ny)-m_index_base,0)=id;	
		m_node_idx(index(nx,ny),0)=id;	
		m_node_idx(id,1)=index(nx,ny);	
	} // ii 
 	// the oor values all map to zero right now, add the others 
	// on second thought, do all of them first 
}

// this needs to deal with negative values... 
// this takes a geo location and returns the geo id.
int index(const int nx, const int ny)  const
{ 
	int rc=bad();
	const IdxTy ooc=oor(nx,ny);
	switch (ooc)
	{
	case 0 : { rc= nx+m_n_x*ny;  break; }
	// below is easy 
	case 4:  // nx<0, return a negative index that is index into peri 
	case 5:  // nx<0, return a negative index that is index into peri 
	case 6:  // nx<0, return a negative index that is index into peri 
			// m_index_base<0 
		//	{rc= -m_index_base-( 1+ nx + ny*m_peri_x);	break; } 

	// above is easy 
	case 8:  // nx<0, return a negative index that is index into peri 
	case 9:  // nx<0, return a negative index that is index into peri 
	case 10:  // nx<0, return a negative index that is index into peri 
			// m_index_base<0 
		//	{ rc= -( 1+m_peri_x - nx+(m_peri_y-ny-1)*m_peri_x);break; }	
	// for y in range, two x cases
	case 1: 
// {rc= m_peri_base+(m_layers-nx)+(ny-m_layers)*2*m_layers; break; } 
	case 2: 
	//{ rc=m_peri_base+(m_layers-nx-m_n_x)+(ny-m_layers)*2*m_layers;  break; }
	const int xzed=nx+m_layers;
	const int yzed=ny+m_layers;
	const int raw= xzed+yzed*(m_n_x+2*m_layers);
	if (ny<0) {rc=-1-raw;}
	else if (ny>=m_n_y) { rc=-1-raw+m_n_x*m_n_y; }
	else { rc=-1-raw+m_n_x*(ny); if (nx>0) rc+=m_n_x;}
	}
if (rc==bad()) { MM_ERR(" bad value "<<MMPR4(ooc,nx,ny,m_n_y)) }
if (ooc!=0) if (rc>0) { MM_ERR(" bad value  for oor rc > 0 "<<MMPR(rc)<<MMPR4(ooc,nx,ny,m_n_y)) }
//if (rc==0) { MM_MSG(" map to zero "<<MMPR4(nx,ny,ooc,rc))} 
return rc;
}
int oor(const int nx, const int ny)  const
{
int rc=0;
if (nx<0) rc|=1;
if (nx>=m_n_x) rc|=2;
if (ny<0) rc|=4;
if (ny>=m_n_y) rc|=8;

return rc;
}
// set up the indexing to include layers of nodes 
// around a core of nx x ny nodes
int make_peripheral_nodes(const IdxTy nx, const IdxTy ny, const IdxTy layers,
const IdxTy nstate, const IdxTy nfixed)
{
	m_layers=layers;
	const IdxTy l2=layers*2;
	// 1 -> 2(nx+ny) + 4, 2-> 4(nx+ny)+16
	const int size=(nx+l2)*(ny+l2)-nx*ny; //l2*((nx+ny)+ l2);
	m_index_base=-size;
	m_peri_x=m_n_x+l2;
	m_peri_y=m_n_y+l2;
	m_n_xl=m_n_x+m_layers;
	m_n_yl=m_n_y+m_layers;
	m_peri_base=m_peri_x*m_layers;
	m_peri_state.resize(size,nstate);
	m_peri_fixed.resize(size,nfixed);
	// this really only needs to be the size of the real grid. 
	m_node_idx.resize(nx*ny+size,3);
	index_nodes();
}

	void set_positions(const IdxTy ix,const IdxTy iy,const D & mesh_xmin
		,const D & mesh_ymin,const D & dx,const D & dy)
{
typedef std::map<int,std::vector<int> > HitMap;
HitMap hm;
const IdxTy sz=m_peri_fixed.dim(0);
for (int i=-m_layers; i<m_n_yl; ++i)
{
const bool inmesh=(i>=0)&&(i<m_n_y);
for (int j=-m_layers; j<m_n_xl; ++j)
{
if (inmesh) if (j==0) j=m_n_x; 
int idx=index(j,i);
//const int eff=idx-m_index_base;
const int eff=-1-idx; // -m_index_base;
if (hm.find(eff)==hm.end())
	{ std::vector<int> loc; loc.push_back(i); loc.push_back(j);  hm[eff]=loc; } 
else 
{
MM_MSG(" duplicate index "<<MMPR3(i,j,eff)<<" vs "<<MMPR3(hm[eff][0],hm[eff][1],eff))
}
if ((eff<0)||(eff>=sz)) {MM_MSG(" bad pos "<<MMPR4(i,j,eff,sz)<<MMPR2(m_index_base,idx))}
m_peri_fixed(eff,ix)=mesh_xmin+dx*j;
m_peri_fixed(eff,iy)=mesh_ymin+dy*i;
} // j 

} // i 



} // set_positions


MyMesh * m_mesh;
FlpTy m_flp;
fl_branches m_branches;
D m_dx,m_dy,m_minx,m_miny;
int m_n_x, m_n_y;
int  m_n_xl, m_n_yl;
NayItor m_ri;
NodeIndex m_node_idx;
int m_layers;
int m_index_base,m_peri_x,m_peri_y,m_peri_base; // minimum value of the peripheral pixels
MyBlock m_peri_state; // same as m_state except it is only the peripheral values that do not come from solver 
MyBlock m_peri_fixed; 
}; // mesh_indexing_new




////////////////////////////////////////////////////////////////

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
static void interrupt()  { raise(SIGINT); } 
// this needs a way to keep track of levels and only enable
// soo deep 
int m_level;
bool m_enable;
const Myt * m_mom;
}; //debug_levels

#if 0 

class debug_and_diagnose_status
{
/*
Maintain a list of interesting points in the simulation
for debug dumps or other output and examination. May be used
for mesh refinement at some point too. 
*/
// map of node index to debug flags. 
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::map<IdxTy, IdxTy> Flags;
typedef std::map<IdxTy, StrTy> Bits;
typedef std::map<IdxTy, StrTy> Stats;

public:

// reset to ctor state
void clear() {}
// set up for next iteration, may just call clear
void iteration() {}
void flag(const IdxTy idx, const IdxTy bits)
{ m_flagged[idx]|=bits; }
bool flagged(const IdxTy idx, const IdxTy bits=~0)
{ return (m_flagged[idx]&bits)!=0; }
void log(const IdxTy idx, const StrTy & lbl , const D & val)
{


}
 
const IdxTy &  flags(const IdxTy idx) { return m_flagged[idx]; } 
void dump() 
{
typedef Stats::const_iterator Ci;
for (Ci ii=m_stats.begin(); ii!=m_stats.end(); ++ii)
{
const IdxTy idx=(*ii).first;

MM_MSG(idx<<" "<<m_flagged[idx]<<" "<<(*ii).second)

} 


}
// dynamically flagged 
Flags m_flagged;
// always flagged
Flags m_defaults;
// create names for the bits 
Bits m_bits;
// text description of interesting points 
Stats m_stats;
}; // debug_and_diagnose_status


////////////////////////////////////////////////////////////////

#endif


class mjm_fick_lanczos_2d // petsc_log_fd :public petsc_fd_base //  public Petnl::jr_base
{
/*
Simulate fick's law with various velocity profiles. This 
leads to lancozs derivatives and appears to give
results that look right but quantitative comparisons lacking. 
The far FD and FE code was removed but easy to add back.
Right now with linar g(v) the abrupt PNJ gives
good results on rho with 1 element depletion region
but lightly doped ( lsmaller slope ) side wall
has oscillations. lol.  



*/


class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef mjm_fick_lanczos_2d Myt;
//typedef petsc_fd_base Super;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;

//typedef debug_and_diagnose_status Dads;

//	typedef  quadratic_fick_lanc_element QFLE;
//typedef QFLE::drift_flags drift_flags;
typedef mjm_defect_levels Defects;
//typedef quadratic_fick_lanc_element_sf QFLESF;
//	typedef  quadratic_fick_lanc_element_sf QFLE;
//	typedef diffuse_iterator Di;

// this makes no real sense now but see GR code 
//	typedef diffuse_iterator_sf Di;

typedef  mjm_fick_es EsTy; // : public EquationSystems
typedef  mjm_fick_jr NonLinTy; // : public NonlinearImplicitSystem,


typedef libMesh::LibMeshInit InitTy;
typedef libMesh::SerialMesh MyMesh;
typedef libMesh::MeshBase::node_iterator NodeItor;
typedef std::map<libMesh::dof_id_type,Real> DiriMap;
typedef std::map<int,D> BCMap;

//typedef libMesh::QUAD4 ElemTy;
libMesh::ElemType ElemTy=libMesh::ElemType::QUAD4;
typedef mesh_indexing_new MeIdx; 
#if 0
 MyMesh mesh(init.comm());
  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
  MeshTools::Generation::build_square (mesh,
                       //                15, 15,
                                       121, 25,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD9);

  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  mesh.print_info();

#endif


public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
mjm_fick_lanczos_2d():m_old(0),m_points(3),m_save_sol(false),m_exit(false),m_init(0),m_mesh(0),m_es(0),m_jr(0) {Init();}
// make this m_points now, wtf
mjm_fick_lanczos_2d(const IdxTy s):m_old(0), m_points(s), m_save_sol(false),m_exit(false),m_init(0),m_mesh(0),m_es(0),m_jr(0) {Init();}
mjm_fick_lanczos_2d(int argc,char **args, InitTy & init) : m_old(0),m_save_sol(false),m_exit(false),m_init(0),m_mesh(0), m_es(0),m_jr(0)
{
m_init=&init;
m_size=1;
m_points=2000; // this really should not have a default.. 
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
const StrTy s=StrTy(args[i]);
if (s==StrTy("-interrupt")) debug_levels::interrupt(); 
m_tree.config("-tree",i,argc,args);
m_flp.config("-params",i,argc,args);
configi(m_points,"-points",i,argc,args);
cmdcmd("-cmd",i,argc,args);
source("-source",i,argc,args);
m_flp.config_set("-set-param",  i,  argc, args);
m_tree.config_set("-set-branch",  i,  argc, args);
if (i==istart) { MM_ERR(" did nothing with i="<<i<< " " <<args[i]) {++i; } } 

}
if (!m_exit) Init();
}
~mjm_fick_lanczos_2d() { delete m_old;  delete m_mesh; delete m_es; 
// owned vy es afaict 
//delete m_jr;
} 
const bool exit() const { return m_exit; } 
////////////////////////////////////////////////////////
// command block

// this should be in the parameters map, nothing special here... 
 void cmdcmd( const StrTy & cmd, int  & i, int argc, char ** args)
{
// this is missing para 
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
if (argc<=i) return; 
const StrTy nm=StrTy(args[i]);
MM_ERR(" executing  "<<nm<<"  ")
command_mode(nm);
++i; // consume param and cmd
}
}

 void source( const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
if (argc<=i) return; 
const StrTy nm=StrTy(args[i]);
MM_ERR(" sourcng  "<<nm<<"  ")
source(nm);
++i; // consume param and cmd
}
}


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
//typedef BranchesVar BraTy;
typedef fl_branches BraTy;
//typedef FickLancParamsVar FlpTy;
typedef fl_params FlpTy;

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

void source(const StrTy & fn)
{

CommandInterpretter li;
li.source(fn,true);
// so then this clverly connects to cin... doh 
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
if (cmd=="source") { if (li.cmd_ok(3))  li.source(li.word(1));  continue; } 
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
if (cmd=="initmesh") { InitMesh(); /*  m_cm.clear(); */  continue; } 
if (cmd=="dump-fixed") { dump_fixed(std::cout,local_label,3); /*  m_cm.clear(); */  continue; } 
if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="push-state") { InitPushOld();  continue; } 
if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="test") { test(li);  continue; } 
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd=="exit") { clean_up(); m_exit=true; return; } 
if (cmd=="interrupt") { debug_levels::interrupt();  continue; } 
if (cmd=="dump-lut") { m_idx.dump_lut(local_label);  continue; } 

MM_ERR(" unrecignized command "<<li.dump())

}




} //command_mode
// when exiting the interpretter
void clean_up()
{


} 
// tests, document what this crap should do lol

// the command object needs a shift function or something
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
#if 0
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
if (cmd=="mbflow") { test_mbflow(li); li.pop(); return ; } 
if (cmd=="mbintegral") { test_mbintegral(li); li.pop(); return ; } 
//void test_gr_integrator( CommandInterpretter & li )
#endif
MM_ERR(" unrecignized TEST command "<<li.dump())

li.pop();
} // test
//////////////////////////////////////////////////////////////




// real code approaching ,move the pixelt crap somewhere

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
#if 0 
gerrymander gm;
MM_ERR(" attempting to refine with thr="<<thr)
MM_MSG(" attempting to refine with thr="<<thr)
//(MyBlock & newpoints, const MyBlock & fixed, const MyBlock & sol, const IdxTy ix, const IdxTy iy, const D & thr)
MyBlock np(3); // this should also take an std::vector 
//gm.max_delta_refine(np,m_fixed,m_state,m_x,m_u,thr,~0); 
gm.max_delta_refine2(np,m_fixed,m_state,m_x,m_n,thr,~0); 
MM_ERR(" attempt to reinit with pointsize = "<<np.size()) 
MM_MSG(" attempt to reinit with pointsize = "<<np.size()) 
Init(0,np);
#endif
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
 	{	MM_STATUS(" "<<ival<<" "<<m_state(m_points-2+0,m_u)<<status_string()) }
//	{	MM_ERR(" iteration "<<ival<<" "<<m_state(m_points-2+0,m_u)<<status_string()) }
	if (m_tree.always_dump_state())  dump_state_etc();
 	else if (m_tree.dump_state_everyfew()&&everyfew)  dump_state_etc();

	if (m_tree.always_dump_to_file()) dump_to_default_file();
	
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
#if 1 
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
ss<<"                     "; // guard, although the status macro should erase the line doh
#endif
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

void dump_to_default_file()
{

	const StrTy fn=m_flp.dump_file();
	const StrTy label=m_flp.dump_label();
	const IdxTy flags=m_flp.dump_flags();
 	dump_to_file( fn ,  label, flags);
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
MM_ONCE(" qoi is fudded ",)
	if (column_names) { os<<"iters"<<""<< " "<< m_state.col_names()<<" "<<qoi.col_names()<<CRLF; }
	os<<m_state.to_string_cols_R(ival,12,12,qoi," ",false);   

}


template <class Os> void dump_fixed(Os & os, const StrTy & label="fick",const IdxTy flags=0)
{

	//const StrTy ival=(m_cm).format_entry("iterations");
	const StrTy lbl=label; // +ival;
	const bool column_names=(1==(flags&2));
	const bool human_readable=(0==(flags&1));
	if (column_names) { os<<"cols"<<lbl<< " "<< m_fixed.col_names()<<CRLF; }
	if (human_readable)  { os<<m_fixed.to_string_cols(lbl,2,8,MyBlock()," ");  }
	else { os<<m_fixed.to_string_cols(lbl,12,12,MyBlock()," ",false);  } 

}
template <class Os> void dump_state_etc(Os & os, const StrTy & label="fick",const IdxTy flags=0)
{
	//const StrTy lbl="fick"+ival;
//#if 0
MM_ONCE(" qoi is still fudded ",)
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
//#endif
}

// move this down into real code 
void compute_qoi( MyBlock & qoi)
{
#if 1 
	const IdxTy nqoi=3;
	qoi.resize(m_points,nqoi);
	qoi.col_name("current",0);
	qoi.col_name("x",1);
	qoi.col_name("y",2);

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
		const D phi_n=s(i,m_v_nx)*s(i,m_n) -m_Vt*m_mun*gradN;
		//const D phi_p=-m_mup*s(i,m_p)*(gradU)-m_Vt*m_mup*gradP;
		const D phi_p=s(i,m_v_px)*s(i,m_p)-m_Vt*m_mup*gradP;
		// this is charge current, the phis are particle current
		const D jc=phi_n-phi_p;
		qoi(i,0)=jc;
		qoi(i,1)=m_fixed(i,m_x); // no easy way to include location lol 
		qoi(i,2)=m_fixed(i,m_y); // no easy way to include location lol 
	} // i 
	// kluge
	qoi(0,0)=qoi(1,0);
	qoi(sz,0)=qoi(sz-1,0);
	qoi(0,1)=m_fixed(0,m_x); // no easy way to include location lol 
	qoi(sz,1)=m_fixed(sz,m_x); // no easy way to include location lol 

#endif

} // compute qoi 

void step(const IdxTy iter ) // const step_param & sp)
{
	// the ordering here is not physical but handy for testing... 
	const FlpTy & flp= m_flp;
	const bool enable_trapping=m_tree.enable_trapping();
	const bool enable_gr=m_tree.enable_gr();
	const bool enable_flmb=m_tree.enable_flmb();
	const bool sanity_check=!true;
	const D dt=flp.time_step(); 
	m_cm.mark("step");
	m_cm.mark("find_pot");
	// this can be FD or FEM but needs diri BC. 
	solve_for_potential(m_state);
	update_drift(m_state); 
	m_cm.cumdel("pottotal","find_pot");
	if (enable_flmb)
	{
		m_cm.mark("mb");
		flmb_flow(dt,m_mun,m_n,m_v_nx,m_v_ny,flp);
		flmb_flow(dt,m_mup,m_p,m_v_px,m_v_py,flp);

		m_cm.cumdel("mbtotal","mb");
		if (sanity_check) check_carriers("mbfl")	;
	}
	impose_bc();
	if (enable_gr)
	{
		m_cm.mark("find_gr");
//		gr_attributing_sf(m_state, dt, m_n,m_p, m_points);
		gr_point_by_point(m_state, dt, m_n,m_p, m_points);
		m_cm.cumdel("grtotal","find_gr");
	}
	if (sanity_check) check_carriers("gr")	;
//	enforce_bc(m_state); // ????? 
	if (enable_trapping){
	m_cm.mark("find_traps");
	const bool do_sf_trap_step= !true; 
	if (do_sf_trap_step)
	{
// this itor gone 
	//	Di din(m_points,m_x,m_n,m_state,m_fixed);
//		Di dip(m_points,m_x,m_p,m_state,m_fixed);
//	 	m_defects.step( din, dip,m_fixed, m_state ,dt,m_u, m_x);
	} else { m_defects.step(m_fixed,m_state,dt,m_points,m_n,m_p,m_u,m_x); }
	m_defects.net_charge(m_state,m_trap,m_points,1);
	m_cm.cumdel("traptotal","find_traps");
} // enable_trapping

	impose_bc();
	m_cm.cumdel("steptotal","step");


} // step 
void oldstep(const IdxTy iter ) // const step_param & sp)
{
const FlpTy & flp= m_flp;
const D t=flp.time_step(); // fl_t(); //1e-11; // sp.dtau;
// this should not be used now 
//const D dxmax=flp.fl_dmax(); // m_h*35; // sp.dxmax;
//const IdxTy pmax=m_points-0;
//for (IdxTy point=0; point<pmax; ++point)
{
	m_cm.mark("step");
	// examine field to come up with velocity distributions
	// as drift shifts this at each point
	m_cm.mark("find_pot");
	// if the shotting fixed the initial e field it would be ok 
} // point


} // step

// this is probably not really needed or makes confusing results 
void enforce_bc(MyBlock & sol)
{
#if 0

#endif
}

// needs nancheck
IdxTy count_negs(MyBlock & state, const IdxTy col)
{
	//typedef mjm_numb Ck;
	IdxTy k=0; 
#if 1 

	const IdxTy sz=m_points;
	for (IdxTy i=0; i<sz; ++i) {
		 if ( state(i,col)<0) {++k;
		MM_MSG(" negative carrier "<<i<<" col "<<col<<" value "<<state(i,col))

} } // for

#endif
return k;

}

// make this a member of MyBlock and do something with FILE and LINE lol 
IdxTy count_nans(MyBlock & state, const IdxTy col)
{
	typedef mjm_numb Ck;
	IdxTy k=0; 
#if 1 
	const IdxTy sz=m_points;
	for (IdxTy i=0; i<sz; ++i) 
	{ 
		if (Ck::denan(state(i,col) ,__FILE__,__LINE__," fudd "))
		{ ++k; 
			MM_MSG(" nan found  for carrier  "<<MMPR3(i,col,state(i,col)))
		}
	}
#endif
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

D rho_fixed(const IdxTy node)
{
return m_fixed(node,m_nd)-m_fixed(node, m_na)+m_state(node,m_trap);
}
D total_q(const IdxTy node)
{
//MM_MSG(MMPR3(m_state(node,m_p),m_state(node,m_n),rho_fixed(node)))
return rho_fixed(node)+m_state(node,m_p)-m_state(node,m_n);
}


// passing state ignores the rhs and boundaries. 
void solve_for_potential(MyBlock & state)
{
	mjm_special_qrule mats;
	MyBlock b(m_nodes);
//	b.resize(m_nodes);
	MySparse A;
	MyBlock mat; // FE coefficients for laplacian 
	mat.resize(4,4); // doh, surprised it ran at all wtf 
	// iterate over the node, needing nearest neighbors to define del^2
	// need fixed charge and Dirichlet BC's.
 	// or just do this like real FEM and integrate over elements... 
	DofVec dof_indices;
//    (*m_es).print_info();
	const IdxTy nvar=1;
	const IdxTy nnode=4;
	const IdxTy nv=nnode*nvar;
 	dof_info di =dof_info((*m_jr).get_dof_map(),nvar); 	
IdxTy serial=0;
//    dof_map.dof_indices (elem, dof_indices);
//    for (IdxTy var=0; var<nvar; var++)
	//dof_map.dof_indices (elem, dof_indices_var[var], var);
	// mesh element itor, get the laplacian code and JxW etc. 
	MeshCItor      el     =(*m_mesh).active_local_elements_begin();
    const MeshCItor end_el = (*m_mesh).active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
		const Elem * elem=*el;
		di.re_init(elem); // FUDD 
		const Node & ll=*((*elem).node_ptr(0));
		const Node & ur=*((*elem).node_ptr(2));
//	MM_MSG(" shot "<<shot.n_dofs(0,0)<<" "<<shot.id())
		//const D dx= (*(*elem).node_ptr(0))(0)-(*(*elem).node_ptr(2))(0);
		const D dx= ll(0)-ur(0);
		const D dy= ll(1)-ur(1);
		//const D dy= (*(*elem).node_ptr(0))(1)-(*(*elem).node_ptr(2))(1);
//MM_MSG(" fudd "<<MMPR2(dx,dy))
		//const D b= (*elem).node_ptr(0)->(1)-(*elem).node_ptr(2)->(1);
//		di.dof_map.dof_indices(elem, dof_indices);
		di.dof_map.dof_indices(elem, dof_indices, 0);
//MM_MSG(" fudd ")
		// this is incredibly dumb, most of a and b are the same, just do this once 
		mats.lap_mat(mat,dx,dy);
		di.sum_dof_mats(A, mat, dof_indices, nv);
		mats.lap_mat_y(mat,dx,dy);
		di.sum_dof_mats(A, mat, dof_indices, nv);
//MM_MSG(" fudd ")
		// needs to be integrated, scaled, and mult b q/eps
		for (IdxTy i=0; i<nnode; ++i)
		{
//			const Node & no=*((*elem).node_ptr(i));
			const IdxTy dofudd=dof_indices[i];	
			const D rhoqe=m_qe*total_q(dofudd)*.25*dx*dy;
	//		MM_MSG(MMPR4(serial,i,dofudd,rhoqe) )
			b(dofudd)+=rhoqe;
		}	
	++serial;
	}

// rhs is the integrated net charges

// get dirichlets 
	impose_diri( &A,&b,1.0);
	impose_bc(); // should be elsewhere. 
	(*m_jr).set_system(A,b);
	(*m_jr).mjm_linear_solve();
	(*m_jr).copy_solution(state,m_u);

}


void impose_diri( MySparse * A, MyBlock * rhs , const D & scale)
{
// moved to bc code 
//const D vd=m_flp.vd();
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(i,true)]=0;
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;

auto di=m_diri.begin();
auto de=m_diri.end();
while (di!=de)
{
const IdxTy dof=(*di).first;
const D val=(*di).second;
if (A) { A->zero_row(dof); (*A)(dof,dof)=scale; }
if (rhs) { (*rhs)(dof)=val*scale; } 
++di;
}


}

void  velocity( D & vx, D & vy, const D & mu
	, const D & vsat, const D & ex, const D & ey, const D & e) const
{
//const D e=sqrt(ex*ex+ey*ey);
const D mueff=mu/(1.0+mu*e/vsat);
vx=mueff*ex;
vy=mueff*ey;
} 


void update_drift(MyBlock & state)
{
    // this is bad with zero points, fwiw
//    state(0,m_u)=m_v0; state(sz,m_u)=m_vd; state(1,m_u)=m_v0;
    const D vsat=m_flp.vsat(); // 1e5;
    const IdxTy iu=m_u;
    const IdxTy sz=m_points;
//    const bool do_2_vel=true;
//    if (do_2_vel) { MM_ONCE(" doing split velocity ",) }
    IdxTy iinit=0;
	const IdxTy bad=m_idx.bad();
	D vx=0, vy=0; 
//    if (do_2_vel ) iinit=0;
	// these are node id's and need to first become location x+nxmesh*y
	// translated and converted back to node id. 
    for (IdxTy i=iinit; i<sz;  ++i)
    { // this really does need to stay on grid for now 
		// and use the lut 
		const IdxTy north=m_idx.north(i);
		const IdxTy south=m_idx.south(i);
		const IdxTy east=m_idx.east(i);
		const IdxTy west=m_idx.west(i);

		const IdxTy top=(north==bad)?i:north;
		const IdxTy bot=(south==bad)?i:south;
		const IdxTy left=(west==bad)?i:west;
		const IdxTy right=(east==bad)?i:east;
		// this probably fails because the first 2 rows are in wrong arrangement,
		// check and skip 
		const D denx= (m_fixed(right,m_x)-m_fixed(left,m_x));
		const D ex=(denx==0)?0: ((m_state(right,m_u)-m_state(left,m_u))/(denx));
		const D deny=(m_fixed(top,m_y)-m_fixed(bot,m_y));
		const D ey= (deny==0)?0:((m_state(top,m_u)-m_state(bot,m_u))/(deny));

		const D e=sqrt(ex*ex+ey*ey);
		velocity( vx, vy, m_mup, vsat, ex, ey,e);
		m_state(i,m_v_px)=-vx;
		m_state(i,m_v_py)=-vy;
		velocity( vx, vy, m_mun, vsat, ex, ey,e);
		m_state(i,m_v_nx)=vx;
		m_state(i,m_v_ny)=vy;

        // for the heck of it track rho too
		// need to get traps etc 
        state(i,m_q)=total_q(i); // rho_fixed(i)+state(i,m_p)-state(i,m_n);
    } // i 
}

// redundant ( ie dx, dy, r2 ) param sets are to save computations 
void kernel_gauss(D * gni, D * gnorm, const D & beta, const D & r2, const D & v_x_t,
const D & dx, const D & dy, const D & vx, const D & vy, const D & t,
const D & n0, const D & ni )
{// this means the first part mus be calculated twice, fudd 
	const D q2=(-beta*(r2+v_x_t*v_x_t));
	const D vexp=exp(q2-2.0*beta*(vx*t*dx+vy*t*dy));
if (gnorm!=0)
{
	const D  gauss=vexp; // exp(q2);
//	if (gnorm!=0) 
	*gnorm+=gauss;
}
	if (gni==0) return; 
	// wt ASSFUDD is wrong with this fudd er 
//	const D vexp=exp(q2+2.0*beta*(vx*t*dx+vy*t*dy));
	// well duh, need to have outbound since that is the relevant v 
	*gni+=(n0*vexp); // -ni/vexp);

}

void flmb_flow(const D & dt,const D & mu, const IdxTy in, const IdxTy ivx, const IdxTy ivy,const FlpTy & flp)
{
 //   { MM_ONCE(" doing flmb ",) }
 //   { MM_ERR(" doing flmb ") }
	const bool dump1=false;
	const bool dump2=false; 
	const bool dump3=false; 
	const bool dump4=m_tree.enable_dump4(); 
	MyBlock next(m_nodes);
	typedef MeIdx::hood_iterator Hi; 
 static const D bereal_cgs=1e-4*9.1e-31/4.11e-21;
    const D vmax=m_flp.vmax();
    const D tau=dt; // this is just the time of flow. ilonger more diffuse
    // the integrand is exp(-beta*(x-x'-v*t)^2)
    const D beta= bereal_cgs*mu/2000.0/vmax/vmax/tau/tau; // 1e-4;   // m/(2kT) lol 

	const D eta=5e-2;
	const D radius=::sqrt(-log(eta)/beta);
	MM_ONCE(MMPR4(radius,beta,bereal_cgs,mu)<<MMPR2(vmax,tau),)
	IdxTy serial=0;
	const NodeItor ie=(*m_mesh).nodes_end();
	for (NodeItor ii=(*m_mesh).nodes_begin(); ii!=ie; ++ii)
	{
		const Node & np=*(*ii);
//	 the lut takes the geo location and returns this id, 
// need to have opposite too... 
		const IdxTy id=np.id();
		const int idgeo=m_idx.lut(id,false);
		const D n0=m_state(id,in);
		const D vx=m_state(id,ivx); // this node it the src node
		const D vy=m_state(id,ivy); 
		const D tvx=vx*dt;
		const D tvy=vy*dt;
		const D vxt=::sqrt(tvx*tvx+tvy*tvy);
		// the iterator is based on geo coords
		Hi hi=m_idx.hood_itor(radius,idgeo,tvx,tvy);
		// now loop over all nodes that can communicate and make the changes
		// need to keep track of stuff for normalization 		
		D gni=0;
		D gnorm=0;
		hi.start();
		while (hi.ok())
		{	// the itor is over geographical nodes 
			const int meshx=hi.meshx();
			const int meshy=hi.meshy();
			const IdxTy dnode=m_idx.index(meshx,meshy); // hi.node();
			// off-mesh nodes are legitimate at the +/- y ends right now. 
			// the walk off increases norm and hence amount transfered to self later. 
			int dnodeloc=0;
			if (dnode<0) 
			{ 
				const IdxTy oorcode=m_idx.oor(meshx,meshy);
				// the oor layers at the y ends are the bulk media
				if (( oorcode|12)==0){ hi.inc(); continue; }
				// this has to find the right one to get the doping level although
				// really now it is just location as it is only outflow 
				dnodeloc=dnode; // m_idx.index(meshx,meshy);
			} 
			// geo coords are turned into mesh id 
			else  dnodeloc=m_idx.lut(dnode,true);
//			if (dnodeloc<0) { hi.inc(); continue; } 
			const D ni=(dnodeloc<0)?
				m_idx.m_peri_state(-dnodeloc-1,in):m_state(dnodeloc,in);	
			//kernel_gauss(gni, gnorm, beta, r2, v_x_t, dx, dy, vx, vy, t, n0, ni);
			if (dump2)	MM_MSG(hi.to_string())
			if (dump2)		
				MM_MSG(MMPR4(gnorm,beta,hi.r2(),vxt)<<MMPR4(hi.delx(),hi.dely(),vx,vy)<<MMPR3(dt,n0,ni))
			// the delx() are DEST-SRC opposite of this convention 
			kernel_gauss(0, &gnorm, beta, hi.r2(), vxt, -hi.delx(), -hi.dely(), vx, vy, dt, n0, ni);
			hi.inc();
		}
		hi.start();
		while (hi.ok())
		{
			//const IdxTy dnode=hi.node();
			const IdxTy dnode=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
			if (dnode<0) { hi.inc(); continue; } 
			const int dnodeloc=m_idx.lut(dnode,true);
/*
		const IdxTy id=dnodeloc;
		const D vx=m_state(id,ivx); // this node it the src node
		const D vy=m_state(id,ivy); 
		const D tvx=vx*dt;
		const D tvy=vy*dt;
		const D vxt=::sqrt(tvx*tvx+tvy*tvy);
*/


			if (dump1) MM_MSG(hi.to_string())
			if (dump1) MM_MSG(MMPR4(hi.meshx(),hi.meshy(),dnode,dnodeloc))
			if (dnodeloc<0) { hi.inc(); continue; } 
//			const D ni=m_state(dnodeloc,in);	
			const D ni=(dnodeloc<0)?
				m_idx.m_peri_state(-dnodeloc-1,in):m_state(dnodeloc,in);	
			gni=0;
			//kernel_gauss(gni, gnorm, beta, r2, v_x_t, dx, dy, vx, vy, t, n0, ni);
			// del are dest - src 
			kernel_gauss(&gni, 0, beta, hi.r2(), vxt, -hi.delx(), -hi.dely(), vx, vy, dt, n0, ni);
			// these are OUTBOUNMD doh 
			//if (gnorm!=0) next(id)=gni/gnorm;
			// this is OUTFLOW from the current node to neighboring nodes.a
			// flow is ADDED to destinations including self. Kluged norm conserves lo Kluged norm conserves loll 
			if (dump1) MM_MSG(MMPR4(gni,gnorm,hi.meshx(),hi.meshy())<<MMPR2(dnodeloc,gni/gnorm))

			if (gnorm!=0) next(dnodeloc)+=gni/gnorm;
	//		MM_MSG(hi.to_string()<<MMPR3(vxt,vx,vy)<<MMPR4(gni,gnorm,ni,n0))
			hi.inc();
		}

		++serial; 	
	}
	// this also needs to loop over the peripheral nodes for outflow. 
	// for now just inflow at y limits to avoid issues at the junction.
	typedef mjm_rectangle_itor ReI;
	const D bcf=1;
	ReI rei(-m_layers,0,m_n_x,m_n_x+m_layers, -m_layers,0,m_n_y,m_n_y+m_layers);
	const int base= -m_layers-m_layers*(m_n_x+2*m_layers); 
	while (rei.ok())
	{
		// this should be negtive
		int node=m_idx.index(rei.x(),rei.y());
		if (dump4) MM_MSG(MMPR3(rei.x(),rei.y(),node))
		//MM_MSG(MMPR2(node,rei.to_string()))
		if (node>=0)
		{
			MM_ERR(" node is i nmesh "<<MMPR3(node,rei.x(),rei.y()))
		}
		const D n0=(node<0)?
		m_idx.m_peri_state(-node-1,in):m_state(node,in);	
		D gni=0;
		D gnorm=0;
		// geo_index is zero based including excludinged area 
		// so the (0,0) hood_itor uses needs to be tanslated 
		if (dump4) MM_MSG(MMPR3(rei.geo_idx(),base,radius))
		// -base makes the input wrong- geo_idx returns mesh coordinates, which is right  
		Hi hi=m_idx.hood_itor(radius,rei.geo_idx()-base,0,0,true);
		if (dump4) MM_MSG(MMPR(hi.to_string()))
		hi.start();
		while (hi.ok())
		{	// the itor is over geographical nodes 
			const int dnode=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
if (dnode>0)
{
		const IdxTy id=m_idx.lut(dnode,true); // dnodeloc;
		const D vx=m_state(id,ivx); // this node it the src node
		const D vy=m_state(id,ivy); 
		const D tvx=vx*dt;
		const D tvy=vy*dt;
		const D vxt=::sqrt(tvx*tvx+tvy*tvy);



			const D ni=(dnode<0)?
				m_idx.m_peri_state(-dnode-1,in):m_state(id,in);	
			//kernel_gauss(gni, gnorm, beta, r2, v_x_t, dx, dy, vx, vy, t, n0, ni);
			kernel_gauss(0, &gnorm, beta, hi.r2(), vxt*bcf, -hi.delx(), -hi.dely(), vx*bcf, vy*bcf, dt, n0, ni);
} else
			kernel_gauss(0, &gnorm, beta, hi.r2(), 0, -hi.delx(), -hi.dely(), 0, 0, dt, n0, 0);



			hi.inc();
		}

		hi.start();
		while (hi.ok())
		{	// the itor is over geographical nodes 
		// this can diffuse into multiple locations, need to get norm right. 	
			const int  dnode=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
			// the normalization is supposed to save this but ill not work here 
			if (dnode<0) { hi.inc(); continue; } 

		const IdxTy id=m_idx.lut(dnode,true); // dnodeloc;
		const D vx=m_state(id,ivx); // this node it the src node
		const D vy=m_state(id,ivy); 
		const D tvx=vx*dt;
		const D tvy=vy*dt;
		const D vxt=::sqrt(tvx*tvx+tvy*tvy);



			const D ni=(dnode<0)?
				m_idx.m_peri_state(-dnode-1,in):m_state(m_idx.lut(dnode,true),in);	

			gni=0;
		// skip while thing if dnone<0
			kernel_gauss(&gni, 0, beta, hi.r2(), vxt*bcf, -hi.delx(), -hi.dely(), vx*bcf, vy*bcf, dt, n0, ni);
			if (dnode>=0) if (gnorm!=0) next(id)+=gni/gnorm;
		if (dump4) 	MM_MSG(" bc "<<MMPR4(hi.meshx(), hi.meshy(),gni/gnorm,hi.r2())<<MMPR4(hi.delx(),hi.dely(), n0,ni)<<MMPR4(dnode,dt,gni,gnorm)) 
			hi.inc();
		}


		rei.inc();
	}

	for (IdxTy i=0; i<m_nodes; ++i)  m_state(i,in)=next(i);

//    { MM_ERR(" done flmb ") }
} // flmb_flow


// probably just integrate this, it should not cross equ or collect 200 dollars lol 
// can be done as a function of position however physics becomes
// a question- maybe it smears out anyway and nothing is added 

/////////////////////////////////////////////////////////
// calculate interpolated gr and attribute to each notde
// this needs to be log n but not now... 
// do ELEMENTS, not NDES
// log interpolation helped nothing. 
// this probably does not do geometry right FIXME
// for quantiatitve agreement
void gr_point_by_point(MyBlock & sol, const D & dt, const IdxTy in, const IdxTy ip, const IdxTy points)
{
	const FlpTy & flp= m_flp;
	const bool debug_gr=true;
	const D fudge= flp.gr_xsection(); // 1e-3;
	const D tau= flp.gr_tau();
	const D phitau=fudge/tau;
	MyBlock diffn(points),diffp(points);
	for (IdxTy i=0; i<points; ++i)
	{
		const D n=m_state(i,in);
		const D p=m_state(i,ip);

		D dn=gr_b2b_point((n),(p),m_ni2,phitau,dt);
		diffn(i)=dn;
		diffp(i)=dn;
	}
	for (IdxTy i=0; i<points; ++i) 
	{ 
		m_state(i,in)+=diffn(i); m_state(i,ip)+=diffp(i); 
	}

}

#if 0
void gr_attributing_sf(MyBlock & sol, const D & dt, const IdxTy in, const IdxTy ip, const IdxTy points)

#endif


////////////////////////////////////////////////////////////////////////

D gr_b2b_point(const D & n , const D & p, const D & ni2
	, const D & phitau, const D & t)
{
D dn=0;
const D expmax=300;
const D expmin=-300;
IdxTy oflow=0;
IdxTy uflow=0;
const D c=-phitau; 
const bool cond=(p<n); // true; // (p<n);
const D nmin=(cond)?p:n;
const D nmax=(cond)?n:p;

const D d0=nmax-nmin; // p-n; 
// FIXME d==0 still allows partial fractions
// but the expressions are diferent 
const D Di= ::sqrt(d0*d0+4.0*ni2);
const D disc=.5*Di;
const D nplus=-.5*d0+.5*Di;
const D nminus=-.5*d0-.5*Di;
const D avg= .5*(nmin+nmax);
//const D ctime=c/d0;
const D tau=Di*c*t;
if (false) MM_MSG(MMPR3(tau,c,t)<<MMPR4(nplus,nminus,d0,disc))
//if (tau>expmax) { return xminus-nmin; } 
if (tau>expmax) { return nminus; } 
//if (tau<expmin) { return xplus-nmin; } 
if (tau<expmin) { return nplus; } 
const D e=exp(tau);
//const D den=(nmin-xplus)*e-(nmin-xminus);
const D den=(avg-disc)*e-(avg+disc);
//const D num=xminus*(nmin-xplus)*e-xplus*(nmin-xminus);
const D num=(nminus-nmin)*(avg-disc)*e-(nplus-nmin)*(avg+disc);
if (false) MM_MSG(MMPR4(den,num,e,nplus)<<MMPR(nminus))
dn=num/den; // -nmin;
if ( true ) return dn;
}
////////////////////////////////////////////////////////////////////

// this is a copy of the code in the point wise gr method
// which can be removed once this seems to work and
// oor handling is fixed 
#if 0
D gr_b2b_point_old(const D & n , const D & p, const D & ni2
	, const D & phitau, const D & t)
#endif

#if 0
#endif
void impose_bc()
{
impose(m_n,m_bc_n);
impose(m_p,m_bc_p);


}
template <class Ty> void impose( const IdxTy col, Ty & m)
{
for(auto i=m.begin(); i!=m.end(); ++i)
 m_state((*i).first,col)=(*i).second;
}


void InitBC()
{
	m_bc_n.clear();
	m_bc_p.clear();

for ( IdxTy i=0; i<m_n_x; ++i)
{
const int idx=m_idx.lut(i,true);
const int idxe=m_idx.lut(m_nodes-i-1,true);
 m_bc_n[idx]=m_n_p;
 m_bc_p[idx]=m_p_p;
 m_bc_n[idxe]=m_n_n;
 m_bc_p[idxe]=m_p_n;


}


	m_diri.clear(); // should be in init	
const D vd=m_flp.vd();
for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(i,true)]=0;
for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;

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
//	const D fac=1.0*point/(pmax-1);
		const bool lr= (m_fixed(point,m_y)>m_y_junction) ;
//{  m_fixed(i,m_nd)=Rev(Nd); m_fixed(i,m_na)=0;	} 
	if (point<m_mid) { sol(point,iu)= m_v0; } else { sol(point,iu)= m_vd; }
	if (!lr) { sol(point,in)= m_n_p; sol(point,ip)=m_p_p; } 
	else { sol(point,in)= m_n_n; sol(point,ip)=m_p_n; }
	sol(point,iz)=1;
}
	if (m_tree.enable_trapping())
	m_defects.net_charge(m_state,m_trap,m_points);

for (IdxTy point=0; point<(-m_idx.m_index_base); ++point)
{
if (m_idx.m_peri_fixed(point,m_y)<=m_y_junction)
{
m_idx.m_peri_state(point,in)=m_n_p; 
m_idx.m_peri_state(point,ip)=m_p_p; 
}
else
{
m_idx.m_peri_state(point,in)=m_n_n; 
m_idx.m_peri_state(point,ip)=m_p_n; 
}

}
	//m_idx.make_peripheral_nodes( nx, ny, m_layers, m_n_state, m_n_fixed);


}

void SetCols()
{
m_fixed.col_name("x",m_x);
m_fixed.col_name("y",m_y);
m_fixed.col_name("Na",m_na);
m_fixed.col_name("Nd",m_nd);


m_state.col_name("voltage",m_u);
m_state.col_name("z",m_z);
m_state.col_name("h",m_p);
m_state.col_name("n",m_n);
m_state.col_name("v_hx",m_v_px);
m_state.col_name("v_hy",m_v_py);
m_state.col_name("v_ex",m_v_nx);
m_state.col_name("v_ey",m_v_ny);
m_state.col_name("rho",m_q);

}
#if 0
#endif
void InitBase()
{
// TODO merge fixed and state and include this crap in a new obkect
IdxTy pos=0; // enum? lol 
 	m_u=pos; ++pos; // 1
	m_z=pos; ++pos; //2
	m_p=pos;  ++pos; // 3
	m_n=pos; ++pos; // 4
	m_v_px=pos; ++pos; // 5
	m_v_py=pos; ++pos; // 6
	m_v_nx=pos; ++pos; // 6
	m_v_ny=pos; ++pos; // 6
	m_q=pos; ++pos; //9 
	m_trap=pos; ++pos; //10 
	m_x=0; m_y=1;m_na=2; m_nd=3;
    //m_best_l2=-1;
    //m_last_l2=-1;
    m_vars=10; // pos??? 
    //m_points=m_size/m_vars;
    m_size=m_points*m_vars;
    m_iter=0;
	m_exit=false;
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

	//m_fixed.resize(m_points,3);
	//m_h=m_flp.total_length()/D(m_points-1); // 1e-6; // good value
	// this should be in the config file 
	m_y_junction=.5*(m_flp.mesh_ymin()+m_flp.mesh_ymax()); // 1e-6; // good value
	m_points=m_flp.mesh_nx()*m_flp.mesh_ny();
    m_size=m_points*m_vars;
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
		if (m_fixed(i,m_y)>m_y_junction) {  m_fixed(i,m_nd)=Rev(Nd); m_fixed(i,m_na)=0;	} 
		else {ilast =i; m_mid=i+1;  m_fixed(i,m_na)=Rev(Na); m_fixed(i,m_nd)=0;	} 
	}
	if (m_tree.enable_trapping())
	{
		m_defects.set(1,m_points);

		typedef Defects::defect_level Level;
		Level level(m_flp.trap_speed());
		//m_defects.add(1e13);
		//m_defects.add(1e13,level);
		m_defects.add(m_flp.trap_N(),level);
m_defects.set_to_equ( m_n_n, m_p_n, m_ni2, 0,  Nd,m_points-1);
m_defects.set_to_equ( m_n_p, m_p_p, m_ni2, Na,  0,0);
//MM_MSG("traps bc "<< MMPR(m_n_n)<<MMPR( m_p_n)<<MMPR(m_p_p)<<MMPR(m_n_p)<<MMPR(m_x_junction))
//void set_to_equ(const IdxTy first, const IdxTy end, const D & ni2, const D & Na, const D & Nd)
m_defects.set_to_equ(0, ilast+1, m_ni2, Na, 0);
m_defects.set_to_equ(ilast+1, m_points, m_ni2, 0,  Nd);
//MM_MSG(MMPR(ilast));
// does nothing as it is reset when initial_guess is done 
// well it does fudd everything up since m_state size is not right yet 
//	m_defects.net_charge(m_state,m_trap,m_points);

	}
}
#if 0
void ProjectOldCrap( const MyBlock & osol, const MyBlock & ofixed, const Defects & odefects, const IdxTy opoints)
{
	const IdxTy sz=m_points; 
	const IdxTy olim=opoints-1;
	MM_MSG(" Project "<<MMPR(sz)<<MMPR(olim))
	m_state.resize(m_points,m_vars);
	if (opoints<2) return ; 
	const bool do_traps=m_tree.enable_trapping();
	IdxTy i1=0;
	IdxTy i2=1;
	for (IdxTy i=0; i<sz; ++i)
	{
		const D x=m_fixed(i,m_x);
		D x2=ofixed(i2,m_x);

		while (i2<olim) { if (x2>=x) break ; ++i2; x2=ofixed(i2,m_x); }
		i1=i2-1;
		D x1=ofixed(i1,m_x);
		D f2=( x-x1)/(x2-x1);
		D f1=1.0-f2;
		// really this needs to use all the old m_x for robustness 
		m_state(i,m_u)= f1*osol(i1,m_u)+f2*osol(i2,m_u);
		m_state(i,m_n)= f1*osol(i1,m_n)+f2*osol(i2,m_n);
		m_state(i,m_p)= f1*osol(i1,m_p)+f2*osol(i2,m_p);
	//	MM_MSG(" Project i="<<i<<MMPR3(m_state(i,m_u),m_state(i,m_n),m_state(i,m_p))<<MMPR4(x1,x2,x1,f2))
		m_state(i,m_z)=1;
		// in reality the other Init may have provided more accurate
		// "N" values as resolution improves. 
		if (do_traps)
		{
//void project_from( const Myt & old, const IdxTy point,const IdxTy  i1, const IdxTy i2, const D & f1, const D & f2)
			m_defects.project_from( odefects, i, i1,  i2, f1, f2,false);

		}
	}
		m_state(0,m_u)= osol(0,m_u);
		m_state(0,m_n)= osol(0,m_n);
		m_state(0,m_p)= osol(0,m_p);

		m_state(sz-1,m_u)= osol(olim,m_u);
		m_state(sz-1,m_n)= osol(olim,m_n);
		m_state(sz-1,m_p)= osol(olim,m_p);

		if (do_traps) { 
			m_defects.project_from( odefects, 0, 0,  0, 1, 0,false);
			m_defects.project_from( odefects, sz-1, olim,  olim, 0, 1,false);

//			MM_ERR(" traps not projected ") 
			MM_ONCE(" traps being projected danger will robinson ",) 
//			m_defects.project_from(odefects,
		} 
//	m_defects.net_charge(m_state,m_trap,m_points);
}

#endif
void InitSol()
{
	if (!m_save_sol) initial_guess(m_state);
	else
	{
	if (m_old==0) { 	MM_ERR(" m_old needs to be non null" ) } 
		Myt & old = *m_old;
		//ProjectOldCrap(old.m_state,old.m_fixed, old.m_defects, old.m_points);
		if (m_tree.enable_trapping())
			m_defects.net_charge(m_state,m_trap,m_points);
	//		if (!m_state.form_is(m_points,m_vars))
//			{
//				MM_ERR(" resizing m_state from"<<m_state.form().to_string())	
//				m_state.resize(m_points,m_vars);
				// there needs to be a projection operation som where. 
				// as well as saved location info 
//			}

	}
} // InitSol


void InitVerify()
{
// this needs to be part of banner and have more for dev crap 
MM_MSG( MMPR(m_n_n)<<MMPR( m_p_n)<<MMPR(m_p_p)<<MMPR(m_n_p)<<MMPR(m_y_junction))
// points monotonic, no nans or negatives etc. 
// perhaps even hmax/hmin diagnose and critique

}
//#endif
void InitEs()
{
//(*m_mesh).prepare_for_use(true);
delete m_es;
//delete m_jr;
const StrTy name="poission";
m_es= new EsTy(*m_mesh);
//m_jr = new  NonLinTy(*m_es, name,1);
 (*m_es).add_system<NonLinTy> (name);
  NonLinTy & system=  (*m_es).get_system<NonLinTy>(name);
m_jr=&system; 
(*m_jr).project_solution_on_reinit()=!true;
  //system_type & system=  * equation_systems.add_system<system_type>(name);
  system.add_variable("u", libMesh::FIRST);
MM_MSG(" variable added ")
//(*m_jr).init();
// assemble linear sol by hand 
// system.attach_assemble_function (assemble_bc);
//    system.attach_init_function (init_cd);
//    system.project_solution_on_reinit()=true;

        (*m_es).init();
//        (*m_es).print_info();

m_jr->get_dof_map().reinit(*m_mesh);
(*m_mesh).prepare_for_use(true);
//(*m_mesh).print_info();
(*m_es).reinit();
}


void InitMesh()
{
delete m_mesh;
  m_mesh= new MyMesh(m_init->comm());

const FlpTy& flp= m_flp; // right now this is being SET in Init not being used to set m_h
const IdxTy mesh_nx=flp.mesh_nx(); // 
const IdxTy mesh_ny=flp.mesh_ny(); // 
const D mesh_xmin=flp.mesh_xmin(); //
const D mesh_ymin=flp.mesh_ymin(); // 
const D mesh_xmax=flp.mesh_xmax(); // 
const D mesh_ymax=flp.mesh_ymax(); // 
MM_MSG(" making mesh "<<MMPR4(mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax)<<MMPR2(mesh_nx,mesh_ny))
  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
// these counts are ELEMENTS NOT NODES... 
  libMesh::MeshTools::Generation::build_square (*m_mesh,
                   //                    mesh_nx, mesh_ny,
                                       mesh_nx-1, mesh_ny-1,
                                       mesh_xmin, mesh_xmax,
                                       mesh_ymin, mesh_ymax,
                                       //QUAD4);
                                       ElemTy);
	m_nodes=m_mesh->n_nodes();
	m_points=m_nodes;
	m_n_x=mesh_nx;
	m_n_y=mesh_ny;
	m_dx=(mesh_xmax-mesh_xmin)/(m_n_x-1);
	m_dy=(mesh_ymax-mesh_ymin)/(m_n_y-1);

  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  m_mesh->print_info();

	//m_diri.clear(); // should be in init	
}
void InitFromMesh()
{
// later????
(*m_mesh).prepare_for_use(true);

 	m_n_state=10;
	m_n_fixed=4;
	m_layers=4;
// these should be one thing 
	m_fixed.resize(m_nodes,m_n_fixed);
	m_state.resize(m_nodes,m_n_state);
	const NodeItor ne=m_mesh->nodes_end();
	for (NodeItor ni=m_mesh->nodes_begin(); ni!=ne; ++ni)
	{
		const Node & n=*(*ni);
		const IdxTy node_n=n.id();
		if (node_n>=m_nodes) { MM_MSG(" bad node "<<MMPR2(node_n,m_nodes))} 
		m_fixed(node_n,m_x)=n(0);
		m_fixed(node_n,m_y)=n(1);



	} // ni

const D mesh_xmin=m_flp.mesh_xmin(); //
const D mesh_ymin=m_flp.mesh_ymin(); // 
//const D mesh_xmax=flp.mesh_xmax(); // 
//const D mesh_ymax=flp.mesh_ymax(); // 

	m_idx=MeIdx(*m_mesh,m_flp,m_tree);
	m_idx.load_regulars();
	m_idx.make_peripheral_nodes( m_n_x, m_n_y, m_layers, m_n_state, m_n_fixed);
	// now need to locate and popylate the fakes layers 
	m_idx.set_positions(m_x,m_y,mesh_xmin,mesh_ymin,m_dx,m_dy);

}

void Init(const IdxTy points=0, const MyBlock & x= MyBlock(1) ) {
	InitBase();
	InitConsts();
	InitMesh();
	InitFromMesh();
	InitFixeds();
	InitBC();
	InitEs();
	InitSol();
	SetCols(); // guess resizes trashing cols possible. 
	InitVerify();
#if 0
	if (m_save_sol)
	{
	// user must do this 
	//	InitPushOld();
	}

#endif

}
#if 0
// note that this does not always get the old config
void InitPushOld()
{
// note that this does NOT get all the old config at
// the user has already reset that 
const bool stacking=false;
if (stacking)
{
m_old = new Myt();
*m_old=*this;
// m_old->m_old is still valid and the stack is a linked list
// useful for retarded time later too. Deleting one of them
// should delete all the ancestors as m_old is deleted in dtor 
//m_old->m_old=0; // lol 
return; 
}

delete m_old;
m_old = new Myt();
*m_old=*this;
// the old m_old was already deleted unless we really are stacking
m_old->m_old=0; // lol 

}

#endif

D Fwd(const D & x) const { return (x); } 
D Rev(const D & x) const { return (x); } 

Myt * m_old;
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
D m_Vt,m_qe;
//D m_ni2,m_tau,m_mup, m_muz, m_mun;
D m_ni2,m_mup, m_muz, m_mun;
D m_n_n, m_p_n,m_p_p,m_n_p,m_v0,m_vd,m_y_junction;
IdxTy m_mid;
//MyBlock m_Na, m_Nd;
int m_u,m_z,m_p,m_n,m_v_nx,m_v_ny,m_v_px,m_v_py,m_q,m_trap;
int m_na,m_nd,m_x,m_y;
int m_n_state,m_n_fixed,m_layers;
bool m_save_sol;

Defects m_defects;
bool m_exit;
InitTy * m_init;
MyMesh * m_mesh;
IdxTy m_nodes;
DiriMap m_diri; // Dirichlet map from dof to value. 
BCMap m_bc_n,m_bc_p;
EsTy *m_es; // : public EquationSystems
NonLinTy *m_jr; // : public NonlinearImplicitSystem,
MeIdx m_idx;
D m_dx, m_dy;
IdxTy m_n_x,m_n_y;

//ZZDads m_dads;
}; // ficklanc




/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////

#ifdef  TEST_FICK__


static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
int main(int argc,char **args)
{
//int  m_ierr = PetscInitialize(&argc,&args,(char*)0,help);if (m_ierr) return m_ierr;
typedef mjm_fick_lanczos_2d Myt;
 // Initialize libraries, like in example 2.
    libMesh::LibMeshInit init (argc, args);

Myt x(argc,args,init);

//x.solve();
if (!x.exit()) x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

