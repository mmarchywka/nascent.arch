#ifndef MJM_FEM_INDEXING_H__
#define MJM_FEM_INDEXING_H__

//#include "mjm_fick_equations.h"
#include "mjm_globals.h"
//#include "mjm_qfle.h"
//#include "mjm_logic_base.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
//#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_diffuse_iterators.h"
//#include "mjm_defect_levels.h"
//#include "mjm_gauss.h"


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
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mjm_libmesh.h"
// fick 
#pragma GCC diagnostic ignored "-Wunused-parameter"


/*
2017-08-13 removed from mjm_fick_lanczos

*/


template<class Tbr, class Tparams>
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
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::MyIntBlock  MyIntBlock;
typedef typename Tr::MyuIntBlock  MyuIntBlock;
typedef typename Tr::MySparse MySparse;
typedef int intt;

//typedef fl_params FlpTy;
typedef Tparams  FlpTy;
typedef Tbr  BraTy;
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
//const D & _dy, const IdxTy padding)
const D & _dy, const IdxTy padding_x,const IdxTy padding_y)
: d(distance),nx(_nx),ny(_ny),dx(_dx),dy(_dy), layers_x(padding_x),layers_y(padding_y) {Init(); } 
// these are in mesh corrdinates
void center(const int cx, const int cy) { cnx=cx; cny=cy; } 
void reference(const int cx, const int cy) { refx=cx; refy=cy; } 
bool ok() { return m_ok; } 
void inc()
{
if (!ok()) return;
++posx;
if ( posx>nmaxx) 
{ 
	posx=nminx; 
	if (posy==nminy) { m_ok=false; }
	else --posy; // unsigned fick again 
}
}
// posx is relative to total grid, NOT just real grod
int centerx() const { return cnx; }
int centery() const { return cny; }
int srcx() const { return refx; }
int srcy() const { return refy; }

int relx() const { return int(posx)-int(cnxl); } 
int rely() const { return int(posy)-int(cnyl); } 
D delx() const { return dx*relx(); } 
D dely() const { return dy*rely(); } 
int meshx() const { return int(posx)-layers_x; }
int meshy() const { return int(posy)-layers_y; }
D r2() const { return relx()*relx()*dx*dx+rely()*rely()*dy*dy; } 
// these are all medh parameters and this will return the right geo node
// for on mesh nodes. However oor makes this wrong 
//int node() const { return (posx-m_layers)+nx*(posy-m_layers); } 
// this is a geo node starting with zero at padding  lower left.
//int node() const { return (posx)+nx*(posy); } 

void start() {
// zero on the mesh is layers,layers on the composite grid 
cnxl=cnx+layers_x;
cnyl=cny+layers_y;
// these are already padded with 2layers
nxl=nx; //+layers;
nyl=ny;//+layers;
// this is signed??? no
// these limits need to extend by layers on each edge
// this can cause int overflow for limiting values... doh 
D testdi=d/dx;
if (testdi>(1<<30)) testdi=(1<<30);
intt xx=intt(testdi+.5)+1;

//if (xx>cnxl) nminx=-layers; else nminx=cnx-xx;
if (xx>cnxl) nminx=0; else nminx=cnxl-xx;
nmaxx=cnxl+xx;
if (nmaxx>= nxl ) nmaxx=nxl-1;

testdi=d/dy;
if (testdi>(1<<30)) testdi=(1<<30);
intt yy=intt(testdi+.5)+1;
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
//void Init() { layers_x=layers; layers_y=layers; m_ok=false; } 
void Init() {  m_ok=false; } 
StrTy to_string() const
{
Ss ss;
ss<<MMPR4(nx,ny,cnx,cny);
ss<<MMPR2(cnxl,cnyl);
ss<<MMPR2(refx,refy);
ss<<MMPR4(nminx,nminy,nmaxx,nmaxy);
ss<<MMPR3(dx,dy,d);
ss<<MMPR3(posx,posy,m_ok);
ss<<MMPR2(layers_x,layers_y);

return ss.str();

}
int nx,ny,nminx,nminy,nmaxx,nmaxy;
int cnx,cny,refx,refy;
int cnxl,cnyl,nxl,nyl;
D dx,dy,d;
int posx,posy;
bool m_ok;
//int layers,layers_x,layers_y;
int layers_x,layers_y;
}; // round_itor

typedef round_itor NayItor;
class peripheral_nodes
{

public:


}; // peripheral_nodes

mesh_indexing_new(MyMesh & m, const FlpTy & p, const BraTy & b) 
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
	if (x< -m_layers_x) return bad();
	int y0=((node/m_n_x));
	int y=(y0+dy);
	if (y>=m_n_yl) return bad();
	if (y<-m_layers_y) return bad();
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
//if (bc ) ll+=2*m_layers;
if (bc ) ll+=2*m_layers_x;
int src_y=(node/ll);
int src_x=node-src_y*ll; // faster than mod? 
//int cx=(node%ll)+ int(sx*offx/m_dx+.5)*sx;
int cx=src_x+ int(sx*offx/m_dx+.5)*sx;
//int cy=(node/ll)+ int(sy*offy/m_dy+.5)*sy;
int cy=src_y+ int(sy*offy/m_dy+.5)*sy;
// now put back on grid if needed
//if (bc ) { cx-=m_layers; cy-=m_layers;}
if (bc ) { cx-=m_layers_x; cy-=m_layers_y;}
// the iterator however should use full sizes
//round_itor ri(r, m_n_x, m_n_y, m_dx,m_dy,m_layers);
//round_itor ri(r, m_n_x+2*m_layers_x, m_n_y+2*m_layers_y, m_dx,m_dy,m_layers);
round_itor ri(r, m_n_x+2*m_layers_x, m_n_y+2*m_layers_y, m_dx,m_dy,m_layers_x,m_layers_y);
ri.center( cx,cy);
ri.reference( src_x,src_y);
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
	const int rcbad=bad();
	int rc=rcbad;
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
	const int xzed=nx+m_layers_x;
	const int yzed=ny+m_layers_y;
	const int raw= xzed+yzed*(m_n_x+2*m_layers_x);
	if (ny<0) {rc=-1-raw;}
	else if (ny>=m_n_y) { rc=-1-raw+m_n_x*m_n_y; }
	else { rc=-1-raw+m_n_x*(ny); if (nx>0) rc+=m_n_x;}
	}
if (rc==rcbad) { MM_ERR(" bad value "<<MMPR4(ooc,nx,ny,m_n_y)) }
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
//int make_peripheral_nodes(const IdxTy nx, const IdxTy ny, const IdxTy layers,
int make_peripheral_nodes(const IdxTy nx, const IdxTy ny, const IdxTy layers_x,const IdxTy layers_y
,const IdxTy nstate, const IdxTy nfixed)
{
//	m_layers=layers;
	m_layers_x=layers_x;
	m_layers_y=layers_y;
	//const IdxTy l2=layers*2;
	const IdxTy l2x=m_layers_x*2;
	const IdxTy l2y=m_layers_y*2;
	// 1 -> 2(nx+ny) + 4, 2-> 4(nx+ny)+16
	const int size=(nx+l2x)*(ny+l2y)-nx*ny; //l2*((nx+ny)+ l2);
	m_index_base=-size;
	m_peri_x=m_n_x+l2x;
	m_peri_y=m_n_y+l2y;
	m_n_xl=m_n_x+m_layers_x;
	m_n_yl=m_n_y+m_layers_y;
//	m_peri_base=m_peri_x*m_layers;
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
for (int i=-m_layers_y; i<m_n_yl; ++i)
{
const bool inmesh=(i>=0)&&(i<m_n_y);
for (int j=-m_layers_x; j<m_n_xl; ++j)
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
BraTy m_branches;
D m_dx,m_dy,m_minx,m_miny;
int m_n_x, m_n_y;
int  m_n_xl, m_n_yl;
NayItor m_ri;
NodeIndex m_node_idx;
//int m_layers;
int m_layers_x,m_layers_y;
int m_index_base,m_peri_x,m_peri_y; //,m_peri_base; // minimum value of the peripheral pixels
MyBlock m_peri_state; // same as m_state except it is only the peripheral values that do not come from solver 
MyBlock m_peri_fixed; 
}; // mesh_indexing_new




////////////////////////////////////////////////////////////////
#endif
