#ifndef MJM_FICK_LANCZOS_2D_H__
#define MJM_FICK_LANCZOS_2D_H__

#include "mjm_flow_kernels.h"
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
#include "mjm_gauss.h"

#include "mjm_fem_indexing.h"
#include "mjm_convergence_tracking.h"
#include "mjm_fick_equations.h"

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
// fick 
#pragma GCC diagnostic ignored "-Wunused-parameter"





// the text has a comment ender lol 
#if 0
 2710  ./test_fick.tex -compilefick2d 
 2711  echo -e "set-param gr_tau 1e-12\nset-param mesh_nx 10\nset-tree enable_dump4 0\nset-param time_step 1e-13\nset-param vd .45\nset-param output_freq 1000\ninit\nbanner\nset-tree enable_flmb 1\nset-tree enable_ss 0\nsteps 100\ndump\nexit"  | ./a.out > xxx
 2712  cat xxx | grep fickiter | sed -e 's/.*=//' | sed -e 's/[() ] */ /g' > yyy
 2713  R
#endif

/*

2017-07-19 The bounndaries are a problem due to the stupid normalization to make diffusion
positive defininte. Since it is normalized, the number of dimensions matter and effectively
there is reflection from the missing elements. The terminals are probably ok but the 
sides need to have the diffusion normalized right or artificial edge effects are created.  
However, fixing the gross problems appears to have helped a lot .

2017-07-08 the quadrant thing sort of works with just one quad, can probably
rotate them as it is probably 16x faster than doing all of them. Tge solver 
seems to have confused indexing but may work when debuffed.

2017-07-06 The flow with impules function seems to sort of work although the
steady state solution made from self consistent u-np solution semes to 
have a memory corruption or node index problem. I'd guess the system matrix is
not indexed right. 

Start to make functor code for kernel and usage of transfer elements.


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
// genearated code could alos make setters to go wtih getters. 
//IdxTy output_freq() const { return m_map.get_uint("output_freq",100); } // // 100;
//D total_length() const { return m_map.get_double("total_length",.25e-6*2000); }
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
D time_step() const { return m_map.get_double("time_step",1e-13); }
IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
IdxTy mesh_ny() const { return m_map.get_uint("mesh_ny",100); } // // 100;
// cgs mostly 
D mesh_xmin() const { return m_map.get_double("mesh_xmin",0); }
D mesh_ymin() const { return m_map.get_double("mesh_ymin",0); }
D mesh_xmax() const { return m_map.get_double("mesh_xmax",1e-4); }
D mesh_ymax() const { return m_map.get_double("mesh_ymax",5e-4); }
int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
int bulk_layers_x() const { return m_map.get_int("bulk_layers_x",bulk_layers()); }
int bulk_layers_y() const { return m_map.get_int("bulk_layers_y",bulk_layers()); }
D bcf() const { return m_map.get_double("bcf",1.0); }
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
IdxTy gr_flags() const { return m_map.get_uint("gr_flags",0); }


D vsat() const { return m_map.get_double("vsat",1e7); } // 1e7 for Si

D Vt() const { return m_map.get_double("Vt",.0259); }
// silicon now 
D qe() const { return m_map.get_double("qe",1.6e-19/8.854e-14/11.68); } // later div 11.68 lol 
D ni2() const { return m_map.get_double("ni2",1.5e10*1.5e10); }
//D mun() const { return m_map.get_double("mun",1000); }
D mun() const { return m_map.get_double("mun",100); }
D mup() const { return m_map.get_double("mun",100); }
D muz() const { return m_map.get_double("mun",10); }
D f() const { return m_map.get_double("f",.5); }
//D f() const { return m_map.get_double("f",50.0); }
D Nd() const { return m_map.get_double("Nd",f()*1e15); }
D Na() const { return m_map.get_double("Na",f()*1e14); }
D v0() const { return m_map.get_double("v0",0); }
D vd() const { return m_map.get_double("vd",.45); }
D trap_speed() const { return m_map.get_double("trap_speed",1.0); }
D trap_N() const { return m_map.get_double("trap_N",1e13); }
D vmax() const { return m_map.get_double("vmax",1e0); }
D fl_cutoff() const { return m_map.get_double("fl_cutoff",1e-4); }
IdxTy quad_pattern() const { return m_map.get_uint("quad_pattern",0); }
D np_alpha() const { return m_map.get_double("np_alpha",.1); }
D test_vpx() const { return m_map.get_double("test_vpx",0); }
D test_vpy() const { return m_map.get_double("test_vpy",0); }
D test_vnx() const { return m_map.get_double("test_vnx",0); }
D test_vny() const { return m_map.get_double("test_vny",0); }
StrTy convergence_tracker() const { return m_map.get_string("convergence_tracker",""); }
StrTy convergence_file() const { return m_map.get_string("convergence_file",""); }
StrTy convergence_label() const { return m_map.get_string("convergence_label",""); }
D converge_tol() const { return m_map.get_double("converge_tol",5e-3); }
IdxTy converge_iters() const { return m_map.get_uint("converge_iters",100); }


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
ss<<"bulk_layers="<<bulk_layers()<<sep;
ss<<"bulk_layers_x="<<bulk_layers_x()<<sep;
ss<<"bulk_layers_y="<<bulk_layers_y()<<sep;
ss<<"bcf="<<bcf()<<sep;
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
ss<<"gr_flags="<<gr_flags()<<sep;
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
ss<<"fl_cutoff="<<fl_cutoff()<<sep;
ss<<"quad_pattern="<<quad_pattern()<<sep;
ss<<"np_alpha="<<np_alpha()<<sep;
ss<<"test_vpx="<<test_vpx()<<sep;
ss<<"test_vpy="<<test_vpy()<<sep;
ss<<"test_vnx="<<test_vnx()<<sep;
ss<<"test_vny="<<test_vny()<<sep;
ss<<"convergence_tracker="<<convergence_tracker()<<sep;
ss<<"convergence_file="<<convergence_file()<<sep;
ss<<"convergence_label="<<convergence_label()<<sep;

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
bool show_iter_status() const { return m_map.get_bool("show_iter_status",!true); }
bool show_iter_status_everyfew() const { return m_map.get_bool("show_iter_status_everyfew",true); }
bool always_dump_state() const { return m_map.get_bool("always_dump_state",!true); }
bool always_dump_to_file() const { return m_map.get_bool("always_dump_to_file",!true); }
bool dump_state_everyfew() const { return m_map.get_bool("dump_state_everyfew",true); }
bool dump_state_atend() const { return m_map.get_bool("dump_state_atend",true); }
bool dump_counter_everyfew() const { return m_map.get_bool("dump_counter_everyfew",true); }
bool cerr_counter_everyfew() const { return m_map.get_bool("cerr_counter_everyfew",true); }
bool debug_xfer() const { return m_map.get_bool("debug_xfer",!true); }
bool debug_xfer2() const { return m_map.get_bool("debug_xfer2",!true); }
bool enable_trapping() const { return m_map.get_bool("enable_trapping",!true); }
bool enable_potential() const { return m_map.get_bool("enable_potential",true); }
//bool enable_gr() const { return m_map.get_bool("enable_gr",!enable_ss()); }
bool enable_gr() const { return m_map.get_bool("enable_gr",true); }
bool debug_gr() const { return m_map.get_bool("debug_gr",!true); }
bool enable_flmb() const { return m_map.get_bool("enable_flmb",true); }
bool enable_overlap() const { return m_map.get_bool("enable_overlap",!true); }
bool enable_ss() const { return m_map.get_bool("enable_ss",!true); }
bool disable_kernle() const { return m_map.get_bool("disable_kernle",!true); }
bool enable_dump4() const { return m_map.get_bool("enable_dump4",!true); }
bool enable_convergence_track() const { return m_map.get_bool("enable_convergence_track",!true); }
bool zero_fd_field() const { return m_map.get_bool("zero_fd_field",!true); }
bool zero_right_field() const { return m_map.get_bool("zero_right_field",!true); }

StrTy to_string(const StrTy & sep=" ") const
{
// should be geneated code, do not edit  wtf ?
Ss ss;
ss<<"show_iter_status="<<show_iter_status()<<sep;
ss<<"show_iter_status_everyfew="<<show_iter_status_everyfew()<<sep;
ss<<"always_dump_state="<<always_dump_state()<<sep;
ss<<"always_dump_to_file"<<always_dump_to_file()<<sep;
ss<<"dump_state_everyfew="<<dump_state_everyfew()<<sep;
ss<<"dump_state_atend="<<dump_state_atend()<<sep;
ss<<"dump_counter_everyfew="<<dump_counter_everyfew()<<sep;
ss<<"cerr_counter_everyfew="<<cerr_counter_everyfew()<<sep;
ss<<"debug_xfer="<<debug_xfer()<<sep;
ss<<"debug_xfer2="<<debug_xfer2()<<sep;
ss<<"enable_trapping="<<enable_trapping()<<sep;
ss<<"enable_potential="<<enable_potential()<<sep;
ss<<"enable_gr="<<enable_gr()<<sep;
ss<<"debug_gr="<<debug_gr()<<sep;
ss<<"enable_flmb="<<enable_flmb()<<sep;
ss<<"enable_overlap="<<enable_overlap()<<sep;
ss<<"enable_ss="<<enable_ss()<<sep;
ss<<"disable_kernle="<<disable_kernle()<<sep;
ss<<"enable_dump4="<<enable_dump4()<<sep;
ss<<"enable_convergence_track="<<enable_convergence_track()<<sep;
ss<<"zero_fd_field="<<zero_fd_field()<<sep;
ss<<"zero_right_field="<<zero_right_field()<<sep;

//ss<<"vmax="<<vmax()<<sep;
return ss.str();
}

}; // fl_params
//////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////

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

// the hardcoded bersions had no benefit except maybe for speed but it was just
// extra clutter to maintain and the map is low overhead. 
// the thing needs a consistent means between config file, cmd line, and command interpretter
//typedef BranchesVar BraTy;
typedef fl_branches BraTy;
//typedef FickLancParamsVar FlpTy;
typedef fl_params FlpTy;




//typedef gauss_kernel2<fl_params,fl_branches> FlowKernel;
typedef gauss_kernel2<FlpTy,BraTy> FlowKernel;
//typedef mesh_indexing_new MeIdx; 
typedef mesh_indexing_new<fl_branches,fl_params>  MeIdx; 
typedef convergence_track_base ConvTrack;
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
class boundary_values
{
public:
D m_n_n, m_p_n,m_p_p,m_n_p,m_v0,m_vd; // ,m_y_junction;
DiriMap m_diri; // Dirichlet map from dof to value. 
BCMap m_bc_n,m_bc_p;
int m_layers_x,m_layers_y;

}; // boundary_values

public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
mjm_fick_lanczos_2d():m_old(0),m_save_sol(false),m_exit(false),m_init(0),m_mesh(0),m_nodes(3),m_es(0),m_es_np(0),m_jr(0),m_jr_np(0),m_del_np(false),m_ct(0) {Init();}
mjm_fick_lanczos_2d(const bool x ):m_old(0),m_save_sol(false),m_exit(false),m_init(0),m_mesh(0),m_nodes(3),m_es(0),m_es_np(0),m_jr(0),m_jr_np(0),m_del_np(false),m_ct(0) {if (x) Init();}

// make this m_points now, wtf
//mjm_fick_lanczos_2d(const IdxTy s):m_old(0), m_save_sol(false),m_exit(false),m_init(0),m_mesh(0),m_nodes(3),m_es(0),m_es_np(0),m_jr(0),m_jr_np(0),m_del_np(false) {Init();}
mjm_fick_lanczos_2d(int argc,char **args, InitTy & init) : m_old(0),m_save_sol(false),m_exit(false),m_init(0),m_mesh(0), m_nodes(3),m_es(0),m_es_np(0),m_jr(0),m_jr_np(0),m_del_np(false),m_ct(0)
{
m_init=&init;
//m_size=1;
//m_points=2000; // this really should not have a default.. 
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
const StrTy s=StrTy(args[i]);
if (s==StrTy("-interrupt")) debug_levels::interrupt(); 
m_tree.config("-tree",i,argc,args);
m_flp.config("-params",i,argc,args);
//configi(m_points,"-points",i,argc,args);
cmdcmd("-cmd",i,argc,args);
source("-source",i,argc,args);
m_flp.config_set("-set-param",  i,  argc, args);
m_tree.config_set("-set-branch",  i,  argc, args);
if (i==istart) { MM_ERR(" did nothing with i="<<i<< " " <<args[i]) {++i; } } 

}
if (!m_exit) Init();
}
~mjm_fick_lanczos_2d() { delete m_old;  delete m_mesh; delete m_es;  if (m_del_np) delete m_es_np;
delete m_ct;
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
if (false) MM_ERR(" processing "<<li.dump())
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="solve") { solve(); continue; } 
if (cmd=="source") { if (li.cmd_ok(3))  li.source(li.word(1));  continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
if (cmd=="potential") { step_potential(); 
		m_cm.inc("iterations");	
		m_cm.inc("potential-steps");	
		//per_step_outputs(iter);

 continue; } 
if (cmd=="flows") { if (li.cmd_ok(2)) { IdxTy iter;  step_flow_only(iter,li.n(1));  } continue; } 
if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
if (cmd=="iter-status") { show_iter_status();  continue; } 
if (cmd=="iter-status-cout") 
	{ if (li.cmd_ok(2)) show_iter_status(4,li.word(1)); 
		else show_iter_status(4,local_label);  continue; } 
//if (cmd=="dump-counter") { m_cm.dump("m_cm-dump "  , std::cout); continue; } 
if (cmd=="dump-counter") { dump_counter(StrTy("dump-counter-cmd")+StrTy(CRLF)); continue; } 
if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
if (cmd=="save") { if (li.cmd_ok(2)) save_to_file(li.word(1),local_label,1);  continue; } 
if (cmd=="load") { if (li.cmd_ok(2)) load_from_file(li.word(1),local_label,1);  continue; } 
if (cmd=="legend") { dump_legend(std::cout,local_label,1);  continue; } 
if (cmd=="dump-human") { dump_state_etc(std::cout,local_label,0);  continue; } 
if (cmd=="init") { Init(); /*  m_cm.clear(); */  continue; } 
if (cmd=="initmesh") { InitMesh(); /*  m_cm.clear(); */  continue; } 
// after zero out the state add back BC for testing 
if (cmd=="initbc") { InitBC();   continue; } 
if (cmd=="tracker") { InitConvergence();   continue; } 
if (cmd=="impulse") { initial_impulse(m_state);   continue; } 
if (cmd=="plates") { 
	if (li.cmd_ok(5))  initial_plates(li.wordf(1),li.wordf(2),li.wordf(3),li.wordi(4));  
	else MM_ERR(" not enough params for plates "<<li.dump() ) continue; } 
if (cmd=="image-pair") { initial_image(m_state);   continue; } 
if (cmd=="zero") { initial_zero(m_state);   continue; } 
if (cmd=="dump-fixed") { dump_fixed(std::cout,local_label,3); /*  m_cm.clear(); */  continue; } 
if (cmd=="save-solution") { m_save_sol=true;  continue; } 
if (cmd=="push-state") { InitPushOld();  continue; } 
if (cmd=="clear-cm") {  m_cm.clear();  continue; } 
if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="test") { test(li);  continue; } 
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(3))  refine(li.wordf(1),li.wordf(2));  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd=="exit") { clean_up(); m_exit=true; return; } 
if (cmd=="interrupt") { debug_levels::interrupt();  continue; } 
if (cmd=="dump-lut") { m_idx.dump_lut(local_label);  continue; } 

if (cmd.length()>1) if (cmd.c_str()[0]!='#')
		MM_ERR(" unrecignized command "<<li.dump())

}




} //command_mode
// when exiting the interpretter
void clean_up()
{


} 

template <class Ts> void dump_counter(const Ts & lbl )
{
 	m_cm.dump(lbl , std::cout); 
 	m_cm.dump(lbl , std::cerr); 
	MM_MSG(lbl<<CRLF<<m_cm.time_table("solver_times"))
	MM_ERR(lbl<<CRLF<<m_cm.time_table("solver_times"))
}


// tests, document what this crap should do lol

// the command object needs a shift function or something
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fick 
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

const bool enable_convergence_track=m_tree.enable_convergence_track();
IdxTy iter=0;
for (iter=0; iter<n; ++iter)
{
	m_cm.inc("iterations");	
	step(iter);
	per_step_outputs(iter);	
	//if (m_ct!=0) if (enable_convergence_track) if ((*m_ct).converged()) break ; 
	if (m_ct!=0) if (enable_convergence_track) if ((*m_ct).exit()) break ; 
} // iter
MM_MSG(" updating the drif tand charge info without updating u")
MM_ERR(" updating the drif tand charge info without updating u")
update_drift(m_state);

}


void refine(const D & xfac,const D & yfac )
{
// sinful old mesh right now 
const IdxTy nx=m_flp.mesh_nx();
const IdxTy ny=m_flp.mesh_ny();
// these should be rats 
const IdxTy newnx=nx*xfac;
const IdxTy newny=ny*yfac;
// now we have the new value, 
InitPushOld(); // this pushes the old solution /m
// reborn  all is new 
m_flp.set("mesh_nx",newnx);
m_flp.set("mesh_ny",newny);
const IdxTy nx2=m_flp.mesh_nx();
const IdxTy ny2=m_flp.mesh_ny();
MM_MSG(" regined sizes are "<<MMPR2(nx2,ny2))
// these should be rationals lol s 
const bool x=m_save_sol; // we are doing this here to avoid dpulicat init
m_save_sol=true; // we are doing this here to avoid dpulicat init
Init();
m_save_sol=x;
MM_ONCE(" preserving m_save_sol on refine just in case later init was supposed to purge",)
}


void config_banner()
{
MM_INC_MSG(m_cm,"test test" ) 
MM_MSG(" configuration "<<m_flp.to_string()<<" points "<<m_nodes)
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
void show_iter_status(const IdxTy flags=0,const StrTy & lbl=StrTy(" show_iter_status") ) const
{	
// 	const StrTy ival=(m_cm).format_entry("iterations");
 	const IdxTy ival=(m_cm).get_count("iterations");
if ((flags==0)||((flags&1)!=0))	
	{ MM_STATUS(" "<<ival<<" "<<m_state(m_nodes-m_n_x-m_n_x*.5+0,m_u)<<status_string())  }
if (((flags&2)!=0))
	{ MM_ERR(lbl<<" "<<ival<<" "<<m_state(m_nodes-m_n_x-m_n_x*.5+0,m_u)<<status_string())  } 
if (((flags&4)!=0))
	{ MM_MSG(lbl<<" "<<ival<<" "<<m_state(m_nodes-m_n_x-m_n_x*.5+0,m_u)<<status_string()) }
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
	if (m_tree.show_iter_status())  { show_iter_status(); } 
	else if (everyfew) if (m_tree.show_iter_status_everyfew())  { show_iter_status(7); } 
 	//{	MM_STATUS(" "<<ival<<" "<<m_state(m_points-2+0,m_u)<<status_string()) }
// 	{	MM_STATUS(" "<<ival<<" "<<m_state(m_nodes-m_n_x-m_n_x*.5+0,m_u)<<status_string()) }
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
D q=0;
D umax=-1e30;
D umin=1e30;

	for (IdxTy i=0; i<m_nodes; ++i)
		{
			D np=m_state(i,m_n)*m_state(i,m_p)/m_ni2;
			if (np>nimax) nimax=np;
			if (np<nimin) nimin=np;
			q+=total_q(i);
			D u=m_state(i,m_u);
			if (u>umax) umax=u;
			if (u<umin) umin=u;
		}
ss<<MMPR(nimax);
ss<<MMPR(nimin);
ss<<MMPR3(q,umax,umin);
ss<<"                     "; // guard, although the status macro should erase the line doh
#endif
return ss.str();

}
//void dump_legend() { dump_state_etc(std::cout,"fick",1); } 
template <class Os> void dump_legend(Os & os, const StrTy & label="legend",const IdxTy flags=0)
{

// the extraction code must include the iter info to R
//  and it needs to be called at the tight time or the iter count is fuiked. 
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
MM_ONCE(" qoi is ficked ",)
MM_MSG(" dump is from this location may have problems in combined output ")
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
MM_ONCE(" qoi is still ficked ",)
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


void save_to_file(const StrTy & fn , const StrTy & label, const IdxTy flags=0)
{
std::ofstream fout(fn);
MM_ERR(" saving to "<<fn);

Myt reload=*this; // hopefully this work without segfault lol
// as a nop this should still compare ok. lol 
reload.load_from_file(fn,label,flags); 
compare(*this,reload);
}
void load_from_file(const StrTy & fn , const StrTy & label, const IdxTy flags=0)
{
//std::ofstream fout(fn);
MM_ERR(" loadin from no real code  "<<fn);

}
bool  compare(Myt & m1, const Myt & m2)
{

return false; // skeleton alert 
} // compare


//////////////////////////////////////////////////////



/////////////////////////////////////////////////////
int  bad() const { return ~ int((~0U)>>1); } 
// move this down into real code 
void compute_qoi( MyBlock & qoi)
{
#if 1 
	const IdxTy nqoi=3;
	qoi.resize(m_nodes,nqoi);
	qoi.col_name("current",0);
	qoi.col_name("x",1);
	qoi.col_name("y",2);

	const IdxTy bad=m_idx.bad();
	const IdxTy bc=0;
	if (m_nodes<bc) return; 
	// really obsolete right now 
	const D Vt=m_flp.Vt(); // .0259;
	const MyBlock & s=m_state;
	const IdxTy df=1; const IdxTy db=1;
	const IdxTy sz=m_nodes-bc;
	for (IdxTy i=bc; i<sz; ++i)
	{

//		const D ih=1.0/(m_fixed(i+df,m_x)-m_fixed(i-db,m_x));

		const IdxTy north=m_idx.north(i);
		const IdxTy south=m_idx.south(i);
		const IdxTy east=m_idx.east(i);
		const IdxTy west=m_idx.west(i);

		const IdxTy top=(north==bad)?i:north;
		const IdxTy bot=(south==bad)?i:south;
		const IdxTy left=(west==bad)?i:west;
		const IdxTy right=(east==bad)?i:east;

		const D denx= (m_fixed(right,m_x)-m_fixed(left,m_x));
	//	const D ex=(denx==0)?0: ((s(right,m_u)-s(left,m_u))/(denx));
		const D deny=(m_fixed(top,m_y)-m_fixed(bot,m_y));
		const D ey= (deny==0)?0:((s(top,m_u)-s(bot,m_u))/(deny));
		const D ih=deny;
		const D gradU= (deny==0)?0:((s(top,m_u)-s(bot,m_u))/(deny));
		const D gradP= (deny==0)?0:((s(top,m_p)-s(bot,m_p))/(deny));
		const D gradN= (deny==0)?0:((s(top,m_n)-s(bot,m_n))/(deny));
		//const D gradU=ih*(s(i+df,m_u)-s(i-db,m_u));
		//const D gradP=ih*(s(i+df,m_p)-s(i-db,m_p));
		//const D gradN=ih*(s(i+df,m_n)-s(i-db,m_n));
		// this probably fails because the first 2 rows are in wrong arrangement,
		// check and skip 


		// this needs to use the calculated velocitiies 
		//const D phi_n=m_mun*s(i,m_n)*(gradU) -m_Vt*m_mun*gradN;
		const D phi_n=s(i,m_v_nx)*s(i,m_n) -Vt*m_mun*gradN;
		//const D phi_p=-m_mup*s(i,m_p)*(gradU)-m_Vt*m_mup*gradP;
		const D phi_p=s(i,m_v_px)*s(i,m_p)-Vt*m_mup*gradP;
		// this is charge current, the phis are particle current
		const D jc=phi_n-phi_p;
		qoi(i,0)=jc;
		qoi(i,1)=m_fixed(i,m_x); // no easy way to include location lol 
		qoi(i,2)=m_fixed(i,m_y); // no easy way to include location lol 
	} // i 
	// kluge
//	qoi(0,0)=qoi(1,0);
//	qoi(sz,0)=qoi(sz-1,0);
//	qoi(0,1)=m_fixed(0,m_x); // no easy way to include location lol 
//	qoi(sz,1)=m_fixed(sz,m_x); // no easy way to include location lol 

#endif

} // compute qoi 

void step_potential( ) // 
{
	m_cm.mark("find_pot");
	// this can be FD or FEM but needs diri BC. 
	solve_for_potential(m_state);
	update_drift(m_state); 
	m_cm.cumdel("pottotal","find_pot");
}
// largely for testing, remove all drift. 
// iter is a global count 
void step_flow_only(IdxTy &  iter, const IdxTy nsteps ) // const step_param & sp)
{
	const FlpTy & flp= m_flp;
	const D dt=flp.time_step(); 
	const D vpx=flp.test_vpx(); 
	const D vpy=flp.test_vpy(); 
	const D vnx=flp.test_vnx(); 
	const D vny=flp.test_vny(); 
	// px, py,nx,ny
	update_drift_zero(m_state,vpx,vpy,vnx,vny); 
	for (IdxTy i=0; i<nsteps; ++i)
	{
		m_cm.mark("mb");
		flmb_flow_kernle(dt,m_mun,m_n,m_v_nx,m_v_ny,flp,true);
		flmb_flow_kernle(dt,m_mup,m_p,m_v_px,m_v_py,flp,false);
		impose_bc();
		m_cm.cumdel("mbtotal","mb");
		// this needs to be ther other to get the stupif output right 
		m_cm.inc("iterations");	
		m_cm.inc("flow-iterations");	
		per_step_outputs(iter);
		++iter;
	}


} // step_flow_only 

void step(const IdxTy iter ) // const step_param & sp)
{
	// the ordering here is not physical but handy for testing... 
	const FlpTy & flp= m_flp;
	const bool enable_trapping=m_tree.enable_trapping();
	const bool enable_potential=m_tree.enable_potential();
	const bool enable_gr=m_tree.enable_gr();
	const bool enable_flmb=m_tree.enable_flmb();
	const bool enable_ss=m_tree.enable_ss();
	const bool disable_kernle=m_tree.disable_kernle();
	const bool enable_convergence_track=m_tree.enable_convergence_track();
	const bool sanity_check=!true;
	const bool recover_flmb=true;	
	const D dt=flp.time_step(); 
	bool failed=false;
	m_cm.mark("step");
	if (enable_potential ) { step_potential();
//		m_cm.mark("find_pot");
		// this can be FD or FEM but needs diri BC. 
//		solve_for_potential(m_state);
//		update_drift(m_state); 
//		m_cm.cumdel("pottotal","find_pot");
	}

	if (enable_ss)
	{
		m_cm.mark("mbss");
		failed|=!solve_for_np(dt,m_mun,m_n,m_v_nx,m_v_ny,m_mup,m_p,m_v_px,m_v_py,flp); 
		if (failed)  { MM_INC_MSG(m_cm, "ss_solver_fails" )}

		m_cm.cumdel("mbsstotal","mbss");
		impose_bc();
	}

	bool try_flmb=enable_flmb||(failed&&recover_flmb);
	if (try_flmb)
	{
		if (failed) { MM_MSG(" recogering ") } 
		failed=false;
		m_cm.mark("mb");
		if (disable_kernle){
		MM_MSG(" old flow disabled but should compile")
		//flmb_flow(dt,m_mun,m_n,m_v_nx,m_v_ny,flp);
		//flmb_flow(dt,m_mup,m_p,m_v_px,m_v_py,flp);
		} else {
		flmb_flow_kernle(dt,m_mun,m_n,m_v_nx,m_v_ny,flp,true);
		flmb_flow_kernle(dt,m_mup,m_p,m_v_px,m_v_py,flp,false);
		}


		m_cm.cumdel("mbtotal","mb");
		if (sanity_check) check_carriers("mbfl")	;
	
		impose_bc();
	}

	if (enable_gr)
	{
		m_cm.mark("find_gr");
		const IdxTy gr_flags=flp.gr_flags(); 
		const bool  debug_gr=m_tree.debug_gr(); 
		const bool do_old_gr=((gr_flags&1)==0);
		const bool do_int_gr=((gr_flags&2)!=0);
		if (debug_gr) { MM_MSG(" gr path "<<MMPR3(gr_flags,do_old_gr,do_int_gr))}
//		gr_attributing_sf(m_state, dt, m_n,m_p, m_points);
		if (do_old_gr)gr_point_by_point(m_state, dt, m_n,m_p, m_nodes);
//	MM_MSG(" integrating gr hardcoded")
		//MM_ONCE(" integrating gr hardcoded out now ",)

		if (do_int_gr) gr_point_by_point_integrating(m_state,  dt, m_n, m_p, m_nodes ,5,5);
		m_cm.cumdel("grtotal","find_gr");
	}
	if (sanity_check) check_carriers("gr")	;
//	enforce_bc(m_state); // ????? 
	if (enable_trapping){
	m_cm.mark("find_traps");
	const bool do_sf_trap_step= !true; 
	MM_ONCE(" 2d trap code has not been tested ",)
	if (do_sf_trap_step)
	{
// this itor gone 
	//	Di din(m_points,m_x,m_n,m_state,m_fixed);
//		Di dip(m_points,m_x,m_p,m_state,m_fixed);
//	 	m_defects.step( din, dip,m_fixed, m_state ,dt,m_u, m_x);
	} else { m_defects.step(m_fixed,m_state,dt,m_nodes,m_n,m_p,m_u,m_x); }
	m_defects.net_charge(m_state,m_trap,m_nodes,1);
	m_cm.cumdel("traptotal","find_traps");
} // enable_trapping
	impose_bc();
	m_cm.cumdel("steptotal","step");
	if (enable_convergence_track)
	{
		m_cm.mark("conv_tests");
		if (m_ct==0) { MM_ONCE(" convergence tracking on but tracker is null",) }
		else
		{
			(*m_ct).failed(failed);
			(*m_ct).add_iteration(m_state,m_fixed,m_nodes,&m_cm);
		}
		m_cm.cumdel("convtotal","conv_tests");
	} // enable_convergence_track
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

	const IdxTy sz=m_nodes;
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
	const IdxTy sz=m_nodes;
	for (IdxTy i=0; i<sz; ++i) 
	{ 
		if (Ck::denan(state(i,col) ,__FILE__,__LINE__," fick "))
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

D rho_fixed(const IdxTy node) const
{
return m_fixed(node,m_nd)-m_fixed(node, m_na)+m_state(node,m_trap);
}
D total_q(const IdxTy node) const
{
//MM_MSG(MMPR3(m_state(node,m_p),m_state(node,m_n),rho_fixed(node)))
return rho_fixed(node)+m_state(node,m_p)-m_state(node,m_n);
}

class id_dof_itor
{
typedef MeshCItor ItorTy;
typedef std::vector<IdxTy>  IdDofVec;
typedef std::vector<D>  ValVec;

public:
id_dof_itor(const MyMesh * mesh) :m_el((*mesh).active_local_elements_begin()),
    m_end_el((*mesh).active_local_elements_end()),m_id_dofs(4),m_vals(4) {Init(); }

typedef IdDofVec dof_vec_type;

	const ItorTy & el() const { return m_el; }
	bool ok() const { return m_ok; } 
	void inc() { ++m_el; Ok(); if (ok()) Update();  }
	//const IdxTy * id_dofs() const { return &m_id_dofs[0];} 
	const IdDofVec & id_dofs() const { return m_id_dofs;} 
	const D dx() const { return m_dx; } 
	const D dy() const { return m_dy; } 

 	const D interpolate(MyBlock & sol,const IdxTy in,const D & fx,const D & fy)
	{
		const IdDofVec & dv=id_dofs();
		const D val=sol(dv[0],in)*(1.0-fx)*(1.0-fy)
					+sol(dv[1],in)*(fx)*(1.0-fy)
					+sol(dv[2],in)*(fx)*(fy)
					+sol(dv[3],in)*(1.0-fx)*(fy);

		return val;
	}
	// user should know if cached values are vald
 	const D interpolate(const D & fx,const D & fy)
	{
		const D val=m_vals[0]*(1.0-fx)*(1.0-fy)
					+m_vals[1]*(fx)*(1.0-fy)
					+m_vals[2]*(fx)*(fy)
					+m_vals[3]*(1.0-fx)*(fy);

		return val;
	}


	void vals(MyBlock & sol, const IdxTy in)
	{
		const IdDofVec & dv=id_dofs();
		for (IdxTy i=0; i<4; ++i) { m_vals[0]=sol(dv[0],in); }

	}
	void allocate(MyBlock & diffn,const D & dn,const D & fx,const D & fy)
	{
		const IdDofVec & dv=id_dofs();
		diffn(dv[0])+=dn*(1.0-fx)*(1.0-fy);
		diffn(dv[1])+=dn*(fx)*(1.0-fy);
		diffn(dv[2])+=dn*(fx)*(fy);
		diffn(dv[3])+=dn*(1.0-fx)*(fy);
	}
	void allocate_log(MyBlock & diffn,const D & dn,const D & fx,const D & fy)
	{
		const IdDofVec & dv=id_dofs();
		diffn(dv[0])+=dn*(1.0-fx)*(1.0-fy);
		diffn(dv[1])+=dn*(fx)*(1.0-fy);
		diffn(dv[2])+=dn*(fx)*(fy);
		diffn(dv[3])+=dn*(1.0-fx)*(fy);
	}

	void allocate_frac(MyBlock & diffn,const D & dn,const D & fx,const D & fy)
	{
		const IdDofVec & dv=id_dofs();
		diffn(dv[0])+=dn*sf(0,fx,fy); // 1.0-fx)*(1.0-fy);
		diffn(dv[1])+=dn*(fx)*(1.0-fy);
		diffn(dv[2])+=dn*(fx)*(fy);
		diffn(dv[3])+=dn*(1.0-fx)*(fy);
	}

	D sf(const IdxTy i, const D & fx, const D & fy)
	{
		switch (i)
		{
			case 0:{return (1.0-fx)*(1.0-fy); break ; } 
			case 1:{return (fx)*(1.0-fy); break ; } 
			case 2:{return (fx)*(fy); break ; } 
			case 3:{return (1.0-fx)*(fy); break ; } 
		} // i 
		MM_ERR(" bad shape "<<MMPR3(i,fx,fy))
	return 0; 
	}


private:
void Init()
{
Ok();  if (ok()) Update(); 
}
void Ok() { m_ok=(m_el!=m_end_el); } 

void Update()
{
		const Elem * elem=*el();
		const Node & ll=*((*elem).node_ptr(0));
		const Node & lr=*((*elem).node_ptr(1));
		const Node & ur=*((*elem).node_ptr(2));
		const Node & ul=*((*elem).node_ptr(3));
		m_id_dofs[0]=ll.id();
		m_id_dofs[1]=lr.id();
		m_id_dofs[2]=ur.id();
		m_id_dofs[3]=ul.id();
		// these are NEGATVE
		m_dx= -(ll(0)-ur(0));
		m_dy= -(ll(1)-ur(1));
}

ItorTy m_el, m_end_el;
bool m_ok; 
//IdxTy m_id_dofs[4];
IdDofVec  m_id_dofs;
D m_dx,m_dy;
ValVec m_vals;
}; // id_dof_itor;

// passing state ignores the rhs and boundaries. 
// this is really doing FEM and the const matrix could be saved 
bool solve_for_potential(MyBlock & state)
{
	mjm_special_qrule mats;
	MyBlock b(m_nodes);
	const IdxTy nnode=4; // vertexes per node 
	const IdxTy maxiters=5000;
	const D rtol=1e-8;
	MySparse local_mat;
	const bool fix_voltage=false;	
	const D qe=m_flp.qe(); ///11.68; // 1.6e-19/8.854e-14;
	const bool some_e_bounds_nodes = !true;
	const bool do_field_eqn = !some_e_bounds_nodes;
	const bool zero_fd_field=m_tree.zero_fd_field(); // &&!some_e_bounds_nodes; // !false;	 // until can lift FEM code for neuman 
	const bool zero_fem_field=false;	 // until can lift FEM code for neuman 
	const bool zero_right_field=m_tree.zero_right_field(); //false;	 // until can lift FEM code for neuman 
	// this also depends on the sizes not changing 
	const bool use_cached=false; // fix_voltage&&!zero_fem_field; // false;	
	bool update_mat= (m_pot_mat.size()==0)||!use_cached;
	MySparse & A=use_cached?m_pot_mat:local_mat;
	MyBlock mat; // FE coefficients for laplacian 
	mat.resize(nnode,nnode); // doh, surprised it ran at all wtf 
	// this is needed to insure that the computed values are not stale
	const D vd=m_flp.vd();
	if (m_bv.m_vd!=vd)
	{
		MM_MSG(" voltage changed "<<MMPR2(vd,m_bv.m_vd))
		m_bv.m_vd=vd;
		InitBC(); // m_diri is the potential bc although n/p similar 
	}

	DofVec dof_indices;
	const IdxTy nvar=1; // enable_ss?3:1;
	const IdxTy nv=nnode*nvar;
 	dof_info di =dof_info((*m_jr).get_dof_map(),nvar); 	
	IdxTy serial=0;
	id_dof_itor dof_itor(m_mesh);
	while (dof_itor.ok())
    {
		const Elem * elem=* dof_itor.el(); // *el;
		di.re_init(elem);  
		di.dof_map.dof_indices(elem, dof_indices, 0);
		if (update_mat)
		{ // these are positive as negative sign from inte by parts
			// is included. Nomralized such that  int sf = .25*a*b 
			bool normal_element=true;
			if (some_e_bounds_nodes )
			{
				IdxTy bno=0;
				D val=0;
				std::vector<IdxTy> dum;
				for (IdxTy i=0; i<nnode; ++i)
				{
					const IdxTy dofick=dof_itor.id_dofs()[i]; // 	
					// this should save the find result and use that 
					const bool is_bc_node= 
						(m_bv.m_diri.find(dofick)!=m_bv.m_diri.end());
					if (is_bc_node)
					{ val=m_bv.m_diri[dofick]; dum.push_back(i);  bno|=(1<<i); }
					//else bno<<=1; 
				}

				const bool zf=(bno!=0)&&(((val==0)&&!zero_right_field) 
					||((val!=0)&&zero_right_field)); 
				if (zf ) {
					mats.der_mat_y(mat,dof_itor.dx(),dof_itor.dy());
					di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
					mats.lap_mat(mat,dof_itor.dx(),dof_itor.dy());
				 	for (IdxTy i=0; i<dum.size(); ++i) mat.zero_row(dum[i]); 	
					di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
					mats.lap_mat_y(mat,dof_itor.dx(),dof_itor.dy());
				 	for (IdxTy i=0; i<dum.size(); ++i) mat.zero_row(dum[i]); 	
					di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
					normal_element=false;
				}
			}
			if (normal_element) {
			mats.lap_mat(mat,dof_itor.dx(),dof_itor.dy());
			di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
			mats.lap_mat_y(mat,dof_itor.dx(),dof_itor.dy());
			di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
			}
		}
		// needs to be integrated, scaled, and mult b q/eps
		const D afac=.25*dof_itor.dx()*dof_itor.dy();
		if (!false ){ for (IdxTy i=0; i<nnode; ++i)
		{
			const IdxTy dofick=dof_itor.id_dofs()[i]; // dof_indices[i];	
			//const D rhoqe=-m_qe*total_q(dofick)*afac;
			const D rhoqe=-qe*total_q(dofick)*afac;
			b(dofick)+=rhoqe;
		}}	
		++serial;
		dof_itor.inc();
	} // while 
	if (fix_voltage) { impose_diri( &A,&b,1.0); } 
	if (zero_fd_field) { impose_zero_field( &A,&b,zero_right_field,do_field_eqn); } 
	(*m_jr).set_system(A,b);
	(*m_jr).mjm_linear_solve(false,rtol,maxiters);
	const int cr=(*m_jr).get_converged_reason();
	const StrTy cs=(*m_jr).converged_string(cr);
	const bool div=(*m_jr).diverged();
//	if (div ) {     MM_MSG(" potential solve failed  see log for count " <<MMPR3(cr,cs,div)) }
	if (div ) {     MM_ONCE(" potential solve failed  see log for count " <<MMPR3(cr,cs,div),) }
	if (div)  { MM_INC_MSG(m_cm, StrTy("potential_solver_fails:")+cs )}
	if (!div)  { MM_INC_MSG(m_cm, "potential_solver_ok" )}
	// go ahead and use it anyway 
	(*m_jr).copy_solution(state,m_u);
	return !div;
}
#if 0

void solve_for_potential_old(MyBlock & state)
{
	mjm_special_qrule mats;
	MyBlock b(m_nodes);
//	b.resize(m_nodes);
	bool update_mat= (m_pot_mat.size()==0);
	MySparse &  A=m_pot_mat;
	MyBlock mat; // FE coefficients for laplacian 
	mat.resize(4,4); // doh, surprised it ran at all wtf 
	// iterate over the node, needing nearest neighbors to define del^2
	// need fixed charge and Dirichlet BC's.
 	// or just do this like real FEM and integrate over elements... 
	DofVec dof_indices;
//	std::vector<IdxTy> id_dofs(4);
//    (*m_es).print_info();
	//const bool enable_ss=m_tree.enable_ss();
	const IdxTy nvar=1; // enable_ss?3:1;
	const IdxTy nnode=4;
	const IdxTy nv=nnode*nvar;
 	dof_info di =dof_info((*m_jr).get_dof_map(),nvar); 	
	IdxTy serial=0;
//    dof_map.dof_indices (elem, dof_indices);
//    for (IdxTy var=0; var<nvar; var++)
	//dof_map.dof_indices (elem, dof_indices_var[var], var);
	// mesh element itor, get the laplacian code and JxW etc. 
	id_dof_itor dof_itor(m_mesh);

//	MeshCItor      el     =(*m_mesh).active_local_elements_begin();
//    const MeshCItor end_el = (*m_mesh).active_local_elements_end();
//    for ( ; el != end_el ; ++el)
	while (dof_itor.ok())
    {
		//const Elem * elem=*el;
		const Elem * elem=* dof_itor.el(); // *el;
		di.re_init(elem); // FUDD 

/*
		const Node & ll=*((*elem).node_ptr(0));
		const Node & lr=*((*elem).node_ptr(1));
		const Node & ur=*((*elem).node_ptr(2));
		const Node & ul=*((*elem).node_ptr(3));
		id_dofs[0]=ll.id();
		id_dofs[1]=lr.id();
		id_dofs[2]=ur.id();
		id_dofs[3]=ul.id();

//	MM_MSG(" shot "<<shot.n_dofs(0,0)<<" "<<shot.id())
		//const D dx= (*(*elem).node_ptr(0))(0)-(*(*elem).node_ptr(2))(0);
		const D dx= ll(0)-ur(0);
		const D dy= ll(1)-ur(1);

*/
		//const D dy= (*(*elem).node_ptr(0))(1)-(*(*elem).node_ptr(2))(1);
//MM_MSG(" fick "<<MMPR2(dx,dy))
		//const D b= (*elem).node_ptr(0)->(1)-(*elem).node_ptr(2)->(1);
//		di.dof_map.dof_indices(elem, dof_indices);
		di.dof_map.dof_indices(elem, dof_indices, 0);
//MM_MSG(" fick ")
		// this is incredibly dumb, most of a and b are the same, just do this once 
		// or better this is const for the const mesh doh 
		if (update_mat)
		{
			//mats.lap_mat(mat,dx,dy);
			mats.lap_mat(mat,dof_itor.dx(),dof_itor.dy());
//		di.sum_dof_mats(A, mat, dof_indices, nv);
			//di.sum_dof_mats(A, mat, id_dofs, nv);
			di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
			//mats.lap_mat_y(mat,dx,dy);
			mats.lap_mat_y(mat,dof_itor.dx(),dof_itor.dy());
		//di.sum_dof_mats(A, mat, dof_indices, nv);
			di.sum_dof_mats(A, mat, dof_itor.id_dofs(), nv);
		}
//MM_MSG(" fick ")
		// needs to be integrated, scaled, and mult b q/eps
		const D afac=.25*dof_itor.dx()*dof_itor.dy();
		if (!false ){ for (IdxTy i=0; i<nnode; ++i)
		{
//			const Node & no=*((*elem).node_ptr(i));
		//	const IdxTy dofick=id_dofs[i]; // dof_indices[i];	
			const IdxTy dofick=dof_itor.id_dofs()[i]; // dof_indices[i];	
			//const D rhoqe=-m_qe*total_q(dofick)*.25*dx*dy;
			const D rhoqe=-m_qe*total_q(dofick)*afac;
	//		MM_MSG(MMPR4(serial,i,dofick,rhoqe) )
			// this is wrong except it is FEM not FD, NOT + but just = 
			b(dofick)+=rhoqe;
//			b(dofick)=rhoqe;
		}}	
	++serial;
	dof_itor.inc();
	}
//} // m_pot_mat.size()==0
//for (IdxTy i=0; i<m_nodes; ++i ) 
// rhs is the integrated net charges

// get dirichlets 
	impose_diri( &m_pot_mat,&b,1.0); // this should use "A" 
	impose_bc(); // should be elsewhere. 
	(*m_jr).set_system(m_pot_mat,b); // this should use "A" 
	(*m_jr).mjm_linear_solve();
	(*m_jr).copy_solution(state,m_u);

}
#endif

void impose_zero_field( MySparse * A, MyBlock * rhs, const bool zero_right_field , const bool do_field_eqn)
{
// moved to bc code 
//const D vd=m_flp.vd();
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(i,true)]=0;
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;
const D scale=1e4;
MM_ONCE(" changing zero field scale to "<<scale,)
const IdxTy bad=m_idx.bad();
auto di=m_bv.m_diri.begin();
auto de=m_bv.m_diri.end();
while (di!=de)
{
const IdxTy dof=(*di).first;
const IdxTy north=m_idx.north(dof);
const IdxTy south=m_idx.south(dof);
IdxTy dofe=dof;
if (north==bad) dofe=south; else dofe=north;
const D val=(*di).second;
if (A) { // A->zero_row(dof); 
const bool zf=((val==0)&&!zero_right_field) ||((val!=0)&&zero_right_field); 
// there is a terrible offset lol 
if (zf) { 
if (do_field_eqn) { A->zero_row(dof); 

(*A)(dof,dof)=scale; (*A)(dof,dofe)=-scale; }}
else{
 A->zero_row(dof); 
(*A)(dof,dof)=scale; 
(*A)(dof,dofe)=0; // hopefully get petsc to open up sparsity pattern for later switch 
//A->zero_row(dofe); 
//(*A)(dofe,dofe)=scale; 
//(*A)(dofe,dof)=-scale; 
//(*A)(dof,dofe)=-scale; 
//if (rhs) { (*rhs)(dofe)=0; } 

}
 }
if (rhs) { (*rhs)(dof)=0; } 
++di;
}
}



void impose_diri( MySparse * A, MyBlock * rhs , const D & scale)
{
// moved to bc code 
//const D vd=m_flp.vd();
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(i,true)]=0;
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;

auto di=m_bv.m_diri.begin();
auto de=m_bv.m_diri.end();
while (di!=de)
{
const IdxTy dof=(*di).first;
const D val=(*di).second;
if (A) { A->zero_row(dof); (*A)(dof,dof)=scale; }
if (rhs) { (*rhs)(dof)=val*scale; } 
++di;
}
}




// this is used by the np solver and the dof calculation is  a question 
template <class Tm> void impose_diri( MySparse * A, MyBlock * rhs , Tm & diri, const IdxTy pitch, const IdxTy offset, const D & scale)
{
// moved to bc code 
//const D vd=m_flp.vd();
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(i,true)]=0;
//for ( IdxTy i=0; i<m_n_x; ++i) m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;

auto di=diri.begin();
auto de=diri.end();
while (di!=de)
{
const IdxTy dof=offset+pitch*(*di).first;
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

void update_drift_zero(MyBlock & state,const D & vpx=0, const D & vpy=0, const D & vnx=9, 
const D & vny=0)
{
    const IdxTy iinit=0;
    const IdxTy sz=m_nodes;
    for (IdxTy i=iinit; i<sz;  ++i)
	{
		state(i,m_v_px)=vpx;
		state(i,m_v_py)=vpy;
		state(i,m_v_nx)=vnx;
		state(i,m_v_ny)=vny;
		// usually irrelevant but fwiw
        state(i,m_q)=total_q(i);  
	} 
}

void update_drift(MyBlock & state)
{
    // this is bad with zero points, fwiw
//    state(0,m_u)=m_v0; state(sz,m_u)=m_vd; state(1,m_u)=m_v0;
    const D vsat=m_flp.vsat(); // 1e5;
    const IdxTy iu=m_u;
    const IdxTy sz=m_nodes;
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

#if 0
// redundant ( ie dx, dy, r2 ) param sets are to save computations 
void kernel_gauss(D * gni, D * gnorm, const D & beta, const D & r2, const D & v_x_t,
const D & dx, const D & dy, const D & vx, const D & vy, const D & t,
const D & n0, const D & ni )
{// this means the first part mus be calculated twice, fick 
	const D q2=(-beta*(r2+v_x_t*v_x_t));
	const D vexp=exp(q2-2.0*beta*(vx*t*dx+vy*t*dy));
if (gnorm!=0)
{
	const D  gauss=vexp; // exp(q2);
//	if (gnorm!=0) 
	*gnorm+=gauss;
}
	if (gni==0) return; 
	// wt  is wrong with this fick er 
//	const D vexp=exp(q2+2.0*beta*(vx*t*dx+vy*t*dy));
	// well duh, need to have outbound since that is the relevant v 
	*gni+=(n0*vexp); // -ni/vexp);

}
#endif

/*
Destination functors. These right now either peform the flows on an update
matrix to be sent back to m_state or matrix eleemnts to find the
steady state solution self consistently with u. 

*/

class flmb_func : public MyBlock
{
public:
typedef MyBlock Super;
flmb_func(const IdxTy n): Super(n) {}
// add to an equation for n and p 
// if from_id<0, this is an rhs term but in thise case they are all rhs, 
void gni(const int to_id, const int from_id, const D & coef_to, const D & coef_from, const
D & nto, const D & nfrom)
{
if (to_id>=0) (*this)(to_id)+=coef_from*nfrom;

}
}; // flmb_func


class flmb_ss 
{ // may be close but not quite, matrix indexes probabably boteched 
public:
// this needs to know the PITCH as it interleaves n and p dofs
// ficking factor of ficking 2 
flmb_ss(const IdxTy n): b(n<<1),m_odd(true) {}
void odd() { m_odd=true; }
void even() { m_odd=!true; }
// add to an equation for n and p 
// if from_id<0, this is an rhs term but in thise case they are all rhs, 
void gni(const int  to_id, const int from_id, const D & coef_to, const D & coef_from, const
D & nto, const D & nfrom)
{
const int off=(m_odd?1:0);
const int i=to_id*2+off;
const int j=from_id*2+off;
//if (i==j) return; // do not count self to self for now 
const IdxTy ty=((to_id<0)?1:0)+((from_id<0)?2:0);
switch (ty)
{
case 0:  // xfer from "from" to "to
{
// const D x=coef_from*nfrom; 
// this is dni / dt partial controbution from j 
//A(i,j)+=coef_from; 
if( j!=i) A(i,j)+=coef_from; 
else { A(j,j)=  coef_from-1; }
//A(j,j)=1;
//if (j!=i) 
// c1 
//{ A(j,j)+=  coef_from; }
//A(i,i)+= -  coef_from;
// note that this makes no contrib to b
break;
}
case 1: // outflow to the boundary
{
// a real loss occurs
//c1 //A(j,j)+=-  coef_from; 
//A(j,j)+=  coef_from; 
// but this is no incoming exitmenet 
break;
}
case 2: // inflow from boundary
{
b(i)-=coef_from*nfrom; 

break;
}
case 3: // irrelevant 
{

break; 
}

} ;// switch 
}

MySparse A;
MyBlock b;
bool m_odd;
}; // flmb_ss

/* 
 This is supposed to find the n and p steady states that correspeond to u
utlimately solving self consistently together. It seems it could work but
also patterns in output suggest corruption or wrong indexes.
Los mixing coefficient is needed too - probably better to just use the
nonlinear solver and write the derivatives. 
*/

bool  solve_for_np(const D & dt,const D & mun, const IdxTy in, const IdxTy ivxn, const IdxTy ivyn
, const D & mup,const IdxTy ip, const IdxTy ivxp, const IdxTy ivyp
,const FlpTy & flp)
{
	const bool show_rhs=false;
	const D alpha=flp.np_alpha(); // .01;
// except for the boundaries, n==p==0 is the solution.
// adding gr would help too 
	flmb_ss Ab(m_nodes);
//	MM_MSG(" start solve_for_np")
	Ab.even();
//	MM_ERR(" gauss kernel used")
	//gauss_kernel kernp(m_flp,dt,mup,m_dx,m_dy,m_tree);
	//gauss_kernel2 kernp(m_flp,dt,mup,m_dx,m_dy,m_tree);
	FlowKernel kernp(m_flp,dt,mup,m_dx,m_dy,m_tree);
	kernp.set_bc(m_bv.m_p_p,m_bv.m_p_n);
	//flmb_do(Ab,  dt, mup, ip, ivxp, ivyp, flp);
	flmb_do(Ab,  kernp,dt, mup, ip, ivxp, ivyp, flp);
	Ab.odd();
//	MM_MSG("  ODD start solve_for_np")
	//gauss_kernel kernn(m_flp,dt,mun,m_dx,m_dy,m_tree);
	//gauss_kernel2 kernn(m_flp,dt,mun,m_dx,m_dy,m_tree);
	FlowKernel  kernn(m_flp,dt,mun,m_dx,m_dy,m_tree);
	kernn.set_bc(m_bv.m_n_p,m_bv.m_n_n);
	//flmb_do(Ab,  dt, mun, in, ivxn, ivyn, flp);
	flmb_do(Ab,  kernn,dt, mun, in, ivxn, ivyn, flp);
	if (show_rhs) {
	for (IdxTy i=0; i<2*m_nodes; ++i)
	{
		if (false) if (Ab.b(i)!=0) { MM_MSG(" excitation "<<MMPR2(i,Ab.b(i))); } 
	}}
	//MM_MSG("  impose  start solve_for_np")
 	impose_diri( &Ab.A, &Ab.b ,m_bv.m_bc_n,2, 1, 1);
 	impose_diri( &Ab.A, &Ab.b ,m_bv.m_bc_p,2, 0, 1);
	//MM_MSG("    gr start solve_for_np")
	// this uses the linear approx which is not likely to help
	// for equillibrium, it is easy enough to find d for (n-d)*(p-d)= ni^2
	if (false) gr_point_by_point(Ab.A,Ab.b,.5*dt,in,ip,m_nodes);
	// this depletes everything 
	if (!false) gr_point_by_point_equ(Ab.A,Ab.b,.5*dt,in,ip,m_nodes);
	//MM_MSG(" call petsc  solve_for_np")
	(*m_jr_np).set_system(Ab.A,Ab.b,true);
	(*m_jr_np).mjm_linear_solve();
	const int cr=(*m_jr_np).get_converged_reason();
	const StrTy cs=(*m_jr_np).converged_string(cr);
	const bool di=(*m_jr_np).diverged();
	if (di ) { MM_MSG(MMPR3(cr,cs,di)) }
	if (di)  { MM_INC_MSG(m_cm, StrTy("np_solver_fails:")+cs )}
	if (!di)  { MM_INC_MSG(m_cm, StrTy("np_solver_ok:")+cs )}
	// m_n MUST = m_p+1, also do linear combs here 
	if (!di||(cr== -3 )) {
		MM_MSG( " taking best fixed point soln")
		 (*m_jr_np).blend_solution(m_state,m_p,alpha,1,true);
		return true;
		}
	else  {MM_MSG(" skipping solution due to divergence" ) }
	//MM_MSG(" return frm  petsc  solve_for_np")
	return !di;
}
#if 0
//template <class Ty> 
void flmb_flow(const D & dt,const D & mu, const IdxTy in, const IdxTy ivx, const IdxTy ivy,const FlpTy & flp)
{
	MM_ERR(" OLD CODE NO  kernel used set to ZERO ")
	flmb_func next(m_nodes);
//	flmb_do(next,  dt, mu, in, ivx, ivy, flp);
	for (IdxTy i=0; i<m_nodes; ++i)  m_state(i,in)=next(i);
}
#endif

void flmb_flow_kernle(const D & dt,const D & mu, const IdxTy in, const IdxTy ivx, const IdxTy ivy,const FlpTy & flp,const bool electrons)
{
	//MM_ERR(" ext kernel  gauss kernel used")
	flmb_func next(m_nodes);
	//gauss_kernel kernle(m_flp,dt,mu,m_dx,m_dy,m_tree);
	//gauss_kernel2 kernle(m_flp,dt,mu,m_dx,m_dy,m_tree,&m_cm);
	FlowKernel kernle(m_flp,dt,mu,m_dx,m_dy,m_tree,&m_cm);
	if( electrons) kernle.set_bc(m_bv.m_p_p,m_bv.m_p_n);
	else kernle.set_bc(m_bv.m_n_p,m_bv.m_n_n);
//	kernle.set(&m_cm); // this needs in ctor call to check radius size 
	flmb_do(next, kernle, dt, mu, in, ivx, ivy, flp);
	for (IdxTy i=0; i<m_nodes; ++i)  m_state(i,in)=next(i);
}



/* obsolete code that may have worked but corrpt not sure 
*/

/*
// new code with objexts for input and output
This iterates over all mesh nodes and exchanges flux with mesh nodes
or bulk nodes then iterates over bulk nodes flowing back to mesh.
*/
// find index for any node real mesh or bulk  this index is based 
D composite_state(const int node, const IdxTy col)
{
return(node<0)? m_idx.m_peri_state(-node-1,col):m_state(node,col);	
}
// this is dangerout as it lets user update the off mesh nodes 
D & composite_state2(const int node, const IdxTy col)
{
return(node<0)? m_idx.m_peri_state(-node-1,col):m_state(node,col);	
}

template <class Ty,class Tx> 
void flmb_do(Ty & dest, Tx & kernel, const D & dt,const D & mu, const IdxTy in, const IdxTy ivx, const IdxTy ivy,const FlpTy & flp)
{
 //   { MM_ERR(" doing flmb ") }
	const bool debug_xfer=m_tree.debug_xfer();
	const bool debug_xfer2=m_tree.debug_xfer2();
	const bool do_bc=!true;
	const bool do_fixed_bulk=true;
	typedef MeIdx::hood_iterator Hi; 
	IdxTy serial=0;
	const NodeItor ie=(*m_mesh).nodes_end();
	for (NodeItor ii=(*m_mesh).nodes_begin(); ii!=ie; ++ii)
	{
		const Node & np=*(*ii);
		const IdxTy idi=np.id();
		const int idigeo=m_idx.lut(idi,false);
		const D n0=m_state(idi,in);
		if (n0==0) { // this is now possible with impulse distro
			MM_ONCE(" carrier counts zero see log ",)
		//	MM_MSG(" n0 is zero "<<MMPR4(n0,idi,in,idigeo))
  continue; }	
		const D vx=m_state(idi,ivx); // this node it the src node
		const D vy=m_state(idi,ivy); 
		kernel.point(idi,n0,vx,vy,m_fixed(idi,m_x),m_fixed(idi,m_y));
		// the iterator is based on geo coords
		Hi hi=m_idx.hood_itor(kernel.radius,idigeo,kernel.tvx,kernel.tvy);
		// now loop over all nodes that can communicate and make the changes
		// need to keep track of stuff for normalization 		
		hi.start();
		while (hi.ok())
		{	// the itor is over geographical nodes 
			const int meshx=hi.meshx();
			const int meshy=hi.meshy();
			const int  idjgeo=m_idx.index(meshx,meshy); // hi.node();
				if (false) if ((meshx<0)||(meshy<0)) { 
				const IdxTy oorcode=m_idx.oor(meshx,meshy);

				MM_MSG(" oorcode "<<MMPR4(idjgeo,oorcode,meshx,meshy)) } 
			// off-mesh nodes are legitimate at the +/- y ends right now. 
			// the walk off increases norm and hence amount transfered to self later. 
//			int idj=0;
			const int idj=m_idx.lut(idjgeo,true);
			if (idj<0) 
			{ 
				const int  oorcode=m_idx.oor(meshx,meshy);
				if (false) MM_MSG(" oorcode "<<MMPR4(idjgeo,oorcode,meshx,meshy))
				// the oor layers at the y ends are the bulk media
				if (( oorcode&12)==0){ hi.inc(); continue; }
				// anything with wrong x is not important. 
				if (( oorcode&3)!=0){ hi.inc(); continue; }
			} 
			const D vx=composite_state(idj,ivx); // m_state(idj,ivx); // this node it the src node
			const D vy=composite_state(idj,ivy); // m_state(idj,ivy); 
			const D ni=composite_state(idj, in);
			// this needs to calculate a flow out into the y bulk 
			// but it is later not updated. 
			//kernel.calculate( idj,hi.r2(), -hi.delx(), -hi.dely(), ni,vx,vy,1.0);
			kernel.calculate( idj,hi.r2(), -hi.delx(), -hi.dely(), ni);
			hi.inc();
		}
		kernel.compile();
		hi.start();
		while (hi.ok())
		{
			const int idjgeo=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
//			if (dnode<0) { hi.inc(); continue; } 
			const int idj=m_idx.lut(idjgeo,true);
			D ni=composite_state(idj, in);
			D gni=0;

	// TODO FIXME Danger Will Robinson
	// This is where the stupid normalization creates a problem. The bulk pixels
	// normalize wrong since all the flow is going one way. 

			const int rc=kernel.result(gni,idj);
			// it is up to the dest to reject out flow into bulk off mesh nods
			if (rc==0) dest.gni(idj, idi, 1,gni/n0,ni,n0); 
			// want hi.relx() and rely() 
			if (debug_xfer) { MM_MSG( " moving "<<MMPR4(hi.relx(),hi.rely(),idi,idj)<<MMPR4(gni,gni/n0,ni,n0)<<hi.to_string()<<MMPR(rc)) }
			if (debug_xfer2) { MM_MSG( " moving "<<MMPR4(hi.centerx(),hi.centery(),hi.meshx(),hi.meshy())<<MMPR4(hi.relx(),hi.rely(),idi,idj)<<MMPR4(gni,gni/n0,ni,n0)<<MMPR4(rc,in,hi.srcx(),hi.srcy())) }
			hi.inc();
	}
	// the boundary pixels do not deplete, that is the values are fixed not the
// total amoun of charge and with in infite time steps they flood the 
// system. This may be more accurate even in regions of good parameter values
if (do_fixed_bulk)
{ 
	// for the given node that has just dispersed its 
	//contents to self and others on and off mesh,
	// add back some amount from a static bulk
	// there are carriers on both ends, need to add them together
	// and pick the right carriers n or p 
const int idbc=-1; //  dummy value <0 to be flagged as bulk of off grid 
			D gni=0;
			// it is up to the kernel to figure out what to do 
			const int rc= kernel.bulk(gni);
//			const int rc=kernel.result(gni,idi);
// the dest can not be confused calling gni with mixed up values, 
			//const int rc=kernel.result(gni,idj);
//	if (rc==0) dest.gni(idj, idi, 1,gni/n0,ni,n0); 
	// this is going on the rhs so only the net (gni/no)*no needs
	// to mean anything 
	const D ni=9;
	if (rc==0) dest.gni(idi, idbc, 1,gni/n0,ni,n0); 

}


		++serial; 	
	}



// This does nothing at the border nodes- the grid points are set equal to bulk with diri 
// and there is no field so they have net zero xfer.a
// but with penetration there can bea  net flow. 
if ( do_bc) {
	// this also needs to loop over the peripheral nodes for outflow. 
	// for now just inflow at y limits to avoid issues at the junction.
	typedef mjm_rectangle_itor ReI;
	// 100 may help 
	const D bcf=m_flp.bcf(); // 100;
	ReI rei(-m_bv.m_layers_x,0,m_n_x,m_n_x+m_bv.m_layers_x, -m_bv.m_layers_y,0,m_n_y,m_n_y+m_bv.m_layers_y);
	const int base= -m_bv.m_layers_x-m_bv.m_layers_y*(m_n_x+2*m_bv.m_layers_x); 
	// TODO FIXME Danger Will Robinson
	// This is where the stupid normalization creates a problem. The bulk pixels
	// normalize wrong since all the flow is going one way. 
	while (rei.ok())
	{
		// this should be negtive
		const int meshx=rei.x();
		const int meshy=rei.y();
		// if meshx and meshy are in range, then idigeo is ok for the nonpadded dims 
		const int idigeo=m_idx.index(meshx,meshy);
		const int  idi=m_idx.lut(idigeo,true); // this does nothing for off grid 
		if (idi>=0) { MM_ERR(" node is i nmesh "<<MMPR4(idi,idigeo,rei.x(),rei.y())) }
		// at least get rid of the ones which cannot contribute worry about terminal
		// normalization later
		const int  oorcode=m_idx.oor(meshx,meshy);
		// the oor layers at the y ends are the bulk media
		if (( oorcode&12)==0){ rei.inc(); continue; }
		if (( oorcode&3)!=0){ rei.inc(); continue; }

		D n0=composite_state(idi, in);
		MM_ONCE(" have a off grid node "<<n0,)
		if (n0==0) { rei.inc(); continue; }	
		MM_ONCE(" have a off grid node that passes  "<<n0,)
		const bool above=(meshy>0);
		const int idrefgeo=m_idx.index(meshx,above?(m_n_y-1):0);
		//const int idrefgeo=m_idx.index(meshx,above?(m_n_y-2):1);
		const int  idref=m_idx.lut(idrefgeo,true); // this does nothing for off grid 
		D vxref=composite_state(idref,ivx); // m_state(idj,ivx); //  src node
		D vyref=composite_state(idref,ivy); // m_state(idj,ivx); //  src node
		const D fac=bcf/D(1<<((above)?(meshy-m_n_y+1):(-meshy)));
		//kernel.point(idi,n0,0,0);
		kernel.point(idi,n0,vxref*fac,vyref*fac,m_fixed(idi,m_x),m_fixed(idi,m_y));
		//kernel.point(idi,n0,vx,vy,m_fixed(idi,m_x),m_fixed(idi,m_y));
		// geo_index is zero based including excludinged area 
		// so the (0,0) hood_itor uses needs to be tanslated 
		// -base makes the input wrong- geo_idx returns mesh coordinates, which is right  
		// this kluge made it worse lol 
		//Hi hi=m_idx.hood_itor(3.0*kernel.radius,rei.geo_idx()-base,0,0,true);
		//Hi hi=m_idx.hood_itor(kernel.radius,rei.geo_idx()-base,0,0,true);
		Hi hi=m_idx.hood_itor(kernel.radius,rei.geo_idx()-base,kernel.tvx,kernel.tvy,true);
		//Hi hi=m_idx.hood_itor(kernel.radius,idigeo,kernel.tvx,kernel.tvy);
		hi.start();
		while (hi.ok())
		{ 	// this can diffuse into multiple locations, need to get norm right. 	
			const int idjgeo=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
			const int idj=m_idx.lut(idjgeo,true); // dnodeloc;
			// in theory now we have the normalization worked out and only need this
			// for dynamic not static nodes
			// note this eliminates the source node 
			if (idj<0) { hi.inc(); continue; } 
		//	if (idjgeo>0)
				// this can't work as it is not a const of the center point
				// and will always be doomed at high fields. 
				const D vx=composite_state(idj,ivx); // m_state(idj,ivx); //  src node
				const D vy=composite_state(idj,ivy); // m_state(idj,ivy); 
				const D ni=composite_state(idj, in);
				// this calculates a flow to anything - including the bulk x pixels 
				//kernel.calculate( idj,hi.r2(), -hi.delx(), -hi.dely(), ni,vx,vy,bcf);
				kernel.calculate( idj,hi.r2(), -hi.delx(), -hi.dely(), ni);
			hi.inc();
		}
		// the kernel needs to normalize and conserve
		// since the source node is bulk, this really does nothing. 
		kernel.compile();
		hi.start();
		while (hi.ok())
		{	// the itor is over geographical nodes 
			const int idjgeo=m_idx.index(hi.meshx(),hi.meshy()); // hi.node();
			const int idj=m_idx.lut(idjgeo,true); // dnodeloc;
			// the normalization is supposed to save this but ill not work here 
			// this should be up to the dest to sort out however, the kernel should see all of these 
			if (idj<0) { hi.inc(); continue; } 
			const D ni=composite_state(idj, in);
			D gni=0;
			const int rc=kernel.result(gni,idj); // , hi.r2(),  -hi.delx(), -hi.dely(),0);
	 		if (rc==0) dest.gni(idj, idi, 1,gni/n0,ni,n0); 
//			MM_MSG( " movingbc  "<<MMPR4(gni,gni/n0,idi,idj)<<hi.to_string())
			if (debug_xfer) { MM_MSG( " movingbc "<<MMPR4(hi.relx(),hi.rely(),idi,idj)<<MMPR4(gni,gni/n0,ni,n0)<<hi.to_string()) } 
			if (debug_xfer2) { MM_MSG( " moving "<<MMPR4(hi.centerx(),hi.centery(),hi.meshx(),hi.meshy())<<MMPR4(hi.relx(),hi.rely(),idi,idj)<<MMPR4(gni,gni/n0,ni,n0)<<MMPR4(rc,in,hi.srcx(),hi.srcy())<<" rei") }
			hi.inc();
		}
		rei.inc();
	}
} // bc
//    { MM_ERR(" done flmb ") }
} // flmb_flow


///////////////////////////////////////////////////////////////////////////////////


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
void gr_point_by_point(MySparse & A, MyBlock & b, const D & dt, 
 const IdxTy in, const IdxTy ip, const IdxTy points)
{
	const FlpTy & flp= m_flp;
	const bool debug_gr=true;
	const D fudge= flp.gr_xsection(); // 1e-3;
	const D tau= flp.gr_tau();
	const D phitau=fudge/tau;
	for (IdxTy i=0; i<points; ++i)
	{
		const IdxTy idxn=1+2*i;
		const IdxTy idxp=2*i;
		const D n=m_state(i,in);
		const D p=m_state(i,ip);
		// phitau*(np-ni^2) =dn
		b(idxn)+=phitau*m_ni2*dt;
		b(idxp)+=phitau*m_ni2*dt;
		A(idxn,idxp)-=phitau*n*dt;
		A(idxp,idxn)-=phitau*p*dt;
	
	//	D dn=gr_b2b_point((n),(p),m_ni2,phitau,dt);

	}

}
void gr_point_by_point_equ(MySparse & A, MyBlock & b, const D & dt, 
 const IdxTy in, const IdxTy ip, const IdxTy points)
{
	const FlpTy & flp= m_flp;
	const bool debug_gr=true;
	const D fudge= flp.gr_xsection(); // 1e-3;
	const D tau= flp.gr_tau();
	const D phitau=fudge/tau;
	for (IdxTy i=0; i<points; ++i)
	{
		const IdxTy idxn=1+2*i;
		const IdxTy idxp=2*i;
		const D n=m_state(i,in);
		const D p=m_state(i,ip);
		const D sum=p+n;
		const D disc=::sqrt(sum*sum+4.0*m_ni2);
		if ((n*p)>m_ni2)
		{
			const D dsum=.5*sum;
			const D ddisc=-.5*disc;
			const D d=.5*sum-.5*disc;
			b(idxn)-=ddisc;
			b(idxp)-=ddisc;
			A(idxn,idxn)+=.5;
			A(idxp,idxp)+=.5;
	//		b(idxn)+=d;
//			b(idxp)+=d;
			continue;
		}
			const D dsum=.5*sum;
			const D ddisc=+.5*disc;
		const D d=.5*sum+.5*disc;
			b(idxn)+=ddisc;
			b(idxp)+=ddisc;
			A(idxn,idxn)+=.5;
			A(idxp,idxp)+=.5;
//		b(idxn)+=d;
//		b(idxp)+=d;

		// phitau*(np-ni^2) =dn
	//	b(idxn)+=phitau*m_ni2*dt;
//		b(idxp)+=phitau*m_ni2*dt;
//		A(idxn,idxp)-=phitau*n*dt;
//		A(idxp,idxn)-=phitau*p*dt;
	
	//	D dn=gr_b2b_point((n),(p),m_ni2,phitau,dt);

	}

}




void gr_point_by_point_integrating(MyBlock & sol, const D & dt, const IdxTy in, const IdxTy ip, const IdxTy points
, const IdxTy qpx, const IdxTy qpy)
{
	const FlpTy & flp= m_flp;
	const bool debug_gr=m_tree.debug_gr(); // true;
	const D fudge= flp.gr_xsection(); // 1e-3;
	const D tau= flp.gr_tau();
	const D phitau=fudge/tau;
	const D nx=qpx;
	const D ny=qpy;
	const D fac=1.0/(nx*ny);
	//MyBlock diffn(points),diffp(points);
	MyBlock diffn(points); // ,diffp(points);
	id_dof_itor dof_itor(m_mesh);
	while (dof_itor.ok())
    {
		const Elem * elem=* dof_itor.el(); 
		const id_dof_itor::dof_vec_type & dofs= dof_itor.id_dofs();
		for (IdxTy i=0; i<nx; ++i)
		{const D fx=D(i)/D(nx-1);
		for(IdxTy j=0; j<ny; ++j)
{
		const D fy=D(j)/D(ny-1);
		const D n=dof_itor.interpolate(sol,in,fx,fy);
		const D p=dof_itor.interpolate(sol,ip,fx,fy);
		D dn=gr_b2b_point((n),(p),m_ni2,phitau,dt);
		dof_itor.allocate_log(diffn,dn,fx,fy);
		//dof_itor.allocate(diffp,dp,fx,fy);
		//diffn(i)+=dn; diffp(i)+=dn;
		}}

		dof_itor.inc();
	} // while
	MM_ONCE(" klugin underflow until interp works",)
	for (IdxTy i=0; i<points; ++i) 
	{ 
		sol(i,in)+=diffn(i)*fac; sol(i,ip)+=diffn(i)*fac; 
	//	if (sol(i,in)<0) sol(i,in)=0;
	//	if (sol(i,ip)<0) sol(i,ip)=0;
	
	}

} // gr_point_by_point_integrating

void gr_point_by_point(MyBlock & sol, const D & dt, const IdxTy in, const IdxTy ip, const IdxTy points)
{
	const FlpTy & flp= m_flp;
	const bool debug_gr=m_tree.debug_gr(); // true;
	const D fudge= flp.gr_xsection(); // 1e-3;
	const D tau= flp.gr_tau();
	const D phitau=fudge/tau;
	MyBlock diffn(points),diffp(points);
	for (IdxTy i=0; i<points; ++i)
	{
		const D n=m_state(i,in);
		const D p=m_state(i,ip);

		D dn=gr_b2b_point((n),(p),m_ni2,phitau,dt);
		if (debug_gr) MM_MSG(" debug_gr "<<MMPR4(i,n,p,(n*p/m_ni2))<<MMPR3(dn,phitau,dt)<<MMPR2((n+dn),(p+dn)))
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
// 2017-07-21 changed the comments below not sure why ??? 
if (tau>expmax) { return nminus-nmin; } 
//if (tau>expmax) { return nminus; } 
if (tau<expmin) { return nplus-nmin; } 
//if (tau<expmin) { return nplus; } 

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
impose(m_n,m_bv.m_bc_n);
impose(m_p,m_bv.m_bc_p);


}
template <class Ty> void impose( const IdxTy col, Ty & m)
{
for(auto i=m.begin(); i!=m.end(); ++i)
 m_state((*i).first,col)=(*i).second;
}


void InitBC()
{
	m_bv.m_bc_n.clear();
	m_bv.m_bc_p.clear();

for ( IdxTy i=0; i<m_n_x; ++i)
{
const int idx=m_idx.lut(i,true);
const int idxe=m_idx.lut(m_nodes-i-1,true);
 m_bv.m_bc_n[idx]=m_bv.m_n_p;
 m_bv.m_bc_p[idx]=m_bv.m_p_p;
 m_bv.m_bc_n[idxe]=m_bv.m_n_n;
 m_bv.m_bc_p[idxe]=m_bv.m_p_n;


}
InitDiri();
}

void InitDiri()
{

	m_bv.m_diri.clear(); // should be in init	
const D vd=m_flp.vd();
for ( IdxTy i=0; i<m_n_x; ++i) m_bv.m_diri[m_idx.lut(i,true)]=0;
for ( IdxTy i=0; i<m_n_x; ++i) m_bv.m_diri[m_idx.lut(m_nodes-i-1,true) ]=vd;

}

// setup initial conditions for testing not reasl sim 
// pairs of impulses to check compare propogation
void initial_image(   MyBlock & sol, const bool z=true) {
if (z) initial_zero(sol);
initial_image(sol,1,m_u);
initial_image(sol,1,m_n);
initial_image(sol,1,m_p);

}
/* 
This is supposed to be a delta line charge or plate in 3D.
The rhs of laplac has the charge multipled by sf integral which
could be eliminated to make a special case but may just be
easier to scale here. 
*/
void initial_plates(const D & x1, const D & x2, const D & q, const int dir)
{
	MyBlock & sol=m_state;
	initial_zero(sol);
	zero_dopants();
	InitDiri(); // this is needed for the potential solver 
	const FlpTy& flp= m_flp; // 
	const IdxTy mesh_nx=flp.mesh_nx(); // 
	const IdxTy mesh_ny=flp.mesh_ny(); // 
	const D mesh_xmin=flp.mesh_xmin(); //
	const D mesh_ymin=flp.mesh_ymin(); // 
	const D mesh_xmax=flp.mesh_xmax(); // 
	const D mesh_ymax=flp.mesh_ymax(); // 
	const IdxTy in =m_na;
	const D dy=D(mesh_ny-1)/(mesh_ymax-mesh_ymin);
	const D x1nd=(x1-mesh_ymin)*dy; // /(mesh_ymax-mesh_ymin)*mesh_ny;
	const D x2nd=(x2-mesh_ymin)*dy; // /(mesh_ymax-mesh_ymin)*mesh_ny;
	int x1n=int(x1nd);
	if (x1n<0) x1n=0; 
	if (x1n>=mesh_ny ) x1n=mesh_ny-1; 
	int x2n=int(x2nd);
	if (x2n<0) x2n=0; 
	if (x2n>=mesh_ny ) x2n=mesh_ny-1; 
	int x1nf=x1n+1;
	int x2nf=x2n+1;
	if (x1nf>=mesh_ny ) x1nf=mesh_ny-1; 
	if (x2nf>=mesh_ny ) x2nf=mesh_ny-1; 

	const D f1=x1nd-x1n;
	const D f2=x2nd-x2n;
	// this is where unsigned saves if you can risk wrsap
	// it si never less than zed but will wrap oor
	MM_MSG(" plat planes "<<MMPR2(x1,x2))
	const bool dirx=false; // make field in y dir 
	const IdxTy nn=(dirx)?mesh_ny:mesh_nx;
	// the FEM solver integrates over the SF, this is a line charge or plate
	const D qplate= q*dy;
	for (IdxTy point=0; point<mesh_nx; ++point)
	{
		const IdxTy p1=m_idx.index(point,x1n);
		// this iss afe but use is not  
		const IdxTy p1f=m_idx.index(point,x1nf);
		const IdxTy p2=m_idx.index(point,x2n);
		const IdxTy p2f=m_idx.index(point,x2nf);

		{
			sol(p1,in)=-qplate*(1.0-f1); 
			sol(p1f,in)+=-qplate*(f1); 
			sol(p2,in)=qplate*(1.0-f2); 
			sol(p2f,in)+=qplate*(f2); 

		} 

	}

}

void initial_impulse(   MyBlock & sol, const bool z=true) {
if (z) initial_zero(sol);
initial_impulse(sol,1,m_u);
initial_impulse(sol,1,m_n);
initial_impulse(sol,1,m_p);


}

// place one charge on and on off the grid 
void initial_image(   MyBlock & sol,const D & v, const IdxTy col) {
if (m_n_x<4)
{
MM_ONCE(" grid sizes too small for default image test ",)
}
if (m_bv.m_layers_y<1)
{
MM_ONCE(" y layers  sizes too small for default image test ",)

}
int idgeo=m_idx.index( m_n_x-2, 0);
int id=m_idx.lut(idgeo,true);
composite_state2(id,col)=v;
idgeo=m_idx.index( 2, -1);
id=m_idx.lut(idgeo,true);
composite_state2(id,col)=-v;

}

void initial_impulse(   MyBlock & sol,const D & v, const IdxTy col) {
const int idgeo=m_idx.index( m_n_x>>1, m_n_y>>1);
const int id=m_idx.lut(idgeo,true);
sol(id,col)=v;
}

void zero_dopants()
{
const IdxTy pmax=m_nodes; // m_points-0;
for (IdxTy point=0; point<pmax; ++point)
{
m_fixed(point,m_na)=0;
m_fixed(point,m_nd)=0;
m_state(point,m_trap)=0;
}

}

void initial_zero(   MyBlock & sol) {

const IdxTy iu=m_u;
const IdxTy ip=m_p;
const IdxTy in=m_n;
const IdxTy iz=m_z;
const IdxTy pmax=m_nodes; // m_points-0;
for (IdxTy point=0; point<pmax; ++point)
{
	sol(point,iu)=0;
	sol(point,iz)=0;
	sol(point,ip)=0;
	sol(point,in)=0;

} // for
for (IdxTy point=0; point<(-m_idx.m_index_base); ++point)
{
m_idx.m_peri_state(point,iu)=0; 
m_idx.m_peri_state(point,iz)=0; 
m_idx.m_peri_state(point,ip)=0; 
m_idx.m_peri_state(point,in)=0; 
}

	m_bv.m_bc_n.clear();
	m_bv.m_bc_p.clear();
	m_bv.m_diri.clear(); // added later 
} // initial_zero

///////////////////////////////////////
void initial_guess(   MyBlock & sol) {
//MM_ERR(" right initial_guess callled doh ")
//sol.resize(size());
sol.resize(m_nodes,m_vars);
/*
// available if needed 
typedef mjm_interpolants Miu;
//typedef mjm_exp_interp Mi;
*/
const IdxTy iu=m_u;
const IdxTy ip=m_p;
const IdxTy in=m_n;
const IdxTy iz=m_z;
const IdxTy pmax=m_nodes; // m_points-0;
for (IdxTy point=0; point<pmax; ++point)
{
//	const D fac=1.0*point/(pmax-1);
		const bool lr= (m_fixed(point,m_y)>m_y_junction) ;
//{  m_fixed(i,m_nd)=Rev(Nd); m_fixed(i,m_na)=0;	} 
	if (point<m_mid) { sol(point,iu)= m_bv.m_v0; } else { sol(point,iu)= m_bv.m_vd; }
	if (!lr) { sol(point,in)=m_bv.m_n_p; sol(point,ip)=m_bv.m_p_p; } 
	else { sol(point,in)= m_bv.m_n_n; sol(point,ip)=m_bv.m_p_n; }
	sol(point,iz)=1;
}
	if (m_tree.enable_trapping())
	m_defects.net_charge(m_state,m_trap,m_nodes);

for (IdxTy point=0; point<(-m_idx.m_index_base); ++point)
{
if (m_idx.m_peri_fixed(point,m_y)<=m_y_junction)
{
m_idx.m_peri_state(point,in)=m_bv.m_n_p; 
m_idx.m_peri_state(point,ip)=m_bv.m_p_p; 
}
else
{
m_idx.m_peri_state(point,in)=m_bv.m_n_n; 
m_idx.m_peri_state(point,ip)=m_bv.m_p_n; 
}

}
	//m_idx.make_peripheral_nodes( nx, ny, m_layers, m_n_state, m_n_fixed);


}

// this may not work, the output never shows up 
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
	// these orders must bu preserved 
	m_p=pos;  ++pos; // 3
	m_n=pos; ++pos; // 4
	// these should be consecutive x,y
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
//    m_size=m_points*m_vars;
    m_iter=0;
	m_exit=false;
} // init

D neutral(const D & N, const D & ni2, bool plus)
{
const D a=.5*N;
//const D b=.5*::sqrt(N*N+4.0*m_ni2);
const D b=.5*::sqrt(N*N+4.0*ni2);
if (plus) return a+b ;

return a-b;
}


void InitConsts()
{

	//m_fixed.resize(m_points,3);
	//m_h=m_flp.total_length()/D(m_points-1); // 1e-6; // good value
	// this should be in the config file 
	m_y_junction=.7*(m_flp.mesh_ymin()+m_flp.mesh_ymax()); // 1e-6; // good value
	m_nodes=m_flp.mesh_nx()*m_flp.mesh_ny();
//	m_nodes=m_points;
    //m_size=m_points*m_vars;
	//m_Vt=m_flp.Vt(); // .0259;
	// silicon assumed
	//m_qe=m_flp.qe()/11.68; // 1.6e-19/8.854e-14;
	m_ni2=m_flp.ni2(); // 1.5e10*1.5e10;
	m_mup=m_flp.mup(); m_mun=m_flp.mun(); m_muz=m_flp.muz();
	const IdxTy midd=0+(m_nodes>>1);
	// this needs to be eliminated 
	m_mid=midd;
	m_bv.m_v0=m_flp.v0(); // 0;
	m_bv.m_vd=m_flp.vd(); // .5-.1*.5;
	//m_vd=.7;
}

void InitFixeds()
{
	const D Nd=m_flp.Nd(); // f*1e15; // + 
	const D Na=m_flp.Na(); // f*1e14; // - 
	// doh, these need to differ from Nx for neutrality. 
 	m_bv.m_n_n= Fwd(neutral(Nd,m_ni2,true )); // Fwd(1e15*f);
 	m_bv.m_p_n=Fwd(m_ni2/Rev(m_bv.m_n_n));
	m_bv.m_p_p=Fwd(neutral(Na,m_ni2, true)); // Fwd(1e14*f);
	m_bv.m_n_p=Fwd(m_ni2/Rev(m_bv.m_p_p)); //1e14;

	//const D midx=
	IdxTy ilast=0;
	for ( IdxTy i=0; i<m_nodes; ++i)
	{
		if (m_fixed(i,m_y)>m_y_junction) {  m_fixed(i,m_nd)=Rev(Nd); m_fixed(i,m_na)=0;	} 
		else {ilast =i; m_mid=i+1;  m_fixed(i,m_na)=Rev(Na); m_fixed(i,m_nd)=0;	} 
	}
	if (m_tree.enable_trapping())
	{
		m_defects.set(1,m_nodes);

		typedef Defects::defect_level Level;
		Level level(m_flp.trap_speed());
		//m_defects.add(1e13);
		//m_defects.add(1e13,level);
		m_defects.add(m_flp.trap_N(),level);
m_defects.set_to_equ( m_bv.m_n_n, m_bv.m_p_n, m_ni2, 0,  Nd,m_nodes-1);
m_defects.set_to_equ( m_bv.m_n_p, m_bv.m_p_p, m_ni2, Na,  0,0);
//MM_MSG("traps bc "<< MMPR(m_n_n)<<MMPR( m_p_n)<<MMPR(m_p_p)<<MMPR(m_n_p)<<MMPR(m_x_junction))
//void set_to_equ(const IdxTy first, const IdxTy end, const D & ni2, const D & Na, const D & Nd)
m_defects.set_to_equ(0, ilast+1, m_ni2, Na, 0);
m_defects.set_to_equ(ilast+1, m_nodes, m_ni2, 0,  Nd);
//MM_MSG(MMPR(ilast));
// does nothing as it is reset when initial_guess is done 
// well it does fick everything up since m_state size is not right yet 
//	m_defects.net_charge(m_state,m_trap,m_points);

	}
}
void ProjectOldCrap( const MyBlock & osol, const MyBlock & ofixed, const Defects & odefects, Myt& old)
{
const bool debug_project=false;
	const bool do_traps=m_tree.enable_trapping();
// this crap somes up alot lol 
const FlpTy& flp= old.m_flp; // right now this is being SET in Init not being used to set m_h
const IdxTy mesh_nx=flp.mesh_nx(); // 
const IdxTy mesh_ny=flp.mesh_ny(); // 
const D mesh_xmin=flp.mesh_xmin(); //
const D mesh_ymin=flp.mesh_ymin(); // 
const D mesh_xmax=flp.mesh_xmax(); // 
const D mesh_ymax=flp.mesh_ymax(); // 
const D dx=old.m_dx;
const D dy=old.m_dy;

	const NodeItor ne=m_mesh->nodes_end();
	for (NodeItor ni=m_mesh->nodes_begin(); ni!=ne; ++ni)
	{
		const Node & n=*(*ni);
		const IdxTy node_n=n.id();
		if (node_n>=m_nodes) { MM_MSG(" bad node "<<MMPR2(node_n,m_nodes))} 
		// location info should already be updated, we really do not need to use node itor
		const D x=m_fixed(node_n,m_x);
		const D y=m_fixed(node_n,m_y);
		// one less, in the case of a denser grid this should work ok.
		// if the new grid is coarser, this still works but all the detail
		// is lost  
		const IdxTy xold=(x-mesh_xmin)/dx;
		const IdxTy xold1=(xold>=(mesh_nx-1))?(mesh_nx-1):(xold+1);
		const IdxTy yold=(y-mesh_ymin)/dy;
		const IdxTy yold1=(yold>=(mesh_ny-1))?(mesh_ny-1):(yold+1);
		const IdxTy i1=old.m_idx.lut(old.m_idx.index(xold,yold),true);
		const IdxTy i2=old.m_idx.lut(old.m_idx.index(xold1,yold),true);
		const IdxTy i3=old.m_idx.lut(old.m_idx.index(xold1,yold1),true);
		const IdxTy i4=old.m_idx.lut(old.m_idx.index(xold,yold1),true);
		const D f1=(dx-(x-xold*dx))/(dx)*(dy-(y-yold*dy))/dy;
if (debug_project)
{		MM_MSG(MMPR4(dx,x,xold,f1)<<MMPR3(dy,y,yold)) }
		const D f2=((x-xold*dx))/(dx)*(dy-(y-yold*dy))/dy;
		const D f3=(x-xold*dx)/(dx)*(y-yold*dy)/dy;
		const D f4=(dx-(x-xold*dx))/(dx)*((y-yold*dy))/dy;
		m_state(node_n,m_u)= f1*osol(i1,m_u)+f2*osol(i2,m_u)+f3*osol(i3,m_u)+f4*osol(i4,m_u);
		m_state(node_n,m_n)= f1*osol(i1,m_n)+f2*osol(i2,m_n)+f3*osol(i3,m_n)+f4*osol(i4,m_n);
		m_state(node_n,m_p)= f1*osol(i1,m_p)+f2*osol(i2,m_p)+f3*osol(i3,m_p)+f4*osol(i4,m_p);
		if (do_traps)
		{

			m_defects.project_from( odefects, node_n, i1,  i2,i3,i4, f1, f2,f3,f4,false);

		}
if (debug_project)
	{	MM_MSG(MMPR4(f1,f2,f3,f4)<<MMPR4(m_state(node_n,m_n),osol(i1,m_n),osol(i2,m_n),osol(i3,m_n)))}
	} // node 
	// no need to fix bc's at this point hopefully that is enforced in the dirichlet maps.
	// not really needed here  except maybe debg 
if (do_traps)	m_defects.net_charge(m_state,m_trap,m_nodes);
// note that traps not done here. 
}

void InitSol()
{
	if (!m_save_sol) initial_guess(m_state);
	else
	{
	if (m_old==0) { 	MM_ERR(" m_old needs to be non null" ) } 
		Myt & old = *m_old;
		ProjectOldCrap(old.m_state,old.m_fixed, old.m_defects, old);
		if (m_tree.enable_trapping())
			m_defects.net_charge(m_state,m_trap,m_nodes);
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
MM_MSG( MMPR(m_bv.m_n_n)<<MMPR( m_bv.m_p_n)<<MMPR(m_bv.m_p_p)<<MMPR(m_bv.m_n_p)<<MMPR(m_y_junction))
// points monotonic, no nans or negatives etc. 
// perhaps even hmax/hmin diagnose and critique

}
//#endif

void InitNPEs()
{

if (m_del_np) delete m_es_np;
const StrTy name="np";
MM_ONCE(" will bomb on delete as es=es_np",)
m_es_np= (!m_del_np)?m_es:( new EsTy(*m_mesh));
(*m_es_np).add_system<NonLinTy> (name);
NonLinTy & system=  (*m_es_np).get_system<NonLinTy>(name);
m_jr_np=&system; 
(*m_jr_np).project_solution_on_reinit()=!true;
system.add_variable("p", libMesh::FIRST);
system.add_variable("n", libMesh::FIRST);
MM_MSG(" np variables added ")
if (m_del_np) (*m_es_np).init();
m_jr_np->get_dof_map().reinit(*m_mesh);
(*m_mesh).prepare_for_use(true);
(*m_es_np).reinit();
(*m_es_np).print_info();



}

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
(*m_es).print_info();
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
	// this needs to be changed to just one var
	m_nodes=m_mesh->n_nodes();
	// m_points=m_nodes;
	m_n_x=mesh_nx;
	m_n_y=mesh_ny;
	m_dx=(mesh_xmax-mesh_xmin)/(m_n_x-1);
	m_dy=(mesh_ymax-mesh_ymin)/(m_n_y-1);
 m_mesh->allow_renumbering(false);
 m_mesh->skip_partitioning(true);
  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  m_mesh->print_info();

	//m_diri.clear(); // should be in init	
}
void InitFromMesh()
{
	//cached valued cleared
	m_pot_mat.clear();
	// later????
	(*m_mesh).prepare_for_use(true);

 	m_n_state=10;
	m_n_fixed=4;
//	m_layers=m_flp.bulk_layers();
	m_bv.m_layers_x=m_flp.bulk_layers_x();
	m_bv.m_layers_y=m_flp.bulk_layers_y();
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
	//m_idx.make_peripheral_nodes( m_n_x, m_n_y, m_layers, m_n_state, m_n_fixed);
	m_idx.make_peripheral_nodes( m_n_x, m_n_y, m_bv.m_layers_x,m_bv.m_layers_y, m_n_state, m_n_fixed);
	// now need to locate and popylate the fakes layers 
	m_idx.set_positions(m_x,m_y,mesh_xmin,mesh_ymin,m_dx,m_dy);

}
void InitConvergence()
{
	delete m_ct;
	m_ct=0;
	const bool enable_convergence_track=m_tree.enable_convergence_track();
	if (!enable_convergence_track ) return;
	const D tol=m_flp.converge_tol();
	const StrTy lbl=m_flp.convergence_label();
	const StrTy ttype=m_flp.convergence_tracker();
	const IdxTy diter=m_flp.converge_iters();
	
	m_ct = new track_voltage( tol,diter,lbl);
	(*m_ct).set_locations(m_u, m_n,m_p,m_x, m_y);

}

//void Init(const IdxTy points=0, const MyBlock & x= MyBlock(1) ) {
void Init() {
// save whatever we have let use do this 
// although save and project flags could differ
	//if (m_save_sol) { InitPushOld(); }

	InitBase();
	InitConsts();
	InitMesh();
	InitFromMesh();
	InitFixeds();
	InitBC();
	InitEs();
	// this is probably superfluous as the first one does everything 
	if (m_tree.enable_ss()) InitNPEs();
	InitSol();
	SetCols(); // guess resizes trashing cols possible. 
	InitConvergence();
	InitVerify();
#if 0
	if (m_save_sol)
	{
	// user must do this 
	//	InitPushOld();
	}

#endif

}
// note that this does not always get the old config
// save these for "old"
// this is a lot easier than writing a copy ctor lol 
void InitPointers()
{
 m_mesh=0;
m_es=0; // : public EquationSystems
m_es_np=0; // : public EquationSystems
m_jr=0; // : public NonlinearImplicitSystem,
m_jr_np=0; // : public NonlinearImplicitSystem,
m_ct=0;
}

void InitPushOld()
{
// note that this does NOT get all the old config at
// the user has already reset that 
const bool stacking=false;
if (stacking)
{
m_old = new Myt(false);
*m_old=*this;
InitPointers();
// m_old->m_old is still valid and the stack is a linked list
// useful for retarded time later too. Deleting one of them
// should delete all the ancestors as m_old is deleted in dtor 
//m_old->m_old=0; // lol 
return; 
}

delete m_old;
m_old = new Myt(false);
*m_old=*this;
InitPointers();
// the old m_old was already deleted unless we really are stacking
m_old->m_old=0; // lol 

}


D Fwd(const D & x) const { return (x); } 
D Rev(const D & x) const { return (x); } 

Myt * m_old;
BraTy m_tree;
FlpTy m_flp; // right now this is being SET in Init not being used to set m_h
MyBlock m_state;
MyBlock m_fixed;
IdxTy m_iter;
IdxTy m_vars; 
CounterMap m_cm;
// misc stuff
//D m_Vt
//D m_qe;
// these are cached from flp and not reeally needed and make consistenty problems 
D m_ni2,m_mup, m_muz, m_mun;
boundary_values m_bv;
//D m_n_n, m_p_n,m_p_p,m_n_p,m_v0,m_vd
D m_y_junction;
IdxTy m_mid;
int m_u,m_z,m_p,m_n,m_v_nx,m_v_ny,m_v_px,m_v_py,m_q,m_trap;
int m_na,m_nd,m_x,m_y;
int m_n_state,m_n_fixed; //,m_layers_x,m_layers_y;
bool m_save_sol;
Defects m_defects;
bool m_exit;
InitTy * m_init;
MyMesh * m_mesh;
IdxTy m_nodes;
//DiriMap m_diri; // Dirichlet map from dof to value. 
//BCMap m_bc_n,m_bc_p;
EsTy *m_es;  EsTy *m_es_np; 
NonLinTy *m_jr;  NonLinTy *m_jr_np;  
bool m_del_np;
MeIdx m_idx;
D m_dx, m_dy;
IdxTy m_n_x,m_n_y;
MySparse m_pot_mat;
ConvTrack * m_ct;
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

