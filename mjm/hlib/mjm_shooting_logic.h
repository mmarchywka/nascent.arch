#ifndef MJM_SHOOTING_LOGIC_H__
#define MJM_SHOOTING_LOGIC_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
//#include "mjm_interpolation.h"
#include "mjm_instruments.h"
//#include "mjm_diffuse_iterators.h"

//#include <algorithm>
#include <string>
//#include <streams>
#include <sstream>

// it uses the fudding defaults fudd 
#include "mjm_petsc_util.h"
#include "mjm_petsc_fd_base.h"


/*
mjm_shooting8.h decision and param classes with unified interface 
*/

/*

while (i<argc)
{
const int istart=i;
config(m_tree,"-tree",i,argc,args);
config(m_flp,"-params",i,argc,args);
configi(m_points,"-points",i,argc,args);
config_set(  i,  argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
Init();
}
////////////////////////////////////////////////////////
// command block
void config_set( int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-set-param")
{
++i;
if (argc<=i) return; 
const StrTy nm=StrTy(args[i]);
m_flp.set(nm);
MM_ERR(" setting "<<nm<<" based on "<<s)
++i;
}

if (s=="-set-branch")
{
++i;
if (argc<=i) return; 
const StrTy nm=StrTy(args[i]);
m_tree.set(nm);
MM_ERR(" setting "<<nm<<" based on "<<s)
++i;
}



}
template <class Tx> void config(Tx & dest, const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
const StrTy nm=StrTy(args[i]);
MM_ERR(" loading file "<<nm<<" based on "<<s)
dest.load(nm);
++i; // consume param and cmd
}

}

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
*/

class ShootingValues
{
public:
class Tr
{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef std::string StrTy;
typedef std::stringstream Ss;

};

typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef ReadWriteMap RWMap;
ShootingValues( const StrTy & nm) : m_map(nm) {}
ShootingValues() {} 

void load(const StrTy & nm) { m_map.load(nm); } 
void set(const StrTy & nm) { m_map.set_line(nm); } 
StrTy  get(const StrTy & nm) { return m_map.get_string(nm," "); } 

//template <class Tx> void config(Tx & dest, const StrTy & cmd, int  & i, int argc, char ** args)
 void config( const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
const StrTy nm=StrTy(args[i]);
MM_ERR(" loading file "<<nm<<" based on "<<s)
load(nm);
++i; // consume param and cmd
}

}
//void config_set( int  & i, int argc, char ** args)
template <class Tcmd> void config_set(const Tcmd & cmd,  int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
//if (s=="-set-param")
if (s==StrTy(cmd))
{
++i;
if (argc<=i) return;
const StrTy nm=StrTy(args[i]);
set(nm);
MM_ERR(" setting "<<nm<<" based on "<<s)
++i;
}

}


protected:
RWMap m_map;

}; // ShootingValues


// all the hokey testing logic changes here easier to recompile than
// make a config file 
class BranchesVar :public ShootingValues
{
typedef ShootingValues Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
//BranchesVar( const StrTy & nm) : m_map(nm) {}
BranchesVar( const StrTy & nm) : Super(nm) {}
BranchesVar() {} 
//void load(const StrTy & nm) { m_map.load(nm); } 
//void set(const StrTy & nm) { m_map.set_line(nm); } 
 // DO NOT EDIT, generate from Branches using ./test_shoot.tex  -bramapgen
// start_brasigvar for text based extration 


bool show_iter_status() const { return m_map.get_bool("show_iter_status",true); }
bool always_dump_state() const { return m_map.get_bool("always_dump_state",true); }
bool dump_state_everyfew() const { return m_map.get_bool("dump_state_everyfew",true); }
bool dump_state_atend() const { return m_map.get_bool("dump_state_atend",true); }
bool dump_counter_everyfew() const { return m_map.get_bool("dump_counter_everyfew",true); }
bool cerr_counter_everyfew() const { return m_map.get_bool("cerr_counter_everyfew",true); }
bool check_for_negs() const { return m_map.get_bool("check_for_negs",true); }
bool shoot_itor_voltage() const { return m_map.get_bool("shoot_itor_voltage",true); }
bool do_reverse_shoot() const { return m_map.get_bool("do_reverse_shoot",false); }
bool shoot_voltage() const { return m_map.get_bool("shoot_voltage",true); }
bool integrate_gr() const { return m_map.get_bool("integrate_gr",true); }
bool integrate_sf_gr() const { return m_map.get_bool("integrate_sf_gr",true); }
bool enable_gr() const { return m_map.get_bool("enable_gr",true); }
bool track_drift_delta() const { return m_map.get_bool("track_drift_delta",true); }
bool track_fl_delta() const { return m_map.get_bool("track_fl_delta",true); }
bool enable_fl_flow() const { return m_map.get_bool("enable_fl_flow",true); }
bool enable_subpixel_fl_flow() const { return m_map.get_bool("enable_subpixel_fl_flow",!true); }
bool enable_sf_fl_flow() const { return m_map.get_bool("enable_sf_fl_flow",true); }
bool enable_drifting() const { return m_map.get_bool("enable_drifting",true); }
// changed default for itor 
bool enable_sub_drift() const { return m_map.get_bool("enable_sub_drift",!true); }
bool enable_sub_drift_itor() const { return m_map.get_bool("enable_sub_drift_itor",!true); }
bool enable_sf_drift() const { return m_map.get_bool("enable_sf_drift",true); }
bool defer_dd_updates() const { return m_map.get_bool("defer_dd_updates",!true); }
bool output_after_drift() const { return m_map.get_bool("output_after_drift",!true); }
bool enable_trapping() const { return m_map.get_bool("enable_trapping",!true); }
bool enable_mb_flow() const { return m_map.get_bool("enable_mb_flow",true); }
bool enable_sf_shooting() const { return m_map.get_bool("enable_sf_shooting", true); } 

StrTy to_string(const StrTy & sep=" ") const
{
// geneated code, do not edit 
Ss ss;

// ./test_shoot.tex  -bradumpgen
ss<<"show_iter_status="<<show_iter_status()<<sep;
ss<<"always_dump_state="<<always_dump_state()<<sep;
ss<<"dump_state_everyfew="<<dump_state_everyfew()<<sep;
ss<<"dump_state_atend="<<dump_state_atend()<<sep;
ss<<"dump_counter_everyfew="<<dump_counter_everyfew()<<sep;
ss<<"cerr_counter_everyfew="<<cerr_counter_everyfew()<<sep;
ss<<"check_for_negs="<<check_for_negs()<<sep;
ss<<"shoot_itor_voltage="<<shoot_itor_voltage()<<sep;
ss<<"do_reverse_shoot="<<do_reverse_shoot()<<sep;
ss<<"shoot_voltage="<<shoot_voltage()<<sep;
ss<<"integrate_gr="<<integrate_gr()<<sep;
ss<<"integrate_sf_gr="<<integrate_sf_gr()<<sep;
ss<<"enable_gr="<<enable_gr()<<sep;
ss<<"track_drift_delta="<<track_drift_delta()<<sep;
ss<<"track_fl_delta="<<track_fl_delta()<<sep;
ss<<"enable_fl_flow="<<enable_fl_flow()<<sep;
ss<<"enable_drifting="<<enable_drifting()<<sep;
/*
ss<<"show_iter_status="<<show_iter_status()<<sep;
ss<<"always_dump_state="<<always_dump_state()<<sep;
ss<<"dump_state_everyfew="<<dump_state_everyfew()<<sep;
ss<<"dump_state_atend="<<dump_state_atend()<<sep;
ss<<"dump_counter_everyfew="<<dump_counter_everyfew()<<sep;
ss<<"cerr_counter_everyfew="<<cerr_counter_everyfew()<<sep;
ss<<"check_for_negs="<<check_for_negs()<<sep;
ss<<"shoot_itor_voltage="<<shoot_itor_voltage()<<sep;
ss<<"do_reverse_shoot="<<do_reverse_shoot()<<sep;
ss<<"shoot_voltage="<<shoot_voltage()<<sep;
ss<<"integrate_gr="<<integrate_gr()<<sep;
ss<<"enable_gr="<<enable_gr()<<sep;
ss<<"track_drift_delta="<<track_drift_delta()<<sep;
ss<<"track_fl_delta="<<track_fl_delta()<<sep;
ss<<"enable_fl_flow="<<enable_fl_flow()<<sep;
ss<<"enable_subpixel_fl_flow="<<enable_subpixel_fl_flow()<<sep;
ss<<"enable_sf_fl_flow="<<enable_sf_fl_flow()<<sep;
ss<<"enable_drifting="<<enable_drifting()<<sep;
*/
ss<<"enable_sub_drift="<<enable_sub_drift()<<sep;
ss<<"enable_sub_drift_itor="<<enable_sub_drift_itor()<<sep;
ss<<"enable_sf_drift="<<enable_sf_drift()<<sep;
ss<<"defer_dd_updates="<<defer_dd_updates()<<sep;
ss<<"output_after_drift="<<output_after_drift()<<sep;
ss<<"enable_trapping="<<enable_trapping()<<sep;
ss<<"enable_mb_flow="<<enable_mb_flow()<<sep;
ss<<"enable_sf_shooting="<<enable_sf_shooting()<<sep;
return ss.str();
}


// end_brasigvar for text based extration 

//ReadWriteMap m_map;




}; // BranchVar

class FickLancParamsVar : public ShootingValues
{
typedef ShootingValues Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
//FickLancParamsVar( const StrTy & nm) : m_map(nm) {}
FickLancParamsVar( const StrTy & nm) : Super(nm) {}
FickLancParamsVar() : Super() {}
//void load(const StrTy & nm) { m_map.load(nm); } 
//void set(const StrTy & nm) { m_map.set_line(nm); } 

//StrTy  get(const StrTy & nm) { return m_map.get_string(nm," "); } 

 // DO NOT EDIT, generate from Branches using ./test_shoot.tex  -bramapgen
// start_flparamvar for text based extration 
// ./test_shoot.tex  -paramgen
D vmax() const { return m_map.get_double("vmax",1e6); }
D norm() const { return m_map.get_double("norm",.5/vmax()); }
D fl_t() const { return m_map.get_double("fl_t",1e-13); }
// if this is too wide, the bC forcer makes current spike 
D shooting_vd_tol_max() const { return m_map.get_double("shooting_vd_tol_max",.00001*1e-3); }
D shooting_vd_tol_min() const { return m_map.get_double("shooting_vd_tol_min",-.00001*1e-3); }
// need to make this depend on spacing or something,. 
D shooting_vi_min() const { return m_map.get_double("shooting_vi_min",-.01); } // .01;
D shooting_vi_max() const { return m_map.get_double("shooting_vi_max",.01); } // .01;
IdxTy shooting_maxiter() const { return m_map.get_uint("shooting_maxiter",100); } //100;
D drift_solver_rtol() const { return m_map.get_double("drift_solver_rtol",PETSC_DEFAULT); }
D drift_solver_abstol() const { return m_map.get_double("drift_solver_abstol",PETSC_DEFAULT); }
D drift_solver_dtol() const { return m_map.get_double("drift_solver_dtol",1e40); }
IdxTy drift_solver_maxits() const { return m_map.get_uint("drift_solver_maxits",1000); } //100;
D vsat() const { return m_map.get_double("vsat",1e7); } // 1e7 for Si
D gr_xsection() const { return m_map.get_double("gr_xsection",1e-12); } // 1e-3;
D gr_tau() const { return m_map.get_double("gr_tau",1e-8); }
IdxTy max_iterations() const { return m_map.get_uint("max_iterations",10000); }
IdxTy output_freq() const { return m_map.get_uint("output_freq",100); } // // 100;
IdxTy output_offset() const { return m_map.get_uint("output_offset",output_freq()-1); } // freq-1;
//D nominal_h() const { return m_map.get_double("nominal_h",.25e-6); }
D total_length() const { return m_map.get_double("total_length",.25e-6*2000); }
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
StrTy output_label() const { return m_map.get_string("output_label","fick"); }

//marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ 
// ./test_shoot.tex  -paramdumpgen
StrTy to_string(const StrTy & sep=" ") const
{
// geneated code, do not edit 
Ss ss;
ss<<"vmax="<<vmax()<<sep;
ss<<"norm="<<norm()<<sep;
ss<<"fl_t="<<fl_t()<<sep;
ss<<"shooting_vd_tol_max="<<shooting_vd_tol_max()<<sep;
ss<<"shooting_vd_tol_min="<<shooting_vd_tol_min()<<sep;
ss<<"shooting_vi_min="<<shooting_vi_min()<<sep;
ss<<"shooting_vi_max="<<shooting_vi_max()<<sep;
ss<<"shooting_maxiter="<<shooting_maxiter()<<sep;
ss<<"drift_solver_rtol="<<drift_solver_rtol()<<sep;
ss<<"drift_solver_abstol="<<drift_solver_abstol()<<sep;
ss<<"drift_solver_dtol="<<drift_solver_dtol()<<sep;
ss<<"drift_solver_maxits="<<drift_solver_maxits()<<sep;
ss<<"vsat="<<vsat()<<sep;
ss<<"gr_xsection="<<gr_xsection()<<sep;
ss<<"gr_tau="<<gr_tau()<<sep;
ss<<"max_iterations="<<max_iterations()<<sep;
ss<<"output_freq="<<output_freq()<<sep;
ss<<"output_offset="<<output_offset()<<sep;
//ss<<"nominal_h="<<nominal_h()<<sep;
ss<<"total_length="<<total_length()<<sep;
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
ss<<"output_label="<<output_label()<<sep;

return ss.str();
} 

// end_flparamvar for text based extration 

//ReadWriteMap m_map;
};

#endif

