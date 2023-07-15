#ifndef MJM_TIMeLINE_LOGIC_H__
#define MJM_TIMELINE_LOGIC_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_instruments.h"

#include <string>
#include <sstream>


/*
mjm_shooting8.h decision and param classes with unified interface 
*/

/*
*/

class TimelineValues
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
TimelineValues( const StrTy & nm) : m_map(nm) {}
TimelineValues() {} 
enum {BAD=~0};
void load(const StrTy & nm) { m_map.load(nm); } 
void set(const StrTy & nm) { m_map.set_line(nm); } 
StrTy  get(const StrTy & nm) { return m_map.get_string(nm," "); } 

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
class TimelineBranches :public TimelineValues
{
typedef TimelineValues Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
//BranchesVar( const StrTy & nm) : m_map(nm) {}
TimelineBranches( const StrTy & nm) : Super(nm) {}
TimelineBranches() {} 


//bool show_iter_status() const { return m_map.get_bool("show_iter_status",true); }

StrTy to_string(const StrTy & sep=" ") const
{
// geneated code, do not edit 
Ss ss;

// ./test_shoot.tex  -bradumpgen
//ss<<"show_iter_status="<<show_iter_status()<<sep;
return ss.str();
}


// end_brasigvar for text based extration 


}; // BranchVar

class TimelineParams : public TimelineValues
{
typedef TimelineValues Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
TimelineParams( const StrTy & nm) : Super(nm) {}
TimelineParams() : Super() {}

// THIS MUST BE GENETATED OR TYPOS WIL TAKE FOREVEROTO FIND 
 // DO NOT EDIT, generate from Branches using ./test_shoot.tex  -bramapgen
// start_flparamvar for text based extration 
// ./test_shoot.tex  -paramgen
//D vmax() const { return m_map.get_double("vmax",1e6); }
//IdxTy drift_solver_maxits() const { return m_map.get_uint("drift_solver_maxits",1000); } //100;
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }

IdxTy ordinal_column() const { return m_map.get_uint("ordinal_column",BAD); } //100;
IdxTy time_column() const { return m_map.get_uint("time_column",0); } //100;
//marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ 
// ./test_shoot.tex  -paramdumpgen
StrTy to_string(const StrTy & sep=" ") const
{
// geneated code, do not edit 
Ss ss;
//ss<<"vmax="<<vmax()<<sep;
ss<<"ordinal_column="<<ordinal_column()<<sep;
ss<<"time_column="<<time_column()<<sep;
return ss.str();
} 

// end_flparamvar for text based extration 

};

#endif

