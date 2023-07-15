#ifndef MJM_CONVERGENCE_TRACKING_H__
#define MJM_CONVERGENCE_TRACKING_H__

//#include "mjm_fem_indexing.h"
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
//#pragma GCC diagnostic ignored "-Wunused-parameter"
//#include "mjm_libmesh.h"
// fick 
//#pragma GCC diagnostic ignored "-Wunused-parameter"





/*

2017-08-13


*/
//////////////////////////////////////////////////////////////////

class convergence_track_base
{

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<int> MyIntBlock;
typedef mjm_block_matrix<IdxTy> MyuIntBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 
public:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MyIntBlock  MyIntBlock;
typedef Tr::MyuIntBlock  MyuIntBlock;
typedef Tr::MySparse MySparse;
typedef int intt;
typedef D fom_type;
typedef D result_type;
typedef std::vector<fom_type> FomHistory;
typedef std::vector<result_type> ResHistory;
typedef mjm_rational RatTy;
convergence_track_base(): m_valid(false),m_failed(false),m_converged(false),m_abort(false) {Init(); }
// this needs a lot more crap for traps etc
virtual void failed(const bool fail) { m_failed=fail; } 
virtual void add_iteration(const MyBlock & state, const MyBlock & fixed, const IdxTy nodes, const CounterMap * cm)
{

}
virtual bool converged() const { return m_converged; }
virtual bool abort() const { return m_abort; }
virtual bool exit() const { return m_abort||m_converged; }
void set_locations(const IdxTy ui, const IdxTy ni, const IdxTy pi, const IdxTy xi, const IdxTy yi)
{ m_u=ui; m_n=ni; m_p=pi; m_x=xi; m_y=yi; }

protected:
void Converged(const bool x ) { m_converged=x; }
void Abort(const bool x ) { m_abort=x; }
bool Failed( ) const { return m_failed; }
void Init()
{
const IdxTy bad=~0U;
m_count=0;
m_u=bad;
m_n=bad;
m_p=bad;
m_x=bad;
m_y=bad;
m_valid=true;
}
void Result(const result_type & fom) { m_results.push_back(fom); }
StrTy Approach(const fom_type & final)
{
Ss ss; 
const IdxTy foms=m_results.size();
for (IdxTy i=0; i<foms; ++i)
{
const D  & f= m_results[foms-1-i];
D norm= (f-final)/final;
RatTy r= norm;
ss<<MMPR4(i,r,f,norm);
}
return ss.str();
}
bool m_valid,m_failed,m_converged,m_abort;
IdxTy m_u,m_n,m_p;
IdxTy m_x,m_y;
IdxTy m_count;
ResHistory m_results;

}; // convergence_track_base


class track_voltage : public convergence_track_base
{
typedef track_voltage Myt;
typedef convergence_track_base Super;
typedef Super Tr;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MyIntBlock  MyIntBlock;
typedef Tr::MyuIntBlock  MyuIntBlock;
typedef Tr::MySparse MySparse;
typedef int intt;


//convergence_track_base(): m_valid(false),m_converged(false),m_abort(false) {Init(); }
// this needs a lot more crap for traps etc
public:
track_voltage( const D & tol, const IdxTy diter,const StrTy & lbl )
: m_tol(tol),m_iters(diter),m_label(lbl) 
{
MM_MSG(" making track_voltage "<<MMPR3(lbl,tol,diter))
Init(); }
virtual void add_iteration(const MyBlock & state, const MyBlock & fixed, const IdxTy nodes, const CounterMap * cm)
{
if (Failed()) { Abort(true); } 
++m_iter;
if ((m_iter%m_iters)!=0) return; 
D umin=1e100;
D umax=-1e100;
for (IdxTy i=0; i<nodes; ++i)
{
const D & v=state(i,m_u);
if (v>umax) umax=v;
if (v<umin) umin=v;
}
const D d=umax-umin;
const D del=d-m_last;
const D avg=.5*(d+m_last);
MM_MSG(MMPR4(m_iter,m_last,del,avg)<<MMPR(d))
MM_ERR(MMPR4(m_iter,m_last,del,avg)<<MMPR(d))
if (avg==0) { Converged(del==0); return; }
const D fom=del/avg;
if (fom<m_tol) if (fom> -m_tol) {Converged(true); }
MM_MSG(MMPR3(m_iter,d,fom))
MM_ERR(MMPR3(m_iter,d,fom))
Result(d);
if (converged()) 
{
MM_MSG("CONVERTED "<<m_label<<" "<<MMPR3(m_iter,d,fom))
MM_ERR("CONVERTED "<<m_label<<" "<<MMPR3(m_iter,d,fom))
MM_MSG(" APPROACH "<<Approach(d))
MM_ERR(" APPROACH "<<Approach(d))
}
m_last=d;
}

protected:
void Init()
{
m_iter=0;
m_last=1e100;
}

D m_tol;
IdxTy m_iters;
StrTy m_label;
IdxTy m_iter;
D m_last;

}; // track_voltage


#endif

