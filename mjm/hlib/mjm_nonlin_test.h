#ifndef MJM_NONLIN_TEST_H__
#define MJM_NONLIN_TEST_H__
#include "mjm_globals.h"
#include "mjm_libmesh.h"
#include "mjm_libmesh_boundaries.h"
#include "mjm_assembler_nlog.h"
#include <mjm_block_matrix.h>
#include <mjm_sparse_matrix.h>
#include <mjm_instruments.h>

// LibMesh includes
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/serial_mesh.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <math.h>

// go ahead and include here for now
#include <mjm_electrode_geometry.h>
#include <mjm_system_params.h>
// https://sourceforge.net/p/libmesh/mailman/message/34518120/
#include "/home/marchywka/d/petsc/petsc-3.7.3/include/petscmat.h"
//#include "/home/marchywka/d/petsc/petsc-3.7.3/include/PETScMatrix.h"
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/petsc_matrix.h>

using libMesh::EquationSystems;
using libMesh::Gradient;
using libMesh::NonlinearImplicitSystem;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Parameters;
using libMesh::Point;
using libMesh::SparseMatrix;
using libMesh::System;
using libMesh::UnstructuredMesh;
using libMesh::NonlinearImplicitSystem;

//template <class Ty> 
class mjm_nonlin_es : public EquationSystems
{
typedef EquationSystems Super;
typedef mjm_nonlin_es Myt;


public:

  	static void Create(Myt ** b, const libMesh::Parallel::Communicator & comm)
{
MM_ERR(" Create mjm_nonlin_es")

}
  	static void Destroy(Myt ** b)
	{
MM_ERR(" Destroy  mjm_nonlin_es")

}

	//mjm_nonlin_es(UnstructuredMesh & m)
template <class FUDD> 	mjm_nonlin_es(FUDD & m)
: Super(m)
{
MM_ERR(" Ctor  mjm_nonlin_es")
//this fails because one system has to be added ibefore
//this->EquationSystems::init();
}
  ~mjm_nonlin_es()
  {
MM_ERR(" Dtor  mjm_nonlin_es")
    // delete _meshRefinement;
  }

  // Public interface functions
  void viewParameters()
{
MM_ERR(" view  mjm_nonlin_es")

}
// this has to be acalled after a system is added 
  void init(){
MM_ERR(" init  mjm_nonlin_es")
this->EquationSystems::init();
}
  void step(const Real & dt = -1.0){}
  void output(int timestep, const Real & t, Real & o_t, bool force = false)
{
MM_ERR(" output  mjm_nonlin_es")

}
  void run(){

MM_ERR(" run  mjm_nonlin_es")

}



}; // mjm_nonlin_es



template <class Ty> 
class mjm_numerical_jr : public NonlinearImplicitSystem,
                       public NonlinearImplicitSystem::ComputeResidualandJacobian,
                       public NonlinearImplicitSystem::ComputeBounds,
                       public System::Initialization
{
typedef mjm_numerical_jr Myt;
typedef Ty assembler_type;
typedef Ty MLATy;
typedef mjm_numb CkTy;
typedef unsigned int IdxTy;
typedef mjm_block_matrix<Real> MyBlock;
typedef mjm_sparse_matrix<Real> MySparse;
typedef mjm_block_matrix<IdxTy> MyIdxBlock;

typedef std::string StrTy;
//typedef EquationSystems EsTy;
typedef mjm_nonlin_es EsTy;
public:
typedef libMesh::dof_id_type Dof;
typedef std::map<Dof,Real> DiMap;
typedef DiMap::const_iterator DiMapItor;

typedef ReadWriteMap RWMap;
typedef electrode_geometry EGeo;
typedef system_params ConditionsTy;

typedef std::map<IdxTy,Real> ZeroEqn;
typedef mjm_libmesh_util::NEVec SurfVec;
typedef mjm_special_qrule MqTy;
typedef std::map<const Elem *, MqTy  > SpecialElements;
// slowly convert to my local sparses... 
	typedef mjm_sparse_matrix<Real> Tmjm_sparse;

mjm_numerical_jr(EquationSystems & eqSys,
                   const StrTy  & name_in,
                 const IdxTy number_in) :
	NonlinearImplicitSystem(eqSys, name_in, number_in),
	m_map("nonlinjr.config"), m_dirichlet(0),m_zero(0), m_surfaces(0)
	,m_es(dynamic_cast<EsTy &>(eqSys)), m_mla(0)
	,m_transuranics(0)
	, m_cm(&m_cm_default)
	, m_iteration(0)
	,m_ass_flags(0)
	,m_eqn_factors(MyBlock(1))
{ ctor(); }
mjm_numerical_jr(EquationSystems & eqSys,
                   const StrTy  & name_in,
                 const IdxTy number_in, CounterMap & cm) :
	NonlinearImplicitSystem(eqSys, name_in, number_in),
	m_map("nonlinjr.config"),m_dirichlet(0),m_zero(0),m_surfaces(0)
	, m_es(dynamic_cast<EsTy &>(eqSys)),
	m_mla(0)
	,m_transuranics(0)
	, m_cm(&cm)
	, m_iteration(0)
	,m_ass_flags(0)
	,m_eqn_factors(MyBlock(1))
{ ctor(); }

~mjm_numerical_jr() { MM_ERR(" fudding dtor caled fon fudding system")} 

void set_iteration(const IdxTy i) { m_iteration=i; } 
enum ASSEMBLY_FLAGS { IMPOSE_ELECTRODE_BC=1, SUPPRESS_TRANSURANIACS=2};
void set_assembly_flags(const IdxTy i) 
	{ m_ass_flags=i; MM_ERR(" setting assembly flag "<<m_ass_flags)   } 
bool impose_electrode_bc() const 
	{ return (0!=(m_ass_flags&IMPOSE_ELECTRODE_BC)); } 
bool suppress_transuraniacs() const 
	{ 
MM_ONCE(" suppress_transuraniacs hard coded to true", return true; )
//MM_ONCE(" suppress_transuraniacs hard coded to false", return !true; )

return (0!=(m_ass_flags&SUPPRESS_TRANSURANIACS)); } 

void mjm_linear_solve()
{
// system matrix, precondition, sol, rhs, tol, iterationr
// this just calls the base which is now obsolete or null 
//assemble();
if (true) {MM_ONCE(" avoding stupid  linear solve ",return; ) } 
MM_ONCE(" doing linear solve ",)
 NumericVector<Number> * null_r=0;
 SparseMatrix<Number> * null_j=0;

assembly(null_r,null_j,*solution);
// no pre
get_linear_solver()->solve(*matrix,0,*solution,*rhs,1e-18,5000);

}

//void FUKKK()
void DumpSolution()
{
// this has no fuddng effect, the only valid thing is a zero fudd 
//*solution=*current_local_solution;
//mjm_libmesh_util::DumpSolutionNodes(std::cout,mesh, *g_es,*g_system,nvar);
//mjm_libmesh_util::DumpSolutionNodes(std::cout,m_mla->mesh, m_es,*this,4);
mjm_libmesh_util::DumpVectorHiRes(std::cout, m_mla->mesh, (*solution), "SolutionDump");
mjm_libmesh_util::DumpActivePolygons(std::cout, m_mla->mesh, (*solution), "Poly");
}


void fudd(CounterMap & _cm) { m_cm=&_cm; } 
void ctor()
{
MM_ERR( " nl ctor ")
//assemble();

// it seems to compute but dump fails on dof fudd???
get_dof_map(); //???????
  // Give the system an object to compute the initial state.
  	attach_init_object(*this);
	project_solution_on_reinit()=true;
  	// Attache the R & J calculation object
  	nonlinear_solver->residual_and_jacobian_object = this;
  	// Attach the bounds calculation object
  	nonlinear_solver->bounds_object = this;
// this only fudding workse where it will get wiped out later FUDD 
//project_solution(Myt::init_value<MLATy>, libmesh_nullptr, m_es.parameters);
// this ASSFUDD compiles but the one that shouldwork does not fudd
// NonlinearImplicitSystem::init_data();
//m_mla->wtf(false);  
} // ctor

void another_fudd()
{
MM_ERR( " reproject fudd  ")
MM_ERR( " doing the init_value call  ")

project_solution(Myt::init_value<MLATy>, libmesh_nullptr, m_es.parameters);

}

// must be called after adding to systems
void super_init_fudd()
{
MM_ERR( " super_init_fudd ")
 NonlinearImplicitSystem::init_data();
}

// called when es init is called
void initialize()
{
MM_ERR( " nl init called ")
// this fudds up the fudding doff map. 
//ImplicitSystem::init_data();
//if (m_mla) get_dof_map().reinit(m_mla->mesh); //???????
get_dof_map(); //???????
if (m_mla)
{
// get_dof_map().reinit(m_mla->mesh); //???????
if ((&(m_mla->mesh)) !=(&(m_es.get_mesh()))) MM_ERR(" fudd mesh mismatch shot")


}
//get_dof_map().reinit(m_es.get_mesh()); //???????
//m_es.get_mesh().renumber_nodes_and_elements();
//get_dof_map().reinit(m_es.get_mesh()); //???????
// NonlinearImplicitSystem::initialize();
// this ass fudd has no fudding effect everythihng is zero fudd 
project_solution(Myt::init_value<MLATy>, libmesh_nullptr, m_es.parameters);
MM_ERR( " nl init data calling  ")

//ImplicitSystem::init_data();
//ImplicitSystem::init_data();
// NonlinearImplicitSystem::init_data();
//NonlinearImplicitSystem::init_data();
//ExplicitSystem::init_data();
// fudding priate 
//add_system_matrix();
//assemble();
// this should project a solution
/*
Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var )
{
    // we could just make a global instance 
    return MLATY::init_value<Number,Parameters>(p,parameters,system,var,g_cond);
}

void init_cd (EquationSystems & es,
              const std::string & system_name)
{
    MM_MSG(" init_cd")
    MLATY::init_cd(es,system_name);
    //g_system = new mjm_numerical_jr(es,system_name,0);

// hopefully initialize gives us better appraoch 
//      system.project_solution(init_value, libmesh_nullptr, es.parameters);





*/

} // initialize
  /**
   * The residual and Jacobian assembly function for the Biharmonic system.
   */
  void residual_and_jacobian(const NumericVector<Number> & u,
                             NumericVector<Number> * R,
                             SparseMatrix<Number> * J,
                             NonlinearImplicitSystem & nlis)
{
//PetscMatrix A=PetscMatrix(matrix);
libMesh::PetscMatrix<Real>*  A=(libMesh::PetscMatrix<Real>*)(matrix);
MatSetOption(A->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE)   ;
const bool have_real_j=true;
MM_ERR( " r and j called  "<<u.size())
if (!have_real_j) default_diags(u,R,J,nlis) ;
if(R)(*R).close();
assembly(R,J,u);
if(R)(*R).close();
//MM_ERR(" ASSFUDD SHOT FUDD")
MM_ONCE(" disable boundary enforce on J and R, should be done on my_sparse",)
//handle_boundaries(u,R,J,nlis) ;
if(R)(*R).close();

//MM_ERR(" ASSFUDD SAGAIN FUDD HIT FUDD")
if (R){ (*R).close(); } 

MM_ERR(" end of res ad j calc maybe print residuals  ")

//if (R) (*R).close();
//if (J) (*J).close();
    if (false) {
 {MM_ONCE(" printing reisudals lots of output  ",) } 
if(R) {
    const StrTy ival=(*m_cm).format_entry("iterations");
mjm_libmesh_util::DumpVector(std::cerr, m_mla->mesh, (*R), StrTy("res")+ival);
}
} else {MM_ONCE(" not printing reisudals here but could be in assembly  ",) } 
} // r and j

  void default_diags(const NumericVector<Number> & u,
                             NumericVector<Number> * R,
                             SparseMatrix<Number> * J,
                             NonlinearImplicitSystem & nlis)
{
MM_ERR( "default digas "<<u.size())
if (J) { 
for (Dof dof=0; dof<u.size(); ++dof)
	 (*J).set(dof,dof,1);
}

} // default_diags

// handle elements with boundary issues found during assembly 
template <class Tvec, class Telems>
  void handle_boundaries(const Telems & be, const Tvec & sol, Tvec  * R, Tmjm_sparse * J)
{
MM_ERR(" boundary size is "<<be.size())
IdxTy cnt=0;
for (auto ii=be.begin(); ii!=be.end(); ++ii)
{
//	MM_ERR(" boundary count "<<cnt)
	++cnt;
	const Elem * elem=*ii;
	const IdxTy bsurface=return_surface_code(elem);
	if (bsurface==0) { MM_ONCE(" have zero surface in boundary list ERROR",) }
	assemble_boundary_elem(elem,sol,R,J,bsurface);
/*	const IdxTy surf=(_bsurface>>16)-1;
	if (surf==0) continue;
	const IdxTy nodes=(*elem).n_nodes();
    for (IdxTy side=0; side<nodes; ++side)
    {
        const IdxTy bit=(1<<side);
        if ((bsurface&bit)==0) continue;
			
	} // side 
*/
} // ii 

}
//	typedef mjm_sparse_matrix<Real> Tmjm_sparse;
/// should implement dirichlet and arbitrary on given
// lists 
template <class Tvec>
  void handle_boundaries(const Tvec & sol, Tvec  * R, Tmjm_sparse * J)
{
	if (J==0) if (R==0) return; 
	if (m_dirichlet==0) return;
	std::vector<Dof> dofs;
	const Real fac=1e2;
	MM_ONCE(" enforcing NEW diri with huge fac "<<fac,)
	DiMapItor di=(*m_dirichlet).begin();
	DiMapItor df=(*m_dirichlet).end();
	while (di!=df)
	{
		const Dof dof=(*di).first;
		const Real v=(*di).second;
		if (J) dofs.push_back(dof);
		//if (R) (*R).set(dof,fac*(u(dof)-v));	
		if (R) (*R)(dof)=-fac*(sol(dof)-v);	
		++di;
	}
	// this needs to zero each row in dofs and put fac on diagonal
	if (J){  (*J).zero_rows(dofs,-1*fac); }

}

// implements dirichlet on given list
  void handle_boundaries(const NumericVector<Number> & u,
                             NumericVector<Number> * R,
                             SparseMatrix<Number> * J,
                             NonlinearImplicitSystem & nlis)
{
MM_ONCE(" old boundary handler called doh",)
if (J==0) if (R==0) return; 
if (m_zero!=0)
{
} // m_zero

//MM_ERR( "handle bounaries "<<u.size())
if (m_dirichlet==0) return;
std::vector<Dof> dofs;
const Real fac=1e2;
MM_ONCE(" enforcing diri with huge fac "<<fac,)
DiMapItor di=(*m_dirichlet).begin();
DiMapItor df=(*m_dirichlet).end();
while (di!=df)
{
	const Dof dof=(*di).first;
	const Real v=(*di).second;
//	if (J) (*J).set(dof,dof,1);
	if (J) dofs.push_back(dof);
	if (R) (*R).set(dof,fac*(u(dof)-v));	
	++di;
}
// petsc complains this needs to be asembled first wtf?
if (J){ (*J).close(); (*J).zero_rows(dofs,1*fac); }

} //handle_boundaries 


  /**
   * Function defining the bounds of the Biharmonic system.
   */
  void bounds(NumericVector<Number> & XL,
              NumericVector<Number> & XU,
              NonlinearImplicitSystem &)
{
MM_ERR( " vounds called  "<<XL.size()<<" "<<XU.size())



} //bound

//template <class Tt,class Tv, class Tx, class Ty,class Tc >
 //       void assemble_given(Tt & Je_var, Tv &  Ge_var, Tc & fe,
 //   ,const Tx & values,const Ty & grads, const IdxTy qp const bool do_lhs)

template <class Tt,class Tv, class Tx, class Tyz,class Tc,class Ts >
        void assemble_some(Tt & Je_var, Tv &  Ge_var
		, Tx & lhs_block, Tyz & rhs_block, Tc & fe, const Ts & sol, const IdxTy qp, const bool do_lhs,
	const bool is_side, const IdxTy side, const Elem * elem)
{
		//fes.fe->reinit(elem,kluge);
		if (is_side) fe.fe->reinit(elem,side);
	// this SHOULD have been done once for the element 
   //     else fe.fe->reinit(elem);
        const IdxTy qpsz=fe.qrule.n_points();
		if (qp>=qpsz ) return ;
		NumVec values; GradVec grads;
   		mjm_libmesh_util::MakeXandGradXVector( values, grads,(*m_mla).di.dof_indices_var, fe , qp, elem, (*m_mla).nvar,sol );
 		(*m_mla).assemble_given(Je_var,  Ge_var, fe , values, grads,  qp , do_lhs);
 		(*m_mla).assemble_given_block(lhs_block, rhs_block, fe , values, grads,  qp , do_lhs,sol);
		}


template < class Tx, class Tyz,class Tc,class Ts >
        void assemble_some(
		 Tx & lhs_block, Tyz & rhs_block, Tc & fe, const Ts & sol, const IdxTy qp, const bool do_lhs,
	const bool is_side, const IdxTy side, const Elem * elem)
{
		//fes.fe->reinit(elem,kluge);
		if (is_side) fe.fe->reinit(elem,side);
// should have been done ince 
  //      else fe.fe->reinit(elem);
        const IdxTy qpsz=fe.qrule.n_points();
		if (qp>=qpsz ) return ;
		NumVec values; GradVec grads;
		// now this has to do it for every qp 
   		mjm_libmesh_util::MakeXandGradXVector( values, grads,(*m_mla).di.dof_indices_var, fe , qp, elem, (*m_mla).nvar,sol );
 	//	(*m_mla).assemble_given(Je_var,  Ge_var, fe , values, grads,  qp , do_lhs);
	// I think that now assemble_given_block is the only one used 
 		(*m_mla).assemble_given_block(lhs_block, rhs_block, fe , values, grads,  qp , do_lhs,sol,elem);
		}




template <class Tt,class Tv, class Tx, class Tyz,class Tc,class Ts >
        void assemble_points(Tt & Je_var, Tv &  Ge_var
		, Tx & lhs_block, Tyz & rhs_block,  const Ts & sol, const Tc * elem, const bool do_lhs,
	const bool is_side, const IdxTy side)
{
        const IdxTy qpsz=(*m_mla).fev.qrule.n_points();
        for (IdxTy qp=0; qp<qpsz; qp++)
        {
         assemble_some(Je_var, Ge_var, lhs_block, rhs_block
			, (*m_mla).fev, sol,qp, do_lhs,false, 0,elem);
		} //qp
		for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
        	if (elem->neighbor(kluge)==libmesh_nullptr)
        	{   
        		const IdxTy qpssz=(*m_mla).fes.qrule.n_points();
        		for (IdxTy qp=0; qp<qpssz; qp++)
         			assemble_some(Je_var, Ge_var, lhs_block, rhs_block, 
						(*m_mla).fes, sol, qp, do_lhs,true, kluge,elem);
        	}

} // points

// user has updated mesh or something
void new_iteration()
{
m_eqn_factors = MyBlock(1);

}

void set_transuranics(SpecialElements * qmap) { m_transuranics=qmap; }

typedef  mjm_libmesh_util::NandE Path;
// this can be done well before they are needed... 
// this is called on each mesh update, use it for an iteration indicator now
void update_transuranics_from_nevec()
{
new_iteration();
// a lot of this was figured out during path tracing, kind of dumb
// to do it here agin
// use m_surfaces[1] and 2, 
// right now turning OFF transuraniacs enforces BC's on electrodes
// and it is easier just to let this run 
//if (!true){static bool crap=false; if (!crap) {  MM_ERR(" skipping transuranics  ") crap=true; }   return ;  }
if (suppress_transuraniacs()) { MM_ONCE("skipping transuraniacs ", return; ) }
MM_ONCE(" seem to be updating transuraniacs so that is enabled but later models may be disabled  ",)
if (m_surfaces==0){  MM_ERR(" no  surfaces ")  return ;  }
if (m_surfaces->size()<3){  MM_ERR(" notenough  surfaces ") return ; } 

SpecialElements & semap=(*m_transuranics);
semap.clear();
for (IdxTy surface=1; surface<3; ++surface)
{
//ZZtypedef mjm_special_qrule MqTy;
const Path  & path=(*m_surfaces)[surface];
typedef std::map<const Elem *,IdxTy> Vim;
Vim visited;
Vim corners;
const bool all_same=false;
const IdxTy sz=path.size();
for (IdxTy i=0; i<sz; ++i)
{
    //const Node * n=path.nv[i];
    const Elem * e=path.ev[i];
	// we may have already marked this due to being part of a singular group
	if (semap.find(e)!=semap.end()) continue; 
   // if (visited.find(e)!=visited.end()) continue;
   // visited[e]+=1;
// determine if this is a singular boundary 
	if (!path.singular(i))	
	{ semap[e]= MqTy(); continue; } 
// something is difficult about this point make a special thing
	if (all_same) { semap[e]= MqTy();  continue; } 
	corners[e]=i;

} // fori 
if (all_same) continue;
// make instance for the points, could also aggregate others for reverse p adaptatation.
make_corner_cases(corners,path);

} // surface
}
template <class Tc> 
void make_corner_cases(Tc  & corners, const Path & path)
{
typedef std::map<const Elem *, IdxTy> Ev;
typedef std::map<const Node * , Ev> NtoE;
typedef std::map<const Node * , IdxTy> NoMap;

SpecialElements & semap=(*m_transuranics);
NtoE groups;
// get the elements around each singular node
for (auto ii=corners.begin(); ii!=corners.end(); ++ii) 
{
const Elem * e=(*ii).first;
const IdxTy sz=(*e).n_nodes();
// this is a collection of all elements that need to be modified due
// to the singular node key
for (IdxTy i=0; i<sz; ++i) 
{
	const Node * n=(*e).get_node(i);
	//if (path.is_singular((*e).get_node(i))) groups[(*e).get_node(i)][e]=1; 
	if (path.is_singular(n)) 
	{
		groups[n][e]=1; 
		// can also add neighbors now
		const IdxTy ized=(*ii).second;
		const IdxTy psz=path.ev.size();
		const IdxTy code=path.bc_codes[ized]; // try to avoid wrap of bc
		if (code !=0)
		{ MM_ERR(" corner code is not zero "<<code) }
		IdxTy iplus=(int(psz+ized)+1)%psz;
		IdxTy iminus=(int(psz+ized)-1)%psz;
		const IdxTy codeminus=path.bc_codes[iminus]; // try to avoid wrap of bc
		const IdxTy codeplus=path.bc_codes[iplus]; // try to avoid wrap of bc

		// this could easily wrap around corners, need to restrict for that 
		for (int  del=-30; del<30; ++del) 
		{
			IdxTy imore=(int(psz+ized)+del)%psz;
			const IdxTy codedel=path.bc_codes[imore]; // try to avoid wrap of bc
			if (del<0) if (codedel!=codeminus) continue;
			if (del>0) if (codedel!=codeplus) continue;

			const Elem * en=path.ev[imore];
    		groups[n][en]=1;
			// now erode out a level this should probably be done in steps
			// 4 levels takes forever lol  
			// not sure that acutally works right 
			// This also needs to be restricted on wrap
			mjm_libmesh_util::n_touchers(groups[n],en,3);
			//MM_ERR(" dildate groud size is now "<<groups[n].size())
/*			typedef std::set<const Elem * > ElSet;
			ElSet toucherset;
			// no idea how slow this is but more like to be robust than my code
			en->find_point_neighbors(toucherset);
			for ( auto ii=toucherset.begin(); ii!=toucherset.end(); ++ii) 
			{ groups[n][(*ii)]=1; }
			*/

		}
	}

}
// for now just analyze do not try to change 
//	semap[e]= MqTy();  

}
// fit to analytical part of integral

for (auto ii=groups.begin(); ii!=groups.end(); ++ii) 
{
// this node is part of each element in g.
const Node * n=(*ii).first;
Ev & g=(*ii).second;
// for the 5 non boundary nodes, fit an
// analytical approximation to integrate thne 
// just do the residuals. The shape functions
// 
using namespace  mjm_libmesh_util;
NoMap nmap,allnmap;
// get all the nodes with no bc's
get_all_nodes(nmap,g,(m_mla->mesh),0); 
get_all_nodes(allnmap,g,(m_mla->mesh),1); 
const IdxTy nmsz=nmap.size();
const IdxTy allnmsz=allnmap.size();
//if (nmsz!=5) 
{ MM_ERR(" nmap sizes "<<nmsz<<" and "<<allnmsz<<" with "<<g.size()) } 
// we can either fit here or in the assembly routine with given solution
// for now just punt
MyBlock coefs(nmsz);
MyBlock residuals(nmsz);
// right now this just prints the fit nothing more 
MqTy & mq=semap[((*g.begin()).first)];
mq.fit_analytic_part( coefs,  residuals ,n,nmap, *solution) ;
//semap[((*g.begin()).first)].fit_analytic_part( coefs,  residuals ,n,nmap, *solution) ;
// give all members of the group a custom qrule 
for (auto jj=g.begin(); jj!=g.end(); ++jj) semap[(*jj).first]=mq;


} // groups
// make MqTy instance for each

// once we have special MqTy instances for the special
// points add them to the transuranics map 


}
template < class Tx, class Tyz,class Tc,class Ts >
        IdxTy  assemble_points_sides( Tx & lhs_block, Tyz & rhs_block,  const Ts & sol
		, const Tc * elem, const bool do_lhs
	)// , const bool is_side, const IdxTy side)
{
		IdxTy sides=0;
		for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
        	if (elem->neighbor(kluge)==libmesh_nullptr)
        	{    ++sides;
        		const IdxTy qpssz=(*m_mla).fes.qrule.n_points();
        		for (IdxTy qp=0; qp<qpssz; qp++)
         			assemble_some( lhs_block, rhs_block, 
						(*m_mla).fes, sol, qp, do_lhs,true, kluge,elem);
        	}
	return sides;
}

// this needs to return the source or figure out otherwise
IdxTy return_surface_code(const Elem * elem)
{
	// this should have source info telling us ty-e 
 	BoundaryInfo & bi = m_mla->mesh.get_boundary_info();
	// this should be a bit string of missing sides, bit 0 is side 0 etc
	const IdxTy side_code=mjm_libmesh_util::bc_code(elem,bi,true);
	IdxTy bsrc=0; 
	if (m_surfaces!=0) 
	{
		for (IdxTy i=0; i<m_surfaces->size(); ++i)
		{ 
			if ((((*m_surfaces)[i])).may_have(elem)) 
				{bsrc=i+1; break; }
		} //i 	
	}else {MM_ONCE(" could get the boundary src from bi",) }

return side_code|(bsrc<<16);
}

const bool use_trans(const Elem * elem, const bool skip)
{
	if (skip) return false;
	if (m_transuranics==0)  return false;
	return  (m_transuranics->find(elem)!=m_transuranics->end());
} // use_trans



template <class Tvec>
  void assemble_boundary_elem
	(const Elem *e, const Tvec & sol, Tvec  * R, Tmjm_sparse * J, const IdxTy surface)
{
	const bool check_if_transuranics_exist= !suppress_transuraniacs(); 
	MqTy mq;
	MqTy * pmq=&mq;
	if (use_trans(e,!check_if_transuranics_exist))
	{
		{ MM_ONCE(" at least one special transuraniac boundaary exists ",) }
		pmq=  &(*m_transuranics)[e];
	}
	// actually this should be able to do it all pass surface code... 
	(*pmq).assemble_boundary_elem( e,sol,R,J,*cond(),surface);
}

template < class Tx, class Tyz,class Tc,class Ts, class Td >
        void assemble_points( Tx & lhs_block, Tyz & rhs_block,  const Ts & sol
		, const Tc * elem, const bool do_lhs
	, Td & my_dof_block // , const IdxTy bcode, const IdxTy sidex)
	)// , const bool is_side, const IdxTy side)
{
// these may update even if not requested yet
	const bool check_if_transuranics_exist= !suppress_transuraniacs(); 
	const bool add_bc_to_specials=false;
	// there is no other code 
	// do all of these here even if not special 
	//const bool do_all_trans=true;
	//const bool disable_bc_on_interiors =(m_iteration!=0);
	const bool disable_bc_on_interiors =!impose_electrode_bc();
	// (surface<<16)|sides 
	const IdxTy bsurface=return_surface_code(elem);
	bool is_special=false;
	const IdxTy sidex=1;
	// only need if bsurface!=0
	const bool do_bcs=false; // ((sidex==0)&&(((bsurface>(1<<15)))||((bsurface!=0)&&!disable_bc_on_interiors)));
   	const bool use_finer_bc = false;
	MqTy mq;
	MqTy * pmq=&mq;
	if (use_trans(elem,!check_if_transuranics_exist))
	{
		{ MM_ONCE(" at least one special transuraniac exists ",) }
		pmq=  &(*m_transuranics)[elem];
		is_special=true;
	}
	// actually this should be able to do it all pass surface code... 
	(*pmq).assemble_given_block(lhs_block, rhs_block, (*m_mla).fev 
		, do_lhs,sol,*cond(),elem,bsurface,sidex);
	if (1==1)
	{
		MM_ONCE(" removing any bc assemble from nonliner ",)
		return; 
	}
	// this ignores both the electrodes AND the periphery
	// probably just better to check if it is on the periphery 
	// but not a big deal 
 	if (do_bcs){  
// FIXME TODO this is using the old gross appraoximationas and libmesh code
     	if (use_finer_bc) 
		{
			(*pmq).assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);}
     	else
		{
			 assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);
		}
 		{ MM_ONCE(" assembling at least one bc "<<m_iteration,) }
	} else {MM_ONCE(" skipping at least one bc "<<is_special,) } 
} // assemble_points ( new )


///////////////////////////////////////////////

template < class Tx, class Tyz,class Tc,class Ts >
        void assemble_points_hybrid( Tx & lhs_block, Tyz & rhs_block,  const Ts & sol
		, const Tc * elem, const bool do_lhs
	)// , const bool is_side, const IdxTy side)
{
 	//BoundaryInfo & bi = m_mla->mesh.get_boundary_info();
// these may update even if not requested yet
const bool check_if_transuranics_exist=!suppress_transuraniacs(); //true;
const bool add_bc_to_specials=false;
// do all of these here even if not special 
const bool do_all_trans=true;
//const bool disable_bc_on_interiors =(m_iteration!=0);
const bool disable_bc_on_interiors =!impose_electrode_bc();

if (check_if_transuranics_exist){ if (m_transuranics!=0) 
{
	if (m_transuranics->find(elem)!=m_transuranics->end())
	{
		// ths should punt to the m_mla for what to do but for now
		// I guess just hard code it... 
		// each element has its own INSTANCE it may as well be inited
		// that way lol. 
//		MqTy mq_dummy;
		{ MM_ONCE(" at least one special transuraniac exists ",) }
		MqTy & mq=  (*m_transuranics)[elem];
	//MM_ERR(" asseembling special asif normal i no bc ")
// 		mq_dummy.assemble_given_block(lhs_block, rhs_block, (*m_mla).fev , do_lhs,sol,*cond(),elem);
 		mq.assemble_given_block(lhs_block, rhs_block, (*m_mla).fev , do_lhs,sol,*cond(),elem);
		// enforcing on the periphery will help make solution more sane 
		// this is using lib mesh not the new code... 
		// and value and grads are gross
  		if (add_bc_to_specials)      
			mq.assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);
			//assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);
		else { MM_ONCE(" assembling SPECIAL WITHOUT bc",) }
		return;
	}
}} else {MM_ONCE("disabled check for special transuraniacs",) }
if (do_all_trans)
{
		MqTy mq;
 		mq.assemble_given_block(lhs_block, rhs_block, (*m_mla).fev , do_lhs,sol,*cond(),elem);
// this ignores both the electrodes AND the periphery
// probably just better to check if it is on the periphery but not a big deal 
bool do_bcs=!false;
if (disable_bc_on_interiors){ 
	if (m_surfaces!=0) if (m_surfaces->size()>0) 
	{ if (m_surfaces->size()>1) if ((((*m_surfaces)[1])).may_have(elem)) do_bcs=false; 
	 if (m_surfaces->size()>2) if ((((*m_surfaces)[2])).may_have(elem)) do_bcs=false; 
	}
	} else {MM_ONCE(" surface BC DISABLED "<<m_iteration,) }
 	if (do_bcs){  
// FIXME TODO this is using the old gross appraoximationas and libmesh code
     	mq.assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);
     	//assemble_points_sides( lhs_block, rhs_block, sol, elem,  do_lhs);
 		{ MM_ONCE(" assembling at least one bc "<<m_iteration,) }
	}else { MM_ONCE(" assembling WITHOUT electrode bc",) }
	return; 
}
{MM_ONCE(" fell through trans doing at least one element as normal ",) }
// residuals don't look different but more R^2 drop on solve?
const bool dups_OK=!false;
/*
        const IdxTy qpsz=(*m_mla).fev.qrule.n_points();
        for (IdxTy qp=0; qp<qpsz; qp++)
        { assemble_some( lhs_block, rhs_block
 , (*m_mla).fev, sol,qp, do_lhs,false, 0,elem); } //qp
*/
//		const bool ignore_bc=!true;
//		if ( ignore_bc) return; 
		IdxTy sides=0;
		for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
        	if (elem->neighbor(kluge)==libmesh_nullptr)
        	{    ++sides;
        		const IdxTy qpssz=(*m_mla).fes.qrule.n_points();
        		for (IdxTy qp=0; qp<qpssz; qp++)
         			assemble_some( lhs_block, rhs_block, 
						(*m_mla).fes, sol, qp, do_lhs,true, kluge,elem);
        	}


if (dups_OK|| (sides==0)) {
        const IdxTy qpsz=(*m_mla).fev.qrule.n_points();
        for (IdxTy qp=0; qp<qpsz; qp++)
        { assemble_some( lhs_block, rhs_block
 , (*m_mla).fev, sol,qp, do_lhs,false, 0,elem); } //qp
}


} // points
template <class Tv> 
const bool Val(const Tv & x) { return mjm_libmesh_util::Val(x); } 
// adds a conservation equation that really needs to be scaled and
// is obliterated with home made dirichlets that wipe out the equation.
// zres is the residue at sol, zjake is the set of deriviatves. 
template <class Tj, class Tr, class Tz>
		void add_zero(Tj & J,Tr & R, Tz & zjake, const Real & zres
		,const IdxTy _zdof, const Real scale=1, const IdxTy flags=0)
{
bool do_zdof=((Val(_zdof))&&(flags==0));
IdxTy cntmax=0;
if (R) R->close();
if (J) J->close();
if ((flags&1)!=0)  cntmax=_zdof;
if (do_zdof)
{
MM_ERR(" only adding one zero at eqn "<<_zdof)
if (R) R->add(_zdof,zres);
if (J) { for (auto ii=zjake.begin(); ii!=zjake.end(); ++ii) 
	{ J->add(_zdof,(*ii).first,(*ii).second); } // ii 
}
return ;
}
IdxTy cnt=0;
for (auto jj=zjake.begin(); jj!=zjake.end(); ++jj) 
{
if (Val(cntmax)) if (cnt>=cntmax) break; 
if (R) R->add((*jj).first,zres*scale);
if (J) { for (auto ii=zjake.begin(); ii!=zjake.end(); ++ii) 
	{ J->add((*jj).first,(*ii).first,scale*(*ii).second); } // ii 

}
++cnt;
} 

MM_ERR(" added "<<cnt<<" zeroes  with eqn dof  "<<_zdof)

}


template <class Te , class Tp>
const Real eqn_at_p( const Te & eq, const Tp & sol)
{
	Real x=0;
	for (auto ii=eq.begin(); ii!=eq.end(); ++ii)
 		{ x+=(*ii).second*(sol((*ii).first)); }
	return x; 
}
template <class Tv, class Ts,class Tfudd> 
	void zero_assembly( Tv & residual, Ts * jacobian, const Tfudd & sol)
{
	(*m_cm).mark("assbegin");

	if (m_zero==0) MM_ERR(" no zero equn")
	if (m_surfaces==0) MM_ERR(" no surfaces ")
	if ((m_zero==0)||(m_surfaces==0)) return; 
	if (m_surfaces->size()<3) MM_ERR(" notenough  surfaces ")
 	BoundaryInfo & bi = m_mla->mesh.get_boundary_info();
	auto solp=sol.clone();
	m_zero->clear();
	ZeroEqn zero_del;
//	ZeroEqn  m_zero_dup;
    using mjm_libmesh_util::CurrentMatch;
    using mjm_libmesh_util::CurrentMatchOld;
    using mjm_libmesh_util::CurrentMatchAvg;
    // using mjm_libmesh_util::CurrentMatchCompare;
    // this relies on boundary first etc 
    CurrentMatchAvg(*m_zero, (*m_surfaces)[1],(*m_surfaces)[2],bi,
		m_mla->di, (*m_mla).nvar, m_mla,sol);
//    CurrentMatch(m_zero_dup, (*m_surfaces)[1],(*m_surfaces)[2],bi,
//		m_mla->di, (*m_mla).nvar, m_mla,sol);

//	const IdxTy iiiii=m_zero_dup.size();
//	MM_ERR(" print "<<iiiii)
//	CurrentMatchCompare(*m_zero,m_zero_dup);	
	// we picked residual as rhs minus lhs, stick with it for now 
	Real r0=-eqn_at_p(*m_zero,sol);
	residual=r0;
	if (jacobian==0) return; 
	const IdxTy sz=(*m_zero).size();
	for (auto ii=(*m_zero).begin(); ii!=(*m_zero).end(); ++ii)
	{
		const IdxTy dof=(*ii).first;	
		Real xzed=(*solp)(dof);
		Real deltaeff=xzed*.000000001;
		if (deltaeff<0) deltaeff=-deltaeff;
		if (deltaeff<1e-28) deltaeff=1e-28;
		solp->add(dof,deltaeff);
    	CurrentMatchAvg(zero_del, (*m_surfaces)[1],(*m_surfaces)[2],bi,
			m_mla->di, (*m_mla).nvar, m_mla,*solp);
		Real rd=-eqn_at_p(zero_del,*solp);
		Real dr_ddof=(rd-r0)/deltaeff;
		(*jacobian)[dof]+=dr_ddof;		
		solp->add(dof,-deltaeff);

	} // ii 


} // zero_assembly 

template <class Tb> void DumpRHS(const Tb & rhs)
{    
StrTy ival=(*m_cm).format_entry("iterations");
// rhssiterations=2
if (ival.length()==11) ival=StrTy("iterations=1");
mjm_libmesh_util::DumpVector(std::cerr, m_mla->mesh, rhs, StrTy("rhss")+ival,true);
}

template <class Tb> void DumpRHS(const Tb & rhs_zed,const Tb & rhs)
{    
StrTy ival=(*m_cm).format_entry("iterations");
// rhssiterations=2
if (ival.length()==11) ival=StrTy("iterations=1");
mjm_libmesh_util::DumpVectorPair(std::cerr, m_mla->mesh, rhs_zed,rhs, StrTy("rhss")+ival,true);
}


class assembly_branches
{
public:
// this is called on each assembly so it can not 
// re-load a config file  each time and need to cache lol 
template <class Tx, class Tz>
assembly_branches(const Tz * J, const Tx * R):
	ddump(!false)
 	,doing_j(true&&(J)) 
 	,doing_residual((true)&&(R))
	,dump_jacobian((ddump)&&doing_j)
	//,dump_rhs((true)&&doing_j)
	,dump_rhs((ddump)&&doing_residual) // before this would dump zed and non-zed
	,dump_residuals_here((ddump)&&doing_residual)
	,doing_linear_ass((!doing_j)&&(!doing_residual))
 	,do_lhs(true )
 	,do_mjm_sparse(true)
	,do_zero_add((!true)&&!doing_linear_ass)
	// this does have an effect 
	,do_eqn_factors(!false) // putting this in disabled changed results? 
//	,dump_rhs_with_mesh(true)
	,dump_sol(ddump)
	,zero_factor(1e-6)
{
// deltaeff???
// MM_ONCE(" putting in disabled normalization changed things wtf",)
}


class dumps
{
public:
dumps(): jacobian(false),rhs(false),residuals(false),sol(false) {}
dumps(const ReadWriteMap & cf ): 
	jacobian(cf.get_bool("dump_jacobian",true))
	,rhs(cf.get_bool("dump_rhs",true))
	,residuals(cf.get_bool("dump_residuals",true))
	,sol(cf.get_bool("dump_solution",true))
{}
bool jacobian,rhs,residuals,sol;

}; // dumps

template <class Tx, class Tz>
assembly_branches(const dumps & d, const Tz * J, const Tx * R):
	ddump(false)
 	,doing_j(true&&(J)) 
 	,doing_residual((true)&&(R))
	,dump_jacobian((d.jacobian)&&doing_j)
	//,dump_rhs((true)&&doing_j)
	,dump_rhs((d.rhs)&&doing_residual) // before this would dump zed and non-zed
	,dump_residuals_here((d.residuals)&&doing_residual)
	,doing_linear_ass((!doing_j)&&(!doing_residual))
 	,do_lhs(true )
 	,do_mjm_sparse(true)
	,do_zero_add((!true)&&!doing_linear_ass)
	,do_eqn_factors(false) // putting this in disabled changed results? 
//	,dump_rhs_with_mesh(true)
	,dump_sol(d.sol)
	,zero_factor(1e-6)
{
// deltaeff???
// MM_ONCE(" putting in disabled normalization changed things wtf",)
}


/*

template <class Tx, class Tz>
assembly_branches(const Tz * J, const Tx * R,const ReadWriteMap & cf):
// this does not probably work at all lol 
 	doing_j(cf.get_bool("doing_j",true)&&(J)) 
 	,doing_residual(cf.get_bool("doing_residual",true)&&(R))
	,dump_jacobian(cf.get_bool("dump_jacobian",true)&&doing_j)
	//,dump_rhs((true)&&doing_j)
	,dump_rhs((true)&&doing_residual) // before this would dump zed and non-zed
 	,dump_rhs(cf.get_bool("dump_rhs",true)&&(doing_residual))
	,dump_residuals_here((true)&&doing_residual)
	,doing_linear_ass((!doing_j)&&(!doing_residual))
 	,do_lhs(true )
 	,do_mjm_sparse(true)
	,do_zero_add((!true)&&!doing_linear_ass)
	,do_eqn_factors(false) // putting this in disabled changed results? 
//	,dump_rhs_with_mesh(true)
	,dump_sol(true)
{
// deltaeff???
// MM_ONCE(" putting in disabled normalization changed things wtf",)
}

*/

	const bool ddump;
 	const bool doing_j; // =true&&(J); 
 	const bool doing_residual; // =(true)&&(R); 
	const bool dump_jacobian;// =(!true)&&doing_j;
	//const bool dump_rhs=true&&doing_residual;
	//const bool dump_rhs=true&&doing_residual&&doing_j;
	const bool dump_rhs; // =(!true)&&doing_j;
	const bool dump_residuals_here; // this is also abailabe in later code... 
	const bool doing_linear_ass; // =(!doing_j)&&(!doing_residual);
 	const bool do_lhs; // =true; 
 	const bool do_mjm_sparse; //=true; 
	const bool do_zero_add;//(!true)&&!ab.doing_linear_ass;
	const bool do_eqn_factors; // (false)
//	const bool dump_rhs_with_mesh;
	const bool dump_sol;
	const Real zero_factor;
}; // assembly_branches




template <class Tx, class Tyy>
void check_node_config(Tx & elem, Tyy & mla)
{
		const IdxTy nnodes_test=elem->n_nodes();
		const IdxTy nvdofs_test=mla.di.n_var_dofs;
		if (nvdofs_test!=4) {MM_ERR(" non 4 dofs "<<nvdofs_test)}
		if (nvdofs_test!=4) {MM_ONCE("  see errors non 4 dofs "<<nvdofs_test,)}
		if (nnodes_test!=4) {MM_ERR(" non 4 nodes "<<nnodes_test)}
		for (IdxTy in=0; in<nnodes_test; ++in)
			if (!elem->is_vertex(in))  { 
			{ MM_ONCE(" at least one non vertex found "<<in<<" of "<<nnodes_test,) }
			MM_ERR(" non vertex "<<in<<" of "<<nnodes_test)}
		// this order is not right see the addLHS2D for example 
} // check_node_config
template <class Tx,class Tyy, class Tz>
void make_deltaeff( Real & deltaeff,Real & my_xzed, Tx & my_sol, const IdxTy dofj,
Tyy & lhs_block_zed, Tz & rhs_block_zed, const IdxTy var, const IdxTy j )
{
	my_xzed=my_sol(dofj);
	const Real deltafac=1.3e-9;
//	const bool is_u=((dofj&&3)==0);
	deltaeff=my_xzed*deltafac; // .00000013; // *1e-9;
	//if (is_u) deltaeff=my_xzed*deltafac; // .00000013; // *1e-9;
//	else deltaeff=my_xzed*deltafac*1e-4; // .00000013; // *1e-9;
	MM_ONCE("deltaeff is now "<<deltafac<<" and now "<<deltaeff,)
	if (deltaeff<0) deltaeff=-deltaeff;
	const Real deltamin=1e-28;
	//if (deltaeff<1e-25) deltaeff=1e-25;
	if (deltaeff<deltamin){
		//const Real x1=lhs_block_zed(var,j,var,j);		
		const Real x1=lhs_block_zed(j,var,j,var);		
		//const Real x2=rhs_block_zed(var,j);	
		const Real x2=rhs_block_zed(j,var);	

		if( 1e-16<(x1*x1)) deltaeff=x2/x1*1e-10;
		MM_ERR(" invoking delta min"<<x1<<" "<<x2<<" "<<deltaeff)
		if (deltaeff<0) deltaeff=-deltaeff;	
		if (deltaeff<deltamin) 	 deltaeff=deltamin;
	}
	my_sol(dofj)+=deltaeff;

} // make_delta_eff
void	set_blocks(MyBlock & rhs_block,MyBlock & rhs_block_zed
		, MyBlock &lhs_block,MyBlock & lhs_block_zed
		, const IdxTy nvars,const IdxTy nvdofs)
{


}

template <class Tv, class Ts,class Tfudd> void assembly( Tv * R, Ts * J, const Tfudd & sol)
{
	(*m_cm).mark("assbegin");
    static int ass_count=0;
    MM_MSG(" ass_count="<<ass_count)
    ++ass_count;
	typedef std::vector<const Elem *> BoundElems;
	BoundElems boundary_elements;
	// this can not be static ...  but could load from m_map
	assembly_branches ab(J,R);
	const bool find_eqn_factors
		=ab.do_eqn_factors&&(m_eqn_factors.size()<2)&&!ab.doing_linear_ass;
	typedef mjm_sparse_matrix<Real> Tmjm;
	if (ab.doing_linear_ass)
	{
		MM_ONCE(" doing a linear assembly ",)
		matrix->zero();
		rhs->zero();
	}
	MLATy &   mla = * m_mla;
    (*m_cm).set("nodes", mla.mesh.n_nodes());
	Tmjm my_sparse;
	const IdxTy solsz=sol.size();
	MyBlock rhss(ab.dump_rhs?solsz:1), rhss_zed(ab.dump_rhs?solsz:1);
	MyBlock my_sol(solsz), my_sol_zed(solsz)
		//,my_res((ab.doing_residual||find_eqn_factors)?(*R).size():1);
		,my_res((ab.doing_residual||ab.do_eqn_factors)?solsz:1);
	if (ab.doing_j||find_eqn_factors) {my_sol.copy_from(sol); }
	my_sol_zed.copy_from(sol); 
    MeshCItor      el     = mla.mesh.active_local_elements_begin();
    const MeshCItor end_el = mla.mesh.active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
        const Elem * elem = *el;
        mla.re_init(elem);
		const IdxTy nvdofs=4;
		const IdxTy nvars=mla.nvar;
	// recalculated many teims 
	const IdxTy bsurface=return_surface_code(elem);
	if (bsurface!=0)
	{
		// surface code can be recompuated later
		boundary_elements.push_back(elem);

	}
		const IdxTy ntot=nvdofs*nvars;
		check_node_config(elem,mla);
		MyIdxBlock my_dof_block(nvdofs,nvars);
		typedef MyBlock::location_type LocTy;
		// this could just be done once since all the dims are the same. 
		//set_blocks(rhs_block,rhs_block_zed, lhs_block,lhs_block_zed,nvars,nvdofs);
		// these indexes are backwards
		const IdxTy dms[]={nvdofs,nvars,nvdofs,nvars};
		const IdxTy dmr[]={nvdofs,nvars};
		LocTy rhsdims=MyBlock::loc(dmr,2);
		LocTy lhsdims=MyBlock::loc(dms,4);
		MyBlock rhs_block(rhsdims),rhs_block_zed(rhsdims);
		MyBlock lhs_block(lhsdims),lhs_block_zed(lhsdims);
		mjm_libmesh_util::get_dof_mat(my_dof_block,elem, nvars);
 		assemble_points( lhs_block_zed, rhs_block_zed, my_sol_zed, elem, ab.do_lhs,my_dof_block);  
		if (ab.doing_j||find_eqn_factors)
		{
			for (IdxTy var=0; var<mla.nvar ; ++var)
				for (IdxTy j=0; j<nvdofs ; ++j)
				{// dof
				rhs_block.zero();
				lhs_block.zero();
                //const IdxTy dofj=my_dof_block(var,j);  
                const IdxTy dofj=my_dof_block(j,var);  
				Real xzed=0,deltaeff=0; // =my_sol(dofj);
				make_deltaeff(deltaeff,xzed,my_sol, dofj,lhs_block_zed,rhs_block_zed, var,j);
				assemble_points( lhs_block, rhs_block, my_sol, elem
					, ab.do_lhs, my_dof_block);  
       			mla.di.my_jake(my_sparse,my_sol_zed,lhs_block_zed ,rhs_block_zed,&my_sol, 
					lhs_block, rhs_block, 1.0/deltaeff,dofj  , my_dof_block);
				my_sol(dofj)=xzed;
				} //for dof and var 
		} // doing_j

     	if (ab.doing_residual||find_eqn_factors) 
			//mla.di.residual(R,my_sol_zed,lhs_block_zed,rhs_block_zed ,my_dof_block);
			mla.di.my_residual(my_res,my_sol_zed,lhs_block_zed,rhs_block_zed ,my_dof_block);
		if (ab.dump_rhs) mla.di.add_block(rhss_zed,rhs_block_zed,my_dof_block);
//		if (ab.dump_rhs) mla.di.add_block(rhss,rhs_block,my_dof_block);

		if (ab.doing_linear_ass)
			mla.di.assemble_jr (matrix, rhs, lhs_block_zed, rhs_block_zed ,my_dof_block);

	} // el element loop
	if (ab.doing_j&&ab.doing_residual)
		MM_MSG("NL doing both now ...")
	if (find_eqn_factors)
	{
		m_eqn_factors= MyBlock(solsz);
		for (IdxTy i=0; i<solsz; ++i)
		{
			//const Real max_row=my_sparse(i,i); // my_sparse.max_abs_row(i);
			//const Real max_row= my_res(i); // my_sparse.max_abs_row(i);
			const Real max_row=  my_sparse.max_abs_row(i);
			if (max_row!=0) m_eqn_factors(i)=1.0/max_row;
			else {
				MM_ERR(" row max is zero "<<i) 	
				 m_eqn_factors(i)=1.0;
			}
		//	if ((i&3)==0) m_eqn_factors(i)*=1e-5;
//			if ((i&3)==2) m_eqn_factors(i)*=100;
			MM_ONCE(" finding first eqn factor "<<m_eqn_factors(i),)
		}
	}
//	if (!false) 
	if (ab.do_eqn_factors)
	{
	  { MM_ONCE(" do_eqn_factors is true normalize  ",) }
		MySparse::sequential_iterator si(my_sparse);
		while (si)
		{
			const IdxTy eqn=si.i();
			const Real fac=m_eqn_factors(eqn);
			my_sparse(eqn,si.j())*=fac;
			if (eqn>=my_res.size()) 
				{ MM_ERR(" residue size is wrong "<<eqn<<" vs "<<my_res.size()) }
			MM_ONCE(" using first eqn factor "<<fac,)
			++si;
		}
		if (m_eqn_factors.size()<my_res.size())
			MM_ERR(" not enough factors for residue "<<m_eqn_factors.size()<<" bs " << my_res.size())
		for (IdxTy i=0; i<my_res.size(); ++i)  my_res(i)*=m_eqn_factors(i);

	} else { MM_ONCE(" do_eqn_factors is false ",) }


  //void handle_boundaries(const Tvec & sol, TVec  * R, Tmjm_sparse * J)
   handle_boundaries(my_sol_zed, R?&my_res:0, J?&my_sparse:0);
	// these are things like getting the particle fluxes to work out
   handle_boundaries(boundary_elements,my_sol_zed, R?&my_res:0, J?&my_sparse:0);
// need to fix everything before this...	
	if (ab.doing_j) {  my_sparse.to_libmesh(*J); }
	if (ab.doing_residual) {  my_res.copy_to_libmesh(*R); }

	if (ab.do_zero_add) 
	{ 
		MM_ONCE(" adding the zero eqns for current eaulity",) 
		Real zres=0;
		ZeroEqn zjake;
		zero_assembly( zres, &zjake, sol);
		IdxTy zdof=~0;	
		// this seems to work better than 1e-3
		add_zero(J,R,zjake,zres,zdof,ab.zero_factor,1);
	}
	else MM_ONCE(" zero add disable",)
	// these do not include the zero add results
	do_dump(ab.dump_residuals_here,my_res,"reszed");
	do_dump(ab.dump_rhs,rhss_zed,"rhszed");
// dumping this is misleading.. 
			//mjm_libmesh_util::DumpVector(std::cerr, m_mla->mesh, rhss, StrTy("rhss")+ival);
	do_dump(ab.dump_sol,my_sol_zed,"solzed");
	
	do_dump(ab.dump_jacobian,my_sparse,"jacobian");
	// it may be more beneficial to dump it this way too 
	if (!false) if (ab.dump_jacobian)
	{
		MM_ONCE(" jacobian dumping enabled here ",) 
//template < class Os,class Tv > void dump_terms(Os & os, const IdxTy flags, const Tv & vector, const StrTy & label)
   		StrTy ival=(*m_cm).format_entry("iterations");
		if (ival.length()==11) ival=StrTy("iterations=0");
		//mjm_libmesh_util::DumpVector(std::cerr, m_mla->mesh, x, StrTy(label)+ival);
		my_sparse.dump_terms(std::cerr,0 ,my_sol_zed,StrTy("jacobianterms")+ival);
// this is so fudding slow.. 
//		mjm_libmesh_util::DumpSparseMatrix(std::cerr,*J,"jacobian",1);
	}

	assembly_exit_cm();

} // assembly

template <class To, class Tl>
void do_dump(const bool do_or_no, const To & x, const Tl&  label)
{
	if (!do_or_no) return;
	MM_ONCE("at least one dumpers "<<label<<" is generatingmasive output ",)
	// need to add zero ... 
   	StrTy ival=(*m_cm).format_entry("iterations");
	if (ival.length()==11) ival=StrTy("iterations=0");
	mjm_libmesh_util::DumpVector(std::cerr, m_mla->mesh, x, StrTy(label)+ival);
}




void assembly_exit_cm()
{
    (*m_cm).lap("ass","assbegin");
    (*m_cm).cum("asstotal","assbegin");
    (*m_cm).dump("EndAssemble", std::cerr); std::cerr<<std::endl;
    (*m_cm).inc("iterations");
    (*m_cm).clear("currentnew");
    (*m_cm).clear("areanew");

}


// so now what to do with m_zero lol 
//if (m_zero!=0)
template <class Tv, class Ts,class Tfudd> void assembly_converting_to_mm( Tv * R, Ts * J, const Tfudd & sol)
{
	(*m_cm).mark("assbegin");
    static int ass_count=0;
    MM_MSG(" ass_count="<<ass_count)
    ++ass_count;
 	const bool doing_j=true&&(J); 
 	const bool doing_residual=(true)&&(R); 
	const bool dump_jacobian=(!true)&&doing_j;
	//const bool dump_rhs=true&&doing_residual;
	//const bool dump_rhs=true&&doing_residual&&doing_j;
	const bool dump_rhs=(!true)&&doing_j;
	const bool doing_linear_ass=(!doing_j)&&(!doing_residual);
 	const bool do_lhs=true; 
 	const bool do_mjm_sparse=true; 
	typedef mjm_sparse_matrix<Real> Tmjm;
	Tmjm my_sparse;
	if (doing_linear_ass)
	{
		MM_ONCE(" doing a linear assembly ",)
		matrix->zero();
		rhs->zero();
	}
	MLATy &   mla = * m_mla;
    (*m_cm).set("nodes", mla.mesh.n_nodes());
//	using mjm_libmesh_util::MakeXandGradXVector;
	auto solp=sol.clone();
	//MyBlock rhss(dump_rhs?(*R).size():1);
	MyBlock rhss(dump_rhs?(*solp).size():1);
	MyBlock rhss_zed(dump_rhs?(*solp).size():1);
	MyBlock my_sol((*solp).size());
	MyBlock my_sol_zed((*solp).size());
	if (do_mjm_sparse&&doing_j) {my_sol.copy_from(*solp); }
	if (do_mjm_sparse) {my_sol_zed.copy_from(*solp); }
 	//BoundaryInfo & bi = mla.mesh.get_boundary_info();
    MeshCItor      el     = mla.mesh.active_local_elements_begin();
    const MeshCItor end_el = mla.mesh.active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
        const Elem * elem = *el;
        //if (elem->has_children()) continue; 
        mla.re_init(elem);
		const IdxTy nvdofs_test=mla.di.n_var_dofs;
		if (nvdofs_test!=4) {MM_ERR(" non 4 dofs "<<nvdofs_test)}
		if (nvdofs_test!=4) {MM_ONCE("  see errors non 4 dofs "<<nvdofs_test,)}
		const IdxTy nvdofs=4;
		const IdxTy nvars=mla.nvar;
		const IdxTy ntot=nvdofs*nvars;
		const IdxTy nnodes_test=elem->n_nodes();
		if (nnodes_test!=4) {MM_ERR(" non 4 nodes "<<nnodes_test)}
		for (IdxTy in=0; in<nnodes_test; ++in)
			if (!elem->is_vertex(in))  { 
			{ MM_ONCE(" at least one non vertex found "<<in<<" of "<<nnodes_test,) }
			MM_ERR(" non vertex "<<in<<" of "<<nnodes_test)}
		// this order is not right see the addLHS2D for example 
		//const IdxTy dms[]={mla.nvar,mla.nvar,mla.di.n_var_dofs,mla.di.n_var_dofs};
		const IdxTy dms[]={nvars,nvdofs,nvars,nvdofs};
		//const IdxTy dmsr[]={mla.nvar,mla.di.n_var_dofs};
		const IdxTy dmsr[]={nvars,nvdofs};
		// make our own from the nodes... eliminate fe stuff wtf
		// FIXME TODO this is transposed, no idea why 
		MyIdxBlock my_dof_block(nvdofs,nvars);
		MyBlock rhs_block(MyBlock::loc(dmsr,2));
		MyBlock rhs_block_zed(MyBlock::loc(dmsr,2));
		MyBlock lhs_block(MyBlock::loc(dms,4));
		MyBlock lhs_block_zed(MyBlock::loc(dms,4));
		// mystery transpose, need to look at reform 
		mjm_libmesh_util::get_dof_mat(my_dof_block,elem, nvars);
		// this uses both now... 
 		//assemble_points( lhs_block_zed, rhs_block_zed, sol, elem, do_lhs); // ,false, 0);
		if (do_mjm_sparse)
 			assemble_points( lhs_block_zed, rhs_block_zed, my_sol_zed, elem, do_lhs,my_dof_block);  
 			else assemble_points( lhs_block_zed, rhs_block_zed, sol, elem, do_lhs,my_dof_block); // ,false, 0);


	//	if (rhs_block_zed(0,0)!=0) 
	//		{ MM_MSG(" rhs_block_zed(0,0)= "<<rhs_block_zed(0,0)) } 

		if (doing_j)
		{
			//const Real delta=.1; 
			//const Real delta=.001; 
			//const Real delta=.0001; 
			//const Real delta=.01; 
 			// const Real rdelta=1.0/delta;
			//const IdxTy n_var_dofs=(*m_mla).di.n_var_dofs;
			for (IdxTy var=0; var<mla.nvar ; ++var)
				//for (IdxTy j=0; j<n_var_dofs ; ++j)
				for (IdxTy j=0; j<nvdofs ; ++j)
				{// dof
				rhs_block.zero();
				lhs_block.zero();
				// doen't this change with fe parames?
                //const IdxTy dofj=(mla.di.dof_indices_var[var][j]);
                const IdxTy dofj=my_dof_block(var,j); // (mla.di.dof_indices_var[var][j]);
				Real xzed=(*solp)(dofj);
				Real my_xzed=my_sol(dofj);
				//Real deltaeff=xzed*.01;
				//Real deltaeff=xzed*.0001;
				//Real deltaeff=xzed*.00001;
				//Real deltaeff=xzed*.000001;
				//Real deltaeff=xzed*.0000001;
				// had been using this for most tests but the init conds are remaining
				//Real deltaeff=xzed*.000001;
				//Real deltaeff=xzed*.0000003;
				Real deltaeff=my_xzed*.00000013*1e-9;
				// this made worse, 
			//	Real deltaeff=xzed*.0001;
				//Real deltaeff=delta; // xzed*(0);
				if (deltaeff<0) deltaeff=-deltaeff;
				const Real deltamin=1e-28;
				//if (deltaeff<1e-25) deltaeff=1e-25;
				if (deltaeff<deltamin){
					const Real x1=lhs_block_zed(var,j,var,j);		
					const Real x2=rhs_block_zed(var,j);	
					if( 1e-16<(x1*x1)) deltaeff=x2/x1*1e-10;
					MM_ERR(" invoking delta min"<<x1<<" "<<x2<<" "<<deltaeff)
					if (deltaeff<0) deltaeff=-deltaeff;	
					if (deltaeff<deltamin) 	 deltaeff=deltamin;

				}
				//solp->add(dofj,delta);

				const Real oldv=(*solp)(dofj); // there was some reason not to 
				solp->add(dofj,deltaeff);
				my_sol(dofj)+=deltaeff;
				const Real fuddv=(*solp)(dofj); // there was some reason not to 

 				//assemble_points( lhs_block, rhs_block, *solp, elem, do_lhs); // ,false, 0);
 				if (do_mjm_sparse)
				{
					assemble_points( lhs_block, rhs_block, my_sol, elem
						, do_lhs, my_dof_block);  
       				//mla.di.jake(J,my_sol_zed,lhs_block_zed ,rhs_block_zed,&my_sol, 
       				mla.di.my_jake(my_sparse,my_sol_zed,lhs_block_zed ,rhs_block_zed,&my_sol, 
						lhs_block, rhs_block, 1.0/deltaeff,dofj  , my_dof_block);


				}
else
{	
				assemble_points( lhs_block, rhs_block, *solp, elem, do_lhs, my_dof_block); // ,false, 0);
//		if (rhs_block(0,0)!=0) { MM_MSG(" rhs_block(0,0)= "<<rhs_block(0,0)) } 

				// \delta(Ax-b)/delta is added to column mapped to dofj in jacobian
        		mla.di.jake(J,sol,lhs_block_zed
					//,rhs_block_zed,solp, lhs_block, rhs_block, rdelta,dofj);
					,rhs_block_zed,solp, lhs_block, rhs_block, 1.0/deltaeff,dofj  , my_dof_block);

}
		


				//solp->add(dofj,-delta);
				//solp->add(dofj,-deltaeff);
				solp->add(dofj,-fuddv);
				solp->add(dofj,oldv);
				my_sol(dofj)=my_xzed;
				// whty the FUDD does this not fudding work if add fudding workds
				//(*solp).set(dofj,oldv);
//				solp->el(dofj)=oldv;
// wrong state wtf	
//			solp->set(dofj,xzed);
				}
		} // doing_j
	//	if (rhs_block_zed(0,0)!=0) { MM_MSG(" rhs_block_zed(0,0)= "<<rhs_block_zed(0,0)) } 

//		if (doing_j)
		const bool let_di_reform=true;
		//lhs_block_zed.push_form(16,16);
		if (!let_di_reform) 
			{ lhs_block_zed.push_form(ntot,ntot); rhs_block_zed.push_form(ntot);}
//		MM_ERR(lhs_block_zed.diffsq(Ke,true))
//		lhs_block_zed.pop_form();
        //mla.di.residual(R,sol,Ke,Fe);
   
     	if (doing_residual&&do_mjm_sparse) mla.di.residual(R,my_sol_zed,lhs_block_zed,rhs_block_zed ,my_dof_block);
     	if (doing_residual&&!do_mjm_sparse) mla.di.residual(R,sol,lhs_block_zed,rhs_block_zed ,my_dof_block);
	//	if (rhs_block_zed(0,0)!=0) 
	//		{ MM_MSG(" rhs_block_zed(0,0)= "<<rhs_block_zed(0,0)) } 

//void assemble_jr (Tm* system, Tx* rhs, Ty & lhs_block, Tx & rhs_block)
		if (dump_rhs) mla.di.add_block(rhss_zed,rhs_block_zed,my_dof_block);
		if (dump_rhs) mla.di.add_block(rhss,rhs_block,my_dof_block);
		//if (dump_rhs) mla.di.add_block(rhss,rhs_block,rhs_block_zed,my_dof_block);


		if (doing_linear_ass)
			mla.di.assemble_jr (matrix, rhs, lhs_block_zed, rhs_block_zed ,my_dof_block);


		if (!let_di_reform) {rhs_block_zed.pop_form(); lhs_block_zed.pop_form();}

	} // el
		if (doing_j&&do_mjm_sparse) {my_sparse.to_libmesh(*J); }

		Real zres=0;
		ZeroEqn zjake;
		// always do this to get the locations to add for now 
		//zero_assembly( zres, doing_j?(&zjake):0, sol);
// this may be causing memory corrpution 
const bool do_zero_add=(!true)&&!doing_linear_ass;
if (do_zero_add) { MM_ONCE(" adding the zero eqns for current eaulity",) }
if(do_zero_add)		zero_assembly( zres, &zjake, sol);
		IdxTy zdof=~0;	
		//,const IdxTy _zdof, const Real scale=1, const IdxTy flags=0)
		//add_zero(J,R,zjake,zres,zdof);
		// max of 10 eqns each multipleid by 10 
		//add_zero(J,R,zjake,zres,10,1,1);

// so far the currents have been very close without this crap 
//MM_ERR(" zero add enabled but seems better without it ")
//if(do_zero_add)		add_zero(J,R,zjake,zres,zdof,1e-3,1);
// this seems to work better than 1e-3
//if(do_zero_add)		add_zero(J,R,zjake,zres,zdof,1e-5,1);
//if(do_zero_add)		add_zero(J,R,zjake,zres,zdof,1e-4,1);
if(do_zero_add)		add_zero(J,R,zjake,zres,zdof,1e-1,1);
else
		MM_ONCE(" zero add disable",)

if (dump_rhs)
{
MM_ONCE(" dumping rhs lots of output ",)
DumpRHS(rhss_zed,rhss); 

}
if (dump_jacobian)
{
MM_ONCE(" jacobian dumping enabled here ",) 
mjm_libmesh_util::DumpSparseMatrix(std::cerr,*J,"jacobian",1);

}


    (*m_cm).lap("ass","assbegin");
    (*m_cm).cum("asstotal","assbegin");
    (*m_cm).dump("EndAssemble", std::cerr); std::cerr<<std::endl;
    (*m_cm).inc("iterations");
    (*m_cm).clear("currentnew");
    (*m_cm).clear("areanew");

} // assembly

template <class Tv, class Ts,class Tfudd> void assembly_working( Tv * R, Ts * J, const Tfudd & sol)
{
	(*m_cm).mark("assbegin");
    static int ass_count=0;
    MM_MSG(" ass_count="<<ass_count)
    ++ass_count;
 	const bool doing_j=true&&(J); 
 	const bool do_lhs=true; 
	MLATy &   mla = * m_mla;
	// make own sparse and dense for access and hierarchy
    DMN  Ke; DVN Fe;
    DMN  Je; DVN Ge;
    DSM Ke_var[mla.nvar][mla.nvar] =
      {
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)}
      };
    DSV Fe_var[mla.nvar] = {DSV(Fe), DSV(Fe), DSV(Fe), DSV(Fe) };
    DSM Je_var[mla.nvar][mla.nvar] =
      {
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)}
      };
    DSV Ge_var[mla.nvar] = {DSV(Ge), DSV(Ge), DSV(Ge), DSV(Ge) };
    (*m_cm).set("nodes", mla.mesh.n_nodes());
	using mjm_libmesh_util::MakeXandGradXVector;
	auto solp=sol.clone();
 	//BoundaryInfo & bi = mla.mesh.get_boundary_info();
    MeshCItor      el     = mla.mesh.active_local_elements_begin();
    const MeshCItor end_el = mla.mesh.active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
        const Elem * elem = *el;
        //if (elem->has_children()) continue; 
        mla.re_init(elem);
//   const Real  xfuddshot=point_value(0,*(*elem).get_node(0),*elem);
//		MM_ERR(" fudd the dof ass fudd shot fuddkkkkk "<<xfuddshot)
        mjm_libmesh_util::ReInitMats(Ke_var, Fe_var,Ke,Fe,mla.di);
		const IdxTy dms[]={mla.nvar,mla.nvar,mla.di.n_var_dofs,mla.di.n_var_dofs};
		const IdxTy dmsr[]={mla.nvar,mla.di.n_var_dofs};
		MyBlock rhs_block(MyBlock::loc(dmsr,2));
		MyBlock rhs_block_zed(MyBlock::loc(dmsr,2));
		MyBlock lhs_block(MyBlock::loc(dms,4));
		MyBlock lhs_block_zed(MyBlock::loc(dms,4));
        if (doing_j) mjm_libmesh_util::ReInitMats(Je_var, Ge_var,Je,Ge,mla.di);
 		assemble_points(Ke_var, Fe_var, lhs_block_zed, rhs_block_zed, sol, elem, do_lhs,false, 0);
		if (doing_j)
		{
			const Real delta=.01; 
 			const Real rdelta=1.0/delta;
			const IdxTy n_var_dofs=(*m_mla).di.n_var_dofs;
			for (IdxTy var=0; var<mla.nvar ; ++var)
				for (IdxTy j=0; j<n_var_dofs ; ++j)
				{// dof
				rhs_block.zero();
				lhs_block.zero();
				// doen't this change with fe parames?
                const IdxTy dofj=(mla.di.dof_indices_var[var][j]);
				solp->add(dofj,delta);
 		assemble_points(Je_var, Ge_var, lhs_block, rhs_block, *solp, elem, do_lhs,false, 0);
				// \delta(Ax-b)/delta is added to column mapped to dofj in jacobian
        		mla.di.jake(J,sol,lhs_block_zed
					,rhs_block_zed,solp, lhs_block, rhs_block, rdelta,dofj);
				solp->add(dofj,-delta);
				}
		} // doing_j
		lhs_block_zed.push_form(16,16);
		rhs_block_zed.push_form(16);
		MM_ERR(lhs_block_zed.diffsq(Ke,true))
//		lhs_block_zed.pop_form();
        //mla.di.residual(R,sol,Ke,Fe);
        mla.di.residual(R,sol,lhs_block_zed,rhs_block_zed);
		rhs_block_zed.pop_form();
		lhs_block_zed.pop_form();

	} // el

    (*m_cm).lap("ass","assbegin");
    (*m_cm).cum("asstotal","assbegin");
    (*m_cm).dump("EndAssemble", std::cerr); std::cerr<<std::endl;
    (*m_cm).inc("iterations");
    (*m_cm).clear("currentnew");
    (*m_cm).clear("areanew");

} // assembly

//////////////////////////////////////////
//template <class Ts> void assembly( Ts & rhs, Ts & lhs, const Ts & sol)
template <class Tv, class Ts,class Tfudd> void assembly_old( Tv * R, Ts * J, const Tfudd & sol)
{
	(*m_cm).mark("assbegin");
    static int ass_count=0;
    MM_MSG(" ass_count="<<ass_count)
    ++ass_count;
 	const bool doing_j=true; 
	MLATy &   mla = * m_mla;
	// make own sparse and dense for access and hierarchy
    DMN  Ke; DVN Fe;
    DMN  Je; DVN Ge;
    DSM Ke_var[mla.nvar][mla.nvar] =
      {
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)},
        { DSM(Ke), DSM(Ke), DSM(Ke), DSM(Ke)}
      };
    DSV Fe_var[mla.nvar] = {DSV(Fe), DSV(Fe), DSV(Fe), DSV(Fe) };

    DSM Je_var[mla.nvar][mla.nvar] =
      {
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)},
        { DSM(Je), DSM(Je), DSM(Je), DSM(Je)}
      };
    DSV Ge_var[mla.nvar] = {DSV(Ge), DSV(Ge), DSV(Ge), DSV(Ge) };

    (*m_cm).set("nodes", mla.mesh.n_nodes());

	using mjm_libmesh_util::MakeXandGradXVector;


	auto solp=sol.clone();

 	BoundaryInfo & bi = mla.mesh.get_boundary_info();
    MeshCItor      el     = mla.mesh.active_local_elements_begin();
    const MeshCItor end_el = mla.mesh.active_local_elements_end();
    for ( ; el != end_el ; ++el)
    {
        const Elem * elem = *el;
        //if (elem->has_children()) continue; 
        mla.re_init(elem);
//   const Real  xfuddshot=point_value(0,*(*elem).get_node(0),*elem);
//		MM_ERR(" fudd the dof ass fudd shot fuddkkkkk "<<xfuddshot)
        mjm_libmesh_util::ReInitMats(Ke_var, Fe_var,Ke,Fe,mla.di);
		const IdxTy dms[]={mla.nvar,mla.nvar,mla.di.n_var_dofs,mla.di.n_var_dofs};
		const IdxTy dmsr[]={mla.nvar,mla.di.n_var_dofs};
		MyBlock rhs_block(MyBlock::loc(dmsr,2));
		MyBlock rhs_block_zed(MyBlock::loc(dmsr,2));
		MyBlock lhs_block(MyBlock::loc(dms,4));
		MyBlock lhs_block_zed(MyBlock::loc(dms,4));
        if (doing_j) mjm_libmesh_util::ReInitMats(Je_var, Ge_var,Je,Ge,mla.di);
        const IdxTy qpsz=mla.fev.qrule.n_points();
        for (IdxTy qp=0; qp<qpsz; qp++)
        {
            NumVec values; GradVec grads;
			// this requires a current solution
			// and these differ for fev and fes... doh
   			MakeXandGradXVector( values, grads,mla.di.dof_indices_var, mla.fev , qp, elem, mla.nvar,sol );
            mla.assemble_elem(Ke_var,Fe_var,values,grads,qp,bi,elem);
			// see if numerical ders will work, rhs should be most of problem. 
            if (doing_j)
			{	const bool do_lhs=!false;
				//const Real delta=.0001; 
				const Real delta=.0001; 
				const IdxTy szdof=mla.nvar*mla.di.n_var_dofs;
   				const IdxTy n_var_dofs = mla.di.n_var_dofs; 
				// di.dof_indices_var[0].size();
				for (IdxTy var=0; var<mla.nvar ; ++var)
				for (IdxTy j=0; j<n_var_dofs ; ++j)
				{// dof
                const IdxTy dofj=(mla.di.dof_indices_var[var][j]);
				solp->add(dofj,delta);
				NumVec valuesp; GradVec gradsp;
            	MakeXandGradXVector( valuesp, gradsp,mla.di.dof_indices_var, mla.fev , qp, elem, mla.nvar,*solp );
				mla.assemble_J(Je_var,Ge_var,valuesp,gradsp,qp,bi,elem,do_lhs);
        		mla.assemble_Jblock(lhs_block,rhs_block,valuesp,gradsp,qp,bi,elem,do_lhs);
				//sum_diff(Je,Ke,Ge,Fe,delta*.0001,j*4+var,szdof);
				//sum_diff(Je,Ke,Ge,Fe,-delta*1e-2,j*4+var,szdof);
				sum_diff(Je,Ke,Ge,Fe,-delta,j*4+var,szdof);
				// need to reform here 
				lhs_block.push_form(szdof,szdof);
				rhs_block.push_form(szdof);
				sum_diff(lhs_block,Ke,rhs_block,Fe,-delta,j*4+var,szdof);
				MM_ERR(lhs_block.diffsq(Je,!true))
				lhs_block.pop_form();
				rhs_block.pop_form();
// needs a liner one
//				MM_ERR(rhs_block.diffsq(Ge))
				//solp(dofj)-=delta;
				solp->add(dofj,-delta); //fudd this shot
				} //k 
			}
        //    if (doing_j) mla.assemble_Jblock(lhs_block,rhs_block,values,grads,qp,bi,elem);
        } //qp
        // this apssed a couple of tests so far
//        if (impose_g_diri) mla.di.impose(Ke,Fe,g_dirichlet);
		// this can be immediately summed 
        //mla.di.assembly(Ke,Fe,mla.system);
        //mla.di.sum(rhs,lhs,sol,Ke,Fe);
        mla.di.residual(R,sol,Ke,Fe);
        if (doing_j) mla.di.jacobian(J,sol,Je,Ge);




    } // for el 

    (*m_cm).lap("ass","assbegin");
    (*m_cm).cum("asstotal","assbegin");
    (*m_cm).dump("EndAssemble", std::cerr); std::cerr<<std::endl;
    (*m_cm).inc("iterations");
    (*m_cm).clear("currentnew");
    (*m_cm).clear("areanew");

} // assemble_bc

///////////////////////////////////
template <class Tm, class Tv,class Tmo, class Tvo>
				void sum_diff(Tm & Je,const Tmo & Ke,const Tv & Ge,const Tvo& Fe, const Real & delta, const IdxTy dofj, const IdxTy sz)
{
 const Real rdelta=1.0/delta;
for (IdxTy i=0; i<sz; ++i)
{// this should not need a dof map for i anyway? 
Je(i,dofj)+=(Ge(i)-Fe(i))*rdelta;
//const Real ji=(Ge(i)-Fe(i))*rdelta;
//Real jic=ji;
// this bombs but may just be bad guess from solver
//	CkTy::denan(jic," jacobian ",StrTy(__FILE__),(__LINE__));	

//MM_ERR(" jacobiran "<<i<<" "<<dofj<<" "<<ji)
} } //sum_diff
// this has to be static so the libmesh can call it 
// this only works if you call it ASSFUDD FUDD 
template <class ASSFUDD> 
static Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var )
{
// ASSFUDD 
//Number ASSFUDD=-111;
//return ASSFUDD; 
//MM_ERR(" assfudd init fudding called ")
    // we could just make a global instance 
    //return 0; 
// none of this fudding shot compiles fudd
//	return MLATy::init_value<Number,Parameters>(p,parameters,system,var,*cond());
    // we could just make a global instance 
 //return ASSFUDD::init_value<libMesh::Number,libMesh::Parameters>(p,parameters,system,var,ConditionsTy());
 return ASSFUDD::ASSFUDD(p,parameters,system,var,*cond());
    //return m_mla->init_value<Number,Parameters>(p,parameters,system,var,ConditionsTy());


//	return 
//		m_mla->init_value<Number,Parameters>(p,parameters,system,var); 
	// ,system,var,g_cond);
}

void dirichlet_map( DiMap & d) { m_dirichlet=&d; } 
void zero_map( ZeroEqn & d) { m_zero=&d; } 
void nevec( SurfVec & d) { m_surfaces=&d; } 
// no idea if these copy
Myt & configure(MLATy * pmla) { m_mla=pmla; return *this; } 
static void configure( const EGeo & _geo, const ConditionsTy & _cond)
{
	geo(&_geo);
	cond(&_cond);
}

static  CounterMap&  counter() { static CounterMap x;  return x;}
static const EGeo * geo(const EGeo * g=0) { static const EGeo * x=0; if (g!=0) x=g; return x;}
static const ConditionsTy * cond(const ConditionsTy * g=0) 
{ static const ConditionsTy * x=0; if (g!=0) x=g; return x;}
template <class Tx> static Tx * g(Tx * g=0) { static Tx * x=0; if (g!=0) x=g; return x;}

RWMap m_map; // ("ecmlanl.config");
DiMap * m_dirichlet;
ZeroEqn * m_zero;
SurfVec * m_surfaces;
EsTy & m_es;
MLATy * m_mla;

SpecialElements* m_transuranics ;


//EGeo m_geo;
//ConditionsTy m_cond;
CounterMap*  m_cm;
CounterMap m_cm_default;
IdxTy m_iteration;
IdxTy m_ass_flags;
MyBlock m_eqn_factors;

}; // mjm_numerical_jr;



#endif


