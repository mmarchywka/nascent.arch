#ifndef MJM_LIBMESH_ASSEMBLER_ZED__
#define MJM_LIBMESH_ASSEMBLER_ZED__

// assume we do not need more than is included here 
#include "mjm_libmesh/mjm_assembler_foo.h"
#include "mjm_libmesh/mjm_instruments.h"

class ZedTag {};


template<class SysTy,class CTy, class ETy> class mjm_libmesh_assembler<SysTy,CTy,ETy,ZedTag>
{
ZedTag m_tag;
typedef  mjm_libmesh_assembler<SysTy,CTy,ETy,ZedTag> Myt;
public:
enum { nvar=4};
typedef SysTy system_type;
mjm_libmesh_assembler
    ( const StrTy & _system_name, EquationSystems & _es, const CTy & _cond, const ETy & _geo)
:system_name(_system_name),es(_es), mesh(es.get_mesh())
  ,dim (mesh.mesh_dimension())
  ,system(es.get_system<SysTy>(system_name))
, u_var (system.variable_number ("u"))
, mplus_var(system.variable_number ("lmplus"))
, mzed_var(system.variable_number ("lmzed"))
, eminus_var(system.variable_number ("leminus"))
,di(system.get_dof_map(),nvar)
,fev( di.dof_map,  dim, dim-0, libMesh::FIFTH,0)
,fes( di.dof_map,  dim, dim-1, libMesh::FIFTH,0)
,cond(_cond)
,geo(_geo)
,is_side(false)
,jfactor(0)
,jelec(0)
{
  libmesh_assert_equal_to (system_name, "Poisson");
//MM_MSG(" var no "<<u_var<<" "<<mplus_var<<" "<<mzed_var<<" "<<eminus_var)
//modify("assembler_log_config");
}
// want to move to static and unspecial template def

  void re_init(const Elem * elem) { di.re_init(elem); fev.fe->reinit (elem); }

template <class Tt,class Tv, class Tx, class Ty >
        void assemble(Tt & Ke_var, Tv &  Fe_var
    ,const Tx & values,const Ty & grads, const IdxTy qp, const bool is_side=false)

   {
	(*this).is_side=is_side;
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1,1,1,1};
   	const Real bf=0; // ::exp(values[u_var]/cond.kT);
  	const Real d2u= 0; // -cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf);
// add the laplacian
 ForceZed( Ke_var,  fev,  qp ,di.n_var_dofs,coefs);
 AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs,coefs);

// then the lower order 
//AddLHSD1(  Ke_var,  fev,  qp, mplus_var, di.n_var_dofs,  -grads[u_var] );
// finally the coefs on the zeroth derivative
//AddLHSD0( Ke_var,  fev,  qp,  cond.mplusv ,di.n_var_dofs,  d2u );
//AddLHSD1(  Ke_var,  fev,  qp, eminus_var, di.n_var_dofs,  grads[u_var] );
//AddLHSD0( Ke_var,  fev,  qp,  cond.eminusv ,di.n_var_dofs,  d2u );
   AddRHS(Fe_var,fev,fes,di,(*this),qp,d2u,values,grads);
// I thought this would work, it seemed to before????
//        if ( internal_electrodes) ImposeVoltages( Ke_var,Fe_var, fev, di,cond,qp, elem);
	jelec=0;
// this is distored terribly but we can apply DD now. 
}

// we can call ourself for rhs independet of cond
    template<class Tx, class Ty, class Tp, class Tz>
    void rhs(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & values,
    const Tz & grads,const Tp &qpoint) const
{ 

    Real jelectrode=0; 
	if (jelec!=0) jelectrode=jelec;
    f_vec(u_var) = 0;
    f_vec(mzed_var) = 0; 
    f_vec(mplus_var) = 0;
    f_vec(eminus_var) = 0;
// MM_ERR(" rhs "<<u_var<<" "<<mzed_var<<" "<<mplus_var<<" "<<eminus_var )
}


template <class Tt,class Tv, class Te >
void side(Tt Ke_var[nvar][nvar],Tv Fe_var[],const Te * elem,const IdxTy side)
{
    fes.fe->reinit(elem, side);
    const bool internal_electrodes=!false;
// this is only a qp loop, may as well move it here... 
//  ImposeBCN(Ke_var,Fe_var,fes,di,cond,elem,side,internal_electrodes); 
    ImposeCallback(Ke_var,Fe_var,fes,di,(*this),elem,side,internal_electrodes);


}

template<class Tx, class Ty, class Tp>
    void sides(Tx & Ke_var, Ty & Fe_var,const Tp * elem,const IdxTy side
,const bool skip_electrodes,const Point & p,const Real & penalty,const IdxTy qp)
{
// do nothing here  default to slope is zero 
}

//template <typename ES, class Ty > static void init_cd (ES & es,
template < class Ty > static void init_cd (EquationSystems & es,
              const Ty & system_name)
{
//  libmesh_assert_equal_to (system_name, "Poisson");
  auto & system = es.get_system<Myt::system_type>(system_name);
//    es.get_system<LinearImplicitSystem>("Poisson");
  // Project initial conditions at time 0
//  es.parameters.set<Real> ("time") = system.time = 0;
//  system.project_solution(exact_value, libmesh_nullptr, es.parameters);
// max_iter
  //es.parameters.set<unsigned int>("linear solver maximum iterations") =3000;
  es.parameters.set<IdxTy>("linear solver maximum iterations") =3000;

}


template< class Number,class Parameters > static Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var ,
					const CTy & cond )
{
return 0;


}

static system_type &  setupES(EquationSystems & equation_systems,const StrTy & name="Poisson")
{
  equation_systems.add_system<system_type> (name);
  system_type & system=  equation_systems.get_system<system_type>(name);
  system.add_variable("u", libMesh::FIRST);
  system.add_variable("lmzed", libMesh::FIRST);
  system.add_variable("lmplus", libMesh::FIRST);
  system.add_variable("leminus", libMesh::FIRST);

return system;
}

//private:
    StrTy system_name;
    EquationSystems & es;
    const MeshBase & mesh;
    const IdxTy dim;
    //  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");
    SysTy & system;
    const VarTy  u_var; //  = system.variable_number ("u");
    const VarTy  mplus_var; //  =g_using_sd?(~0): system.variable_number ("mplus");
    const VarTy  mzed_var ; // = system.variable_number ("mzed");
    const VarTy  eminus_var ; // = g_using_sd?(~0):system.variable_number ("eminus");


    dof_info di ; //=dof_info(system.get_dof_map(),nvar); 
    fe_info fev ;// = fe_info( di.dof_map,  dim, dim-0, FIFTH,0);
    fe_info fes ;// = fe_info( di.dof_map,  dim, dim-1, FIFTH,0);
    const CTy & cond;
    const ETy & geo;
	
	bool is_side;
	Real jfactor,jelec;	
	static const ETy & default_geo() { static ETy x;  return x; } 

	
}; // mjm_libmesh_assembler




#endif

