#ifndef MJM_LIBMESH_ASSEMBLER_FOO__
#define MJM_LIBMESH_ASSEMBLER_FOO__


#include <mjm_globals.h>
#include <mjm_species.h>
#include <mjm_geometry.h>
#include <mjm_libmesh.h>
#include <mjm_instruments.h>
// C++ include files that we need
#include <iostream>
#include <string>
#include <vector>
#include <math.h>


class FooTag {};
//class LogTag {};


// make a non templated base class?
template<class SysTy,class CTy, class Ety, class TagTy> class mjm_libmesh_assembler
{
protected:
typedef ReadWriteMap ConfigMap;
/*void modify(const StrTy & nm) { ReadWriteMap rwm(nm); //rwm.get("jfactor",jfactor); //MM_ERR(" jfactor "<<jfactor); } */
static ConfigMap & config() { static ConfigMap rwm; return rwm; } 

};

// this needs a few prototypes althogh those will move into the headers above
// the specialization seem to be allowed very late

template<class SysTy,class CTy, class ETy> class mjm_libmesh_assembler<SysTy,CTy,ETy,FooTag>
{
FooTag m_tag;
typedef  mjm_libmesh_assembler<SysTy,CTy,ETy,FooTag> Myt;
public:
enum { nvar=4};
typedef SysTy system_type;
mjm_libmesh_assembler
    ( const StrTy & _system_name, EquationSystems & _es, const CTy & _cond, const ETy & _geo)
:system_name(_system_name),es(_es), mesh(es.get_mesh())
  ,dim (mesh.mesh_dimension())
  ,system(es.get_system<SysTy>(system_name))
 // ,system(es.get_system<SysTy>("Poisson"))
, u_var (system.variable_number ("u"))
, mplus_var(system.variable_number ("mplus"))
, mzed_var(system.variable_number ("mzed"))
, eminus_var(system.variable_number ("eminus"))
//,di(dof_info(system.get_dof_map(),nvar)) 
//,fev(fe_info( di.dof_map,  dim, dim-0, FIFTH,0))
//,fes(fe_info( di.dof_map,  dim, dim-1, FIFTH,0))
,di(system.get_dof_map(),nvar)
,fev( di.dof_map,  dim, dim-0, libMesh::FIFTH,0)
,fes( di.dof_map,  dim, dim-1, libMesh::FIFTH,0)
,cond(_cond)
,geo(_geo)
{
  libmesh_assert_equal_to (system_name, "Poisson");
//MM_MSG(" var no "<<u_var<<" "<<mplus_var<<" "<<mzed_var<<" "<<eminus_var)

}

  void re_init(const Elem * elem) { di.re_init(elem); fev.fe->reinit (elem); }

//template <class Tt,int NVAR,class Tv, class Tx, class Ty > 
//        void assemble(Tt Ke_var[NVAR][NVAR], Tv Fe_var[]    
//  ,const Tx & values,const Ty & grads, const IdxTy qp)
//template <class Tt,class Tv, class Tx, class Ty > 
//        void assemble(Tt Ke_var[nvar][nvar], Tv Fe_var[nvar]    
//  ,const Tx & values,const Ty & grads, const IdxTy qp)

template <class Tt,class Tv, class Tx, class Ty >
        void assemble(Tt & Ke_var, Tv &  Fe_var
    ,const Tx & values,const Ty & grads, const IdxTy qp)

   {const Real kT=cond.kT; // .0259;
   const Real coefs[]={1.0,kT,kT,kT};
   const Real bf=::exp(values[cond.uv]/cond.kT);
   const Real d2u=  -cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf);
// hardcoded
//   AddLHS(  Ke_var, fev,cond,  qp ,di.n_var_dofs,coefs,d2u,values,grads);

// add the laplacian
 AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs,coefs);

// then the lower order 
AddLHSD1(  Ke_var,  fev,  qp, cond.mplusv, di.n_var_dofs,  grads[cond.uv] );
// finally the coefs on the zeroth derivative
AddLHSD0( Ke_var,  fev,  qp,  cond.mplusv ,di.n_var_dofs,  d2u );
AddLHSD1(  Ke_var,  fev,  qp, cond.eminusv, di.n_var_dofs,  grads[cond.uv] );
AddLHSD0( Ke_var,  fev,  qp,  cond.eminusv ,di.n_var_dofs,  d2u );


   // This has included some LHS terms
//   this punts to cond
//   AddForcingFunctionsCond(Fe_var,fev,fes,di,cond,qp,d2u,values,grads);
// we can get called back and do it here 
   AddRHS(Fe_var,fev,fes,di,(*this),qp,d2u,values,grads);
// I thought this would work, it seemed to before????
//        if ( internal_electrodes) ImposeVoltages( Ke_var,Fe_var, fev, di,cond,qp, elem);
}

// we can call ourself for rhs independet of cond
    template<class Tx, class Ty, class Tp, class Tz>
    void rhs(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & values,
    const Tz & grads,const Tp &qpoint) const
{ // for now just keep punting 
    //cond.rhs(f_vec,di.nvar,d2u,values,grads,qpoint);
    const Real gr_plus=
        (cond.kf*values[cond.mzedv]-cond.kr*values[cond.mplusv]*values[cond.eminusv]);
    const Real gr_zed=-gr_plus;
    const Real gr_e=gr_plus;
    const Real jelectrode=geo.j(qpoint(0),qpoint(1));
    f_vec(cond.uv) = d2u;
    f_vec(cond.mzedv) =  cond.umeinstein*gr_zed;
    f_vec(cond.mplusv) = cond.ummplus*gr_plus;
    f_vec(cond.eminusv) = cond.umeminus*gr_e+jelectrode;


}


template <class Tt,class Tv, class Te >
void side(Tt Ke_var[nvar][nvar],Tv Fe_var[],const Te * elem,const IdxTy side)
{
    fes.fe->reinit(elem, side);
    const bool internal_electrodes=false;
// this is only a qp loop, may as well move it here... 
//  ImposeBCN(Ke_var,Fe_var,fes,di,cond,elem,side,internal_electrodes); 
    ImposeCallback(Ke_var,Fe_var,fes,di,(*this),elem,side,internal_electrodes);


}

template<class Tx, class Ty, class Tp>
    void sides(Tx & Ke_var, Ty & Fe_var,const Tp * elem,const IdxTy side
,const bool skip_electrodes,const Point & p,const Real & penalty,const IdxTy qp)
{
        const Real xf = p(0);
        const Real yf = p(1);
const fe_info & fe=fes;
if (!skip_electrodes)
{
Real value;
if (geo.epoint(value,true,xf,yf,0))
    ImposeN( Ke_var, Fe_var, fe, di,qp, value,penalty,cond.uv);
}

if ( xf>(-.5)) if (xf<.5)
{
 ImposeN( Ke_var, Fe_var, fe, di,qp, cond.mpluseq,penalty,cond.mplusv);
 ImposeN( Ke_var, Fe_var, fe, di,qp, cond.mzedeq,penalty,cond.mzedv);
 ImposeN( Ke_var, Fe_var, fe, di,qp, cond.eeq,penalty,cond.eminusv);
}

}

//template <typename ES, class Ty > static void init_cd (ES & es,
template < class Ty > static void init_cd (EquationSystems & es,
              const Ty & system_name)
{
  libmesh_assert_equal_to (system_name, "Poisson");
  auto & system = es.get_system<Myt::system_type>("Poisson");
//    es.get_system<LinearImplicitSystem>("Poisson");
  // Project initial conditions at time 0
//  es.parameters.set<Real> ("time") = system.time = 0;
//  system.project_solution(exact_value, libmesh_nullptr, es.parameters);
// max_iter
  //es.parameters.set<unsigned int>("linear solver maximum iterations") =3000;
  es.parameters.set<IdxTy>("linear solver maximum iterations") =300;

}


template< class Number,class Parameters > static Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var ,
					const CTy & cond )
{
if( var=="mplus") return cond.mpluseq;
if( var=="mzed") return cond.mzedeq;
if( var=="eminus") return cond.eeq;
if( var=="u")
{
if (p(0)>.99) if ( p(1)*p(1)<.4) return default_geo().vplus();
if (p(0)<-.99) if ( p(1)*p(1)<.4) return default_geo().vminus();
return 0;
}
MM_MSG(" unknown var "<<var)

return 0;


}

static system_type &  setupES(EquationSystems & equation_systems)
{
  equation_systems.add_system<system_type> ("Poisson");
  system_type & system=  equation_systems.get_system<system_type>("Poisson");
  system.add_variable("u", libMesh::FIRST);
  system.add_variable("mzed", libMesh::FIRST);
  system.add_variable("mplus", libMesh::FIRST);
  system.add_variable("eminus", libMesh::FIRST);

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

static const ETy & default_geo() { static ETy x;  return x; } 
}; // mjm_libmesh_assembler


#endif

