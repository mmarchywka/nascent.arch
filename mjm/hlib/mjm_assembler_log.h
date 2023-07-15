#ifndef MJM_LIBMESH_ASSEMBLER_LOG__
#define MJM_LIBMESH_ASSEMBLER_LOG__

// assume we do not need more than is included here 
#include "mjm_libmesh/mjm_assembler_foo.h"
#include "mjm_libmesh/mjm_instruments.h"
#include "mjm_libmesh/mjm_libmesh_boundaries.h"

class LogTag {};


template<class SysTy,class CTy, class ETy> class mjm_libmesh_assembler<SysTy,CTy,ETy,LogTag>
{
LogTag m_tag;
typedef mjm_libmesh_assembler<SysTy,CTy,ETy,LogTag> Myt;
typedef mjm_numb CkTy;
public:
enum { nvar=4};
enum { DECOUPLE=0}; // unFIXME
typedef SysTy system_type;
mjm_libmesh_assembler
    ( const StrTy & _system_name, EquationSystems & _es, const CTy & _cond, const ETy & _geo,
	 const bool add_diri=false)
:fuddinit(deb()),system_name(_system_name),es(_es), mesh(es.get_mesh())
  ,dim (mesh.mesh_dimension())
  ,system(es.get_system<SysTy>(system_name))
, u_var (system.variable_number ("u"))
, mplus_var(system.variable_number ("lmplus"))
, mzed_var(system.variable_number ("lmzed"))
, eminus_var(system.variable_number ("leminus"))
,di(system.get_dof_map(),nvar)
,cond(_cond)
,geo(_geo)
,checkpoint(wtf(add_diri))
,fev( di.dof_map,  dim, dim-0, libMesh::FIFTH,0)
,fes( di.dof_map,  dim, dim-1, libMesh::FIFTH,0)
//,cond(_cond)
//,geo(_geo)
,is_side(false)
,have_diri(add_diri)
,jfactor(1)
,jelec(0)
{
  libmesh_assert_equal_to (system_name, "Poisson");
//MM_MSG(" var no "<<u_var<<" "<<mplus_var<<" "<<mzed_var<<" "<<eminus_var)
modify("assembler_log_config");
// add_dirichlet();
}
// want to move to static and unspecial template def
void modify(const StrTy & nm)
{
ReadWriteMap rwm(nm);
rwm.get("jfactor",jfactor);
MM_ERR(" jfactor "<<jfactor);
}
IdxTy deb()
{
MM_ERR(" fudding ctor called ")
return 0; 
}
IdxTy wtf(const bool add_diri)
{
if (add_diri) add_dirichlet();
return 1; 
}
void add_dirichlet()
{
// eventually want to subset the paths, tracehere?
//const BoundaryInfo & bi=mesh.get_boundary_info();


// find the electrode and "infinity" conditions and apply them to dof maps
//DirichletBoundary db ;
//DirichletMaker(const IdxTy src, const D & value):m_src(src),m_value(value),m_func(m_value),m_db(0) {}
//const Db & boundary() const { return *m_db;}
//Db & boundary()  { return *m_db;}
//void add_var(const IdxTy v) {m_vars.push_back(v); }
DirichletMaker eplus(1, geo.voltage(0)) ;
DirichletMaker eminus(2, geo.voltage(1)) ;
eplus.add_var(0);
eminus.add_var(0);
eplus.init();
eminus.init();

//di.dof_map.add_dirichlet_boundary(eplus.boundary());
//di.dof_map.add_dirichlet_boundary(eminus.boundary());

//system.get_dof_map().add_dirichlet_boundary(eplus.boundary());

//if ((yf>.9)||(yf<-.9)) if ( xf>(-.5)) if (xf<.5)
//    mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.mzedeq),penalty,mzed_var);
//    mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.mpluseq),penalty,mplus_var);
//    mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.eeq),penalty,eminus_var);

// these now do entire periphery... 
DirichletMaker szed(0, fwd(cond.mzedeq)) ;
DirichletMaker splus(0, fwd(cond.mpluseq)) ;
DirichletMaker sminus(0, fwd(cond.eeq)) ;
szed.add_var(mzed_var);
splus.add_var(mplus_var);
sminus.add_var(eminus_var);
szed.init();
splus.init();
sminus.init();


//di.dof_map.add_dirichlet_boundary(szed.boundary());
//di.dof_map.add_dirichlet_boundary(splus.boundary());
//di.dof_map.add_dirichlet_boundary(sminus.boundary());


if ( !have_diri) { MM_ERR(" doing fudding init of es fudd ") 
es.init();
 }
else{ MM_ERR(" more complicaed  fudding es.init fudl ")
//di.dof_map.reinit(mesh);
//di.dof_map.create_dof_constraints(mesh);
//di.dof_map.process_constraints(mesh);
//di.dof_map.reinit(mesh);
//system.init();
es.init();
}
//
have_diri=true; 
}

  void re_init(const Elem * elem) { di.re_init(elem); fev.fe->reinit (elem); }

//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////

template <class Tt,class Tv, class Tx, class Ty >
        void assemble_elem(Tt & Ke_var, Tv &  Fe_var
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi, const Elem * elem)

   {
	 bool is_side=false; // kluge for noe
	IdxTy bsrc=mjm_libmesh_util::bsrc(bi,elem);
	is_side=mjm_libmesh_util::Val(bsrc);
	if (is_side) 
	{is_side=is_side;}
	else
	{is_side=is_side;}
	
	switch (bsrc)
	{
		case 0: { break; } // boundary
		case 1: { break; } // boundary
		case 2: { break; } // boundary
		case mjm_libmesh_util::BAD: { break; } // boundary
		default: { break; } // not known now. 

	}; // switch 

	(*this).is_side=is_side;
//	if (have_diri) if (is_side) return;  
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
   	//const Real coefs[]={-1.0,kT,kT,kT}; 
   	const Real bf=::exp(values[u_var]/cond.kT);
//  	CkTy::denan(bf," bf ",StrTy(__FILE__),(__LINE__));	
	//MM_ERR(" fudd u "<<values[u_var])
	// const Real d2u=((DECOUPLE==0)?1.0:0.0)*( - cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf)); //  unFIXME
	const Real d2u=((DECOUPLE==0)?1.0:0.0)*( - cond.q_over_e*(rev(values[mplus_var])-rev(values[eminus_var]))); //  unFIXME
 	CkTy::denan(d2u," d2u ",StrTy(__FILE__),(__LINE__));	

// not doing sides causes nan exposion 
//	if (!is_side) 
	mjm_libmesh_util::AddLHSD2( Ke_var,  qp,fev, di,coefs);

 const IdxTy qpsz=fes.qrule.n_points();

if (false)if (qp<qpsz) {	 // this generates oor from denorm etc. 
	for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
		if (elem->neighbor(kluge)==libmesh_nullptr)
		{	fes.fe->reinit(elem,kluge); 
			mjm_libmesh_util::AddLHSD2( Ke_var,  qp,fes, di,coefs);
if (DECOUPLE==0) { // unFIXME
	// change grad u signds 
// unFIXME
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fes,di,  qp, mplus_var, 1.0*( - grads[u_var]) ); // 
// finally the coefs on the zeroth derivative
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fes,di,  qp, eminus_var,  1.0*( grads[u_var]) ); //  
   	mjm_libmesh_util::AddRHS(Fe_var,fev,fes,di,(*this),qp,d2u,values,grads,true);

} 
// I thought this would work, it seemed to before????


		}
	fev.fe->reinit(elem);
}


if (DECOUPLE==0) { // unFIXME
	// change grad u signds 
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, mplus_var, 1.0*( - grads[u_var]) ); // 
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, mplus_var, cond.kT*( - grads[mplus_var]) ); // 
// finally the coefs on the zeroth derivative
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, eminus_var,  1.0*( grads[u_var]) ); //  
	mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, eminus_var,  cond.kT*( grads[eminus_var]) ); //  



   	mjm_libmesh_util::AddRHS(Fe_var,fev,fes,di,(*this),qp,d2u,values,grads,false);

} 
// I thought this would work, it seemed to before????
//        if ( internal_electrodes) ImposeVoltages( Ke_var,Fe_var, fev, di,cond,qp, elem);
	Real value;
	const Point p=fev.q_point[qp];
	jelec=0;
// this is distored terribly but we can apply DD now. 
//	if (geo.contains(value,is_side,p(0),p(1)))
	if (!have_diri) 
	if (bsrc>0) if (is_side) 
	{
		geo.voltage(value,bsrc-1,p(0),p(1));
		const Real penalty=1e10;
//		const Real valuei=0;
// with new mesh, handled by sides thing. 
	//MM_ERR(" lhs impose skipped "<<p(0)<<" "<<p(1)<<" "<<value) 
    if (false) 	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fev, di,qp, value,penalty,u_var);
// probably just impose zero normal slope, cut the grid and make them "sides" 
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,mzed_var);
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,mplus_var);
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,eminus_var);
		if (value>0) jelec=-1e28; else jelec=1e28;
		if (p(0)>0) jelec=jelec*cond.mueminus/cond.mumplus;
		//MM_ERR(" setting "<<p(0)<<" "<<p(1)<<" val "<<value)	
	}
}


////////////////////////////////////////////////////////////////////////




// we can call ourself for rhs independet of cond
    template<class Tx, class Ty, class Tp, class Tz>
    void rhs(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & values,
    const Tz & grads,const Tp &qpoint) const
{ // for now just keep punting 
    //cond.rhs(f_vec,di.nvar,d2u,values,grads,qpoint);
    //const Real gr_plus= (cond.kf*values[mzed_var]-cond.kr*values[mplus_var]*values[eminus_var]);
//    const Real gr_plus= (cond.kf*rev(values[mzed_var]-cond.kr*values[mplus_var]*values[eminus_var]);
 if (DECOUPLE!=0) return; // unFIXME 
// const Real gr_zed=-gr_plus;
 //   const Real gr_e=gr_plus;
// this needs to be fixed
  //  const Real jelectrode=geo.j(qpoint(0),qpoint(1));
    Real jelectrode=0; // (is_side)? geo.j(qpoint(0),qpoint(1)):0;
	if (jelec!=0) jelectrode=jelec;
    f_vec(u_var) = d2u;
	 Real pl=values[mplus_var]; //if (pl>200) pl=200;
	 Real mi=values[eminus_var]; //if (mi>200) mi=200;
	 Real ne=values[mzed_var]; //if (ne>200) ne=200;
	const Real rpl=rev(pl);
	const Real rmi=rev(mi); // going oor
	const Real rne=rev(ne);
//  	CkTy::denan(rpl," rpl ",StrTy(__FILE__),(__LINE__));	
//  	CkTy::denan(rmi," rmi ",StrTy(__FILE__),(__LINE__));	
//  	CkTy::denan(rne," rne ",StrTy(__FILE__),(__LINE__));	

    //f_vec(mzed_var) = -cond.kr+cond.kf*cond.umeinstein *(rev(pl+mi-ne));
    f_vec(mzed_var) = -cond.kr*cond.umeinstein+cond.kf*cond.umeinstein *(rpl*rmi/rne);
    //f_vec(mplus_var) = d2u+cond.ummplus*(cond.kr*rev(ne-pl) -cond.kf*rev(mi));
	// negative signs may work better here ... 
    f_vec(mplus_var) = d2u+cond.ummplus*(cond.kr*rne/rpl -cond.kf*rmi); // unFIXME sign
    //f_vec(eminus_var) = -d2u+cond.umeminus*(cond.kr*rev(ne-mi) -cond.kf*rev(pl)) 
    f_vec(eminus_var) =- d2u+cond.umeminus*(cond.kr*rne/rmi -cond.kf*rpl)  //  unFIXME sign 
			//+jelectrode*1e+6/rmi;
			+jelectrode*jfactor/rmi; // best so far 1.121e-17
//MM_ERR(" rhs "<<f_vec(0)<<" "<<f_vec(1)<<" "<<f_vec(2)<<" "<<f_vec(3) )
// this seems to be the right order
// MM_ERR(" rhs "<<u_var<<" "<<mzed_var<<" "<<mplus_var<<" "<<eminus_var )
}

template <class Tt,class Tv, class Te >
void side(Tt Ke_var[nvar][nvar],Tv Fe_var[],const Te * elem,const IdxTy side
, const BoundaryInfo & bi )
{
    fes.fe->reinit(elem, side);
    const bool internal_electrodes=!false;
	IdxTy bsrc=mjm_libmesh_util::bsrc(bi,elem);
// this is only a qp loop, may as well move it here... 
//  ImposeBCN(Ke_var,Fe_var,fes,di,cond,elem,side,internal_electrodes); 
    //mjm_libmesh_util::ImposeCallback(Ke_var,Fe_var,fes,di,(*this),elem,side,internal_electrodes);
    mjm_libmesh_util::ImposeCallback(Ke_var,Fe_var,fes,di,(*this),elem,bsrc,internal_electrodes);


}

template <class Tt,class Tv, class Te >
void side(Tt Ke_var[nvar][nvar],Tv Fe_var[],const Te * elem,const IdxTy side)
{
    fes.fe->reinit(elem, side);
    const bool internal_electrodes=!false;
// this is only a qp loop, may as well move it here... 
//  ImposeBCN(Ke_var,Fe_var,fes,di,cond,elem,side,internal_electrodes); 
    mjm_libmesh_util::ImposeCallback(Ke_var,Fe_var,fes,di,(*this),elem,side,internal_electrodes);


}

template<class Tx, class Ty, class Tp>
    void sides(Tx & Ke_var, Ty & Fe_var,const Tp * elem,const IdxTy side
,const bool skip_electrodes,const Point & p,const Real & penalty,const IdxTy qp)
{
        const Real xf = p(0);
        const Real yf = p(1);
const fe_info & fe=fes;
// eh with the new mesh the interior electrode surfaces are now
// sides as far as libmesh is concerned... doh 
//if (!false) if (!skip_electrodes)
{
Real value;
//	if (geo.epoint(value,true,xf,yf,0))
	//if (geo.epoint2(value,true,p(0),p(1),0))
//	if (geo.contains(value,true,p(0),p(1)))
if (!have_diri) if (  geo.voltage(value,side-1,p(0),p(1)))
{
//if (yf*yf<.7*.7) if (xf*xf<.7*.7)
		//MM_ERR(" actually imposing "<<p(0)<<" "<<p(1)<<" "<<value)
    	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, value,penalty,u_var);
		}
}

/*
const Real jelectrode=geo.j(p(0),p(1));
if (jelectrode!=0)
{
	const Real mi=values[eminus_var];
	const Real ne=values[mzed_var];
	const Real rpl=rev(pl);
	const Real rmi=rev(mi);
}
*/

// also need a y gate for internal electrodes
// only electrodes now 
 if (!have_diri) 
if ((yf>.9)||(yf<-.9))
if ( xf>(-.5)) if (xf<.5)
{
	//MM_ERR(" imposing bc "<<xf<<" "<<yf) 
	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.mzedeq),penalty,mzed_var);
 	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.mpluseq),penalty,mplus_var);
 	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fe, di,qp, fwd(cond.eeq),penalty,eminus_var);
//MM_ERR(" imposing on other planes "<<fwd(cond.mzedeq))
}

}

template<class Ty, class Tx> Real j(const Ty & values, const Tx & grads, const IdxTy & dir) const
{
	// note that the carrier grad is log.... 
	Real pl=values[mplus_var]; //if (pl>200) pl=200;
	Real mi=values[eminus_var]; //if (mi>200) mi=200;
//	Real ne=values[mzed_var]; //if (ne>200) ne=200;

// unused 
//	const Real rpl=rev(pl);
//	const Real rmi=rev(mi); // going oor


//	const Real rne=rev(ne);
	Real jv=0;
	jv+= cond.q*cond.mumplus*pl*(grads[u_var](dir)-cond.kT*grads[mplus_var](dir));
	jv+= cond.q*cond.mueminus*mi*(grads[u_var](dir)+cond.kT*grads[eminus_var](dir));

return jv;
}

//template <typename ES, class Ty > static void init_cd (ES & es,
template < class Ty > static void init_cd (EquationSystems & es,
              const Ty & system_name)
{
//  libmesh_assert_equal_to (system_name, "Poisson");

// unused but ok 
//  auto & system = es.get_system<Myt::system_type>(system_name);


//    es.get_system<LinearImplicitSystem>("Poisson");
  // Project initial conditions at time 0
//  es.parameters.set<Real> ("time") = system.time = 0;
//  system.project_solution(exact_value, libmesh_nullptr, es.parameters);
// max_iter
  //es.parameters.set<unsigned int>("linear solver maximum iterations") =3000;
  es.parameters.set<IdxTy>("linear solver maximum iterations") =1000;

}


template< class Number,class Parameters > static Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var ,
					const CTy & cond )
{
if( var=="lmplus") return fwd(cond.mpluseq);
if( var=="lmzed") return fwd(cond.mzedeq);
if( var=="leminus") return fwd(cond.eeq);
// these apparently need to be changed with electrode configuration
if( var=="u")
{
//if (p(0)>.99) if ( p(1)*p(1)<.4) return default_geo().vplus();
//if (p(0)<-.99) if ( p(1)*p(1)<.4) return default_geo().vminus();
///Real value;
//if (default_geo().epoint2(value,false,p(0),p(1),0)) return value; 


//if (default_geo().contains(value,!false,p(0),p(1),0)) return value; 
Real v= default_geo().guess_init_value(p(0),p(1));

return v;
}
MM_MSG(" unknown var "<<var)

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

//	const bool have_diri=true;
//	MLATY  mla(system_name,es,g_cond,g_geo);
//	if (have_diri) mla.add_dirichlet();
// 	add_dirichlet();
return system;
}

static Number fwd(const Number & x) { return ::log(x); } 
static Number rev(const Number &  x) {Number y=x; //if ( y>200){ y=200; MM_ERR(" x = "<<x) }
 return ::exp(y); } 
// this does not work as the rhs was divied by exp 
//static Number fwd(const Number & x) { return (x); } 
//static Number rev(const Number & x) { return (x); } 



//private:
	const IdxTy fuddinit;
    StrTy system_name;
    EquationSystems & es;
	// now using this for refinement fixes non-const
   //  const 
	MeshBase & mesh;
	//ZZSerialMesh & mesh;
    const IdxTy dim;
    //  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");
    SysTy & system;
    const VarTy  u_var; //  = system.variable_number ("u");
    const VarTy  mplus_var; //  =g_using_sd?(~0): system.variable_number ("mplus");
    const VarTy  mzed_var ; // = system.variable_number ("mzed");
    const VarTy  eminus_var ; // = g_using_sd?(~0):system.variable_number ("eminus");


    dof_info di ; //=dof_info(system.get_dof_map(),nvar); 
    const CTy & cond;
    const ETy & geo;
	
	const IdxTy checkpoint;

    fe_info fev ;// = fe_info( di.dof_map,  dim, dim-0, FIFTH,0);
    fe_info fes ;// = fe_info( di.dof_map,  dim, dim-1, FIFTH,0);
   // const CTy & cond;
   // const ETy & geo;
	
	bool is_side;
	bool have_diri;
	Real jfactor,jelec;	
	// this will fail to get the snapped geometries ???????
	static const ETy & default_geo() { static ETy x;  return x; } 

	
}; // mjm_libmesh_assembler


#if 0
//template <class Tt,int NVAR,class Tv, class Tx, class Ty > 
//        void assemble(Tt Ke_var[NVAR][NVAR], Tv Fe_var[]    
//  ,const Tx & values,const Ty & grads, const IdxTy qp)
//template <class Tt,class Tv, class Tx, class Ty > 
//        void assemble(Tt Ke_var[nvar][nvar], Tv Fe_var[nvar]    
//  ,const Tx & values,const Ty & grads, const IdxTy qp)

template <class Tt,class Tv, class Tx, class Ty >
        void assemble(Tt & Ke_var, Tv &  Fe_var
    ,const Tx & values,const Ty & grads, const IdxTy qp, const bool is_side=false)

   {
	(*this).is_side=is_side;
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
   	const Real bf=::exp(values[u_var]/cond.kT);
//  	CkTy::denan(bf," bf ",StrTy(__FILE__),(__LINE__));	
	//const Real d2u=(bf<1e-200)?0:(  -cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf));
	const Real d2u=(  -cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf));
 	CkTy::denan(d2u," d2u ",StrTy(__FILE__),(__LINE__));	
// hardcoded
//   AddLHS(  Ke_var, fev,cond,  qp ,di.n_var_dofs,coefs,d2u,values,grads);

// add the laplacian
//AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs,coefs);
mjm_libmesh_util::AddLHSD2( Ke_var,  qp,fev, di,coefs);
// AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs(0),coefs);
 
//AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs(0),coefs[0],0);
// AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs(1),coefs[1],1);
// AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs(2),coefs[2],2);
// AddLHSD2( Ke_var,  fev,  qp ,di.n_var_dofs(3),coefs[3],3);

// then the lower order 
//AddLHSD1(  Ke_var,  fev,  qp, mplus_var, di.n_var_dofs(mplus_var),  -grads[u_var] );
//AddLHSD1(  Ke_var,  fev,  qp, mplus_var, di.n_var_dofs,  -grads[u_var] );
mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, mplus_var,   -grads[u_var] );
// finally the coefs on the zeroth derivative
//AddLHSD0( Ke_var,  fev,  qp,  cond.mplusv ,di.n_var_dofs,  d2u );
//AddLHSD1(  Ke_var,  fev,  qp, eminus_var, di.n_var_dofs(eminus_var),  grads[u_var] );
//AddLHSD1(  Ke_var,  fev,  qp, eminus_var, di.n_var_dofs,  grads[u_var] );
mjm_libmesh_util::AddLHSD1(  Ke_var,  fev,di,  qp, eminus_var,   grads[u_var] );
//AddLHSD0( Ke_var,  fev,  qp,  cond.eminusv ,di.n_var_dofs,  d2u );


   // This has included some LHS terms
//   this punts to cond
//   AddForcingFunctionsCond(Fe_var,fev,fes,di,cond,qp,d2u,values,grads);
// we can get called back and do it here 
   mjm_libmesh_util::AddRHS(Fe_var,fev,fes,di,(*this),qp,d2u,values,grads);
// I thought this would work, it seemed to before????
//        if ( internal_electrodes) ImposeVoltages( Ke_var,Fe_var, fev, di,cond,qp, elem);
	Real value;
	const Point p=fev.q_point[qp];
	jelec=0;
// this is distored terribly but we can apply DD now. 
// this is now for he sides handler... 
// includes boarder
	//if (geo.epoint2(value,is_side,p(0),p(1),0))
	if (geo.contains(value,is_side,p(0),p(1)))
	{
		const Real penalty=1e10;
//		const Real valuei=0;
// with new mesh, handled by sides thing. 
	//MM_ERR(" lhs impose skipped "<<p(0)<<" "<<p(1)<<" "<<value) 
    if (false) 	mjm_libmesh_util::ImposeN( Ke_var, Fe_var, fev, di,qp, value,penalty,u_var);
// probably just impose zero normal slope, cut the grid and make them "sides" 
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,mzed_var);
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,mplus_var);
//    	ImposeN( Ke_var, Fe_var, fev, di,qp, valuei,penalty,eminus_var);
		if (value>0) jelec=-1e28; else jelec=1e28;
		if (p(0)>0) jelec=jelec*cond.mueminus/cond.mumplus;
		//MM_ERR(" setting "<<p(0)<<" "<<p(1)<<" val "<<value)	
	}
}




#endif // 0 


#endif // guard

