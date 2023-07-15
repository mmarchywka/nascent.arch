#ifndef MJM_LIBMESH_ASSEMBLER_NLOG__
#define MJM_LIBMESH_ASSEMBLER_NLOG__

// assume we do not need more than is included here 
#include "mjm_libmesh/mjm_assembler_foo.h"
#include "mjm_libmesh/mjm_instruments.h"
#include "mjm_libmesh/mjm_libmesh_boundaries.h"

#include "mjm_block_matrix.h"

// derived from LogTag try to get all the terms right apparently
// including nl term for now just iterate to fixed point, 
// calculate jacobians etc later depending on how this goes....

class NLogTag {};


template<class SysTy,class CTy, class ETy> class mjm_libmesh_assembler<SysTy,CTy,ETy,NLogTag>
{
NLogTag m_tag;
typedef mjm_libmesh_assembler<SysTy,CTy,ETy,NLogTag> Myt;
typedef mjm_numb CkTy;
typedef mjm_block_matrix<Real> MyBlock;
public:
enum { nvar=4};
enum { DECOUPLE=0}; // unFIXME
typedef SysTy system_type;
typedef SysTy nl_system_type;
/*
mjm_libmesh_assembler
    ( const StrTy & _system_name, EquationSystems & _es, const CTy & _cond, const ETy & _geo,
	 const bool add_diri=false)
:fuddinit(deb()),system_name(_system_name),es(_es), mesh(es.get_mesh())
  ,dim (mesh.mesh_dimension())
  ,ass_fudd_system(es.get_system<SysTy>(system_name))
, u_var (ass_fudd_system.variable_number ("u"))
, mplus_var(ass_fudd_system.variable_number ("lmplus"))
, mzed_var(ass_fudd_system.variable_number ("lmzed"))
, eminus_var(ass_fudd_system.variable_number ("leminus"))
,di(ass_fudd_system.get_dof_map(),nvar)
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
//  libmesh_assert_equal_to (system_name, "Poisson");
//MM_MSG(" var no "<<u_var<<" "<<mplus_var<<" "<<mzed_var<<" "<<eminus_var)
modify("assembler_log_config");
// add_dirichlet();
}
*/
template <class Fu> mjm_libmesh_assembler
    ( Fu & system, EquationSystems & _es, const CTy & _cond, const ETy & _geo,int ass,
	 const bool add_diri=false)
:fuddinit(deb()),system_name("ASSFUDD"),es(_es), mesh(es.get_mesh())
  ,dim (mesh.mesh_dimension())
  ,ass_fudd_system(es.get_system<SysTy>(system_name))
  //,system((SysTy)&sys)
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
// wth !true, residuals failed to decrease and u was never updated however
// the species distro was more interestng 
,d2u_on_right(true)
{
//  libmesh_assert_equal_to (system_name, "Poisson");
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
 MM_ERR(" Diri crap removed from nlog ass ") 
es.init();
//
have_diri=true; 
}

  void re_init(const Elem * elem) { di.re_init(elem); fev.fe->reinit (elem); }


////////////////////////////////////////////////////////////////////////////
void debug_branches(bool is_side, const IdxTy bsrc)
{
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

} // debug_branches

template <class Tt,class Tv, class Tx, class Ty,class Tf >
        void assemble_fex(Tt & Ke_var, Tv & Fe_var, Tf & fe
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi
, const Real & d2u, const Real * coefs, const IdxTy bsrc, const bool is_side, const bool do_lhs=true)
{
using namespace mjm_libmesh_util;
	//mjm_libmesh_util::
if (do_lhs)		AddLHSD2( Ke_var, qp,fe, di,coefs);

if (DECOUPLE==0) { // unFIXME
	const Real kT=cond.kT; // .0259;
	//mjm_libmesh_util::
//		AddLHSD1( Ke_var, fe,di, qp, mplus_var, 1.0*(  grads[u_var]) ); // 
if (do_lhs) {
AddLHSD1( Ke_var, fe,di, qp, mplus_var, kT*(  grads[mplus_var]) ); // 
AddLHSD1( Ke_var, fe,di, qp, mplus_var, ( -grads[u_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, eminus_var,kT*( grads[eminus_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, eminus_var, ( grads[u_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, mzed_var,kT*( grads[mzed_var]) ); //  
}
AddRHS(Fe_var,fe,fe,di,(*this),qp,d2u,values,grads,is_side);
} // decouple

} // assemble_rhs

template <class Tt,class Tv, class Tx, class Ty,class Tf >
        void assemble_fex_block(Tt & Ke_block, Tv & Fe_block, Tf & fe
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi
, const Real & d2u, const Real * coefs, const IdxTy bsrc, const bool is_side, const bool do_lhs=true)
{
using namespace mjm_libmesh_util;
	//mjm_libmesh_util::
if (do_lhs)		AddLHSD2Block( Ke_block, qp,fe, di,coefs);

if (DECOUPLE==0) { // unFIXME
	const Real kT=cond.kT; // .0259;
	//mjm_libmesh_util::
//		AddLHSD1( Ke_var, fe,di, qp, mplus_var, 1.0*(  grads[u_var]) ); // 
if (do_lhs) {
AddLHSD1Block( Ke_block, fe,di, qp, mplus_var, kT*(  grads[mplus_var]) ); // 
AddLHSD1Block( Ke_block, fe,di, qp, mplus_var, ( -grads[u_var]) ); //  
AddLHSD1Block( Ke_block, fe,di, qp, eminus_var,kT*( grads[eminus_var]) ); //  
AddLHSD1Block( Ke_block, fe,di, qp, eminus_var, ( grads[u_var]) ); //  
AddLHSD1Block( Ke_block, fe,di, qp, mzed_var,kT*( grads[mzed_var]) ); //  
}
AddRHSBlock(Fe_block,fe,fe,di,(*this),qp,d2u,values,grads,is_side);
// this needs the solution vector 
//AddRHSBlockMake(Fe_block,fe,fe,di,(*this),qp,d2u,values,grads,is_side);
} // decouple

} // assemble_rhs







template <class Tt,class Tv, class Tx, class Ty >
        void assemble_elem(Tt & Ke_var, Tv &  Fe_var
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi, const Elem * elem)

   {
	IdxTy bsrc=0;
	const Real d2u=find_d2u(bsrc,values,grads,qp,bi,elem);
	
/*  bool is_side=false; // kluge for noe
	IdxTy bsrc=mjm_libmesh_util::bsrc(bi,elem);
	is_side=mjm_libmesh_util::Val(bsrc);
	debug_branches(is_side,bsrc);
	(*this).is_side=is_side;
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
   	//const Real coefs[]={-1.0,kT,kT,kT}; 
   	const Real bf=::exp(values[u_var]/cond.kT);
//  	CkTy::denan(bf," bf ",StrTy(__FILE__),(__LINE__));	
	//MM_ERR(" fudd u "<<values[u_var])
	 Real d2ux=((DECOUPLE==0)?1.0:0.0)*( - cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf)); //  unFIXME
//	 Real d2ux=((DECOUPLE==0)?1.0:0.0)*( - cond.q_over_e*(rev(values[mplus_var])-rev(values[eminus_var]))); //  unFIXME
	CkTy::dnz(d2ux); // slow??? 
	const Real d2u=d2ux;
	CkTy::denan(d2u," d2u ",StrTy(__FILE__),(__LINE__));	
*/

	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
// not doing sides causes nan exposion 
//	if (!is_side) 
assemble_fex(Ke_var,Fe_var,fev,values,grads,qp,bi,d2u,coefs,bsrc,is_side);
//	mjm_libmesh_util::AddLHSD2( Ke_var,  qp,fev, di,coefs);

 const IdxTy qpsz=fes.qrule.n_points();

if (!false)if (qp<qpsz) {	 // this generates oor from denorm etc. 
	for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
		if (elem->neighbor(kluge)==libmesh_nullptr)
		{	fes.fe->reinit(elem,kluge); 
			assemble_fex(Ke_var,Fe_var,fes,values,grads,qp,bi,d2u,coefs,bsrc,is_side);
		}
	fev.fe->reinit(elem);
} // qp 

// I thought this would work, it seemed to before????
//        if ( internal_electrodes) ImposeVoltages( Ke_var,Fe_var, fev, di,cond,qp, elem);
	jelec=0;
// this is distored terribly but we can apply DD now. 
//	if (geo.contains(value,is_side,p(0),p(1)))
//	if (bsrc>0) if (is_side) 
} // assemble 
/////////////////////////////////////////////////////////////////////

// assemble using given fe and values
template <class Tt,class Tv, class Tx, class Ty,class Tc >
        void assemble_given(Tt & Ke_var, Tv &  Fe_var, Tc & fe
    ,const Tx & values,const Ty & grads, const IdxTy qp, const bool do_lhs)
   {
	const Real d2u=find_d2u_only(values,grads);
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
using namespace mjm_libmesh_util;
if (do_lhs)		AddLHSD2( Ke_var, qp,fe, di,coefs);
if (DECOUPLE==0) { // unFIXME
	const Real kT=cond.kT; // .0259;
if (do_lhs) {
AddLHSD1( Ke_var, fe,di, qp, mplus_var, kT*(  grads[mplus_var]) ); // 
AddLHSD1( Ke_var, fe,di, qp, mplus_var, ( -grads[u_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, eminus_var,kT*( grads[eminus_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, eminus_var, ( grads[u_var]) ); //  
AddLHSD1( Ke_var, fe,di, qp, mzed_var,kT*( grads[mzed_var]) ); //  
}
AddRHS(Fe_var,fe,fe,di,(*this),qp,d2u,values,grads,is_side);
} // decouple
} // given

// this should be the only one in use now, the dense vecs not used
// and the jacobian is numerical 
template <class Tt,class Tv, class Tx, class Ty,class Tc , class Tsol>
        void assemble_given_block(Tt & Ke_var, Tv &  Fe_block, Tc & fe
    ,const Tx & values,const Ty & grads, const IdxTy qp, const bool do_lhs, const Tsol & sol, const Elem * e=0 )
   {
	//const Real d2u=find_d2u_only(values,grads);
	// const Real d2uold =find_d2u_only(values,grads);
	const Real d2u=	find_d2u_sol(fe,di,qp,sol);

	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
using namespace mjm_libmesh_util;
if (do_lhs)		AddLHSD2Block( Ke_var, qp,fe, di,coefs);
if (!d2u_on_right) {
if (do_lhs)		AddLHSD2CrossBlock( Ke_var, qp,fe, di,1.0,eminus_var,u_var);
if (do_lhs)		AddLHSD2CrossBlock( Ke_var, qp,fe, di,-1.0,mplus_var,u_var);
}

const bool  dump_quadrature_info=!true;
if (dump_quadrature_info)  AddLHSD2BlockDump( Ke_var, qp , fe, di, e);

// try to couple nstead of kluging for once lol
// maybe this will converge a lot better 
//if (do_lhs)		AddLHSD2CrossBlock( Ke_var, qp,fe, di,1,eqn,var);
if (DECOUPLE==0) { // unFIXME
	const Real kT=cond.kT; // .0259;

//AddLHSD1SqBlock( Ke_var,qp,fe,di,coefs,eqn,var,sol)

// I changed ther order in the dofs here see what happens
// see comments in mjm_libmesh.h where this is impemented
// seesm to be an issue. 
if (do_lhs) {
AddLHSD1Block( Ke_var, fe,di, qp, mplus_var, kT*(  grads[mplus_var]) ); // 
//AddLHSD1SqBlock( Ke_var,qp,fe,di,kT,mplus_var,mplus_var,sol);
AddLHSD1Block( Ke_var, fe,di, qp, mplus_var, ( -grads[u_var]) ); //  
//AddLHSD1SqBlock( Ke_var,qp,fe,di,-1,mplus_var,u_var,sol);
AddLHSD1Block( Ke_var, fe,di, qp, eminus_var,kT*( grads[eminus_var]) ); //  
//AddLHSD1SqBlock( Ke_var,qp,fe,di,kT,eminus_var,eminus_var,sol);
AddLHSD1Block( Ke_var, fe,di, qp, eminus_var, ( grads[u_var]) ); //  
//AddLHSD1SqBlock( Ke_var,qp,fe,di,1,eminus_var,u_var,sol);
AddLHSD1Block( Ke_var, fe,di, qp, mzed_var,kT*( grads[mzed_var]) ); //  
//AddLHSD1SqBlock( Ke_var,qp,fe,di,kT,mzed_var,mzed_var,sol);

}
// this is not right with integration of exponential 
//AddRHSBlock(Fe_var,fe,fe,di,(*this),qp,d2u,values,grads,is_side);
AddRHSBlockMake(Fe_block,fe,fe,di,(*this),qp,d2u,values,grads,is_side,sol);
} // decouple

} // given








// return the jacobian 
template <class Tt,class Tv, class Tx, class Ty >
        void assemble_J(Tt & Je_var, Tv &  Ge_var
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi, const Elem * elem, const bool do_lhs)

   {
	IdxTy bsrc=0;
	const Real d2u=find_d2u(bsrc,values,grads,qp,bi,elem);
	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
assemble_fex(Je_var,Ge_var,fev,values,grads,qp,bi,d2u,coefs,bsrc,is_side, do_lhs);

 const IdxTy qpsz=fes.qrule.n_points();
if (!false)if (qp<qpsz) {	 // this generates oor from denorm etc. 
	for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
		if (elem->neighbor(kluge)==libmesh_nullptr)
		{	fes.fe->reinit(elem,kluge); 
			assemble_fex(Je_var,Ge_var,fes,values,grads,qp,bi,d2u,coefs,bsrc,is_side,do_lhs);
		}
	fev.fe->reinit(elem);
} // qp 

	jelec=0;
} // assemble 
/////////////////////////////////////////////////////////////////////

template < class Tx, class Ty >
        Real find_d2u( IdxTy & bsrc,
    const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi, const Elem * elem)
{
//	 bool is_side=false; // kluge for noe
	//IdxTy bsrc=mjm_libmesh_util::bsrc(bi,elem);
	bsrc=mjm_libmesh_util::bsrc(bi,elem);
	bool is_side=mjm_libmesh_util::Val(bsrc);
	debug_branches(is_side,bsrc);
	(*this).is_side=is_side;
	return find_d2u_only(values,grads);
}

//	const Real d2u=	find_d2u_sol(fe,di,qp,sol);
template < class Tx, class Ty,class Tz, class Tw > 
	Real find_d2u_sol( const Tx & fe,const Ty & di, const Tz & qp, const Tw & sol)
{
	const IdxTy n_var_dofs = di.dof_indices_var[0].size();
	Real d2ux=0;
	// true gives charge shielding but is very slow 
	const bool eq_d2u=!true;
            for (IdxTy  j=0; j<n_var_dofs; j++)
            {
				const Real phij=fe.phi[j][qp];
                const auto plj=sol(di.dof_indices_var[mplus_var][j]);
                const auto mij=sol(di.dof_indices_var[eminus_var][j]);
                const auto uj=sol(di.dof_indices_var[u_var][j]);
	if ( eq_d2u)
	{
		// this is now way off, need a factor of q here. 
   		const Real bf=::exp(uj/cond.kT);
		d2ux+=((DECOUPLE==0)?1.0:0.0)*phij*
			( - cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf));
	}
	else
	{
	 	d2ux+=((DECOUPLE==0)?1.0:0.0)*phij*
			( - cond.q_over_e*(rev(plj)-rev(mij))); 
	}
} // j

	CkTy::dnz(d2ux); // slow??? 
	const Real d2u=d2ux;
	CkTy::denan(d2u," d2u ",StrTy(__FILE__),(__LINE__));	

 return d2u;
}



template < class Tx, class Ty > 
	Real find_d2u_only( const Tx & values,const Ty & grads)
{
//	const Real kT=cond.kT; // .0259;
//   	const Real coefs[]={1.0,kT,kT,kT};
	Real d2ux=0;
	// true gives charge shielding but is very slow 
	const bool eq_d2u=!true;
	if ( eq_d2u)
	{
   		const Real bf=::exp(values[u_var]/cond.kT);
		//  	CkTy::denan(bf," bf ",StrTy(__FILE__),(__LINE__));	
		//MM_ERR(" fudd u "<<values[u_var])
		d2ux=((DECOUPLE==0)?1.0:0.0)*
			( - cond.q_over_e*(cond.mpluseq/bf-cond.eeq*bf)); //  unFIXME
	}
	else
	{
	 	d2ux=((DECOUPLE==0)?1.0:0.0)*
			( - cond.q_over_e*(rev(values[mplus_var])-rev(values[eminus_var]))); //  unFIXME
	}





	CkTy::dnz(d2ux); // slow??? 
	const Real d2u=d2ux;
	CkTy::denan(d2u," d2u ",StrTy(__FILE__),(__LINE__));	
 return d2u;

}


//	 bool is_side=false; // kluge for noe
	//IdxTy bsrc=mjm_libmesh_util::bsrc(bi,elem);



template <class Tt,class Tv, class Tx, class Ty >
        void assemble_Jblock(Tt & Je_block, Tv &  Ge_block
    ,const Tx & values,const Ty & grads, const IdxTy qp, const BoundaryInfo & bi, const Elem * elem, const bool do_lhs)

   {
	IdxTy bsrc=0;
	const Real d2u=find_d2u(bsrc,values,grads,qp,bi,elem);

	const Real kT=cond.kT; // .0259;
   	const Real coefs[]={1.0,kT,kT,kT};
	assemble_fex_block(Je_block,Ge_block,fev,values,grads,qp,bi,d2u,coefs,bsrc,is_side, do_lhs);

 const IdxTy qpsz=fes.qrule.n_points();
if (!false)if (qp<qpsz) {	 // this generates oor from denorm etc. 
	for(IdxTy kluge=0; kluge<elem->n_sides(); ++kluge)
		if (elem->neighbor(kluge)==libmesh_nullptr)
		{	fes.fe->reinit(elem,kluge); 
			assemble_fex_block(Je_block,Ge_block,fes,values,grads,qp,bi,d2u,coefs,bsrc,is_side,do_lhs);
		}
	fev.fe->reinit(elem);
} // qp 


} // assemble 










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
//	const Real kT=cond.kT; // .0259
	const Real vf=1; // rev(values[u_var]/kT);
//  	CkTy::denan(rpl," rpl ",StrTy(__FILE__),(__LINE__));	
//  	CkTy::denan(rmi," rmi ",StrTy(__FILE__),(__LINE__));	
//  	CkTy::denan(rne," rne ",StrTy(__FILE__),(__LINE__));	

    //f_vec(mzed_var) = -cond.kr+cond.kf*cond.umeinstein *(rev(pl+mi-ne));
    f_vec(mzed_var) = +cond.kr*cond.umeinstein-cond.kf*cond.umeinstein *(rpl*rmi/rne);
    //f_vec(mplus_var) = d2u+cond.ummplus*(cond.kr*rev(ne-pl) -cond.kf*rev(mi));
	// negative signs may work better here ... 
if (d2u_on_right) {
    f_vec(mplus_var) = d2u-cond.ummplus*(cond.kr*rne/rpl -cond.kf*rmi)*vf;
    //f_vec(eminus_var) = -d2u+cond.umeminus*(cond.kr*rev(ne-mi) -cond.kf*rev(pl)) 
    f_vec(eminus_var) =- d2u-cond.umeminus*(cond.kr*rne/rmi -cond.kf*rpl)/vf   
			//+jelectrode*1e+6/rmi;
			+jelectrode*jfactor/rmi; // 
}else{
    f_vec(mplus_var) = -cond.ummplus*(cond.kr*rne/rpl -cond.kf*rmi); 
    f_vec(eminus_var) =-cond.umeminus*(cond.kr*rne/rmi -cond.kf*rpl)  
			+jelectrode*jfactor/rmi; // 
}
//MM_ERR(" rhs "<<f_vec(0)<<" "<<f_vec(1)<<" "<<f_vec(2)<<" "<<f_vec(3) )
// this seems to be the right order
// MM_ERR(" rhs "<<u_var<<" "<<mzed_var<<" "<<mplus_var<<" "<<eminus_var )
}

    template<class Tx, class Ty, class Tp, class Tz,class Tsol>
    void rhsmake(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & fe,
    const Tz & dof_indices_var,const Tp &qp, const Tsol & sol) const
{ // for now just keep punting 
 if (DECOUPLE!=0) return; // unFIXME 
//    Real jelectrode=0; // (is_side)? geo.j(qpoint(0),qpoint(1)):0;
//	if (jelec!=0) jelectrode=jelec;
//    f_vec(u_var) = d2u;
	

	const IdxTy n_var_dofs = dof_indices_var[0].size();
//    const IdxTy var_i=var;
			Real rplrmirne=0;
			Real rnermi=0;
			Real rnerpl=0;
			Real rmi=0;
			Real rpl=0;
	           // no jacobian here, this is just averaging. 
            for (IdxTy  j=0; j<n_var_dofs; j++)
            {
				const Real phij=fe.phi[j][qp];
                const auto plj=sol(dof_indices_var[mplus_var][j]);
                const auto mij=sol(dof_indices_var[eminus_var][j]);
                const auto nej=sol(dof_indices_var[mzed_var][j]);
				const Real rplj=rev(plj);
				const Real rmij=rev(mij); // going oor
				const Real rnej=rev(nej);
    //            u_ += fe.phi[j][qp]*x;
                rplrmirne += phij*rplj*rmij/rnej;
                rnermi += phij*rnej/rmij;
                rnerpl += phij*rnej/rplj;
                rmi += phij*rmij;
                rpl += phij*rplj;
   //             grad_u.add_scaled (fe.dphi[j][qp], x);
            }

//	const Real kT=cond.kT; // .0259
	const Real vf=1; // rev(values[u_var]/kT);
	const Real kf=cond.kf;
	const Real kr=cond.kr;
	Real d2ueff=d2u;
	d2ueff=-cond.q_over_e*(rpl-rmi);
    f_vec(u_var) = d2ueff;
    f_vec(mzed_var) = +kr*cond.umeinstein-kf*cond.umeinstein *(rplrmirne);
	// negative signs may work better here ... 
if (d2u_on_right) {
    f_vec(mplus_var)= d2ueff-cond.ummplus*(kr*rnerpl -kf*rmi)*vf;
    f_vec(eminus_var)=- d2ueff-cond.umeminus*(kr*rnermi -kf*rpl)/vf   ;
}else{
    f_vec(mplus_var)= -cond.ummplus*(kr*rnerpl -kf*rmi); 
    f_vec(eminus_var)= -cond.umeminus*(kr*rnermi -kf*rpl)  ;
}
} //rhs new




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
//        const Real xf = p(0);
//        const Real yf = p(1);
//const fe_info & fe=fes;
// eh with the new mesh the interior electrode surfaces are now
// sides as far as libmesh is concerned... doh 
//if (!false) if (!skip_electrodes)

}

Real jcoef(const IdxTy var, const IdxTy grad) const
{
	// the eqn former needs to know that these are exponentiated
	//jv+= cond.q*cond.mumplus*pl*(grads[u_var](dir)-cond.kT*grads[mplus_var](dir));
	if ((var==mplus_var)&&(grad==u_var)) 
		return cond.q*cond.mumplus;
	if ((var==mplus_var)&&(grad==mplus_var)) 
		return +cond.q*cond.mumplus*cond.kT;
	//jv+= cond.q*cond.mueminus*mi*(grads[u_var](dir)+cond.kT*grads[eminus_var](dir));
	if ((var==eminus_var)&&(grad==u_var)) 
		return cond.q*cond.mueminus;
	if ((var==eminus_var)&&(grad==eminus_var)) 
		return -cond.q*cond.mueminus*cond.kT;

	return 0;
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
	// wtf? need to get signs etc right.. 
	jv+= cond.q*cond.mumplus*pl*(grads[u_var](dir)+cond.kT*grads[mplus_var](dir));
	jv+= cond.q*cond.mueminus*mi*(grads[u_var](dir)-cond.kT*grads[eminus_var](dir));

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
  //es.parameters.set<IdxTy>("linear solver maximum iterations") =1000;
  es.parameters.set<IdxTy>("linear solver maximum iterations") =10;

}


template< class Parameters >
 static Number ASSFUDD (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var , const CTy & cond)
{ 
//static bool once=false;  MM_ERR(" init called ") once=true; 
// this one seems to trigger
static bool once=false;  if (!once) {MM_ERR(" init called ")   once=true;  }
// return 0; 
// the electrode geometry snapped to grid is not
// the same as the default so this is fudded up 
Real v= default_geo().guess_init_value(p(0),p(1));
Real bf=1; // exp(v/cond.kT);
if( var=="lmplus") return fwd(cond.mpluseq/bf);
if( var=="lmzed") return fwd(cond.mzedeq);
if( var=="leminus") return fwd(cond.eeq*bf);
if( var=="u")
{
//Real v= default_geo().guess_init_value(p(0),p(1));
return v;
}
MM_MSG(" unknown var "<<var)

return 0;


} 

// THISASSFUDD WILL NOT FUDDING COMPUIER IN CALL FROM FUDDING NL SYSTEM FUDD 

template< class Number,class Parameters > static Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var ,
					const CTy & cond )
{
static bool once=false;  if (!once) {MM_ERR(" init called ")   once=true;  }
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

#if ASSFUDD
template< class Number,class Parameters >  Number init_value (const Point & p,
                    const Parameters & parameters,
                    const std::string & system,
                    const std::string & var  )
 {
//static bool once=false;  MM_ERR(" init called ") once=true; 
static bool once=false;  if (!once) {MM_ERR(" init called ")   once=true;  }


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
Real v= geo.guess_init_value(p(0),p(1));

return v;
}
MM_MSG(" unknown var "<<var)

return 0;

}

#endif


static system_type &  setupES(EquationSystems & equation_systems,const StrTy & name="ASSFUDD")
{
   equation_systems.add_system<system_type> (name);
  system_type & system=  equation_systems.get_system<system_type>(name);
  //system_type & system=  * equation_systems.add_system<system_type>(name);
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
/*
  // Set default parameters
  // These were chosen to match the Petsc defaults
  es.parameters.set<Real>        ("linear solver tolerance") = 1e-5;
  es.parameters.set<Real>        ("linear solver minimum tolerance") = 1e-5;
  es.parameters.set<unsigned int>("linear solver maximum iterations") = 10000;

  es.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 50;
  es.parameters.set<unsigned int>("nonlinear solver maximum function evaluations") = 10000;

  es.parameters.set<Real>("nonlinear solver absolute residual tolerance") = 1e-35;
  es.parameters.set<Real>("nonlinear solver relative residual tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver absolute step tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver relative step tolerance") = 1e-8;




*/


template <class Ty> static Ty &  setupNLES(EquationSystems & equation_systems,const StrTy & name="ASSFUDD")
{
  //equation_systems.add_system<Ty> (name);
  //Ty & system=  equation_systems.get_system<Ty>(name);

  // Here we specify the tolerance for the nonlinear solver and
  // the maximum of nonlinear iterations.
  //equation_systems.parameters.set<Real>         ("nonlinear solver tolerance")          = 1.e2;
  equation_systems.parameters.set<Real>         ("nonlinear solver tolerance")          = 1.e6;
  equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations") = 10;
  equation_systems.parameters.set<unsigned int> ("linear solver maximum iterations") = 10;


  //Ty & system=   
	equation_systems.add_system<Ty>(name);
  Ty & system=  equation_systems.get_system<Ty>(name);
  system.add_variable("u", libMesh::FIRST);
  system.add_variable("lmzed", libMesh::FIRST);
  system.add_variable("lmplus", libMesh::FIRST);
  system.add_variable("leminus", libMesh::FIRST);
MM_ERR("setupNLES   ")
// this crashes 
//	system.init();
//MM_ERR("assmeble    ")
//	system.assemble();
//MM_ERR("dun ass   ")
// called in diri crap near ctor fudd 
//equation_systems.init();

	// this will make a copy, ZZ
// 	system.configure(geo,cond);

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
    SysTy&  ass_fudd_system;
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
	static const ETy & default_geo(ETy * y=0 ) 
	{ static ETy x;static ETy *xx=&x; if (y!=0) xx=y;   return *xx; } 

bool d2u_on_right;
	
}; // mjm_libmesh_assembler



#endif // guard

