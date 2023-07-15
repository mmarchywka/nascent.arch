#ifndef MJM_SHOOTING_H__
#define MJM_SHOOTING_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"

#define USING_PETSC
#ifdef USING_PETSC
// deviate from shooting to matrix 
#include "mjm_petsc_util.h"
#endif
/*


g++ -DTEST_SHOOTING__ -Wall  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -Wno-unused-function -I.. -x c++ mjm_shooting.h


Petsc is not needed for shooting alone, remove abouve inludes etc,

g++ -DTEST_SHOOTING__ -Wall  -s=gnu++11 -gdwarf-3 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -Wno-unused-function -Wno-non-template-friend -x c++ mjm_shooting.h


g++ -DTEST_PETSC_CLASS_UTIL_MAIN__ -Wall  -std=gnu++11 -I..  -I/home/marchywka/d/petsc/petsc-3.7.3/include -I/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/include -L/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -lpetsc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib/gcc -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib  -Wl,-rpath,/home/marchywka/d/petsc/petsc-3.7.3/arch-linux2-c-debug/lib -x c++ mjm_petsc_util.h 


*/



#ifdef USING_PETSC
typedef mjm_petsc_nl_solve Petnl;
// use the petsc nl solver to solve the normal linear fd
// equations for 1D continuous current flow with reaction.
class petsc_lin_fd : public Petnl::jr_base
{
/*
AFAICT, this solves F(x)=0 and J is dF(x)/dx
This test for sqrt does seem to work and botching J will make it faile
*/
typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
public :
petsc_lin_fd(): m_size(0) {Init();}
petsc_lin_fd(const IdxTy s): m_size(s) {Init();}

class params 
{
typedef unsigned int IdxTy;
typedef double D;
public:
params(): iu(0), iz(1), ip(2), in(3) {}
int  iu,iz,ip,in;

}; // params 

void bc_simple( MySparse * pJ, MyBlock & rhs
	, const MyBlock & sol, const bool do_j, const bool do_f, const params p=params())
{
// "R" is called rhs inaccurately but is lhs-rhs.
int d=m_vars;

int  iu=p.iu; IdxTy iz=p.iz; IdxTy ip=p.ip; IdxTy in=p.in;
if (do_f)
{
// this makes J>0
rhs(iu)=sol(iu)-0;
//rhs(iz)=sol(iz)-1e18;
rhs(ip)=sol(ip)-m_p_0;
rhs(in)=sol(in)-m_n_0;
rhs(iz)=sol(iz)-m_kf/m_kr*sol(ip)*sol(in);
}
if (do_j)
{
	MySparse & J =*pJ;
	J(iu,iu)=1;
	J(ip,ip)=1;
	J(in,in)=1;
	J(iz,iz)=1; 
	J(iz,ip)=-m_kf/m_kr*sol(in);
	J(iz,in)=-m_kf/m_kr*sol(ip);
} // j 

//	iu+=d; iz+=d; ip+=d; in+=d;
	iu=m_size-d; iz=iu+1; ip=iu+2; in=iu+3;
if (do_f)
{
rhs(iu)=sol(iu)-1;
//rhs(iz)=sol(iz)-1e18;
rhs(ip)=sol(ip)-m_p_d;
rhs(in)=sol(in)-m_n_d;

rhs(iz)=sol(iz)-m_kf/m_kr*sol(ip)*sol(in);
}
if (do_j)
{
	MySparse & J =*pJ;
	J(iu,iu)=1;
	J(iz,iz)=1;
	J(ip,ip)=1;
	J(in,in)=1;
	J(iz,ip)=-m_kf/m_kr*sol(in);
	J(iz,in)=-m_kf/m_kr*sol(ip);
} // j 

}


void jr( MySparse & J, MyBlock & R, MyBlock & sol)
{
//MM_ERR(" right jr called "<<m_size)
//for (IdxTy i=0; i<m_size; ++i) J(i,i)=2.0*sol(i);
// for the linear case,these are just the coefficients. 
// u,z,p,n
const IdxTy pmax=m_points-1;
const D qe=m_qe;
const D h=m_h;
const D ih2=1.0/(h*h);
const D h2=h*h;
const D ih5=.5/h;
const D ih52vt=ih5*ih5/m_Vt;
const D ih2vt=1.0/(h*h*m_Vt);
const int d=m_vars;
//const int dz=0; // -m_vars;
//const int pz1=pmax/2;
//const int pz2=pz1+1;
int  iu=0; IdxTy iz=1; IdxTy ip=2; IdxTy in=3;
iu+=d; iz+=d; ip+=d; in+=d;
for (IdxTy point=1; point<pmax; ++point)
{
	// laplacian u = -q/e(p-n)

//	rhs(iu)+=lapG*m_Vt-qe*(sol(in)-sol(ip));

	J(iu,iu+d)=ih2; J(iu,iu)=-2*ih2; J(iu,iu-d)=ih2; J(iu,ip)=qe; J(iu,in)=-qe;
	
	const D gradG=ih5*(sol(iu+d)-sol(iu-d))/m_Vt;
	const D lapG=(sol(iu+d)+sol(iu-d)-2.0*sol(iu))*ih2/m_Vt;

	//	rhs(iz)+= lapZ+fgen/m_muz;
	J(iz,iz)=-(2.0)*ih2;
	J(iz,iz+d)=1*ih2; 
	J(iz,iz-d)=1*ih2; 
	J(iz,iz)+=m_kr/m_muz/m_Vt;
	J(iz,in)+=-m_kf*sol(ip)/m_muz/m_Vt;
	J(iz,ip)+=-m_kf*sol(in)/m_muz/m_Vt;	

//	rhs(ip)=lapP+gradP*gradG +sol(ip)*lapG-fgen/m_mup;	
	const D imupvt=1.0/(m_mup*m_Vt);
	J(ip,ip)=-2.0*ih2+lapG;
	J(ip,ip+d)=ih2+ih5*gradG; 
	J(ip,ip-d)=ih2-ih5*gradG; 
	J(ip,iz)=-m_kr*imupvt;
	J(ip,in)=m_kf*sol(ip)*imupvt;
	J(ip,ip)+=m_kf*sol(in)*imupvt;	
	J(ip,iu+d)=+ih52vt*sol(ip+d);
	J(ip,iu+d)+=-ih52vt*sol(ip-d);
	J(ip,iu-d)=-ih52vt*sol(ip+d);
	J(ip,iu-d)+=+ih52vt*sol(ip-d);

	J(ip,iu)+=-2.0*sol(ip)*ih2vt;
	J(ip,iu+d)+=sol(ip)*ih2vt;
	J(ip,iu-d)+=sol(ip)*ih2vt;

//	rhs(in)=lapN-gradN*gradG -sol(in)*lapG-fgen/m_mun;	
	J(in,in)=-2.0*ih2-lapG;
	J(in,in+d)=ih2-ih5*gradG; 
	J(in,in-d)=ih2+ih5*gradG;
	const D imunvt=1.0/(m_mun*m_Vt);
 	J(in,iz)=-m_kr*imunvt;
	J(in,in)+=m_kf*sol(ip)*imunvt;
	J(in,ip)+=m_kf*sol(in)*imunvt;	
	J(in,iu+d)=-ih52vt*sol(in+d);
	J(in,iu+d)+=+ih52vt*sol(in-d);
	J(in,iu-d)=+ih52vt*sol(in+d);
	J(in,iu-d)+=-ih52vt*sol(in-d);


	J(in,iu)+=2.0*sol(in)*ih2vt;
	J(in,iu+d)+=-sol(in)*ih2vt;
	J(in,iu-d)+=-sol(in)*ih2vt;

	iu+=d; iz+=d; ip+=d; in+=d;
}

 bc_simple( &J,  R,  sol, true, false);
MM_ONCE(" printing J disbled ",return; )
MM_MSG(" J ")
J.dump(std::cout,1,"J"); // row per line 
if (true) return; 
//J(m_size-3,m_size-3)=1;
//rhs(m_size-3)=sol(iz)-m_kf/m_kr*sol(in)*sol(ip);
J(m_size-3,m_size-3)=1;
J(m_size-3,m_size-2)=-m_kf/m_kr*sol(m_size-1);
J(m_size-3,m_size-1)=-m_kf/m_kr*sol(m_size-2);

J(1,1)=1;
J(1,2)=-m_kf/m_kr*sol(3);
J(1,3)=-m_kf/m_kr*sol(2);

J(1,1)=1;
J(0,0)=1; // this voltage is fixed;
J(m_size-4,m_size-4)=1; // this voltage is fixed;
//const D c1=-m_Vt/(sol(m_vars)-sol(0)-m_Vt);

	//rhs(2)=(sol(2)*gradG+m_muz/m_mup*gradZ+gradP);
	//rhs(ip)=(sol(ip)*gradGe+m_muz/m_mup*gradZe+gradPe);
if (false) {
	const D sf=1e-0;
J(2,2)=(sol(m_vars)-sol(0))/h/m_Vt- 1.0/h;   // p 
J(2,1)=-m_muz/m_mup/h;
J(2,1+m_vars)=m_muz/m_mup/h;
J(2,2+m_vars)=1.0/h;
J(2,0)=-sol(2)/h/m_Vt;
J(2,4)=sol(2)/h/m_Vt;
J(2).multiply(sf);
const int asee=m_size-m_vars;
const int ase=m_size-2*m_vars;
J(asee+2,ase+2)=(sol(ase+m_vars)-sol(ase+0))/h/m_Vt- 1.0/h;   // p 
J(asee+2,ase+1)=-m_muz/m_mup/h;
J(asee+2,ase+1+m_vars)=m_muz/m_mup/h;
J(asee+2,ase+2+m_vars)=1.0/h;
J(asee+2,ase+0)=-sol(ase+2)/h/m_Vt;
J(asee+2,ase+4)=sol(ase+2)/h/m_Vt;
} // false 
J(2,2)=1;
J(3,3)=1;
//J(m_size-3,m_size-3)=1; //  zthis voltage is fixed;
//J(m_size-2,m_size-2)=1; // this voltage is fixed;



// fix n with current equality 
const D gradg0=-(sol(0)-sol(m_vars))/m_Vt/h;
const D gradp0=-(sol(2)-sol(2+m_vars))/h;
const D gradn0=-(sol(3)-sol(3+m_vars))/h;
const D gradgd=-(sol(m_size-2*m_vars)-sol(m_size-m_vars))/m_Vt/h;
//const D gradpd=(sol(m_size-2*m_vars+2)-sol(m_size-m_vars+2));
//const D gradnd=(sol(m_size-2*m_vars+3)-sol(m_size-m_vars+3));
J(m_size-1,m_size-1)=1; //n 
J(m_size-2,m_size-2)=1; //n 
//`J(m_size-3,m_size-3)=1; //n 
// J anode + J cathode = 0 

   // currents are the same 
//    const D jc=-gradG*(m_mup*sol(2)+m_mun*sol(3))-m_mup*gradP+m_mun*gradN;
// -([0]-[4])*(m_mup*[2] + m_un[3])-m_up*([2]-[6])+m_un([3]-[7])
// -([m_size-8]-[m_size-4])*(m_mup*[m_size-2] + m_un[m_size-1])-m_up*([m_size-6]-[m_size-2])
//      +m_un([m_size-5]-[m_size-1])
//    const D ja=-gradGe*(m_mup*sol(ip)+m_mun*sol(in))-m_mup*gradPe+m_mun*gradNe;

	//const D jc=-gradG*(m_mup*sol(2)+m_mun*sol(3))-m_mup*gradP+m_mun*gradN;
	//const D ja=-gradGe*(m_mup*sol(ip)+m_mun*sol(in))-m_mup*gradPe+m_mun*gradNe;
const D jf=1.0e-20;
J(3,3)=1;
//    rhs(3)=jc+ja;
if (false) {
J(3,3)=(-gradg0*m_mun+m_mun)*jf;
J(3,0)=(-m_mup*sol(2)+m_mun*sol(3)/m_Vt)*jf;
J(3,4)=(m_mup*sol(2)+m_mun*sol(3))/m_Vt*jf;
J(3,2)=(-m_mup)*jf; 
J(3,6)=(m_mup)*jf;
J(3,7)=(-m_mun)*jf;

const int base=m_size-8;
const int basee=0; // m_size-4;
J(basee+3,base+3)=(-gradgd*m_mun+m_mun)*jf;
J(basee+3,base+0)=(-m_mup*sol(base+2)+m_mun*sol(base+3))/m_Vt*jf;
J(basee+3,base+4)=(m_mup*sol(base+2)+m_mun*sol(base+3))/m_Vt*jf;
J(basee+3,base+2)=(-m_mup)*jf; 
J(basee+3,base+6)=(m_mup)*jf;
J(basee+3,base+7)=(-m_mun)*jf;

}

//MyBlock search(m_size);
//J.multiply(search,sol);
//MM_MSG(" presumed search direction "<<search.to_string())
//search.catalog();
MM_MSG(" J ")
//J.dump(std::cout,0,"J"); // R format 
J.dump(std::cout,1,"J"); // row per line 

/*
J(3,3)=-gradg0*m_mun+m_mun; //n 
J(3,3+m_vars)=-m_mun; //n 
J(3,2)=-gradg0*m_mup-m_mup; //n 
J(3,2+m_vars)=+m_mup; //n 

J(3,m_size-m_vars-1)=+m_mun; //n 
J(3,m_size-1)=-gradgd*m_mun-m_mun; //n 
J(3,m_size-2-m_vars)=-m_mup; //n 
J(3,m_size-1-m_vars)=-gradg0*m_mup+m_mup; //n 
*/

}
 void initial_guess(   MyBlock & sol) {
//MM_ERR(" right initial_guess callled doh ")
sol.resize(size());
//for (IdxTy i=0; i<m_size; ++i) sol(i)=i*10+1;
//for (IdxTy i=0; i<m_size; i+=2) sol(i)=-sol(i);
// this needs to fit boundary conditions
for (IdxTy i=0; i<m_size; ++i) sol(i)=0;
const IdxTy pmax=m_points-0;
int  iu=0; IdxTy iz=1; IdxTy ip=2; IdxTy in=3;
for (IdxTy point=0; point<pmax; ++point)
{
const int d=m_vars;
sol(iu)= 1.0*point/(pmax-1);
sol(in)=m_n_b; // 1.0*1e8; // point/pmax;
//sol(iz)=1.0e8;
sol(ip)=m_p_b;// 1.00000e8;
//sol(iz)=1e18; // m_kf/m_kr*sol(in)*sol(ip);;
sol(iz)=m_kf/m_kr*sol(in)*sol(ip);

	iu+=d; iz+=d; ip+=d; in+=d;
}
//sol(2)=0; sol(3)=0; sol(94)=0;
sol(3)=m_n_0; // 1e6;
sol(2)=m_p_0; // 1e10;
sol(m_size-2)=m_p_d;
sol(m_size-1)=m_n_d;
sol(1)=m_kf/m_kr*sol(2)*sol(3);
sol(m_size-3)=m_kf/m_kr*sol(m_size-3)*sol(m_size-2);
}
// actually EVERYTHING is on the right, lhs-rhs=0
 void rhs( const MyBlock & sol, MyBlock & rhs) {
//MM_ERR(" right rhs callled doh ")
rhs.resize(size());
//for (IdxTy i=0; i<m_size; ++i) rhs(i)=sol(i)*sol(i)-(i);
//MySparse J;  // not exactly J but close 
const IdxTy pmax=m_points-1;
const D h=m_h;
const D qe=m_qe;
const D ih2=1.0/(h*h);
const D h2=h*h;
const D ih5=.5/h;
const D ihvt=1.0/h/m_Vt;
const D ih=1.0/h;
const int d=m_vars;

int dz=0; // -m_vars;
const int pz1=pmax/2;
const int pz2=pz1+1;

int  iu=0; IdxTy iz=1; IdxTy ip=2; IdxTy in=3;
//++iu; ++iz; ++ip; ++in;
	iu+=d; iz+=d; ip+=d; in+=d;
for (IdxTy point=1; point<pmax; ++point)
{

	// laplacian u = -q/e(p-n)
//	J(iu,iu+d)=ih2; J(iu,iu)=-2*ih2; J(iu,iu-d)=ih2; J(iu,ip)=qe; J(iu,in)=-qe;
	
	const D gradG=ih5*(sol(iu+d)-sol(iu-d))/m_Vt;
	const D lapG=(sol(iu+d)+sol(iu-d)-2.0*sol(iu))*ih2/m_Vt;
	const D gradZ=ih5*(sol(iz+d)-sol(iz-d));
	const D lapZ=(sol(iz+d)+sol(iz-d)-2.0*sol(iz))*ih2;
	const D gradP=ih5*(sol(ip+d)-sol(ip-d));
	const D lapP=(sol(ip+d)+sol(ip-d)-2.0*sol(ip))*ih2;
	const D gradN=ih5*(sol(in+d)-sol(in-d));
	const D lapN=(sol(in+d)+sol(in-d)-2.0*sol(in))*ih2;

//	const D den=1.0-.5*h*gradG;
//	const D num=1.0+.5*h*gradG;
//	bool fixz=((point==(pz1))||(point==pz2));
//	bool bump=((point==pz2));
	const D fgen=(sol(iz)*m_kr-m_kf*sol(ip)*sol(in))/m_Vt;

	rhs(iu)+=lapG*m_Vt-qe*(sol(in)-sol(ip));
	rhs(iz)+= lapZ+fgen/m_muz;
	rhs(ip)=lapP+gradP*gradG +sol(ip)*lapG-fgen/m_mup;	
	rhs(in)=lapN-gradN*gradG -sol(in)*lapG-fgen/m_mun;	

	iu+=d; iz+=d; ip+=d; in+=d;
}

 bc_simple( 0,  rhs,  sol, !true, !false);
 check_new_sol( rhs,  sol);
if (true) return; 
// fix voltage
rhs(0)=sol(0)-0; rhs(m_size-4)=sol(m_size-4)-1;
rhs(3)=sol(3)-1e8;
rhs(2)=sol(2)-1e8;
rhs(1)=sol(1)-m_kf/m_kr*sol(2)*sol(3);

rhs(m_size-1)=sol(m_size-1)-1e8;
rhs(m_size-2)=sol(m_size-2)-1e8;
rhs(m_size-3)=sol(m_size-3)-m_kf/m_kr*sol(m_size-1)*sol(m_size-2);

//rhs(1+m_vars*pz1)=sol(1+m_vars*pz1)-1e8;
//rhs(1+m_vars*pz2)=sol(1+m_vars*pz2)-1e8;
if (false)
{
	iu=0; iz=1; ip=2; in=3;
	const D gradG=ihvt*(sol(iu+d)-sol(iu));
	//const D lapG=(sol(iu+d)+sol(iu-d)-2.0*sol(iu))*ih2;
	const D gradZ=ih*(sol(iz+d)-sol(iz));
	//const D lapZ=(sol(iz+d)+sol(iz-d)-2.0*sol(iz))*ih2;
	const D gradP=ih*(sol(ip+d)-sol(ip));
	//const D lapP=(sol(ip+d)+sol(ip-d)-2.0*sol(ip))*ih2;
	const D gradN=ih*(sol(in+d)-sol(in));
	//const D lapN=(sol(in+d)+sol(in-d)-2.0*sol(in))*ih2;

	iu=m_size-2*m_vars; iz=iu+1; ip=iu+2; in=iu+3;
	const D gradGe=ihvt*(sol(iu+d)-sol(iu));
	//const D lapGe=(sol(iu+d)+sol(iu-d)-2.0*sol(iu))*ih2;
	const D gradZe=ih*(sol(iz+d)-sol(iz));
	//const D lapZe=(sol(iz+d)+sol(iz-d)-2.0*sol(iz))*ih2;
	const D gradPe=ih*(sol(ip+d)-sol(ip));
	//const D lapPe=(sol(ip+d)+sol(ip-d)-2.0*sol(ip))*ih2;
	const D gradNe=ih*(sol(in+d)-sol(in));
	//const D lapNe=(sol(in+d)+sol(in-d)-2.0*sol(in))*ih2;
	const D sf=1e-0;
	// the p and z things needs to balance as currnet flow 
	rhs(2)=sf*(sol(2)*gradG+m_muz/m_mup*gradZ+gradP);
	rhs(ip+m_vars)=sf*(sol(ip)*gradGe+m_muz/m_mup*gradZe+gradPe);

	// currents are the same 
	const D jc=-gradG*(m_mup*sol(2)+m_mun*sol(3))-m_mup*gradP+m_mun*gradN;
	const D ja=-gradGe*(m_mup*sol(ip)+m_mun*sol(in))-m_mup*gradPe+m_mun*gradNe;

	rhs(3)=sol(3)-1e8; // -(jc+ja)*1e-20;

rhs(m_size-3)=sol(iz)-m_kf/m_kr*sol(in)*sol(ip);
rhs(1)=sol(1)-m_kf/m_kr*sol(2)*sol(3);
rhs(m_size-1)=sol(m_size-1)-1e8;

}


//rhs(2)=0; rhs(3)=0; rhs(94)=0;
MM_MSG(" rhs cata log "<< rhs.to_string())
MM_MSG(" sol cata log "<< sol.to_string())
const D l2=rhs.norm();
MM_MSG(" rhs l2 norm  "<< l2<<" "<<m_best_l2)
if ((l2<m_best_l2)||(m_best_l2<0) )
{
m_best_l2=l2;
if (sol.size()==m_best_sol.size())
{
MyBlock xxx=sol-m_best_sol;
MM_MSG(" diff sol "<<xxx.to_string())
}
m_best_sol=sol;
}
//rhs.catalog();

}

void check_new_sol(const MyBlock & rhs, const MyBlock & sol)
{
MM_ONCE(" printing disabled for rhs and sol ", )
if (false) {
MM_MSG(" rhs cata log "<< rhs.to_string())
MM_MSG(" sol cata log "<< sol.to_string())
}

const D l2=rhs.norm();
MM_MSG(" rhs l2 norm  "<< l2<<" best  "<<m_best_l2<<" diff "<<(l2-m_best_l2)<<" iter "<<m_iter)
++m_iter;
if ((l2<m_best_l2)||(m_best_l2<0) )
{
m_best_l2=l2;
if(false) if (sol.size()==m_best_sol.size())
{
MyBlock xxx=sol-m_best_sol;
MM_MSG(" diff sol "<<xxx.to_string())
}
m_best_sol=sol;
}


}


MyBlock & solution() { return m_best_sol; }
// this is the total number of eqns
virtual IdxTy size() const { return m_size; }

private:
void Init() {
	m_best_l2=-1;
	m_vars=4;
	m_points=m_size/m_vars;
	m_h=1e-6;
	 m_Vt=(300.0*1.3807e-16/1.6e-19);
//	 m_Vt=.0259;
	m_qe=1.6e-19/8.854e-14;
	m_kf=100;
	m_kr=1e18*m_kf;
	m_mup=100;
	m_mun=1000;
	m_muz=10;
 	m_n_0=1e12;
 	m_n_d=1e16;
	m_p_0=1e16;
	m_p_d=1e12;
	m_n_b=1e14;
	m_p_b=1e14;
	m_iter=0;
 }

IdxTy m_size;
MyBlock m_best_sol;
D m_best_l2;
// misc stuff
IdxTy m_vars, m_points;
D m_Vt,m_h,m_qe,m_kf,m_kr;
D m_mup, m_muz, m_mun;
D m_n_0, m_n_d,m_p_0,m_p_d,m_n_b,m_p_b;
IdxTy m_iter=0;
};





#endif // USING_PETSC


class mjm_shooting
{

typedef mjm_shooting Myt;

public:
class Tr{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
//typedef LocTy location_type;
//typedef mjm_rational RatTy;
}; // Tr

typedef Tr::D D;
typedef Tr::SsTy SsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<IdxTy> MyBlockInt;
typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
typedef mjm_generic_itor<2> LoopItor2;
typedef mjm_generic_itor<3> LoopItor3;
typedef mjm_generic_itor<4> LoopItor4;

// really the log and conc versions should be subclasses
// with same method names doh 
mjm_shooting () {}

~mjm_shooting () {}

class state
{
public:
// h must be set or it bombs 
//  f=0; // these vary between equations, zero for now 
// in both cases, f begins as generation rate and in log case
// is divided by carrier concentration, exp(u)
 state():
	 grad_g(0),  lap_g(0), Vt(300*1.3807e-16/1.6e-19),  f(0), u_(0), u_minus_1(0), h(0),q(1)
{}
void step(const D & u)
{
u_minus_1=u_;
u_=u;
}
 D  grad_g,  lap_g, Vt,  f, u_, u_minus_1, h,q ;
}; //state

void  u_n_plus_1_orig( D & p1,  D & p2, const state & s )
{ u_n_plus_1_orig(p1,p2, s.grad_g,  s.lap_g, s.Vt,  s.f, s.u_, s.u_minus_1, s.h, s.q) ; }

// actual carrier concentration with g as u/kT and f is the rhs.
// answer is p1/p2
void  u_n_plus_1_orig( D & p1,  D & p2, const D & grad_g, const D & lap_g, const D & Vt, const D & f,
  const D & u_, const D & u_minus_1, const D & h , const D & q=1)
{
const D h2=.5*h;
const D hh=h*h;
p2=1.0-h2*grad_g*q;
const D cn=2*(1.0+h2*lap_g*q);
const D cn1=-(1.0+h2*grad_g*q);
p1=(cn1*u_minus_1+cn*u_+hh*f/Vt);
}
// residual of normal FD equation 
// note that the "f" is different in the two cases 
D  conc_residual(const D & sol , const D & grad_g, const D & lap_g, const D & Vt, const D & f,
  const D & u_, const D & u_minus_1, const D & h )
{
const D lhs= sol + u_minus_1-2.0*u_-.5*h*(sol-u_minus_1)*(grad_g)
			- h*h*u_*lap_g;
//sol*(hhi+.25*hhi*sol-.5*hhi*u_minus_1-hgg) + u_*(-2*hhi) +u_minus_1*(hhi+.25*hhi*u_minus_1+hgg);
const D rhs=f/Vt;
return rhs-lhs;
}


// log of carrier concentration with g as u/kT and f is the rhs.
// answer is p1 +/- p2 assuming p2 is real.  
void  u_n_plus_1( D & p1,  D & p2, const D & grad_g, const D & lap_g, const D & Vt, const D & f,
  const D & u_, const D & u_minus_1, const D & h, const D & q=1 )
{
const D h2g=2.0*h*grad_g*q;
const D b=4.0-2*u_minus_1 - h2g;
// in my hotes this is 4ch^2, 
const D c = -8.0*u_+u_minus_1*(4.0+u_minus_1+h2g) -4.0*h*h*((f/Vt)+lap_g*q);

 p1=-b*.5;
 p2=.5*::sqrt(b*b-4.0*c);

}
void  u_n_plus_1( D & p1,  D & p2, const state & s )
{ u_n_plus_1(p1,p2, s.grad_g,  s.lap_g, s.Vt,  s.f, s.u_, s.u_minus_1, s.h,s.q) ; }

D  log_residual( const D & sol, const state & s )
{return  log_residual(sol, s.grad_g,  s.lap_g, s.Vt,  s.f, s.u_, s.u_minus_1, s.h) ; }


// verify the solutions are right for log FD eqn
D  log_residual(const D & sol , const D & grad_g, const D & lap_g, const D & Vt, const D & f,
  const D & u_, const D & u_minus_1, const D & h )
{
const D hhi=1.0/(h*h);
const D hgg=.5*grad_g/h;
const D lhs=sol*(hhi+.25*hhi*sol-.5*hhi*u_minus_1-hgg)
		+ u_*(-2*hhi)
		+u_minus_1*(hhi+.25*hhi*u_minus_1+hgg);
const D rhs=f/Vt+lap_g;
return rhs-lhs;
}
} ;// mjm_shooting

#ifdef  TEST_SHOOTING__
// actual system
typedef double D;
typedef unsigned int IdxTy;
template <class Ts> void set_sims( Ts * states, const IdxTy n,
	const D & h, const D & gg, const D & lap_g)
{
for (IdxTy i=0; i<n; ++i)
{
states[i]->h=h;
states[i]->grad_g=gg;
states[i]->lap_g=lap_g;
}
}
template <class Ts> void set_t( Ts * states, const IdxTy n,
	const D & tv)
{
for (IdxTy i=0; i<n; ++i)
{
states[i]->Vt=tv;
}
}


template <class Ts> void set_both( Ts & lins, Ts & log, const D & u0, 
const D & uminus)
{

log.u_=u0; log.u_minus_1=uminus;
lins.u_=::exp(log.u_); lins.u_minus_1=::exp(log.u_minus_1);

}

template <class Ts> void setup_system(int argc, char **args,
	  Ts &  pn, Ts & pl,Ts & nn, Ts & nl,Ts & zn,Ts & zl,Ts & phi,
	 D & mup,  D & muz, D & mun, D & kf, D & kr)
{
const IdxTy n=7;
Ts * states[n];
states[0]=&pn; states[1]=&pl; 
states[2]=&nn; states[3]=&nl; 
states[4]=&zn; states[5]=&zl; 
states[6]=&phi; 
pn.q=1; pl.q=1;
zn.q=0; zl.q=0;
nn.q=-1; nl.q=-1;

D h=1e-9;
D grad_g=0;
D lap_g=0;
D Vt=300*1.3807e-16/1.6e-19;
// cleverly, the na and p and z are LOG concentration 
D n0=19-1e-3; D n1=19;
D p0=14+1.00001e-3; D p1=14;
D z0=log(kf/kr*exp(n0+p0));
D z1=log(kf/kr*exp(n1+p1));
// D z1=10;
D phi1=0;
// phi1 set by constraint 
for (IdxTy i=1; i<IdxTy(argc); i+=2)
{ 
	if (strcmp(args[i],"-h")==0)
	{ 	h=atof(args[i+1]); MM_MSG(" setting h to "<<h) }

}
// this is not consistent with the grad phi needed to keep stochiometry
// at electrode 
set_sims(states,n,h,grad_g,lap_g);
set_t(states,n,Vt);
// the "1" is minus 1
set_both( pn, pl, p0,p1); 
set_both( nn, nl, n0,n1); 
set_both( zn, zl, z0,z1); 
phi.u_minus_1=phi1;
// fluxes need to match reaction rate 

phi.u_= phi.u_minus_1+ Vt*(muz/mup*exp(z0-p0)*(z0-z1)+(p0-p1));
}

// first petsc example


//class petsc_lin_fd : public Petnl::jr_base
static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
int main(int argc,char **args)
{

typedef mjm_petsc_nl_solve Myt;
//PetscInitialize(&argc,&args,(char*)0,help);
int  m_ierr = PetscInitialize(&argc,&args,(char*)0,help);if (m_ierr) return m_ierr;
typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
//  PetscErrorCode ierr;
//  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
const IdxTy r=100;

 petsc_lin_fd jr(r); 
//test_nl jr(r);
Myt x( &jr);
x.set(500000,Myt::MAXIT);
x.set(10000000,Myt::MAXF);
x.set(1e-19,Myt::RTOL);
x.set(1e-39,Myt::STOL);
//x.init();
//typedef double D;
//typedef mjm_block_matrix<D> MyBlock;
//typedef mjm_sparse_matrix<D> MySparse;
x.solve();
MyBlock & zx=x.solution();
MyBlock & y=jr.solution();
for (IdxTy i=0; i<y.size(); ++i) MM_MSG(" sol "<<i<<" "<<y(i)<<" "<<zx(i))
return 0;

} //petsc main


// last shooting example

int main3(int argc,char **args)
{
typedef  mjm_shooting Ty;
typedef  mjm_shooting::state Tystate;
typedef Ty::IdxTy IdxTy;
typedef Ty::D D;
Ty x;
//Tystate sn,sl;
Tystate pn,pl,nn,nl,zn,zl,phi;
// these are not const as can be changed by params
D mup=100;
D muz=10;
D mun=1000;
const D qe=1.6e-19/8.854e-14 ; // /1.6e-19;
D kf=100e3;
D kr=10e3;

setup_system(argc, args,  pn,pl,nn,nl,zn,zl,phi,mup,muz,mun,kf,kr);
// argh, repeated code 
const IdxTy n=7;
Tystate * states[n];
states[0]=&pn; states[1]=&pl; 
states[2]=&nn; states[3]=&nl; 
states[4]=&zn; states[5]=&zl; 
states[6]=&phi; 

D generation_rate=0;

D p1=0;
D p2=0;
D p1_exp=0;
D p2_exp=0;
const D dist=1;
const IdxTy iend=dist/nl.h;
for ( IdxTy i=0; i<iend; ++i)
{
	x.u_n_plus_1( p1,  p2,  nl );
	x.u_n_plus_1_orig ( p1_exp, p2_exp, nn );
	nl.step(p1+p2);
	// this needs to be divided by mobility etc
	nl.f=generation_rate*::exp(-nl.u_)/mun;
	nn.step(p1_exp/p2_exp);
	nn.f=generation_rate/mun;

	x.u_n_plus_1( p1,  p2,  pl );
	x.u_n_plus_1_orig ( p1_exp, p2_exp, pn );
	pl.step(p1+p2);
	pl.f=generation_rate*::exp(-pl.u_)/mup;
	pn.step(p1_exp/p2_exp);
	pn.f=generation_rate/mup;

	// we can zero out the g terms or make a special one...
	x.u_n_plus_1( p1,  p2,  zl );
	x.u_n_plus_1_orig ( p1_exp, p2_exp, zn );
	zl.step(p1+p2);
	zl.f=-generation_rate*::exp(-zl.u_)/muz;
	zn.step(p1_exp/p2_exp);
	zn.f=-generation_rate/muz;
	//const D lap_g=qe*(pn.u_-nn.u_);
	const D lap_g=qe*(exp(pl.u_)-exp(nl.u_));
 D 	phi_next=2.0*phi.u_-phi.u_minus_1-phi.h*phi.h*lap_g;
MM_MSG("soln "<<(pn.h*i)<<" "<<phi.u_<<" "<<exp(zl.u_)<<" "<<exp(pl.u_)<<" "<<exp(nl.u_)<<" q "<<(lap_g/qe))
generation_rate=kf*exp(nl.u_)*exp(pl.u_)-exp(zl.u_)*kr;
// 

// this should have a version without h or make it track fwd and back h
//set_sims(states,n,h,grad_g,lap_g);
set_sims(states,n,phi.h,.5/phi.h*(phi_next-phi.u_)/phi.Vt,lap_g/phi.Vt);
	phi.step(phi_next);

} // i 

return 0;

}


// debug gthe linear and log thing 
int main2(int argc,char **args)
{
typedef  mjm_shooting Ty;
typedef  mjm_shooting::state Tystate;
typedef Ty::IdxTy IdxTy;
typedef Ty::D D;
Ty x;
Tystate sn,sl;
sl.h=.00001;
sl.grad_g=-5e6*1e0;
sl.lap_g=0;
sn=sl;
sl.u_=30;
sl.u_minus_1=sl.u_-1e-3;
sn.u_=::exp(sl.u_);
sn.u_minus_1=::exp(sl.u_minus_1);
if (argc>1)
{
sl.h=atof(args[1]);
sn.h=sl.h;
MM_MSG(" setting h to "<<sl.h)
}

D generation_rate=0;
D p1=0;
D p2=0;
D p1_exp=0;
D p2_exp=0;
const D dist=1;
const IdxTy iend=dist/sl.h;
for ( IdxTy i=0; i<iend; ++i)
{
	x.u_n_plus_1( p1,  p2,  sl );
	const D r_plus=x.log_residual( p1+ p2, sl );
	const D r_minus=x.log_residual( p1- p2,  sl );

	x.u_n_plus_1_orig ( p1_exp, p2_exp, sn );

	MM_MSG(( sl.h*i )<<" "<<sl.u_<<" exp "<<::exp(sl.u_)<<" orig "<<(sn.u_)<<" del "<<(::exp(sl.u_)-sn.u_)<<" "<<sl.u_minus_1<<" "<<p1<<" +/- "<<p2<<" u- "<<(p1-p2)<<" u+ "<<(p1+p2)<<" R "<<r_plus<<" "<<r_minus)
	sl.step(p1+p2);
	sl.f=generation_rate*::exp(-sl.u_);
	sn.step(p1_exp/p2_exp);
	sn.f=generation_rate;

} // for

return 0;
}  //main

// initial testing 
int main1(int argc,char **args)
{
typedef  mjm_shooting Ty;
//typedef  mjm_shooting::state Tystate;
typedef Ty::IdxTy IdxTy;
typedef Ty::D D;
Ty x;

D p1=0;
D p2=0;
D p1_exp=0;
D p2_exp=0;

const D Vt=300*1.3807e-16/1.6e-19;
D u_=30;
D u_minus_1=u_-1e-3;
D u_exp=::exp(u_);
D u_minus_1_exp=::exp(u_minus_1);

D grad_g=5e4;
D lap_g=0;
D f=0; // these vary between equations, zero for now 
D h=.00001;


for ( IdxTy i=0; i<100; ++i)
{
	x.u_n_plus_1( p1,  p2,  grad_g,  lap_g,  Vt,  f, u_, u_minus_1, h );
	const D r_plus=x.log_residual( p1+ p2,  grad_g,  lap_g,  Vt,  f, u_, u_minus_1, h );
	const D r_minus=x.log_residual( p1- p2,  grad_g,  lap_g,  Vt,  f, u_, u_minus_1, h );

	x.u_n_plus_1_orig
		( p1_exp, p2_exp, grad_g, lap_g,  Vt,  f, u_exp, u_minus_1_exp,h );

	MM_MSG(i<<" "<<u_<<" exp "<<::exp(u_)<<" orig "<<(u_exp)<<" del "<<(::exp(u_)-u_exp)<<" "<<u_minus_1<<" "<<p1<<" +/- "<<p2<<" u- "<<(p1-p2)<<" u+ "<<(p1+p2)<<" R "<<r_plus<<" "<<r_minus)
	u_minus_1=u_;
	u_minus_1_exp=u_exp;
	u_=p1+p2;
	u_exp=p1_exp;
//	u_=p1-p2;


} //i

return 0; 
}


#endif // TEST_SHOOTING__


#endif

