#ifndef MJM_DD_BASE2_H__
#define MJM_DD_BASE2_H__
 
#include "mjm_text_data.h"
#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
//#include <mjm_templates.h>
// needed for file io crap
#include <mjm_csv_ops.h>
// finally added rule hits
//#include <mjm_sequence_hits.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>
extern "C" {
void dgtsv_(int * N,int * NRHS, double * subdiag, double * diag,
 double * superdiag, double * solution,int * LDB,int * INFO);

void dgbsv_(int * N,int * KL, int * KU, int * NRHS, double * AB, int * LDAB,
int * IPIV, double * b,int * LDB,int * INFO);


};


namespace mjm_dd_base 
{

class dummy2;
typedef dummy Told;
typedef dummy2 Tnew;

const bool verbose_crap=!false;

//typedef <class Ty> bool isfinite(const Ty & x) { return std::isfinite(x); }
 bool isfinite(const double & x) { return std::isfinite(x); }

class material_state
{
public:
material_state(const Grid & g) : v(g),n(g),p(g),Nteff(g) {}
Values  v,n,p,Nteff;
void load(const StrTy & fn) { Load(fn); }
void save(const StrTy & fn) const { Save(fn); }


void resample()
{
const Grid & g=n.grid();
Grid gn(g.min(),g.max(),g.size()<<1);
Values x(gn);
x=v.resample(gn); v=x;
x=n.resample(gn); n=x;
x=p.resample(gn); p=x;
x=Nteff.resample(gn); Nteff=x;



} // resample 
private:
typedef mjm_csv::mjm_csv_ops CsvTy;
typedef mjm_csv::LineTy LineTy;
static const LineTy  header()  
{
LineTy x;
x.push_back("v");
x.push_back("n");
x.push_back("p");
x.push_back("Nteff");

return x; 
}
static void Load(Value & d, const StrTy & v)
{ SsTy ss(v.c_str()); ss>>d; }

// use a csv format with lots of precision
// for now just do it by hand
void Load(const StrTy & fn) {
const IdxTy sz=v.grid().size();
// get header line
const IdxTy maxline=(1<<16);
const IdxTy maxchar=maxline-2;
ChTy line[maxline];
IsTy * is;
bool dodel=false;
CsvTy::open_input(is,dodel,fn);
LineTy head ;
bool x=CsvTy::load_line_from_string(head, is, line, maxchar);
// now need to sort fields into members
if ( head!=header())
{
MM_MSG("Warning loading "<<fn<<" header starting "<<((head.size()!=0)?head[0]:StrTy(" ") )<<"does notmatch prototype "<<x)
}
IdxTy loaded=0;
while (CsvTy::ok(*is))
{
LineTy lt;
// this should check for comments etc
bool x=CsvTy::load_line_from_string(lt, is, line, maxline-2);
if ( !x) break;
// shoudl have some easy way to make all doubles
if ( loaded<sz)
{
const IdxTy lsz=lt.size();
if ( lsz>0) Load(v[loaded],lt[0]);
if ( lsz>1) Load(n[loaded],lt[1]);
if ( lsz>2) Load(p[loaded],lt[2]);
if ( lsz>3) Load(Nteff[loaded],lt[3]);
}


++loaded;
} // ok
if ( loaded!=v.grid().size())
{
MM_MSG("Warning loading "<<fn<<" size" <<loaded <<" not match grid "<<v.grid().size())
}
if ( dodel) delete is;


 } //load


void Save(const StrTy & fn) const {
const StrTy sep=" ";
//const IdxTy maxline=(1<<16);
//const IdxTy maxchar=maxline-2;
//ChTy line[maxline];
OsTy * os;
bool dodel=false;
CsvTy::open_output(os,dodel,fn);
LineTy head =header();
const IdxTy hsz=head.size();
if ( head.size()!=0) (*os)<<head[0];
for ( IdxTy i=1; i<hsz;  ++i) (*os)<<sep<<head[i];
(*os)<<CRLF;
const IdxTy sz=v.grid().size();
const IdxTy prec=32;
for ( IdxTy i=0; i<sz; ++i)
{
SsTy ss;
// this should probably dump the grid into too and a comment etc
ss.precision(prec);
ss<<v[i]<<sep;
ss.precision(prec);
ss<<n[i]<<sep;
ss.precision(prec);
ss<<p[i]<<sep;
ss.precision(prec);
ss<<Nteff[i]<<CRLF;

(*os)<<ss.str();


} //i 


if ( dodel) delete os;
 } // save


}; //state

class material_properties
{
public:

material_properties() :  Nc(0),Nv(0),ni(0),ni2(0),egap(0),epsilon(0),mu(0){}
Value Nc,Nv,ni,ni2,egap,epsilon,mu;
void Si()
{
 Nc=2.8e19;
 Nv=1.83e19;
 ni=1.4e10;
 ni2=ni*ni;
 egap=1.1;
epsilon=12;
mu=100;
}



} ; // properties





class JunkBin
{
public:
JunkBin () : iter(0),terminate_code(0),last_v(0) {}
IdxTy  iter; 
IdxTy terminate_code;
Value last_v; // for checking rate of V change
}; //JunkBin

class sim_flags
{
public:
sim_flags():vfold(.1),nfold(.1), sep(" ")
,resample_interval(0), resample_max(0),fixed_start(3),fixed_end(2)
//,solver(0), non_linear_fac(.01)
//,solver(0), non_linear_fac(2e7)
,solver(0), non_linear_fac(1),maxstep(20),minstep(-20),maxnn(1e20),minnn(1e-2)
 {}
bool resample( const IdxTy iter) const 
{ 
if ( resample_max<=iter ) return false;
if ( (iter%resample_interval) == (resample_interval-1)) return true; 
return false;

 } 
Value vfold; // amount of old voltage to mix in with updates
Value nfold; // amount of old carriers to mix in with updates
StrTy sep; // default io word separater
IdxTy resample_interval,resample_max; // 2* resample times
// number of samples at each end to exclude from size and remain fixed
IdxTy fixed_start, fixed_end; //sf=
IdxTy solver;
// scaling of amount to progress along grad due to failure to
// include e field to maintain band matrix,,,,,
Value non_linear_fac;
Value maxstep,minstep,maxnn,minnn;
}; // sim_flags

typedef material_state MatSt;
typedef material_properties MatPr;
typedef GeneralParameters GP;

class current_components
{

public:
current_components() {}
// eventually this iwill be a collection of carrier types,
// n and p ok for now 
Values n_drift,n_diff, p_drift,p_diff,j_totx,j_toty;  

 void compute( const MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
Values en(g),gradn(g),gradp(g),totn(g),totp(g);
Values xn(g),xp(g);
Told::D(en,g,m.v); // negative of field
Told::D(gradn,g,m.n);
Told::D(gradp,g,m.p);
n_drift=Grid(g);
n_diff=Grid(g);
p_drift=Grid(g);
p_diff=Grid(g);
Told::D(totn,g,m.v*m.n);
Told::D(totp,g,m.v*m.p);
Told::D(xn,g,m.n.ln());
Told::D(xp,g,m.p.ln());
// better approx to grad n 
j_totx=totn-gradn*m.v-(((xn*m.n)*(gp.vt())));
j_totx=j_totx+totp-gradp*m.v-(((xp*m.p)*(0-gp.vt())));

j_toty=totn-gradn*(m.v+gp.vt());
j_toty=j_toty+totp-gradp*(m.v-gp.vt());


//j_totx=totn-((gradn*(m.v+gp.vt())));
//j_totx=j_totx+totp-((gradp*(m.v-gp.vt())));


// the field is negative 
const Value a1=-gp.q()*pr.mu;
const Value a2=gp.q()*gp.vt()*pr.mu;
for ( IdxTy i=0; i<(sz); ++i)
{
n_drift[i]=en[i]*m.n[i]*a1;
p_drift[i]=en[i]*m.p[i]*a1;

n_diff[i]=a2*gradn[i];
//n_diff[i]=a2*xn[i]*m.n[i];
p_diff[i]=-a2*gradp[i];
//p_diff[i]=-a2*xp[i]*m.p[i];
j_totx[i]=j_totx[i]*a1;
j_toty[i]=j_toty[i]*a1;

} // i 


} //compute

}; // current_components


class dummy2
{

public:
// make a single N type section that has rhs intrinsic, see if 1/2
// section integrates right 
static void make_one_side(MatSt & m, const MatPr & pr, const Value & Nd, const Value & Na, const Value junction, const GP & gp )
{

const Value ni=pr.ni;
const Value ni2=pr.ni2;
const Value Nc=pr.Nc;
const Value Nv=pr.Nv;
const Value Ec=pr.egap;
const Value Ev=0;
const Value eps=pr.epsilon*gp.epsilon_vacu();
const Value kT=gp.vt();

// the N material is on the left, 
const Value efn=Told::fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efp=Told::fermi_level( Na, Ev+.001,  0,0 ,  Nc,  Nv,  kT , Ec,Ev);
const Value builtin=(efn-efp);  

const Grid & g=m.n.grid();
const Value jun=junction;
const Value sz=g.size();
//const Value dx=g.min_dx();
//const Value dxf=dx*dx/eps*gp.q();
const Value xp=::sqrt(2.0*builtin*Nd/(Na*(Na+Nd))*eps/gp.q());
const Value xn=::sqrt(2.0*builtin*Na/(Nd*(Na+Nd))*eps/gp.q());
MM_MSG("builtin= "<<builtin<<" efp= "<<efp<<" efn= "<<efn<<" xp= "<<xp<<" xn= "<<xn)

const Value Nd2=Nd/2;
const Value Na2=Na/2;
const Value Ndeff=Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Npdeff=-Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Naeff=Na2+::sqrt(Na2*Na2+ni2);
const Value Nnaeff=-Na2+::sqrt(Na2*Na2+ni2);
bool last=true;
for ( IdxTy i=0; i<(sz); ++i)
{
const bool left=(g[i]<junction);
Value Nt=(left)?(-Nd):(Na);
m.Nteff[i]=Nt;
if ( (!left)&&last) m.Nteff[i]=0;
last=left;
} // i 
Values & p=m.p;
Values & n=m.n;
Values & v=m.v;
Value qini=0;
bool was_depleted=false;
for ( IdxTy i=0; i<sz; ++i)
{
const bool left=(g[i]<jun);
const Value dn=jun-xn-g[i];
const Value dp=jun+xp-g[i];
const bool depleted= ( left)?(g[i]>(jun-xn)):(g[i]<(jun+xp));
if ( !depleted){
if ( left) { n[i]=Ndeff; p[i]=Npdeff;v[i]=0; qini+=n[i]-p[i]-Nd;   }
else { p[i]=Naeff; n[i]=Nnaeff; v[i]=-builtin; qini+=n[i]-p[i]+Na;}
} // depleted 
if ( depleted) { n[i]=ni; p[i]=ni; 
if ( left  ) {v[i]=-.5*gp.q()*Nd/eps*dn*dn; qini-=Nd; } 
else { v[i]=-builtin+.5*gp.q()*Na/eps*dp*dp; qini+=Na; } 
} // depleted
// dump the charge here for now 
if ( !depleted&&was_depleted)
{
// simplyusing these values kept charge at 1e12 instead of creeping
// to 1e16 in a few iters. 
const Value Na2x=(Na+3e-13*qini)/2.0;
const Value Naeffx=Na2x+::sqrt(Na2x*Na2x+ni2);
const Value Nnaeffx=-Na2x+::sqrt(Na2x*Na2x+ni2);
p[i-1]=Naeffx; n[i-1]=Nnaeffx;
//p[i]+=qini;
} // neutralize
was_depleted=depleted;
} //i 


} // make_one_side



// junction is now a coordinate, not an index into g
static void make_abrupt(MatSt & m, const MatPr & pr, const Value & Nd, const Value & Na, const Value junction, const GP & gp )
{
const bool zero_dep=!true;
const Value ni=pr.ni;
const Value ni2=pr.ni2;
const Value Nc=pr.Nc;
const Value Nv=pr.Nv;
const Value Ec=pr.egap;
const Value Ev=0;
const Value eps=pr.epsilon*gp.epsilon_vacu();
const Value kT=gp.vt();

// the N material is on the left, 
const Value efn=Told::fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efp=Told::fermi_level( Na, Ev+.001,  0,0 ,  Nc,  Nv,  kT , Ec,Ev);
const Value builtin=(efn-efp);  

const Grid & g=m.n.grid();
const Value jun=junction;
const Value sz=g.size();
//const Value dx=g.min_dx();
//const Value dxf=dx*dx/eps*gp.q();
const Value xp_depl=::sqrt(2.0*builtin*Nd/(Na*(Na+Nd))*eps/gp.q());
const Value xn_depl=::sqrt(2.0*builtin*Na/(Nd*(Na+Nd))*eps/gp.q());

const Value xp=(zero_dep)?0:xp_depl;
const Value xn=(zero_dep)?0:xn_depl;


MM_MSG("builtin= "<<builtin<<" efp= "<<efp<<" efn= "<<efn<<" xp= "<<xp<<" xn= "<<xn<<" "<<zero_dep)

const Value Nd2=Nd/2;
const Value Na2=Na/2;
const Value Ndeff=Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Npdeff=-Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Naeff=Na2+::sqrt(Na2*Na2+ni2);
const Value Nnaeff=-Na2+::sqrt(Na2*Na2+ni2);
for ( IdxTy i=0; i<(sz); ++i)
{
const bool left=(g[i]<junction);
Value Nt=(left)?(-Nd):(Na);
m.Nteff[i]=Nt;
} // i 
Values & p=m.p;
Values & n=m.n;
Values & v=m.v;
Value qini=0;
bool was_depleted=false;
for ( IdxTy i=0; i<sz; ++i)
{
const bool left=(g[i]<jun);
const Value dn=jun-xn-g[i];
const Value dp=jun+xp-g[i];
const bool depleted= ( left)?(g[i]>(jun-xn)):(g[i]<(jun+xp));
if ( !depleted){
if ( left) { n[i]=Ndeff; p[i]=Npdeff;v[i]=0; qini+=n[i]-p[i]-Nd;   }
else { p[i]=Naeff; n[i]=Nnaeff; v[i]=-builtin; qini+=n[i]-p[i]+Na;}
} // depleted 
if ( depleted) { n[i]=ni; p[i]=ni; 
if ( left  ) {v[i]=-.5*gp.q()*Nd/eps*dn*dn; qini-=Nd; } 
else { v[i]=-builtin+.5*gp.q()*Na/eps*dp*dp; qini+=Na; } 
} // depleted
// dump the charge here for now 
if ( !depleted&&was_depleted)
{
// simplyusing these values kept charge at 1e12 instead of creeping
// to 1e16 in a few iters. 
const Value Na2x=(Na+3e-13*qini)/2.0;
const Value Naeffx=Na2x+::sqrt(Na2x*Na2x+ni2);
const Value Nnaeffx=-Na2x+::sqrt(Na2x*Na2x+ni2);
p[i-1]=Naeffx; n[i-1]=Nnaeffx;
//p[i]+=qini;
} // neutralize
was_depleted=depleted;
} //i 

}  //make_abrupt




static void dump_crap( OsTy & os, const MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
const StrTy sep=flags.sep; // flags?
const StrTy label="dump_crap";
const IdxTy iter=junk.iter;
os<<label<<sep<<"iter"<<sep<<"i"<<sep<<"x"<<sep<<"voltage"<<sep<<"rho"<<sep<<"netQ"<<sep<<"n"<<sep<<"h";
os<<sep<<"ndrift"<<sep<<"ndiff"<<sep<<"pdrift"<<sep<<"pdiff"<<sep<<"jtot"<<sep<<"jtotx"<<sep<<"jtoty";
os<<CRLF;
Value qnet=0;
current_components cc;
cc.compute(m,pr,gp,flags,junk);
for (IdxTy i=0; i<sz; ++i)
{
const Value rho=-m.Nteff[i]-m.n[i]+m.p[i];
qnet+=rho;
const Value jtot=cc.n_drift[i]+cc.n_diff[i]+cc.p_drift[i]+cc.p_diff[i];
os<<label<<sep<<iter<<sep<<i<<sep<<g[i]<<sep<<m.v[i]<<sep<<rho<<sep<<qnet<<sep<<m.n[i]<<sep<<m.p[i];
os<<sep<<cc.n_drift[i]<<sep<<cc.n_diff[i]<<sep<<cc.p_drift[i]<<sep<<cc.p_diff[i]<<sep<<jtot<<sep<<cc.j_totx[i]<<sep<<cc.j_toty[i];
os<<CRLF;

}


} // dump_crap
static void measure( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{



}

// sum rho dot V and D dot E with arbitrary V origin using nx and px
// instead of n and p in m.
static Value energy( MatSt & m, const Values & nx, const Values & px, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
Values  rho(g),e(g),v(g),temp(g);
// use the existing v? no for now  
rho=nx-px+m.Nteff;
Told::I(e,g,rho);
Told::I(v,g,e);
Told::I(temp,g,e*e);
const Value d_dot_e=temp[sz-1];
Told::I(temp,g,rho*v);
const Value v_dot_rho=temp[sz-1];

// cgs to ev? Something still not right i- factor of 1e7,
//return gp.q()*gp.q()*(d_dot_e-v_dot_rho)/pr.epsilon;
const Value eps=pr.epsilon*gp.epsilon_vacu();
return gp.q()*(d_dot_e+v_dot_rho)/eps*1e-7;

}

static Value energy2( Values * grade, const MatSt & m, const Values & nx, const Values & px, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
Value energy=0;

const Grid & g=m.n.grid();
const Value sz=g.size();
const Value eps=pr.epsilon*gp.epsilon_vacu();
const Value ei=1.0/eps; // pr.epsilon;
Values  rho(g),e(g),v(g),temp(g);
// use the existing v? no for now  
rho=nx-px+m.Nteff;
Told::I(e,g,rho); // his is negatice 
Told::I(v,g,e); // v is right since rho is negative jj

if (grade!=0 )
{
Values & gr= *grade;
const Value vref=v[sz-1];
Values irhox(g),irho(g);
//Told::I(irho,g,rho);
Told::I(irhox,g,rho*g);

// this is the grad wrt RHO not n
for (IdxTy i=0; i<sz; ++i)
{
// this is neatgive but so is rho get the signs consistent
// this should be grad rho but div charge, positve for increasing "n"
// mks
gr[i]=gp.q()*(v[i]-vref)*ei*100.0; // wtf, eps is farads/cm,  
//gr[i]=(v[i]-vref);
// distance etc are cgs, 
//gr[i]+=1e-7*gp.q()*(v[i]+ei*(irhox[sz-1]-irhox[i]-g[i]*(e[sz-1]-e[i])));
gr[i]+=100.0*gp.q()*(v[i]*ei+ei*(irhox[sz-1]-irhox[i]-g[i]*(e[sz-1]-e[i])));

}

return energy;
}
 
Told::I(temp,g,e*e);
const Value d_dot_e=.5*temp[sz-1];
Told::I(temp,g,rho*v);
const Value v_dot_rho=-temp[sz-1];

// cgs to ev? Something still not right i- factor of 1e7,
//return gp.q()*gp.q()*(d_dot_e-v_dot_rho)/pr.epsilon;
//energy=gp.q()*(d_dot_e+v_dot_rho)*ei*1e-7;
energy=gp.q()*(d_dot_e*gp.q()+v_dot_rho*gp.q())*ei*100.0;
return energy;
}



static Value entropy( MatSt & m, const Values & nx, const Values & px, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
const Value ni2=pr.ni2;
Value s=0;
Value sp=0;
Value N=0;
for ( IdxTy i=0; i<(sz); ++i)
{
const Value nii=nx[i];
if ( nii == 0) continue;
N+=nii;
s-=nii*(::log(nii)-1.0);
const Value np=ni2/nii;
sp-=np*(::log(np)-1.0);

}
return s+sp+0*N*(::log(N)-1.0); 
}

// as above, but more compute der etc too.
// use nx and px, NOT the values in m but presume ni2 is usable
// and n*p == ni2 
// note that grad RHO is NOT grad n as np product is maintained.

static Value entropy2( Values * gradrho, const MatSt & m, const Values & nx, const Values & px, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
const Value dxx=1; // g.min_dx();
const Value dxxi=1.0/dxx;
const Value ni2=pr.ni2;
// this needs to belog P, not log of number of microstates
// so it needs an accurate denomiantor.
const IdxTy si=flags.fixed_start;
const IdxTy sf=sz-flags.fixed_end;
const IdxTy szeff=sf-si;
const Value normalize=::log(1.0*szeff);
//const Value qi=-1.0/gp.q();
if ( gradrho!=0)
{
Values & gr= * gradrho;
// this is the derivative wrt rho not n this way the constant rho
// constraint can be enforced and use chain rule back to grad n
Value ntot=0; Value ptot=0;

for ( IdxTy i=si; i<(sf); ++i) { ntot+=nx[i]; ptot+=px[i]; } 
// entropy is log P, need the total number of configs... 
if ( ntot!=0 ) ntot=::log(ntot)-normalize;
if ( ptot!=0 ) ptot=::log(ptot)-normalize;

for ( IdxTy i=si; i<(sf); ++i)
{
Value ns=0;
Value ps=0;
const Value n=nx[i];
if (n!=0) ns=::log(n);
const Value p=px[i];
if (p!=0) ps=::log(p);
// check the sign here...
// just make grad_rho to mean net particule count, div by charge
//gr[i]=qi*(ns*(n*n/(n*n+ni2))-ps*(p*p/(p*p+ni2)));
gr[i]=-dxxi*(((-ntot+ns)/(n*n/(n*n+ni2)))-(ni2/(n*n)*(ptot-ps)/(p*p/(p*p+ni2))));


} // i 
return 0; // later we want to combine and save the logs for both s and grad.
}


Value s=0;
Value sp=0;
Value N=0;
Value P=0;
for ( IdxTy i=si; i<(sf); ++i)
{

const Value pii=px[i];
if ( pii!=0)
{
P+=pii;
s-=pii*(::log(pii)-1.0);

}
const Value nii=nx[i];
if ( nii == 0) continue;
N+=nii;
s-=nii*(::log(nii)-1.0);
//const Value np=ni2/nii;
//sp-=np*(::log(np)-1.0);


} // i 
if ( N==0) N=1;
if ( P==0) P=1;
// this should probably have the const log M or log (sz) somewhere
// since it is multipled by N and P right? 
//return dxxi*(s+sp+N*(::log(N)-1.0)+P*(::log(P)-1.0)); 
return dxxi*(s+sp+N*(::log(N/szeff)-1.0)+P*(::log(P/szeff)-1.0)); 
//return s+sp; 

} // entropy2 

// change the grad and remove components that change total rho
static void  fix_grad( Values & gradn, const Values & gradrho, const MatSt & m, const Values & nx, const Values & px, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk,const bool do_grad_n, const bool do_log_der )
{
const Grid & g=m.n.grid();
const Value sz=g.size();
if ( sz==0) return; 
const IdxTy si=flags.fixed_start;
const IdxTy sf=sz-flags.fixed_end;
const IdxTy szeff=sf-si;
const Value ni2=pr.ni2;
//const Value qi=-gp.q();
Value sum=0;
for ( IdxTy i=si; i<(sf); ++i) { sum+=gradrho[i]; } //i 
sum=sum/szeff;
// no obvious effect?
//sum=0;

// at zero T, this seems to have a problem with low N region,
// could either do grad log N or "greater" 
// grad log N may be better idea 
if ( !do_grad_n)
{
for ( IdxTy i=si; i<(sf); ++i) 
{

const Value pp2=px[i]*px[i];
if ( pp2!=0) { 
const Value drhodn=pp2/(pp2+ni2);
gradn[i]=(gradrho[i]-sum)*drhodn;
} else gradn[i]=0;

} // i 

return;
} 

for ( IdxTy i=si; i<(sf); ++i) 
{ 
const Value nn2=nx[i]*nx[i];
if ( nn2!=0)
{
// this should change derivative from -rho/q to actual n
// which is offset byt the change in ni2/n == p 
// this should probvably be Positive ? 
const Value drhodn=-nn2/(nn2+ni2);
gradn[i]=(gradrho[i]-sum)*drhodn;
} else gradn[i]=0;

} //i 





}
// if a is too big, negtive values can not be fixed... 
static Value update_from_grad(Values & nnew, Values & pnew, const Value a, const Values & gradF, const Value sz, const MatSt & m, const Value ni2, const sim_flags & flags, const bool do_neutral)
{
Value rhotot=0;
Value logrhotot=0;
Value S=0;
Value R=0;
Value N=0;
const IdxTy si=flags.fixed_start; // 3;
const IdxTy sf=sz-flags.fixed_end; // 2;
// this worked great with entropy turned off, need a better way to
// distribute charge
//const bool  attempte_to_neutralize_log=!true;
const bool  attempte_to_neutralize_log=do_neutral;
// dealing with underflow fudded this up 
const bool  attempte_to_neutralize=!true;
const bool  attempte_log_update=true;
const Value def=::sqrt(::sqrt(ni2));
if ( !attempte_log_update)
{
// def is a miserable idea, look at junction lol 
for ( IdxTy i=0; i<(sz); ++i)
{
const bool var_region= ((i>=si)&&(i<sf)) ;
if ( true) nnew[i]=m.n[i]+a*gradF[i]; else nnew[i]=m.n[i];
//if ( nnew[i]<0) nnew[i]=def;
// this is lkely makiing a mess.. 
if ( nnew[i]<0) nnew[i]=m.n[i];

pnew[i]=ni2/nnew[i];
if (var_region) { S+=nnew[i]; R+=pnew[i]; N+=m.Nteff[i];
Value qnet=nnew[i]-pnew[i]+m.Nteff[i]; rhotot+=qnet;
if ( qnet!=0) { if ( qnet>0) logrhotot+=::log(qnet);
else logrhotot-=::log(-qnet);
} 
}



} // i 

} //! attempte_log_update


if ( attempte_log_update)
{
for ( IdxTy i=0; i<(sz); ++i)
{
const bool var_region= ((i>=si)&&(i<sf)) ;
if ( var_region) nnew[i]=m.n[i]*::exp(1e-70*a*gradF[i]*m.n[i]); else nnew[i]=m.n[i];
pnew[i]=ni2/nnew[i];

if (var_region) { S+=nnew[i]; R+=pnew[i]; N+=m.Nteff[i];
Value qnet=nnew[i]-pnew[i]+m.Nteff[i]; rhotot+=qnet;
if ( qnet!=0) { if ( qnet>0) logrhotot+=::log(qnet);
else logrhotot-=::log(-qnet);
} 
}




} // i 


} //! attempte_log_update



if ( attempte_to_neutralize_log)
{
Value rhototnew=0;
//Value logrhototnew=0;
//Value logrhoxxx=logrhotot/(sf-si);
Value x1=-N/2.0/S; 
Value x2=.5/S*::sqrt(N*N+4.0*R*S);
Value fac=x1+x2; // this must be > 0
for ( IdxTy i=si; i<(sf); ++i)
{
//Value qnet=nnew[i]-pnew[i]+m.Nteff[i];
//if ( qnet>0) nnew[i]*=fac; else nnew[i]/=fac;
nnew[i]*=fac;
pnew[i]=ni2/nnew[i];
rhototnew+=nnew[i]-pnew[i]+m.Nteff[i];

}
rhotot=rhototnew;


}

// thisis more fudded up than when I started with defalt. 
// probably this needs to look at the sign and realize
// n varies many rders of magnitude.. 
if ( attempte_to_neutralize)
{
for ( IdxTy j=0; j<4; ++j)
{
Value rhototnew=0;
Value rhoxxx=rhotot/(sf-si);
Value okmin=rhoxxx;
// in this casae we can make n less than zero 
if ( rhoxxx>0) { 
for ( IdxTy i=si; i<(sf); ++i)
{

const Value nn2=nnew[i]*nnew[i];
if ( nn2==0) { okmin=0; break; } 
const Value dr=(nn2+ni2)/nn2;
const Value thismin=nnew[i]*dr;
if ( thismin<okmin) okmin=thismin;
//if ( okmin==0) okmin=thismin;
}
// ipostiive  amonts?
if ( okmin<rhoxxx) rhoxxx=okmin*.9;
} // positive

for ( IdxTy i=si; i<(sf); ++i)
{
const Value nn2=nnew[i]*nnew[i];
//nnew[i]-=rhoxxx*nn2/(nn2+ni2);
const Value dr=rhoxxx*nn2/(nn2+ni2);
// this is making a huge mess, just return false or something,. 
//if ( nnew[i]<0) { nnew[i]=def; } 
if ( nnew[i]<dr) { nnew[i]=def; } 
else {  nnew[i]-=dr; } 
pnew[i]=ni2/nnew[i];
rhototnew+=nnew[i]-pnew[i]+m.Nteff[i];
}
rhotot=rhototnew;
} // j 
} // neutralize


return rhotot;
}

// try to follow gradient to minimize F gnerally based on jdftx
// more careful analysis than last time but still missing things
// like actual energy level of carriers- just treat as rho time voltage
// regardless of conduction band level etc 
static void update_carriers_grad_2( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{

const bool dump_a_search=false;
const bool dump_at_dun=!false;
const bool dump_grad_mag=false;


const Value ni2=pr.ni2;
//const Value ni=pr.ni;
//const Value epx=1.0/pr.epsilon;
const Value kt=gp.vt();
//const Value sf=0; // gp.q()*kt;
// q is wrong? is this cgs units or ev? 
//const Value sf= gp.q()*kt*1e7; // entropy needs to be rechcked
const Value sf= gp.q()*kt; // wtf mks now... entropy needs to be rechcked
// I think this is MKS now? 
//const Value sf= gp.q()*kt; // entropy needs to be rechcked
//const Value sf= kt; // entropy needs to be rechcked
const Grid & g=m.n.grid();
// shoudl assert uniform
//const Value dx=g.min_dx();
//const Value dxf=dx*dx;
const Value sz=g.size();
Values nnew(g),pnew(g),x(g),y(g),gradF(g);
Values irhox(g),irho(g),gradrhos(g),gradrhoe(g),graden(g);
Values rhot=m.n+m.Nteff-m.p;
Told::I(irho,g,rhot);
Told::I(irhox,g,rhot*g);

static IntTy siter=-1;
++siter;

const Value siz=sf*entropy2(0,m,m.n,m.p,pr,gp,flags,junk);
const Value eiz=energy2(0,m,m.n,m.p,pr,gp,flags,junk);
const Value fiz=eiz-siz; 

// this should return grad_{net_count} or effective rho time q. 
entropy2(&gradrhos,m,m.n,m.p,pr,gp,flags,junk);
energy2(&gradrhoe,m,m.n,m.p,pr,gp,flags,junk);
// ideally I guess we could try grad log n ?
// this should remove mean in rho and change to grad_{n} 
fix_grad(graden,gradrhoe-gradrhos*sf,m,m.n,m.p,pr,gp,flags,junk,true,false);
//fix_grad(graden,gradrhoe-gradrhos*0,m,m.n,m.p,pr,gp,flags,junk);
// garden should be zero drho mean ( drho=0 ) with n now.
// F can be either linear or log in n???

const bool do_neutral=true; // (siter<500);
// rezero 
Value rhotot=update_from_grad(m.n,m.p,0,graden,sz,m,ni2,flags,do_neutral);
const Value si=sf*entropy2(0,m,m.n,m.p,pr,gp,flags,junk);
const Value ei=energy2(0,m,m.n,m.p,pr,gp,flags,junk);
const Value fi=ei-si; 

const Value mag_graden=::sqrt(graden.sum_squares());
if ( dump_grad_mag)
{
// metrics
const Value mag_gradrhoe=::sqrt(gradrhoe.sum_squares());
const Value mag_gradrhos=sf*::sqrt(gradrhos.sum_squares());

MM_MSG(" mags graden= "<<mag_graden<<" gradrhoe= "<<mag_gradrhoe<<" gradrhos= "<< mag_gradrhos<<" rezed= "<<(fi-fiz)<<" ezed= "<<(ei-eiz)<<" szed= "<<(si-siz)<<" rhototx= "<<rhotot)
}


bool dun=false;
Value iter=0;
// if the energies change scale, this needs to change too
// although perhaps scale by norm of grad or something?
// added 1e-10 after attempting to get to mks, worked for while
// then suddently exploded with higher a. 
//Value a=-1e-5*1e-10;
const Value ff=::pow(1.2,600);
Value a=-1e-5*1e-14*ff;
if ( mag_graden>1) a=mag_graden*1e-21*ff;
//const Value afac=1.2;
const Value afac=1.2;
//Value fbest=fi;
Value delfmax=0;
const IdxTy iter_max=6000;
while (iter<iter_max)
{ // the less than zero carrier con can not really be fixed,
// that should set a max a... 
Value rhotot=update_from_grad(nnew,pnew,a,graden,sz,m,ni2,flags,do_neutral);
const Value sa=sf*entropy2(0,m,nnew,pnew,pr,gp,flags,junk);
const Value ea=energy2(0,m,nnew,pnew,pr,gp,flags,junk);
const Value fa=ea-sa; 
const Value delF=fa-fi;
if ( dump_a_search)
{
MM_MSG(siter<<" "<<iter<<" a= "<<a<<" fi= "<<fi<<" ei= "<<ei<<" si= "<<si<<" fa= "<<fa<<" ea= "<<ea<<" sa= "<<sa<<" df= "<<(fa-fi)<<" ds = "<<(sa-si)<<" de= "<<(ea-ei)<<" rhotot= "<<rhotot)
} // dump

if ( dun ) break;
if ( delfmax!=0 ){  if ( delF>delfmax)  { a=a/afac; dun=true; }} 
if ( delF<0){   if ( delF<delfmax) delfmax=delF; }
if ( !dun)  a=a*afac;
++iter;
if ( delF<0) if ( iter==iter_max) dun=true; 
if ( dun&&dump_at_dun)
{
MM_MSG("dump_at_dun "<< siter<<" "<<iter<<" a= "<<a<<" fi= "<<fi<<" ei= "<<ei<<" si= "<<si<<" fa= "<<fa<<" ea= "<<ea<<" sa= "<<sa<<" df= "<<(fa-fi)<<" ds = "<<(sa-si)<<" de= "<<(ea-ei)<<" rhotot= "<<rhotot)

}


} // iter

if ( !dun)
{

MM_MSG(" failed to improve so break ")
junk.terminate_code=1; 
return;

}

// this destroys the pn=ni2 thing.. 
const Value fold=0; // flags.nfold;
const IdxTy szi=flags.fixed_start; // 3;
const IdxTy szf=sz-flags.fixed_end; // 2;
//m.n.mix(nnew,fold,2,sz-2);
m.n.mix(nnew,fold,szi,szf);
//m.p.mix(pnew,fold,2,sz-2);
m.p.mix(pnew,fold,szi,szf);

}



// try to follow gradient to minimize F gnerally based on jdftx
static void update_carriers_grad( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Value ni2=pr.ni2;
const Value ni=pr.ni;
const Value epx=1.0/pr.epsilon;
const Value kt=gp.vt();
const Value sf=gp.q()*kt;
const Grid & g=m.n.grid();
// shoudl assert uniform
//const Value dx=g.min_dx();
//const Value dxf=dx*dx;
const Value sz=g.size();
Values nnew(g),pnew(g),x(g),y(g),gradF(g);
Values irhox(g),irho(g);
Values rhot=m.n+m.Nteff-m.p;
Told::I(irho,g,rhot);
Told::I(irhox,g,rhot*g);

//const Value c=gp.q()/(gp.vt()*pr.epsilon*gp.epsilon_vacu());
static Value siter=0;
const Value si=sf*entropy(m,m.n,m.p,pr,gp,flags,junk);
const Value ei=energy(m,m.n,m.p,pr,gp,flags,junk);
const Value fi=ei-si; 
MM_MSG(siter<<" ini e= "<<ei<<" si= "<<si<<" F= "<<(fi))
Value dN=0;
Value fmag=0;
//Value sni2=::log(ni2);
for ( IdxTy i=0; i<(sz); ++i)
{
Value n2=m.n[i]*m.n[i];
// in theory anyway n2 should not be zero...
if (n2==0 ) n2=1.0/ni2;
const Value drhodn=(n2-ni2)/n2;
// this hsould be the rho v term not D dot E 
const Value drhov= drhodn*(-0*m.v[i]+ gp.q()*epx*(irhox[sz-1]-irhox[i]-g[i]*(irho[sz-1]-irho[i])));
const Value edote=drhodn*2.0*(m.v[sz-1]-m.v[i]);

//const Value dn=(drhov+edote)-kt*(::log(m.n[i]));
const Value pp=ni2/m.n[i];
const Value dn=(drhov+edote)-kt*(::log(m.n[i])-(pp/m.n[i])*::log(pp));
gradF[i]=dn;
dN+=dn;
fmag+=dn*dn;

} // i 
fmag=::sqrt(fmag);

Value aiter=0;
// should project out changes in N?
Value amin=-.001;
Value amax=-1e16;

Value a=amin;
Value delfmax=0;
bool dun=false;
// yes, this makes sense
while (a>amax)
{
// do a search through "a" and try to find mini F
// and verify it makes sense...
Value rhotot=0;
for ( IdxTy i=0; i<(sz); ++i)
{
nnew[i]=m.n[i]+a*gradF[i];
if ( nnew[i]<0) nnew[i]=::sqrt(ni);
pnew[i]=ni2/nnew[i];
rhotot+=nnew[i]-pnew[i]+m.Nteff[i];
}


Value sx=sf*entropy(m,nnew,pnew,pr,gp,flags,junk);
Value ex=energy(m,nnew,pnew,pr,gp,flags,junk);
const Value fx=ex-sx; 
const Value delF=fx-fi; 
MM_MSG(" "<<dun<<" "<<siter<<" "<<aiter<<" iter ini a= "<<a<< " e= "<<ex<<" si= "<<sx<<" delF= "<<(delF)<<" delE= "<<(ex-ei)<<" delS= "<<(sx-si)<<" rtot= "<<rhotot<<" dN= "<<dN<<" fmag= "<<fmag)
++aiter;
if ( dun ) break;
if ( delfmax!=0 ){  if ( delF>delfmax)  { a=a/2; dun=true; }} 
if ( delF<0){   if ( delF<delfmax) delfmax=delF; }
if ( !dun)  a=a*2;
} // while a 

if ( !dun)
{

MM_MSG(" failed to improve so break ")
junk.terminate_code=1; 
return;

}

const Value fold=flags.nfold;
m.n.mix(nnew,fold,2,sz-2);
m.p.mix(pnew,fold,2,sz-2);

}

typedef double BlasTy;

class blas_diags
{

typedef BlasTy Dd;
public:

blas_diags( const IdxTy szn_, const IdxTy KL_, const IdxTy KU_)
{
szn=szn_;
KL=KL_;
KU=KU_;
szwork=2*KL+KU+1;
szdiag=szwork*szn;
rhs =new BlasTy[szn];
// thjis is a fudding transpose.. 
diags  =new BlasTy[szdiag];
diagsfudd=0;
divisors=0;
// only for diagnostics right now 
//::memset(diags,0,sizeof(BlasTy)*szdiag);
for ( IdxTy i=0; i<IdxTy(szdiag); ++i) diags[i]=0;
for ( IdxTy i=0; i<IdxTy(szn); ++i) rhs[i]=0;
set_pointers ();




}
void set_pointers()
{
diag=diags+(KL+KU)*szn;
// the first n elements not used on uppers, lowers leave tailsjj
if ( KU>=1) u1diag=diags+(KL+KU-1)*szn+1; else u1diag=0;
if ( KU>=2) u2diag=diags+(KL+KU-2)*szn+2; else u2diag=0;
if ( KU>=3) u3diag=diags+(KL+KU-3)*szn+3; else u3diag=0;
//u3diag=0;
if ( KL>=1) d1diag=diags+(KL+KU+1)*szn; else d1diag=0;
if ( KL>=2) d2diag=diags+(KL+KU+2)*szn; else d2diag=0;
if ( KL>=3) d3diag=diags+(KL+KU+3)*szn; else d3diag=0;

}

// the fixed names for diagonals rather than subscripts 
// are making a mess but ok for now. 
void add_upper(Dd * newdiag)
{
IdxTy KUnew=KU+1;
IdxTy szworknew=2*KL+KUnew+1;
IdxTy  szdiagnew=szworknew*szn;
Dd* diagsnew  =new BlasTy[szdiagnew];
const IdxTy udate=szn;
for ( IdxTy i=0; i<IdxTy(szdiag); ++i) diagsnew[i+udate]=diags[i];

szdiag=szdiagnew;
szwork=szworknew;
KU=KUnew;
delete [] diags;
diags=diagsnew;
set_pointers();

u3diag=diags+(KL+KU-3)*szn+3;
for ( IdxTy i=0; i<(szn-3); ++i) u3diag[i]=newdiag[i];

}


~blas_diags() { delete [] rhs; delete [] diags;
delete[] diagsfudd; delete [] divisors;  } 

Dd *  reform() 
{
 diagsfudd  =new BlasTy[szdiag];
for ( IdxTy ifudd=0; ifudd<IdxTy(szwork); ++ifudd)
{
const IdxTy base=szn*ifudd;
for ( IdxTy jfudd=0; jfudd<IdxTy(szn); ++jfudd)
{
const IdxTy fuk=jfudd*szwork+ifudd;
diagsfudd[fuk]=diags[base+jfudd];
// MM_MSG(ifudd<<" "<<jfudd<<" "<<diagsfudd[fuk])
}
}
return diagsfudd;
} // reform
// precondition here
void assign(Dd & d, const Dd * x) { if (x!=0) d=*x; } 
void assign(Dd * d, const Dd & x) { if (d!=0) *d=x; } 
void precondition(const IdxTy starter, const IdxTy end2)
{
//if ( true ) return; 
// this seems important now for getting abrupt right 
const bool  add_and_sub_n_and_p=!false;

divisors =new BlasTy[szn];
// the boundary conditions need to be fixed lol 
for(IdxTy i=0; i<IdxTy(starter); ++i){ diag[i]=1; rhs[i]=0;}
for(IdxTy i=IdxTy(end2); i<IdxTy(szn); ++i){ diag[i]=1; rhs[i]=0;}
// no longer part of solve
//for ( IdxTy eqn=0; eqn<IdxTy(szn); ++eqn)
for ( IdxTy eqn=starter; eqn<end2; ++eqn)
{
divisors[eqn]=* brute_force(eqn,0); // rhs[eqn];
if ( divisors[eqn]==0) continue; 
//rhs[eqn]=1;
rhs[eqn]/=divisors[eqn];
for ( int di=-KL; di<=KU; ++di)
{
Dd * idx=brute_force(eqn,di);
if ( idx!=0 ) (*idx)/=divisors[eqn];
} // di 
} // eqn
if ( add_and_sub_n_and_p)
{
// also note that these alternate as a-0-(-2a)-0-a and
// a-a or a-(-a)  at least with huge GR term . 
const IdxTy coefs=KU+KL+1;
Dd t1[coefs],t2[coefs];
for ( IdxTy eqn=0; eqn<IdxTy(szn); eqn+=2)
{
// the lowest diagonal for "y" is zero, adding it to the "x"
// equation beneath does nothing and similarly for 
// the "p" equation adding back to the above eqn
for ( int di=-KL; di<=KU; ++di) 
{assign(t1[KL+di],brute_force(eqn,di)); 
assign(t2[KL+di],brute_force(eqn+1,di));  } 


//assign(brute_force(eqn+1,KU),t1[KU]+t2[KU-1]); 
for ( int di=-(KL-1); di<=(KU-1); ++di) 
{assign(brute_force(eqn,di),t1[KL+di]+t2[KL+di-1]); 
assign(brute_force(eqn+1,di),t1[KL+di+1]-t2[KL+di]); 
} // assign back 

const Dd r1=rhs[eqn];
const Dd r2=rhs[eqn+1];
rhs[eqn]=r1+r2;
rhs[eqn+1]=r1-r2;
// MM_MSG(" newy "<<dump_eqn(eqn)<<" rhs "<<rhs[eqn])
// MM_MSG(" newx "<<dump_eqn(eqn+1)<<" rhs "<<rhs[eqn+1])

} // eqn
 } // add and subtract n and p 

// fix the boundaries agagin...

}

// after solving, rhs should have solution, now fix.. 
// at this point, the rhs needs to bekkkkkkk
void decondition()
{



}
// note d is SIGNED wih main diag as zed
// actually eqn is unsigned 
Dd * brute_force(const IntTy eqn, const IntTy d)
{
switch (d)
{
case -3 : { if ( d3diag!=0) if ( eqn>=-d) return &d3diag[eqn+d]; }
case -2 : { if ( d2diag!=0) if ( eqn>=-d) return &d2diag[eqn+d]; }
case -1 : { if ( d1diag!=0) if ( eqn>=-d) return &d1diag[eqn+d]; }
case 0 : { return &diag[eqn]; }
case 1 : { if ( u1diag!=0) return &u1diag[eqn]; }
case 2 : { if ( u2diag!=0) return &u2diag[eqn]; }
case 3 : { if ( u3diag!=0) return &u3diag[eqn]; }
default : { MM_MSG(" blas_diag default idx "<<eqn<<" "<<d) } 

}; // swich

return 0; 
}

StrTy dump_eqn(const IdxTy jfudd) const
{
SsTy ss;
if ( jfudd>2) if ( d3diag!=0) ss<<" "<<d3diag[jfudd-3];
if ( jfudd>1) ss<<" "<<d2diag[jfudd-2];
if ( jfudd>0) ss<<" "<<d1diag[jfudd-1];
ss<<" "<<diag[jfudd];
if ( jfudd<(szn-1)) ss<<" "<<u1diag[jfudd];
if ( jfudd<(szn-2)) ss<<" "<<u2diag[jfudd];
if ( jfudd<(szn-3)) if ( u3diag!=0) ss<<" "<<u3diag[jfudd];

//MM_MSG(siter<<" "<<jfudd<<" eq "<<ss.str()<<" rhs "<<rhs[jfudd]) 

return ss.str();
}

void dump_msg(const IdxTy siter)
{

 for ( int jfudd=0; jfudd<szn; ++jfudd) 
{  
SsTy ss;
if ( jfudd>1) ss<<d2diag[jfudd-2];
if ( jfudd>0) ss<<" "<<d1diag[jfudd-1];
ss<<" "<<diag[jfudd];
if ( jfudd<(szn-1)) ss<<" "<<u1diag[jfudd];
if ( jfudd<(szn-2)) ss<<" "<<u2diag[jfudd];
if ( u3diag!=0) if ( jfudd<(szn-3)) ss<<" "<<u3diag[jfudd];

MM_MSG(siter<<" "<<jfudd<<" eq "<<ss.str()<<" rhs "<<rhs[jfudd]) 

} 
} // msg
/////////////////////////////////////////
void solve(const bool dump_solution,const IdxTy siter)
{
double * rhsx= 0;
// let user call if needed
//precondition();
BlasTy * diagsfudd  =reform();
int NRHS=1;
int INFO=99999;
int N = szn; 
if ( dump_solution)
{rhsx= new double[szn];
::memcpy(rhsx,rhs,szn*sizeof(double));
}
// pivot does not appear tobe needed to sort out results lol
int * ipiv = new int[szn] ; 
 dgbsv_(&N,&KL, & KU, &NRHS, diagsfudd, &szwork,ipiv, rhs, &szn ,&INFO);
if (INFO!=0) { { MM_MSG(siter<<" DGBSV error code is "<<INFO) } 
if ( INFO>0) {  MM_MSG(" info egn "<<dump_eqn(INFO-1)) } exit(-1); }
delete [] ipiv; 
// we need to know what, if anything to do here
// or just let user call 
//decondition();
if ( dump_solution)
{ for (IdxTy i=0; i<IdxTy(szn); ++i) { 
MM_MSG("solutx "<<siter<<" "<<i<<dump_eqn(i)<<" rhx "<<rhsx[i]<<" = "<<rhs[i])
} //i 
delete [] rhsx;
}
} // solve

////////////////////////////////////////////
int KL,KU,szn;
int  szwork,szdiag;
Dd * rhs;
Dd * diags;
Dd * diagsfudd;
Dd * diag;
Dd * u1diag;
Dd * u2diag;
Dd * u3diag;
Dd * d1diag;
Dd * d2diag;
Dd * d3diag;

Dd * newdiag;
// first shot at conditioning...
Dd * divisors;

}; // blas_diags

// I think that adding the e-field derivative can be reduced
// to banded form by exploiting the product form of the e-field
// elements. Have to verify this however, esp with n and p both
// non-band elements should be q/eps*(n,p)[i]*grady[j] and gauss reduce
// out to leave banded form 
// consider the ones below diag positive, those above negative.
// "e" goes across, that makes grady the "f" term 

// eliminate alt equations which should be same carriers
// hopefully more stable jjjjjjj
static void reduce_efield_car (  blas_diags & bd,  int & KU,
const double *  grade, const double *  grady, const IdxTy  modifiers)
{
const bool assert_finite=true;
static IdxTy  siter=0;
//const Grid & g=grade.grid();
//const IdxTy sz=g.size();
//if ( sz==0) return; 

BlasTy * rhs=bd.rhs;
BlasTy * diag=bd.diag;
BlasTy * u1diag=bd.u1diag;
BlasTy * u2diag=bd.u2diag;
BlasTy * d1diag=bd.d1diag;
BlasTy * d2diag=bd.d2diag;

const IdxTy szn=bd.szn;
BlasTy * newdiag=new BlasTy[szn]; // bd.newdiag;
for ( IdxTy i=0; i<szn; ++i) newdiag[i]=0;
//MM_MSG(" test xxxxx " ) 

// do row i and i-1, last eqn needs to be fixed again... 
for ( IdxTy i=2; i<(szn-2); ++i)
{
// multiply by ratio of "f" or row indexed factor and subtract
// to zero out everything except things crossing the diags. 
// if f is zero, then what? 
// the abs carrier concentrations altnerate but these are logs
// do this for all of them until substracdted out?
diag[i]+=-grady[i]*grade[i];
if ( i>0) u1diag[i]+=-grady[i]*grade[i+1];
if( grady[i+1]==0){
MM_MSG(" facfac "<<siter<<" "<<i<<" "<<"xxx"<<" "<<grady[i]<<" "<<grady[i+1])
 continue;
}

const Value fac=-grady[i]/grady[i+1];
MM_MSG(" facfac "<<siter<<" "<<i<<" "<<fac<<" "<<grady[i]<<" "<<grady[i+1])

const Value limx=::sqrt(fac*fac);
if ( limx>100) continue;
if ( limx<.001) continue;

// eh, this was stupid and not right we only need to do the DIAGS
// not all the zero elements lol 
// rightnow the diags are transposed from the way blas ultiamtely needs them,
// this would be easier to do later... 
//const IdxTy idx=i+j*sz+othercrap; 
//const IdxTy idx1=i+j*sz+othercrap; 
//diag[idx]+=fac*diags[idx1];
const IdxTy ii=i;
// the off diag subscripting is wierd, should put this in object
// none of the efield terms have been included yet, put them here
//diag[ii]+=fac*d1diag[ii]-grady[i]*grade[ii];
diag[ii]+=fac*d1diag[ii];
d1diag[ii-1]+=fac*d2diag[ii-1];
// diags are zed initially, so need to add the field term . 
//if ( ii>0) u1diag[ii]+=fac*diag[ii+1]+grady[i]*grade[ii+1];
if ( ii>0) u1diag[ii]+=fac*diag[ii+1]-grady[i]*grade[i+1];;
if ( ii>0) u2diag[ii]+=fac*u1diag[ii+1];
if ( ii>0) newdiag[ii]+=fac*u2diag[ii+1];
// now do rhs too
rhs[i]+=fac*rhs[i+1];
if ( assert_finite)
{
bool isok=isfinite(diag[ii])&&isfinite(d1diag[ii-1])&&isfinite(u1diag[ii])
&&isfinite(u2diag[ii])&&isfinite(newdiag[ii]);
if (!isok)
{
MM_MSG( " notfinite "<<ii<<" "<< fac<<" "<<d1diag[ii]<<" "<<d2diag[ii-1]<<" "<<
diag[ii+1]<<" u1= "<<u1diag[ii+1]<< " u2= "<<u2diag[ii+1]<<" "<<grady[i]<<" "<<grade[i+1]<<" new "<<newdiag[ii]<<" "<<u2diag[ii+1])

}
} // assert_finite 
} // i 
// need to fix first eqn too, shold have saved some crpa first. 
bd.add_upper(newdiag);
KU+=2;
delete [] newdiag; 
++siter; 
}


// see comments in eliminate routine that calls this...
// this eliminates adjacent equations that use different carriers
static void reduce_efield (  blas_diags & bd,  int & KU,
const double *  grade, const double *  grady, const IdxTy  modifiers)
{
const bool assert_finite=true;
static IdxTy  siter=0;
//const Grid & g=grade.grid();
//const IdxTy sz=g.size();
//if ( sz==0) return; 

BlasTy * rhs=bd.rhs;
BlasTy * diag=bd.diag;
BlasTy * u1diag=bd.u1diag;
BlasTy * u2diag=bd.u2diag;
BlasTy * d1diag=bd.d1diag;
BlasTy * d2diag=bd.d2diag;

const IdxTy szn=bd.szn;
BlasTy * newdiag=new BlasTy[szn]; // bd.newdiag;
for ( IdxTy i=0; i<szn; ++i) newdiag[i]=0;
//MM_MSG(" test xxxxx " ) 

// do row i and i-1, last eqn needs to be fixed again... 
for ( IdxTy i=10; i<(szn-10); ++i)
{
// multiply by ratio of "f" or row indexed factor and subtract
// to zero out everything except things crossing the diags. 
// if f is zero, then what? 
// the abs carrier concentrations altnerate but these are logs
// do this for all of them until substracdted out?
diag[i]+=-grady[i]*grade[i];
if ( i>0) u1diag[i]+=-grady[i]*grade[i+1];
 if( grady[i+1]==0){
if ( verbose_crap) 
{MM_MSG(" facfac "<<siter<<" "<<i<<" "<<"xxx"<<" "<<grady[i]<<" "<<grady[i+1]) } 
// continue;
}
const bool zedd=(grady[i+1]==0);
// assume that adjacent grads are opposite or -1 to maintain np. 
Value fac=zedd?(1):(-grady[i]/grady[i+1]);
if ( verbose_crap) 
{ MM_MSG(" facfac "<<siter<<" "<<i<<" "<<fac<<" "<<grady[i]<<" "<<grady[i+1]) } 

const Value limx=::sqrt(fac*fac);
//if ( limx>10) continue;
if ( limx>10) { fac=(fac<0)?(-10):(10);  }
//if ( limx<.1) continue;
if ( limx<.1) { fac=(fac<0)?(-.1):(.1);  } 

// eh, this was stupid and not right we only need to do the DIAGS
// not all the zero elements lol 
// rightnow the diags are transposed from the way blas ultiamtely needs them,
// this would be easier to do later... 
//const IdxTy idx=i+j*sz+othercrap; 
//const IdxTy idx1=i+j*sz+othercrap; 
//diag[idx]+=fac*diags[idx1];
const IdxTy ii=i;
// the off diag subscripting is wierd, should put this in object
// none of the efield terms have been included yet, put them here
//diag[ii]+=fac*d1diag[ii]-grady[i]*grade[ii];
diag[ii]+=fac*d1diag[ii];
d1diag[ii-1]+=fac*d2diag[ii-1];
// diags are zed initially, so need to add the field term . 
//if ( ii>0) u1diag[ii]+=fac*diag[ii+1]+grady[i]*grade[ii+1];
if ( ii>0) u1diag[ii]+=fac*diag[ii+1]-grady[i]*grade[i+1];;
if ( ii>0) u2diag[ii]+=fac*u1diag[ii+1];
if ( ii>0) newdiag[ii]+=fac*u2diag[ii+1];
// now do rhs too
rhs[i]+=fac*rhs[i+1];
if ( assert_finite)
{
bool isok=isfinite(diag[ii])&&isfinite(d1diag[ii-1])&&isfinite(u1diag[ii])
&&isfinite(u2diag[ii])&&isfinite(newdiag[ii]);
if (!isok)
{
MM_MSG( " notfinite "<<ii<<" "<< fac<<" "<<d1diag[ii]<<" "<<d2diag[ii-1]<<" "<<
diag[ii+1]<<" u1= "<<u1diag[ii+1]<< " u2= "<<u2diag[ii+1]<<" "<<grady[i]<<" "<<grade[i+1]<<" new "<<newdiag[ii]<<" "<<u2diag[ii+1])

}
} // assert_finite 
} // i 
// need to fix first eqn too, shold have saved some crpa first. 
bd.add_upper(newdiag);
++KU;
delete [] newdiag; 
++siter; 
}







class GRFunction
{
typedef double D;
// the default is for gr/n and gr/p to match
// the current stupid needs. jj
public:
GRFunction() : m_g(0),m_n(0),m_p(0),m_ni(0),m_ni2(0),
m_gr(0),m_grx(0),m_gry(0),m_dyy(0),m_dxx(0),m_dyx(0),m_dxy(0) {}

GRFunction(const D & g, const D & ni) : m_g(g),m_n(0),m_p(0),
m_ni(ni),m_ni2(ni*ni),
m_gr(0),m_grx(0),m_gry(0),m_dyy(0),m_dxx(0),m_dyx(0),m_dxy(0) {}


// for wierd things like field depenednece, add more sigs etc
//template <IdxTy flags> 
void set( const D & n, const D & p, const IdxTy flags )
{ m_n=n; m_p=p;

// user should divide by x or y
const D d=m_n+m_p+2.0*m_ni;
m_gr=m_g*(m_ni2-m_n*m_p)/(d);

if ( m_p!=0) m_grx=m_gr/m_p; else m_grx=0;
if ( m_n!=0) m_gry=m_gr/m_n; else m_gry=0;

//  doh, the errors at each location in clude p and n errors
// and both have x and y gr derivatices,,,
// first numerator then denominator 
if ( m_p!=0) m_dxx=-m_g*m_ni2/m_p/d-m_gr/d;
else m_dxx=0;
// the y derivative of gr div x, or exp(ln(n))
if ( m_n!=0) m_dxy=-m_g*m_n/d-m_grx*m_n/d;
else m_dxy=0;
// the ders of the der y thing, or p term 
if ( m_n!=0) m_dyy=-m_g*m_ni2/m_n/d-m_gr/d;
else m_dyy=0;

if ( m_p!=0) m_dyx=-m_g*m_p/d-m_gry*m_p/d;
else m_dyx=0; 
}
// for the newer version, need another method... 
// the gr term has been multipled bu more, for now have const 
// denominateor too... 
void set_alt( const D & n, const D & p, const IdxTy flags )
{
 m_n=n; m_p=p;
const Value d=1;
m_gr=m_g*(m_ni2-m_n*m_p)/(d);
m_grx=m_g*(m_ni2+m_ni2/m_p*m_n-m_n*m_n-m_n*m_p);
m_gry=m_g*(m_ni2+m_ni2/m_n*m_p-m_p*m_p-m_n*m_p);
// dAy/dy
m_dyy=m_g*(-m_ni2/m_n*m_p-m_n*m_p);
m_dyx=m_g*(m_ni2/m_n*m_p-2.0*m_p*m_p-m_n*m_p);
// dAx/x
m_dxx=m_g*(-m_ni2/m_p*m_n-m_n*m_p);
m_dxy=m_g*(m_ni2/m_p*m_n-2.0*m_n*m_n-m_n*m_p);

} // set_alt
// the convention changed, these are all used as POSITVIE
// byu caller.... 
void set_alt2( const D & n, const D & p, const IdxTy flags )
{
 m_n=n; m_p=p;
const Value d=1;
m_gr=m_g*(m_ni2-m_n*m_p)/(d);
m_grx=m_g*(m_ni2/m_p-m_n);
m_gry=-m_g*(m_ni2/m_n-m_p);
// dAy/dy
m_dyy=m_g*(m_ni2/m_n);
m_dyx=m_g*(m_p);
// dAx/x
m_dxx=-m_g*(m_ni2/m_p);
m_dxy=-m_g*(m_n);

} // set_alt



// private:
D m_g,m_n,m_p,m_ni,m_ni2;

// gr is actual gr rate for given n and p, 
// dy and dx are derivaives wrt n=exp(y) and 
// p=exp(x) 
D m_gr, m_grx, m_gry, m_dyy,m_dxx,m_dyx,m_dxy;

}; // GRFunction  
static void make_old_crap(blas_diags & bd, const BlasTy * rhs,  const Values & nx, const Values & px,
const Values & Nteff,const Values & y, const Values & x, const Values & e, GRFunction & grfunc
, const Value & dxi2,const Value & vtd2, const Value & qe,
const Value & c,const IdxTy & sz, const bool assert_finite
)
{
static IdxTy siter=0;
const bool  efield = true; 
const bool  make_rhs_dy_dt =true;

const Value vtd22=vtd2*2.0;
const IdxTy szn=sz*2;
double * & diag=bd.diag;
double * & u1diag=bd.u1diag;
double * & u2diag=bd.u2diag;
double * & d1diag=bd.d1diag;
double * & d2diag=bd.d2diag;




IdxTy j=2;
for (IdxTy i=1; i<(sz-1); ++i)
{
// 1 flag should calc all ders
grfunc.set(nx[i],px[i],1);

// this can not work for the ( make_rhs_dy_dt ) case yet  
const Value xterm= 0; //    np*(qe-fi);
const Value eid=e[i]*dxi2;
const Value dy=.5*(y[i+1]-y[i-1]);

if ( make_rhs_dy_dt )
{
// the gr term is surely not right yet. 
// grn was needed to remove the "G" from the rhs term, so not needed
// here 
//diag[j]= (-vtd22-qe*nx[i])+grn- gr_part2;
//const Value dg=fi*(ni2/(nx[i]*nx[i])+gr);
const Value dg=grfunc.m_dyy; // -gr_part2+ fi*(ni2/(nx[i]));
diag[j]= (-vtd22-qe*nx[i])+dg;
MM_MSG(siter<< " diagterm "<<i<<" "<<vtd22<<" "<<(qe*nx[i])<<" "<<dg<<" "<<grfunc.m_gry<<" "<<((nx[i]*px[i])/grfunc.m_ni2)<<" "<<grfunc.m_gr<<" rhs "<<rhs[j])

d2diag[j-2]=(-eid+vtd2*(1.0-dy));
u2diag[j]=(eid+vtd2*(1.0+dy));
d1diag[j-1]=0;
// for the n term, this is the p derivative 
u1diag[j]=grfunc.m_dyx; // xterm/nx[i];
}
else{
// this needs to be fixed, just to get it to compile
diag[j]= nx[i]*(-vtd22-qe*nx[i])+0; // grn+rhs[j]- gr_part2*nx[i];
d2diag[j-2]=nx[i]*(-eid+vtd2*(1.0-dy));
u2diag[j]=nx[i]*(eid+vtd2*(1.0+dy));
d1diag[j-1]=0;
// for the n term, this is the p derivative 
u1diag[j]=xterm;


}//   ( make_rhs_dy_dt )


// inclue the e field terms that are on the digaonals, the
// others are ignored although this could go either way
// when neutral. The charge total must be zero..
// as field is zero at both sides 
if ( efield&& make_rhs_dy_dt )
{ // n terms are negative, 
const Value fudd=qe*dy; // integral and der deltas cancner 
d2diag[j-2]+=(-fudd);
d1diag[j-1]+=fudd;
// I guess this is debatable, 1/2 maybe?
// probably exactly zed, 
//diag[j]+=.5*(-fudd)*nx[i];
u1diag[j]-=fudd;
u2diag[j]-=(-fudd);
}
else if ( efield&& !make_rhs_dy_dt )
{
const Value fudd=qe*dy; // integral and der deltas cancner 
d2diag[j-2]+=(-fudd)*nx[i-1];
d1diag[j-1]+=fudd*px[i-1];
// I guess this is debatable, 1/2 maybe?
// probably exactly zed, 
diag[j]+=.5*(-fudd)*nx[i];
//u1diag[j]-=fudd*px[i];

}

// now for p part
++j;

const Value dyx=.5*(x[i+1]-x[i-1]);

// gr not right yet 
if ( make_rhs_dy_dt )
{
// need to fix gr term yet... 
//const Value dg=fi*(ni2/(px[i]*px[i])+gr);
const Value dg=grfunc.m_dxx; // -gr_part2+ fi*(ni2/(px[i]));
//diag[j]= (vtd22+qe*px[i])-grn+gr_part2;
diag[j]= (vtd22+qe*px[i])-dg;
MM_MSG(siter<< " diagtermp "<<i<<" "<<vtd22<<" "<<(qe*px[i])<<" "<<dg<<" "<<grfunc.m_grx<<" "<<((nx[i]*px[i])/grfunc.m_ni2)<<" "<<grfunc.m_gr<<" rhs "<<rhs[j])
// this needs a check if loop runs to limits 
// this is the -1, which is now off minus 2
d2diag[j-2]=(-eid-vtd2*(1.0-dyx));
u2diag[j]=(eid-vtd2*(1.0+dyx));
d1diag[j-1]=-grfunc.m_dxy; // -xterm/px[i];
//d1diag[j-1]=   - np*(qe+fi);
u1diag[j]=0;
}
else
{
diag[j]= px[i]*(vtd22+qe*px[i])-0; // grn+rhs[j]+gr_part2*px[i];
d2diag[j-2]=px[i]*(-eid-vtd2*(1.0-dyx));
u2diag[j]=px[i]*(eid-vtd2*(1.0+dyx));
d1diag[j-1]=-xterm;
//d1diag[j-1]=   - np*(qe+fi);
u1diag[j]=0;

}//  ( make_rhs_dy_dt )


if ( efield&& make_rhs_dy_dt )
{
const Value fudd=-qe*dyx; // integral and der deltas cancner 
d2diag[j-2]+=(-fudd);
d1diag[j-1]+=fudd;
//diag[j]+=.5*(-fudd)*px[i];
u1diag[j]-=fudd;
u2diag[j]-=(-fudd);


}
else if ( efield&& !make_rhs_dy_dt )
{ // n terms are negative, 
const Value fudd=-qe*dyx; // integral and der deltas cancner 
d2diag[j-2]+=(-fudd)*px[i-1];
d1diag[j-1]+=fudd*nx[i];
diag[j]+=.5*(-fudd)*px[i];
//u1diag[j]-=fudd*nx[i+1];
//u2diag[j]-=(-fudd)*px[i+1];
}
++j;
} // i diagonals...
if ( true)
{
// this should be similar to the next diags for numerical issues
// although I guess it should decouple anyway.
diag[0]=diag[2]; diag[1]=diag[3];
u1diag[0]=0;
u2diag[0]=0;
u1diag[1]=0;
u2diag[1]=0;
// wtf? 
d1diag[0]=0;
// now lets see if the others work lol
const IdxTy endfudd=szn-1;
diag[endfudd]=diag[endfudd-2];
diag[endfudd-1]=diag[endfudd-3];
d1diag[endfudd-1]=0;
d1diag[endfudd-2]=0;
u1diag[endfudd-1]=0;
d2diag[endfudd-2]=0;
d2diag[endfudd-3]=0;


} // boundary 1


} // make_old_crap



static void make_bd(blas_diags & bd, const BlasTy * rhs,  const Values & nx, const Values & px,
const Values & Nteff,const Values & y, const Values & x, const Values & e, GRFunction & grfunc
, const Value & dxi2,const Value & vtd2, const Value & qe,
const Value & c,const IdxTy & sz, const bool assert_finite
)
{
static IdxTy siter=0;
const Value vtd22=vtd2*2.0;
const IdxTy szn=sz*2;
double * & diag=bd.diag;
double * & u1diag=bd.u1diag;
double * & u2diag=bd.u2diag;
double * & d1diag=bd.d1diag;
double * & d2diag=bd.d2diag;


IdxTy j=2;
for (IdxTy i=1; i<(sz-1); ++i)
{
// 1 flag should calc all ders
grfunc.set(nx[i],px[i],1);
const Value eid=e[i]*dxi2;
{
const Value dy=.5*(y[i+1]-y[i-1]);

const Value dg=grfunc.m_dyy; 
diag[j]= (-vtd22-qe*nx[i])+dg;
if ( verbose_crap){
MM_MSG(siter<< " diagterm "<<i<<" "<<eid<<" "<<dy<<" "<<vtd22<<" "<<(qe*nx[i])<<" getc "<<dg<<" "<<grfunc.m_gry<<" "<<((nx[i]*px[i])/grfunc.m_ni2)<<" "<<grfunc.m_gr<<" rhs "<<rhs[j])
}

d2diag[j-2]=(-eid+vtd2*(1.0-dy));
u2diag[j]=(eid+vtd2*(1.0+dy));
d1diag[j-1]=0;
// for the n term, this is the p derivative 
u1diag[j]=grfunc.m_dyx +qe*px[i]; // xterm/nx[i];

} // n  scoping 
// p part 
++j; 
const Value dyx=.5*(x[i+1]-x[i-1]);
// blocks to scope variables
{
const Value dg=grfunc.m_dxx; // -gr_part2+ fi*(ni2/(px[i]));
//diag[j]= (vtd22+qe*px[i])-grn+gr_part2;
diag[j]= (vtd22+qe*px[i])-dg;
if ( verbose_crap){
MM_MSG(siter<< " diagtermp "<<i<<" "<<eid<<" "<<dyx<<" "<<vtd22<<" "<<(qe*px[i])<<" getc "<<dg<<" "<<grfunc.m_grx<<" "<<((nx[i]*px[i])/grfunc.m_ni2)<<" "<<grfunc.m_gr<<" rhs "<<rhs[j])
}
// this needs a check if loop runs to limits 
// this is the -1, which is now off minus 2
d2diag[j-2]=(-eid-vtd2*(1.0-dyx));
u2diag[j]=(eid-vtd2*(1.0+dyx));
d1diag[j-1]=-grfunc.m_dxy -qe*nx[i]; // -xterm/px[i];
//d1diag[j-1]=   - np*(qe+fi);
u1diag[j]=0;

} // p scope

++j;


} // i
// boundari
diag[0]=diag[2]; diag[1]=diag[3];
u1diag[0]=0;
u2diag[0]=0;
u1diag[1]=0;
u2diag[1]=0;
// wtf? 
d1diag[0]=0;
// now lets see if the others work lol
const IdxTy endfudd=szn-1;
diag[endfudd]=diag[endfudd-2];
diag[endfudd-1]=diag[endfudd-3];
d1diag[endfudd-1]=0;
d1diag[endfudd-2]=0;
u1diag[endfudd-1]=0;
d2diag[endfudd-2]=0;
d2diag[endfudd-3]=0;

++siter;


} // make_bd

// gauss eliminate efield derivatives with only carriers of the same
// type. This adds another diagonal but hopefully works better.

// this should take consecutives and gaus eliminate the grad(e)*grad(y)
// terms but it does not seem to work well. Version above adds another
// a diagonal reducing adjacent eqns for the same carriers 
static void do_eliminate_efield(blas_diags & bd,int & KU, const Values & nx,const Values & px,
const Values & x,const Values & y,const IdxTy & szn,const Value & qef,
const Value & dxi2)
{
double *  grade= new double[szn];
double *  grady= new double[szn];
IdxTy j=1;
for ( IdxTy i=2; i<(szn-2); i+=2)
{
grade[i]=-qef*nx[j];
grade[i+1]=qef*px[j];
const Value dydx=(y[j+1]-y[j-1])*dxi2;
const Value dxdx=(x[j+1]-x[j-1])*dxi2;
grady[i]=dydx;
grady[i+1]=dxdx;
++j;
}
grade[0]=0; grade[1]=0; grade[szn-2]=0; grade[szn-1]=0;
grady[0]=0; grady[1]=0; grady[szn-2]=0; grady[szn-1]=0;
reduce_efield( bd, KU, grade, grady, 0);
delete [] grade;
delete [] grady;

} //do_eliminate_efield



 
static void make_rhs(BlasTy * rhs,  const Values & nx, const Values & px,
const Values & Nteff,const Values & y, const Values & x, const Values & e, 
GRFunction & grfunc
, const Value & dxi2,const Value & dxfi, const Value & qe,
const Value & c,const IdxTy & sz,const IdxTy & ncar
, const bool make_rhs_dy_dt, const bool assert_finite
)


{
static IdxTy siter=0;
const IdxTy szn=sz*ncar;
const IdxTy sii=1;
const IdxTy siinc=ncar*sii;
const IdxTy sff=sz-1;
const IdxTy sffnc=sff*ncar;
IdxTy ptr=siinc; // actual pointer into combo area
for ( IdxTy i=sii; i<(sff); ++i)
{
const IdxTy in=i+1;
const IdxTy ip=i-1;
const Value dydx=(y[in]-y[ip])*dxi2;
const Value lapn=(y[in]+y[ip]-2.0*y[i])*dxfi;
//const Value fi=1.0/fd[i];
//const Value gr= -fi*(ni2-nx[i]*px[i]);
// zero flag could avoid derivative calculation 
grfunc.set(nx[i],px[i],0);
const Value rhoq=(nx[i]-px[i]+Nteff[i])*(-qe);


const Value dxdx=(x[in]-x[ip])*dxi2;
const Value lapp=(x[in]+x[ip]-2.0*x[i])*dxfi;

// eh, either this of change the sign... 
// this is the error, NOT the rhs unless we change sign of
// result... doh
// this is actually the time derivative of the carrier concentration
// but div by nx makes it log which would seem to be what we
// want... 
// hopefullyu compiler can move this crap lol 
if ( verbose_crap){
MM_MSG(" rhsrhs "<<siter<<" "<<i<<" "<<e[i]<<" "<<dydx<<" "<<dxdx<<" "<<rhoq<<" gr "<<grfunc.m_gry<<" "<<grfunc.m_grx<<" la "<<lapn<<" "<<lapp<<" "<<nx[i]<<" "<<px[i])
}

if ( make_rhs_dy_dt )
{
if ( nx[i]!=0) rhs[ptr]=(e[i]*dydx +c*(lapn+dydx*dydx)+rhoq) +grfunc.m_gry; // gr/nx[i];  
else rhs[ptr]=(e[i]*dydx +c*(lapn+dydx*dydx)+rhoq) ;

if ( assert_finite) { if ( !isfinite(rhs[ptr])) {
MM_MSG("notfinite "<<siter<<" "<<ptr<<" "<<rhs[ptr]<<" "<<e[i]<<" "<<dydx<<" "<<" "<<lapn<<" "<<rhoq<<" "<<grfunc.m_gry<<" np= "<<nx[i]<<","<<px[i]<<" yy "<<y[in]<<" "<<y[ip]<<" "<<nx[in]<<" "<<nx[ip])
} } // assert_finite
  
++ptr;
if ( px[i]!=0) rhs[ptr]=(e[i]*dxdx -c*(lapp+dxdx*dxdx)+rhoq) - grfunc.m_grx; // gr/px[i]; 
else rhs[ptr]=(e[i]*dxdx -c*(lapp+dxdx*dxdx)+rhoq) ; 

if ( assert_finite) { if ( !isfinite(rhs[ptr])) {
MM_MSG("notfinitep "<<siter<<" "<<ptr<<" "<<rhs[ptr]<<" "<<e[i]<<" "<<dydx<<" "<<e[i]<<" "<<lapn<<" "<<rhoq<<" "<<grfunc.m_gry<<" "<<px[i])
} } // assert_finite
  


++ptr;
}
else{
rhs[ptr]=nx[i]*(e[i]*dydx +c*(lapn+dydx*dydx)+rhoq) +grfunc.m_gr; // gr;  
++ptr;
rhs[ptr]=px[i]*(e[i]*dxdx -c*(lapp+dxdx*dxdx)+rhoq) - grfunc.m_gr; // gr; 
++ptr;



} // if ( make_rhs_dy_dt ) 



} //i 
for ( IdxTy i=0; i<siinc; ++i) {rhs[i]=0;  } 
for ( IdxTy  i=sffnc; i<IdxTy(szn); ++i) { rhs[i]=0; } 

++siter;

} // make_rhs

// different set of ewns
// uses solution's electric field, a function of the constant 
// current flux. 
static void solve_difference_eq_eflux3(Values & dy, Values & dxc, const MatPr & pr, const GP & gp, const Values & y, const Values & x
,const Values & nx, const Values & px, const Values &Nteff
, const Values & e,  GRFunction & grfunc, const Value & c, 
const sim_flags & flags, JunkBin & junk, const IdxTy solverflags)
{

static IdxTy siter=0;
const Grid & g=y.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 
bool flagx=!true;
if ( (solverflags&1)!=0)  flagx=true; 
const bool dump_terms=flagx;
const bool dump_solution=flagx;

const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
const Value dxfi2=2.0*dxfi;
// for alts ohly
//const Value dxi=1.0/dx;
const Value dxi2=.5/dx;
//const Value vtd2=c*dxfi; // kt/q/del^2
int ncar=2;
int szn=ncar*sz;
int KL=3;
int  KU=3;
//MM_MSG(" using the eflux thing")
blas_diags bd(szn,KL,KU);

BlasTy * rhs =bd.rhs;
BlasTy * subdiag =bd.d1diag;
BlasTy * sub2diag =bd.d2diag;
BlasTy * sub3diag =bd.d3diag;
BlasTy * diag =bd.diag;
BlasTy * superdiag =bd.u1diag;
BlasTy * super2diag =bd.u2diag;
BlasTy * super3diag =bd.u3diag;

// old test with flux useless, try again 
const Value flux=0;
const Value fluxdxi2=flux*dxi2;

const IdxTy bound=3;
const IdxTy start=0+bound;
const IdxTy end=sz-bound;
const IdxTy starter=ncar*start;
IdxTy j=starter;
const IdxTy end2=szn-j;

for(IdxTy i=start; i<end; ++i)
{
const Value n=nx[i];
const Value p=px[i];
// this probably works with the old one 
//grfunc.set_alt(n,p,1);
grfunc.set_alt2(n,p,1);
const Value sum=n+p;
const Value sumi=1.0/sum; // should never be zero... 
const Value sumi2=sumi*sumi; // should never be zero... 
const Value prosumi2=n*p*sumi2;
const Value cdxfisum=c*dxfi*sum;
const Value cdxi2=c*dxi2;
const Value cdxfi2sum=c*dxfi2*sum;
// this could benefit from f(x,y) vs -f(y,x)
const  Value rhoterm=qe*(p-n-Nteff[i]);
// in equillibrubium, gradx=-grady. 
const Value gradx=dxi2*(x[i+1]-x[i-1]);
const Value grady=dxi2*(y[i+1]-x[i-1]);

const Value clapx=(x[i+1]+x[i-1]-2.0*x[i])*dxfi*c;
const Value clapy=(y[i+1]+y[i-1]-2.0*y[i])*dxfi*c;

const Value gradsum=gradx+grady;
// check grfunc, change convention these are all POSTIVIE
const Value enp=grfunc.m_gry+flux*grady*sumi+c*p*gradsum*grady*sumi
+clapy+rhoterm;
const Value enp_dy=grfunc.m_dyy-flux*grady*n*sumi2
-c*prosumi2*gradsum*grady -c*dxfi2-qe*n;
const Value enp_dyp=dxi2*flux*sumi+dxi2*c*p*sumi*(gradx+2.0*grady)+c*dxfi;
const Value enp_dym=-dxi2*flux*sumi-dxi2*c*p*sumi*(gradx+2.0*grady)+c*dxfi;

const Value enp_dx=grfunc.m_dyx-flux*grady*p*sumi2
+c*prosumi2*gradsum*grady+qe*p;
const Value enp_dxp=dxi2*c*p*sumi*grady;
const Value enp_dxm=-enp_dxp;


const Value epp=grfunc.m_grx+flux*gradx*sumi-c*n*gradsum*gradx*sumi
-clapx+rhoterm;
const Value epp_dx=grfunc.m_dxx-flux*gradx*p*sumi2
+c*prosumi2*gradsum*gradx+c*dxfi2+qe*p;
const Value epp_dxp=dxi2*flux*sumi-dxi2*c*n*sumi*(grady+2.0*gradx)-c*dxfi;
const Value epp_dxm=-dxi2*flux*sumi+dxi2*c*n*sumi*(grady+2.0*gradx)-c*dxfi;
const Value epp_dy=grfunc.m_dxy-flux*gradx*n*sumi2
-c*prosumi2*gradsum*gradx-qe*n;
const Value epp_dyp=-dxi2*c*n*sumi*gradx;
const Value epp_dym=-epp_dyp;

/*
if ( dump_terms)
{
MM_MSG("pieces "<<siter<<" "<<j<<" An "<<An<<" Bn "<<Bn<<" Bndxminus "<<Bndxminus<<" Bndyminus "<<Bndyminus<<" Andy "<<Andy<<" Bndy "<<Bndy<<" Bndx "<<Bndx<<" Andx "<<Andx <<" BndyPlus "<<Bndyplus<<" Bndxplus "<<Bndxplus)
 
} // dump_terms
*/
// y is first
 rhs[j] =enp;
 subdiag[j-1] =enp_dxm; // Bndxminus; // x-1
 sub2diag[j-2] =enp_dym; // Bndyminus; // y-1
 diag[j] =enp_dy; // Andy+Bndy; // y 
 superdiag[j] =enp_dx; // Bndx+Andx; // x
 super2diag[j] =enp_dyp; // Bndyplus; // y+1
 super3diag[j] =enp_dxp; // Bndxplus; // x+1

// second half
++j;
/*
if ( dump_terms)
{
MM_MSG(" pieces "<<siter<<" "<<j<<" Ap "<<Ap<<" Bp "<<Bp<<" Bpdy "<<Bpdy<<" Apdy "<<Apdy<<" Bpdxminus "<<Bpdxminus<<" Bpdyminus "<<Bpdyminus<<" Apdx "<<Apdx<<" Bpdx "<<Bpdx <<" BpdyPlus "<<Bpdyplus<<" Bpdxplus "<<Bpdxplus)
}
*/
 rhs[j] =epp; // x
 subdiag[j-1] =epp_dy; // Bpdy+Apdy; // y 
 sub2diag[j-2] =epp_dxm; // Bpdxminus; // x-1
sub3diag[j-3] =epp_dym; // Bpdyminus; // y-1
 diag[j] =epp_dx; // Apdx+Bpdx; // x
 superdiag[j] =epp_dyp; // Bpdyplus; // y+1
 super2diag[j] =epp_dxp; // Bpdxplus; // x+1


++j;
} //i 

// fix the stupid numerical issues note this dies if starter<szn
// this also needs to make derivative zero ZZ
for(IdxTy i=0; i<IdxTy(starter); ++i){ diag[i]=diag[starter]; rhs[i]=0;}
for(IdxTy i=IdxTy(end2); i<IdxTy(szn); ++i){ diag[i]=diag[end2-1]; rhs[i]=0;}
// not part of solve now
bd.precondition(starter,end2);
bd.solve(dump_solution,siter);
bd.decondition();

// there is another adjusment later lol... 
// this had been checking for finte too
// this needs to project out the thing that kills
// net charge==0-
 j=0;
for (IdxTy i=0; i<(sz); ++i) { 
dy[i]=- fixnp(rhs[j]); ++j; 
dxc[i]=-fixnp(rhs[j]); ++j; 
} // i

++siter;

}
static void update(Values  & nx, Values & px,
	const Values & dy, const Values & dx
	, const Values & Nteff , const sim_flags & flags, JunkBin & junk
)
{



} // update


// new and improved eflux



// new and improved eflux

static void solve_difference_eq_eflux2(Values & dy, Values & dxc, const MatPr & pr, const GP & gp, const Values & y, const Values & x
,const Values & nx, const Values & px, const Values &Nteff
, const Values & e,  GRFunction & grfunc, const Value & c, 
const sim_flags & flags, JunkBin & junk, const IdxTy solverflags)
{

static IdxTy siter=0;
const Grid & g=y.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 
bool flagx=!true;
if ( (solverflags&1)!=0)  flagx=true; 
const bool dump_terms=flagx;
const bool dump_solution=flagx;

const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
const Value dxfi2=2.0*dxfi;
// for alts ohly
//const Value dxi=1.0/dx;
const Value dxi2=.5/dx;
//const Value vtd2=c*dxfi; // kt/q/del^2
int ncar=2;
int szn=ncar*sz;
int KL=3;
int  KU=3;
//MM_MSG(" using the eflux thing")
blas_diags bd(szn,KL,KU);

BlasTy * rhs =bd.rhs;
BlasTy * subdiag =bd.d1diag;
BlasTy * sub2diag =bd.d2diag;
BlasTy * sub3diag =bd.d3diag;
BlasTy * diag =bd.diag;
BlasTy * superdiag =bd.u1diag;
BlasTy * super2diag =bd.u2diag;
BlasTy * super3diag =bd.u3diag;

// old test with flux useless, try again 
const Value flux=0;
const Value fluxdxi2=flux*dxi2;

const IdxTy bound=3;
const IdxTy start=0+bound;
const IdxTy end=sz-bound;
const IdxTy starter=ncar*start;
IdxTy j=starter;
const IdxTy end2=szn-j;

for(IdxTy i=start; i<end; ++i)
{
const Value n=nx[i];
const Value p=px[i];
grfunc.set_alt(n,p,1);
const Value sum=n+p;
const Value cdxfisum=c*dxfi*sum;
const Value cdxi2=c*dxi2;
const Value cdxfi2sum=c*dxfi2*sum;
// this could benefit from f(x,y) vs -f(y,x)
const  Value rhoterm=qe*(p*p-n*n-Nteff[i]*(sum));
const Value rhotermx=qe*(2.0*p*p-Nteff[i]*p);
const Value rhotermy=qe*(-2.0*n*n-Nteff[i]*n);
// in equillibrubium, gradx=-grady. 
const Value gradx=dxi2*(x[i+1]-x[i-1]);
const Value grady=dxi2*(y[i+1]-x[i-1]);

const Value clapx=(x[i+1]+x[i-1]-2.0*x[i])*dxfi*c;
const Value clapy=(y[i+1]+y[i-1]-2.0*y[i])*dxfi*c;

const Value gradsum=gradx+grady;
const Value An=rhoterm-grfunc.m_gry;
const Value Ap=rhoterm+grfunc.m_grx;
const Value Bn=flux*grady+c*p*gradsum*grady+clapy*sum;
const Value Bp=flux*gradx-c*n*gradsum*gradx-clapx*sum;
const Value Andy=rhotermy-grfunc.m_dyy;
const Value Andx=rhotermx-grfunc.m_dyx;
const Value Apdy=rhotermy+grfunc.m_dxy;
const Value Apdx=rhotermx+grfunc.m_dxx;
const Value Bndx=+c*p*gradsum*grady+clapy*p;
// this also has a term due to the y part of der product
const Value Bndy=-cdxfi2sum+clapy*n;
const Value Bpdy=-c*n*gradsum*gradx-clapx*n;
const Value Bpdx=cdxfi2sum-clapx*p;
const Value Bndypm=+cdxi2*p*(gradsum+grady);
const Value Bndyplus=fluxdxi2+Bndypm+cdxfisum;
const Value Bndyminus=-fluxdxi2-Bndypm+cdxfisum;

const Value Bndxpm=cdxi2*p*grady;
const Value Bndxplus=Bndxpm; // cdxi2*p*grady;
const Value Bndxminus=-Bndxpm; // -cdxi2*p*grady;


const Value Bpdxpm=+cdxi2*n*(gradsum+gradx);
const Value Bpdxminus=-fluxdxi2+Bpdxpm-cdxfisum;
const Value Bpdxplus=fluxdxi2-Bpdxpm-cdxfisum;
const Value Bpdypm=cdxi2*n*gradx;
const Value Bpdyplus=-Bpdypm;
const Value Bpdyminus=Bpdypm; // cdxi2*n*gradx;


if ( dump_terms)
{
MM_MSG("pieces "<<siter<<" "<<j<<" An "<<An<<" Bn "<<Bn<<" Bndxminus "<<Bndxminus<<" Bndyminus "<<Bndyminus<<" Andy "<<Andy<<" Bndy "<<Bndy<<" Bndx "<<Bndx<<" Andx "<<Andx <<" BndyPlus "<<Bndyplus<<" Bndxplus "<<Bndxplus)
 
} // dump_terms

// y is first
 rhs[j] =An+Bn;
 subdiag[j-1] =Bndxminus; // x-1
 sub2diag[j-2] =Bndyminus; // y-1
 diag[j] =Andy+Bndy; // y 
 superdiag[j] =Bndx+Andx; // x
 super2diag[j] =Bndyplus; // y+1
 super3diag[j] =Bndxplus; // x+1


// second half
++j;
if ( dump_terms)
{
MM_MSG(" pieces "<<siter<<" "<<j<<" Ap "<<Ap<<" Bp "<<Bp<<" Bpdy "<<Bpdy<<" Apdy "<<Apdy<<" Bpdxminus "<<Bpdxminus<<" Bpdyminus "<<Bpdyminus<<" Apdx "<<Apdx<<" Bpdx "<<Bpdx <<" BpdyPlus "<<Bpdyplus<<" Bpdxplus "<<Bpdxplus)
}
 rhs[j] =Ap+Bp; // x
 subdiag[j-1] =Bpdy+Apdy; // y 
 sub2diag[j-2] =Bpdxminus; // x-1
sub3diag[j-3] =Bpdyminus; // y-1
 diag[j] =Apdx+Bpdx; // x
 superdiag[j] =Bpdyplus; // y+1
 super2diag[j] =Bpdxplus; // x+1


++j;
} //i 

// fix the stupid numerical issues note this dies if starter<szn
for(IdxTy i=0; i<IdxTy(starter); ++i){ diag[i]=diag[starter]; rhs[i]=0;}
for(IdxTy i=IdxTy(end2); i<IdxTy(szn); ++i){ diag[i]=diag[end2-1]; rhs[i]=0;}
// not part of solve now
bd.precondition(starter,end2);
bd.solve(dump_solution,siter);
bd.decondition();

// there is another adjusment later lol... 
// this had been checking for finte too
 j=0;
for (IdxTy i=0; i<(sz); ++i) { 
dy[i]=- fixnp(rhs[j]); ++j; 
dxc[i]=-fixnp(rhs[j]); ++j; 
} // i

++siter;

}


// create n-diagonal by using the target current/flux to 
// remove e-field. 
// shoudl be obolsted by eflux2 above
static void solve_difference_eq_eflux(Values & dy, Values & dxc, const MatPr & pr, const GP & gp, const Values & y, const Values & x
,const Values & nx, const Values & px, const Values &Nteff
, const Values & e,  GRFunction & grfunc, const Value & c, 
const sim_flags & flags, JunkBin & junk, const IdxTy solverflags)
{

static IdxTy siter=0;
const Grid & g=y.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 
bool flagx=!true;
if ( (solverflags&1)!=0)  flagx=true; 
const bool dump_terms=flagx;
const bool dump_solution=flagx;

// check for bad values right awy
// they should now be prevented
//const bool assert_finite=!true;

#define MJM_INCLUDE_ALT_TERMS 0
#if MJM_INCLUDE_ALT_TERMS 
// this creates an oscillation at edge of depletion region,
// over shoot e and h but still has total charge...  
const bool use_alt_terms=!true;

#endif

const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
const Value dxfi2=2.0*dxfi;
// for alts ohly
#if MJM_INCLUDE_ALT_TERMS 
const Value dxi=1.0/dx;
#endif
const Value dxi2=.5/dx;
//const Value vtd2=c*dxfi; // kt/q/del^2
int ncar=2;
int szn=ncar*sz;
int KL=3;
int  KU=3;
//MM_MSG(" using the eflux thing")
blas_diags bd(szn,KL,KU);

BlasTy * rhs =bd.rhs;
BlasTy * subdiag =bd.d1diag;
BlasTy * sub2diag =bd.d2diag;
BlasTy * sub3diag =bd.d3diag;
BlasTy * diag =bd.diag;
BlasTy * superdiag =bd.u1diag;
BlasTy * super2diag =bd.u2diag;
BlasTy * super3diag =bd.u3diag;

// old test with flux useless, try again 
const Value flux=0;

const IdxTy bound=3;
const IdxTy start=0+bound;
const IdxTy end=sz-bound;
const IdxTy starter=ncar*start;
IdxTy j=starter;
const IdxTy end2=szn-j;
// this is a numerical mess lol see at end it is "fixed" 
//const Value diagb=1e28;
//for(IdxTy i=0; i<starter; ++i){ diag[i]=diagb; rhs[i]=0;}
//for(IdxTy i=end2; i<IdxTy(szn); ++i){ diag[i]=diagb; rhs[i]=0;}

for(IdxTy i=start; i<end; ++i)
{
const Value n=nx[i];
const Value p=px[i];
grfunc.set_alt(n,p,1);
const Value sum=n+p;
const  Value rhoterm=qe*(p*p-n*n-Nteff[i]*(sum));
const Value rhotermx=qe*(2.0*p*p-Nteff[i]*p);
const Value rhotermy=qe*(-2.0*n*n-Nteff[i]*n);

const Value gradx=dxi2*(x[i+1]-x[i-1]);
const Value grady=dxi2*(y[i+1]-x[i-1]);

#if MJM_INCLUDE_ALT_TERMS
// the grad terms appears as squares, it may be helpful to write
// these as fwd*back 
// not well documented in text, possibly helping but I may 
// have also reduced increment limits preventing getting
// into basin of spurious roots
// this causes pre-echo like oscillations at depletion edge
const Value gradxf=dxi*(x[i+1]-x[i]);
const Value gradyf=dxi*(y[i+1]-x[i]);
const Value gradxb=dxi*(x[i]-x[i-1]);
const Value gradyb=dxi*(y[i]-x[i-1]);
const Value gradxfb=gradxf*gradxb;
const Value gradyfb=gradyf*gradyb;
const Value gradxy=gradx*grady;
#endif
// end of alt

const Value lapx=(x[i+1]+x[i-1]-2.0*x[i])*dxfi;
const Value lapy=(y[i+1]+y[i-1]-2.0*y[i])*dxfi;

const Value gradsum=gradx+grady;
const Value An=rhoterm-grfunc.m_gry;
const Value Ap=rhoterm+grfunc.m_grx;
const Value Bn=flux*grady+c*p*gradsum*grady+c*lapy*sum;
const Value Bp=flux*gradx-c*n*gradsum*gradx-c*lapx*sum;
const Value Andy=rhotermy-grfunc.m_dyy;
const Value Andx=rhotermx-grfunc.m_dyx;
const Value Apdy=rhotermy+grfunc.m_dxy;
const Value Apdx=rhotermx+grfunc.m_dxx;
const Value Bndx=+c*p*gradsum*grady+c*lapy*p;
// this also has a term due to the y part of der product
const Value Bndy=-c*dxfi2*sum+c*lapy*n;
const Value Bpdy=-c*n*gradsum*gradx-c*lapx*n;
const Value Bpdx=c*dxfi2*sum-c*lapx*p;
const Value Bndyplus=flux*dxi2+c*p*dxi2*(gradsum+grady)+c*dxfi*sum;
const Value Bndyminus=-flux*dxi2-c*p*dxi2*(gradsum+grady)+c*dxfi*sum;
// this sign seems wrong on gradyb
const Value Bndxplus=c*p*dxi2*grady;
const Value Bndxminus=-c*p*dxi2*grady;

const Value Bpdxminus=-flux*dxi2+c*n*dxi2*(gradsum+gradx)-c*dxfi*sum;
const Value Bpdxplus=flux*dxi2-c*n*dxi2*(gradsum+gradx)-c*dxfi*sum;
const Value Bpdyplus=-c*n*dxi2*gradx;
const Value Bpdyminus=c*n*dxi2*gradx;

#if MJM_INCLUDE_ALT_TERMS
const Value Bn_alt=flux*grady+c*p*(gradxy+gradyfb)+c*lapy*sum;
const Value Bp_alt=flux*gradx-c*n*(gradxy+gradxfb)-c*lapx*sum;
//const Value Bndy_alt=-c*dxfi2*sum+c*lapy*n+c*p*lapy;
const Value Bndy_alt=c*sum*(lapy-dxfi2);
//const Value Bpdx_alt=c*dxfi2*sum-c*lapx*p-c*lapx*n;
const Value Bpdx_alt=-c*sum*(lapx-dxfi2);
const Value Bndx_alt=+c*p*(gradxy+gradyfb+lapy);
const Value Bpdy_alt=-c*n*(gradxy+gradxfb+lapx);
//const Value Bndyplus_alt=flux*dxi2+c*p*dxi2*(gradsum+gradyb*2.0)+c*dxfi*sum;
const Value Bndyplus_alt=flux*dxi2+c*p*dxi2*(gradx+2*gradyb)+c*dxfi*sum;
//const Value Bpdxplus_alt=flux*dxi2-c*n*dxi2*(gradsum+gradxf*2.0)-c*dxfi*sum;
const Value Bpdxplus_alt=flux*dxi2-c*n*dxi2*(grady+2*gradxb)-c*dxfi*sum;
//const Value Bndyminus_alt=-flux*dxi2-c*p*dxi2*(gradsum-gradyb*2.0)+c*dxfi*sum;
const Value Bndyminus_alt=-flux*dxi2-c*p*dxi2*(gradx+2*gradyf)+c*dxfi*sum;
//const Value Bpdxminus_alt=-flux*dxi2+c*n*dxi2*(gradsum-gradxb*2.0)-c*dxfi*sum;
const Value Bpdxminus_alt=-flux*dxi2+c*n*dxi2*(grady+2*gradxf)-c*dxfi*sum;

#endif

if ( dump_terms)
{
MM_MSG("pieces "<<siter<<" "<<j<<" An "<<An<<" Bn "<<Bn<<" Bndxminus "<<Bndxminus<<" Bndyminus "<<Bndyminus<<" Andy "<<Andy<<" Bndy "<<Bndy<<" Bndx "<<Bndx<<" Andx "<<Andx <<" BndyPlus "<<Bndyplus<<" Bndxplus "<<Bndxplus)
 
} // dump_terms

// y is first
 rhs[j] =An+Bn;
 subdiag[j-1] =Bndxminus; // x-1
 sub2diag[j-2] =Bndyminus; // y-1
 diag[j] =Andy+Bndy; // y 
 superdiag[j] =Bndx+Andx; // x
 super2diag[j] =Bndyplus; // y+1
 super3diag[j] =Bndxplus; // x+1

#if MJM_INCLUDE_ALT_TERMS
if ( use_alt_terms)
{
if ( dump_terms)
{
MM_MSG("piecesalt "<<siter<<" "<<j<<" Bn_alt "<<Bn_alt<<" Bndyminus_alt "<<Bndyminus_alt<<" Bndy_alt "<<Bndy_alt<<" Bndx_alt "<<Bndx_alt<<" Bndyplus_alt "<<Bndyplus_alt )   
}

 rhs[j] =An+Bn_alt;
 subdiag[j-1] =Bndxminus; // x-1
 sub2diag[j-2] =Bndyminus_alt; // y-1
 diag[j] =Andy+Bndy_alt; // y 
 superdiag[j] =Bndx_alt+Andx; // x
 super2diag[j] =Bndyplus_alt; // y+1
 super3diag[j] =Bndxplus; // x+1
}
#endif

// second half
++j;
if ( dump_terms)
{
MM_MSG(" pieces "<<siter<<" "<<j<<" Ap "<<Ap<<" Bp "<<Bp<<" Bpdy "<<Bpdy<<" Apdy "<<Apdy<<" Bpdxminus "<<Bpdxminus<<" Bpdyminus "<<Bpdyminus<<" Apdx "<<Apdx<<" Bpdx "<<Bpdx <<" BpdyPlus "<<Bpdyplus<<" Bpdxplus "<<Bpdxplus)
}
 rhs[j] =Ap+Bp; // x
 subdiag[j-1] =Bpdy+Apdy; // y 
 sub2diag[j-2] =Bpdxminus; // x-1
sub3diag[j-3] =Bpdyminus; // y-1
 diag[j] =Apdx+Bpdx; // x
 superdiag[j] =Bpdyplus; // y+1
 super2diag[j] =Bpdxplus; // x+1

#if MJM_INCLUDE_ALT_TERMS
if ( use_alt_terms)
{
if ( dump_terms)
{
MM_MSG("piecesalt "<<siter<<" "<<j<<" Bp_alt "<<Bp_alt<<" Bpdy_alt "<<Bpdy_alt<<" Bpdxminus_alt "<<Bpdxminus_alt<<" Bpdx_alt "<<Bpdx_alt<<" Bpdxplus_alt "<<Bpdxplus_alt )   
}

 rhs[j] =Ap+Bp_alt; // x
 subdiag[j-1] =Bpdy_alt+Apdy; // y 
 sub2diag[j-2] =Bpdxminus_alt; // x-1
sub3diag[j-3] =Bpdyminus; // y-1
 diag[j] =Apdx+Bpdx_alt; // x
 superdiag[j] =Bpdyplus; // y+1
 super2diag[j] =Bpdxplus_alt; // x+1



}
#endif


++j;
} //i 

// fix the stupid numerical issues note this dies if starter<szn
for(IdxTy i=0; i<IdxTy(starter); ++i){ diag[i]=diag[starter]; rhs[i]=0;}
for(IdxTy i=IdxTy(end2); i<IdxTy(szn); ++i){ diag[i]=diag[end2-1]; rhs[i]=0;}
// no longer part of solve
bd.precondition(starter,end2);
bd.solve(dump_solution,siter);
bd.decondition();

/*
double * rhsx= 0;

{ // blas scoping 
bd.precondition();
// this still needs 5 diags 
BlasTy * diagsfudd  =bd.reform();
// this method does likely not work 
//if ( verbose_crap ) { bd.dump_msg(siter); }
int NRHS=1;
int INFO=99999;
int N = szn; 
if ( dump_solution)
{rhsx= new double[szn];
::memcpy(rhsx,rhs,szn*sizeof(double));
}
// this can move to become a bd method
// pivot does not appear tobe needed to sort out results lol
int * ipiv = new int[szn] ; // this does not appear to be needed
 dgbsv_(&N,&KL, & KU, &NRHS, diagsfudd, &bd.szwork,ipiv, rhs,  &bd.szn ,&INFO);
if (INFO!=0) { { MM_MSG(siter<<" DGBSV error code is "<<INFO) } 
if ( INFO>0) {  MM_MSG(" info egn "<<bd.dump_eqn(INFO-1)) } exit(-1); }
delete [] ipiv; 
bd.decondition();
} // blas scoping 
if ( dump_solution)
{ for (IdxTy i=0; i<IdxTy(szn); ++i) { 
MM_MSG("solutx "<<siter<<" "<<i<<bd.dump_eqn(i)<<" rhx "<<rhsx[i]<<" = "<<rhs[i])
} //i 
delete [] rhsx;
}
*/

// there is another adjusment later lol... 
 j=0;
for (IdxTy i=0; i<(sz); ++i) { 
dy[i]=- fixnp(rhs[j]); ++j; 
dxc[i]=-fixnp(rhs[j]); ++j; 
} // i

++siter;

} //solve_difference_eq_eflux

// avoid going to spurious root right away
static Value fixnp(const Value & dx)
{
// iirc, as G gets larger the spurious min/max are closer 
const Value fudge=.8;
const Value mmax=.5;
const Value mmin=-.5;
Value x=dx;
if ( true)
{
// this helped when there was a missing line in the code
// that added dx,y to x and y but no idea if it helps now
// right now the inverse blows up with high GR
//Value xlog=::log(1.0+::fabs(dx));
Value xlog=(0.0+::fabs(dx));
if ( dx>0) x=xlog;
else x= -xlog;
}
if ( x>mmax) x=mmax;
if (x<mmin ) x=mmin;

return x*fudge; 
}
// as below, but attempt to do n and p together to make
// a band diagonal bigger than tridiagonal and 2x in size.
// this avoids solving n and p seprately and still
// ignore the off-digamonal electric field.
// format is interleaved, n0,p0,n1,p1 .. 

static void solve_difference_eq_together(Values & dy, Values & dxc, const MatPr & pr, const GP & gp, const Values & y, const Values & x
,const Values & nx, const Values & px, const Values &Nteff
, const Values & e,  GRFunction & grfunc, const Value & c, const sim_flags & flags, JunkBin & junk)
{
static IdxTy siter=0;
const Grid & g=y.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 
// probably the value at n has coeff zero and the off diags
// should go all the way on both sides, need a way to solve
// they are small but there are a lot of them. 
//const bool efield=true;

const bool  eliminate_efield=true;
// originally I divved by exp(y) or nx but on first attempt
// this failed although it looks like it was better
// this means removing exp(y) from rhs and backing out
// partials for matrix elements. 
// todo : xterm, efield, gr derivative but so far so good
// the xterm may be ok, should get boundary conditions although
// probably not a big deal right now. 
const bool make_rhs_dy_dt=true; 

// check for bad values right awy
const bool assert_finite=true;

//const Value ni2=pr.ni2;
// the fudding gr term is an assfudd
//const Value ni=pr.ni;

const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
//const Value dxi=1.0/dx;
const Value dxi2=.5/dx;
const Value vtd2=c*dxfi; // kt/q/del^2

//MM_MSG(" test "<<sz<<" "<<dxfi<<" "<<dxf)
//const IdxTy si=flags.fixed_start;
//const IdxTy sf=sz-flags.fixed_end;
// there is no point in making grids, just use contig memory

int ncar=2;
int szn=ncar*sz;
// 2 upper and lowe diags,+/-2 are the adjacents, +1 issame p for n
// -1 is same n for p   
int NRHS=1;
int KL=ncar;
int  KU=ncar;
blas_diags bd(szn,KL,KU);

BlasTy * rhs =bd.rhs;

make_rhs(rhs,nx,px,Nteff,y,x,e,grfunc,dxi2,dxfi,qe,c,sz,ncar,make_rhs_dy_dt,assert_finite );


// this needs a graceful way to handle gr and boundary conditions 
// skip KL unused rows then the KU slots
const bool use_old_crap=false;

if  (!use_old_crap)
{
make_bd(bd,rhs,nx,px,Nteff,y,x,e,grfunc,dxi2,vtd2,qe,c,sz,assert_finite);

}
// leaving out the xy coupling and the boundary terms apparently
// is singular, either decouple or nan in unused area lol 
else {  // use_old_crap
make_old_crap(bd,rhs,nx,px,Nteff,y,x,e,grfunc,dxi2,vtd2,qe,c,sz,assert_finite); }

if ( eliminate_efield)
{
const Value qef=qe*dx;
do_eliminate_efield(bd,KU,nx,px,x,y,szn,qef,dxi2);
//++KU;

}


BlasTy * diagsfudd  =bd.reform();

if ( verbose_crap ) { bd.dump_msg(siter); }

int INFO=99999;
int N = szn; 
// pivot does not appear tobe needed to sort out results lol
int * ipiv = new int[szn] ; // this does not appear to be needed
 dgbsv_(&N,&KL, & KU, &NRHS, diagsfudd, &bd.szwork,ipiv, rhs,  &bd.szn ,&INFO);
if (INFO!=0) { { MM_MSG(siter<<" DGBSV error code is "<<INFO) } 
if ( INFO>0) {  MM_MSG(" info egn "<<bd.dump_eqn(INFO-1)) } exit(-1); }
const bool FUDD=false; 
if ( FUDD)
{ for (IdxTy i=0; i<IdxTy(szn); ++i) {  
if ( (ipiv[i]-1)!=i) 
{const Value fu=rhs[i]; rhs[i]=rhs[ipiv[i]-1]; rhs[ipiv[i]-1]=fu;  } } }

IdxTy j=0;
for (IdxTy i=0; i<(sz); ++i) { 
 // MM_MSG(i<<" "<<j<<" "<<rhs[j]<<" "<<rhs[j+1]<<" ipiv "<<ipiv[j])
// we never changed the sign of the rhs, so the solution points
// in the wrong dir doh
//if ( (ipiv[j]-1)!=j) {MM_MSG(" ipiv_diff "<<siter<<" j= "<<j<<" ipiv= "<<(ipiv[j]-1)) } 
dy[i]=- rhs[j]; ++j; 
if ( assert_finite) { if ( !isfinite(rhs[j])) { 
MM_MSG(" notfinite "<<siter<<" "<<j<<" "<<rhs[j]<<" "<<bd.dump_eqn(j)) 
}
} // assert_tinite
//if ( (ipiv[j]-1)!=j) {MM_MSG(" ipiv_diff "<<siter<<" j= "<<j<<" ipiv= "<<(ipiv[j]-1)) } 
dxc[i]=-rhs[j]; ++j; 
if ( assert_finite) { if ( !isfinite(rhs[j])) { 
MM_MSG(" notfinite "<<siter<<" "<<j<<" "<<rhs[j]<<" "<<bd.dump_eqn(j)) 
}
} // assert_tinite
} // i  

delete [] ipiv;

++siter;

} // together


// for the non-zed current case, try to iterate n,p, and electric
// field to avoid lower triangular matrix
// solve for dy with x being the fixed "other" although 
// it would be nice to try lower diagonal 2x matrix solving
// field, eletroncs, and holes at same time. 
// "c" is kt/q with sign change etc 
// fd is not the fermi but demoninator of generation rate assumed
// to be np-ni2 in numerator 
// note that nx and px cn be interchanged n and p 
static void solve_difference_eq(Values & dy, Values & dxc, const MatPr & pr, const GP & gp, const Values & y, const Values & x
,const Values & nx, const Values & px, const Values &Nteff
, const Values & e, const Values & fd, const Value & c, const sim_flags & flags, JunkBin & junk)
{
static IdxTy siter=0;
const Grid & g=y.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 
const Value ni2=pr.ni2;
// the fudding gr term is an assfudd
const Value ni=pr.ni;

const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
//const Value dxi=1.0/dx;
const Value dxi2=.5/dx;

const Value vtd2=c*dxfi; // kt/q/del^2

//MM_MSG(" test "<<sz<<" "<<dxfi<<" "<<dxf)
//const IdxTy si=flags.fixed_start;
//const IdxTy sf=sz-flags.fixed_end;

const IdxTy sii=1;
const IdxTy sff=sz-1;

//const IdxTy szeff=sf-si;
// try different n and p values and see how similar they are 
Values nnew(g),pnew(g),errn(g),errp(g);

// calculate the current non-solution ( error )
for ( IdxTy i=sii; i<(sff); ++i)
{
const IdxTy in=i+1;
const IdxTy ip=i-1;
const Value dydx=(y[in]-y[ip])*dxi2;
const Value lapn=(y[in]+y[ip]-2.0*y[i])*dxfi;
const Value fi=1.0/fd[i];
// with fi<0, gr is positive when depleted. 
 const Value gr= -fi*(ni2-nx[i]*px[i]);


const Value rhoq=(nx[i]-px[i]+Nteff[i])*(-qe);
// swapping nx and px will fudd up rho.... 
// if nx is zero, this blows up 
//errn[i]=e[i]*dydx +c*(lapn+dydx*dydx)+(nx[i]-px[i]+Nteff[i])*(-qe)
// -fi*(ni2/nx[i]-px[i]);
errn[i]=nx[i]*(e[i]*dydx +c*(lapn+dydx*dydx)+rhoq)
+gr; //  -fi*(ni2-nx[i]*px[i]);
// MM_MSG(" i= "<<i<<" errn= "<<errn[i]<<" dydx= "<<dydx<<" lapn= "<<lapn <<" rhoq= " << rhoq<<" gr= "<<gr)
const Value dxdx=(x[in]-x[ip])*dxi2;
const Value lapp=(x[in]+x[ip]-2.0*x[i])*dxfi;
//errp[i]=e[i]*dxdx -c*(lapp+dxdx*dxdx)+(nx[i]-px[i]+Nteff[i])*(-qe)
// -fi*(ni2/px[i]-nx[i]);
errp[i]=px[i]*(e[i]*dxdx -c*(lapp+dxdx*dxdx)+rhoq)
- gr; //  -fi*(ni2/-px[i]*nx[i]);


// MM_MSG(" i= "<<i<<" errn= "<<errn[i]<<" errp= "<<errp[i])

}
for ( IdxTy i=0; i<sii; ++i) {errn[i]=0; errp[i]=0; } 
for ( IdxTy i=sff; i<sz; ++i) { errn[i]=0; errp[i]=0;} 


int N=sz;
int NRHS=1; // number of rhs columns
int LDB=N;
int INFO=0;
// solution must be loaded with rhs first 
double * subdiag= new double[sz];
double * diag= new double[sz];
double * superdiag= new double[sz];
double * solution= new double[sz];



//const Value dxx=-2.0*dxfi;
const Value vtd22=vtd2*2.0;
for (IdxTy i=1; i<(sz-1); ++i)
{
solution[i]=-errn[i];
//diag[i]= -vtd22+ ni2/(fd[i]*nx[i])-qe*nx[i];
// this is confusing but the gr term seems right,
// errn include both pars, the ni2 term here removes
// the equ part .
const Value fi=1.0/fd[i];
// with fi<0, gr is positive when depleted. 
 const Value gr= -fi*(ni2-nx[i]*px[i]);
const Value gr_part2= gr/(nx[i]+px[i]+2*ni); // this assume form of fd... doh 
diag[i]= nx[i]*(-vtd22-qe*nx[i])+ni2/fd[i]+errn[i]- gr_part2*nx[i];
const Value eid=e[i]*dxi2;
// this needs a check if loop runs to limits 
const Value dy=.5*(y[i+1]-y[i-1]);
//subdiag[i]=-eid+vtd2*(1.0-dy);
subdiag[i]=nx[i]*(-eid+vtd2*(1.0-dy));
//superdiag[i]=eid+vtd2*(1.0+dy);
superdiag[i]=nx[i]*(eid+vtd2*(1.0+dy));

}

const IdxTy szi=sz-1;
{
// unlike before, now we need to get both limits fixed here
//including in diag 
Value fi=1.0/fd[0];
// with fi<0, gr is positive when depleted. 
Value gr= -fi*(ni2-nx[0]*px[0]);
 Value gr_part2= gr/(nx[0]+px[0]+2*ni); // this assume form of fd... doh 
diag[0]= nx[0]*(-vtd22 -qe*nx[0])+(ni2/fd[0])+errn[0]-gr_part2*nx[0];
fi=1.0/fd[szi];
// with fi<0, gr is positive when depleted. 
 gr= -fi*(ni2-nx[szi]*px[szi]);
 gr_part2= gr/(nx[szi]+px[szi]+2*ni); // this assume form of fd... doh 
diag[szi]= nx[szi]*(-vtd22-qe*nx[szi])+ni2/fd[szi]+errn[szi]-gr_part2*nx[szi];
// also need to get the missing ones, not just those forced to zed

subdiag[sz-2]=0;
subdiag[sz-1]=0;
superdiag[0]=0;
subdiag[0]=0;
superdiag[sz-2]=0;
}

// doh, these get over written so second call fudded up 
dgtsv_(&N,&NRHS,&subdiag[0],&diag[0],&superdiag[0],&solution[0],&LDB,&INFO);

if (INFO!=0) { MM_MSG(" DGTSV error code is "<<INFO) }


for (IdxTy i=0; i<(sz); ++i) dy[i]=solution[i];

for (IdxTy i=1; i<(sz-1); ++i)
{
solution[i]=-errp[i];

const Value fi=1.0/fd[i];
// with fi<0, gr is positive when depleted. 
 const Value gr= -fi*(ni2-nx[i]*px[i]);
const Value gr_part2= gr/(nx[i]+px[i]+2*ni); // this assume form of fd... doh 
// drift-diffusion changes sign on kt/q BUT also 
// continuity equation, so just change sign of "g" LOL
diag[i]= px[i]*(vtd22+qe*px[i])-(ni2/fd[i])+errp[i]+gr_part2*px[i];
const Value eid=e[i]*dxi2;
// this needs a check if loop runs to limits 
const Value dy=.5*(x[i+1]-x[i-1]);
subdiag[i]=px[i]*(-eid-vtd2*(1.0-dy));
superdiag[i]=px[i]*(eid-vtd2*(1.0+dy));

}
{
Value  fi=1.0/fd[0];
// with fi<0, gr is positive when depleted. 
 Value gr= -fi*(ni2-nx[0]*px[0]);

 Value gr_part2= gr/(nx[0]+px[0]+2*ni); // this assume form of fd... doh 
diag[0]= px[0]*(vtd22+qe*px[0])-(ni2/fd[0])+errp[0]+gr_part2*px[0];
 fi=1.0/fd[szi];
 gr= -fi*(ni2-nx[szi]*px[szi]);
 gr_part2= gr/(nx[0]+px[0]+2*ni); // this assume form of fd... doh 
diag[szi]= px[szi]*(vtd22+qe*px[szi])-(ni2/fd[szi])+errp[szi]+gr_part2*px[szi];

subdiag[sz-2]=0;
subdiag[sz-1]=0;
subdiag[0]=0;
superdiag[0]=0;
superdiag[sz-2]=0;
}

// doh, these get over written so second call fudded up 
dgtsv_(&N,&NRHS,&subdiag[0],&diag[0],&superdiag[0],&solution[0],&LDB,&INFO);
if (INFO!=0) { MM_MSG(" DGTSV error code is "<<INFO) }
for (IdxTy i=0; i<(sz); ++i) dxc[i]=solution[i];



delete [] diag;  delete [] subdiag;
delete [] superdiag; delete [] solution;


++siter;

}

// non-equ approach right now using fixed e field and form for the GR term.
// fixed e field preserves tri diagonal form for linearized eqn to solve.
// interesting to see if lower triangular form is tractable with 
// 1 M data points. 
static void update_carriers_4( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk, const IdxTy idx, const IdxTy solverflags )
{
const bool assert_finite=true;
const Grid & g=m.n.grid();
const IdxTy sz=g.size();
if ( sz==0) return; 

const Value ni=pr.ni;
Values dy(g), dx(g),x(g),y(g),fd(g),e(g),rho(g);
// using a fixed e field is dumb, but it keeps
// matrix tri diagnoanl although lower triangular may be ok 
const Value qe=gp.q()/(pr.epsilon*gp.epsilon_vacu());
rho=(m.p-m.n-m.Nteff)*(qe);

if ( assert_finite)
{
for (IdxTy i=0; i<sz; ++i)
{ 
bool isok=isfinite(rho[i]);
if (!isok)
{
MM_MSG( " notfinite "<< junk.iter<<" "<<i<<" "<<m.n[i]<<" "<<m.p[i]<<" "<<m.Nteff[i] )

} // isok

}

}


// this is negative and units scaled.. 
//fd=(m.p+m.n+(2*ni))*(-1e1)*(1e-26);
//fd=(m.p+m.n+(2*ni))*(-1e1)*(1e-20);
//fd=(m.p+m.n+(2*ni))*(-1e1)*(1e-18);
//fd=(m.p+m.n+(2*ni))*(-1e1)*(1e-16);

//GRFunction grfunc(1e12,ni);
// did not work with step inits but ok with zero bias smooth start
//GRFunction grfunc(1e9,ni);
// on first iter, e and h both decrease, maybe this will help 
// wow on one iter very close, good debye tail satr etc
//GRFunction grfunc(1e12,ni);
// if some is good more is better 
// looks better after 1 iter but still oscialtes
//GRFunction grfunc(1e13,ni);

//GRFunction grfunc(1e6,ni);
// see if I am backwards, 
//GRFunction grfunc(1e9,ni);
// likely creating more pre-echo 
//GRFunction grfunc(1e8,ni);
// instant nan
//GRFunction grfunc(1e7,ni);
// with 100 iters and 300k pointd, only 2 clean lobes in pre-echo
//GRFunction grfunc(1e11,ni);
// 3 loops, when others die seems the matrix solver returns nan/inf
GRFunction grfunc(1e8,ni);
// this is higher amplitude crap
//GRFunction grfunc(1e13,ni);
// nan from matrix solver
//GRFunction grfunc(1e15,ni);


// after adding more predociioners, still wrong signs with 1e6
// near boarders maybe gr too hight 
// this dies on attempt to solve
//GRFunction grfunc(1e1,ni);
// this may be slightly betteralthough h still goes wrong way 
// no good with 1k points 
//GRFunction grfunc(1e4,ni);
// turn off +/- preconditioner
//GRFunction grfunc(1e-3,ni);


//GRFunction grfunc(1e14,ni);
// never see a debye tail
//GRFunction grfunc(1e15,ni);
// with limits on the dx this sort of works for e but h is wrong
// on n side and stble on p side
//GRFunction grfunc(1e3,ni);
//GRFunction grfunc(1e14,ni);
// with 1e9 on steps, the rhs terms fluctuate a lot, numerical issue?
// GRFunction grfunc(1e-8,ni);

// if rho is infinte, this blows up
Told::I(e,g,rho);

//x=m.p.ln();
//y=m.n.ln();

x=m.p.ln(1,0);
y=m.n.ln(1,0);


// mks 
const Value c=gp.vt()*gp.q();
// swapping n and p will not work, 
if ( idx==4 ) 
{fd=(m.p+m.n+(2*ni))*(-1e1)*(1e-2);
solve_difference_eq(dy,dx,pr,gp,y,x,m.n,m.p,m.Nteff,e,fd,c,flags,junk);
}
else if ( idx==5)
solve_difference_eq_together(dy,dx,pr,gp,y,x,m.n,m.p,m.Nteff,e,grfunc,c+0*gp.vt(),flags,junk);
else if ( idx==6)solve_difference_eq_eflux2(dy,dx,pr,gp,y,x,m.n,m.p,m.Nteff,e,grfunc,0*c+gp.vt(),flags,junk,solverflags);
else if ( idx==7)solve_difference_eq_eflux(dy,dx,pr,gp,y,x,m.n,m.p,m.Nteff,e,grfunc,0*c+gp.vt(),flags,junk,solverflags);
else if ( idx==8){ solve_difference_eq_eflux3(dy,dx,pr,gp,y,x,m.n,m.p,m.Nteff,e,grfunc,0*c+gp.vt(),flags,junk,solverflags);
//continue; 
}  // do thisin themethod
else {MM_MSG("bad solver "<<idx)}


// skipped for some solvers

Values & nnew=m.n;
Values & pnew=m.p;
// if we include everything, this should be able to 
// be close to 1 just requiring more iterations
// but leaving out e field may make a mess, no idea yet
//const Value fac=.1;
const Value fac=flags.non_linear_fac;
for ( IdxTy i=0; i<sz; ++i)
{
double fn=dy[i]*fac;
if ( fn>flags.maxstep) fn=flags.maxstep;
if ( fn<flags.minstep) fn=flags.minstep;

nnew[i]=::exp(y[i]+fn);
if ( nnew[i]>flags.maxnn) nnew[i]=flags.maxnn;
if ( nnew[i]<flags.minnn) nnew[i]=flags.minnn;

// what the FUDD? This was missing???
fn=dx[i]*fac;

if ( fn>flags.maxstep) fn=flags.maxstep;
if ( fn<flags.minstep) fn=flags.minstep;

pnew[i]=::exp(x[i]+fn);
if ( pnew[i]>flags.maxnn) pnew[i]=flags.maxnn;
if ( pnew[i]<flags.minnn) pnew[i]=flags.minnn;
if ( assert_finite)
{
bool isok=isfinite(nnew[i])&&isfinite(pnew[i]);
if (!isok)
{
MM_MSG( " notfinite "<< junk.iter<<" "<<i<<" yx "<<y[i]<<" "<<dy[i]<<" "<<x[i]<<" "<<dx[i]<<" new "<<nnew[i]<<" "<<pnew[i] )

} 

} // assert_finnite

//MM_MSG(" i= "<<i<<" yi= "<<y[i]<<" dy= "<<dy[i]<<" nnew= "<<nnew[i] <<" xi= "<<x[i]<<" dx= "<<dx[i]<<" pnew= "<<pnew[i])
 
}

} // update_carriers_4



// try to linearize and solve at onece

static void update_carriers_3( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
static IdxTy siter=0;
const Value ni2=pr.ni2;
const Value ni=pr.ni;
const Grid & g=m.n.grid();
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value dxfi=1.0/dxf;
const IdxTy sz=g.size();
//MM_MSG(" test "<<sz<<" "<<dxfi<<" "<<dxf)
if ( sz==0) return; 
const IdxTy si=flags.fixed_start;
const IdxTy sf=sz-flags.fixed_end;


//const IdxTy szeff=sf-si;
Values nnew(g),pnew(g),x(g),y(g),err(g),err2(g);
x=m.p.ln();
y=m.n.ln();
const Value c=gp.q()/(gp.vt()*pr.epsilon*gp.epsilon_vacu());
// first calculate the error at each location with this
/// set of values
// I think that the error at beginning and end is zero alrad
for ( IdxTy i=0; i<si; ++i) err[i]=0;
for ( IdxTy i=0; i<si; ++i) err2[i]=0;
for ( IdxTy i=sf; i<sz; ++i) err[i]=0;
for ( IdxTy i=sf; i<sz; ++i) err2[i]=0;
// this should guard the der calc
//for ( IdxTy i=si; i<sf; ++i)
for ( IdxTy i=1; i<(sz-1); ++i)
{
const Value rhoerr=-c*(m.n[i]-m.p[i]+m.Nteff[i]);
const  Value dererr=dxfi*(y[i+1]+y[i-1]-2.0*y[i]);
const  Value dererr2=dxfi*(x[i+1]+x[i-1]-2.0*x[i]);
// this should be zero, this is the error we want to fix
err[i]=rhoerr+dererr;
err2[i]=-rhoerr+dererr2;
} //i 

err[0]=0; err[sz-1]=0; err[1]=0; err[sz-2]=0;
err2[0]=0; err2[sz-1]=0; err2[1]=0; err2[sz-2]=0;
{
int N=sz;
int NRHS=1; // number of rhs columns
int LDB=N;
int INFO=0;
const bool do_solution2=!false; 
// solution must be loaded with rhs first 
// probably a stack size limit lol 
//double subdiag[sz-1*0],diag[sz],diag2[sz],superdiag[sz-1*0], solution[sz],solution2[sz];
double * subdiag= new double[sz];
double * diag= new double[sz];
double * superdiag= new double[sz];

double * diag2= new double[sz];
double * solution= new double[sz];
double * solution2= new double[sz];



const Value dxx=-2.0*dxfi;
for (IdxTy i=0; i<sz; ++i)
{
solution[i]=-err[i];
solution2[i]=-err2[i];
const Value dfdy=-c*(m.n[i]+m.p[i]);
//const Value dfdy=0; // need tomake this work 
diag[i]=dxx+dfdy;
diag2[i]=dxx+dfdy;
// this is sloppy but easier..ja
 subdiag[i]=dxfi;
 superdiag[i]=dxfi;

}
subdiag[sz-2]=0;
subdiag[sz-1]=0;
superdiag[0]=0;
//MM_MSG(" calling ")
// doh, these get over written so second call fudded up 
dgtsv_(&N,&NRHS,&subdiag[0],&diag[0],&superdiag[0],&solution[0],&LDB,&INFO);
if ( do_solution2)
{
for (IdxTy i=0; i<sz; ++i) {subdiag[i]=dxfi; superdiag[i]=dxfi;}
subdiag[sz-2]=0;
subdiag[sz-1]=0;
superdiag[0]=0;
dgtsv_(&N,&NRHS,&subdiag[0],&diag2[0],&superdiag[0],&solution2[0],&LDB,&INFO);

} // soluion 2

if (INFO!=0)
{
MM_MSG(" DGTSV error code is "<<INFO)
}
const Value fac=1; // dx*1000; 
Value mag_sol=0;
Value mag_err=0;
Value max_sol=0;

for ( IdxTy i=0; i<sz; ++i)
{
mag_err+=err[i]*err[i];
// MM_MSG(" i= "<<i<<" x= "<<y[i]<<" soln= "<<solution[i]<<" n= "<<exp(y[i]+solution[i]))
const bool pbranch=(do_solution2)&&(m.p[i]>m.n[i]);
if (!pbranch)//if ( nnew[i]>pnew[i])
{
const Value sx=solution[i]*solution[i];
mag_sol+=sx;
if ( sx>max_sol) max_sol=sx;
solution[i]*=fac;
nnew[i]=::exp(y[i]+solution[i]);
if ( nnew[i]!=0) pnew[i]=ni2/nnew[i]; else pnew[i]=ni;
}
else
{
const Value sx=solution2[i]*solution2[i];
mag_sol+=sx;
if ( sx>max_sol) max_sol=sx;
solution2[i]*=fac;
pnew[i]=::exp(x[i]+solution2[i]);
if ( pnew[i]!=0) nnew[i]=ni2/pnew[i]; else nnew[i]=ni;


}
}
if (( siter& 511)==0) 
{
MM_MSG(" siter= "<<siter<<" errfig max_sol= "<<::sqrt(max_sol)<<" mag_sol= "<<::sqrt(mag_sol)<<" mag_err= "<<::sqrt(mag_err))
}
delete [] diag; delete [] diag2; delete [] subdiag;
delete [] superdiag; delete [] solution; delete [] solution2;


// messing up ni2 constraint in some cases
const Value fold=0*flags.nfold;
m.n.mix(nnew,fold,2,sz-2);
//m.n.mix(nnew,fold,si,sf);
m.p.mix(pnew,fold,2,sz-2);
//m.p.mix(pnew,fold,si,sf);
} // blas call 

++siter;


} // update_carriers_3


// various integrators for diff eq
static void update_carriers_2( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const bool do_backward=true;
const bool do_original=true;
const bool do_enhanced=!true;
const bool do_smoothed=true;
const Value ni2=pr.ni2;
const Value ni=pr.ni;
const Grid & g=m.n.grid();
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value sz=g.size();
Values nnew(g),pnew(g),x(g),y(g);
x=m.p.ln();
y=m.n.ln();
const Value c=gp.q()/(gp.vt()*pr.epsilon*gp.epsilon_vacu());
// shoaoting?
//Value test=y[1];
//Value hi=test+.0000001;
//Value lo=test-.001; 
//IdxTy iix=0;
//static IdxTy siter=0;
//while ( true)
//{

// y-only with p from ni2 works ok but current has problems
// try to do both and pick.. 
nnew[0]=y[0];
nnew[1]=y[1];
nnew[sz-1]=y[sz-1];
nnew[sz-2]=y[sz-2];


pnew[0]=x[0];
pnew[1]=x[1];
pnew[sz-1]=x[sz-1];
pnew[sz-2]=x[sz-2];

// this is now asymmetrical and may help convergence
if ( do_original) { 
for ( IdxTy i=2; i<(sz-2); ++i)
{
const Value rhsn=c*(m.n[i]-m.p[i]+m.Nteff[i]);
//nnew[i]=.5*(-dxf*rhsn+y[i-1]+y[i+1]);
nnew[i]=.5*(-dxf*rhsn+nnew[i-1]+y[i+1]);
//pnew[i]=.5*(dxf*rhsn+x[i-1]+x[i+1]);
pnew[i]=.5*(dxf*rhsn+pnew[i-1]+x[i+1]);
} // fwd
} // normal 

if ( do_smoothed) { 
//Values & ys=y;
Values & ys=nnew;
Values & xs=pnew;
for ( IdxTy i=3; i<(sz-3); ++i)
{
// this is obviously positive
Value drhsn=c*(m.n[i]+m.p[i]);
Value drhsnb=c*(m.n[i-1]+m.p[i-1]);
Value drhsnf=c*(m.n[i+1]+m.p[i+1]);

const Value nn2=.5*(ys[i+2]-2.0*ys[i+1]-2.0*ys[i-1]+ys[i-2]);
nnew[i]=-ys[i+1]-ys[i-1]+nn2-dxf*.5*(.5*ys[i-1]*(drhsn+drhsnb)+.5*ys[i+1]*(drhsn+drhsnf));
const Value den=-3.0-.25*dxf*(drhsn*2.0+drhsnb+drhsnf);
nnew[i]=nnew[i]/den;

const Value pp2=.5*(xs[i+2]-2.0*xs[i+1]-2.0*xs[i-1]+xs[i-2]);
pnew[i]=-xs[i+1]-xs[i-1]+pp2+dxf*.5*(.5*xs[i-1]*(drhsn+drhsnb)+.5*xs[i+1]*(drhsn+drhsnf));

const Value dep=-3.0+.25*dxf*(drhsn*2.0+drhsnb+drhsnf);
pnew[i]=pnew[i]/dep;



} //i 

} // smoothed


if ( do_enhanced) { 
const Value dxf3=dxf/3.0;
Value rhsnlast=c*(m.n[1]-m.p[1]+m.Nteff[1]);
Value rhsn=c*(m.n[2]-m.p[2]+m.Nteff[2]);
for ( IdxTy i=2; i<(sz-2); ++i)
{
Value rhsnnext=c*(m.n[i+1]-m.p[i+1]+m.Nteff[i+1]);
Value drhsn=c*(m.n[i]+m.p[i]);
Value drhsp=c*(-m.n[i]+m.p[i]);
//nnew[i]=.5*(-dxf*rhsn+y[i-1]+y[i+1]);
nnew[i]=.5*(-dxf3*(rhsn+rhsnnext+rhsnlast+drhsn*(2.0*y[i]-y[i+1]-nnew[i-1]))+nnew[i-1]+y[i+1]);
//pnew[i]=.5*(dxf*rhsn+x[i-1]+x[i+1]);
//pnew[i]=.5*(dxf*rhsn+x[i-1]+x[i+1]);
pnew[i]=.5*(dxf3*(rhsn+rhsnnext+rhsnlast-drhsp*(2.0*x[i]-x[i+1]-pnew[i-1]))+pnew[i-1]+x[i+1]);


rhsnlast=rhsn;
rhsn=rhsnnext;
} // 
} //  enhanced 




if ( do_backward)
{
for ( IdxTy i=(sz-3); i>=2; --i)
{
const Value rhsn=c*(m.n[i]-m.p[i]+m.Nteff[i]);
nnew[i]=.5*(-dxf*rhsn+nnew[i-1]+nnew[i+1]);
pnew[i]=.5*(dxf*rhsn+pnew[i-1]+pnew[i+1]);
}
}





for ( IdxTy i=1; i<sz; ++i)
{

if ( nnew[i]>pnew[i])
{
nnew[i]=::exp(nnew[i]);
if ( nnew[i]!=0) pnew[i]=ni2/nnew[i]; else pnew[i]=ni;
}
else
{
pnew[i]=::exp(pnew[i]);
if ( pnew[i]!=0) nnew[i]=ni2/pnew[i]; else nnew[i]=ni;


}


}
// messing up ni2 constraint in some cases
const Value fold=0*flags.nfold;
m.n.mix(nnew,fold,1,sz-1);
m.p.mix(pnew,fold,1,sz-1);


} // update_2


static void update_carriers( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Value ni2=pr.ni2;
const Value ni=pr.ni;
const Grid & g=m.n.grid();
// shoudl assert uniform
const Value dx=g.min_dx();
const Value dxf=dx*dx;
const Value sz=g.size();
Values nnew(g),pnew(g),x(g),y(g);
x=m.p.ln();
y=m.n.ln();
const Value c=gp.q()/(gp.vt()*pr.epsilon*gp.epsilon_vacu());
// shoaoting?
Value test=y[1];
Value hi=test+.0000001;
Value lo=test-.001; 
IdxTy iix=0;
static IdxTy siter=0;
while ( true)
{
nnew[0]=y[0];
nnew[1]=test;
for ( IdxTy i=2; i<sz; ++i)
{
// this looks like it worked on the first pass, but then
// the nexttime all fudded up 
//const Value yi=y[i];
//const Value rhs=c*(::exp(yi)-ni2*::exp(-yi)+m.Nteff[i]);
//const Value rhs=c*(m.n[i]-m.p[i]+m.Nteff[i]);
const Value rhsn=c*(m.n[i-1]-m.p[i-1]+m.Nteff[i-1]);
//const Value ne=exp(nnew[i-1]);
//const Value rhsn=c*(ne-ni2/ne+m.Nteff[i-1]);
nnew[i]=dxf*rhsn+2.0*nnew[i-1]-nnew[i-2];
MM_MSG(" siter "<<siter<<" "<<iix<<" i= "<<i<<" n "<<nnew[i]<<" "<<(test-y[0])<<" "<<y[sz-1])
} //i 
const Value final=::exp(nnew[sz-1]);
// this can go inf.. 
//if ( final>m.n[sz-1]) hi=test; else lo=test;
// the search range is not monotonic lol 
if ( nnew[sz-1]>y[sz-1]) hi=test; else lo=test;


test=(hi+lo)/2.0;
MM_MSG(" "<<iix<<" final= "<<final<<" "<<nnew[sz-1]<<" tgt= "<<m.n[sz-1]<<" hi= "<<hi<<" lo= "<< lo << " test= "<<( test-y[0]) ) 
++iix;
//if ((hi/lo)<1.000000001) {
if ( iix>100) {
break;
}

} // true


for ( IdxTy i=1; i<sz; ++i)
{
nnew[i]=::exp(nnew[i]);
if ( nnew[i]!=0) pnew[i]=ni2/nnew[i]; else pnew[i]=ni;

}
const Value fold=flags.nfold;
m.n.mix(nnew,fold,1,sz-1);
m.p.mix(pnew,fold,1,sz-1);

++siter;

}


static void update_voltage( MatSt & m, const MatPr & pr, const GP & gp, const sim_flags & flags, JunkBin & junk )
{
const Value fold=flags.vfold; 
const Grid & g=m.n.grid();
const Value sz=g.size();
const Value eps=pr.epsilon*gp.epsilon_vacu();
Values ee(g),rhoeps(g),qn(g),vx(g);
Value qnet=0;
for ( IdxTy i=0; i<(sz); ++i) {
const Value Nt=m.Nteff[i]; 
const Value nq=Nt+m.n[i]-m.p[i];
qn[i]=nq;
// this starts out ok and drifts off by Nx times dx
qnet+=nq;
rhoeps[i]=nq*gp.q()/eps;
}
Told::I(ee,g,rhoeps);
Told::I(vx,g,ee);
//const Value fnew=1.0-fold;
// for ( IdxTy i=2; i<sz; ++i) m.v[i]=fold*v[i]+fnew*vx[i];
// probably want to return vx??? 
m.v.mix(vx,fold,2,sz);

} // update_voltage


class one_point_info
{
public:
one_point_info(const Value & aa, const Value & bb, const Value & cc, const Value & NN)
: m_a(aa),m_b(bb),m_c(cc),m_N(NN),m_x(0),m_y(0),m_en(0),m_ep(0),
m_dendx(0),m_dendy(0),m_depdx(0), m_depdy(0) {}

void set(const Value & x, const Value & y) { m_x=x; m_y=y; } 
Value e(const Value x) const { return ::exp(x); } 
Value rho() const { return m_a*(e(m_x)-e(m_y)-m_N );}
Value g(const Value & x, const Value & y) const
{ return m_b*(e(x)-m_c*(e(-y)));}
Value dd(const Value & x) const
{ return m_a*e(x)+m_b*m_c*(e(-x));}
Value dd2(const Value & x) const
{ return m_a*e(x)-m_b*m_c*(e(-x));}


Value en() const { return rho()-g(m_x,m_y);}
Value ep() const { return rho()+g(m_y,m_x);}
Value dendx() const { return (m_a-m_b)*e(m_x); } 
Value dendy() const { return -dd(m_y); } 
Value depdy() const { return -(m_a-m_b)*e(m_y); } 
Value depdx() const { return dd(m_x); } 



void update()
{
m_en=en();
m_ep=ep();
m_dendx=dendx();
m_dendy=dendy();
m_depdx=depdx();
m_depdy=depdy();
}
void solve_1()
{
const Value A=(-m_a+m_b)*e(m_x)/dd(m_x);
const Value B=(-m_a+m_b)*e(m_y)/dd(m_y);
//m_dy=-(m_en+A*m_ep)/(m_dendy*m_depdx-m_dendx*m_depdy);
m_dy=-(m_en+A*m_ep)/(m_dendy+A*m_depdy);
m_dx=-(m_ep+B*m_en)/(m_depdx+B*m_dendx);
const Value test_en=m_dx*m_dendx+m_dy*m_dendy;
const Value test_ep=m_dx*m_depdx+m_dy*m_depdy;
const Value rho=e(m_x)-e(m_y)-m_N;
MM_MSG(" comps dxy "<<m_dx<<" "<<m_dy<<" en= "<<m_en<<" "<<test_en<<" "<<(m_en+test_en)<<" ep= "<<m_ep<<" "<<test_ep<<" "<<(m_ep+test_ep)<<" rho "<<rho)



}
void update_dy_dx()
{
m_x+=m_dx;
m_y+=m_dy;
MM_MSG(" x = "<<m_x<<" "<<e(m_x)<<" m_y= "<<m_y<<" "<<e(m_y)<<" dd2 "<<dd2(m_x)<<" "<<(-dd2(m_y)))
}
void analytical()
{
const Value fx=.5; // (m_N<0)?(.5):(-.5);
const Value ddel=4.0*m_c/(m_N*m_N);
const Value rad=m_N*::sqrt(1.0+ddel);
const Value rada=::sqrt(m_N*m_N+4.0*m_c);
Value _y=-m_N/2.0+fx*rada;
Value dif=-m_N/2.0-.5*rad;
Value sum=-m_N/2.0+.5*rad;
Value _x=m_c/_y;
m_x_anal=::log(_x);
m_y_anal=::log(_y);
MM_MSG(" analytical rad= "<<rad<<" "<<(m_N/2.0)<<" "<<dif<<" "<<sum<<" x= "<<_x<<" "<<_y<<" del "<<(_x-m_N)<<" "<<(_y+m_N)<<" ddel "<<ddel<<" rho "<<(_x-_y-m_N))


}

Value m_a,m_b,m_c,m_N;
Value m_x,m_y;
Value m_en,m_ep;
Value m_dendx,m_dendy,m_depdx,m_depdy;
Value m_x_anal, m_y_anal;
Value m_dx,m_dy;

} ; // one_point_info


static void test_one_point( const Value&  a, const Value & b, const Value
& c, const Value & N  ,const IdxTy iter_count=4000
,const StrTy & load_from="",const StrTy & save_to="",const IdxTy flag=0)
{

one_point_info opi(a,b,c,N);
const Value div=1.0/(::log(::exp(1)));
//if N>0, these are ACCEPTORS and NEGATIVE 
// else N<0 POSITIVE DONOTRS
const Value x1=::log(::fabs(N))*div;
const Value y1=::log(c/x1)*div;
const Value x=(N>0)?x1:y1; // ::log(N)*div;
const Value y=(N>0)?y1:x1; // ::log(c/x)*div;
opi.analytical();
MM_MSG(" a= "<<a<<" b= "<<b<<" c= "<<c<<" N= "<<N<<" x= "<<x<<" y= "<<y<<" "<<opi.m_x_anal<<" "<<opi.m_y_anal)
opi.set(x,y);
for ( IdxTy i=0; i<iter_count; ++i)
{
opi.update();
opi.solve_1();
opi.update_dy_dx();
}


} // test_one_point




// try to get the debye tails for an abrupt junction
// using depletion approx as starting point
// note that "flag" is really the solver, but 
static void test_abrupt2( const IdxTy npts, const Value & dist, const Value
& Eext_, const Value tol ,const IdxTy iter_count=4000
,const StrTy & load_from="",const StrTy & save_to="",const IdxTy flag=0,
const IdxTy solverflag=0)
{
const bool load_materials_state=(load_from.length()!=0);
const bool save_materials_state=(save_to.length()!=0);
Grid g(0,dist,npts);
IdxTy sz=npts;
GeneralParameters gp;
material_properties si; si.Si();
material_state state(g);
const Value Nd=3e16;
const Value Na=3e16;
const Value junction=g[npts>>1];
// this work on quick check, current still bad at deplt
// see if this works better 
//const Value junction=g[npts-1];
make_abrupt(state,si,Nd,Na,junction,gp);
//make_one_side(state,si,Nd,Na,junction,gp);

sim_flags flags;
flags.solver=flag;
// return various diagnostics and intermediates j
JunkBin jb;
if ( load_materials_state) { state.load(load_from); }


while ( true)
{
const bool do_resample=flags.resample(jb.iter); 
if ( do_resample)
{
state.resample(); sz=state.n.size();
}
// these update the state, the old crap is then lost
// which makes combined things hard. For now just make
// copies, this is stupid but probably ok for now. 
update_voltage( state, si, gp,flags,  jb );

switch (flags.solver)
{
// various integrators to smooth differeence eqn
// slow convergence but probably ok 
case 0 : { update_carriers_2( state, si, gp,flags,  jb );
break; }

// linearizes and solves ssytems of equations generated
// by difference equations
// seems to converge fast but not sure about residual
// error and stability in general case 
case 1 : { update_carriers_3( state, si, gp,flags,  jb ); break; }

// Finds an energy and sort-of entropy and attempts to minize
// not work well 
case 2 : { 
update_carriers_grad_2( state, si, gp,flags,  jb );
update_voltage( state, si, gp,flags,  jb );
break; }
case 4 : { update_carriers_4( state, si, gp,flags,  jb,4,solverflag ); break; }
case 5 : { update_carriers_4( state, si, gp,flags,  jb,5,solverflag ); break; }
case 6 : { update_carriers_4( state, si, gp,flags,  jb,6,solverflag ); break; }
case 7 : { update_carriers_4( state, si, gp,flags,  jb,7,solverflag ); break; }
case 8 : { update_carriers_4( state, si, gp,flags,  jb,8,solverflag ); break; }


default : { {MM_MSG(" bad solver "<<flags.solver) } return; } 
};


// added solver flag to sim flags
if ( false )
{

// sofar this no longer seems to cycle through bizarre voltages
// with cleaned up size term to entropy and  and more careful
// size denomintors accounting for fixed ends
// however, it seems to inhbit converges and leaves tails. 
// bool br2=   ((jb.iter&63)!=0);
//  bool br2=   ((jb.iter&1023)!=0);
  bool br2= false; //   ((jb.iter&(1<<12))!=0);
//const IdxTy iter_swap=iter_count-20;
//bool br= (jb.iter<iter_swap); 
bool  br=!false;
//if ( br) update_carriers_2( state, si, gp,flags,  jb );
if ( br){
// smoothed intgeration
if ( br2) update_carriers_2( state, si, gp,flags,  jb );
//if ( br) 
else 
{ // solve difference equation system using blas
update_carriers_3( state, si, gp,flags,  jb );
}

//update_carriers_grad( state, si, gp,flags,  jb );
} // br
else{ // should be dun afterward... 
update_carriers_grad_2( state, si, gp,flags,  jb );
update_voltage( state, si, gp,flags,  jb );
}
} // false


measure( state, si, gp,flags,  jb );
if ( (jb.iter%255) == 0 ) {
MM_MSG(" iter= "<<jb.iter<<" sz= "<<sz<<" v= " <<state.v[sz-1]<<" dv= "<<(state.v[sz-1]-jb.last_v))
jb.last_v=state.v[sz-1]; 

}
//dump_crap(std::cout, state, si, gp,flags,  jb );
if (jb.iter>=iter_count) break;
if ( jb.terminate_code!=0) break; 
++jb.iter;

} // true

if ( save_materials_state) { state.save(save_to); }
dump_crap(std::cout, state, si, gp,flags,  jb );


/*
const bool dump_each=false;
// put the junction in the center, go back and forth 
Grid g(0,dist,npts);
// these are now changed due to resample
IdxTy junction=npts>>1;
IdxTy sz=npts;
g.check();
GP gp;
//Values v(g),vx(g), n(g),p(g),rho(g),n0g(g),p0g(g),nnew(g),pnew(g);
Values v(g),vx(g), n(g),p(g),nnew(g),pnew(g),Nteff(g);
//Values e(g),gradj(g),gradn(g);
//Values temp(g);
//Value Eext=Eext_;
// put Na or p side on the right 
//const Value Na=1e16;
// forgot dielectic const, lol
const Value Na=3e16;
//const Value Nd=3e15;
const Value Nd=Na;
const Value kT=gp.vt();
const Value eps=gp.epsilon();
const Value Nc=2.8e19;
const Value Nv=1.83e19;
const Value ni=1.4e10;
const Value ni2=ni*ni;
const Value Ev=0;
const Value egap=1.1;
const Value Ec=Ev+egap;
// the efn is measured from confuction band
// wtf, jkust make them both from Ec and figure it out later lol
//const Value efn=kT*::log(Nd/Nc);
// relative to Ec...
//const Value efn=Ev-Ec+fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efn=fermi_level( 0, 0,  Nd,  Ec-.001,  Nc,  Nv,  kT , Ec,Ev);
const Value efp=fermi_level( Na, Ev+.001,  0,0 ,  Nc,  Nv,  kT , Ec,Ev);
// from valuence band 
//const Value efp=kT*::log(Na/Nv);
//const Value vlim=egap;
//const Value Ev=-egap;
//const Value builtin=::fabs(efn+efp+egap);
// the efn is now refernced to Ec just like the efp so just subtact 
// this is positive, as efn>efp 
const Value builtin=(efn-efp); // it looks forward biases lol 

//const Value phi_target=0; // applied voltage measured

Value phi=0; // applied voltage measured
// should assert uniform grid
// but now we resample 
Value dx=g.min_dx();
Value dxf=dx*dx/eps*gp.q();
// try to guess the depltion region
const Value xp=::sqrt(2.0*builtin*Nd/(Na*(Na+Nd))*eps/gp.q());
const Value xn=::sqrt(2.0*builtin*Na/(Nd*(Na+Nd))*eps/gp.q());
Value Ezed= .5*gp.q()*Na*xp/eps;
 Ezed+= .5*gp.q()*Nd*xn/eps;
//Value Ezed_max=Ezed*10;
//Value Ezed_min=-Ezed_max;
MM_MSG("builtin= "<<builtin<<" efp= "<<efp<<" efn= "<<efn<<" xp= "<<xp<<" xn= "<<xn)
//const Value phitol=::fabs(builtin)*1e-7;
IdxTy cnt=0;
const Value Nd2=Nd/2;
const Value Na2=Na/2;
const Value Ndeff=Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Npdeff=-Nd2+::sqrt(Nd2*Nd2+ni2);
const Value Naeff=Na2+::sqrt(Na2*Na2+ni2);
const Value Nnaeff=-Na2+::sqrt(Na2*Na2+ni2);
for ( IdxTy i=0; i<(sz); ++i)
{
const bool left=(i<junction);
// this is negaive of "q"
Value Nt=(left)?(-Nd):(Na);
Nteff[i]=Nt;
}

// this integrates out to about the builtin voltage but it
// seems to not sum the charge exactly to zero so the end drifts
// back but this goes away with finer grid, 
Value qini=0;
bool was_depleted=false;
for ( IdxTy i=0; i<sz; ++i)
{
const bool left=(i<junction);
const Value dn=g[junction]-xn-g[i];
const Value dp=g[junction]+xp-g[i];
const bool depleted= ( left)?(g[i]>(g[junction]-xn)):(g[i]<(g[junction]+xp));
if ( !depleted){
if ( left) { n[i]=Ndeff; p[i]=Npdeff;v[i]=0; qini+=n[i]-p[i]-Nd;   }
else { p[i]=Naeff; n[i]=Nnaeff; v[i]=-builtin; qini+=n[i]-p[i]+Na;}
}

if ( depleted) { n[i]=ni; p[i]=ni; 
if ( left  ) {v[i]=-.5*gp.q()*Nd/eps*dn*dn; qini-=Nd; } 
else { v[i]=-builtin+.5*gp.q()*Na/eps*dp*dp; qini+=Na; } 
} // depleted
// dump the charge here for now 
if ( !depleted&&was_depleted)
{
// simplyusing these values kept charge at 1e12 instead of creeping
// to 1e16 in a few iters. 
const Value Na2x=(Na+3e-13*qini)/2.0;
const Value Naeffx=Na2x+::sqrt(Na2x*Na2x+ni2);
const Value Nnaeffx=-Na2x+::sqrt(Na2x*Na2x+ni2);
p[i-1]=Naeffx; n[i-1]=Nnaeffx;
//p[i]+=qini;
if (false) if ( p[i]!=0)
{
Value qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;
qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;
qxx=n[i]-ni2/p[i];
n[i]=ni2/p[i];
p[i]+=qxx;

}


//n[i]-=0*qini/2;
} 
was_depleted=depleted;
} //i init
// if this did not achieve neutrality, dump more charge somewhere. 
const bool some_resamples=!true;
IdxTy rescnt=0;
Grid * gnew=NULL;
while ( true) // phi-rho loop 
{
const bool do_resample=some_resamples&&(rescnt<6)&&((cnt%1000)==999);
if ( do_resample)
{
const IdxTy sznew=g.size()*2;
MM_MSG(" begin resample sz="<<sznew<<" "<<g.size())
//Grid g(0,dist,npts);
Grid * gold= gnew;
gnew= new  Grid(0,dist,sznew);
MM_TRACE
Grid & gn= * gnew;
// we can not give ctor address of a temp lol 
// except now the next time through, wtf... 
Values xx(gn);
xx=v.resample2x(gn); v=xx;
xx=vx.resample2x(gn); vx=xx;
xx=p.resample2x(gn); p=xx;
xx=pnew.resample2x(gn); pnew=xx;
xx=n.resample2x(gn); n=xx;
xx=nnew.resample2x(gn); nnew=xx;
xx=Nteff.resample2x(gn); Nteff=xx;
junction=junction*2;
sz=gn.size();
g=gn;
dx=g.min_dx();
dxf=dx*dx/eps*gp.q();
// try to guess the depltion region
++rescnt;
delete gold;
MM_MSG(" done resampl "<<rescnt)
} // resamples

v[0]=0;
v[1]=0;
Value qtot=0;
// this does not gie good asymptotes because total charge is not zero.
// that should occur automatically. 
for ( IdxTy i=0; i<sz; ++i) vx[i]=v[i];

const bool udate_finite_diff=false;
if ( udate_finite_diff)
{
for ( IdxTy i=2; i<(sz-1); ++i)
//for ( IdxTy i=1; i<(sz-2); ++i)
{
// this code is notusednow anyway... 
//const bool left=(i<junction);
//Value Nt=(left)?(-Nd):(Na);
Value Nt=Nteff[i];
const Value nq=Nt+n[i]-p[i];
vx[i]=.5*(v[i-1]+v[i+1]-dxf*(nq));
//vx[i+1]=2*v[i]-v[i-1]+dxf*(nq);

qtot+=nq;
if ( dump_each) MM_MSG(cnt<<" "<<left<<" depl "<<i<<" "<<g[i]<<" v= "<<v[i]<<" n= "<<n[i]<<" p= "<<p[i]<<" n+p= "<<(n[i]+p[i])<<" "<<(n[i]*p[i]/ni2)<<" qtot= "<<qtot<<" nq= "<<nq)
// update n and p for zero J constraint 
}
} // update+finite
else
{
const bool force_voltage=!true;
const bool force_voltage_bound=!true;
const bool force_neutral=!true;
Values ee(g),rho(g),qn(g);
Value qnet=0;
for ( IdxTy i=0; i<(sz); ++i) {
//const bool left=(i<junction);
//Value Nt=(left)?(-Nd):(Na);
const Value Nt=Nteff; 
const Value nq=Nt+n[i]-p[i];
qn[i]=nq;
// this starts out ok and drifts off by Nx times dx
qnet+=nq;
rho[i]=nq*gp.q()/eps;
}
if ( force_neutral)
{
Value qnet2=0;
const Value nq=-qnet/sz;
//MM_MSG("ASSFUDD "<<nq<<" qnet= "<<qnet)
for ( IdxTy i=0; i<(sz); ++i) {
if ( nq<0) p[i]-=nq; else n[i]+=nq;
rho[i]+=nq*gp.q()/eps;
qnet2+=rho[i];
}
//MM_MSG(" qnet was "<<qnet<<" now is "<<((eps/gp.q())*qnet2))
} // force_neutral

I(ee,g,rho);
I(vx,g,ee);
// boundary conditions should create field of zer at both ends..
// if not there is a net charge
Value a=vx[0];
Value b=(vx[sz-1]-vx[0]+builtin)/(sz-1.0);
if ( (cnt&255) == 0 ) 
{MM_MSG(cnt<<" qnet= "<<qnet<<" eerror= "<<b<<" "<<vx[sz-1])}

if (force_voltage) for ( IdxTy i=0; i<(sz); ++i) { vx[i]-=a+b*i; } 
if (force_voltage_bound) { vx[sz-1]=-builtin; vx[sz-2]=vx[sz-1];} 


//vx[0]=v[0]; vx[1]=v[1];
//vx[0]=v[0]; vx[1]=v[1];
// this is slower than the computations, not surprising
if (false) for ( IdxTy i=0; i<(sz); ++i) {
MM_MSG(" "<<cnt<<" dumpin i= "<<i<<" ee= "<<ee[i]<<"  rho= "<<rho[i]<<" vx= "<<vx[i]<<" qnet= "<<qn[i])
}


}
const bool do_backwards_too=!true;
if ( do_backwards_too)
{
for ( IdxTy i=(sz-2); i>1; --i)
{
const bool left=(i<junction);
Value Nt=(left)?(-Nd):(Na);
const Value nq=Nt+n[i]-p[i];
vx[i]=.5*vx[i]+.5*(.5*(v[i-1]+v[i+1]-dxf*(nq)));

} // i 

} // backwards

const Value fold=.1;
const Value fnew=1.0-fold;
for ( IdxTy i=2; i<sz; ++i) v[i]=fold*v[i]+fnew*vx[i];

const bool fixed_current=!true; 
const bool fixed_current_too=!true; 
const bool use_equ=false; // !fixed_current; 
const bool fixed_x=false; // !fixed_current; 
const bool use_srh=false; // !fixed_current; 
const bool use_sum_and_diff=!false; // !fixed_current; 
if ( use_sum_and_diff){

const Value G=1e14;
const Value mu=gp.mu();
const Value q=gp.q();

 //sum_and_diff(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
// sum_and_diff2(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
 //sum_and_diff3(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0.0,Ndeff*mu);
 // sum_and_diff4(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0,Ndeff*mu);
 sum_and_diff5(nnew,pnew,v,n,p,ni2,G,mu,dx,sz,kT,q,0,Ndeff*mu, Nteff);
}

if ( fixed_x)
{
//const Value G=400000000;
// negative values explode, positive sign is right AFAICT
// 1e2 with the finit diff, 1e13 with the integrator lol 
//const Value G=1e2;
//const Value G=1e13;
//const Value G=1e15;
const Value G=1e12;
//const Value G=00;
const Value mu=gp.mu();
const Value q=gp.q();
Values rhodiveps(g),e(g);
const bool need_garbage=false;
if ( need_garbage)
{
D(e,g,v); // sign is wrong here. 
for ( IdxTy i=0; i<sz; ++i) {
const bool left=(i<junction);
Value Nt=(left)?(-Nd):(Na);
rhodiveps[i]=-q*(Nt+n[i]-p[i])/gp.epsilon();
}
} // garbage
// q negative"? 
fixed_j_total(nnew,pnew,v,n,p,rhodiveps,e,ni2,G,mu,dx,sz,kT,q);
}

if ( use_srh)
{
const Value G=1e17+0*1e2;
const Value Gsrh=1e18;
const Value mu=gp.mu();
const Value q=gp.q();
Values junk;

srh_total(nnew,pnew,v,n,p,junk,junk,ni2,G,Gsrh,mu,dx,sz,kT,q);
} // srh 

// this causes junction creep 
if ( fixed_current) fixed_j(v,n,p,nnew,pnew,gp,sz,kT);
// efn and efp referenced to Ec and Ev, 
if ( use_equ ) equill(v,n,p,nnew,pnew,gp,sz,kT,junction,Na,Nd,Nc,Nv,efn,efp,egap);
//Value d1=0;
//Value d2=0;
//for ( IdxTy i=0; i<(sz); ++i) d1+=(n[i]-nnew[i])*(n[i]-nnew[i]);

for ( IdxTy i=2; i<(sz-1); ++i)
{ n[i]=(fold*n[i]+fnew*nnew[i]); p[i]=(fold*p[i]+fnew*pnew[i]); } 


if ( fixed_current_too) { fixed_j(v,n,p,nnew,pnew,gp,sz,kT);
for ( IdxTy i=2; i<(sz-1); ++i)
{ n[i]=(.75*n[i]+.25*nnew[i]); 
p[i]=(.75*p[i]+.25*pnew[i]);
//if ( n[i]>0) p[i]=ni2/n[i];
 } 
} 

++cnt;
phi=v[0]-v[sz-1];
//const Value dphi=phi-builtin;
const bool converged=true; // (::fabs(dphi)<phitol);
if ( (cnt%31) ==0 ) MM_MSG("cnt "<<cnt<<" phi "<<phi<<" v0= "<<v[0]<<" vd= "<<v[sz-1]<<" builtin= "<<builtin)

if (converged&&(cnt>iter_count))
{ 
//static void D(Values & d, const Grid & g, const Values & n)
Values e(g), grn(g), grp(g);
D(e,g,v); D(grn,g,n); D(grp,g,p);
const Value f1=gp.q()*gp.mu();
for ( IdxTy i=0; i<sz; ++i){
const Value Jn=f1*(-n[i]*e[i]+kT*grn[i]);
const Value Jp=f1*(-p[i]*e[i]-kT*grp[i]);
MM_MSG("cnt= "<<cnt<<" i= " <<i<<" "<<g[i]<< " phi= "<<v[i]<<" n= "<<n[i]<<" p= "<<p[i]<<" Jp= "<<(Jp)<<" Jn= "<<(Jn)<<" Jt= "<<(Jn+Jp)<<" ndrift= "<<(-f1*n[i]*e[i])<<" "<<(f1*kT*grn[i])) 

}
break;
}
} // phi-rho loop 
delete gnew;
*/

} // test_abrupt

}; // dummy2 ?

}; // ns  mjm_dd_base


#endif


