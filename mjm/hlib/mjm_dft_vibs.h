#ifndef MJM_TCL_DFT_VIBS_BASE_H__
#define MJM_TCL_DFT_VIBS_BASE_H__


#include <stdlib.h>
//#include <tcl.h>
//#include "mjm_tcl_base.h"
#include "mjm_solvers.h"
#include "mjm_jdftx_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
// this should be moved and unified.
//#include <cblas.h>
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>

//#include "mjm_sample_spectrum.h"
//#include "mjm_io.h"
// the crap that needs this shold be moved.
#include "prob_current.h"
#include "mjm_interchanges.h"



class pert_class : public mjm_jdftx_typedefs
{
public:
pert_class( ) {  } 
template <class IsTy> pert_class(IsTy * _ions, IsTy *_fis ) {
ions.load(_ions); 
 forces.load(_fis);
// zomewhat of a kluge?
// no obvious effect
//forces.zero_means(); 
 } 
mjm_jdftx_ionpos ions;
mjm_jdftx_forces forces;

};


class perts_class : public  std::vector<pert_class> , mjm_jdftx_typedefs
{
typedef std::map<StrTy,D> Masses;
typedef mjm_jdftx_ionpos IonsTy;
typedef IonsTy::iterator IonItor;
typedef mjm_jdftx_forces ForcesTy;
typedef mjm_2D_data Data2DTy;
typedef mjm_2D_data_uint CountMatTy;

Masses m_mass_map;

void pput(const char * fuch, D wtf)
{
m_mass_map[StrTy(fuch)]=wtf;
}

void load_default_masses()
{
pput("h",1.0078);
pput("c",12);
pput("o",16);
pput("fe",55);

}

D mass(const IdxTy atom, const IonsTy  & ions)
{
const StrTy & species=ions.get_label(atom,1);
return m_mass_map[species];

}
public:
perts_class()
{

load_default_masses();

}
void geometry_matrix(Data2DTy & dm, const IonsTy & ionx)
{
IonItor ii(ionx),jj(ionx);
while ( ii) { 
while ( jj) { 
D x=ionx.distance(ii,jj);
info()<<MM_MARK<< " in geo  "<<*ii<<" "<<*jj<<" "<<x<<" "<<CRLF;
if (x!=0)
{
D dalpha=(ionx.get(ii.atom(),ii.coord())-ionx.get(jj.atom(),ii.coord()));
D dbeta=(ionx.get(ii.atom(),jj.coord())-ionx.get(jj.atom(),jj.coord()));
D d3=x*x*x;
if (ii.coord()==jj.coord()) x=(1.0-3.0*dalpha*dbeta/x/x)/d3;
else x=-3.0*dalpha*dbeta/d3/x/x;
dm(ii,jj)=x;
}
 ++jj; }
jj.clear();
++ii;
}



}
IdxTy index_zeff(const IonItor & ii, const IonItor & jj, const IdxTy sz)
{
IdxTy i=ii.atom();
IdxTy j=jj.atom();
// for now, just do the simple thing and make too many
if ( i<j) return  i*sz+j;
return j*sz+i; 

}
template <class Tx, class Ty > void zero( Ty  d, const Tx nsz ) 
{ ::memset(d,0,nsz); } 
void find_zeff(D* zpairs, const Data2DTy & c, const Data2DTy & a, const IdxTy nzeff2, const IdxTy sz)
{
info()<<MM_MARK<< " finding zeff  "<<nzeff2<<" "<<CRLF;
D num[nzeff2],den[nzeff2];
//IdxTy cnt[nzeff2];
zero(&num[0],nzeff2*sizeof(D)); 
zero(&den[0],nzeff2*sizeof(D)); 
//zero(&cnt[0], nzeff2*sizeof(int));
IonItor ii(sz),jj(sz);
while ( ii) { 
jj.clear();
while ( jj) { 
IdxTy idx=index_zeff(ii,jj,sz);
info()<<MM_MARK<< " finding zeff idx "<<idx<<" "<<*ii<<" "<<*jj<<" "<<CRLF;
num[idx]+=a(ii,jj)*c(ii,jj);
den[idx]+=a(ii,jj)*a(ii,jj);
++jj;
}
++ii;
}
// this leaves unused slots but it makes the index calculations
// simpler for now
for (IdxTy i=0; i<nzeff2; ++i) 
if (den[i]!=0) zpairs[i]=0.0-1.0*num[i]/den[i];
else zpairs[i]=0;


}


void kluge_matrix( Data2DTy & dest, const IonsTy & ionx)
{
IonItor ii(ionx),jj(ionx);
while ( ii) { 
jj.clear();
while ( jj) { 
if (ii.atom()!=jj.atom())
{
//const IdxTy idx=index_zeff(ii,jj,szi);
//info()<<MM_MARK<< " finding zeff idx "<<idx<<" "<<*ii<<" "<<*jj<<" "<<CRLF;
//dest(*ii,*jj)=dm(*ii,*jj)*Zeff[idx];
if(*ii>*jj)
{
dest(*ii,*jj)=.5*(dest(*ii,*jj)+dest(*jj,*ii));
dest(*jj,*ii)=dest(*ii,*jj);
}

}
++jj;
}
++ii;
}

info()<<MM_MARK<< " dun kluging  "<<" "<<CRLF;
}


/*
void make_masses( Data2DTy & dest, const IonsTy & ionx )
{
typedef IonsTy::iterator Itor;
Itor ii(ionx),jj(ionx);
IdxTy szi=ionx.size();
IdxTy modes=szi*3;
dest.set_size(modes,modes);
while ( ii)
{



}
*/

void clean_matrix( Data2DTy & dest, const IonsTy & ionx)
{
const IdxTy szi=ionx.size();
const IdxTy modes=szi*3;
Data2DTy dm;
dm.set_size(modes,modes);
// set up the relative loations matrix
geometry_matrix(dm,ionx);
info()<<MM_MARK<< " made it past geo mat "<<" "<<CRLF;
//IdxTy nzeff2=szi*(szi-1)/2;
// to make the index generator simpler, just make too manuy
IdxTy nzeff2=szi*(szi-0);
// for off diag terms, find teh Z to minimize errors, then
// reset terms to match
D Zeff[nzeff2];
find_zeff(Zeff, dest,dm,nzeff2,szi);
info()<<MM_MARK<< " zeff ios now "<<" ";
for (IdxTy fu =0 ;fu<nzeff2; ++fu) info()<<Zeff[fu]<<"~"<<(::sqrt(::fabs(Zeff[fu])/4.0/3.14159))<<" ";
info()<<CRLF;

// we should be able to reuse these but they are cheap 
IonItor ii(ionx),jj(ionx);
while ( ii) { 
jj.clear();
while ( jj) { 
if (ii.atom()!=jj.atom())
{
const IdxTy idx=index_zeff(ii,jj,szi);
info()<<MM_MARK<< " finding zeff idx "<<idx<<" "<<*ii<<" "<<*jj<<" "<<CRLF;
dest(*ii,*jj)=dm(*ii,*jj)*Zeff[idx];
}
++jj;
}
++ii;
}

info()<<MM_MARK<< " dun cleaning  "<<" "<<CRLF;
}
void make_masses( Data2DTy & dest, const IonsTy & ionx )
{
typedef IonsTy::iterator Itor;
Itor ii(ionx),jj(ionx);
IdxTy szi=ionx.size();
IdxTy modes=szi*3;
dest.set_size(modes,modes);
while ( ii)
{
D m1=mass(ii.atom(),ionx);
while (jj)
{
D m2=mass(jj.atom(),ionx);
dest(ii,jj)=1.0/::sqrt(m2*m1);
++jj;
}
jj.clear();
++ii;
}

}

void find_perturbations_non_cart( Data2DTy & dest, const IonsTy & eq_ions, 
const ForcesTy  & eq_forces, const StrTy & param)
{
//const bool doing_cart=(param!="natural");
IdxTy szi=eq_ions.size();
const IdxTy modes=szi*3; // leave in translation rotation etc. 
dest.set_size(modes,modes);
CountMatTy counts(modes,modes);
Data2DTy prior(modes,modes),fprior(modes,modes),dprior(modes,modes);
typedef IonsTy::iterator Iitor;
Iitor ii(eq_ions);
const IdxTy persz=this->size();
// for each perturbtion
for ( IdxTy i=0; i<persz; ++i)
{
// get ions and forces for that displacement
mjm_jdftx_ionpos & iox=(*this)[i].ions;
mjm_jdftx_forces & fx =(*this)[i].forces;
// find out which atom moved, this shoudl be more general
IdxTy thisatom=0;
D dx=0;

for ( IdxTy ai=0; ai<szi; ++ai)
{
D dthis=eq_ions.distance2(ai,iox,ai);
if ( dthis>dx) { thisatom=ai; dx=dthis; } 

} // ai atom itor 
IdxTy thiscoord=eq_ions.get_coord(thisatom,iox);
IdxTy mode=thisatom*3+thiscoord;
ii.clear();
while (ii)
{
const D bg_force=eq_forces.get(ii.atom(),ii.coord());
D res=0;
// fc,dist,2bg,3observed  force dot bg force, 4observer dot dir, 5bg dot dir
enum { FC=0, DIST=1, BG=2, F_DOT_BG=3, RAW_DOT_DIR=4, BG_DOT_DIR=5};
std::vector<D> rev=eq_ions.resolve(thisatom,thiscoord,ii.atom(),ii.coord(),iox,fx,eq_forces);
res=rev[FC];
IdxTy  mcnt= counts(*ii,mode);
if ( mcnt==0) { dest(*ii,mode)=res; prior(*ii,mode)=rev[F_DOT_BG];
fprior(*ii,mode)=rev[RAW_DOT_DIR]; dprior(*ii,mode)=rev[DIST]; 
} 
else{
const D good=.8;
const D bad =1.0- good;
const D dthis=rev[DIST]; 
const D dold=dprior(*ii,mode);
// these are distnaces
//const D dthis2=dthis*dthis;
//const D dold2=dold*dold;
// this is crap
//const D dthis3=dthis*dthis*dthis;
//const D dold3=dold*dold*dold;
// this is crap
//const D coef_b=1.-(dthis3+dold3)/(dthis3-dold3);
//D rhs=(rev[RAW_DOT_DIR]+fprior(*ii,mode)-2.0*rev[BG_DOT_DIR])/(dthis-dold);
//rhs-=(rev[RAW_DOT_DIR]-fprior(*ii,mode))*(dthis3+dold3)/((dthis3-dold3)*(dthis-dold));
//const D bb=rhs/coef_b;
// different appraoch
// ok, this is the x^3 coef LOL
const D tforce=(rev[RAW_DOT_DIR]-rev[BG_DOT_DIR]);
const D oforce=(fprior(*ii,mode)-rev[BG_DOT_DIR]);
//const D bb3= (dthis*(oforce) -dold*(tforce))/(dthis*dold3-dold*dthis3);

//const D bb=(oforce-bb3*dthis3)/dthis*30.0/36*2.0/3.0*3.0/4.0;
//const D bb=(oforce-bb3*dthis3)/dthis*30.0/51.0*30.0/51.0;
//const bool bb_usable=false; // (::fabs(dthis-dold)>.011);
const D exp=-2;
const D p1=::pow(dthis,exp); // *::sqrt(::fabs(dthis));
const D p2=::pow(dold,exp); // *::sqrt(::fabs(dold));

//const D bb=(dthis2*(oforce)-dold2*tforce)/(dthis2*dold-dold2*dthis); 
const D bb=(p1*(oforce)-p2*tforce)/(p1*dold-p2*dthis); 


const bool bb_usable=true; // (::fabs(dthis-dold)>.011);
//const bool bb_usable=false; // (::fabs(dthis-dold)>.011);
if (bb_usable) dest(*ii,mode)=bb;
else
{
// add more weight if they point opposite
if ((rev[3]<0) &&(prior(*ii,mode)>0)) 
{ dest(*ii,mode)=good*dest(*ii,mode)+bad*res; }
else if ((rev[3]>0) &&(prior(*ii,mode)<0)) 
{ dest(*ii,mode)=bad*dest(*ii,mode)+good*res; }
else
 dest(*ii,mode)=.5*(dest(*ii,mode)+res);
}
// may help slightly
// if ( ::fabs(dest(*ii,mode))<1e-3) dest(*ii,mode)=0;
}
 
++counts(*ii,mode);

D fom=res/bg_force;
if (0)if ( ::fabs(fom)<100)
{ dest(*ii,mode)=0;
info()<<MM_MARK<<" DANGER WILL ROBINSN SETTING COEF TO ZERO KLUGE KLUG  "<<CRLF;
}
//if ( dest(mode,j)*dest(mode,j)<1e-5) dest(mode,j)=0;
info()<<MM_MARK<<" FORCES RESOLUT  atom="<<ii.atom()<<" coord="<<ii.coord()<<" k="<<res<<" "
<<dest(*ii,mode)<<" bg="<<bg_force<<" f="<<fx.get(ii.atom(),ii.coord())
<<" pct="<<fom
<<" j="<<*ii<<" mode="<<mode<<CRLF;
++ii;
} // ii is mode


} // i


}



void find_perturbations( Data2DTy & dest, const IonsTy & eq_ions, 
const ForcesTy  & eq_forces, const StrTy & param)
{
const bool doing_cart=(param!="natural");
if ( !doing_cart)
{
find_perturbations_non_cart( dest,eq_ions,eq_forces, param);
return;
}
IdxTy szi=eq_ions.size();
IdxTy ffi=eq_forces.size();
if ( szi!=ffi) {
 info()<<MM_MARK<<"ERROR ion size is "<<szi<<" not eq forces "<<ffi<<CRLF;
if ( ffi<szi) szi=ffi;
}
const IdxTy persz=this->size();
const IdxTy modes=szi*3; // leave in translation rotation etc. 
dest.set_size(modes,modes);
for ( IdxTy i=0; i<persz; ++i)
{

const IdxTy atom_in_theory=i/3;
//const IdxTy coord_in_theory=i%3;
mjm_jdftx_ionpos & iox=(*this)[i].ions;
mjm_jdftx_forces & fx =(*this)[i].forces;
if ( iox.size()!=eq_ions.size())
{ 
info()<<MM_MARK<<" skipping entry due to ion size problem "<<iox.size()<<CRLF; 
continue; 
}
if ( fx.size()!=eq_forces.size())
{ 
info()<<MM_MARK<<" skipping entry due to force size problem "<<fx.size()<<CRLF; 
continue; 
}


D dmax=0;
D dmaxv=0;
//D totd=0;
IdxTy atom=0;
IdxTy coord=0;
for ( IdxTy ai=0; ai<szi; ++ai)
{

for ( IdxTy ci=0; ci<3; ++ci)
{
D dmv=iox.get(ai,ci)-eq_ions.get(ai,ci);
D dm=dmv*dmv;
if ( dm>dmax) {
//info()<<MM_MARK<<" setting now "<<ai<<" "<<ci<<" "<<dmv<<CRLF; 
 dmax=dm; dmaxv=dmv; atom=ai; coord=ci; } 
}
} // ai
IdxTy coord_cartesian=coord;
if  ( !doing_cart) coord=eq_ions.get_coord(atom,iox);

 IdxTy mode=atom*3+coord;
// these can be missing s things fill in, should just keep track of max etc.
if ( mode!=i) info()<<MM_MARK<<" the best fit is not where it should be "
<< " assuming naural coodes "<<CRLF;
if ( atom!=atom_in_theory)  info()<<MM_MARK<<" atoms not match "<<
atom<<" implied by i is "<<atom_in_theory<<CRLF;

//mode=i ; // 
//info()<<MM_MARK<< " !!! setting mode to sequence not empirical for natural " <<CRLF;
//totd=iox.distance(atom,eq_ions,atom);

info()<<MM_MARK<<" dmax "<<dmaxv<<" at "<<atom<<" "<<coord<< "(was "<<coord_cartesian<<")"<<" mode="<<mode<<CRLF;
const IdxTy thisatom=atom;
const IdxTy thiscoord=coord;
atom=0;
coord=0;
for ( IdxTy j=0; j<modes; ++j)
{
// need the masses at some point. 
// by mobing this one atom, these are the forces each other
// atom, the mass then is the other atom mass.
//dest(j,mode)=(fx.get(atom,coord)-eq_forces.get(atom,coord))/(dmaxv)/mass(atom,iox);
const D bg_force=eq_forces.get(atom,coord);
//dest(mode,j)=(fx.get(atom,coord)-eq_forces.get(atom,coord))/(dmaxv)/mass(atom,iox);
D res=0;
//if (coord==thiscoord) 
//res=(fx.get(atom,coord)-bg_force)/(dmaxv)/(D)mass(atom,iox);
//D m1=mass(atom,iox);
//D m2=mass(thisatom,iox);
//if ( atom==0) m1=2;
//if ( thisatom==0) m2=2;
//res=(fx.get(atom,coord)-bg_force)/(dmaxv)/::sqrt(m1*m1);

if ( doing_cart) res=(fx.get(atom,coord)-bg_force)/(dmaxv); // /::sqrt(m1*m1);
else res=eq_ions.resolve(thisatom,thiscoord,atom,coord,iox,fx,eq_forces)[0];

// eh, tjhis needs to find the direction now and the force too

//res=(fx.get(atom,coord)-bg_force)/(totd); // /::sqrt(m1*m1);


//info()<<MM_MARK<<" MASSES "<<mass(atom,iox)<<" "<<mass(thisatom,iox)<<CRLF;
//res=(fx.get(atom,coord)-bg_force)/(dmaxv)/::sqrt((D)12);
dest(j,mode)=res;
D fom=res/bg_force;
if (0)if ( ::fabs(fom)<100)
{ dest(j,mode)=0;

info()<<MM_MARK<<" DANGER WILL ROBINSN SETTING COEF TO ZERO KLUGE KLUG  "<<CRLF;
}
//if ( dest(mode,j)*dest(mode,j)<1e-5) dest(mode,j)=0;
info()<<MM_MARK<<" FORCES RESOLUT  atom="<<atom<<" coord="<<coord<<" k="<<res<<" "
<<dest(j,mode)<<" bg="<<bg_force<<" f="<<fx.get(atom,coord)
<<" pct="<<fom
<<" j="<<j<<" mode="<<mode<<CRLF;
//info()<<MM_MARK<<" mass "<<mass(atom,iox)<<" "<<iox.species(atom)<<CRLF;
++coord;
if ( coord==3) { coord=0; ++atom; } 
}


}


}

};

#if 0 
};
/*
template <class SigClass> class mjm_tcl_base


*/
#endif
#endif
