#ifndef MJM_SPECIES_H__
#define MJM_SPECIES_H__

#include <string>
#include "mjm_globals.h"

class mjm_species
{

typedef double R;
typedef std::string StrTy;
public:
mjm_species(): q(1),mu(10), qunits(1.6e-19) {}
mjm_species(const R & _q, const R & _mu): q(_q),mu(_mu), qunits(1.6e-19) {}

// this is flux or particle current density not charge
// using mu for D may be dumb here
// note this iactually kT/q or .0259 
R phidd(const R& n, const R& gradn, const R& E, const R& kT ) const
{

if ( q==0) return -mu*kT*gradn;
// the firld points in dir of positive charge force
if ( q>0) return mu*n*E-mu*kT*gradn;
return -mu*n*E-mu*kT*gradn;


}

R Jdd(const R& n, const R& gradn, const R& E, const R& kT ) const
{
// electrons use plus sign 
if ( q<0) return qunits*q*mu*(n*E+kT*gradn);
return qunits*q*mu*(n*E-kT*gradn);

}


const StrTy  name( const int nidx=0, const StrTy & sep=" ") const
{

switch (nidx)
{
case 0: return  m_name;
case 1: return  m_variable_name;
case 2: return  m_chem_name;
case 3: return  m_notes;


}

return name(0)+sep+name(1)+sep+name(2)+sep+name(3);

}

template<class Os> Os& dump(Os & os, const int nidx=0) const
{
	const StrTy sep=" ";
	os<<name(nidx)<<sep;
	os<<"q= "<<q<<sep;
	os<<"mu= "<<mu<<sep;
//	os<<"qunits= "<<qunits<<sep;

}


R q,mu,qunits;

private:

StrTy m_name;
StrTy m_variable_name;
StrTy m_chem_name;
StrTy m_notes;


}; // mjm_species



#endif
