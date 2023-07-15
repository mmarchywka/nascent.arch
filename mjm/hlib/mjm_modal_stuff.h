#ifndef MJM_MODAL_STUFF_H__
#define MJM_MODAL_STUFF_H__


#include <stdlib.h>
#include <tcl.h>
//#include "mjm_tcl_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
//#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

// #include <cblas.h>
#include "mjm_conflicts.h"
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

//#include "mjm_sample_spectrum.h"
#include "mjm_io.h"
//#include "prob_current.h"
#include "mjm_interchanges.h"

extern "C"
{
// eigenvalue solver, real coeefs
void dgeev_(char * JOBVL, char * JOBVR, int * N, double * A,
int * LDA, double * WR, double * WL, double * VL,int * LDVL,
double * VR, int * LDVR, double * WORK, int * LWORK, int * INFO); 
// want to add svd for det finding etc

};


class mjm_modal_typedefs
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;


static OsTy & info() { return std::cout; } 

};

class mjm_modal_mode : public mjm_modal_typedefs
{




};


class mjm_modal_stuff : public mjm_modal_typedefs
{
protected:
typedef std::vector<StrTy> Labels;
typedef std::vector<D> Data;
typedef std::pair<Labels,Data> Line;
typedef std::vector<Line> Lines;

Lines m_data;
public:
IdxTy size() const { return m_data.size(); } 
const D & get(const IdxTy atom, const IdxTy coord) const
{ return m_data[atom].second[coord]; }
const StrTy & get_label(const IdxTy atom, const IdxTy coord) const
{ return m_data[atom].first[coord]; }

void perturb(const IdxTy atom, const IdxTy coord, const D delta)
{
m_data[atom].second[coord]+=delta;

}
} ;


void find_modes(double * array, const IdxTy n)
{
/*
void dgeev_(char * JOBVL, char * JOBVR, int * N, double * A,
int * LDA, double * WR, double * WL, double * VL,int * LDVL,
double * VR, int * LDVR, double * WORK, int * LWORK, int * INFO); 
*/
char  JOBVL='N';
char  JOBVR='V';
int N=n; // m_port.size();
double * A =array; // m_port.raw();
int LDA=N;
double WR[N],WI[N];
double * VL=0;
// complex eigenvalues return a real,imag pair of vectors LOL
double * VR= new double[N*N];
int  LDVL=1;
int  LDVR=N;
int LWORK=N*N;
double WORK[LWORK];
int INFO=0;
// these are actually omega SQUARED
dgeev_(&JOBVL,&JOBVR,&N, A, &LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK, &LWORK
, &INFO);
info()<<MM_MARK<<" info= "<<INFO<<CRLF;
// Hartrees/au^2/am -> mks finally s-1 to ev 
// amu-> grms/mole / 6e23 -> grams/atom * 1e-3 -> kg/atom
const D uc=6.023e23/1e-3/.529e-10/.529e-10*4.36e-18;
const D hbar_ev=4.136e-15;
for (int i=0; i<N; ++i)
{
const D ev=::sqrt(::fabs(WR[i]*uc))*hbar_ev;
const D cminverse=ev*8065.5;
 info()<<MM_MARK<<" "<<i<<" "<<ev<<" ( "<<cminverse<<" ) "<<WR[i]<<" "<<WI[i]<<CRLF;
}

//delete [] A;
delete [] VR;
}





};

#endif


