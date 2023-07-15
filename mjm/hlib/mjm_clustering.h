#ifndef MJM_CLUSTERING_H__
#define MJM_CLUSTERING_H__


/***

A collection of methods for dealing with points in space such as
atoms in molecule.

*/

#include "mjm_globals.h"
#include "mjm_conflicts.h" //  #include <cblas.h>
#include <vector>

// lapack sigs
extern "C"
{       
/*
void zheevr_(char* JOBZ, char* RANGE, char* UPLO, int * N, CoTy* A, int * LDA,
                double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
                double* W, CoTy* Z, int* LDZ, int* ISUPPZ, CoTy* WORK, int* LWORK,
                double* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO);
*/

//void cgesv_(int * N, int * NRHS, CoTy * A, int* LDA, int * IPIV, CoTy * B, int * LDB, int* INFO);

// DGEBRD real double SVD
void dgesvd_(char * JOBU, char * JOBVT,int * M, int * N, double * A,
 int *LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK, 
int * INFO);


}; // lapacks


//typedef mjm_generic_traits mjm_node_traits;
namespace mjm_clustering
{
typedef double DefaultD;
typedef unsigned int IdxTy;
class density1d 
{

typedef density1d Myt;
typedef DefaultD D;
class mypair
{public: D x,y; mypair(const D & xx, const D & yy, ): x(xx),y(yy) {}
}
typedef mypair S;
typedef std::vector<S> m_strength;
enum { DIM=1};
// just dump in other crap too LOL.
D m[DIM];
D m_min[DIM],m_max[DIM];
IdxTy m_n;
 
public:
density1d() { 
//for (unsigned int i=0; i<DIM; ++i) m[i]=0; 
//for (unsigned int i=0; i<DIM; ++i) m_min[i]=1e100; 
//for (unsigned int i=0; i<DIM; ++i) m_max[i]=-1e100; 
}

~density1d() {  } 

template < class Ts> Myt & operator()(const Ts & loc, const Ts & val)
{
m_strength.push_back(S(loc,val));
return *this;
}

D f( const D & sample, const D & grid)
{ return 1e-1+::fabs(sample-grid); } 
template < class Td > void solve(Td & result)
{
// should not throw on zero n
// compute density at each point on grid, differentiated peak detect
const IdxTy pts= m_strength.size();
D maxx=0;
// these are values
for ( IdxTy i=0; i<pts; ++i)  if (m_strength[i].x>maxx ) maxx=m_strength[i].x;


const D step=range/gpts;
for ( IdxTy j=0; j<gpts; ++j) 
{
// uniform could just add but expect to make non-uniform
const D pt=j*step;
for ( IdxTy i=0; i<pts; ++i) 
{
D rho=0;

for ( IdxTy i=0; i<pts; ++i) { const S&  v=m_strength[i]; rho+=v.y/f(v.x,pt); }



}
}  

}
private:

}; // average














}; // clusterers  ns
#endif

