#ifndef MJM_CONFORMALS_H__
#define MJM_CONFORMALS_H__

#include "mjm_globals.h"
#include "mjm_generic_iterators.h"
#include "mjm_rational.h"
#include "mjm_block_matrix.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <math.h>
#include <string.h>


/*
conformal mapping tests 
vi mjm_block_matrix.h 
// #ifdef TEST_CONFORMALS_MAIN__
g++ -DTEST_CONFORMALS_MAIN__  -std=gnu++11 -Wall -I.. -x c++ mjm_conformals.h 
 2331  ./a.out



https://arxiv.org/pdf/physics/0604033


*/
class mjm_conformals 
{
typedef mjm_conformals Myt;
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
typedef mjm_generic_itor<2> LoopItor2;
typedef LocTy location_type;
typedef mjm_rational RatTy;


static D  rat_power_map( const RatTy & x, const RatTy & y, const RatTy & theta, const RatTy & p, const RatTy & xz,const RatTy & yz) 
{

D sol=0;
RatTy xe=x; xe.subtract(xz);
RatTy ye=y; ye.subtract(yz);
 D theta1=std::atan2(ye.approx(),xe.approx());
if (theta1<0) theta1=2*M_PI+theta1; // 0 to 2pi 
// if this is inside the wedge the result is zero 
const D thetaeff=(theta1-theta.approx());

if (thetaeff>M_PI/p.approx()) return sol;
// this may not work well lol if this and that are the same. 
const D r=xe.mult(xe).approx()+ye.mult(ye).approx();
const D s=std::sin(p.approx()*(thetaeff));
const D result=s*std::pow(r,.5*p.approx());

sol=(result);

return sol; 

}



} ; // mjm_block_matrix

#ifdef TEST_CONFORMALS_MAIN__






int main(int argc, char ** argv)
{
typedef mjm_conformals T;
typedef mjm_generic_itor<2> Itor;
typedef mjm_rational RatTy;
typedef unsigned int IdxTy;
//MM_ERR("location")
RatTy x,y,theta,p,xz,yz;
RatTy theta2,xz2,yz2;
Itor itor(100,100);
x.set(-50);
y.set(-50);
p.set(2.0/3.0);
theta2.set(M_PI);
xz2.set(25);
yz2.set(25);
if ( itor.done()) return 0 ;
IdxTy bump=1;
while (true)
{
switch (bump)
{ 
case 2: return 0;
case 1: x.set(int(itor[0])-50);
case 0: y.set(int(itor[1])-50);
double d=T::rat_power_map( x, y, theta,  p,  xz, yz) ;
double f=T::rat_power_map( x, y, theta2,  p,  xz2, yz2) ;
MM_MSG(" c0nform "<<x.approx()<<" "<<y.approx()<<" "<<d<<" "<<f<<" "<<(d+f))
} // swith 
bump=itor.inc();

} // true
//MM_ERR("ctor")

return 0;
}


#endif // main



#endif

