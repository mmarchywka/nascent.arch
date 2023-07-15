#ifndef MJM_GENERIC_POTENTIALS_H__
#define MJM_GENERIC_POTENTIALS_H__

#include <math.h>
#include <map>
#include <string>
#include <sstream>


namespace generic_potentials
{

typedef double D;
typedef int I ;
typedef std::string StrTy ;
typedef unsigned int IdxTy ;
typedef std::map<StrTy, D> LutTy;

template<class Ty> D dist2(const Ty & v, const I sz)
{ D x=0; for ( I i=0; i<sz; ++i)
{
//MM_MSG(" in loop "<<i)
 x+=v[i]*v[i];
}

 return x; }
template<class Ty> D dist(const Ty & v, const I sz) { return ::sqrt(dist2(v,sz)); } 
template<class Ty> D r3(const Ty & v, const I sz) { D x= (dist(v,sz)); return x*x*x; } 
template<class Ty> D r5(const Ty & v, const I sz) { D x= (dist(v,sz)); return x*x*x*x*x; } 

// this is the derivative of dist or r^2 
template<class Ty> D dxr2(const Ty & v, const I sz,const I d) { return 2.0*v[d]; }
// second der wrt r^2
template<class Ty> D d2xr2(const Ty & v, const I sz,const I d1,const I d2) { return (d1==d2)?2.0:0; }
// r= ( x*x+y*y)^(1/2)
// dr/dx=x/r
// derivative of 1/r in direction d -> -2(1/2)x/r^2 
template<class Ty> D dxr(const Ty & v, const I sz,const I d) { return v[d]/dist(v,sz); }

// d( 1/r)/dx =  -x/r^(3)
template<class Ty> D dx_over_r(const Ty & v, const I sz,const I d) { return -v[d]/r3(v,sz); }

// d^2( 1/r)/(dxdy) = //   -x/r^(3)
template<class Ty> D d2x_over_r(const Ty & v, const I sz,const I d1,const I d2) 
{ 
const bool diag=(d1==d2);
//D ddr=dist(v,sz);
D od=-3*v[d1]*v[d2]/r5(v,sz);
//return -d2x(v,sz,d1,d2)/r3(v,sz) +2.0*dx(v,sz,d1)*dx(v,sz,d2)/ddr;
if ( diag) { return -1.0/r3(v,sz)+ od; }else { return od; }
}

LutTy & charge_map()
{
static LutTy m;
static bool inited=false;
if ( !inited)
{
m["c"]=4;
m["C"]=4;
m["H"]=1;
m["h"]=1;


inited=true;
}
return m; 
}

template <class Td,class Tl, class Tm  > void charges( Td *  is,  Tl & lbl,  Tm & mmap)
{
for (IdxTy i=0; i<lbl.size(); ++i) is[i]=mmap[StrTy(lbl[i][0])];  
for (IdxTy i=0; i<lbl.size(); ++i) MM_MSG(" charge "<<is[i]<<" "<<mmap[StrTy(lbl[i][0])])  

}

// load form a input stream using something like jdftx ionpos format 
template <class Ty,class Tv,class Tl > void load( Ty * is, Tv & v, Tl & lbl)
{
typedef std::stringstream Ss;
typedef typename Tl::value_type LblTy;
typedef typename Tl::value_type::value_type StrTy;
typedef typename Tv::value_type LocTy;
typedef typename Tv::value_type::value_type PosTy;
StrTy nm;
PosTy xx;
while ( !is->eof()&&is->good())
{
LocTy lo;
LblTy lb;
const int LEN=1<<10;
char line[LEN];
(*is).getline(line,LEN-1);
if ( !(!is->eof()&&is->good())) break;
Ss ss(line);
bool add_this_one=true;
if ( add_this_one)
{
StrTy ion,species;
(ss)>>ion;
if ( ion!="ion")
{
std::cout<<MM_MARK<<" ignoring "<<ion<<" which should be ion..."<<CRLF;
continue;
}
//(*is)>>species>>x>>y>>z>>flag;
ss>>nm; lb.push_back(nm);
lbl.push_back(lb);
ss>>xx; lo.push_back(xx);
ss>>xx; lo.push_back(xx);
ss>>xx; lo.push_back(xx);
v.push_back(lo);


} // add
}
}



/* fudd 
template <class Ty,int sz> Ty distance( const Ty & v, const int sz)
{
return v[sz-1]*v[sz-1]+distance(v,sz-1); 
}

template <> Ty distance<Ty,0>( const Ty & v, const int sz ) { return 0; } 
*/


}; // namespace



#endif //MJM_GENERIC_POTENTIALS_H__
