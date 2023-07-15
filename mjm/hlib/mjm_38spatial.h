#ifndef MJM_38SPATIAL_H__
#define MJM_38SPATIAL_H__
/***

A collection of methods for dealing with points in space such as
atoms in molecule.

*/

#include "mjm_globals.h"
#include "mjm_matrix_ops.h"
//#include "mjm_strings.h"
//#include "mjm_io.h"
#include "mjm_graphs.h"
#include "mjm_templates.h"
//#include "mjm_timing.h"
//#include "mjm_config.h"
// sort 
#include <algorithm>
#include <vector>
#include <map>
#include <sstream>

//typedef mjm_generic_traits mjm_node_traits;


class mjm_default_38spatial_traits : public mjm_generic_traits
{
private:
typedef mjm_default_38spatial_traits  Myt;

typedef mjm_generic_traits Tr;
public:

typedef unsigned int IdxTy;
typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;

typedef Tr::StrTy FieldTy;
typedef Tr::StrTy StrTy;
typedef std::vector<FieldTy> LineTy;

typedef IdxTy NicTy; // nic name type should be efficient
typedef NicTy NodeNic;
typedef NicTy EdgeNic;

typedef std::vector<NodeNic> NeighTy;
//typedef std::map<NodeNic,NeighTy> AdjacencyTy;
typedef std::vector<EdgeNic> EdgeList;
// directed, from key to value[i] without regard to edge 
typedef std::map<NodeNic, NeighTy> AdjacencyTy;
typedef std::map<NodeNic, EdgeList> ConnectTy;
typedef std::map<FieldTy, NodeNic> NodeIdx;
typedef NodeIdx::iterator NodeIdxItor;
typedef std::map<FieldTy, EdgeNic> EdgeIdx;

typedef double CoordTy;
typedef std::vector<CoordTy> LocationTy;


//typedef std::vector<NodeTy> NodeList;
//typedef NodeList::iterator NodeListItor;
//typedef NodeList::const_iterator NodeListConstItor;


//typedef std::vector<EdgeTy> EdgeList;


}; //traits


namespace mjm_38spatial
{

typedef mjm_default_38spatial_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::FieldTy FieldTy;
typedef Tr::LineTy LineTy;
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;
typedef Tr::LocationTy LocationTy;
// in traits? 
//typedef unsigned int FlagTy;
//typedef unsigned int IdxTy;
typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;

typedef double DefaultFloat;
// line

template <class Tv > void ZeroLocation( Tv & x)
{
x.clear();
x.push_back(0);
x.push_back(0);
x.push_back(0);



}
template < class Td, class Ts, class Tstr  > void subset_equals( Td & d, const Ts & s, const Tstr & x)
{
//const bool check_x=(x.length()>0);
//const bool check_y=(y.length()>0);
const IdxTy size=s.size();
for ( IdxTy i =0; i<size; ++i)
{
const Tstr & id= Tstr(s[i]); 
if (id==x) d.push_back(i);
}


}
// humour me, "-" for location type may be order sensitive ... 
template < class Ts > DefaultFloat *  exhaustive_distance_matrix( const Ts & s)
{
// this should come from the parameters. 
typedef LocationTy Tl;
typedef DefaultFloat D;
const IdxTy size=s.size();
D  * p= new D[size*size];
D * q=p;
// should use stl conventions doooo
for ( IdxTy i =0; i<size; ++i)
// nodes only need to be converted into a location type 
// hopefully this eliminated the copy ctor call... 
{const Tl & xi= Tl(s[i]); 

std::cout<<MM_MARK<<" size is  "<<xi.size()<<CRLF;
const bool ok= ( xi.size()==3);
{
for (IdxTy j=0; j<size; ++j)
{
const Tl & yi= Tl(s[j]); 
// this is dumb ,needs a selective iteor or other thing
if (ok&&( yi.size()==3))
{
const D d=mjm::distance2(xi,yi);
std::cout<<MM_MARK<<" d is "<<d<<CRLF;
*(q++)=d;
} else *(q++)=0;
}}}
return p;
}


// anayze collection "s" of spatial things and return result in to"d"
// which supports a few basic operations

// want things like average and principal axes etc.
// then can translate etc
// bounds, center of mass, moments. 
class analyze_param
{
typedef DefaultFloat D;
typedef StrTy TestTy;
typedef std::vector<TestTy> TestList;
typedef std::vector<StrTy> StackTy;
typedef std::vector<D> VectorResults;
typedef typename std::map<StrTy, D> ResultStack;
typedef typename std::map<StrTy, VectorResults> VectorResultStack;
typedef typename VectorResultStack::iterator VItor;


TestList m_tests;
StrTy m_current;

public:
template <class Ts, class Tn> analyze_param(Ts & strings, Tn & numbers)
{




}

const TestTy   test(const IdxTy i )  const
{
// this can't return a ref to a default ?
if ( i>=m_tests.size()) return TestTy();
return m_tests[i];

}

void add_commentary( const StrTy & s) {m_comments.push_back(s); }
void current(const StrTy &nm ) { m_current=nm;
if ( m_vectors.find(nm)==m_vectors.end()) ZeroLocation(m_vectors[nm]);

 } 
// this will return a zero length one on first hit... 
operator VectorResults&() { return m_vectors[m_current]; } 
void store(const StrTy & nm ,const VectorResults & r) { m_vectors[nm]=r; } 
void dump(OsTy & os)
{

DumpMap(os);


}
private:

StackTy m_comments;
ResultStack m_scalars;
VectorResultStack m_vectors;
void DumpMap(OsTy & os)
{
const StrTy lbl="results";
VItor itor= m_vectors.begin();
while ( itor!=m_vectors.end())
{
StrTy nm=(*itor).first;
VectorResults v=(*itor).second;
os<<lbl<<" "<<nm;
for (IdxTy i=0; i<v.size(); ++i) os<<" "<<v[i];
os<<CRLF;

++itor;
}


}


}; 

template < class Tm, class Ti,class Tf  > void  do_something( Tm & results, Ti & itor, Tf & f )
{

OsTy & os =std::cout;
typedef LocationTy LoTy;
typedef LoTy::value_type N;
while ( itor)
{
//os<<MM_MARK<<" deref  matrix "<<""<<CRLF;
const LoTy & lo = LoTy( *itor);
//os<<MM_MARK<<" loading functor "<<lo.size()<<CRLF;
// this really needs to just avoid "comment" atoms.
if ( lo.size()==3) f(lo);
//const N n1=lo[0]; const N n2=lo[1]; const N n3=lo[2];
++itor;
}
os<<MM_MARK<<" solving f  "<<""<<CRLF;
f.solve(results);

}
// make a matrix of sum of moments up to dim n. 
// results MATRIX only need be double subscrippted
template < class Tm, class Ti  > void  moment_matrix( Tm & results, Ti & itor )
{
// normally these all take iterators, not objects. 
// I could doe the stl itor scheme or my own verswion using bool for end
// for now assume it is of my scheme type.
// must support *,++ and bool. 
//const IdxTy sz=s.size();
OsTy & os =std::cout;
typedef LocationTy LoTy;
typedef LoTy::value_type N;
typedef mjm_functors::moments Mo;
LoTy means; 
ZeroLocation(means);
typedef mjm_functors::average Av;
Av av;
Ti itor2=itor; // must be able to copy 
do_something(means,itor2,av);
results.current("means");
(LoTy&)results=means;
os<<MM_MARK<<" means are "<<means[0]<<" "<<means[1]<<" "<<means[2]<<""<<CRLF;
Mo mof(means[0],means[1],means[2]);
os<<MM_MARK<<" loading matrix"<<CRLF;
while ( itor)
{
//os<<MM_MARK<<" deref  matrix "<<""<<CRLF;
const LoTy & lo = LoTy( *itor);
//os<<MM_MAkkkkkkkkkkkkkkkkRK<<" loading matrix "<<lo[0]<<CRLF;
// this really needs to just avoid "comment" atoms.
if ( lo.size()==3) mof(lo);
//const N n1=lo[0]; const N n2=lo[1]; const N n3=lo[2];

++itor;
} // itor;

os<<MM_MARK<<" solvbinb matrix "<<""<<CRLF;
results.current("svd");
mof.solve((LoTy&)results);
os<<MM_MARK<<" solvbinb matrix "<<""<<CRLF;


}
// bounds, average, moments, name counts. 
// should also find closest approach from unit cell neibhbors

template < class Tm, class Tp  > void  test_summaries( Tp & results, Tm & s)
{
// this needs a mat class, blas? 


}

template < class Tm, class Tp  > void  misc_stats( Tp & results, Tm & s)
{



}

template < class Ta, class Ts > void  analyze_1( Ta & d, const Ts & s)
{

typedef LocationTy Tl;
typedef DefaultFloat D;
const IdxTy size=s.size();
typedef std::vector<IdxTy> QualTy;
QualTy quals;
// first see where they are. 
D * distances=exhaustive_distance_matrix(s);
// find histograms of C-C distances, for long molecules the max
// could be large and make bins too coarse, need to just get
// the low end of the curve. 
subset_equals(quals,s,"C");
// quals should probably just be an iterator rather than a collection.
// want first peak in histo but need to pick a bin size fine enough to fine it.





}


// this uses idiosyncratic iterators....
template < class Ts > void test_exhaustive_distance_matrix( const Ts & s)
{
typedef typename  Ts::const_iterator It;
typedef typename  Ts::const_iterator::target_name Tt;
// this should come from the parameters. 
typedef LocationTy Tl;
const IdxTy size=s.size();
// should use stl conventions doooo
It itor= s.const_itor(); 
while (itor)
{
const Tl & x= Tl(*itor); 

++itor;
} // while



}





};



#endif

