/***
mjm_neighbors.h : find nearby things 
mjm_graphs.h : graph classes
mjm_tree.cpp : tree util 
mjm_perf.cpp: launch a bunch of threads and dump timing stats.
derived from, used copies of headers etc:
logs_db_lut.cpp : integrate the various log and othere tasks inot one place

duphus_agg.cpp: derioved from 
x-ui-util.cpp : Phluant x-windows manipulator for qa and other
image needs. 

Mike Marchywka May 2011
May 30, 2011 : ran fine overnight but mem diag is zed,
start adding "Break" after load LOL, input escapes, temp files
on too big etc


*/

#include "mjm_graphs.h"
// sort 
#include <algorithm>
#include <vector>
#include <map>



class mjm_neighbors
{
private:
typedef mjm_neighbors  Myt;

typedef mjm_node_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::FieldTy FieldTy;
typedef Tr::LineTy LineTy;
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;
typedef Tr::NodeList NodeList;
typedef Tr::NodeListConstItor  NodeListConstItor;
typedef Tr::EdgeList EdgeList;
typedef Tr::NeighTy NeighTy;
typedef Tr::AdjacencyTy  AdjacencyTy;
typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;
// line

// we never own this so do not delete
const NodeList * m_nl;
// this maps node nics in DG 
AdjacencyTy m_ad;
StrTy m_desc;
// 
template <class tag, int param=2 > coord
{
typedef coord<tag,param> Myt;
public:
typedef double location_type;
typedef double distance_type;
enum {DIMENSION=param};
// this is NOT initialized, caveat emptor 
location_type x[DIMENSION];  

void locate(const NodeTy & n)
{// for geo data really need exact conversion of make sure that atof is ok


}
distance_type operator(const Myt & a, const Myt & b)
{
location_type dx=a.x[0]-b.x[0];
location_type dy=a.x[1]-b.x[1];
return dx*dx+dy*dy;
}


};


template <class tag, class coord_tag > class distance
{

typedef coord<coord_tag> CoordTy;

typedef std::map<NodeNic,CoordTy > NodeLocations;
NodeLocations m_nl;
CoordTy m_d;
public:
typedef double distance_type;
void set(const NodeList & nl)
{
NodeListConstItor ii=nl.begin();
NodeListConstItor ie=nl.end();
CoordTy d=m_d;
while (ii!=ie)
{
const Node & n = *ii;
d.locate(n);
// this makes a copy I hope LOL 
m_nl[n.nic()]=d;
++ii;
} 


}
distance_type operator() (const NodeNic & a, const NodeNic & b)
{
 return m_d(m_nl[a],m_nl[b]);
}


} 
// this is an approach which will NOT work 
/*
template <class tag > distance
{
// this will implement euclid for illustration.
typedef std::vector<IdxTy> FieldList;
FieldList m_list;
IdxTy m_sz;
public:
distance():m_sz(0) {}
void add_field(const IdxTy i) { m_list.push_back(i); m_sz=m_list.size(); } 
// specializations can do whatever they want but
// users expect a distance_type and this operator
typedef double distance_type;
distance_type operator() ( const NodeTy & a, const NodeTy b) 
{
distance_type sse=;
for ( IdxTy i=0; i<m_sz; ++i)
{
// now this is where the strings are a mess... 



}
return sse;

}

};
*/

public:


mjm_neighbors():m_nl(0),m_ad(),m_desc("")  { }
mjm_node(const StrTy & desc,  const NodeList * nl )
: m_nl(nl), m_desc(desc)
{  } 

const StrTy & name() const { return m_desc; }







}; //neighbors















/*
class mjm_csv
{
private:
typedef mjm_netwerk  Myt;

typedef mjm_generic_traits Tr;



*/






