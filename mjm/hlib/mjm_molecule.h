#ifndef MJM_MOLECULE_H__
#define MJM_MOLECULE_H__
/***

This needs to make a collection of "atoms" into a sort of graph
with location awareness. Special node that has location and
fuzzy bonds or connections. 

mjm_molecule.h : based on graph classes
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

#include "mjm_globals.h"
//#include "mjm_strings.h"
//#include "mjm_io.h"
#include "mjm_graphs.h"
#include "mjm_templates.h"
#include "mjm_38spatial.h"
//#include "mjm_timing.h"
//#include "mjm_config.h"
// sort 
#include <algorithm>
#include <vector>
#include <map>
#include <sstream>

//typedef mjm_generic_traits mjm_node_traits;


class mjm_molecule_traits : public mjm_generic_traits
{
private:
typedef mjm_molecule_traits  Myt;

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



class mjm_atom : public mjm_node
{
private:
typedef mjm_atom  Myt;

typedef mjm_molecule_traits Tr;

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
// line

//LineTy m_node;
// if this is not zero, then nodes can exist without a valid name 
//IdxTy m_name;
//NodeNic m_nic;

StrTy m_species;
LocationTy m_location;
enum { DEFAULT=~0U, REAL_ATOM=0, COMMENT=1};

IdxTy m_class; // real, comment, etc
// full line used to make this item, mostly for comments
StrTy m_full_line;


public:
mjm_atom():m_species(""), m_location(), m_class(DEFAULT),m_full_line() {

Clear();

}
void set(const LocationTy & there, const StrTy & species,const StrTy & line)
{
m_location=there;
m_species=species;
if ( m_species.length()!=0) m_class=REAL_ATOM;
m_full_line=line;

} 
// constness not an issue now?  The stupid [] is not const... 
const LocationTy & location() const { return m_location; } 
LocationTy & location()  { return m_location; } 
const StrTy &  species() const { return m_species; } 
const StrTy & full() const { return m_full_line; }
bool is_atom() const { return m_class==REAL_ATOM; } 
operator LocationTy() { return m_location; }
operator const LocationTy() const { return m_location; }
operator const StrTy() const { return m_species; }
 
/*
mjm_atom():m_name(0),m_nic(0) {append(""); }
mjm_atom(const FieldTy & name,  const NodeNic & nic )
: m_name(0), m_nic(nic)
{ m_node.push_back(name); } 

mjm_atom(const LineTy & that, const IdxTy n, const NodeNic & nic )
	:m_name(n),m_nic(nic) {m_node=that; }
void append(const StrTy & v) { m_node.push_back(v); } 
Myt & set(const LineTy & node) { m_node=node; return *this; }
Myt & set(const NodeNic & nic ) { m_nic=nic; return *this; } 
// for non-zed names, should check size or segfault1
const StrTy & name() const { return m_node[m_name]; }
const StrTy & value(const IdxTy & i) const { return m_node[i]; }
IdxTy size() const { return m_node.size(); }
const NodeNic & nic() const { return m_nic; }
*/
private:
void Clear()
{
m_location.clear();
m_location.push_back(0);
m_location.push_back(0);
m_location.push_back(0);

}


}; //node



class mjm_relation : public mjm_edge
{
private:
typedef mjm_relation  Myt;

typedef mjm_node_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::LineTy LineTy;
typedef Tr::NicTy NicTy;

typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;
// line

LineTy m_node;
// if this is not zero, then nodes can exist without a valid name 
IdxTy m_name;
NicTy m_nic,m_from,m_to;
public:
/*
mjm_relation():m_name(0),m_nic(0),m_from(0),m_to(0)  {append(""); }
mjm_relation(const LineTy & that, const IdxTy n, const NicTy & nic, const NicTy & from, const NicTy & to  )
	:m_name(n),m_nic(nic),m_from(from), m_to(to) {m_node=that; }
void append(const StrTy & v) { m_node.push_back(v); } 
Myt & set(const LineTy & node) { m_node=node; return *this; }
Myt & set(const NicTy & nic ) { m_nic=nic; return *this; } 
Myt & set(const NicTy & from, const NicTy & to ) { m_from=from; m_to=to;  return *this; } 
// for non-zed names, should check size or segfault1
const StrTy & name() const { return m_node[m_name]; }
const StrTy & value(const IdxTy & i) const { return m_node[i]; }
IdxTy size() const { return m_node.size(); }
const NicTy & nic() const { return m_nic; }
const NicTy & from() const { return m_from; }
const NicTy & to() const { return m_to; }
*/

}; //edge/



// this does not derive from graph
class mjm_molecule
{
private:
typedef mjm_molecule  Myt;

typedef mjm_molecule_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::StrTy StrTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;

typedef mjm_atom AtomTy;
/*
typedef Tr::LineTy LineTy;
typedef Tr::FieldTy FieldTy;
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;
typedef Tr::NodeIdx NodeIdx;
typedef Tr::NodeIdxItor NodeIdxItor;
typedef Tr::EdgeIdx EdgeIdx;
//typedef Tr::NodeList NodeList;
//typedef Tr::EdgeList EdgeList;


typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;

typedef mjm_io IoTy;
typedef mjm_node NodeTy;
typedef mjm_edge EdgeTy;

typedef std::vector<NodeTy> NodeList;
typedef NodeList::iterator NodeListItor;
typedef NodeList::const_iterator NodeListConstItor;

typedef std::vector<EdgeTy> EdgeList;



// line

NodeList m_nodes;
NodeIdx m_node_idx;
EdgeList m_edges;
// I guess the edges can have repeated names? 
// EdgeIdx m_edge_idx;

LineTy m_graph;
// if this is not zero, then nodes can exist without a valid name 
IdxTy m_name;

// for analysis
typedef Myt SubGraph;
typedef NicTy SgNic;
// and idx would map names to nic of self this maps to graph
typedef std::map<NodeNic, SgNic> SgMap;
typedef SgMap::iterator SgMapItor;
typedef std::vector<SubGraph> SgList;
// rather than move things, just bitmap the obsolste
typedef std::vector<IdxTy> SgStatus;
// a vector of subgraphs
SgList m_sgl;
// maps node nic to subgraph
// node nic was assigned when it was added to the main graph
// and that points to a sub graph in m_sgl
SgMap m_sgm;
// should match length of sgl, 1= used 0= obsolste
SgStatus m_sgs;
*/
typedef std::vector<AtomTy> VectorTy;
VectorTy m_atoms;
class myitor
{

Myt & m;
IdxTy m_ptr,m_size;
public: 
myitor(Myt & mol) : m(mol),m_ptr(0),m_size(m.size()) { }
AtomTy & operator*() { return m.atom(m_ptr); }  
typedef AtomTy target_name;
myitor & operator++() { ++m_ptr; return *this; } 
 operator bool() { return m_ptr<m_size; } 
}; 

class myconst_itor
{

const Myt & m;
IdxTy m_ptr,m_size;
public: 
myconst_itor(const Myt & mol) : m(mol),m_ptr(0),m_size(m.size()) { }
const AtomTy & operator*() { return m.atom(m_ptr); }  
myconst_itor & operator++() { ++m_ptr; return *this; } 
typedef const AtomTy target_name;
 operator bool() { return m_ptr<m_size; } 
}; 




public:

typedef AtomTy atom_type;

typedef myitor iterator;
typedef myconst_itor const_iterator;
void add(const AtomTy & a) { m_atoms.push_back(a); } 
IdxTy atoms() const { return m_atoms.size(); } 
IdxTy size() const { return atoms(); } 

const AtomTy & atom(const IdxTy i) const { return m_atoms[i]; } 
AtomTy & atom(const IdxTy i) { return m_atoms[i]; } 
AtomTy & operator[](const IdxTy i) { return m_atoms[i]; } 
const AtomTy & operator[](const IdxTy i) const { return m_atoms[i]; } 

iterator itor()   {return iterator(*this); }
const_iterator const_itor() const  {return const_iterator(*this); }
//iterator end() const  {}

/*mjm_molecule() : m_name(0) {append(""); }
void append( const StrTy & v) {m_graph.push_back(v); }

// add a "line" in the form that it represents
// an edge with 2 nodes defined only by name.
// this is fine when you don't have const nics, 
void add(const LineTy & l, const IdxTy & edge, const IdxTy & from, const IdxTy & to)
{
const FieldTy & fr_node=l[from];
const FieldTy & to_node=l[to];
const NicTy & fr_nic=addNode(fr_node);
const NicTy & to_nic=addNode(to_node);
const NicTy & edge_nic=m_edges.size();

m_edges.push_back(EdgeTy(l,edge,edge_nic,fr_nic,to_nic));

} // add
// line is complete node spec, n is index for field containing name
void add_node(const LineTy & l, const IdxTy & n, const bool flag_dups)
{
	addNode(l,n,flag_dups); 
}
void add(const EdgeTy & e)
{ m_edges.push_back(e); }
void add(const Myt & g)
{
const EdgeList & el=g.m_edges;
const IdxTy n=el.size();
for ( IdxTy i=0; i<n; ++i)
{ add(el[i]); }

}
NodeList & node_list() { return m_nodes; } 
*/

/*
// this needs to take a template param pointer
// for a escaping function 
void dump_subgraphs(OsTy & os)
{
// const IdxTy nic=m_nodes.nic();
const IdxTy n=m_nodes.size();
for (IdxTy i=0; i<n; ++i)
{
os<<m_sgm[i];
const NodeTy & no= m_nodes[i];
const IdxTy l=no.size();
for ( IdxTy j=0; j<l; ++j)
{
os<<" "<<no.value(j);
}
os<<CRLF;
}
os.flush();


}

void dump_edges_by_subgraph(OsTy & os)
{
dump_edges_by_subgraph(os,0);
}
template <mjm_strings::EscPtr p=& mjm_strings::add_nothing > void dump_edges_by_subgraph(OsTy & os, const IdxTy & flag)
{
enum { SHOW_SG=1, STUPUD=~0};
const bool show_sg=(flag&1);

const IdxTy sz=m_edges.size();
for (IdxTy i=0; i<sz; ++i)
{
const EdgeTy & e= m_edges[i];
const NodeNic from=e.from();
const NodeNic to=e.to();
const IdxTy sg= (m_sgm.find(from)!=m_sgm.end())?m_sgm[from]:STUPUD;
if ( m_sgm[to]!=sg) 
{osx()<<MM_MARK<<" ASSERT FAILS FOR SUBGRAPH AT EDGE "<<i<<" "<<to<<" from "<<from<<CRLF; osx().flush(); }
const IdxTy es=e.size();
if ( show_sg) { if ( sg==(IdxTy)STUPUD ) os<<"-1"<<" "; else os<<sg<<" ";} 
if (es>0U) os<<" "<<(*p)(e.value(0).c_str(),0);
for (IdxTy j=1; j<es; ++j) os<<" "<<(*p)(e.value(j).c_str(),0);
}
os<<CRLF;

os.flush();

}
*/



private:
// this messes up merges, where the nic refers
// to the graph index. 
/*
const NicTy & addNode(const FieldTy & n)
{
NodeIdxItor ni=m_node_idx.find(n);
if ( ni==m_node_idx.end())
{
// tempted to index from 1 so zero is invliad?
const NicTy & nic = m_nodes.size();
m_nodes.push_back( NodeTy(n,nic));  
m_node_idx[n]=nic;
return nic;
}
// this is kind of dum, we just did this a minute ago
return (*ni).second; //  m_node_idx[n];

}

// these should beloaded apriori, before edges but dups
// could just add details to node stubs
const NicTy & addNode(const LineTy & line, const IdxTy n, const bool flag_dups )
{
const FieldTy & name = line[n];
NodeIdxItor ni=m_node_idx.find(name);
if ( ni==m_node_idx.end())
{
// tempted to index from 1 so zero is invliad?
const NicTy & nic = m_nodes.size();
m_nodes.push_back( NodeTy(line,n,nic));  
m_node_idx[name]=nic;
return nic;
}
if (flag_dups) { osx()<<MM_MARK<<" duplicate node name "<<name<<CRLF; 
osx().flush(); } 


// this is kind of dum, we just did this a minute ago
return (*ni).second; //  m_node_idx[n];

}

*/






};



class mjm_molecule_io
{

private:
typedef mjm_molecule_io Myt;
typedef mjm_molecule_traits Tr;
typedef Tr::ChTy ChTy; 
typedef Tr::StrTy StrTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::LineTy LineTy;
typedef Tr::FieldTy FieldTy;
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;
typedef Tr::NodeIdx NodeIdx;
typedef Tr::NodeIdxItor NodeIdxItor;
typedef Tr::EdgeIdx EdgeIdx;
typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;
typedef Tr::LocationTy LocationTy;
typedef Tr::CoordTy  D;
// should be traits
typedef std::stringstream SsTy;
typedef std::vector<StrTy>  VecTy;

// these are part of api should be public 
typedef mjm_molecule MolTy;
typedef MolTy::atom_type AtomTy;
//typedef MolTy::location_type LocationTy;
typedef IdxTy StatusTy;
public:
// for now,just make default ctor and then load etc
mjm_molecule_io():m_debug_os(& std::cout) {}
StatusTy  load(MolTy & d, const StrTy & s) { return Load(d,s); }
StatusTy  store(const MolTy & d, OsTy  & os) { return Store(d,os); }


private:
class fileDesc
{
}; // fileDesc
typedef fileDesc AdHocFileDesc;

StatusTy  Store(const MolTy & d, OsTy  & os)
{
const StrTy sep=" ";
MolTy::const_iterator itor = d.const_itor();
while ( itor)
{
const AtomTy & a= *itor;
if ( a.is_atom())
{
const LocationTy & al=a.location();

os<<a.species()<<sep<<al[0]<<sep<<al[1]<<sep<<al[2];

} else os<<a.full();
os<<CRLF;

++itor;
}
return 0; 
}


StatusTy  Load(MolTy & d, const StrTy & s)
{

IsTy * is = FindStream(s);
const StatusTy & x= LoadAdHoc(d,is, AdHocFileDesc());
CloseStream(is);
return x; 
}
StatusTy LoadAdHoc(MolTy &d, IsTy * is, const AdHocFileDesc & ahf )
{
const IdxTy size=1<<16;
ChTy buf[size];

while ( ok(is))
{
(*is).getline(buf,size);
StrTy x = StrTy(buf);
AddLine(d,x,ahf);

}


return 0; 
}
static D af(const StrTy & s) { return ::atof(s.c_str()); }

StatusTy AddLine(MolTy & d, const StrTy & line, const AdHocFileDesc &ahf)
{
VecTy v;
Parse(v,line,ahf);
// lines not part of molecule should be kept in separate location erally....
AtomTy atom;
LocationTy l;
StrTy species="";
// if possible extract info else just insert as a comment.
if ( v.size()==4)
{
species=v[0];
l.push_back(af(v[1]));
l.push_back(af(v[2]));
l.push_back(af(v[3]));

}
atom.set(l,species,line);


d.add(atom);

return 0; 
}
void Parse(VecTy & v, const StrTy & line, const AdHocFileDesc & ahf)
{
SsTy ss(line);
while (ok(&ss)) { StrTy x ; ss>>x; v.push_back(x); } 

}


static bool ok(IsTy * is ) { return ( (*is).good()&&!(*is).eof()); } 

IsTy * FindStream(const StrTy & s)
{
IsTy * is = mjm_io::get_in(s.c_str(), *m_debug_os);

return is;

}





void CloseStream(IsTy * is)
{
delete is;
}


private:

OsTy  * m_debug_os;

} ; //molecule_io













#endif

