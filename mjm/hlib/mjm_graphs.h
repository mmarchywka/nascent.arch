#ifndef MJM_GRAPHS_H__
#define MJM_GRAPHS_H__
/***
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
#include "mjm_strings.h"
#include "mjm_io.h"
#include "mjm_timing.h"
#include "mjm_config.h"
// sort 
#include <algorithm>
#include <vector>
#include <map>

//typedef mjm_generic_traits mjm_node_traits;


class mjm_node_traits : public mjm_generic_traits
{
private:
typedef mjm_node_traits  Myt;

typedef mjm_generic_traits Tr;
public:
/*
typedef Tr::ChTy ChTy; 
typedef Tr::ErTy ErTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
*/

typedef unsigned int IdxTy;
typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;

typedef Tr::StrTy FieldTy;
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


//typedef std::vector<NodeTy> NodeList;
//typedef NodeList::iterator NodeListItor;
//typedef NodeList::const_iterator NodeListConstItor;


//typedef std::vector<EdgeTy> EdgeList;


}; //traits



class mjm_node
{
private:
typedef mjm_node  Myt;

typedef mjm_node_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::FieldTy FieldTy;
typedef Tr::LineTy LineTy;
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;

// in traits? 
//typedef unsigned int FlagTy;
//typedef unsigned int IdxTy;
typedef mjm_io IoTy;
typedef IoTy::OsTy OsTy;
typedef IoTy::IsTy IsTy;
// line

LineTy m_node;
// if this is not zero, then nodes can exist without a valid name 
IdxTy m_name;
NodeNic m_nic;
public:
mjm_node():m_name(0),m_nic(0) {append(""); }
mjm_node(const FieldTy & name,  const NodeNic & nic )
: m_name(0), m_nic(nic)
{ m_node.push_back(name); } 

mjm_node(const LineTy & that, const IdxTy n, const NodeNic & nic )
	:m_name(n),m_nic(nic) {m_node=that; }
void append(const StrTy & v) { m_node.push_back(v); } 
Myt & set(const LineTy & node) { m_node=node; return *this; }
Myt & set(const NodeNic & nic ) { m_nic=nic; return *this; } 
// for non-zed names, should check size or segfault1
const StrTy & name() const { return m_node[m_name]; }
const StrTy & value(const IdxTy & i) const { return m_node[i]; }
IdxTy size() const { return m_node.size(); }
const NodeNic & nic() const { return m_nic; }

}; //node



class mjm_edge
{
private:
typedef mjm_edge  Myt;

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
mjm_edge():m_name(0),m_nic(0),m_from(0),m_to(0)  {append(""); }
mjm_edge(const LineTy & that, const IdxTy n, const NicTy & nic, const NicTy & from, const NicTy & to  )
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


}; //edge/




class mjm_graph
{
private:
typedef mjm_graph  Myt;

typedef mjm_node_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
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


public:
mjm_graph() : m_name(0) {append(""); }
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
// not equal? No idea what unsigned gt is here.. 
if (es>0U) os<<" "<<(*p)(e.value(0).c_str(),0);
for (IdxTy j=1; j<es; ++j) os<<" "<<(*p)(e.value(j).c_str(),0);
}
os<<CRLF;

os.flush();

}

// split into disconnected regions
void make_subgraphs(const IdxTy & flags )
{
enum { VERBOSE_FLAG=1,WARN_CYCLE=2};

const bool verbose=(flags&VERBOSE_FLAG);
const bool warn_cycle=(flags&WARN_CYCLE);

// a vector of subgraphs
SgList & sgl=m_sgl;
// maps node nic to subgraph
SgMap & sgm=m_sgm;
// should match length of sgl, 1= used 0= obsolste
SgStatus & sgs=m_sgs;


const IdxTy sz=m_edges.size();
for (IdxTy i=0; i<sz; ++i)
{
const EdgeTy & e= m_edges[i];
const NodeNic & f=e.from();
const NodeNic & t=e.to();
SgMapItor smf=sgm.find(f);
SgMapItor smt=sgm.find(t);
if (verbose) { osx()<<MM_MARK<<" edge "<<i<<" connects "<<f<<" "<<m_nodes[f].name()<<" to "<<t<<" "<<m_nodes[t].name()<<CRLF; osx().flush(); }
// new node
const bool new_from= (smf==sgm.end()) ;
const bool new_to= (smt==sgm.end()) ;
// create new subgraph and add edge
if ( new_from&&new_to)
{
const IdxTy sg_nic=sgl.size();
 SubGraph  sg= SubGraph();
// this messes up the node nics, we don't need anyway
// sg.add(e);
// add both nodes to table
sgm[f]=sg_nic;
sgm[t]=sg_nic;
if (verbose) { osx()<<MM_MARK<<" edge "<<i<<" both new sq is  "<<sg_nic<<CRLF; osx().flush(); }
// indicate this is here and in use
sgl.push_back(sg);
sgs.push_back(1);
continue; 

}
// add the from node to the to graph list
if ( new_from&&!new_to) { sgm[f]=sgm[t]; continue;  }
if ( !new_from&&new_to) { sgm[t]=sgm[f]; continue;  }


if ( !new_from &&! new_to)
{
// we have the itors above, faster to dref
const IdxTy & from_sg=sgm[f];
const IdxTy & to_sg=sgm[t];
if (verbose) { osx()<<MM_MARK<<" edge "<<i<<" both old from  "<<from_sg<<" to "<<to_sg<<CRLF; osx().flush(); }
//SubGraph & fg= sgl[from_sg];
// SubGraph & tg= sgl[to_sg];

// want to move smaller but this ok for now,
// presume that the smaller idx is larger
// just merge the node associations, not the empty graphs
//if ( from_sg<to_sg) { fg.add(tg); } else { tg.add(fg); } 
// update the node tables
// this is brure force, merges should not be real common? LOL
SgMapItor ii=sgm.begin();
SgMapItor ie=sgm.end();
if ( from_sg==to_sg)
{
if (warn_cycle) { osx()<<MM_MARK<<" edge "<<i<<" cycle  "<<m_nodes[f].name() <<" to "<<m_nodes[t].name()<<CRLF; osx().flush(); }
 continue; // ????? cycle? 
}
if ( from_sg<to_sg)
 {
sgs[to_sg]=0; // mark the to subgraph as unused, 
// repoint the nodes in "to_sg" to from_sg
// this is pretty stpid, but worry about it lateer
while (ii!=ie) { if ((*ii).second==to_sg) (*ii).second=from_sg; ++ii;}  
  } 
else 
{  

sgs[from_sg]=0; // mark the to subgraph as unused, 
// repoint the nodes in "to_sg" to from_sg
// this is pretty stpid, but worry about it lateer
while (ii!=ie) { if ((*ii).second==from_sg) (*ii).second=to_sg; ++ii;}  
} 
continue;
} //if else  merge


} // loop 



} //subgraph 



private:
// this messes up merges, where the nic refers
// to the graph index. 
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







};














/*
class mjm_csv
{
private:
typedef mjm_netwerk  Myt;

typedef mjm_generic_traits Tr;



*/





#endif

