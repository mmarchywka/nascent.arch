#ifndef MJM_MOM_H__
#define MJM_MOM_H__
/***
mjm_mom.h : map of maps, can generalize but just a 2 lvel for now
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
#include "mjm_graphs.h"
//typedef mjm_generic_traits mjm_node_traits;





class mjm_mom
{
private:
typedef mjm_mom  Myt;

typedef mjm_node_traits Tr;

typedef Tr::ChTy ChTy; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::ErTy ErTy;
typedef Tr::LineTy LineTy;
typedef Tr::FieldTy FieldTy;
// could just use generic traists? 
/*
typedef Tr::NicTy NicTy;
typedef Tr::NodeNic NodeNic;
typedef Tr::NodeIdx NodeIdx;
typedef Tr::NodeIdxItor NodeIdxItor;
typedef Tr::EdgeIdx EdgeIdx;
//typedef Tr::NodeList NodeList;
//typedef Tr::EdgeList EdgeList;
*/

typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;

typedef mjm_io IoTy;
/*
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
public:   class branch_class
{
private:
typedef branch_class Myt;
typedef FieldTy KeyTy;
typedef LineTy TerminalTy;
typedef std::map<KeyTy,Myt> LowerTy;
typedef LowerTy UpperTy;
// shouldhae cont here. 
typedef LowerTy::iterator LoItor;
// these have to point to instance owned elsewhere. 
TerminalTy const *  m_terminal; 
// these are owned by us
LowerTy * m_lower;
UpperTy * m_mom;
// internal temps
typedef std::map<StrTy,int> FreqTy;

public:
branch_class():m_terminal(0),m_lower(0),m_mom(0) {}
branch_class( TerminalTy & t): m_terminal(&t),m_lower(0),m_mom(0){}
branch_class(LowerTy & l) : m_terminal(0),m_lower(&l),m_mom(0) {}
~branch_class(){ delete m_lower; }
// put line into hierarchy defined by idarray. Each element
// of idarray is an index that specifies a field element in d
// to us as next level key. This is teraminated with ~0. 
// note that this is predicated on unique 
void put(const IdxTy * idarray,  const TerminalTy  *   d)
{
if (*idarray==(~0U)) { m_terminal=d; return; }
const KeyTy & k= (*d)[*idarray];
if (m_lower==NULL) m_lower= new LowerTy();
(*m_lower)[k].put((++idarray),d);

}
// make one line for each first level key and substract seconds, 
template <class Ti> void  walk_group_second(OsTy & os, const Ti & iit, const Ti & iet,const char * opt  )
{
const bool skip_uid=(::strrchr(opt,'u')!=0);
const bool skip_time=(::strrchr(opt,'t')!=0);
const bool make_freq=(::strrchr(opt,'f')!=0);
const bool make_domain=(::strrchr(opt,'d')!=0);

const StrTy sep=" "; 
LoItor li=(*m_lower).begin();
LoItor ie=(*m_lower).end();

while (li!=ie)
{
if (!skip_uid) os<<(*li).first<<sep;
LoItor li2=(*li).second.m_lower->begin();
LoItor ie2=(*li).second.m_lower->end();

if (make_freq)
{
FreqTy fq;
while (li2!=ie2)
{
const KeyTy & k=(*li2).first;
const Myt & v=(*li2).second;
const TerminalTy & t=*(v.m_terminal);
Ti x= iit;
std::stringstream ss;
while (x!=iet)
{
//os<<" x is "<<(*x)<<"for "<< t[*x];
// reall should hve templates to help comiler remove this crap 
if (!make_domain) ss<< t[*x]<<"+";
else 
{
FieldTy u=t[*x];
const char * uc=u.c_str();
char cf[strlen(uc)+10];
strcpy(cf,uc);
char * cfe=strrchr(cf,'.');
if ( cfe!=0) 
{
*cfe=0;
cfe=strrchr(cf,'.');
if (cfe==0) ss<<cf<<"+";
else ss<<(cfe+1)<<"+";
}
else ss<<cf<<"+";
} // make_domainZZhead 
++x;
}
FieldTy kf;
// this is space delimiteed causing problems. 
ss>>kf;
++fq[kf];
++li2;
} // uid records
FreqTy::iterator fi= fq.begin();
os<<"sz="<<fq.size()<<" ";
while (fi!=fq.end())
{
os<<(*fi).first<<"["<<(*fi).second<<"] ";
++fi;
}

}
else
{ // don't make freq
IdxTy n=0;
typedef mjm_strings::my64int rint64;
rint64 t1=0;
rint64 t2=0;
while (li2!=ie2)
{
const KeyTy & k=(*li2).first;
if (!skip_time) 
{
if (n==0) t1=mjm_strings::str_to_time(k.c_str());
t2=mjm_strings::str_to_time(k.c_str());
os<<" "<<(t2-t1)<<sep; 
 } else os<<sep;
const Myt & v=(*li2).second;
const TerminalTy & t=*(v.m_terminal);
Ti x= iit;
while (x!=iet)
{
//os<<" x is "<<(*x)<<"for "<< t[*x];
os<< t[*x];
++x;
if ( x!=iet) os<<sep; 
} // x
++n; 
++li2;
}// li2
} // making freq
os<<CRLF;

++li;
} // itor



}


};


typedef branch_class DocTy;
IdxTy * m_idarray;
IdxTy m_max;
DocTy m_doc;
//typedef std::map<FieldTy,LineTy> LeafTy;
//typedef std::map<FieldTy,LeafTy> TopTy;

public:
mjm_mom() : m_idarray(0),m_max(0)  { }
~mjm_mom() { delete[] m_idarray; } 
// copy key order into an array which we own.
void setOrder(const IdxTy * p, const IdxTy n)
{
if (m_idarray!=0) delete[] m_idarray;
m_idarray= new IdxTy[n+1];
m_max=0;
for(IdxTy i=0; i<n; ++i) 
{ const IdxTy x=p[i]; if (x>m_max) m_max=x; m_idarray[i]=x;}
m_idarray[n]=~0; 

}
// add a "line" in the form that it represents
// an edge with 2 nodes defined only by name.
// this is fine when you don't have const nics, 
void add(const LineTy & l)
{
//const FieldTy & fr_node=l[from];
if ( m_max>=l.size()) throw new std::runtime_error("line is too short missing keys");
if (m_idarray==0) throw new std::runtime_error("no order set in map of maps");
m_doc.put(m_idarray,&l);

} // add

void dump(OsTy & os, IdxTy mode, IdxTy * idxx,IdxTy sz, const char * opt)
{
m_doc.walk_group_second(os,idxx,idxx+sz,opt);
}

private:



};














/*
class mjm_csv
{
private:
typedef mjm_netwerk  Myt;

typedef mjm_generic_traits Tr;



*/





#endif

