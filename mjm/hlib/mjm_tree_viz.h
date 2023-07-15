#ifndef MJM_TREE_VIZ_H__
#define MJM_TREE_VIZ_H__
 
#include "mjm_globals.h"
// TODO this does not belong here but easiest palce to add for now. 
#include "mjm_pheno_notes.h"
#include "mjm_rendering_temp_data.h"
#include "mjm_latex_writer.h"
// for the presence absence vector 
#include "mjm_taxon_tools.h"
#include "mjm_char_mat.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include "mjm_cli_ui.h"
//#include "../mjm_fasta_ii.h"
#include "mjm_fasta_ii.h"


#include "mjm_collections.h"
#include "mjm_svg_writer.h"


//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
//#include "mjm_biom_hdf5.h"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>
#include <stdint.h>

/*
TODO FIXME
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 



////////////////////////////////////////////////////////////////

class tree_viz_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
tree_viz_params( const StrTy & nm) : Super(nm) {}
tree_viz_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
IdxTy maxdepth() const { return m_map.get_uint("maxdepth",3); } // // 100;
bool print_counts() const { return m_map.get_bool("print_counts",!true); }
bool print_haves() const { return m_map.get_bool("print_haves",!true); }
bool print_havenots() const { return m_map.get_bool("print_havenots",!true); }
bool print_if_have() const { return m_map.get_bool("print_if_have",!true); }
bool suppress_vector() const { return m_map.get_bool("suppress_vector",!true); }
bool add_level() const { return m_map.get_bool("add_level",true); }
bool print_hit() const { return m_map.get_bool("print_hit",true); }
IdxTy human_digits() const { return m_map.get_uint("human_digits",1); }
IdxTy data_digits() const { return m_map.get_uint("data_digits",12); }
StrTy human_sep() const { return m_map.get_string("human_sep","\t"); }
StrTy data_sep() const { return m_map.get_string("data_sep"," "); }
StrTy pheno_notes_name() const { return m_map.get_string("pheno_notes_name","z"); }
StrTy pheno_notes_field() const { return m_map.get_string("pheno_notes_field","anaerob"); }
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"accmode="<<accmode()<<sep;
ss<<"maxdepth="<<maxdepth()<<sep;
ss<<"print_counts="<<print_counts()<<sep;
ss<<"print_haves="<<print_haves()<<sep;
ss<<"print_havenots="<<print_havenots()<<sep;
ss<<"print_if_have="<<print_if_have()<<sep;
ss<<"suppress_vector="<<suppress_vector()<<sep;
ss<<"add_level="<<add_level()<<sep;
ss<<"print_hit="<<print_hit()<<sep;
ss<<"human_digits="<<human_digits()<<sep;
ss<<"data_digits="<<data_digits()<<sep;
ss<<"human_sep="<<human_sep()<<sep;
ss<<"data_sep="<<data_sep()<<sep;
ss<<"pheno_notes_name="<<pheno_notes_name()<<sep;
ss<<"pheno_notes_field="<<pheno_notes_field()<<sep;


//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace tree_viz_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
//typedef unsigned int KeyCode;
//typedef unsigned __int128 KeyCode;
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
typedef std::vector<IdxTy > Locations;



//typedef mjm_sparse_matrix<D> MySparse;
}; // 


class layout_consts {
typedef  ::tree_viz_traits::Tr  Tr;
typedef layout_consts Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef mjm_taxon_tools Mtt;
typedef Mtt::ignore_t ignore_t;

public:
layout_consts(const IdxTy _node_count):
 node_count(_node_count),
nodew(fscale(node_count)),
 w(node_count*nodew),
 hscale(3.0/4.0*.5),
 h(w*hscale ), 
 ymargin(.15),
 ylegend(0.02*h ),
 ytaxa0(.3*h ),
ytaxa1(.8*h),
 xoff(0),
 yoff(0 ),
// itaxa(0),
 linesz(.002*h),
 nodesz(rmin(.8,2.5/nodew)),
 visthick(.001*h),
 taxaszconst(bignode(.1*h*.75*nodew,nodew)),
 taxasz(bignode(rmin(5.0*nodew,.04*h*nodew),nodew)),
tw(1.05*w),th((1.0+ymargin)*h)
 {}
static D rmin(const D & x1, const D & x2) { return (x1<x2)?x1:x2;}
static D rmax(const D & x1, const D & x2) { return (x1>x2)?x1:x2;}
D bignode(const D & x1, const D & x2) { return (nodew>1)?x2:x1;}
static IdxTy fscale(const IdxTy & nc) { if (nc>100) return 1;  return 10; }

//ow(x1,y1,np.x1,lc.taxlocedge(right,alt));
D taxlocedge(const bool right, const bool alt ) const
{
if (right&&alt) return .01*h;
if (right&&!alt) return .11*h;
if (!right&&alt) return .99*h;
if (!right&&!alt) return .89*h;
return 0;
}

const ignore_t & ignores() const { return m_mtt.ignores(); }

template <class Tt> 
StrTy informative_name(const  Mtt::tax_t & tv, const Tt & tt) const
{ return  m_mtt.informative_name(tv,tt) ; }
private:
public:
const D node_count,nodew, w , hscale, h , ymargin, ylegend,  ytaxa0, ytaxa1, xoff, yoff;
const D linesz, nodesz, visthick, taxaszconst, taxasz;
const D tw,th;

mjm_taxon_tools m_mtt;

}; // layout_consts

typedef Tr::D D; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::StrTy StrTy  ;  

class wpos { public: wpos():x(0),y(0),color("black") {} D x,y; StrTy color; StrTy nm; };
class npos { public: npos():x0(0),x1(0) {} npos(const D & i ) : x0(i),x1(i+1),ntaxon(0){} D x0,x1; IdxTy ntaxon;  };

class legmap : public std::map<StrTy,wpos>
{
public: 
void add(const StrTy & n, const wpos & w) { (*this)[n]=w;}
class compc{ public: 
static bool cmp(const wpos & x, const wpos & y) { return x.x<y.x; }
bool comp(const wpos & x, const wpos & y) { return x.x<y.x; }
template <class Tfudd > bool operator()(const Tfudd  & x, const Tfudd  & y) { return x.y<y.y; }

};
typedef std::vector<wpos> sorted;
sorted sort()
{
sorted r; 
compc compcfudd;
MM_LOOP(ii,(*this)) { (*ii).second.nm=(*ii).first;  r.push_back((*ii).second); } 
std::sort(r.begin(),r.end(),compcfudd);
return r; 
} 

}; // legmap 

class layout_state
{
typedef  ::tree_viz_traits::Tr  Tr;
typedef layout_state Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;


typedef legmap LegendMap; //  legend_map;


public:
class node_layout : public  std::map<IdxTy,npos>  // NodeLayout;
{
public:
typedef std::vector<IdxTy> Locs;
typedef std::vector<IdxTy> Cursor;
typedef std::vector<Locs> Breaks;
typedef std::vector<IdxTy> Vals;
static const IdxTy & bad() { static const IdxTy v= ~0U;  return v; } 
IdxTy breaks(const IdxTy & loc) const
{
return bad(); 

}
const Cursor&  cursor(const IdxTy & i )
{
return m_breaks[i];
}
// this is a LINEAGE,the first value is lowest level.
void new_values(const Vals & _v, const IdxTy & loc)
{
Vals v;
reverse(v,_v);
const IdxTy newsz=v.size();
const IdxTy oldsz=m_vals.size();
if (newsz>oldsz) add_lvl(bad(),loc,newsz-oldsz); 
MM_SZ_LOOP(i,v,vsz) { 
// lowest level change invalidates lower levels 
if (v[i]!=m_vals[i]) {add_break(i,loc); break; }
}
// this only need copy those that changed... 
m_vals=v;
}
void reverse(Vals & v, const Vals & _v)
{ MM_SZ_LOOP(i,v,sz) { v.push_back(_v[sz-i-1]); } }

void add_lvl(const IdxTy & val,  const IdxTy& loc, const IdxTy n=1 )
{
Locs nv;
if (loc!=bad()) nv.push_back(loc); 
for (IdxTy i=0; i<n; ++i)
{ m_vals.push_back(val);  m_breaks.push_back(nv); }
}

void add_break(const IdxTy&  lvl, const IdxTy& loc)
{
Locs nv;
//nv.push_back(loc);
while (m_breaks.size()<=lvl){ add_lvl(bad(),bad()); }
for (IdxTy i=lvl; i<m_breaks.size(); ++i) m_breaks[i].push_back(loc); 

}

Breaks m_breaks;
Vals m_vals;
}; // node_layout

// this is all in layout_state

typedef std::map<StrTy,wpos> SampleCursor;
//typedef std::map<IdxTy,npos> NodeLayout;
typedef node_layout NodeLayout;
typedef std::map<StrTy,D> SampleFig;
layout_state(): itaxa(0) {}
typedef mjm_svg_writer sw_type;
mjm_svg_writer sw;
typedef NodeLayout nl_type;
typedef NodeLayout::Cursor  LayoutCursor;
NodeLayout nl;
//LayoutCursor nlc;
typedef SampleCursor sc_type;
SampleCursor sc;
typedef SampleFig sv_type;
SampleFig sv;
LegendMap legend_map;
IdxTy itaxa;
}; //layout_state


class node_values  // : public  std::map<IdxTy, NodeValue> NodeValues;
{

typedef node_values Myt;
typedef  ::tree_viz_traits::Tr  Tr;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy,D> NodeValue;
typedef std::vector<StrTy> SampleOrder;
typedef std::vector<IdxTy> NodeOrder;
typedef std::map<IdxTy, NodeValue> NodeValues;
typedef NodeValues::const_iterator nv_iterator;

public:
typedef NodeValue node_value;
typedef std::vector<D> values_t;
typedef SampleOrder sample_order_t;
const node_value & operator[](const IdxTy n) { return m_nv[n];} 
void add(const IdxTy node, const StrTy &sample, const D&  x)
{
m_nv[node][sample]+=x;
if (m_samples.find(sample)==m_samples.end()) m_order.push_back(sample); 
m_samples[sample]+=x;
}
D cvalue(const IdxTy node, const StrTy &sample) const 
{
const auto x=m_nv.find(node);
if (x==m_nv.end())
{
 MM_ERR(" missing node "<<MMPR3(node,sample,m_nv.size()))
 return 0;
}
const auto y=(*x).second.find(sample);
if ( y==(*x).second.end()) return 0;

return (*y).second;
}

IdxTy nsamples() const { return m_samples.size(); } 
const sample_order_t & order() const { return m_order; }  
// for removing non terminal nodes.. 
void move(const IdxTy nn,const IdxTy nt)
{
m_nv[nn]=m_nv[nt];
m_nv[nt].clear();
}
void erase(const IdxTy nt)
{
m_nv.erase(m_nv.find(nt));
}

nv_iterator find(const IdxTy n) { return m_nv.find(n); } 
nv_iterator begin() { return m_nv.begin(); } 
nv_iterator end() { return m_nv.end(); } 

values_t values(const IdxTy n) const
{
values_t v;
auto jj=m_nv.find(n);
// all of nothing, otherwise no diea which are missing 
if (jj==m_nv.end()) return v;
auto&  snv=(*jj).second;
//MM_LOOP(ii,m_samples)
MM_LOOP(oo,m_order)
{
auto ii=m_samples.find(*oo);
auto kk=snv.find((*ii).first);
if (kk==snv.end()) { v.push_back(0); } 
else { v.push_back((*kk).second/(*ii).second) ; } 
//v.push_back(((*jj).second)[(*ii).first]/(*ii).second) ;
}

return v; 
}

void print_values(Ss & ss,const IdxTy & node,const IdxTy & digits
		,const StrTy & sep,const IdxTy & flags)
{
std::vector<D> values=(*this).values(node);
ss.precision(digits);
MM_LOOP(ii,values) {
if ((*ii)==0) ss<<sep<<" 0 "; 
 else { 
const D v=(*ii);
ss<<std::scientific<<sep<<v; }}
}

class sorted_nodes
{
public:
sorted_nodes(): m_values(0),m_ptr(0), m_sample() {}
class sv
{
public:
template <class Ty> 
bool operator()(const Ty & a, const Ty & b) const
{ // TODO FIXME these values need to be normalized for cross sample tests 

auto ii=nvr.find(a);
// TODO not sure this is exactly right but ok for now 
if (ii==nvr.end()) return false; 
auto jj=nvr.find(b);
if (jj==nvr.end()) return true; 
auto kk=(*ii).second.find(m_sample);
if (kk==((*ii).second.end())) return false;
auto ll=(*jj).second.find(m_sample);
if (ll==((*jj).second.end())) return true;

return ((*kk).second<((*ll).second)); 
//return ( nvr[a][m_sample]<nvr[b][m_sample]); 

}
StrTy m_sample;
const NodeValues & nvr;
}; // sv
bool is_valid() const { return m_ptr<m_node_order.size(); } 
NodeValues * m_values;
IdxTy m_ptr;
StrTy m_sample;
NodeOrder m_node_order;

}; // sorted_nodes

NodeValues m_nv;
std::map<StrTy,D> m_samples;
SampleOrder m_order;
}; // node_values


}; // tree_viz_traits



class mjm_tree_viz 
{
typedef  tree_viz_traits::Tr  Tr;
typedef mjm_tree_viz Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef tree_viz_params Logic;
typedef mjm_logic_base VariableStore;


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;



typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;

typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;

// TODO FIXME this needs to be a more comlicated object 
//typedef std::map<StrTy,D> NodeValue;
//typedef std::map<IdxTy, NodeValue> NodeValues;
typedef tree_viz_traits::node_values  NodeValues;
typedef NodeValues::node_value NodeValue;

typedef std::vector<StrTy> TaxVec;

typedef tree_viz_traits::layout_consts LayoutConst;

typedef tree_viz_traits::wpos  wpos;
typedef tree_viz_traits::npos  npos;
typedef tree_viz_traits::legmap  LegendMap;

public:
// TODO FIXME put this some where in FRP ... 
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }
// TODO FIXME find FRP for this.. 
static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

public :
bool done() const  { return m_done; } 
public :
mjm_tree_viz():m_dmel(new Dmel()) {Init();}
mjm_tree_viz(int argc,char **_args) : m_dmel(new Dmel())
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<IdxTy(argc); ++i) args[i]=_args[i];
for (IdxTy i=argc; i<ikluge; ++i) args[i]=&dummy[0];
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
//m_tree.config("-tree",i,argc,args);
//m_flp.config("-params",i,argc,args);
//configi(m_points,"-points",i,argc,args);
//m_flp.config_set("-set-param",  i,  argc, args);
//m_tree.config_set("-set-branch",  i,  argc, args);
cmdlcmd( i, argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
}
~mjm_tree_viz()
{
//clear_handlers();
delete m_dmel;
}
////////////////////////////////////////////////////////
// command block

// this should be in the parameters map, nothing special here... 
 void configi(IdxTy & dest, const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
//const StrTy nm=StrTy(args[i]);
dest=myatoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-source") { ++i; command_modef(args[i]); ++i; }
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { ++i; clean_up(); }
if (s=="-about") { ++i; about(); }
if (confirm) {}
} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 
void start_date(const StrTy & d) { m_flp.set("start_date",d); }
void end_date(const StrTy & d) { m_flp.set("end_date",d); }

void command_modef(const char * fn)
{ std::ifstream fin(fn); CommandInterpretter li(&fin); command_mode(li); }
void command_mode() { CommandInterpretter li(&std::cin); command_mode(li); }
void command_mode(const StrTy & cmd) 
{ CommandInterpretter li; li.set(cmd,1); command_mode(li); }
 void cli_cmd( CliTy::list_type & choices,  const char * frag)
{
//MM_ERR("cli_cmd"<<MMPR(frag))
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}
 void cli_param( CliTy::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
const StrTy cmd=CliTy::word(StrTy(_cmd),0);
auto ii=m_comp_map.find(cmd);
if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag); 
}


//void cmd_write_svg(Cip & cip , LocalVar & lv ) 
//{
//m_char_mat.cmd_write_svg(cip,lv);
//}

void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))

} 

}

//////////////////////////////////

void xxx_cmd_hit_or_miss_fasta(Cip & cip , LocalVar & lv ) 
{ 
//typedef std::map<StrTy, IdxTy> HaveMap;
//typedef  mjm_string_index_collection  Fidc;
//typedef std::vector< string_indexer> Fid;
const StrTy deflin=StrTy("nolineage");
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
} // cmd_hit_or_miss



///////////////////////////////////////////////////////////////



////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
TaxTree & tt = m_tax_tree;
StrTy local_label="tat";
m_cli.set_target(*this);
m_cli.set_command_handler(&Myt::cli_cmd);
m_cli.set_param_handler(&Myt::cli_param);
m_cli.activate();
LocalVar mloc;
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
if (m_flp.log_commands()) 
		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd.c_str()[0]=='#' ) continue; 

if (m_cmd_map.find(cmd)!=m_cmd_map.end())
{
 CommandInterpretterParam  cip(li); 
(this->*m_cmd_map[cmd])(cip,mloc);
continue;

}
const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
if (cmd=="about") { about();  continue; } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="print") { MM_MSG(li.line())  continue; } 
if (cmd=="err") { MM_ERR(li.line())  continue; } 
if (cmd=="status") { MM_STATUS(li.line())  continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 

//m_tax_tree.standard_commnds(cmd,p1,p2,li);
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }

/*
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
*/

if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");
if (m_flp.exit_on_err())
{
MM_ERR(" quiting "<<MMPR(m_flp.exit_on_err()))
clean_up();
return; 

}
}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
void about()
{
Ss ss;
ss<<" mjm_trees_and_tables "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
ss<<" Code to read various files related to 16S rRNA analyses "<<CRLF;
ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;
ss<<"database link provided by Bradley.stevenson@ou.edu"<<CRLF;
ss<<"http://www.earthmicrobiome.org/data-and-code/ "<<CRLF;
ss<<"Sample data provided by echen@zymoresearch.com "<<CRLF;
ss<<"http://www.isppweb.org/smc_files/bull%20et%20al.%202010%20jpp%20list.pdf"<<CRLF;
ss<<"http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.690.647&rep=rep1&type=pdf"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/Pathogenic_bacteria"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/List_of_bacteria_genera"<<CRLF;

std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();

}

// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
// use CommandInterpretterParam 
//const StrTy & cmd=li.word(0);

//if (cmd=="int") { test_int(li);} //   continue; } 
//else if (cmd=="diff") { test_diff(li);} //   continue; } 
//else if (cmd=="fl") { test_fl(li);} //   continue; } 
//void test_gr_integrator( CommandInterpretter & li )
if (false) {}else { MM_ERR(" unrecignized TEST command "<<li.dump()) } 

li.pop();
} // test

void dump_cm()
{
 m_cm.dump("solve_step"  , std::cout);
 MM_MSG("solve_times"<<CRLF<<m_cm.time_table("solver_times"))
}

void config_banner()
{
MM_INC_MSG(m_cm,"test test" )
MM_MSG(" configuration "<<m_flp.to_string())
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general
/*

StrTy tod_from_string(const StrTy & w)
{
if (w.length()!=6) return StrTy("BAD")+w;
return w.substr(0,4);
}

*/



void parse_dmel( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
// do nothing with this for now 
		if (m_dmel!=0) 
		{ (*m_dmel).event("userdmelentry" ,name, d,w); } 

}
////////////////////////////////////////////////////////////////////
void add_node_value( const TaxVec & v, const StrTy & sample, const D & val)
{
TaxTree & tt=m_tax_tree;

NodeValues& nv=  m_values;

const IdxTy node=tt.add_new_nodes(v,0);
//nv[node][sample]+=val;
nv.add(node,sample,val);

}
void tt_operation(const StrTy & cmd, const StrTy & p1) //const
{
m_tax_tree.basic_operations(cmd,p1);
}
void write_tax_tree(OsTy & os,const IdxTy flags=0)
{
m_tax_tree.write_composite(os);
}

// right now only terminal nodes have values 
template <class Tsv,class Tnl, class Tnv>
void  check_nonterminal_nodes(IdxTy & node_count, Tsv & sv, Tnl & nl, Tnv & nv,  TaxTree & tt, const IdxTy&  _node, const IdxTy & taxo )
{
IdxTy node=_node;
auto aii=tt.begin_obj();
auto aee=tt.end_obj();
std::map<IdxTy,IdxTy> nt_nodes;
while (aii!=aee)
{
auto & n=(*aii).second;
if (n.size()==0){ ++aii;  continue; } 
auto ii=nv.find((n.id()));
MM_LOOP(jj,(*ii).second)
{
D dh=(*jj).second;
if (dh!=0)
{
MM_ERR(" non terminal node has counts "<<MMPR4(n.id(),(*jj).first,dh,tt.node_info(n.id(),257)))
MM_ERR(" non terminal node can add fake kids or distribute over phyla   ")
++nt_nodes[n.id()];
}
}
++aii;
}
MM_LOOP(ii,nt_nodes) { 
IdxTy nt=(*ii).first;
IdxTy nn=tt.add_generations(nt,3,StrTy("MISSING"));
//auto ii=nv.find((n.id()));
nv.move(nn,nt);
//nv[nn]=nv[nt];
//nv[nt].clear();

 }
} // check_nonterminal_nodes







template <class Tsv,class Tnl, class Tnv>
void  find_small_groups(IdxTy & ntaxon, IdxTy & node_count, Tsv & sv, Tnl & nl, Tnv & nv,  TaxTree & tt, const IdxTy&  _node, const IdxTy & taxo )
//find_small_groups(node_count,ls.sv,ls.nl,nv,  tt, node,taxo);
{
std::vector<IdxTy> tv;
tt.lineage(tv,_node);
IdxTy tvsz=tv.size();
// see if all DESCENDANTS count is right hit sis kluge TODO FIXME 
IdxTy nkids=(tvsz>=taxo)?tt.nkids(tv[tvsz-taxo]):0;
if (nkids<2) ntaxon=1;
IdxTy ndesc=tv[tvsz-taxo];
while ( nkids==1)
{
ndesc=tt.child(ndesc,0);
//tv.clear(); // std::vector<IdxTy> tv;
//tt.lineage(tv,ndesc);
//tvsz=tv.size();
nkids=tt.nkids(ndesc);
}
if (nkids>1) if ( nkids!=(~0U)) ntaxon=nkids;
}



//template <class Tsv,class Tnl, class Tnv>
template <class Tls, class Tnv>
//void  global_stats(IdxTy & node_count, Tsv & sv, Tnl & nl, Tnv & nv,  TaxTree & tt, const IdxTy&  _node, const IdxTy & taxo )
void  global_stats(IdxTy & node_count, Tls & ls, Tnv & nv,  TaxTree & tt, const IdxTy&  _node, const IdxTy & taxo )
{
IdxTy node=_node;
typename Tls::sv_type & sv=ls.sv;
typename Tls::nl_type & nl=ls.nl;
check_nonterminal_nodes(node_count,sv,nl,nv,  tt, node,taxo);
#if 0
auto aii=tt.begin_obj();
auto aee=tt.end_obj();
std::map<IdxTy,IdxTy> nt_nodes;
while (aii!=aee)
{
auto & n=(*aii).second;
if (n.size()==0){ ++aii;  continue; } 
auto ii=nv.find((n.id()));
MM_LOOP(jj,(*ii).second)
{
D dh=(*jj).second;
if (dh!=0)
{
MM_ERR(" non terminal node has counts "<<MMPR4(n.id(),(*jj).first,dh,tt.node_info(n.id(),257)))
MM_ERR(" non terminal node can add fake kids or distribute over phyla   ")
++nt_nodes[n.id()];
}
}
++aii;
}
MM_LOOP(ii,nt_nodes) { 
IdxTy nt=(*ii).first;
IdxTy nn=tt.add_generations(nt,3,StrTy("MISSING"));
//auto ii=nv.find((n.id()));
nv[nn]=nv[nt];
nv[nt].clear();

 }

#endif


TaxTree::tree_terminal_iterator tis(tt,node);
while (tis.valid())
{ 
node=tis.node();
const IdxTy nkids=tt.nkids(node);
if (nkids!=0) { MM_ERR(" nonterminal node "<<MMPR2(node,tt.node_info(node,256))) } 
// count taxon members, 
IdxTy ntaxon=0;
const bool  find_small_groups=true;
const bool  find_breaks=true;
if ( find_small_groups)
{
(*this).find_small_groups(ntaxon,node_count,sv,nl,nv,  tt, node,taxo);
#if 0 
std::vector<IdxTy> tv;
tt.lineage(tv,node);
IdxTy tvsz=tv.size();
// see if all DESCENDANTS count is right hit sis kluge TODO FIXME 
IdxTy nkids=(tvsz>=taxo)?tt.nkids(tv[tvsz-taxo]):0;
if (nkids<2) ntaxon=1;
IdxTy ndesc=tv[tvsz-taxo];
while ( nkids==1)
{
ndesc=tt.child(ndesc,0);
//tv.clear(); // std::vector<IdxTy> tv;
//tt.lineage(tv,ndesc);
//tvsz=tv.size();
nkids=tt.nkids(ndesc);
}
if (nkids>1) if ( nkids!=(~0U)) ntaxon=nkids;
#endif

}

if (find_breaks)
{
std::vector<IdxTy> tv;
tt.lineage(tv,node);
nl.new_values(tv, node_count);

}


nl[node].ntaxon=ntaxon;
nl[node]=npos(node_count);
 
auto ii=nv.find((node));

MM_LOOP(jj,(*ii).second)
{
D dh=(*jj).second;
sv[(*jj).first]+=dh;
} // jj 
tis.inc();
++node_count;
}
} // global_stats 

IdxTy rooter(IdxTy & taxo, TaxTree & tt)
{
// this does not have a root if bacteria and archea have no common parent 
//const auto rootstt=tt.roots_unify();
//MM_ERR(" trying toots unify again danger will robinson")
const auto rootstt=tt.roots();
if (rootstt.size()==0) { MM_ERR(" no roots ") return bad() ; }
taxo=2;
IdxTy ridx=0;
MM_SZ_LOOP(i,rootstt,szr)
{
//if (tt.node_name(rootstt[i])=="bacteria") { ridx=i; break; } 
// non-parents point here lol 
// but does not propoagate lol 
if (rootstt[i]==bad()) { ridx=i; taxo=2; break; } 
}
IdxTy node=rootstt[ridx];
MM_ERR("Setting roo node "<<MMPR4(node,ridx,tt.node_name(node),tt.node_info(node,256)))
return node;
}


//template < class Np, class Os, class Ow, class Sw,class Lc   >
template < class Np, class Os, class Ow, class Ls,class Lc   >
void mark_taxa_change(const StrTy & t2,Ls & ls,Os & os,Ow& ow,Np & np,const Lc & lc, const bool right)
{
IdxTy & itaxa=ls.itaxa;
typename Ls::sw_type & sw=ls.sw;
const bool alt=((itaxa&1)==0);
const bool mark_central=false;
const bool mark_edge=!false;
D x1,y1;
D tsz=lc.taxasz;
if (np.ntaxon>0) if (np.ntaxon<2) tsz=lc.nodesz*lc.nodew;
if (np.ntaxon>0) if (np.ntaxon==2) tsz=2*lc.nodesz*lc.nodew;
//if (np.ntaxon>0) 
if (np.ntaxon>5) tsz=lc.taxaszconst;
if (mark_edge)
{
//D taxloc=alt?lc.ytaxa0:lc.ytaxa1;
ow(x1,y1,np.x1,lc.taxlocedge(right,alt));
if (right){ 
if (alt) os<<sw.vtext_text(t2,x1,y1,tsz);
else os<<sw.vtext_text(t2,x1,y1,tsz);
}
else{
if (alt)  os<<sw.vtextb_text(t2,x1,y1,tsz);
else   os<<sw.vtextb_text(t2,x1,y1,tsz);
}



}

if (mark_central) {
D taxloc=alt?lc.ytaxa0:lc.ytaxa1;
if (!right) taxloc=lc.h-taxloc;
//if (alt) ow(x1,y1,np.x1,ytaxa0); else ow(x1,y1,np.x1,ytaxa1);
ow(x1,y1,np.x1,taxloc);
//os<<sw.vtext_text(t2,x1,y1,3.0*(x1-x0));
if (right){ 
if (alt) os<<sw.vtext_text(t2,x1,y1,tsz);
else os<<sw.vtextb_text(t2,x1,y1,tsz);
}
else{
if (alt)  os<<sw.vtextb_text(t2,x1,y1,tsz);
else   os<<sw.vtext_text(t2,x1,y1,tsz);
}

} // mark_central



++itaxa;
}




template <class Tv, class Np, class Os, class Ow, class Sw ,class Lc  >
 StrTy  taxon_splits(StrTy & t2,Tv & tv,Tv & tvold,Np & np, Ow & ow, Os & os,Sw & sw,TaxTree & tt,const IdxTy node,const IdxTy taxo, const D & h, const Lc & lc  )
{ 
StrTy taxon=tt.node_name(node);
tv.clear();
tt.lineage(tv,node);
// remove uninformative words
const bool find_bon_mots=true;
if (find_bon_mots)
{
//static std::map<StrTy,IdxTy> ignored;
//typename Lc::ignore_t  ignored=lc.ignores();
#if 0
static bool init=false;
if (!init)
{
ignored[StrTy("uncultured bacterium")]=1;
ignored[StrTy("unculturedbacterium")]=1;
ignored[StrTy("unidentified")]=1;
ignored[StrTy("unclassified")]=1;
ignored[StrTy("blank")]=1;
// these should hve been removed before 
ignored[StrTy("g__")]=1;
ignored[StrTy("s__")]=1;
ignored[StrTy("o__")]=1;
ignored[StrTy("f__")]=1;
ignored[StrTy("c__")]=1;
// place holder for non-terminal counts
ignored[StrTy("MISSING")]=1;
ignored[StrTy("na")]=1;

init=true; 
}
#endif

taxon=lc.informative_name(tv,tt);
//taxon="";
/*
const StrTy sep=" ";
for (int  i=0; i<tv.size(); ++i) 
{
const StrTy & w=tt.node_name(tv[i]);
if (ignored.find(w)==ignored.end()) 
	{if (taxon.length()!=0) { taxon=w+sep+taxon; break; } taxon=w; } 

}
*/

}


const IdxTy szlim=taxo;
const IdxTy sztv=tv.size();
const IdxTy sztvold=tvold.size();
if (!find_bon_mots) if (sztv>1) if (sztvold>1) if (tv[1]!=tvold[1]) taxon=taxon+StrTy(" ")+tt.node_name(tv[1]);
const bool case1=(sztv>szlim) && (sztvold>szlim) &&(tv[sztv-taxo]!=tvold[sztvold-taxo]);
const bool case2=(sztv>szlim) && (sztvold<=szlim);
if (case1||case2)
{
D x0,y0,x1,y1;
ow(x0,y0,np.x0,0);
ow(x1,y1,np.x0,1.0*h);

os<<sw.line_text(x0,y0,x1,y1, StrTy("blue"),.1*lc.nodew); os<<CRLF;
 //taxon=taxon+tt.node_name(tv[tv.size()-2]);
 t2=tt.node_name(tv[tv.size()-taxo]);
}
return taxon;
}
template <class Sc, class Lm > 
wpos & get_cursor(const StrTy & sname,IdxTy & snamed,Sc & sc,Lm & legend_map,IdxTy & node,IdxTy & nidx,const IdxTy &node_count)
{ 
//const StrTy & sname=(*jj).first;
wpos & cp=sc[sname];
if (cp.x==0)
{
const IdxTy ci=snamed%5;
if (ci==1) cp.color="red";
if (ci==2) cp.color="green";
if (ci==3) cp.color="blue";
if (ci==3) cp.color="white";
legend_map.add(sname,cp);
//if (snamed>20) { MM_ONCE(" breaking due to too many samples ",) break; }
++snamed; 
}
// add a line 
if (nidx>(node_count>>1)) if  (legend_map[sname].x==0) { legend_map[sname]=cp; } 
return cp;
}

template < class Np, class Os, class Ow, class Lc,  class Ls  >
//void label_node( Os & os,Sw& sw,Ow& ow,Lc & lc,const bool right,Np & np,const StrTy & taxon, Ls & ls)
void label_node( Os & os,Ow& ow,Lc & lc,const bool right,Np & np,const StrTy & taxon, Ls & ls, const IdxTy nidx)
{  
typename Ls::sw_type & sw=ls.sw;
const typename Ls::LayoutCursor & laycur=ls.nl.cursor(0);

const bool label_on_empty_edge=false;
const bool label_middle=true;
D x0,y0,x1,y1;
// TODO FIXME the alignement crap is wrong??? 
// kluge TODO FIXME this is an offset problem... 
D dx=(np.x1-np.x0)/lc.nodew;
MM_ONCE(" check fing bseline f "<<dx,)
if (label_on_empty_edge)
{
D yi=(right)?lc.ylegend:(lc.h-lc.ylegend);
ow(x0,y0,np.x0+dx,yi);
ow(x1,y1,np.x1+dx,yi);
//os<<sw.vtext_text(taxon,np.x0,ylegend,(np.x1-np.x0));
if (right) { os<<sw.vtext_text(taxon,x0,y0,lc.nodesz*(x1-x0)); } 
else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));} 
}
if (label_middle)
{
D yimax=lc.ylegend;
D yimin=(lc.h-lc.ylegend);
D yimaxc=yimin;
D yiminc=yimax;
const D mid=.5*(yimax+yimin);

D yi=0;
D ysum=0; IdxTy n=0;
MM_LOOP(ii,ls.sc)
{
if (ls.sv[(*ii).first]==0)
{
MM_ERR(" ho hits "<<(*ii).first) 
continue;
}
const D y=(*ii).second.y;
ysum+=y;
++n;
if (y>yimaxc) yimaxc=y;
if (y<yiminc) yiminc=y;
//MM_ERR(MMPR4(y,yimaxc,yiminc,yimax))
}
if (yimaxc>yimax) yimaxc=yimax;
if (yiminc<yimin) yiminc=yimin;
yi=ysum/n; // .5*(yimaxc+yiminc);
const bool align_top=(yi<mid);
const D  len=taxon.length()*lc.nodesz;
if (align_top) if (len>(lc.h-yi)) yi=.01*lc.h;
if (!align_top) if (len>(yi)) yi=.98*lc.h;
ow(x0,y0,np.x0+dx,yi);
ow(x1,y1,np.x1+dx,yi);
if (align_top) { os<<sw.vtext_text(taxon,x0,y0,lc.nodesz*(x1-x0)); } 
else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));} 

}

ow(x0,y1,np.x0,lc.h);
ow(x0,y0,np.x0,0);

os<<sw.line_text(x0,y0,x0,y1, StrTy("blue"),.01*lc.nodew); os<<CRLF;
MM_ONCE(" unused highlight_column dusabed ",)
const bool highlight_column =false; // (np.x0>30)&&(np.x0<50);
if (highlight_column)
{
ow(x0,y0,np.x1,0);
// x,y,w,h, doh !
os<<sw.shade_rect_text(x0,y0,x0,y1, StrTy("red"),.5); os<<CRLF;
}



}

//typedef  ::tree_viz_traits::Tr  Tr;
typedef tree_viz_traits::layout_state LayoutState;


template <class Os, class Ow, class Np, class Lc, class Ls>
void plateau_and_slope(Os & os, Ow& ow,Np & np, const StrTy & sn, const D & counts,Lc & lc,Ls & ls,wpos &cp)
{
// plateau and slope 
const D totals=ls.sv[sn];
if (totals==0) { MM_ERR(" total count zero for "<<sn) return; }
D dh=1.0*lc.node_count*counts/totals*lc.hscale*lc.nodew;
D i=cp.x-lc.xoff; // nidx;
D j=cp.y-lc.yoff; // sc[(*jj).first];
D i1=np.x0-lc.xoff;
D j1=j; // -yoff;
D fi=np.x1-lc.xoff;
D fj=j+dh;
cp.x=fi;
cp.y=fj;
//nv[node][sample]+=val;
D x0,y0,x1,y1,x2,y2;
ow(x0,y0,i,j);
ow(x1,y1,i1,j1);
ow(x2,y2,fi,fj);
//MM_ERR(MMPR((*jj).first)<<MMPR4(i,j,fi,fj)<<MMPR4(x0,y0,x1,y1))
//plateau 
os<<ls.sw.line_text(x0,y0,x1,y1, cp.color,lc.linesz); os<<CRLF;
// slope
os<<ls.sw.line_text(x1,y1,x2,y2, cp.color,lc.linesz); os<<CRLF;
}

void fix_layout(LayoutState & ls,const LayoutConst & lc)
{
MM_LOOP(ii,ls.nl) { (*ii).second.x1*=lc.nodew; (*ii).second.x0*=lc.nodew; }


}

class page_replay
{
public:
template <class Ty> 
page_replay( const Ty & v, const IdxTy n, const IdxTy d)
:values(v),node(n),depth(d) {}
page_replay() {}
std::vector<D> values;
IdxTy node;
IdxTy depth;
}; // page_replay

void write_latex(std::ostream & os, const IdxTy flags)
{
const bool include_docs=true;
mjm_latex_writer lw;
NodeValues& nv=  m_values;
TaxTree & tt=m_tax_tree;
tt.sort_for_ui();
IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
IdxTy nidx=0;
if (include_docs) { os<<lw.start_text();}
const StrTy caption="test table";
const StrTy label="test:table";
const IdxTy depth=10;
const IdxTy nsamp=nv.nsamples();
std::vector<page_replay> pb(depth+1);
const auto & order=nv.order();
os<<lw.tree_table_start_text( caption,  label,  depth, nsamp);
os<<lw.tree_table_header_text(StrTy("taxa"),order,depth,nsamp,4);
IdxTy mind=depth;
IdxTy nlines=0;
TaxTree::tree_hierarch_iterator ti(tt,node);
while (ti.valid())
{ 
++nlines;
node=ti.node();
//Ss ss;
// get node values

std::vector<IdxTy> tv;
tt.lineage(tv,node);
IdxTy tvsz=tv.size();
const IdxTy bumps=ti.bumps();
//StrTy indent=" ";
//for( IdxTy i=0; i<ti.depth(); ++i) indent+=StrTy("    ");
//ss<<ti.dump();
std::vector<D> values=nv.values(node);
if (nidx>2)// if (ti.depth()<4) { 
if ((nlines>30) ||(ti.depth()<4)) { 
os<<lw.tree_table_end_text();
os<<lw.tree_table_start_text( caption,  label,  depth, nsamp);
if (!false){ 
os<<lw.tree_table_header_text(StrTy("taxa"),order,depth,nsamp,4);
for(IdxTy i=mind; i<ti.depth(); ++i)
os<<lw.tree_table_text(lw.escape_dash(tt.node_name(pb[i].node)),pb[i].values,pb[i].depth,depth,nsamp,4);
}
nlines=ti.depth()-mind;
}
const StrTy nn=tt.node_name(node);
if (nn!="!") { //  continue; 
pb[ti.depth()]=page_replay(values,node,ti.depth());
if (ti.depth()<mind) mind=ti.depth();
//MM_LOOP(ii,values) { ss<<" "<<(*ii); } 
//os<<indent<<tt.node_name(node)<<" "<<ss.str()<<CRLF;
os<<lw.tree_table_text(lw.escape_dash(nn),values,ti.depth(),depth,nsamp,4);
}
// StrTy tree_table_text(const StrTy & w, const Ty & v, const IdxTy level,const IdxTy depth, const IdxTy vals, const IdxTy precision) const

++nidx;
ti.inc();
} // ti
os<<lw.tree_table_end_text();
if (include_docs) { os<<lw.end_text();}
} // write_latex

void write_txt(std::ostream & os, const IdxTy flags)
{
Logic & flp=m_flp;
StrTy spacerstring=" ";
StrTy hierlabel=" Hierarchy   ";
const bool human=true;
StrTy sep=human?flp.human_sep():flp.data_sep(); // "\t";
IdxTy digits=human?flp.human_digits():flp.data_digits(); // 1; // :8;

NodeValues& nv=  m_values;
TaxTree & tt=m_tax_tree;
tt.sort_for_ui();
IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
IdxTy nidx=0;
TaxTree::tree_hierarch_iterator ti(tt,node);
const auto & order=nv.order();
{
Ss ss;
ss<<hierlabel;
MM_LOOP(ii,order) ss<<" "<<(*ii);
os<<ss.str()<<CRLF;
}

while (ti.valid())
{ 
node=ti.node();
Ss ss;
std::vector<IdxTy> tv;
tt.lineage(tv,node);
IdxTy tvsz=tv.size();
const IdxTy bumps=ti.bumps();
StrTy indent=" ";
for( IdxTy i=0; i<ti.depth(); ++i) indent+=spacerstring; // StrTy(" ");
//ss<<ti.dump();
nv.print_values(ss,node,digits,sep,flags);
/*
std::vector<D> values=nv.values(node);
ss.precision(digits);
MM_LOOP(ii,values) {
if ((*ii)==0) ss<<sep<<" 0 "; 
 else { ss<<std::scientific<<sep<<(*ii); }}
*/
 
StrTy nname=tt.node_name(node);
while (nname.length()<8) {nname=nname+StrTy(" "); } 
os<<indent<<nname<<" "<<ss.str()<<CRLF;

++nidx;
ti.inc();
} // ti
} // write_txt

template <class Tdnc, class Tt,class Tnv, class Toracle,class Lc >
void collect_nodes(Tdnc & dnc, Tt & tt,Tnv & nv, Toracle &oracle, Lc & lc  )
{IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
// this really should get ALL nodes and let user prune thgm 
TaxTree::tree_terminal_iterator ti(tt,node);
while (ti.valid())
{ 
node=ti.node();
auto ii=nv.find((node));
// samples need to be there. 
if ((*ii).second.size()==0) { ti.inc(); continue; } 
std::vector<StrTy> tv,toracle;
tt.lineage(tv,node);
oracle.lookup(toracle,tv);
// this is a fudding map from sample to value fudd 
//MM_LOOP(jj,(*ii).second) { } // jj 
//dnc.add(tv,(*ii).second);
MM_ONCE(" adding notes even if there are none",)
if (toracle.size()!=0)
{ MM_ERR(MMPR4(node, toracle.size(),toracle[0],tt.node_name(node))) } 
//dnc.add(tv,(*ii).second,toracle);
if ( toracle.size()==0) 
 dnc.add(tv,(*ii).second,toracle); // ,(toracle.size()!=0),0);
// #error  fix this here 

MM_SZ_LOOP(it,toracle,szto)
{
if (it==0)  dnc.add(tv,(*ii).second,toracle,(toracle.size()!=0),lc.flag_number(toracle[0]));
else dnc.flag_last_node(lc.flag_number(toracle[it]));
}


dnc.set_last_src(dnc.st()(tt.st()(tt.node(node).src())));
ti.inc();
} // ti 

} // collect_nodes
/*
template <class Tlay, class Tdnc>
void analyze_layout(Tlay & lay, Tdnc & dnc)
{
// typedef typename Tdnc::iterator Ii;
MM_LOOP(ii,dnc)
{


} // ii 

} // analyze_layout 
*/

template < class Tnv, class Lc >
void  check_nonterminal_nodes_ii(Tnv & nv,  TaxTree & tt, Lc & lc )
{
auto aii=tt.begin_obj(); auto aee=tt.end_obj();
std::map<IdxTy,IdxTy> nt_nodes;
while (aii!=aee)
{
auto & n=(*aii).second;
if (n.size()==0){ ++aii;  continue; } 
auto ii=nv.find((n.id()));
MM_LOOP(jj,(*ii).second)
{
D dh=(*jj).second;
if (dh!=0)
{
if (n.id()!=bad())
{
if (false)
{MM_ERR(" non terminal node has counts "<<MMPR4(n.id(),(*jj).first,dh,tt.node_info(n.id(),257)))
MM_ERR(" non terminal node can add fake kids or distribute over phyla   ")
} //  false
MM_ONCE(" bad nodes ",)
++nt_nodes[n.id()];
} 
else { MM_ONCE(" non terminal bad node ignored ",) }
}
}
++aii;
}
MM_LOOP(ii,nt_nodes) { 
IdxTy nt=(*ii).first;
//IdxTy nn=tt.add_generations(nt,3,StrTy("MISSING"));
IdxTy nn=tt.add_generations(nt,1,lc.terminal_name());
nv.move(nn,nt);
//nv.erase((nt));
 }
}



template <class Os, class Ow, class Sw, class La, class Lc>
void svg_ii_start(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc)
{
D xs,ys;
ow(xs,ys,lc.total_width(),lc.total_height());
os<<sw.start_text(lc.doc_name(),xs,ys);
os<<CRLF;
os<<sw.frame_text(lc.bg_color(),xs,ys);
os<<CRLF;
os<<sw.stroke_text("red",1000);
os<<CRLF;

}
template <class Os, class Ow, class Sw, class La, class Lc>
void svg_iia_start(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc)
{
D xs,ys;
ow(xs,ys,lc.total_iia_width(),lc.total_iia_height());
os<<sw.start_text(lc.doc_name(),xs,ys);
os<<CRLF;
os<<sw.frame_text(lc.bg_color(),xs,ys);
os<<CRLF;
os<<sw.stroke_text("red",1000);
os<<CRLF;
}
template <class Os, class Ow, class Sw, class La, class Lc>
void svg_hm_start(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc)
{
D xs,ys;
ow(xs,ys,lc.total_hm_width(),lc.total_hm_height());
os<<sw.start_text(lc.doc_name(),xs,ys);
os<<CRLF;
os<<sw.frame_text(lc.bg_color(),xs,ys);
os<<CRLF;
os<<sw.stroke_text("red",1000);
os<<CRLF;
}





template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_flag_legend(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
typedef rendering_temp_data::legend_layout Ll;
Ll ll(lc.flag_colors(), lc.flag_names());
D rsize=.1*lc.bottom_margin();
D x=lc.plot_x();
D y=lc.plot_y()+lc.plot_height()+lc.bottom_margin()-rsize;
D w=lc.plot_width();
D h=rsize;// lc.bottom_margin();
ll.layout(x,y,w,h);
const IdxTy sz=lc.nflags();
for(IdxTy i=0; i<sz; ++i)
{
StrTy text,color,size;
D angle;
//void label(StrTy & text, StrTy & size, D & angle, const IdxTy i) const
//void symbol( D & x, D & y, D & w , D & h, StrTy &color, const IdxTy i ) const
ll.symbol( x, y, w , h, color, i );
os<<sw.shade_rect_text(x,y,w,h, color,1); os<<CRLF;
ll.label(x,y, text, size, angle, i); 
//os<<sw.htext_text(text,x,y,atof(size.c_str()),StrTy("white"),"middle"); 
os<<sw.htext_text(text,x,y,atof(size.c_str()),StrTy("white"),"left"); 
os<<CRLF;
//os<<sw.shade_rect_text(x,y,w,h, c,1); os<<CRLF;
//os<<sw.htext_text(r.m_name,xpos,y0,10,StrTy("white"),"middle"); os<<CRLF;
} 
} // svg_ii_flag_legend
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_flag_legend(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
os<<sw.comment_text(" iia_flag_lengend" ) ; 
//if (true) return; 
//os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.node_iia_color(ii.i()),1); os<<CRLF;
typedef rendering_temp_data::data_node_collection Nodes;
typedef Nodes::iterator node_iterator;
typedef rendering_temp_data::legend_layout Ll;
typedef std::vector<StrTy> Nv;
Nv node_colors,node_names;
const node_iterator ee=dnc.end();
for (node_iterator ii=dnc.begin(); ii!=ee; ++ii) 
{
//os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.node_iia_color(ii.i()),1); os<<CRLF;
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
node_names.push_back(taxon);
node_colors.push_back(lc.node_iia_color(ii.i()));
}

//Ll ll(lc.flag_colors(), lc.flag_names());
Ll ll(node_colors, node_names);
D rsize=.1*lc.bottom_margin();
D x=lc.plot_x();
D y=lc.plot_y()+lc.plot_height()+lc.bottom_margin()-rsize;
D w=lc.plot_width();
D h=rsize;// lc.bottom_margin();
ll.layout(x,y,w,h);
const IdxTy sz=lc.nflags();
for(IdxTy i=0; i<sz; ++i)
{
StrTy text,color,size;
D angle;
//void label(StrTy & text, StrTy & size, D & angle, const IdxTy i) const
//void symbol( D & x, D & y, D & w , D & h, StrTy &color, const IdxTy i ) const
ll.symbol( x, y, w , h, color, i );
if (h<1) h=1;
if ( w<1) w=1;
os<<sw.shade_rect_text(x,y,w,h, color,1); os<<CRLF;
ll.label(x,y, text, size, angle, i); 
D tsz=atof(size.c_str());
if (tsz<1) tsz=1;
//os<<sw.htext_text(text,x,y,atof(size.c_str()),StrTy("white"),"middle"); 
os<<sw.htext_text(text,x,y,tsz,StrTy("white"),"left"); 
os<<CRLF;
//os<<sw.shade_rect_text(x,y,w,h, c,1); os<<CRLF;
//os<<sw.htext_text(r.m_name,xpos,y0,10,StrTy("white"),"middle"); os<<CRLF;
} 
} // svg_iia_flag_legend





template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_cats(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
 typedef  typename La::Runs Ru;
 typedef  typename La::run R;
const Ru & rup=la.m_cat_runs["phylum"];
IdxTy ij=0;
MM_LOOP(ii,rup)
{
const R & r=(*ii);
if (r.len()>30)
{
const IdxTy i=r.m_first; 
const IdxTy e=r.m_last; 
D  mean=la.m_means[i];
D x0=lc.plot_x()+lc.node_size()*i;
D x1=lc.plot_x()+lc.node_size()*(e+1);
//D x1=x0; // lc.plot_x()+lc.node_size()*(ii.i()+1);
//D y0=lc.plot_y();
D ypos=.1;
if (mean<.5) ypos=1.0-ypos;
D y0=lc.plot_y()+lc.plot_height()*ypos;
D y1=lc.plot_y()+lc.plot_height()*ypos;
ow(x0,y0);
ow(x1,y1);
D xpos=x0+(-x0+x1)/2.0;
const StrTy color=(ij&1)?"red":"blue";
os<<sw.line_text(x0,y0,x1,y1, color,2+.5*lc.node_size()); os<<CRLF;
os<<sw.htext_text(r.m_name,xpos,y0,10,StrTy("white"),"middle"); os<<CRLF;
++ij;
}


} // ii 
const Ru & rup2=la.m_cat_runs["pseudogenus"];
MM_LOOP(ii,rup2)
{
const R & r=(*ii);
if (r.len()>10)
{
const IdxTy i=r.m_first; 
const IdxTy e=r.m_last; 
D  mean=la.m_means[i];
D  maxx=la.m_maxxs[i];
D  minn=la.m_minns[i];
D x0=lc.plot_x()+lc.node_size()*i;
D x1=lc.plot_x()+lc.node_size()*(e+1);
//D x1=x0; // lc.plot_x()+lc.node_size()*(ii.i()+1);
//D y0=lc.plot_y();
D ypos=.2;
if (mean<.5) ypos=1.0-ypos;
if (minn>(1.0-maxx)) { ypos=.75*minn; }
else { ypos=1.0-.75*(1.0-maxx); }
if ( ypos<.2) ypos=.2;
if (ypos>.8) ypos=.8;
D y0=lc.plot_y()+lc.plot_height()*ypos;
D y1=lc.plot_y()+lc.plot_height()*ypos;
ow(x0,y0);
ow(x1,y1);
D xpos=x0+(-x0+x1)/2.0;
const StrTy color=(ij&1)?"orange":"yellow";
os<<sw.line_text(x0,y0,x1,y1, color,2+.5*lc.node_size()); os<<CRLF;
//os<<sw.htext_text(r.m_name,xpos,y0,10,StrTy("white")); os<<CRLF;
os<<sw.vtext_text(r.m_name,xpos,y0,5,color, "middle"); os<<CRLF;
++ij;
}

} // ii 

}


template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_curves(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

typedef rendering_temp_data::data_node_collection Nodes;
typedef Nodes::iterator node_iterator;

MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
/// this sample is an m_st() index value NOT an ordingal 
const IdxTy sample=la.m_sample_order[j];
node_iterator ee=dnc.end();
std::vector<D> xv,yv;
const D xorg=lc.plot_x();
const D yorg=lc.plot_y(); 
xv.push_back(xorg);
yv.push_back(yorg);
std::vector<IdxTy> plateau;
const StrTy plab=lc.sample_name(j);
const IdxTy pspace=2*plab.length()/lc.node_size()+2;
const IdxTy plen=plab.length()/lc.node_size();
for (node_iterator ii=dnc.begin(); ii!=ee; ++ii) 
{
D y=la.m_int[sample][ii.i()]*lc.plot_height()+yorg; 
D x=xorg+lc.node_size()*(ii.i()+1);
ow(x,y);
xv.push_back(x);
yv.push_back(y);
bool next_plateau=(plateau.size()==0);
if (!next_plateau) next_plateau=(ii.i()>(plateau.back()+pspace));
if (next_plateau) if (yv.size()>(plen+5)) if (yv[ii.i()]==yv[ii.i()-plen]) plateau.push_back(ii.i()-plen);
} // node 
const bool label_plateau=true;
if (label_plateau) { 
MM_LOOP(ii,plateau)
{
const IdxTy idx=(*ii);
os<<sw.htext_text(plab,xv[idx],yv[idx],2,lc.sample_contrast_color(j)); os<<CRLF;
}
}

//MM_ERR(MMPR2(xv.size(),yv.size()))
//os<<sw.vector_text(xv,yv,lc.sample_color(sample),lc.curve_width());
// make partially transporaent 
os<<sw.vector_text(xv,yv,lc.sample_color(j),lc.curve_width(),.5);
os<<CRLF;
} // j sample...

} // curves
/////////////////////////////////////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_rect(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
os<<sw.comment_text(" iia_rect" ) ; 
typedef rendering_temp_data::data_node_collection Nodes;
typedef Nodes::iterator node_iterator;

MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
node_iterator ee=dnc.end();
//std::vector<D> xv,yv;
const D xorg=lc.sample_iia_x(j); // plot_x();
const D dxorg=lc.sample_iia_filldx(j); // plot_x();
const D yorg=lc.plot_y(); 
D x0=xorg;
D y0=yorg;
for (node_iterator ii=dnc.begin(); ii!=ee; ++ii) 
{
D y1=la.m_int[sample][ii.i()]*lc.plot_height()+yorg; 
D y1a=y1;
x0=xorg; 
D x1=xorg+dxorg;
ow(x0,y0);
//D x1=xorg+lc.node_size()*(ii.i()+1);
if (y1>y0) {
ow(x1,y1);
os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.node_iia_color(ii.i()),1); os<<CRLF;
}
y0=y1a;
} // node 
} // j sample...
} // _iia_rect  ex curves
/////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_hm_rect(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
os<<sw.comment_text(" hm_rect" ) ; 
typedef rendering_temp_data::data_node_collection Nodes;
typedef Nodes::iterator node_iterator;
const D w=lc.node_size(); // 10;
const D h=lc.plot_height()/la.m_sample_order.size(); // 10;
D zmax=-1;
D zmin=-1;
std::vector<D> xv,yv;
std::vector<StrTy> cv;
for (IdxTy pass=0; pass<2; ++pass)
{
MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
const IdxTy jold=j; // sample;
node_iterator ee=dnc.end();
//std::vector<D> xv,yv;
const D xorg=lc.sample_iia_x(jold); // plot_x();
const D dxorg=lc.sample_iia_filldx(jold); // plot_x();
const D yorg=lc.plot_y(); 
//D x0=xorg;
//D y0=yorg;
D z0=0;
for (node_iterator ii=dnc.begin(); ii!=ee; ++ii) 
{
D z1=la.m_int[sample][ii.i()]; 
if (pass==0)
{
D dz=z1-z0;
z0=z1;
if (dz==0) continue;
if (zmax<0) zmax=dz;
if (zmax<dz) zmax=dz;
if (zmin<0) zmin=dz;
if (zmin>dz) zmin=dz;

continue;
}
//D z1a=z1;
//x0=xorg; 
//D x1=xorg+dxorg;
//ow(x0,y0);
//D x1=xorg+lc.node_size()*(ii.i()+1);
if (pass==1)if (z1>z0) {

D dz=z1-z0;
z0=z1;
D val=-log(dz)/log(10);
/*
D rho=(val-zmin)/(zmax-zmin);
D red=1.0-2.0*rho;
if (red<0) red=0;
D blue=2.0*rho-1.0;
if (blue<0) blue=0;
D green=(rho-.75)*(rho-.25);
if (green>0) green=0;
else { green=sqrt(-16.0*green); } 
color(col,red,green,blue,0);
*/
StrTy col= color(val,zmin,zmax);


D y1= lc.hm_sample_y(jold); // j*h+lc.plot_y();
D xsugg=lc.node_start_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
D x1=xsugg; // ii.i()*w+lc.plot_x();

ow(x1,y1);
//os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.node_iia_color(ii.i()),1); os<<CRLF;
xv.push_back(x1);
yv.push_back(y1);
cv.push_back(col);

}
//z0=z1a;
} // node 
} // j sample...
if (pass==0) { zmin=-log(zmin)/log(10); } 
if (pass==0) { zmax=-log(zmax)/log(10); } 
} // pass

os<<sw.sparse_shade_rect_array_text(xv,yv,w,h, cv,1); os<<CRLF;

xv.clear();
yv.clear();
cv.clear();
D yleg=lc.plot_height()+lc.plot_y()+10;
D xleg=lc.plot_x();
D wleg=lc.plot_width();
const IdxTy nblocks=10;
D legp=wleg/(nblocks);
D hleg=.25*lc.bottom_margin()+5; // lc.w();
D ssz=15;
MM_ONCE(" fudding legend is fudded ",) 
for (IdxTy i=0; i<nblocks; ++i)
{
xv.push_back(xleg+legp*i);
yv.push_back(yleg);
const D val=1.0*i/(nblocks-1)*(zmax-zmin)+zmin ; // (1.0*i/16 -zmin)/(zmax-zmin);
StrTy col= color(val,zmin,zmax);
cv.push_back(col);
Ss ss;
ss<<std::setprecision(3)<<(val); // (val*(zmax-zmin)+zmin); 
os<<sw.gtext_text(ss.str(),xv.back()+.5*legp,yv.back()+hleg+10,ssz,cv.back(),StrTy("left"),45); os<<CRLF;

}
//os<<sw.sparse_shade_rect_array_text(xv,yv,wleg/nblocks,hleg, cv,1); os<<CRLF;
os<<sw.sparse_shade_rect_array_text(xv,yv,legp,hleg, cv,1); os<<CRLF;


os<<sw.gtext_text(StrTy("marchywka test code"),xleg,yleg+.5*hleg,hleg,StrTy("black"),StrTy("left"),0); os<<CRLF;

} // _iia_rect  ex curves

StrTy color(const D & val, const D & zmin, const D & zmax) const
{
D rho=(val-zmin)/(zmax-zmin);
D red=1.0-2.0*rho;
if (red<0) red=0;
D blue=2.0*rho-1.0;
if (blue<0) blue=0;
D green=(rho-.75)*(rho-.25);
if (green>0) green=0;
else { green=sqrt(-16.0*green); } 
StrTy col="";
color(col,red,green,blue,0);
return col;

}

const char  lut(const IdxTy p) const
{
static char c[16];
static bool init=false;
if (!init)
{
for(IdxTy i=0; i<10; ++i) c[i]='0'+i;
for(IdxTy i=0; i<6; ++i) c[i+10]='A'+i;

}
return c[p];
}
StrTy myhex(const D & x) const 
{
Ss ss;
const IdxTy n=IdxTy(x*255);
ss<<(lut((n>>4)&15));
ss<<(lut((n>>0)&15));
return ss.str();
}
void color(StrTy & col,const D & red,const D & green,const D & blue,const IdxTy & flags)
const 
{
Ss ss;
ss<<"#0";
ss<<myhex(red);
ss<<myhex(green);
ss<<myhex(blue);
col=ss.str();
}


/////////////////////////////////////////
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_lines(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
typedef rendering_temp_data::data_node_collection Nodes;
typedef Nodes::iterator node_iterator;
const IdxTy sosz=la.m_sample_order.size();
for (IdxTy j=1;j<sosz; ++j)
{
const IdxTy sample0=la.m_sample_order[j-1];
const IdxTy sample1=la.m_sample_order[j];
D v0last=0;
D  vflast=0;
node_iterator ee=dnc.end();
for (node_iterator ii=dnc.begin(); ii!=ee; ++ii) 
{
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
const D xorg=lc.sample_iia_x(j-1); // plot_x();
D x1=lc.sample_iia_x(j); // plot_x();
const D dxorg=lc.sample_iia_filldx(j-1); // plot_x();
const D yorg=lc.plot_y(); 
D x0=xorg +dxorg;
const D v0=la.m_int[sample0][ii.i()];
const D vf=la.m_int[sample1][ii.i()];
const bool brk=( la.m_la_break_levels[ii.i()]<3);
// otherwise they overlaps and the colors are confusing 
const bool too_small=( ((vf-vflast)<.02) && ((v0-v0last)<.02) ) ; // continue;
vflast=vf;
v0last=v0;
if (too_small&&!brk) continue;

D y0=v0*lc.plot_height()+yorg; 
D y1=vf*lc.plot_height()+yorg; 
const D dx=x1-x0;
const D dy=y1-y0;
const D len=sqrt(dx*dx+dy*dy);
D xmid=x0+.1*(dx);
D ymid=y0+.1*(dy);
D angle=180.0/M_PI*atan2(dy,dx);
xmid+=dy/len;
ymid-=dx/len;
{
ow(x0,y0);
ow(x1,y1);
ow(xmid,ymid);
D fsz=.5*lc.node_size();
if (brk) fsz=fsz*5;
const StrTy colr= lc.node_iia_color(ii.i());
//os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.node_iia_color(ii.i()),1); os<<CRLF;
//os<<sw.line_text(x0,y0,x1,y1, lc.node_iia_color(ii.i()),fsz); os<<CRLF;
os<<sw.line_text(x0,y0,x1,y1, colr,fsz); os<<CRLF;
if (brk) fsz=fsz*.5;
os<<sw.gtext_text(taxon,xmid,ymid,5*fsz,colr,StrTy("left"),angle); os<<CRLF;
}
//y0=y1a;
} // node 
} // j sample...
} // _iia_lines




/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////

//svg_ii_legend(os,ow,sw,la,lc,dnc);
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_legend(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
typedef std::vector<D> Coord;
Coord x,y;
lc.sample_legend_pos(x,y,la.m_sample_order);
MM_SZ_LOOP(i,la.m_sample_order,sosz)
{
D x0=x[i];
D y0=y[i];
ow(x0,y0);
const IdxTy j=la.m_sample_order[i];
D x1=lc.leg_x(j);
D y1=lc.leg_y(j);
ow(x1,y1);
// TODO the subscript is ordinal NOT m_st()
const StrTy & colr=lc.sample_color(i);
if (false) { MM_ERR(" svg_ii_legend "<<MMPR4(i,j,x0,y0)<<MMPR3(x1,y1,colr))}
//const StrTy & sname=lc.sample_name(i);
//os<<sw.htext_text((*ii).first,x1,y1,fitsz,colr); os<<CRLF;
//os<<sw.htext_text(sname,x0,y0,lc.sample_legend_font(),colr); os<<CRLF;
os<<sw.faintline_text(x0,y0,x1,y1, colr,lc.sample_ptr_size()); os<<CRLF;
} // i 

MM_SZ_LOOP(i,la.m_sample_order,sosz2)
{
D x0=x[i];
D y0=y[i];
ow(x0,y0);
// TODO the subscript is ordinal NOT m_st()
const StrTy & colr=lc.sample_color(i);
const StrTy & sname=lc.sample_name(i);
//os<<sw.htext_text((*ii).first,x1,y1,fitsz,colr); os<<CRLF;
os<<sw.htext_text(sname,x0,y0,lc.sample_legend_font(),colr); os<<CRLF;
//os<<sw.faintline_text(x0,y0,x1,y1, colr,lc.sample_ptr_size()); os<<CRLF;
} // i 
} // svg_ii_legend


template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_hm_legend(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
typedef std::vector<D> Coord;
Coord x,y;
const IdxTy nsamples=lc.samples();
const IdxTy nnodess=lc.nodes();
lc.hm_sample_legend_pos(x,y,la.m_sample_order);
IdxTy maxlen=0;
D hspace=lc.right_legend_width();
D vspace=lc.plot_height();

D bestsz=20;
D bestangle=60;
for(IdxTy i=0; i<nsamples; ++i)
{
const StrTy & sname=lc.sample_name(i);
if (sname.length()>maxlen) maxlen=sname.length(); 
}
D hfac=1.0*maxlen/hspace;
D vfac=1.0*nsamples/vspace;
D iangle= 180.0/M_PI*atan2(hfac,vfac)*.5;
if (iangle<0) iangle=0;
if (iangle>75) iangle=75;
 bestangle=iangle;

D maxvsz=1.0/(vfac*(1.0+sin(bestangle*M_PI/180.0)));
D maxhsz=1.0/(hfac*cos(bestangle*M_PI/180.0));
//D oksz=(maxhsz>maxvsz)?maxhsz:maxvsz;
D oksz=.25*(maxhsz+maxvsz);
bestsz=oksz;
//if (oksz<1) bestsz=1;
//if (oksz>20) bestsz=1.0/(hfac*cos(bestangle*M_PI/180.0));

if (bestsz<1) bestsz=1;
if (bestsz>20) bestsz=20;
// also include splitting option 
const D fsz=bestsz;
const D angle=bestangle;
MM_ERR(MMPR4(maxlen,hspace,nsamples,vspace)<<MMPR2(maxvsz,maxhsz)<<MMPR4(hfac,vfac,iangle,oksz)<<MMPR3(bestsz,fsz,angle))

for(IdxTy i=0; i<nsamples; ++i)
{
D x0=x[i];
D y0=y[i];
ow(x0,y0);
// TODO the subscript is ordinal NOT m_st()
const StrTy & colr=lc.sample_color(i);
const StrTy & sname=lc.sample_name(i);
//os<<sw.htext_text((*ii).first,x1,y1,fitsz,colr); os<<CRLF;
//os<<sw.htext_text(sname,x0,y0,lc.sample_legend_font(),colr); os<<CRLF;
os<<sw.gtext_text(sname,x0,y0,fsz,colr,StrTy("left"),angle); os<<CRLF;
//os<<sw.faintline_text(x0,y0,x1,y1, colr,lc.sample_ptr_size()); os<<CRLF;
} // i 
} // svg_hm_legend





/////////////////////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_legend(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
os<<sw.comment_text(" iia_legend" ) ; 
typedef rendering_temp_data::name_splitter  Ns;
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
//typedef std::vector<D> Coord;
//Coord x,y;
//lc.sample_legend_pos(x,y,la.m_sample_order);
Ns ns;
Ns::splitted_type splitted;
Ns::input_type snames;
MM_SZ_LOOP(j,la.m_sample_order,soszxx)
{
const StrTy & sname=lc.sample_name(j);
snames.push_back(sname);
}

ns.split(splitted, snames);

MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
//node_iterator ee=dnc.end();
//std::vector<D> xv,yv;
const D xorg=lc.sample_iia_x(j); // plot_x();
const D dxorg=lc.sample_iia_filldx(j); // plot_x();
const D yorg=lc.plot_y()+lc.plot_height(); 
D x0=xorg;
D y0=yorg;
ow(x0,y0);
// TODO the subscript is ordinal NOT m_st()
const StrTy & colr=lc.sample_color(j);
//const StrTy & sname=lc.sample_name(j);
//os<<sw.htext_text((*ii).first,x1,y1,fitsz,colr); os<<CRLF;

const Ns::input_type & slines=splitted[j];
const D fsx=lc.sample_legend_font();
D dxeff=dxorg-fsx;
if (dxeff<0) dxeff=0;
MM_SZ_LOOP(i,slines,szs)
{
const StrTy & sname=slines[i];
const D xoff=.5*fsx+dxeff/szs*(szs-1-i);
//os<<sw.vtext_text(sname,x0,y0,lc.sample_legend_font(),colr); os<<CRLF;
os<<sw.vtext_text(sname,x0+xoff,y0,fsx,colr); os<<CRLF;

}

//os<<sw.faintline_text(x0,y0,x1,y1, colr,lc.sample_ptr_size()); os<<CRLF;
} // i 





} // svg_iia_legend
/////////////////////////////////////////////////////////


///////////////////////////////////////////////////////
//svg_ii_node_labels(os,ow,sw,la,lc,dnc);
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_node_labels(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
D yposlast=0;
MM_LOOP(ii,dnc)
{
 D  mean=la.m_means[ii.i()];
//ow(x0,y1,np.x0,lc.h);
//ow(x0,y0,np.x0,0);
const bool brk=( la.m_la_break_levels[ii.i()]<3);
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
const StrTy taxongs=taxon;
const StrTy taxon2=dnc.informative_name(ii.serial(),-3);
//MM_ERR(" ASS FUDD "<<MMPR4(taxon2, ii.i(),ii.serial(),la.m_la_break_levels[ii.serial()] ))
if (brk) taxon=taxon2+StrTy(" ... ")+taxon;
//m_cat_runs
 typedef  typename La::Runs Ru;
 typedef  typename La::run R;


D x0=lc.node_label_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
D y0=lc.node_label_y(mean); // lc.plot_y()+lc.plot_height()*mean;
const D & w=lc.node_label_width();
const D & w2=5*lc.node_label_width();

if (brk)
{
D ypos=0;
IdxTy lastsz=0;

if ( lc.node_size()<10)
{
const Ru & rup=la.m_cat_runs["phylum"];
IdxTy ij=0;
MM_LOOP(jj,rup)
{
const R & r=(*jj);
if (r.includes(ii.i()))
{
if ( lastsz<5) if (yposlast==0) ypos=.25;
break;
}
lastsz=r.len();
}
} // node_width

const D ph=lc.plot_height();
if (mean>.5) 
	{ os<<sw.vtext_text(taxon2,x0,lc.plot_y()+ypos*ph,w2);}
else 
	{ os<<sw.vtextb_text(taxon2,x0,lc.plot_y()+ph*(1.0-ypos),w2);}
yposlast=ypos;
os<<CRLF;
} // brk 
 
//const bool align_top=(yi<mid);
//const D  len=taxon.length()*lc.nodesz;
//if (align_top) if (len>(lc.h-yi)) yi=.01*lc.h;
//if (!align_top) if (len>(yi)) yi=.98*lc.h;
//ow(x0,y0,np.x0+dx,yi);
//ow(x1,y1,np.x1+dx,yi);
//if (align_top) { os<<sw.vtext_text(taxon,x0,y0,lc.nodesz*(x1-x0)); }
//else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));}



D len=taxon.length()*w*.5;
if (len>lc.plot_height()*.5) { y0=lc.plot_y(); mean=0; } 
ow(x0,y0);
if (mean>.5) { os<<sw.vtextb_text(taxon,x0,y0,w);}
else { os<<sw.vtext_text(taxon,x0,y0,w);}
if (false) if (brk)
{

if (mean>.5) { os<<sw.vtext_text(taxon2,x0,lc.plot_y(),w2);}
else { os<<sw.vtextb_text(taxon2,x0,lc.plot_y()+lc.plot_height(),w2);}

} // brk


os<<CRLF;
const bool label_big_moves=true;
//if (label_big_moves&&ii.next_ok())
if (label_big_moves&&(ii.i()!=0))
{
// TODO this is not consistent with coords above doh 
//const D xorg=lc.plot_x();
const D yorg=lc.plot_y(); 
MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
// -1 ? TODO FIXME WF I think this is a map but may be a vector 
if (ii.i()<1) { MM_ONCE(" fix this crap ",) } 
const IdxTy iii=ii.i();
D vi=(iii>0)?la.m_int[sample][ii.i()-1]:0; 
D vf=la.m_int[sample][ii.i()]; 
const bool nostraddle=((vf-mean)*(vi-mean)>.01);
if (nostraddle&&((vf-vi)>.05)){
// TODO make this consistent 
D x0=lc.node_label_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
D y=.5*(vi+vf)*lc.plot_height()+yorg; 
D x=x0; // xorg+lc.node_size()*(ii.i()+1);
ow(x,y);

//os<<sw.vtext_text(r.m_name,xpos,y0,5,color, "middle"); os<<CRLF;
	{ os<<sw.vtext_text(taxongs,x,y,1,"white","middle");}
break; // only label one 
} // if 
} // j 

} // big moves

} // ii 

} // svg_ii_node_labels


////////////////////////////////////////////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_hm_node_labels(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
// TODO FIXME what the fudd is wrong here? 
const IdxTy node_kluge=0;

const typename Lc::SORT_KEY sk=lc.sort_key();
const bool sortgroups= (sk== Lc::SORT_ALPHA_TREE) // :{MM_ERR(" alpha sort ")  dnc.sort_alpha_lineage();  break; } 
||(sk== Lc::SORT_ALPHA_GS) // :{ MM_ERR(" alpha gs sort ") dnc.sort_alpha_gs();  break; } 
;
const IdxTy nsamples=lc.samples();
const IdxTy nnodes=lc.nodes();
const bool dobrks=sortgroups;
const bool label_big_moves=false; // true;
IdxTy nmany=0;
if (nsamples>4) if (nnodes>300) if (lc.node_size()<5) nmany=(nsamples>>1);
if (nmany>4) nmany=4;
const bool label_mean=(nmany==0);
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
D yposlast=0;
MM_LOOP(ii,dnc)
{
 D  mean=la.m_means[ii.i()];
//ow(x0,y1,np.x0,lc.h);
//ow(x0,y0,np.x0,0);
const bool brk=dobrks&&( la.m_la_break_levels[ii.i()]<3);
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
const StrTy taxongs=taxon;
const StrTy taxon2=dnc.informative_name(ii.serial(),-3);
//MM_ERR(" ASS FUDD "<<MMPR4(taxon2, ii.i(),ii.serial(),la.m_la_break_levels[ii.serial()] ))
if (brk) taxon=taxon2+StrTy(" ... ")+taxon;
//m_cat_runs
 typedef  typename La::Runs Ru;
 typedef  typename La::run R;


D x0=lc.node_label_x(ii.i()+node_kluge); // ii.i()*lc.node_size()+lc.plot_x();
D y0=lc.node_label_y(mean); // lc.plot_y()+lc.plot_height()*mean;
const D & w=lc.node_label_width();
const D & w2=5*lc.node_label_width();

if (brk)
{
D ypos=0;
IdxTy lastsz=0;

if ( lc.node_size()<10)
{
const Ru & rup=la.m_cat_runs["phylum"];
IdxTy ij=0;
MM_LOOP(jj,rup)
{
const R & r=(*jj);
if (r.includes(ii.i()))
{
if ( lastsz<5) if (yposlast==0) ypos=.25;
break;
}
lastsz=r.len();
}
} // node_width

const D ph=lc.plot_height();
if (mean>.5) 
	{ os<<sw.vtext_text(taxon2,x0,lc.plot_y()+ypos*ph,w2);}
else 
	{ os<<sw.vtextb_text(taxon2,x0,lc.plot_y()+ph*(1.0-ypos),w2);}
yposlast=ypos;
os<<CRLF;
} // brk 
 
//const bool align_top=(yi<mid);
//const D  len=taxon.length()*lc.nodesz;
//if (align_top) if (len>(lc.h-yi)) yi=.01*lc.h;
//if (!align_top) if (len>(yi)) yi=.98*lc.h;
//ow(x0,y0,np.x0+dx,yi);
//ow(x1,y1,np.x1+dx,yi);
//if (align_top) { os<<sw.vtext_text(taxon,x0,y0,lc.nodesz*(x1-x0)); }
//else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));}



D len=taxon.length()*w*.5;
if (label_mean ){
if (len>lc.plot_height()*.5) { y0=lc.plot_y(); mean=0; } 
ow(x0,y0);
if (mean>.5) { os<<sw.vtextb_text(taxon,x0,y0,w);}
else { os<<sw.vtext_text(taxon,x0,y0,w);}
} // mean 
//const IdxTy nmany=5;
for(IdxTy imany=0; imany<nmany; ++imany)
{
const D mpos=1.0*imany/nmany;
D x0a=lc.node_label_x(ii.i()+node_kluge); // ii.i()*lc.node_size()+lc.plot_x();
D y0a=lc.node_label_y(mpos); // lc.plot_y()+lc.plot_height()*mean;
ow(x0a,y0a);
{ os<<sw.vtext_text(taxon,x0a,y0a,w);}
}

if (false) if (brk)
{

if (mean>.5) { os<<sw.vtext_text(taxon2,x0,lc.plot_y(),w2);}
else { os<<sw.vtextb_text(taxon2,x0,lc.plot_y()+lc.plot_height(),w2);}

} // brk




os<<CRLF;
//if (label_big_moves&&ii.next_ok())
if (label_big_moves&&(ii.i()!=0))
{
// TODO this is not consistent with coords above doh 
//const D xorg=lc.plot_x();
const D yorg=lc.plot_y(); 
MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
// -1 ? TODO FIXME WF I think this is a map but may be a vector 
if (ii.i()<1) { MM_ONCE(" fix this crap ",) } 
const IdxTy iii=ii.i();
D vi=(iii>0)?la.m_int[sample][ii.i()-1]:0; 
D vf=la.m_int[sample][ii.i()]; 
const bool nostraddle=((vf-mean)*(vi-mean)>.01);
if (nostraddle&&((vf-vi)>.05)){
// TODO make this consistent 
D x0=lc.node_label_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
D y=.5*(vi+vf)*lc.plot_height()+yorg; 
D x=x0; // xorg+lc.node_size()*(ii.i()+1);
ow(x,y);

//os<<sw.vtext_text(r.m_name,xpos,y0,5,color, "middle"); os<<CRLF;
	{ os<<sw.vtext_text(taxongs,x,y,1,"white","middle");}
break; // only label one 
} // if 
} // j 

} // big moves

} // ii 

} // svg_hm_node_labels










/////////////////////////////////////////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_node_labels(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

os<<sw.comment_text(" iia_node_labels" ) ; 
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
const D yorg=lc.plot_y(); 
D yposlast=0;
D lastvi=0;
D lastvf=0;
MM_LOOP(ii,dnc)
{
const IdxTy sample0=la.m_sample_order[0];
const IdxTy last=la.m_sample_order.size()-1;
const IdxTy samplef=la.m_sample_order[last];
D vi=la.m_int[sample0][ii.i()]; 
D vf=la.m_int[samplef][ii.i()]; 
const D xorg=lc.sample_iia_x(last); // plot_x();
const D dxorg=lc.sample_iia_filldx(last); // plot_x();
//const D yorg=lc.plot_y(); 
D xf=xorg+dxorg;
D dyi= -(vi-lastvi)*.5;
D dyf= -(vf-lastvf)*.5;
D yi=(vi+dyi)*lc.plot_height()+yorg; 
D yf=(vf+dyf)*lc.plot_height()+yorg; 
lastvi=vi;
lastvf=vf;
/// D  mean=la.m_means[ii.i()];
//ow(x0,y1,np.x0,lc.h);
//ow(x0,y0,np.x0,0);
const bool brk=( la.m_la_break_levels[ii.i()]<3);
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
const StrTy taxongs=taxon;
const StrTy taxon2=dnc.informative_name(ii.serial(),-3);
//MM_ERR(" ASS FUDD "<<MMPR4(taxon2, ii.i(),ii.serial(),la.m_la_break_levels[ii.serial()] ))
if (brk) taxon=taxon2+StrTy(" ... ")+taxon;
//m_cat_runs
 //typedef  typename La::Runs Ru;
 //typedef  typename La::run R;
D x0=lc.plot_width()*.5; // lc.node_label_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
const D fsz=1;
D y0=yposlast+fsz; // lc.node_label_y(mean); // lc.plot_y()+lc.plot_height()*mean;
//const D & w=lc.node_label_width();
//const D & w2=5*lc.node_label_width();
if (!true) {D _x=x0; D _y=y0; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
if (true) {D _x=xf; D _y=yf; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
if (true) {D _x=0; D _y=yi; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
//else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));}
yposlast=y0;
os<<CRLF;
} // ii 
} // svg_oo_node_labels

/////////////////////////////////////////////////////////////////////////////

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_node_overlay(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

os<<sw.comment_text(" iia_node_overlays" ) ; 
//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
const D yorg=lc.plot_y(); 
//D lastvi=0;
//D lastvf=0;
MM_LOOP(ii,dnc)
{
//const IdxTy sample0=la.m_sample_order[0];
//const IdxTy last=la.m_sample_order.size()-1;
//const IdxTy samplef=la.m_sample_order[last];
//D vi=la.m_int[sample0][ii.i()]; 
//D vf=la.m_int[samplef][ii.i()]; 
//const D xorg=lc.sample_iia_x(last); // plot_x();
//const D dxorg=lc.sample_iia_filldx(last); // plot_x();
//const D yorg=lc.plot_y(); 
//D xf=xorg+dxorg;
//D dyi= -(vi-lastvi)*.5;
//D dyf= -(vf-lastvf)*.5;
//D yi=(vi+dyi)*lc.plot_height()+yorg; 
//D yf=(vf+dyf)*lc.plot_height()+yorg; 
//lastvi=vi;
//lastvf=vf;
/// D  mean=la.m_means[ii.i()];
//ow(x0,y1,np.x0,lc.h);
//ow(x0,y0,np.x0,0);
const bool brk=( la.m_la_break_levels[ii.i()]<3);
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
const StrTy taxongs=taxon;
const StrTy taxon2=dnc.informative_name(ii.serial(),-3);
//MM_ERR(" ASS FUDD "<<MMPR4(taxon2, ii.i(),ii.serial(),la.m_la_break_levels[ii.serial()] ))
if (brk) taxon=taxon2+StrTy(" ... ")+taxon;
MM_SZ_LOOP(j,la.m_sample_order,sosz)
{
const IdxTy sample=la.m_sample_order[j];
const auto & vm=(*ii).values();
const auto & dvi=vm.find(sample);
D val=0;
if (dvi!=vm.end()) val=(*dvi).second;
else continue;
if (val==0) continue;
const D xorg=lc.sample_iia_x(j); // plot_x();
const D dxorg=lc.sample_iia_filldx(j); // plot_x();
const D yorg=lc.plot_y(); 
// -1 ? TODO FIXME WF I think this is a map but may be a vector 
if (ii.i()<1) { MM_ONCE(" fix this crap ",) } 
const IdxTy iii=ii.i();
D vi=(iii>0)?la.m_int[sample][ii.i()-1]:0; 
D vf=la.m_int[sample][ii.i()]; 
D dyi=(vf-vi)*.75;
D y0=(vi+dyi)*lc.plot_height()+yorg; 
D ybase=(vf)*lc.plot_height()+yorg; 
if ( ybase<(y0+1)) continue;
//m_cat_runs
 //typedef  typename La::Runs Ru;
 //typedef  typename La::run R;
D x0=xorg+.2*dxorg; // lc.plot_width()*.5; // lc.node_label_x(ii.i()); // ii.i()*lc.node_size()+lc.plot_x();
const D fsz=1;
// +1 is for contrast right now no specific function call 
StrTy concol=lc.node_iia_contrast_color(ii.i());
//D y0=yposlast+fsz; // lc.node_label_y(mean); // lc.plot_y()+lc.plot_height()*mean;
//const D & w=lc.node_label_width();
//const D & w2=5*lc.node_label_width();
//if (!true) {D _x=x0; D _y=y0; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
//if (!true) {D _x=xf; D _y=yf; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
//if (!true) {D _x=0; D _y=yi; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz); }
if (true) 
	{D _x=x0; D _y=y0; ow(_x,_y);  os<<sw.htext_text(taxon,_x,_y,fsz,concol); os<<CRLF;}
//else { os<<sw.vtextb_text(taxon,x0,y0,lc.nodesz*(x1-x0));}
//yposlast=y0;
} // j 

//os<<CRLF;
} // ii 
} // svg_oo_node_labels

/////////////////////////////////////////////////////////////////////////////









/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//svg_ii_node_outlines(os,ow,sw,la,lc);
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_node_outlines(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

//typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
MM_LOOP(ii,dnc)
{
//ow(x0,y1,np.x0,lc.h);
//ow(x0,y0,np.x0,0);
D x0=lc.plot_x()+lc.node_size()*ii.i();
D x1=lc.plot_x()+lc.node_size()*(ii.i()+1);
//D x1=x0; // lc.plot_x()+lc.node_size()*(ii.i()+1);
D y0=lc.plot_y();
D y1=lc.plot_y()+lc.plot_height();
ow(x0,y0);
ow(x1,y1);
const IdxTy brlevel=( la.m_la_break_levels[ii.i()]) ;


//if ( la.m_la_break_levels[ii.i()]<2) 
if ( brlevel<3) 
{
os<<sw.line_text(x0,y0,x0,y1, StrTy("blue"),.05*lc.node_size()); os<<CRLF;
}
else
{

if ( brlevel>4) 
{os<<sw.line_text(x0,y0,x0,y1, StrTy("blue"),.02*lc.node_size()); os<<CRLF;}
else {os<<sw.line_text(x0,y0,x0,y1, StrTy("red"),.02*lc.node_size()); os<<CRLF;}


}
//const bool highlight_column =(np.x0>30)&&(np.x0<50);
const bool highlight_columns =true; // (np.x0>30)&&(np.x0<50);
if (highlight_columns)
{
const bool flag_node=(*ii).flagged();

if (flag_node){ 
StrTy taxon=dnc.informative_name(ii.serial(),lc.name_levels());
for (IdxTy i=0; i<(*ii).nflags(); ++i)
{
const IdxTy & flag= (*ii).get_flag(i);
MM_ERR(" flagged node "<<MMPR4(flag,ii.i(),ii.serial(),taxon)<<MMPR2(lc.flag_color(flag), lc.flag_name(flag))) 
// x,y,w,h doh !
os<<sw.shade_rect_text(x0,y0,x1-x0,y1-y0, lc.flag_color(flag),.5); os<<CRLF;
} // i 
} // flag_node 

} // highlight_columns


} // ii 

} // node_outlines

template <class Os, class Ow, class Sw, class La, class Lc>
void svg_ii_margins(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc)
{
D x1=lc.plot_x(); //,ys;
D y1=lc.plot_y(); //,ys;
D x2=lc.plot_width()+x1; //,ys;
D y2=lc.plot_height()+y1; //,ys;
ow(x1,y1); 
ow(x2,y2); 

//os<<sw.hstripes(x1,y1,x2,y2,10,"black",lc.visthick);
os<<sw.hstripes(x1,x2,y1,y2,10,lc.stripe_color(),lc.stripe_thickness());
os<<CRLF;
//draw_sample_legend(os,ow,ls.sw, ls.legend_map,lc,  node_count,  lc.h);
///os<<ss.str();
x1=lc.sign_x(); // 3.0*D(lc.nodew*node_count/4);
y1=lc.sign_y(); // (1+.9*lc.ymargin)*lc.h-lc.taxasz;
ow(x1,y1);
os<<sw.htext_text(lc.sign(),x1,y1,lc.sign_height(),lc.sign_color());
os<<CRLF;

x1=lc.ylab_x(); // 3.0*D(lc.nodew*node_count/4);
y1=lc.ylab_y(); // (1+.9*lc.ymargin)*lc.h-lc.taxasz;
ow(x1,y1);
//os<<sw.vtextb_text(lc.ylab(),x1,y1,lc.ylab_size()); 
os<<sw.vtext_text(lc.ylab(),x1,y1,lc.ylab_size()); 
os<<CRLF;

for(IdxTy i=0; i<=10; ++i)
{
x1=.5*lc.plot_x(); // 3.0*D(lc.nodew*node_count/4);
y1=lc.plot_y()+lc.plot_height()*(D(i)/10.0); // (1+.9*lc.ymargin)*lc.h-lc.taxasz;
ow(x1,y1);
Ss ss; ss<<(D(i)/10.0); 
//os<<sw.vtextb_text(lc.ylab(),x1,y1,lc.ylab_size()); 
os<<sw.vtext_text(ss.str(),x1,y1,.02*lc.plot_height()); 
os<<CRLF;





}



}

typedef std::map<IdxTy, StrTy> SuggOrd;
SuggOrd m_suggested_order;
// this is a map from order to string name 
void initial_sample_order(const SuggOrd & m)
{
m_suggested_order=m;

}

typedef std::map<StrTy,StrTy> param_hash_t;
void write_svg_ii(std::ostream & os, const IdxTy flags, mjm_pheno_notes & pn
, param_hash_t & ph )
{
const StrTy proname=m_flp.pheno_notes_name(); // "z";
const StrTy field=m_flp.pheno_notes_field(); // "aerob";
mjm_pheno_notes::hbpheno_oracle oracle(pn,proname,field);
write_svg_ii(os,flags,oracle,ph);

}

void write_svg_ii(std::ostream & os, const IdxTy flags
, param_hash_t & ph )
{
mjm_pheno_notes::null_oracle oracle;
write_svg_ii(os,flags,oracle,ph);
}
template <class To > 
void write_svg_ii(std::ostream & os, const IdxTy flags, To & oracle 
, param_hash_t & ph )
{
typedef rendering_temp_data::data_node_collection Nodes;
//typedef Nodes::iterator node_iterator;
typedef Nodes::analysis La;
typedef Nodes::consts Lc;
Lc lc(ph);
std::vector<StrTy> vsugg;
MM_LOOP(ii,m_suggested_order) { vsugg.push_back((*ii).second); } 
lc.m_suggested_initial_sample_order=vsugg;
const IdxTy style=flags;
const bool ignore_notes=false;
//typedef rendering_temp_data::data_node Node;
Nodes dnc;
NodeValues& nv=  m_values;
TaxTree & tt=m_tax_tree;
MM_ERR(" checking nonterminal nodes")
check_nonterminal_nodes_ii( nv, tt,lc );
MM_ERR(" collectin node data")
if (ignore_notes) {
mjm_pheno_notes::null_oracle noracle;
collect_nodes(dnc,tt,nv,noracle,lc);

//MM_ERR(" danger will robinson saorting by gs")
//if (false) dnc.sort_alpha_lineage(); 
//dnc.sort_alpha_gs(); 
}
else { 
//MM_ERR(" checking nonterminal nodes")
//check_nonterminal_nodes_ii( nv, tt );
collect_nodes(dnc,tt,nv,oracle,lc);
//dnc.sort_alpha_notes(); 
//MM_ONCE(" sort order fudded up ",)
//MM_ERR(" danger will robinson saorting by gs")
//if ( false) dnc.sort_alpha_lineage(); 

//dnc.sort_alpha_gs(); 
}

const Lc::SORT_KEY sort_keys=lc.sort_key();
MM_ERR(" sorting  node data "<<MMPR(sort_keys))
switch (sort_keys)
{
case Lc::SORT_ALPHA_TREE :{MM_ERR(" alpha sort ")  dnc.sort_alpha_lineage();  break; } 
case Lc::SORT_ALPHA_GS :{ MM_ERR(" alpha gs sort ") dnc.sort_alpha_gs();  break; } 
case Lc::SORT_NOTES  :{ MM_ERR(" notes sort ") dnc.sort_alpha_notes();  break; }
case Lc::SORT_VALUES  :{ 
// these do not exist until the analysis is done... 
// the nodes have been collected by dnc but that is all 
MM_ERR(" thes amepls look up is fudded fudd shot this will bomb shot ")
StrTy sname=""; // lc.sample_name(1);
StrTy smod=""; // lc.sample_name(1);
//MM_LOOP(ff,nv){ MM_LOOP(gg,(*ff).second) {sname=(*gg).first; break; }}
sname=ph["sort-sample"];
smod=ph["sort-modifiers"];
MM_LOOP(ii,ph) { MM_ERR(" plot-hash "<<MMPR2((*ii).first,(*ii).second)) }
//MM_ERR(" values sort "<<sname ) dnc.sort_values(sname);  break; }
//MM_ERR(" values sort "<<sname ) dnc.sort_values(sname,nv);  break; }
MM_ERR(" values sort "<<sname ) dnc.sort_values(sname,smod);  break; }
default: MM_ERR(" bad sort key "<<MMPR(sort_keys)) 
} // switch 

// now this needs a list of all sample names that may not be on each node 
MM_ERR(" analyzing layout")
La la;
dnc.analyze(la,lc);
MM_ERR(" starting writing svg")
//analyze_layout(la,dnc);
mjm_svg_writer sw;
mjm_output_warp ow;
//Lc lc(la);
lc.do_analysis(la);
switch (style)
{
case 0: { svg_ii_do_write(os,ow,sw,la,lc,dnc); break; } 
case 1: { svg_iia_do_write(os,ow,sw,la,lc,dnc); break; } 
case 2: {//lc.hm_sample_sort(true); 
//m_suggested_order=m;
if (m_suggested_order.size()!=0) lc.hm_sample_sort(m_suggested_order,la,dnc.st()); 
// #error fudding la sample_order fuc 
// the la sample order needs to match FUDD 
 svg_hm_do_write(os,ow,sw,la,lc,dnc); break; } 

} // style 

} // write_svg_ii

template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_ii_do_write(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{
svg_ii_start(os,ow,sw,la,lc);
svg_ii_flag_legend(os,ow,sw,la,lc,dnc);
// draw background first including stripes and node flagging 
// (os,sw,la,dnc);
svg_ii_cats(os,ow,sw,la,lc,dnc);
svg_ii_margins(os,ow,sw,la,lc);
svg_ii_legend(os,ow,sw,la,lc,dnc);
svg_ii_node_outlines(os,ow,sw,la,lc,dnc);
svg_ii_node_labels(os,ow,sw,la,lc,dnc);
// add legend lines so the curves need to be calculated first
// add the curves
svg_ii_curves(os,ow,sw,la,lc,dnc);
// finally the node names which should not be obscured
MM_ERR(" writing end text ")
os<<sw.end_text();
os<<CRLF;
} // 


template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_iia_do_write(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

svg_iia_start(os,ow,sw,la,lc);
svg_iia_flag_legend(os,ow,sw,la,lc,dnc);
// draw background first including stripes and node flagging 
// (os,sw,la,dnc);
//svg_ii_cats(os,ow,sw,la,lc,dnc);
svg_ii_margins(os,ow,sw,la,lc);
svg_iia_legend(os,ow,sw,la,lc,dnc);
//svg_ii_node_outlines(os,ow,sw,la,lc,dnc);
// add legend lines so the curves need to be calculated first

// add the curves
svg_iia_rect(os,ow,sw,la,lc,dnc);
svg_iia_lines(os,ow,sw,la,lc,dnc);

svg_iia_node_labels(os,ow,sw,la,lc,dnc);
svg_iia_node_overlay(os,ow,sw,la,lc,dnc);

// finally the node names which should not be obscuresnamed
MM_ERR(" writing end text ")
os<<sw.end_text();
os<<CRLF;
} // 
template <class Os, class Ow, class Sw, class La, class Lc,class Dnc>
void svg_hm_do_write(Os&os,Ow&ow,Sw&sw,La&la,Lc&lc, Dnc& dnc)
{

svg_hm_start(os,ow,sw,la,lc);

// add the curves first as these can be overwritten 
svg_hm_rect(os,ow,sw,la,lc,dnc);

///svg_iia_flag_legend(os,ow,sw,la,lc,dnc);
// draw background first including stripes and node flagging 
// (os,sw,la,dnc);
svg_ii_cats(os,ow,sw,la,lc,dnc);
///svg_ii_margins(os,ow,sw,la,lc);
///svg_iia_legend(os,ow,sw,la,lc,dnc);
svg_hm_legend(os,ow,sw,la,lc,dnc);
//svg_ii_node_outlines(os,ow,sw,la,lc,dnc);
// add legend lines so the curves need to be calculated first
///svg_iia_lines(os,ow,sw,la,lc,dnc);

svg_hm_node_labels(os,ow,sw,la,lc,dnc);
//svg_iia_node_overlay(os,ow,sw,la,lc,dnc);

// finally the node names which should not be obscured
MM_ERR(" writing end text ")
os<<sw.end_text();
os<<CRLF;
} // 


void write_svg(std::ostream & os, const IdxTy flags)
{

NodeValues& nv=  m_values;

typedef std::vector<IdxTy> Tv;
Tv tv,tvold;
//mjm_output_warp ow(1.0,1.0,1,1);
mjm_output_warp ow;
//mjm_circular_warp ow(xpitch,ypitch,w,h);
IdxTy node_count=0; 
TaxTree & tt=m_tax_tree;
tt.sort_for_ui();

IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
LayoutState ls;
//IdxTy  global_stats(Tsv & sv, Tnl & nl, Tnv & nv,  TaxTree & tt, const IdxTy&  _node)
//global_stats(node_count,ls.sv,ls.nl,nv,  tt, node,taxo);
global_stats(node_count,ls,nv,  tt, node,taxo);

LayoutConst lc(node_count);
//IdxTy itaxa=0;
fix_layout(ls,lc);

D xs=0; // 
D ys=0; //
//ow(xs,ys,1.05*lc.w,(1.0+lc.ymargin)*lc.h);
ow(xs,ys,lc.tw,lc.th);
os<<ls.sw.start_text(" test foo",xs,ys);
os<<ls.sw.frame_text("#00808080",xs,ys);
os<<ls.sw.stroke_text("red",1000);
os<<CRLF;
Ss ss;

node=best_root; // tt.roots()[ridx];
TaxTree::tree_terminal_iterator ti(tt,node);
IdxTy nidx=0;
//MM_LOOP(ii,nv)
IdxTy snamed=0;
while (ti.valid())
{ 
node=ti.node();
const npos & np=ls.nl[node];
const bool right=((np.x0)>(lc.nodew*(node_count>>1)));
//MM_ERR(MMPR2(node,tt.node_name(node)))
auto ii=nv.find((node));
IdxTy sami=0;
MM_LOOP(jj,(*ii).second)
{
if (false) if (sami>20)
{ { MM_ONCE(" breaking due to too many samples ",) break; } }
const StrTy & sn=(*jj).first;

wpos & cp=get_cursor(sn,snamed,ls.sc,ls.legend_map,node,nidx,node_count);
plateau_and_slope(os,ow,np,sn,(*jj).second,lc,ls,cp);
++sami;
} // jj
StrTy t2="";
// caption 
 StrTy taxon=taxon_splits(t2,tv,tvold,np,ow,os,ls.sw,tt,node,taxo,lc.h,lc );
// plateaus still overwrite 
 //label_node( ss,ls.sw,ow,lc,right,np,taxon,ls);
 label_node( ss,ow,lc,right,np,taxon,ls,nidx);

if (t2.c_str()[0]!=0) 
	//{ mark_taxa_change(t2,ls.itaxa,os,ls.sw,ow,np,lc.ytaxa0,lc.ytaxa1,lc.taxasz,lc.h,right); }
	//{ mark_taxa_change(t2,ls.itaxa,os,ls.sw,ow,np,lc,right); }
	{ mark_taxa_change(t2,ls,os,ow,np,lc,right); }

tvold=tv;
ti.inc();
++nidx;
} // ti 
os<<ls.sw.hstripes(0,lc.w,0,lc.h,10,"black",lc.visthick);
draw_sample_legend(os,ow,ls.sw, ls.legend_map,lc,  node_count,  lc.h);
os<<ss.str();
ss.str()="";
{
D x1,y1;
ow(x1,y1,3.0*D(lc.nodew*node_count/4),(1+.9*lc.ymargin)*lc.h-lc.taxasz);
os<<ls.sw.htext_text(StrTy( "Marchywka   test code"),x1,y1,lc.taxasz,"black");
os<<CRLF;
}


os<<ls.sw.end_text();
os<<CRLF;
}


template <class Tos, class Tw, class Ts, class Tm,class Lc  > 
void draw_sample_legend(Tos & os, Tw & ow, Ts & sw,Tm & legend_map, Lc & lc, const IdxTy & node_count, const D & h)
{
const D uymargin=.1;
const IdxTy n=legend_map.size();
IdxTy maxsz=0;
MM_LOOP(ii,legend_map) { const IdxTy len=(*ii).first.length(); if (len>maxsz) maxsz=len; } 
const D A=lc.nodew*node_count*uymargin*h;
const D cc=n*maxsz;
if (cc==0) return; 
const D fitsz=.5*sqrt(A/cc);
const D dy=1.1*fitsz;
D fudd1=lc.nodew*node_count;
if (fudd1>maxsz) fudd1-=maxsz;
IdxTy div=IdxTy(1*0+(fudd1/(maxsz+2)/fitsz));
if (div==0) div=1;
MM_ERR(" levend map "<<MMPR4(maxsz,fitsz,dy,div))
IdxTy legs=0;
typename Tm::sorted x=legend_map.sort();
//MM_LOOP(ii,legend_map)
const IdxTy nt=1+(n/div);
MM_LOOP(ii,x)
{
const IdxTy legst=legs;

//const D ykey=1.02*h+dy*(legst/div);
//const D xkey=1.0*(legst%div)*(node_count)/div;

D ykey=1.02*h+dy*(legs%nt);
D xkey=1.0*(IdxTy((legs/nt)))*(1.0*lc.nodew*node_count)/D(div);
if (xkey>(lc.nodew*(node_count>>1)))
{
xkey=1.0*(IdxTy(div+(div>>1)-(legs/nt)))*(1.0*node_count*lc.nodew)/D(div);
ykey=1.02*h+dy*((legs+nt-1)%nt);

}


//MM_ERR(MMPR4(legs,nt,div,ykey)<<MMPR2(xkey,node_count))

const auto & cp=(*ii); // (*ii).second;
D x0,y0,x1,y1;
ow(x1,y1,2.0+xkey,ykey+.5*fitsz);
ow(x0,y0,cp.x,cp.y);
const StrTy & colr=cp.color;
//os<<sw.htext_text((*ii).first,x1,y1,fitsz,colr); os<<CRLF;
os<<sw.htext_text((*ii).nm,x1,y1,fitsz,colr); os<<CRLF;
os<<sw.faintline_text(x0,y0,x1,y1, colr,.001*node_count*lc.nodew); os<<CRLF;
++legs;
}


}



/////////////////////////////////////////////////////////////////////
private:
void Init()
{
//delete m_dmel;
//m_dmel= new Dmel();
m_done=false;
//load_literals();
//load_shares();
//load_handlers();
}

/*
void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<" "<<code<<CRLF;
    ss.str(StrTy(""));
}

*/

void DMel(const StrTy & e)
{
MM_ERR(e)
if (m_dmel!=0) {m_dmel->event("wtf",e); }
}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
	if (m_dmel!=0) {m_dmel->event("wtf",ss.str()); }
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; }
	if (m_dmel!=0) {m_dmel->event(code,ss.str()); }
    ss.str(StrTy(""));
}
void DumpDMel() // (OsTy & os)
{
if (m_dmel!=0) { MM_ERR(MMPR((m_dmel->string(1)))) } 
else { MM_ERR(" m_dmel is null ") }

}

CmdMap m_cmd_map;
CompMap m_comp_map;


bool m_done;

Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
CharMat m_char_mat;
MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 

NodeValues m_values;


CounterMap m_cm;
CliTy m_cli;

}; //mjm_tree_viz



/////////////////////////////////////////////////////////

#ifdef  TEST_STRING_SEQ__
int main(int argc,char **args)
{
typedef mjm_string_seq  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

