#ifndef MJM_LINC_GRAPH_H__
#define MJM_LINC_GRAPH_H__
 
#include "mjm_globals.h"

#include "mjm_config_resolve.h"
#include "mjm_string_picker.h"
#include "mjm_rule_list.h"
#include "mjm_coefficients.h"
#include "mjm_file_name_gen.h"
#include "mjm_rag_splitter.h"
// units
#include "mjm_unit_crap.h"

// #define USE_RATIONALS
#ifdef USE_RATIONALS
#define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
#include "mjm_bigint_mmpr.h"
//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
// add for the multiplicity tests
#include "mjm_integrals.h"
//#include "../copied_code/gmp_factorize.h"
#include "mjm_rational.h"
#endif // USE_RATIONALS

// for the presence absence vector 
//#include "mjm_ordering.h"
#include "mjm_char_mat.h"
#include "mjm_part_iterators.h"
#include "mjm_iterator_base.h"

#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
//#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_svg_writer.h"

#include "mjm_canned_methods.h"
#include "mjm_dog_plots.h"

//#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
#include "mjm_fasta_ii.h"

//#include "mjm_collections.h"
//#include "mjm_svg_writer.h"

#ifdef USE_RATIONALS
#include "mjm_rational_things.h"
#include "mjm_table_managers.h"
#endif

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


2018-08-02 copy from mjm_linc_graph
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 
/*

*/


////////////////////////////////////////////////////////////////

class linc_graph_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
linc_graph_params( const StrTy & nm) : Super(nm) {}
linc_graph_params() : Super() {}
// the defaults are only loaded into the map on default instance.
// but still supplied in live instance. 
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
//IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
StrTy datemin() const { return m_map.get_string("datemin",""); }
StrTy datemax() const { return m_map.get_string("datemax",""); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//IdxTy maxdepth() const { return m_map.get_uint("maxdepth",3); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }
//bool print_haves() const { return m_map.get_bool("print_haves",!true); }
//bool print_havenots() const { return m_map.get_bool("print_havenots",!true); }
//bool print_if_have() const { return m_map.get_bool("print_if_have",!true); }
//bool suppress_vector() const { return m_map.get_bool("suppress_vector",!true); }
//bool add_level() const { return m_map.get_bool("add_level",true); }
//bool print_hit() const { return m_map.get_bool("print_hit",true); }
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  

void defaults()
{
set("log_commands","0");
set("exit_on_err","0");
set("exit_on_err","0");
set("ragged_params","."); 
set("maxcnt",100); 
//StrTy datemin() const { return m_map.get_string("datemin",""); }
//StrTy datemax() const { return m_map.get_string("datemax",""); }

}

StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
//ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
//ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;
ss<<"datemin"<<datemin()<<sep;
ss<<"datemax"<<datemax()<<sep;

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace linc_graph_traits
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
//typedef std::vector<IdxTy > Locations;
//typedef mjm_string_ordering Ordering;

//typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


#if 0 

#endif

}; // trees_and_tables_traits
///////////////////////////////////////////////////////////////


#if 1 
class mjm_digraph
{
typedef  linc_graph_traits::Tr  Tr;
//typedef  linc_graph_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_digraph Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

typedef Tr::MyBlock  MyBlock;
typedef linc_graph_params ParamGlob;

typedef string_tokenizer St;
typedef std::map<StrTy,StrTy> SeqMap;
class node_entry
{
public:
typedef node_entry Myt;
node_entry(): m_v(0) {}
node_entry(const D & x): m_v(x) {}
operator double () const   { return m_v; } 
Myt & operator+=( const D & x) { m_v+=x; return *this; } 
D m_v;
}; // node_entry

typedef IdxTy Nono;
typedef node_entry Novo;
typedef std::map<Nono,Novo> Dmap;
typedef std::map<Nono,Dmap> Gr;

typedef std::map<Nono,IdxTy> Hm;
typedef std::map<Nono,D> Score;

//typedef std::map<IdxTy,IdxTy> TermMap;
//typedef std::map<IdxTy,TermMap> SpeciesMap;
//typedef std::map<IdxTy,SpeciesMap> GenusMap;

static int myatoi(const StrTy & s )  { return myatoi(s.c_str()); } 
static int myatoi(const char * c)  { return ::strtol(c,0,0); }

public:
mjm_digraph(): m_read_only(false)
,m_min_words(bad()),m_max_words(0),m_fields(0)

 {}
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
const IdxTy  size() const { return m_gr.size(); }

void load(const StrTy & fn )
{
std::ifstream  ifs(fn);  
load(ifs); 
}
void load(std::istream & is )
{
 IdxTy lines=0;
CommandInterpretter li(&is);
Setup(li,m_flp);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if (sz<3) 
		{
			MM_ERR(" bad line "<<li.line())
			continue;
		}
		const Words & w=li.words();
		const D score=atof(w[0].c_str());
		const IdxTy dest=m_st(w[1]);
		const IdxTy src=m_st(w[2]);
		m_gr[src][dest]+=score;
	//	m_table.add(li.words());
		++lines;
    } // nextok()
//MM_ERR(" done loading "<<MMPR3(m_min_words,m_max_words,m_fields)<< MMPR3(lines,size(),m_order.dump(0)))
//MM_ERR(MMPR(lines)<< stats())
MM_ERR(MMPR(lines))
}
void subs(const IdxTy flags) 
{
const bool dump_walk=((flags&(1<<0))==0);
const bool max2=((flags&(1<<1))!=0);
//IdxTy ml=0;
//IdxTy lev_lim=10;
//const IdxTy max_level=ml;
MM_LOOP(ii,m_gr)
{
Hm  hits;
Hm  level;
Score score;
const IdxTy & i=(*ii).first;
const StrTy & nm=m_st(i);
++hits[i];
const Dmap & dm= (*ii).second;
walk(m_gr,dm,hits,level,1,score,1);
if (dump_walk) { 
IdxTy max=0;
if (max2)
{
MM_LOOP(jj, hits) {
	IdxTy j=(*jj).first; 
	const IdxTy & lvl=level[j];
if ( lvl>max) max=lvl;
//if (max>lev_lim) break; 	
 } // j 
//if (max>lev_lim) continue; 	
} // max2

MM_LOOP(jj, hits) 
{	
	IdxTy j=(*jj).first; 
	const StrTy & nmj=m_st(j); 
	const IdxTy & lvl=level[j];
	const D  & s=score[j];
//	if ((max_level==0)||(lvl<max_level)) 
	{ MM_ERR(MMPR4(nm, lvl,nmj,s)) }

} 
}

} // ii 


} // subs

void walk(const Gr & gr, const Dmap & dm, Hm & hits, Hm & level, const IdxTy d, Score & s,  const D & scorein )
{

MM_LOOP(ii, dm)
{
const IdxTy & i=(*ii).first;
//if (hits[i]!=0) continue; // return; 
// this needs to get to the end... 
//if (hits.find(i)!=hits.end()) continue; // return; 
const StrTy & nm=m_st(i);
if (false) {MM_ERR(MMPR(nm)) } 
++hits[i];
IdxTy&  dh=level[i];
const D scoreout=scorein*(*ii).second;
if (( dh==0)||( dh>d)) {  dh=d;  s[i]=scoreout;  } 
if (hits[i]>1) continue; // return; 
const auto ig= gr.find(i);
const Dmap & dmi= (*ig).second;
walk(gr,dmi,hits,level,d+1,s,scoreout);
} // ii 

} // walk




private:
void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
//li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,';');
li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,' ');
}

ParamGlob m_flp;
//SeqMap m_map;
Gr m_gr;
//Ordering m_order;
St m_st;
bool m_read_only;
IdxTy m_min_words,m_max_words,m_fields;

}; // mjm_digraph

#endif




class mjm_linc_graph 
{
typedef linc_graph_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_linc_graph Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef mjm_string_picker<Tr> Picker;

typedef mjm_rule_list<Tr> ConfigRules;
typedef mjm_coefficients<Tr> Coefficients;
typedef mjm_coefficients_map<Tr> CoefficientsMap;


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;


typedef mjm_canned_methods Canned;
typedef mjm_rag_splitter<Tr> RagSplit ;
typedef mjm_file_name_gen<Tr> Names; 

typedef linc_graph_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

//typedef mjm_clustalw_aln Aln;
//typedef std::map<StrTy, Aln> AlnMap;

typedef mjm_digraph Dig;
typedef std::map<StrTy, Dig> DigMap;


//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;

typedef mjm_config_resolve<Tr> CoRes;
typedef mjm_logic_resolve<ParamGlob,Tr> LoRes;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

typedef mjm_unit_crap<Tr> Units;
typedef Units::qty_ty qty_ty;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;
typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;

#ifdef USE_RATIONALS
typedef mjm_rational_things  MulRat;
typedef  MulRat::big_coef_type big_coef_type;
typedef mjm_table_base<big_coef_type> TableMan;
#endif



public:

static int myatoi(const StrTy & s )  { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }

public :
mjm_linc_graph():m_flp(m_lo.flp()), m_flp_def(m_lo.flp_def()), m_dmel(new Dmel()) {Init();}
mjm_linc_graph(int argc,char **_args) : 
m_flp(m_lo.flp()), m_flp_def(m_lo.flp_def()) , m_dmel(new Dmel())
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
~mjm_linc_graph()
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



CmdMap m_cmd_map;
CompMap m_comp_map;

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
/*

static void stream_edit_fasta(const StrTy & d,const StrTy & s,const IdxTy flags)
static void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv )
static void cmd_read_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
static void cmd_write_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
static void cmd_add_to_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
static void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv, FastaMap & m, RaggedMap & rm  )
static StrTy zymo_condense(const StrTy & s ) // const
static void cmd_transpose_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
static void cmd_read_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
static void cmd_write_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
static void cmd_add_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
template <class Ttm > static void cmd_tt(Cip & cip , LocalVar & lv, Ttm & ttm ) 



*/




void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_stream_edit_fasta( cip ,  lv ); }
void cmd_read_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_read_fasta(cip ,  lv,  m_fasta_map  ); }
void cmd_write_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_write_fasta(cip ,  lv,  m_fasta_map  ); }


class dog_day_stat
{
public:
// 2021-09-23 changed default to output units 
dog_day_stat(): m_dog(),m_name(), //  m_invalid(true),  
m_n(0),m_ntot(0), m_tot(0),m_max(0), m_avg(0),m_wavg(0), m_avg_dose(0)
,m_freq(0),m_weight(0)
,m_invalid(true),m_mixed_units(false),m_output_units(!false)
{

//IdxTy m_n,m_ntot;
//D m_tot, m_max, m_avg, m_avg_dose,m_freq;
//bool m_invalid;

}
const StrTy & key_name() const { return m_name; } 
void key_name(const StrTy & x)  {  m_name=x; } 
// the n here is the non-zero terms, ntot is the
// toal time period 
void calculate(const IdxTy ntot)
{
if (ntot==0) return ;
m_ntot=ntot;
m_avg=m_tot/ntot;
//if (m_weight==0) m_wavg=m_avg; else m_wavg=m_tot/m_weight;
m_wavg=m_weight;
m_freq=1.0*m_n/ntot;
m_invalid=false;

}
// this no longer works right with binning etc.
// it relies on zero entries for all days. 
void add_dose(const qty_ty & _v, const IdxTy & serial, const D & w=1 )
{
const D v=_v.qty();
const StrTy sfx=_v.sfx();
MM_ERR(" stubbed out "<<MMPR2(v,sfx))
if (m_n==0)
{
m_units=sfx;

} // n==0
else if (m_units!=sfx) { m_mixed_units=true;
MM_ERR(" mixed units "<<MMPR2(m_units,sfx))
}
++m_n;
++m_unit_check[sfx];
m_tot+=v;
m_weight+=w*v;
if (v>m_max) m_max=v;
m_invalid=!false;
if (v!=0) ++m_non_zed[serial];
++m_serials[serial];
} // add_dose


void add_dose_xx(const D & v, const IdxTy & serial, const D & w=1 )
{
++m_n;
m_tot+=v;
m_weight+=w*v;
if (v>m_max) m_max=v;
m_invalid=!false;
if (v!=0) ++m_non_zed[serial];
++m_serials[serial];
}
IdxTy fields() const { return 0; } 
StrTy to_ssv(const IdxTy flags =0,const IdxTy field=0 ) const
{

// case ? lol 
//const bool one=Bit(flags,0);
Ss ss;
/*
//ss<<std::setprecision(8);
IdxTy fn=0;
bool next_ok=(!one||(field==fn));
if (next_ok)ss<<m_avg;
++fn; bool next_ok=(!one||(field==fn));
if ( m_output_units&&next_ok) { 
if (ss.str().length()) ss<<" ";
if (m_units=="") ss<<"-"; else ss<<m_units; } 

++fn; bool next_ok=(!one||(field==fn));
if (next_ok){
if (ss.str().length()) ss<<" ";
ss<<m_max;
}
++fn; bool next_ok=(!one||(field==fn));
if (next_ok){
if (ss.str().length()) ss<<" ";
ss<<(m_non_zed.size());
}
++fn; bool next_ok=(!one||(field==fn));
if (next_ok){
if (ss.str().length()) ss<<" ";
ss<<(m_ntot);
}
*/
return ss.str();
} // to_ssv
StrTy to_string() const
{
Ss ss;
if (m_max!=0){
//	ss <<std::setw(6+1);
//if (m_max>(1e-3)) ss <<std::fixed;else  ss <<std::scientific;
ss<<std::setprecision(2);
if (m_avg>=100) ss<<std::setprecision(3);
if (m_avg>=1000) ss<<std::setprecision(4);
if (m_avg<=(1e-2)) { ss <<std::scientific<<std::setprecision(2); } 

ss<<m_avg;
if ( m_output_units) { ss<<" "; ss<<m_units; } 

//<<"(<="
ss<<";";
ss.unsetf(std::ios_base::floatfield); 
//ss <<std::scientific<<std::setprecision(2);}
ss<<std::setprecision(2);
if (m_max>=100) ss<<std::setprecision(3);
if (m_max>=1000) ss<<std::setprecision(4);
//if (m_max>(1e3))  {ss <<std::fixed<<std::setprecision(3);}
//if (m_max>(1e6))  {ss.unsetf(std::ios_base::floatfield); ss <<std::scientific<<std::setprecision(2);}

//<<m_max<<",\\%="
//<<m_max<<",\\%="
ss<<m_max<<";"
//<<std::setprecision(2)
//<<(m_freq)<<")";
//<<(m_n)<<"/"<<m_ntot<<")";

//<<(m_n)<<"/"<<m_ntot<<"";
//<<(m_non_zed.size())<<"/"<<(m_serials.size())<<"";
<<(m_non_zed.size())<<"/"<<(m_ntot)<<"";


}
return ss.str();

}
StrTy dump(const IdxTy flags=0, const StrTy&  sep=StrTy(" ")) const
{
Ss ss;
if (flags==0) ss<<MMPR2(m_n,m_ntot)<<MMPR4(m_tot,m_max,m_avg,m_avg_dose)<<MMPR3(m_freq,m_wavg,m_invalid);
if (flags==1) ss<<MMPR4S(m_n,m_ntot,m_tot,m_max)<<MMPR4S(m_avg,m_avg_dose,m_freq,m_wavg)<<MMPRS(m_invalid);
if (m_output_units) ss<<MMPR2(m_units,m_mixed_units);
//if (flags==1) ss<<MMPR4S(m_n,m_ntot,m_tot,m_max)<<MMPR4S(m_avg,m_avg_dose,m_freq,m_invalid);
return ss.str();

}

StrTy m_dog, m_name;
IdxTy m_n,m_ntot;
D m_tot, m_max, m_avg, m_wavg, m_avg_dose,m_freq,m_weight;
bool m_invalid,m_mixed_units,m_output_units;
StrTy m_units;
std::map<StrTy,IdxTy> m_unit_check;
std::map<IdxTy,IdxTy> m_serials,m_non_zed;
}; // dog_day_stat
typedef dog_day_stat Stat;



// key is dietary component  name
typedef std::map<StrTy, Stat> Dstat;
// key is dog name 
typedef std::map<StrTy, Dstat> DDstat;

class table_decorators
{
public:

void decorate_y(const StrTy & nm, StrTy & name) {
latex_blocks::add_superscripts(name,superscripts[name]);
// this is the old approach where post hoc units are in the confi file
// but now they are in the data file doh... 
// appending units should prevent these from being found... 
if (units.find(nm)!=units.end()) {name=name+"("+units[nm]+")"; }
}

void decorate_y(const StrTy & nm, StrTy & name, const StrTy & key,
const IdxTy flags=0) {
latex_blocks::add_superscripts(name,superscripts[key],flags);
// this is the old approach where post hoc units are in the confi file
// but now they are in the data file doh... 
// appending units should prevent these from being found... 
if (units.find(key)!=units.end()) {name=name+"("+units[key]+")"; }
}




template <class Ty> 
void concat(StrTy & x, Ty & line,  const IdxTy first, const IdxTy last)
{
 for (IdxTy j=first; j<last; ++j) 
{
	if ( x.length()!=0) x+=StrTy(" ");  
	x+=line[j] ;  
} 

}


const bool check_x(const StrTy & name) 
{
//if (!_check_x) return true;
//if (subdogs.find(name)==subdogs.end()) return false; 
//return  (subdogs.found(name)); //  return false; 
return  (subdogs.test(name)); //  return false; 
//return true;
}
//if (check_foods) { if (subnames.find(key)==subdogs.end()) continue; } 
const bool check_y(const StrTy & name) 
{
//if (!_check_y) return true;
//if (subnames.find(name)==subnames.end()) return false; 
return  (subnames.test(name)); //  return false; 
//return true;

}

std::map<StrTy, std::vector<StrTy > >  superscripts;
//std::map<StrTy, IdxTy> subdogs,subnames;
//std::map<StrTy, IdxTy> subnames;
Picker subdogs ,subnames,subdates;
std::map<StrTy, StrTy> units;
typedef std::map<StrTy,StrTy> TableParams;
TableParams tparam,vars,footnotes; //,varnames;
bool scope_block // ,_check_x
//,_check_y
;
StrTy comment;
 IdxTy dsz;
 IdxTy imin;
 IdxTy imax;
 IdxTy serialmin;
 IdxTy serialmax;
 StrTy datemin;
StrTy datemax;

IdxTy period,filtering;
CoefficientsMap m_coefs_map;


void decorate_latex(OsTy & os)
{
if (scope_block) os<<"}"<<CRLF;
}
void decorate_latex(OsTy & os,const StrTy & datefmin,const StrTy & datefmax)
{
if (scope_block) os<<"{"<<CRLF;

if ( vars.find("datemin")!=vars.end()) 
	latex_blocks::variable(os,vars ["datemin"], datefmin);
if ( vars.find("datemax")!=vars.end()) 
	latex_blocks::variable(os,vars ["datemax"], datefmax);
if ( vars.find("footnotes")!=vars.end()) 
{
StrTy fn="";
MM_LOOP(ii,footnotes) {fn+=StrTy("{\\bf ")+(*ii).first+StrTy(") } ")+(*ii).second+StrTy(".");}
	latex_blocks::variable(os,vars ["footnotes"], fn);
}
if ( comment.length()!=0) latex_blocks::comment(os,comment,true);
}

template <class Ty>
void decorate(const Ragged & params, const IdxTy&  ds, Ty & flp)
{
 dsz=ds;// ata.size();
imin=0;
 imax=dsz;
 serialmin=0;
 serialmax=0;
 datemin=flp.datemin(); // "";
datemax=flp.datemax(); // "";
period=1;
filtering=1;
tparam["headers"]="1";
tparam["caption"]=" ";
tparam["maxsz"]="30";
tparam["period"]="1";
//tparam["filtering"]="1";
const StrTy & pf=flp.get("period");
if (pf.c_str()[0]!=0) tparam["period"]=pf; // "1";
const StrTy tpp=tparam["period"];
MM_ERR(" period info "<<MMPR2(pf,tpp))
period=myatoi(tpp.c_str());
if (period==0) { period=1; tparam["period"]="1";}
comment="";
scope_block=false;
for(IdxTy i=0; i<params.size(); ++i)
{
const Ragged::Line & line=params.line(i);
const IdxTy ll=line.size();
if (ll<2) continue;
const StrTy&  scmd=line[0];
const StrTy&  skey=line[1];
tparam[scmd]=skey; // added catch all
MM_ONCE(" parameters may shadow with recent devafilat see above ",)
if (scmd=="superscript")
{
for (IdxTy j=2; j<ll; ++j) superscripts[skey].push_back(line[j]);
}

//if (scmd=="dogs") { for (IdxTy j=1; j<ll; ++j) subdogs[line[j]]=subdogs.size() ; }
if (scmd=="dogs") { for (IdxTy j=1; j<ll; ++j) subdogs.add(line[j]) ; }
//if (scmd=="foods") { for (IdxTy j=1; j<ll; ++j) subnames[line[j]]=subnames.size(); }
//if (scmd=="output-foods") { for (IdxTy j=1; j<ll; ++j) subnames[line[j]]=subnames.size(); }
if (scmd=="output-foods") { for (IdxTy j=1; j<ll; ++j) subnames.add(line[j]); }

if (scmd=="units") { for (IdxTy j=2; j<ll; ++j) units[line[j]]=skey ; }
if (scmd=="caption") { tparam["caption"]="";
 concat( tparam["caption"],  line, 2, ll);
// for (IdxTy j=2; j<ll; ++j) tparam["caption"]+=StrTy(" ")+line[j] ; }
}
if (scmd=="plot-food") 
{ tparam["plot-food"]=""; concat( tparam["plot-food"],  line, 2, ll); }

if (scmd=="filter") { if ( ll>3)  m_coefs_map[skey].load(line,2);  }
if (scmd=="variable") { concat( vars[skey],  line, 2, ll); }
if (scmd=="footnote") { 
 concat( footnotes[skey],  line, 2, ll);
// varnames[skey]=line[2] ; 
// for (IdxTy j=2; j<ll; ++j) footnotes[skey]+=StrTy(" ")+line[j] ; 
}


if (scmd=="comment") { 
 for (IdxTy j=2; j<ll; ++j) comment+=StrTy(" ")+line[j] ; }

if (scmd=="scope") {scope_block=true; }  
if (scmd=="nscope") {scope_block=!true; }  
if (scmd=="period") {if (pf.length()==0) period=myatoi(skey);
MM_ERR(" period found "<<MMPR3(pf,period,skey))
 }  
if (scmd=="datemin") {if (flp.datemin()=="") datemin=(skey); }  
if (scmd=="datemax") {if (flp.datemax()=="") datemax=(skey); }  

if (scmd=="serialmin") {serialmin=myatoi(skey); }  
if (scmd=="serialmax") {serialmax=myatoi(skey); }  
if (scmd=="filtering") {filtering=myatoi(skey); }  

}
// _check_x=subdogs.selecting(); // (subdogs.size()!=0);
// _check_y=subnames.selecting(); // (subnames.size()!=0);

} // decorate

bool oor( const IdxTy x, const IdxTy x0, const IdxTy xm)
{
if (x<x0) return true;
if (xm!=0) if (x>=xm) return true;
return false;

}
bool serial_oor(const IdxTy serial){
return oor(serial,serialmin,serialmax); 
}
bool date_oor(const StrTy & date){
if (date<datemin) return true; // continue;
//if (datemax.c_str()[0]!=0) if (date>=datemax) continue;
if (datemax.c_str()[0]!=0) if (date>=datemax) return true;
return false;
}


}; // table_decorators
typedef table_decorators::TableParams TableParams;
class bin_meta_data
{
public:
typedef std::map<StrTy,IdxTy  > ExMap;
bin_meta_data() {Init(); }

void Init() { m_days=0; }

void data(const IdxTy serial, const StrTy & dog, const IdxTy days)
{
if (m_days==0) m_days=days;
if (m_days!=days)
{
MM_ERR(" day change "<<MMPR4(serial,dog,days,m_days))
}
++bindaze[serial];
++bindogdaze[dog][serial];
// this is really stpid as bindogdaze should be sorted. 
extreme(bindogmax, serial,true,dog);
extreme(bindogmin, serial,!true,dog);

}
const bool entries(const StrTy & d)   { return  bindogdaze[d].size()!=0; }
const IdxTy size(  const StrTy & d)   { return  bindogdaze[d].size(); }
const IdxTy first(const StrTy & d)   { return  bindogmin[d]; }
const IdxTy last(const StrTy & d)   { return  bindogmax[d]; }

void extreme(ExMap & m, const IdxTy s, const bool p, const StrTy & d)
{
auto ii=m.find(d);
bool found=(ii!=m.end());
if ( ! found ) { m[d]=s; return; } 
if ( p ) if (m[d]<s) { m[d]=s; return;} // max if true 
if (m[d]>s) m[d]=s; 
}

std::map<IdxTy,IdxTy  > bindaze;
ExMap bindogmin,bindogmax;
std::map< StrTy, std::map<IdxTy,IdxTy > >  bindogdaze;
IdxTy m_days;

}; // bin_meta_data

// this is the format coming from mjm_snacks. Roughly date, dog, and food/dose pairs
// reorganized before handing to dog_plots
class snacks_format
{
typedef qty_ty Val;
//typedef std::map<StrTy, D> KVTy;
typedef std::map<StrTy, Val> KVTy;
public:

snacks_format(){ Init(); Reset();}
bool use_units() const { return have_units; } 
void use_units(const bool tf) { have_units=tf; del_idx=tf?3:2; } 
template <class Ty, class Td> bool line(const Ty & line, Td & td )
{
Reset();
nf=line.size();
if (nf<nfmax) { return false; }
serial=myatoi(line[serial_idx]);
// this mayhave been meant as the end but AFAICT the serial format is daily 
//serial=myatoi(line[serial_idx]);
if (td.serial_oor(serial)) { return false;}
 date=(line[date_idx]);
if (td.date_oor(date)) {return false; } ;

catalog_date_serial(date,serial);



  dog=(line[name_idx]);
 //if (!td.check_x(dog)) { return false; } 
 if (!td.subdogs.test(dog)) { return false; } 
for(IdxTy j=fields_idx; j<nf; j+=del_idx) { 

// the unit is being ignored with atof, add to keyu 
#if 0
const StrTy & key=line[j];
//if (!td.check_y(key)) continue;   
 if (!td.subnames.test(key)) { return false; } 
StrTy v=line[j+1];
if (have_units) v=v+line[j+2];
MM_ERR(" stubbed "<<MMPR(v))
#endif

StrTy  key=line[j];
if (have_units) key=key+"("+line[j+2]+")";
//if (!td.check_y(key)) continue;   
 if (!td.subnames.test(key)) { return false; } 
StrTy v=line[j+1];
//if (have_units) v=v+line[j+2];
MM_ERR(" stubbed "<<MMPR2(key,v))


//const D & v=atof(line[j+1].c_str());
kv[key]+=v;

}
return true; 
}


template <class Ty> bool line(const Ty & line )
{
Reset();
nf=line.size();
if (nf<nfmax) { return false; }
serial=myatoi(line[serial_idx]);
if (serial==0) 
{MM_ERR(" zero zerial "<<MMPR2(line.size(),line[serial_idx]))  
return false; 
}
 date=(line[date_idx]);
catalog_date_serial(date,serial);
/*if (date_to_serial.find(date)!=date_to_serial.end())
{
if (date_to_serial[date]!=serial)
{
MM_ERR(" dates not const "<<MMPR3(serial, date_to_serial[date],date))
}
}
date_to_serial[date]=serial;
serial_to_date[serial]=date;
*/

  dog=(line[name_idx]);
for(IdxTy j=fields_idx; j<nf; j+=del_idx) { 
// the unit is being ignored with atof, add to keyu 
#if 0
const StrTy & key=line[j];
StrTy v=line[j+1];
if (have_units) v=v+line[j+2];
#endif
 StrTy key=line[j];
StrTy v=line[j+1];
if (have_units) key=key+"("+line[j+2]+")";;

// const D & v=atof(line[j+1].c_str());
kv[key]+=v;
}
return true; 
}


void catalog_date_serial(const StrTy & date,const IdxTy & serial)
{

if (date_to_serial.find(date)!=date_to_serial.end())
{
if (date_to_serial[date]!=serial)
{
MM_ERR(" dates not const "<<MMPR3(serial, date_to_serial[date],date))
}
}
date_to_serial[date]=serial;
serial_to_date[serial]=date;

}




void Init()
{
 serial_idx=0;
 date_idx=1;
 name_idx=2;
 fields_idx=3;
have_units=false;
del_idx=2;
 nfmax=3;
date_to_serial.clear();
serial_to_date.clear();
}
void Reset()
{
nf=0;
serial=0;
date="";
dog="";
kv.clear();

}
 IdxTy serial_idx;
 IdxTy date_idx;
 IdxTy name_idx;
 IdxTy fields_idx;
bool have_units;
IdxTy del_idx;
 IdxTy nfmax;
std::map<StrTy, IdxTy> date_to_serial;
std::map<IdxTy,StrTy> serial_to_date;

 IdxTy nf,serial;
StrTy date,dog;
KVTy kv;

}; //  snacks_format

class day_plot_data_table
{
public:
typedef std::map<StrTy, DDstat> Bins;
typedef std::map<StrTy, bin_meta_data> Bmd;
day_plot_data_table(): m_serial_min(~0), m_serial_max(0){}
void have_units(const bool tf) {sf.use_units(tf); } 
void load(const Ragged & data) // , const IdxTy period)
{
MM_ERR(" loading dog_plot_data_table"<<MMPR2(td.imin,td.imax))
const IdxTy period=td.period;
lines=0;
for(IdxTy i=td.imin; i<td.imax; ++i)
{
const Ragged::Line & line=data.line(i);
if (!sf.line(line,td) ) continue ;
++lines;
++dates[sf.date];
const IdxTy & serial=sf.serial;
if (serial>m_serial_max) m_serial_max=serial;
if (serial<m_serial_min) m_serial_min=serial;
++serials[serial];
const StrTy & dog=sf.dog;
MM_ERR( "loading "<<MMPR3(dog,sf.date,sf.serial))
++dogs[dog];
//Ss ss;
//ss<<IdxTy(serial/period);
const StrTy bin=make_bin(serial,period); // ss.str();
bmd[bin].data(sf.serial, dog,period);
//++bindaze[bin][sf.serial];
//++bindogdaze[bin][dog][sf.serial];
MM_LOOP(ii,sf.kv) 
{ 
++foods[(*ii).first]; // moved from generate???? 
bins[bin][dog][(*ii).first].add_dose((*ii).second,sf.serial); }
 } // i 
MM_ERR(" parsing "<<MMPR3(lines, dates.size(),dogs.size()))
} // load 

StrTy make_bin(const IdxTy serial, const IdxTy period)
{
Ss ss;
ss<<IdxTy(serial/period);
return ss.str();
}


//ConfigRules rl; // xxxxx
void load_filtered(ConfigRules & rl, const Ragged & data) // , const IdxTy period)
{
//if (scmd=="filter") { if ( ll>3)  m_coefs_map[skey].load(line,2);  }
MM_ERR(" loading dog_plot_data_table"<<MMPR2(td.imin,td.imax))
const IdxTy period=td.period;
const bool do_filtering=(td.filtering!=0);
lines=0;
for(IdxTy i=td.imin; i<td.imax; ++i)
{
const Ragged::Line & line=data.line(i);
if (!sf.line(line,td) ) continue ;
++lines;
++dates[sf.date];
const IdxTy & serial=sf.serial;
if (serial>m_serial_max) m_serial_max=serial;
if (serial<m_serial_min) m_serial_min=serial;
++serials[serial];
const StrTy & dog=sf.dog;
MM_ERR( "loading "<<MMPR3(dog,sf.date,sf.serial))
++dogs[dog];
const StrTy bin=make_bin(serial,period); // ss.str();
if (!do_filtering) bmd[bin].data(sf.serial, dog,period);
MM_LOOP(ii,sf.kv) 
{ 
 StrTy food=(*ii).first;
D scale=1.0;
//foodalias
bool apply_food_xform=true;
if ( apply_food_xform)
{
std::vector<StrTy> desc;
if ( rl.apply(desc,"foodalias",food))
{
const IdxTy dsz=desc.size();
if ( dsz>0) food=desc[0];
if ( dsz>1) { scale=atof(desc[1].c_str()); } 
} // rl.apply
// stop now or it gets included
if (scale==0) continue; 
} // apply

++foods[food]; // moved from generate???? 
//++foods[(*ii).first]; // moved from generate???? 
const StrTy name=dog+StrTy("-")+food;

if ( do_filtering) {
StrTy fname;
rl.apply(fname,"path-filter",name);
//rl.apply(m_attmap["color"], "color",name());
 Coefficients & f=td.m_coefs_map[fname];
MM_ERR(MMPR4(fname,f.dump(),td.m_coefs_map.dump(),rl.dump()))
const bool hard_code=false;
if (hard_code)
{
const int lim=3;
const D w=1.0/D(2*lim+1);
for(int k=-lim; k<=lim; ++k)
{ 
const IdxTy serialk=int(serial)+k;
const StrTy bink=make_bin(serialk,period); // ss.str();
bmd[bink].data(serialk, dog,period);
bins[bink][dog][(*ii).first].add_dose((*ii).second,serialk,w);

} //k
} // hard_code
else
{
if (period!=1) { MM_ERR(" non unity period odd results") } 
auto iic=f.begin();
auto eec=f.end();
while (iic!=eec) 
{

const IdxTy serialk=int(serial)+((*iic).first);
const StrTy bink=make_bin(serialk,period); // ss.str();
bmd[bink].data(serialk, dog,period);
//bins[bink][dog][(*ii).first].add_dose((*ii).second,serialk,(*iic).second);
//bins[bink][dog][food].add_dose(scale*(*ii).second,serialk,(*iic).second);
bins[bink][dog][food].add_dose((*ii).second*scale,serialk,(*iic).second);

++iic;
} // iic

} // !hard_code

} // do_filtering

//else { bins[bin][dog][(*ii).first].add_dose((*ii).second,sf.serial);}
else { bins[bin][dog][food].add_dose((*ii).second,sf.serial);}

}


 } // i 


MM_ERR(" parsing "<<MMPR3(lines, dates.size(),dogs.size()))
} // load_filtered



void fill_missing_points()
{
const IdxTy period=td.period;
MM_ERR(" adding zero for missing day/dog/food points")
for(IdxTy serial=m_serial_min; serial<=m_serial_max; ++serial)
{
const StrTy bin=make_bin(serial,period); // ss.str();
auto & b=bins[bin];
MM_LOOP(jj,dogs)
{
const StrTy dog=(*jj).first;
bmd[bin].data(serial, dog,period);
auto & b2=b[dog];
MM_LOOP(kk,foods)
{
const StrTy food=(*kk).first;
//bins[bin][dog][food].add_dose(0,serial); }
MM_ERR(" filling missing doses over writes good ones too ")
b2[food].add_dose(0,serial); 
} // food
} // dogs
} // serial 
}

void calculate(const bool use_sample_days=false)
{
MM_ONCE(" day persiod "<<MMPR(use_sample_days),)
MM_LOOP(ii,bins)
{ 
//const IdxTy ndays=bindaze[(*ii).first].size(); // ?? property of bin or data?? 
MM_LOOP(jj,(*ii).second)
{
//const IdxTy ndogdays=bindogdaze[(*ii).first][(*jj).first].size(); // ?? property of bin or data?? 
const IdxTy ndogdays=bmd[(*ii).first].bindogdaze[(*jj).first].size(); // ?? property of bin or data?? 
const IdxTy ncalendardays=bmd[(*ii).first].m_days; // ?? property of bin or data?? 
MM_LOOP(kk,(*jj).second)
{
//(*kk).second.calculate(ndays);
//++foods[(*kk).first]; // moved to load???? 
if (use_sample_days) (*kk).second.calculate(ndogdays);
else (*kk).second.calculate(ncalendardays);
}
} // jj 
} /// ii 
} // calculate 
void generate(OsTy & os, Ragged & ro,const bool dump_all,const bool write_svg)
{

if (dump_all||write_svg)
{
const StrTy sep=" ";
MM_LOOP(ii,bins)
{
const StrTy & bin=(*ii).first;
const IdxTy day=myatoi(bin)*td.period;
if (false) MM_ERR(MMPR2(bin,day))

MM_LOOP(jj,(*ii).second)
{
const StrTy & dog=(*jj).first;
IdxTy  pereff= bmd[bin].size(dog);

 IdxTy dayfirst= bmd[bin].first(dog);
 IdxTy daylast= bmd[bin].last(dog);
StrTy lbl=sf.serial_to_date[daylast];
StrTy lblfirst=sf.serial_to_date[dayfirst];
Ss ss;
//ss<<sf.serial_to_date[dayfirst]<<" ("<<pereff<<" days) ";
//ss<<dayfirst<<" "<<daylast<<" "<<lbl<<" "<<pereff<<" ";
//ss<<dayfirst<<" "<<daylast<<" "<<lblfirst<<" "<<lbl<<" "<<pereff<<" ";
//                                                                ^^^^^ space shift 
ss<<dayfirst<<" "<<daylast<<" "<<lblfirst<<" "<<lbl<<" "<<pereff;
 StrTy  binlab=ss.str();

MM_LOOP(kk,foods)
{
const StrTy & food=(*kk).first;
const Stat & s= bins[bin][dog][food];
{
Ss ss;
// FIXME this is dumb, the dump does not need to be consistent with ro
// it shoudl derive output from that unless intent is to contrast.
if (dump_all){
ss<<s.dump();
os<<bin<<sep<<binlab<<sep<<dog<<sep<<food<<sep<<ss.str()<<CRLF;
}
if (write_svg)
{
ss<<s.dump(1);
Ss sv;
sv<<bin<<sep<<binlab<<sep<<dog<<sep<<food<<sep<<ss.str()<<CRLF;
ro.add(sv.str());
}


} // 

} // kk 
} // jj 
} // ii 
} // dump_all || write_svg
} // generate


snacks_format sf;
table_decorators td;
std::map<StrTy, IdxTy> dogs,dates,foods;
std::map<IdxTy, IdxTy> serials;
Bmd bmd;
Bins bins;
IdxTy lines;
IdxTy m_serial_min, m_serial_max;


}; // day_plot_data_table



static bool flag_bit(const IdxTy & flags,const IdxTy & n,const bool t)
{
const bool f=((flags&(1<<n))!=0); 
return t?f:!f;
}


void cmd_snacks_txt_svg_i(Cip & cip , LocalVar & lv )
{

const StrTy dest=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy dname=cip.wif(3);
const StrTy pname=cip.wif(4);
const StrTy pmarkname=cip.wif(5);
Ragged & data=m_ragged_map[dname];
const Ragged & params=m_ragged_map[pname];
const bool dump_all=flag_bit(flags,0,true);
const bool write_svg=flag_bit(flags,1,true);
const bool fill_missing_points=flag_bit(flags,2,true);
const bool make_from_scratch=flag_bit(flags,3,!true);
const bool write_ro=flag_bit(flags,4,true);
const bool dump_dl=flag_bit(flags,5,true);
const bool loading_filtered=flag_bit(flags,6,true);
const bool split_foods=flag_bit(flags,7,true);
const bool have_units=flag_bit(flags,8,true);

if ((cip.help()||(dest.c_str()[0]==0)))
{
Ss ss;
ss<<" cmd_snacks_txt_svg_i dest flags dname pname "<<CRLF;
ss<< " make various svg plots of dose stats per dog over various intervals"<<CRLF;
ss<< " dest : destination file name for output"<<CRLF;
ss<< " flags : bits set for various actions:"<<CRLF;
ss<<"        2:svg,"<<CRLF;
ss<<"        8:input is  parsed, skip parsing step conerting from multipe foods to single dog/food/day per line   "<<CRLF;
ss<<"       16: write parsed. requies 2 to work "<<CRLF;
ss<<"       64:moving average enable   "<<CRLF;
ss<<"      128:make one plot per food, most useful with al foods enabled    "<<CRLF;

ss<<"  ./mjm_linc_graph.out -cmd \"read-ragged in $SE_DATA_EXCHANGE/snacks_collated.ssv\" -cmd \"read-ragged p pin2.txt\"  -cmd \"snacks-txt-svg-i happy2.svg 2 in p\" -quit " <<CRLF;
ss<<" loads a data file derived from mjm_snacks and creates an intermediate format "<<CRLF;
ss<<" with one dog-time-dose point per line instead of the multiple doses on each dog-day line"<<CRLF;
ss<<" . With moving avg set period to one "<<CRLF;
cip.help(ss.str());
return; 
}
//snacks_format sf;
Ragged rtemp;
Ragged &  ro=(make_from_scratch)?rtemp:data;
MM_ERR( "cmd_snacks_txt_svg "<<MMPR4(dest,flags,dname,pname)<<MMPR2(data.size(),params.size())<<MMPR3(dump_all,write_svg,fill_missing_points)<<MMPR2(make_from_scratch,write_ro)<<MMPR(have_units))
// the normal route unless the dog-period-food per line file has been genrated already. 
if ( make_from_scratch)
{
std::ofstream ofs(dest);
std::ostream & os= ofs;// std::cout;
day_plot_data_table dpdt;
dpdt.have_units(have_units);

ConfigRules rl; // xxxxx
rl.add_lines(params.begin(), params.end()); //adfasdfadfs


dpdt.td.decorate(params,data.size(),m_flp);
if( loading_filtered) dpdt.load_filtered(rl, data);
else dpdt.load(data);
if (fill_missing_points) dpdt.fill_missing_points();
dpdt.calculate();
// ro spaces vary depending on all foods or not
dpdt.generate(os,ro,dump_all,write_svg);
}
ro.write_file("xxx.kluge",0);
// alternative output to svg but svg needs to be specified in flags too 
if (write_ro)
{
//std::ofstream ofs(dest);
//std::ostream & os= ofs;// std::cout;
//Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ;
const IdxTy wflags=0;
ro.sep(" ");
//if (tabsep) r.sep("\t");
//if (csv) r.sep(",");
//if (ignore_hash) r.ignore_hash(ignore_hash);
MM_ERR(" writing dat file "<<MMPR2(dest,ro.size()))
ro.write_file(dest,wflags);
}

if (write_svg&&!write_ro)
{

Ragged p2=params;
// markups globbed into params but should be separate file 
Ragged pmark=params;
if (pmarkname.length()!=0) pmark+=m_ragged_map[pmarkname];
MM_ERR(MMPR2(p2.size(),params.size()))
// the period here is not consistent anymore. 
ReadWriteMap rwm;
// this can be hierarchially resolved, except for the things like the markups
// first make this map into a temp 
p2.to_map(rwm);
// starting with m_flp_def, then add the p2 map and finally the m_flp contents;
//CoRes cores;
//rwm.map()=cores.def_config_cmd(m_flp_def.map().map(),rwm.map(),m_flp.map().map());
rwm.map()=m_lo.def_config_cmd(m_flp_def.map().map(),rwm.map(),m_flp.map().map());
//rwm.map()=cores.def_config_cmd(m_flp_def.map().map(),p2,m_flp.map().map());
Names fng;
RagSplit rs;
if (split_foods)
{
fng.set_file(dest);
rs.field(7);
rs.split(ro);
}
auto ii=rs.begin();
do { 
if (split_foods) if (ii==rs.end()) break;
StrTy desteff=dest;
if (split_foods) desteff=fng.sfx((*ii).first);
MM_MSG(MMPR(desteff));
std::ofstream ofs(desteff);
std::ostream & os= ofs;// std::cout;
Ragged & roe=(split_foods)?((*ii).second) : ro;
if (split_foods) ++ii; 
mjm_dog_plots lb;
//lb.setup_day_plot(ro,p2);
// this works but ideally should be distinct from the parameters and configuration 
//lb.set_markups(p2);
lb.set_markups(pmark);
lb.setup_day_plot(roe,rwm,p2);
//lb.find_day_statistics(ro);
//lb.sort(ro);
mjm_svg_writer sw;
//std::ofstream  os(fn);
lb.write_day_svg(os,sw);
if (dump_dl) MM_ERR(MMPR(lb.dump_dl())) 
if (!split_foods) break; 
} while (true);
} // write_svg


} // cmd_snacks_txt_svg_i

void cmd_snacks_txt_svg(Cip & cip , LocalVar & lv )
{

const StrTy dest=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy dname=cip.wif(3);
const StrTy pname=cip.wif(4);
const Ragged & data=m_ragged_map[dname];
const Ragged & params=m_ragged_map[pname];
const bool dump_all=flag_bit(flags,0,true);
const bool write_svg=flag_bit(flags,1,true);

//snacks_format sf;

std::ofstream ofs(dest);
std::ostream & os= ofs;// std::cout;
MM_ERR( "cmd_snacks_txt_svg "<<MMPR4(dest,flags,dname,pname)<<MMPR2(data.size(),params.size()))
day_plot_data_table dpdt;
//table_decorators td;
dpdt.td.decorate(params,data.size(),m_flp);
//typedef std::map<StrTy, DDstat> Bins;
//typedef std::map<StrTy, bin_meta_data> Bmd;

//const IdxTy period=dpdt.td.period;
//IdxTy lines=0;
dpdt.load(data);

dpdt.calculate();

Ragged ro;
dpdt.generate(os,ro,dump_all,write_svg);

if (write_svg)
{
layout_blocks lb;
Ragged p2=params;
MM_ERR(MMPR2(p2.size(),params.size()))
lb.setup_day_plot(ro,p2);
lb.find_day_statistics(ro);
//lb.sort(ro);
mjm_svg_writer sw;
//std::ofstream  os(fn);
lb.write_day_svg(os,sw);
} // write_svg


} // cmd_snacks_txt_svg

// unused, skeleton for making tables for particular amino acids or ineractions etc 
void cmd_make_latex_table(Cip & cip , LocalVar & lv )
{
MM_ONCE(" skeleton code unused yet ",)
const StrTy dest=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy dname=cip.wif(3);
const StrTy pname=cip.wif(4);
const Ragged & data=m_ragged_map[dname];
const Ragged & params=m_ragged_map[pname];

const IdxTy serial_idx=0;
const IdxTy date_idx=1;
const IdxTy name_idx=2;
const IdxTy fields_idx=3;
const IdxTy nfmax=3;

std::ofstream ofs(dest);
std::ostream & os= ofs;// std::cout;
MM_ERR( "cmd_make_latex_table "<<MMPR4(dest,flags,dname,pname)<<MMPR2(data.size(),params.size()))
table_decorators td;
td.decorate(params,data.size(),m_flp);

//const bool check_dogs=td.subdogs.selecting(); // (td.subdogs.size()!=0);
//const bool check_foods=td.subnames.selecting(); // (td.subnames.size()!=0);

DDstat dds;
std::map<StrTy, IdxTy> dogs,names;
IdxTy lines=0;
 StrTy datefmin="z";
 StrTy datefmax="";


for(IdxTy i=td.imin; i<td.imax; ++i)
{
const Ragged::Line & line=data.line(i);
const IdxTy nf=line.size();
if (nf<nfmax)
{

continue;
}
const IdxTy serial=myatoi(line[serial_idx]);

if (td.serial_oor(serial)) continue;

const StrTy&  date=(line[date_idx]);
if (td.date_oor(date)) continue;

if (date<datefmin) datefmin=date;
if (date>datefmax) datefmax=date;

const StrTy&  name=(line[name_idx]);
//if (check_dogs) { if (!td.check_x(name)) continue;  } 
//if (check_dogs) 
{ if (!td.subdogs.test(name)) continue;  } 
++dogs[name];
++lines;
for(IdxTy j=fields_idx; j<nf; j+=2) { 
const StrTy & key=line[j];
//if (check_foods) { if (!td.check_y(key)) continue;  } 
{ if (!td.subnames.test(key)) continue;  } 
const D & v=atof(line[j+1].c_str());
dds[name][key].add_dose(v,serial);
++names[key];
}
} // i 
// parsing phase 
const IdxTy tlines=names.size();
const IdxTy twords=dogs.size();
Ragged rtable(tlines+1,twords+1);

IdxTy iname=0;
IdxTy idog=0;
MM_LOOP(jj,dogs)
{rtable.ncline(iname)[idog]=(*jj).first;
++idog;
}
++iname;
MM_LOOP(ii,names)
{
StrTy  name=(*ii).first;
StrTy  nm=name;

td.decorate_y(nm,name);

IdxTy idog=0;
rtable.ncline(iname)[idog]=name;
++idog;
MM_LOOP(jj,dogs)
{

const StrTy & dog=(*jj).first;
// cleverly this went from drug to dog "name" lol. 
Stat &   q= dds[dog][nm];
//q.calculate(dogs[name]);
q.calculate((*jj).second);
StrTy x=q.to_string();
rtable.ncline(iname)[idog]=x;

++idog;
} // jj 

++iname;
} // ii 
Ss cs;
cs<<" From "<<datefmin<<" to "<<datefmax<<" ";
td.tparam["caption"]=cs.str()+td.tparam["caption"];

MM_ERR(" writing latex table "<<MMPR3(cs.str(),rtable.size(),td.tparam.size()))

td.decorate_latex(os,datefmin,datefmax);


latex_blocks::ragged_table(os,rtable,td.tparam);
td.decorate_latex(os);


} // cmd_make_latex_table


void cmd_snacks_rags(Cip & cip , LocalVar & lv )
{

const StrTy dest=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy dname=cip.wif(3);
const StrTy pname=cip.wif(4);
const Ragged & data=m_ragged_map[dname];
const Ragged & params=m_ragged_map[pname];

const bool have_units=Bit(flags,0);
if (cip.help()||(dest.length()==0)) {
Ss ss; 
ss<<"  cmd_snacks_rags : make latex table for various snack dog food groups"<<CRLF;
ss<< "snacks-rags dest flags dname pname data params "<<CRLF;
ss<<" dest : output file name "<<CRLF;
ss<<" flags : 0, have units  "<<CRLF;
ss<<" dname : data rags name  "<<CRLF;
ss<<" pname : parameter rags  name "<<CRLF;
ss<<" ./mjm_linc_graph.out -cmd \"read-ragged in $SE_DATA_EXCHANGE/snacks_collated.ssv\" -cmd \"read-ragged p happy.txt\" -cmd \"snacks-rags htable 0 in p\" -quit"<<CRLF;

cip.help(ss.str());
}

const IdxTy serial_idx=0;
const IdxTy date_idx=1;
const IdxTy name_idx=2;
const IdxTy fields_idx=3;
const IdxTy nfmax=3;

std::ofstream ofs(dest);
std::ostream & os= ofs;// std::cout;
MM_ERR( "cmd_snacks_rags "<<MMPR4(dest,flags,dname,pname)<<MMPR3(have_units,data.size(),params.size()))
table_decorators td;

ConfigRules rl; // xxxxx
rl.add_lines(params.begin(), params.end()); //adfasdfadfs

td.decorate(params,data.size(),m_flp);
DDstat dds;
IdxTy lines=0;

 StrTy datefmin="z";
 StrTy datefmax="";


for(IdxTy i=td.imin; i<td.imax; ++i)
{
const Ragged::Line & line=data.line(i);
const IdxTy nf=line.size();
if (nf<nfmax) { continue; }
//typedef std::map<StrTy, D> Qmap;
//Qmap m;
const IdxTy serial=myatoi(line[serial_idx]);
// in theory this can compute the start but
// not a big deal 
//if (serial<serialmin) continue;
// this should break here but can flag and check data 
//if (serialmax!=0) if (serial>=serialmax) continue;

if (td.serial_oor(serial)) continue;

const StrTy&  date=(line[date_idx]);
if (td.date_oor(date)) continue;

if (date<datefmin) datefmin=date;
if (date>datefmax) datefmax=date;

const StrTy&  name=(line[name_idx]);
if ( !td.subdogs.test(name)) continue;
++lines;
//const bool have_units=true;
MM_ONCE(" have_units is now flag bit 0 "<<MMPR(have_units),)
const IdxTy field_inc=(have_units)?3:2;
//for(IdxTy j=fields_idx; j<nf; j+=2) { 
for(IdxTy j=fields_idx; j<nf; j+=field_inc) { 
StrTy key=line[j];
//if ( !td.subnames.test(key)) continue;
bool apply_food_xform=true;
if ( apply_food_xform)
{
D scale=1.0;
std::vector<StrTy> desc;
if ( rl.apply(desc,"foodalias",key))
{
const IdxTy dsz=desc.size();
if ( dsz>0) key=desc[0];
if ( dsz>1) { scale=atof(desc[1].c_str()); } 
} // rl.apply
// stop now or it gets included
if (scale==0) continue; 
if ( !td.subnames.test(key)) continue;
// adding for the new units logic... 
if ((j+1)>=nf)
{
Ss ss;
MM_LOOP(ee,line) { ss<<(*ee)<<" "; } 
MM_ERR(" probably units error "<<MMPR4(j,nf,key,ss.str()))
}
else {
 D  v=scale*atof(line[j+1].c_str());
dds[name][key].add_dose(v,serial);
} // units logic error test 

} // if apply
else { 
if ( !td.subnames.test(key)) continue;
const D & v=atof(line[j+1].c_str());
dds[name][key].add_dose(v,serial);
} // no xform 

}
} // i 
// parsing phase 
const IdxTy tlines=td.subnames.size(); // names.size();
const IdxTy twords=td.subdogs.size(); // dogs.size();
Ragged rtable(tlines+1,twords+1);

IdxTy iname=0;
IdxTy idog=0;
//rtable.ncline(iname)[idog]="Ingredient";
MM_LOOP(jj,td.subdogs)
{rtable.ncline(iname)[idog]=(*jj).first;
++idog;
}
++iname;
MM_LOOP(ii,td.subnames)
{
StrTy  name=(*ii).first;
StrTy  nm=name;

td.decorate_y(nm,name);

IdxTy idog=0;
rtable.ncline(iname)[idog]=name;
++idog;
MM_LOOP(jj,td.subdogs)
{

const StrTy & dog=(*jj).first;
// cleverly this went from drug to dog "name" lol. 
Stat &   q= dds[dog][nm];
q.calculate((*jj).second);
StrTy x=q.to_string();
rtable.ncline(iname)[idog]=x;

++idog;
} // jj 

++iname;
} // ii 
Ss cs;
cs<<" From "<<datefmin<<" to "<<datefmax<<" ";
td.tparam["caption"]=cs.str()+td.tparam["caption"];

MM_ERR(" writing latex table "<<MMPR3(cs.str(),rtable.size(),td.tparam.size()))

td.decorate_latex(os,datefmin,datefmax);
latex_blocks::ragged_table(os,rtable,td.tparam);
td.decorate_latex(os);

} // cmd_snacks_rags

//////////////////////////////////////////////////////////////////////////
#if 0
ls ../mikemail/dev/out/
dog_daily.txt  dog_glob.txt  dog_used.txt
marchywka@happy:/home/documents/cpp/proj/linc_graph$ cat  ../mikemail/dev/out/ *.txt
18551 2020-10-16 Beauty sugar 900 tsp salt 1 tsp

18551 2020-10-16 9 Beauty sugar - 900 tsp
18551 2020-10-16 9 Beauty salt - 1 tsp

sugar 900 tsp
salt 1 tsp
#endif
//////////////////////////////////////////////////////////////////
typedef std::map<StrTy , std::map<IdxTy, IdxTy > >  DopType;


class misc_format_params
{
public:
misc_format_params() { Init(); } 

void Init()
{

//datefmin datefmax
lines=0;
ncol=0;
//have_units append_units append_units_after
serial_idx=0; date_idx=1; name_idx=2; fields_idx=3; nfmax=3;
split_columns=(ncol!=0);
ditem=0;
group_idx=0;
} // Init 
IdxTy ncol,serial_idx, date_idx, name_idx, fields_idx,nfmax,group_idx;
StrTy datefmin, datefmax;
IdxTy lines,ditem; 
bool split_columns, have_units,  append_units,  append_units_after; 


}; // misc_format_params


typedef misc_format_params MiscForm;
typedef std::map< StrTy , std::map<StrTy, int> > GroupMap;


void parse_daily_rag( DopType & dop, DDstat & dds, table_decorators & td
,GroupMap &gm, MiscForm & mf,  ConfigRules & rl, const Ragged & data ) // const
{
//const bool split_columns=(ncol!=0);


for(IdxTy i=td.imin; i<td.imax; ++i)
{
const Ragged::Line & line=data.line(i);
const IdxTy nf=line.size();
if (nf<mf.nfmax) { continue; }
//typedef std::map<StrTy, D> Qmap;
//Qmap m;
const IdxTy serial=myatoi(line[mf.serial_idx]);
// in theory this can compute the start but
// not a big deal 
//if (serial<serialmin) continue;
// this should break here but can flag and check data 
//if (serialmax!=0) if (serial>=serialmax) continue;

if (td.serial_oor(serial)) continue;

const StrTy&  date=(line[mf.date_idx]);
if (td.date_oor(date)) continue;
StrTy datebin;
if ( !rl.apply(datebin,"datealias",date)) continue; 
if ( !td.subdates.test(datebin)) continue;

if (date<mf.datefmin) mf.datefmin=date;
if (date>mf.datefmax) mf.datefmax=date;



const StrTy&  name=(line[mf.name_idx]);
if ( !td.subdogs.test(name)) continue;
MM_MSG(" binning "<<MMPR4(name,date,datebin,serial))
++dop[datebin][serial];
++mf.lines;
const IdxTy fields_del=mf.ditem+(mf.have_units?3:2);
// just check once doh
int nfd=int(nf)-int(mf.fields_idx);
int left_over=nfd%fields_del;
if (left_over)
{
Ss ss; MM_LOOP(zz,line) ss<<(*zz)<<" ";
MM_ERR("line size wrong check units flats "<<MMPR4(left_over,mf.fields_idx,nf,ss.str()))
} // left_over
StrTy key=" ";
const IdxTy nfeff=(mf.group_idx==0)?nf:mf.group_idx;
for(IdxTy j=mf.fields_idx; j<nfeff; j+=fields_del) { 
//for(IdxTy j=fields_idx; j<nf; j+=2) { 
//18551 2020-10-16 9 Beauty sugar - 900 tsp
//StrTy key=line[j];
key=line[j];
if (mf.append_units ) 
{
const StrTy & u=line[j+2+mf.ditem];
if (u!="-") key=key+StrTy("(")+u+StrTy(")");
} // append_units

//if ( !td.subnames.test(key)) continue;
if ((j+1)>=nf)
{
Ss ss; MM_LOOP(zz,line) ss<<(*zz)<<" ";
MM_ERR(" will bomb now, may need to set flag bit zero for units "<<MMPR3(j,nf,ss.str()))

}
bool apply_food_xform=true;
if ( apply_food_xform)
{
D scale=1.0;
const StrTy keyorig=key;
std::vector<StrTy> desc;
if ( rl.apply(desc,"foodalias",key))
{
const IdxTy dsz=desc.size();
if ( dsz>0) key=desc[0];
if ( dsz>1) { scale=atof(desc[1].c_str()); } 
} // rl.apply
// stop now or it gets included
if (scale==0) continue; 
const StrTy keyo=key;
if (mf.append_units_after ) 
{
const StrTy & u=line[j+2+mf.ditem];
if (u!="-") key=key+StrTy("(")+u+StrTy(")");
} // append_units_after
//if ( !td.subnames.test(key)) continue;
if ( !td.subnames.test(keyorig)) continue;
 D  v=scale*atof(line[j+1+mf.ditem].c_str());
dds[datebin][key].add_dose(v,serial);
td.subnames.key_name(key,keyo);
//dds[datebin][key].key_name(keyo);
} // if apply
else { 
if ( !td.subnames.test(key)) continue;

const D & v=atof(line[j+1+mf.ditem].c_str());
dds[datebin][key].add_dose(v,serial);
} // no xform 
} // j ??? 

MM_ERR(" WTF "<<MMPR3(nf,left_over,mf.group_idx))
if ( key.length())
{
//for(IdxTy j=(nf-left_over); j<mf.group_idx; ++j) { 
for(IdxTy j=(mf.group_idx); j<nf; ++j) { 
MM_ERR(" adding group "<<MMPR4(j,nf,left_over,mf.group_idx) <<MMPR2(key,line[j]))
++gm[key][line[j]]; } // j 
}// key.length() // mf.group_idx




} // i


} // parse_daily_rag

class char_defs_class
{

public:
char_defs_class()  {}
char_defs_class(const StrTy & s) : m_name(s),m_latex(s) {}
operator StrTy() const { return m_name; } 
operator IdxTy() const { return myatoi(m_value); } 
const StrTy&  latex () const { return m_latex; } 
const StrTy&  value () const { return m_value; } 
const StrTy&  latex (const StrTy & x )  {m_latex=x;  return m_latex; } 
const StrTy&  value (const StrTy & x )  {m_value=x;  return m_value; } 
StrTy m_name;
StrTy m_latex;
StrTy m_value;

}; // char_defs_class

typedef char_defs_class CharDefsTy;
//typedef std::map<StrTy, StrTy> CharDefMap;
typedef std::map<StrTy, CharDefsTy> CharDefMap;
typedef std::map<StrTy, IdxTy> CharMap;
typedef std::vector<IdxTy> PermTy;


void load_char_defs(CharDefMap & char_defs, const Ragged & r) const
{
MM_LOOP(ii,r)
{
const Ragged::Line & l=(*ii);
if (l.size()<3) continue;
if (l[0]=="grouping")
{
const StrTy& nm=l[1];
const StrTy& v=l[2];
CharDefsTy x(nm);
x.value(v);
if (l.size()>3)x.latex(l[3]); 

//char_defs[l[1]]=(l[2]);
// lookup either by name of str value...
char_defs[nm]=x; // (l[2]);
char_defs[v]=x; // (l[2]);
//char_defs[l[2]]=(l[1]);
} // grouping 
} // ii

} // load_char_defs

//////////////////////////////////////////////////////////////////////////

void cmd_snacks_time(Cip & cip , LocalVar & lv )
{

const StrTy dest=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy dname=cip.wif(3);
const StrTy pname=cip.wif(4);
const IdxTy ncol=myatoi(cip.wif(5));
const Ragged & data=m_ragged_map[dname];
const Ragged & params=m_ragged_map[pname];
// these are the data file units added wth Muqed but lacking
// in the original hard coded version 
const bool have_units=Bit(flags,0);
// append units soon, any modifiers in the config files need to know the units
const bool append_units=Bit(flags,1);
// add the units later so config file decoratiosn apply to name not units
const bool append_units_after=Bit(flags,2);
const bool parse_glob_file=Bit(flags,3);
const bool enable_grouping=!Bit(flags,4);

if (cip.help()||(dest.length()==0)) {
Ss ss; 
ss<<"  cmd_snacks_time : make latex table for various snack dog food groups"<<CRLF;
ss<< " snacks-rags dest flags dname pname data params "<<CRLF;
ss<<" dest : output file name "<<CRLF;
ss<<" flags : 0, have_units   "<<CRLF;
ss<<" dname : data rags name  "<<CRLF;
ss<<" pname : parameter rags  name "<<CRLF;
ss<<" ncol : if non zero split columns  "<<CRLF;
ss<<" ./mjm_linc_graph.out -cmd \"read-ragged in $SE_DATA_EXCHANGE/snacks_collated.ssv\" -cmd \"read-ragged p happy.txt\" -cmd \"snacks-rags htable 0 in p\" -quit"<<CRLF;

cip.help(ss.str());
}

const IdxTy serial_idx=0;
const IdxTy date_idx=1;
const IdxTy name_idx=2;
const IdxTy fields_idx=3;
const IdxTy nfmax=3;
const bool split_columns=(ncol!=0);
std::ofstream ofs(dest);
std::ostream & os= ofs;// std::cout;
MM_ERR( "cmd_snacks_time "<<MMPR4(dest,flags,dname,pname)<<MMPR3(ncol,data.size(),params.size()))
table_decorators td;
ConfigRules rl; // xxxxx
rl.add_lines(params.begin(), params.end()); //adfasdfadfs
CharDefMap char_defs;
load_char_defs(char_defs,params);

td.decorate(params,data.size(),m_flp);
DDstat dds;
IdxTy lines=0;

// StrTy datefmin="z";
// StrTy datefmax="";
// days in period, date bin, serial, count
//std::map<StrTy , std::map<IdxTy, IdxTy > > dop;
GroupMap gm; 
DopType dop;
//void parse_daily_rag( DopType & dop, DDstat & dds, table_decorators & td
//,MiscForm & mf,  ConfigRules & rl, const Ragged & data ) // const
MiscForm mf;
mf.ncol=ncol;
mf.serial_idx=serial_idx;
mf.date_idx=date_idx;
mf.name_idx=name_idx;
mf.fields_idx=fields_idx;
mf.nfmax=nfmax;
mf.datefmin="z"; //  datefmin;
mf.datefmax=""; //  datefmax;
mf.lines=lines; 
mf.split_columns=split_columns;
mf.have_units= have_units; mf.append_units=  append_units;
mf.append_units_after=  append_units_after; 
if ( parse_glob_file)
{
mf.ditem=1;
++mf.name_idx;
++mf.fields_idx;
mf.group_idx=mf.fields_idx+4;
}
parse_daily_rag(dop,dds,td,gm,mf,rl,data);
#if 0 
} // i
#endif

//write_latex_table(os,td,dds);
Ragged rtable= write_latex_table(td.subdates,td.subnames,td,dds,dop,&gm,  split_columns?1:0, "Name");

const auto & dogs=td.subdogs;
const IdxTy ndogs=dogs.size();
Ss ds; MM_LOOP(ii,dogs)  { ds<<(*ii).first<<" "; } 
Ss cs;
cs<<" Events Summary for "<<ds.str()<<"  from "<<mf.datefmin<<" to "<<mf.datefmax;
//cs<<" From "<<mf.datefmin<<" to "<<mf.datefmax<<" "<<ds.str();
td.tparam["caption"]=cs.str()+td.tparam["caption"];

// noun to class
CharMap char_map;
if (enable_grouping)
{
load_char_map(char_map, char_defs, gm );
} // enable_grouping

MM_ERR(" writing latex table "<<MMPR4(split_columns,cs.str(),rtable.size(),td.tparam.size()))

if (enable_grouping)
{
organize(rtable,td,char_map,char_defs);
} // enable_grouping 
///////////////////////////////////////////////////


if (split_columns)
{
RagSplit rs;
// the this param is the header kluge 
// which is now fixed in the ragged_Table thing
//rs.split_columns(rtable, ncol,1);
rs.split_columns(rtable, ncol+1,0);
td.decorate_latex(os,mf.datefmin,mf.datefmax);
MM_LOOP(ii,rs)
{
// this is cleverly putting headers into tparam split from table
// the top header is a kluge, shifted by one and 
// so is not quite right. Flag 1 should fix headers  
latex_blocks::ragged_table(os,(*ii).second,td.tparam,1);
}
td.decorate_latex(os);
}
else 
{
td.decorate_latex(os,mf.datefmin,mf.datefmax);
latex_blocks::ragged_table(os,rtable,td.tparam);
td.decorate_latex(os);
}




} // cmd_snacks_time
void load_char_map(CharMap& char_map, CharDefMap & char_defs, GroupMap & gm, const IdxTy flags=0 )
{
const bool take_first=Bit(flags,0);
const bool alias_order_check=Bit(flags,1);
const bool highest=Bit(flags,2);
const bool lowest=Bit(flags,3);
MM_LOOP(ii,gm)
{
const StrTy & n=(*ii).first;
IdxTy v=~0;
const auto & x=(*ii).second;
MM_LOOP(jj,x) 
{
MM_ERR("S" << MMPR3(n,(*jj).first,(*jj).second))
const StrTy & x=(*jj).first;
const auto ii=char_defs.find(x);
//if (ii!=char_defs.end()) v=myatoi((*ii).second);
// doh lol 
if (ii!=char_defs.end())
{
 v=(*ii).second; /// myatoi((*ii).second);
if (take_first) break; 
}
//if ((*jj).first=="food") v=1;
//else if ((*jj).first=="amino") v=2;
} 
char_map[n]=v;

} // ii 
} //load_char_map

void sort_perm(PermTy &perm,const Ragged & rtable, CharMap & char_map, TableParams & tparam, const IdxTy xxnhdr )
{
#if 1 
const IdxTy nhdr=myatoi(tparam["headers"]);
MM_ERR(MMPR2(nhdr,xxnhdr))
std::sort(perm.begin(), perm.end(),
[&](const IdxTy& a, const IdxTy& b)
{  
// headers
// these go oor if sort is not consistent apparently 
//MM_ERR(MMPR2(a,b))
//if (a<=nhdr) { return a<b;}
if (a<nhdr) { return a<b;}
//if (b<=nhdr) { return a>b;}
if (b<nhdr) { return !true;}
//const StrTy & k1=rtable[perm[a]][0]; 
const StrTy & k1=rtable[a][0]; 
//const StrTy & k2=rtable[perm[b]][0];
const StrTy & k2=rtable[b][0];
const IdxTy c1=char_map[k1];
const IdxTy c2=char_map[k2];
if (c1==c2) return (k1<k2); return (c1<c2);  });
#endif

} // sort_perm


void organize(Ragged& rtable, table_decorators & td, CharMap & char_map, CharDefMap & char_defs)
{
// rtable can be sorted if the headers tracked ok... 
/////////////////////////////////////////////////
const IdxTy tlines=td.subnames.size(); // names.size();
const IdxTy twords=td.subdates.size(); // dogs.size();
const IdxTy nhdr=myatoi(td.tparam["headers"]);
Ragged rtable2(tlines+1,twords+1);
PermTy  perm(tlines+1);
MM_ERR(" sorting "<<MMPR4(tlines,perm.size(),rtable2.size(),rtable.size()))
{ MM_SZ_LOOP(i,perm,sz) { perm[i]=i; }  } 
Ss ss; MM_LOOP(ii,perm) { ss<< (*ii)<<" ";} MM_ERR(MMPR(ss.str()))
//std::sort(perm.begin(), perm.end(), &Myt::ASUCK);
//sort_perm(PermTy &perm,const Ragged & rtable, CharMap & char_map );
sort_perm(perm,rtable, char_map,td.tparam,nhdr );
// doh, this is about the same as number of nouns lol.... 
//IdxTy chars=char_map.size();
IdxTy chars=char_defs.size();
MM_ERR(" sort done") 
//MM_SZ_LOOP(i,perm,sz) { rtable2[i]=rtable[perm[i]]; } // i 
MM_SZ_LOOP(i,perm,sz) { rtable2.set(rtable[perm[i]],i); } // i 
MM_ERR(" ser done") 
IdxTy pchar=~1; // ~0 is taken as is zero.... 
rtable=Ragged(tlines+chars+1,twords+1);
IdxTy j=0;
Ragged::Line hdr(twords+1);
MM_SZ_LOOP(i,rtable2,szr)
{
const auto& line=rtable2[i];
const IdxTy ch=char_map[line[0]];
MM_ERR(" chars "<<MMPR4(i,pchar,ch,line[0]))
if ((ch!=pchar)&&(i>=nhdr))
{
Ss ss; ss<<ch;
//hdr[0]=StrTy("{\\bf ")+StrTy(char_defs[ss.str()])+"}"; // " new group ";
hdr[0]=(char_defs[ss.str()].latex()); // " new group ";
rtable.set(hdr,j); ++j;
pchar=ch;
} // diff 
//if (line[0]=="") 
MM_ERR(" maybe blank input line "<<MMPR4(j,i,line.size(),line[0])<<MMPR3(tlines,chars,rtable2.size()))
rtable.set(line,j); ++j;
}
while (j<(tlines+chars)) { hdr[0]=""; rtable.set(hdr,j);  ++j;}
//rtable=rtable2;
}

#if 0 
static bool ASUCK(const IdxTy& a, const IdxTy& b)
{  
// headers
MM_ERR(MMPR2(a,b))
if (a<=1) { return a<b;}
if (b<=1) { return a>b;}
const StrTy & k1=rtable[perm[a]][0]; 
const StrTy & k2=rtable[perm[b]][0];
const IdxTy c1=char_map[k1];
const IdxTy c2=char_map[k2];
if (c1==c2) return (k1>k2); 
return (c1>c2); 
}
#endif




bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }

template <class Tp,class Tdop>
Ragged write_ragged_table
(Tp & subx, Tp & suby, table_decorators & td,   DDstat &  dds, Tdop & dop, const IdxTy flags=0,
const StrTy & name=""  )
{
Ragged rtable;
IdxTy iname=0;
MM_LOOP(ii,suby)
{
StrTy  name=(*ii).first;
StrTy  nm=name;
//td.decorate_y(nm,name,suby.key_name(name));
IdxTy idog=0;
Ragged::Line x;
//rtable.ncline(iname)[idog]=name;
x.push_back(name);
++idog;
MM_LOOP(jj,subx)
{

const StrTy & dog=(*jj).first;
// cleverly this went from drug to dog "name" lol. 
Stat &   q= dds[dog][nm];
//q.calculate((*jj).second);
q.calculate(dop[dog].size());
//StrTy xxx=q.to_string();
//rtable.ncline(iname)[idog]=x;
const IdxTy fields=q.fields();
for(IdxTy i=0; i<fields; ++i)  x.push_back(q.to_ssv(0,i));
++idog;
} // jj 

++iname;
} // ii 


return rtable;
} // write_ragged_table

//void write_latex_table(OsTy & os, table_decorators &td,  DDstat &  dds)
template <class Tp,class Tdop>
Ragged write_latex_table
(Tp & subx, Tp & suby, table_decorators & td,   DDstat &  dds, Tdop & dop, GroupMap * gm, const IdxTy flags=0,
const StrTy & name=""  )
{

const bool fix_header_shift=Bit(flags,0); 
// parsing phase 
const IdxTy tlines=suby.size(); // names.size();
const IdxTy twords=subx.size(); // dogs.size();
Ragged rtable(tlines+1,twords+1);

IdxTy iname=0;
IdxTy idog=0;
if (fix_header_shift) { rtable.ncline(iname)[idog]=name; ++idog; } 
//rtable.ncline(iname)[idog]="Ingredient";
MM_LOOP(jj,subx)
{rtable.ncline(iname)[idog]=(*jj).first;
++idog;
}
++iname;
MM_LOOP(ii,suby)
{
StrTy  name=(*ii).first;
StrTy  nm=name;

//td.decorate_y(nm,name);
// add flat 1 to always enclose in math mode for sort order
td.decorate_y(nm,name,suby.key_name(name),1);
if ( gm ) (*gm)[name]=(*gm)[nm];
IdxTy idog=0;
rtable.ncline(iname)[idog]=name;
++idog;
MM_LOOP(jj,subx)
{

const StrTy & dog=(*jj).first;
// cleverly this went from drug to dog "name" lol. 
Stat &   q= dds[dog][nm];
//q.calculate((*jj).second);
q.calculate(dop[dog].size());
StrTy x=q.to_string();
rtable.ncline(iname)[idog]=x;

++idog;
} // jj 

++iname;
} // ii 
/*Ss cs;
//cs<<" From "<<datefmin<<" to "<<datefmax<<" ";
td.tparam["caption"]=cs.str()+td.tparam["caption"];

MM_ERR(" writing latex table "<<MMPR3(cs.str(),rtable.size(),td.tparam.size()))

//td.decorate_latex(os,datefmin,datefmax);
latex_blocks::ragged_table(os,rtable,td.tparam);
td.decorate_latex(os);
*/

return rtable;
} // write_latex_table


////////////////////////////////////////////////////////////////////
void cmd_zymo_rags(Cip & cip , LocalVar & lv )
{

const StrTy ragin=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const Ragged & qio=m_ragged_map[ragin];
const Ragged & params=m_ragged_map[cmd2];
const IdxTy nval=myatoi(cmd1);
const IdxTy maxcnt=m_flp.maxcnt();
const bool output_latex=(((1<<0)&flags)!=0);
const bool  output_rank_table=(((1<<1)&flags)!=0);
const bool  output_ranks=(((1<<2)&flags)!=0)&&!output_rank_table;
const bool output_summary=(((1<<3)&flags)!=0);

const bool output_dump=!output_latex&&!output_rank_table&&!output_ranks
&&!output_summary;
MM_ERR( "cmd_zymo_rags "<<MMPR4(nval,maxcnt,ragin,cmd2)<<MMPR4(flags,output_latex,output_ranks,output_summary))
std::ostream & os= std::cout;
std::vector<StrTy> val_names(nval);
const IdxTy psz=params.size();
for(IdxTy i=0; i<nval; ++i)
{
const bool have=(i<psz);
val_names[i]=have?params.line(i)[0]:StrTy("x");

}

MyBlock vals;
std::vector<StrTy> names;
const IdxTy sz=qio.size();
vals.resize(sz,nval);
std::vector<D> merit(sz);
std::vector<IdxTy> idxs(sz);
std::vector<IdxTy> idxs2(sz);
//IdxTy ranks[sz*nval];
std::vector<std::map<IdxTy, IdxTy> > rank_maps(nval);

for(IdxTy i=0; i<sz; ++i)
{
idxs[i]=i;
const Ragged::Line & l=qio.line(i);
const IdxTy lsz=l.size();
//MM_ERR(MMPR2(i,lsz))
for(int j=0; j<int(nval); ++j)
{
int diff=nval-j-1;
if (diff>=int(lsz)) continue;
//MM_ERR(MMPR4(i,lsz,diff,l[lsz-diff-1]))
const D x=atof(l[lsz-diff-1].c_str());
vals(i,diff)=x;
merit[i]+=x;
}

// ) ; // sort
StrTy name="";
IdxTy jmax=lsz;
if (lsz>nval) jmax=lsz-nval;
//MM_ERR(MMPR3(i,lsz,jmax))
for(IdxTy j=0; j<jmax; ++j)
{
const StrTy & w=l[j];

 name=name+w;
}

names.push_back(name);
} // i 
if (!false) { 
std::sort(idxs.begin(), idxs.end(),
//std::sort(xord.begin(),xord.end(),
[&](const IdxTy& a, const IdxTy& b)
{ return merit[a]>merit[b]; });

for(IdxTy j=0; j<(nval); ++j)
{
for(IdxTy i=0; i<sz; ++i) {  idxs2[i]=i; } 
std::sort(idxs2.begin(), idxs2.end(),
//std::sort(xord.begin(),xord.end(),
[&](const IdxTy& a, const IdxTy& b)
{ return vals(a,j)>vals(b,j); });
auto & m=rank_maps[j];
if (output_rank_table){ 
for(IdxTy i=0; i<sz; ++i) {  m[idxs2[i]]=i; }  }
if (output_ranks){ 
for(IdxTy i=0; i<sz; ++i) {  m[i]=idxs2[i]; }  }
} // j 
}  // false
if ( output_summary)
{
IdxTy counts[nval];
D H[nval];
D norm[nval];
const D scale=1.0/log(2);
latex_blocks::table_header(os, nval, val_names,0);
for(IdxTy q=0; q<(nval); ++q)
{
counts[q]=0;
H[q]=0;
norm[q]=0;
for(IdxTy _i=0; _i<sz; ++_i) {  
const IdxTy i=_i; // rank_maps[q][_i];
const D & v=vals(i,q);
norm[q]+=v;
}
if (norm[q]==0) {  MM_ERR(" danger will robinson zero norm "<<MMPR3(norm[q],sz,q))}
const D nv=1.0/norm[q]; // 
for(IdxTy _i=0; _i<sz; ++_i) {  
const IdxTy i=_i; // rank_maps[q][_i];
const D & v=vals(i,q)*nv;
const bool count= (v!=0);
if (count)
{
++counts[q];
H[q]-=v*log(v)*scale;

}// count

} // _i

} // q

latex_blocks::table_line(os,nval,"counts",&counts[0]);
latex_blocks::table_line(os,nval,"H",&H[0]);

/*
os<<"counts";
for(IdxTy q=0; q<nval; ++q)
{
os<<"&"<<counts[q];
}
os<<"\\\\"<<CRLF;
*/


latex_blocks::table_end(os,0);
return;
} // output_summary

if (output_ranks)
{
for(IdxTy q=0; q<(nval); ++q)
{

os<<"\\begin{table}[H] \\centering"<<CRLF;
os<<"\\begin{tabular}{|";
 os<<"l|";
for(IdxTy p=1; p<=(nval); ++p) os<<"r|";
os<<"}"<<CRLF;
os<<"\\hline"<<CRLF;
os<<"Name";
for(IdxTy p=0; p<(nval); ++p) os<<"&"<<val_names[p];
os<<"\\\\"<<CRLF;
os<<"\\hline"<<CRLF;

for(IdxTy _i=0; _i<sz; ++_i) {  
if ( _i>=maxcnt) break; 
const IdxTy i=rank_maps[q][_i];
Ss ss;
ss<<fixslash(shorter(names[i]));
//for(IdxTy l=0; l<nval; ++l) { ss<<"&"<<std::setprecision(3)<<vals(i,l); } 

for(IdxTy l=0; l<nval; ++l) { 
Ss sv;
if ( vals(i,l)!=0) sv<<std::setprecision(3)<<(vals(i,l));

ss<<"&"<<sv.str(); } //  (rank_maps[l][i]+1);} 
os<<ss.str()<<"\\\\"<<CRLF;
} // i  

os<<"\\hline"<<CRLF;
os<<"\\end{tabular}"<<CRLF;
os<<"\\end{table}"<<CRLF;

} // q 
MM_ERR(" done output ranks ")
return; 
} // output_ranks


if ( output_rank_table)
{
os<<"\\begin{table}[H] \\centering"<<CRLF;
os<<"\\begin{tabular}{|";
 os<<"l|";
for(IdxTy p=1; p<=(nval); ++p) os<<"r|";
os<<"}"<<CRLF;
os<<"\\hline"<<CRLF;
os<<"Name";
for(IdxTy p=0; p<(nval); ++p) os<<"&"<<val_names[p];
os<<"\\\\"<<CRLF;
os<<"\\hline"<<CRLF;
for(IdxTy _i=0; _i<sz; ++_i) {  
if ( _i>=maxcnt) break; 
const IdxTy i=idxs[_i];
Ss ss;
ss<<fixslash(shorter(names[i]));
//for(IdxTy l=0; l<nval; ++l) { ss<<"&"<<std::setprecision(3)<<vals(i,l); } 

for(IdxTy l=0; l<nval; ++l) { 
Ss sv;
if ( vals(i,l)!=0) sv<<(rank_maps[l][i]+1);

ss<<"&"<<sv.str(); } //  (rank_maps[l][i]+1);} 
os<<ss.str()<<"\\\\"<<CRLF;
} // i  

os<<"\\hline"<<CRLF;
os<<"\\end{tabular}"<<CRLF;
os<<"\\end{table}"<<CRLF;

return; 
} // output_rank_table


MM_ERR(" tabulated "<<MMPR(sz))

IdxTy combmax=(1<<nval);
for(IdxTy i=0; i<combmax; ++i)
{
Ss sl;
IdxTy cnt=0;
for(IdxTy _j=0; _j<sz; ++_j)
{
const IdxTy j=idxs[_j];
//MM_ERR(MMPR2(j,_j))
IdxTy k=0;
for (; k<nval; ++k)
{
const bool bit=(((1<<k)&i)!=0);
if ( bit ) if (vals(j,k)==0) break;
if ( !bit ) if (vals(j,k)!=0) break;
} // k 
if ( k==nval)
{
 ++cnt;
if ( cnt<maxcnt) { 
Ss ss;
if (output_dump)  {
for(IdxTy l=0; l<nval; ++l) { ss<<vals(j,l)<<" "; } 
 MM_ERR(MMPR2(ss.str(),(shorter(names[j])))) } 
if (output_latex) { 
for(IdxTy l=0; l<nval; ++l) { ss<<"&"<<std::setprecision(3)<<vals(j,l); } 
//MM_ERR(MMPR2(ss.str(),shorter(names[j]))) 
sl<<fixslash(shorter(names[j]))<<""<<ss.str()<<"\\\\"<<CRLF;
} //output_latex 

}}
} // j 

Ss ss;
if (output_dump) { 
for(IdxTy l=0; l<nval; ++l) { 
const bool bit=(((1<<l)&i)!=0);
if (bit) ss<<val_names[l]<<" ";
} // l 
ss<< std::hex<<i<<" "<<std::dec<<cnt<<CRLF;
MM_ERR(ss.str()) } 
if (output_latex) { 
const IdxTy padctc=nval+1-3;
IdxTy padct=padctc;
StrTy pad="";
while (padct!=0){ pad+="&"; --padct; } 
for(IdxTy l=0; l<nval; ++l) { 
const bool bit=(((1<<l)&i)!=0);
if (bit) ss<<val_names[l]<<" ";
} // l 
ss<<"&"<< std::hex<<i<<"& "<<std::dec<<cnt<<pad<<"\\\\"<<CRLF;
os<<"\\begin{table}[H] \\centering"<<CRLF;
os<<"\\begin{tabular}{|";
 os<<"l|";
for(IdxTy p=1; p<=(nval); ++p) os<<"r|";
os<<"}"<<CRLF;
os<<"\\hline"<<CRLF;
os<<ss.str(); // <<CRLF;
os<<"\\hline"<<CRLF;
os<<"Name";
for(IdxTy p=0; p<(nval); ++p) os<<"&"<<val_names[p];
os<<"\\\\"<<CRLF;
os<<"\\hline"<<CRLF;
os<<sl.str();
os<<"\\hline"<<CRLF;
os<<"\\end{tabular}"<<CRLF;
os<<"\\end{table}"<<CRLF;


}
 

} // i 
} // zymo_rags
StrTy fixslash(const StrTy & s)
{
const char * p =s.c_str();
const IdxTy sz=s.length();
char d[sz+1];
d[sz]=0;
for(IdxTy i=0; i<sz; ++i)
{
char c=p[i];
if (c=='_') c='.';
d[i]=c;
}

return StrTy(d);

}

StrTy shorter(const StrTy & s)
{
const IdxTy sz=s.length();
if (sz<4) return s;
const char * posi=s.c_str();
const char * p=posi;
while ( *(p+2) !=0)
{
if ( *p=='g') if (*(p+1)=='_') if (*(p+2)=='_') 
{
const char c1=*(p+3);
if (c1==0) return  s;
const char c2=*(p+4);
if (c2==0) return  s;
if ((c1=='N')&&(c2=='A'))
{
IdxTy b=0;
IdxTy pos=p-s.c_str();
IdxTy bc=0;
for(IdxTy jx=0; jx<pos; ++jx)
{
const char c1=*(posi+pos-jx+1);
const char c2=*(posi+pos-jx);
if ((c1=='_')&&(c2=='_') ) 
{
++bc;
if ( bc>1) return s.substr(pos-jx-1,jx);

}


}
 return s.substr(b,pos-b); 
}

return StrTy(p); 
}


++p;
}

return s; 
}
typedef mjm_partition_itor PartItor;
void cmd_linc_graph(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
//const StrTy sfasta=cip.p2; // wif(3);
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy srag=cip.wif(3);
const StrTy qiime=cip.wif(4);
// was going to take a raggeed input 
Dig & sr=m_dig_map[ragin];
//const Ragged & qiiq=m_ragged_map[qiime];
sr.subs(flags);

} // cmd_linc_graph



typedef mjm_cursor_pair_itor CpiTy; // (m_i, l,  highest);
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
typedef std::vector<D> Dpv;
typedef std::vector< Dpv > Dp;

#ifdef HAVE_OLD_PROBS

#endif // #ifdef HAVE_OLD_PROBS

void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
{
Canned::cmd_add_to_fasta(cip ,  lv,  m_fasta_map  );
#if 0 
const StrTy name=cip.p1;
const StrTy seq=cip.p2;
Fasta & f=m_fasta_map[name];
StrTy sn=cip.wif(3);
if (sn.length()==0)
{ Ss ss ; ss<<"seq"<<(f.size()+1); sn=ss.str(); }
sn=StrTy(">") + sn;
f.add(sn,seq); // load(fn);
MM_ERR("adding to  "<<MMPR4(name,m_fasta_map[name].size(),sn,seq))
#endif
}

void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
{
Canned::cmd_zymo_merge_fasta(cip ,  lv,  m_fasta_map, m_ragged_map  );
#if 0
//const StrTy name=cip.p1;
const StrTy dfasta=cip.p1;
// file with sequence numbers only 
const StrTy sfasta=cip.p2; // wif(3);
// annotation in zymo format
const StrTy srag=cip.wif(3);
const StrTy qiime=cip.wif(4);
//CommandInterpretter li;
Fasta & df=m_fasta_map[dfasta];
const Fasta & sf=m_fasta_map[sfasta];
const Ragged & sr=m_ragged_map[srag];
const Ragged & qiiq=m_ragged_map[qiime];
// xlate_map xlate_field_map(const IdxTy key, const IdxTy val)
Ragged::xlate_double_map qii_map=qiiq.xlate_field_max_map(8,3);
MM_ERR(" SHTCK "<<MMPR(qii_map.size()))
std::map<StrTy, StrTy> kv;
const IdxTy szsr=sr.size();
for (IdxTy i=0; i<szsr; ++i)
{
//MM_ERR(" ASUCK "<<MMPR2(i,sr.line(i).size()))
const IdxTy fick=sr.line(i).size();
if (fick!=1)
{
MM_ERR(" the ragged file needs to be loaded with flag bit 0 set (1) "<<MMPR3(i,fick,srag))
return;
}
const char * c=sr.line(i)[0].c_str();
StrTy k="";
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; k=StrTy(d);
while (c[j]!=0) { if ((c[j]!='\t') &&(c[j]!=' ')) break; ++j; }
kv[k]=StrTy(c+j);

}

MM_ERR("zumo merge fasta fn into name "<<MMPR3(dfasta,sfasta,srag)<<MMPR3(df.size(),sf.size(),sr.size()))
const IdxTy sz=sf.size();
StrTy sep=" ";
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & seq=sf.seq(i);
StrTy sn=sf.name(i);
StrTy sqname="";
const char * c=sn.c_str();
if (c[0]!=0) { 
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; 
StrTy k=StrTy(d+1);

StrTy annot=kv[k];
const char * shtfick=sn.c_str();
sqname=StrTy(shtfick+1); // sn; // k;
//MM_ERR(MMPR2(k,annot))
StrTy sfx=""; 
auto im=qii_map.find(sqname);
//MM_ERR(" looking up "<<MMPR(sqname))
if (im!=qii_map.end()) { Ss sf; sf<<" qzscore="; sf<<(*im).second; sfx=sf.str(); } 
else { 
// MM_ERR(" ASSFUCL "<<MMPR2(sqname,(*(qii_map.begin())).first))

 } 
sn=sn+sep+zymo_condense(annot)+sfx;
}
df.add(sn,seq); // load(fn);
}

#endif


}

#if 0 
StrTy zymo_condense(const StrTy & s ) const
{
const char * c=s.c_str();
const IdxTy sz=s.length();
//IdxTy start=0;
StrTy genus="";
StrTy species="";
if (sz<3) return StrTy();
for (IdxTy i=0; i<(sz-3); ++i)
{
if ((c[i]=='g')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
//  buserror? wtf memcpy up? 
// the bus error got me, I think the length j is just wrong doh... 
char d[j+1-i]; //memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
genus=StrTy(d+3);
continue;
}
if ((c[i]=='s')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
// memcpy probably ok just had wrong length.. 
char d[j+1-i]; // memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
species=StrTy(d+3);
continue;
}



} // i 
return genus+StrTy(" ")+species;
} 

#endif

//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
void cmd_write_svg_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const StrTy prag=(cip.wif(4));

Ragged & r=m_ragged_map[name];
Ragged & pr=m_ragged_map[prag];

ReadWriteMap rwm;
pr.to_map(rwm);
mjm_svg_writer sw;
std::ofstream  os(fn);

IdxTy text_idx=0;rwm.get("text_idx",text_idx);
IdxTy jmax=3;rwm.get("jmax",jmax);
IdxTy jfirst=1;rwm.get("jfirst",jfirst);

const IdxTy min_line=jfirst+jmax;
r.min_line_size(min_line);
const IdxTy cats=r.size();
IdxTy tm=50+30; rwm.get("tm",tm);
IdxTy bm=330;rwm.get("bm",bm);
IdxTy lm=200;rwm.get("lm",lm);
IdxTy rm=50;rwm.get("rm",rm);
D catpitch=30;rwm.get("catpitch",catpitch);
D ysz=600;rwm.get("ysz",ysz);
D ynames=.5*D(tm); rwm.get("ynames",ynames);
D namesz=10; rwm.get("namesz",namesz);

//D ymin=tm;
//D ymax=ymin+ysz;
D yleg=tm+ysz+30;
D szleg=catpitch*.9;
StrTy colorleg="red";
IdxTy xs=IdxTy( catpitch*cats+.5)+lm+rm;
IdxTy ys=ysz+tm+bm;
D wr=catpitch;
D hr=.01*ysz;
D xmin=1e-4;
D vmin=4;
D vmax=vmin-1;

IdxTy ygrid=9;
IdxTy xz=lm;
IdxTy xf=lm+cats*catpitch;
IdxTy yz=tm;
IdxTy yf=tm+ysz;
D xleg=.9*lm;
StrTy ylab= StrTy(" Range limited log ratio to mean");

std::vector<StrTy> colors(jmax);
std::vector<StrTy> names(jmax);
MM_SZ_LOOP(i,names,xxx)
{
Ss ss,sc;
StrTy name="name";
StrTy col="color";
ss<<name<<i;
sc<<col<<i;
name=ss.str();
//MM_ERR(MMPR2(rwm.size(),name))
rwm.get(name,name);
names[i]=name;
StrTy c="";
rwm.get(sc.str(),c);
colors[i]=c;
} // names 

MM_ERR(" cmd_write_svg_ragged to "<<MMPR4(cmd,fn,name,flags)<<MMPR3(r.size(),cats,prag))


//ow.bounds(x,y,lc.m_level_sz);
os<<sw.start_text(" test foo",xs,ys);
//os<<sw.frame_text("#00808080",xs,ys);
os<<sw.frame_text("#FFFFFF",xs,ys);
os<<CRLF;

std::vector<D> merit(cats);
std::vector<D> avg(cats);
for(IdxTy i=0; i<cats; ++i)
{
D sum=0; 
D n=0;
D mmax=0;
D mmin=0;
const Ragged::Line line=r.line(i);
if (line.size()>=4) { 
for (IdxTy j=0; j<jmax; ++j)
{
//MM_ERR("FUR RRRR  "<<MMPR3(i,j,line.size()))
const StrTy & vs= line[jfirst+j];
const D dv=atof(vs.c_str());
if (dv==0) continue;
if (n==0) { mmax=dv; mmin=dv; } 
if (dv>mmax) mmax=dv;
if (dv<mmin) mmin=dv;
sum+=dv;
n+=1;
} //j 
} // if 
if ( n!=0) avg[i]=sum/n; else avg[i]=0;
if ( n==jmax) merit[i]=mmax/mmin;
else merit[i]=-mmax;
} //i 

//std::sort(idxs.begin(), idxs.end(),
//[&](const IdxTy& a, const IdxTy& b)
//{ return merit[a]>merit[b]; });

Ragged::sort_order so=r.ordering();
r.sort(so,[&](const IdxTy& a, const IdxTy& b) -> bool
{ 
if (merit[a]<0) 
{
if (merit[b]<0) return merit[a]<merit[b]; 
return !false;
}
if (merit[b]<0) 
{
//if (merit[a]<0) return merit[a]>merit[b]; 
return false;
}


return merit[a]>merit[b]; 

});


for (IdxTy j=0; j<jmax; ++j)
{
const IdxTy m=(j%3);
StrTy c="black";
if (m==0) c="red";
else if (m==1) c="green";
else if (m==2) c="blue";

if (colors[j].length()==0) colors[j]=c;
}

for(IdxTy pass=0; pass<2; ++pass)
{
for (IdxTy j=0; j<jmax; ++j)
{
const StrTy & cj=colors[j];
for(IdxTy i=0; i<cats; ++i)
{
const IdxTy is=so[i];
const D x=lm+i*catpitch;

const Ragged::Line line=r.line(is);
//MM_ERR("  RRRR  "<<MMPR3(i,j,line.size()))
if (line.size()<=(jfirst+j)) continue; 
const StrTy & vs= line[jfirst+j];
const D vraw=atof(vs.c_str());
D v=xmin;
if (vraw!=0) v=log(vraw/avg[is])/log(10);
//else  v=log(v/avg[i])/log(10);
if (pass==0)
{
if (vraw!=0) { 
if ( vmax<vmin) { vmax=v; vmin=v; }
else {
if ( v>vmax) vmax=v;
if (v<vmin) vmin=v;
} // max<min
 } // vraw
 
} // pass==0
if (pass==1) { 
v=(v-vmin)/(vmax-vmin);
// v is between zero and 1 
const D y=tm+(1.0-v)*ysz;
// TODO overlaps ???? 
//if ( vraw!=0) os<<sw.shade_rect_text(x,y,wr,hr,cj,.5);
if ( vraw!=0) os<<sw.shade_rect_text(x,y,wr,hr,cj,.8);
else
{
const D y=tm+(1.0-0)*ysz;
 //os<<sw.shade_rect_text(x,y,wr,hr,cj,.5);
 os<<sw.shade_rect_text(x,y,wr,hr,cj,.8);
os<<sw.gtext_text(StrTy("*"),x+(D(j+1)/(jmax+2))*wr , y, szleg, cj,  StrTy("start") ,0); 

}
os<<CRLF;
} // pass
} //i  
} // j 
if (pass==0)
{
Ss ss,sr;
ss<<std::setprecision(2)<<(vmax*1.1);
ss>>vmax;
// this is negative need to make it  more negative  
sr<<std::setprecision(2)<<(vmin*1.1);
sr>>vmin;


}

} // pass



os<<sw.gtext_text(ylab, szleg+10, yf, szleg, colorleg,  StrTy("end") ,90); 
//StrTy hstripes(const D & xz,const D & xf, const D & yz, const D & yf,const IdxTy & n, const StrTy& col, const D & thick)
os<<sw.hstripes( xz, xf, yz,  yf, ygrid, StrTy("black"), 1);
os<<sw.vstripes( xz, xf, yz,  yf, cats, StrTy("black"), 1);


for(IdxTy i=0; i<=ygrid; ++i)
{
const D f=D(i)/ygrid;
const D y=tm+ysz*(1.0-f);
Ss ss;
ss<<std::setprecision(3)<<((vmin+(vmax-vmin)*f));
const StrTy & text=ss.str();
// this  does not anchor on  right or  levt 
os<<sw.gtext_text(text, xleg, y, szleg, colorleg,  StrTy("end") ,0); 
os<<CRLF;

}

for(IdxTy i=0; i<cats; ++i)
{
const IdxTy is=so[i];
const Ragged::Line line=r.line(is);
if (line.size()<=text_idx) continue; //  
const StrTy & text= line[text_idx];
const D x=lm+i*catpitch+szleg;
os<<sw.gtext_text(text, x, yleg, szleg, colorleg,  StrTy("end") ,-45); 
os<<CRLF;
} // i 

{
const D xnames=lm+10;
D x=xnames;
for (IdxTy j=0; j<jmax; ++j)
{
const StrTy & cj=colors[j];
const StrTy & text=names[j];
const D y=ynames;
os<<sw.shade_rect_text(x,y,wr,hr,cj,.8);
os<<sw.gtext_text(text, x+wr, y, namesz, cj,  StrTy("start") ,0); 
os<<CRLF;



 x+=+wr+10+text.length()*namesz; 
} // j
} // names scoping


os<<sw.end_text();
os<<CRLF;

} // write_svg_ragged
void cmd_transpose_ragged(Cip & cip , LocalVar & lv ) 
{
Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ;
#if 0 
const StrTy cmd=cip.cmd();
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
Ragged & r=m_ragged_map[sname];
Ragged & d=m_ragged_map[dname];
MM_ERR(" cmd_transpose_ragged "<<MMPR3(dname,sname,flags)<<MMPR2(r.size(),d.size()))
r.transpose(d);
MM_ERR(" cmd_transpose_ragged "<<MMPR3(dname,sname,flags)<<MMPR2(r.size(),d.size()))
#endif
} // transpose

void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{
Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ;
#if 0 
const StrTy cmd=cip.cmd();
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const bool load_string=(cmd=="string-ragged");
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
const bool csv=((flags&16)!=0);
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_read_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR2(debug,csv)<<MMPR4(load_string, load_parsed,tabsep,ignore_hash))
Ragged & r=m_ragged_map[name];
r.sep(" ");
if (tabsep) r.sep("\t");
if (csv) r.sep(",");
if (ignore_hash) r.ignore_hash(ignore_hash);
if (load_string) r.load_from_string(fn,debug);
else 
if (load_parsed) r.load(fn,debug);
else r.load_lines(fn,debug);

MM_ERR(MMPR2(r.size(),name))
#endif

}


void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{
Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ;
#if 0
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const bool load_string=(cmd=="string-ragged");
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
const bool csv=((flags&16)!=0);
const IdxTy wflags=0;
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_write_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR2(debug,csv)<<MMPR4(load_string, load_parsed,tabsep,ignore_hash)<<MMPR(wflags))
Ragged & r=m_ragged_map[name];
r.sep(" ");
if (tabsep) r.sep("\t");
if (csv) r.sep(",");
if (ignore_hash) r.ignore_hash(ignore_hash);
//if (load_string) r.load_from_string(fn,debug);
//else 
//if (load_parsed) r.load(fn,debug);
//else r.load_lines(fn,debug);
r.write_file(fn,wflags);

MM_ERR(MMPR2(r.size(),name))
#endif

}



void cmd_app_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy dest=cip.p1;
const StrTy src=cip.p2;
//const StrTy src=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(3));
if ((cip.help()||(src.c_str()[0]==0)))
{
Ss ss;
ss<<" cmd_app_ragged dest src flags "<<CRLF;
ss<< " append the ragged table src to dest "<<CRLF;
ss<< " dest : destination ragged name t"<<CRLF;
ss<< " src : sourc ragged name t"<<CRLF;
ss<< " flags : none  "<<CRLF;
cip.help(ss.str());
return; 
}
Ragged & d= m_ragged_map[dest];
const Ragged & s= m_ragged_map[src];
MM_ERR(" appending raggeds "<<MMPR3(dest,src,flags))
d+=s;

}

void cmd_add_ragged(Cip & cip , LocalVar & lv ) 
{
Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ;
#if 0 
const StrTy cmd=cip.cmd();
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
//const bool load_string=(cmd=="string-ragged");
//const bool load_parsed=((flags&1)==0);
//const bool tabsep=((flags&2)!=0);
//const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_add_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR(debug))
Ragged & r=m_ragged_map[name];
r.add(cip.words2());

MM_ERR(MMPR2(r.size(),name))
#endif

}



void cmd_read_dig(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
//const bool load_parsed=((flags&1)==0);
//const bool tabsep=((flags&2)!=0);
MM_ERR(" cmd_read_dig from "<<fn<<" to "<<name<<" flags "<<flags)
 m_dig_map[name].load(fn);
//if (tabsep) m_dig_map[name].sep("\t");
//if (load_parsed) m_dig_map[name].load(fn);
//else m_dig_map[name].load_lines(fn);

MM_ERR(MMPR2(m_dig_map[name].size(),name))
}






void cmd_query_aln(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy query=cip.p2;
const StrTy field=cip.wif(3);
//const Aln & pro=m_aln_map[name];
MM_ERR(" query_protrait "<<MMPR3(name,query,field))
//auto w=pro.query(query,field);
//MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR4(name, query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))

}

void cmd_tt(Cip & cip , LocalVar & lv ) 
{
Canned::cmd_tt( cip ,  lv, m_tax_trees ) ;
#if 0 
const StrTy cmd=cip.p1;
const StrTy name=cip.p2;
const StrTy p1=cip.wif(3);
TaxTree & tt=m_tax_trees[name];
MM_ERR(" cmd_tt "<<MMPR3(cmd,name,p1))
do { 
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }
} while (false); 

#endif

}
#if 0
// this needs a rational version 
D binom(const IdxTy top, const IdxTy bot,const D * logtable)
{
D fck= exp(logtable[top]-logtable[bot]-logtable[(top-bot)]);
//MM_ERR(MMPR4(fck,top,bot,exp(logtable[top])))
return fck;
}

#ifdef USE_RATIONALS
typedef mjm_rational RatTy;

#if DUPLICATE_RAT_CODE

std::vector<RatTy> m_fact_table;
RatTy binom(const IdxTy top, const IdxTy bot)
{
const bool use_table=true;
if (use_table)
{
if ( m_fact_table.size()<=top)
{
if (m_fact_table.size()==0) m_fact_table.push_back(1);
while ( m_fact_table.size()<=top) 
{ m_fact_table.push_back(m_fact_table.back()*IdxTy(m_fact_table.size())); } 
}
// this may be slower than the multiplication lol 
return m_fact_table[top]/(m_fact_table[bot]*m_fact_table[top-bot]); 
}
RatTy c=1;
for (IdxTy i=top; i>(top-bot); --i) c=c*(i);
for (IdxTy i=2; i<=bot; ++i) c=c/(i);
//MM_ERR(MMPR3(top,bot,c))
return c; 
}

typedef std::vector<RatTy> RatCo;
RatCo sum_coefs(const IdxTy n, const IdxTy flags=0)
{
//std::vector<RatTy> 
RatCo coefs;
for (IdxTy j=0; j<=n; ++j)
{
RatTy sum=0;
for (IdxTy i=n+2-j; i<=(n+1); ++i)
{
sum+=coefs[n+1-i]*binom(i,n-j);

} // i 
RatTy an1j=(binom(n,j)-sum)/binom(n+1-j,n-j);
coefs.push_back(an1j);
} // j 
if ((flags&1)==0) { 
RatTy total=0;
MM_SZ_LOOP(i,coefs,sz)
{
const IdxTy p=sz-i;
const RatTy & x=coefs[i];
total+=x; // coefs[i];
MM_MSG(MMPR4(i,x,p,total)<<MMPR(D(x)))
} }
return coefs;
} // sum_coefs

StrTy  format_poly(const IdxTy i, const IdxTy n, const RatCo & x,  const IdxTy style)
{
Ss ss;
if (style==0) {
ss<<i;
for(IdxTy j=0; j<x.size(); ++j)
{
ss<<"&"<<x[j];
}
for(IdxTy j=x.size(); j<=n; ++j) ss<<"&";

ss<<"\\\\";
} // 0
else if (style==1)
{
ss<<"y"<<i<<"=";
for(IdxTy j=0; j<x.size(); ++j)
{
ss<<"+("<<x[j]<<")*x^("<<(i+1-j)<<")";
}
ss<<";"<<CRLF;

} // 1

return ss.str();
}
void nesting(RatCo & d)
{
RatCo  y;
for (IdxTy i=0; i<d.size(); ++i)
{
RatCo x=sum_coefs(d.size()-i,1);
while ( x.size()>y.size() ) y.push_back(0);
for (IdxTy j=0; j<x.size(); ++j)
{
const IdxTy xpow=x.size()-j-1;
const IdxTy yc=y.size()-j-1;
y[yc]=y[yc]+d[i]*x[xpow];
} // j 

} // i 
d=y;
}

void cmd_sum_nest(Cip & cip , LocalVar & lv )
{
const IdxTy n=myatoi(cip.p1);
const IdxTy style=myatoi(cip.p2);
RatCo nest;
nest.push_back(1);
for(IdxTy i=1; i<=n; ++i)
{
nesting(nest);
//std::cout<<ss.str()<<CRLF;
std::cout<<format_poly(i,n,nest,style)<<CRLF;
} //i 
}


void cmd_coef_table(Cip & cip , LocalVar & lv )
{
const IdxTy n=myatoi(cip.p1);
const IdxTy style=myatoi(cip.p2);
for(IdxTy i=1; i<=n; ++i)
{
RatCo x=sum_coefs(i,1);
//std::cout<<ss.str()<<CRLF;
std::cout<<format_poly(i,n,x,style)<<CRLF;
} //i 
}


#else

// TODO needs rational numbers 
#endif // DUPLICATE_RAT_CODE

#endif

std::map<IdxTy, IdxTy> find_multiplicities(const LocTy & c)
{
typedef mjm_ragged_card_itor Mrc;
Mrc mrc(c);
//std::vector<IdxTy> c;
//c.push_back(
std::map<IdxTy, IdxTy> cnts;
for(mrc.reset(); mrc.ok(); ++mrc)
{
//MM_MSG(mti.to_string());
++cnts[mrc.sum()];
}
MM_LOOP(ii,cnts) { MM_MSG(MMPR2((*ii).first,(*ii).second)) }
MM_LOOP(ii,c) { MM_MSG("C "<<MMPR((*ii))) }
return cnts;
}

void cmd_find_mult(Cip & cip , LocalVar & lv )
{
LocTy c;
//for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
//const IdxTy n=myatoi(cip.p1);
IdxTy x=0;
IdxTy i=2;
while ((x=myatoi(cip.wif(i)))!=0) { c.push_back(x); ++i; } 
find_multiplicities(c);
}

#endif




void cmd_write_svg(Cip & cip , LocalVar & lv ) 
{
m_char_mat.cmd_write_svg(cip,lv);
}

void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))

} 

}

void cmd_list(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_dig_map) { MM_MSG("m_dig_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();
/*
ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
ProtraitMap m_protrait_map;
PhenoMap m_pheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;
*/


}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("app-ragged")]=&Myt::cmd_app_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;
m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
m_cmd_map[StrTy("zymo-rags")]=&Myt::cmd_zymo_rags;
m_cmd_map[StrTy("snacks-rags")]=&Myt::cmd_snacks_rags;
m_cmd_map[StrTy("snacks-time")]=&Myt::cmd_snacks_time;
//void cmd_make_latex_table(Cip & cip , LocalVar & lv )
m_cmd_map[StrTy("make-latex-table")]=&Myt::cmd_make_latex_table;
m_cmd_map[StrTy("snacks-txt-svg")]=&Myt::cmd_snacks_txt_svg;
m_cmd_map[StrTy("snacks-txt-svg-i")]=&Myt::cmd_snacks_txt_svg_i;

m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;
#ifdef HAVE_OLD_PROBS
#endif // #ifdef HAVE_OLD_PROBS

#ifdef  DUPLICATE_RAT_CODE
//m_cmd_map[StrTy("string-prob13")]=&Myt::cmd_string_prob13;
//m_cmd_map[StrTy("string-prob14")]=&Myt::cmd_string_prob14;
//void cmd_find_mult(Cip & cip , LocalVar & lv )
m_cmd_map[StrTy("find-mult")]=&Myt::cmd_find_mult;
//void cmd_sum_code(Cip & cip , LocalVar & lv )
//m_cmd_map[StrTy("sum-code")]=&Myt::cmd_sum_code;
#ifdef USE_RATIONALS
m_cmd_map[StrTy("coef-table")]=&Myt::cmd_coef_table;
m_cmd_map[StrTy("sum-nest")]=&Myt::cmd_sum_nest;
#endif

#endif // DUPLICATE_RAT_CODE

m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
TaxTree & tt = m_tax_tree;
StrTy local_label="tat";
//typedef void (Tsrc::* TargCmd)( ListTy & choices,  const char * frag);
//typedef void (Tsrc::* TargParam)( ListTy & choices, const char *  cmd, const char * frag);
m_cli.set_target(*this);
//void set_command_handler(TargCmd * p ) { m_targ_cmd=p; }
//void set_param_handler(TargParam * p ) { m_targ_param=p; }
m_cli.set_command_handler(&Myt::cli_cmd);
m_cli.set_param_handler(&Myt::cli_param);
//std::vector<StrTy> some;
//some.push_back(StrTy("load-tax-nodes"));
//m_cli.set_list(some);
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
os<<ss.str();

}
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



// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code..  
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
bool done() const  { return m_done; } 

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




/////////////////////////////////////////////////////////////////////
private:
void Init()
{
m_flp_def.defaults();
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



bool m_done;

LoRes m_lo;
ParamGlob& m_flp; // m_flp_def;
ParamGlob& m_flp_def;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
FastaMap m_fasta_map;
RaggedMap m_ragged_map;
//CoefficientsMap m_coefs_map;
DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
Units m_units;
//MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;

#ifdef USE_RATIONALS
TableMan m_tables;
#endif


CounterMap m_cm;
CliTy m_cli;

}; //mjm_linc_graph



/////////////////////////////////////////////////////////

#ifdef  TEST_LINC_GRAPH__
int main(int argc,char **args)
{
typedef mjm_linc_graph  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

