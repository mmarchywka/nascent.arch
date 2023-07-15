#ifndef MJM_VIOLIN_H__
#define MJM_VIOLIN_H__
 
#include "mjm_globals.h"


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

//#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
#include "../mjm_fasta_ii.h"

//#include "mjm_collections.h"
#include "mjm_tokenized_collections.h"
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

class violin_box_etc_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
violin_box_etc_params( const StrTy & nm) : Super(nm) {}
violin_box_etc_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
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
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace violin_box_etc_traits
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
typedef  violin_box_etc_traits::Tr  Tr;
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
typedef violin_box_etc_params ParamGlob;

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

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

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

}; // mjm_digraph end

#endif




class mjm_violin_box_etc 
{
typedef violin_box_etc_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_violin_box_etc Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;


typedef mjm_canned_methods Canned;


typedef violin_box_etc_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;



typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

//typedef mjm_clustalw_aln Aln;
//typedef std::map<StrTy, Aln> AlnMap;

typedef mjm_digraph Dig;
typedef std::map<StrTy, Dig> DigMap;


//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

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

int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }

public :
//class mjm_violin_box_etc 
mjm_violin_box_etc():m_dmel(new Dmel()) {Init();}
mjm_violin_box_etc(int argc,char **_args) : m_dmel(new Dmel())
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
if (i==istart) { MM_ERR(" did nothing with "<<args[i]) ++i;  } 

}
}
~mjm_violin_box_etc()
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
if (s=="-quit") { MM_ERR(" FUDD "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
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



void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_stream_edit_fasta( cip ,  lv ); 
}

void cmd_read_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_read_fasta(cip ,  lv,  m_fasta_map  ); }

void cmd_write_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_write_fasta(cip ,  lv,  m_fasta_map  );

}
bool flag_bit(const IdxTy flag, const IdxTy bit, const bool pol=true)
{
const bool x=((flag&(1<<bit))!=0);
return pol?x:!x;
}

typedef mjm_group_sample_map GsMap;
class gsm_loader
{
typedef TokRagged InpTy;
public:
gsm_loader() { Init(); } 

void setup(TokRagged & r, Ragged & pr)
{
ReadWriteMap&  rwm= m_rwm;
pr.to_map(rwm);
rwm.get("group_idx",group_idx);
rwm.get("sample_idx",sample_idx);
rwm.get("angroup_idx",angroup_idx);
rwm.get("otu_idx",otu_idx);
rwm.get("otu_long_idx",otu_long_idx);
rwm.get("val_idx",val_idx);
}

void load(GsMap & gsm , TokRagged & qio)
{
const IdxTy sz=qio.size();
for(IdxTy i=0; i<sz; ++i)
{
const InpTy::Line & l = qio.line(i);
// yes tokenizing a flot is stupid 
//gsm.add(l[group_idx],l[sample_idx],l[otu_idx],atof(qio.st(l[val_idx]).c_str()));
gsm.add(l[group_idx],l[sample_idx],l[angroup_idx],l[otu_idx],atof(qio.st(l[val_idx]).c_str()));
} // i 
}

private:

void Init()
{
group_idx=0;
sample_idx=1;
angroup_idx=2;
otu_long_idx=3;
otu_idx=4;
val_idx=5;
}

ReadWriteMap m_rwm;
 IdxTy group_idx;
 IdxTy sample_idx;
 IdxTy angroup_idx;
 IdxTy otu_long_idx;
 IdxTy otu_idx;
 IdxTy val_idx;

}; // gsm_loader

void cmd_group_stats(Cip & cip , LocalVar & lv )
{

typedef TokRagged InpTy;
typedef Ragged ParTy;
const StrTy ragin=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const StrTy fn=cip.wif(5);
const bool just_dump_rag=((flags&(1<<16))!=0);
const bool just_dump_otu=false;
const bool write_svg_file=true;
InpTy & qio=m_tokragged_map[ragin];
ParTy & params=m_ragged_map[cmd2];
//const IdxTy nval=myatoi(cmd1);
MM_ERR("cmd_group_stats "<<MMPR4(ragin,cmd1,cmd2,flags)<<MMPR(fn))
GsMap gsm(&qio.st());
gsm_loader gl;
gl.setup(qio,params);
gl.load(gsm,qio);
//gsm.calculate();
if (just_dump_otu) { gsm.calculate_otu();}
if (write_svg_file)
{
Ragged rt;
gsm.calculate_otu(rt);
Ragged  rout; // =m_ragged_map[name];
if (just_dump_rag)
{
Ragged& r=  rt; // =m_ragged_map[name];
r.sep(" ");
//if (tabsep) r.sep("\t");
//if (csv) r.sep(",");
//if (ignore_hash) r.ignore_hash(ignore_hash);
//if (load_string) r.load_from_string(fn,debug);
//else 
//if (load_parsed) r.load(fn,debug);
//else r.load_lines(fn,debug);
const IdxTy wflags=0;
r.write_file(fn,wflags);
return;
} // just_dump_rag


rt.transpose(rout);

Ragged & pr=m_ragged_map[cmd1];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
layout_blocks lb;

lb.setup(rout,pr);
lb.find_statistics(rout);
lb.sort(rout);
mjm_svg_writer sw;
std::ofstream  os(fn);
lb.write_svg(os,sw,flags);
} // write_svg_file

} // cmd_group_stats


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
latex_blocks::table_header(os, nval, val_names);
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


latex_blocks::table_end(os);
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

void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
{ Canned::cmd_add_to_fasta(cip ,  lv,  m_fasta_map  ); }

void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
{ Canned::cmd_zymo_merge_fasta(cip ,  lv,  m_fasta_map, m_ragged_map  ); }

void cmd_write_svg_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const StrTy prag=(cip.wif(4));
Ragged & r=m_ragged_map[name];
Ragged & pr=m_ragged_map[prag];
MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
layout_blocks lb;

lb.setup(r,pr);
lb.find_statistics(r);
lb.sort(r);
mjm_svg_writer sw;
std::ofstream  os(fn);
lb.write_svg(os,sw);

}



//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
void cmd_write_svg_ragged_orig(Cip & cip , LocalVar & lv ) 
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
os<<sw.frame_text("#FFFFFFFF",xs,ys);
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
//MM_ERR(" SHOT FURCJK RRRR FUDD "<<MMPR3(i,j,line.size()))
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
//MM_ERR(" SHOT FURCJK RRRR FUDD "<<MMPR3(i,j,line.size()))
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
// this is negative need to make it fudding more negative fudd 
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
// this ASSFUDD does not anchor on fudding right or fudding levt 
os<<sw.gtext_text(text, xleg, y, szleg, colorleg,  StrTy("end") ,0); 
os<<CRLF;

}

for(IdxTy i=0; i<cats; ++i)
{
const IdxTy is=so[i];
const Ragged::Line line=r.line(is);
if (line.size()<=text_idx) continue; // fudd 
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
{ Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ; } // transpose
void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_add_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ; }


void cmd_transpose_tragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_transpose_ragged(cip ,  lv, m_tokragged_map  ) ; } // transpose
void cmd_read_tragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_tokragged_map  ) ; }
void cmd_write_tragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_tokragged_map  ) ; }
void cmd_add_tragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_add_ragged( cip ,  lv, m_tokragged_map  ) ; }




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
{ Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }




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
MM_LOOP(ii,m_tokragged_map) { MM_MSG("m_tokragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
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
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;

m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
m_cmd_map[StrTy("zymo-rags")]=&Myt::cmd_zymo_rags;
m_cmd_map[StrTy("group-stats")]=&Myt::cmd_group_stats;

//m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;

m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
//TaxTree & tt = m_tax_tree;
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

/*

if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }

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
os<<ss;

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

ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
FastaMap m_fasta_map;
RaggedMap m_ragged_map;
TokRaggedMap m_tokragged_map;
DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
//MiiMap m_luts;
//TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;

#ifdef USE_RATIONALS
TableMan m_tables;
#endif


CounterMap m_cm;
CliTy m_cli;

}; //mjm_violin 



/////////////////////////////////////////////////////////

#ifdef  TEST_VIOLIN_BOX__
int main(int argc,char **args)
{
typedef mjm_violin_box_etc Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

