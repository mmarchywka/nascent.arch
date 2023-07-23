
#ifndef MJM_string_tax_H__
#define MJM_string_tax_H__
 
#include "mjm_globals.h"

#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"

#include "mjm_canned_methods.h"


#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"


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
made from wizard script on Wed Mar  2 07:11:16 EST 2022
string_tax
*/


////////////////////////////////////////////////////////////////

class string_tax_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
string_tax_params( const StrTy & nm) : Super(nm) {}
string_tax_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
// samples, 
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }

// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
//ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;
return ss.str();
}


}; // string_tax_params


namespace string_tax_traits
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
typedef IdxTy FlagTy;
typedef  data_model_error_log Dmel; 
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::FlagTy FlagTy;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // string_tax_traits
///////////////////////////////////////////////////////////////




class mjm_string_tax 
{
typedef string_tax_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_string_tax Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::FlagTy FlagTy;
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
typedef string_tax_params ParamGlob;
typedef mjm_logic_base VariableStore;
typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;

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


public:
// FIXME doh put this somwhere lol 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }

public :
mjm_string_tax():m_dmel(new Dmel()) {Init();}
mjm_string_tax(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_string_tax()
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
if (s=="-quit") { MM_ERR(" wtf "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
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
void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_stream_edit_fasta( cip ,  lv ); 
}

void cmd_read_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_read_fasta(cip ,  lv,  m_fasta_map  ); }

void cmd_write_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_write_fasta(cip ,  lv,  m_fasta_map  );

}
*/

bool flag_bit(const IdxTy flag, const IdxTy bit, const bool pol=true)
{
const bool x=((flag&(1<<bit))!=0);
return pol?x:!x;
}


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

} // zymo_rags

//void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_add_to_fasta(cip ,  lv,  m_fasta_map  ); }

//void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_zymo_merge_fasta(cip ,  lv,  m_fasta_map, m_ragged_map  ); }

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

}


void cmd_transpose_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ; } // transpose
void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_dump_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_dump_ragged( std::cout, cip ,  lv, m_ragged_map  ) ; }



void cmd_add_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ; }
//void cmd_transpose_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_transpose_ragged(cip ,  lv, m_tokragged_map  ) ; } // transpose
//void cmd_read_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_read_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_write_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_write_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_add_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_add_ragged( cip ,  lv, m_tokragged_map  ) ; }



void cmd_tt(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }


void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_quit(Cip & cip , LocalVar & lv )
{ clean_up(); return;  }
void cmd_cm(Cip & cip , LocalVar & lv )
{ dump_cm(); return;  }
void cmd_banner(Cip & cip , LocalVar & lv )
{ config_banner(); return;  }

void cmd_set_param(Cip & cip , LocalVar & lv )
{  //kkkkkkkkkkk 
//const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
 if (cip.li().cmd_ok(3))  m_flp.set(cip.li().cmd_set());
        return;  }
void cmd_get_param(Cip & cip , LocalVar & lv )
{
 if (cip.li().cmd_ok(2))
        std::cout<<lv["local_label"]<<" "<<cip.li().word(1)<<" "
                        <<m_flp.get(cip.li().word(1))<<CRLF;
         }





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
//MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_dig_map) { MM_MSG("m_dig_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();


}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;

m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("quit")]=&Myt::cmd_quit;
m_cmd_map[StrTy("cm")]=&Myt::cmd_cm;
m_cmd_map[StrTy("banner")]=&Myt::cmd_banner;
m_cmd_map[StrTy("set-param")]=&Myt::cmd_set_param;
m_cmd_map[StrTy("get-param")]=&Myt::cmd_get_param;




m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("dump-ragged")]=&Myt::cmd_dump_ragged;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;

//m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
//m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
//m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
//m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
//m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
//m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
//m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
//m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
//m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
//m_cmd_map[StrTy("zymo-rags")]=&Myt::cmd_zymo_rags;
//m_cmd_map[StrTy("group-stats")]=&Myt::cmd_group_stats;

//m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;

//m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
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
mloc["local_label"]=local_label;
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
if ( m_done) return;
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
ss<<" mjm_string_tax built  "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com started  Wed Mar  2 07:11:16 EST 2022 "<<CRLF;
ss<<" Empty skeleton code add credits like this,   "<<CRLF;
ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;

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
m_done=false;
}
bool Bit(const FlagTy flags, const IdxTy b) { return (flags&(1<<b))!=0; }

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

// MEMBERS members

bool m_done;

ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
TokRaggedMap m_tokragged_map;
//DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
//CharMat m_char_mat;
//MiiMap m_luts;
//TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;

}; //mjm_string_tax 

/////////////////////////////////////////////////////////

#ifdef  TEST_string_tax__
int main(int argc,char **args)
{
typedef mjm_string_tax Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif
