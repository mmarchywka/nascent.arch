
#ifndef MJM_muqed_util_H__
#define MJM_muqed_util_H__


#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_canned_methods.h"
#include "mjm_cli_ui.h"
#include "mjm_tokenized_collections.h"

#include "mjm_globals.h"
#include "mjm_generic_recipe.h"
#include "mjm_compiled_noun.h"
#include "mjm_thread_util.h"
#include "mjm_ragged_forms.h"
#include "mjm_snacks_accumulator.h"
#include "mjm_constrained_file.h"
#include "mjm_unit_crap.h"
// faster than calling date in bash doh 
#include "mjm_calendar.h"
#include "mjm_diary_parse_state.h"

#include "mjm_diet_diary_form.h"



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
made from wizard script on Sun Oct 25 08:33:20 EDT 2020
muqed_util
*/


////////////////////////////////////////////////////////////////

class muqed_util_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
muqed_util_params( const StrTy & nm) : Super(nm) {}
muqed_util_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
//StrTy daily= m_flp.dog_daily(); "dog_daily.txt";
StrTy dog_daily() const { return m_map.get_string("dog_daily","dog_daily.txt"); }
//StrTy glob=m_flp.glob(); // "dog_glob.txt";
StrTy glob() const { return m_map.get_string("glob","dog_glob.txt"); }
//StrTy usedd=m_flp.dog_used(); // "dog_used.txt";
StrTy dog_used() const { return m_map.get_string("dog_used","dog_used.txt"); }
// samples, 
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
StrTy start_date() const { return m_map.get_string("start_date",""); }
StrTy end_date() const { return m_map.get_string("end_date",""); }
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
ss<<"dog_daily"<<dog_daily()<<sep;
ss<<"glob"<<glob()<<sep;
ss<<"dog_used"<<dog_used()<<sep;
ss<<"start_date"<<start_date()<<sep;
ss<<"end_date"<<end_date()<<sep;
return ss.str();
}


}; // muqed_util_params


namespace muqed_util_traits
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


}; // muqed_util_traits
///////////////////////////////////////////////////////////////




class mjm_muqed_util 
{
typedef muqed_util_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_muqed_util Myt;
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
typedef muqed_util_params ParamGlob;
typedef mjm_logic_base VariableStore;
typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
//typedef mjm_tax_tree TaxTree;
//typedef std::map<StrTy, TaxTree> TaxTrees;
//typedef std::map<IdxTy,IdxTy> Mii;
//typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;

typedef mjm_diet_diary_form<Tr>  DiaryForm;
typedef DiaryForm::parse_settings PO;
typedef DiaryForm::parse_state PS;

typedef mjm_loo_parsing<Tr> Loo;

public:
// FIXME doh put this somwhere lol 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }

public :
mjm_muqed_util():m_dmel(new Dmel()) {Init();}
mjm_muqed_util(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_muqed_util()
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
if (s=="-quit") { MM_ERR(" quit  "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
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



//void cmd_tt(Cip & cip , LocalVar & lv ) { Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }


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
// use the loo 
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tokragged_map) { MM_MSG("m_tokragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();


}
/////////////////////////////////////////////////////////


void cmd_about(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_about "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("about" ) } 

about();
} // cmd_about

void cmd_dates(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_dates "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("dates" ) } 
//m_startdate=cip.p1; m_enddate=cip.p2; 
start_date(cip.p1); end_date(cip.p2); 


} // cmd_dates
#if 0 
void cmd_quit(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_quit "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("quit" ) } 


} // cmd_quit
#endif

void cmd_dump(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_dump "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("dump" ) } 
MM_ERR(m_diary.dump())

} // cmd_dump

void cmd_parse(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_parse "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("parse" ) } 


} // cmd_parse

void cmd_parse_week(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_parse_week "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("parse_week" ) } 


} // cmd_parse_week

void cmd_blank(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_blank "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("blank" ) } 


} // cmd_blank

void cmd_markup(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_markup "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("markup" ) } 
//if (cmd=="markup")
{PO po(1<<PO::MARKUP); Myt::PS ps; ps.debug_os=&std::cerr;
const StrTy startdate=m_flp.start_date(); // m_startdate; // cip.wif(3);
const StrTy enddate=m_flp.end_date(); // m_enddate; // cip.wif(4);
po.date_range(startdate,enddate);
m_diary.parse(ps,po);
//StrTy fn=cmdp1;
//typedef std::stringstream Ss;
MM_ERR(" dump of markuyp ")
//Ss ss; x.save(ss,*(ps.markup)); MM_ERR(ss.str()) 
if ((fn=="-")||(fn=="")) {Ss ss; m_diary.save(ss,*(ps.markup)); MM_ERR(ss.str()) }
else
{ std::ofstream ofs(fn); m_diary.save(ofs,*(ps.markup)); MM_ERR(" saved to "<<MMPR(fn))}

} // markup

if (debugg()) {  MM_ERR("end markup" ) } 

} // cmd_markup

void cmd_markup_week(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_markup_week "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("markup_week" ) } 

//if (cmd=="markup")
PO po(1<<PO::MARKUP); Myt::PS ps; ps.debug_os=&std::cerr;
m_diary.parse_week(ps,po);
//StrTy fn=cmdp1;
//typedef std::stringstream Ss;
MM_ERR(" dump of markuyp ")
//Ss ss; x.save(ss,*(ps.markup)); MM_ERR(ss.str()) 
if ((fn=="-")||(fn=="")) {Ss ss; m_diary.save(ss,*(ps.markup)); MM_ERR(ss.str()) }
else
{ std::ofstream ofs(fn); m_diary.save(ofs,*(ps.markup)); MM_ERR(" saved to "<<MMPR(fn))}


} // cmd_markup_week
#if 0

if (cmd=="commit")
{
x.clear_form();
x.load_form();
StrTy fn=cmdp1;
std::ifstream ifs(fn); x.load_form(ifs);
x.save_form();
MM_ERR(" form updated")

} // commit 
if (cmd=="template")
{
StrTy fn=cmdp1;
PO po(1<<PO::BLANK); Myt::PS ps; ps.debug_os=&std::cerr;
//x.parse(ps,po);
x.parse_week(ps,po); // may not do whole anyway lol 
typedef std::stringstream Ss;
if ((fn=="-")||(fn=="")) {Ss ss; x.save(ss,*(ps.blank)); MM_ERR(ss.str()) }
else
{ std::ofstream ofs(fn); x.save(ofs,*(ps.blank)); MM_ERR(" saved to "<<MMPR(fn))}
} // template

#endif

void cmd_commit(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_commit "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("commit" ) } 
m_diary.clear_form();
m_diary.load_form();
std::ifstream ifs(fn); m_diary.load_form(ifs);
m_diary.save_form();
MM_ERR(" form updated")



} // cmd_commit

void cmd_template(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_template "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("template" ) } 
PO po(1<<PO::BLANK); Myt::PS ps; ps.debug_os=&std::cerr;
//x.parse(ps,po);
m_diary.parse_week(ps,po); // may not do whole anyway lol 
typedef std::stringstream Ss;
if ((fn=="-")||(fn=="")) {Ss ss; m_diary.save(ss,*(ps.blank)); MM_ERR(ss.str()) }
else
{ std::ofstream ofs(fn); m_diary.save(ofs,*(ps.blank)); MM_ERR(" saved to "<<MMPR(fn))}


} // cmd_template

void cmd_eval(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_eval "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("eval" ) } 


} // cmd_eval

void cmd_eval_week(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.cmd();
IdxTy flags=myatoi(cip.p1.c_str());
IdxTy flagout=myatoi(cip.p2.c_str());
const StrTy startdate=m_flp.start_date(); // m_startdate; // cip.wif(3);
const StrTy enddate=m_flp.end_date(); // m_enddate; // cip.wif(4);
StrTy daily= m_flp.dog_daily(); //  "dog_daily.txt";
StrTy glob=m_flp.glob(); // "dog_glob.txt";
//StrTy usedd="dog_used.txt";
StrTy usedd=m_flp.dog_used(); // "dog_used.txt";
const bool dump_dmglob=Bit(flagout,0);
const bool dump_rest=Bit(flagout,1);
const bool save_dmglob=Bit(flagout,2);
const bool save_daily=Bit(flagout,3);
const bool save_used=Bit(flagout,4);

if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_eval_week "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("eval_week"<<MMPR4(cmd,flags,flagout,startdate )<<MMPR4(enddate,dump_dmglob,dump_rest,save_dmglob)<<MMPR2(save_daily,save_used)) } 
IdxTy idx=3;
if (save_dmglob)
{ if (cip.wif(idx)!="") if (cip.wif(idx)!="-") glob=cip.wif(idx); ++idx;  }
if (save_daily)
{ if (cip.wif(idx)!="") if (cip.wif(idx)!="-") daily=cip.wif(idx); ++idx;  }
if (save_used)
{ if (cip.wif(idx)!="") if (cip.wif(idx)!="-") usedd=cip.wif(idx); ++idx;  }

const bool week=(cmd=="eval-week");

MM_ERR(
MMPR4(flags,flagout,daily,glob)<<MMPR4(usedd,dump_dmglob,dump_rest, save_dmglob)<<MMPR2(save_daily,save_used))
//MM_FAULT
PO po(1<<PO::EVAL);
PS ps; ps.debug_os=&std::cerr;
if (!week)
{
po.date_range(startdate,enddate);
MM_ERR(MMPR2(startdate,enddate))
m_diary.parse(ps,po);
}else m_diary.parse_week(ps,po);
typedef std::stringstream Ss;
if (dump_dmglob) {
Ss ss; ss<<ps.dmglob->dump();
//x.save(ss,*(ps.blank));
MM_ERR("%%%%% dump of dmglob ") MM_ERR(ss.str())
} // dump_dmglob
std::vector<StrTy> no;

typedef mjm_ragged_table Ragged;
Ragged d,used;
// void dog_daily(Ragged & d, Ragged * used,noun_order &no, const IdxTy flags)
ps.dmglob->dog_daily(d, & used, no,flags);
if (dump_rest) {
Ss rr;
rr<<"%%%%% dump of d or daily dog report "<<CRLF ;
d.dump_os(rr,3);
rr<<"%%%%% dump of used or inventory report  "<<CRLF ;
used.dump_os(rr,3);
MM_ERR(rr.str())
 } // dump_rest

if (save_dmglob) {std::ofstream osf(glob.c_str());
osf<<(*(ps.dmglob)).dump(); }
if (save_daily) {std::ofstream osf(daily.c_str()); d.dump_os(osf,3); }
if (save_used) {std::ofstream osf(usedd.c_str()); used.dump_os(osf,3); }



} // cmd_eval_week
#if 0 
void cmd_load(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load" ) } 


} // cmd_load

if (cmd=="load-form") { x.set_file(cip.p1); x.load_form();
// MM_ERR("loading "<<MMPR(x.dump()))  
}
if (cmd=="load-units") { x.load_units(cip.p1);
// MM_ERR("loading "<<MMPR(x.dump()))  
}
if (cmd=="load-nouns") { x.load_nouns(cip.p1);
//MM_ERR("loading "<<MMPR(x.dump()))  
}
if (cmd=="load-ignores") { x.load_ignores(cip.p1);
//MM_ERR("loading "<<MMPR(x.dump()))  
}
if (cmd=="load-reserveds") { x.load_reserveds(cip.p1);
// MM_ERR("loading "<<MMPR(x.dump())) 
 }
if (cmd=="load-recipes") { x.load_recipes(cip.p1);  MM_ERR("loading "<<MMPR(x.dump()))  }
if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="nouns") {x.index();  MM_ERR(x.dump_nouns()) }





#endif






void cmd_load_form(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_form "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_form" ) } 

m_diary.set_file(cip.p1); 
m_diary.load_form();

} // cmd_load_form

void cmd_load_units(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_units "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_units" ) } 
m_diary.load_units(cip.p1);


} // cmd_load_units

void cmd_load_finder(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_finder "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_finder" ) } 

m_diary.load_finder(cip.p1);

} // cmd_load_finder


void cmd_load_nouns(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_nouns "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_nouns" ) } 

m_diary.load_nouns(cip.p1);

} // cmd_load_nouns

void cmd_load_ignores(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_ignores "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_ignores" ) } 

m_diary.load_ignores(cip.p1);

} // cmd_load_ignores

void cmd_load_reserveds(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_reserveds "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_reserveds" ) } 

m_diary.load_reserveds(cip.p1);

} // cmd_load_reserveds

void cmd_load_recipes(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_recipes "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load_recipes" ) } 

m_diary.load_recipes(cip.p1);  //MM_ERR("loading "<<MMPR(x.dump()))  }

} // cmd_load_recipes
#if 0 
void cmd_dump(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_dump "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("dump" ) } 


} // cmd_dump
#endif

void miss_template(Ss & ss, const StrTy &no, const IdxTy t ) const
{
switch (t)
{
// "Vitamin D" canon D-3 expand 1 rpart ",D-3,1"
case 1 :{
ss<<"\""<<no<<"\" canon \""<<no<<"\" flags syn runits mg ";
break; }  // case

default: // infinite loop ... 
ss<<"\""<<no<<"\" canon \""<<no<<"\" expand 1 flags syn part \","<<no<<",1\"";


} // switch 
} // miss_template

void cmd_missing(Cip & cip , LocalVar & lv )
{
const IdxTy t=myatoi(cip.p1);
//const IdxTy flags=myatoi(cip.p1);
//const bool dump=!Bit(flags,0);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_missing "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("missing "<<MMPR(t) ) } 
m_diary.index();  
const auto  m=m_diary.missing();
// "Vitamin D" canon D-3 expand 1 rpart ",D-3,1"
MM_LOOP(ii,m) { MM_ERR(MMPR2((*ii).first,(*ii).second)) } // m 
std::ostream &os=std::cerr;
MM_LOOP(ii,m) {Ss ss; miss_template(ss,(*ii).first,t); os<<ss.str()<<CRLF; } // m 


//Loo::Dump(std::cerr,"missing",m);
//if (dump) { MM_ERR(m_diary.dump_nouns())  } 

} // cmd_missing


void cmd_nouns(Cip & cip , LocalVar & lv )
{
//const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p1);
const bool dump=!Bit(flags,0);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_nouns "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("nouns" ) } 
m_diary.index();  
if (dump) { MM_ERR(m_diary.dump_nouns())  } 

} // cmd_nouns

void cmd_expand(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_expand "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("expand" ) } 
MM_ERR(" expand "<<MMPR3(cip.line(),cip.p1, cip.p2));
m_diary.set_expansion(cip.p1,myatoi(cip.p2.c_str()));
} // cmd_expand

void cmd_xxx(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_xxx "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("xxx" ) } 


} // cmd_xxx

void cmd_load(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_load "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("load" ) } 


} // cmd_load

void cmd_clear(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_clear "<<CRLF;

cip.help(ss.str());
return;
}
if (debugg()) {  MM_ERR("clear" ) } 


} // cmd_clear

///////////////////////////////////////////////////////

static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 
bool debugg() const { return true; } 
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
//m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;
/////////////////////////////////////////////////////////////////

m_cmd_map[StrTy("about")]=&Myt::cmd_about;
m_cmd_map[StrTy("dates")]=&Myt::cmd_dates;
m_cmd_map[StrTy("quit")]=&Myt::cmd_quit;
m_cmd_map[StrTy("dump")]=&Myt::cmd_dump;
m_cmd_map[StrTy("parse")]=&Myt::cmd_parse;
m_cmd_map[StrTy("parse-week")]=&Myt::cmd_parse_week;
m_cmd_map[StrTy("blank")]=&Myt::cmd_blank;
m_cmd_map[StrTy("markup")]=&Myt::cmd_markup;
m_cmd_map[StrTy("markup-week")]=&Myt::cmd_markup_week;
m_cmd_map[StrTy("commit")]=&Myt::cmd_commit;
m_cmd_map[StrTy("template")]=&Myt::cmd_template;
m_cmd_map[StrTy("eval")]=&Myt::cmd_eval_week;
m_cmd_map[StrTy("eval-week")]=&Myt::cmd_eval_week;
m_cmd_map[StrTy("load")]=&Myt::cmd_load;
m_cmd_map[StrTy("load-form")]=&Myt::cmd_load_form;
m_cmd_map[StrTy("load-units")]=&Myt::cmd_load_units;
m_cmd_map[StrTy("load-nouns")]=&Myt::cmd_load_nouns;
m_cmd_map[StrTy("load-finder")]=&Myt::cmd_load_finder;
m_cmd_map[StrTy("load-ignores")]=&Myt::cmd_load_ignores;
m_cmd_map[StrTy("load-reserveds")]=&Myt::cmd_load_reserveds;
m_cmd_map[StrTy("load-recipes")]=&Myt::cmd_load_recipes;
m_cmd_map[StrTy("dump")]=&Myt::cmd_dump;
m_cmd_map[StrTy("missing")]=&Myt::cmd_missing;
m_cmd_map[StrTy("nouns")]=&Myt::cmd_nouns;
m_cmd_map[StrTy("expand")]=&Myt::cmd_expand;
m_cmd_map[StrTy("xxx")]=&Myt::cmd_xxx;
m_cmd_map[StrTy("load")]=&Myt::cmd_load;
m_cmd_map[StrTy("clear")]=&Myt::cmd_clear;



//////////////////////////////////////////////////////////////
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
ss<<" mjm_muqed_util built  "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com started  Sun Oct 25 08:33:20 EDT 2020 "<<CRLF;
ss<<" Empty skeleton code add credits like this,   "<<CRLF;
//ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;

std::ostream & os=std::cout;
os<<ss.str();

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

// members MEMBERS

bool m_done;

ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
TokRaggedMap m_tokragged_map;

DiaryForm m_diary;
//StrTy m_startdate,m_enddate;
//DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
//CharMat m_char_mat;
//MiiMap m_luts;
//TaxTree m_tax_tree; // now have sequences ID's to taxon 
//TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;

}; //mjm_muqed_util 

/////////////////////////////////////////////////////////

#ifdef  TEST_muqed_util__
int main(int argc,char **args)
{
typedef mjm_muqed_util Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif
