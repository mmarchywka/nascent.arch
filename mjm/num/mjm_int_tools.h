#ifndef MJM_INT_TOOLS_H__
#define MJM_INT_TOOLS_H__

#include "mjm_globals.h"

#include "mjm_globals.h"
#define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
#include "mjm_bigint_mmpr.h"
//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
#include "mjm_integrals.h"
#include "../copied_code/gmp_factorize.h"

#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"

// at least for algorithm testing 
#include "mjm_data_model_error_log.h"
// some of these pickup the math libs... 
// fraction input support 
// for now just code here, not too hard 
// need this for Alfredo lol
//#include "mjm_rational.h"
 // #include "mjm_generic_iterators.h"
// not really used but tyedefed
//#include "mjm_block_matrix.h"
//#include "mjm_interpolation.h"
//#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
// try to use the home made version of erf(x1)-erf(x2)
//#include "mjm_integrals.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>

/*
2017-10-27
copied from mjm_snacks.h
copied pollard rho from demo progam 
https://gmplib.org/manual/Demonstration-Programs.html


*/

/*
g++ -DTEST_INT_TOOLS__ -std=gnu++11 -I gmp/gmp-6.1.2 -Lgmp/gmp-6.1.2/.libs -lgmp -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -gdwarf-3 -O0 -I.. -Wall -Wno-unused-variable -Wno-unused-function -Wno-non-template-friend -x c++ mjm_int_tools.h -o mjm_int_tools.out

*/

////////////////////////////////////////////////////////////////

class int_tool_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
int_tool_params( const StrTy & nm) : Super(nm) {}
int_tool_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  

//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
//bool always_dump_to_file() const { return m_map.get_bool("always_dump_to_file",!true); }
/*
bool skip_old() const { return m_map.get_bool("skip_old",true); }
bool print_dog_days() const { return m_map.get_bool("print_dog_days",!true); }
bool accumulate_dog_days() const { return m_map.get_bool("accumulate_dog_days",true); }
*/
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
/*
StrTy snack_log() const { return m_map.get_string("snack_log","snacks_log"); }
		//if (skip_old) if (date==StrTy("2017-04-22")) skipping=false;
StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
*/
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

//ss<<"time_step="<<time_step()<<sep;
/*
ss<<"skip_old="<<skip_old()<<sep;
ss<<"print_dog_days="<<print_dog_days()<<sep;
ss<<"accumulate_dog_days="<<accumulate_dog_days()<<sep;
*/
ss<<"log_commands"<<log_commands()<<sep;
/*
ss<<"snack_log="<<snack_log()<<sep;
ss<<"start_date="<<start_date()<<sep;
ss<<"end_date="<<end_date()<<sep;
*/
return ss.str();
}


}; // snack_params





class mjm_int_tool
{
/*

*/

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef mjm_block_matrix<D> MyBlock;
//typedef mjm_sparse_matrix<D> MySparse;
}; // 



typedef mjm_int_tool Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

//typedef mjm_logic_base Logic;
typedef int_tool_params Logic;
typedef mjm_logic_base VariableStore;
//typedef std::map<StrTy,StrTy > DateInfo;
//typedef std::map<StrTy,StrTy > Canonical;
//typedef std::vector<StrTy> CanonFails;
//typedef std::vector<StrTy> DumpOrder;
//typedef std::map<StrTy,IdxTy > FreqInfo;
//typedef std::map<StrTy,IdxTy > Literals;

//typedef literal_map LiteralClass;
//typedef LiteralClass::LITERAL_TYPE LITERAL_TYPE;
//typedef literal_map LiteralClass;

// collection of things for each date 
//typedef std::map<StrTy, DateInfo > MappedInfo;
//typedef std::map<StrTy, D > ShareFactor;
//typedef std::map<StrTy, ShareFactor > ShareFactorMap;

typedef data_model_error_log Dmel;

/// need to have dog_day def available. 
// MACRO is indexed by name and then date 

public :
mjm_int_tool():m_dmel(new Dmel()) {Init();}
mjm_int_tool(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_int_tool()
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
dest=::atoi(args[i]);
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
//if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
//if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { ++i; clean_up(); }
} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 

//void start_date(const StrTy & d) { m_flp.set("start_date",d); }
//void end_date(const StrTy & d) { m_flp.set("end_date",d); }

void command_modef(const char * fn)
{
std::ifstream fin(fn);
CommandInterpretter li(&fin);
command_mode(li);


}
void command_mode()
{
//LineIterator 
CommandInterpretter li(&std::cin);
command_mode(li);

}
void command_mode(const StrTy & cmd)
{

CommandInterpretter li;
li.set(cmd,1);
//li.set(cmd,0);
command_mode(li);
}


void command_mode(CommandInterpretter & li)
{
StrTy local_label="fick";
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
if (m_flp.log_commands()) 
		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz<1) continue;
const StrTy cmd=li.word(0);
const StrTy p1=(sz>1)?li.word(1):StrTy("");
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
//if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
//if (cmd=="parse") { parse();  continue; } 
//if (cmd=="pair") { pair_stats(li.word(1),li.word(2),li.word(3));  continue; } 
//if (cmd=="pair-all") { pair_stats(li.word(1));  continue; } 
//if (cmd=="dose-vectors") { if (sz>1) { dose_vectors(li.word(1));} else { MM_ERR( "need dogname ") }  continue; } 
//void pair_stats(const StrTy & x, const StrTy & y)
/*
if (cmd=="dump-drugs") { collect_all_drugs(); continue; } 
if (cmd=="dump-drugs-cerr") { collect_all_drugs(2); continue; } 
if (cmd=="list-drugs-cerr") { collect_all_drugs(3); continue; } 
if (cmd=="list-dogs-cerr") { collect_all_drugs(4); continue; } 
if (cmd=="dump-counts") { stream_counts(std::cout); continue; } 
if (cmd=="dump-canon-fails") { dump_canon_fails(std::cout); continue; } 
if (cmd=="clear-canon-fails") { clear_canon_fails(); continue; } 
if (cmd=="dump-dmel") { dump_dmel(std::cout); continue; } 
if (cmd=="dump-dmel-cerr") { dump_dmel(std::cerr); continue; } 
if (cmd=="dump-dog-days") { dump_dog_days(std::cout,p1); continue; } 
if (cmd=="n-canon-fails") { MM_MSG(MMPR(m_canon_fails.size())) continue; } 

if (cmd=="clear-order") { clear_order(); continue; } 
if (cmd=="dump-order") { dump_order(std::cout); continue; } 
if (cmd=="push-order") 
		//{ int i=1; while ( li.cmd_ok(i+1)){ push_order(li.word(i)); ++i; } continue; } 
		{ int i=1; while ( i<sz){ push_order(li.word(i)); ++i; } continue; } 
if (cmd=="remove-order") 
		{ int i=1; while ( i<sz){ remove_order(li.word(i)); ++i; } continue; } 

if (cmd=="load-order") { if( li.cmd_ok(2)) load_order(li.word(1)); continue; } 
if (cmd=="order-all") { collect_all_drugs(1); continue; } 
*/
//if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
//if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
//if (cmd=="process") { if (li.cmd_ok(2)) process_file(li.word(1));  continue; } 
//if (cmd=="add") { if (li.cmd_ok(2)) add_timeline(li.word(1));  continue; } 
if (cmd=="init") { Init();  continue; } 
//if (cmd=="hlatex") { std::cout<<m_chart.to_latex_h();  continue; } 
//if (cmd=="hssv") { std::cout<<m_chart.to_ssv(true);  continue; } 
//if (cmd=="vssv") { std::cout<<m_chart.to_ssv(!true);  continue; } 
//if (cmd=="add-event") { if (li.cmd_ok(3)) add_event(li.word(1),li.word(2));  continue; } 
//if (cmd=="add-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
//if (cmd=="set-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
//if (cmd=="reset-sum") { if (li.cmd_ok(2)) reset_sum(li.word(1));  continue; } 
//if (cmd=="make") { make_chart();  continue; } 
//if (cmd=="unused") { dump_unused();  continue; } 
//if (cmd=="clear-unused") { m_unused.clear();  continue; } 
//if (cmd=="legend") { dump_legend(std::cout,local_label,1);  continue; } 
//if (cmd=="dump-human") { dump_state_etc(std::cout,local_label,0);  continue; } 
//if (cmd=="init") { Init();   continue; } 
//if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
/*
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
*/
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");

}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 


// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
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
//MM_MSG(" configuration "<<m_flp.to_string())
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}
bool done() const  { return m_done; } 
////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
private:
void Init()
{
//delete m_dmel;
//m_dmel= new Dmel();
m_done=false;
}
bool m_done;
Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;
CounterMap m_cm;

}; //mjm_int_tools



/////////////////////////////////////////////////////////

#ifdef  TEST_INT_TOOLS__
int main(int argc,char **args)
{
typedef mjm_int_tool  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

