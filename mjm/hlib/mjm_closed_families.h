#ifndef MJM_CLOSED_FAMILIES_H__
#define MJM_CLOSED_FAMILIES_H__

#include "mjm_globals.h"

#include "mjm_globals.h"
//#define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
//#include "mjm_bigint_mmpr.h"


//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
#include "mjm_integrals.h"
// this will take forever to fix... 
//#include "../copied_code/gmp_factorize.h"

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

2023-03-12 I think these are Hermites? lol 

*/

/*
g++ -DTEST_CLOSED_FAMILIES__ -std=gnu++11 -I gmp/gmp-6.1.2 -Lgmp/gmp-6.1.2/.libs -lgmp -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -gdwarf-3 -O0 -I.. -Wall -Wno-unused-variable -Wno-unused-function -Wno-non-template-friend -x c++ mjm_closed_families.h -o mjm_closed_families.out

*/

////////////////////////////////////////////////////////////////

class closed_families_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
closed_families_params( const StrTy & nm) : Super(nm) {}
closed_families_params() : Super() {}
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

class mjm_closed_families
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



typedef mjm_closed_families Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

//typedef mjm_logic_base Logic;
typedef closed_families_params Logic;
typedef mjm_logic_base VariableStore;

typedef data_model_error_log Dmel;

public :
mjm_closed_families():m_dmel(new Dmel()) {Init();}
mjm_closed_families(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_closed_families()
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
StrTy local_label="closed_families";
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
if (cmd=="parse") { parse_mode(li);  continue; } 
//void add_single_term(const IdxTy  expn, const IdxTy p, const IdxTy flag, const value_type & c)
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



/////////////////////////////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////////////

/*
 A sum over products of a polynomial \sum_n  P(x)exp(-bnx*x) for later use in rational 
quotients that are then closed under division and multiplication and
allowing non-integer powers to be differnetiated 
This is DENSE in both n and the polynomial coefficients although latter should defer
to typedef of poly 
This is closed under multiplication but not division, need the exp_poly_ratio for that 

*/
template <class Tv=double> class poly_norm
{

class Tr{
public:
//typedef double D;
typedef Tv  D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::D D;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::StrTy StrTy;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<IdxTy> MyBlockInt;
typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
typedef mjm_generic_itor<2> LoopItor;

typedef poly_norm<Tv>  Myt;
typedef Tv value_type;
typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;

typedef SFi::simple_polynomial<Tv> Tpoly;
typedef std::vector<Tpoly> PolyVec;
public:
poly_norm():m_beta(1) { Init(); Default(); }
//poly_norm(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S) 
//:m_P(P),m_Q(Q),m_R(R),m_S(S)
//{ Init();  }
poly_norm(const value_type & x): m_beta(1)  {  from_const(x); } 

void from_const(const value_type & x)
{
//Mono(m_P,0,0); Mono(m_Q,0,x); Mono(m_R,0,0); Mono(m_S,0,1);
}
void normal()
{
//Mono(m_P,0,1); Mono(m_Q,0,0); Mono(m_R,0,0); Mono(m_S,0,1);
}
const IdxTy size() const { return m_vec.size(); } 
const Tpoly& operator[](const IdxTy i) const { return m_vec[i]; } 
void push_back(const Tpoly & x) 
{
m_vec.push_back(x); 
}

bool is_zero(const Tpoly & x) const 
{
const IdxTy sz=x.size();
for (IdxTy i=0; i<sz; ++i) if (x[i]!=0) return false;
return true;

}



bool operator==(const Myt & that) const
{
const IdxTy sz=m_vec.size();
const IdxTy szt=that.m_vec.size();
for(IdxTy i=sz; i<szt; ++i) { if (!is_zero(that.m_vec[i])) return false; }
for(IdxTy i=szt; i<sz; ++i) { if (!is_zero(m_vec[i])) return false; }
const IdxTy szz=(sz<szt)?sz:szt;
for(IdxTy i=0; i<szz; ++i)
{
if (m_vec[i]==that.m_vec[i]) continue;
if (!is_zero(m_vec[i])) return false;
if (!is_zero(that.m_vec[i])) return false;

}
return true; 
}

Myt operator*( const value_type & that) const 
{ Tpoly p; p.push_back(that); return (*this)*p; }

Myt operator*( const Tpoly & that) const 
{
Myt y;
const IdxTy sz=m_vec.size();
for (IdxTy i=0; i<sz; ++i)
{
Tpoly p;
mi.multiply_polynomials(p,m_vec[i],that);
y.m_vec.push_back(p);

}
y.trim();
return y;
}

Myt operator*( const Myt & that) const 
{
Myt y;
const IdxTy sz=m_vec.size();
const IdxTy szt=that.m_vec.size();
const int sznet=sz+szt-1;
if (sznet<1) return y;
y.m_vec= PolyVec(sznet);
for (IdxTy i=0; i<sz; ++i)
for (IdxTy j=0; j<szt; ++j)
{
Tpoly p;
// this stpid mi thing needs to be fixed, it could at least be mutable for const crap lol TODO 
y.mi.multiply_polynomials(p,m_vec[i],that.m_vec[j]);
y.mi.accumulate_polynomial(y.m_vec[i+j],p);
//y.m_vec[i+j]+=p;

}
y.trim();
return y;
}

Myt operator+( const Myt & that) const 
{
Myt y;
Tpoly pzero;
const IdxTy sz=m_vec.size();
const IdxTy szt=that.m_vec.size();
const IdxTy szmax=(sz>szt)?sz:szt;
if (szmax<1) return y;
y.m_vec= PolyVec(szmax);
for (IdxTy i=0; i<szmax; ++i)
{
Tpoly p;
const Tpoly & px=(i<sz)?m_vec[i]:pzero;
const Tpoly & py=(i<szt)?that.m_vec[i]:pzero;
mi.add_polynomials(p,px,py);
y.m_vec[i]=p;

}
y.trim();
return y;
}

Myt operator-( const Myt & that) const 
{
Myt y;
Tpoly pzero;
const IdxTy sz=m_vec.size();
const IdxTy szt=that.m_vec.size();
const IdxTy szmax=(sz>szt)?sz:szt;
if (szmax<1) return y;
y.m_vec= PolyVec(szmax);
for (IdxTy i=0; i<szmax; ++i)
{
Tpoly p;
const Tpoly & px=(i<sz)?m_vec[i]:pzero;
const Tpoly & py=(i<szt)?that.m_vec[i]:pzero;
mi.subtract_polynomials(p,px,py);
y.m_vec[i]=p;

}
y.trim();
return y;
}






void trim()
{
const IdxTy sz=m_vec.size();
if (sz<2) return; 
IdxTy nonzed=sz;
for (IdxTy i=0; i<sz; ++i)
{
const IdxTy j=sz-1-i;
Tpoly &  p=m_vec[j];
mi.trim_polynomial(p);
if (is_zero(p)) if (nonzed==(j+1)) --nonzed;


}
if (nonzed!=sz) m_vec.resize(nonzed);

}
// this is only really good for when P(0)==0 
void div_exp()
{
if (m_vec.size()!=0) m_vec.erase(m_vec.begin()); 

}
IdxTy  lowzed() const 
{
const IdxTy sz=m_vec.size();
const IdxTy zed=0;
for (IdxTy i=0; i<sz; ++i)
{
if (!is_zero(m_vec[i])) return zed;
++zed;
}
}


value_type min_coef() const { return extreme_coef(0); }
value_type max_coef() const { return extreme_coef(1); }

// hopefully this casts ok lol 
IdxTy min_xcoef() const { return IdxTy(extreme_xcoef(3)); }


IdxTy extreme_xcoef(const IdxTy flag) const
{
const IdxTy sz=m_vec.size();
IdxTy extre=0;
bool have_value=false;
for (IdxTy i=0; i<sz; ++i)
{
const Tpoly & p=m_vec[i];
const IdxTy szp=p.size();
IdxTy szpi=(have_value)?extre:szp;
if (szpi>szp) szpi=szp;
for (IdxTy j=0; j<szpi; ++j)
{
bool hits=(p[j]!=0);
if ((j==0)&&(hits)) return 0; // no common x factors 
if (hits) if ((j<=extre)||(extre==0) ) {extre=j; have_value=true;  break; } 
}
}	 
return extre;
}

value_type extreme_coef(const IdxTy flag) const
{
const IdxTy sz=m_vec.size();
value_type extre=0;
bool have_value=false;
for (IdxTy i=0; i<sz; ++i)
{
const Tpoly & p=m_vec[i];
const IdxTy szp=p.size();
for (IdxTy j=0; j<szp; ++j)
{
bool hits=false;
switch (flag)
{
case 0: { hits=(p[j]<extre); break; } 
case 1: { hits=(p[j]>extre); break; } 
//case 3: { hits=(p[j]>extre); break; } 
default: { MM_ERR(" bad extreme_coef op "<<flag); } 
}

if ((hits)||!have_value) extre=p[j];
}	 

}

return extre;
}






void  scale_coefs(const IdxTy flag, const value_type & v=value_type(1))
{
const IdxTy sz=m_vec.size();
value_type extre=0;
for (IdxTy i=0; i<sz; ++i)
{
 Tpoly & p=m_vec[i];
const IdxTy szp=p.size();
for (IdxTy j=0; j<szp; ++j)
{
switch (flag)
{
case 0: { p[j]=p[j]*v; break; } 
case 1: { p[j]=p[j]/v; break; } 
case 2: { if (j!=0) p[j-1]=p[j]; break; } 
default: { MM_ERR(" bad extreme_coef op "<<flag); } 
}

} // j	 
switch (flag)
{
case 2: { p[szp-1]=0; break; } 
//default: { MM_ERR(" bad extreme_coef op "<<flag); } 
}


}
if (flag==2) trim();
}




/*
Myt operator*(const value_type & that) //const
{
Myt x;
Tpoly x1;
x1.push_back(that);
mi.multiply_polynomials(x.m_Q,m_Q,x1); // zero E numerator
mi.multiply_polynomials(x.m_P,m_P,x1); // zero E numerator
//MM_MSG(MMPR3(mi.print_polynomial(m_Q),mi.print_polynomial(x1),mi.print_polynomial(x.m_Q)))
//MM_MSG(MMPR3(mi.print_polynomial(m_P),mi.print_polynomial(x1),mi.print_polynomial(x.m_P)))

x.m_R=m_R;
x.m_S=m_S;
mi.trim_polynomial(x.m_P);
mi.trim_polynomial(x.m_Q);
mi.trim_polynomial(x.m_R);
mi.trim_polynomial(x.m_S);
return x; 
}
*/
// derivative and put the higher powers into e2 for top and bottom 
/*
TODO FIXME FUDD see if this should work cancelling exp terms lol 
mjm_burmann.h350  bcc.to_string()=(1)*exp(-x*x)+0
----------------------------------------------------------------
(-2x)*exp(-x*x)+

mjm_burmann.h352  f.der(bcc).to_string()=()*exp(-x*x)+
----------------------------------------------------------------
()*exp(-x*x)+

*/

void add_single_term(const IdxTy  expn, const IdxTy p, const IdxTy flag, const value_type & c)
{
Tpoly po=Tpoly();
while (m_vec.size()<=expn) m_vec.push_back(po); 
Tpoly & pi=m_vec[expn];
while (pi.size()<=p) pi.push_back(0); 
pi[p]+=c;

}

 Myt  der( )
{
Myt y;
der(y);

return y; 
}
 void der( Myt & y)
{
y.clear();
const IdxTy sz=m_vec.size();
for (IdxTy i=0; i<sz; ++i)
{
Tpoly p;
mi.next_gauss_der_beta_polynomial(p,m_vec[i],m_beta*i);
y.m_vec.push_back(p);

}

}

 Tv evaluate( const Tv & x) {return evaluate(x,exp(-m_beta*x*x)); }
 Tv evaluate( const Tv & x, const Tv & normv)
{
Tv res=Tv(0);
Tv expv=Tv(1);
const IdxTy sz=m_vec.size();
for (IdxTy i=0; i<sz; ++i)
{
const Tpoly & p=m_vec[i];
res+=mi.evaluate_polynomial(p,x)*expv;
expv=expv*normv;
}
return res;
}


/* 
 void der( Myt & y, Tpoly & e2n, Tpoly & e2d)
{
// keep within this constraint, punt the higher order polys
// (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x)) ->
// (N'D- D'N)/D^2  from quotient rule.
// N' = P'E-2xPE+Q', D'=R'E-2xRE+S'
Tpoly n1,d1,n2,d2,t1,t2,t3,t4;
mi.next_gauss_der_polynomial(n1,m_P);
mi.next_gauss_der_polynomial(d1,m_R);
mi.differentiate_polynomial(n2,m_Q);
mi.differentiate_polynomial(d2,m_S);
mi.multiply_polynomials(y.m_S,m_S,m_S);
mi.multiply_polynomials(t1,m_R,m_S);
mi.add_polynomials(y.m_R,t1,t1,1.0,1.0);
mi.multiply_polynomials(e2d,m_R,m_R);
// N' is n1*E+n2 and D' is d1*E+d2
// N = PE+Q, D=RE+S
// N'D= E(R*n2+S*n1)+n2*S+EER*n1 
// ND'= EE*(d1*P) + E*(d1*Q+d2*P) +Q*d2  
// N'D- ND'
mi.multiply_polynomials(t1,d1,m_P);
mi.multiply_polynomials(t2,n1,m_R);
mi.add_polynomials(e2n,t2,t1,1.0,-1.0);
mi.multiply_polynomials(t1,d1,m_Q);
mi.multiply_polynomials(t2,d2,m_P);
mi.add_polynomials(t3,t1,t2); //  ND' E term 
mi.multiply_polynomials(t1,n1,m_S);
mi.multiply_polynomials(t2,n2,m_R);
mi.add_polynomials(t4,t1,t2); //  N'D E term 
mi.add_polynomials(y.m_P,t3,t4,-1.0,1.0); //  N'D E term 

mi.multiply_polynomials(t1,d2,m_Q);
//mi.multiply_polynomials(y.m_Q,d2,m_Q);
mi.multiply_polynomials(t2,n2,m_S);
mi.add_polynomials(y.m_Q,t1,t2,-1.0,1.0); //  N'D E term 
if ( is_zero(y.m_S) && is_zero(y.m_Q))
{
y.m_S= y.m_R; y.m_Q=y.m_P;
y.m_R=e2d; y.m_P=e2n;
e2d.clear();
e2n.clear();
}
if ( !is_zero(e2n) ||!is_zero(e2d))
{
MM_MSG(" der genearated e 2 terms "<<MMPR2(e2n.to_string(),e2d.to_string()))
}
//template<class Td>  void polynomial_list_reduce(Td & dest, const Td & v1)
std::vector<Tpoly> dest,v1;
v1.push_back(y.m_P);
v1.push_back(y.m_Q);
v1.push_back(y.m_R);
v1.push_back(y.m_S);
mi.polynomial_list_reduce( dest,  v1);
y.m_P=dest[0];
y.m_Q=dest[1];
y.m_R=dest[2];
y.m_S=dest[3];

}
*/

// this is going to get evaluated at zero's of the denominator,
// however l'hopital can save it 
// (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x)) ->
/*
Tv eval(const Tpoly & p, const Tpoly & q, const Tv & x, const Tv & ex)
{
return mi.evaluate_polynomial(p,x)*ex+
			mi.evaluate_polynomial(q,x);
}
Tv evaluate(const Tv & x) //const
{
const Tv ex=exp(-x*x);
return evaluate(x,ex);
}

Tv evaluate(const Tv & x, const Tv & ex) //const
{
const Tv n=eval(m_P,m_Q,x,ex); 
const Tv d=eval(m_R,m_S,x,ex); 
// ok this is not right lol as small non zero prevent l'hopital 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( d!=0) return (n/d);
Tpoly P,Q,R,S;
Tpoly P2,Q2,R2,S2;
P=m_P; Q=m_Q; R=m_R; S=m_S;
P2=P; Q2=Q; R2=R; S2=S;

while (true)
{
mi.next_gauss_der_polynomial(P2,P);
mi.differentiate_polynomial(Q2,Q);
mi.next_gauss_der_polynomial(R2,R);
mi.differentiate_polynomial(S2,S);
const Tv n=eval(P2,Q2,x,ex); 
const Tv d=eval(R2,S2,x,ex); 
//	MM_MSG(MMPR(print(P,Q,R,S)))	
//	MM_MSG(MMPR4(print(P2,Q2,R2,S2),x,n,d))	
if (d!=0)  return (n/d);
// TODO decide if you wan to do this FIXME 
//if ((n!=0)||(d!=0)) return (n/d);
// return the nan apropos for Tv 
// not sure cycles are possible with the exp poly 
if ((P2==P) &&(Q2==Q) && (R2==R)&&(S2==S)) return (n/d);
P=P2; Q=Q2; R=R2; S=S2;
} // while 
// not reachable 
return 0; 

} // evaluate
*/
/*
Tv evaluate_at_zero() //const
{
// FIXME TODO this will bomb if the sizes are zero... 
const IdxTy szp=m_P.size();
const IdxTy szq=m_Q.size();
const IdxTy szr=m_R.size();
const IdxTy szs=m_S.size();
const Tv p0=(szp>0)?m_P[0]:0;
const Tv q0=(szq>0)?m_Q[0]:0;
const Tv r0=(szr>0)?m_R[0]:0;
const Tv s0=(szs>0)?m_S[0]:0;
const Tv n=p0+q0;  
const Tv d=r0+s0; // eval(m_R,m_S,x,ex); 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( d!=0) return (n/d);
//Tpoly expn,expd,temp,neff,deff;
// TODO FIXME need to find real limits ... 
// int nterms=int(szq)-int(szp);
//	if (nterms<10) nterms=10;
// int dterms=int(szs)-int(szr);
//	if (dterms<10) dterms=10;
//mi.gauss_taylor_polynomial(expn, nterms);
//mi.gauss_taylor_polynomial(expd, dterms);
//MM_MSG(MMPR(mi.print_polynomial(expn)))
//MM_MSG(MMPR(mi.print_polynomial(expd)))
//mi.multiply_polynomials(temp,expn,m_P);
//mi.add_polynomials(neff,temp,m_Q);
//mi.multiply_polynomials(temp,expd,m_R);
//mi.add_polynomials(deff,temp,m_S);
// TODO FIXME KLUGE fix the limit lol 
const IdxTy szdeff=szp+szq+10 ; // neff.size();
const IdxTy szneff=szr+szs+10  ; // deff.size();

IdxTy i=0;
while (true)
{
Tv n=0;
Tv d=0; 
if (i>=szneff) if (i>=szdeff)
{
MM_ERR(" no limit found "<<MMPR4(szp,szq,szr,szs)<<MMPR4(szneff,szdeff,mi.print_polynomial(m_P),mi.print_polynomial(m_R)))
 return (n/d);
}
d=mi.gauss_term_polynomial( m_R,i );
const Tv si=(szs>i)?m_S[i]:0;
d+=si;
// this is not really right as 1/0 is a problem however that should
// not occur but the user needs to be careful.
// Danger Will Robinson !!!!!! FIXED 
// and it is no faster anyway 
//if (d!=0)
{
n=mi.gauss_term_polynomial( m_P,i );
const Tv qi=(szq>i)?m_Q[i]:0;
n+=qi;
}

//if (i<szneff) n=neff[i];
//if (i<szdeff) d=deff[i];
if (d!=0) return (n/d);
++i;
} // true 
MM_ERR(" should not get here ")
return 0; // unreadable
} // evaluate_at_zero
*/

/*StrTy print(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S)
{
SsTy ss;
const StrTy exps="exp(-x*x)";
ss<<"("<<mi.print_polynomial(P)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(Q)<<CRLF;
ss<<"----------------------------------------------------------------"<<CRLF;
ss<<"("<<mi.print_polynomial(R)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(S)<<CRLF;

return ss.str();
}
*/

void clear() { m_vec.clear(); } 
//StrTy to_string()  { return print(m_P,m_Q,m_R,m_S); }

StrTy exp_string(const IdxTy n )const   
{
SsTy ss;
if (n==0) return ss.str(); 
if (n==1) ss<<"exp(-("<<m_beta<<")*x^2)";
else ss<<"exp(-("<<n<<"*"<<m_beta<<")*x^2)";

return ss.str();
}
StrTy to_string() const   
{ 
SsTy ss;
const IdxTy sz=m_vec.size();
if (sz==0) { ss<<"0"; return ss.str(); }
IdxTy padd=0;
for (IdxTy i=0; i<sz; ++i)
{
const IdxTy idx=sz-1-i;
Tpoly  p=m_vec[idx];
if (is_zero(p)) continue;
const StrTy e=exp_string(idx);
if (padd!=0) ss<<"+";
ss<<"("<<mi.print_polynomial(p)<<")"<<e;
++padd;
}
if (padd==0) { ss<<"0"; return ss.str(); }
return ss.str();
 }

private:

void  Init()
{
m_vec.clear();
//Mono(m_default,0,1);
//return m_default;
}
void Default() { 
//Default(m_P); Default(m_Q); Default(m_R); Default(m_S);
 }
void Mono(Tpoly & p, const IdxTy n, const Tv & v)
{
//p = Tpoly(n+1);
//p[n]=v;
}

//void Default(Tpoly & x) { x=m_default; }
PolyVec m_vec;
D m_beta;
//Tpoly m_default;
// rational function (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x))
//Tpoly m_P,m_Q,m_R,m_S;
// TODO wtf 
mutable Mi mi;

}; // poly_norm
////////////////////////////////////////////////////////////////////////

template <class Tv=double> class exp_poly_ratio
{

class Tr{
public:
//typedef double D;
typedef Tv  D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr
public:
typedef typename Tr::D D;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::StrTy StrTy;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_block_matrix<IdxTy> MyBlockInt;
typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
typedef mjm_generic_itor<2> LoopItor;

typedef exp_poly_ratio<Tv>  Myt;
typedef Tv value_type;
typedef mjm_integrals Mi;
typedef shape_function_integrals SFi;

typedef poly_norm<Tv> PolyPoly;
typedef SFi::simple_polynomial<Tv> Tpoly;
typedef std::vector<Tpoly> PolyVec;
public:
exp_poly_ratio() { Init(); Default(); }
//poly_norm(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S) 
//:m_P(P),m_Q(Q),m_R(R),m_S(S)
//{ Init();  }
exp_poly_ratio(const value_type & x) {  from_const(x); } 
exp_poly_ratio(const StrTy & x) {  from_string(x); } 

void from_const(const value_type & x)
{
//Mono(m_P,0,0); Mono(m_Q,0,x); Mono(m_R,0,0); Mono(m_S,0,1);
}
void poly(const Tpoly & p)
{
m_n.clear();
m_d.clear();
m_n.push_back(p);
Tpoly one;
one.push_back(1);
m_d.push_back(one);
}
void poly(const Tpoly & pn, const Tpoly & pd)
{
m_n.clear();
m_d.clear();
m_n.push_back(pn);
m_d.push_back(pd);
}


void from_string(const StrTy & x)
{
Tpoly one, zero,expp;
zero.push_back(0);
one.push_back(1);
expp.push_back(0);
expp.push_back(1);
// this really needs an expression parser
if (x=="one") { m_n.push_back(one); m_d.push_back(one); return ; } 
if (x=="exp") { m_n.push_back(zero); m_n.push_back(one); m_d.push_back(one); return ; } 


}

void normal()
{
//Mono(m_P,0,1); Mono(m_Q,0,0); Mono(m_R,0,0); Mono(m_S,0,1);
Tpoly one, zero;
zero.push_back(0);
one.push_back(1);
 m_n.push_back(zero); m_n.push_back(one); m_d.push_back(one); 


}
bool is_zero(const Tpoly & x)
{
const IdxTy sz=x.size();
for (IdxTy i=0; i<sz; ++i) if (x[i]!=0) return false;
return true;

}
void add_single_term(const IdxTy  expn, const IdxTy p, const IdxTy flag, const value_type & c)
{
PolyPoly & x=((flag&1)==0)?m_n:m_d;
x.add_single_term(expn,p,flag,c);
}

Myt operator+(const Myt & that) const
{
Myt x;
if (m_d==that.m_d)
{
x.m_d=m_d;
x.m_n=m_n+that.m_n;
return x; 
}
x.m_n=m_n*that.m_d+ that.m_n*m_d;
x.m_d=m_d*that.m_d;
return x; 
}
Myt operator-(const Myt & that) const
{
Myt x;
if (m_d==that.m_d)
{
x.m_d=m_d;
x.m_n=m_n-that.m_n;
return x; 
}
x.m_n=m_n*that.m_d- that.m_n*m_d;
x.m_d=m_d*that.m_d;
return x; 
}






Myt operator*(const Myt & that) const
{
Myt x;
x.m_n=m_n*that.m_n;
x.m_d=m_d*that.m_d;
return x; 
}



Myt operator*(const value_type & that) const
{
Myt x;
x.m_n=m_n*that;
x.m_d=m_d;
return x; 
}




Myt operator/(const Myt & that) const
{
Myt x;
if (m_d==that.m_d)
{
x.m_n=m_n;
x.m_d=that.m_n;
return x;
}
x.m_n=m_n*that.m_d;
x.m_d=m_d*that.m_n;
x.reduce();
//MM_MSG(MMPR4(x.m_P.size(),x.m_Q.size(),x.m_R.size(),x.m_S.size())<<MMPR(x.to_string()))

return x;
}

// derivative and put the higher powers into e2 for top and bottom 
/*
TODO FIXME FUDD see if this should work cancelling exp terms lol 
mjm_burmann.h350  bcc.to_string()=(1)*exp(-x*x)+0
----------------------------------------------------------------
(-2x)*exp(-x*x)+

mjm_burmann.h352  f.der(bcc).to_string()=()*exp(-x*x)+
----------------------------------------------------------------
()*exp(-x*x)+

*/

Myt der() //const
{Myt y;

der(y);
return y;

}

 void der( Myt & y)
{
// the ratio of the PolyPoly types should be cloased now
PolyPoly nprime,dprime;
nprime=m_n.der();
dprime=m_d.der();
y.m_n=m_d*nprime-m_n*dprime;
y.m_d=m_d*m_d;
y.reduce();
}
void reduce()
{
reduce_exp();
reduce_scale();
reduce_x();
}
void reduce_x()
{
IdxTy  minn=m_n.min_xcoef();
IdxTy  mind=m_d.min_xcoef();
IdxTy cnt=(minn>mind)?mind:minn;
MM_MSG("REDUCE "<<MMPR3(cnt,minn,mind))
while (cnt!=0)
{
m_n.scale_coefs(2);
m_d.scale_coefs(2);
--cnt; 
} 


}
void reduce_scale()
{
// this need order of magniute abs etc 
auto maxn=m_n.max_coef();
auto maxd=m_d.max_coef();
auto minn=m_n.min_coef();
auto mind=m_d.min_coef();
auto f=maxn;
if (f!=0)
{
	MM_MSG(" scaling "<<MMPR4(maxn,maxd,minn,mind)<<MMPR(f))
	m_n.scale_coefs(1,f);
	m_d.scale_coefs(1,f);
}



}


// remove the common exp terms
void reduce_exp()
{
const IdxTy szn=m_n.size();
const IdxTy szd=m_d.size();
const IdxTy szmin=(szn<szd)?szn:szd;
IdxTy comfac=0;
for ( IdxTy i=0; i<szmin; ++i)
{
if (is_zero(m_n[i])&&is_zero(m_d[i])) { ++comfac; } else {break;}

}
while (comfac>0)
{
m_n.div_exp();
m_d.div_exp();

--comfac;
}

}






// this is going to get evaluated at zero's of the denominator,
// however l'hopital can save it 
// (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x)) ->
/*
Tv eval(const Tpoly & p, const Tpoly & q, const Tv & x, const Tv & ex)
{
return mi.evaluate_polynomial(p,x)*ex+
			mi.evaluate_polynomial(q,x);
}
*/
Tv evaluate(const Tv & x) //const
{
const Tv ex=exp(-x*x);
return evaluate(x,ex);
}

Tv evaluate(const Tv & x, const Tv & ex) //const
{
 Tv ne=m_n.evaluate(x,ex); 
 Tv de=m_d.evaluate(x,ex); 
//MM_MSG(MMPR4(x,ex,ne,de)<<MMPR2(m_n.to_string(),m_d.to_string()))
//const Tv d=eval(m_R,m_S,x,ex); 
// ok this is not right lol as small non zero prevent l'hopital 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( de!=0) return (ne/de);
PolyPoly nprime=m_n;
PolyPoly dprime=m_d;
PolyPoly n=m_n;
PolyPoly d=m_d;

while (true)
{
n.der(nprime);
d.der(dprime);
const Tv dprimev=dprime.evaluate(x,ex);
const Tv nprimev=nprime.evaluate(x,ex);
//	MM_MSG(MMPR(print(P,Q,R,S)))	
//	MM_MSG(MMPR4(print(P2,Q2,R2,S2),x,n,d))	
if (dprimev!=0)  return (nprimev/dprimev);
// TODO decide if you wan to do this FIXME 
//if ((n!=0)||(d!=0)) return (n/d);
// return the nan apropos for Tv 
// not sure cycles are possible with the exp poly 
//if ((P2==P) &&(Q2==Q) && (R2==R)&&(S2==S)) return (n/d);
if ((nprime==n)&&(dprime==d)) return (nprimev/dprimev);
n=nprime;
d=dprime;
//P=P2; Q=Q2; R=R2; S=S2;
} // while 
// not reachable 
return 0; 

} // evaluate

Tv evaluate_at_zero() //const
{
Tv n=m_n.evaluate_at_zero();
Tv d=m_d.evaluate_at_zero();
if ( d!=0) return (n/d);

/*
// FIXME TODO this will bomb if the sizes are zero... 
const IdxTy szp=m_P.size();
const IdxTy szq=m_Q.size();
const IdxTy szr=m_R.size();
const IdxTy szs=m_S.size();
const Tv p0=(szp>0)?m_P[0]:0;
const Tv q0=(szq>0)?m_Q[0]:0;
const Tv r0=(szr>0)?m_R[0]:0;
const Tv s0=(szs>0)?m_S[0]:0;
const Tv n=p0+q0;  
const Tv d=r0+s0; // eval(m_R,m_S,x,ex); 
//		MM_MSG(MMPR4(print(m_P,m_Q,m_R,m_S),x,n,d))	
if ( d!=0) return (n/d);
// TODO FIXME KLUGE fix the limit lol 
const IdxTy szdeff=szp+szq+10 ; // neff.size();
const IdxTy szneff=szr+szs+10  ; // deff.size();

IdxTy i=0;
while (true)
{
Tv n=0;
Tv d=0; 
if (i>=szneff) if (i>=szdeff)
{
MM_ERR(" no limit found "<<MMPR4(szp,szq,szr,szs)<<MMPR4(szneff,szdeff,mi.print_polynomial(m_P),mi.print_polynomial(m_R)))
 return (n/d);
}
d=mi.gauss_term_polynomial( m_R,i );
const Tv si=(szs>i)?m_S[i]:0;
d+=si;
// this is not really right as 1/0 is a problem however that should
// not occur but the user needs to be careful.
// Danger Will Robinson !!!!!! FIXED 
// and it is no faster anyway 
//if (d!=0)
{
n=mi.gauss_term_polynomial( m_P,i );
const Tv qi=(szq>i)?m_Q[i]:0;
n+=qi;
}

//if (i<szneff) n=neff[i];
//if (i<szdeff) d=deff[i];
if (d!=0) return (n/d);
++i;
} // true 
MM_ERR(" should not get here ")

*/

return 0; // unreadable
} // evaluate_at_zero

/*StrTy print(const Tpoly & P, const Tpoly & Q, const Tpoly & R, const Tpoly & S)
{
SsTy ss;
const StrTy exps="exp(-x*x)";
ss<<"("<<mi.print_polynomial(P)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(Q)<<CRLF;
ss<<"----------------------------------------------------------------"<<CRLF;
ss<<"("<<mi.print_polynomial(R)<<")*"<<exps;
ss<<"+"<<mi.print_polynomial(S)<<CRLF;

return ss.str();
}
*/

//StrTy to_string()  { return print(m_P,m_Q,m_R,m_S); }

StrTy to_string()const   
{ 
SsTy ss;
ss<<"("<<m_n.to_string()<<")/("<<m_d.to_string()<<")";
return ss.str();
}

private:

void  Init()
{
//m_vec.clear();
//Mono(m_default,0,1);
//return m_default;
}
void Default() { 
//Default(m_P); Default(m_Q); Default(m_R); Default(m_S);
 }
void Mono(Tpoly & p, const IdxTy n, const Tv & v)
{
//p = Tpoly(n+1);
//p[n]=v;
}

//void Default(Tpoly & x) { x=m_default; }
//PolyVec m_vec;
PolyPoly m_n,m_d;
//Tpoly m_default;
// rational function (P(x)exp(-x*x)+Q(x))/(R(x)exp(-x*x)+S(x))
//Tpoly m_P,m_Q,m_R,m_S;
Mi mi;

}; // exp_poly_ratio
////////////////////////////////////////////////////////////////////////
template <class St> void stack_op(St& sstack,const IdxTy op, const IdxTy p= 0 )
	{ 	const IdxTy sz=sstack.size(); 
		auto  x1=sstack[sz-1];
		auto  x2=sstack[sz-2];
		if ((p&1)!=0) 		
		{ sstack.pop_back(); sstack.pop_back(); }
		switch (op)
		{
		case 0:{sstack.push_back(x1+x2); break; }
		case 1: { sstack[sz-1]=x2; sstack[sz-2]=x1;	 break; } 
		case 2:{sstack.push_back(x1*x2); break; }

		default: MM_ERR(" bad case stack op "<<MMPR2(op,p))
		} //  switch 
	} 

void parse_mode(CommandInterpretter & li)
{
StrTy local_label="closed_families";
typedef exp_poly_ratio<> VarTy;
typedef VarTy::value_type Tvs;
std::vector<VarTy> sstack;
bool echo_cmd=false;
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
if (m_flp.log_commands()) 
		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz<1) continue;
const StrTy cmd=li.word(0);
const StrTy p1=(sz>1)?li.word(1):StrTy("");
if (echo_cmd) {MM_MSG(MMPR(li.line())) } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

if (cmd=="echo") { echo_cmd=true;  continue; } 
if (cmd=="no-echo") { echo_cmd=!true;  continue; } 
if (cmd=="push") { sstack.push_back(VarTy(p1));  continue; } 
if (cmd=="size") { MM_MSG(MMPR(sstack.size()));  continue; } 
if (cmd=="pop") { sstack.pop_back();  continue; } 
if (cmd=="dup") { sstack.push_back(sstack[sstack.size()-1]);  continue; } 
//template <class St> void stack_op(St& sstack,const IdxTy op, const IdxTy p  )
if (cmd=="swap") 
	{ 	stack_op(sstack,1);
//		const IdxTy sz=sstack.size();
//		VarTy x1=sstack[sz-1];
//		VarTy x2=sstack[sz-2];
//		sstack[sz-1]=x2; sstack[sz-2]=x1;	
		continue; 
	} 
if (cmd=="print") { MM_MSG(MMPR(sstack.back().to_string()));  continue; } 
if (cmd=="eval") { MM_MSG(MMPR(sstack.back().evaluate(Tvs(atof(p1.c_str())))));  continue; } 
if (cmd=="stack") { for (IdxTy i=0; i<sstack.size(); ++i ) { MM_MSG(MMPR2(i,sstack[i].to_string()));}  continue; } 
if (cmd=="add") 
	{ //	const IdxTy sz=sstack.size(); 
//		VarTy  x1=sstack[sz-1];
//		VarTy  x2=sstack[sz-2];
//		sstack.pop_back(); sstack.pop_back();
//		sstack.push_back(x1+x2);  
	 	stack_op(sstack,0);
		continue; 
	} 
if (cmd=="mult") 
	{ 	stack_op(sstack,2);
//const IdxTy sz=sstack.size(); sstack.push_back(sstack[sz-1]*sstack[sz-2]);  

	continue; } 
// TODO FIXME Danger Will robinson the push_Back is not executed if nested
// corruption ???????????????
if (cmd=="der") 
	{ const IdxTy sz=sstack.size(); 
		VarTy y=sstack[sz-1].der();
 		sstack.push_back(y);  
	continue; } 
if (cmd=="div") 
	{ const IdxTy sz=sstack.size(); sstack.push_back(sstack[sz-1]/sstack[sz-2]);  continue; } 

//void add_single_term(const IdxTy  expn, const IdxTy p, const IdxTy flag, const value_type & c)
if (cmd=="set-term") { if (!li.cmd_ok(5)) continue;
const IdxTy expn= atoi(li.word(1).c_str());
const IdxTy p =atoi(li.word(2).c_str());
const IdxTy flag =atoi(li.word(3).c_str());
const D  c =atof(li.word(4).c_str());
MM_MSG(MMPR4(expn,p,flag,c))
sstack.back().add_single_term(expn,p,flag,c);
continue;
 } 


//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
if (cmd=="quit") { MM_MSG(" exit parser ")  return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized parse command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad parse command ",li.line(),"NOTDATA");

}
} //parse_mode













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

}; //mjm_closed_families



/////////////////////////////////////////////////////////

#ifdef  TEST_CLOSED_FAMILIES__
int main(int argc,char **args)
{
typedef mjm_closed_families  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

