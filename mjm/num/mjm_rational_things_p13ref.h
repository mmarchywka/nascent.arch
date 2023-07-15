#ifndef MJM_SHOT_FUDD_H__
#define MJM_SHOT_FUDD_H__
 
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
//#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
#include "../mjm_fasta_ii.h"


//#include "mjm_collections.h"
//#include "mjm_svg_writer.h"


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
/*

g++ -DTEST_RATIONAL_THINGS__ -std=gnu++11 -DUSE_RATIONALS -I gmp/gmp-6.1.2 -Lgmp/gmp-6.1.2/.libs -lgmp -lreadline -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -O0 -gdwarf-3 -I.. -I../json/rapidjson-master/include/ -Wall -x c++ mjm_rational_things.h -o mjm_rational_things.out



g++ -DTEST_RATIONAL_THINGS__ -std=gnu++11 -DUSE_RATIONALS -I gmp/gmp-6.1.2 -Lgmp/gmp-6.1.2/.libs -lgmp -lreadline -O3 -I.. -I../json/rapidjson-master/include/ -Wall -x c++ mjm_rational_things.h -o mjm_rational_things.out


*/


////////////////////////////////////////////////////////////////

class rational_things_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
rational_things_params( const StrTy & nm) : Super(nm) {}
rational_things_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
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

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace rational_things_traits
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
typedef mjm_block_matrix<IdxTy> MyUBlock;
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


}; // trees_and_tables_traits
///////////////////////////////////////////////////////////////


class mjm_rational_things 
{
typedef rational_things_traits::Tr  Tr;
typedef mjm_rational_things Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MyUBlock  MyUBlock;

typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;


typedef rational_things_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;
public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
mjm_rational_things():m_dmel(new Dmel()) {Init();}
mjm_rational_things(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_rational_things()
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



// with lt[0]=0, lt[i] = sum_{i==j} ( log(j))
void log_table( D * lt, const IdxTy sz, const bool summ=true) const
{
lt[0]=0;
for(IdxTy i=1; i<sz; ++i)
{
if (summ) lt[i]=lt[i-1]+log(D(i));
else  lt[i]=log(D(i));
}
}
// log(N-i)
void slog_table( D * lt, const IdxTy sz ) const
{
for(IdxTy i=0; i<=sz; ++i)
{
lt[i]=log(D(sz-i));
}
}





////////////////////////////////////////////////////

typedef mjm_cursor_pair_itor CpiTy; // (m_i, l,  highest);
//////////////////////////////////////////////////////////////


void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
if (tabsep) m_ragged_map[name].sep("\t");
if (load_parsed) m_ragged_map[name].load(fn);
else m_ragged_map[name].load_lines(fn);

MM_ERR(MMPR2(m_ragged_map[name].size(),name))
}



// this needs a rational version 
D binom(const IdxTy top, const IdxTy bot,const D * logtable) const
{
D fudd= exp(logtable[top]-logtable[bot]-logtable[(top-bot)]);
//MM_ERR(MMPR4(fudd,top,bot,exp(logtable[top])))
return fudd;
}

#ifdef USE_RATIONALS
typedef mjm_integrals Mi;
typedef mjm_rational RatTy;
typedef std::vector<RatTy> RatCo;
mutable std::vector<RatTy> m_fact_table;



void make_nest(RatCo & d, const IdxTy n ) const
{
Mi mi;
RatCo nest;
nest.push_back(1);
for(IdxTy i=1; i<=n; ++i) { nesting(nest); } //i 
RatCo x=nest;
// missing const term of zero 
x.push_back(0);
mi.reverse_polynomial(d,x);
} // sum_nest
//////////////////////////////////////////////////

// mjm_strict_triangle_itor(const LocTy & lengths ) :  m_lengths(lengths)

void recurse_triangle(int * mvec, const LocTy & lv, const int &  ilim
, const IdxTy & cnt, const IdxTy depth, const IdxTy start, const int l)
const
{
for(IdxTy j=(start); j<cnt; ++j)
{
const int l2=l+lv[j];
if (int(l2+depth)>ilim) return; 
const bool dp=((depth&1)!=0);
if (dp) --mvec[l2+depth]; else ++mvec[l2+depth];
recurse_triangle(mvec,lv,ilim,cnt,depth+1,j+1,l2);
} // j 

}
///////////////////////////////////////////////////////////////////

//typedef long long int big_coef_type;
typedef long int big_coef_type;

big_coef_type * make_F_table(const IdxTy cnt, const IdxTy sz) const
{
Mi mi;
const IdxTy lvl=cnt-2;
RatCo tp;
// this is the only function needed
make_nest(tp,lvl);
// now make a table 
const IdxTy table_sz=sz; // ilim+2;
big_coef_type * pctable= new big_coef_type[table_sz];
for(IdxTy i=0; i<table_sz; ++i)
{
RatTy ni=mi.evaluate_polynomial(tp,RatTy(i));
ni.simplify();
if (big_coef_type(ni.get_d()) !=1) 
{MM_ERR(" njot int "<<MMPR3(ni,i,lvl))   ; delete [] pctable; return 0; }
pctable[i]=big_coef_type(ni.get_n());
if (pctable[i]<0)
{
MM_ERR(" possible overflow wrap "<<MMPR4(i,cnt,sz,pctable[i]))
}
} // i 
return pctable;
} // make_F_table



void make_multiplicities_table_lens(LocTy &  m, const LocTy & lv
, const big_coef_type * _ctable=0, const IdxTy flags=0) const
{
const bool print_res=((flags&1)==0); 
const bool print_tab=((flags&2)==0); 
Mi mi;
const IdxTy cnt=lv.size();
if (cnt<1) return; 
IdxTy ns=0;
const IdxTy ltsz=cnt;
for(IdxTy i=0; i<ltsz; ++i)  { ns+=lv[i]; } 
const int  ilim=ns>>1;
int mvec[ilim+1];
memset(&mvec[0],0,(ilim+1)*sizeof(int));
recurse_triangle(mvec,lv,ilim,cnt,0,0,0);
const IdxTy lvl=cnt-2;

const big_coef_type * ctable=_ctable;
if (ctable==0) { 
ctable= make_F_table(cnt, ilim+2);
#if 0  
RatCo tp;
// this is the only function needed
make_nest(tp,lvl);
// now make a table 
const IdxTy table_sz=ilim+2;
big_coef_type * pctable= new big_coef_type[table_sz];
for(IdxTy i=0; i<table_sz; ++i)
{
RatTy ni=mi.evaluate_polynomial(tp,RatTy(i));
ni.simplify();
if (big_coef_type(ni.get_d()) !=1) {MM_ERR(" njot int "<<MMPR3(ni,i,lvl))}
pctable[i]=big_coef_type(ni.get_n());
} // i 
ctable=pctable;
#endif

}
std::vector<int> pos,npos;
//IdxTy ptr=0;
//IdxTy tptr=1;
int i=0;
m=LocTy(ns+1);
if (print_tab) { MM_ERR(MMPR4(m.size(),ns,lvl,ilim)) } 
for(; i<=ilim; ++i)  
{ 
const IdxTy ie=i+1;
big_coef_type ni=ctable[ie];
//RatTy ni=mi.evaluate_polynomial(tp,RatTy(ie));
for(IdxTy k=0; k<pos.size(); ++k)  
// ni-=RatTy(npos[k])*mi.evaluate_polynomial(tp,RatTy(ie-pos[k]));
ni-=npos[k]*ctable[ie-pos[k]];
if ( print_res) { MM_ERR(MMPR2(i,ni)) } 
m[i]=ni; // (IdxTy(D(ni)));
m[ns-i]=m[i];
if ( mvec[i]!=0) { pos.push_back(ie); npos.push_back(mvec[i]);  }  
} // i  
if (_ctable==0) {  delete [] ctable; } 
} 


///////////////////////////////////////////////////////////////////////
void make_multiplicities_with_lens(LocTy &  m, const LocTy & lv, const IdxTy flags=0) const
{
const bool print_res=((flags&1)==0); 
Mi mi;
const IdxTy cnt=lv.size();
if (cnt<1) return; 
IdxTy ns=0;
const IdxTy ltsz=cnt;
for(IdxTy i=0; i<ltsz; ++i)  { ns+=lv[i]; } 
const int  ilim=ns>>1;
int mvec[ilim+1];
memset(&mvec[0],0,(ilim+1)*sizeof(int));
recurse_triangle(mvec,lv,ilim,cnt,0,0,0);

const IdxTy lvl=cnt-2;
RatCo tp;
make_nest(tp,lvl);
std::vector<int> pos,npos;
//IdxTy ptr=0;
//IdxTy tptr=1;
int i=0;
m=LocTy(ns+1);
MM_ERR(MMPR4(m.size(),ns,lvl,ilim))
for(; i<=ilim; ++i)  
{ 
const IdxTy ie=i+1;
RatTy ni=mi.evaluate_polynomial(tp,RatTy(ie));
for(IdxTy k=0; k<pos.size(); ++k)  
 ni-=RatTy(npos[k])*mi.evaluate_polynomial(tp,RatTy(ie-pos[k]));
if ( print_res) { MM_ERR(MMPR2(i,ni)) } 
m[i]=(IdxTy(D(ni)));
m[ns-i]=m[i];
if ( mvec[i]!=0) { 
//MM_ERR(" fudd "<<MMPR2(i,mvec[i]))

pos.push_back(ie); npos.push_back(mvec[i]);  }  
//m.push_back(ni);
} // i  
} 



//////////////////////////////////////////////////////////////////

void make_multiplicities_with_less_old_lens(LocTy &  m, const LocTy & lv, const IdxTy flags=0) const
{
const bool print_res=((flags&1)==0); 
Mi mi;
const IdxTy cnt=lv.size();
if (cnt<1) return; 
IdxTy ns=0;
const IdxTy ltsz=cnt;
for(IdxTy i=0; i<ltsz; ++i)  { ns+=lv[i]; } 
const int  ilim=ns>>1;
int mvec[ilim+1];
memset(&mvec[0],0,(ilim+1)*sizeof(int));
//for(int  i=0; i<=ilim;  ++i) mvec[i]=0;
for(IdxTy i=0; i<cnt; ++i)
{
const int l=lv[i];
if (l>ilim) break; 
++mvec[l]; // =x+1;
for(IdxTy j=(i+1); j<cnt; ++j)
{
const int l2=l+lv[j];
if (l2>ilim) break; 
//MM_ERR(MMPR4(l2,ilim,l,mvec[l])<<MMPR(mvec[l2]))
--mvec[l2+1];
//mvec[l2]=mvec[l2]-1;
for(IdxTy k=(j+1); k<cnt; ++k)
{
const int l3=l2+lv[k];
if (l3>ilim) break; 
//MM_ERR(MMPR4(l2,ilim,l,mvec[l])<<MMPR(mvec[l2]))
++mvec[l3+2];
//mvec[l2]=mvec[l2]-1;


} // k 


} // j
} // i 
// mjm_strict_triangle_itor sti(lv) ;
int mvecsti[ilim+1];
memset(&mvecsti[0],0,(ilim+1)*sizeof(int));
recurse_triangle(mvecsti,lv,ilim,cnt,0,0,0);
//sti.make_mult_array(mvecsti);
for(int i=0; i<=ilim;  ++i) {int d= mvec[i]-mvecsti[i];
MM_ERR(MMPR4(i,mvec[i],mvecsti[i],d))}
//for(IdxTy i=0; i<(cnt-1); ++i)  { l=lv[i]; ns+=l; lt[i+1]=l+lt[i]; } 
for(int i=0; i<=ilim;  ++i) { mvec[i]=mvecsti[i]; } 

const IdxTy lvl=cnt-2;
RatCo tp;
make_nest(tp,lvl);
std::vector<int> pos,npos;
//IdxTy ptr=0;
//IdxTy tptr=1;
int i=0;
m=LocTy(ns+1);
for(; i<=ilim; ++i)  
{ 
const IdxTy ie=i+1;
RatTy ni=mi.evaluate_polynomial(tp,RatTy(ie));
for(IdxTy k=0; k<pos.size(); ++k)  
 ni-=RatTy(npos[k])*mi.evaluate_polynomial(tp,RatTy(ie-pos[k]));
if ( print_res) { MM_ERR(MMPR2(i,ni)) } 
m[i]=(IdxTy(D(ni)));
m[ns-i]=m[i];
if ( mvec[i]!=0) { 
MM_ERR(" fudd "<<MMPR2(i,mvec[i]))

pos.push_back(ie); npos.push_back(mvec[i]);  }  
//m.push_back(ni);
} // i  
} 


/////////////////////////////////////////////////

void make_multiplicities_with_lens_old(LocTy &  m, const LocTy & lv, const IdxTy flags=0) const
{
const bool print_res=((flags&1)==0); 
Mi mi;
const IdxTy cnt=lv.size();
if (cnt<1) return; 
IdxTy ns=0;
const IdxTy ltsz=cnt;
IdxTy lt[ltsz];
IdxTy l=lv[0]; ns+=l; lt[0]=l+1;  
for(IdxTy i=1; i<ltsz; ++i)  { l=lv[i]; ns+=l; lt[i]=l+1+lt[i-1]; } 

//for(IdxTy i=0; i<(cnt-1); ++i)  { l=lv[i]; ns+=l; lt[i+1]=l+lt[i]; } 
const IdxTy lvl=cnt-2;
RatCo tp;
make_nest(tp,lvl);
std::vector<int> pos,npos;
IdxTy ptr=0;
IdxTy tptr=1;
IdxTy i=0;
m=LocTy(ns+1);
const IdxTy ilim=ns>>1;
for(; i<=ilim; ++i)  
{ 
const IdxTy ie=i+1;
RatTy ni=mi.evaluate_polynomial(tp,RatTy(ie));
for(IdxTy k=0; k<pos.size(); ++k)  
 ni-=RatTy(npos[k])*mi.evaluate_polynomial(tp,RatTy(ie-pos[k]));
if ( print_res) { MM_ERR(MMPR2(i,ni)) } 
m[i]=(IdxTy(D(ni)));
m[ns-i]=m[i];
while (ptr<lv.size()) { if (lv[ptr]!=i) break; pos.push_back(ie); npos.push_back(1); ++ptr; }  
//while (tptr<ltsz) { { MM_ERR(MMPR4(i,tptr,ltsz,lt[tptr])) } if (lt[tptr]!=ie) break; pos.push_back(ie); npos.push_back(-1); ++tptr; }  
while (tptr<ltsz) { 
const IdxTy lb=lv[0]+lv[tptr]+2; 
{ MM_ERR(MMPR4(i,tptr,ltsz,lt[tptr])<<MMPR(lb)) }  

if (lb!=ie) break; pos.push_back(ie);
int neq=tptr+1;
while (neq<int(ltsz)) {if (lv[neq]!=lv[tptr]) break; ++neq; } 
int neqd=neq-tptr;
int neqdf=1;
for(int k=1; k<=neqd; ++k) neqdf*=k;
//++tptr;
tptr=neq;
 npos.push_back(-neqdf); 
 }  
//m.push_back(ni);
} // i  
} 




RatTy binom(const IdxTy top, const IdxTy bot) const 
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

RatCo sum_coefs(const IdxTy n, const IdxTy flags=0) const
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
void nesting(RatCo & d) const
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
} // sum_nest


void cmd_nest_table(Cip & cip , LocalVar & lv )
{
typedef mjm_integrals Mi;
Mi mi;
const IdxTy n=myatoi(cip.p1);
const IdxTy m=myatoi(cip.p2);
RatCo nest;
nest.push_back(1);
for(IdxTy i=1; i<=n; ++i)
{
nesting(nest);
RatCo x=nest;
// missing const term of zero 
x.push_back(0);
RatCo y;
mi.reverse_polynomial(y,x);
for(IdxTy j=0; j<m; ++j)
{
RatTy v= mi.evaluate_polynomial(y,RatTy(j));
MM_MSG(MMPR3(i,j,v))

} 
//std::cout<<ss.str()<<CRLF;
//std::cout<<format_poly(i,n,nest,style)<<CRLF;
} //i 
} // sum_nest





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
#error needs rationals 
#endif

std::map<IdxTy, IdxTy> find_multiplicities(LocTy & res, const LocTy & c, const IdxTy flags=0)
{
typedef mjm_ragged_card_itor Mrc;
Mrc mrc(c);
const bool print_res=((flags&1)==0);
//std::vector<IdxTy> c;
//c.push_back(
std::map<IdxTy, IdxTy> cnts;
for(mrc.reset(); mrc.ok(); ++mrc)
{
//MM_MSG(mti.to_string());
++cnts[mrc.sum()];
}
if ( print_res) {
int  dn=0;
int  last=0;
int  last2=0;
int  dn2=0;
MM_LOOP(ii,cnts) { 
const IdxTy i=(*ii).first;
const IdxTy m=(*ii).second;
dn=m-last;
dn2=m-last2;
last2=last;
last=m;
MM_MSG(MMPR4(i,m,dn,dn2)) }
MM_LOOP(ii,c) { MM_MSG("C "<<MMPR((*ii))) }
} // print_res
for(IdxTy i=0; i<cnts.size(); ++i) res.push_back(cnts[i]); 
return cnts;
}

D  count_spaces( const IdxTy L, const LocTy & c, const IdxTy flags=0)
{
Mi mi;
IdxTy sum=0;
MM_LOOP(ii,c) sum+=(*ii);
const IdxTy total=L-sum;
const IdxTy n=c.size()+1;
RatCo d; 
make_nest( d,  n-2 ); //  const
return ((n)*D(mi.evaluate_polynomial(d,RatTy(total+1))));
}
template <class Tx, class Ty> 
void counts1d(Tx& cache, Ty & p, Ty & P , const IdxTy L, const LocTy & c, IdxTy flags=0)
{
Mi mi;
IdxTy sum=0;
const IdxTy offset=2;
MM_LOOP(ii,c) sum+=(*ii);
const IdxTy total=1+L-sum;
const IdxTy n=c.size()+1;
// p_spaces01 calculates the higher order poly so it can get the 
// normalization with few evals but we are doing all of them 
// anyway so just use lower order
const IdxTy highest=c.size()+2;
//std::vector<RatCo>  cache(highest+1);
if (cache.size()<(highest+1)) cache= Tx(highest+1);
for(IdxTy i=0; i<cache.size(); ++i) 
	cache[i]=(simplify(fix_f(sum_coefs(i,1))));

const RatCo & d=cache[n-offset];
RatTy Pl=0;
for(int  i=0; i<=total; ++i)
{
 RatTy pi=(RatTy(n)*(mi.evaluate_polynomial(d,RatTy(total-i))));
Pl=Pl+pi;
pi.simplify();
Pl.simplify();
p.push_back(pi);
P.push_back(Pl);
}

}

void p_spaces01( D & p0, D & p1, const IdxTy L, const LocTy & c
, D * some, const IdxTy nsome, const IdxTy flags=0)
{
Mi mi;
IdxTy sum=0;
// changing to 1 small effect but sums a little worse 
const IdxTy offset=2;
MM_LOOP(ii,c) sum+=(*ii);
//const IdxTy total=L-sum;
MM_ONCE(" test kluge forgot plus one ",)
const IdxTy total=1+L-sum;
const IdxTy n=c.size()+1;
static std::vector<RatCo> cache;
while (cache.size()<=(n-offset))
{
//MM_ONCE(" creating static junk bin mem leak suspect",)
//MM_ERR(" making "<<cache.size())
RatCo d; 
make_nest( d,  cache.size() ); //  const
for(IdxTy i=0; i<d.size(); ++i) d[i].simplify();
cache.push_back(d);
}

//make_nest( d,  n-2 ); //  const
//const RatCo & d=cache[n-2];
const RatCo & d=cache[n-offset];
const int newo=0;
//const D t= ((n)*D(mi.evaluate_polynomial(d,RatTy(total+1))));
const D t= ((n)*D(mi.evaluate_polynomial(d,RatTy(total+newo))));
const D it=1.0/t;
const D t1= ((n)*D(mi.evaluate_polynomial(d,RatTy(total+newo-1))));
const D t2= ((n)*D(mi.evaluate_polynomial(d,RatTy(total+newo-2))));
// this is hte p of getting space length of exactly zero. 
p0=(t-t1)/t;
p1=(t1-t2)/t;

if (some!=0)
{
D prior=t;
for(IdxTy i=0; i<nsome; ++i)
{
// TODO verify preceisely what these mean.. 
//some[i]= ((n)*D(mi.evaluate_polynomial(d,RatTy(total-i))))/t;
const D p= ((n)*D(mi.evaluate_polynomial(d,RatTy(total-i-1))));
// integraed probability, space lt i 
some[i]=p*it;
// point probability 
some[i+nsome]=(prior-p)*it;
prior=p;
}

}

}


RatCo fix_f(const RatCo & x)
{
Mi mi;
RatCo y=x;
y.push_back(0);
RatCo d; 
mi.reverse_polynomial(d,y);
return d; 
}

std::map<IdxTy, IdxTy> make_spaces(LocTy & res, const IdxTy L, const LocTy & c, const IdxTy flags=0)
{
Mi mi;
IdxTy sum=0;
MM_LOOP(ii,c) sum+=(*ii);
const IdxTy total=L-sum;
const IdxTy n=c.size()+1;
res=LocTy(total+1);
RatCo x= sum_coefs(n-3);
//void make_nest(RatCo & d, const IdxTy n ) const
// missing const term of zero 
x.push_back(0);
RatCo d; 
mi.reverse_polynomial(d,x);
make_nest( d,  n-3 ); //  const
for(IdxTy i=0; i<=total; ++i)
{
res[i]=(n)*D(mi.evaluate_polynomial(d,RatTy(total+1-i)));

}
std::map<IdxTy, IdxTy> cnts;

return cnts;

}


std::map<IdxTy, IdxTy> find_spaces(LocTy & res, const IdxTy L, const LocTy & c, const IdxTy flags=0)
{

//mjm_total_spaces_itor(const LocTy & l, const IdxTy L )
//:  m_l(l), m_L(L), m_ls(Sum(m_l)), m_ss(m_L-m_ls)

typedef mjm_total_spaces_itor Mrc;
Mrc mrc(c,L);
const bool print_res=((flags&1)==0);
//std::vector<IdxTy> c;
//c.push_back(
std::map<IdxTy, IdxTy> cnts;
for(mrc.reset(); mrc.ok(); ++mrc)
{
//MM_MSG(mti.to_string());
MM_SZ_LOOP(i,mrc.spaces(),msz) { ++cnts[mrc.spaces()[i]]; } 

}
if ( print_res) {
int  dn=0;
int  last=0;
int  last2=0;
int  dn2=0;
MM_LOOP(ii,cnts) { 
const IdxTy i=(*ii).first;
const IdxTy m=(*ii).second;
dn=m-last;
dn2=m-last2;
last2=last;
last=m;
MM_MSG(MMPR4(i,m,dn,dn2)) }
MM_LOOP(ii,c) { MM_MSG("C "<<MMPR((*ii))) }
} // print_res
for(IdxTy i=0; i<cnts.size(); ++i) res.push_back(cnts[i]); 
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
LocTy brute;
find_multiplicities(brute,c);
}

//void make_multiplicities_with_lens(LocTy &  m, const LocTy & lv) const
void cmd_model_mult(Cip & cip , LocalVar & lv )
{
LocTy c;
//for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
//const IdxTy n=myatoi(cip.p1);
IdxTy x=0;
IdxTy i=2;
while ((x=myatoi(cip.wif(i)))!=0) { c.push_back(x); ++i; } 
LocTy res;
MM_ERR(" model_mult "<<MMPR(c.size()))
make_multiplicities_with_lens( res, c); 
//find_multiplicities(c);
}
void cmd_hello(Cip & cip , LocalVar & lv )
{
typedef long  int Ti;
 Ti x=(Ti('h')<<(8*0));
 x|=(Ti('h')<<(8*1));
 x|=(Ti('e')<<(8*2));
 x|=(Ti('l')<<(8*3));
 x|=(Ti('l')<<(8*4));
 x|=(Ti('o')<<(8*5));
 x|=(Ti('x')<<(8*7));
x=8647033796017088616L;
// x|=('l'<<(8*4));
// x|=('l'<<(8*3));
std::string s((char*)&x);
MM_MSG(MMPR2(x,s))
}
RatCo&  simplify( RatCo && x)
{
for (IdxTy i=0; i<x.size(); ++i) x[i].simplify(); 
return x;
}
void cmd_value2dn_sum(Cip & cip , LocalVar & lv )
{
const IdxTy L=myatoi(cip.p1);
LocTy c;
IdxTy i=3;
IdxTy x=0;
while ((x=myatoi(cip.wif(i)))!=0) {MM_ERR(" length "<<x)  c.push_back(x); ++i; } 
MM_ERR(" cmd_value2dn_sum "<<MMPR(c.size()))
std::vector<RatCo> cache;
std::vector<RatTy>  p,P;
counts1d( cache, p,  P ,  L, c, 0);
for(IdxTy i=0; i<10; ++i)
{
MM_ERR(MMPR3(i,p[i],P[i]))
}
// P(o)+p(o)*(l-o-1)
const IdxTy l1=c[0];
const IdxTy l2=c[1];
for(IdxTy o1=1; o1<l1; ++o1)
{
for(IdxTy o2=1; o2<l2; ++o2)
{
const IdxTy o=o1+o2;
//const RatTy po1=P[o1]+(l1-o1-1)*p[o1];
//const RatTy po2=P[o2]+(l2-o2-1)*p[o2];
//MM_ERR(MMPR4(o,o1,o2,po1)<<MMPR4(po2,(po1*po2),P[o1],p[o1]))
}
}


}

void cmd_value2d_sum(Cip & cip , LocalVar & lv ){

Mi mi;
const IdxTy n=myatoi(cip.p1);
const IdxTy k=myatoi(cip.p2);
MM_ERR(MMPR2(n,k))
RatCo res(2+4*n);
const IdxTy highest=2*n+2;
std::vector<RatCo>  cache(highest+1);
for(IdxTy i=0; i<cache.size(); ++i) 
	cache[i]=(simplify(fix_f(sum_coefs(i,1))));
MM_ERR("done making cache")
const RatCo& an=cache[n-1];
// note this should be the nested G functions but not doing cross term ye 
const RatCo& bn=cache[n-1];

for(IdxTy j=1; j<=n; ++j)  { 
const RatTy & b=bn[j];
MM_ERR(MMPR2(j,b))
for(IdxTy k=0; k<=j; ++k)  { 
const int c=((k&1)==0)?1:(-1);
const RatTy bin=binom(j,k)*RatTy(c);
for(IdxTy i=1; i<=n; ++i)  {
const RatTy & a=an[i];
const RatTy abb= a*b*bin;
const IdxTy pm=i+k-1;
if (pm>=cache.size()) { MM_ERR(" oor "<<MMPR3(i,k,cache.size())) } 
const RatCo & f= cache[pm];
for(IdxTy p=1; p<=pm; ++p)  {
res[j-k+p]+=abb*bin*f[p];
} // p 
res[j+i]-=abb*bin;
} // i 
} // k
} // j

MM_MSG(mi.print_polynomial(res));
RatTy total=0;
for(IdxTy j=0; j<res.size(); ++j) total+=res[j];
total.simplify();
MM_MSG(MMPR(total))
for(IdxTy j=0; j<res.size(); ++j)
{
RatTy v=mi.evaluate_polynomial(res,RatTy(j));
v.simplify();
MM_MSG(MMPR3(j,v,D(v)))
} 


} // value2d

void cmd_fvalue_sum(Cip & cip , LocalVar & lv )
{
Mi mi;
const IdxTy n=myatoi(cip.p1);
const IdxTy k=myatoi(cip.p2);
MM_ERR(MMPR2(n,k))
RatCo res;
RatCo x= sum_coefs(n,1);
RatCo y=fix_f(x);
// now y is in normal order, index is power of variable
const IdxTy sz=y.size();
for(IdxTy i=0; i<sz; ++i)
{

RatCo f= fix_f(sum_coefs(i+k,1));
RatCo term;
mi.multiply_polynomials(term,f,y);
mi.accumulate_polynomial(res, term);
for(IdxTy j=0; j<res.size(); ++j) res[j].simplify();

} // i 

MM_MSG(mi.print_polynomial(res));
RatTy total=0;
for(IdxTy j=0; j<res.size(); ++j) total+=res[j];
total.simplify();
MM_MSG(MMPR(total))
for(IdxTy j=0; j<res.size(); ++j)
{
RatTy v=mi.evaluate_polynomial(res,RatTy(j));
v.simplify();
MM_MSG(MMPR3(j,v,D(v)))
} 

}

void cmd_both_mult(Cip & cip , LocalVar & lv )
{
LocTy c;
//for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
//const IdxTy n=myatoi(cip.p1);
IdxTy x=0;
IdxTy i=2;
while ((x=myatoi(cip.wif(i)))!=0) { c.push_back(x); ++i; } 
LocTy res,rest,brute;
MM_ERR(" both_mult "<<MMPR(c.size()))
make_multiplicities_with_lens( res, c,1); 
MM_ERR(" both_mult "<<MMPR(res.size()))
make_multiplicities_table_lens(rest,c,0,1); 
MM_ERR(" both_mult "<<MMPR(rest.size()))
//make_multiplicities_with_less_old_lens(res,c,1); 
//find_multiplicities(brute, c,1);
MM_ERR(" both_mult "<<MMPR(brute.size()))
//if (true) return;
IdxTy tmod=0;
IdxTy tmodt=0;
IdxTy tloop=0;
for(IdxTy i=0; i<res.size(); ++i) { 
const int mod=res[i];
const int modt=res[i];
const int loop=0 ; // brute[i];
//const int d=mod-loop;
const int d=mod-modt;
tmod+=mod;
tmodt+=modt;
tloop+=loop;
MM_MSG(MMPR4(i,mod,modt,loop)<<MMPR4(d,tmod,tmodt,tloop))}
//find_multiplicities(c);
}


void cmd_both_spaces(Cip & cip , LocalVar & lv )
{
LocTy c;
//for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
const IdxTy L=myatoi(cip.p1);
IdxTy x=0;
IdxTy i=2;
while ((x=myatoi(cip.wif(i)))!=0) { c.push_back(x); ++i; } 
LocTy res,rest,brute;
MM_ERR(" both_spaces "<<MMPR(c.size()))
make_spaces( res, L,c,1); 
MM_ERR(" both_spaces "<<MMPR(res.size()))
//make_multiplicities_table_lens(rest,c,0,1); 
MM_ERR(" both_spaces "<<MMPR(rest.size()))
//make_multiplicities_with_less_old_lens(res,c,1); 
find_spaces(brute, L,c,1);
MM_ERR(" both_spaces "<<MMPR(brute.size()))
//if (true) return;
IdxTy tmod=0;
//IdxTy tmodt=0;
IdxTy tloop=0;
for(IdxTy i=0; i<brute.size(); ++i) { 
const int mod=res[i];
//const int modt=res[i];
const int loop= brute[i];
//const int d=mod-loop;
//const int d=mod-modt;
tmod+=mod;
//tmodt+=modt;
tloop+=loop;
//MM_MSG(MMPR4(i,mod,modt,loop)<<MMPR4(d,tmod,tmodt,tloop))}
MM_MSG(MMPR3(i,loop,mod)<<MMPR2(tloop,tmod))}
//find_multiplicities(c);
const IdxTy cs=count_spaces(L,c,1);
D p0,p1;
p_spaces01( p0, p1, L,  c,0,0);
MM_MSG(MMPR3(cs,p0,p1))

}





void cmd_sum_code(Cip & cip , LocalVar & lv )
{
const IdxTy n=myatoi(cip.p1);
sum_coefs(n);
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
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
//void cmd_find_mult(Cip & cip , LocalVar & lv )
m_cmd_map[StrTy("find-mult")]=&Myt::cmd_find_mult;
m_cmd_map[StrTy("model-mult")]=&Myt::cmd_model_mult;
m_cmd_map[StrTy("both-mult")]=&Myt::cmd_both_mult;
m_cmd_map[StrTy("both-spaces")]=&Myt::cmd_both_spaces;
//void cmd_sum_code(Cip & cip , LocalVar & lv )
m_cmd_map[StrTy("sum-code")]=&Myt::cmd_sum_code;
#ifdef USE_RATIONALS
m_cmd_map[StrTy("coef-table")]=&Myt::cmd_coef_table;
m_cmd_map[StrTy("sum-nest")]=&Myt::cmd_sum_nest;
m_cmd_map[StrTy("fvalue-sum")]=&Myt::cmd_fvalue_sum;
m_cmd_map[StrTy("value2d-sum")]=&Myt::cmd_value2d_sum;
m_cmd_map[StrTy("value2dn-sum")]=&Myt::cmd_value2dn_sum;
m_cmd_map[StrTy("hello")]=&Myt::cmd_hello;
m_cmd_map[StrTy("nest-table")]=&Myt::cmd_nest_table;
#endif

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
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
RaggedMap m_ragged_map;
MyUBlock m_soi_table;

CounterMap m_cm;
CliTy m_cli;
}; //mjm_trees_and_tables



/////////////////////////////////////////////////////////

#ifdef  TEST_RATIONAL_THINGS__
int main(int argc,char **args)
{
typedef mjm_rational_things  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

