
#ifndef MJM_mjsonu_H__
#define MJM_mjsonu_H__
 
#include "mjm_globals.h"

#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_string_array_formatter.h"

// TODO FIXME this or similar needs fftw3 to compile wtf3? lol
#include "mjm_canned_methods.h"


#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"

//#include "rapidjson/reader.h"
//#include "rapidjson/istreamwrapper.h"
#include "rapidjsonalt/rapidjson/reader.h"
#include "rapidjsonalt/rapidjson/istreamwrapper.h"

#include <iostream>
#include <sstream>





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

2019-02-16 copy from violin plot to make template 

g++ -std=gnu++11 -DTEST_mjsonu__ -gdwarf-3 -O0 -I. -I/home/documents/cpp/mjm/hlib/ -I/home/documents/cpp/mjm/num -I/home/documents/cpp/mjm/ -I/home/documents/cpp/mjm/tcl -I/home/documents/cpp/mjm/include -I/home/documents/cpp/pkg/include -I/home/documents/cpp/pkg -Wall -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -Wno-misleading-indentation -x c++ ./mjsonu.h -o mjsonu.out -lreadline




*/


////////////////////////////////////////////////////////////////

class mjsonu_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
mjsonu_params( const StrTy & nm) : Super(nm) {}
mjsonu_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
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
ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;
return ss.str();
}


}; // mjsonu_params


namespace mjsonu_traits
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
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // mjsonu_traits
///////////////////////////////////////////////////////////////
// need to eliminate crap 
using namespace rapidjson;
//using namespace std;

// copied rapidjson simplepullreader
struct MyHandler {

	typedef unsigned int IdxTy; 
	typedef std::stringstream Ss; 
	typedef std::string StrTy; 
typedef mjsonu_traits::Tr Tr; 
typedef mjm_string_array_formatter<Tr> Formatter;

// If you can require C++11, you could use std::to_string here
template <typename T> std::string stringify(T x) {
// TODO FIXME this needs more preceisions for floats to this is a mess. 
    	std::stringstream ss;
    ss << x;
    return ss.str();
}

class object_state
{
	public:
object_state(): m_kind(0),m_index(0),m_key() {}
object_state(const StrTy & k, const IdxTy n ): m_kind(n),m_index(0),m_key(k) {}
object_state(const StrTy & k, const IdxTy n, const IdxTy idx )
: m_kind(n),m_index(idx),m_key(k) {}
		IdxTy m_kind;
		IdxTy m_index;
	StrTy m_key;
}; // object_state

typedef mjm_ragged_table Ragged;

    const char* type;
    std::string data;
std::vector<object_state> m_stack;
StrTy m_key,m_filter;
bool m_dump_all;
IdxTy m_out_mode;
Formatter m_formatter;
Ragged * m_table;
bool m_no_stdout;
IdxTy m_global; // 2022-01-07
void set_table(Ragged * p ) { m_table=p;}
void no_stdout(const bool x ) { m_no_stdout=x;}

void new_parse() { Init(); } 
void Init()
{
if ( m_stack.size())
{
//IdxTy idx=m_stack.back().m_index;
//++idx;
m_stack.clear();
/// 2022-07-01 july move inc first
++m_global;

	   m_stack.push_back(object_state(StrTy("GLOBAL"),0,m_global)); 
//++m_global;
}
else{
// clear is non functiona; 
m_stack.clear();
// 2022-07-03 now it will not start at zero but everything should work.. 
++m_global;
	   m_stack.push_back(object_state(StrTy("GLOBAL"),0,m_global)); 
	  // m_stack.push_back(object_state(StrTy("GLOBAL"),0)); 
}

type=0;
m_table=NULL;
m_no_stdout=false;
} // Init
   //MyHandler(): m_filter(),m_dump_all(true)  {} 

    MyHandler() : type(), data() ,m_filter(),m_dump_all(true),m_out_mode(0) ,
m_global(0)  
{
Init(); //m_stack.push_back(object_state(StrTy("GLOBAL"),0)); 
    }
    MyHandler(const IdxTy ng) : type(), data() ,m_filter(),m_dump_all(true),m_out_mode(0) ,
m_global(ng)  
{
Init(); //m_stack.push_back(object_state(StrTy("GLOBAL"),0)); 
    }


IdxTy global() const { return m_global; }

void clean_output_mode(const bool x ) {if (x) m_out_mode=1; else m_out_mode=0; }  
void out_mode(const IdxTy m) { m_out_mode=m;} 
bool stack_ok()
{
const bool x=  ( m_stack.size()!=0);
if ( !x) MM_ERR(" json stack underflow ")

return x; 
}

bool new_item() 
{
switch (m_out_mode)
{
case 1: return new_item_2();
case 0:
default : return new_item_1();
} // switch

} // new_item
// historical brackets and semicolons etc
bool new_item_1() 
{
if (stack_ok())
{
object_state & x = m_stack.back();
++x.m_index;
const StrTy sep=";";
StrTy y="";
Ss ss;
//MM_LOOP(ii,m_stack) { y=y+sep+(*ii).m_key; }
MM_LOOP(ii,m_stack) { 
const StrTy & k=(*ii).m_key;
//if (k!="") 
{ 	ss<<sep+k<<"["<<(*ii).m_index<<"]";  } 

}
y=ss.str();
if (m_key=="")
{
MM_SZ_LOOP(i,m_stack,sz)
{
auto & x=m_stack[sz-1-i];
if (x.m_key!="") { m_key=x.m_key; break; }
}

}
if (m_dump_all) { MM_MSG(y<<":"<<m_key<<"="<<data) } 
else 
{
if (m_filter==m_key)
{
MM_MSG(MMPR(m_key))
	std::cout<<data<<CRLF;

}

}
m_key="";
}
	return true; 
} // new_item_1



// better format for parsing. 
bool new_item_2() 
{
typedef std::vector<StrTy> Line;
Line l;
if (stack_ok())
{
object_state & x = m_stack.back();
// 2022-01-7
//++x.m_index;
const StrTy sep=";";
StrTy y="";
//Ss ss;
//MM_LOOP(ii,m_stack) { y=y+sep+(*ii).m_key; }
MM_LOOP(ii,m_stack) { 
const StrTy & k=(*ii).m_key;
//if (k!="") 
{ 
// this looks backwards... 
l.push_back(k);
Ss ss;
ss<<(*ii).m_index;
l.push_back(ss.str());
//	ss<<sep+k<<"["<<(*ii).m_index<<"]";  

} 


}
//y=ss.str();
if (m_key=="")
{
MM_SZ_LOOP(i,m_stack,sz)
{
auto & x=m_stack[sz-1-i];
if (x.m_key!="") { m_key=x.m_key; break; }
}
}
l.push_back(m_key);
l.push_back("=");
l.push_back(data);
//{ MM_MSG(y<<":"<<m_key<<" = "<<data) } 
if (!m_no_stdout ) MM_MSG(m_formatter.format(l)) 
if (m_table) (*m_table).add(l);
#if 0 
if (m_dump_all) { MM_MSG(y<<":"<<m_key<<"="<<data) } 
else 
{
if (m_filter==m_key)
{
MM_MSG(MMPR(m_key))
	std::cout<<data<<CRLF;

}

}
#endif

m_key="";
}
	return true; 
} // new_item_2




void pop_object()
{
if ( stack_ok())
{
object_state & x = m_stack.back();
if ( x.m_kind!=1)
{
MM_ERR(" unkind "<<MMPR2(x.m_kind,1))
}
m_stack.pop_back();
}
else MM_ERR(" json stack underflow ")
}

void pop_array()
{
if ( m_stack.size()!=0)
{
object_state & x = m_stack.back();
if ( x.m_kind!=2)
{
MM_ERR(" unkind "<<MMPR2(x.m_kind,2))
}
m_stack.pop_back();
}
else MM_ERR(" json stack underflow ")

}
bool new_thing(const IdxTy n)
{
if (stack_ok())
{
object_state & x = m_stack.back();
// 2022-01-07
//++x.m_index;
}
// 2022-07-01 wtf 2022-01-07? now it IS July lol... 
//m_stack.push_back(object_state(m_key,n)); 
// this is picking up a stale object numer move increment first doh
++m_global;
m_stack.push_back(object_state(m_key,n,m_global)); 
//++m_global;
m_key="";
return true;
}
//typedef rapidjson::Null Jnull;
//typedef rapidjson::Boll Jbool;
//typedef rapidjson::Int  Jint;
//typedef rapidjson::Int64  Jint64;
//typedef rapidjson::Uint  Juint;
//typedef rapidjson::Uint64  Juint64;
//typedef rapidjson::Double  Jdouble;
typedef rapidjson::SizeType  Jsize;

    bool Null() { type = "Null"; data.clear();  return new_item(); }
    bool Bool(bool b) { type = "Bool:"; data = b? "true": "false"; return new_item(); }
    bool Int(int i) { type = "Int:"; data = stringify(i); return new_item(); }
    bool Uint(unsigned u) { type = "Uint:"; data = stringify(u); return new_item(); }
    bool Int64(int64_t i) { type = "Int64:"; data = stringify(i); return new_item(); }
    bool Uint64(uint64_t u) { type = "Uint64:"; data = stringify(u); return new_item(); }
    	// arrrgrgh wtf??? precisions gones .... doh 
    bool Double(double d) { type = "Double:"; data = stringify(d); return new_item(); }
    bool RawNumber(const char* str, Jsize length, bool) { type = "Number:"; data = std::string(str, length); return new_item(); }
    bool String(const char* str, SizeType length, bool) { type = "String:"; data = std::string(str, length); return new_item(); }
    bool StartObject() { type = "StartObject"; data.clear(); 
	return new_thing(1); } 
 //	    m_stack.push_back(object_state(m_key,1)); 
	    
	   // return true; }
    bool Key(const char* str, SizeType length, bool) { type = "Key:"; data = std::string(str, length); m_key=data;  return true; }
    bool EndObject(SizeType memberCount) { type = "EndObject:"; data = stringify(memberCount); 
	   pop_object(); 
	    return true; }
    bool StartArray() { type = "StartArray"; data.clear();
	   
	return new_thing(2); } 
//	   m_stack.push_back(object_state(m_key,2)); 
//	    return true; }
    bool EndArray(SizeType elementCount) { type = "EndArray:"; data = stringify(elementCount); pop_array();  return true; }

   void select(const StrTy & f) { m_filter=f; m_dump_all=(f==""); } 
    private:
    MyHandler(const MyHandler& noCopyConstruction);
    MyHandler& operator=(const MyHandler& noAssignment);
};
//////////////////////////////////////////////////////////////////////////







class mjm_mjsonu 
{
typedef mjsonu_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_mjsonu Myt;
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
typedef mjsonu_params ParamGlob;
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
mjm_mjsonu():m_dmel(new Dmel()) {Init();}
mjm_mjsonu(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_mjsonu()
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
if (s=="-quit") { MM_ERR(" quitting  "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
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

typedef std::map<StrTy, std::vector<StrTy> > VecMap;
typedef std::map<StrTy, VecMap > TreeMap;

StrTy any_kv(const Ragged & r, const StrTy & k, const IdxTy n, const IdxTy flags)
{
IdxTy hits=0;
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l =r[i];
const IdxTy len=l.size();
if (len<3) continue;
if (l[len-3]==k) 
{
if (hits==n) return l[len-1];
++hits;
}
} // i 


return StrTy();
}// any_kv

TreeMap index( const Ragged & r, const IdxTy flags)
{
TreeMap tm;
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l =r[i];
const IdxTy len=l.size();
if (len<3) continue;
StrTy key="";
for(IdxTy j=0; j<len-3; ++j)
{
key=key+l[j];

} // j 
tm[key][l[len-3]].push_back(l[len-1]);
} // i 
return  tm;
} // index

StrTy kvm( const TreeMap & tm, const StrTy & k1, const StrTy & v1, const StrTy & k2)
{
MM_LOOP(ii,tm)
{
auto &m =(*ii).second;
auto jj=m.find(k1);
if (jj==m.end()) continue;
auto v=(*jj).second;
if ( v[0]==v1)
{
jj=m.find(k2);
if (jj==m.end()) continue;
return (*jj).second[0];
} // found 

}
return StrTy();
} // kvm
 


//typedef mjm_ragged_table Ragged;
class parse_result_
{
public:
parse_result_() {Init();}
void Init() { count=0; } 
StrTy text;
Ragged res;
IdxTy count;
}; // parse_result_

typedef parse_result_ parse_result;
//bool Bit(const IdxTy f, const IdxTy b) const { return (f&(1<<b))!=0; } 
parse_result dirty_parse(IsTy & is, const IdxTy flags)
{
parse_result pr;
const bool clean_output_mode=true;
const bool dump_here=!true;
const bool reset_handler=true;
const bool debug_out=!true;
const bool do_nostd =true;
const bool good=is.good();
const bool eof=is.eof();
const StrTy cmd1="";
//MM_ERR(MMPR4(cip.cmd(),good,eof,use_stdin)<<MMPR4(select,flags,cmd1,cmd2))
MM_ERR(MMPR4(__FUNCTION__,flags,clean_output_mode,dump_here)<<MMPR4(reset_handler, debug_out, do_nostd,good)<<MMPR(eof))
//  json reader and  zero length file TODO FIXME?  
if (eof||!good)
{
MM_ERR(" null input to json parser "<<MMPR2(eof,good))
return pr; 
}
const bool abort=false;
    MyHandler handler;
    handler.select(cmd1);
	if (clean_output_mode) handler.clean_output_mode(true);
if (do_nostd) handler.no_stdout(true);
handler.set_table(&pr.res);
    Reader reader;
    //StringStream ss(json);
    BasicIStreamWrapper<std::istream>  ss(is);
    reader.IterativeParseInit();
    while (!reader.IterativeParseComplete()) {
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good;
if (bad)
{
MM_ERR("  json dies before parse "<<MMPR2(eof,good))
pr.count=handler.global();
return pr;
}
 } // ParseComplete  
// this logic is now wrong, 
const bool mjm_plod=true;
//bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler)||mjm_plod;
bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler);
if (mjm_plod&&!ok) { 
if ( debug_out) { MM_ERR(" reseting ok "<<MMPR4(ok,is.good(),is.eof(),ss.Tell())<<MMPR(int(ss.Peek()))) }
    ss.Take();  reader.IterativeParseInit(); ok=true; 
// needs the  again... 
if (reset_handler){
 handler.new_parse();
if( do_nostd) handler.no_stdout(true);
handler.set_table(&pr.res);

}
// I thought this had to reset the stack in the hanlder 
} 
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good||!ok;
if (bad)
{
const IdxTy n=ss.Tell();
IdxTy nmax=(n>20)?20:n;
char cc[nmax+2];
for(IdxTy i=0; i<nmax; ++i) cc[i]=ss.Take();
cc[nmax]=0;
MM_ERR("  json dies before parse  "<<MMPR4(ok,eof,good,n)<<MMPR(cc))
pr.count=handler.global();

return pr;
}
 } 
if ( dump_here) {     std::cout << handler.type << handler.data << std::endl; } 
    } // IterativeParseComplete
pr.count=handler.global();
return pr;
} // dirty_parse;

//////////////////////////////////////////////////////////////////////

parse_result better_parse(IsTy & is, const IdxTy firstn, const IdxTy flags)
{
parse_result pr;
pr.count=firstn;
const bool clean_output_mode=!Bit(flags,0); // true;
// FIXME Danger Will Robinson setting this will cause abort??? 
const bool dump_here=Bit(flags,1); // !true;
if (dump_here)
{ MM_ERR(" setting dump_here causes abort wtf")
}
const bool reset_handler=!Bit(flags,2); // true;
const bool debug_out=Bit(flags,3); // !true;
const bool do_nostd =!Bit(flags,4); // true;
const bool good=is.good();
const bool eof=is.eof();
const StrTy cmd1="";
//MM_ERR(MMPR4(cip.cmd(),good,eof,use_stdin)<<MMPR4(select,flags,cmd1,cmd2))
MM_ERR(MMPR4(__FUNCTION__,flags,clean_output_mode,dump_here)<<MMPR4(reset_handler, debug_out, do_nostd,good)<<MMPR2(eof,firstn))
//  json reader for zero length file


if (eof||!good)
{
MM_ERR(" null input to json parser "<<MMPR2(eof,good))
return pr; 
}
const bool abort=false;
    MyHandler handler(firstn);
    handler.select(cmd1);
	if (clean_output_mode) handler.clean_output_mode(true);
if (do_nostd) handler.no_stdout(true);
handler.set_table(&pr.res);
    Reader reader;
    //StringStream ss(json);
    BasicIStreamWrapper<std::istream>  ss(is);
    reader.IterativeParseInit();
    while (!reader.IterativeParseComplete()) {
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good;
if (bad)
{
MM_ERR("  json dies before parse "<<MMPR2(eof,good))
pr.count=handler.global();
return pr;
}
 } // ParseComplete  
// this logic is now wrong, 
const bool mjm_plod=true;
//bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler)||mjm_plod;
bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler);
if (mjm_plod&&!ok) { 
if ( debug_out) { MM_ERR(" reseting ok "<<MMPR4(ok,is.good(),is.eof(),ss.Tell())<<MMPR(int(ss.Peek()))) }
    ss.Take();  reader.IterativeParseInit(); ok=true; 
// needs the  again... 
if (reset_handler){
 handler.new_parse();
if( do_nostd) handler.no_stdout(true);
handler.set_table(&pr.res);

}
// I thought this had to reset the stack in the hanlder 
} 
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good||!ok;
if (bad)
{
const IdxTy n=ss.Tell();
IdxTy nmax=(n>20)?20:n;
char cc[nmax+2];
for(IdxTy i=0; i<nmax; ++i) cc[i]=ss.Take();
cc[nmax]=0;
MM_ERR("  json dies before parse  "<<MMPR4(ok,eof,good,n)<<MMPR(cc))
pr.count=handler.global();

return pr;
}
 } 
if ( dump_here) {     std::cout << handler.type << handler.data << std::endl; } 
    } // IterativeParseComplete
pr.count=handler.global();
return pr;
} // better_parse;










/* for page scraping, json is as dirty as html for which sed/grep
work where the real parser does not on at least one occasion.
Its important to either get everything or just the proper stuff
and flip modes

*/
void cmd_sax_json(Cip & cip , LocalVar & lv )
{
// better than typing but aaarrrfj
using namespace rapidjson;
using namespace std;

const StrTy ragin=cip.p1;
// this is both a flag and mask string? WTF
//const StrTy select=cip.p2;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const StrTy select=cmd2;
const bool use_stdin=(ragin=="-")||(ragin.length()==0);
const bool dump_here=false; 
const bool clean_output_mode=Bit(flags,0);
// this changes default behavior
const bool reset_handler=!Bit(flags,1);
std::ifstream rin(ragin.c_str());
std::istream&  is=use_stdin?(std::cin):(rin);
const bool good=is.good();
const bool eof=is.eof();
MM_ERR(MMPR4(cip.cmd(),good,eof,use_stdin)<<MMPR4(select,flags,cmd1,cmd2))
// or zero length file 
if (eof||!good)
{
MM_ERR(" null input to json parser "<<MMPR2(eof,good))
return; 
}
const bool abort=false;
    MyHandler handler;
    //handler.select(select);
    handler.select(cmd1);
	if (clean_output_mode) handler.clean_output_mode(true);
    Reader reader;
    //StringStream ss(json);
    BasicIStreamWrapper<std::istream>  ss(is);
    reader.IterativeParseInit();
    while (!reader.IterativeParseComplete()) {
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good;
if (bad)
{
MM_ERR("  json dies before parse "<<MMPR2(eof,good))
return;
}
 } 
// this logic is now wrong, 
const bool mjm_plod=true;
//bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler)||mjm_plod;
bool ok= reader.IterativeParseNext<kParseDefaultFlags>(ss, handler);
if (mjm_plod&&!ok) { 
MM_ERR(" reseting ok "<<MMPR4(ok,is.good(),is.eof(),ss.Tell())<<MMPR(int(ss.Peek())))
    ss.Take();  reader.IterativeParseInit(); ok=true; 
if (reset_handler) handler.new_parse();
// I thought this had to reset the stack in the hanlder 
} 
{const bool good=is.good();
const bool eof=is.eof();
const bool bad=eof||!good||!ok;
if (bad)
{
const IdxTy n=ss.Tell();
IdxTy nmax=(n>20)?20:n;
char cc[nmax+2];
for(IdxTy i=0; i<nmax; ++i) cc[i]=ss.Take();
cc[nmax]=0;
MM_ERR("  json dies before parse  "<<MMPR4(ok,eof,good,n)<<MMPR(cc))

return;
}
 } 
if ( dump_here) {     cout << handler.type << handler.data << std::endl; } 
    } // IterativeParseComplete
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
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;

//m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
//m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
//m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
//m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("sax-json")]=&Myt::cmd_sax_json;
//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
//m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
//m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
//m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
//m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
//m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
//m_cmd_map[StrTy("sax-json")]=&Myt::cmd_sax_json;
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
ss<<" mjm_mjsonu "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com Mon Feb 18 15:01:06 UTC 2019 "<<CRLF;
ss<<" A simple utility program to convert json into various line oriented formats of varying precisions.    "<<CRLF;
ss<<" It is mostly a front end for rapidjson which does the parsing "<<CRLF;
ss<<" RapidJSON  Copyright (C) 2015 THL A29 Limited, a Tencent company, and Milo Yip. ";
ss<< "RAPIDJSON_VERSION_STRING="; 
ss<< QUOTE(RAPIDJSON_VERSION_STRING); 
ss<<CRLF;
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
bool Bit(const IdxTy f, const IdxTy b) const { return (f&(1<<b))!=0; } 
void Init()
{
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

}; //mjm_mjsonu 

/////////////////////////////////////////////////////////

#ifdef  TEST_mjsonu__
int main(int argc,char **args)
{
typedef mjm_mjsonu Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif
