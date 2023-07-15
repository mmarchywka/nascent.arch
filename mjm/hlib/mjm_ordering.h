#ifndef MJM_ORDERING_H__
#define MJM_ORDERING_H__
 
#include "mjm_globals.h"
// for the presence absence vector 
//#include "mjm_char_mat.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
//#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
//#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
//#include "mjm_strings.h"
//#include "mjm_string_index.h"

//#include "mjm_cli_ui.h"
//#include "../mjm_fasta_ii.h"

// ragged table is here 
#include "mjm_collections.h"
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



////////////////////////////////////////////////////////////////

class ordering_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
ordering_params( const StrTy & nm) : Super(nm) {}
ordering_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
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

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace ordering_traits
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
//typedef mjm_block_matrix<D> MyBlock;
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
}; // trees_and_tables_traits

// derived from string_index class mjm_1_to_v
template <class Tscalar, class Tval, class Tvec >
class mjm_1_to_some
{
typedef mjm_1_to_some Myt;
protected:
typedef ordering_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
//typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;
//typedef std::map<char, IdxTy > Wcmap;
public:
typedef Tscalar key_t;
typedef Tval value_t;
//typedef std::vector<value_t> vector_t;
typedef Tvec vector_t;
typedef std::map<key_t,vector_t> object_t;
typename object_t::iterator begin() { return m_obj.begin(); }
typename object_t::iterator end() { return m_obj.end(); }
typename object_t::const_iterator begin() const { return m_obj.begin(); }
typename object_t::const_iterator end() const { return m_obj.end(); }

IdxTy size() const { return m_obj.size(); } 
void add(const key_t & k, const value_t & v) { m_obj[k].push_back(v); }
const vector_t & operator[]( const key_t& k) const
{
auto ii=m_obj.find(k); if (ii!=m_obj.end()) { return (*ii).second; }
return m_null;
}
const vector_t starting_with( const key_t& k) const
{
auto ii=m_obj.lower_bound(k); if (ii!=m_obj.end()) {
vector_t v;
const IdxTy len=k.length();
// TODO does not remove uniques etc 
while (strncmp((*ii).first.c_str(),k.c_str(),len)==0)
{
 append(v,(*ii).second);
++ii;
}

return v; 
 }
return m_null;
}

void append( vector_t & d, const vector_t & s)  const
{ MM_LOOP(ii,s) { d.push_back((*ii)); } 
}

void clear() { m_obj.clear(); }
private:
object_t m_obj;
vector_t m_null;
};  // mjm_1_to_some



class mjm_string_ordering 
{
typedef  ordering_traits::Tr  Tr;
typedef mjm_string_ordering Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;


typedef ordering_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef IdxTy OrdTy;
typedef StrTy TgtTy;

typedef std::vector<TgtTy> FwdTy;
typedef std::vector<OrdTy> OrderVec;
typedef FwdTy::const_iterator CItor;

typedef mjm_1_to_some<TgtTy,OrdTy, OrderVec> RevMap;

public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
typedef OrdTy order_t;
typedef OrderVec order_vector_t;
typedef TgtTy target_t;
typedef CItor const_iterator;

mjm_string_ordering():m_dmel(new Dmel()) {Init();}
~mjm_string_ordering() { delete m_dmel; }

template <class Ty> void add(const Ty & v) { MM_LOOP(ii,v) { Next((*ii)); } } // add
// this fails when in and out are the same types... 
const TgtTy & operator[](const OrdTy i ) const  { return m_fwd[i]; } 

const order_vector_t  & operator[](const TgtTy &s) const  { return m_rev[s]; } 
const order_vector_t   starting_with(const TgtTy &s) const  { return m_rev.starting_with(s); } 

StrTy dump(const IdxTy flags ) const  { return to_string(flags); } 
StrTy to_string(const IdxTy flags ) const 
{
Ss ss; 
switch (flags)
{
default:
MM_SZ_LOOP(i, m_fwd,sz) { ss<<" "<<i<<":"<<m_fwd[i]; } 
}; // switch 
return ss.str();
}

////////////////////////////////////////////////////////
// command block

static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
// when exiting the interpretter
void clean_up() { m_done=true; } void about()
{
Ss ss;
ss<<" mjm_trees_and_tables "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();
}
void config_banner() { MM_MSG(" configuration "<<m_flp.to_string()) }
bool done() const  { return m_done; } 
void  dump_dmel(OsTy & os )  const
{
if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}
/*
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
*/

private:
void Init()
{
m_done=false;
}
template <class Ty> void Next(const Ty & val)
{
const IdxTy ord=m_fwd.size();
m_fwd.push_back(val);
m_rev.add(val,ord);
}
/*
void DMel(const StrTy & e) { }
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
Dmel * m_dmel;
Ss m_ss;
RaggedMap m_ragged_map;

FwdTy m_fwd;
RevMap m_rev;

}; //mjm_string_ordering


////////////////////////////////////////////////////////////////////////

class mjm_uint_ordering 
{
typedef  ordering_traits::Tr  Tr;
typedef mjm_string_ordering Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;


typedef ordering_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef IdxTy OrdTy;
typedef IdxTy TgtTy;

typedef std::vector<TgtTy> FwdTy;
typedef std::vector<OrdTy> OrderVec;
typedef FwdTy::const_iterator CItor;

typedef mjm_1_to_some<TgtTy,OrdTy, OrderVec> RevMap;

public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
typedef OrdTy order_t;
typedef OrderVec order_vector_t;
typedef TgtTy target_t;
typedef CItor const_iterator;

mjm_uint_ordering():m_dmel(new Dmel()) {Init();}
~mjm_uint_ordering() { delete m_dmel; }

template <class Ty> void add(const Ty & v) { MM_LOOP(ii,v) { Next((*ii)); } } // add
template <class Ty> void add_unique(const Ty & v) { MM_LOOP(ii,v) { 
if (m_rev[*ii].size()==0) 
Next((*ii)); } } // add

void add_unique_value(const TgtTy & v) {  if (m_rev[v].size()==0) Next(v);  } // add


void next() { Next(); } 
// this fails when in and out are the same types... 
const TgtTy & operator[](const OrdTy i ) const  { return m_fwd[i]; } 
// int to int ambugious 
//const order_vector_t  & operator[](const TgtTy &s) const  { return m_rev[s]; } 
const order_vector_t  & operator()(const TgtTy &s) const  { return m_rev[s]; } 
//const order_vector_t   starting_with(const TgtTy &s) const  { return m_rev.starting_with(s); } 


template <class Ts> void sort(Ts & s)
{
FwdTy x=m_fwd;
std::sort(x.begin(),x.end(),s);
MM_ERR(" NEED TO FIX REV MAP")
//m_fwd.clear();
clear();
MM_LOOP(ii,x) Next(*ii);
}
void clear() { m_fwd.clear(); m_rev.clear(); } 

StrTy dump(const IdxTy flags ) const  { return to_string(flags); } 
StrTy to_string(const IdxTy flags ) const 
{
Ss ss; 
switch (flags)
{
default:
MM_SZ_LOOP(i, m_fwd,sz) { ss<<" "<<i<<":"<<m_fwd[i]; } 
}; // switch 
return ss.str();
}

////////////////////////////////////////////////////////
// command block

 IdxTy size()const  {  return m_fwd.size(); } 
static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
// when exiting the interpretter
void clean_up() { m_done=true; } void about()
{
Ss ss;
ss<<" mjm_trees_and_tables "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();
}
void config_banner() { MM_MSG(" configuration "<<m_flp.to_string()) }
bool done() const  { return m_done; } 
void  dump_dmel(OsTy & os )  const
{
if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}
/*
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
*/

private:
void Init()
{
m_done=false;
}

template <class Ty> void Next(const Ty & val)
{
const IdxTy ord=m_fwd.size();
m_fwd.push_back(val);
m_rev.add(val,ord);
}
 void Next()
{
const IdxTy ord=m_fwd.size();
m_fwd.push_back(ord);
m_rev.add(ord,ord);
}



/*
void DMel(const StrTy & e) { }
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
Dmel * m_dmel;
Ss m_ss;
RaggedMap m_ragged_map;

FwdTy m_fwd;
RevMap m_rev;

}; //mjm_string_ordering












#undef MM_DMEL 
#endif

