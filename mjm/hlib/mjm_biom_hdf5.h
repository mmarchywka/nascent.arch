#ifndef MJM_BIOM_HDF5_H__
#define MJM_BIOM_HDF5_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
// some of these pickup the math libs... 
// fraction input support 
// for now just code here, not too hard 
// need this for Alfredo lol
#include "mjm_rational.h"
 // #include "mjm_generic_iterators.h"
// not really used but tyedefed
#include "mjm_block_matrix.h"
//#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
#include "mjm_calendar.h"
// get the taxa names conforming 
#include "mjm_strings.h"
#include "mjm_collections.h"
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


#define HAVE_RAPIDJSON
#ifdef HAVE_RAPIDJSON
#include "rapidjson/reader.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/error/en.h"
#include <cstdio>



#endif

/*

2065  wget -O xxx.biom  -S -v "ftp://ftp.microbio.me/emp/release1/otu_tables/picrust/emp_cr_gg_13_8.normalized.biom"

2098  h5dump xxx.biom | grep "[a-zA-Z]" > xxx.data
 2099  more xxx.data
 2100  vi xxx.data
 2101  grep -v "__"  xxx.data | more
 2102  grep -v "__"  xxx.data | grep -v "([0-9]*): "
 2103  grep -v "__"  xxx.data | more
 2104  grep  "__"  xxx.data | more
 2105  grep -v "__"  xxx.data | more
 2106  grep -v "__"  xxx.data | grep -i soil
 2107  more xxx.data
 2108  h5dump xxx.biom | grep -B 5 -A 100 indptr

HDF5 "xxx.biom" {
GROUP "/" {
   ATTRIBUTE "creation-date" {
      DATATYPE  H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
      DATASPACE  SCALAR
      DATA {
      (0): "2016-08-23T16:21:22.893878"
   ATTRIBUTE "format-url" {
      DATATYPE  H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
      DATASPACE  SCALAR
      DATA {
      (0): "http://biom-format.org"
   ATTRIBUTE "format-version" {
      DATATYPE  H5T_STD_I64LE
      DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
--More--(0%)


g++ -DTEST_SNACK__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_snacks.h
*/

////////////////////////////////////////////////////////////////

class biom_hdf5_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
biom_hdf5_params( const StrTy & nm) : Super(nm) {}
biom_hdf5_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  

//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
//bool always_dump_to_file() const { return m_map.get_bool("always_dump_to_file",!true); }

//bool skip_old() const { return m_map.get_bool("skip_old",true); }
//bool print_dog_days() const { return m_map.get_bool("print_dog_days",!true); }
//bool accumulate_dog_days() const { return m_map.get_bool("accumulate_dog_days",true); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
//bool print_conflicts() const { return m_map.get_bool("print_conflicts",!true); }

//StrTy snack_log() const { return m_map.get_string("snack_log","snacks_log"); }
		//if (skip_old) if (date==StrTy("2017-04-22")) skipping=false;
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }

IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

//ss<<"time_step="<<time_step()<<sep;
//ss<<"skip_old="<<skip_old()<<sep;
//ss<<"print_dog_days="<<print_dog_days()<<sep;
//ss<<"accumulate_dog_days="<<accumulate_dog_days()<<sep;
ss<<"log_commands"<<log_commands()<<sep;
//ss<<"print_conflicts"<<print_conflicts()<<sep;
//ss<<"snack_log="<<snack_log()<<sep;
//ss<<"start_date="<<start_date()<<sep;
//ss<<"end_date="<<end_date()<<sep;
ss<<"otu_format="<<otu_format()<<sep;
return ss.str();
}


}; // biom_hdf5_params


namespace biom_hdf5_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // biom_hdf5_traits


#if 0 


#endif


#ifdef HAVE_RAPIDJSON
namespace { 

using namespace rapidjson;
using namespace std;
class biom_json 
{
typedef biom_json Myt;
typedef biom_hdf5_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef Tr::Dmel  Dmel; 
public:
biom_json() {}

struct MyBlankHandler {
    bool Null() { return true; }
    bool Bool(bool b) {  return true; }
    bool Int(int i) {  return true; }
    bool Uint(unsigned u) { return true; }
    bool Int64(int64_t i) {  return true; }
    bool Uint64(uint64_t u) {  return true ; }
    bool Double(double d) {  return true; }
    bool RawNumber(const char* str, SizeType length, bool copy) { return true; }
    bool String(const char* str, SizeType length, bool copy) { return true; }
    bool StartObject() { return true; }
    bool Key(const char* str, SizeType length, bool copy) { return true; }
    bool EndObject(SizeType memberCount) {  return true; } 
	bool StartArray() { return true; }
    bool EndArray(SizeType elementCount) {  return true; }
}; // MyBlankHandler


class MyVirtHandler {
public:
virtual ~MyVirtHandler() {}

virtual bool Null() {return true;}
virtual bool Bool(bool b) {return true;}
virtual bool Int(int i) {return true;}
virtual bool Uint(unsigned u) {return true;}
virtual bool Int64(int64_t i) {return true;}
virtual bool Uint64(uint64_t u) {return true ;}
virtual bool Double(double d) {return true;}
virtual bool RawNumber(const char* str, SizeType length, bool copy) {return true;}
virtual bool String(const char* str, SizeType length, bool copy) {return true;}
virtual bool StartObject() {return true;}
virtual bool Key(const char* str, SizeType length, bool copy) {return true;}
virtual bool EndObject(SizeType memberCount) {return true;}
virtual bool StartArray() {return true;}
virtual bool EndArray(SizeType elementCount) {return true;}


}; // MyVirtHandler


class MyPtrHandler : public MyVirtHandler  {

// this needs some concept of state 
typedef std::vector<StrTy> StateVec;
public:
MyPtrHandler() :m_debug(true), m_ptr(0),m_array(0),m_obj(0) {}
virtual ~MyPtrHandler() {}
//MyPtrHandler() : m_ptr(this),m_array(0),m_obj {}
MyPtrHandler(MyVirtHandler * p) : m_ptr(p),m_array(0),m_obj(0) {}
virtual bool Null(){m_ptr->Null(); return true;}
virtual bool Bool(bool b){m_ptr->Bool(b); return true;}
virtual bool Int(int i){m_ptr->Int(i); return true;}
virtual bool Uint(unsigned u){m_ptr->Uint(u); return true;}
virtual bool Int64(int64_t i){m_ptr->Int64(i); return true;}
virtual bool Uint64(uint64_t u){m_ptr->Uint64(u); return true ;}
virtual bool Double(double d){m_ptr->Double(d); return true;}
virtual bool RawNumber(const char* str, SizeType length, bool copy){m_ptr->RawNumber(str,length,copy); return true;}
virtual bool String(const char* str, SizeType length, bool copy){m_ptr->String(str,length,copy); return true;}
virtual bool StartObject(){CheckPush(); ++m_obj; m_ptr->StartObject(); return true;}
virtual bool Key(const char* str, SizeType length, bool copy){NewKey(str,length,copy);  return true;}
virtual bool EndObject(SizeType memberCount){m_ptr->EndObject(memberCount); --m_obj;  CheckPop();  return true;}
virtual bool StartArray(){CheckPush(); ++m_array; m_ptr->StartArray(); return true;}
virtual bool EndArray(SizeType elementCount){m_ptr->EndArray(elementCount); --m_array; CheckPop();  return true;}

void NewKey ( const char * str, SizeType l, bool c)
{
m_key=str;
//CheckPop();
m_ptr->Key(str,l,c);
}
void CheckPush()
{

}
void CheckPop()
{
if (m_stack.size()==0) return;
auto & s=m_stack.back();
if (s.array==m_array) if (s.obj==m_obj) {

if (m_debug) {MM_ERR(" popping handler ptr  ")}

 Load(s); m_stack.pop_back(); } 
}

class St { public: MyVirtHandler * p; StrTy key; IdxTy array,obj; };
void push(MyVirtHandler *p ) { Push(); m_ptr=p; } 
void Push()
{
St st;
st.p=m_ptr;
st.key=m_key;
st.array=m_array;
st.obj=m_obj;
m_stack.push_back(st);
}
void Load(const St & st)
{
m_ptr=st.p;
m_key=st.key;
m_array=st.array;
m_obj=st.obj;
}
bool m_debug;
MyVirtHandler * m_ptr;
//std::vector< MyVirtHandler * > m_stack;
std::vector< St > m_stack;
StrTy m_key;
StateVec m_state;
IdxTy m_array,m_obj;
}; // MyPtrHandler

class biom_rows : public MyVirtHandler
{
public:
typedef std::vector<StrTy> Taxa;
typedef string_tokenizer St;
typedef std::vector<IdxTy> TokenTaxa;
//class  row  { public: void clear() { taxa.clear(); } StrTy id; Taxa taxa; } ; 
class  row  { public: 
	void clear() { taxa.clear(); id="uninitialized"; } 
	IdxTy size() const { return taxa.size(); }
Taxa taxav( St & st) const { Taxa x; st.tokenize_vector(x,taxa); return x; }   
	 StrTy id; 

// fudd private:
TokenTaxa taxa; 


} ; 
typedef std::vector<row> Data; // m_rows;
biom_rows():m_debug(true),m_array(0),m_obj(0),m_st(0) {}
biom_rows(St * st ):m_debug(true),m_array(0),m_obj(0),m_st(st) {}
virtual ~biom_rows() {}

IdxTy size() const { return m_rows.size(); }
Data & data() { return m_rows; }
virtual bool Null() {return true;}
virtual bool Bool(bool b) {return true;}
virtual bool Int(int i) {return true;}
virtual bool Uint(unsigned u) {return true;}
virtual bool Int64(int64_t i) {return true;}
virtual bool Uint64(uint64_t u) {return true ;}
virtual bool Double(double d) {return true;}
virtual bool RawNumber(const char* str, SizeType length, bool copy) {return true;}
virtual bool String(const char* str, SizeType length, bool copy) {ValueS(str); return true;}
virtual bool StartObject() {++m_obj; return true;}
virtual bool Key(const char* str, SizeType length, bool copy) {m_key=str; return true;}
virtual bool EndObject(SizeType memberCount) {--m_obj; return true;}
virtual bool StartArray() {++m_array; return true;}
virtual bool EndArray(SizeType elementCount) {--m_array; CheckPush(); return true;}
void CheckPush()
{
if (m_array==1) { 
if ((m_rows.size()%10000)==0) 
	{ MM_ERR(" pushing row "<<MMPR3(m_rows.size(),m_id,m_taxa.size()))}
row r;
if (m_id=="alreadypushed")
{
MM_ERR(" parsing a stale id "<<MMPR2(m_rows.size(),m_taxa.size()))
}
r.id=m_id;
if (r.id.length()<1) 
{ 
	StrTy info=""; 
	if (m_taxa.size()>0) info=m_taxa[m_taxa.size()-1]; 
	MM_ERR(" blank row id detected "<<MMPR3(m_rows.size(),m_taxa.size(),info))
}
//r.taxa=m_taxa;
(*m_st).tokenize_vector(r.taxa,m_taxa);

m_rows.push_back(r); m_taxa.clear(); m_id="alredypushed"; 
} 

}
void ValueS(const char * s) 
{
if (m_key=="id") {m_id=s; } 
else if (m_key=="taxonomy" )  { m_taxa.push_back(s); } 
else {MM_ERR(" bad key in rowy data "<<MMPR2(m_key,s)) } 
}
bool m_debug;
IdxTy m_array,m_obj;
StrTy m_key;
StrTy m_id;
Taxa m_taxa;
//std::vector<row> 
Data m_rows;
St*  m_st;
} ; // biom_rows


/////////////////////////////////////////////////////////////j
class biom_columns : public MyVirtHandler
{
//typedef std::vector<StrTy> Taxa;
//class  row  { public: void clear() { taxa.clear(); } StrTy id; Taxa taxa; } ; 
public:
typedef std::map<StrTy, StrTy > Hash;
class col { public: void clear() { map.clear(); } Hash map; }  ;
typedef std::vector<col>  Data;
biom_columns():m_debug(true),m_array(0),m_obj(0) {}
virtual ~biom_columns() {}

IdxTy size() const { return m_cols.size(); }
Data & data() { return m_cols; } 
virtual bool Null() {Value(StrTy("null")); return true;}
virtual bool Bool(bool b) {Value(b); return true;}
virtual bool Int(int i) {Value(i); return true;}
virtual bool Uint(unsigned u) {Value(u); return true;}
virtual bool Int64(int64_t i) {Value(i); return true;}
virtual bool Uint64(uint64_t u) {Value(u); return true ;}
virtual bool Double(double d) {Value(d); return true;}
virtual bool RawNumber(const char* str, SizeType length, bool copy) 
	{Value(str); return true;}
virtual bool String(const char* str, SizeType length, bool copy) 
	{Value(str); return true;}
virtual bool StartObject() {++m_obj; return true;}
virtual bool Key(const char* str, SizeType length, bool copy) {m_key=str; return true;}
virtual bool EndObject(SizeType memberCount) {--m_obj; CheckPush(); return true;}
virtual bool StartArray() {++m_array; return true;}
virtual bool EndArray(SizeType elementCount) {--m_array; CheckPush();
 return true;}
void CheckPush()
{
//MM_ERR(MMPR3(m_array,m_obj,m_col.map.size()))
if (m_obj==~0)
{
MM_ERR(" should fix the convention "<<MMPR(m_cols.size()))
MM_LOOP(ii,m_cols[0].map)
{
MM_ERR(MMPR2((*ii).first,(*ii).second))
}
}
if (m_obj==0)
{
	if ((m_cols.size()%5000)==0) 
	{
		MM_ERR(" columns "<<MMPR3(m_cols.size(),m_col.map["id"],m_col.map.size()))
		MM_LOOP(ii,m_col.map)
		{
			MM_ERR(MMPR2((*ii).first,(*ii).second))
		}
	}
	// this shoould be the test instead of object???
	if (m_array==1)
	{
		if (m_cols.size()==0) { MM_ERR(" pushing null biom column "<<MMPR(m_cols.size()) )  } 
		m_cols.push_back(m_col);
		m_col.clear();
	}
}
/*
if (m_array==1) { 
if ((m_rows.size()%10000)==0) 
	{ MM_ERR(" pushing row "<<MMPR3(m_rows.size(),m_id,m_taxa.size()))}
row r;
r.id=m_id;
r.taxa=m_taxa;
m_rows.push_back(r); m_taxa.clear(); } 
*/
}
template <class Tv > void Value( const Tv & v)
{
	Ss ss;
	ss<<v;
	if (m_col.map.find(m_key)!=m_col.map.end()) 
		{ MM_ONCE(" appending col or otu data to map "<<MMPR(m_key), ) } 
	//m_col.map[m_key]=ss.str();
	m_col.map[m_key]+=ss.str();
}
void ValueS(const char * s) 
{
/* if (m_key=="id") {m_id=s; } 
else if (m_key=="taxonomy" )  { m_taxa.push_back(s); } 
else {MM_ERR(" bad key in rowy data "<<MMPR2(m_key,s)) } 
*/
}
bool m_debug;
IdxTy m_array,m_obj;
StrTy m_key;
//StrTy m_id;
//Taxa m_taxa;
//std::vector<row> m_rows;
col m_col;
Data m_cols;
} ; // biom_rows


//////////////////////////////////////////////////////////////
// sparse
class biom_data : public MyVirtHandler
{
public:
typedef std::ostream Dest;
class Tuple { public: IdxTy i,j; D v; };
typedef std::vector<Tuple> Data;
biom_data():m_debug(true), m_array(0),m_idx(0),m_discard_data(false),m_dest(0) {}
virtual ~biom_data() {}
void discard_data() { m_discard_data=true; } 
void dest() { m_dest=0; }
void dest(Dest * p ) { m_dest=p;  }

IdxTy size() const { return m_data.size(); }
virtual bool Null() {return true;}
virtual bool Bool(bool b) {return true;}
virtual bool Int(int i) {Value(i); return true;}
virtual bool Uint(unsigned u) {Value(u); return true;}
virtual bool Int64(int64_t i) {Value(i); return true;}
virtual bool Uint64(uint64_t u) {Value(u); return true ;}
virtual bool Double(double d) {Value(d); return true;}
virtual bool RawNumber(const char* str, SizeType length, bool copy) {
if (m_debug) { MM_ERR(" string for num "<<MMPR3(m_t.i,m_t.j,str)) } 
Value(atof(str));
return true;}
virtual bool String(const char* str, SizeType length, bool copy) {return true;}
virtual bool StartObject() {return true;}
virtual bool Key(const char* str, SizeType length, bool copy) {return true;}
virtual bool EndObject(SizeType memberCount) {return true;}
virtual bool StartArray() {++m_array; SetIdx(); return true;}
virtual bool EndArray(SizeType elementCount) {--m_array; CheckDun();  return true;}
template <class Tv> 
void Value(const Tv & x )
{
switch ( m_idx)
{
case 0:{ m_t.i=x; break;  } 
case 1:{ m_t.j=x; break; } 
case 2:{ m_t.v=x; break; } 
default:{ MM_ERR(" bad fild in biom data "<<MMPR4(m_idx,m_t.i,m_t.j,m_t.v)) } 
};
++m_idx;

}
void Dispatch()
{
m_idx=0;
if (m_dest==0){    m_data.push_back(m_t); return; } 
const StrTy sep=" ";
(*m_dest)<<m_t.i<<sep<<m_t.j<<sep<<m_t.v<<CRLF;

}
void CheckDun()
{
if (m_array==1){ if (m_idx==3){
if (!m_discard_data) Dispatch();//   m_data.push_back(m_t);
else {MM_ONCE(" data discard activated ... ",)
MM_DUNCE(" discard ",1000000)
}
if (m_debug&&!m_discard_data&&(m_dest==0))if (( m_data.size()%1000000)==0) 
{ MM_ERR( " push back "<<MMPR4(m_data.size(), m_t.i,m_t.j,m_t.v)) } 
}
else { MM_ERR(" incomplete "<<MMPR4(m_array, m_idx,m_t.i,m_t.j)<<MMPR(m_t.v)) } 
}
}
void SetIdx()
{
if (m_array==2) m_idx=0;
}
bool m_debug;
IdxTy m_array,m_idx;
Tuple m_t;
Data m_data;
bool m_discard_data;
Dest* m_dest;
} ; // biom_data




class biom_handler : public  MyVirtHandler 
{
typedef std::map<StrTy, StrTy> InfoMap;

typedef std::map<StrTy, MyVirtHandler * > VectorMap;
typedef std::ostream Dest;

typedef string_tokenizer St;
public:
biom_handler():m_debug(true), m_virt(new MyPtrHandler(this)),m_rows(&m_st),m_dest(0) {Init(); }
virtual ~biom_handler() { delete m_virt;
if (m_vft.size()!=0) {MM_ERR(" result stil in m_vft lol") } 
}
void discard_data() { m_data.discard_data(); } 
void dest() { m_dest=0;m_data.dest(); }
void dest(Dest * p ) { m_dest=p; m_data.dest(p);  }

virtual bool Key(const char* str, SizeType length, bool copy) {CheckVector(str); return true;}
virtual bool Null() {return true;}
virtual bool Bool(bool b) {return Value(b);}
virtual bool Int(int i) {return Value(i);}
virtual bool Uint(unsigned u) {return Value(u);}
virtual bool Int64(int64_t i) {return Value(i);}
virtual bool Uint64(uint64_t u) {return Value(u);}
virtual bool Double(double d) {return Value(d);}
virtual bool RawNumber(const char* str, SizeType length, bool copy) {return Value(str);}
virtual bool String(const char* str, SizeType length, bool copy) {return Value(str);}
virtual bool StartObject() {return true;}
virtual bool EndObject(SizeType memberCount) {return true;}
virtual bool StartArray() {return true;}
virtual bool EndArray(SizeType elementCount) {return true;}
MyPtrHandler * ptr()  { return m_virt; } 
StrTy string() const
{
Ss ss;
MM_LOOP(ii, m_map)
{
ss<<( MM_STR_LOC ) <<MMPR2((*ii).first,(*ii).second)<<CRLF;
} // ii 
ss<<MMPR(m_data.size())<<CRLF;
ss<<MMPR(m_rows.size())<<CRLF;
ss<<MMPR(m_cols.size())<<CRLF;

return ss.str();
}
typedef biom_data::Data data_type;
typedef biom_rows::Data rows_type;
typedef biom_columns::Data cols_type;
biom_data::Data & data() { return m_data.m_data; } 
biom_rows::Data & rows() { return m_rows.m_rows; } 
biom_columns::Data & cols() { return m_cols.m_cols; } 
St & tokenizer() { return m_st; } 
private:
void Init()
{
m_vft["rows"] = &m_rows; // new biom_rows();
m_vft["columns"] = &m_cols; // new biom_columns();
//m_map["columns"] = new biom_columns();
m_vft["data"] = &m_data; // new biom_data();

}
template <class Tv > bool Value(const Tv & v) 
{
Ss ss;
ss<<v;
if (m_debug) { MM_ERR(" mapping value "<<MMPR2((m_virt->m_key),(ss.str()))) } 
auto & x=m_map[m_virt->m_key];
if (x.length()!=0) x+=StrTy(" ")+ss.str(); else x=ss.str();

return true;
}
// a new key 
void CheckVector(const char * str)
{
auto p = m_vft.find(str);
if (p== m_vft.end()) return;
if (m_debug) {MM_ERR(" pushing "<<MMPR(str))}
m_virt->push((*p).second);  

}
bool m_debug;
MyPtrHandler * m_virt;
InfoMap m_map;
St m_st;
VectorMap m_vft;
biom_rows m_rows;
biom_columns m_cols;
biom_data m_data;
Dest* m_dest;

}; // biom_handler


class MyInspectHandler {
	public:
MyInspectHandler(): m_array(0),m_obj(0) {}
    bool Null() { out("null"); return true; }
    bool Bool(bool b) {out("bool",b);   return true; }
    bool Int(int i) { if (!silent()){ out("int",i);}   return true; }
    bool Uint(unsigned u) { if (!silent()){ out("uint",u);}  return true; }
    bool Int64(int64_t i) { if (!silent()) { out("i64",i);} return true; }
    bool Uint64(uint64_t u) { if (!silent()){ out("ui64",u);}  return true ; }
    bool Double(double d) { if (!silent()) {  out("double",d);}  return true; }
    bool RawNumber(const char* str, SizeType length, bool copy) { out("rn",str,length); return true; }
    bool String(const char* str, SizeType length, bool copy) {out("str",str,length);  return true; }
    bool StartObject() { ++m_obj;out("startobj"); return true; }
    bool Key(const char* str, SizeType length, bool copy) {out("key",str,length);  return true; }
    bool EndObject(SizeType memberCount) {  out("endobj",memberCount,m_obj); --m_obj; return true; } 
	bool StartArray() { ++m_array; if (! silent()) { out("startarray");  } 
return true; }
    bool EndArray(SizeType elementCount) { if (!silent()){  out("endarray",elementCount); }  
 --m_array; return true; }

template <class Tx,class Ty,class Tz> 
void out(const Ty & x, const Tz & y,const Tx & z) 
{
Ss ss;
ss<<x<<" "<<y<<" "<<z<<" arrays "<<m_array;
outs(ss);
}

template <class Tx,class Ty> 
void out(const Ty & x, const Tx & y) 
{
Ss ss;
ss<<x<<" "<<y;
outs(ss);
}

template <class Ty> 
void out(const Ty & x) 
{
Ss ss;
ss<<x;
outs(ss);
}

bool silent() const { return !true; } // return m_array!=0; } 
void outs(Ss & ss)
{
std::ostream & os=std::cout;
//if (!silent())
 os<<ss.str()<<CRLF;

}

IdxTy m_array;
IdxTy m_obj;

}; // MyInspectHandler
//////////////////////////////////////////////////////////

class JSON_acceptor
{
public:
    virtual bool Bool(bool b) {   return true; }
    virtual bool Int(int i) {    return true; }
    virtual bool Uint(unsigned u) {  return true; }
    virtual bool Int64(int64_t i) { return true; }
    virtual bool Uint64(uint64_t u) {  return true ; }
    virtual bool Double(double d) {   return true; }
    virtual bool RawNumber(const char* str, SizeType length, bool copy) { return true; }
    virtual bool String(const char* str, SizeType length, bool copy) { return true; }


}; // JSAON_acceptor


class MyAccHandler {
	public:
MyAccHandler(): m_array(0),m_obj(0) {}
    bool Null() { out("null"); return true; }
    bool Bool(bool b) {out("bool",b);   return m_ptr->Bool(b);; }
    bool Int(int i) { if (!silent()){ out("int",i);}   return m_ptr->Int(i); }
    bool Uint(unsigned u) { if (!silent()){ out("uint",u);}  return m_ptr->Uint(u); }
    bool Int64(int64_t i) { if (!silent()) { out("i64",i);} return m_ptr->Int64(i);; }
    bool Uint64(uint64_t u) { if (!silent()){ out("ui64",u);}  return m_ptr->Uint64(u); }
    bool Double(double d) { if (!silent()) {  out("double",d);}  return m_ptr->Double(d);  }
    bool RawNumber(const char* str, SizeType length, bool copy) { out("rn",str,length); return m_ptr->RawNumber(str,length,copy) ; }
    bool String(const char* str, SizeType length, bool copy) {out("str",str,length);  return m_ptr->String(str,length,copy); }
    bool StartObject() { ++m_obj;out("startobj"); return true; }
    bool Key(const char* str, SizeType length, bool copy) {out("key",str,length);  return true; }
    bool EndObject(SizeType memberCount) {  out("endobj",memberCount,m_obj); --m_obj; return true; } 
	bool StartArray() { ++m_array; if (! silent()) { out("startarray");  } 
return true; }
    bool EndArray(SizeType elementCount) { if (!silent()){  out("endarray",elementCount); }  
 --m_array; return true; }

template <class Tx,class Ty,class Tz> 
void out(const Ty & x, const Tz & y,const Tx & z) 
{
Ss ss;
ss<<x<<" "<<y<<" "<<z<<" arrays "<<m_array;
outs(ss);
}

template <class Tx,class Ty> 
void out(const Ty & x, const Tx & y) 
{
Ss ss;
ss<<x<<" "<<y;
outs(ss);
}

template <class Ty> 
void out(const Ty & x) 
{
Ss ss;
ss<<x;
outs(ss);
}

bool silent() const { return !true; } // return m_array!=0; } 
void outs(Ss & ss)
{
std::ostream & os=std::cout;
//if (!silent())
 os<<ss.str()<<CRLF;

}

IdxTy m_array;
IdxTy m_obj;
JSON_acceptor * m_ptr;

}; // MyAccHandler






/////////////////////////////////////////////////////////////////////////

struct MyHandler {
    bool Null() { cout << "Null()" << endl; return true; }
    bool Bool(bool b) { cout << "Bool(" << boolalpha << b << ")" << endl; return
 true; }
    bool Int(int i) { cout << "Int(" << i << ")" << endl; return true; }
    bool Uint(unsigned u) { cout << "Uint(" << u << ")" << endl; return true; }
    bool Int64(int64_t i) { cout << "Int64(" << i << ")" << endl; return true; }
    bool Uint64(uint64_t u) { cout << "Uint64(" << u << ")" << endl; return true
; }
    bool Double(double d) { cout << "Double(" << d << ")" << endl; return true; 
}
    bool RawNumber(const char* str, SizeType length, bool copy) { 
        cout << "Number(" << str << ", " << length << ", " << boolalpha << copy 
<< ")" << endl;
        return true;
    }
    bool String(const char* str, SizeType length, bool copy) { 
        cout << "String(" << str << ", " << length << ", " << boolalpha << copy 
<< ")" << endl;
        return true;
    }
    bool StartObject() { cout << "StartObject()" << endl; return true; }
    bool Key(const char* str, SizeType length, bool copy) {
        cout << "Key(" << str << ", " << length << ", " << boolalpha << copy << 
")" << endl;
        return true;
    }
    bool EndObject(SizeType memberCount) { cout << "EndObject(" << memberCount << ")" << endl; return true; } bool StartArray() { cout << "StartArray()" << endl; return true; }
    bool EndArray(SizeType elementCount) { cout << "EndArray(" << elementCount << ")" << endl; return true; }
}; // MyHandler




int main() {
    const char json[] = " { \"hello\" : \"world\", \"t\" : true , \"f\" : false, \"n\": null, \"i\":123, \"pi\": 3.1416, \"a\":[1, 2, 3, 4] } ";

    MyHandler handler;
    Reader reader;
    StringStream ss(json);
    reader.Parse(ss, handler);
return 0; 
}

int parse(const StrTy & fn, const IdxTy ty, const char * dest=0 ) {
    Reader reader;
  	FILE *fp = fopen(fn.c_str(), "rb");
//  	FILE *fdest =0;
	std::ostream * ofs=0;
	if ( dest!=0){
		 ofs= new std::ofstream(StrTy(dest)); // fdest= fopen(dest, "wb");
		MM_ERR(" opening "<<dest)	
	}
	const IdxTy bsize=1<<16;
	char b[bsize];
	FileReadStream  ss(fp,b,bsize);
switch (ty)
{
   case 0:{ MyBlankHandler handler; reader.Parse(ss, handler); break; } 
   case 1:{ MyInspectHandler handler; reader.Parse(ss, handler); break; } 
   case 2:{ MyHandler handler; reader.Parse(ss, handler); break; } 
   case 3:{ MyAccHandler handler; reader.Parse(ss, handler); break; } 
   case 4:{ MyVirtHandler handler; reader.Parse(ss, handler); break; } 
   case 5:{ //biom_handler handler; 
		biom_handler & handler= m_bh;
		handler.dest(ofs);
		reader.Parse(ss, *handler.ptr());
		handler.dest();
		//std::cout<<handler.string()<<CRLF; 
		MM_ERR(handler.string()) 
 		break; 
		} 
   case 6:
		{ 
			biom_handler & handler= m_bh;
			m_bh.discard_data();
			reader.Parse(ss, *handler.ptr());
			MM_ERR(handler.string()) 
 			break; 
		} 


default : {MM_ERR(" bad case "<<MMPR(ty)) } 
}
fclose(fp);
//if (fdest!=0) { fclose(fdest); } 
delete ofs;
return 0; 
}
typedef biom_handler::data_type data_type;
typedef biom_handler::rows_type rows_type;
typedef biom_handler::cols_type cols_type;
data_type & data() { return m_bh.data(); } 
rows_type & rows() { return m_bh.rows(); } 
cols_type & cols() { return m_bh.cols(); } 
StrTy info() { return m_bh.string(); } 
// this needs to be typedefed some fing wehre 
string_tokenizer  & tokenizer() { return m_bh.tokenizer(); } 
biom_handler m_bh;



}; //biom_json

}; // namespace 

#endif



class mjm_biom_hdf5 
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



typedef mjm_biom_hdf5 Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

//typedef mjm_logic_base Logic;
typedef biom_hdf5_params Logic;
typedef mjm_logic_base VariableStore;

typedef std::vector<StrTy> DumpOrder;



typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

typedef mjm_calendar CalTy;

////////////////////////////////////////////////////////
//typedef otu_collection OtuCol;
//typedef std::map<StrTy,OtuCol> OtuMap; 
//typedef seq_sample Sample;
//typedef std::map<StrTy,Sample> SampleMap; 
//typedef std::map<StrTy,SampleMap> SampleGroup; 

////////////////////////////////////////////////////

public :
mjm_biom_hdf5():m_dmel(new Dmel()) {Init();}
mjm_biom_hdf5(int argc,char **_args) : m_dmel(new Dmel())
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<argc; ++i) args[i]=_args[i];
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
~mjm_biom_hdf5()
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
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { ++i; clean_up(); }
if (s=="-about") { ++i; about(); }
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


void dump_unused()
{
Ss ss;
//for (auto ii=m_unused.begin(); ii!=m_unused.end(); ++ii)
{
//ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
MM_MSG("unused:"<<CRLF<<ss.str())

}

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
if (cmd=="about") { about();  continue; } 
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
if (cmd=="parse-biom-json") { parse_biom_json(li.word(1),li.word(2));  continue; } 
if (cmd=="parse-in") { parse_biom(li.word(1),li.word(2));  continue; } 


/*
if (cmd=="dump-otu") { dump_otu(std::cout);  continue; } 
if (cmd=="catalog-otu") { catalog_otu(std::cout);  continue; } 

// this works for  format but notin the frontiers paper
if (cmd=="parse-reads") { parse_reads(li.word(1),li.word(2));  continue; } 
// this is 
if (cmd=="parse-front-reads") 
	{ parse_front_reads(li.word(1),li.word(2));  continue; } 
if (cmd=="dump-reads") { dump_reads(std::cout);  continue; } 
if (cmd=="sum-samples") { sum_samples(li.word(1),atoi(li.word(2).c_str()));  continue; } 
if (cmd=="dump-group-raw") { dump_group(li.word(1),true,false);  continue; } 
if (cmd=="dump-group-collate") { dump_group(li.word(1),!true,!false);  continue; } 

if (cmd=="dump-otu-frac") {  dump_collated_group(std::cout, StrTy("dummy")); continue; } 
*/


//if (cmd=="pair-all") { pair_stats(li.word(1));  continue; } 
//if (cmd=="dose-vectors") { if (sz>1) { dose_vectors(li.word(1));} else { MM_ERR( "need dogname ") }  continue; } 
//void pair_stats(const StrTy & x, const StrTy & y)
//if (cmd=="dump-drugs") { collect_all_drugs(); continue; } 
//if (cmd=="dump-drugs-cerr") { collect_all_drugs(2); continue; } 
//if (cmd=="list-drugs-cerr") { collect_all_drugs(3); continue; } 
//if (cmd=="list-dogs-cerr") { collect_all_drugs(4); continue; } 
//if (cmd=="dump-counts") { stream_counts(std::cout); continue; } 
//if (cmd=="dump-canon-fails") { dump_canon_fails(std::cout); continue; } 
//if (cmd=="clear-canon-fails") { clear_canon_fails(); continue; } 
if (cmd=="dump-dmel") { dump_dmel(std::cout); continue; } 
if (cmd=="dump-dmel-cerr") { dump_dmel(std::cerr); continue; } 
//if (cmd=="dump-dog-days") { dump_dog_days(std::cout,p1); continue; } 
//if (cmd=="dump-dog-times") { dump_dog_times(std::cout,p1); continue; } 
//if (cmd=="n-canon-fails") { MM_MSG(MMPR(m_canon_fails.size())) continue; } 

//if (cmd=="clear-order") { clear_order(); continue; } 
//if (cmd=="dump-order") { dump_order(std::cout); continue; } 
//if (cmd=="push-order") 
		//{ int i=1; while ( li.cmd_ok(i+1)){ push_order(li.word(i)); ++i; } continue; } 
//		{ int i=1; while ( i<sz){ push_order(li.word(i)); ++i; } continue; } 
//if (cmd=="remove-order") 
//		{ int i=1; while ( i<sz){ remove_order(li.word(i)); ++i; } continue; } 

//if (cmd=="load-order") { if( li.cmd_ok(2)) load_order(li.word(1)); continue; } 
//if (cmd=="order-all") { collect_all_drugs(1); continue; } 

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
void about()
{
Ss ss;
ss<<" mjm_biom_hdf5 "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
ss<<" Code to read various files related to 16S rRNA analyses "<<CRLF;
ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;
ss<<"database link provided by Bradley.stevenson@ou.edu"<<CRLF;
ss<<"http://www.earthmicrobiome.org/data-and-code/ "<<CRLF;
ss<<"Sample data provided by echen@zymoresearch.com "<<CRLF;
ss<<" uses rapidjson header library for biom reading "<<CRLF;
ss<<" needs h5dump for the hdf5 reads  "<<CRLF;
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

///////////////////////////////////////////////////////////////////



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

#if 0 
void dump_otu(std::ostream & os)
{
MM_LOOP(ii,m_otus)
{
// each one is an otu
os<<(*ii).second.dump(0,(*ii).first)<<CRLF;
} // ii 
}

void catalog_otu(std::ostream & os)
{
const StrTy sep=" ";
std::map<StrTy, IdxTy> m;
MM_LOOP(ii,m_otus) { (*ii).second.catalog(m); } // ii 

MM_LOOP(ii,m)
{
os<<(*ii).first<<sep<<(*ii).second<<CRLF;
}

}
#endif

enum {  BAD=~0,
ATTRIBUTE=0, CSET=1, CTYPE=2, DATA=3, DATASET=4, DATASPACE=5, 
DATATYPE=6, GROUP=7, H5T_C_S1=8, H5T_CSET_ASCII=9, H5T_IEEE_F64LE=10, 
H5T_STD_I32LE=11, H5T_STD_I64LE=12, H5T_STRING=13, 
H5T_STR_NULLTERM=14, H5T_VARIABLE=15, HDF5=16, SCALAR=17, 
SIMPLE=18, STRPAD=19, STRSIZE=20,LEFT=21, RIGHT=22 };
// keywords finder
template <class Tm > 
IdxTy  find_biom_kw(StrTy & w, IdxTy&  loc, IdxTy & pc, const Tm & m, const Words & words)
{
//MM_LOOP(i,words,sz)
const IdxTy sz=words.size();
for (IdxTy i=pc; i<sz; ++i)
{
++pc;
const auto & v=m.find(words[i]);
if (v!=m.end()) { w=words[i]; loc=i; return (*v).second; } 
}

return BAD ; 
}

template <class Tm > 
void biom_kw_map(  Tm & m)
{
IdxTy seq=0;
m["ATTRIBUTE"]=seq; ++seq; 
m["CSET"]=seq; ++seq; 
m["CTYPE"]=seq; ++seq; 
m["DATA"]=seq; ++seq; 
m["DATASET"]=seq; ++seq; 
m["DATASPACE"]=seq; ++seq; 
m["DATATYPE"]=seq; ++seq; 
m["GROUP"]=seq; ++seq; 
m["H5T_C_S1"]=seq; ++seq; 
m["H5T_CSET_ASCII"]=seq; ++seq; 
m["H5T_IEEE_F64LE"]=seq; ++seq; 
m["H5T_STD_I32LE"]=seq; ++seq; 
m["H5T_STD_I64LE"]=seq; ++seq; 
m["H5T_STRING"]=seq; ++seq; 
m["H5T_STR_NULLTERM"]=seq; ++seq; 
m["H5T_VARIABLE"]=seq; ++seq; 
m["HDF5"]=seq; ++seq; 
m["SCALAR"]=seq; ++seq; 
m["SIMPLE"]=seq; ++seq; 
m["STRPAD"]=seq; ++seq; 
m["STRSIZE"]=seq; ++seq; 
m["{"]=seq; ++seq; 
m["}"]=seq; ++seq; 

}
template <class Ts1, class Ts2> void pop(Ts1 & s1, Ts2 & s2, const IdxTy & n)
{
for(IdxTy i=0; i<n; ++i)
{
s1.pop_back();
s2.pop_back();

}

}
bool isdig(const char & si )
{return (  (si>='0') && (si<='9'))  ; }

StrTy  get_number(const StrTy & w)
{
const IdxTy sz=w.length();
char c[sz+2];
const char * s=w.c_str();
IdxTy pc=0;
for(IdxTy i=0; i<sz; ++i)
{
	const char si=s[i];
if (isdig(si)) { 	c[pc]=si; ++pc; } 
	
}
c[pc]=0;
return StrTy(c);
}
StrTy dequote(const StrTy & w,const bool decomma=false)
{
const IdxTy sz=w.length();
char c[sz+2];
const char * s=w.c_str();
IdxTy pc=0;
for(IdxTy i=0; i<sz; ++i)
{
	const char si=s[i];
	if (si!='"') {c[pc]=si; ++pc; }
	if (decomma) if (si==',') { --pc; } 
}
c[pc]=0;
return StrTy(c);
}

// ( a,b,c,... ) 
template <class Td > void parse_ntuple( Td & v, const StrTy & w)
{
const IdxTy sz=w.length();
char c[sz+2];
const char * s=w.c_str();
IdxTy pc=0;
for(IdxTy i=0; i<sz; ++i)
{
	const char si=s[i];
if (isdig(si)) { continue; } 
if ( i==pc) { pc=i+1; continue;  } 
const IdxTy n=atoi(s+pc);
v.push_back(n);
pc=i+1;	
}
if ( pc<sz)
{
const IdxTy n=atoi(s+pc);
v.push_back(n);


}

} // parse_ntuple
	typedef std::vector< IdxTy> BIV;
	typedef std::vector< StrTy> BRV;
	typedef std::vector< std::vector< StrTy> > BSV;
	typedef std::vector< std::vector< IdxTy > > BTV;

void assign(BRV & d, const StrTy & s, const std::vector<IdxTy > & loc, 
string_tokenizer * ps=0)
{
const IdxTy l=loc[0];
// const IdxTy x=atoi( get_number(s).c_str());
if ( d.size()<=l) d.resize(l+1);
d[l]=s;
}

void assign(BIV & d, const StrTy & s, const std::vector<IdxTy > & loc, 
string_tokenizer * ps=0)
{
const IdxTy l=loc[0];
 const IdxTy x=atoi( get_number(s).c_str());
if ( d.size()<=l) d.resize(l+1);
d[l]=x;
}



void assign(BSV & d, const StrTy & s, const std::vector<IdxTy > & loc, string_tokenizer * ps=0)
{
const IdxTy l=loc[0];
// const IdxTy x=atoi( get_number(s).c_str());
if ( d.size()<=l) d.resize(l+1);
auto & v=d[l];
const IdxTy l2=loc[1];
if ( v.size()<=l2) v.resize(l2+1);
//if (ps==0) 
{ v[l2]=s; }// else {}
}

void assign(BTV & d, const StrTy & s, const std::vector<IdxTy > & loc, string_tokenizer * ps=0)
{
const IdxTy l=loc[0];
// const IdxTy x=atoi( get_number(s).c_str());
if ( d.size()<=l) d.resize(l+1);
auto & v=d[l];
const IdxTy l2=loc[1];
if ( v.size()<=l2) v.resize(l2+1);
//if (ps==0) 
{ v[l2]=(*ps)(s); }// else {}
}


template <class Td > 
void parse_biom_line (Td & d , const Words & w, string_tokenizer * ps=0)
{
std::vector<IdxTy> loc;
MM_SZ_LOOP(i,w,sz)
{
const StrTy & s= w[i];
const IdxTy len=s.length();
if ( len==0) continue;
//MM_ERR(" fudd "<<MMPR2(s,i))
if (s.c_str()[0]=='(') { loc.clear(); parse_ntuple(loc,s); }  
else 
{
// const D x=atof( get_number(s));
if (loc.size()==0) loc.push_back(0);
//assign(d,x,loc);
assign(d,s,loc,ps);
//MM_ERR(" assigning "<<MMPR3(s,loc[0],loc.size())) 
++loc[loc.size()-1];
}

}

} // parse_biom_line

void  dump_biom_tax(const StrTy & fn, BTV & ds,string_tokenizer * st)
{
//std::ostream&  os=std::cout;
std::ofstream  os(fn);
const StrTy sep=" ";
IdxTy i=0; 
MM_LOOP(ii,ds)
{
Ss ss;
auto dsii=(*ii);
ss<<i;
MM_LOOP(jj,dsii)
{
ss<<sep<<dequote((*st)(*jj),true);
} // jj 
os<<ss.str()<<CRLF;
++i; 
} // ii 

} // dump_viom_tax

void  dump_biom_sam(const StrTy & fn, BRV & ds,string_tokenizer * st)
{
//std::ostream&  os=std::cout;
std::ofstream  os(fn);
const StrTy sep=" ";
IdxTy i=0; 
MM_LOOP(ii,ds)
{
Ss ss;
auto dsii=(*ii);
ss<<i<<sep<<dequote(dsii,true);
os<<ss.str()<<CRLF;
++i; 
} // ii 

} // dump_viom_tax



// this must be piped in now from h5dump or it defeats the
// purpose. 
void parse_biom(const StrTy &tax, const StrTy & sample)
{

	typedef std::map<StrTy, IdxTy> Wma;
	//typedef std::vector< std::vector< StrTy> > BSV;
	typedef std::map<StrTy, BIV > DsIntMap;
	//typedef std::map<StrTy, BSV > DsTaxMap;
	typedef std::map<StrTy, BTV > DsTaxMap;
	typedef std::map<StrTy, BRV > DsSamMap;
	DsIntMap ds_int;
	DsTaxMap ds_tax;
	DsSamMap ds_sam;
	string_tokenizer st; 
	const bool debug=!false;
	if (debug) { MM_ERR(" parse_biom "<<MMPR2(tax,sample)) } 
	// FUDD SHOT std::ifstream  isn(m_flp.snack_log().c_str());
	//std::ifstream  isn(fn.c_str());
	std::istream  & isn= std::cin; //  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	//const bool skip_old=m_flp.skip_old();
	//const bool print_dog_days=m_flp.print_dog_days();
	//const bool accumulate_dog_days=m_flp.accumulate_dog_days();
	Wma m;
	biom_kw_map(m);
	StrTy kww;
	StrTy dataset,group;
	CommandInterpretter li(is);
	IdxTy cntt=0;
	std::vector<IdxTy> ps;
	std::vector<StrTy> pv;
	IdxTy sample_at=0;
//	auto & otumap= m_otus[name];
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		IdxTy loc=0;
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		IdxTy pc=0;
		while (pc<sz) { 
		const IdxTy kw= find_biom_kw( kww, loc,pc,m,li.words());
		if (kw!=BAD) { ps.push_back(kw); pv.push_back(kww);
				if (debug) {MM_ERR("pushing "<<kww) } 
 } 
//		if (debug) { MM_ERR(MMPR4(kw,kww,pc,loc))}
		switch (kw)
		{

			case DATASET: { dataset=dequote(li.word(loc+1)); 
			if (debug)  { MM_ERR(" parsing dataset "<<MMPR(dataset))}
				++pc; break; } // DATASET 
// this is hierarchail fck 
			case GROUP: { group=dequote(li.word(loc+1)); 
			if (group=="sample") sample_at=ps.size(); 
			if (debug)  { MM_ERR(" parsing group "<<MMPR(group))}
				++pc; break; } // DATASET 


			case DATA: {  break; } // DATASET 
			case LEFT: { 
// this needs to use up all the data 
//			if (ps.back()==DATA) 
//				{ parse_biom_line(ds_int[dataset], li.line()); pc=sz; break; } 
 break; } // DATASET 
			case RIGHT: {
			if (ps[ps.size()-3]==DATA){ // { parse_biom_data(); break; } 
			// for now get rid of the last one to save memory 
			if (dataset=="taxonomy") { dump_biom_tax(tax,ds_tax[dataset],&st);
				ds_tax[dataset].clear(); dataset="";  } 
			if ((sample_at!=0)&&(dataset=="ids")) 
				{ dump_biom_sam(sample,ds_sam[dataset],&st);
				ds_sam[dataset].clear(); dataset="";  } 

			} // data
 		pop(ps,pv,3); 
	if (ps.size()<sample_at) sample_at=0; 

break; } // DATASET 
			default:
			if (ps[ps.size()-2]==DATA) // { parse_biom_data(); break; } 
			{ 
				if ( dataset=="taxonomy" ) 
				{parse_biom_line(ds_tax[dataset], li.words(),&st); pc=sz;
			if ((cntt%1000)==0) {MM_STATUS(" dataset  "<<dataset<<" "<<ds_tax[dataset].size()) } 
			++cntt;
 break; }
			if ((dataset=="ids") && (sample_at!=0))
			{	
				parse_biom_line(ds_sam[dataset], li.words(),&st); pc=sz;
			if ((cntt%1000)==0) {MM_STATUS(" dataset  "<<dataset<<" "<<ds_sam[dataset].size()) } 
			++cntt;
			break;
			}
			// ignore for now 
			
			if ((cntt%100000)==0) {MM_STATUS(" skipping "<<cntt<<" "<<li.line()) } 
			++cntt;
			//	parse_biom_line(ds_int[dataset], li.words(),&st); 
				pc=sz; break; 

			} 
				if (debug) { MM_ERR(" bad case  "<<MMPR2(kw,kww)<<MMPR((li.line())))} 
		} ; // switch kw 
		} // while pc 
	//	if (debug) { MM_ERR(" fudd "<<MMPR((li.line())))} 
		//otumap.add(li.words());
	} // li.next
} // parse_otu

#if 0 
void dump_reads(std::ostream & os)
{
const StrTy sep=" ";
MM_LOOP(jj,m_sample_groups)
{
// each one is an samle
auto & samples=(*jj).second;
const StrTy & nmm=(*jj).first; 
MM_LOOP(ii,samples)
{ 
//os<<nmm<<sep<<(*ii).second.dump(0,nm+sep+(*ii).first)<<CRLF;
os<<(*ii).second.dump(0,nmm+sep+(*ii).first); // <<CRLF;

} // ii 
} // jj 


}
#endif

void parse_biom_json(const StrTy &tax, const StrTy & sample)
{
#ifdef HAVE_RAPIDJSON
biom_json x;
x.parse(tax,atoi(sample.c_str()));
#else
MM_ERR(" skipping parse_biom_json not compiled with json parser")
#endif

}

#if 0 
// froniers paper format samples on first line 
void parse_front_reads(const StrTy &name, const StrTy &fn)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_reads "<<MMPR2(name,fn)) } 
	// FUDD SHOT std::ifstream  isn(m_flp.snack_log().c_str());
	std::ifstream  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	CommandInterpretter li(is);
	std::vector<StrTy> snames; // names;
	auto & samples= m_sample_groups[name];
	// the first line is a bunch of  sample names  
	// field as a name name  
	if  (li.nextok())
	{
		const IdxTy sz=li.size();
		for(IdxTy i=0; i<sz; ++i) snames.push_back(li.word(i));

	}
	// each following line is a otu  name folled by read counts 
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		//const StrTy  sn=li.word(0);
		const StrTy  otun=li.word(0);
		if (debug) { MM_ERR(" fudd "<<MMPR((li.line())))} 
		//const IdxTy sz=li.words().size();
		for(IdxTy i=1; i<sz; ++i)
		//samples[names[i]]+=::atof(li.word(i).c_str());
		samples[snames[i]].add(otun,::atof(li.word(i).c_str()));
	////	samples.add(li.words());
	} // li.next
} // parse_otu

//  format sequnces in first line 
void parse_reads(const StrTy &name, const StrTy &fn)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_reads "<<MMPR2(name,fn)) } 
	// FUDD SHOT std::ifstream  isn(m_flp.snack_log().c_str());
	std::ifstream  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	//const bool skip_old=m_flp.skip_old();
	//const bool print_dog_days=m_flp.print_dog_days();
	//const bool accumulate_dog_days=m_flp.accumulate_dog_days();
	CommandInterpretter li(is);
	std::vector<StrTy> names;
	// the format is a line of sequence names
	// follow by sample name and read coutn line for eahc sample 
	auto & samples= m_sample_groups[name];
	// the first line is a bunch of otu or sequence names with first 
	// field as a name name  
	if  (li.nextok())
	{
		const IdxTy sz=li.size();
		for(IdxTy i=0; i<sz; ++i) names.push_back(li.word(i));

	}
	// each following line is a sample name folled by read counts 
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		const StrTy  sn=li.word(0);
		//otumap.add(li.line());
		if (debug) { MM_ERR(" fudd "<<MMPR((li.line())))} 
		//const IdxTy sz=li.words().size();
		for(IdxTy i=1; i<sz; ++i)
		//samples[names[i]]+=::atof(li.word(i).c_str());
		samples[sn].add(names[i],::atof(li.word(i).c_str()));
	////	samples.add(li.words());
	} // li.next
} // parse_otu





void sum_samples(const StrTy &group, const IdxTy & level)
{


	auto & samples= m_sample_groups[group];
	// this iterates over all samples within group group 
	MM_LOOP(ii,samples)
	{
//	samples[sn].add(names[i],::atof(li.word(i).c_str()));
	} // ii 
}

//#if 0 
void dump_group(const StrTy &group,const bool do_raw, const bool do_collate, std::ostream  & os=std::cout)
{
	const bool debug=false; 
	const StrTy sep=" ";
	const IdxTy otu_format=m_flp.otu_format();
//	const bool do_raw=false;
//	const bool do_collate=true;
	auto & samples= m_sample_groups[group];
	auto & otumap= m_otus[group];
	MM_ERR(MMPR2(samples.size(),otumap.size()))
	// this iterates over all samples within group group 
	MM_LOOP(ii,samples)
	{
		const StrTy & sn=(*ii).first;
		const  Sample & s= (*ii).second;
		typedef std::map<StrTy, D> OtuCollated;
		OtuCollated collate;
		D hit_total=0;
		MM_LOOP(jj,s)
		{

			const StrTy & otu=(*jj).first;
			// this really needs to sum over the same flats 
			const auto & otuii=(*jj).second; // .flat();
			const auto & otui=otumap[otuii.otu()]; // .flat();
			//const StrTy & flat=otui.flat();
			const StrTy & flat=otui.string(otu_format);
			if (debug) MM_ERR(MMPR3(otuii.otu(),otui.size(),flat))
			const D hits=otuii.hits();
			hit_total+=hits;
			if (hits!=0)
			{ 	
				if (do_collate) collate[flat]+=hits;
				if (do_raw) {
				os<<group<<sep<<sn<<sep<<otuii.otu();
				os<<sep<<hits<<sep<<flat<<CRLF;}

			}

		}
		// if anything was collated dump it now 
		MM_LOOP(kk,collate)
		{
			const auto & flat=(*kk).first; 
			auto  hits  =(*kk).second; 
			if (hit_total!=0) hits=hits/hit_total;
			os<<group<<sep<<sn;
			os<<sep<<hits<<sep<<flat<<CRLF;
		} // collate
//	samples[sn].add(names[i],::atof(li.word(i).c_str()));


	} // ii 
}
#endif

#if 0 
///////////////////////////////////////////////////////
template <class Tv, class Tm> 
void invert_map(Tv & v, const Tm & m)
{
v = Tv(m.size());
MM_LOOP(ii,m)
{
v[(*ii).second]=(*ii).first; 
}
}

template <class Tvv> 
void all_samples( Tvv & v)
{
	MM_LOOP(ii,m_sample_groups)
	{
		std::vector<StrTy>  vv;
		vv.push_back((*ii).first);
		MM_LOOP(jj, (*ii).second)
		{
			vv.push_back((*jj).first);
		}
		v.push_back(vv);	
	} // ii 
} // all_samples 
template <class Tvv> 
StrTy gs_header( const Tvv & v)
{
Ss ss;
StrTy gr="";
const StrTy sep=" ";
const StrTy concat="_";
ss<<"serial"<<sep<<"catagory";
MM_LOOP(ii,v)
{
MM_LOOP(jj,(*ii))
{
if (jj==(*ii).begin()) { gr=*jj; continue; } 
ss<<sep<< gr<<concat<<(*jj);
//sep=" ";
}

}
return ss.str();
}
//#if 0
void dump_collated_group(std::ostream & os,  const StrTy  &group)
{
	const StrTy sep=" ";
	typedef std::map<StrTy, IdxTy> Tcat;
	typedef std::map<IdxTy, std::map< StrTy, std::map< StrTy, D > > >  Tcompos;
	typedef std::vector< std::vector< StrTy > > GSVec;
	GSVec allgs;
	all_samples(allgs);	
	Ss ss;
	ss<<gs_header(allgs);
	ss<<CRLF;	
	Tcompos d;
	Tcat m;
 	phyla_seq( m);
	std::vector<StrTy> groups,imap;
	invert_map(imap,m);
	groups.push_back(group);
	//collate_groups(d,m,groups);
	collate_groups(d,m,allgs);
	MM_LOOP(ii,d)
	{
		const IdxTy&  otu_class= (*ii).first;
		ss<<otu_class<<sep<<imap[otu_class];
		MM_LOOP(jj,allgs)
		{
			const auto&  gr=(*jj);
			StrTy grs="xxx";
			MM_LOOP(kk,gr)
			{
				const StrTy  & sample=(*kk);
				// this is kind of a kluge doh 
				if ( kk==gr.begin()) { grs=*kk; continue; } 
				const D & hits=(*ii).second[grs][sample];
				//ss<<sep<<grs<<sample<<sep<<hits;	
				ss<<sep<<hits;	

			} // kk or sample 

		} // jj or group
		ss<<CRLF;
 	} // ii or class 
os<<ss.str();
}


template <class Tmap, class Tcat, class Tg>
void collate_groups(Tmap & d, const Tcat & m, const Tg &groups)
{
	const bool debug=false; 
	const StrTy sep=" ";
	const IdxTy otu_format=m_flp.otu_format();
// each vector in groups has a group as first element followed
// by samples 
	MM_LOOP(kk,groups)
{
	//const StrTy & group=(*kk);
	const StrTy & group=(*kk)[0];
	auto & samples= m_sample_groups[group];
	auto & otumap= m_otus[group];
	for (auto ll=(++(*kk).begin()); ll!=(*kk).end(); ++ ll)
	{
//	MM_ERR(MMPR2(samples.size(),otumap.size()))
	// this iterates over all samples within group group 
	//MM_LOOP(ii,samples)
	{
		const StrTy & sn=(*ll); // (*ii).first;
		const  Sample & s= samples[sn] ; // (*ii).second;
		typedef std::map<StrTy, D> OtuCollated;
		OtuCollated collate;
		D hit_total=0;
		MM_LOOP(jj,s)
		{

			const StrTy & otu=(*jj).first;
			// this really needs to sum over the same flats 
			const auto & otuii=(*jj).second; // .flat();
			const auto & otui=otumap[otuii.otu()]; // .flat();
//			const StrTy & flat=otui.string(otu_format);
			const IdxTy otu_class=otui.classify(m);
//			if (debug) MM_ERR(MMPR3(otuii.otu(),otui.size(),flat))
			const D hits=otuii.hits();
			hit_total+=hits;
			if (hits!=0)
			{ 	
				d[otu_class][group][sn]+=hits;

			}

		}
		// if anything was collated dump it now 
		for (IdxTy otu_class=0; otu_class<m.size(); ++otu_class)
		{
			if (hit_total!=0) //hits=hits/hit_total;
				d[otu_class][group][sn]/=hit_total;
		} // collate

} 
	} // ii 
} // kk 
} // collate_groups 




/////////////////////////////////////////////////////
// for cats that are not at same level, archaea for example,
// the otu needs to figure out which cat it is in based
// on a hit somehwere in the hieratchy 
template <class Tmap> void phyla_seq(Tmap & m)
{
IdxTy seq=0;
m["blank"]=seq; ++seq; 
m["acidobacteria"]=seq; ++seq; 
m["actinobacteria"]=seq; ++seq; 
m["bacteroidetes"]=seq; ++seq; 
m["candidate"]=seq; ++seq; 
m["candidatedivisionjs1"]=seq; ++seq; 
m["candidatedivisionop3"]=seq; ++seq; 
m["candidatedivisionop8"]=seq; ++seq; 
m["candidatedivisionop9"]=seq; ++seq; 
m["chlorobi"]=seq; ++seq; 
m["chloroflexi"]=seq; ++seq; 
m["cyanobacteria"]=seq; ++seq; 
m["deferribacteres"]=seq; ++seq; 
m["euryarchaeota"]=seq; ++seq; 
m["fibrobacteres"]=seq; ++seq; 
m["firmicutes"]=seq; ++seq; 
m["fusobacteria"]=seq; ++seq; 
m["gemmatimonadetes"]=seq; ++seq; 
m["hydrogenedentes"]=seq; ++seq; 
m["lentisphaerae"]=seq; ++seq; 
m["planctomycetes"]=seq; ++seq; 
m["proteobacteria"]=seq; ++seq; 
m["saccharibacteria"]=seq; ++seq; 
m["spirochaetae"]=seq; ++seq; 
m["synergistetes"]=seq; ++seq; 
m["tenericutes"]=seq; ++seq; 
m["thermodesulfobacteria"]=seq; ++seq; 
m["thermotogae"]=seq; ++seq; 
m["unidentified"]=seq; ++seq; 
m["verrucomicrobia"]=seq; ++seq; 
}

#endif

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
bool m_done;
//BraTy m_tree;
//FlpTy m_flp; // right now this is being SET in Init not being used to set m_h

//TimeLines m_timelines;
//KeyCounts m_unused;

//MyChart m_chart;
//time_map m_time_map;
//IdxTy m_ord_col; // input column containing a time ordinal (such as day number)
//IdxTy m_t_col; // input col with absolute time string such as date.
Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;

//OtuMap m_otus;
//SampleGroup m_sample_groups;


CounterMap m_cm;

}; //mjm_biom_hdf5



/////////////////////////////////////////////////////////

#ifdef  TEST_BIOM_HDF5__
int main(int argc,char **args)
{
typedef mjm_biom_hdf5  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


/////////////////////////////////////////////////////
#ifdef  TEST_BF_FLOW__
int main(int argc,char **args)
{
typedef mjm_gauss  Myt;
typedef double D;
typedef unsigned int IdxTy;
Myt x(argc,args);
const D k=atof(args[1]);
const D l=1;
const IdxTy npts=6;
const unsigned int n=atoi(args[2]); // 100;
std::vector<D> ni(npts);
MM_MSG(MMPR(k))
ni[npts/2]=1;
x.do_brute_force_sim( k,  l,  ni,  n );


//if (!x.done()) x.command_mode();

return 0;
}

#endif


#ifdef  TEST_GAUSS__

int main(int argc,char **args)
{
typedef mjm_gauss  Myt;

Myt x(argc,args);

if (!x.done()) x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

