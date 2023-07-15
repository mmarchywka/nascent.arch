
#ifndef MJM_bibtex_entry_H__
#define MJM_bibtex_entry_H__
 
#include "mjm_globals.h"

#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
#include "mjm_strings.h"

//#include "mjm_canned_methods.h"


//#include "mjm_cli_ui.h"

//#include "mjm_tokenized_collections.h"
// uses ragged table for comment parsing 
//#include "mjm_collections.h"

#include "mjm_pawnoff.h"
#include "mjm_blob.h"
#include "mjm_read_buffer.h"
#include "mjm_misc_parse.h"
#include "mjm_bibtex_parse.h"

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
made from wizard script on Fri Sep 13 14:50:58 EDT 2019
toobib
*/


////////////////////////////////////////////////////////////////

namespace bibtex_entry_traits
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


}; // toobib_traits
///////////////////////////////////////////////////////////////

/*


typedef toobib_traits::Tr  Tr;
typedef mjm_toobib Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;

typedef Tr::MyBlock  MyBlock;


//typedef std::vector<StrTy> Words;
//typedef data_model_error_log Dmel;

////////////////////////////////////////////////////


//typedef mjm_pawnoff<Tr>  Handler;
//typedef mjm_blob<Tr> Blob;
//typedef mjm_blob_map<Tr>  BlobMap;
//typedef mjm_misc_parse<Tr> ParseTable;
//typedef mjm_read_buffer<Tr> RdBuf;


*/


/*

@www{$name,
  author =   {}, title =     {$title}, journal =    {}, year = {$year},
  volume =   {}, number =    {}, pages =     {}, note =  {},
urldate={$date},
url={$uin}

*/
template <class Tr> 
class temp_progress
{
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_read_buffer<Tr> RdBuf;
typedef CommandInterpretter Ci;

typedef std::map<StrTy, StrTy> Mmap;
public:

Blob & blob(const StrTy & n) { return m_blobs[n]; } 
// I guess it could find some way to make a ref 
StrTy  get(const StrTy & n) { return StrTy(m_blobs[n]); } 
// this is only from original page... 
Mmap & meta() { return m_meta; } 

void dump(OsTy & os)  const
{
MM_LOOP(ii,m_blobs) 
{ os<<"TempPRogress "<<MMPR2((*ii).first, StrTy((*ii).second))<<CRLF; }
MM_LOOP(ii,m_meta) 
{ os<<"TempPRogressMeta "<<MMPR2((*ii).first, (*ii).second)<<CRLF; }



}

private:

BlobMap m_blobs;
Mmap m_meta;

}; 




template <class Tr> 
class mjm_bib_name_generator
{
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_read_buffer<Tr> RdBuf;



typedef std::map<StrTy, StrTy> Map;
typedef std::vector<StrTy> V;
public:

void append_interesting(StrTy & n, const V & v){ AppendInteresting( n, v); } 
void append_interesting_url(StrTy & n, const StrTy & u){ AppendInterestingU( n, u); } 


void make_words(V & v, const StrTy & s) { MakeWords( v,  s); } 
private:
void MakeWords(V & v, const StrTy & s)
{
IdxTy len=s.length();
char c[len+2];
IdxTy ptr=0;
const IdxTy mask= ParseTable::UC|ParseTable::LC|ParseTable::DIGIT; //  m_clut
for(IdxTy i=0; i<len; ++i)
{
const char cp=s.c_str()[i];
IdxTy lu= m_clut.lut(cp);
if ((lu&mask)!=0) { c[ptr]=cp; ++ptr; }
else{
if (ptr!=0) { c[ptr]=0; StrTy x= StrTy(c); v.push_back(x); 
ptr=0; } 
}
}
if (ptr!=0) { c[ptr]=0; StrTy x= StrTy(c); v.push_back(x); 
}

} //MakeWords
bool Test(const char *  ps, const IdxTy i, const IdxTy mask, const bool p )
{
const char cp=ps[i];
IdxTy lu= m_clut.lut(cp);
const bool x= ((lu&mask)!=0); 
return p?x:!x;
}
void MakeWords(V & v, V & sep, const StrTy & s)
{
IdxTy len=s.length();
char c[len+2];
IdxTy ptr=0;
bool  mode=true;
const IdxTy mask= ParseTable::UC|ParseTable::LC|ParseTable::DIGIT; //  m_clut
const char * ps=s.c_str();
if (len>0) mode=Test(ps,0,mask,true);
for(IdxTy i=0; i<len; ++i)
{
//const char cp=s.c_str()[i];
//IdxTy lu= m_clut.lut(cp);
//if ((lu&mask)!=0) { c[ptr]=cp; ++ptr; }
if (Test(ps,i,mask,mode)) { c[ptr]=ps[i]; ++ptr; }
else{ mode=!mode;
if (ptr!=0) { c[ptr]=0; StrTy x= StrTy(c); (mode?sep:v).push_back(x); 
ptr=1; c[0]=ps[i]; } 
}
}
if (ptr!=0) { c[ptr]=0; StrTy x= StrTy(c); (mode?v:sep).push_back(x); 
}

} //MakeWords



// this needs a routine to mince the word with the Parse thing, 
//enum { BAD=~0,BITS=8*sizeof(Ch), CARD=(1<<BITS),LC =(1<<0), DIGIT=(1<<1),
//WHITE=(1<<2), CTRL=(1<<3), UC=(1<<4), PUNC=(1<<5),EOL=(1<<6),START=(1<<7),
//STOP=(1<<8), OPER=(1<<9) } ;
void AppendInteresting(StrTy & n, const V & v)
{
IdxTy i=0;
IdxTy j=0;
while (i<v.size()) {
if (v[i].length()>3) 
{ n=n+v[i]; ++j; } ++i; if (i>3) break; } ;
} //AppendInteresting 
void AppendInterestingU(StrTy & n, const StrTy & u)
{
V w,s;
MakeWords(w,s,u);
{ MM_SZ_LOOP(i,w,sz) { MM_ERR(MMPR2(i,w[i]))  }}
{ MM_SZ_LOOP(i,s,sz) { MM_ERR(MMPR2(i,s[i]))  }}

{ const IdxTy sz=w.size();
IdxTy i=0; 
for(; i<sz; ++i)
 { 
if (s[i]=="/") { IdxTy l=i-1; if ( l<=w.size()) n=n+w[l]; break; }
 }
if ( i==sz){ IdxTy l=i-1; if ( l<=w.size()) n=n+w[l]; }
}

//AppendInteresting(n,w); 


//IdxTy i=0;
//IdxTy j=0;
//while (i<t.size()) {
//if (t[i].length()>3) 
//{ m_name=m_name+t[i]; ++j } ++i; if (i>3) break; } ;

} //AppendInterestingU 



ParseTable m_clut;
// this should be static or moved to a name generator object
Map m_exclude;
}; // name_generator

template<class Tr>
class mjm_bibtex_entry
{

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;



typedef std::map<StrTy, StrTy> Map;
typedef std::map<StrTy, StrTy> Ord;
typedef std::vector<StrTy> V;
typedef mjm_bib_name_generator<Tr> Ngen;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_read_buffer<Tr> RdBuf;
typedef temp_progress<Tr>  TempProg;
typedef CommandInterpretter Ci;
bool Bit(const IdxTy f, const IdxTy b) { return ((1<<b)&f)!=0; }
bool Mask(const IdxTy f, const IdxTy m) { return ((m)&f)!=0; }


public:
mjm_bibtex_entry() {Init(); } 
IdxTy bad() const  { return ~0; } 
void clear() { Init(); } 
void serial(const IdxTy & x) { m_serial=x;}
IdxTy serial()const  {return  m_serial;}
const StrTy source() const { return m_source_file; } 
void source(const StrTy & s )  { m_source_file=s; } 
void type(const StrTy & x) { m_type=x;}
void name(const StrTy & x) { m_name=x;}
const StrTy &  type() const  { return  m_type;}
const StrTy &  name() const  { return  m_name;}
void comment( const StrTy & key, const StrTy & v) { m_comments[key]=v; } 
// repeated calls now overwrite
void comment(  const StrTy & v) 
//{Ss ss; ss<<m_comments.size();  m_comments[ss.str()]=v; } 
{Ss ss;   m_comments[ss.str()]=v; } 
StrTy comment()
{
if (m_comments.size()==0) return StrTy();
return (*(m_comments.begin())).second;

}
void set_order_vector(const V & order_vector) { m_order_vector=order_vector; } 
void push_order_vector(const StrTy & k) { m_order_vector.push_back(k); } 
void set( const StrTy & key, const StrTy & v) { m_map[key]=v; } 
void seto( const StrTy & key, const StrTy & v) {
if (missing(key)) push_order_vector(key);set(key,v); } 
const StrTy & get( const StrTy & key) {auto ii=m_map.find(key);
if (ii==m_map.end()) return m_blank;
return (*ii).second ; } 
template <class Tm> void add(Tm &m ) 
{ MM_LOOP(ii,m) { m_map[(*ii).first]=(*ii).second; } }
template <class Tm> void count(Tm &m ) 
{ MM_LOOP(ii,m_map) { ++m[(*ii).first]; } }
void make_name(TempProg & tp) { MakeName(tp); } 
StrTy format( ) { Ss ss; Format(ss); return ss.str(); }
void save(OsTy & ss,const IdxTy flags) { Save(ss,flags); } 
StrTy dump(const IdxTy flags ) { Ss ss; Dump(ss,flags); return ss.str(); }
bool missing(const StrTy & f) { return m_map.find(f)==m_map.end(); } 
IdxTy make_keys_uniform(const IdxTy flags=0) { return MakeKeysUniform(flags); } 
private:
StrTy Canonical(const StrTy & s)
{
const IdxTy slen=s.length();
const char *p =s.c_str();
char dc[slen+1];
for(IdxTy i=0; i<slen; ++i)
{
char c=p[i];
//int  v=m_clut.lut(c);
//if ((v&(ParseTable::UC))!=0)  c=char(int(c)+(int('a')-int('A')));
// FIXME need to come up with translate version of lut lol
// Want to do lut based parser too... 
if ((c>='A')&&(c<='Z'))  c=char(int(c)+(int('a')-int('A')));
dc[i]=c;
}
dc[slen]=0;
return StrTy(dc); 
}
StrTy CanonicalValue(const StrTy & s)
{
StrTy sc=ParseTable::remove_mask(s,ParseTable::EOL);
sc=ParseTable::remove_mask(sc,ParseTable::WHITE,1);
sc=ParseTable::remove_mask(sc,ParseTable::WHITE,2);
return sc;
}

IdxTy  MakeKeysUniform(const IdxTy flags) 
{
IdxTy i=0;
Map m;
MM_LOOP(ii,m_map) 
{
const StrTy &  k=(*ii).first;
const StrTy & v=(*ii).second;
const StrTy  vc=CanonicalValue(v);
const StrTy d=Canonical(k);
//if (k!=d) { ++i; }
if (m.find(d)!=m.end()) { ++i;  MM_ERR(" keys overlap "<<MMPR4(k,d,v,m[d])) } 
m[d]=vc;

} // m_map
m_map=m;
MM_LOOP(ii,m_order_vector) { (*ii)=Canonical((*ii)); }
return i;
}


void MakeName(TempProg & tp) { 
V t,a,b; // u;
//MM_LOOP(ii,m_map) { MM_MSG(MMPR2((*ii).first,(*ii).second)) }
StrTy yr=m_map["year"];
m_gen.make_words(t,m_map["title"]);
m_gen.make_words(a,m_map["authors"]);
m_gen.make_words(b,m_map["author"]);
//m_gen.make_words(u,m_map["url"]);
m_gen.append_interesting(m_name,t);
m_gen.append_interesting(m_name,a);
m_gen.append_interesting(m_name,b);
m_gen.append_interesting_url(m_name,m_map["url"]);

m_name=m_name+yr;

} 

void Dump(OsTy & ss,const IdxTy flags)
{
const StrTy sep=" ";
const StrTy lead=m_name+sep+m_type;
MM_LOOP(ii,m_map) { ss<<lead<<sep<<(*ii).first<<sep<<(*ii).second<<CRLF; } 

}
void Save(OsTy & ss,const IdxTy flags) { 
//Save(ss,flags); 
const bool use_map_order=(Bit(flags,0));
const bool include_comments=(Bit(flags,1));
//MM_ERR("  "<<MMPR2(flags,use_map_order))
const StrTy spaces="    ";
const IdxTy osz=m_order_vector.size();
const IdxTy kvsz=m_map.size();
if (osz!=0) if (osz!=kvsz) { 
ss<<"% FIXME autogenerated  map and order differ "<<MMPR2(osz,kvsz)<<CRLF ; 
MM_ERR(" size error "<<MMPR2(m_order_vector.size(),m_map.size())) 
IdxTy vi=0; 
MM_LOOP(ii,m_order_vector) { if (m_map.find(*ii)==m_map.end())
{ ss<<"% map is missing "<<(*ii)<<" "<<vi<<CRLF;  } 
++vi;
}
} // size error 
// this is up to the caller
//ss<<CRLF;
if (include_comments) 
{ MM_LOOP(ii,m_comments) { ss<< ParseTable::remove_run((*ii).second,ParseTable::EOL) <<CRLF; } }
ss<<"@"<<m_type<<"{"<<m_name; // <<","<<CRLF;
const IdxTy sz=m_map.size();
const bool use_order=(osz==kvsz)&&!use_map_order;
if (use_order)
{
MM_LOOP(ii,m_order_vector) 
{ ss<<","<<CRLF<<spaces; Format(ss,*ii,m_map[*ii]); }
} // use_order
else
{
MM_LOOP(ii,m_map) 
{ss<<","<<CRLF<<spaces;  Format(ss,(*ii).first,(*ii).second); }

} // use_order 
if (kvsz>0) ss<<CRLF;
ss<<"}"<<CRLF;

} // save

void Format(OsTy & ss)
{
MM_LOOP(ii,m_comments) { ss<<"%" << (*ii).second<<CRLF; } 
ss<<"@"<<m_type<<"{"<<m_name; // <<CRLF;
const IdxTy sz=m_map.size();
IdxTy i=0; 
if ( i!=sz) ss<<","; ss<<CRLF;  
MM_LOOP(ii,m_map) 
{
	++i; Format(ss,(*ii).first,(*ii).second);
	if ( i!=sz) ss<<","; ss<<CRLF;  
}
ss<<"}"<<CRLF;
}

void Format(OsTy & ss, const StrTy & key, const StrTy & v)
{ ss<<key<<" = {"<<v<<"}"; }

void Init()
{
m_type="";
m_name="";
m_map.clear();
m_order.clear();
m_comments.clear();
m_serial=bad();
}
IdxTy  Parse(const char * p, const IdxTy len)
{
typedef mjm_bibtex_parse<Tr> Bp;
Bp bib_parse;

return bib_parse.parse(m_map,p,len,0);
}

StrTy m_type, m_name,m_lexi_added_date,m_source_file;
Map m_map;
Ord m_order;
V m_order_vector;
Map m_comments;
Ngen m_gen;
IdxTy m_serial;
StrTy m_blank;
}; // mjm_bibtex_entry


template <class Tr> 
class mjm_bibtex_fixer 
{

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;

typedef mjm_bibtex_entry<Tr> Tgt;
typedef std::vector<Tgt> TgtV;
typedef std::map<StrTy, TgtV> TgtM;
typedef mjm_blob<Tr> Blob;
typedef mjm_read_buffer<Tr> RdBuf;
typedef mjm_var_parse<Tr> ParseTable;
typedef CommandInterpretter Ci;
public:
mjm_bibtex_fixer() { Init(); } 

private:

void Init() {}

}; // mjm_bibtex_fixer


template <class Tr> 
class mjm_bibtex_entry_map 
{

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;

typedef mjm_bibtex_entry<Tr> Tgt;
typedef std::vector<Tgt> TgtV;
typedef std::map<StrTy, TgtV> TgtM;
typedef mjm_blob<Tr> Blob;
typedef mjm_read_buffer<Tr> RdBuf;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_pawnoff<Tr>  Handler;
typedef CommandInterpretter Ci;
public:
mjm_bibtex_entry_map() { Init(); } 

// these rely on have no values with zero length vectors.. 
class iter
{
public:
iter(TgtM & m ) : m_map(m) {
m_ii=m_map.begin();
m_i=0;
}
Tgt & operator*() { return(*m_ii).second[m_i]; }
iter & operator++() { ++m_i; 
if (m_i>=(*m_ii).second.size()){  ++m_ii; m_i=0; }
return *this;
}
 operator bool () {  return m_ii!=m_map.end(); } 
TgtM & m_map;
IdxTy m_i;
typename TgtM::iterator m_ii;

};

typedef std::vector<StrTy> OrdVec;
class iter_ord
{
public:
iter_ord(TgtM & m, const OrdVec & o  ) : m_map(m), m_ord(o) {
m_sz=m_ord.size();
m_o=0;
if ( m_sz>0) m_ii=m_map.find(m_ord[0]);
m_i=0;
}
Tgt & operator*() { return(*m_ii).second[m_i]; }
iter_ord & operator++() { ++m_i; 
if (m_i>=(*m_ii).second.size()){ ++m_o ;  
if (m_o<m_sz)  m_ii=m_map.find(m_ord[m_o]);
 m_i=0; }
return *this;
}
 operator bool () {  return (m_o<m_sz); } 
TgtM & m_map;
const OrdVec  & m_ord;
IdxTy m_i,m_o,m_sz;
typename TgtM::iterator m_ii;

};


void add(const Tgt & t) { Add(t); } 
void parse( const StrTy nm)  { Parse(nm); } 
void dump(OsTy & os, const IdxTy flags) { Dump(os,flags); } 
void save(const StrTy  & fn, const IdxTy flags) { Save(fn,flags); } 
void command_mode(Ci & li) { CommandMode( li) ; } 
IdxTy cmd(const StrTy & p1, const StrTy & p2 )  { return Cmd(p1,p2); } 
void try_to_fix(Ci & li) {  TryToFix( li); } 
void check_comments(OsTy & os, const IdxTy flags) 
 {  CheckComments( os,  flags) ;} 
void missing(OsTy & os, const StrTy & f, const IdxTy flags) { Missing(os,f,flags); } 
private:
bool Bit(const IdxTy f, const IdxTy b) { return ((1<<b)&f)!=0; }
bool Mask(const IdxTy f, const IdxTy m) { return ((m)&f)!=0; }
IdxTy Cmd(const StrTy & p1, const StrTy & p2 )  { 

if (p1=="isc") IntegrateStructuredComments();


return 0; 

} 
void CommandMode(Ci & li)
{
while (li.nextok())
{
const IdxTy sz=li.size();
CommandInterpretterParam  cip(li);

if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd.c_str()[0]=='#' ) continue;

const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
//if (cmd=="about") { about();  continue; }
//if (cmd=="print") { MM_MSG(li.line())  continue; }

} // nextok

} //CommandMode


void Missing(OsTy & os, const StrTy & f, const IdxTy flags) 
{ 
 // Missing(os,f,flags); 
iter x(m_map);
while (x) { if ((*x).missing(f)) { os<<(*x).type()<<" "<<(*x).name()<<CRLF;}  ++x; } // iter

} 

void Save(const StrTy  & fn, const IdxTy flags) 
{ 
std::ofstream os(fn.c_str());
Save(os,flags); 
} 

void Save(OsTy & os, const IdxTy flags) 
{
Handler h;
StrTy today=ParseTable::remove_mask(h.today(),ParseTable::WHITE);
os<<"% programmatically fixed probably bu toobib"<<CRLF;
os<<"% loaded from "<<m_loaded_from<<" written on "<<today<<CRLF;
IdxTy serial=0;
iter_ord x(m_map,m_ord);
while (x) {  
os<<"%"<< serial<<" prior "<<(*x).serial()<<CRLF;
(*x).save(os,flags);  ++x; 
++serial; 
} // iter

}

void Dump(OsTy & os, const IdxTy flags) 
{ 
//Dump(os,flags); 
const bool all_entries=(!Bit(flags,0));
const bool summary=(Bit(flags,1));
const bool count=(Bit(flags,2));
const bool lines=(Bit(flags,3));
if ( all_entries) {
iter x(m_map);
while (x) { os<<(*x).format(); ++x; } // iter
}
if (summary)
{
os<<MMPR(m_map.size())<<CRLF;
MM_LOOP(ii,m_map) 
{ os<<"name="<<(*ii).first<<" size="<<(*ii).second.size()<<CRLF; }
}
if (count)
{
typedef std::map<StrTy, IdxTy> Cm;
Cm m;
IdxTy n=0;
iter x(m_map);
while (x) { (*x).count(m); ++m[StrTy("@")+(*x).type()]; ++x; }
MM_LOOP(ii,m) { os<<" field="<<(*ii).first<<" "<<(*ii).second<<CRLF; ++n; } 
} // count
if (lines)
{
IdxTy n=0;
iter x(m_map);
while (x) { os<<(*x).dump(0); ++x; }
//MM_LOOP(ii,m) { os<<" field="<<(*ii).first<<" "<<(*ii).second<<CRLF; } 
} // count



} // Dump

 
void Parse( const StrTy nm)
{
typedef std::map<StrTy,StrTy> M;
typedef mjm_bibtex_parse<Tr> Bp;
Bp bib_parse;
Blob x;
m_loaded_from=nm;
x.load(nm);
IdxTy len=x.size();
IdxTy pos=0;
while (pos<len)
{
M m;
Tgt t;
IdxTy dx= bib_parse.parse(m,x.ptr()+pos,len-pos,0);
pos+=dx;
// may get spurs at end or elsewhere
const bool has_type=(m.end()!=m.find("type_name"));
const bool has_name=(m.end()!=m.find("name_name"));
const IdxTy msz=m.size();
// the map size can't be zero if it has a name or type
// just for thought 
if (msz==0) continue;
//if (!has_name) continue;
// remove name, type, and comments then add to map
t.name(m["name_name"]);
t.type(m["type_name"]);
t.comment(m["junk_name"]);
m.erase(m.find("name_name"));
m.erase(m.find("type_name"));
m.erase(m.find("junk_name"));
t.serial(m_entries);
t.add(m);
std::vector<StrTy> v;
bib_parse.get_last_keys(v);
t.set_order_vector(v);
//IdxTy dx=t.parse(m,x.ptr()+pos,len-pos,0);
//os<<MMPR3(pos,dx,len);
//MM_LOOP(ii,m) { os<<MMPR2((*ii).first, (*ii).second)<<CRLF; }
add(t);

}
} // Parse

void TryToFix(Ci & li)
{

iter x(m_map);
while (x) {IdxTy rc=TryToFix(*x,li);if (rc==1) break;  ++x; } // iter
// Duplicates
MM_LOOP(ii,m_map)
{ if ((*ii).second.size()!=1) ResolveDuplicates(li,(*ii).first); }

}
void ResolveDuplicates(Ci &li, const StrTy &nm)
{
TgtV orig=m_map[nm];
MM_LOOP(ii,orig)
{ Tgt & x=(*ii); 
const StrTy & authors= x.get("authors");
MM_MSG(MMPR4(x.serial(),x.name(),x.type(),authors)) }


} // ResolveDuplicates


StrTy Word(const char *p, const IdxTy len,  const IdxTy  first)
{ return ParseTable::word(p,len,first); } 
#if 0
StrTy Word(const char *p, const IdxTy len, const IdxTy first)
{
IdxTy ifirst=~0;
IdxTy ilast=~0;
IdxTy state=0;
for(IdxTy i=first; i<len; ++i)
{
const char c=p[i];
if (c==0) { ilast=i; break; }
typename ParseTable::Iv v=m_clut.lut(c);
if (Mask(v,ParseTable::WHITE)) { 
if (state==0) continue;
if (state==1) { ilast=i; break ; } 
} //mask
else
{ if (state==0){  ifirst=i; state=1; } }
}// for 
if (ifirst==~0) return StrTy(); 
if (ilast==~0) ilast=len; 

const IdxTy sz=ilast-ifirst;
char ca[sz+1];
memcpy(ca,p+ifirst,sz);
ca[sz]=0;
return StrTy(ca);
} // Word
#endif

//typedef mjm_ragged_table Ragged;
IdxTy IntegrateStructuredComments()
{
MM_MSG(" doing Integrate Structured Comments")
IdxTy olaps=0;
iter x(m_map);
while (x) {
const IdxTy olap=(*x).make_keys_uniform();
olaps+=olap;
const IdxTy rc=IntegrateStructuredComments(*x);if (rc!=0) return rc;  ++x; } 
if (olaps!=0) { MM_ERR(" some keys overalpped "<<MMPR(olaps)) } 
return 0; 
}

//void  CommentsToRags( gg
IdxTy IntegrateStructuredComments(Tgt & bib)
{
const bool missing_url=bib.missing("url");
const bool missing_srcurl=bib.missing("srcurl");
const bool missing_citeurl=bib.missing("citeurl");
const bool missing_medtobib=bib.missing("medtobib");
const bool missing_mjmdate=bib.missing("mjmdate");
const StrTy comment=bib.comment();
// a string is taken as a file name
//Ci cparse(comment);
Ss ss; ss<<comment;
Ss snew;
Ci cparse(&ss);
IdxTy mods=0; 
//Ragged ragged;
//ragged.load_lines(ss);
while (cparse.nextok())
{
const StrTy & l=cparse.line();
const char * p=l.c_str();
const IdxTy len=l.length();
//ParseTable::word(
StrTy val=Kvals(p,len,"srcurl");
if (val.length()!=0) {++mods;  bib.seto("srcurl",val); if (missing_url) bib.seto("url",val); } 
else
{
StrTy val=Kvals(p,len,"citeurl");
if (val.length()!=0) {++mods;  bib.seto("citeurl",val);} 
else{
StrTy val=Kvals(p,len,"med2bib",true);
if (val.length()!=0) { ++mods;  bib.seto("medtobib",val);} 
else{
StrTy val=Kvals(p,len,"date",true);
if (val.length()!=0) { 
Handler h;
//MM_ERR(" trying to get lexi date for "<<val)
val=h.dates(val);
val=ParseTable::remove_mask(val,ParseTable::WHITE);
++mods;
 bib.seto("mjmdate",val);

} 
else snew<<l<<CRLF;
}
}
} 
} // while 
if (mods!=0) bib.comment(snew.str());
return 0; 
}

StrTy Kvals(const char * p, const IdxTy len, const StrTy & key,bool rol=false)
{
const IdxTy pos=Blob::find_blob(p,len,0,key.c_str());
if (pos==~0) return StrTy();
StrTy val= Word(p, len, pos);
if ( !rol) { 
return Word(p,len,pos+val.length()); 
}
val=StrTy(p+pos+val.length());
return val; 
}


IdxTy TryToFix(Tgt & bib, Ci & li)
{
// get the comments, fix, ignore, mark FIXME 
// if there is no url, see if src url exists in comments
const bool missing_url=bib.missing("url");
const bool update_comment=true;
const StrTy comment=bib.comment();
const char * pp=comment.c_str();
const IdxTy slen=comment.length();
const IdxTy pos=Blob::find_blob(pp,slen,0,"srcurl");
const StrTy link=(pos!=~0) ?Word(pp,slen,pos+7):StrTy("");

if (missing_url)
{ Ss ss;
ss<<" url is missing ";
//if (link.length()!=0) 
ss<<bib.serial()<<" "<<bib.name()<<" "<<bib.get("author")<< "  srcurl is "<<link;
//else ss<<" no srcurl "; 
MM_MSG(ss.str())
MM_MSG(" ignore fix mark comment title quit ")
while (li.nextok()){
const IdxTy sz=li.size();
//CommandInterpretterParam  cip(li);
//if (sz<1) continue;
const StrTy cmd=(sz<1)?StrTy():li.word(0);
// comment should be updated when fixed 
if (cmd=="f") {  bib.set("url",link); bib.push_order_vector("url");  break ; } 
if (cmd=="m") { bib.set("url","FIXME");  bib.push_order_vector("url");  break; } 
if (cmd=="c") {  MM_MSG(comment)  } 
if (cmd=="t") {  MM_MSG(bib.get("title") )  } 
if (cmd=="a") {  MM_MSG(bib.dump(0) )  } 
if (cmd=="q")  return 1;  
if (cmd=="i")  break;  
MM_MSG(" one of f m c t q i")
} // nextok

} // missing_url 
// fix comment  
if (update_comment)
{
char nc[3*slen+2];
IdxTy  pc=0;
IdxTy mods=0;
IdxTy state=0;
for(IdxTy i=0; i<slen; ++i)
{
const char c=pp[i];
typename ParseTable::Iv v=m_clut.lut(c);
// if at beginning of line, ignore white put in percent 
switch (state)
{
case 0 :  { 
if (Mask(v,ParseTable::WHITE)) {break; }  
if (Mask(v,ParseTable::EOL)) {break; }  
if ( c!='%') { nc[pc]='%'; ++pc; ++mods; }
 nc[pc]=c; ++pc;  state=1;
 break;
}
case 1: { nc[pc]=c; ++pc;
if (Mask(v,ParseTable::EOL)) {state=0; break; }  
break;
}


} // switch 

} // i 
if (mods!=0)
{
nc[pc]=0;
const StrTy ncom=StrTy(nc);
MM_MSG( " old comment "<<comment<<CRLF);
MM_MSG( " new comment "<<nc<<CRLF);
MM_MSG( " accept reject  "<<CRLF);
{(li.nextok());
const IdxTy sz=li.size();
const StrTy cmd=(sz<1)?StrTy():li.word(0);
// comment should be updated when fixed 
if (cmd=="a") {  bib.comment(nc);   } 
}
} // mods 


} // update_comment




return 0; 
}

enum { OKB4=ParseTable::WHITE, BEOL=ParseTable::EOL, ASCII=ParseTable::ASCII };

void CheckComment(StrTy & fixed, StrTy & input, Ci& li,  const IdxTy flags) 
{



}
void CheckComments(OsTy & os, const IdxTy flags) 
{ 
//Dump(os,flags); 
const bool all_entries=(!Bit(flags,0));
const bool summary=(Bit(flags,1));
const bool count=(Bit(flags,2));
const bool lines=(Bit(flags,3));
if ( all_entries) {
iter x(m_map);
//while (x) { os<<(*x).format(); ++x; } // iter
while (x) {StrTy p=(*x).comment();  CheckComment((*x).name(),p.c_str(),p.length());  ++x; } // iter
}
}

IdxTy CheckComment(const StrTy & lbl,const char * p, const IdxTy len, const IdxTy flags=0)
{
RdBuf b(10);
b.start_new();
bool comment=false;
bool line_bad=false;

for(IdxTy i=0; i<len; ++i)
{
bool flagged=false;
const char c=p[i];
typename ParseTable::Iv v=m_clut.lut(c);
if (Mask(v,BEOL)) {b.start_new(); 
if (line_bad)
{
MM_ERR(" bad line "<<MMPR2(lbl,b.next_last_string()))

}
line_bad=false;
 comment=false; continue;  } 
if (!Mask(v,ASCII)) { line_bad=true; flagged=true; MM_ERR(" non ASCII in comment "<<std::hex<<MMPR3(lbl,IdxTy(c),c))  } 
if (c=='%') comment=true;
if (!comment) if (!Mask(v,OKB4)) 
{
flagged=true;
line_bad=true;
}
if (flagged || !flagged) b.append(c);

}
b.cap();
if (line_bad)
{
MM_ERR(" bad line "<<MMPR2(lbl,b.next_last_string()))
}
return 0;
}

void Init()
{
m_entries=0;
}


void Add(const Tgt & t) { 
// also can use lexi date when added 
m_map[t.name()].push_back(t);
m_ord.push_back(t.name()); 
++m_entries;
 } 



TgtM m_map;
IdxTy m_entries;
OrdVec m_ord;
ParseTable m_clut;
StrTy m_loaded_from;

}; // mjm_bibtex_entry_map



#endif
