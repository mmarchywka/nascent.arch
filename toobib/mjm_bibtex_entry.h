
#ifndef MJM_bibtex_entry_H__
#define MJM_bibtex_entry_H__
 
#include "mjm_globals.h"

#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_bibtex_fields.h"

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
#include "mjm_wovdb.h"

#include <algorithm>
#include <map>
#include <ctype.h>
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
//typedef mjm_block_matrix<D> MyBlock;
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
{ if (n.length()) n=n+"_"; n=n+v[i]; ++j; } ++i; if (i>3) break; } ;
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
if (s[i]=="/") { IdxTy l=i-1; if ( l<=w.size())if(w[l].length()>4)  n=n+"_"+w[l]; break; }
 }
if ( i==sz){ IdxTy l=i-1; if ( l<=w.size()) if(w[l].length()>4)  n=n+"_"+w[l]; }
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



typedef std::map<StrTy, StrTy> MapTy;
typedef std::vector<StrTy> LiTy;
// this is the idiotic map for the key value paris
// instead of DB etc.. 
typedef std::map<StrTy, LiTy> MmapTy;
typedef std::map<StrTy, StrTy> Ord;
typedef std::vector<StrTy> V;
typedef mjm_bib_name_generator<Tr> Ngen;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_read_buffer<Tr> RdBuf;
typedef temp_progress<Tr>  TempProg;
typedef mjm_wovdb<Tr,StrTy> WovDB;

typedef mjm_bibtex_fields<Tr>  Fields;

typedef typename WovDB::vector_type Vt;
typedef CommandInterpretter Ci;
static bool Bit(const IdxTy f, const IdxTy b) { return ((1<<b)&f)!=0; }
static bool Mask(const IdxTy f, const IdxTy m) { return ((m)&f)!=0; }


public:
typedef MapTy kvp_type;
mjm_bibtex_entry() {Init(); } 
IdxTy bad() const  { return ~0; } 
StrTy operator[](const StrTy & nm) const
{
auto ii=m_map.find(nm);
if (ii==m_map.end()) return StrTy();
auto & jj=(*ii).second;
if (jj.size()==0)  return StrTy();
return jj[0];
}

// this was never changed to a dv, just a str to vec map doh 
const MmapTy & map() const { return m_map; }
void clear() { Init(); } 
void serial(const IdxTy & x) { m_serial=x;}
IdxTy serial()const  {return  m_serial;}
IdxTy size() const { return m_map.size(); } 
IdxTy errors() const { return m_errors.size(); } 
IdxTy observations() const { return m_observations.size(); } 
const StrTy source() const { return m_source_file; } 
void source(const StrTy & s )  { m_source_file=s; } 
void type(const StrTy & x) { m_type=x;}
void type(const Vt *  x,WovDB &m) {if (x)if((*x).size())  m_type=m[(*x)[0]];}
void clean_name() { CleanName();}
void name(const StrTy & x) { m_name=x; CleanName();}
void verbatim_name(const StrTy & x) { m_name=x; }
void name(const Vt *  x, WovDB & m ) {if (x)if((*x).size())  m_name=m[(*x)[0]];}
StrTy missing_type() const { return "missing";}
StrTy missing_name() const { return "NONE";}
const bool has_good_name() const { return (m_name!="")&&(m_name!=missing_name()); } 
const bool has_good_type() const { return (m_type!="")&&(m_type!=missing_type()); } 
const StrTy & type() const { return m_type;}
const StrTy & name() const { return m_name;}
const StrTy & text() const { return m_original_text;} 
void text(const char * p, const IdxTy s, const IdxTy sz)
{ m_original_text=StrTy(p+s, sz); m_buf_start=s; m_buf_size=sz; }

void comment( const StrTy & key, const StrTy & v) { m_comments[key]=v; } 
// repeated calls now overwrite
void comment(  const StrTy & v) 
//{Ss ss; ss<<m_comments.size();  m_comments[ss.str()]=v; } 
{Ss ss;   m_comments[ss.str()]=v; } 

void comment(const Vt *  x, WovDB &m) 
{if (x)if((*x).size())  m_comments[""]=m[(*x)[0]];}

StrTy comment()
{
if (m_comments.size()==0) return StrTy();
return (*(m_comments.begin())).second;

}
void set_order_vector(const V & order_vector) { m_order_vector=order_vector; } 
void push_order_vector(const StrTy & k) { m_order_vector.push_back(k); } 
//void set( const StrTy & key, const StrTy & v) { m_map[key]=v; } 
void set( const StrTy & key, const StrTy & v) { m_map[key].push_back(v); } 
void seto( const StrTy & key, const StrTy & v) {
//if (missing(key)) push_order_vector(key);set(key,v); } 
push_order_vector(key);set(key,v); } 
const StrTy & get( const StrTy & key,const IdxTy i=0) {
auto ii=m_map.find(key);
if (ii==m_map.end()) return m_blank;
const auto & v=(*ii).second;
if (v.size()<=i ) return m_blank;
//return (*ii).second ; 
return v[i] ; 
} 
// the get( sig is used in nor-const  for str return  
//const LiTy & get(const StrTy & key) const { 
const LiTy & get_list(const StrTy & key) const { 
const auto& ii=m_map.find(key);
if (ii!=m_map.end()) return (*ii).second;
static LiTy x;
return x; 
} // get

const IdxTy count_values( const StrTy & key) {
auto ii=m_map.find(key);
if (ii==m_map.end()) return 0;
const auto & v=(*ii).second;
return v.size();
} 
IdxTy non_x(const std::map<StrTy, IdxTy> & rejects) const 
{
IdxTy n=0;
MM_LOOP(ii,m_map) {
const StrTy & x=(*ii).first;
const char *  p=x.c_str();
bool no_sp=false;
for(IdxTy i=0; i<x.length(); ++i) if (p[i]!=' '){ no_sp=true; break; }
if (!no_sp )continue ;
if (( p[0]==0)) continue;
if (( p[0]=='x')||(p[0]=='X')) if (p[1]=='_') continue;
// X-UA-Compatible access pagetitle publisher robots theme-color title viewport
if ( rejects.find(x)!=rejects.end()) continue;
const auto & v=(*ii).second;
if (v.size()==0) continue;
IdxTy sz=0;
MM_LOOP(jj,v) sz+=(*jj).length();
if (sz==0) continue;
++n;
} // ii 
return n;
} // non_x
StrTy key_str() const
{
Ss ss;
MM_LOOP(ii,m_map) {
const StrTy & x=(*ii).first;
ss<<" "<<x;
} // ii 


return ss.str();
} // key_str

const StrTy & blank() const { return m_blank; } 
//template <class Tm> void add(Tm &m ) 
//{ MM_LOOP(ii,m) { m_map[(*ii).first]=(*ii).second; } }
void clear_map() { m_map.clear(); } 
void add(const MapTy &m ) 
{ MM_LOOP(ii,m) { m_map[(*ii).first].push_back((*ii).second); } }
void add(const StrTy & k, const StrTy & v ) 
{  { m_map[k].push_back(v); } }
// set was alredy used for push.. 
void only(const StrTy & k, const StrTy & v ) 
{  { m_map[k].clear(); m_map[k].push_back(v); } }


void add(MmapTy &m ) 
{ MM_LOOP(ii,m) { m_map[(*ii).first]=((*ii).second); } }
void replace(const StrTy& key, const StrTy & v, const IdxTy pos) 
{
auto ii=m_map.find(key);
if (ii==m_map.end()) { m_map[key].push_back(v); return ; }
auto& ve=(*ii).second;
if (ve.size()==0) { ve.push_back(v); return ; }
ve[pos]=v;

} // replace

IdxTy remove_long( const StrTy & key, const IdxTy len,const IdxTy flags=0) 
{
IdxTy rc=0;
const bool note_dropped=!Bit(flags,0);
const bool first_few=!Bit(flags,1);
auto ii=m_map.find(key);
//if (ii==m_map.end()) { m_map[key].push_back(v); return ; }
auto& ve=(*ii).second;
auto vnew=ve;
vnew.clear();
MM_LOOP(jj,ve) { if ((*jj).length()<len) vnew.push_back((*jj));
else
{
++rc;
Ss ss;
if (note_dropped) 
{
ss<<" deleted for space";
} // note_dropped
if (first_few) 
{
Ss ss;
const IdxTy nfirst=20;
const IdxTy sz=(len<nfirst)?len:nfirst;
char * p= new char[sz+1];
memcpy(p,(*jj).c_str(),sz);
p[sz]=0;
ss<<StrTy(p);
} // note_dropped


if (ss.str().length()) vnew.push_back(ss.str());


} // else

 }  // jj 
//if (ve.size()==0) { ve.push_back(v); return ; }
//ve[pos]=v;
m_map[key]=vnew;
return rc;
} // remove_long



IdxTy remove(const StrTy& key,const IdxTy flags) 
{
auto ii=m_map.find(key);
if (ii==m_map.end()) { return 1  ; }
m_map.erase(ii);
return 0; 
} //remove 
void add_db(WovDB &m ) 
{ 
auto ii=m.begin("");
//if (ii) MM_ERR(MMPR(ii.key()))
//else MM_ERR(" ii at end ")
while (ii) {
//MM_ERR(MMPR2(bool(ii),ii.key()))
 m_map[ii.key()].push_back((*ii)); ++ii; } 
//MM_ERR(" exit add_db")
add_something("error",m,m_errors);
add_something("observe",m,m_observations);
/*
{ auto ii=m.begin("error");
while (ii) {
//MM_ERR(MMPR2(bool(ii),ii.key()))
 m_errors[ii.key()].push_back((*ii)); ++ii; } 
}
*/
} // add_db
// TODO wov const iter 
void add_something(const char * nm, /*const */  WovDB & m,MmapTy & map) const
{
 auto ii=m.begin(nm);
while (ii) {
 map[ii.key()].push_back((*ii)); ++ii; } 

} // add_something 
template <class Tm> void count(Tm &m ) 
//{ MM_LOOP(ii,m_map) { ++m[(*ii).first]; } }
{ MM_LOOP(ii,m) { ++m[(*ii).first]; } }

void make_name() {TempProg tp;  MakeName(tp); } 
void make_name(TempProg & tp) { MakeName(tp); } 
StrTy format()const  { Ss ss; Format(ss); return ss.str(); }
template <class Trag > void  format_rag(Trag  & r,const IdxTy flags) const 
{ return FormatRag( r,flags );} 
void save(OsTy & ss,const IdxTy flags)const { Save(ss,flags); } 
IdxTy  save(const StrTy  & fn,const IdxTy flags)const 
{ std::ofstream of(fn); Save(of,flags); return !of.good();  } 
StrTy dump(const IdxTy flags ) { Ss ss; Dump(ss,flags); return ss.str(); }
StrTy dump_errors(const IdxTy flags ) const  { Ss ss; DumpErrors(ss,flags); return ss.str(); }
bool missing(const StrTy & f) { return m_map.find(f)==m_map.end(); } 
//const StrTy d=Canonical(k);
bool missing_canonical(const StrTy & f) { return missing(Canonical(f)); } 
// TODO this may not work right. 
IdxTy make_keys_uniform(const IdxTy flags=0) { return MakeKeysUniform(flags); } 

IdxTy fix_keys()
{
MM_ERR(MMPR(__FUNCTION__))
const bool annotate_changes=!true;
ParseTable clut;
IdxTy nfix=0;
auto assfck=m_map;
// changes modify the iterator doh 
MM_LOOP(ii,assfck)
{
const StrTy & kold=(*ii).first;
MM_ERR(MMPR2(kold,(*ii).second.size()))
const char * p=kold.c_str();
IdxTy i=0;
bool mods=false;
IdxTy d=0;
IdxTy ed=0;
char cc[kold.length()+1];
char ecc[3*kold.length()+3];
while (p[i]!=0)
{
const char c=p[i];
bool bad=(c=='*')||(c=='$')||(c=='@')|| (c=='\\')||(c=='\'')||(c=='"')||(c=='{')||(c=='}') ||(c==',')|| (c=='-');
const IdxTy iv=clut.lut(c);
bad|=((iv&(ParseTable::DARK))==0) ;
bad|= ( c&128) ;
bool ok=!bad;
if ( ok ) { cc[d]=p[i]; ++d; ecc[ed]=c; ++ed; } else
{
 mods=true;
++nfix;
ecc[ed]='\\'; ++ed;
ecc[ed]=c; ++ed;
}
++i;
}
if (i==0)
{
mods=true;
memcpy(cc,"X_BLANK_KEY",12);
d=strlen(ecc);
}
cc[d]=0;
ecc[ed]=0;
if (mods)
{
MM_ERR(" fixing key "<<MMPR4(name(),StrTy(cc),kold,nfix))
if ( annotate_changes){
// this risks a problem in the value field though... doh 
Ss ss; ss<<MMPR3(StrTy(cc),StrTy(ecc),nfix);
 add(StrTy("x_fix_key"),ss.str()); 
}
 move(StrTy(cc),kold,2);

}
} // ii 
return nfix;
} // fix_keys

// renmae is common from various sources oftern for case or plural 
IdxTy move(const StrTy & knew, const StrTy & kold,const IdxTy flags)
{
const bool only_exists=Bit(flags,0);
const bool append=Bit(flags,1);
const bool actually_copy=Bit(flags,2);
// any existing one has to go anyway.
if (!only_exists) { 
auto ii=m_map.find(knew);
if ( ii!=m_map.end()) m_map.erase(ii);
}
auto jj=m_map.find(kold);
if (jj!=m_map.end())
{
if (!append)  m_map[knew]=(*jj).second;
else {
MM_LOOP(kk,(*jj).second) m_map[knew].push_back(*kk); 
}
if (!actually_copy) { 
auto kk=m_map.find(kold);
m_map.erase(kk); 
 } // actually_copy

}// jj 
return 0;
} // move
// make the keys all uniform... BUT ALSO REPLACE...  
IdxTy make_keys_canonical(const IdxTy flags)
{
//typedef std::map<StrTy, LiTy> MmapTy;
MmapTy newmap;
// unfortunately this should be backed with order
// preserving WovDB and then eliminated... 
MM_LOOP(ii,m_map)
{
const StrTy& k=(*ii).first;
const auto & v=(*ii).second;
const StrTy d=Canonical(k);
auto & vn=newmap[d];
MM_LOOP(jj,v) { vn.push_back(*jj); } 

}
m_map=newmap;
return 0; 
} // make_keys_canonical

IdxTy trim_values(const IdxTy flags)
{
//typedef std::map<StrTy, LiTy> MmapTy;
MmapTy newmap;
// unfortunately this should be backed with order
// preserving WovDB and then eliminated... 
MM_LOOP(ii,m_map)
{
const StrTy& k=(*ii).first;
const auto & v=(*ii).second;
const StrTy d=k; // Canonical(k);
auto & vn=newmap[d];
MM_LOOP(jj,v) { vn.push_back(mjm_strings::trim(*jj)); } 

}
m_map=newmap;
return 0; 
} // make_keys_canonical



typedef mjm_strings StrUtil;
template<class Tkm>
IdxTy matching_keys(Tkm & jm, const StrTy & s, const IdxTy flags)
{

MM_LOOP(ii,m_map)
{
const StrTy& k=(*ii).first;
const int  pos=StrUtil::indexOf(k.c_str(),s.c_str());
if (pos!=(-1)) jm[k]=pos;
} // ii 
return -1;
} // matching_keys
template<class Tkm>
IdxTy all_keys(Tkm & v,  const IdxTy flags=0)
{

MM_LOOP(ii,m_map)
{
const StrTy& k=(*ii).first;
v.push_back(k);
} // ii 
return 0;
} // matching_keys



private:
void CleanName()
{
// need a char-class thing but this is ok for now...
const IdxTy len=m_name.length();
const char * p=m_name.c_str();
// this is always a risk with data dependent sizes unchecked...
char cp[len+1];
IdxTy j=0;
for(IdxTy i=0; i<len; ++i)
{
const char c=p[i];
if ( c&128) continue;
if (!isalnum(c)&&(c!='_')&&(c!='/')) continue;
cp[j]=c; ++j;
} //i 
cp[j]=0; // doh WTF 2021-06-02
m_name=StrTy(cp);

}  // CleanName



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
// mjm_strings::fancy_to_lower?
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
//Map m;
MmapTy m;
MM_LOOP(ii,m_map) 
{
const StrTy &  k=(*ii).first;
//const StrTy & v=(*ii).second;
const IdxTy szv=(*ii).second.size(); // all keys must have one value 
const StrTy & v=(*ii).second[0]; // all keys must have one value 
const StrTy d=Canonical(k);
//if (k!=d) { ++i; }
//if (m.find(d)!=m.end()) { ++i;  MM_ERR(" keys overlap "<<MMPR4(k,d,v,m[d])) } 
if (m.find(d)!=m.end()) { ++i;  MM_ERR(" keys overlap "<<MMPR4(k,d,v,m[d][0])) } 
for(IdxTy i=0; i<szv; ++i) { 
const StrTy & v=(*ii).second[i];  
const StrTy  vc=CanonicalValue(v);
//m[d]=vc;
m[d].push_back(vc);
}
} // m_map
m_map=m;
MM_LOOP(ii,m_order_vector) { (*ii)=Canonical((*ii)); }
return i;
}
StrTy GetLastNames(StrTy & name, V & blast, const StrTy & s, const IdxTy flags)
{
typedef typename Fields::parsed_names_class Pnc;
Pnc pnc;
Fields x;
//StrTy name;
//const StrTy ali=get("author");
const StrTy ali=s; // get("author");
x.parse_names(pnc,ali,0);
IdxTy ia=0;
if (pnc.size()) name="";
MM_LOOP(ii,pnc)
{
const StrTy & n=(*ii).last();
if (n.length()==0) continue;
blast.push_back(n);
if (ia<3) 
{
if (m_name.length()) name+="_";
name+=n;

} // ia 
++ia;
} // ii 
return name; 
} // GetLastNames

// TODO this is making invalid names when there are not usable entries FIXME
void MakeName(TempProg & tp) { 
V t,a,b,blast; // u;
//MM_LOOP(ii,m_map) { MM_MSG(MMPR2((*ii).first,(*ii).second)) }
// adding zero will add spurs with zero length 
//StrTy yr=m_map["year"][0];
//StrTy yr=m_map["year"];
StrTy yr=get("year");
//m_gen.make_words(t,m_map["title"]);
m_gen.make_words(t,get("title"));
//m_gen.make_words(a,m_map["authors"]);
m_gen.make_words(a,get("authors"));
//m_gen.make_words(b,m_map["author"]);
m_gen.make_words(b,get("author"));
Fields f;
//StrTy GetLastNames(StrTy & name, V & blast, const StrTy & s, const IdxTy flags)
StrTy sname=get("author");
GetLastNames(m_name,blast,sname,0);
sname=get("editor");
GetLastNames(m_name,blast,sname,0);

/*

{
typedef typename Fields::parsed_names_class Pnc;
Pnc pnc;
Fields x;
const StrTy ali=get("author");
x.parse_names(pnc,ali,0);
IdxTy ia=0;
if (pnc.size()) m_name="";
MM_LOOP(ii,pnc)
{
const StrTy & n=(*ii).last();
if (n.length()==0) continue;
blast.push_back(n);
if (ia<3) 
{
if (m_name.length()) m_name+="_";
m_name+=n;

} // ia 
++ia;
} // ii 
} // scoping
*/

//m_gen.make_words(u,m_map["url"]);
m_gen.append_interesting(m_name,t);
// this will ignore short names
//m_gen.append_interesting(m_name,blast);
if (blast.size()==0)
{
m_gen.append_interesting(m_name,a);
m_gen.append_interesting(m_name,b);
}

//m_gen.append_interesting_url(m_name,m_map["url"]);
if (m_name.length()<10) m_gen.append_interesting_url(m_name,get("url"));
m_name=m_name+"_"+yr;
MM_ERR(MMPR3(__FUNCTION__,yr,m_name))
MM_LOOP(ii,t) { MM_ERR(MMPR2(__FUNCTION__,(*ii))) } 
// need a char-class thing but this is ok for now...
const IdxTy len=m_name.length();
const char * p=m_name.c_str();
// this is always a risk with data dependent sizes unchecked...
char cp[len+1];
IdxTy j=0;
// TODO make static and add a bit for this 
ParseTable clut;
for(IdxTy i=0; i<len; ++i)
{
const char c=p[i];
if (c=='\\') continue;
if (c=='\'') continue;
if (c=='"') continue;
if (c=='{') continue;
if (c=='}') continue;
if (c==',') continue;
if (c=='-') continue;
if (c=='_') { cp[j]=c; ++j; continue; } 
const IdxTy iv=clut.lut(c);
if ((iv&(ParseTable::DARK))==0) continue;
if ((iv&ParseTable::ALPHA)) { cp[j]=c; ++j; continue; } 
if ((iv&ParseTable::DIGIT)) { cp[j]=c; ++j; continue; } 
if ( c&128) continue;
//if (!isalnum(c)) continue;
//cp[j]=c; ++j;
} //i 
cp[j]=0; // doh WTF 2021-06-02
m_name=StrTy(cp);
if (m_name.length()<3) {m_name+="short_name"; } 
}  // MakeNAme


void Dump(OsTy & ss,const IdxTy flags)
{
const StrTy sep=" ";
const StrTy lead=m_name+sep+m_type;
//MM_LOOP(ii,m_map) { ss<<lead<<sep<<(*ii).first<<sep<<(*ii).second<<CRLF; } 
MM_LOOP(ii,m_map) { 
const auto & v=(*ii).second;
const IdxTy sz=v.size();
if (sz==0)
{
MM_ERR(" zero length value is wrong, "<<MMPR2(lead, (*ii).first))
}
IdxTy i=0;
MM_LOOP(jj,v)
{ if (sz==1) { ss<<lead<<sep<<(*ii).first<<sep<<(*jj)<<CRLF; }
else {  ss<<lead<<sep<<(*ii).first<<"["<<i<<"]"<<sep<<(*jj)<<CRLF; }
++i;  }
} 

}

void DumpErrors(OsTy & ss,const IdxTy flags) const 
{
const StrTy sep=" ";
const StrTy lead=m_name+sep+m_type;
const bool info_if_errors=Bit(flags,0);
const bool info_too=Bit(flags,1);
const bool dump_info=info_too||(info_if_errors&&(m_errors.size()));
MM_ERR(" dumping rrrors "<<MMPR4(flags, lead,info_if_errors,info_too))
MM_ERR(" dumping rrrors "<<MMPR4(flags, m_errors.size(),m_observations.size(),info_too))
//MM_LOOP(ii,m_map) { ss<<lead<<sep<<(*ii).first<<sep<<(*ii).second<<CRLF; } 
MM_LOOP(ii,m_errors) { 
const auto & v=(*ii).second;
const IdxTy sz=v.size();
if (sz==0)
{
MM_ERR(" zero length value is wrong, "<<MMPR2(lead, (*ii).first))
}
IdxTy i=0;
MM_LOOP(jj,v)
{ if (sz==1) { ss<<lead<<sep<<(*ii).first<<sep<<(*jj)<<CRLF; }
else {  ss<<lead<<sep<<(*ii).first<<"["<<i<<"]"<<sep<<(*jj)<<CRLF; }
++i;  }
} 

if (dump_info)
{
MM_LOOP(ii,m_observations)
{
const auto & v=(*ii).second;
MM_LOOP(jj,v)
{
ss<<(*ii).first<<" "<<(*jj)<<CRLF; 

} // jj 

} // ii 

} // dump_info


} // DumpErrors


void Save(OsTy & ss,const IdxTy flags) const  { 
//Save(ss,flags); 
const bool use_map_order=(Bit(flags,0));
const bool include_comments=(Bit(flags,1));
//MM_ERR(" FCK "<<MMPR2(flags,use_map_order))
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
// { ss<<","<<CRLF<<spaces; Format(ss,*ii,m_map[*ii]); }
{
// making this const fcks the sht fck brakcet fck 
const auto kk=m_map.find(*ii);
if (kk!=m_map.end())
{
MM_LOOP(jj,(*kk).second)
{ss<<","<<CRLF<<spaces;  Format(ss,(*ii),(*jj)); }
} // kk 
//pTgtV * p=m_wov_map.find("name",(*ii));
//if (p!=0)
//{
//MM_LOOP(jj,*p)
//{ss<<","<<CRLF<<spaces;  Format(ss,(*ii),(*jj)); }

//} // p 

} // ii 
} // use_order
else
{
MM_LOOP(ii,m_map) 
{//{ss<<","<<CRLF<<spaces;  Format(ss,(*ii).first,(*ii).second); }
MM_LOOP(jj,(*ii).second)
{ss<<","<<CRLF<<spaces;  Format(ss,(*ii).first,(*jj)); }
}

} // use_order 
if (kvsz>0) ss<<CRLF;
ss<<"}"<<CRLF;

} // save

void FormatOld(OsTy & ss) const 
{
//MM_ERR(" start format ")
MM_LOOP(ii,m_comments) {
//MM_ERR(" comment "<<(*ii).second)
 //ss<<"%" << (*ii).second<<CRLF;
// now comments are multilkine with percent already added 
 ss<< (*ii).second<<CRLF;
 } 
ss<<"@"<<m_type<<"{"<<m_name; // <<CRLF;
const IdxTy sz=m_map.size();
IdxTy i=0; 
if ( i!=sz) ss<<","; ss<<CRLF;  
MM_LOOP(ii,m_map) 
{
//	++i; Format(ss,(*ii).first,(*ii).second);
//	if ( i!=sz) ss<<","; ss<<CRLF;  

// wth was this?? 
//++ii;
++i;


const auto & v=(*ii).second;
const IdxTy szj=v.size();
IdxTy j=0;
MM_LOOP(jj,v)
{
++j; Format(ss,(*ii).first,(*jj));
if (( i!=sz)||(j!=szj)) ss<<","; 
ss<<CRLF;  
} // jj 

}
ss<<"}"<<CRLF;
//MM_ERR(" done format ")
}
////////////////////////////////////////////////

template <class Trag> 
void FormatRag(Trag & r, const IdxTy flags) const 
{
const bool parse=Bit(flags,0);
typedef typename Trag::Line Line;
//r.code_split(0);
r.splitter(7);
r.dynamic_split();
Ss ss,rr;
MM_LOOP(ii,m_comments) { 
Line l;
l.push_back( (*ii).second); 
r.add(l);
} 
// the parsers have a problem if these are blacnk....
{
const StrTy t=m_type; // (m_type.length()==0)?missing_type():m_type;
const StrTy n=m_name; // (m_name.length()==0)?missing_name():m_name;
//ss<<"@"<<m_type<<"{"<<m_name; // <<CRLF;
//ss<<"@"<<t<<"{"<<n; // <<CRLF;
Line l;
l.push_back("@"); l.push_back(t); l.push_back(n);
r.add(l);
}

MM_LOOP(ii,m_map) 
{
const StrTy & key=(*ii).first;
// must always have char zero... 
//const bool ex=(key.c_str()[0]=='%'); 
const auto & v=(*ii).second;
//Ss & qq=ex?rr:ss;
MM_LOOP(jj,v)
{
//if (qq.str().length()!=0) { qq<<","<<CRLF; } 
// Format(qq,key,(*jj)); //  qq<<CRLF;  
if (parse)
{
StrTy l=key+" "+(*jj);
r.add(l);
} // parse 
else
{

Line l;
l.push_back(key); l.push_back(*jj);
r.add(l);
} // not parse

} // jj 
} // ii 
//ss<<CRLF<<"}"<<CRLF;
//_ss<<rr.str();
//_ss<<CRLF;
//_ss<<ss.str();
} // FormatRag


///////////////////////////////////////////////
void Format(OsTy & _ss) const 
{
Ss ss,rr;
MM_LOOP(ii,m_comments) { rr<< (*ii).second<<CRLF; } 
// the parsers have a problem if these are blacnk....

const StrTy t=(m_type.length()==0)?missing_type():m_type;
const StrTy n=(m_name.length()==0)?missing_name():m_name;
//ss<<"@"<<m_type<<"{"<<m_name; // <<CRLF;
ss<<"@"<<t<<"{"<<n; // <<CRLF;


MM_LOOP(ii,m_map) 
{
const StrTy & key=(*ii).first;
// must always have char zero... 
const bool ex=(key.c_str()[0]=='%'); 
const auto & v=(*ii).second;
Ss & qq=ex?rr:ss;
MM_LOOP(jj,v)
{
if (qq.str().length()!=0) { qq<<","<<CRLF; } 
 Format(qq,key,(*jj)); //  qq<<CRLF;  

} // jj 

} // ii 
ss<<CRLF<<"}"<<CRLF;
_ss<<rr.str();
_ss<<CRLF;
_ss<<ss.str();

} // Format

///////////////////////////////////////////////////
void Format(OsTy & ss, const StrTy & key, const StrTy & v) const 
{ ss<<key<<" = {"<<v<<"}"; }
//{ ss<<key<<" = {\""<<v<<"\"}"; }

void Init()
{
m_type="";
m_name="";
m_map.clear();
m_order.clear();
m_comments.clear();
m_serial=bad();
m_buf_start=~0;
m_buf_size=~0;
m_ampersand=~9;
m_insert_point=~0;

}
IdxTy  Parse(const char * p, const IdxTy len)
{
typedef mjm_bibtex_parse<Tr> Bp;
Bp bib_parse;

return bib_parse.parse(m_map,p,len,0);
}

StrTy m_type, m_name,m_lexi_added_date,m_source_file;
MmapTy m_map;
MmapTy m_errors,m_observations;
Ord m_order;
V m_order_vector;
MapTy m_comments;
Ngen m_gen;
IdxTy m_serial;
StrTy m_blank;
StrTy m_original_text;
IdxTy m_buf_start,m_buf_size, m_ampersand,m_insert_point;
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

typedef char Ch;
typedef mjm_bibtex_entry<Tr> Tgt;
typedef std::vector<Tgt> TgtV;
typedef std::map<StrTy, TgtV> TgtM;
typedef typename TgtM::const_iterator Mcitor;
typedef mjm_blob<Tr> Blob;
typedef mjm_read_buffer<Tr> RdBuf;
typedef mjm_var_parse<Tr> ParseTable;
typedef mjm_pawnoff<Tr>  Handler;
typedef mjm_wovdb<Tr,StrTy> WovDB;
typedef mjm_wovdb<Tr,Tgt> WovBibDB;
typedef typename WovBibDB::iter itor;
typedef typename WovBibDB::citer const_itor;
typedef typename WovBibDB::vector_type pTgtV;
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
// TODO this needs to allow zero size entries doh 
if (m_i>=(*m_ii).second.size()){  ++m_ii; m_i=0; }
return *this;
}
 operator bool () {  return m_ii!=m_map.end(); } 
void see_ok()
{
if (m_ii==m_map.end()) return; 
if ((*m_ii).second.size()==0)
MM_ERR(" aero lengtrh vector "<<(*m_ii).first);

}
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

}; // iter_ord
class iter_wov_ord
{
public:
iter_wov_ord(WovBibDB  & m, const OrdVec & o  ) : m_db_map(m), m_ord(o) {
m_sz=m_ord.size();
m_o=0;
bump();
//if ( m_sz>0) m_pv=m_db_map.find("name",m_ord[0]);

m_i=0;
}
const Tgt & operator*() { return m_db_map[(*m_pv)[m_i]]; }
iter_wov_ord & operator++() { ++m_i; 
if (m_i>=(*m_pv).size()){ ++m_o ;  
bump();
//while (m_o<m_sz){   m_pv=m_db_map.find("name",m_ord[m_o]);
//if (m_pv!=0) if ((*m_pv).size()!=0) break;
//++m_o; } 

 m_i=0; }
return *this;
}
 operator bool () {  return (m_pv!=0)&&(m_o<m_sz); } 

void bump()
{
while (m_o<m_sz){   m_pv=m_db_map.find("name",m_ord[m_o]);
if (m_pv!=0) if ((*m_pv).size()!=0) break;
++m_o; } 

}

WovBibDB& m_db_map;
const OrdVec  & m_ord;
IdxTy m_i,m_o,m_sz;
const pTgtV*  m_pv;

}; // iter_wov_ord


//Mcitor begin() const  { return m_map.begin(); } 
//Mcitor end() const  { return m_map.end(); } 

const_itor begin() const  { return m_wov_map.begin("name"); } 
//const_itor end() const  { return m_wov_map.end("name"); } 

//IdxTy size(const StrTy & s) const { auto ii=m_map.find(s);
//if (ii==m_map.end()) return 0; return (*ii).second.size(); }  
IdxTy size(const StrTy & s) const {const  pTgtV* p =m_wov_map.find("name",s);
if (p==0) return 0; return (*p).size(); }  
// typically there will only be one... 
const Tgt * operator[](const StrTy & nm) const  
{ const pTgtV * p=m_wov_map.find("name",nm);
if (p==0) return 0; if ((*p).size())  return & m_wov_map[(*p)[0]]; return 0; }

 Tgt * operator[](const StrTy & nm)   
{ const pTgtV * p=m_wov_map.find("name",nm);
if (p==0) return 0; if ((*p).size())  return & m_wov_map[(*p)[0]]; return 0; }

Tgt & operator[](const IdxTy i) { return m_wov_map[i]; } 
const Tgt & operator[](const IdxTy i) const  { return m_wov_map[i]; } 

void add(const Tgt & t) { Add(t); } 
void parse( const StrTy nm)  { Parse(nm); } 
void parse( const StrTy nm, const Ch * ptr, const IdxTy len){Parse_db(nm,ptr,len); } 
void dump(OsTy & os, const IdxTy flags) { Dump(os,flags); } 
StrTy dump(const IdxTy flags) { Ss ss; dump(ss,flags); return ss.str(); } 
void save(const StrTy  & fn, const IdxTy flags) { Save(fn,flags); } 
void command_mode(Ci & li) { CommandMode( li) ; } 
IdxTy cmd(const StrTy & p1, const StrTy & p2 )  { return Cmd(p1,p2); } 
void try_to_fix(Ci & li) {  TryToFix( li); } 
void check_comments(OsTy & os, const IdxTy flags) 
 {  CheckComments( os,  flags) ;} 
void missing(OsTy & os, const StrTy & f, const IdxTy flags) {Missing(os,f,flags); } 
// TODO this is not really right as it is number of unique names 
//IdxTy size() const { return m_map.size(); }
IdxTy size() const { return m_wov_map.size("name"); }
// this is now easy with wovdb
IdxTy entries() const { return m_wov_map.size(); }
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
//iter x(m_map);
//while (x) { if ((*x).missing(f)) { os<<(*x).type()<<" "<<(*x).name()<<CRLF;}  ++x; } // iter
MM_WOV_LOOP(ii,m_wov_map,"name")
{ if ((*ii).missing(f)) { os<<(*ii).type()<<" "<<(*ii).name()<<CRLF;}  } // iter

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
#if 0
iter_ord x(m_map,m_ord);
while (x) {  
os<<"%"<< serial<<" prior "<<(*x).serial()<<CRLF;
(*x).save(os,flags);  ++x; 
++serial; 
} // iter
#endif
iter_wov_ord ii(m_wov_map,m_ord);
while (ii) {  
os<<"%"<< serial<<" prior "<<(*ii).serial()<<CRLF;
(*ii).save(os,flags);  ++ii; 
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
//MM_ERR(" dumping "<<MMPR4(all_entries,summary,count,lines))
if ( all_entries) {
//iter x(m_map);
//while (x) {os<<CRLF;  os<<(*x).format(); os<<CRLF;  ++x; } // iter
MM_LOOP(ii,m_wov_map) {os<<CRLF;  os<<(*ii).format(); os<<CRLF;  } // iter

}
if (summary)
{
//os<<MMPR(m_map.size())<<CRLF;
os<<MMPR(m_wov_map.size())<<CRLF;
//MM_LOOP(ii,m_map) 
//{ os<<"name="<<(*ii).first<<" size="<<(*ii).second.size()<<CRLF; }
//MM_LOOP(ii,m_wov_map.find("name")) 
//for(auto ii=m_wov_map.begin("name"); (ii); ii.next())
MM_WOV_KLOOP(ii,m_wov_map,"name")
{ os<<"name="<<(ii).key()<<" size="<<ii.vsize()<<CRLF; }


}
if (count)
{
typedef std::map<StrTy, IdxTy> Cm;
Cm m;
IdxTy n=0;
//iter x(m_map);
//while (x) { (*x).count(m); ++m[StrTy("@")+(*x).type()]; ++x; }
//for(auto ii=m_wov_map.begin("name"); (ii); ++ii)
MM_WOV_LOOP(ii,m_wov_map,"name")
{ (*ii).count(m); ++m[StrTy("@")+(*ii).type()];  }

MM_LOOP(ii,m) { os<<" field="<<(*ii).first<<" "<<(*ii).second<<CRLF; ++n; } 
} // count
if (lines)
{
IdxTy n=0;
//iter x(m_map);
//while (x) { os<<(*x).dump(0); ++x; }
//for(auto ii=m_wov_map.begin("name"); (ii); ++ii)
MM_WOV_LOOP(ii,m_wov_map,"name")
 {os<<(*ii).dump(0); }


//MM_LOOP(ii,m) { os<<" field="<<(*ii).first<<" "<<(*ii).second<<CRLF; } 
} // count



} // Dump

 
void Parse( const StrTy nm)
{
Blob x;
x.load(nm);
IdxTy len=x.size();
//MM_MSG(MMPR2(nm,len));
Parse_db(nm,x.ptr(),len);
}
// TODO  now to make this vector valued 
class asgn {
typedef std::map<StrTy,StrTy> M;
typedef WovDB Mdb;
public:
asgn(M & m) : m_map(m) {}
void kv(const StrTy & k, const StrTy & v, int stage) {m_map[k]=v; }
M &  m_map;

}; // asgn
class asgndb {
typedef std::map<StrTy,StrTy> M;
typedef WovDB Mdb;
public:
asgndb(Mdb & m) : m_map(m) {}
void kv(const StrTy & k, const StrTy & v, int stage) {
/*
grep m.kv mjm_bibtex_parse.h | sed -e 's/.*m\.kv/m.kv/'
m.kv(leading_name,rdb.string(),0); 
m.kv(type_name,rdb.string_sane(0),1);     return pc;  } 
m.kv(type_name,rdb.string_sane(0),2);     return pc;  } 
m.kv(type_name,rdb.string_sane(0),3); 
m.kv(type_name,rdb.string_sane(1),4); return pc; } 
m.kv(name_name,rdb.string_sane(1),5);
m.kv(rdb.string_sane(kvc-2),rdb.string_sane(kvc-1),6); 
*/
//m_map[k]=v; 
IdxTy loc=m_map.add(v);
if (stage==~0)  { m_map.index(loc,"error",k); return ; }
if (stage==~1)  { m_map.index(loc,"observe",k); return ; }
const StrTy idxname=(stage==6)?"":"special";
m_map.index(loc,idxname,k);

}
Mdb &  m_map;

}; // asgndb


void Parse_old( const StrTy nm, const Ch * ptr, const IdxTy len)
{
typedef std::map<StrTy,StrTy> M;
typedef WovDB Mdb;
typedef mjm_bibtex_parse<Tr> Bp;
Bp bib_parse;
IdxTy pos=0;
m_loaded_from=nm;
//MM_ERR(" parsing "<<MMPR2(nm,len))
while (pos<len)
{
// TODO FIXME this needs to be a vector valued map like
// the db although just adding an = val to be a push would work. 
M m;
Mdb mdb;
Tgt t;
//IdxTy dx= bib_parse.parse(m,x.ptr()+pos,len-pos,0);
asgn x(m);
asgndb xdb(mdb);
IdxTy dx= bib_parse.parse(x,ptr+pos,len-pos,0);
pos+=dx;
if (dx==0)
{
MM_ERR(" infininte loop possible "<<MMPR3(pos,len,(ptr+pos)))
break;

}
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


IdxTy  Parse_db( const StrTy nm, const Ch * ptr, const IdxTy len)
{
typedef WovDB Mdb;
typedef typename  Mdb::vector_type Vt;
typedef mjm_bibtex_parse<Tr> Bp;
Bp bib_parse;
IdxTy pos=0;
m_loaded_from=nm;
//MM_ERR(" parsing "<<MMPR2(nm,len))
while (pos<len)
{
// TODO FIXME this needs to be a vector valued map like
// the db although just adding an = val to be a push would work. 
Mdb mdb;
Tgt t;
//IdxTy dx= bib_parse.parse(m,x.ptr()+pos,len-pos,0);
asgndb xdb(mdb);
IdxTy dx= bib_parse.parse(xdb,ptr+pos,len-pos,0);
t.text(ptr,pos,dx);
pos+=dx;
if (dx==0)
{
MM_ERR(" infininte loop possible "<<MMPR3(pos,len,(ptr+pos)))
break;
}
// may get spurs at end or elsewhere

const Vt * vtype=(mdb.find("special","type_name"));
//const bool has_type=(vtype);
//const bool has_name=(mdb.find("special","name_name"));
const Vt * vname=(mdb.find("special","name_name"));
//const bool has_name=(vname);
const Vt * vlead=(mdb.find("special","junk_name"));
const IdxTy msz=mdb.size();
// the map size can't be zero if it has a name or type
// just for thought 
if (msz==0) continue;
//if (!has_name) continue;
// remove name, type, and comments then add to map
t.name(vname,mdb); // t.name(m["name_name"]);
t.type(vtype,mdb); // t.type(m["type_name"]);
t.comment(vlead,mdb); // t.comment(m["junk_name"]);
//MM_MSG(MMPR4(t.name(),pos,len,size()))
//MM_ERR(MMPR2((*vtype)[0],(*vname)[0]))
//m.erase(m.find("name_name"));
//m.erase(m.find("type_name"));
//m.erase(m.find("junk_name"));
t.serial(m_entries);
t.add_db(mdb);
std::vector<StrTy> v;
bib_parse.get_last_keys(v);
t.set_order_vector(v);
t.source(nm);
//IdxTy dx=t.parse(m,x.ptr()+pos,len-pos,0);
//os<<MMPR3(pos,dx,len);
//MM_LOOP(ii,m) { os<<MMPR2((*ii).first, (*ii).second)<<CRLF; }
add(t);
}
return pos;
} // Parse_db

void TryToFix(Ci & li)
{

//for(auto ii=m_wov_map.begin("name"); (ii); ii.next())
MM_WOV_LOOP(ii,m_wov_map,"name")
{IdxTy rc=TryToFix(*ii,li);if (rc==1) break;  } // iter

//iter x(m_map);
//while (x) {IdxTy rc=TryToFix(*x,li);if (rc==1) break;  ++x; } // iter
// Duplicates
//MM_LOOP(ii,m_map)
//{ if ((*ii).second.size()!=1) ResolveDuplicates(li,(*ii).first); }

MM_WOV_KLOOP(ii,m_wov_map,"name")
{ if (ii.vsize()!=1) ResolveDuplicates(li,ii.key()); }

}
void ResolveDuplicates(Ci &li, const StrTy &nm)
{
//TgtV orig=m_map[nm];
//MM_LOOP(ii,orig)
//{ Tgt & x=(*ii); 
//const StrTy & authors= x.get("authors");
//MM_MSG(MMPR4(x.serial(),x.name(),x.type(),authors)) }
//typedef typename WovBibDB::vector_type pTgtV;
// ptr could be null except only valid names passed 
const pTgtV &   orig=*(m_wov_map.find("name",nm));
MM_LOOP(ii,orig)
{ Tgt & x=m_wov_map[(*ii)]; 
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
#if 0 
iter x(m_map);
while (x) {
const IdxTy olap=(*x).make_keys_uniform();
olaps+=olap;
const IdxTy rc=IntegrateStructuredComments(*x);if (rc!=0) return rc;  ++x; } 
#endif


//iter x(m_map); while (x) 
MM_WOV_LOOP(ii,m_wov_map,"name")
{
const IdxTy olap=(*ii).make_keys_uniform();
olaps+=olap;
const IdxTy rc=IntegrateStructuredComments(*ii);if (rc!=0) return rc;  } 

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
if ( c!='%') {
MM_ERR(" adding percent")
 nc[pc]='%'; ++pc; ++mods; }
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
#if 0
iter x(m_map);
//while (x) { os<<(*x).format(); ++x; } // iter
while (x) {StrTy p=(*x).comment();  CheckComment((*x).name(),p.c_str(),p.length());  ++x; } // iter
#endif
MM_WOV_LOOP(ii,m_wov_map,"name")
{StrTy p=(*ii).comment();  CheckComment((*ii).name(),p.c_str(),p.length());  } // iter

} // all_entries

} // CheckComments

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
//m_map[t.name()].push_back(t);
m_ord.push_back(t.name()); 

const IdxTy loc=m_wov_map.add(t);
m_wov_map.index(loc,"name",t.name());

++m_entries;
 } 



//TgtM m_map;
WovBibDB m_wov_map;
IdxTy m_entries;
OrdVec m_ord;
ParseTable m_clut;
StrTy m_loaded_from;

}; // mjm_bibtex_entry_map



#endif
