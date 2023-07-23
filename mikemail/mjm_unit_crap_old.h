#ifndef MJM_UNIT_CRAP_H__
#define MJM_UNIT_CRAP_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

// should just be ragged table but allows more saving. 
//#include "mjm_ragged_forms.h"
#include "mjm_collections.h"
#include "mjm_string_tokenizer.h"
#include "mjm_read_buffer.h"
#include "mjm_indexed_map.h"
//#include "mjm_constrained_file.h"

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>


// Sun Aug 30 10:55:15 EDT 2020
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_unit_crap   
// g++ -Wall -std=gnu++11 -DTEST_MJM_UNIT_CRAP -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_unit_crap.h  -lpthread -lreadline





template <class Tr>
class mjm_unit_crap 
{
 typedef mjm_unit_crap Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef string_tokenizer St;
typedef St::coded_type Co;
typedef mjm_read_buffer<Tr>  Rb;

typedef std::map<StrTy,int> Powers;
class _parse_rational_expression
{

typedef  _parse_rational_expression Myt;
public:
_parse_rational_expression() { }
// TODO m_sfx botched with operations 
_parse_rational_expression(const char * exp) { m_sfx=exp;  Parse(exp); }
_parse_rational_expression(const StrTy &  exp) { Parse(exp.c_str()); }
Myt operator/(const Myt & that) const { return Div(that); } 
Myt operator*(const Myt & that) const { return Div(that,1); } 
int  ratio(const Myt & that) const { return Ratio(that); } 
const StrTy&  sfx() const { return m_sfx; }
StrTy  usfx() const { return Sfx(); }
StrTy  dump() const { return Dump(); }
bool equiv(const Myt & that) const { return Equiv(that); } 
private:
StrTy Dump() const
{
//const Powers & p=m_pow;
Ss ss; 
ss<<usfx(); 
//MM_LOOP(ii,p){ss<<m_st((*ii).first) <<"^"<<(*ii).second;  }
//MM_LOOP(ii,p){ss<<((*ii).first) <<"^"<<(*ii).second;  }
return ss.str(); 
} // DUmp
bool Equiv(const Myt & that) const 
{ return this->E(that)&&that.E(*this); }  // Equiv
bool E(const Myt & that) const { 
const Powers & p=m_pow;
const Powers & q=that.m_pow;
//bool e=true;
MM_LOOP(ii,p) { int c=(*ii).second; 
if (!c) continue;
const auto jj=q.find((*ii).first);
if (jj==q.end()) return false;
if (c!=(*jj).second) return false;
} // ii
return true;
} // E

// keep positives on left, negatives precede with slash. 
// no negative exponents 
StrTy Sfx() const
{
const Powers & p=m_pow;
std::map<int, std::vector< StrTy > > omap;
MM_LOOP(ii,p){omap[-(*ii).second].push_back((*ii).first);  }
Ss ss; 
IdxTy ip=0;
MM_LOOP(ii,omap)
{
int n=-(*ii).first;
if (n==0) continue;
if (n==1) {MM_LOOP(jj,(*ii).second) 
	{  if (ip!=0) ss<<"-"; ss<<(*jj) ; ++ip; } continue;}
if (n>1) {MM_LOOP(jj,(*ii).second) 
	{  if (ip!=0) ss<<"-"; ss<<(*jj)<<"^"<<n ; ++ip; } continue;}
if (n== -1) {MM_LOOP(jj,(*ii).second) 
	{  if (ip==0) ss<<"1"; ss<<"/"<<(*jj) ; ++ip; } continue;}
{MM_LOOP(jj,(*ii).second) 
{  if (ip==0) ss<<"1"; ss<<"/"<<(*jj)<<"^"<<(-n ); ++ip; } continue;}
}// ii  
return ss.str(); 
} // Sfx



class ParseTree
{
public:
ParseTree():m_op(0),l(0),r(0) {}
ParseTree(const char * x ):m_op(0),l(0),r(0) {m_p[x]=1; }
ParseTree(const StrTy &  x ):m_op(0),l(0),r(0) {m_p[x]=1; }
ParseTree(const StrTy &  x,ParseTree * _l,ParseTree * _r )
:m_op(0),l(_l),r(_r) {m_p[x]=1; }
void newleft() { l = new ParseTree(); }
void freeleft() {delete  l ; l=0;  }
void newright() { r = new ParseTree(); }
void freeright() {delete  r ; r=0;  }
void addleft() {MM_LOOP(ii,l->m_p) { m_p[(*ii).first]+=(*ii).second; }  }
void addright() {MM_LOOP(ii,r->m_p) { m_p[(*ii).first]+=(*ii).second; }  }
void subleft() {MM_LOOP(ii,l->m_p) { m_p[(*ii).first]-=(*ii).second; }   }
void subright() {MM_LOOP(ii,r->m_p) { m_p[(*ii).first]-=(*ii).second; }   }
void mult(const int p ) {MM_LOOP(ii,m_p) { m_p[(*ii).first]*=p; }  }
void add(const StrTy & n, const int d) { m_p[n]+=d; }
const IdxTy size() const { return m_p.size(); }
int  val() const 
{
auto ii=m_p.begin();
if (ii!=m_p.end())  return atoi((*ii).first.c_str()); 
return 0;
}
char m_op;
Powers m_p;
ParseTree *l, *r;

}; // ParseTree

class ParseState
{
typedef ParseState Myt;
public:
ParseState(): e(0),i(0),len(0) {}
ParseState(const char * e) {this->e=e;  i=0; len=strlen(e); }
ParseState(ParseState * ps) { e=ps->e; i=ps->i; len=ps->len; }
const char c() const { return e[i]; } 
Myt & operator++() { ++i; return *this; } 
operator bool() const { return i<len; }   
const char * e;
IdxTy i,len;
}; // ParseState
// treat - as * 
// recognize -,*,/,^ should use charclass lut however. 
// needs to accept () 

void Add(ParseTree * pt, Rb & nm, const int d=1)
{
nm.cap();
if (nm[0][0]!=0) pt->add(StrTy(nm[0]),d); 
nm.clear();
}
int number( const char * exp, ParseState & ps )
{
int j=0;
bool neg=false;
int i=0;
while (ps)
{
const char c=ps.c(); 
if (c==0) return neg?(-j):j;
if (i==0) if (c=='-'){ ++ps;  ++i;  neg=true; continue; } 
if (c<='9') if (c>='0') { ++ps;  ++i; j=j*10+(c-'0'); continue; } 
return neg?(-j):j;

} // true;
return neg?(-j):j;
} // number 
void Partial(ParseTree * pt, const char * exp, ParseState & ps )
{
Rb rb;
rb.clear();
rb.mark();
while (ps)
{
const char c=ps.c();
++ps;
if (c==' ') {  continue; } 
if (c=='(') { Add(pt,rb);  pt->newleft(); Partial(pt->l,exp,ps); 
pt->addleft();
pt->freeleft();  }
else if (c==')' )  {Add(pt,rb); } 
else if (c=='-')  { Add(pt,rb); } 
else if (c=='*')  { Add(pt,rb); } 
else if (c=='^')  { 
Add(pt,rb,1);
//pt->newright(); 
int d=number(exp,ps);
//Partial(pt->r,exp,ps);  // this needs to eval to number
//int d=pt->r->val();
pt->mult(d);
//Add(pt,rb,d); 
//pt->freeright(); 
}  // ^
else if (c=='/')  { 
Add(pt,rb);
pt->newright(); Partial(pt->r,exp,ps);  // this needs to eval to number
pt->subright();
pt->freeright(); }  // ^
else rb.append(c);
} // ps 
Add(pt,rb);
if ((pt->l)|| (pt->r)) MM_ERR(" still have kids ")

} // Partial 
void Parse(const char * exp)
{
ParseState ps(exp);
ParseTree pt;
Partial(&pt,m_exp.c_str(),ps);
Ss ss; 
MM_LOOP(ii,(pt.m_p)) {  ss<<(*ii).first<<"^"<<(*ii).second; }
//MM_ERR("parsing "<<MMPR(ss.str()))
m_pow=pt.m_p;
}

Myt Div(const Myt & that, const int p=-1 ) const 
{ 
Myt x=*this;
x.m_exp="/";
//MM_LOOP(ii,that.m_pow) { x.m_pow[(*ii).first]-=(*ii).second; }
MM_LOOP(ii,that.m_pow) { x.m_pow[(*ii).first]+=(*ii).second*p; }
return x; 
} // Div

int  Ratio(const Myt & that) const 
{ 
Myt x=*this; // avoids inserting zero into map 
int low=0, high=0,rrx=0;
MM_LOOP(ii,that.m_pow) {
const int p1=(*ii).second;
if (p1==0) continue;
// this creates a zero if none before , using a copy
const int p2= x.m_pow[(*ii).first]; 
if (p2==0) continue;
const int ap1=(p1<0)?(-p1):p1;
const int ap2=(p2<0)?(-p2):p2;
const int r=(ap1>ap2)?(p1/p2):(p2/p1);
if (r>high) high=r;
if (r<low) low=r;
}
if(high!=0) if (low==0) rrx=high;
if(low!=0) if (high==0) rrx=low;
return rrx; 
} // Div


StrTy m_sfx; // questionable value, make the power map more readable.  
StrTy m_exp;
StrTy m_nanme;
public:
Powers m_pow;

}; // _parse_rational_expression
public:
typedef _parse_rational_expression ratio_parse;
protected:
class _unit_desc
{
typedef _unit_desc Myt;
typedef std::map<Co, Co > StrMap;
typedef std::map<StrTy, IdxTy > Aliases;
typedef std::map<Co, D > ConvMap; // qty???
typedef std::map<StrTy,int> Umap;

public:


typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef StrTy Word;
typedef string_tokenizer St;
// BAD refers to the quality of the abbreviation not the
// computer issues 
enum { ERROR, MULTIPLIER, BAD, EXTENSIVE, INTENSIVE };
_unit_desc() {}
_unit_desc(const Line & l, St & st) { Init(l, st); } 
const IdxTy & name() const { return m_name; }
const Aliases & labels() const { return m_aliases; } 
bool error() const { return m_error; } 
const StrTy &  error_string() const { return m_error_string; } 
bool unit() const { return false; }
bool multiplier() const { return false; }
bool acronym() const { return false; }
bool mfg() const { return false; }
// caller has tokenizer if needed 
IdxTy type() const { return m_type; }
const StrTy & expansion() const  { return m_def; } 
template <class Tmap > 
bool bridge(D & v, const Myt & that, const Tmap & m ) const { return Bridge(v,that,m); } 
void dump(Ss & ss , const IdxTy flags=0) { Dump(ss,flags); } 
void dump(Ss & ss , St &st , const IdxTy flags=0) { Dump(ss,st,flags); } 
private:
//StrTyZZ

bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
bool U(IdxTy & x, const Line & l, IdxTy & i, const IdxTy sz,St& st)
{ ++i; if (i>=sz) return false; x=st(l[i]); return true; }
bool U(StrTy & x, const Line & l, IdxTy & i, const IdxTy sz)
{ ++i; if (i>=sz) return false; x=l[i]; return true; }
bool UMap( const Line & l, IdxTy & i, const IdxTy sz, St& st)
{ ++i; if (i>=sz) return false; const StrTy & w=l[i]; 
Rb rb;
std::vector<Word> v;
rb.vectorize(v,w,',');
if (v.size()!=2)
{
m_error=true;
m_error_string=" convert needs a name,value pair not "+w; 
return false;
} // ==2
m_conv[st(v[0])]=atof(v[1].c_str()); 
return true; }
bool Labels( const Line & l, IdxTy & i, const IdxTy sz, St& st)
{ ++i; if (i>=sz) return false; const StrTy & w=l[i]; 
Rb rb;
std::vector<Word> v;
rb.vectorize(v,w,',');
MM_LOOP(ii,v) { m_aliases[(*ii)]=1;  } 
return true; }

bool UnitsRatio( const Line & l, IdxTy & i, const IdxTy sz, St& st)
{ ++i; if (i>=sz) return false; const StrTy & w=l[i]; 
ratio_parse rp=ratio_parse(w);
m_units=rp.m_pow;
return true; }




template <class Tmap > 
bool Bridge(D & v, const Myt & that, const Tmap & m ) const { 

// first see if they have a common base unit and types
// are compatible
if (m_base==that.m_base) 
{v=m_conv[that.m_base]/that.m_conv[m_base]; return true; } 
if (m_type!=that.m_type) return false; 
// try to find the conversion things in the map...


return false; 

} 

Co BadCo() { return ~0; } 

void Init(const Line & l, St & st ) 
{ 
m_error=false;
m_type=BadCo();
m_base=BadCo();
m_up=BadCo();
m_down=BadCo();
Rb m_scratch;
const IdxTy sz=l.size();
if (sz==0) return ; 
m_name=st(l[0]);
for(IdxTy i=0; i<sz; ++i)
{
if (m_error) break; 
const Word & w=l[i];
if (w=="text"){  U(m_text,l,i,sz); continue;}
if (w=="definiition"){  U(m_def,l,i,sz); continue;}
// 
if (w=="type"){  U(m_type,l,i,sz,st); continue;}
if (w=="system"){  U(m_system,l,i,sz); continue;}
if (w=="base"){  U(m_base,l,i,sz,st); continue;}
if (w=="up"){  U(m_up,l,i,sz,st); continue;}
if (w=="down"){  U(m_down,l,i,sz,st); continue;}
if (w=="units"){  U(m_units_string,l,i,sz); continue;}
if (w=="convert"){  UMap(l,i,sz,st); continue;}
if (w=="aliases"){  Labels(l,i,sz,st); continue;}
if (w=="units"){  Labels(l,i,sz,st); continue;}

} // i 

}  // Init 
void Dump(Ss & ss , const IdxTy flags=0)
{
ss<<MMPR4(m_name,m_text,m_def,m_type);
ss<<MMPR4(m_sfx_map.size(),m_conv.size(),m_error,m_error_string); 

}
void Dump(Ss & ss , St & st, const IdxTy flags=0)
{
ss<<MMPR4(st(m_name),m_text,m_def,st(m_type));
ss<<MMPR4(m_sfx_map.size(),m_conv.size(),m_error,m_error_string); 

}



// members 
Co m_name;
Aliases m_aliases;
StrTy m_text,m_def,m_units_string, m_system;
Co m_type,m_base,m_up,m_down ;
StrMap m_sfx_map;
ConvMap m_conv; // conversion factors for various names 
bool m_error;
StrTy m_error_string;
Umap  m_units;

friend mjm_unit_crap;

}; // _unit_desc

class _qty_dim
{
public:
// uggh thought I could avoid these doh
_qty_dim() : m_d(0) {}
_qty_dim(const D & v ) : m_d(v) {}
_qty_dim(const IdxTy & v ) : m_d(v) {}
_qty_dim(const int  & v ) : m_d(v) {}

_qty_dim(const char * exp)
{ char * last; m_d=strtod(exp,&last);if (last==exp) m_d=1;  m_rp=ratio_parse(last); }

_qty_dim(const D & v,const char * exp) { m_d=v; m_rp=ratio_parse(exp); } 
_qty_dim(const D & v,const StrTy &  exp) { m_d=v; m_rp=ratio_parse(exp); } 
_qty_dim(const StrTy &  exp)
{ char * last; m_d=strtod(exp.c_str(),&last); 
if (exp.c_str()==last) m_d=1;
m_rp=ratio_parse(last); }

D & qty() { return m_d; }
const D & qty() const { return m_d; }
StrTy sfx() const { return m_rp.sfx(); } 
StrTy dump() const
{ Ss ss; ss<<m_d<<m_rp.dump(); return ss.str(); } // dump

D m_d;
ratio_parse m_rp;

}; // _qty_dim

public:
typedef _unit_desc unit_desc;
typedef _qty_dim qty_dim;
//typedef std::map<StrTy, unit_desc> units_map;
//typedef std::map<StrTy, unit_desc> units_map;
typedef mjm_indexed_map<Tr,unit_desc> units_map;
mjm_unit_crap() {}
~mjm_unit_crap() {}
void fn(const StrTy & fn ) { m_fn=fn; }
IdxTy load() { return Load();} 
IdxTy save() { return Save();} 
IdxTy parse() { return Parse();} 
IdxTy size() const { return m_map.size(); }
void clear() { m_map.clear(); m_defs.clear();  }
qty_dim base(const StrTy & s ) { qty_dim x(s); return ToBase(x); } 
IdxTy to_base(qty_dim & x) const { return ToBase(x); } 
void intensive (const StrTy & s )
{
ratio_parse rp= ratio_parse(s.c_str());

}
// this needs to look up abbreviations 
// as well as whole names. 
bool have(const StrTy & s) { return Have(s); }
IdxTy  to(D & r, const StrTy & to, const StrTy & from, const StrTy & v)
{ return To(r, to,from,v); } 
IdxTy  to(D & r, const char * to, const char *  from)
{ return To(r, to,from); } 

IdxTy  convert(D & d, const D & q, const StrTy &  to, const StrTy &  from,
const StrTy & rho)
{ return  Convert(d, q, to, from,rho); }
void dump(Ss & ss,const IdxTy flags=0) {  Dump(ss,flags); }
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss; 
Dump(ss,flags); 
 return ss.str(); }
void Dump(Ss & ss, const IdxTy flags=0)
{
MM_ERR(MMPR(m_map.size()))
MM_LOOP(ii,m_map) {//MM_ERR(" dumping "); 

 (*ii).second.dump(ss,m_st); ss<<CRLF;  }
}

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
IdxTy Load() { 
m_defs.load(m_fn);
return 0;

} 
IdxTy Save() { 
Ss ss;
ss<<m_defs.dump(3);
std::ofstream os(m_fn.c_str());
os<<ss.str();

return 0;

}
IdxTy Add(const Line & l)
{
unit_desc x(l,m_st);
if (x.error())
{
// if length is zero no error possible
Ss ss;
ss<<" error in line :";
MM_LOOP(ii,l) { ss<<" "<<(*ii);  }
ss<<CRLF;
if (l.size()>0)
MM_ERR(" bad units abbrv "<<MMPR2(l[0],x.error_string()))
else
MM_ERR(" error with zero length line wtf")
}
//m_map[m_st(x.name())]=x;
m_map.add_multi(m_st(x.name()),x.labels(),x);

} // Add
 
IdxTy Parse() { 
Ragged & r= m_defs;
IdxTy lno=0;
MM_LOOP(ii,r)
{
const Line & l =(*ii);
if (l.size()==0)  { ++lno; continue; } 
if (l[0].c_str()[0]=='#') { ++lno; continue; } 
//MM_ERR(" parsing line "<<MMPR2(lno,l.size()))
unit_desc x(l,m_st);
if (x.error())
{
// if length is zero no error possible
if (l.size()>0)
MM_ERR(" bad units abbrv "<<MMPR3(lno,l[0],x.error_string()))
else
MM_ERR(" error with zero length line wtf")
}
//m_map[m_st(x.name())]=x;
m_map.add_multi(m_st(x.name()),x.labels(),x);
++lno;
} // ii 
return 0;

} // Parse
// FIXME this is dumb, if it is not zero then what?
bool Have(const StrTy & s) { return (m_map[s]!=0); }
IdxTy To(D & r, const StrTy & to, const StrTy & from, const StrTy & v)
{
unit_desc * in=m_map[from];
unit_desc * out=m_map[to];
if ( !in) return 1;
if (!out) return 2;
const D vin=atof(v.c_str());
const IdxTy ib=(*in).m_base;
const IdxTy  ob=(*out).m_base;
if (ib!=ob) {MM_ERR(" not chasing convert  chains "<<MMPR4(to,from,ib,ob)) } 
// the tokenizer is generally a waste here but
// may be useful in other places, with the right typedefs
// a dummy can be used 
//const D f=(*in).m_conv[((*in).m_base)]/(*out).m_conv[((*out).m_base)];
const D f=(*in).m_conv[ib]/(*out).m_conv[ob];
r=f*vin;
return 0;
}

IdxTy To(D & r, const StrTy & to, const StrTy & from)
{
unit_desc * in=m_map[from];
unit_desc * out=m_map[to];
if ( !in) return 1;
if (!out) return 2;
const IdxTy ib=(*in).m_base;
const IdxTy ob=(*out).m_base;
if (ib!=ob) {MM_ERR(" not chasing convert  chains "<<MMPR4(to,from,m_st(ib),m_st(ob))) } 
//const D f=(*in).m_conv[((*in).m_base)]/(*out).m_conv[((*out).m_base)];
const D f=(*in).m_conv[ib]/(*out).m_conv[ob];
r=f;
return 0;
}
StrTy  Rat(const Powers & p)
{
Ss ss; 
//MM_LOOP(ii,p){ss<<m_st((*ii).first) <<"^"<<(*ii).second;  }
MM_LOOP(ii,p){ss<<((*ii).first) <<"^"<<(*ii).second;  }
return ss.str(); 
}
// Convert each entry into its "base" unit or common
// one to make harmonization easier 
IdxTy ToBase(D & v, Powers & p)
{
Powers pb;
MM_LOOP(ii,p)
{
unit_desc * descp=m_map[(*ii).first];
if (!descp) 
{MM_ERR(" unit not found die now "<<MMPR2((*ii).first,(*ii).second))
return ~0;  } 
// doh, this needs the power of the exponent... 
const int pi=(*ii).second;
if (pi==0) continue; // should clean the map 
const IdxTy ib=(*descp).m_base;
D vc=(*descp).m_conv[ib];
int pii=(pi<0)?(-pi):pi;
if (pi<0) vc=1.0/vc;
// note for high int power, sqrt may be best 
while  (pii>0) { v=v*vc; --pii; } 
//if (pi==1) v=v*vc; else if (pi== -1 ) v=v/vc; else if (pi==0) {}
//else v=v*::pow(vc,pi);
// the input units should condense to fewer base units. 
pb[m_st(ib)]+=pi; // (*ii).second;
} // ii 
p=pb;
return 0;
} // ToBase

IdxTy ToBase(qty_dim & qd) { return  ToBase(qd.m_d, qd.m_rp.m_pow); } 

// gg/dump
IdxTy  Convert(D & d, const D & q, const StrTy &  to, const StrTy &  from,
const StrTy & rho)
{
MM_ERR(MMPR4(q,to,from,rho))
if (to==from) { d=q; return 0; } 
if (to.length()==0) { d=q; return 0; } 
// need a member element for caching, can't use static 
// with multiple instances 
qty_dim toqd(1,to);
const IdxTy rc4=ToBase(toqd);
// NB q added here doh 
qty_dim fromqd(q,from);
//MM_ERR(MMPR(fromqd.dump()))
const IdxTy rc2= ToBase(fromqd);
MM_ERR(MMPR2(toqd.dump(),fromqd.dump()))
//if (fromqd.m_rp.m_pow==toqd.m_rp.m_pow) 
if (fromqd.m_rp.equiv(toqd.m_rp))
{  // cache c here 
const D c=fromqd.m_d/toqd.m_d;
d=c;   MM_ERR(MMPR4(d,q,toqd.m_d,fromqd.m_d)) 
return 0; } 
if (rho.length()==0) return 1; 
qty_dim rhoqd(rho.c_str());
MM_ERR(MMPR(rhoqd.dump()))
const IdxTy rc3= ToBase(rhoqd);
MM_ERR(MMPR2(toqd.m_d,rhoqd.dump()))
ratio_parse x=toqd.m_rp/fromqd.m_rp;
int rr=x.ratio(rhoqd.m_rp);
ratio_parse y;
D c=1.0/toqd.m_d*fromqd.m_d;
const D d1=c;
//MM_ERR(MMPR4(d1,q,toqd.dump(),fromqd.dump())) 
if (rr==1) {  y=x/rhoqd.m_rp; c=c*rhoqd.m_d; }
else if (rr== -1) {  y=x*rhoqd.m_rp; c=c/rhoqd.m_d; }
// cache c here 
d=c;
//MM_ERR(MMPR2(d,c)<<MMPR4(y.dump(),v,rr,rhoqd.m_d))
MM_ERR(MMPR2(d,c)<<MMPR4(y.dump(),toqd.m_d,rr,rhoqd.m_d))
return rc2|rc3|rc4; 
}

// members
St m_st;
StrTy m_fn;
Ragged m_defs;
units_map m_map;

}; // mjm_unit_crap

//////////////////////////////////////////////

template <class Tr>
class mjm_unit_crap_map : public std::map<typename Tr::StrTy, mjm_unit_crap< Tr > >  
{
 typedef mjm_unit_crap_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_unit_crap< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_unit_crap_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }

private:

void Init()
{
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;


}
return ss.str();
// return Dump(flags); 

}




private:

}; // mjm_unit_crap_map




////////////////////////////////////////////
#ifdef  TEST_MJM_UNIT_CRAP
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_unit_crap <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_UNIT_CRAP "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_unit_crap<Tr>  Myt;
//Myt x(argc,args);
Myt x;

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="fn") { x.fn(cip.p1); MM_ERR(x.dump()) }
if (cmd=="load") { x.load(); MM_ERR(x.dump()) }
if (cmd=="save") { x.save(); MM_ERR(x.dump()) }
if (cmd=="parse") { x.parse(); MM_ERR(x.dump()) }
if (cmd=="eval") { x.intensive(cip.p1); MM_ERR(x.dump()) }
if (cmd=="to") { double v=0;  IdxTy rc=x.to(v,cip.p1,cip.p2,cip.wif(3)); 
MM_ERR(MMPR4(rc,v,cip.p1,cip.p2))
MM_ERR(x.dump()) }
//IdxTy  to(D & r, const StrTy & to, const StrTy & from, const StrTy & v)
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_UNIT_CRAP_H__ 
