#ifndef MJM_TABULATED_FUNCTIONS_H__
#define MJM_TABULATED_FUNCTIONS_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_integrals.h"
#include "mjm_closed_families.h"
#include "mjm_rods_array.h"
#include "mjm_collections.h"
#include "mjm_string_tokenizer.h"

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>

// copied from : 

//g++ -std=gnu++11 -DTEST_MJM_TABULATED_FUNCTIONS  -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_tabulated_functions.h  -lpthread -lreadline


// Sat Jan  4 18:20:12 EST 2020
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_legendre_stuff   
// g++ -std=gnu++11 -DTEST_MJM_LEGENDRE_STUFF -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_legendre_stuff.h  -lpthread -lreadline




template <class Tr, class _CoefTy=double>
class mjm_eval_var_cache 
{
 typedef mjm_eval_var_cache Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef _CoefTy Coef; // want to allow rationals and arb precision int 
typedef std::vector<Coef> Pows;
// these may be slower than tree map 
typedef std::unordered_map<IdxTy, int > VarPowMap;
typedef std::unordered_map<IdxTy, Pows > VarValsMap;
const bool Bit(const IdxTy f, const IdxTy b) const { return ((1<<b)&f)!=0; }
public:
mjm_eval_var_cache(string_tokenizer* stp): m_stp(stp),m_lmax(~0) {}
string_tokenizer * st() { return  m_stp;}
void st(string_tokenizer * st) {  m_stp=st;}

// note that the string_tokenizer is not cleared, that is up to 
// owner... 
void clear() { m_map.clear(); } 
IdxTy lmax() { return m_lmax; } 
void lmax(const IdxTy l ) {  m_lmax=l; } 
// need var values AND scale factors for each (l,m)
template <class Tx, class Ty> 
//void load(const Ty & vars, const VarPowMap & vpm)
// vars is the user supplied variable values, vpm
// the maximum values needed. 
void load(const Ty & vars, const Tx & vpm, const IdxTy flags)
{
const bool dump=Bit(flags,0);
MM_ERR(" loading vars "<<MMPR(dump))
MM_LOOP(ii,vars)
{
const IdxTy v=(*m_stp)((*ii).first);
const Coef & c=(*ii).second;
Coef cp=1;
//MM_ERR(MMPR3((*ii).first,v,c))
const auto jj= vpm.find(v);
if (jj!=vpm.end())
{
const int pmax=(*jj).second;
Pows & pv=m_map[v];
//MM_ERR(MMPR3((*ii).first,pmax,c))
while (pv.size()<=pmax) { pv.push_back(cp); cp=cp*c; } 
}// jj 
else
{
// it is not an error to have extra vars  FIXME
MM_ERR(MMPR3(v,(*ii).first,c)<<" not found in max power map  ")
//MM_ERR(MMPR3(v,(*ii).first,c)<<" not found in max power map which has ")
//MM_LOOP(kk,vpm){ MM_ERR(" vpm has "<<MMPR3((*kk).first,((*m_stp)((*kk).first)),(*kk).second))
//}
} // not found in vpm
} // ii 
// make sure all vars are there however... 
MM_LOOP(kk,vpm){ 
const IdxTy v=(*kk).first;
const int p=(*kk).second;
if (m_map.find(v)==m_map.end())
{
MM_ERR(" no value for vpm listing  "<<MMPR3(v,((*m_stp)(v)),p))
}

} // kk
if (dump)
{
MM_LOOP(kk,vpm){ 
const IdxTy v=(*kk).first;
const int p=(*kk).second;
Coef c=Coef(0);
const auto ff=vars.find((*m_stp)(v));
const bool found= (ff!=vars.end()) ;
if (found) c=(*ff).second; 
if (found) MM_ERR(" "<<MMPR4(v,((*m_stp)(v)),c,p))
else MM_ERR(" missing "<<MMPR3(v,((*m_stp)(v)),p))
} // kk


} // dump 

} // load


const Coef & varpow(const IdxTy v, const int p) 
{ 
auto ii=m_map.find(v);
if (ii==m_map.end()) {

 MM_ERR(" not found "<<MMPR4(v,((*m_stp)(v)),p,m_map.size()))
MM_LOOP(jj,m_map) { 
MM_ERR(" have "<<MMPR3((*ii).first,((*m_stp)((*ii).first)),(*ii).second.size()))
}
 } 
if (p>=m_map[v].size()) { MM_ERR(" power too big "<<MMPR4(p,v,((*m_stp)(v)),m_map[v].size())) } 
return m_map[v][p]; 

}  
virtual Coef scale_factor(const IdxTy l, const int m,const int gh ) { return Coef(1); } //  =0;


VarValsMap m_map;
string_tokenizer * m_stp;
IdxTy m_lmax;
}; // mjm_eval_var_cache

template <class Tr, class _CoefTy=double>
class mjm_tabulated_functions 
{
 typedef mjm_tabulated_functions Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef _CoefTy Coef; // want to allow rationals and arb precision int 


typedef std::vector<Coef> Poly;
typedef mjm_integrals Mi;
typedef std::unordered_map<IdxTy, int > VarPowMap;
public:
typedef mjm_ragged_table Ragged;
private:
typedef mjm_eval_var_cache <Tr,Coef> EvalVarsBase;
template <class Td,int N=2> class ntri_index{
public:
// ordering is m=l, m=l-1, etc for ease of indexing... 
//typedef Desc data_type;
typedef Td data_type;
// FIXME this is not right but seems to work leave gaps
// although l==0 seems to find something valid wtf?
static IdxTy size(const IdxTy lmax)  { return (N*(lmax+1)*(lmax+2))>>1; } 
static IdxTy index(const IdxTy l,const int  m)  
{ return  size(l)+(~m); }  // signed arthimetic, bacward store offset lol. 
}; // ntri_index



/* Finally a table based system in which integrals can be computer
exactly and multiplied by variable parameters. Right now this is
the integral of a schmidt-quasi normalized function l' m' left
unnormalized, multiplied by l 1 also unnormalized. These are done
almost exactly depending on how pi is eventually handled.

This usage is not right though, it needs an expandable storage
like std::vector with the caching used for indexing. The cache
points to the index entries which point to the internal vector  


*/
class tab_term
{

//typedef std::map<IdxTy, int > VarPowMap;
public:
enum { BAD=~0};
IdxTy bad() const { return BAD; } 
tab_term() : m_l(BAD),m_m(BAD),m_lp(BAD),m_mp(BAD),m_gh(BAD),m_c(0)
,m_prior(BAD), m_here(BAD),m_next(BAD) {}
tab_term(const IdxTy l, const IdxTy m, const IdxTy lp, const IdxTy mp,
const IdxTy gh,
const Coef c) : m_l(l),m_m(m),m_lp(lp),m_mp(mp),m_gh(gh),m_c(c) 
,m_prior(BAD), m_here(BAD),m_next(BAD) {}
IdxTy l() { return m_l; } 
IdxTy m() { return m_m; } 
void var(const IdxTy v, const int  p) { m_map[v]+=p; }
int  pow(const IdxTy v) { return m_map[v]; } 

Coef evaluate( EvalVarsBase * ev )
{
static Coef cz=Coef(0);
//static Coef cu=Coef(1);
if (m_c==cz) return cz;
Coef c=m_c;
MM_LOOP(ii,m_map)
{
const IdxTy v=(*ii).first;
const int p = (*ii).second;
if (p!=0)  { 
const Coef & cp=ev-> varpow( v,  p) ;
// these are sparse enough 
if (cp==cz) return cz; 
c=c*cp ; // ev-> varpow( v,  p) ;
//MM_ERR(MMPR4(c,m_c,cp,v)<<MMPR2(p,ev->varpow(v,1)))
}
//if (c==cz) return c; 
}
//if (c!=cz) {  MM_ERR(" non zero result "<<MMPR(c))}
return c;
} // evaluate


template <class Td>const  Coef & extract(Td &  d,  EvalVarsBase * ev ) const 
{
static Coef cz=Coef(0);
//static Coef cu=Coef(1);
//if (m_c==cz) return cz;
Coef c=m_c;
MM_LOOP(ii,m_map)
{
const IdxTy v=(*ii).first;
const int p = (*ii).second;
if (p!=0)  { d[(*(*ev).st())(v)]=p; } 
//const Coef & cp=ev-> varpow( v,  p) ;
// these are sparse enough 
//if (cp==cz) return cz; 
//c=c*cp ; // ev-> varpow( v,  p) ;
//MM_ERR(MMPR4(c,m_c,cp,v)<<MMPR2(p,ev->varpow(v,1)))
//}
//if (c==cz) return c; 
}
return m_c;
//if (c!=cz) {  MM_ERR(" non zero result "<<MMPR(c))}
//return c;
} // evaluate




void update_max_pow(VarPowMap & m) const 
{
MM_LOOP(ii,m_map)
{
IdxTy v=(*ii).first;
int x= (*ii).second;
// this is dumb but can divide for now if needed... 
if (x>m[v])  m[v]=x;
else if ((-x)>m[v])  m[v]=-x;
}

}
// links integrated into data type...
IdxTy next() const { return m_next; } 
StrTy dump() const
{
 Ss ss;
ss<<MMPR4(m_l,m_m,m_lp,m_mp)<<MMPR(m_gh);
ss<<MMPR4(m_c,int(m_prior),int(m_here),int(m_next)); 
MM_LOOP(ii,m_map) { ss<<MMPR2((*ii).first, (*ii).second); } 
return ss.str(); 
}
StrTy dump(string_tokenizer & st ) const
{
 Ss ss;
ss<<MMPR4(m_l,m_m,m_lp,m_mp)<<MMPR(m_gh);
ss<<MMPR4(m_c,int(m_prior),int(m_here),int(m_next)); 
MM_LOOP(ii,m_map) {
 //ss<<MMPR2(st((*ii).first), (*ii).second);
 ss<<" "<<(st((*ii).first))<<"^"<< (*ii).second;
 } 
return ss.str(); 
}




// keep all info even if indexed for later asserts. 
int  m_l,m_m,m_lp,m_mp,m_gh;
// coef and powers map
Coef m_c;

VarPowMap m_map;
IdxTy m_prior,m_here,m_next;

}; // tab_term 

typedef tab_term TabTerm;
class tab_list_index_entry
{
enum { BAD=~0};
public:
tab_list_index_entry():m_first(BAD),m_last(BAD),m_lp(BAD),m_mp(BAD),m_gh(BAD) {}
IdxTy bad() const { return BAD; } 

IdxTy add_link(const IdxTy loc, const IdxTy lp, const int mp, const int gh )
{

if (m_lp==BAD) m_lp=lp;
if (m_mp==BAD) m_mp=mp;
if (m_gh==BAD) m_gh=gh;
if ( (m_lp==lp ) && (m_mp==mp ) && (m_gh==gh ))
{
}else MM_ERR(" index mismatch "<<MMPR3(m_lp,m_mp, m_gh)<<MMPR4(lp,mp,gh,loc))
const IdxTy rv=m_last;
if (m_first==BAD){  m_first=loc; m_last=loc; } 
else { m_last=loc; } 
return rv; 

}
StrTy dump()
{
Ss ss;
ss<<MMPR3(m_lp,m_mp,m_gh);
ss<<MMPR2(int(m_first),int(m_last));

return ss.str();

}
IdxTy m_first, m_last;
// m_lp should be positive but not necesarily others 
int  m_lp,m_mp,m_gh;

}; // tab_list_index_entry
typedef tab_list_index_entry  TabListIndexEntry;


typedef ntri_index<TabListIndexEntry,2> TabTermIndex;
typedef mjm_rods_array<Tr,TabTermIndex> TabTermBase;



class tab_term_cache : public TabTermBase
{
typedef TabTermBase Super;
typedef std::vector<TabTerm> TabTermVector;
typedef ntri_index<TabListIndexEntry,2> Tindex;
public:
typedef TabTerm tab_term; 
// why the FUDD
typedef typename Super::data_type data_type; // FUDD
tab_term_cache() : m_state(0),m_lpmax(0) {}
VarPowMap & vpm() {  return m_max_pow; }
IdxTy lpmax() const { return m_lpmax;}
IdxTy index_size() const { return Super::size();}
void add(const tab_term & t) {
if (m_state!=0) { MM_ERR(" already compiled "<<t.dump()) }
 m_vec.push_back(t); 
if (t.m_lp>m_lpmax) m_lpmax=t.m_lp;
t.update_max_pow(m_max_pow);
//MM_ERR(" term added "<<t.dump(&m_st))
} 
bool  next( tab_term & t) 
{
//MM_ERR(" next for "<<t.dump(m_st))
IdxTy n =t.next();
if (n==t.bad()) return false;
t=m_vec[n];
return true; 
}
IdxTy entries( const IdxTy lp, const int mp, const int gh )
{
IdxTy n=0;
tab_term t;
if(!(*this)(t,lp,mp,gh)) return 0;
do { ++n; } while (next(t));
return n;
}
IdxTy bad() const  { return ~0; } 
bool  operator()(tab_term & t, const IdxTy lp, const int mp, const int gh )
//{return (*this)(Tindex::index(l,m)); } 
{
if (m_state==0) { MM_ERR(" not compiled "<<MMPR3(lp,mp,gh)<<t.dump()) }
//const IdxTy x=(Tindex::index(lp,gh+(mp<<1))); 
const IdxTy x=(*this)(lp,mp,gh); // (Tindex::index(lp,gh+(mp<<1))); 
if (x>=index_size()) return false; 
if (x==bad()) return false; 
data_type & idx= ((Super*)(this))->operator()(x); 
//data_type & idx= ((Super*)(this))->operator()(Tindex::index(lp,gh+(mp<<1))); 
if ( idx.m_first==idx.bad()) return false;
t= m_vec[idx.m_first];
return true; 
} 
//IdxTy operator()( const IdxTy lp, const int mp, const int gh )
//{ return ((Super*)(this))->operator()(Tindex::index(lp,gh+(mp<<1)));  } 
IdxTy operator()( const IdxTy lp, const int mp, const int gh )
{ 
if (lp==bad()) return bad();
if (mp==bad()) return bad();
if (mp<0) return bad();
if (gh==bad()) return bad();
if (gh<0) return bad();

const IdxTy x= Tindex::index(lp,gh+(mp<<1));  
if (x>=index_size()) return bad(); 
return x;
} 


data_type &  operator()(const IdxTy n )
{
return ((Super*)(this))->operator()(n);  
} 
void compile()
{
const IdxTy sz=m_vec.size();
((Super*)(this))->size(Tindex::size(m_lpmax),0);
for(IdxTy i=0; i<sz; ++i)
{
tab_term & t= m_vec[i];
IdxTy loc=(*this)(t.m_lp,t.m_mp,t.m_gh);
data_type & ie=(*this)(loc);
const IdxTy last=ie.add_link(i,t.m_lp,t.m_mp,t.m_gh);
t.m_here=i;
if (last!=ie.bad())
{
m_vec[last].m_next=i; 
t.m_prior=last;
t.m_next=t.bad();
}

} // i 
m_state=1;

} // compile


StrTy dump(const IdxTy flags)
{
Ss ss;
ss<< "Called the right code "<<Super::used()<<CRLF;
// why the FUDD can't it find this shot 
const IdxTy sz=Super::used();
for(IdxTy i=0; i<sz; ++i)
{
ss<<i<<" "<<(*this)(i).to_string()<<CRLF;
}
return ss.str();
}
TabTermVector m_vec;
IdxTy m_state;
IdxTy m_lpmax;
VarPowMap m_max_pow;
}; // tab_term_cache 

typedef tab_term_cache TabTermCache;

typedef string_tokenizer St;

public:
St * st() { return & m_st; } 
template <class Td, class Tc > void extract(Td & d, Tc & coefs,  EvalVarsBase * ev, const IdxTy lp, const int  mp,const int gh, const IdxTy flags)
{  Extract( d, coefs, ev, lp,  mp, gh,  flags); } 

Coef evaluate(EvalVarsBase * ev, const IdxTy lp, const int  mp,const int gh, const IdxTy flags)
{ return Evaluate(ev,lp,mp,gh,flags); } 
StrTy list_terms(const IdxTy lp, const IdxTy mp,const IdxTy gh)
{ return ListTerms(lp,mp,gh); } 
/*
TabTerm t; // = m_ttc.
bool ok= m_ttc( t,  lp,  mp,  gh );
const IdxTy k= m_ttc.entries(  lp,  mp,  gh );
IdxTy n=0;
if (ok)
{
do{
MM_ERR(MMPR3(n,k,t.dump(m_st))) 
++n;
} while (m_ttc.next(t));

} // if
return "";
} // list_terms
*/
IdxTy bad() const { return ~0; } 
void load(Ragged & s,const IdxTy flags ) { Load( s,flags); } 
VarPowMap & vpm() { return m_ttc.vpm(); } 
StrTy dump_index(const IdxTy flags=0) { return DumpIndex(flags); }
StrTy dump(const IdxTy flags=0) { return Dump(flags); }

Coef Evaluate(EvalVarsBase * ev, const IdxTy lp, const int  mp,const int gh, const IdxTy flags)
{
Coef c=Coef(0);

TabTerm t; // = m_ttc.
bool ok= m_ttc( t,  lp,  mp,  gh ); 
const IdxTy k= m_ttc.entries(  lp,  mp,  gh );
IdxTy n=0;
if (ok)
{
do{
if (t.l()<=ev->lmax()) { // MM_ERR(MMPR3(n,k,t.dump(m_st)))
Coef tc=t.evaluate(ev); 
if (tc!=0) {
tc=tc*ev->scale_factor(t.l(),t.m(),gh);
c+=tc; }
}
++n;
} while (m_ttc.next(t));
//c=c*ev->scale_factor(l,m,gh);

}
return c; 
}
template <class Td,class Tc > void Extract(Td & d, Tc & coefs, EvalVarsBase * ev, const IdxTy lp, const int  mp,const int gh, const IdxTy flags)
{
ev->st(&m_st);
typedef typename Td::value_type Te;
TabTerm t;
bool ok= m_ttc( t,  lp,  mp,  gh ); 
const IdxTy k= m_ttc.entries(  lp,  mp,  gh );
IdxTy n=0;
if (ok)
{
do{
Te v;
v.clear();
if (t.l()<=ev->lmax()) { // MM_ERR(MMPR3(n,k,t.dump(m_st)))
//Coef tc=t.evaluate(ev); 
const Coef & c= t.extract(v,ev); 
//v["\\coef"]=c;
d.push_back(v);
//Ss ss;
//MM_LOOP(kk,v) { ss<<(*kk).first<<" "<<(*kk).second<<" "; }
coefs.push_back(c);
//MM_ERR(MMPR2(c,ss.str()))
// now this needs a way to add a scale factor to the extracted term 
//tc=tc*ev->scale_factor(t.l(),t.m(),gh);
}
++n;
} while (m_ttc.next(t));
}

} // Extract

StrTy ListTerms(const IdxTy lp, const IdxTy mp,const IdxTy gh)
{
TabTerm t; // = m_ttc.
bool ok= m_ttc( t,  lp,  mp,  gh ); 
const IdxTy k= m_ttc.entries(  lp,  mp,  gh );
IdxTy n=0;
if (ok)
{
do{
MM_ERR(MMPR3(n,k,t.dump(m_st)))
++n;
} while (m_ttc.next(t));

} // if
return "";

} // ListTerms
/*
// need to profile and re-run without " ana..." space and change h/g to number

 analyze_vars 1 1 1 1 h coef approx divpi scaled approxs salpha sbeta
analyze_terms 1 1 1 1 h -884279719003555/422212465065984 -2.094395102393195 -0.66
66666666666666 -1 -1 1 1

 analyze_vars 1 1 1 1 g coef approx divpi scaled approxs salpha cbeta
analyze_terms 1 1 1 1 g -884279719003555/422212465065984 -2.094395102393195 -0.66
66666666666666 -1 -1 1 1

 analyze_vars 1 1 2 0 h coef approx divpi scaled approxs
*/
bool Bit(const IdxTy f, const IdxTy b) const { return ((1<<b)&f)!=0; } 
void Load(Ragged & s,const IdxTy flags )
{
const bool debug=Bit(flags,0);
const IdxTy sz=s.size();
typedef std::map<StrTy, IdxTy> Vm;
typedef std::map< IdxTy,StrTy> Vim;
Vim m;
Vm im;
IdxTy pos=0;
IdxTy approxs=0;
int l=bad(),m_=bad(),lp=bad(),mp=bad(),gh=bad();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line line=s.line(i);
const IdxTy len=line.size();
if (debug) { MM_ERR(MMPR(len));
for(IdxTy i=0; i<len; ++i ) MM_ERR(MMPR2(i,line[i])); 
}
if (len<2) continue;
if (line[0]=="analyze_vars")
{
m.clear();
im.clear();
for(IdxTy j=1; j<len; ++j) 
{const StrTy & w=line[j]; im[w]=j;   m[j]=w;

if (w=="coef") pos=j;   
if (w=="approxs") approxs=j;   

} 
//{const StrTy & w=line[j]; im[w]=j;   m[j]=w;if (w=="approx") pos=j;   } 
l=atoi(line[1].c_str());
m_=atoi(line[2].c_str());
lp=atoi(line[3].c_str());
mp=atoi(line[4].c_str());
gh=atoi(line[5].c_str());
if(line[5]=="h") gh=1;

} // if 
else if (line[0]=="analyze_terms")
{
TabTerm t;

//tab_term(const IdxTy l, const IdxTy m, const IdxTy lp, const IdxTy mp, const IdxTy gh,
//const Coef c) : m_l(l),m_m(m),m_lp(lp),m_mp(mp),m_gh(gh),m_c(c) {}
//void var(const IdxTy v, const IdxTy p) { m_map[v]+=p; }

for(IdxTy j=1; j<len; ++j) 
{
const StrTy & w=line[j];
if (j<pos) // two loops? wtf
{
if (w!=m[j]) MM_ERR("mismatch in input "<<MMPR4(i,line[0],line[j],m[j]))
}// j<pos
else if (j==pos)  {}
else if (j==(pos+1)) 
{Ss ss; ss<<w; Coef c; Ss sc(ss.str()); sc>>c ;  t= TabTerm(l,m_,lp,mp,gh,c); }
else
{
// these should be ints... 
const int p=atoi(w.c_str());
//t.var(m_st(m[j]),atoi(w.c_str()));
if ( j>approxs) if( p!=0) t.var(m_st(m[j]),p);

} // j>=pos
} // j 
if (len>=(pos+1)) { 
//MM_ERR(" adding "<<t.dump(m_st))
m_ttc.add(t);
}
} else {MM_ERR( " unkonwn line "<<MMPR3(i,line[0],line[1])) } 

} // i 
// do not really want to do this here. 
m_ttc.compile();
} // Load

private:

StrTy Dump(const IdxTy flags)
{
Ss ss;


return ss.str();

}

StrTy DumpIndex(const IdxTy flags)
{
Ss ss;
for(IdxTy i=0; i<m_ttc.index_size(); ++i)
{
const auto & x=m_ttc(i);
const IdxTy k= m_ttc.entries(  x.m_lp,  x.m_mp,  x.m_gh );
ss<<MMPR(k)<<" "<<m_ttc(i).dump();
ss<<CRLF;
}
return ss.str();

} // DumpIndex



public: // ASSFUDD 
St m_st;

private:
TabTermCache m_ttc;

}; // mjm_tabulated_functions
 
////////////////////////////////////////////
#ifdef  TEST_MJM_TABULATED_FUNCTIONS
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

class tester {
typedef mjm_cli_ui<tester> Cli;
public:
 void cli_cmd( Cli::list_type & choices,  const char * frag)
{
/*const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
*/
}

 void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

 }; // tester
typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_TABULATED_FUNCTIONS "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_tabulated_functions<Tr>  Myt;
//Myt x(argc,args);
Myt x;
mjm_ragged_table r;
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
if (cmd=="load") {
const IdxTy flag=atoi(cip.p2.c_str());
 r.load(cip.p1,flag!=0 ); x.load(r,flag);  } 
if (cmd=="list") {
const IdxTy lp=atoi(cip.p1.c_str());
const IdxTy mp=atoi(cip.p2.c_str());
const IdxTy gh=atoi(cip.wif(3).c_str());
x.list_terms(lp,mp,gh);
} 
if (cmd=="index") {
MM_ERR(x.dump_index())
}
if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_LEGENDRE_STUFF_H__ 
