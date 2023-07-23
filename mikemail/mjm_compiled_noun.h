#ifndef MJM_COMPILED_NOUN_H__
#define MJM_COMPILED_NOUN_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_loo_parsing.h"
//#include "mjm_ragged_forms.h"
//#include "mjm_snacks_accumulator.h"
//#include "mjm_constrained_file.h"
#include "mjm_unit_crap.h"
#include "mjm_collections.h"
#include "mjm_misc_parse.h"
// faster than calling date in bash doh 
#include "mjm_calendar.h"
#include "mjm_diary_parse_error.h"
#include "mjm_string_tokenizer.h"

#include "mjm_muqed_phrase.h"



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

/*

2020-12-30 or so enum
2021-01-02 allow rational or fraction input  of form N/M with or
without units but not pxq where p or q are fractions.  


*/


//template <class Tr> class mjm_compiled_noun_collection<Tr>;
template <class Tr>
class mjm_compiled_noun //   : public mjm_ragged_forms<Tr> 
{
 typedef mjm_compiled_noun Myt;
//typedef  mjm_ragged_forms<Tr>  Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::OsTy Os;
 typedef typename Tr::Ofs Ofs;
// why the fck wont' these inherit fck 
//typedef typename Super::Ragged Ragged;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line;
typedef StrTy Word;
typedef std::vector<Word> Words;
typedef std::map<StrTy, IdxTy> DictTy;
typedef mjm_var_parse<Tr> CharClass;
//typedef typename Super::CharClass CharClass; 
//typedef mjm_constrained_file<Tr> BackingFile;
//ZZZZBackingFile m_file;
typedef typename CharClass::read_buffer Rb;

typedef mjm_unit_crap<Tr> Units;
typedef typename Units::dimmed_value_type dimmed_value_type;
typedef typename Units::conversion_type conversion_type;
typedef typename Units::qty_dim Qdim;
typedef std::map<StrTy, Qdim> SfxMap;
// eventually derive from Units 
typedef std::vector<StrTy> LocalFactors;


typedef mjm_loo_parsing<Tr> Loo;
typedef std::map<StrTy, int > ModMap;
typedef std::map<StrTy,Qdim> PieceMap;
typedef std::map<StrTy,D> AdjMap;

// these need to be tokenized
// the value is meaningless these should be sets.. 
typedef std::map<StrTy,int> GroupMap;
typedef std::map<StrTy,int> SelectMap;
//typedef std::map<StrTy,StrTy> GroupMap;
typedef std::map<StrTy,int> ConflictMap;
//typedef std::map<StrTy,StrTy> ConflictMap;

typedef std::vector< Line> LineParts;
typedef std::map<StrTy,PieceMap> PieceMaps;
typedef std::map<StrTy,LineParts> LineMaps;


typedef std::map<StrTy,Line> IndustryCodes;
typedef std::map<StrTy,Line> MiscMap;

typedef std::map<StrTy,StrTy> EnumMap;
typedef mjm_muqed_phrase<Tr> Phrase;
typedef typename Phrase::file_finder_type Finder;


typedef std::map<StrTy,StrTy> PropMap;



class noun_evaluation
{
typedef noun_evaluation Myt;
typedef mjm_compiled_noun Cn;
//typedef parse_error Error;
typedef mjm_diary_parse_error<Tr> Error;
typedef std::vector<Error> Errors;
public:
typedef  PropMap properties_type;
//noun_evaluation():m_val(-1),m_of(-1) {}
class eval_unit;
noun_evaluation():m_eflags(0), m_cn(0),m_convert_status(0),m_rqty(-999),m_failed(false) {}
noun_evaluation(const StrTy & a ):m_eflags(0),m_cn(0),m_convert_status(0),m_alias(a),m_rqty(-999) ,m_failed(false){}
noun_evaluation(const StrTy & a,const StrTy & u ):m_eflags(0),m_cn(0),m_convert_status(0),m_alias(a),
m_runits(u) , m_rqty(-999) ,m_failed(false){}
noun_evaluation(Cn * cn ):m_eflags(0),m_cn(cn),m_convert_status(0),m_alias(cn->canon()),
m_runits(cn->runits()) , m_rqty(-999) ,m_failed(false){}

// called during explosion expansion. 
// evaluates to the amount of cn implied by evaluation eup 
// using the contents amount qd.
noun_evaluation(Cn * cn,Myt * eup, const Qdim & qd, Units * pu  ):
m_eflags(eup->m_eflags),
m_ate(eup->m_ate),m_offered(eup->m_offered),
m_cn(cn),m_convert_status(0&&(eup->m_convert_status)),m_alias(cn->canon()),
m_runits(cn->runits()) ,m_rqty(-999),m_failed(false)
{
Scale(cn,eup,qd,pu);
} 


void noun( Cn * cn ) { m_cn=cn; } 
Cn & noun() {

if (m_cn==0){
MM_MSG(" request for null noun() in noun_evaluation doh ")
MM_ERR(" request for null noun() in noun_evaluation doh ")
}
 return *m_cn; }
const Cn & noun()const  {
if (m_cn==0){
MM_MSG(" request for const null noun() in noun_evaluation doh ")
MM_ERR(" request for const null noun() in noun_evaluation doh ")
}
 return *m_cn; }
class eval_unit
{
typedef eval_unit Myt;
public:
enum { A=0, B=1, SFX=2};
eval_unit(): m_a(1),m_b(1),m_q(1),m_state(0),m_sfx() {}
void clear() { m_a=1; m_b=1; m_state=9; m_sfx=""; } 
const StrTy & canon() const { return m_alias; } 
void count(const Qdim& ns) {
MM_ERR("count in eval "<<MMPR3(have_sfx(),ns.qty(),ns.sfx()))
n(ns.qty()); sfx(ns.sfx());   }
Myt  operator*(const D & p) { Myt x=(*this); x.m_q*=p; return x; } 
bool has_qty_data() const
{// zero specified is ok, if not zero needs to be by default 
if (m_q==0) return false;
if (Bit(m_state,A)) return true;
return Bit(m_state,B); 
}
// this is a dumb place for updating m_q... 
void a(const D & _a) { m_a=_a; m_q=m_a*m_b; m_state|=(1<<A);}
void b(const D & _b) { m_b=_b; m_q=m_a*m_b; m_state|=(1<<B);}
bool n(const D & _n) { 
if (Bit(m_state,A)) 
{
if (!Bit(m_state,B)) b(_n);
else return false;
}
 else { a(_n);}  return true;  }
bool n(const char * s) {return n(atof(s)); }  
void sfx(const StrTy & _s) { m_sfx=_s; m_state|=(1<<SFX);}
void sfx(const char *  _s) { m_sfx=StrTy(_s); m_state|=(1<<SFX);}
bool have_sfx() const { return Bit(m_state,SFX); } 
const StrTy & sfx()const {return  m_sfx;}
const D& qty() const { return m_q; }
const D& hi() const { return m_q; }
const D& lo() const { return m_q; }
void set_pct(const Myt & x) { m_sfx=x.m_sfx; m_q=x.m_q*m_q*.01; } 
StrTy dump() const { Ss ss; return dump(ss).str(); }
Ss &   dump(Ss & ss) const { 
ss<<MMPR(Bit(m_state,SFX)); ss<<MMPR(Bit(m_state,A)); 
ss<<MMPR2(Bit(m_state,B),m_q); ss<<MMPR4(m_a,m_b,m_state,m_sfx);
return ss; } 

D m_a,m_b,m_q;
IdxTy m_state;
StrTy m_sfx;
friend noun_evaluation;
}; // eval_unit
// so the noun_evaluation creates maybe 2 numbers of eaten and offered
// which may be adapted to outcomes or other non-consumed entries
// Results are entered most directly in terms of how they were measured
// usually volume kitchen units. Reporting units usually metric weights.
// A units conversion table converts among fixed generic units
// that can be metric, english or local quirky units.
// Private units specific to the noun may be provided too.  
// Density or iu etc are commonly expected. 
// Conversion may not succeed for many reasons though. 
// TODO check API but conversion success may be validated
// during markup but in any case reporting has to handle it. 
// back to noun_evaluation
const StrTy & alias() const { return m_alias; } 
bool ok() const { return m_errors.size()==0; }
enum EVAL_FLAGS  {STRING=0  }; 
const IdxTy flags() const { return m_eflags; }
void is_str(const bool x) {  Bit(m_eflags,STRING,x); } 
const bool  is_str() { return Bit(m_eflags,STRING); } 
// non-const as it may be computed later
const StrTy& string_value() { return StringValue(); }
void  string_value(const StrTy& s ) {m_string_value=s; is_str(true); }
const D & qty() const { return m_ate.qty();}
const StrTy & sfx() const { return m_ate.sfx();}
const D & of() const { return m_offered.qty();}
const StrTy & of_sfx() const { return m_offered.sfx();}
const bool  has_of() const { return m_offered.has_qty_data();}
const D pct() const { if (of()!=0) return qty()/of()*100.0; return -1;  } 
// the noun evaluation could be vector valued or at least a pair
// of offered and eaten. Conversion goes with each number doh 
// this has to account for failure and missing and "of" too
enum { RUSEFUL=0 }; 
const bool ruseful() const { return Bit(m_convert_status,RUSEFUL); } 
//const D rqty() const {if (m_rqty<0) return qty();  return m_rqty;}
const D rqty() const {if (ruseful()) return m_rqty;  return qty();}
//const StrTy & rsfx() const {if (m_rqty<0) return sfx(); return m_runits;}
const StrTy & rsfx() const {if (ruseful()) return m_runits; return sfx(); }
const StrTy & rrsfx() const { return m_runits;}
bool failed() const { return m_failed; } 
// this is easier than a ref and gets a failure easier... 
const StrTy * property(const StrTy & k ) const { return Property(k); } 
bool  property(const StrTy & k,const StrTy & v ) { return Property(k,v); } 
bool set(StrTy & d,const StrTy & s ) { return Set(d,s); } 
// crap
IdxTy calculate(Units & u, const LocalFactors & lcl )
{ return Calculate( u, lcl ); } 
//const D & accepted() const {if (m_of<=0) return 1;  return m_val/m_of;}
const StrTy & canon() const { return m_alias; } 
//const StrTy & units() const { return m_units; } 
Errors & errors() { return m_errors; } 
void add_error(const StrTy & w, const IdxTy pos, Ss & ss)
{ add_error(w,pos,ss.str().c_str()); } 
void add_error(const StrTy & w, const IdxTy pos, const char * msg)
{ 
MM_ERR(MMPR3(w,pos,msg))
m_errors.push_back(Error(w,pos,msg)); }
StrTy dump() const { Ss ss; return dump(ss).str(); }
Ss& dump(Ss & ss) const 
{
ss<<MMPR(m_alias);
ss<<"m_ate "; m_ate.dump(ss);
ss<<"m_offered "; m_offered.dump(ss);
MM_LOOP(ii,m_errors) { (*ii).dump(ss); } 
return ss; 
}

bool has_qty_data() const 
{ return m_ate.has_qty_data()||m_offered.has_qty_data(); } 
Myt operator*(const D eta) const 
{ Myt x=(*this); x.m_ate.m_q*=eta;
 x.m_offered.m_q*=eta;
return x;
}

// AFCK
//template <class Tos> static  Tos& operator<<(Tos & os, const Qdim & q) { os<<q.dump(); return os; }



void scale(const D offered, const D ate )
{
// this is operating on a macro def, not right
D frac=(offered>0)?(ate/offered):(-1);
//if (frac>=0) m_offered.m_q*=m_ate.m_q/frac;
// need to tell it it has real data not force m_q 
m_ate.m_q*=ate;
if (frac>0) m_offered.n(m_ate.m_q/frac);
//MM_ERR(" scaling "<<MMPR4(offered,ate,frac,m_offered.m_q))
if (m_rqty>0) { m_rqty=m_rqty*ate; } // m_rqty>0
} // scale 
private:

typedef typename Units::power_key_type Pk;
typedef typename Qdim::powers_map_type Qpow;
typedef typename Qdim::ratio_parse Qratio;
StrTy Nuid(const Cn & noun) const { return noun.name(); } 
// the noun density may have been over ridden during evalu 
// t1 is parent serving info, t2 is the evaluation amount/unit

// make a noun evaluation for cn that is a component of 
// the noun evaluated for eup and scale by its quantity
// with qd as the component amount. 
void Scale(Cn * cn,Myt * eup, const Qdim & qd, Units * pu  )
{
// so, the reporting output for the parent is not here AFAICT
// doh... 
Myt & eu= *eup;
// amount of this noun  per serving of this noun  
Qdim f1=cn->serving();
// if the serving size is zero, just report in qd units
// can escape here but done later for now... 
// serving size of parent
Qdim f2=eup->noun().serving();
//const StrTy * densityp= eup->property("density");
//const StrTy density=(densityp?(*densityp):StrTy());
// this does not have "thing" but implies parent/parent
const StrTy density=eup->noun().density();
MM_ERR(" ebval amount  "<<MMPR3(eu.noun().canon(),eu.qty(),eu.sfx()))
MM_ERR(" parent  serving size "<<MMPR3(eu.noun().canon(),f2.qty(),f2.sfx()))
MM_ERR(" component  serving size "<<MMPR3(cn->canon(),f1.qty(),f1.sfx()))
// 2021-09-09 escape with no component serving size
//if (f1.qty()==0)
{
// this is returning the originally entered amount for the parent.. 
// because the reporting amount is not this doh... 
// go ahead and convert cleanly... 
// eu.sfx() and eq.qty() are the  parent amounts
// want eu*qd/f2 in reporting units for cn;
MM_ERR(" FCK "<<MMPR3(qd.dump(),f2.dump(),density))
D sv;
//IdxTy purc1=pu->convert(sv,f2.qty(),eu.sfx(),f2.sfx(),density);
// instead, convert eu to same units as serving f2, 
IdxTy purc1=pu->convert(sv,eu.qty(),f2.sfx(),eu.sfx(),density);
D nsvngs=sv/f2.qty();
MM_ERR(MMPR3(purc1,sv,nsvngs))
// now convert to cn reporting of serving units if they exist

eval_unit x;
x.sfx(eu.sfx());
x.n(1);
x.n(eu.qty());
eval_unit y;
y.sfx(qd.sfx());
y.n(nsvngs);
y.n(qd.qty());
m_ate=y; // I guess we are the object "cn"? 
MM_ERR(" escaping  same? " << MMPR3(y.dump(),this,cn))

return ; 
} // no serving 
} // Scale

bool Debugg() const { return false; }

// calculate(
/*
Given "x_v" of units  "x_u" 
and a set of conversion factors a_{i,j},
find "y_v" in units of "y_u" 
The units cluster into subgraphs ( length, volume, mass, etc)
connected by object specific factors. A conversion may
exist between two objects as with contents .  

// create "r" according to conversion specified by ct
// using stored units and anything in vc that is needed
// return requested value of anything possible
IdxTy Convert(dimmed_value_type & r, const conversion_type & ct
    , const Vc & vc, const IdxTy flags) const
{

_conversion_entry() : m_factor(0) {}
_conversion_entry( const StrTy &  from_unit, const StrTy & from_thing,
const StrTy & to_unit, const StrTy & to_thing, const F & f)
:m_from_unit(from_unit), m_from_thing(from_thing),
m_to_unit(to_unit), m_to_thing(to_thing),m_factor(f) {}
_conversion_entry( const StrTy &  from_unit,
const StrTy & to_unit,  const F & f)
:m_from_unit(from_unit),
m_to_unit(to_unit), m_factor(f) {}


*/
// Factors are stringgs with units for conversion 
// include density or serving 

// the noun_evaluation consists of a direct parsing
// of up to two factors and a unit. These can be arbitrary
// and noun specific. The noun or other specification
// has a reporting unit that is derived from the fixed
// units and noun-specific of local factors.
// Impl is now mostly hidden although currently
// ambiguous handling of likely imperfect cases. 
IdxTy Calculate(Units & u, const LocalFactors & lcl )
{
// we are in a noun_evaluztion, 
if (Debugg()) 
{ MM_ERR(" calculaet1 "<<MMPR(qty())<<MMPR3(m_rqty,rrsfx(),sfx())) } 
//D x=0;
MM_ERR(MMPR(alias()));
//MM_LOOP(ii,lcl) { MM_ERR(" Local Factors "<<MMPR((*ii).dump())) }
MM_LOOP(ii,lcl) { MM_ERR(" Local Factors "<<MMPR((*ii))) }


//if (rrsfx.length()==0) { m_rqty=qty(); retrun 0; } 
if (rrsfx().length()==0) { return 0; } 
dimmed_value_type result;
const StrTy dname="";
StrTy in_unit=sfx();
// this is in servings if that exists else a "count"
// if there is no serving and the qty has no units 
// it is just taken as a count or factor for the expansion. 
if (in_unit.length()==0){ in_unit="serving";    } 
//conversion_type ct(sfx(),dname,rrsfx(),dname,qty());
// the specified input unit is either the suffix entered
// by the user or "serving" . The "thing" values are 
// left blank for now although during an explosion
// they are important.  
// the noun alias and local factor "things" need to be
// considered 
// this should have a ctor with order reversed led by qty. 
conversion_type ct(in_unit,dname,rrsfx(),dname,qty());
//IdxTy rc=0; // u.convert(x, qty(),rrsfx(),sfx(),rho);
// in the case of a blank "to" unit the input unit should be adopted
// although its not clear that happens... 
IdxTy rc= u.convert(result, ct,lcl,0);
// TODO FIXME if the units differ the reporting unit needs to be fixed doh
const bool all_ok=((rc==0)||(rc==1))&&(result.u()==ct.to_unit());
const bool catastrophe=(qty()!=0)&&(result.v()==0);
if (Debugg()||!all_ok||catastrophe) 
{ Ss ss; 
ss<<" convert failed? "<<MM_MARK<<" "<<MMPR4(rc,result.dump(),ct.dump(),m_rqty);
ss<<MMPR4(qty(),rsfx(),rrsfx(),sfx())<<MMPR(catastrophe); 
ss<<MMPR(Loo::Dump(lcl));
//if (!all_ok) throw std::logic_error(ss.str());
// this should not happen... 
//if (all_ok&&catastrophe) throw std::logic_error(ss.str());
MM_ERR(ss.str());
//ss<<" convert failed "<<MMPR4(rc,m_rqty,qty(),rsfx())<<MMPR3(rrsfx(),sfx(),rho);
if (!all_ok||catastrophe) add_error(alias(),0,ss); 
}
if (all_ok&&!catastrophe) {   m_rqty=result.v(); m_convert_status|=(1<<RUSEFUL); } 
return rc;
} // Calculate 

const StrTy * Property(const StrTy & k ) const 
{ 
const auto & ii = m_props.find(k);
if (ii!=m_props.end()) return &(*ii).second;
return 0;
} 

const bool Property(const StrTy & k,const StrTy & v )  
{ 
const auto & ii = m_props.find(k);
const bool hit= (ii!=m_props.end());
 m_props[k]=v; 
return hit;
} 

const bool Set( StrTy & d,const StrTy & s )  
{ 
const auto & ii = m_props.find(s);
const bool hit= (ii!=m_props.end());
if (hit)d=(*ii).second;
return hit;
} 

const StrTy & StringValue() { return m_string_value; }

Errors m_errors;
IdxTy m_eflags;
StrTy m_string_value;
eval_unit m_ate,m_offered;
//D m_val; D m_of; 
Cn * m_cn;
IdxTy m_convert_status;
StrTy m_alias; 
StrTy m_runits;
D m_rqty;
bool m_failed;
PropMap m_props;
//StrTy m_units;
friend mjm_compiled_noun;
}; // noun_evaluation 


enum { SEL_GROUP=1,SEL_ADJ=2,SEL_MISC=4, SEL_HUGE=1<<16}; 

class _noun_score 
{
typedef _noun_score Myt;
//class mjm_compiled_noun<Tr>;
typedef mjm_compiled_noun<Tr> Cn;
public:
//nsp[in]=NounScore(ps.admap,ps.dates, &xvec[i]);
_noun_score(): m_score(~0) {}
// the adjectives are organized in phrase, date is in d
// and noun to analyze is pointed to by np 
// quantifiers for eaten and of , override modifiers
// and selectors.
// mfg should select, but vs fresh/frozen? etc 
_noun_score(const Phrase & p, const StrTy & d, Cn * np)
{ Score(p,d,np); }

_noun_score(const ModMap & mm, const StrTy & d, Cn * np)
{
Score(mm,d,np);

} 
IdxTy score() const { return m_score; } 
private:

// for a collection of selecting values in c, if the propery map p
// contains an entry for key and that entry is in c return v; 
template <class Tm,class Tc > 
int LutScore(Tm & p, const StrTy & key, const Tc & c, const int v)
{
const auto ii=p.find(key);
if (ii==p.end()) return 0;
// note that "c" must not contain null or invalid keys. 
const auto jj=c.find((*ii).second);
if (jj==c.end()) return 0;
return v;
} // LutScore

void Score(const Phrase & p, const StrTy & d, Cn * np)
{
m_score=0;
const bool date_ok=(*np).date_ok(d);
if (!date_ok) return; 
const bool looks_ok=(*np).has_right_modifiers(p.admap());
if (!looks_ok) return; 
m_score+=8;

// the value here is a position but the score
// is a bit array. 
auto sm=p.admap();
Cn::select_map_type smm;
MM_LOOP(ii,sm) smm[(*ii).first]=0; 
// dates - no definitive but close 
np->scores(smm);
IdxTy s=m_score;
MM_LOOP(ii,smm) s+=(*ii).second; 
m_score=s;

const auto & prop=p.prop();
// m_upc hopefully is used with number as the key with value the code system(s)
// wth is access allowed to m_upc??? 
int s1= LutScore(prop, "UPC", np->m_upc, SEL_HUGE);
m_score+=s1;

// select words from mfg, groups, maybe modifiers 
// ideally user should make a different noun for each
// combination but things come up . 
// specific values like UPC codes definitive but maybe flag
// if dates bad etc. 

// properties - can override or select. 

// ignore others like density which override not select


//MM_DIE(" no code")
//Score(p.admap(),d,np); 

} // Score
void Score(const ModMap & mm, const StrTy & d,  Cn * np)
{

const bool looks_ok=(*np).has_right_modifiers(mm);
const bool date_ok=(*np).date_ok(d);
m_score=(looks_ok?1:0)+(date_ok?2:0);
const bool slow=true;
if (slow)  m_have=mm;
} // Score 

StrTy m_date,m_startdate,m_enddate;
typedef std::map<StrTy,int> Map;
Map m_have,m_need,m_not;
IdxTy m_ds,m_score;

}; // _noun_score



// note this is a BAD number NOT a bad word 
enum {BAD=~0 }; 
public:
mjm_compiled_noun() {Init(Line()); }
mjm_compiled_noun(const Line & l) {Init(l); }
mjm_compiled_noun(const Line & l,const StrTy & name, const IdxTy f ) 
{Init(l,name,f); }
enum NOUN_FLAGS { IGNORE=0,SYNONYM=1, DUMMY=2, NEED_UNITS=3, DMEL=4, WIN_TIE=5, EXP_EVAL=6 ,SHRIMP_OK=7};
enum NOUN_VALUE_FLAGS {FREE=0, RESOURCE=1,ENUM=2,VALUE_IS_VALUE=3 };
typedef _noun_score noun_score;
typedef noun_evaluation evaluation;
typedef typename evaluation::eval_unit Eu;
typedef LineParts line_parts_type;
typedef PieceMaps parts_type;
typedef PieceMap partss_type;

const PieceMaps &   pieces() const { return m_alt_parts; } 

typedef SelectMap select_map_type;
const select_map_type  & selectors() const { return m_selects; }
void scores(SelectMap & sm) const  { Scores(sm); } 
bool ignore() const { return  (Bit(m_flags,IGNORE)); } 
bool def() const { return  (Bit(m_flags,WIN_TIE)); } 
bool shrimp_ok() const { return  (Bit(m_flags,SHRIMP_OK)); } 
// TODO these can both be parsed during init now doh 
const StrTy & density() const { return m_density;}
const StrTy & runits() const { return m_report_units;}
const StrTy & name() const { return m_name; } 
const StrTy & canon() const 
{if (m_canon.length()!=0)  return m_canon; return m_name; }
bool has_right_modifiers(const ModMap & m) { return HasRightModifiers(m); } 
bool date_ok(const StrTy & d) {
if (m_enddate!="") if (d>m_enddate) return false;  
return d>=m_startdate; } 

evaluation evaluate(const StrTy & name, Phrase & phrase,  CharClass * pcc, Units & units)
{ return Evaluate(name, phrase,pcc,units); } 

bool is_private_sfx(const char *  a) const{return is_private_sfx(StrTy(a)); } 
bool is_private_sfx(const StrTy & a)const 
{ return (m_private_sfxs.find(a)!=m_private_sfxs.end()); } 
const StrTy & rconv(const StrTy & x) const {  return Rconv(x); } 
bool expand() const { return Bit(m_expand,0); }  // 0
bool expand_lines() const { return Bit(m_expand,3); }  // 3
IdxTy expand_rank() const { return (m_expand>>1)&3; } //1,2 
PieceMap & parts() { return  m_alt_parts[""]; } 
const PieceMap & parts()const  { return  m_alt_parts[""]; } 
PieceMap & parts(const StrTy & nm) { return  m_alt_parts[nm]; } 
PieceMap  parts(const StrTy & nm,const IdxTy fix)const  { return  Parts(nm,fix); } 

const PieceMap & parts(const StrTy & nm) const  { return  m_alt_parts[nm]; } 
LineParts & lparts() { return m_alt_lines[StrTy()]; } 
const LineParts & lparts()const  { return m_alt_lines[StrTy()]; } 
LineParts & lparts(const StrTy & nm ) { return m_alt_lines[nm]; } 
const LineParts & lparts(const StrTy & nm )const  { return m_alt_lines[nm]; } 

bool have_parts(const StrTy & nm) { return  m_alt_parts.find(nm)!=m_alt_parts.end(); } 
bool have_line_parts(const StrTy & nm) { return  m_alt_lines.find(nm)!=m_alt_lines.end(); } 

//const PieceMap & parts(const StrTy & nm) const  { return  m_alt_parts[nm]; } 
const Qdim & serving() const { return m_default; } 
const GroupMap & groups() const { return m_groups ; } 
void dump(Os & os, const IdxTy flags)const  { Dump(os,flags); } 
private:

void Scores(SelectMap & sm) const 
{ MM_LOOP(ii,sm) { auto jj=m_selects.find((*ii).first);
if (jj!=m_selects.end()) (*ii).second=(*jj).second;  }
 }


const StrTy & Rconv(const StrTy & sfx) const { 
const auto ii=m_private_sfxs.find(sfx);
// only thing possible 
if (ii==m_private_sfxs.end()) return density();
return (*ii).second;
} // Rconv 
static bool Bit(const IdxTy f, const IdxTy b)   { return  ((f>>b)&1)!=0; }
static void  Bit(IdxTy& f, const IdxTy b,const bool x)   
{ if (x) f|=(1<<b); else f&= ~(1<<b);   }

// probably strip all to non-digit/decimal place from right
// to get sfx etc. Need catagorires like manufacturer or other
// matching adjective TODO FIXME 
bool HasRightModifiers(const ModMap & m) { 

MM_LOOP(ii,m_mods) // HasRightModifiers
{
//const StrTy & a=(*ii).first;

// FIXME this logic no good 
const bool must_have=((*ii).second>0);
const bool have=(m.find((*ii).first)!=m.end());
const bool violation=(must_have!=have);
if ( violation ) return false;
} // ii 
MM_LOOP(ii,m_nots)
{
if (m_nots.find((*ii).first)!=m_nots.end()) return false;
}
return  true; 

} 



// 0 used, 1 unused, -1 bad
IdxTy  EvaluateAdj(Eu & eu, evaluation & rc,  const StrTy & a, CharClass * pcc, const IdxTy & i,Units & units  )
{
IdxTy irc=0;
const StrTy & name=m_name;
 Rb rb;
rb.clear();
(*pcc).parse_groups(rb,a.c_str());
const IdxTy n=rb.string_count();
//MM_ERR(" parsing result "<<MMPR2(a,n))
if (n==0) return 1; // continue;
if (n==2)
{
if ((rb[0][0]==CharClass::TINT)||(rb[0][0]==CharClass::TFLOAT))
{
if (!eu.n(rb[1])) {  rc.add_error(a,i,"too many numbers"); irc=BAD;  } 
} // int or float
else { Ss ss; ss<<"bad adjective for "<<MMPR(name);   
	rc.add_error(a,i,ss.str().c_str()); irc=BAD;irc=~0; } 
return irc; // continue;
} // n==2
if (n==4) // number sfx
{
const bool one_no= ((rb[0][0]==CharClass::TINT)||(rb[0][0]==CharClass::TFLOAT));
const bool two_sfx= ((rb[2][0]==CharClass::TALPHA));
//MM_ERR(" parsing result "<<MMPR4(a,n,one_no,two_sfx))

if ((one_no)&&(two_sfx))
{	
if (!eu.n(rb[1])) {irc=BAD;  rc.add_error(a,i,"too many numbers");  } 
// drop pct code, include with sfx
const bool is_pct=(strcmp(rb[3],"pct")==0) ;
const bool is_unit=units.have(rb[3])||is_private_sfx(rb[3]); // m_private
if (is_unit||is_pct)
{
if (eu.have_sfx()){  
irc=BAD;
Ss ss; ss<<" multiple sfx or pct  "<<MMPR2(a,eu.sfx()); 
rc.add_error(a,i,ss);  } 
//MM_ERR(" adding suffix "<<MMPR(rb[3]))
eu.sfx(rb[3]);
} // units  and pct 
else{
irc=BAD;
Ss ss;
ss<<MM_MARK<< " unknown suffix "<<rb[3]<<" from "<<a;
 rc.add_error(a,i, ss.str().c_str());  }

return irc; // continue;
} // n+ sfc 
else { irc=BAD; rc.add_error(a,i,"bad fields ");  } 
} // n==4

// TODO 0,5 is a typo and may have utility later gets mess up here though 
// this xhould also support n-x-m units dettachd
else if (n==6) //  added "/" for rats, "x", shrimp only convert to mg 
{
const bool one_no=((rb[0][0]==CharClass::TINT)||(rb[0][0]==CharClass::TFLOAT));
const bool three_no= ((rb[4][0]==CharClass::TINT));
const bool isx=(strcmp(rb[3],"x")==0);
const bool israt=(strcmp(rb[3],"/")==0);
const bool iscomma=(strcmp(rb[3],",")==0);
const bool isshrimp=(one_no&&three_no&&isx);
if ((one_no&&three_no&&israt))
{
// this should be ints but be flexible as hex not allowed loll 
const D n=atof(rb[1]);
const D d=atof(rb[5]);
if (d==0) rc.add_error(a,i,"zerp denominator in rational input ");
const D v=n/d;
if (!eu.n(v)){
irc=BAD;
Ss ss; ss<<" toomany numbers "<<MMPR4(rb[1],rb[3],rb[5],v);
rc.add_error(a,i,ss);  
return irc;
}
//if (!eu.n(rb[5])) rc.add_error(a,i,"rb[5] too many numbers");
//if (eu.have_sfx()){  rc.add_error(a,i," n8 multiple sfx ");  } 
//eu.sfx(rb[7]);
return irc;
} // israt
if (!isx)
{
irc=BAD;
Ss ss; ss<<" syntax not supported "<<MMPR(iscomma);
rc.add_error(a,i,ss);  
return irc;

}
//if (!(one_no&&three_no&&isx))
if (!(isshrimp))
{ irc=BAD;  rc.add_error(a,i," bad alphanum adj"); } 
if ((isshrimp&&!shrimp_ok()))
{ irc=BAD;  rc.add_error(a,i," attach sfx to form MxN to avoid shrimp syntax "); } 


const D q=atof(rb[1]);
const IdxTy sz1=atoi(rb[5])/100;
// this can seg fault... 
const IdxTy lns=strlen(rb[5]);
if (lns<4)
{
irc=BAD;
Ss ss; ss<<" this looks like shrimp but not ribht and no other impl ";
rc.add_error(a,i,ss);  
return irc;
}
const IdxTy sz2=(lns>2)?atoi(rb[5]+2):0;
//MM_ERR(MMPR2(sz1,sz2))
//MM_FAULT
D avg=D(sz1+sz2)/2;
D oz=D(q)/avg;
// this returns the conversion factor, 
const IdxTy rcr=units.to(avg,"g","pounds");
oz=oz*avg;
if (rcr!=0) { irc=BAD;
Ss ss; ss<< " missing units for shrump oz to grams "<<MMPR2(rcr,oz);
 rc.add_error(a,i,ss) ;} 
//avg=avg*oz;
//avg=D(q)/avg*28.0;
eu.sfx("g");
if (!eu.n(oz)){ irc=BAD;  rc.add_error(a,i,"too many numbers"); } 
return irc; // continue;
} // n==6 

else if (n==8) // number x number units, units must be attached  
//else if (n==8) // also want  number / number units,  
{

const bool one_no= ((rb[0][0]==CharClass::TINT)||(rb[0][0]==CharClass::TFLOAT));
const bool two_x= ((rb[3][0]=='x'))||((rb[3][1]==0));
const bool two_rat= ((rb[3][0]=='/'))||((rb[3][1]==0));
const bool three_no= ((rb[4][0]==CharClass::TINT)||(rb[4][0]==CharClass::TFLOAT));
const bool ok=one_no&&two_x&&three_no;
const bool rat=one_no&&two_rat&&three_no;
if (rat)
{
const D n=atof(rb[1]);
const D d=atof(rb[5]);
const D v=n/d;
if ((v==0)|| !eu.n(v)){
irc=BAD;
Ss ss; ss<<" too many or bad numbers "<<MMPR4(rb[1],rb[3],rb[5],v);
rc.add_error(a,i,ss);  
return irc;
}
if (eu.have_sfx()){  rc.add_error(a,i," n8 rat multiple sfx ");  } 
eu.sfx(rb[7]);

return irc;
} // rat 
if (!ok) rc.add_error(a,i,"unknown adjective sequence");
else { 
if (!eu.n(rb[1])) rc.add_error(a,i,"rb[1] too many numbers");
if (!eu.n(rb[5])) rc.add_error(a,i,"rb[5] too many numbers");
if (eu.have_sfx()){  rc.add_error(a,i," n8 multiple sfx ");  } 
eu.sfx(rb[7]);
}
return 0; // continue;

} // 8
// TODO variations of rational inputs with "X" and units

irc=BAD;
Ss ss; ss<<" bad adjective too many fields "<<MMPR4(a,i,n,rb[n-1]);
rc.add_error(a,i,ss);
return irc;
} // EvaluateAdj  6 param no pct

// 0 used, 1 unused, -1 bad

#if 0 
#endif

// alt
// find "of" or pct, split
// remove adj used for calculations
void EvaluateOne(evaluation & rc,Eu & eu, const StrTy & name,  Words & adj, CharClass * pcc, Units & units) // const 
{
Words unused;
Word a="",last="";
MM_SZ_LOOP(i,adj,sz)  { 
last=a;
a=adj[i];
const char * ap=a.c_str();
//if (a=="of")
const bool is_pct=(a=="pct");
const bool private_sfx=is_private_sfx(a);
if (units.have(ap)||is_pct||private_sfx) { 
if (eu.have_sfx()){ 
Ss ss; ss<<" multiple suffixes "<<MMPR2(a,eu.sfx());
 rc.add_error(a,i,ss);  
} 
eu.sfx(a);
continue; } // sfx

// not important in this evaluation but may give
if (m_adjs.find(a)!= (m_adjs.end())) {unused.push_back(a);  continue; }

IdxTy  eac=EvaluateAdj(eu, rc, a,  pcc,i,units );
//if (eac==0) continue; 
if (eac==1) unused.push_back(a);

} // i
adj=unused; 
} // EvaluateOne
IdxTy EvaluateSpecial10(evaluation & rc, Words & a1, Words & a2, Phrase & phrase, CharClass * pcc, Units & units, const StrTy & name)
{
// insert as value into rc,
if(Bit(m_vflags,ENUM)) { 
auto ii=m_enum.find(a1[0]);
if (ii==m_enum.end()) 
{
Ss ss; 
ss<<a1[0]<< " not in enum    ";
MM_LOOP(ii,m_enum) { ss<<" "<<(*ii).second; } 
rc.add_error(a1[0],0,ss);  
return ~0; 
}
a1[0]=(*ii).second; 
// otherwise the mapped word will just be evaluated as normally. 
if(Bit(m_vflags,VALUE_IS_VALUE)) { rc.string_value(a1[0]); }
return 0; 
} 
if(Bit(m_vflags,VALUE_IS_VALUE)) { rc.string_value(a1[0]); }
if(Bit(m_vflags,FREE)) { a1.clear(); return 0; } 
Finder * pf= phrase.finder();
if(!Bit(m_vflags,RESOURCE)) { a1.clear(); return 0; } 
if (!pf) 
{
Ss ss; ss<<" no resource finder  "<<MMPR(name);
rc.add_error(a1[0],0,ss);  
return ~0;
}
const IdxTy cnt=(*pf).count(a1[0]);
if (cnt!=1)
{
Ss ss; ss<<"resource count( non-canonicalized )  not 1  "<<MMPR2(cnt,a1[0]);
rc.add_error(a1[0],0,ss);  
return ~0;
}
a1.clear();

return 0;

}// EvaluateSpecial10

IdxTy EvaluateSpecial(evaluation & rc, Words & a1, Words & a2, Phrase & phrase, CharClass * pcc, Units & units, const StrTy & name)
{
IdxTy rcv=0;
//const bool viv=Bit(m_vflags,VALUE_IS_VALUE);
const IdxTy sz1=a1.size();
const IdxTy sz2=a2.size();
if (sz1==1) if (sz2==0) 
{
// this MAY fall through and get a numerical value equal
// to the count or other rules 
const IdxTy evalrc= EvaluateSpecial10(rc,a1,a2,phrase,pcc,units,name); 
return evalrc;
} // one word 



return rcv;
} 
 evaluation Evaluate(const StrTy & name, Phrase & phrase, CharClass * pcc, Units & units)  
{
evaluation rc(canon(),runits());
rc.noun(this);
// TODO photo name and path and description
// TODO adjectives that translate into numeric ranges etc
// can still have a value except the numeric parsers need to ignore
// this stuff. Perhaps take it our from a1 first and
// then parse a1 and a2 normally if anything is left.
// these now come from the Phrase but can be depleted
// not sure this is good now 
//MM_ERR(" phrase evalulate")
Words a1=phrase.a1();
Words a2=phrase.a2();
// non-quant data first
// may want to report text as a suffix for enumerated events
// or not as for tihngs like file names. 
// paths are a problem for the latex compatibility. 
if (Bit(m_flags,EXP_EVAL))
{
// set some flags in rc based on phrases
const bool viv=Bit(m_vflags,VALUE_IS_VALUE);
IdxTy rces= EvaluateSpecial(rc, a1, a2, phrase, pcc, units,name);
// I guess to validate a resource the info should be in phrase?
// We don't know our container and maybe it can vary 
if (viv) if (!rces) return rc;

}// EXP_EVAL
// fall through although this will reset qty and sfx maybe not 

//{ Ss ss; MM_SZ_LOOP(i,a1,sz) { ss<<MMPR(a1[i]); } MM_ERR(ss.str()) }  
EvaluateOne(rc,rc.m_ate,name,a1,pcc,units);
if (a2.size()) EvaluateOne(rc,rc.m_offered,name,a2,pcc,units);
//MM_ERR(MMPR2(rc.m_ate.sfx(),rc.m_ate.qty()))
if (rc.m_ate.sfx()=="pct")
{ // get units etc from offered, 
rc.m_ate.set_pct(rc.m_offered);
}
// check for missing units then use the default or serving size
// which right now is just a numbers DIMENSIONLESS
if (rc.m_ate.sfx().length()==0)
{
// some nouns may require a unit...
if (Bit(m_flags,NEED_UNITS))
{
Ss rr;
MM_LOOP(zz,a1) { rr<<(*zz)<<" "; } 
Ss ss; 
ss<<" noun requires unit but none here  "<<MMPR3(rc.m_ate.qty(),rc.m_ate.sfx(),rr.str());
rc.add_error(name,0,ss); 
} // requre unit
if (m_report_units!="serving")
{
if (m_default.qty()!=0) {
// m_serving m_default
rc.m_ate.count(m_default);
if(a2.size()) rc.m_offered.count(m_default);
//MM_FAULT
}  // m_default
//else MM_ERR(" need default serving info doh "<<MMPR(m_name) )
else MM_ONCE(" VERIFY THIS IS OK need default serving info doh "<<MMPR(m_name), )
} // serving

} // sfx leng
// a1 and a2 contain unused adjectives
//if (rho.length()==0) rho=density(); 
StrTy rho=density(); 
rc.set(rho,"density");
LocalFactors lf; // this can be a noun member
if (m_density.length()) lf.push_back(m_density); 
// this already has been parsed
//lf.push_back(m_default+StrTy("/serving"));  // serving
if (m_default.qty()!=0) lf.push_back(m_default.dump()+StrTy("/serving"));  // serving
//MM_LOOP(ii,m_private_sfxs){lf.push_back((*ii).second+StrTy("/")+(*ii).first);} 
MM_LOOP(ii,m_private_sfxs){lf.push_back((*ii).second);} 
// finally calculate rc ate in reporting units. 
//rqty
//const IdxTy val=rc.calculate(units,rho); // find the reporting amount. 
const IdxTy val=rc.calculate(units,lf); // find the reporting amount. 
//MM_ERR(" fixing "<<MMPR4(m_name,rho,val,rc.qty())<<MMPR(rc.sfx()))

if (val!=0)
{
MM_ERR(" fixing "<<MMPR4(m_name,rho,val,rc.qty())<<MMPR(rc.sfx()))
Ss ss; 
ss<<" unable to fix units, maybe missing density "<<MMPR3(val,rc.m_ate.sfx(),density());
rc.add_error(name,0,ss); 
}
if (ignore())  // (Bit(m_flags,IGNORE))
{
if (rc.has_qty_data())
{ Ss ss; ss<<"ignore placeholder noun has meaninful entry  "<<canon(); 
 rc.add_error(m_name,0,ss); 
}
} // ignore
//MM_ERR(MMPR(rc.dump()))

return rc;
} // Evalute phrase



#if 0 
#endif

typedef std::map<StrTy, IdxTy> CmdMap;
enum 
{ PART=1,CANON=2,IGNOREN=3,DEFAULTQTY=4,MODIFIER=5
,NOT=6,ADJ=7,ADJECTIVE=8,xxCONTAINS=9,xxPIECEOF=10,TEXT=11
,GROUP=12,CONFLICT=13,STARTDATE=14,MFG=15,xxxQTY=16,UNITS=17,
DENSITY=18,URL=19,BIBTEX=20,SERVING=21, ENDDATE=22,RUNITS=23,
PRIVATEUNIT=24,EXPAND=25,BOMTEX=26, RPART=27,LPART=28,
//FLAGS=29,  
UPC=30, SETFLAGS=31, MISC=32,
ACTIVE=33, INGREDIENTS=34, SETVFLAGS=35,ENUMM=36
};
static CmdMap  Lut() 
{
CmdMap x;
x[StrTy("canon")]=CANON; // some systematic generic name 
x[StrTy("flags")]=SETFLAGS; // noun attribute flags 
x[StrTy("vflags")]=SETVFLAGS; // attributes of value if non-numeric/units 
x[StrTy("ignore")]=IGNOREN; // a placehold noun should have zero qty
x[StrTy("serving") ]=SERVING; // the noun has a default "size" 
x[StrTy("default-qty") ]=DEFAULTQTY;
x[StrTy("punit") ]=PRIVATEUNIT;
x[StrTy("units") ]=UNITS;
x[StrTy("runits") ]=RUNITS;
x[StrTy("density") ]=DENSITY;
x[StrTy("modifier")]=MODIFIER;
x[StrTy("adj")]=ADJ;
x[StrTy("adjective") ]=ADJECTIVE;
x[StrTy("not")]=NOT;
x[StrTy("startdate") ]=STARTDATE; // first date used for disambiguation 
x[StrTy("enddate") ]=ENDDATE; // last date this is valid 
x[StrTy("text") ]=TEXT; // human reable description 
x[StrTy("group") ]=GROUP; // user defined group memberships aminoacid, etc 
x[StrTy("conflict") ]=CONFLICT; // dietary design rule conclift list 
x[StrTy("UPC") ]=UPC; // disambiguation may still not be complete. 
x[StrTy("mfg") ]=MFG; // maker 
x[StrTy("url") ]=URL;
x[StrTy("bibtex") ]=BIBTEX;
x[StrTy("bomtex") ]=BOMTEX;
x[StrTy("misc") ]=MISC;
x[StrTy("active") ]=ACTIVE;
x[StrTy("part")]=PART; // various components into which the noun may expand
x[StrTy("rpart")]=RPART;
x[StrTy("lpart")]=LPART;
x[StrTy("ingredients") ]=INGREDIENTS;
x[StrTy("expand")]=EXPAND;
x[StrTy("enum")]=ENUMM;



return x;

}
static CmdMap & Llut() 
{
static CmdMap  map=Lut();
return map;
}
// fixing this means trashing the reference 
// some are all ranks, others may contain rank.
// likely all rank except for one or two actives
// need a constraint on rank amounts. Serving maybe
// should total to 1-active.   
PieceMap  Parts(const StrTy & nm,const IdxTy fix)  const
{ 
const IdxTy cner=expand_rank();
auto xx=m_alt_parts.find(nm);
if (xx==m_alt_parts.end()) return PieceMap();
PieceMap x=(*xx).second;
const IdxTy sz=x.size();
std::vector<IdxTy> ranks;
MM_LOOP(ii,x)
{
//const StrTy & wp=(*ii).first;
Qdim  qdp=(*ii).second;
const StrTy qsfx=qdp.sfx();
const bool is_rank=(qsfx=="rank");
if (is_rank) ranks.push_back(IdxTy(qdp.qty()));
} // ii 
if (ranks.size()==0) return x;
if (ranks.size()==sz) {}
// should remove the non-ranks but usually small or not imporant
D sum=0;
MM_LOOP(ii,x)
{
Qdim & y=(*ii).second;
Qdim z;
Rescale(z,y,cner,1);
sum+=z.qty();
} // ii 
if (sum==0) sum=1; // just do nothing to the scales...  
MM_LOOP(ii,x)
{
Qdim & y=(*ii).second;
Rescale(y,y,cner,(1.0/sum));
} // ii 

return  x;  
} 

void Rescale(Qdim & qout, const Qdim& qin,const IdxTy cner, const D df) const
{
switch (cner)
{
case 0 : {break; } // just use the rank as a qty, not intuitive though 
case 1 : {qout=Qdim(df,""); break; } // just count  
case 2 : {qout=Qdim(df/qin.qty(),""); break; } // more intuitive   
case 3 : {const bool cutoff=(qin.qty()>5); qout=Qdim(cutoff?0:df,""); break; } // cutoff   
//default: { MM_ERR(" bad expand case "<<MMPR3(w,wp,cn.expand_rank())) } 
default: { MM_ERR(" bad expand case "<<MMPR2(qin.sfx(),cner)) }
} // switch 

} // rescale 


static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); }
static int myatoi(const char * c) { return ::strtol(c,0,0); }

// template<class Tp >
//static void SetFlags (IdxTy & d, const Tp & m, IdxTy & i, const Line & l, const IdxTy sz)
#define NCASEFLAG(x,y,z) case x : { Loo::SetFlags(y,z,i,l,szl); break; } 

typedef std::map<StrTy,IdxTy> FlagMap;
//  non-quant values were added after the fact and kluged in
// these can be filenames to imges or supporting files other than bibtex
// or maybe free form text yet to be reported as other than "text" 
//enum NOUN_VALUE_FLAGS { };
static void BitVFlagMap(FlagMap & m)
{
// resource and not free means must exist 
// all files must have unique names regardless of dir 
m["free"]=1<<FREE; // free form human readable text  
m["file"]=1<<RESOURCE; // the value is a data file that needs to exist  
m["string"]=1<<VALUE_IS_VALUE; // eval to string rather than "1" or rest of phrase 
m["enum"]=1<<ENUM; // value must be in given list although two need notbe consistent  


} // BitVFlagMap
//enum NOUN_FLAGS
static void BitFlagMap(FlagMap & m)
{
m["ignore"]=1<<IGNORE; // mnemomnic or placeholder must have zero qty. 
m["syn"]=1<<SYNONYM; // accomodates a synonum from foreign ingredient list 
m["synonym"]=1<<SYNONYM;
m["dummy"]=1<<DUMMY; // machine generated to replace a missing one. 
m["units"]=1<<NEED_UNITS; // there is no concept of a serving must have units 
m["dmel"]=1<<DMEL; // something was wrong with this before earlier entries dmel 
m["corrected"]=1<<DMEL; // something was wrong with this before earlier entries dmel 
m["default"]=1<<WIN_TIE; // pick this silently in case of a tie 
m["alpha"]=1<<EXP_EVAL; // evaluation is exceptional ( not a number )  
m["text"]=1<<EXP_EVAL; // evaluation is exceptional ( not a number )  
m["shrimp"]=1<<SHRIMP_OK; // evaluation is exceptional ( not a number )  
// text here is part of flag string, other text is a description
// of the noun in human format


}
StrTy FlagBits(const IdxTy f) const 
{
Ss ss;
if (f!=0) ss<<" ";
FlagMap m;
BitFlagMap(m);
MM_LOOP(ii,m)
{
if ((f&((*ii).second))!=0)
{
// this picks up synonyms but should be obvious... 
// could reset bits as printed but that can masks non 2^n 
if (ss.str().length()!=0) ss<<"|";
 ss<<(*ii).first;
}
} // ii 
return ss.str();
}


#define NCASEEQ(x,y) case x : {++i; if (i<szl) y=l[i]; break; } 
#define NCASEEQI(x,y) case x : {++i; if (i<szl) y=myatoi(l[i]); break; } 
#define NCASEINC(x,y) case x : {++i; if (i<szl) y[l[i]]+=1; break; } 
#define NCASEDEC(x,y) case x : {++i; if (i<szl) y[l[i]]-=1; break; } 
#define NCASESPLIT(x,y) case x : {Loo::Split(i,y,l,szl);  break; } 
#define NCASEENUM(x,y) case x : {Loo::SplitEnum(i,y,l,szl);  break; } 
#define NCASESPLITPART(x,y,z) case x : {Loo::SplitPart(i,y,z,l,szl);  break; } 
#define NCASESPLITPARTXY(x,y) case x : {Loo::SplitPartFull(i,y,l,szl,canflags);  break; } 
#define NCASESPLITPARTXYLBL(x,y) case x : {Loo::SplitPartFullLbl(i,y,l,szl,m_name,canflags);  break; } 
#define NCASESPLITPARTRANKN(x,y) case x : {Loo::SplitPartRankName(i,y,l,szl,canflags);  break; } 
//#define NCASESPLITPARTLINE(x,y) case x : {Loo::SplitPartLine(m_line_parts,i,y,l,szl);  break; } 
#define NCASESPLITALTLINE(x,y) case x : {Loo::SplitFullPartLine(y,i,l,szl);  break; } 
//m_line_parts[w]=r.line(0);


/* the concern for hidden components motivates
machine readable comments and ingredient lists. 
Human readable ingredient lists for later grepping and maybe
parsing LATER.
Masses that reflect components and volumes that reflect whole.
Pehaps an "active" tag but consier metal salts where weight
are amount of metal or even a capsule content. Volume
needs to reflect entire formulation. 

*/

void Init(const Line & l) 
{
if (l.size()<1) return; 
m_name=l[0];

Init(l,m_name,1);
} // Init
void Init(const Line & _l,const StrTy & name, const IdxTy f) 

{
FlagMap fm,fmv;
BitFlagMap(fm);
BitVFlagMap(fmv);
m_flags=0;
m_vflags=0;
m_expand=0;
//m_default=1;
m_pcc=0;
// take out parents or pre-process here so params etc
// can have comma deferred or comment text
int skipping=0;
Line l;
 Ss sse;
const IdxTy sz_l=_l.size();
for(IdxTy i=0; i<sz_l; ++i)
{
const Word & w=_l[i];
sse<<w<<" ";
if (w.c_str()[0]=='#') { break; } 
if ( w=="(") { ++skipping; continue; } 
if ( w==")") { --skipping;
if (skipping<0) 
{
MM_ERR(" imbalanced parens in noun line "<<MMPR4(i,skipping,l.size(),sse.str()))
}
 continue; } 
if (skipping<=0) l.push_back(w);
}
if (skipping!=0) 
{
MM_ERR(" imbalanced parens in noun line "<<MMPR3(skipping,l.size(),sse.str()))
}

const IdxTy canflags=~0;
IdxTy szl=l.size();
//for(IdxTy i=0; i<szl; ++i)
for(IdxTy i=0; i<szl; ++i)
{
//MM_ERR(" noun input dump "<<MMPR2(i,l[i]))
// TODO this is dumb since the Ragged parse will expose this if quoted
if (l[i].c_str()[0]=='#') { szl=i; break; } 
}
// TODO  <=  f? 
if (szl<1) return; 

//m_name=l[0];
m_name=name;
CmdMap & m=Llut();
for(IdxTy i=f; i<szl; ++i)
{
// Ragged parse should group with quotes  
const Word & w=l[i];
// in args, the parens get soaked up so ok here 

const IdxTy code=m[w];
if (code==CANON)
{
if (m_alt_parts.size())
{
MM_ERR(" m_alt_parts using canon for name so define canon first ..."<<MMPR2(m_name,m_alt_parts.size())) 
} 
} // canon 
switch (code) { 
// parts need to work with adjectives... 
NCASESPLITPARTXYLBL(PART,m_alt_parts) 
NCASESPLITPARTRANKN(RPART,m_alt_parts) 
//NCASESPLITPARTLINE(LPART,m_line_parts) 
NCASESPLITALTLINE(LPART,m_alt_lines) 
NCASEEQ(CANON,m_canon) 
NCASEFLAG(SETFLAGS,m_flags,fm) 
NCASEFLAG(SETVFLAGS,m_vflags,fmv) 
// tokenized tree? too short still 
case IGNOREN: { m_flags|=(1<<IGNORE); break; } 
// reference value for contents and for quantity with no units. 
case SERVING:  
case DEFAULTQTY: { ++i; if (i<szl) {  m_default=Qdim(l[i]);
//MM_ERR(" setting serving size "<<MMPR3(l[i],m_default.qty(),m_default.sfx()))
 }
 break; } 
case PRIVATEUNIT: { ++i; 
if (i<szl)
{
Rb rb(l[i]);
rb.split_and_mark(',');
const IdxTy n=rb.string_count();
if (n&1) 
{
MM_ERR(" bad private unit "<<MMPR4(l[0],l[i-1],l[i],rb[n-1]))
//break;
}
for(IdxTy j=0; j<(n&(~1)); j+=2)
{
const StrTy sfx=rb[j];// =Qdim(l[i]); 
const StrTy cfac=StrTy(rb[j+1]);
//m_private_sfxs[sfx]=Qdim(rb[j+1]);
m_private_sfxs[sfx]=Qdim(cfac+"/"+sfx);
} // i 
} // i<sz

break; } 
NCASEINC(MODIFIER,m_mods);
NCASEDEC(NOT,m_nots);
NCASEINC(MFG,m_adjs);
NCASESPLIT(UNITS,m_adjs)
NCASESPLIT(ADJ,m_adjs)
NCASESPLIT(ADJECTIVE,m_adjs)
NCASEEQ(TEXT,m_text) 
NCASEEQI(EXPAND,m_expand) 
NCASESPLIT(GROUP,m_groups)
NCASEENUM(ENUMM,m_enum)
NCASEINC(CONFLICT,m_conflicts)
NCASEEQ(STARTDATE,m_startdate)
NCASEEQ(ENDDATE,m_enddate)
NCASEEQ(ACTIVE,m_active)
NCASEEQ(INGREDIENTS,m_ingredients)
//NCASEEQ(QTY,m_canon_qty) 
NCASEEQ(RUNITS,m_report_units) 
NCASEEQ(URL,m_url) 
NCASEEQ(DENSITY,m_density) 
NCASEEQ(BIBTEX,m_bibtex) 
NCASEEQ(BOMTEX,m_bomtex) 
case UPC : {Loo::SplitUPC(i,m_upc,l,szl);  break; } 
case MISC : {Loo::SplitUPC(i,m_misc,l,szl);  break; } 
default:
MM_ERR(" unknown noun mod "<<MMPR3(m_name,w,i))
} // switch 
} // i 
SetSelectors();
} // Init

#undef NCASEEQ
#undef NCASEEQI
#undef NCASEINC
#undef NCASEDEC
#undef NCASESPLIT
#undef NCASESPLITPART
#undef NCASESPLITPARTXY
#undef NCASESPLITPARTXYLBL
#undef NCASESPLITPARTRANKN
//#undef NCASESPLITPARTLINE
#undef NCASESPLITALTLINE





void SetSelectors()
{
SelectMap & m=m_selects;
m.clear();
IdxTy mask=SEL_GROUP;
MM_LOOP(ii,m_groups) { m[(*ii).first]|=mask; } 
mask=SEL_ADJ;
MM_LOOP(ii,m_adjs) { m[(*ii).first]|=mask; } 
mask=SEL_MISC;
MM_LOOP(ii,m_misc) { MM_LOOP(jj,(*ii).second) m[(*jj)]|=mask; } 
// add very high value for unique

 

} // SetSelectors


//void Di(
void Dump(Os & os, const IdxTy flags)  const
{ 
//Dump(os,flags); 
os<<MMZS(m_name); os<<MMZS(m_canon);
os<<MMN(m_expand);
os<<MMN(m_flags)<<FlagBits(m_flags);
os<<MMN(m_vflags); // <<FlagBits(m_flags);
os<<MMZS(m_startdate);
os<<MMZS(m_enddate);
os<<MMZS(m_density);
os<<MMZS(m_active);
os<<MMZS(m_report_units);
os<<MMZS(m_url);
os<<MMZS(m_bibtex);
os<<MMZS(m_bomtex); 
os<<MMZS(m_text);
os<<MMZS(m_expansion_name);
os<<MMZS(m_ingredients);

// this needs to be a string and will create a zero leading to printing
os<<MMZS(m_default); ///  os<<MMZS(m_canon);
//Loo::Dump(os,"m_parts",m_parts,0);
//Loo::Dump(os,"m_line_parts",m_line_parts,0);
Loo::Dump(os,"m_alt_parts",m_alt_parts,0);
Loo::Dump(os,"m_alt_lines",m_alt_lines,0);
Loo::Dump(os,"m_adjs",m_adjs,0);
Loo::Dump(os,"m_mods",m_mods,0);
Loo::Dump(os,"m_nots",m_nots,0);
Loo::Dump(os,"m_private_sfxs",m_private_sfxs,0);
Loo::Dump(os,"m_groups",m_groups,0);
Loo::Dump(os,"m_conflicts",m_conflicts,0);
Loo::Dump(os,"m_upc",m_upc,0);
Loo::Dump(os,"m_misc",m_misc,0);
Loo::Dump(os,"m_enum",m_enum,0);
}

// MEMBERS members 

StrTy m_name,m_canon; //,m_canon_qty;
//typedef std::map<StrTy, Line> LineParts;
//PieceMap m_parts;
//LineParts m_line_parts;
// these store qdim so you need naes available/m_alt_parts
PieceMaps m_alt_parts;
LineMaps m_alt_lines;
IdxTy m_expand;
AdjMap m_adjs;
ModMap m_mods;
ModMap m_nots;
IdxTy m_flags,m_vflags;
//D m_default;
Qdim  m_default;
// these are actually conversion factors
SfxMap m_private_sfxs;
CharClass * m_pcc;
MiscMap m_misc;
EnumMap m_enum;
// these are user define.  
GroupMap m_groups; // TODO FIXME need a member map 
SelectMap m_selects; // TODO FIXME need a member map 
ConflictMap m_conflicts;
IndustryCodes m_upc;
// density is a special case of conversion factors or 
// actually similar to content - say energy kcal/cup same as vitamin/kg
// want to tokenize too /to
StrTy m_startdate,m_enddate,m_density,m_report_units,m_url,m_bibtex,m_bomtex;
StrTy m_text,m_expansion_name;
StrTy m_active,m_ingredients;
//friend  mjm_compiled_noun_collection<Tr>;
}; // mjm_compiled_noun

template <class Tr> 
class mjm_compiled_noun_collection 
{
typedef mjm_compiled_noun_collection Myt;
typedef mjm_compiled_noun<Tr> Cnoun;
//typedef  mjm_ragged_forms<Tr>  Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::OsTy Os;
 typedef typename Tr::Ofs Ofs;

typedef std::vector<Cnoun> Cnouns; 
typedef std::map<StrTy, Cnouns> CnounMap;
typedef typename CnounMap::iterator Itor;
typedef typename CnounMap::const_iterator CItor;

typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line;
typedef StrTy Word;
typedef std::vector<const Cnoun*> PtrVec;
typedef mjm_var_parse<Tr> CharClass;
typedef typename CharClass::read_buffer Rb;

public:

//typedef std::map<StrTy, PtrVec> MissingMap;
typedef std::map<StrTy, int> MissingMap;
const IdxTy size() const { return m_map.size(); }
Cnouns & operator[](const StrTy & k) { return m_map[k];}
const Cnouns & operator[](const StrTy & k)const  { return m_map[k];}
const Cnouns * operator()(const StrTy & k)const  
{auto ii=m_map.find(k); if (ii==m_map.end()) return 0;  return &  (*ii).second;}
Cnouns * operator()(const StrTy & k)  
{auto ii=m_map.find(k); if (ii=m_map.end()) return 0;  return (*ii);}
MissingMap missing() const { return Missing(); } 
Itor find(const StrTy & k) { return m_map.find(k); } 
CItor find(const StrTy & k)const { return m_map.find(k); } 

Itor begin() { return m_map.begin(); }
CItor begin()const  { return m_map.begin(); }

Itor end()  { return m_map.end(); }
CItor end()const  { return m_map.end(); }
void clear() { m_map.clear(); } 

void load(const Ragged & r) { Load( r); }
template <class Tm> void index_groups(Tm & m ) { IndexGroups(m); } 
template <class Tm> void count_groups(Tm & m ) { CountGroups(m); } 
void dump(Os & os, const IdxTy flags) const { Dump(os,flags); } 
private:
MissingMap Missing() const 
{ 
MissingMap m;
MM_LOOP(ii,(*this))
{
const Cnouns & c= (*ii).second;
const StrTy & noun=(*ii).first;
MM_LOOP(qq,c)
{
const typename Cnoun::parts_type & p=(*qq).pieces();
MM_LOOP(jj,p)
{
const StrTy & basis=(*jj).first;
MM_LOOP(kk,(*jj).second)
{
const StrTy & s=(*kk).first;
const Cnouns * cpp=(*this)(s);
if (cpp==0) ++m[s];

}

} // jj 
} // qq 
} /// ii 

return m; 
} // Missing
 
template <class Tm> void IndexGroups(Tm & m )
{
MM_LOOP(ii,m_map)
{
//const StrTy  & nm=(*ii).first;
const Cnouns & vec=(*ii).second;
MM_LOOP(jj,vec)
{
const Cnoun & cn=(*jj);
const Cnoun * p=&cn;
MM_LOOP(kk,cn.groups())
{
//m[nm].push_back(p);
m[(*kk).first].push_back(p);
} // kk 
} // jj 
} // ii 

} // IndexGroups
 
template <class Tm> void CountGroups(Tm & m )
{
MM_LOOP(ii,m_map)
{
//const StrTy  & nm=(*ii).first;
const Cnouns & vec=(*ii).second;
MM_LOOP(jj,vec)
{
const Cnoun & cn=(*jj);
const Cnoun * p=&cn;
MM_LOOP(kk,cn.groups())
{
//m[nm].push_back(p);
++m[(*kk).first];
} // kk 
} // jj 
} // ii 

} // CountGroups



void Load(const Ragged & r)
{
MM_LOOP(ii,r)
{
//IdxTy i=(*ii).second;
const Line & l =(*ii); // m_nouns[i]; 
if (l.size()==0) continue;
// comment needs to begin a word, ZZ
if ( l[0].c_str()[0]=='#') continue; 
//m_cnouns[i]=Cnoun(l);
// TODO FIXME allow l[0] to be a comma sep list of equal names
// should work for units too. 
//m_map[l[0]].push_back(Cnoun(l));
Rb rb(l[0]);
rb.split_and_mark(',');
const IdxTy n=rb.string_count();
for(IdxTy j=0; j<n; ++j)
{
const StrTy name=rb[j];// =Qdim(l[i])
m_map[name].push_back(Cnoun(l,name,1));
} // j 

} // ii 

} // Load

void Dump(Os & os, const IdxTy flags) const 
{ 
MM_LOOP(ii,m_map) { 
const StrTy nm=(*ii).first;
const auto & v =(*ii).second;
MM_LOOP(jj,v)
{
os<<nm<<" ";
(*jj).dump(os,flags);
os<<CRLF; 
} // jj 
} // ii
//Dump(os,flags); 

} // Dump  


CnounMap m_map;

}; // mjm_compiled_noun_collection


#endif

