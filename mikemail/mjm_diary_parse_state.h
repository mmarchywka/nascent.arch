#ifndef MJM__DIARY_PARSE_STATE_H__
#define MJM__DIARY_PARSE_STATE_H__

#include "mjm_globals.h"
#include "mjm_compiled_noun.h"
#include "mjm_thread_util.h"
#include "mjm_ragged_forms.h"
#include "mjm_snacks_accumulator.h"
#include "mjm_constrained_file.h"
#include "mjm_unit_crap.h"
// faster than calling date in bash doh 
#include "mjm_calendar.h"

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



template <class Tr>

class mjm_diary_parse_state
{
//class mjm_diet_diary_form  : public mjm_ragged_forms<Tr> 
//{
// typedef mjm_diet_diary_form Myt;
//typedef  mjm_ragged_forms<Tr>  Super;
typedef mjm_diary_parse_state Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// why the fck wont' these inherit fck 
//typedef typename Super::Ragged Ragged;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line;
typedef StrTy Word;
typedef std::vector<Word> Words;
typedef std::map<StrTy, IdxTy> DictTy;
//typedef typename Super::CharClass CharClass; 
//typedef mjm_constrained_file<Tr> BackingFile;
//ZZZZBackingFile m_file;
//typedef typename CharClass::read_buffer Rb;

typedef mjm_unit_crap<Tr> Units;
typedef typename Units::qty_dim Qdim;
typedef std::map<StrTy, int > ModMap;

//typedef compiled_noun Cnoun;
typedef mjm_compiled_noun<Tr> Cnoun;
typedef mjm_diary_parse_error<Tr> parse_error;
//typedef std::vector<Cnoun> CnounVec;
typedef std::vector<Cnoun> Cnouns; 
typedef std::map<StrTy, Cnouns> CnounMap;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
//public:
//#define MACRONOUNS QUOTE(MACRONOUN)
// need an "eaten vs served percentage and prepared ahead of time?
//enum { IGNORE,LPAREN,RPAREN,RESERVED,NOUN,TIME,DATE,MACRO,MACROMOD,NAME,AMPERSAND,BSLASH, BAD=Super::BAD};
//enum { MACRONOUN=0};
//MM_ERR("class "<<MMPR3(w,c,class_enumval(c)))
//enum { MACROMASK=CharClass::UC|CharClass::DIGIT, DIGIT=CharClass::DIGIT};
//enum { LC=CharClass::LC};
//enum { UC=CharClass::UC};


typedef ::parse_settings<Tr> parse_settings;
typedef parse_settings PO;

typedef mjm_muqed_phrase<Tr> Phrase;

typedef typename Phrase::file_finder_type Finder;





typedef parse_error Error;
typedef std::vector<Error> Errors;
enum {DATED=1, CONSUMER=(1<<1),RESERVED=(1<<2),
NOUNED=(1<<3),LATEX=(1<<4), ADJ=(1<<5), MACROEVAL=(1<<6) };
// MACRO def must be at DATED&!CONSUMER&!RESERVED&!NOUNED 
//Test(DATED)&&!Test(CONSUMER|RESERVED|NOUNED)
bool Test( const IdxTy mask) const { return (state&mask)!=0; } 
bool Test(const IdxTy state, const IdxTy mask) const { return (state&mask)!=0; } 
bool TestAll( const IdxTy mask) const { return (state&mask)==mask; } 
bool TestAll(const IdxTy state, const IdxTy mask) const 
{ return (state&mask)==mask; } 

//typedef typename compiled_noun::evaluation E; 
typedef typename Cnoun::evaluation E; 
// date,consumer,noun,time, qty (  , units )
typedef std::map<D,E> Emap;
typedef std::map<StrTy,Emap> NounMap;
typedef std::map<StrTy, NounMap> ConsumerMap;
typedef std::map<StrTy,ConsumerMap> DateMap;
typedef std::vector<E> Evec;
//ne=dm[cdates][consumer][canon][tmft()]; 
public:
typedef mjm_snacks_accumulator<Tr> DateMap2;
mjm_diary_parse_state():form(0),line(0),nparen(0), d(0),j(0),j0(0),state(0)
,m_date(0), m_old_date(0),time(0),depth(0),
markup(0),blank(0),canonical(0),mom(0),root(this),
per_noun_os(0),debug_os(0),
dmlocal(0),dmglob(0),exp(0),m_top(true),m_pu(0),m_pfinder(0)  {} 
// macro def needs to get the line numer 
mjm_diary_parse_state(Myt & p ):form(p.form),line(p.line),nparen(0), d(0),j(0),j0(0),state(0)
,m_date(p.m_date), m_old_date(p.m_old_date),time(p.time),depth(p.depth+1),
markup(p.markup),blank(p.blank),canonical(p.canonical)
,mom(&p),root((p.root==&p)?mom:p.root) ,
per_noun_os(0),debug_os(0),dmlocal(0),
dmglob(0),exp(1),m_top(false),m_pu(p.m_pu),m_pfinder(p.m_pfinder)   {
//MM_ERR(" copy ctor called ")
m_phrase.finder(m_pfinder); 
} 
~mjm_diary_parse_state() { delete dmlocal; if (m_top) free();  }
void free()
{
delete markup;
delete blank;
delete canonical;
delete dmglob;
delete dmlocal; dmlocal=0;
}
bool warn_depth() const { return (depth>30); }
bool stop_depth() const { return (depth>40); }
void keep_blank() { blank=0; } 
void finder(Finder * p ) { m_pfinder=p; m_phrase.finder(p);  }
Finder *  finder() { return m_pfinder; }
const Finder *  finder() const  { return m_pfinder; }
void adj(const Word & w ) 
{ state|=ADJ;m_phrase.add_word(w,j);  
//admap[w]=j;  adjectives.push_back(w); 
}
enum { IGNORE=1 };
bool skipping() const { return nparen!=0; } 
const StrTy & consumer_name() const { return consumer; } 
// need to verify this is in a name location and NOT a adjective
void have_name(const StrTy & w) {
if (!have_date()) add_error(w,j," need date before name "); 
//MM_ERR( " should ba name "<<MMPR2(w,Test(state,CONSUMER)))
if ( Test(state,CONSUMER))  adj(w); 
else { //if (!expanding()) 
consumer=w; 
state|=CONSUMER; } 
}
// comments, note , outcome, DMEL, AKA obsolete etc
// aggregate and do something 
void have_reserved(const StrTy & w,const Line & l ) {
if (!have_date()) add_error(w,j," need date before reserved word "); 
if (!Test(state,CONSUMER)){ consumer=w; state|=RESERVED; 
// TODO FIXME should not do parsing or actions here but easy for now. jjjj
//MiscTypes m_dmel, m_comments; MiscConTypes m_notes;
if (w=="NOTE") if (l.size()>(j+1) ) m_notes[l[j+1]].push_back(l);
if (w=="DMEL") m_dmel.push_back(l);
if (w=="COMMENT") m_dmel.push_back(l);
/////

}  // have_reserved
else // if we already have a consumer, this is not a reserved
// word but rather bad choice of macro or adjective
{
add_error(w,j,"can not use reserved word here   ");
}

}
void have_latex(const StrTy & w) {
state|=LATEX;

}
enum {BAD=~0 } ;

//void organize_adjectives( ) { organize_adjectives(m_phrase); } 
void organize_adjectives( ) { m_phrase.organize_adjectives(*this); } 
void organize_adjectives( Phrase &phrase)
{
//organize_adjectives(phrase.prop(),phrase.a1(),phrase.a2()); 

} // organize_adjectives phrase
#if 0 
#endif


D tfmt() const { return int(time*60)/60; } 
// a new data point for date, time, consumer, noun, qty 
void have_noun(E & ne) {
const bool ignore=ne.noun().ignore()||(!Test(state,CONSUMER));
// could just save ne's in big collection 
// date,consumer,noun,time, qty (  , units )
// really want a tokenized ragged with float columns .... 
state|=NOUNED;
if (!ignore) { 
if (depth==0) // top level so evaluate 
{
//if (dm==0) dm= new DateMap2();
if (dmglob) dmglob->add(*this,ne,(ne.noun()));
if (false) if (debug_os!=0)  
	(*debug_os)<<" DEBUG NOUN "<<MMPR4(cdates,consumer,ne.canon(),tfmt())
		<<MMPR2(ne.qty(),ne.sfx())<<CRLF; 
if (tfmt()==0)
	add_error(ne.canon(),0," time of zero likely an error ");
if (m_date==0)
	add_error(ne.canon(),0," date is still zero at this point  ");


// MM_ERR(" FICL "<<MMPR((debug_os!=0)))
} 
else { evec.push_back(ne); } // in a recurser
} // ignore
m_phrase.clear();
//adjectives.clear();
//admap.clear();
} // have_noun 1 param 

//  ne is based on the adjectives for the macro which was expanded into ps2
// this needs to allow for variables etc. 
void have_noun(E & ne, Myt & ps2) 
{
state|=NOUNED;
if (ps2.time>time) time=ps2.time;
MM_LOOP(ii,ps2.evec)
{
// offered and ate should be in ne not in ps2 entries, 
// just multiply ne terms by each ps2 term
E neq=(*ii);
//neq=neq*ne.offered();
neq.scale(ne.has_of()?ne.of():0,ne.qty());
//operator
have_noun(neq);

} // ii 
m_phrase.clear();
//adjectives.clear();
//admap.clear();
// ps2 not reusable will go out of scope soon 
} // have_noun 2 param 
template <class Tv> 
void have_noun4(E & ne, Myt & ps2, const Tv & minus, const Tv & plus) 
{
MM_ONCE(" need to implement macro modifiers",) 
state|=NOUNED;
if (ps2.time>time) time=ps2.time;
MM_LOOP(ii,ps2.evec)
{
// offered and ate should be in ne not in ps2 entries, 
// just multiply ne terms by each ps2 term
E neq=(*ii);
//neq=neq*ne.offered();
neq.scale(ne.has_of()?ne.of():0,ne.qty());
//operator
have_noun(neq);

} // ii 
m_phrase.clear();
//adjectives.clear();
//admap.clear();
// ps2 not reusable will go out of scope soon 
} // have_noun 4 param 


void have_time( const StrTy & w, const D & val)
{
if (val>time) time=val;  
if (!have_date()) add_error(w,j," no date prior to time ");

}
const StrTy date() const { return cdates; } 
bool have_date() const { return Test(DATED); } 
void reset_old_date(const IdxTy old=0 ) { m_old_date=old; } 
void have_date(const IdxTy ival,const StrTy & w) 
{ 
// if we already have a date, need to punt to adjectives
if (Test(DATED)) add_error(w,j," alreadhy have a date");
//dates=w; cdates=w;  m_date=ival; state|=DATED;
//if (m_date<m_old_date)
if (ival<m_old_date)
{
Ss ss;
ss<<MM_MARK<<" new date less than old date "<<MMPR4(m_date,m_old_date,cdates,w)
<<MMPR2(have_date(),(state&DATED));
//add_error(w,j," new date less than old date");
add_error(w,j,ss.str().c_str());
}
dates=w; cdates=w;  m_date=ival; state|=DATED;
} // have_date

bool have_macro_def(const StrTy & w, DictTy & dmacros )
{ // true of def, false if requires expansion 
tname=w; // stupid kluge should be saved ... 
// in this state, it should be a definition not a usage 
const bool want_def=Test(DATED)&&!Test(CONSUMER|RESERVED|NOUNED);
//auto ii=dmacros.find(w);
// these need to be defined everyday for a macro recipe. 
//bool found=(ii!=dmacros.end());
// with prefix adjectives, the expansion already has the prefactors
// for evaluation .. 
return want_def;

} // have_macro_def
void rparen() {unparen(); }
void paren() {++nparen; }
void unparen() {
// verify not zero ... 
// TODO this may not work if the left is before the consumer name... 
if (nparen==0) 
	add_error("",j," too many right parens ");
//else 
// 2020-05-06 let it go negative anyway... 
	--nparen; 

}
//StrTy dump() const { Ss ss; ss<<MMPR4(d,j,state,date)<<MMPR2(time,adjectives.size());  return ss.str(); }
StrTy dump() const { Ss ss; ss<<MMPR4(d,j,state,m_date)<<MMPR2(time,m_phrase.words());  return ss.str(); }
// at end of processing a line
void check_eol() { 
// do someething with dangling adjectives 
//if (adjectives.size()!=0) add_error("",j," unused adjectives at end of line  ");
if (m_phrase.words()!=0) add_error("",j," unused adjectives at end of line  ");
// check and reset paren count 
//if (nparen!=0) add_error("",j," dangling  parens ");
if (nparen>0) add_error("",j," dangling  left  parens ");
if (nparen<0) add_error("",j," dangling  right parens ");
// 2021-06-29 line starting with bad date failed to make error message
if (!have_date()) add_error("",j," no date found for line ");

} 
void new_line() {
check_eol();
nparen=0;
old_err[line]=errors;
errors.clear();
m_old_date=m_date;
d=0; j=0; state=0; m_date=0; time=0;
m_phrase.clear();
//adjectives.clear();
//admap.clear();
m_notes.clear();
m_dmel.clear();
m_comments.clear();
dates=""; cdates=""; entry_type="";
consumer="";
++line; }
IdxTy errors_size() const { return errors.size(); } 
void add_error(const StrTy & w, const IdxTy pos, Ss & ss)
{add_error(w,pos,ss.str().c_str()); } 
void add_error(const StrTy & w, const IdxTy pos, const char * msg)
{ 
MM_ERR(MMPR3(line,dates,consumer)<<MMPR3(w,pos,msg))
errors.push_back(Error(w,pos,msg)); }
void add_errors(Errors & e) { MM_LOOP(ii,e) errors.push_back(*ii); e.clear(); }
void add_errors(Myt & that) { Errors & e=that.errors;  
Ss ss;
ss<<tname<<"("<<j<<")";
const StrTy pfx=ss.str(); // (that.evec.size()>0)?that.evec[0].canon():"";
MM_LOOP(ii,e)
{
//auto ee=(*ii);
 errors.push_back((*ii).scope(pfx));
}
 e.clear(); }
bool expanding() const { return exp!=0; } 


Errors errors; 
typedef std::map<IdxTy,Errors> OldErr;
OldErr old_err;
void markup_form(const IdxTy flags=0 )
{
Ragged & r=*form;
Ragged & d=*markup;

MM_SZ_LOOP(i,r,szr)
{
const Line & l=r[i];
d.add(l);
auto ii=old_err.find(i);
if (ii==old_err.end()) continue;
MM_LOOP(jj,(*ii).second)
{
Line x;
x.push_back("##");
x.push_back((*jj).dump());
d.add(x);
} // jj 
} // i 

} // markup




Ragged * form;

IdxTy line,nparen;
IdxTy d,j,j0;
IdxTy state;
IdxTy m_date,m_old_date;
// date string as found and canonical version 
StrTy dates,cdates,consumer;
StrTy entry_type;
D time; // wtf
IdxTy depth; // recursion for recipe expansion 
Phrase m_phrase;
//Words adjectives;
//ModMap admap;
Ragged * markup;
Ragged * blank;
Ragged * canonical;
Myt * mom,* root;

OsTy * per_noun_os;
OsTy * debug_os;
Evec evec;
DateMap*  dmlocal; // don't really want this here for recursion.. 
DateMap2*  dmglob; // don't really want this here for recursion.. 
// TODO FIXME NOTE, DMEL, COMMNENTS...
typedef std::vector<  Line > MiscTypes;
typedef std::map<StrTy , std::vector< Line > >  MiscConTypes;
MiscTypes m_dmel, m_comments;
MiscConTypes m_notes;
/////
IdxTy exp;
StrTy tname;
bool m_top;
Units * m_pu;
Finder * m_pfinder;
}; // parse_state 
 




#endif // MJM_..._H__ 
