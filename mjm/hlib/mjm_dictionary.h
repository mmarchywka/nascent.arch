#ifndef MJM_DICTIONARY_H__
#define MJM_DICTIONARY_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>

/*
*/

/*
g++ -DTEST_DICTIONARY__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_dictionary.h
*/

////////////////////////////////////////////////////////////////

class dict_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
dict_params( const StrTy & nm) : Super(nm) {}
dict_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  

//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy output_label() const { return m_map.get_string("output_label","fick"); }
//bool always_dump_to_file() const { return m_map.get_bool("always_dump_to_file",!true); }

bool skip_old() const { return m_map.get_bool("skip_old",true); }
bool print_dog_days() const { return m_map.get_bool("print_dog_days",!true); }
bool accumulate_dog_days() const { return m_map.get_bool("accumulate_dog_days",true); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }

StrTy snack_log() const { return m_map.get_string("snack_log","snacks_log"); }
		//if (skip_old) if (date==StrTy("2017-04-22")) skipping=false;
StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }

// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

//ss<<"time_step="<<time_step()<<sep;
ss<<"skip_old="<<skip_old()<<sep;
ss<<"print_dog_days="<<print_dog_days()<<sep;
ss<<"accumulate_dog_days="<<accumulate_dog_days()<<sep;
ss<<"log_commands"<<log_commands()<<sep;
ss<<"snack_log="<<snack_log()<<sep;
ss<<"start_date="<<start_date()<<sep;
ss<<"end_date="<<end_date()<<sep;
return ss.str();
}


}; // snack_params

#if 0 
// record exceptions and warnings for data or code problems diagnosis 
class data_model_error_log 
{
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef double D;
typedef unsigned int IdxTy;
//typedef std::map<StrTy, D > DoseMap;
typedef std::vector<StrTy > Words;
class entry
{
public:
StrTy dump() const
{
Ss ss;
ss<<MMPR3(word,date,problem);
for (IdxTy i=0; i<adjectives.size(); ++i) ss<<" "<<adjectives[i];
return ss.str();
}
StrTy word;
Words adjectives;
StrTy date;
StrTy problem;
 
}; // entry 
typedef entry Event;
typedef std::vector<Event > EventList;
typedef std::map<StrTy, IdxTy> ProblemCount;
public:
void event(const StrTy & p , const StrTy & word)
{ Event x; x.problem=p; x.word=word; event(x); }
void event(const StrTy & p , const StrTy & word, const StrTy & d)
{ Event x; x.problem=p; x.word=word; x.date=d; event(x); }
void event(const StrTy & p , const StrTy & word, const StrTy & d,const Words & words)
{ Event x; x.problem=p; x.word=word; x.date=d; x.adjectives=words; event(x); }

void event(const Event & x)
{
	m_events.push_back(x);
	++m_problems[x.problem];
}

StrTy string() const
{
Ss ss;
const IdxTy sz=m_events.size();
for (IdxTy i=0; i<sz; ++i)
{
ss<<i<<" "<<m_events[i].dump()<<CRLF;
}

return ss.str();

}

EventList m_events;
ProblemCount m_problems;
}; // data_model_error_log
#endif
/*
namespace  nested_map_iterator
{

typedef std::string StrTy;
typedef std::stringstream Ss;
typedef double D;
typedef unsigned int IdxTy;
typedef std::map<StrTy, D > TerminalMap;
typedef TerminalMap::iterator TerminalMapItor;
typedef StrTy Word;
typedef std::vector<Word> Words;
//Words m_keys;
//public:

template <typename  Tx, typename Ty> Tx next(Words & w, Ty & ii) 
	//{ m_keys.push_back((*ii).first); return (*ii).second; }
	{ w.push_back((*ii).first); return next((*ii).second); }
template < > D next(Words & w, TerminalMapItor & ii) 
	{ w.push_back((*ii).first); return (*ii).second; }


}; //mested_map_iterator
*/

/*
class dog_day
{
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef double D;
typedef unsigned int IdxTy;
typedef std::map<StrTy, D > DoseMap;
public:
dog_day()  {} // only for stl crap 
dog_day(const StrTy & n, const StrTy & d) : m_date(d),m_name(n) {}
// canonical name and amount 
void inc(const StrTy & k, const D & v) { m_doses[k]+=v; }  
void set_dose(const StrTy & n,const D & v ) {m_doses[n]=v;}
void add(const dog_day & that, const D & fac)
{
for (auto ii=that.m_doses.begin(); ii!=that.m_doses.end(); ++ii)
	{
auto key=(*ii).first;
m_doses[key]+=fac*(*ii).second;
//MM_MSG(MMPR4(m_name,key,fac,(*ii).second)<<MMPR2(that.m_name,m_doses[key]))
}
}
template <class Ty> void add(const dog_day & that, const D & fac, Ty & handlers)
{
for (auto ii=that.m_doses.begin(); ii!=that.m_doses.end(); ++ii)
	{
auto key=(*ii).first;
auto hh=handlers.find(key);
if (hh==handlers.end()) { m_doses[key]+=fac*(*ii).second; continue; }
//MM_MSG(MMPR4(m_name,key,fac,(*ii).second)<<MMPR2(that.m_name,m_doses[key]))
const auto * ha=(*hh).second;
//void share_doses(D & d, const D & f, const D & amt) { d+=f*amt; }
ha->share_doses(m_doses[key],fac,(*ii).second);
}


}

 StrTy string( ) const
{
const StrTy sep=" ";
Ss ss;
ss<<m_date<<sep;
ss<<m_name;
ss<<dose_string();
if (true ) return ss.str(); 
for(auto ii=m_doses.begin(); ii!=m_doses.end(); ++ii)
{
ss<<sep<<(*ii).first<<"="<<(*ii).second;
}

return ss.str();

}

 StrTy dose_string( ) const
{
const StrTy sep=" ";
Ss ss;
for(auto ii=m_doses.begin(); ii!=m_doses.end(); ++ii)
{
ss<<sep<<(*ii).first<<"="<<(*ii).second;
}

return ss.str();

}



template <class Ty>  StrTy string(const Ty & v ) const
{
if (v.size()==0) return string();
const StrTy sep=" ";
Ss ss;
ss<<m_date<<sep;
ss<<m_name;
ss<<dose_string(v);
return ss.str();
}
template <class Ty>  StrTy dose_string(const Ty & v ) const
{//for(auto ii=m_doses.begin(); ii!=m_doses.end(); ++ii)
if (v.size()==0) return string();
const StrTy sep=" ";
Ss ss;
for(IdxTy i=0; i<v.size(); ++i)
{
	auto ii=m_doses.find(v[i]);
	if (ii!=m_doses.end())	
		{ ss<<sep<<(*ii).first<<"="<<(*ii).second; } 
	else { ss<<sep<<v[i]<<"=0"; }
}

return ss.str();

}



StrTy m_date;
StrTy m_name;
DoseMap m_doses;

}; // dog_day
*/

// reserved words or specific things with meanings different
// from canonical names .

class literal_map
{
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
public:
// TIME etc are not included here although "noon" and "midnight" etc
// may be useful 
enum LITERAL_TYPE { ADJECTIVE=1, QUANTIFIER=2, NOUN=3, LMACRO=5,IGNORE=4, NONE=0 };
typedef std::map<StrTy,LITERAL_TYPE > Literals;
typedef std::map<StrTy,IdxTy > FreqInfo;
typedef std::map<int,FreqInfo > FreqClassInfo;
template <class Ty> void add_adjective(const  Ty & nm) { m_lit[StrTy(nm)]=ADJECTIVE; }
template <class Ty> void add_noun(const  Ty & nm) { m_lit[StrTy(nm)]=NOUN; }
template <class Ty> void add_ignore(const  Ty & nm) { m_lit[StrTy(nm)]=IGNORE; }
template <class Ty> void add_adjectives(const  Ty & v) 
	{ for (auto ii=v.begi(); ii!=v.end(); ++ii) {m_lit[StrTy(*ii)]=ADJECTIVE; }}

bool is_adjective(const StrTy & a )  const 
{ 
	auto ii=map().find(a);
	if (ii==map().end() ) return false; 
	return  ((*ii).second== ADJECTIVE ); 
}  

bool is_noun(const StrTy & a ) 
{ 
	if (map().find(a)==map().end() ) return false; 
	return  (map()[a]== NOUN ); 
}  

Literals & m() { return m_lit;}
Literals & map() { return m_lit;}
const Literals & map() const { return m_lit;}
		void inc_word(const IdxTy j,const IdxTy & name_class,const StrTy & f
	, const StrTy & cn)
//		++m_words[j][name_class][f];
{		++m_words[j][name_class][f]; m_class_names[name_class]=cn;}

void stream_counts(OsTy & os)
{
	const IdxTy sz=5;
	for (IdxTy pos=0; pos<sz; ++pos)
	{
		const FreqClassInfo & fcmap=m_words[pos];
		for( auto kk=fcmap.begin(); kk!=fcmap.end(); ++kk)
		{
			const int name_class=(*kk).first;
			const FreqInfo & fmap=(*kk).second; // m_words[i];
			IdxTy j=9; 
			for (auto ii=fmap.begin(); ii!=fmap.end(); ++ii)
			{
				const StrTy word=(*ii).first;
				const  IdxTy count=(*ii).second;	
				os<<MMPR(name_class)<<" "<<name_class_name(name_class)
					<<" "<<MMPR4(pos,j,word,count)<<CRLF;
				++j;
			}
		}
	}
}

// fix this crap lol 
StrTy name_class_name(const IdxTy i ) { return m_class_names[i]; }



Literals m_lit;
FreqClassInfo m_words[5];
std::map<IdxTy, StrTy> m_class_names;

};

class mjm_dictionary_entry
{
typedef mjm_dictionary_entry Myt;
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef StrTy Word;
typedef StrTy GrammarType;

// parts of speech this can be 
typedef std::map<GrammarType,D> UsageProb;
typedef std::map<Word,D> Synonyms;
public:

private:
Word m_name;
UsageProb m_types;
Synonyms m_synonyms;

}; // mjm_dictionary_entry


/*
Currently a handler goes with a given canonical name which also
is the aggrgation unit for the parser. Ideally many food and drug
items fall into specific types such as salmon or vitamin supplements.
However, the multi B for example has ingreidents of several types
to aggregate individually. 
 So making a handler name different from aggregation name may be useful.

n-nouns handled by on hanndler ( cname ).
each non cound contribute to several aggregation names ( aname). 
The handler needs to maintain composition info for each noun
it handles  and be able to update all anames. 

This requires composition info on the trivial names / nouns either separetely
or within the handler.
noun/trivial -> handler cname -> various aggregation anames  

*/

/*
Handlers for a group of ingredients with same canonical name. Handler
sorts out trivial names based on value and date and other context to
create a quantity for tha cname- usually mg but could be anything.
Eventually each trivial name will map to several canonical names
as ingredients or components map out variously. 

Also needs to some error bars or quality algebra. 

*/
/*

class drugs
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
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

typedef  data_model_error_log Dmel; 
typedef std::map<StrTy, D> BrandQty;
typedef std::map<StrTy, D> IgnoreMap;
typedef literal_map LiteralClass;
typedef LiteralClass::LITERAL_TYPE LIT_TYPE;
enum  LITERAL_TYPE { ADJECTIVE=LIT_TYPE::ADJECTIVE, QUANTIFIER=LIT_TYPE::QUANTIFIER
		, NOUN=LIT_TYPE::NOUN, LMACRO=LIT_TYPE::LMACRO
		,IGNORE=LIT_TYPE::IGNORE, NONE=LIT_TYPE::NONE };

typedef std::map<StrTy, D> SfxMap;
typedef std::vector<IdxTy> Used;

public:
typedef std::vector<StrTy> Words;
drugs(): m_name("none"),m_qtu(1),m_ignore_share(false) {}
drugs(const StrTy & n ): m_name(n),m_qtu(1)
		,m_ignore_share(false) {}
drugs(const D & qtu): m_name("none"),m_qtu(qtu)
		,m_ignore_share(false) {}
drugs(const StrTy & n, const D & qtu): m_name(n)
		,m_qtu(qtu),m_ignore_share(false) {}
virtual ~drugs() {}
// if special ajective analysis is needed for this group 
virtual bool custom_qty() const { return false; }
// this should return quality indication for last call to get quantity 
enum QUALITY_EST {DEFAULTS=0};
virtual QUALITY_EST quality() const 
{
return DEFAULTS;
}
virtual void share_doses(D & d, const D & f, const D & amt) const  {
if (m_ignore_share) { d+=amt; } else { d+=f*amt;} }
virtual void ignore_share(const bool tf) {m_ignore_share=tf; }
// this is always needed as it relates COUNT to quantity such as mg. 
//virtual D qty_to_units() const { return m_qtu; }
virtual D qty_to_units(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{ 
auto ii=m_qtu_map.find(w);
if (ii!=m_qtu_map.end()) return (*ii).second; 
return m_qtu; 
}
// add the mg/count value for a given trivial name to this canonical handler 
// this is only for simple cases as trivial names could change over time
// or with adjectives ( NOW vs GNC foo)
// for variant name w potency or mg/count is this s
virtual void add_potency(const StrTy & w, const D & s) { m_qtu_map[w]=s; } 
// () const? 
// this can be ad hoc logic to increase  val based on all context
// although cname should not be updated, it was used to find this handler  
//virtual bool safe_to_ignore(const StrTy & a) { return false;}
virtual bool safe_to_ignore(const StrTy & a) { return (m_ignores.find(a)!=m_ignores.end());}


virtual D qty_to_units2(Used & adj_used, D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
	// see if the trvial name w implies a count-to-qty number 
	auto ii=m_qtu_map.find(w);
	if (ii!=m_qtu_map.end()) return (*ii).second; 
	// othewise sort through the things  although should be 
	// zero if map above is used  let adj_used sort it out 
	const IdxTy asz=adjectives.size();
	for(IdxTy i=0; i<asz; ++i)
	{
		const StrTy & a= adjectives[i];	
		auto ii=m_qtu_map.find(a);
		if (ii!=m_qtu_map.end()) { ++adj_used[i]; return (*ii).second;  } 
		// a suffix means a size 
		const StrTy sfx=suffix(a);
		const D conv=sfx_map()[sfx];
		if (sfx!="x") if (conv<=0) continue; // if this is NOT a suffix then it is a count 
		if ((sfx==a)&&(i>0)) // theere is white space between number and sfx 
		{
			++adj_used[i];
			++adj_used[i-1];
			// need to check this is a NUMBER 
			D vall=atof(adjectives[i-1].c_str());
			return vall*conv;	
		}
		D vall=q_x_sz(a,a.length()-sfx.length()); // this calculates the a x v syntax result 
		if (vall>=0)
		{
			++adj_used[i];
			return vall*conv;
		}
	}
	return m_qtu; 
}


virtual D count(Used & adj_used, D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
	D explicit_count=-1;
	const IdxTy asz=adjectives.size();
//	bool db= (cname=="Mg");
//	bool db2=false;
	for(IdxTy i=0; i<asz; ++i)
	{
		const StrTy & a= adjectives[i];
//		if (a==".2") db2=true;
//		if (db&&db2) MM_ERR(" debugg "<<MMPR3(a,i,asz))	
		const char * s=a.c_str();
		// pure number with decimal or fraction, "x" handled in other thing 
		char clast=firstnondigit(s);
		//if (containsnondigit(s)) continue;
		if ((clast!=0)&&(clast!='x')) continue;
		if (clast=='x') if (adj_used[i]>0) continue;
		// make sure the the NEXT adjective is not a suffix for this qty
		if ((i+1)<asz)
		{ // if the next item is a valid suffix skip this as qt 
			const StrTy & b=adjectives[i+1];
			const StrTy sfx=suffix(b);
			const D conv=sfx_map()[sfx];
			// is it is 'x' then this is still a count 
			if (sfx==b) if (sfx!="x") if (conv>0) continue; // 
		}
		++adj_used[i];	
		const D v=atof(s);
		if (explicit_count<0) explicit_count=v;
		else explicit_count+=v;	
		//if (db&&db2) MM_ERR(MMPR3(s,v,explicit_count))
	}
	// it COULD explicitly be ZERO 
	if (explicit_count>=0) return explicit_count;
//	if (db&&db2) MM_ERR(MMPR(explicit_count))
	// assume 1 if nother explicit 
	return 1; 
}
// see what was used and check for possible errors
virtual IdxTy  errcount(Used & adj_used, const D & upc, const D & qty, D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
	IdxTy er=0;
	const IdxTy asz=adjectives.size();
	for(IdxTy i=0; i<asz; ++i)
	{
		const StrTy & a=adjectives[i];
		if (safe_to_ignore(a)) continue;
		// this needs to distinguish adjectives specifying which type from 
		if (adj_used[i]!=1)
		{ // in some cases an adjective can be unused because it is
		// redundant with the explict size. Check if it is in the map
			const StrTy & a=adjectives[i];
			auto ii=m_qtu_map.find(a);
			if (ii!=m_qtu_map.end()) { if (upc== (*ii).second) continue;  } 
		// a suffix means a size 
			Ss ss; ss<<" adj "<<a<<" used "<<adj_used[i];
			exp_log(__FILE__,__LINE__,"badadjuse",
				ss.str().c_str(),w,d,adjectives,dmel);
			continue;
		}
	}
return er;
}

virtual void qty_from_adj2(D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
// default impl
// qty_from_adj(cname,val,w,d,adjectives,dmel);
	const IdxTy asz=adjectives.size();
	Used adj_used(asz);
	IdxTy aused=0;
	D val_added=0;
	D units_per_count=qty_to_units2(adj_used,val,cname,w,d,adjectives,lit,dmel);
	D quantity=count(adj_used,val,cname,w,d,adjectives,lit,dmel);
	IdxTy erc=errcount(adj_used,units_per_count,quantity,val,cname,w,d,adjectives,lit,dmel);
	val_added=units_per_count*quantity;
	val+=val_added;
	if (true) return; 
return ; 
} //  _adj2


void exp_log(const char * f, const IdxTy line, const char * prob, const char * msg, 
const StrTy & w, const StrTy & d, const Words & adjectives, Dmel * dmel) const
{
	std::cerr<<f<<":"<<line<<MMPR2(prob,msg)<<MMPR2(w,d)<<CRLF;
	if (dmel!=0) (*dmel).event(prob,w,d,adjectives);
	m_cm.inc_location(f,line,prob );
}

D sfx_convert(	 const StrTy & sfx)
{
	SfxMap & s=sfx_map();
	if (s.find(sfx)==s.end()) return -1;
	const D conv=sfx_map()[sfx];
	return conv;
}
SfxMap & sfx_map() // const StrTy & sfx)
{
static SfxMap sfxv;
static bool init=false;
if (!init)
{
sfxv["mg"]=1;
sfxv["g"]=1000;
sfxv["mcg"]=.001;
sfxv["x"]=0;
sfxv[""]=0;


init=true;
}
return sfxv;
}


#if 0 
// obsolete 
virtual void qty_from_adj(StrTy & cname, D& val,const StrTy & w
		,const StrTy & d, const Words & adjectives, Dmel * dmel=0)
{

// this should call qty_t_units 
D mg_per_qty=qty_to_units(w,d,adjectives,dmel); // 100; // need lut for w and maybe date 

const IdxTy asz=adjectives.size();
if (asz==0) { val=1; return ; }
for (IdxTy i=0; i<asz; ++i)
{
// this can not be used in the general case as literals
// have been gated out 
const StrTy & a=adjectives[i];
// in the other code this is included 
//if (m_literals[a]== ADJECTIVE ) continue;
//if (!parse_mg(a,val)) val+=atof(a.c_str());
const StrTy sfx=suffix(a);
const bool is_num=(sfx.length()==0);
const bool maybe_time=(a.length()==6);
D vall=q_x_sz(a);
if (is_num) { val+=vall*mg_per_qty; }
else if (maybe_time&&(sfx=="PM")) { val+=0;} 
else if (maybe_time&&(sfx=="AM")) { val+=0;} 
else if (sfx=="mg") { val+=vall;} 
else if (sfx=="x") { val+=mg_per_qty*atof(a.c_str());} 
else if (sfx=="mcg") { val+=vall*.001;} 
else if (sfx=="g") { val+=vall*1000.;} 
else { val+=atof(a.c_str()); { 
if (dmel!=0) (*dmel).event("unknown suffix ",w,d);

MM_ERR(" unknown suffix "<<a<<" sfx "<<sfx<<" date "<<d) }{  MM_INC_MSG(m_cm,"qtysuffox")}  } 

}
}

#endif

bool contains(const char * s, const char c)
{
IdxTy i=0; 
while (s[i]!=0) { if (s[i]==c) return true; ++i; }
return false;
}
bool containsdigit(const char * s)
{
IdxTy i=0; 
while (s[i]!=0) { if (digitt(s[i])) return true; ++i; }
return false;
}
bool containsnondigit(const char * s)
{
IdxTy i=0; 
while (s[i]!=0) { const char c=s[i]; if (!digitt(c)&&(c!='.')&&(c!='/')) return true; ++i; }
return false;
}
char firstnondigit(const char * s)
{
IdxTy i=0; 
while (s[i]!=0) { const char c=s[i]; if (!digitt(c)&&(c!='.')&&(c!='/')) return c; ++i; }
return 0;
}



bool digitt(const char c ) const { return (c>='0')&&(c<='9'); } 
bool lc(const char c ) const { return (c>='a')&&(c<='z'); } 
bool uc(const char c ) const { return (c>='A')&&(c<='Z'); } 
StrTy suffix(const StrTy & a)
{
const IdxTy sz=a.length();
const char * s=a.c_str();
StrTy sfx="";
int  i=sz-1;
while (i>=0) { if (digitt(s[i])) break; --i; } 
sfx=s+i+1;
return sfx;
}

D q_x_sz(const StrTy & a, const IdxTy lim=0)
{
const IdxTy sz=a.length();
if (sz==0) return 0;
D val=0;
D factor=1;
const char * s=a.c_str();
for (IdxTy i=0; i<sz; ++i)
{
if (lim!=0) if (i>=lim) break; // skip suffix for valir 
if (s[i]=='x')
{
const D f=atof(s);
const D v=atof(s+i+1);
val=f*v*factor;
return val;
} else if ((s[i]!='.')&&!digitt(s[i])) return -1;

}
val=factor*atof(s); // just return the value of numeric part 
return val;
}

bool parse_mg(const StrTy & a, D & val)
{
const IdxTy sz=a.length();
if (sz==0) return false;
const char * s=a.c_str();
D factor=1;
if (s[sz-1]!='g') return false;
if (s[sz-2]=='m') { factor=1;}
else { if (s[sz-2]!='c') { return false;}
 if (s[sz-3]!='m') { return false;}
factor=1e-3;
}
for (IdxTy i=0; i<sz; ++i)
{
if (s[i]=='x')
{
const D f=atof(s);
const D v=atof(s+i+1);
val+=f*v*factor;
return true;
}

}
val+=factor*atof(s); // just return the value of numeric part 
return true; 
}
void ignore(const StrTy & a) {m_ignores[a]=1; }
void ignore(const char *  a) {m_ignores[StrTy(a)]=1; }
private:
StrTy m_name;
D m_qtu;
BrandQty m_qtu_map;
IgnoreMap m_ignores;
mutable CounterMap m_cm;
bool m_ignore_share;

}; // drugs

*/



class mjm_dictionary
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



typedef mjm_dictionary Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

//typedef mjm_logic_base Logic;
typedef dict_params Logic;
typedef mjm_logic_base VariableStore;



typedef std::vector<StrTy> DumpOrder;

typedef literal_map LiteralClass;
//typedef LiteralClass::LITERAL_TYPE LITERAL_TYPE;
//typedef literal_map LiteralClass;

// collection of things for each date 
typedef std::map<StrTy, DateInfo > MappedInfo;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

/// need to have dog_day def available. 
// MACRO is indexed by name and then date 


public :
mjm_dictionary():m_dmel(new Dmel()) {Init();}
mjm_dictionary(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_snacks()
{
clear_handlers();
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
//if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
//if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { ++i; clean_up(); }
} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 

//void start_date(const StrTy & d) { m_flp.set("start_date",d); }
//void end_date(const StrTy & d) { m_flp.set("end_date",d); }

#if 0
void dump_unused()
{
Ss ss;
//for (auto ii=m_unused.begin(); ii!=m_unused.end(); ++ii)
{
//ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
MM_MSG("unused:"<<CRLF<<ss.str())

}
#endif

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
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="dump-dmel") { dump_dmel(std::cout); continue; } 
if (cmd=="dump-dmel-cerr") { dump_dmel(std::cerr); continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
//if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
#if 0
if (cmd=="parse") { parse();  continue; } 
if (cmd=="pair") { pair_stats(li.word(1),li.word(2),li.word(3));  continue; } 
if (cmd=="pair-all") { pair_stats(li.word(1));  continue; } 
if (cmd=="dose-vectors") { if (sz>1) { dose_vectors(li.word(1));} else { MM_ERR( "need dogname ") }  continue; } 
//void pair_stats(const StrTy & x, const StrTy & y)
if (cmd=="dump-drugs") { collect_all_drugs(); continue; } 
if (cmd=="dump-drugs-cerr") { collect_all_drugs(2); continue; } 
if (cmd=="list-drugs-cerr") { collect_all_drugs(3); continue; } 
if (cmd=="list-dogs-cerr") { collect_all_drugs(4); continue; } 
if (cmd=="dump-counts") { stream_counts(std::cout); continue; } 
if (cmd=="dump-canon-fails") { dump_canon_fails(std::cout); continue; } 
if (cmd=="clear-canon-fails") { clear_canon_fails(); continue; } 
if (cmd=="dump-dog-days") { dump_dog_days(std::cout,p1); continue; } 
if (cmd=="n-canon-fails") { MM_MSG(MMPR(m_canon_fails.size())) continue; } 

if (cmd=="clear-order") { clear_order(); continue; } 
if (cmd=="dump-order") { dump_order(std::cout); continue; } 
if (cmd=="push-order") 
		//{ int i=1; while ( li.cmd_ok(i+1)){ push_order(li.word(i)); ++i; } continue; } 
		{ int i=1; while ( i<sz){ push_order(li.word(i)); ++i; } continue; } 
if (cmd=="remove-order") 
		{ int i=1; while ( i<sz){ remove_order(li.word(i)); ++i; } continue; } 

if (cmd=="load-order") { if( li.cmd_ok(2)) load_order(li.word(1)); continue; } 
if (cmd=="order-all") { collect_all_drugs(1); continue; } 

#endif

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
if (cmd=="test") { test(li);  continue; } 
/*
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
*/
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrec0gnized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");

}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

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
//MM_MSG(" configuration "<<m_flp.to_string())
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}
bool done() const  { return m_done; } 

void clear_order()
{
m_order.clear();
}

#if 0 
void dump_order(OsTy & os)
{
for (IdxTy i=0; i<m_order.size(); ++i)
	os<<i<<" "<<m_order[i]<<CRLF;
}
void push_order(const StrTy & w)
{
if (w!="") m_order.push_back(w);
}
// TODO have this take a map of things to delete 
void remove_order(const StrTy & w)
{
const IdxTy sz=m_order.size();
Words wv(sz);
for (IdxTy i=0; i<sz; ++i) 
{
const StrTy & wo=m_order[i];
if (w!=wo ) wv.push_back(wo);
}
m_order=wv; 
}


void load_order(const StrTy & fn )
{
	std::ifstream  isn(fn.c_str());
	while (isn.good()&&!isn.eof()) 
		{ StrTy x=""; isn>>x; if (x.length()!=0) m_order.push_back(x); }
}
#endif

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general
// these should have been moved 
#if 0 
bool digitt(const char c ) const { return (c>='0')&&(c<='9'); } 
bool lc(const char c ) const { return (c>='a')&&(c<='z'); } 
bool uc(const char c ) const { return (c>='A')&&(c<='Z'); } 

bool date(const char * s) const 
{
if (strlen(s)!=10) return false;
if (s[4]!='-') return false;
if (s[7]!='-') return false;
if (!digitt(s[0])) return false;
if (!digitt(s[1])) return false;
if (!digitt(s[2])) return false;
if (!digitt(s[3])) return false;
return true;
}
bool dog_name(const char * s) const
{
const IdxTy sz=strlen(s);
if (sz<2) return false;
if (!uc(s[0])) return false;
for (IdxTy i=1; i<sz; ++i) if (!lc(s[i])) return false;
return true;
}
bool macro_name(const char * s) const
{
const IdxTy sz=strlen(s);
if (sz<2) return false;
//if (!uc(s[0])) return false;
for (IdxTy i=0; i<sz; ++i) if (!uc(s[i])&&!digitt(s[i])) return false;
return true;
}


///////////////////////////////////////////////////////////////////



// after the dates, the first word is either a dog name or
// a definition name or reserved command word
// MACRO contains no lower case, dog name is Normal capitalizaiton 
int parse_name(const StrTy & name)
{
const IdxTy sz=name.length();
if (sz<2) return INVALID;
const char * cstr=name.c_str();
if (name=="AKA") return ALIAS;
if (name=="NOTE") return NOTE;
if (name=="COMMENT") return NOTE;
if (dog_name(cstr)) return DOG_NAME;
if (date(cstr)) return DATE;
//if (cstr[2]>='Z') return MACRO;
if (sz==6) if (cstr[5]=='M') 
	if (digitt(cstr[0])&&digitt(cstr[1])&&digitt(cstr[2])&&digitt(cstr[3]))
		return TIME;
if (macro_name(cstr)) return MACRO;
return INVALID;
}




// determine if the given word is a macro ref possible modifed
bool macro_ref(const StrTy & w, const StrTy & d, const IdxTy lt )
{
if (lt==NONE)  return !notmacro(w.c_str());
if (lt==LMACRO) return true;
return false; 
}
// numeric in some form 
bool contains(const char * s, const char c)
{
IdxTy i=0; 
while (s[i]!=0) { if (s[i]==c) return true; ++i; }
return false;
}
bool containsdigit(const char * s)
{
IdxTy i=0; 
while (s[i]!=0) { if (digitt(s[i])) return true; ++i; }
return false;
}


// return true if this can NOT be a macro but
// false is not informaive 
bool notmacro(const char * s)
{
IdxTy i=0; 
bool found_minus=false;
bool found_uc=false;
while (s[i]!=0) { 
	if (s[i]=='-') found_minus=true;
	// modifiers allowed 
	if (!found_minus) if (lc(s[i])) return true; 
	if (uc(s[i])) found_uc= true; 
	if (s[i]=='.') return true;
++i; }
return !found_uc; // false;
}

bool adjective_ref(const StrTy & w, const StrTy & d,const IdxTy lt)
{
//enum LITERAL_TYPE { ADJECTIVE=1, QUANTIFIER=2, NOUN=3, LMACRO=5,IGNORE=4, NONE=~0 };
if (lt==ADJECTIVE) return true;
if (lt==QUANTIFIER) return true;
if (lt!=NONE) return !true;

const char * c=w.c_str();
const IdxTy sz=w.length();
if (sz==0) return false;
//if (m_literals.find(w)!=m_literals.end()) return true;
//if (c[0]='.') return true;
//if (contains(c,'.')) return true;
if (notmacro(c)) if (containsdigit(c)) return true;
//if (sz>2) if (c[sz-2]=='m') if (c[sz-1]=='g') if (digitt(c[sz-3])) return true;

return false;
}

// each macro gets parsed into a "dog_day" item 
dog_day macro_lut(const StrTy & w, const StrTy & d)
{
// find base name,
std::vector<StrTy> words;
//const char * c=w.c_str();

const IdxTy sz=w.length();
char c[sz+1];
memcpy(c,w.c_str(),sz);
c[sz]=0;
 IdxTy ptr=0;
for (IdxTy i=0; i<sz; ++i)
{
if (c[i]=='-')
{
c[i]=0;
words.push_back(StrTy(c+ptr));
ptr=i+1;
}

}
if (ptr<sz) words.push_back(StrTy(c+ptr));
if (words.size()==0) words.push_back(StrTy(""));
auto ii=m_macro.find(words[0]);
if (ii==m_macro.end())
{
MM_MSG(" no macro ever for "<<words[0]<<" from "<<w<<" "<<d)
MM_INC_MSG(m_cm, "nomacroever")
dog_day dd;
return dd;
}
auto jj=(*ii).second.find(d);
if (jj==(*ii).second.end())
{

MM_MSG(" no macro  for "<<words[0]<<" "<<d<<" from "<<w<<" "<<d)
MM_INC_MSG(m_cm, "nomacrodate")
dog_day dd;
return dd;

}
// so this returns a COPY of 
dog_day dd=m_macro[words[0]][d];
const IdxTy wsz=words.size();
// now subtract off all the entries with base names in minus list 
for (IdxTy i=1; i<wsz; ++i) { dd.set_dose(words[i],0); }

return dd; 

}

// for a given name w on date d, find the canonical name cname and value
// from adjective list 
void	qty_from_adj(StrTy & cname, D& val,const StrTy & w,const StrTy & d)
{
	if (m_canon.find(w)==m_canon.end()) 
	{ 
		exp_log(__FILE__,__LINE__,"addcanon","adding canonical self name for ",w,d);
		m_canon[w]=w;
		m_canon_fails.push_back(w);
		m_canon_fails_date.push_back(d);
	}	// always create some cname 
	cname=m_canon[w];
	Noun default_handler;
	Noun * h=& default_handler;
	auto ii=m_handlers.find(cname);
	if (ii!=m_handlers.end()) { h=(*ii).second;  }
	(*h).qty_from_adj2(val,cname,w, d, m_adjectives,m_literals,m_dmel);

}
#endif

void exp_log(const char * f, const IdxTy line, const char * prob, const char * msg, 
const StrTy & w, const StrTy & d)
{
	std::cerr<<f<<":"<<line<<MMPR4(f,line,prob,msg)<<MMPR2(w,d)<<CRLF;
	if (m_dmel!=0) (*m_dmel).event(prob,w,d,m_adjectives);
	m_cm.inc_location(f,line,prob );
}

#if 0 
void	qty_from_adj_old(StrTy & cname, D& val,const StrTy & w,const StrTy & d)
{
// look up the cname first and get any rules for making the quantity value
// specific to that 
const bool dump_qty_parse=false; // (w=="pantothenate");
if (m_canon.find(w)==m_canon.end()) 
{
MM_ERR(" adding canonical self name for "<<w<<" "<<d);
if (m_dmel!=0) (*m_dmel).event("addcanon",w,d);
MM_INC_MSG(m_cm,"addcanon" )
m_canon[w]=w;
m_canon_fails.push_back(w);
m_canon_fails_date.push_back(d);
}
cname=m_canon[w];
Noun * h=0;
auto ii=m_handlers.find(cname);
if (ii!=m_handlers.end())
{
if (dump_qty_parse) { MM_MSG(" found a handler for "<<cname) }
// it needs the adjective lst too 
//h=((*ii).second.qty_from_adj(cname,  val,w, d, m_adjectives));
h=(*ii).second; // .qty_from_adj(cname,  val,w, d, m_adjectives));
if ((*h).custom_qty()) {  (*h).qty_from_adj(cname,  val,w, d, m_adjectives); return; } 
}
// this does not work .. 
//static Noun dgeneric;
//return dgeneric.second.qty_from_adj(cname,  val,w, d);
// this is not needed and may FAIL unless there is a count instead of amount 
//D mg_per_qty=(h==0)?1:((*h).qty_to_units(w,d,m_adjectives,m_dmel)); // need lut for w and maybe date 
//if (dump_qty_parse) {MM_MSG(MMPR(mg_per_qty)) }
const IdxTy asz=m_adjectives.size();
IdxTy aused=0;
// see aused logic at bottom 
if (false) if (asz==0) { 
D mg_per_qty=(h==0)?1:((*h).qty_to_units(w,d,m_adjectives,m_dmel)); // need lut for w and maybe date 
if (dump_qty_parse) {MM_MSG(MMPR(mg_per_qty)) }

val+=mg_per_qty; return ; }

D val_added=0;
for (IdxTy i=0; i<asz; ++i)
{
const StrTy & a=m_adjectives[i];
// this needs to distinguish adjectives specifying which type from 
// quantifieers 
// FUDD THESE FUDDING NUM FUDDK SHOT 
//if (m_literals.map()[a]== ADJECTIVE ) 
if (m_literals.is_adjective(a)) 
{
	MM_ERR(" ignoring literal adjective "<<a<<" for "<<w<<" "<<d);
	if (m_dmel!=0) (*m_dmel).event("badlitadj",w,d,m_adjectives);
	MM_INC_MSG(m_cm,"badlitadj" )
	continue;
}
++aused;
//if (!parse_mg(a,val)) val+=atof(a.c_str());
const StrTy sfx=suffix(a);
const bool is_num=(sfx.length()==0);
D vall=q_x_sz(a);
if (is_num) { 
D mg_per_qty=(h==0)?1:((*h).qty_to_units(w,d,m_adjectives,m_dmel)); // need lut for w and maybe date 
if (dump_qty_parse) {MM_MSG(MMPR(mg_per_qty)) }
val_added+=vall*mg_per_qty;
 }
else if (sfx=="mg") { val_added+=vall;
} 
else if (sfx=="x") { 
D mg_per_qty=(h==0)?1:((*h).qty_to_units(w,d,m_adjectives,m_dmel)); // need lut for w and maybe date 
if (dump_qty_parse) {MM_MSG(MMPR(mg_per_qty)) }
val_added+=mg_per_qty*atof(a.c_str());

} 
else if (sfx=="mcg") { val_added+=vall*.001;} 
else if (sfx=="g") { val_added+=vall*1000.;} 
else { val_added+=atof(a.c_str()); { 
if (m_dmel!=0) (*m_dmel).event("badsuffix",a,d);
MM_ERR(" unknown suffix "<<a<<" sfx "<<sfx<<" date "<<d) }{  MM_INC_MSG(m_cm,"qtysuffox")}  } 

}
if (aused==0) { 
D mg_per_qty=(h==0)?1:((*h).qty_to_units(w,d,m_adjectives,m_dmel)); // need lut for w and maybe date 
if (dump_qty_parse) {MM_MSG(MMPR(mg_per_qty)) }

val+=mg_per_qty; return ; }

val+=val_added;
} // qty_from_adj_old

#endif


#if 0 
StrTy suffix(const StrTy & a)
{
const IdxTy sz=a.length();
const char * s=a.c_str();
StrTy sfx="";
int  i=sz-1;
while (i>=0) { if (digitt(s[i])) break; --i; } 
sfx=s+i+1;
return sfx;
}

D q_x_sz(const StrTy & a)
{
const IdxTy sz=a.length();
if (sz==0) return 0;
D val=0;
D factor=1;
const char * s=a.c_str();
for (IdxTy i=0; i<sz; ++i)
{
if (s[i]=='x')
{
const D f=atof(s);
const D v=atof(s+i+1);
val=f*v*factor;
return val;
}

}
val=factor*atof(s); // just return the value of numeric part 
return val;
}

bool parse_mg(const StrTy & a, D & val)
{
const IdxTy sz=a.length();
if (sz==0) return false;
const char * s=a.c_str();
D factor=1;
if (s[sz-1]!='g') return false;
if (s[sz-2]=='m') { factor=1;}
else { if (s[sz-2]!='c') { return false;}
 if (s[sz-3]!='m') { return false;}
factor=1e-3;
}
for (IdxTy i=0; i<sz; ++i)
{
if (s[i]=='x')
{
const D f=atof(s);
const D v=atof(s+i+1);
val+=f*v*factor;
return true;
}

}
val+=factor*atof(s); // just return the value of numeric part 
return true; 
}


#endif








// the amount of macro w shared by dog n on date d
void catalog_entries ( const int name_class, const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	for (IdxTy i=0; i<sz; ++i)
	{
		const StrTy f=li.word(i);
		const IdxTy j=(i<4)?i:4;
		m_literals.inc_word(j,name_class,f,name_class_name(name_class));
//		++m_words[j][name_class][f];
	} 

}
void clear_canon_fails()
{
m_canon_fails.clear();
m_canon_fails_date.clear();
}
void dump_canon_fails(OsTy & os)
{
const IdxTy sz=m_canon_fails.size();
for (IdxTy i=0; i<sz; ++i) os<<m_canon_fails[i]<<" "<<m_canon_fails_date[i]<<CRLF;

}

void stream_counts(OsTy & os)
{
	m_literals.stream_counts(os); 
}


void parse()
{
	//IsTy * is=&std::cin;
	std::ifstream  isn(m_flp.snack_log().c_str());
	IsTy * is=& isn; // &std::cin;
	const bool skip_old=m_flp.skip_old();
	const bool print_dog_days=m_flp.print_dog_days();
	const bool accumulate_dog_days=m_flp.accumulate_dog_days();
	const StrTy start_date=m_flp.start_date();
	const StrTy end_date=m_flp.end_date();
	bool flush_olds=!true;
	bool skipping=skip_old;
	StrTy last_date="";
//	MM_ERR(MMPR2(skip_old,skipping))
	// LineIterator would work too
	CommandInterpretter li(is);
	while (li.nextok())
	{
	//	catalog_entries(li);
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<4) continue;
		const StrTy date=li.word(0);
		if (!this->date(date.c_str()))
		{ // TODO FIXME use the err_log thing 
			MM_INC_MSG(m_cm,"baddate" )
			if (m_dmel!=0) (*m_dmel).event("baddate",date,date);
			{MM_ERR(" bad date on line "<<li.line()); }
		}
		if (last_date!=date)
		{
			end_of_day_check();
			if (flush_olds) clear_date_cache();
			last_date=date;	
		}
		if (skipping)
		{
//	MM_ERR(MMPR3(skip_old,skipping, date))
			//if (skip_old) if (date==StrTy("2017-04-22")) skipping=false;
			if (skip_old) if (date>=start_date) skipping=false;
			if (skipping) continue;
		}
		if (date>end_date) break;
		//MM_MSG(li.line());
		const StrTy  amper=li.word(1); // should be ampersand
		if (amper!=StrTy("&")) {
	// TODO FIXME use the exp_log thing 
MM_INC_MSG(m_cm,"badline" )
if (m_dmel!=0) (*m_dmel).event("missingamp",amper,date);

MM_ERR(" bad line "<<li.line()); }
		const StrTy name=li.word(2);
		// this is either a dog name, or a macro for a meal or a NOTE or EVENT or AKA
		const int name_class=parse_name(name);
		catalog_entries(name_class, li);
//enum NAME_CLASS { DOG_NAME=0, ALIAS=1, MACRO=2, TIME=3, DATE=4, NOTE=5, INVALID=~0 };
		switch (name_class)
		{
			case DOG_NAME : { dog_day dd=parse_dog_day(date,name,3,li);
				if (false)	MM_MSG(dd.string())
				if (print_dog_days)	MM_MSG(dd.string(m_order))
				if (accumulate_dog_days){ m_dog_days[date][name]=dd;  }	
 				break; }
			case MACRO : { parse_macro(date,name,3,li); break; }
			case NOTE : { parse_note(date,name,3,li); break; }
			case ALIAS: { break; } // this is defined as TEXT not like MACRO 
			default:
			{ 
			// TODO FIXME use the exp_log 
				MM_INC_MSG(m_cm,"badclassparse" )
				if (m_dmel!=0) (*m_dmel).event("badclassparse",name,date);
				MM_ERR(" bad name_class "<<MMPR4(name_class,name_class_name(name_class),name,date)) }
		} // name_class
	}
	end_of_day_check();

}
void end_of_day_check()
{
// make sure macro usages make sense 

}
// TODO FIXME this really needs an iteroatr 
// better off just to make itor 
template <typename Ty > void traverse_all_dog_days(Ty & f)
{
Words  w(3);
const DogDays & dd= m_dog_days;
for(auto ii=dd.begin(); ii!=dd.end(); ++ii)
{
const auto & day=(*ii).first;
w[0]=day;
const auto & dm=(*ii).second;
for(auto jj=dm.begin(); jj!=dm.end(); ++jj)
{
const auto & dog=(*jj).first;
w[1]=dog;
const auto & dom=(*jj).second.m_doses;
for(auto kk=dom.begin(); kk!=dom.end(); ++kk)
{
const StrTy drug=(*kk).first;
w[2]=drug;
const D dose=(*kk).second;
f(w,dose);
} // kk

} // jj

} // ii 

}

// give each dose combination a vector number in order of
// first usage making them indications of recency novelty etc
void dose_vectors( const StrTy & thedog)
{
const bool dump_ts=true;
const bool dump_h=true;
const bool dump_freq=true;

Words w(3);
typedef std::map<StrTy, IdxTy >  HMap;
HMap hmap,vmap,tsmap;
IdxTy n=0;
const IdxTy sz= m_order.size();
const DogDays & dd= m_dog_days;
for(auto ii=dd.begin(); ii!=dd.end(); ++ii)
{
const auto & day=(*ii).first;
w[0]=day;
const auto & dm=(*ii).second;
for(auto jj=dm.begin(); jj!=dm.end(); ++jj)
{
const auto & dog=(*jj).first;
if (dog!=thedog) continue;
w[1]=dog;
const auto & dom=(*jj).second;
Ss ssx;
for (IdxTy i=0; i<sz; ++i)
{
//const StrTy drug=(*kk).first;
const StrTy & key=m_order[i];
const auto kkx=dom.m_doses.find(key);
bool dxe=(kkx==dom.m_doses.end()) ;
const D dosex=dxe?0:(*kkx).second;
if (i!=0) ssx<<",";
ssx<<key<<"="<<dosex; 
//f(w,dose);
} // for 
const StrTy & k=ssx.str();
if (vmap.find(k)==vmap.end()) vmap[k]=vmap.size()+1;
tsmap[w[0]]=vmap[k];
++hmap[k];
++n;
} // jj
} // ii 
if (dump_ts){
for(auto ii=tsmap.begin(); ii!=tsmap.end(); ++ii)
{
const StrTy & d=(*ii).first;
const IdxTy cw=(*ii).second;
MM_MSG(MMPR2(d,cw))
}
}
D h=0;
if (n==0) return;
D den=1.0/n;
D scale=1.0/log(2.0);
IdxTy c=0;
for(auto ii=hmap.begin(); ii!=hmap.end(); ++ii)
{
const auto & cw=(*ii).first;
//w[0]=day;
const auto & dm=(*ii).second;
const D p=dm*den;
const D dh=p*log(p)*scale;
h-=dh;
if (dump_freq) { MM_MSG("hmapv "<<MMPR4(cw,dm,dh,vmap[cw])) } 
++c;
} // ii 
if (dump_h) { MM_MSG("hmap disorder "<<MMPR3(h,n,c)) } 
} // dose_vectors

void pair_stats( const StrTy & thedog)
{
const auto & v=m_order;
const IdxTy sz=v.size();
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & x=v[i];
for(IdxTy j=i; j<sz; ++j)
{
const StrTy & y =v[j];
pair_stats( x, y, thedog);
} // j
} // i

}
void pair_stats(const StrTy & x, const StrTy & y, const StrTy & thedog)
{
Words  w(3);
typedef std::map<StrTy, std::map< StrTy, IdxTy > >  HMap;
HMap hmap;
IdxTy n=0;
const DogDays & dd= m_dog_days;
for(auto ii=dd.begin(); ii!=dd.end(); ++ii)
{
const auto & day=(*ii).first;
w[0]=day;
const auto & dm=(*ii).second;
for(auto jj=dm.begin(); jj!=dm.end(); ++jj)
{
const auto & dog=(*jj).first;
if (dog!=thedog) continue;
w[1]=dog;
const auto & dom=(*jj).second.m_doses;
const auto kkx=dom.find(x);
const auto kky=dom.find(y);
//if (kkx==dom.end()) continue;
bool dxe=(kkx==dom.end()) ;
//if (kky==dom.end()) continue;
bool dye=(kky==dom.end()) ;

Ss ssx, ssy;
//const StrTy drug=(*kk).first;
const D dosex=dxe?0:(*kkx).second;
const D dosey=dye?0:(*kky).second;
ssx<<dosex; ssy<<dosey;
++hmap[ssx.str()][ssy.str()];
//f(w,dose);
++n;
} // jj

} // ii 
D h=0;
if (n==0) return;
D den=1.0/n;
D scale=1.0/log(2.0);
IdxTy c=0;
for(auto ii=hmap.begin(); ii!=hmap.end(); ++ii)
{
const auto & dosex=(*ii).first;
//w[0]=day;
const auto & dm=(*ii).second;
for(auto jj=dm.begin(); jj!=dm.end(); ++jj)
{
const auto & dosey=(*jj).first;
const auto & cnt=(*jj).second;
const D p=D(cnt)*den;
MM_MSG(MMPR4(dosex,dosey,cnt,p)<<MMPR(n))
h=h-p*log(p)*scale;
++c;
} // jj 
} // ii 
MM_MSG(MMPR4(x,y,h,c)<<MMPR(n))
}





void collect_all_drugs(const IdxTy flags=0)
{
typedef std::map<StrTy, IdxTy> AllDrugs;
AllDrugs ad;
AllDrugs adogs;
Ss ss;
//typedef std::map< StrTy,  std::map< StrTy, dog_day> > DogDays;
const DogDays & dd= m_dog_days;
for(auto ii=dd.begin(); ii!=dd.end(); ++ii)
{
const auto & day=(*ii).first;
const auto & dm=(*ii).second;
for(auto jj=dm.begin(); jj!=dm.end(); ++jj)
{
const auto & dog=(*jj).first;
++adogs[dog];
const auto & dom=(*jj).second.m_doses;
for(auto kk=dom.begin(); kk!=dom.end(); ++kk)
{
if ((*kk).second!=0) ++ad[(*kk).first];
} // kk

} // jj 
} // ii 
//if (flags==0) {
IdxTy i=0; 
for(auto ii=ad.begin(); ii!=ad.end(); ++ii)
{
const StrTy key=(*ii).first;
const IdxTy count=(*ii).second;
switch (flags)
{
case 1:{ m_order.push_back(key); break ; } 
case 2:{ MM_ERR(MMPR3(i,key,count)) break; } 
case 3:{ ss<<" "<<key; break; } 
case 4:{ break ; } 
default:
MM_MSG(MMPR3(i,key,count))
}

++i;
} // ii 
if (flags==3) { MM_ERR(ss.str())}
if (flags==4) { 
IdxTy i=0;
for (auto ii=adogs.begin(); ii!=adogs.end(); ++ii)
{
const StrTy dog=(*ii).first;
const IdxTy count=(*ii).second;
MM_ERR(MMPR3(i,dog,count))
++i;
}

}
return;
//} // flags==0


}


void load_literals()
{
auto & litmap=m_literals.map();
auto & lit=m_literals;

lit.add_adjective("NOW");
lit.add_adjective("GNC");
lit.add_adjective("sundown");
// ??? TODO WTF 
lit.add_adjective("adjective");
lit.add_adjective("delta");
lit.add_adjective("extra");
// ignore?
//lit.add_adjective("and");
lit.add_adjective("small");
lit.add_adjective("some");
lit.add_adjective("wild");
lit.add_adjective("keta");
lit.add_adjective("walmart");
lit.add_adjective("atlantic");
lit.add_adjective("unknown");
lit.add_adjective("bulk");

/*
m_literals[StrTy("NOW")]=ADJECTIVE;
m_literals[StrTy("GNC")]=ADJECTIVE;
m_literals[StrTy("sundown")]=ADJECTIVE;
m_literals[StrTy("adjective")]=ADJECTIVE;
m_literals[StrTy("delta")]=ADJECTIVE;
m_literals[StrTy("extra")]=ADJECTIVE;
m_literals[StrTy("and")]=ADJECTIVE;
m_literals[StrTy("small")]=ADJECTIVE;
m_literals[StrTy("some")]=ADJECTIVE;
m_literals[StrTy("wild")]=ADJECTIVE;
m_literals[StrTy("keta")]=ADJECTIVE;
m_literals[StrTy("walmart")]=ADJECTIVE;
m_literals[StrTy("atlantic")]=ADJECTIVE;
m_literals[StrTy("unknown")]=ADJECTIVE;
m_literals[StrTy("bulk")]=ADJECTIVE;
*/

/*
m_literals[StrTy("etc")]=IGNORE;
m_literals[StrTy("")]=IGNORE;
// only appear in PHRASES
m_literals[StrTy("Diamond")]=IGNORE;
m_literals[StrTy("Naturals")]=IGNORE;
m_literals[StrTy("Extreme")]=IGNORE;
m_literals[StrTy("Athlete")]=IGNORE;
m_literals[StrTy("Fromm")]=IGNORE;
m_literals[StrTy("Puppy")]=IGNORE;
m_literals[StrTy("Gold")]=IGNORE;


m_literals[StrTy("Flint")]=IGNORE;
m_literals[StrTy("River")]=IGNORE;
m_literals[StrTy("Ranch")]=IGNORE;
m_literals[StrTy("Salmon")]=IGNORE;
m_literals[StrTy("And")]=IGNORE;
m_literals[StrTy("Potato")]=IGNORE;
m_literals[StrTy("Adult")]=IGNORE;

m_literals[StrTy("Chicken")]=IGNORE;
m_literals[StrTy("Rice")]=IGNORE;
m_literals[StrTy("American")]=IGNORE;
m_literals[StrTy("Journey")]=IGNORE;
m_literals[StrTy("Lamb")]=IGNORE;
m_literals[StrTy("Sweet")]=IGNORE;
m_literals[StrTy("Wellness")]=IGNORE;

m_literals[StrTy("nopatothenate")]=IGNORE;
m_literals[StrTy("nospinach")]=IGNORE;

// error report  awked 
m_literals[StrTy("noshrimp")]=IGNORE;
m_literals[StrTy("noCu")]=IGNORE;
m_literals[StrTy("nolecithin")]=IGNORE;
m_literals[StrTy("noMgcitrate")]=IGNORE;
m_literals[StrTy("nod-serine")]=IGNORE;
m_literals[StrTy("nolipoicacid")]=IGNORE;
m_literals[StrTy("nochili")]=IGNORE;
m_literals[StrTy("nonwkelp")]=IGNORE;
m_literals[StrTy("nosalmonbroth")]=IGNORE;
m_literals[StrTy("noarginine")]=IGNORE;
m_literals[StrTy("nooliveoil")]=IGNORE;
m_literals[StrTy("nopantothenic")]=IGNORE;
m_literals[StrTy("nobiotin")]=IGNORE;
m_literals[StrTy("noFleaAway")]=IGNORE;
m_literals[StrTy("noctbrothbs")]=IGNORE;
//m_literals[StrTy("unknown")]=IGNORE;
m_literals[StrTy("nopantothenate")]=IGNORE;
m_literals[StrTy("nogarlic")]=IGNORE;
m_literals[StrTy("noFleAway")]=IGNORE;
*/

lit.add_ignore("etc");
lit.add_ignore("");
lit.add_ignore("Diamond");
lit.add_ignore("Naturals");
lit.add_ignore("Extreme");
lit.add_ignore("Athlete");
lit.add_ignore("Fromm");
lit.add_ignore("Puppy");
lit.add_ignore("Gold");
lit.add_ignore("Flint");
lit.add_ignore("River");
lit.add_ignore("Ranch");
lit.add_ignore("Salmon");
lit.add_ignore("And");
lit.add_ignore("and");
lit.add_ignore("Potato");
lit.add_ignore("Adult");
lit.add_ignore("Chicken");
lit.add_ignore("Rice");
lit.add_ignore("American");
lit.add_ignore("Journey");
lit.add_ignore("Lamb");
lit.add_ignore("Sweet");
lit.add_ignore("Wellness");
lit.add_ignore("nopatothenate");
lit.add_ignore("nospinach");
lit.add_ignore("noshrimp");
lit.add_ignore("noCu");
lit.add_ignore("nolecithin");
lit.add_ignore("noMgcitrate");
lit.add_ignore("nod-serine");
lit.add_ignore("nolipoicacid");
lit.add_ignore("nochili");
lit.add_ignore("nonwkelp");
lit.add_ignore("nosalmonbroth");
lit.add_ignore("noarginine");
lit.add_ignore("nooliveoil");
lit.add_ignore("nopantothenic");
lit.add_ignore("nobiotin");
lit.add_ignore("noFleaAway");
lit.add_ignore("noctbrothbs");
// bad syntax modified wrong thing in source 
lit.add_ignore("unknown");
lit.add_ignore("nopantothenate");
lit.add_ignore("nogarlic");
lit.add_ignore("noFleAway");





// each key maps to a more limited set of canonical name synonym

m_canon[StrTy("Mgcitrate")]=StrTy("Mg");
m_canon[StrTy("Ivermectin")]=StrTy("Ivermectin");
m_canon[StrTy("Mn")]=StrTy("Mn");
m_canon[StrTy("LESI")]=StrTy("Iodine");
m_canon[StrTy("nwkelp")]=StrTy("Iodine");
m_canon[StrTy("kelp")]=StrTy("Iodine");
m_canon[StrTy("Cu")]=StrTy("Cu");
m_canon[StrTy("HMCU")]=StrTy("HMCU");
m_canon[StrTy("HMCu")]=StrTy("HMCU");
m_canon[StrTy("K2")]=StrTy("K2");
m_canon[StrTy("noK2")]=StrTy("noK2");
m_canon[StrTy("11KC")]=StrTy("11KC");
m_canon[StrTy("01KC")]=StrTy("01KC");
m_canon[StrTy("KCl")]=StrTy("10KC");
m_canon[StrTy("10KC")]=StrTy("10KC");
// these look like macro modified if not definee although hypen number is not
// going to exist although macro-macro is possible 
m_canon[StrTy("B-1")]=StrTy("B-1");
m_canon[StrTy("B-2")]=StrTy("B-2");
m_canon[StrTy("B-3")]=StrTy("B-3");
//m_canon[StrTy("SNB3")]=StrTy("B-3");
m_canon[StrTy("SNB-3")]=StrTy("B-3");
//m_canon[StrTy("NSB-3")]=StrTy("B-3");
m_canon[StrTy("B-6")]=StrTy("B-6");
m_canon[StrTy("B-12")]=StrTy("B-12");
//m_canon[StrTy("B-12")]=StrTy("B-12");
m_canon[StrTy("B-100")]=StrTy("B-100");
m_canon[StrTy("garlic")]=StrTy("garlic");
m_canon[StrTy("b20ngnc")]=StrTy("b20ngnc");
m_canon[StrTy("nob20ngnc")]=StrTy("nob20ngnc");
m_canon[StrTy("b15ngnc")]=StrTy("b15ngnc");
m_canon[StrTy("nob15ngnc")]=StrTy("nob15ngnc");
m_canon[StrTy("carrot")]=StrTy("carrot");
m_canon[StrTy("carrots")]=StrTy("carrot");
m_canon[StrTy("biotin")]=StrTy("biotin");
m_canon[StrTy("lipoicacid")]=StrTy("lipoicacid");
m_canon[StrTy("pantothenate")]=StrTy("pantothenate");
m_canon[StrTy("metasilicate")]=StrTy("metasilicate");
m_canon[StrTy("thiosulfate")]=StrTy("thiosulfate");
m_canon[StrTy("arginine")]=StrTy("arginine");
m_canon[StrTy("d-serine")]=StrTy("d-serine");
m_canon[StrTy("SnAgOx")]=StrTy("SnAgOx");
m_canon[StrTy("salmonbrothbs")]=StrTy("salmon");
m_canon[StrTy("ketasalmonbrothbs")]=StrTy("salmon");
m_canon[StrTy("salmonbroths")]=StrTy("salmon");
m_canon[StrTy("salmonbrothb")]=StrTy("salmon");
m_canon[StrTy("salmonbroth")]=StrTy("salmon");
m_canon[StrTy("salmon")]=StrTy("salmon");
m_canon[StrTy("wildsalmon")]=StrTy("salmon");
m_canon[StrTy("ctbroth")]=StrTy("ctbroth");
m_canon[StrTy("ctbrothbs")]=StrTy("ctbroth");
m_canon[StrTy("ctbroths")]=StrTy("ctbroth");
m_canon[StrTy("ctbrothb")]=StrTy("ctbroth");
m_canon[StrTy("ctskin")]=StrTy("ctskin");
m_canon[StrTy("coffee")]=StrTy("coffee");
m_canon[StrTy("oliveoil")]=StrTy("oliveoil");


m_canon[StrTy("spinach")]=StrTy("spinach");
m_canon[StrTy("chili")]=StrTy("chili");
m_canon[StrTy("Kibble")]=StrTy("kibble");
m_canon[StrTy("kibble")]=StrTy("kibble");
m_canon[StrTy("KibbleMixDiamondJourney")]=StrTy("kibble");
m_canon[StrTy("shrimp")]=StrTy("shrimp");
m_canon[StrTy("FleaAway")]=StrTy("FleaAway");
m_canon[StrTy("taurine")]=StrTy("taurine");
m_canon[StrTy("swings")]=StrTy("swings");
m_canon[StrTy("eggo3")]=StrTy("eggo3");
m_canon[StrTy("egg")]=StrTy("eggo3"); // for now ignore o3 status
m_canon[StrTy("noeggo3")]=StrTy("noeggo3");
m_canon[StrTy("yeggo3")]=StrTy("yeggo3");
m_canon[StrTy("weggo3")]=StrTy("weggo3");
m_canon[StrTy("lecithin")]=StrTy("lecithin");
m_canon[StrTy("pantothenic")]=StrTy("pantothenate"); // anachronism 
// this needs to be a type that checkes for over writes 
//m_canon[StrTy("nwkelp")]=StrTy("nwkelp");
m_canon[StrTy("worm")]=StrTy("worm");
//m_canon[StrTy("mgcitrate")]=StrTy("mgcitrate");


add_canonicals_to_literals();
}
void add_canonicals_to_literals()
{
// FUDD these FUDDING enums FUDD 
for (auto ii=m_canon.begin(); ii!=m_canon.end(); ++ii)
	//m_literals.map()[(*ii).first]=NOUN;
	m_literals.add_noun((*ii).first);


}
void load_handlers()
{
// TODO FIXME danger will robinsom memory leak esp if assigning
// over existing one. Need to subclass and delete existing  
clear_handlers();
m_handlers["pantothenate"] = new Noun("pantothenate",500);
m_handlers["biotin"] = new Noun("biotin",5);
m_handlers["taurine"] = new Noun("taurine",1000);
m_handlers["arginine"] = new Noun("arginine",500);
m_handlers["Mg"] = new Noun("Mg",200);
m_handlers["B-1"] = new Noun("B-1",300);
m_handlers["B-2"] = new b2_handler(); // new Noun("B-6",100);
//m_handlers["B-2"] = new Noun("B-2",100);
m_handlers["B-3"] = new Noun("B-3",100);
//m_handlers["B-3"] = new Noun("NSB-3",100);
//m_handlers["B-3"] = new Noun("SNB-3",100);
m_handlers["B-6"] = new b6_handler(); // new Noun("B-6",100);
m_handlers["B-12"] = new b12_handler(); // new Noun("B-6",100);
m_handlers["Mn"] = new mn_handler(); // new Noun("B-6",100);
m_handlers["lecithin"] = new lecithin_handler(); // new Noun("B-6",100);
m_handlers["10KC"] = new kcl_handler(); // new Noun("B-6",100);
m_handlers["salmon"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["shrimp"] = new shrimp_handler(); // new Noun("B-6",100);
m_handlers["d-serine"] = new Noun("d-serine",700);

m_handlers["Iodine"] = new kelp_handler();
m_handlers["lipoicacid"] = new lipoicacid_handler();
m_handlers["metasilicate"] = new metasilicate_handler();
m_handlers["thiosulfate"] = new thiosulfate_handler();
//m_handlers["Iodine"] = new Noun("Iodine",.300);
//m_handlers["Iodine"]->add_potency("LESI",1.0);
//m_handlers["Iodine"]->add_potency("nwkelp",.3);
//m_handlers["Iodine"]->add_potency("kelp",.325);

}
void clear_handlers()
{
for ( auto ii=m_handlers.begin(); ii!=m_handlers.end(); ++ii)
{ // do not change teh map yet lol 
	delete (*ii).second;
}	
m_handlers.clear();
}
void load_shares()
{
m_share_factor[StrTy("Spicey")]=.05;
m_share_factor[StrTy("Hershey")]=.05;
m_share_factor[StrTy("Dexter")]=.05;
m_share_factor[StrTy("Peapod")]=.1;
m_share_factor[StrTy("Beauty")]=.25;
m_share_factor[StrTy("Moe")]=.25;
m_share_factor[StrTy("Greta")]=.25;


m_share_factor[StrTy("Spicey")]=.05;
m_share_factor[StrTy("Hershey")]=.05;
m_share_factor[StrTy("Dexter")]=.05;
m_share_factor[StrTy("Peapod")]=.2;
m_share_factor[StrTy("Beauty")]=.25;
m_share_factor[StrTy("Moe")]=.2;
m_share_factor[StrTy("Greta")]=.2;


ShareFactor & sf=m_share_factors[StrTy("PPHOME")];
sf[StrTy("Spicey")]=.0;
sf[StrTy("Hershey")]=.0;
sf[StrTy("Dexter")]=.0;
sf[StrTy("Peapod")]=1.;
sf[StrTy("Beauty")]=.0;
sf[StrTy("Moe")]=.0;
sf[StrTy("Greta")]=.0;

}
// things specific to older dates 
void  clear_date_cache()
{
m_adjectives.clear(); // should be line specific
m_macro.clear(); // old macros except alias ignored 

}

/////////////////////////////////////////////////////////////////////
private:
void Init()
{
//delete m_dmel;
//m_dmel= new Dmel();
m_done=false;
load_literals();
load_shares();
load_handlers();
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
// make exhaustive list of everything found in different places for QC
//FreqClassInfo m_words[5];
// push back adjectives until discovering what they modify 
Words m_adjectives;
MacroMap m_macro;
//Literals m_literals;
LiteralClass m_literals;
Canonical m_canon;
CanonFails m_canon_fails,m_canon_fails_date;
DumpOrder m_order;
ShareFactor m_share_factor;
ShareFactorMap m_share_factors;
NounHandler m_handlers;
DogDays m_dog_days;
CounterMap m_cm;

}; //mjm_timeline 



/////////////////////////////////////////////////////////

#ifdef  TEST_SNACK__
int main(int argc,char **args)
{
typedef mjm_snacks  Myt;
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

