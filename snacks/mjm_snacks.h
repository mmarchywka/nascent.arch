#ifndef MJM_SNACKS_H__
#define MJM_SNACKS_H__

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

/*
TODO FIXME probably should parse all modifiers first instead
of lr . Or, even entire dog day entry but modifier precedes
noun ias about only rule lol
TODO FIXME need to map modifiers-noun to cname then distinguish cname
from an aggregation name or aname for output to user. cname
handles common issues like trival quantitiy 
quirks (see 3x4150 for all shrimp lol).

Parse semi stuctured manual input log entries to get total dose info
as function of day. This relies on some prepartins,


~/d/latex$ cat cases.tex | test_timeline -date-lines > $cpp/mjm_libmesh/snacks_log

creating entires like this,

2017-09-06 & Greta AMSNACK 0930AMSNACK 500mg pantothenate  
2017-09-06 & Beauty AMSNACK 500mg pantothenate 0930AMSNACK 500mg pantothete  
2017-09-06 & Moe AMSNACK 0930AMSNACK  
2017-09-06 & Hershey AMSNACK  
2017-09-06 & Peapod AMSNACK .5 B-1 2mg Cu 500mg pantothenate 0930AMSNACK  
2017-09-06 & Spicey AMSNACK  




TODO FIXME
*/

/*
g++ -DTEST_SNACK__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_snacks.h
*/

////////////////////////////////////////////////////////////////

class snack_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
snack_params( const StrTy & nm) : Super(nm) {}
snack_params() : Super() {}
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
bool print_conflicts() const { return m_map.get_bool("print_conflicts",!true); }

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
ss<<"print_conflicts"<<print_conflicts()<<sep;
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

class dog_day
{
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef double D;
typedef unsigned int IdxTy;
typedef std::map<StrTy, D > DoseMap;
typedef std::map<StrTy, DoseMap > TimeDoseMap;
typedef std::vector<StrTy> Def;
public:
dog_day():m_parse_tod(false),m_open(false)  {} // only for stl crap 
dog_day(const StrTy & n, const StrTy & d) : m_date(d),m_name(n),m_parse_tod(false),m_open(false) {}
// canonical name and amount 
void inc(const StrTy & k, const D & v) { m_doses[k]+=v;if (m_parse_tod) m_times[m_tod][k]+=v;  }  
void set_dose(const StrTy & n,const D & v ) {m_doses[n]=v;}
void set_tod(const StrTy & tod ) {m_tod=tod;m_parse_tod=true;}
void reparse( const Def & words)
{
m_def=words; m_open=true; 
MM_ERR(" found open macro "<<MMPR2(m_date,m_name))
}
const Def & def() const { return  m_def;}
const bool open() const { return m_open; } 
void add(const dog_day & that, const D & fac)
{
// FIXME the input file semantics may not reflect this 
if (m_parse_tod) {if (that.m_tod>m_tod) m_tod=that.m_tod; }
for (auto ii=that.m_doses.begin(); ii!=that.m_doses.end(); ++ii)
	{
auto key=(*ii).first;
m_doses[key]+=fac*(*ii).second;
// really adopting that tod is quesitonable TODO 
if (m_parse_tod) { m_times[m_tod][key]+=fac*(*ii).second;}
//MM_MSG(MMPR4(m_name,key,fac,(*ii).second)<<MMPR2(that.m_name,m_doses[key]))
}
}
template <class Ty> void add(const dog_day & that, const D & fac, Ty & handlers)
{
if (m_parse_tod) {if (that.m_tod>m_tod) m_tod=that.m_tod; }
for (auto ii=that.m_doses.begin(); ii!=that.m_doses.end(); ++ii)
	{
auto key=(*ii).first;
auto hh=handlers.find(key);
if (hh==handlers.end()) {
 m_doses[key]+=fac*(*ii).second; 
if (m_parse_tod)  m_times[m_tod][key]+=fac*(*ii).second; 
continue; }
//MM_MSG(MMPR4(m_name,key,fac,(*ii).second)<<MMPR2(that.m_name,m_doses[key]))
const auto * ha=(*hh).second;
//void share_doses(D & d, const D & f, const D & amt) { d+=f*amt; }
ha->share_doses(m_doses[key],fac,(*ii).second);
if (m_parse_tod) ha->share_doses(m_times[m_tod][key],fac,(*ii).second);
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

} // string

 StrTy time_string( ) const
{
const StrTy sep=" ";
Ss ss;
//if (true ) return ss.str(); 
for(auto ii=m_times.begin(); ii!=m_times.end(); ++ii)
{
const StrTy & time=(*ii).first;
const auto & doses=(*ii).second;
for(auto jj=doses.begin(); jj!=doses.end(); ++jj)
{
ss<<m_date<<sep<<time<<sep;
ss<<m_name;
ss<<sep<<(*jj).first<<sep<<(*jj).second<<CRLF;
}
}

return ss.str();

} // string

 template <class Tomap > StrTy time_string(const Tomap & om, const IdxTy & flags=0 ) const
{
const StrTy sep=" ";
Ss ss;
const bool conflict_count=((flags&1)!=0);
//if (true ) return ss.str(); 
for(auto ii=m_times.begin(); ii!=m_times.end(); ++ii)
{
const StrTy & time=(*ii).first;
const auto & doses=(*ii).second;
IdxTy count=0;
for(auto jj=doses.begin(); jj!=doses.end(); ++jj)
{
if (om.find((*jj).first)==om.end()) continue; 
++count;
ss<<m_date<<sep<<time<<sep;
ss<<m_name;
ss<<sep<<(*jj).first<<sep<<(*jj).second<<CRLF;
}
if (conflict_count) if (count>1) { ss<<m_date<<sep<<time<<sep<<m_name<<sep<<"conflict_count"<<sep<<count<<CRLF;}
}

return ss.str();

} // string


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
TimeDoseMap m_times;
StrTy m_tod;
bool m_parse_tod;
bool m_open;
Def m_def;

}; // dog_day

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




class adjectives
{


}; // adjectives
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
drugs(): m_name("none"),m_qtu(1),m_ignore_adj(false),m_ignore_share(false) {}
drugs(const StrTy & n ): m_name(n),m_qtu(1)
		,m_ignore_adj(false),m_ignore_share(false) {}
drugs(const D & qtu): m_name("none"),m_qtu(qtu)
		,m_ignore_share(false) {}
drugs(const StrTy & n, const D & qtu): m_name(n)
		,m_qtu(qtu),m_ignore_adj(false),m_ignore_share(false) {}
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
virtual void ignore_adj(const bool tf) {m_ignore_adj=tf; }
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
	if (m_ignore_adj) return 1; 
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
	if (m_ignore_adj) return 1; 
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
		if (m_ignore_adj) continue; // obviously this whole loop can be skpped, TODO 
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
/*	for (IdxTy i=0; i<asz; ++i)
	{
		const StrTy & a=adjectives[i];
		if (safe_to_ignore(a)) continue;
		// this needs to distinguish adjectives specifying which type from 
		if (false) if (lit.is_adjective(a)) 
		{
			exp_log(__FILE__,__LINE__,"badlitadj",
				" ignoring literal adjectice",w,d,adjectives,dmel);
			continue;
		}
		++aused;
		const StrTy sfx=suffix(a);
		const bool is_num=(sfx.length()==0);
		D vall=q_x_sz(a); // this calculates the a x v syntax result 
		//if (is_num) { 
		const D conv=sfx_map()[sfx];
		if (conv>0) val_added+=conv*vall; 
		else if (conv==0)
		{
			D mg_per_qty=qty_to_units(w,d,adjectives,dmel); //  
			val_added+=mg_per_qty*vall;
		}
	} // for 
	if (aused==0)
	{
		D mg_per_qty=qty_to_units(w,d,adjectives,dmel); //  
		val_added+=mg_per_qty;
	}
	val+=val_added;
*/
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
sfxv["oz"]=28350;
sfxv["mg"]=1;
sfxv["kg"]=1e6; // for safety let user change units lol 
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
bool m_ignore_adj, m_ignore_share;

}; // drugs


class kelp_handler : public drugs
{
typedef drugs Super;
public:
kelp_handler():Super(StrTy("kelp"), 99999) {
add_potency("LESI",1.0);
add_potency("nwkelp",.3);
add_potency("NOW",.325);
}

// this is only needed for the harder part
//virtual bool custom_qty() const { return !false; }
virtual D qty_to_units
	(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{
if (w=="LESI") { return 1.0;}
if (w=="nwkelp") { return .3;}
if (w=="kelp") { if (adj.size()!=0) if (adj[adj.size()-1]=="NOW") return .325;}
//m_handlers["Iodine"] = new Noun("Iodine",.300);
//m_handlers["Iodine"]->add_potency("LESI",1.0);
//m_handlers["Iodine"]->add_potency("nwkelp",.3);
//m_handlers["Iodine"]->add_potency("kelp",.325);
//if (dmel!=0) (*dmel).event("kelpfails",w,d,adj);
exp_log(__FILE__,__LINE__,"kelpfails"," qty_to_units",w,d,adj,dmel);
return Super::qty_to_units(w,d,adj,dmel);

}

}; // kelp_handler


class b12_handler : public drugs
{
typedef drugs Super;
public:
b12_handler():Super(StrTy("B-12"), 1.5) {
add_potency("NOW",1.5);
add_potency("GNC",1.5);
}
}; // b12_handler
class b2_handler : public drugs
{
typedef drugs Super;
public:
b2_handler():Super(StrTy("B-2"), 100) {
add_potency("NOW",100);
}
}; // b2_handler
class b3_handler : public drugs
{
typedef drugs Super;
public:
b3_handler():Super(StrTy("B-3"), 100) {
//add_potency("NOW",100);
ignore("bs"); 
}
}; // b3_handler
class d3_handler : public drugs
{
typedef drugs Super;
public:
d3_handler():Super(StrTy("D-3"), 1) {
//add_potency("NOW",100);
ignore("kroger"); 
ignore("iu"); 
}
}; // b3_handler


class k1_handler : public drugs
{
typedef drugs Super;
public:
k1_handler():Super(StrTy("K1"), 100000) {
//add_potency("NOW",100);
ignore("bs"); 
}
}; // b3_handler



class b1_handler : public drugs
{
typedef drugs Super;
public:
b1_handler():Super(StrTy("B-1"), 300) {
ignore("nutricost");
ignore("NOW");
}
}; // b2_handler



class mn_handler : public drugs
{
typedef drugs Super;
public:
mn_handler():Super(StrTy("Mn"), 8) {
add_potency("small",1);
}
}; // mn_handler



class b6_handler : public drugs
{
typedef drugs Super;
public:
b6_handler():Super(StrTy("B-6"), 99999) {
add_potency("NOW",100);
add_potency("GNC",100);
add_potency("sundown",50);
}
// this is only needed for the harder part
//virtual bool custom_qty() const { return !false; }
// this is really stupic 
virtual bool safe_to_ignore(const StrTy & a) { return (a=="sundown")||(a=="GNC")||(a=="NOW");}
// this code should never be called 
virtual D qty_to_units
	(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{
//if (w=="LESI") { return 1.0;}
//if (w=="nwkelp") { return .3;}
//if (w=="kelp") { if (adj.size()!=0) if (adj[adj.size()-1]=="NOW") return .325;}
exp_log(__FILE__,__LINE__,"b6obsoletecode"," qty_to_units",w,d,adj,dmel);
// must have adjective
const IdxTy as=adj.size();
if (as==0)
{
//if (dmel!=0) (*dmel).event("b6fails",w,d,adj);
exp_log(__FILE__,__LINE__,"b6fails"," qty_to_units",w,d,adj,dmel);
return 0; 
}
const StrTy at=adj[as-1];
if (at=="NOW") return 100;
if (at=="GNC") return 100;
if (at=="sundown") return 50;
//m_handlers["Iodine"] = new Noun("Iodine",.300);
//m_handlers["Iodine"]->add_potency("LESI",1.0);
//m_handlers["Iodine"]->add_potency("nwkelp",.3);
//m_handlers["Iodine"]->add_potency("kelp",.325);
exp_log(__FILE__,__LINE__,"b6fails"," qty_to_units",w,d,adj,dmel);
//if (dmel!=0) (*dmel).event("b6fails",w,d,adj);
return Super::qty_to_units(w,d,adj,dmel);

}

}; // b6_handler


class lecithin_handler : public drugs
{
typedef drugs Super;
public:
lecithin_handler():Super(StrTy("lecithin"), 1) {}
// this is only needed for the harder part
//virtual bool custom_qty() const { return !false; }

virtual bool safe_to_ignore(const StrTy & a) {
// needed now to make fake quantity appear right 
//return false; 
 return (a=="bulk");
}
virtual void qty_from_adj2(D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
const IdxTy sz=adjectives.size();
//for (IdxTy i=0; i<sz; ++i) { if (adjectives[i]=="bulk") { val=1200; return; }}
Super::qty_from_adj2(val,cname,w,d,adjectives,lit,dmel);
for (IdxTy i=0; i<sz; ++i) { if (adjectives[i]=="bulk") { val=val*1200; return; }}
}

// not used 
virtual D qty_to_units
	(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{
//if (w=="LESI") { return 1.0;}
//if (w=="nwkelp") { return .3;}
//if (w=="kelp") { if (adj.size()!=0) if (adj[adj.size()-1]=="NOW") return .325;}
// must have adjective
const IdxTy as=adj.size();
if (as==0) { return 1; }

const StrTy at=adj[as-1];
//if (at=="bulk") return 1;
// this is just to make it look better when the 1200mg capsules ere
// added for dental andto use them up 
if (at=="bulk") return 1200;

//	Ss ss; ss<<" adj "<<a<<" used "<<adj_used[i];
exp_log(__FILE__,__LINE__,"lecithinfails"," qty_to_units",w,d,adj,dmel);
//if (dmel!=0) (*dmel).event("lecithinfails",w,d,adj);


return Super::qty_to_units(w,d,adj,dmel);

}
}; // lecithin_handler



class Cu_handler : public drugs
{
typedef drugs Super;
public:
Cu_handler():Super(StrTy("Cu"), 1) 
{
ignore("carlson");
ignore("glycinate");
}

}; //  Cu_handler
/// magnesium 

class mg_handler : public drugs
{
typedef drugs Super;
public:
mg_handler():Super(StrTy("Mg"), 200) {
// qty should still be right 
//ignore("bulk");
}

virtual bool safe_to_ignore(const StrTy & a) {
// needed now to make fake quantity appear right 
//return false; 
 return (a=="bulk");
}
virtual void qty_from_adj2(D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
const IdxTy sz=adjectives.size();
//for (IdxTy i=0; i<sz; ++i) { if (adjectives[i]=="bulk") { val=1200; return; }}
Super::qty_from_adj2(val,cname,w,d,adjectives,lit,dmel);
}

// not used 
virtual D qty_to_units
	(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{
MM_ERR( "should not be called")
// must have adjective
const IdxTy as=adj.size();
if (as==0) { return 1; }

const StrTy at=adj[as-1];
//if (at=="bulk") return 1;
// this is just to make it look better when the 1200mg capsules ere
// added for dental andto use them up 
if (at=="bulk") return 1200;
exp_log(__FILE__,__LINE__,"lecithinfails"," qty_to_units",w,d,adj,dmel);
return Super::qty_to_units(w,d,adj,dmel);
}
}; // mg_handler












class salmon_handler : public drugs
{
typedef drugs Super;
public:
salmon_handler():Super(StrTy("salmon"), 1) 
{
ignore("keta");
ignore("coho");
ignore("pink");
ignore("wild");
ignore("atlantic");
ignore("walmart");
ignore("kroger");
ignore("sockeye");
ignore("pacific");
ignore("farmed");
ignore("chile");
}

}; //  salmon_handler
class tuna_handler : public drugs
{
typedef drugs Super;
public:
tuna_handler():Super(StrTy("tuna"), 1) 
{
ignore("pink");
ignore("wild");
ignore("atlantic");
ignore("walmart");
ignore("kroger");
ignore("pacific");
ignore("farmed");
}

}; //  tuna_handler



class vitamina_handler : public drugs
{
typedef drugs Super;
public:
vitamina_handler():Super(StrTy("vitamina"), 1) 
{
ignore("iu");
ignore("bronson");
}

}; //  salmon_handler



class garlic_handler : public drugs
{
typedef drugs Super;
public:
garlic_handler():Super(StrTy("garlic"), 1) 
{
ignore("minced");
ignore("pressed");
ignore("chopped");
}

}; //  garlic_handler

class leucine_handler : public drugs
{
typedef drugs Super;
public:
leucine_handler():Super(StrTy("leucine"), 625) 
{
ignore("nutricost");
ignore("bulk");
}

}; //  leucine_handler

class lysinehcl_handler : public drugs
{
typedef drugs Super;
public:
lysinehcl_handler():Super(StrTy("lysinehcl"), 650) 
{
ignore("nutricost");
ignore("bulk");
}

}; //  lysinehcl_handler


class taurine_handler : public drugs
{
typedef drugs Super;
public:
taurine_handler():Super(StrTy("taurine"), 650) 
{
ignore("nutricost");
ignore("bulk");
}

}; //  taurine_handler




class biotin_handler : public drugs
{
typedef drugs Super;
public:
biotin_handler():Super(StrTy("biotin"), 5) 
{
ignore("nutricost");
}

}; //  biotin_handler


class pantothenate_handler : public drugs
{
typedef drugs Super;
public:
pantothenate_handler():Super(StrTy("pantothenate"), 500) 
{
ignore("nutricost");
}

}; //  pantothenate_handler




class spinach_handler : public drugs
{
typedef drugs Super;
public:
spinach_handler():Super(StrTy("spinach"), 1) 
{
ignore("frozen");
ignore("raw");
ignore("birdseye");
}

}; //  bor_l_immune_handler


class bor_l_immune_handler : public drugs
{
typedef drugs Super;
public:
bor_l_immune_handler():Super(StrTy("bor-l-immune"), 1) 
{
ignore("drops");
}

}; //  bor_l_immune_handler



class lipoicacid_handler : public drugs
{
typedef drugs Super;
public:
lipoicacid_handler():Super(StrTy("lipoicacid"), 1) 
{
ignore_share(!true); // TODO rerun check this
//ignore_adj(!true);
ignore("best");
ignore("naturals");
ignore("natural");
}

}; //  lipoicacid_handler




class thiosulfate_handler : public drugs
{
typedef drugs Super;
public:
thiosulfate_handler():Super(StrTy("thiosulfate"), 1) 
{
ignore_share(true);
}
}; //  thiosulfate_handler
class metasilicate_handler : public drugs
{
typedef drugs Super;
public:
metasilicate_handler():Super(StrTy("metasilicate"), 1) 
{
ignore_share(true);
}
}; //  metasilicate_handler



class kcl_handler : public drugs
{
typedef drugs Super;
public:
kcl_handler():Super(StrTy("kcl"), 1) 
{
ignore("extra");
ignore("small");
}

}; //  salmon_handler




class shrimp_handler : public drugs
{
typedef drugs Super;
public:
//shrimp_handler():Super(StrTy("B-6"), 99999) {}
shrimp_handler():Super(StrTy("shrimp"), 99999) {} // TODO WTF 2018-12-15???
// this is only needed for the harder part
virtual bool custom_qty() const { return !false; }
// this is NOT called if we have a custom _from_units
virtual D qty_to_units
	(const StrTy & w, const StrTy & d, const Words & adj, Dmel * dmel=0) const 
{
	Ss ss; 
//ss<<" adj "<<a<<" used "<<adj_used[i];
ss<<" not implemente in shrimp handler "<<MMPR2(w,d);
	exp_log(__FILE__,__LINE__,"shrimp_not_impl",
		ss.str().c_str(),w,d,adj,dmel);
//MM_ERR(" not implemente in shrimp handler "<<MMPR2(w,d))
//if (dmel!=0) (*dmel).event("shrimp_not_impl",w,d,adj);

return 0; 

}

virtual void qty_from_adj2(D& val, const StrTy & cname,const StrTy & w
		,const StrTy & d, const Words & adjectives, const LiteralClass & lit, Dmel * dmel=0)
{
// default impl
StrTy fcname=cname; //  cons t
 qty_from_adj(fcname,val,w,d,adjectives,dmel);
}
// really this is obsolete but new sig calls it 
virtual void qty_from_adj(StrTy & cname, D& val,const StrTy & w
		,const StrTy & d, const Words & adj, Dmel * dmel=0)
{

const IdxTy asz=adj.size();
if (asz==0) {val+= gm_per_count(50,41);  return; } 
 
if (asz>1)
{
	MM_ERR(" more than one shrimp modifiied "<<MMPR4(w,d,adj[0],adj[1]))
	if (dmel!=0) (*dmel).event("shrimpfails",w,d,adj);
}
// the dmel should have a pushable context and know the date or line being processed
const D sc=parse_shrimp(adj[0]);
if (sc<0)
{
	MM_ERR(" bad  shrimp modifiied "<<MMPR3(w,d,adj[0]))
	if (dmel!=0) (*dmel).event("shrimpmodfails",w,d,adj);
	return; 
}
val+=sc;
} // qty_from_adj
D parse_shrimp(const StrTy & w)
{
const IdxTy sz=w.length();
if (sz<6) return -1; 
const char * s=w.c_str();
//const IdxTy cnt=atoi(s);
// this MAY be finding a "radix character" for 0x wtf lol
D cnt=atof(s);
if (*s=='0') if (s[1]=='x') cnt=0;
const IdxTy maxn=(s[sz-1]-'0')+(s[sz-2]-'0')*10;
const IdxTy minn=(s[sz-3]-'0')+(s[sz-4]-'0')*10;
return cnt*gm_per_count(maxn,minn);
}
D gm_per_count(const D & maxn, const D & minn)
{
const D av=(maxn+minn)*.5; // argue about < 1/x > lol 
return 16.0/av*g_to_oz(); 
}
//static D g_to_oz() { return 2.20462*16; } 
static D g_to_oz() { return 28.3495; } 
}; // shrimp_handler






class mjm_snacks
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



typedef mjm_snacks Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

//typedef mjm_logic_base Logic;
typedef snack_params Logic;
typedef mjm_logic_base VariableStore;
// the ~0 values were failing probably due to sign problems wtf 
enum NAME_CLASS { DOG_NAME=0, ALIAS=1, MACRO=2, TIME=3, DATE=4, NOTE=5, HYPOTHESIS=6,DMEL=7, COMMENT=DMEL+1, INVALID=0 };
// list of things for a given date
typedef std::map<StrTy,StrTy > DateInfo;


typedef std::map<StrTy,StrTy > Canonical;
typedef std::vector<StrTy> CanonFails;


typedef std::vector<StrTy> DumpOrder;
//typedef std::map<StrTy,IdxTy > FreqInfo;
//typedef std::map<StrTy,IdxTy > Literals;

typedef literal_map LiteralClass;
//typedef LiteralClass::LITERAL_TYPE LITERAL_TYPE;
//typedef literal_map LiteralClass;
typedef LiteralClass::LITERAL_TYPE LIT_TYPE;
enum  LITERAL_TYPE { ADJECTIVE=LIT_TYPE::ADJECTIVE, QUANTIFIER=LIT_TYPE::QUANTIFIER
		, NOUN=LIT_TYPE::NOUN, LMACRO=LIT_TYPE::LMACRO
		,IGNORE=LIT_TYPE::IGNORE, NONE=LIT_TYPE::NONE };


// collection of things for each date 
typedef std::map<StrTy, DateInfo > MappedInfo;
typedef std::map<StrTy, D > ShareFactor;
typedef std::map<StrTy, ShareFactor > ShareFactorMap;

typedef std::vector<StrTy> Words;
typedef drugs Noun;
typedef std::map<StrTy, Noun*> NounHandler;
typedef data_model_error_log Dmel;

/// need to have dog_day def available. 
// MACRO is indexed by name and then date 
typedef std::map<StrTy , std::map< StrTy, dog_day> > MacroMap;
typedef std::map< StrTy,  std::map< StrTy, dog_day> > DogDays;
// generic map indxed by dog then day and containing a string /
typedef std::map< StrTy,  std::map< StrTy, std::vector<StrTy> > > DogDayMap;

typedef mjm_calendar CalTy;

StrTy name_class_name(const int i)
{
static std::map<int, StrTy > x;
static bool init=false;
if (!init)
{
x[DOG_NAME]=StrTy("DOG_NAME");
x[ALIAS]=StrTy("ALIAS");
x[MACRO]=StrTy("MACRO");
x[TIME]=StrTy("TIME");
x[DATE]=StrTy("DATE");
x[NOTE]=StrTy("NOTE");
x[COMMENT]=StrTy("COMMENT");
x[HYPOTHESIS]=StrTy("HYPOTHESIS");
x[DMEL]=StrTy("DMEL");
x[INVALID]=StrTy("INVALID");
init=true;
}
if (x.find(i)==x.end()) { 
if (m_dmel!=0) { Ss ss; ss<<i; (*m_dmel).event("bad class ",ss.str()); } 
MM_INC_MSG(m_cm,"badclass" )
MM_ERR(" bad class "<<i); return StrTy("NOTFOUND"); } 
return x[i];
}




public :
mjm_snacks():m_dmel(new Dmel()) {Init();}
mjm_snacks(int argc,char **_args) : m_dmel(new Dmel())
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
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]<<" "<<i) } 

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
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
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
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
//if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
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
if (cmd=="dump-dmel") { dump_dmel(std::cout); continue; } 
if (cmd=="dump-dmel-cerr") { dump_dmel(std::cerr); continue; } 
if (cmd=="dump-dog-days") { dump_dog_days(std::cout,p1); continue; } 
if (cmd=="dump-dogs") { dump_dogs(std::cout,p1); continue; } 
if (cmd=="dump-dog-times") { dump_dog_times(std::cout,p1); continue; } 
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


// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. 
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

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general

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
if (name=="HYPOTHESIS") return HYPOTHESIS;
if (name=="DMEL") return DMEL;
	// TODO FIXME this is not right, a NOTE si structured 
//if (name=="COMMENT") return NOTE;
if (name=="COMMENT") return COMMENT;
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
template <class Tv> 
void new_toke(Tv & v,  char * c, IdxTy& ptr, const IdxTy i)
{ c[i]=0; v.push_back(StrTy(c+ptr)); ptr=i+1; }  
// each macro gets parsed into a "dog_day" item 
// this returns an evaluated hacked up set of defininte quantities
dog_day macro_lut(const StrTy & w, const StrTy & d, const Words & adj )
{
// find base name,
typedef std::vector<StrTy> V;
V words;
V parameters;
//const char * c=w.c_str();
// go into param parsing 
IdxTy state=0;
const IdxTy sz=w.length();
char c[sz+1];
memcpy(c,w.c_str(),sz);
c[sz]=0;
 IdxTy ptr=0;
for (IdxTy i=0; i<sz; ++i)
{
const  char cp =c[i];
if (cp=='-') { c[i]=0; words.push_back(StrTy(c+ptr)); ptr=i+1; }
else if (cp=='[' ) { ++state; new_toke(words,c,ptr,i); } // if 
else if (cp==']' ) { --state; new_toke(parameters,c,ptr,i); } // if 
else if (cp==',' ) {  c[i]=0; new_toke(parameters,c,ptr,i); } // if 

} // for i

if ( state!=0)
{
MM_MSG(" macro param error  "<<words[0]<<" from "<<w<<" "<<state<<" "<<d)
exp_log(__FILE__,__LINE__,"nomacrodef","bad  macro param def ",w,d);
MM_INC_MSG(m_cm, "nomacroever")

}

if (ptr<sz) words.push_back(StrTy(c+ptr));
if (words.size()==0) words.push_back(StrTy(""));
auto ii=m_macro.find(words[0]);
if (ii==m_macro.end())
{
MM_MSG(" no macro ever for "<<words[0]<<" from "<<w<<" "<<d)
exp_log(__FILE__,__LINE__,"nomacrodef","missing macro def ",w,d);
MM_INC_MSG(m_cm, "nomacroever")
dog_day dd;
return dd;
}
auto jj=(*ii).second.find(d);
if (jj==(*ii).second.end())
{

MM_MSG(" no macro  for "<<words[0]<<" "<<d<<" from "<<w<<" "<<d)
exp_log(__FILE__,__LINE__,"nocurrentmacrodef","missing applicable macro def ",w,d);
MM_INC_MSG(m_cm, "nomacrodate")
dog_day dd;
return dd;

}
// so this returns a COPY of 
dog_day dd=m_macro[words[0]][d];
// if there is more than one adjective, use the others as variables

const IdxTy wsz=words.size();
// now subtract off all the entries with base names in minus list 
for (IdxTy i=1; i<wsz; ++i) { dd.set_dose(words[i],0); }
//MM_ERR(" found a macro "<<MMPR4(w,d,dd.open(),sz))
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
// w is the macro names
D share_factor(const StrTy & n,const StrTy & d,const StrTy & w)
{
const IdxTy sz=m_adjectives.size();
if (sz!=0)
{
	const char * s=m_adjectives[0].c_str();
	// DOH TODO zero had been disallowed and assumed nonumberic 
	const D f=atof(s);
//	const D f=atof(m_adjectives[0].c_str());
	//if ((f<=0)||(sz>1))
	//if ((f<=0)||(f>1)) // now variables are allowed  and f>1 ok ???
	if ((f<=0)) // now variables are allowed  and f>1 ok ???
	{
		const IdxTy len=strlen(s);
		if (len>0)
		{
			char c=*s;
			IdxTy i=1;
			while (c==' ') {c=s[i]; ++i; } 
			if (c=='0') return f; 
			if (c=='.') return f; 

		}
		if (m_dmel!=0) { (*m_dmel).event("badmacroshare" ,w, d,m_adjectives); } 
	}
	else return f; 
}
auto ii=m_share_factors.find(w);
if (ii!=m_share_factors.end()) return (*ii).second[n];
if (d>StrTy("2018-09")) 
if (m_dmel!=0) { (*m_dmel).event("obsoletesharefactor" ,w, d,m_adjectives); } 

return m_share_factor[n];


}
StrTy tod_from_string(const StrTy & w)
{
if (w.length()!=6) return StrTy("BAD")+w;
return w.substr(0,4);
}

// returning an instance hopefully works right.. 
dog_day parse_dog_day( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const Words & w=li.words();
	return parse_dog_day(d,n,start,w);
}

dog_day parse_dog_day( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li, const Words & vars)
{
	Words  w=li.words();
	return parse_dog_day(d,n,start,w,vars);
}

IdxTy is_var( const StrTy & w, const IdxTy len )
{
	IdxTy st=0;
	const char * p=w.c_str();
//	const char c=*p; // may be null but it shold be there
	if (p[0]=='$') st=1;
	if ( len>2) if (p[0]=='\\') if ( p[1]=='$') st=2;
	return st;
}
// make a COPY of the word 
dog_day parse_dog_day( const StrTy & d, const StrTy & n, const IdxTy start,   Words  w, const  Words & vars)
{


	const IdxTy sz=w.size();
	const IdxTy szv=vars.size();
	std::vector<IdxTy>  used(szv);
//MM_ERR(" looking for subs "<<MMPR4(d,n,w[0],vars.size()))
	for(IdxTy i=0; i<sz; ++i)
	{
		const IdxTy len=w[i].length();
		if ( len<2) continue; 
		const char * p=w[i].c_str();
		const char c=*p; // may be null but it shold be there
		//MM_ERR(" check "<<MMPR3(c,p,i))
		IdxTy st=is_var(w[i],len) ; // 0;
//		if (c=='$') st=1;
//		if ( len>2) if (c=='\\') if ( *(p+1)=='$') st=2;
		if ( st!=0) 
		{
			const IdxTy v=::atoi(p+st);
			if ( v<szv) {  w[i]=vars[v]; ++used[v]; MM_ERR(MMPR4(v,w[i],vars[v],(p+1))) }
			else
			{
				if (m_dmel!=0) (*m_dmel).event("badvariables",w[0],d);
				MM_ERR(" bad variables "<<MMPR4(w[0],v,sz,szv)<<MMPR(d))
			}
		}
	}
	MM_SZ_LOOP(i,used,szu)
	{
		if ( i==0) continue; // qty 
		if ( (used[i])==0){
				if (m_dmel!=0) (*m_dmel).event("missingvariables",w[0],d);
				MM_ERR(" missing  variables "<<MMPR4(w[0],i,sz,szv)<<MMPR(d))
			}
	}
	return parse_dog_day(d,n,start,w);
}
// this may be called with vars, 
dog_day parse_dog_day( const StrTy & d, const StrTy & n, const IdxTy start,  const Words & _w)
{
MM_ONCE(" backed up parse_dog_day check results ",)

	const bool parse_time_of_day=!false;
	dog_day dd(n,d);
	const IdxTy sz=_w.size(); // li.size();
	// FIXME tthis member screws up recursion 
	m_adjectives.clear();
	StrTy tod="0000";
	tod="0430";
	//StrTy tods="tod";
	if (parse_time_of_day) dd.set_tod(tod);
	for (IdxTy i=start; i<sz; ++i)
	{
		const StrTy w=_w[i]; // li.word(i);
		if (w=="") continue; 
		IdxTy st=is_var( w, w.length());
		if ( st!=0) { dd.reparse(_w) ; return dd; } 
//		const char c=w[0];
//		if ( c=='$') { dd.reparse(_w) ; return dd; } 
//		if ( c=='\\') {if (w.length()>2) if ( w[1]=='$') { dd.reparse(_w) ; return dd;} } 
		const bool dump_qty_parse=false; // (w=="pantothenate");
		if (dump_qty_parse) {MM_MSG(MMPR4(d,n,i,w))}
//enum LITERAL_TYPE { ADJECTIVE=1, QUANTIFIER=2, NOUN=3, LMACRO=5,IGNORE=4, NONE=~0 };
		//LITERAL_TYPE lt=NONE;
		IdxTy lt=NONE; //  THESE ENUNS  
		if (m_literals.map().find(w)!=m_literals.map().end())
		{
			lt=m_literals.map()[w];
			if (lt==IGNORE)
			{	// do NOT ignore adjectives just nouns and macro 
				m_adjectives.clear();
				continue;
			}
		}
		// could be adjective, noun, or macro. Macro is self contained, no leaders
		// but maybe trailing adjectives macro-missing  
		//  FIXME recursive adjectives needs to be passed not a memebnr
		if (macro_ref(w,d,lt))
		{
			if (parse_name(w)==TIME)
			{
				if (!parse_time_of_day) continue;	
				// in theory this just requires getting time words and keeping track except that
				// PPHOME etc rely on the PMDINNER tod info. 
				tod=tod_from_string(w);
				dd.set_tod(tod);
				continue; 
			}
	//		MM_ERR( "looking up "<<MMPR2(w,d))
			// NB this is not a const ref 
			dog_day macro= macro_lut(w,d,m_adjectives);
			
			if ( macro.open()) 
			{
				MM_ERR( " making var subgs   ") 
				 // SO NOW WHAT? 
//			 macro= macro_lut(w,d,m_adjectives);
//dog_day parse_dog_day( const StrTy & d, const StrTy & n, const IdxTy start,  const Words & w, const Words & vars)
			 Words x=m_adjectives;
			 macro= parse_dog_day(d,n,start,macro.def(),m_adjectives);
			m_adjectives=x;
			}

			if (w=="PPHOME") { if (macro.m_tod<tod) macro.set_tod(StrTy("1300")); }
			if (w=="0930AMSNACK") { if (macro.m_tod<StrTy("0730AM") ) macro.set_tod(StrTy("0930")); }

			// rhia only uses the FIRST adjective, others are macro sub
			// this now recognizes leading adjectives or quantifiers
			const D macfac=share_factor(n,d,w);
			//MM_MSG(MMPR4(w,macfac,n,d))
			dd.add(macro,macfac,m_handlers);
			m_adjectives.clear();
		} else
		// if adjective just push, otherwise use noun's rules for
		// popping the adjectives 
		if (adjective_ref(w,d,lt))
		{ 	// these may be quantifiers or modifiers
			// are times included here?
			m_adjectives.push_back(w);
		} else
		// this could be "everything ELSE except should drop times in adjectives
		// noun needs to use up all the adjectives 
		{
			StrTy cname=w;
			D val=0;
//			qty_from_adj_old(cname, val,w,d);
			qty_from_adj(cname, val,w,d);
MM_ERR(MMPR4(n,d,cname,val)<<MMPR(w))
			//if (dump_qty_parse) {MM_MSG(MMPR4(w,cname,val,m_adjectives.size()))}
			//if ((cname=="Mg")&&(m_adjectives.size()>0)) {MM_ERR(MMPR4(w,cname,val,m_adjectives.size()))}
			dd.inc(cname,val);
			if (parse_time_of_day) 
			{ // this does not pickup the macro cogtrbitionbs. 
//			std::cout<<d<<" "<<MMPR4(tods,n,tod,cname)<<MMPR(val)<<CRLF;
			}	
			m_adjectives.clear();
		}

	}
if (false)	MM_MSG(dd.string())
	return dd;
}
// this defines a possibly partial set of quantities for later evaluation
void parse_macro( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	//a macro is just like the dog entry except it is a definition for that date
	// may return an incoimpletely evaluated thing now.  
	dog_day dd=parse_dog_day( d, n, start, li);
	//MM_MSG(" no macro spurious "<<MMPR2(n,d))
	m_macro[n][d]=dd;
	

}
///////////////////////////////////
void parse_alias( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	// need to check the data before doing this could be a mess lol 
	//an AKA alias defines a new noun for now   
	//MM_MSG(" no macro spurious "<<MMPR2(n,d))
// Kibble typica
//m_canon[StrTy("KibbleMixDiamondJourney")]=StrTy("kibble");
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
//	for (IdxTy i=start; i<sz; ++i)
	//if (m_canon.find(name)!=m_canon.end())
	if (m_literals.is_noun(name))
	{
		if (m_dmel!=0) 
		{ (*m_dmel).event("akanounchanged" ,name, d,w); } 

	}
// will wreck some things 
	//m_canon[n]=w[start];
// the AKA simply is defininf a canonical name 
//	m_canon[name]=name; // w[start];
 //add_canonical_noun(name);
 m_literals.add_noun(name);
//MM_MSG(" setting aka  "<<MMPR2(n,w[start]))
//MM_ERR(" setting aka  "<<MMPR2(n,w[start]))
	

}



////////////////////////////////////
void parse_comment( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
	StrTy entry=(li.rest_of_line(start-1));
	static CalTy  cal;
	// TODO FIXME this does not work ... 
	StrTy dayno=""; // cal.xlate(w[0]);
	entry=dayno+StrTy(" ")+entry;
	// insert day no 
	m_notes[name][d].push_back(entry );
// do nothing with this for now 
}


////////////////////////////////////
void parse_note( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
	StrTy entry=(li.rest_of_line(start-1));
	static CalTy  cal;
	// TODO FIXME this does not work ... 
	StrTy dayno=""; // cal.xlate(w[0]);
	entry=dayno+StrTy(" ")+entry;
	// insert day no 
	m_notes[name][d].push_back(entry );
// do nothing with this for now 
}


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
///////////////////////////////////////////////

void dump_dogs(OsTy & os, const StrTy & d=StrTy())
{// this is not ORDERED but ok for now 
//std::map<StrTy, IdxTy> vmap;
std::map<StrTy , IdxTy >  dogmap;
for (auto ii=m_dog_days.begin(); ii!=m_dog_days.end(); ++ii)
{
	const StrTy dii=(*ii).first;
	const auto & dog_for_day=(*ii).second;
	for (auto jj=dog_for_day.begin(); jj!=dog_for_day.end(); ++jj)
	{
		const StrTy dog=(*jj).first;
		++dogmap[dog];
//		os<<res<<MMPR(vmap[cw])<<CRLF;
	} //  jj
} // ii 
MM_LOOP(ii,dogmap)
{
os<<" dump_dogs "<<MMPR2((*ii).first, (*ii).second) <<CRLF;

}
} // dump_dogs

////////////////////////////////////////////////

void dump_dog_days(OsTy & os, const StrTy & d=StrTy())
{// this is not ORDERED but ok for now 
//std::map<StrTy, IdxTy> vmap;
std::map<StrTy , std::map<StrTy, IdxTy> >  vmapmap;
for (auto ii=m_dog_days.begin(); ii!=m_dog_days.end(); ++ii)
{
	const StrTy dii=(*ii).first;
	const auto & dog_for_day=(*ii).second;
	for (auto jj=dog_for_day.begin(); jj!=dog_for_day.end(); ++jj)
	{
		const StrTy dog=(*jj).first;
		if (d!="") if (dog!=d) continue;
		auto & vmap = vmapmap[dog];
		const StrTy res=(*jj).second.string(m_order);
		const StrTy cw=(*jj).second.dose_string(m_order);
		if (vmap.find(cw)==vmap.end()) vmap[cw]=vmap.size();
//		MM_ERR(MMPR3(dog,vmap.size(),cw))
		os<<res<<MMPR(vmap[cw])<<CRLF;
	} //  jj
} // ii 
} // dump_dog_days

template <class Timap> void index_order(Timap & omap)
{
for (auto ii=m_order.begin(); ii!=m_order.end(); ++ii) ++omap[(*ii)];

}
void dump_dog_times(OsTy & os, const StrTy & d=StrTy())
{// this is not ORDERED but ok for now 
//std::map<StrTy, IdxTy> vmap;
//std::map<StrTy , std::map<StrTy, IdxTy> >  vmapmap;
typedef std::map<StrTy, int> Omap;
Omap omap;
index_order(omap);
const bool print_conflicts=m_flp.print_conflicts();
const IdxTy flags=(print_conflicts?1:0);
for (auto ii=m_dog_days.begin(); ii!=m_dog_days.end(); ++ii)
{
	const StrTy dii=(*ii).first;
	const auto & dog_for_day=(*ii).second;
	for (auto jj=dog_for_day.begin(); jj!=dog_for_day.end(); ++jj)
	{
		const StrTy dog=(*jj).first;
		if (d!="") if (dog!=d) continue;
//		auto & vmap = vmapmap[dog];
		const StrTy res=(*jj).second.time_string(omap,flags);
//		const StrTy cw=(*jj).second.dose_string(m_order);
//		if (vmap.find(cw)==vmap.end()) vmap[cw]=vmap.size();
//		MM_ERR(MMPR3(dog,vmap.size(),cw))
		os<<res; // <<CRLF;
	} //  jj
} // ii 
} // dump_dog_times



void stream_counts(OsTy & os)
{
	m_literals.stream_counts(os); 
/* 
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
*/
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
	li.set_split(4,' '); // 2019-05-14
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
				if (dd.open())
				{

				if (m_dmel!=0) (*m_dmel).event("missing valuses ",name,date);
				MM_ERR(" bad variables "<<MMPR3(date,name,dd.string())) 

				}
				if (print_dog_days)	MM_MSG(dd.string(m_order))
				// there shold only be one entry per dog per day but
				// actually allowing more with sum  here should work 
				if (accumulate_dog_days){ m_dog_days[date][name]=dd;  }	
 				break; }
			case MACRO : { parse_macro(date,name,3,li); break; }
			case NOTE : { parse_note(date,name,3,li); break; }
			case COMMENT : { parse_comment(date,name,3,li); break; }
			case HYPOTHESIS : {  break; }
			case DMEL : { parse_dmel(date,name,3,li); break; }
			case ALIAS: { parse_alias(date,name,3,li); break; } // this is defined as TEXT not like MACRO 
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
lit.add_adjective("bs");
// ??? TODO WTF 
lit.add_adjective("adjective");
lit.add_adjective("delta");
lit.add_adjective("extra");
// ignore?
//lit.add_adjective("and");
lit.add_adjective("small");
lit.add_adjective("minced");
lit.add_adjective("some");
lit.add_adjective("wild");
lit.add_adjective("keta");
lit.add_adjective("coho");
lit.add_adjective("walmart");
lit.add_adjective("atlantic");
lit.add_adjective("chile");
lit.add_adjective("farmed");
lit.add_adjective("unknown");
lit.add_adjective("bulk");
lit.add_adjective("best");
lit.add_adjective("natural");
lit.add_adjective("naturals");
lit.add_adjective("drop");
lit.add_adjective("drops");
lit.add_adjective("crystal"); // used with thiosulfate
lit.add_adjective("nutricost"); // 
lit.add_adjective("canned"); // 
lit.add_adjective("cooked"); // 


lit.add_adjective("raw"); // 
lit.add_adjective("frozen"); // 
lit.add_adjective("birdseye"); // 

lit.add_adjective("kroger"); // 
lit.add_adjective("sockeye"); // 
lit.add_adjective("pink"); // 
lit.add_adjective("carlson"); // 
lit.add_adjective("glycinate"); // 
lit.add_adjective("pacific"); // 
lit.add_adjective("bronson"); // 
lit.add_adjective("iu"); // 

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
lit.add_ignore("novaline");
lit.add_ignore("nothreonine");
lit.add_ignore("nomethionine");
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
// this is a noun rather than adjective as it refers to a general "thing"
m_canon[StrTy("cooked")]=StrTy("cooked");

m_canon[StrTy("Mgcitrate")]=StrTy("Mg");
m_canon[StrTy("Ivermectin")]=StrTy("Ivermectin");
m_canon[StrTy("ivermectin")]=StrTy("Ivermectin");
m_canon[StrTy("Mn")]=StrTy("Mn");
m_canon[StrTy("LESI")]=StrTy("Iodine");
m_canon[StrTy("nwkelp")]=StrTy("Iodine");
m_canon[StrTy("kelp")]=StrTy("Iodine");
m_canon[StrTy("Cu")]=StrTy("Cu");
m_canon[StrTy("HMCU")]=StrTy("HMCU");
m_canon[StrTy("HMCu")]=StrTy("HMCU");
m_canon[StrTy("K2")]=StrTy("K2");
m_canon[StrTy("K1bs")]=StrTy("K1");
m_canon[StrTy("noK2")]=StrTy("noK2");
m_canon[StrTy("11KC")]=StrTy("11KC");
m_canon[StrTy("11PC")]=StrTy("11PC");
m_canon[StrTy("10PC")]=StrTy("10PC");
m_canon[StrTy("01PC")]=StrTy("01PC");
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
m_canon[StrTy("D-3")]=StrTy("D-3");
m_canon[StrTy("garlic")]=StrTy("garlic");
m_canon[StrTy("b7ngnc")]=StrTy("b7ngnc");
m_canon[StrTy("b8ngnc")]=StrTy("b9ngnc");
m_canon[StrTy("b10ngnc")]=StrTy("b10ngnc");
m_canon[StrTy("b20ngnc")]=StrTy("b20ngnc");
m_canon[StrTy("b25ngnc")]=StrTy("b25ngnc");
m_canon[StrTy("bssngnc")]=StrTy("bssngnc");
m_canon[StrTy("b20ng")]=StrTy("b20ng");
m_canon[StrTy("b15ng")]=StrTy("b15ng");
m_canon[StrTy("nob20ngnc")]=StrTy("nob20ngnc");
m_canon[StrTy("b15ngnc")]=StrTy("b15ngnc");
m_canon[StrTy("nob15ngnc")]=StrTy("nob15ngnc");
m_canon[StrTy("p20ngnc")]=StrTy("p20ngnc");
m_canon[StrTy("carrot")]=StrTy("carrot");
m_canon[StrTy("carrots")]=StrTy("carrot");
m_canon[StrTy("canned")]=StrTy("canned");
m_canon[StrTy("biotin")]=StrTy("biotin");
m_canon[StrTy("bor-l-immune")]=StrTy("bor-l-immune");
m_canon[StrTy("lipoicacid")]=StrTy("lipoicacid");
m_canon[StrTy("pantothenate")]=StrTy("pantothenate");
m_canon[StrTy("metasilicate")]=StrTy("metasilicate");
m_canon[StrTy("thiosulfate")]=StrTy("thiosulfate");
m_canon[StrTy("arginine")]=StrTy("arginine");
m_canon[StrTy("lysine")]=StrTy("lysine");
m_canon[StrTy("leucine")]=StrTy("leucine");
m_canon[StrTy("isoleucine")]=StrTy("isoleucine");
m_canon[StrTy("valine")]=StrTy("valine");
m_canon[StrTy("methionine")]=StrTy("methionine");
m_canon[StrTy("lysinehcl")]=StrTy("lysinehcl");
m_canon[StrTy("d-serine")]=StrTy("d-serine");
m_canon[StrTy("SnAgOx")]=StrTy("SnAgOx");
m_canon[StrTy("SnAg")]=StrTy("SnAg");

//#if 0 
m_canon[StrTy("tuna")]=StrTy("tuna");
m_canon[StrTy("salmonbrothbs")]=StrTy("salmon");
m_canon[StrTy("ketasalmonbrothbs")]=StrTy("salmon");
m_canon[StrTy("salmonbroths")]=StrTy("salmon");
m_canon[StrTy("salmonbrothb")]=StrTy("salmon");
m_canon[StrTy("salmonbroth")]=StrTy("salmon");
m_canon[StrTy("salmon")]=StrTy("salmon");
m_canon[StrTy("salmons")]=StrTy("salmon");
m_canon[StrTy("wildsalmon")]=StrTy("salmon");
//#endif


m_canon[StrTy("ctbroth")]=StrTy("ctbroth");
m_canon[StrTy("ctbrothbs")]=StrTy("ctbroth");
m_canon[StrTy("ctbroths")]=StrTy("ctbroth");
// not exact match FIXME  
m_canon[StrTy("clbrothbs")]=StrTy("ctbroth");
m_canon[StrTy("cbbrothbs")]=StrTy("ctbroth");
m_canon[StrTy("cbbroth")]=StrTy("ctbroth");

m_canon[StrTy("ctbrothb")]=StrTy("ctbroth");
m_canon[StrTy("ctskin")]=StrTy("ctskin");
m_canon[StrTy("ctskinb")]=StrTy("ctskinb");
m_canon[StrTy("coffee")]=StrTy("coffee");
m_canon[StrTy("oliveoil")]=StrTy("oliveoil");
m_canon[StrTy("coconutoil")]=StrTy("coconutoil");


//m_canon[StrTy("spinach")]=StrTy("spinach");
m_canon[StrTy("chili")]=StrTy("chili");
m_canon[StrTy("Kibble")]=StrTy("kibble");
m_canon[StrTy("kibble")]=StrTy("kibble");
m_canon[StrTy("KibbleMixDiamondJourney")]=StrTy("kibble");
m_canon[StrTy("KibbleEliteSeriesBeef")]=StrTy("kibble");
m_canon[StrTy("shrimp")]=StrTy("shrimp");
m_canon[StrTy("vitamina")]=StrTy("vitamina");
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
m_canon[StrTy("diroban")]=StrTy("diroban");
//m_canon[StrTy("mgcitrate")]=StrTy("mgcitrate");
// 2018-02-18
m_canon[StrTy("clindamycin")]=StrTy("clindamycin");
m_canon[StrTy("prednisone")]=StrTy("prednisone");
m_canon[StrTy("cephalexin")]=StrTy("cephalexin");
m_canon[StrTy("doxycycline")]=StrTy("doxycycline");
m_canon[StrTy("Temaril-P")]=StrTy("Temaril-P");
m_canon[StrTy("metronidazole")]=StrTy("metronidazole");
m_canon[StrTy("ketoconazole")]=StrTy("ketoconazole");
m_canon[StrTy("serrano")]=StrTy("serrano");
m_canon[StrTy("habanero")]=StrTy("habanero");
m_canon[StrTy("jalapeno")]=StrTy("jalapeno");
m_canon[StrTy("teeth")]=StrTy("teeth");
m_canon[StrTy("brush")]=StrTy("brush");
m_canon[StrTy("tryptophan")]=StrTy("tryptophan");
m_canon[StrTy("tyrosine")]=StrTy("tyrosine");
m_canon[StrTy("phenylalanine")]=StrTy("phenylalanine");


m_canon[StrTy("threonine")]=StrTy("threonine");
m_canon[StrTy("histidinehcl")]=StrTy("histidinehcl");
m_canon[StrTy("cerenia")]=StrTy("cerenia");



add_canonicals_to_literals();
}
void add_canonical_noun(const StrTy & nm)
{
m_canon[nm]=nm;
m_literals.add_noun(nm);
}
void add_canonicals_to_literals()
{
//  enums  
for (auto ii=m_canon.begin(); ii!=m_canon.end(); ++ii)
	//m_literals.map()[(*ii).first]=NOUN;
	m_literals.add_noun((*ii).first);


}
void load_handlers()
{
// TODO FIXME danger will robinsom memory leak esp if assigning
// over existing one. Need to subclass and delete existing  
clear_handlers();
//m_handlers["pantothenate"] = new Noun("pantothenate",500);
m_handlers["pantothenate"] = new pantothenate_handler(); // Noun("biotin",5);
m_handlers["biotin"] = new biotin_handler(); // Noun("biotin",5);
m_handlers["bor-l-immune"] = new bor_l_immune_handler(); // Noun("biotin",5);
m_handlers["spinach"] = new spinach_handler(); // Noun("biotin",5);
//m_handlers["taurine"] = new Noun("taurine",1000);
m_handlers["taurine"] = new taurine_handler(); // Noun("biotin",5);
m_handlers["arginine"] = new Noun("arginine",500);
m_handlers["lysine"] = new Noun("lysine",500);
//m_handlers["lysinehcl"] = new Noun("lysinehcl",0);
m_handlers["leucine"] = new leucine_handler(); // Noun("biotin",5);
m_handlers["lysinehcl"] = new lysinehcl_handler(); // Noun("biotin",5);
//m_handlers["Mg"] = new Noun("Mg",200);
//  so the number on the front is per capsule not per made up serving 
//m_handlers["Mg"] = new Noun("Mg",200); // the  number on front is NOT per  SERVING 
m_handlers["Mg"] = new mg_handler(); // Noun("biotin",5);
m_handlers["Cu"] = new Cu_handler(); // Noun("biotin",5);
//m_handlers["B-1"] = new Noun("B-1",300);
m_handlers["B-1"] = new b1_handler(); // Noun("B-1",300);
m_handlers["B-2"] = new b2_handler(); // new Noun("B-6",100);
m_handlers["B-3"] = new b3_handler(); // new Noun("B-6",100);
m_handlers["D-3"] = new d3_handler(); // new Noun("B-6",100);
m_handlers["K1"] = new k1_handler(); // new Noun("B-6",100);
//m_handlers["B-2"] = new Noun("B-2",100);
//m_handlers["B-3"] = new Noun("B-3",100);
//m_handlers["B-3"] = new Noun("NSB-3",100);
//m_handlers["B-3"] = new Noun("SNB-3",100);
m_handlers["B-6"] = new b6_handler(); // new Noun("B-6",100);
m_handlers["B-12"] = new b12_handler(); // new Noun("B-6",100);
m_handlers["Mn"] = new mn_handler(); // new Noun("B-6",100);
m_handlers["lecithin"] = new lecithin_handler(); // new Noun("B-6",100);
m_handlers["10KC"] = new kcl_handler(); // new Noun("B-6",100);


m_handlers["tuna"] = new tuna_handler(); // new Noun("B-6",100);
m_handlers["salmon"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["salmonbrothbs"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["salmonbroths"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["salmonbrothb"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["salmonbroth"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["ketasalmonbrothbs"] = new salmon_handler(); // new Noun("B-6",100);
m_handlers["wildsalmon"] = new salmon_handler(); // new Noun("B-6",100);


//m_canon[StrTy("salmonbrothbs")]=StrTy("salmon");
//m_canon[StrTy("ketasalmonbrothbs")]=StrTy("salmon");
//m_canon[StrTy("salmonbroths")]=StrTy("salmon");
//m_canon[StrTy("salmonbrothb")]=StrTy("salmon");
//m_canon[StrTy("salmonbroth")]=StrTy("salmon");
//m_canon[StrTy("salmon")]=StrTy("salmon");
//m_canon[StrTy("wildsalmon")]=StrTy("salmon");



m_handlers["garlic"] = new garlic_handler(); // new Noun("B-6",100);
m_handlers["shrimp"] = new shrimp_handler(); // new Noun("B-6",100);
m_handlers["vitamina"] = new vitamina_handler(); // new Noun("B-6",100);
//  this is per  SERVING 2017-12-28ZZ
//m_handlers["d-serine"] = new Noun("d-serine",700);
m_handlers["d-serine"] = new Noun("d-serine",700.0/4.0);

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
// the main data indexed 
DogDays m_dog_days;
// dog,date, note map 
DogDayMap m_notes;
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

