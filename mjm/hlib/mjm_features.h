#ifndef MJM_FEATURES_H__
#define MJM_FEATURES_H__

#include "mjm_text_data.h"
#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
// needed for file io crap
#include <mjm_csv_ops.h>
// finally added rule hits
#include <mjm_sequence_hits.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

class features_typedefs
{
public:
typedef unsigned int IdxTy;
typedef  int IntTy;
typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::stringstream SsTy;
};

namespace mjm_features 
{
typedef features_typedefs Tr;

typedef Tr::IdxTy IdxTy;
typedef  Tr::IntTy IntTy;
typedef Tr::ChTy ChTy;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;
enum { BAD=~0U}; 

// right now these are only reverse-complement hits
class feature_index_entry
{

public:
feature_index_entry() : m_pos(0),m_mate(0),m_length(0),m_dist(0) {}

// should assert length if it is there 
void read(IsTy & is) { is>>m_pos>>m_mate>>m_hit; Fix();
// is should be a line input stream, reading more should be ok////

 } 
StrTy dump(const StrTy&  sep=" ") const 
{
SsTy x;
x<<m_pos<<sep<<m_mate<<sep<<m_hit<<sep<<m_dist<<sep<<m_length;
return x.str();
}

// fudding unsigned fudd 
// als a fudd when they overlap or are is a pallindrome 
IntTy spacing() const { 
if ( m_mate>m_pos) return IntTy(m_mate)-IntTy(m_pos+m_length); 
return IntTy(m_pos-m_mate)-IntTy(m_length); 
}


// expected distance between hits in random sequence 
IdxTy rna_dist() const { return 1<<(m_length<<1); }
//  return true if not closer together than n time chance distance 
// avoid division 
bool noise(const IdxTy n ) const { return spacing()>(n*rna_dist()); } 
// position of hit, mate and cacilated distance and length
IdxTy m_pos, m_mate,m_length;
IntTy m_dist;
// actual hit contents
StrTy m_hit;

private:
void Fix()
{

m_dist=IdxTy(m_mate)-IdxTy(m_pos) ; // +0*m_length; // note that there is a convention issue here
m_length=m_hit.length(); // sometimes read in from the file but this is better

}



}; // feature_index_entry


class match_criterion
{
public:


void read(IsTy & is) 
{ is>>m_type>>m_first>>m_second>>m_min_dist>>m_max_dist; 
is>>m_min_size>>m_max_size;
is>>m_min_mate>>m_max_mate;
is>>m_symbol;
Fix();

}
StrTy dump(const StrTy&  sep=" ") const 
{
SsTy x;
x<<m_type<<sep<<m_first<<sep<<m_second<<sep<<m_min_dist<<sep<<m_max_dist;
x<<sep<<m_min_size<<sep<<m_max_size<<sep<<m_min_mate<<sep<<m_max_mate;
return x.str();

}
// this really should 
IdxTy m_first,m_second;
IntTy  m_min_dist,m_max_dist;

IdxTy m_min_size,m_max_size;
// nomrally these probably are not used but good for
// the rule in/out
int  m_min_mate,m_max_mate;
StrTy m_type,m_symbol;
private:
void Fix()
{


}


 }; // rule to meet to make sequence a hit  

class concise_rule
{

public:
IdxTy m_start, m_rule;
//IdxTy m_end; 
StrTy m_hit;

}; // consise_rule


// this is a list of criteria which form one regular expression
class criteria : public  std::vector<match_criterion> 
{
typedef value_type Criterion;
typedef criteria Criteria;

public:
void make_groups() const { ItemGroups(m_list,m_last,*this); }
// these are the same 
const IdxTy type(const IdxTy i ) const { return m_list[i]; }
const IdxTy list(const IdxTy i ) const { return m_list[i]; }
const IdxTy last(const IdxTy i ) const { return m_last[i]; }
private:
typedef std::vector<IdxTy> Lists;
// well, need to fix this junk
mutable Lists m_list,m_last;
// probably should merge f and r before... 
// now mx can point to either f or r, need to figure out from criteria
template <class Ty> static void ItemGroups(Ty & list,  Ty & lasts, const Criteria & c )
{
// arrgh
list.clear();
lasts.clear();

// well if there was stuff in there,,,
const IdxTy sz=c.size();
const IdxTy lsz=list.size();
std::map<StrTy , IdxTy> xlate;
// these become zed and unity 
xlate["rule"]=2;
xlate["rc"]=1;
if ( lasts.size()<xlate.size()) lasts.resize(xlate.size()); 

for (IdxTy i=0; i<sz; ++i)
{
const Criterion & ci=c[i];
const IdxTy group=ci.m_first;
// should be consistent
const StrTy & t=ci.m_type;
const IdxTy val=xlate[t];
const IdxTy lse=list.size()-lsz;
if (group==lse)
{
if ( val==0)  {MM_MSG(" invalid criterion tyhpe "<<t<<" see source code")}
const IdxTy vale=val-1;
list.push_back(vale);
lasts[vale]=group+1; //make zero a valid one too sand just comp < size
} else  if ( lse<group) 
{ MM_MSG(" missing criteria at "<<i<<" first "<<group<<" type "<<t) } 
else if  ( (val-1)!=list[group+lsz]) 
{MM_MSG(" inconsistent type "<<t<<" at "<<i<<" for "<<group)}
while ( list.size()<group) list.push_back(0);
}
}



}; // criteria

// this is a collection of criteria which can be scanned to 
// restrict fatures and rules for later testing.
// the idea here is to make a better way to find things
// like mating pairs without having to find each one on each
// pass, perhaps iteratively compiling features and criteria

class criteria_collection : public std::vector<criteria> 
{




}; // criteria_collection




// contains the info on each hit...
// this should be able to exist without later knowledge of
// criteria or features but it is still easier to use with them.
// the matching classes that just kept indicies had to look at
// criteria and features/rules to figure out the matching text.
class hit_element
{
typedef hit_element Myt;
// if this element needs to be popped, resume scan from this state;
public:
typedef std::vector<IdxTy> StateVector; 
hit_element():m_group(BAD),m_start(BAD),m_length(BAD),m_end(BAD) {}
hit_element(const IdxTy i ):m_group(i),m_start(BAD),m_length(BAD),m_end(BAD) {}
hit_element(const IdxTy g, const IdxTy s, const IdxTy len)
	:m_group(g),m_start(s),m_length(len),m_end(s+len) {}

hit_element(const IdxTy g, const IdxTy s, const IdxTy len, const StrTy & hit,
const StrTy & sym)
	:m_group(g),m_start(s),m_length(len),m_end(s+len),m_hit(hit),m_symbol(sym)
 {}


Myt & operator=( const feature_index_entry & fi)
{
m_start=fi.m_pos; 
m_hit=fi.m_hit;
m_length=m_hit.length(); 

return *this;

}

Myt & operator=( const concise_rule & cr)
{

m_start=cr.m_start; 
m_hit=cr.m_hit;
m_length=m_hit.length(); 
return *this;

}


IdxTy m_group; 
IdxTy m_start, m_length,m_end;
StrTy m_hit;
StrTy m_symbol;

// scanner state when added, scanner defined.
StateVector m_state;


}; // hit_element


class hit_struct  : public std::vector<hit_element>
{

typedef value_type E;

public:

int operator()(const E & r, const E & l) { return r.m_start<l.m_start; } 
// the only reason to do was to find extent, now a waste lol. j


void dump(StrTy & pieces, StrTy & groups,const StrTy & label, const StrTy &clabel="piecepos ", const StrTy & glabel="groupsym ") const 
{
//const StrTy label="rclabel ";
const IdxTy len=extent();
const IdxTy sz=size();
MM_MSG(" printing "<<len<<" "<<sz)
if (sz==0 ) return ;
const StrTy sep=" ";
const IdxTy base=start(); // (*this)[0].m_start;
ChTy c[len+1],s[len+1];
c[len]=0; s[len]=0;
for (IdxTy i=0; i<len; ++i) {c[i]='.'; s[i]='.'; }
IdxTy last=0;
//IdxTy offsets[sz];
for (IdxTy i=0; i<sz; ++i)
{
const E & e=(*this)[i];
const IdxTy ix=e.m_start-base;
//offsets[i]=ix;
const IdxTy lxx=e.m_length;
const IdxTy lf=lxx*sizeof(ChTy);
// no longer sorted by start order. 
MM_MSG(" info "<<label<<sep<<i<<sep<<e.m_group<<sep<<e.m_start<<sep<<(IntTy(ix)-IntTy(last))<<sep<<(e.m_start%3)<<sep<<ix<<sep<<e.m_length<<sep<<base)
if ((ix+lxx)>len) { MM_MSG(" overrun "<<ix<<" "<<lxx<<" "<<len)}
if (IntTy(ix)<0) { MM_MSG(" underrrun "<<IntTy(ix)<<" "<<lxx<<" "<<len)}
::memcpy(c+ix, e.m_hit.c_str(),lf);
// do not fill with zed.. 
if ( e.m_symbol.length()!=0)
	::memset(s+ix, e.m_symbol.c_str()[0],lf);

last=ix+e.m_length;
} //i 
pieces=clabel+StrTy(c);
groups=glabel+StrTy(s);



}

void sort()
{
//MM_MSG("sort stopped bad idea")
//std::sort(begin(),end(),*this); 
if ( size()==0) { m_start=0; m_extent=0; return ; }
const IdxTy s=(*this)[0].m_start;
IdxTy ss=s;
IdxTy e=s;
for ( IdxTy i=0; i<size(); ++i) 
{
const E & f=(*this)[i];
const IdxTy ssi=f.m_start;
if ( ssi<ss) ss=ssi;
const IdxTy ei=f.m_start+f.m_length;
if ( ei>e) e=ei;
//MM_MSG(" SHOTFUDD "<<i<<" "<<f.m_start<<" "<<f.m_length<<" "<<f.m_hit<<" "<<e<<" "<<s)

}
m_start=ss;
m_extent=1+e-ss;
//MM_MSG(" FFFUUUCCCKKK AAAS FFFUKCKC "<<e<<" "<<s)
}
// this does not work.. 
//IdxTy extent() const { if (size()==0) return 0; 
//return (*this)[size()-1].m_end-(*this)[0].m_start;}
IdxTy start() const { return m_start; } 
IdxTy extent() const { return m_extent; } 
IdxTy m_start,m_extent;

}; // hit_struct;

typedef hit_element HitElement;
// this is essentially one match to a regex
typedef hit_struct HitStruct;

// this is a collection of matches to a given regex ( criteria set)
class hit_list : public std::vector<hit_struct>
{



} ; // hit_list

class hit_format
{
public:
hit_format() : m_sep(" ") {}
const StrTy & sep () const { return m_sep;}
const StrTy & label () const { return m_label;}

StrTy m_sep;
StrTy m_label;
}; // hit_format

typedef hit_format HitFormat;

typedef HitElement::StateVector ScanState;

// fudding names space wont fudding tyedef
class dummy {
public:

typedef concise_rule ConciseRule;
typedef std::vector<ConciseRule> ConciseRules;

typedef feature_index_entry Feature;
typedef match_criterion Criterion;

// right now just fit all in memory 
typedef std::vector<Feature> Features;
// these should be sorted in order of first
// candidate they qualify so that it can be validated
// in a contiguous sequence. 
//typedef std::vector<Criterion> Criteria;
typedef criteria  Criteria;



// indexes into Features, probabluy should just copy features but
// for now this may be easier.   
typedef std::vector<IdxTy> Matches;

// collection of Match sets that meet criteria
typedef std::vector<Matches> Matching;
//  reorganize rules to fit index. Usually they will occur
// as distance between features, change to range of dists
// between mates first as these are indexed.

typedef hit_list HitList;

typedef mjm_csv::FileTy RegexRules;
typedef mjm_csv::LineTy RegexRule;

static void make_concise_rule_hits(ConciseRules & cr, const RegexRules & f, const IdxTy seq)
{

const IdxTy sz=f.size();
const IdxTy col_seq=f.column(StrTy("sequence")); 
const IdxTy col_start=f.column(StrTy("start")); 
const IdxTy col_hit=f.column(StrTy("hit")); 
const IdxTy col_rule=f.column(StrTy("rule")); 
const bool check_seq=(col_seq!=~0U);
const IdxTy min=col_rule+1;
MM_MSG(" using columns "<<col_seq<<" "<<col_start<<" "<<col_hit<<" "<<col_rule)

for ( IdxTy i=0; i<sz; ++i)
{
const RegexRule&  rr=f[i];
if (rr.size()<min) MM_MSG(" have a rule too short "<<i<<" size is "<<rr.size());

if ( check_seq) if ( IdxTy(::atoi(rr[col_seq].c_str()))!=seq) continue;
ConciseRule r;
r.m_start=::atoi(rr[col_start].c_str());
r.m_rule=::atoi(rr[col_rule].c_str());
r.m_hit=(rr[col_hit]);

cr.push_back(r);


} //i 
MM_MSG(" retainred "<< cr.size()<<" concise rules")

}


static void reorganize_rules(Criteria & d, const Criteria & s)
{



}
// make sure that rules as stated obey assumptions on order 
static void assert_rules(const Criteria & s)
{



}

template<class Ty> static void load_things(Ty & f, IsTy & is,const bool dump=true)
{
typedef mjm_csv::mjm_csv_ops Tx;
enum { MAXLINE=1<<16 } ;
char line[MAXLINE];
while  ( Tx::ok(is))
{
// need to move to a csv class that returns string stream lines 
is.getline(line,MAXLINE-1);
SsTy ss(line);
if (ss.str().c_str()[0]=='#') continue;
// this init crap always bothers me...
typename Ty::value_type fi;
//if ( Ty::ok(is))
fi.read(ss);
if (!Tx::ok(is)) continue;
 f.push_back(fi);
if ( dump ) MM_MSG(fi.dump());
}

}



enum { NONE=~0U};
// lines of pos,mate,hit without dist or lengths
// the csv or other reader should return a stringstream per line
static void load_features_from_string_test(Features & f, IsTy & is, const IdxTy dump=0, const IntTy max_mate=NONE, const IdxTy n_chance=NONE)
{
IdxTy skipped=0;
typedef mjm_csv::mjm_csv_ops Ty;
// if bit 1 is set, the following parameters remove hits of length
// less than mate at distance n_chance or more 
enum { MAXLINE=1<<16,DUMP=1,LEN_DIST_LIM=2 } ;
const bool do_dump=((dump&DUMP)!=0);
MM_MSG("fudding dump "<<dump<<" "<<do_dump<<" "<<max_mate<<" "<<n_chance)
char line[MAXLINE];
while  ( Ty::ok(is))
{
// these need to be somewhere for common processing like comment
// removeal
is.getline(line,MAXLINE-1);
// for some reason this was taken out? 
if  (! Ty::ok(is)) break;
SsTy ss(line);
if (ss.str().c_str()[0]=='#') continue;
// this init crap always bothers me...
Feature fi;
//if ( Ty::ok(is))
fi.read(ss);
// thiswas another assfudd from unsigned fudd shot
if ( max_mate!=NONE) if ( fi.spacing()>max_mate){ 
//MM_MSG(" AFF FUDD "<<fi.spacing()<<" shot "<<max_mate)
++ skipped;  continue;}
if ( max_mate!=NONE) if ( fi.spacing()<(-IntTy(max_mate))){ 
//MM_MSG(" AFF FUDD  shut "<<fi.spacing()<<" shot "<<(-max_mate))

++ skipped;  continue;}
if ( n_chance!=NONE) if ( fi.noise(n_chance)){ ++ skipped;  continue;}

 
f.push_back(fi);
if ( do_dump) MM_MSG(fi.dump());
}

if ( skipped ) MM_MSG(" skipped "<<skipped<<" features due to restrictions")

}





enum { FITS=0, TRY_AGAIN=1, GIVE_UP=~0U,MATCHES=2, PASSES=3};
// this msut return fit, try again, give up 
static IdxTy fits(const Matches & mx, const Criteria & c,  IdxTy  & cptr, IdxTy & fptr, const Features & f)
{
const IdxTy pos=mx.size();
const IdxTy csz=c.size();
// need to do something 
if ( csz==0 ) return GIVE_UP;
while ( cptr<csz)
{
const Criterion & cn=c[cptr];
if  ( cn.m_first!=pos) { return FITS;}
const Feature & fi=f[fptr];
// second has to be less than or equal, if equal point to self
const bool self_test=(cn.m_first==cn.m_second);
const Feature & fi2=self_test?fi:(f[mx[cn.m_second]]);
// rules are ordered 
// SGINED
const int  dis=(int(fi.m_pos)-int(fi2.m_pos));
//MM_MSG(" distance "<<dis<<" must be from "<<cn.m_min_dist<<" TO "<<cn.m_max_dist)
if ( dis>int(cn.m_max_dist)) return GIVE_UP;
// the coding here is really dumb, hopefully compiler can fix this..
if ( dis<int(cn.m_min_dist)) return TRY_AGAIN;
// then need to make sure it is the right kind of feature
const IdxTy len=fi.m_length;
// in general this will be the most frequent killer
if ( len<cn.m_min_size) return TRY_AGAIN;
if ( len>cn.m_max_size) return TRY_AGAIN;
// should know what this mates with
const IntTy mate_dist=fi.m_dist;
if ( mate_dist<cn.m_min_mate) return TRY_AGAIN;
if ( mate_dist>cn.m_max_mate) return TRY_AGAIN;
++cptr;
}

// all the crieria were met so we are done
return MATCHES;
}

static IdxTy fit(const Criterion & cn, const ConciseRule & ri, const ConciseRule & ri2)
{
//MM_MSG("rule vs rule")
const int  dis=(int(ri.m_start)-int(ri2.m_start));
//MM_MSG("rule vs rule"<<dis)
if ( dis>int(cn.m_max_dist)) return GIVE_UP;
if ( dis<int(cn.m_min_dist)) return TRY_AGAIN;
// right now ignore length although there could be other criteria.

// there is no mate 
const IdxTy rule=ri.m_rule;
//MM_MSG(" testing "<<rule<<" vs "<<cn.m_min_mate<<" "<<cn.m_max_mate)
if ( rule<(cn.m_min_mate)) return TRY_AGAIN;
if ( rule>(cn.m_max_mate)) return TRY_AGAIN;
//MM_MSG(" passes rule")
return PASSES;
}

static IdxTy fit(const Criterion & cn, const ConciseRule & ri, const Feature & fi2, const bool nswap )
{

 int  dis=(int(ri.m_start)-int(fi2.m_pos));
if ( !nswap ) dis=-dis;
//if ( !nswap) MM_MSG(" f-r testing rx="<<ri.m_start<<" fi2="<<fi2.m_pos<<" "<<dis<<" vs "<<cn.m_min_dist<<" "<<cn.m_max_dist)
if ( dis>int(cn.m_max_dist)) return GIVE_UP;
if ( dis<int(cn.m_min_dist)) return TRY_AGAIN;
// do not bother with mate or other stuff, self rule should catch that 
// unless the self rule was omitted ...
// probably better to add rule range check here than require self rule.

//MM_MSG(" tested rx="<<ri.m_start<<" fi2="<<fi2.m_pos<<" "<<dis<<" vs "<<cn.m_min_dist<<" "<<cn.m_max_dist)

return PASSES;
}




// in a sefl-test, fi2 is fi so they have zero offset
static IdxTy fit(const Criterion & cn, const Feature & fi, const Feature & fi2)
{
const int  dis=(IntTy(fi.m_pos)-IntTy(fi2.m_pos));
//MM_MSG(" distance "<<dis<<" must be from "<<cn.m_min_dist<<" TO "<<cn.m_max_dist)
// for self test, zed distance  assuerd that fi2 occurs not after fi,
// use neg distance as a code...
const IntTy  mx=cn.m_max_dist;
const IntTy  mn=cn.m_min_dist;
// special code for this, these two
// must be maters
if (mn<0)
{
// they need to point to each other and have same string.
//MM_MSG(" comparing "<<fi.m_pos<<" "<<fi.m_mate<<" "<<fi2.m_mate<<" "<<fi.m_hit)
if ( fi.m_pos!=fi2.m_mate) return TRY_AGAIN;
if ( fi2.m_pos!=fi.m_mate) return TRY_AGAIN;
if ( fi2.m_length!=fi.m_length) return TRY_AGAIN;
// crap, these are rc, just check len lol 
//if (fi.m_hit!=fi2.m_hit) return TRY_AGAIN;
}
else
{
if ( dis>mx) return GIVE_UP;
// the coding here is really dumb, hopefully compiler can fix this..
if ( dis<mn) return TRY_AGAIN;
}

// then need to make sure it is the right kind of feature
const IdxTy len=fi.m_length;
//MM_MSG(" wtf len "<<len<<" "<<cn.m_min_size<<" "<<cn.m_max_size)
// in general this will be the most frequent killer
if ( len<cn.m_min_size) return TRY_AGAIN;
if ( len>cn.m_max_size) return TRY_AGAIN;
// should know what this mates with
const IntTy  mate_dist=fi.m_dist;
//MM_MSG(" wtf mate "<<mate_dist<<" "<<cn.m_min_mate<<" "<<cn.m_max_mate)
if ( mate_dist<cn.m_min_mate) return TRY_AGAIN;
if ( mate_dist>cn.m_max_mate) return TRY_AGAIN;
return PASSES;
}

// go though all criteria that need to be met to accept these
// element into hs 
// the rule type , feature pointer, rule pointer, and criteria point
// define the scan state
enum { RT=0, FP=1, RP=2, CP=3};
static IdxTy fits(const HitStruct & hs, ScanState & s, const Criteria & c,  const Features & f, const ConciseRules & r)
{
const IdxTy pos=hs.size();
const IdxTy csz=c.size();
//MM_MSG(" fitt pos="<<pos<<" csz="<<csz<<" s[CP]="<<s[CP])
// need to do something 
if ( csz==0 ) return GIVE_UP;

while ( s[CP]<csz)
{
const Criterion & cn=c[s[CP]];
// these are the groups to which the criteria apply,
// second must be in the list alrady 
const IdxTy g1=cn.m_first;
//MM_MSG(" fitt "<<pos<<" "<<g1<<" "<<s[CP])
if  ( g1!=pos) { return FITS;}
const IdxTy g2=cn.m_second;
const bool self_test=(g1==g2);
// this may not work rtight 
//const IdxTy type1=s[RT]; // c.list(g1); // tv[g1];
const IdxTy type1= c.list(g1); 
const IdxTy type2= c.list(g2); // self test  fudd hs[g2].m_state[RT];  // this SHOUOD work lol 
const IdxTy br=type1+(type2<<1);
//MM_MSG(" fitt "<<g1<<" "<<g2<<" "<<s[CP]<<" "<<self_test)

switch (br)
{
case 0:{
const Feature & fi=f[s[FP]];
const Feature & fi2=self_test?fi:(f[hs[g2].m_state[FP]]);
const IdxTy fitt=fit(cn,fi,fi2);
if ( fitt!=PASSES) return fitt;
break;
}
case 1:
{
//MM_MSG(" case 1, fptr="<<s[FP]<<" rptr="<<s[RP])
if (self_test) MM_MSG( " self test is true for "<<s[CP])
const ConciseRule & ri=r[s[RP]];
const Feature & fi2=(f[hs[g2].m_state[FP]]);
const IdxTy fitt=fit(cn,ri,fi2,true);
if ( fitt!=PASSES) return fitt;
//MM_MSG(" case 1 passes "<<s[FP]<<" "<<s[RP]<<" "<<s[CP]<<" "<<csz)
break;
}

case 2:
{
//MM_MSG(" case 2, fptr="<<s[FP]<<" rptr="<<s[RP])
//if (self_test) MM_MSG( " self test is true for "<<s[CP])
const Feature & fi=f[s[FP]];
const ConciseRule & ri2=(r[hs[g2].m_state[RP]]);
const IdxTy fitt=fit(cn,ri2,fi,false); // saves a signature, 
if ( fitt!=PASSES) return fitt;
//MM_MSG(" case 2 passes "<<s[FP]<<" "<<s[RP]<<" "<<s[CP]<<" "<<csz)
break;
}

case 3:
{
//MM_MSG(" case 3 "<<s[RP]<<" "<<r.size()<<" "<<self_test)
const ConciseRule & ri=r[s[RP]];
const ConciseRule & ri2=self_test?ri:(r[hs[g2].m_state[RP]]);
//MM_MSG(" case 3 "<<rptr<<" "<<self_test)
const IdxTy fitt=fit(cn,ri,ri2);
//MM_MSG(" case 3 "<<rptr<<" "<<self_test)
if ( fitt!=PASSES) return fitt;
break;

}

default : MM_MSG(" bad case "<<br )

} // switch 

// rules are ordered 
// SGINED
++s[CP];

} // cptr

// all the crieria were met so we are done
return MATCHES;
}


// I suppose that since the list is ordered, a binary search would
// make sense but it is not likely to heplp much esp with memory
// thrashing although ti depends on getting end pts right . 
static IdxTy fits(const Matches & mx, const Criteria & c,  IdxTy  & cptr, IdxTy & fptr, const Features & f, IdxTy & rptr, const ConciseRules & r)
{
const IdxTy pos=mx.size();
const IdxTy csz=c.size();
// need to do something 
if ( csz==0 ) return GIVE_UP;

// this needs to be part of Criteria 
typedef std::vector<IdxTy> Types;
// return the last location of each type in lasts

// terminating will be a little trickier
while ( cptr<csz)
{
const Criterion & cn=c[cptr];
// this criteria is wrong, we should probably keep checking for other
// but assume at least one for each and in order
const IdxTy g1=cn.m_first;
if  ( g1!=pos) { return FITS;}
const IdxTy g2=cn.m_second;
const bool self_test=(g1==g2);
//if (mx.size()>2) MM_MSG(cptr<<" groups g1="<<g1<<" "<<"g2="<<g2)
const IdxTy type1=c.list(g1); // tv[g1];
const IdxTy type2=c.list(g2); // tv[g2];
const IdxTy br=type1+(type2<<1);

// second has to be less than or equal, if equal point to self
switch (br)
{
case 0:{
const Feature & fi=f[fptr];
const Feature & fi2=self_test?fi:(f[mx[g2]]);
const IdxTy fitt=fit(cn,fi,fi2);
if ( fitt!=PASSES) return fitt;
break;
}
case 1:
{
//MM_MSG(" case 1, fptr="<<fptr<<" rptr="<<rptr)
//if (self_test) MM_MSG( " self test is true for "<<cptr)
const ConciseRule & ri=r[rptr];
const Feature & fi2=(f[mx[g2]]);
const IdxTy fitt=fit(cn,ri,fi2,true);
if ( fitt!=PASSES) return fitt;
break;
}

case 2:
{
//MM_MSG(" case 2, fptr="<<fptr<<" rptr="<<rptr)
//if (self_test) MM_MSG( " self test is true for "<<cptr)
const Feature & fi=f[fptr];
const ConciseRule & ri2=(r[mx[g2]]);
const IdxTy fitt=fit(cn,ri2,fi,false); // saves a signature, 
if ( fitt!=PASSES) return fitt;
break;
}





case 3:
{
//MM_MSG(" case 3 "<<rptr<<" "<<r.size()<<" "<<self_test)
const ConciseRule & ri=r[rptr];
const ConciseRule & ri2=self_test?ri:(r[mx[g2]]);
//MM_MSG(" case 3 "<<rptr<<" "<<self_test)
const IdxTy fitt=fit(cn,ri,ri2);
//MM_MSG(" case 3 "<<rptr<<" "<<self_test)
if ( fitt!=PASSES) return fitt;
break;

}

default : MM_MSG(" bad case "<<br )

} // switch 



// rules are ordered 
// SGINED
++cptr;
}

// all the crieria were met so we are done
return MATCHES;
}

static bool dunn(const ScanState & s, const Criteria & c, const IdxTy szm, const IdxTy szf, const IdxTy szr)
{
const bool fp_dun=(s[FP]>=szf);
const bool rp_dun=(s[RP]>=szr);
// if we need more but are out of then dun 
const bool dun=((szm<c.last(0))&&(fp_dun))||((szm<c.last(1))&&(rp_dun));
return dun;

}

// I suppose a state iterator rather than limit on inc?
static bool inc(ScanState & s)
{
if ( s[RT]==0) { ++s[FP];} else if (s[RT]==1) { ++s[RP];}
else {{MM_MSG(" FITS errors in rt "<<s[RT]<<" "<<" ") } return false;}
return true;
} //
static bool assign(HitElement & he, const ScanState &s, const IdxTy cptr_start,const Criteria & c, const Features & f, const ConciseRules & r )
{
const IdxTy hsz=he.m_group;
if ( s[RT]==0) {he=f[s[FP]]; }
else if ( s[RT]==1) { he=r[s[RP]]; } 
else {{MM_MSG(" type error "<<s[RT]<<" "<<hsz); } return false; } 
he.m_symbol=c[cptr_start].m_symbol;
ScanState &  si=he.m_state;
si=s;
si[CP]=cptr_start;

return true;
}

static void match(HitList & hl,  const Criteria & c, const Features & f, const ConciseRules & r )
{


IdxTy szf=f.size();
IdxTy szr=r.size();
if ( c.size()==0) return ; 
// fix this crap
c.make_groups();
std::vector<ScanState> restart;
// this should take a ctor with sizes and then have
// inc method. 
ScanState s;
s.resize(4); // rt,fp,rp,cptr
s[RT]=c.list(0);  // this only needs to be valid for inc but...
HitStruct hs;
IdxTy next=TRY_AGAIN; // need to check give up first.. 
while ( true)
{

const IdxTy szm=hs.size();
//const bool fp_dun=(s[FP]>=szf);
//const bool rp_dun=(s[RP]>=szr);
// if we need more but are out of then dun 
//const bool dun=((szm<c.last(0))&&(fp_dun))||((szm<c.last(1))&&(rp_dun));
const bool dun=dunn(s,c,szm,szf,szr);
//MM_MSG("trying dun "<< next<<" "<<s[FP]<<" "<<s[RP]<<" "<<dun<<" "<<szm)
if (( next==GIVE_UP) ||(dun))
{
//MM_MSG("give up test"<<hs.size())
if ( szm==0) {{ MM_MSG(" scan completes on "<<next<<" "<<dun)} break;} // nothing left to try 
next=TRY_AGAIN; // unless popping pops
s= hs.back().m_state;
//MM_MSG(" szm = "<<szm<<" restart size "<<restart.size())
//s=restart[szm-1]; 
// this elements requires these and we are popping so no need to check 
// could use inc but then need to figure out how to handle boundary.
// putting szf etc into ScanState ctor may work. 
if (!inc(s)) break;
// only one of these tests is needed, should add logic to state
// but only needed in a few places, perhaps iterator?
if (s[FP]>=szf) next=GIVE_UP;
if (s[RP]>=szr) next=GIVE_UP;
hs.pop_back();
} // GIVE_UP
IdxTy cptr_start=s[CP]; // arrgghhhh fix this 
if ( next!=GIVE_UP) next= ( fits(hs,s,c,f,r)) ;

if ( next==FITS) 
{ 
// really all of this except for state pushing can be defered
// until we decide to keep a complete hit, but it is
// a possible debug help
const IdxTy hsz=hs.size();
s[RT]=c.list(hsz); // the number of criter caan not be zed
HitElement he(hsz);
if ( !assign(he,s,cptr_start,c,f,r)) break ;
hs.push_back(he);
if ( restart.size()<=hs.size()) restart.resize(hs.size()+1);
// the older size, could be zed
restart[hsz]=s;

//if ( hs.size()>3) MM_MSG(" have size now "<<hs.size()<<" cp="<<fudd<<" "<<c.size())

// this inc is bad  in theory anyway same element of differnt stack etc. 
//inc(s);

s[RT]=c.list(hs.size()); // the number of criter caan not be zed
} //fits
// we have a complete set, save and restart scanning by
// only popping the last element although the whole
// stack could be cleared out. 
if ( next==MATCHES) 
{ 
const IdxTy hsz=hs.size();
s[RT]=c.list(hs.size()); // the number of criter caan not be zed
HitElement he(hsz);
if ( !assign(he,s,cptr_start,c,f,r)) break ;
hs.push_back(he);
if ( restart.size()<=hs.size()) restart.resize(hs.size()+1);
// the older size, could be zed
restart[hsz]=s;
// probably not matter here 
s[RT]=c.list(hs.size()); // the number of criter caan not be zed
hs.sort();
hl.push_back(hs); 
//MM_MSG(" have size now for match "<<hs.size()<<" and "<< hl.size())
next=GIVE_UP; // take this one and try to resume scanning 
} // MATCHES
// this element did not pass all criteria, but there are 
//more candidates 
if ( next==TRY_AGAIN) 
{ 
// to be cleaner, this should only be done on init and when size changes
//s[RT]=c.list(hs.size()); // the number of criter caan not be zed
// try again for this level also means too early so
// it can be ruled out for next pass,
// this is off one level doh.
inc(s);
s[CP]=cptr_start;  

}

//MM_MSG("trying loop "<< next<<" rt="<<s[0]<<" cp="<<s[CP])
} // true

} // match hitstruct




static void match(Matching & m, const Criteria & c, const Features & f, const ConciseRules & r )
{
const bool olap=false;
//MM_MSG("olap turned off for now note polarity os backwards")
IdxTy szf=f.size();
IdxTy szr=r.size();
IdxTy fp=0; IdxTy rp=0;
IdxTy cptr=0; // critera_pointer
// assumed later
if ( c.size()==0) return ; 
// right now this is only ints need criteria to know if these
//Types tv,lasts;
// sort the elements into being defined by features or rules
// this is really stpud should just merge the features and ruls
//item_groups(tv,lasts,c);
c.make_groups();

// are rules or rc
Matches mx; // tentative matches so far 


Matches cptrs;
// this needs a struct lol 
Matches optrs;
IdxTy x=TRY_AGAIN; // need to check give up first.. 
while ( true)
{
//if (mx.size()>2) MM_MSG(" another pass "<<mx.size())
// the give_up test needs to be done first, and teh value of x retained...
{
// this is probably missing cases since the test for 
// szr should not be needed, probably needed for sxf too 
// danger red alert todo 
//const bool fp_dun=(fp>=szf)&&(szf!=0);
const bool fp_dun=(fp>=szf);
// something was wrong, need to find out why the test is needed here. 
//const bool rp_dun=(rp>=szr)&&(szr!=0);
const bool rp_dun=(rp>=szr);
const IdxTy szm=mx.size();
// if we need and have exhauseted this is dun;
// if there are none, the lasts will be zed, will this work right?
// the last() value is now designed to work this way 
// see if this is now missing cases at the end.. 
//const bool dun=((szm<=c.last(0))&&(fp_dun))||((szm<=c.last(1))&&(rp_dun));
// if the sz is greater than both, we need nothing and should be dun lol
// but then give_up should be set. 
const bool dun=((szm<c.last(0))&&(fp_dun))||((szm<c.last(1))&&(rp_dun));

//if ( szm>2) MM_MSG(" almost thre "<<szm<<" "<<c.last(0)<<" "<<c.last(1) <<" "<<fp_dun<<" "<<rp_dun<<" "<<fp<<" "<<rp<<" "<<x)

if (( x==GIVE_UP) ||(dun))
{
x=TRY_AGAIN; // unless popping pops
if ( szm==0)
{
MM_MSG(" breaking out now")
 break; // nothing left to try 

}
 //MM_MSG(" giving up no "<<szm<<" "<<dun<<" "<<rp<<" "<<fp)
// so how do we back up now with the other ptr? 
IdxTy rt=c.list(mx.size()-1); // the number of criter caan not be zed
if ( rt==0) {
//MM_MSG(" popping feature "<<szm)
fp=mx.back()+1; 
if (fp>=szf) x=GIVE_UP;
rp=optrs.back();
}
else if ( rt==1){
//MM_MSG(" popping rule "<<szm)
 rp=mx.back()+1;
if (rp>=szf) x=GIVE_UP;
fp=optrs.back(); 
}else 
MM_MSG(" popping ERROR  "<<szm<<" "<<rt)

mx.pop_back(); optrs.pop_back(); cptr=cptrs.back(); cptrs.pop_back();

} // give up 



} // give up block 
IdxTy cptr_start=cptr; // arrgghhhh fix this 

// this will fail if sizes are zero, fp and rp need to be checked..
 if ( x!=GIVE_UP) x= ( fits(mx,c,cptr,fp,f,rp,r)) ;
//if ( x!=1) MM_MSG(" fits returns "<<x<<" mx size "<<mx.size()<<" "<<fp<<" "<<rp)
//if ( mx.size()>2) MM_MSG(" fits returns "<<x<<" mx size "<<mx.size()<<" "<<fp<<" "<<rp)
if ( x==FITS) 
{ 
//MM_MSG(" a fit at "<<mx.size())
IdxTy rt=c.list(mx.size()); // the number of criter caan not be zed
if ( rt==0 ) { mx.push_back(fp); ++fp; optrs.push_back(rp); }
else if ( rt==1 ) { //MM_MSG(" have hit rule ") 
mx.push_back(rp); ++rp; optrs.push_back(fp);}
else MM_MSG(" bad rt="<<rt)
cptrs.push_back(cptr_start); 
} //fits
if ( x==MATCHES) 
{ 
IdxTy rt=c.list(mx.size()); // the number of criter caan not be zed
if ( rt==0 ) { mx.push_back(fp); optrs.push_back(rp); }
else if ( rt==1 ) { mx.push_back(rp); optrs.push_back(fp); }
// matter in olap case
cptrs.push_back(cptr_start);
m.push_back(mx); 


if (olap)
{ 
// eh we needed to push this new one doh.
IdxTy rtz=c.list(0); // the number of criter caan not be zed
if ( mx.size()!=0)
{
// this causes hits with the same first location to only
// pickup the first one, should only bump on failure of 
// higher ups, that is all the others should GIVE_UP befoire
// comding down here,just remove the last one and keeptrying . 
if ( rtz==0)  {fp=mx[0]+1; rp=optrs[0];}
else if ( rtz==1) { rp=mx[0]+1; fp=optrs[0]; } 

}
 else if ( rtz==0) ++fp; else ++rp;


   mx.clear();cptr=0; cptrs.clear(); optrs.clear();
} // olap
else x=GIVE_UP; // just pop the last one 

 }  //MATCHES

 
if ( x==TRY_AGAIN) 
{ 
IdxTy rt=c.list(mx.size()); // the number of criter caan not be zed
if ( rt==0) ++fp;
else if ( rt==1) ++rp;

cptr=cptr_start;  
}



} // true


} // match


// givena set of match criteria, scan a list of hits an assemble 
// subsets which fullfill all the criteria


// givena set of match criteria, scan a list of hits an assemble 
// subsets which fullfill all the criteria
static void match(Matching & m, const Criteria & c, const Features & f )
{
// need a list ofe features and a ptr
MM_MSG(" matching without rules")
IdxTy sz=f.size();
IdxTy fp=0; // featuer_pointer
IdxTy cptr=0; // critera_pointer
// really the matches object should push back all this
// crap at once including ctpr if that can vary 
Matches mx; // tentative matches so far 
// in reality, these are fixed and known ahead of time but for now
// leave variable 
Matches cptrs;
while ( true)
{

//const Criterion & ci=f[fp];
// read through the features, for each one see if it can be the n-th
// of a new match.
IdxTy cptr_start=cptr; // arrgghhhh fix this 
IdxTy x= ( fits(mx,c,cptr,fp,f)) ;
//MM_MSG(" at "<<fp<<" mxsz="<<mx.size()<<" code="<<x)
if ( x==FITS) { mx.push_back(fp); ++fp; cptrs.push_back(cptr_start); } 
// the size here can not be zero, should check.. 
if ( x==MATCHES) 
{ 
mx.push_back(fp);
//cptrs.push_back(cptr_start);
m.push_back(mx); 
// eh we needed to push this new one doh.
if ( mx.size()!=0) fp=mx[0]+1; else ++fp;
   mx.clear();cptr=0; cptrs.clear();
 } 
// the use of incremented fp means no switch statement here. 
// we need to pop the cptr stack...
if ( x==TRY_AGAIN) { ++fp;cptr=cptr_start;  } 
if (( x==GIVE_UP) ||(fp>=sz))
{
if ( mx.size()==0) break; // nothing left to try 
fp=mx.back()+1; 
mx.pop_back(); 
cptr=cptrs.back();
cptrs.pop_back();

}

// if the first one is possible, then start looking for the 
// second one knowing where first and mates are

// contine for each i-th features or go back to looking for a new
// first one.

} // fp 

}

// these could be const I suppose
static void find_features( Matching & m , Features & f, Criteria & c,ConciseRules & cr)
{
//Fy::Matching m;
// this must use the rules one to pickup the type info 
// I changed that now.. 
//Fy::match(m,c,f);
//if ( cr.size()!=0) 
MM_MSG(" forcing match with rules as better even if no rules present")
match(m,c,f,cr);
//else match(m,c,f);
//print_features(m,c,f,cr);
}


static void print_ex_features( OsTy & os, const HitList & hl, const HitFormat & hf=HitFormat())
{
const StrTy & sep= hf.sep();
const StrTy & rclabel=hf.label();
const IdxTy sz=hl.size();
for ( IdxTy i=0; i<sz; ++i)
{
const HitStruct & hs= hl[i];
if ( hs.size()==0) continue;
SsTy ss;
ss<<rclabel<<sep<<i<<sep;
StrTy p; StrTy g;
hs.dump(p,g,ss.str());
//dump(os, rclabel,mx, f, cr,c);
os<<p<<CRLF;
os<<g<<CRLF;
os<<" offsets "<<i<< "dxj: ";
for ( IdxTy j=1; j<hs.size(); ++j) 
os<<sep<<IntTy(hs[j].m_start-hs[j-1].m_start);
os<<CRLF;

} // i 
} //print_ex


static void print_ex_features( OsTy & os, Matching & m, Features & f, Criteria & c,ConciseRules & cr, const StrTy & rclabel="rclabel")
{
const StrTy sep=" ";
const IdxTy msz=m.size();
MM_MSG(" matches count is "<<msz)
for ( IdxTy i=0; i<msz; ++i)
{
const Matches & mx=m[i];
if ( mx.size()==0) continue;
HitStruct hs;
make_hit_struct(hs,mx,f,cr,c);
SsTy ss;
ss<<rclabel<<sep<<i<<sep;
StrTy p; StrTy g;
hs.dump(p,g,ss.str());
//dump(os, rclabel,mx, f, cr,c);
os<<p<<CRLF;
os<<g<<CRLF;
os<<" offsets "<<i<< "dxj: ";
for ( IdxTy j=1; j<mx.size(); ++j) os<<sep<<IntTy(hs[j].m_start-hs[j-1].m_start);

os<<CRLF;
}

}



static void print_features( OsTy & os, Matching & m, Features & f, Criteria & c,ConciseRules & cr, const StrTy & rclabel="rclabel")
{
//OsTy & os = std::cout;
const IdxTy msz=m.size();
MM_MSG(" matches count is "<<msz)
for ( IdxTy i=0; i<msz; ++i)
{
const Matches & mx=m[i];
dump(os, rclabel,mx, f, cr,c);
os<<CRLF;
}
}

// should use this everywhere early on... 
static void make_hit_struct(HitStruct & hs ,  const Matches & mx, const Features & f, const ConciseRules & r, const Criteria & c)
{

const IdxTy sz=mx.size();
const IdxTy rsz=r.size();
const IdxTy csz=c.size();
const bool include_rules=(rsz!=0);
const bool debug=!true;
if ( debug) MM_MSG(" nclude rules  "<<include_rules)
IdxTy cptr=0;
for ( IdxTy i=0; i<sz; ++i)
{
// only the distance to mate is indicated, not the mate index 
// but need to fix that... 
IdxTy ss=0;
IdxTy lens=0;
StrTy hits="";
if (include_rules)
{
if ( debug) MM_MSG(" checking type for "<<i)
if ( c.type(i)==0)
{
const Feature & fi=f[mx[i]];
ss=fi.m_pos;
lens=fi.m_length;
hits=fi.m_hit;
}

else if ( c.type(i)==1)
{
if ( debug) MM_MSG(" doing a rule at "<<i)
const ConciseRule & cr= r[mx[i]];
ss=cr.m_start;
// really the read needs to be fixed
lens=cr.m_hit.length(); // cr.m_end-cr.m_start+1;
hits=cr.m_hit;
} // 1
else MM_MSG(" bad type "<<c.type(i)<<" at "<<i)
} //include rules
else
{
const Feature & fi=f[mx[i]];
ss=fi.m_pos;
lens=fi.m_length;
hits=fi.m_hit;
}
while ( cptr<csz) 
{
if ( c[cptr].m_first==i) break;
if ( c[cptr].m_first>i) break;
++cptr;
}
HitElement e(i,ss,lens);
e.m_hit=hits;
if ( cptr<csz) if ( c[cptr].m_first==i)  {e.m_symbol=c[cptr].m_symbol;}


hs.push_back(e);
} // fo
//MM_MSG(" sorting hits "<<sz)
// sort based on pos...
hs.sort(); 
// calculate extent , should be dun in sort 

//MM_MSG(" sorting hits done "<<sz<<" "<<hs.extent())

} //make_hit

// dump the matches by looking up what their indexes point at in featuers
// and maybe relate to criteria
// this should be expanded with a fasta files and maybe genbank too
// this should also now allow for out-or-order hits even though
// the scanning right now only supports forward, cross list hits
// can differ in sequence 
static void dump(OsTy & os, const StrTy & label, const Matches & mx, const Features & f, const ConciseRules & r, const Criteria & c)
{
const IdxTy sz=mx.size();
const IdxTy csz=c.size();
const bool debug=false; 
if ( sz==0) return;
const StrTy sep=" ";
const IdxTy rsz=r.size();
const bool include_rules=(rsz!=0);
if ( debug) MM_MSG(" nclude rules  "<<include_rules)
SsTy bases;
SsTy groups;
// size is at least one 
IdxTy start=0;
if ( !include_rules) start=f[mx[0]].m_pos;
else
{
// this should not happen, we should check the size but easy to do.
c.make_groups();
if ( debug) MM_MSG(" dun making ")

start=(c.type(0)==0)?f[mx[0]].m_pos:r[mx[0]].m_start;
}

IdxTy ptr=0;
char cs[2];
IdxTy cptr=0;
for ( IdxTy i=0; i<sz; ++i)
{
// only the distance to mate is indicated, not the mate index 
// but need to fix that... 
IdxTy ss=0;
IdxTy lens=0;
StrTy hits="";
if (include_rules)
{
if ( debug) MM_MSG(" checking type for "<<i)
if ( c.type(i)==0)
{
const Feature & fi=f[mx[i]];
ss=fi.m_pos;
lens=fi.m_length;
hits=fi.m_hit;
}

else if ( c.type(i)==1)
{
if ( debug) MM_MSG(" doing a rule at "<<i)
const ConciseRule & cr= r[mx[i]];
ss=cr.m_start;
// really the read needs to be fixed
lens=cr.m_hit.length(); // cr.m_end-cr.m_start+1;
hits=cr.m_hit;
} // 1
else MM_MSG(" bad type "<<c.type(i)<<" at "<<i)
} //include rules
else
{
const Feature & fi=f[mx[i]];
ss=fi.m_pos;
lens=fi.m_length;
hits=fi.m_hit;
}

MM_MSG(" info "<<label<<sep<<i<<sep<<ss<<sep<<(ss%3)<<sep<<(ss-start)<<sep<<lens<<sep<<hits<<sep<<start)

cs[0]='1'+((i)%10);
cs[1]=0;
while ( cptr<csz) 
{
if ( c[cptr].m_first==i) break;
if ( c[cptr].m_first>i) break;
++cptr;
}
if ( cptr<csz) if ( c[cptr].m_first==i) if ( c[cptr].m_symbol.length()!=0)
{
cs[0]=c[cptr].m_symbol.c_str()[0];

}

 
StrTy sym= StrTy(cs);
// features are in order but mates?
// this works until you have overlapping hits for riboswitch lol 
const IdxTy dif=ss-start;
while ( ptr<dif) { groups<<"."; bases<<"."; ++ptr; }
const IntTy dx=IntTy(dif)-IntTy(ptr);
// this needs to belimited to the lack of overlap
IntTy dd=lens;
if ( dx<0){ dd=lens+dx;
for (IntTy ii=-dx; ii<IntTy(lens); ++ii) bases<<hits.c_str()[ii];
}
else bases<<hits; 

ptr+=dd;

//for ( IdxTy j=0; j<fi.m_length; ++j) groups<<(((i+1)%10));
for ( IntTy j=0; j<dd; ++j) groups<<sym;


}
os<<"piecepos"<<sep<<label<<sep<<start<<sep<<bases.str()<<CRLF;
os<<"groupsym"<<sep<<label<<sep<<start<<sep<<groups.str()<<CRLF;
}


} ; // assfycj how the hell do you typedef a fudding name space?
}; // ns  mjm_features


#endif


