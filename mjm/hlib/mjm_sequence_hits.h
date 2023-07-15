#ifndef MJM_SEQUENCE_HITS_H__
#define MJM_SEQUENCE_HITS_H__

#include "mjm_text_data.h"
#include <stdlib.h>
#include <tcl.h>
//#include "mjm_tcl_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
#include <mjm_csv_ops.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
#include <algorithm>

class sequence_hits_typedefs
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

namespace mjm_sequence_hits 
{
typedef sequence_hits_typedefs Tr;

typedef Tr::IdxTy IdxTy;
typedef  Tr::IntTy IntTy;
typedef Tr::ChTy ChTy;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;

// right now these are just handled as csv files,
// will want to decouple from that 
class hit_class
{
public:

void set(const IdxTy seq, const IdxTy ru, const IdxTy s, const IdxTy e)
{ m_seq=seq; m_rule=ru; m_start=s; m_end=e; }
void set(const StrTy & h, const StrTy &  desc) { m_hit=h; m_desc=desc;}


private:
// these should be const really, 
IdxTy m_seq, m_rule,m_start, m_end;
// hit is the actual matching text.
StrTy m_hit,m_desc;


}; // hit_class
typedef hit_class HitClass;
// need to examine relationshipt between hits to sort out regex
// and debug rules files etc
class hit_pair
{
typedef hit_class Ru; // the hits, could be csv file lines too 
public:

private:

Ru m_x, m_y;
// stats- relative to start of x;
IdxTy m_start_overlap, m_end_overlap;
IdxTy m_overlaps, m_x_singles,m_y_singles;

bool m_x_in_y, m_y_in_x;




}; // hit_pair;


typedef mjm_csv::LineTy LineTy;
// note this is SIGNED 
typedef mjm_csv::IntLineTy IntLineTy;


typedef mjm_csv::FileTy FileTy;
typedef mjm_csv::IndexMapTy  IndexMapTy;
typedef mjm_csv::BinarySiteMapTy BinarySiteMapTy;
typedef mjm_csv::LineIndexMapTy LineIndexMapTy;
typedef mjm_csv::InverseLineIndexMapTy InverseLineIndexMapTy;
typedef mjm_csv::StrIndexMapTy StrIndexMapTy;



// these deal with locations of rule hits in various sequences
// and relations among them  Generally rules are regular expressions
// and hits are locations, rule number, and subsequences that match
class mjm_sequence_hits
{
//
typedef LineTy::const_iterator LI;
typedef FileTy::const_iterator FI;
typedef IndexMapTy::const_iterator II;
typedef LineIndexMapTy::const_iterator LII;
typedef InverseLineIndexMapTy::const_iterator ILII;


public:
typedef mjm_csv::mjm_csv_ops FockTy;
static void dump_map(OsTy & os, const InverseLineIndexMapTy & s, const StrTy & sep)
{
ILII ii=s.begin();
while ( ii!=s.end())
{
FockTy::dump_line(os,(*ii).first,sep,1);
os<<sep<<"|"; 
const FileTy & f=(*ii).second;
FockTy::dump_file(os,f,sep,1);

os<<CRLF;
++ii;
}


}

static void dump_map(OsTy & os, const IndexMapTy & s, const StrTy & sep)
{
II i=s.begin();
while ( i!=s.end())
{
os<<((*i).first)<<sep<<"|";
const LineTy & line= (*i).second;
FockTy::dump_line(os,line,sep); //LI li=line.begin(); while (li!=line.end()) { os<<sep<<(*li); ++li; }
os<<CRLF;
++i; 
}

}


// include the forward map for easer, 
static void make_binary_site_map(BinarySiteMapTy & bm, const InverseLineIndexMapTy & im,const LineIndexMapTy & idx)
{
//first go through each line-key and compile a list of all sites
typedef std::vector< ILII> IV;
typedef IV::const_iterator IVII;
typedef std::map<IntTy,IV >  LoMap;
typedef LoMap::const_iterator LMI;
// this is for something stupid....
const IdxTy seqs=idx.size();
IdxTy io=0;
IdxTy order[seqs];
for ( IdxTy i=0; i<seqs; ++i) order[i]=~0;  // this shouuld be somethign invald
LoMap lomap;
ILII ilii=im.begin();
//MM_MSG("")
while ( ilii!=im.end())
{
const LineTy & line=(*ilii).first;
LI li = line.begin();
while ( li!=line.end())
{
IntTy val=::atoi((*li).c_str()); 
lomap[val].push_back(ilii); // this relies on map being stable and copy of value.. 
++li; }
// as long as we are going through making the keys again, wth just make
// the ordering based on the already existin lists too LOL
const FileTy & fi=(*ilii).second;
FI fii=fi.begin();
while ( fii!=fi.end())
{
// terribly specialized
const IntTy seq=::atoi((*fii)[0].c_str());
//MM_MSG(" seq is "<<seq<<" and sv is "<<sv);
order[io]=seq;
++io;
++fii; }



 
++ilii; }
std::vector<IntTy> site_values;
LMI lmi=lomap.begin();
while ( lmi!=lomap.end()) { site_values.push_back((*lmi).first); ++lmi; } 
std::sort(site_values.begin(),site_values.end());
MM_MSG("")

const IdxTy sites=lomap.size();
typedef char CodeTy;
CodeTy codes[seqs][sites];
::memset(codes,'0' ,seqs*sites*sizeof(CodeTy));
//LII lii=idx.begin(); while ( lii!=idx.end()) { const LineTy & key =(*lii).first; ++lii; } 
std::vector<IntTy>::const_iterator  svi=site_values.begin();
//MM_MSG("")
IdxTy siteno=0;

while ( svi!=site_values.end())
{
MM_MSG(" site "<<siteno<<" value "<<(*svi))
++svi;++siteno;
}
svi=site_values.begin();
siteno=0;
// this is dumb, the allowed keyw are in teh map originally do
// not need to do exhasutive search again
while ( svi!=site_values.end())
{
//MM_MSG(" processing site "<<(*svi))
const IdxTy sv=(*svi);
IVII ivii=lomap[sv].begin();
IVII ivie=lomap[sv].end();
while (ivii!=ivie)
{ // the first part of this is the key, which contains sv at some point and thenthesecond part
// is a list of sequences that need serial numbers. 
const FileTy & fi = (*(*ivii)).second;
// this is specialized to a single numeric field as the key will fail else...
FI fii=fi.begin();
while ( fii!=fi.end())
{
// terribly specialized
const IntTy seq=::atoi((*fii)[0].c_str());
//MM_MSG(" seq is "<<seq<<" and sv is "<<sv);

codes[seq][siteno]='1';
++fii; }

++ivii; }
//MM_MSG("")

++siteno;
++svi; }

// printing code, the marix needs its own object
StrTy sep=" ";
OsTy & os=std::cout;
if ( (seqs*sites)<1) return; 
CodeTy * plast =codes[order[0]];
std::vector<IdxTy> sames;
for ( IdxTy i=0; i<seqs; ++i)
{
const IdxTy oi=order[i];
if ( oi==~0U ){  MM_MSG(" invalid at "<<i<<" with "<<((int)oi)) continue; } 
CodeTy * p =codes[oi];
IdxTy diffs=0,ones=0;;
for ( IdxTy j=0; j<sites; ++j){ if ( p[j]!='0') ++ones; if ( plast[j]!=p[j]) ++diffs; } 
const bool print_this=(diffs>0)||(i==0);
if ( print_this) {
os<<sep<<sames.size()<<CRLF;
if ( sames.size()!=0)
{
os<<MM_MARK<<" same as ";
for ( IdxTy ss=0; ss<sames.size(); ++ss) os<<sep<<sames[ss];
sames.clear();
os<<CRLF;

 } //sames

os<<i<<sep<<order[i]<<sep<<diffs<<sep<<ones<<"\t";
for ( IdxTy j=0; j<sites; ++j)
{ 

if ( p==plast) os<<p[j]; 
else {
const IdxTy pcode=((p[j]=='0')?0:1)+((plast[j]=='0')?0:2);
// if (plast[j]==p[j]) os<<" "; else os<<p[j];
switch ( pcode) {
case 0: os<<'0'; break;
case 1: os<<'+'; break;
case 2: os<<'-'; break;
case 3: os<<'1'; break;
}; 


 } 
}
// this is no delayed until sames sizeis abaiable but this is dumb
//os<<CRLF;

} else sames.push_back(order[i]);
plast=p;
}

if ( sames.size()!=0)
{
os<<sep<<sames.size()<<CRLF;
os<<MM_MARK<<" same as ";
for ( IdxTy ss=0; ss<sames.size(); ++ss) os<<sep<<sames[ss];
sames.clear();
os<<CRLF;

 } //sames
else os<<sep<<0<<CRLF;

//MM_MSG("")


// this is memory problem for large number of sites and sequences but usualyh not a big deal
// for hundred sites and thousand sequences.  ( about 1 Meg or so )

// these need to be numerical, and sorted.



}


// given a list of key fields, make a key and value and then dump by value
static void index_by_key(StrIndexMapTy & d, const FileTy & s, const IdxTy key, const StrTy & sep)
{
const IdxTy sz=s.size();
const StrTy def=StrTy("");
for ( IdxTy i=0; i<sz; ++i)
{
const LineTy & line=s[i];
const IdxTy len=line.size();
const StrTy & k=(key<len)?line[key]:def;
d[k]+=FockTy::rol(line,key,sep);

} //i 

}

static void index_by_key(LineIndexMapTy & d, const FileTy & s, const IdxTy key, const StrTy & sep)
{
const IdxTy sz=s.size();
const StrTy def=StrTy("");
for ( IdxTy i=0; i<sz; ++i)
{
const LineTy & line=s[i];
const IdxTy len=line.size();
const StrTy & k=(key<len)?line[key]:def;
LineTy lkey;
lkey.push_back(k);
d[lkey]=(FockTy::rolv(line,key,sep));

} //i 

}





template < class Ty> static bool ok(const Ty & is) { return is.good()&&!is.eof(); } 

// this is specialized to the rule hits,
/*
slightly processed into sequence, start, stop, hit, "rule", rule number , comments

 clustals -lihitssn oxx/vp40/vp40s 13   | more
0 0 0 M >rule 14    hydrophobic, http://en.wikipedia.org/wiki/Amino_acid
0 1 1 R >rule 15    basic , http://en.wikipedia.org/wiki/Amino_acid
0 2 2 R >rule 15    basic , http://en.wikipedia.org/wiki/Amino_acid
0 3 3 V >rule 14    hydrophobic, http://en.wikipedia.org/wiki/Amino_acid
*/
// make a bunch of common stats let user sort with grep... 
template <class Ty> static void make(LineTy & d, const LineTy & s, Ty & f, const IdxTy sz)
{
for (IdxTy i=0; i<sz; ++i) d.push_back(s[f[i]]); 
}
template <class Ty> static IdxTy summary(OsTy& os, const FileTy & f,const Ty &param= Ty())
{
const StrTy sep=" ";
// these could be overridden
const IdxTy col_seq=f.column(StrTy("sequence")); 
const IdxTy col_start=f.column(StrTy("start")); 
//const IdxTy col_hit=f.column(StrTy("hit")); 
const IdxTy col_rule=f.column(StrTy("rule")); 
const IdxTy col_desc=f.column(StrTy("description")); 
//const IdxTy col_start=f.column(StrTy("rule_literal")); 
// show total hits by sequence and rule in format seq rule rule rule etc.
typedef std::vector<IdxTy> ListSeq;
typedef std::map<IdxTy,IdxTy> KeyCount;
typedef std::map<IdxTy,ListSeq> KeyList;
typedef std::map<IdxTy,StrTy> DescMap;
typedef std::map<IdxTy, KeyCount> Count2Map;
typedef std::map<IdxTy, KeyList> Count3Map;

//typedef std::map<LineTy, IdxTy> CountMap;
//CountMap cm;
Count2Map cm2;
IdxTy maxseq=0;
std::vector<IdxTy> fields;
DescMap dm;
fields.push_back(col_seq);
fields.push_back(col_rule);
for( auto i=f.begin(); i!=f.end(); ++i)
{
IdxTy seq=::atoi(((*i)[col_seq]).c_str());
if (seq>maxseq) maxseq=seq;
}
++maxseq; // this is a bit pathological if there are none... 
for( auto i=f.begin(); i!=f.end(); ++i)
{
//LineTy key;
//make(key,*i,fields,fields.size());
//++cm[key];
const IdxTy rule=::atoi((*i)[col_rule].c_str());
//++cm2[::atoi((*i)[col_rule].c_str())][::atoi((*i)[col_seq].c_str())];
++cm2[rule][::atoi((*i)[col_seq].c_str())];
dm[rule]=(*i)[col_desc]; // redundant, need this elsewhere... 
}

summarize(os,sep,"summary", cm2,maxseq,dm);
// as above but for each hit type at each site
fields.clear();
fields.push_back(col_rule);
fields.push_back(col_start);
cm2.clear();
// this is really backwards, want a lits of site occupancy for each seq. 
Count3Map cm3;
// this can not use cm2, needs a vector of values. 
for( auto i=f.begin(); i!=f.end(); ++i)
{
const IdxTy rule=::atoi((*i)[col_rule].c_str());
const IdxTy pos=::atoi((*i)[col_start].c_str());
const IdxTy seq=::atoi((*i)[col_seq].c_str());
// in theory anyway the vector should now be sorted in sequence order...
auto & vec=cm3[rule][pos];
const IdxTy vsz=vec.size();
if ( vsz>0) if (seq<=vec[vsz-1]) MM_MSG(" hits need to be sorted "
<< rule<<" "<<" seqis "<< seq<<" after "<<vec[vsz-1]<<" at "<<pos)
cm3[rule][pos].push_back(seq);

}

site_survey(os,sep,"site-survey", cm3,maxseq);

return 0; 
}


template <class Ty> static void site_survey(OsTy &os, const StrTy & sep, const StrTy& label, const Ty &cm3, const IdxTy maxseq)
{


}

//  for each rule number, show number of hits in each sequence along
// with min max etx
template <class Ty,class Tx > static void summarize(OsTy &os, const StrTy & sep, const StrTy& label, const Ty &cm2, const IdxTy maxseq,  Tx & dm)
{
const bool header=true;
const bool print_desc=true;
const bool skip_identical=false;
const bool skip_different=false;
const bool short_identical=false;
if ( header) { 
os<<"# label"<<sep<<"rule"<<sep<<"min"<<sep<<"max"<<sep<<"diffs"<<
sep<<"seqmin"<<sep<<"seqmax";
for (IdxTy k=0; k<maxseq; ++k) os<<sep<<" seq"<<k;
os<<CRLF;
}
for( auto i=cm2.begin(); i!=cm2.end(); ++i)
{
const IdxTy rule=(*i).first;
const StrTy & desc=dm[rule];
//const KeyCount & kc=(*i).second;
const auto & kc=(*i).second;
// there is no assurance the maps traverse in the same order...
IdxTy c[maxseq];
IdxTy crsmin=~0U;
IdxTy crsmax=0;
IdxTy seqmax=0;
IdxTy seqmin=0;

for (IdxTy k=0; k<maxseq; ++k) c[k]=0;
for( auto j=kc.begin(); j!=kc.end(); ++j)
{
IdxTy crs=(*j).second; 
const IdxTy seq=(*j).first; 
c[seq]=crs; 
if (( crs<crsmin)||(crsmin==~0U)) { seqmin=seq; crsmin=crs;} 
if ( crs>crsmax) {seqmax=seq; crsmax=crs;}
}
// probably mean this?
for (IdxTy k=0; k<maxseq; ++k) if(c[k]==0){ seqmin=k; crsmin=0; break; } 
const bool differ=(crsmax-crsmin);
const bool do_print=(!skip_identical&&!differ) ||( !skip_different&&differ);
const bool do_full=do_print&&!(!differ&&short_identical);
if (do_print)
{
if ( do_full)
os<<label<<sep<<rule<<sep<<crsmin<<sep<<crsmax<<sep<<(crsmax-crsmin)
<<sep<<seqmin<<sep<<seqmax;
else os<<label<<sep<<rule<<sep<<crsmin;

} // do print 
for (IdxTy k=0; k<maxseq; ++k) os<<sep<<c[k];
if ( print_desc&&do_full) os<<sep<<">"<<sep<< desc;
if (do_print) os<<CRLF;

}



}
// load a translation table to use with main rules hits 
// for sequence name matching.  This should work gracefully with
// tokenized or compressed csv 
static IdxTy load_index(FileTy & f, const StrTy & nm)
{


return 0;
}

// check that the file is likely to be a rules hit file
static IdxTy validate(FileTy & f, const IdxTy flags)
{
IdxTy errs=0;
IdxTy line=0;
// make sure that each rule literal is "rule" and that all fields
// are there, numerics are valid etc.
for ( auto i=f.begin(); i!=f.end(); ++i)
{
const LineTy & l=(*i);
if ( l.size()<6)  { ++errs;MM_MSG("bad line "<<line<<" "<<FockTy::dump(l))  } 
else if ( l[4]!=">rule") 
{ ++errs; 
MM_MSG("line "<<line<<" instead of >rule have "<<l[4] <<" fow line "<< FockTy::dump(l)); } 
++line;

}



return errs;

}

// right now these are passed around as csv files,
// should have struct


static IdxTy ai(const StrTy & x) { return ::atoi(x.c_str()); } 
template < class Ty> static void load_struct(Ty & d, const FileTy & f)
{
const IdxTy col_seq=f.column(StrTy("sequence")); 
const IdxTy col_start=f.column(StrTy("start")); 
const IdxTy col_end=f.column(StrTy("end")); 
const IdxTy col_hit=f.column(StrTy("hit")); 
const IdxTy col_rule=f.column(StrTy("rule")); 
const IdxTy col_desc=f.column(StrTy("description")); 
const IdxTy col_lit=f.column(StrTy("rule_literal")); 
const bool load_redundant_crap=false;
const bool check_rule_literal=false;


for ( auto i=f.begin(); i!=f.end(); ++i)
{
const LineTy & l=(*i);
if ( check_rule_literal) if ( l[col_lit]!=">rule") {}
HitClass hc;
hc.set(ai(l[col_seq]),ai(l[col_rule]),ai(l[col_start]),ai(l[col_end]));
if ( load_redundant_crap) {hc.set(l[col_hit],l[col_desc]);  }
else {hc.set(l[col_hit],"");  }


d.push_back(hc);


} //i


} // load_struct


// numeric conversion should be part of csv stuff 
static IdxTy load(FileTy & f, const StrTy & nm)
{
const IdxTy fields=7;
FockTy::load_ragged(f,nm,fields,0);
IdxTy b=0;
f.set_column(StrTy("sequence"),b); ++b;
f.set_column(StrTy("start"),b); ++b;
f.set_column(StrTy("end"),b); ++b;
f.set_column(StrTy("hit"),b); ++b;
f.set_column(StrTy("rule_literal"),b); ++b;
f.set_column(StrTy("rule"),b); ++b;
f.set_column(StrTy("description"),b); ++b;

f.make_index();




return 0; 
/*
const IdxTy max=1<<16;
const bool dodel=(nm!=StrTy("-"));
 IsTy * is = dodel?(new std::ifstream(nm.c_str())):(&std::cin);
ChTy * line = new ChTy[max+2];
while (ok(*is))
{
is->getline(line,max);
if ( !ok(*is)) break;
SsTy ss((line));
LineTy x;
while ( ok(ss)) 
{
StrTy y; ss>>y;
if ( ok(ss)) x.push_back(y);
}
f.push_back(x);
}
// mem leak on throw? 
delete [] line;
if ( dodel) delete is;
return 0; 
*/

}
 


}; // mjm_ ?? 

}; // mjm_sequence_hist ns

#endif


