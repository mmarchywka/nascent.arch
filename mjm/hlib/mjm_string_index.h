#ifndef MJM_STRING_INDEX_H__
#define MJM_STRING_INDEX_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
#include "mjm_generic_iterators.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_block_matrix.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_collections.h"
#include "mjm_svg_writer.h"
#include "mjm_color_table.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>

// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
//#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
//#define MM_DMELF(fm,code,x) f.DMel(code, f.m_ss<<MM_STR_LOC<<" "<<x); 

namespace string_index_traits
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
//typedef mjm_sparse_matrix<D> MySparse;
}; // 
}; // biom_hdf5_traits

template<class Tscalar, class Tval, class Tvec=std::vector<Tval> > 
class mjm_1_to_v 
{
typedef mjm_1_to_v Myt;
protected:
typedef string_index_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;
//typedef std::map<char, IdxTy > Wcmap;
public:
typedef Tscalar key_t;
typedef Tval value_t; 
//typedef std::vector<value_t> vector_t;
typedef Tvec vector_t;
typedef std::map<key_t,vector_t> object_t;
typename object_t::iterator begin() { return m_obj.begin(); }
typename object_t::iterator end() { return m_obj.end(); }
typename object_t::const_iterator begin() const { return m_obj.begin(); }
typename object_t::const_iterator end() const { return m_obj.end(); }

void add(const key_t & k, const value_t & v) { m_obj[k].push_back(v); }
const vector_t & operator[]( const key_t& k) const
{
auto ii=m_obj.find(k); if (ii!=m_obj.end()) { return (*ii).second; }
return m_null; 
} 
void clear() { m_obj.clear(); }
private:
object_t m_obj;
vector_t m_null;
};  // mjm_1_to_v

//template <class Tobj, class Traits >
class mjm_string_index 
{
typedef mjm_string_index Myt;
protected:
typedef string_index_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

typedef std::map<char, IdxTy > Wcmap;

class int_vector : public std::vector<int>
{
public:

int sum() const {int s=0;  MM_LOOP(ii,(*this)) { s+=(*ii); } return s; } 
//int maxkey() const {int sum=0;  MM_LOOP(ii,(*this)) { s+=(*ii); } return s; } 


}; // int_Vector

typedef  mjm_1_to_v <int,int, int_vector> One2N;

public:
mjm_string_index(): m_dmel(0) {dna_map(); }
virtual ~mjm_string_index() {}

IdxTy comp_class( const char c1, const char c2, Wcmap & m )
{
const IdxTy& t1=m[c1];
const IdxTy& t2=m[c2];
if ((t1!=1) ||(t2!=1)) return 2;
if (c1==c2) return 0;
return 1;
}

void base_map()
{
Wcmap&  m= m_wcmap;
m.clear();
m['A']=1;
m['C']=1;
m['G']=1;
m['T']=1;


}
class Summary
{
public:
Summary(const IdxTy e, const IdxTy d, const IdxTy n ): eq(e),diff(d),ncmp(n) {}
IdxTy eq,diff,ncmp;

};
typedef Summary summary_type;



class align_type {
typedef align_type Myt;
public:
align_type ( const IdxTy ii, const IdxTy jj, const StrTy & ss): i(ii),j(jj),s(ss),del(int(i)-int(j)),len(s.length()) {}
const IdxTy &  loc(const bool one_not_two) const { return one_not_two? i:j; } 
const int  delta(const bool one_not_two) const { return one_not_two?(del):(-del);  } 
const IdxTy & size() const { return len; } 
 IdxTy olaps(const Myt & that) const
{
bool nolap1=(i>(that.i+that.len)||((i+len)<that.i));
bool nolap2=(j>(that.j+that.len)||((j+len)<that.j));
return (nolap1?0:1) + (nolap2?0:2) ; }
bool greater_than(const Myt & that ) const
{ return ((i>(that.i+that.len))&&(j>(that.j+that.len))); }
bool less_than(const Myt & that ) const
{ return that.greater_than(*this); }
bool  crosses(const Myt & that ) const
{ return !(greater_than(that)||less_than(that));  }

class start_sort
{
public:
template <class Tx, class Ty> 
bool operator()(Tx & x, Ty  & y ) const { return x.i<y.i; }
};

// fudding sort fudd 
//const IdxTy i,j;
//const StrTy s;
//const int del;
//const IdxTy len;
 IdxTy i,j;
 StrTy s;
 int del;
 IdxTy len;

}; // align_type

class space_type {
public:
space_type ( const IdxTy ii, const IdxTy jj): x(ii),len(jj) {}
const IdxTy x,len;
//const StrTy s;
};
// TODO fix conventions 
typedef std::vector<align_type> HiVec;
typedef std::map<int,IdxTy> offset_hist_type;



summary_type comp_base_sequences(const char * s1, const char * s2 , const IdxTy & flags)
{
const char * p1=s1;
const char * p2=s2;
IdxTy eq=0;
IdxTy diff=0;
IdxTy ncmp=0;
while ((*p2!=0)&&(*p1!=0))
{
IdxTy cmp=comp_class(*p1,*p2,m_wcmap);
switch (cmp)
{
case 0: { ++eq;   break; }
case 1: { ++diff; break; }
case 2: { ++ncmp; break; }

default : {  MM_ERR(" bad comp "<<MMPR4(cmp,p1,p2,s1)<<MMPR(s2)) } 
}; // case

++p1;
++p2;
}

return summary_type(eq,diff,ncmp); 
}


/*

Setting 100 perc_identity and word_size equal length fails on non-AGCT chars.
Identity needs to eliminate them, (bases-pseudo)/bases and word size does
not seem to help speed or anything so just set to a default. 
COMMAND blastn -evalue 1e-3 -ungapped -word_size 4 -perc_identity 95 -num_alignments 3 -remote -db nr

https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
		A  adenosine          C  cytidine             G  guanine
		T  thymidine          N  A/G/C/T (any)        U  uridine 
		K  G/T (keto)         S  G/C (strong)         Y  T/C (pyrimidine) 
		M  A/C (amino)        W  A/T (weak)           R  G/A (purine)        
		B  G/T/C              D  G/A/T                H  A/C/T      
		V  G/C/A              -  gap of indeterminate length
*/


void dna_map()
{
Wcmap&  m= m_wcmap;
Wcmap&  p= m_wcsubs;
m.clear();
p.clear();
m['K']=2;
m['M']=2;
m['B']=3;
m['V']=3;
m['N']=4;
m['S']=2;
m['W']=2;
m['D']=3;
m['Y']=2;
m['R']=2;
m['H']=3;

const char * chars= "GTACGTCGCAAGCTGCATGATTCGAACT";
::memcpy(m_subs,chars,strlen(chars));

p['K']=0;
p['M']=p['K']+m['K'];
p['B']=p['M']+m['M'];
p['V']=p['B']+m['B'];
p['N']=p['V']+m['V'];
p['S']=p['N']+m['N'];
p['W']=p['S']+m['W'];
p['D']=p['W']+m['W'];
p['V']=p['D']+m['D'];
p['R']=p['V']+m['V'];
p['H']=p['R']+m['R'];

}
template < class Td> 
void enumerate_wilds( Td & v, const char * s, const IdxTy cardmax=0)
{
typedef mjm_generic_itor<1> Itor;
typedef Itor::location_type location_type;
location_type loc;
std::vector<IdxTy> cloc;
Wcmap&  wcmap= m_wcmap;
Wcmap&  wcsubs= m_wcsubs;
const IdxTy slen=strlen(s);
char c[slen+1];
for (IdxTy i=0; i<slen; ++i)
{
const char & sc=s[i];
if (wcmap.find(sc)!=wcmap.end()) {loc.push_back(wcmap[sc]);  cloc.push_back(i); }
c[i]=sc;

}
c[slen]=0;

IdxTy sz=loc.size();
if (sz==0) { v.push_back(StrTy(c)); return; } 
IdxTy card=1;
// TODO FIXME this can easily oflow doh doh doh 
MM_LOOP(ll,loc) { card*=(*ll);if (cardmax!=0) if (card>cardmax) return;  } 
MM_ERR(" expanding "<<MMPR2(c,card))
Itor ii(loc);

while (!ii.done())
{
IdxTy idx=0; 
MM_LOOP(j,cloc)
{
c[*j]=m_subs[wcsubs[idx]+ii[idx]];

++idx;
}
//MM_ERR(MMPR(c))
v.push_back(StrTy(c));
ii.inc(); 
} // ii 

} 
//////////////////////////////////////////////////
// count number of chars thought to be offset by that amount 
void offset_histogram( offset_hist_type & oh, const HiVec & v)
{
MM_LOOP(ii,v) { oh[(*ii).del]+=(*ii).len; } // v

} // offset_histogram

class discovered_offsets
{
typedef discovered_offsets Myt;
typedef One2N Pt;
typedef std::vector<Pt> Pvec;
public:
typedef Pt offsets_type;
typedef Pvec offsets_v_type;
discovered_offsets(const IdxTy sz):m_sz(sz),m_p(m_sz)  {}
void locations_and_lengths(const HiVec & v, const IdxTy flags)
{
//const bool one_not_two=((flags&1)!=0);
MM_LOOP(ii,v)
{
//const auto & vi=(*ii);
//const IdxTy loc=vi.loc(one_not_two); // ? vi.i:vi.j;
//const int del=vi.delta(one_not_two); // ?(vi.del):(-vi.del); 
//const IdxTy len=vi.size();
//const IdxTy loce=loc+len;


}

}
offsets_type::vector_t  best_offset(const offsets_type & ot, const IdxTy flags)
{
 offsets_type freq_inverse;
IdxTy i=0;
MM_LOOP(ii,ot)
{
// this maps the score into the list of possible offsets
freq_inverse.add(-(*ii).second.sum(),(*ii).first);
MM_ERR(" freq_inverse score and offset "<<MMPR3(i,-(*ii).second.sum(),(*ii).first))
++i;
} // ot  
auto ii=freq_inverse.begin();
//MM_LOOP(jj,freq_inverse) { MM_ERR(" freq_inverse "<<MMPR2((*jj).first,(*jj).second.)) } 
if (ii!=freq_inverse.end()) return (*ii).second;

return offsets_type::vector_t(); 
} // best_offset


// for each location in the target sequence, produce a 
// figure ot merit for alignment preference at that location.
// This is a map from offset del  to a vector of strengths
void analyze(const HiVec & v, const IdxTy flags)
{
// this is sequence 1 if set, align_t.i points to us  
const bool one_not_two=((flags&1)!=0);
MM_LOOP(ii,v)
{
const auto & vi=(*ii);
const IdxTy loc=vi.loc(one_not_two); // ? vi.i:vi.j;
const int del=vi.delta(one_not_two); // ?(vi.del):(-vi.del); 
const IdxTy len=vi.size();
const IdxTy loce=loc+len;
const IdxTy inc=1*len;
MM_ERR(" hit mark "<<MMPR4(loc,del,vi.i,vi.j)<<MMPR3(loc,vi.s,loce))
for (IdxTy i=loc; i<loce; ++i)
{
// this is trying to create a historgram of offsets at a given
// point but the utlimate goal is a freuqncy key into offsets 
if ( i>=m_sz) 
{
MM_ERR(" dnager will robinson oob "<<MMPR4(m_sz,i,loce,loc)<<MMPR4(len,vi.i,vi.j,del))
}
else m_p[i].add(del,inc);

} //i 

} // v 


} // analyze 
const offsets_type&  offsets_at(const IdxTy & i ) { return m_p[i]; } 
private:
const IdxTy m_sz;
Pvec m_p;
}; // discovered_offsets

/* 
Two extreme cases- unque strings occur once in each of the pair. Alt, small frequent
strings or single chars create signal to correlate. Right now it looks like just
eliminating the impossible will work- anitcrossing for example. 
*/
static const IdxTy & bad()  { static const IdxTy i=~0U; return i; } 
void anti_crossing(HiVec & d, const HiVec & s)
{
std::map<IdxTy,IdxTy> m;
IdxTy i=0;
IdxTy mac=0;
IdxTy idxm=bad();
MM_LOOP(ii,s)
{
IdxTy j=0; 
const auto & vi=(*ii);
MM_LOOP(jj,s)
{
if (ii==jj) continue; 
const auto & vj=(*jj);
if ( vi.crosses(vj)) {++m[i]; ++m[j]; 
if (m[i]>mac) {mac=m[i]; idxm=i; }
if (m[j]>mac) {mac=m[j]; idxm=j; }

 }
++j;
} // jj 
++i;
} // ii 
MM_SZ_LOOP(k,s,sz)
{
if (k!=idxm) d.push_back(s[k]);
else {MM_ERR(" anticorss reject "<<MMPR4(k,mac,s[k].i,s[k].j)) }
}

}


void unique(HiVec & d, const HiVec & s, const HiVec & s2)
{
// in the case of one strig inside another, take the logest lol. 
HiVec t;
const bool check_unique=true;
const bool debug_unique=!true;
if (check_unique) { 
std::map<StrTy, IdxTy> m;
MM_LOOP(ii,s) { ++m[(*ii).s]; }
MM_LOOP(ii,s2) { ++m[(*ii).s]; }
//MM_LOOP(ii,s) {MM_ERR(MMPR2((*ii).s,m[(*ii).s])) if (m[(*ii).s]==2) d.push_back((*ii)); }  
MM_LOOP(ii,s) {

if (debug_unique ) { MM_ERR(MMPR2((*ii).s,m[(*ii).s]))  } 

if (m[(*ii).s]==2) t.push_back((*ii)); }  
}
if (true)
{
std::map<IdxTy, IdxTy> mloc;
const HiVec & _s=(check_unique?t:s);
MM_LOOP(ii,_s)
{
bool dokeep=true;
for (auto jj=(_s.begin()); jj<_s.end(); ++jj)
{
if (ii==jj) continue;
if ( (*jj).olaps((*ii))!=0) if ((*ii).len<(*jj).len) {

if (debug_unique ) { MM_ERR(" reject "<<MMPR((*ii).del)<<MMPR4((*ii).i,(*ii).j,(*ii).len, (*ii).s)<<MMPR((*jj).del)<<MMPR4((*jj).i,(*jj).j,(*jj).len, (*jj).s)) } 

 dokeep =false; break;  } 

} // jj
if (dokeep) {
d.push_back((*ii)); 
if (debug_unique) { MM_ERR(" keep "<<MMPR((*ii).del)<<MMPR4((*ii).i,(*ii).j,(*ii).len, (*ii).s)) } 
}
} // ii 

} // true 



} // unique 
/*
The strings need to be ONLY between the two being aligned lol 
mjm_string_index.h387  (*ii).s=AAGTCTGAT m[(*ii).s]=2
mjm_string_index.h387  (*ii).s=AAGTCTGATGTGAAAATGC m[(*ii).s]=2
mjm_string_index.h387  (*ii).s=AAGTCTGAT m[(*ii).s]=2
mjm_string_index.h387  (*ii).s=AAGTCTGATGTGAAAATGC m[(*ii).s]=2
mjm_string_index.h362  hit mark  loc=358 del=-198 vi.i=358 vi.j=556 loc=358 vi.s=AAGTCTGAT loce=367
mjm_string_index.h362  hit mark  loc=550 del=-6 vi.i=550 vi.j=556 loc=550 vi.s=AAGTCTGATGTGAAAATGC loce=569
mjm_string_index.h362  hit mark  loc=556 del=198 vi.i=358 vi.j=556 loc=556 vi.s=AAGTCTGAT loce=565
mjm_string_index.h362  hit mark  loc=556 del=6 vi.i=550 vi.j=556 loc=556 vi.s=AAGTCTGATGTGAAAATGC loce=575



*/

// try go minimize errors 
void maxmatches (  char * c1, char * c2, const char * p1, const char * p2 
,int&  ps1, int&  ps2, int&  pd1, int&  pd2
, const int & loc1, const int &  loc2, const int & len, const char space)
{
const bool debug_dist=true;
// this is the number of SOURCE characters to use up 
int rsz=loc1-ps1;
// pd2f>=pd1f
int pd1f=pd1+loc1-ps1;
int pd2f=pd2+loc2-ps2;
// start slot in destination 
int start=pd1;
// the last dest pointer needs to match 2 
int end= pd2+loc2-ps2;
//int end=pd2f-pd1f+ pd2+loc2-ps2;
// number of spaces to write
int spaces=pd2f-pd1f; // end-start-rsz;
// number of chars to write 
int chars=rsz;
const int ssz=chars;
if (debug_dist) { MM_ERR(MMPR4(start,end,spaces,chars)) } 
mjm_block_matrix<IdxTy> cumoff(spaces+1,chars);

for(int osett=0; osett<=spaces; ++osett)
{
Ss ss;
for(int i=0; i<chars; ++i)
{
const bool agree=(*(p2+ps2+i+osett)==*(p1+ps1+i));
cumoff(osett,i)=(agree?1:0);
ss<<cumoff(osett,i);
//if (i!=0) cumoff(osett,chars-1-i)+=cumoff(osett,chars-i);
}
MM_ERR(MMPR2(osett,ss.str()))
}

for(int osett=0; osett<=spaces; ++osett)
{
Ss ss;
if ( chars!=0) { ss<<cumoff(osett,0)<<" "; } 
for(int i=1; i<chars; ++i)
{
//cumoff(osett,chars-1-i)+=cumoff(osett,chars-i);
cumoff(osett,i)+=cumoff(osett,i-1);
//ss<<cumoff(osett,chars-1-i)<<" ";
ss<<cumoff(osett,i)<<" ";
}
MM_ERR(MMPR2(osett,ss.str()))
}



int  osetp=0;
for(int i=0; i<(chars); ++i)
{
while (true) { 
bool insert_space=(osetp<spaces);
if (insert_space) insert_space=(cumoff(osetp+1,i)>cumoff(osetp,i));
if (!insert_space) { *(c1+pd1+i+osetp) =*(p1+ps1+i); break; } 
else {
if (debug_dist) { MM_ERR(" space "<<MMPR4(i,osetp,cumoff(osetp+1,i),cumoff(osetp,i))) } 
*(c1+pd1+i+osetp) =space; ++osetp; }  
}
} // i 
while (osetp<spaces)  {*(c1+pd1+chars+osetp) =space; ++osetp; }  


// now copy the short one 
ps1+=loc1-ps1; // yset this should be the loc1 lol 
ps1=loc1; // this is right but above should be right too 
pd1+=ssz+pd2f-pd1f; // pd2-loc1+ps1;
rsz=loc2-ps2;
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2=loc2; // ps2+=rsz;
if (debug_dist)
{
char ca[rsz+1];
::memcpy(ca,c1+start,rsz);
ca[rsz]=0;
MM_ERR("1 "<<ca);
::memcpy(ca,c2+pd2-rsz,rsz);
MM_ERR("2 "<<ca);
}


} // maxmatches


// the ends are most likely spots of difference just put the gaps there 
void distribute (  char * c1, char * c2, const char * p1, const char * p2 
,int&  ps1, int&  ps2, int&  pd1, int&  pd2
, const int & loc1, const int &  loc2, const int & len, const char space)
{
const bool debug_dist=true;
// this is the number of SOURCE characters to use up 
int rsz=loc1-ps1;

// this is the dest pointer final value after writinge up to loc1 ( 
int pd1f=pd1+loc1-ps1;
int pd2f=pd2+loc2-ps2;
// start slot in destination 
int start=pd1;
// the last dest pointer needs to match 2 
int end= pd2+loc2-ps2;
//int end=pd2f-pd1f+ pd2+loc2-ps2;
// number of spaces to write
int spaces=pd2f-pd1f; // end-start-rsz;
// number of chars to write 
int chars=rsz;
if (debug_dist) { MM_ERR(MMPR4(start,end,spaces,chars)) } 
int cko=0; 
int cke=0; 
int idxs=0;
int idxs2=0;
int idxd=0;
int placed=0;
int iter=0;
const int ssz=chars;
while ((chars>0)&&(spaces>0))  // (chars>spaces)
{
const bool left= ((iter&1)==0);
if (left)
{
idxd=start+(iter>>1)-pd1;
idxs=idxd-cke; // ps1+cke-ps1;
idxs2=idxd-cke; // ps1+cke-ps1;
} // left
else
{
idxd=end-(iter>>1)-pd1-1;
idxs=idxd+cko; // ps1+ssz-cko-ps1;
idxs2=idxd+cko-(+pd2f-pd1f); // ps1+cke-ps1;

} // right 
const bool achar=(*(p2+ps2+idxs2)==*(p1+ps1+idxs));
if (debug_dist){ 
IdxTy b2=ps2+idxd;
IdxTy b1=ps1+idxs; 
MM_ERR(MMPR2(b1,b2)<<MMPR4( achar,p2[b2],p1[b1],idxs)<<MMPR4(idxd,chars,spaces,placed))

 } 
if (achar) {
// these are the number of chars inserted to offset 
 *(c1+pd1+idxd) =*(p1+ps1+idxs2); --chars; ++placed; }
else { 
if (left) {++cke;} else {++cko; } 

*(c1+pd1+idxd) =space; --spaces;}
++iter;
} // chars>spaces
// the remnants need to be placed at only available locations hopefully near middle
idxd=start+((iter+1)>>1)-pd1;
idxs=idxd+ps1-cke-ps1;
while (spaces>0) {*(c1+pd1+idxd)=space; ++idxd;     --spaces; }
while (chars>0) {  *(c1+pd1+idxd)=*(p1+ps1+idxs); ++idxd; ++idxs;    --chars; }
// 
ps1+=loc1-ps1; // yset this should be the loc1 lol 
ps1=loc1; // this is right but above should be right too 
pd1+=ssz+pd2f-pd1f; // pd2-loc1+ps1;
rsz=loc2-ps2;
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2=loc2; // ps2+=rsz;
if (debug_dist)
{
char ca[rsz+1];
::memcpy(ca,c1+start,rsz);
ca[rsz]=0;
MM_ERR("1 "<<ca);
::memcpy(ca,c2+pd2-rsz,rsz);
MM_ERR("2 "<<ca);
}


} // distribute


void make_alignment( StrTy & d1, StrTy&  d2, const StrTy s1, const StrTy & s2, const HiVec &  _vf, const HiVec & _vr) 
{
const bool debug_align=!false;
IdxTy sz1=s1.length();
IdxTy sz2=s2.length();
IdxTy sumzed=sz1+sz2+1;
typedef std::vector<space_type> SpVec;
if (debug_align) { 
MM_ERR(" aligning "<<MMPR3(_vf.size(),sz1,s1))
MM_ERR(" aligning "<<MMPR3(_vr.size(),sz2,s2))
}
SpVec s1s,s2s;
discovered_offsets dd1(sz1),dd2(sz2);
// each loaction has a map of offset to a vector of hint strengths 
// first remove hits that have a string occuring in multiple locations ...
HiVec vf,vr;
unique(vf,_vf,_vr);
unique(vr,_vr,_vf);

while (true) { 
HiVec vf2=vf;
vf.clear(); // doh 
anti_crossing( vf,  vf2);
if ( vf.size()==vf2.size()) break;
} // while (vf.size()<vf2.size());
if (debug_align) { MM_ERR(" done anticrossing "<<MMPR(vf.size())) } 
const bool direct_write=true;
if (direct_write)
{
// these must be in order 
align_type::start_sort so;
// eh these need to be sorted lol
std::sort(vf.begin(),vf.end(), so);
const char space='.';
char c1[10*sumzed], c2[10*sumzed];
::memset(&c1[0],'x',10*sumzed);
::memset(&c2[0],'x',10*sumzed);
const char * p1=s1.c_str();
const char * p2=s2.c_str();

int ps1=0;
int  ps2=0;
int  pd1=0;
int  pd2=0;
int  sz=0;
MM_LOOP(ii,vf)
{
const auto & vi=(*ii);
const int loc1=vi.i; // vi.loc(one_not_two); // ? vi.i:vi.j;
const int loc2=vi.j; // vi.loc(one_not_two); // ? vi.i:vi.j;
const int del=loc1-loc2; // ?(vi.del):(-vi.del); 
const int len=vi.size();
//const int loce=loc+len;
//const int inc=1*len;
if (debug_align) { 
MM_ERR(" hit mark "<<MMPR4(loc1,del,vi.i,vi.j)<<MMPR3(len,loc2,vi.s))
}
// write out whateverr exists between last write and these coordinates
/// note these are not matching and oculd be interspersed on short stting 
#if 0
int rsz=loc1-ps1;
::memcpy(c1+pd1,p1+ps1,rsz); pd1+=rsz; ps1+=rsz;
rsz=loc2-ps2;
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2+=rsz;
#endif

//int rsz=loc1-ps1;
int pd1f=pd1+loc1-ps1;
int pd2f=pd2+loc2-ps2;
// fill in with gaps in the discontiguous sequence
// just string in longer one, add spacer in shorter one
//const bool 
const bool distribute=true;
if (pd1f<pd2f)
{
if (distribute)
{
//this->distribute(c1,c2,p1,p2,ps1,ps2,pd1,pd2,loc1,loc2,len,space);
this->maxmatches(c1,c2,p1,p2,ps1,ps2,pd1,pd2,loc1,loc2,len,space);
} // distribute 
else { 
int rsz=loc1-ps1;
//MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR2(pd2,pd1))
::memcpy(c1+pd1,p1+ps1,rsz); pd1+=rsz; ps1+=rsz;
rsz=loc2-ps2;
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2+=rsz;

// these can be distributed but both ends, except beginning
// are mismatch and in fact both could be spaced alternatively but
// ignore that for now  
const IdxTy sz= pd2-pd1;
::memset(c1+pd1,space,sz); pd1+=sz;
 } // distribute 
}
else 
{ // pd1>=pd2
if (distribute)
{
//this->distribute(c2,c1,p2,p1,ps2,ps1,pd2,pd1,loc2,loc1,len,space);
this->maxmatches(c2,c1,p2,p1,ps2,ps1,pd2,pd1,loc2,loc1,len,space);
} // distribute 
else { 
int rsz=loc1-ps1;
//MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR2(pd2,pd1))
::memcpy(c1+pd1,p1+ps1,rsz); pd1+=rsz; ps1+=rsz;
rsz=loc2-ps2;
//MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR3(pd2,pd1,rsz))
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2+=rsz;

const IdxTy sz= pd1-pd2;
//MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR3(pd2,pd1,sz))
::memset(c2+pd2,space,sz); pd2+=sz;
} // distribute 
}

if (pd1!=pd2)
{
MM_ERR(" outputs should be same dwr "<<MMPR4(pd1,pd2,ps1,ps2))
}
if (debug_align) { MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR2(pd2,pd1)) } 
// write out common string to both 
if (ps1!=loc1)
{
MM_ERR(" ps1!=loc1 "<<MMPR4(ps1,loc1,pd1,ps2)<<MMPR2(loc2,pd2))
}

if (debug_align) { MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR3(pd2,pd1,len)) }
::memcpy(c1+pd1,p1+ps1,len); pd1+=len; ps1+=len;
if (ps2!=loc2)
{
MM_ERR(" ps2!=loc2 "<<MMPR4(ps1,loc1,pd1,ps2)<<MMPR2(loc2,pd2))
}
if (debug_align) { MM_ERR(MMPR4(loc1,ps1,loc2,ps2)<<MMPR3(pd2,pd1,len))}
::memcpy(c2+pd2,p2+ps2,len); pd2+=len; ps2+=len;

//::memcpy(c1+pd1,p1+ps1,rsz); pd1+=rsz; ps1+=rsz;
//sz=(*ii).len;
//::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2+=rsz;
//sz=(*ii).len;

} // ii 
if (debug_align) { MM_ERR(MMPR4(ps1,ps2,sz1,sz2)<<MMPR2(pd2,pd1)) } 

sz=sz1-ps1;
{ ::memcpy(c1+pd1,p1+ps1,sz); pd1+=sz; ps1+=sz; }
sz=sz2-ps2;
{ ::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz; }

c1[pd1]=0;
c2[pd2]=0;
d1=StrTy(c1);
d2=StrTy(c2);
const bool sanity_check=true;
if (sanity_check)
{
verify_align_sane(c1,p1);
verify_align_sane(c2,p2);


} // sanity_check


return;
} // direct_write
// more dubious code... 

dd1.analyze(vf, 1);
//dd2.analyze(vf, 0);
IdxTy szmin=(sz1<sz2)?sz1:sz2;
int  coff=0;
// TODO probably can terminate early 
for (int  i=0; IdxTy(i)<(szmin); ++i)
{
if (IdxTy(i)>=sz1) break; // no further offsets needed
// get the map between offsets and strength vectors 
const auto & oa= dd1.offsets_at(i); 
// this is a vector of offsets with the best strength  
const auto & best1= dd1.best_offset(oa,0);
if ((i+coff)>=int(sz2)) break; // no further offsets needed
if (coff+i<0) continue;
// get the similar distribution of hints at the correpsonding current offset 
const auto & ob= dd2.offsets_at(i+coff); 
const auto & best2= dd2.best_offset(ob,0);
std::map<int,int> histo;
Ss ss; 
MM_LOOP(ii,best1) { ss<<" best1 "<<(*ii)<<" ";  ++histo[*ii]; }
MM_LOOP(ii,best2) { ss<<" best2 "<<(*ii)<<" "; ++histo[-(*ii)]; }

int most_best=coff;
int streng=0;
MM_LOOP(ii,histo) { if ((*ii).second>streng) { streng=(*ii).second; most_best=(*ii).first; }}
int diff=most_best-coff;
if (diff==0) continue;
MM_ERR(MMPR(ss.str())<<MMPR4(i,coff,most_best,streng)<<MMPR(diff))
// these are maps of offsets and a vector of values  
// this is BACKWARD looking not FWD 
if (diff>0) { MM_ERR( " SPACE! "<<MMPR2(i,diff))  space_type sp(i,diff); s1s.push_back(sp); } 
else { MM_ERR(" SPACE2 "<<MMPR2(i,-diff)) space_type sp(i+0*coff,-diff); s2s.push_back(sp); } 
coff=most_best;

}



/*
int offset=0;
int sp1=0;
int sp2=0;
MM_LOOP(ii,vf)
{
IdxTy loc1=(*ii).i;
IdxTy loc2=(*ii).j;
MM_ERR(" hit "<<MMPR3(loc1,loc2,(*ii).s))
int offnew=loc2-loc1;
int del=offnew-offset;
if (del>0)
{
if (loc1<=sp1) continue;
sp1=loc1;
s1s.push_back(space_type(sp1,del));
sp1+=del;
}
else if (del<0)
{
if (loc2<=sp2) continue;
sp2=loc2;
s2s.push_back(space_type(sp2,-del));
sp2-=del;

}
offset=offnew;
} // vf 
*/

const char space='.';
char c1[10*sumzed], c2[10*sumzed];
const char * p1=s1.c_str();
const char * p2=s2.c_str();

IdxTy ps1=0;
IdxTy ps2=0;
IdxTy pd1=0;
IdxTy pd2=0;
IdxTy sz=0;
MM_LOOP(ii,s2s)
{
//sz=(*ii).x-ps1;
sz=(*ii).len; // (*ii).x-ps2;
const IdxTy loc=(*ii).x;
const bool anticross_violation=(ps1>loc);
const IdxTy rsz=loc-ps1;
MM_ERR(" space record "<<MMPR4((*ii).x,ps1,pd1,sz)<<MMPR2(rsz,anticross_violation))
::memcpy(c1+pd1,p1+ps1,rsz); pd1+=rsz; ps1+=rsz;
sz=(*ii).len;
::memset(c1+pd1,space,sz); pd1+=sz;
}
bool anticross_violation=(ps1>sz1);
//if (sz1<ps1)
{
MM_ERR(" sizes maybe wrong "<<MMPR4(sz1,ps1,pd1,sz)<<MMPR2(sz2,anticross_violation))
}
sz=sz1-ps1;
//MM_ERR(" fudd "<<MMPR4(pd1,ps1,si

//MM_ERR(" fudd "<<MMPR4(pd1,ps1,sz,sz1))
if (!anticross_violation) 
{::memcpy(c1+pd1,p1+ps1,sz); pd1+=sz; ps1+=sz; }


MM_LOOP(ii,s1s)
{
sz=(*ii).len; // (*ii).x-ps2;
const IdxTy loc=(*ii).x;
const bool anticross_violation=(ps2>loc);
const IdxTy rsz=loc-ps2;
MM_ERR(" space record "<<MMPR4((*ii).x,ps2,pd2,sz)<<MMPR2(rsz,anticross_violation))
::memcpy(c2+pd2,p2+ps2,rsz); pd2+=rsz; ps2+=rsz;
sz=(*ii).len;
::memset(c2+pd2,space,sz); pd2+=sz;
}
anticross_violation=(ps2>sz2);
MM_ERR(" sizes maybe wrong "<<MMPR4(sz2,ps2,pd2,sz)<<MMPR2(sz1,anticross_violation))
sz=sz2-ps2;
//MM_ERR(" fudd "<<MMPR4(pd2,ps2,sz,sz2))
if (!anticross_violation) 
::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz;
//while (true)
//{
//::memset(c2+pd2,space,sz); pd2+=sz;
//::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz;
//break; 
//} // while 
c1[pd1]=0;
c2[pd2]=0;
d1=StrTy(c1);
d2=StrTy(c2);
//MM_ERR(" str "<<MMPR2(d1.length(),d2.length()))
}


IdxTy verify_align_sane( const char * c1, const char * p1) const
{
IdxTy ps=0; IdxTy pd=0; IdxTy cnt=0;
while (c1[pd]!=0)
{
if (c1[pd]!='.') 
{
if (p1[ps]!=c1[pd]) { ++cnt; MM_ERR(" out of place "<<MMPR4(ps,pd,p1[ps],c1[pd])) }
++ps;
}
++pd;
}
return cnt;
}



/////////////////////////////////////////////



protected:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }

void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}


Dmel * m_dmel;
Ss m_ss;
//To m_obj;
Wcmap m_wcmap,m_wcsubs;
char m_subs[1<<8];

};  // mjm class_

class mjm_align_collection
{
typedef mjm_align_collection Myt;
typedef mjm_string_index Msi;
protected:
typedef string_index_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

typedef mjm_ragged_table Ragged;

typedef Msi::align_type align_type;
typedef Msi::HiVec hit_vector;

typedef std::map<IdxTy, hit_vector> HvMap;
typedef std::map<IdxTy, HvMap>  HvijMap;

typedef mjm_color_table Mct;
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }


class seq_layout
{
public:
seq_layout()
	:m_x(0),m_y(0),m_pitch(0),m_serial(~0),m_class(~0) {}

seq_layout(const D& x, const D & y, const D & p, const IdxTy & s,const IdxTy & c)
	:m_x(x),m_y(y),m_pitch(p), m_serial(s),m_class(c) {}

D base_x( const IdxTy i) const { return m_x+i*m_pitch; } 
D base_y( const IdxTy i) const { return m_y; } 
const IdxTy seq() const { return m_serial;}
const IdxTy sclass() const { return m_class;}
private:
D  m_x, m_y, m_pitch;
IdxTy m_serial,m_class;
}; 

typedef std::map<IdxTy, seq_layout> SlMap;
typedef std::map<IdxTy, SlMap>  SlijMap;

typedef std::map<IdxTy, StrTy> Names;

public:
mjm_align_collection() {}
~mjm_align_collection() {}
const IdxTy&  bad() const { static IdxTy f= ~0U; return f;  } 
void add(const IdxTy & i, const IdxTy & j, const hit_vector & hv)
{ m_map[i][j]=hv; }

// TODO this is dumb the calleer uses j and krep mixing up i and j doh 
void add(const IdxTy & i, const IdxTy & j, const hit_vector & hv, const StrTy & iname, const StrTy & jname)
{ m_map[i][j]=hv; m_jnames[i]=iname; m_knames[j]=jname; }

void add(const IdxTy & i, const IdxTy & j, const hit_vector & hv, const StrTy & iname, const StrTy & jname, const StrTy & kseq, const StrTy & jseq)
{ m_map[i][j]=hv; m_jnames[i]=iname; m_knames[j]=jname; 
m_jseqs[i]=kseq; m_kseqs[j]=jseq;
}




void set_color_table( const Mct & ct) { m_ct=ct;  }

IdxTy size() const { return m_map.size(); }



void write_drift(const StrTy & fn , const IdxTy flags) 
{
std::ofstream os(fn);
 write_drift(os,flags);
}


void set_annotation(const Ragged & r)
{
m_annotation=r; 
}
void write_svg(const StrTy & fn , const IdxTy flags) 
{
std::ofstream os(fn);
 write_svg(os,flags);
}


template <class Tl > 
void layout_sequences(SlijMap&  seqs, HvijMap& hvij, Tl & lo )
{
IdxTy seqi=0;
MM_LOOP(ii,hvij)
{
const IdxTy & i=(*ii).first;
const HvMap & hm= (*ii).second;
const D x=lo.x(seqi);
const D y=lo.y(seqi);
const D pitch=lo.pitch(seqi);
seqs[i][bad()]= seq_layout(x,y,pitch,seqi,0);
++seqi;
MM_LOOP(jj,hm)
{
const IdxTy & j=(*jj).first;
if (seqs[i].find(j)!=seqs[i].end()) continue; 
//const hit_vector & hv=(*jj).second;
const D x=lo.x(seqi);
const D y=lo.y(seqi);
const D pitch=lo.pitch(seqi);
seqs[i][j]= seq_layout(x,y,pitch,seqi,1);

++seqi;
} // jj 
} // ii 


} // layout
class layout_info
{

public:
layout_info() { Init(); }
// 2022-09
void need_width(const D & w) { m_xs=w*m_pitch+m_lm+m_rm+10; }
void need_seqs(const D & h) { m_ys=h*m_ny+m_tm+m_bm+10; }
void need_height(const D & h) { m_ys=h+m_tm+m_bm+10; }
void width(const D & w) { m_xs=w; }
void height(const D & w) { m_ys=w; }

D xs( ) const { return m_xs; }
D ys( ) const { return m_ys; }
D x(const IdxTy i ) const { return i*m_nx+m_lm; }
D y(const IdxTy i ) const { return i*m_ny+m_tm; }
D pitch(const IdxTy & i ) const { return m_pitch; } 
D seq_name_x(const IdxTy & i) { return ( x(i)-.25*m_lm*(3-(i%3))); } 
const D & lm() const { return m_lm;} 
private:
void Init()
{
m_xs=2000;
m_ys=2000;
m_lm=100;
m_tm=200;
m_rm=100;
m_bm=300;
m_nx=1;
m_ny=250;
m_pitch=1;
} 
D m_xs,m_ys;
D m_lm,m_tm, m_nx,m_ny,m_rm,m_bm;
D m_pitch;
}; // layout_info;

#if 1 
template<class Td> 
void coverage_map(Td & m, const HvMap &hvm,const bool one)
{
MM_LOOP(ii,hvm)
{
//const IdxTy j=(*jj).first;
const hit_vector & hvv=(*ii).second;
MM_LOOP(jj,hvv)
{
const IdxTy base=(*jj).loc(one); 
for(IdxTy i=0; i<(*jj).len; ++i) 
{
++m[base+i];
} // i 
} // jj 
} // ii 
} // converage_map

//HvijMap m_map;
template <class Tv> 
void scan_hits( Tv & xx, Tv & yy, const HvijMap & m)
{

MM_LOOP(ii,m)
{
const IdxTy i=(*ii).first;
const HvMap & hvm=(*ii).second;
xx[i]=0;
yy[i]=0;

MM_LOOP(jj,hvm)
{
//const IdxTy j=(*jj).first;
const  hit_vector & hv=(*jj).second;
MM_LOOP(kk,hv)
{
const IdxTy extx=(*kk).j+(*kk).len-1;
const IdxTy exty=(*kk).i+(*kk).len-1;
if (extx>xx[i]) xx[i]=extx;
if (exty>yy[i]) yy[i]=exty;

} // kk 

} // jj 

} // ii 

} // scan_hits

// write out in format suitable for determining drift from
// the known sequence 
void write_drift(std::ofstream & os , const IdxTy flags) 
{
const bool offset_histo=(((1<<0)&flags)!=0);
const bool zig_zag=(((1<<1)&flags)!=0);
const bool list_strings=(((1<<2)&flags)!=0);
const bool xy_histo=(((1<<3)&flags)!=0);
const bool cover_histo=(((1<<4)&flags)!=0);
const int ylim=10;
MM_ERR("write_drift "<<MMPR4(flags,offset_histo,zig_zag,list_strings)<<MMPR3(xy_histo,cover_histo,ylim))
// this should be all of the alignment records
// organized by probe string
// there should only be one... 
typedef std::vector<int> Dr;
typedef std::map<Dr,IdxTy> DrMap;
typedef std::vector<int> Xy;
typedef std::map<Xy,Xy> XyMap;
// TODO FIXME these are in the wrong palce neeed to be inside the i loop 
XyMap xy,cover;
DrMap dr;
IdxTy in=0;
const IdxTy incsz=5;
MM_LOOP(ii,m_map)
{
//const IdxTy i=(*ii).first;
const IdxTy i=(*ii).first;
const StrTy target_name= m_jnames[i].substr(1,30);
const StrTy&  target_seq=m_jseqs[i];
const HvMap & hm=(*ii).second;
// put in entire string with dummy counts 
if (cover_histo){
const char * p = target_seq.c_str();
IdxTy pos=0;
while (*p!=0)
{
Xy loc(2),inc(incsz);
loc[0]=pos;
loc[1]=0;
inc[0]=IdxTy(*p);
cover[loc]=inc;
++pos;
++p;
} 


}

Ss ss;
//IdxTy jc=0;
// for test  sting j find all the alignment records 
IdxTy jn=0;
MM_LOOP(jj,hm)
{
const IdxTy j= (*jj).first;
const StrTy probe_name= m_knames[j].substr(1,30);
// these are the matches between i and j at all locations 
const hit_vector & hvv=(*jj).second;
// first find the mode
int minflat=0;
int mode=0;
std::map<int, int> dd;
MM_LOOP(kk,hvv) { dd[(*kk).delta(true)]+=(*kk).len; }
MM_LOOP(kk,dd)
{
if ( (*kk).second>minflat) {  minflat=(*kk).second; mode=(*kk).first; }
}
if (offset_histo)
{
MM_LOOP(kk,hvv) { 
// this is the location on the target
const int x=(*kk).j;
//const int extx=(*kk).j+(*kk).len-1;
const int len=(*kk).len;
// this is the location on the probe
const int y=(*kk).i;
//const int exty=(*kk).i+(*kk).len-1;
for (int p=0; p<len; ++p) 
{ Dr v(2); v[0]=x+p; v[1]=y+p-mode-v[0]; ++dr[v]; } // p 
}// kk 
} // offset_histo
if ( list_strings)
{
const IdxTy hvsz=hvv.size();
for(IdxTy ih=0; ih<hvsz; ++ih)
{
const auto & h=hvv[ih];
const int x=h.j;
const int len=h.len;
const int y=h.i-x-mode;
const StrTy & s=h.s;
os<<"list_strings "<<MMPR3(in,jn,ih)<<MMPR4(x,y,len,s)<<MMPR3(target_name,probe_name,mode)<<CRLF;
} // ih
} // list_strings

if ( xy_histo||cover_histo)
{
const IdxTy hvsz=hvv.size();
std::vector<IdxTy> ord(hvsz);
{MM_SZ_LOOP(i,ord,ordsz) { ord[i]=i; } }
std::sort(ord.begin(),ord.end(),
[&](const int& a, const int& b)
{ return hvv[a].j<hvv[b].j; });
int lasty=0;
int lastgood=0;
for(IdxTy ih=0; ih<hvsz; ++ih)
{
const IdxTy io=ord[ih];
const IdxTy iidx=io;
const auto & h=hvv[iidx];
const int x=h.j;
const int len=h.len;
const int y=h.i-x-mode;
//const StrTy & s=h.s;
Xy loc(2), inc(incsz);
loc[0]=x;
loc[1]=y;
if (cover_histo)
{
for(int  k=0; k<len; ++k)
{
if (cover.find(loc)!=cover.end()) { inc=cover[loc]; }
else inc=Xy(incsz);
inc[0]=IdxTy(target_seq.c_str()[x+k]);
++inc[1];
cover[loc]=inc;
++loc[0];
}
} // cover_histo
if (xy_histo) {

inc[0]=1;
inc[1]=len;
inc[2]=(y!=lasty)?1:0;
inc[1]=1;
const bool good_shift=(y!=lastgood);
const bool ygood= ((y>-ylim)&&(y<ylim));
inc[3]=(good_shift&&ygood)?1:0;
if (xy.find(loc)!=xy.end()) {
const auto & p=xy[loc];
inc[0]+=p[0]; inc[1]+=p[1]; inc[2]+=p[2]; inc[3]+=p[3]; inc[4]+=p[4]; }
xy[loc]=inc;
inc=Xy(5);
inc[4]=1;
for(int  k=1; k<len; ++k)
{
++loc[0];
if (xy.find(loc)==xy.end()) { xy[loc]=inc;}
else {xy[loc][4]+=inc[4]; } 

} // k 
 

//os<<"list_strings "<<MMPR3(in,jn,ih)<<MMPR4(x,y,len,s)<<MMPR2(target_name,probe_name)<<CRLF;
lasty=y;
if (ygood)  lastgood=y;

} // xy_histo


} // ih
} //xy_histo 





if (zig_zag)
{
const IdxTy hvsz=hvv.size();
std::vector<IdxTy> ord(hvsz);
{MM_SZ_LOOP(i,ord,ordsz) { ord[i]=i; } }
std::sort(ord.begin(),ord.end(),
[&](const int& a, const int& b)
{ return hvv[a].j<hvv[b].j; });

IdxTy lastx=0;
IdxTy lasty=0;
{MM_SZ_LOOP(i,ord,ordsz) { 
const auto & h=hvv[ord[i]];
const int x=h.j;
const int len=h.len;
const int y=h.i-x-mode;
const int xend=x+len-1;
const int yend=y;

//{ Dr v(2); v[0]=x+p; v[1]=y+p-mode-v[0]; ++dr[v]; } // p 


//if (lastx!=0) 
{os<<"zig_zag "<<MMPR3(in,jn,i)<<MMPR2(x,y)<<CRLF;}
os<<"zig_zag "<<MMPR3(in,jn,i)<<MMPR2(xend,yend)<<CRLF;

lastx=xend;
lasty=yend;

 } } // i 

if (false ) MM_ERR(MMPR2(lastx,lasty))
} // zig_zag


++jn;
} // jj 
++in;
} // ii 
 // TODO FIXME these are in the wrong place but work ok as long
// as there is only ONE target sequence. 

if (cover_histo)
{
MM_LOOP(ii,cover)
{
const Xy & loc=(*ii).first;
const Xy & cnt=(*ii).second;
const IdxTy & x =(loc[0]);
const int  & y=loc[1];
const char c =char( cnt[0]);
const IdxTy & n=cnt[1];
//const IdxTy & shift=cnt[2];
//const IdxTy & shiftlim=cnt[3];
//const IdxTy & covers=cnt[4];

//os<<"xy_histo "<<MMPR4(x,y,n,nlen)<<MMPR3(shift,shiftlim,covers)<<CRLF;
os<<"cover_histo "<<MMPR4(x,y,c,n)<<CRLF;

} // ii 




} // cover_histo
if (offset_histo) { 
MM_LOOP(ii,dr)
{
const auto & v=(*ii).first;
const auto & cnt=(*ii).second;
os<<MMPR3(v[0],v[1],cnt)<<CRLF;
}
} // offset_histo
if ( xy_histo)
{
MM_LOOP(ii,xy)
{
const Xy & loc=(*ii).first;
const Xy & cnt=(*ii).second;
const IdxTy & x=loc[0];
const int  & y=loc[1];
const IdxTy & n=cnt[0];
const IdxTy & nlen=cnt[1];
const IdxTy & shift=cnt[2];
const IdxTy & shiftlim=cnt[3];
const IdxTy & covers=cnt[4];

os<<"xy_histo "<<MMPR4(x,y,n,nlen)<<MMPR3(shift,shiftlim,covers)<<CRLF;

} // ii 


} // xy_histo


} // drift


typedef  Ragged::Line NoteLine;
typedef std::map<IdxTy,StrTy> NoteMap;
void add_to(  NoteMap & m, const IdxTy & szl, const NoteLine & l, const IdxTy flags)
{
if (flags==0) 
{
IdxTy base=0;
IdxTy cnt=0;
StrTy col="";
if (szl>1) base=myatoi(l[1]);
if (szl>2) cnt=myatoi(l[2]);
if (szl>3) col=(l[3]);
for(IdxTy i=0; i<cnt; ++i) { m[base+i]=col ; } 
}



}

class _write_svg_params
{

public:
_write_svg_params() { Init(); }


void Init()
{
display_target_bases=true;
sort_by_goodness=true;
label_names=true;
show_baseline=true;
write_annotation=true;
} // Init

bool display_target_bases;
bool sort_by_goodness;
bool label_names;
bool show_baseline;
bool write_annotation;


}; // _write_svg_params
typedef _write_svg_params  write_svg_params;

void write_svg(std::ofstream & os , const IdxTy flags) 
{
write_svg_params wsp;
write_svg(os,wsp,flags);
}

void write_svg(std::ofstream & os , const write_svg_params & wsp, const IdxTy flags) 
{
layout_info li;
SlijMap seqs;
const IdxTy dpflags=flags;
MM_FLAG(remove_diag,4,1);

const bool display_target_bases=wsp.display_target_bases;
const bool sort_by_goodness=wsp.sort_by_goodness;
const bool label_names=wsp.label_names;
const bool show_baseline=wsp.show_baseline;
const bool write_annotation= wsp.write_annotation;


NoteMap  base_colors;
NoteMap  reference_lines; 

if (write_annotation)
{
const IdxTy sz=m_annotation.size();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line & l=m_annotation.line(i);
const IdxTy szl=l.size();
if (szl<1) continue;
MM_ERR(MMPR2(l[0],szl))
const StrTy & cmd=l[0];
if (cmd=="base-color") add_to(base_colors,szl,l,0);
if (cmd=="reference-line") add_to(reference_lines,szl,l,0);

} // i 

} // write_annotation


std::map<IdxTy,IdxTy> maxx,maxy;
scan_hits(maxx,maxy,m_map);

layout_sequences(seqs,m_map,li);
// 2022-09
// y is more complicated, need all the sequences but large offset
// hits can be a problem 
IdxTy maxX=0;
IdxTy maxY=0;
MM_LOOP(ii,maxx) { if ((*ii).second>maxX) maxX=(*ii).second; } // ii 
MM_LOOP(ii,maxy) { if ((*ii).second>maxY) maxY=(*ii).second; } // ii 
const bool set_active_area=true;
if (set_active_area)
{
li.need_width(maxX);
// this is probably the probe size not the realones doh 
// ASSFUICK 
MM_ERR(MMPR2(seqs.size(),maxY))
li.need_seqs(seqs.size());
li.need_height(maxY);

}



//const bool omit_seg_labels=((flags&1)!=0);


if (m_ct.size()==0)
{
//2071  ./mjm_string_seq.out -cmd "set-param mflags 2"  -cmd "read-fasta knowns limis_k2.fasta" -cmd "index-fasta knowns 8"  -cmd "set-param maxscores 15" -cmd "set-param save_aligns 1"  -cmd "mt-explore-streaming-hits knowns ku.fasta" -cmd "cmd-align default some_colors" -cmd "write-align xxx.svg default"  -quit # | tee linear_combs_cross._txt



MM_ERR(" color table is blanks, consider a cmd like cmd-align default \\\"some_colors\\\"   ")

MM_ERR(" setting to default, to grey out all set one to grey... ")
m_ct.some_colors();
}






const IdxTy xs=li.xs(); // 100;
const IdxTy ys=li.ys(); // 100;
mjm_svg_writer sw;
os<<sw.start_text(" test foo",xs,ys);
os<<sw.frame_text("#00808080",xs,ys);
os<<CRLF;
Ss sstemp;
MM_LOOP(ii,m_map)
{
const IdxTy i=(*ii).first;
const HvMap & hm=(*ii).second;
Ss ss;
ss<<"j"<<i<< m_jnames[i].substr(1,30);
const D xx=li.seq_name_x(seqs[i][bad()].seq()); // base_x(0);
const D yy=seqs[i][bad()].base_y(0);
os<<sw.gtext_text(ss.str(),xx,yy,20,StrTy("red"),StrTy("left"),90); os<<CRLF;

{
const IdxTy full=maxy[i];
std::map<IdxTy, IdxTy> m;
coverage_map(m, hm,!true);
IdxTy last=0;
IdxTy l2=0;
IdxTy lbase=0;
MM_LOOP(mm,m)
{
const IdxTy base=(*mm).first;
const IdxTy cnt=(*mm).second;
const D xi=seqs[i][bad()].base_x(base);
const D yi=seqs[i][bad()].base_y(base);
{os<<sw.line_text(xi,yi,xi,yi-cnt, StrTy("blue"),1,1); os<<CRLF; 
Ss ss; ss<<cnt;
//const D kluge=0.5;
const D kluge=-0.5;
MM_ONCE(" viewer dependent kluge set to "<<MMPR(kluge),)
// This is probably then "off by one " 
// #error check offset 
//sstemp<<sw.label_line_text(ss.str(),xi-kluge,yi,xi-kluge,yi-cnt,1,StrTy("red"),"start",.1);
sstemp<<sw.label_line_text(ss.str(),xi-kluge,yi,xi-kluge,yi-cnt,1,StrTy("white"),"middle",.1);
//{os<<sw.line_text(xi,yi,xi,yi-cnt, StrTy("white"),th,opaq); os<<CRLF; }
}
D xlo=.5;
if (base!=(lbase+1))
{ 
const D filldefi=.05;
{os<<sw.line_text(xi+filldefi,yi,xi+filldefi,yi+full, StrTy("black"),1.0-2*filldefi,.9); os<<CRLF; } 
MM_ONCE(" new labels for black lines",)
Ss ss;
ss<<MMPR(base);
os<<sw.label_line_text(ss.str(),xi+xlo,yi,xi+xlo,yi+full,1,StrTy("white"),"start",.001);

}
if (cnt>0) if (last>cnt) if ( l2<last)
{ 
{os<<sw.line_text(xi,yi,xi,yi+full, StrTy("red"),1,.9); os<<CRLF; }
Ss ss;
ss<<MMPR(base);
os<<sw.label_line_text(ss.str(),xi+xlo,yi,xi+xlo,yi+full,1,StrTy("white"),"start",.001);
 }
lbase=base;
l2=last;
last=cnt;

}
}

const bool do_reference=true;
if (do_reference)
{
const IdxTy full=maxy[i];
const StrTy & seq=m_jseqs[i];
//const char * p =seq.c_str();
const IdxTy sz=seq.length();
for(IdxTy ip=0; ip<sz; ++ip)
{
const D xi=seqs[i][bad()].base_x(ip);
const D yi=seqs[i][bad()].base_y(ip);
const auto iirl=reference_lines.find(ip);
if (iirl!=reference_lines.end())
{
const D filldefi=.0;
os<<sw.line_text(xi+filldefi,yi,xi+filldefi,yi+full, (*iirl).second,1.0-2*filldefi,.9); 
os<<CRLF;  
}
} // ip 
} // do_reference

// in the blus histogram at top 
if ( display_target_bases)
{
const StrTy & seq=m_jseqs[i];
const char * p =seq.c_str();
const IdxTy sz=seq.length();
MM_ERR( " writing target bases "<<MMPR(sz))
for(IdxTy ip=0; ip<sz; ++ip)
{
const D xi=seqs[i][bad()].base_x(ip);
const D yi=seqs[i][bad()].base_y(ip);
//{os<<sw.line_text(xi,yi,xi,yi-cnt, StrTy("blue"),1,1); os<<CRLF; 
Ss ss; ss<<p[ip];
const D kluge=0;
StrTy bcol=base_colors[ip];
if(bcol.c_str()[0]==0) bcol="white";
//sstemp<<sw.label_line_text(ss.str(),xi-kluge,yi,xi-kluge+3,yi-cnt,1,StrTy("white"),"start",.1);
sstemp<<sw.gtext_text(ss.str(),xi-kluge,yi,1,bcol,StrTy("middle"),0); os<<CRLF;

if (( ip%10)==0) sstemp<<sw.gtext_text("0",xi-kluge,yi+3,1,StrTy("white"),StrTy("middle"),0); os<<CRLF;


} // ip 


} // display_target_bases

/////////////////////////////////////////////////////

std::vector<IdxTy> sord;
std::map<IdxTy, IdxTy> sohits;
if (sort_by_goodness) { 
MM_LOOP(jj,hm) { sord.push_back((*jj).first);
std::map<IdxTy, IdxTy> szm;
const hit_vector & hvv=(*jj).second;
MM_LOOP(kk,hvv)
{
for(IdxTy k=(*kk).j; k<((*kk).j+(*kk).len); ++k) ++szm[k];
}
 sohits[(*jj).first]=(szm.size()); 
} 


std::sort(sord.begin(),sord.end(),
[&](const int& a, const int& b)
{ return sohits[a]>sohits[b]; });


}
sohits.clear();
MM_SZ_LOOP(l,sord,sosz) sohits[sord[l]]=l;

////////////////////////////////////////////////////
IdxTy jc=0;
//const bool xj=( sort_by_goodness)&& (remove_diag);
MM_LOOP(jj,hm)
{

const IdxTy j= (*jj).first;
//const IdxTy jeff=(xj)?sord[(*jj).first]: (*jj).first;


const hit_vector & hvv=(*jj).second;
if (!remove_diag)
{
Ss ss;
ss<<"krep"<<j<< m_knames[j].substr(1,30);
const D xx=li.seq_name_x(seqs[i][j].seq()); // base_x(0);
const D yy=seqs[i][j].base_y(0);
///ss<<std::setprecision(3)<<(val); // (val*(zmax-zmin)+zmin); 
//os<<sw.gtext_text(ss.str(),xv.back()+.5*legp,yv.back()+hleg+10,ssz,cv.back(),StrTy("left"),45); os<<CRLF;
if (false) 
{ os<<sw.gtext_text(ss.str(),xx,yy,10,StrTy("green"),StrTy("left"),45); 
os<<CRLF; }

os<<sw.gtext_text(ss.str(),xx,yy,10,StrTy("green"),StrTy("left"),90); os<<CRLF;


}
D minflat=0;
if (remove_diag)
{
D  mode=0;
std::map<int, int> dd;
MM_LOOP(kk,hvv) { dd[(*kk).delta(true)]+=(*kk).len; }
MM_LOOP(kk,dd)
{
if ( (*kk).second>minflat) {  minflat=(*kk).second; mode=(*kk).first; } 
}
if( sort_by_goodness)
{

minflat=50.0-mode+10.0*sohits[j]; // sord[j];
MM_ERR(MMPR4(minflat,j,jc,sord[j])<<MMPR2(sohits[j],sord[jc]))
}


else minflat=50.0-mode+10.0*j;

if ( label_names)
{
const StrTy col=m_ct(j); // 
const D xz=seqs[0][bad()].base_x(0)-.5*li.lm();
const D yz=seqs[0][bad()].base_y(0)+minflat+mode-2;
Ss ss;
if ( show_baseline) 
{ os<<sw.fainterline_text(xz+.5*li.lm(),yz,xz+xs,yz, col,.5); os<<CRLF; } 
//{ os<<sw.line_text(xz+.5*li.lm(),yz,xz+xs,yz, col,3); os<<CRLF; } 
ss<<"#"<<sohits[j]<<"k"<<j<< m_knames[j].substr(1,60);
//const D xx=li.seq_name_x(seqs[i][j].seq()); // base_x(0);

//const D offseti=(remove_diag)?(minflat-hv.j):0;
//const D offsetf=(remove_diag)?(len):0;
/// TODO FIXME  this needs to use the pitch or something similar 
//const D  x0=xz+hv.j;
//const D  y0=yz+hv.i+offseti;
// these are line ends NOT recangle stuff
//const D x1=x0+len-0;
//const D y1=y0+len-0-offsetf;
//const D opaq=.7;
//const D th=3;
os<<sw.gtext_text(ss.str(),xz,yz,4,col,StrTy("start"),0); os<<CRLF;

} // labeal_names
}


MM_LOOP(kk,hvv)
{
const align_type & hv=(*kk);
if (false) { connect(sstemp, os, sw,  seqs, hv,  i, j ); } 
dotplot(sstemp, os, sw,  seqs, hv,  i, j,minflat,dpflags );





} // kk 
++jc;
} // jj 
} // ii 

//os<<sw.label_line_text(ss.str(),x,y,x1,y1,30,StrTy("white"),"left",.3);
//os<<CRLF;
//if (w<.9) {os<<sw.line_text(x,y,x1,y1, StrTy("blue"),th); os<<CRLF; }
//else if (w<.95){ os<<sw.line_text(x,y,x1,y1, StrTy("green"),th); os<<CRLF;}
//else { os<<sw.line_text(x,y,x1,y1, StrTy("red"),th); os<<CRLF; }
//if (!omit_seg_labels) 
{ os<<sstemp.str(); } 
os<<sw.end_text();
os<<CRLF;
} // write_svg ?


template <class Oss,class Os, class Sw,class Se > 
void dotplot(Oss & sstemp, Os & os, Sw & sw, Se & seqs,const align_type & hv,  const IdxTy i, const IdxTy j,const D & minflat, const IdxTy flags )
{
MM_FLAG(label_seg_perp,0,0)
MM_FLAG(show_seg,1,0)
MM_FLAG(label_seg,2,1)
MM_FLAG(label_seg_seq,3,1);
MM_FLAG(remove_diag,4,1);
const bool include_id=true;
//const D minflat=200;
const IdxTy len=hv.len;
const D offseti=(remove_diag)?(minflat-hv.j):0;
const D offsetf=(remove_diag)?(len):0;
//const IdxTy len2=len>>1;
//const IdxTy jmid=hv.i+(len2);
//const IdxTy imid=hv.j+(len2);
const D xz=seqs[0][bad()].base_x(0);
const D yz=seqs[0][bad()].base_y(0);
/// TODO FIXME  this needs to use the pitch or something similar 
const D  x0=xz+hv.j;
const D  y0=yz+hv.i+offseti;
// these are line ends NOT recangle stuff
const D x1=x0+len-0;
const D y1=y0+len-0-offsetf;
const D opaq=.7;
const D th=3;
const StrTy col=m_ct(j); // 
//const bool show_seg=((flags&(1<<0))==0);
//const bool label_seg=((flags&(1<<1))==0);
//const bool label_seg_perp=((flags&(1<<2))==0);
//const bool label_seg_seq=((flags&(1<<3))==0);
if (show_seg)
{os<<sw.line_text(x0,y0,x1,y1, col,th,opaq); os<<CRLF; }
if (label_seg) {
Ss ss;
ss<<"("<<i<<","<<j<<","<<len<<","<<hv.delta(true)<<")";
sstemp<<sw.label_line_text(ss.str(),x0,y0,x1,y1,3,StrTy("white"),"left",.01);
sstemp<<CRLF; }
if (label_seg_seq) { 
Ss ss;
if ( include_id) { ss<<MMPR(j)<<"d"<<hv.delta(true); } 
ss<<hv.s;
sstemp<<sw.label_line_text(ss.str(),x0,y0,x1,y1,1,StrTy("white"),"left",.01);
sstemp<<CRLF; }


if (label_seg_perp) {
Ss ss;
ss<<"("<<i<<","<<j<<","<<len<<","<<hv.delta(true)<<")";

sstemp<<sw.gtext_text(ss.str(),x1,y1,3,StrTy("white"),StrTy("start"),-45); 
sstemp<<CRLF;}



} // dotplot

template <class Oss,class Os, class Sw,class Se > 
void connect(Oss & sstemp, Os & os, Sw & sw, Se & seqs,const align_type & hv,  const IdxTy i, const IdxTy j )
{
const IdxTy len=hv.len;
const IdxTy len2=len>>1;
const IdxTy jmid=hv.i+(len2);
const IdxTy imid=hv.j+(len2);
const D x0=seqs[i][bad()].base_x(imid);
const D y0=seqs[i][bad()].base_y(imid);
const D x0l=seqs[i][bad()].base_x(hv.j);
const D y0l=seqs[i][bad()].base_y(hv.j);
const D x0r=seqs[i][bad()].base_x(hv.j+len);
const D y0r=seqs[i][bad()].base_y(hv.j+len);


const D x1=seqs[i][j].base_x(jmid);
const D y1=seqs[i][j].base_y(jmid);
const D x1l=seqs[i][j].base_x(hv.i);
const D y1l=seqs[i][j].base_y(hv.i);
const D x1r=seqs[i][j].base_x(hv.i+len);
const D y1r=seqs[i][j].base_y(hv.i+len);


// these are presentation or view issues not properties of the
// alignment markers. 
const D dx=x1-x0;
const D dy=y1-y0;
const bool retro=(dx*dy<0);
D opaq=(len>60)?.1:.2;
D th=len;
Ss ss;
ss<<hv.len;
 StrTy col=(retro)?StrTy("red"):StrTy("blue");
// wide lines will have angled bottoms, need rect or array of
// stakes... 
if (len<10){os<<sw.line_text(x0,y0,x1,y1, col,th,opaq); os<<CRLF; }
else
{
col="white";
th=1; opaq=.6;
{os<<sw.line_text(x0r,y0r,x1r,y1r, col,th,opaq); os<<CRLF; }
{os<<sw.line_text(x0l,y0l,x1l,y1l, col,th,opaq); os<<CRLF; }
{os<<sw.line_text(x0l,y0l,x0r,y0r, col,th,opaq); os<<CRLF; }
{os<<sw.line_text(x1l,y1l,x1r,y1r, col,th,opaq); os<<CRLF; }


}

{sstemp<<sw.label_line_text(ss.str(),x1,y1,x0,y0,3,StrTy("white"),"left",.01);
sstemp<<CRLF; }



} // connect 

void cmd(const StrTy & cmd, const IdxTy flags)
{ CommandInterpretter li; li.set(cmd,1);

command_mode(li);

}

void command_mode(CommandInterpretter &li)
{
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="some_colors") { m_ct.some_colors(); continue; }
if (cmd=="clear") { m_ct.clear(); continue; }
if (sz==2)
{
if (cmd=="add") { m_ct.add(li.word(1)); continue; }
if (cmd=="load") { m_ct.load(li.word(1),0); continue; }
} // 1 
if (sz==3)
{
if (cmd=="set") { m_ct.set(myatoi(li.word(1)),li.word(2)); continue; }
if (cmd=="load") { m_ct.load(li.word(1),myatoi(li.word(2))); continue; }
} // 1 


}

}






#endif

private:

HvijMap m_map;
Names m_jnames,m_knames;
Names m_jseqs,m_kseqs;
Mct m_ct;
Ragged m_annotation;



}; // mjm_align_collection

// }; //mjm _



// #undef MM_DMEL

#endif

