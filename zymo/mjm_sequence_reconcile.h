#ifndef MJM_SEQUENCE_RECONCILE_H__
#define MJM_SEQUENCE_RECONCILE_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_fasta_ii.h"
#include "mjm_instruments.h"
#include "mjm_strings.h"
#include "mjm_collections.h"
#include "mjm_kv_vector.h"

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>


// Mon Aug 15 14:01:37 EDT 2022
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_sequence_reconcile   
// g++  -Wall -std=gnu++11 -DTEST_MJM_SEQUENCE_RECONCILE -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_sequence_reconcile.h  -o mjm_sequence_reconcile.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_sequence_reconcile("mjm_sequence_reconcile"
, "  ");

template <class Tr>
class mjm_sequence_reconcile 
{
 typedef mjm_sequence_reconcile Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef Ragged::Line Line;
typedef std::vector<Line> Lines;
typedef mjm_kv_vector<StrTy,D,Tr> SumVec;
typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy , Fasta> FastaMap;

// a set of fastas 
class _Entry
{
public:
StrTy dump() const
{
Ss ss;
ss<<MMPR2(fasta,name);
MM_LOOP(ii,tax) ss<<" "<<(*ii); 
return ss.str();
}
StrTy fasta,name;
std::vector<StrTy> tax;
}; // _Entry

typedef _Entry Entry;

typedef std::vector<Entry>  EntryVector;
typedef std::map<StrTy,Entry>  EntryMap;
// the key is a sequence full length 
typedef std::map<StrTy, EntryVector> MapTy;
typedef std::map<StrTy, EntryMap> MapMapTy;

typedef std::map< StrTy, std::map< StrTy, IdxTy > > IndexTy;
typedef std::vector<StrTy> AssBaseTy;

class _AssTy : public AssBaseTy
{
public:
_AssTy() : AssBaseTy() {}
_AssTy(const AssBaseTy & x): AssBaseTy(x) {}
const StrTy & fasta() const { return m_fasta;}
const StrTy & fasta(const StrTy & x) {m_fasta=x;  return m_fasta;}
StrTy m_fasta;

};
typedef _AssTy AssTy;

typedef std::map< StrTy,  AssTy > AssMapTy;

typedef std::map<StrTy, AssMapTy> AssMapMap;

//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_sequence_reconcile() {}
~mjm_sequence_reconcile() {}
// load and index a fasta file. 
// nm is the fasta map name and fn is the file, append
// to whatever is there. index by nm
IdxTy load(const StrTy & nm, const StrTy & fn, const IdxTy flags)
{ return Load(nm,fn,flags); } 
IdxTy load_assignments(const StrTy & nm, const StrTy & fn, const IdxTy flags)
{ return LoadAss(nm,fn,flags); } 

IdxTy dump_fasta(Ragged & r,  const IdxTy flags) { return DumpFasta(r,flags); }


typedef AssTy assign_type;
AssTy reconcile(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
{ return Reconcile(f,seqno,ass,flags); }

AssTy lookup(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
{ return Lookup(f,seqno,ass,flags); }
template <class LiTy> IdxTy find_olaps( Ragged & r, const LiTy & bl, const IdxTy redlev, const IdxTy flags)
{ return  FindOlaps(r, bl, redlev, flags); }



StrTy dump_eq(const IdxTy flags=0) { return DumpEq(flags); }
StrTy dump_ass(const IdxTy flags=0) { return DumpAss(flags); }
StrTy dump_ass_fasta(const IdxTy flags=0) { return DumpAssFasta(flags); }
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {
Ss ss;  
MM_LOOP(ii,m_map)
{
const StrTy & nm=(*ii).first;
ss<<MMPR2(nm,(*ii).second.size())<<CRLF;
MM_LOOP(jj,(*ii).second)
{
ss<<MMPR2((*jj).fasta,(*jj).name);
ss<<CRLF;
} // jj 
} // m_map
MM_LOOP(ii,m_fasta_map)
{
const StrTy & nm=(*ii).first;
const auto & f=(*ii).second;
ss<<MMPR3(nm,f.name(),f.size())<<CRLF;

} // ii 
return ss.str(); 

} // Dump
StrTy DumpEq(const IdxTy flags=0) {
Ss ss;  
MM_LOOP(ii,m_map)
{
const StrTy & nm=(*ii).first;
ss<<MMPR2(nm,(*ii).second.size())<<CRLF;
ss<<" eqsize= "<<(*ii).second.size()<<" ";
MM_LOOP(jj,(*ii).second)
{
StrTy& seq=(*jj).name;
const StrTy & fasta=(*jj).fasta;
ss<<MMPR2(fasta,seq);
} // jj 
ss<<CRLF;
} // m_map
return ss.str(); 
} // DumpEq
IdxTy DumpFasta(Ragged & r,  const IdxTy flags) 
{ 
IdxTy rc=0;

MM_LOOP(ii,m_fmap)
{
const StrTy & seq=(*ii).first;
Line lname(1);
lname[0]=">recon ";
StrTy & x=lname[0];
MM_LOOP(jj,(*ii).second)
{
const StrTy & fname=(*jj).first;
const Entry & e=(*jj).second;
x=x+e.fasta+" "+e.name+" ";
MM_LOOP(kk,e.tax) x+=(*kk)+" ";
x+=";";
} // jj 
const auto ll=m_ass.find(seq);
if (ll!=m_ass.end())
{
Ss ss;
ss<<" ASSIGNED "<<(*ll).second.fasta();
MM_LOOP(kk,(*ll).second) ss<<(*kk)<<" ";
lname.push_back(ss.str());
}
r.add(lname);
Line s;
s.push_back(seq);
r.add(seq);
} // ii 


#if 0
// sequence should only appear once per fasta... 
Entry & ev=m_fmap[seq][f];
//StrTy fasta,name; std::vector<StrTy> tax; }; // _Entry
ev.fasta=f; ev.name=seqno; ev.tax=ass;
// rbegin finds the most recent if the names start lexi chrono
//const AssTy &  best=(*(m_fmap[seq].rbegin())).second.tax;
AssTy   best=(*(m_fmap[seq].rbegin())).second.tax;
best.fasta(f);
const bool best_seq=IsSeq(best);
const bool best_tax=HasTax(best);
if ( best.size()) { m_ass[seq]=best; return best ; }
#endif 




return rc; 
}// DumpFasta

StrTy DumpAssFasta(const IdxTy flags=0) {
Ss ss;
AssTy ASSFUC;
//Fasta & x=m_fasta_map[nm];
//const AssTy  & ass=m_lut[fasta][seq];
MM_LOOP(ii,m_fasta_map)
{
const StrTy & key=(*ii).first;
const auto & asslut=m_lut[key];
Fasta & x=(*ii).second;
const IdxTy sz=x.size();
for(IdxTy i=0; i<sz; ++i)
{
// this  name has multiple  fields wtf 
const StrTy & ni=x.name(i);
const char * p = ni.c_str();
StrTy asfck;
char cass[ni.length()+1];
IdxTy sp=0;
IdxTy dp=0;
if (*p=='>') sp=1;
while (p[sp]!=0)
{
if (p[sp]==' ') break;
if (p[sp]=='\t') break;
cass[dp]=p[sp];
++dp;
++sp;
} // while
cass[dp]=0;

StrTy nf=StrTy(cass);
const StrTy & seq=x.seq(i);
//auto ff=asslut.find((*p=='>')?StrTy(p+1):ni);
auto ff=asslut.find(nf);
const AssTy & av=(ff==asslut.end())?ASSFUC:(*ff).second;
MM_ERR(" WTF "<<MMPR(ni)<<MMPR4(nf,asslut.size(),key,x.name()))
Ss rr;
rr<<x.name()<<" ";
MM_LOOP(jj,av) { rr<<(*jj)<<" "; }
Ss tt;
//ss<<ni<<" "<<rr.str()<<CRLF;
if (*p != '>') tt<<">";
tt<<ni<<" "<<rr.str(); // <<CRLF;
// wtf clustalw chokes on spaces in thef cuking name  
// this pointer to the  string does not seem to 
// work? Its not perm? WTF 
//const char * p2 = tt.str().c_str();
const StrTy shit=tt.str();
const char * p2 = shit.c_str();
//char cassf[tt.str().length()+1];
char cassf[shit.length()+2];
MM_ERR(MMPR(shit.length()))
sp=0;
dp=0;
//if (*p2=='>') sp=1;
while (p2[sp]!=0)
{
if (p2[sp]==' ') { cassf[dp]=';';   } 
else cassf[dp]=p2[sp];
//MM_ERR(MMPR4(sp,dp,p2[sp],cassf[dp]))
++dp;
++sp;
} // while
cassf[dp]=0;

ss<<StrTy(cassf)<<CRLF;
//ss<<tt.str()<<CRLF;

ss<<seq<<CRLF;
} // i 

} // ii 
return ss.str();
} // DumpAssFasta


StrTy DumpAss(const IdxTy flags=0) {
const bool one_liner=Bit(flags,0);
const bool two_liner=Bit(flags,1);
IdxTy gseqno=0;
Ss ss;  
// the m_map  is indexed by sequence  
//m_map[seq].push_back(x);
// while the lut is by project
//auto & at=m_lut[nm];
const bool idcrlf=!two_liner&&!one_liner;
const bool mcr= (two_liner||!one_liner);
MM_ERR(MMPR4(one_liner,two_liner,idcrlf,mcr)<<MMPR(flags))
MM_LOOP(ii,m_map)
{
const StrTy & nm=(*ii).first;
ss<<MMPR3(gseqno,nm,(*ii).second.size());
if (mcr) ss<<CRLF;
if (two_liner) ss<<MMPR(gseqno);
if ( two_liner ) { ss<<" eqsize= "<<(*ii).second.size()<<" "; } 
std::vector<AssTy> past;
std::vector<StrTy> pastid;
past.clear();
Ss fs;
MM_LOOP(jj,(*ii).second)
{
StrTy& seq=(*jj).name;
const StrTy & fasta=(*jj).fasta;
fs<<fasta;
if (idcrlf) {  ss<<MMPR(gseqno); } 
ss<<MMPR2(fasta,seq);
const AssTy  & ass=m_lut[fasta][seq];
if ( past.size()) if (ass!=past.back())
{ 
Ss rr; MM_LOOP(kk,ass) rr<<(*kk)<<" ";
rr<<" was  ";
MM_LOOP(kk,past.back()) rr<<(*kk)<<" ";
MM_ERR( " assignment changed  was "<<pastid.back()<<" now "<<MMPR2(fasta,seq)<<" "<<rr.str()) 
 } 
past.push_back(ass);
{ Ss rr;  rr<<MMPR2(fasta,seq)<<" ";
pastid.push_back(rr.str()); }

MM_LOOP(kk,ass) { ss<<" "<<(*kk); }
if (idcrlf ) { ss<<CRLF; }  
if (mcr ) { ss<<CRLF; }  
} // jj 
if (!idcrlf) ss<<CRLF;
ss<<CRLF<<MMPR(fs.str())<<CRLF;
++gseqno;
} // m_map
return ss.str(); 
} // DumpAss



typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

IdxTy Load(const StrTy & nm, const StrTy & fn, const IdxTy flags)
{ 
// this does not allow on the fly indexing.. 
Fasta & x=m_fasta_map[nm];
const IdxTy szi=x.size();
x.load(fn);
// the name will be known at look up, set the name to fn
x.name(fn);
const IdxTy szf=x.size();
for(IdxTy i=szi; i<szf; ++i)
{
Add(nm,i,x.seq(i),x.name(i),flags);
} // i 
MM_ERR(MMPR4(__FUNCTION__,nm,fn,x.name())<<MMPR3(flags,x.preamble(),x.size()))
return 0; 
}  // Load
// name is the sequence name, seq is the actual sequence, nm is the fasta name, and
// idx is the seuences position in the fasta. 
void Add(const StrTy & nm, const IdxTy idx, const StrTy & seq, const StrTy & name,  const IdxTy flags)
{
// zymo convention >seqxxx count=nnn
const bool first_word=Bit(flags,0);
Entry x;
x.name=name;
if (first_word)
{

char c[name.length()+1];
char * p=c;
memcpy(c,name.c_str(),name.length()+1);
while ( *p!=0)
{
if (*p==' ') break;
++p;
} // p 
*p=0;
x.name=StrTy(c);
if (c[0]=='>') x.name=StrTy(c+1);
} // first_word;

x.fasta=nm;
m_map[seq].push_back(x);
// now update index for lookup
MM_ERR(MMPR4(first_word,nm,name,x.name)<<"FUCL")
// fasta name and sequence name point to location in the fasta 
m_idx[nm][x.name]=idx;

// optionally maintain a contains data base... 


// create a null assignment... 
//Reconcile(nm,x.name,AssTy(),0);
} // Add

// extracted from Reconcile
AssTy Lookup(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
{
const bool append_seq=Bit(flags,2);
// the "flag" is actually a pass count and not updating on second pass
//const bool update=!Bit(flags,0);
//const bool presume_lexi_chrono=true;
// this contaminates the index  if seqno is not found...
// The fasta or sequence may be missing.
//}

// optionally check to see if f/seqno contains a smaller sequence
// that is the root due to inconsistent trimming. 
/// in that case, the better assignment of this or the smaller needs to be considered
// or munged. 
//const auto cc=m_contains.find(Entry(f,seqno));

bool found=false;
auto ii=m_idx.find(f);
bool fasta_found=(ii!=m_idx.end());
if (!fasta_found)
{MM_ERR(" no index found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno)) return AssTy(); }
auto kk=m_fasta_map.find(f);
if (kk==m_fasta_map.end())
{MM_ERR(" no fasta found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno)) return AssTy();
}
const Fasta & ff=(*kk).second; // m_fasta_map[f];
const IdxTy ffsz=ff.size();
// apparently m_idx finds the seuqnece index for fasta f
const auto & sm=(*ii).second;
// this apparentl  returns the sequence location for the name seqno 
auto jj=sm.find(seqno);
if (jj!= sm.end()) found=true;
if (!found)
{MM_ERR(" no sequence found "<<MMPR2(found,fasta_found)<<MMPR4(f,ff.name(),ffsz,seqno))
return AssTy();
}
const IdxTy loc=(*jj).second; // m_idx[f][seqno];
//const StrTy  seq=m_fasta_map[f].seq(loc);
const StrTy  seq=(*kk).second.seq(loc);
AssTy old=m_ass[seq];
// this makes it too long lol 
if (append_seq)  // old.push_back(seq);
{
if (old.size()) old.back()+=" "+seq;
else old.push_back(seq); 

}
//if (!update) return old;
return old;
} // Lookup
class Floc { public:
 Floc( ): f(),n(~0),dn(0) {}
 Floc(const StrTy _f, IdxTy _n ): f(_f),n(_n),dn(0) {}
 Floc(const StrTy _f, IdxTy _n, IdxTy _dn ): f(_f),n(_n),dn(_dn){ }
bool same_seq(const Floc & that) const { return (f==that.f)&&(n==that.n); } 
 StrTy f; IdxTy n; IdxTy dn; 
OsTy & dump(OsTy & os, const IdxTy flags) const 
{ os<<MMPR3(f,n,dn); } // dump
Ss & dump(Ss & os, const IdxTy flags) const 
{ os<<MMPR3(f,n,dn); return os; } // dump



StrTy dump(const IdxTy flags)const { Ss ss; return dump(ss,flags).str(); } // dump
bool operator<(const Floc & that ) const 
{ return (f<that.f)||((f==that.f)&&(n<that.n)); }
}; // floc
typedef std::vector<Floc>  Flocs;
typedef std::map<int,Flocs> LenLocs;
typedef std::map<Floc,Flocs> OlapIdx;
// TODO this is probably better indexed on exact full length
// sequence but the trimming combination are n1*n2 as each 
// end is indepenent and tested alone. 
template <class LiTy> IdxTy IndexByLength( LenLocs& lenidx, const LiTy & bl, const IdxTy flags)
{ 
MM_LOOP(ii,bl)
{
const Fasta & ff=m_fasta_map[(*ii)]; 
const IdxTy sz=ff.size();
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & nm=ff.name(i);
const StrTy & seq=ff.seq(i);
const IdxTy len=seq.length();
Floc fl((*ii),i);
lenidx[-len].push_back(fl);
} // i 
} // ii 

return 0; 
} // IndexByLength

StrTy DumpOlaps(const OlapIdx & olaps, const IdxTy df)
{
Ss ss;
MM_LOOP(ii,olaps)
{
ss<<MMPR((*ii).first.dump(df));
const int len1=m_fasta_map[(*ii).first.f].seq((*ii).first.n).length();
MM_LOOP(jj,(*ii).second)
{
const int len2=(m_fasta_map[(*jj).f].seq((*jj).n).length());
ss<<" : "<<MMPR2((*jj).dump(df),(len1-len2));
} // jj 
ss<<CRLF;

} // ii 
//MM_ERR(MMPR(ss.str()))
return ss.str();
} // DumpOlaps
IdxTy TrimFastas(const OlapIdx & olaps, const IdxTy df)
{

MM_LOOP(ii,olaps)
{
// make al the other ones this.. 
const StrTy & x=m_fasta_map[(*ii).first.f].seq((*ii).first.n);
MM_LOOP(jj,(*ii).second)
{
m_fasta_map[(*jj).f].seq((*jj).n,x);
} // jj 

} // ii 
//MM_ERR(MMPR(ss.str()))


return 0;
} // TrimFastas

template <class LiTy> IdxTy FindOlaps( Ragged & r, const LiTy & bl, const IdxTy redlev, const IdxTy flags)
{ 
OlapIdx olaps;
return FindOlaps(olaps,bl,redlev,flags);
} // FindOlaps

template <class LiTy> IdxTy FindOlaps( OlapIdx & olaps, const LiTy & bl, const IdxTy redlev, const IdxTy flags)
{
const bool index_resolver=Bit(flags,0);
const bool trim_fasta=Bit(flags,1);
const IdxTy df=0;
const IdxTy tf=0;

LenLocs lenidx;
IndexByLength(lenidx,bl,0);
MM_LOOP(ii,bl)
{
const Fasta & ff=m_fasta_map[(*ii)]; 
const IdxTy sz=ff.size();
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & nm=ff.name(i);
const StrTy & seq=ff.seq(i);
const int len=seq.length();
Floc smaller((*ii),i);
MM_LOOP(jj,lenidx)
{
const int clen=-(*jj).first;
if (clen<len) break;
const int dlen=clen-len; //  >0
const Flocs & fl=(*jj).second;
MM_LOOP(kk,fl)
{
const StrTy & f2=(*kk).f;
const IdxTy n2=(*kk).n;
// find nm/seq inside f2[n2]a
std::cout<<". ";
//MM_ERR(MMPR4(f2,n2,m_fasta_map[f2].size(),dlen))
const StrTy & seq2=m_fasta_map[f2].seq(n2);
for(IdxTy io=0; io<=dlen; ++io)
{
//MM_ERR(MMPR4(clen,len,dlen,io))
const int sc=::strncmp((seq.c_str()),seq2.c_str()+io,len);
// overlap found 
if(sc==0)
{
Floc larger(f2,n2,io);
if (!smaller.same_seq(larger))
{
 if (index_resolver) olaps[larger].push_back(smaller);
 else olaps[smaller].push_back(larger);
} 

} // sc 
} // io 

} // kk 
std::cout<<CRLF;
} // jj 

} // i 
} // ii 
MM_ERR(MMPR(DumpOlaps(olaps,df)));
if (trim_fasta) TrimFastas(olaps,tf);
return 0; //   FindOlaps(r, bl, redlev, flags); 

} // FindOlaps





// this is invoked from the abundance table assignments which in theory
// should just come from the other biom's which have been sorted out
// and each sequence assignment can be done once from m_map entry
// for each sequence.
//MapMapTy m_fmap;

////////////////////////////////////////////// 

// flags>>8 from user 
AssTy Reconcile(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
{
// the "flag" WAS actually a pass count and not updating on second pass

const bool trim_to_match=Bit(flags,0);
const bool prefer_no_na=Bit(flags,1);

// lookup pass moved to Lookup  
//const bool update=!Bit(flags,0);
const bool update=true; // !Bit(flags,0);
const bool presume_lexi_chrono=true;
// this contaminates the index  if seqno is not found...
// The fasta or sequence may be missing.
//}

// optionally check to see if f/seqno contains a smaller sequence
// that is the root due to inconsistent trimming. 
/// in that case, the better assignment of this or the smaller needs to be considered
// or munged. 
//const auto cc=m_contains.find(Entry(f,seqno));

bool found=false;
auto ii=m_idx.find(f);
bool fasta_found=(ii!=m_idx.end());
if (!fasta_found)
{MM_ERR(" no index found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno)) return AssTy(); }
auto kk=m_fasta_map.find(f);
if (kk==m_fasta_map.end())
{MM_ERR(" no fasta found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno)) return AssTy();
}
const Fasta & ff=(*kk).second; // m_fasta_map[f];
const IdxTy ffsz=ff.size();
// apparently m_idx finds the seuqnece index for fasta f
const auto & sm=(*ii).second;
// this apparentl  returns the sequence location for the name seqno 
auto jj=sm.find(seqno);
if (jj!= sm.end()) found=true;
if (!found)
{MM_ERR(" no sequence found "<<MMPR2(found,fasta_found)<<MMPR4(f,ff.name(),ffsz,seqno))
return AssTy();
}
const IdxTy loc=(*jj).second; // m_idx[f][seqno];
//const StrTy  seq=m_fasta_map[f].seq(loc);
const StrTy  seq=(*kk).second.seq(loc);
const AssTy old=m_ass[seq];
if (!update) return old;
const StrTy fold=old.fasta();
const bool old_seq=IsSeq(old);
const bool old_tax=HasTax(old);
// so now we hae the sequence named seqno in f at location loc 
// see if it has already been read elsewhere. 
// note this is not a reference and it can modify m_map. 
// but that should be stable now after a prior loading pass
// multi is a vector of fasta/name pairs 
// this was updated during "Add" and should be stable but
// calling "[]" is not const 
// no assignment files used.... doh 
// these did not came from assignments files which should be lexi=chrono order
// taking the most recent may be best. 
//MapMapTy m_fmap;

// sequence should only appear once per fasta... 
Entry & ev=m_fmap[seq][f];
//StrTy fasta,name; std::vector<StrTy> tax; }; // _Entry
ev.fasta=f; ev.name=seqno; ev.tax=ass;
//ev.tax.fasta(f); 
// hopefully the longer sequencs come in later....
StrTy trimmed;
if (trim_to_match)
{
const IdxTy sl=seq.length();
if (sl>4){
char * pnew= new char[sl+2];
memcpy(pnew,seq.c_str()+1,sl);
pnew[sl-2]=0;
StrTy ptrim=StrTy(pnew);
const auto ff=m_fmap.find(ptrim);
if (ff!=m_fmap.end())
{
trimmed=ptrim;
MM_ERR( "trimmed hit "<<MMPR(ev.dump())<<" vs "<<(*(*ff).second.rbegin()).second.dump())
}
delete [] pnew; 
} // gt 4
} // trim_to_match
// for a trim match, the trimmed sequence entry may need to be updated
// but they need to be identified as equivalent in any case. 
// the assignments for trimmed and seq should match but the
// overall history and fasta needs some thought. 
#if 0 
for(IdxTy pass=0; pass<2; ++pass)
{
// rbegin finds the most recent if the names start lexi chrono
const auto &  ms=m_fmap[seq];
auto xx=ms.rbegin();
if (xx==ms.rend())
{
MM_ERR(" logic error for sequence "<<MMPR(ev.dump()))
}
while (xx!=ms.rend())
{
//const AssTy &  best=(*(m_fmap[seq].rbegin())).second.tax;
//const Entry&  c=(*(m_fmap[seq].rbegin())).second;
const Entry&  c=(*xx).second;
//AssTy   best=(*(m_fmap[seq].rbegin())).second.tax;
AssTy best=c.tax;
// this ddid not work  reliably
best.fasta(c.fasta);

const bool best_seq=IsSeq(best);
// this removes things like __unidenfied or __NA
//const bool best_tax=HasTax(best);
const bool best_tax=HasUnderScores(best);
bool acceptable=best.size()&&best_tax;
if ((pass==0)&&prefer_no_na&&acceptable)
{
MM_LOOP(ww,best){ const char * p=(*ww).c_str(); 
if ((*ww).length()<4) {acceptable=false; break;}
if (strcmp(p+1,"__NA")==0) { acceptable=false; break;}
if (strcmp(p+1,"__unidentified")==0) { acceptable=false; break;}
} // ww 
if (!acceptable) MM_ERR( " not accepted on first pass "<<MMPR(ev.dump()))
} // pass ==0 

if (acceptable) { m_ass[seq]=best; return best ; } 
++xx;
} // xx
} // pass
#endif
const auto &  ms=m_fmap[seq];
auto xx=ms.rbegin();
if (xx==ms.rend())
{
MM_ERR(" logic error for sequence "<<MMPR(ev.dump()))
}
const IdxTy bf=0;
while (xx!=ms.rend())
{
const Entry&  c=(*xx).second;
AssTy best=c.tax;
best.fasta(c.fasta);
IdxTy aq=AssQuality(best,0,bf);
const bool acceptable=(aq==0);
if (acceptable) { m_ass[seq]=best; return best ; } 
++xx;
} // xx
// now either return the ass or 
AssTy _old=old;
//StrTy ssrc=f+":;"+seqno;
StrTy first="X__sequences";
//bool added_new=false;

const bool use_seq=true; 
if (use_seq)
{
AssTy novel;
novel.push_back(first);

//if (old.size()>1) ssrc=ssrc+","+old[1];
StrTy ssrc;
StrTy tseq;
MM_LOOP(xx,ms)
{
const Entry&  c=(*xx).second;
const auto & ff=m_fasta_map[c.fasta];
auto ii=m_idx.find(c.fasta);
const auto & sm=(*ii).second;
auto jj=sm.find(c.name);
const IdxTy loc=(*jj).second; // m_idx[f][seqno];
tseq=ff.seq(loc);
if (ssrc.length()) ssrc=ssrc+",";
// wth? 
 //ssrc+=c.fasta+":;"+seqno;
 ssrc+=c.fasta+":;"+c.name;
//ev.fasta=f; ev.name=seqno; ev.tax=ass;
//AssTy best=c.tax;
//best.fasta(c.fasta);
} // xx
novel.push_back(ssrc);
m_ass[seq]=novel;
m_ass[seq].fasta(StrTy());
return novel;
}
//} // found really dum now 
m_ass[seq]=ass;
m_ass[seq].fasta(f);
return ass;

} // Reconcile 

#if 0

const IdxTy szass=ass.size();
bool use_seq=false;
if (szass==0) use_seq=true;
else {
//StrTy  x=ass[0]; 
StrTy  x=ass[ass.size()-1]; 
if(strncmp(x.c_str(),"X__",3)==0) use_seq=true; 
if (x.length()>3) if(x.c_str()[2]=='_') x=x.substr(3);
if (x=="blank") use_seq=true;
if (x=="None") use_seq=true;
if (x=="") use_seq=true;
if (x=="-") use_seq=true;
if (x=="unidentified") use_seq=true;
if (x=="s__unidentified") use_seq=true;
if (x=="NA") use_seq=true;
}
if (use_seq)
{
//const StrTy val=f+" "+seqno;
//StrTy val=StrTy("X__")+f+"_"+seqno;
//if (old_seq) val=val+"_"+old[old.size()-1];
//if (old_seq) if(val!=old[old.size()-1]) val=val+"_"+old[old.size()-1];
// one separator or the other lol, easier than escaping etc lol 
AssTy _old=old;
const StrTy ssrc=f+":;"+seqno;
StrTy first="X__sequences";
bool added_new=false;
AssTy novel;
const IdxTy oldsz=_old.size();
if (oldsz)
{
if (_old[0]==first)
{
if (oldsz>1){  _old[1]+=","+ssrc; added_new=true; } 
} // first

} // old 
//MM_LOOP(ii,old) novel.push_back(*ii);
// great so now this is a bunch of refs
if (!added_new)
{
novel.push_back(first);
novel.push_back(ssrc);
// TODO FIXME this could grow beyond the expected hierarchy size
// and nt sure if code is good lol 
} // added_new 
MM_LOOP(ii,_old) novel.push_back(*ii);
// this seems great but is hard to read
/*
novel.push_back(seqno);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
*/

m_ass[seq]=novel;
return novel;
}
//} // found really dum now 
m_ass[seq]=ass;
m_ass[seq].fasta(f);
return ass;

//}
//else MM_ERR(" no sequence found "<<MMPR(loc)<<MMPR4(f,ff.name(),ffsz,seqno))
//return AssTy();



} // Reconcile

#endif


// to save pointless tests, stop at "threshold" 
IdxTy AssQuality(const AssTy & ass, const IdxTy threshold, const IdxTy flags)
{
IdxTy rc=0;
const IdxTy sz=ass.size();
if ( sz==0) return ~0;
MM_LOOP(ii,ass)
{
StrTy  x=(*ii); 
if(strncmp(x.c_str(),"X__",3)==0) return 1; 
if (x.length()>3) if(x.c_str()[2]=='_') x=x.substr(3);
if (x=="blank") return 2;
if (x=="None") return 2;
if (x=="") return 2;
if (x=="-") return 2;
if (x=="unidentified") return 2;
//if (x=="bacterium") return 2;
if (x=="s__unidentified") return 2;
if (x=="NA") return 2;
if (x=="uous_taxa") return 2;

} // ii 

return rc;
} // AssQuality 



////////////////////////////////////////////////////////////j
////////////////////////////////////////////////////////////j



////////////////////////////////////////////// 
AssTy ReconcileOld(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
{
const bool update=!Bit(flags,0);
const bool presume_lexi_chrono=true;
// this contaminates the index  if seqno is not found...
// The fasta or sequence may be missing.
//}

// optionally check to see if f/seqno contains a smaller sequence
// that is the root due to inconsistent trimming. 
/// in that case, the better assignment of this or the smaller needs to be considered
// or munged. 
//const auto cc=m_contains.find(Entry(f,seqno));

bool found=false;
auto ii=m_idx.find(f);
bool fasta_found=(ii!=m_idx.end());
if (!fasta_found)
{MM_ERR(" no sequence found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno))
return AssTy();
}
auto kk=m_fasta_map.find(f);
if (kk==m_fasta_map.end())
{MM_ERR(" no fasta_map found "<<MMPR2(found,fasta_found)<<MMPR2(f,seqno))
return AssTy();
}
const Fasta & ff=(*kk).second; // m_fasta_map[f];
const IdxTy ffsz=ff.size();
// apparently m_idx finds the seuqnece index for fasta f
const auto & sm=(*ii).second;
// this apparentl  returns the sequence location for the name seqno 
auto jj=sm.find(seqno);
if (jj!= sm.end()) found=true;
//else MM_ERR(" no sequence found "<<MMPR(loc)<<MMPR4(f,ff.name(),ffsz,seqno))
//return AssTy();
if (!found)
{MM_ERR(" no sequence found "<<MMPR2(found,fasta_found)<<MMPR4(f,ff.name(),ffsz,seqno))
return AssTy();
}
const IdxTy loc=(*jj).second; // m_idx[f][seqno];
const StrTy  seq=m_fasta_map[f].seq(loc);
// so now we hae the sequence named seqno in f at location loc 
// see if it has already been read elsewhere. 
// note this is not a reference and it can modify m_map. 
// but that should be stable now after a prior loading pass
// multi is a vector of fasta/name pairs 
// this was updated during "Add" and should be stable but
// calling "[]" is not const 
// no assignment files used.... doh 
// these did not came from assignments files which should be lexi=chrono order
// taking the most recent may be best. 
const EntryVector multi=m_map[seq];
const IdxTy multiplicity=multi.size();
MM_ERR(MMPR4(f,seqno,loc,multiplicity))
if (multiplicity==0)
{ MM_ERR(" no entries in m_map "<<MMPR(seq)) }
// find the current assignment for this sequence 
const AssTy old=m_ass[seq];
if (old.size()==0)
{
MM_ERR(" existing assignment is zero length "<<MMPR4(update,f,seqno,multiplicity))
}
// if not updating, this should be stable now. 
if (!update) return m_ass[seq];
// otherwise pick the best and updated. 
const bool old_seq=IsSeq(old);
const bool new_seq=IsSeq(ass);
const bool old_has_tax=HasTax(old)&&!old_seq;
const bool new_has_tax=HasTax(ass)&&!new_seq;
if (old_has_tax&&new_has_tax&&(old!=ass))
{
MM_ERR(" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA ")
// the tax assignments are not in the fasta files which just 
// reqlates sequence to numer.
// in fact, the tax assignment files are never referenced in the command line
MM_ERR(MMPR4(f,seqno,loc,multiplicity))
Ss ss; 
for(IdxTy i=0; i<multi.size(); ++i )  { ss<<MMPR3(i,multi[i].fasta,multi[i].name)<<" "; }
MM_ERR(ss.str())
MM_ERR(" reconciling "<<MMPR2(old[old.size()-1],ass[ass.size()-1]))
{ Ss ss; MM_LOOP(ii,old) { ss<<(*ii)<<" "; }MM_ERR(" recon oold "<<MMPR(ss.str()))  }
{ Ss ss; MM_LOOP(ii,ass) { ss<<(*ii)<<" "; }MM_ERR(" recon ass "<<MMPR(ss.str()))  }

MM_ERR(" ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ ")
}
// the important thing was to unit unknowns 
// however, this fails as "Add" comes from the fasta with zero tax
// info. This needs to use assignment files 
if(presume_lexi_chrono)
{
StrTy recent="";
Entry ex;
bool sfound=false;
MM_LOOP(ii,multi)
{
const Entry & e=(*ii);
//StrTy fasta,name; std::vector<StrTy> tax;
// again, presumed chronological paths taking the most recent
if (e.fasta>recent){ recent=e.fasta; ex=e; }
if (e.tax==ass) sfound=true;
} // multi
if (!sfound)
{
MM_ERR(" the included assignment was not found in Added assignments ")
}

if (ex.tax.size()) { m_ass[seq]=ex.tax; return m_ass[seq]; }
else
{
MM_ERR(" no good tax "<<MMPR4(recent,ex.tax.size(),f,seqno))
} // size 

} // presume_lexi_chrono
// keep the old one 
if (old_has_tax && ! new_has_tax) return m_ass[seq]; 
// this implies new has tax too now... 
if ( old_has_tax) { if ( old.size()>ass.size()) return m_ass[seq]; } 
// new_has_tax but either old does not have tax or is shorter 
// fall through is ok 
//if ( !old_has_tax) { m_ass[seq]=ass; return ass; } 
// now old has tax but it is shorter...

const IdxTy szass=ass.size();
bool use_seq=false;
if (szass==0) use_seq=true;
else {
//StrTy  x=ass[0]; 
StrTy  x=ass[ass.size()-1]; 
if (x.length()>3) if(x.c_str()[2]=='_') x=x.substr(3);
if (x=="blank") use_seq=true;
if (x=="None") use_seq=true;
if (x=="") use_seq=true;
if (x=="-") use_seq=true;
if (x=="unidentified") use_seq=true;
if (x=="s__unidentified") use_seq=true;
if (x=="NA") use_seq=true;
}
if (use_seq)
{
//const StrTy val=f+" "+seqno;
StrTy val=StrTy("X__")+f+"_"+seqno;
//if (old_seq) val=val+"_"+old[old.size()-1];
if (old_seq) if(val!=old[old.size()-1]) val=val+"_"+old[old.size()-1];
AssTy novel;
novel.push_back(seqno);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
novel.push_back(val);
m_ass[seq]=novel;
return novel;
}
//} // found really dum now 
m_ass[seq]=ass;
return ass;

//}
//else MM_ERR(" no sequence found "<<MMPR(loc)<<MMPR4(f,ff.name(),ffsz,seqno))
//return AssTy();



} // ReconcileOld
////////////////////////////////////////////////////////////j



bool HasTax(const AssTy & old)
{
bool old_has_tax=false;
if ( old.size()) 
{
const StrTy last= old[old.size()-1];
if ( last.length()>3) if (last.c_str()[1]=='_') old_has_tax=true;
if (old_has_tax) if (strcmp(last.c_str()+1,"__NA")==0) old_has_tax=false;
if (old_has_tax) if (strcmp(last.c_str()+1,"__unidentified")==0) old_has_tax=false;
if (old_has_tax) if (strcmp(last.c_str()+1,"__None")==0) old_has_tax=false;

} // old.size
// remove x__NA from consideration.. 
if (old_has_tax) //if (old.size()>3)
{
//const char * p=o


}

return old_has_tax;
} // HasTax


bool HasUnderScores(const AssTy & old)
{
bool old_has_tax=false;
if ( old.size()) 
{
old_has_tax=true;
MM_LOOP(ii,old)
{
const StrTy & x=(*ii);
const char * p=x.c_str();
if (x.length()<4) return false;
if (p[0]=='X') return false;
if (p[0]=='x') return false;
if (p[1]!='_') return false;
if (p[2]!='_') return false;
} // ii 
//const StrTy last= old[old.size()-1];
//if ( last.length()>3) if (last.c_str()[1]=='_') old_has_tax=true;
//if (old_has_tax) if (strcmp(last.c_str()+1,"__NA")==0) old_has_tax=false;
//if (old_has_tax) if (strcmp(last.c_str()+1,"__unidentified")==0) old_has_tax=false;
//if (old_has_tax) if (strcmp(last.c_str()+1,"__None")==0) old_has_tax=false;

} // old.size
return old_has_tax;
} // HasUnderScores





bool IsSeq(const AssTy & old)
{
bool old_has_tax=false;
if ( old.size()) 
{
const StrTy last= old[old.size()-1];
if ( last.c_str()[0]=='X') if (last.c_str()[1]=='_') old_has_tax=true;
} // old.size
return old_has_tax;
} // IsSeq 

IdxTy LoadAss(const StrTy & nm, const StrTy & fn, const IdxTy flags)
{ 
std::ifstream fis(fn);
auto & at=m_lut[nm];
while ( fis.good()&&!fis.eof())
{
StrTy s="";
AssTy ass;
StrTy seq;
std::getline(fis,s);
const IdxTy slen=s.length();
char c[slen+1];
memcpy(c,s.c_str(),slen+1);
char * p =c;
char * pstart=p;
// there is a tab at the beginning but then more crap later
//possibly 
/*
==> ./data/2018-03-09/tax_assignments.txt <==
seq1	k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__perfoetens	1	1
seq2	k__Bacteria;p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Veillonellaceae;g__Megamonas;s__funiformis	1	1
*/

while (*p)
{
// fortunately since there is no pstart this continues to 
// work by accident lol. But it fails to get
// the s__ term....  
//if ((*p=='\t')||(*p==' '))
if ((*p=='\t')&&(ass.size()==0))
{ *p=0; seq=StrTy(c); ++p; pstart=p; continue;  }
if ( *p==';')
{ *p=0; ass.push_back(StrTy(pstart)); ++p; pstart=p; continue;  }
if ( *p=='\t')
{ *p=0; ass.push_back(StrTy(pstart)); ++p; pstart=p; break;  }
++p;
} // p
if (pstart!=p)  ass.push_back(StrTy(pstart));
at[seq]=ass;
} // fis
MM_ERR(MMPR4(__FUNCTION__,nm,fn,at.size()))
return 0;
} // LoadAss 


// MEMBERS
RaggedMap m_ragged_map;
FastaMap m_fasta_map;
IndexTy m_idx;
AssMapTy m_ass;
MapTy m_map;
MapMapTy m_fmap;
AssMapMap m_lut; // assignments for each fasta of same name... 
}; // mjm_sequence_reconcile

//////////////////////////////////////////////

template <class Tr>
class mjm_sequence_reconcile_map : public std::map<typename Tr::StrTy, mjm_sequence_reconcile< Tr > >  
{
 typedef mjm_sequence_reconcile_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_sequence_reconcile< Tr> >   Super;
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
mjm_sequence_reconcile_map() {}
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

}; // mjm_sequence_reconcile_map




////////////////////////////////////////////
#ifdef  TEST_MJM_SEQUENCE_RECONCILE
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
typedef tester_< mjm_sequence_reconcile <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_SEQUENCE_RECONCILE "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_sequence_reconcile<Tr>  Myt;
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
if (cmd=="dump-eq") { MM_ERR(x.dump_eq()) }
if (cmd=="dump-ass") { MM_ERR(x.dump_ass(atoi(cip.p1.c_str()))) }
if (cmd=="dump-ass-fasta") { 
MM_ERR(MMPR3(cmd,cip.p1,cip.p2))
StrTy s= x.dump_ass_fasta(atoi(cip.p2.c_str()));
std::ofstream fox(cip.p1);
fox<<s;

 }
else if (cmd=="load") { x.load(cip.p1,cip.p2,atoi(cip.wif(3).c_str())); }
else if (cmd=="load-ass") { x.load_assignments(cip.p1,cip.p2,atoi(cip.wif(3).c_str())); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_SEQUENCE_RECONCILE_H__ 
