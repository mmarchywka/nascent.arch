#ifndef MJM_TAXON_TOOLS_H__
#define MJM_TAXON_TOOLS_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
#include "mjm_generic_iterators.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_block_matrix.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_collections.h"

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

namespace taxon_tools_traits
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

class mjm_taxon_tools 
{
typedef mjm_taxon_tools Myt;
protected:
typedef taxon_tools_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;


public:
typedef std::map<StrTy,IdxTy> ignore_t;
typedef std::vector<IdxTy> TaxVec;
typedef std::vector<StrTy> TaxNameVec;
typedef TaxVec tax_t;
typedef TaxNameVec tax_name_t;
class tax_split { public: tax_split() {} tax_split(const tax_name_t & _t, const D & f) : t(_t),frac(f) {} tax_name_t t; D frac; };
typedef std::vector<tax_t> Taxes;
typedef Taxes taxes_t;

typedef std::vector<StrTy> Words;
typedef std::vector<tax_split> tax_splits_t;

mjm_taxon_tools(): m_dmel(0),m_info_flags(0),m_unknown_conformed("UNK") {Init(); }
mjm_taxon_tools(const IdxTy flags ): m_dmel(0),m_info_flags(flags),
m_unknown_conformed("UNK") {Init(); }
virtual ~mjm_taxon_tools() {}
// yes, I know... 
static bool is_lower(const char c) { return (c>='a')&&(c<='z'); } 
Words split(const StrTy & w) const
{
Words v;
const char * c=w.c_str();
const IdxTy sz=w.length()-1;
char d[sz+2];
::memcpy(d,c,sz+1);
d[sz+1]=0;
IdxTy sp=0;
for(IdxTy i=1; i<sz; ++i)
{
if (d[i]=='-') if (is_lower(d[i-1])&&is_lower(d[i+1]))
{
d[i]=0;
v.push_back(d+sp);
sp=i+1;
}
} // i 
if (sp<(sz+1)) v.push_back(d+sp);
return v;
}

tax_splits_t split(const tax_name_t tv) const
{
tax_splits_t t;
//tax_t ts;

MM_LOOP(ii,tv)
{
// FIXME realistically split should only occur ate one level.. 
const Words w=split((*ii));
const IdxTy alts=w.size(); 
if (alts==0) 
{
MM_ERR(" zero size "<<MMPR2(tv.size(),(*ii))) 
}
if (t.size()==0) t.push_back(tax_split());
const IdxTy tsz=t.size();
for(IdxTy j=1; j<alts; ++j)
for (IdxTy i=0; i<tsz; ++i) { t.push_back(t[i]); }

for(IdxTy j=0; j<alts; ++j)
for (IdxTy i=0; i<tsz; ++i) 
{
t[j*tsz+i].t.push_back(w[j]); }

} // i 

MM_LOOP(ii,t) { (*ii).frac=1.0/t.size(); } 
return t;
}


const ignore_t & ignores() const { return m_ignored; } 
const ignore_t & default_ignores() const 
{
static ignore_t ignored;
static bool init=false;
if (!init)
{
ignored[StrTy("uncultured bacterium")]=1;
ignored[StrTy("unculturedbacterium")]=1;
ignored[StrTy("unidentified")]=1;
ignored[StrTy("unclassified")]=1;
ignored[StrTy("blank")]=1;
// these should hve been removed before 
ignored[StrTy("g__")]=1;
ignored[StrTy("s__")]=1;
ignored[StrTy("o__")]=1;
ignored[StrTy("f__")]=1;
ignored[StrTy("c__")]=1;
// place holder for non-terminal counts
ignored[StrTy("MISSING")]=1;
ignored[StrTy("na")]=1;
init=true;
}

return ignored;
}

template <class Tt>
StrTy informative_name(const TaxVec & tv, const Tt & tt) const
{
const ignore_t & ignored=ignores();
StrTy taxon="";
const StrTy sep=" ";
for (int  i=0; i<tv.size(); ++i)
{
const StrTy & w=tt.node_name(tv[i]);
if (ignored.find(w)==ignored.end())
    {if (taxon.length()!=0) { taxon=w+sep+taxon; break; } taxon=w; }

}
return taxon;
}



tax_name_t  informative_vector(const tax_name_t & tv) const
{
tax_name_t taxon;
const ignore_t & ignored=ignores();
const StrTy sep=" ";
for (int  i=0; i<tv.size(); ++i)
{
const StrTy & w=(tv[i]);
if (ignored.find(w)==ignored.end()) taxon.push_back(w);
else { break; }  // these should not be LINEAGES, bacteria first entery all after uninforative bad. 
}
return taxon;
}

StrTy  informative_n_name(const tax_name_t & tv, const int  n) const
{
StrTy taxon;
const ignore_t & ignored=ignores();
StrTy sep="";
IdxTy j=0; 
for (int  i=0; i<tv.size(); ++i)
{
const StrTy & w=(tv[i]);
if (ignored.find(w)==ignored.end()){ ++j;  const bool ok=((n<0)&&((tv.size()-i)==(-n-1)))||(n>=0);

//if (n<0) { MM_ERR(" ASSSFCUUUUUUUUK "<<MMPR4(w,i,n,j)<<MMPR4(tv.size(),(tv.size()-i),(-n-1),ok)) } 

  if (ok) { taxon=w+sep+taxon; sep=" ";} }
if (n>0) if (j>=n) break; 
//else { break; }  // these should not be LINEAGES, bacteria first entery all after uninforative bad. 
}
return taxon;
}




protected:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }

void Init()
{
 m_ignored = default_ignores();
}



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
IdxTy m_info_flags;
StrTy m_unknown_conformed;
ignore_t m_ignored;
};  // mjm class_


#endif

