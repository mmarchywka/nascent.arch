#ifndef MJM_ZYMO_H__
#define MJM_ZYMO_H__

#include "mjm_globals.h"
#include "mjm_pheno_notes.h"
#include "mjm_taxon_tools.h"
#include "mjm_tree_viz.h"
#include "mjm_data_model_error_log.h"
#include "mjm_zymo_composition_table.h"
#include "mjm_column_plot.h"
#include "mjm_kv_vector.h"
#include "mjm_fasta_ii.h"
#include "mjm_sequence_reconcile.h"

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
// get the taxa names conforming 
#include "mjm_strings.h"

// use the ncbi taxa DB dump 
// ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt

#include "mjm_collections.h"
#include "mjm_tokenized_collections.h"

//  2409  echo -e "read-comp x1 data/2022-02-18/zr5958_l7.txt 8\nread-comp x2 junk/yyy 8\nsync-comp x1\nsync-comp x2\nmerge-comp x2 x1\nsave-comp-notes-ssv x2 foo"| ./mjm_zymo.out 2>fack


//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
#include "mjm_biom_hdf5.h"

#include "mjm_biom_to_rag.h"
#include "mjm_rag_reduce.h"

#ifdef PYTHON_BOOST_BUILD
#include <boost/python.hpp>
#endif




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

TODO FIXME
*/

/*
Adding ncbi tax tree,
2584  echo -e "load-tax-tree /home/marchywka/junk/tax\nsave-tax " | ./mjm_trees_and_tables.out    > xxx
 2585  mv xxx tax-info 
 2586  echo -e "load-tax tax-info \nsave-tax " | ./mjm_trees_and_tables.out    > xxx
 2606  run_trees_and_tables -opt -compile 
 2610  echo -e "load-tax tax-info \nbest-taxon \"cellular organisms\" bacteria proteobacteria alphaproteobacteria  " | ./mjm_trees_and_tables.out    > xxx
 2616  echo -e "load-tax tax-info \nbest-indexed-taxon \"cellular organisms\" bacteria proteobacteria alphaproteobacteria  " | ./mjm_trees_and_tables.out    > xxx
 2617  echo -e "load-tax tax-info \nbest-indexed-taxon helicobacter  " | ./mjm_trees_and_tables.out    > xxx
 2618  echo -e "load-tax tax-info \nbest-indexed-taxon pylori  " | ./mjm_trees_and_tables.out    > xxx
 2619  echo -e "load-tax tax-info \nbest-indexed-taxon asdfasdf  " | ./mjm_trees_and_tables.out    > xxx
=====================================================
 


otf_angles, note that conform_table.txt was added 

4777  run_zymo -opt   -compile 
 4778  ./mjm_zymo.out -source z8.txt -quit  > zzz8
 4779  grep unclassif taxa
 4780  history | grep zz8
 4781  echo dump-otf-class-stats zzz8 | ./mjm_zymo.out  2>xxx | mjm eq | awk '{print "dummy group "$1" "$6" "$8" "$2" "$NF}' > ppp
 4782  history | grep ppp
 4783  echo otf-angles ppp | ./mjm_zymo.out  2>xxx  > fff

==========================================================================================
2018-01-23 make angle plots between all samples from differen tsources,

4584  ./mjm_zymo.out -source z8.txt -quit  > zzz8
 4585  history | grep zzz8
 4586  echo dump-otf-class-stats zzz8 | ./mjm_zymo.out  2>xxx | mjm eq | awk '{print "dummy group "$1" "$6" "$8" "$2" "$NF}' > ppp
 4587  history | grep ppp
 4588  echo otf-angles ppp | ./mjm_zymo.out  2>xxx  > fff

cat z8.txt | grep -v "^#"
read-catagory-map phyla phyla.txt
read-catagory-map genera genera_plus.txt
dump-catagory-map phyla
set-param catagory_map genera 
set-param use_sample_filter 0 
parse-biom-json emp /home/marchywka/junk/emp.json   
copy-biom emp foo 
clear-biom emp  

parse-biom-json zymo /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 
copy-biom zymo zymo 
clear-biom zymo   
set-param otu_format 0x000000
classify-sample-group foo /home/marchywka/junk/sample_envo
update-order-format foo 0x0C 
update-order-format zymo 0x0C 
sort-order
dump-group-n foo  0x01E    
dump-group-n zymo 0x00E    
=================================================================================================


 run_zymo   -compile 
g++ -DTEST_ZYMO__ -std=gnu++11 -gdwarf-3 -O0 -I.. -I../json/rapidjson-master/include/ -Wall -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_zymo.h -o mjm_zymo.out

*/





////////////////////////////////////////////////////////////////

class zymo_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
zymo_params( const StrTy & nm) : Super(nm) {}
zymo_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
//	const StrTy sample_name_class_map=m_flp.adhoc_classify_name("xlate");
StrTy adhoc_classify_name() const { return m_map.get_string("adhoc_classify_name","xlate"); }
StrTy groups_to_classify() const { return m_map.get_string("groups_to_classify","grclass"); }
StrTy plot_hash_name() const { return m_map.get_string("plot_hash_name","plot-hash"); }
StrTy otu_lumps_name() const { return m_map.get_string("otu_lumps_name","otu-hash"); }
IdxTy otu_depth_limit() const { return m_map.get_uint("otu_depth_limit",0); } // // 100;
StrTy allowed_otu_name() const { return m_map.get_string("allowed_otu_name","otu-names"); }

bool crp_have_header() const { return m_map.get_bool("crp_have_header",false); }
bool crp_include_header() const { return m_map.get_bool("crp_include_header",false); }
bool crp_fix_offset() const { return m_map.get_bool("crp_fix_offset",false); }
bool crp_add_abundances() const { return m_map.get_bool("crp_add_abundances",!false); }
bool crp_add_identicals() const { return m_map.get_bool("crp_add_identicals",false); }
IdxTy crp_saf() const { return m_map.get_int("crp_saf",3); } // // 100;
int crp_osize() const { return m_map.get_int("crp_osize",5); } // // 100;
int crp_n1() const { return m_map.get_int("crp_n1",1); } // // 100;
int crp_n2() const { return m_map.get_int("crp_n2",2); } // // 100;
int crp_ngroup() const { return m_map.get_int("crp_ngroup",5); } // // 100;
int crp_p() const { return m_map.get_int("crp_p",-1); } // // 100;



// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands"<<log_commands()<<sep;
ss<<"exit_on_err"<<exit_on_err()<<sep;
ss<<"otu_format="<<otu_format()<<sep;
ss<<"catagory_map="<<catagory_map()<<sep;
ss<<"sample_filter="<<sample_filter()<<sep;
ss<<"use_sample_filter="<<use_sample_filter()<<sep;
ss<<"adhoc_classify_name="<<adhoc_classify_name()<<sep;
ss<<"groups_to_classify="<<groups_to_classify()<<sep;
ss<<"plot_hash_name="<<plot_hash_name()<<sep;
ss<<"otu_lumps_name="<<otu_lumps_name()<<sep;
ss<<"otu_depth_limit="<<otu_depth_limit()<<sep;
ss<<"allowed_otu_name="<<allowed_otu_name()<<sep;

ss<<crp_have_header();
ss<<"crp_include_header="<<crp_include_header();
ss<<crp_fix_offset(); 
ss<<crp_add_abundances(); 
ss<<"crp_add_identicals="<<crp_add_identicals(); 
ss<<"crp_saf="<<crp_saf();
ss<<crp_osize();
ss<<crp_n1();
ss<<crp_n2();
ss<<crp_ngroup();
ss<<crp_p() ;




return ss.str();
}


}; // zymo_params


// reserved words or specific things with meanings different
// from canonical names .



namespace zymo_traits
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
}; // zymo_traits

class sparse_connectivity
{
typedef sparse_connectivity Myt;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef  Tr::Dmel  Dmel; 

typedef string_tokenizer St;


//typedef std::vector<IdxTy> Dvec;
typedef std::map<IdxTy,D> Dvec;
typedef std::map<IdxTy, Dvec> Cmap;

public:
void connect(const StrTy&  r, const StrTy & c , const D & v) 
{ m_connections[m_xlate(r)][m_xlate(c)]+=v; }
void connect(const IdxTy&  r, const IdxTy & c , const D & v) 
{ m_connections[r][c]+=v; }

void normalize(const IdxTy samples)
{
// rows, first index is OTU and second one is sample 
const IdxTy sz=m_connections.size(); // this should be the OTU count
IdxTy i=0; 
std::vector<D> totals(samples);
MM_LOOP(ii,m_connections)
{
const auto & c=(*ii).second;
MM_LOOP(jj,c)
{
totals[(*jj).first]+=(*jj).second;
} // jj
++i;
}
MM_LOOP(ii,m_connections)
{
auto & c=(*ii).second;
MM_LOOP(jj,c)
{
const D sum=totals[(*jj).first];
if (sum!=0)  { (*jj).second/=sum; } 
} // jj 
} // ii

}


void angles(const IdxTy sz)
{
MyBlock m(sz,sz);
MM_LOOP(ii,m_connections)
{
const auto & c=(*ii).second;
MM_LOOP(jj,c)
{
const IdxTy i=(*jj).first;
const D vj=(*jj).second;
MM_LOOP(kk,c)
{
const IdxTy j=(*kk).first;
const D vk=(*kk).second;
m(i,j)+=vk*vj;
} // kk 
} // jj

} // ii 
const StrTy sep=" ";
//std::cout<<m.to_string();
std::ostream & os =std::cout;
for (IdxTy i=0; i<sz; ++i)
{
for (IdxTy j=0; j<sz; ++j)
{
//os<<m(i,j)<<sep;
const D norm=m(i,i)*m(j,j);
if (norm>0) os<<i<<sep<<j<<sep<<m(i,j)/sqrt(m(i,i)*m(j,j))<<CRLF;
else os<<i<<sep<<j<<sep<<m(i,j)<<CRLF;
}
//os<<CRLF; 
}

} // angles

private:

St m_xlate;
Cmap m_connections;


}; // sparse_connectivity


/*
The zymo taxa info looks like this, conforming means remove level pfx and make sure
it is in right place. Flag the XIII thing for now. 
zymo seq111 k__Bacteria p__Actinobacteria c__Actinobacteria o__Corynebacteriales f__Mycobacteriaceae g__Mycobacterium s__aurum
-austroafricanum-chubuense-gilvum-iranicum-parafortuitum-vanbaalenii
zymo seq112 k__Bacteria p__Firmicutes c__Clostridia o__Clostridiales f__Family XIII g__Family-[Eubacterium] s__unidentified
zymo seq113 k__Bacteria p__Actinobacteria c__Actinobacteria o__Corynebacteriales f__Corynebacteriaceae g__Corynebacterium s__m
ucifaciens-pilbarense-ureicelerivorans

*/
// needless to say, this needs a dmel object wtf 
class otu_conform
{
typedef otu_conform Myt;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef  Tr::Dmel  Dmel; 
enum{ BAD=~0,ZYMO_TAXSZ=6};
typedef std::map<char,IdxTy> TaxMap;
typedef std::map<char,IdxTy> BadCharMap;
// may have levels 
typedef std::map<StrTy, StrTy > XlateMap;
typedef std::map<StrTy, StrTy > ExpMap;
typedef std::map<IdxTy,IdxTy > TaxaHist;


public:

typedef std::map<StrTy, IdxTy > CanonMap;
typedef std::map< IdxTy,CanonMap > CanonLevelMap;

otu_conform():m_dmel(0) ,m_serial(0),m_canon(0)
{
Init();

}
void canon_map( CanonLevelMap * m ) { m_canon=m; } 
void Init()
{
	char c[ZYMO_TAXSZ+1]; 
	ZymoTaxc(&c[0]);
	LoadTaxMap(m_m,& c[0]);
	MM_ERR(" inited otu_conform") 
}
// iterate over something convertible into a StrTy and put into
// subscriptable container v
template <class Tv, class Ti, class Te >
void zymo_to_vector( Tv & v, Ti  ii, const Te & ee, const IdxTy flags=0)
{
++m_serial;
TaxMap & m=m_m;
IdxTy level=0;
const bool track_used=((flags&1)!=0);
const bool track_hier_size=((flags&2)!=0);
while ( ii!=ee)
{
	const auto & s=(*ii);
	if (s.length()<4)
	{
		DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false);
		++ii;
		++level;
		m_exp[s]="length"; 
		continue;
	}
	const  IdxTy olevel=v.size();
	IdxTy pfxloc=3;
	const IdxTy slevel=Ilevel(pfxloc,s,m,olevel);
	if (slevel!=olevel) // (slevel==BAD)
	{
		DMel(m_ss<<MM_STR_LOC<<MMPR4(level,olevel,slevel,s));
		m_exp[s]="level"; 
		++ii;
++level;
	if ( olevel<slevel)	continue; 
	while ( slevel>v.size()) {v.push_back("missing"); } 
	} 
	StrTy sconformed=Scanforbads(mjm_strings::to_lower(s.substr(pfxloc)),level,s);
	// check canon and xlation table
	typedef XlateMap::const_iterator Xi;
	Xi xi=m_xlate.find(sconformed);
	// TODO m_xlate is NOT the fing canon lut doh 
	if (xi!=m_xlate.end()) sconformed=(*xi).second; 
	else 
	{
//		DMel(m_ss<<MM_STR_LOC<<" canon miss "<<MMPR4(sconformed,olevel,slevel,s));
//		m_exp[s]="missing"; 
	}

	if (!IsCanon(sconformed,slevel))
	{
		DMel(m_ss<<MM_STR_LOC<<" canon fail "<<MMPR4(sconformed,olevel,slevel,s));
		m_exp[s]="canon"; 
	}
if (track_used)	++m_used[v.size()][sconformed];
	v.push_back(sconformed);
++level;
++ii;
} // ii 
if ( track_hier_size) ++m_taxa_hist[v.size()];
} // zymo_to_vector
IdxTy  unique_exceptions() const { return m_exp.size(); } 
void exception_dump ( OsTy & os, const StrTy & label="exceptions" ) const
{
MM_LOOP(ii,m_exp)
{
os<<MM_STR_LOC<<" "<<label<<" "<<(*ii).second<<" "<<(*ii).first<<CRLF;
} // ii

} // exception_dumo

void used_dump ( OsTy & os, const StrTy & label="usedcanon" ) const
{
MM_LOOP(ii,m_used)
{
const IdxTy level=(*ii).first;
MM_LOOP(jj,(*ii).second)
{
os<<MM_STR_LOC<<" "<<label<<" "<<level<<" "<<(*jj).first<<" "<<(*jj).second<<CRLF;
} // jj 
} // ii

} // exception_dumo


private:
void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss, const bool print=true)
{
	Ss& ss=dynamic_cast<Ss& >(_ss);
	if (print) { std::cerr<<ss.str()<<CRLF; } 
	ss.str(StrTy(""));

}
bool IsCanon(const StrTy & sconformed,const IdxTy & slevel)
{
if (m_canon==0) return true;
if (m_canon->size()==0) return true;
MM_ONCE(" doing at least one canon check "<<MMPR2(sconformed,slevel),) 
CanonMap & m=(*m_canon)[slevel];
return (m.find(sconformed)!=m.end()); 
}

StrTy Scanforbads(const StrTy & s, const IdxTy level, const StrTy & sorg)
{
Ss ss,bad;
const char * sc=s.c_str();
const IdxTy sz=s.length();
for (IdxTy i=0; i<sz; ++i)
{
 char  c=sc[i];
if (c==' ') c='_';
else if (c=='\t') c='_';

if ((c<'a')||(c>'z')) if ((c<'0')||(c>'9')) if (c!='_')
{
		bad<<c;
		continue;
}

ss<<c;
} // i 
if (bad.str().length()!=0) { 
		DMel(m_ss<<MM_STR_LOC<< " bad chars "<<MMPR3(s,level, bad.str()),false);
		m_exp[sorg]="badchars"; 
} 
return ss.str();
}

IdxTy Ilevel( IdxTy & pfxloc, const StrTy & s, const TaxMap & m, const IdxTy & def) 
{
const char c=s.c_str()[0];
const char c1=s.c_str()[1];
const char c2=s.c_str()[2];
if ( c1!='_'){ pfxloc=0;  return def; } 
if ( c2!='_'){ pfxloc=0;  return def; }
auto ii= m.find(c);
if ( ii==m.end()) return BAD;
return (*ii).second;

}
void ZymoTaxc(char * c)
{
 IdxTy p=0; 
c[p]='k';  ++p; 
c[p]='p';  ++p; 
c[p]='c';  ++p; 
c[p]='o';  ++p; 
c[p]='f';  ++p; 
c[p]='g';  ++p; 
c[p]='s';  ++p;
c[p]=0;  ++p;
}
void LoadTaxMap( TaxMap & m, const char * c)
{ const char * p=c; IdxTy i=0; while (*(p+i)!=0) {m[*(p+i)]=i; ++i; } }
void LoadBadChars()
{
IdxTy pos=1;
m_bad[' ']=pos; ++pos;
m_bad['\t']=pos; ++pos;
m_bad['"']=pos; ++pos;
m_bad['[']=pos; ++pos;
m_bad[']']=pos; ++pos;


}
Dmel * m_dmel;
IdxTy m_serial;
Ss m_ss;
TaxMap m_m;
BadCharMap m_bad;
XlateMap m_xlate;
CanonLevelMap*  m_canon;
CanonLevelMap m_used;
ExpMap m_exp; // this belongs in Dmel ... 
TaxaHist m_taxa_hist;
}; // otu_conform


/*
cyanobacteria | 1 | _;_bacteria;_terrabacteria_group;_cyanobacteria/melainabacteria_group AKA _cyanophyceae _cyanophycota _cya
nophyta _oxygenic_photosynthetic_bacteria _oxyphotobacteria
deferribacteres | 2 | _;_bacteria _;_bacteria;_deferribacteres AKA _deferribacteraeota _flexistipes_group
dictyoglomi | 2 | _;_bacteria;_dictyoglomi _;_bacteria AKA _dictyoglomi _dictyoglomaeota _dictyoglomus_group
e
*/

class synonym_canon_lut 
{
typedef synonym_canon_lut Myt;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef  Tr::Dmel  Dmel; 
enum{ BAD=~0};
// may have levels 
// don't use tokenizer here, O(0)
typedef IdxTy Level;
typedef StrTy TrivialName;
typedef StrTy CanonicalName;

typedef std::map<TrivialName, CanonicalName> CanonLUT;
typedef std::map<Level, CanonLUT> CanonMap;

typedef std::vector<TrivialName> TrivialNames;
typedef std::map<CanonicalName, TrivialNames> SynonymLUT;
typedef std::map<Level, SynonymLUT> SynonymMap;

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public:
IdxTy size() const { return m_fwd.size(); } 
template <class Tv> 
	void conform_vector(Tv & taxa)
{
MM_SZ_LOOP(i,taxa,sz)
{
auto & v=taxa[i];
// TODO nonconst wtf 
if (m_fwd.find(i)!=m_fwd.end()) 
{
auto & mfi=m_fwd[i];
 if ( mfi.find(taxa[i])!=mfi.end())
v=m_fwd[i][v];
if (v!=m_fwd[i][v])
{
	//DMel("non_idempotent_map", m_ss<<MM_STR_LOC<<MMPR3(v,m_fwd[i][v],i));
	//MM_ONCE("non_idempotent_map adding "<< m_ss<<MM_STR_LOC<<MMPR3(v,m_fwd[i][v],i),)
	MM_ONCE("non_idempotent_map adding "<< m_ss.str()<<MM_STR_LOC<<MMPR3(v,m_fwd[i][v],i),)
	m_fwd[i][v]=v;
}	
}

} // i 



}
// synonyms need to be lumped to the right lineage 
// 5265  run_zymo -canon-check emp_canon_file.txt 2>qqq > canon_check.1
// 5267  run_zymo -canon-check emp_canon_file.txt 2>qqq | tee canon_check.1
// dictyoglomi | 2 | _;_bacteria;_dictyoglomi _;_bacteria AKA _dictyoglomi _dictyoglomaeota _dictyoglomus_group
template <class Ty> 
void add_tax_line(const Ty & v, const StrTy & line, const IdxTy & flags=0 )
{

	const IdxTy sz=v.size();
	if (sz<5)
	{
		DMel("too few words", m_ss<<MM_STR_LOC<<MMPR2(sz,line));
	return; // 	continue;
	}
	IdxTy pc=0;
	const StrTy target=v[pc]; ++pc;
	const IdxTy hits=myatoi(v[pc+1]);
	const IdxTy akaloc=pc+3+hits;
 	const bool fixed_loc=(v[pc]=="|")&&(v[pc+2]=="|")&&(v[akaloc]=="AKA");	
	if (!fixed_loc)
	{
 		const bool fixed_loc1=(v[pc]=="|");
		const bool fixed_loc2=(v[pc+2]=="|");
		const bool fixed_loc3=(v[akaloc]=="AKA");	
		DMel("fixed_locations_wrong", m_ss<<MM_STR_LOC<<MMPR4(fixed_loc1,fixed_loc2,fixed_loc3,line)<<MMPR(akaloc));
		return; // continue;

	}
	pc+=3;
	// right now all syns work on all levels 
	std::vector<IdxTy> levels;
	std::vector<StrTy> syns;
	while (pc<akaloc)
	{
		const StrTy hier=v[pc];
		const IdxTy len=hier.length();
		IdxTy level=0;
		for(IdxTy i=0; i<len; ++i) { if (hier.c_str()[i]==';') ++level; }
		levels.push_back(level); // note that this is off by one AFAICT 
		++pc;
	}
	++pc;
	while (pc<sz)
	{
		syns.push_back(v[pc].substr(1));
		++pc;
	}
	MM_SZ_LOOP(i,levels,szl)
	{
		IdxTy level=levels[i];
		m_fwd[level][target]=target;
		// TODO subtract one? 
		MM_LOOP(ii,syns)
		{
			// TODO need to check that target maps to target and that this
			// is not already there ... 
			m_fwd[level][(*ii)]=target;
		} // syns

	} // i 
}

void assemble_canon_map(const StrTy & str)
{
	IdxTy cnt=0;
	std::ifstream is (str);
	CommandInterpretter li(&is);
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		 add_tax_line(li.words(), li.line(), 0 );
		++cnt;
		//MM_ERR(" processing "<<li.dump())
	} // li
	MM_ERR(" done loadiing "<<MMPR2(str,cnt))
} // assemble_canon_map


//typedef std::map<TrivialName, CanonicalName> CanonLUT;
//typedef std::map<Level, CanonLUT> CanonMap;

void load_canon_map(const StrTy & str) { std::ifstream is (str); load_canon_map(is); }
void load_canon_map(IsTy & is)
{	
	const IdxTy base=0; 
	IdxTy cnt=0;
	IdxTy skip=0;
	CommandInterpretter li(&is);
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<(base+3)) { ++skip; continue; }
		const IdxTy level=myatoi(li.word(base));
		const StrTy & t=li.word(base+1);
		const StrTy & c=li.word(base+2);
		m_fwd[level][t]=c;
		++cnt;
	} // li 
	MM_ERR(" finish loading canon "<<MMPR2(cnt,skip))
} // load 
void dump_canon_map(OsTy & os)
{
MM_ERR(" dumping canon map ")
const StrTy sep=" ";
MM_LOOP(ii,m_fwd)
{
	const IdxTy level=(*ii).first;
	const CanonLUT & m=(*ii).second;
	MM_LOOP(jj,m)
	{
		const TrivialName & t=(*jj).first;
		const CanonicalName & c=(*jj).second;
		os<<level<<sep<<t<<sep<<c<<CRLF;
	} // jj 

} // ii 


} // dump_canon_map

private:
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
void DMel(const StrTy & code, OsTy & _ss)
{
	Ss& ss=dynamic_cast<Ss& >(_ss);
	std::cerr<<ss.str()<<" "<<code<<CRLF;
	ss.str(StrTy(""));
}






Dmel * m_dmel;
IdxTy m_serial;
Ss m_ss;
CanonMap m_fwd;
SynonymMap m_rev;

}; // synonym_canon_lut 


class otu_struct 
{
typedef otu_struct Myt;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;

//typedef  data_model_error_log Dmel; 
typedef  Tr::Dmel  Dmel; 
typedef StrTy HierTy;
//typedef std::vector<HierTy> AsgnVec;
typedef std::vector<IdxTy> AsgnVec;

// TODO FIXME danger will robinson memory leak here as garbage bin can not be cleared 
typedef string_tokenizer St;
static St&  st()
{
static St s;
return s;
}
//typedef std::map<StrTy, IdxTy> Ttok;
//static Ttok & tokens() 
//{
//static Ttok x;

//return x; 

//}
public:
otu_struct() {}
otu_struct(const StrTy & name ):m_otu(name) {}
const IdxTy size() const { return m_hier.size(); }
//const AsgnVec & tax() const { return m_hier; } 
//const StrTy & flat() const { return m_flat;} 
const StrTy string(const IdxTy flags) const
{
Ss ss;
const StrTy sep=" ";
IdxTy i=0;
const IdxTy sel_mask=256;
const bool to_lower=((flags&1)==0);
const bool remove_pfx=((flags&2)==0); // true;
const bool make_single=((flags&4)==0); // true;
MM_LOOP(ii,m_hier)
{
	const bool exclude=((flags&(sel_mask<<i))!=0);
	//MM_ERR(MMPR4(exclude,flags,i,sel_mask))
	if (exclude) { ++i; continue; }
	const HierTy & hi=st()(*ii);
	StrTy s=hi;
//	MM_ERR(" adding "<<MMPR(s)) 
	const IdxTy len=hi.length();	
	s=cuname(s,to_lower,remove_pfx,make_single);
	//if ( to_lower) s=mjm_strings::to_lower(s);
	//if (remove_pfx)  s=RemovePfx(s);
	//if ( make_single) mjm_strings::make_single(s);
	ss<<s<<sep;
	++i;
} // ii 

return ss.str();
} //string
// this is obsoleete but still callled somewhere 
template <class Tm> IdxTy classify(const Tm & m) const
{
// get highest( lowest?)  level hit that matches anything 
MM_ERR(" dont use this ")
MM_LOOP(ii,m_hier)
{
	const bool to_lower=true;
	const bool remove_pfx=true;
	const bool make_single=true;
const HierTy & hi=st()(*ii);
	StrTy s=hi;
	const IdxTy len=hi.length();	
	s=cuname(s,to_lower,remove_pfx,make_single);
	//if ( to_lower) s=mjm_strings::to_lower(s);
	//if (remove_pfx)  s=RemovePfx(s);
	//if ( make_single) mjm_strings::make_single(s);
	const auto mf=m.find(s);
	if (mf!=m.end()) { return (*mf).second; } 

} // ii 
return 0; 
}

template <class Tm> StrTy  classify_name(const Tm & m) const
{
// get highest(lwest?)  level hit that matches anything 
//MM_LOOP(ii,m_hier)
const IdxTy sz=m_hier.size();
for(IdxTy i=0; i<sz; ++i)
{
	const IdxTy ii=sz-i-1;
	const bool to_lower=true;
	const bool remove_pfx=true;
	const bool make_single=true;
const HierTy & hi=st()(m_hier[ii]); // (*ii);
	StrTy s=hi;
	const IdxTy len=hi.length();	
	s=cuname(s,to_lower,remove_pfx,make_single);
//	if ( to_lower) s=mjm_strings::to_lower(s);
//	if (remove_pfx)  s=RemovePfx(s);
//	if ( make_single) mjm_strings::make_single(s);
	//MM_ERR(" finding "<<MMPR2(s,m.size()))
	const auto mf=m.find(s);
	// a set would work here but still 
	if (mf!=m.end()) { //MM_ERR(" returning a hit "<<MMPR((*mf).first))  

return (*mf).first; } 

} // ii 
const StrTy flat=string(0);
if (flat.length()==0) { MM_ONCE(" blanks in otu detected ",) } 
else { MM_ERR(" no match for "<< string(0))  } 
return StrTy("nomatch"); 
}

static StrTy cuname(const StrTy & s,const bool to_lower,const bool remove_pfx,const bool make_single, const bool conform=true) 

{
	StrTy  r=s;
	if ( to_lower) r=mjm_strings::to_lower(r);
	if (remove_pfx)  r=RemovePfx(r);
	if ( make_single) mjm_strings::make_single(r);
	if (conform) conform_map(r);
return r; 
}
StrTy  uname(const StrTy & s) const  { StrTy r= RemovePfx(mjm_strings::to_lower(s)); mjm_strings::make_single(r); return conform_map(r); } 
typedef std::map<StrTy, StrTy>  ConformMap ;  
static ConformMap &  conform_map () { static ConformMap m; return m; } 
static StrTy  conform_map (const StrTy & s) { auto ii=conform_map().find(s); if (ii==conform_map().end()) return s;
else return (*ii).second;  } 
StrTy level(const IdxTy n, const bool conform)  const 
{
const IdxTy sz=m_hier.size();
if (n>=sz) return StrTy("");
return  cuname(st()(m_hier[n]),true,true,true,conform);
}
//const IdxTy size() const { return m_hier.size(); } 
// a nested map  of genus and species with classification
template <class Tm> StrTy  gs_classify(const Tm & m) const
{
IdxTy locs=6;
IdxTy locg=5;
const StrTy nomatch="nomatch";
const IdxTy sz=m_hier.size();
MM_ERR( " classify gs"<<MMPR2(sz,locs)) 
if (sz>locs)
{
 	const StrTy species= uname(st()(m_hier[locs]));
	const StrTy genus= uname(st()(m_hier[locg]));
	MM_ERR(MMPR2(species,genus))
	const auto mf=m.find(genus);
	if (mf!=m.end())
	{
		const auto & ms=(*mf).second;
		const auto & msi=ms.find(species);
		// each genus must have a non-blank default 
		if ( msi==ms.end()) { return (*(ms.find(StrTy("")))).second; }
		return (*msi).second;
	}	

}

return nomatch;
}

Myt & push(const HierTy & t ) { // m_flat=m_flat+StrTy(" ")+t; 
	m_hier.push_back(st()(t)); return *this;  } 

Myt & pop( IdxTy n=1 ) { // m_flat=m_flat+StrTy(" ")+t; 
 while (n>0 ) { m_hier.pop_back(); --n; } return *this;  } 

StrTy dump(const IdxTy flags=0) const
{
const StrTy sep=" ";
Ss ss;
ss<<m_otu;
const IdxTy sz=size();
for (IdxTy i=0; i<sz; ++i) { ss<<sep<<(st()(m_hier[i])); } 
return ss.str();
}

std::vector<StrTy > taxon(const IdxTy flags=0) const
{
	std::vector<StrTy> v;
if (flags==0) { 	MM_SZ_LOOP(i,m_hier,sz) { v.push_back(st()(m_hier[i])); }}
else
{
	MM_SZ_LOOP(i,m_hier,sz)
	{
	StrTy s= st()(m_hier[i]);
	const bool to_lower=true,remove_pfx=true,make_single=true;
	s=cuname(s,to_lower,remove_pfx,make_single);
	v.push_back(s);
}
}
	return v;
}

/*
virtual bool conform(const IdxTy flags=0) 
{
bool ok=true;
for (auto ii=m_hier.begin(); ii!=m_hier.end(); ++ii)
{
	HierTy & hi=st()(*ii);
	const IdxTy len=hi.length();	

}

return ok;
} */

template <class Tmap> void catalog( Tmap & m) const
{
MM_LOOP(ii,m_hier)
{
m[st()(*ii)]+=1;
}
}

private:
// this does leave in things like s__ which is ok for marking now. 
static StrTy RemovePfx( const StrTy & s) // const
{
const IdxTy sz=s.length();
if (sz<3) return s; 
StrTy d=s;
const char * cs=s.c_str();
for (IdxTy i=0; i<(sz-2); ++i )
{
if (cs[i]=='_') if (cs[i+1]=='_') return s.substr(i+2); 

} // i 
return d; 
}
const StrTy & name() const { return m_otu; } 
private:
StrTy m_otu;
AsgnVec m_hier;
// StrTy m_flat;
// should have the "reads" that map to this along with parsed
// quality issues


 // DB featuers - pointers to other structs, data about composition
// abundance not related to otu per se. 


}; // otu_struct


class otu_collection : public std::map<StrTy, otu_struct>
{
typedef otu_collection Myt;
typedef  std::map<StrTy, otu_struct> Super;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef  Tr::Dmel  Dmel; 

//typedef string_tokenizer St;

public:
otu_collection() {}
otu_collection(const StrTy & name ):m_name(name) {}
const StrTy & name() const { return m_name; } 
// this is a thing with first element as otu name followed by taxa defs
template <class Ty> Myt & add( const Ty & v) 
{
if (v.size()<1) return *this;
auto ii=v.begin();
const StrTy name=(*ii);
otu_struct os= otu_struct(name);
++ii;
for (; ii!=v.end(); ++ii)
//MM_LOOP(ii,v)
{
//MM_ERR(" pushing "<<(*ii))
os.push(*ii);
//const auto & line=(*ii);
//for (auto jj=line.begin(); jj!=line.end(); ++jj)
//MM_LOOP(jj,line)
//{

//} // jj 

}
if (find(name)!=end())
{
MM_ERR(" adding another ostu named "<<name)
}
(*this)[name]=os;
return *this;
}

template <class Ty> Myt & add( const StrTy & id, const Ty & v) 
{
otu_struct os= otu_struct(id);
MM_LOOP(ii,v) { os.push(*ii); }
if (find(id)!=end())
{ MM_ERR(" adding another ostu named "<<id) }
(*this)[id]=os;
return *this;
}

bool has(const StrTy & id) const
{ return (find(id)!=end()); }

//const IdxTy size() const { return size(); }
StrTy dump(const IdxTy flags,const StrTy & nm) const
{
Ss ss;
const StrTy sep=" ";
for (auto ii=begin(); ii!=end(); ++ii)
{
ss<<nm<<sep;
ss<<(*ii).second.dump()<<CRLF;

}

return ss.str();
}

template <class Tmap> void catalog( Tmap & m) const
{
MM_LOOP(ii,(*this))
{
(*ii).second.catalog(m);
}
}
private:



private:
StrTy m_name;


}; // otu_collection

/*
 A sample is a thing that was analyzedand the analysis attributes are included
in that. A class number comes from groups different from processing groups.

*/

class seq_sample 
{
typedef seq_sample Myt;
typedef zymo_traits::Tr Tr;
protected:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::MyBlock  MyBlock;
typedef Tr::Dmel  Dmel; 
// TODO FIXME memory leak when the collection gets deleted need to make this non-static 
typedef string_tokenizer St;
/*
static St & st() 
{
static St st;
return st;
}
*/

public:
class otu_info
{
typedef otu_info Myt;
public:
otu_info() : m_otu(),m_hits(0) {}
//otu_info(const StrTy & nm) : m_otu(st()(nm)),m_hits(0) {}
otu_info(const IdxTy  & nt) : m_otu(nt),m_hits(0) {}
//const StrTy & otu() const { return st()(m_otu);}
const StrTy & otu(St & st ) const { return st(m_otu);}
const IdxTy & otu( ) const { return (m_otu);}
//void  otu(const StrTy & otu )  {  m_otu=st()(otu);}
void  otu(St & st, const StrTy & otu )  {  m_otu=st(otu);}
void  otu( const IdxTy & otu )  {  m_otu=(otu);}
const D & hits() const { return m_hits; } 
bool operator ==( const D & x ) { return (m_hits==x); } 
Myt & operator +=( const D & x ) {m_hits+=x; return *this; } 
//otu_struct m_otu;
//StrTy  m_otu;
IdxTy   m_otu;
D m_hits;
}; // otu_info 
//typedef std::map<StrTy, otu_info> DataMap;
typedef std::map<IdxTy, otu_info> DataMap;
typedef DataMap::const_iterator CItor;
public:
//seq_sample(): m_class(~0)  {}
seq_sample(): m_class(~0)  {}
CItor  begin() const  { return  m_map.begin(); }
CItor  end() const { return  m_map.end(); }
//const IdxTy & operator()( const StrTy & x)   { return m_st(x); } 
//const StrTy & operator()( const IdxTy & x) const  { return m_st(x); } 

//string_tokenizer & operator() { return m_st; } 
// operator string_tokenizer() { return m_st; } 
#if 0
void set_name(const StrTy & nm)
{
const IdxTy & t=m_st(nm);
//m_map[nm].otu(nm);
m_map[t].otu(t);
}
#endif
void set_class( const IdxTy & c) { m_class=c;}
const IdxTy & get_class() const { return m_class; }  

// this sets an OTU name, 
//void set_name(const IdxTy & t) { m_map[t].otu(t); }
void set_sample_name(const StrTy  & t) { m_name=t; }
void set_otu_token(const IdxTy & t) { m_map[t].otu(t); }

//Myt & operator +=( const D & x ) {m_hits+=x; return *this; } 
#if 0
Myt & add( const StrTy & nm , const D & v)
{
const auto & t=m_st(nm);
auto & mm=m_map[t];
//auto & mm=m_map[nm];
//if (mm==0) mm.otu(nm);
if (mm==0) mm.otu(t);
//m_map[nm]+=v;
mm+=v;
return *this; 
}
#endif

Myt & add( const IdxTy & t , const D & v)
{
//const auto & t=m_st(nm);
auto & mm=m_map[t];
if (mm==0) mm.otu(t);
mm+=v;
return *this; 
}


StrTy dump(const St & st,  const IdxTy flags=0, const StrTy&  label=StrTy(" ")) const
{
Ss ss;
const StrTy sep=" ";
MM_LOOP(ii,m_map)
{
const D & hit=(*ii).second.hits();
// right now nameis ntused 
//if (hit!=0) { 	ss<<label<<sep<<m_name<<sep<<(*ii).first<<sep<<hit<<CRLF; } 
//if (hit!=0) { 	ss<<label<<sep<<(m_st((*ii).first))<<sep<<hit<<CRLF; } 
if (hit!=0) { 	ss<<label<<sep<<(st((*ii).first))<<sep<<hit<<CRLF; } 
} // ii
return ss.str();
}

private:
	StrTy m_name;
	DataMap m_map;
	IdxTy m_class;
//	St m_st;
}; // seq_sample


class mjm_zymo 
{
typedef  zymo_traits::Tr  Tr;
typedef mjm_zymo Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

//typedef Tr::MySparse MySparse;

typedef zymo_params Logic;
typedef mjm_logic_base VariableStore;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

typedef string_tokenizer Tokenizer;

typedef mjm_column_plot<Tr> ColPlot;

//typedef mjm_calendar CalTy;

////////////////////////////////////////////////////////
typedef otu_collection OtuCol;
typedef std::map<StrTy,OtuCol> OtuMap; 
typedef seq_sample Sample;

typedef mjm_tax_tree TaxTree;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef Ragged::Line Line;
typedef std::vector<Line> Lines;


typedef mjm_pheno_notes PhenoNotes;

typedef mjm_kv_vector<StrTy,D,Tr> SumVec;


typedef mjm_biom_to_rag<Tr> BiomToRag;
typedef mjm_rag_reduce<Tr> RagReduce;

typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy , Fasta> FastaMap;
typedef mjm_sequence_reconcile<Tr> Recon;
//typedef std::map<StrTy,Sample> SampleMap; 
//class sample_map : public std::map<StrTy, Sample>
class sample_map : public std::map<IdxTy, Sample>
{
typedef sample_map Myt;
typedef std::map<IdxTy, Sample>  Super;
typedef Super::iterator Itor;
public:
typedef string_tokenizer St;
St & st() { return m_st;}
// FIXME TODO need to move this in with iteror and overload etc 
const St & st() const { return m_st;}
const IdxTy & operator()( const StrTy & x)   { if (x.c_str()[0]==0) return m_st(StrTy("blank")); return m_st(x); } 
const StrTy & operator()( const IdxTy & x) const  { return m_st(x); } 
// infinite loop wtf 
//Sample & operator[](const StrTy & s) { return Super::operator[(m_st(s))]; } 
//Sample & operator[](const StrTy & s) { return Super(*this)[(m_st(s))]; } 
//Itor find(const StrTy & s) { return Super::find(m_st(s)); } 
//CItor find(const StrTy & s)const  { return Super::find(m_st(s)); } 

string_tokenizer m_st;

}; // sample_map




//typedef std::map<StrTy,Sample> SampleMap; 
typedef sample_map  SampleMap; 
typedef std::map<StrTy,SampleMap> SampleGroup; 
//typedef mjm_biom_hdf5::biom_json BiomFile;
typedef biom_json BiomFile;
typedef std::map<StrTy, BiomFile> BiomFiles;

typedef std::map<StrTy, IdxTy > CatagoryMap;
typedef std::map<StrTy, CatagoryMap> CatagoryMaps;
//typedef std::map<StrTy, IdxTy > OrderMap;
typedef collate_order_map OrderMap;

typedef std::map<StrTy, OrderMap> CollateOrders;


typedef std::map<StrTy, IdxTy > SampleFilter;
typedef std::map<StrTy, SampleFilter> SampleFilters;

typedef std::map<StrTy, StrTy> SpeciesMap;
typedef std::map<StrTy, SpeciesMap> GenusMap;

typedef mjm_zymo_composition_table<Tr> CompTable;
typedef std::map<StrTy,CompTable> CompMap;




class genus_species_map
{
public:
typedef GenusMap Gm;

Gm m_map;

}; // genus_species_map
typedef genus_species_map GSMap;

typedef std::map<StrTy, GSMap > GSMaps;


////////////////////////////////////////////////////


public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const {
//::strtol(s.c_str(),0,0);
return ::strtol(c,0,0);
}

public :
mjm_zymo():m_dmel(new Dmel()) {Init();}
mjm_zymo(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_zymo()
{
//clear_handlers();
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
dest=myatoi(args[i]);
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
if (s=="-about") { ++i; about(); }
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
StrTy xwords(const IdxTy first,const std::vector<StrTy> & w) 
{
StrTy s;
IdxTy i=first;
if (i<w.size()) s=w[i];
++i;
for (; i<w.size(); ++i) s=s+StrTy(" ")+w[i];
return s;
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
if (cmd=="") continue;
if (cmd.c_str()[0]=='#' ) continue; 
const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
const StrTy p3=(sz>3)?li.word(3):StrTy("");
const StrTy p4=(sz>4)?li.word(4):StrTy("");
const StrTy p5=(sz>5)?li.word(5):StrTy("");
const StrTy p6=(sz>6)?li.word(6):StrTy("");
if (cmd=="about") { about();  continue; } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="print") { MM_MSG(li.line())  continue; } 
if (cmd=="err") { MM_ERR(li.line())  continue; } 
if (cmd=="source") 
	{ { MM_ERR(" sourcing from "<<p1) } command_modef(p1.c_str()); continue;  }
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
//if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
//if (cmd=="parse") { parse();  continue; } 
//if (cmd=="pair") { pair_stats(li.word(1),li.word(2),li.word(3));  continue; } 
if (cmd=="pn") { MM_ERR(" enter pheno_notes ") m_pheno_notes.command_mode(li); MM_ERR(" return to zymo ")   continue; } 
if (cmd=="pheno-notes") {  m_pheno_notes.command_mode(p1);  continue; } 
if (cmd=="parse-otu") { parse_otu(li.word(1),li.word(2));  continue; } 
if (cmd=="dump-otu") { if (p1.length()>0)  {Ofs ofs(p1); dump_otu(ofs,p1); } else dump_otu(std::cout,"std::cout");  continue; } 
if (cmd=="dump-sample-group") { if (sz>2) { Ofs ofs(p2); dump_sample_group(p1,ofs); }  else dump_sample_group(p1);   continue; } 
if (cmd=="assemble-canon-lut") { m_canon_lut.assemble_canon_map(p1);  continue; } 
if (cmd=="load-canon-lut") { m_canon_lut.load_canon_map(p1);  continue; } 
if (cmd=="dump-canon-lut") { m_canon_lut.dump_canon_map(std::cout);  continue; } 
if (cmd=="read-canon-map") { read_canon_map(p1);  continue; } 
if (cmd=="read-ragged") { m_ragged_map[p1].load(p2);  continue; } 
if (cmd=="clear-ragged") { m_ragged_map[p1].clear();  continue; } 
if (cmd=="add-ragged") {  m_ragged_map[p1].add(xwords(2,li.words()));  continue; } 
if (cmd=="sum-ragged") { sum_vec(m_ragged_map[p1],m_ragged_map[p2],myatoi(p3),myatoi(p4),myatoi(p5));  continue; } 
if (cmd=="dump-ragged") { MM_MSG(m_ragged_map[p1].dump()) ;  continue; } 
if (cmd=="ssv-ragged") { std::ofstream ofs(p1); ofs<<m_ragged_map[p2].dump_ssv();  continue; } 
if (cmd=="read-taxa-conform") { read_taxa_conform(p1,p2);  continue; } 
if (cmd=="dump-taxa-conform") { dump_taxa_conform(p1,p2);  continue; } 
if (cmd=="catalog-taxa-names") { catalog_taxa_names();  continue; } 
if (cmd=="col-plot") { 
// needs to set ppt.m_nkeys if it is used... 
ColPlot cp;
ColPlot::plot_params_type ppt;
ppt.m_colx=myatoi(p3);
ppt.m_coly=myatoi(p4);
ppt.m_colc=myatoi(p5);
const IdxTy flags=myatoi(p6);
ppt.m_xx=flags&3;
ppt.m_xy=(flags>>2)&3;

IdxTy rc=cp.r_plot(m_ragged_map[p1],m_ragged_map[p2],ppt,flags>>4);
MM_ERR(MMPR(ppt.dump())<<MMPR4(m_ragged_map[p1].size(), m_ragged_map[p2].size(),p1,p2))

 continue; } 

if (cmd=="col-bins") { 
ColPlot cp;
IdxTy rc=cp.r_stat(m_ragged_map[p1],m_ragged_map[p2],myatoi(p3),myatoi(p4));
MM_ERR(MMPR4(m_ragged_map[p1].size(), m_ragged_map[p2].size(),myatoi(p3),myatoi(p4)))
continue;
} // col-bins

if (cmd=="col-stuff") { 
ColPlot cp;
//  vals,flags
IdxTy rc=cp.r_stuff(m_ragged_map[p1],m_ragged_map[p2],myatoi(p3),myatoi(p4));
MM_ERR(MMPR4(m_ragged_map[p1].size(), m_ragged_map[p2].size(),myatoi(p3),myatoi(p4)))
continue;
} // col-stuff
/*
 Convert 1 or more biom files into ssv seq,sample,count,otu0,otu1.... 
2544  ./mjm_biom_to_rag.out  "b2r b2.biom 0 b2.ssv" quit
 2545  ./mjm_biom_to_rag.out  "b2r b3.biom 0 b3.ssv" quit
 2546  . srecord
Sum counts by sample and n key fields
 2559  ./mjm_rag_reduce.out "reduce s3.ssv b3.ssv S0S2K4-9K1" quit 2>&1 
 2560  ./mjm_rag_reduce.out "reduce s2.ssv b2.ssv S0S2K4-9K1" quit 2>&1 
 2561  ./mjm_rag_reduce.out "reduce s1.ssv b1.ssv S0S2K4-9K1" quit 2>&1 

Put all sample counts on one line followed by otu0,otu1... 

 2562  ./mjm_rag_reduce.out "merge m13.ssv s1.ssv s3.ssv" quit 2>&1 
 2563  more m13.ssv | grep Porp


*/

if (cmd=="bioms-to-ssv") { 
CommandInterpretterParam  cip(li);
const IdxTy rc= cmd_bioms_to_ssv(li);
MM_ERR(MMPR(rc))
 continue; } // 
if (cmd=="reconciled") { 
MM_ERR(MMPR(cmd))
CommandInterpretterParam  cip(li);
const IdxTy rc= cmd_reconciled(li);
MM_ERR(MMPR(rc))
 continue; } // 
// do common things to a count-key format ragged
if (cmd=="count-key-select") { 
MM_ERR(" fack ")
CommandInterpretterParam  cip(li);
const IdxTy rc= cmd_count_key_select(li);
MM_ERR(MMPR(rc))
 continue; } // 
if (cmd=="load-recon") { 
MM_ERR(" fack ")
CommandInterpretterParam  cip(li);
const IdxTy rc= cmd_load_recon(li);
MM_ERR(MMPR(rc))
 continue; } // 
if (cmd=="find-olaps") { 
MM_ERR(" fack ")
CommandInterpretterParam  cip(li);
const IdxTy rc= cmd_find_olaps(li);
MM_ERR(MMPR(rc))
 continue; } // 


if (cmd=="quiet") { _quiet(); continue; } 
if (cmd=="nquiet") { _nquiet(); continue; } 


if (cmd=="list") { list(std::cout);  continue; } 
if (cmd=="catalog-otu") { catalog_otu(std::cout);  continue; } 
if (cmd=="read-catagory-map") { load_catagory_map(li.word(1),li.word(2));  continue; } 
if (cmd=="dump-catagory-map") { dump_catagory_map(li.word(1));  continue; } 
if (cmd=="add-catagory") { add_catagory(li.word(1),li.word(2));  continue; } 
if (cmd=="clear-orders") { m_orders.clear();  continue; } 
if (cmd=="dump-orders") { dump_orders();  continue; } 
//if (cmd=="dump") { MM_ERR(x.dump()) }
//else if (cmd=="load") { x.load(cip.p1,atoi(cip.p2.c_str())); }
//else if (cmd=="load_notes") { x.load_notes(cip.p1,atoi(cip.p2.c_str())); }
//else if (cmd=="save") { x.save(cip.p1,cip.p2,0); }
//else if (cmd=="save_gs") { x.save_gs(cip.p1,atoi(cip.p2.c_str())); }

if (cmd=="dump-comp") { MM_MSG(m_comps_map[p1].dump(myatoi(p2)));  continue; } 
if (cmd=="save-comp") { m_comps_map[p1].save(p2,(p3),myatoi(p4));  continue; } 
if (cmd=="merge-comp") { m_comps_map[p1].merge(m_comps_map[p2],myatoi(p3));  continue; } 
if (cmd=="sync-comp") { m_comps_map[p1].sync_notes(myatoi(p2));  continue; } 
if (cmd=="save-comp-ssv") { m_comps_map[p1].save_ssv(p2);  continue; } 
if (cmd=="save-comp-notes-ssv") { m_comps_map[p1].save_notes_ssv(p2);  continue; } 
if (cmd=="save-comp-gs") { m_comps_map[p1].save_gs(p2,myatoi(p3));  continue; } 
if (cmd=="read-comp") { m_comps_map[p1].load(p2,myatoi(p3));  continue; } 
if (cmd=="read-comp-notes") { m_comps_map[p1].load_notes(p2,myatoi(p3));  continue; } 

// this should call li.push and then pass the thing to it 
if (cmd=="tree-viz") { m_tvz.command_mode(p1);   continue; } 
if (cmd=="tv") { MM_ERR(" enter tvz") m_tvz.command_mode(li);  MM_ERR(" exit tvz" )   continue; } 

if (cmd=="otf-angles") { otf_angles(p1,p2);  continue; } 
if (cmd=="dump-otf-class-stats") { dump_otf_class_stats(p1);  continue; } 
//void dump_otf_class_stats(const StrTy & fn)

if (cmd=="read-sample-filter") { read_sample_filter(li.word(1),li.word(2));  continue; } 
//void  add_gs_map(const StrTy & mapname, const StrTy & fn) 
if (cmd=="read-gs-map") { add_gs_map(p1,p2);  continue; } 
//biom_json x;
//x.parse(tax,atoi(sample.c_str()));
if (cmd=="parse-biom-json") { parse_biom_json(p1,p2);  continue; } 
if (cmd=="parse-biom-json-no-data") { parse_biom_json(p1,p2,true);  continue; } 
if (cmd=="parse-biom-json-to-file") { parse_biom_json(p1,p2,!true,li.word(3));  continue; } 
if (cmd=="dot-sparse-biom") { dot_sparse_biom(li.word(1),li.word(2));  continue; } 
if (cmd=="copy-biom") { copy_biom(p1,p2,0);  continue; } 
if (cmd=="copy-biom-canon") { copy_biom(p1,p2,7);  continue; } 
if (cmd=="copy-biom-conform") { copy_biom(p1,p2,6);  continue; } 
if (cmd=="move-biom") { move_biom(p1,p2);  continue; } 
if (cmd=="clear-biom") { clear_biom(p1);  continue; } 
if (cmd=="clear-all-biom") { clear_biom(StrTy(""));  continue; } 

// this works for zymo format but notin the frontiers paper
if (cmd=="parse-reads") { parse_reads(p1,p2);  continue; } 
// this is 
if (cmd=="parse-front-reads") 
	{ parse_front_reads(li.word(1),li.word(2));  continue; } 
if (cmd=="dump-reads") { dump_reads(std::cout);  continue; } 
if (cmd=="sum-samples") { sum_samples(li.word(1),myatoi(li.word(2).c_str()));  continue; } 

	//const bool do_raw= ((settings&1)!=0);
	//const bool do_collate= ((settings&2)!=0);
	//const bool update_order= ((settings&4)!=0);
	//const bool use_map= ((settings&8)!=0);

if (cmd=="add-order") { m_orders[""].add_at(p1,myatoi(p2));  continue; } 
if (cmd=="dump-select-n") { dump_select(p1,myatoi(p2),p3);  continue; } 
//void adhoc_tree(const StrTy &list,const IdxTy & settings,  std::ostream  & os=std::cout)
if (cmd=="adhoc-tree") { adhoc_tree(p1,myatoi(p2),p3);  continue; } 
if (cmd=="adhoc-rags") { adhoc_rags(p1,myatoi(p2),p3);  continue; } 


//if (cmd=="dump-group-raw") { dump_group(li.word(1),true,false);  continue; } 
if (cmd=="dump-group-n") { dump_group(p1,myatoi(p2));  continue; } 
if (cmd=="dump-group-raw") { dump_group(p1,13);  continue; } 
// this appears now to make a decent bar plto with coll.R but not
// checked for accuracy beyond spot check with Zymo plot. 
if (cmd=="dump-group-collate") { dump_group(p1,14);  continue; } 
if (cmd=="dump-group-format") { dump_group(p1,6);  continue; } 
if (cmd=="dump-class-format") { dump_group(li.word(1),(6|16));  continue; } 
if (cmd=="update-order-format") {
		MM_ERR(" not right needs format doh include 8,16, or 32 ") ;  
		IdxTy n=(sz>2)?myatoi(p2):4;dump_group(p1,n);  
			continue; } 
if (cmd=="update-order-map") { dump_group(li.word(1),12);  continue; } 
if (cmd=="sort-order") { sort_order(p1);  continue; } 

if (cmd=="classify-sample-group") { classify_group(p1,p2);  continue; } 

if (cmd=="dump-otu-frac") {  dump_collated_group(std::cout, StrTy("dummy")); continue; } 
//m_tax_tree.standard_commnds(cmd,p1,p2,li);

if (cmd=="load-tax-nodes") { m_tax_tree.read_ncbi_nodes_dmp(p1);  continue; }
if (cmd=="load-tax-tree") { m_tax_tree.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { m_tax_tree.dump(p1);  continue; }
if (cmd=="write-tax-single") { m_tax_tree.write_single(p1);  continue; }
if (cmd=="save-tax") { m_tax_tree.write_composite(p1);  continue; }
if (cmd=="load-tax") {m_tax_tree.clear();  m_tax_tree.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { m_tax_tree.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { m_tax_tree.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { m_tax_tree.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { m_tax_tree.traverse_full_tree(p1,4);  continue; }
if (cmd=="best-taxon") {
std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  m_tax_tree.best_taxon(ww,4);  continue;
}
if (cmd=="best-indexed-taxon") {
std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  m_tax_tree.best_indexed_taxon(ww,4);  continue;
}
if (cmd=="node-info") { local_label=m_tax_tree.node_info(myatoi(p1),myatoi(p2)); MM_MSG((local_label))   continue; }
if (cmd=="lineage") { StrTy ll=m_tax_tree.lineage(myatoi(p1)); MM_MSG(MMPR2(myatoi(p1),ll))   continue; }
if (cmd=="node-no") { auto x=m_tax_tree.node_numbers(p1); Ss ss ; MM_LOOP(ii,x) { ss<<" "<<(*ii); } MM_MSG(ss.str())   continue; }

if (cmd=="reconcile-n") { reconcile(p1,myatoi(p2));  continue; } 
if (cmd=="group-tax-tree") { group_tax_tree(p1,myatoi(p2));  continue; } 






//if (cmd=="pair-all") { pair_stats(li.word(1));  continue; } 
//if (cmd=="dose-vectors") { if (sz>1) { dose_vectors(li.word(1));} else { MM_ERR( "need dogname ") }  continue; } 
//void pair_stats(const StrTy & x, const StrTy & y)
//if (cmd=="dump-drugs") { collect_all_drugs(); continue; } 
//if (cmd=="dump-drugs-cerr") { collect_all_drugs(2); continue; } 
//if (cmd=="list-drugs-cerr") { collect_all_drugs(3); continue; } 
//if (cmd=="list-dogs-cerr") { collect_all_drugs(4); continue; } 
//if (cmd=="dump-counts") { stream_counts(std::cout); continue; } 
//if (cmd=="dump-canon-fails") { dump_canon_fails(std::cout); continue; } 
//if (cmd=="clear-canon-fails") { clear_canon_fails(); continue; } 
if (cmd=="dump-dmel") { dump_dmel(std::cout); continue; } 
if (cmd=="dump-dmel-cerr") { dump_dmel(std::cerr); continue; } 
//if (cmd=="dump-dog-days") { dump_dog_days(std::cout,p1); continue; } 
//if (cmd=="dump-dog-times") { dump_dog_times(std::cout,p1); continue; } 
//if (cmd=="n-canon-fails") { MM_MSG(MMPR(m_canon_fails.size())) continue; } 

//if (cmd=="clear-order") { clear_order(); continue; } 
//if (cmd=="dump-order") { dump_order(std::cout); continue; } 
//if (cmd=="push-order") 
		//{ int i=1; while ( li.cmd_ok(i+1)){ push_order(li.word(i)); ++i; } continue; } 
//		{ int i=1; while ( i<sz){ push_order(li.word(i)); ++i; } continue; } 
//if (cmd=="remove-order") 
//		{ int i=1; while ( i<sz){ remove_order(li.word(i)); ++i; } continue; } 

//if (cmd=="load-order") { if( li.cmd_ok(2)) load_order(li.word(1)); continue; } 
//if (cmd=="order-all") { collect_all_drugs(1); continue; } 

//if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
//if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
//if (cmd=="process") { if (li.cmd_ok(2)) process_file(li.word(1));  continue; } 
//if (cmd=="add") { if (li.cmd_ok(2)) add_timeline(li.word(1));  continue; } 
if (cmd=="init") { Init();  continue; } 
//if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");
if (m_flp.exit_on_err())
{
MM_ERR(" quiting "<<MMPR(m_flp.exit_on_err()))
clean_up();
return; 

}
else
{
MM_ERR("about print err set-param get-param parse-otu dump-otu dump-sample-group assemble-canon-lut load-canon-lut dump-canon-lut read-canon-map read-taxa-conform dump-taxa-conform catalog-taxa-names list catalog-otu read-catagory-map dump-catagory-map add-catagory clear-orders dump-orders otf-angles dump-otf-class-stats read-sample-filter read-gs-map parse-biom-json parse-biom-json-no-data parse-biom-json-to-file dot-sparse-biom copy-biom copy-biom-canon copy-biom-conform move-biom clear-biom clear-all-biom parse-reads parse-front-reads dump-reads sum-samples dump-group-n dump-group-raw dump-group-collate dump-group-format dump-class-format update-order-format update-order-map sort-order classify-sample-group dump-otu-frac load-tax-nodes load-tax-tree dump-tax-tree write-tax-single save-tax load-tax dump-lineages dump-normal-lineages check-normal-lineages traverse-full-tree best-taxon best-indexed-taxon node-info reconcile-n group-tax-tree dump-dmel dump-dmel-cerr init banner cm test quit")


} // help

}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
void about()
{
Ss ss;
ss<<" mjm_zymo "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
ss<<" Code to read various files related to 16S rRNA analyses "<<CRLF;
ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;
ss<<"database link provided by Bradley.stevenson@ou.edu"<<CRLF;
ss<<"http://www.earthmicrobiome.org/data-and-code/ "<<CRLF;
ss<<"Sample data provided by echen@zymoresearch.com "<<CRLF;
ss<<"http://www.isppweb.org/smc_files/bull%20et%20al.%202010%20jpp%20list.pdf"<<CRLF;
ss<<"http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.690.647&rep=rep1&type=pdf"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/Pathogenic_bacteria"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/List_of_bacteria_genera"<<CRLF;

#if 0
2522  more xxx | grep "^<td><i>" | sed -e 's/wiki/\nwiki/g' | grep wiki | awk '{print $1}'  | sed -e 's/.....//' | tr "[A-Z]" "[a-z]"| sed -e 's/_/ /g' | sed -e 's/"//g' | sed -e 's/m. /mycobacterium /g' | sed -e 's/n. /neisseria /g' | sed -e 's/b. /borrelia /g' | sed -e 's/(.*)//g' | sort | uniq | awk '{print $0" | wikipathogen" }' > wiki_pathogens.txt 
 2523  history | grep txt
 2524  pdftohtml -stdout xxx.pdf | sed -e 's/<i><b>/\n<i><b>/g' | grep "<b>"| sed -e 's/<i><b>//' | sed -e 's/<\/[ib]>.*//' | sed -e 's/&[^;]*;//g' | grep -v "[<>=]" | more | sort | uniq > plant_pathogens.txt
 2528  pdftohtml -stdout xxx.pdf | sed -e 's/<i><b>/\n<i><b>/g' | grep "<b>"| sed -e 's/<i><b>//' | sed -e 's/<\/[ib]>.*//' | sed -e 's/&[^;]*;//g' | grep -v "[<>=]" | tr "[A-Z]" "[a-z]" | sort | uniq | awk '{print $0" | plantrelated" }'  > plant_pathogens.txt 
#endif


std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();

}

// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fack 
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
MM_MSG(" configuration "<<m_flp.to_string())
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}
bool done() const  { return m_done; } 

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general
/*

StrTy tod_from_string(const StrTy & w)
{
if (w.length()!=6) return StrTy("BAD")+w;
return w.substr(0,4);
}

*/
// levels starts at 7 the number to sum 
void sum_vec(Ragged & d, const Ragged & s, const IdxTy nlevels, const IdxTy nsum, const IdxTy flags)
{
SumVec sv;
const IdxTy nx=nlevels-nsum;
MM_LOOP(ii,s)
{
const Line & l=(*ii);
Line k;
for(IdxTy i=0; i<nx; ++i) k.push_back(l[i]); 
std::vector<D> v;
for(IdxTy i=nlevels; i<l.size(); ++i) v.push_back(atof(l[i].c_str())); 
sv.add(k,v);
} 
sv.to_ragged(d,0);
} // sum_vec

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


#if 0
#endif


void list(std::ostream & os)
{

MM_MSG(" list ")
os<<MMPR(m_otus.size()) <<CRLF;
MM_LOOP(ii,m_otus)
{
//MM_LOOP(jj,(*ii).second)
{
//os<<MMPR3((*ii).first, (*jj).first, (*jj).second.size()) <<CRLF;
os<<MMPR3((*ii).first, (*ii).second.name(), (*ii).second.size()) <<CRLF;
}

}

os<<"SampleGroups"<<CRLF;
MM_LOOP(ii,m_sample_groups)
{
//MM_LOOP(jj,(*ii).second)
//{
//os<<MMPR2((*ii).first,(*jj).first)<<CRLF;
os<<MMPR2((*ii).first,(*ii).second.size())<<CRLF;
//}
} // ii 


os<<"Raggeds"<<CRLF;
MM_LOOP(ii,m_ragged_map) { os<<MMPR2((*ii).first,(*ii).second.size())<<CRLF; }//ii/

MM_LOOP(ii,m_fasta_map) { os<<MMPR2((*ii).first,(*ii).second.size())<<CRLF; }//ii/

os<<"BiomFiles"<<CRLF;
MM_LOOP(ii,m_biom_files)
{
//MM_LOOP(jj,(*ii).second)
//{
//os<<MMPR2((*ii).first,(*jj).first)<<CRLF;
auto & bh=(*ii).second;
os<<MMPR4((*ii).first,bh.data().size(),bh.rows().size(),bh.cols().size())<<CRLF;
//}
} // ii 


MM_MSG(" end of list ")
}


void dump_otu(std::ostream & os, const StrTy & comment="")
{
MM_ERR(" dumping otu ... "<<comment)
MM_LOOP(ii,m_otus)
{
// each one is an otu
os<<(*ii).second.dump(0,(*ii).first)<<CRLF;
} // ii 
}

void catalog_otu(std::ostream & os)
{
const StrTy sep=" ";
std::map<StrTy, IdxTy> m;
// for each otu catalog it into m 
MM_LOOP(ii,m_otus) { (*ii).second.catalog(m); } // ii 

MM_LOOP(ii,m)
{
os<<(*ii).first<<sep<<(*ii).second<<CRLF;
}

}

void clear_biom(const StrTy &name)
{
MM_ERR(" clearing "<<MMPR2(name,m_biom_files[name].info()))
if (name.length()!=0) m_biom_files.erase(name); // =x;
else  m_biom_files.clear(); // =x;


}

void dot_sparse_biom(const StrTy &name, const StrTy &fn)
{
sparse_connectivity sc;
BiomFile&  x= m_biom_files[name]; // =x;
BiomFile::data_type & d=x.data();
BiomFile::rows_type & r=x.rows();
BiomFile::cols_type & c=x.cols();
MM_LOOP(ii,d)
{
auto & v=(*ii);
sc.connect(v.i,v.j,v.v);
}
sc.normalize(c.size());
sc.angles(c.size());

}



// (cmd=="parse-biom-json-to-file") parse_biom_json(p1,p2,!true,li.word(3)); } 

void parse_biom_json(const StrTy &name, const StrTy &fn, const bool discard_data=false, const StrTy & dest="")
{
//biom_json x;
//x.parse(tax,atoi(sample.c_str()));
MM_ERR(" parsing json biom "<<MMPR4(name,fn,discard_data,dest))
BiomFile&  x= m_biom_files[name]; // =x;
if (discard_data) x.parse(fn,6);
else if ( dest=="") x.parse(fn,5);
else x.parse(fn,5,dest.c_str());
// the values in x went out of scope...
}

void move_biom(const StrTy &biom, const StrTy &group)
{ 
	copy_biom(biom, group,0);
}

void copy_biom(const StrTy &biom, const StrTy &group, const IdxTy & track_flags)
{
MM_ERR(" copy rows "<<biom)
 copy_biom_row(biom, group,track_flags);
MM_ERR(" copy cols "<<biom)
 copy_biom_col(biom, group);
MM_ERR(" dun copy cols "<<biom)

}





void copy_biom_row(const StrTy &biom, const StrTy &group, const IdxTy track_flags=0)
{	const bool debug=false;
	if (debug) { MM_ERR(" parse_otu "<<MMPR2(biom,group)) } 
	otu_conform oc;
	const bool track_used=((track_flags&1)!=0);
	const bool apply_conform=((track_flags&4)!=0);
	const bool track_canon=((track_flags&8)!=0);
	//const bool track_canon=true;
	const bool use_canon_lut=(m_canon_lut.size()!=0);
	MM_ERR(" copy_biom_row "<<MMPR4(track_used,apply_conform, track_canon,use_canon_lut))
//	const bool track_hier_size=((track_flags&2)!=0);

		//const StrTy & map_name= m_flp.catagory_map();
	//	if (map_name==".") phyla_seq(m);
//		else m=m_catagory_maps[map_name];

	//oc.canon_map(& (m_catagory_maps["genera"] )); 
	if (track_canon)
	{
		 oc.canon_map(& (m_clm)); 
	}
	//const bool skip_old=m_flp.skip_old();
	BiomFile&  x= m_biom_files[biom]; // =x;
	BiomFile::rows_type & r=x.rows();
	auto & otumap= m_otus[group];
	MM_LOOP(ii,r )
	{
		auto & row=(*ii);	
		if (row.id.length()<1){ MM_ERR(" ignoring blank row id")  continue; } 
		//if (debug) { MM_ERR(" fack "<<MMPR2(row.id,row.taxa.size()))} 
		if (debug) { MM_ERR(" fack "<<MMPR2(row.id,row.size()))} 
		if (apply_conform) { 
		const auto &  intaxa=row.taxav(x.tokenizer());
		std::vector<StrTy> conformedtaxa;
		//void zymo_to_vector( Tv & v, Ti & ii, Te & ee)
		oc.zymo_to_vector(conformedtaxa,intaxa.begin(),intaxa.end(),track_flags);
		if ( use_canon_lut ) m_canon_lut.conform_vector(conformedtaxa);
	//	MM_ERR(MMPR(conformedtaxa.size()))
		//otumap.add(row.id,row.taxa);
		otumap.add(row.id,conformedtaxa);
		} 
		//MM_ERR("added ok " )
		else { otumap.add(row.id,row.taxav(x.tokenizer())); } 
	} //ii 
if (track_used) oc.exception_dump(std::cout);
if (track_used) oc.used_dump(std::cout);
} // copy_biom  
// the sample counts are not coming out right 
void copy_biom_col(const StrTy &biom, const StrTy &group)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_otu "<<MMPR2(biom,group)) } 
	//const bool skip_old=m_flp.skip_old();
	BiomFile&  x= m_biom_files[biom]; // =x;
	BiomFile::rows_type & r=x.rows();
	BiomFile::cols_type & c=x.cols();
	// this needs to just make an otu index and then add the data points with index lut
	BiomFile::data_type & d=x.data();
	SampleMap  & samples= m_sample_groups[group];
	const IdxTy oldsamples=samples.size();
	MM_SZ_LOOP(i,d,sz)
	{
		const auto & p=d[i];
		const StrTy & sname=c[p.j].map["id"];
		const IdxTy snametok=samples.st()(sname);
		const StrTy & otun=r[p.i].id;
		const IdxTy otutok=samples.st()(otun);
		samples[snametok].add(otutok,p.v);
		samples[snametok].set_sample_name(sname);
		samples[snametok].set_otu_token(otutok);
	}
		const IdxTy newsamples=samples.size();
		const IdxTy biomsamples=c.size();
		MM_ERR(" sample stats "<<MMPR4(newsamples,oldsamples,(newsamples-oldsamples),biomsamples))
		if ((newsamples-oldsamples)!=biomsamples){
		MM_LOOP(ii,c)
		{
			// this returns a SAMPLE name as the columns aer samples 
			const StrTy & bn=(*ii).map["id"];
			const IdxTy snametok=samples.st()(bn);
			// this sets the OTU NUMER to be consistent ... doh 
			//samples[samples.st()(bn)].set_name(samples.st()(bn));
			samples[snametok].set_sample_name((bn));

		}
		const IdxTy newsamples2=samples.size();
		MM_ERR(" sample stats later "<<MMPR4(newsamples2,oldsamples,(newsamples2-oldsamples),biomsamples))
		}
} // copy_biom  



void parse_otu(const StrTy &name, const StrTy &fn)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_otu "<<MMPR2(name,fn)) } 
	otu_conform oc;
	const bool try_to_conform=!true;
	const IdxTy track_flags=7;
 //	std::vector<StrTy> conformedtaxa;
        //void zymo_to_vector( Tv & v, Ti & ii, Te & ee)
  //      oc.zymo_to_vector(conformedtaxa,intaxa.begin(),intaxa.end(),track_flags);
   //     if ( use_canon_lut ) m_canon_lut.conform_vector(conformedtaxa);



	// FICK CRAP std::ifstream  isn(m_flp.snack_log().c_str());
	std::ifstream  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	//const bool skip_old=m_flp.skip_old();
	//const bool print_dog_days=m_flp.print_dog_days();
	//const bool accumulate_dog_days=m_flp.accumulate_dog_days();
	CommandInterpretter li(is);
	auto & otumap= m_otus[name];
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		const StrTy otu_name=li.word(0);
		//otumap.add(li.line());
		if (debug) { MM_ERR(" fack "<<MMPR((li.line())))} 
		if (try_to_conform)
		{
 			std::vector<StrTy> conformedtaxa;
			// need to remove pfx etc 
			MM_LOOP(ii,(li.words())) {conformedtaxa.push_back(mjm_strings::to_lower(*ii)); } 
//        	oc.zymo_to_vector(conformedtaxa,li.words().begin(),li.words().end(),track_flags);
//			if (conformedtaxa.size()<1) conformedtaxa.push_back(StrTy());
			conformedtaxa[0]=li.words()[0];
			otumap.add(conformedtaxa);

		}
		else otumap.add(li.words());
	} // li.next
} // parse_otu

void dump_reads(std::ostream & os)
{
const StrTy sep=" ";
MM_LOOP(jj,m_sample_groups)
{
// each one is an samle
SampleMap & samples=(*jj).second;
const StrTy & nmm=(*jj).first; 
MM_LOOP(ii,samples)
{ 
//os<<nmm<<sep<<(*ii).second.dump(0,nm+sep+(*ii).first)<<CRLF;
os<<(*ii).second.dump(samples.st(),0,nmm+sep+(samples.st()((*ii).first))); // <<CRLF;

} // ii 
} // jj 


}

// froniers paper format samples on first line 
void parse_front_reads(const StrTy &name, const StrTy &fn)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_reads "<<MMPR2(name,fn)) } 
	// FICK CRAP std::ifstream  isn(m_flp.snack_log().c_str());
	std::ifstream  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	CommandInterpretter li(is);
	std::vector<StrTy> snames; // names;
	SampleMap & samples= m_sample_groups[name];
	// the first line is a bunch of  sample names  
	// firt entry is otu NAME and last is "taxonomy" unassigned 
	// field as a name name  
	if  (li.nextok())
	{
		const IdxTy sz=li.size();
		for(IdxTy i=0; i<(sz-1); ++i) snames.push_back(li.word(i));

	}
	const IdxTy szsamp=snames.size();
	// each following line is a otu  name folled by read counts 
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		//const StrTy  sn=li.word(0);
		const StrTy  otun=li.word(0);
		// TODO the stupid OTU file is not complete, add any stragglers
	const bool kluge_otu_fack=true;
	if (kluge_otu_fack)
	{	
		auto & otumap= m_otus[name];
		std::vector<StrTy> nv(2);
		nv[0]=li.word(0);
		nv[0]=li.word(sz-1);
		if (!otumap.has(nv[0])) {  otumap.add(nv); } 
		}

		if (debug) { MM_ERR(" fack "<<MMPR((li.line())))} 
		//const IdxTy sz=li.words().size();
		const IdxTy szeff=(sz<szsamp)?sz:szsamp;
		for(IdxTy i=1; i<szeff; ++i)
		//samples[names[i]]+=::atof(li.word(i).c_str());
		samples[samples.st()(snames[i])].add(samples.st()(otun),::atof(li.word(i).c_str()));
	////	samples.add(li.words());
	} // li.next
} // parse_otu

// zymo format sequnces in first line 
void parse_reads(const StrTy &name, const StrTy &fn)
{
	const bool debug=false;
	if (debug) { MM_ERR(" parse_reads "<<MMPR2(name,fn)) } 
	// FICK CRAP std::ifstream  isn(m_flp.snack_log().c_str());
	std::ifstream  isn(fn.c_str());
	IsTy * is=& isn; // &std::cin;
	//const bool skip_old=m_flp.skip_old();
	//const bool print_dog_days=m_flp.print_dog_days();
	//const bool accumulate_dog_days=m_flp.accumulate_dog_days();
	CommandInterpretter li(is);
	std::vector<StrTy> names;
	// the zymo format is a line of sequence names
	// follow by sample name and read coutn line for eahc sample 
	SampleMap & samples= m_sample_groups[name];
	// the first line is a bunch of otu or sequence names with first 
	// field as a name name  
	if  (li.nextok())
	{
		const IdxTy sz=li.size();
		for(IdxTy i=0; i<sz; ++i) names.push_back(li.word(i));

	}
	// each following line is a sample name folled by read counts 
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		//MM_ERR(" processing "<<li.dump())
		if (sz<1) continue;
		const StrTy  sn=li.word(0);
		//otumap.add(li.line());
		if (debug) { MM_ERR(" fack "<<MMPR((li.line())))} 
		//const IdxTy sz=li.words().size();
		for(IdxTy i=1; i<sz; ++i)
		//samples[names[i]]+=::atof(li.word(i).c_str());
		samples[samples.st()(sn)].add(samples.st()(names[i]),::atof(li.word(i).c_str()));
	////	samples.add(li.words());
	} // li.next
} // parse_otu





void sum_samples(const StrTy &group, const IdxTy & level)
{


	SampleMap  & samples= m_sample_groups[group];
	// this iterates over all samples within group group 
	MM_LOOP(ii,samples)
	{
//	samples[sn].add(names[i],::atof(li.word(i).c_str()));
	} // ii 
}
// make a single bash or R word out of this 
void make_single(StrTy & f)
{ Ss ss;
//MM_ERR(MMPR(f))
ss<<f;
const StrTy forig=f;
f=StrTy("");
StrTy g=f;
// FICK 
while (!ss.eof() && ss.good()) {g="";  ss>>g; 
//if (!ss.good()) break;
//if (ss.eof()) break; 
if ( f!="" ) f=f+StrTy("_"); 
f=f+g;
 } 
//MM_ERR(MMPR3(forig,f,g))
if (f==StrTy("")) {f="blank"; return ; }
//if (true ) return;
const IdxTy sz=f.length();
if (sz==0) { f="blank"; return; } 
char c[sz+1];
const char * p=f.c_str();
IdxTy pc=0;
for(IdxTy i=0; i<sz; ++i)
{
char d=p[i];
if (d==' ') d='_';
else if (d=='\t') d='_';
else if (d==',') d='_';

c[pc]=d; // p[i];
++pc;
}
c[pc]=0;
f=&c[0];

}
void dump_sample_group(const StrTy &group, std::ostream  & os=std::cout)
{
	IdxTy i=0;
	const StrTy sep=" ";
	SampleMap  & samples= m_sample_groups[group];
	MM_ERR(" dumping "<<MMPR2(group,samples.size()))
	MM_LOOP(ii,samples)
	{
		const StrTy & sn=samples.st()((*ii).first);
		os<<i<<sep<<group<<sep<<sn<<CRLF;
		i=i+1;
	} // ii 

}
void dump_orders(std::ostream & os = std::cout) const
{

MM_ERR(" dump_orders ")
os<<" dump_orders"<<CRLF;
const StrTy sep=" ";
MM_LOOP(ii,m_orders)
{
const StrTy&  name=(*ii).first;
StrTy pfx=name;
//const auto & gr = (*ii).second;
//MM_LOOP(jj,gr)
{
//const StrTy & oname=(*jj).first;
//const StrTy pfx2=pfx+sep+oname; 
//(*jj).second.dump(os,pfx,0);
(*ii).second.dump(os,pfx,0);


} // gr

} // groups




}
void sort_order(const StrTy &name)
{
MM_ERR(" sorting "<<name)
	OrderMap &  order =(name==".")? m_orders[""]:m_orders[name];
//	std::sort(order.begin(),order.end());
	order.sort();
MM_ERR(" done sorting "<<name<<" "<<order.size())
}

// a lits of sample name and catgory names
// note it isi not stored now 
void classify_group(const StrTy &group,const StrTy & fn, const IdxTy & settings=0)
{
	// TODO note that groups not named will have  class of zero which is translated
	// by the string tokenizer to whatever lol 
	SampleMap  & samples= m_sample_groups[group];
	std::ifstream fin(fn);
	CommandInterpretter li(&fin);
	IdxTy i=0; 
	while (li.nextok())
	{
		if (li.size()<2) continue;
		const StrTy & w=li.word(0);
	 	StrTy  w2=li.word(1);
		if (w.length()==0) continue;
		if (w2.length()==0) w2="blankgroup";
		const IdxTy &  p=samples.st()[w];	
		// if the sample is not classified just allow user to output real name later 
		if ( p!=samples.st().badidx()) samples[p].set_class(samples.st()(w2));	
//		else samples[p].set_class(~0);	
		++i;
	} // li 

}
void adhoc_rags(const StrTy &list,const IdxTy & settings,  const StrTy & fn)
{
if (fn.length()!=0 ) { std::ofstream ofs(fn); adhoc_rags(list,settings,ofs); return;}
adhoc_rags(list,settings,std::cout);

}

void adhoc_tree(const StrTy &list,const IdxTy & settings,  const StrTy & fn)
{
if (fn.length()!=0 ) { std::ofstream ofs(fn); adhoc_tree(list,settings,ofs); return;}
adhoc_tree(list,settings,std::cout);

}

void all_samples(Ragged & sglist)
{
		MM_LOOP(ii,m_sample_groups)
		{
			const StrTy group=(*ii).first;
			SampleMap & samples=(*ii).second;
			MM_LOOP(jj,samples)
			{
				const StrTy sample=samples.st()((*jj).first);
				std::vector<StrTy> v;
				v.push_back(group);
				v.push_back(sample);
			//	MM_ERR(" adding "<<MMPR2(group,sample))
				sglist.add(v);
			} // jj 
		} // ii 	
} // all_samples

bool fbit( const IdxTy f, const IdxTy bit) { return (f&(1<<bit))!=0; } 

// can write to rag file but memory tight at this point 
void adhoc_rags(const StrTy &list,const IdxTy & settings,  std::ostream  & os=std::cout)
{
	const IdxTy flags=settings;
	const bool informative_names=fbit(flags,0);	
	const bool split_ambiguous_taxa=fbit(flags,1);
	const bool skip_unclassified=!(fbit(flags,17));
	const bool lump_otus=(fbit(flags,19));
	const bool debug=(fbit(flags,20));
	const bool debug_otu=(fbit(flags,20));
	const bool dump_config_only=(fbit(flags,21));
	const StrTy allowed_otu_name=m_flp.allowed_otu_name(); // ("xlate");
	const StrTy sample_name_class_map=m_flp.adhoc_classify_name(); // ("xlate");
//	const StrTy plot_hash_name=m_flp.plot_hash_name(); // ("xlate");
	const StrTy otu_lumps_name=m_flp.otu_lumps_name(); // ("otu-hash");
	const IdxTy otu_depth_limit=m_flp.otu_depth_limit();

	mjm_taxon_tools mtt;
//	Ragged & allowed_otu_list=m_ragged_map[allowed_otu_name];
//	Ragged::xlate_map otuxlate= allowed_otu_list.xlate_field_map(0,1);

	Ragged & xlate_list=m_ragged_map[sample_name_class_map];
	Ragged::xlate_map mxlate= xlate_list.xlate_field_map(0,1);
	//Ragged::xlate_map grxlate= grclass_list.xlate_field_map(0,0);

	typedef std::map<StrTy,StrTy> Lumps;
	Lumps otu_lumps=m_ragged_map[otu_lumps_name].vkmap(0);

	Ragged & sglist=m_ragged_map[list];
	if (sglist.size()==0){ sglist.load(list,true); }
	if (sglist.size()==0) { all_samples(sglist); } // size==0
	if (debug)
	{
		MM_ERR(MMPR4(list,flags,informative_names,split_ambiguous_taxa)<<MMPR4(skip_unclassified, lump_otus,debug,dump_config_only)<<MMPR4(allowed_otu_name, sample_name_class_map,otu_lumps_name,otu_depth_limit))

	}

	if (dump_config_only)
	{
		//map_dump(sglist,"sglist");
		MM_ERR(sglist.dump());
		map_dump(otu_lumps,"otu_lumps");
		map_dump(mxlate,"mxlate");
//		map_dump(otuxlate,"otuxlate");

		MM_ERR(MMPR4(list,flags,informative_names,split_ambiguous_taxa)<<MMPR4(skip_unclassified, lump_otus,debug,dump_config_only)<<MMPR4(allowed_otu_name, sample_name_class_map,otu_lumps_name,otu_depth_limit))
		return; 
	} // dump_config_only
	typedef mjm_tokenized_vector_map OutMap; 
	OutMap om;
	MM_SZ_LOOP(i,sglist,szlist)
	{
		const auto & w=sglist.line(i);
		if (w.size()<2) continue;
		const StrTy & group = w[0];
//		const bool class_this_group=true; // (grxlate.find(group)!=grxlate.end());
		StrTy  sname = w[1];
		const StrTy & comment =(w.size()>2)? w[2]:StrTy("");
		StrTy snfack=sname;
		auto mi=mxlate.find(sname);
		const bool foundmx=(mi!=mxlate.end());
		if (debug) { MM_ERR(MMPR4(i,group,sname,foundmx))  }
		if (foundmx) snfack=(*mi).second;	
		else if (skip_unclassified ) continue; 

		StrTy sn=group+StrTy("::")+snfack;
		if (comment.length()>0) sn=sn+StrTy(",")+comment;
		SampleMap  & samples= m_sample_groups[group];
		OtuCol  & otumap= m_otus[group];
		const auto ii = samples.find(samples.st()(sname));
		if (ii==samples.end())
		{	
			MM_ERR( " ok, not found "<<MMPR4(sn,sname,snfack,group))
			continue;
		}
		if (debug) { MM_ERR(MMPR4(i,group,sname,snfack))  }
		const Sample & s= (*ii).second;
		MM_LOOP(jj,s)
		{
			const StrTy & otu=samples.st()((*jj).first);
			const Sample::otu_info & otuii=(*jj).second; // .flat();
			const StrTy otukey=otu; // samples.st()(otuii.otu());
			const D hits=otuii.hits();
			typedef std::vector<StrTy> OtuNm;
			//const otu_struct & otui=not_found?(dummy_otu):otumap[otukey]; 
			const otu_struct & otui=otumap[otukey]; 
			OtuNm tv=otui.taxon(7);
			if (informative_names) { tv=mtt.informative_vector(tv); }
			if (lump_otus)
			{
				StrTy nm=StrTy("Other");
				MM_SZ_LOOP(o,tv,tvszz) 
				{ 
				auto it=otu_lumps.find(tv[tvszz-o-1]);
				if (it!=otu_lumps.end()) { nm= (*it).second;  break; }}
				tv.clear();
				tv.push_back(nm);
			}
			if (otu_depth_limit!=0)	
			{
				const IdxTy sz=tv.size();
				if (sz>otu_depth_limit)
				{
					OtuNm tvx(otu_depth_limit);
					for(IdxTy oi=0; oi<otu_depth_limit; ++oi) tvx[oi]=tv[oi];
					tv=tvx;
				
				}

			}
			if (split_ambiguous_taxa)
			{
				auto tvaa=mtt.split(tv);	
				MM_SZ_LOOP(ia,tvaa,tvaasz)
				{
					const auto & tva=tvaa[ia];
//					tvz.add_node_value( tva.t, sn, hits*tva.frac);    
					//otu_select(tva.t,sn,sname,group,hits*tva.frac,otuxlate);
					otu_select(tva.t,sn,sname,group,hits*tva.frac,otu_lumps,om,debug_otu);
				}	
			}
			// have the ignored list here... 
			else { //tvz.add_node_value( tv, sn, hits);   
				//otu_select(tv,sn,sname,group,hits,otuxlate);
				otu_select(tv,sn,sname,group,hits,otu_lumps,om,debug_otu);

			} 
		} // jj sample 
	} // i, sglist 
om.dump(os);
} // adhoc_rags
template<class Tv,class Tx>
void map_dump(Tv & m, const Tx& lbl)
{
	MM_LOOP(ii,m) {MM_ERR(lbl<<" "<<MMPR2((*ii).first,(*ii).second))}
}
template<class Tv, class Tx,class To>
void otu_select(Tv & tv, const StrTy & sn, const StrTy & sname, const StrTy & group, const D & hits, Tx & otuxlate, To & om, const bool debug=true ) 
{
bool ok=false;
StrTy tax="";
StrTy xlate="";
const auto nf =otuxlate.end();
static const StrTy wc="*";
StrTy other="";
auto ii=otuxlate.find(wc);
if ( ii!=nf)  { other=(*ii).second; } 
bool do_other=(other.c_str()[0]!=0);
IdxTy matches=0;
const IdxTy sz=tv.size();
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & n=tv[i];
auto ii=otuxlate.find(n);
// TODO FIXME this must allow multiple hits and
// make normalization work out with "other" ok 
if (ii!=nf){ ++matches;  ok=true; xlate=(*ii).second; } 
tax=tax+StrTy(";")+n;
if (debug) {MM_ERR(MMPR4(i,ok,tax,n))}
} // i 
if (!ok) if ( do_other) { ok=true; xlate=other;  tax=other; } 
if (!ok) { return; } 
if (matches>1)
{
MM_ERR(" need to output all amtches "<<MMPR2(matches,tax)) 
}
//MM_MSG(MMPR4(group,sname,sn,tax)<<MMPR2(xlate,hits))
std::vector<StrTy> k;
k.push_back(group);
k.push_back(sname);
k.push_back(sn);
k.push_back(tax);
k.push_back(xlate);

om.add(k,hits);
} // otu_select


// TODO implement manipulation of sample data prior to getting here 
// 
void adhoc_tree(const StrTy &list,const IdxTy & settings,  std::ostream  & os=std::cout)
{

	const IdxTy flags=settings;
	const bool informative_names=fbit(flags,0);	
	const bool split_ambiguous_taxa=fbit(flags,1);
	const bool write_svg=fbit(flags,2);
	const bool write_latex=fbit(flags,3);
	const bool write_txt=fbit(flags,4);
	const bool write_tax=fbit(flags,5);
	const bool write_svg_ii=fbit(flags,6);
	const bool write_svg_iia=fbit(flags,7);
	const bool write_svg_hm=fbit(flags,8);
	const bool use_sample_name_classes=fbit(flags,16);
	const bool skip_unclassified=!(fbit(flags,17));
	const bool classify_as_group=(fbit(flags,18));
	const bool lump_otus=(fbit(flags,19));
	const StrTy sample_name_class_map=m_flp.adhoc_classify_name(); // ("xlate");
	const StrTy plot_hash_name=m_flp.plot_hash_name(); // ("xlate");
	const StrTy otu_lumps_name=m_flp.otu_lumps_name(); // ("otu-hash");
	const StrTy groups_to_classify=m_flp.groups_to_classify(); // ("xlate");
	const IdxTy otu_depth_limit=m_flp.otu_depth_limit();
	MM_ERR(" enter adhoc_tree "<<MMPR4(flags,list,informative_names,split_ambiguous_taxa)<<MMPR4(write_svg,write_latex,write_txt,write_tax)<<MMPR4(write_svg_ii,use_sample_name_classes,skip_unclassified, classify_as_group)<<MMPR4(sample_name_class_map,groups_to_classify,lump_otus,otu_depth_limit)) 
	MM_ERR(MMPR3(sample_name_class_map,plot_hash_name,groups_to_classify))
	std::map<IdxTy,StrTy> suggested;
	typedef std::map<StrTy,StrTy> TvzParam;
	typedef std::map<StrTy,StrTy> Lumps;
	TvzParam ph;
	Lumps otu_lumps=m_ragged_map[otu_lumps_name].vkmap(0);
	Ragged & phkv=m_ragged_map[plot_hash_name];
	ph= phkv.xlate_field_map(0,1);
	
	otu_struct dummy_otu=otu_struct(StrTy("missing"));
	MM_ONCE(" not using the member tvz", ) 
	mjm_tree_viz kluge ; // m_tvz;
	mjm_tree_viz&  tvz= kluge; // m_tvz;
	mjm_taxon_tools mtt;
	IdxTy lost=0, found=0;
	Ragged & xlate_list=m_ragged_map[sample_name_class_map];
	Ragged & grclass_list=m_ragged_map[groups_to_classify];
	Ragged::xlate_map mxlate= xlate_list.xlate_field_map(0,1);
	Ragged::xlate_map grxlate= grclass_list.xlate_field_map(0,0);
	Ragged & sglist=m_ragged_map[list];
	if (sglist.size()==0){
			 sglist.load(list,true);
		MM_ERR(" facker should have echoed some facking lines")
	}
	// just dump everything. 
	if (sglist.size()==0)
	{
		all_samples(sglist);

	} // size==0
	MM_ERR("adhoc_tree "<<MMPR3(settings,list,sglist.size())<<MMPR4(informative_names,split_ambiguous_taxa,write_svg, write_latex)<<MMPR(write_txt))
	MM_SZ_LOOP(i,sglist,szlist)
	{
		const auto & w=sglist.line(i);
		if (w.size()<2) continue;
		const StrTy & group = w[0];
		const bool class_this_group=(grxlate.find(group)!=grxlate.end());
		StrTy  sname = w[1];
		const StrTy & comment =(w.size()>2)? w[2]:StrTy("");
		StrTy snfack=sname;
		if (classify_as_group) { snfack="groupavg"; } 
		else if (use_sample_name_classes&&class_this_group)
		{
			auto mi=mxlate.find(sname);
			//MM_ERR(" xlate "<<MMPR2(mxlate.size(),sname))
			if (mi!=mxlate.end()) snfack=(*mi).second;	
			else if (skip_unclassified ) continue; 
		}
		StrTy sn=group+StrTy("::")+snfack;
		if (comment.length()>0) sn=sn+StrTy(",")+comment;
		suggested[i]=sn;
		SampleMap  & samples= m_sample_groups[group];
		OtuCol  & otumap= m_otus[group];
		const auto ii = samples.find(samples.st()(sname));
	//	MM_ERR(" sample in fo "<<MMPR4(sn,(samples.st()(sname)),group,sname))
		if (ii==samples.end())
		{	
			MM_ERR( " not found "<<MMPR4(sn,sname,snfack,group))
			continue;
		}
		const Sample & s= (*ii).second;
		MM_LOOP(jj,s)
		{
			const StrTy & otu=samples.st()((*jj).first);
			const Sample::otu_info & otuii=(*jj).second; // .flat();
			const StrTy otukey=otu; // samples.st()(otuii.otu());
			const bool not_found=(otumap.find(otukey)==otumap.end());
			if (not_found)
			{ //// TODO FIXME wtf
				const StrTy otukey2= samples.st()(otuii.otu());
				MM_ONCE(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2),)	
				++lost;	
			}
			else ++found;
//otu_struct os= otu_struct(name);
			const otu_struct & otui=not_found?(dummy_otu):otumap[otukey]; 
			const D hits=otuii.hits();
			typedef std::vector<StrTy> OtuNm;
			OtuNm tv=otui.taxon(7);
			if (informative_names)
			{
				tv=mtt.informative_vector(tv);
			}
			if (lump_otus)
			{
				StrTy nm=StrTy("Other");
				MM_SZ_LOOP(o,tv,tvszz) 
				{ 
				auto it=otu_lumps.find(tv[tvszz-o-1]);
				if (it!=otu_lumps.end()) { nm= (*it).second;  break; }}
				tv.clear();
				tv.push_back(nm);
			}
			if (otu_depth_limit!=0)	
			{
				const IdxTy sz=tv.size();
				if (sz>otu_depth_limit)
				{
					OtuNm tvx(otu_depth_limit);
					for(IdxTy oi=0; oi<otu_depth_limit; ++oi) tvx[oi]=tv[oi];
					tv=tvx;
				
				}

			}
			if (split_ambiguous_taxa)
			{
				auto tvaa=mtt.split(tv);	
				MM_SZ_LOOP(ia,tvaa,tvaasz)
				{
					const auto & tva=tvaa[ia];
					tvz.add_node_value( tva.t, sn, hits*tva.frac);    
				}	
			}
			// have the ignored list here... 
			else { tvz.add_node_value( tv, sn, hits);   } 
		} // jj
	} // i
	tvz.initial_sample_order(suggested);	
 	if (write_svg) tvz.write_svg(os,flags); 
 	if (write_svg_ii) tvz.write_svg_ii(os,0,m_pheno_notes,ph); 
 	if (write_svg_iia) tvz.write_svg_ii(os,1,m_pheno_notes,ph); 
 	//if (write_svg_hm) tvz.write_svg_ii(os,2,m_pheno_notes,ph); 
 	if (write_svg_hm) tvz.write_svg_ii(os,2,ph); 
 	if (write_latex) tvz.write_latex(os,flags); 
 	if (write_txt) tvz.write_txt(os,flags); 
 	if (write_tax) tvz.write_tax_tree(os,flags); 
	MM_ERR(" done adhoc_tree "<<MMPR2(lost,found))
} // adhoc_tree

		typedef std::map<StrTy, D> HitAcc;

void dump_select(const StrTy &group,const IdxTy & settings,  const StrTy & fn)
{
// this needs to find the leng
//if (fn.length()==0) 
if (fn.c_str()[0]==0) { dump_select(group,settings); return; }
MM_ERR(" dump_select "<<MMPR3(group,settings,fn))
std::ofstream fos(fn);
dump_select(group,settings,fos);

}
void dump_select(const StrTy &group,const IdxTy & settings,  std::ostream  & os=std::cout)
{
	SampleMap  & samples= m_sample_groups[group];
	OtuCol  & otumap= m_otus[group];
	otu_struct dummy_otu=otu_struct(StrTy("missing"));
	const bool use_hitacc=((settings&1)==0);
	const bool use_tvz=((settings&2)!=0);
	const StrTy sep=" ";
	mjm_tree_viz tvz;
	typedef OrderMap  Torder;
	IdxTy sampi=0;
	const IdxTy nsamp=samples.size();
	//typedef std::map<StrTy, IdxTy > Tlut;
	Torder&  order =m_orders[""];
	IdxTy lost=0, found=0;
	MM_ERR(" dump_select "<<MMPR3(settings,use_hitacc,use_tvz)<<MMPR4(group,samples.size(),otumap.size(),order.size()))
	MM_LOOP(ii,samples)
	{
		Ss ss;
		HitAcc hitacc;
		D hittot=0;
		const IdxTy wtf=(*ii).first;
		const StrTy & sn=samples.st()(wtf);
		if (sn.length()==0)
		{ MM_ERR(" sample name is blank ignore ") continue; }
		const Sample & s= (*ii).second;
		const bool include_sample=true;
		if ( include_sample) { 
		MM_LOOP(jj,s)
		{
			const StrTy & otu=samples.st()((*jj).first);
			const Sample::otu_info & otuii=(*jj).second; // .flat();
			const StrTy otukey=otu; // samples.st()(otuii.otu());
			const StrTy otukey2= samples.st()(otuii.otu());
			const bool not_found=(otumap.find(otukey)==otumap.end());
			if (not_found)
			{ //// TODO FIXME wtf
			if (false) 	MM_ERR(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2))	
				MM_ONCE(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2),)	
				++lost;	
			}
			else ++found;
//otu_struct os= otu_struct(name);
			const otu_struct & otui=not_found?(dummy_otu):otumap[otukey]; 
			const D hits=otuii.hits();
			std::vector<StrTy> tv=otui.taxon(7);

			//void add_node_value( const TaxVec & v, const StrTy & sample, const D & val)
			if (use_tvz) { tvz.add_node_value( tv, sn, hits); } 
			if (use_hitacc) { 
			MM_SZ_LOOP(i,tv,tvsz)
			{
//				MM_ERR(" compare "<<MMPR3(i,tv[i],order(tv[i])))
				if (order(tv[i])!=order.bad()) 
				{hitacc[tv[i]]+=hits; // break;
 				}
			} // tv
			} // use_hitacc
			hittot+=hits;
		} // jj
		// TODO FIXME the order thing is probably not right,
		// need to check this u look at outpu at elast 
		if (use_hitacc) { 
		const IdxTy sclass=s.get_class();
		const bool valid_class=(sclass!=bad());
		const StrTy&  clv=((valid_class?samples.st()(sclass):"noclass"));
		ss<<group<<sep<<sn<<sep<<clv;
		MM_SZ_LOOP(kk,order,orsz)
		{
			const StrTy &key=order(kk); // (*kk).first;
			if (key.c_str()[0]==0)
			{ 
				MM_ONCE(" warnin bkanls zero in order "<<MMPR2(kk,orsz),)
				continue;
 } 
			D v=0;
			auto mf=hitacc.find(key);
			if (mf!=hitacc.end()) v= (*mf).second;
			if (hittot!=0) v=v/hittot;
			ss<<sep<<key<<sep<<v;
		}
		os<<ss.str()<<CRLF;
		} // use_hitacc
		if (false) if (use_tvz) if (sampi>20)
		{
			MM_ONCE(" only doing 20 samples for tvz",)
			break;
		}
		} // include_sample 
		++sampi;
		if ((sampi&31)==0) {
			 MM_STATUS(" dump_select "<<MMPR4(sampi,nsamp,hitacc.size(),hittot)) }
	} // ii 
	if ( use_tvz) {tvz.write_svg(os,0); }
MM_ERR(MMPR2(lost,found))
} // dump_select


/*
 Convert 1 or more biom files into ssv seq,sample,count,otu0,otu1.... 
2544  ./mjm_biom_to_rag.out  "b2r b2.biom 0 b2.ssv" quit
 2545  ./mjm_biom_to_rag.out  "b2r b3.biom 0 b3.ssv" quit
 2546  . srecord
Sum counts by sample and n key fields
 2559  ./mjm_rag_reduce.out "reduce s3.ssv b3.ssv S0S2K4-9K1" quit 2>&1 
 2560  ./mjm_rag_reduce.out "reduce s2.ssv b2.ssv S0S2K4-9K1" quit 2>&1 
 2561  ./mjm_rag_reduce.out "reduce s1.ssv b1.ssv S0S2K4-9K1" quit 2>&1 

Put all sample counts on one line followed by otu0,otu1... 

 2562  ./mjm_rag_reduce.out "merge m13.ssv s1.ssv s3.ssv" quit 2>&1 
 2563  more m13.ssv | grep Porp


*/



typedef std::vector<StrTy> BiomList;
typedef std::vector<Ragged> ReducedList;
bool Bit(const IdxTy flag, const IdxTy bit) const {return (0!=((flag>>bit)&1)); } 
// called with use flags>>8
IdxTy reconcile(ReducedList & rl, const BiomList & bl, const IdxTy flags) 
{
IdxTy rc=0;
// the sequid is not here FICK 
const IdxTy seqcol=3;
const IdxTy blsz=bl.size();
MM_ERR(MMPR(m_recon.dump_eq()))
for(IdxTy pass=0; pass<2; ++pass)
{
for(IdxTy i=0; i<blsz; ++i)
{
Ragged & r=rl[i];
MM_LOOP(jj,r)
{
Line & l=(*jj);
const StrTy & seqid=l[seqcol];
//MM_ERR(MMPR(seqid))
// AssTy reconcile(const StrTy & f, const StrTy & seqno, const AssTy & ass, const IdxTy flags)
Recon::assign_type ai,af;
for(IdxTy j=seqcol+1; j<l.size(); ++j) ai.push_back(l[j]); 
if (pass==0) af=m_recon.reconcile(bl[i],seqid,ai ,flags);
else af=m_recon.lookup(bl[i],seqid,ai ,flags);
// only needed on second pass really... 
if (af!=ai)
{
{
Ss ss;
MM_LOOP(kk,af) { ss<<(*kk)<<":"; }
Ss rr;
MM_LOOP(kk,ai) { rr<<(*kk)<<":"; }
MM_ERR( " reconciling "<<MMPR3(r.name(),ss.str(),rr.str()))
}
MM_ERR(MMPR2(seqcol,af.size()))
Line l2(seqcol+1+af.size());
for(IdxTy j=0; j<=seqcol; ++j) l2[j]=l[j]; 
for(IdxTy j=seqcol+1; j<l2.size(); ++j) l2[j]=af[j-seqcol-1]; 
l=l2;
}

} // jj

} // i 
} // pass 
return rc;
} // reconcile 

IdxTy bioms_to_ssv(Ragged & rdest, const BiomList & bl, const IdxTy redlev, const IdxTy flags)
{
const bool dump_ssv=Bit(flags,0);
const bool dump_reduced=Bit(flags,1);
// this has to match merge usage... 
const bool include_header=Bit(flags,2);
const bool do_reconcile=Bit(flags,3);
MM_ERR(MMPR4(bl.size(),redlev,__FUNCTION__,do_reconcile)<<MMPR3(dump_ssv,dump_reduced,include_header))
if (!do_reconcile) return bioms_to_ssv_old(rdest,bl,redlev,flags); 
IdxTy rc=0;
BiomToRag btor;
RagReduce rred;
//Ragged rdest;
const IdxTy blsz=bl.size();
//std::vector<Ragged> reduceds(blsz);
ReducedList expands(blsz);
ReducedList reduceds(blsz);
//const StrTy redstr="S0S2K4-9K1"; 
Ss ss; ss<<"S0S2K4-"; ss<<(10-redlev); ss<<"K1";
Ss rr;  rr<<redlev; 

const StrTy redstr=ss.str(); // "S0S2K4-9K1"; 
const StrTy redlevs=rr.str();
MM_ERR(MMPR2(redstr,redlevs))
IdxTy i=0;
MM_LOOP(ii,bl)
{
	const StrTy & fn=(*ii);
	MM_ERR(MMPR(fn))
//	Ragged r,d; 
	Ragged r; 
 	btor.load_ragged(r,fn,0);
	r.name(fn);
	btor.parse_to_ragged(expands[i],r,0);
	if (dump_ssv)
	{
		MM_MSG(expands[i].dump());
		const StrTy fnout=fn+".ssv";
		std::ofstream ofs(fnout);
		ofs<<expands[i].dump_ssv();
	}
++i;
} // ii 
// these still contain sequence numbers... 
if (do_reconcile) reconcile(expands,bl,flags>>8);
i=0;
MM_LOOP(ii,bl)
{

 //2561  ./mjm_rag_reduce.out "reduce s1.ssv b1.ssv S0S2K4-9K1" quit 2>&1 
// add reconciler to merge sequences not taxo ass 
	IdxTy rc=rred.reduce(reduceds[i],expands[i],redstr,0);
	//reduceds[i].name(r.name()); 
	reduceds[i].name((*ii)); 
	if (dump_reduced)
	{
		MM_ERR(reduceds[i].dump())
		const StrTy fnout=(*ii)+".reduced"+redlevs;
		std::ofstream ofs(fnout);
		ofs<<reduceds[i].dump_ssv();
	}

++i;
} // ii 
// at this point, reduced[] contained the expanded text
// as sample, count, sequence, tax[]
// and two pass recon can be done now...
// this is too late, they have been added and sequid no meaning 
//if (do_reconcile) reconcile(reduceds,bl,0);
// this needs all the files at once to get sample count.. 
	IdxTy rcmerge=rred.merge(rdest,reduceds,flags);

return rc; 
} // bioms_to_ssv

// facked up reconcile no sequid 
IdxTy bioms_to_ssv_old (Ragged & rdest, const BiomList & bl, const IdxTy redlev, const IdxTy flags)
{
const bool dump_ssv=Bit(flags,0);
const bool dump_reduced=Bit(flags,1);
// this has to match merge usage... 
const bool include_header=Bit(flags,2);
// code is in wrong palce 
const bool do_reconcile=false; // Bit(flags,3);
MM_ERR(MMPR4(bl.size(),redlev,__FUNCTION__,do_reconcile)<<MMPR3(dump_ssv,dump_reduced,include_header))
IdxTy rc=0;
BiomToRag btor;
RagReduce rred;
//Ragged rdest;
const IdxTy blsz=bl.size();
//std::vector<Ragged> reduceds(blsz);
ReducedList reduceds(blsz);
//const StrTy redstr="S0S2K4-9K1"; 
Ss ss; ss<<"S0S2K4-"; ss<<(10-redlev); ss<<"K1";
Ss rr;  rr<<redlev; 

const StrTy redstr=ss.str(); // "S0S2K4-9K1"; 
const StrTy redlevs=rr.str();
MM_ERR(MMPR2(redstr,redlevs))
IdxTy i=0;
MM_LOOP(ii,bl)
{
	const StrTy & fn=(*ii);
	MM_ERR(MMPR(fn))
	Ragged r,d; 
 	btor.load_ragged(r,fn,0);
	r.name(fn);
	btor.parse_to_ragged(d,r,0);
	if (dump_ssv)
	{
		MM_MSG(d.dump());
		const StrTy fnout=fn+".ssv";
		std::ofstream ofs(fnout);
		ofs<<d.dump_ssv();
	}
 //2561  ./mjm_rag_reduce.out "reduce s1.ssv b1.ssv S0S2K4-9K1" quit 2>&1 
// add reconciler to merge sequences not taxo ass 
	IdxTy rc=rred.reduce(reduceds[i],d,redstr,0);
	reduceds[i].name(r.name()); 
	if (dump_reduced)
	{
		MM_ERR(reduceds[i].dump())
		const StrTy fnout=fn+".reduced"+redlevs;
		std::ofstream ofs(fnout);
		ofs<<reduceds[i].dump_ssv();
	}

++i;
} // ii 
// at this point, reduced[] contained the expanded text
// as sample, count, sequence, tax[]
// and two pass recon can be done now...
// this is too late, they have been added and sequid no meaning 
if (do_reconcile) reconcile(reduceds,bl,0);
// this needs all the files at once to get sample count.. 
	IdxTy rcmerge=rred.merge(rdest,reduceds,flags);

return rc; 
} // bioms_to_ssv_old



IdxTy cmd_bioms_to_ssv( const CommandInterpretterParam  &  cip)
{
IdxTy rc=0;
//IdxTy bioms_to_ssv(const BiomList & bl, const IdxTy redlev, const IdxTy flags)
//const bool dump_ssv=Bit(flags,0);
//const bool dump_reduced=Bit(flags,1);
const StrTy fout=cip.p1;
const StrTy dname=cip.p2;
//const IdxTy redlev=myatoi(cip.p2);
const IdxTy redlev=myatoi(cip.wif(3));
// change dest to named... 
const IdxTy flags=myatoi(cip.wif(4));
const IdxTy fn1=5;
BiomList  bl;
for(IdxTy i=fn1; cip.wif(i).length()!=0; ++i) bl.push_back(cip.wif(i)); 
// the trimmer expects shorter ones first which should be old
std::sort(bl.begin(),bl.end());

//Ragged r;
Ragged&  r=m_ragged_map[dname];

{
Ss ss; MM_LOOP(ii,bl) ss<<(*ii)<<" ";
MM_ERR(MMPR4(__FUNCTION__,fout,dname,r.size())<<MMPR3(redlev,flags,ss.str()))
}

// start passign flags>> 8 ? 
rc=bioms_to_ssv(r,bl,redlev,flags);
if (fout.length())
{
std::ofstream fos(fout);
fos<<r.dump_ssv();
}
else MM_ERR(" not saving since fout is empty "<<MMPR2(__FUNCTION__,fout))

return rc;
} // cmd_bioms_to_ssv


// outpu reconciliation data such as sequence assignment fasta etc 
IdxTy cmd_reconciled( const CommandInterpretterParam  &  cip)
{
IdxTy rc=0;
//IdxTy bioms_to_ssv(const BiomList & bl, const IdxTy redlev, const IdxTy flags)
//const bool dump_ssv=Bit(flags,0);
//const bool dump_reduced=Bit(flags,1);
const StrTy fout=cip.p1;
const StrTy dname=cip.wif(3);
//const IdxTy redlev=myatoi(cip.p2);
//const IdxTy redlev=myatoi(cip.wif(3));
// change dest to named... 
const IdxTy flags=myatoi(cip.wif(2));
const bool just_ass_fasta=(flags==0);
MM_ERR(MMPR4(fout,dname,flags,just_ass_fasta))
//const IdxTy fn1=5;
//Ragged r;
// start passign flags>> 8 ? 
if (just_ass_fasta)
{
if (fout!="-")
{
std::ofstream fos(fout);
fos<< m_recon.dump_ass_fasta(flags);
}
else
{
std::cout<< m_recon.dump_ass_fasta(flags);
}
return rc;
}
Ragged&  r=m_ragged_map[dname];
rc=m_recon.dump_fasta(r,flags);
if (fout.length())
{
std::ofstream fos(fout);
// this changes spaces to underscore? Should escap maybe
fos<<r.dump_ssv();
}
else MM_ERR(" not saving since fout is empty "<<MMPR2(__FUNCTION__,fout))

return rc;
} // cmd_reconciled




// supposed to make plotting easier.... 
// this needs to pick two integer counts and two keys for
// label and "color" such as a genus and phylum 
// output is sorted by log ratio. Exceptionas
// are zero in one count but none if bothzer.
// need to make a ratio for the thexceptions 
// for uitable plotting. 
class _count_ratio_params
{
public:

 _count_ratio_params() { Init(); } 
 _count_ratio_params(const Logic & flp) {  
Init(); 
have_header=flp.crp_have_header(); // false;
include_header=flp.crp_include_header(); // false;
fix_offset=flp.crp_fix_offset(); // !true;
add_abundances=flp.crp_add_abundances(); // !true;
add_identicals=flp.crp_add_identicals(); // !true;
saf=flp.crp_saf(); // 3;
osize=flp.crp_osize(); // 2+2+1;
n1=flp.crp_n1(); // 1;
n2=flp.crp_n2(); // 2;
ngroup=flp.crp_ngroup(); // 5;
p= flp.crp_p(); // -1; 



} 


void Init()
{
have_header=false;
include_header=false;
fix_offset=!true;
add_abundances=!true;
add_identicals=!true;
saf=3;
osize=2+2+1;
n1=1;
n2=2;
ngroup=5;
p=-1; 

} // Init 

StrTy dump() const
{
Ss ss;
ss<<MMPR4(have_header, include_header, fix_offset, saf);
ss<<MMPR4(osize,n1,n2,ngroup);
ss<<MMPR3(add_identicals,add_abundances,p); 
return ss.str();
} // dump

 bool have_header; // =false;
 bool include_header; // =false;
 bool fix_offset; // =!true;
 bool add_abundances; // =!true;
 bool add_identicals; // =!true;
 IdxTy saf; // =3;
 int osize; // =2+2+1;
 int n1; // =1;
 int n2; // =2;
 int ngroup; // =5;
 int p; // =-1; 


}; // _count_ratio_params
typedef _count_ratio_params count_ratio_param_type;

IdxTy count_ratios(Ragged & d,const Ragged & _s,const count_ratio_param_type & crp, const IdxTy flags)
{
IdxTy rc=0;
Ragged s=_s;
const bool have_header=crp.have_header; // false;
const bool include_header=crp.include_header; // false;
const bool fix_offset=crp.fix_offset; // !true;
const bool add_abundances=crp.add_abundances;
const bool add_identicals=crp.add_identicals; // true;
const IdxTy fix_fp_flags=1;
const IdxTy saf=crp.saf; // 3;
const int osize=crp.osize; // 2+2+1;
const int n1=crp.n1; // 1;
const int n2=crp.n2; // 2;
const int ngroup=crp.ngroup; // 5;
const int p=crp.p; // -1; 

// in some cases, specie and genus are the same but things like 
// family changed in pipeline... in particular ventriculi 
// this should be fixed in reconcile 
// this needs to put them all together and then make ok and zed's again
// doh.. 
const IdxTy istart=have_header?1:0;
if (add_identicals)
{
Ragged q;
typedef std::map< std::vector<StrTy>, IdxTy> Pp;
Pp pp;
const IdxTy szs=s.size();
if (istart) if (szs) q.add(s[0]);
for(IdxTy i=istart;i<szs; ++i)
{
const Line & l=s[i];
const int len=l.size();
if ((len<=n1) ||(len<=n2)) { q.add(l); continue; } 

const StrTy k1=((len+p)>0)?l[len+p]:StrTy("-");
const StrTy k2=(len>ngroup)?l[ngroup]:StrTy("-");
if ((k1=="-") ||(k2=="-")) { q.add(l); continue; } 
std::vector<StrTy> kk;
kk.push_back(k1);
kk.push_back(k2);
const auto ii=pp.find(kk);
if (ii==pp.end()){   pp[kk]=q.size(); q.add(l); continue; } 
const IdxTy idx=(*ii).second;
// just keep the first one and add the two important columns.. 
const IdxTy c1=(len>n1)?myatoi(l[n1]):0;
const IdxTy c2=(len>n2)?myatoi(l[n2]):0;
const IdxTy qlen=q[idx].size();
const IdxTy c1q=(qlen>n1)?myatoi(q[idx][n1]):0;
const IdxTy c2q=(qlen>n2)?myatoi(q[idx][n2]):0;
Line  l2=s[i];
{ Ss ss; ss<<(c1+c1q); l2[n1]=ss.str(); } 
{ Ss ss; ss<<(c2+c2q); l2[n2]=ss.str(); } 
q[idx]=l2;

} // i

s=q;

} // add_identicals

int lenlim=n1;
if (n2>lenlim) lenlim=n2;
// this is not good for the NA/ stuff do not want to exclude.. 
if (ngroup>lenlim) lenlim=ngroup;

if ((-p)>lenlim) lenlim=(-p);
int sn1=0, sn2=0, max1=0,max2=0,min1=0, min2=0;
Lines ok,zed1,zed2; // ,except;
const IdxTy szs=s.size();
for(IdxTy i=istart;i<szs; ++i)
{
const Line & l=s[i];
const int len=l.size();
// I'm not sure these were being dropped but no idea how included
// lol 
//if (len<lenlim) { except.push_back(l); continue; } // lenlim 
Line lout(osize);
const IdxTy c1=(len>n1)?myatoi(l[n1]):0;
const IdxTy c2=(len>n2)?myatoi(l[n2]):0;
sn1+=c1; sn2+=c2;
lout[0]=(len>n1)?l[n1]:StrTy("0");
lout[1]=(len>n1)?l[n2]:StrTy("0");
// in some cases these seem to be blank confusing R.. 
lout[3]=((len+p)>0)?l[len+p]:StrTy("-");
lout[4]=(len>ngroup)?l[ngroup]:StrTy("-");

const bool c1z=(c1==0);
const bool c2z=(c2==0);
if (c1z&&c2z) continue;
if(!c1z){if (min1==0) min1=c1; if (min1>c1) min1=c1; if ( c1>max1) max1=c1;}
if(!c2z){if (min2==0) min2=c2; if (min2>c2) min2=c2; if ( c2>max2) max2=c2;}

if (c1z) {lout[2]=l[n2]; zed1.push_back(lout); continue; }  
if (c2z) {lout[2]=l[n1]; zed2.push_back(lout); continue; }  
//Ss ss; ss<<(D(c1)/D(c2));
//lout[2]=ss.str(); 
lout[2]=fix_fp_fmt((D(c1)/D(c2)),fix_fp_flags);
ok.push_back(lout);
} // i 
// this is in the wrong place it needs to work on each
// read or sequence doh 
if (false&&fix_offset)
{
MM_LOOP(ii,ok)
{
Line & l=(*ii);
const D c1=myatoi(l[0]);
const D c2=myatoi(l[1]);
Ss ss; ss<<((c1+2*min1)/(c2+2*min2));
l[2]=ss.str(); 
} // ii 
} // fix_offset
sort_add(d,zed1,1e-2*D(sn2),saf);
sort_add(d,ok,D(sn1)/D(sn2),saf);
sort_add(d,zed2,1e-2*D(sn1),saf);
if( add_abundances )
{
Ragged dx; // problems with refs and iterators...
const D f1=1.0/D(sn1);
const D f2=1.0/D(sn2);

MM_LOOP(ii,d)
{
Line l=(*ii); // modify this one...  
l.push_back(div(l[0],f1,fix_fp_flags));
l.push_back(div(l[1],f2,fix_fp_flags));

dx.push_back(l);
} // dx
d=dx; // wtf 
} // add_abundances

return rc; 
} // count_ratios

StrTy fix_fp_fmt(const D & v, const IdxTy flags )
{
const bool def=Bit(flags,0);
Ss ss;
if (!def)
{
ss<<std::setprecision(3);
if (v<.1) ss<<std::scientific;
} // def
ss<<v;
return ss.str();
} // fix_fp_fmt;
// this is really dumb as the denominator can be done
// just once and even turned into a multiply... 
//StrTy div(const StrTy & n, const StrTy & d) 
//{ return atof(n.c_str())/atof(d.c_str()); } // div

StrTy div(const StrTy & n, const D & f,const IdxTy flags) 
//{ Ss ss; ss<<( f*atof(n.c_str())) ; return ss.str(); } // div
{return fix_fp_fmt( f*atof(n.c_str()),flags); } // div


void sort_add(Ragged & d, Lines & x, const D& norm,const IdxTy flags)
{
const bool normalize=Bit(flags,0);
const bool fill=Bit(flags,1);
//std::sort(idxs.begin(), idxs.end(), [&](const IdxTy& a, const IdxTy& b) { return merit[a]>merit[b]; });
std::sort(x.begin(),x.end(), [&](const Line& a, const Line& b) { return Myt::merit(a)>Myt::merit(b); });
MM_LOOP(ii,x) { 
if (normalize)
{
D v=merit(*ii);
Ss ss; ss<<(v/norm); (*ii)[2]=ss.str(); 
} // normalize 
if (fill)
{// should be char for char non-white check 
MM_LOOP(jj,(*ii))
{
StrTy  x=(*jj);
bool blank=true;
const char * xp=x.c_str();
while (*xp) {if (*xp!=' ') { blank=false; break; }     ++xp;} 
if (blank) (*jj)="blank";
} // jj 
} // fill 
d.push_back(*ii); } 
//MM_LOOP(ii,ok) { d.push_back(*ii); } 
//MM_LOOP(ii,zed2) { d.push_back(*ii); } 

} // sort 

static D merit(const Line & a) 
{ return atof(a[2].c_str()); }

IdxTy count_key_select(Ragged & d,const Ragged & s,const count_ratio_param_type & crp, const IdxTy flags)
{
IdxTy rc=0;
return count_ratios(d,s,crp,flags);


//std::sort(idxs.begin(), idxs.end(), [&](const IdxTy& a, const IdxTy& b) { return merit[a]>merit[b]; });

return rc;
} // count_key_select

IdxTy cmd_count_key_select( const CommandInterpretterParam  &  cip)
{
IdxTy rc=0;
//const bool dump_ssv=Bit(flags,0);
//const bool dump_reduced=Bit(flags,1);
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
//const IdxTy redlev=myatoi(cip.p2);
const IdxTy flags=myatoi(cip.wif(3));
const bool load_s=Bit(flags,0);
const bool save_d=Bit(flags,1);
count_ratio_param_type crp(m_flp); 
Ragged & s=m_ragged_map[sname];
Ragged & d=m_ragged_map[dname];
if (load_s) s.load(sname);

MM_ERR(MMPR4(__FUNCTION__,sname,dname,d.size())<<MMPR2(s.size(),flags)<<MMPR3(load_s,save_d,crp.dump()))

rc=count_key_select(d,s,crp, flags>>2);
if (save_d) { std::ofstream fos(dname); fos<<d.dump_ssv(); }


return rc;
} //cmd_count_key_select

void _quiet() { mjm_global_flags::mm_err_enable=!true; }
void _nquiet() { mjm_global_flags::mm_err_enable=true; }

IdxTy cmd_load_recon( const CommandInterpretterParam  &  cip)
{
IdxTy rc=0;
// load into fasta named f from file fn 
const StrTy f=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_recon.load(f,fn,flags);
return 0;
} // cmd_load_recon
IdxTy cmd_find_olaps( const CommandInterpretterParam  &  cip)
{
IdxTy rc=0;
//IdxTy bioms_to_ssv(const BiomList & bl, const IdxTy redlev, const IdxTy flags)
//const bool dump_ssv=Bit(flags,0);
//const bool dump_reduced=Bit(flags,1);
const StrTy fout=cip.p1;
const StrTy dname=cip.p2;
//const IdxTy redlev=myatoi(cip.p2);
const IdxTy redlev=myatoi(cip.wif(3));
// change dest to named... 
const IdxTy flags=myatoi(cip.wif(4));
const IdxTy fn1=5;
BiomList  bl;
for(IdxTy i=fn1; cip.wif(i).length()!=0; ++i) bl.push_back(cip.wif(i)); 
// the trimmer expects shorter ones first which should be old
std::sort(bl.begin(),bl.end());

//Ragged r;
Ragged&  r=m_ragged_map[dname];

{
Ss ss; MM_LOOP(ii,bl) ss<<(*ii)<<" ";
MM_ERR(MMPR4(__FUNCTION__,fout,dname,r.size())<<MMPR3(redlev,flags,ss.str()))
}

// start passign flags>> 8 ? 
rc=m_recon.find_olaps(r,bl,redlev,flags);
//rc=bioms_to_ssv(r,bl,redlev,flags);
if (fout.length())
{
std::ofstream fos(fout);
fos<<r.dump_ssv();
}
else MM_ERR(" not saving since fout is empty "<<MMPR2(__FUNCTION__,fout))


return rc;
} // cmd_find_olaps



static const IdxTy & bad() { static const IdxTy i=~0U; return i; } 

enum { RAW=1, COLLATE=2, UPDATE_ORDER=4, USE_MAP=8,OUTPUT_CLASS=16, GS_MAP=32, TAX_TREE=64, DEBUG_INFO=(1<<16)};

// The main dump and collate routine 
//void dump_group(const StrTy &group,const bool do_raw, const bool do_collate, const bool use_map=true, std::ostream  & os=std::cout)
void dump_group(const StrTy &group,const IdxTy & settings,  std::ostream  & os=std::cout)
{
enum { RAW=1, COLLATE=2, UPDATE_ORDER=4, USE_MAP=8,OUTPUT_CLASS=16, GS_MAP=32, TAX_TREE=64, DEBUG_INFO=(1<<16)};
	// dump as encountered 
	const bool do_raw= ((settings&RAW)!=0);
	// collate the entries before dumping 
	const bool do_collate= ((settings&COLLATE)!=0);
	// and output order map is maintained so a consistent cumulative value can be created
	// and something needs to be done with new taxa types 
	const bool update_order= ((settings&UPDATE_ORDER)!=0);
	//  map the taxa into catagories  
	const bool use_map= ((settings&USE_MAP)!=0);
	// instead of labelling with sample name, label with class to which sample belongs
	// or varius actions if no class found either output sample name or default class 
	const bool output_class= ((settings&OUTPUT_CLASS)!=0);
	const bool use_gs_map= ((settings&GS_MAP)!=0);
	const bool use_tax_tree= ((settings&TAX_TREE)!=0);

	const bool debug=((settings&(DEBUG_INFO))!=0); 
	IdxTy flatten_scheme=0;		
	if (use_map) flatten_scheme=0;
	else if ( use_gs_map) flatten_scheme=1;
	else if (use_tax_tree) flatten_scheme =3;
	else flatten_scheme =2;
	MM_ERR(" dump_group "<<MMPR3(settings,debug,flatten_scheme)<<MMPR4(do_raw,do_collate,update_order,use_map)<<MMPR2(output_class,use_gs_map))



	const StrTy sep=" ";
	const IdxTy otu_format=m_flp.otu_format();
	const bool filter_samples_from_map=m_flp.use_sample_filter();
	const StrTy & sfname= m_flp.sample_filter();
	const StrTy gs_mapname="blank";
	GSMap & gs_map= m_gs_maps[gs_mapname];
	SampleMap  & samples= m_sample_groups[group];
	auto & otumap= m_otus[group];
	typedef OrderMap  Torder;
	typedef std::map<StrTy, IdxTy > Tlut;
	Torder&  order =m_orders[""];
	Tlut m;
	Tlut & sample_filter=m_sample_filters[sfname];
	//Tokenizer cst;
	null_string_tokenizer  cst;
	// this needs to get one from the collection 
	if (use_map) 
	{ 
		const StrTy & map_name= m_flp.catagory_map();
		if (map_name==".") phyla_seq(m);
		else m=m_catagory_maps[map_name];
		MM_ERR(" catagoy mapping with "<<MMPR2(map_name,m.size()))
		order.append(m); 
 	} 
	MM_ERR(MMPR2(samples.size(),otumap.size()))
	// this iterates over all samples within group group 
	MM_LOOP(ii,samples)
	{
		const IdxTy wtf=(*ii).first;
		//MM_ERR(MMPR(wtf))
		const StrTy & sn=samples.st()(wtf);
	//	MM_ERR(MMPR(sn))
		if (sn.length()==0)
		{
			MM_ERR(" sample name is blank ignore ")
			continue; 
		}
		// apply filer map here
		if (filter_samples_from_map)
		{
			if ( sample_filter.find(sn)==sample_filter.end()) continue;
		}
	//	const  
		const Sample & s= (*ii).second;
		// if the name is not found, it outputs zero which the st() happily
		// translates into nonsensen
		const IdxTy sclass=s.get_class();
		const bool valid_class=(sclass!=(~0));
		//const StrTy&  sout=(output_class&&valid_class)?(samples.st()(sclass)):(sn);
		const StrTy&  sout=(output_class)?((valid_class?samples.st()(sclass):"noclass")):(sn);
		typedef std::map<StrTy, D> OtuCollated;
		//typedef std::map<IdxTy, D> OtuCollated;
		OtuCollated collate;
		D hit_total=0;
		MM_LOOP(jj,s)
		{

			const StrTy & otu=samples.st()((*jj).first);
			// this really needs to sum over the same flats 
			const Sample::otu_info & otuii=(*jj).second; // .flat();
			const StrTy otukey=otu; // samples.st()(otuii.otu());
			const StrTy otukey2= samples.st()(otuii.otu());
			if (otumap.find(otukey)==otumap.end())
			{ //// TODO FIXME wtf
				MM_ERR(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2))	
				MM_ONCE(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2),)	
			}
			const auto & otui=otumap[otukey]; // otumap[samples.st()(otuii.otu())]; // .flat();
			//const StrTy & flat=otui.flat();
			// GSMap & m= m_gs_maps[mapname];
//template <class Tm> StrTy  gs_classify(const Tm & m) const
			StrTy flat;
			// TODO FIXME reorg and include best_indexed_taxon() with a string
			// mashed name and taxon node id. 
//			if (use_map) flat=otui.classify_name(m);
//			else { if  (use_gs_map) { flat=otui.gs_classify(gs_map.m_map);}
//			else { flat=otui.string(otu_format); //MM_ERR("fack"<<MMPR2(flat,otu_format))
//}}
			switch (flatten_scheme)
			{
				case 0: { flat=otui.classify_name(m); break; } 
				case 1: { flat=otui.gs_classify(gs_map.m_map); break; }
		 		case 2: { flat=otui.string(otu_format);  break; } 
				case 3: { 	
							std::vector<StrTy> tv=otui.taxon();
							flat=m_tax_tree.best_indexed_taxon_str(tv,4);  
							break; 
						}
				default : {MM_ONCE(" bad default "<<MMPR(flatten_scheme),) }
			} // flatten_scheme


			// R can not handle ragged right FUYCK this 
			make_single(flat);	
			//MM_ERR(MMPR(flat))
			if (flat.length()==0) flat="blank"; // TODO needs to get all whitespace strings 
			if (debug) MM_ERR(MMPR3(samples.st()(otuii.otu()),otui.size(),flat))
			const D hits=otuii.hits();
			hit_total+=hits;
			if (hits!=0)
			{ 	
				// moved order.add here from separate loop not sure if good idea. 
				if (do_collate) {  collate[cst(flat)]+=hits;}
				if (update_order) {  order.add(flat);  }
				if (do_raw) {
				//os<<group<<sep<<sn<<sep<<otuii.otu();
				//os<<sep<<hits<<sep<<flat<<CRLF;}
				MM_MSG(group<<sep<<sn<<sep<<(samples.st()(otuii.otu()))<<sep<<hits<<sep<<flat) }
			}
		}
		// if anything was collated dump it now 
		D cum=0;
		if ( do_collate) { 
		// this was originally a simple itor over collate now much alower 
		// but allows use to calculate cumulative amounts and impose order 
		for(IdxTy i=0; i<order.size(); ++i)
	//	MM_LOOP(kko,ordv)
		{
			//const auto kk=collate.find(*kko);
			const StrTy flat=cst(order(i));
			//const auto kk=collate.find(cst(order(i)));
			const auto kk=collate.find(flat);
			// still want to output it, just not add zero to the maps 
			bool zed= (kk==collate.end());
//>>>			if (kk==collate.end()) continue;
// >>>			const auto & flat=cst((*kk).first); 
//			if ( order.find(flat)==order.end()) order[flat]=order.size();
			//const IdxTy serial=order[flat];
			const IdxTy serial=order(flat);
			auto  hits  =zed?0:(*kk).second; 
			if (hit_total!=0) hits=hits/hit_total;
			// this makes no sense now as the order has changed 
			cum+=hits;
			// os<<group<<sep<<sn; 
			// os<<sep<<hits<<sep<<flat<<CRLF;
			//MM_MSG(group<<sep<<sn<<sep<<hits<<sep<<cum<<sep<<flat<<sep<<serial) 
			MM_MSG(group<<sep<<sout<<sep<<hits<<sep<<cum<<sep<<flat<<sep<<serial) 
		} // kko  
} // colalte
//	samples[sn].add(names[i],::atof(li.word(i).c_str()));


	} // ii 
}
template <class Tm > void unknown_words(Tm & m) const
{
m[StrTy("unknown")]=1;
m[StrTy("unidentified")]=1;

}
void reconcileold(const StrTy &group,const IdxTy & settings,  std::ostream  & os=std::cout)
{
	SampleMap  & samples= m_sample_groups[group];
	auto & otumap= m_otus[group];
	MM_ERR(MMPR2(samples.size(),otumap.size()))
	null_string_tokenizer st;
	std::map<StrTy,IdxTy> exclude;
	unknown_words(exclude);
	// this iterates over all samples within group group 
	MM_LOOP(ii,samples)
	{
		const IdxTy wtf=(*ii).first;
		const StrTy & sn=samples.st()(wtf);
	//	MM_ERR(MMPR(sn))
		if (sn.length()==0)
		{
			MM_ERR(" sample name is blank ignore ")
			continue; 
		}
		const Sample & s= (*ii).second;
		MM_LOOP(jj,s)
		{

			const StrTy & otu=samples.st()((*jj).first);
			// this really needs to sum over the same flats 
			const Sample::otu_info & otuii=(*jj).second; // .flat();
			const StrTy otukey=otu; // samples.st()(otuii.otu());
			const StrTy otukey2= samples.st()(otuii.otu());
			if (otumap.find(otukey)==otumap.end())
			{ //// TODO FIXME wtf
				MM_ERR(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2))	
				MM_ONCE(" missing otu "<<MMPR3(otuii.otu(),otukey,otukey2),)	
			}
			// non const due to fixing 	
			auto & otui=otumap[otukey]; // otumap[samples.st()(otuii.otu())]; // .flat();
		
			std::vector<StrTy>  tv=otui.taxon();
			const bool remove_unknowns=true;
			if ( remove_unknowns)
			{
				IdxTy j=0;
				MM_SZ_LOOP(i,tv,sz)
				{
					if (exclude.find(tv[sz-1-i])==exclude.end()) break;
					++j;
				}
				otui.pop(j);
				if (j!=0)   tv=otui.taxon();
			}
				m_tax_tree.descend(tv);
		} // jj 

	} // ii 
} // reconcileold

void group_tax_tree(const StrTy &group,const IdxTy & settings,  std::ostream  & os=std::cout)
{
	auto & otumap= m_otus[group];
	MM_ERR("group_tax_tree trying to mak tree from "<<MMPR2(group,otumap.size()))
	if (otumap.size()==0)
	{
		MM_ERR(" group size is zero perhaps try these... ")
		MM_LOOP(ii,m_otus) { MM_ERR(MMPR2((*ii).first,(*ii).second.size())) } 
	}
	MM_ERR(" starting tax tree "<<MMPR(m_tax_tree.size()))
	null_string_tokenizer st;
	std::map<StrTy,IdxTy> exclude;
	//unknown_words(exclude);
	// this iterates over all samples within group group 
	MM_LOOP(ii,otumap)
	{
		// non const due to fixing 	
		const auto &  otukey=(*ii).first;
		auto & otui=otumap[otukey]; // otumap[samples.st()(otuii.otu())]; // .flat();
		std::vector<StrTy>  tv=otui.taxon();
		m_tax_tree.add_new_nodes(tv);
	} // ii 
	MM_ERR(" ending tax tree "<<MMPR(m_tax_tree.size()))

} // group_tax_tree
///////////////////////////////////////////////////////////////

void reconcile(const StrTy &group,const IdxTy & flags,  std::ostream  & os=std::cout)
{
	enum { REMOVE_UNKNOWN=1, MAKE_TAX_TREE=2, DUMP_TAXON=4};
	const bool remove_unknowns=((REMOVE_UNKNOWN&flags)!=0);
	const bool make_tax_tree=((MAKE_TAX_TREE&flags)!=0);
	const bool dump_taxon=((DUMP_TAXON&flags)!=0);
	auto & otumap= m_otus[group];
	MM_ERR(MMPR4(flags,remove_unknowns,make_tax_tree,dump_taxon)<<MMPR(otumap.size()))
	null_string_tokenizer st;
	std::map<StrTy,IdxTy> exclude;
	unknown_words(exclude);
	// this iterates over all samples within group group 
	MM_LOOP(ii,otumap)
	{
		// non const due to fixing 	
		const auto &  otukey=(*ii).first;
		auto & otui=otumap[otukey]; // otumap[samples.st()(otuii.otu())]; // .flat();
		std::vector<StrTy>  tv=otui.taxon();
		if (dump_taxon)
		{ 	Ss ss;
			ss<<otukey;
			MM_LOOP(jj,tv) { ss<<" "<<(*jj); }
			os<<MM_STR_LOC<<" "<<ss.str()<<CRLF;
		}
		if ( remove_unknowns)
		{
			IdxTy j=0;
			MM_SZ_LOOP(i,tv,sz)
			{
				if (exclude.find(tv[sz-1-i])==exclude.end()) break;
				++j;
			}
			otui.pop(j);
			if (j!=0)   tv=otui.taxon();
		}
		if (make_tax_tree) { m_tax_tree.descend(tv); } 

	} // ii 
} // reconcile





///////////////////////////////////////////////////////



void read_taxa_conform(const StrTy & fn,const StrTy & p2="") 
{
MM_ERR(" loading maps from adhoc names to confomed names ")
	std::ifstream fin(fn);
	CommandInterpretter li(&fin);
	IdxTy i=0; 
	auto & m =otu_struct::conform_map();
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		if ( sz<1) continue;
		StrTy keyin=li.word(0);
		if (keyin=="") continue;
		// the map only operates on conformed values so make the key conform 
		StrTy key=otu_struct::cuname(keyin,true,true,true,false);
		if ( key!=keyin ) { MM_ERR(" loading taxa key conformed "<<MMPR2(keyin, key)) } 	
		if (sz==1) m[key]="";
		else { ++i; m[key]=li.word(1); } 
	}
MM_ERR(" done loading maps from adhoc names to confomed names loaded "<<MMPR(i))
} // read_taxa_conform 

void dump_taxa_conform(const StrTy & fn,const StrTy & p2="") 
{
MM_ERR(" dumping  maps from adhoc names to confomed names ")
//	std::ifstream fin(fn);
	IdxTy i=0; 
	auto & m =otu_struct::conform_map();
	MM_LOOP(ii,m)
	{
		const StrTy & trivial=(*ii).first;
		const StrTy & conformed=(*ii).second;
		MM_MSG(MMPR3(i,trivial,conformed))
		++i;	
	}
MM_ERR(" done dumping  maps from adhoc names to confomed names loaded "<<MMPR(i))
} // dump_taxa_conform 

void catalog_taxa_names()
{
typedef std::map<StrTy,IdxTy>  CntMap;
typedef std::map<IdxTy, CntMap> LvlMap;
typedef std::map<StrTy,LvlMap> GrMap; 
GrMap ctmap,ctmapc;
LvlMap globmap,globmapc;
IdxTy maxlvels=0;
MM_LOOP(ii,m_sample_groups)
{
	const StrTy & group=(*ii).first;
	LvlMap & gmap=ctmap[group];
	LvlMap & gmapc=ctmapc[group];
	auto & otumap= m_otus[group];
	SampleMap  & samples= (*ii).second; // m_sample_groups[group];
	MM_LOOP(jj,samples)
	{
		string_tokenizer & st=samples.st();
		const StrTy & sn=st((*jj).first);
		const Sample & s= (*jj).second;
		MM_LOOP(kk,s)
		{
			const StrTy & otu=st((*kk).first);
			const Sample::otu_info & otuii=(*kk).second; // .flat();
			const otu_struct & otui=otumap[st(otuii.otu())]; // .flat();
			MM_SZ_LOOP(lvel,otui,sz)
			{	
				const StrTy nm= otui.level(lvel,false);
				const StrTy nmc= otui.level(lvel,!false);
				++gmap[lvel][(nm)];	
				++gmapc[lvel][(nmc)];	
				++globmap[lvel][(nm)];
				++globmapc[lvel][(nmc)];	
			} // lvel 
			if (sz>=maxlvels) maxlvels=sz;

		}  // kk 
	} // jj 
} // ii 
const StrTy sep=" ";
const StrTy lbl=MM_STR_LOC;
std::ostream & os=std::cout;
for (IdxTy i=0; i<maxlvels; ++i)
{

MM_LOOP(ii,globmap[i])
{
 //	MM_LOOP(jj,(*ii).second)
	{
		Ss ss;
		const StrTy nm=(*ii).first;
		const IdxTy gcnt=(*ii).second;
		ss<<lbl<<sep<<MMPR3(i,nm,gcnt);
		MM_LOOP(gg,m_sample_groups)
		{
			const StrTy & group=(*gg).first;
			// non const so [] works ... 
			LvlMap & gmap=ctmap[group];
//			LvlMap & gmapc=ctmapc[group];
			ss<<MMPR(group)<<MMPR(gmap[i][nm]);
		} // gg 
		os<<ss.str()<<CRLF;	
	} // jj

} // ii 


} // i 


} // catalog_taxa_names


// should be obsolete
///////////////////////////////////////////////////////
template <class Tv, class Tm> 
void invert_map(Tv & v, const Tm & m)
{
v = Tv(m.size());
MM_LOOP(ii,m)
{
v[(*ii).second]=(*ii).first; 
}
}

template <class Tvv> 
void all_samples( Tvv & v)
{
	MM_LOOP(ii,m_sample_groups)
	{
		std::vector<StrTy>  vv;
		vv.push_back((*ii).first);
		SampleMap & samples=(*ii).second;
		//MM_LOOP(jj, (*ii).second)
		MM_LOOP(jj, samples)
		{
			vv.push_back(samples.st()((*jj).first));
		}
		v.push_back(vv);	
	} // ii 
} // all_samples 
template <class Tvv> 
StrTy gs_header( const Tvv & v)
{
Ss ss;
StrTy gr="";
const StrTy sep=" ";
const StrTy concat="_";
ss<<"serial"<<sep<<"catagory";
MM_LOOP(ii,v)
{
MM_LOOP(jj,(*ii))
{
if (jj==(*ii).begin()) { gr=*jj; continue; } 
ss<<sep<< gr<<concat<<(*jj);
//sep=" ";
}

}
return ss.str();
}
void dump_collated_group(std::ostream & os,  const StrTy  &group)
{
	const StrTy sep=" ";
	typedef std::map<StrTy, IdxTy> Tcat;
	typedef std::map<IdxTy, std::map< StrTy, std::map< StrTy, D > > >  Tcompos;
	typedef std::vector< std::vector< StrTy > > GSVec;
	GSVec allgs;
	all_samples(allgs);	
	Ss ss;
	ss<<gs_header(allgs);
	ss<<CRLF;	
	Tcompos d;
	Tcat m;
 	phyla_seq( m);
	std::vector<StrTy> groups,imap;
	invert_map(imap,m);
	groups.push_back(group);
	//collate_groups(d,m,groups);
	collate_groups(d,m,allgs);
	MM_LOOP(ii,d)
	{
		const IdxTy&  otu_class= (*ii).first;
		ss<<MM_MARK<<otu_class<<sep<<imap[otu_class];
		MM_LOOP(jj,allgs)
		{
			const auto&  gr=(*jj);
			StrTy grs="xxx";
			MM_LOOP(kk,gr)
			{
				const StrTy  & sample=(*kk);
				// this is kind of a kluge doh 
				if ( kk==gr.begin()) { grs=*kk; continue; } 
				const D & hits=(*ii).second[grs][sample];
				//ss<<sep<<grs<<sample<<sep<<hits;	
				ss<<sep<<hits;	

			} // kk or sample 

		} // jj or group
		ss<<CRLF;
 	} // ii or class 
os<<ss.str();
}

// don't remembe what this does 
template <class Tmap, class Tcat, class Tg>
void collate_groups(Tmap & d, const Tcat & m, const Tg &groups)
{
	const bool debug=false; 
	const StrTy sep=" ";
	const IdxTy otu_format=m_flp.otu_format();
// each vector in groups has a group as first element followed
// by samples 
	MM_LOOP(kk,groups)
{
	//const StrTy & group=(*kk);
	const StrTy & group=(*kk)[0];
	SampleMap & samples= m_sample_groups[group];
	auto & otumap= m_otus[group];
	for (auto ll=(++(*kk).begin()); ll!=(*kk).end(); ++ ll)
	{
//	MM_ERR(MMPR2(samples.size(),otumap.size()))
	// this iterates over all samples within group group 
	//MM_LOOP(ii,samples)
	{
		const StrTy & sn=(*ll); // (*ii).first;
		// TODO the sting tokenizer is not const correct wtf
		const  
			Sample & s= samples[samples.st()(sn)] ; // (*ii).second;
		typedef std::map<StrTy, D> OtuCollated;
		OtuCollated collate;
		D hit_total=0;
		MM_LOOP(jj,s)
		{

			const StrTy & otu=samples.st()((*jj).first);
			// this really needs to sum over the same flats 
			const auto & otuii=(*jj).second; // .flat();
			const auto & otui=otumap[samples.st()(otuii.otu())]; // .flat();
//			const StrTy & flat=otui.string(otu_format);
			const IdxTy otu_class=otui.classify(m);
//			if (debug) MM_ERR(MMPR3(otuii.otu(),otui.size(),flat))
			const D hits=otuii.hits();
			hit_total+=hits;
			if (hits!=0)
			{ 	
				d[otu_class][group][sn]+=hits;

			}

		}
		// if anything was collated dump it now 
		for (IdxTy otu_class=0; otu_class<m.size(); ++otu_class)
		{
			if (hit_total!=0) //hits=hits/hit_total;
				d[otu_class][group][sn]/=hit_total;
		} // collate

} 
	} // ii 
} // kk 
} // collate_groups 

 void add_catagory(const StrTy & name, const StrTy & cat)
{
m_catagory_maps[name][cat]=(m_catagory_maps[name]).size();

}
 void dump_catagory_map(const StrTy & name)
{
std::ostream & os=std::cout;
auto & m=m_catagory_maps[name];
MM_LOOP(ii,m)
{
MM_ERR(" cat map "<<MMPR3(name,(*ii).first,(*ii).second))
}
}
 void read_sample_filter(const StrTy & name, const StrTy & fn)
{

read_catagory_map(m_sample_filters[name],fn);
MM_ERR(" sample filter "<<MMPR2(name,m_sample_filters[name].size()))

}

 void load_catagory_map(const StrTy & name, const StrTy & fn)
{
read_catagory_map(m_catagory_maps[name],fn);

}

template <class Tmap> void read_catagory_map(Tmap & m, const StrTy & fn)
{

std::ifstream fin(fn);
CommandInterpretter li(&fin);
IdxTy i=0; 
while (li.nextok())
{
if (li.size()<1) continue;
const StrTy & w=li.word(0);
if (w.length()==0) continue;
const StrTy wcon= otu_struct::cuname(w,true,true,true,true);
if (w!=wcon) { MM_ERR(" coercing read_catagory_map key "<<MMPR2(w,wcon)) } 
m[wcon]=i;
++i;
} // li 

}
 void read_canon_map( const StrTy & fn)
{
MM_ERR(" reading canon map file "<<fn)
std::ifstream fin(fn);
CommandInterpretter li(&fin);
const IdxTy base=0;
IdxTy i=0; 
while (li.nextok())
{
if (li.size()<(3+base)) continue;
const StrTy & w=li.word(base); // level 
const StrTy & y=li.word(base+1); // canonized 
const StrTy & z=li.word(base+2); // zero for now just a flag 
if (w.length()==0) continue;
const IdxTy level=myatoi(w);
const IdxTy flag=myatoi(z);
m_clm[level][y]=flag;
++i;
} // li 
MM_ERR(" done with m_clm "<<MMPR2(m_clm.size(),i))
}







/////////////////////////////////////////////////////
// for cats that are not at same level, archaea for example,
// the otu needs to figure out which cat it is in based
// on a hit somehwere in the hieratchy 
template <class Tmap> void phyla_seq(Tmap & m)
{
IdxTy seq=0;
m["blank"]=seq; ++seq; 
//m["arcobacter"]=seq; ++seq; 
//m["helicobacteraceae"]=seq; ++seq; 
m["acidobacteria"]=seq; ++seq; 
m["actinobacteria"]=seq; ++seq; 
m["bacteroidetes"]=seq; ++seq; 
m["candidate"]=seq; ++seq; 
m["candidatedivisionjs1"]=seq; ++seq; 
m["candidatedivisionop3"]=seq; ++seq; 
m["candidatedivisionop8"]=seq; ++seq; 
m["candidatedivisionop9"]=seq; ++seq; 
m["chlorobi"]=seq; ++seq; 
m["chloroflexi"]=seq; ++seq; 
m["cyanobacteria"]=seq; ++seq; 
m["deferribacteres"]=seq; ++seq; 
m["euryarchaeota"]=seq; ++seq; 
m["fibrobacteres"]=seq; ++seq; 
m["firmicutes"]=seq; ++seq; 
m["fusobacteria"]=seq; ++seq; 
m["gemmatimonadetes"]=seq; ++seq; 
m["hydrogenedentes"]=seq; ++seq; 
m["lentisphaerae"]=seq; ++seq; 
m["planctomycetes"]=seq; ++seq; 
m["proteobacteria"]=seq; ++seq; 
m["saccharibacteria"]=seq; ++seq; 
m["spirochaetae"]=seq; ++seq; 
m["synergistetes"]=seq; ++seq; 
m["tenericutes"]=seq; ++seq; 
m["thermodesulfobacteria"]=seq; ++seq; 
m["thermotogae"]=seq; ++seq; 
m["unidentified"]=seq; ++seq; 
m["verrucomicrobia"]=seq; ++seq; 
}


void  add_gs_map(const StrTy & mapname, const StrTy & fn) 
{
	MM_ERR(" add_gs_map enter only blank supported now  "<<MMPR2(mapname,fn) ) 
	typedef string_tokenizer St;
	std::ifstream is(fn);
	GSMap & m= m_gs_maps[mapname];
	St st;
	CommandInterpretter li(&is);
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		// this requires 3 words : at least a genus , a "|" , and a value 
		IdxTy j=0;	
		StrTy genus="";
		StrTy species="";
		StrTy value="";
		for (IdxTy i=0; i<sz; ++i)
		{
			const StrTy & w=li.word(i);
			if (w.length()<1) continue;
			if (w=="|") { j=3; continue; } 
			if (j==0)  genus=w; 
			else if (j==1) species=w;
			else if ( j==3) { value=w;	 break; } 
			++j;
		} // i 
		if (value.length()==0) value=StrTy("blank_")+genus+StrTy("_")+species;
		m.m_map[genus][species]=value;

	} // li 
	MM_ERR(" add_gs_map exit size now "<<MMPR(m.m_map.size()) ) 
} // add_gs_map


/*
Read a file generated like this from z5.txt,

classify-sample-group foo /home/marchywka/junk/sample_envo
update-order-format foo
sort-order
dump-class-format foo



mjm_zymo.h1721 foo ENVO:forest 6.11303e-06 0.000409573 archaea_euryarchaeota_thermoplasmata_e2_f___g___s___ 135
mjm_zymo.h1721 foo ENVO:forest 0.0269768 0.0273864 bacteria_ 140
mjm_zymo.h1721 foo ENVO:forest 0.00440138 0.0317878 bacteria_acidobacteria_ 169
mjm_zymo.h1721 foo ENVO:forest 0.000843598 0.0326314 bacteria_acidobacteria_[chloracidobacteria]_11-24_f___g___s___ 171
mjm_zymo.h1721 foo ENVO:forest 8.55824e-05 0.0327169 bacteria_acidobacteria_[chloracidobacteria]_ellin7246_f___g___s___ 175
mjm_zymo.h1721 foo ENVO:forest 3.05651e-05 0.0327475 bacteria_acidobacteria_[chloracidobacteria]_o___f___g___s___ 176
mjm_zymo.h1721 foo ENVO:forest 9.78085e-05 0.0328453 bacteria_acidobacteria_[chloracidobacteria]_pk29_f___g___s___ 177

and output stats of class and otu <%>, <cum>,min, maxetc

echo dump-otf-class-stats zzz7 | ./mjm_zymo.out 2>xxx | mjm eq > aaa7



*/
class otf_stats { 

public :
otf_stats():n(0),pct(0),cum(0),minpct(1),maxpct(0),mincum(1),maxcum(0),pct2(0),cum2(0) {}
//	sm[st(clas)][st(otu)].add(atof(pct.c_str()), atof(cum.c_str()));	
void add( const D & _pct, const D & _cum)
{
++n; pct+=_pct; cum+=_cum;
maxx(maxpct,_pct);
minn(minpct,_pct);
maxx(maxcum,_cum);
minn(mincum,_cum);
pct2+=_pct*_pct;
cum2+=_cum*_cum;
}
void maxx(D & x, const D & y) { if ( y>x ) x=y; } 
void minn(D & x, const D & y) { if ( y<x ) x=y; } 

void dump(std::ostream & os, const IdxTy & flags=0)
{
switch ( flags)
{
case 0 :{ os<<MMPR(n)<<MMPR4(pct,cum,minpct,maxpct)<<MMPR4(mincum,maxcum,pct2,cum2); break; }
case 1 :{ 
const D apct=pct/n; 
const D acum=cum/n;
const D sdpct=sqrt((pct2/n-apct*apct));
const D sdcum=sqrt((cum2/n-acum*acum));
os<<MMPR(n)<<MMPR4(apct,acum,minpct,maxpct)<<MMPR4(mincum,maxcum,sdpct,sdcum); break; }
case 2 :{ 
const D apct=pct; 
const D acum=cum;
const D sdpct=sqrt((pct2-apct*apct));
const D sdcum=sqrt((cum2-acum*acum));
os<<MMPR(n)<<MMPR4(apct,acum,minpct,maxpct)<<MMPR4(mincum,maxcum,sdpct,sdcum); break; }

// TODO  dmel... 
default: { MM_ERR(" bad basae "<<flags) }
}  //switch 
} // dump


IdxTy n;
D pct,cum,minpct,maxpct,mincum,maxcum,pct2,cum2;

 };

void dump_otf_class_stats(const StrTy & fn)
{
std::ifstream is(fn);
 otf_class_stats(std::cout,  is) ;
}
/*
z8.txt

read-catagory-map phyla phyla.txt
dump-catagory-map phyla
set-param catagory_map phyla
set-param use_sample_filter 0
#read-sample-filter aa /home/marchywka/junk/aa 
parse-biom-json emp /home/marchywka/junk/emp.json
#dot-sparse-biom zymo foo 
copy-biom emp foo
clear-biom emp

parse-biom-json zymo /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom
copy-biom zymo zymo
clear-biom zymo
#dump-otu
#catalog-otu
set-param otu_format 0x000000
classify-sample-group foo /home/marchywka/junk/sample_envo
# RAW=1, COLLATE=2, UPDATE_ORDER=4, USE_MAP=8,OUTPUT_CLASS=16, GS_MAP=32, DEBUG_INFO=(1<<16)
update-order-format foo 0x0C
update-order-format zymo 0x0C
sort-order
dump-group-n foo  0x01E      
dump-group-n zymo 0x00E      


to reorganize output for taxaplot2.R use this, 

4490  echo dump-otf-class-stats zzz8 | ./mjm_zymo.out  2>xxx | mjm eq | awk '{print "dummy group "$1" "$6" "$8" "$2" "$NF}' > ppp


*/

//mjm_zymo.h1721 foo ENVO:forest 9.78085e-05 0.0328453 bacteria_acidobacteria_[chloracidobacteria]_pk29_f___g___s___ 177
// the outpu stream  after mjm eq to remove = 

//ENVO:podzol nomatch  n 4 apct 0.980886 acum 0.980886 minpct 0.954457 maxpct 0.997024 mincum 0.954457 maxcum 0.997024 sdpct 0.016795 sdcum 0.016795
//ENVO:podzol plantrelated  n 4 apct 0.00293145 acum 0.983818 minpct 0.000537316 maxpct 0.00473277 mincum 0.959081 maxcum 0.997561 sdpct 0.00180642 sdcum 0.0153096

void  otf_class_stats(std::ostream & os, std::istream & is) 
{
	MM_ERR(" staart otf_class_stats") 
	typedef string_tokenizer St;
	typedef otf_stats Ac;
	typedef std::map<IdxTy, Ac> Amap;
	typedef std::map<IdxTy, Amap> Statmap;
	typedef std::map<IdxTy, D> Pct;
	typedef std::map<IdxTy,  std::vector<IdxTy> > SeqMap;
//	typedef std::map<IdxTy, IdxTy> Ordermap;

	Statmap sm;
	St st;
	Pct pctm,cumm;
SeqMap seq;
//	Ordermap omap,imap;
	OrderMap iomap; //collate_order_map
	CommandInterpretter li(&is);
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		if (sz!=7)
		{
			// use Dmel
			MM_ERR(" bad input line "<<li.line())
			continue;
		}
		const Words & w=li.words();
		const StrTy & clas =  w[2];
		const StrTy & pct=w[3];
		const StrTy & cum=w[4];
		const StrTy & otu=w[5];
		const StrTy & order=w[6];
		const IdxTy stc=st(clas);
		const D dpct=atof(pct.c_str());
//		omap[st(otu)]=myatoi(order);
//		imap[myatoi(order)]=st(otu);
		const IdxTy oi=myatoi(order);
		iomap.add_at(otu,oi);
		seq[stc].push_back(oi);
//		omap.add(otu);
// TODO FIXME these cumulative things make no sense at this point 
	//sm[st(clas)][st(otu)].add(dpct, atof(cum.c_str()));	
	sm[stc][st(otu)].add(dpct, 0);	
	pctm[stc]+=dpct; // atof(pct.c_str());
	} // li
	MM_LOOP(ii,seq) { std::sort((*ii).second.begin(), (*ii).second.end()); } 
	const StrTy sep=" ";
	MM_LOOP(ii,sm)
	{
		auto & l1=(*ii).second;
		const IdxTy classno=(*ii).first; 
		const StrTy n1= st(classno);
		// this is a loop over the map and the cum order is now wrong 
		//MM_LOOP(jj,l1)
//		D cumnew=0;
		// too slow and does not work wtf 
		//MM_SZ_LOOP(j,iomap,sz)
		IdxTy jlast=~0;
		MM_LOOP(kk,seq[classno])
		{
			const IdxTy j=(*kk);
			if ( j==jlast) continue;
			jlast=j;
			const StrTy n2=iomap[j]; //  st((*jj).first);
			const IdxTy order=j ; // omap[(*jj).first];
			auto  jj=l1.find(st(n2));
			if (jj!=l1.end())
			{
				auto & stat= (*jj).second;
				const D op=stat.pct;
				const D newpct=op/pctm[classno];
				stat.pct=newpct;
				cumm[classno]+=stat.pct; // /stat.n;
//				cumnew+=newpct;	
				stat.cum=cumm[classno];
				// the stat members are public, just fix the pct and cum 
				os<<n1<<sep<<n2<<sep;
				stat.dump(os,1+1);
				os<<sep<<order;
				os<<CRLF;
			}
			else { MM_ERR(" missing entry "<<MMPR3(n1,n2,j)) } 
		} //  kk

	} // ii
	MM_ERR(" end otf_class_stats") 
} // otf_class_stats

// typical input to taxaplot2.R "ppp 
// dummy group ENVO:forest 0.378233 0.378233 acidobacteria 1
// this is backwards to allow output defalt 
void  otf_angles(const StrTy & fin, const StrTy & fout)
{
if (fin!="") {std::ifstream is(fin); otf_angles(std::cout,is); return;  }
otf_angles();
}
 
void  otf_angles(std::ostream & os=std::cout, std::istream & is=std::cin) 
{
	MM_ERR(" start  otf_angles")

	typedef string_tokenizer St;
	St st;
//	MySparse m;
	typedef std::map<IdxTy, std::map< IdxTy, D >> Sparse;
	typedef std::map<IdxTy, D> MeanMap;
	MeanMap mm;
	Sparse m;
	CommandInterpretter li(&is);
	while (li.nextok())
	{
		const IdxTy sz=li.size();
		if (sz!=7)
		{
			// use Dmel
			MM_ERR(" bad input line "<<li.line())
			continue;
		}
		const Words & w=li.words();
		const StrTy & clas =  w[2];
		const StrTy & pct=w[3];
		const StrTy & cum=w[4];
		const StrTy & otu=w[5];
		const StrTy & order=w[6];
		const IdxTy coord=myatoi(order);
		const IdxTy vec=st(clas);
		const D val=atof(pct.c_str());
		mm[coord]+=val;
		m[vec][coord]+=val ; // atof(pct.c_str());
	//sm[st(clas)][st(otu)].add(atof(pct.c_str()), atof(cum.c_str()));	
	} // li
	D sum=0;
	MM_LOOP(ii,mm) { sum+=(*ii).second * (*ii).second; } 
const D neff=mm.size();
	D origin=-sum/neff/neff+2.0/neff; 
	const IdxTy ms=m.size();
	const StrTy sep= " ";
	MyBlock p(ms,ms);
	MM_LOOP(ii,m)	
	{
		const IdxTy i=(*ii).first;
		MM_LOOP(jj,m)
		{
			const IdxTy j=(*jj).first;
			MM_LOOP(kk,(*jj).second) 
			{ 	const IdxTy coord=(*kk).first;  
				if ( m[i].find(coord)!=m[i].end()) p(i,j)+=m[i][coord]*m[j][coord]; 
			} 
		}

	} // ii 
	MM_ERR( " subtracting compuated origin of "<<MMPR(origin))
//	for(IdxTy i=0; i<ms; ++i) for(IdxTy j=0; j<ms; ++j) p(i,j)-=origin;
	for(IdxTy i=0; i<ms; ++i)
	{
		for(IdxTy j=0; j<ms; ++j)
		{
			const D & v=p(i,j);
			const D & d=sqrt(p(i,i)*p(j,j));
			if ( v!=0) os<<"otf_angles "<<st(i)<<sep<<st(j)<<sep<<(v/d)<<CRLF;
			else {MM_ONCE(" orthogonal groups not printed "<<MMPR2(st(i),st(j)),) } 
		} // j 
	} // i 

	MM_ERR(" done otf_angles")

} // otf_angles



// stub in progress TODO FIXME 
void  otu_stats(const StrTy & group, const IdxTy & settings) 
{
	//typedef string_tokenizer St;
	//typedef otf_stats Ac;
	//typedef std::map<IdxTy, Ac> Amap;
//	typedef std::map<IdxTy, Amap> Statmap;
//	typedef std::vector<IdxTy> Taxv;

	//typedef std::map<Taxv,Ac> StatMap;

	const bool filter_samples_from_map=false; 

	SampleMap  & samples= m_sample_groups[group];
	MM_LOOP(ii,samples)
	{
		const StrTy & sn=samples.st()((*ii).first);
		if (sn.length()==0)
		{
			MM_ERR(" sample name is blank ignore ")
			continue; 
		}
		// apply filer map here
		if (filter_samples_from_map)
		{
//			if ( sample_filter.find(sn)==sample_filter.end()) continue;
		}
	// facking tokenizer not const correct 
	//	const  
		const Sample & s= (*ii).second;
		//const StrTy&  sout=output_class?(samples.st()(s.get_class())):(sn);
		//typedef std::map<IdxTy, D> OtuCollated;
		//OtuCollated collate;
		D hit_total=0;
		MM_LOOP(jj,s)
		{
			// still needs to normalize this crap... uggghh 
			const StrTy & otu=samples.st()((*jj).first);

		} // jj 
	} // ii 


} // otu_stats


/*

void end_of_day_check()
{
// make sure macro usages make sense 

}
// TODO FIXME this really needs an iteroatr 
// better off just to make itor 

*/


/////////////////////////////////////////////////////////////////////
private:
void Init()
{
//delete m_dmel;
//m_dmel= new Dmel();
m_done=false;
//load_literals();
//load_shares();
//load_handlers();
}


// MEMBERS

bool m_done;

Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;

OtuMap m_otus;
SampleGroup m_sample_groups;
BiomFiles m_biom_files;
CatagoryMaps m_catagory_maps;
CollateOrders m_orders;
SampleFilters m_sample_filters;
GSMaps m_gs_maps;
otu_conform::CanonLevelMap m_clm;

synonym_canon_lut  m_canon_lut;

TaxTree m_tax_tree;
RaggedMap m_ragged_map;

mjm_tree_viz m_tvz;
PhenoNotes m_pheno_notes;

CounterMap m_cm;
//typedef std::map<StrTy,CompTable> CompMap;
CompMap m_comps_map;
FastaMap m_fasta_map;
Recon m_recon;
}; //mjm_zymo

#ifdef PYTHON_BOOST_BUILD
#warning building boost python
using namespace boost::python;
//namespace {
typedef mjm_zymo Tg;
static Tg  * single() {
static Tg * p = new Tg();
return p;
}
class python_foo {
public:
static void cmd(const StrTy & s) { (*single()).command_mode(s); } 
}; // python_foo

//}; // namespace 

// wtf lol, 
// strings mjm_zymo.so | grep PyInit
//PyInit_Tg
//BOOST_PYTHON_MODULE(Tg)
BOOST_PYTHON_MODULE(mjm_zymo)
{
        //.def("command_mode", &fack_python::FICKBOOSTs)
      //  boost::python::def("command_mode", &fack_python::FICKBOOSTs) ;
        boost::python::def("cmd", &python_foo::cmd) ;
}


#endif



/////////////////////////////////////////////////////////

#ifdef  TEST_ZYMO__
int main(int argc,char **args)
{
typedef mjm_zymo  Myt;
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
const unsigned int n=Myt::myatoi(args[2]); // 100;
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

