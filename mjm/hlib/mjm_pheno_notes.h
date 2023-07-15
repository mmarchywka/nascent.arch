#ifndef MJM_PHENO_NOTES_H__
#define MJM_PHENO_NOTES_H__
 
#include "mjm_globals.h"
// for the presence absence vector 
#include "mjm_ordering.h"
#include "mjm_char_mat.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
//#include "../mjm_fasta_ii.h"
#include "mjm_fasta_ii.h"


#include "mjm_collections.h"
#include "mjm_svg_writer.h"


//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
//#include "mjm_biom_hdf5.h"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>
#include <stdint.h>

/*
TODO FIXME
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#ifndef MM_DMEL
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 
#endif
/*

 ./mjm_pheno_notes.out -cmd "tt load-tax zymo zymo-tt-info" -cmd "load-hbpheno foo gs_16s_attributes.txt" -cmd "hbpheno-annotate zymo foo  0x0100"
 3184  ./mjm_pheno_notes.out -source z1pheno.txt -quit 2>xxx 
 3185  ./mjm_pheno_notes.out -cmd "tt load-tax zymo zymo-tt-info" -cmd "load-hbpheno foo gs_16s_attributes.txt" -cmd "hbpheno-annotate zymo foo  0x0100" 2>xxx
 3186  ./mjm_pheno_notes.out -cmd "tt load-tax zymo zymo-tt-info" -cmd "load-hbpheno foo gs_16s_attributes.txt" -cmd "hbpheno-annotate zymo food  0x0100"
 3187  ./mjm_pheno_notes.out -cmd "tt load-tax zymo zymo-tt-info" -cmd "load-hbpheno foo gs_16s_attributes.txt" -cmd "hbpheno-annotate zymo food  0x0100" -quit 2>xxx
 3188  ./mjm_pheno_notes.out -cmd "tt load-tax zymo zymo-tt-info" -cmd "load-hbpheno foo gs_16s_attributes.txt" -cmd "hbpheno-annotate zymo food  0x0100" -quit 2>xxx > gmp16.txt
 3793  vi mjm_pheno_notes.h


cat z1pheno.txt 
load-protrait foo /home/marchywka/junk/pheno/IJSEM_pheno_db_v1.0.txt
#query-protrait foo "helicobacter" "spore production"
tt load-tax zymo zymo-tt-info
#annotate zymo foo "spore production"
#annotate zymo foo "oxygen preference"
#annotate zymo foo "Sole carbon substrate use"
#load-pheno foo /home/marchywka/junk/pheno/ProTraits_precisionScores.txt  
load-hbpheno foo gs_16s_attributes.txt
#pheno-annotate zymo foo pathogen asdf qer
#hbpheno-annotate zymo foo pathogen asdf qer
hbpheno-annotate zymo foo  0x0100 



*/


////////////////////////////////////////////////////////////////

class pheno_notes_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
pheno_notes_params( const StrTy & nm) : Super(nm) {}
pheno_notes_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//IdxTy maxdepth() const { return m_map.get_uint("maxdepth",3); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }
//bool print_haves() const { return m_map.get_bool("print_haves",!true); }
//bool print_havenots() const { return m_map.get_bool("print_havenots",!true); }
//bool print_if_have() const { return m_map.get_bool("print_if_have",!true); }
//bool suppress_vector() const { return m_map.get_bool("suppress_vector",!true); }
//bool add_level() const { return m_map.get_bool("add_level",true); }
//bool print_hit() const { return m_map.get_bool("print_hit",true); }
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace pheno_notes_traits
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
//typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
//typedef unsigned int KeyCode;
//typedef unsigned __int128 KeyCode;
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
typedef std::vector<IdxTy > Locations;
typedef mjm_string_ordering Ordering;

//typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;

class util
{
public:
static void Split( Words & w, const StrTy & s) // const
{ //Traits::Split(w,s); 
const IdxTy sz=s.length();
const char * sc=s.c_str();
char c[sz+1];
//::memcpy(c,sc,sz);
IdxTy dptr=0;
IdxTy sptr=0;
while ( sc[sptr] != 0 )
{
const char cp=sc[sptr];  ++sptr;
const bool ender= ( cp==',')||(cp=='(')||(cp==')')||(cp=='"')||(cp=='|');
if (ender) 
{
if (dptr!=0){ c[dptr]=0; w.push_back(StrTy(c)); dptr=0;}}
else { c[dptr]=cp; ++dptr; } 
 
} // while sc sptr
if (dptr!=0) { c[dptr]=0; w.push_back(StrTy(c)); }

} //split



static void Split2( Words & w, const StrTy & s) // const
{ //Traits::Split(w,s); 
const IdxTy sz=s.length();
const char * sc=s.c_str();
char c[sz+1];
//::memcpy(c,sc,sz);
IdxTy dptr=0;
IdxTy sptr=0;
while ( sc[sptr] != 0 )
{
const char cp=sc[sptr];  ++sptr;
const bool ender= (cp==' ')||( cp==',')||(cp=='(')||(cp==')')||(cp=='"')||(cp=='|');
if (ender) 
{
if (dptr!=0){ c[dptr]=0; w.push_back(StrTy(c)); dptr=0;}}
else { c[dptr]=cp; ++dptr; } 
 
} // while sc sptr
if (dptr!=0) { c[dptr]=0; w.push_back(StrTy(c)); }

} //split

 template <class Ty>
static void reduce(std::ostream &ss, const Ty & w , const IdxTy flags)
{
//MM_ERR(" old reduce ")
//if ( w.size()>0) {reduce(ss,w,0); } 
switch (flags)
{
case 1:
{
std::map<StrTy,IdxTy> m;
std::map<D,StrTy> ms;
MM_LOOP(ii,w) { ++m[(*ii)]; }
MM_ONCE(" danger will robinson fix this code or get unreliable oracle ",)
MM_LOOP(ii,m) { ms[-((*ii).second + 1.0/(ms.size()+10))]=(*ii).first;  }
//MM_LOOP(ii,m) { ss<<" "<<(*ii).second<<":"<<(*ii).first; } 
MM_LOOP(ii,ms) { ss<<" "<<(*ii).second<<":"<<IdxTy(-(*ii).first); }
break;
}
default:
{ MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  }
}
}

 template <class Tm, class Ty>
static void reduce_unique(Tm & wr, const Ty & w , const IdxTy flags)
{

std::map<StrTy,IdxTy> m;
//std::map<D,StrTy> ms;
MM_LOOP(ii,w) { ++m[(*ii)]; }
MM_ONCE(" danger will robinson fix this code or get unreliable oracle ",)
//MM_LOOP(ii,m) { ms[-((*ii).second + 1.0/(ms.size()+10))]=(*ii).first;  }
MM_LOOP(ii,m) { wr.push_back((*ii).first); } 
}

static bool fast_grep(const StrTy &s, const StrTy & p)
{
const char * pp=p.c_str();
const char * sp=s.c_str();
const IdxTy psz=p.length();
const IdxTy ssz=s.length();
if (ssz<psz) return false; 
const IdxTy sd=ssz-psz; 
for(IdxTy i=0; i<sd; ++i)
{
if (strncmp(sp+i,pp,psz)==0) return true; 

}
return false; 
}


 template <class Ty>
static void reduce(std::ostream &ss, const Ty & _w , const IdxTy flags,const StrTy & field)
{
MM_ERR(" field  reduce ")
Ty w;
MM_LOOP(ii,_w) { if ( fast_grep((*ii),field)) w.push_back((*ii)); }
//if ( w.size()>0) {reduce(ss,w,0); } 
switch (flags)
{
case 1:
{
std::map<StrTy,IdxTy> m;
std::map<D,StrTy> ms;
MM_LOOP(ii,w) { ++m[(*ii)]; }
MM_LOOP(ii,m) { ms[-((*ii).second + 1.0/(ms.size()+10))]=(*ii).first;  }
//MM_LOOP(ii,m) { ss<<" "<<(*ii).second<<":"<<(*ii).first; } 
MM_LOOP(ii,ms) { ss<<" "<<(*ii).second<<":"<<IdxTy(-(*ii).first); }
break;
}
default:
{ MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  }
}
}






}; // util

}; // trees_and_tables_traits
///////////////////////////////////////////////////////////////

class mjm_home_brew_pheno
{
/*
cat ProTraits_precisionScores.txt  | more
Tax_ID;Organism_name;Phenotype;Minority;Original_label;Hamap_+;Hamap_-;MicrobeWi
ki_+;MicrobeWiki_-;Wikipedia_+;Wikipedia_-;Other_text_+;Other_text_-;PMC_article
s_+;PMC_articles_-;Pubmed_abstracts_+;Pubmed_abstracts_-;Proteome_composition_+;
Proteome_composition_-;Phyletic_profiles_+;Phyletic_profiles_-;Coexisting_microb
es_+;Coexisting_microbes_-;Gene_neighbourhood_+;Gene_neighbourhood_-;Codon_usage
_bias_+;Codon_usage_bias_-;Integrated_score_+;Integrated_score_-
7;Azorhizobium caulinodans;alkaline_phosphatase;-;?;0.58;0.033;?;?;?;?;?;?;?;?;0
.76;0.198;1.0;0.128;0.89;0.219;0.648;0.159;0.859;0.161;0.866;0.149;0.89;0.198
7;Azorhizobium caulinodans;arbutin;+;?;0.287;0.33;?;?;?;?;?;?;?;?;?;?;0.274;0.23
3;0.608;0.639;?;?;0.078;0.307;0.225;0.508;0.287;0.508
7;Azorhizobium caulinodans;cellobiose;-;?;0.008;1.0;?;?;?;?;?;?;?;?;?;?;0.156;0.
656;0.345;0.445;0.253;0.088;0.422;0.367;0.325;0.116;0.345;0.656
7;Azorhizobium caulinodans;d-arabinose;+;?;0.532;0.025;?;?;?;?;?;?;?;?;?;?;0.863
;0.301;0.837;0.251;?;?;0.76;0.151;0.72;0.231;0.837;0.251


*/
typedef  pheno_notes_traits::Tr  Tr;
typedef  pheno_notes_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_home_brew_pheno Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

//typedef Tr::MyBlock  MyBlock;
typedef pheno_notes_params ParamGlob;

typedef string_tokenizer St;
typedef std::map<IdxTy,IdxTy> TermMap;
typedef std::map<IdxTy,TermMap> SpeciesMap;
typedef std::map<IdxTy,SpeciesMap> GenusMap;

//static 
class nested_map_iterator {


typedef GenusMap::const_iterator Gi;
typedef SpeciesMap::const_iterator Si;
typedef TermMap::const_iterator Ti;

public:
nested_map_iterator(const GenusMap & m, const St & st) : m_m(m),m_st(st) {Init();}
// TODO normally bool operator nad pre inc support but keeping these sep for now 
bool valid() const { return !m_done; } 
bool done() const { return m_done; } 
const StrTy & genus() const { return m_st((*m_gi).first); } 
const StrTy & species() const { return m_st((*m_si).first); } 
const StrTy & term() const { return m_st((*m_ti).first); } 
const IdxTy & count() const { return (*m_ti).second; } 

void inc()
{
BumpTI();

}
void seek(const StrTy & g, const StrTy & s , const StrTy & t) {Seek( g, s , t); }

private:
void Init()
{
m_done=false;
m_gi=m_m.begin();
if (m_gi==m_m.end()) {m_done=true; return; }
m_si=(*m_gi).second.begin();
if (m_si==(*m_gi).second.end()) {BumpGI(); return; }
m_ti=(*m_si).second.begin();
if (m_ti==(*m_si).second.end()) {BumpSI(); return; }

}
void Seek(const StrTy & g, const StrTy & s , const StrTy & t)
{
m_done=false;
if (g=="" ) { Init(); return; } 
m_gi=m_m.find(m_st[g]);
if (m_gi==m_m.end()) {m_done=true; return; }
if (s=="" ) { m_si= (*m_gi).second.begin(); } 
else m_si=(*m_gi).second.find(m_st[s]);
if (m_si==(*m_gi).second.end()) {m_done=true; return; }
if (t=="" ) { m_ti= (*m_si).second.begin(); } 
else m_ti=(*m_si).second.find(m_st[t]);
if (m_ti==(*m_si).second.end()) {m_done=true;  return; }

}




void BumpGI()
{
++m_gi;
if (m_gi==m_m.end()) {m_done=true; return; }
m_si=(*m_gi).second.begin();
if (m_si==(*m_gi).second.end()) {BumpGI(); return; }
m_ti=(*m_si).second.begin();
if (m_ti==(*m_si).second.end()) {BumpSI(); return; }

}
void BumpSI()
{
++m_si;
if (m_si==(*m_gi).second.end()) {BumpGI(); return; }
m_ti=(*m_si).second.begin();
if (m_ti==(*m_si).second.end()) {BumpSI(); return; }
}
void BumpTI()
{
++m_ti;
if (m_ti==(*m_si).second.end()) {BumpSI(); return; }
}

void Carry(const IdxTy lvl)
{
switch (lvl)
{
case 0: 
if (m_gi==m_m.end()) {m_done=true; return; }
m_si=(*m_gi).second.begin();
case 1:
if (m_si==(*m_gi).second.end()) {BumpGI(); return; }
m_ti=(*m_si).second.begin();
case 2:
if (m_ti==(*m_si).second.end()) {BumpSI(); return; }

} // switch 


} // carry 




bool m_done;
IdxTy m_bumps;

const GenusMap&  m_m;
const St & m_st;
Gi m_gi;
Si m_si;
Ti m_ti;
 
}; // nested_map_iterator;




public:
mjm_home_brew_pheno(): m_read_only(false)
,m_min_words(bad()),m_max_words(0),m_fields(0)

 {}
typedef  nested_map_iterator iterator; 
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
const IdxTy  size() const { return m_map.size(); }
void merge(const Myt & that ) { MergeMap(that); }
void load(const StrTy & fn ) { std::ifstream ifs(fn); load(ifs); }
void load(std::istream & is )
{
 IdxTy lines=0;
IdxTy f1=0; // bad();
IdxTy f2=1; // bad();
CommandInterpretter li(&is);
Setup(li,m_flp);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
//		MM_ERR(MMPR3(lines,sz,li.word(0)))
	//	m_table.add(li.words());
        if (m_fields==0) m_min_words=sz;
		m_fields+=sz;
        if (sz<m_min_words) m_min_words=sz;
        if (sz>m_max_words) m_max_words=sz;
		if ((sz>f1)&&(sz>f2))
		{
//			const StrTy key= mjm_strings::to_lower(li.word(f1));
//			Words gs;
//			Util::Split2(gs,key);
//			const IdxTy gssz=gs.size();
			const StrTy genus=li.word(f1);//(gssz>0)?gs[0]:StrTy("");
			const StrTy species=li.word(f2); // (gssz>1)?gs[1]:StrTy("");
//			if (gssz>2) MM_ONCE( " more words "<<MMPR2(gssz,key),)
//			const StrTy value= mjm_strings::to_lower(li.word(f2));
//			gs.clear();
//			Util::Split(gs,value);
			if ((f2+1)>=sz)
			{
			const StrTy & s=StrTy("blank"); // li.word(i);
			add(genus,species,s);
			//MM_ERR(" adding only "<<MMPR3(genus,species,s))
			//++m_map[m_st(genus)][m_st(species)][m_st(s)]; // =1;
			//++m_map[m_st(genus)][m_st(StrTy())][m_st(s)]; // =1;
			//++m_map[m_st(genus)][m_st(StrTy("*"))][m_st(s)]; // =1;
			}
			for (IdxTy i=(f2+1); i<sz; ++i) { 
			const StrTy & s=li.word(i);
			if (s.length()!=0) { add(genus,species,s); } 
	//		MM_LOOP(ii,gs){
//			MM_ERR(" adding "<<MMPR3(genus,species,s))
//			++m_map[m_st(genus)][m_st(species)][m_st(s)]; // =1;
//			++m_map[m_st(genus)][m_st(StrTy())][m_st(s)]; // =1;
//			++m_map[m_st(genus)][m_st(StrTy("*"))][m_st(s)]; // =1;
			}
//			m_index.add(key,lines);
		}
 //       m_lines.push_back(li.words());
		++lines;
    } // nextok()
MM_ERR(" done loading "<<MMPR3(m_min_words,m_max_words,m_fields)<< MMPR3(lines,size(),m_order.dump(0)))
}
void add(const StrTy & genus, const StrTy & species, const StrTy & s, const IdxTy flags=0)
{
	MM_ERR(" adding "<<MMPR3(genus,species,s))
	++m_map[m_st(genus)][m_st(species)][m_st(s)]; // =1;
	++m_map[m_st(genus)][m_st(StrTy())][m_st(s)]; // =1;
	++m_map[m_st(genus)][m_st(StrTy("*"))][m_st(s)]; // =1;
}


Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
{
Words w;
IdxTy gidx=m_st[genus];
IdxTy sidx=m_st[species];
MM_ERR(" find  entries for "<<MMPR4(gidx,sidx,genus,species))
//MM_ERR(" made it here ok ")
const auto & sp=m_map.find(gidx);
if (sp==m_map.end()) return w;
const auto & ip=(*sp).second.find(sidx);
if (ip==(*sp).second.end()) return w;
MM_ERR(" found entries for "<<MMPR2(genus,species))
StrTy s;
//MM_LOOP(ii,(*ip).second) {s=s+m_st((*ii).first); } w.push_back(s);
MM_LOOP(ii,(*ip).second) {for (IdxTy i=0; i<(*ii).second; ++i) w.push_back(m_st((*ii).first)); } 
return w;
}

void format_values(OsTy & os, const TermMap & m, const IdxTy flags) const 
{
//const bool print_kv=((flags&2)==0);
const bool skip_blank=((flags&4)==0);
const StrTy sep=" ";
MM_LOOP(ii,m)
{
const StrTy k=m_st((*ii).first);
if (skip_blank) if (k=="blank") continue; 
const IdxTy cnt=(*ii).second;
if (cnt!=0) os<<sep<<k;
if (cnt!=1) os<<":"<<cnt;
}
}

void dump(OsTy & os, const IdxTy flags=0) const 
{
const StrTy sep=" ";
const bool print_some_kv=((flags&1)==0);
MM_LOOP(ii,m_map)
{
const StrTy genus=m_st((*ii).first);
MM_LOOP(jj,(*ii).second)
{
const StrTy species=m_st((*jj).first);
os<<genus<<sep<<species;
if (print_some_kv) { format_values(os, (*jj).second,flags); }
os<<CRLF;
} // jj 
} // ii 
} // dump 

private:
void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,' ');
}
void MergeMaps(  const Myt  & b, const Myt & s)//  const
{
//if (d.size()!=0) { MM_ONCE(" may need to clear d before calling ",)} 
///m_map.clear();
// NB a simple equals operator would probably work for now...
MergeMap(b);
MergeMap(s);


}
void MergeMap(const Myt & that)
{

iterator mni(that.m_map,that.m_st);
while (mni.valid())
{
const StrTy & genus=mni.genus();
const StrTy & species=mni.species();
const StrTy & term=mni.term();
add(genus,species,term);

mni.inc();
} // mni


} // MergeMap



ParamGlob m_flp;
GenusMap m_map;
Ordering m_order;
// TODO FIXME make this a pointer and have a static shared for everyone 
St m_st;
bool m_read_only;
IdxTy m_min_words,m_max_words,m_fields;

}; // mjm_pro_pheno






///////////////////////////////////////////////////////////////




class mjm_pro_pheno
{
/*
cat ProTraits_precisionScores.txt  | more
Tax_ID;Organism_name;Phenotype;Minority;Original_label;Hamap_+;Hamap_-;MicrobeWi
ki_+;MicrobeWiki_-;Wikipedia_+;Wikipedia_-;Other_text_+;Other_text_-;PMC_article
s_+;PMC_articles_-;Pubmed_abstracts_+;Pubmed_abstracts_-;Proteome_composition_+;
Proteome_composition_-;Phyletic_profiles_+;Phyletic_profiles_-;Coexisting_microb
es_+;Coexisting_microbes_-;Gene_neighbourhood_+;Gene_neighbourhood_-;Codon_usage
_bias_+;Codon_usage_bias_-;Integrated_score_+;Integrated_score_-
7;Azorhizobium caulinodans;alkaline_phosphatase;-;?;0.58;0.033;?;?;?;?;?;?;?;?;0
.76;0.198;1.0;0.128;0.89;0.219;0.648;0.159;0.859;0.161;0.866;0.149;0.89;0.198
7;Azorhizobium caulinodans;arbutin;+;?;0.287;0.33;?;?;?;?;?;?;?;?;?;?;0.274;0.23
3;0.608;0.639;?;?;0.078;0.307;0.225;0.508;0.287;0.508
7;Azorhizobium caulinodans;cellobiose;-;?;0.008;1.0;?;?;?;?;?;?;?;?;?;?;0.156;0.
656;0.345;0.445;0.253;0.088;0.422;0.367;0.325;0.116;0.345;0.656
7;Azorhizobium caulinodans;d-arabinose;+;?;0.532;0.025;?;?;?;?;?;?;?;?;?;?;0.863
;0.301;0.837;0.251;?;?;0.76;0.151;0.72;0.231;0.837;0.251


*/
typedef  pheno_notes_traits::Tr  Tr;
typedef  pheno_notes_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_pro_pheno Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

//typedef Tr::MyBlock  MyBlock;
typedef pheno_notes_params ParamGlob;

typedef string_tokenizer St;
typedef std::map<IdxTy,IdxTy> TermMap;
typedef std::map<IdxTy,TermMap> SpeciesMap;
typedef std::map<IdxTy,SpeciesMap> GenusMap;

public:
mjm_pro_pheno(): m_read_only(false)
,m_min_words(bad()),m_max_words(0),m_fields(0)

 {}
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
const IdxTy  size() const { return 0; }

void load(const StrTy & fn )
{
std::ifstream  ifs(fn);  
load(ifs); 
}
void load(std::istream & is )
{
 IdxTy lines=0;
IdxTy f1=bad();
IdxTy f2=bad();
CommandInterpretter li(&is);
Setup(li,m_flp);
    if  (li.nextok()) m_order.add(li.words()); 
	const auto & col1= m_order[StrTy("Organism_name")];
	if (col1.size()!=0) f1=col1[0];	
	const auto & col2= m_order[StrTy("Phenotype")];
	if (col2.size()!=0) f2=col2[0];	
    while (li.nextok())
    {
        const IdxTy sz=li.size();
//		MM_ERR(MMPR3(lines,sz,li.word(0)))
	//	m_table.add(li.words());
        if (m_fields==0) m_min_words=sz;
		m_fields+=sz;
        if (sz<m_min_words) m_min_words=sz;
        if (sz>m_max_words) m_max_words=sz;
		if ((sz>f1)&&(sz>f2))
		{
			const StrTy key= mjm_strings::to_lower(li.word(f1));
			Words gs;
			Util::Split2(gs,key);
			const IdxTy gssz=gs.size();
			const StrTy genus=(gssz>0)?gs[0]:StrTy("");
			const StrTy species=(gssz>1)?gs[1]:StrTy("");
			if (gssz>2) MM_ONCE( " more words "<<MMPR2(gssz,key),)
			const StrTy value= mjm_strings::to_lower(li.word(f2));
			gs.clear();
			Util::Split(gs,value);
			MM_LOOP(ii,gs){
			//MM_ERR(" adding "<<MMPR2(genus,species))
			++m_map[m_st(genus)][m_st(species)][m_st(*ii)]; // =1;
			++m_map[m_st(genus)][m_st(StrTy())][m_st(*ii)]; // =1;
			}
//			m_index.add(key,lines);
		}
 //       m_lines.push_back(li.words());
		++lines;
    } // nextok()
MM_ERR(" done loading "<<MMPR3(m_min_words,m_max_words,m_fields)<< MMPR3(lines,size(),m_order.dump(0)))
}

Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
{
Words w;
IdxTy gidx=m_st[genus];
IdxTy sidx=m_st[species];
MM_ERR(" find  entries for "<<MMPR4(gidx,sidx,genus,species))
MM_ERR(" made it here ok ")
const auto & sp=m_map.find(gidx);
if (sp==m_map.end()) return w;
const auto & ip=(*sp).second.find(sidx);
if (ip==(*sp).second.end()) return w;
MM_ERR(" found entries for "<<MMPR2(genus,species))
StrTy s;
//MM_LOOP(ii,(*ip).second) {s=s+m_st((*ii).first); } w.push_back(s);
MM_LOOP(ii,(*ip).second) {for (IdxTy i=0; i<(*ii).second; ++i) w.push_back(m_st((*ii).first)); } 
return w;
}


private:
void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,';');
}

ParamGlob m_flp;
GenusMap m_map;
Ordering m_order;
St m_st;
bool m_read_only;
IdxTy m_min_words,m_max_words,m_fields;

}; // mjm_pro_pheno


class mjm_protrait
{
/*

INSTRUCTIONS
IJSEM_pheno_db.txt is the raw data.

If you are interested in contributing time to expanding the IJSEM bacterial phen
otype database, please contact Stuart Jones at sjones20@nd.edu.

SOME CURATION STEPS IN R:
i

~/junk/pheno$ ls -al
total 164836
drwxrwxr-x  2 marchywka marchywka      4096 2018-01-29 15:17 .
drwxr-xr-x 15 marchywka marchywka      4096 2018-03-20 10:28 ..
-rw-rw-r--  1 marchywka marchywka   1958534 2018-01-29 10:19 4272392.zip
-rw-rw-r--  1 marchywka marchywka   1955556 2018-01-29 10:14 IJSEM_pheno_db_v1.0.txt
-rw-rw-r--  1 marchywka marchywka     93660 2018-01-29 13:37 Protraits_names_only.txt
-rw-rw-r--  1 marchywka marchywka 164469969 2016-10-14 05:41 ProTraits_precisionScores.txt
-rw-rw-r--  1 marchywka marchywka      2662 2018-01-29 10:14 README.txt


*/
//typedef pheno_notes_traits Traits;
typedef  pheno_notes_traits::Tr  Tr;
typedef  pheno_notes_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_protrait Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

//typedef Tr::MyBlock  MyBlock;
typedef pheno_notes_params ParamGlob;
//typedef std::vector<StrTy > Words;

typedef mjm_ragged_table Ragged;

//template <class Tscalar, class Tval, class Tvec > //class mjm_1_to_some
typedef mjm_1_to_some<StrTy, IdxTy, std::vector< IdxTy > > IndexTy;

public:
mjm_protrait(): m_min_words(0),m_max_words(0),m_fields(0) {}
void load(const StrTy & fn )
{
std::map<StrTy,IdxTy> m;
//std::istream & is =  std::cin;
std::ifstream  is(fn); //  =  std::cin;

load(is);
}
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
//const Line & line(const IdxTy i ) const { return m_lines[i];}
//const IdxTy  size() const { return m_lines.size(); }
const IdxTy  size() const { return m_table.size(); }
//void add(const Line & li ) { m_lines.push_back(li); }
void load(std::istream & is )
{
 IdxTy lines=0;
IdxTy f1=bad();
IdxTy f2=bad();
// IdxTy skip=0;
//const  IdxTy interval=100000;
//const IdxTy w=2;
CommandInterpretter li(&is);
Setup(li,m_flp);
    if  (li.nextok()) m_order.add(li.words()); 
	const auto & col1= m_order[StrTy("Genus name")];
	if (col1.size()!=0) f1=col1[0];	
	const auto & col2= m_order[StrTy("species name")];
	if (col2.size()!=0) f2=col2[0];	

    while (li.nextok())
    {
        const IdxTy sz=li.size();
//		MM_ERR(MMPR3(lines,sz,li.word(0)))
		m_table.add(li.words());
        if (m_fields==0) m_min_words=sz;
		m_fields+=sz;
        if (sz<m_min_words) m_min_words=sz;
        if (sz>m_max_words) m_max_words=sz;
		if ((sz>f2)&&(sz>f1))
		{
			const StrTy key= mjm_strings::to_lower(li.word(f1)+StrTy(" ")+li.word(f2));
			m_index.add(key,lines);
		}
 //       m_lines.push_back(li.words());
		++lines;
    } // nextok()
MM_ERR(" done loading "<<MMPR4(m_min_words,m_max_words,m_fields,m_index.size())<< MMPR3(lines,size(),m_order.dump(0)))
} // load 

StrTy columns(const IdxTy f=0) const { return m_order.dump(f); }
template <class Td>
void enumerate_column(Td & m, const StrTy & nm) const
{
	IdxTy f2=bad();
	const auto & col2= m_order[nm];
	if (col2.size()!=0) f2=col2[0];	

MM_SZ_LOOP(i,m_table,tsz)
{
const auto & x=m_table.line(i); // (*ii);
const IdxTy sz=x.size();
if (sz>f2) ++m[x[f2]];

}

}

Words query( const StrTy & q, const StrTy & field,const IdxTy flags=0) const
{
	Words w;
	const bool just_results=((flags&1)!=0);
	const bool split_words=((flags&2)!=0);
	//const auto & col= m_order[q];
//	IdxTy f1=bad();
	IdxTy f2=bad();
//	if (col1.size()!=0) f1=col1[0];	
	const auto & col2= m_order[field];
	if (col2.size()!=0) f2=col2[0];	
	else { MM_ERR(" no columns for "<<MMPR2(q,field) ) return w;} 
	std::vector<IdxTy> locs= m_index[q];
//	MM_ERR(" query results "<<MMPR(locs.size()))
	if (locs.size()==0)
	{
		IdxTy fg=bad(); IdxTy fs=bad();
	 	locs= m_index.starting_with(q);
		const auto & colg= m_order[StrTy("Genus name")];
		if (colg.size()!=0) fg=colg[0];	
		const auto & cols= m_order[StrTy("species name")];
		if (cols.size()!=0) fs=cols[0];	
		const StrTy sep=" ";
	//	MM_ERR(" doing bounds search "<<MMPR4(locs.size(),fg,fs,f2))
		MM_LOOP(ii,locs)
		{
			const IdxTy line=(*ii);
			const auto & lv=m_table.line(line);
			StrTy x;
			const IdxTy sz=lv.size();
			if ( sz>fg) x=x+lv[fg];
			if ( sz>fs) x=x+sep+lv[fs];
			if ( sz>f2) x=x+sep+lv[f2];
			if (just_results){ 

 if (split_words) 
{
	Words w2;
	Split(w2,lv[f2]);
	MM_LOOP(jj,w2) { w.push_back(*jj); } 
}
else w.push_back(lv[f2]);

}
			else w.push_back(x);
		}
		return w;
	}
	MM_LOOP(ii,locs)
	{
		const auto & x=m_table.line((*ii));
		if (x.size()>=f2) {

//void Split( Words & w, const StrTy & s) const
 if (split_words) 
{
	Words w2;
	Split(w2,x[f2]);
	MM_LOOP(jj,w2) { w.push_back(*jj); } 
}
 else { w.push_back(x[f2]); } 

//MM_ERR(MMPR3(w.size(),x[f2],f2))
 }
		else w.push_back(StrTy());
	} // ii 
	return w; 
} // query 


void dump_words( const IdxTy flags) const  { dump_words(std::cerr,flags); } 
void dump_words(OsTy & os, const IdxTy flags) const 
{
const bool only_words=false;
const bool words_and_pos=!false;

MM_SZ_LOOP(i,m_table,sz)
{
const auto & line=m_table.line(i); // m_table[i];
IdxTy j=0;
MM_LOOP(jj,line)
{
std::vector<StrTy> split;
Split(split,(*jj));
MM_LOOP(kk,split)
{
if (only_words) { os<<(*kk)<<CRLF; } 
if (words_and_pos){ os<<j<<" "<<m_order[j]<<" "<<(*kk)<<CRLF; } 
} // kk 
++j;
} // jj 

} //i
}


private:

void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
li.readline_ok(false); li.use_alt_eol('\r',!false); li.set_split(1,'\t');
}
void Split( Words & w, const StrTy & s) const
{ 
 Util::Split(w,s); 
/*
const IdxTy sz=s.length();
const char * sc=s.c_str();
char c[sz+1];
//::memcpy(c,sc,sz);
IdxTy dptr=0;
IdxTy sptr=0;
while ( sc[sptr] != 0 )
{
const char cp=sc[sptr];  ++sptr;
const bool ender= ( cp==',')||(cp=='(')||(cp==')')||(cp=='"')||(cp=='|');
if (ender) 
{
if (dptr!=0){ c[dptr]=0; w.push_back(StrTy(c)); dptr=0;}}
else { c[dptr]=cp; ++dptr; } 
 
} // while sc sptr
if (dptr!=0) { c[dptr]=0; w.push_back(StrTy(c)); }
*/

}

ParamGlob m_flp;
Ordering m_order;
Ragged m_table;
IdxTy m_min_words,m_max_words,m_fields;
IndexTy m_index;
};   // mjm_protrait










class mjm_pheno_notes 
{
typedef  pheno_notes_traits::Tr  Tr;
typedef  pheno_notes_traits::util  Util;
typedef mjm_pheno_notes Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;


typedef pheno_notes_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_protrait Protrait;
typedef std::map<StrTy, Protrait> ProtraitMap;

typedef  mjm_pro_pheno ProPheno;
typedef std::map<StrTy, ProPheno> PhenoMap;

typedef  mjm_home_brew_pheno HbPheno;
typedef std::map<StrTy, HbPheno> HbMap;



//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;
typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
mjm_pheno_notes():m_dmel(new Dmel()) {Init();}
mjm_pheno_notes(int argc,char **_args) : m_dmel(new Dmel())
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<IdxTy(argc); ++i) args[i]=_args[i];
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
~mjm_pheno_notes()
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
if (confirm) {}
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

void command_modef(const char * fn)
{ std::ifstream fin(fn); CommandInterpretter li(&fin); command_mode(li); }
void command_mode() { CommandInterpretter li(&std::cin); command_mode(li); }
void command_mode(const StrTy & cmd) 
{ CommandInterpretter li; li.set(cmd,1); command_mode(li); }



CmdMap m_cmd_map;
CompMap m_comp_map;

 void cli_cmd( CliTy::list_type & choices,  const char * frag)
{
//MM_ERR("cli_cmd"<<MMPR(frag))
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);

}

}
 void cli_param( CliTy::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
const StrTy cmd=CliTy::word(StrTy(_cmd),0);
auto ii=m_comp_map.find(cmd);
if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag); 
}



void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_ragged_map[name].load(fn);
MM_ERR(MMPR2(m_ragged_map[name].size(),name))
}

void cmd_load_protrait(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_protrait_map[name].load(fn);
MM_ERR(MMPR2(m_protrait_map[name].size(),name))
}
void cmd_load_pheno(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_pheno_map[name].load(fn);
MM_ERR(MMPR2(m_pheno_map[name].size(),name))
}
void cmd_load_hbpheno(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_hbpheno_map[name].load(fn);
MM_ERR(MMPR2(m_hbpheno_map[name].size(),name))
}
void cmd_write_hbpheno(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
if (name.length()!=0) { 
std::ofstream os(name);
m_hbpheno_map[name].dump(os);
}
else 
m_hbpheno_map[name].dump(std::cout);

MM_ERR(MMPR2(m_hbpheno_map[name].size(),name))
}






void cmd_merge_hbpheno(Cip & cip , LocalVar & lv ) 
{
const StrTy dname=cip.p1;
const StrTy base=cip.p2;
const StrTy updates=cip.wif(3);
HbPheno & dp= m_hbpheno_map[dname];
HbPheno & bp= m_hbpheno_map[base];
HbPheno & up= m_hbpheno_map[updates];
MM_ERR("merge_hbpheno writing to "<<MMPR4(dname, base,updates,bp.size())<<MMPR(up.size()))

dp.merge(bp);
dp.merge(up);

} // cmd_merge_hbpheon






void cmd_query_pheno(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy genus=cip.p2;
const StrTy species=cip.wif(3);
ProPheno &   pm=m_pheno_map[name];//.load(fn);
//Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
const StrTy field="";
Words w=pm.query(genus,species,field);

MM_LOOP(ii,w) { MM_ERR((*ii)) }

MM_ERR(MMPR4(pm.size(),name,genus,species))
//MM_ERR(MMPR2(m_pheno_map[name].size(),name))
}

void cmd_query_hbpheno(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy genus=cip.p2;
const StrTy species=cip.wif(3);
HbPheno &   pm=m_hbpheno_map[name];//.load(fn);
//Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
const StrTy field="";
Words w=pm.query(genus,species,field);

MM_LOOP(ii,w) { MM_ERR((*ii)) }

MM_ERR(MMPR4(pm.size(),name,genus,species))
//MM_ERR(MMPR2(m_pheno_map[name].size(),name))
}




void cmd_query_protrait(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy query=cip.p2;
const StrTy field=cip.wif(3);
const Protrait & pro=m_protrait_map[name];
MM_ERR(" query_protrait "<<MMPR3(name,query,field))
auto w=pro.query(query,field);
MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR4(name, query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))

}
void cmd_enumerate_protrait(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy query=cip.p2;
const StrTy field=cip.wif(3);
const Protrait & pro=m_protrait_map[name];
MM_ERR(" enumerate_protrait "<<MMPR3(name,query,field))
std::map<StrTy,IdxTy> m;
pro.enumerate_column(m,field);
IdxTy i=0;
MM_LOOP(ii,m)
{
MM_ERR(MMPR3(i,(*ii).first,(*ii).second))
++i;
}
//MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR4(name, query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))

} // cmd_enumerate_protrait


void cmd_annotate(Cip & cip , LocalVar & lv ) 
{
const StrTy ttname=cip.p1;
const StrTy pname=cip.p2;
const StrTy field=cip.wif(3);
MM_ERR(" cmd_annotate "<<MMPR3(ttname,pname,field))
//void annotate_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
//annotate_tree(std::cout, ttname, pname, field, 0);
annotate_pro_tree(std::cout, ttname, pname, field, 0);

MM_ERR(" exist cmd_annotate "<<MMPR3(ttname,pname,field))
}

void cmd_pheno_annotate(Cip & cip , LocalVar & lv ) 
{
const StrTy ttname=cip.p1;
const StrTy pname=cip.p2;
const StrTy field=cip.wif(3);
MM_ERR(" cmd_annotate "<<MMPR3(ttname,pname,field))
//void annotate_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
//annotate_tree(std::cout, ttname, pname, field, 0);
annotate_pheno_tree(std::cout, ttname, pname, field, 0);

MM_ERR(" exist cmd_annotate "<<MMPR3(ttname,pname,field))
}

void cmd_hbpheno_annotate(Cip & cip , LocalVar & lv ) 
{
const StrTy ttname=cip.p1;
const StrTy pname=cip.p2;
const StrTy field=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(4));
MM_ERR(" cmd_annotate "<<MMPR4(ttname,pname,field,flags))
//void annotate_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
//annotate_tree(std::cout, ttname, pname, field, 0);
annotate_hbpheno_tree(std::cout, ttname, pname, field, flags);

MM_ERR(" exist cmd_annotate "<<MMPR3(ttname,pname,field))
}




void cmd_protrait_words(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const Protrait & pro=m_protrait_map[name];
MM_ERR("protrait_words"<< MMPR2(pro.size(),name))
pro.dump_words(0);
}

void cmd_tt(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.p1;
const StrTy name=cip.p2;
const StrTy p1=cip.wif(3);
TaxTree & tt=m_tax_trees[name];
MM_ERR(" cmd_tt "<<MMPR3(cmd,name,p1))
do { 
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }
} while (false); 


}






void cmd_write_svg(Cip & cip , LocalVar & lv ) 
{
m_char_mat.cmd_write_svg(cip,lv);
}

void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))

} 

}

void cmd_list(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_protrait_map) { MM_MSG("m_protrait_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();
/*
ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
ProtraitMap m_protrait_map;
PhenoMap m_pheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;
*/


}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("load-protrait")]=&Myt::cmd_load_protrait;
m_cmd_map[StrTy("load-pheno")]=&Myt::cmd_load_pheno;
m_cmd_map[StrTy("load-hbpheno")]=&Myt::cmd_load_hbpheno;
m_cmd_map[StrTy("write-hbpheno")]=&Myt::cmd_write_hbpheno;
m_cmd_map[StrTy("merge-hbpheno")]=&Myt::cmd_merge_hbpheno;
m_cmd_map[StrTy("query-protrait")]=&Myt::cmd_query_protrait;
m_cmd_map[StrTy("query-pheno")]=&Myt::cmd_query_pheno;
m_cmd_map[StrTy("query-hbpheno")]=&Myt::cmd_query_hbpheno;
m_cmd_map[StrTy("enumerate-protrait")]=&Myt::cmd_enumerate_protrait;
m_cmd_map[StrTy("protrait-words")]=&Myt::cmd_protrait_words;
m_cmd_map[StrTy("annotate")]=&Myt::cmd_annotate;
m_cmd_map[StrTy("pheno-annotate")]=&Myt::cmd_pheno_annotate;
m_cmd_map[StrTy("hbpheno-annotate")]=&Myt::cmd_hbpheno_annotate;
m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
TaxTree & tt = m_tax_tree;
StrTy local_label="tat";
//typedef void (Tsrc::* TargCmd)( ListTy & choices,  const char * frag);
//typedef void (Tsrc::* TargParam)( ListTy & choices, const char *  cmd, const char * frag);
m_cli.set_target(*this);
//void set_command_handler(TargCmd * p ) { m_targ_cmd=p; }
//void set_param_handler(TargParam * p ) { m_targ_param=p; }
m_cli.set_command_handler(&Myt::cli_cmd);
m_cli.set_param_handler(&Myt::cli_param);
//std::vector<StrTy> some;
//some.push_back(StrTy("load-tax-nodes"));
//m_cli.set_list(some);
m_cli.activate();
LocalVar mloc;
li.set_split(1,' ');
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

if (m_cmd_map.find(cmd)!=m_cmd_map.end())
{
 CommandInterpretterParam  cip(li); 
(this->*m_cmd_map[cmd])(cip,mloc);
continue;

}
const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
if (cmd=="about") { about();  continue; } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="print") { MM_MSG(li.line())  continue; } 
if (cmd=="err") { MM_ERR(li.line())  continue; } 
if (cmd=="status") { MM_STATUS(li.line())  continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 

//m_tax_tree.standard_commnds(cmd,p1,p2,li);


if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }

/*
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
*/

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
ss<<" mjm_trees_and_tables "<<__DATE__<<" "<<__TIME__<<CRLF;
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

std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();

}
IdxTy rooter(IdxTy & taxo, TaxTree & tt)
{
// this does not have a root if bacteria and archea have no common parent 
//const auto rootstt=tt.roots_unify();
//MM_ERR(" trying toots unify again danger will robinson")
const auto rootstt=tt.roots();
if (rootstt.size()==0) { MM_ERR(" no roots ") return bad() ; }
taxo=2;
IdxTy ridx=0;
MM_SZ_LOOP(i,rootstt,szr)
{
//if (tt.node_name(rootstt[i])=="bacteria") { ridx=i; break; } 
// non-parents point here lol 
// but does not propoagate lol 
if (rootstt[i]==bad()) { ridx=i; taxo=2; break; }
}
IdxTy node=rootstt[ridx];
MM_ERR("Setting roo node "<<MMPR4(node,ridx,tt.node_name(node),tt.node_info(node,256)))
return node;
}

#if 0 
template <class Ty>
void reduce(std::ostream &ss, const Ty & w , const IdxTy flags)
{

//if ( w.size()>0) {reduce(ss,w,0); } 
switch (flags)
{
case 1: 
{
std::map<StrTy,IdxTy> m;
std::map<D,StrTy> ms;
MM_LOOP(ii,w) { ++m[(*ii)]; } 
MM_LOOP(ii,m) { ms[-((*ii).second + 1.0/(ms.size()+10))]=(*ii).first;  } 
//MM_LOOP(ii,m) { ss<<" "<<(*ii).second<<":"<<(*ii).first; } 
MM_LOOP(ii,ms) { ss<<" "<<(*ii).second<<":"<<IdxTy(-(*ii).first); } 
break;
}
default:
{ MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
}

}
#endif
class pro_oracle
{
typedef  pheno_notes_traits::util  Util;
public:
pro_oracle(Myt & x, const StrTy & proname,const StrTy & f):  pro(x.m_protrait_map[proname]),field(f) {}
void lookup(std::ostream & ss,const std::vector<IdxTy> & tv,const std::vector<StrTy> & proq)
{
const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
//MM_ERR(" query_protrait "<<MMPR2(query,field))
auto w=pro.query(query,field,3);
if ( w.size()>0) {Util::reduce(ss,w,1); } // MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
else
{ // non-terminal nodes non GS should just use one name... doh 
const StrTy query2=proq[0];
w=pro.query(query2,field,3);
Util::reduce(ss,w,1);
//{MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
}
}
const Protrait & pro; // =m_protrait_map[proname];
const StrTy field;
}; //pro_oracle






class pheno_oracle
{
typedef  pheno_notes_traits::util  Util;
public:
pheno_oracle(Myt & x, const StrTy & proname,const StrTy & f):  pro(x.m_pheno_map[proname]),field(f) {}
void lookup(std::ostream & ss,const std::vector<IdxTy> & tv,const std::vector<StrTy> & proq)
{

//Words w=pm.query(genus,species,field);
//MM_LOOP(ii,w) { MM_ERR((*ii)) }
//MM_ERR(MMPR4(pm.size(),name,genus,species))

const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
MM_ERR(" query_protrait "<<MMPR3(query,field,proq.size()))
if ( proq.size()==0) return; 
const StrTy genus=(proq.size()>1)?proq[1]:StrTy();
const StrTy species=proq[0];
//Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
auto w=pro.query(genus,species,field,3);
if (!true) { 
if ((w.size())==0)
{
const StrTy blank="";
auto w=pro.query(species,blank,field,3);
ss<<w.size(); return ; } 

ss<<w.size(); return ; } 
//if (false) 
if (field.length()!=0)
{
MM_ERR(" looking for field "<<MMPR(field))
if ( w.size()>0) {Util::reduce(ss,w,1,field); }   
else
{  const StrTy query2=proq[0]; Util::reduce(ss,w,1,field); } 
}

if (!true){
if ( w.size()>0) {Util::reduce(ss,w,1); } // MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
else
{ // non-terminal nodes non GS should just use one name... doh 
const StrTy query2=proq[0];
//w=pro.query(query2,field,3);
Util::reduce(ss,w,1);
//{MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
}
}
}
const ProPheno & pro; // =m_protrait_map[proname];
const StrTy field;
}; //pheno_oracle

////////////////////////////////////////////

class null_oracle
{
typedef std::vector<StrTy> Sv;
typedef std::vector<IdxTy> Iv;
public:
null_oracle(Myt & x, const StrTy & proname,const StrTy & f) {}
null_oracle() {}
void lookup(Sv & tr,const Sv & proq){}
void lookup(std::ostream & ss,const Iv & tv,const Sv & proq){}



}; // null_oracle


class hbpheno_oracle
{
typedef  pheno_notes_traits::util  Util;
public:
hbpheno_oracle(Myt & x, const StrTy & proname,const StrTy & f):  pro(x.m_hbpheno_map[proname]),field(f) {}


void lookup(std::vector<StrTy> & tr,const std::vector<StrTy> & proq)
{
const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
MM_ERR(" query_homebrew  "<<MMPR3(query,field,proq.size()))
if ( proq.size()==0) return; 
const StrTy genus=(proq.size()>1)?proq[1]:StrTy();
const StrTy species=proq[0];
//Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
auto _w=pro.query(genus,species,field,3);
std::vector<StrTy> w;
MM_LOOP(ii,_w)
{
if (Util::fast_grep(*ii,field)) w.push_back(*ii);
}
MM_ERR(" query_homebrew  "<<MMPR4(w.size(),genus,species,field))
if (w.size()!=0) MM_ERR( MMPR(w[0]))
Util::reduce_unique(tr,w,0);

}

void lookup(std::ostream & ss,const std::vector<IdxTy> & tv,const std::vector<StrTy> & proq)
{

//Words w=pm.query(genus,species,field);
//MM_LOOP(ii,w) { MM_ERR((*ii)) }
//MM_ERR(MMPR4(pm.size(),name,genus,species))

const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
MM_ERR(" query_homebrew  "<<MMPR3(query,field,proq.size()))
if ( proq.size()==0) return; 
const StrTy genus=(proq.size()>1)?proq[1]:StrTy();
const StrTy species=proq[0];
//Words query( const StrTy & genus, const StrTy & species, const StrTy & field,const IdxTy flags=0) const
auto w=pro.query(genus,species,field,3);
MM_ERR(" query_homebrew  "<<MMPR4(w.size(),genus,species,field))
if (w.size()!=0) MM_ERR( MMPR(w[0]))
if (!true) { 
if ((w.size())==0)
{
const StrTy blank="";
auto w=pro.query(species,blank,field,3);
ss<<w.size(); return ; } 

ss<<w.size(); return ; } 
//if (false) 
if (field.length()!=0)
{
MM_ERR(" looking for field "<<MMPR(field))
if ( w.size()>0) {Util::reduce(ss,w,1,field); }   
else
{  const StrTy query2=proq[0]; Util::reduce(ss,w,1,field); } 
return; 
}

if (true){
MM_ERR(" wild card look")
if ( w.size()>0) {Util::reduce(ss,w,1); } // MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
else
{ // non-terminal nodes non GS should just use one name... doh 
const StrTy query2=proq[0];
//w=pro.query(query2,field,3);
Util::reduce(ss,w,1);
//{MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
}
}
}
const HbPheno & pro; // =m_protrait_map[proname];
const StrTy field;
}; //pheno_oracle


/////////////////////////////////////



void annotate_pro_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
{
//ParamGlob & flp=m_flp;
pro_oracle oracle(*this,proname,field);
annotate_tree(os,ttname,oracle,flags);
}
void annotate_pheno_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
{
//ParamGlob & flp=m_flp;
pheno_oracle oracle(*this,proname,field);
annotate_tree(os,ttname,oracle,flags);
}

void annotate_hbpheno_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
{
//ParamGlob & flp=m_flp;
hbpheno_oracle oracle(*this,proname,field);
annotate_tree(os,ttname,oracle,flags);
}





template <class To > 
void annotate_tree(std::ostream & os, const StrTy & ttname, To & oracle, const IdxTy flags)
{
const bool gs_only = ((flags&(1<<8))!=0);
StrTy spacerstring=" ";
StrTy hierlabel=" Hierarchy   ";
StrTy sep="\t"; // human?flp.human_sep():flp.data_sep(); // "\t";
TaxTree & tt=m_tax_trees[ttname]; // m_tax_tree;
IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
IdxTy nidx=0;
TaxTree::tree_hierarch_iterator ti(tt,node);
MM_ERR(" begin tree annotation ") 
while (ti.valid())
{
node=ti.node();
if (gs_only ) if (tt.nkids(node)!=0) { ti.inc(); continue; } 
Ss ss;
std::vector<IdxTy> tv;
std::vector<StrTy> proq;
tt.lineage(tv,node);
if (tv.size()==0 ) { MM_ERR(" node is not valid "<<MMPR(node))  ti.inc(); continue; } 
//IdxTy tvsz=tv.size();
MM_LOOP(ii,tv) { proq.push_back(tt.node_name(*ii)); } 
//MM_ERR(MMPR(proq.size()))
oracle.lookup(ss,tv,proq);
#if 0 
const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
MM_ERR(" query_protrait "<<MMPR2(query,field))
auto w=pro.query(query,field,3);
if ( w.size()>0) {reduce(ss,w,1); } // MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
else
{ // non-terminal nodes non GS should just use one name... doh 
const StrTy query2=proq[0];
w=pro.query(query2,field,3);
reduce(ss,w,1);
//{MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
}

#endif

//MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR3( query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))
//const IdxTy bumps=ti.bumps();
if (!gs_only)
{
StrTy indent=" ";
for( IdxTy i=0; i<ti.depth(); ++i) indent+=spacerstring; // StrTy(" ");
//nv.print_values(ss,node,digits,sep,flags);
StrTy nname=tt.node_name(node);
while (nname.length()<8) {nname=nname+StrTy(" "); }
os<<indent<<nname<<" "<<ss.str()<<CRLF;
}
else
{
const StrTy g=(proq.size()>1)? proq[1]:proq[0];
const StrTy sp=(proq.size()>1)?proq[0]:StrTy("*");
os<<g<<" "<<sp<<" "<<ss.str()<<CRLF; 

}


++nidx;
ti.inc();
} // ti
} // write_txt





void annotate_tree(std::ostream & os, const StrTy & ttname, const StrTy & proname, const StrTy & field, const IdxTy flags)
{
//ParamGlob & flp=m_flp;
StrTy spacerstring=" ";
StrTy hierlabel=" Hierarchy   ";
//const bool human=true;
StrTy sep="\t"; // human?flp.human_sep():flp.data_sep(); // "\t";
//IdxTy digits=1; // human?flp.human_digits():flp.data_digits(); // 1; // :8;

//NodeValues& nv=  m_values;
TaxTree & tt=m_tax_trees[ttname]; // m_tax_tree;
const Protrait & pro=m_protrait_map[proname];
//tt.sort_for_ui();
IdxTy taxo=2;
const IdxTy best_root=rooter(taxo,tt);
IdxTy node=best_root;
IdxTy nidx=0;
TaxTree::tree_hierarch_iterator ti(tt,node);
//const auto & order=nv.order();
//{ Ss ss; ss<<hierlabel; MM_LOOP(ii,order) ss<<" "<<(*ii); os<<ss.str()<<CRLF; }

while (ti.valid())
{
node=ti.node();
Ss ss;
std::vector<IdxTy> tv;
std::vector<StrTy> proq;
tt.lineage(tv,node);
// if the idea is no valid, it returns nothing
if (tv.size()==0 ) { MM_ERR(" node is not valid "<<MMPR(node))  ti.inc(); continue; } 
//IdxTy tvsz=tv.size();
MM_LOOP(ii,tv) { proq.push_back(tt.node_name(*ii)); } 
MM_ERR(MMPR(proq.size()))
const StrTy query=(proq.size()>1)?(proq[1]+StrTy(" ")+proq[0]):proq[0];
MM_ERR(" query_protrait "<<MMPR2(query,field))
auto w=pro.query(query,field,3);
if ( w.size()>0) {Util::reduce(ss,w,1); } // MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 
else
{ // non-terminal nodes non GS should just use one name... doh 
const StrTy query2=proq[0];
w=pro.query(query2,field,3);
Util::reduce(ss,w,1);
//{MM_SZ_LOOP(j,w,wsz){ ss<<" "<<w[j]; }  } 

}
//MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR3( query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))
//const IdxTy bumps=ti.bumps();
StrTy indent=" ";
for( IdxTy i=0; i<ti.depth(); ++i) indent+=spacerstring; // StrTy(" ");
//nv.print_values(ss,node,digits,sep,flags);
StrTy nname=tt.node_name(node);
while (nname.length()<8) {nname=nname+StrTy(" "); }
os<<indent<<nname<<" "<<ss.str()<<CRLF;
++nidx;
ti.inc();
} // ti
} // write_txt



// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
// use CommandInterpretterParam 
//const StrTy & cmd=li.word(0);

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

/*
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

*/

void DMel(const StrTy & e)
{
MM_ERR(e)
if (m_dmel!=0) {m_dmel->event("wtf",e); }
}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
	if (m_dmel!=0) {m_dmel->event("wtf",ss.str()); }
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; }
	if (m_dmel!=0) {m_dmel->event(code,ss.str()); }
    ss.str(StrTy(""));
}
void DumpDMel() // (OsTy & os)
{
if (m_dmel!=0) { MM_ERR(MMPR((m_dmel->string(1)))) } 
else { MM_ERR(" m_dmel is null ") }

}



bool m_done;

ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
ProtraitMap m_protrait_map;
PhenoMap m_pheno_map;
HbMap m_hbpheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;

}; //mjm_trees_and_tables



/////////////////////////////////////////////////////////

#ifdef  TEST_PHENO_NOTES__
int main(int argc,char **args)
{
typedef mjm_pheno_notes  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

