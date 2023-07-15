#ifndef MJM_STRING_TAX_H__
#define MJM_STRING_TAX_H__
 
#include "mjm_globals.h"
// for the presence absence vector 
//#include "mjm_ordering.h"
#include "mjm_char_mat.h"
//#include "mjm_part_iterators.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
//#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"
//#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
#include "../mjm_fasta_ii.h"

//#include "mjm_collections.h"
#include "mjm_svg_writer.h"

//#include "mjm_biom_hdf5.h"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <random>
// rand()
#include <stdlib.h>
#include <stdint.h>

/*
TODO FIXME
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 
/*


the fudding eog viewer crashed and I had no history records for
any of this cool sorting shot wtf FUDD 


2026  ./mjm_string_tax.out  -cmd "load-genus-scores foo linear_combs_nohits._txt 1" -cmd "set-param xlast 100" -cmd "set-param ylast 100" -cmd "write-svg-hm foo.svg foo 58" -quit
 2027  eog foo.svg


for loreiopsis, it buts the thing on right due to missing diagonal
connection 

 ./mjm_string_tax.out  -cmd "load-genus-strings-fasta foo xxxg" -cmd "set-param xlast 0" -cmd "set-param ylast 0" -cmd "write-svg-hm foo.svg foo 58" -quit

 ./mjm_string_tax.out  -cmd "load-genus-scores foo linear_combs_nohits._txt 1" -cmd "set-param xlast 100" -cmd "set-param ylast 100" -cmd "organism-pairs foo.svg foo 3" -cmd "write-svg-hm foo.svg foo 0"  -quit



*/


////////////////////////////////////////////////////////////////

class string_tax_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
string_tax_params( const StrTy & nm) : Super(nm) {}
string_tax_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy xfirst() const { return m_map.get_uint("xfirst",0); } // // 100;
IdxTy yfirst() const { return m_map.get_uint("yfirst",0); } // // 100;
IdxTy ylast() const { return m_map.get_uint("ylast",0); } // // 100;
IdxTy xlast() const { return m_map.get_uint("xlast",0); } // // 100;
IdxTy sort_iters() const { return m_map.get_uint("sort_iters",1); } // // 100;
IdxTy moment_sort_iters() const { return m_map.get_uint("moment_sort_iters",100); } // // 100;
StrTy run_file() const { return m_map.get_string("run_fule","mjm_string_tax.dump"); }
StrTy center_name() const { return m_map.get_string("center_name",""); }
StrTy highlight_name() const { return m_map.get_string("highlight_name",""); }
StrTy anote_name() const { return m_map.get_string("anote_name",""); }
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
ss<<"xfirst="<<xfirst()<<sep;
ss<<"yfirst="<<yfirst()<<sep;
ss<<"xlast="<<xlast()<<sep;
ss<<"ylast="<<ylast()<<sep;
ss<<"sort_iters="<<sort_iters()<<sep;
ss<<"moment_sort_iters="<<moment_sort_iters()<<sep;
ss<<"run_file="<<run_file()<<sep;
ss<<"center_name="<<center_name()<<sep;
ss<<"highlight_name="<<highlight_name()<<sep;
ss<<"anote_name="<<anote_name()<<sep;

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace string_tax_traits
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
//typedef std::vector<IdxTy > Locations;
//typedef mjm_string_ordering Ordering;

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

class mjm_string_graph
{
typedef  string_tax_traits::Tr  Tr;
typedef  string_tax_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_string_graph Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

//typedef Tr::MyBlock  MyBlock;
typedef string_tax_params ParamGlob;

typedef string_tokenizer St;
typedef std::map<StrTy,StrTy> SeqMap;
//typedef std::map<IdxTy,IdxTy> TermMap;
//typedef std::map<IdxTy,TermMap> SpeciesMap;
//typedef std::map<IdxTy,SpeciesMap> GenusMap;

typedef std::vector<D> ScoreVec;
typedef std::map<IdxTy,ScoreVec> ToScores;
typedef std::map<IdxTy, ToScores> FromScores;

typedef  std::map<IdxTy, IdxTy> Dm; 
typedef Dm::const_iterator DmCi;
typedef std::map < IdxTy, Dm  > Bm;
typedef Bm::const_iterator BmCi;
typedef std::map<IdxTy,IdxTy>  ToFromMap;




int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public:
mjm_string_graph(): m_read_only(false)
,m_min_words(bad()),m_max_words(0),m_fields(0)

 {}
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
const IdxTy  size() const { return m_scores.size(); }

void load(const StrTy & fn, const IdxTy &flags  )
{
std::ifstream  ifs(fn);  
load(ifs,flags); 
}


/*
mjm_string_seq.h4368  bests  j=9 n1=>NR_026093.1:Atopobium:parvulum bests|  1,>NR_026093.1:Atopobium:parvulum 1,>NR_113159.1:Atopobium:parvulum 1,>NR_102936.1:Atopobium:parvulum 0.941395,>NR_036819.1:Atopobium:rimae 0.941395,>NR_113038.1:Atopobium:rimae 0.905045,>NR_116936.1:Olsenella:umbonata 0.901335,>NR_036821.1:Olsenella:profusa
*/
bool is_digit_dot(const char & c) const { return ((c>='0')&&(c<='9'))||(c=='.'); } 
bool score_split(D& score, StrTy & rec, StrTy & genus, StrTy & species, const StrTy & s) const
{
const char score_sep=',';
const char name_sep=':';
const char kn_sep='=';
const IdxTy sz=s.length();
char c[sz+1];
memcpy(c,s.c_str(),sz);
c[sz]=0;
IdxTy i=0;
IdxTy sp=0;
IdxTy pc=0;
while (true)
{
char & ci=c[i];
if (ci==0)
{
			if (pc==3) { species=StrTy(c+sp); sp=i; ++pc;   } 
 			return (sp!=0);
}
if (ci==score_sep)
{
ci=0;
score=atof(c+sp);
sp=i+1;
++pc;
}
else 
{
// rec is ok 
//	if ((sp==0)&&(!is_digit_dot(ci))) return false; 
	if (ci==name_sep)
	{
		ci=0;
		switch (pc)
		{
			case 0: { rec=StrTy(c+sp); sp=i+2; ++pc; break; } 
			case 1: { rec=StrTy(c+sp); sp=i+1; ++pc; break; } 
			case 2: { genus=StrTy(c+sp); sp=i+1; ++pc; break; } 
			case 3: { species=StrTy(c+sp); sp=i+1; ++pc; break; } 
			default: { MM_ERR(" bad case "<<MMPR4(s,sp,pc,i)) } 
		} // pc
	}
	else if (ci==kn_sep) // not valid here but a different type 
	{
		ci=0;
	++pc;	
		rec=StrTy(c+sp);
		score=atof(c+i+1); // do not update anything just return with false later	
		genus=StrTy(c+i+1); // do not update anything just return with false later	
	} 

}

++i;
}



return true;
}


void load(std::istream & is, const IdxTy & flags  )
{
const bool get_from=((flags&1)!=0);
 IdxTy lines=0;
CommandInterpretter li(&is);
Setup(li,m_flp);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		const Words & w=li.words();
	//	MM_ERR(MMPR3(lines,sz,li.line()))
//		MM_ERR(MMPR3(lines,sz,li.word(0)))
	//	m_table.add(li.words());
        if (m_fields==0) m_min_words=sz;
		m_fields+=sz;
        if (sz<m_min_words) m_min_words=sz;
        if (sz>m_max_words) m_max_words=sz;
		if (sz<3) { ++lines; continue; }
		StrTy from="";
		for (IdxTy j=0; j<sz; ++j){
		StrTy rec="";
		StrTy genus="";
		StrTy species="";
		D score=0;
//bool score_split(D& score, StrTy & rec, StrTy & genus, StrTy & species, const StrTy & s) const
 if (score_split(score, rec, genus,  species,w[j]))
{
if ( from.length()==0){
//MM_ERR(" setting from "<<MMPR4(w[j],genus,rec,species))
//if (genus==StrTy("")) 
{ MM_ONCE(" consider setting flags  name format  "<<MMPR4(w[j],j,rec,species)<<MMPR(genus),) } 
 from= genus;
} else 
{m_scores[m_st(from)][m_st(genus)].push_back(score); 
//MM_ERR(MMPR3(m_st(from),m_scores.size(),flags)<<MMPR4(from,genus,species,score))
}}
else if ( rec==StrTy("j")) {
//MM_ERR(" j value "<<MMPR3(rec,score,species))

}
else {
//MM_ERR(" shotfudd  "<<MMPR(w[j])<<MMPR4(rec,score,species,genus))
if ( get_from)  if ( rec==StrTy("n1")) { from=genus;
MM_ONCE(" getting from from n1 entry  "<<MMPR(from),)
//MM_ERR(" getting from from n1 entry "<<MMPR2(rec,from))
 }
if (!get_from) { MM_ONCE(" getting name from parsed field",) }  
//MM_ERR(" unknown "<<MMPR2(j,w[j])<<  MMPR(rec)<<MMPR4(from,genus,species,score)) 
} 
} // j 
		++lines;
    } // nextok()
//MM_ERR(" done loading "<<MMPR3(m_min_words,m_max_words,m_fields)<< MMPR3(lines,size(),m_order.dump(0)))

//MM_ERR(MMPR(lines)<< stats())
MM_ERR(MMPR3(lines,m_scores.size(),size()))
} // load

/////////////////////////////////////////////////////////
// >genuseek2467770 10  Acaryochloris:3 Cuniculiplasma:1 Cyanobium:2 Defluviicoccus:1 Gloeobacter:2 Methanimicrococcus:2 Methanobrevibacter:1 Paenibacillus:1 Prochlorococcus:1 Synechococcus:1 Tetragenococcus:2

void load_genus_strings_fasta(const StrTy  & fn )
{
std::ifstream is(fn);
load_genus_strings_fasta(is);
}
void load_genus_strings_fasta(std::istream & is )
{
 IdxTy lines=0;
CommandInterpretter li(&is);
Setup(li,m_flp);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		const Words & w=li.words();
		if (sz<4) { ++lines; continue; }
		if (w[0].c_str()[0]!='>') continue;
		IdxTy lenhit=atof(w[1].c_str());
		std::vector<StrTy> genera;
		std::vector<D> scores;
		for (IdxTy j=2; j<sz; ++j){
		const IdxTy slen=w[j].length();
		char c[slen+1];
		memcpy(c,w[j].c_str(),slen+1);
		StrTy genus="";
		D score=0;
		for(IdxTy k=0; k<slen; ++k)
		{
			if (c[k]==':')
			{
				c[k]=0;
				genus=StrTy(c);
				score=atof(c+k+1)/(sz-2);
				break;
			}

		}
			if (genus!="") 
		//{ genera.push_back(genus); scores.push_back(score*lenhit); } 
		{ genera.push_back(genus); scores.push_back(lenhit*(1+0*score)); } 
} // j 
MM_LOOP(ii,genera)
{
IdxTy ix=0;
MM_LOOP(jj,genera)
{
if ((*ii)==(*jj)){ ++ix;  continue;}
//MM_ERR(" FUDD SHOT "<<MMPR2((*ii),(*jj)))
auto & v=m_scores[m_st(*ii)][m_st(*jj)];//.push_back(scores[ix]); 
if (v.size()==0) v.push_back(scores[ix]);
else v.back()+=scores[ix];
++ix;
}
}


		++lines;
    } // nextok()
MM_ERR(MMPR3(lines,m_scores.size(),size()))
} // load




////////////////////////////////////////////////////////

class placement_state
{
public:
class point { public: point(const IdxTy x, const IdxTy y ):m_x(x),m_y(y) {}
const IdxTy m_x,m_y;
}; // poingt
typedef std::vector<point> History; 
placement_state(const IdxTy nx, const IdxTy ny)
 :m_nx(nx),m_ny(ny),m_nz(1) ,m_sz(m_nx*m_ny*m_nz)
{
Init();
}
~placement_state() 
{
delete [] m_used;
}
const IdxTy&  w() const { return m_nx; } 
const IdxTy&  h() const { return m_ny; } 
const point & place_next()
{

IdxTy x=0;
IdxTy y=0;
FindAvail(x,y);
m_h.push_back(point(x,y));
return m_h.back();

}

const point & place_near(const IdxTy & tx, const IdxTy & ty)
{
IdxTy x=tx;
IdxTy y=ty;
FindAvailNear(x,y);
m_h.push_back(point(x,y));
return m_h.back();

}
private:
void Rmax(int * b, const int  r) const 
{
const IdxTy r2=2*r;
IdxTy i=0;
int  curse=r;
for (; i<=r2; ++i) b[i]=curse;
while (curse>(-r)) { --curse ;  b[i]=curse;++i; } 
for (; i<=(3*r2); ++i) b[i]=curse;
while (curse<(r)) { ++curse ;  b[i]=curse; ++i; } 
//Ss ss;
//ss<<MMPR(r);
//for (i=0; i<IdxTy(8*r); ++i) ss<<MMPR2(i,b[i]);
//MM_ERR(ss.str())
}
void FindAvailNear(IdxTy &x, IdxTy &y) const
{
const int  xc=x;
const int  yc=y;
int dx=int(x)-xc;
int dy=int(y)-yc;
// this could start at 0 but too cluttered 
int r=1; int pc=0;
int N=8*r;
pc=N-1;
int rmax=100;
int * rt = new int[8*rmax];
Rmax(rt,r);
IdxTy idx=Idx(x,y);
int off=2*r;
while (!Avail(idx))
{
// try to create a swirl with r offset 
if (pc==(N-1)) { r=r+1; off=2*r; if (x>(3*m_nx/4)) off+=4*r;
if ( x<(m_nx/4)) off=0;
if (y<(m_ny/4)) off+=2*r;
if (y>(3*m_ny/4)) off+=r;


if (r>=rmax) { delete [] rt; rmax=2*r; rt = new int[8*rmax]; } 
 N=8*r; pc=0; dx=r; dy=0; Rmax(rt,r); }
else {++pc;  dx=rt[(off+pc+0)%N]; dy=rt[(off+pc+(2*r))%N];
 } 

//MM_ERR(MMPR4(x,y,dx,dy)<<MMPR3(r,pc,N))
Clip(x,dx,m_nx); 
Clip(y,dy,m_ny); 
 idx=Idx(x,y);
}
Mark(idx);
delete [] rt;
}

void Clip(IdxTy &z, const int dx, const IdxTy lim) const
{
if ((int(z)+dx)>=int(lim)) z=(lim-1);
else {if ((int(z)+dx)<0) z=0;
else z=int(z)+dx;}

}
void FindAvail(IdxTy &x, IdxTy &y) const
{
static std::mt19937_64 mt;
 x=mt()%m_nx;
 y=mt()%m_ny;

IdxTy idx=Idx(x,y);
while (!Avail(idx))
{
 x=mt()%m_nx;
 y=mt()%m_ny;
 idx=Idx(x,y);
}
Mark(idx);

}

void Init()
{
m_used= new IdxTy[m_sz];
memset(m_used,0,m_sz*sizeof(IdxTy));
m_state=0;

}
IdxTy Idx(const IdxTy x, const IdxTy y) const { return x+y*m_nx; } 
bool Avail(const IdxTy x, const IdxTy y) const
{
return Avail(Idx(x,y)); 
}
bool Avail(const IdxTy idx) const
{
return (m_used[idx]==0); 
}
void  Mark(const IdxTy idx) const
{ m_used[idx]=1; }



// TODO FIXME use a vector, iterator cursor for this crap 
const IdxTy m_nx,m_ny,m_nz;
const IdxTy m_sz;
IdxTy m_state;
IdxTy * m_used;
History m_h;
}; // placement_state

class node_location
{
public:
IdxTy m_id, m_nx, m_ny; 
std::vector<IdxTy> m_conn;
std::vector<D> m_w;
}; // node_location

typedef std::map<IdxTy,node_location> NodePos;
#if 0 
void allocate(NodePos & np, const IdxTy & idx, IdxTy * used, const IdxTy nx, const IdxTy ny, IdxTy & x, IdxTy & y , IdxTy & dir  )
{
node_location nl;
static std::mt19937_64 mt;
if (np.find(idx)!=np.end()){ return; }  //  continue;
if (np.find(idx)==np.end()){ //  continue;
nl.m_id=idx;
IdxTy loc=y*nx+x;
while (used[loc]!=0) 
{
if (false) { if ( dir==0) { ++x; if (x>=nx) { x=0; ++y; } dir=1;}
else if ( dir==1) { ++y; if (y>=ny) { y-=2; }dir=2; }
else if ( dir==2) { if (x>0) --x; else ++x; ++y; if (y>=ny) { y-=2; }dir=3; }
else if ( dir==3) { if (y>0) --y; else ++y; ++x; if (x>=nx) { x-=2; }dir=0; }
}
D dnx=2.0*(D(mt()%nx)-nx/2)/nx;
dnx=dnx*dnx*dnx*dnx*dnx;
D dnxnx=nx*dnx;
if ((dnxnx+D(x))<0) x=0;
else {if ((dnxnx+D(x))>=nx) x=nx-1;
else x=IdxTy(dnxnx+x);
}
D dny=2.0*(D(mt()%ny)-ny/2)/ny;
dny=dny*dny*dny*dny*dny;
D dnyny=ny*dny;
if ((dnyny+D(y))<0) y=0;
else {if ((dnyny+D(y))>=ny) y=ny-1;
else y=IdxTy(dnyny+y); }
//y=mt()%ny;


loc=y*nx+x;
// this is a dummy insert as it inserts a COPY and we need to add to it but tell
// callees it is in there 
}
nl.m_nx=x;
nl.m_ny=y;
used[loc]=bad(); 
np[idx]=nl;
}
nl=np[idx];
MM_LOOP(jj,m_scores[idx])
{
IdxTy nidx=(*jj).first;
nl.m_conn.push_back(nidx);
allocate(np,nidx,used,nx,ny,x,y,dir);
D w=0;
MM_LOOP(kk,(*jj).second) { w+=(*kk); } 
nl.m_w.push_back(w);
} // nn 

np[idx]=nl;
//++x; if (x>=nx) { x=0; ++y; } 

}

void xxxposition_nodes(NodePos & np,const IdxTy nx, const IdxTy ny) 
{
IdxTy x=nx/2;
IdxTy y=ny/2;
IdxTy dir=0; 
IdxTy * used = new IdxTy[nx*ny];
::memset(used,0,nx*ny*sizeof(IdxTy));
MM_LOOP(ii,m_scores)
{
const IdxTy idx=(*ii).first;
allocate(np,idx,used,nx,ny,x,y,dir);

} // ii 
delete [] used; 
} // position_nodes

#endif

//////////////////////////////////////////////////////////

void place(NodePos & np, node_location & nl, const IdxTy idx, const IdxTy x, const IdxTy y)
{
nl.m_id=idx;
nl.m_nx=x;
nl.m_ny=y;
np[idx]=nl;

}
void allocate(NodePos & np,const IdxTy idx, placement_state & ps, const bool recurse) 
{
node_location nl;
// the connections are cyclic NOT parent child 
//if (np.find(idx)!=np.end()){ return; }  //  continue;
if (np.find(idx)==np.end()){ 
// see if it canbe placed near a placed component 
IdxTy putnear=bad();
IdxTy putnear2=bad();
D locmax=0;
D locmax2=0;
IdxTy tx=0;
IdxTy tx2=0;
IdxTy ty=0;
IdxTy ty2=0;
MM_LOOP(jj,m_scores[idx])
{
// the cyclic part, two nodes is a cycle doh 
IdxTy nidx=(*jj).first;
if (np.find(nidx)==np.end()){  continue; } 
D w=0;
MM_LOOP(kk,(*jj).second) { w+=(*kk); } 
//if (w>locmax) { tx=np[nidx].m_nx; ty=np[nidx].m_ny;  locmax=w; putnear=nidx; } 
if (w>locmax) {tx2=tx; ty2=ty; putnear2=putnear; locmax2=locmax;  tx=np[nidx].m_nx; ty=np[nidx].m_ny;  locmax=w; putnear=nidx; } 
else {if (w>locmax2)  tx2=np[nidx].m_nx; ty2=np[nidx].m_ny;  locmax2=w; putnear2=nidx; } 

}
const bool anywhere=(putnear==bad());
const bool anywhere2=(putnear2==bad());
int txc=.5*(tx+tx2);
int tyc=.5*(ty+ty2);
const int dx=int(tx2)-int(tx);
const int dy=int(ty2)-int(ty);
//if (dy==0) tyc+=.2*dx;
//else 
{ //D slope=-dx/dy; 
if (txc>int(ps.w()/2)) {txc+=3; } else { txc-=3; } 
if (txc<=int(ps.w()/2)) {txc-=3; } else { txc+=3; } 
if (tyc>int(ps.h()/2)) {tyc+=3; } else { tyc-=3; } 
if (tyc<=int(ps.h()/2)) {tyc-=3; } else { tyc+=3; } 

txc+=.2*dy; if (txc<0) txc=0;  if (txc>=int(ps.w())) txc=ps.w();
tyc-=.2*dx; if (tyc<0) tyc=0; if (tyc>int(ps.h())) tyc=ps.h();
}
const placement_state::point & psp= (anywhere)?ps.place_next():
(anywhere2? ps.place_near( tx,  ty):ps.place_near(txc,tyc));
place(np,nl,idx,psp.m_x,psp.m_y);
//nl.m_id=idx;
//nl.m_nx=psp.m_x;
//nl.m_ny=psp.m_y;
//np[idx]=nl;
} // put into a new location 
if (!recurse) return; 

nl=np[idx];
MM_LOOP(jj,m_scores[idx])
{
// the cyclic part, two nodes is a cycle doh 
IdxTy nidx=(*jj).first;
if (idx==nidx)
{
MM_ERR(" connection loop wtf "<<MMPR2(idx,nidx))
continue;
}
nl.m_conn.push_back(nidx);
allocate(np,nidx,ps,!true);
D w=0;
MM_LOOP(kk,(*jj).second) { w+=(*kk); } 
nl.m_w.push_back(w);
} // nn 

np[idx]=nl;
//++x; if (x>=nx) { x=0; ++y; } 


}

void position_nodes(NodePos & np,const IdxTy nx, const IdxTy ny) 
{
placement_state ps(nx,ny);
IdxTy x=nx/4;
IdxTy y=ny/4;
const D xcc=3*nx/8;
const D ycc=ny/2;
const IdxTy jmax=5;
const D rxpos=nx/4;
const D rypos=ny/4;
const D phi=2.0*M_PI/sqrt(D(jmax));
std::map<IdxTy,IdxTy> placed;
for (IdxTy j=0; j<jmax; ++j)
{
//x=((j&1)+1)*nx/4; 
//y=((j&2)+1)*ny/4; 
if (false)
{
switch (j)
{
case 0:{ x=nx/8+nx/2-1; y=ny/8+ny/2-1; break; } 
case 1:{ x=nx/8+nx/8; y=ny/8+ny/8+2; break; } 
case 2:{ x=7*nx/8+2; y=ny/8; break; } 
case 3:{ x=nx/8+nx/2-1; y=7*ny/8; break; } 
}
} //
{
x=IdxTy(xcc+rxpos*cos(phi*sqrt(D(j))));
y=IdxTy(ycc+rypos*sin(phi*sqrt(D(j))));


}


IdxTy cidx=bad();
MM_LOOP(ii,m_scores)
{
// TODO FIXME this is putting a bad eleent in doh.. 
// and could mess up itor too lol 
IdxTy thatsz=0;
auto jj=m_scores.find(cidx);
if (jj!=m_scores.end()) thatsz=(*jj).second.size();
//if ((*ii).second.size()>m_scores[cidx].size())
if ((*ii).second.size()>thatsz)
{ IdxTy i=(*ii).first;
	if (placed.find(i)==placed.end()) cidx=i; 
}
}
if (cidx==bad())
{
MM_ERR(" breaking at "<<MMPR(j))
 break;
}
node_location nl;
const placement_state::point & psp= ps.place_near(x,y);
MM_ERR(" placing "<<MMPR3(cidx,x,y))
place(np,nl,cidx,psp.m_x,psp.m_y);
placed[cidx]=1;
}

//allocate(np,cidx,ps,!true);
MM_LOOP(ii,m_scores)
{
const IdxTy idx=(*ii).first;
allocate(np,idx,ps,true);

} // ii 

} // position_nodes



void write_svg(const StrTy & fn)
{
std::ofstream os(fn);
write_svg(os);
}
void write_svg(std::ostream & os)
{
NodePos np;
const IdxTy xsz=60;
const IdxTy ysz=60;
position_nodes(np,xsz,ysz);

mjm_svg_writer sw;
const IdxTy xpitch=120;
const IdxTy ypitch=90;
const  D lblsz=2.0*.5*60;
const IdxTy tm=200;
const IdxTy lm=500;
const IdxTy bm=50;
const IdxTy rm=300;
const IdxTy wr=xsz; // 30;
const IdxTy hr=ysz; // 30;
//mjm_output_warp ow(xpitch,ypitch,w,h);
//mjm_circular_warp ow(xpitch,ypitch,w,h);
//ow_bounds(x,y,maxlevel,maxsz);
IdxTy xs=lm+rm+xpitch*wr;
IdxTy ys=bm+tm+ypitch*hr;
//ow.bounds(x,y,lc.m_level_sz);
os<<sw.start_text(" test foo",xs,ys);
os<<sw.frame_text("#00808080",xs,ys);
os<<CRLF;
D norm=0;
MM_LOOP(ii,np) {MM_LOOP(ww,(*ii).second.m_w) { if (norm<(*ww)) norm=(*ww);  }} 
for(IdxTy pass=0; pass<2; ++pass)
{

MM_LOOP(ii,np)
{
const StrTy nm=m_st((*ii).first);
const auto &no=(*ii).second;
D x=no.m_nx*xpitch+lm;
D y=no.m_ny*ypitch+tm;
const StrTy colr="red";
const D angle=30;
if (pass==1) 
{ os<<sw.gtext_text(nm,x,y,lblsz,colr,StrTy("middle"),angle); os<<CRLF; } 
auto ww=no.m_w.begin();
MM_LOOP(kk,no.m_conn)
{
const node_location & dest=np[*kk];
if (no.m_id==dest.m_id) continue;
D x1=dest.m_nx*xpitch+lm;
D y1=dest.m_ny*ypitch+tm;
//MM_ERR(" FUDD "<<MMPR(nm)<<MMPR4(no.m_nx,no.m_ny,dest.m_nx,dest.m_ny))
const D w=(*ww)/norm;
Ss ss;
ss<<w;
const D th=10;
//StrTy label_line_text(const StrTy &text, const D& x0, const D& y0, const D & x1, const D& y1,const D & sz, const StrTy & color="#00ffffff", const StrTy & anchor="left",const D & frac=.5) const
if (pass==0) {
os<<sw.label_line_text(ss.str(),x,y,x1,y1,30,StrTy("white"),"left",.3);
os<<CRLF; 
}
if (pass==0) { 
if (w<.9) {os<<sw.line_text(x,y,x1,y1, StrTy("blue"),th); os<<CRLF; }
else if (w<.95){ os<<sw.line_text(x,y,x1,y1, StrTy("green"),th); os<<CRLF;}
else { os<<sw.line_text(x,y,x1,y1, StrTy("red"),th); os<<CRLF; } 
}

++ww;
}
} // pass


} // np
os<<sw.end_text();
os<<CRLF;


}



/*
xv.push_back(x1);
yv.push_back(y1);
cv.push_back(col);
if (pass==0) { zmin=-log(zmin)/log(10); }
if (pass==0) { zmax=-log(zmax)/log(10); }
os<<sw.sparse_shade_rect_array_text(xv,yv,w,h, cv,1); os<<CRLF;
*/
const char  lut(const IdxTy p) const
{
static char c[16];
static bool init=false;
if (!init)
{
for(IdxTy i=0; i<10; ++i) c[i]='0'+i;
for(IdxTy i=0; i<6; ++i) c[i+10]='A'+i;

}
return c[p];
}

StrTy myhex(const D & x) const
{
Ss ss;
const IdxTy n=IdxTy(x*255);
ss<<(lut((n>>4)&15));
ss<<(lut((n>>0)&15));
return ss.str();
}

void color(StrTy & col,const D & red,const D & green,const D & blue,const IdxTy & flags)
const
{
Ss ss;
ss<<"#0";
ss<<myhex(red);
ss<<myhex(green);
ss<<myhex(blue);
col=ss.str();
}

StrTy color(const D & val, const D & zmin, const D & zmax) const
{
D rho=(val-zmin)/(zmax-zmin);
D red=1.0-2.0*rho;
if (red<0) red=0;
D blue=2.0*rho-1.0;
if (blue<0) blue=0;
D green=(rho-.75)*(rho-.25);
if (green>0) green=0;
else { green=sqrt(-16.0*green); }
StrTy col="";
color(col,red,green,blue,0);
return col;

}





template <class Tfudd> void ASSFUDD(FromScores & _scores,  FromScores & rev_map, Tfudd & lmap_from,Tfudd & lmap_to, D & zmin, D & zmax, const bool max_score)
{
MM_LOOP(ii,_scores)
{const auto &no=(*ii).second;
const IdxTy src=(*ii).first;
// NOTE size starts at 1 maybe? Sure? 
if (lmap_from.find(src)==lmap_from.end()) lmap_from[src]=lmap_from.size();

MM_LOOP(jj,no)
{D score=0;
const IdxTy dest=(*jj).first;
if (max_score)
{MM_LOOP(kk,(*jj).second) 
{  rev_map[dest][src].push_back(*kk); 
if ((*kk)>score) score=(*kk); } 
}
else{
MM_LOOP(kk,(*jj).second) {  rev_map[dest][src].push_back(*kk); score+=(*kk); } 
}
if (score==0) continue; 
//if (lmap.find(dest)==lmap.end()) lmap[dest]=lmap.size();
if (lmap_to.find(dest)==lmap_to.end()) lmap_to[dest]=lmap_to.size();

if (zmax<0) zmax=score;
if (zmax<score) zmax=score;
if (zmin<0) zmin=score;
if (zmin>score) zmin=score;
}
}
}

//template <class Tm>
static bool sort_map_connects(const FromScores  & m, const int & a, const int & b) 
{
auto ia=m.find(a);
auto ib=m.find(b);
if (ia==m.end()) { MM_ERR(" sort problem "<<MMPR2(a,b)) }
if (ib==m.end()) { MM_ERR(" sort problem "<<MMPR2(a,b)) }
//return (m_scores[a].size() > m_scores[b].size()); 
return ((*ia).second.size() > ((*ib).second.size()));

}



static bool sort_map_amp(const FromScores  & m, const bool max_score, const int & a, const int & b) 
{
auto ia=m.find(a);
auto ib=m.find(b);
if (ia==m.end()) { MM_ERR(" sort problem "<<MMPR2(a,b)) }
if (ib==m.end()) { MM_ERR(" sort problem "<<MMPR2(a,b)) }
D scorea=0;
//MM_LOOP(ii,m_scores[a]) { 
MM_LOOP(ii,(*ia).second) { D sx=0;
{MM_LOOP(kk,(*ii).second){ if (max_score) {if ((*kk)>sx) sx=(*kk); } else  scorea+=(*kk);} } scorea+=sx; }

D scoreb=0;
//MM_LOOP(ii,m_scores[b]) { 
MM_LOOP(ii,(*ib).second) { D sx=0;
//{MM_LOOP(kk,(*ii).second)  scoreb+=(*kk); }  } 
{MM_LOOP(kk,(*ii).second){ if (max_score) {if ((*kk)>sx) sx=(*kk); } else  scoreb+=(*kk);} } scoreb+=sx;  }

return (scorea > scoreb); }





typedef std::vector<IdxTy> XYord; //  xord,yord;
void  moments(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const IdxTy sz=xord.size();
//std::vector<D> mox(sz);
std::map<IdxTy,D> mox; 
MM_SZ_LOOP(i,xord,xordsz)
{
const bool normalize=true;
D mo=0; D norm=0;
const auto sii=s.find(xord[i]);
const auto & mv=(*sii); // s[xord[i]];
// TODO this needs the maps not tohe orders as this is sparse and nwo slow 
MM_SZ_LOOP(j,yord,yordsz)
{
const auto  jf=mv.second.find(yord[j]);
if ((jf)!=mv.second.end()){
D str=0;
const bool just_count=true;
if (just_count) str=1; else { 
MM_LOOP(jj,(*jf).second) { str+=(*jj); } }
// this looks better 
mo+=str/(j+1);
// mo+=str/(yord[j]+1);
norm+=str;
}
} // jj 
if (normalize&&((norm!=0))) mox[xord[i]]=mo/norm; else mox[xord[i]]=mo;
} // ii 
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ return mox[a]>mox[b]; });
} // moments 
/////////////////////////////////////
void  centroid(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
const bool rev=((flags&1)!=0);
//const IdxTy sz=xord.size();
//std::vector<D> mox(sz);
std::map<IdxTy,D> mox; 
MM_SZ_LOOP(i,xord,xordsz)
{
D mo=0; D norm=0;
const auto sii=s.find(xord[i]);
const auto & mv=(*sii); // s[xord[i]];
// TODO this needs the maps not tohe orders as this is sparse and nwo slow 
MM_SZ_LOOP(j,yord,yordsz)
{
const auto  jf=mv.second.find(yord[j]);
if ((jf)!=mv.second.end()){
D str=0;
const bool just_count=true;
if (just_count) str=1; else{
MM_LOOP(jj,(*jf).second) { str+=(*jj); }  }

// original 
 mo+=str*(j+1);
//if (!rev) mo+=str*(yordsz-j);
//else mo+=str*(j+1);


// mo+=str*(yord[j]+1);
norm+=str;
}
} // jj 
if (norm!=0) mox[xord[i]]=mo/norm; else mox[xord[i]]=0;
} // ii 
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
// FUDDING deos not matter FUDD 	
	//  return mox[a]<mox[b];
	if (rev) return mox[a]<mox[b];
	return mox[a]>mox[b];

 });
} // centroid 



/////////////////////////////////
void  runlen(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool rev=((flags&1)!=0);
//const IdxTy sz=xord.size();
//std::vector<D> mox(sz);
std::map<IdxTy,D> mox; 
MM_SZ_LOOP(i,xord,xordsz)
{
IdxTy rpos=0;
IdxTy rlen=0;
IdxTy rbest=0;
IdxTy pbest=0;
const auto sii=s.find(xord[i]);
const auto & mv=(*sii); // s[xord[i]];
// TODO this needs the maps not tohe orders as this is sparse and nwo slow 
MM_SZ_LOOP(j,yord,yordsz)
{
const auto  jf=mv.second.find(yord[j]);
if ((jf)!=mv.second.end()){
//if (rlen==0) rpos=yord[j];
if (rlen==0) rpos=j;
++rlen;
} else {
if (rlen>rbest) { rbest=rlen; pbest=rpos; }
rlen=0; 
}
} // jj 

mox[xord[i]]=pbest;
//if (norm!=0) mox[xord[i]]=mo/norm; else mox[xord[i]]=0;
} // ii 
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
//if (rev) return mox[a]>mox[b]; 

return mox[a]<mox[b]; }


);
} // centroid 



////////////////////////////////



/////////////////////////////////////////////////////
void  ranks(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
/*
const IdxTy sz=xord.size();
//std::vector<D> mox(sz);
std::map<IdxTy,D> mox; 
MM_SZ_LOOP(i,xord,xordsz)
{
D mo=0;
const auto sii=s.find(xord[i]);
const auto & mv=(*sii); // s[xord[i]];
// TODO this needs the maps not tohe orders as this is sparse and nwo slow 
MM_SZ_LOOP(j,yord,yordsz)
{
const auto  jf=mv.second.find(yord[j]);
if ((jf)!=mv.second.end()){
D str=0;
MM_LOOP(jj,(*jf).second) { str+=(*jj); } 
 mo+=str/(j+1);
}
} // jj 
mox[xord[i]]=mo;
} // ii */
const bool lt=true;
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
// a and b are values in xord to be compared by their maps entries for yord
//MM_ONCE(" danger will robinson kluge lt eq ",)
//if (a==b) return true; 
const auto sia=s.find(a);
const auto sib=s.find(b);
const bool sia_end=(sia==s.end());
const bool sib_end=(sib==s.end());
if (sia_end||sib_end)
{
MM_ERR(" sorting error "<<MMPR2(a,b))
}
const auto & ma=(*sia).second;
const auto & mb=(*sib).second;
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const auto dib=mb.find(v);
const bool dia_end=(dia==ma.end());
const bool dib_end=(dib==mb.end());
// there is an entry for a not b it is earlier 
if (!dia_end&&dib_end) return lt;
if (dia_end&&!dib_end) return !lt;
if (!dia_end&&!dib_end) 
{
D s1=0; D s2=0;
MM_LOOP(jj,(*dia).second) { s1+=(*jj); } 
MM_LOOP(jj,(*dib).second) { s2+=(*jj); } 
if (s1==s2) continue; // return false;
if (s1>s2) return lt;
return !lt;
}

}

// these are equal..
return false; 



});
} // moments 

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void  maxscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool lt=true;
std::map<IdxTy,D> vals;
MM_LOOP(ii,xord)
{
D s1=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
if( !dia_end) { MM_LOOP(jj,(*dia).second) { if ((*jj)>s1)s1=(*jj); }  } 

} // j 
vals[(*ii)]=s1;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
return vals[a]>vals[b];
#if 0 
// a and b are values in xord to be compared by their maps entries for yord
//MM_ONCE(" danger will robinson kluge lt eq ",)
//if (a==b) return true; 
const auto sia=s.find(a);
const auto sib=s.find(b);
const bool sia_end=(sia==s.end());
const bool sib_end=(sib==s.end());
if (sia_end||sib_end)
{
MM_ERR(" sorting error "<<MMPR2(a,b))
}
D s1=0; D s2=0;
const auto & ma=(*sia).second;
const auto & mb=(*sib).second;
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const auto dib=mb.find(v);
const bool dia_end=(dia==ma.end());
const bool dib_end=(dib==mb.end());
// there is an entry for a not b it is earlier 
{
if( !dia_end) { MM_LOOP(jj,(*dia).second) { if ((*jj)>s1)s1=(*jj); }  } 
if (!dib_end) { MM_LOOP(jj,(*dib).second) { if ((*jj)>s2)s2=(*jj); }  } 
}

}
return (s1>s2); 
// these are equal..
//return false; 
#endif

});



} // maxscores




void  maxscoresxxxhhhhhhhhhhh(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool lt=true;
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
// a and b are values in xord to be compared by their maps entries for yord
//MM_ONCE(" danger will robinson kluge lt eq ",)
//if (a==b) return true; 
const auto sia=s.find(a);
const auto sib=s.find(b);
const bool sia_end=(sia==s.end());
const bool sib_end=(sib==s.end());
if (sia_end||sib_end)
{
MM_ERR(" sorting error "<<MMPR2(a,b))
}
D s1=0; D s2=0;
const auto & ma=(*sia).second;
const auto & mb=(*sib).second;
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const auto dib=mb.find(v);
const bool dia_end=(dia==ma.end());
const bool dib_end=(dib==mb.end());
// there is an entry for a not b it is earlier 
//if (!dia_end&&dib_end) return lt; // TODO kluge order? 
//if (dia_end&&!dib_end) return !lt;
//if (!dia_end&&!dib_end) 
{
if( !dia_end) { MM_LOOP(jj,(*dia).second) { if ((*jj)>s1)s1=(*jj); }  } 
if (!dib_end) { MM_LOOP(jj,(*dib).second) { if ((*jj)>s2)s2=(*jj); }  } 
//if (s1==s2) continue; // return false;
//if (s1>s2) return lt; // kluge
//return !lt;
}

}
return (s1>s2); 
// these are equal..
//return false; 

});
} // maxscores

void  totscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool lt=true;
std::map<IdxTy,D> vals;
MM_LOOP(ii,xord)
{
D s1=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 

} // j 
vals[(*ii)]=s1;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
return vals[a]>vals[b];
#if 0 
#endif

});

} // maxscores

/////////////////////////////////////////

void  moscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
const bool fwdmo=((flags&1)==0);
//const bool lt=true;
std::map<IdxTy,D> vals;
MM_LOOP(ii,xord)
{
D s1=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
if (fwdmo)
{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(yordsz-j)*(*jj); }  } 
} else { 
if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(j)*(*jj); }  } 
}
} // j 
vals[(*ii)]=s1;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
return vals[a]>vals[b];
#if 0 
#endif

});

} // maxscores


void  peakscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool fwdmo=((flags&1)==0);
//const bool lt=true;
std::map<IdxTy,D> vals;
MM_LOOP(ii,xord)
{
D peak=0;
IdxTy loc=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
D s1=0;
{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 
if (s1>peak) { peak=s1; loc=j; } 
//if (fwdmo)
//{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(yordsz-j)*(*jj); }  } 
//} else { 
//if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(j)*(*jj); }  } 
}
} // j 
vals[(*ii)]=loc;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
return vals[a]>vals[b];
#if 0 
#endif

});

} // peakscores


void  crossscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
const bool fwdmo=((flags&1)==0);
//const bool lt=true;
std::map<IdxTy,D> vals;
const D ycent=.5*yord.size();
MM_LOOP(ii,xord)
{
D sum=0;
//IdxTy loc=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
//D s1=0;
//{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 
{if( !dia_end) { sum+=(j-ycent);  } 
//if (s1>peak) { peak=s1; loc=j; } 
//if (fwdmo)
//{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(yordsz-j)*(*jj); }  } 
//} else { 
//if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(j)*(*jj); }  } 
}
} // j 
vals[(*ii)]=sum;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
if (fwdmo) return vals[a]>vals[b];
return vals[a]<vals[b]; // orig cdoe
#if 0 
#endif

});

} // crossscores


void  firstscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool fwdmo=((flags&1)==0);
//const bool lt=true;
std::map<IdxTy,D> vals;
//const D ycent=.5*yord.size();
MM_LOOP(ii,xord)
{
D sum=0;
//IdxTy loc=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const bool dia_end=(dia==ma.end());
//D s1=0;
//if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 
if( !dia_end) { sum=j; break;  } 
//if( !dia_end) { sum=s1; break;  } 
//if (s1>peak) { peak=s1; loc=j; } 
//if (fwdmo)
//{if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(yordsz-j)*(*jj); }  } 
//} else { 
//if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(j)*(*jj); }  } 

} // j 
vals[(*ii)]=sum;
} // ii

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
//if (fwdmo) return vals[a]>vals[b];
//return vals[a]>vals[b]; // orig cdoe
return vals[a]<vals[b]; // orig cdoe
#if 0 
#endif

});

} // firstscores

class cluster
{
public:
static const IdxTy bad() { return ~0; } 
cluster(): m_1(bad()),m_2(bad()) {}
cluster(const IdxTy i)
: m_1(bad()),m_2(bad()) { add(i);  }

cluster(const IdxTy i, const IdxTy j)
: m_1(bad()),m_2(bad()) { add(i); add(j); }

void add(const IdxTy i) { m_set.push_back(i); }
void bridge(const IdxTy i) 
{ if (m_1==bad()) m_1=i; else if (m_2==bad()) m_2=i; }

IdxTy m_1,m_2;
std::vector<IdxTy> m_set;

}; // cluster
void  cprobscores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
std::map<IdxTy,IdxTy> ycnts,himap,xhisto;
const IdxTy szx=xord.size();
const IdxTy szy=yord.size();
const IdxTy maxma=300;
std::vector<IdxTy> yprob(szy);
std::map<IdxTy,IdxTy> hitbits;
MM_LOOP(ii,xord)
{
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// ycnts is the frequnece of hits to the index or y value 
MM_LOOP(kk,ma) { ++ycnts[(*kk).first]; } 
} // xord
MM_ERR("ycnts made "<<MMPR3(szx,szy,ycnts.size()))
//  yprob is yvalues appearoing in ycnts but this could be a problem !!
 IdxTy i=0; MM_LOOP(ii,ycnts) { yprob[i]=(*ii).first; ++i;  } 
// in any case then this ends up as y values, not yord idx, in descending cnt order 
std::sort(yprob.begin(),yprob.end(),
[&](const int& a, const int& b) { return ycnts[a]>ycnts[b]; });

MM_ERR("yprob sort ")
const IdxTy jmax=(szy<8)?szy:8;
for(IdxTy j=0; j<jmax; ++j) 
{
MM_ERR(MMPR4(j,yprob[j],ycnts[yprob[j]],m_st(yprob[j])))

}
D frac=0;
const bool make_pairs=!false&&(((szx|szy)&(~0x0FFFF))==0);
std::map<IdxTy,IdxTy> pp;
std::vector<IdxTy> ppv;
for ( i=0; i<szx; ++i)
{
IdxTy & hi=hitbits[xord[i]];
hi=0; 
const auto sia=s.find(xord[i]);
const auto & ma=(*sia).second;
++xhisto[ma.size()];
frac+=ma.size();
for(IdxTy j=0; j<jmax; ++j) 
	if (ma.find(yprob[j]) != ma.end()) hi|=(1<<j);
++himap[hi];

if (make_pairs&&(ma.size()<maxma))
{
const IdxTy sp=16;
IdxTy ikk=0;
MM_LOOP(kk,ma)
{
ikk=(*kk).first;
IdxTy ikki=ikk+1;
for(auto  kkp=kk; kkp!=ma.end(); ++kkp)
{
if (kkp==kk)continue; //  ++kkp; 
// make sure score maps are non-zero
ikki=(*kkp).first;

const IdxTy pploc=(ikk<<sp)|ikki;
++pp[pploc];
++ikki;
}
++ikk;
} // kk

} // make_pairs

} // i 
frac=frac/(szx*szy);
MM_ERR("himap made  "<<MMPR2(himap.size(),frac))
if (make_pairs)
{
const IdxTy sp=16;
const IdxTy mask=(1<<sp)-1;
const IdxTy nmin=10;
const IdxTy nlim=szy>>1;
MM_LOOP(ip,pp) {
const IdxTy & pploc=(*ip).first;
const IdxTy i=pploc>>sp; IdxTy j=pploc&(mask);  
const IdxTy n=(*ip).second; 
if ((n>=nmin)&&(n<nlim)) ppv.push_back(pploc);
// MM_ERR("ppairs " << MMPR3(i,j,n)<<MMPR2(m_st(i),m_st(j)))
 } 
// ppv now contains the codes for some pairs 
// that can be extendned
typedef std::vector<IdxTy> Sh;
typedef std::vector<Sh> Shv;
// the key is  xord  vale
typedef std::map<IdxTy, Shv> ShMap;
ShMap shmap;
IdxTy maxsh=0;
IdxTy maxshlim=21000;
MM_ERR("make shmap "<<MMPR(ppv.size()))
MM_SZ_LOOP(iv,ppv,ppvsz)
{
const IdxTy & pploc=ppv[iv];
const IdxTy i=pploc>>sp; IdxTy j=pploc&(mask);  
for (IdxTy jv=iv+1; jv<ppvsz; ++jv)
{
const IdxTy & pplocp=ppv[jv];
const IdxTy ip=pplocp>>sp; IdxTy jp=pplocp&(mask);  
for ( IdxTy k=0; k<szx; ++k)
{
const IdxTy idx=xord[k];
// this makes a mess and is slow not sure it is useful. 
//if (shmap[idx].size()>1000) { MM_ERR(" too big shmap entry "<<MMPR(m_st(idx)))
//continue; } 
const auto sia=s.find(idx);
const auto & ma=(*sia).second;
auto & shidx= shmap[idx];
if (shidx.size()>maxshlim) { continue; }  
if ((ma.size()>=maxma)) continue; 
if (ma.find(i)==ma.end()) continue;
if (ma.find(j)==ma.end()) continue;
if (ma.find(ip)==ma.end()) continue;
if (ma.find(jp)==ma.end()) continue;
shidx.push_back( Sh(4)); 
if (shidx.size()>maxsh) { maxsh=shidx.size(); }  
auto & d= shidx.back(); 
d[0]=i; d[1]=j; d[2]=ip; d[3]=jp;
} // k 
} // j 
} // i 
MM_ERR("done with 4word  "<<MMPR(maxsh))
const bool make_longer=!true;
if ( make_longer)
{
for(IdxTy longp=0; longp<4; ++longp)
{
MM_ERR("make longer "<<MMPR(longp))
for ( IdxTy k=0; k<szx; ++k)
{
const IdxTy xval=xord[k];
if ( shmap.find(xval) == shmap.end()) continue;
const auto sia=s.find(xval);
const auto & ma=(*sia).second;

Shv shv;
auto & s=shmap[xval];

if (s.size()>=szy)
{
MM_ERR(" too big already "<<MMPR4(xval,m_st(xval),k,s.size()))
continue;
}

MM_LOOP(ii,s)
{
//orginal destroyed 
// this either needs to pick up reorderd dups or give up lol 
MM_SZ_LOOP(iv,ppv,ppvsz)
{
Sh  existing=(*ii);
const IdxTy & pploc=ppv[iv];
const IdxTy i=pploc>>sp; IdxTy j=pploc&(mask);  
if (j<=(*ii).back()) continue; 
if (ma.find(i)==ma.end()) continue;
if (ma.find(j)==ma.end()) continue;
bool add=true;
MM_LOOP(jj,existing)
{
if ((*jj)==i) { add=false; break; } 
if ((*jj)==j) { add=false; break; } 

} // jj 
if (add)
{
existing.push_back(i);
existing.push_back(j);
shv.push_back(existing); 
}
} // iv
} // ii 

s=shv;
} // k 
} // longp
} // make_longer



std::map<IdxTy, IdxTy > xmerit,ymerit;
MM_LOOP(ii, shmap)
{
const IdxTy xval=(*ii).first;
//if (((*ii).second.size()>=maxma))
{
// continue; 
}
MM_LOOP(jj,(*ii).second)
{
xmerit[xval]+=(*jj).size();
MM_LOOP(kk,(*jj))
{
++ymerit[*kk];
} // kk 
 	
} // jj 

const auto sia=s.find(xval);
const auto & ma=(*sia).second;
if ((ma.size()>=maxma)) xmerit[xval]=~0; 


} // ii 

std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
const IdxTy ma=xmerit[a];
const IdxTy mb=xmerit[b];
if (ma==mb)
{
if (ma==0) return false;
const Shv & va= shmap[a];
const Shv & vb= shmap[b];
const IdxTy vasz=va.size();
const IdxTy vbsz=vb.size();
const IdxTy minsz=(vasz<vbsz)?vasz:vbsz;
for(IdxTy ivv=0; ivv<minsz; ++ivv)
{
const Sh & vva=va[ivv];
const Sh & vvb=vb[ivv];
const IdxTy szvva=vva.size();
const IdxTy szvvb=vvb.size();
const IdxTy minvvsz=(szvva<szvvb)?szvva:szvvb;
for(IdxTy jvv=0; jvv<minvvsz; ++jvv)
{
const IdxTy & av=vva[jvv];
const IdxTy & bv=vvb[jvv];
if (av!=bv) return (av>bv); 

} // jvv
if ( szvva!=szvvb) return (szvva>szvvb);

} // ivv 
if ( vasz!=vbsz) return vasz<vbsz;
return false; // ia>ib;
}
return xmerit[a]>xmerit[b]; 

});


std::sort(yord.begin(),yord.end(),
[&](const int& a, const int& b) { return ymerit[a]>ymerit[b]; });



return ; 
} // make_pairs


MM_LOOP(ii,xhisto) 
{
MM_ERR("xhisto" << MMPR2((*ii).first,(*ii).second))
}
MM_LOOP(ii,himap)
{
Ss ss;
const IdxTy cnt=(*ii).second;
const IdxTy bits=(*ii).first;
ss<<std::hex<<MMPR(bits);
MM_ERR(ss.str()<<MMPR(cnt))
}

// sort in order of frequency of the hitvector 
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ return himap[hitbits[a]]>himap[hitbits[b]]; });

std::sort(yord.begin(),yord.end(),
[&](const int& a, const int& b) { return ycnts[a]>ycnts[b]; });

} // crpobscores

void  anglescores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool fwdmo=((flags&1)==0);
const bool tie=((flags&2)==0);
//const bool lt=true;
std::map<IdxTy,D> vals;
MM_LOOP(ii,xord)
{
D peak=0;
//IdxTy loc=0;
const auto sia=s.find((*ii));
const auto & ma=(*sia).second;
// FIXME this is still slow and no idea why not just itor over ma??? 
//MM_SZ_LOOP(j,yord,yordsz)
MM_LOOP(kk,ma)
{
//const IdxTy v=yord[j];
//const auto dia=ma.find(v);
//const bool dia_end=(dia==ma.end());
D s1=0;
//if( !dia_end) 
//{ MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 
{ MM_LOOP(jj,(*kk).second) { s1+=(*jj); }  } 
peak+=s1*s1;
} // kk 
vals[(*ii)]=1.0/sqrt(peak);
} // ii

//MyBlock angles;
const IdxTy sz=xord.size();
IdxTy slog=1;
IdxTy szl=sz;
//const IdxTy compmax=sz>>4;
const IdxTy compmax=sz;
while (szl!=0) { szl>>=1; ++slog; } 
const IdxTy fakesz=(1<<slog);  
MM_ERR(" ficl lll "<<MMPR4(sz,szl,slog,fakesz))
//angles.resize(sz,sz);
std::map<D,IdxTy> dotp;
IdxTy i=0;
MM_LOOP(ii,xord)
{
const auto ma=(*(s.find((*ii)))).second;
IdxTy j=i+1;
const IdxTy ibase=i<<slog; // fakesz;
D dot=0;
IdxTy comp=0;
for (auto ip=ii+1; ip!=xord.end(); ++ip)
{
const auto mb=(*(s.find((*ip)))).second;
MM_LOOP(ll,mb)
{
if (comp>compmax) { dot=0; break; } 
auto x=ma.find((*ll).first);
if (x==ma.end()) continue;
D s1=0; D s2=0;
MM_LOOP(mm,(*ll).second) { s1+=*mm; }
MM_LOOP(mm,(*x).second) { s2+=*mm; }
dot+=s1*s2;
++comp;
}
//angles(i,j)=dot*vals[i]*vals[j]; // j>i
D score=dot*vals[i]*vals[j]; // j>i
//D score=-dot; // *vals[i]*vals[j]; // j>i
if (score!=0) {
while (dotp.find(-score)!=dotp.end()) score*=1.0+1e-9;
dotp[-score]=ibase+j; // j>i
}
//else { MM_ERR(" zero score "<<MMPR3(i,j,score)) } 
//dotp[-dot*vals[i]*vals[j]]=ibase+j; // j>i
// this is a MAP VALUE NOT and index 
//dotp[score]=ibase+j; // j>i
++j;
} // ip 
++i;
} // ii 
// find the closest hits and move them
// together somehow... 
const IdxTy mask=fakesz-1; 
std::map<IdxTy,IdxTy> clmap;
std::vector<cluster> clusters;
MM_LOOP(ii,dotp)
{
const IdxTy & idx=(*ii).second;
const D & val=(*ii).first; 
const IdxTy i=idx>>slog;
const IdxTy j=idx&mask;
// MM_ERR(MMPR4(val,idx,i,j))
const bool inew= (clmap.find(i)==clmap.end());
const bool jnew= (clmap.find(j)==clmap.end());
// these should be unique, one or other must be missing 
// actually this is ok later as these are just ones not
// of any relevance... 
if (!inew&&!jnew)
{
const IdxTy jl=clmap[j];
const IdxTy il=clmap[i];
if (jl==il) continue; 
//MM_ERR( " bridging "<<MMPR4(i,j,il,jl))
//if ( jl<il) clusters[jl].bridge(clmap[il]);
//else clusters[il].bridge(clmap[jl]);
if ( jl<il) clusters[jl].bridge(il);
else clusters[il].bridge(jl);


// 
// MM_ERR(" duplicate "<<MMPR4(i,j,clmap[i],clmap[j]))
continue;
}
if (inew&&jnew)
{
clmap[i]=clusters.size();
clmap[j]=clusters.size();
clusters.push_back(cluster(i,j));
//if (vals[i]>vals[j])clusters.push_back(cluster(i,j));
//xelse clusters.push_back(cluster(j,i));
continue;
}
// bridges 
if (inew)
{
//clmap[i]=clusters.size();
//clusters.push_back(cluster(i));
//clusters[clmap[j]].bridge(clmap[i]);

clmap[i]=clmap[j];
clusters[clmap[j]].add(i);
continue;
} // inew 
if (jnew)
{
clmap[j]=clmap[i];
clusters[clmap[i]].add(j);
//clmap[j]=clusters.size();
//clusters.push_back(cluster(j));
//clusters[clmap[i]].bridge(clmap[j]);


} // jnew

} // dotp 
IdxTy ptr=0;
IdxTy dptr=0;
std::map<IdxTy,IdxTy> visited;
XYord x2=xord;
while (true)
{
auto& cp=clusters[ptr];
MM_LOOP(ii,cp.m_set)
{
if (xord.size()<=dptr)
{
MM_ERR(" dptr too big "<<MMPR2(dptr,xord.size()))
break; 
}
//if ((*ii)!=(~0U)) 
{ x2[dptr]=xord[(*ii)]; ++dptr; } 
} // ii 
visited[ptr]=1;
//MM_ERR(" done with "<<MMPR(ptr))
++ptr;
if (ptr==clusters.size()) break; 
/*
if ( cp.m_1==(~0U)) break;
if (visited[cp.m_1]==0) { ptr=cp.m_1; continue; } 
if ( cp.m_2==(~0U)) break;
if (visited[cp.m_2]==0) { ptr=cp.m_2; continue; } 
break;
*/
}
if (xord.size()!=dptr)
{
MM_ERR(" did not come out right "<<MMPR4(xord.size(),dptr,ptr,clusters.size()))
}
xord=x2;
//std::sort(xord.begin(),xord.end(),
//[&](const int& a, const int& b) 
//{ return vals[a]>vals[b]; }); 
} // anglescores







///////////////////////////////////////



void  xxxxxtotores(const FromScores & s, const FromScores & r, XYord &xord, XYord & yord , const IdxTy flags )
{
//const bool lt=true;
std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ 
// a and b are values in xord to be compared by their maps entries for yord
//MM_ONCE(" danger will robinson kluge lt eq ",)
//if (a==b) return true; 
const auto sia=s.find(a);
const auto sib=s.find(b);
const bool sia_end=(sia==s.end());
const bool sib_end=(sib==s.end());
if (sia_end||sib_end)
{
MM_ERR(" sorting error "<<MMPR2(a,b))
}
D s1=0; D s2=0;
const auto & ma=(*sia).second;
const auto & mb=(*sib).second;
MM_SZ_LOOP(j,yord,yordsz)
{
const IdxTy v=yord[j];
const auto dia=ma.find(v);
const auto dib=mb.find(v);
const bool dia_end=(dia==ma.end());
const bool dib_end=(dib==mb.end());
// there is an entry for a not b it is earlier 
//if (!dia_end&&dib_end) return lt; // TODO kluge order? 
//if (dia_end&&!dib_end) return !lt;
//if (!dia_end&&!dib_end) 
{
if( !dia_end) { MM_LOOP(jj,(*dia).second) { s1+=(*jj); }  } 
if (!dib_end) { MM_LOOP(jj,(*dib).second) { s2+=(*jj); }  } 
//if (s1==s2) continue; // return false;
//if (s1>s2) return lt; // kluge
//return !lt;
}

}
return (s1>s2); 
// these are equal..
//return false; 

});
} // maxscores





///////////////////////////////////////////////////////////


typedef mjm_ragged_table Ragged;

void write_svg_hm(const StrTy & fn, const IdxTy flags, Ragged & fudd) 
{ std::ofstream os(fn); write_svg_hm(os,flags,fudd); }



void score_trim( FromScores & to, const FromScores & from)
{
MM_LOOP(ii,from)
{
auto sc=(*ii).second;
MM_LOOP(jj,sc)
{
D score=0;
MM_LOOP(kk,(*jj).second) {  score+=(*kk); }  
//MM_ERR(MMPR(score))
if (score==0) continue; 
to[(*ii).first][(*jj).first]=(*jj).second;
} // jj 
} // ii 
} // score_trim;

//m_scores_map[nm].organisms_paris(fn,flags);

void write_svg_hm(std::ostream & os, const IdxTy flags, Ragged & anotes)
{
D zmin=-1;
D zmax=-1;
const IdxTy _xfirst=m_flp.xfirst();
const IdxTy _yfirst=m_flp.yfirst();
const IdxTy _xlast=m_flp.xlast();
const IdxTy _ylast=m_flp.ylast();
const StrTy center_name=m_flp.center_name();
const StrTy anote_name=m_flp.anote_name();
const StrTy highlight_name=m_flp.highlight_name();
const IdxTy moment_sort_iters=m_flp.moment_sort_iters();
const IdxTy sort_iters=m_flp.sort_iters();
const StrTy run_file=m_flp.run_file();
// make a copy or dummy 
//const bool annotate=(anote_name.c_str()[0]!=0);
//Ragged  anote=annotate?m_ragged_map[anote_name]:Ragged();
std::vector<D> xv,yv;
std::vector<StrTy> cv;
std::map<IdxTy,IdxTy> lmap_from,lmap_to;
FromScores rev_scores;
const bool max_score=true;
const bool alt_loc=(center_name.c_str()[0]!=0);
ASSFUDD(m_scores,rev_scores,lmap_from,lmap_to,zmin,zmax,max_score);
zmin=log(zmin)/log(10);
zmax=log(zmax)/log(10);
MM_ERR(MMPR4(moment_sort_iters,sort_iters,zmax,zmin)<<MMPR4(_xfirst,_yfirst,_xlast,_ylast)<<MMPR2(center_name,anote_name))

FromScores r_scores_trim,f_scores_trim;
score_trim( f_scores_trim, m_scores);
score_trim( r_scores_trim, rev_scores);
m_scores=f_scores_trim;
rev_scores=r_scores_trim;
MM_ERR(" trimmed score matrix ") 
/*
MM_LOOP(ii,m_scores)
{
auto sc=(*ii).second;
MM_LOOP(jj,sc)
{
D score=0;
MM_LOOP(kk,(*jj).second) {  score+=(*kk); }  
if (score==0) continue; 
f_scores_trim[(*ii).first]=(*jj).first;
} // jj 
} // ii 
*/


//std::vector<IdxTy> xord,yord;
XYord xord,yord;
// the first index into scores is x 
MM_LOOP(ii,lmap_from) { xord.push_back((*ii).first); }
MM_LOOP(ii,lmap_to) { yord.push_back((*ii).first); }
const bool sort_size_x=((flags&1)!=0);
const bool sort_conn_x=((flags&2)!=0);
const bool sort_size_y=((flags&4)!=0);
const bool sort_conn_y=((flags&8)!=0);
const bool sort_moments=((flags&16)!=0);
const bool sort_ranks=((flags&32)!=0);
const bool sort_alpha_x=((flags&64)!=0);
const bool sort_alpha_y=((flags&128)!=0);
const bool find_blocks=((flags&256)!=0); // 0x0100
const bool sort_centroid=((flags&(1<<9))!=0); // 0x0200
const bool sort_runlen=((flags&(1<<10))!=0); // 0x0400
const bool sort_maxscores=((flags&(1<<11))!=0); // 0x0800 
const bool sort_totscores=((flags&(1<<12))!=0); // 0x01000
const bool sort_moscores=((flags&(1<<13))!=0); // 8192 0x02000
const bool sort_peakscores=((flags&(1<<14))!=0); //  0x04000
const bool sort_crossscores=((flags&(1<<15))!=0); //  0x08000
const bool dump_sort_order=((flags&(1<<16))!=0); // 10000
const bool sort_firstscores=((flags&(1<<17))!=0); //  0x20000
const bool sort_anglescores=((flags&(1<<18))!=0); //  0x40000
const bool print_scores=((flags&(1<<19))!=0); //  0x80000
const bool sort_cprobscores=((flags&(1<<20))!=0); //  0x100000
//const bool dump_sort_order=((flags&2048)!=0);
//const bool print_scores=true;
const bool bg_highlight=true;

MM_ERR(MMPR(sort_size_x) <<MMPR(sort_conn_x) <<MMPR(sort_size_y) <<MMPR(sort_conn_y) <<MMPR(sort_moments) <<MMPR(sort_ranks) <<MMPR(sort_alpha_x) <<MMPR(sort_alpha_y) <<MMPR(find_blocks) <<MMPR(sort_centroid) <<MMPR(sort_runlen)<<MMPR(sort_maxscores)<<MMPR(sort_totscores)<<MMPR(sort_moscores)<<MMPR(sort_peakscores)<<MMPR(sort_crossscores)<<MMPR(sort_firstscores)<<MMPR(sort_anglescores)<<MMPR(sort_cprobscores) <<MMPR(dump_sort_order) )


for( IdxTy dummy=0; dummy<sort_iters; ++dummy )
{
if (sort_alpha_x) { std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) { return m_st(a)<m_st(b); }); }
if (sort_alpha_y) { std::sort(yord.begin(),yord.end(),
[&](const int& a, const int& b) { return m_st(a)<m_st(b); }); }
if (sort_size_x) { std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ return Myt::sort_map_connects(m_scores,a,b); }); }

if (sort_size_y) { std::sort(yord.begin(),yord.end(),
[&](const int& a, const int& b) { 
return Myt::sort_map_connects(rev_scores,a,b); }); }


if (sort_conn_x) { std::sort(xord.begin(),xord.end(),
[&](const int& a, const int& b) 
{ return Myt::sort_map_amp(m_scores,max_score,a,b); }); }
if (sort_conn_y) { std::sort(yord.begin(),yord.end(),
[&](const int& a, const int& b) 
{ return Myt::sort_map_amp(rev_scores,max_score,a,b); }); }

if (sort_moments)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
moments(m_scores,rev_scores, xord, yord ,0) ;
moments(rev_scores,m_scores, yord, xord ,1) ;
}
}

if (sort_centroid)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
centroid(m_scores,rev_scores, xord, yord ,0) ;
centroid(rev_scores,m_scores, yord, xord ,1) ;
}
}


if (sort_runlen)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
runlen(m_scores,rev_scores, xord, yord ,0) ;
runlen(rev_scores,m_scores, yord, xord ,1) ;
}
}




if (sort_ranks)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
ranks(m_scores,rev_scores, xord, yord ,0) ;
ranks(rev_scores,m_scores, yord, xord ,0) ;
}
}

if (sort_maxscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
maxscores(m_scores,rev_scores, xord, yord ,0) ;
maxscores(rev_scores,m_scores, yord, xord ,0) ;
}
}


if (sort_totscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
totscores(m_scores,rev_scores, xord, yord ,0) ;
totscores(rev_scores,m_scores, yord, xord ,0) ;
}
}


if (sort_moscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
moscores(m_scores,rev_scores, xord, yord ,0) ;
moscores(rev_scores,m_scores, yord, xord ,1) ;
}
}

if (sort_peakscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
peakscores(m_scores,rev_scores, xord, yord ,0) ;
peakscores(rev_scores,m_scores, yord, xord ,1) ;
}
}


if (sort_crossscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
crossscores(m_scores,rev_scores, xord, yord ,0) ;
crossscores(rev_scores,m_scores, yord, xord ,1) ;
}
}

if (sort_firstscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
firstscores(m_scores,rev_scores, xord, yord ,0) ;
firstscores(rev_scores,m_scores, yord, xord ,1) ;
}
}

if (sort_anglescores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
anglescores(m_scores,rev_scores, xord, yord ,2*(i&1)) ;
anglescores(rev_scores,m_scores, yord, xord ,1+ 2*(i&1)) ;
}
}


if (sort_cprobscores)
{
for (IdxTy i=0; i<moment_sort_iters; ++i)
{
cprobscores(m_scores,rev_scores, xord, yord ,2*(i&1)) ;
}
}







} // dummy
//yord=xord;
if (dump_sort_order)
{
std::ofstream os(run_file);
MM_ERR(" writing sort order to "<<MMPR(run_file))
MM_SZ_LOOP(i,xord,xsz) { os<<MMPR3(i,xord[i],m_st(xord[i]))<<CRLF; }
MM_SZ_LOOP(i,yord,ysz) { os<<MMPR3(i,yord[i],m_st(yord[i]))<<CRLF; }

}

IdxTy __xfirst=0;
IdxTy __yfirst=0;
IdxTy __xlast=0;
IdxTy __ylast=0;
if (alt_loc)
{
const IdxTy xt=m_st(center_name);
const IdxTy yt=m_st(center_name);

MM_SZ_LOOP(i,xord,xsz) 
{if(xord[i]==xt){ if (i>10) __xfirst=i-10; else __xfirst=0;   break; }  }
__xlast=__xfirst+_xlast;
MM_SZ_LOOP(i,yord,ysz) 
{if(yord[i]==yt){if (i>10) __yfirst=i-10; else __yfirst=0;   break; }  }
__ylast=__yfirst+_ylast;


}


const IdxTy xfirst=(alt_loc)?__xfirst:_xfirst; // =m_flp.xfirst();
const IdxTy yfirst=(alt_loc)?__yfirst:_yfirst; // =m_flp.yfirst();
 IdxTy xlast=(alt_loc)?__xlast:_xlast;//=m_flp.xlast();
 IdxTy ylast=(alt_loc)?__ylast:_ylast; // =m_flp.ylast();
if (xlast>xord.size()) xlast=xord.size();
if (ylast>yord.size()) ylast=yord.size();

MM_ERR(MMPR4(xfirst,xlast,yfirst,ylast))

std::map<IdxTy,IdxTy> xfudd,yfudd;
const IdxTy xordsz=(xlast!=0)?xlast:xord.size();
const IdxTy yordsz=(ylast!=0)?ylast:yord.size();
for(IdxTy i=xfirst; i<xordsz; ++i ) {xfudd[xord[i]]=i;} 
for(IdxTy i=yfirst; i<yordsz; ++i ) {yfudd[yord[i]]=i;} 

const IdxTy xorg=xfirst;
const IdxTy yorg=yfirst;

//MM_SZ_LOOP(i,yord,yordsz) {yfudd[xord[i]]=i;} 

const IdxTy nodesx=xfudd.size(); // lmap.size();
const IdxTy nodesy=yfudd.size(); // lmap.size();

const IdxTy xsz=xfudd.size()+1; // lmap.size()+1;
const IdxTy ysz=yfudd.size()+1; // xsz;
//std::sort(indexes.begin(),indexes.end(),
//[&](const int& a, const int& b) { return (misses[a] < misses[b]); }



mjm_svg_writer sw;
IdxTy xpitch=120;
IdxTy ypitch=90;
if (xsz>100) xpitch=xpitch/10;
if (ysz>100) ypitch=ypitch/10;

 D lblsz=1.5*2.0*.5*60;
D lblszx=lblsz;
D lblszy=lblsz;
if (xsz>100) lblszy=lblszy/10;
if (ysz>100) lblszx=lblszx/10;

const IdxTy tm=2*300; // 15*xpitch; // 200;
const IdxTy lm=2*600; // 18*ypitch; // 500;
const IdxTy bm=2*800; // 41*xpitch; // 500;
const IdxTy rm=2*500; // 30*ypitch; // 300;
const IdxTy wr=xsz; // 30;
const IdxTy hr=ysz; // 30;
//mjm_output_warp ow(xpitch,ypitch,w,h);
//mjm_circular_warp ow(xpitch,ypitch,w,h);
//ow_bounds(x,y,maxlevel,maxsz);
IdxTy xs=lm+rm+xpitch*wr;
IdxTy ys=bm+tm+ypitch*hr;
//ow.bounds(x,y,lc.m_level_sz);
os<<sw.start_text(" test foo",xs,ys);
os<<sw.frame_text("#00808080",xs,ys);
os<<CRLF;
os<<sw.hstripes(lm, lm+xpitch*nodesx, tm, tm+ypitch*nodesy, nodesy, StrTy("black"),1);
os<<sw.vstripes(lm, lm+xpitch*nodesx, tm, tm+ypitch*nodesy, nodesx, StrTy("black"),1);

if (bg_highlight)
{
const IdxTy xt=m_st(highlight_name);
const IdxTy yt=m_st(highlight_name);

{
if (xfudd.find(xt)!=xfudd.end())
{
const IdxTy i=((*(xfudd.find(xt))).second);
D x=(i-xorg)*xpitch+lm;
os<<sw.shade_rect_text( x, tm,  xpitch,ypitch*hr  , "orange", .5)<<CRLF ; 
}
} // xfudd

{
if (yfudd.find(yt)!=yfudd.end())
{
const IdxTy i=((*(yfudd.find(yt))).second);
D y=(i-yorg)*ypitch+tm;
os<<sw.shade_rect_text( lm, y,  xpitch*wr,ypitch  , "orange", .5)<<CRLF ; 
}
} // yfudd
MM_ERR(" Checking annotations ")
const IdxTy asz=anotes.size();
for(IdxTy i=0; i<asz; ++i)
{
const  Ragged::Line & l=  anotes.line(i);
const IdxTy lsz=l.size();
MM_ERR(" Checking annotations "<<MMPR(lsz))
if (lsz<4) continue;
if (l[0]=="mark")
{
MM_ERR(" Checking annotations "<<MMPR2(lsz,l[1]))
if ( l[1]=="x")
{
MM_ERR(" Checking annotations "<<MMPR2(lsz,l[2]))
const auto ii=xfudd.find(m_st(l[2]));
if (ii==xfudd.end()) continue;
const IdxTy i=(*ii).second;
D x=(i-xorg)*xpitch+lm;
MM_ERR(" Checking annotations "<<MMPR3(lsz,l[3],x))
os<<sw.shade_rect_text( x, tm,  xpitch,ypitch*hr  , l[3], .5)<<CRLF ; 

} // x
else
{
const auto ii=yfudd.find(m_st(l[2]));
if (ii==yfudd.end()) continue;
const IdxTy i=(*ii).second;
D y=(i-yorg)*ypitch+tm;
os<<sw.shade_rect_text( lm, y,  xpitch*wr,ypitch  , l[3], .5)<<CRLF ; 


} // y 
continue;
} // mark

} // i 

}  // bg_highlight 


Ss ss;

for(IdxTy pass=0; pass<1; ++pass)
{

MM_LOOP(ii,m_scores)
{
const IdxTy src=(*ii).first;
// NOTE size starts at 1 maybe? Sure? 
//if (lmap.find(src)==lmap.end()) lmap[src]=lmap.size();
//const StrTy nm=m_st((*ii).first);
const auto &no=(*ii).second;
//D x=lmap[src]*xpitch+lm;
auto si= xfudd.find(src);
if (si==xfudd.end()) { continue; } 
D x=((*si).second-xorg)*xpitch+lm;
MM_LOOP(jj,no)
{
D score=0;
MM_LOOP(kk,(*jj).second) { if (max_score) {if ((*kk)>score) score=(*kk); } else score+=(*kk); } 
//MM_ERR(MMPR(score))
if (score==0) continue; 
const IdxTy dest=(*jj).first;
//if (lmap.find(dest)==lmap.end()) lmap[dest]=lmap.size();

auto sj= yfudd.find(dest);
if (sj==yfudd.end()) {
MM_ERR(" skipping "<<MMPR2(m_st(dest),m_st(src)))
 continue; } 
if (print_scores) 
{ 
	const IdxTy x= (*si).second; 
	const IdxTy y= (*sj).second; 
	MM_ERR(MMPR3(x,y,score)<<MMPR2(m_st((*si).first),m_st((*sj).first))<<" scores")

}
//D y=lmap[dest]*ypitch+tm;
//D y=yfudd[dest]*ypitch+tm;
D y=((*sj).second-yorg)*ypitch+tm;
score=log(score)/log(10);
const StrTy col=color(score,zmin,zmax);
xv.push_back(x);
yv.push_back(y);
cv.push_back(col);

} // ii 

} //jj 
} // pass
os<<sw.sparse_shade_rect_array_text(xv,yv,xpitch,ypitch, cv,1); os<<CRLF;
IdxTy in_gridsx=4;
IdxTy in_gridsy=4;
if (nodesy<150) in_gridsy=3;
if (nodesy<100) in_gridsy=2;
if (nodesy<30) in_gridsy=1;
if (nodesx<150) in_gridsx=3;
if (nodesx<100) in_gridsx=2;
if (nodesx<30) in_gridsx=1;


MM_LOOP(ii,xfudd)
{
const IdxTy src=(*ii).first;
const StrTy & namex=m_st(src);
//IdxTy ylbl=tm+ypitch*(2+lmap.size());
IdxTy ylbl=tm+ypitch*(2+nodesy);
// the rectangle is top left and the anchor is bottom lol 
IdxTy x=lm+xpitch*((*ii).second-xorg);
{const StrTy colr="red";
const D angle=60; os<<sw.gtext_text(namex,x,ylbl,lblszy,colr,StrTy("left"),angle); os<<CRLF; } 
const IdxTy ass=in_gridsy;
for (IdxTy grf=0; grf<ass; ++grf)
//{IdxTy ylbl=tm+ypitch*(2+lmap.size())*grf/ass;
{IdxTy ylbl=tm+ypitch*(nodesy)*grf/ass-ypitch;
const StrTy colr="red";
{D angle=90; if (grf==0) angle=60;  os<<sw.gtext_text(namex,x,ylbl,lblszy,colr,

(grf==0)?StrTy("end"):StrTy("left")

//StrTy("left")

,angle); os<<CRLF; } 
}
} // xord 

MM_LOOP(ii,yfudd)
{
const IdxTy src=(*ii).first;
const StrTy & namex=m_st(src);
//IdxTy xlbl=lm+xpitch*(2+lmap.size());
IdxTy xlbl=lm+xpitch*(2+nodesx);
// the rectangle is top left and the anchor is bottom lol 
IdxTy y=tm+ypitch*(1+(*ii).second-yorg);
{const StrTy colr="blue";
const D angle=30; os<<sw.gtext_text(namex,xlbl,y,lblszx,colr,StrTy("left"),angle); os<<CRLF; } 
const IdxTy ass=in_gridsx;
for (IdxTy grf=0; grf<ass; ++grf)
//{IdxTy xlbl=lm+xpitch*(2+lmap.size())*grf/ass;
{IdxTy xlbl=lm+xpitch*(nodesx-1)*grf/ass-xpitch;
const StrTy colr="blue";
{ D angle=0; if (grf==0) angle=30;  os<<sw.gtext_text(namex,xlbl,y,lblszx,colr,
(grf==0)?StrTy("end"):StrTy("left")
,angle); os<<CRLF; } 
}
} // xord 

const D yleg=.5*bm+tm+ypitch*(2+nodesy);
const D legp=xpitch*D(nodesx)/16;
const D wleg=xpitch*nodesx;
const D hleg=bm/4;
const D ssz=legp/2;
hm_leg(os, sw, zmax, zmin,lm, yleg, legp, wleg, hleg, ssz);
os<<sw.end_text();
os<<CRLF;



if (find_blocks)
{
Bm blocks;
const int nblocks=40;
const IdxTy deltamax=30;
const bool find_diag_blocks=true;
const bool ignore_diag=true;
output_corrs(os,blocks,m_scores,lmap_from,lmap_to,find_diag_blocks,!true,nblocks);
(*this).find_blocks( blocks, xord,yord, deltamax, ignore_diag);


} // find_blocks 




}
template <class Os, class Sw> 
void hm_leg(Os& os, Sw & sw,const D & zmax, const D & zmin,const D & xleg, const D & yleg, const D & legp, const D & wleg, const D & hleg, const D & ssz )
{
std::vector<D> xv,yv;
std::vector<StrTy> cv;
std::map<IdxTy,IdxTy> lmap;

for (IdxTy i=0; i<16; ++i)
{
xv.push_back(xleg+legp*i);
yv.push_back(yleg);
const D val=1.0*i/16*(zmax-zmin)+zmin ; // (1.0*i/16 -zmin)/(zmax-zmin);
StrTy col= color(val,zmin,zmax);
cv.push_back(col);
Ss ss;
ss<<std::setprecision(3)<<(val); // (val*(zmax-zmin)+zmin); 
os<<sw.gtext_text(ss.str(),xv.back()+.5*legp,yv.back()+hleg+10,ssz,cv.back(),StrTy("left"),45); os<<CRLF;

}
//os<<sw.sparse_shade_rect_array_text(xv,yv,wleg/16,hleg, cv,1); os<<CRLF;
os<<sw.sparse_shade_rect_array_text(xv,yv,wleg/16,hleg, cv,1); os<<CRLF;


//os<<sw.gtext_text(StrTy("marchywka test code"),xleg,yleg+.5*hleg,hleg,StrTy("black"),StrTy("left"),0); os<<CRLF;
os<<sw.gtext_text(StrTy("marchywka test code"),xleg,yleg+.5*hleg,ssz*4,StrTy("black"),StrTy("left"),0); os<<CRLF;






} //hm_leg





void organism_pairs(const StrTy & fn, const IdxTy flags) { std::ofstream os(fn); organism_pairs(os,flags); }
void organism_pairs(std::ostream & os, const IdxTy flags)
{
//const IdxTy xfirst=m_flp.xfirst();
//const IdxTy yfirst=m_flp.yfirst();
//const IdxTy xlast=m_flp.xlast();
//const IdxTy ylast=m_flp.ylast();
//const IdxTy moment_sort_iters=m_flp.moment_sort_iters();
//const IdxTy sort_iters=m_flp.sort_iters();
const int nblocks=30;

const bool write_os=((flags&1)==0);
const bool write_into_scores=((flags&2)!=0);
const bool find_diag_blocks=((flags&4)!=0);
//std::vector<D> xv,yv;
//std::vector<StrTy> cv;

Bm  blocks;
ToFromMap lmap_from,lmap_to;

XYord xord,yord;
if (find_diag_blocks)
{
FromScores rev_scores;
const bool max_score=true;
D zmin=-1;
D zmax=-1;
ASSFUDD(m_scores,rev_scores,lmap_from,lmap_to,zmin,zmax,max_score);
zmin=log(zmin)/log(10);
zmax=log(zmax)/log(10);

//std::vector<IdxTy> xord,yord;
// the first index into scores is x 
MM_LOOP(ii,lmap_from) { xord.push_back((*ii).first); }
MM_LOOP(ii,lmap_to) { yord.push_back((*ii).first); }
}




FromScores corr;
MM_LOOP(ii,m_scores)
{
MM_LOOP(jj,(*ii).second)
{
const IdxTy & second=((*jj).first);
MM_LOOP(ll,(*ii).second)
{
const IdxTy & first=((*ll).first);
D scorej=0;
D scorel=0;
MM_LOOP(kk,(*jj).second) {if ((*kk)>scorej) scorej=(*kk); }
MM_LOOP(kk,(*ll).second) {if ((*kk)>scorel) scorel=(*kk); }
//MM_LOOP(kk,(*ll).second) {}
auto & vvv= corr[first][second];
const D tscore=scorel*scorej;
if (vvv.size()==0) vvv.push_back(tscore);
else vvv.back()+=tscore; // scorel*scorej;
} // ll 
} // jj 

} // ii 
if (write_os)
{


output_corrs(os,blocks,corr,lmap_from,lmap_to,find_diag_blocks,true,nblocks);
#if 0

MM_LOOP(ii,corr)
{
// indexing these things by order of occurence is almost
// clustering as block diagonals appear... 
const StrTy & first=m_st((*ii).first);
const IdxTy & firstn=((*ii).first);
auto fno=lmap_from.find(firstn);
const IdxTy firsto=(*fno).second;
MM_LOOP(jj,(*ii).second)
{
const StrTy & second=m_st((*jj).first);
const D v=(*jj).second[0];
if (find_diag_blocks)
{
const IdxTy & secondn=((*jj).first);
auto tno=lmap_to.find(secondn);
const IdxTy secondo=(*tno).second;
const int offdiag=int(secondn)-int(firstn);
const int offdiago=int(secondo)-int(firsto);
//if ((offdiago<nblocks)&&(offdiago>(-nblocks))) 
if ((offdiag<nblocks)&&(offdiag>(-nblocks))) 
		blocks[firstn][secondn]=1;
		//blocks[firsto][secondo]=1;

os<<MMPR3(secondo,firsto,offdiago)<<MMPR3(firstn,secondn,offdiag);
}
os<< MMPR3(first,second,v)<<CRLF;
}
} // ii 

#endif






if (find_diag_blocks)
{
FromScores rev_scores;
const bool max_score=true;
D zmin=-1;
D zmax=-1;
//m_scores.clear(); 
lmap_from.clear(); lmap_to.clear();
ASSFUDD(corr,rev_scores,lmap_from,lmap_to,zmin,zmax,max_score);
zmin=log(zmin)/log(10);
zmax=log(zmax)/log(10);
xord.clear(); yord.clear();
//std::vector<IdxTy> xord,yord;
// the first index into scores is x 
MM_LOOP(ii,lmap_from) { xord.push_back((*ii).first); }
MM_LOOP(ii,lmap_to) { yord.push_back((*ii).first); }
}


if (find_diag_blocks)
{

const IdxTy deltamax=30;
find_blocks(blocks,xord,yord,deltamax);
#if 0 

MM_SZ_LOOP(i,xord,xsz)
{
//const IdxTy from=(*ii).first;
//const std::map<IdxTy, IdxTy> & x= (*ii).second;
IdxTy delta=0;
while (delta<deltamax)
{
bool block_ok=true;
for(IdxTy j=0; j<=delta; ++j)
{
BmCi  bm1=blocks.find(xord[i+delta]);
BmCi  bm2=blocks.find(xord[i+delta-j]);
//BmCi  bm1=blocks.find(i+delta);
//BmCi  bm2=blocks.find(i+delta-j);


if (bm1==blocks.end()) { block_ok=false; break; } 
if (bm2==blocks.end()) { block_ok=false; break; } 
DmCi dm1=(*bm1).second.find(yord[i+delta-j]);
DmCi dm2=(*bm2).second.find(yord[i+delta]);
//DmCi dm1=(*bm1).second.find(i+delta-j);
//DmCi dm2=(*bm2).second.find(i+delta);



if (dm1==(*bm1).second.end()) { block_ok=false; break; } 
if (dm2==(*bm2).second.end()) { block_ok=false; break; } 

//const bool findx=x.find(from);
}
if (!block_ok) break;
++delta;
}
if (delta>5)
{
for(IdxTy  j=0; j<delta; ++j)
{
MM_ERR(" block "<<MMPR4(i,j,m_st(xord[i+j]),m_st(yord[i+j])))
//MM_ERR(" block "<<MMPR4(i,j,m_st(i+j),m_st(i+j)))
}

}
if (delta>0) { i+=delta-1; } 
} // i 
#endif


}






} // write_os
// TODO a simple equals assinment should work here. 
if (write_into_scores)
{
m_scores=corr;
#if 0
m_scores.clear();
MM_LOOP(ii,corr)
{
//const StrTy & first=m_st((*ii).first);
const IdxTy & first=((*ii).first);
MM_LOOP(jj,(*ii).second)
{
//const StrTy & second=m_st((*jj).first);
const IdxTy & second=((*jj).first);
const D v=(*jj).second[0];
m_scores[first][second].push_back(v);
//os<<MMPR3(first,second,v)<<CRLF;
}
} // ii 
#endif

} // write_os


} // organism_paris

void output_corrs(std::ostream & os, Bm & blocks,FromScores & corr,ToFromMap & lmap_from,ToFromMap & lmap_to, const bool find_diag_blocks, const bool output_corr, const int nblocks)
{
MM_LOOP(ii,corr)
{
// indexing these things by order of occurence is almost
// clustering as block diagonals appear... 
const StrTy & first=m_st((*ii).first);
const IdxTy & firstn=((*ii).first);
auto fno=lmap_from.find(firstn);
const IdxTy firsto=(*fno).second;
MM_LOOP(jj,(*ii).second)
{
const StrTy & second=m_st((*jj).first);
const D v=(*jj).second[0];
if (find_diag_blocks)
{
const IdxTy & secondn=((*jj).first);
auto tno=lmap_to.find(secondn);
const IdxTy secondo=(*tno).second;
const int offdiag=int(secondn)-int(firstn);
const int offdiago=int(secondo)-int(firsto);
//if ((offdiago<nblocks)&&(offdiago>(-nblocks))) 
if ((offdiag<nblocks)&&(offdiag>(-nblocks))) 
		blocks[firstn][secondn]=1;
		//blocks[firsto][secondo]=1;

if (output_corr) 
{ os<<MMPR3(secondo,firsto,offdiago)<<MMPR3(firstn,secondn,offdiag); }

}
if (output_corr) { os<< MMPR3(first,second,v)<<CRLF; } 
}
} // ii 

} // output_corr

void find_blocks(const Bm & blocks,const XYord & xord,const XYord & yord,const IdxTy deltamax, const bool ignore_diag=false)
{

MM_SZ_LOOP(i,xord,xsz)
{
//const IdxTy from=(*ii).first;
//const std::map<IdxTy, IdxTy> & x= (*ii).second;
//IdxTy deltamin=ignore_diag?1:0;
const IdxTy deltamin=0;
IdxTy delta=deltamin;
while (delta<deltamax)
{
bool block_ok=true;
for(IdxTy j=deltamin; j<=delta; ++j)
{
// for many maps x and y will differ but elements
// joingin sames have been zeroed so we ignore those 
const IdxTy xod=xord[i+delta];
const IdxTy xodj=xord[i+delta-j];
const IdxTy yod=yord[i+delta];
const IdxTy yodj=yord[i+delta-j];
//BmCi  bm1=blocks.find(xord[i+delta]);
if (!ignore_diag||(xod!=yodj))
{
BmCi  bm1=blocks.find(xod);
if (bm1==blocks.end()) { block_ok=false; break; } 
DmCi dm1=(*bm1).second.find(yodj);
if (dm1==(*bm1).second.end()) { block_ok=false; break; } 
}

if (!ignore_diag||(xodj!=yod))
{
//BmCi  bm2=blocks.find(xord[i+delta-j]);
BmCi  bm2=blocks.find(xodj);
if (bm2==blocks.end()) { block_ok=false; break; } 
DmCi dm2=(*bm2).second.find(yod);
if (dm2==(*bm2).second.end()) { block_ok=false; break; } 
}
//BmCi  bm1=blocks.find(i+delta);
//BmCi  bm2=blocks.find(i+delta-j);
//DmCi dm1=(*bm1).second.find(yord[i+delta-j]);
//DmCi dm2=(*bm2).second.find(yord[i+delta]);
//DmCi dm1=(*bm1).second.find(i+delta-j);
//DmCi dm2=(*bm2).second.find(i+delta);
//const bool findx=x.find(from);
} // for j 
if (!block_ok) break;
++delta;
}
if (delta>5)
{
for(IdxTy  j=0; j<delta; ++j)
{
MM_ERR(" block "<<MMPR4(i,j,m_st(xord[i+j]),m_st(yord[i+j])))
//MM_ERR(" block "<<MMPR4(i,j,m_st(i+j),m_st(i+j)))
}
}
if (delta>0) { i+=delta-1; } 
} // i 
//}

} // find_blocks

#if 0

IdxTy dummy( const StrTy & s) const
{
IdxTy zed=0;
const char * c=s.c_str(); 
const IdxTy sz=s.length();
for (IdxTy i=0; i<sz; ++i) {if (c[i]=='-') ++zed; }  

return zed;
}
StrTy stats()
{
Ss ss;
IdxTy i=0;
MM_LOOP(ii,m_map)
{
const StrTy & name=(*ii).first;
const StrTy & seq=(*ii).second;
//ss<<MMPR(i)<<MMPR4(name,seq,seq.length(),dummy(seq))<<CRLF;
ss<<MMPR(i)<<MMPR3(name,seq.length(),dummy(seq))<<CRLF;
++i;
}

return ss.str();
}

#endif


#if 0

void  compare(const StrTy & nm) const
{
auto is=m_map.find(nm);
if (is==m_map.end())
{
MM_ERR(" exact match not found "<<nm) 
MM_LOOP(ii,m_map)
{
if (strncmp((*ii).first.c_str(),nm.c_str(),nm.length())==0)
{
MM_ERR(" using "<<nm) 
is=ii;
break;
}
}
if (is==m_map.end()) return;
}
std::vector<IdxTy> lens;
MM_LOOP(ii,m_map){ lens.push_back((*ii).second.length()); }
const StrTy & unk=(*is).second;
const IdxTy len=unk.length();
const IdxTy lenreal=len-dummy(unk);
const char * s2=unk.c_str();

std::vector<IdxTy> misses,indexes;

std::vector<StrTy> names;
IdxTy j=0;
MM_LOOP(ii,m_map)
{
const char * s1=(*ii).second.c_str();
const IdxTy sz=(len<lens[j])?len:lens[j];
IdxTy miss=0;
for (IdxTy i=0; i<sz; ++i)
{
if ( s2[i]!='-') if ( s2[i]!=s1[i]) ++miss;

} // i 
misses.push_back(miss);
indexes.push_back(j);
names.push_back((*ii).first);
//const D pct=100.0*(lenreal-miss)/lenreal;
//MM_MSG(MMPR4(j,(*ii).first, (*is).first,miss)<<MMPR(pct))
++j;
} // ii 
j=0;
std::sort(indexes.begin(),indexes.end(),
[&](const int& a, const int& b) { return (misses[a] < misses[b]); }

);
for (IdxTy i=0; i<indexes.size(); ++i)
{
const IdxTy idx=indexes[i];
IdxTy miss=misses[idx];
const D pct=100.0*(lenreal-miss)/lenreal;
MM_MSG(MMPR4(j,names[idx], (*is).first,miss)<<MMPR(pct))
++j;

}
} //compare
#endif
void fudd_params(const ParamGlob & fudd)
{ m_flp=fudd; } 
private:
void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
//li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,';');
li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,' ');
}

ParamGlob  m_flp;
//SeqMap m_map;
FromScores m_scores;
//Ordering m_order;
St m_st;
bool m_read_only;
IdxTy m_min_words,m_max_words,m_fields;

}; // mjm_pro_pheno




class mjm_string_tax
{
typedef string_tax_traits::Tr  Tr;
typedef string_tax_traits::util  Util;
typedef mjm_string_tax Myt;
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


typedef string_tax_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

typedef mjm_string_graph Aln;
typedef std::map<StrTy, Aln> AlnMap;

//typedef mjm_align_collection Agc;
//typedef std::map<StrTy, Agc> AgcMap;


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
mjm_string_tax():m_dmel(new Dmel()) {Init();}
mjm_string_tax(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_string_tax()
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



void cmd_read_fasta(Cip & cip , LocalVar & lv )
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_fasta_map[name].load(fn);
MM_ERR("reading fasta fn into name "<<MMPR3(fn,m_fasta_map[name].size(),name))
}
void cmd_write_fasta(Cip & cip , LocalVar & lv )
{
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
MM_ERR("writing fasta fn into name "<<MMPR3(fn,m_fasta_map[name].size(),name))
Fasta & fasta=m_fasta_map[name];
std::ofstream os(fn);
MM_SZ_LOOP(i,fasta,fsz)
{
os<<fasta.name(i)<<CRLF;
os<<fasta.seq(i)<<CRLF;


}

}



///////////////////////////////////////////////////////////////////////




void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
{
const StrTy name=cip.p1;
const StrTy seq=cip.p2;
Fasta & f=m_fasta_map[name];
StrTy sn=cip.wif(3);
if (sn.length()==0)
{ Ss ss ; ss<<"seq"<<(f.size()+1); sn=ss.str(); }
sn=StrTy(">") + sn;
f.add(sn,seq); // load(fn);
MM_ERR("adding to  "<<MMPR4(name,m_fasta_map[name].size(),sn,seq))
}

void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
{
//const StrTy name=cip.p1;
const StrTy dfasta=cip.p1;
// file with sequence numbers only 
const StrTy sfasta=cip.p2; // wif(3);
// annotation in zymo format
const StrTy srag=cip.wif(3);
const StrTy qiime=cip.wif(4);
//CommandInterpretter li;
Fasta & df=m_fasta_map[dfasta];
const Fasta & sf=m_fasta_map[sfasta];
const Ragged & sr=m_ragged_map[srag];
const Ragged & qiiq=m_ragged_map[qiime];
// xlate_map xlate_field_map(const IdxTy key, const IdxTy val)
Ragged::xlate_double_map qii_map=qiiq.xlate_field_max_map(8,3);
MM_ERR(" SHTFUDD "<<MMPR(qii_map.size()))
std::map<StrTy, StrTy> kv;
const IdxTy szsr=sr.size();
for (IdxTy i=0; i<szsr; ++i)
{
//MM_ERR(" ASSFUDD "<<MMPR2(i,sr.line(i).size()))
const IdxTy fudd=sr.line(i).size();
if (fudd!=1)
{
MM_ERR(" the ragged file needs to be loaded with flag bit 0 set (1) "<<MMPR3(i,fudd,srag))
return;
}
const char * c=sr.line(i)[0].c_str();
StrTy k="";
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; k=StrTy(d);
while (c[j]!=0) { if ((c[j]!='\t') &&(c[j]!=' ')) break; ++j; }
kv[k]=StrTy(c+j);

}

MM_ERR("zumo merge fasta fn into name "<<MMPR3(dfasta,sfasta,srag)<<MMPR3(df.size(),sf.size(),sr.size()))
const IdxTy sz=sf.size();
StrTy sep=" ";
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & seq=sf.seq(i);
StrTy sn=sf.name(i);
StrTy sqname="";
const char * c=sn.c_str();
if (c[0]!=0) { 
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; 
StrTy k=StrTy(d+1);

StrTy annot=kv[k];
const char * shotfudd=sn.c_str();
sqname=StrTy(shotfudd+1); // sn; // k;
//MM_ERR(MMPR2(k,annot))
StrTy sfx=""; 
auto im=qii_map.find(sqname);
//MM_ERR(" looking up "<<MMPR(sqname))
if (im!=qii_map.end()) { Ss sfudd; sfudd<<" qzscore="; sfudd<<(*im).second; sfx=sfudd.str(); } 
else { 
// MM_ERR(" ASSFUCL "<<MMPR2(sqname,(*(qii_map.begin())).first))

 } 
sn=sn+sep+zymo_condense(annot)+sfx;
}
df.add(sn,seq); // load(fn);
}



}


StrTy zymo_condense(const StrTy & s ) const
{
const char * c=s.c_str();
const IdxTy sz=s.length();
//IdxTy start=0;
StrTy genus="";
StrTy species="";
if (sz<3) return StrTy();
for (IdxTy i=0; i<(sz-3); ++i)
{
if ((c[i]=='g')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
// fudding buserror? wtf memcpy fudded up? 
// the bus error got me, I think the length j is just wrong doh... 
char d[j+1-i]; //memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
genus=StrTy(d+3);
continue;
}
if ((c[i]=='s')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
// memcpy probably ok just had wrong length.. 
char d[j+1-i]; // memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
species=StrTy(d+3);
continue;
}



} // i 
return genus+StrTy(" ")+species;
} 


void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
if (tabsep) m_ragged_map[name].sep("\t");
if (load_parsed) m_ragged_map[name].load(fn);
else m_ragged_map[name].load_lines(fn);

MM_ERR(MMPR2(m_ragged_map[name].size(),name))
}

//void load_genus_strings_fasta(const StrTy  & fn )
void cmd_load_genus_strings_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_scores_map[name].load_genus_strings_fasta(fn);
MM_ERR(" loading genus strings "<<MMPR3(fn,m_scores_map[name].size(),name))
}


void cmd_load_genus_scores(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_scores_map[name].load(fn,flags);
MM_ERR(MMPR4(fn,m_scores_map[name].size(),name,flags))
}
void cmd_cmp_aln(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy sn=cip.p2;
//m_scores_map[name].compare(sn);
MM_ERR(MMPR3(sn,m_scores_map[name].size(),name))
}





void cmd_query_aln(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy query=cip.p2;
const StrTy field=cip.wif(3);
const Aln & pro=m_scores_map[name];
MM_ERR(" query_protrait "<<MMPR3(name,query,field))
//auto w=pro.query(query,field);
//MM_SZ_LOOP(i,w,sz) {MM_MSG(MMPR4(name, query, field, w[i])<<MMPR(i)) } 
//MM_ERR(MMPR2(m_protrait_map[name].size(),name))

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


void cmd_organism_pairs(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy nm=cip.p2;
//const StrTy field=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(" organism pairs   "<<MMPR3(fn,nm,flags))
m_scores_map[nm].fudd_params(m_flp);
m_scores_map[nm].organism_pairs(fn,flags);
MM_ERR(" done organism paris  "<<MMPR3(fn,nm,flags))
}


StrTy help_write_svg_hm()
{
Ss ss;
ss<<" write_svg_hm fn source flags"<<CRLF; 
ss<<" sort_size_x=((flags&1)!=0);"<<CRLF; 
ss<<" sort_conn_x=((flags&2)!=0);"<<CRLF; 
ss<<" sort_size_y=((flags&4)!=0);"<<CRLF; 
ss<<" sort_conn_y=((flags&8)!=0);"<<CRLF; 
ss<<" sort_moments=((flags&16)!=0);"<<CRLF; 
ss<<" sort_ranks=((flags&32)!=0);"<<CRLF; 
ss<<" sort_alpha_x=((flags&64)!=0);"<<CRLF; 
ss<<" sort_alpha_y=((flags&128)!=0);"<<CRLF; 
ss<<" find_blocks=((flags&256)!=0);"<<CRLF; 
ss<<" sort_centroid=((flags&512)!=0);"<<CRLF; 
ss<<" sort_runlen=((flags&1024)!=0);"<<CRLF; 
ss<<" dump_sort_order=((flags&2048)"<<CRLF;

return ss.str();
}

void cmd_write_svg_hm(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy nm=cip.p2;
//const StrTy field=cip.wif(3);
const StrTy ragged=cip.wif(4);
const IdxTy flags=myatoi(cip.wif(3));
if (fn=="-help")
{
MM_MSG(help_write_svg_hm())
return;
}

MM_ERR(" wrting svg hm  "<<MMPR3(fn,nm,flags))
m_scores_map[nm].fudd_params(m_flp);
m_scores_map[nm].write_svg_hm(fn,flags,m_ragged_map[ragged]);
MM_ERR(" done wrting svg hm "<<MMPR3(fn,nm,flags))
}






void cmd_write_svg(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy nm=cip.p2;
const StrTy field=cip.wif(3);
MM_ERR(" wrting svg "<<MMPR2(fn,nm))
m_scores_map[nm].write_svg(fn);
MM_ERR(" done wrting svg "<<MMPR2(fn,nm))
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
MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_scores_map) { MM_MSG("m_scores_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_ONCE(" something wrong here need to check svn",)
//MM_LOOP(ii,m_aligns_map) { MM_MSG("m_aligns_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//AgcMap m_aligns_map;
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();

}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("load-genus-scores")]=&Myt::cmd_load_genus_scores;

m_cmd_map[StrTy("load-genus-strings-fasta")]=&Myt::cmd_load_genus_strings_fasta;
//void load_genus_strings_fasta(const StrTy  & fn )

m_cmd_map[StrTy("compare-aln")]=&Myt::cmd_cmp_aln;
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
m_cmd_map[StrTy("write-svg")]=&Myt::cmd_write_svg;
m_cmd_map[StrTy("write-svg-hm")]=&Myt::cmd_write_svg_hm;
m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
m_cmd_map[StrTy("organism-pairs")]=&Myt::cmd_organism_pairs;
//m_cmd_map[StrTy("string-prob")]=&Myt::cmd_string_prob;


/*
m_cmd_map[StrTy("load-pheno")]=&Myt::cmd_load_pheno;
m_cmd_map[StrTy("load-hbpheno")]=&Myt::cmd_load_hbpheno;
m_cmd_map[StrTy("write-hbpheno")]=&Myt::cmd_write_hbpheno;
m_cmd_map[StrTy("merge-hbpheno")]=&Myt::cmd_merge_hbpheno;
*/
m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
//m_cmd_map[StrTy("query-pheno")]=&Myt::cmd_query_pheno;
//m_cmd_map[StrTy("query-hbpheno")]=&Myt::cmd_query_hbpheno;
//m_cmd_map[StrTy("enumerate-protrait")]=&Myt::cmd_enumerate_protrait;
//m_cmd_map[StrTy("protrait-words")]=&Myt::cmd_protrait_words;
//m_cmd_map[StrTy("annotate")]=&Myt::cmd_annotate;
//m_cmd_map[StrTy("pheno-annotate")]=&Myt::cmd_pheno_annotate;
//m_cmd_map[StrTy("hbpheno-annotate")]=&Myt::cmd_hbpheno_annotate;
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
os<<ss;

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

#if 0 

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

#endif

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
FastaMap m_fasta_map;
RaggedMap m_ragged_map;
AlnMap m_scores_map;
//AgcMap m_aligns_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
CharMat m_char_mat;
//MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;

}; //mjm_trees_and_tables



/////////////////////////////////////////////////////////

#ifdef  TEST_STRING_TAX__
int main(int argc,char **args)
{
typedef mjm_string_tax  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

