#ifndef MJM_RENDERING_TEMP_DATA_H__
#define MJM_RENDERING_TEMP_DATA__H__
 
#include "mjm_globals.h"
#include "mjm_logic_base.h"
#include "mjm_ordering.h"
#include "mjm_string_tokenizer.h"
#include "mjm_taxon_tools.h"
/*

#include "mjm_latex_writer.h"
// for the presence absence vector 
#include "mjm_char_mat.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include "mjm_cli_ui.h"
#include "../mjm_fasta_ii.h"


#include "mjm_collections.h"
#include "mjm_svg_writer.h"
*/

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

////////////////////////////////////////////////////////////////

class rendering_temp_data_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
rendering_temp_data_params( const StrTy & nm) : Super(nm) {}
rendering_temp_data_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
/*
IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
IdxTy maxdepth() const { return m_map.get_uint("maxdepth",3); } // // 100;
bool print_counts() const { return m_map.get_bool("print_counts",!true); }
bool print_haves() const { return m_map.get_bool("print_haves",!true); }
bool print_havenots() const { return m_map.get_bool("print_havenots",!true); }
bool print_if_have() const { return m_map.get_bool("print_if_have",!true); }
bool suppress_vector() const { return m_map.get_bool("suppress_vector",!true); }
bool add_level() const { return m_map.get_bool("add_level",true); }
bool print_hit() const { return m_map.get_bool("print_hit",true); }
IdxTy human_digits() const { return m_map.get_uint("human_digits",1); }
IdxTy data_digits() const { return m_map.get_uint("data_digits",12); }
StrTy human_sep() const { return m_map.get_string("human_sep","\t"); }
StrTy data_sep() const { return m_map.get_string("data_sep"," "); }
*/
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
/*
ss<<"accmode="<<accmode()<<sep;
ss<<"maxdepth="<<maxdepth()<<sep;
ss<<"print_counts="<<print_counts()<<sep;
ss<<"print_haves="<<print_haves()<<sep;
ss<<"print_havenots="<<print_havenots()<<sep;
ss<<"print_if_have="<<print_if_have()<<sep;
ss<<"suppress_vector="<<suppress_vector()<<sep;
ss<<"add_level="<<add_level()<<sep;
ss<<"print_hit="<<print_hit()<<sep;
ss<<"human_digits="<<human_digits()<<sep;
ss<<"data_digits="<<data_digits()<<sep;
ss<<"human_sep="<<human_sep()<<sep;
ss<<"data_sep="<<data_sep()<<sep;
*/

//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace rendering_temp_data /// _traits
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
//typedef  data_model_error_log Dmel; 
//typedef unsigned int KeyCode;
//typedef unsigned __int128 KeyCode;
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
typedef std::vector<IdxTy > Locations;



//typedef mjm_sparse_matrix<D> MySparse;
}; // 

class point { public : 
typedef  ::rendering_temp_data::Tr  Tr;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
point(const D & xx, const D & yy): x(xx),y(yy) {}
point(): x(0),y(0) {}
 operator StrTy()
{ Ss ss;
ss<<MMPR2(x,y); 
return ss.str();
}
  D x,y;};


class name_splitter  {
typedef  ::rendering_temp_data::Tr  Tr;
typedef name_splitter Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef StrTy Sym;
public:

typedef std::vector<StrTy> input_type;
typedef std::vector<input_type> splitted_type;

void split( splitted_type & st, const input_type & it, const IdxTy flags=0) const
{

MM_LOOP(ii,it)
{
input_type x;
// gradually push pieces of *ii onto x.
IdxTy p=0;
IdxTy q=0;
const char * _s=(*ii).c_str();
const IdxTy len=(*ii).length();
char s[len+1];
memcpy(s,_s,len+1);
//s[len]=0;
while ( true)
{
char&  sq=s[q];
if (sq==0) break; 
if (sq==':')
{
sq=0;
// push back anyway to let user sort it out? Null scoping strings? 
if ( p<q) x.push_back(StrTy(s+p));
p=q+1;

}
else
{ }
++q;
} 
if (s[p]!=0)  x.push_back(StrTy(s+p)); 
st.push_back(x);
} // it

}


private:


}; // name_splitter


class legend_layout  {
typedef  ::rendering_temp_data::Tr  Tr;
typedef legend_layout Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef StrTy Sym;
typedef StrTy Label;
class location
{
public:
void rect( D & x, D & y, D & w , D & h, StrTy &color) const
{
x=m_x; y=m_y; w=m_w; h=m_h;
color=m_color;

}
void text(D& x, D&y, StrTy & text, StrTy & size, D & angle) const 
{
x=m_x; y=m_y; // w=m_w; h=m_h;
text=m_text; size=m_size; angle=m_angle;
}
void symbol(const StrTy & s) { m_color=s; } 
void label(const StrTy & s) { m_text=s; } 
void locate(const D & x,  const D& y, const StrTy & size,const D & angle)
{ m_x=x; m_y=y; m_size=size;m_angle=angle; }
void locate( const D & x, const D & y, const D & w , const D & h) // const
{
m_x=x; m_y=y; m_w=w; m_h=h;
//m_color=color;

}

//private:
D m_x,m_y,m_w,m_h;
StrTy m_color;
StrTy m_font;
StrTy m_text;
StrTy m_size;
D m_angle;

};

class sort_alpha
{
public:
sort_alpha(Myt & m) : myt(m) {}
bool operator()(const IdxTy a, const IdxTy b) const
{ return myt.m_label[a]<myt.m_label[b];} 
Myt & myt;
}; // alpha
typedef std::vector<Sym> SymbolArray ;
typedef std::vector<Label> LabelArray ;
typedef std::vector<location> SymLocations;
//typeded std::vector<location> SymLocations;

public:
legend_layout(const SymbolArray & s, const LabelArray & l)
: m_sym(s), m_label(l) { Init(); }  

void label(D & x, D & y, StrTy & text, StrTy & size, D & angle, const IdxTy i) const
{
const location & fudd=m_label_placed[m_order[i]];
fudd.text(x,y,text, size, angle);
}

void symbol( D & x, D & y, D & w , D & h, StrTy &color, const IdxTy i ) const
{

const location & z=m_sym_placed[m_order[i]];
z.rect( x, y, w , h, color); 

}
void layout( const D & x, const D & y, const D & w, const D &h)
{
IdxTy n=m_label.size();
IdxTy chars=0;
IdxTy maxchar=0;
MM_LOOP(ii, m_label) 
{
const IdxTy sz=(*ii).length();
if (sz>maxchar) maxchar=sz;
chars+=sz;

}
const D amax=5.0*maxchar*n;
const D area=w*h;
const D scale=sqrt(area/amax);

const IdxTy dx=(maxchar+2)*scale;
const IdxTy dy=5*scale;
const IdxTy cols=w/dx;
const IdxTy rows=h/dy;
MM_SZ_LOOP(i,m_sym,sz) 
{
IdxTy j=m_order[i];
D xi=x+(i%cols)*dx;
D yi=y+(i/cols)*dy;
D wi=2*scale;
D hi=2*scale;
D sizei=scale;
Ss ss; ss<<sizei;
m_sym_placed[j].locate(xi,yi,wi,hi);
// just put the text over the block 
// this alsmost works except for white on white 
//xi+=0*wi; yi+=.5*hi;
xi+=wi; yi+=.5*hi;
m_label_placed[j].locate(xi,yi,ss.str(),0);
} 


}

private:
void Init()
{
MM_LOOP(ii, m_sym) 
	{m_sym_placed.push_back(location()); m_sym_placed.back().symbol(*ii); }
MM_LOOP(ii, m_label) 
	{m_label_placed.push_back(location()); m_label_placed.back().label(*ii); }
MM_SZ_LOOP(i,m_sym,sz) {m_order.push_back(i); } 
MM_MSG(MMPR4(m_sym_placed.size(),m_label_placed.size(),m_sym.size(),m_label.size()))
MM_LOOP(ii,m_sym) {MM_MSG("m_sym" << MMPR(*ii)) }
MM_LOOP(ii,m_label) {MM_MSG("m_label" << MMPR(*ii)) }

}
friend sort_alpha;

SymbolArray m_sym;
LabelArray m_label;
SymLocations m_sym_placed;
SymLocations m_label_placed;
std::vector<IdxTy> m_order;
}; // legend_layout



class data_node  {
typedef  ::rendering_temp_data::Tr  Tr;
typedef data_node Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::MyBlock  MyBlock;

//typedef mjm_taxon_tools Mtt;
//typedef Mtt::ignore_t ignore_t;
typedef std::vector<IdxTy> OriginHier;
typedef std::vector<IdxTy> DisplayHier;
typedef std::vector<IdxTy> Flags;
typedef std::vector<IdxTy> Notes;
//typedef std::vector<D> Values;
typedef std::map<IdxTy, D> Values;
public:
data_node( const OriginHier & h, const Values & v) : m_hier(h),m_names(h),m_values(v),
m_flagged(false), m_src(0)  {
//MM_LOOP(ii,m_values) { MM_ERR(" fuddi g"<<MMPR2((*ii).first,(*ii).second))}
}
data_node( const OriginHier & h, const Values & v, const Notes & notes) 
	: m_hier(h),m_names(h),m_values(v),m_notes(notes) , m_flagged(false) , m_src(0)  {
//MM_LOOP(ii,m_values) { MM_ERR(" fuddi g"<<MMPR2((*ii).first,(*ii).second))}
}

const IdxTy & src() const { return m_src; } 
const IdxTy & src(const IdxTy src)  { m_src=src; return m_src; } 

const bool flagged() const { return m_flagged; } 
//const IdxTy flag() const { return m_flag; }  
const IdxTy nflags( ) const { return m_flags.size(); }  
const IdxTy get_flag(const IdxTy i ) const { return m_flags[i]; }  
//const IdxTy flag(const IdxTy f )  { m_flag=f; m_flagged=true;  return m_flag; }  
const Myt &  set_flag(const IdxTy f )  { m_flags.push_back(f) ; m_flagged=true;  return *this; }  
const Values & values() const { return m_values; } 
const D value(const IdxTy & sample) const
{
const auto ii=m_values.find(sample);
if (ii==m_values.end()) return 0;
return (*ii).second;
}
const D value() const
{
D sum=0;
//MM_LOOP(ii,m_values) { sum+=(*ii).second; } 
MM_LOOP(ii,m_values) { const D x=(*ii).second; if (x!=0) sum+=log(x); 
else MM_ERR(" danger will robinson a ero value "<<MMPR((*ii).first))  } 
return sum;
}

const Notes & notes() const { return m_notes; } 
typedef DisplayHier lineage;
const lineage &  names() const { return m_names; } 
Myt & operator+=(const Myt & that)
{
MM_ONCE(" combining values but trust user on equality lol",)
combine_values(that.m_values);
return *this;
}
void combine_values(const Values & x) 
{
MM_LOOP(ii,x) { m_values[(*ii).first]+=(*ii).second; } 

}
void value_lineage()
{
if (m_values.size()>0)
{
if (m_hier.size()>0) m_hier[0]=m_values[0];
else m_hier.push_back(m_values[0]);
if (m_names.size()>0) m_names[0]=m_values[0];
else m_names.push_back(m_values[0]);



}

}

private:

OriginHier m_hier;
DisplayHier m_names;
Flags m_flags;
Values m_values;
Notes m_notes;
bool m_flagged;
//std::vector<IdxTy> m_flags;
//IdxTy m_flag;
IdxTy m_src;

}; // data_node


class data_node_collection  {
typedef  ::rendering_temp_data::Tr  Tr;
typedef data_node_collection Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;

typedef string_tokenizer St;
typedef mjm_uint_ordering  Or;
typedef data_node Node;
typedef std::vector<Node>  Nodes;
typedef std::vector<StrTy>  Notes;
//typedef Tr::MyBlock  MyBlock;

typedef std::vector<StrTy >  Lineage;


//typedef std::vector<D >  Values;
// this is a fcuking map 
typedef std::map<StrTy,D >  Values;


//typedef mjm_taxon_tools Mtt;
//typedef Mtt::ignore_t ignore_t;

//class sort_alpha_lineage_comp
class sort_base_lineage_comp
{
public:
//sort_alpha_lineage_comp( const St & st, const Nodes & n) : m_st(st), m_nodes(n) {}
sort_base_lineage_comp( const St & st, const Nodes & n) : m_st(st), m_nodes(n) {}
// note that LINEAGE starts with the node in question NOT the kingdom 

/*
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{
const auto & x=m_nodes[_x].names();
const auto & y=m_nodes[_y].names();
return vector_lt(x,y);
}
*/
template < class Tv> 
bool vector_lt(const Tv & x, const Tv & y) const
{
IdxTy i=0; 
IdxTy szx=x.size(); 
IdxTy szy=y.size(); 
while (true)
{
bool s1= (i<szx);
bool s2= (i<szy);
// let optimizer fix this lol 
if ( s1&&!s2) return false;
if ( (!s1)&&s2) return !false;
if ( (!s1)&&(!s2)) return false;
const auto &  vx=x[szx-1-i];
const auto &  vy=y[szy-1-i];
if (vx!=vy)  return m_st(vx)<m_st(vy); 
++i;
}
// unreachable 
return false; 
}
protected:

const St & m_st;
const Nodes & m_nodes;

};


class sort_alpha_lineage_comp : public  sort_base_lineage_comp
{
typedef sort_base_lineage_comp Super;
public:
sort_alpha_lineage_comp( const St & st, const Nodes & n) : Super(st,n) {}
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{
const auto & x=m_nodes[_x].names();
const auto & y=m_nodes[_y].names();
return vector_lt(x,y);
}

}; // sort_slpha_lineage_comp



template <class NodeValues>
class sort_values_fick_comp : public  sort_base_lineage_comp
{
typedef sort_base_lineage_comp Super;
public:
//sort_values_comp( const St & st, const Nodes & n, const Values & v) : Super(st,n),m_v(v) {}
//sort_values_comp( const St & st, const Nodes & n) : Super(st,n) {}
sort_values_fick_comp( const St & st, const Nodes & n, const StrTy &nm, const NodeValues & nv) 
: Super(st,n),m_sn(nm),m_nv(nv) {}
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{

const auto & x=m_nv.cvalue(_x,m_sn); // m_nv[_x][m_sn]; // m_nodes[_x].values();
const auto & y=m_nv.cvalue(_y,m_sn); // m_nodes[_y].values();
MM_ERR(" value sort "<<MMPR4(_x,_y,m_sn,x)<<MMPR(y))
return x<y;
//const auto fudd1=x.find(m_st[m_sn]); // x.begin(); // find(_x);
//const auto fudd2=y.find(m_st[m_sn]); // y.begin(); // find(_x);
//return (*fudd1).second<(*fudd2).second; 
//D a=(fudd1!=x.end())? (*fudd1).second:0; // <(*fudd2).second; 
//D b=(fudd2!=y.end())? (*fudd2).second:0; // <(*fudd2).second; 
//return (*fudd1).second<(*fudd2).second; 

//return a<b;
//const auto & x=m_nodes[_x].names();
//const auto & y=m_nodes[_y].names();
//const IdxTy s=0;
//return m_v[s][_x]<m_v[s][_y];
//const auto fudd2=m_v.find(_y);
//if (fudd1==m_v.end()) { MM_ERR(" fudd1 missing "<<MMPR2(_x,_y))}
//if (fudd2==m_v.end()) { MM_ERR(" fudd2 missing "<<MMPR2(_x,_y))}

//return (*fudd1).second<(*fudd2).second; 
// this is not const
//return m_v[_x]<m_v[_y] ;
//return vector_lt(x,y);
}
//const Values & m_v;
const StrTy m_sn;
const NodeValues & m_nv;
}; // sort_slpha_lineage_comp

/////////////////////////////////////////////

class sort_values_comp : public  sort_base_lineage_comp
{
typedef sort_base_lineage_comp Super;
public:
//sort_values_comp( const St & st, const Nodes & n, const Values & v) : Super(st,n),m_v(v) {}
//sort_values_comp( const St & st, const Nodes & n) : Super(st,n) {}
sort_values_comp( const St & st, const Nodes & n, const StrTy &nm, const StrTy &mods) 
: Super(st,n),m_sn(nm),m_mods(mods) {}
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{
if(m_sn.length()!=0)
{
// at this point, the nodes have been collected by dnc but nothing else is
// done. The string tokenizer should be valid 
const auto & x= m_nodes[_x].value(m_st[m_sn]);
const auto & y= m_nodes[_y].value(m_st[m_sn]);

if ((y+x)>0)return x<y;
const auto & xv=m_nodes[_x].values();
const auto & yv=m_nodes[_y].values();
MM_LOOP(ii,xv)
{
const auto & jj=yv.find((*ii).first);
if (jj!=yv.end()) return ((*ii).second<(*jj).second); }

}
const auto & x= m_nodes[_x].value();
const auto & y= m_nodes[_y].value();

//dncvalue(_x,m_sn); // m_nv.cvalue(_x,m_sn); // m_nv[_x][m_sn]; // m_nodes[_x].values();
//const auto & y=dncvalue(_y,m_sn); // m_nv.cvalue(_y,m_sn); // m_nodes[_y].values();
//MM_ERR(" value sort "<<MMPR4(_x,_y,m_sn,x)<<MMPR(y))
return x<y;
//const auto fudd1=x.find(m_st[m_sn]); // x.begin(); // find(_x);
//const auto fudd2=y.find(m_st[m_sn]); // y.begin(); // find(_x);
//return (*fudd1).second<(*fudd2).second; 
//D a=(fudd1!=x.end())? (*fudd1).second:0; // <(*fudd2).second; 
//D b=(fudd2!=y.end())? (*fudd2).second:0; // <(*fudd2).second; 
//return (*fudd1).second<(*fudd2).second; 

//return a<b;
//const auto & x=m_nodes[_x].names();
//const auto & y=m_nodes[_y].names();
//const IdxTy s=0;
//return m_v[s][_x]<m_v[s][_y];
//const auto fudd2=m_v.find(_y);
//if (fudd1==m_v.end()) { MM_ERR(" fudd1 missing "<<MMPR2(_x,_y))}
//if (fudd2==m_v.end()) { MM_ERR(" fudd2 missing "<<MMPR2(_x,_y))}

//return (*fudd1).second<(*fudd2).second; 
// this is not const
//return m_v[_x]<m_v[_y] ;
//return vector_lt(x,y);
}
//const Values & m_v;
const StrTy m_sn,m_mods;
//const NodeValues & m_nv;
}; // sort_slpha_lineage_comp

/////////////////////////////////////////////




////////////////////////////////////////////


class sort_alpha_notes_comp : public  sort_base_lineage_comp
{
typedef sort_base_lineage_comp Super;
public:
sort_alpha_notes_comp( const St & st, const Nodes & n) : Super(st,n) {}
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{
const auto & x=m_nodes[_x].notes();
const auto & y=m_nodes[_y].notes();
return vector_lt(x,y);
}
bool equals(const Node & x, const Node & y) const
{
return x.values()==y.values(); 
}

}; // sort_slpha_lineage_comp

class sort_alpha_gs_comp : public  sort_base_lineage_comp
{
typedef sort_base_lineage_comp Super;
public:
sort_alpha_gs_comp( const St & st, const Nodes & n, mjm_taxon_tools & tt) : Super(st,n), m_tt(tt)  {}
bool operator()( const IdxTy & _x, const IdxTy & _y) const
{
const auto & x=m_nodes[_x].names();
const auto & y=m_nodes[_y].names();
std::vector<StrTy>  inx= m_tt.informative_vector(m_st(x));
std::vector<StrTy>  iny= m_tt.informative_vector(m_st(y));
while (inx.size()<2) { inx.push_back("MISSING"); }
while (iny.size()<2) { iny.push_back("MISSING"); }
if (inx[1]>iny[1]) return false;
if (inx[1]<iny[1]) return !false;
return inx[0]<iny[0];
//return vector_lt(x,y);
}

private:
//St & m_st;
mjm_taxon_tools& m_tt;

}; // sort_slpha_lineage_comp






class node_iterator 
{
typedef node_iterator Myt;
public:
node_iterator(St & s, Nodes & n, Or & o ) : m_st(s),m_nodes(n),m_order(o)
,m_ptr(0),m_size(m_order.size()) {}
node_iterator end() { m_ptr=m_size; return *this;  }
node_iterator begin() { m_ptr=0; return *this;  }
Node&  operator*() { return m_nodes[m_order[m_ptr]]; }  
bool valid()const  { return  m_ptr<m_size; }
bool next_ok()const  { return  (m_ptr+1)<m_size; }
IdxTy serial() const { return m_order[m_ptr]; } 
IdxTy i() const { return m_ptr; } 
void inc() { ++m_ptr; }
// this should compare the other members too I guess lol 
bool operator==( const node_iterator & that ) { return m_ptr==that.m_ptr; } 
bool operator!=( const node_iterator & that ) { return m_ptr!=that.m_ptr; } 
Myt & operator++() { ++m_ptr;  return *this; } 
operator bool() const { return valid(); } 
private:
St&  m_st;
Nodes&  m_nodes; 
const Or& m_order;
IdxTy m_ptr,m_size;
}; // node iterat0r;

public:
void add( const Lineage & h, const Values & v) {m_order.next();  m_nodes.push_back(Node(m_st(h),m_st(v))); }
void add( const Lineage & h, const Values & v, const Notes & notes) 
{m_order.next();  m_nodes.push_back(Node(m_st(h),m_st(v),m_st(notes))); }

void add( const Lineage & h, const Values & v, const Notes & notes,const bool flag, const IdxTy f=0) 
{m_order.next();  m_nodes.push_back(Node(m_st(h),m_st(v),m_st(notes)));
if (!flag) return; 
m_nodes.back().set_flag(f);
 }

const D dncvalue(const IdxTy node, const StrTy & sample) const
{
// apparently this is a VECTOR wtf??? 
return m_nodes[node].value(m_st[sample]);
}
void flag_last_node(const IdxTy f ) { m_nodes.back().set_flag(f); }


void set_last_src(const IdxTy src) { m_nodes.back().src(src); } 

void sort_alpha_lineage()
{
sort_alpha_lineage_comp salc(m_st,m_nodes);
m_order.sort(salc);
if (false ){MM_LOOP(ii,(*this)){MM_ERR(MMPR3(ii.i(),ii.serial(),informative_name(ii.serial(),0)))}}

}
void sort_alpha_notes()
{
sort_alpha_notes_comp salc(m_st,m_nodes);
m_order.sort(salc);
if (false ){MM_LOOP(ii,(*this)){MM_ERR(MMPR3(ii.i(),ii.serial(),informative_name(ii.serial(),0)))}}
combine_eqal_nodes(m_nodes,salc,begin(),end());
// updaye the ordering 
m_order.clear();
MM_LOOP(ii,m_nodes) { m_order.next(); } 
sort_alpha_notes_comp salc2(m_st,m_nodes);
m_order.sort(salc);

}


//class sort_values_comp : public  sort_base_lineage_comp
// this needs to use the dnc node values NOT the values object 
// that was the source for making the tree... 
//{m_order.next();  m_nodes.push_back(Node(m_st(h),m_st(v),m_st(notes))); }
template <class TFock>
void sort_values_fick(const StrTy & nm,const TFock & xnv)
{
//sort_values_comp salc(m_st,m_nodes,m_values[0]);
sort_values_fick_comp<TFock>  salc(m_st,m_nodes,nm,xnv); // ,m_values[0]);
m_order.sort(salc);
if (!false )
{MM_LOOP(ii,(*this)){MM_ERR(MMPR3(ii.i(),ii.serial(),informative_name(ii.serial(),0)))}}
}

void sort_values(const StrTy & nm, const StrTy &mods)
{
//sort_values_comp salc(m_st,m_nodes,m_values[0]);
sort_values_comp  salc(m_st,m_nodes,nm,mods); // ,m_values[0]);
m_order.sort(salc);
if (false )
{MM_LOOP(ii,(*this)){MM_ERR(MMPR3(ii.i(),ii.serial(),informative_name(ii.serial(),0)))}}
}







void sort_alpha_gs()
{
sort_alpha_gs_comp sagc( m_st, m_nodes, m_tt) ;
m_order.sort(sagc);

} 


template <class Te >
void combine_eqal_nodes(Nodes & d, Te & eq, node_iterator ii, const node_iterator  ee )
{
Nodes x;
if (ii!=ee){ x.push_back(*ii); ++ii; } 
while (ii!=ee)
{
Node & last=x.back();
if (eq.equals(last,*ii)) { last+=(*ii); } 
else { x.push_back(*ii); x.back().value_lineage(); } 
++ii;
}
d=x;

}


class layout_analysis
{

typedef mjm_uint_ordering  Sr;
typedef std::vector<IdxTy> LocArray;
typedef std::vector<D> ValArray;
typedef std::vector<point> PtArray;
typedef std::map<IdxTy, LocArray> Marks;
typedef std::map<IdxTy, ValArray> Vals;
typedef std::map<IdxTy, point>  MarkPts;
typedef Node::lineage Lg;

public:
class run
{
public:
run ( const StrTy & name, const IdxTy first) : m_name(name), m_first(first),m_last(first) {}
void include(const IdxTy end ) { m_last=end; } 
IdxTy len() const { return m_last+1-m_first; } 
bool includes(const IdxTy i ) const { return  (i>=m_first)&&(i<=m_last); } 
const StrTy m_name;
const IdxTy m_first;
// build as it goes along or call at the end...
IdxTy m_last;
}; 
typedef std::vector<run> Runs;
typedef std::map<StrTy, Runs> RunHier;

RunHier m_cat_runs;

static IdxTy bad() { return ~0; } 
//la.cat_runs(in,inold,(ii).i());
void cat_runs(const std::vector<StrTy>  & in,const std::vector<StrTy> & inold,const IdxTy idx)
{
Runs & rp= m_cat_runs["phylum"];
const IdxTy phloc=2;
const IdxTy loc=in.size()-phloc;
//const IdxTy locold=inold.size()-phloc;
if (in.size()>=phloc)
{
const StrTy & phylum=in[loc];
if (rp.size()==0) rp.push_back(run(phylum,idx));
run & ru=rp.back();
if ( ru.m_name==phylum) ru.include(idx);
else  rp.push_back(run(phylum,idx));
}
Runs & rp2= m_cat_runs["pseudogenus"];
const IdxTy pgloc=1;
if (in.size()>pgloc)
{
const StrTy & pg=in[pgloc];
if (rp2.size()==0) rp2.push_back(run(pg,idx));
run & ru=rp2.back();
if ( ru.m_name==pg) ru.include(idx);
else  rp2.push_back(run(pg,idx));

}


}

// lineage starts at terminal node and ends at root 
void mark_breaks(const Lg & lg,const Lg & lgold,const IdxTy idx)
{
IdxTy i=0;
m_la_break_levels[idx]=0; 
const IdxTy lgtop=lg.size()-1;
const IdxTy lgoldtop=lgold.size()-1;
if ( lgoldtop==bad()) {mark_idx(idx,i);  return ; } 
if ( lgtop==bad()) {mark_idx(idx,i); return ; } 
while (true)
{
if (( i>lgtop)||(i>lgoldtop)) {mark_idx(idx,i); return; }
const IdxTy v=lg[lgtop-i];
const IdxTy vold=lgold[lgoldtop-i];

if ( v!=vold) { mark_idx(idx,i); return; }
++i;
}
}
void mark_idx(const IdxTy idx, const IdxTy i)
{
m_la_break_levels[idx]=i+1; m_la_breaks[i].push_back(idx);

}

template <class Ti > 
void leg_points(Nodes & nodes, Ti _ii, Ti ee)
{
D ylowest=0;
D ymin=.1;
std::map<IdxTy,IdxTy> attached;
const IdxTy nosz=nodes.size();
MM_SZ_LOOP(j,m_sample_order,sosz) 
//for  (; ii!=ee; ++ii )
{
IdxTy spn= D(sosz)/(.2*D(nosz));
if (spn<1) spn=1;
const IdxTy sample=m_sample_order[j];
D xmin=.4*nosz;
D xmax=.75*nosz;
D ybest=0;
D xbest=0;
//MM_LOOP(ii,nodes){ // nodes 
for (auto ii=_ii; ii!=ee; ++ii) {
D x=ii.i(); // (ii).i(); // this is a node number not a location 
D y=m_int[sample][(ii).i()]; // /=norms[sample];
if (x<xmax) if (attached[ii.i()]>=spn ) continue;
if (x>xmax) if ( ybest==0) 
{
xbest=x; ybest=y; 
break;
}
if (x<xmin) continue;
if (y<ymin) continue;
if (y>ylowest){  xbest=x; ybest=y; ylowest=y; break;  } 
} 
if (ybest!=0 ) 
{  
	m_leg_pts[sample] = point(xbest,ybest); 
	++attached[IdxTy(xbest)];
}

}
//MM_SZ_LOOP(j,m_sample_order,sosz)
// {
//const IdxTy sample=la.m_sample_order[j];
//}

} // leg_points

const D & leg_x(const IdxTy& sample )  { return m_leg_pts[sample].x; }
const D & leg_y(const IdxTy& sample )  { return m_leg_pts[sample].y; }


const IdxTy nsamples() const { return m_sample_order.size(); } 
StrTy sample_color(const IdxTy & sample) { return StrTy("red"); } 

std::vector<StrTy> m_sample_names;
IdxTy m_n,m_max_sample,m_n_sample_chars;
Sr m_sample_order;
Marks m_la_breaks;
std::map<IdxTy,IdxTy> m_la_break_levels;
Vals m_int;
MarkPts m_leg_pts;
ValArray m_means;
ValArray m_minns;
ValArray m_maxxs;
}; // layout_analysis
class layout_consts
{

typedef std::map<IdxTy, point>  MarkPts;
typedef std::vector<StrTy> StrVec;
public:
layout_consts(  ) {Init(); } 
layout_consts(const std::map<StrTy,StrTy>  & ph  ) 
{
	m_ph=ph;
	Init(); 

} 
void do_analysis( layout_analysis & la ) {Init(la); } 
layout_consts( layout_analysis & la ) {Init(); Init(la); } 
const StrTy& doc_name() const { return m_doc_name; }
const StrTy& sign() const { return m_doc_name; }
const StrTy& bg_color() const { return m_bg_color; }
const D& total_width() const { return m_total_width; }
const D& total_iia_width() const { return m_total_width; }
const D& total_hm_width() const { return m_total_hm_width; }
const D& total_height() const { return m_total_height; }
const D& total_iia_height() const { return m_total_height; }
const D& total_hm_height() const { return m_total_height; }

const D& plot_width() const { return m_plotw; }
const D& plot_height() const { return m_ploth; }

const D& plot_x() const { return m_l_margin; }
const D& plot_y() const { return m_t_margin; }
const D& bottom_margin() const { return m_b_margin; }
const D& top_margin() const { return m_t_margin; }
const D& right_legend_width() const { return m_r_margin; }

const D& node_size() const { return m_node_width; }  
const D& node_label_width() const { return m_node_label_width; }  
const IdxTy& nodes() const { return m_nodes; }  
const IdxTy& samples() const { return m_samples; }  

const D& stripe_thickness() const { return m_stripe_thickness;}
const StrTy& stripe_color() const { return m_stripe_color;}
const D& sign_height() const { return m_sign_height;}
const StrTy& sign_color() const { return m_sign_color;}
const D& sign_x() const { return m_sign_x; }
const D& sign_y() const { return m_sign_y; }


const StrTy& ylab() const { return m_ylab;}
const D& ylab_x() const { return m_ylab_x; }
const D& ylab_y() const { return m_ylab_y; }
const D& ylab_size() const { return m_ylab_size; }
const D& curve_width() const { return m_curve_width; }

const D&  sample_legend_font() const { return m_leg_font; } 
// TODO FIXME this is faster in a loop 
//const D leg_x(const IdxTy i ) const { return 0; }
//const D leg_y(const IdxTy i ) const { return 0; }
//MM_SZ_LOOP(j,la.m_sample_order,sosz) {
// This will be a different sample ordering but still have an iterator 
template<class Tm,class To> 
void sample_legend_pos(Tm & vx, Tm & vy, const To & so)
{
IdxTy ncols=m_n_leg_cols;
IdxTy nrows=m_n_leg_rows; 
MM_SZ_LOOP(i,so,sosz) 
{
//const IdxTy col= i/nrows;
//const IdxTy row= i%nrows;

const IdxTy col= i%ncols;
const IdxTy row= i/ncols;

D x= m_leg_x+ m_leg_col*col;
D y= m_leg_y+m_leg_row*row+m_leg_font;
vx.push_back(x);
vy.push_back(y);
}
}
template<class Tm,class To> 
void hm_sample_legend_pos(Tm & vx, Tm & vy, const To & so)
{
IdxTy ncols=1; // m_n_leg_cols;
IdxTy nrows=m_samples; // m_n_leg_rows; 
const D xpos=m_plotw+m_l_margin+5;
MM_SZ_LOOP(i,so,sosz) 
{
//const IdxTy col= i/nrows;
//const IdxTy row= i%nrows;
//const IdxTy col= i%ncols;
//const IdxTy row= i/ncols;
//D x= m_leg_x+ m_leg_col*col;
//D y= m_leg_y+m_leg_row*row+m_leg_font;
// this fails with the m_hm_sample_order crap
//D y=hm_sample_y(i)+.5*(hm_sample_y(i+1)-hm_sample_y(i));
D y=hm_sample_y(i)+.5*(hm_sample_pitch());
//D y=hm_sample_y(so[i])+.5*(hm_sample_pitch());
vx.push_back(xpos);
vy.push_back(y);
}
}

const D hm_sample_pitch() const
{
D h=D(plot_height())/m_samples; // la.m_sample_order.size();
return h;
} 

std::map<IdxTy,IdxTy>& ma= m_hm_sample_order;
const D  hm_sample_y( const IdxTy i) const
{
//D h=lc.plot_height()/la.m_sample_order.size();
D h=D(plot_height())/m_samples; // la.m_sample_order.size();
return plot_y()+h*i;
//return plot_y()+h*m_hm_sample_order[i];

}

const D & sample_ptr_size() const { return m_sample_ptr_size; } 
const StrTy & sample_color(const IdxTy i) const { 
if (i>=m_sample_colors.size()) { MM_ERR(" sample out of range "<<MMPR2(i,m_sample_colors.size())) }
//MM_ERR(MMPR2(m_sample_colors.size(),i))
return m_sample_colors[i] ; } 
const StrTy & sample_contrast_color(const IdxTy i) const { 
if (i>=m_sample_contrast_colors.size()) { MM_ERR(" sample out of range "<<MMPR2(i,m_sample_contrast_colors.size())) }
//MM_ERR(MMPR2(m_sample_colors.size(),i))
return m_sample_contrast_colors[i] ; } 

const IdxTy & flag_number(StrTy & name ) //const 
{ 
if (m_flag_numbers.find(name)==m_flag_numbers.end())
{
add_flag(name,default_colors(nflags()));
}
return m_flag_numbers[name]; } 



const StrVec & flag_colors() const { return m_flag_colors; } 
const StrVec & flag_names() const { return m_flag_names; } 
const StrTy & flag_color(const IdxTy i) const { 
if (i>=m_flag_colors.size()) { MM_ERR(" sample out of range "<<MMPR2(i,m_flag_colors.size())) }
return m_flag_colors[i%(m_flag_colors.size())];
} 
const StrTy  flag_name(const IdxTy i) const { 
const IdxTy sz=m_flag_names.size();
if (i>=sz) { MM_ERR(" sample out of range "<<MMPR2(i,m_flag_names.size())) }
if (sz==0) return StrTy();
return m_flag_names[i%(sz)];
} 
const StrTy & node_iia_color(const IdxTy i) const { 
if (i>=m_flag_colors.size()) 
{ MM_ERR(" iia node  out of color range "<<MMPR2(i,m_flag_colors.size())) }
return m_flag_colors[i%(m_flag_colors.size())];
} 
const StrTy & node_iia_contrast_color(const IdxTy i) const { 
const auto & v=m_flag_colors;
const IdxTy cyclo=v.size();
IdxTy comp=contrast(i,cyclo); // (cyclo-1+cyclo*i-i)%(cyclo);
//if (comp==i) comp=(comp+1)%cyclo;
if (i>=cyclo) 
{ MM_ERR(" iia node  out of color range "<<MMPR2(i,m_flag_colors.size())) }
return v[comp];
} 

IdxTy contrast(const IdxTy& i, const IdxTy& cyclo) const
{
IdxTy comp=(cyclo-1+cyclo*i-i)%(cyclo);
if (comp==i) comp=(comp+1)%cyclo;
return comp;
}



const D  sample_iia_x(const IdxTy i ) const
{
return plot_x()+i*m_sample_iia_dx;
}
const D & sample_iia_dx(const IdxTy i ) const
{return m_sample_iia_dx; } 
const D  sample_iia_filldx(const IdxTy i ) const
{return m_sample_iia_dx*.5; } 



const IdxTy&  nflags() const { return m_nflags; }
//const IdxTy ntflag=lc.add_flag("assugnment","green");
const IdxTy add_flag(const StrTy & name, const StrTy & color) // "assugnment","green");
{
if (m_flag_numbers.find(name)!=m_flag_numbers.end()) { return m_flag_numbers[name]; } 
const StrTy _color=(color.length()==0)?default_colors(nflags()):color;
const IdxTy flag= m_nflags;
m_flag_numbers[name]=m_flag_names.size();
m_flag_colors.push_back(_color);
m_flag_names.push_back(name);

++m_nflags;
return flag;
}

//std::map<IdxTy,IdxTy>& ma= m_hm_sample_order;
const StrTy & sample_name(const IdxTy i) const 
{ return m_sample_names[m_hm_sample_order[i]] ; } 
//{ return m_sample_names[i] ; } 
const D & leg_x(const IdxTy& sample )  { return m_leg_pts[sample].x; }
const D & leg_y(const IdxTy& sample )  { return m_leg_pts[sample].y; }
const IdxTy & name_levels( ) const  { return m_name_levels; }


const StrTy & terminal_name() const { return m_terminal_name; } 

//D x0=ii.i()*lc.node_size()+lc.plot_x();
//D y0=lc.plot_y()+lc.plot_height()*mean;
// fudding +1 is a svg kluge 
const D node_label_x(const IdxTy i) const { return node_size()*i+plot_x()+1+IdxTy(node_size()/4); } 
const D node_start_x(const IdxTy i) const { return node_size()*i+plot_x(); } 
const D node_label_y(const D&  mean ) const { return plot_height()*mean +plot_y(); } 
enum  SORT_KEY { SORT_ALPHA_TREE=0, SORT_ALPHA_GS=1, SORT_NOTES=2,SORT_VALUES=3 } ; 
//const SORT_KEY sort_key() const { return SORT_ALPHA_TREE; } 
const SORT_KEY sort_key() const { return SORT_VALUES; } 
const StrTy sort_sample() const 
{
const auto ii=m_ph.find("sort_sample");
if (ii==m_ph.end()) return StrTy("");
return (*ii).second; 

}
private:


const StrTy & default_colors(const IdxTy i) 
{
static std::vector<StrTy> cs;
static bool init=false;
if (!init)
{
cs.push_back("white");
cs.push_back("red");
cs.push_back("green");
cs.push_back("orange");
cs.push_back("blue");
cs.push_back("black");
init=true;
}
return cs[i%cs.size()];
}



void Init()
{

m_doc_name="Marchywka Test";
m_bg_color="#0808080";
// not eliminated 
m_terminal_name="DUMMY";
// eliminated by ignored lsit
//m_terminal_name="MISSING";
m_sign_color="black";

m_nflags=0;
}
public:

// so is sequence maybe with skips and name
void hm_sample_sort(const std::map<IdxTy, StrTy >&  so, layout_analysis & la, St & st)
{
// name=x[idx] 
std::map<IdxTy,IdxTy>& ma= m_hm_sample_order;
std::vector<StrTy> no;
std::map<StrTy,IdxTy> mord,haves;
MM_SZ_LOOP(i,m_sample_names,snsz) haves[m_sample_names[i]]=i;
MM_LOOP(ii,so) 
{
const bool ha= (haves.find((*ii).second)!=haves.end()); 
	MM_ERR(" suggested order "<<MMPR3((*ii).first, (*ii).second,ha))
	if (ha)	 no.push_back((*ii).second);
}

MM_SZ_LOOP(i,no,xsz) mord[no[i]]=i;
if (m_samples!=no.size())
{
MM_ERR(" sample count wrong "<<MMPR3(so.size(),no.size(),m_samples))
}
for (IdxTy i=0; i<m_samples; ++i)
{
const StrTy & nm=m_sample_names[i];
const IdxTy suggloc=mord[m_sample_names[i]];
// this should worl FICL 
ma[suggloc]=i;
//ma[i]=suggloc;
MM_ERR(" sample order "<<MMPR4(i,suggloc,nm,la.nsamples()))
} 
std::vector<IdxTy> nr;
//MM_LOOP(ii,ma) nr.push_back((*ii).second);
MM_LOOP(ii,ma)
{
 nr.push_back(st(m_sample_names[(*ii).second]));
MM_ERR(" final ord er "<<MMPR(m_sample_names[(*ii).second]))
}
//MM_LOOP(ii,no) nr.push_back((*ii));
la.m_sample_order.clear();
la.m_sample_order.add(nr);
MM_ERR(" sample order  backout ")
//hm_sample_sort(false); 
MM_ERR(" sample order fini ")
}
void hm_sample_sort(const bool sugo=false)
{
MM_ERR( " sorting ghm sample sort ")
std::map<StrTy,IdxTy> mord;
//{ return m_sample_names[m_hm_sample_order[i]] ; } 
// name=x[idx] 
StrVec & x= m_suggested_initial_sample_order;
MM_SZ_LOOP(i,x,xsz) mord[x[i]]=i;
std::map<IdxTy,IdxTy>& ma= m_hm_sample_order;
for (IdxTy i=0; i<m_samples; ++i)
{
if (!sugo) { ma[i]=i; continue; } 
//const StrTy & nm=m_sample_names[i];
//if (mord.find(nm)==mord.end()) continue;
//const IdxTy suggloc=mord[m_sample_names[i]];
//ma[i]=suggloc;
//MM_ERR(" hm_sample_sort "<<MMPR3(i,suggloc,x[i]))
}
MM_ERR( " done sorting ghm sample sort ")
}
private:
void Init( layout_analysis & la ) 
{

m_nodes=la.m_n;
m_samples=la.nsamples();
hm_sample_sort();
m_aspect_ratio=1.0/3.0;
m_node_width=1;
m_name_levels=2;
//m_node_size=1;
// try to keep size between about 500 and 2000
const D tgtns=1000.0/(m_nodes+1);
if (m_nodes<100) m_node_width=10;
m_node_width=IdxTy(tgtns+1);
if (m_node_width<1) m_node_width=1;
if (m_node_width>10) m_node_width=10;

m_node_label_width=.5+.5*m_node_width;
if ( m_node_label_width>6) m_node_label_width=6;
m_plotw=m_nodes*m_node_width;
m_ploth=IdxTy(m_plotw*m_aspect_ratio);
m_l_margin=10;
m_r_margin=20;
m_t_margin=20;
m_b_margin=20+.1*m_nodes*m_node_width;
m_total_width=m_plotw+m_l_margin+m_r_margin;
m_hm_sample_width=100;
m_total_hm_width=m_total_width+m_hm_sample_width; // m_plotw+m_l_margin+m_r_margin;
m_total_height=m_ploth+m_b_margin+m_t_margin;

m_leg_x=m_l_margin;
m_leg_y=m_t_margin+m_ploth+2;
m_leg_w=m_plotw;
m_leg_h=m_b_margin-4;

D areamax=m_samples*la.m_max_sample;
if (areamax==0) areamax=1;
D layarea=m_leg_w*m_leg_h;
D fsz=1+.5*sqrt(layarea/(areamax+1));
MM_ERR(MMPR(fsz))
if (fsz>m_leg_h) fsz=.5*m_leg_h;
m_leg_col=fsz*la.m_max_sample;
if (m_leg_col<1) m_leg_col=1;
if (m_leg_col>3) m_leg_col-=1;
m_leg_row=1.2*fsz;
m_n_leg_rows=IdxTy(m_leg_h/fsz);
if (m_n_leg_rows<1) m_n_leg_rows=1;
m_n_leg_cols=IdxTy(m_leg_w/m_leg_col);
if (m_n_leg_cols<1) m_n_leg_cols=1;
m_leg_font=fsz;

m_stripe_color="black";
m_stripe_thickness=.3;
m_sign_y=m_ploth+m_t_margin+.8*m_b_margin;
m_sign_w=.5*m_plotw;
m_sign_h=.2*m_b_margin;
m_sign_height=.5*m_sign_h; 
IdxTy doc_len= m_doc_name.length();
if (doc_len*m_sign_height> m_sign_w) m_sign_height=1+IdxTy(m_sign_w/doc_len);
m_sign_x=m_l_margin+m_plotw-doc_len*m_sign_height; // .25*m_plotw+m_l_margin;
if (m_sign_x<0) m_sign_x=0;

m_ylab=" Cumulative Hits Fraction Explained";
m_ylab_x=m_plotw+m_l_margin+.1*m_r_margin;
m_ylab_y=.2*m_ploth+m_t_margin;
m_ylab_size=.6*m_r_margin;
doc_len= m_ylab.length();
if (doc_len*m_ylab_size> .6*m_ploth) m_ylab_size=1+IdxTy(.6*m_ploth/doc_len);
m_curve_width=1.9; // *m_node_width;
m_sample_ptr_size=.2;
m_sample_names=la.m_sample_names;
std::vector<StrTy> cs;
cs.push_back("white");
cs.push_back("red");
cs.push_back("green");
cs.push_back("orange");
cs.push_back("blue");
cs.push_back("black");
// TODO FIXME temp for test 
//m_flag_colors=cs;
const IdxTy cycle=cs.size();

MM_SZ_LOOP(i,m_sample_names,snsz)
{
m_sample_colors.push_back(cs[i%cycle]);
m_sample_contrast_colors.push_back(cs[(cycle+i*cycle-i-1)%cycle]);
}

m_leg_pts=la.m_leg_pts;
MM_LOOP(ii,m_leg_pts)
{
point & p=(*ii).second;
p.x=m_l_margin+m_node_width*(p.x+1);
p.y=m_t_margin+m_ploth*p.y;

}


 m_sample_iia_dx=1.0*m_plotw/m_samples; 

MM_ERR(dump())
{MM_SZ_LOOP(i,la.m_sample_order,szs)  { MM_ERR(MMPR2(i,(la.m_sample_order[i]))<<" ") } }
} 

StrTy dump()
{
Ss ss;

ss<<MMPR4( m_doc_name, m_bg_color,m_stripe_color,m_sign_color); 
ss<<MMPR2(m_terminal_name,m_nflags);
ss<<MMPR4(m_node_width,m_node_label_width, m_nodes,m_samples);
ss<<MMPR4(m_ylab, m_nodes,m_samples,m_total_width);
ss<<MMPR4(m_total_height,m_node_width, m_aspect_ratio,m_l_margin);
ss<<MMPR4( m_r_margin,m_t_margin, m_b_margin, m_plotw);
ss<<MMPR4( m_ploth ,m_leg_x,m_leg_y,m_leg_w);
ss<<MMPR4(m_leg_h,m_leg_font,m_leg_col,m_leg_row)<<MMPR2(m_n_leg_rows,m_n_leg_cols)<<MMPR3( m_stripe_thickness,m_sign_height,m_sign_x);
ss<<MMPR4(m_sign_y, m_ylab_x,m_ylab_y,m_ylab_size);
ss<<MMPR3(m_curve_width,m_sample_ptr_size,m_name_levels);
{MM_SZ_LOOP(i,m_sample_names,szs)  { ss<<MMPR2(i,m_sample_names[i])<<" "; } }
{MM_SZ_LOOP(i,m_sample_colors,szs)  { ss<<MMPR2(i,m_sample_colors[i])<<" "; } }
{MM_SZ_LOOP(i,m_sample_contrast_colors,szs)  { ss<<MMPR2(i,m_sample_contrast_colors[i])<<" "; } }
{MM_LOOP(ii,m_leg_pts)  { ss<<MMPR2((*ii).first,StrTy(m_leg_pts[(*ii).first]))<<" "; } }
return ss.str();
}

StrTy m_doc_name, m_bg_color,m_stripe_color,m_sign_color,m_ylab,m_terminal_name;
IdxTy m_nodes,m_samples,m_n_leg_rows,m_n_leg_cols,m_name_levels,m_nflags;

D m_total_width,m_total_height,m_node_width, m_node_label_width, m_aspect_ratio,m_l_margin, m_r_margin,m_t_margin, m_b_margin;
D m_total_hm_width,m_hm_sample_width;
D m_plotw, m_ploth; // , m_node_size;
D m_leg_x,m_leg_y,m_leg_w,m_leg_h,m_leg_font,m_leg_col,m_leg_row;
D m_stripe_thickness,m_sign_height,m_sign_x,m_sign_y,m_sign_w,m_sign_h;
D m_ylab_x,m_ylab_y,m_ylab_size;
D m_curve_width,m_sample_ptr_size;

StrVec m_sample_names,m_sample_colors,m_sample_contrast_colors,m_flag_colors,m_flag_names;
std::map<StrTy, IdxTy> m_flag_numbers;
 MarkPts m_leg_pts;

D  m_sample_iia_dx; 
public:
StrVec m_suggested_initial_sample_order;
mutable std::map<IdxTy,IdxTy>  m_hm_sample_order;

std::map<StrTy,StrTy> m_ph;
}; // layout_consts

public:

typedef layout_analysis analysis;
typedef layout_consts consts;

void analyze( analysis  & la, consts & lc)
{
typedef Node::lineage Lg;
Lg lgold;
la.m_n=m_nodes.size();

std::vector<StrTy>  inold;



std::map<IdxTy,IdxTy> mso;
std::map<IdxTy,D> norms;
// this multi pass thing is stupid but easy to code...

// this is a string vector
bool used_order= 
	(lc.m_ph["manual-sample-order"]!="") &&(lc.m_suggested_initial_sample_order.size()!=0);
if (used_order)
{
MM_LOOP(ii,lc.m_suggested_initial_sample_order)
{
// this needs to eliminated dead samples or later normalization will bomb...
const IdxTy isamp=m_st(*ii);
MM_LOOP(jj,(*this) ) { if ((*jj).values().find(isamp)!=(*jj).values().end()){ 
la.m_sample_order.add_unique_value(isamp); break; } } 
} // ii 

}
else{ 
MM_LOOP(ii,(*this)) { MM_LOOP(jj,(*ii).values()) { la.m_sample_order.add_unique_value((*jj).first); }}
}

const IdxTy ntflag=lc.add_flag("assugnment","");
// this is the wrong fudding string tokenizer, nees to use tax tre. 
const IdxTy fakeno=m_st(StrTy("fake"));
MM_LOOP(ii,(*this)) // nodes 
{
 auto & node=(*ii);
if (false) 
	{ MM_ERR(" ASSFUDD "<<MMPR4(fakeno,ii.i(),ii.serial(),int(node.src()))) } 
if (node.src()==fakeno) { node.set_flag(ntflag); }
} // ii 


MM_LOOP(ii,(*this)) // nodes 
{
const auto & node=(*ii);
const auto & v=node.values();
Lg lg=node.names();
std::vector<StrTy>  in= m_tt.informative_vector(m_st(lg));
//la.mark_breaks(lg,lgold,(ii).serial());
la.mark_breaks(lg,lgold,(ii).i());
// TODO FIXME his needs to assemble on the fly or be called at end to get the termina ones
la.cat_runs(in,inold,(ii).i());
lgold=lg;
inold=in;
/// find hieratchy breaks, flags, and integrate samples. 
MM_SZ_LOOP(j,la.m_sample_order,sosz) {
const IdxTy sample=la.m_sample_order[j];
auto fudd=v.find(sample);
D dv=(fudd==v.end())?0:(*fudd).second; // v[sample];
auto &val=la.m_int[sample];
const D vold=(val.size()>0)?val.back():0;
const D vnew=vold+dv;
val.push_back(vnew);
norms[sample]=vnew;
} // values
} // ii 

MM_LOOP(ii,(*this)){ // nodes 
D sum=0; IdxTy n=0;
D minn=1e10; D maxx=0;
MM_SZ_LOOP(j,la.m_sample_order,sosz) {
const IdxTy sample=la.m_sample_order[j];
la.m_int[sample][(ii).i()]/=norms[sample];
//MM_ERR("NORMALIZ FUICK "<<MMPR4(n, la.m_int[sample][(ii).i()],norms[sample],m_st(sample) )<<MMPR(ii.i()))
//++n; sum+=la.m_int[sample][ii.i()];
++n; 
D s=la.m_int[sample][ii.i()];
sum+=s;
if ( s<minn ) minn=s;
if ( s>maxx) maxx=s;
}
const D avg=(n==0)?0:(sum/n);
la.m_means.push_back(avg);
la.m_minns.push_back(minn);
la.m_maxxs.push_back(maxx);
}

la.m_max_sample=0;
la.m_n_sample_chars=0;
MM_SZ_LOOP(j,la.m_sample_order,sosz) {
const IdxTy sample=la.m_sample_order[j];
const StrTy&  sn= m_st(sample);
const IdxTy snlen=sn.length();
la.m_n_sample_chars+=snlen;
if (snlen>la.m_max_sample) la.m_max_sample=snlen;
}


la.leg_points(m_nodes,begin(),end());

// TODO this needs to scope the fing size fudd 
MM_SZ_LOOP(j,la.m_sample_order,sosz2) {
const IdxTy sample=la.m_sample_order[j];
MM_ERR(MMPR4(la.m_sample_names.size(),j,sample,m_st(sample)))
la.m_sample_names.push_back(m_st(sample));
}


} // analyze


typedef node_iterator iterator;
node_iterator begin() { return node_iterator(m_st,m_nodes,m_order); }
node_iterator end() { return node_iterator(m_st,m_nodes,m_order).end(); }
StrTy informative_name(const IdxTy i,const IdxTy lvl ) { return m_tt.informative_n_name(m_st(m_nodes[i].names()),lvl); } 
St & st() { return m_st; } 
private:
St m_st;
Nodes m_nodes; 
Or m_order;
mjm_taxon_tools m_tt;
}; // data_node_collection 





}; // ns 




#if 0 
class layout_consts {
typedef  ::tree_viz_traits::Tr  Tr;
typedef layout_consts Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef mjm_taxon_tools Mtt;
typedef Mtt::ignore_t ignore_t;

public:
layout_consts(const IdxTy _node_count):
 node_count(_node_count),
nodew(fscale(node_count)),
 w(node_count*nodew),
 hscale(3.0/4.0*.5),
 h(w*hscale ), 
 ymargin(.15),
 ylegend(0.02*h ),
 ytaxa0(.3*h ),
ytaxa1(.8*h),
 xoff(0),
 yoff(0 ),
// itaxa(0),
 linesz(.002*h),
 nodesz(rmin(.8,2.5/nodew)),
 visthick(.001*h),
 taxaszconst(bignode(.1*h*.75*nodew,nodew)),
 taxasz(bignode(rmin(5.0*nodew,.04*h*nodew),nodew)),
tw(1.05*w),th((1.0+ymargin)*h)
 {}
static D rmin(const D & x1, const D & x2) { return (x1<x2)?x1:x2;}
static D rmax(const D & x1, const D & x2) { return (x1>x2)?x1:x2;}
D bignode(const D & x1, const D & x2) { return (nodew>1)?x2:x1;}
static IdxTy fscale(const IdxTy & nc) { if (nc>100) return 1;  return 10; }

//ow(x1,y1,np.x1,lc.taxlocedge(right,alt));
D taxlocedge(const bool right, const bool alt ) const
{
if (right&&alt) return .01*h;
if (right&&!alt) return .11*h;
if (!right&&alt) return .99*h;
if (!right&&!alt) return .89*h;
return 0;
}

const ignore_t & ignores() const { return m_mtt.ignores(); }

template <class Tt> 
StrTy informative_name(const  Mtt::tax_t & tv, const Tt & tt) const
{ return  m_mtt.informative_name(tv,tt) ; }
private:
public:
const D node_count,nodew, w , hscale, h , ymargin, ylegend,  ytaxa0, ytaxa1, xoff, yoff;
const D linesz, nodesz, visthick, taxaszconst, taxasz;
const D tw,th;

mjm_taxon_tools m_mtt;

}; // layout_consts

typedef Tr::D D; 
typedef Tr::IdxTy IdxTy; 
typedef Tr::StrTy StrTy  ;  

class wpos { public: wpos():x(0),y(0),color("black") {} D x,y; StrTy color; StrTy nm; };
class npos { public: npos():x0(0),x1(0) {} npos(const D & i ) : x0(i),x1(i+1),ntaxon(0){} D x0,x1; IdxTy ntaxon;  };

class legmap : public std::map<StrTy,wpos>
{
public: 
void add(const StrTy & n, const wpos & w) { (*this)[n]=w;}
class compc{ public: 
static bool cmp(const wpos & x, const wpos & y) { return x.x<y.x; }
bool comp(const wpos & x, const wpos & y) { return x.x<y.x; }
template <class Tfudd > bool operator()(const Tfudd  & x, const Tfudd  & y) { return x.y<y.y; }

};
typedef std::vector<wpos> sorted;
sorted sort()
{
sorted r; 
compc compcfudd;
MM_LOOP(ii,(*this)) { (*ii).second.nm=(*ii).first;  r.push_back((*ii).second); } 
std::sort(r.begin(),r.end(),compcfudd);
return r; 
} 

}; // legmap 


class layout_state
{
typedef  ::tree_viz_traits::Tr  Tr;
typedef layout_state Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef legmap LegendMap; //  legend_map;


public:
class node_layout : public  std::map<IdxTy,npos>  // NodeLayout;
{
public:
typedef std::vector<IdxTy> Locs;
typedef std::vector<IdxTy> Cursor;
typedef std::vector<Locs> Breaks;
typedef std::vector<IdxTy> Vals;
static const IdxTy & bad() { static const IdxTy v= ~0U;  return v; } 
IdxTy breaks(const IdxTy & loc) const
{
return bad(); 

}
const Cursor&  cursor(const IdxTy & i )
{
return m_breaks[i];
}
// this is a LINEAGE,the first value is lowest level.
void new_values(const Vals & _v, const IdxTy & loc)
{
Vals v;
reverse(v,_v);
const IdxTy newsz=v.size();
const IdxTy oldsz=m_vals.size();
if (newsz>oldsz) add_lvl(bad(),loc,newsz-oldsz); 
MM_SZ_LOOP(i,v,vsz) { 
// lowest level change invalidates lower levels 
if (v[i]!=m_vals[i]) {add_break(i,loc); break; }
}
// this only need copy those that changed... 
m_vals=v;
}
void reverse(Vals & v, const Vals & _v)
{ MM_SZ_LOOP(i,v,sz) { v.push_back(_v[sz-i-1]); } }

void add_lvl(const IdxTy & val,  const IdxTy& loc, const IdxTy n=1 )
{
Locs nv;
if (loc!=bad()) nv.push_back(loc); 
for (IdxTy i=0; i<n; ++i)
{ m_vals.push_back(val);  m_breaks.push_back(nv); }
}

void add_break(const IdxTy&  lvl, const IdxTy& loc)
{
Locs nv;
//nv.push_back(loc);
while (m_breaks.size()<=lvl){ add_lvl(bad(),bad()); }
for (IdxTy i=lvl; i<m_breaks.size(); ++i) m_breaks[i].push_back(loc); 

}

Breaks m_breaks;
Vals m_vals;
}; // node_layout

// this is all in layout_state

typedef std::map<StrTy,wpos> SampleCursor;
//typedef std::map<IdxTy,npos> NodeLayout;
typedef node_layout NodeLayout;
typedef std::map<StrTy,D> SampleFig;
layout_state(): itaxa(0) {}
typedef mjm_svg_writer sw_type;
mjm_svg_writer sw;
typedef NodeLayout nl_type;
typedef NodeLayout::Cursor  LayoutCursor;
NodeLayout nl;
//LayoutCursor nlc;
typedef SampleCursor sc_type;
SampleCursor sc;
typedef SampleFig sv_type;
SampleFig sv;
LegendMap legend_map;
IdxTy itaxa;
}; //layout_state


class node_values  // : public  std::map<IdxTy, NodeValue> NodeValues;
{

typedef node_values Myt;
typedef  ::tree_viz_traits::Tr  Tr;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy,D> NodeValue;
typedef std::vector<StrTy> SampleOrder;
typedef std::vector<IdxTy> NodeOrder;
typedef std::map<IdxTy, NodeValue> NodeValues;
typedef NodeValues::const_iterator nv_iterator;

public:
typedef NodeValue node_value;
typedef std::vector<D> values_t;
typedef SampleOrder sample_order_t;
const node_value & operator[](const IdxTy n) { return m_nv[n];} 
void add(const IdxTy node, const StrTy &sample, const D&  x)
{
m_nv[node][sample]+=x;
if (m_samples.find(sample)==m_samples.end()) m_order.push_back(sample); 
m_samples[sample]+=x;
}

IdxTy nsamples() const { return m_samples.size(); } 
const sample_order_t & order() const { return m_order; }  
void move(const IdxTy nn,const IdxTy nt)
{
m_nv[nn]=m_nv[nt];
m_nv[nt].clear();
}

nv_iterator find(const IdxTy n) { return m_nv.find(n); } 
nv_iterator begin() { return m_nv.begin(); } 
nv_iterator end() { return m_nv.end(); } 

values_t values(const IdxTy n) const
{
values_t v;
auto jj=m_nv.find(n);
// all of nothing, otherwise no diea which are missing 
if (jj==m_nv.end()) return v;
auto&  snv=(*jj).second;
//MM_LOOP(ii,m_samples)
MM_LOOP(oo,m_order)
{
auto ii=m_samples.find(*oo);
auto kk=snv.find((*ii).first);
if (kk==snv.end()) { v.push_back(0); } 
else { v.push_back((*kk).second/(*ii).second) ; } 
//v.push_back(((*jj).second)[(*ii).first]/(*ii).second) ;
}

return v; 
}

void print_values(Ss & ss,const IdxTy & node,const IdxTy & digits
		,const StrTy & sep,const IdxTy & flags)
{
std::vector<D> values=(*this).values(node);
ss.precision(digits);
MM_LOOP(ii,values) {
if ((*ii)==0) ss<<sep<<" 0 "; 
 else { 
const D v=(*ii);
ss<<std::scientific<<sep<<v; }}
}

class sorted_nodes
{
public:
sorted_nodes(): m_values(0),m_ptr(0), m_sample() {}
class sv
{
public:
template <class Ty> 
bool operator()(const Ty & a, const Ty & b) const
{ // TODO FIXME these values need to be normalized for cross sample tests 

auto ii=nvr.find(a);
// TODO not sure this is exactly right but ok for now 
if (ii==nvr.end()) return false; 
auto jj=nvr.find(b);
if (jj==nvr.end()) return true; 
auto kk=(*ii).second.find(m_sample);
if (kk==((*ii).second.end())) return false;
auto ll=(*jj).second.find(m_sample);
if (ll==((*jj).second.end())) return true;

return ((*kk).second<((*ll).second)); 
//return ( nvr[a][m_sample]<nvr[b][m_sample]); 

}
StrTy m_sample;
const NodeValues & nvr;
}; // sv
bool is_valid() const { return m_ptr<m_node_order.size(); } 
NodeValues * m_values;
IdxTy m_ptr;
StrTy m_sample;
NodeOrder m_node_order;

}; // sorted_nodes

NodeValues m_nv;
std::map<StrTy,D> m_samples;
SampleOrder m_order;
}; // node_values


}; // tree_viz_traits



class mjm_tree_viz 
{
typedef  tree_viz_traits::Tr  Tr;
typedef mjm_tree_viz Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef tree_viz_params Logic;
typedef mjm_logic_base VariableStore;


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;



typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;

typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;

// TODO FIXME this needs to be a more comlicated object 
//typedef std::map<StrTy,D> NodeValue;
//typedef std::map<IdxTy, NodeValue> NodeValues;
typedef tree_viz_traits::node_values  NodeValues;
typedef NodeValues::node_value NodeValue;

typedef std::vector<StrTy> TaxVec;

typedef tree_viz_traits::layout_consts LayoutConst;

typedef tree_viz_traits::wpos  wpos;
typedef tree_viz_traits::npos  npos;
typedef tree_viz_traits::legmap  LegendMap;

public:
// TODO FIXME put this some where in FRP ... 
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }
// TODO FIXME find FRP for this.. 
static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

public :
bool done() const  { return m_done; } 
#endif


#undef MM_DMEL

#endif

