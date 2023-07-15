#ifndef MJM_DOG_PLOTS_H__
#define MJM_DOG_PLOTS_H__
 
#include "mjm_globals.h"
#include "mjm_canned_methods.h"
// add day number to notes 
#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_svg_writer.h"
#include "mjm_shapes.h"
#include "mjm_color_generator.h"
#include "mjm_rule_list.h"

//#include "mjm_collections.h"

#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>

/*


2018-10-18 copy from mjm_canned_methods
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 
/*

*/


////////////////////////////////////////////////////////////////


namespace dog_plots_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ostream Os;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
//typedef mjm_string_ordering Ordering;

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }
static D myatof(const StrTy & s ) { return atof(s.c_str()); } 
static bool mybit(const IdxTy & flags, const IdxTy & bit, const bool pol=true)
{
const bool b=(((1<<bit)&flags)!=0);
return pol?b:!b;
}
static IdxTy  bad() {  return ~0; } 
//typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::OsTy Os;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // _traits
///////////////////////////////////////////////////////////////


class mjm_dog_plots 
{
typedef  dog_plots_traits::Tr  Tr;
typedef mjm_dog_plots Myt;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::OsTy Os;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock MyBlock;
//typedef Tr::Ordering Ordering;
typedef mjm_ragged_table Ragged;
typedef std::vector<StrTy> Strings;
typedef std::vector<D> Numbers;
typedef std::map< IdxTy , std::map< IdxTy, D > > SparseData;

typedef mjm_color_generator<Tr> Colors;

typedef std::map<StrTy, StrTy> AttMapBase;

class fudd: public AttMapBase
{
typedef AttMapBase Super;
public:
fudd():  Super(),m_def("") {} 
template <class Tfudd>
StrTy & operator[](const Tfudd& c) { return Super::operator[](StrTy(c)); }
template <class Tfudd>
const StrTy &   operator[](const Tfudd& c)const
 {
const StrTy k=StrTy(c);
auto ii=find(k);
if (ii==end()) return m_def; 
  return (*ii).second;

 }

StrTy m_def;
};

typedef fudd AttMap;

//typedef std::map<const char * , StrTy> AttMap;

typedef mjm_svg_writer Sw;

class group_inverse
{
public:
group_inverse(): m_j(~0), m_seq(~0) {}
void add_location(const IdxTy j) { m_j=j; } 
void set_seq(const IdxTy j) { m_seq=j; } 
IdxTy index() const { return m_seq;}
IdxTy j() const { return m_j;}

IdxTy m_j, m_seq;
} ; // group_inverse

typedef group_inverse Jinv;

typedef std::map<StrTy, Jinv> JinvMap;

public:

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }
static D myatof(const StrTy & s ) { return atof(s.c_str()); } 
static bool mybit(const IdxTy & flags, const IdxTy & bit, const bool pol=true)
{ return Tr::mybit(flags,bit,pol); } 
static IdxTy  bad() {  return ~0; } 

class path_param
{
public:
path_param() { clear(); }
void set(const StrTy & color, const D & thick, const StrTy & sym, const IdxTy lbl)
{
//m_color=color;
m_attmap["color"]=color;
m_thick=thick;
m_opaq=1;
m_symbol_size=thick;
//m_symbol=sym;
m_attmap["symbol"]=sym;
m_label_segs=lbl;
}
void opaq(const D & o) { m_opaq=o; } 
// void set_color(const StrTy & s) { 
//m_attmap["color"]=s;
//m_color=s; 
//} 
// void set_symbol(const StrTy & s)
//{ 
//m_attmap["symbol"]=s;
//m_symbol=s;
//BuildSymbol();
//}
void set_fudder(const StrTy & nm, const StrTy & s)
{
m_attmap[nm]=s;
//if ( s=="symbol") m_symbol=s;

}
bool plotted() const { return m_plotted; } 
bool have_polygon() const { return m_x.size()!=0; } 
bool have_lines() const { return m_lines_x.size()!=0; } 
bool custom() const { return have_polygon()||have_lines(); } 
 void set_symbol_size(const StrTy & s){ m_symbol_size=atof(s.c_str());
//BuildSymbol();
 }
// TODO why the FUDD does this compile? lol 
 //void set_y_func(const StrTy & s){ m_y_func=atof(s.c_str());
 void set_y_func(const StrTy & s){ m_y_func=(s.c_str());
 }

//template <class Tv> void set_dashes(const Tv& dashes){ m_stroke_array=dashes; }
 void set_dashes(const std::vector<IdxTy> & dashes){ m_stroke_array=dashes; }
void clear()
{

set("red",1,"",1);
m_clear_x=0;
m_clear_y=0;
m_stroke_array.clear();
 m_traverse_order="2";
m_use_ps=!false; // TODO WTF  rescaling fudd 
m_plotted=false;
}
bool use_ps() const { return m_use_ps; } 
void begin_plotting()
{
BuildSymbol();
BeginPlotting();
m_plotted=true;
}

template <class Trule> 
void setup(Trule & rl, const StrTy & dog, const StrTy & food)
{ // xxxxx
m_name=dog+StrTy("-")+food;
rl.apply(m_attmap["color"], "color",name());
rl.apply(m_thick, "thick",name());
rl.apply(m_opaq, "opaq",name());
rl.apply(m_attmap["symbol"], "symbol",name());
rl.apply(m_symbol_size, "symbolsz",name());
rl.apply(m_label_segs, "labels",name());
// eventually the api will append now assign 
if (rl.apply(m_stroke_array, "dashes",name()))
{
m_stroke_array.clear();
rl.apply(m_stroke_array, "dashes",name());
}

const StrTy & attcol=m_attmap["color"];
MM_ERR(" setting color to "<<MMPR4(dog,food,attcol,m_name))
m_attmap["seg_label"]=m_name; // dog+StrTy("-")+food;
//rl.apply(m_label_segs, "labels",name());
m_seg_label_color=color();
m_seg_label_size=thick();

}

void setup(ReadWriteMap & rwm, const StrTy & base, const StrTy & dog, const StrTy & food)
{ // xxxxx
m_name=dog+StrTy("-")+food;
rwm.get(base+StrTy("color"),m_attmap["color"]);
rwm.get(base+StrTy("thick"),m_thick);
rwm.get(base+StrTy("opaq"),m_opaq);
//rwm.get(base+StrTy("symbol"),m_symbol);
rwm.get(base+StrTy("symbol"),m_attmap["symbol"]);
rwm.get(base+StrTy("symbolsz"),m_symbol_size);
rwm.get(base+StrTy("labels"),m_label_segs);
const StrTy & attcol=m_attmap["color"];
MM_ERR(" setting color to "<<MMPR4(dog,food,attcol,base))

//pp.setup(m_dl.m_rwm,base,dog,food);
//m_seg_label=dog+StrTy("-")+food;
m_attmap["seg_label"]=m_name; // dog+StrTy("-")+food;
m_seg_label_color=color();
m_seg_label_size=thick();


}
//const StrTy & color() const { return m_color;}
const StrTy & color() const { return m_attmap["color"];}
const D & thick() const { return m_thick; } 
const D & opaq() const { return m_opaq; } 
const StrTy&  dog()const  { return m_dog; } 
const StrTy&  dog(const StrTy & x) {m_dog=x;  return m_dog; } 
const StrTy&  food() const  { return m_food; } 
const StrTy&  food(const StrTy & x) {m_food=x;  return m_food; } 
const StrTy&  name() const  { return m_name; } 
const StrTy&  name(const StrTy & x) {m_name=x;  return m_name; } 


StrTy dump()
{
Ss ss;
//ss<<MMPR4(m_color,m_thick,m_symbol,m_label_segs);
const StrTy & symbol=m_attmap["symbol"];
const StrTy & color=m_attmap["color"];
//ss<<MMPR4(m_color,m_thick,symbol,m_label_segs);
ss<<MMPR4(color,m_thick,symbol,m_label_segs);
ss<<MMPR(m_opaq);
return ss.str();
}
typedef mjm_shape_generator ShapeGen;
typedef CommandInterpretter::Words Cw;
void BeginPlotting()
{

if ( m_y_func.length()==0) return; 
CommandInterpretter li;
li.set_split(1,' ');
MM_ERR("BeginPlot " << MMPR2(m_y_func.length(),m_y_func))
li.nextline(m_y_func.c_str());
const Cw& line=li.words();
 const IdxTy len=line.size();
MM_ERR("BeginPlot " << MMPR(len))
const StrTy& sym=line[0];
MM_ERR("BeginPlot " << MMPR(sym))
if ( sym=="pp") { m_use_ps=true; }

}
void BuildSymbol()
{
m_x.clear();
m_y.clear();
m_lines_x.clear();
m_lines_y.clear();
StrTy & s=m_attmap["symbol"];
//if ( m_symbol.length()==0) return; 
if ( s.length()==0) return; 
CommandInterpretter li;
li.set_split(1,' ');
//li.nextline(m_symbol.c_str());
li.nextline(s.c_str());
const Cw& line=li.words();
 const IdxTy len=line.size();
const StrTy& sym=line[0];
if ( sym=="triangle") { ShapeGen::ngon(m_x,m_y,m_symbol_size,3,0,0,0,1); }
else if ( sym=="lines")
{
for(IdxTy i=1; i<(len-1); i+=2)
{
m_lines_x.push_back(atof(line[i].c_str()));
m_lines_y.push_back(atof(line[i+1].c_str()));
}

}
else if ( sym=="ngon")
{
D sz=IfSet( m_symbol_size, line, 1, len);
IdxTy n=IfSet( 3, line, 2, len);
D phi =IfSet( 0, line, 3, len);
IdxTy m=IfSet( 1, line, 4, len);
if (len>1) { D szx=atof(line[1].c_str()); if (szx>0) sz=szx; } 
 ShapeGen::ngon(m_x,m_y,sz,n,0,0,phi,m); 
}
else if ( sym=="square")
{
//static void ngon( Tv & x, Tv & y, const D & r, const IdxTy & n
//, const D & xz=0, const D & yz=0, const D & phiz=0, const IdxTy & m=1)
//StrTy polygon_text( const std::vector<D>  x, const std::vector<D> &  y
//,const StrTy & stroke
//, const StrTy & fill="#00ffffff",const D & sw=1) const
ShapeGen::ngon(m_x,m_y,m_symbol_size,4);
}
} // BuildSymbol 
D  IfSet(  const D & def, const Cw & line, const IdxTy loc, const IdxTy len)
{
D d=def;
if (loc>=len) return d;
 D szx=atof(line[loc].c_str()); 
if (szx>0) d=szx;  
return d; 
}

//StrTy m_color;
D m_thick,m_opaq;
//StrTy m_symbol;
D m_symbol_size;
IdxTy m_label_segs;
//StrTy m_seg_label;
StrTy m_seg_label_color;
D m_seg_label_size;
D m_clear_x,m_clear_y;
std::vector<IdxTy> m_stroke_array;
std::vector<D> m_x,m_y,m_lines_x,m_lines_y;
StrTy m_traverse_order,m_y_func;
AttMap m_attmap;
bool m_use_ps;
bool m_plotted;
StrTy m_dog,m_food,m_name;
}; // path_param

class plot_box
{
typedef std::vector<D> Pos;
typedef std::vector<StrTy> Names;
mutable mjm_calendar m_cal;
public:
// TODO needs clipping etc. 
D utodx(const D & xin) const { return m_lm+m_xpitch*xin; }
D xend() const { return m_xend; }
D utodx(const StrTy & s, const D & xbase) const 
{
//const char * c=s.c_str();
D x=myatoi(s);
if ( x>10000) return utodx(x-xbase);
//x=D(myatoi(m_cal.xlate(s).c_str()))-xbase;
x=D((m_cal.date_to_number(s)))-xbase;
return utodx(x);

}

D utody(const D & yin) const { return m_tm+m_ypitch*(1.0-yin); }
D legx(const D &  ip)const { return  m_lm+ip*m_xpitch+m_szleg;}
D legy(const D &  ip)const { return  m_yleg;}
const D & min_x_space() const  { return m_min_x_space;}


typedef std::vector<D> LegSt; // xv,yv;
void legend_strokes(LegSt & xv, LegSt & yv, const D & _x, const D & _y, const D & wr, const D & hr
,  const D& th,  const D& len )  
//, const  Names & labels)
{
//std::vector<D> xv,yv;
xv.push_back(_x);
yv.push_back(_y+hr+th);
xv.push_back(_x+wr+len*m_leg_sz);
yv.push_back(_y+hr+th);
}



//m_pb.legend_box(x,y,wr,hr,m_pb.m_leg_sz,w,h,m_pb.m_leg_x, m_pb.m_leg_y,m_names);
void legend_box(Pos & _x, Pos & _y, D & wr, D & hr
,  D& w,  D& h 
, const  Names & labels)
{
legend_box(_x,_y,wr,hr,m_leg_sz,w,h,m_leg_x,m_leg_y,labels);
}

void legend_box(Pos & _x, Pos & _y, D & wr, D & hr,const D&  sz
, const D& w, const D& h , const D & xz, const D & yz
, const  Names & labels)
{

D pn=0;
const IdxTy nlabels=labels.size();
for (IdxTy j=0; j<nlabels; ++j)
{
const StrTy & text=labels[j];
 const D p=text.length()*sz;
if  ( p>pn) pn=p;
}
pn+=10;
D x=0;
D y=0;
D yinc=sz+2;
hr=sz;
wr=sz;
for (IdxTy j=0; j<nlabels; ++j)
{
const StrTy & text=labels[j];
if (x>=w) { x=0; y+=yinc; }
_x.push_back(x+xz);
_y.push_back(y+yz);
if ( pn==0)  x+=wr+10+text.length()*sz;
else  x+=wr+pn; // 
} // j


}

void limits(const D & maxx, const D & minn)
{
m_current_max=maxx;
m_current_min=minn;
Compute();

}
StrTy y_label(const D & f)
{
Ss ss;
ss<<std::setprecision(3);
// TODO fix this division crap 
//D q= m_scale*(f-m_current_min); // /(maxx-minn); // m_dl.descale(f); // ((m_dl.m_vmin+(m_dl.m_vmax-m_dl.m_vmin)*f));
//D descale(const D & f){ return  ((m_vmin+(m_vmax-m_vmin)*f)); } 
D q=  ((m_current_min+(m_scale)*f));  
//if ((q>=1e3)&&(q<1e5)) ss<<std::fixed<<std::setprecition(0);
if ((q>=1e3)&&(q<1e5)) ss<<int(q);
else ss<<q;
const StrTy & text=ss.str();
return text;
}


void setup_geometry(ReadWriteMap&  rwm, const D & range)
{

m_title="Set a title" ; rwm.get("title",m_title);

m_tm=50+30; rwm.get("tm",m_tm);
m_bm=330;rwm.get("bm",m_bm);
m_lm=200;rwm.get("lm",m_lm);
m_rm=50;rwm.get("rm",m_rm);

m_xpitch=30;rwm.get("xpitch",m_xpitch);
m_histo_fill=1;rwm.get("histo_fill",m_histo_fill);
// wtf this period is not consistent with the effective one, need
// to make organized precedence ... 
MM_ONCE(" periods dont match",)
IdxTy period=1; rwm.get("period",period);
m_histo_fill=m_histo_fill*D(period);
m_ysz=600;rwm.get("ysz",m_ysz);
m_ypitch=m_ysz;
m_xend=m_xpitch*range+m_lm;
m_xsz=m_xpitch*range+1;rwm.get("xsz",m_xsz);
m_leg_x_frac=1; rwm.get("leg_x_frac",m_leg_x_frac);
m_xs=IdxTy( m_xsz+1)+m_lm+m_rm;
rwm.get("xs",m_xs);
m_ys=m_ysz+m_tm+m_bm;
rwm.get("ys",m_ys);

 m_szleg=m_xpitch; rwm.get("szleg",m_szleg);
 m_yleg=m_tm+m_ysz+2*m_szleg; rwm.get("yleg",m_yleg);
 m_colorleg="red"; rwm.get("colorleg",m_colorleg);
m_min_x_space=1.4*m_szleg; rwm.get("min_x_space",m_min_x_space);

m_xleg=.9*m_lm; rwm.get("xleg",m_xleg);

//m_ylab= StrTy(" Range limited log ratio to mean");
m_ylab= StrTy("Normalized amounts "); rwm.get("ylab",m_ylab);
m_ylab_sz=.2*m_lm; rwm.get("ylab_sz",m_ylab_sz);
m_ylab_x=.2*m_lm; rwm.get("ylab_x",m_ylab_x);
m_ylab_y=.5*m_ysz+m_tm+m_ylab_sz*m_ylab.length()*.5*.5; rwm.get("ylab_y",m_ylab_y);

m_xz=m_lm;
//m_xf=m_lm+m_cats*m_xpitch;
m_yz=m_tm;
m_yf=m_tm+m_ysz;

m_date_y_loc=m_ys-m_bm-20;
D delyd=0;
rwm.get("deldatey",delyd);
m_date_y_loc+=delyd;
rwm.get("datey",m_date_y_loc);
MM_ERR(" ASSFUDD " <<MMPR4(m_date_y_loc,delyd,m_ys,m_bm))


m_leg_x=m_lm; rwm.get("leg_x",m_leg_x);
m_leg_y=.15*m_tm+10; rwm.get("leg_y",m_leg_y);
m_leg_sz=10; rwm.get("leg_sz",m_leg_sz);
m_xf=m_lm+range*m_xpitch;
//MM_ERR(dump())
m_x_rule=0; rwm.get("x_rule",m_x_rule);
m_current_max=1;
m_current_min=0;
Compute();
}

StrTy dump()
{
Ss ss;
ss<<MMPR(m_title);
ss<<MMPR(m_tm);
ss<<MMPR(m_bm);
ss<<MMPR(m_lm);
ss<<MMPR(m_rm);
ss<<MMPR(m_xpitch);
ss<<MMPR(m_histo_fill);
ss<<MMPR(m_xsz);
ss<<MMPR(m_leg_x_frac);
ss<<MMPR(m_ypitch);
ss<<MMPR(m_ysz);
ss<<MMPR(m_date_y_loc);
ss<<MMPR(m_xs);
ss<<MMPR(m_ys);
ss<<MMPR(m_szleg);
ss<<MMPR(m_yleg);
ss<<MMPR(m_min_x_space);
ss<<MMPR(m_xleg);
ss<<MMPR(m_colorleg);
ss<<MMPR(m_ylab);
ss<<MMPR(m_ylab_x);
ss<<MMPR(m_ylab_y);
ss<<MMPR(m_ylab_sz);
ss<<MMPR(m_leg_x);
ss<<MMPR(m_leg_y);
ss<<MMPR(m_leg_sz);
ss<<MMPR(m_xz);
ss<<MMPR(m_xf);
ss<<MMPR(m_yf);
ss<<MMPR(m_yz);
ss<<MMPR(m_x_rule);
ss<<MMPR(m_xend);
ss<<MMPR(m_current_max);
ss<<MMPR(m_current_min);

return ss.str();
}

void Compute()
{
m_scale=(m_current_max-m_current_min);
}
StrTy m_title;
IdxTy m_tm,m_bm, m_lm,m_rm;
D m_xpitch,m_histo_fill, m_xsz,m_leg_x_frac;
D m_ypitch, m_ysz,m_date_y_loc;
D m_xs, m_ys;
D m_szleg,m_yleg,m_min_x_space,m_xleg;
StrTy m_colorleg;
StrTy m_ylab;
D  m_ylab_x, m_ylab_y, m_ylab_sz;
D m_leg_x, m_leg_y,m_leg_sz;
D m_xz,m_xf,m_yf,m_yz;
D m_x_rule;
D m_xend;

D m_current_max, m_current_min, m_scale;
}; // plot_box

class data_loader
{
public:
class path_stats
{
public:
path_stats(): m_max(-1e300), m_min(1e300),m_scale(1) {}
void point(const D & x, const D & y) 
{
if (y>m_max) m_max=y;
if (y<m_min) m_min=y;
}
void point(const D & x, const D & y,const IdxTy serial_end, const StrTy & date_end) 
{
if (y>m_max) m_max=y;
if (y<m_min) m_min=y;
//MM_ERR(" ASS FUDD "<<MMPR4(x,y,serial_end,date_end))
m_serial_end_map[x]=serial_end;
m_date_end_map[x]=date_end;
}
StrTy dump() const
{
Ss ss;
ss<<MMPR4(m_max,m_min,m_serial_end_map.size(),m_date_end_map.size());
return ss.str();
}
void begin_plotting(const bool minzed=false)
{
if (minzed) m_min=0;
D d=m_max-m_min;
if (d<=0){
MM_ERR(" max and min wrong "<<MMPR2(m_max,m_min))
 d=1;
}

m_scale=1.0/d;

}
D rescale( const D & v)
{
return (v-m_min)*m_scale;
}
D scale_value() const { return m_scale;} 
std::map<D,IdxTy> m_serial_end_map;
std::map<D,StrTy> m_date_end_map;
std::map<D,StrTy> m_errors_etc;
D m_max,m_min,m_scale;
}; // path_stats ?

public:
typedef std::vector<StrTy> PlotKey;
typedef std::vector<IdxTy> Dashes;
//typedef std::vector<D> PathStats;
typedef path_stats PathStats;
typedef PathStats Ps;
typedef std::map<PlotKey,Dashes> DashMap;
typedef std::map<PlotKey,StrTy> SymbolMap;
typedef std::map<PlotKey,D> ParamMap;
typedef std::map<PlotKey,PathStats> PathStatsMap;
// for allowing descriptions later... 
typedef std::map<PlotKey,StrTy> SymbolSizeMap;

typedef std::map<StrTy, SymbolMap> MapCollection;
typedef std::map<StrTy, SymbolSizeMap> SizeMapCollection;

template<class Td> 
void concat(Td & d, const Ragged::Line& line, const IdxTy n,const IdxTy m )
{
const IdxTy sz=line.size();
if (sz<(n+1)) return; 
PlotKey pk(n-m);
for(IdxTy i=m; i<n; ++i) pk[i-m]=line[i];
//StrTy 
typedef typename Td::mapped_type Tm;
Tm s;
s=line[n];
for(IdxTy j=(n+1); j<sz; ++j){
s=s+StrTy(" ")+(line[j]);
MM_ERR(" adding   "<<MMPR4(line[0],line[j],pk[0],pk[1]))
}
d[pk]=s;
} // dashes


template<class Td> 
void rconcat(Td & d, const Ragged::Line& line, const IdxTy n,const IdxTy m )
{
const IdxTy sz=line.size();
if (sz<(n+1)) return; 
//StrTy 
typedef typename Td::mapped_type Tm;
Tm s=line[m];
PlotKey pk(n);
IdxTy i=m+1;
while ((i+n)<=sz)
{
for(IdxTy iz=0; iz<n; ++iz) pk[iz]=line[i+iz];
MM_ERR(" adding   "<<MMPR4(s,line[0],pk[0],pk[1]))
d[pk]=s;
i+=n;
}
//for(IdxTy i=m; i<n; ++i) pk[i-m]=line[i];
//for(IdxTy j=(n+1); j<sz; ++j){
//s=s+StrTy(" ")+(line[j]);
//}
//d[pk]=s;
MM_ERR(" fudd added   "<<MMPR4(s,line[0],pk[0],pk[1]))


} // dashes




void setup(Ragged & r, Ragged & pr)
{

MM_ERR("setups day plot ")
ReadWriteMap&  rwm = m_rwm;
pr.to_map(rwm);
setup(r,rwm,pr);
}

void setup(Ragged & r, ReadWriteMap & rwm, Ragged & pr)
{
m_rwm=rwm; // if called from above, this is stupid 
m_bin_idx=0;rwm.get("bin_idx",m_bin_idx);
// TODO the serial index has to match the date used in the liner entry
// or else it will create offset date labels in the plot 
m_serial_idx=1;rwm.get("serial_idx",m_serial_idx);
m_serial_end_idx=2;rwm.get("serial_end_idx",m_serial_end_idx);
m_date_idx=3;rwm.get("date_idx",m_date_idx);
m_date_end_idx=4;rwm.get("date_end_idx",m_date_end_idx);
// kluge moving from emachines beaver to here something changed... 
m_sample_idx=6+0;rwm.get("sample_idx",m_sample_idx);
m_food_idx=7+0;rwm.get("food_idx",m_food_idx);
// eithe pick a stat field number or 
m_value_idx=(1+16);rwm.get("value_idx",m_value_idx);
// use a name and the value is in the following column. 
m_value_idx_name="" ;rwm.get("value_idx_name",m_value_idx_name);

//m_plot_thick=5;rwm.get("plot_thick",m_plot_thick);
m_dates_label=" All dates approximate";
//rwm.get("dates_label",m_dates_label);

const IdxTy prsz=pr.size();
for(IdxTy i=0; i<prsz; ++i)
{
const Ragged::Line& line=pr.line(i);
const IdxTy szl=line.size();
if (szl<1) continue;
const StrTy & cmd=line[0];
MM_ERR(MMPR2(cmd,szl))
//MM_ERR(MMPR(cmd))
if (cmd=="foods") {
for(IdxTy j=1; j<szl; ++j){
 ++m_foods_to_plot[line[j]]; 
MM_ERR(" adding food "<<MMPR2(line[j],m_foods_to_plot.size()))
}
}
if (cmd=="nfoods") {
for(IdxTy j=1; j<szl; ++j){
//  ++m_foods_to_plot[line[j]]; 
const StrTy & lx=line[j];
auto ii=m_foods_to_plot.find(lx);
if (ii!=m_foods_to_plot.end())
{
m_foods_to_plot.erase(ii);
MM_ERR(" removing food "<<MMPR2(lx,m_foods_to_plot.size()))
}
else if ( lx=="*")
{
m_foods_to_plot.clear();
MM_ERR(" removing all foods "<<MMPR2(lx,m_foods_to_plot.size()))
} // star
else MM_ERR(" removing nothing from  foods "<<MMPR2(lx,m_foods_to_plot.size()))

} // j 
} // nfoods


if (cmd=="normalization") {
MM_ERR(" normalization  "<<MMPR(line[1])) 
concat(m_norms_map,line,3,1); } 

if (cmd=="dashes") {
if (szl<3) continue;
const StrTy& dog=line[1];
const StrTy& food=line[2];
//typedef std::vector<StrTy> PlotKey;
//typedef std::vector<IdxTy> Dashes;
//typedef std::map<PlotKey,Dashes> DashMap;
PlotKey pk(2);
pk[0]=dog;
pk[1]=food;
Dashes da(szl-3);
for(IdxTy j=3; j<szl; ++j){
da[j-3]=myatoi(line[j]);
MM_ERR(" adding dashes "<<MMPR3(line[j],food,dog))
}
m_dash_map[pk]=da;
} // dashes
if (cmd=="rsymbol") {
MM_ERR(" rsumbol "<<MMPR(line[1])) 
//rconcat(m_symbol_map,line,2,1); } 
rconcat(m_maps["symbol"],line,2,1); } 
if (cmd=="rcolor") {
MM_ERR(" rcolor "<<MMPR(line[1])) 
rconcat(m_color_map,line,2,1); } 


if (cmd=="rsymbolsz") {
MM_ERR(" rsumbolsz "<<MMPR(line[1])) 
rconcat(m_symbol_size_map,line,2,1); } 

if (cmd=="ryfunc") {
MM_ERR(" ryfunc "<<MMPR(line[1])) 
rconcat(m_yfunc_map,line,2,1); } 


if (cmd=="symbol") {
//if (szl<3) continue;
MM_ERR(" sumbol "<<MMPR(line[1])) 
//concat(m_symbol_map,line,3,1);
concat(m_maps["symbol"],line,3,1);
/* const StrTy& dog=line[1];
const StrTy& food=line[2];
PlotKey pk(2);
pk[0]=dog;
pk[1]=food;
StrTy s;
if (szl>3) s=line[3];
for(IdxTy j=4; j<szl; ++j){
s=s+StrTy(" ")+(line[j]);
MM_ERR(" adding symbol  "<<MMPR3(line[j],food,dog))
}
m_symbol_map[pk]=s;
*/
} // dashes





} // i<prsz
m_day_max=0;
m_day_min=0;
const bool all_dogs=(m_dogs_to_plot.size()==0);
const bool all_foods=(m_foods_to_plot.size()==0);
const IdxTy rsz=r.size();
if (rsz==0)
{
MM_ERR("********************************** no data pounts ")
}
for(IdxTy i=0; i<rsz; ++i)
{
const Ragged::Line& line=r.line(i);
const IdxTy szl=line.size();
if (false) { Ss ss; for(IdxTy k=0; k<szl; ++k) { ss<<MMPR2(k,line[k]); } MM_ERR(ss.str()) }
const StrTy& bin=(line[m_bin_idx]);
const IdxTy serial=myatoi(line[m_serial_idx]);
const StrTy& date=(line[m_date_idx]);
// TODO FIXME thre is a spurious space here moving some fields 
const StrTy& sample=(line[m_sample_idx]);
const StrTy & dog=sample;
const StrTy& food=(line[m_food_idx]);
MM_ONCE(MMPR4(bin,serial,date,sample)<<MMPR(food),)
if (m_day_max<serial) { m_day_max=serial; if (m_day_min==0) m_day_min=serial; }
if (m_day_min>serial) m_day_min=serial;
m_date_to_serial[date]=serial;
m_serial_to_date[serial]=date;
if (all_foods) ++m_foods_to_plot[food];
if (all_dogs) ++m_dogs_to_plot[dog];
if (m_foods_to_plot.find(food)==m_foods_to_plot.end())
{
MM_ERR(" ignoring food "<<MMPR(food))
{ Ss ss; for(IdxTy k=0; k<szl; ++k) { ss<<MMPR2(k,line[k]); } MM_ERR(ss.str()) }

 continue;
}
//if (m_dogs_to_plot.find(dog)==m_dogs_to_plot.end()) continue;
const D value=find_value(line,m_value_idx,m_value_idx_name); // atof(line[m_value_idx].c_str());
//const D value=atof(line[m_value_idx].c_str());
if (i==0) { m_y_max=value; m_y_min=value; }
if (value>m_y_max) m_y_max=value;
if (value<m_y_min) m_y_min=value;
//MM_ERR(MMPR2(serial,date))
++m_samples[sample];
++m_bins[bin];
} 
m_range=m_day_max-m_day_min+1;
MM_ERR(MMPR3(m_range,m_samples.size(),m_bins.size()))
} // setup_day_plot


void load(Ragged & r)
{
//const IdxTy sz=m_samples.size();
const IdxTy rsz=r.size();
//MM_LOOP(ii,m_foods_to_plot)
//{ m_dense_map[(*ii).first].resize(m_range,sz); }
IdxTy i=0;
i=m_day_min;
for(; i<=m_day_max; ++i)  { m_point_order[i]=i-m_day_min;  } 

m_vmin=0;
m_vmax=0;
m_norm=0;
for(IdxTy i=0; i<rsz; ++i)
{

const Ragged::Line& line=r.line(i);
if (false)
{Ss ss;
MM_SZ_LOOP(_i,line,_sz)
{
ss<<line[_i]<<" ";
} // _i
MM_ERR(MMPR(ss.str()))
}

//const IdxTy szl=line.size();
const IdxTy serial=myatoi(line[m_serial_idx]);
const IdxTy serial_end=myatoi(line[m_serial_end_idx]);
//const StrTy& date=(line[m_date_idx]);
const StrTy& sample=(line[m_sample_idx]);
const StrTy & dog=sample;
const StrTy& date=(line[m_date_idx]);
const StrTy& date_end=(line[m_date_end_idx]);
const StrTy& food=(line[m_food_idx]);
if (m_foods_to_plot.find(food)==m_foods_to_plot.end()) continue;
const D _value=find_value(line,m_value_idx,m_value_idx_name); // atof(line[m_value_idx].c_str());
//const D _value=atof(line[m_value_idx].c_str());
StrTy xxx;
std::basic_string<char>& ( std::basic_string<char>::  * _p) (const std::basic_string<char> &);
_p=&StrTy::operator=;

// TODO WTF  does nothing.. 
hierarchy_set(_p , xxx,m_norms_map,dog,food);


//D normxi=atof(m_norms_map[dog][food].c_str());
D normxi=atof(xxx.c_str());
m_norm=normxi;
const D value= (normxi!=0)?(_value/normxi):_value;
MM_ERR(MMPR3(serial,date,sample)<<MMPR4(food,value,_value,normxi))
if (i==0) { m_vmax=value; m_vmin=value; }
if (value>m_vmax) m_vmax=value;
if (value<m_vmin) m_vmin=value;

if (m_sample_order.find(sample)==m_sample_order.end())
{
MM_ERR(" not found sample "<<MMPR4(date,food,sample,m_sample_order.size())) 
}
const IdxTy n=m_sample_order[sample];
IdxTy serial_eff=serial;
if (m_point_order.find(serial_eff)==m_point_order.end())
{
MM_ERR(" mising point order "<<MMPR4(serial,sample,date,food)<<MMPR3(m_day_min,m_day_max,m_point_order.size()))
// TODO NOTE first noted on end missing points, may be a problem here
// too. 

while (serial_eff<=serial_end)
{
if  (m_point_order.find(serial_eff)!=m_point_order.end()) break;
++serial_eff;
} // while
const bool found=(m_point_order.find(serial_eff)!=m_point_order.end());
MM_ERR(" kluged start date adjust "<<MMPR3(serial_eff,serial,found))
}// if

IdxTy serial_end_eff=serial_end;
if (m_point_order.find(serial_end_eff)==m_point_order.end())
{
MM_ERR(" mising point order END  "<<MMPR4(serial_end,sample,date,food)<<MMPR3(m_day_min,m_day_max,m_point_order.size()))
while (serial_end_eff>=serial_eff)
{
if (m_point_order.find(serial_end_eff)!=m_point_order.end()) break;
--serial_end_eff;
if (serial_end_eff>serial_end) break; // roll over 
} // while 

MM_ERR(" Using botched end "<<MMPR4(serial,serial_eff,serial_end,m_point_order[serial_eff])<<MMPR3(serial_end,serial_end_eff,m_point_order[serial_end_eff]))

}// if 

if (0==m_point_order[serial_end_eff])
{
MM_ERR(" probably should not be zero ")
MM_ERR(" Using botched end BAD ZERO MAYBE  "<<MMPR4(serial,serial_eff,serial_end,m_point_order[serial_eff])<<MMPR3(m_point_order[serial_end_eff],serial_end_eff,m_point_order[serial_end_eff]))
}

// this just removes the time or x origin and AFAICT could be a float
//const IdxTy x=m_point_order[serial];
// doh, serial_end may not exist... 
const IdxTy x=(m_point_order[serial_eff]+m_point_order[serial_end_eff])>>1;
// don't user serial_end as will put spurious value if none altready there 
MM_ERR(MMPR(serial_eff)<<MMPR4(x,serial,serial_end,m_point_order[serial_eff])<<MMPR3(m_point_order[serial_end_eff],serial_end_eff,m_point_order[serial_end_eff]))


std::vector<StrTy> pk(2);
pk[0]=dog;
pk[1]=food;
//m_path_stats_map[dog][food].point(x,value);
// FIXME serial_end_eff???
m_path_stats_map[pk].point(x,value,serial_end,date_end);
// this needs to be sparse or find invalid like negative 
//m_dense_map[food](x,n)+=value; // atof(line[m_value_idx]);
// FIXME this is the total over the allowed time period but 
// no pderiod length denominator so I guess ok 
m_sparse_map[food][n][x]+=value;
MM_ERR(MMPR2(n,x))
}
MM_ERR(" load done "<<MMPR4(m_vmin,m_vmax,m_day_min,m_day_max))
} //load 


IdxTy range() const { return m_range; } 

 D rescale(const D & v){ return  (v-m_vmin)/(m_vmax-m_vmin);}
//D q=m_dl.descale(f); // ((m_dl.m_vmin+(m_dl.m_vmax-m_dl.m_vmin)*f));
D descale(const D & f){ return  ((m_vmin+(m_vmax-m_vmin)*f)); } 
template <class Tm > StrTy dump(Tm & m, const StrTy & lbl )
{ Ss ss;
MM_LOOP(ii,m) { ss<<MMPR3(lbl,(*ii).first,(*ii).second); } 
return ss.str();
}
StrTy dump()
{
Ss ss;
ss<<dump(m_dogs_to_plot,"dogs");
ss<<dump(m_foods_to_plot,"foods");
ss<<MMPR2(m_y_max,m_y_min);
ss<<MMPR2(m_day_max,m_day_min);
ss<<MMPR2(m_vmax,m_vmin);
return ss.str();
}

IdxTy update_sample_order()
{
m_sample_order.clear();
m_isample_order.clear();
IdxTy i=0;
MM_LOOP(ii,m_samples) { 
m_sample_order[(*ii).first]=i;
m_isample_order[i]=(*ii).first;

 ++i; } 

return i;
}

Ps & find_path_stats(const StrTy & dog, const StrTy & food)
{
std::vector<StrTy> pk(2); pk[0]=dog; pk[1]=food;
return find_path_stats(pk); 
}

Ps & find_path_stats(std::vector<StrTy> & pk)
{ 
data_loader::PathStatsMap& psm= m_path_stats_map;
auto psi =psm.find(pk); // m_dl.m_path_stats_map.find(pk);
const bool found_ps=(psi!=psm.end()); // m_dl.m_path_stats_map.end());
return found_ps?((*psi).second):m_psdef; // m_dl.m_path_stats_map[pk];
}


IdxTy find_in_sample_order(const StrTy & dog) const
{
auto ii=m_sample_order.find(dog);
if ((ii)==m_sample_order.end())
{ MM_ERR(" not found "<<MMPR3(0,dog, m_sample_order.size())) return bad(); }
IdxTy i=(*ii).second; // m_dl.m_sample_order[(*is).first];
return i; 
}

D find_value(const Ragged::Line  & line,const IdxTy value_idx,const StrTy & value_idx_name) // atof(line[m_value_idx].c_str());
{
if (value_idx!=0) return atof(line[value_idx].c_str());
// NB MINUS ONE doh... 
const IdxTy sz=line.size()-1;
if (sz==~0) return 0;
for (IdxTy i=0; i<sz; ++i)
{

if (line[i]==value_idx_name)  return atof(line[i+1].c_str());
}
return 0;
} // find_value


StrTy  find_string(const Ragged::Line  & line,const StrTy & value_idx_name) // atof(line[m_value_idx].c_str());
{
// NB MINUS ONE doh... 
const IdxTy sz=line.size()-1;
if (sz==~0) return "" ;
for (IdxTy i=0; i<sz; ++i)
{

if (line[i]==value_idx_name)  return (line[i+1].c_str());
}
return"" ;
} // find_string



ReadWriteMap m_rwm;
IdxTy m_bin_idx,m_serial_idx, m_serial_end_idx, m_date_idx,m_date_end_idx, m_sample_idx, m_food_idx,m_value_idx;
IdxTy m_range;
StrTy m_value_idx_name;
//std::map<StrTy, MyBlock> m_dense_map;
std::map<StrTy, SparseData> m_sparse_map;
StrTy m_dates_label;
D m_y_max,m_y_min;
IdxTy m_day_min,m_day_max;
D m_vmin, m_vmax;
D m_norm;
typedef std::map<StrTy,IdxTy> IdxMap;
typedef std::map<IdxTy,StrTy> StrMap;
IdxMap m_bins;
StrMap m_serial_to_date;
IdxMap m_date_to_serial;
IdxMap m_samples;
IdxMap m_sample_order;
StrMap m_isample_order;
std::map<IdxTy,IdxTy> m_point_order;
IdxMap m_foods_to_plot;
IdxMap m_dogs_to_plot;
//SymbolMap m_symbol_map;
SymbolMap m_color_map;
SymbolMap m_yfunc_map;
MapCollection m_maps;
SymbolSizeMap m_symbol_size_map;
DashMap m_dash_map;
// TODO FIXME this is a String value for easy of use and flexibility
// but repated use may get slow. 
SymbolMap m_norms_map;

PathStatsMap m_path_stats_map;
Ps m_psdef;
}; // data_loader


class data_rounding
{
public:

//TODO FIXME want to force a grid line at zero sometimes 
static IdxTy  adjust_grid(D & mx, D & mn, const IdxTy  & ygrid)
//static D  adjust_grid(D & mx, D & mn, const IdxTy  & ygrid)
{
const D vmaxi=mx; // m_dl.m_vmax;
const D vmini=mn; // m_dl.m_vmin;
const IdxTy ygridi=ygrid;
D e= 1.1;  
exp_range(mx,true,e);
exp_range(mn,!true,e);
D del=1;
IdxTy ibest=ygridi;
for(IdxTy i=5; i<16; ++i)
{

//D diff=m_dl.m_vmax-m_dl.m_vmin;
D diff=mx-mn;
D idelg=D(i-0)/diff;
idelg=idelg-IdxTy(idelg);
if (idelg<del) {del=idelg; ibest=i; }
}
D ygrido=ibest;
if ( ibest==ygridi) ygrido=5; // 2019-05
//MM_ERR(MMPR4(m_dl.m_vmax,vmaxi,m_dl.m_vmin,vmini)<<MMPR2(m_ygrid,ygridi))
MM_ERR(MMPR4(mx,vmaxi,mn,vmini)<<MMPR2(ygrido,ygridi))
return ygrido;
}

static void  exp_range( D & x, const bool top, const D&  e )
{
const bool f=(top&&(x>0))||!(top&&(x<0));
D y=x*((f)?e:(1.0/e));
Ss ss;
ss<<std::setprecision(2)<<y;
ss>>x;
}

static D  trunc( const D & x, const IdxTy n  )
{
Ss ss;
ss<<std::setprecision(n)<<x;
D y;
ss>>y;
return y;
}

}; // data_rounding 
/*
void default_rainbod(Strings & rainbow)
{
rainbow.push_back("red");
rainbow.push_back("blue");
rainbow.push_back("green");
rainbow.push_back("orange");
rainbow.push_back("black");

}
*/

void colors_and_day_names(ReadWriteMap & rwm)
{
IdxTy i=m_dl.update_sample_order();

MM_ERR(MMPR3(i,m_dl.m_sample_order.size(),m_xlabels.size()))
//m_names.resize(i);
//m_colors.resize(i);
m_xlabels.resize(m_dl.range());
MM_SZ_LOOP(i,m_xlabels,xxx)
{
Ss ss;
StrTy xlabel="xlabel";
ss<<m_dl.m_serial_to_date[(i+m_dl.m_day_min)];
StrTy & n=m_xlabels[i];
if (n.length()==0) xlabel=ss.str();
else xlabel=n;
//MM_ERR(MMPR2(rwm.size(),name))
rwm.get(xlabel,xlabel);
m_xlabels[i]=xlabel;
} // i 
MM_LOOP(ii,m_dl.m_samples)
{
const StrTy & nameauto=(*ii).first;
IdxTy i=m_dl.m_sample_order[nameauto];
#if 0 
{
Ss sc;
StrTy col="color";
sc<<col<<i;
// TODO m_default_color is blank at this point 
StrTy c=m_default_color; // rainbow[i%rainbow.size()];
//StrTy c= rainbow[i%rainbow.size()];
rwm.get(sc.str(),c);
MM_ERR(" setting color "<<MMPR2(i,c))
m_colors[i]=c;
}

Ss sc;
m_names[i]=nameauto;
StrTy col="name";
sc<<col<<":"<<nameauto;
rwm.get(sc.str(),m_names[i]);

#endif

++i;
} // ii  samples 

} // colors_and_day_names


// xy plots with headers and limits derived from data scan
// 5855 17567 17567 2018-02-05 1  Greta B-6  m_n=1 m_ntot=1 m_tot=12.5 m_max=12.5 m_avg=12.5 m_avg_dose=0 m_freq=1 m_invalid=0
// ./mjm_linc_graph.out -cmd "read-ragged in xxx2a" -cmd "add-ragged p period 3" -cmd "read-ragged p lys.txt"  -cmd "snacks-txt-svg srxxx 1 in p" -quit

// ss<<dayfirst<<" "<<daylast<<" "<<lbl<<" "<<pereff<<" ";
// StrTy  binlab=ss.str();
// os<<bin<<sep<<binlab<<sep<<dog<<sep<<food<<sep<<ss.str()<<CRLF;
//StrTy lbl=sf.serial_to_date[daylast]



typedef data_loader::path_stats Ps;

void write_svg_day_best_labels(OsTy & os, Sw& sw, const IdxTy flags)
{
//const IdxTy sz=m_dl.range();
auto & sm=m_dl.m_sparse_map;
MM_LOOP(ii,sm)
{
const StrTy & food=(*ii).first;
auto & sd=(*ii).second;
MM_LOOP(jj,sd)
{
const StrTy & dog=m_dl.m_isample_order[(*jj).first];
// this is now a sparse datamap 
auto & p=(*jj).second; // x and y values 
MM_LOOP(iip,sm)
{
auto & sdp=(*iip).second;
MM_LOOP(jjp,sdp)
{
auto & pp=(*jjp).second; // x and y values 
MM_ERR(MMPR2(p.size(),pp.size()))
} // jjp 

} // iip

path_param & pp = m_pm[dog][food];
pp.m_clear_x=0;
pp.m_clear_y=0;

} // jj 

} // ii 


}

StrTy dump_dl() { return m_dl.dump(); } 
// pure histo with biggest plotted first 
void write_svg_day_data_histo(OsTy & os, Sw& sw, const IdxTy flags)
{
const bool noverlap=true;
// this is now in display units or pixels I guess 
//const D bin_val=1.0/.02; // resolve vertical equality bin size 
const D bin_val=20; // resolve vertical equality bin size 
Ss ss;
std::map<IdxTy, std::map< D, std::vector< path_param * > > > pm;
std::map<IdxTy, std::map< D, std::vector< D >  > > pmrealy;
MM_LOOP(ff,m_dl.m_foods_to_plot)
{
const StrTy & food=(*ff).first; // 
auto & sd=m_dl.m_sparse_map[food];
MM_LOOP(is,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*is).first;
IdxTy i=m_dl.find_in_sample_order(dog);
auto  ii=sd.find(i);
if (ii==sd.end()) { ++i; continue;}
auto & s=(*ii).second;
//IdxTy i=m_dl.find_in_sample_order(dog);
//auto  ii=sd.find(i);
//if (ii==sd.end()) { ++i; continue;}
path_param & pp = m_pm[dog][food];
pp.dog(dog);
pp.food(food);
Ps & ps=m_dl.find_path_stats( dog, food);
pp.begin_plotting();
ps.begin_plotting(true);
MM_ERR(MMPR2(pp.dump(),ps.dump()))


MM_LOOP(jj,s)
{
const IdxTy  ip=(*jj).first;
D v=(*jj).second; // sdf[ip];
//const StrTy & xtext= m_xlabels[ip]; // line[text_idx];
// this one being called 
//MM_ERR("PLOTDATA2 "<<MMPR4(dog,food,xtext,ip)<<MMPR(v))
if (pp.use_ps()) v=ps.rescale(v);
else v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
//assert_normal(v,sd,ii,jj,m_dl);
//const D y=m_pb.utody(v);
const D yreal=m_pb.utody(v);
D y=yreal;
// TODO when the things are similar but not equal one obscures
// the other, want to split in that case... 
// this is in pixels or a display unit, want a const amt doh 
if ( bin_val!=0)
{
D binx=(y/bin_val)+.5;
D binf=(y*bin_val);
if (binf!=0)
{
D binfx=(y/binx)+.5;
//D binfx=(y/binf)+.5;
D yy=int(binx)*bin_val;
//D yy=int(binfx)*binf;
MM_ERR(MMPR2(binf,binfx)<<MMPR4(binx,yy,y,bin_val))
y=yy;
}
} // bin_val
//plot_seg2(ss,sw,xv,yv,y,iplast,ip,ipnext,s,pp,ps);
pm[ip][y].push_back(&pp);
pmrealy[ip][y].push_back(yreal);
}
/*
const IdxTy sz=m_dl.range();
for (IdxTy ip=0; ip<sz; ++ip)
{
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
auto ff=sdf.find(ip);
const bool found=(ff!=sdf.end());
D v=(!found)?0:(*ff).second;
v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
const D y=m_pb.utody(v);
*/

} // is 

} // ff

const D xpitch=m_pb.m_xpitch; // m_lm+ip*m_catpitch;
const D xfill=xpitch*m_pb.m_histo_fill; // m_lm+ip*m_catpitch;
const D yzed=m_pb.utody(0);
D ylast=yzed;
MM_LOOP(iip,pm)
{
const IdxTy ip=(*iip).first;
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
MM_LOOP(vv,(*iip).second)
{
// these should really be binned better... 
const D y=(*vv).first;
D yf=yzed;
if (noverlap) { 
auto vnext=(vv); ++vnext;
//if ( vnext!=(*iip).second.end()) yf=(*vnext).first;
if ( vnext!=(*iip).second.end()) yf=pmrealy[ip][(*vnext).first][0];
}

// should only have one entry 
const IdxTy dsplit=(*vv).second.size();
//if (dsplit>1) yf=yzed;

//https://stackoverflow.com/questions/13069446/simple-fill-pattern-in-svg-diagonal-hatching

const D effhg= y-yf;
const bool new_way=true;
if ( new_way)
{
std::vector<StrTy>  colors; // (dsplit);
 D opaq=1;
// this has a shadow in pmrealy... 
const auto &  yreal=pmrealy[ip][y];
MM_LOOP(vvv,(*vv).second) { path_param & pp= *(*vvv); colors.push_back(pp.color()); opaq=pp.opaq(); } 
// TODO the non-histo plot is still wrong this uses the start date not end date 
/*

template <class Tv,class Tc> 
StrTy sparse_shade_rect_array_text( const Tv& xv, const Tv& yv, const Tv &xw, const Tv  & xh, const Tc & colorv, const D & opaq=1.0) const
{
ZZ
*/
const bool old_bin=false;
if (! old_bin){
std::vector<D> xv,yv,wv,hv;
D del=xfill/D(colors.size());
for(IdxTy l=0; l<colors.size(); ++l)
{
xv.push_back(x+D(l)*del);
const auto &  yreal=pmrealy[ip][y][l];
yv.push_back(yreal);
wv.push_back((del));
hv.push_back(yreal-yf);

} // l 
os<<sw.sparse_shade_rect_array_text(xv,yv,wv,hv,colors,opaq);

}else{
os<<sw.shade_rect_text(x-0*xfill,y,xfill,effhg,colors,opaq);
} // old_bin
//  try to remove bins but still multuple colors, needs a y vector  
//os<<sw.shade_rect_text(x-0*xfill,yreal[0],xfill,effhg,colors,opaq);
os<<CRLF;
D sz=1.0+.05*xfill;
D yeff=y;
if (yeff<(yzed-sz)) yeff+=sz; 
IdxTy seg=0;
const IdxTy nsegs=(*vv).second.size();
MM_LOOP(vvv,(*vv).second) { path_param & pp= *(*vvv);  
//os<<sw.gtext_text(pp.food(), x+.5*xfill, yeff, sz, "white",  StrTy("middle") ,0);
os<<sw.gtext_text(pp.food(), x+xfill/D(nsegs)*(.5+D(seg)), yeff, sz/D(nsegs), "white",  StrTy("middle") ,0);
os<<CRLF;
++seg;
} // vvv
}
else { // old
MM_LOOP(vvv,(*vv).second)
{
path_param & pp= *(*vvv);
const StrTy & cj=pp.color(); // m_colors[i];
//const StrTy & text=m_names[i];
D effh= y-yf;
//if ( effh<1) { yf=yzed; effh=y-yf; } 
// ASSFUDD 
//os<<sw.shade_rect_text(x-.5*xpitch,yf,xpitch,-y+yf,cj,.5);
//os<<sw.shade_rect_text(x-.5*xfill,y,xfill,y-yf,cj,.5);
// it looks like the dates are the end dates, need to check alignement 
os<<sw.shade_rect_text(x-xfill,y,xfill,effh,cj,.5);
//ss<<sw.gtext_text(symbol, x, y+0*pp.m_symbol_size, pp.m_symbol_size, pp.color(),  StrTy("middle") ,0);
D sz=1.0+.05*xfill;
D yeff=y;
if (yeff<(yzed-sz)) yeff+=sz; 
os<<sw.gtext_text(pp.food(), x-.5*xfill, yeff, sz, "white",  StrTy("middle") ,0);
// seg2
}
} // old
ylast=y;

} // vv
} // iip 


} // write_svg_day_data_histo


// pass first 10 then 101 aka 2 and  5
void write_svg_day_data(OsTy & os, Sw& sw, const IdxTy flags)
{

Ss ss;
MM_LOOP(ff,m_dl.m_foods_to_plot)
{
const StrTy & food=(*ff).first; // 
auto & sd=m_dl.m_sparse_map[food];
MM_ERR(" writing food "<<MMPR3(sd.size(),food, m_dl.m_foods_to_plot.size()))
MM_LOOP(is,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*is).first;
IdxTy i=m_dl.find_in_sample_order(dog);
auto  ii=sd.find(i);
if (ii==sd.end()) { ++i; continue;}
auto & s=(*ii).second;
path_param & pp = m_pm[dog][food];
MM_ERR(" begin  plotting "<<MMPR3(dog,food,pp.color()))
Ps & ps=m_dl.find_path_stats( dog, food);
pp.begin_plotting();
ps.begin_plotting(true);
MM_ERR(MMPR2(pp.dump(),ps.dump()))
const IdxTy traverse_order=myatoi(pp.m_traverse_order);
std::vector<D> xv,yv;
assert_found((*ii).first,sd,dog,food);
auto & sdf=sd[(*ii).first];
if (traverse_order==2)
{
IdxTy iplast=bad();
MM_LOOP(jj,s)
{
const IdxTy  ip=(*jj).first;
auto jnext=jj; ++jnext;
IdxTy ipnext=(jnext==s.end())?bad():(*jnext).second;
assert_found(ip,sdf,dog,food);
D v=sdf[ip];
const StrTy & xtext= m_xlabels[ip]; // line[text_idx];
// this one being called 
MM_ERR("PLOTDATA2 "<<MMPR4(dog,food,xtext,ip)<<MMPR(v))
if (pp.use_ps()) v=ps.rescale(v);
else v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
assert_normal(v,sd,ii,jj,m_dl);
const D y=m_pb.utody(v);
plot_seg2(ss,sw,xv,yv,y,iplast,ip,ipnext,s,pp,ps);

iplast=ip;
} // jj 


} // 2
else if (traverse_order==0)
{
MM_LOOP(jj,s)
{
const D ip=(*jj).first-0; // m_day_min;
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
assert_found(ip,sdf,dog,food);
D v=sdf[ip];
v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
assert_normal(v,sd,ii,jj,m_dl);
const D y=m_pb.utody(v);
MM_ERR("PLOTDATA1 "<<MMPR4(dog,food,ip,v))
D ipnext=ip+1; // m_day_min;
if ( m_histogram_like) { 
auto jjn=jj; ++jjn;
if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
}
plot_seg(ss,sw,xv,yv,x,y,ip,ipnext,s,pp,ps);
//plot_seg(ss,sw,xv,yv,x,y,ip,jj,s,pp);
} // jjj, data points
} // traverser_order==0
else if (traverse_order==1)
{
const IdxTy sz=m_dl.range();
for (IdxTy ip=0; ip<sz; ++ip)
{
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
auto ff=sdf.find(ip);
const bool found=(ff!=sdf.end());
D v=(!found)?0:(*ff).second;
MM_ERR("PLOTDATA0 "<<MMPR4(dog,food,ip,v))
v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
const D y=m_pb.utody(v);
int pse=(int(ps.m_serial_end_map[ip])-int(m_dl.m_day_min))-int(ip);
MM_ERR(MMPR4( pse,ps.m_serial_end_map[ip],m_dl.m_day_min,x)<<MMPR(ip))
if (pse<0) pse=1;
plot_seg(ss,sw,xv,yv,x,y,ip,ip+pse,s,pp,ps);

} // ip 
} // ==1
else MM_ERR(" bad traversing "<<MMPR3(traverse_order,dog,food))
if (xv.size()!=0) { 
const D yz=m_pb.utody(0);
if (m_histogram_like) { xv.push_back(xv.back()); yv.push_back(yz); } 
if (pp.m_stroke_array.size()==0) 
	os<<sw.vector_text(xv,yv,pp.color(),pp.thick(),.7); 
else{
MM_ERR(" using stroke array "<<MMPR(pp.m_stroke_array.size()))
 os<<sw.vector_text(xv,yv,pp.m_stroke_array, pp.color(),pp.thick(),.7); 
}
os<<CRLF;
 }

} // dogs

} // foods
if (m_label_foods) os<<ss.str();

if (false)  write_svg_day_best_labels( os, sw, flags);

} // write_svg_day_data


void write_svg_day_data_old(OsTy & os, Sw& sw, const IdxTy flags)
{

Ss ss;
//const bool label_foods=true;
//const bool histogram_like=!true;
MM_LOOP(ff,m_dl.m_foods_to_plot)
{
const StrTy & food=(*ff).first; // 
auto & sd=m_dl.m_sparse_map[food];
MM_ERR(" writing food "<<MMPR3(sd.size(),food, m_dl.m_foods_to_plot.size()))
MM_LOOP(is,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*is).first;
//IdxTy find_in_sample_order() const
IdxTy i=m_dl.find_in_sample_order(dog);
//if (m_dl.m_sample_order.find((*is).first)==m_dl.m_sample_order.end())
//{ MM_ERR(" not found "<<MMPR3(0,(*is).first, m_dl.m_sample_order.size())) }
//IdxTy i=m_dl.m_sample_order[(*is).first];


auto  ii=sd.find(i);
//if (ii==m_sparse_data.end()) { ++i; continue;}
if (ii==sd.end()) { ++i; continue;}
auto & s=(*ii).second;
MM_ERR(" plotting "<<MMPR2(dog,food))
path_param & pp = m_pm[dog][food];
std::vector<StrTy> pk(2);
pk[0]=dog;
pk[1]=food;
// FIXME this will load up the map with non existing ones 
data_loader::PathStatsMap& psm= m_dl.m_path_stats_map;
auto psi =psm.find(pk); // m_dl.m_path_stats_map.find(pk);
const bool found_ps=(psi!=psm.end()); // m_dl.m_path_stats_map.end());
Ps psdef;
Ps &  ps =found_ps?((*psi).second):psdef; // m_dl.m_path_stats_map[pk];
//data_loader::path_stats&  ps =found_ps?((*psi).second):Ps(); // m_dl.m_path_stats_map[pk];

MM_ERR(MMPR(pp.dump()))
const IdxTy traverse_order=myatoi(pp.m_traverse_order);
std::vector<D> xv,yv;
assert_found((*ii).first,sd,dog,food);
//if ( sd.find((*ii).first)==sd.end())
//{
//MM_ERR(" missing sparse data "<<MMPR2(dog,food))
//}
auto & sdf=sd[(*ii).first];
if (traverse_order==0)
{
MM_LOOP(jj,s)
{
const D ip=(*jj).first-0; // m_day_min;
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
assert_found(ip,sdf,dog,food);

//if ( sdf.find(ip)==sdf.end())
//{ MM_ERR(" missing sparse data "<<MMPR2(dog,food)) }

D v=sdf[ip];
v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
assert_normal(v,sd,ii,jj,m_dl);
//if ((v>1)||(v<0))
//{
//D vfudd=sd[(*ii).first][(*jj).first];
//MM_ERR(" v oor "<<MMPR4((*ii).first, (*jj).first, v,vfudd)<<MMPR2(m_dl.m_vmax,m_dl.m_vmin))
//}
const D y=m_pb.utody(v);


D ipnext=ip+1; // m_day_min;
if ( m_histogram_like) { 
auto jjn=jj; ++jjn;
if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
}
plot_seg(ss,sw,xv,yv,x,y,ip,ipnext,s,pp,ps);
//plot_seg(ss,sw,xv,yv,x,y,ip,jj,s,pp);
} // jjj, data points
} // traverser_order==0
else if (traverse_order==1)
{
const IdxTy sz=m_dl.range();
for (IdxTy ip=0; ip<sz; ++ip)
{
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
auto ff=sdf.find(ip);
const bool found=(ff!=sdf.end());
D v=(!found)?0:(*ff).second;
v=m_dl.rescale(v); // (v-m_dl.m_vmin)/(m_dl.m_vmax-m_dl.m_vmin);
//assert_normal(v,sd,ii,jj,m_dl);
const D y=m_pb.utody(v);
//m_range=m_day_max-m_day_min+1;
int pse=(int(ps.m_serial_end_map[ip])-int(m_dl.m_day_min))-int(ip);
MM_ERR(MMPR4( pse,ps.m_serial_end_map[ip],m_dl.m_day_min,x)<<MMPR(ip))
if (pse<0) pse=1;
plot_seg(ss,sw,xv,yv,x,y,ip,ip+pse,s,pp,ps);

} // ip 
} // ==1
else MM_ERR(" bad traversing "<<MMPR3(traverse_order,dog,food))
if (xv.size()!=0) { 
if (pp.m_stroke_array.size()==0) 
	os<<sw.vector_text(xv,yv,pp.color(),pp.thick(),.7); 
else{
MM_ERR(" using stroke array "<<MMPR(pp.m_stroke_array.size()))
 os<<sw.vector_text(xv,yv,pp.m_stroke_array, pp.color(),pp.thick(),.7); 
}
os<<CRLF;
 }

} // dogs

} // foods
if (m_label_foods) os<<ss.str();

if (false)  write_svg_day_best_labels( os, sw, flags);

} // write_svg_day_data


//assert_normal(v,sd,ii,jj,m_dl);
template <class Ta, class Tb,class Tc, class Td  > 
static bool assert_normal(const D & v,Ta & sd, Tb & ii, Tc & jj,Td & dl)
{
const bool oor= ((v>1)||(v<0));
if (oor)
{
D vfudd=sd[(*ii).first][(*jj).first];
MM_ERR(" v oor "<<MMPR4((*ii).first, (*jj).first, v,vfudd)<<MMPR2(dl.m_vmax,dl.m_vmin))
}
return !oor;
} 


template <class Tc, class Tm,class Ts > 
static bool assert_found( const Tc & key, const Tm & m, const Ts & s1, const Ts & s2 )
{
const bool found = ( m.find(key)!=m.end());
if (!found)
{
MM_ERR(" missing "<<MMPR3(key,s1,s2))
} 
return found; 
}

int xendf(const IdxTy ip, const Ps & ps)
{ 
 auto const  & sem=ps.m_serial_end_map;
auto const  & semii=sem.find(ip);
if (semii==sem.end()) return 0; 
return m_pb.utodx(int((*semii).second)-int(m_dl.m_day_min)); }

template<class Tj, class Ts>
void plot_seg2(OsTy& ss, Sw& sw
, std::vector<D> & xv, std::vector<D> & yv, const D & y
, const IdxTy iplast, const IdxTy ip , const Tj & ipnext, Ts & s 
, const path_param & pp ,const Ps & ps ) 
{
const StrTy & llabel=pp.m_attmap["seg_label"]; // pp.m_seg_label;
//const StrTy & symbol=pp.m_symbol;
const StrTy & symbol=pp.m_attmap["symbol"];
//const D xbar=m_pb.utodx(ipnext); // m_lm+(ipnext)*m_catpitch;
// TODO FIXME danger will robinson this is all floating point eq tests lol 
const D x=m_pb.utodx(ip); // m_lm+ip*m_catpitch;
const D xend=xendf(ip,ps);
const D xendx=xend+m_pb.m_xpitch;
const D xmid=(x+xend+m_pb.m_xpitch )*.5;
//const D xlast=m_pb.utodx(iplast); // m_lm+ip*m_catpitch;
//const D xlastend=xendf(iplast,ps);
const D xlastendx=xendf(iplast,ps)+m_pb.m_xpitch;
//const D xlastmid=(xlast+xlastend+m_pb.m_xpitch )*.5;
//const D xnext=m_pb.utodx(ipnext); // m_lm+ip*m_catpitch;
//const D xnextend=xendf(ipnext,ps);
const bool plot_symbol=(symbol.length()!=0);
const bool have_prior=(xv.size()!=0);
//const bool have_next=(ipnext!=bad());
const D yz=m_pb.utody(0);

if ( m_histogram_like&&!have_prior) { xv.push_back(x); yv.push_back(yz); }

if (!have_prior) { xv.push_back(x); yv.push_back(y);}
else {
D xprior=xv.back();
D yprior=yv.back();
D dx=(x-xprior-m_pb.m_xpitch)/m_pb.m_xpitch;
const bool dxdiff=(dx*dx>1e-2);
const bool gaphisto=((ip-iplast)>1); // ((dx*dx)>1e-2); //(x!=xprior);
//if ((x!=xprior)||(y!=yprior)) { 
if ((dxdiff)||(y!=yprior)) { 
if (m_histogram_like) { 
if (gaphisto) {

MM_ERR("DUXK"<<MMPR4(x,xprior,ip,iplast)<<MMPR3(gaphisto,yz,dx))
xv.push_back(xprior); yv.push_back(yz);
xv.push_back(x); yv.push_back(yz);
}
else 
 {xv.push_back(xprior); yv.push_back(y);}

xv.push_back(x); yv.push_back(y);  
D xendxr=xendx;
if (xendxr>m_pb.xend()) xendxr=m_pb.xend(); 
xv.push_back(xendxr); yv.push_back(y);  
 }
else {  // !histogram_like
// need to force to zero if missing points... 
//if (xprior!=xlastmid) 
if (x!=xlastendx) 
//if (gap) 
{
xv.push_back(xprior); yv.push_back(yz);  

xv.push_back(xmid); yv.push_back(yz);  
}
xv.push_back(xmid); yv.push_back(y);  
} // !histogram_like

MM_ERR(MMPR4(m_label_foods,ip,m_label_mod,m_label_res))
if ((m_label_foods) &&((int(ip)%m_label_mod)==m_label_res)){
if  (m_histogram_like) 
ss<<sw.label_line_text(llabel,xprior,yprior,xprior,y,2,"white","center",0) ; 
else
ss<<sw.label_line_text(llabel,xprior,yprior,xmid,y,2,"white","center",0) ; 
ss<<CRLF;
}
}
}

if (plot_symbol) { 
MM_ERR(" plotting symbol "<<MMPR4(xmid,y,symbol,pp.thick()))
plot_pp_symbol( ss, sw, pp, symbol,xmid ,y);
} // plot_symbol 

}

template<class Tj, class Ts>
void plot_seg_older(OsTy& ss, Sw& sw
, std::vector<D> & xv, std::vector<D> & yv, const D & x, const D & y
, const IdxTy ip
//, const bool label_foods, const bool histogram_like
//,const StrTy & llabel
//,  Tj & jj, Ts & s 
, const Tj & ipnext, Ts & s 
, const path_param & pp
,const Ps & ps ) 

{
//const StrTy & llabel=pp.m_seg_label;
const StrTy & llabel=pp.m_attmap["seg_label"]; // pp.m_seg_label;
//const StrTy & symbol=pp.m_symbol;
const StrTy & symbol=pp.m_attmap["symbol"];;
const bool plot_symbol=(symbol.length()!=0);
if (xv.size()>0){ 
D xlast=xv.back();
D ylast=yv.back();
if ((x!=xlast)||(y!=ylast)) { xv.push_back(x); yv.push_back(y);// }
MM_ERR(MMPR4(m_label_foods,ip,m_label_mod,m_label_res))
if ((m_label_foods) &&((int(ip)%m_label_mod)==m_label_res)){
ss<<sw.label_line_text(llabel,xlast,ylast,x,y,2,"white","center",0) ; 
ss<<CRLF;
}}}
if (plot_symbol) { 
MM_ERR(" plotting symbol "<<MMPR4(x,y,symbol,pp.thick()))
plot_pp_symbol( ss, sw, pp, symbol,x ,y);
} // plot_symbol 

if (ip==0) { xv.push_back(x); yv.push_back(y);}

//D ipnext=ip+1; // m_day_min;
if ( m_histogram_like) { 
//auto jjn=jj; ++jjn;
//if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
//}
// FIXME this makes an odd looking plot and may extend past the end
}
const D xbar=m_pb.utodx(ipnext); // m_lm+(ipnext)*m_catpitch;
xv.push_back(xbar);
yv.push_back(y);
}


template<class Tj, class Ts>
void plot_seg(OsTy& ss, Sw& sw
, std::vector<D> & xv, std::vector<D> & yv, const D & x, const D & y
, const IdxTy ip , const Tj & ipnext, Ts & s 
, const path_param & pp ,const Ps & ps ) 
{
//const StrTy & llabel=pp.m_seg_label;
const StrTy & llabel=pp.m_attmap["seg_label"]; // pp.m_seg_label;
//const StrTy & symbol=pp.m_symbol;
const StrTy & symbol=pp.m_attmap["symbol"];;
const bool plot_symbol=(symbol.length()!=0);
if (xv.size()>0){ 
D xlast=xv.back();
D ylast=yv.back();
if ((x!=xlast)||(y!=ylast)) { xv.push_back(x); yv.push_back(y);// }
MM_ERR(MMPR4(m_label_foods,ip,m_label_mod,m_label_res))
if ((m_label_foods) &&((int(ip)%m_label_mod)==m_label_res)){
ss<<sw.label_line_text(llabel,xlast,ylast,x,y,2,"white","center",0) ; 
ss<<CRLF;
}}}
if (plot_symbol) { 
MM_ERR(" plotting symbol "<<MMPR4(x,y,symbol,pp.thick()))
plot_pp_symbol( ss, sw, pp, symbol,x ,y);
} // plot_symbol 

if (ip==0) { xv.push_back(x); yv.push_back(y);}

//D ipnext=ip+1; // m_day_min;
if ( m_histogram_like) { 
//auto jjn=jj; ++jjn;
//if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
//}
// FIXME this makes an odd looking plot and may extend past the end
}
const D xbar=m_pb.utodx(ipnext); // m_lm+(ipnext)*m_catpitch;
xv.push_back(xbar);
yv.push_back(y);
}

template<class Tj, class Ts>
void plot_seg_old(OsTy& ss, Sw& sw
, std::vector<D> & xv, std::vector<D> & yv, const D & x, const D & y
, const IdxTy ip
//, const bool label_foods, const bool histogram_like
//,const StrTy & llabel
//,  Tj & jj, Ts & s 
, const Tj & ipnext, Ts & s 
, const path_param & pp
,const Ps & ps ) 

{
//const StrTy & llabel=pp.m_seg_label;
const StrTy & llabel=pp.m_attmap["seg_label"]; // pp.m_seg_label;
//const StrTy & symbol=pp.m_symbol;
const StrTy & symbol=pp.m_attmap["symbol"];;
const bool plot_symbol=(symbol.length()!=0);
if (xv.size()>0){ 
D xlast=xv.back();
D ylast=yv.back();
if ((x!=xlast)||(y!=ylast)) { xv.push_back(x); yv.push_back(y);// }
MM_ERR(MMPR4(m_label_foods,ip,m_label_mod,m_label_res))
if ((m_label_foods) &&((int(ip)%m_label_mod)==m_label_res)){
ss<<sw.label_line_text(llabel,xlast,ylast,x,y,2,"white","center",0) ; 
ss<<CRLF;
}}}
if (plot_symbol) { 
MM_ERR(" plotting symbol "<<MMPR4(x,y,symbol,pp.thick()))
plot_pp_symbol( ss, sw, pp, symbol,x ,y);
} // plot_symbol 

if (ip==0) { xv.push_back(x); yv.push_back(y);}

//D ipnext=ip+1; // m_day_min;
if ( m_histogram_like) { 
//auto jjn=jj; ++jjn;
//if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
//}
// FIXME this makes an odd looking plot and may extend past the end
}
const D xbar=m_pb.utodx(ipnext); // m_lm+(ipnext)*m_catpitch;
xv.push_back(xbar);
yv.push_back(y);
}


void plot_pp_symbol( Os& ss, Sw & sw, const path_param & pp,  const StrTy & symbol
, const D x, const D y )
{

if (pp.have_polygon())
{
ss<<sw.polygon_text(pp.m_x,pp.m_y,x,y,pp.color(),pp.color(),1); 
//StrTy polygon_text( const std::vector<D>  x, const std::vector<D> &  y
//,const StrTy & stroke
//, const StrTy & fill="#00ffffff",const D & sw=1) const
}
else if (symbol=="circle")
{
ss<<sw.circle_text(x,y,pp.m_symbol_size,pp.color(),pp.color(),1); 
} else 
{
ss<<sw.gtext_text(symbol, x, y+0*pp.m_symbol_size, pp.m_symbol_size, pp.color(),  StrTy("middle") ,0);
}
ss<<CRLF;

} // plot_pp_symbol





void write_svg_grid(OsTy & os, Sw& sw)
{

os<<sw.gtext_text(m_pb.m_ylab, m_pb.m_ylab_x, m_pb.m_ylab_y, m_pb.m_ylab_sz, m_pb.m_colorleg,  StrTy("end") ,90);
//StrTy hstripes(const D & xz,const D & xf, const D & yz, const D & yf,const IdxTy & n, const StrTy& col, const D & thick)
os<<sw.hstripes( m_pb.m_xz,m_pb.m_xf,m_pb.m_yz,m_pb.m_yf,m_ygrid,StrTy("black"), 1);
const bool  normal_x_grid = (m_pb.m_x_rule==0);
MM_ONCE(" using new vstripe method",)
if ( normal_x_grid) 
{os<<sw.vstripes( m_pb.m_xz,m_pb.m_xf,m_pb.m_yz,m_pb.m_yf,m_dl.range(), StrTy("black"), 1);
}else {
os<<sw.vstripes_del( m_pb.m_xz,m_pb.m_xf,m_pb.m_yz,m_pb.m_yf,m_pb.m_xpitch*m_pb.m_x_rule, StrTy("black"), 1);
}
}
void write_svg_ylabels(OsTy & os, Sw& sw)
{
const bool norm=(m_dl.m_yfunc_map.size()!=0);
if (!norm) { m_pb.limits(m_dl.m_vmax,m_dl.m_vmin);} 
for(IdxTy i=0; i<=m_ygrid; ++i)
{
const D f=D(i)/m_ygrid;
const D y=m_pb.utody(f); // m_tm+m_ysz*(1.0-f);
StrTy text= m_pb.y_label( f) ; // ,m_dl.m_vmax, m_dl.m_vmin);
#if 0
const D y=m_pb.utody(f); // m_tm+m_ysz*(1.0-f);
Ss ss;
ss<<std::setprecision(3);
D q=m_dl.descale(f); // ((m_dl.m_vmin+(m_dl.m_vmax-m_dl.m_vmin)*f));
//if ((q>=1e3)&&(q<1e5)) ss<<std::fixed<<std::setprecition(0);
if ((q>=1e3)&&(q<1e5)) ss<<int(q);
else ss<<q;
const StrTy & text=ss.str();
#endif

// this ASSFUDD does not anchor on fudding right or fudding levt 
os<<sw.gtext_text(text, m_pb.m_xleg, y, m_pb.m_szleg, m_pb.m_colorleg,  StrTy("end") ,0);
os<<CRLF;
}
}

void write_svg_day_xlabels(OsTy & os, Sw& sw)
{
D xlast=-1e100;
const D & min_dx= m_pb.min_x_space();
const IdxTy ll=m_xlabels.size();
MM_ERR(" XLABEL "<<MMPR2(min_dx,ll))
for(IdxTy i=0; i<ll; ++i)
{
const IdxTy ip=i; // -m_first_line;
const IdxTy is=i; // m_so[i];
//const Ragged::Line line=r.line(is);
//if (line.size()<=text_idx) continue; // fudd 
const StrTy & text= m_xlabels[is]; // line[text_idx];
//MM_ERR(MMPR3(text,is,i))
const D x=m_pb.legx(ip); // m_lm+ip*m_catpitch+m_szleg;
MM_ERR(" XLABEL "<<MMPR4(x,x-xlast,text,i))
bool occludes=(text.c_str()[0]!=0);
// min_dx==0 just for clarity although I guess x need not be monontinc..  
if (occludes) if ((min_dx==0)||((x-xlast)>=min_dx))
{
xlast=x;
const D y=m_pb.legy(ip); // m_lm+ip*m_catpitch+m_szleg;
//os<<sw.gtext_text(text, x, m_yleg, m_szleg, m_colorleg,  StrTy("end") ,-45);
os<<sw.gtext_text(text, x, y, m_pb.m_szleg, m_pb.m_colorleg,  StrTy("end") ,-45);
os<<CRLF;
} // dx

} // i 

//os<<sw.gtext_text(m_dl.m_dates_label,m_pb.utodx(m_dl.range()>>1) , m_pb.m_yleg+10*m_pb.m_szleg, m_pb.m_szleg, m_pb.m_colorleg,  StrTy("center") ,0);
os<<sw.gtext_text(m_dl.m_dates_label,m_pb.utodx(m_dl.range()>>1) , m_pb.m_date_y_loc, m_pb.m_szleg, m_pb.m_colorleg,  StrTy("center") ,0);
os<<CRLF;

}
void write_svg_legend(OsTy & os, Sw& sw, const IdxTy flags)
{
std::vector<D> x,y;
std::vector<StrTy> names;
D w,h,wr,hr;
w=m_pb.m_xsz*m_pb.m_leg_x_frac; // m_pb.m_leg_sz;
h=m_pb.m_tm; // m_pb.m_leg_sz;
{
//m_names.clear();
//m_colors.clear();
std::vector<path_param*> ppp;
MM_LOOP(ff,m_dl.m_foods_to_plot)
{
const StrTy & food=(*ff).first; // 
MM_LOOP(is,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*is).first;
// 2021-09-22 try to include scale
Ps & ps=m_dl.find_path_stats( dog, food);
MM_ERR(MMPR3(dog,food,ps.scale_value()))
D psc=ps.scale_value();
if (psc!=0) psc=1.0/psc;
//StrTy nc=dog+StrTy("-")+food;
path_param & pp=m_pm[dog][food];
ppp.push_back(&pp);
//m_names.push_back(nc);
const bool include_scale=true;
if (include_scale)
{
StrTy namescale=pp.name();

Ss ss; ss<<std::setprecision(3)<< (psc); 
MM_ERR(MMPR3(dog,food,ss.str()))
namescale+=StrTy("( s=")+ss.str()+")";
names.push_back(namescale);
}else
names.push_back(pp.name());



/////////////////m_colors.push_back(pp.m_color);
//m_colors.push_back(pp.m_attmap["color"]);

} // is
} // ff

// the names here are used to size the boxes, 
//m_pb.legend_box(x,y,wr,hr,w,h,m_names);
m_pb.legend_box(x,y,wr,hr,w,h,names);
MM_ERR(MMPR4(w,h,wr,hr))
// fudd
//MM_SZ_LOOP(i,m_names,sz)
MM_SZ_LOOP(i,ppp,sz)
{
const path_param & pp=* ppp[i];
// this is now linearized, why the fudd not names[i]???
const StrTy & _text=pp.name(); // m_names[i];
const StrTy & text= names[i];
MM_ERR(MMPR2(_text,text))
#if 0 
// 2021-09-22 try to include scale
StrTy text=_text;
const bool include_scale=true;
if (include_scale)
{
Ps & ps=m_dl.find_path_stats( dog, food);
MM_ERR(MMPR3(dog,food,ps.scale_value()))
D psc=ps.scale_value();
if (psc!=0) psc=1.0/psc;
Ss ss; ss<< (psc); 
MM_ERR(MMPR3(dog,food,ss.str()))
text+=StrTy("( s=")+ss.str()+")";
}
#endif

if (!ppp[i]->plotted()){
MM_ERR(" skipping non plotted entry "<<MMPR(text))
 continue;
}
{ plot_pp_symbol( os, sw, *ppp[i], ppp[i]->m_attmap["symbol"],x[i] ,y[i]); }
std::vector<D> xv,yv;
//m_pb.legend_strokes(xv,yv,x[i],y[i],wr,hr,pp.thick(),text.length());
m_pb.legend_strokes(xv,yv,x[i],y[i],wr,hr,pp.thick(),.25*wr);

if (pp.m_stroke_array.size()!=0) 
{os<<sw.vector_text(xv,yv,pp.m_stroke_array, pp.color(),pp.thick(),.7+.3);}
	else os<<sw.vector_text(xv,yv,pp.color(),pp.thick(),.7); 
  os<<CRLF;
os<<sw.gtext_text(text, x[i]+wr, y[i]+hr , m_pb.m_leg_sz, pp.color(),  StrTy("start") ,0);
os<<CRLF;

} // names

} // else

} // write_svg_legend


#if 0 
void write_svg_legend_old(OsTy & os, Sw& sw, const IdxTy flags)
{
std::vector<D> x,y;
D w,h,wr,hr;
w=m_pb.m_xsz; // m_pb.m_leg_sz;
h=m_pb.m_tm; // m_pb.m_leg_sz;

{
std::vector<path_param*> ppp;
std::vector<StrTy> names;
MM_LOOP(ff,m_dl.m_foods_to_plot)
{
const StrTy & food=(*ff).first; // 
MM_LOOP(is,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*is).first;
path_param & pp=m_pm[dog][food];
ppp.push_back(&pp);
names.push_back(pp.name());

} // is
} // ff

m_pb.legend_box(x,y,wr,hr,w,h,names);
MM_ERR(MMPR4(w,h,wr,hr))

//MM_SZ_LOOP(i,m_names,sz)
MM_SZ_LOOP(i,ppp,sz)
{
const path_param & pp=* ppp[i];
const StrTy & text=pp.name(); // m_names[i];
if (!ppp[i]->plotted()){
MM_ERR(" skipping non plotted entry "<<MMPR(text))
 continue;
}
 plot_pp_symbol( os, sw, *ppp[i], ppp[i]->m_attmap["symbol"],x[i] ,y[i]); }
const bool  plot_segment =true;
if ( plot_segment) 
{
std::vector<D> xv,yv;
m_pb.legend_strokes(xv,yv,x[i],y[i],wr,hr,pp.thick(),text.length());

//xv.push_back(x[i]);
//yv.push_back(y[i]+hr+pp.thick());
//xv.push_back(x[i]+wr+text.length()*m_pb.m_leg_sz);
//yv.push_back(y[i]+hr+pp.thick());
 
 
if (pp.m_stroke_array.size()!=0) 
os<<sw.vector_text(xv,yv,pp.m_stroke_array, pp.color(),pp.thick(),.7+.3); 

}
os<<sw.gtext_text(text, x[i]+wr, y[i]+hr , m_pb.m_leg_sz, pp.color(),  StrTy("start") ,0);
os<<CRLF;

}

} // else





} // write_svg_legend_old

#endif


void setup_day_plot(Ragged & r, ReadWriteMap & rwm, Ragged & pr)
{
m_dl.setup(r,rwm,pr);
setup_day_plot2(r,rwm,pr);
}
void setup_day_plot(Ragged & r, Ragged & pr)
{
// data r has spurious space and field shift
m_dl.setup(r,pr);
ReadWriteMap & rwm= m_dl.m_rwm;
m_markups=pr;
setup_day_plot2(r,rwm,pr);
}
void set_markups(Ragged &pr) { m_markups=pr; } 
void setup_day_plot2(Ragged & r, ReadWriteMap & rwm, Ragged & pr )
{
colors_and_day_names(rwm);
//m_markups=pr;
m_dl.load(r);
// uses period 
m_pb.setup_geometry(rwm,m_dl.range());

//m_default_color="red"; rwm.get("default_color",m_default_color);
setif(m_default_color,"red",rwm,"default_color");
setif(m_default_thick,2,rwm,"default_thick");

//m_default_thick=2; rwm.get("default_thick",m_default_thick);
m_default_sym=""; rwm.get("default_sym",m_default_sym);
m_default_lbl=0; rwm.get("default_lbl",m_default_lbl);
IdxTy label_foods=1; rwm.get("label_foods",label_foods); m_label_foods=(label_foods!=0);
 IdxTy hl=0; rwm.get("histogram_like",hl); m_histogram_like=(hl!=0);
m_histogram_code=hl;
m_label_mod=5; rwm.get("label_mod",m_label_mod);
m_label_res=1; rwm.get("label_res",m_label_res);

MM_ERR(MMPR(m_histogram_code))
MM_ERR(MMPR(m_pb.dump()))
setup_paths(r,pr);
}
template <class Tx, class Ty,class Tz > 
void setif(Ty & d, const Tz & def, ReadWriteMap & rwm, const Tx & c) { d=def; rwm.get(c,d); }

void setup_paths_old(Ragged & r, Ragged & pr)
{
//Strings rainbow;
//default_rainbod( rainbow);
Colors colors;
//const IdxTy cyc=rainbow.size();
StrTy color=m_default_color;
const D thick=m_default_thick;
 D opaq=1.0;
m_dl.m_rwm.get("opaq",opaq);
const StrTy sym=m_default_sym;
const IdxTy lbl=m_default_lbl;
IdxTy col=0;
MM_LOOP(ii,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*ii).first;
path_param def;
def.opaq(opaq);
const StrTy base=StrTy("path-")+dog;
def.setup(m_dl.m_rwm,base,StrTy(""),StrTy(""));
MM_LOOP(jj,m_dl.m_foods_to_plot)
{
const StrTy & food=(*jj).first;
path_param & pp=m_pm[dog][food];
//const StrTy base=StrTy("path-")+dog+StrTy("-")+food+StrTy("-");
const StrTy base=StrTy("path-")+pp.name()+StrTy("-");
//if (cyc>0) { color=rainbow[col%cyc]; ++col; } 
color=colors.next_contrast();
def.set(color,thick,sym,lbl);
pp=def;
pp.setup(m_dl.m_rwm,base,dog,food);
//pp.m_seg_label=dog+StrTy("-")+food;
//pp.m_seg_label_color=pp.color();
//pp.m_seg_label_size=pp.thick();

auto & dm=m_dl.m_dash_map;
hierarchy_set(&path_param::set_dashes,pp,dm,dog,food);
hierarchy_set(&path_param::set_symbol_size,pp,m_dl.m_symbol_size_map,dog,food);
//hierarchy_set(&path_param::set_symbol,pp,m_dl.m_symbol_map,dog,food);
//hierarchy_set(&path_param::set_symbol,pp,m_dl.m_maps["symbol"],dog,food);
hierarchy_set(&path_param::set_fudder,pp,m_dl.m_maps,dog,food, StrTy("symbol"));
//hierarchy_set(&path_param::set_color,pp,m_dl.m_color_map,dog,food);
hierarchy_set(&path_param::set_fudder,pp,m_dl.m_maps,dog,food, StrTy("color"));
hierarchy_set(&path_param::set_y_func,pp,m_dl.m_yfunc_map,dog,food);


} //jj
} // ii

} // setup_paths_old


void setup_paths(Ragged & r, Ragged & pr)
{
typedef mjm_rule_list<Tr> ConfigRules;
ConfigRules rl; // xxxxx
rl.add_lines(pr.begin(), pr.end()); //adfasdfadfs
//MM_MSG(rl.dump())
Colors colors;
StrTy color=m_default_color;
const D thick=m_default_thick;
 D opaq=1.0;
m_dl.m_rwm.get("opaq",opaq);
const StrTy sym=m_default_sym;
const IdxTy lbl=m_default_lbl;
IdxTy col=0;
MM_LOOP(ii,m_dl.m_dogs_to_plot)
{
const StrTy & dog=(*ii).first;
path_param def;
def.opaq(opaq);
const StrTy base=StrTy("path-")+dog;
def.setup(m_dl.m_rwm,base,StrTy(""),StrTy(""));
MM_LOOP(jj,m_dl.m_foods_to_plot)
{
const StrTy & food=(*jj).first;
path_param & pp=m_pm[dog][food];
color=colors.next_contrast();
// the name is not valid yet as base params not known 
//rule_list.apply(color,"color",pp.name());// xxxxxxxx
def.set(color,thick,sym,lbl);
pp=def;
pp.setup(rl,dog,food);

} //jj
} // ii

} // setup_paths



template <class Tx > static IdxTy xx_string( std::vector<Tx> & x)
{
return x.size();
}
static D xx_string( const D  & x) { return x; } 
static StrTy xx_string( const StrTy  & x)
{
return x;
}

template < class Tf,class Tp, class Td,class Tv >
static void hierarchy_set( Tf  f, Tp & pp, Td & dm, Tv & pk )
{
if (dm.find(pk)!=dm.end()) {
Ss ss;
for (IdxTy i=0; i<pk.size(); ++i) ss<<pk[i]<<" ";
MM_ERR(" seting  "<<MMPR2(ss.str(),xx_string(dm[pk])))
((&pp)->*f)(dm[pk]);}
}
template < class Tf,class Tp, class Td >
static void hierarchy_set( Tf  f, Tp & pp, Td & dm, const StrTy & dog, const StrTy & food )
{
data_loader::PlotKey pk(2);
pk[0]=StrTy("*"); pk[1]=StrTy("*");
hierarchy_set(f,pp,dm,pk);
//if (dm.find(pk)!=dm.end()) {
//MM_ERR(" seting  "<<MMPR3(dog,food,xx_string(dm[pk])))
//((&pp)->*f)(dm[pk]);}

pk[0]=StrTy("*");
pk[1]=food;
hierarchy_set(f,pp,dm,pk);
//if (dm.find(pk)!=dm.end()) {
//MM_ERR(" seting  "<<MMPR3(dog,food,xx_string(dm[pk])))
//((&pp)->*f)(dm[pk]);}

//} //jj
pk[0]=dog;
pk[1]=StrTy("*");
hierarchy_set(f,pp,dm,pk);
//if (dm.find(pk)!=dm.end()) {
//MM_ERR(" seting  "<<MMPR3(dog,food,xx_string(dm[pk])))
//((&pp)->*f)(dm[pk]);}
pk[0]=dog;
pk[1]=food;
hierarchy_set(f,pp,dm,pk);

//if (dm.find(pk)!=dm.end()) {
//MM_ERR(" seting  "<<MMPR3(dog,food,xx_string(dm[pk])))
//((&pp)->*f)(dm[pk]);}


}

//////////////////////////////////////////////////////////


template < class Tf,class Tp, class Td,class Tv >
static void hierarchy_set( Tf  f, Tp & pp, Td & dm, Tv & pk ,const StrTy & nm)
{
auto & dmi=dm[nm];
if (dmi.find(pk)!=dmi.end()) {
Ss ss;
for (IdxTy i=0; i<pk.size(); ++i) ss<<pk[i]<<" ";
MM_ERR(" seting  "<<MMPR2(ss.str(),xx_string(dmi[pk])))
((&pp)->*f)(nm,dmi[pk]);}
}
template < class Tf,class Tp, class Td >
static void hierarchy_set( Tf  f, Tp & pp, Td & dm, const StrTy & dog, const StrTy & food, const StrTy & nm  )
{
data_loader::PlotKey pk(2);
pk[0]=StrTy("*"); pk[1]=StrTy("*");
hierarchy_set(f,pp,dm,pk,nm);

pk[0]=StrTy("*");
pk[1]=food;
hierarchy_set(f,pp,dm,pk,nm);

pk[0]=dog;
pk[1]=StrTy("*");
hierarchy_set(f,pp,dm,pk,nm);

pk[0]=dog;
pk[1]=food;
hierarchy_set(f,pp,dm,pk,nm);


}




void xy(D & x, D & xend, D & y, D & yend, const Ragged::Line & line
, IdxTy & pos, const IdxTy flags=~0)
{
if (mybit(flags,0)){ x= m_pb.utodx(line[pos],m_dl.m_day_min); ++pos; }
if (mybit(flags,1)) {  xend=  m_pb.utodx(line[pos],m_dl.m_day_min); ++pos; }
if (mybit(flags,2)){ y= m_pb.utody(myatof(line[pos])); ++pos;}
if (mybit(flags,3)) {  yend= m_pb.utody(myatof(line[pos])); ++pos; }

}

class markup_state
{
public:
void set(const StrTy & k, const StrTy & v)
{
m_kv[k]=v;

}

StrTy get(const StrTy & k, const StrTy & d)
{
auto ii=m_kv.find(k);
if (ii==m_kv.end()) return d;
return (*ii).second;
}
D get(const StrTy & k, const D & d)
{
auto ii=m_kv.find(k);
if (ii==m_kv.end()) return d;
return atof(((*ii).second).c_str());
}

void set(const StrTy & k, const D & d)
{
Ss ss;
ss<<std::setprecision(16)<<d;
m_kv[k]=ss.str();
}

std::map<StrTy, StrTy> m_kv;

}; // markup_state
// wrong  line 2018-09-07 1 2018-11-01 0 blue 20 .2
typedef markup_state MarkVar; 
void markup_line_wtf(Os & ss, Sw& sw, const Ragged::Line & line, 
const IdxTy szl, IdxTy & pos, const StrTy & p1, const StrTy & p2, const StrTy & cmd, MarkVar & mv)
{
D x,xend,y,yend;
if ((szl+pos)<6)
{
MM_ERR(" line command too short "<<MMPR4(szl,cmd,p1,p2))
}
// TDO 15??? 
xy(x,xend,y,yend, line, pos,5);
const StrTy& cj= (line[pos]); ++pos;
const D stroke= myatof(line[pos]); ++pos;
const D opaq= myatof(line[pos]); ++pos;
//const StrTy anchor = (line[pos]); ++pos;
//StrTy text = (line[pos]); ++pos;
MM_ERR(" markup_line command  "<<MMPR4(x,y,xend,yend)<<MMPR3(cj,stroke,opaq))
ss<<sw.line_text(x, y, xend,yend, cj, stroke ,opaq);
}

void markup_line(Os & ss, Sw& sw, const Ragged::Line & line, 
const IdxTy szl, IdxTy & pos, const StrTy & p1, const StrTy & p2, const StrTy & cmd, MarkVar & mv)
{
D x,xend,y,yend;
if ((szl+pos)<6)
{
MM_ERR(" line command too short "<<MMPR4(szl,cmd,p1,p2))
}
xy(x,xend,y,yend, line, pos,15);
const StrTy& cj= (line[pos]); ++pos;
const D stroke= myatof(line[pos]); ++pos;
const D opaq= myatof(line[pos]); ++pos;
//const StrTy anchor = (line[pos]); ++pos;
//StrTy text = (line[pos]); ++pos;
MM_ERR(" markup_linef command  "<<MMPR4(x,y,xend,yend)<<MMPR3(cj,stroke,opaq))
ss<<sw.line_text(x, y, xend,yend, cj, stroke ,opaq);
}



// condition date name text 
void markup_event(Os & ss, Sw& sw, const Ragged::Line & line, 
const IdxTy szl, IdxTy & pos, const StrTy & p1, const StrTy & p2, const StrTy & cmd, MarkVar & mv)
{
D x,xend,y,yend;
if ((szl+pos)<6)
{
MM_ERR(" line command too short "<<MMPR4(szl,cmd,p1,p2))
}
const StrTy& name= (line[pos]); ++pos;
xy(x,xend,y,yend, line, pos,1); // flag 1 means only get one post 
 y= m_pb.utody(0); 
 yend= m_pb.utody(1);

StrTy type= (line[pos]); ++pos;
//while (pos<szl) { const StrTy & x= line[pos]; ++pos; if (x=="#") break; text=text+StrTy(" ")+x; } 
//y=.5;
D size= 1*m_pb.m_xpitch;
const StrTy cj= line[pos]; ++pos; // "red";
const D stroke=size;
const D opaq=1;
//const StrTy anchor="start";
//const D angle=90;
//const D stroke= myatof(line[pos]); ++pos;
//const D opaq= myatof(line[pos]); ++pos;
//const StrTy anchor = (line[pos]); ++pos;
//StrTy text = (line[pos]); ++pos;
//ss<<sw.line_text(x, y, xend,yend, cj, stroke ,opaq);
//ss<<sw.gtext_text(text, x, y, size, cj, anchor ,angle);
ss<<sw.line_text(x, y, x,yend, cj, stroke ,opaq);
} // markup_event




// condition date name text 
void markup_condition(Os & ss, Sw& sw, const Ragged::Line & line, 
const IdxTy szl, IdxTy & pos, const StrTy & p1, const StrTy & p2, const StrTy & cmd, MarkVar & mv)
{
D x,xend,y,yend;
if ((szl+pos)<6)
{
MM_ERR(" line command too short "<<MMPR4(szl,cmd,p1,p2))
}
const StrTy& name= (line[pos]); ++pos;
xy(x,xend,y,yend, line, pos,1);
StrTy text= (line[pos]); ++pos;
while (pos<szl) { const StrTy & x= line[pos]; ++pos; if (x=="#") break; text=text+StrTy(" ")+x; } 
y=.5;
D size= m_pb.m_xpitch;
const StrTy cj="red";
const StrTy anchor="start";
const D angle=90;
//const D stroke= myatof(line[pos]); ++pos;
//const D opaq= myatof(line[pos]); ++pos;
//const StrTy anchor = (line[pos]); ++pos;
//StrTy text = (line[pos]); ++pos;
//ss<<sw.line_text(x, y, xend,yend, cj, stroke ,opaq);
ss<<sw.gtext_text(text, x, y, size, cj, anchor ,angle);
}



void markup(Os & ss, Sw& sw, Ragged & r)
{
D x,xend,y,yend;
MarkVar mv;
const D top=m_pb.utody(1);
const D bottom=m_pb.utody(0);

const IdxTy rsz=r.size();
for(IdxTy i=0; i<rsz; ++i)
{

const Ragged::Line& line=r.line(i);
const IdxTy szl=line.size();
if (szl==0) continue;
if (false)
{Ss ss; MM_SZ_LOOP(_i,line,_sz) { ss<<line[_i]<<" "; } // _i
MM_ERR(MMPR(ss.str()))
}
const StrTy & cmd=line[0];
StrTy p1,p2;
if ( szl>1) p1=line[1];
if ( szl>2) p2=line[2];

if (cmd=="box")
{
IdxTy pos=1;
xy(x,xend,y,yend, line, pos);
const D w= xend-x; // 
const D h= yend-y; // 
const StrTy& cj= (line[pos]); ++pos;
const D op= myatof(line[pos]); ++pos;
MM_ERR(" writing box "<<MMPR4(x,y,w,h)<<MMPR2(cj,op))
ss<<sw.shade_rect_text(x,y,w,h,cj,op);
} // box
if (cmd=="condition")
{
IdxTy pos=1;
markup_condition(ss,sw,line,szl,pos,p1,p2,cmd,mv);
ss<<CRLF; 
continue;
}

if (cmd=="event")
{
IdxTy pos=1;
markup_event(ss,sw,line,szl,pos,p1,p2,cmd,mv);
ss<<CRLF; 
continue;
}


//StrTy line_text( const D& x0, const D& y0, const D & x1, const D & y1, const StrTy & color="#00ffffff",const D & sw=10, const D opaq=1.0) const
if (cmd=="line")
{
IdxTy pos=1;
markup_line(ss,sw,line,szl,pos,p1,p2,cmd,mv);
ss<<CRLF; 
continue;
}

if (cmd=="text")
{
IdxTy pos=1;
if (szl<7)
{
MM_ERR(" text command too short "<<MMPR4(szl,cmd,p1,p2))
}
xy(x,xend,y,yend, line, pos,5);
const StrTy& cj= (line[pos]); ++pos;
const D size= myatof(line[pos]); ++pos;
const D angle= myatof(line[pos]); ++pos;
const StrTy anchor = (line[pos]); ++pos;
StrTy text = (line[pos]); ++pos;
while (pos<szl){ text=text+StrTy(" ") +line[pos]; ++pos; } 
MM_ERR(" writing text "<<MMPR2(x,y)<<MMPR3(cj,angle,text))
ss<<sw.gtext_text(text, x, y, size, cj, anchor ,angle);

ss<<CRLF; 
continue;




} //text 
if (cmd=="bad-data")
{
IdxTy pos=1;
xy(x,xend,y,yend, line, pos,3);
const StrTy& cj= "red"; // (line[pos]); ++pos;
const D angle= 90; // myatof(line[pos]); ++pos;
const StrTy anchor = "middle"; // (line[pos]); ++pos;
StrTy text = "bad data "; // + (line[pos]); ++pos;
while (pos<szl){ text=text+StrTy(" ") +line[pos]; ++pos; } 
const D op=.2;
MM_ERR(" writing bad text  "<<MMPR2(x,y)<<MMPR3(cj,angle,text))
const D w=xend-x;
const D y=bottom;
const D h=top-bottom;
ss<<sw.shade_rect_text(x,bottom,w,top-bottom,cj,op);
const D size= h*.01;
ss<<sw.gtext_text(text, x+.5*w, y+.5*h, size, "black", anchor ,angle);

ss<<CRLF; 
continue;



} //text 
if (cmd=="note-box")
{
IdxTy pos=1;
xy(x,xend,y,yend, line, pos,~0);
const StrTy& cj= (line[pos]); ++pos;
// D angle=  myatof(line[pos]); ++pos;
const StrTy anchor = (line[pos]); ++pos;
StrTy text = (line[pos]); ++pos;
while (pos<szl){ text=text+StrTy(" ") +line[pos]; ++pos; }
const D op=.2;
const D w=xend-x;
const D h=yend-y;
ss<<sw.shade_rect_text(x,y,w,h,cj,op);
D angle=0;
D size= w*.01;
if (w>h){  angle=90; size=h*.01; }
ss<<sw.gtext_text(text, x+.5*w, y+.5*h, size, "black", anchor ,angle);
MM_ERR(" writing bad text  "<<MMPR2(x,y)<<MMPR3(cj,angle,text))

ss<<CRLF; 
continue;


} //text 






#if 0 
os<<sw.shade_rect_text(x[i],y[i],wr,hr,cj,.8);
return m_pb.utodx(int((*semii).second)-int(m_dl.m_day_min)); }
m_pb.utodx(s,m_dl.m_day_min);
ss<<sw.label_line_text(llabel,xprior,yprior,xprior,y,2,"white","center",0) ; 
ss<<sw.gtext_text(symbol, x, y+0*pp.m_symbol_size, pp.m_symbol_size, pp.color(),  StrTy("middle") ,0);
 ShapeGen::ngon(m_x,m_y,sz,n,0,0,phi,m); 
	os<<sw.vector_text(xv,yv,pp.color(),pp.thick(),.7); 
#endif


} //i

} // makrup 

void markup_end(Os & ss, Sw& sw, Ragged & r)
{
/*
D x,xend,y,yend;
const D top=m_pb.utody(1);
const D bottom=m_pb.utody(0);

const IdxTy rsz=r.size();
for(IdxTy i=0; i<rsz; ++i)
{
const Ragged::Line& line=r.line(i);
const IdxTy szl=line.size();
if (szl<2) continue; 
if (false)
{Ss ss; MM_SZ_LOOP(_i,line,_sz) { ss<<line[_i]<<" "; } // _i
MM_ERR(MMPR(ss.str()))
}
const StrTy & level=line[0];
if (level!=StrTy("mend")) continue;
const StrTy & cmd=line[1];
IdxTy pos=2;
//void xy(D & x, D & xend, D & y, D & yend, const Ragged::Line & line
//, IdxTy & pos, const IdxTy flags=~0)
//if (mybit(flags,0)){ x= m_pb.utodx(line[pos],m_dl.m_day_min); ++pos; }
D x,y,xend,yend;
if (cmd=="")
{
// get x and y from line at pos and pos+1
xy(x,xend,y,yend,line,pos,5);

} // 

} // end

*/

} // markup_end





//////////////////////////////////////////////////////////


void write_day_svg(OsTy & os, Sw& sw)
{
MM_ERR("write_svg")
os<<sw.start_text(" test foo",m_pb.m_xs,m_pb.m_ys);
os<<sw.frame_text("#FFFFFF",m_pb.m_xs,m_pb.m_ys);
os<<CRLF;
// need this to get the max and min ut shold draw laer 
//write_svg_day_data(os, sw,2);
MM_ERR(" writing grid"<<MMPR3(m_dl.m_vmax,m_dl.m_vmin,m_ygrid))
m_ygrid=data_rounding::adjust_grid(m_dl.m_vmax,m_dl.m_vmin,m_ygrid);
markup(os, sw, m_markups);
write_svg_grid(os,sw);
MM_ERR(" writing data")
// need to call begin_plotting to generate smbols first 
if ( m_histogram_code>255){ write_svg_day_data_histo(os, sw,5);}
else {write_svg_day_data(os, sw,5);}
MM_ERR(" writing xlabels")
write_svg_day_xlabels(os,sw);
MM_ERR(" writing ylabels")
write_svg_ylabels(os,sw);
write_svg_legend(os,sw,0);
markup_end(os, sw, m_markups);
MM_ERR(" writing done")

os<<sw.end_text();
os<<CRLF;

} // write_day_svg



protected:

ReadWriteMap   m_rwm;
plot_box m_pb;
data_loader m_dl;

typedef std::map<StrTy, std::map< StrTy, path_param> > PathMap;
PathMap m_pm;

StrTy m_default_color;
D m_default_thick;
StrTy  m_default_sym;
IdxTy m_default_lbl;

//std::vector<StrTy> m_names,m_xlabels; // , m_colors;
std::vector<StrTy> m_xlabels; // , m_colors;
bool m_label_foods, m_histogram_like;
IdxTy m_histogram_code;
//StrTy  m_colorleg;
IdxTy m_label_mod,m_label_res,m_ygrid;

Ragged m_markups;

};  //layout_blocks




/////////////////////////////////////////////////////////

#ifdef  TEST_DOG_PLOTS__
int main(int argc,char **args)
{
//typedef mjm_linc_graph  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif

#undef MM_DMEL 

#endif
