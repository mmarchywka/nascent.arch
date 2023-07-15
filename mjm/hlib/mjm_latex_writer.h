#ifndef MJM_LATEX_WRITER_H__
#define MJM_LATEX_WRITER_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
//#include "mjm_rational.h"
 // #include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include "mjm_cli_ui.h"


#include "mjm_collections.h"


//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
//#include "mjm_biom_hdf5.h"
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
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 



////////////////////////////////////////////////////////////////

class latex_writer_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
latex_writer_params( const StrTy & nm) : Super(nm) {}
latex_writer_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",true); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands"<<log_commands()<<sep;
ss<<"exit_on_err"<<exit_on_err()<<sep;
//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace latex_writer_traits
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
}; // trees_and_tables_traits



class mjm_latex_writer 
{


typedef  latex_writer_traits::Tr  Tr;
typedef mjm_latex_writer Myt;

typedef mjm_cli_ui<Myt> CliTy;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef latex_writer_params Logic;
typedef mjm_logic_base VariableStore;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;
typedef std::map<StrTy,TaxTree> TaxTrees;

public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
mjm_latex_writer():m_dmel(new Dmel()) {Init();}
mjm_latex_writer(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_latex_writer()
{
//clear_handlers();
// this is often shared, wtf TODO FIXME
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


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

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

void cmd_help(Cip & cip , LocalVar & lv ) 
{ MM_LOOP(ii,m_cmd_map) { MM_ERR(MMPR((*ii).second)) } }


StrTy start_text() const
{
Ss ss;
ss<<"\\documentclass{article}"<<CRLF;
ss<<"\\usepackage{multirow}"<<CRLF;
ss<<"\\begin{document}"<<CRLF;
return ss.str();
}
StrTy end_text() const
{
Ss ss;
ss<<"\\end{document}"<<CRLF;
return ss.str();
}

void line(std::ostream & ss, const StrTy & x) const  { ss<<x<<CRLF;} 
//https://tex.stackexchange.com/questions/131867/using-multicolumn-in-latex
StrTy tree_table_start_text(const StrTy & caption, const StrTy & label, const IdxTy depth, const IdxTy vals) const
{
Ss ss;

line(ss,"\\begin{table}");
line(ss,"\\centering");
line(ss,StrTy("\\caption{")+caption+StrTy("}"));
line(ss,StrTy("\\label{")+label+StrTy("}"));
StrTy cx="{";
for (IdxTy i=0; i<depth; ++i) cx=cx+StrTy("r|");
//for (IdxTy i=0; i<vals; ++i) cx=cx+StrTy("r@{.}l|");
for (IdxTy i=0; i<vals; ++i) cx=cx+StrTy("l|");
cx=cx+StrTy("}");
line(ss,StrTy("\\begin{tabular}[\\textwidth]")+StrTy(cx)); // {c|c|c|c|c|c|c|c|c|}
return ss.str();
}

StrTy tree_table_end_text() const
{
Ss ss;

line(ss,"\\end{tabular}");
line(ss,"\\end{table}");
return ss.str();
}
StrTy escape_dash(const StrTy & s) const
{
const IdxTy sz=s.length();
char c[sz*3];
const char * sp=s.c_str();
IdxTy dptr=0;
for (IdxTy i=0; i<sz; ++i)
{if (sp[i]=='_') { c[dptr]='\\'; ++dptr; } 
c[dptr]=sp[i]; ++dptr;
}
c[dptr]=0;
return StrTy(c);

}
StrTy help_split(const StrTy & s) const
{
const IdxTy sz=s.length();
char c[sz*3];
const char * sp=s.c_str();
IdxTy dptr=0;
for (IdxTy i=0; i<sz; ++i)
{
if (sp[i]==',') { 
c[dptr]='\\'; ++dptr;
c[dptr]='\\'; ++dptr;

 } 
c[dptr]=sp[i]; ++dptr;
if (sp[i]>='a')if (sp[i]<='z') if(sp[i+1]>='0') if (sp[i+1]<='9')  {
 c[dptr]='\\'; ++dptr;
 c[dptr]='\\'; ++dptr;

 } 
}
c[dptr]=0;
return StrTy(c);

}



template <class Ty> 
StrTy tree_table_text(const StrTy & w, const Ty & v, const IdxTy level, const IdxTy depth, const IdxTy vals, const IdxTy precision) const
{
const IdxTy space=level;
const IdxTy taxon=depth-level;
Ss ss;
ss.precision(precision);
StrTy spacer="";
for (IdxTy i=0; i<space; ++i) spacer=spacer+StrTy(".");
ss<<"\\multicolumn{"<<space<<"}{|l|}{"<<spacer<<"}";
ss<<"&\\multicolumn{"<<taxon<<"}{|l|}{"<<w<<"}"; // \\\\"<<CRLF;
MM_SZ_LOOP(i,v,sz) {Ss ns; ns.precision(precision); ns<<std::scientific<<v[i];  ss<<"&"<<ns.str(); }
ss<<"\\\\"<<CRLF;

return ss.str();
}
template <class Ty> 
StrTy tree_table_header_text(const StrTy & w, const Ty & v, const IdxTy depth, const IdxTy vals, const IdxTy precision) const
{
const IdxTy level=1; 
const IdxTy space=level;
const IdxTy taxon=depth-level;
IdxTy rows=1;
MM_LOOP(ii,v) { IdxTy ri=1+((*ii).length()>>3); if (ri>rows) rows=ri; } 
Ss ss;
ss.precision(precision);
StrTy spacer="";
for (IdxTy i=0; i<space; ++i) spacer=spacer+StrTy(".");
ss<<"\\multicolumn{"<<space<<"}{|l|}{\\multirow{"<<rows<<"}{*}{"<<spacer<<"}}";
ss<<"&\\multicolumn{"<<taxon<<"}{l}{\\multirow{"<<rows<<"}{*}{"<<w<<"}}"; // \\\\"<<CRLF;
MM_SZ_LOOP(i,v,sz) {Ss ns; ns.precision(precision); 
ns<<std::scientific<<v[i];  
//ss<<"&"<<escape_dash(ns.str()); 
ss<<"&\\multicolumn{1}{|p{.5in}|}{\\multirow{"<<rows<<"}{*}{\\parbox[t]{.5in}{"<<escape_dash(help_split(ns.str()))<<"}}}"; 

}
ss<<"\\\\"<<CRLF;

return ss.str();
}


#if 0 

StrTy start_text(const StrTy &desc, const D & w, const D & h) const
{
Ss ss;
ss<<"<svg style=\"overflow: hidden; -moz-user-select: none; cursor: default;";
ss<<" position: relative; background-color: rgb(255, 255, 255);\"" ;
ss<<" xmlns=\"http://www.w3.org/2000/svg\" width=\"";
ss<<w<<"\" version=\"1.1\" height=\""<<h<<"\">";
ss<<"<desc>";
ss<<desc;
ss<<"</desc>";

return ss.str();
}
StrTy frame_text(const StrTy &color, const D & w, const D & h) const
{
Ss ss;

ss<<"<rect fill-opacity=\"1\" stroke-width=\"1\" stroke-opacity=\"1\" style=\"\" stroke=\"";
ss<<color<<"\" fill=\""<<color<<"\" ry=\"0\" rx=\"0\" r=\"0\"";
ss<<" height=\""<<h<<"\" width=\""<<w<<"\" y=\""<<0<<"\" x=\""<<0<<"\">";
ss<<"</rect>";

return ss.str();
}

StrTy stroke_text(const StrTy &color, const D & w, const StrTy & fill="none") const
{
Ss ss;

ss<<" <g fill=\""<<fill<<"\" stroke=\""<<color<<"\" stroke-width=\""<<w<<"\" />";

return ss.str();
}

StrTy line_text( const D& x0, const D& y0, const D & x1, const D & y1, const StrTy & color="#00ffffff",const D & sw=10) const
{
Ss ss;
//const StrTy fill="transparent";
//ss<<"<path d=\"M"<<x0<<" "<<y0<<" L"<<x1<<" "<<y1<<"\"  fill=\""<<fill<< "\" stroke=\""<<color<<"\"/>";
ss<<"<path d=\"M"<<x0<<" "<<y0<<" L"<<x1<<" "<<y1
<<"\" stroke=\""<<color<<"\""
<<" stroke-width=\""<<sw<<"\"/>";

return ss.str();
}
StrTy faintline_text( const D& x0, const D& y0, const D & x1, const D & y1, const StrTy & color="#00ffffff",const D & sw=10) const
{
Ss ss;
//const StrTy fill="transparent";
//ss<<"<path d=\"M"<<x0<<" "<<y0<<" L"<<x1<<" "<<y1<<"\"  fill=\""<<fill<< "\" stroke=\""<<color<<"\"/>";
ss<<"<path d=\"M"<<x0<<" "<<y0<<" L"<<x1<<" "<<y1<<"\""
<<add_kv("stroke",color)
<<add_kv("stroke-width",sw)
<<add_kv("stroke-dasharray",StrTy("1 .5 "))
//<<"\" stroke=\""<<color<<"\""
//<<" stroke-width=\""<<sw<<"\"/>";
<<"/>";

return ss.str();
}



StrTy text_text(const StrTy &text, const D& x, const D& y, const D & w, const D & h, const StrTy & color="#00ffffff") const
{
Ss ss;
//const StrTy color="00ff00ff";
if (false) { 
ss<<"<rect  fill=\"#"<<color<<"\" ry=\"0\" rx=\"0\" r=\"0\"";
ss<<" height=\""<<h<<"\" width=\""<<w<<"\" y=\""<<y<<"\" x=\""<<x<<"\">";
ss<<"</rect>";
}
ss<<"<text transform=\"matrix(1,0,0,1,0,0)\" font-weight=\"bold\"";
ss<<" font-size=\"8px\" font-family=\"Verdana\" fill=\"";
ss<<""<<color<<"\" stroke=\"none\" text-anchor=\"middle\"";
ss<<" y=\""<<y<<"\" x=\""<<x<<"\" style=\"text-anchor: left;";
ss<<" font-family: Verdana; font-size: 12px; font-weight: bold;\">";
ss<<text;
ss<<"</text>";

return ss.str();
}
StrTy add_kv( const StrTy & k, const StrTy & v)const { return StrTy(" ")+k+StrTy("=\"")+v+StrTy("\" "); } 
StrTy add_kv( const StrTy & k, const D & v)const {Ss ss; ss<<v;  return StrTy(" ")+k+StrTy("=\"")+ss.str()+StrTy("\" "); } 
StrTy vtext_text(const StrTy &text, const D& x, const D& y, const D & sz, const StrTy & color="#00ffffff") const
{
Ss ss;
//ss<<"<text transform=\"rotate(90,"<<x<<","<<y<<")\" font-weight=\"bold\"";
ss<<"<text transform=\"rotate(90,"<<x<<","<<y<<")\" font-weight=\"normal\"";
ss<<" font-size=\""<<sz<<"\" font-family=\"Verdana\" fill=\"";
ss<<""<<color<<"\""; //  stroke=\"none\" text-anchor=\"left\"";
ss<<" y=\""<<y<<"\" x=\""<<x<<"\" style=\"text-anchor: left;\"";
//ss<<" font-family: Verdana; font-size: 12px; font-weight: bold;\">";
// ss<<add_kv("dominant-baseline","alphabetic");
ss<<">";
ss<<text;
ss<<"</text>";
return ss.str();
}
StrTy vtextb_text(const StrTy &text, const D& x, const D& y, const D & sz, const StrTy & color="#00ffffff") const
{
Ss ss;
//ss<<"<text transform=\"rotate(90,"<<x<<","<<y<<")\" font-weight=\"bold\"";
ss<<"<text transform=\"rotate(90,"<<x<<","<<y<<")\" font-weight=\"normal\"";
ss<<" font-size=\""<<sz<<"\" font-family=\"Verdana\" fill=\"";
ss<<""<<color<<"\""; //  stroke=\"none\" text-anchor=\"left\"";
ss<<" y=\""<<y<<"\" x=\""<<x<<"\" style=\"text-anchor: end;";
//ss<<" font-family: Verdana; font-size: 12px; font-weight: bold;\">";
ss<<"\">";
ss<<text;
ss<<"</text>";
return ss.str();
}




StrTy htext_text(const StrTy &text, const D& x, const D& y, const D & sz, const StrTy & color="#00ffffff") const
{
Ss ss;
//ss<<"<text transform=\"rotate(90,"<<x<<","<<y<<")\" font-weight=\"bold\"";
ss<<"<text transform=\"rotate(0,"<<x<<","<<y<<")\" font-weight=\"normal\"";
ss<<" font-size=\""<<sz<<"\" font-family=\"Verdana\" fill=\"";
ss<<""<<color<<"\""; //  stroke=\"none\" text-anchor=\"left\"";
ss<<" y=\""<<y<<"\" x=\""<<x<<"\" style=\"text-anchor: left;";
//ss<<" font-family: Verdana; font-size: 12px; font-weight: bold;\">";
ss<<"\">";
ss<<text;
ss<<"</text>";

return ss.str();
}





StrTy rect_text( const D& x, const D& y, const D & w, const D & h, const StrTy & color="#00ffffff") const
{
Ss ss;
//const StrTy color="00ff00ff";
if (!false) { 
ss<<"<rect  fill=\""<<color<<"\" ry=\"0\" rx=\"0\" r=\"0\"";
ss<<" height=\""<<h<<"\" width=\""<<w<<"\" y=\""<<y<<"\" x=\""<<x<<"\">";
ss<<"</rect>";
}
return ss.str();
}

StrTy shade_rect_text( const D& x, const D& y, const D & w, const D & h, const StrTy & color="#00ffffff", const D & opaq=1.0) const
{
Ss ss;
//const StrTy color="00ff00ff";
if (!false) { 
ss<<"<rect  fill=\""<<color<<"\" ry=\"0\" rx=\"0\" r=\"0\"";
ss<<add_kv("fill-opacity",opaq); 
ss<<" height=\""<<h<<"\" width=\""<<w<<"\" y=\""<<y<<"\" x=\""<<x<<"\">";
ss<<"</rect>";
}
return ss.str();
}




StrTy end_text() const
{
Ss ss;
ss<<"</svg>";

return ss.str();
}


//os<<sw.hstripes(0,w,0,h,10,"black");
StrTy hstripes(const D & xz,const D & xf, const D & yz, const D & yf,const IdxTy & n, const StrTy& col, const D & thick)
{
Ss ss; 
D y=yz;
const D del=(yf-yz)/D(n);
for (IdxTy i=0; i<=n; ++i)
{
ss<< line_text(  xz,  y, xf,  y, col,thick) ;
y=y+del;
}


return ss.str();
}


#endif




void cmd_start(Cip & cip , LocalVar & lv ) 
{
OsTy & os=*cip.m_os;
//os<<start_text() ; // 



}

//void cmd_navigate(Cip & cip , LocalVar & lv,TaxTree & tt, IdxTy node ) 


void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("write-start")]=&Myt::cmd_start;
//m_cmd_map[StrTy("ncurses")]=&Myt::cmd_ncurses;
//m_cmd_map[StrTy("load-tax-nodes")]=&Myt::cmd_load_tax_nodes;
//m_cmd_map[StrTy("navigate")]=&Myt::cmd_navigate;
//m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;
//m_cmd_map[StrTy("tt-bin")]=&Myt::cmd_tt_bin;


}
#if 0 

bool tax_tree_commands(CommandInterpretter & li, const StrTy & cmd,TaxTree & tt,const StrTy & p1 ,const StrTy & p2)
{
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
if (cmd=="best-taxon") { 
std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  tt.best_taxon(ww,4);  continue; 
} 
if (cmd=="best-indexed-taxon") { 
std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  tt.best_indexed_taxon(ww,4);  continue; 
} 
// TODO FIXME has no way to return this now, should have local vars etc
if (cmd=="node-info") { StrTy local_label=tt.node_info(myatoi(p1),myatoi(p2)); MM_MSG((local_label))   continue; } 
return false; 
} while (false);
return true;
}
#endif

 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
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
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 

//m_tax_tree.standard_commnds(cmd,p1,p2,li);
//if ( tax_tree_commands(li,cmd,m_tax_tree,p1,p2)) continue ;
//if (cmd=="load-tax-nodes") { m_tax_tree.read_ncbi_nodes_dmp(p1);  continue; } 
//if (cmd=="load-tax-tree") { m_tax_tree.read_ncbi_dmp(p1);  continue; } 
//if (cmd=="dump-tax-tree") { m_tax_tree.dump(p1);  continue; } 
//if (cmd=="write-tax-single") { m_tax_tree.write_single(p1);  continue; } 
//if (cmd=="save-tax") { m_tax_tree.write_composite(p1);  continue; } 
//if (cmd=="load-tax") { m_tax_tree.read_composite(p1);  continue; } 
//if (cmd=="dump-lineages") { m_tax_tree.dump_lineages(p1,2);  continue; } 
//if (cmd=="dump-normal-lineages") { m_tax_tree.dump_lineages(p1,3);  continue; } 
//if (cmd=="check-normal-lineages") { m_tax_tree.dump_lineages(p1,4);  continue; } 
//if (cmd=="traverse-full-tree") { m_tax_tree.traverse_full_tree(p1,4);  continue; } 
//if (cmd=="best-taxon") { 
//std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  m_tax_tree.best_taxon(ww,4);  continue; 
//} 
//if (cmd=="best-indexed-taxon") { 
//std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  m_tax_tree.best_indexed_taxon(ww,4);  continue; 
//} 
//if (cmd=="node-info") { local_label=m_tax_tree.node_info(myatoi(p1),myatoi(p2)); MM_MSG((local_label))   continue; } 


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
ss<<" mjm_svg_writer "<<__DATE__<<" "<<__TIME__<<CRLF;
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
//MM_MSG(" configuration "<<m_flp.to_string())
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


/////////////////////////////////////////////////////////////////////
private:
void Init()
{
m_done=false;
}


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

Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//TaxTree m_tax_tree;
//TaxTrees m_tax_trees;

CounterMap m_cm;
CliTy m_cli;


CmdMap m_cmd_map;
CompMap m_comp_map;


}; //mjm_trees_and_tables



/////////////////////////////////////////////////////////

#ifdef  TEST_SVG_WRITER__
int main(int argc,char **args)
{
typedef mjm_svg_writer  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

