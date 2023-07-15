#ifndef MJM_TREES_AND_TABLES_H__
#define MJM_TREES_AND_TABLES_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
#include "mjm_rational.h"
 // #include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include "mjm_cli_ui.h"


#include "mjm_collections.h"
#include "mjm_svg_writer.h"


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
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 



////////////////////////////////////////////////////////////////

class trees_and_tables_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
trees_and_tables_params( const StrTy & nm) : Super(nm) {}
trees_and_tables_params() : Super() {}
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



namespace trees_and_tables_traits
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



class mjm_trees_and_tables 
{


typedef  trees_and_tables_traits::Tr  Tr;
typedef mjm_trees_and_tables Myt;

typedef mjm_cli_ui<Myt> CliTy;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef trees_and_tables_params Logic;
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
mjm_trees_and_tables():m_dmel(new Dmel()) {Init();}
mjm_trees_and_tables(int argc,char **_args) : m_dmel(new Dmel())
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
~mjm_trees_and_tables()
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


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

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

//void start() { initscr(); }
//void finish() { endwin(); }
//char show() { refresh(); return getch(); }
//void branch(const StrTy & base, const std::vector<StrTy> & kids)

class node_nav_state { public : 
node_nav_state(const IdxTy _i, const IdxTy _iold, const IdxTy fc, const IdxTy ks):i(_i),iold(_iold),firstchild(fc),kidsz(ks) {}

node_nav_state():i(0),iold(0),firstchild(0),kidsz(0) {}
void operator()(IdxTy& _i, IdxTy& _iold, IdxTy& fc, IdxTy& ks)
{_i=i; _iold=iold; fc=firstchild; ks=kidsz;}
IdxTy i,iold,firstchild,kidsz ; };

void cmd_tt(Cip & cip , LocalVar & lv ) 
{
const StrTy & cmd=cip.p1;
const StrTy & name=cip.p2;
const StrTy p1=cip.wif(3);
const StrTy p2=cip.wif(4);
if (cmd=="list")
{
MM_LOOP(ii,m_tax_trees) { MM_ERR(MMPR2((*ii).first, (*ii).second.size())) } 
return; 
}

TaxTree & tt=m_tax_trees[name];
MM_ERR(" cmd_tt "<<MMPR2(name,cmd)<<MMPR3(p1,p2,tt.size()))
if ( tax_tree_commands(cip.m_li,  cmd,tt, p1 , p2))
{
	MM_ERR(" cmd_tt "<<MMPR2(name,tt.size()))
 	return;
}
MM_ERR(" did nothing with "<<MMPR2(cmd,name))
} // cmd_tt

void cmd_tt_bin(Cip & cip , LocalVar & lv ) 
{
const StrTy & cmd=cip.p1;
const StrTy & name1=cip.p2;
const StrTy name2=cip.wif(3);
const StrTy nodeq=cip.wif(4);
//const StrTy p2=cip.wif(5);

TaxTree & tt1=m_tax_trees[name1];
TaxTree & tt2=m_tax_trees[name2];
MM_ERR(" cmd_tt_bin "<<MMPR3(name1,tt1.size(),tt2.size()))
do
{
if (cmd=="eqnode") // tt eqnode q r, find best node in r that matches node in q
{
translate_tt1_to_tt2(tt1,tt2);

continue;
}


MM_ERR(" did nothing with "<<MMPR3(cmd,name1,name2))
} while (false); 
} // cmd_tt_bin


void cmd_navigate(Cip & cip , LocalVar & lv ) 
{
const StrTy& name=cip.p2;
TaxTree & tt=(name.size()==0)?m_tax_tree:m_tax_trees[name];
IdxTy node=myatoi(cip.p1);
cmd_navigate(cip,lv,tt,node);
}

void cmd_navigate(Cip & cip , LocalVar & lv,TaxTree & tt, IdxTy node ) 
{ 
typedef mjm_ncurses_links<> NCurse; 
typedef NCurse::ScreenLocations Sl;
std::map<IdxTy, node_nav_state> nostate;
Sl sl; 
NCurse foo; 
IdxTy row,col,i,iold, firstchild,kidsz;
i=0;
firstchild=0;
iold=i;
IdxTy done=0;
foo.start();
foo.clear();
tt.sort_for_ui();
bool draw_again=true;
while (done==0)
{
const StrTy base=tt.node_name(node);
//MM_ERR(MMPR3(cip.p1,base,node))
if (draw_again) {
foo.clear();
sl.clear();
std::vector<StrTy> kids= tt.kids(node);
kidsz=kids.size();
foo.branch(base,kids,firstchild,&sl);
}
const IdxTy szlist=sl.size();
if (i>=szlist) i=szlist-1;
const IdxTy overhang=(kidsz>(szlist-1))?(kidsz+1-szlist):0;
row=sl[i].y;
col=sl[i].x;
//if (i>0) row=row-firstchild;
move(row,col);
// TODO char here was not flagged 
int  c=foo.show();
//MM_ERR(MMPR(int(c)))
switch (c)
{
case 'j':
case KEY_DOWN: 
{ i=(i+1)%szlist; 
if ((i==0)&&( overhang>firstchild)) 
{ ++firstchild;i=szlist-1;   
if (overhang!=0){  firstchild=(firstchild%(overhang));
if (firstchild==0) i=1; 
} 
else {i=1; firstchild=0; } 
}
if (i==0) i=1;
  break; }
case 'k':
case KEY_UP: { i=(i+szlist-1)%szlist; 
if ((i==1)&&(firstchild!=0)) { --firstchild; i=1; } 

break; }
case 'h':
case KEY_LEFT: {
if (i==0) { 
nostate[node]=node_nav_state(i,iold,firstchild,kidsz);
node=tt.parent(node);  
nostate[node](i,iold,firstchild,kidsz);
foo.clear(); break; }
iold=i; i=0; break; }
case 'l':
case KEY_RIGHT: {
if (i!=0) { 
nostate[node]=node_nav_state(i,iold,firstchild,kidsz);
node=tt.child(node,i-1+firstchild); 
nostate[node](i,iold,firstchild,kidsz);
foo.clear(); break; }
if (i==0 )  i=iold;if(i==0) i=(1%szlist);  break; }
case '\n':
case '\r': {
nostate[node]=node_nav_state(i,iold,firstchild,kidsz);
if(i==0) node=tt.parent(node); else  node=tt.child(node,i-1+firstchild); 
nostate[node](i,iold,firstchild,kidsz);
foo.clear(); break; }
default: {done=1;}
}; // sxitch 
} // done
foo.finish();

 }

void cmd_write_svg(Cip & cip , LocalVar & lv ) 
{
const StrTy& name=cip.p2;
TaxTree & tt=(name.size()==0)?m_tax_tree:m_tax_trees[name];
TaxTree::tree_layout lc;
typedef TaxTree::tree_layout::position Pos;
tt.traverse_full_tree(lc, StrTy(""), 0);

mjm_svg_writer sw;
OsTy & os=*cip.m_os; // 
if (true)
{
const IdxTy xpitch=200;
const IdxTy w=xpitch;
const IdxTy ypitch=40;
const IdxTy h=ypitch;
//mjm_output_warp ow(xpitch,ypitch,w,h);
mjm_circular_warp ow(xpitch,ypitch,w,h);
IdxTy x=0;
IdxTy y=0;
/*
IdxTy maxsz=0;
IdxTy maxlevel=0;
MM_LOOP(ii,lc.m_level_sz)
{
if ((*ii).first>maxlevel) maxlevel=(*ii).first;
if ((*ii).second>maxsz) maxsz=(*ii).second;
}
*/

//ow_bounds(x,y,maxlevel,maxsz);
ow.bounds(x,y,lc.m_level_sz);
os<<sw.start_text(" test foo",x,y);
os<<sw.frame_text("#00808080",x,y);
os<<CRLF;
MM_LOOP(ii,lc.m_layout)
{
const StrTy name=tt.node_name((*ii).first);
const Pos & p=(*ii).second;
ow(x,y,p.x,p.y);
//const IdxTy x=xpitch*(*ii).second.x+100;
//const IdxTy y=ypitch*(*ii).second.y+100;
os<<sw.text_text(name,x,y,w,h);
os<<CRLF;

}
MM_LOOP(ii,lc.m_parent)
{ // the index is the node and the value is the parent 
const Pos &  p0=lc.m_layout[(*ii).first];
const Pos &  p1=lc.m_layout[(*ii).second];
IdxTy x0,x1,y0,y1;
ow(x0,y0,p0.x,p0.y);
ow(x1,y1,p1.x,p1.y);
//const IdxTy x0=xpitch*p0.x+100;
//const IdxTy x1=xpitch*p1.x+100;
//const IdxTy y0=ypitch*p0.y+100;
//const IdxTy y1=ypitch*p1.y+100;
os<<sw.line_text(x0,y0,x1,y1);
os<<CRLF;
}

os<<sw.end_text();
os<<CRLF;
}


if (false)
{
os<<sw.start_text(" test foo",100,100);
os<<CRLF;
os<<sw.frame_text("#00808080",400,400);
os<<CRLF;
os<<sw.text_text(" foo doo do ",10,10,20,20);
os<<CRLF;
os<<sw.end_text();
os<<CRLF;
}


} // cmd_write_svg


void cmd_ncurses(Cip & cip , LocalVar & lv ) 
{ mjm_ncurses_links<> foo; foo.main(); }

void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).second))

} 

}

void cmd_load_tax_nodes(Cip & cip , LocalVar & lv ) 
{ m_tax_tree.read_ncbi_nodes_dmp(cip.p1); } 


void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("ncurses")]=&Myt::cmd_ncurses;
m_cmd_map[StrTy("load-tax-nodes")]=&Myt::cmd_load_tax_nodes;
m_cmd_map[StrTy("write-svg")]=&Myt::cmd_write_svg;
m_cmd_map[StrTy("navigate")]=&Myt::cmd_navigate;
m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;
m_cmd_map[StrTy("tt-bin")]=&Myt::cmd_tt_bin;


}
bool tax_tree_commands(CommandInterpretter & li, const StrTy & cmd,TaxTree & tt,const StrTy & p1 ,const StrTy & p2)
{
do {
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; } 
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; } 
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; } 
if (cmd=="save-tax") { tt.write_composite(p1);  continue; } 
if (cmd=="load-tax") { tt.clear(); tt.read_composite(p1);  continue; } 
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
if ( tax_tree_commands(li,cmd,m_tax_tree,p1,p2)) continue ;
//if (cmd=="load-tax-nodes") { m_tax_tree.read_ncbi_nodes_dmp(p1);  continue; } 
if (cmd=="load-tax-tree") { m_tax_tree.read_ncbi_dmp(p1);  continue; } 
if (cmd=="dump-tax-tree") { m_tax_tree.dump(p1);  continue; } 
if (cmd=="write-tax-single") { m_tax_tree.write_single(p1);  continue; } 
if (cmd=="save-tax") { m_tax_tree.write_composite(p1);  continue; } 
if (cmd=="load-tax") { m_tax_tree.read_composite(p1);  continue; } 
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


if (cmd=="lineage") 
{// local_label=m_tax_tree.node_info(myatoi(p1),myatoi(p2)); MM_MSG((local_label))   continue; } 
std::vector<StrTy> ww= li.words(); ww.erase(ww.begin());  
IdxTy taxon=m_tax_tree.best_indexed_taxon(ww,4); //  continue; 
std::vector<StrTy> lineage;
m_tax_tree.lineage(lineage,taxon);
Ss ss; 
MM_LOOP(ii,lineage)  { ss<<(*ii)<<" "; }
MM_MSG(ss.str())
continue;
}


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

void translate_tt1_to_tt2(TaxTree & tt1,TaxTree & tt2)
{
// go through each node in tt1 and find a good march in tt2
typedef std::map<IdxTy, IdxTy> Tree2Tree;
typedef std::vector<IdxTy> NodeList;
Tree2Tree t2t;
IdxTy good_ones=0;
IdxTy mult_ones=0;
IdxTy miss_ones=0;
IdxTy resol_ones=0;
MM_ERR(" translate_tt1_to_tt2")
// nned to move this out  this only makes rev not fwd 
mjm_index_lut splits=tt2.split_words();
MM_ERR(" translate_tt1_to_tt2 done word split ")

// for each node in the new tree, find equivalent node in the
// ncbi tree presumed in tt2
for(auto ii=tt1.begin_node();  ii!=tt1.end_node(); ++ii)
{
//const StrTy & name=(*ii).first;
const IdxTy&  sample_node=(*ii).first;
//const NodeList & snodes=(*ii).second;
const StrTy & name=(*ii).second;
MM_ERR(MMPR(name))
IdxTy vsz=0; // v.size();
// Try to find "name" in the ncbi tree tt2ZZ 
const auto jj=tt2.find_synonym(name);
 // this may be invalid 
const auto & v=(*jj).second;
if (jj!=tt2.end_synonym()){  vsz=v.size(); }
if (vsz==1) {++good_ones;  t2t[sample_node]=v[0]; continue; }
if (vsz>1)
{ // hopefully one lineage matches this one better lol. 
MM_DMEL("multiples",MMPR2(name,vsz))
++mult_ones;
}
if (vsz==0)
{ // go word by word, trying to find a unique match... 
MM_DMEL("nothingexact",MMPR(name))
++miss_ones;
const IdxTy namelen=name.length();
char c[namelen+1];
for (IdxTy i=0; i<namelen; ++i) 
if (name.c_str()[i]!='_' ) c[i]=name.c_str()[i]; else c[i]=' ';
c[namelen]=0;
StrTy name2(c);

Ss ss(name2);
std::vector<StrTy> pieces; 
std::vector<IdxTy> sizes; 
while (true) { StrTy y; ss>>y; if (y.length()==0) break; 
pieces.push_back(y);
} // true
IdxTy ones=0;
IdxTy sum=0;
IdxTy max=0;
MM_LOOP(jj,pieces) {
 const IdxTy hits=(splits((*jj)).size());
 sizes.push_back(hits);
if ( hits>max) max=hits;
sum+=hits;
if (hits==1) ++ones;
} // jj pieces 
if ( ones==1)
{
++resol_ones;
MM_LOOP(jj,pieces) {
 const IdxTy hits=(splits((*jj)).size());
if (hits==1) {t2t[sample_node]=(*jj)[0]; 
MM_DMEL("singeuniq",MMPR4(name,sum,max,(*jj)[0]))
break; } 
}
} // ones uniq
else
{
MM_DMEL("zedsizefail",MMPR4(name,sum,max,ones))

}
} // vsz==0
else
{
MM_DMEL("multisizefail",MMPR2(name,vsz))

}

} //ii 

 DumpDMel();
MM_ERR(MMPR4(good_ones,miss_ones,mult_ones,resol_ones))

} // translate_tt1_to_tt2


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

Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
TaxTree m_tax_tree;
TaxTrees m_tax_trees;

CounterMap m_cm;
CliTy m_cli;

}; //mjm_trees_and_tables



/////////////////////////////////////////////////////////

#ifdef  TEST_TREES_AND_TABLES__
int main(int argc,char **args)
{
typedef mjm_trees_and_tables  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

