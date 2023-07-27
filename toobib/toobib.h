/*

<one line to give the program's name and a brief idea of what it does.>

toobib - a collection of hacks to obtain bibtex descirbing a work
at a given url with or without advanced knowledge of the page's
contents based on domain or other clues. 
In its current form, it writes significant debugging information
to temporary files in the current directory. It is currently
only suggested for developers running in a sandbox or protected
area to avoid over writing other files on the user's computer.
It represents a collection of code but not a finished program
of any kind. 
Conceived and written by Mike Marchywka from 2019 to present.
See dates in individual code pieces as they were 
generated from my wizards. 
2023-07-13;
Copyright (C) <year> <name of author>

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@software{mmarchywka-toobib-2023-07-13,
  author = {Michael J Marchywka},
  title = {TooBib: an Impromptu Discovery tool for Bibtex Entries},
abstract={ A collection of hacks and empirical rules to search for bibtex given a url or possibly a local file},
institution={ Self },
license={ Knowledge sir should be free to all },
publisher={Mike Marchywka},
email={marchywka@hotmail.com},
authorid={orcid.org/0000-0001-9237-455X},
  filename = {toobib.h },
  url = {},
  version = {0.0.0},
  date-started = {2019-09-13},
  date-public = {2023-07-13}
  date = {2023-07-13}
}


*/


/*

toobib is best run and compiles with scripts. A recent working 
command line is as follows, assuming the headers are in the include paths,

g++ -std=gnu++11 -DTEST_toobib__ -gdwarf-3 -O0 -MMD -MF toobib.deps -I. -I/home/documents/cpp/mjm/hlib/ -I/home/documents/cpp/mjm/num -I/home/documents/cpp/mjm/  -I/home/documents/cpp/mjm/include -I/home/documents/cpp/pkg/include -I/home/documents/cpp/pkg -I /usr/include/libxml2 -Wall -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -Wno-misleading-indentation -Wno-unused-but-set-variable -x c++ ./toobib.h -o toobib.out -lreadline -lxml2 -lpthread
 2158  ./toobib.out quiet  -about -quit




*/

#ifndef MJM_toobib_H__
#define MJM_toobib_H__

//#define MIXED_CERR_STREAMS 


#include "mjm_globals.h"

#include "mjm_bibtex_find.h"
#include "mjm_med2bib_guesses.h"
// the dmel is not really needed or used but also include elaswhereh
#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"

#include "mjm_canned_methods.h"


#include "mjm_cli_ui.h"

/*
extern "C" {
#define TESTHTML_HEADER 
#include "junk/testHTML.c"
}
*/

#include "mjm_xml2_interface.h"
#include "mjm_tokenized_collections.h"
//#include "mjm_taxonomy_stuff.h"

#include "mjm_pawnoff.h"
#include "mjm_blob.h"
#include "mjm_misc_parse.h"
#include "mjm_bibtex_entry.h"
#include "mjm_bibtex_parse.h"
#include "mjm_ext_things.h"
#include "mjm_loo_parsing.h"

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
made from wizard script on Fri Sep 13 14:50:58 EDT 2019
toobib
*/


////////////////////////////////////////////////////////////////

class toobib_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
toobib_params( const StrTy & nm) : Super(nm) {}
toobib_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
// samples, 
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }

// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;
return ss.str();
}


}; // toobib_params


namespace toobib_traits
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
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // toobib_traits
///////////////////////////////////////////////////////////////




class mjm_toobib 
{
typedef toobib_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_toobib Myt;
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

typedef mjm_canned_methods Canned;
typedef toobib_params ParamGlob;
//typedef mjm_logic_base VariableStore;
typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;
typedef mjm_ragged_cursor RaggedCursor;



typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
//typedef mjm_tax_tree TaxTree;
//typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;

typedef mjm_pawnoff<Tr>  Handler;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
//typedef mjm_misc_parse<Tr> ParseTable;
typedef mjm_var_parse<Tr> ParseTable;
typedef typename ParseTable::Iv Iv;
typedef mjm_read_buffer<Tr> RdBuf;

typedef mjm_ext_things<Tr> ExtThings;
typedef mjm_bibtex_find<Tr> FindBib;

typedef temp_progress<Tr>  TempProg;


/*
not sure what this was lol, 
@www{$name,
  author =   {}, title =     {$title}, journal =    {}, year = {$year},
  volume =   {}, number =    {}, pages =     {}, note =  {},
urldate={$date},
url={$uin}
*/

typedef mjm_bibtex_entry<Tr> BibEntry;
typedef mjm_bibtex_entry_map<Tr> BibMap;
typedef std::map<StrTy,BibMap> BibMapMap;
typedef mjm_xml2_interface<Tr> MyXml;
typedef std::map<StrTy, MyXml> XmlMap;

typedef  mjm_med2bib_guesses <Tr>  Guesser;

typedef mjm_loo_parsing<Tr> Loo;


public:
// FIXME doh put this somwhere lol 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }

// API 
public :
mjm_toobib():m_dmel(new Dmel()) {Init();}
mjm_toobib(int argc,char **_args) : m_dmel(new Dmel())
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
if (i==istart) { MM_ERR(" did nothing with "<<args[i]) ++i;  } 

}
}
~mjm_toobib()
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
enum { LEGACY=1, DEBUG=2};
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-source") { ++i; command_modef(args[i]); ++i; }
if (s=="-cred") { ++i; m_guesser.set_chromate_launch_flags(1);  }
if (s=="-chromate") { ++i; m_guesser.set_chromate_launch_flags(myatoi(args[i]));  ++i; }
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-legacy") {m_mode|=(1<<LEGACY);  ++i; }
if (s=="-debug") {_loud(); _save_temps(); m_mode|=(1<<DEBUG);  ++i; }
if (s=="-devel") {_loud(); _save_temps();   ++i; }
if (s=="-quiet") {_quiet();  ++i; }
if (s=="quiet") {_quiet();  ++i; }
if (s=="-quit") { //  MM_ERR(" FK "<<MMPR4(i,argc,args[i],s)) 
 ++i; clean_up(); }
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
//{ CommandInterpretter li; li.set(cmd,1); command_mode(li); }
{ CommandInterpretter li; li.set(cmd,3); command_mode(li); }

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



bool flag_bit(const IdxTy flag, const IdxTy bit, const bool pol=true)
{
const bool x=((flag&(1<<bit))!=0);
return pol?x:!x;
}


void cmd_map_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy oname=cip.p1;
const StrTy iname=cip.p2;
const StrTy mname=cip.wif(3);
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
Ragged & rin=m_ragged_map[iname];
Ragged & rout=m_ragged_map[oname];
Ragged & rm=m_ragged_map[mname];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
map_ragged(rout,rin,rm);

MM_ERR("cmd_map_ragged "<<MMPR4(cmd,oname,iname,mname))
}

void cmd_bib_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy dest=cip.p1;
const StrTy rname=cip.p2;
//const StrTy mname=cip.wif(3);
//const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
MM_ERR("cmd_bib_ragged "<<MMPR4(cmd,dest,rname,flags))

Ragged & r=m_ragged_map[rname];
std::ofstream os(dest);
bib_ragged(os,r,flags);

}

#if 0

// mjm_med2bib_guesses commnands


if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="select") { x.select(cip.p1,atoi(cip.p2.c_str()));  } // Select
if (cmd=="skip") { x.select(cip.p1,1);  } // Select
if (cmd=="only") { x.select(cip.p1,2);  } // Select
if (cmd=="force") { x.select(cip.p1,4);  } // Select
if (cmd=="clip") {
MM_ERR(" in clip ")
MM_ERR(" done clip ")

}
if (cmd=="?") { 
MM_MSG(" quit dump add select skip only force guess clip ? h")

continue ;

#endif

IdxTy string_to_flag(const StrTy & s,const IdxTy t) const
{
const bool leg=Bit(m_mode,LEGACY);
const bool deb=Bit(m_mode,DEBUG);
if (t==3) // ref 
{

}
else if(t==2) // select
{
if  (s=="skip") return 1; 
else if (s=="only") return 2; 
else if (s=="force") return 4; 
} // t==2;
else if (t==1) // add
{
if (s=="legacy") return 1; // pmc stems 
else if (s=="force") return 12; // append anything  
else if (s=="paste") return 16; //   
else if (s=="")  
{ IdxTy rv=0;
if (leg) rv|=1;
if (!leg) rv|=16;
return rv;
}

} // t==1;
else if (t==0) // clip guess
{
//0 : show decisions  not passed to In
//1 : show trial not passed to In
//2 : keep going after first found or defer,  not paassed to in 
//3 : collect_all hits in each handler,   becomes Bit 2 in In flags 
//4 : paste citation to clipboard at end   becomes Bit 3 in In flags
if (s=="legacy") return  32+16+0+2; 
if (s=="all") return  8+4+2+1; 
if (s=="")  
{ IdxTy rv=0;
//if (leg) rv|= 32+16+2+1; 
// turn off PMC/PMID punting as this is a better converter
// than legacy code but if need speed can still punt 
if (leg) rv|= 16+2+0; 
if (deb) rv|= 8+4+2+1; 
return rv;
}
} // t==0
return myatoi(s);

}
void cmd_select_regex(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
const StrTy&  p2=cip.p2;
const IdxTy flags=string_to_flag(p2,2); // myatoi(p2);
Guesser & x= m_guesser;
 //x.clip(cip.p1.c_str(),strtod(cip.p2.c_str(),0)); MM_ERR(x.dump())
if (cmd=="select-regex") x.select(p1,flags); 
else if (cmd=="skip") x.select(p1,1); 
else if (cmd=="only") x.select(p1,2); 
else if (cmd=="force") x.select(p1,4); 
// MM_ERR(x.dump()) 
 } // cmd_select.

// format references for copy/paste into text 
void cmd_ref(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
const StrTy&  p2=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flags=string_to_flag(p2,3); // myatoi(p2);
Guesser & x= m_guesser;
 //x.clip(cip.p1.c_str(),strtod(cip.p2.c_str(),0)); MM_ERR(x.dump())
//if (cmd=="add") { x.add(cip.p1,cip.p2.c_str(),atoi(cip.wif(3).c_str()));  }
MM_MSG(MMPR4(cmd,p1,p2,flags)<<MMPR(x.flag_to_string(flags,3)))
x.ref(p1,flags); 
// MM_ERR(x.dump()) 
 } // cmd_add


//const bool stem_only=Bit(flags,0); // use old "non_pmc" prefix etc
//const bool show_dups=Bit(flags,1); // if there is a dup name show it 
//const bool append_dups=Bit(flags,2); // update file even if dup 
//const bool append_bad=Bit(flags,3); /

void cmd_add(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
const StrTy&  p2=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flags=string_to_flag(cip.wif(3),1); // myatoi(p2);
Guesser & x= m_guesser;
 //x.clip(cip.p1.c_str(),strtod(cip.p2.c_str(),0)); MM_ERR(x.dump())
//if (cmd=="add") { x.add(cip.p1,cip.p2.c_str(),atoi(cip.wif(3).c_str()));  }
MM_MSG(MMPR4(cmd,p1,p2,flags)<<MMPR(x.flag_to_string(flags,1)))
x.add(p1,p2.c_str(),flags); 
// MM_ERR(x.dump()) 
 } // cmd_add


void cmd_clip(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
const StrTy&  p2=cip.p2;
//const IdxTy flags=myatoi(p2);
const IdxTy flags=string_to_flag(p2,0); // myatoi(p2);
Guesser & x= m_guesser;
MM_MSG(MMPR4(cmd,p1,p2,flags)<<MMPR(x.flag_to_string(flags,0)))
 //x.clip(cip.p1.c_str(),strtod(cip.p2.c_str(),0)); MM_ERR(x.dump())
x.clip(p1,flags); 
// MM_ERR(x.dump()) 
 } // cmd_clip

/*
//enum { SHOW_DEC=0, SHOW_TRY=1, KEEP_GOING=2, GET_ALL=3};
const bool collect_all=Bit(flags,3);
const bool paste=Bit(flags,4);
// TODO once one caller sets these they were  set, use a temp  doh
IdxTy mif=m_in_flags;
if (collect_all) mif|=(1<<InpTy::COLLECT_ALL);
if (paste) mif|=(1<<InpTy::PASTE);
// punt_pmc similarly etc 
InpTy in(this,mif,flags);

//////////////////////////////
//enum { SHOW_DEC=0, SHOW_TRY=1, KEEP_GOING=2, GET_ALL=3};
const bool show_decisions=Bit(flags,0);
const bool show_trial=Bit(flags,1);
const bool keep_going=Bit(flags,2);
enum { SKP=0, ONE=1, COLLECT_ALL=2, PASTE=3, PUNT_PMC=4};


enum { SKP=0, ONE=1, COLLECT_ALL=2, PASTE=3, PUNT_PMC=4};

*/
// param 3 is m_rflags not used for anything. 
// InpTy in(this,m_in_flags,flags);
/*

Bit
0 : show decisions  not passed to In
1 : show trial not passed to In
2 : keep going after first found or defer,  not paassed to in 
3 : collect_all hits in each handler,   becomes Bit 2 in In flags 
4 : paste citation to clipboard at end   becomes Bit 3 in In flags
5 : punt PMC instead of collectings now, Bit 4 in In flags
legacy =32+16+2+!
*/


void cmd_guess(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
const StrTy&  p2=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flags=string_to_flag(cip.wif(3),0); // myatoi(p2);
Guesser & x= m_guesser;
MM_MSG(MMPR4(cmd,p1,p2,flags)<<MMPR(x.flag_to_string(flags,0)))
//if (cmd=="guess") { x.guess(cip.p1,cip.p2.c_str(),atoi(cip.wif(3).c_str())); MM_ERR(x.dump()) }
x.guess(p1,p2.c_str(),flags); 
// MM_ERR(x.dump()) 
 } // cmd_guess

void cmd_config_guesser(Cip & cip , LocalVar & lv ) 
{
const StrTy&  cmd=cip.cmd();
const StrTy&  p1=cip.p1;
//const StrTy&  p2=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
Guesser & x= m_guesser;
// may have been read with the code reader... 
const Ragged & r =m_ragged_map[p1];

//typedef mjm_ragged_cursor RaggedCursor;
RaggedCursor rc(r);
//if (cmd=="guess") { x.guess(cip.p1,cip.p2.c_str(),atoi(cip.wif(3).c_str())); MM_ERR(x.dump()) }
x.configure(rc); 
// MM_ERR(x.dump()) 
 } // cmd_config_guesser



////////////////////////////
//void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_add_to_fasta(cip ,  lv,  m_fasta_map  ); }

//void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_zymo_merge_fasta(cip ,  lv,  m_fasta_map, m_ragged_map  ); }

void cmd_link(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
//Ragged & r=m_ragged_map[name];
//Ragged & pr=m_ragged_map[prag];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
//const StrTy x=m_test_hand.get_output("","myclip -echo");
ExtThings et;
//const StrTy x=m_test_hand.get_clipboard();
const StrTy x=et.get_clipboard();
MM_ERR(MMPR2(cmd,x))
//ExtThings et;
const StrTy y=et.get_url(x,1);
MM_ERR(MMPR2(x,y))
}

// find bib(s) for the link or clipboard 
void cmd_find(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.p2);
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
//Ragged & r=m_ragged_map[name];
//Ragged & pr=m_ragged_map[prag];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
//const StrTy x=m_test_hand.get_output("","myclip -echo");
ExtThings et;
//const StrTy x=m_test_hand.get_clipboard();
const StrTy x=fn.length()?fn:et.get_clipboard();
if (flags!=0) { MM_ERR(MMPR2(x,flags)) } 
// TODO the finder uses the crossref xform which is bad ... 
const StrTy  bib=m_finder.find(x,flags);
//MM_ERR(MMPR2(x,bib))
//MM_ERR(MMPR2(cmd,x))
//ExtThings et;
//const StrTy y=et.get_url(x,1);
//MM_ERR(MMPR2(x,y))
}

void cmd_out_bib(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flags=myatoi(cip.p2);
//const StrTy prag=(cip.wif(4));
//Ragged & r=m_ragged_map[name];
//Ragged & pr=m_ragged_map[prag];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
//const StrTy x=m_test_hand.get_output("","myclip -echo");
//ExtThings et;
//const StrTy x=m_test_hand.get_clipboard();
if (fn=="-")
{
// TODO the finder uses the crossref xform which is bad ... 
const IdxTy x =m_finder.output(std::cout,flags);
}
else
{
std::ofstream fos(fn.c_str());
// TODO the finder uses the crossref xform which is bad ... 
const IdxTy x =m_finder.output(fos,flags);
}

//MM_ERR(MMPR2(x))
//MM_ERR(MMPR2(cmd,x))
//ExtThings et;
//const StrTy y=et.get_url(x,1);
//MM_ERR(MMPR2(x,y))
}




void cmd_read_blob(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy blob=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(cmd,blob,fn,flags))
m_blob_map[blob].load(fn);
}
void cmd_write_blob(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy blob=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(cmd,fn,blob,flags))
m_blob_map[blob].save(fn);
}


void cmd_gather_html(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy nx=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flagsh=myatoi(cip.wif(4));
const bool print_html = Bit(flags,0);
const bool print_html_names = Bit(flags,1);
MM_ERR(MMPR4(cmd,fn,nx,flags)<<MMPR3(flagsh,print_html,print_html_names))
MyXml& x=m_xml_map[nx];
typedef std::map<StrTy, StrTy> Map;
Map d;
x.html_kv(d,fn,0);
// this should be a map of all k-v pairs for each tag, sa meta tag
// but instead singles out say name and content for each pair
// need a more general thing.. 
MM_LOOP(ii,d) { MM_ERR(MMPR2((*ii).first,(*ii).second)) } 
//if (flags==0) { MM_ERR(" parse done but no output, flags==0") } 
//if (print_html) x.err_dump_html();
//if (print_html_names) x.print_html_names(flagsh);
}


void cmd_parse_html(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy nx=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const IdxTy flagsh=myatoi(cip.wif(4));
const bool print_html = Bit(flags,0);
const bool print_html_names = Bit(flags,1);
MM_ERR(MMPR4(cmd,fn,nx,flags)<<MMPR3(flagsh,print_html,print_html_names))
//htmlDocPtr doc = NULL;
Blob & bp =m_blob_map[fn];
const char * p=bp.ptr();
const IdxTy sz=bp.size();
MM_ERR(MMPR3(cmd,fn,sz))
MyXml& x= m_xml_map[nx];
x.parse_html(p,sz);
if (flags==0) { MM_ERR(" parse done but no output, flags==0") } 
if (print_html) x.err_dump_html();
if (print_html_names) x.print_html_names(flagsh);
}


void cmd_parse(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy fn=cip.p1;
if (fn.length()==0) fn=m_curr_bib_map;
else m_curr_bib_map=fn;
m_bib_map[m_curr_bib_map].parse(fn);
} // cmd_parse
void cmd_select(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy fn=cip.p1;
m_curr_bib_map=fn;
MM_ERR(MMPR2(cmd,m_curr_bib_map))
}

void cmd_map_dump(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
const IdxTy n=myatoi(cip.p1);
m_bib_map[m_curr_bib_map].dump(std::cout,n);
} // cmd_parse

void cmd_inter_fix(Cip & cip , LocalVar & lv ) 
{
//const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
//const IdxTy n=myatoi(cip.p1);
//m_bib_map.command_mode(cip.li());
m_bib_map[m_curr_bib_map].try_to_fix(cip.li());
} // cmd_parse
void cmd_save(Cip & cip , LocalVar & lv ) 
{
//const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
//m_bib_map.command_mode(cip.li());
m_bib_map[m_curr_bib_map].save(fn,flags);
MM_ERR(" saving with "<<MMPR2(fn,flags))
} // cmd_parse

void cmd_shell(Cip & cip , LocalVar & lv ) 
{
//const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const IdxTy flags=myatoi(cip.p2);
//m_bib_map.command_mode(cip.li());
//m_bib_map[m_curr_bib_map].save(fn,flags);
//MM_ERR(" saving with "<<MMPR2(fn,flags))
if (fn.length()==0) 
{
MM_ERR(" no shell command ")
return;
}
Blob d,out,err;
const StrTy cmd=cip.line().c_str()+2;
//IdxTy rc= m_test_hand.fileio(out,err,d,cmd);
IdxTy rc= m_test_hand.raw(cmd);
MM_ERR(MMPR4(StrTy(d),StrTy(out),StrTy(err),rc))

} // cmd_shell

void cmd_flags(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
//if (s=="-legacy") {m_mode|=(1<<LEGACY);  ++i; }
//if (s=="-debug") {m_mode|=(1<<DEBUG);  ++i; }
if (cmd=="legacy") { m_mode=(1<<LEGACY); }
else if (cmd=="debug") { m_mode=(1<<DEBUG); }
else MM_ERR(" logic error in cmd_flags "<<MMPR(cmd))
MM_ERR("cmd_flags "<<MMPR2(cmd,m_mode))
//const StrTy fn=cip.p1;
//const IdxTy flags=myatoi(cip.p2);
//m_bib_map.command_mode(cip.li());
//m_bib_map[m_curr_bib_map].save(fn,flags);
//MM_ERR(" saving with "<<MMPR2(fn,flags))
//if (fn.length()==0) 
//{
//MM_ERR(" no shell command ")
//return;
//}
//Blob d,out,err;
//const StrTy cmd=cip.line().c_str()+2;
//IdxTy rc= m_test_hand.fileio(out,err,d,cmd);
//IdxTy rc= m_test_hand.raw(cmd);
//MM_ERR(MMPR4(StrTy(d),StrTy(out),StrTy(err),rc))

} // cmd_flags




void cmd_check_comments(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
const IdxTy n=myatoi(cip.p1);
m_bib_map[m_curr_bib_map].check_comments(std::cout,n);
} // cmd_parse
void cmd_bib_cmd(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
//const IdxTy n=myatoi(cip.p1);
m_bib_map[m_curr_bib_map].cmd(cip.p1,cip.p2);
} // cmd_parse



void cmd_missing(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy f=cip.p1;
const IdxTy n=myatoi(cip.p2);
m_bib_map[m_curr_bib_map].missing(std::cout,f,n);
} // cmd_parse




void cmd_bib(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
//Ragged & r=m_ragged_map[name];
//Ragged & pr=m_ragged_map[prag];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
Bib();
MM_ERR(MMPR(m_bib.format()))
}

void cmd_paste(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
//Ragged & r=m_ragged_map[name];
//Ragged & pr=m_ragged_map[prag];
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
PasteCitation();
MM_ERR(MMPR(m_bib.format()))
}






void cmd_write_svg_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const StrTy prag=(cip.wif(4));
Ragged & r=m_ragged_map[name];
Ragged & pr=m_ragged_map[prag];
MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))

}


void cmd_transpose_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ; } // transpose
void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_dump_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_dump_ragged(std::cout,  cip ,  lv, m_ragged_map  ) ; }


void cmd_add_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ; }
//void cmd_transpose_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_transpose_ragged(cip ,  lv, m_tokragged_map  ) ; } // transpose
//void cmd_read_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_read_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_write_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_write_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_add_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_add_ragged( cip ,  lv, m_tokragged_map  ) ; }


//void cmd_tt(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }


void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_quit(Cip & cip , LocalVar & lv )
{ clean_up(); return;  }
void cmd_quiet(Cip & cip , LocalVar & lv ) { _quiet(); return;  }
void cmd_loud(Cip & cip , LocalVar & lv ) {_loud();  return;  }
void _loud() { mjm_global_flags::mm_err_enable=true; }
void _quiet() { mjm_global_flags::mm_err_enable=!true; }

void _save_temps() { mjm_global_flags::mm_delete_temps=!true;  }
void cmd_save_temps(Cip & cip , LocalVar & lv )
{ 
_save_temps();
return;  }




void cmd_cm(Cip & cip , LocalVar & lv )
{ dump_cm(); return;  }
void cmd_banner(Cip & cip , LocalVar & lv )
{ config_banner(); return;  }

void cmd_set_param(Cip & cip , LocalVar & lv )
{  //kkkkkkkkkkk 
//const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
 if (cip.li().cmd_ok(3))  m_flp.set(cip.li().cmd_set());
        return;  }
void cmd_get_param(Cip & cip , LocalVar & lv )
{
 if (cip.li().cmd_ok(2))
        std::cout<<lv["local_label"]<<" "<<cip.li().word(1)<<" "
                        <<m_flp.get(cip.li().word(1))<<CRLF;
         }





void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{ // TODO it should be possible to call the f ptr and get a msg back 
MM_ERR(MMPR((*ii).first))

} 

} // help 

void cmd_list(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_bib_map) { MM_MSG("m_bib_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tokragged_map) { MM_MSG("m_tokragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_blob_map) { MM_MSG("m_blob_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_xml_map) { MM_MSG("m_xml_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_dig_map) { MM_MSG("m_dig_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(m_guesser.dump())
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
m_cmd_map[StrTy("skip")]=&Myt::cmd_select_regex;
m_cmd_map[StrTy("only")]=&Myt::cmd_select_regex;
m_cmd_map[StrTy("force")]=&Myt::cmd_select_regex;
m_cmd_map[StrTy("select-regex")]=&Myt::cmd_select_regex;
m_cmd_map[StrTy("ref")]=&Myt::cmd_ref;
m_cmd_map[StrTy("add")]=&Myt::cmd_add;
m_cmd_map[StrTy("clip")]=&Myt::cmd_clip;
m_cmd_map[StrTy("guess")]=&Myt::cmd_guess;
m_cmd_map[StrTy("configure-guesser")]=&Myt::cmd_config_guesser;
m_cmd_map[StrTy("link")]=&Myt::cmd_link;
m_cmd_map[StrTy("find")]=&Myt::cmd_find;
m_cmd_map[StrTy("out-bib")]=&Myt::cmd_out_bib;
m_cmd_map[StrTy("bib")]=&Myt::cmd_bib;
m_cmd_map[StrTy("parse")]=&Myt::cmd_parse;
m_cmd_map[StrTy("select")]=&Myt::cmd_select;
m_cmd_map[StrTy("parse-html")]=&Myt::cmd_parse_html;
m_cmd_map[StrTy("gather-html")]=&Myt::cmd_gather_html;
m_cmd_map[StrTy("read-blob")]=&Myt::cmd_read_blob;
m_cmd_map[StrTy("write-blob")]=&Myt::cmd_write_blob;
m_cmd_map[StrTy("map-dump")]=&Myt::cmd_map_dump;
m_cmd_map[StrTy("check-comments")]=&Myt::cmd_check_comments;
m_cmd_map[StrTy("bib-cmd")]=&Myt::cmd_bib_cmd;
m_cmd_map[StrTy("inter-fix")]=&Myt::cmd_inter_fix;
m_cmd_map[StrTy("missing")]=&Myt::cmd_missing;
m_cmd_map[StrTy("paste")]=&Myt::cmd_paste;
m_cmd_map[StrTy("save")]=&Myt::cmd_save;
m_cmd_map[StrTy("!")]=&Myt::cmd_shell;
m_cmd_map[StrTy("legacy")]=&Myt::cmd_flags;
m_cmd_map[StrTy("debug")]=&Myt::cmd_flags;
m_cmd_map[StrTy("all")]=&Myt::cmd_flags;

m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("quit")]=&Myt::cmd_quit;
m_cmd_map[StrTy("quiet")]=&Myt::cmd_quiet;
m_cmd_map[StrTy("loud")]=&Myt::cmd_loud;
m_cmd_map[StrTy("save-temp")]=&Myt::cmd_save_temps;
m_cmd_map[StrTy("cm")]=&Myt::cmd_cm;
m_cmd_map[StrTy("banner")]=&Myt::cmd_banner;
m_cmd_map[StrTy("set-param")]=&Myt::cmd_set_param;
m_cmd_map[StrTy("get-param")]=&Myt::cmd_get_param;

// needs a nwords parameter
m_cmd_map[StrTy("read-code")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("dump-ragged")]=&Myt::cmd_dump_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;

//m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
//m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
//m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
//m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("map-ragged")]=&Myt::cmd_map_ragged;
m_cmd_map[StrTy("bib-ragged")]=&Myt::cmd_bib_ragged;
//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
//m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
//m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
//m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
//m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
//m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
//m_cmd_map[StrTy("group-stats")]=&Myt::cmd_group_stats;

//m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;

//m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
//m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
//TaxTree & tt = m_tax_tree;
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
mloc["local_label"]=local_label;
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
if ( m_done) return;
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

MM_MSG_CRLF

} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;
// this ends with cmd prompt and no crlf putting new prompt on same line
MM_MSG_CRLF

} 


void about()
{
Ss ss;
ss<<" toobib built  "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com started  Fri Sep 13 14:50:58 EDT 2019 "<<CRLF;
ss<<" Primarily now a replacement for med2bib a script for bibtex scraping "<<CRLF;
ss<<" Given a url or filename on clipboard, find possible bibtex entries  "<<CRLF;
ss<<" for the content featured on web page or in file   "<<CRLF;
ss<<" Empty skeleton code add credits like this,   "<<CRLF;
ss<<" Yet another program to go TO Organized BIBtex "<<CRLF;
ss<<" Originally intended to construct bibtex from webpages  "<<CRLF;
ss<<" lacking it based on scraping but expanded   "<<CRLF;
ss<<" to fix script based bibtex collections   "<<CRLF;
ss<<" With some input from texhax@tug.org "<<CRLF;
ss<<" Uses external libraries for xml and json parsing as well as many   "<<CRLF;
ss<<" text and web getting tools. readline for user interface    "<<CRLF;

ss<<" html parser  downloaded from https://dev.w3.org/XML/testHTML.c "<<CRLF;
ss<<" html parser   Daniel.Veillard@w3.org " <<CRLF;
ss<<" json parser : https://rapidjson.org/ "<<CRLF;
ss<<" links to -lreadline -lxml2 -lpthread "<<CRLF;
ss<<" Web interface : wget and headless chrome "<<CRLF;
ss<<" Nodejs wscat until websocket lib integrated  "<<CRLF;
ss<<" Generally specific sites like crossref are credited when used  "<<CRLF;
ss<<" from header files :"<<CRLF;
ss<<mjm_global_credits::about();

std::ostream & os=std::cout;
os<<ss.str();

}
#if 0 
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

#endif


// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fck 
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

#if 1 
StrTy  bib_name_ragged( const Ragged & r,const IdxTy flags  )
{
TempProg tp;
BibEntry _bib;
_bib.clear();
IdxTy sz=r.size();
if ( sz==0 ) return StrTy(); 
for(IdxTy i=0; i<sz;++i)
{
Ragged::Line l=r.line(i);
const IdxTy szz=l.size();
if (szz==0) continue;
StrTy v="";
if (szz>1) v=l[1];
for(IdxTy j=2; j<szz; ++j) v=v+" "+l[j]; 
//MM_MSG(" add "<<MMPR2(l[0],v))
_bib.set(l[0],v);
}
_bib.make_name(tp);
return _bib.name();
/*
m_bib.set("url",link);
const StrTy date=RemoveCRLF(m_test_hand.today());
m_bib.set("urldate",date);
m_bib.comment(date);
m_bib.comment("autogenerated by toobib " );
Blob & bt=tp.blob("title");
FindHtmlTitle(bt,bp);
m_bib.set("title",RemoveCRLF(StrTy(bt)));
//std::map<StrTy,StrTy> m;
FindMeta(tp.meta(),bp);
MM_LOOP(ii,tp.meta()) { MM_ERR(MMPR2((*ii).first,(*ii).second)) } 
m_bib.make_name(tp);
Blob & bhtml=tp.blob("rendered");
et.get_url(bhtml,link,1);
//et.render_html(bhtml,bp,link);
*/

}

#endif


template <class Os> 
void bib_ragged(Os & os, const Ragged & r,const IdxTy flags  )
{
const bool last_name=Bit(flags,0);
const bool no_name=Bit(flags,1);
const bool indent=Bit(flags,2);
const bool make_name=Bit(flags,3);
IdxTy sz=r.size();
if ( sz==0 ) return; 
const StrTy leader=indent?"    ":"";
IdxTy i=sz-1;
IdxTy start=0;
IdxTy end=sz;
StrTy  newname= bib_name_ragged(  r,0 );
//MM_MSG(MMPR(newname))
if (!no_name)
{
if (last_name){  i=sz-1; end=i; } 
else { i=0; start=1;}
Ragged::Line l=r.line(i);
if (!make_name)
{
if (l.size()>0) { os<<"@"<<l[0]<<" { "; }
if (l.size()>1) { os<<l[1]<<","; }
if (l.size()>0) { os<<CRLF; }
}
if (make_name)
{
if (l.size()>0) { os<<"@"<<l[0]<<" { "; }
if (l.size()>0) { os<<newname<<","; }
if (l.size()>0) { os<<CRLF; }
} // make_name
} // no_name

for (i=start; i<end; ++i)
{
Ragged::Line l=r.line(i);
const IdxTy szz=l.size();
if (szz==0) continue;
StrTy v="";
if (szz>1) v=l[1];
for(IdxTy j=2; j<szz; ++j) v=v+" "+l[j]; 
os<<leader<<l[0]<<"={"<<v<<"}";
if (i!=(end-1)) os<<",";
os<<CRLF;
} 
os<<"}"<<CRLF;

}
void map_ragged(Ragged & rout,  const Ragged & rin, const Ragged & rm )
{
ReadWriteMap rwm;
rm.to_map(rwm,0);
const auto & m = rwm.map();
MM_LOOP(ii,rin)
{
Ragged::Line l=(*ii);
if ( l.size()>0)
{
MM_ERR(" reading line "<<MMPR2(l.size(),l[0]))
StrTy & k=l[0];
const auto ii=m.find(k);
if (ii!=m.end()) k=(*ii).second; 
}
rout.add(l);
} // ii 


} // map_ragged



/////////////////////////////////////////////////////////////////////
private:
static bool Bit(const IdxTy f, const IdxTy b) { return ((f>>b)&1)!=0; } 
void Init()
{
m_done=false;
m_mode=0;
}








IdxTy FindBlob(const char * p, const IdxTy sz, const IdxTy first, const char * s)
{
const IdxTy slen=strlen(s);
for (IdxTy i=first; i<sz; ++i)
{
if (strncmp(p+i,s,slen)==0) return i;
}
return ~0;
}
//<meta name="SERVER" content="PDSVCMS01" />
template <class Tkv> 
void FindMeta(Tkv & m, const Blob & b)
{
RdBuf x(10); // make 1k default block size for now 
x.start_new();
IdxTy kv=0;
const char * p=b.ptr();
const IdxTy sz=b.size();
IdxTy start=0;
IdxTy end=0;
while (start!= ~0) { 
start= FindBlob(p,sz,end,"<meta name=");
if (start<sz) start=start+12;
end=FindBlob(p,sz,start,"\"");
if(start<sz) if ( end<sz)if (end>=start)
{  x.append(p+start,end-start);
x.start_new();
start= FindBlob(p,sz,end,"content=");
if (start<sz) start=start+9;
end=FindBlob(p,sz,start,"\"");
if(start<sz) if ( end<sz)if (end>=start)
{  x.append(p+start,end-start);
x.start_new();
kv+=2;
 continue ; }}
break;

} // while
for (IdxTy i=0; i<kv; i+=2 ) m[x.string(i+1)]=x.string(i+2); 
}


void FindHtmlTitle(Blob & bt, const Blob & bp)
{
RdBuf x(10); // make 1k default block size for now 
const char * p=bp.ptr();
const IdxTy sz=bp.size();
IdxTy start=FindBlob(p,sz,0,"<title>");
if (start<sz) start=start+7;
IdxTy end=FindBlob(p,sz,start,"</title>");;
if(start<sz) if ( end<sz)if (end>=start)  x.append(p+start,end-start);
//x.load(fn);
//Release();
x.append(char(0));
bt.adopt(x.ptr(),x.size(),x.alloc());
x.release();

}

void Bib()
{
TempProg&  tp = m_tp;
ExtThings et;
//const StrTy link= et.get_clipboard();
Blob & b=tp.blob("uin");
et.get_clipboard(b);
const StrTy link= StrTy(b); // tp.et.get_clipboard();
//MM_ERR(MMPR2(cmd,x))
Blob & bp=tp.blob("page");
//const StrTy y=et.get_url(link);
et.get_url(bp,link,0);
const StrTy y=StrTy(bp); // et.get_url(link);
//MM_ERR(MMPR2(x,y))
m_bib.clear();
m_bib.type("www");
m_bib.set("url",link);
const StrTy date=RemoveCRLF(m_test_hand.today());
m_bib.set("urldate",date);
m_bib.comment(date);
m_bib.comment("autogenerated by toobib " );
Blob & bt=tp.blob("title");
FindHtmlTitle(bt,bp);
m_bib.set("title",RemoveCRLF(StrTy(bt)));
//std::map<StrTy,StrTy> m;
FindMeta(tp.meta(),bp);
MM_LOOP(ii,tp.meta()) { MM_ERR(MMPR2((*ii).first,(*ii).second)) } 
m_bib.make_name(tp);
Blob & bhtml=tp.blob("rendered");
et.get_url(bhtml,link,1);
//et.render_html(bhtml,bp,link);
//MM_ERR(MMPR(StrTy(bhtml)))
}


StrTy RemoveCRLF( const StrTy s)
{
const IdxTy len=s.length();
char c[len+2];
IdxTy ptr=0;
const IdxTy mask= ParseTable::EOL; //  m_clut
for(IdxTy i=0; i<len; ++i)
{
const char cp=s.c_str()[i];
IdxTy lu= m_clut.lut(cp);
if ((lu&mask)==0) { c[ptr]=cp; ++ptr; } 

}
c[ptr]=0;
return StrTy(c); 
}

void Paste(const StrTy & x)
{
ExtThings et;
et.put_clipboard(x); 
}

void Paste()
{
ExtThings et;
//et.put_clipboard(Citation(m_bib.name())); 
et.put_clipboard(m_bib.name()); 
}

void PasteCitation()
{
Paste(Citation(m_bib.name())); 

}

StrTy Citation(const StrTy & n) 
{ return StrTy("\\cite{")+n+StrTy("}"); } 


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


// MEMBERS members 

bool m_done;
IdxTy m_mode;
ParamGlob m_flp;
//VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
TokRaggedMap m_tokragged_map;
//DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
//CharMat m_char_mat;
//MiiMap m_luts;
//TaxTree m_tax_tree; // now have sequences ID's to taxon 
//TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;
// not used except for test
Handler m_test_hand;
BlobMap m_blob_map;
BibEntry m_bib;
ParseTable m_clut;
TempProg  m_tp; // probably obsolete 
BibMapMap m_bib_map;
StrTy m_curr_bib_map;
 XmlMap m_xml_map;
// TODO the finder uses the crossref xform which is bad ... 
FindBib m_finder;

Guesser m_guesser;


}; //mjm_toobib 

/////////////////////////////////////////////////////////

#ifdef  TEST_toobib__
int main(int argc,char **args)
{
typedef mjm_toobib Myt;
///typedef double D;
//typedef unsigned int IdxTy;
if (strcmp(args[0],"-quiet")==0) { mjm_global_flags::mm_err_enable=!true; }
if (strcmp(args[0],"quiet")==0) { mjm_global_flags::mm_err_enable=!true; }
if (argc) { 
if (strcmp(args[1],"-quiet")==0) { mjm_global_flags::mm_err_enable=!true; }
if (strcmp(args[1],"quiet")==0) { mjm_global_flags::mm_err_enable=!true; }
}

Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif
