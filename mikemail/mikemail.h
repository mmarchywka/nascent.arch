
#ifndef MJM_mikemail_H__
#define MJM_mikemail_H__

// moving this down turns timeout() into stdscr or something from ncurses? 
//#include "MailCore/MailCore.h"
#include <dirent.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <typeinfo>
#include <stdio.h>
#include <string>

#include "mjm_globals.h"
//#define HAVE_MIKE_MUTT_SERVER
#ifdef HAVE_MIKE_MUTT_SERVER
#warning including HAVE_MIKE_MUTT_SERVER code 
#include "mjm_mutt_server.h"
#include "mjm_neomutt_interface.h"
//#undef QUOTE
#endif

#include <mjm_thread_util.h>
#include <mjm_line_istreams.h>
#include <mjm_ascii_tree.h>
// must come before mjm_globalsh
// mor now moved to build file as flag 
//#define LINE_LOCK_MM_ERR

#include <mjm_libetpan_readmsg.h>


#include "mjm_pawnoff.h"
#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_canned_methods.h"
#include "mjm_cli_ui.h"
#include "mjm_tokenized_collections.h"

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

// these do not have their own includes and need to be down here 
#include <mjm_generic_message.h>
#include <mjm_message_store.h>
#include <mjm_mbox_formats.h>
#include <mjm_paranoid_ofstream.h>
#include <mjm_mailproc/mjm_mailproc.h>

/*


1)Download files from IMAP server for archival and analysis. 
2) Mail processing,
for example, 
 cat mboxdiet.txt
load-mbox foo /var/mail/marchywka
#view-src foo 1 3 Received:
#mmp-proc mproc foo
read-ragged foorag diet_rules.txt
list
mmp-config mproc foorag
mmp-config mproc foorag 1
#list
mmp-proc mproc foo

cat diet_rules.txt 
handler default desc crap script "cat -" dump 1 
handler diet desc "dietformhandler" noun-file diet/canon_nouns.txt ignores-file diet/ignores.txt reserveds-file diet/reserveds.txt
# ls diet canon_nouns.txt  diary.txt  ignores.txt  reserveds.txt
handler default desc doo script ls dump 1 
handler default desc moo script df dump 1 
rule n 0 str marchywka location from: code 2
rule n 0 str diet location subject: code 1 
pair 0 0
pair 1 1 
pair 0 0
pair 0 0
vars foo doo dee duu
vars doh rey me fa

=============================================
=============================================
Alternatively, for the yahoo or an imap server using the neomutt
library,  may work with libetpan too not sure code is there,

./hist_doing: 1877  ./mikemail.out -source muttdiet.txt -quit 
./hist_doing: 1882  vi muttdiet.txt 
./hist_doing: 1883  cp muttdiet.txt muttdiet_spam.txt 
./hist_doing: 1953  ./mikemail.out -source muttdiet.txt -quit 
./hist_doing: 1756  cat muttdiet.txt 
marchywka@happy:/home/documents/cpp/proj/mikemail$ cat muttdiet.txt
read-ragged foorag diet_rules.txt
mmp-config mproc foorag
mmp-config mproc foorag 1
mutt-proc proc_yahoo 1 mproc




*/



////////////////////////////////////////////////////////////////

class mikemail_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
public:
mikemail_params( const StrTy & nm) : Super(nm) {}
mikemail_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
//IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
StrTy mail_clip_file() const { return m_map.get_string("mail_clip_file","mail_clip_file.txt"); }
//IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
// samples, 
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }
bool buffer_msg() const { return m_map.get_bool("buffer_msg",!true); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;
ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"ragged_params="<<ragged_params()<<sep;
ss<<"buffer_msg="<<buffer_msg()<<sep;
ss<<"mail_clip_file="<<mail_clip_file()<<sep;
return ss.str();
}
}; // mikemail_params


namespace mikemail_traits
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
typedef uint64_t KeyCode;
enum { BAD=~0};
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef std::vector<StrTy> Words;

typedef mjm_libetpan_readmsg<Tr> LibetpanReader;
typedef mjm_folder_downloader<Tr> FolderDownloader;
 

}; // mikemail_traits
///////////////////////////////////////////////////////////////

class mjm_mikemail 
{
typedef mikemail_traits::Tr  Tr;
typedef mjm_mikemail Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
enum {BAD=Tr::BAD};

typedef mjm_thread_util<Tr>::mutex_vector MutexVector;

typedef mikemail_traits::LibetpanReader LibetpanReader;
typedef mikemail_traits::FolderDownloader FolderDownloader;
typedef mjm_generic_message<Tr>  MyMsg;
typedef MyMsg::alternative_ranks AltRank;
typedef std::vector<MyMsg> MyMsgVector;

typedef mjm_message_store<Tr,MyMsg> MyMsgMap;
typedef std::map<StrTy, MyMsgMap> MsgMapMap;

typedef mjm_mailproc<Tr>  MailProc;
typedef mjm_mailproc_map<Tr>  MailProcMap;
typedef MailProc::proc_io_type proc_io_type;

typedef std::vector<IdxTy> MyMsgSortOrder;
typedef mjm_pawnoff<Tr> Hand;
typedef mjm_paranoid_ofstream<Tr> SafeOsTy;
typedef mjm_mbox_formats<Tr> FormatterTy;
typedef mjm_ascii_tree<Tr> TreeTy;

typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

typedef mjm_canned_methods Canned;
typedef mikemail_params ParamGlob;
typedef mjm_logic_base VariableStore;
typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;

typedef std::map<StrTy, LibetpanReader> LibetMap;
typedef std::map<StrTy, FolderDownloader> DownloaderMap;

#ifdef HAVE_MIKE_MUTT_SERVER
typedef mjm_neomutt_interface<Tr> MuttInterface; 
#endif

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

class mbox_cursor
{
public:
void inc() { m_n=m_n+1; }
IdxTy ptr() { return m_n; }
IdxTy ptr(const IdxTy p) { m_n=p; return m_n; }
IdxTy m_n;
}; 
typedef std::map<StrTy, mbox_cursor> MboxCursorMap;
// MboxCursorMap m_cursors;

class message_view
{
public:
message_view() {Init();}
~message_view() {}
bool Bit(const IdxTy f, const IdxTy b) { return (f&(1<<b))!=0; } 
const Words & header_fields() const { return m_hf; } 
const Words & tree_fields() const { return m_tf; } 
const Words & header_lines() const { return m_hl; } 
const IdxTy mode() const { return m_mode;} 
void  mode(const IdxTy m )  {  m_mode=m ;} 
AltRank & alt_prefs() const  { return m_alts; } 
private:
void Init()
{
{
Words& w= m_hf;
w.push_back("From:");
w.push_back("Date:");
w.push_back("Content-Length:");
w.push_back("Content-Type:");
//w.push_back("boundary");
w.push_back("mime");
w.push_back("Subject:");
}
Words& w= m_hl;
w.push_back("From:");
w.push_back("To:");
w.push_back("Subject:");

Words&  v=m_tf;
v.push_back("mime");
v.push_back("Content-Type:");
v.push_back("Content-Disposition:");
v.push_back("Content-Transfer-Encoding:");

m_alts["text/plain"]=1;
m_alts["text/html"]=2;


m_mode=2;
}
Words m_hf,m_hl; // header fields 
Words m_tf;  // tree fields 
mutable AltRank m_alts;
IdxTy m_mode;

}; // message_view

typedef message_view  ViewParam ;



//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
// not currently used, intended for sorting later. 
typedef mjm_tax_tree TaxTree;
typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

//typedef std::vector<IdxTy> LocTy;

public:
// FIXME doh put this somwhere lol  - canned is a big junk bin 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }
bool Bit(const IdxTy f, const IdxTy b) { return ((f&(1<<b))!=0); } 
public :
mjm_mikemail():m_dmel(new Dmel()) {Init();}
mjm_mikemail(int argc,char **_args) : m_dmel(new Dmel())
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
cmdlcmd( i, argc, args);
if (i==istart) { MM_ERR(" did nothing with "<<args[i]) ++i;  } 

}
}
~mjm_mikemail()
{
//pthread_mutex_destroy(mutex1());
//pthread_mutex_destroy(mutex2());
delete m_dmel;
}
////////////////////////////////////////////////////////
// command block
/*
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
*/
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-source") { ++i; command_modef(args[i]); ++i; }
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { MM_ERR(" quitting "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
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
// these are not really used here. 
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
////////////////////////////////////////////////////
// real code 
void load_mbox(MyMsgMap & mmm,const StrTy &fn)
{
std::ifstream ifs(fn);
//typedef mjm_line_istream<Tr> Ty;
#if 0 
typedef mjm_loaded_istream<Tr> Ty;
Ty ifs(fn.c_str());
ifs.launch();
#endif
LineIterator li(&ifs);
li.set_split(1,' ');
li.remove_crlf(true);
li.use_alt_eol('\n',true);
IdxTy i=0; 
while ( ifs.good()&&!ifs.eof())
{
//MyMsg m=m_formatter.next_mbox_msg(ifs);
// this is where the stupid assignment is, not
// sure wat it returns but Assign was wrong so I know
// it was  using it... 
MyMsg m=m_formatter.next_mbox_msg(li);
MM_STATUS(MMPR2(m.verbatim().length(),i)<<"        ")
if (m.verbatim().length()!=0) mmm.append(m);
++i;
}
MM_ERR(" loaded mbox messages "<<MMPR3(i,ifs.good(),ifs.eof())<<"    " )
// cross fingers ifs goes out of scope well lol 
} // load_mbox

void save_attachments(MyMsgMap & mmm)
{
const IdxTy sz=mmm.size();
MM_ERR(" saving attachments "<<MMPR(sz))
for (IdxTy i=0; i<sz; ++i)
{
// non const as subparts will populate 
MyMsg & msg=mmm[i];
m_formatter.decode_parts(msg);
} // i 

} // save_attachments


//const StrTy & mime() const { return m_mime; }
//const IdxTy parts() const { return m_parts.size(); }
//const IdxTy depth() const { return m_depth; }
//const IdxTy depth(const IdxTy i ) { m_depth=i;  return m_depth; }
// ref should return as it is just a ptr 
//Myt& part(const IdxTy i )  { return *m_parts[i]; }

//Content-Type: image/jpeg; name="20190409_205520.jpg"^M
//Content-Disposition: attachment; filename="20190409_205520.jpg"^M
//Content-Transfer-Encoding: base64^M



typedef mjm_ascii_tree<Tr> AsciiTree;
typedef  AsciiTree::Node AtNode;
//AsciiTree at;
void  att_tree(AsciiTree & at, MyMsg & m, const  Words & w, AtNode * p )
{
AtNode * pn=0;
StrTy name=StrTy("->");// +m.hval_lc("Content-Type:");
// name+=m.hval_lc("Content-Disposition:");
// name+=m.hval_lc("Content-Transfer-Encoding:");
MM_LOOP(ii,w) { name+=m.hval_lc(*ii); } 
//MM_ERR(MMPR(name))
if (p==0) pn=&at.add_root(name);
else pn=&at.add_kid(*p,name);

const IdxTy sz=m.parts();
for(IdxTy i=0; i<sz; ++i)
{
att_tree(at,m.part(i),w,pn);
}

} // att_tree

StrTy tree_attachments( MyMsg & msg, const Words & w)
{
AsciiTree at;
att_tree(at,msg,w,0);
return at.ascii();
}

void  tree_attachments(MyMsgMap & mmm)
{
Words w;
w.push_back("From:");
w.push_back("Date:");
w.push_back("Content-Length:");
w.push_back("Content-Type:");
w.push_back("boundary");
w.push_back("Subject:");

Words v;
v.push_back("Content-Type:");
v.push_back("Content-Disposition:");
v.push_back("Content-Transfer-Encoding:");

const IdxTy sz=mmm.size();
for (IdxTy i=0; i<sz; ++i)
{
// non const as subparts will populate 
MyMsg & msg=mmm[i];
const StrTy h= msg.header_string(w);
// this does not do anything that the ctor did not do 
//m_formatter.decode_parts(msg);
//AsciiTree at;
//att_tree(at,msg,0);
const StrTy x=tree_attachments(msg,v); //  at.ascii();
//if ( msg.parts()>0) 
//if ( Bit(flags,0)==0) 
{ MM_ERR(MMPR2(i,h)); std::cerr<<x<<CRLF;  } 
// m_formatter.decode_parts(msg);
} // i 
//return at;
} // tree_attachments

void select_range(std::vector<IdxTy> & selvec,const StrTy & fn, mbox_cursor & cu)
{
const IdxTy sz=fn.length();
const char * px=fn.c_str();
char pa[sz+1];
if (px[0]=='+')  { cu.inc(); selvec.push_back(cu.ptr()); return ; } 
memcpy(pa,px,sz+1);
char * p=&pa[0];
IdxTy start=0;
IdxTy prior=BAD;
bool r=false;
for(IdxTy i=0; i<=sz; ++i)
{
const char c=p[i];
if ((c==0)||(c==',')) { 
p[i]=0;
IdxTy n=myatoi(p+start); 
if (r) {
for(IdxTy j=prior; j<=n; ++j) selvec.push_back(j); 
MM_ERR(" selecting "<<MMPR2(prior,n))
r=false;
}
else
{ selvec.push_back(n); MM_ERR(" selecting "<<MMPR(n)) }
start=i+1;
prior=n;
} // comma  
if (c=='-' )
{ p[i]=0; IdxTy n=myatoi(p+start); r=true; prior=n; start=i+1; } // dash
}
const IdxTy szs=selvec.size();
if (szs!=0){ cu.ptr(selvec.back()); } 
} // select_range



StrTy view(MyMsg & msg , const ViewParam & vp=ViewParam())
{
Ss ss;

ss<<msg.header_string(vp.header_fields()); //  at.ascii();
ss<<CRLF;
MM_LOOP(ii,vp.header_lines())
{
const StrTy v=msg.hval_lc((*ii));
ss<<(*ii)<<" "<<v<<CRLF;
}
ss<<tree_attachments(msg,vp.tree_fields()); //  at.ascii();

AltRank & alts=vp.alt_prefs();
const IdxTy szalts=alts.size();
switch(vp.mode())
{
case 0: { ss<<msg.text(alts); break; } 
case 1: { ss<<msg.body(); break; } 
case 2: { ss<<msg.text(m_hand,alts); break; } 
case 3: { ss<<msg.verbatim(); break; } 

} ; // mode
const IdxTy sz2=alts.size();
if (sz2!=szalts)
{
ss<<" new alt found "<<CRLF;
MM_LOOP(ii,alts)
{
const StrTy t=(*ii).first;
const IdxTy r=(*ii).second;
ss<<MMPR2(t,r)<<CRLF;
} // ii 

}
return ss.str();
}
 



// could be polymorphic 
class thread_entry_param
{
public:
thread_entry_param()
: p(0),opcode(0),nthreads(0),name(),flags(0)
{}
void load_folder(Myt * _p, const StrTy & _name) { p=_p; name=_name; opcode=1; }

void download_folder(Myt * _p, const StrTy & _name, const StrTy & _name2, const StrTy & _dest, const IdxTy _nthreads)
{ p=_p; name=_name; name2=_name2; dest=_dest; opcode=2; nthreads=_nthreads; }

void download_all_folders(Myt * _p, const StrTy & _name, const StrTy & _name2, const StrTy & _dest, const IdxTy _nthreads)
{ p=_p; name=_name; name2=_name2; dest=_dest; opcode=3; nthreads=_nthreads; }



Myt * p;
IdxTy opcode,nthreads;
StrTy name;
StrTy name2;
IdxTy flags;
StrTy dest;
}; // thread_entry_param

class df_thread_param
{
public:
df_thread_param() :p(0),fd(0),x(0),serial(0),buffer_msg(true) {}

Myt * p;
FolderDownloader* fd;   
LibetpanReader * x; 
StrTy dest;
IdxTy serial;
bool buffer_msg;
}; // df_thread_param

IdxTy  cancel_cancel()
{
enter_serial();
m_exit_bg=!true; // cancel this is user initiated another 
IdxTy signal=m_signal; 
exit_serial();
return signal;
}

bool signal_or_die(bool & signalled , IdxTy & signal)
{
enter_serial();
// to notify user, set a flag do NOT do IO holding the mutex. 
//if (m_exit_bg) break;  
if (m_exit_bg){ exit_serial(); return true; } //   break;   } 
if (m_signal!=signal){ //MM_ERR(MMPR2(signal,p->m_signal))  
signalled=true; signal=m_signal; } 
exit_serial();
return false;
}


void write_to_file( const StrTy & dest,  MyMsg & msg, IdxTy & qc ) // nothrow
{
//MM_ERR(" start write")
//const StrTy v=m_formatter.convert_to_mbox(msg);
StrTy v;
// this is only one thread
IdxTy res=m_formatter.to_mbox(v,msg);

// ned to move this crap into the "safe os"
enter_serial2(); // TODO FIXME regardless of the name it does nothing doh. 
SafeOsTy sos(dest);
sos<<v; // need to be nothrow now 
sos.finish();
exit_serial2();
//MM_ERR(" end  write")
if (!sos.good()) { qc=2; }

} // write_to_file


// it would be easier to just call a member function asap
static void*  df_thread_entry(void * _msgp)
{
typedef df_thread_param Pt;
Myt * p= ((Pt*)_msgp) ->p;
Pt  msg=* (Pt*)_msgp;
StrTy dest="";
StrTy a="";
try { 
p->thread_start();
delete (Pt*) _msgp;
IdxTy signal=p->cancel_cancel();

FolderDownloader& fd=*(msg.fd);   
LibetpanReader & x=*(msg.x); 
x.set_verbosity(p->m_verbosity);
// make the conns in parllel lol
x.post_assignment_fix();
a=x.account(); // 
dest= msg.dest;
const bool save_to_file=(dest.length()!=0);
const bool buffer_msg=msg.buffer_msg;// true;
const bool track_headers=true;
IdxTy i=0; 
IdxTy tries=0,fails=0;
// may be other updating 
// this interface is not good for multi threads, ok for now. 
while  (fd.more_to_download()&&x.worth_trying())
{
++tries;
bool signalled=false;
//const IdxTy m=n-i; // get from the reader. 
if (p->signal_or_die(signalled , signal)) break;

if (!x.worth_trying()) { MM_ERR(MMPR(msg.serial)<< " not worth trying ") }

const IdxTy m= fd.next_message(); 
//MM_ERR(" FCK "<<MMPR(m))
// assert? this was redundant with while predicate 
if (m==BAD) { break; }
x.set_message(m);
x.get_message();
const StrTy & h= x.header();
const StrTy & b= x.body();
const StrTy & v= x.verbatim();
const StrTy & sfs= x.string_flags_string();
//MyMsg nmsg=MyMsg(h,b,v,sfs,a);
//MyMsg nmsg=MyMsg(v,a);
MyMsg nmsg=MyMsg(v,a,0,fd.uuid(m));
nmsg.flag_string(sfs);
//MM_ERR(" fck "<<MMPR2(m,fd.uuid(m)))
// too late, the kids have already been made 
// and it would be dumb to propogate 
//nmsg.server_uuid(fd.uuid(m));
//MM_ERR(" fck "<<MMPR3(m,nmsg.server_uuid(),fd.uuid(m)))
Ss ssd;
ssd<<MMPR3(m,x.account(),fd.message_data(m));
// #error need to include uuid and make unique file for save
nmsg.my_data(ssd.str());
if (v.length()!=0)
{
if (buffer_msg) { 
// the message store is suppposed to be safe and maybe this caused deadlock 
//p->enter_serial2(); //   
//MyMsg msg=MyMsg(h,b,a);
//MM_ERR(" on msg  "<<MMPR3(i,msg.header().length(),b.length()))
p->m_messages.append(nmsg); 
//p->exit_serial2(); //   
} // buffer
if (track_headers)
{
// keep track of the number of times various header keys come up fwiw 

} // track_headers
IdxTy qc=0;
if (save_to_file)
{
//p->write_to_file( qc,  dest,v);  // nothrow
//p->write_to_file( qc,  dest,v,sfs);  // nothrow
p->write_to_file(dest, nmsg, qc);  // nothrow
if ( qc!=0) { MM_ERR(" no good write  "<<MMPR3(m,dest,a)) ++fails; }
}
fd.last_message_quality(qc,m);
} else
{
// doing nothig AFAICT will retry excet the first zero entry is bad 
// somay infintie loop 
fd.last_message_quality(1,m);
++fails; 
}
// these cause deadlock of unkonwn cause everyone waiting for enter_Serial2
//p->enter_serial2(); //   
if (signalled){  MM_ERR(a<< "at msg "<<MMPR4(msg.serial,m,tries,fails)<< MMPR3(i,h.length(),b.length())) ; 
//p->exit_serial2(); //   
// TODO this does not wotk if the lead thread dies first lol 
if (msg.serial==0){  MM_ERR(MMPR(fd.status())) } 

 } // signalled  
++i; 
} // while
// only the last to die needs this. 

    } catch (...) {
//MM_ERR(" throwed something ")
p->exit_serial2(); // not sure why this occurs but this did not fix it. 
p->exit_serial();
p->thread_done();
if (p->verbosity(1)) { 
   // https://stackoverflow.com/questions/315948/c-catching-all-exceptions
  std::exception_ptr pe = std::current_exception();
        MM_ERR((pe ? pe.__cxa_exception_type()->name() : "null")) 
MM_ERR( " thread throws "<<MMPR2(msg.serial,a))
}
        throw;
}
p->thread_done();
if ( p->verbosity(10) ) { MM_ERR( " thread terminated normally "<<MMPR2(msg.serial,a)) } 

return 0; 
} // df_thread_entry
void download_all_folders_mt( const StrTy & name, const StrTy & name2, const IdxTy flags, const thread_entry_param tep)
{
typedef  df_thread_param Tp;
const IdxTy nthreads=tep.nthreads;
//if (nthreads<2){
//  download_folder( x,  name2, flags, tep);
//return; } 
try { 
thread_start();
LibetpanReader & x= m_libet_rdr_map[name]; //

const bool buffer_msg = m_flp.buffer_msg(); // true; 
const IdxTy n=x.count();
const StrTy a=x.account(); //  x.m_user+StrTy(";")+host+StrTy(";")+folder;
const StrTy & dest=tep.dest;
bool save_to_file=(dest.length()!=0);
MM_ERR("download folder   "<<MMPR4(a,n,name2,dest))
FolderDownloader fd; // putting this into the map requires synchronization. 
fd.set_verbosity(m_verbosity);
x.set_index(fd.idx());
Ragged & fnames= m_ragged_map[name2];
//enter_serial();
//IdxTy signal=m_signal;
//exit_serial();
// this should allow flags for unseen etc. 
MM_ERR(" not passing flags to get index as want all now, "<<MMPR(flags))
x.get_index();
fd.download_mailbox(name2); // get index and compare 
//IdxTy i=0; 
// these all stay in scope until join so all must exist.
// a tree of disconnects to children is ok. 
LibetpanReader lpr[nthreads];
Tp *  msgp[nthreads];
for(IdxTy  i=0; i<nthreads; ++i)
{
msgp[i]= new Tp();
Tp & ms=*msgp[i];
ms.p=(this);
ms.fd=&fd;
lpr[i]=x; // this wastes x
// this is serializing all the conns wtf
//lpr[i].post_assignment_fix();
ms.x= &lpr[i];
ms.dest=dest;
ms.serial=i;
ms.buffer_msg=buffer_msg;
} // i 
// there is also lap and set 
	m_cm.mark("download_all_mt");

   auto pfunc= Myt::df_thread_entry;
launch(nthreads, pfunc ,(void**)&msgp[0]);
fd.save_historical_data(name2);
// this should be in "finally" block but just measurig now
	m_cm.cumdel("download_all_mt_total","download_all_mt");
D td=m_cm.get_time("download_all_mt");
D ts=m_cm.get_time("download_all_mt_total");
MM_ERR(MMPR3(ts,td,fd.status()))

    } catch (...) {
thread_done();
        throw;
}
thread_done();
} // download_all_folders_mt


void download_folder_mt( LibetpanReader & x, const StrTy & name2, const IdxTy flags, const thread_entry_param tep)
{
typedef  df_thread_param Tp;
const IdxTy nthreads=tep.nthreads;
if (nthreads<2){
  download_folder( x,  name2, flags, tep);
return; } 
try { 
thread_start();
const bool buffer_msg = m_flp.buffer_msg(); // true; 
const IdxTy n=x.count();
const StrTy a=x.account(); //  x.m_user+StrTy(";")+host+StrTy(";")+folder;
const StrTy & dest=tep.dest;
bool save_to_file=(dest.length()!=0);
MM_ERR("download folder   "<<MMPR4(a,n,name2,dest))
FolderDownloader fd; // putting this into the map requires synchronization. 
fd.set_verbosity(m_verbosity);
x.set_index(fd.idx());
//enter_serial();
//IdxTy signal=m_signal;
//exit_serial();
MM_ERR(" not passing flags to get_index for now want all "<<MMPR(flags))
x.get_index();
fd.download_mailbox(name2); // get index and compare 
//IdxTy i=0; 
// these all stay in scope until join so all must exist.
// a tree of disconnects to children is ok. 
LibetpanReader lpr[nthreads];
Tp *  msgp[nthreads];
for(IdxTy  i=0; i<nthreads; ++i)
{
msgp[i]= new Tp();
Tp & ms=*msgp[i];
ms.p=(this);
ms.fd=&fd;
lpr[i]=x; // this wastes x
// this is serializing all the conns wtf
//lpr[i].post_assignment_fix();
ms.x= &lpr[i];
ms.dest=dest;
ms.serial=i;
ms.buffer_msg=buffer_msg;
} // i 
// there is also lap and set 
	m_cm.mark("download_mt");

   auto pfunc= Myt::df_thread_entry;
launch(nthreads, pfunc ,(void**)&msgp[0]);
fd.save_historical_data(name2);
// this should be in "finally" block but just measurig now
	m_cm.cumdel("download_mt_total","download_mt");
D td=m_cm.get_time("download_mt");
D ts=m_cm.get_time("download_mt_total");
MM_ERR(MMPR3(ts,td,fd.status()))

    } catch (...) {
thread_done();
        throw;
}
thread_done();
} // download_folder_mt


// tep is by value although need to fix the threading. 
void download_folder( LibetpanReader & x, const StrTy & name2, const IdxTy flags, const thread_entry_param tep)
{
try { 
thread_start();
const bool buffer_msg = true; 
const IdxTy n=x.count();
const StrTy a=x.account(); //  x.m_user+StrTy(";")+host+StrTy(";")+folder;
const StrTy & dest=tep.dest;
bool save_to_file=(dest.length()!=0);
MM_ERR("download folder   "<<MMPR4(a,n,name2,dest))
FolderDownloader fd; // putting this into the map requires synchronization. 
fd.set_verbosity(m_verbosity);
x.set_index(fd.idx());
enter_serial();
IdxTy signal=m_signal;
exit_serial();
x.get_index();
fd.download_mailbox(name2); // get index and compare 
IdxTy i=0; 
IdxTy tries=0,fails=0;
// may be other updating 
// this interface is not good for multi threads, ok for now. 
while  (fd.more_to_download())
{
++tries;
bool signalled=false;
//const IdxTy m=n-i; // get from the reader. 
enter_serial();
// to notify user, set a flag do NOT do IO holding the mutex. 
//if (m_exit_bg) break;  
if (m_exit_bg){ exit_serial();  break;   } 
if (m_signal!=signal){  signalled=true; signal=m_signal; } 
exit_serial();
const IdxTy m= fd.next_message(); 
if (m==BAD) {
MM_ERR(" bad next message wtf ")
}
x.set_message(m);
x.get_message();
const StrTy & h= x.header();
const StrTy & b= x.body();
const StrTy & v= x.verbatim();
if (v.length()!=0)
{
IdxTy qc=0;
if (save_to_file)
{
//MM_ERR(" begin save "<<MMPR(dest))
//MM_ERR(" begin save "<<MMPR(v.length()))
std::ofstream ofs(dest,std::ios_base::app);
ofs<<v;
if (!ofs.good()) 
{
MM_ERR(" no good write e "<<MMPR2(dest,a))
qc=2;
}
//MM_ERR(" end save "<<MMPR(dest))
}
fd.last_message_quality(qc,m);

} else
{
// doing nothig AFAICT will retry excet the first zero entry is bad 
// somay infintie loop 
fd.last_message_quality(1,m);
++fails; 
}
//y.view(h,b);
if (signalled) MM_ERR(a<< " downloading on msg "<<MMPR2(tries,fails)<< MMPR3(i,h.length(),b.length()))
if (buffer_msg) { 
MyMsg msg=MyMsg(h,b,a);
//MM_ERR(" on msg  "<<MMPR3(i,msg.header().length(),b.length()))
m_messages.append(msg); } // buffer
++i; 
} // while
// only the last to die needs this. 
fd.save_historical_data(name2);
MM_ERR(" done downloading folder"<<MMPR(a))
    } catch (...) {
thread_done();
        throw;
}

thread_done();


} // download_folder

//typedef std::map<StrTy, LibetpanReader> LibetMap;
//LibetpanReader & x= m_libet_rdr_map[sname]; //
void load_folder( LibetpanReader & x, const IdxTy flags, const thread_entry_param tep)
{
thread_start();
try{ 
//x.all_at_once(user, pass, host, port,  msg,  folder);
//x.set_account_info(rwm,folder);
//x.get_index();
const IdxTy n=x.count();
const StrTy a=x.account(); //  x.m_user+StrTy(";")+host+StrTy(";")+folder;
//typedef  mjm_pawnoff_mail Hand;
//Hand y;
MM_ERR("load folder   "<<MMPR2(a,n))
enter_serial();
IdxTy signal=m_signal;
exit_serial();
for (IdxTy i=0; i<n; ++i)
{
bool signalled=false;
const IdxTy m=n-i;
enter_serial();
// to notify user, set a flag do NOT do IO holding the mutex. 
if (m_exit_bg){ exit_serial();  break;   } 
if (m_signal!=signal){  signalled=true; signal=m_signal; } 
exit_serial();
x.set_message(m);
x.get_message();
const StrTy & h= x.header();
const StrTy & b= x.body();
//y.view(h,b);
if (signalled) MM_ERR(a<< " loading  on msg  "<<MMPR3(i,h.length(),b.length()))
MyMsg msg=MyMsg(h,b,a);
//this is the hokey index thing again... 
// this is even worse in multhreaded environ need lock arond
// thjis atomic operation.. or nuclear lol 
//const IdxTy idx=m_messages.size();
//m_messages[idx]=msg;
//MM_ERR(" on msg  "<<MMPR3(i,msg.header().length(),b.length()))
m_messages.append(msg);
}
MM_ERR(" done loading folder"<<MMPR(a))
    } catch (...) {
        //throw;
}

thread_done();
} // load_folder
void thread_start() { enter_serial(); ++m_thread_total;  ++m_thread_count; exit_serial(); } 
void thread_done() { enter_serial(); --m_thread_count; exit_serial(); } 



////////////////////////////////////////////////////
//static void mt_shrink_group(mt_shrink_group_param * msgp)
static void*  thread_entry_func(void * _msgp)
{
thread_entry_param * msgp= (thread_entry_param*)_msgp;
Myt * p= msgp->p;
//MM_ERR(" in thread func  ")
p->enter_serial();
p->m_exit_bg=!true; // cancel this is user initiated another 
p->exit_serial();
thread_entry_param newtep=*msgp;
const IdxTy opcode=msgp->opcode;
const StrTy  name= msgp->name; // lol refs no good doh
const StrTy  name2= msgp->name2; // lol refs no good doh
const IdxTy flags=msgp->flags;
delete msgp;
if (p->verbosity(20)) { MM_ERR(" in  tep "<<MMPR3(name,opcode,flags)) } 
if (opcode==1)
{
//MM_ERR(" calling tep "<<MMPR(name))
LibetpanReader & x= p->m_libet_rdr_map[name]; //
p->load_folder( x, flags, newtep);
return 0; 
}
if (opcode==2) // download_folder
{
//MM_ERR(" calling tep "<<MMPR(name))
LibetpanReader & x= p->m_libet_rdr_map[name]; //
//p->download_folder( x, name2,flags,newtep);
p->download_folder_mt( x, name2,flags,newtep);
return 0; 
}
if (opcode==3) // download_all_folders
{
//MM_ERR(" calling tep "<<MMPR(name))
//LibetpanReader & x= p->m_libet_rdr_map[name]; //
//p->download_folder( x, name2,flags,newtep);
p->download_all_folders_mt( name, name2,flags,newtep);
return 0; 
}



return 0; 
} // thread_entry_func

// TODO FIXME this double pointer needs to be fixed in earlier code 
//void launch(const IdxTy nthreads, void *(*thread_function)(void* ) ,thread_entry_param *msgp)
void launch(const IdxTy nthreads, void *(*thread_function)(void* ) ,void **msgp)
{
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, thread_function, msgp[i] ); }
   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}
void fire_and_forget(thread_entry_param *msgp)
{
//MM_ERR(" firing")
   auto pfunc= Myt::thread_entry_func;
   IdxTy nthreads=1;
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, pfunc, msgp+i ); }
//   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}

#ifdef HAVE_MIKE_MUTT_SERVER
void cmd_mutt_server(Cip & cip , LocalVar & lv )
{
const StrTy config=cip.p1;
//const StrTy folder=cip.p2;
//const StrTy progfile=cip.wif(3); // can store in rags map here ... 
//const StrTy dest=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
//const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_mutt_server "<<CRLF;
//ss<<" Example: -cmd \"download-folder foo INBOX oldcrap mboxx 25\""<<CRLF;
//ss<<" ragin folder progfile dest nthreads "<<CRLF;
//ss<<" ragin : the session or account rag name  "<<CRLF;
//ss<<" folder : folder known by server  "<<CRLF;
//ss<<" progfil : progress file, journal not implemnted no need yet   "<<CRLF;
//ss<<"  dest : location to store downloaded mail in mbox format   "<<CRLF;
//ss<<"  nthreads : connetion count, extras dies without consequence   "<<CRLF;
cip.help(ss.str());
return; 
}
mutt_main_server(config.c_str(),0);

} // cmd_mutt_server

void cmd_mutt_load(Cip & cip , LocalVar & lv )
{
const StrTy config=cip.p1;
const IdxTy flags=myatoi(cip.p2);
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_mutt_load "<<CRLF;

cip.help(ss.str());
return; 
}
MyMsgMap & m=m_messages;
MuttInterface mi;
mutt_set_init_only(true);
int i=0;
int rc=mutt_main_server(config.c_str(),0);
if (rc!=0)
{
MM_ERR(" failed to open mutt "<<MMPR(rc))
return;
}
struct Email ** elist= mutt_email_list(flags);
if (elist!=0)
{
while (elist[i]!=0)
{

MM_ERR("adding "<< MMPR(i))
MyMsg msg;
mi.convert(msg,elist[i]);

MM_ERR("adding2  "<< MMPR(i))
m[elist[i]->index]=msg;
MM_ERR("adding3  "<< MMPR(i))
++i;
} // i 

}
MM_ERR("free "<<MMPR(i))
mutt_mjm_free(elist);
shutdown_mutt();
MM_ERR(MMPR(i))
} // cmd_mutt_load


void cmd_mutt_proc(Cip & cip , LocalVar & lv )
{
// this is a neomutt style file pass to neomutt library
const StrTy config=cip.p1;
const IdxTy flags=myatoi(cip.p2);
const StrTy fn=cip.wif(3);
const IdxTy mflags=myatoi(cip.wif(4)); // 0;

// TODO should have been fixed but code not tested yet 
// xTODO xFIXME this will hang waiting on a prompt if the
// server is valud but not the fold as it wants to 
// ask if you want to make it... 
// mutt-proc proc_yahoo 9 mproc 0 imaps://imap.mail.yahoo.com:993/old
const StrTy AFCK=(cip.wif(5)); // 0;
if (cip.help())
{
Ss ss;
ss<<cip.cmd()<<" cmd_mutt_proc "<<CRLF;

cip.help(ss.str());
return; 
}


const bool unread_only=Bit(flags,0);
const bool reset_on_fail=Bit(flags,1);
const bool server_update=!Bit(flags,2);
const bool move_to_archive=Bit(flags,3);
const bool skip_setting_read=Bit(flags,4);
MM_ERR(MMPR4(config,flags,mflags,fn))
MM_ERR(MMPR(skip_setting_read))
MM_ERR(MMPR4(unread_only,reset_on_fail,server_update,move_to_archive))

MyMsgMap & m=m_messages;


struct Email ** elist= open_get_mutt( m,  config, unread_only);

MM_ERR(" returned messages "<<MMPR2(m.size(),(elist!=0)))
if (elist!=0) {
//MyMsgMap & m=m_messages; // def?m_messages:m_messages_map[mbox];
MailProc & mp= m_proc_map[fn];
mp.respond_clear();
mp.respond_imap();
//proc_msg_map(mp,mmm);
MM_SZ_LOOP(k,m,szm)
{
const int idx=elist[k]->index;
const MyMsg & msg=m[idx];
//MM_ERR(MMPR2(i,msg.verbatim()))
//MM_ERR(MMPR4(szm,k,idx,msg.verbatim().size()))
MM_ERR(MMPR4(szm,k,idx,msg.from_address(0).size()))
// processor or handler should keep a journal and log
// but can configure here too. 
proc_io_type pio;
const IdxTy prc=mp.proc(pio,msg,mflags);
MM_ERR(MMPR(prc))
if (server_update) { 
//MM_ERR("ignore fcking read flag  try  server update ")
if (!skip_setting_read){ 
if (prc==0) set_read_flag(elist[k]);
else if (reset_on_fail) reset_read_flag(elist[k]);
} // skip_setting_read
MM_ERR(" server update arch "<<MMPR2(move_to_archive, flags) )
if (move_to_archive)
{
// setting flag to 3 should allow creation of folders
// with proper domain prefix ... 
const int mvrc= mjm_mutt_move_msg(elist[k], AFCK.c_str(), 1 );
MM_ERR(" server update done"<<MMPR(mvrc) )
} 

MM_ERR(" server update done")
}
/*
char * subj="subect";
char * text=" message text";
char ** lto= new char*[2];
lto[0]="marchywka@hotmail.com";
lto[1]=0;
// the reply function fails but sending new messag ok 
//int rcs=mjm_send_via_mutt(lto,subj, text,0,elist[k]);
int rcs=mjm_send_via_mutt(lto,subj, text,0,0);

MM_ERR(MMPR(rcs))
*/
} // k loop 

} // elist !=0 
//MM_ERR("free "<<MMPR(i))
if (elist) mutt_mjm_free(elist);
shutdown_mutt();
//MM_ERR(MMPR(i))
} // cmd_mutt_proc

struct Email ** open_get_mutt(MyMsgMap & m,  const StrTy & mutt_config, const bool unread_only)
{
m.clear();
MuttInterface mi;
mutt_set_init_only(true);
int i=0;
//int rc=mutt_main_server(config.c_str(),0);
int rc=mutt_main_server(mutt_config.c_str(),0);
if (rc!=0)
{
MM_ERR(" failed to open mutt "<<MMPR(rc))
return 0 ;
}
struct Email ** elist= mutt_email_list(unread_only?1:0);
if (elist!=0)
{
std::map<IdxTy,IdxTy> serial;
while (elist[i]!=0)
{
MyMsg msg;
// StrTy account(){ return  m_user+StrTy(";")+m_host+StrTy(";")+m_folder; }

mi.convert(msg,elist[i]);
serial[elist[i]->index]=i;
m[elist[i]->index]=msg;
++i;
} // i 
} // elist 

return elist; 

} // open_get_mutt


void cmd_mutt_fck(Cip & cip , LocalVar & lv )
{
const StrTy config=cip.p1;
const IdxTy flags=myatoi(cip.p2);
const bool unread_only=Bit(flags,0);
const bool reset_on_fail=Bit(flags,1);
const bool server_update=!Bit(flags,2);
const StrTy fn=cip.wif(3);
MM_ERR(MMPR3(config,flags,fn))
if (cip.help())
{

}
MuttInterface mi;
mutt_set_init_only(true);
int i=0;
int rc=mutt_main_server(config.c_str(),0);
if (rc!=0)
{
MM_ERR(" failed to open mutt "<<MMPR(rc))
return;
}
// need to fix mjm_send signature 
//const 
const char * subj="subect";
//const 
const char * text=" message text";
//const 
const char * * lto= new const char*[2];
lto[0]="marchywka@hotmail.com";
lto[1]=0;
const char *  * aaa= new const char*[3];
aaa[0]="xxx.txt";
aaa[1]=0;
// the reply function fails but sending new messag ok 
//int rcs=mjm_send_via_mutt(lto,subj, text,0,elist[k]);
// this signature needs to be fixed... 
int rcs=mjm_send_via_mutt(lto,subj, text,aaa,0);

MM_ERR(MMPR(rcs))

shutdown_mutt();

} // cmd_mutt_fck 





#endif

void cmd_download_folder(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
const StrTy folder=cip.p2;
const StrTy progfile=cip.wif(3); // can store in rags map here ... 
const StrTy dest=cip.wif(4); // can store in rags map here ... 
const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(nthreads==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_download_folder "<<CRLF;
ss<<" Example: -cmd \"download-folder foo INBOX oldcrap mboxx 25\""<<CRLF;
ss<<" ragin folder progfile dest nthreads "<<CRLF;
ss<<" ragin : the session or account rag name  "<<CRLF;
ss<<" folder : folder known by server  "<<CRLF;
ss<<" progfil : progress file, journal not implemnted no need yet   "<<CRLF;
ss<<"  dest : location to store downloaded mail in mbox format   "<<CRLF;
ss<<"  nthreads : connetion count, extras dies without consequence   "<<CRLF;
cip.help(ss.str());
return; 
}

StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
   const StrTy user=(rwm.get_string("user","missing"));
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
LibetpanReader & x= m_libet_rdr_map[sname]; //
x.set_verbosity(m_verbosity);
MM_ERR(" loading in background "<<MMPR2(a,dest))
x.set_account_info(rwm,folder);
MM_ERR(" loading in background "<<MMPR2(a,dest))
// must be defered until we get a folder downloader or that should be created here before lauching linked 
// getters 
//x.get_index();
//const IdxTy n=x.count();
//MM_ERR(" downloading in background "<<MMPR(n))
//thread_entry_param tep;
//tep.load_folder(this,sname);
// except that tep goes out of scope immediately... 
//fire_and_forget(&tep);
thread_entry_param*  tepp= new thread_entry_param();
tepp->download_folder(this,sname,progfile,dest,nthreads);
// could wait on lock but then can hang. If called deletes that is ok too
// and if not something else is wrong an leak no big deal 
// keep it around for debugger 
fire_and_forget(tepp);

} // cmd_download_folder


void cmd_test_imap(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
const StrTy folder=cip.p2;
const StrTy uids=cip.wif(3); // can store in rags map here ... 
const StrTy sflags=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
const IdxTy uididx=myatoi(uids);
const IdxTy rflags=myatoi(sflags);

if (cip.help()||(ragin.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_test_imap "<<CRLF;
ss<<cip.cmd()<<" test-imap ragin folder uids sflags "<<CRLF;
ss<<" list and modify flags for an imap mailbox folder.  "<<CRLF;
ss<<" ragin : the session or account rag name  "<<CRLF;
ss<<" folder : folder or mbox name known to server   "<<CRLF;
ss<<" progfil : uid or index of message to modify    "<<CRLF;
ss<<"  sflags :  0 , set_seen; 1, not_uid   "<<CRLF;
cip.help(ss.str());
return; 
}
StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
   const StrTy user=(rwm.get_string("user","missing"));
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
LibetpanReader & x= m_libet_rdr_map[sname]; //
x.set_verbosity(m_verbosity);
//MM_ERR(" ... "<<MMPR2(a,dest))
MM_ERR(" ... "<<MMPR4(ragin,folder,uididx,rflags))
x.set_account_info(rwm,folder);

//IdxTy fast_mark2( IdxTy uuid=0, const IdxTy flags=0)
//const bool  set_seen=Bitt(flags,0);
//const bool  not_uid=Bitt(flags,1);

x.fast_mark2(uididx, rflags);


}

void cmd_download_all_folders(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
const StrTy folder=cip.p2;
const StrTy progfile=cip.wif(3); // can store in rags map here ... 
const StrTy dest=cip.wif(4); // can store in rags map here ... 
const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(nthreads==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_download_folder "<<CRLF;
ss<<" Example: -cmd \"download-folder foo INBOX oldcrap mboxx 25\""<<CRLF;
ss<<" ragin folder progfile dest nthreads "<<CRLF;
ss<<" ragin : the session or account rag name  "<<CRLF;
ss<<" folders : ragged holding folder names  "<<CRLF;
ss<<" progfil : progress root file, journal not implemnted no need yet   "<<CRLF;
ss<<"  dest : root location to store downloaded mail in mbox format   "<<CRLF;
ss<<"  nthreads : connetion count, extras dies without consequence   "<<CRLF;
cip.help(ss.str());
return; 
}

StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
   const StrTy user=(rwm.get_string("user","missing"));
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
LibetpanReader & x= m_libet_rdr_map[sname]; //
x.set_verbosity(m_verbosity);
MM_ERR(" loading in background "<<MMPR2(a,dest))
x.set_account_info(rwm,folder);
MM_ERR(" loading in background "<<MMPR2(a,dest))


thread_entry_param*  tepp= new thread_entry_param();
tepp->download_all_folders(this,sname,progfile,dest,nthreads);
// could wait on lock but then can hang. If called deletes that is ok too
// and if not something else is wrong an leak no big deal 
// keep it around for debugger 
fire_and_forget(tepp);

} // cmd_download_all_folders







void cmd_start_session(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
const StrTy folder=cip.p2;
StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
   const StrTy user=(rwm.get_string("user","missing"));
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
LibetpanReader & x= m_libet_rdr_map[sname]; //
x.set_verbosity(m_verbosity);
x.set_account_info(rwm,folder);
x.get_index(1);
const IdxTy n=x.count();
MM_ERR(" session connected and found index info   "<<MMPR2(n,a))

} // cmd_start_session


void cmd_list_folders(Cip & cip , LocalVar & lv )
{
const StrTy sname=cip.p1;
// this must by blank to work when session is inbox or lower 
const StrTy root=""; // cip.p2;
//const StrTy root= cip.p2;
MM_ERR(" ADCASDCSADCAS may not work with non-blank mbox "<<MMPR(root))
const StrTy ragout=cip.p2; // cip.wif(3);

if (cip.help()||(sname.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_list_folders sname [ ragout ] "<<CRLF;
ss<<" list folders in whatever box session sname was started  "<<CRLF;
ss<<" ragout : optional name of ragged file to store mbox list for processing.   "<<CRLF;
cip.help(ss.str());
return; 
}
//const bool def=(mbox=="-");



LibetpanReader & x= m_libet_rdr_map[sname]; //
MM_ERR(" getting folder info for "<<x.account()<<MMPR(root) )

typedef LibetpanReader::FolderList Fl ;

Fl * pfl=x.get_folder_status(root);
if ( pfl==0) return; 
Fl & pflr=*pfl;
MM_ERR(" returned ok   "<< MMPR(pflr.size()))
MM_LOOP(ii,pflr) { MM_ERR((*ii).dump()); }
if (ragout.length()!=0)
{
Ragged & x = m_ragged_map[ragout];
x.append_list(pflr);

}
delete pfl;

} // cmd_list_folders






void cmd_load_folder(Cip & cip , LocalVar & lv )
{
const StrTy ragin=cip.p1;
const StrTy folder=cip.p2;
const StrTy sflags=cip.wif(3);
// bit 0 : do not fecth list, non sense here
// bit 1 : get envelopes too
MM_ERR(MMPR3(ragin,folder,sflags));
 
const IdxTy flags=myatoi(sflags);
StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
   const StrTy user=(rwm.get_string("user","missing"));
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
LibetpanReader & x= m_libet_rdr_map[sname]; //
x.set_verbosity(m_verbosity);
x.set_account_info(rwm,folder);
x.get_index(flags);
const IdxTy n=x.count();
MM_ERR(" loading in background "<<MMPR2(n,a))
//thread_entry_param tep;
//tep.load_folder(this,sname);
// except that tep goes out of scope immediately... 
//fire_and_forget(&tep);
thread_entry_param*  tepp= new thread_entry_param();
tepp->load_folder(this,sname);
// could wait on lock but then can hang. If called deletes that is ok too
// and if not something else is wrong an leak no big deal 
// keep it around for debugger 
fire_and_forget(tepp);

} // cmd_load_folder


void cmd_stop_loads(Cip & cip , LocalVar & lv )
{
enter_serial();
m_exit_bg=true;
exit_serial();

}
void cmd_signal(Cip & cip , LocalVar & lv )
{
// in theory a thread could miss exactly 1<<64 of these 
MM_MSG(" message start ")

enter_serial();
MM_MSG(" message start entered  ")
++m_signal; 
IdxTy signal=m_signal;
IdxTy threads=m_thread_count;
IdxTy total_threads=m_thread_total;
exit_serial();
MM_ERR(" signal count now "<<MMPR3(signal,threads,total_threads))
}

MyMsgMap & mmm_lut(const StrTy & mbox) //  const
{
const bool def=(mbox=="-");
if (def){
MM_ERR(" wiating for default box to be done") 
 wait_to_finish(); 
}
return def?m_messages:m_messages_map[mbox];
}
void cmd_cpy_mbox(Cip & cip , LocalVar & lv )
{

const StrTy cmd=cip.cmd();
const StrTy d=cip.p1;
const StrTy s=cip.p2;
if (cip.help()||(d.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_cpy_mbox "<<CRLF;
ss<<" mbox filename-in "<<CRLF;
cip.help(ss.str());
return; 
}
MyMsgMap & md=mmm_lut(d);
MyMsgMap & ms=mmm_lut(s);
md=ms;
if (cmd=="move-mbox") ms.clear();
MM_ERR(MMPR3(cmd,d,s))
} // cmd_cpy_mbox
void cmd_load_mbox(Cip & cip , LocalVar & lv )
{

const StrTy mbox=cip.p1;
const StrTy fn=cip.p2;
//const StrTy progfile=cip.wif(3); // can store in rags map here ... 
//const StrTy dest=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
//const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_load_mbox "<<CRLF;
ss<<" mbox filename-in "<<CRLF;
cip.help(ss.str());
return; 
}
//const bool def=(mbox=="-");
//list_thing(m_messages_map,"m_messages_map");
//MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
MyMsgMap & mmm=mmm_lut(mbox);
load_mbox(mmm,fn);

} // cmd_load_mbox




void cmd_save_att(Cip & cip , LocalVar & lv )
{

const StrTy mbox=cip.p1;
const StrTy fn=cip.p2;
//const StrTy progfile=cip.wif(3); // can store in rags map here ... 
//const StrTy dest=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
//const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_save_att "<<CRLF;
ss<<" mbox filename-in "<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=(mbox=="-");
//list_thing(m_messages_map,"m_messages_map");
MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
//load_mbox(mmm,fn);
save_attachments(mmm);
//void tree_attachments(MyMsgMap & mmm)
} // cmd_save_att

void cmd_tree_att(Cip & cip , LocalVar & lv )
{

const StrTy mbox=cip.p1;
const StrTy fn=cip.p2;
//const StrTy progfile=cip.wif(3); // can store in rags map here ... 
//const StrTy dest=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
//const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_tree_att "<<CRLF;
ss<<" mbox filename-in "<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=(mbox=="-");
//list_thing(m_messages_map,"m_messages_map");
MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
//load_mbox(mmm,fn);
tree_attachments(mmm);
//void tree_attachments(MyMsgMap & mmm)
} // cmd_save_att


IdxTy proc_msg_map(MailProc & mp, MyMsgMap & mmm)
{
IdxTy mflags=0;
MM_SZ_LOOP(i,mmm,szm)
{
const MyMsg & msg=mmm[i];
//MM_ERR(MMPR2(i,msg.verbatim()))
MM_ERR(MMPR2(i,msg.verbatim().size()))
mp.proc(msg,mflags);
}
return 0;
} // proc_msg_map

void cmd_mmp_proc(Cip & cip , LocalVar & lv )
{
//const StrTy mbox=cip.p1;
const StrTy fn=cip.p1;
const StrTy mbox=cip.p2;
const StrTy smode=cip.wif(3);
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_mmp_proc "<<CRLF;
ss<<" mmp mbox mode "<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=(mbox=="-");
MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
MailProc & mp= m_proc_map[fn];
mp.respond_clear();
mp.respond_mbox();
proc_msg_map(mp,mmm);
} // cmd_mmp_proc


void cmd_mmp_config(Cip & cip , LocalVar & lv )
{
//const StrTy mbox=cip.p1;
const StrTy fn=cip.p1;
const StrTy rname=cip.p2;
const StrTy smode=cip.wif(3);
const IdxTy flags=myatoi(smode);
const bool just_dump=Bit(flags,0);
MM_ERR("cmd_mmp_config "<<MMPR4(fn,rname,smode,flags)<<MMPR(just_dump))
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_mmp_config "<<CRLF;
ss<<" mmp mbox mode "<<CRLF;
cip.help(ss.str());
return; 
}
Ragged & r= m_ragged_map[rname];
MailProc & mp= m_proc_map[fn];
//mp.config_rules(r);
if (!just_dump) mp.config(r);
if (just_dump) MM_ERR(mp.dump());
MM_ERR( MMPR3(rname,fn,r.size()))

} // cmd_mmp_config


void cmd_view_src(Cip & cip , LocalVar & lv )
{
const StrTy mbox=cip.p1;
const StrTy fn=cip.p2;
const StrTy smode=cip.wif(3);
const StrTy hdr=cip.wif(4);
const IdxTy flags=myatoi(smode);
const bool just_dump=Bit(flags,0);
const bool show_hdr=Bit(flags,1);

MM_ERR(MMPR2(just_dump,show_hdr))
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_view_src "<<CRLF;
ss<<" mbox msgno smode [ header ]  "<<CRLF;
ss<<" smode 1<<0 just_dump  "<<CRLF;
ss<<" smode 1<<1 show_hdr  "<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=(mbox=="-");
MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
const IdxTy mno=myatoi(fn);
const MyMsg & msg=mmm[mno];

if (just_dump) {  MM_MSG(msg.verbatim()); } 
else m_hand.view(msg.verbatim());
if ( show_hdr)
{
const StrTy hval= (msg.hval_lc(hdr));
MM_MSG(MMPR2(hdr,hval))
}


} // cmd_view_src

void cmd_view(Cip & cip , LocalVar & lv )
{

const StrTy mbox=cip.p1;
const StrTy fn=cip.p2;
const StrTy smode=cip.wif(3);
//const StrTy progfile=cip.wif(3); // can store in rags map here ... 
//const StrTy dest=cip.wif(4); // can store in rags map here ... 
//const StrTy nthreadss=cip.wif(5); // can store in rags map here ... 
//const IdxTy nthreads=myatoi(nthreadss);
// can be after parse 
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_tree_att "<<CRLF;
ss<<" mbox msgno "<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=(mbox=="-");
//list_thing(m_messages_map,"m_messages_map");
MyMsgMap & mmm=def?m_messages:m_messages_map[mbox];
const IdxTy mno=myatoi(fn);
const bool save_att=(fn.length()>0)?(fn.c_str()[fn.length()-1]=='s'):false;
MyMsg & msg=mmm[mno];
//load_mbox(mmm,fn);
//tree_attachments(mmm);
ViewParam vp;
if (smode.length()!=0) vp.mode(myatoi(smode));
if (save_att) m_formatter.decode_parts(msg);
StrTy x= view(msg,vp);
m_hand.view(x);
// save attachments 
//MyMsg & msg=mmm[i];
//void tree_attachments(MyMsgMap & mmm)
} // cmd_view


void cmd_v(Cip & cip , LocalVar & lv )
{

//const StrTy mbox=cip.p1;
const StrTy fn=cip.p1;
const StrTy smode=cip.p2;
const StrTy & flagss=cip.wif(3);
const StrTy & fw= m_flp.mail_clip_file();
const StrTy & box="-";
if (cip.help()||(fn.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_v "<<CRLF;
ss<<" v range mode flags "<<CRLF;
ss<<" flags: 0, view; 1: save to mail_clip_file"<<CRLF;
cip.help(ss.str());
return; 
}
const bool def=true; // (mbox=="-");
const IdxTy flags= myatoi(flagss);

bool view_interactive=!Bit(flags,0);
bool write_to_file=Bit(flags,1);
MM_ERR(MMPR4(fn,flagss,view_interactive,write_to_file))
auto & cursor=m_cursors[box];
std::vector<IdxTy> selvec;
//select_range(selvec,fn,m_cursors["-"]);
select_range(selvec,fn,cursor);
//list_thing(m_messages_map,"m_messages_map");
MyMsgMap & mmm=m_messages; // :m_messages_map[mbox];
Ss ss;
const IdxTy msz=mmm.size();
if (msz==0)
{
MM_ERR(" box  is empty "<<MMPR(box))
return; 
}
MM_LOOP(ii,selvec)
//for (IdxTy i=0; i<idx; ++i)
{
const IdxTy i=(*ii);
if (i>=msz) {cursor.ptr(msz-1);   continue; } 
//const IdxTy mno=myatoi(fn);
const bool save_att=false; // (fn.length()>0)?(fn.c_str()[fn.length()-1]=='s'):false;
MyMsg & msg=mmm[i];
ViewParam vp;
if (smode.length()!=0) vp.mode(myatoi(smode));
if (save_att) m_formatter.decode_parts(msg);
ss<<CRLF<<" ***************** MESSAGE "<<i<<CRLF;
ss<<view(msg,vp);
ss<<CRLF;
} // i 
const StrTy&  x=ss.str(); //  view(msg,vp);

if (view_interactive ) m_hand.view(x);
if (write_to_file)
{
//SafeOsTy os(fw,ios_base::app);
SafeOsTy os(fw);
os<<x<<CRLF;
}
// save attachments 
//MyMsg & msg=mmm[i];
//void tree_attachments(MyMsgMap & mmm)
} // cmd_v




void cmd_escape(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
//const StrTy fn=cip.p2;
//const StrTy smode=cip.wif(3);
if (cip.help()||(cmd.length()==0))
{
Ss ss;
ss<<cip.cmd()<<" cmd_escape "<<CRLF;
ss<<" execute shell escape "<<CRLF;
cip.help(ss.str());
return; 
}
MM_ERR(" executing "<<MMPR(cmd))
//m_hand.cmd(cmd);
m_hand.cmd(StrTy(cmd.c_str()+1));
} // cmd_escape


void cmd_headers(Cip & cip , LocalVar & lv )
{
//#define MM_STR_LOC StrTy( __FILE__)+StrTy(QUOTE(__LINE__))
//#define MM_STR_TIMESTAMP StrTy( __DATE__)+StrTy(QUOTE(__TIME__))

std::vector<StrTy> w;
//w.push_back("Date:");
w.push_back("From:");
w.push_back("Content-Length:");
w.push_back("Subject:");
/*
Date: Thu, 11 Apr 2019 19:27:32 +0000 (UTC)
From: Yahoo <no-reply@cc.yahoo-inc.com>
Reply-To: Yahoo <no-reply@cc.yahoo-inc.com>
To: marchywka@yahoo.com
Message-ID: <97835104.4652.1555010852151@be115.member.bf1.yahoo.com>
Subject: Security settings changed on your Yahoo account
MIME-Version: 1.0
Content-Type: text/html; charset="utf-8"
Content-Transfer-Encoding: 7bit
Content-Length: 7273
*/


Ss ss;
const IdxTy idx=m_messages.size();
for (IdxTy i=0; i<idx; ++i)
{
// the crlf needs to be handled elsewere to get the efficiency goals wtf
const MyMsg & msg=m_messages[i];
// now hval_lc and drop nocrlf
//const StrTy d= nocrlf(msg.hval("Date:"));
//const StrTy f= nocrlf(msg.hval("Frome:"));
//const StrTy su= nocrlf(msg.hval("Subject:"));
//const StrTy len= nocrlf(msg.hval("Content-Length:"));


const StrTy d= (msg.hval_lc("Date:"));
const StrTy f= (msg.hval_lc("From:"));
const StrTy su= (msg.hval_lc("Subject:"));
const StrTy len= (msg.hval_lc("Content-Length:"));


//Hand hand;
Hand&  hand=m_hand;
// for speed want to do all at once which is a design goal but not notable on new machines 
const StrTy ds= nocrlf(hand.dates(d));
// ds has a crlf at end 
//MM_ERR(i<<" "<<ds<< msg.header_string(w))
// MM_ERR(i<<" "<<f<<" "<<ds<<" "<<len<<" "<<su)
ss<<i<<" "<<f<<" "<<ds<<" "<<len<<" "<<su<<CRLF;
}
//Hand x;
Hand&  x=m_hand;
x.view(ss.str());
}
StrTy nocrlf( const StrTy & v) const
{
const IdxTy sss=v.length();
char vx[v.length()+2];
IdxTy ptr=0;
for(IdxTy j=0; j<sss; ++j)
{
char c=v.c_str()[j];
if (c!='\r') if (c!='\n') { vx[ptr]=c; ++ptr;}
}
vx[ptr]=0;
return StrTy(vx); 

}

//void all_at_once(const StrTy & user, const StrTy & pass, const StrTy & host, const IdxTy port, IdxTy n, const StrTy & folder)
void cmd_read_folder(Cip & cip , LocalVar & lv )
//static void testMessageBuilder1(String * path)
{
const StrTy ragin=cip.p1;
//const IdxTy msg=myatoi(cip.p2); // wif(3);
const StrTy folder=cip.p2;
const StrTy sflags=cip.wif(3);
// bit 0 : do not fecth list, non sense here
// bit 1 : get envelopes too

const IdxTy flags=myatoi(sflags);

//const StrTy cmd2=cip.wif(4);
StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
//const Ragged & params=m_ragged_map[cmd2];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
//    const StrTy sport=rwm.get_string("port","missing");
//    const IdxTy port= atoi(sport.c_str());
   const StrTy user=(rwm.get_string("user","missing"));
//   const StrTy pass=rwm.get_string("password","missing");
//MM_ERR(MMPR4(user,pass,host,port)<<MMPR2(msg,folder))
MM_ERR(MMPR3(user,host,folder))
const StrTy a= user+StrTy(";")+host+StrTy(";")+folder;
//typedef std::map<StrTy, LibetpanReader> LibetMap;
LibetpanReader & x= m_libet_rdr_map[sname]; //
//x.all_at_once(user, pass, host, port,  msg,  folder);
x.set_verbosity(m_verbosity);
x.set_account_info(rwm,folder);
x.get_index(flags);
const IdxTy n=x.count();
//Hand y;
Hand&  y=m_hand;
for (IdxTy i=0; i<n; ++i)
{
const IdxTy m=n-i;
x.set_message(m);
x.get_message();
const StrTy & h= x.header();
const StrTy & b= x.body();
IdxTy rc=y.view(h,b);
MM_ERR(" on msg  "<<MMPR2(i,rc))
if ( rc!=0 ) break ; 
MM_ERR(" outdated MyMsg ctor")
MyMsg msg=MyMsg(h,b,a);
//this is the hokey index thing again... 
const IdxTy idx=m_messages.size();
m_messages[idx]=msg;

}


//typedef  mjm_pawnoff_mail Hand;
//Hand y;
//y.view(h,b);


}

void cmd_browse(Cip & cip , LocalVar & lv )
{

std::vector<StrTy> w;
w.push_back("Date:");
w.push_back("From:");
w.push_back("Content-Length:");

IdxTy n=m_messages.size();
//typedef  mjm_pawnoff_mail Hand;
//Hand y;
Hand&  y=m_hand;
for (IdxTy i=0; i<n; ++i)
{
// a ref is supposed to be safe lol 
MyMsg m=m_messages[i];
const StrTy s=m.header_string(w);
MM_ERR(s)
//char c=(char)0; 
//std::cin>> c;
//if ( c=='q' ) break;
y.view(m.header(),m.body());
y.view(m.body(),m.body());
}

}

//void all_at_once(const StrTy & user, const StrTy & pass, const StrTy & host, const IdxTy port, IdxTy n, const StrTy & folder)
void cmd_read_lib_msg(Cip & cip , LocalVar & lv )
//static void testMessageBuilder1(String * path)
{
const StrTy ragin=cip.p1;
const IdxTy msg=myatoi(cip.p2); // wif(3);
const StrTy folder=cip.wif(3);
const StrTy sflags=cip.wif(4);
const IdxTy flags=myatoi(sflags);
StrTy sname=ragin;
const Ragged & server=m_servers_map[ragin];
//const Ragged & params=m_ragged_map[cmd2];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy host=rwm.get_string("hostname","missing"); //
    const StrTy sport=rwm.get_string("port","missing");
    const IdxTy port= atoi(sport.c_str());
   const StrTy user=(rwm.get_string("user","missing"));
   const StrTy pass=rwm.get_string("password","missing");
MM_ERR(MMPR4(user,pass,host,port)<<MMPR2(msg,folder))
//auto & x= m_libet_rdr_map[sname].all_at_once(user, pass, host, port,  msg,  folder);
auto & x= m_libet_rdr_map[sname]; // .all_at_once(user, pass, host, port,  msg,  folder);
x.all_at_once(user, pass, host, port,  msg,  folder,flags);

const StrTy & h= x.header();
const StrTy & b= x.body();
//typedef  mjm_pawnoff<Tr>  Hand;
Hand&  y=m_hand;
y.view(h,b);
//MM_ERR(h);
//MM_ERR(b);

}





void cmd_read_msg(Cip & cip , LocalVar & lv )
//static void testMessageBuilder1(String * path)
{
const StrTy ragin=cip.p1;
const IdxTy msg=myatoi(cip.p2); // wif(3);
const StrTy folder=cip.wif(3);
//const StrTy cmd2=cip.wif(4);
//const Ragged & server=m_servers_map[ragin];
//const Ragged & params=m_ragged_map[cmd2];
//ReadWriteMap rwm;
//server.to_map( rwm);
StrTy sname=ragin;
//SessTy sess;
// this makes a copy. 
//m_sessions_map[sname]=SessTy();
#if HAVE_MAIL_CORE
m_sessions_map[sname].get_message( folder, msg);
#endif

}


//void to_map(ReadWriteMap & rwm) const
void cmd_setup_session(Cip & cip , LocalVar & lv )
//static void testMessageBuilder1(String * path)
{
const StrTy ragin=cip.p1;
//const IdxTy flags=myatoi(cip.p2); // wif(3);
//const StrTy cmd1=cip.wif(3);
//const StrTy cmd2=cip.wif(4);
const Ragged & server=m_servers_map[ragin];
//const Ragged & params=m_ragged_map[cmd2];
ReadWriteMap rwm;
server.to_map( rwm);
StrTy sname=ragin;
//SessTy sess;
// this makes a copy. 

#if HAVE_MAIL_CORE
m_sessions_map[sname]=SessTy();
m_sessions_map[sname].setup_session(rwm);
#endif

}

//void to_map(ReadWriteMap & rwm) const
void cmd_fetch_test(Cip & cip , LocalVar & lv )
//static void testMessageBuilder1(String * path)
{
const StrTy ragin=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const Ragged & server=m_servers_map[ragin];
//const Ragged & params=m_ragged_map[cmd2];
ReadWriteMap rwm;
server.to_map( rwm);
#if HAVE_MAIL_CORE
StrTy sname=ragin;
//SessTy sess;
// this makes a copy. 
m_sessions_map[sname]=SessTy();
m_sessions_map[sname].run(rwm);

#endif

}



 


bool flag_bit(const IdxTy flag, const IdxTy bit, const bool pol=true)
{
const bool x=((flag&(1<<bit))!=0);
return pol?x:!x;
}


void cmd_zymo_rags(Cip & cip , LocalVar & lv )
{

const StrTy ragin=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const Ragged & qio=m_ragged_map[ragin];
const Ragged & params=m_ragged_map[cmd2];
const IdxTy nval=myatoi(cmd1);
const IdxTy maxcnt=0; // m_flp.maxcnt();
const bool output_latex=(((1<<0)&flags)!=0);
const bool  output_rank_table=(((1<<1)&flags)!=0);
const bool  output_ranks=(((1<<2)&flags)!=0)&&!output_rank_table;
const bool output_summary=(((1<<3)&flags)!=0);

} // zymo_rags
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
void cmd_dump_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy ragn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
Ragged & r=m_ragged_map[ragn];
MM_ERR(MMPR3(cmd,ragn,r.size()))
MM_ERR(MMPR(r.dump()))
//MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))
}



void cmd_transpose_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ; } // transpose
void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_read_server(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_servers_map  ) ; }


void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ; }
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


void cmd_tt(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }


void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_quit(Cip & cip , LocalVar & lv ) { clean_up(); return;  }
void cmd_join(Cip & cip , LocalVar & lv ) { wait_to_finish(); return;  }
void cmd_cm(Cip & cip , LocalVar & lv )
{ dump_cm(); return;  }
void cmd_banner(Cip & cip , LocalVar & lv )
{ config_banner(); return;  }

void cmd_set_param(Cip & cip , LocalVar & lv )
{  
 if (cip.li().cmd_ok(3))  m_flp.set(cip.li().cmd_set());
        return;  }

void cmd_get_param(Cip & cip , LocalVar & lv )
{
 if (cip.li().cmd_ok(2))
        std::cout<<lv["local_label"]<<" "<<cip.li().word(1)<<" "
                        <<m_flp.get(cip.li().word(1))<<CRLF;
         }


void cmd_verbosity(Cip & cip , LocalVar & lv )
{

const StrTy lvl=cip.p1;
IdxTy v=myatoi(lvl);
m_verbosity=v;
if (verbosity(1)) { MM_ERR(" verbosity set to "<<MMPR(m_verbosity))
}
}





void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))
} 
/*
MM_ERR()
MM_ERR()
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
*/


}

//RaggedMap m_servers_map;
//RaggedMap m_accounts_map;
//RaggedMap m_configs_map;
template <class Ty> 
void list_thing ( Ty & x, const StrTy & nm)
{
MM_LOOP(ii,x) { MM_MSG(nm<<" "<<MMPR2((*ii).first,(*ii).second.size())) }
}
template <class Ty> 
void list_thing_name ( Ty & x, const StrTy & nm)
{
MM_LOOP(ii,x) { MM_MSG(nm<<" "<<MMPR2((*ii).first,(*ii).second.name())) }
}

void cmd_list(Cip & cip , LocalVar & lv ) 
{
	list_thing(m_servers_map,"m_servers_map");
	list_thing(m_accounts_map,"m_accounts_map");
	list_thing(m_configs_map,"m_configs_map");
	list_thing(m_messages_map,"m_messages_map");
	list_thing(m_proc_map,"m_proc_map");
//MailProcMap m_proc_map;
//	list_thing_name(m_sessions_map,"m_sessions_map");
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tokragged_map) { MM_MSG("m_tokragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_dig_map) { MM_MSG("m_dig_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(MMPR2(m_messages.size(),m_thread_count))
MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();


}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("!")]=&Myt::cmd_escape;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;

m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("quit")]=&Myt::cmd_quit;
m_cmd_map[StrTy("join")]=&Myt::cmd_join;
m_cmd_map[StrTy("cm")]=&Myt::cmd_cm;
m_cmd_map[StrTy("banner")]=&Myt::cmd_banner;
m_cmd_map[StrTy("set-param")]=&Myt::cmd_set_param;
m_cmd_map[StrTy("get-param")]=&Myt::cmd_get_param;
m_cmd_map[StrTy("verbosity")]=&Myt::cmd_verbosity;




m_cmd_map[StrTy("dump-ragged")]=&Myt::cmd_dump_ragged;
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("read-server")]=&Myt::cmd_read_server;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;
m_cmd_map[StrTy("build-msg")]=&Myt::cmd_fetch_test;
m_cmd_map[StrTy("fetch-test")]=&Myt::cmd_fetch_test;
m_cmd_map[StrTy("read-msg")]=&Myt::cmd_read_msg;
m_cmd_map[StrTy("read-lib-msg")]=&Myt::cmd_read_lib_msg;
m_cmd_map[StrTy("read-folder")]=&Myt::cmd_read_folder;
m_cmd_map[StrTy("load-folder")]=&Myt::cmd_load_folder;
m_cmd_map[StrTy("start-session")]=&Myt::cmd_start_session;
m_cmd_map[StrTy("list-folders")]=&Myt::cmd_list_folders;
m_cmd_map[StrTy("download-folder")]=&Myt::cmd_download_folder;
m_cmd_map[StrTy("test-imap")]=&Myt::cmd_test_imap;
//void cmd_test_imap(Cip & cip , LocalVar & lv )
m_cmd_map[StrTy("headers")]=&Myt::cmd_headers;
#ifdef HAVE_MIKE_MUTT_SERVER
m_cmd_map[StrTy("mutt-server")]=&Myt::cmd_mutt_server;
m_cmd_map[StrTy("mutt-load")]=&Myt::cmd_mutt_load;
m_cmd_map[StrTy("mutt-proc")]=&Myt::cmd_mutt_proc;
m_cmd_map[StrTy("mutt-fck")]=&Myt::cmd_mutt_fck;
//void cmd_mutt_server(Cip & cip , LocalVar & lv )
#endif
m_cmd_map[StrTy("view")]=&Myt::cmd_view;
m_cmd_map[StrTy("view-src")]=&Myt::cmd_view_src;
m_cmd_map[StrTy("mmp-proc")]=&Myt::cmd_mmp_proc;
m_cmd_map[StrTy("mmp-config")]=&Myt::cmd_mmp_config;
m_cmd_map[StrTy("v")]=&Myt::cmd_v;
m_cmd_map[StrTy("load-mbox")]=&Myt::cmd_load_mbox;
m_cmd_map[StrTy("move-mbox")]=&Myt::cmd_cpy_mbox;
m_cmd_map[StrTy("copy-mbox")]=&Myt::cmd_cpy_mbox;
m_cmd_map[StrTy("save-att")]=&Myt::cmd_save_att;
m_cmd_map[StrTy("tree-att")]=&Myt::cmd_tree_att;
m_cmd_map[StrTy("signal")]=&Myt::cmd_signal;
m_cmd_map[StrTy("browse")]=&Myt::cmd_browse;
m_cmd_map[StrTy("stop-loads")]=&Myt::cmd_stop_loads;
m_cmd_map[StrTy("setup-session")]=&Myt::cmd_setup_session;
	;
//void cmd_build_msg(Cip & cip , LocalVar & lv )

//m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
//m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
//m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
//m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
//m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
//m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
//m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
//m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
//m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
//m_cmd_map[StrTy("zymo-rags")]=&Myt::cmd_zymo_rags;
//m_cmd_map[StrTy("group-stats")]=&Myt::cmd_group_stats;

//m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;

//m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

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


} //command_mode
// when exiting the interpretter
void clean_up()
{
if (verbosity(20)) 
{ MM_ERR(" clean up enters serial") } 
enter_serial();
m_done=true;
m_exit_bg=true;
exit_serial();
if (verbosity(20)) 
{MM_ERR(" clean up exits fitst serial")}
while (true )
{
enter_serial();
IdxTy tc=m_thread_count;
IdxTy tctotal=m_thread_total;
exit_serial();
if (tc==0){
if (verbosity(20)) 

{MM_ERR(" clean_up returns") }  return; } 
MM_ERR(" ctrl-C to stop,  waiting for conns to close "<<MMPR2(tc,tctotal))
sleep(1);
}

} 

void signal_check()
{
enter_serial();
MM_MSG(" message start entered  ")
++m_signal;
IdxTy signal=m_signal;
IdxTy threads=m_thread_count;
IdxTy total_threads=m_thread_total;
exit_serial();
MM_ERR(" signal count now "<<MMPR3(signal,threads,total_threads))
}


// this can return before the tc increases and that is not good. 
void wait_to_finish()
{
IdxTy cnt=0; 
while (true )
{
signal_check();
enter_serial();
IdxTy tc=m_thread_count;
IdxTy tctotal=m_thread_total;
exit_serial();
if (tc==0){
if (verbosity(20)) 

{MM_ERR(" wait_to_finish") }  
if ((cnt<12)&&( tctotal<3)) { MM_ERR(" tctotal too small to exit wait more ...") } else
return; } 
MM_ERR(" ctrl-C to stop, lose download no journals yet  "<<MMPR2(tc,tctotal))
sleep(10);
++cnt;
}



}



void about()
{
Ss ss;
ss<<" mjm_mikemail "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com Sun Apr 14 17:49:38 UTC 2019 "<<CRLF;
ss<<" May use libetpan, neomutt library, other sources adapted,   "<<CRLF;
ss<<" Empty skeleton code add credits like this,   "<<CRLF;
ss<<" Uses readline, libetpan, boost, and neomutt among others "<<CRLF;
ss<<"  -lreadline -lz -L /home/ubuntu/dev/libetpan/libetpan/src/.libs -letpan -lxml2 -lsasl2 -ltidy -lctemplate -lpthread -luuid -licudata -licui18n -licuio -liculx -licutest -licutu -licuuc -lssl -lcrypto -lglib-2.0 -lboost_filesystem -lboost_system -L /home/ubuntu/dev/neomutt/neomutt-20200501 -lneomutt" <<CRLF;

//ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;

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
 m_cm.dump("response_times"  , std::cout);
 MM_MSG("resp_times"<<CRLF<<m_cm.time_table("resp_times"))
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


bool verbosity( const IdxTy lvl ) const
{
return (m_verbosity>=lvl); 

}


/////////////////////////////////////////////////////////////////////
private:
void Init()
{
m_done=false;
m_verbosity=0;
MuInit();
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



void enter_serial() { m_mutex_vector.enter_serial(0); }
void exit_serial() { m_mutex_vector.exit_serial(0); }
void enter_serial2() { m_mutex_vector.enter_serial(1); }
void exit_serial2() { m_mutex_vector.exit_serial(1); }
void enter_serial(const IdxTy i ) { m_mutex_vector.enter_serial(i); }
void exit_serial(const IdxTy i ) { m_mutex_vector.exit_serial(i); }


void MuInit()
{
m_thread_count=0;
m_thread_total=0;
m_signal=0; // not strictly needed ... 
m_exit_bg=false; 
m_mutex_vector = MutexVector(4); // want to add enums for usaged

}

// MEMBERS

MutexVector m_mutex_vector;
volatile IdxTy  m_verbosity;
volatile bool m_done; // not intended for killing previously 
volatile bool m_exit_bg; // only terminate bg threads. 
volatile IdxTy m_thread_count,m_thread_total, m_signal; // In theory this can be unitted lol. 
ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
RaggedMap m_ragged_map;
RaggedMap m_servers_map;
RaggedMap m_accounts_map;
RaggedMap m_configs_map;
//SessMap m_sessions_map;
LibetMap m_libet_rdr_map;
DownloaderMap m_downloader_map;

TokRaggedMap m_tokragged_map;
MsgMapMap m_messages_map;
MyMsgMap m_messages;

TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;
FormatterTy m_formatter;
Hand m_hand;
MboxCursorMap m_cursors;
MailProcMap m_proc_map;
}; //mjm_mikemail 

/////////////////////////////////////////////////////////

#ifdef  TEST_mikemail__
int main(int argc,char **args)
{
typedef mjm_mikemail Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif
