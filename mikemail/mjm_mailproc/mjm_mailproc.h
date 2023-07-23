#ifndef MJM_MAILPROC_H__
#define MJM_MAILPROC_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
#include "mjm_collections.h"
#include "mjm_canned_methods.h"
#include "mjm_msg_rules.h"
#include "mjm_msg_jrnl.h"
// include and registero handlers 
#include "mjm_msg_responder.h"
#include "mjm_global_mailusers.h"
#include "mjm_mailproc_handler.h"
#include "mjm_mailproc_diet_handler.h"
#include "mjm_mailproc_toobib_handler.h"
#include "mjm_mailproc_csv_mmp_handler.h"
#include "mjm_mailproc_file_req_handler.h"
#include "mjm_mailproc_user_handler.h"

#include "../mjm_generic_message.h"

#include "mjm_pawnoff.h"
// order is temporary get required includes in righ file 
#include "mjm_cpp_regex.h"

#ifdef  TEST_MJM_MAILPROC

#include "../mjm_message_store.h"
#endif

/*
// Sat Jun  6 10:11:50 EDT 2020
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_mailproc   

//g++ -std=gnu++11 -DTEST_MJM_MAILPROC -I. -I.. -I../../../mjm/hlib -I../../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_mailproc.h  -lpthread -lreadline


// g++ -std=gnu++11 -DTEST_MJM_MAILPROC -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_mailproc.h  -lpthread -lreadline

*/

template <class Tr>
class mjm_mailproc 
{
 typedef mjm_mailproc Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

typedef mjm_generic_message<Tr>  MyMsg;
typedef mjm_msg_rules<Tr>  MyRules;
typedef mjm_msg_jrnl<Tr>  MyJrnl;
typedef  mjm_mailproc_handler<Tr> MyMailHand;
typedef  mjm_mailproc_diet_handler<Tr> MyMailHandDiet;
typedef  mjm_mailproc_toobib_handler<Tr> MyMailHandToobib;
typedef  mjm_mailproc_csv_mmp_handler<Tr> MyMailHandCsv;
typedef  mjm_mailproc_file_req_handler<Tr> MyMailHandFile;
typedef  mjm_mailproc_user_handler<Tr> MyMailHandUser;
typedef typename   MyMailHand::ptr MyMailHandMap;
typedef typename MyMailHand::dispatch_vars DisVars;
typedef typename MyMailHand::ctor_vars CtorVars;
//typedef typename   MyMailHand::ctor_vars Cv;
typedef std::vector< DisVars> DisVarsVec;
//typedef mjm_mailproc_handler_map<Tr> MyMailHandMap;
//typedef std::vector< MyMailHandPtr>  MyMailHandVec;
typedef typename MyMailHand::vector  MyMailHandVec;
typedef mjm_pawnoff<Tr> Hand;
typedef mjm_msg_responder<Tr> MsgResponder;
typedef mjm_global_mailusers<Tr> UserList;
typedef typename UserList::user_desc Ud;

typedef std::regex Regex;
typedef mjm_ragged_table Ragged;
typedef  Ragged::Line Rline;
//typedef MyMsg::alternative_ranks AltRank;
#ifdef  TEST_MJM_MAILPROC
typedef std::vector<MyMsg> MyMsgVector;
typedef mjm_message_store<Tr,MyMsg> MyMsgMap;
typedef std::map<StrTy, MyMsgMap> MsgMapMap;

#endif

// FIXME doh put this somwhere lol  - canned is a big junk bin 
int myatoi(const StrTy & s ) const { return mjm_canned_methods::myatoi(s.c_str()); }
int myatoi(const char * c) const { return mjm_canned_methods::myatoi(c); }


class disp_pair
{
public:
disp_pair( ): m_hand(0), m_disp(0) {} 
disp_pair(const IdxTy hand, const IdxTy disp ): m_hand(hand), m_disp(disp) {} 
IdxTy m_hand, m_disp;
StrTy dump() const { Ss ss; ss<<MMPR2(m_hand,m_disp); return ss.str(); } 
};

typedef std::map<IdxTy, disp_pair> DispPairMap;
typedef std::vector< disp_pair> DispPairVec;



static int  Register()
{

//typedef  mjm_mailproc_handler<Tr> MyMailHand;
//typedef typename   MyMailHand::ptr MyMailHandMap;
MyMailHandMap::name(StrTy("default"),  ( (MyMailHand::factory)));
MyMailHandMap::name(StrTy("default"),  MyMailHand::factoryc);

MyMailHandMap::name(StrTy("diet"),  ( (MyMailHandDiet::factory)));
MyMailHandMap::name(StrTy("diet"),  MyMailHandDiet::factoryc);
MyMailHandMap::name(StrTy("toobib"),  ( (MyMailHandToobib::factory)));
MyMailHandMap::name(StrTy("toobib"),  MyMailHandToobib::factoryc);


MyMailHandMap::name(StrTy("csvform"),  ( (MyMailHandCsv::factory)));
MyMailHandMap::name(StrTy("csvform"),  MyMailHandCsv::factoryc);

MyMailHandMap::name(StrTy("file"),  ( (MyMailHandFile::factory)));
MyMailHandMap::name(StrTy("file"),  MyMailHandFile::factoryc);

MyMailHandMap::name(StrTy("user"),  ( (MyMailHandUser::factory)));
MyMailHandMap::name(StrTy("user"),  MyMailHandUser::factoryc);


return 0; 
}

static int  _Register()
{
static int x=Register();
return x;
}

enum { MSG_MU=0 , MU_SZ};

public:
//typedef std::map<StrTy, StrTy> proc_io_type;
typedef typename MyMailHand::proc_io_type proc_io_type;

mjm_mailproc():m_mutex_vector(MU_SZ) {_Register();}
~mjm_mailproc() {}
enum { RC_NOACTION=0 };
// right now needed for listing, no meaning yet 
IdxTy size() const { return 0; }
void config_handlers(const Ragged & r) { ConfigHandlers(r,0); } 
void config_rules(const Ragged & r) { m_rules.config(r,0,m_index); } 

// Handlers need to be set before rules now for lut to work. 
void config(const Ragged & r) 
{ ConfigHandlers(r,0); m_rules.config(r,0,m_index);  } 

void respond_clear() { m_responder.clear_response_server();  }
void respond_imap() { m_responder.add_response_server(MsgResponder::MUTTINT);  }
void respond_mbox() { m_responder.add_response_server(MsgResponder::SENDMAILEXT);  }

IdxTy  proc(const MyMsg & m, const IdxTy flags) {return Proc(m,flags); } 
IdxTy  proc(proc_io_type & pio, const MyMsg & m, const IdxTy flags) 
{return Proc(pio,m,flags); } 


StrTy dump(const IdxTy flags=0) { return Dump(flags); }


#ifdef  TEST_MJM_MAILPROC

#endif

private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  
//Hand m_hand; MyRules m_rules; MyJrnl m_jrnl; MyMailHandVec m_mail_hands;
//DisVarsVec m_dis; DispPairVec m_dp_vec;
ss<<m_rules.dump()<<CRLF;
MM_LOOP(ii,m_mail_hands) {ss<<(*(*ii)).dump()<<CRLF; } 
//MM_LOOP(ii,m_dp_vec) {ss<<((*ii)).dump()<<CRLF; } 



return ss.str(); }

IdxTy Proc(const MyMsg & m, const IdxTy _flags)
{
proc_io_type pio;
return Proc(pio,m,_flags);

} // Proc 2 param 


IdxTy  Proc(proc_io_type & pio, const MyMsg & m, const IdxTy _flags) 
{
IdxTy rc=RC_NOACTION;
const IdxTy flags=0;
// the crlf needs to be handled elsewere to get the efficiency goals wtf
const bool do_nothing=Bit(_flags,0);
const MyMsg & msg=m; // m_messages[i];
const StrTy d= (msg.hval_lc("Date:"));
const StrTy f= (msg.hval_lc("From:"));
const StrTy su= (msg.hval_lc("Subject:"));
const StrTy len= (msg.hval_lc("Content-Length:"));
const StrTy from=msg.from_address(0);
StrTy udir=msg.from_address(1);
MM_ERR(MMPR4(d,f,su,len)<<MMPR(from))
DisVars u_vars;
u_vars.add("udir",udir);
// make a copy for updating etc. 
Ud  ud=m_user_list.look_up(from);
if (ud.user()=="") { ud.user(from); } 
MM_ERR(MMPR(ud.dump()))
if (true) { udir="shared"; u_vars.add("udir",udir); } 
 
//return mjm_cpp_regex::stuff::has(s,m_r,0);
EnterSerial(MSG_MU);
// in theory could have multiple hits but should not... 
const IdxTy code=m_rules.code(m,ud,0);
// create journal entry this needs hierarchy
// and maybe human interaction 
IdxTy rcs=m_jrnl.start(m,code,flags);
// unless journal says already done, 
const IdxTy rch=0;
// invoke handler
// handler could success synchronously, queue, punt, or fail 
// responses include interacting with sender and needs a facility
// for that
// this add_lut thing is probably never needed now 
if (!do_nothing) { 
m_responder.add_lut(m.server_tag(),m.server_uuid());
const IdxTy rcf=Dispatch(pio,m,ud,u_vars,flags,code);
}


// indicate journal done
IdxTy rcd=m_jrnl.done(m,rch,flags);
if (!do_nothing) m_user_list.set(ud);
ExitSerial(MSG_MU);
MM_ERR(MMPR4(rc,rcd,rcs,int(code)))
return rc;
}
IdxTy Dispatch(proc_io_type & pio, const MyMsg & m, Ud & ud, DisVars & u_vars, const IdxTy flags,const IdxTy code)
{
IdxTy rc=0;
//MM_ERR(" dispatching "<<MMPR2(code,m_dp_vec.size()))
MM_ERR(" dispatching "<<MMPR2(code,m_mail_hands.size()))
const StrTy from=m.from_address(0);
MM_ONCE(" obsolete pairs map ignored ",)
/*
if (code>=m_dp_vec.size()) return ~0; 
const disp_pair & dp=m_dp_vec[code];
if ((dp.m_hand>=m_mail_hands.size())|| (dp.m_disp>=m_dp_vec.size()))
{ MM_ERR(" bad code will bomb "<<MMPR(code)<<MMPR4(dp.m_hand, m_mail_hands.size(), dp.m_disp, m_dp_vec.size())) }
MM_ERR(" dispatching "<<MMPR(dp.m_hand))
*/
//const IdxTy sel=dp.m_disp;
const IdxTy sel=code;
//rc=m_mail_hands[dp.m_hand]().handle(m,m_dis[dp.m_disp]);
DisVars dv;
if (sel<m_dis.size()) dv=m_dis[sel];
MM_ERR(MMPR(dv.dump()))
dv+=u_vars;
MM_ERR(MMPR(dv.dump()))
//rc=(m_mail_hands[dp.m_hand])->handle(m,m_dis[dp.m_disp]);
if (sel>=m_mail_hands.size())
{
MM_ERR(" bad select code "<<MMPR2(sel,m_mail_hands.size()))
return rc; 
} // sel >=size 
rc=(m_mail_hands[sel])->handle(pio,m,ud,dv);
// may have to open a conversation or ticket


return rc;
}

void ConfigHandlers(const Rline & rline,const IdxTy n) 
{
IdxTy i=n;
const IdxTy sz=rline.size();
if (sz<=i) return;
const  StrTy nm=rline[i]; ++i;
if (sz<=i) return;
const  StrTy ty=rline[i]; ++i;
CtorVars cv(rline,3);
MyMailHand * mh=MyMailHandMap::factory(ty,cv);
if (mh==0)
{
MM_ERR(" unable to make handler no name "<<MMPR2(nm,ty))
Ss ss;
MM_LOOP(ii,rline) ss<<(*ii)<<" ";
MM_ERR("bad handler line :"<<ss.str())
return; 
}
mh->responder(&m_responder);
mh->name(nm);
if (m_index.find(nm)!=m_index.end())
{
MM_ERR(" duplicate handler name "<<MMPR3(nm,m_index[nm],m_mail_hands.size()))

}
m_index[nm]=m_mail_hands.size();
m_mail_hands.push_back(mh);

} // ConfigHandlers

void ConfigUsers(const Rline & rline,const IdxTy n) 
{
//IdxTy i=n;
const IdxTy sz=rline.size();
//if (sz<=i) return;
for(IdxTy i=n; i<sz; ++i)
{ 
const  StrTy nm=rline[i];///  ++i;
if (nm=="fn"){ ++i ;  if (i<sz) {
MM_ERR(" setting user file to "<<MMPR(rline[i])) 
 m_user_list.fn(rline[i]); } continue;  }
 
} // i 
m_user_list.load();
//CtorVars cv(rline,2);
} // ConfigUsers

bool U(StrTy & x, const Rline & l, IdxTy & i, const IdxTy sz)
{ ++i; if (i>=sz) return false; x=(l[i]); return true; }
bool U(IdxTy & x, const Rline & l, IdxTy & i, const IdxTy sz)
{ ++i; if (i>=sz) return false; x=myatoi(l[i]); return true; }



void ConfigDspv(const Rline & rline,const IdxTy n) 
{
IdxTy i=n;
const IdxTy sz=rline.size();
if (sz<=i) return;
const  StrTy nm=rline[i]; ++i;
DisVars cv(rline,2);
m_dis_index[nm]=m_dis.size();
m_dis.push_back(cv);
MM_ERR(" VARSSSSSS "<<MMPR4(nm,m_dis_index[nm],cv.dump(),n))
} // ConfigHandlers


/*
void ConfigPair(const Rline & rline,const IdxTy n) 
{
IdxTy i=n;
const IdxTy sz=rline.size();
if (sz<=(i+1)) return;
const  StrTy nm=rline[i]; ++i;
const IdxTy h=myatoi(nm);
const  StrTy nm2=rline[i]; ++i;
const IdxTy v=myatoi(nm2);
MM_ERR(MMPR4(nm,h,nm2,v))
disp_pair codepair(h,v);
m_dp_vec.push_back(codepair);

} // ConfigHandlers
*/


void ConfigHandlers(const Ragged & r,const IdxTy st) 
{ 
const IdxTy sz=r.size();
for(IdxTy i=st; i<sz; ++i)
{
const Rline & line=r.line(i);
if (line.size()<1) continue;
if (line[0]=="handler") { ConfigHandlers(line,1); } // 
//if (line[0]=="pair") { ConfigPair(line,1); } // 
if (line[0]=="vars") { ConfigDspv(line,1); } // 
if (line[0]=="users") { ConfigUsers(line,1); } // 
} // for i 

// The Dispatch Vars need to be ordered to match the handler numbers,


} 


typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

// MEMBERS members
Hand m_hand;
MyRules m_rules;
MyJrnl m_jrnl;
MyMailHandVec m_mail_hands;
typedef std::map<StrTy,IdxTy> HandlerIndex;
HandlerIndex m_index;
DisVarsVec m_dis;
HandlerIndex m_dis_index;
//DispPairVec m_dp_vec;
MsgResponder m_responder;
UserList m_user_list;
}; // mjm_mailproc

//////////////////////////////////////////////

template <class Tr>
class mjm_mailproc_map : public std::map<typename Tr::StrTy, mjm_mailproc< Tr > >  
{
 typedef mjm_mailproc_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_mailproc< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_mailproc_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
StrTy list(const IdxTy flags=0) { return List(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }
IdxTy size() const { return Super::size(); } 
private:

void Init()
{
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;
}
return ss.str();
// return Dump(flags); 
}

StrTy List(const IdxTy flags=0)
{
Ss ss;
ss<<MMPR((size()))<<CRLF;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;
}
return ss.str();
// return Dump(flags); 
}





private:

}; // mjm_mailproc_map




////////////////////////////////////////////
#ifdef  TEST_MJM_MAILPROC
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
// there is a tpedef  FCK 
//typedef Cli::list_type Choices;
typedef std::vector<StrTy> Choices;
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_mailproc<Tr> > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_MAILPROC "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_mailproc<Tr>  Myt;
//Myt x(argc,args);
Myt x;

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_MAILPROC_H__ 
