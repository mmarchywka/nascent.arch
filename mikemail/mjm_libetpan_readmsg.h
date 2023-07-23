#ifndef MJM_LIBETPAN_READMSG_H__
#define MJM_LIBETPAN_READMSG_H__

/*
Copied from the libetpan readmsg.c file for testing without MailCore 

*/

#include "mjm_globals.h"

#include "mjm_instruments.h"
#include "mjm_collections.h"


#include <mjm_ts_map.h>

#include <mjm_thread_util.h>
#include <mjm_paranoid_ofstream.h>


#include <stdexcept>
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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <stdlib.h>
#ifndef WIN32
#include <unistd.h>
#endif

#include <libetpan/charconv.h>
#include <libetpan/libetpan.h>

#include <sys/stat.h>
#ifndef WIN32
#       include <sys/mman.h>
#       include <unistd.h>
#endif
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>

//class END_MJM_STR {};
//END_MJM_STR MM_EOS;

//std::string operator<<(std::stringstream & s , const END_MJM_STR & ) { return s.str(); } 
//std::stringstream &  operator<<(std::stringstream & s , const char * c ) { ((std::ostream)s)<<c;  return s; } 

template <class Tr>
class mjm_folder_status
{
//typedef mikemail_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_folder_status Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
typedef std::vector<StrTy> Line;
public:
mjm_folder_status() { Init();}
mjm_folder_status(const StrTy & nm) : m_name(nm)  { Init();}
const StrTy & name() const { return m_name; } 
//int mailsession_status_folder(mailsession * session, const char * mb,
//    uint32_t * result_messages, uint32_t * result_recent,
//    uint32_t * result_unseen);
Myt & set ( const IdxTy m,const IdxTy r, const IdxTy u)
{
m_messages=m; m_recent=r; m_unread=u; 
return *this; }
StrTy dump() const
{
Ss ss;
const StrTy sep=" ";
ss<<MMPR2(m_name.length(),m_name)<<"      "<<sep<<MMPR3(m_messages,m_recent,m_unread);
return ss.str(); 
}
IdxTy line_size() const { return 4; } 
void line( Line & w) const 
{
w.push_back(m_name);
{ Ss ss; ss<< m_messages; w.push_back(ss.str()); } 
{ Ss ss; ss<< m_recent; w.push_back(ss.str()); } 
{ Ss ss; ss<< m_unread; w.push_back(ss.str()); } 


}

private:
void Init() 
{

m_messages=0; m_recent=0; m_unread=0; 

}

StrTy m_name;
IdxTy m_messages, m_recent, m_unread; 

}; // mjm_folder_status


template <class Tr>
class mjm_folder_downloader
{
//typedef mikemail_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_folder_downloader Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
//typedef typename Tr::Mutex Mutex;

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;
typedef mjm_ragged_table Ragged;
typedef  mjm_paranoid_ofstream<Tr> Jfs;
enum { BAD=~0}; 
enum { MU0=0, MU1=1, MU_J=2, MU_SZ=3};

public:
// this is needed because the msg number is not really a number doh 
// attempts etc could be included here but for long index and
// small machines just reading separate the download attempt history 
class folder_index
{
public:
folder_index(const IdxTy s, const IdxTy n, const StrTy & uuid)
: m_serial(s), m_number(n), m_uuid(uuid) {}
folder_index(const IdxTy s, const IdxTy n, const char *  uuid)
: m_serial(s), m_number(n), m_uuid(uuid) {}
folder_index()
: m_serial(0), m_number(0), m_uuid() {}


const IdxTy & serial() const { return m_serial; }
const IdxTy & n() const { return m_number; }
const StrTy &  uuid() const { return m_uuid; }
StrTy dump()  const
{
Ss ss;
ss<<MMPR3(m_serial,m_number,m_uuid);
return ss.str();
}

private:
// serial is reundant for consistency check and compare later.  
IdxTy m_serial, m_number;
StrTy m_uuid;

}; // folder_index
private:
class download_history 
{
public:
// n and uuid should match a folder index but serial adjustablea
// for indexing. 
download_history(const IdxTy s, const IdxTy n, const StrTy & uuid)
: m_serial(s), m_number(n), m_uuid(uuid), m_downloaded(false),m_partial(false),m_failed(0),m_state(0) {}
download_history()
: m_serial(0), m_number(0), m_uuid(), m_downloaded(false),m_partial(false),m_failed(0),m_state(0)  {}
const IdxTy & serial() const { return m_serial; }
const IdxTy & n() const { return m_number; }
const StrTy &  uuid() const { return m_uuid; }
void  set_uid(const StrTy & u )  {  m_uuid=u; }
bool did_download() const { return m_downloaded; } 
void downloaded_ok()  {  m_downloaded=true;m_partial=false; } 
void partial()  {  m_partial=!false; } 
bool  dirty()  {  return m_partial; } 
void upstate(const IdxTy s) { if (s>m_state) m_state=s; m_partial=(s==1); m_downloaded=(s==2); } 
void download_failed()  { ++m_failed;} 
void faileds(const IdxTy f )  { m_failed=f;} 
const IdxTy & fails() const   { return m_failed;} 
void state( const IdxTy s) { m_state=s; } 
StrTy dump() const
{
Ss ss;
ss<<MMPR3(m_serial,m_number,m_uuid)<<MMPR2(m_downloaded,m_failed); 
return ss.str();
}

StrTy journal() const
{
Ss ss;
ss<<MMPR3(m_serial,m_number,m_uuid)<<MMPR(m_state); 
return ss.str();
}

private:
// serial is reundant for consistency check and compare later.  
IdxTy m_serial, m_number;
StrTy m_uuid;
bool m_downloaded,m_partial; 
IdxTy m_failed;
IdxTy m_state;
}; //download_history 



public: // part of API
// this does not need to be thread safe, read one at start IIRC.
typedef mjm_ts_map<IdxTy, folder_index,Tr > FoldIdx;
//typedef std::map<IdxTy, folder_index > FoldIdx;
//typedef std::map<<StrTy, download_history > DownIdx;
typedef mjm_ts_map<StrTy, download_history,Tr > DownIdx;

private:


// Take the input ragged table and create a download history
// and maybe modify the folder index. 
void ragged(FoldIdx & f, DownIdx & d, const Ragged & s )
{
const IdxTy sz= s.size();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line &  w= s.line(i);
const IdxTy sz=w.size();
if ( sz<3)
{
MM_ERR(" wrong size "<<MMPR2(w.size(),i) );
}
const IdxTy idx=i;
const IdxTy serial= atoi(w[0].c_str());
const IdxTy n= atoi(w[1].c_str());
const StrTy&  uid=w[2];
f[idx]=folder_index(idx,n,uid);
if (sz>3) {
const char * c=w[3].c_str();
if (strncmp(c,"failed",5)==0) d[uid].faileds(atoi(c+6));
if (w[3]=="ok")  d[uid].downloaded_ok(); 
d[uid].set_uid(uid);
}
//MM_MSG(MMPR3(w[0],w[1],w[2]))
}

} // ragged
// take the journal data in s and update the loads
// verify no dirty states 
void journal(FoldIdx & f, DownIdx & d, const Ragged & s )
{
const IdxTy sz= s.size();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line &  w= s.line(i);
const IdxTy sz=w.size();
if ( sz<3)
{
MM_ERR(" wrong size "<<MMPR2(w.size(),i) );
}
const IdxTy idx=i;
const IdxTy serial= atoi(w[0].c_str());
const IdxTy n= atoi(w[1].c_str());
const StrTy&  uid=w[2];
f[idx]=folder_index(idx,n,uid);
if (sz>3) {
const char * c=w[3].c_str();
const IdxTy state=atoi(w[3].c_str());
const IdxTy DONE=2;
auto & x=d[uid];
if (state==DONE) x.downloaded_ok();
if (state>0) if (state!=DONE) x.partial();
x.upstate(state);
//if (strncmp(c,"failed",5)==0) d[uid].faileds(atoi(c+6));
//if (w[3]=="ok")  d[uid].downloaded_ok(); 
x.set_uid(uid);
}
//MM_MSG(MMPR3(w[0],w[1],w[2]))
}

} // journal 






// Convert the combined folder index and download history into
// a table ( ragged ) for editing or io. 
Ragged ragged( const FoldIdx & f, const DownIdx & d, const IdxTy flags=0)
{
Ragged r; 
MM_LOOP(ii,f)
{
const folder_index & fi= (*ii).second;
std::vector<StrTy> w;
Ss ss; 
ss<<fi.serial()<<" "; // should equal ii first. 
ss<<fi.n()<<" ";
ss<<fi.uuid()<<" ";
auto jj= d.find(fi.uuid());
if (jj!=d.end())
{
const auto & dh=(*jj).second; 
 if (dh.did_download()) ss<<"ok "; 
 if (dh.fails()!=0) ss<<"failed"<<dh.fails()<<" "; 

}
// letting it reparse is dumb but need to conver int to str 
r.add(ss.str());
}
return r;
}
void Init()
{
m_idx_ptr=BAD;
m_count=0; 
m_mutex_vector = MutexVector(MU_SZ);
m_idx_invalid=true; 
m_advance_on_result=false;
m_verbosity_ptr=0;
}



void enter_serial(const IdxTy n) { m_mutex_vector.enter_serial(n); }
void exit_serial(const IdxTy n) { m_mutex_vector.exit_serial(n); }

MutexVector m_mutex_vector; 



FoldIdx m_idx;
DownIdx m_down;
IdxTy m_idx_ptr,m_count;
bool  m_idx_invalid;
bool m_advance_on_result;
std::map<StrTy , IdxTy> m_uid_fails;
std::map<IdxTy , IdxTy> m_q_histo;
volatile IdxTy * m_verbosity_ptr;
StrTy m_journal;
//IdxTy m_jflag;

public:
Myt & set_verbosity( IdxTy * p ) {  m_verbosity_ptr=p; return *this;   } 
Myt & set_verbosity( volatile  IdxTy & p ) {  m_verbosity_ptr=& p; return *this;   } 

void verbosity( const IdxTy * p) { m_verbosity_ptr=p; } 
bool verbosity(const int  n) { return verbosity((IdxTy)n, true); } 
//bool verbosity(const int  n) { return verbosity((IdxTy)n, def); } 

bool verbosity(const IdxTy n,const bool def=true) 
{
if ( m_verbosity_ptr==0 ) return def;
return  (*m_verbosity_ptr>=n); 
}
void journal(const StrTy & fn )  { m_journal=fn; } 


mjm_folder_downloader() {Init(); }
~mjm_folder_downloader() {

}
FoldIdx * idx() { return &m_idx; } 
void load_historical_data(const StrTy & fn)
{
Ragged r;
r.load(fn); // use space sep and no quote or escape - need to specifiy alt 
ragged(m_idx,m_down,r);


}
void save_historical_data(const StrTy & fn)
{
Ragged r= ragged(m_idx,m_down,0);
r.write_file(fn,3); // space sep  disable serial numbers 

}
void journal( const download_history & dh)
{
if (m_journal.c_str()[0]==0) return; 
// TODO FIXME never hold a lock duing IO... 
enter_serial(MU_J);
Jfs jfs(m_journal);
jfs<<dh.journal()<<CRLF;
exit_serial(MU_J);
}
// this really needs to be done by the caller, bad location 
// just make a list of things to download???? 
void download_mailbox(const StrTy & fn)
{

FoldIdx oldf;
DownIdx oldd;
Ragged r;
r.load(fn); // use space sep and no quote or escape - need to specifiy alt 
//void ragged(FoldIdx & f, DownIdx & d, const Ragged & s )
m_count=m_idx.size();
ragged(oldf,oldd,r);
//if ( m_idx_invalid) { get_index(); }
m_idx_invalid=(m_idx.size()==0);  // TODO FIXME see line below 
if ( m_idx_invalid||(m_idx.size()==0) ) { m_idx_ptr=~0; return;  }
m_down=oldd; // .clear();
m_idx_ptr=0;
//f[idx]=folder_index(idx,n,uid);
reconcile_past();
next_nonloaded_msg();
} // download_mailbox
bool more_to_download(const IdxTy flags=0)
{
enter_serial(0);
IdxTy i=m_idx_ptr;
exit_serial(0);
// make sure valid to look for next
if ( i== ~0) return false; 
return true; 
}
StrTy status()
{
Ss ss;
// should just chunk these up into locals and then stream 
enter_serial(0);
const IdxTy sz=m_uid_fails.size();
const IdxTy ptr=m_idx_ptr;
StrTy ui[sz]; IdxTy cnt[sz];
IdxTy i=0;
MM_LOOP(ii,m_uid_fails)
{ ui[i]=(*ii).first; cnt[i]=(*ii).second;
++i;
}
exit_serial(0);
ss<<MMPR(ptr);
for(IdxTy i=0; i<sz; ++i) ss<<"uuid="<<ui[i]<<" fails "<<cnt[i]<<";";
return ss.str();
}
const StrTy uuid(const IdxTy n)
{
const folder_index & fi =m_idx[n];
const StrTy & iud=fi.uuid();
//MM_ERR(" uuid in folder downloader "<<MMPR2(iud,n))
return iud;
}

StrTy message_data( const IdxTy n )
{
const folder_index & fi =m_idx[n];
const StrTy & iud=fi.uuid();
Ss ss;
//#error this needs toi be synchronized rwm on map doh 
ss<<fi.dump()<<m_down[iud].dump();
return ss.str(); // fi.dump();

}
// this needs to be synchronized somwhwere.. 
void last_message_started( const IdxTy n, const IdxTy state=1)
{
if (m_journal.length()==0) return;
const folder_index & fi =m_idx[n];
const StrTy & iud=fi.uuid();
enter_serial(0);
auto & x= m_down[iud];
x.set_uid(iud);
x.upstate(state);
//x.partial();
// and then this gets ANOTHER lock... 
journal(x);
exit_serial(0);

}


void last_message_quality(const IdxTy q, const IdxTy n)
{
// accept it as good 
//const folder_index & fi =m_idx[m_idx_ptr];
const folder_index & fi =m_idx[n];
const StrTy & iud=fi.uuid();
enter_serial(0);
++m_q_histo[q];
if (q!=0)
{
++m_uid_fails[iud];
}
m_down[iud].set_uid(iud);
// this is stupid, the default ctor will make the instance uuid inconsistent with
// accurate map key, need to fix TODO FIXME 
if ( q==0) {  m_down[iud].downloaded_ok();
// TODO FIXME these will deadlock if they need to hold locks during parameter evaluation
// and the err line lock is on 
if ( verbosity(20)) {MM_ERR(" should indicate success "<<MMPR3(fi.uuid(),m_down.size(), m_down[fi.uuid()].did_download())) } 
 }
else  {  m_down[iud].download_failed(); 
if ( verbosity(1)) { MM_ERR(" should indicate FAIL "<<MMPR3(fi.uuid(),m_down.size(), m_down[fi.uuid()].did_download())) } 

}
// note that this need not be sequential for this thread.
if ( m_advance_on_result) { // =false;
if (m_idx_ptr!=BAD) ++m_idx_ptr; // do not try again ... 
 next_nonloaded_msg(); } 
exit_serial(0);
}
// now this needs to increment on launch not result. 
IdxTy next_message() {

enter_serial(0);
if ( m_advance_on_result){ IdxTy p=m_idx_ptr; exit_serial(0);   return p; } 
const IdxTy p=m_idx_ptr;
 next_nonloaded_msg();  
// don't want to wait for that lol
if (m_idx_ptr!=BAD)  ++m_idx_ptr;
if (m_idx_ptr>=m_count) m_idx_ptr=BAD; 
exit_serial(0);
const StrTy puuid=m_idx[p].uuid();
if (puuid.c_str()[0]==0) 
{
MM_ERR("blank uuid "<<MMPR3(p,m_idx_ptr,m_idx[p].n()));
}
return p; 
 } 
void reconcile_past()
{
IdxTy already_done=0;
IdxTy down_entries=0;
// probably faster to just scan the m_down entries
// however the m_idx is not indexed by uuid. 
if (verbosity(20)) { MM_ERR(MMPR2(m_down.size(),m_idx.size())) } 
for (IdxTy i=0; i<m_count; ++i)
{
const folder_index & fi =m_idx[i];
auto ff=m_down.find(fi.uuid()); 
if (ff!=m_down.end())
{
//MM_ERR(MMPR3(fi.uuid(),(*ff).first,(*ff).second.uuid()))
++down_entries;
if ( m_down[fi.uuid()].did_download() )   
{ ++already_done; } // did_download
}
} // i 

if (verbosity(0)) { MM_ERR(MMPR4(m_down.size(),m_count,already_done,down_entries))  } 
} // reconcile_past

// not always serialized although this is slow should be done without a lock. ... 
IdxTy next_nonloaded_msg()
{
// actually count could equal this but realistically not... 
if (m_idx_ptr==BAD) return BAD;
while ( m_idx_ptr<m_count )
{
if (m_idx.find(m_idx_ptr)==m_idx.end()) break;
const folder_index & fi =m_idx[m_idx_ptr];
auto ff=m_down.find(fi.uuid()); 
if (ff!=m_down.end())
{
if ( m_down[fi.uuid()].did_download() )   
{ ++m_idx_ptr; continue;
} // did_download
} // found
// found one.
// returning member ? on stac I guess  
//m_n=m_idx_ptr;
return m_idx_ptr; 
} // while  m_idx_ptr
m_idx_ptr=BAD; 
return m_idx_ptr; 
} // next_nonloaded_msg




}; // mjm_folder_downloader



////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <class Tr>
class mjm_libetpan_readmsg
{
//typedef mikemail_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_libetpan_readmsg Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy,StrTy> FlagMap;

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

typedef  mjm_folder_downloader<Tr> FolderDownloader;
//typedef  typename  mjm_folder_downloader<Tr>::folder_index folder_index;
typedef  typename  FolderDownloader::folder_index folder_index;

//typedef pthread_mutex_t Mutex;
//typedef Tr::Mutex Mutex;

//typedef std::map<IdxTy, folder_index> FoldIdx;
typedef typename FolderDownloader::FoldIdx FoldIdx; 

public:
typedef mjm_folder_status<Tr> FolderStatus;
typedef std::vector<FolderStatus> FolderList;

#if 0 

#endif




mjm_libetpan_readmsg() { Init(); } 
mjm_libetpan_readmsg( IdxTy * p ) { Init(); m_verbosity_ptr=p;  } 



~mjm_libetpan_readmsg() 
{
// probably wrong heap.. 
//delete m_storage;
if ( m_storage!=0) {mailstorage_free(m_storage); m_storage=0; } 
if ( m_f!=0) {mailfolder_free(m_f); m_f=0; } 
//delete m_f;
}
enum { BAD=~0}; 
enum { FETCH_HEADER=(1<<0), FETCH_BODY=(1<<1), FETCH_VERBATIM=(1<<2), FETCH_SPLIT=(1<<3), FETCH_FLAGS=(1<<4) }  ; 
const StrTy& header() const { return m_hdr;}
const StrTy& body() const { return m_body;}
const StrTy& verbatim() const { return m_verbatim;}
const IdxTy  bit_flags() const { return m_bit_flags; } 
const FlagMap &  string_flags() const { return m_string_flags; } 
const StrTy  string_flags_string() const 
{StrTy x;  StrTy sep="";  MM_LOOP(ii, m_string_flags) {x=x+sep+(*ii).first;  sep=" ";   } return x;  } 
// TODO want to replace with uuid
void set_message(const IdxTy n ) { m_n=n; }
Myt & set_verbosity( IdxTy * p ) {  m_verbosity_ptr=p; return *this;   } 
Myt & set_verbosity( volatile  IdxTy & p ) {  m_verbosity_ptr=& p; return *this;   } 
void disconnect() { Disconnect(); } 
void fetch_flags(const IdxTy n ) { m_fetch_flags=n; }
void get_message() { GetAMessage(); }
bool worth_trying()  {if (m_r!=0) { m_worth_trying=false; }  return m_worth_trying; } 
// TODO generally recognized as a problem but ok for now, uuid should be better
// all the info is here for that. 
IdxTy message_number() { return m_n; }
//const StrTy a=x.account(); //  x.m_user+StrTy(";")+host+StrTy(";")+folder;
StrTy account(){ return  m_user+StrTy(";")+m_host+StrTy(";")+m_folder; } 
void set_index(FoldIdx * p ) {  m_idxp=p; } 
// should take an account info thing...
void set_account_info(const ReadWriteMap & rwm,const StrTy folder)
{
m_host=rwm.get_string("hostname","missing"); //
const StrTy sport=rwm.get_string("port","missing");
m_port= atoi(sport.c_str());
 m_user=(rwm.get_string("user","missing"));
 m_pass=rwm.get_string("password","missing");
m_folder=folder;
m_path=m_folder;
MM_ERR(MMPR4(m_host,m_port,m_user,m_pass))
}

const IdxTy count() const { return m_count; } 

// flags : 0=return without getting list 
// flags : 1=get envelopes 

void get_index(const IdxTy flags=0)
{
m_idx_invalid=true;
if ( m_idxp==0) set_index(m_fd.idx());
//m_idx.clear();
m_idxp->clear();
InitStorage();
if (Die(m_r," init storage fails in get_index")) { return ; }
//MM_ERR(" storage worked next folders ")
InitFolders();
if (Die(m_r," init folders fails in get_index")) { return ; }
//if (m_r!=0) { return ; }
MM_ERR(" getting index   ")
if (Bits(flags,1)) return;
if (m_idx_invalid) GetList((flags>>1));
if ( m_idx_invalid)
{ MM_ERR(" folder index did not initialize") } 

}

// need to write copy ctor etc 
void post_assignment_fix()
{
m_storage=0;
m_f=0;
m_worth_trying=false;
InitStorage();
if (Die(m_r, " reallocate storag ",2)) return;
InitFolders();
if (Die(m_r, " reallocate folder ",2)) return;

}
void all_at_once(const StrTy & user, const StrTy & pass, const StrTy & host, const IdxTy port,
IdxTy n, const StrTy & folder, const IdxTy flags )
{
m_user=user;
m_pass=pass;
m_host=host;
m_folder=folder;
m_path=m_folder;
m_port=port;
m_n=n;
InitStorage();
if (m_r!=0) { return ; }
if (verbosity(10)) { MM_ERR(" storage worked next folders ") } 
InitFolders();
if (m_r!=0) { return ; }
if ( verbosity(10)) { MM_ERR(" connection my have worked  ") } 

if (m_idx_invalid) GetList(flags);
if ( m_idx_invalid)
{ MM_ERR(" folder index did not initialize") } 
else GetAMessage();

}


FolderList*  get_folder_status(const StrTy & root)
{
if (m_f==0){  if (m_r==0) m_r=BAD; Die(m_r," m_F is null" ) ; return 0 ; } 
if (m_f->fld_session==0){  if (m_r==0) m_r=BAD; Die(m_r," m_F->session  is null" ) ; return 0 ; } 
if (  Die(m_r, "no valid conn for folder status"  )) return 0 ; 
//m_r = mailsession_get_messages_list(m_f->fld_session, &res);
//int mailsession_status_folder(mailsession * session, const char * mb,
//    uint32_t * result_messages, uint32_t * result_recent,
//    uint32_t * result_unseen);
mail_list * res;
//int mailsession_list_folders(mailsession * session, const char * mb,
//                             struct mail_list ** result);
MM_ERR(MMPR(root))
m_r=mailsession_list_folders(m_f->fld_session, root.c_str(),&res);
if (Die(m_r, " unable to get folder list ")) return 0 ; 
//struct mail_list { clist * mb_list; /* elements are (char *) */
FolderList&  flist = *( new FolderList());

clist * pf=res->mb_list; // clist *
clistiter* ii= clist_begin(pf);
clistiter* ie= clist_end(pf);
while (ii!=ie)
{
const char * f=(const char *)clist_content(ii);
if ( f!=0){
flist.push_back(FolderStatus(f));
// MM_STATUS(MMPR(f)<<"       ")
 MM_ERR(MMPR(f)<<"       ")
}
else { MM_ERR(" folder mb name ") } 
ii=clist_next(ii);
}
mail_list_free(res);
MM_LOOP(ii,flist)
{
FolderStatus & f=(*ii);
uint32_t m,r,u;
const char * mb=f.name().c_str();
 MM_ERR(" getting "<<MMPR2(strlen(mb),mb)<<"   FUC " )
IdxTy rr= mailsession_status_folder(m_f->fld_session,  mb, &m, &r,&u); 
Ss ss;
//if (Die(m_r,(ss<<" getting folder info "<<MM_EOS))) return; 
//if (Die(m_r," getting folder info ")) return; 
//Warn(rr, " getting folder info"); 
f.set(m,r,u);
//MM_STATUS(f.dump()<<"                    " )
// the fcking status was fcking up the fcking output FCK... 
MM_ERR(f.dump()<<"                    " )
}
//int mailsession_get_envelopes_list(mailsession * session,
//                                   struct mailmessage_list * result);

return & flist;
} // get_folder_status


Myt &  operator=(const Myt & x) 
{
//return mjm_libetpan_readmsg(x); 
assign(x);
return *this;
}
mjm_libetpan_readmsg(const Myt & x)//  { Init(); (*this)=x;  } 
//Myt  operator=(const Myt & x) 
{
assign(x); 
}

void assign(const Myt & x) 
{
// free our stuff
//if ( m_storage!=0) { mailstorage_free(m_storage); m_storage=0; } 
//if ( m_f!=0) { mailfolder_free(m_f); m_f=0; } 
m_storage=0;
m_f=0;
m_worth_trying=false;
m_verbosity_ptr=x.m_verbosity_ptr;
m_user=x.m_user;
m_pass=x.m_pass;
m_host=x.m_host;
m_folder=x.m_folder;
m_port=x.m_port;
m_fetch_flags=x.m_fetch_flags;
m_n=x.m_n;
m_r=x.m_r;
m_conn=x.m_conn;
m_auth=x.m_auth;
m_xoauth2=x.m_xoauth2;
m_path=x.m_path;
m_cache_dir=x.m_cache_dir;
m_flags_dir= x.m_flags_dir;
m_idxp=x.m_idxp;
m_count=x.m_count;
m_idx_invalid=x.m_idx_invalid;
m_hdr=x.m_hdr;
m_body=x.m_body;
m_verbatim=x.m_verbatim;
m_bit_flags=x.m_bit_flags;
m_string_flags=x.m_string_flags;
//return *this;
}

StrTy dump()
{
Ss ss;
ss<<MMPR4(m_user,m_pass,m_host,m_folder);
ss<<MMPR4(m_port,m_fetch_flags,m_n,m_r);
ss<<MMPR4(m_conn,m_auth,m_xoauth2,m_path);
ss<<MMPR2(m_count,m_idx_invalid);

return ss.str();
}
void verbosity( const IdxTy * p) { m_verbosity_ptr=p; } 
bool verbosity(const int  n) { return verbosity((IdxTy)n, true); } 
//bool verbosity(const int  n) { return verbosity((IdxTy)n, def); } 

bool verbosity(const IdxTy n,const bool def=true) 
{
if ( m_verbosity_ptr==0 ) return def;
return  (*m_verbosity_ptr>=n); 
}

private:

mailstorage * m_storage;
mailfolder * m_f;
bool m_worth_trying;
volatile IdxTy * m_verbosity_ptr;
StrTy  m_user,m_pass, m_host, m_folder; 
IdxTy  m_port,m_fetch_flags,m_n,m_r,m_conn,m_auth;
bool m_xoauth2;
StrTy m_path,m_cache_dir, m_flags_dir;
FolderDownloader m_fd; // putting this into the map requires synchronization. 
//x.set_index(fd.idx());

FoldIdx*  m_idxp;
#if 0 
FoldIdx m_idx;
DownIdx m_down;
IdxTy m_idx_ptr;
#endif

//IdxTy m_count,m_idx_ptr;
IdxTy m_count;
bool m_idx_invalid;
StrTy m_hdr,m_body,m_verbatim;
IdxTy m_bit_flags;
FlagMap m_string_flags;
// copied from the libetpan examples 
#include <libetpan-readmsg-common.h>



void Init()
{
m_idxp=0;
m_storage=0;
m_f=0; 
m_worth_trying=false;
m_verbosity_ptr=0;
//m_folder=0;
m_conn=CONNECTION_TYPE_TRY_STARTTLS;
m_conn=CONNECTION_TYPE_TLS;
m_auth=0;
m_port=0;
m_auth=IMAP_AUTH_TYPE_PLAIN;
m_n=0;
m_count=0;
#if 0
m_idx_ptr=~0;
#endif

m_r=0;
m_xoauth2=!true;
m_path=m_folder; // "+INBOX";
m_cache_dir=".";
m_flags_dir=".";
//m_fetch_flags=FETCH_VERBATIM|FETCH_SPLIT;
m_fetch_flags=FETCH_VERBATIM|FETCH_SPLIT|FETCH_FLAGS;
m_idx_invalid=true;
//m_idxp->clear();
InitMsg();
}
void InitStorage()
{
if ( m_storage!=0) { mailstorage_free(m_storage); m_storage=0; } 
//if ( m_f!=0) mailfolder_free(m_f); m_f=0; } 
//delete m_f;
  m_storage = mailstorage_new(NULL);
  if (m_storage == NULL) { m_r=1;
     printf("error initializing storage\n"); rte(" storage failes to init");        } 
else
{
/* the 1 is for imap, 
static struct storage_name storage_tab[] = {
  {POP3_STORAGE, "pop3"},
  {IMAP_STORAGE, "imap"},
  {NNTP_STORAGE, "nntp"},
  {MBOX_STORAGE, "mbox"},
  {MH_STORAGE, "mh"},
  {MAILDIR_STORAGE, "maildir"},
  {FEED_STORAGE, "feed"},
};

this works, 
 ./readmsg -s imap.mail.yahoo.com -p 993 -u marchywka@yahoo.com -v "Generic3!
" -d imap -t -l INBOX 10



*/
//m_r=init_storage(m_storage,1,m_host.c_str(),m_port, m_conn,m_user.c_str(),
//m_pass.c_str(), m_auth, m_xoauth2,m_path.c_str(),m_cache_dir.c_str(),m_flags_dir.c_str());
//m_pass.c_str(), m_auth, m_xoauth2,m_path.c_str(),0,0);
//m_pass.c_str(), m_auth, m_xoauth2,0,0,0);

 if (m_xoauth2) {
      m_r = imap_mailstorage_init_sasl(m_storage, m_host.c_str(), m_port, m_path.c_str(), m_conn,
          "xoauth2", NULL, NULL, NULL, NULL, m_user.c_str(), m_pass.c_str(), NULL, false, NULL);
    } else {
      m_r = imap_mailstorage_init(m_storage, m_host.c_str(), m_port, m_path.c_str(), m_conn,
      //m_r = imap_mailstorage_init(m_storage, m_host.c_str(), m_port, 0, m_conn,
          IMAP_AUTH_TYPE_PLAIN, m_user.c_str(), m_pass.c_str(), false,NULL);
    }
// this should throw but return anyway in case 
if (  Die(m_r, "unable to connec"  )) return; 

}
}
void rte(const char * s)
{
// 2023-05-30 new compiler 
//throw std::runtime_error::runtime_error(s); 
throw std::runtime_error(s); 
}
void InitFolders()
{
//if ( m_storage!=0) mailstorage_free(m_storage); m_storage=0; } 
if ( m_f!=0) { mailfolder_free(m_f); m_f=0; } 
//delete m_f;
  m_f = mailfolder_new(m_storage, m_path.c_str(), NULL);
  //m_f = mailfolder_new(m_storage, NULL, NULL);
  if (m_f == NULL) {
//m_r=1;
    if  ( verbosity(1)) { MM_ERR("error initializing folder"<<MMPR(m_r)) } 
	rte("error initializing folder");
    return; // goto free_storage;
  }

  m_r = mailfolder_connect(m_f);
  if (m_r != MAIL_NO_ERROR) {
 if ( verbosity(1)) {    MM_ERR("error initializing folder2"<<MMPR(m_r)) } 
	rte("error initializing folder2");
	return; 
    //goto free_folder;
  }


	if (m_f->fld_session==0) {m_r=999;  Die(m_r,"null session"); return ; } 
m_worth_trying=true;
SetProgressCallback();
} // InitFolders

void Disconnect()
{
m_worth_trying=false;
if (m_f==NULL ) return;
  m_r = mailfolder_disconnect(m_f);

 if ( verbosity(10)) {    MM_ERR("disconnected "<<MMPR2(m_r,account())) } 
if (Die(m_r," disconnecting")) return;
//  if (m_r != MAIL_NO_ERROR) {
// if ( verbosity(1)) {    MM_ERR("error initializing folder2"<<MMPR(m_r)) } 
//	rte("error initializing folder2");
//	return; 
    //goto free_folder;
  //}
}


#if 0

struct mailimf_fields {
  clist * fld_list; /* list of (struct mailimf_field *), != NULL */
};

struct mailmessage {
  mailsession * msg_session;
  mailmessage_driver * msg_driver;
  uint32_t msg_index;
  char * msg_uid;

  size_t msg_size;
  struct mailimf_fields * msg_fields;
  struct mail_flags * msg_flags;

  int msg_resolved;
  struct mailimf_single_fields msg_single_fields;
  struct mailmime * msg_mime;

  /* internal data */

  int msg_cached;
  void * msg_data;

 /*
   msg_folder field :
   used to reference the mailfolder, this is a workaround due
   to the problem with initial conception, where folder notion
   did not exist.
 */
  void * msg_folder;
  /* user data */
  void * msg_user_data;
};
struct mailimf_field {
  int fld_type;
  union {
    struct mailimf_return * fld_return_path;              /* can be NULL */
    struct mailimf_orig_date * fld_resent_date;    /* can be NULL */
    struct mailimf_from * fld_resent_from;         /* can be NULL */
    struct mailimf_sender * fld_resent_sender;     /* can be NULL */
    struct mailimf_to * fld_resent_to;             /* can be NULL */
    struct mailimf_cc * fld_resent_cc;             /* can be NULL */
    struct mailimf_bcc * fld_resent_bcc;           /* can be NULL */
    struct mailimf_message_id * fld_resent_msg_id; /* can be NULL */
    struct mailimf_orig_date * fld_orig_date;             /* can be NULL */
    struct mailimf_from * fld_from;                       /* can be NULL */
    struct mailimf_sender * fld_sender;                   /* can be NULL */
    struct mailimf_reply_to * fld_reply_to;               /* can be NULL */
    struct mailimf_to * fld_to;                           /* can be NULL */
    struct mailimf_cc * fld_cc;                           /* can be NULL */
    struct mailimf_bcc * fld_bcc;                         /* can be NULL */
    struct mailimf_message_id * fld_message_id;           /* can be NULL */
    struct mailimf_in_reply_to * fld_in_reply_to;         /* can be NULL */
    struct mailimf_references * fld_references;           /* can be NULL */
    struct mailimf_subject * fld_subject;                 /* can be NULL */
    struct mailimf_comments * fld_comments;               /* can be NULL */
    struct mailimf_keywords * fld_keywords;               /* can be NULL */
    struct mailimf_optional_field * fld_optional_field;   /* can be NULL */
  } fld_data;
};
enum {
  MAILIMF_FIELD_NONE,           /* on parse error */
  MAILIMF_FIELD_RETURN_PATH,    /* Return-Path */
  MAILIMF_FIELD_RESENT_DATE,    /* Resent-Date */
  MAILIMF_FIELD_RESENT_FROM,    /* Resent-From */
  MAILIMF_FIELD_RESENT_SENDER,  /* Resent-Sender */
  MAILIMF_FIELD_RESENT_TO,      /* Resent-To */
  MAILIMF_FIELD_RESENT_CC,      /* Resent-Cc */
  MAILIMF_FIELD_RESENT_BCC,     /* Resent-Bcc */
  MAILIMF_FIELD_RESENT_MSG_ID,  /* Resent-Message-ID */
  MAILIMF_FIELD_ORIG_DATE,      /* Date */
  MAILIMF_FIELD_FROM,           /* From */
  MAILIMF_FIELD_SENDER,         /* Sender */
  MAILIMF_FIELD_REPLY_TO,       /* Reply-To */
  MAILIMF_FIELD_TO,             /* To */
  MAILIMF_FIELD_CC,             /* Cc */
  MAILIMF_FIELD_BCC,            /* Bcc */
  MAILIMF_FIELD_MESSAGE_ID,     /* Message-ID */
  MAILIMF_FIELD_IN_REPLY_TO,    /* In-Reply-To */
  MAILIMF_FIELD_REFERENCES,     /* References */
  MAILIMF_FIELD_SUBJECT,        /* Subject */
  MAILIMF_FIELD_COMMENTS,       /* Comments */
  MAILIMF_FIELD_KEYWORDS,       /* Keywords */
  MAILIMF_FIELD_OPTIONAL_FIELD  /* other field */
};

#endif
// void (*)(size_t, size_t, void*) 

static void Handler(mailimap_msg_att * att,  void * context)
{ Myt & x=(Myt&) context;
//MM_ERR(MMPR2(l1,l2))
mailimap_msg_att_free(att);
x.Handler();
}


static void BodyProgress(size_t l1, size_t l2, void * context)
{ Myt & x=(Myt&) context;
MM_ERR(MMPR2(l1,l2))
x.BodyProgress();
}

static void ItemProgress(size_t l1, size_t l2, void * context)

{ Myt & x=(Myt&) context;
MM_ERR(MMPR2(l1,l2))
x.ItemProgress(); 
}

void BodyProgress()
{
const char * c=mailimap_read_line((mailimap*)m_f->fld_session);
MM_ERR(" body "<<MMPR(c))
}
void ItemProgress()
{
const char * c=mailimap_read_line((mailimap*)m_f->fld_session);
MM_ERR(" item "<<MMPR(c))

}
void Handler()
{
const char * c=mailimap_read_line((mailimap*)m_f->fld_session);
MM_ERR(" hand "<<MMPR(c))

}

// probably need to look at mailstream in libetpan
void SetProgressCallback()
{
if (false) { 
mailimap_set_progress_callback((mailimap*)m_f->fld_session,  ( &Myt::BodyProgress), ( &Myt::ItemProgress),this);
mailimap_set_progress_callback((mailimap*)m_storage->sto_session,  ( &Myt::BodyProgress), ( &Myt::ItemProgress),this);
mailimap_set_msg_att_handler((mailimap*)m_f->fld_session,(&Myt::Handler), this);
mailimap_set_msg_att_handler((mailimap*)m_storage->sto_session,(&Myt::Handler), this);
} 
//m_r=r;
//if (Die(m_r," setting callback ")) return; 

}

void GetSessEnvelopes(mailsession * sess)
{
mailmessage_list  *res; // kk =mailmessage_list_new(0);
int afck_r = mailsession_get_messages_list(sess, &res);
MM_ERR(" sht fck ")
int _r = mailsession_get_envelopes_list(sess, res);
    if (_r != MAIL_NO_ERROR) {
      MM_ERR("** message  not found ** - "<< MMPR2(_r, maildriver_strerror(_r)))
      return; 
} // error 
MM_ERR(" ass sht fcker  ")
carray * p= res->msg_tab;
IdxTy cnt= p->len;
m_count=cnt;
MM_ERR(" ass adcacascad  ")
MM_ERR("envelope count "<< MMPR(cnt))
  for(IdxTy i = 0 ; i < carray_count(p) ; i ++) {
    mailmessage * msg;
    msg = (mailmessage*)carray_get(p, i);
 char *  uid=msg->msg_uid;
IdxTy idx=msg->msg_index;
const int flags= msg->msg_flags->fl_flags;
const char * subj="none";
const char * mess="none";
StrTy sdate;
const char * odate="none";
  //struct mailimf_fields * msg_fields;
   mailimf_fields * mf=msg->msg_fields;
clist * cl= mf->fld_list;;
clistiter * ci=clist_begin(cl);
while (ci!=clist_end(cl))
{
//MM_ERR(" FUKSIJNG SHT SHTIS FCK SHT")
mailimf_field * pf=(mailimf_field*)clist_content(ci); 
const char * doh =**(char***)((&pf->fld_data)); 
//if (pf->fld_type==MAILIMF_FIELD_MESSAGE_ID) 
if (pf->fld_type==MAILIMF_FIELD_ORIG_DATE) {sdate=dt_format(doh);  odate= sdate.c_str(); } 
else if (pf->fld_type==MAILIMF_FIELD_MESSAGE_ID) { mess = doh; } 
else if (pf->fld_type==MAILIMF_FIELD_SUBJECT) 
{ 
subj =**(char***)((&pf->fld_data)); 
// AFCK SHT FCK 
//char * afck =((mailimf_message_id*)(pf->fld_data))->mid_value; 
//char * afck =pf->mailimf_message_id-> mid_value; 
//char * afck =pf->fld_data.mailimf_message_id->mid_value; 
//char * afck =**(char**)(void**)(&(pf->fld_data)); 
//char * afck =**(char***)((&pf->fld_data)); 
//MM_ERR(MMPR3(msg->msg_index,msg->msg_uid,afck)) 
}
ci=clist_next(ci);
} // ci 
MM_ERR(MMPR4(idx,uid,flags,subj)<<MMPR2(mess,odate))
} // for i 
// not sure this is right as res is an auto does it free mem?
mailmessage_list_free(res);
} // GetSessEnvelopes

StrTy dt_format(const void * _dt)
{
const mailimf_date_time * dt=(const mailimf_date_time *)_dt;
Ss ss;
ss<<(dt->dt_year)<<"-";
ss<<(dt->dt_month)<<"-";
ss<<(dt->dt_day)<<":";
ss<<(dt->dt_hour)<<":";
ss<<(dt->dt_min)<<":";
ss<<(dt->dt_sec)<<"+";
ss<<(dt->dt_zone);

return ss.str();
}
/*struct mailimf_date_time {
  int dt_day;
  int dt_month;
  int dt_year;
  int dt_hour;
  int dt_min;
  int dt_sec;
  int dt_zone;
};

*/

///////////////////////////////////////////////////
#ifdef stupid_afck_sht
void GetEnvelopes()
{

mailmessage_list* res;
m_r = mailsession_get_envelopes_list(m_f->fld_session, &res);
    if (m_r != MAIL_NO_ERROR) {
      MM_ERR("** message  not found ** - "<< MMPR2(m_n, maildriver_strerror(m_r)))
      return; ;
}
carray * p= res[0].msg_tab;
IdxTy cnt= p->len;
m_count=cnt;
MM_ERR("envelope count "<< MMPR(cnt))
  for(IdxTy i = 0 ; i < carray_count(p) ; i ++) {

/*

struct mail_flags {
  uint32_t fl_flags;
  clist * fl_extension; //
};
struct mailmessage {
  mailsession * msg_session;
  mailmessage_driver * msg_driver;
  uint32_t msg_index;
  char * msg_uid;
  size_t msg_size;
  struct mailimf_fields * msg_fields;
  struct mail_flags * msg_flags;
  int msg_resolved;
  struct mailimf_single_fields msg_single_fields;
  struct mailmime * msg_mime;
  int msg_cached;
  void * msg_data;
*/


    mailmessage * msg;
    msg = (mailmessage*)carray_get(p, i);
const int flags= msg->msg_flags->fl_flags;
	(*m_idxp)[i]=folder_index(i,msg->msg_index,msg->msg_uid);
//	MM_ERR(MMPR4(i,msg->msg_index,msg->msg_uid,msg->msg_size))
}

mailmessage_list_free(res);



} // GetEnvelopes
#endif
public:
       // MMPR2(m_n, maildriver_strerror(m_r)))
IdxTy fast_mark2( IdxTy uuid=0, const IdxTy flags=0)
{
const bool  set_seen=Bitt(flags,0);
const bool  not_uid=Bitt(flags,1);

 mailstorage * storage= mailstorage_new(0);
IdxTy _r=0;
      MM_ERR(MMPR3( m_host.c_str(), m_port, m_path.c_str())<<MMPR2( m_user.c_str(), m_pass.c_str()))
const bool fck=!true;
 if (fck) {
      _r = imap_mailstorage_init_sasl(storage, m_host.c_str(), m_port, m_path.c_str(), m_conn,
          "xoauth2", NULL, NULL, NULL, NULL, m_user.c_str(), m_pass.c_str(), NULL, false, NULL);
    } else {
      _r = imap_mailstorage_init(storage, m_host.c_str(), m_port, m_path.c_str(), m_conn,
          IMAP_AUTH_TYPE_PLAIN, m_user.c_str(), m_pass.c_str(), false,NULL);
    }
MM_ERR(MMPR2(_r,maildriver_strerror(_r)))
// this should throw but return anyway in c
_r=mailstorage_connect(storage);
MM_ERR(MMPR2(_r,maildriver_strerror(_r)))

mailfolder * mf;
 mf = mailfolder_new(storage, m_path.c_str(), NULL);
mailfolder_connect(mf); 
//m_r = mailsession_get_messages_list(m_f->fld_session, &res);
uint32_t nmsg,nrecent, nus;
_r=mailsession_status_folder(mf->fld_session, m_path.c_str(), &nmsg, &nrecent, &nus);
MM_ERR(MMPR2(_r,maildriver_strerror(_r)))
MM_ERR(MMPR4(_r, nmsg,nrecent, nus))
GetSessEnvelopes(mf->fld_session);

if (set_seen) {
    struct mailimap_set *set;
    set = mailimap_set_new_single(uuid);
MM_ERR("ADCASDCA")
//void SessMarkAsRead(mailsession * sess,  mailimap_set * set, const bool index_set, const IdxTy mflags,  const IdxTy flags)
	SessMarkAsRead(mf->fld_session,  set, not_uid, 0, 0);
MM_ERR("RETWERTEWT")
    mailimap_set_free(set);
}

mailfolder_disconnect(mf); 
mailfolder_free(mf);
mailstorage_disconnect(storage);
mailstorage_free(storage);
return 0;
}




IdxTy fast_mark( IdxTy uuid)
{

//mailstream * mailstream_new(mailstream_low * low, size_t buffer_size);
mailstream * s=0; // mailstream_new(low,  buffer_size);
//mailstream * s=0;
mailsession_driver * sess_driver=imap_session_driver;

mailsession * mysess= mailsession_new( sess_driver);
int rc= mailsession_connect_stream(mysess,  s); 
MM_ERR(MMPR(rc))
//int mailsession_starttls(mailsession * session);
rc= mailsession_starttls(mysess);
MM_ERR(MMPR(rc))
//int mailsession_login(mailsession * session, const char * userid, const char * password);
rc= mailsession_login(mysess, m_user.c_str(), m_pass.c_str());
MM_ERR(MMPR(rc))
//int mailsession_noop(mailsession * session);
//int mailsession_status_folder(mailsession * session, const char * mb, uint32_t * result_messages, uint32_t * result_recent, uint32_t * result_unseen);

//int mailsession_unseen_number(mailsession * session, const char * mb, uint32_t * result);
//int mailsession_subscribe_folder(mailsession * session, const char * mb);


//int mailsession_unsubscribe_folder(mailsession * session, const char * mb);

//int mailsession_logout(mailsession * session);
rc= mailsession_logout(mysess);
MM_ERR(MMPR(rc))
mailsession_free(mysess);
//int mailstream_close(mailstream * s);
mailstream_close(s);
return 0;
} // fast_mark

private:
// https://stackoverflow.com/questions/16208781/how-to-set-email-imap-seen-flag-through-mailcore-libetpan/16217009
void MsgMarkAsRead(mailmessage * msg, const IdxTy flags)
{

    struct mailimap_set *set;
    set = mailimap_set_new_single(msg->msg_index);
	const int mflags= msg->msg_flags->fl_flags;
    MarkAsRead(set,true,mflags,flags); 
    mailimap_set_free(set);
}

void IdxMarkAsRead(const IdxTy idx, const bool index_set,const IdxTy mflags, const IdxTy flags)
{

    struct mailimap_set *set;
    set = mailimap_set_new_single(idx);
    MarkAsRead(set,index_set,mflags,flags); 
    mailimap_set_free(set);
}

void SessMarkAsRead(mailsession * sess,  mailimap_set * set, const bool index_set, const IdxTy mflags,  const IdxTy flags)
{
    struct mailimap_flag_list*flist;
    struct mailimap_store_att_flags * store_flags;
    int err;
    flist = mailimap_flag_list_new_empty();
	mailimap_flag* nfl;
// afck nfl.fl_data=(void*)0;
//nfl.fl_type=mflags;
MM_ERR(" SHTTTT")
nfl=mailimap_flag_new_seen();
MM_ERR(" xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    //mailimap_flag_list_add(flist,mailimap_flag_new_seen()|mflags);
    mailimap_flag_list_add(flist,nfl);
MM_ERR("wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww ")
    store_flags = mailimap_store_att_flags_new_set_flags(flist);
MM_ERR("iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii ")
if ( index_set)
    err = mailimap_store((mailimap*)sess, set, store_flags);
else    err = mailimap_uid_store((mailimap*)sess, set, store_flags);
MM_ERR("rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr ")


MM_ERR(" fcking shtif cuk double frre fck ")
  //  mailimap_store_att_flags_free(store_flags);
  //  mailimap_flag_list_free(flist);

    if (err != MAILIMAP_NO_ERROR) {
m_r=err;
      MM_ERR("** could notup date d ** - "<< MMPR2(m_n, maildriver_strerror(m_r)))
    

}



}
void MarkAsRead(struct mailimap_set * set, const bool index_set, const IdxTy mflags,  const IdxTy flags)
{ //struct mailmessage *messageStruct = [self messageStruct];
    struct mailimap_flag_list*flist;
    struct mailimap_store_att_flags * store_flags;
    int err;
    flist = mailimap_flag_list_new_empty();
    // TODO: ensure that we're not overwriting original flags of message
//	const int mflags= msg->msg_flags->fl_flags;

/*
struct mailimap_flag {
  int fl_type;
  union {
    char * fl_keyword;   // can be NULL 
    char * fl_extension; // can be NULL 
  } fl_data;
*/


	mailimap_flag nfl;
//nfl.fl_keyword=0;
nfl.fl_data=0;
nfl.fl_type=mflags;
// wtf is sht here 
nfl=*mailimap_flag_new_seen();;
    //mailimap_flag_list_add(flist,mailimap_flag_new_seen()|mflags);
    mailimap_flag_list_add(flist,&nfl);
    store_flags = mailimap_store_att_flags_new_set_flags(flist);
if ( index_set)
    err = mailimap_store((mailimap*)m_f->fld_session, set, store_flags);
else    err = mailimap_uid_store((mailimap*)m_f->fld_session, set, store_flags);

//    mailimap_set_free(set);
    mailimap_store_att_flags_free(store_flags);

    if (err != MAILIMAP_NO_ERROR) {
m_r=err;
      MM_ERR("** could notup date d ** - "<< MMPR2(m_n, maildriver_strerror(m_r)))
    

}

}
//////////////////////////////////////////////////
// flags 0= get_envelopes 
void GetList(const IdxTy flags)
{

//SetProgressCallback();
//	mailmessage * msg=0;
//struct mailmime * mime=0;    
mailmessage_list* res;
m_r = mailsession_get_messages_list(m_f->fld_session, &res);
    if (m_r != MAIL_NO_ERROR) {
      MM_ERR("** message  not found ** - "<< MMPR2(m_n, maildriver_strerror(m_r)))
      return; ;
    }
bool get_envelopes=Bits(flags,1); // false;
if (get_envelopes)
{
// fck 
m_r = mailsession_get_envelopes_list(m_f->fld_session, res);
    if (m_r != MAIL_NO_ERROR) {
      MM_ERR("** envelopes   not found ** - "<< MMPR2(m_n, maildriver_strerror(m_r)))
      return; ;
    }

}
else MM_ERR(" not getting envdelope "<<MMPR(flags))

//MM_ERR(" something ")
carray * p= res[0].msg_tab;
IdxTy cnt= p->len;
m_count=cnt;
IdxTy actual=0;
//while (res[cnt]!=0) ++cnt;
MM_ERR("message count "<< MMPR(cnt))
  for(IdxTy i = 0 ; i < carray_count(p) ; i ++) {
    mailmessage * msg;

    msg = (mailmessage*)carray_get(p, i);
bool index=true;
if (get_envelopes)
{
const bool want_seen=Bits(flags,2); // false;
const bool check_seen=Bits(flags,4); // false;
const int mflags= msg->msg_flags->fl_flags;
const bool seen=(Bits(mflags,MAIL_FLAG_SEEN)) ;
//index=! (Bit(flags,MAIL_FLAG_SEEN)) ;
index=!check_seen ||( seen==want_seen);
MM_ERR(MMPR4(want_seen,check_seen,flags,seen)<<MMPR2(index,mflags))
/*

if (Bit(m_bit_flags,MAIL_FLAG_NEW)) Flag("new");
if (Bit(m_bit_flags,MAIL_FLAG_SEEN)) Flag("read");
if (Bit(m_bit_flags,MAIL_FLAG_FLAGGED)) Flag("flagged");
if (Bit(m_bit_flags,MAIL_FLAG_DELETED)) Flag("deleted");
if (Bit(m_bit_flags,MAIL_FLAG_ANSWERED)) Flag("replied");
if (Bit(m_bit_flags,MAIL_FLAG_FORWARDED)) Flag("forwarded");
*/


}
	//m_idx[i]=folder_index(i,msg->msg_index,msg->msg_uid);
//if (index) 	(*m_idxp)[i]=folder_index(i,msg->msg_index,msg->msg_uid);
if (index){
 	(*m_idxp)[actual]=folder_index(actual,msg->msg_index,msg->msg_uid);
++actual; }
//	MM_ERR(MMPR4(i,msg->msg_index,msg->msg_uid,msg->msg_size))
// or maybe not now with get envelopes FCK 
// apparently flags can be non-zero junk 
// dump(msg,i);
}
m_count=actual;
MM_ERR("message count "<< MMPR2(cnt,actual))
mailmessage_list_free(res);

	m_idx_invalid=!true;

}
void dump(mailmessage * msg, const IdxTy i )
{
	if (msg->msg_flags!=0) MM_ERR(MMPR3(i,msg->msg_flags->fl_flags,msg->msg_flags->fl_extension))
   mailimf_fields * x=msg->msg_fields;
if (x!=0) 
{
clist * pf=x->fld_list; // clist *
clistiter* ii= clist_begin(pf);
clistiter* ie= clist_end(pf);
while (ii!=ie)
{
mailimf_field * clc=(mailimf_field*)clist_content(ii);
const char * pray = (const char *)(clc+sizeof(int)); // lol 
if (pray!=0) MM_ERR(MMPR(pray))
++ii;
}
}
} // dump

StrTy ErrText(const IdxTy r)
{ Ss ss;
 ss<< maildriver_strerror(r);
ss<<" "<<MMPR(r);
return ss.str();
}
bool Die(const IdxTy r, const char * txt )
{
if (r==0)  return false; 
m_worth_trying=false;
Ss ss;
ss<<MMPR3(txt,r,ErrText(r));
MM_ERR(ss.str())
MM_ERR(dump())
//MM_ERR(MMPR3(txt,r,ErrText(r)))
rte(ss.str().c_str());
return true; 
}
bool Die(const IdxTy r, const char * txt, const IdxTy verb )
{
if (r==0)  return false; 
m_worth_trying=false;
Ss ss;
ss<<MMPR3(txt,r,ErrText(r));
if (verbosity(verb)) { 
MM_ERR(ss.str())
MM_ERR(dump()) } 
//MM_ERR(MMPR3(txt,r,ErrText(r)))
rte(ss.str().c_str());
return true; 
}




bool Warn(const IdxTy r, const char * txt )
{
if (r==0)  return false; 
Ss ss;
ss<<MMPR3(txt,r,ErrText(r));
MM_ERR(ss.str())
MM_ERR(dump())
//MM_ERR(MMPR3(txt,r,ErrText(r)))
//rte(ss.str().c_str());
return true; 
}



void InitMsg()
{
m_hdr="";
m_body="";
m_verbatim="";
m_bit_flags=0;
m_string_flags.clear();

}
// logs don't OR lol
bool Bits(const IdxTy flags, const IdxTy  b)
{
const IdxTy v= flags&(b);
return (v!=0); 
}
bool Bitt(const IdxTy flags, const IdxTy  b)
{
const IdxTy v= flags&(1<<b);
return (v!=0); 
}


void FetchStuff(mailmessage * msg, const IdxTy flags)
{
StrTy hdr, body;
InitMsg();
char * res=0; // I think this needs to be free'd 
size_t  len=0;
if (Bits(m_fetch_flags, FETCH_HEADER))
{m_r=mailmessage_fetch_header(msg,&res,&len);
if (Die(m_r, "fetch_header" )) return; 
hdr=StrTy(res);
mailmessage_fetch_result_free(msg,res);
}
if (Bits(m_fetch_flags, FETCH_BODY)){ 
m_r=mailmessage_fetch_body(msg,&res,&len);
if (Die(m_r, "fetch_body" )) return; 
body=StrTy(res);
mailmessage_fetch_result_free(msg,res);
}
if (Bits(m_fetch_flags, FETCH_VERBATIM|FETCH_SPLIT)) 
{
// should have blank line sep 
m_r=mailmessage_fetch(msg,&res,&len);
if (Die(m_r, "fetch_both" )) return; 
if (Bits(m_fetch_flags, FETCH_VERBATIM)) { m_verbatim=StrTy(res); } 
if (Bits(m_fetch_flags, FETCH_SPLIT)) 
{
IdxTy ptr=0;
while (ptr<len)
{
char c=res[ptr];
if (c==(char)0) break;
char c2=res[ptr+1];
bool crlf=(( c=='\r') || (c=='\n'));
bool crlf2=(( c2=='\r') || (c2=='\n'));
if (crlf&&crlf2) 
{
char c3=res[ptr+2];
bool crlf3=(( c3=='\r') || (c3=='\n'));
if ( c==c2) { res[ptr+1]=0; ptr+=2; break; } 
if ( c==c3) { res[ptr+2]=0; 
//char cx=res[ptr+3];
ptr+=3; break; } 

}
++ptr;
}
hdr=StrTy(res);
body=StrTy(res+ptr);
}
mailmessage_fetch_result_free(msg,res);
}


if (Bits(m_fetch_flags, FETCH_FLAGS)) 
{
mail_flags * resf;
m_r=mailmessage_get_flags(msg,&resf);
if (Die(m_r, "gettingflags" )) return; 
m_bit_flags=resf->fl_flags;

if (Bits(m_bit_flags,MAIL_FLAG_NEW)) Flag("new");
if (Bits(m_bit_flags,MAIL_FLAG_SEEN)) Flag("read");
if (Bits(m_bit_flags,MAIL_FLAG_FLAGGED)) Flag("flagged");
if (Bits(m_bit_flags,MAIL_FLAG_DELETED)) Flag("deleted");
if (Bits(m_bit_flags,MAIL_FLAG_ANSWERED)) Flag("replied");
if (Bits(m_bit_flags,MAIL_FLAG_FORWARDED)) Flag("forwarded");
// apparently this is dones when the message is deleted.
//mail_flags_free(resf);

}

/*
  MAIL_FLAG_NEW       = 1 << 0,
  MAIL_FLAG_SEEN      = 1 << 1,
  MAIL_FLAG_FLAGGED   = 1 << 2,
  MAIL_FLAG_DELETED   = 1 << 3,
  MAIL_FLAG_ANSWERED  = 1 << 4,
  MAIL_FLAG_FORWARDED = 1 << 5,
  MAIL_FLAG_CANCELLED = 1 << 6
};

  mail_flags is the value of a flag related to a message.
  
  - flags is the standard flags value

  - extension is a list of unknown flags for libEtPan!

struct mail_flags {
  uint32_t fl_flags;
  clist * fl_extension; // elements are (char *) 
};



*/
//MM_ERR(MMPR2(hdr,body))
m_body=body;
m_hdr=hdr;

}
void Flag(const StrTy f) { m_string_flags[f]=f;} 

void GetAMessage()
{
// must be done eraly to avoid stale data 
InitMsg();
	mailmessage * msg=0;
struct mailmime * mime=0;    
const IdxTy midx=(*m_idxp)[m_n].n();
//m_r = mailsession_get_message(m_f->fld_session, m_n, &msg);
if ( m_f==0) rte(" fails sanity check");
if ( m_f->fld_session==0) rte(" fails second sanity check");
m_r = mailsession_get_message(m_f->fld_session, midx, &msg);
    if (m_r != MAIL_NO_ERROR) {
      MM_ERR("** message %i not found ** - "<< MMPR3(midx,m_n, maildriver_strerror(m_r)))
      return; ;
    }

    m_r = mailmessage_get_bodystructure(msg, &mime);
    if (m_r != MAIL_NO_ERROR) {
      MM_ERR("** message  not found ** - "<<MMPR(m_r)<< MMPR4(midx,m_n, msg->msg_size,maildriver_strerror(m_r)))
if ( msg->msg_data!=0) 	MM_ERR("data"<<(const char *) msg->msg_data);
if ( msg->msg_user_data!=0) //	MM_ERR("data"<<(const char *) msg->msg_data);
	MM_ERR("udata"<<(const char *) msg->msg_user_data);
      //printf("** message %i not found - %s **\n", m_n,
       //   maildriver_strerror(m_r));
      mailmessage_free(msg);
	return;     
}
FetchStuff(msg,0);

//dump( msg, 0 );

//    m_r = etpan_render_mime(stdout, msg, mime);

//dump( msg, 0 );

    mailmessage_free(msg);


}


}; // mjm_libetpan_readmsg




#endif // MJM_LIBETPAN_READMSG_H__ 
