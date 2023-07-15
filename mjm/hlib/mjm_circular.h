
#ifndef MJM_CIRCULAR_H__
#define MJM_CIRCULAR_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pthread.h>
#include "mjm_globals.h"
#include "mjm_io.h"

/* 
Buffer input stream for and maybe concat files too
*/
void * mjm_circular_thread(void *);


class mjm_circular
{
private:
typedef mjm_circular Myt;

public:
typedef mjm_generic_traits Tr;
typedef Tr::ChTy ChTy;
typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;
typedef Tr::StrTy StrTy;
typedef unsigned int FlagTy;
typedef Tr::IdxTy  IdxTy;
typedef Tr::ErTy  ErTy;

typedef std::map<StrTy, StrTy> ConfigTy;
typedef ConfigTy::iterator CoItor;
typedef ConfigTy::const_iterator CocItor;

typedef std::vector<StrTy> FileListTy;
typedef FileListTy::const_iterator FileListItor;

const bool verbose;



typedef IdxTy MacU; // a machine sized unsiged int
typedef IdxTy MachU; // a machine sized unsiged int
typedef char MachX; // some unit that should be efficient on the machine
typedef char MachTy; // some unit that should be efficient on the machine
typedef pthread_t ThreadTy;

private:
// I have to check semantics, these are meant for threading
// but this may just be const vio? 
volatile MacU m_rd_ptr, m_wr_ptr;
const MacU m_sz;
// api for stream means this should be char
 MachX * m_tgt;




// consider templated
// is power of 2
// the only real benefit now is using unconditinal
// mask versus mod or test. Normally increments only
// allow one option(>) and then subtract works. 
bool m_pot;

ThreadTy m_loader;

//StrTy m_fn;
FileListTy m_file_list;


enum { MU_WAIT_4_DATA=10000}; 
private:
pthread_mutex_t m_ptr_lock;
pthread_cond_t m_cond;
void get_ptr_lock()
{
pthread_mutex_lock(&m_ptr_lock);
}
void free_ptr_lock()
{
pthread_mutex_unlock(&m_ptr_lock); 

}
IdxTy m_last_read;
bool m_good, m_eof;
// init these so that presumumption is that loader will
// activate if not already active. 
bool m_good_in, m_eof_in;

public:
mjm_circular(MacU n ): verbose(false), m_rd_ptr(0), m_wr_ptr(0),m_sz(n), 
m_tgt( new MachX[m_sz]),m_pot((m_sz&&(m_sz-1))==0), m_loader(0),
m_ptr_lock(PTHREAD_MUTEX_INITIALIZER),m_cond(PTHREAD_COND_INITIALIZER),
m_last_read(0),
m_good(true),m_eof(false),m_good_in(true), m_eof_in(false)  {}


~mjm_circular()
{
delete [] m_tgt;
m_tgt=NULL;
// find pthread delete thingy
pthread_mutex_destroy(&m_ptr_lock);

//delete  m_loader;
//m_loader=NULL;

}
// this doesn't like const
/*
// const StrTy & get(const StrTy &key) const {  return m_config[key]; }
const StrTy & get(const StrTy &key) const 
{  
CocItor ci=m_config.find(key);
if ( ci!=m_config.end()) return (*ci).second;
return m_default; 
//return m_config[key]; 

}
void put(const StrTy & key, const StrTy &v) { m_config[key]=v; }
*/
void load(const StrTy & fn) { load(fn.c_str()); }
void start(const ChTy * fn)
{
//m_fn=fn;
m_file_list.push_back(StrTy(fn));

int rc=pthread_create(&m_loader, NULL, mjm_circular_thread ,(void*) this);
osx()<<MM_MARK<<" pthread create returns "<<rc<<CRLF;

}

void start( const ChTy * const * fn, IdxTy ii, IdxTy imax )
// void start( const ChTy * * fn, IdxTy ii, IdxTy imax )
//void start( ChTy * * fn, IdxTy ii, IdxTy imax )
{
for (IdxTy i=ii; i<imax; ++i)
{
m_file_list.push_back(StrTy(fn[i]));
}

int rc=pthread_create(&m_loader, NULL, mjm_circular_thread ,(void*) this);
osx()<<MM_MARK<<" pthread create returns "<<rc<<CRLF;

}


IdxTy available()
{
IdxTy d=0; 
get_ptr_lock();
// save locks...
// this stupid thing needs to know that the index things
// are UNSIGNED and hopefully temps never less than zed??? 
if ( m_wr_ptr>=m_rd_ptr) d= m_wr_ptr-m_rd_ptr;
else d=m_sz+m_wr_ptr-m_rd_ptr;

free_ptr_lock();

return d;

} 
void load(const ChTy * fn , const bool is_really_the_end )
{
IsTy * is = mjm_io::get_in(fn);
if ( is==NULL) throw  ErTy("not found");
load(*is,is_really_the_end);
delete is;

} 
// this is called by the loading thread after everything is set up.
void load()
{

osx()<<MM_MARK<<" pthread loader running  "<<" " <<CRLF;

// load(m_fn);

FileListItor fi=m_file_list.begin();
FileListItor fe=m_file_list.end();

while (fi!=fe)
{
const char * fnx=(*fi).c_str();
load(fnx,false);
osx()<<MM_MARK<<" pthread loader new file   "<<fnx<<" " <<CRLF;
osx().flush();
++fi;

}


m_eof_in=true;
get_ptr_lock();
pthread_cond_signal(&m_cond);
free_ptr_lock();


}
void load(IsTy & is, const bool is_really_the_end)
{


//ChTy line[m_max_line];
//ChTy * fields[m_max_fields];
if ( verbose) osx()<<MM_MARK<<" loading cbuffer "<<CRLF;

get_ptr_lock();
IdxTy wr=m_wr_ptr;
IdxTy rd=m_rd_ptr;
free_ptr_lock();
while (!is.eof()&&is.good())
{
get_ptr_lock();
//wr=m_wr_ptr;
rd=m_rd_ptr;
free_ptr_lock();
// this points to the invalid location
// neother sz nor rd should be written to
// but the test is when? 


IdxTy space=m_sz-wr;
if ( rd==0) --space;
if ( rd>wr) space=rd-wr-1;


// this hangs, probably due to inability to fill last niche
// if ( space>0) --space; 
// small writes into large buffer seem to be reliable
// probably is likely the wrap only. 
//if ( space>100) space=100;
if ( space==0)
{
if (verbose) { osx()<<MM_MARK<<" publishing due to no space "
<<rd<<" "<<wr<<CRLF;osx().flush(); }
get_ptr_lock();
//pthread_cond_signal(&m_cond);
// should recalc here once code is good

free_ptr_lock();

}
IdxTy got=0; 
//IdxTy got=(space==0)?0:is.readsome(m_tgt+wr ,space);
//IdxTy ffff=0;
if (space!=0)
{ 

//is.read(m_tgt+wr,space);
//is.readsome(m_tgt+wr,space);
//got=is.gcount();
IdxTy a=is.rdbuf()->in_avail();
if ( a>space) a=space;
if ( a==0 ) a=space-(space>>1);
if (got==0){
 is.read(m_tgt+wr,a);
got=is.gcount();


}
}



// This doesnt FCKNG find an EOF even in s file stream,
// maybe usefulfro files being updated etc. 
// IdxTy got=(space==0)?0:is.readsome(m_tgt+wr ,space);
if ( got==0) { 
//  may publish pointeres here 
if (0)
{
get_ptr_lock();
if ( m_rd_ptr!=m_wr_ptr) pthread_cond_signal(&m_cond); 
free_ptr_lock();
}

if ( verbose)
{
osx()<<MM_MARK<<" fck got is zed available == "<<available()
<<" space is "<<space<<" wr "<<m_wr_ptr<<" "
<<m_rd_ptr<<" is eof="<< is.eof()<<" "<<is.good()<<CRLF;
osx().flush();
}
usleep(MU_WAIT_4_DATA*1); continue;
// usleep(1); continue;



 } 
if ( verbose) 
{
osx()<<MM_MARK<<" loading data wr="<<wr<<" got "<< got
<<" available is now "<<available()<<CRLF; osx().flush(); } 


wr+=got;
//
//IdxTy had_before;
// atomic
// if ( wr>m_sz) throw ErTy( "size too big");
if ( wr==m_sz) wr=0; // -=m_sz;
get_ptr_lock();
m_wr_ptr=wr;
//free_ptr_lock();
//get_ptr_lock();
// this broadcase is only needed if someone else was waiting
// for data. When the other gy blocks, should publish current
// pointers so we can check for abail == 0 
pthread_cond_signal(&m_cond);
free_ptr_lock();
//if ( had_before==0)

if ( verbose) {
osx()<<MM_MARK<<" LOADER got= "<< got
<<" available is now "<<available()<<" rd="<<m_rd_ptr<<" "<<m_wr_ptr<<CRLF;
osx().flush();
} // verbose

} // while 

if ( verbose) 
{
osx()<<MM_MARK<<"XXXXXXXX dun loading file  "<<CRLF; osx().flush(); }

if ( is_really_the_end)
{
m_eof_in=true;

get_ptr_lock();
pthread_cond_signal(&m_cond);

free_ptr_lock();
} // real end

}







bool good() {

if ( verbose) {osx()<<MM_MARK<<" good "<<m_good<<CRLF; }

return m_good; }
bool eof() {if ( available()==0) m_eof|=m_eof_in;

if ( verbose) {osx()<<MM_MARK<<" eof "<<m_eof<<CRLF; }
 return m_eof;  }
IdxTy gcount() {return m_last_read; }





IdxTy getlinex( MachTy* const bstar, const MachU maxlen)
{
IdxTy len_read=0;
MachTy * b=bstar;
if ( verbose) 
{
osx()<<MM_MARK<<" get linemax "<<maxlen<<"  with avauk="<<available()<<CRLF; 
osx().flush(); }


if ( maxlen==0) return 0;
// make this len of string, tack on zed 
const MachTy * bend=b+maxlen-1;
// for strict equal nothing is there 
// this pushed out possibly invalid valued

get_ptr_lock();
IdxTy wr=m_wr_ptr;
IdxTy rd=m_rd_ptr;
free_ptr_lock();

while ( b<bend)
{
// in this state, we should reload 
if ( wr==rd)
{

get_ptr_lock();
wr=m_wr_ptr;
//m_rd_ptr=rd;
free_ptr_lock();
if ( rd==wr)
{
if (eof()) break;
// need notify or wakeup etc.j
if ( verbose) { osx()<<MM_MARK<<"reader waiting on data  "<<CRLF; osx().flush();  } 
get_ptr_lock();
// above is a faster but hazard test, 
if ( m_wr_ptr==m_rd_ptr) if ( !m_eof_in ) pthread_cond_wait(&m_cond,&m_ptr_lock);


free_ptr_lock();
if ( verbose) { osx()<<MM_MARK<<"signalled ok "<<CRLF;osx().flush();  } 
if ( available()==0) 
{
if ( verbose) { osx()<<MM_MARK<<"read dies on no availab le wake"<<CRLF; osx().flush(); } 
break; // woke up for nothing.
} //break
}
}


MachTy c= m_tgt[rd++];
if ( rd==m_sz) rd=0;
// turning this off seems to make things worse
//if (0)
{
get_ptr_lock();
m_rd_ptr=rd;
free_ptr_lock();
}

++len_read;
if ( c=='\n')
{
//ZZif (0) 
if ( verbose) osx()<<MM_MARK<<" n found rd="<<m_rd_ptr<<" w="<<m_wr_ptr<<CRLF;
 break; // read the new line and exit. 
}
// dos files come a mess. 
if ( c=='\r')
{
osx()<<MM_MARK<<"Found superfluous cr in data do not use DOS format "<<CRLF;
osx().flush();
 continue;
}
*(b++)=c;
} // while
// we SHOULD also remove all following '\r' etc 
*b=0;
get_ptr_lock();

m_rd_ptr=rd;


m_last_read=len_read;
free_ptr_lock(); 
if ( verbose)
{
osx()<<MM_MARK<<" DUN READING    rd="<<m_rd_ptr<<" w="<<m_wr_ptr<<
" len="<<len_read<< CRLF; osx().flush();
}

// the zero goes at b 
return len_read;


}



IdxTy getline( MachTy* const bstar, const MachU maxlen)
{
IdxTy len_read=0;
MachTy * b=bstar;

if ( maxlen==0) return 0;
get_ptr_lock();
IdxTy  rd=m_rd_ptr; 
IdxTy  wr=m_wr_ptr; 
//volatile IdxTy & rd=m_rd_ptr; 
free_ptr_lock();
// make this len of string, tack on zed 
const MachTy * bend=b+maxlen-1;
// for strict equal nothing is there 
// this pushed out possibly invalid valued


while ( b<bend)
{
// in this state, we should reload 
if ( wr==rd)
{ 
// need notify or wakeup etc.j
get_ptr_lock();
// above is a faster but hazard test, 
wr=m_wr_ptr;
if ( wr==rd) if ( !m_eof_in )
{
m_rd_ptr=rd;
 pthread_cond_wait(&m_cond,&m_ptr_lock);

}
free_ptr_lock();
if ( available()==0) 
{
if ( verbose) { osx()<<MM_MARK<<"read dies on no availab le wake"<<CRLF; osx().flush(); } 
break; // woke up for nothing.
} //break
} // fast check


// get_ptr_lock();
MachTy c= m_tgt[rd++];

 if (rd==m_sz)
{
get_ptr_lock();
 rd=0;
m_rd_ptr=rd; 
free_ptr_lock();
}

// turning this off seems to make things worse
//if (0)
// free_ptr_lock();

// ++len_read;
if ( c=='\n')
{
len_read=b-bstar+1; // len_read;
//ZZif (0) 
if ( verbose) osx()<<MM_MARK<<" n found rd="<<m_rd_ptr<<" w="<<m_wr_ptr<<CRLF;
 break; // read the new line and exit. 
}
// dos files come a mess. 
if ( c=='\r')
{
osx()<<MM_MARK<<"Found superfluous cr in data do not use DOS format "<<CRLF;
osx().flush();
 continue;
}
*(b++)=c;
} // while
// we SHOULD also remove all following '\r' etc 
*b=0;

get_ptr_lock();
m_rd_ptr=rd; 
free_ptr_lock();

//m_last_read=b-bstar+1; // len_read;
m_last_read= len_read;

// the zero goes at b 
return len_read;


}




/*

void dump(OsTy & os)
{
CoItor ii= m_config.begin();
CoItor ie= m_config.end();
while(ii!=ie)
{

os<<(*ii).first<<" "<<(*ii).second<<CRLF;
++ii; 

} // while




}
*/






};  // mjm_io

void * mjm_circular_thread(void * x )
{

mjm_circular * t= (mjm_circular*)x;
t->load();
return 0;

}

#endif

