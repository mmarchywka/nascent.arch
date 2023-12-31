#ifndef MJM_XDVIK_INTERFACE_H__
#define MJM_XDVIK_INTERFACE_H__

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
#include <signal.h>


//#ifdef HAVE_XDVIK


#warning including xdvik
extern "C" {


#include <X11/extensions/XTest.h>
#include "xdvi-config.h"

#include <locale.h>

#include <X11/Intrinsic.h>

#include "kpathsea/proginit.h"
#include "kpathsea/expand.h"

#include "xdvi.h"
#include "util.h"
#include "x_util.h"
#include "sfSelFile.h"
#include "my-snprintf.h"
#include "dvi-init.h"
#include "filehist.h"
#include "mag.h"
#include "message-window.h"


//#include "xdvi_code.h"
#include "xdvi_code_main.h"

extern int mjm_dy_font();
extern int mjm_fix_y();
extern void mjm_bbox(int x, int y, int w, int h);
extern void mjm_set_string_block (const char * p, int * _x, int * _y, int * _w, int * _h  );

 //texk/xdvik/mjm_xdvi_stuff.h
#define interface_only
#include "xdvik/mjm_xdvi_stuff.h"


};

//#endif


#define MM_SET_GET(x,y) const IdxTy & x() const { return y; } const IdxTy & x(const IdxTy v ) {y=v;  return y; } 



// Tue Jul  9 20:06:35 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_xdvik_interface   
// g++ -std=gnu++11 -DTEST_MJM_XDVIK_INTERFACE -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_xdvik_interface.h  -lpthread -lreadline

template <class Tr>
class mjm_xdvik_interface 
{
 typedef mjm_xdvik_interface Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
enum { USR_EVENT_MU=0, PAGE_LOAD=1 } ; 
// should use "bad" but due to preamble zero is not called. 
mjm_xdvik_interface():m_prior_vis(1),m_blockes_since_bop(0),m_inter(0),m_last_pc(0), m_max_pc(0),m_active(0), m_flags(0),
m_mutex_vector(1)
 {

mjm_block_callback=&Myt::SpecialsHandler;
//extern void (*mjm_select_callback)(int x, int y, int w, int h);
mjm_select_callback=&Myt::XSelectHandler;

mjm_block_callback_obj=(this);

}
~mjm_xdvik_interface() {}
void signal() { Signal(); } 
static void handle_usr2(int signo) { HandleUsr2(signo);}
bool running() {return ( mjm_made_it_to_mjm_do_pages!=0);  }
bool failed() {return ( mjm_failed_to_init!=0);  }
void launch_page(const StrTy & fn)
{ LaunchPage(fn); } 

// TODO FIXME threading issue here, should set event flag or something for GUI to do . 
void redraw() {
globals.ev.flags |= EV_RELOAD;
MM_ERR(" sending fake event")
// ::redraw_page();
if (!running() || failed()) 
{
MM_ERR(" ignoring XEvent "<<MMPR2(running(),failed()))
return;

}
//XTestFakeMotionEvent(DISP,0,0,0,0);
XClientMessageEvent ce;
XEvent xe;
xe.xclient=ce;
XSendEvent(DISP,0,0,0,&xe);


 } 
bool wait_eop(const IdxTy n=0) { bool b= WaitEOP(n);MM_ERR(" wait return ") return b;   } 

IdxTy show(const StrTy & n,const IdxTy flags) { return Show(n,flags); }
IdxTy hide(const StrTy & n,const IdxTy flags) { return Hide(n,flags); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:

void Set(const IdxTy mask, const IdxTy mu) { EnterSerial(mu);  m_flags|=mask; ExitSerial(mu); }
void Set(const IdxTy mask) {   m_flags|=mask; }
void Reset(const IdxTy mask, const IdxTy mu) { EnterSerial(mu);  m_flags&=~mask; ExitSerial(mu); }
void Reset(const IdxTy mask) {   m_flags&=~mask;  }
bool Flags(const IdxTy mask, const IdxTy mu)
{
bool b=false; 
EnterSerial(mu);
 b=((m_flags&mask)!=0);
ExitSerial(mu);
return b;
}





bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const IdxTy flags=0) {Ss ss;  
ss<<m_vb_map.dump(flags); 

return ss.str(); }

class ViewingBlock
{
typedef ViewingBlock Myt;
public:
ViewingBlock(const StrTy & n, const IdxTy order )
	:m_order(order), m_name(n),m_x(bad()), m_y(bad()), m_xe(bad()), m_ye(bad()),m_visibility(1),m_start(0),m_end(~0) {} 
ViewingBlock(const StrTy & n ):m_order(0), m_name(n),m_x(bad()), m_y(bad()), m_xe(bad()), m_ye(bad()),m_visibility(1),m_start(0),m_end(~0) {} 
ViewingBlock(): m_order(0),m_x(bad()), m_y(bad()), m_xe(bad()), m_ye(bad()),m_visibility(1),m_start(0),m_end(~0) {} 
//ViewingBlock(const Myt & n ): m_name(n.m_name),m_x(n.m_x), m_y(n.m_y), m_w(n.m_w), m_h(n.m_h),m_visibility(n.m_visibility) {} 

static IdxTy bad() { return ~0; } 
const StrTy & name() const { return m_name; } 
const StrTy & name(const StrTy & n )  {m_name=n;  return m_name; } 
const IdxTy & visibility() const { return m_visibility; } 
const IdxTy & visibility(const IdxTy v ) {m_visibility=v;  return m_visibility; } 
const IdxTy & flip(const IdxTy v ) {m_visibility^=v;  return m_visibility; } 

MM_SET_GET(start,m_start)
MM_SET_GET(end,m_end)
const void update_start(const IdxTy x, const IdxTy y) { m_x=x; m_y=y; } 
const void update_start_y( const IdxTy y) {  m_y=y; } 
const void update_end(const IdxTy x, const IdxTy y) { m_xe=x; m_ye=y; } 
const bool is_in(const IdxTy x, const IdxTy y) { return  (m_x<=x)&&(m_xe>=x)&&(m_y<=y)&&(m_ye>=y); } 
// return true if any part of this is in the rectangle in parameter 
const bool is_in(const IdxTy x, const IdxTy y,const IdxTy w, const IdxTy h) 
{ 
const bool yy=includes(y,y+h,m_y,m_ye); // is_in(m_y,m_ye,y)||is_in(m_y,m_ye,y+h)||includes(m_y,m_ye,y,y+h);
if (!yy) return false; 
// if the text is visible the entire line is available 
if ((m_visibility==1) ) return !false; 
const bool xx=includes(x,x+w,m_x,m_xe); // 
//const bool xx=is_in(m_x,m_xe,x)||is_in(m_x,m_xe,x+w)||includes(m_x,m_xe,x,x+w);
return xx;

} 
// if x between x0 and x1
const bool is_in(const IdxTy x0, const IdxTy x1,const IdxTy x) { return  ( x0<=x)&&(x1>=x); }  
const bool includes(const IdxTy x0, const IdxTy x1,const IdxTy xi, const IdxTy xf) 
{ return  is_in(x0,x1,xi)||is_in(x0,x1,xf)||is_in(xi,xf,x0)||is_in(xi,xf,x1); }  
const bool includesxx(const IdxTy x0, const IdxTy x1,const IdxTy xi, const IdxTy xf) { return  ( x0>=xi)&&(x1<=xf); }  

const bool in(const IdxTy p ) const { 
// for now INCLUDE missing end, see ctor defaults
bool bin= (p>=m_start)&&(p<=m_end);
//MM_ERR(MMPR4(m_name,p,m_start,m_end)<<MMPR(bin))
return bin;
 } 

void outline() 
{ 
int dy=(mjm_visibility==1)?0:(mjm_dy_font()>>1);
	 mjm_bbox( m_x, m_y+(0), m_xe-m_x,m_ye+dy-m_y); 

} 

StrTy dump() const {
Ss ss;
ss<<MMPR4(m_name,m_x,m_y,m_visibility);
return ss.str();
}
StrTy geo() const {
Ss ss;
ss<<MMPR2(m_name,m_visibility);
ss<<MMPR4(m_x,m_xe,m_y,m_ye);
return ss.str();
} // geo 
IdxTy y() const { return m_y;}
private:
StrTy m_name;
IdxTy m_order;
int m_x,m_y,m_xe,m_ye;
IdxTy m_visibility;
IdxTy m_start,m_end;

}; // ViewingBlock

class ViewingBlockMap 
{
typedef ViewingBlockMap Myt;
typedef ViewingBlock Tgt;
typedef typename std::map<StrTy, Tgt > Map;
typedef typename std::map<IdxTy, Tgt* >   CoverMap;
public:
typedef typename Map::iterator itor;

bool have(const StrTy & n) const { return m_map.find(n)!=m_map.end(); } 
// this should only be called in response to reading DVI file so that
// if not seen it is added in the order it is found. 
Tgt & get(const StrTy & n, const IdxTy pos, const IdxTy flags ) { return GetOrMake(n,pos,flags); } 
Tgt & get(const StrTy & n, const IdxTy pos ) { return GetOrMake(n,pos); } 
Tgt & get(const StrTy & n ) { return GetOrMake(n); } 
Tgt & get(const IdxTy  n ) { return GetOrMake(m_ord[n]); } 
Tgt * lb(const IdxTy i ) { return Lb(i); } 
IdxTy size() const { return m_ord.size(); } 
StrTy name(const IdxTy i) { return m_ord[i];}
// very slow way to do prior
Tgt * prior(const IdxTy i ) { if (i==0) return 0; return &m_map[m_ord[i-1]]; } 
Tgt * next(const IdxTy i ) { if (i>=size()) return 0; return &m_map[m_ord[i+1]]; } 
itor begin() { return m_map.begin();}
itor end() { return m_map.end();}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
// TODO FIXME the threading here is a mess, GUI and CL conflict

Tgt & GetOrMake(const StrTy & n ) { 
itor ii= m_map.find(n);
if (ii!=m_map.end()) return (*ii).second;
// should not occur unless pages views ooo 
m_ord.push_back(n); 
m_map[n]=Tgt(n,m_map.size());
return m_map[n]; 
} 


Tgt & GetOrMake(const StrTy & n, const IdxTy pos  ) { 
itor ii= m_map.find(n);
if (ii!=m_map.end()) return (*ii).second;
m_ord.push_back(n); 
m_map[n]=Tgt(n,m_map.size());
return m_map[n]; 
} 


Tgt & GetOrMake(const StrTy & n, const IdxTy pos, const IdxTy flags  ) { 
itor ii= m_map.find(n);
if (ii==m_map.end()){
//  return (*ii).second;
m_ord.push_back(n); 
m_map[n]=Tgt(n,m_map.size());
}
Tgt& x=m_map[n];
const bool is_start=Bit(flags,0);
const bool is_end=Bit(flags,1);
if (is_start) { m_cover[~pos]=&x;  x.start(pos); } 
if (is_end) { x.end(pos); } 

return x; 
} 

Tgt * Lb(const IdxTy i ) { 

auto ii=m_cover.lower_bound(~i); 
if (ii==m_cover.end()) return 0; 
// it could be that the "end" of this block has not been seen yet
// in which case the block should be visible for now and set active to zero
Tgt* p=(*ii).second;
MM_ERR(MMPR2(i,p->name()))
if (p->in(i)) return p; 
return 0; 
} 



StrTy Dump(const IdxTy flags=0) { 
Ss ss; 
MM_LOOP(ii,m_map)  ss<<(*ii).second.dump()<<CRLF;
return ss.str();
}
bool Bit(const IdxTy x, const IdxTy b) { return ((x&(1<<b))!=0); }

Map m_map;
typedef std::vector<StrTy> Ordering;
Ordering m_ord;
CoverMap m_cover;

}; // ViewingBlockMap



class ThreadParam
{
public: 
StrTy m_fn;
sigset_t ss;
};
static void * ThreadFunction(void * _p)
{
const ThreadParam * p = (const ThreadParam*)_p;
const StrTy fn=p->m_fn;
const char * args[2];
args[0]="./a.out";
args[1]=fn.c_str();
const IdxTy argc=2;
 (void)sigdelset(&(p->ss), SIGUSR2);

sigprocmask(SIG_SETMASK,&(p->ss),0);
delete p;
main_xdvi(argc,args);
return 0;
}
void LaunchPage(const StrTy & fn )
{
Reset(PAGE_LOAD, USR_EVENT_MU);
typedef mjm_thread_util<Tr> Threads;
ThreadParam*  p= new ThreadParam();
p->m_fn=fn;
sigset_t fs;
sigfillset(&fs);
 (void)sigaddset(&fs, SIGUSR2);
// so right now all sigs disabled, may want to check that lol 
sigprocmask(SIG_BLOCK,&fs,&(p->ss));
 (void)sigdelset(&(p->ss), SIGUSR2);
Threads::fire_and_forget( 1, &Myt::ThreadFunction, p);

}

void Signal() { 
//::signal(SIGUSR2, &Myt::handle_usr2); 
//kill(getpid(),SIGUSR2); 
kill(getpid(),SIGUSR1); 
}

static void HandleUsr2(int signo) {

// HandleUsr2(signo);
MM_ERR(" handling usrsig2")

}

// this needs to be fast however... 
class special_parse { 
typedef std::vector<StrTy> Psp;
enum { MAX_LEN=(1<<10), MAX=MAX_LEN-3 } ;
public:
special_parse(const char * cp): m_group(0),m_cmd(0),m_block(0)  {m_data[0]=0;  Parse(cp); }
const bool is_us() const { return strcmp(group(),"::mjm")==0; } 
const bool is(const char * _cmd ) const { return strcmp(cmd(),_cmd)==0; } 
const char * group() { return m_data+m_group; } 
const char * cmd() { return m_data+m_cmd; } 
const char * block() { return m_data+m_block; } 
const StrTy dump() const { Ss ss; ss<<MMPR3(group(),cmd(),block()); return ss.str(); }
private:
static void mjm_word(const char * cp, int start, int * end)
{
int p=start;
while ( ( cp[p]!=' ')&&(cp[p]!=0)) ++p;
*end=p;
}
static void Parse( Psp& v, const char * cp)
{
int start=0;
int end=0;
int * pend=&end;
mjm_word(cp,start,pend);
if (cp[end]==0) return;
start=end+1;

}
void Parse(const char * cp)
{
IdxTy ptr=0;
IdxTy pc=0;
// leading spaces
while (cp[ptr]==' ') ++ptr;
m_group=pc;
while (cp[ptr]!=' ')
{
m_data[pc]=cp[ptr]; 
if (cp[ptr]==0) break;
 ++ptr;
++pc;
if ( pc==MAX) {m_data[pc]=0;  MM_ERR(" too long "<<cp ); return ; } 
}
m_data[pc]=0; ++pc;

while (cp[ptr]==' ') ++ptr;
m_cmd=pc;
while (cp[ptr]!=' ')
{
m_data[pc]=cp[ptr]; 
if (cp[ptr]==0) break;
 ++ptr;
++pc;
if ( pc==MAX) {m_data[pc]=0;  MM_ERR(" too long "<<cp ); return ; } 
}
m_data[pc]=0; ++pc;
while (cp[ptr]==' ') ++ptr;
m_block=pc;
while (true)
{
m_data[pc]=cp[ptr]; 
if (cp[ptr]==0) break;
 ++ptr;
++pc;
if ( pc==MAX) {m_data[pc]=0;  MM_ERR(" too long "<<cp ); return ; } 

}

} // Parse


char m_data[MAX_LEN];
IdxTy m_group, m_cmd,m_block;

}; // special_parse

bool WaitEOP( IdxTy n )
{
MM_ERR(" WaitEOP")
while ( true )
{
if (Flags(PAGE_LOAD, USR_EVENT_MU)) return true; 
--n;  // note that zero is effecitvely forever 
if ( n==0) return false; 
sleep(1);

MM_ERR(" waiting to load "<<MMPR(m_flags))
}
}

static void* XSelectHandler(int x, int y, int w, int h,int b)
{

MM_ERR("XSelectHandler "<<MMPR4(x,y,w,h)<<MMPR(b))


return ((Myt*)mjm_block_callback_obj)->ThisXSelect(x,y,w,h,b);

}
 void*  ThisXSelect(int x, int y, int w, int h,int b)
{
IdxTy hits=0;
EnterSerial(USR_EVENT_MU);
bool flip_last=(2==b);
ViewingBlock * la=0;
const IdxTy sz=m_vb_map.size();
for(IdxTy i=0; i<sz; ++i)
//MM_LOOP(ii,m_vb_map)
{
//ViewingBlock & vb=m_vb_map.get(bname,m_last_pc,1);
ViewingBlock & vb=m_vb_map.get((i)) ; // (*ii).second;
//MM_ERR(" testing "<<MMPR4(x,y,w,h)<<MMPR(vb.geo()))
if (vb.is_in(x,y,w,h) )
{  MM_ERR(" selection includes "<<MMPR(vb.name())<<MMPR4(x,y,w,h)<<MMPR(vb.geo())) 
vb.flip(1);
if (flip_last) { la= m_vb_map.next(i);  if (la!=0) la->flip(1); } 
++hits;
}

la=&vb;
} // ii 

ExitSerial(USR_EVENT_MU);
if (hits!=0)  redraw();
return NULL;
}

static void SpecialsHandler(void * obj, char * cp, int x, int y)
{
//if (true ) return; 
//MM_ERR(" SpecialHandler "<<MMPR3(cp,x,y))
Myt * us=(Myt*)obj;
if (us==0){
MM_ERR(" us iz null, making a memory leak ")
 us = new Myt(); 
mjm_block_callback_obj=(us);
}
// need to figure out how to get the longjump target and rethrow lol. 
if (cp==~0) { MM_ERR(" danger unwinding with longjmp ")  us->ExitSerial(USR_EVENT_MU);  return; }

us->EnterSerial(USR_EVENT_MU);
 us-> OurSpecialsHandler(us,  cp,  x, y);
us->ExitSerial(USR_EVENT_MU);
}
#define MYBOP 0x08B 
#define MYEOP 0x08C 


//int _x,_y,_w,_h;
//mjm_set_string_block (NULL, & _x, & _y, & _w, & _h  );
//vb.update_start(_x,_y-_h);
//m_y_active

void UpdateYActive()
{
int _x,_y,_w,_h;
mjm_set_string_block (NULL, & _x, & _y, & _w, & _h  );
//vb.update_start(_x,_y-_h);
m_y_active=_y;


}

void UpdateActive()
{
m_active=m_vb_map.lb(m_last_pc);
MM_ERR(" update active "<<mactive(m_active))
if (m_active!=0) NewVisibility(m_active->visibility());
else ToVisible();


}
StrTy mactive(ViewingBlock * a)
{
if (a==0) return "NULL";
return a->name();
}

// default page breaks may leave dangling blocks. 
// these may either straddle page breaks and the end has not been found yet
// or the end may already be known. In the first case visibility is 
// assumed but in the second may be invisible. Normally wait until
// end to get right font but want one on this page in whatever
// font is being used. May have tree structure but for now
// flat although also overlapping is possible. Just check all for now
//   

void CheckOpenBlocks()
{
if (m_active==0) return; 
ViewingBlock & vb=*m_active;
const bool prior_vis=(vb.visibility());
const bool set_last=true&&(!prior_vis);
if (set_last)
{
mjm_y_onlys(vb.y()+mjm_dy_font());
int _x,_y,_w,_h;
char dumm=0;
mjm_set_string_block (&dumm, & _x, & _y, & _w, & _h  );
mjm_set_string_block (vb.name().c_str(), & _x, & _y, & _w, & _h  );
vb.update_start(_x,_y-_h);
vb.update_end(_x+_w,_y);
vb.outline();
MM_ERR(" updating block on this page "<<MMPR(vb.geo()))
}
}
int GetY()
{
// I'm guessing the font is not set to null before the code starts lol. 
return mjm_y_only();
}
 void PCAdvance( char * cp, int x, int y)
{
 //UpdateYActive();
// this was working better when stuck wtf
//m_last_pc=87; // y; 
m_last_pc= y; 
if (y>m_max_pc) m_max_pc=y;
//MM_ERR(" setting PC "<<MMPR4(x,(x==MYBOP),m_last_pc,m_max_pc))
// setting visibility to 1 messes up if a block crosses
// page boundaries. But if these is none, it should start
// visible 
//if ( m_last_pc>100) MM_MSG(" uck "<<MMPR3(std::hex,m_last_pc,GetY()));

if (x==MYBOP) {m_blockes_since_bop=0; mjm_visibility=1;  m_active=0; MM_ERR(" BOP found "<<mactive(m_active))   Reset(PAGE_LOAD); }
// already locked now 
//if (x==MYEOP) Set(PAGE_LOAD, USR_EVENT_MU);
else if (x==MYEOP){CheckOpenBlocks();  MM_ERR(" EOP found "<<mactive(m_active))  Set(PAGE_LOAD); } 
// if not already in one, wait to encounter block start or stop
// but this does not work at top of page.  
else if (m_active==0){   return; } 
if (m_active!=0) if (m_active->in(m_last_pc)){ UpdateYActive();  return;  } 
// one was active now we are oor. This may be due to a page change
// but need to look for one. 
// 
UpdateActive();
//if (m_active!=0 ) 
UpdateYActive();
}
void Change(ViewingBlock * old)
{
if (m_active!=old) 
{
//MM_ERR(" active changed from "<<MMPR3(mactive(old),mactive(m_active),m_last_pc));

}

}
 void OurSpecialsHandler(void * obj, char * cp, int x, int y)
{

ViewingBlock* act=m_active; 
if (cp==0) { 
// if we had a non-vis as the last thing, issue a crlf
// this puts in a spur at top of page but ok for now 
// also this can be a follow on special, 
if (m_inter==0) {
// kluge but ok for now although check dvi code.. 
if ((x<239)||(x>242)) mjm_xdvi_xr(1);
}
++m_inter; PCAdvance(cp,x,y);Change(act);  return; } 
special_parse sp(cp);
if (y<0)
{
MM_ERR(" fixing bad y "<<MMPR4(cp,x,y,m_last_pc))
mjm_fix_y();
int _x,_y,_w,_h;
char dumm=0;
 mjm_set_string_block (&dumm, & _x, & _y, & _w, & _h  );
}

//{ MM_MSG("  special "<<MMPR3(m_last_pc,cp,sp.dump()))  } 
if (!sp.is_us()) { MM_MSG(" wrong special "<<MMPR2(cp,sp.dump())) return ; } 
if (sp.is("blockstart")) {
++m_blockes_since_bop;
StrTy bname=sp.block();
//MM_MSG(" block start "<<MMPR4(sp.dump(),cp,y,mjm_dy_font()))
ViewingBlock & vb=m_vb_map.get(bname,m_last_pc,1);
vb.update_start(x,y-mjm_dy_font());
NewVisibility(vb.visibility());
m_active=&vb;
//MM_MSG(MMPR4(bname, mjm_visibility,x,y))
// if we are already invisible, just continue on 
// otherwise insert a crlf and start on a new line. 
// push and pop are problems though 
const bool interspersed=(m_inter!=1);
if ((mjm_visibility==0))
{

if (interspersed||!true) mjm_xdvi_xr(0);
// this one causes errors if omited lol. 
//	MM_ERR(" trying to typset "<<MMPR(bname))
//	 mjm_set_string(bname.c_str());
// won't work since we don't have font yet 
const bool set_first=false;
if (set_first)
{
int _x,_y,_w,_h;
char dumm=0;
// interspersed free text?
//if (m_blockes_since_bop==1) 
mjm_set_string_block (&dumm, & _x, & _y, & _w, & _h  );
mjm_set_string_block (bname.c_str(), & _x, & _y, & _w, & _h  );
vb.update_start(_x,_y-_h);
vb.update_end(_x+_w,_y);
} // set_first


} else
{
// each block muyst start on left anywey 
mjm_xdvi_xr(0);
int _x,_y,_w,_h;
mjm_set_string_block (NULL, & _x, & _y, & _w, & _h  );
vb.update_start(_x,_y-_h);


} 
Change(act);
//MM_ERR(" done with block start "<<MMPR2(sp.dump(),cp))
return ;
} // blockstart

if (sp.is("blockend")) {
 mjm_visibility=1;
//MM_MSG(" block end "<<MMPR2(sp.dump(),cp))
StrTy bname=sp.block();
ViewingBlock & vb=m_vb_map.get(bname,m_last_pc,2);
// this seems to work now 
//MM_ERR(" end of block "<<MMPR4(vb.name(),m_last_pc,vb.start(),vb.end()))
m_prior_vis=(vb.visibility());
// if this was visible the block extends the entire text length
static const int xmax=10000; 
// if we are setting at the end and not visible do that now
// with latest font 
const bool set_last=true&&(!m_prior_vis);
if (set_last)
{
int _x,_y,_w,_h;
char dumm=0;
// interspersed free text?
//if (m_blockes_since_bop==1) 
mjm_set_string_block (&dumm, & _x, & _y, & _w, & _h  );
mjm_set_string_block (bname.c_str(), & _x, & _y, & _w, & _h  );
vb.update_start(_x,_y-_h);
vb.update_end(_x+_w,_y);
} // set_last

// although ideally the is_in checks would just drop the y comp if visible 
//if (m_prior_vis) vb.update_end(m_prior_vis?xmax:x,y);
// this does not help due to pop before calling us. 
// if (m_prior_vis) vb.update_end(xmax,y);
 if (m_prior_vis) vb.update_end(xmax,m_y_active);
//if( !m_prior_vis) vb.update_start_y(y-mjm_dy_font());
if (!m_prior_vis) vb.outline();
m_active=0;
//MM_ERR(" end of block "<<MMPR4(vb.geo(),m_last_pc,vb.start(),vb.end()))

Change(act);
m_inter=0;
return; 
} // blockend 
if (sp.is("visible")) { ToVisible(); 

Change(act);

return; }
Change(act);

MM_MSG(" bad special block  "<<MMPR2(sp.dump(),cp))
} // OurSpecialHandler 

#if 0 
void OldOurSpecialsHandler(void * obj, char * cp, int x, int y)
{

// the dvi interpretter calls us for each instruction executed but
// there is no assurance that pages are viewed in order 
if (cp==0) { PCAdvance(cp,x,y); return; } 
i// otherwise this is a special with a string to parse 
// could use the command parser but write low level here for speed
int start=0;
int end=0;
int * pend=&end;
mjm_word(cp,start,pend);
if (memcmp(cp+start, "::mjm", end-start) != 0) {
//printf(" got wrong special %s end %d start %d \n ",cp,end,start);
MM_ERR(" wrong special "<<MMPR3(cp,start,end))
return; 
} // ::mjm
if (start==end) return;
if (cp[end]==0) return;
start=end+1;
mjm_word(cp,start,pend);
// this now points to the type of word 
// but can not yet extract the names etc 

if (memcmp(cp+start, "blockstart", end-start) == 0) {
//printf(" block start  %s end %d start %d \n ",cp,end,start);
// don't allow null names 
if (cp[end]==0) return;
// eveything else is the block name 
//mjm_word(cp,start,pend);
StrTy bname=cp+end+1;
MM_ERR(" block start "<<MMPR4(cp,start,end,bname))
ViewingBlock & vb=m_vb_map.get(bname,m_last_pc,1);
//m_prior_vis=m_vp_map.prior_vis(bname);
//if (bname=="section2 ") vb.visibility(0);
//mjm_visibility=vb.visibility();
NewVisibility(vb.visibility());
m_active=&vb;
MM_ERR(MMPR2(bname, mjm_visibility))
if (mjm_visibility==0)
{
	MM_ERR(" trying to typset "<<MMPR(bname))
	 mjm_set_string(bname.c_str());
} else
{
// each block muyst start on left anywey 
mjm_xdvi_xr(0);

} 

MM_ERR("done with  block start "<<MMPR4(cp,start,end,bname))

return;
} // start 
// if typeset, it should have 
if (memcmp(cp+start, "blockend", end-start) == 0) {
//printf(" block end   %s end %d start %d \n ",cp,end,start);
//ToVisible(); 
 mjm_visibility=1;
MM_ERR(" block end "<<MMPR3(cp,start,end))
// don't allow null names 
if (cp[end]==0) return;
// eveything else is the block name 
//mjm_word(cp,start,pend);
StrTy bname=cp+end+1;
//ViewingBlock & vb=m_vb_map.get(bname);
ViewingBlock & vb=m_vb_map.get(bname,m_last_pc,2);
m_prior_vis=(vb.visibility());
m_active=0;

return;
} // end 
if (memcmp(cp+start, "visible", end-start) == 0) {
MM_ERR(" visible"<<MMPR3(cp,start,end))
ToVisible(); 
// mjm_visibility=1;
return; 
}
MM_ERR(" bad mjm special  "<<MMPR3(cp,start,end))

}

#endif



/*
Each block is a paragraph. 
Prior      Next       Action 
plain     visible     none 
plain     invisible   print crlf caption, push last
visible   plain       none 
visible   visible     none 
visible   invisible   print crlf caption, push last 
invisible plain       print crlf, then text
invisible visible     print crlf , then text 
invisible invisible   print caption, push last  

plain : current cursor not equal last. 
*/


void ToVisible()
{

//mjm_xdvi_xr(2) ;
if ( m_prior_vis==1 ) { mjm_visibility=1; return; } 
// last one was not displayed, either have prior as invsi or plain
const bool m=mjm_moved();
//if ( m) 
mjm_xdvi_xr(2) ;
//if ( mjm_visibility==0) mjm_xdvi_xr() ;
mjm_visibility=1;
}
void ToInVisible()
{
//mjm_xdvi_xr(2) ;
if ( m_prior_vis==0 ) { mjm_visibility=0; return; } 
//m_prior_vis=m_visibility;
//const bool m=mjm_moved();
//if ( mjm_visibility==1) 
//if (mjm_moved()) 
mjm_xdvi_xr(2) ;

mjm_visibility=0;
}

void NewVisibility( const IdxTy v )
{
if (v==0) ToInVisible();
else ToVisible();

}
bool MatchSet(ViewingBlock& vb, const StrTy & n, const IdxTy v,const bool pol )
{

//MM_ERR(MMPR3(vb.name().c_str(),n.c_str(),n.length()));  
const bool x= strncmp(vb.name().c_str(),n.c_str(),n.length())==0; 
if ( x&&pol || ( !x&&!pol)  ) vb.visibility(v);
return x; 
} 


IdxTy Show(const StrTy & n,const IdxTy flags ) { 
const bool pol=Bit(flags,0);
IdxTy cnt=0;
EnterSerial(USR_EVENT_MU);
MM_LOOP(ii,m_vb_map) { if ( MatchSet((*ii).second , n , 1,pol)) ++cnt;  } 
ExitSerial(USR_EVENT_MU);
return cnt; 

}
IdxTy Hide(const StrTy & n,const IdxTy flags) { 
const bool pol=Bit(flags,0);
IdxTy cnt=0;
EnterSerial(USR_EVENT_MU);
MM_LOOP(ii,m_vb_map) { if ( MatchSet((*ii).second , n , 0,pol)) ++cnt;  } 
ExitSerial(USR_EVENT_MU);

return cnt; }




IdxTy m_prior_vis,m_blockes_since_bop,m_inter; 

ViewingBlockMap m_vb_map;
IdxTy m_last_pc, m_max_pc;
ViewingBlock * m_active;
IdxTy m_flags,m_y_active;

}; // mjm_xdvik_interface

//////////////////////////////////////////////

template <class Tr>
class mjm_xdvik_interface_map : public std::map<typename Tr::StrTy, mjm_xdvik_interface< Tr > >  
{
 typedef mjm_xdvik_interface_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_xdvik_interface< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_xdvik_interface_map() {}






typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


StrTy dump(const IdxTy flags=0) { return Dump(flags); }

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




private:

}; // mjm_xdvik_interface_map




////////////////////////////////////////////
#ifdef  TEST_MJM_XDVIK_INTERFACE
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

class tester {
typedef mjm_cli_ui<tester> Cli;
public:
 void cli_cmd( Cli::list_type & choices,  const char * frag)
{
/*const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
*/
}

 void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

 }; // tester
typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_XDVIK_INTERFACE "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_xdvik_interface<Tr>  Myt;
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

#endif // MJM_XDVIK_INTERFACE_H__ 
