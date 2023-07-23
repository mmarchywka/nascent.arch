#ifndef MJM_VIEW_STRUCT_ALPHA_H__
#define MJM_VIEW_STRUCT_ALPHA_H__

// MM_ERR messages swallowed when c code uses stderr
#define MIXED_CERR_STREAMS

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_alt_display.h"
#include "mjm_cursing_viewer.h"
//#include "mjm_2d_states.h"

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


#ifdef HAVE_XDVIK
#warning including xdvik

#include "mjm_xdvik_interface.h"
#if 0 
extern "C" { 
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


};
#endif

#endif



// Wed Jun  5 13:11:54 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_view_struct_alpha   
// g++ -std=gnu++11 -DTEST_MJM_VIEW_STRUCT_ALPHA -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_view_struct_alpha.h  -lpthread -lreadline -lncurses
/*
dvixxx="/home/ubuntu/dev/xdvik/xdvik-22.87.03"
dvia="$dvixxx/texk/xdvik/libxdvi.a"
dvib="$dvixxx/texk/kpathsea/.libs/libkpathsea.a"
dvic="$dvixxx/libs/freetype2/libfreetype.a"
La=" -L$dvixxx/texk/xdvik "
la=" -lxdvi "
Lb=" -L$dvixxx/texk/kpathsea/.libs -lkpathsea"
Lc=" -L$dvixxx/libs/freetype2 -lfreetype"
h1="-I$dvixxx/texk/xdvik"
h2="-I$dvixxx/texk"
h3="-I$dvixxx/texk/xdvik/gui"
h4="-I$dvixxx/libs/freetype2/freetype2"
# permissive for the stupid main.c copy
g++ -std=c++11  -fpermissive  -DHAVE_XDVIK -DTEST_MJM_VIEW_STRUCT_ALPHA -I. -I../../mjm/hlib -I../../mjm/num $h1 $h2 $h3 $h4  -gdwarf-3 -O0   $dvib $dvic  -x c++ mjm_view_struct_alpha.h -x $dvia  $La -lxdvi  $la $Lb $Lc   -lpthread -lreadline -lncurses -lXaw -lXmu -lXt -lSM -lICE -lXi -lXext -lXpm -lX11 -lm

*/




template <class Tr>
class mjm_view_struct_alpha 
{
 typedef mjm_view_struct_alpha Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef mjm_alt_display<Tr> AltDisp;
typedef mjm_cursing_viewer<Tr> Viewer;
//typedef mjm_2d_states<Tr,IdxTy> ScreenPos;
typedef typename AltDisp::El El;

// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_view_struct_alpha():m_dflags(0) {}
~mjm_view_struct_alpha() {}
void view_file(const StrTy & fn) { ViewFile(fn); } 
void view_file_interactive(const StrTy & fn) { ViewFileInteractive(fn); } 
StrTy make_string() { return MakeString(); } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const IdxTy flags=0) {Ss ss; 

ss<<m_v.dump();
ss<<m_ad.dump();

 return ss.str(); }
void ViewFile(const StrTy & fn) 
{ 
 m_ad.load_file(fn); 
// with current states, get the view
const StrTy x=m_ad.display();
Ss ss;
ss<<x;
m_v.display(ss);
} // ViewFile 

int NextNotOf( const void * elp, const El & el, IdxTy & cxu, IdxTy & cyu, IdxTy posi, const bool er )
{
const bool exclude_root=er; // true;
IdxTy pe=0;
const IdxTy els=el.size();
IdxTy kluge=posi; 
// doh, stupid kluge assumed contiguous things. Needs to start at current location 
for(; kluge<els; ++kluge) { if (el[kluge]==elp) 
{
if (exclude_root) {
while  ((el[kluge]==elp)||(el[kluge]==&m_ad)){  ++kluge; if (kluge==els) { kluge=0; break; }} 
} else {
while  ((el[kluge]==elp)){  ++kluge; if (kluge==els) { kluge=0; break; }} 
}

IdxTy pe=kluge;
//IdxTy cxu=0;
//IdxTy cyu=0;
return (m_v.cursor_pos(cxu,cyu,pe)); // { //kcx=cxu; cy=cyu; } 
//break;
}}
return ~0;
} 

StrTy MakeString()
{
El el;
const StrTy x=m_ad.display(el,m_dflags);
return x;
}
void ViewFileInteractive(const StrTy & fn) 
{ 
 m_ad.load_file(fn); 
int  xi=0;
int cx=0;
int cxx=0;
int cy=0;
int cyx=0;
IdxTy BAD=~0;
IdxTy epos=BAD;
IdxTy pos=0;
void * elp=0;
bool er=true;
//IdxTy dflags=0;
IdxTy&  dflags=m_dflags;
while (true) { 
El el;
// with current states, get the view
// MVC lol? 
// maintain list of element to string pos tables
//const StrTy x=m_ad.display(el,xi);
const StrTy x=m_ad.display(el,dflags);
// get string pos vs (x,y) 
int xi2=  m_v.display_and_curse(x,xi,cxx,cyx);
if (epos!=BAD)
{
IdxTy cxu=0;
IdxTy cyu=0;
// if this fails, the content needs to scroll 
if ( 0==NextNotOf(elp,el,cxu,cyu,pos,er) ) { cx=cxu; cy=cyu; } 
// after redisplay this needs to be fixed, 
else {xi=xi2; cy=0;
el.clear();
// now x is not consistent 
StrTy xs=
m_ad.display(el,xi);
xi2=  m_v.display_and_curse(xs,xi,cxx,cyx);
if ( 0==NextNotOf(elp,el,cxu,cyu,pos,er) ) { cx=cxu; cy=cyu; } 
else {MM_ERR(" scroll fails to work")}
} 
// set the cursor to string position for element pos 
epos=BAD;
} // epos!=BAD
// keypress on string pos translates to element
// get char and pos
IdxTy ch;; // ,pos;
const IdxTy rc= m_v.cursing_pos_tab(ch,pos,cx,cy,0);
//MM_ERR(" returns with "<<MMPR3(ch,pos,rc))
if (pos<el.size())
{
//MM_ERR(MMPR(el[pos]->name())) MM_ERR(MMPR(x[pos]))
// allowing this for m_ad means repated - moves cursor backwards ... 
//if ( el[pos]!=&m_ad)
{
if (ch=='+') el[pos]->set_level("1");
// collapsing this alone removes crlf 
if ( el[pos]!=&m_ad)
if (ch=='-') {elp=el[pos];  el[pos]->set_level("0"); epos=pos; er=true; } 
}
if (ch=='>') el[pos]->cmd("cycle");
if (ch=='%') m_ad.cmd("cycle");
if (ch=='=') el[pos]->cmd("cycle_same");
if (m_v.maskac(rc)==2) {elp=el[pos];  if ( el[pos]->level()=="0")  epos=pos; er=false; } 
}
if (ch=='!'){  el[pos]->cmd("cycle_all"); m_ad.set_level("1"); } 
if (ch=='1'){   m_ad.cmd("++"); } 
if (ch=='0'){   m_ad.cmd("--"); } 
if (ch=='@'){   el[pos]->cmd("next"); } 

if (ch=='q') break;
if (ch=='t'){ xi=0; cx=0; cy=0;  } 
if (ch=='D'){ dflags^=1;  } 
if (ch=='T'){ dflags^=2;  } 

// this pages, want to scroll 
//if ((rc&4)!=0) { xi=xi2;  cy=0;  }
if ((rc&4)!=0) { int xinew=m_ad.next_line(x,xi+1); if (xinew!=BAD) xi=xinew;   }


//if ((rc&8)!=0) { int xinew=xi-80; if (xinew<0) xinew=0;  xi=xinew; cx=0; cy=0;  }
// started on line break need to go back 
if ((rc&8)!=0) { int xinew=m_ad.prior_line(x,xi-2-2); if (xinew<0) xinew=0; if( xinew==BAD) xi=0;   xi=xinew; cx=cx; cy=0;  }

//
} // while 

} // ViewFile 


AltDisp m_ad;
Viewer m_v;
IdxTy m_dflags;
}; // mjm_view_struct_alpha

//////////////////////////////////////////////

template <class Tr>
class mjm_view_struct_alpha_map : public std::map<typename Tr::StrTy, mjm_view_struct_alpha< Tr > >  
{
 typedef mjm_view_struct_alpha_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_view_struct_alpha< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_view_struct_alpha_map() {}
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

}; // mjm_view_struct_alpha_map




////////////////////////////////////////////
#ifdef  TEST_MJM_VIEW_STRUCT_ALPHA
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
ss<<" MJM_VIEW_STRUCT_ALPHA "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_view_struct_alpha<Tr>  Myt;
//Myt x(argc,args);
Myt x;
#ifdef HAVE_XDVIK
typedef mjm_xdvik_interface<Tr> Dvi;
Dvi dvi;
#endif

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
Ss ss;
for(IdxTy i=1; i<argc; ++i) ss<<args[i]<<CRLF;
li.push(&ss, false);
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
if (cmd=="view") { x.view_file_interactive(cip.p1); }
if (cmd=="save") { MM_ERR(" writing to "<<MMPR(cip.p1)) std::ofstream os(cip.p1);  os<<x.make_string(); }


#ifdef HAVE_XDVIK
if (cmd=="good-dvi") { MM_ERR(is_good_dvi_file(cip.p1.c_str(),false)) }
if (cmd=="launch-dvi") { dvi.launch_page(cip.p1); ;  }
if (cmd=="start-dvi") { main_xdvi(argc,args);  }
if (cmd=="sig") {dvi.signal(); MM_ERR(" signalled ")  }
if (cmd=="hide") {IdxTy n=dvi.hide(cip.p1,atoi(cip.p2.c_str())); dvi.redraw(); MM_ERR(" hid "<<n)  }
if (cmd=="show") {IdxTy n=dvi.show(cip.p1,atoi(cip.p2.c_str())); dvi.redraw(); MM_ERR(" show  "<<n)  }
if (cmd=="x") {StrTy x=dvi.dump(0); MM_ERR("dumping "<< x)  }
if (cmd=="redraw") { dvi.redraw(); MM_ERR(" redraw ")  }
if (cmd=="wait-load") { 
MM_ERR(" wait-load")  
dvi.wait_eop(); 
MM_ERR(" load wait  return  ")  
// could redef macro lol 
fprintf(stderr," printf version,  load wait  return     \n") ;
}
if (cmd=="p") { mjm_paging=atoi(cip.p1.c_str()); MM_ERR(" paging "<<MMPR(mjm_paging))  }
#endif


//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }
MM_ERR(" done with "<<cmd)
} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_VIEW_STRUCT_ALPHA_H__ 