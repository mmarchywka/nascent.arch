#ifndef MJM_COLOR_TABLE_H__
#define MJM_COLOR_TABLE_H__

#include "mjm_globals.h"
//#include "mjm_data_model_error_log.h"
//#include "mjm_generic_iterators.h"
#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
//#include "mjm_block_matrix.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_collections.h"
//#include "mjm_svg_writer.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
// #define COLOR_OUT typedef char Ch; const Ch  e1[]={Ch(27),'[','4','4','m',0}; const Ch  e1a[]={Ch(27),'[','4','2','m',0}; const Ch  e1b[]={Ch(27),'[','4','1','m',0}; const Ch  e2[]={Ch(27),'[','0','m',0}; StrTy _blue=StrTy(e1); StrTy _green=StrTy(e1a); StrTy _red=StrTy(e1b); StrTy _end=StrTy(e2);


namespace color_table_traits
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
//typedef  data_model_error_log Dmel;
}; // 
}; // color_table_traits
//template <class Tobj, class Traits >
class mjm_color_table 
{
typedef mjm_color_table Myt;
protected:
typedef color_table_traits::Tr  Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
//typedef typename Tr::MyBlock  MyBlock;
//typedef typename Tr::Dmel  Dmel;

typedef std::map<IdxTy, StrTy  > Cm;
typedef std::map<IdxTy, IdxTy  > Fm;
typedef Cm::const_iterator CmItor;
typedef Fm::const_iterator FmItor;
typedef std::vector<IdxTy> Strats;

class Env
{
public:
Env(const Cm & c,const Fm & fm,  const IdxTy b, const IdxTy m, const StrTy& d, const
Strats &s): m_colmap(c),m_flagmap(fm), m_bounds(b),m_mod(m),m_default(d), m_states(s) {}

Cm m_colmap;
Fm m_flagmap;
IdxTy m_bounds;
IdxTy m_mod;
StrTy m_default;
Strats m_states;
}; // Env


public:
mjm_color_table():m_bounds(0)  , m_mod(0),m_default("grey")  {Init(); }
virtual ~mjm_color_table() {}


//#define COLOR_OUT typedef char Ch; const Ch  e1[]={Ch(27),'[','4','4','m',0}; const Ch  e1a[]={Ch(27),'[','4','2','m',0}; const Ch  e1b[]={Ch(27),'[','4','1','m',0}; const Ch  e2[]={Ch(27),'[','0','m',0}; StrTy _blue=StrTy(e1); StrTy _green=StrTy(e1a); StrTy _red=StrTy(e1b); StrTy _end=StrTy(e2);
// VT100
static StrTy vt100_lut(const StrTy & x )
{
static bool loaded=false;
static std::map<StrTy,StrTy> vt100lut;
if (!loaded)
{
typedef char Ch; //const Ch  e1[]={Ch(27),'[','4','4','m',0}; 
//const Ch  e1a[]={Ch(27),'[','4','2','m',0}; 
//const Ch  e1b[]={Ch(27),'[','4','1','m',0}; 
const Ch  e2[]={Ch(27),'[','0','m',0}; 
//StrTy _blue=StrTy(e1); 
//StrTy _green=StrTy(e1a); 
//StrTy _red=StrTy(e1b); 
StrTy _end=StrTy(e2);
//vt100lut["red"]=_red;
//vt100lut["green"]=_green;
//vt100lut["blue"]=_blue;
//vt100lut["end"]=_end;
// FUDDING strings now the zero is included doh 
vt100lut["end"]=_end=StrTy({Ch(27),'[','0','m'}); 
vt100lut["blackb"]=StrTy({Ch(27),'[','4','0','m',0}); 
//vt100lut["redb"]=StrTy({Ch(27),'[','4','1','m',0}); 
vt100lut["redb"]=StrTy({Ch(27),'[','4','1','m'}); 
vt100lut["greenb"]=StrTy({Ch(27),'[','4','2','m',0}); 
vt100lut["yellowb"]=StrTy({Ch(27),'[','4','3','m',0}); 
vt100lut["blueb"]=StrTy({Ch(27),'[','4','4','m',0}); 
vt100lut["magentab"]=StrTy({Ch(27),'[','4','5','m',0}); 
vt100lut["cyanb"]=StrTy({Ch(27),'[','4','6','m',0}); 
vt100lut["whiteb"]=StrTy({Ch(27),'[','4','7','m',0}); 

vt100lut["blackf"]=StrTy({Ch(27),'[','3','0','m',0}); 
vt100lut["redf"]=StrTy({Ch(27),'[','3','1','m'}); 
vt100lut["greenf"]=StrTy({Ch(27),'[','3','2','m',0}); 
vt100lut["yellowf"]=StrTy({Ch(27),'[','3','3','m',0}); 
vt100lut["bluef"]=StrTy({Ch(27),'[','3','4','m',0}); 
vt100lut["magentaf"]=StrTy({Ch(27),'[','4','5','m',0}); 
vt100lut["cyanf"]=StrTy({Ch(27),'[','3','6','m',0}); 
vt100lut["whitef"]=StrTy({Ch(27),'[','3','7','m',0}); 




loaded=true;
} // loaded 
//return vt100lut[x]; // this does pollute it with 
auto ii=vt100lut.find(x);
if (ii==vt100lut.end()) return StrTy();
return (*ii).second;

} //vt100_lut

static StrTy vt100(const StrTy & x, const StrTy & c, const IdxTy flags)
{
const StrTy & pfx=vt100_lut(c);
const StrTy & end=vt100_lut("end");
if (!pfx.length()) 
{
MM_ERR(" no vt100 code for "<<MMPR3(x,c,flags))
return x;
}
return pfx+x+end;
}// vt100 


const StrTy & operator()(const IdxTy i )  const
{
CmItor ci=m_colmap.find(i); if (ci!=m_colmap.end()) return (*ci).second; 
MM_LOOP(ii,m_states)
{
switch (*ii)
{
case 0:{ if (m_mod==0) continue; CmItor ci=m_colmap.find(i%m_mod); 
		if (ci!=m_colmap.end()) return (*ci).second;  break ; }
//case 1:{ if (m_mod==0) continue; CmItor ci=m_colmap.find(ialt); 
//		if (ci!=m_colmap.end()) return (*ci).second;  break ; }
default: break; 
} // switch
} // ii 
return m_default;
}

// TODO return a ref except for stupid zero at tail 
const IdxTy  flags(const IdxTy i )  const
{
FmItor ci=m_flagmap.find(i); if (ci!=m_flagmap.end()) return (*ci).second; 
MM_LOOP(ii,m_states)
{
switch (*ii)
{
case 0:{ if (m_mod==0) continue; FmItor ci=m_flagmap.find(i%m_mod); 
		if (ci!=m_flagmap.end()) return (*ci).second;  break ; }
//case 1:{ if (m_mod==0) continue; CmItor ci=m_colmap.find(ialt); 
//		if (ci!=m_colmap.end()) return (*ci).second;  break ; }
default: break; 
} // switch
} // ii 
return 0; // m_default;
}


void load(const StrTy & fn, const IdxTy flags ) 
{ std::ifstream ifs(fn); load(ifs,flags);
MM_ERR(MMPR3(fn,m_colmap.size(),flags))

 } 

void load(std::istream  & is, const IdxTy flags) 
{
    CommandInterpretter li(&is);
//    li.set_split(1,'|');
	IdxTy line=0;
    while (li.nextok())
    {
		MM_ERR(" fudd "<<li.line())
        const IdxTy sz=li.size();
        if (sz<1) continue;
		const IdxTy mx=m_colmap.size();
 //       if (li.line().size()!=0) ts.push_back(li.word(1));
		if (sz==1) m_colmap[mx]=li.word(0);
		if (sz==2) m_colmap[myatoi(li.word(0))]=li.word(1);
		if (sz==3){
			const IdxTy loc=myatoi(li.word(0));
			 m_colmap[loc]=li.word(1);
			m_flagmap[loc]=myatoi(li.word(2));
			}
		++line; 
   } // nextok()


} // load 

StrTy color(const D & val, const D & zmin, const D & zmax) const
{
D rho=(val-zmin)/(zmax-zmin);
D red=1.0-2.0*rho;
if (red<0) red=0;
D blue=2.0*rho-1.0;
if (blue<0) blue=0;
D green=(rho-.75)*(rho-.25);
if (green>0) green=0;
else { green=sqrt(-16.0*green); }
StrTy col="";
color(col,red,green,blue,0);
return col;

}

static const char  lut(const IdxTy p) // const
{
static char c[16];
static bool init=false;
if (!init)
{
for(IdxTy i=0; i<10; ++i) c[i]='0'+i;
for(IdxTy i=0; i<6; ++i) c[i+10]='A'+i;

}
return c[p];
}
static void myhex(Ss & ss, const char c)
{
ss<<lut((c>>4)&15);
ss<<lut(c&15);
}
// wth is this 255 times a float??}
StrTy myhex(const D & x) const
{
Ss ss;
const IdxTy n=IdxTy(x*255);
ss<<(lut((n>>4)&15));
ss<<(lut((n>>0)&15));
return ss.str();
}



void color(StrTy & col,const D & red,const D & green,const D & blue,const IdxTy & flags)
const
{
Ss ss;
ss<<"#0";
ss<<myhex(red);
ss<<myhex(green);
ss<<myhex(blue);
col=ss.str();
}
void clear() { Init();}
void set(const IdxTy i, const StrTy & col) { m_colmap[i]=col;}

void add(const StrTy & col)
{
	const IdxTy mx=m_colmap.size();
	set(mx,col);
}

void push()
{ m_stack.push_back(Env(m_colmap,m_flagmap, m_bounds, m_mod, m_default, m_states)); 
}
void pop()
{
const Env & x=m_stack.back();
m_colmap=x.m_colmap;
m_flagmap=x.m_flagmap;
m_bounds=x.m_bounds;
m_default=x.m_default;
m_states=x.m_states;
m_stack.pop_back();
}
void some_colors() { SomeColors(); }
IdxTy size() const { return m_colmap.size(); } 



void cmd(const StrTy & cmd, const IdxTy flags)
{ CommandInterpretter li; li.set(cmd,1); 

command_mode(li); 

}
void command_mode(CommandInterpretter &li)
{
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="some_colors") { some_colors(); continue; }
if (cmd=="clear") { clear(); continue; }
if (sz==1)
{
if (cmd=="add") { add(li.word(1)); continue; }
if (cmd=="load") { load(li.word(1),0); continue; }
} // 1 
if (sz==3)
{
if (cmd=="set") { set(myatoi(li.word(1)),li.word(2)); continue; }
if (cmd=="load") { load(li.word(1),myatoi(li.word(2))); continue; }
} // 1 


}

}


protected:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }
void Init()
{
m_colmap.clear();
m_states.clear();
m_states.push_back(0);
}
void SomeColors()
{
add("white");
add("red");
add("green");
add("blue");
add("black");
add("yellow");
m_mod=m_colmap.size();
}




Cm m_colmap;
Fm m_flagmap;
IdxTy m_bounds;
IdxTy m_mod;
StrTy m_default;
Strats m_states;

std::vector<Env> m_stack;



}; // mjm_align_collection


#endif

