#ifndef MJM_SCATTER_SHATTER_H__
#define MJM_SCATTER_SHATTER_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_collections.h"
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


// Fri Feb 12 14:04:24 EST 2021
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_scatter_shatter   
// g++  -Wall -std=gnu++11 -DTEST_MJM_SCATTER_SHATTER -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_scatter_shatter.h  -lpthread -lreadline

template <class Tr>
class mjm_scatter_shatter 
{
 typedef mjm_scatter_shatter Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;
class _point
{

public:
_point():
m_x(0),m_y(0),m_r(0),m_drdx(0),m_drdy(0),m_drdxdy(0),m_drdx2(0),m_drdy2(0) {}

_point(const D & x, const D & y ):
m_x(x),m_y(y),m_r(0),m_drdx(0),m_drdy(0),m_drdxdy(0),m_drdx2(0),m_drdy2(0) {}

D curve() const { return sqrt(m_drdx2*m_drdx2+m_drdy2*m_drdy2); } 
D grad() const { return sqrt(m_drdx*m_drdx+m_drdy*m_drdy); } 
D val() const { return m_r; } 
D mix() const { return m_drdxdy; } 

StrTy dump() const
{
Ss ss;
ss<<MMPR4(m_x,m_y,m_r,m_drdx)<<MMPR4(m_drdy,m_drdx2,m_drdy2,m_drdxdy);
return ss.str();
} // dump
D m_x,m_y,m_r,m_drdx,m_drdy,m_drdxdy,m_drdx2,m_drdy2;

}; // _point
typedef _point point;


class _gauss_point
{
public:
_gauss_point(const D & x, const D & y, const D & w, const D & a)
: m_x(x),m_y(y),m_w(w),m_a(a) {Init(); } 
_gauss_point()
: m_x(0),m_y(0),m_w(0),m_a(0) {Init(); } 
StrTy dump(const IdxTy flags=0)
{
Ss ss;
ss<<MMPR4(m_x,m_y,m_w,m_a);
return ss.str();
}
D & rho(const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
return m_w*exp(-m_a*(dx*dx+dy*dy));
} // rho
D & drdx(const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
return -2*m_a*(dx)*m_w*exp(-m_a*(dx*dx+dy*dy));
} //drdx 
D & drdy(const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
return -2*m_a*(dy)*m_w*exp(-m_a*(dx*dx+dy*dy));
} //drdy 
D & drdx2(const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
return -2*m_a*(1.0-2*m_a*dx*dx)*m_w*exp(-m_a*(dx*dx+dy*dy));
} //drdx2
D & drdy2(const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
return -2*m_a*(1.0-2*m_a*dy*dy)*m_w*exp(-m_a*(dx*dx+dy*dy));
} //drdx 

D & dall(D & drdx, D & drdy, D & drdx2, D& drdy2, const D & x, const D & y) const
{
const D dx=x-m_x;
const D dy=y-m_y;
const D p=m_w*exp(-m_a*(dx*dx+dy*dy));
const D p2=-2.0*m_a*p;
drdx=p2*dx;
drdy=p2*dy;
drdx2=p2*(1.0-2*m_a*dx*dx);
drdy2=p2*(1.0-2*m_a*dy*dy);
return p;
} // dall

template <class Tp> IdxTy  dall(Tp & xx) // , const D & x, const D & y) const
{
const D x=xx.m_x;
const D y=xx.m_y;
const D dx=x-m_x;
const D dy=y-m_y;
const D p=m_w*exp(-m_a*(dx*dx+dy*dy));
const D p2=-2.0*m_a*p;
xx.m_drdx=p2*dx;
xx.m_drdy=p2*dy;
xx.m_drdxdy=p2*(-2.0*m_a)*dx*dy;
xx.m_drdx2=p2*(1.0-2*m_a*dx*dx);
xx.m_drdy2=p2*(1.0-2*m_a*dy*dy);
xx.m_r= p;
return 0;
} // dall


//template <class Tp> IdxTy  sall(Tp & xx, const D & x, const D & y) const
template <class Tp> IdxTy  sall(Tp & xx) const
{
const D x=xx.m_x;
const D y=xx.m_y;
const D dx=x-m_x;
const D dy=y-m_y;
const D p=m_w*exp(-m_a*(dx*dx+dy*dy));
const D p2=-2.0*m_a*p;
xx.m_drdx+=p2*dx;
xx.m_drdy+=p2*dy;
xx.m_drdx2+=p2*(1.0-2*m_a*dx*dx);
xx.m_drdy2+=p2*(1.0-2*m_a*dy*dy);
xx.m_drdxdy+=p2*(-2.0*m_a)*dx*dy;
xx.m_r+= p;
return 0;
} // sall






private:
void Init()
{

} // Init

D m_x,m_y,m_w,m_a;

}; // _gauss_point

typedef _gauss_point gauss_point;

typedef std::vector<gauss_point> PointVec;

public:
mjm_scatter_shatter() {Init(); }
~mjm_scatter_shatter() {}
void grid( const IdxTy nx, const IdxTy ny, const D & dx, const D & dy, const D & xi, const D & yi)
//void grid( const IdxTy nx, const IdxTy ny, const D & xi, const D & yi)
{ Grid(nx,ny,dx,dy,xi,yi); }

void add_point(const D & x, const D & y, const D & w, const D & a)
{Point(x,y,w,a); } 
void solve(const IdxTy flags) { Solve(flags); } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  
MM_LOOP(ii,m_points) ss<<(*ii).dump(flags)<<CRLF;

return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
void Solve(const IdxTy flags) { 
const D xc=m_xt/m_n;
const D yc=m_yt/m_n;
const D meana=m_a/m_n;
D x=xc, y=yc;
D xbox=m_xmax-m_xmin;
D ybox=m_ymax-m_ymin;
point d(x,y);

Eval(d,0);
MM_ERR(MMPR(d.dump()))
if (d.grad()<1e-8)
{
 x=x+.01*meana;
 //x=x+.001*xbox;
// y=y+.001*ybox;
 y=y+.01*meana;
}
else // move perpendicular to grad OR find the max... 
{
// move along grad until ok
D f=.001;
x=x+d.m_drdx*f;
y=y+d.m_drdy*f;

}
d=point(x,y);
Eval(d,0);
MM_ERR(MMPR(d.dump()))
// get back in the saddle 
while ( true)
{
D f=.00001;
D gv=d.grad();
D dxp=d.m_drdx*f/gv;
D dyp=d.m_drdy*f/gv;
//d=point(x,y);
point dplus(x+dxp,y+dyp);
point dminus(x-dxp,y-dyp);
Eval(dplus,0);
Eval(dminus,0);
MM_ERR(MMPR(dplus.dump()))
MM_ERR(MMPR(dminus.dump()))
MM_ERR(MMPR4(x,y,(dplus.val()-d.val()),(dminus.val()-d.val())))
if ( dplus.val()<d.val()) if (dminus.val()<d.val()) break;
if (dplus.val()>d.val()) { x=x+dxp; y=y+dyp;}
else if (dminus.val()>d.val()) { x=x-dxp; y=y-dyp;}
d=point(x,y);
Eval(d,0);
MM_ERR(MMPR(d.dump()))
} // top of ridge 
for(IdxTy iwalk=0; iwalk<0; ++iwalk)
{
// now go perpendicular to the gradient and keep verifying the position...
D gv=d.grad();
D f=.03;
D dxp=-d.m_drdy*f/gv;
D dyp=d.m_drdx*f/gv;
x=x+dxp; y=y+dyp;
d=point(x,y);
Eval(d,0);
MM_ERR("XXXXXXXXXXXXXXXXXXXXXX " << MMPR(d.dump()))
} // iwalk 

} // Solve
void Grid( const IdxTy nx, const IdxTy ny, const D & dx, const D & dy, const D & xi, const D & yi)
{ 
for(IdxTy i=0; i<nx; ++i)
{
for(IdxTy j=0; j<ny; ++j)
{
point d(xi+i*dx, yi+j*dy);
Eval(d,0);
MM_ERR(MMPR(d.dump()));
} // j

} // i 


} // Grid

 

void Eval( point & d, const IdxTy flags)
{
MM_LOOP(ii,m_points) { (*ii).sall(d); }

}
void Point(const D & x, const D & y, const D & w, const D & a)
{
if (m_n==0) { m_xmax=x; m_xmin=y; m_ymin=y; m_ymax=y; } 
++m_n;
m_xt+=x;
m_yt+=y;
m_wt+=w;
m_a+=a;
if (x>m_xmax){ m_xmax=x;}
if (y>m_ymax){ m_ymax=y;}
if (x<m_xmin){ m_xmin=x;}
if (y<m_ymin){ m_ymin=y;}

m_points.push_back(gauss_point(x,y,w,a));

}
void Init()
{
m_xt=0;
m_yt=0;
m_wt=0;
m_a=9;
m_n=0;
 m_xmax=0;
m_xmin=0;
m_ymin=0;
m_ymax=0;
} // Init


// MEMBERS
D m_xt,m_yt,m_wt,m_a;
D m_xmax,m_xmin,m_ymin,m_ymax;
IdxTy m_n;
Ragged m_config;
PointVec m_points;

}; // mjm_scatter_shatter

//////////////////////////////////////////////

template <class Tr>
class mjm_scatter_shatter_map : public std::map<typename Tr::StrTy, mjm_scatter_shatter< Tr > >  
{
 typedef mjm_scatter_shatter_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_scatter_shatter< Tr> >   Super;
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
mjm_scatter_shatter_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
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

}; // mjm_scatter_shatter_map




////////////////////////////////////////////
#ifdef  TEST_MJM_SCATTER_SHATTER
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
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_scatter_shatter <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_SCATTER_SHATTER "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

typedef double D;
int main(int argc,char **args)
{
about();
typedef mjm_scatter_shatter<Tr>  Myt;
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
if (cmd=="new") {x=Myt();  MM_ERR(x.dump()) }
if (cmd=="solve") {x.solve(0); MM_ERR(x.dump()) }
if (cmd=="point") { 
const D xx=atof(cip.wif(1).c_str()); 
const D y=atof(cip.wif(2).c_str()); 
const D w=atof(cip.wif(3).c_str()); 
const D a=atof(cip.wif(4).c_str()); 
x.add_point(xx,y,w,a);
MM_ERR(x.dump()) }
//void grid( const IdxTy nx, const IdxTy ny, const D & dx, const D & dy, const D & xi, const D & yi)
if (cmd=="grid") { 
const IdxTy nx=atof(cip.wif(1).c_str()); 
const IdxTy ny=atof(cip.wif(2).c_str()); 
const D dx=atof(cip.wif(3).c_str()); 
const D dy=atof(cip.wif(4).c_str()); 
const D xi=atof(cip.wif(5).c_str()); 
const D yi=atof(cip.wif(6).c_str()); 
x.grid(nx,ny,dx,dy,xi,yi);
MM_ERR(x.dump()) 
}


//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_SCATTER_SHATTER_H__ 