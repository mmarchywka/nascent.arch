#ifndef MJM_GAUSS_H__
#define MJM_GAUSS_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// try to use the home made version of erf(x1)-erf(x2)
#include "mjm_integrals.h"
#include "mjm_sparse_2nomial.h"

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
Likely the origin shifts are not right between x and xp. Check all of this

Ths point of the class is to do 2D gaussian integrals integrating flow from source
to dest areas with Maxwell Boltzman velocity profile. Everything else is just for
testing and could be ifdeffed out after testing. 
*/

/*

For missing CXXABI, on my system it need to point to the right location
when using the cool stuff in mjm_integrals ( see its compile line)

https://github.com/FoldingAtHome/fah-issues/issues/1147

g++ -DTEST_GAUSS__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -x c++ mjm_gauss.h 




g++ -DTEST_BF_FLOW__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_gauss.h
g++ -DTEST_GAUSS__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_gauss.h
*/

////////////////////////////////////////////////////////////////

class mjm_gauss
{
/*

*/

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef mjm_gauss Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;

typedef mjm_logic_base Logic;

public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
mjm_gauss() {Init();}
mjm_gauss(int argc,char **_args) // : m_save_sol(false)
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<IdxTy(argc); ++i) args[i]=_args[i];
for (IdxTy i=argc; i<IdxTy(ikluge); ++i) args[i]=&dummy[0];

//m_size=1;
//m_points=2000; // this really should not have a default.. 
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
dest=::atoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
source("-source",i,argc,args);
quit("-quit",i,argc,args);

//if ( s==StrTy("-harvest")) { ++i;
//dest=::atoi(args[i]);
//const StrTy fn=StrTy(args[i]);
//process_file(fn);
//if (confirm) MM_ERR(" done harvesting "<<fn)
//++i; // consume param and cmd
//return;
//}
//if ( s==StrTy("-htext")) { ++i;
//if(confirm) MM_ERR(" output as htext ") // <<s<<" to "<<dest)
//StrTy s=m_chart.to_ascii_h();
//std::cout<<s<<CRLF;
//++i; // consume param and cmd
//return;
//}

#if 0 
if ( s==StrTy("-vtext")) { ++i;
//if(confirm) MM_ERR(" setting "<<s<<" to "<<dest)
if(confirm) MM_ERR(" output as vtext ") // <<s<<" to "<<dest)
StrTy s=m_chart.to_ascii_v();
std::cout<<s<<CRLF;
//++i; // consume param and cmd
return;
}

if ( s==StrTy("-cmd")) { ++i;
const StrTy cmd=StrTy(args[i]);
if(confirm) MM_ERR(" command mode exec of  "<<cmd) // <<s<<" to "<<dest)
command_mode(cmd); 
++i; // consume param and cmd
return;
}

// tl 
if ( s==StrTy("-hlatex")) { arg_cmd(i,args,0,"hlatex",confirm); return; }
if ( s==StrTy("-banner")) { arg_cmd(i,args,0,"banner",confirm); return; }
#endif

} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 

 void source( const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
if (argc<=i) return;
const StrTy nm=StrTy(args[i]);
MM_ERR(" sourcng  "<<nm<<"  ")
source(nm);
++i; // consume param and cmd
}
}
 void quit( const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return;
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
clean_up();
//if (argc<=i) return;
//const StrTy nm=StrTy(args[i]);
//MM_ERR(" sourcng  "<<nm<<"  ")
//source(nm);
//++i; // consume param and cmd
}
}


void dump_unused()
{
Ss ss;
//for (auto ii=m_unused.begin(); ii!=m_unused.end(); ++ii)
{
//ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
MM_MSG("unused:"<<CRLF<<ss.str())

}
void source(const StrTy & fn)
{

CommandInterpretter li;
li.source(fn,true);
// so then this clverly connects to cin... doh 
command_mode(li);
}


void command_mode()
{
//LineIterator 
CommandInterpretter li(&std::cin);
command_mode(li);

}
void command_mode(const StrTy & cmd)
{

CommandInterpretter li;
li.set(cmd,1);
//li.set(cmd,0);
command_mode(li);
}


void command_mode(CommandInterpretter & li)
{
StrTy local_label="fick";
while (li.nextok())
{
const IdxTy sz=li.size();
// MM_ERR(" processing "<<li.dump())
if (sz<1) continue;
const StrTy cmd=li.word(0);
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
//if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
//if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
//if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
//if (cmd=="process") { if (li.cmd_ok(2)) process_file(li.word(1));  continue; } 
//if (cmd=="add") { if (li.cmd_ok(2)) add_timeline(li.word(1));  continue; } 
if (cmd=="init") { Init();  continue; } 
//if (cmd=="hlatex") { std::cout<<m_chart.to_latex_h();  continue; } 
//if (cmd=="hssv") { std::cout<<m_chart.to_ssv(true);  continue; } 
//if (cmd=="vssv") { std::cout<<m_chart.to_ssv(!true);  continue; } 
//if (cmd=="add-event") { if (li.cmd_ok(3)) add_event(li.word(1),li.word(2));  continue; } 
//if (cmd=="add-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
//if (cmd=="set-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
//if (cmd=="reset-sum") { if (li.cmd_ok(2)) reset_sum(li.word(1));  continue; } 
//if (cmd=="make") { make_chart();  continue; } 
//if (cmd=="unused") { dump_unused();  continue; } 
//if (cmd=="clear-unused") { m_unused.clear();  continue; } 
//if (cmd=="legend") { dump_legend(std::cout,local_label,1);  continue; } 
//if (cmd=="dump-human") { dump_state_etc(std::cout,local_label,0);  continue; } 
//if (cmd=="init") { Init();   continue; } 
//if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
//if (cmd=="banner") { config_banner();  continue; } 
//if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
/*
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
*/
if (cmd=="quit") { clean_up(); return; } 
if (cmd!="") 
{
if (cmd.c_str()[0]!='#')
	MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
}

}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
/*
note the noise floor, 
jm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=77 beta=0.01 vxt=47 res=2.70887e-10 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=78 beta=0.01 vxt=48 res=1.00684e-10 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=79 beta=0.01 vxt=49 res=3.66356e-11 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=80 beta=0.01 vxt=50 res=1.30882e-11 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=81 beta=0.01 vxt=51 res=4.59011e-12 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=82 beta=0.01 vxt=52 res=1.56319e-12 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=83 beta=0.01 vxt=53 res=4.9738e-13 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=84 beta=0.01 vxt=54 res=1.13687e-13 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=85 beta=0.01 vxt=55 res=5.68434e-14 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=86 beta=0.01 vxt=56 res=1.13687e-13 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=87 beta=0.01 vxt=57 res=0 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=88 beta=0.01 vxt=58 res=0 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=89 beta=0.01 vxt=59 res=-7.10543e-14 s1=0 s0=1 c1=0 c0=1
mjm_gauss.h261  xpl=0 xpu=10 xl=10 xu=20 i=90 beta=0.01 vxt=60 res=0 s1=0 s0=1 c1=0 c0=1



*/
void test_int( CommandInterpretter & li )
{
D beta=.01;
D xpu=10;
D xpl=0;
D xu=20;
D xl=10;
D vxt=-30;
D s1=0;
D s0=1;
D c1=0;
D c0=1;
for (IdxTy i=0; i<100; ++i)
{
D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
MM_MSG(MMPR4(xpl,xpu,xl,xu)<<MMPR4(i,beta,vxt,res)<<MMPR4(s1,s0,c1,c0))
vxt+=1;

}


}


//D evaluateq( const D & xpc, const D & ypc, const D & xc, const D & yc
//, const D & vxt, const D & vyt, const D & a, const D & beta,const IdxTy quad, const IdxTy quad) 

void test_quad( CommandInterpretter & li )
{
D beta=1;
D xpc=0;
D ypc=0;
D xc=0;
D yc=0;
D vxt=-5;
D vyt=0;
D a=1;
D b=1;
const IdxTy quad=0;
const IdxTy quadp=0;

for (IdxTy i=0; i<100; ++i)
{
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad) ;
MM_MSG(MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(i,res, quadp, quad)) 
vxt+=.1;

}


}

D xrand(const D & mmax, const D & mmin) const
{
const IdxTy mask=(1<<16)-1;
const D de=((D)(rand()&mask))/(1.0*mask+1.0);
const D res= (mmax-mmin)*(de)+mmin;
// the IdxTy was getting fudding with -n being huge positive fudd 
//MM_MSG(" ASS FUDD "<<MMPR4(de,mmax,mmin,res));
return res; 
}

void test_plot( CommandInterpretter & li )
{
D beta=.01; D xpc=0; D ypc=0; D xc=0;
D yc=0; D vxt=0; D vyt=0; D a=1; D b=1;
const IdxTy quad=1; const IdxTy quadp=1;
const D rn=30;
for (IdxTy bi=0; bi<10; ++bi)
{
beta=1000;
for (IdxTy ib=0; ib<bi; ++ib) beta=beta/10;
for (IdxTy i=0; i<1000; ++i)
{
xc=xrand(rn,-rn);
xpc=vxt+ xrand(rn,-rn);
yc= xrand(rn,-rn);
ypc=vyt + xrand(rn,-rn);
const D rad=(xc+vxt-xpc)*(xc+vxt-xpc)+(yc+vyt-ypc)*(yc+vyt-ypc);
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
//D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad,false) ;
bool olf=!false;
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta,1 , 1,olf||false) ;
D res2= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 2, 2,olf||false) ;
D res3= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 3, 3,olf||false) ;
D res4= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 4, 4,olf||false) ;
MM_MSG(MMPR4(i,(xc-xpc),(yc-ypc),res)<<MMPR4(res2,res3,res4,rad))
//MM_MSG(MMPR(i)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
} // i 
} // bi 
} // test_plot

void test_fl( CommandInterpretter & li )
{

D a=1;
D b=1;
D xpc=0;
D ypc=0;
//D xc=a;
D yc=0;
D vxt=0;
D vyt=0;
D beta=li.wordf(1);
D dest=D(li.wordf(2))/2.0;
D xpu=1+dest;
D xpl=0+dest;
D xu=1-dest;
D xl=0-dest;
D vt=0;
bool qpneg=!false;
bool qneg=!false; // change for test
// mb_fl_diff_poly<0  ress=-9.16713e-12 res=-9.16713e-12 iz=0 beta=0.01 _xpu=-23.0081 _xpl=-24.0081 _xu=27.8284 _xl=26.8284 vt=0 qpneg=1 qneg=1 base=-1.91013
 beta=0.01; xpu=-23.0081; xpl=-24.0081;xu=27.8284;xl=26.8284;//vt=0 qpneg=1 qneg=1 base=-1.91013

//IdxTy pattern=0;
//StrTy x=m_flp.get("pattern");
//if (x.length()!=0) pattern=::atoi(x.c_str());
//StrTy x=m_flp.get("xc");
//if (x.length()!=0) xc=::atof(x.c_str());
//MM_MSG(MMPR4(x,int(x.c_str()[0]),xc,pattern))

// fl_mb_diff does not appear to work yet 
// D res=fl_mb_diff( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);
 D reses=esimple( beta, xpu, xpl, xu, xl ,vt, qpneg, qneg);
 D respoly=fl_mb_diff_poly( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);
 D rescomp=fl_mb_diff_comp( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);
 D rescomp2=fl_mb_diff_comp2( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);
 D reslanczos=fl_mb_diff_lanczos( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);

//D res=fl_mb_diff(beta,xpc,ypc,xc,yc,vxt,vyt,a,b,beta,pattern,noverlap);
//MM_MSG(MMPR2(beta,dest)<<MMPR4((respoly-reses),respoly,res,reses))
MM_MSG(MMPR4(beta,dest,rescomp,rescomp2)<<MMPR3((respoly-reses),respoly,reses))
MM_MSG(MMPR3(beta,dest,(rescomp-respoly))<<MMPR3((reslanczos-reses),reslanczos,reses))

} // test_fl
void test_mappoly( CommandInterpretter & li )
{
D beta=.01; D xpc=0; D ypc=0; D xc=0;
D yc=0; D vxt=0; D vyt=0; D a=1; D b=1;
const IdxTy quad=1; const IdxTy quadp=1;
const D rn=30;
for (IdxTy i=0; i<10000; ++i)
{
xc=xrand(rn,-rn);
xpc=vxt+ xrand(rn,-rn);
yc= xrand(rn,-rn);
ypc=vyt + xrand(rn,-rn);
const D rad=(xc+vxt-xpc)*(xc+vxt-xpc)+(yc+vyt-ypc)*(yc+vyt-ypc);
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
//D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad,false) ;
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta,1 , 1,false) ;
D res2= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 2, 2,false) ;
D res3= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 3, 3,false) ;
D res4= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 4, 4,false) ;
//D respoly=fl_mb_diff_poly( beta, xpu, xpl, xu, xl, vt, qpneg, qneg);
MM_MSG(MMPR4(i,(xc-xpc),(yc-ypc),res)<<MMPR4(res2,res3,res4,rad))
//MM_MSG(MMPR(i)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
} // i 


} // test_mappoly



void test_diff( CommandInterpretter & li )
{

D a=1;
D b=1;
D xpc=0;
D ypc=0;
D xc=a;
D yc=0;
D vxt=0;
D vyt=0;
D beta=li.wordf(1);
IdxTy pattern=0;
StrTy x=m_flp.get("pattern");
if (x.length()!=0) pattern=::atoi(x.c_str());
x=m_flp.get("xc");
if (x.length()!=0) xc=::atof(x.c_str());
MM_MSG(MMPR4(x,int(x.c_str()[0]),xc,pattern))


bool noverlap=false;
D res=evaluate_node_pair(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,pattern,noverlap);
MM_MSG(MMPR(res))

} // test_diff



////////////////////////////////////////////////////

void test_beta( CommandInterpretter & li )
{
//D beta=.01; D xpc=.1; D ypc=.20; D xc=xpc+2;
//D yc=xpc+2; D vxt=0; D vyt=0; D a=.1; D b=.1;

D beta=.01; D xpc=20; D ypc=20; D xc=xpc;
D yc=xpc; D vxt=0; D vyt=0; D a=.1; D b=.1;


const IdxTy quad=1; const IdxTy quadp=1;
//const D rn=30;
for (IdxTy i=0; i<20000; ++i)
{
//const D beta=log(1.0+1.0*(1+i)/10000*300*3*100*100000);
const D beta=log(1.0+0.25*(1+i)*1e-5);
const D rad=(xc+vxt-xpc)*(xc+vxt-xpc)+(yc+vyt-ypc)*(yc+vyt-ypc);
const D f=beta*rad;
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
//D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad,false) ;
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta,1 , 1,false) ;
D res2=0; //  evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 2, 2,false) ;
D res3=0; //  evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 3, 3,false) ;
D res4= 0; // evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 4, 4,false) ;
//MM_MSG(MMPR4(i,(xc-xpc),(yc-ypc),res)<<MMPR4(res2,res3,res4,rad)<<MMPR2(beta,f))
MM_MSG(MMPR4(i,(xc-xpc),(yc-ypc),res)<<MMPR(rad)<<MMPR2(beta,f))
//MM_MSG(MMPR(i)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
} // i 

} // test_beta



//////////////////////////////////////////////////

void test_origin( CommandInterpretter & li )
{
D beta=.01; D xpc=0; D ypc=0; D xc=0;
D yc=0; D vxt=-5; D vyt=0; D a=1; D b=1;
const IdxTy quad=1; const IdxTy quadp=1;

::srand(123);
D reslast=0;
// This is an ASSFUDD it fist converts -20 to an usngiedn fudding int fudd 
//const IdxTy rn=20;
const D rn=20;
IdxTy fails=0;
// this passes but only by brute force setting of the origins. 
for (IdxTy i=0; i<10000; ++i)
{
xc=xrand(rn,-rn);
xpc= xrand(rn,-rn);
yc= xrand(rn,-rn);
ypc= xrand(rn,-rn);
const D rad=(xc+vxt-xpc)*(xc+vxt-xpc)+(yc+vyt-ypc)*(yc+vyt-ypc);
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad,false) ;
MM_MSG(MMPR(i)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
if (false) if (res!=reslast)
{
MM_MSG(MMPR2(i,(res-reslast))<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
if (i!=0) { ++fails; }
}

reslast=res;
//vxt+=.1;

}
MM_MSG("test 1 " << MMPR2((reslast),fails)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR2( quadp, quad)) 

fails=0;
reslast=0;

for (IdxTy i=0; i<000; ++i)
{
vxt= xrand(rn,-rn);
xc= xrand(rn,-rn);
xpc= xc+vxt; // xrand(20,-20);
yc= xrand(rn,-rn);
vyt= xrand(rn,-rn);
ypc=4+yc+vyt; // xrand(20,-20);
const D rad=(xc+vxt-xpc)*(xc+vxt-xpc)+(yc+vyt-ypc)*(yc+vyt-ypc);
//D res=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1,  s0,  c1,  c0);
// zero is probably ok but the new base is 1 to give 1-4 wtf 
D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, 1, 1,false) ;
D res3= evaluateq( xpc+1, ypc+1, xc, yc, vxt, vyt, a, b,beta, 3, 1,false) ;
D res3p= evaluateq( xpc, ypc, xc+1, yc+1, vxt, vyt, a, b,beta, 1, 3,false) ;
//MM_MSG(MMPR(i)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
if (res!=reslast)
{
MM_MSG(MMPR(i)<<MMPR4(res,(res3-res3p),(res-res3),(res-reslast))<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR4(rad,res, quadp, quad)) 
if (i!=0) { ++fails; }
}
reslast=res;
} // i 
MM_MSG("test 2 "<< MMPR2((reslast),fails)<<MMPR4( xpc, ypc, xc, yc)<<MMPR4( vxt, vyt, a, beta)<<MMPR2( quadp, quad)) 



}



// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
if (cmd=="int") { test_int(li);} //   continue; } 
else if (cmd=="quad") { test_quad(li);} //   continue; } 
else if (cmd=="origin") { test_origin(li);} //   continue; } 
else if (cmd=="plot") { test_plot(li);} //   continue; } 
else if (cmd=="beta") { test_beta(li);} //   continue; } 
else if (cmd=="diff") { test_diff(li);} //   continue; } 
else if (cmd=="fl") { test_fl(li);} //   continue; } 
else if (cmd=="mappoly") { test_mappoly(li);} //   continue; } 
//void test_gr_integrator( CommandInterpretter & li )
else { MM_ERR(" unrecignized TEST command "<<li.dump()) } 

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

////////////////////////////////////////////////
// fopiced from mjm_qfle.h 2017-07-07
// check math and conventions on decary factor etc 
D K2K1K0Ke2Ke1Ke0(const D& yf, const D& yi,const D& a
, const D & c3, const D& c2, const D & c1, const D& c0
, const D& ce2, const D & ce1, const D& ce0) const
{

static D ps=::sqrt(M_PI);
static D third=1.0/3.0;
const D azi=a*yi;
const D azi2=azi*azi;
const D azf=a*yf;
const D azf2=azf*azf;

const D cerfi=erf(azi);
const D cerff=erf(azf);
const D cexpi=exp(-azi2); // /(ps*a);
const D cexpf=exp(-azf2); // /(ps*a);
const D expd=1.0/(ps*a);



const D k0=c0*( yf*cerff+cexpf*expd-cexpi*expd-yi*cerfi);
const D k1=c1*(.5*yf*yf*cerff-.5*yi*yi*cerfi- .25/a/a*(cerff-cerfi)+.5*expd*(yf*cexpf-yi*cexpi));
const D k2=c2*third*(((1.0+azf2)*cexpf*expd-(1.0+azi2)*cexpi*expd) /(a*a)
        +yf*yf*yf*cerff-yi*yi*yi*cerfi);

D k3=0;
const D  a3=a*a*a;
const D  a4=a3*a;
const D z3f=yf*yf*yf;
const D z3i=yi*yi*yi;
const D z4f=z3f*yf;
const D z4i=z3i*yi;
// http://www.wolframalpha.com/input/?i=integrate+x^3+erf%28a*x%29
k3= c3*(4.0*a3*(z3f*cexpf-z3i*cexpi)+6.0*(cexpf*azf-cexpi*azi)
    +ps*a4*4.0*(cerff*z4f-cerfi*z4i)-3.0*ps*(cerff-cerfi))/(16*ps*a4);

const D erfn=.5*ps;
const D ke2=ce2*(-.5/(a*a)*(yf*cexpf-(cerff-cerfi)*erfn/a-yi*cexpi));   // Ke2
const D ke1=ce1*( -.5/(a*a)*(cexpf-cexpi)); // Ke1 
const D ke0=ce0*(((cerff-cerfi)*erfn/a)); // Ke0 

return k3+k2+k1+k0+ke2+ke1+ke0;
}
// difference of two erf's and two normals
// erf(xf)-erf(xi) and exp(xf2)-exp(xi^2)
int doen( D & derf, D & dnorm, const D & xf, const D & xi)
{

return 0; 

}

// evaluate the integral over a line from xpl to xpu away from a source at x
// with a collecting function Q(x') as below 
// the src integral is just exp(-beta*(x')^2)
static  D const_line_source( const D & beta, const D & _xpu,const D & _xpl, const D & _x 
, const D & vt, const bool qpneg, const bool limneg)
{

//mi.doen_split(derf,dnorm, xf,xi,nn,true);
//mi.doen_split(derf,dnorm, xf,xi,0,false);
	static mjm_integrals mi;
	const bool simple_erf=!false;
	D dnorm=0;
	IdxTy branch=0;
	const bool dump_coefs=!true;
	const bool dump_js=!true;
	const bool dump_limits=!true;
	const bool check_small_bx=!true; // check the small value expression for comparison to "exact"
	const bool canonize=true; // remove offsets from coordinates 
	bool small_ok=true; // ok to use expansion for small parameters 
	static const  D ps=::sqrt(M_PI);
	D base=0;
	if (canonize)
	{	// this is not right but seemed to work ??? TODO FIXME 
		//base=-(_xpu+_xpl+_x+vt)/3.0;
//		base=-(_xpu+_xpl+_x+vt)/3.0;
		base=-(_xpu+_xpl)*.5;
		MM_ONCE(" canonicalization code needs to be sign checked",)
	}
	const D xpu=_xpu+base;
	const D xpl=_xpl+base;
	const D x=_x+base;
	const D dxp=xpu-xpl;
	const D q1=(qpneg?1:(-1.0))/dxp;
	const D q0=(qpneg?(-xpl):xpu)/dxp;
	if (beta==0)
	{
		D res=xpu-xpl;
		res=(.5*q1*res+q0)*res;
		branch=1;
		MM_MSG(" beta is zero, danger will robinson "<<MMPR4(xpu,xpl,x,res))
		return res;
	}
	const D bs=::sqrt(beta);  
	const D bsyf=bs*(xpu-x-vt);
	const D bsyi=bs*(xpl-x-vt);
	// see comments for exp(-26) likely 1+delta problem 
	static const D maxexp=5.5-.5; // sqrt(100);
	static const D minexp=-maxexp;

// this logic is wrong doh
	if ( bsyf>maxexp) if (bsyi>maxexp) return 0;
	if ( bsyf<minexp) if (bsyi<minexp) return 0;

	const D small_lim=2e-2;
	if (small_ok) {small_ok&=(bsyf<small_lim)&&(bsyf>-small_lim);
	if (small_ok){ small_ok&=(bsyi<small_lim)&&(bsyi>-small_lim);}}
	D ressbx=0;
	const D b2xv=2.0*beta*(x+vt);
	// in the beta -> zero limit this looks more stable the other one starts to oscillate
	if (check_small_bx||small_ok)
	{ // taylor expand the exp and do polynomial integration
		const D c1=1.0-beta*(x+vt)*(x+vt);
		// this needs to handle the overlap conditions 
		D xpue=xpu;
		D xple=xpl;
		D xcu=x+vt;
		// may work except for partial overlap 
		if (!limneg&&(xple<xcu)) { xple=xcu; }
		if (limneg&&(xpue>xcu)) { xpue=xcu; }
		if (xple!=xpl)  { ressbx-=xpl-xcu;}
		if (xpue!=xpu)  { ressbx+=xpu-xcu;}
		 bool skip = ( !limneg&&(xpue<xcu));
		 skip |= ( limneg&&(xple>xcu));

		if (!skip)
		{
		const D xpu2=xpue*xpue;
		const D xpl2=xple*xple;
		const D t3=(xpu2+xpl2)*(xpu2-xpl2)*.25;
		const D t2=(xpue*xpue*xpue-xple*xple*xple)/3.0;
		const D t1=(xpue*xpue-xple*xple)*.5;
		const D t0=(xpue-xple);
		ressbx+=t3*(-beta*q1)+t2*(-beta*q0+q1*b2xv)+t1*(q1*c1+q0*b2xv)+t0*q0*c1;
		}	
		//if (dump_js) 
			MM_MSG(" small "<<MMPR4(ressbx,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR(bsyi))
		branch=2;
		if (small_ok) 
		{
			if (ressbx<0) { MM_MSG(" res<0 "<<MMPR4(ressbx,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR2(bsyi,branch))}
			return ressbx;
		}
	}
	D res=0;
	// for negative limit, just check lowest
	if (!limneg&&(bsyi<0))
	{
		if (bsyf>0) { 	
	if (simple_erf) { 	res+=.5*ps*q0/bs*(erf(bsyf)-erf(0)); } 
else { 		mi.doen_split(res,dnorm, bsyf,0,0,false); res*=.5*ps*q0/bs; }
		res-=q1*.5/beta*(exp(-bsyf*bsyf)-(1.0)); // doh sign 
		res-=(xpl-x-vt);
		} else res+=xpu-xpl;
		if (dump_js) MM_MSG(" limneg "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR(bsyi))
		branch=3;
	}
	else if ((limneg)&&(bsyf>0)) // positive limit
	{
		if (bsyi<0)
		{
		if (simple_erf) { res+=.5*ps*q0/bs*(erf(0)-erf(bsyi)); } 
		else { mi.doen_split(res,dnorm, 0,bsyi,0,false); res*=.5*ps*q0/bs; } 
		res-=q1*.5/beta*(1.0-exp(-bsyi*bsyi)); // doh sign // dreaded 1-episilon 
		res+=xpu-x-vt;
		} else res+=xpu-xpl;
		if (dump_js) MM_MSG(" limpos "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyi)<<MMPR(bsyf))
		branch=4;
	}	
	else 
	{
		// this is just a simple integral of (q1*xp+q0)*exp(-beta*(xp-x-vt)^2)
		if ( simple_erf) { res+=.5*ps*q0/bs*(erf(bsyf)-erf(bsyi));}
		else { mi.doen_split(res,dnorm, bsyf,bsyi,0,false); res*=.5*ps*q0/bs; } 
		//res+=q1*.5/beta*(exp(-bsyf*bsyf)-exp(-bsyi*bsyi));
		res-=q1*.5/beta*(exp(-bsyf*bsyf)-exp(-bsyi*bsyi)); // doh sign 
		if (dump_js) MM_MSG(" res "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyi)<<MMPR(bsyf))
		branch=5;
	}
if (check_small_bx) { MM_MSG("taylor vs integral "<< MMPR4(res,ressbx,beta,bsyf)<<MMPR(bsyi)) } 
	if (res<0) { 
//MM_MSG(" res<0 "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR4(bsyi,q1,q0,branch))
if (res<-1e-4) MM_MSG(" res<0 "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR4(bsyi,q1,q0,branch))
MM_ONCE(" res<0 rest suppressed "<<MMPR4(res,beta,x,vt)<<MMPR4(xpu,xpl,base,bsyf)<<MMPR4(bsyi,q1,q0,branch),)
MM_ONCE(" setting res=0",)
res=0;
}
	return res; 
} // const_line_source

// fl_mb_diff class for sf 
class fl_sf_auto
{
public:
fl_sf_auto( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg, const D & base)
{
	xpu=_xpu+base;
	xpl=_xpl+base;
	xu=_xu+base;
	xl=_xl+base;
	const D dx=xu-xl;
	const D dxp=xpu-xpl;
	q1=(qpneg?1:(-1.0))/dxp;
	q0=(qpneg?(-xpl):xpu)/dxp;
	p1=(qneg?1:(-1.0))/dx;
	p0=(qneg?(-xl):xu)/dx;
	bs=::sqrt(beta);  
	bsyff=bs*(xpu-xu-vt);
	// changed convention, now first idx is for x not xp 
	bsyfi=bs*(xpl-xu-vt);
	bsyif=bs*(xpu-xl-vt);
	bsyii=bs*(xpl-xl-vt);
	
	bsyff2=bsyff*bsyff;
	bsyfi2=bsyfi*bsyfi;
	bsyii2=bsyii*bsyii;
	bsyif2=bsyif*bsyif;


	overu=(xu<xpu)?xu:xpu;
	overl=(xl>xpl)?xl:xpl;
	//if (dump_limits ) { MM_MSG(MMPR4(xpu,xpl,xu,xl)<<MMPR4(beta,vt,qpneg,qneg)<<MMPR4(q1,q0,p1,p0))}

}
bool olap() const { return overu>overl; }

void msg() const
{
MM_MSG(MMPR4( xpu,xpl,xu,xl)<<MMPR4(q1,q0,p1,p0)<<MMPR2(overu,overl))
MM_MSG(MMPR4(bs,bsyff,bsyfi,bsyif)<<MMPR(bsyii))
MM_MSG(MMPR4(bsyff2,bsyfi2,bsyif2,bsyii2))

}

D xpu,xpl,xu,xl,q1,q0,p1,p0,overu,overl;
D bs,bsyff,bsyfi,bsyif,bsyii;
D bsyff2,bsyfi2,bsyif2,bsyii2;

} ; // fl_sf_auto
// see comments for esimple, do the normalized diference for flow
// note that the x-x' does not contain a vt term as this is the zero time
// case so only the difference is noted. 
// \int dx' Q(x') \int dx P(x)(  exp(-beta*(x'-x-vt)^2)*A - \delta(x-x'))
// http://www.wolframalpha.com/input/?i=integrate+y*erf%28sqrt%28b%29*y%29 
static  D fl_mb_diff( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	static mjm_integrals mi;
	D erff,normf,erfi,normi;
	// this should work for all ranges 
	// TODO IIRC this is the wrong integral again, getting 
	// mixed up between double integrals 
	mi.doen_split(erff,normf, sf.bsyff,sf.bsyfi,0,false); 
	mi.doen_split(erfi,normi, sf.bsyif,sf.bsyii,0,false); 
	// the origi effects the shape functions NOT the exponents, 
	const D qp2=sf.q1*sf.p1;	
	const D qp1=sf.q0*sf.p1+sf.q1*sf.p0;
	const D qp0=sf.q0*sf.p0;
	const D bs=sqrt(beta);
	const D A=bs/ps;
	const D cxf=sf.q0+sf.q1*sf.xu+sf.q1*vt;
	const D cxi=sf.q0+sf.q1*sf.xl+sf.q1*vt;
	const D c0=sf.p0-sf.p1*vt;
	const D M=.5/beta;
	const D N=.5*ps/bs;
	const D t2=qp2*M*M*(erff-erfi)*A;
	const D t1=sf.p1*N*M*(cxf*erff-cxi*erfi)*A;
	const D t=A*c0*(cxf*(N*erff+M*normf)-cxi*(N*erfi+M*normi));
	// the use of C0 is not right, need to ao algebra on other terms 
	// unstable, this could be large - large	
	// miscellaneous terms of form y^n*e(y)
	const D fexf=exp(-sf.bsyff2)+exp(-sf.bsyfi2);
	const D fexi=exp(-sf.bsyif2)+exp(-sf.bsyii2);
	const D ferf=erf(sf.bsyff)+erf(sf.bsyfi);
	const D feri=erf(sf.bsyif)+erf(sf.bsyii);

	const D dyf=sf.bsyff-sf.bsyfi;
	const D dxx=sf.xu-sf.xl;
	const D dxx2=sf.xu*sf.xu-sf.xl*sf.xl;
	const D syf=sf.bsyff+sf.bsyfi;
	const D sxx=sf.xu+sf.xl;
	const D sxx2=sf.xu*sf.xu+sf.xl*sf.xl;
	const D dyi=sf.bsyif-sf.bsyii;
	const D syi=sf.bsyif+sf.bsyii;
	const D syf2=sf.bsyff2+sf.bsyfi2;
	const D dyf2=sf.bsyff2-sf.bsyfi2;
	const D dyi2=sf.bsyif2-sf.bsyii2;
	const D syi2=sf.bsyif2+sf.bsyii2;
	// this is just x^n*e(x) - y^n*e(y) reorganized 
	const D derff=.5*(dyf*ferf+syf*erff);
	const D derfi=.5*(dyi*feri+syi*erfi);
	const D derff2=.5*(dyf2*ferf+syf2*erff);
	const D derfi2=.5*(dyi2*feri+syi2*erfi);
	const D dexpf=.5*(dyf*fexf+syf*normf);
	const D dexpi=.5*(dyi*fexi+syi*normi);

	const D I01f=A*N*(.5*derff2-.5*M*dexpf-.5*M*erff);
	const D I01i=A*N*(.5*derfi2-.5*M*dexpi-.5*M*erfi);
	
	const D I00f=(N*derff+M*normf)*A;
	const D I00i=(N*derfi+M*normi)*A;

MM_MSG(MMPR4(I01f,I01i,derff2,derfi2)<<MMPR4(dexpf,dexpi,erff,erfi))
	const D t0=-c0*sf.q1*(I01f-I01i);

	const D tkluge1=-c0*(cxf*I00f-cxi*I00i); // -c0*sf.q1*(I01f-I01i); 
MM_MSG(MMPR4(cxf,I00f,cxi,I00i))

	if (sf.olap())
	{
		const D y1=(sf.overu-sf.overl);
		const D y2=(sf.overu*sf.overu-sf.overl*sf.overl);
		const D y3=(sf.overu*sf.overu*sf.overu-sf.overl*sf.overl*sf.overl);
		iz=y3*qp2/3.0 + .5*y2*qp1+y1*qp0;	
	}
	res=-iz+t2+t1+t0+t+tkluge1;
 	MM_MSG(MMPR4( beta,_xpu,_xpl,_xu)<<MMPR3(_xl,vt,base))
	sf.msg();
	MM_MSG(MMPR4(res,iz,t2,t1)<<MMPR3(t0,t,tkluge1))
	return res; 
} //fl_mb_diff


/*

http://www.wolframalpha.com/input/?i=integrate+%28x^0%29+exp%28-%28b%29x^2%29dx
integral x^0 exp(-b x^2) dx = (sqrt(π) erf(sqrt(b) x))/(2 sqrt(b)) + constant

http://www.wolframalpha.com/input/?i=integrate+%28x^1%29+exp%28-%28b%29x^2%29dx
integral x^1 exp(-b x^2) dx = -e^(-b x^2)/(2 b) + constant

http://www.wolframalpha.com/input/?i=integrate+%28x^2%29+exp%28-%28b%29x^2%29dx

integral x^2 exp(-b x^2) dx = (sqrt(π) erf(sqrt(b) x))/(4 b^(3/2)) - (x e^(-b x^2))/(2 b) + constant


http://www.wolframalpha.com/input/?i=integrate+%28x^3%29+exp%28-%28b%29x^2%29dx
integral x^3 exp(-b x^2) dx = -(e^(-b x^2) (b x^2 + 1))/(2 b^2) + constant

http://www.wolframalpha.com/input/?i=integrate+%28x^4%29+exp%28-%28b%29x^2%29dx
integral x^4 exp(-b x^2) dx = (3 sqrt(π) erf(sqrt(b) x))/(8 b^(5/2)) - (x e^(-b x^2) (2 b x^2 + 3))/(4 b^2) + constant


http://www.wolframalpha.com/input/?i=integrate+%28x^5%29+exp%28-%28b%29x^2%29dx

integral x^5 exp(-b x^2) dx = -(e^(-b x^2) (b^2 x^4 + 2 b x^2 + 2))/(2 b^3) + constant

// http://www.wolframalpha.com/input/?i=integrate+y*erf%28sqrt%28b%29*y%29 

*/

// see comments for fl_mb_diff  impleented using general poly 
// math for finding numercal problems
// , do the normalized diference for flow
// note that the x-x' does not contain a vt term as this is the zero time
// case so only the difference is noted. 
// MAY BE NORMALIZED \int dx' Q(x') \int dx P(x)(  exp(-beta*(x'-x-vt)^2)*A - \delta(x-x'))
// NOTE check usage of const D A=bs/ps; // 1.0; // ps/bs;
// with A==1 or disabled and iz or overlap OFF this should 
// match esimple. 
// TODO FIXME the integrators add high order terms with zero coeffs wtf 
static  D fl_mb_diff_poly( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=!false;
	const bool debug_step_2=!false;
	const bool debug_step_3=!false;
	const bool debug_pieces=!false;
	const bool debug_ltz=!false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);
	const D bsyfi=bs*yfi; // (sf.xpl-sf.xu-vt);
	const D bsyii=bs*yii; // (sf.xpl-sf.xl-vt);

if (debug_step_1){ 
	MM_MSG(MMPR4(sf.xu,sf.xl,sf.xpu,sf.xpl)<<MMPR(vt))
	MM_MSG(MMPR4(bsyff,bsyfi,bsyif,bsyii))
}

	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	//typedef shape_function_integrals::simple_polynomial<D> Po;
	typedef Pol::simple_polynomial<D> Po;
	typedef Pol::dense_2_polynomial<D> D2;
	Po P,Q,pxferf,pxierf,pxfexp,pxiexp,pulim,pllim;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	if (debug_step_1) { MM_MSG(MMPR2(P.to_string(),Q.to_string())) }
	D2 yint(P,1.0,-1.0,-vt);
	if (debug_step_1) { MM_MSG(MMPR(yint.to_string())) } 
	// integrate multiplied by exp(-b*y*y)
	const IdxTy erfsz=P.size()+1;
	const IdxTy expsz=P.size();

	D2 yerf(erfsz),yexp(expsz);
	// beta,factor, dir 
	//D2::integrate_e(yexp,yerf,yint,beta,1,0);
	// TODO FIXME use _new here 
	D2::integrate_e_old(yexp,yerf,yint,beta,1,0);
	if (debug_step_1) { MM_MSG(MMPR(yexp.to_string())) } 
	if (debug_step_1) { MM_MSG(MMPR(yerf.to_string())) } 
	pulim.push_back(-sf.xu-vt); pulim.push_back(1);
	pllim.push_back(-sf.xl-vt); pllim.push_back(1);
// evaluate( Td & dest, const Myt & x, const value_type & v, const IdxTy dir)
	// this is not right, these need to be evaluated at y which creates
	// higher power terms in x'. These need to be evaluaated at a polynomail  
//	D2::evaluate( pxferf, yerf, sf.xu, 1);
	const IdxTy edir=1;
	D2::evaluate( pxferf, yerf, pulim, edir);
	D2::evaluate( pxierf, yerf, pllim, edir);
	D2::evaluate( pxfexp, yexp, pulim, edir);
	D2::evaluate( pxiexp, yexp, pllim, edir);
	if ( debug_step_2) {
	MM_MSG(MMPR(pxferf.to_string()))
	MM_MSG(MMPR(pxierf.to_string()))
	MM_MSG(MMPR(pxfexp.to_string()))
	MM_MSG(MMPR(pxiexp.to_string()))
	 } // debug_step_2
	// convert y to x'
	Po T,yfp,yip,pip,pfp;  

	pfp.push_back((vt+sf.xu));
	pip.push_back((vt+sf.xl));
	pfp.push_back(1);
	pip.push_back(1);

	// these are now polynomials after being evaluated at x.
	// multiply by the Q polynomial, do the shifts
	// and integrate x'
	Ss ss; 
	D rf=0, ri=0;
	const IdxTy erffac=0;
	const IdxTy expfac=1;
	//D res=0;
	// this is now a polynomial in x' NOT y
	mi.multiply_polynomials(T,pxferf,Q);
	if (debug_step_3) { MM_MSG("pxferf "<<MMPR(T.to_string())) }
	piece( rf, ri, mi, T,beta, yff,yfi,  bsyff, bsyfi,erffac,pfp,debug_pieces);
	res+=rf-ri;
	if (debug_step_3) { ss<<" p1 "<<MMPR4(res,rf,ri,(rf-ri)); } 
	mi.multiply_polynomials(T,pxfexp,Q);
	if (debug_step_3) { MM_MSG("pxfexp "<<MMPR(T.to_string()))}
	piece( rf, ri, mi, T,beta,  yff,yfi,  bsyff, bsyfi,expfac,pfp,debug_pieces);
	res+=rf-ri;
	if (debug_step_3) { ss<<" p2 "<<MMPR4(res,rf,ri,(rf-ri)); }
	mi.multiply_polynomials(T,pxierf,Q);
	if (debug_step_3) { MM_MSG("pxierf "<<MMPR(T.to_string())) } 
	piece( rf, ri, mi, T,beta,  yif,yii,  bsyif, bsyii,erffac,pip,debug_pieces);
	res-=rf-ri;
	if (debug_step_3) { ss<<" p3 "<<MMPR4(res,rf,ri,(rf-ri)); } 
	mi.multiply_polynomials(T,pxiexp,Q);
	if (debug_step_3) { MM_MSG("pxiexp "<<MMPR(T.to_string())) } 
	piece( rf, ri, mi, T,beta,  yif,yii,  bsyif, bsyii,expfac,pip,debug_pieces);
	res-=rf-ri;
	if (debug_step_3) { ss<<" p4 "<<MMPR4(res,rf,ri,(rf-ri)); } 

	// now evaluate at xf and xi, then multiply the polynomials
	// by Q and integreate the exp and erf polynomials
	// the x polynomial is now in y=x'-x-vt which is then rearranged
	// in terms of x' for the second integral
	// return d = g(f(x))
	// mi.composite_polynomial(d,g,f);
	if (sf.olap())
	{
		Po R,S;
		mi.multiply_polynomials(R,P,Q);
		mi.integrate_polynomial(S,R);
		// efficiency wtf 
		iz=mi.evaluate_polynomial(S,sf.overu)
		 -mi.evaluate_polynomial(S,sf.overl);
	}
	// the whole thing ignored a -dy in the change of variable 
	res=-res;
	const D ress=res;
	//res=A*res-0*iz; // +t2+t1+t0+t+tkluge1;
//	res=1.0*res-0*iz; // +t2+t1+t0+t+tkluge1;
	res=1.0*res-1.0*iz; // +t2+t1+t0+t+tkluge1;
	if (debug_step_3) { ss<<MMPR(iz); MM_MSG(ss.str()) } 
	if (debug_ltz) if (ress<0)
	{
		MM_MSG(" mb_fl_diff_poly<0 "<<MMPR4(ress,res,iz,beta)<<MMPR4( _xpu, _xpl,  _xu, _xl)<<MMPR4( vt, qpneg, qneg,base))

	}
	return res;
} // fl_mb_diff_poly


//template <class Tm, class Sf,class Tp > 
template <class Tm, class Tp > 
static		void piece( D & rf, D & ri, Tm & mi, const Tp & tin  
	,const D & beta, const D & yf, const D & yi,  const D & bsyf, const D & bsyi, 
	const IdxTy factor,const Tp & pf, const bool debug_pieces)
{
	Tp T,Texp,Terf;
	// TODO make this lol 
// void integrate_polynomial_e(derf, dexp, v1, b, factor, init_dest=true )
// right now the polynomial is in x' ( aka x ) but the factor
// is in y evaluted at xf or xi 
	mi.composite_polynomial(T,tin,pf);
	mi.integrate_polynomial_e(Terf,Texp,T,beta,factor);
	if (debug_pieces) {MM_MSG(" "<<MMPR2(T.size(),T.to_string(0,StrTy("y"))))
	MM_MSG(" "<<MMPR2(Terf.size(),Terf.to_string(0,StrTy("y"))))
	MM_MSG(" "<<MMPR2(Texp.size(),Texp.to_string(0,StrTy("y"))))
}
	// So now the polynomials are functions of y 
	rf=mi.evaluate_polynomial(Texp,yf)*exp(-bsyf*bsyf)
			+mi.evaluate_polynomial(Terf,yf)*erf(bsyf);
	ri=mi.evaluate_polynomial(Texp,yi)*exp(-bsyi*bsyi)
			+mi.evaluate_polynomial(Terf,yi)*erf(bsyi);
}
/////////////////////////////////////////////////////////
static D aerfc(const D & x) 
{
const bool pos=(x>0);
return (pos)?erfc(x):-erfc(-x);
}
// add these two togeterh to get result 
static void aerfc(D & i, D & v, const D & x) 
{
if (x>.5) { v=-erfc(x); i=1; return ; } 
if (x<-.5) { v=erfc(-x); i=-1; return ; } 
i=0;
v=erf(x);
}


static  D fl_mb_diff_comp( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=false;
	const bool debug_step_2=false;
	const bool debug_step_2b=false;
	const bool debug_step_3=!false;
	const bool debug_step_3r=!false;
	const bool debug_resvar=!false;
//	const bool debug_pieces=!false;
	const bool debug_ltz=!false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);
	const D bsyfi=bs*yfi; // (sf.xpl-sf.xu-vt);
	const D bsyii=bs*yii; // (sf.xpl-sf.xl-vt);

if (debug_step_1){ 
	MM_MSG(MMPR4(sf.xu,sf.xl,sf.xpu,sf.xpl)<<MMPR(vt))
	MM_MSG(MMPR4(bsyff,bsyfi,bsyif,bsyii))
}

	Ss ss; 
	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	//typedef shape_function_integrals::simple_polynomial<D> Po;
	typedef Pol::simple_polynomial<D> Po;
	typedef Pol::dense_2_polynomial<D> D2;
	Po P,Q,pxferf,pxierf,pxfexp,pxiexp,pulim,pllim;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	if (debug_step_1) { MM_MSG(MMPR2(P.to_string(),Q.to_string())) }
	// TODO FIXME check the integration variable has to be sqrt(b)y NOT y
	// this only works for beta==1 unless care with erf expansion later 
	D2 yint(P,-1.0,1.0,-vt);
	D2 qint(Q,1.0,1.0,vt);
	D2 yy(1.0,1.0,vt);
	D2 yyx(1.0,-1.0,-vt);
	if (debug_step_1) { MM_MSG(MMPR(yint.to_string(0,"y","x'"))) } 
	if (debug_step_1) { MM_MSG(MMPR(yy.to_string(0,"y","x'"))) } 
	// integrate multiplied by exp(-b*y*y)
	const IdxTy erfsz=P.size()+1;
	const IdxTy expsz=P.size();

	D2 yerf(erfsz),yexp(expsz),yexp2,yerf2,yexp3,yerf3;
	// beta,factor, dir 
//integrate_e( dexp,derf,  x,  b,factor,  dir)
	const bool use_rerf=!false;
	if (use_rerf) D2::integrate_e_rnew(yexp,yerf,yint,beta,1,0);
	else D2::integrate_e_new(yexp,yerf,yint,beta,1,0);
	if (debug_step_2) { MM_MSG(" the exp and erf parts of x integral" ) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp.to_string(0,"y","x'"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf.to_string(0,"y","x'"))) } 
	// now change x' to x
	D2::composite(yexp2,yexp,yy,1);
	D2::composite(yerf2,yerf,yy,1);
	if (debug_step_2) { MM_MSG(" x integral as function of x and y " ) } 
	if (debug_step_2) { MM_MSG(MMPR(yy.to_string(0,"x","y"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	D2 temp,temp2;
	// do not evaluate yet, change variables and evaluate at end
	temp=qint.transpose();
	if (debug_step_2) { MM_MSG(MMPR(qint.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(temp.to_string(0,"x","y"))) } 
qint=temp;
	// the variable "x' " in yexp and yerf needs to be chagend to y
	D2::multiply(yexp,yexp2,qint);
	D2::multiply(yerf,yerf2,qint);
	if (debug_step_2b) { MM_MSG(" final integrand as f of y and x " ) } 
	if (debug_step_2b) { MM_MSG(MMPR(yexp.to_string(0,"y","x"))) } 
	if (debug_step_2b) { MM_MSG(MMPR(yerf.to_string(0,"y","x"))) } 
	

/*	D2::composite(yexp2,yexp,yy,0);
	D2::composite(yerf2,yerf,yy,0);
	D2::multiply(temp,yexp2,Q,1);
	D2::multiply(temp2,yerf2,Q,1);
*/
	// now the matrix should be in x an x' again 
	// multiply by Q(x') and integrate replacing x' with y
	const IdxTy isz=yexp.get_size()+3;
	yexp2.set_size(isz);
	yerf2.set_size(isz);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yexp,beta,1,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
	if (debug_step_3) { MM_MSG("  exp integrand  parts " ) } 
	if (debug_step_3) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// add the erf terms to those generated from exp 
	// note that this may be a numerical problem will need to check
	// priot to evaluation 
	//D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yerf,beta,0,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
	if (debug_step_3r) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3r) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// covert the erf polynomial form y and x to x' and x 
	D2::composite(yerf3,yerf2,yyx,0);
	if (debug_step_3r) { MM_MSG(MMPR(yerf3.to_string(0,"(x')","x"))) } 
//	for (D aoff=-5; aoff<6; aoff+=1)
//	{
		D eff,eif,efi,eii, cff,cif,cfi,cii;
		aerfc(cff,eff,bsyff);
		aerfc(cii,eii,bsyii);
		aerfc(cfi,efi,bsyfi);
		aerfc(cif,eif,bsyif);
		// needs rerf term 
		D ssp=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,
	 		erf(bsyii),erf(bsyff), erf(bsyif),erf(bsyfi));
		D sspc=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,
	 		//aerfc(bsyii),aerfc(bsyff), aerfc(bsyif),aerfc(bsyfi));
	 		eii,eff, eif,efi);
		//D sspone=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,1.0,1.0,1.0,1.0);
		D sspone=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,cii,cff,cif,cfi);
		MM_MSG(MMPR2(ssp,sspone))
//	}
	yexp=yexp2;//.transpose();
	yerf=yerf2;//.transpose();
	const D rerff=(use_rerf)?((ps/bs*.5)):1.0;
	//const D yif=(sf.xpu-sf.xl-vt);
	const D r1=yexp.evaluate(yff,sf.xu)*exp(-beta*yff*yff);
	// the erf term is likely the problem. 
	const D r2=yerf.evaluate(yff,sf.xu)*rerff*erf(bsyff);
	const D r2alt=yerf3.evaluate(sf.xpu,sf.xu)*rerff*erf(bsyff);
	const D r2altp=yerf3.evaluate(sf.xpu,sf.xu);
	const D r3=yexp.evaluate(yii,sf.xl)*exp(-beta*yii*yii);
	const D r4=yerf.evaluate(yii,sf.xl)*rerff*erf(bsyii);
	const D r4alt=yerf3.evaluate(sf.xpl,sf.xl)*rerff*erf(bsyii);
	const D r4altp=yerf3.evaluate(sf.xpl,sf.xl);

	const D r5=yexp.evaluate(yfi,sf.xu)*exp(-beta*yfi*yfi);
	const D r6=yerf.evaluate(yfi,sf.xu)*rerff*erf(bsyfi);
	const D r6alt=yerf3.evaluate(sf.xpl,sf.xu)*rerff*erf(bsyfi);
	const D r6altp=yerf3.evaluate(sf.xpl,sf.xu);
	const D r7=yexp.evaluate(yif,sf.xl)*exp(-beta*yif*yif);
	const D r8=yerf.evaluate(yif,sf.xl)*rerff*erf(bsyif);
	const D r8alt=yerf3.evaluate(sf.xpu,sf.xl)*rerff*erf(bsyif);
	const D r8altp=yerf3.evaluate(sf.xpu,sf.xl);
	
	///D resalt= r1+r3-r5-r7+rerff*ssp;	
	//D resalt= r1+r3-r5-r7-rerff*sspc+rerff*sspone;	
	D resalt= r1+r3-r5-r7+rerff*sspc+rerff*sspone;	
	D resold= r1+r2+r3+r4-r5-r6-r7-r8;	
	res+=resalt;
	//D resalt= r1+r2alt+r3+r4alt-r5-r6alt-r7-r8alt;	
	// very obvious at importnant parameter values, now try the alts 
	if( debug_resvar) { 
	D sspexp= r1+r3-r5-r7;	
	MM_MSG(MMPR4(res,sspexp,ssp*rerff,sspc*rerff)<<MMPR4(yff,yii,yfi,yif)<<MMPR4(sf.xpu,sf.xpl,sf.xu,sf.xl))
	MM_MSG(MMPR(resold)<<MMPR4(r1,r2,r3,r4)<<MMPR4(r5,r6,r7,r8))
	MM_MSG(MMPR(resalt)<<MMPR4(r1,r2alt,r3,r4alt)<<MMPR4(r5,r6alt,r7,r8alt))
	MM_MSG(MMPR((r2altp+r4altp-r6altp-r8altp))<<MMPR4(r2altp,r4altp,r6altp,r8altp))
	}

	if (sf.olap())
	{
		Po R,S;
		mi.multiply_polynomials(R,P,Q);
		mi.integrate_polynomial(S,R);
		// efficiency wtf 
		iz=mi.evaluate_polynomial(S,sf.overu)
		 -mi.evaluate_polynomial(S,sf.overl);
	}
	// the whole thing ignored a -dy in the change of variable 
	res=-res;
	const D ress=res;
	//res=A*res-0*iz; // +t2+t1+t0+t+tkluge1;
//	res=1.0*res-0*iz; // +t2+t1+t0+t+tkluge1;
	res=1.0*res-1.0*iz; // +t2+t1+t0+t+tkluge1;
	if (debug_step_3) { ss<<MMPR(iz); MM_MSG(ss.str()) } 
	if (debug_ltz) if (ress<0)
	{
		MM_MSG(" mb_fl_diff_poly<0 "<<MMPR4(ress,res,iz,beta)<<MMPR4( _xpu, _xpl,  _xu, _xl)<<MMPR4( vt, qpneg, qneg,base))

	}
	return res;
} // fl_mb_diff_comp


//////////////////////////////////////////////////////////////////////////////
// this looks non-negative down to about 1e-300
// TODO FIXME use of dense 2d polynomial things is stupid, max 2-3 terms out of 9 or so lol 
static  D fl_mb_diff_comp2( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=false;
	const bool debug_step_2=false;
	const bool debug_step_2b=false;
	const bool debug_step_3=false;
	const bool debug_step_3r=false;
	const bool debug_resvar=false;
//	const bool debug_pieces=!false;
	const bool debug_ltz=!false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);
	const D bsyfi=bs*yfi; // (sf.xpl-sf.xu-vt);
	const D bsyii=bs*yii; // (sf.xpl-sf.xl-vt);

if (debug_step_1){ 
	MM_MSG(MMPR4(sf.xu,sf.xl,sf.xpu,sf.xpl)<<MMPR(vt))
	MM_MSG(MMPR4(bsyff,bsyfi,bsyif,bsyii))
}

	Ss ss; 
	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	//typedef shape_function_integrals::simple_polynomial<D> Po;
	typedef Pol::simple_polynomial<D> Po;
	typedef Pol::dense_2_polynomial<D> D2;
	Po P,Q,pxferf,pxierf,pxfexp,pxiexp,pulim,pllim;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	if (debug_step_1) { MM_MSG(MMPR2(P.to_string(),Q.to_string())) }
	// TODO FIXME check the integration variable has to be sqrt(b)y NOT y
	// this only works for beta==1 unless care with erf expansion later 
	D2 yint(P,-1.0,1.0,-vt);
	D2 qint(Q,1.0,1.0,vt);
	D2 yy(1.0,1.0,vt);
	D2 yyx(1.0,-1.0,-vt);
	if (debug_step_1) { MM_MSG(MMPR(yint.to_string(0,"y","x'"))) } 
	if (debug_step_1) { MM_MSG(MMPR(yy.to_string(0,"y","x'"))) } 
	// integrate multiplied by exp(-b*y*y)
	const IdxTy erfsz=P.size()+1;
	const IdxTy expsz=P.size();

	D2 yerf(erfsz),yexp(expsz),yexp2,yerf2,yexp3,yerf3;
	// beta,factor, dir 
//integrate_e( dexp,derf,  x,  b,factor,  dir)
	const bool use_rerf=!false;
	if (use_rerf) D2::integrate_e_rnew(yexp,yerf,yint,beta,1,0);
	else D2::integrate_e_new(yexp,yerf,yint,beta,1,0);
	if (debug_step_2) { MM_MSG(" the exp and erf parts of x integral" ) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp.to_string(0,"y","x'"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf.to_string(0,"y","x'"))) } 
	// now change x' to x
	D2::composite(yexp2,yexp,yy,1);
	D2::composite(yerf2,yerf,yy,1);
	if (debug_step_2) { MM_MSG(" x integral as function of x and y " ) } 
	if (debug_step_2) { MM_MSG(MMPR(yy.to_string(0,"x","y"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	D2 temp,temp2;
	// do not evaluate yet, change variables and evaluate at end
	temp=qint.transpose();
	if (debug_step_2) { MM_MSG(MMPR(qint.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(temp.to_string(0,"x","y"))) } 
qint=temp;
	// the variable "x' " in yexp and yerf needs to be chagend to y
	D2::multiply(yexp,yexp2,qint);
	D2::multiply(yerf,yerf2,qint);
	if (debug_step_2b) { MM_MSG(" final integrand as f of y and x " ) } 
	if (debug_step_2b) { MM_MSG(MMPR(yexp.to_string(0,"y","x"))) } 
	if (debug_step_2b) { MM_MSG(MMPR(yerf.to_string(0,"y","x"))) } 

	// now the matrix should be in x an x' again 
	// multiply by Q(x') and integrate replacing x' with y
	const IdxTy isz=yexp.get_size()+3;
	yexp2.set_size(isz);
	yerf2.set_size(isz);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yexp,beta,1,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
	if (debug_step_3) { MM_MSG("  exp integrand  parts " ) } 
	if (debug_step_3) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// add the erf terms to those generated from exp 
	// note that this may be a numerical problem will need to check
	// priot to evaluation 
	//D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yerf,beta,0,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
	if (debug_step_3r) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3r) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// covert the erf polynomial form y and x to x' and x 
	D2::composite(yerf3,yerf2,yyx,0);
	if (debug_step_3r) { MM_MSG(MMPR(yerf3.to_string(0,"(x')","x"))) } 
//	for (D aoff=-5; aoff<6; aoff+=1)
//	{
		// needs rerf term 
//	D ssp=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu, erf(bsyii),erf(bsyff),
//		erf(bsyif),erf(bsyfi));
		D eff,eif,efi,eii, cff,cif,cfi,cii;
		aerfc(cff,eff,bsyff);
		aerfc(cii,eii,bsyii);
		aerfc(cfi,efi,bsyfi);
		aerfc(cif,eif,bsyif);
		IdxTy pflags;
		D sspc=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,
	 		//aerfc(bsyii),aerfc(bsyff), aerfc(bsyif),aerfc(bsyfi));
	 		eii,eff, eif,efi,&pflags);
		//D sspone=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,1.0,1.0,1.0,1.0);
	// I'm only sure this is right when there are no x*x' terms TODO FIXME 
		bool skip_const=(pflags!=0)&&(cii==cff)&&(cii==cif)&&(cii==cfi);
		D sspone=skip_const?0:yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,cii,cff,cif,cfi);
	//	MM_MSG(MMPR2(ssp,sspone))

	if (debug_step_3r){		MM_MSG(MMPR2(sspc,sspone)) } 
//	}
	yexp=yexp2;//.transpose();
	yerf=yerf2;//.transpose();
	const D rerff=(use_rerf)?((ps/bs*.5)):1.0;
	//const D yif=(sf.xpu-sf.xl-vt);
	const D r1=yexp.evaluate(yff,sf.xu)*exp(-beta*yff*yff);
	// the erf term is likely the problem. 
//	const D r2=yerf.evaluate(yff,sf.xu)*rerff*erf(bsyff);
//	const D r2alt=yerf3.evaluate(sf.xpu,sf.xu)*rerff*erf(bsyff);
//	const D r2altp=yerf3.evaluate(sf.xpu,sf.xu);
	const D r3=yexp.evaluate(yii,sf.xl)*exp(-beta*yii*yii);
//	const D r4=yerf.evaluate(yii,sf.xl)*rerff*erf(bsyii);
//	const D r4alt=yerf3.evaluate(sf.xpl,sf.xl)*rerff*erf(bsyii);
//	const D r4altp=yerf3.evaluate(sf.xpl,sf.xl);

	const D r5=yexp.evaluate(yfi,sf.xu)*exp(-beta*yfi*yfi);
//	const D r6=yerf.evaluate(yfi,sf.xu)*rerff*erf(bsyfi);
//	const D r6alt=yerf3.evaluate(sf.xpl,sf.xu)*rerff*erf(bsyfi);
//	const D r6altp=yerf3.evaluate(sf.xpl,sf.xu);
	const D r7=yexp.evaluate(yif,sf.xl)*exp(-beta*yif*yif);
//	const D r8=yerf.evaluate(yif,sf.xl)*rerff*erf(bsyif);
//	const D r8alt=yerf3.evaluate(sf.xpu,sf.xl)*rerff*erf(bsyif);
//	const D r8altp=yerf3.evaluate(sf.xpu,sf.xl);
	
	//D resalt= r1+r3-r5-r7+rerff*ssp;	
	//D resalt= r1+r3-r5-r7-rerff*(sspc-sspone);	
	D resalt= r1+r3-r5-r7+rerff*(sspc+sspone);	
//	D resold= r1+r2+r3+r4-r5-r6-r7-r8;	
	res+=resalt;
	//D resalt= r1+r2alt+r3+r4alt-r5-r6alt-r7-r8alt;	
	// very obvious at importnant parameter values, now try the alts 
//	if( debug_resvar) { 
//	MM_MSG(MMPR(res)<<MMPR4(yff,yii,yfi,yif)<<MMPR4(sf.xpu,sf.xpl,sf.xu,sf.xl))
//	MM_MSG(MMPR(resold)<<MMPR4(r1,r2,r3,r4)<<MMPR4(r5,r6,r7,r8))
//	MM_MSG(MMPR(resalt)<<MMPR4(r1,r2alt,r3,r4alt)<<MMPR4(r5,r6alt,r7,r8alt))
//	MM_MSG(MMPR((r2altp+r4altp-r6altp-r8altp))<<MMPR4(r2altp,r4altp,r6altp,r8altp))
//	}
// this needs to be fixed, it is NOT positive defininte as
// IIRC more can flow out that in . 
	if (false) if (sf.olap())
	{
		Po R,S;
		mi.multiply_polynomials(R,P,Q);
		mi.integrate_polynomial(S,R);
		// efficiency wtf 
		iz=mi.evaluate_polynomial(S,sf.overu)
		 -mi.evaluate_polynomial(S,sf.overl);
		iz=iz*A;
	}
	// the whole thing ignored a -dy in the change of variable 
	res=-res;
	const D ress=res;
	//res=A*res-0*iz; // +t2+t1+t0+t+tkluge1;
//	res=1.0*res-0*iz; // +t2+t1+t0+t+tkluge1;
	res=1.0*res-1.0*iz; // +t2+t1+t0+t+tkluge1;
	// TODO FIXME on quick testing the only negatives were e-320 or so
	// that is tolerable for now
	if (res<0) if (res>-1.0e-200 ) res=0; 
	if (debug_step_3) { ss<<MMPR(iz); MM_MSG(ss.str()) } 
	if (debug_ltz) if ((ress<0)||(res<0))
	{
		MM_MSG(" mb_fl_diff_comp2<0 "<<MMPR4(ress,res,iz,beta)<<MMPR4( _xpu, _xpl,  _xu, _xl)<<MMPR4( vt, qpneg, qneg,base)<<MMPR4(sspc,sspone,(_xpu-_xu),(_xpl-_xl)))

	}
	return res;
} // fl_mb_diff_comp2


//////////////////////////////////////////////////////////////////////////////

static  D fl_mb_diff_comp_sparse( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=false;
	const bool debug_step_2=false;
	const bool debug_step_2b=false;
	const bool debug_step_3=false;
	const bool debug_step_3r=false;
	const bool debug_resvar=false;
//	const bool debug_pieces=!false;
	const bool debug_ltz=!false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);
	const D bsyfi=bs*yfi; // (sf.xpl-sf.xu-vt);
	const D bsyii=bs*yii; // (sf.xpl-sf.xl-vt);

if (debug_step_1){ 
	MM_MSG(MMPR4(sf.xu,sf.xl,sf.xpu,sf.xpl)<<MMPR(vt))
	MM_MSG(MMPR4(bsyff,bsyfi,bsyif,bsyii))
}

	Ss ss; 
	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	//typedef shape_function_integrals::simple_polynomial<D> Po;
	typedef Pol::simple_polynomial<D> Po;
	//typedef Pol::dense_2_polynomial<D> D2;
	typedef sparse_2nomial<D> D2;
	Po P,Q,pxferf,pxierf,pxfexp,pxiexp,pulim,pllim;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	if (debug_step_1) { MM_MSG(MMPR2(P.to_string(),Q.to_string())) }
	// TODO FIXME check the integration variable has to be sqrt(b)y NOT y
	// this only works for beta==1 unless care with erf expansion later 
	D2 yint(P,-1.0,1.0,-vt);
	D2 qint(Q,1.0,1.0,vt);
	D2 yy(1.0,1.0,vt);
	D2 yyx(1.0,-1.0,-vt);
	if (debug_step_1) { MM_MSG(MMPR(yint.to_string(0,"y","x'"))) } 
	if (debug_step_1) { MM_MSG(MMPR(yy.to_string(0,"y","x'"))) } 
	// integrate multiplied by exp(-b*y*y)
	const IdxTy erfsz=P.size()+1;
	const IdxTy expsz=P.size();

	D2 yerf(erfsz),yexp(expsz),yexp2,yerf2,yexp3,yerf3;
	// beta,factor, dir 
//integrate_e( dexp,derf,  x,  b,factor,  dir)
	const bool use_rerf=!false;
	if (use_rerf) D2::integrate_e_rnew(yexp,yerf,yint,beta,1,0);
	else D2::integrate_e_new(yexp,yerf,yint,beta,1,0);
	if (debug_step_2) { MM_MSG(" the exp and erf parts of x integral" ) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp.to_string(0,"y","x'"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf.to_string(0,"y","x'"))) } 
	// now change x' to x
	D2::composite(yexp2,yexp,yy,1);
	D2::composite(yerf2,yerf,yy,1);
	if (debug_step_2) { MM_MSG(" x integral as function of x and y " ) } 
	if (debug_step_2) { MM_MSG(MMPR(yy.to_string(0,"x","y"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	D2 temp,temp2;
	// do not evaluate yet, change variables and evaluate at end
	temp=qint.transpose();
	if (debug_step_2) { MM_MSG(MMPR(qint.to_string(0,"y","x"))) } 
	if (debug_step_2) { MM_MSG(MMPR(temp.to_string(0,"x","y"))) } 
qint=temp;
	// the variable "x' " in yexp and yerf needs to be chagend to y
	D2::multiply(yexp,yexp2,qint);
	D2::multiply(yerf,yerf2,qint);
	if (debug_step_2b) { MM_MSG(" final integrand as f of y and x " ) } 
	if (debug_step_2b) { MM_MSG(MMPR(yexp.to_string(0,"y","x"))) } 
	if (debug_step_2b) { MM_MSG(MMPR(yerf.to_string(0,"y","x"))) } 

	// now the matrix should be in x an x' again 
	// multiply by Q(x') and integrate replacing x' with y
	const IdxTy isz=yexp.get_size()+3;
	yexp2.set_size(isz);
	yerf2.set_size(isz);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0);
	//D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yexp,beta,1,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yexp,beta,1,0,false);
	if (debug_step_3) { MM_MSG("  exp integrand  parts " ) } 
	if (debug_step_3) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// add the erf terms to those generated from exp 
	// note that this may be a numerical problem will need to check
	// priot to evaluation 
	//D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
if (use_rerf)	D2::integrate_e_rnew(yexp2,yerf2,yerf,beta,0,0,false);
else 	D2::integrate_e_new(yexp2,yerf2,yerf,beta,0,0,false);
	if (debug_step_3r) { MM_MSG(MMPR(yexp2.to_string(0,"y","x"))) } 
	if (debug_step_3r) { MM_MSG(MMPR(yerf2.to_string(0,"y","x"))) } 
	// covert the erf polynomial form y and x to x' and x 
	D2::composite(yerf3,yerf2,yyx,0);
	if (debug_step_3r) { MM_MSG(MMPR(yerf3.to_string(0,"(x')","x"))) } 
//	for (D aoff=-5; aoff<6; aoff+=1)
//	{
		// needs rerf term 
//	D ssp=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu, erf(bsyii),erf(bsyff),
//		erf(bsyif),erf(bsyfi));
		D eff,eif,efi,eii, cff,cif,cfi,cii;
		aerfc(cff,eff,bsyff);
		aerfc(cii,eii,bsyii);
		aerfc(cfi,efi,bsyfi);
		aerfc(cif,eif,bsyif);
		IdxTy pflags;
		D sspc=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,
	 		//aerfc(bsyii),aerfc(bsyff), aerfc(bsyif),aerfc(bsyfi));
	 		eii,eff, eif,efi,&pflags);
		//D sspone=yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,1.0,1.0,1.0,1.0);
	// I'm only sure this is right when there are no x*x' terms TODO FIXME 
		bool skip_const=(pflags!=0)&&(cii==cff)&&(cii==cif)&&(cii==cfi);
		D sspone=skip_const?0:yerf3.evaluate_four(sf.xpl,sf.xpu,sf.xl,sf.xu,cii,cff,cif,cfi);
	//	MM_MSG(MMPR2(ssp,sspone))

	if (debug_step_3r){		MM_MSG(MMPR2(sspc,sspone)) } 
//	}
	yexp=yexp2;//.transpose();
	yerf=yerf2;//.transpose();
	const D rerff=(use_rerf)?((ps/bs*.5)):1.0;
	//const D yif=(sf.xpu-sf.xl-vt);
	const D r1=yexp.evaluate(yff,sf.xu)*exp(-beta*yff*yff);
	// the erf term is likely the problem. 
//	const D r2=yerf.evaluate(yff,sf.xu)*rerff*erf(bsyff);
//	const D r2alt=yerf3.evaluate(sf.xpu,sf.xu)*rerff*erf(bsyff);
//	const D r2altp=yerf3.evaluate(sf.xpu,sf.xu);
	const D r3=yexp.evaluate(yii,sf.xl)*exp(-beta*yii*yii);
//	const D r4=yerf.evaluate(yii,sf.xl)*rerff*erf(bsyii);
//	const D r4alt=yerf3.evaluate(sf.xpl,sf.xl)*rerff*erf(bsyii);
//	const D r4altp=yerf3.evaluate(sf.xpl,sf.xl);

	const D r5=yexp.evaluate(yfi,sf.xu)*exp(-beta*yfi*yfi);
//	const D r6=yerf.evaluate(yfi,sf.xu)*rerff*erf(bsyfi);
//	const D r6alt=yerf3.evaluate(sf.xpl,sf.xu)*rerff*erf(bsyfi);
//	const D r6altp=yerf3.evaluate(sf.xpl,sf.xu);
	const D r7=yexp.evaluate(yif,sf.xl)*exp(-beta*yif*yif);
//	const D r8=yerf.evaluate(yif,sf.xl)*rerff*erf(bsyif);
//	const D r8alt=yerf3.evaluate(sf.xpu,sf.xl)*rerff*erf(bsyif);
//	const D r8altp=yerf3.evaluate(sf.xpu,sf.xl);
	
	//D resalt= r1+r3-r5-r7+rerff*ssp;	
	//D resalt= r1+r3-r5-r7-rerff*(sspc-sspone);	
	D resalt= r1+r3-r5-r7+rerff*(sspc+sspone);	
//	D resold= r1+r2+r3+r4-r5-r6-r7-r8;	
	res+=resalt;
	//D resalt= r1+r2alt+r3+r4alt-r5-r6alt-r7-r8alt;	
	// very obvious at importnant parameter values, now try the alts 
//	if( debug_resvar) { 
//	MM_MSG(MMPR(res)<<MMPR4(yff,yii,yfi,yif)<<MMPR4(sf.xpu,sf.xpl,sf.xu,sf.xl))
//	MM_MSG(MMPR(resold)<<MMPR4(r1,r2,r3,r4)<<MMPR4(r5,r6,r7,r8))
//	MM_MSG(MMPR(resalt)<<MMPR4(r1,r2alt,r3,r4alt)<<MMPR4(r5,r6alt,r7,r8alt))
//	MM_MSG(MMPR((r2altp+r4altp-r6altp-r8altp))<<MMPR4(r2altp,r4altp,r6altp,r8altp))
//	}
// this needs to be fixed, it is NOT positive defininte as
// IIRC more can flow out that in . 
	if (false) if (sf.olap())
	{
		Po R,S;
		mi.multiply_polynomials(R,P,Q);
		mi.integrate_polynomial(S,R);
		// efficiency wtf 
		iz=mi.evaluate_polynomial(S,sf.overu)
		 -mi.evaluate_polynomial(S,sf.overl);
		iz=iz*A;
	}
	// the whole thing ignored a -dy in the change of variable 
	res=-res;
	const D ress=res;
	//res=A*res-0*iz; // +t2+t1+t0+t+tkluge1;
//	res=1.0*res-0*iz; // +t2+t1+t0+t+tkluge1;
	res=1.0*res-1.0*iz; // +t2+t1+t0+t+tkluge1;
	// TODO FIXME on quick testing the only negatives were e-320 or so
	// that is tolerable for now
	if (res<0) if (res>-1.0e-200 ) res=0; 
	if (debug_step_3) { ss<<MMPR(iz); MM_MSG(ss.str()) } 
	if (debug_ltz) if ((ress<0)||(res<0))
	{
		MM_MSG(" mb_fl_diff_comp2<0 "<<MMPR4(ress,res,iz,beta)<<MMPR4( _xpu, _xpl,  _xu, _xl)<<MMPR4( vt, qpneg, qneg,base)<<MMPR4(sspc,sspone,(_xpu-_xu),(_xpl-_xl)))

	}
	return res;
} // fl_mb_diff_comp_sparse





////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
// lanczos integral by derivative doh 
static  D fl_mb_diff_lanczos( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=false;
	const bool debug_step_2=false;
	const bool debug_step_3=false;
	const bool debug_pieces=false;
	const bool debug_ltz=false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);

	const D dxp=sf.xpu-sf.xpl;
	const D dx=sf.xu-sf.xl;
	const D xphat=.5*(sf.xpu+sf.xpl);
	const D xhat=.5*(sf.xu+sf.xl);
	const D yhat=xphat-xhat-vt;
	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	typedef Pol::simple_polynomial<D> Po;
	typedef Pol::dense_2_polynomial<D> D2;
	//typedef std::vector<Po> Pstack;
	typedef std::vector<D> Dstack;
	typedef std::vector<D2> D2stack;
	typedef std::vector< D2stack > Ystack;
	Po P,Q;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	D2 pqy(P,Q);
	D2 yint(-1.0,1.0,-vt);
	Ystack ystack;
	ystack.push_back( D2stack());
	ystack[0].push_back(pqy);
//	MM_MSG("FUDDKKKK "<<MMPR3(ystack[0].size(), pqy.to_string(),ystack[0][0].to_string()))
	if (debug_step_1) 
{ MM_MSG(MMPR2(P.to_string(),Q.to_string())) }
//	Pstack pstack,qstack;	
//	pstack.push_back(P); qstack.push_back(Q);
	Dstack dxs,dxps;
	dxs.push_back(dx);
	dxps.push_back(dxp);
	D term=0;
	D fac=1;
	IdxTy level=0;
	IdxTy nlevel=1;
	const D c=exp(-beta*yhat*yhat);
	while (true)
	{
		for(IdxTy i=1; i<=level; i+=2)
		{
		const IdxTy j=level-i;
		while (i>=dxs.size()) { dxs.push_back(dx*dxs[dxs.size()-1]); }
		while (j>=dxps.size()) { dxps.push_back(dxp*dxps[dxps.size()-1]); }
		while (i>=ystack.size())
		{	
			const IdxTy ix=ystack.size()-1;
			const D2 & d2=ystack[ix][0];
//			MM_MSG(MMPR3(d2.to_string(),ystack[ix].size(),ix))
			D2 pq,temp;
			D2::differentiate(temp,d2,1);
//			MM_MSG(MMPR2(temp.to_string(),d2.to_string()))
			D2::multiply(pq,d2,yint,-2.0*beta);
			pq.accumulate(temp);
			// TODO convert d2 into pq doh 	
			ystack.push_back(D2stack());
			ystack[ix+1].push_back(pq);
//			MM_MSG(MMPR(ystack[ix+1][0].to_string()))
//			MM_MSG(MMPR4(ystack.size(),pq.to_string(),temp.to_string(),yint.to_string()))
		}
		D2stack & d2s=ystack[i-1];
		while (j>=d2s.size())
		{
			const IdxTy ix=d2s.size()-1;
			const D2 & d2=d2s[ix]; // [iy];
			D2 pq,temp;
			D2::differentiate(temp,d2,0);
//			MM_MSG(MMPR2(temp.to_string(),d2.to_string()))
			D2::multiply(pq,d2,yint,2.0*beta);
			pq.accumulate(temp);
			d2s.push_back(pq);
//			MM_MSG(MMPR4(d2s.size(),pq.to_string(),temp.to_string(),yint.to_string()))
		}

// TODO 		fac = 1/i!/j! 
			fac=4.0;
			for (IdxTy crap=1; crap<=i; ++crap) fac=.5*fac/D(crap);
			for (IdxTy crap=1; crap<=j; ++crap) fac=.5*fac/D(crap);
			//for (IdxTy crap=1; crap<=(j+i); ++crap) fac=fac/D(crap);

		const D & dxi=dxs[i-1];
		const D & dxpi=dxps[j-1];
		// these can be cached too 
		//const D & ymat=ystack[i][j].evaluate(xhat,yhat);
	if (debug_step_2) { 
		MM_MSG(MMPR(d2s[j-1].to_string()))}
		const D & ymat=d2s[j-1].evaluate(xhat,xphat);
		//term=fac*c*Pi*Qi*dxi*dxpi;
		term=fac*c*ymat*dxi*dxpi;
		res+=term;
	if (debug_step_2) { 
		MM_MSG(MMPR4(i,j,term,res)<<MMPR4(dxi,dxpi,ymat,xhat)<<MMPR2(xphat,fac))
}	

	} // i 
		level+=2;
		nlevel+=1;
		if (nlevel>7) break; 
	} // true 

	return res;

} //   fl_mb_diff_lanczos

//////////////////////////////////////////

// expand integrand around peak and approximate that way 
// may be similar to steepest descents  
static  D fl_mb_diff_peak( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
	const bool debug_step_1=!false;
	const bool debug_step_2=!false;
	const bool debug_step_3=!false;
	const bool debug_pieces=!false;
	const bool debug_ltz=!false;
	D res=0;
	const bool canonize=true; // remove offsets from coordinates 
	static const  D ps=::sqrt(M_PI);
	D iz=0;

	D base=0;	
	if (canonize)
	{
		base=-.5*(_xu+_xl);
		base=-.25*(_xu+_xl+_xpu+_xpl);
		MM_ONCE(" canonicalization code needs to be sign checked "<< base,)
	} // canonize
 	const fl_sf_auto sf( beta, _xpu, _xpl,  _xu, _xl, vt, qpneg, qneg,base);
	const D bs=sqrt(beta);
	const D A=bs/ps; // 1.0; // ps/bs;
	// xFIXME xTODO not the subscript convention here is NOT consistent 
	const D yff=(sf.xpu-sf.xu-vt);
	const D yif=(sf.xpu-sf.xl-vt);
	const D yfi=(sf.xpl-sf.xu-vt);
	const D yii=(sf.xpl-sf.xl-vt);
	const D bsyff=bs*yff; // (sf.xpu-sf.xu-vt);
	const D bsyif=bs*yif; // (sf.xpu-sf.xl-vt);

	const D dxp=sf.xpu-sf.xpl;
	const D dx=sf.xu-sf.xl;
	const D xphat=.5*(sf.xpu+sf.xpl);
	const D xhat=.5*(sf.xu+sf.xl);
	const D yhat=xphat-xhat-vt;
	static mjm_integrals mi;
	typedef shape_function_integrals Pol;
	typedef Pol::simple_polynomial<D> Po;
	typedef Pol::dense_2_polynomial<D> D2;
	typedef std::vector<D> Dstack;
	//typedef std::vector<D2> D2stack;
	//typedef std::vector< D2stack > Ystack;
	//typedef std::vector<Po> Pstack;
	Po P,Q;
	P.push_back(sf.p0);
	P.push_back(sf.p1);
	Q.push_back(sf.q0);
	Q.push_back(sf.q1);
	D2 pqy(P,Q);
	D2 yint(-1.0,1.0,-vt);
	Dstack dxs,dxps;
	dxs.push_back(dx);
	dxps.push_back(dxp);
	D term=0;
	D fac=1;
	IdxTy level=0;
	IdxTy nlevel=1;
	const D c=exp(-beta*yhat*yhat);
	while (true)
	{
		for(IdxTy i=1; i<=level; i+=2)
		{
		const IdxTy j=level-i;
		while (i>=dxs.size()) { dxs.push_back(dx*dxs[dxs.size()-1]); }
		while (j>=dxps.size()) { dxps.push_back(dxp*dxps[dxps.size()-1]); }
		const D & dxi=dxs[i-1];
		const D & dxpi=dxps[j-1];
		} // i 
		level+=2;
		nlevel+=1;
		if (nlevel>7) break; 
	} // true 

	return res;

} //   fl_mb_diff_peak



/////////////////////////////////////////

// rederivce for simple case with xp as destination outer integral
// and x inner integral of source. Polynomaials P(x) for source and
// Q(x') foe xp destination. iHard coded for now, fix when numerical
// problems occur lol
// \int dx' Q(x') \int dx P(x) exp(-beta*(x'-x-vt)^2)
// I0= \int exp(-beta*y*y) = cf2*erf(y*sqrt(beta))
// I1= \int y exp(-beta*y*y) = cf1*exp(-beta*y*y)
// coefs:
// A5 : y*y*I0, A4 : y*I1 , A3: y*I0, A2 : I1, A1: I0
// http://www.wolframalpha.com/input/?i=integrate+y*erf%28sqrt%28b%29*y%29 
static  D esimple( const D & beta, const D & _xpu,const D & _xpl, const D & _xu, const D & _xl 
, const D & vt, const bool qpneg,const bool qneg)
{
//const D & p1m=1.0;
// p0m=1.0, const D & q1m=1.0, const D & q0m=1.0 ) 
	const bool dump_coefs=!true;
	const bool dump_js=!true;
	const bool dump_limits=!true;
	const bool check_small_bx=!true; // check the small value expression for comparison to "exact"
	// turning these on does eliminate a lot of spurious negative values that I THOUGHT
	// were in the limit as beta*x goes to infinity NOT zero although the OFFSETS
	// are a huge problem and maybe the normalize PLUS expansion works. 
	const bool canonize=true; // remove offsets from coordinates 
	bool small_ok=true; // ok to use expansion for small parameters 
	bool big_ok=!true; // ok to use big approximatin 
	static const  D ps=::sqrt(M_PI);

	D base=0;	
	// optionally make variable origins canonical to reduce noise and numeric crap
	if (canonize)
	{	// this is not right but seemed to work ??? TODO FIXME 
		base=-.25*(_xpu+_xpl+_xu+_xl);
		MM_ONCE(" canonicalization code needs to be sign checked",)
	}
	const D xpu=_xpu+base;
	const D xpl=_xpl+base;
	const D xu=_xu+base;
	const D xl=_xl+base;

	// for quadrant zero to zero, the source x and destination d are at lower end of range
	const D dx=xu-xl;
	const D dxp=xpu-xpl;
// variable change is now buried 	
	// negative means pgoressing from upper limit hence sign is right 
	const D q1=(qpneg?1:(-1.0))/dxp;
	// q0 depends on xmin or xmax... 
	const D q0=(qpneg?(-xpl):xpu)/dxp;
	const D p1=(qneg?1:(-1.0))/dx;
//	const D p1=-1.0/dx;
	const D p0=(qneg?(-xl):xu)/dx;

	if (dump_limits ) { MM_MSG(MMPR4(xpu,xpl,xu,xl)<<MMPR4(beta,vt,qpneg,qneg)<<MMPR4(q1,q0,p1,p0))}
	if (beta==0)
	{
		// this is not right as the limit, that is there needs to be a normalization
		// as beta -> zero. 
		D res=(.5*p1*(xu+xl)+p0)*(xu-xl);
		res*=(.5*q1*(xpu+xpl)+q0)*(xpu-xpl);
		MM_MSG(" beta is zero, danger will robinson "<<MMPR4(xpu,xpl,xu,xl)<<MMPR(res))
		return res;
	}

	const D bs=::sqrt(beta);  
	const D bsyff=bs*(xpu-xu-vt);
	const D bsyfi=bs*(xpu-xl-vt);
	const D bsyif=bs*(xpl-xu-vt);
	const D bsyii=bs*(xpl-xl-vt);
	// this needs to be more accurate limiting the ranges better
	// note that this is exp(-36) which is not that small 
	// right at the resolution limit, 1+delta need to fix this crap 
	//> exp(-36)
	//[1] 2.319523e-16
	const D big_lim=5.5-.5;	 
	if (big_ok)
	{
		const bool not_cross_zed=((bsyfi>0)&&(bsyif>0))||((bsyfi<0)&&(bsyif<0));
		if (!not_cross_zed) big_ok=false;
	}
	if (big_ok) {big_ok&=(bsyff>big_lim)||(bsyff<-big_lim);
	if (big_ok) {big_ok&=(bsyfi>big_lim)||(bsyfi<-big_lim);
	if (big_ok) {big_ok&=(bsyif>big_lim)||(bsyif<-big_lim);
	if (big_ok) {big_ok&=(bsyii>big_lim)||(bsyii<-big_lim);
}}}}
	if (big_ok) return 0; 

	const D small_lim=2e-2;
	if (small_ok) {small_ok&=(bsyff<small_lim)&&(bsyff>-small_lim);
	if (small_ok){ small_ok&=(bsyfi<small_lim)&&(bsyfi>-small_lim);
	if (small_ok){ small_ok&=(bsyif<small_lim)&&(bsyif>-small_lim);
	if (small_ok){ small_ok&=(bsyii<small_lim)&&(bsyii>-small_lim);
}} }}
	const D xl2=xl*xl;
	const D xu2=xu*xu;
	D ressbx=0;
	// in the beta -> zero limit this looks more stable the other one starts to oscillate
	if (check_small_bx||small_ok)
	{ // taylor expand the exp and do polynomial integration
		//\int dxp \int dx  (q1*xp+q0)*(p1*x+p0)*(1-b*xp*xp-b*(x+vt)*(x+vt)+2*b*xp*(x+vt))
		const D xv= xu-xl; // :)
		const D xv2= xu2-xl2;
		const D xv3=xu2*xu-xl2*xl;
		const D xv4=xu2*xu2-xl2*xl2;
		const D w2=-beta;
		const D Xp0=.25*p1*w2*xv4+xv3/3.0*(p0*w2+p1*2.0*vt*w2)+.5*xv2*(p1-p1*beta*vt*vt-2.0*p0*vt*beta)+p0*xv-p0*xv*vt*vt*beta;	
		const D Xp1=p1/3.0*2*beta*xv3+.5*xv2*(p0*2.0*beta+2.0*p1*beta*vt)+ 2.0*xv*p0*beta*vt;
		const D Xp2=.5*xv2*p1*w2+p0*xv*w2;
		const D xpu2=xpu*xpu;
		const D xpl2=xpl*xpl;
		const D t3=(xpu2+xpl2)*(xpu2-xpl2)*.25;
		const D t2=(xpu*xpu*xpu-xpl*xpl*xpl)/3.0;
		const D t1=(xpu*xpu-xpl*xpl)*.5;
		const D t0=(xpu-xpl);
//   (t3*(Xp2*q1))=-914.149 (t2*(Xp2*q0+q1*Xp1))=2741.67 (t1*(q1*Xp0+q0*Xp1))=-384853 
// (t0*(q0*Xp0))=382068 t3=794.02 t2= 39.8003 t1=1.995 t0=0.1
if (check_small_bx)
{		MM_MSG(MMPR4((t3*(Xp2*q1)),(t2*(Xp2*q0+q1*Xp1)),(t1*(q1*Xp0+q0*Xp1)),(t0*(q0*Xp0)))<<MMPR4(t3,t2,t1,t0)) ; }	 
		ressbx=t3*(Xp2*q1)+t2*(Xp2*q0+q1*Xp1)+t1*(q1*Xp0+q0*Xp1)+t0*(q0*Xp0);	 
		if (small_ok) 
		{
			if (ressbx<0) { MM_MSG(" res<0 "<<MMPR4(xpu,xpl,xu,xl)<<MMPR4(base,beta,vt,ressbx)<<MMPR4(q1,q0,p1,p0)) }
			return ressbx;
		}
	}


	const D cf1=-.5/beta;
	const D cf2=.5*ps/bs;

	const D A5=-p1*q1;
	const D A4=p1*q1;
	//const D A30=2.0*vt*p1*q1+q1*p0-q1*p1*vt+q0*p1;
	//const D A30=vt*p1*q1+q1*p0+q0*p1;
	const D A30=-vt*p1*q1-q1*p0-q0*p1;
	//const D A31=2.0*p1*q1;
	const D A31=-2.0*p1*q1;
	//const D A20=-q1*p1*vt-q0*p1;
	const D A20=q1*p1*vt+q0*p1;
	const D A21=q1*p1;
	//const D A1c=q0*p1+q1*p0-q1*p1*vt;
	//const D A12=p1*q1;
	////const D A12=p1*q1;
	const D A12=-p1*q1;
	//const D A11=2.0*vt*p1*q1+A1c;
	////const D A11=q1*p0+p1*(q0+q1*vt);
	const D A11=-(q1*p0+p1*(q0+q1*vt));
	//const D A10=p1*q1*vt*vt+vt*A1c+q0*p0-q0*p1*vt;
	////const D A10=q0*p0+p0*q1*vt;
	const D A10=-(q0*p0+p0*q1*vt);

	if (dump_coefs ) { MM_MSG(MMPR4(A5,A4,A30,A31)<<MMPR4(A20,A21,A12,A11)<<MMPR(A10))}


	const D bsyff3=bsyff*bsyff*bsyff;
	const D bsyfi3=bsyfi*bsyfi*bsyfi;
	const D bsyif3=bsyif*bsyif*bsyif;
	const D bsyii3=bsyii*bsyii*bsyii;

	const D eyff=exp(-bsyff*bsyff);
	const D eyfi=exp(-bsyfi*bsyfi);
	const D eyif=exp(-bsyif*bsyif);
	const D eyii=exp(-bsyii*bsyii);
	const D eryff=erf(bsyff);
	const D eryfi=erf(bsyfi);
	const D eryif=erf(bsyif);
	const D eryii=erf(bsyii);

	const D b2v=bsyff*bsyff*eryff-bsyfi*bsyfi*eryfi-bsyif*bsyif*eryif+bsyii*bsyii*eryii;
	const D b2vx=bsyff*bsyff*eryff*xu-bsyfi*bsyfi*eryfi*xl-bsyif*bsyif*eryif*xu+bsyii*bsyii*eryii*xl;
	const D expb2v=bsyff*bsyff*eyff-bsyfi*bsyfi*eyfi-bsyif*bsyif*eyif+bsyii*bsyii*eyii;
	const D expv=(eyff-eyfi-eyif+eyii);
	const D expvx=(eyff*xu-eyfi*xl-eyif*xu+eyii*xl);
	const D expvxx=(eyff*xu2-eyfi*xl2-eyif*xu2+eyii*xl2);
//	const D expvq=(eyff*q0u-eyfi*q0l-eyif*q0u+eyii*q0l);
//else { 		mi.doen_split(res,dnorm, bsyf,0,0,false); res*=.5*ps*q0/bs; }

	const D erfv=(eryff-eryfi-eryif+eryii);
	const D erfvx=(eryff*xu-eryfi*xl-eryif*xu+eryii*xl);
	const D erfvxx=(eryff*xu2-eryfi*xl2-eryif*xu2+eryii*xl2);
//	const D erfvq=(eryff*q0u-eryfi*q0l-eryif*q0u+eryii*q0l);
	const D expyv=(bsyff*eyff-bsyfi*eyfi-bsyif*eyif+bsyii*eyii);
	const D expyvx=(bsyff*eyff*xu-bsyfi*eyfi*xl-bsyif*eyif*xu+bsyii*eyii*xl);
	const D expyv3=(bsyff3*eyff-bsyfi3*eyfi-bsyif3*eyif+bsyii3*eyii);
	const D erfyv3=(bsyff3*eryff-bsyfi3*eryfi-bsyif3*eryif+bsyii3*eryii);
	const D erfyv=(bsyff*eryff-bsyfi*eryfi-bsyif*eryif+bsyii*eryii);
	const D erfyvx=(bsyff*eryff*xu-bsyfi*eryfi*xl-bsyif*eryif*xu+bsyii*eryii*xl);
	const D erfyvxx=(bsyff*eryff*xu2-bsyfi*eryfi*xl2-bsyif*eryif*xu2+bsyii*eryii*xl2);
// Asymptote seems to be J5, J3x , and J1xx 
	//const D erfyvq=(bsyff*eryff*q0u-bsyfi*eryfi*q0l-bsyif*eryif*q0u+bsyii*eryii*q0l);
//J5: e^(-b y^2)/(3 sqrt(π) b^(3/2)) + 1/3 y^3 erf(sqrt(b) y) + (y^2 e^(-b y^2))/(3 sqrt(π) sqrt(b)) + constant
const D  J5v=cf2*( expv/(3*ps*beta*bs)+erfyv3/(3.0*beta*bs)+expb2v/(3*ps*bs*beta));
//J4, 1.0/(4*beta^2) * exp (-beta*y*Y)
//const D J4v=cf1*(.25/(beta*beta)*expv);
const D J4v=(.25/(beta*beta)*expv);
// J3 1/2 y^2 erf(sqrt(b) y) - (erf(sqrt(b) y))/(4 b) + (y e^(-b y^2))/(2 sqrt(π) sqrt(b)) + constant
const D  J3v=cf2*( .5*b2v/beta-.25/beta*erfv+ expyv/(2.0*ps*beta));
const D  J3vx=cf2*( .5*b2vx/beta-.25/beta*erfvx+ expyvx/(2.0*ps*beta));
//J2, -1/(2*beta)*sqrt(pi)/2/sqrt(beta)*erf(sqrt(b)*y)
const D J2v=cf1*(cf2*erfv);
const D J2vx=cf1*(cf2*erfvx);
// J1, //integral erf(sqrt(b) y) dy = y erf(sqrt(b) y) + e^(-b y^2)/(sqrt(π) sqrt(b)) + constant
const D J1v= cf2*(erfyv/bs+expv/(ps*bs)); 
const D J1vx= cf2*(erfyvx/bs+expvx/(ps*bs)); 
const D J1vxx= cf2*(erfyvxx/bs+expvxx/(ps*bs)); 

	D res=A5*J5v;
	res+=A4*J4v;
	res+=A30*J3v;
	res+=A31*J3vx;
	res+=A20*J2v;
	res+=A21*J2vx;
	res+=A10*J1v;
	res+=A11*J1vx;
	res+=A12*J1vxx;
	const D vxcc=A31*J3vx+A21*J2vx+A11*J1vx+A12*J1vxx;
	if (dump_js ) { MM_MSG(MMPR4(J5v,J4v,J3v,J3vx)<<MMPR4(J2v,J2vx,J1v,J1vx)<<MMPR3(J1vxx,res,vxcc))}
	// this is only the ix or iy component not the res returned to the caller after ix*iy 	
	if (check_small_bx)
	{MM_MSG("smallbx "<<MMPR4(res, ressbx,(ressbx-res),((ressbx-res)/(ressbx+res)))<<MMPR4(bsyff,bsyfi,bsyif,bsyii)<<MMPR(beta))	  }
	if (res<0)
	{
	static int bad_cnt=0;
			if ((bad_cnt%(1<<12))==0)  { MM_MSG(" res<0 "<<MMPR(bad_cnt)<<MMPR4(xpu,xpl,xu,xl)<<MMPR4(base,beta,vt,res)<<MMPR4(q1,q0,p1,p0)<<MMPR4(bsyff,bsyfi,bsyif,bsyii)) }

++bad_cnt;
	}
	return res; 

}


// evaluate double integral in 1D,  (c1*xp+c0)*(s1*x+s0)exp(-beta(xp-x-vxt)^2)
D evaluate1( const D & beta, const D & xpu,const D & xpl, const D & xu, const D & xl 
, const D & vxt , const D & s1, const D & s0, const D & c1, const D & c0, const D & xpbase=0,
const D & xbase=0) 
{
static const D N=1;
static const  D ps=::sqrt(M_PI);

//const D xbase=0; // try to fix numerical errors 
//const D xpbase=0; // try to fix numerical errors 
// TODO check the origins, this was a big deal in old code and shifts the
// shape function coefficients
// FIXME fudding sign mistake  and power of coef mistake  
const D vt0=-vxt; // v0*tau;
// this should have been fixed, done by caller 
// on a quick test, this is wrong. We need to square root beta
//const D aeff=beta; // not sure if this is beta or beta to a power lol 
const D aeff=::sqrt(beta); // not sure if this is beta or beta to a power lol 
// do any xbase translate here... 
const D a=0;// linear doping gradient B;
//const D b=C-a*xbase;
const D b=1.0; // doping is a CONSTANT deal with it later C; // -a*xbase;
const D c=c1;
const D d=c0;

// keep these numbers smalll 
const D xf=xu-xbase;
const D xi=xl-xbase;
//MM_MSG(" carrieres nf="<<(a*xf+b)<<" ni="<<(a*xi+b))
// for now these work best with same base 
const D xpf=xpu-xpbase;
const D xpi=xpl-xpbase;
// copied from old code, check these 

const D yff=xf-xpf-vt0;
const D yii=xi-xpi-vt0;
const D yif=xi-xpf-vt0;
const D yfi=xf-xpi-vt0;
// this is not right 
// fudd, c is the SOURCE and s is the DEST wtf 
const D syf=-s0-s1*(xf-vt0);
const D syi=-s0-s1*(xi-vt0);

// fix later lol 
const D ac=a*c;
const D ca=a*c;

const D cy1f=(-2.0*xf*ac-d*a-b*c);
const D cy1i=(-2.0*xi*ac-d*a-b*c);
const D cy0f=(d*b+ac*xf*xf+d*a*xf+b*c*xf);
const D cy0i=(d*b+ac*xi*xi+d*a*xi+b*c*xi);

const D jf0=ps*.5/aeff;
const D c0i3f=s1*ac;
const D c0i3i=s1*ac;
const D c0i2f=ac*syf+s1*cy1f;
const D c0i2i=ac*syi+s1*cy1i;
const D c0i1f=s1*(cy0f)+syf*(cy1f);
const D c0i1i=s1*(cy0i)+syi*(cy1i);
const D c0i0f=syf*cy0f;
const D c0i0i=syi*cy0i;

//const bool do_old_int=false;
//const bool do_alt_int=!false; 
//D J0alt=0;

const D jf1=-.5/(aeff*aeff);
const D c1i2f=s1*(-2.0*ca);
const D c1i2i=s1*(-2.0*ca);
const D c1i1f=s1*(-cy1f)+syf*(-2.0*ca);
const D c1i1i=s1*(-cy1i)+syi*(-2.0*ca);
const D c1i0f=syf*(-cy1f);
const D c1i0i=syi*(-cy1i);



const D jf2=1; // .5*a*c/(aeff*aeff);
const D a3=(1.0/(aeff*aeff*aeff));
const D a2=(1.0/(aeff*aeff));
const D c2i3f=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i3i=-s1*ac*.5*a2; // /aeff/aeff;
const D c2i2f=-syf*ac*.5*a2; // /aeff/aeff;
const D c2i2i=-syi*ac*.5*a2; // /aeff/aeff;
const D c2i1f=.25*ps*a3*s1*ac;
const D c2i1i=.25*ps*a3*s1*ac;
const D c2i0f=.25*ps*a3*syf*ac;
const D c2i0i=.25*ps*a3*syi*ac;

D  alt1=K2K1K0Ke2Ke1Ke0(yff,  yfi, aeff ,  c0i3f*jf0,c0i2f*jf0,  c0i1f*jf0+c2i1f*jf2,  c0i0f*jf0+c2i0f*jf2 ,
      c1i2f*jf1+c2i3f*jf2,  c1i1f*jf1+c2i2f*jf2,  c1i0f*jf1);
D  alt2= K2K1K0Ke2Ke1Ke0(yii,  yif, aeff ,c0i3i*jf0,  c0i2i*jf0,  c0i1i*jf0+c2i1i*jf2,  c0i0i*jf0+c2i0i*jf2 ,
      c1i2i*jf1+c2i3i*jf2,  c1i1i*jf1+c2i2i*jf2,  c1i0i*jf1);
D res=alt1+alt2 ; // +j0k3ai*jf0;
return res;
}





/////////////////////////////////////////////////////////////////////
// evaluate double integral in 1D,  (cxp+d)*(ax+b)exp(-beta(xp-x-vxt)^2)
// evaluate double integral in 1D,  (c1*p+c0)*(s1*x+s0)exp(-beta(xp-x-vxt)^2)
D evaluate( const D & beta, const D & xpu,const D & xpl, const D & xu, const D & xl 
, const D & vxt , const D & a, const D & b, const D & c, const D & d) 
{
//  change variables on gthe x integral to y=xp-x-vt, dy=-dx
// (a*(xp-y-vt)+b)exp(-beta*y^2))(-dy), ymin=xp-xmin-vt etc.

return 0; 
}

D evaluate_line( const D & xpc, const D & ypc, const D & xc, const int dir 
, const D & vxt, const D & vyt, const D & a, const D & b, const D & beta,const IdxTy pattern) 
{
D res=0;
//const D quads=1; 
const D quads=4; 
const IdxTy qbase=1;
switch (pattern)
{
case 256:
case 257:
case 0:
{
	for (IdxTy i=0; i<quads; ++i)
		res+= evaluate_ls( xpc, ypc, xc, dir, vxt, vyt, a, b,beta, i+qbase) ;
	break;
}
case 1:{res+=evaluate_ls(xpc,ypc,xc,dir,vxt,vyt,a,b,beta,0+qbase);break;}
case 2:{res+=evaluate_ls(xpc,ypc,xc,dir,vxt,vyt,a,b,beta,1+qbase);break;}
case 3:{res+=evaluate_ls(xpc,ypc,xc,dir,vxt,vyt,a,b,beta,2+qbase);break;}
case 4:{res+=evaluate_ls(xpc,ypc,xc,dir,vxt,vyt,a,b,beta,3+qbase);break;}

} // switch

return res;
}

//	const bool dumpep=false;
// evaluate the transfer from (x-a,a+a) to (xp+a,xp-a) with velocity distro
// exp (-beta*x^2) and vt displacement from x. Same origins for all x and y 
// evaluate_node_pair(  xpc,  ypc, xc, yc ,  vxt,  vyt,  a, b, beta)
D evaluate_node_pair( const D & xpc, const D & ypc, const D & xc, const D & yc
, const D & vxt, const D & vyt, const D & a, const D & b, const D & beta,const IdxTy pattern, const bool noverlap=true) 
{
//	const bool dumpep=false;
D res=0;
//const D quads=1; 
const D quads=4; 
const IdxTy qbase=1;
switch (pattern)
{
case 0:
{
	// this was zero based now changed 
	for (IdxTy i=0; i<quads; ++i)
	for (IdxTy ip=0; ip<quads; ++ip)
	res+= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, ip+qbase, i+qbase,noverlap) ;
	break;
}
case 1:{res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,0+qbase,0+qbase,noverlap);break;}
case 2:{res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,1+qbase,1+qbase,noverlap);break;}
case 3:{res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,2+qbase,2+qbase,noverlap);break;}
case 4:{res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,3+qbase,3+qbase,noverlap);break;}
// coded single ones
case (16+0): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,0+qbase,0+qbase,noverlap);break;}
case (16+1): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,0+qbase,1+qbase,noverlap);break;}
case (16+2): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,0+qbase,2+qbase,noverlap);break;}
case (16+3): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,0+qbase,3+qbase,noverlap);break;}
case (16+4): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,1+qbase,0+qbase,noverlap);break;}
case (16+5): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,1+qbase,1+qbase,noverlap);break;}
case (16+6): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,1+qbase,2+qbase,noverlap);break;}
case (16+7): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,1+qbase,3+qbase,noverlap);break;}
case (16+8): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,2+qbase,0+qbase,noverlap);break;}
case (16+9): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,2+qbase,1+qbase,noverlap);break;}
case (16+10): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,2+qbase,2+qbase,noverlap);break;}
case (16+11): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,2+qbase,3+qbase,noverlap);break;}
case (16+12): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,3+qbase,0+qbase,noverlap);break;}
case (16+13): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,3+qbase,1+qbase,noverlap);break;}
case (16+14): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,3+qbase,2+qbase,noverlap);break;}
case (16+15): {res+=evaluateq(xpc,ypc,xc,yc,vxt,vyt,a,b,beta,3+qbase,3+qbase,noverlap);break;}



// if noverlap, ignore quadrants that may be identical 
case 256:
{
IdxTy f[16];
// for now tolerances should not be a big deal a centers should be assigments
// from sames sources and not calculations 
if (noverlap)
{
const bool same=(xc==xpc)&&(yc==ypc);
//if (same){ for (IdxTy i=0; i<16; ++i)  f[i]=0; }
if (same){ for (IdxTy i=0; i<16; ++i)  if ((i&3)==(i>>2)) f[i]=0; else f[i]=1; }
else  // maybe just remove some quadrants
{
const D dx=xpc-xc;
const D dy=ypc-yc;
// these should be quantized, so now this is not real trick 
const D xtol=1.5*a;
const D ytol=1.5*b;
const bool clearx=(dx>xtol)||(dx<-xtol);
const bool cleary=(dy>ytol)||(dy<-ytol);
for (IdxTy i=0; i<16; ++i)  f[i]=1; 
if ((!clearx)&&!cleary)
{
const bool xgt=(xpc>(xc+.5*a));
const bool xlt=(xpc<(xc-.5*a));
const bool xeq=!(xgt||xlt); // risky but can fix later 

const bool ygt=(ypc>(yc+.5*b));
const bool ylt=(ypc<(yc-.5*b));
const bool yeq=!(ygt||ylt); // risky but can fix later 
 // tick tac toe now
if (xgt)
{
if (ygt){ // both gt, xp is NE , top 2 bits qp quadrant 
// skip qp-3 from q-1
//for (IdxTy i=0; i<4; ++i) {  f[i]=0;  f[i]=0; } 
f[2*4]=0;
} else if (ylt)  // q-4 and qp-2
{ f[3+4]=0;
} else // q-1 olap with qp-2 and q-4 olap qp-3
{
f[4]=0; f[3+(4*2)]=0; 

} // yeq


} // xgt
else if (xlt)
{
if (ygt){ // qp-4 and q-2
f[3*4+1]=0;
} else if (ylt) // qp-1 and q-3
{
f[2]=0;
} else // qp1 olap q2 and qp4 olap q3
{
f[1]=0; f[3*4+2]=0;
} // yeq


} // xlt
else // xeq
{
if (ygt){ // qp-3 olap q2 and qp4 olap q1
f[8+1]=0; f[4*3]=0;
} else if (ylt) // qp2 olap q3 and qp1 olap q4
{
f[4+2]=0; f[3]=0;
} else
{
MM_ONCE(" should not get here "<<MMPR4(xc,xpc,yc,ypc)<<MMPR2(a,b),)
MM_MSG(" should not get here "<<MMPR4(xc,xpc,yc,ypc)<<MMPR2(a,b))

} // yeq

} // xeq 

} // not clear in x or y 

} // same

} // noverlap
else {
for (IdxTy i=0; i<16; ++i)  f[i]=1; // user did not want to exclude olap 
}
for (IdxTy i=0; i<16; ++i) 
{
	const IdxTy qp=(i>>2)+qbase;
	const IdxTy q=(i&3)+qbase;
	if (f[i]!=0) 
		res+= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, qp, q,false) ;
}
break;
}

case (257):{
const D xf=xc+vxt;
const D yf=yc+vyt;
IdxTy qp=0;
IdxTy q=0;
if ((xf-xpc>a)&&(yf-ypc>b)) {qp=0; q=2;}
else if ((xf-xpc>a)&&(yf-ypc>0)) {qp=0; q=1;}
else if ((xf-xpc>a)&&(yf-ypc<-b)) {qp=3; q=1;}
else if ((xf-xpc>a)&&(yf-ypc<0)) {qp=3; q=2;}
else if ((xf-xpc>0)&&(yf-ypc>b)) {qp=0; q=2;}
else if ((xf-xpc>0)&&(yf-ypc>0)) {qp=0; q=1;}
else if ((xf-xpc>0)&&(yf-ypc<-b)) {qp=3; q=1;}
else if ((xf-xpc>0)&&(yf-ypc<0)) {qp=3; q=2;}
else if ((xf-xpc<(-a))&&(yf-ypc>b)) {qp=0; q=2;}
else if ((xf-xpc>(-a))&&(yf-ypc>0)) {qp=0; q=1;}
else if ((xf-xpc>(-a))&&(yf-ypc<-b)) {qp=3; q=1;}
else if ((xf-xpc>(-a))&&(yf-ypc<0)) {qp=3; q=2;}
else if (yf-ypc>b) {qp=0; q=2;}
else if (yf-ypc>0) {qp=0; q=1;}
else if (yf-ypc<-b) {qp=3; q=1;}
else  {qp=3; q=2;}

res+= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, qp+qbase, q+qbase,noverlap) ; 


break;
}
default: { MM_MSG(" bad quadrant pattern "<<pattern)  
 MM_ONCE(" bad quadrant pattern(s) see cout  "<<pattern,) } 
} // pattern 


//const D quads=16; 
// FIXME TODO 
//MM_ONCE(" only one quadrant being evaluated instead of 4x4",)
return res;

}

D normalize( const D & a, const D & b, const D & ainf, const D & binf 
, const D & beta, const bool noverlap=!true) 
{
D res=0;
// the best thing to do is eliminate the destination shape functions but this is 
// easier 
// get all the self terms
for (IdxTy q=0; q<4; ++q)
for (IdxTy qp=0; qp<4; ++qp)
res+=normalize(a,b,ainf,binf,beta,qp,q,4,noverlap);

for (IdxTy q=0; q<4; ++q)
{
	res+=normalize(a,b,ainf,binf,beta,1,q,0,noverlap);
	res+=normalize(a,b,ainf,binf,beta,2,q,1,noverlap);
	res+=normalize(a,b,ainf,binf,beta,3,q,2,noverlap);
	res+=normalize(a,b,ainf,binf,beta,4,q,3,noverlap);
}

return res;
}

// source regions axb received by ainf x binf 
D normalize( const D & a, const D & b, const D & ainf, const D & binf 
, const D & beta,const IdxTy quadp,const IdxTy quad, const IdxTy node, const bool noverlap=!true) 
{
const bool dumpep=false;
const bool xn=(quad==3)||(quad==2);
const bool xpn=(quadp==3)||(quadp==2);
const bool yn=(quad==3)||(quad==4);
const bool ypn=(quadp==3)||(quadp==4);
// these need to be either self and identical OR the other half
// centered on remove nodes at ainf and binf. 
D xpc_node=0;
D ypc_node=0;
switch (node)
{
// for really large values, the noise floor will make the results spuriously high 
case 0:{ xpc_node=-ainf; ypc_node=-binf; break; } 
case 1:{ xpc_node=ainf; ypc_node=-binf; break; } 
case 2:{ xpc_node=ainf; ypc_node=binf; break; } 
case 3:{ xpc_node=-ainf; ypc_node=binf; break; } 
// use 4 for cetner 
}
const D xpc=xpc_node;
const D ypc=ypc_node;
const D xc=0;
const D yc=0;
// the src quadrants have only one sf, the dest quadrants need both 
// which comes from the dest nodes and the srd node
const D xpu=xpc+(xpn?0:ainf);
const D xpl=xpc+(xpn?(-ainf):0);// -a;
const D xu=xc+(xn?0:a);
const D xl=xc+(xn?(-a):0); // -a;
const D ypu=ypc+(ypn?0:binf);
const D ypl=ypc+(ypn?(-binf):0); // -a;
const D yu=yc+(yn?0:b);
const D yl=yc+(yn?(-b):0); // -a;
if (noverlap)
{// this is where rationals would be nice 
const D xplen=(xpu-xpl);
const D xlen=(xu-xl);
const D yplen=(ypu-ypl);
const D ylen=(yu-yl);
const bool xolap=(xc-xpc)*(xc-xpc)<(xlen*xplen)*(.25-.01);
const bool yolap=(yc-ypc)*(yc-ypc)<(ylen*yplen)*(.25-.01);

if (xolap&&yolap) return 0;
}
D ix=  esimple(  beta,  xpu, xpl,  xu, xl ,0,(quadp==2)||(quadp==3),(quad==2)||(quad==3));
D iy= esimple(  beta,  ypu, ypl,  yu, yl ,0,(quadp==4)||(quadp==3),(quad==4)||(quad==3));
D res=ix*iy;
//if (dbug) { MM_MSG(MMPR3(ix,iy,res)) } 
return res;

} // normalie 
// evaluate the gaussin  distro near a fixed line source. 
static D evaluate_ls( const D & xpc, const D & ypc, const D & xc, const int dir 
, const D & vxt, const D & vyt, const D & a, const D & b, const D & beta,const IdxTy quadp) 
{
const bool xpn=(quadp==3)||(quadp==2);
const bool ypn=(quadp==3)||(quadp==4);
const D xpu=xpc+(xpn?0:a);
const D xpl=xpc+(xpn?(-a):0);// -a;
const D ypu=ypc+(ypn?0:b);
const D ypl=ypc+(ypn?(-b):0); // -a;

const bool xd=(dir==1)||(dir==3);
const bool yd=(dir==0)||(dir==2);
const bool nd= (dir==1)||(dir==2);
//D xdel=0; D ydel=0;
// this needs to check for overlap with the constant region on the "otherside"
// of the line actually a half plane. 
D res=0;
// const_line_source( beta,  _xpu, _xpl,  _x , vt,  qpneg)
// limiting the overlap here will mess up the shape functions which are made
// to span integration span. 
const D ix=xd?const_line_source( beta, xpu, xpl, xc ,vxt,(quadp==2)||(quadp==3),nd):a;
const D iy=yd?const_line_source( beta, ypu, ypl, xc ,vyt,(quadp==4)||(quadp==3),nd):b;
// need to see if it overlaps with the fixed area

if (false) res=-ix*iy; else res=ix*iy;
if (res<0)
{
//MM_MSG(" for ce zer res<0 "<<MMPR4(res,ix,iy,dir)<<MMPR4(beta,xpu,xpl,xc))
if (res<-1e-4) {MM_MSG(" for ce zer res<0 danger willr robinson "<<MMPR4(res,ix,iy,dir)<<MMPR4(beta,xpu,xpl,xc))}
MM_ONCE(" for ce zer res<0 "<<MMPR4(res,ix,iy,dir)<<MMPR4(beta,xpu,xpl,xc),)
res=0;
}
// the stuff on the otherside is in the const bulk however this stull needs
// to work with shape functions
//res+=xdel*b+ydel*a;
return res;
} // evaluate_ls

//////////////////////////////////////////////////////////////////////
// evaluate transfer from node at (x,y) to node at (xp,yp)
// integrating over a square of side 2a at each location.
// the integral is \int c(xp)dxp \int  s(x)exp(-beta*(xp-x-vt) dx
// over a square of size 4a^2 at both locations. c(xp) and s(x) are
// the linear shape functions. This resolves into a product
// but there are some numeric issues with large x 
// the normalization is ignored in the app right now
// but there are still quality issues at larges beta*x as
// erf goes to 1. 
//D res= evaluateq( xpc, ypc, xc, yc, vxt, vyt, a, b,beta, quadp, quad) ;
static D evaluateq( const D & xpc, const D & ypc, const D & xc, const D & yc
, const D & vxt, const D & vyt, const D & a, const D & b, const D & beta,const IdxTy quadp,const IdxTy quad, const bool noverlap=true) 
{
const bool dumpep=false;
const bool dbug=false;
if (dumpep) { MM_MSG(MMPR2(quadp,quad)<<MMPR4(xpc,ypc,xc,yc)<<MMPR4(vxt,vyt,a,beta)) }
D res=0;
// the important thing is teh DIFERENCE between flowed smearer distro
// and \delta(x-x') distrobution which exists with overlap 
const bool subtract_olap=true;
bool have_olap=false;
// this is ZERO based wtf 
const bool xn=(quad==3)||(quad==2);
const bool xpn=(quadp==3)||(quadp==2);
const bool yn=(quad==3)||(quad==4);
const bool ypn=(quadp==3)||(quadp==4);



// the absolute value is a big fudd. Need to thinkhere doh 
// need abs values etc this is a mess lol 
const D xpu=xpc+(xpn?0:a);
const D xpl=xpc+(xpn?(-a):0);// -a;
const D xu=xc+(xn?0:a);
const D xl=xc+(xn?(-a):0); // -a;
const D ypu=ypc+(ypn?0:b);
const D ypl=ypc+(ypn?(-b):0); // -a;
const D yu=yc+(yn?0:b);
const D yl=yc+(yn?(-b):0); // -a;

// do not do overlapping elelemtns
if (noverlap||subtract_olap)
{// this is where rationals would be nice 
const D xplen=(xpu-xpl);
const D xlen=(xu-xl);
const D yplen=(ypu-ypl);
const D ylen=(yu-yl);
// these need to eliminate the case where the boundaries over lap 
//const bool xoolap=(xc<xpu)&(xc>xpl);
//const bool yoolap=(yc<ypu)&(yc>ypl);
//const bool xolap=((xc-xpc)*(xc-xpc))<((xlen*xplen)*(.25-.01));
const bool xolap=((xc-xpc)*(xc-xpc))<((xlen*xplen)*(1.01));
//const bool yolap=((yc-ypc)*(yc-ypc))<((ylen*yplen)*(.25-.01));
const bool yolap=((yc-ypc)*(yc-ypc))<((ylen*yplen)*(1.01));
//const bool qolap=(quad==1)&&(quadp==2)&&(xpc<=xu)&&(xpc>=xc)&&(ypc>=yc)&&(ypc<=yu)
//				||(quad==2)&&(quadp==1)&&(xc<=xpu)&&(xc>=xpc)&&(yc>=ypc)&&(yc<=ypu)
//				||(quad==4)&&(quadp==3)&&(xc<=xpu)&&(xc>=xpc)&&(yc<=ypc)&&(yc>=ypl)

have_olap=(xolap&&yolap); //  return 0;
//MM_MSG("overalsps "<<MMPR4(xolap,yolap,have_olap,subtract_olap))
//if (noverlap) if (xolap&&yolap) return 0;
if (noverlap&&have_olap) return 0;
//if (qolap) return 0;
}

/*
// these acutally need to be abs values and split even more
// also the scale is probably wrong, need 1/a, 1/b etc 
D s1x=-1.0/a;
D s1y=-1.0/b;
const D s0=1;
D c1x=-1.0/a;
D c1y=-1.0/b;
const D c0=1;
if (xn ) { s1x=-s1x;}
if (yn ) { s1y=-s1y;}
if (xpn ) { c1x=-c1x;}
if (ypn ) { c1y=-c1y;}
*/

// fudd, c is the SOURCE and s is the DEST wtf 
// the code should do this if sf is rght 
//s1x=s1y=c1x=c1y=0;
//const D s0x= s0 -s1x*(xpc+0*vxt);
//const D s0y= s0 -s1y*(ypc+0*vyt);
// need to check this 
//const D c0x=c0 -c1x*(xc-vxt);
//const D c0y= c0 -c1y*(yc-vyt);

// the sf origins have to be the center point even when negative now 
// so now the sf is ok but what about the limits lol wtf 
//const D ff=0; // 0; // 1.0;
// translate to make the c shape function right, now just fix s
//const D x1=ff*xpc;
//const D x2=ff*xc;
//const D y1=ff*ypc;
//const D y2=ff*yc;
// change variables on the inner integral leading to limits that are first
// order polynomials ( axp+b) in the outer variable .
//if (dbug) { MM_MSG(MMPR4(beta,vxt,x1,x2)<<MMPR4(xpu,xpl,xu,xl)<<MMPR4(s1x,s0,c1x,c0)) } 
//if (dbug) { MM_MSG(MMPR4(beta,vyt,y1,y2)<<MMPR4(ypu,ypl,yu,yl)<<MMPR4(s1y,s0,c1y,c0)) } 

//D esimple(  beta,  xpu, xpl,  xu, xl ,vt ) 
const bool use_simple=true;
if (use_simple)
{
/*
D ix00= esimple(  beta,  xpu, xpl,  xu, xl ,vxt,0,1,0,1 );
D ix10= esimple(  beta,  xpu, xpl,  xu, xl ,vxt,1,0,0,1 );
D ix11= esimple(  beta,  xpu, xpl,  xu, xl ,vxt,0,1,1,0 );
D ix01= esimple(  beta,  xpu, xpl,  xu, xl ,vxt,1,0,1,0 );
 
D iy00= esimple(  beta,  ypu, ypl,  yu, yl ,vyt,0,1,0,1 );
D iy10= esimple(  beta,  ypu, ypl,  yu, yl ,vyt,1,0,0,1 );
D iy01= esimple(  beta,  ypu, ypl,  yu, yl ,vyt,0,1,1,0 );
D iy11= esimple(  beta,  ypu, ypl,  yu, yl ,vyt,1,0,1,0 );


if (dbug) { MM_MSG(MMPR4(ix00,ix10,ix01,ix11)<<MMPR4(iy00,iy10,iy01,iy11)) } 
res+=iy00*(ix00-ix10-ix01-ix11);
res-=iy10*(ix00-ix10-ix01-ix11);
res-=iy01*(ix00-ix10-ix01-ix11);
res-=iy11*(ix00-ix10-ix01-ix11);
*/
//MM_ONCE(" dropping ix for now",)
//D ix= 1; // esimple(  beta,  xpu, xpl,  xu, xl ,vxt,1,1,1,1 );
D ix=0;
D iy=0;
if (false)
{
D ix=  esimple(  beta,  xpu, xpl,  xu, xl ,vxt,(quadp==2)||(quadp==3),(quad==2)||(quad==3));
D iy= esimple(  beta,  ypu, ypl,  yu, yl ,vyt,(quadp==4)||(quadp==3),(quad==4)||(quad==3));
}else{

//D ix=  fl_mb_diff_poly(  beta,  xpu, xpl,  xu, xl ,vxt,(quadp==2)||(quadp==3),(quad==2)||(quad==3));
//D iy= fl_mb_diff_poly(  beta,  ypu, ypl,  yu, yl ,vyt,(quadp==4)||(quadp==3),(quad==4)||(quad==3));

D ix=  fl_mb_diff_comp2(  beta,  xpu, xpl,  xu, xl ,vxt,(quadp==2)||(quadp==3),(quad==2)||(quad==3));
D iy= fl_mb_diff_comp2(  beta,  ypu, ypl,  yu, yl ,vyt,(quadp==4)||(quadp==3),(quad==4)||(quad==3));


res=ix*iy;
return res; 
}

if (!have_olap) res=ix*iy;
else
{
if (subtract_olap)
{
const bool xsf2=(xn&&xpn)||(!xn&&!xpn);
const bool ysf2=(yn&&ypn)||(!yn&&!ypn);
const D norm=sqrt(beta/M_PI);
const D ixzed=.5*(xsf2?(a/3.0):(a/6.0));
//const D ixzed=xsf2?(1.0/3.0):(1.0/6.0);
const D iyzed=.5*(ysf2?(b/3.0):(b/6.0));
//const D iyzed=ysf2?(1.0/3.0):(1.0/6.0);
const D ixz=ix*norm-ixzed;
const D iyz=iy*norm-iyzed;
res=(ixz)*(iyz);
MM_MSG(" have overlap "<<MMPR3(norm,(ix*norm),(iy*norm))<< MMPR4(ix,iy,ixzed,iyzed)<<MMPR4(ixz,iyz,res,quad)<<MMPR4(quadp,beta,vyt,xpu)<<MMPR4(xpl,xu,xl,ypu)<<MMPR3(ypl,yu,yl)<<MMPR4(xc,yc,xpc,ypc))

if ((res<0) ||(ixz<0)||(iy<0))
{
MM_MSG( MMPR4(ix,iy,ixzed,iyzed)<<MMPR4(ixz,iyz,res,quad)<<MMPR4(quadp,beta,vyt,xpu)<<MMPR4(xpl,xu,xl,ypu)<<MMPR3(ypl,yu,yl))

res=0; 
}
} // subtract 

} // have 
if (dbug) { MM_MSG(MMPR3(ix,iy,res)) } 
}
else
{
/*
//D evaluate1( beta, xpu,xpl, xu, xl , vxt ,  s1, s0, c1,  c0, xpbase=0, xbase=0) 
D ix00=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  0,  s0x,  0,  c0x,x1,x2);
D ix10=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1x,  0,  0,  c0x,x1,x2);
D ix01=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  0,  s0x,  c1x,  0,x1,x2);
D ix11=evaluate1(  beta,  xpu, xpl, xu , xl , vxt ,  s1x,  0,  c1x,  0,x1,x2);

D iy00=evaluate1(  beta,  ypu, ypl, yu , yl , vyt ,  0,  s0y,  0,  c0y,y1,y2);
D iy10=evaluate1(  beta,  ypu, ypl, yu , yl , vyt ,  s1y,  0,  0,  c0y,y1,y2);
D iy01=evaluate1(  beta,  ypu, ypl, yu , yl , vyt ,  0,  s0y,  c1y,  0,y1,y2);
D iy11=evaluate1(  beta,  ypu, ypl, yu , yl , vyt ,  s1y,  0,  c1y,  0,y1,y2);

if (dbug) { MM_MSG(MMPR4(ix00,ix10,ix01,ix11)<<MMPR4(iy00,iy10,iy01,iy11)) } 
res+=iy00*(ix00-ix10-ix01-ix11);
res-=iy10*(ix00-ix10-ix01-ix11);
res-=iy01*(ix00-ix10-ix01-ix11);
res-=iy11*(ix00-ix10-ix01-ix11);
*/
}



return res; 
}

////////////////////////////////////////////////////////////////////

// added for doing overlaps in meaninful way 
// flow by spreading out exp(-k*(x-x')^2) given initial distro in ni
// spacing between nodes is l 
template <class Ty> void do_brute_force_sim(const D & k, const D & l, const Ty & ni, const IdxTy n )
{
Ty res(n),zed(n);
const IdxTy points=ni.size();
const D di=(points-1)*l/(n-1); 
const D A=sqrt(k/M_PI);
// for very large k, this will oscillate as it beats with the other points
const IdxTy nsf=n+100; // 100;
const D dnsf=l/D(1.0*nsf);
for (IdxTy i=0; i<n; ++i) // for each dest point
{
const D x=D(i)*di;
for (IdxTy j=0; j<points; ++j) // each src point 
{
const D xpzed=D(j)*l;
const D dzed=x-xpzed;
const D intf=di;
// this is the total integrated charge at each point summing over the 
// source shape functions 
if (dzed>=0) if (dzed<l) zed[i]+=intf*ni[j]*(l-dzed)/l;
if (dzed<0) if (dzed>-l) zed[i]+=intf*ni[j]*(l+dzed)/l;

for (IdxTy kk=0; kk<nsf; ++kk)
{
const D doff=(D(kk)+.5)*dnsf;
const D xpp=xpzed+doff;
const D xpn=xpzed-doff;
const D dijp=x-xpp;
const D dijn=x-xpn;
// may be a factor of 2 here, but symptotes seem right 
//if ((dij>=0)&&(dij<l)) 
// this is the resulting charge distribution convolved 
res[i]+=di*ni[j]*exp(-k*dijn*dijn)*(l-doff)/l*A*dnsf;
//if ((dij<=0)&&(dij>-l)) 
if (doff!=0) res[i]+=di*ni[j]*exp(-k*dijp*dijp)*(l-doff)/l*A*dnsf;
//MM_MSG(MMPR3(i,j,kk)<<MMPR4(res[i],di,ni[j],k)<<MMPR4(dijp,dijn,l,doff)<<MMPR2(A,dnsf))
}



}


}
for (IdxTy i=0; i<n; ++i) MM_MSG(" bf_flow "<<MMPR4(i,i*di,res[i],zed[i])) 

for (IdxTy j=0; j<points; ++j) // each src point 
{
const D xpzed=D(j)*l;
D sump=0;
D ressp=0;
D zeddp=0;
D sumn=0;
D ressn=0;
D zeddn=0;

const D dif=1;
for (IdxTy i=0; i<n; ++i) // for each dest point
{
const D x=D(i)*di;
const D doff=x-xpzed;
// this now inegrates the original, convolved, and net over the dest shape functions 
if (doff>=0) if (doff<l) sump+=dif*(res[i]-zed[i])*(l-doff)/l;
if (doff<0) if (doff>-l) sumn+=dif*(res[i]-zed[i])*(l+doff)/l;
if (doff>=0) if (doff<l) ressp+=dif*(res[i])*(l-doff)/l;
if (doff<0) if (doff>-l) ressn+=dif*(res[i])*(l+doff)/l;
if (doff>=0) if (doff<l) zeddp+=dif*(zed[i])*(l-doff)/l;
if (doff<0) if (doff>-l) zeddn+=dif*(zed[i])*(l+doff)/l;


} // i 
D sum=sump+sumn;
D ress=ressp+ressn;
D zedd=zeddp+zeddn;
MM_ERR(" sum "<<MMPR4(j,sum,ress,zedd)<<MMPR4(sump,sumn,ressp,ressn))
MM_MSG(" sum "<<MMPR4(j,sum,ress,zedd)<<MMPR4(sump,sumn,ressp,ressn))
//MM_MSG(" sum "<<MMPR4(j,sum,ress,zedd))
} // j 


}

/////////////////////////////////////////////////////////////////////
private:
void Init()
{
//m_ord_col=m_flp.ordinal_column();
//m_t_col=m_flp.time_column();
m_done=false;

}
bool m_done;
//BraTy m_tree;
//FlpTy m_flp; // right now this is being SET in Init not being used to set m_h

//TimeLines m_timelines;
//KeyCounts m_unused;

//MyChart m_chart;
//time_map m_time_map;
//IdxTy m_ord_col; // input column containing a time ordinal (such as day number)
//IdxTy m_t_col; // input col with absolute time string such as date.
Logic m_flp;

CounterMap m_cm;

}; //mjm_timeline 



/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
#ifdef  TEST_BF_FLOW__
int main(int argc,char **args)
{
typedef mjm_gauss  Myt;
typedef double D;
typedef unsigned int IdxTy;
Myt x(argc,args);
const D k=atof(args[1]);
const D l=1;
const IdxTy npts=6;
const unsigned int n=atoi(args[2]); // 100;
std::vector<D> ni(npts);
MM_MSG(MMPR(k))
ni[npts/2]=1;
x.do_brute_force_sim( k,  l,  ni,  n );


//if (!x.done()) x.command_mode();

return 0;
}

#endif


#ifdef  TEST_GAUSS__

int main(int argc,char **args)
{
typedef mjm_gauss  Myt;

Myt x(argc,args);

if (!x.done()) x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

