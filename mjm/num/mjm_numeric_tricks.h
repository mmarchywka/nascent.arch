#ifndef MJM_NUMERIC_TRICKS_H__
#define MJM_NUMERIC_TRICKS_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>

/*
*/

/*
g++ -DTEST_NUMERIC_TRICKS__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_numeric_tricks.h
*/

////////////////////////////////////////////////////////////////



class mjm_numeric_tricks 
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

typedef mjm_numeric_tricks Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;




public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
mjm_numeric_tricks() {Init();}
mjm_numeric_tricks(int argc,char **_args) // : m_save_sol(false)
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<argc; ++i) args[i]=_args[i];
for (IdxTy i=argc; i<ikluge; ++i) args[i]=&dummy[0];

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
if ( s==StrTy("-vssv")) { arg_cmd(i,args,0,"vssv",confirm); return; }
if ( s==StrTy("-hssv")) { arg_cmd(i,args,0,"hssv",confirm); return; }
if ( s==StrTy("-init")) { arg_cmd(i,args,0,"init",confirm); return; }
if ( s==StrTy("-quit")) { arg_cmd(i,args,0,"quit",confirm); return; }
if ( s==StrTy("-add")) { arg_cmd(i,args,1,"add",confirm); return; }
// tl and event
if ( s==StrTy("-add-event")) { arg_cmd(i,args,2,"add-event",confirm); return; }
if ( s==StrTy("-add-sum")) { arg_cmd(i,args,1,"add-sum",confirm); return; }
if ( s==StrTy("-set-sum")) { arg_cmd(i,args,1,"set-sum",confirm); return; }
if ( s==StrTy("-reset-sum")) { arg_cmd(i,args,1,"reset-sum",confirm); return; }
if ( s==StrTy("-make")) { arg_cmd(i,args,0,"make",confirm); return; }
if ( s==StrTy("-unused")) { arg_cmd(i,args,0,"unused",confirm); return; }
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


void dump_unused()
{
Ss ss;
//for (auto ii=m_unused.begin(); ii!=m_unused.end(); ++ii)
{
//ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
MM_MSG("unused:"<<CRLF<<ss.str())

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
MM_ERR(" processing "<<li.dump())
if (sz<1) continue;
const StrTy cmd=li.word(0);
//if (cmd=="solve") { solve(); continue; } 
//if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
//if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
//if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

//if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
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
if (cmd=="init") { Init();   continue; } 
//if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
//if (cmd=="banner") { config_banner();  continue; } 
//if (cmd=="cm") { dump_cm();  continue; } 
/*
if (cmd=="test") { test(li);  continue; } 
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
*/
if (cmd=="quit") { clean_up(); return; } 

MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())

}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
//void test_gr_integrator( CommandInterpretter & li )
MM_ERR(" unrecignized TEST command "<<li.dump())

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
//////////////////////////////////////////////////////////////////////
// compute ca*exp(a)-cb*exp(b) with better precision when the 
// operands are mcch gretaer than the difference . Uses repeated
// difference of two squares to turn difference into a factor
// this now computes ca*exp(a)-cb*exp(b)
// where we expect ca and cb close to unity and a,b near zero and negative
// designed for case where maybe a complementary normal distribution would
// make sense 
// copied from sqdiff.cpp 2017-05-03
//
/*
Copied from expdiff.cpp. Try to find the difference
between two similar terms near 1 with difference of
square dilation of difference

(a^(2^n)-b^(2^n))=(a-b) \pi (a^(2^i)+b^(2^i)) , 0<i<=n

*/

static D sqdiffe(const D & a, const D & b, const D & ca, const D & cb, const IdxTy n)
{
D den=1.0;
D an=a;
D can=ca;
D bn=b;
D cbn=cb;
const IdxTy neff=n;
for (IdxTy i=0; i<neff; ++i)
{
if (an<-50) break;
if (bn<-50) break;
// this is essentially 2^s where s= i+(i-1)+(i-2) ... 
if (true)
{D temp=  (can*exp(an)-cbn*exp(bn))/den;
MM_MSG(MMPR(i)<<MMPR(den)<<MMPR(can)<<MMPR(an)<<MMPR(cbn)<<MMPR(bn)<<MMPR(temp))
}
den=den*(can*exp(an)+cbn*exp(bn));
// these are both about 1
an=an*2;
bn=bn*2;
can=can*can;
cbn=cbn*cbn;
// arbitrary cutoff for now.  This should not go <0 for now
if (den>1e12) break; 
if (den<1e-12) break; 
}
return (can*exp(an)-cbn*exp(bn))/den;

}

// computes ca*exp(a)-cb*exp(b) when the difference is much smaller than one
// works better than the above usually in 1 term for the cases of interest. 
static D sqdiffetaylor(const D & a, const D & b, const D & ca, const D & cb, const IdxTy n)
{
if (n==0)  return ca*exp(a)-cb*exp(b);
D sum=0;
D xn=1;
D yn=1;
D den=1;
// really, an array makes more sense here i this is time criticla code 
std::vector<D> terms;
terms.push_back((ca*xn-cb*yn)/den);
for(IdxTy i=1; i<n; ++i)
{
xn=xn*a;
yn=yn*b;
den=den*i;
terms.push_back((ca*xn-cb*yn)/den);

}
while (terms.size()!=0) {sum=sum+terms.back(); terms.pop_back(); }
return sum;
}
// as above for erf, ca*erf(a)-cb*erf(b)
static D sqdifferftaylor(const D & a, const D & b, const D & ca, const D & cb, const IdxTy n)
{
if (n==0)  return ca*erf(a)-cb*erf(b);
D sum=0;
D xn=a;
D yn=b;
D den=1;
// really, an array makes more sense here i this is time criticla code 
std::vector<D> terms;
terms.push_back((ca*xn-cb*yn)/den);
for(IdxTy i=1; i<n; ++i)
{
xn=xn*a*a;
yn=yn*b*b;
den=den*i;
const D sign=((i&1)==0)?1:-1;
terms.push_back(sign*(ca*xn-cb*yn)/((2*i+1)*den));

}
while (terms.size()!=0) {sum=sum+terms.back(); terms.pop_back(); }
static const D fac=.5/sqrt(M_PI);
return sum*fac;
}
#if 0 
// thia does not work 
// term i of taylor series integrated twice with n poly term
// prefix. 
// n is the power of rht prefeactor, x^n * exp(- beta * y^2) 
// nint is the number of integrations performed, thought to be 2 before
static D mb_taylor(const D & beta,const D & yff,const D & yii,
	const D & yfi,const D & yif,const IdxTy n,const IdxTy nterms, const IdxTy nint)
{
D sum=0;
D yff2=yff*yff; D yii2=yii*yii; D yfi2=yfi*yfi; D yif2=yif*yif;
//D zff=yff2; D zii=yii2; D zfi=yfi2; D zif=yif2;
// include zero but it drops out exactly 
D zff=1; D zii=1; D zfi=1; D zif=1;
// this adds in the terms for hte polynomial prefactor 
for (IdxTy i=0; i<n; ++i)
 {zff=yff*zff; zii=yii*zii; zfi=yfi*zfi; zif=yif*zif;}
// integrations raise the squared thing too 
for (IdxTy i=0; i<nint; ++i)
 {zff=yff2*zff; zii=yii2*zii; zfi=yfi2*zfi; zif=yif2*zif;}
D fac=1.0;
// this better be zero lol  for poly just one 
sum+=fac/((n+1)*(n+2))*(zff+zii-zfi-zif);
for (IdxTy i=1; i<nterms; ++i)
{
fac=fac*(-beta)/D(i);
 zff=yff2*zff; zii=yii2*zii; zfi=yfi2*zfi; zif=yif2*zif;
sum+=fac/((2*i+n+1)*(2*i+n+2))*(zff+zii-zfi-zif);
}
return sum;
}
#endif

// rederive in parking lot, see notes for terminaology.
// specialized to g(v)=exp(-beta*y^2) with polynomial
// (s1*x' + s0 ) \int_x (x1*x+c0)*(ax+b) g(v) dx' dx
// and requiring s0 differ for xi and xf due to change of var
// p_i are the net coeefs of y^2, y^1,y^0 after change of var 
// I forgot the p_i have x and y terms fudd 
// (beta,yff,yii,yfi,yif,s1,syf,syi,p2,p1f,p1i,p1y,p0f,p0i,poyf,p0yi,p0y2,10);
//(beta,yff,yii,yfi,yif,s1,syf,syi,p2,p1f,p1i,p1y,p0f,p0i,p0yf,p0yi,p0y2,10);
static D mb_taylor_old(const D & beta,const D & yff,const D & yii,
	const D & yfi,const D & yif,const D & s1, const D & s0f,const D & s0i
,const D & p2, const D & p1f, const D & p1i, const D & p1y
,const D & p0f, const D & p0i
, const D & p0yf, const D p0yi, const D & p0y2 
,  const IdxTy nterms)
{
D sum=0;
D yff2=yff*yff; D yii2=yii*yii; D yfi2=yfi*yfi; D yif2=yif*yif;
D zff=1; D zii=1; D zfi=1; D zif=1;
 zff=zff*yff2;  zii=zii*yii2;  zfi=zfi*yfi2;  zif=zif*yif2;
D fac=1.0;
D zffj=zff; D ziij=zii; D zfij=zfi; D zifj=zif;
MM_MSG(MMPR(yff)<<MMPR(yii)<<MMPR(yfi)<<MMPR(yif))
for (IdxTy j=0; j<3; ++j)
{
if (j==2)
{
D wij=p2/(D(j+1)*D(j+3));
D zij=p2/((D(j+1)*(D(j+2))));
// this needs to change in j loop and then revert.. fudd 
D sumy=s1*(zffj*yff+ziij*yii-zfij*yfi-zifj*yif);
// NB get fudding subscripts in right order 
D sumyp=(zffj*s0f+ziij*s0i-zfij*s0i-zifj*s0f);
sum+=wij*sumy+zij*sumyp;
MM_MSG("j==2 "<<MMPR(sum)<<MMPR(wij)<<MMPR(sumy)<<MMPR(zij)<<MMPR(sumyp))
}
// arrrghhh
if (j==1)
{
D wij1=1.0/(D(j+1)*D(j+3));
D zij1=1.0/((D(j+1)*(D(j+2))));

D sumy1=s1*(zffj*yff*p1f+ziij*yii*p1i-zfij*yfi*p1i-zifj*yif*p1f);
D sumyp1=(zffj*s0f*p1f+ziij*s0i*p1i-zfij*s0i*p1i-zifj*s0f*p1f);
sum+=wij1*sumy1+zij1*sumyp1;
// then the p1y term
D wij2=p1y/(D(j+1)*D(j+4));
D zij2=p1y/((D(j+1)*(D(j+3))));
D sumy2=s1*(zffj*yff2+ziij*yii2-zfij*yfi2-zifj*yif2);
D sumyp2=(zffj*s0f*yff+ziij*s0i*yii-zfij*s0i*yfi-zifj*s0f*yif);
sum+=wij2*sumy2+zij2*sumyp2;

MM_MSG("j==1 "<<MMPR(sum)<<MMPR(wij1)<<MMPR(sumy1)<<MMPR(zij1)<<MMPR(sumyp1)<<MMPR(wij2)<<MMPR(sumy2)<<MMPR(zij2)<<MMPR(sumyp2))
} // j==1
if (j==0)
{

D wij0=1.0/(D(j+1)*D(j+3));
D zij0=1.0/((D(j+1)*(D(j+2))));

D sumy0=s1*(zffj*yff*p0f+ziij*yii*p0i-zfij*yfi*p0i-zifj*yif*p0f);
D sumyp0=(zffj*s0f*p0f+ziij*s0i*p0i-zfij*s0i*p0i-zifj*s0f*p0f);
sum+=wij0*sumy0+zij0*sumyp0;


D wij2=1.0/(D(j+1)*D(j+4));
D zij2=1.0/((D(j+1)*(D(j+3))));
D sumy2=s1*(zffj*yff2*p0yf+ziij*yii2*p0yi-zfij*yfi2*p0yi-zifj*yif2*p0yf);
D sumyp2=(zffj*s0f*yff*p0yf+ziij*s0i*yii*p0yi-zfij*s0i*yfi*p0yi-zifj*s0f*yif*p0yf);
sum+=wij2*sumy2+zij2*sumyp2;

D wij22=p0y2/(D(j+1)*D(j+5));
D zij22=p0y2/((D(j+1)*(D(j+4))));
// this needs to change in j loop and then revert.. fudd 
D sumy22=s1*(zffj*yff*yff2+ziij*yii*yii2-zfij*yfi*yfi2-zifj*yif*yif2);
// NB get fudding subscripts in right order 
D sumyp22=(zffj*s0f*yff2+ziij*s0i*yii2-zfij*s0i*yfi2-zifj*s0f*yif2);
sum+=wij22*sumy22+zij22*sumyp22;

MM_MSG("j==0 "<<MMPR(sum)<<MMPR(wij0)<<MMPR(sumy0)<<MMPR(zij0)<<MMPR(sumyp0) <<MMPR(wij2)<<MMPR(sumy2)<<MMPR(zij2)<<MMPR(sumyp2)<<MMPR(wij22)<<MMPR(sumy22)<<MMPR(zij22)<<MMPR(sumyp22))



} // j=0=


MM_MSG(MMPR(j)<<MMPR(sum))
 zffj*=yff;  ziij*=yii;  zfij*=yfi; zifj*=yif;
}

for (IdxTy i=1; i<nterms; ++i)
{
fac=fac*(-beta)/D(i);
 zff=zff*yff2;  zii=zii*yii2;  zfi=zfi*yfi2;  zif=zif*yif2;
 zffj=zff;  ziij=zii;  zfij=zfi; zifj=zif;
for (IdxTy j=0; j<3; ++j)
{
if (j==2)
{
D wij=p2/(D(2*i+j+1)*D(2*i+j+3));
D zij=p2/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy=s1*(zffj*yff+ziij*yii-zfij*yfi-zifj*yif);
D sumyp=(zffj*s0f+ziij*s0i-zfij*s0i-zifj*s0f);
sum+=fac*(wij*sumy+zij*sumyp);
}
if (j==1)
{

D wij1=1.0/(D(2*i+j+1)*D(2*i+j+3));
D zij1=1.0/((D(2*i+j+1)*(D(2*i+j+2))));

D sumy1=s1*(zffj*yff*p1f+ziij*yii*p1i-zfij*yfi*p1i-zifj*yif*p1f);
D sumyp1=(zffj*s0f*p1f+ziij*s0i*p1i-zfij*s0i*p1i-zifj*s0f*p1f);
sum+=fac*(wij1*sumy1+zij1*sumyp1);
// then the p1y term
D wij2=p1y/(D(2*i+j+1)*D(2*i+j+4));
D zij2=p1y/((D(2*i+j+1)*(D(2*i+j+3))));
D sumy2=s1*(zffj*yff2+ziij*yii2-zfij*yfi2-zifj*yif2);
D sumyp2=(zffj*s0f*yff+ziij*s0i*yii-zfij*s0i*yfi-zifj*s0f*yif);
sum+=fac*(wij2*sumy2+zij2*sumyp2);

} // j==1
if (j==0)
{

D wij0=1.0/(D(2*i+j+1)*D(2*i+j+3));
D zij0=1.0/((D(2*i+j+1)*(D(2*i+j+2))));

D sumy0=s1*(zffj*yff*p0f+ziij*yii*p0i-zfij*yfi*p0i-zifj*yif*p0f);
D sumyp0=(zffj*s0f*p0f+ziij*s0i*p0i-zfij*s0i*p0i-zifj*s0f*p0f);
sum+=fac*(wij0*sumy0+zij0*sumyp0);


D wij2=1.0/(D(2*i+j+1)*D(2*i+j+4));
D zij2=1.0/((D(2*i+j+1)*(D(2*i+j+3))));
D sumy2=s1*(zffj*yff2*p0yf+ziij*yii2*p0yi-zfij*yfi2*p0yi-zifj*yif2*p0yf);
D sumyp2=(zffj*s0f*yff*p0yf+ziij*s0i*yii*p0yi-zfij*s0i*yfi*p0yi-zifj*s0f*yif*p0yf);
sum+=fac*(wij2*sumy2+zij2*sumyp2);

D wij22=p0y2/(D(2*i+j+1)*D(2*i+j+5));
D zij22=p0y2/((D(2*i+j+1)*(D(2*i+j+4))));
// this needs to change in j loop and then revert.. fudd 
D sumy22=s1*(zffj*yff*yff2+ziij*yii*yii2-zfij*yfi*yfi2-zifj*yif*yif2);
// NB get fudding subscripts in right order 
D sumyp22=(zffj*s0f*yff2+ziij*s0i*yii2-zfij*s0i*yfi2-zifj*s0f*yif2);

sum+=fac*(wij22*sumy22+zij22*sumyp22);
} // j=0=

MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(sum))


 zffj*=yff;  ziij*=yii;  zfij*=yfi; zifj*=yif;
}

}

return sum;
}
/////////////////////////////////////////////////////////////////////

class ztrack
{
typedef double D;
public:
ztrack(const D & yff,const D & yii,const D & yfi,const D & yif)
:m_yff(yff),m_yii(yii),m_yfi(yfi),m_yif(yif)
,m_yff2(yff*yff),m_yii2(yii*yii),m_yfi2(yfi*yfi),m_yif2(yif*yif)
, m_zff(m_yff2), m_zii(m_yii2), m_zfi(m_yfi2),m_zif(m_yif2)
{

}

void zj()
{ m_zffj=m_zff;  m_ziij=m_zii;  m_zfij=m_zfi; m_zifj=m_zif; }
void next_j() 
{
 m_zffj*=m_yff;  m_ziij*=m_yii;  m_zfij*=m_yfi; m_zifj*=m_yif;
}
void next_i(){ // pre ince
 m_zff=m_zff*m_yff2;  m_zii=m_zii*m_yii2;  m_zfi=m_zfi*m_yfi2;  m_zif=m_zif*m_yif2;
}

 D zj_sum(const D & cf, const D & ci) const
{
D sumyp=(m_zffj*cf+m_ziij*ci-m_zfij*ci-m_zifj*cf);
return sumyp;
}

D zjy_sum() const
{
D sumy=(m_zffj*m_yff+m_ziij*m_yii-m_zfij*m_yfi-m_zifj*m_yif);
return sumy;
}
D zjy2_sum() const
{
D sumy=(m_zffj*m_yff2+m_ziij*m_yii2-m_zfij*m_yfi2-m_zifj*m_yif2);
return sumy;
}

D zjy3_sum() const
{
D sumy=(m_zffj*m_yff2*m_yff+m_ziij*m_yii2*m_yii-m_zfij*m_yfi2*m_yfi-m_zifj*m_yif2*m_yif);
return sumy;
}


D zjy_sum(const D & cf, const D & ci) const
{
D sumy=(m_zffj*m_yff*cf+m_ziij*m_yii*ci-m_zfij*m_yfi*ci-m_zifj*m_yif*cf);
return sumy;
}
D zjy2_sum(const D & cf, const D & ci) const
{
D sumy=(m_zffj*m_yff2*cf+m_ziij*m_yii2*ci-m_zfij*m_yfi2*ci-m_zifj*m_yif2*cf);
return sumy;
}




 const D  m_yff, m_yii, m_yfi,m_yif;
 const D  m_yff2, m_yii2, m_yfi2,m_yif2;
 D  m_zff, m_zii, m_zfi,m_zif;
 D  m_zffj, m_ziij, m_zfij,m_zifj;

}; // ztrack

static D mb_taylor(const D & beta,const D & yff,const D & yii,
	const D & yfi,const D & yif,const D & s1, const D & s0f,const D & s0i
,const D & p2, const D & p1f, const D & p1i, const D & p1y
,const D & p0f, const D & p0i
, const D & p0yf, const D p0yi, const D & p0y2 
,  const IdxTy nterms, bool & failed)
{
D sum=0;
ztrack zt( yff,yii,yfi,yif);
D fac=1.0;
zt.zj();
// MM_MSG(MMPR(yff)<<MMPR(yii)<<MMPR(yfi)<<MMPR(yif))
#if 0
for (IdxTy j=0; j<3; ++j)
{
if (j==2)
{
D wij=p2/(D(j+1)*D(j+3)); D zij=p2/((D(j+1)*(D(j+2))));
D sumy=s1*zt.zjy_sum(); D sumyp=zt.zj_sum(s0f,s0i); 
sum+=wij*sumy+zij*sumyp;
//MM_MSG("j==2 "<<MMPR(sum)<<MMPR(wij)<<MMPR(sumy)<<MMPR(zij)<<MMPR(sumyp))
}
if (j==1)
{
D wij1=1.0/(D(j+1)*D(j+3)); D zij1=1.0/((D(j+1)*(D(j+2))));
D sumy1=s1*zt.zjy_sum(p1f,p1i); D sumyp1=zt.zj_sum(s0f*p1f,s0i*p1i); 
sum+=wij1*sumy1+zij1*sumyp1;
// then the p1y term
D wij2=p1y/(D(j+1)*D(j+4)); D zij2=p1y/((D(j+1)*(D(j+3))));
D sumy2=s1*zt.zjy2_sum(); D sumyp2=zt.zjy_sum(s0f,s0i); 
sum+=wij2*sumy2+zij2*sumyp2;
//MM_MSG("j==1 "<<MMPR(sum)<<MMPR(wij1)<<MMPR(sumy1)<<MMPR(zij1)<<MMPR(sumyp1)<<MMPR(wij2)<<MMPR(sumy2)<<MMPR(zij2)<<MMPR(sumyp2))
} // j==1
if (j==0)
{
D wij0=1.0/(D(j+1)*D(j+3)); D zij0=1.0/((D(j+1)*(D(j+2))));
D sumy0=s1*zt.zjy_sum(p0f,p0i); D sumyp0=zt.zj_sum(s0f*p0f,s0i*p0i); 
sum+=wij0*sumy0+zij0*sumyp0;

D wij2=1.0/(D(j+1)*D(j+4)); D zij2=1.0/((D(j+1)*(D(j+3))));
D sumy2=s1*zt.zjy2_sum(p0yf,p0yi); D sumyp2=zt.zjy_sum(s0f*p0yf,s0i*p0yi); 
sum+=wij2*sumy2+zij2*sumyp2;

D wij22=p0y2/(D(j+1)*D(j+5)); D zij22=p0y2/((D(j+1)*(D(j+4))));
D sumy22=s1*zt.zjy3_sum(); D sumyp22=zt.zjy2_sum(s0f,s0i); 
sum+=wij22*sumy22+zij22*sumyp22;
//MM_MSG("j==0 "<<MMPR(sum)<<MMPR(wij0)<<MMPR(sumy0)<<MMPR(zij0)<<MMPR(sumyp0) <<MMPR(wij2)<<MMPR(sumy2)<<MMPR(zij2)<<MMPR(sumyp2)<<MMPR(wij22)<<MMPR(sumy22)<<MMPR(zij22)<<MMPR(sumyp22))
} // j=0=
MM_MSG(MMPR(j)<<MMPR(sum))
zt.next_j();
} // i==0 j loop 
#endif

// pulling i==0 out was for anaaluzis but it does remove the extra mult lol 
const IdxTy _nterms=1000;
std::vector<D> terms;
for (IdxTy i=0; i<_nterms; ++i)
{

#if 0
for (IdxTy j=0; j<3; ++j)
{
if (j==2)
{
D wij=p2/(D(2*i+j+1)*D(2*i+j+3)); D zij=p2/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy=s1*zt.zjy_sum(); D sumyp=zt.zj_sum(s0f,s0i); 
sum+=fac*(wij*sumy+zij*sumyp);
}
if (j==1)
{
D wij1=1.0/(D(2*i+j+1)*D(2*i+j+3)); D zij1=1.0/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy1=s1*zt.zjy_sum(p1f,p1i); D sumyp1=zt.zj_sum(s0f*p1f,s0i*p1i); 
sum+=fac*(wij1*sumy1+zij1*sumyp1);

D wij2=p1y/(D(2*i+j+1)*D(2*i+j+4)); D zij2=p1y/((D(2*i+j+1)*(D(2*i+j+3))));

D sumy2=s1*zt.zjy2_sum(); D sumyp2=zt.zjy_sum(s0f,s0i); 
sum+=fac*(wij2*sumy2+zij2*sumyp2);
} // j==1
if (j==0)
{
D wij0=1.0/(D(2*i+j+1)*D(2*i+j+3)); D zij0=1.0/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy0=s1*zt.zjy_sum(p0f,p0i); 
D sumyp0=zt.zj_sum(s0f*p0f,s0i*p0i); 
sum+=fac*(wij0*sumy0+zij0*sumyp0);

D wij2=1.0/(D(2*i+j+1)*D(2*i+j+4)); D zij2=1.0/((D(2*i+j+1)*(D(2*i+j+3))));
D sumy2=s1*zt.zjy2_sum(p0yf,p0yi); 
D sumyp2=zt.zjy_sum(s0f*p0yf,s0i*p0yi); 
sum+=fac*(wij2*sumy2+zij2*sumyp2);

D wij22=p0y2/(D(2*i+j+1)*D(2*i+j+5)); D zij22=p0y2/((D(2*i+j+1)*(D(2*i+j+4))));
D sumy22=s1*zt.zjy3_sum();
D sumyp22=zt.zjy2_sum(s0f,s0i); 

sum+=fac*(wij22*sumy22+zij22*sumyp22);
} // j=0=

MM_MSG(MMPR(i)<<MMPR(j)<<MMPR(sum))

zt.next_j(); //  zffj*=yff;  ziij*=yii;  zfij*=yfi; zifj*=yif;
}
#endif
D term=0;
IdxTy j=0;
{
D wij0=1.0/(D(2*i+j+1)*D(2*i+j+3)); D zij0=1.0/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy0=s1*zt.zjy_sum(p0f,p0i); 
D sumyp0=zt.zj_sum(s0f*p0f,s0i*p0i); 
term+=fac*(wij0*sumy0+zij0*sumyp0);

D wij2=1.0/(D(2*i+j+1)*D(2*i+j+4)); D zij2=1.0/((D(2*i+j+1)*(D(2*i+j+3))));
D sumy2=s1*zt.zjy2_sum(p0yf,p0yi); 
D sumyp2=zt.zjy_sum(s0f*p0yf,s0i*p0yi); 
term+=fac*(wij2*sumy2+zij2*sumyp2);

D wij22=p0y2/(D(2*i+j+1)*D(2*i+j+5)); D zij22=p0y2/((D(2*i+j+1)*(D(2*i+j+4))));
D sumy22=s1*zt.zjy3_sum();
D sumyp22=zt.zjy2_sum(s0f,s0i); 

term+=fac*(wij22*sumy22+zij22*sumyp22);
} //

++j; zt.next_j(); // j==1
{
D wij1=1.0/(D(2*i+j+1)*D(2*i+j+3)); D zij1=1.0/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy1=s1*zt.zjy_sum(p1f,p1i); D sumyp1=zt.zj_sum(s0f*p1f,s0i*p1i); 
term+=fac*(wij1*sumy1+zij1*sumyp1);

D wij2=p1y/(D(2*i+j+1)*D(2*i+j+4)); D zij2=p1y/((D(2*i+j+1)*(D(2*i+j+3))));

D sumy2=s1*zt.zjy2_sum(); D sumyp2=zt.zjy_sum(s0f,s0i); 
term+=fac*(wij2*sumy2+zij2*sumyp2);
}

++j; zt.next_j(); // j==2
{
D wij=p2/(D(2*i+j+1)*D(2*i+j+3)); D zij=p2/((D(2*i+j+1)*(D(2*i+j+2))));
D sumy=s1*zt.zjy_sum(); D sumyp=zt.zj_sum(s0f,s0i); 
term+=fac*(wij*sumy+zij*sumyp);
}
if (isnan(term))
{

if (false) MM_MSG("my_taylor_nan "<< MMPR4(i,term,sum,beta)<<MMPR4(yii,zt.m_zii,yff,zt.m_zff))
if (!false) MM_MSG("my_taylor_explode "<< MMPR3(i,sum,beta)<<MMPR4(yii,zt.m_zii,yff,zt.m_zff))
failed=true;
break;
}
terms.push_back(term);
sum+=term;
//MM_MSG(MMPR4(i,term,sum,beta))
if ((term*term)<(1e-18*sum*sum )) break;
fac=fac*(-beta)/D(i+1);
zt.next_i(); zt.zj(); 

}
sum=0;
while (terms.size()!=0) { sum+=terms.back(); terms.pop_back(); } 
return sum;
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


CounterMap m_cm;

}; //mjm_timeline 



/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////

#ifdef  TEST_NUMERIC_TRICKS__

int main(int argc,char **args)
{
typedef mjm_numeric_tricks  Myt;

Myt x(argc,args);

if (!x.done()) x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

