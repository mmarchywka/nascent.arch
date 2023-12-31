#ifndef MJM_2D_STATES_H__
#define MJM_2D_STATES_H__

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


// Wed Jun 12 21:08:58 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_2d_states   
// g++ -std=gnu++11 -DTEST_MJM_2D_STATES -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_2d_states.h  -lpthread -lreadline

template <class Tr,class Tstate>
class mjm_2d_states 
{
 typedef mjm_2d_states Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_2d_states():m_ptr(0),m_default(),m_x(0),m_y(0),m_sz(0) {}
mjm_2d_states(const Tstate & def):m_ptr(0),m_default(def),m_x(0),m_y(0),m_sz(0) {}
~mjm_2d_states() {delete [] m_ptr; }
enum {BAD=~0 };
// needs assign and copy doh
private:
Myt & operator=(const Myt & x) { return *this; } 
mjm_2d_states(const Myt & x):m_ptr(0),m_x(0),m_y(0) {}

public:

IdxTy offset(const IdxTy x, const IdxTy y) const  { return Offset(x,y,m_x,m_y); }
const Tstate &  state(const IdxTy x, const IdxTy y) const  { return State(x,y); }
void state(const Tstate & s,const IdxTy x, const IdxTy y) { State(s,x,y); }
void clear() { Clear(m_ptr,m_sz); }
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
IdxTy next_not(IdxTy & x, IdxTy & y, const IdxTy x0,const IdxTy y0,  const IdxTy dx, const IdxTy dy, const Tstate & d) const
{
int xi=x;
int yi=y;
int rc= NextNot(xi,yi,x0,y0,dx,dy,d); 
x=xi; y=yi;
return rc;
}

private:
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
// position after having NOT moved once in dx dy; 
// this is based on there being at least one char on screen and
// not forcing a scroll, 
IdxTy NextNot(int & x,int & y, const IdxTy x0,const IdxTy y0,  const int dx, const int dy, const Tstate & d) const
{
x=x0;
y=y0;
if (m_x==0) return BAD;
if (m_y==0) return BAD;
const IdxTy flags=0;
IdxTy rc=Next(x,y,x,y,dx,dy,flags);
while ( State(x,y)==d)
{
if (rc==0) { 
// just moved into bad area, should not occur going left 
// no need as with vi to be on a new space, 
rc=Next(x,y,x,y,0,1,flags); 

} 
else if (Bit(rc,1)) { rc=Next(x,y,x,y,0,1,flags); } 
// ignore for now, want to fixs though 
else break; // rc=Next(x,y,x,y,dx,dy,flags);

}

return 0;
}

IdxTy Next(int & x, int & y, const IdxTy x0,const IdxTy y0,  const int dx, const int dy, const IdxTy flags) const
{
x=x0+dx;
y=y0+dy;
IdxTy rc=0;
if (m_x<=x){rc|=1; x=m_x-1; } 
if (0>x){rc|=4; x=0; } 
if (m_y<=y){rc|=2; y=m_y-1;  } 
if (0>y){rc|=8;y=0;  } 
//return m_ptr[offset(x,y)]; 

return rc; 
}

template <class Ty> 
IdxTy NeighborHisto(Ty & m, const IdxTy  x, const IdxTy  y) const
{


return 0;
}



void State(const Tstate & s,const IdxTy x, const IdxTy y) 
{ 
	//State(s,x,y); 
	Expand(x,y);
m_ptr[Offset(x,y,m_x,m_y)] = s; 
}
const Tstate &  State(const IdxTy x, const IdxTy y)  const 
{ 
if (m_x<=x) return m_default;
if (m_y<=y) return m_default;
return m_ptr[offset(x,y)]; 
}
void Expand(const IdxTy x, const IdxTy y) 
{
const bool xx=(x>=m_x);
const bool yy=(y>=m_y);
if (xx||yy)
{
const IdxTy xn=xx?(x+1):m_x;
const IdxTy yn=yy?(y+1):m_y;
//IdxTy sz=(xx?x:m_x)*(yy?y:m_y);
IdxTy sz=xn*yn;
Tstate * pnew= new Tstate[sz];
Clear(pnew,sz);
// should just make a new object and assign 
// otherwise we have ANOTHER offset calc that can be inconsistent with other 
for(IdxTy i=0; i<m_x; ++i)
{
for(IdxTy j=0; j<m_y; ++j){  pnew[Offset(i,j,xn,yn)]=m_ptr[Offset(i,j,m_x,m_y)];    }
} // i 
delete [] m_ptr;
m_ptr=pnew;
m_x=xn; // xx?x:m_x;
m_y=yn; // yy?y:m_y;
m_sz=sz;

} // if 

} // Expand

IdxTy Offset(const IdxTy x, const IdxTy y, const IdxTy mx, const IdxTy my ) const 
{ return mx*y+x; }
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
void Clear(Tstate * _ptr, const IdxTy sz )
{
 Tstate * p=_ptr;
const Tstate * pend=p+sz;
while (p!=pend){  *p=m_default; ++p; } 
}


StrTy Dump(const IdxTy flags=0) 
{Ss ss;  
ss<<MMPR3(m_x,m_y,m_sz)<<CRLF;
for(IdxTy i=0; i<m_x; ++i)
for(IdxTy j=0; j<m_y; ++j)
ss<<MMPR3(i,j,m_ptr[offset(i,j)])<<CRLF;

return ss.str(); }

Tstate * m_ptr;
Tstate m_default;
IdxTy m_x,m_y,m_sz;

}; // mjm_2d_states

//////////////////////////////////////////////

template <class Tr,class Tstate>
class mjm_2d_states_map : public std::map<typename Tr::StrTy, mjm_2d_states< Tr,Tstate > >  
{
 typedef mjm_2d_states_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_2d_states< Tr,Tstate> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_2d_states_map() {}
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

}; // mjm_2d_states_map




////////////////////////////////////////////
#ifdef  TEST_MJM_2D_STATES
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
ss<<" MJM_2D_STATES "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_2d_states<Tr,int>  Myt;
//Myt x(argc,args);
Myt x(~0);

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
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
if (cmd=="clear") { x.clear(); MM_ERR(x.dump()) }
if (cmd=="state") {
IdxTy xi=cip.myatoi(cip.p1);
IdxTy y=cip.myatoi(cip.p2);
IdxTy s=cip.myatoi(cip.wif(3));

MM_ERR(MMPR3(xi,y,s))
x.state(s,xi,y);

 MM_ERR(x.dump()) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_2D_STATES_H__ 
