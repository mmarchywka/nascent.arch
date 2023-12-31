#ifndef MJM_CHOICES_ABBRV_H__
#define MJM_CHOICES_ABBRV_H__

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


// Thu Feb 24 19:36:41 EST 2022
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_choices_abbrv   
// g++  -Wall -std=gnu++11 -DTEST_MJM_CHOICES_ABBRV -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_choices_abbrv.h  -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_choices_abbrv("mjm_choices_abbrv"
, "  ");

template <class Tr>
class mjm_choices_abbrv 
{
 typedef mjm_choices_abbrv Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef std::map<StrTy, StrTy> AbbrMap;
public:
mjm_choices_abbrv() {}
~mjm_choices_abbrv() {}
typedef AbbrMap abbreviate_map;
template <class Tv> IdxTy abbreviate(AbbrMap & m, const Tv  & inv, const IdxTy flags )
{ return Abbreviate(m,inv,5,flags); }

template <class Tv> IdxTy abbreviate(AbbrMap & m, const Tv  & inv, const IdxTy lenmin,const IdxTy flags )
{ return Abbreviate(m,inv,lenmin,flags); }


StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
template <class Tv> IdxTy Abbreviate(AbbrMap & m, const Tv  & inv, const IdxTy lenmin,const IdxTy flags )
{
std::map<StrTy, IdxTy> f;
MM_LOOP(ii,inv) { ++f[(*ii)]; }
std::map<IdxTy, std::vector<StrTy> > hm  ;
MM_LOOP(ii,f) {if (m.find((*ii).first)==m.end())  m[(*ii).first]="";
 hm[(*ii).second].push_back((*ii).first); }
AltMap(m,f,lenmin,flags);
if (false) {
MM_LOOP(ii,f)
{
//MM_LOOP(jj,f) if ((*ii).first!=(*jj).first) MinExt(m,(*ii).first,(*jj).first,5,flags); 
MM_LOOP(jj,f) if ((*ii).first!=(*jj).first) MinExt(m,(*ii).first,(*jj).first,lenmin,flags); 
} // ii 
}
return 0;
} // Abbreviate

template <class Tv> IdxTy AltMap(AbbrMap & m, const Tv  & f, const IdxTy lenmin,const IdxTy flags )
{
const bool local_min=true;
IdxTy minext=lenmin;
IdxTy maxlen=0;
MM_LOOP(ii,f)
{
if (local_min )  minext=lenmin;
const auto & is=(*ii).first;
MM_LOOP(jj,f)
{
const auto & js=(*jj).first;
const IdxTy jlen=js.length();
if (jlen>maxlen) maxlen=jlen;
if (is!=js)
{
const char * p1=is.c_str();
const char * p2=js.c_str();
IdxTy i=0;
IdxTy lsame=0;
while (p1[i])
{
if (p1[i]!=p2[i]) break;
++i;
} // i 
lsame=i;
if (lsame>minext) minext=lsame;
} // !=  
} // jj 
if (local_min)
{
++minext;
char c[maxlen+4];
const char * p1=is.c_str();
memcpy(c,p1,is.length()+1);
c[minext]=0;
m[is]=StrTy(c);

} // local_min
} //  ii 
if (!local_min) {
++minext;
char c[maxlen+4];
MM_LOOP(ii,f)
{
const auto & s=(*ii).first;
const char * p1=s.c_str();
memcpy(c,p1,s.length()+1);
c[minext]=0;
m[(*ii).first]=StrTy(c);

} // ii 
} // local_min
return 0; 
} // AltMap
IdxTy MinExt( AbbrMap & m, const StrTy & k1, const StrTy & k2,const IdxTy szmin, IdxTy flags)
{
const char * p1=k1.c_str();
const char * p2=k2.c_str();
char k1c[k1.length()+3+szmin];
char k2c[k2.length()+3+szmin];
memcpy(k1c,p1,strlen(p1)+1);
memcpy(k2c,p2,strlen(p2)+1);

if (m.find(k1)!=m.end()) if (m[k1].length())  p1=m[k1].c_str();
if (m.find(k2)!=m.end()) if (m[k2].length()) p2=m[k2].c_str();
IdxTy i=0;
while (p1[i])
{
//k1c[i]=p1[i]; k2c[i]=p2[i];
if (p1[i]!=p2[i]) break;
++i;
}
if (i>szmin) { k1c[i+1]=0; k2c[i+1]=0; }
else
{  k1c[szmin+1]=0; k2c[szmin+1]=0; } 

if (p1[i]!=p2[i])
{
m[k1]=StrTy(k1c);
m[k2]=StrTy(k2c);

}  
else
{m[k1]=k1; m[k2]=k2;}
//while (true) { StrTy v1=m[k1]; StrTy v2=m[k2]; } // true  

return 0;
}  // MinExt

// MEMBERS



}; // mjm_choices_abbrv

//////////////////////////////////////////////

template <class Tr>
class mjm_choices_abbrv_map : public std::map<typename Tr::StrTy, mjm_choices_abbrv< Tr > >  
{
 typedef mjm_choices_abbrv_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_choices_abbrv< Tr> >   Super;
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
mjm_choices_abbrv_map() {}
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

}; // mjm_choices_abbrv_map




////////////////////////////////////////////
#ifdef  TEST_MJM_CHOICES_ABBRV
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
typedef tester_< mjm_choices_abbrv <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_CHOICES_ABBRV "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_choices_abbrv<Tr>  Myt;
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

#endif // MJM_CHOICES_ABBRV_H__ 
