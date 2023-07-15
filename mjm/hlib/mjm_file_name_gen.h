#ifndef MJM_FILE_NAME_GEN_H__
#define MJM_FILE_NAME_GEN_H__

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


// Mon May 27 05:51:38 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_file_name_gen   
// g++ -std=gnu++11 -DTEST_MJM_FILE_NAME_GEN -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_file_name_gen.h  -lpthread -lreadline

template <class Tr>
class mjm_file_name_gen 
{
 typedef mjm_file_name_gen Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_file_name_gen(): m_serial(0),m_mode(0)  {}
~mjm_file_name_gen() {}
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


const StrTy next() { return Next(); }  

const StrTy sfx(const StrTy & s) { return Sfx(s); }  
void   set_file( const StrTy & fn) {  SetFile(fn); }
const StrTy   set_file( const StrTy & dir, const StrTy & base, const StrTy & ext)
{  SetFile(dir,base,ext); return Next(); } 
StrTy friendly_file( const StrTy & fn)
{ return FriendlyFile(fn); }

StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
void SetFile( const StrTy & fn)
{ 
StrTy dir,base,ext;
const IdxTy sz=fn.length();
char c[sz+1];
memcpy(c,fn.c_str(),sz+1);
IdxTy p=sz;
while (p!= ~0)
{
if (c[p]=='.') {c[p]=0; ext=StrTy(c+p+1); }
if (c[p]=='/') {c[p]=0; base=StrTy(c+p+1); dir=StrTy(c); break; }
--p; 
}
if (p==~0)  base=StrTy(c);
m_dir=dir; m_base=base; m_ext=ext; 

}

void SetFile( const StrTy & dir, const StrTy & base, const StrTy & ext)
{ m_dir=dir; m_base=base; m_ext=ext; }

StrTy & Next()
{
Ss ss; ss<<m_dir; ss<<m_base;
++m_serial;
ss<<m_serial; // m_serial;
ss<<"."; ss<<m_ext; m_fn=ss.str(); return m_fn;
}
const StrTy Sfx(const StrTy & s) { 

Ss ss; ss<<m_dir; ss<<m_base;
ss<<"_";
ss<<s;
ss<<"."; ss<<m_ext; m_fn=ss.str(); return m_fn;
return m_fn; 

}  

StrTy FriendlyFile( const StrTy & fn)
{
//Ss ss;
const IdxTy sz=fn.length();
const char * p=fn.c_str();
char c[sz+1];
for(IdxTy i=0; i<sz; ++i)
{
if ( p[i]!=' ') c[i]=p[i];
else c[i]='_';
}
c[sz]=0;
return StrTy(&c[0]);
//return ss.str();
}



StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }

StrTy m_dir, m_base,m_ext,m_fn;
IdxTy m_serial;
IdxTy m_mode;


}; // mjm_file_name_gen

//////////////////////////////////////////////

template <class Tr>
class mjm_file_name_gen_map : public std::map<typename Tr::StrTy, mjm_file_name_gen< Tr > >  
{
 typedef mjm_file_name_gen_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_file_name_gen< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_file_name_gen_map() {}
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

}; // mjm_file_name_gen_map




////////////////////////////////////////////
#ifdef  TEST_MJM_FILE_NAME_GEN
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


int main(int argc,char **args)
{
typedef mjm_file_name_gen<Tr>  Myt;
//Myt x(argc,args);
Myt x;

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
if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
 if (cmd=="set") { 
if (sz>3) { x.set_file(li.word(1),li.word(2),li.word(3)); } 

 }
 if (cmd=="get") {  MM_MSG(MMPR(x.sfx(li.word(1)))) } 
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_FILE_NAME_GEN_H__ 
