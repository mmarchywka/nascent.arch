#ifndef MJM_CHAR_STAR_MAP_H__
#define MJM_CHAR_STAR_MAP_H__

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


// Sat Nov 30 08:33:08 EST 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_char_star_map   
// g++ -std=gnu++11 -DTEST_MJM_CHAR_STAR_MAP -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_char_star_map.h  -lpthread -lreadline

class mjm_string_less
{
public:
template <class Ty> 
bool operator()(const Ty & x, const Ty & y)  { int z=strcmp(x,y); return (z<0); } 

}; // mjm_string_less

template <class Tr,class K, class V >
class mjm_char_star_map : public std::map<K,V,mjm_string_less> 
{
 typedef mjm_char_star_map Myt;
typedef std::map<K,V,mjm_string_less>  Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_char_star_map() {}
~mjm_char_star_map() {}
//template <int n> V & operator[](const char[n]  k) { return (*this)[K(&k)]; }  
 V & operator[](const char *   k) { return (Super(*this))[K(k)]; }  

StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
}; // mjm_char_star_map


template <class Tr> 
class mjm_string_t_less
{
typedef mjm_read_buffer<Tr>  RdBuf;
public:
mjm_string_t_less(RdBuf ** p) : p_buf(p) {}
//template <class Ty> 
bool operator()(const int  x, const int  y)  
{ int z=strcmp((**p_buf).string(x),(**p_buf).string(y)); 
return (z<0); } 
RdBuf ** p_buf;
}; // mjm_string_less


template <class Tr, class Ch=char> // ,class K, class V >
class mjm_char_star_int_map //  : public std::map<K,V,mjm_string_less> 
{
 typedef mjm_char_star_int_map Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
typedef mjm_read_buffer<Tr> RdBuf;
typedef mjm_string_t_less<Tr>  Tl;
typedef typename std::map<IdxTy , IdxTy, Tl  > Map;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_char_star_int_map():m_tl(Tl(&m_buf)), m_map(m_tl), m_buf(new RdBuf(10)) {}
~mjm_char_star_int_map() { delete m_buf; }
//template <int n> V & operator[](const char[n]  k) { return (*this)[K(&k)]; }  
// V & operator[](const char *   k) { return (*this)[K(k)]; }  
// to find a key , append it to the buffer and do a normal look up.
// then delete it from buffer. 
 // or make a dummy type

StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }

Tl m_tl;
Map m_map;
RdBuf*  m_buf;

}; // mjm_char_star_int_map



//////////////////////////////////////////////





//////////////////////////////////////////////

template <class Tr , class K, class V>
class mjm_char_star_map_map : public std::map<typename Tr::StrTy, mjm_char_star_map< Tr,K,V > >  
{
 typedef mjm_char_star_map_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_char_star_map< Tr,K,V >  >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_char_star_map_map() {}
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

}; // mjm_char_star_map_map




////////////////////////////////////////////
#ifdef  TEST_MJM_CHAR_STAR_MAP
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
ss<<" MJM_CHAR_STAR_MAP "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_char_star_map<Tr>  Myt;
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

#endif // MJM_CHAR_STAR_MAP_H__ 
