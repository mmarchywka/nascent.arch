#ifndef MJM_BLOCK_LUT_H__
#define MJM_BLOCK_LUT_H__

#if 0 
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

#endif

#include "mjm_pickup_missing.h" 

// Tue Oct  8 18:04:17 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_block_lut   
// g++ -std=gnu++11 -DTEST_MJM_BLOCK_LUT -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_block_lut.h  -lpthread -lreadline

template <class Tr,class Tin=char, class Tout=long int >
class mjm_block_lut 
{
 typedef mjm_block_lut Myt;
protected:
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

//typedef typename std::map<Tin,Tout> LutTy;
// typedef typename Tr::MyBlock  MyBlock;
enum { CARD=1<<(8*sizeof(Tin)) } ; 
//typedef Tost[CARD] LutTy;
typedef Tout LutTy;

public:
mjm_block_lut() {zero();}
~mjm_block_lut() {}
const Tout & lut(const Tin & in) { return m_lut[in]; } 
Tout * lut() { return &m_lut[0]; } 
void zero() { FillToEnd(0,0); } 
void fill(const Tin & s, const Tin & e, const Tout & p) { 
if ( e==(CARD-1) ) 
 FillToEnd(s,p); 
else Fill(s,e,p); 

} 

#if 0 
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
#endif


StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
void Fill(const Tin & s, const Tin & e, const Tout & p)
{
for(Tin i=s; i<=e; ++i){
//MM_ERR(MMPR3(int(i),int(s),int(e)))
 m_lut[i]=p; 
}
}
void FillToEnd(const Tin & s, const Tout & p)
{
for(Tin i=s; i<(CARD-1); ++i){
//MM_ERR(MMPR3(int(i),int(s),int(CARD-1)))
 m_lut[i]=p; 
}
m_lut[CARD-1]=p;
}



bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }

LutTy m_lut[CARD];


}; // mjm_block_lut

//////////////////////////////////////////////

template <class Tr,class Tin=char, class Tout=long int >
class mjm_block_lut_map : public std::map<typename Tr::StrTy, mjm_block_lut< Tr > >  
{
 typedef mjm_block_lut_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_block_lut< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_block_lut_map() {}
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

}; // mjm_block_lut_map




////////////////////////////////////////////
#ifdef  TEST_MJM_BLOCK_LUT
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
ss<<" MJM_BLOCK_LUT "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_block_lut<Tr>  Myt;
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

#endif // MJM_BLOCK_LUT_H__ 