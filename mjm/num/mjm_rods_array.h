#ifndef MJM_RODS_ARRAY_H__
#define MJM_RODS_ARRAY_H__

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


// Sun Jan  5 09:38:39 EST 2020
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_rods_array   
// g++ -std=gnu++11 -DTEST_MJM_RODS_ARRAY -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_rods_array.h  -lpthread -lreadline

// read-only data structure array 
// Td must supply and index operation and a target data type
template <class Tr , class Td >
class mjm_rods_array 
{
 typedef mjm_rods_array Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

public:
// this should be a COMPILED type without malloced pts,
// put everything inline eventually no fixed size objects. 
typedef typename  Td::data_type data_type;
enum { EL_SIZE=sizeof(data_type)};
// typedef typename Tr::MyBlock  MyBlock;
mjm_rods_array() {Init(); }
~mjm_rods_array() {Release(); }
//void size(const IdxTy sz,const IdxTy flags) { Resize(Td::size(sz,flags)); }  
void size(const IdxTy sz,const IdxTy flags) { Resize(sz,flags); }  
IdxTy size() const  { return m_sz; }  
IdxTy used() const  { return m_used; }  
const data_type & operator()(IdxTy i) const { 

//MM_ERR( " access "<<MMPR3(i,m_sz,m_used))
return m_ptr[i];}
data_type & operator()(IdxTy i)  {
//MM_ERR( " access "<<MMPR3(i,m_sz,m_used))
 return m_ptr[i];}

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
void Resize(const IdxTy sz,const IdxTy flags)
{
data_type * p = new data_type[sz];
const bool clean=Bit(flags,0);
const bool discard=Bit(flags,1);
const bool mark_as_used=Bit(flags,2);
const IdxTy used=mark_as_used?sz:(discard?0:m_used);
// this would not work because the fudding object dtors get fudded up.
// need to think this through again and use void * through out
// but then the object copy and ass-ign wont' work fudd 
//Copy(p,m_ptr,used);
Release();
m_used=used;
m_ptr=p;
m_sz=sz;
//if (clean) Clear(2); 
}
void Copy(data_type * p,const data_type * ptr, IdxTy used)
{
if (ptr==0) return; 
::memcpy(p,ptr,used*EL_SIZE);
}
// the fudding alloc allocs shot for the invalid fudders doh 
void xxClear(const IdxTy flags)
{
const bool def=Bit(flags,0);
const bool unused=Bit(flags,1);
const IdxTy start=unused?m_used:0;
if (def)
{
data_type d;
for (IdxTy i=start; i<m_sz; ++i)  m_ptr[i]=d;
//for (IdxTy i=start; i<m_sz; ++i)  (data_type*)(m_ptr+i*EL_SIZE)=d;

}
else { ::memset(m_ptr+start,0,(m_sz-start)*EL_SIZE); }
//else { ::memset(m_ptr+start*EL_SIZE,0,(m_sz-start)*EL_SIZE); }
}

void Release()
{
if (m_ptr==0) return; 
delete [] m_ptr; m_ptr=0; 
m_sz=0; m_used=0;
}
void Init()
{
m_ptr=0;
m_sz=0;
m_used=0;
}
// this needs a rods mem manager
// unforuntaly dtor on empty crap here 
data_type * m_ptr;
//void * m_ptr;
IdxTy m_sz;
IdxTy m_used;

}; // mjm_rods_array

//////////////////////////////////////////////

template <class Tr,class Td >
class mjm_rods_array_map : public std::map<typename Tr::StrTy, mjm_rods_array< Tr, Td > >  
{
 typedef mjm_rods_array_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_rods_array< Tr, Td > >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_rods_array_map() {}
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

}; // mjm_rods_array_map




////////////////////////////////////////////
#ifdef  TEST_MJM_RODS_ARRAY
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
ss<<" MJM_RODS_ARRAY "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_rods_array<Tr>  Myt;
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

#endif // MJM_RODS_ARRAY_H__ 
