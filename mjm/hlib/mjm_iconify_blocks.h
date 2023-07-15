#ifndef MJM_ICONIFY_BLOCKS_H__
#define MJM_ICONIFY_BLOCKS_H__

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

// Tue Oct  8 17:56:32 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_iconify_blocks   
// g++ -std=gnu++11 -DTEST_MJM_ICONIFY_BLOCKS -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_iconify_blocks.h  -lpthread -lreadline

template <class Tr>
class mjm_iconify_blocks 
{
 typedef mjm_iconify_blocks Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef int InTy;

class range_desc
{
public:

range_desc() : m_start(bad()), m_end(bad()), m_visibility(bad()) {Init(); }
range_desc(const StrTy & nm) : m_start(bad()), m_end(bad()),m_name(nm), m_visibility(bad()) {Init(); }
range_desc(const StrTy & nm, const IdxTy s) : m_start(s), m_end(bad()),m_name(nm), m_visibility(bad()) {Init(); }
range_desc(const StrTy & nm, const IdxTy s, const IdxTy e) : m_start(s), m_end(e),m_name(nm), m_visibility(1) {Init(); }
const StrTy & name() const { return m_name; } 
void start(const InTy & s) { m_start=s; }
void start(const InTy & x,const InTy & y) { m_start_x=x;m_start_y=y; }
void start(const InTy & x,const InTy & y, const InTy & page) { m_start_x=x;m_start_y=y; m_start_page=page;  }
void end(const IdxTy & e) { m_end=e; }
void end(const InTy & x,const InTy & y) { m_end_x=x;m_end_y=y; }
void end(const InTy & x,const InTy & y, const InTy page) { m_end_x=x;m_end_y=y;m_end_page=page; }
void visibility(const IdxTy v) { m_visibility=v;}
IdxTy visibility()const  { return m_visibility;}
bool visible()const  { return m_visibility!=0;}
void flip() { if (!visible()) m_visibility=1; else m_visibility=0; } 
void show() { m_visibility=1;}
void hide() { m_visibility=0;}
bool valid() const { return (m_start!=bad())&&(m_end!=bad()); } 
// if start had been found but not end, everything after is included but not valid 
bool in(const IdxTy n) const { return ((m_start<=n)&&(m_end>=n)); } 
bool pix_in_y(const InTy y) const { return ((m_start_y<=y)&&(m_end_y>=y)); } 
bool pix_in_x(const InTy x) const { return ((m_start_x<=x)&&(m_end_x>=x)); } 
bool pix_in(const InTy x,const InTy y) const { return pix_in_x(x)&&pix_in_y(y); } 
bool pix_in_v(const InTy x,const InTy y) const { return pix_in_x(x)&&pix_in_y(y) || pix_in_y(y)&&visible(); } 
bool whole(const IdxTy page ) const { return (m_start_page<page)&&(m_end_page>page); } 
bool pix_in(const InTy x,const InTy y,const InTy page) const 
{ 
bool starts= (m_start_page<=page);
bool ends= (m_end_page>=page);
if (!starts&&!ends) return false;
if (starts&&ends) 
{
if (whole(page)) return true;
return pix_in_x(x)&&pix_in_y(y); 
}
if (starts) return (m_start_x<=x)&&(m_start_y<=y);  
if (ends) return (m_end_x>=x)&&(m_end_y>=y);  
MM_ERR(" logic error "<<dump())
}
bool pix_in_v(const InTy x,const InTy y,const InTy page) const 
{ 
bool starts= (m_start_page<=page);
bool ends= (m_end_page>=page);
if (!starts&&!ends) return false;
if (starts&&ends) 
{
if (whole(page)) return true;
return pix_in_x(x)&&pix_in_y(y); 
}
if (starts) return ((visible()||m_start_x<=x))&&(m_start_y<=y);  
if (ends) return ((visible()||m_end_x>=x))&&(m_end_y>=y);  
MM_ERR(" logic error "<<dump())
}




IdxTy bad() { return ~0; }
StrTy dump(const IdxTy flags=0 ) const { return Dump(flags); } 
private:
StrTy Dump(const IdxTy flags ) const { 
Ss ss;
ss<<MMPR4(m_name,m_visibility,m_start,m_end);
ss<<MMPR4(m_start_x,m_start_y,m_end_x,m_end_y);
ss<<MMPR2(m_start_page,m_end_page);
return ss.str(); 
} 
void Init() 
{
m_start_page=bad();
m_end_page=bad();
}
IdxTy m_start, m_end;
// the only one that matters really is y, the x is too hard for now 
InTy m_start_x,m_start_y,m_start_page, m_end_x,m_end_y,m_end_page;
StrTy m_name;
IdxTy m_visibility;

}; // range_desc

typedef range_desc Ra;
// eventually this needs to tokenize strings, make all the maps int keys
typedef std::map<StrTy, Ra> Map;
typedef std::map<StrTy,IdxTy> Blist;
typedef std::map<IdxTy, Blist> MapLoc;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_iconify_blocks(): m_visibility(1), m_bl(0) {Init();}
~mjm_iconify_blocks() {}
IdxTy bad() { return ~0; }

bool next(const IdxTy pc, const IdxTy op ) { return Next(pc,op); } 
#if 0 
if( m_bl==0) return visible(); 
if (m_bl->in(pc)) 
{
// this may change asynchronously although this should only be 
// done at block boundaries 
//m_visibility=m_bl->visibility();
}
else
{ // find a new block containing this - wait for special to execute  
m_bl=0;
m_visibility=1;

} // finding new block 
return visible(); 
} 
#endif

bool visible() { return m_visibility!=0; } 

void start_block(const StrTy & nm, const IdxTy p, const IdxTy flags=0) { StartBlock(nm,p,flags); } // start_block
void start_block(const StrTy & nm, const IdxTy p, const InTy x, const IdxTy y,const IdxTy page, const InTy flags) 
{ StartBlock(nm,p,x,y,page,flags); } // start_block

bool has(const StrTy & nm) const { return m_map.find(nm)!=m_map.end(); } 
StrTy name() const {if (m_bl==0) return StrTy();  return m_bl->name(); } 
void end_block(const StrTy & nm, const IdxTy p) { EndBlock(nm,p); }
void end_block(const StrTy & nm, const IdxTy p,const InTy x, const InTy y,const InTy page) { EndBlock(nm,p,x,y,page); }
IdxTy flip_at(const IdxTy x, const IdxTy y) { return FlipAt(x,y); } 
//IdxTy flip_at(const IdxTy x, const IdxTy y) { return FlipAtDebug(x,y); } 
IdxTy size() const { return m_map.size(); } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
Ra operator[](const StrTy &nm) const { return Block(nm); } 
Ra operator[](const IdxTy n ) const { return Block(m_order_found[n]); } 
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
void Init()
{
 m_mutex_vector=MutexVector(1);
}

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
Ra Block(const StrTy &nm) const { 
auto ii=m_map.find(nm);
if (ii==m_map.end()) return Ra();
return (*ii).second; 
} 
IdxTy  FlipAt(const InTy x, const InTy y) 
{ 
IdxTy n=0;
EnterSerial(MAP_MU);
MM_LOOP(ii,m_map)
{
Ra & ra=(*ii).second;
if (ra.pix_in_v(x,y)){
 ra.flip(); ++n;
//MM_ERR(" flipped "<<ra.dump()) 
}
//if (!ra.visible()) MM_ERR(" hidden "<<ra.dump()) 
//MM_ERR(" block  "<<MMPR2(x,y)<<ra.dump()) 

//else MM_ERR("not  flipped "<<MMPR2(x,y)<<ra.dump()) 
}
ExitSerial(MAP_MU);
return n; 
} 
IdxTy  FlipAtDebug(const InTy x, const InTy y) 
{ 
IdxTy n=0;
EnterSerial(MAP_MU);
MM_LOOP(ii,m_order_found)
{
Ra & ra=m_map[*ii]; // (*ii).second;
if (ra.pix_in_v(x,y)){
 ra.flip(); ++n;
MM_ERR(" flipped "<<ra.dump()) 
}
//if (!ra.visible()) MM_ERR(" hidden "<<ra.dump()) 
MM_ERR(" block  "<<MMPR2(x,y)<<ra.dump()) 

//else MM_ERR("not  flipped "<<MMPR2(x,y)<<ra.dump()) 
}
ExitSerial(MAP_MU);
return n; 
} 





bool Next(const IdxTy pc, const IdxTy op ) { 

if( m_bl==0) return visible(); 
if (m_bl->in(pc)) 
{
// this may change asynchronously although this should only be 
// done at block boundaries 
//m_visibility=m_bl->visibility();
}
else
{ // find a new block containing this - wait for special to execute  
m_bl=0;
m_visibility=1;

} // finding new block 
return visible(); 
} 
void StartBlock(const StrTy & nm, const IdxTy p, const InTy x, const InTy y, const IdxTy page, const IdxTy flags)
{
const IdxTy f2=flags|2;
StartBlock(nm,p,f2);
m_map[nm].start(x,y,page);
ExitSerial(MAP_MU);
}
void StartBlock(const StrTy & nm, const IdxTy p, const IdxTy flags=0)
{
EnterSerial(MAP_MU);
auto ii=m_map.find(nm);
if ( ii!=m_map.end())
{
Ra & b=(*ii).second;
m_visibility=b.visibility();
}
else
{
m_map[nm]= Ra(nm,p);
m_order_found.push_back(nm);
if (Bit(flags,0)) m_map[nm].hide();
else m_map[nm].show();
}
m_bl=&m_map[nm];
++m_s[p][nm];

if (!Bit(flags,1)) ExitSerial(MAP_MU);
} // StartBlock

void EndBlock(const StrTy & nm, const IdxTy p,const InTy x, const InTy y,const InTy page)
{
EndBlock(nm,p);
EnterSerial(MAP_MU);
m_map[nm].end(x,y,page);
ExitSerial(MAP_MU);

}

void EndBlock(const StrTy & nm, const IdxTy p)
{
EnterSerial(MAP_MU);
auto ii=m_map.find(nm);
if ( ii!=m_map.end())
{
Ra & b=(*ii).second;
if ( p!=0) b.end(p);
}
else
{
m_map[nm]= Ra(nm); 
if (p!=0) m_map[nm].end(p);
m_map[nm].show();

}
MM_ONCE(" change  1 here ", )
// TODO FIXME this was in but then inter-block spaced zapped no idea how other code ran 
// taking this out still leaves thing missing and each icon on new line 
// m_bl=0;
// this still swallows up the next text and puts in crlf. 
if (p!=0) {  m_bl=0; m_visibility=1; } 
++m_e[p][nm];
ExitSerial(MAP_MU);

} // EndBlock



StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }


IdxTy m_visibility;
Ra * m_bl;
Map m_map;
MapLoc m_s,m_e;
std::vector<StrTy> m_order_found;
}; // mjm_iconify_blocks

//////////////////////////////////////////////

template <class Tr>
class mjm_iconify_blocks_map : public std::map<typename Tr::StrTy, mjm_iconify_blocks< Tr > >  
{
 typedef mjm_iconify_blocks_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_iconify_blocks< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_iconify_blocks_map() { }

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

}; // mjm_iconify_blocks_map




////////////////////////////////////////////
#ifdef  TEST_MJM_ICONIFY_BLOCKS
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
ss<<" MJM_ICONIFY_BLOCKS "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_iconify_blocks<Tr>  Myt;
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

#endif // MJM_ICONIFY_BLOCKS_H__ 