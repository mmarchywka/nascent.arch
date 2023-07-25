#ifndef MJM_SVG_PLOTTER_BASE_H__
#define MJM_SVG_PLOTTER_BASE_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_floorplan_points.h"

#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"

#include "mjm_canned_methods.h"


#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"

#include "mjm_canned_methods.h"
// add day number to notes 
#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_svg_writer.h"

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


// Sat Apr 30 19:10:56 EDT 2022
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_svg_plotter_base   
// g++  -Wall -std=gnu++11 -DTEST_MJM_SVG_PLOTTER_BASE -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_svg_plotter_base.h  -o mjm_svg_plotter_base.out -lpthread -lreadline

mjm_global_credits::credit __credit__mjm_svg_plotter_base("mjm_svg_plotter_base"
, "  ");

template <class Tr>
class mjm_svg_plotter_base 
{
 typedef mjm_svg_plotter_base Myt;
protected:
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
typedef mjm_canned_methods Canned;
typedef mjm_svg_writer Sw;

typedef mjm_floorplan_points<Tr> FloorPlan;



enum PASSES {DEFS,BODY,POST};
public:
// floating point equality is not assured 
//typedef D Zidx;
typedef IdxTy Zidx;
typedef std::map<D,IdxTy> Zmap;
// a component can have multiple Z-levels, traversal is
// potential mess ... 
class _Zord
{
public:
_Zord() : m_z(0),znext(m_z) {}
_Zord(const IdxTy z ) : m_z(z),znext(m_z) {}
void next(const Zidx n) { if ((znext<n)||(znext==m_z)) znext=n; } 
bool done() { return (m_z==znext); } 
void  next() {  m_z=znext; } 
const Zidx&  z() const { return m_z; }
void dump(Ss  & ss) { ss<<MMPR2(m_z,znext); } 
StrTy  dump() { Ss  ss; dump(ss); return ss.str(); } 
Zidx m_z;
Zidx znext;
}; // _Zord
typedef _Zord Zord;
typedef std::vector<Myt * > Kids;
typedef std::map<Zidx,Kids > KidMap;

class _Point
{

// FIXME doh put this somwhere lol 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); }
int myatoi(const char * c) const { return Canned::myatoi(c); }


public:
_Point() {}
_Point(const StrTy & n) {m_v.push_back(n);}
_Point(const D & n) {m_v.push_back(f(n));}
_Point(const StrTy & n, const StrTy & m ) {m_v.push_back(n); m_v.push_back(m);}
_Point(const D & n, const D & m ) {m_v.push_back(f(n)); m_v.push_back(f(m));}

const IdxTy size() const { return m_v.size(); } 
const StrTy & str(const IdxTy i) const { return m_v[i]; }
const D operator[](const IdxTy i) const { return atof( m_v[i].c_str()); }
const D d(const IdxTy i) const { return atof( m_v[i].c_str()); }
// could appcet hex etc?
const IdxTy  i(const IdxTy i) const { return myatoi( m_v[i]); }
StrTy f(const D & x) const 
{ Ss ss; ss<<std::setprecision(16)<<x; return ss.str(); }
// TODO rationals? Think about this ... 

// slow but verstile.. 
typedef std::vector<StrTy> V;

V m_v; 

}; // _Point


typedef _Point Point;


class _Xform
{

public:
_Xform() {}
_Xform(const D& x, const D & y, const D & max, const D & may):
 m_x(x), m_ax(max), m_y(y), m_ay(may){}
Point xform(const Point x, const IdxTy flags=0)
{

return Point((x.d(0)-m_x)*m_ax, (x.d(1)-m_y)* m_ay);

}
D m_x, m_ax, m_y, m_ay;
}; // _Xform

typedef _Xform Xform;
typedef std::vector<Point> Points;
typedef std::map<StrTy,StrTy> Values;
typedef std::map<StrTy,StrTy> Properties;
public:


class _Params
{
public:

}; // _Params
typedef _Params params_type;




// I guess these could be virtual, don't emember ctor order on hierarchy.. 
mjm_svg_plotter_base() {Init();}
virtual ~mjm_svg_plotter_base() {Dtor();}
IdxTy size() const { return m_kids.size(); } 
virtual void make_layout( FloorPlan & fp ) {}
FloorPlan & fp() { return m_fp;}
virtual IdxTy ctor(const Ragged & r, const StrTy & type, const StrTy & nm, params_type & param, const IdxTy flags) 
{
MM_ERR(" base class emptry ctor ")
return 0;
} // ctor 


IdxTy write(OsTy & os, const IdxTy flags) //const
{
Zord zord;
Sw sw;

write_svg(os,sw,this,zord,DEFS,flags);
write_svg(os,sw,this,zord,BODY,flags);
write_svg(os,sw,this,zord,POST,flags);



return 0; 
}

virtual IdxTy add( Myt * child, const IdxTy flags) // const 
{
m_kids[0].push_back(child);
return 0; 
}
virtual IdxTy write_svg(OsTy & os, Sw & sw, Myt * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{
//Ss ss; ss<<MMPR3(__FILE__,m_name,m_type);
// viewer  omplains about crap before xml  
if(false)
{Ss ss; ss<<"auck "<<MMPR3(m_rtti,m_name,m_type)<<MMPR2(zord.dump(),pass);
 os<<sw.comment_text(ss);
os<<CRLF;
}
write_self_svg(os,sw,this,zord,pass,flags);
write_kids_svg(os,sw,this,zord,pass,flags);
write_self_end_svg(os,sw,this,zord,pass,flags);

return 0;
} // write_svg

// at least override this... 
virtual IdxTy write_self_svg(OsTy & os, Sw & sw, Myt * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{

Ss ss; 
ss<<" should not get here ... "<< MM_MARKF<<" ";
ss<<MMPR3(m_rtti,m_name,m_type)<<MMPR2(zord.dump(),pass);
 os<<sw.comment_text(ss);
os<<CRLF;
return 0; 

} // write_self_svg
virtual IdxTy write_self_end_svg(OsTy & os, Sw & sw, Myt * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{

Ss ss; 
ss<<""<< MM_MARKF<<" ";
ss<<MMPR3(m_rtti,m_name,m_type)<<MMPR2(zord.dump(),pass);
 os<<sw.comment_text(ss);
os<<CRLF;
return 0; 

} // write_self_svg



virtual IdxTy write_kids_svg(OsTy & os, Sw & sw, Myt * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{
//Ss ss; ss<<MMPR3(__FILE__,m_name,m_type);

// draw everyting in the map  at this level 
// in theory preabmble and defs dont need to be in z-ordder... 
auto ii=m_kids.begin();
while (ii!=m_kids.end())
{
if ((*ii).first==zord.z())
{
MM_LOOP(jj,(*ii).second ){  (*jj)->write_svg(os,sw,root,zord,pass,flags); }
++ii;
if  (ii!=m_kids.end()) zord.next((*ii).first);
break;
} // ==

++ii;
} // ii 

return 0;
} // write_kids_svg 


template <class Ty,class Tv > void setup(Ty & d, const Tv & def, const StrTy & nm, ReadWriteMap&  rwm,const IdxTy flags=0)
{

//m_title="Set a title" ; rwm.get("title",m_title);
//m_title="Set a title" ; rwm.get("title",m_title);
d=def; rwm.get(nm,d);
return 0; 
} //setup 



StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
protected:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
void ClearMap()
{
MM_LOOP(ii,m_kids) { MM_LOOP(jj,(*ii).second) { delete (*jj); }}

} 
void Init()
{

} // Init
void Dtor()
{
// The kids are not owned here let them be culled from source
//ClearMap();

} // Dtor 

// MEMBERS
//Zmap m_zmap;
StrTy m_name,m_type,m_rtti;
KidMap m_kids;
Points m_points;
Values m_values;
Properties m_prop;
FloorPlan m_fp;
}; // mjm_svg_plotter_base

//////////////////////////////////////////////

template <class Tr>
class mjm_svg_plotter_base_map : public std::map<typename Tr::StrTy, mjm_svg_plotter_base< Tr > >  
{
 typedef mjm_svg_plotter_base_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_svg_plotter_base< Tr> >   Super;
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
mjm_svg_plotter_base_map() {}
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

}; // mjm_svg_plotter_base_map




////////////////////////////////////////////
#ifdef  TEST_MJM_SVG_PLOTTER_BASE
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
typedef tester_< mjm_svg_plotter_base <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_SVG_PLOTTER_BASE "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_svg_plotter_base<Tr>  Myt;
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

#endif // MJM_SVG_PLOTTER_BASE_H__ 