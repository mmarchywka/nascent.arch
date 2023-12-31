#ifndef MJM_BIBTEX_FORMAT_H__
#define MJM_BIBTEX_FORMAT_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_bibtex_entry.h"
#include "mjm_bibtex_render.h"
#include "mjm_collections.h"

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


// Mon Jun 28 05:43:15 EDT 2021
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_bibtex_format   
// g++  -Wall -std=gnu++11 -DTEST_MJM_BIBTEX_FORMAT -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_bibtex_format.h  -lpthread -lreadline


template <class Tr>
class mjm_bibtex_format_entry 
{
 typedef mjm_bibtex_format_entry Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
 typedef mjm_bibtex_entry<Tr> BibEntry;
 typedef mjm_bibtex_entry_map<Tr> BibMap;
typedef mjm_bibtex_render<Tr> Renderer;
typedef mjm_ragged_table Ragged;
class form_param 
{

public:

private:

}; // form_param
typedef form_param FormParam;
typedef FormParam form_param_type;

public:
mjm_bibtex_format_entry() {}
~mjm_bibtex_format_entry() {}

virtual IdxTy  format(Renderer & r, const BibEntry & e, const FormParam  fp )
{
// finally read part of, 
// http://tug.org/texmf-docs/bibtex/ btxdoc.tex
// It does say include everything in the bibtex which the ld+json can now do.
// But, I think I can generalize
// It does have a croffref concept that I ignored but
// now may generalize with mjm_in to allow entry details specific
// to the work ( one thing produced by identifiable authors )
// enclosed in some larger thing indefinintley 
// the concept of DOI does not exist yet AFAICT
// TODO FIXME 
// FIXME My code does not support "@" strings or related "#"
 
// look at entry type, real_type, and then pick tokens
// to pass to renderer
// composition, person, material
// compositions are normally hierarchial and have specific
// location modifiers ( like page or para or chapter ) from the bibtex entry 
// composition : most interest then locators last then exposed url  
// authors, year, title<texturl>, enclosing_title,publisher or source
//  issue,  location ( page, line, etc) , texturl, infourl  , id (doi)  
/
// person needs to disambiguate, allow contact, and provide relevant
// biography intro for the citation being made 
// person: "author", year ( of information ), location, affiliation or employer, title, background/alumni ed,  id(orcid, ssno, driver etc lol ). 
//
// bom : description ( noun and adjectives, title like), mfg(city), domain, model number , exceptional features, subs allowed, quantity, specs url, purchase url_n 
 
return 0;
} // format 


private:


}; // mjm_bibtex_format_entry


template <class Tr>
class mjm_bibtex_format 
{
 typedef mjm_bibtex_format Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
 typedef mjm_bibtex_entry<Tr> BibEntry;
 typedef mjm_bibtex_entry_map<Tr> BibMap;
typedef mjm_bibtex_render<Tr> Renderer;
typedef mjm_ragged_table Ragged;
public:
mjm_bibtex_format() {}
~mjm_bibtex_format() {}
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


// MEMBERS



}; // mjm_bibtex_format

//////////////////////////////////////////////

template <class Tr>
class mjm_bibtex_format_map : public std::map<typename Tr::StrTy, mjm_bibtex_format< Tr > >  
{
 typedef mjm_bibtex_format_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_bibtex_format< Tr> >   Super;
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
mjm_bibtex_format_map() {}
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

}; // mjm_bibtex_format_map




////////////////////////////////////////////
#ifdef  TEST_MJM_BIBTEX_FORMAT
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
typedef tester_< mjm_bibtex_format <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_BIBTEX_FORMAT "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_bibtex_format<Tr>  Myt;
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

#endif // MJM_BIBTEX_FORMAT_H__ 
