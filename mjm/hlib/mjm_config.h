
#ifndef MJM_CONFIG_H__
#define MJM_CONFIG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


#include "mjm_globals.h"
#include "mjm_io.h"


class mjm_config
{
private:
typedef mjm_config Myt;
typedef std::runtime_error CmdEr;

public:
typedef mjm_generic_traits Tr;
typedef Tr::ChTy ChTy;
typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;
typedef Tr::StrTy StrTy;
typedef unsigned int FlagTy;
typedef Tr::IdxTy  IdxTy;
typedef Tr::ErTy  ErTy;

typedef std::map<StrTy, StrTy> ConfigTy;
typedef ConfigTy::iterator CoItor;
typedef ConfigTy::const_iterator CocItor;
const bool verbose;
ConfigTy m_config;

IdxTy m_max_line;
IdxTy m_max_fields;
StrTy m_default; 
public:
mjm_config(): verbose(true), m_max_line(1<<16),m_max_fields(16),m_default("") {}


~mjm_config()
{
// don't want to really risk killing
// on bad files, 
// std::cout<<MM_MARK<<" dtor called at "<<" " <<" lines/sec="<<CRLF; 
// counting on dtor to flush files and close

}
// this doesn't like const

// const StrTy & get(const StrTy &key) const {  return m_config[key]; }
const StrTy & get(const StrTy &key) const 
{  
CocItor ci=m_config.find(key);
if ( ci!=m_config.end()) return (*ci).second;
return m_default; 
//return m_config[key]; 

}
void put(const StrTy & key, const StrTy &v) { m_config[key]=v; }
void load(const StrTy & fn) { load(fn.c_str()); }
void load(const ChTy * fn )
{
IsTy * is = mjm_io::get_in(fn);
if ( is==NULL) throw  ErTy("not found");
load(*is);


}
void load(IsTy & is)
{
ChTy line[m_max_line];
ChTy * fields[m_max_fields];
if ( verbose) osx()<<MM_MARK<<" loading config "<<CRLF;

while (!is.eof()&&is.good())
{
is.getline(line,m_max_line);
if ( verbose) osx()<<MM_MARK<<" loading config  have line "<<strlen(line)<<" "<< line<<CRLF;
if ( line[0]=='#' ) continue;
if ( line[0]=='/' ) continue;
IdxTy field=0;
const IdxTy len= ::strlen(line);
mjm_strings::parse_fields<true>(line,len, fields,field,m_max_fields);

if ( field<2) continue;
m_config[StrTy(fields[0])]=StrTy(fields[1]);

} // while 

if ( verbose) osx()<<MM_MARK<<" dun loading config  "<<CRLF;

}


void dump(OsTy & os)
{
CoItor ii= m_config.begin();
CoItor ie= m_config.end();
while(ii!=ie)
{

os<<(*ii).first<<" "<<(*ii).second<<CRLF;
++ii; 

} // while




}







};  // mjm_io


#endif

