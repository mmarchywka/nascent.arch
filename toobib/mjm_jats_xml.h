#ifndef MJM_JATS_XML_H__
#define MJM_JATS_XML_H__

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
// copied from bomtex and xref_json
#include "mjm_pawnoff.h"
#include "mjm_collections.h"
#include "mjm_wovdb.h"
#include "mjm_strings.h"
#include "testHTML.h"
//#include <../mjsonu/mjsonu.h>

#include "mjm_hier_two.h"

//#include <../mjsonu/mjsonu.h>

#include "mjm_pawnoff.h"
#include "mjm_collections.h"
#include "mjm_wovdb.h"
#include "mjm_strings.h"




// Tue Feb  8 07:46:50 EST 2022
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_jats_xml   
//g++  -Wall -std=gnu++11 -DTEST_MJM_JATS_XML -I. -I../../mjm/hlib -I../../mjm/num  -I /usr/include/libxml2  -gdwarf-3 -O0  -x c++ mjm_jats_xml.h  -lpthread -lreadline -lxml2 -Wno-unused-variable -Wno-unused-function
// g++  -Wall -std=gnu++11 -DTEST_MJM_JATS_XML -I. -I../../mjm/hlib -I../../mjm/num  -I /usr/include/libxml2  -gdwarf-3 -O0  -x c++ mjm_jats_xml.h  -lpthread -lreadline -lxml2
// g++  -Wall -std=gnu++11 -DTEST_MJM_JATS_XML -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_jats_xml.h  -lpthread -lreadline

// -I /usr/include/libxml2 -I/home/documents/cpp/pkg/include  -gdwarf-3 -O0  -x c++ mjm_bomtex_json.h  -lpthread -lreadline -lxml2


template <class Tr>
class mjm_jats_xml 
{
 typedef mjm_jats_xml Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;


typedef mjm_pawnoff<Tr> Hand;
typedef typename Hand::blob Blob;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef Ragged::Word Word;
typedef mjm_ragged_cursor RaggedCursor;
typedef mjm_wovdb<Tr,StrTy> Tdb;
typedef mjm_hier_two<Tr> HierUtil;


//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_jats_xml() {}
~mjm_jats_xml() {}

StrTy jats_xform( const StrTy & fn, const IdxTy  flags) const
{
Ragged r;
IdxTy rc=load_ragged(r,fn,0);
std::map<StrTy,StrTy> m;
IdxTy rca=assemble(m,r,0);
StrTy s;
IdxTy rcf=format(s,m,0);
MM_ERR(MMPR(fn)<<MMPR4(s,m.size(),r.size(),rc))
return s;
} // crossref_xform


IdxTy load_ragged(Ragged & r, const StrTy & fn, const IdxTy  flags) const
{
//MM_MSG(" using integrated json for xref may not work doh ")
std::ifstream is(fn);
IdxTy rc=testHTML::parse_stream(r,is,0);
MM_ERR(r.dump())
return 0;
}


template <class Tm > IdxTy assemble(Tm & m, const Ragged & r,  const IdxTy  flags) const
{
return Assemble(m,r,flags);
}
template <class Tm > IdxTy format(StrTy & s, Tm & m,  const IdxTy  flags) const
{
return Format(s,m,flags);
}




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
// this seems to be hierarchial and may generate multiple
// entries. The raged could be chopped to one article v and h dirs
// or they can be partsed into a map of Tm;s.

/*

mjm_jats_xml.h340 0|1|html|3|body|4|article|5|xmlns:xsi|6|text|http://www.w3.org
/2001/XMLSchema-instance
1|1|html|3|body|4|article|7|xmlns:mml|8|text|http://www.w3.org/1998/Math/MathML
2|1|html|3|body|4|article|9|xmlns:xlink|10|text|http://www.w3.org/1999/xlink
3|1|html|3|body|4|article|11|dtd-version|12|text|1.1
4|1|html|3|body|4|article|13|xml:lang|14|text|en
5|1|html|3|body|4|article|15|front|16|journal-meta|17|journal-id|18|text|authore
a
6|1|html|3|body|4|article|15|front|16|journal-meta|19|publisher|20|publisher-nam
e|21|text|Authorea
7|1|html|3|body|4|article|15|front|22|article-meta|23|article-id|24|pub-id-type|
25|text|doi
8|1|html|3|body|4|article|15|front|22|article-meta|23|article-id|26|text|10.2254
1/au.159103680.00306295


22|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|46|contrib|
56|address|58|institution|59|text|Universidade Federal Fluminense
23|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|46|contrib|
60|text|
        
24|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|61|contrib|
62|contrib-type|63|text|author
25|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|61|contrib|
64|corresp|65|text|no
tmm[l[kid+3][l[len-4]]=v;

80|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|181|contrib|182|contrib-type|183|text|author
81|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|181|contrib|184|corresp|185|text|no
82|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|181|contrib|186|name|187|surname|188|text|Fava
83|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|181|contrib|186|name|189|given-names|190|text|Claudia Del
84|1|html|3|body|4|article|15|front|22|article-meta|30|contrib-group|181|contrib|191|address|192|text|




*/


template <class Tm > IdxTy Assemble(Tm & _m, const Ragged & r,  const IdxTy  flags) const
{
typedef std::map<StrTy, Tm > Tmx;
Tmx tmx;

const IdxTy aid_pos=4; // aritcle node id pos
const IdxTy kid_pos=11; // key id  pos
const IdxTy sz=r.size();
typedef std::map< StrTy, std::map< StrTy, std::vector<StrTy >  > > TwoMap;
std::map<StrTy, TwoMap>  tmmm;
for(IdxTy i=0; i<sz; ++i)
{
const Line & l=r[i];
const IdxTy len=l.size();
if (len<2) continue; // no kvp 
MM_ERR(MMPR2(len,l[len-1]))
if (len<=kid_pos) continue; // no kvp 
const StrTy & v=l[len-1];
Tm & m=tmx[l[aid_pos]];
TwoMap& tmm=tmmm[l[aid_pos]];;
const StrTy & kid=l[kid_pos];
MM_ERR(MMPR2(l[aid_pos],kid))
if (kid=="journal-id") m["journal"]=(v);
else if (kid=="publisher") m["publisher"]=(v);
//7|1|html|3|body|4|article|15|front|22|article-meta|23|article-id|24|pub-id-type|
//8|1|html|3|body|4|article|15|front|22|article-meta|23|article-id|26|text|10.2254
else if (kid=="contrib-group")
{  
//StrTy q=l[kid_pos+2];
//if (q=="name") q=l[kid_pos+4];
//tmm[l[kid_pos+3]][q].push_back(v);
MM_ERR(MMPR3(l[kid_pos+1],l[len-4],v))
tmm[l[kid_pos+1]][l[len-4]].push_back(v);
}
else if (kid=="article-id")  tmm[l[kid_pos]][l[kid_pos+2]].push_back(v);
// TODO this code may be bombinf if the data is bad... 
//else tmm[l[kid_pos+3]][l[len-4]].push_back(v);
else if ( len>(kid_pos+3)) tmm[l[kid_pos+3]][l[len-4]].push_back(v);
} // i main line loop 
//////////////////////////////////////////
MM_LOOP(oo,tmx)
{
Tm & m=tmx[(*oo).first];
TwoMap& tmm=tmmm[(*oo).first];
bool used=false;
// now dismember the misc map in tmm...
auto& aidm=tmm["article-id"];
if (aidm.size())
{
const StrTy & t=aidm["pub-id-type"][0];
if (t.length()){  m[t]=aidm["text"][0]; } 
} // aidm.size
//if (!used) 
if ( tmm["abstract"].size())
{ // this  doesnt  do see the  later in loop  
auto& aidm=tmm["abstract"];
const StrTy & t=aidm["p"][0];
if (t.length()) m["abstract"]=t; // aidm["text"][0];
} // abstract
//else 
if (!used) { 
MM_LOOP(ii,tmm)
{
auto & em=(*ii).second;
if (em.find("contrib-type")!=em.end())
{
const StrTy ty=em["contrib-type"][0];
if (ty=="author")
{
MM_LOOP(ff,em)
{
MM_ERR(" fuc "<<MMPR2((*ff).first,(*ff).second[0]))
}
Vadd(em,m,"author","surname","given-names",0);
/*
StrTy last;
MM_LOOP(jj,em["surname"]) last+=(*jj);
 StrTy first;
MM_LOOP(jj,em["given-names"]) first+=(*jj);
if (m["author"].length()) m["author"]+=" and ";
m["author"]+=first+" "+last;
*/
} // author
else 
Vc(em,m,0);
//MM_LOOP(jj,em) m[(*jj).first]+=(*jj).second[0];
} // contrib-type
else
{
Vc(em,m,0);
//MM_LOOP(jj,em) m[(*jj).first]+=(*jj).second[0];

}
} // ii
} // else 
m["abstract"]=m["p"];
Mv(m,"title","article-title");
Mv(m,"abstract","p");
m["type-name"]="article";
m["name-name"]="xxx";


} // oo 

MM_ERR(MMPR(tmx.size()))
if (tmx.size()>0) _m=(*tmx.begin()).second;
return 0;
}


template <class Tm,class Te> 
IdxTy Vc(Te & em, Tm & m, const IdxTy flags)
const 
{
MM_LOOP(jj,em)
{
const auto k=(*jj).first;

 const auto ve=(*jj).second;
MM_LOOP(kk,ve) {
if (m[k].length()) m[k]+=", ";
 m[k]+=(*kk);
} // kk 
}
return 0; 
} // Vc


template <class Tm> 
IdxTy Mv( Tm & m, const StrTy & d, const StrTy & s, const IdxTy flags=0)
const 
{
auto ii=m.find(s);
if (ii==m.end()) return 1;
m[d]=m[s];
m.erase(ii);
return 0;
} // Mv


template <class Tm,class Te> 
IdxTy Vadd(Te & em, Tm & m, const StrTy & d, const StrTy & dx,const StrTy & s, const IdxTy flags)
const 
{
StrTy last;
 StrTy first;
MM_LOOP(jj,em[s]) last+=(*jj);
if (dx.length()) { MM_LOOP(jj,em[dx]) first+=(*jj); } 
if (m[d].length()) m[d]+=" and ";
if ( first.length()) m[d]+=first+" ";
m[d]+=last;

return 0;
}


template <class Tm > IdxTy Format(StrTy & s, Tm & m,  const IdxTy  flags) const
{
Ss ss;
ss<<"@"<<m["type-name"]<<"{"<<m["name-name"]; // <<","<<CRLF;
MM_LOOP(ii,m)
{
const StrTy & key=(*ii).first;
if (key=="type-name") continue;
if (key=="name-name") continue;

ss<<","<<CRLF;
const StrTy & v=(*ii).second;
ss<<key<<" = {"<<v<<"} ";

}
ss<<CRLF<<"}"<<CRLF;

s=ss.str();
return 0;
}// Format 

// MEMBERS
mutable Hand m_hand;



}; // mjm_jats_xml

//////////////////////////////////////////////

template <class Tr>
class mjm_jats_xml_map : public std::map<typename Tr::StrTy, mjm_jats_xml< Tr > >  
{
 typedef mjm_jats_xml_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_jats_xml< Tr> >   Super;
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
mjm_jats_xml_map() {}
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

}; // mjm_jats_xml_map




////////////////////////////////////////////
#ifdef  TEST_MJM_JATS_XML
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
typedef tester_< mjm_jats_xml <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;

typedef mjm_ragged_table Ragged;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_JATS_XML "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_jats_xml<Tr>  Myt;
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
else if (cmd=="load") {
// x.load(li.words(),1); 
Ragged r;
IdxTy rc=x.load_ragged(r,cip.p1,0);
std::map<StrTy,StrTy> m;
IdxTy rca=x.assemble(m,r,0);
StrTy s;
IdxTy rcf=x.format(s,m,0);
MM_MSG(r.dump())
MM_MSG(MMPR4(s,m.size(),r.size(),rc))

//IdxTy load_ragged(Ragged & r, const StrTy & fn, const IdxTy  flags) const
}


//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_JATS_XML_H__ 