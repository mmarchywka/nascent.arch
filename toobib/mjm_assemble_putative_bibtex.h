#ifndef MJM_ASSEMBLE_PUTATIVE_BIBTEX_H__
#define MJM_ASSEMBLE_PUTATIVE_BIBTEX_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_pawnoff.h"
#include "mjm_collections.h"
#include "mjm_wovdb.h"
#include "mjm_read_buffer.h"
#include "mjm_misc_parse.h"
#include "mjm_bibtex_entry.h"
#include "mjm_bibtex_fields.h"
#include "mjm_strings.h"



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

// caopy from mjm_scrape_xform  2021-01-05
// g++  -Wall -std=gnu++11 -DTEST_MJM_ASSEMBLE_PUTATIVE_BIBTEX -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_assemble_putative_bibtex.h  -lpthread -lreadline

template <class Tr>
class mjm_assemble_putative_bibtex 
{
 typedef mjm_assemble_putative_bibtex Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;


typedef mjm_pawnoff<Tr> Hand;
typedef typename Hand::blob Blob;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef Ragged::Word Word;
typedef mjm_wovdb<Tr,StrTy> Tdb;
typedef mjm_read_buffer<Tr> Rdbuf;

typedef std::map<StrTy, StrTy> XformMap;

typedef mjm_bibtex_entry<Tr> BibEntry;
typedef mjm_bibtex_entry_map<Tr> BibMap;

typedef mjm_bibtex_fields<Tr> Fields;
//typedef mjm_var_parse<tr> CharTable;
typedef mjm_strings StrUtil;
class _ResultType
{
typedef _ResultType Myt;
public:
_ResultType() {Init(0); }
_ResultType(const IdxTy n) {Init(n); }
Myt & operator=(const IdxTy n ) { m_n=n;  return *this; }
Myt & operator=(const BibEntry &  n ) { m_be=n;  return *this; }
operator IdxTy() const { return m_n; } 
void errors(const IdxTy e ) { m_errors=e; } 
IdxTy  errors( ) const  { return m_errors; } 
BibEntry & be() { return m_be; } 
BibMap & bm() { return m_bm; } 
const BibEntry & bib() const  { return m_be; } 
StrTy & latex_output() { return m_latex_output; } 
const StrTy & rendered()const  { return m_rendered; } 
void  rendered(const StrTy & s) {m_rendered=s; } 
StrTy & validation_fn() { return m_validation_fn; } 
void  validation_fn(const StrTy & s) {m_validation_fn=s; } 


StrTy dump(const IdxTy flags=0) const
{
Ss ss;
ss<<MMPR2(m_n,m_errors);
ss<<MMPR4(m_be.name(),m_be.type(),m_be.size(),m_be.errors());
ss<<MMPR(m_latex_output);
return ss.str();
}

private:

void Init(const IdxTy n) {m_n=n; m_errors=0;}

IdxTy m_n,m_errors;
BibEntry m_be;
BibMap  m_bm;
StrTy m_latex_output;
StrTy m_rendered,m_validation_fn;
}; // _Result_Type

typedef _ResultType ResultType;
typedef ResultType AccessTy;


/*
Examine a newly downloaded candidate bibtex
entry and validate with old and new approaches.
Clean up and insert fields etc. 


*/
class _AssembleOpt
{


}; // _AssembleOpt
typedef _AssembleOpt AssembleOpt;

public:
typedef AccessTy return_code;
typedef Blob blob_type;
typedef AssembleOpt assem_opt;
mjm_assemble_putative_bibtex() {Init();}
mjm_assemble_putative_bibtex(const Ragged & r) {Init(r); }
/*

// add trailing cr, dos2unix, maybe indent insert comments etc 
// add k-v pairs to bibtex OR as leading comments. 
AccessTy final_assembly(Blob & d, const Blob &s,const  Keys & k, const Values & v, const IdxTy flags)
 
AccessTy canonicalized(Blob & d, const Blob & s,  const IdxTy flags)
// not a dup name etc. 
AccessTy ok_for_file(Blob & b, const StrTy & fn, const IdxTy flags)

*/
// normally append but maybe sort or something odd. 
AccessTy add_to_file(Blob & b, const StrTy & fn, const IdxTy flags)
{return AddToFile(b, fn,  flags); }

AccessTy safe_final_assembly(Blob & d, const Blob & s, const assem_opt & ao)
const {
return SafeFinalAssembly(d,s,ao);
}
AccessTy light_clean(const StrTy & fn, const StrTy & fntemp, const IdxTy flags) const
{ return LightClean(fn,fntemp,flags); } // light_clearn 

// completely parse and write back out with changes 
AccessTy rewrite_parse(const StrTy & fn, const StrTy & fntemp, const IdxTy flags) const
{
return  ReWriteParse( fn, fntemp, flags); 
}
IdxTy notate_edit(BibEntry & be, const StrTy & field, const StrTy & s) const
{ return NotateEdit(be,field,s); } 
// old style and try to parse if with new code. 
// this "works" with a blank blob. 
AccessTy validate(Blob & b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const
{ return Validate(b,ofn,temp,flags); } 

AccessTy parse(Blob & b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const
{ return Parse(b,ofn,temp,flags); } 



StrTy dump(const IdxTy flags=0) const  { return Dump(flags); }




private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0)const  {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

AccessTy LightClean(const StrTy & fn, const StrTy & fntemp, const IdxTy flags) const
{

AccessTy rc=0;
StrTy cmd;
IdxTy lc=Cmd(cmd,"lightclean",0);
const StrTy fnq=StrTy("\"") + fntemp+"\" \""+fn+"\"" ;
//cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; cp \""+fntemp+"\" \""+fn+"\"" ;
cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; cp "+fnq ;
IdxTy rt= Execute(0,0, cmd); 
IdxTy lc2=Cmd(cmd,"lightcleanchange",0);
cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; diff "+fnq+" ;  cp "+fnq ;
Blob data,out,err;
IdxTy rt2=m_hand.fileio(out,err,data,cmd,3);
MM_ERR(" light cleaning change err  "<<MMPR2(cmd,StrTy(err)))
MM_ERR(" light cleaning change out  "<<MMPR(StrTy(out)))
rc.rendered(StrTy(out));

rc=rt;
return rc;
}
IdxTy ReWriteKvp(BibEntry & be,  const StrTy & dn, const StrTy & sn,
const StrTy & cmd,const IdxTy flags) const
{
const bool forces=Bit(flags,0);
const bool forced=Bit(flags,1);
const bool nonzed=Bit(flags,2);
StrTy d,s;
s=be.get(sn);
MM_ERR(MMPR4(d,s,sn,dn)<<MMPR3(cmd,forces,forced));
if (!forces) if (s.length()==0) return 1; 
m_hand.cmd_exec(d,s,cmd,4);
MM_ERR(MMPR4(d,s,sn,dn)<<MMPR(cmd));
if (!forced) if (d.length()==0) return 2; 
if (nonzed) if (atoi(d.c_str())==0) return 2; 
be.replace(dn,d,0);
Ss sss;
sss<<"ReWriteKvp"<<MMPR3(dn,sn,flags);
NotateEdit(be,dn,sss.str()); 
//IdxTy NotateEdit(BibEntry & be, const StrTy & field, const StrTy & s) const


return 0; 
} //ReWriteKvp 

IdxTy FixBeKvp(BibEntry & be, const StrTy & dn ,  const StrTy & sn,
const StrTy & cmd,const IdxTy flags) const
{
const bool force=Bit(flags,0);
StrTy d;
const StrTy & s=be.get(sn);
m_hand.cmd_exec(d,s,cmd,4);
Ss ss;
ss<<__FUNCTION__ <<MMPR4(s,cmd,d,dn);
MM_ERR(ss.str())
if (!force) if (d.length()==0) return 1;
be.replace(dn,d,0);
NotateEdit(be,dn,ss.str());
//be.add("X_TooBib", ss.str());
return 0; 
} // FixBeKvp
IdxTy NotateEdit(BibEntry & be, const StrTy & field, const StrTy & s) const
{
be.add("X_TooBib", field+": "+ s);
return 0; 
}

template <class Tv> 
IdxTy CheckInOrder(StrTy & s, StrTy & slin, BibEntry&  be, const Tv & list, const IdxTy flags)
const
{
StrTy cr="\n";
MM_LOOP(ii,list) { 
const StrTy & x=be.get(*ii);
MM_ERR(MMPR2(x,(*ii)))
s+=cr+(x); slin+=" "+(x);
} // ii 

return 0;
} // CheckInOrder


AccessTy ReWriteParse(const StrTy & fn, const StrTy & fntemp, const IdxTy flags) const
{
const bool make_lc=Bit(flags,0);
// see  comments below 
const bool adhoc_date_munge=!true;
const bool canonical_author_list=true;
const bool fix_non_ascii=true;
const bool field_fix_dates=true;
const bool remove_longs=true;
const bool fix_all_names=!true;
AccessTy rc=0;
StrTy cmd;
BibMap bm;
bm.parse(fn);
if (bm.size()==0)
{
MM_ERR(MMPR(bm.size()))
return rc;
}
//BibEntry bef,bei;
// TODO also have not made uniform lower case etc yet. 
BibEntry  be=bm[0];
if (remove_longs)
{
const IdxTy lenlim=1024*2;
std::vector<StrTy>  v;
be.all_keys(v,0);
MM_LOOP(ii,v)
{
be.remove_long(*ii,lenlim);
} // ii 
} // remove_longs
if (make_lc) be.make_keys_canonical(0);
if (make_lc) be.trim_values(0);
Ragged r;
be.format_rag(r,1);
StrTy name=be.name();
// TODO FIXME this needs to be used ISNTEAD of "date" whic
// takes a 4 digit year as a TIME and returns today's  date
//  
if (field_fix_dates)
{
const IdxTy datelen=be.get("date").length();
const IdxTy ylen=be.get("year").length();
const IdxTy mlen=be.get("month").length();
const IdxTy dlen=be.get("day").length();
{ MM_SSF(MMPR4(datelen,ylen,mlen,dlen)) MM_ERR(ss.str()) }

typedef typename Fields::parsed_dates_class Pdc;
//typedef typename Pdc::value_type Date;

Pdc pdc;
Fields x;
StrTy cr="\n";

/*
date = {2016},
date-created = {2016-05-10T17:01:22Z},
date-deposited = {2017-02-022017-02-02T21:52:04Z},
date-indexed = {2021-12-16T14:03:18Z},
date-issued = {2016},
date-published-print = {2016},


*/
StrTy slin="";
StrTy s=be.get("date");
slin+=s;
typedef std::vector<StrTy> Vlist;
Vlist vdate;
// ? This was used somewhere and may have worked???? 
//vdate.push_back("journal-issue");
vdate.push_back("date-journal-issue");
 vdate.push_back("publication_date");
 vdate.push_back("date_published");
 vdate.push_back("date-published-print");
 vdate.push_back("publication-date");
 vdate.push_back("date-posted");
 vdate.push_back("datemodified");
 vdate.push_back("datepublished");
 vdate.push_back("online_date");
 vdate.push_back("date-published-online");
 vdate.push_back("dateaccepted");

//IdxTy CheckInOrder(StrTy & s, StrTy & slin, BibEntry&  be, const Tv & list, const IdxTy flags)
 CheckInOrder( s, slin,  be, vdate, 0);
/*
s+=cr+be.get("journal-issue");
slin+=" "+be.get("journal-issue");

s+=cr+be.get("publication_date");
slin+=" "+be.get("publication_date");

s+=cr+be.get("date_published");
slin+=" "+be.get("date_published");
s+=cr+be.get("date-published-print");
slin+=" "+be.get("date-published-print");

s+=cr+be.get("datemodified");
slin+=" "+be.get("datemodified");

s+=cr+be.get("datepublished");
slin+=" "+be.get("datepublished");


s+=cr+be.get("online_date");
slin+=" "+be.get("online_date");

s+=cr+be.get("date-published-online");
slin+=" "+be.get("date-published-online");

s+=cr+be.get("dateaccepted");
slin+=" "+be.get("dateaccepted");
*/

// this should return a result for each line 
// this is just munving them into one lines doh.. 
// this could be right now.. 
x.parse_dates(pdc,s,0);
MM_LOOP(ii,pdc) { MM_ERR(MMPR((*ii).dump(0))) } 
//StrTy ndate=pdc[0].lexi();
MM_LOOP(ii,pdc)
{
const auto & bd=(*ii); // pdc[0];
MM_ERR(MMPR(bd.dump()))
if (bd.is_zed()) continue;
//StrTy ndate=pdc[0].best();
StrTy ndate=bd.best();
//StrTy sal=x.canon_author_list(pnc,0);
//be.replace("date",ndate,0);
be.only("date",ndate);
// this is f- up because it is a munge of lines... 
// but it is due to slin which is intended.... gets everything. 
be.replace("date_orig",slin,0);
NotateEdit(be,"date",ndate);
// except that does not  work... 
// this is postponed to let the date expander do it, rewrite_ymd
if (ylen==0)  if (bd.has_year()){ be.add("year",bd.year()); 
NotateEdit(be,"year",bd.year()+StrTy(",  infield_fix_dates ")  );
if (mlen==0)  if (bd.has_month()){ be.add("month",bd.month()); 
NotateEdit(be,"month",bd.month()+StrTy(",  infield_fix_dates ")  );

if (dlen==0)  if (bd.has_day()){ be.add("day",bd.day()); 
NotateEdit(be,"day",bd.day()+StrTy(",  infield_fix_dates ")  );


}}


}

break;
//const IdxTy datelen=be.get("date").length();
//const IdxTy ylen=be.get("year").length();
//const IdxTy mlen=be.get("month").length();
//const IdxTy dlen=be.get("day").length();
// first find a good date or creation time
// determine if y,m, and d fields exist and try to populate
} // ii 

} // field_fix_dates



if (adhoc_date_munge)
{
// first want  year which may bemissing, parse the date....
const IdxTy datelen=be.get("date").length();
const IdxTy ylen=be.get("year").length();
const IdxTy mlen=be.get("month").length();
const IdxTy dlen=be.get("day").length();
{ MM_SSF(MMPR4(datelen,ylen,mlen,dlen)) MM_ERR(ss.str()) }

if ((datelen==0)&&(ylen==0))
{
// first try to find anything close, ZZ
std::map<StrTy,IdxTy> mm;
be.matching_keys(mm,"date",0);
if (mm.size()==0)
{
FixBeKvp( be, "date_downloaded" ," ", "date \"+%Y-%m-%d\"",0);
be.matching_keys(mm,"date",0);

} // no dates

MM_LOOP(ii,mm)
{
// should probably check for earliest... 
MM_ERR(MMPR2(__FUNCTION__,(*ii).first))
FixBeKvp( be, "date" ,(*ii).first, "date -f -  \"+%Y-%m-%d\"",0);

} // mm 

} // datelen, ylen == 0

const IdxTy datelen2=be.get("date").length();
const IdxTy ylen2=be.get("year").length();
const bool rewrite_ymd=(datelen2!=0)&&(ylen2==0); 
if (rewrite_ymd)
{
//FixBeKvp( be, "date" ,"", "date \"+%Y-%m-%d\"",0);
// TODO FIXME this is bad with date -f as a 4 digit numbers
// is a TIME and it operates on TODAY's date  
ReWriteKvp(be,"year","date", "date -f - \"+%Y\"",4);
ReWriteKvp(be,"month","date", "date -f - \"+%m\"",7);
ReWriteKvp(be,"day","date", "date -f - \"+%d\"",7);
const IdxTy ylen3=be.get("year").length();
// "date" does not work well with yyyy-mm
if (ylen3==0) if (datelen2==7)
{
IdxTy y4=atoi(be.get("date").c_str());
// may not work well for antiquities ... :)
if (y4>0) { Ss ss; ss<<y4;  be.add("year",ss.str()); } 
y4=atoi(be.get("date").c_str()+5);
if (y4>0) { Ss ss; ss<<y4;  be.add("month",ss.str()); } 
NotateEdit(be,"year","extract from date leng== 7");
} // 
} // rewrite_ymd
} // adhoc_date_munge

//FixBeKvp( be, "date_downloaded" ," ", "date \"+%Y-%m-%d\"",0);
if (be.get("urldate").length()<1)
{
FixBeKvp( be, "urldate" ," ", "date \"+%Y-%m-%d\"",0);
}


if (be.get("author").length()<2)
{
Ss ss; ss<<__FUNCTION__ <<MMPR(be.get("author"));
NotateEdit(be,"author",ss.str());
be.move("author","authors",1+4);
}

if (be.get("author").length()<2)
{
Ss ss; ss<<__FUNCTION__ <<MMPR(be.get("author"));
NotateEdit(be,"author",ss.str());
be.move("author","name",1+4);
}



if (be.get("publisher").length()<2)
{
const StrTy dest="publisher";
const StrTy s="site_name";
Ss ss; ss<<__FUNCTION__ <<MMPR2(be.get(s),be.get(dest));
NotateEdit(be,dest,ss.str());
be.move(dest,s,1+4);
}

if (be.get("journal").length()<2)
{
const StrTy dest="journal";
const StrTy s="publisher";
Ss ss; ss<<__FUNCTION__ <<MMPR2(be.get(s),be.get(dest));
NotateEdit(be,dest,ss.str());
be.move(dest,s,1+4);
}

// this f- up pages like wikipedia if there is a better title... 
const StrTy & head=be.get("headline");
const StrTy & ctitle=be.get("title");
if (ctitle.length()==0) if ( head.length()>10){
NotateEdit(be,"title",be.get("title"));
 be.only("title",head); }

//digitalcommons, ,
//dissertation_institution = {The University of Maine},
//dissertation_name = {Doctor of Philosophy (PhD)},

const StrTy& xdi=be.get("dissertation_institution");
if (xdi.length()>3)
{
const StrTy x=be.get("dissertation_name");
const IdxTy pos=mjm_strings::indexOf(x.c_str(),"PhD");
if (pos!=(~0)) be.type("phdthesis");
// this should happen in guessdigitalcommon but failing.. 
const StrTy xsc=be.get("school");
if (xsc.length()==0) be.add("school",xdi);

} // inst

// single author and year is also a collision concern now.. 
//if (name.length()<5)
#if 0 
const bool use_good_names=true;
if (name.length()<12)
{
MM_ERR(" using klged name thin "<<__FUNCTION__)
StrTy newname="";
if (! use_good_names) {
// make a good name
const StrTy year=be.get("year");
newname=year;
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l=r[i];
if (l.size()<2) continue;
const bool title= (l[0]=="title") ;
const bool author=((l[0]=="author"));
//if ((l[0]=="author") || (l[0]=="title")) 
if (title||author) 
{
for(IdxTy j=1; j<l.size(); ++j) 
{
const StrTy & w=l[j];
if ( w.length()<2) continue;
// should use a low P dictionary 
if (title) if ( w.length()<5) continue;
// with the code splitter this should pick up "and" and comma
// allowing for last names only... 
if (author) if ((w.c_str()[1]<'a')||(w.c_str()[2]>'z')) continue;
if (author) if ((w=="and")) continue;
if (newname.length()==0) newname=w;
else newname=newname+"_"+w;
if (newname.length()>36) break; 
} // j 
} /// author 
if (newname.length()>36) break; 
} // i 
} // use_good_names
//}// length 
#endif
const bool missing=(name=="missing");
const bool bad_name=missing||(name.length()<5);
const bool fix_name_later=true;
if (!fix_name_later)
{
if (fix_all_names||bad_name)
{
be.make_name();
name=be.name();
//if (newname.length()>4) be.name(newname); 
if (name.length()<6)  
{
StrTy doi=be.get("doi");
if ( doi.length()>8) be.name(doi); 
} // newname
} // name 
} // name later

if (fix_non_ascii)
{
Fields f;
IdxTy changes=0;
auto  m=be.map();
auto m2=m; // just to get auto to wok  
m2.clear(); // 
MM_LOOP(ii,m)
{
const StrTy & k =(*ii).first;
const StrTy klc=StrUtil::fancy_to_lower(k);
const bool origk =StrUtil::ends_with(klc.c_str(),"_orig")
 || StrUtil::ends_with(klc.c_str(),"_as_rcvd");
const bool try_to_fix=(klc=="author") ||(klc=="title");
const bool verbatim_key=origk&& !try_to_fix;
const auto orig=(*ii).second;
bool some_changed=false;
MM_LOOP(jj,(*ii).second)
{
//StrTy x;
StrTy x=(*jj);
// this should be done after html fixes... 
//IdxTy n=f.utf8_to_latex(x,(*jj),0);
//MM_ERR("FIX_NON_ASCII "<<MMPR3(n,x,(*jj)))
StrTy y;
// don't do this with names... the braces etc make
// parsing hard...
IdxTy mm=0;
// this makes a mess of many fields 
bool skip_html=verbatim_key; // false;
// https://pdfs.semanticscholar.org/0ee8/33694e13ddc65109305cc0a59c7545f8f9ce.pdf

//skip_html=(strncmp(k.c_str(),"author",6)==0);
skip_html|=(strncmp(k.c_str(),"X_",2)==0);
skip_html|=(strncmp(k.c_str(),"x_",2)==0);
if (skip_html ) {y=x; }
else { mm=f.html_to_latex(y,x,0); }

//IdxTy n=f.utf8_to_latex(y,(*jj),0);
StrTy ynew;
IdxTy n=f.utf8_to_latex(ynew,y,0);
y=ynew;
//IdxTy nll=f.latex_to_latex(ynew,(*jj),0);
IdxTy nll=f.latex_to_latex(ynew,y,0);
if ( nll) y=ynew;
//MM_ERR("FIX_HTML_ASCII "<<MMPR4(n,mm,x,y))
if (n== ~0) continue; // ignore invalids.. 
if (mm== ~0) continue; // ignore invalids.. 
if (n+mm+nll) { changes+=n+mm+nll; (*jj)=y;
some_changed=true;
}
} // jj
if (some_changed)
{
m2[k+"_as_rcvd"]=orig;
} // some_changed

} // ii 
if (changes) { 
MM_LOOP(ii,m2) 
{
MM_LOOP(jj,(*ii).second){ m[(*ii).first].push_back((*jj)); } 
} // ii 
 be.clear_map(); be.add(m); }

} // fix_non_ascii




if (canonical_author_list)
{
typedef typename Fields::parsed_names_class Pnc;
Pnc pnc;
Fields x;
const StrTy ali=be.get("author");
x.parse_names(pnc,ali,0);
StrTy sal=x.canon_author_list(pnc,0);
be.replace("author",sal,0);
be.replace("author_orig",ali,0);
NotateEdit(be,"author",sal);

}
/// refernces normally too bi
StrTy reff=be.get("references");
if (reff.length()>100)
{
be.remove("references",0);
NotateEdit(be,"references","deleted for space " );
}
reff=be.get("reference");
if (reff.length()>100)
{
be.remove("reference",0);
NotateEdit(be,"reference","deleted for space " );
}

if (fix_name_later)
{
if (fix_all_names||bad_name)
{
be.make_name();
name=be.name();
//if (newname.length()>4) be.name(newname); 
if (name.length()<6)  
{
StrTy doi=be.get("doi");
if ( doi.length()>8) be.name(doi); 
} // newname
} // name 
} // name later



// make all the keys uniform lower case 


// change authors to author

// Final date checks 


StrTy x= be.format();
Blob b;
b=x;
b.save(fn);

#if 0 
IdxTy lc=Cmd(cmd,"lightclean",0);
const StrTy fnq=StrTy("\"") + fntemp+"\" \""+fn+"\"" ;
//cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; cp \""+fntemp+"\" \""+fn+"\"" ;
cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; cp "+fnq ;
IdxTy rt= Execute(0,0, cmd); 
IdxTy lc2=Cmd(cmd,"lightcleanchange",0);
cmd="cat \""+fn+"\" | "+cmd + "| tee \"" + fntemp+"\" ; diff "+fnq+" ;  cp "+fnq ;
Blob data,out,err;
IdxTy rt2=m_hand.fileio(out,err,data,cmd,3);
MM_ERR(" light cleaning change err  "<<MMPR2(cmd,StrTy(err)))
MM_ERR(" light cleaning change out  "<<MMPR(StrTy(out)))
rc.rendered(StrTy(out));

rc=rt;
#endif

return rc;
}



AccessTy SafeFinalAssembly(Blob & d, const Blob & s, const assem_opt & ao)
const {
AccessTy rc=0;
d=StrTy(s); 

// dos2unix, insert url, insert comments 

return rc;
}

AccessTy AddToFile(Blob & b, const StrTy & fn, const IdxTy flags)
{
AccessTy rc=0;
b.save(fn,1);
return rc;
}
IdxTy Cmd(StrTy & cmd, const StrTy & xform, const IdxTy flags) const 
{
auto ii=m_map.find(xform);
if (ii==m_map.end()) 
{
MM_ERR(" no xform for "<<MMPR3(xform,cmd,flags))
return ~0;
} 
 cmd=(*ii).second;
return 0;
} // Cmd
IdxTy Cmd2(StrTy & cmd, const StrTy & xform, const IdxTy flags) const 
{
auto ii=m_map_2.find(xform);
if (ii!=m_map_2.end()) { cmd+=(*ii).second; } 
return 0;
} // Cmd2
IdxTy Code(StrTy & cmd, const StrTy & xform, const IdxTy flags=0) const 
{
auto ii=m_code.find(xform);
if (ii!=m_code.end()) { cmd+=(*ii).second; } 
return 0;
} // Cmd2



AccessTy FileToBlob(Blob & d, const StrTy & fn, const StrTy & xform,  const IdxTy flags) const 
{
AccessTy rc=0;
StrTy cmd;
if (Cmd(cmd,xform,flags)) return rc;
//cmd="cat \""+fn+"\" | "+cmd;
cmd=cmd+" "+fn;
Cmd2(cmd,xform,flags);
Blob data,out,err;
IdxTy c=m_hand.fileio(out,err,data,cmd);
MM_ERR(MMPR4(c,StrTy(err),xform,fn))
return rc; 
}
AccessTy BlobToBlob(Blob & d, const Blob & s, const StrTy & xform,  const IdxTy flags) const 
{
AccessTy rc=0;
StrTy cmd;
if (Cmd(cmd,xform,flags)) return rc;
//cmd="cat \""+fn+"\" | "+cmd;
//cmd=cmd+" - ";
Cmd2(cmd,xform,flags);
Blob out,err;
IdxTy c=m_hand.fileio(out,err,s,cmd);
MM_ERR(MMPR3(c,StrTy(err),xform))
return rc; 
}



AccessTy FileToFile(const StrTy  & dfn, const StrTy & fn, const StrTy & xform,  const IdxTy flags) const 
{
AccessTy rc=0;
StrTy cmd;
if (Cmd(cmd,xform,flags)) return rc;
// hopefully ok with the fileio impl else write a new one 
//cmd="cat \""+fn+"\" | "+cmd+" > \""+dfn+"\"";;
//cmd=cmd+" \""+fn+"\"  > \""+dfn+"\"";
cmd=cmd+" \""+fn+"\"";
Cmd2(cmd,xform,flags);
cmd=cmd+"  > \""+dfn+"\"";
//cmd=cmd+" "+fn;
Blob data,b,err;
IdxTy c=m_hand.fileio(data,err,b,cmd);
MM_ERR(MMPR4(c,StrTy(err),xform,fn))
return rc; 
}
IdxTy Execute(Ragged * rout, Ragged * rerr, const StrTy & cmd) const
{
Blob data,out,err;
//IdxTy c=m_hand.lineio(out,err,data,cmd);
// remove trailing nl but also enclose whatever is there in parens,
IdxTy c=m_hand.fileio(out,err,data,cmd,3);
//MM_ERR(MMPR4(c,cmd,StrTy(err),StrTy(out)))
if (rout) LoadRagged(*rout,StrTy(out));
if (rerr) LoadRagged(*rerr,StrTy(err));
return c; 
} // Execute
AccessTy Parse(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const 
{
AccessTy rc=0;
BibMap&  bm=rc.bm();
bm.parse(ofn);
const IdxTy sz=bm.size();
if( sz!=1)
{ MM_ERR(" should have exactly one entry "<<MMPR2(ofn,sz))

}
//if( sz>4 ) { rc=~0;   MM_ERR(" something very wrong with s"<<MMPR2(ofn,sz)) return rc; } 

if (!sz) { rc=1; return rc; }
rc.be()=*(bm.begin()); //  ((*(bm.begin())).second)[0];
const BibEntry & be=rc.be(); //  ((*(bm.begin())).second)[0];
const StrTy & nm=be.name();
const StrTy & ty=be.type();
const IdxTy sze=be.size();
const IdxTy errors=be.errors();
rc.errors(errors);
MM_ERR(MMPR4(nm,ty,sze,errors))
IdxTy rci=0;
if (nm.length()==0)  { rci|=2; }
if (ty.length()==0)  { rci|=4; }
if (sze<3) { rci|=8;}
if (errors) { rci|=16;}
rc=rci;
return rc;

}

AccessTy Validate(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const 
{
AccessTy rc=0;
if (ofn=="")
{
MM_MSG(" the ofn into Validate is blank something is wrong...")
MM_ERR(" the ofn into Validate is blank something is wrong...")
rc=-1;
return rc;
}
StrTy vcode;
Code(vcode, "validate");
Ragged r,rout,rerr;
LoadRagged(r,vcode);
StrTy p=ReplaceRagged(r,ofn, "mybib");
//MM_ERR(MMPR(p))
IdxTy rcl=Execute(&rout,&rerr,p);
//MM_ERR(MMPR(rout.dump(18*0,"")));
//MM_ERR(MMPR(rerr.dump(18*0,"")));
 StrTy fnerr;
Cmd(fnerr,"checkout",0);
rout.save(fnerr,16+2,"");
rerr.append(fnerr,16+2,"");

{
Blob data,out,err;
// TODO move to maps for later modification... 
const StrTy  cmd="cat "+fnerr+" | grep -a  -B 3  -i  \"^l\\.[0-9][0-9]* \\|output\\|error\\|warning\" | grep -a  -v \"^[0-9]*(\" |grep -a  -v \"^rm: \" |  sed -e 's/  */ /g' | grep -v \"^[[:space:]]$\" | grep -v \"/usr/share/texlive/\\|/var/lib/texmf\\|^--$\\|aux[)]$\\|texmf-dist\\|have changed. Rerun to get cross-references\" ";
//MM_ERR(MMPR(cmd))
IdxTy c=m_hand.fileio(out,err,data,cmd,3);
//MM_ERR(MMPR4(c,StrTy(err),StrTy(out),StrTy(data)))
if (StrTy(out).length()!=0) MM_MSG(CRLF<<StrTy(out))
rc.latex_output()=StrTy(out);
rc.validation_fn(fnerr);
}
//MM_ERR(MMPR2(rcl,fnerr))
Blob rendered;
Cmd(fnerr,"rendered",0);
rendered.load(fnerr);
//MM_ERR(MMPR(StrTy(rendered)))
rc.rendered(StrTy(rendered));

return rc; 
}


#if 0 
AccessTy Validate_xxx(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const 
{
AccessTy rc=0;
StrTy vcode;
// remove the old test files 
Code(vcode, "rmtestdoc");
Ragged r;
LoadRagged(r,vcode);
StrTy p=ReplaceRagged(r,ofn, "mybib");
MM_ERR(MMPR(p))
// write_file(fn,18);
Ragged rout, rerr;
Execute(&rout,&rerr,p);
r.clear();
vcode=StrTy();
Code(vcode, "testdoc");
LoadRagged(r,vcode);
MM_ERR(MMPR(r.dump()))
const StrTy oftex=ofn+".tex";
const StrTy ofbib=ofn+".bib";
_b.save(ofbib);
 p=ReplaceRagged(r,ofn, "mybib");
{std::ofstream ofs(oftex); ofs<<p; }

Cmd(p,"latex",0);
p=p+" "+ofn;
IdxTy rcl=Execute(&rout,&rerr,p);
Cmd(p,"bibtex",0);
p=p+" "+ofn;
IdxTy rcf=Execute(&rout,&rerr,p);
MM_ERR(MMPR(rout.dump(18*0,"")));
MM_ERR(MMPR(rerr.dump(18*0,"")));
const StrTy fnerr=m_map["checkout"];
rout.save(fnerr);
rerr.append(fnerr);
MM_ERR(MMPR3(rcl,rcf,fnerr))
return rc; 
}
#endif

void LoadRagged(Ragged & r, const StrTy & s) const
{ // r.set_split(8,' '); 
r.splitter(8);
Ss ss; ss<<s; r.load(ss); }
StrTy ReplaceRagged(Ragged & r, const StrTy & sout, const StrTy & sin) const
{
MM_LOOP(ii,r) MM_LOOP(jj,(*ii)) 
{ if ((*jj)==sin) *jj=sout; } 
Ss ss;
ss<<r.dump(18,"");
return ss.str();
} // ReplaceRagged

//IdxTy c=m_hand.fileio(dest,err,data,cmd);

void Init(const Ragged & r)
{
Init();
const IdxTy sz=r.size();
for(IdxTy i=0; i<sz; ++i)
{


} // i 

} // Init 

/*
#http://tex.stackexchange.com/questions/173621/how-to-validate-check-a-biblatex-bib-file
testbibtex()
{
cat - << ----
\documentclass{article}
\begin{document}
\nocite{*}
\bibliographystyle{plain}
\bibliography{mybib}
\end{document}
----

}

*/

void InitCode() { 
{Ss ss;
ss<<"\\documentclass{article}"<<eol;
ss<<"\\usepackage{hyperref}"<<eol;
ss<<"\\usepackage[paperwidth=200in, paperheight=200in]{geometry}"<<eol;
ss<<"\\begin{document}"<<eol;
ss<<"\\nocite{*}"<<eol;
ss<<"\\bibliographystyle{plainurl}"<<eol;
ss<<"\\bibliography{mybib}"<<eol;
ss<<"\\end{document}"<<eol;
m_code["testdoc"]=ss.str(); 
}
{
Ss ss;
ss<<"rm mybib.xxx"<<eol;
ss<<"rm mybib.pdf"<<eol;
ss<<"rm mybib.aux"<<eol;
ss<<"rm mybib.blg"<<eol;
ss<<"rm mybib.bbl"<<eol;
ss<<"rm mybib.log"<<eol;
ss<<"rm mybib.bib"<<eol;

m_code["rmtestdoc"]=ss.str(); 
} 
{
Ss ss;
ss<<"rm $x.xxx"<<eol;
ss<<"rm $x.pdf"<<eol;
ss<<"rm $x.aux"<<eol;
ss<<"rm $x.blg"<<eol;
ss<<"rm $x.bbl"<<eol;
ss<<"rm $x.log"<<eol;
//ss<<"rm $x.bib"<<eol;

m_code["xrmtestdoc"]=ss.str(); 
} 




{
Ss ss;
ss<<"bibtemp()"<<eol;
ss<<"{"<<eol;
ss<<"cat - << ---xxx---"<<eol;
ss<<m_code["testdoc"];
ss<<"---xxx---"<<eol;
ss<<"}"<<eol;
// this seems to be a bug in bibtex as the suffix is ambiguous
ss<<"x=`basename mybib | sed -e 's/\\.bib$//'`"<<eol;
ss<<"rm  /tmo/$x.pdf /tmp/$x.bib"<<eol;
ss<<"cp mybib /tmp/$x.bib"<<eol;
ss<<"cd /tmp"<<eol;
ss<<"pwd"<<eol;
ss<<m_code["xrmtestdoc"]; // <<eol;
ss<<"bibtemp |tee  $x.tex "<<eol;
// sed -n 'w /tmp/out'
//ss<<m_map["latex"]<<" "<<"$x.tex |tee  "<<m_map["checkout"]<<eol;
// https://superuser.com/questions/954849/standard-utility-to-write-stdin-to-a-file-with-no-stdout
// dump first run as it always has errors and picked up in second 
ss<<m_map["latex"]<<" "<<"$x.tex |sed -n 'w "<<m_map["checkout"]<<"'"<<eol;
ss<<m_map["bibtex"]<<" "<<"$x  | tee -a "<<m_map["checkout"]<<eol;
ss<<m_map["latex"]<<" "<<"$x.tex |tee -a  "<<m_map["checkout"]<<eol;
ss<<"echo latex_rc $?"<<eol;
ss<<"rm \""+m_map["rendered"]+"\""<<eol;
ss<<"rm \""+m_map["rtemp"]+"\""<<eol;
//ss<<"pdftotext \"$x.pdf\" \""+m_map["rendered"]+"\""<<eol;
ss<<"pdftotext \"$x.pdf\" \""+m_map["rtemp"]+"\""<<eol;
ss<<"cat \""+m_map["rtemp"] +"\" | strings -n 1 | tee \""+ m_map["rendered"]+"\""<<eol;
m_code["validate"]=ss.str(); 
}

} // Initcode


void Init() { 
 eol="\n";
m_map["latex"]="pdflatex -file-line-error -halt-on-error "; 
m_map["bibtex"]= "bibtex "; 
m_map["checkout"]="checkbib_test_output.xxx";
m_map["rendered"]="/tmp/checkbib_rendered.txt";
m_map["rtemp"]="/tmp/.checkbib_rendered.txt";
m_map["lightclean"]="dos2unix | grep -a  -v \"^[[:space:]]*$\"";
m_map["lightcleanchange"]="sed -e 's/\\([^\\\\]\\)&/\\1\\\\\\&/g'";
InitCode();
} // Init
// MEMBERS


mutable Hand m_hand;
Tdb m_config;
XformMap m_map;
XformMap m_map_2;

XformMap m_code;
StrTy eol;
//StrTy m_hdoi,m_tdoi,m_mdoi,m_udoi;
//StrTy m_doised2;

}; // mjm_assemble_putative_bibtex

//////////////////////////////////////////////

template <class Tr>
class mjm_assemble_putative_bibtex_map : public std::map<typename Tr::StrTy, mjm_assemble_putative_bibtex< Tr > >  
{
 typedef mjm_assemble_putative_bibtex_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_assemble_putative_bibtex< Tr> >   Super;
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
mjm_assemble_putative_bibtex_map() {}
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

}; // mjm_assemble_putative_bibtex_map




////////////////////////////////////////////
#ifdef  TEST_MJM_ASSEMBLE_PUTATIVE_BIBTEX
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
typedef tester_< mjm_assemble_putative_bibtex <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_ASSEMBLE_PUTATIVE_BIBTEX "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_assemble_putative_bibtex<Tr>  Myt;
//Myt x(argc,args);
Myt x;
typedef Myt::blob_type Blob;
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
//AccessTy Validate(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) const 
else if (cmd=="val") { Blob b; auto rc=x.validate(b,cip.p1,cip.p2,atoi(cip.wif(3).c_str())); MM_ERR(rc.dump()) }
else if (cmd=="parse") { Blob b; auto rc=x.parse(b,cip.p1,cip.p2,atoi(cip.wif(3).c_str())); MM_ERR(MMPR2(rc.dump(),rc.be().dump_errors(1))) }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_ASSEMBLE_PUTATIVE_BIBTEX_H__ 
