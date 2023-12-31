#ifndef MJM_DOI_SCRAPE_H__
#define MJM_DOI_SCRAPE_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_pawnoff.h"
#include "mjm_collections.h"
#include "mjm_wovdb.h"
#include "mjm_read_buffer.h"
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


// Mon Jan  4 13:43:47 EST 2021
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_doi_scrape   
// g++  -Wall -std=gnu++11 -DTEST_MJM_DOI_SCRAPE -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_doi_scrape.h  -lpthread -lreadline


/*

contrary to their comment, the html is usually dirty although it can
be easier, 
// https://unix.stackexchange.com/questions/383150/sed-to-print-only-first-pattern-match-of-the-line

10:29:07 UTC
marchywka@happy:/home/documents/cpp/proj/toobib$ cat /tmp/fileOXmTYg | grep 10 | perl -pe 's/.*?10/10/' | sed -e 's/ /\n/g' | grep "10\."

*/



template <class Tr>
class mjm_doi_scrape 
{
 typedef mjm_doi_scrape Myt;
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
typedef mjm_strings StrUtil; 
class _Scraping
{

public:
_Scraping() { m_flags=0; } 
_Scraping(const StrTy & s, const StrTy & m, const IdxTy f)
:m_scrape(s),m_method(m),m_flags(f) {}
_Scraping(const StrTy & s)
:m_scrape(s),m_method(),m_flags(0) {}
_Scraping & operator=(const StrTy & s) { m_scrape=s; return *this; } 
operator StrTy() const { return m_scrape; }
// wtf 
const char * c_str() const { return m_scrape.c_str(); }
const StrTy method() const { return m_method; } 

StrTy m_scrape;
StrTy m_method;
IdxTy m_flags;


}; // _Scraping

typedef _Scraping Scraping;
typedef std::vector<Scraping> Scrapings;

class _ResultType
{
typedef _ResultType Myt;
public:
_ResultType() {Init(0); }
_ResultType(const IdxTy n) {Init(n); }
Myt & operator=(const IdxTy n ) { Init(n);  return *this; }
operator IdxTy() const { return 0; } 
private:

void Init(const IdxTy n) {}


}; // _Result_Type

typedef _ResultType ResultType;
typedef ResultType AccessTy;


/*
Scrape a file for anything that looks like a doi.
File could be text, html, pdf binary, rendered html,
pdf2text output, exif output, various structured formats.
Do not use various deefined methods or handlers yet. 

*/


public:
typedef std::vector<Scraping> scraped_vector;

typedef AccessTy result_code;
mjm_doi_scrape() {Init();}
mjm_doi_scrape(const Ragged & r) {Init(r); }
StrTy dump(const IdxTy flags=0) { return Dump(flags); }

template <class Td> AccessTy scrape_url(Td & v, const StrTy & url, const IdxTy flags)
{ return ScrapeURL(v,url,flags); }

template <class Td> AccessTy scrape_file(Td & v, const StrTy & fn, const IdxTy flags)
{ return ScrapeFile(v,fn,flags); }

template <class Td> AccessTy scrape_blob(Td & v, const Blob & b, const IdxTy flags)
{ return ScrapeBlob(v,StrTy(),b,flags); }

template <class Td> AccessTy scrape_bin_thing(Td & v, const StrTy & pfx, const Blob & b,const StrTy & type,  const IdxTy flags)
{ return ScrapeBinThing(v,pfx,b,type,flags); } 
template <class Td> AccessTy scrape_bin_file(Td & v, const StrTy & fn,const StrTy & type,  const IdxTy flags)
{ return ScrapeBinThing(v,CatPipe(fn),Blob(),type,flags); } 


template <class Td> AccessTy scrape_bin_blob(Td & v,  const Blob & b,const StrTy & type,  const IdxTy flags)
{ return ScrapeBinThing(v,"",b,type,flags); } 


private:
const StrTy CatPipe(const StrTy & fn) { return "cat \""+fn+"\" |"; }
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

//SortScrapes(td,pfx,b,flags);
template <class Td> AccessTy SortScrapes(Td & v, const StrTy & pfx, const Blob & b,  const IdxTy flags)
{
AccessTy rc=0;
Td vout;
StrTy abs=StrTy(b);
const IdxTy vs=v.size();
std::vector<int> vidx(vs);
IdxTy idx=0;
MM_LOOP(ii,v)
{
IdxTy loc=StrUtil::indexOf(abs.c_str(),(*ii).c_str());
MM_ERR(MMPR3(loc,idx,(*ii).c_str()))
vidx[idx]=loc;
++idx;
}
MM_ERR(" sort the ")
// cant believe c__ does not have permutation order sort 
// https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
std::vector<int> index(vidx.size());
for (int i = 0 ; i < index.size() ; i++) { index[i] = i; }
std::sort(index.begin(), index.end(),
    [&](const int& a, const int& b) {
MM_ERR(MMPR2(a,b))
        return (vidx[a] < vidx[b]);
    }
);
idx=0;
MM_ERR(" done with sort")
MM_LOOP(ii,index)
{
MM_ERR(" "<<MMPR2(idx,(*ii)))
vout.push_back(v[(*ii)]); 
MM_ERR(MMPR3(idx,StrTy(vout.back()),(*ii)))
++idx;
}
MM_ERR(" done with pushing  sort")

MM_ERR(" sorting "<<MMPR2(v.size(),abs.length()))

v=vout;
MM_ERR(" done with voit assign"  )

return rc;
} // SortScrapes


template <class Td> AccessTy ScrapeBlob(Td & v, const StrTy & pfx, const Blob & b, const StrTy & method, const IdxTy flags)
{
AccessTy rc=0;
const bool one_per_customer=Bit(flags,0);
const IdxTy sz0=v.size();
{
Blob out,err;
StrTy cmd=pfx+method;
IdxTy c=m_hand.fileio(out,err,b,cmd,3);
MM_ERR(" scraping "<<MMPR4(c,cmd,StrTy(out),StrTy(err))<<MMPR(StrTy(b)))
out.line_vector(v);
}
const IdxTy sz1=v.size();
const IdxTy dx=sz1-sz0;
if (one_per_customer) if (dx>1) v.resize(sz0+1);
MM_ERR(MMPR4(one_per_customer,dx,sz0,v.size())) 
return rc;
} // ScrapeBlob 


template <class Td> AccessTy ScrapeBlob(Td & v, const StrTy & pfx, const Blob & b, const IdxTy flags)
{
std::vector<StrTy> methods;
methods.push_back(m_ddoi);
methods.push_back(m_mdoi);
methods.push_back(m_tdoi);
methods.push_back(m_hdoi);
methods.push_back(m_arxiv);
return  ScrapeList(v,  pfx, b, methods,   "",   flags);

/*const bool first_result_only=Bit(flags,1);
AccessTy rc=0;
const IdxTy sz0=v.size();
AccessTy rc1=ScrapeBlob(v,pfx,b,m_mdoi,flags);
if (first_result_only) if (v.size()) return rc;
AccessTy rc2=ScrapeBlob(v,pfx,b,m_tdoi,flags);
if (first_result_only) if (v.size()) return rc;
AccessTy rc3=ScrapeBlob(v,pfx,b,m_hdoi,flags);
return rc; 
*/

} // ScrapeBlob 

template <class Td,class Tl > AccessTy ScrapeList(Td & v, const StrTy & pfx, const Blob & b,const Tl & methods, const StrTy & type,  const IdxTy flags)
{
const bool first_result_only=Bit(flags,1);
AccessTy rc=0;
MM_LOOP(ii,methods)
{
rc=ScrapeBlob(v,pfx,b,(*ii),flags);
MM_ERR(" in scrape list "<<MMPR2(v.size(),(*ii)))
if (first_result_only) if (v.size()) return rc;
} // ii 
// now sort by order of occurence, although this relies on html order... 
SortScrapes(v,pfx,b,flags);
return rc;
} // ScrapeList

template <class Td> AccessTy ScrapeBinThing(Td & v, const StrTy & pfx, const Blob & b,const StrTy & type,  const IdxTy flags)
{
const bool first_result_only=Bit(flags,1);
MM_ERR(" ScrapeBinThing")
// this was picking up the wrong one first due to "^doi:" being later in the order although 
// it eventually picked up the right one too, 
// https://cob.silverchair-cdn.com/cob/content_public/journal/jcs/133/8/10.1242_jcs.241976/3/jcs241976.pdf?Expires=1627343651&Signature=rFMJjpRroH8tV86-MCYI6yiCq7V2U8Oxl0u8TuhQGTRPckPhoOPWAEjjqud8ZtoGzCLvx0SiO56R3RhS7xNgFqdEi8qV43VgL-KwwozVe8W1z8JXadZIRWfr3DHa3uITkGApvMUzM77YVnxEJbzWEs19fNiJQBuHLot~FEcVnaqWV~gl2nHcUDlgBiex48ZZ3EHbtiCg7xeoZPn-E5ITka0eMQaOWqO1IYB59EKNsX8BNbwFCngpUBjD2XRNGZ7VMlyziQE8KQKlJ4jTHXQKmYYXivZYPuA3ZVvubFb8co5y02yf2k-aZhKiyaJfRtqgA4X3NLpNUVfeKK7hRdwqkg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

// the URI filter below is due to this, 
 // https://www.jneurosci.org/content/jneuro/37/23/5770.full.pdf
// although the pdftotext should probably be doi scanned before the
// bin the doi is below the first page and the citations apparently
// create doi links FCK 

//AccessTy rc=0;
//const StrTy sfx=" | grep -v \"\\\\\\\\\" |  sed -e 's/[.;)]$//' | sed -e 's/[)<>].*//' | grep \"10\" ";
//const StrTy sfx=" | grep -v \"\\\\\\\\\" |  sed -e 's/[.;)]$//' | sed -e 's/[)<>].*//' |"+m_isolate_10;
const StrTy sfx=" | grep -v \"\\\\\\\\\" |sed -e 's/^doi=//i' |  sed -e 's/[.;)]$//' | sed -e 's/[)<>].*//' |"+m_isolate_10;
//const StrTy c1_doix="`grep -a -i  \"doi\\|10\\..*/\" \""+fn+"\" | sed -e 's/.this version posted.*//' | grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;
//const StrTy c0_doix="strings -n 1 |grep -a -i  \"^doi[ =:]*10\" "+sfx;
//const StrTy c0_doix="strings -n 1 |sed -e 's/ /\\n/g' | grep -a -i  \"^doi[ =:]*10\" "+sfx;
const StrTy c0_doix="strings -n 1 |sed -e 's/ /\\n/g' | grep -a -i  \"^doi[ =:]*10\\|doi.org/\" "+sfx;
const StrTy c1_doix="strings -n 1 | grep -v \"/URI \" |  grep -a -i  \"doi\\|10\\..*/\"  | sed -e 's/.this version posted.*//' | grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;
//const StrTy c2_doix="grep -a -i  \"first published \" \""+fn+"\" |grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;
const StrTy c2_doix="strings -n 1 | grep -v \"/URI \" | grep -a -i  \"first published \" |grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;
//const StrTy c3_doix="strings -n 1 |  sed -e 's/  */\\\\n/g' | grep 10| head -n 1 | sed -e 's/http.*org.//' "+sfx;
const StrTy c3_doix="strings -n 1 |grep -v \"/URI \" |   sed -e 's/  */\\\\n/g' | grep 10| head -n 1 | sed -e 's/http.*org.//' "+sfx;
const StrTy c4_doix="strings -n 1 |  sed -e 's/  */\\\\n/g' | grep -i  \"DOI:\" | head -n 1 | sed -e 's/.*  *//' "+sfx;
const StrTy c5_doix="grep -i DOI | awk '{if (NF==2) print $2;}' | grep \"^10\" "+sfx;

const IdxTy sz0=v.size();
/*
AccessTy rc1=ScrapeBlob(v,pfx,b,c1_doix,flags);
if (first_result_only) if (v.size()) return rc;
AccessTy rc2=ScrapeBlob(v,pfx,b,c2_doix,flags);
if (first_result_only) if (v.size()) return rc;
AccessTy rc3=ScrapeBlob(v,pfx,b,c3_doix,flags);
*/
std::vector<StrTy> methods;
methods.push_back(c0_doix);
methods.push_back(c1_doix);
methods.push_back(c2_doix);
methods.push_back(c3_doix);
methods.push_back(c4_doix);
methods.push_back(c5_doix);
return  ScrapeList( v,  pfx, b, methods,  type,   flags);

} // ScrapeBinThing
  

//blob.line_vector(dest);
//IdxTy c=m_hand.fileio(dest,err,data,cmd);

template <class Td> AccessTy ScrapeFile(Td & v, const StrTy & fn, const IdxTy flags)
{
const StrTy pfx="cat \""+fn+"\" |";
MM_ERR(" scraping file "<<MMPR2(fn,pfx))
// the sorter now needs the blod loaded
Blob b;
b.load(fn);

//return ScrapeBlob(v,pfx,Blob(),flags);
return ScrapeBlob(v,pfx,b,flags);
// turn lines into strings 
//blob.line_vector(dest);
} // ScrapeFile


template <class Td> AccessTy ScrapeURL(Td & v, const StrTy & url, const IdxTy flags)
{
AccessTy rc=0;
const IdxTy szi=v.size();
Blob in,out,err;
StrTy cmd="echo \""+url+"\" | " + m_udoi;
IdxTy c=m_hand.fileio(out,err,in,cmd);
out.line_vector(v);
const IdxTy szf=v.size();
// turn lines into strings 
//blob.line_vector(dest);
return rc;
} // ScrapeUrl


//IdxTy c=m_hand.fileio(dest,err,data,cmd);

void Init(const Ragged & r)
{
Init();
const IdxTy sz=r.size();
for(IdxTy i=0; i<sz; ++i)
{


} // i 

} // Init 

#if 0 
urlfck="$1" 
 doixxx=`echo $urlfck | sed -e 's/.*\(\/10\.[0-9][0-9]*\/[^;&?#]*\).*/\1/' | grep -v "^http://\|^https://" `
# remove common suffixes not like part of doi 
doixxx=`echo $doixxx | sed -e 's/\/full$//'| sed -e 's/^\///' `
echo $doixxx
}


#canhandle=`cat "$fn" | grep -a -i "doi=" | grep -v "data-journal" |  sed -e 's/.*doi=/doi=/ig' |  sed -e 's/["]//g'| sed -e 's/[^/0-9a-zA-Z.=%].*/\n/g'  | grep "doi=" | grep -v "data-journal" | sed -e 's/doi=//'  | head -n 1`


doisfrommeta()
{
sed -e 's/<meta name=/\n<meta name=/g' | sed -e 's/[}>]/&\n/g' | grep -a  citation_doi | sed -e 's/.*content="//' | sed -e 's/".*//' |doised2 | sed -e 's/ .*//g'| grep 10 | head -n 1

}
#endif


#if 0


handledoi()
doised2()
{
sed -e 's/[^-()/a-zA-Z_0-9.].*//'
}

handledoi()
{
echo "  handledoi " | smsg
fn="$1"
uin="$3"
fn2="$temp5"
canhandle=`echo $uin | sed -e 's/.*\(\/10\.[0-9][0-9]*\/[^;&?#]*\).*/\1/' | grep -v "^http://\|^https://" `
if [ "$canhandle" == "" ]
then
canhandle=`echo $uin| sed -e 's/.*doi=//' | sed -e 's/[&=;#].*//' | sed -e 's/^http.*//' `
fi


# remove common suffixes not like part of doi 
canhandle=`echo $canhandle | sed -e 's/\/full$//'| sed -e 's/^\///' `
if [ "$canhandle" != "" ]
then
echo FOO doi found in url $canhandle | smsg
fi

if [ "$canhandle" == "" ]
then
#echo getting was fn again someone is trashing it... | smsg
#$WGET --no-check-certificate --user-agent Firefox -q -O "$fn" "$uin"
#cat "$fn" | grep "doi=" | smsg
#canhandle=`cat "$fn" | grep "doi" `
#canhandle=`cat "$fn" | sed -e 's/<meta name=/\n<meta name=/g' | grep -a  citation_doi | sed -e 's/.*content="//' | sed -e 's/".*//' | sed -e 's/[^/a-zA-Z0-9.].*//' `

canhandle=`cat "$fn" | sed -e 's/<meta name=/\n<meta name=/g' | grep -a  citation_doi | sed -e 's/.*content="//' | sed -e 's/".*//' | sed -e 's/[^-()/a-zA-Z0-9.].*//'| sed -e 's/ .*//g'| grep 10 | head -n 1 `

# original 
#canhandle=`cat "$fn" | sed -e 's/<meta name=/\n<meta name=/g' | sed -e 's/[}>]/&\n/g' | grep -a  citation_doi | sed -e 's/.*content="//' | sed -e 's/".*//' | sed -e 's/[^-()/a-zA-Z0-9.].*//'| sed -e 's/ .*//g'| grep 10 | head -n 1 `
canhandle=`cat "$fn" | sed -e 's/<meta name=/\n<meta name=/g' | sed -e 's/[}>]/&\n/g' | grep -a  citation_doi | sed -e 's/.*content="//' | sed -e 's/".*//' |doised2 | sed -e 's/ .*//g'| grep 10 | head -n 1 `


#canhandle=`cat "$fn" | sed -e 's/<meta name=/\n<meta name=/g' | grep -a  citation_doi | sed -e 's/.*content="//' | grep 10 | head -n 1 `
#canhdnalge=`echo $canhandle |awk '{print $1}' ` 
canhandle=`echo $canhandle |awk '{print $1}' `
if [ "$canhandle" == "" ]
then
echo try doi extract onlyu doi | smsg
#canhandle=`cat "$fn" | sed -e 's/doi=./\ndoi /gi' | sed -e 's/["<>]/\n/g' | grep doi | awk '{print $2}' ` 
canhandle=`cat "$fn" | grep -i "doi=" | sed -e 's/doi=/\ndoi /gi' | sed -e 's/[<>]/\n/g' | sed -e 's/doi [^0-9]*10/doi 10/g' |  sed -e 's/["<>]/\n/g' | grep doi | awk '{print $2}' `

#endif
// doisfromtext
// grep -a -i "doi[=:]"  | sed -e 's/.*doi[:=]/\ndoi=/gi' | sed -e 's/["]//g'| sed -e "s/'//g" |  sed -e 's/[^/()0-9a-zA-Z.=% -].*/\n/g'  |  grep "^doi=" | sed -e 's/doi=//i'  |  awk '{print $1}'  | head -n 1


void Init() { 
//m_isolate_10=" grep 10 | perl -pe 's/.*?10/10/' | sed -e 's/ /\\n/g' | grep \"10\\.\"";
m_isolate_10=" grep 10 | perl -pe 's/.*?10\\./10./' | sed -e 's/ /\\n/g' | grep \"10\\.\"";
// doised2() { sed -e 's/[^-()/a-zA-Z_0-9.].*//' }
//m_doised2=" doised2() { sed -e 's/[^-()/a-zA-Z_0-9.].*//' }";
m_doised2="  sed -e 's/[^-()/a-zA-Z_0-9.].*//' ";
// doifromutl
//m_udoi="strings -n 1 |  sed -e 's/.*\\(\\/10\\.[0-9][0-9]*\\/[^;&?#]*\\).*/\\1/' | grep -v \"^http://\\|^https://\"  | sed -e 's/\\/full$//'| sed -e 's/^\\///'"; 
m_udoi="strings -n 1 |  sed -e 's/.*\\(\\/10\\.[0-9][0-9]*\\/[^;&?#]*\\).*/\\1/' | grep -v \"^http://\\|^https://\"  | sed -e 's/[\\/.]pdf[?].*//' | sed -e 's/[\\/.]pdf$//' |  sed -e 's/\\/full$//'| sed -e 's/^\\///'"; 
// doisfrommeta 
//m_mdoi=def_doised2+ "; sed -e 's/<meta name=/\\n<meta name=/g' | sed -e 's/[}>]/&\\n/g' | grep -a  citation_doi | sed -e 's/.*content=\"//' | sed -e 's/\".*//' |doised2 | sed -e 's/ .*//g'| grep 10 | head -n 1";
//m_mdoi= "strings -n 1 |  sed -e 's/<meta name=/\\n<meta name=/g' | sed -e 's/[}>]/&\\n/g' | grep -a  citation_doi | sed -e 's/.*content=\"//' | sed -e 's/\".*//' |"+m_doised2+" | sed -e 's/ .*//g'| grep 10 | head -n 1";
m_ddoi= "strings -n 1 |  sed -e 's/http/\\nhttp/ig' | sed -e 's/[)}>]/\\n/g' | grep -a  doi\\.org | sed -e 's/.*doi.org.//' | sed -e 's/\".*//' |"+m_doised2+" | sed -e 's/ .*//g'| "+m_isolate_10+ "| head -n 1";
m_mdoi= "strings -n 1 |  sed -e 's/<meta name=/\\n<meta name=/g' | sed -e 's/[}>]/&\\n/g' | grep -a  citation_doi | sed -e 's/.*content=\"//' | sed -e 's/\".*//' |"+m_doised2+" | sed -e 's/ .*//g'| "+m_isolate_10+ "| head -n 1";
// doisfromtext 
//m_tdoi="grep -a -i \"doi[=:]\"  | sed -e 's/.*doi[:=]/\\ndoi=/gi' | sed -e 's/[\"]//g'| sed -e \"s/'//g\" |  sed -e 's/[^/()0-9a-zA-Z.=% -].*/\\n/g'  |  grep \"^doi=\" | sed -e 's/doi=//i'  |  awk '{print $1}'  | head -n 1";
// this should not be done if there is a citation_blah but just in case
// remove it 
m_tdoi="strings -n 1 | grep -v -i \"citation_reference\\|ref-cit\" |   grep -a -i \"doi[=:]\"  | sed -e 's/.*doi[:=]/\\ndoi=/gi' | sed -e 's/[\"]//g'| sed -e \"s/'//g\" |  sed -e 's/[^/()0-9a-zA-Z.=% -].*/\\n/g'  |  grep \"^doi=\" | sed -e 's/doi=//i'  |  awk '{print $1}'| "+m_isolate_10;
// handledoi thing 
 //m_hdoi="strings -n 1 |  grep -v -i \"citation_reference\\|ref-cit\" | grep -i \"doi=\" | sed -e 's/doi=/\\ndoi /gi' | sed -e 's/[<>]/\\n/g' | sed -e 's/doi [^0-9]*10/doi 10/g' |  sed -e 's/[\"<>]/\\n/g' | grep doi | awk '{print $2}'"; 
 m_hdoi="strings -n 1 |  grep -a -v -i \"citation_reference\\|ref-cit\" | grep -a -i \"doi=\" | sed -e 's/doi=/\\ndoi /gi' | sed -e 's/[<>]/\\n/g' | sed -e 's/doi [^0-9]*10/doi 10/g' |  sed -e 's/[\"<>]/\\n/g' | grep doi | awk '{print $2}'| "+m_isolate_10; 
 m_arxiv="strings -n 1 | grep \"arXiv:\" | awk '{print $1}' | head -n 1  "; 

 // canhandle=`cat "$fn" | grep -v "data-journal" |doisfromtext`




} // Init
// MEMBERS


Hand m_hand;
Tdb m_config;
StrTy m_hdoi,m_tdoi,m_mdoi,m_ddoi,m_udoi,m_arxiv;
StrTy m_isolate_10;
StrTy m_doised2;

}; // mjm_doi_scrape

//////////////////////////////////////////////

template <class Tr>
class mjm_doi_scrape_map : public std::map<typename Tr::StrTy, mjm_doi_scrape< Tr > >  
{
 typedef mjm_doi_scrape_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_doi_scrape< Tr> >   Super;
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
mjm_doi_scrape_map() {}
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

}; // mjm_doi_scrape_map




////////////////////////////////////////////
#ifdef  TEST_MJM_DOI_SCRAPE
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
typedef tester_< mjm_doi_scrape <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_DOI_SCRAPE "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_doi_scrape<Tr>  Myt;
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

#endif // MJM_DOI_SCRAPE_H__ 
