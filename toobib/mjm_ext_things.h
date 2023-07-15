
#ifndef MJM_ext_things_H__
#define MJM_ext_things_H__
 
#include "mjm_globals.h"

//#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
//#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_read_buffer.h"
#include "mjm_char_star_map.h"

//#include "mjm_canned_methods.h"


//#include "mjm_cli_ui.h"

//#include "mjm_tokenized_collections.h"

#include "mjm_pawnoff.h"
#include "mjm_blob.h"
//#include "mjm_misc_parse.h"
//#include "mjm_bibtex_entry.h"
//#include "mjm_bibtex_parse.h"

#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>
#include <stdint.h>

/*
made from wizard script on Fri Sep 13 14:50:58 EDT 2019
toobib
*/


////////////////////////////////////////////////////////////////


namespace ext_things_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
//typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // toobib_traits
///////////////////////////////////////////////////////////////



template <class Tr> 
class mjm_ext_things
{

typedef mjm_ext_things<Tr> Myt;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;
//typedef typename Tr::MyBlock  MyBlock;

typedef unsigned char Ch;

typedef mjm_pawnoff<Tr>  Handler;
typedef mjm_blob<Tr> Blob;
typedef mjm_blob_map<Tr>  BlobMap;
//typedef mjm_misc_parse<Tr> ParseTable;
//typedef mjm_read_buffer<Tr> RdBuf;
typedef std::map<StrTy,StrTy> Gmap; 

typedef mjm_read_buffer<Tr> RdBuf;

public:
mjm_ext_things() { Init(); } 
StrTy today(const IdxTy flags=0) { return m_hand.today(flags); } 
StrTy get_clipboard() { return Get(m_rdclip); } 
IdxTy get_clipboard(Blob & b) { Get(b,m_rdclip); return 0;  } 
// TODO this needs to insure the program does not terminate before bg
// proc launches 
IdxTy put_clipboard(const StrTy & p) { Put(p,m_copy); sleep(1);  return 0;  } 
IdxTy paste(const StrTy & p) { Put(p,m_copy); sleep(1);  return 0;  } 

StrTy get_url(const StrTy url,const IdxTy flags ) { return RdURL(url,flags); } 
IdxTy get_url(Blob & b, const StrTy & url,const IdxTy flags) { return RdURL(b,url,flags);} 
IdxTy render_html(Blob & bd, const Blob & bs, const StrTy & base) { return RenderHtml(bd,bs,base);} 

StrTy html_info(const StrTy & s, const StrTy & base) {Blob bd;  HtmlInfo(bd,Blob(s),base); return bd;} 
IdxTy html_info(Blob & bd, const Blob & bs, const StrTy & base) { return HtmlInfo(bd,bs,base);} 

IdxTy exif_info(Blob & bd, const Blob & bs, const StrTy & flags) { return ExifInfo(bd,bs,flags);} 

IdxTy pdftotext(Blob & bd, const Blob & bs, const StrTy & flags) { return PdfToText(bd,bs,flags);} 


IdxTy test_latex(Blob & bd, const Blob & bs, const StrTy & flags) { return TestLatex(bd,bs,flags);} 

IdxTy grep(Blob & bd, const Blob & bs, const StrTy & flags) { return Grep(bd,bs,flags);} 
IdxTy sed(Blob & bd, const Blob & bs, const StrTy & flags) { return Sed(bd,bs,flags);} 

// do not use often makes a mess of tracking deps
IdxTy cmd(Blob & bd, const Blob & bs, const StrTy & flags) { return Cmd(bd,bs,flags);} 




StrTy dump(const IdxTy flags=0) const { return Dump(flags); } 
//StrTy dump(const IdxTy flags=0) const { return Dump(flags); } 
private:
bool Bit(const IdxTy f, const IdxTy b) const { return ((1<<b)&f)!=0; } 
IdxTy RenderHtml(Blob & bd, const Blob & bs, const StrTy & base) {
bd=(m_hand.get_output_with_input(StrTy(bs),m_render+base));
 return 0;

} 

IdxTy HtmlInfo(Blob & bd, const Blob & bs, const StrTy & base) { 
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump));
return  0; // HtmlInfo(bd,bs,base);
} 

IdxTy ExifInfo(Blob & bd, const Blob & bs, const StrTy & flags) { 
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
const StrTy cmd=m_exif+StrTy("  ")+flags; 
//MM_ERR(MMPR(cmd));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
(m_hand.get_output_with_input(bd,bs,cmd));
return  0; // HtmlInfo(bd,bs,base);
} 


IdxTy PdfToText(Blob & bd, const Blob & bs, const StrTy & flags) { 
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
//MM_ERR(" trying pdf to text ")
const StrTy cmd=m_pdftotext+StrTy("  ")+flags; 
//MM_ERR(MMPR(cmd));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
(m_hand.get_output_with_input(bd,bs,cmd));
return  0; // HtmlInfo(bd,bs,base);
} 

IdxTy Grep(Blob & bd, const Blob & bs, const StrTy & flags) { 
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
const StrTy cmd=m_grep+StrTy("  ")+flags; 
//MM_ERR(MMPR(cmd));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
(m_hand.get_output_with_input(bd,bs,cmd));
return  0; // HtmlInfo(bd,bs,base);
} 

IdxTy Sed(Blob & bd, const Blob & bs, const StrTy & flags) { 
return Canned(m_sed,bd,bs,flags); } 

IdxTy Canned(const StrTy & base, Blob & bd, const Blob & bs, const StrTy & flags) { 
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
const StrTy cmd=base+StrTy("  ")+flags; 
//MM_ERR(MMPR(cmd));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
(m_hand.get_output_with_input(bd,bs,cmd));
return  0; // HtmlInfo(bd,bs,base);
} 




IdxTy Cmd(Blob & bd, const Blob & bs, const StrTy & flags) { 
const StrTy cmd=flags; 
(m_hand.get_output_with_input(bd,bs,cmd));
return  0; // HtmlInfo(bd,bs,base);
} 




/*

XL="pdflatex -file-line-error -halt-on-error "
# ampersand bibs are making it threw... 
#$XLBM "$fntest" > /dev/null
$XL "$fntest" &> checkbib_test_output.xxx
rczed=$?
cp "$fn" "$fn.bib"
bibtex "$fntestb" &>> checkbib_test_output.xxx
rcone=$?
$XL "$fntest" &>> checkbib_test_output.xxx
rctwo=$?
# echo rc= $?


*/
// the blob needs to concat... read_buffer???
IdxTy TestLatex(Blob & bd, const Blob & bs, const StrTy & flags) { 
StrTy m_latex,m_bibtex;
m_latex="pdflatex -file-line-error -halt-on-error ";
m_bibtex="bibtex ";
StrTy dir="/tmp";
StrTy fn="asdfqwert";
//MM_ERR(MMPR((m_html_dump+m_info_filter)));
//MM_ERR(" trying to test pdf and bibtex ")
// get a temp dir and file name 
const StrTy xxx=
StrTy("cat - >") +dir+"/"+fn
+";cd "+dir+";rm mybibtemp.bib ;rm "+fn+StrTy(".pdf")
+" ;"+m_latex+" "+fn+" ;"+m_bibtex+" "+fn
  +";"+m_latex+" "+fn
+StrTy(";rc=$? ;echo RETURNCODE= $rc ")
+";cat "+fn+".log ;pdftotext "+fn+".pdf - ";
const StrTy cmd=xxx+flags; 
//MM_ERR(MMPR(cmd));
//bd=(m_hand.get_output_with_input(StrTy(bs),m_html_dump+m_info_filter));
(m_hand.get_output_with_input(bd,bs,cmd));


return  0; // HtmlInfo(bd,bs,base);
} 






void AddHdr(StrTy & x, const Gmap & m)
{
MM_LOOP(ii,m) { 
x=x+StrTy(" --header=\"")+(*ii).first+":"+(*ii).second +"\"";

}
}

StrTy RdURLS(const StrTy uin,const IdxTy flags)
{
if (Bit(flags,0)) { return m_get_render+uin; } 
const bool dbg=(Bit(flags,1)) ;
const bool pcmd=(Bit(flags,2)) ;

const StrTy dbgflags=dbg?"-S -v":"-q";
//StrTy h1="Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8";
//StrTy h2="--user-agent=\"Mozilla\" ";
StrTy h2="--user-agent=\""+m_ua+"\" ";
//StrTy h3="Referer: "+uin;

StrTy cook=" --load-cookies="+m_cc+" --save-cookies="+m_cc+" --keep-session-cookies";
Gmap &  m = m_headers;
//m["Accept"]="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8";
m["Referer"]=uin;
//StrTy x=m_rdurl+ " -O \"-\" -S -v";
StrTy x=m_rdurl+ " -O \"-\" "+dbgflags;
//+StrTy(" --header=\"")+h1 +"\""
AddHdr(x,m);
x=x+cook+" "
//+StrTy(" --header=\"")+h3+"\" "
+h2+" \""+uin+"\"" ; 

if (pcmd) {MM_ERR(MMPR(x)) } 
return(x);
}
IdxTy RdURL(Blob & b, const StrTy & url, const IdxTy flags) // { return RdURL(b,url);} 
{
const StrTy cmd=RdURLS(url,flags);
Get(b,cmd);
return 0;
}
StrTy RdURL(const StrTy u,const IdxTy flags)
{
const StrTy cmd=RdURLS(u,flags);
Blob blob;
Get(blob,cmd);
return StrTy(blob); 
}

void Init()
{
m_debug=0;
m_cc="crapcookie.txt";
m_ua="Mozilla";
m_rdclip="myclip -echo";

// keep these vars 
//m_render="lynx -force_html -dump -width=10000 -useragent=\""+m_ua+"\" -stdin -base ";
m_render="lynx -force_html -dump -width=10000 -useragent=\"Mozilla\" -stdin -base ";
//m_get_render="lynx -force_html -dump -width=10000 -useragent=\"Mozilla\" -accept_all_cookies -cookie_file="+m_cc+" -cookie_save_file="+m_cc+" ";
m_get_render="lynx -force_html -dump -width=10000 -useragent=Mozilla -accept_all_cookies -cookie_file="+m_cc+" -cookie_save_file="+m_cc+" ";


m_rdurl="wget "; 
m_html_dump="testHTML -quiet -mjm - ";
// this terminates when someone else pasts to clipboard- need to
// fork off. 
//m_copy="xclip -verbose -in -l 0  -selection clipboard ";
m_copy="xclip -in -l 0  -selection clipboard &  ";
//m_copy="xclip -verbose -in -l 1 ";
//m_copy="cat -n -";
//m_info_filter="| grep "\( meta \| href.*text \|text 10\.\)";
m_info_filter="| grep \"\\( meta \\| href.*text \\|text 10\\.\\)\"";
m_exif="exifutil -list - ";
m_pdftotext="pdftotext - - ";
m_grep="grep ";
m_sed="sed ";
//m_base="basename";
//m_dir="dirname";

//cat xxx | sed -s 's/^m_/m_lut["/g' | sed -e 's/=/"]=/'
m_lut["cc"]="crapcookie.txt";
m_lut["ua"]="Mozilla";
m_lut["rdclip"]="myclip -echo";
//m_lut["render"]="lynx -force_html -dump -width=10000 -useragent=\""+m_ua+"\" -stdin -base ";
//m_lut["get_render"]="lynx -force_html -dump -width=10000 -accept_all_cookies -cookie_file="+m_cc+" -cookie_save_file="+m_cc+" ";
m_lut["get_render"]="lynx -force_html -dump -width=10000 -useragent=\"Mozilla\" -accept_all_cookies -cookie_file=crapcookie.txt -cookie_save_file=crapcookie.txt ";
m_lut["rdurl"]="wget ";
m_lut["html_dump"]="testHTML -mjm - ";
m_lut["copy"]="xclip -in -l 0  -selection clipboard &  ";
m_lut["info_filter"]="| grep \"\\( meta \\| href.*text \\|text 10\\.\\)\"";
m_lut["exif"]="exifutil -list - ";
m_lut["pdftotext"]="pdftotext - - ";
m_lut["grep"]="grep ";
m_lut["sed"]="sed ";
Gmap&  m=m_headers;
m["Accept"]="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8";
//m["Referer"]=uin;


}
StrTy Get(const StrTy &  x) {return m_hand.get_output("",x);}
StrTy Put(const StrTy&  x, const StrTy & cmd) 
{return m_hand.get_output_with_input_bg(x,cmd);}
void Get(Blob & blob,const StrTy &  x) {m_hand.get_output(blob,"",x);}

StrTy Dump(const IdxTy flags) const
{
Ss ss;
ss<<MMPR4SL(m_cc,m_ua,m_rdclip,m_render);
ss<<MMPR4SL(m_get_render,m_rdurl,m_copy,m_html_dump);
ss<<MMPRSL(m_info_filter);
ss<<MMPRSL(m_exif);
//ss<<MMPRSL(m_base);
//ss<<MMPRSL(m_dir);
ss<<MMPRSL(m_pdftotext);
ss<<MMPRSL(m_grep);
ss<<MMPRSL(m_sed);
MM_LOOP(ii,m_lut) { ss<<" m_lut "<<(*ii).first<<" = "<<(*ii).second<<CRLF; }
return ss.str();
}

Handler m_hand;
IdxTy m_debug;
// this works but is not flexible
StrTy  m_cc,m_ua,m_rdclip,m_render,m_get_render,m_rdurl,m_copy;
StrTy m_html_dump,m_info_filter,m_exif,m_pdftotext,m_grep,m_sed;
//typedef std::map<const Ch * ,const Ch * > Lut;
// this probably works but is still limited
typedef mjm_char_star_map<Tr , const char * ,const char * > Lut;
Lut m_lut;
Gmap m_headers;
// this probably works best with strings in a single buffer and made at runtime
// requiring rebuild when one is changed but this is not a big event. 
typedef std::map<StrTy, IdxTy > BufTokenMap;
BufTokenMap m_btm;
RdBuf m_strings;
//StrTy m_base,m_dir;

}; // ext_things
//typedef ext_things ExtThings;




#endif
