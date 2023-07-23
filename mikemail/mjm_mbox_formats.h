#ifndef MJM_MBOX_FORMATS_H__
#define MJM_MBOX_FORMATS_H__



#include <mjm_generic_message.h>
#include <mjm_paranoid_ofstream.h>
#include <mjm_read_buffer.h>
#include <mjm_instruments.h>


template <class Tr>
class mjm_mbox_formats 
{
 typedef mjm_mbox_formats Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
enum {BAD=Tr::BAD}; 
//typedef typename mjm_generic_message<Tr>  MyMsg;

// these should be completely in Tr wtf
typedef  mjm_generic_message<Tr>  MyMsg;
typedef mjm_pawnoff<Tr> Hand;
typedef mjm_paranoid_ofstream<Tr> Pos;
typedef mjm_read_buffer<Tr> Rb;
public:
mjm_mbox_formats(): m_attach_pfx("./attach/") {}

bool is_email(const StrTy & s)
{
const char * p=s.c_str();
while ( *p!=0)
{
if (*p=='@') return true;
++p;
}
return false;
}
template <class Ty>
StrTy find_address(Ty & v)
{
MM_LOOP(ii,v)
{
if (is_email(*ii)) return (*ii);
}
return "xxx@xxx";
}

IdxTy decode_parts(MyMsg & msg, const IdxTy flags=0)
{
const IdxTy sz=msg.parts();
for(IdxTy i=0; i<sz; ++i)
{
MyMsg & part=msg.part(i);
const StrTy & mime=part.mime();
// MM_ERR(MMPR3(i,mime,flags))
decode_parts(part,flags+1);
/*
Content-Type: image/gif^M
Content-ID: <wrt.1225515037.3587.87814.w125@yahoogroups.com>^M
Content-Transfer-Encoding: base64^M
Content-Disposition: inline; filename="bg_lblue_right.gif"^M
*/
StrTy  uuid=part.server_uuid();
Ss su;
su<<uuid<<"__"<<flags<<"_"<<i;
//MM_ERR("WTF"<<MMPR3(uuid,flags,i))
uuid=su.str();
const StrTy & ctype=part.hval_lc("content-type:");
const StrTy & cte=part.hval_lc("content-transfer-encoding:");
const StrTy & cd=part.hval_lc("content-disposition:");
//const StrTy & disp=part.hval_lc("content-disposition:");
StrTy fn=part.field_in_field("content-disposition:","filename");
if (fn.length()==0)  fn=part.field_in_field("content-disposition:","name");
//fn=fn+uuid;
if ( cte=="base64")
{
//MM_ERR(MMPR4(fn,cte,cd,ctype))
char * result;
IdxTy len;
Pos fos;
if (fn.length()!=0)
{
//const StrTy fn_ =fos.friendly_file(fn);
//const StrTy fn_ =fos.friendly_unique_file(fn+uuid);
const StrTy fn_ =fos.friendly_file(fn+uuid);
//const StrTy fntotal=StrTy("./attach/")+fn_;
//const StrTy fntotal=fos.unique("./attach/",fn_,"");
const StrTy fntotal=fos.unique(m_attach_pfx,fn_,"");
fos.set_file(fntotal);
//MM_ERR(" file set "<<MMPR(fntotal))
}
else
{
const StrTy ext=mime_to_ext(msg.mime());
//fos.set_file("./attach/","unk",StrTy("._")+ext);
//fos.set_file(fos.unique("./attach/",uuid,StrTy("._")+ext));
fos.set_file(fos.unique(m_attach_pfx,uuid,StrTy("._")+ext));

//MM_ERR(" alt file set "<<MMPR(ext))
}
// these writes should not need locks 
m_hand.base64(&result, &len,part.body());
//MM_ERR(MMPR2(fos.fn(),len))
fos.write_nothrow(result,len);
fos.finish();

}


}

return 0;
}
//MyMsg next_mbox_msg(std::istream & ifs)
// can return alt, just as easy tis way
// but exta assign is fcking 
/*
MyMsg next_mbox_msg(LineIterator & li )
{
MyMsg m;
 next_mbox_msg(m,LineIterator & li )
return m;
}
*/
MyMsg next_mbox_msg( LineIterator & li )
{//LineIterator li(&ifs);
MyMsg m;
StrTy s,a,uuid;
IdxTy d=0;
bool skipping=true;
m_rb.clear();
IdxTy i=0; 
while (li.nextok())
{
//Ss ss;
++i; 
const IdxTy sz=li.words().size();
//std::cerr<<MM_MARK<<" "<<MMPR2(i,sz)<<"                        "<<'\r';
if (sz==0) { if (skipping) continue;
else {m_rb.append('\n');   //ss<<'\n'; s=s+ss.str(); 
continue; }   } 
const StrTy & cmd=li.word(0);
if ( cmd=="From")
{
skipping=!skipping;
if ( skipping) 
{
li.save_last_line(true);
// put the from line back. 
// mjm_generic_message(const StrTy & v,  const StrTy & a, const IdxTy d, const StrTy & uuid)
s=m_rb.string();
MyMsg m2(s,a,d,uuid);
s="";
m_rb.clear();
return m2;
//break;
} //skipping 

} // from 
//else { ss<<li.line()<<'\n'; s=s+ss.str(); } 
else { m_rb.append(li.line()); m_rb.append('\n');} 

} // while 
//if ( s!="")
if (m_rb.size()!=0)
{
s=m_rb.string();
MyMsg m2(s,a,d,uuid);
return m2;


}
return m; 
}



IdxTy  to_mbox(StrTy & v, MyMsg & msg, const IdxTy & flags=0 )
{
const bool include_sfs=!Bit(flags,0); // true;
const StrTy _v=msg.verbatim();
StrTy d=msg.hval_lc("Date:");
StrTy from=msg.hval_lc("From:");
StrTy to=msg.hval_lc("To:");
StrTy subj=msg.hval_lc("Subject:");
StrTy sfs=StrTy("Mikemail-Source-Info: ")+msg.string_flags()+StrTy(" ")+ msg.string_data();
//MM_ERR(MMPR4(d,from,to,subj)<<MMPR(sfs))
bool have_from=false;
Hand y;
Ss ss;
ss<<"From: "<<from<<CRLF;
ss<<"To: "<<to<<CRLF;
ss<<"Date: "<<d<<CRLF;
ss<<"Subject: "<<subj<<CRLF;
ss<<_v;
LineIterator li(&ss);
li.set_split(1,' ');
li.remove_crlf(true);
// from 
li.nextok(); 
 from=find_address(li.words());  
li.nextok();  to=li.line();
li.nextok();  IdxTy sz=li.words().size();    d=li.reassemble(1,sz); d=y.date_default(d);   
li.nextok();   subj=li.line(); 
 
Ss sr;
IdxTy modified=0;
while (li.nextok())
{

const IdxTy sz=li.size();
const StrTy cmd=li.word(0);
if (sz>=1){
if (cmd=="From" ) {  sr<<">"; ++modified; } 
} // sz 
sr<<li.line()<<'\n';
} // nextok()
bool bad_parse= ((d.length()==0) || ( from.length()==0));
if (bad_parse) //  ((d.length()==0) || ( from.length()==0))
{
MM_ERR(" date or from mising "<<MMPR2(d,from))
}
//const StrTy ds= nocrlf(hand.dates(d));
//sr<<_v;
Ss sfrom;
//sfrom<<"From "<<from<<" "<<d<<'\n';
sfrom<<"From "<<from<<" "<<d;
//sfrom<<to<<'\n';
//sfrom<<subj<<'\n';
if (include_sfs)  sfrom<<sfs<<'\n';
v=sfrom.str();
if (modified !=0)
{
MM_ERR(" message modified  "<<MMPR3(d,from,modified))
sr<<'\n'<<" [ Modified by mikemail for mbox format "<<modified<<" times. ]\n";
v+=sr.str()+StrTy("\n");
}
else v+=_v+StrTy("\n"); 
// this is really a save parts thing as they have already been decoded 
decode_parts( msg);
return bad_parse?BAD:modified;

} // to_mbox
#if 0 
IdxTy  convert_to_mbox(StrTy & v, MyMsg & msg, const StrTy & flags )
{
const bool include_sfs=!Bit(flags,0); // true;
const StrTy & _v=msg.verbatim();
// TODO FIXME this needs to concatenate headers... 
Ss ss(_v);
Ss sr;
Hand y;
StrTy d="";
StrTy from="";
StrTy to="";
StrTy subj="";
bool have_from=false;
//CommandInterpretter li(ss);
LineIterator li(&ss);
li.set_split(1,' ');
li.remove_crlf(true);
IdxTy modified=0;
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
//if (m_flp.log_commands()) 
//		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz>=1){
const StrTy cmd=li.word(0);

//if (cmd=="From:" ) {   } 
//if (from.length()==0) if (cmd=="From:" ) {   StrTy din="" ; for(IdxTy j=1; j<sz; ++j) din=din+StrTy(" ")+li.word(j); from=din; }  
//if (from.length()==0) if (cmd=="From:" ) {   from=li.reassemble(1,sz);  }  
if (from.length()==0) if (cmd=="From:" ) {   from=find_address(li.words());  }  
if (cmd=="From" ) {  sr<<">"; ++modified; } 
//if (d.length()==0) if (cmd=="Date:" ) {   StrTy din="" ; for(IdxTy j=1; j<sz; ++j) din=din+StrTy(" ")+li.word(j); d=y.date_default(din); } 
if (d.length()==0) if (cmd=="Date:" ) {   d=li.reassemble(1,sz); d=y.date_default(d);  } 
if (to.length()==0) if (cmd=="To:" ) {   to=li.line(); } 
if (subj.length()==0) if (cmd=="Subject:" ) {   subj=li.line(); } 
 }
sr<<li.line();
} // nextok()
bool bad_parse= ((d.length()==0) || ( from.length()==0));
if (bad_parse) //  ((d.length()==0) || ( from.length()==0))
{
MM_ERR(" date or from mising "<<MMPR2(d,from))
}
//const StrTy ds= nocrlf(hand.dates(d));
//sr<<_v;
Ss sfrom;
//sfrom<<"From "<<from<<" "<<d<<'\n';
sfrom<<"From "<<from<<" "<<d;
sfrom<<to<<'\n';
sfrom<<subj<<'\n';
if (include_sfs)  sfrom<<sfs<<'\n';
v=sfrom.str();
if (modified !=0)
{
MM_ERR(" message modified  "<<MMPR3(d,from,modified))
sr<<'\n'<<" [ Modified by mikemail for mbox format "<<modified<<" times. ]\n";
v+=sr.str()+StrTy("\n");
}
else v+=_v+StrTy("\n"); 
return bad_parse?BAD:modified;
} // convert_to_mbox
#endif

 StrTy mime_to_ext(const StrTy & mime)
{
// to_lower
LineIterator li;
li.set_split(3,'/'); // preserve quotes for second pass
li.remove_crlf(true);
li.nextline(mime.c_str());
const auto & w=li.words();
const IdxTy sz=w.size();
StrTy type="";
if (sz==0) return "xxx";
if (sz==1) return w[0];
if (sz==2) return w[1];
return li.reassemble(0,sz,".");

} // mime_to_ext


private:

Hand m_hand; 
Rb m_rb;

StrTy m_attach_pfx;

bool Bit(const IdxTy flag, const IdxTy b ) const { return ((flag&(1<<b))!=0); } 








}; // mjm_mbox_formats

#endif // MJM_MBOX_FORMATS_H__ 
