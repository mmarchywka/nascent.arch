#ifndef MJM_GENERIC_MESSAGE_H__
#define MJM_GENERIC_MESSAGE_H__

// need to have includes here doh 

#include "mjm_pawnoff.h"
#include "mjm_strings.h"
#include "mjm_paranoid_ofstream.h"
#include <mjm_read_buffer.h>
#include <mjm_collections.h>

#include <vector>
#include <map>


template <class Tr>
class mjm_generic_message
{
//typedef mikemail_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_generic_message Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::Ofs Ofs;
//typedef typename Tr::MyBlock  MyBlock;
enum {BAD=Tr::BAD};
//typedef pthread_mutex_t Mutex;
//typedef typename Tr::Mutex Mutex;
typedef std::vector<StrTy > KeySpec;
typedef std::map<StrTy, StrTy> KeyMap;
typedef mjm_ragged_table Ragged;
typedef  std::vector<StrTy> KeyList; 
typedef mjm_pawnoff<Tr> Hand;
typedef mjm_read_buffer<Tr> Rb;

typedef std::vector<Myt*> Parts;
typedef std::vector<Myt*> MsgList;
typedef std::map<StrTy, IdxTy> AltRank;

typedef std::vector<StrTy> Vvalues;
typedef std::map<StrTy,Vvalues> Hvalues;
typedef std::map<StrTy, Hvalues> VectorHeaders;

public:
typedef  AltRank alternative_ranks;
mjm_generic_message() {}
~mjm_generic_message() {Release(); }
//mjm_generic_message(const StrTy & h, const StrTy & b) :m_hdr(h),m_body(b)  {}
mjm_generic_message(const StrTy & h, const StrTy & b,const StrTy & a)
:m_hdr(h),m_body(b)  {Init(); MakeKeys(h);  m_keys["account"]=a; MakeLC(); }
//MyMsg msg=MyMsg(h,b,v,sfs,a);
// changed to niclde verbatim 2020-07-18 may be some reason not to 
mjm_generic_message(const StrTy & h, const StrTy & b,const StrTy & v, const StrTy & sfs, const StrTy & a)
//:m_hdr(h),m_body(b)  {Init(); MakeKeys(h);  m_keys["account"]=a;  MakeLC(); }
:m_hdr(h),m_body(b),m_verbatim(v)  {Init(); MakeKeys(h);  m_keys["account"]=a;  MakeLC(); }
mjm_generic_message(const StrTy & v,  const StrTy & a)
:m_verbatim(v)  {Init(); MakeKeys(v);  m_keys["account"]=a; MakeLC();  SplitPayload();  }
mjm_generic_message(const StrTy & v,  const StrTy & a, const IdxTy d)
:m_verbatim(v)  {Init(); depth(d);  MakeKeys(v);  m_keys["account"]=a; MakeLC();  SplitPayload();  }
// TODO FIXME this needs a boundary stack for hierarchial attachments 
mjm_generic_message(const StrTy & v,  const StrTy & a, const IdxTy d, const StrTy & uuid)
:m_verbatim(v)  {Init(); depth(d);  MakeKeys(v);  m_keys["account"]=a; server_uuid(uuid);  MakeLC();  SplitPayload();  }
// AddPart
//void AddPart(const StrTy & v) { m_parts.push_back( new Myt(v,account(),depth()+1,server_uuid())); }
mjm_generic_message(const Myt & p, const StrTy & v)
:m_verbatim(v)  {Init(); 
depth(p.depth()+1);  
// TODO this logic has been fixed IIRC need to remove junk from here 
MakeKeys(v);  
m_keys["account"]=p.account(); 
server_uuid(p.server_uuid());  
MakeLC();  
SplitPayload();  
}

void remake_headers(const IdxTy flags=0)
{

MakeKeys(m_verbatim);  
MakeLC();  
SplitPayload();  

}

// that's great but the message buffer uses a f'ing map that makes a f'ing copy with 
// that's great but the message buffer uses a f'ing map that makes a f'ing copy with 
// copy ctor and ASSign operator fing stuff

mjm_generic_message( const Myt & x)
{
Assign(*this,x);
}
Myt & operator=(const Myt & x)
{
Assign(*this,x);
return *this;
}

void my_data(const StrTy & s) { m_my_data=s; } 
void flag_string(const StrTy & sfs) { m_sfs=sfs; } 
const StrTy&  string_flags() const  {return  m_sfs; } 
const StrTy&  string_data() const  {return  m_my_data; } 

//MyMsg msg=MyMsg(h,b,v,sfs,a);

const StrTy & header() const { return m_hdr;}
const StrTy & body() const { return m_body;}
const StrTy  decoded_body() const {Ss ss=DecodeBody();  return ss.str();}
const StrTy & verbatim() const { return m_verbatim;}
StrTy from_address(const IdxTy flags ) const { return FromAddress(flags); } 
// should use something with RO [] 
const StrTy & account() { return m_keys["account"]; } 
const StrTy & server_uuid()const  { return m_keys["uuid"]; }
void  server_uuid(const StrTy & val) { add_key("uuid", val); }

const StrTy & mime() const { return m_mime; } 
const IdxTy parts() const { return m_parts.size(); } 
const IdxTy depth() const { return m_depth; } 
const IdxTy depth(const IdxTy i ) { m_depth=i;  return m_depth; } 
// ref should return as it is just a ptr 
Myt& part(const IdxTy i )  { return *m_parts[i]; } 
const Myt& part(const IdxTy i )const   { return *m_parts[i]; } 

const StrTy name() const
{
StrTy type;
StrTy f1=field_in_field("content-disposition:","filename",&type);
if (f1.length()==0) return type;
return f1;
}


// only works for lexical string no notion of numbers. 
//void make_key( const KeySpec & k )
//{

//}
// for making a header line. 
StrTy header_string( const KeyList & w) const
{
Ss ss;
MM_LOOP(ii,w)
{
const StrTy xx=m_keys[(*ii)];
// this looke ok 
//MM_ERR(MMPR4(xx.length(),strlen(xx.c_str()),(*ii),m_keys[(*ii)])) 
//ss<<m_keys[(*ii)]<<"+"; 
ss<<xx.c_str()<<"+"; 
//MM_ERR(MMPR(ss.str()))
} // ii 
// this not ok 
return ss.str();


}


template <class Tv> void  dump_keys(Tv & v) const
{
Ss ss;
MM_LOOP(ii,m_keys)
{
//ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
ss.push_back((*ii).first); // <<" "<<(*ii).second<<CRLF;
}
//return ss.str();
} // dump_keys



StrTy dump_keys() const
{
Ss ss;
MM_LOOP(ii,m_keys)
{
ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
return ss.str();
} // dump_keys


void split_payload()
{
SplitPayload();


}

static const StrTy&  blank_key()  { static StrTy  blank=""; return blank; } 
void add_key(const StrTy & key, const StrTy & val)
{
m_keys[key]=val;
StrTy _n = mjm_strings::fancy_to_lower(key);
m_keys_lc[_n]=val;
m_v_headers_lc[_n][key].push_back(val);
m_v_headers_lc[_n][blank_key()].push_back(val);

}

/*
void MakeLC()
{
MM_LOOP(ii,m_keys)
{
const StrTy k=mjm_strings::fancy_to_lower((*ii).first);
m_keys_lc[k]=(*ii).second;

} // m_keys

}

*/



// this only works for uniq keys 
const StrTy hval(const StrTy & n) const { return m_keys[n]; } 
const StrTy hval_lc(const StrTy & n) const 
{ 
StrTy _n = mjm_strings::fancy_to_lower(n);
return m_keys_lc[_n]; 

} 
StrTy sane_headers() const
{
Ss ss;
MM_LOOP(ii,m_v_headers_lc)
{
const auto &  vv=((*ii).second)[blank_key()];
MM_LOOP(jj,vv) { ss<<(*ii).first<<" "<<(*jj)<<CRLF; } // jj
} // ii 

return ss.str();
} // sane_headers

const Vvalues  &  hval_vec(const StrTy & n) const 
{ 
//Vvalues vv;
StrTy _n = mjm_strings::fancy_to_lower(n);
//return m_keys_lc[_n]; 

const auto & vv=m_v_headers_lc[_n];
//MM_ERR(MMPR3(m_v_headers_lc.size(),vv.size(),_n))
const auto ii=vv.find(blank_key());
if (ii==vv.end())
{
MM_ERR(" found nothing for "<<MMPR(n))
static Vvalues vvv;
return vvv;
}
//return vv[blank_key()]; 
return (*ii).second; 
//return (*(vv.find(blank_key()))).second; 
/*MM_LOOP(ii,m_v_headers_lc)
{
const auto j=(*ii).second;
MM_LOOP(jj,j) { MM_LOOP(kk,(*jj).second) vv.push_back((*kk)); }

}
*/


//return vv;
} 




const Ragged & parsed() const { return m_parsed; } 
StrTy 
field_in_field(const StrTy & hdr, const StrTy & name, StrTy * p=0 )
const 
{
const StrTy & ct=hval_lc(hdr);
return FieldInField(ct,name,p);
} 

void viewable_parts(MsgList & l, AltRank & map )
{
const IdxTy sz=parts();
const StrTy & mt=mime();
if (sz==0) { if ( strncmp(mt.c_str(),"text",4)==0) l.push_back(this); return ; } 
if ( mt=="multipart/alternative")
{
IdxTy besti=BAD;
IdxTy bestr=BAD;
IdxTy defi=BAD;
StrTy def="text/plain";
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & mtp=part(i).mime();
//if (mtp=="text/plain") {l.push_back(&part(i)); return; } 
if (mtp==def) {defi=i; } 
IdxTy rank=map[mtp];
if (rank!=0) if ((rank<bestr)||(bestr==BAD)) { besti=i; bestr=rank; }
} // for i 
if ( bestr!=BAD)  {l.push_back(&part(besti)); return; } 
if ( defi!=BAD)  {l.push_back(&part(defi)); return; } 
l.push_back(&part(0)); 
return;  
} // alternatives  

for(IdxTy i=0; i<sz; ++i) { part(i).viewable_parts(l,map); } // i 

} // viewable_parts 

StrTy text(AltRank & map) 
{
Ss ss;
MsgList l;
viewable_parts(l,map);
MM_LOOP(ii,l)
{
//const StrTy code= (*ii)->hval_lc("Content-Transfer-Encoding:");
ss<<(*ii)->body();
ss<<CRLF;
} // ii 
return ss.str();
} // text

template <class Ty> StrTy text(Ty & hand, AltRank map ) 
{
Ss ss;
MsgList l;
viewable_parts(l,map);
MM_LOOP(ii,l)
{
//Ss sr;
const StrTy mime=(*ii)->mime();

Ss sr=(*ii)->DecodeBody();
//const StrTy code= (*ii)->hval_lc("Content-Transfer-Encoding:");
//if ( code=="base64") { sr<<hand.base64((*ii)->body()); } 
//else if ( code=="quoted-printable") { sr<<hand.qp((*ii)->body()); } 
//else { sr<<(*ii)->body(); } 

if (mime=="text/html")
{ ss<<hand.html(sr.str()); }
else { ss<<sr.str(); } 
ss<<CRLF;
} // ii 
return ss.str();
} // text


void server_tag(void * p ) { m_p_msg_server=p; } 
void *  server_tag( )const  {return  m_p_msg_server; } 


private:

StrTy m_hdr;
StrTy m_body;
StrTy m_verbatim;
StrTy m_sfs;
StrTy m_my_data;
StrTy m_mime;
StrTy m_boundary;
std::vector<StrTy> m_boundary_stack;
Parts m_parts;
IdxTy m_depth;

Ragged m_parsed;
void * m_p_msg_server;

// dumbly extracted from headers AND account info 
//KeyMap m_keys;
mutable ReadWriteMap m_keys,m_keys_lc;

mutable VectorHeaders m_v_headers_lc;
// threading???? 
// made prior to sorting although should be tokeninzed
// to something easy to sort. 
mutable StrTy m_sort_key;

void Init()
{
m_depth=0;
m_p_msg_server=0;
}


static void Assign(Myt & d, const Myt & s )
{
// do nothing really as the parts get parse in pieces. 
d.m_hdr=s.m_hdr;
d.m_body=s.m_body;
d.m_verbatim=s.m_verbatim;
d.m_sfs=s.m_sfs;
d.m_my_data=s.m_my_data;
d.m_mime=s.m_mime;
d.m_boundary=s.m_boundary;
d.m_boundary_stack=s.m_boundary_stack;
d.Release();
DeepCopy(d,s);
d.m_depth=s.m_depth;
d.m_parsed=s.m_parsed;
d.m_keys=s.m_keys;
d.m_keys_lc=s.m_keys_lc;
d.m_sort_key=s.m_sort_key;
//m_v_headers_lc[k][key].push_back(val);
// TODO FIXME someone is making an asignment instead of a reference 
// in the proc mail code
d.m_v_headers_lc=s.m_v_headers_lc;
d.m_p_msg_server=s.m_p_msg_server;
}

// this does not pickup the keys. wtf 
// TODO FIXME hieratchial boundary stack... 
void AddPart(const StrTy & v) { m_parts.push_back( new Myt(v,account(),depth()+1,server_uuid())); }


static bool Bit(const IdxTy flag, const IdxTy b )  { return ((flag&(1<<b))!=0); }

// this is really wasteful. In reality there are no copied needed just handing off the
// pointers would be ok 
void Release()
{
MM_LOOP(ii,m_parts) { delete (*ii); } 
m_parts.clear();
}
static void DeepCopy( Myt & d, const Myt & s)
{
MM_LOOP(ii,s.m_parts)
{
Myt * p= new Myt();
d.m_parts.push_back(p);
Assign(*p,**ii);

}

}

void MakeLC()
{
MM_LOOP(ii,m_keys)
{
const StrTy & key=(*ii).first;
const StrTy k=mjm_strings::fancy_to_lower((*ii).first);
const StrTy & val=(*ii).second; 
m_keys_lc[k]=val; // (*ii).second;
/////////////////////////////////////////////////
m_v_headers_lc[k][key].push_back(val);
m_v_headers_lc[k][blank_key()].push_back(val);
//MM_ERR(" adding headers "<<MMPR3(key,k,val))
////////////////////////////////////////////////

} // m_keys

}


StrTy FromAddress(const IdxTy flags ) const
{
const bool remove_ampersand=Bit(flags,0); 
const StrTy fr= (hval_lc("From:"));
Rb sub(fr);
//sub.split_and_mark(su,' ');
//const IdxTy n=sub.strings();
//if (n<2) return HandleNull(m,dv);
//const char * cmd=sub[1];
Rb s2;
s2.take_between(sub,'<','>',' ');
IdxTy n=s2.strings();
// mayu not be in brackets, look for ampersand 
if (n==0)
{
s2.clear();
s2.append(fr); s2.cap();
s2.split_and_mark(' ');
n=s2.strings();
for(IdxTy i=0; i<n; ++i) { IdxTy xx=s2.csub('@','@',i);
if (xx==1) return StrTy(s2[i]);  
} // i 

} // n==0
if (n!=1){
 MM_ERR(" failued searching for address in "<<MMPR2(fr,n))
return fr; 
}
if (remove_ampersand) s2.csub('_','@',0);
StrTy x=s2[0];
MM_ERR(" searching from returns  "<<MMPR2(fr,x))
return x;
}

void SplitPayload()
{
GetBasics();
 const bool multipart=(m_boundary.c_str()[0]!=0);
// look for lines exactly --boundary making children from the first
// to the end retaining full text in parent.
if (multipart)
{ ParseParts(); 
} // multipart 



}

/*
Should check the speck but this looks like pattern. 
--00000000000064c0cf05862201f7
--00000000000064c0cb05862201f5
--00000000000064c0cb05862201f5
--00000000000064c0cb05862201f5--
--00000000000064c0cf05862201f7
--00000000000064c0cf05862201f7
--00000000000064c0cf05862201f7
--00000000000064c0cf05862201f7--

*/

Ss DecodeBody() const 
{
return Decode(body());
}
Ss  Decode(const StrTy & v) const 
{
Ss sr;
const StrTy code= hval_lc("Content-Transfer-Encoding:");
if ( code=="base64") {Hand hand;  sr<<hand.base64(v); } 
else if ( code=="quoted-printable") { Hand hand;  sr<<hand.qp(v); } 
else { sr<< v; } 

return sr;
}

void ParseParts()
{
Rb rb; 
rb.clear();
const StrTy hh=StrTy("--");
StrTy hit=hh+m_boundary;
StrTy end=hit+hh;
//StrTy pt="";
IdxTy part=BAD;
//Ss vs(m_verbatim);
//Ss vs(Decode(m_body));
// TODO FIXME in general this will not work if body is encoded
Ss vs(DecodeBody().str());
LineIterator li(&vs);
li.set_split(1,' ');
li.remove_crlf(true);
while (li.nextok())
{
const IdxTy sz=li.size();
const auto & w=li.words();
//void AddPart(const StrTy & v) { push_back( new Myt(v,account())); }
if ((sz==1)&&(w[0]==hit))
{
//if (part!=BAD)  AddPart(pt);
if (part!=BAD)  AddPart(rb.string());
rb.clear();  // pt="";
++part;
// this is turning into a huge time sink on large attachments  
//} else { pt=pt+li.line()+StrTy("\n"); } 
} 
else if ((sz==1)&&(w[0]==end))
{
if (part!=BAD)if (rb.size()!=0)  { ++part;  AddPart(rb.string()); rb.clear(); } 

break;  // wait, this is descending and not need to stack doh now the next one needs to know where to start
}
else { rb.append(li.line());  rb.append('\n'); } 

} // nextok
//if (part!=BAD)if (pt.length()!=0)  { ++part;  AddPart(pt); } 

// this should not exist 
/*
and created infinite recursion if found, need to check the encoding  
 v="Content-Type: text/html;\n\tboundary=\"----=_NextPart_000_001E_01D3BFB7.C22B89F0\"; charset=\"UTF-8\"\nContent-Transfer-Encoding: quoted-printable\n\n<html xmlns:v=3D\"urn:schemas-microsoft-com:vml\" xmlns:o=3D\""...) at ./mjm_generic_message.h:291
#162 0x000000000044d19a in mjm_generic_message<mikemail_traits::Tr>::ParseParts (this=0x7fff57502360) at ./mjm_generic_message.h:393
#163 0x000000000043d0e3 in mjm_generic_message<mikemail_traits::Tr>::SplitPayload (this=0x7fff57502360) at ./mjm_generic_message.h:333
#164 0x000000000042c1c4 in mjm_generic_message<mikemail_traits::Tr>::mjm_generic_message (this=0x7fff57502360, 
    v="Content-Type: text/html;\n\tboundary=\"----=_NextPart_000_001E_01D3BFB7.C22B89F0\"; charset=\"UTF-8\"\nContent-Transfer-Encoding: quoted-printable\n\n<html xmlns:v=3D\"urn:schemas-microsoft-com:vml\" xmlns:o=3D\""..., 
    a="marchywka@hotmail.com;imap-mail.outlook.com;BounceFck", d=7079, uuid="158-4346") at ./mjm_generic_message.h:59
#165 0x0000000000456f49 in mjm_generic_message<mikemail_traits::Tr>::AddPart (this=0x7fff574ef580, 
    v="Content-Type: text/html;\n\tboundary=\"----=_NextPart_000_001E_01D3BFB7.C22B89F0\"; charset=\"UTF-8\"\nContent-Transfer-Encoding: quoted-printable\n\n<html xmlns:v=3D\"urn:schemas-microsoft-com:vml\" xmlns:o=3D\""...) at ./mjm_generic_message.h:291
#166 0x000000000044d19a in mjm_generic_message<mikemail_traits::Tr>::ParseParts (this=0x7fff574ef580) at ./mjm_generic_message.h:393
#167 0x000000000043d0e3 in mjm_generic_message<mikemail_traits::Tr>::SplitPayload (this=0x7fff574ef580) at ./mjm_generic_message.h:333
#168 0x000000000042c1c4 in mjm_generic_message<mikemail_traits::Tr>::mjm_generic_message (this=0x7fff574ef580, 
    v="Content-Type: text/html;\n\tboundary=\"----=_NextPart_000_001E_01D3BFB7.C22B89F0\"; charset=\"UTF-8\"\nContent-Transfer-Encoding: quoted-printable\n\n<html xmlns:v=3D\"urn:schemas-microsoft-com:vml\" xmlns:o=3D\""..., 
    a="marchywka@hotmail.c
*/
// make sure it found at least a starting point to chop something off. 
if (part!=BAD)
if (rb.size()!=0)  { MM_ERR("  recursion mess, leftover text "<<server_uuid() )  ++part;  AddPart(rb.string()); } 

//MM_ERR(MMPR(part))




} // ParseParts
StrTy 
FieldInField(const StrTy & ct, const StrTy & name, StrTy * p=0)
const
{
LineIterator li;
//li.set_split(3,' '); // preserve quotes for second pass
li.set_split(5,' '); // preserve quotes and remove all whiet space for second pass

li.remove_crlf(true);

LineIterator lie;
lie.set_split(1,'=');
lie.remove_crlf(true);
// lie.remove_white(true); // this worked but convert to splitter 4

StrTy boundary="";
//const StrTy & ct=hval_lc("Content-Type:");
// MIME-Version: 1.0
//Content-Type: multipart/alternative; boundary="-----------------------------1152365197"

li.nextline(ct.c_str());
const auto & w=li.words();
const IdxTy sz=w.size();
StrTy type="";
if (sz>0) type=w[0];
if (p!=0) *p=type;
MM_LOOP(ii,w)
{
const auto & x=(*ii);
const IdxTy len=x.length();
char c[len+1];
memcpy(c,x.c_str(),len+1);
if (len>0) if ( c[len-1]==';') c[len-1]=0; 
lie.remove_junk(c);
//lie.nextline(x.c_str());
lie.nextline(&c[0]);
const auto & y=lie.words();
const IdxTy sze=y.size();
if ( sze<2) continue;
const StrTy& cmd=y[0];
const StrTy val=y[1];
if (cmd==name) {boundary=val;  //MM_ERR(MMPR(val)) 
} 
} // ii 
return boundary;
} // FieldInField

void GetBasics()
{
StrTy type;
StrTy boundary=field_in_field("content-type:","boundary",&type);
m_mime=type;
m_boundary=boundary;
add_key("boundary",m_boundary);
add_key("mime",m_mime);

//MM_ERR(MMPR3(m_mime,m_boundary,ct))
if (m_mime=="") {  
m_mime="text/plain/presumed";
//MM_ERR(m_verbatim)
//MM_ERR(" non mime write to nomime.xxx"<<MMPR2(m_verbatim.length(),m_verbatim))
//MM_ERR(" non mime write to nomime.xxx"<<MMPR3(depth(),account(), m_verbatim.length()))
// TODO FIXME check the RFC etc. 
//MM_ERR(" null mime assumed text/plain/presumed "<<MMPR3(depth(),account(), m_verbatim.length()))
if (false) //if (m_verbatim!="")
{
mjm_paranoid_ofstream<Tr>  sos(StrTy("nomime.xxx"));
//const StrTy v=m_formatter.convert_to_mbox(msg);
sos<<" non mime write to nomime.xxx"<<MMPR3(depth(),account(), m_verbatim.length()) ;
sos<<m_verbatim;
sos.finish();
} // not blank 
} // non mime
}


void MakeKeys(const StrTy & s)
{


Ragged&  x=m_parsed;
// read them into a ragged table do not remove quotes etc. 
// just make sure the first word is ok. 
x=SplitHeaders(s);
// concatenate tose continued with blanks 
x=MakeNormalHeaders(x);
// finally index for ease of use. 
x.to_map(m_keys,0);
/*
Ss ss;
ss<<m_hdr;
// this is making one big field no idea why 
//x.load_lines(ss);
// this needs to set options for quotes and lf/cr etc. 
x.load(ss,!true);
// headers can have repeated keys so this is not really right 

*/

}

// this needs to stop at top and let others parse the multipart hdr
Ragged SplitHeaders(const StrTy & h, const IdxTy flags=1)
{
Ragged dest;
const bool stop_at_blank=Bit(flags,0);
const bool remake_hdr=Bit(flags,1);
const bool keep_blank_hdr=Bit(flags,2);
StrTy nhdr;
Ss ss(h);
Ss sr;
Hand y;
StrTy d="";
StrTy from="";
StrTy to="";
StrTy subj="";
bool have_from=false;
//CommandInterpretter li(ss);
LineIterator li(&ss);
li.set_split(0,' ');
li.remove_crlf(true);
//IdxTy modified=0;
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
//if (m_flp.log_commands()) 
//              if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if( li.line().size()==0){
if (remake_hdr) m_hdr=nhdr;
else 
if (m_hdr.length()==0) if (!keep_blank_hdr) m_hdr=nhdr;
 if (stop_at_blank) break; 
}
if (remake_hdr) nhdr+=li.line()+"\r"; //// 
dest.add(li.line()); 
} // li.nextok()
// fill in with verbatim presumsed
if (m_body.length()==0)
{
Ss sr;
while (li.nextok()) { sr<<li.line()<<'\n'; }
m_body=sr.str();

}
return dest;
} // split_headers


// note that with repeated keys this messes up
// this is designed to concatenate mail headers 
// the crlf crap can go as it is cleaned up on parsing. 
// the new one needs to pickup the params of hinraw 
 Ragged MakeNormalHeaders(const Ragged & hinraw ) //const
{
Ragged dest;
const IdxTy nidx=0; // key;
const IdxTy vidx=1;
//const IdxTy minlen=((nidx>vidx)?nidx:vidx)+1;
StrTy last="";
typedef Ragged::Line Line;
Line x;
MM_SZ_LOOP(i,hinraw,hinsz)
{
const Line & line=hinraw.line(i);
//MM_ERR(MMPR2(line.size(),minlen))
if (line.size()<1) continue;
const StrTy&   nm=line[nidx];
const bool not_white=(nm.c_str()[0]!='\t');
const bool concat=(nm.length()!=0)&&not_white;
//if (nm.length()!=0) 
if (concat) 
	{if (x.size()!=0)  dest.add(x); x.clear(); } 
MM_SZ_LOOP(j,line,linsz)
{
//if (j!=0) 
x.push_back(line[j]);
} // j
/*
StrTy  v="";
StrTy sep="";
for ( IdxTy j=0; j<line.size(); ++j) if ( j!=key){if (line[j]!="\n") if (line[j]!="\r" )   v=v+sep+ line[j]; sep=StrTy(" "); }

const IdxTy sss=v.length();
char vx[v.length()+2];
IdxTy ptr=0;
for(IdxTy j=0; j<sss; ++j)
{
char c=v.c_str()[j];
if (c!='\r') if (c!='\n') { vx[ptr]=c; ++ptr;}
}
vx[ptr]=0;
if  ( nm.length()>0) rwm[nm]=vx;
else rwm[last]=rwm[last]+StrTy(" ")+vx;
last=nm;
//MM_ERR(" setting key "<<MMPR2(nm,v))
//rwm.set( nm,v,true);

*/

}

if (x.size()!=0) dest.add(x);  

return dest;
} // MakeNormalHeaders




}; //  mjm_generic_message




#endif // MJM_GENERIC_MESSAGE_H__ 
