
////////////////////////////////////////////
// hard coded and hard to change as static with limite
// params, no In thing... 
static IdxTy cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, OutTy & out, const IdxTy flags)
{
const bool always=true; // Bit(flags,0);
const bool file= Bit(flags,1);
const bool notrailingcrlf= Bit(flags,2); // 4 
const bool dfile= Bit(flags,3); /// 8
Blob d,err,cout;
StrTy _cmd;
//if (file) { _cmd="cat \""+s+"\" |"+cmd; }
if (file) { _cmd=cmd; d.load(s); }
else { d=s; _cmd=cmd; } 
IdxTy rc=out.hand().fileio(cout,err,d,_cmd,3);
if (always||(rc!=0)) { MM_ERR(MMPR3(dest,s,dfile)<< MMPR4(file,notrailingcrlf,flags,rc)<<MMPR4(cmd,StrTy(cout),StrTy(err),StrTy(d))) }
if (notrailingcrlf) {
IdxTy ncrlf= cout.no_trailing_crlf();
if (always) { MM_ERR(MMPR(ncrlf)) } 
 }
if (!dfile) dest=StrTy(cout);
else { cout.save(dest); } 
return rc;
}// cmd_exec

static IdxTy parse_hier_html(HierKvp & hkvp, Kvp * m, const StrTy & fntemp, const StrTy & fnin, const InpTy & in, OutTy & out, const IdxTy xflags )
{
const StrTy s1="testHTML --mjmflags hier -mjm \""+fnin+"\"";
MM_ERR(MMPR(s1))
StrTy fntemps=fntemp;
IdxTy rch=cmd_exec(fntemps,StrTy(),s1,out,8);
Ragged r;
r.kvp_splitter(0);
r.load(fntemp); // this should now be parsed with value in the end of line field 
MM_ERR(" loaded "<<MMPR2(fntemp,r.size()))
//HierKvp hkvp;
//Kvp m,me;
IdxTy lhrc=load_hier_obj(hkvp,m,r,in,out,0,xflags);
//Blob b(hkvp.dump());
return lhrc;
} // parse_hier_html 

static IdxTy filter_small_file(const StrTy & dfn, const StrTy & sfn, const StrTy & cmd, OutTy & out, const IdxTy flags)
{
const bool always=true; // Bit(flags,0);
//const bool file= Bit(flags,1);
//const bool notrailingcrlf= Bit(flags,2); // 4 
Blob d,err,cout;
StrTy _cmd="cat \""+sfn+"\" |"+cmd+" | tee \""+dfn+"\"";
IdxTy rc=out.hand().fileio(cout,err,d,_cmd,3);
if (always||(rc!=0)) { MM_ERR(MMPR2(flags,rc)<<MMPR4(cmd,StrTy(cout),StrTy(err),StrTy(d))) }
return rc;
} // filter_small_file

// could also use imagemagick identify 
static IdxTy identify(StrTy & dest, const StrTy & fn, const StrTy &_cmd, OutTy & out, const IdxTy flags)
{
const StrTy EXIF="exiftool -struct -f -x -ee  ";
const StrTy cmd=(_cmd.length())?_cmd:(EXIF+"\""+fn+"\" | grep \"File Type[^a-zA-Z]*:\" | tr \"[A-Z]\" \"[a-z]\" | awk '{printf \"%s\", $NF}'  ") ;
IdxTy rc= cmd_exec(dest,"",cmd,out,0);
// TODO FIXME update a cache in out
// out.file_info(fn,"identify",dest);

MM_ERR(" identifies as "<<MMPR3(dest,fn,cmd))
return rc;
}

//IdxTy rcx=get_xref(out,in,doi,nm,0);
// TODO this is redundant with guessdoi in the funcs header.. 
static  IdxTy get_xref(OutTy & out, const InpTy & in , const StrTy & doi, const StrTy & nm , const IdxTy flags )  { 
// this is a different counter, although neds to be the same
const IdxTy countz=out.id_count(doi,0);
// the caller is incrementing it first... doh... 
// this may be incorrectly reutrning "1" wtf... 
const IdxTy b4z=out.id_count("doi",doi,0);
MM_ERR(MMPR4(doi,countz,b4z,nm))
// in any case, gated out at out.id_count for now ... 
/*
this is fng up -all as later strategies fail to return results,
// should cache maybe... Note the count is global for "out"
 not sure if it descends... doph... 
./mjm_med2bib_handlers.h77  doi=10.1038/s41560-019-0404-4 countz=4 b4z=0 nm=handledoi
./mjm_med2bib_handlers.h81 returning since tried alrady  doi=10.1038/s41560-019-0404-4 countz=4 b4z=0
./mjm_med2bib_handlers.h554  breaking loop as later ones usually are refs or spurs
*/

if ((b4z+countz)>3 )
{
MM_ERR("should be returning since tried alrady "<<MMPR3(doi,countz,b4z))
if (false) return 0; 
}
const IdxTy count=out.id_count(doi,1);
const IdxTy b4=out.id_count("doi",doi,1);
MM_ERR(MMPR3(doi,count,b4))
const StrTy fnbib=out.fn("crossref",1); // in.ZZ
Grc grc=in.getter().crossref(fnbib,doi,0);
const StrTy src=grc.name();
const StrTy srcx=(src.length())?(StrTy("(")+src+")"):StrTy();
// should validate before adding
//parse(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) 
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grc.url(),nm+srcx,0);
return rcf;
} // get_xref
static StrTy isolatedomain()
{return " sed -e 's;\\([^/]\\)/[^/].*;\\1;'"; }

typedef std::map<StrTy, StrTy> KV;
KV kv;
static IdxTy remove_kv( KV & kv, const StrTy & key)
{
auto ii =kv.find(key);
if (ii==kv.end()) return 0;
kv.erase(ii);
return 1;
}
static IdxTy build_kv_map(KV & kv, const StrTy & href, const IdxTy flags)
{
typedef std::map<StrTy,StrTy> Lut;
Lut lut;
lut["author"]="authors";
lut["journal_title"]="journal";
lut["publication_date"]="date";
// buildentry
Ragged r;
Ss ss; ss<<href;
r.load(ss,false);
const IdxTy sz=r.size();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line l=r[i];
const IdxTy len=l.size();
if (len<1) continue;
StrTy key=l[0];
if (lut.find(key)!=lut.end()) key=lut[key];
StrTy& val=kv[key];
StrTy sep=" ";
if (key=="authors") //sep=" and ";
{
if (val.length()) { val=val+" and "; } 
for(IdxTy k=1; k<len; ++k) { if (val.length()) val=val+sep; val=val+l[k]; }}
else
{for(IdxTy k=1; k<len; ++k) { if (val.length()) val=val+sep; val=val+l[k]; }
}
} // i 
return 0;
}



//////////////////////////////////////////////////

 static  IdxTy handlecitmgr(const InpTy & in , OutTy & out , const IdxTy xflags=0 )  { 
const IdxTy cef=0;
const StrTy nm="handlecitmgr"; // +out.partial();;
out.enter(nm);

const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
// this will fail ????  if called from a domain guess. ZZ

//if (out.only(nm,in,"normalget"))  { 
//MM_ERR(" skipping "<<MMPR(nm)<<" due to only failure ") 
//out.exit(nm);return ~0; }


// only if missing and debug  
Grc grc=in.getter().normalget(fn,uin,16+32);
StrTy href="",sfx="",base="";
{
StrTy cmd="cat \""+fn+"\" | grep \"/citmgr.*gca=\"";
cmd_exec(href,"",cmd,out,cef);
if (href.length()==0) { out.exit(nm); return 1; } 
}
{
StrTy cmd=" sed -e 's;\\([^/]\\)/[^/].*;\\1;' "; 
cmd_exec(base,"",cmd,out,cef);
}
{
StrTy cmd=" sed -e 's;\\([^/]\\)/[^/].*;\\1;' "; 
cmd_exec(sfx,"",cmd,out,cef);
}
StrTy url=base+"/"+sfx+"&type-bibtex";
MM_ERR(MMPR3(nm,url,uin))
const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Grc grcb=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcb.url(),nm,0);
#if 0
hascitmgr=`cat "$fn" | grep "/citmgr.*gca=" `
if [ "$hascitmgr" != "" ]
then
guesscitmgr "$uin" "$fn" "$hascitmgr"


guesscitmgr()
{
uin="$1"
fn="$2"
href="$3"
base=`echo "$uin" | sed -e 's;\([^/]\)/[^/].*;\1;' `
sfx=`echo $href | sed  -e 's/.*citmgr/citmgr/' | sed -e 's/\(gca=[^"]*\)".*/\1/'`
url="$base/$sfx&type=bibtex"
# this should call the one that gets it 
lynx -source -dump "$url" > "$fn"

#endif


MM_ERR(" not impl"<<MMPR(nm))  
//MM_MSG(" not impl"<<MMPR(nm))  
out.exit(nm);
return 0; 
} // handlecitmgr


 static  IdxTy handlespring(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlespring"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

Grc grc=in.getter().normalget(fn,uin,16);
StrTy chan,url;
{StrTy cmd="cat \""+fn+"\" | grep \"citation.*format=bibtex&flavour=citation\"";
cmd_exec(chan,"",cmd,out,0);
}
if (chan.length()==0) { out.exit(nm);  return 1; } 
{StrTy cmd="sed -e 's/.*http:/http:/' | sed -e 's/[\"].*//'";
cmd_exec(url,chan,cmd,out,0);
}
const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Grc grcb=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcb.url(),nm,0);

#if 0
canhandle=`cat "$fn" | grep "citation.*format=bibtex&flavour=citation" `
if [ "$canhandle" != "" ]
then
url=`echo $canhandle | isolatelink`
echo found handlespring for $nm  something for $2  | smsg
echo handlespringer was trashing this | smsg
#$WGET --no-check-certificate --user-agent Firefox -q -O "$fn2" "$url"
normalget "$fn2" "$url"
#lynx -source -dump "$url" > "$fn" 
citepastecat  "$fn2" "$uin" "$url" "$2"

isolatelink()
{
sed -e 's/.*http:/http:/' | sed -e 's/["].*//'
}


#endif



MM_ERR(" not impl") 
out.exit(nm);
 return 0; }


// TODO may not be having a comma after doi line altough upper case URL
// fied is confusing the entry canonicallization but it should work  
 static  IdxTy handledowncit(const InpTy & in , OutTy & out , const IdxTy xflags=0)  {

const StrTy nm="handledowncit"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();

const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")){ out.exit(nm);  return ~0; } 

Grc grc=in.getter().normalget(fn,uin,16);
StrTy x,domain,doi;
{StrTy cmd="cat \""+fn+"\" | grep \"/action/showCitFormats\"";
cmd_exec(x,"",cmd,out,0);
}
if (x.length()==0){ out.exit(nm);  return 1; } 
{StrTy cmd="sed -e 's/.*[:]\\/\\/\\([^/#?]*\\).*/\\1/'"; // isolatedomain
cmd_exec(domain,uin,cmd,out,0);
}
if (domain.length()==0) { out.exit(nm);  return 2; } 
StrTy url=StrTy("http://")+domain+"/action/downloadCitation";
{StrTy cmd="urldecode | sed -e 's/.*\\(10\\.[0-9][0-9][0-9][0-9][0-9a-zA-Z.]*\\/[-0-9a-zA-Z.]*\\).*/\\1/'"; 
cmd_exec(doi,uin,cmd,out,0);
}
if (doi.length()==0) { out.exit(nm);  return 3; } 
// at this point calling crossref would work but
// pub may bebetter 
StrTy pd=StrTy("doi=")+doi+StrTy("&include=abs&format=bibtex&direct=&submit=Download");
const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Grc grcb=in.getter().normalpost(fnbib,url,pd,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcb.url(),nm,0);

#if 0
nm="handledowncit"
x=`cat "$fn" | grep "/action/showCitFormats"`

domain=`echo $uin | sed -e 's/.*[:]\/\/\([^/#?]*\).*/\1/'`
echo domain $domain | smsg
#url="http://atsjournals.org/action/downloadCitation"
url="http://$domain/action/downloadCitation"
doi=`echo $uin| urldecode | sed -e 's/.*\(10\.[0-9][0-9][0-9][0-9][0-9a-zA-Z.]*\/[-0-9a-zA-Z.]*\).*/\1/' `
pd="doi=$doi&include=abs&format=bibtex&direct=&submit=Download"

normalpost "$temp5" "$url" "$pd"
x=`cat "$temp5" | grep -i "@article" `
if [ "$x" != "" ]
then
citepastecat  "$temp5" "$uin" "$url $pd" 
#endif




// MM_ERR(" not impl") 

out.exit(nm);
 return 0; } // handledowncit




 static  IdxTy handlehighwire(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handledhighwire"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm) ;  return ~0; } 

Grc grc=in.getter().normalget(fn,uin,16);
StrTy chan,base,bibtext;
{StrTy cmd="cat \""+fn+"\" | sed -e 's/\"/\\n/g' | grep highwire | grep bibtex| sed -e 's/bibtext.*$/bibtext/'";
cmd_exec(chan,"",cmd,out,4);
}
if (chan.length()==0) { out.exit(nm);  return 1; } 
{StrTy cmd=" sed -e 's;\\([^/]\\)/[^/].*;\\1;'"; // isolatedomain
cmd_exec(base,uin,cmd,out,4);
}
if (base.length()==0) { out.exit(nm);  return 2; } 

StrTy url=base+StrTy("/")+chan;

const StrTy temp=out.fn("temp",1); // in.ZZ
const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ

Grc grct=in.getter().normalget(temp,url,0);

{StrTy cmd="cat \""+temp+"\" | awk '{if ($1==\"@article\") {x=\"\"; i=1; while (i<=NF){  x=x$i; i=i+1; } print x;       } else print $0;  }'";
//"cat \""+fn+"\" | sed -e 's/\"/\\n/g' | grep highwire | grep bibtex";
cmd_exec(bibtext,"",cmd,out,0);
}
if (bibtext.length()==0) { out.exit(nm);  return 3; } 
Blob b(bibtext);
b.save(fnbib);

IdxTy rcf= out.good_enough(b,fnbib,in,grct.url(),nm,0);


#if 0

handlehighwire()
{
fn="$1"
uin="$3"
canhandle=`cat  "$fn" | sed -e 's/"/\n/g' | grep highwire | grep bibtex`
#/highwire/citation/10397/bibtex
base=`echo $uin | isolatedomain`
url="$base/$canhandle"
#lynx -source -dump "$url" > "$fn"
$WGET --no-check-certificate -O "$temp5" -q "$url"
cat "$temp5" | awk '{if ($1=="@article") {x=""; i=1; while (i<=NF){  x=x$i; i=i+1; } print x;       } else print $0;  }' > "$temp6"

#endif





//MM_ERR(" not impl")  
out.exit(nm);
return rcf; 

} // handlehighwire

// seen in google books otherwise gets spuriosu hits 
 static  IdxTy handlebiblink(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlebiblink"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

Grc grc=in.getter().normalget(fn,uin,16);
StrTy biblink,base;
// bibtext should take out highwire hits 
{StrTy cmd="cat \""+fn+"\" | sed -e 's/\"/\\n/g' | grep bibtex | grep -v bibtext | head -n 1 | htmldecode"; 
cmd_exec(biblink,"",cmd,out,4);
}
if (biblink.length()==0) { out.exit(nm);  return 1; }
{StrTy cmd=" sed -e 's;\\([^/]\\)/[^/].*;\\1;'"; // isolatedomain
cmd_exec(base,uin,cmd,out,4); // no trailing crlf 
}
MM_ERR(MMPR3(nm,base,biblink))
if (base.length()==0) { out.exit(nm);  return 2; } 
StrTy url=base+StrTy("/")+biblink;

const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Grc grcb=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcb.url(),nm,0);


#if 0

nm=handlebiblink
echo trying $nm $2 | smsg
fn="$1"
uin="$3"
fn2="$temp5"
fn3="$temp6"
biblinks=`cat "$fn" | sed -e 's/"/\n/g' | grep bibtex| htmldecode`
echo found `echo $biblinks | wc ` | smsg
biblinks=`cat "$fn" | sed -e 's/"/\n/g' | grep bibtex| head -n 1 | htmldecode `
rm "$fn3"
echo trying $biblinks | smsg
for f in $biblinks ; do
## if this is a relative link 
xxx=`echo $f | grep http://`
if [ "$xxx" == "" ]
then
basex=`echo $uin | isolatedomain`
f="$basex/$f"
fi
# #echo found something for $2  | smsg
echo trying totl link $f | smsg

normalget "$fn2" "$f"


#endif
out.exit(nm); 

return 0; 

}
 static  IdxTy handledoilink(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
// mutate uin to get doi or related info  
const StrTy nm="handledoilink"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
//const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
typedef std::vector<StrTy> CmdVec;
CmdVec cv;
// first observed on citeseer
const StrTy termcrap="| sed -e 's/[&?].*//'";
//https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0014707
cv.push_back( "sed -e 's/.*=//' | grep \"^10\" | grep \"/\" "+termcrap); // doi=, plos
// https://onlinelibrary.wiley.com/doi/abs/10.1111/sjoe.12341
cv.push_back( "sed -e 's/.*doi//' | sed -e 's/[.=a-zA-Z;\\/]*10/10/' |  grep \"^10\" | grep \"/\" "+termcrap); // doi=, plos

// these 'tards not automation friendly get normalget right... 
//https://iopscience.iop.org/article/10.1088/0034-4885/79/8/084502
cv.push_back( "sed -e 's/.*article\\/10/10/'| sed -e 's/.*pdf\\/10/10/'  | grep \"^10\" |sed -e 's/[?&#].*//' |  grep \"/\" "+termcrap); // doi=, plos

//https://www.recercat.cat/bitstream/handle/2072/368596/10.1038%40s41560-019-0404-4.pdf?sequence=1

// 2023-01-03 lol took 3 urldecodes to get doi lol 
// https://scanr.enseignementsup-recherche.gouv.fr/publication/doi10.1088%25252f1742-5468%25252f2014%25252f09%25252fp09008
// const IdxTy poshttp=StrUtil::indexOf(l[lsz-1].c_str(),"http");
IdxTy decodes=0;
StrTy uind=uin;
while (decodes<5)
{
//const StrTy uind=MutateOnly(uin, "urldecode" ,out);
IdioRL  urlp(uind);
urlp.parse();
StrTy du=urlp.find_doi(0);
if (du.length())
{
Ss ss; ss<<nm+"(idiourl, "<<MMPR(decodes)<<")"; 
//IdxTy rcx=get_xref(out,in,du,nm+"(idiorl)" ,0);
IdxTy rcx=get_xref(out,in,du,ss.str(),0);
IdxTy nnow=out.found();
if (!in.collect_all()) if (nstart!=nnow) { out.exit(nm);  return 0; } 
} // du
 const IdxTy npct=StrUtil::indexOf(uind.c_str(),"%");
 const IdxTy nplus=StrUtil::indexOf(uind.c_str(),"+");
if (npct==~0) if (nplus==~0) break;
uind=MutateOnly(uind, "urldecode" ,out);
++decodes;
} // decodes



// should be handled in xref 
std::map<StrTy,IdxTy> seen;
MM_LOOP(ii,cv)
{
const StrTy& cmd=(*ii);
const StrTy doi=MutateOnly(uin, cmd,out);
MM_ERR(MMPR3(uin,doi,cmd))
if (doi.length()<5 ) continue;
//if (seen.find(doi)!=seen.end()) continue;
++seen[doi];
IdxTy rcx=get_xref(out,in,doi,nm,0);
IdxTy nnow=out.found();
if (!in.collect_all()) if (nstart!=nnow) break;
} // ii 

out.exit(nm); 
return 0; 
} // handledoilink
// TODO this is not really safe if the input is pdf, should use handlepdf 
// first or risk spurious results. 
 static  IdxTy handledoi(const InpTy & in , OutTy & out , const IdxTy xflags=0 )  { 
/*

Crossref returns unicode escapes, 
"given":"Fr\u00e9d\u00e9ric"
"affiliation":[{"name":"Institut de Chimie Mol\u00e9culaire et des mat\u00e9riaux d'Orsay (UMR 8182) Univ Paris Sud; Universit\u00e9 Paris Saclay; 91405 Orsay cedex France"}]}]
"author":"S\u00e9gaud"
"author":"Cant\u00fa\u2005Reinhard"
"author":"Cant\u00fa\u2005Reinhard"

That the json parser turns into this which appears to be ok utf-8
except the e' becomes A in latin-1 wth.. ,
./mjsonu.h312 GLOBAL 0  0 message 1 author 27  58 given = Fr\xc3\xa9d\xc3\xa9ric
./mjsonu.h312 GLOBAL 0  0 message 1 author 27  58 affiliation 59  60 name = Institut de Chimie Mol\xc3\xa9culaire et des mat\xc3\xa9riaux d'Orsay (UMR 8182) Univ Paris Sud; Universit\xc3\xa9 Paris Saclay; 91405 Orsay cedex France
./mjsonu.h312 GLOBAL 0  0 message 1 reference 64  81 author = S\xc3\xa9gaud
./mjsonu.h312 GLOBAL 0  0 message 1 reference 64  103 author = Cant\xc3\xba\xe2\x80\x85Reinhard
./mjsonu.h312 GLOBAL 0  0 message 1 reference 64  105 author = Cant\xc3\xba\xe2\x80\x85Reinhard



*/


const StrTy nm="handledoi"; // +out.partial();;
out.enter(nm);
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 
const IdxTy nstart=out.found();
Grc grc=in.getter().normalget(fn,in.uin(),16);
//std::vector<StrTy> v;
// if this is pdf, use pdftotext -layout FIRST lol. 
Scrapes v;
Src src=in.scraper().scrape_file(v,fn,1);
MM_ERR(" found "<<MMPR(v.size()))
MM_LOOP(ii,v)
{
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
const StrTy  doi=StrTy(*ii);
MM_ERR(" doi = "<<MMPR(StrTy(*ii)))
//MM_MSG(" doi = "<<MMPR2(nm,StrTy(*ii)))
if (StrTy(*ii).length()<5) continue;
// wtf, just let xref do it or allow disable.. 
//const IdxTy b4=out.id_count("doi",doi,1);
//const IdxTy b4=out.id_count("doi",doi,0);
//if (b4) continue; 
//IdxTy rcx=get_xref(out,in,StrTy(*ii),nm,0);
IdxTy rcx=get_xref(out,in,doi,nm,0);
MM_ERR(" breaking loop as later ones usually are refs or spurs")
//if (out.found()) break;
const IdxTy nnow=out.found();
//if (in.collect_all())  { MM_ERR(MMPR(nm)<<" ignoring collect all here") } 
if (!in.collect_all())  
if (nstart!=nnow) break;

} // ii 
out.exit(nm); 
return 0; 


} // handledoi

 static  IdxTy handlebibxml(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlebibxml"; // +out.partial();;
out.enter(nm);
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,in.uin(),16);
 
const StrTy TH="testHTML";
const StrTy THFLAGS=" --mjmflags 0 --mjm  ";

StrTy ispdf,fnp,link,base;
IdxTy rci= identify(ispdf, fn,"",  out, 0);
if (ispdf==StrTy("pdf"))
{
MM_ERR( MMPR(nm) << "skipping as "<<MMPR2(fn, ispdf)<<" is pdf ")
{ out.exit(nm);  return 1; } 
}
{StrTy cmd=TH+ THFLAGS+"\""+fn+"\" | grep -i \" href \" | grep -i bibtex";
cmd_exec(fnp,"",cmd,out,0);
}
if (fnp.length()==0)
{
MM_ERR( MMPR(nm) << "no links in  "<<MMPR(fn)<<" from xml parse ")
out.exit(nm); 
return 2;
}
{StrTy cmd="head -n 1 | awk '{print $NF}'";
//TH+ THFLAGS+"\""+fn+"\" | grep -i \" href \" | grep -i bibtex";
cmd_exec(link,fnp,cmd,out,0);
}
if (link.length()==0)
{
MM_ERR( MMPR(nm) << "no links in  "<<MMPR(fn)<<" from xml parse ")
out.exit(nm);
return 4;
}

{StrTy cmd=" sed -e 's;\\([^/]\\)/[^/].*;\\1;'"; // isolatedomain
cmd_exec(base,uin,cmd,out,0);
}
if (base.length()==0) { out.exit(nm);  return 5; } 
StrTy url=base+StrTy("/")+link;

const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Grc grcb=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcb.url(),nm,0);

#if 0

handlebibxml()
{
handler="$2" # "handlebibxml"
echo "  $handler fin is $1 " | smsg
ispdf=`$EXIF "$fn" | grep "File Type" | grep -i pdf `
if [ "$ispdf" != "" ]
then
echo skipping $handler due to $ispdf | smsg
else
$TH $THFLAGS "$fn" | grep -i " href " | grep -i bibtex  > "$fnp"
echo fnp is | smsg
cat "$fnp" | smsg
link=`cat "$fnp" | head -n 1 |  awk '{print $NF}'`
echo link is $link | smsg
if [ "$link" != "" ]
then
#testcrossref  "$fn" "$uin" "$doi" "$handler" > "$fn2"
base=`echo $uin | sed -e 's/\([^\/]\)\/\([^\/]\).*/\1/' `
url="$base$link"
echo full link is $url | smsg
normalget "$fn2" "$url"
#  ../toobib.out -cmd "parse bbb" -cmd "bib-cmd isc" -cmd "save bbb.bib 1" -quit
$TB -cmd "parse $fn2" -cmd "save $fnp 1" -quit
cat "$fnp" | grep -v "^%" > "$fn"
citepastecat  "$fn" "$uin" "$url" "$nm"  > "$fn2"
fi

#endif


//MM_ERR(" not impl")  

out.exit(nm); 
return rcf; 




}
// TODO see if this is redundant or otherwise can get arxiv or other ID info 
 static  IdxTy handledoixml(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="handledoixml"; // +out.partial();;
out.enter(nm);
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
const StrTy fntemp=out.fn("temp"); //
const IdxTy nstart=out.found();
//const StrTy TH="testHTML";
//const StrTy THFLAGS=" --mjmflags 0 --mjm  ";
StrTy ispdf,fnp,doi;
IdxTy rci= identify(ispdf, fn,"",  out, 0);
if (ispdf==StrTy("pdf"))
{
MM_ERR( MMPR(nm) << "skipping as "<<MMPR2(fn,ispdf)<<" is pdf ")
out.exit(nm); 
 return 1;
}
Blob b;
MM_ERR(" may be slow "<<MMPR2(nm,fn))
// this is obsolete, but should include eq sign.. 
Xrc xrc=in.xformer().html_to_parsed(b,fn,0);
//////////////////////////////////////////////////////

Ragged h,j;
// better but has change of "="
if (!true) 
{
typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
// really don't want to do this.. 
const bool only_json=true; 
const IdxTy fnty=256+8+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rc=bt.parse_both(h,j,fn,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))
} 


//////////////////////////////////////////////////////
b.save(fntemp);
MM_ERR(" may be slow too too "<<MMPR3(nm,fn,fntemp))
// $TH $THFLAGS "$fn" | grep -i "\( doi$\|[0-9] text 10\.[^ ]*$\)"  > "$fnp"
//{StrTy cmd=TH+ THFLAGS+"\""+fn+"\" |grep -i \"\\( doi$\\|[0-9] text 10\\.[^ ]*$\\)\"";
//cmd_exec(fnp,"",cmd,out,0);
//const StrTy cmd= " awk '{print $NF}' | grep \"/\" | grep \"^10\\.[^ ]*$\"";
MM_ERR(" trying html parser for doixml wtf")
//const StrTy cmd= " awk '{print $NF}' |tee xxx1 |  grep \"/\" | tee xxx2 |  grep \"\\(^10\\.[^ ]*$\\|^doi.*\\)\"| tee xxx3 |  sed -e 's/.*[:= ]//' | tee fick 4 |  sed -e 's/^[^1]*//'";
const StrTy cmd= " awk '{print $NF}' |  grep \"/\" |  grep \"\\(^10\\.[^ ]*$\\|^doi.*\\)\"| sed -e 's/.*[:= ]//' |  sed -e 's/^[^1]*//'";
fnp=out.fn("temp"); // fnp="AFCKSHTFCK";
cmd_exec(fnp,fntemp,cmd,out,2+8);
MM_ERR(" the file fnp  "<<MMPR(fnp))
Ragged r;
r.load(fnp);
const IdxTy ndois=r.size();
if (ndois==0)
{
MM_ERR( MMPR(nm) << "no links in  "<<MMPR(fnp)<<" from xml parse ")
out.exit(nm); 
return 2;
}
//doi=`cat "$fnp" | grep -A 1000 -i "doi"  | grep " 10\." | awk '{print $NF}'`

//{StrTy cmd= " grep -A 1000 -i \"doi\"  | grep \" 10\\.\" | head -n 1 | awk '{print $NF}'";
//cmd_exec(doi,fnp,cmd,out,0);
//}
//if (doi.length()==0)
//{
//MM_ERR( MMPR(nm) << "no doi  in  "<<MMPR(fn)<<" from xml parse ")
//out.exit(nm);
//return 3;
//}

// const IdxTy b4=out.id_count("doi",doi,1);
//std::map<StrTy, IdxTy> seen;
IdxTy goods=0;
MM_SZ_LOOP(i,r,szr)
{
const Ragged::Line & l=r[i];
if (l.size()!=1)
{
MM_ERR("pathological line "<<MMPR3(fnp,i,l.size()))
continue;
}
const StrTy & doi=l[0];
//if (seen[doi]!=0) { ++seen[doi]; continue ;} 
const IdxTy b4=out.id_count("doi",doi,0);
MM_ERR(MMPR2(doi,b4))
if (b4) continue; 
// only try one and even this one may be bad... 
if (goods) 
{
MM_ERR(" limiting "+nm+" to one hit "<<MMPR4(doi, b4,szr,i))
break;
}
// this prevents xref from trying it .... 
//out.id_count("doi",doi,1);
++goods;
Ss ss;
ss<<nm<<" ( quality"<<MMPR3(i,szr,goods)<<" )"; 
IdxTy rcx=get_xref(out,in,doi,ss.str(),0);
//IdxTy rcx=get_xref(out,in,doi,nm,0);
const IdxTy nnow=out.found();
if (!in.collect_all()) if (nnow!=nstart) break;
//++seen[doi];
} // i 

#if 0
const StrTy fnbib=out.fn("crossref",1); // in.ZZ
Grc grcc=in.getter().crossref(fnbib,doi,0);
// should validate before adding
//parse(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) 
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcc.url(),nm,0);
#endif

#if 0
handledoixml()
{
handler="$2" # "handledoixml"
echo "  $handler fin is $1 " | smsg
fn="$1"
uin="$3"
fn2="$temp2"
fnp="$temp3"
ispdf=`$EXIF "$fn" | grep "File Type" | grep -i pdf `
if [ "$ispdf" != "" ]
then
echo skipping $handler due to $ispdf | smsg
else
$TH $THFLAGS "$fn" | grep -i "\( doi$\|[0-9] text 10\.[^ ]*$\)"  > "$fnp"
echo fnp is | smsg
cat "$fnp" | smsg
doi=`cat "$fnp" | grep -A 1000 -i "doi"  | grep " 10\." | awk '{print $NF}'`
if [ "$doi" != "" ]
then
testcrossref  "$fn" "$uin" "$doi" "$handler" > "$fn2"
fi

#endif

out.exit(nm); 
//MM_ERR(" not impl")  
return 0; }

 static  IdxTy handlepdfexif(const InpTy & in , OutTy & out , const IdxTy xflags=0 )  { 
const StrTy nm="handlepdfexif"; // +out.partial();;
out.enter(nm);
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 
const IdxTy nstart=out.found();

//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,in.uin(),16);
const StrTy EXIF="exifutil -list ";
//xxxdoi=`exifutil -list "$fn" | grep -i "^doi[ \t]*:" | awk '{print $NF}'| head -n 1 `
StrTy xxxdoi,xxxbib,xxxisbn;
//{StrTy cmd=EXIF+"\""+fn+"\"  | grep -i \"^doi[ \t]*:\" | awk '{print $NF}'| head -n 1 ";
{StrTy cmd=EXIF+"\""+fn+"\"  | grep -i \"^doi[[:space:]]*:\" | awk '{print $NF}'| head -n 1 ";
cmd_exec(xxxdoi,"",cmd,out,0);
}
{StrTy cmd=EXIF+"\""+fn+"\"  | grep -i \"^isbn[[:space:]]*:\" | awk '{print $NF}'| head -n 1 ";
cmd_exec(xxxisbn,"",cmd,out,0);
}
if ( xxxisbn!="")
{
MM_ERR(" have an isbn "<<MMPR(xxxisbn))
}
const StrTy fnbib=out.fn("crossref",1); // in.ZZ
if (xxxdoi.length())
{
IdxTy rcf=get_xref(out,in,xxxdoi,nm,0);

#if 0
Grc grcc=in.getter().crossref(fnbib,xxxdoi,0);
// should validate before adding
//parse(Blob & _b, const StrTy & ofn, const StrTy & temp,const IdxTy flags) 
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grcc.url(),nm,0);
#endif

const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nstart!=nnow) { out.exit(nm);  return rcf; } 

//if (out.found()) { out.exit(nm);  return rcf; } 
} // xxxdoi.

//xxxbib=`exifutil -list "$fn" | grep -i "^x-bib-bibtex[ \t]*:" | awk '{ $1=""; $2="" ;  print $0}'| head -n 1 `

{StrTy cmd=EXIF+"\""+fn+"\" | grep -i \"^x-bib-bibtex[ \t]*:\" | awk '{ $1=\"\"; $2=\"\" ;  print $0}'| head -n 1 ";
cmd_exec(xxxbib,"",cmd,out,0);
}
if (xxxbib.length())
{
// echo "$xxxbib" | sed -e 's/^[ \t][ \t]*//g' | sed -e 's/\}[ \t]*, /\},\n/g'| sed -e 's/\}[ \t]*\}$/\}\n\}/g'  |  sed -e 's/^@\([^}]*\), /@\1,\n/' > "$fnfack"
{StrTy cmd=" sed -e 's/^[ \\t][ \\t]*//g' | sed -e 's/\\}[ \\t]*, /\\},\\n/g'| sed -e 's/\\}[ \\t]*\\}$/\\}\\n\\}/g'  |  sed -e 's/^@\\([^}]*\\), /@\\1,\\n/' ";
cmd_exec(xxxbib,xxxbib,cmd,out,0);
}

const StrTy fnbib2=out.fn("crossref",1); // in.ZZ
if (fnbib!=fnbib2)
{
MM_ERR(" bib names not matching probably lost... "<<MMPR2(fnbib,fnbib2))

}
Blob b(xxxbib);
b.save(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,uin,nm,0);
{
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nstart!=nnow) { out.exit(nm);  return rcf; } 
}

//if (out.found()){ out.exit(nm);  return rcf; } 
} // xxxbib

#if 0

handlepdfexif()
{
nm="handlepdfexif"
echo $nm | smsg
fn="$1"
uin="$3"
ispdf=`strings "$fn" | head -n 4 | grep -i "%PDF"`
if [ "$ispdf" == "" ]
then
echo probably not pdf, skipping `strings $fn | head -n 2 ` | smsg
else
fnfack="$temp5"
xxxdoi=`exifutil -list "$fn" | grep -i "^doi[ \t]*:" | awk '{print $NF}'| head -n 1 `
echo check doi $xxxdoi | smsg
checkdoi "$fnfack" "$uin" "$url" "$nm" "$xxxdoi"
if [ "$handled" == "" ]
then
xxxbib=`exifutil -list "$fn" | grep -i "^x-bib-bibtex[ \t]*:" | awk '{ $1=""; $2="" ;  print $0}'| head -n 1 `
if [ "$xxxbib" != "" ]
then
#echo "$xxxbib" | sed -e 's/^[ \t][ \t]*//g' | sed -e 's/\}[ \t]*, /\},\n/g'| sed -e 's/\}[ \t]*\}$/\}\n\}/g' | sed -e 's/^@\([^\{\}]*/), /@\1,\n/'  > "$fnfack"
#echo "$xxxbib" | sed -e 's/^[ \t][ \t]*//g' | sed -e 's/\}[ \t]*, /\},\n/g'| sed -e 's/\}[ \t]*\}$/\}\n\}/g'   > "$fnfack"
echo "$xxxbib" | sed -e 's/^[ \t][ \t]*//g' | sed -e 's/\}[ \t]*, /\},\n/g'| sed -e 's/\}[ \t]*\}$/\}\n\}/g'  |  sed -e 's/^@\([^}]*\), /@\1,\n/' > "$fnfack"
#cat "$fnfack"
citepastecat  "$fnfack" "$uin" "$url" "$nm" > "$temp6"
#echo  adfasfasfs  




#endif




///MM_ERR(" not impl")  
out.exit(nm); 
return 0; 

}
 static  IdxTy handlepdf(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlepdf"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
// first hit is from the pdf thing doh 
// http://su.diva-portal.org/smash/record.jsf?dswid=-8098&pid=diva2%3A189901&c=1&searchType=SIMPLE&language=en&query=189901&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all

if (true) if (type!="pdf")
{

MM_ERR(" type is not pdf so abort this now "<<MMPR4(uin,type,fn,tci))
out.exit(nm); 
return 0;

} // type != pdf 
// move this to the scraper... 
Scrapes v;
Src src=in.scraper().scrape_bin_file(v,fn,type,1);
MM_ERR(" found "<<MMPR(v.size()))
//std::map<StrTy,IdxTy> seen;
MM_LOOP(ii,v)
{
const StrTy doi=StrTy(*ii);
// xref should do this.. 
//const IdxTy b4=out.id_count("doi",doi,1);
const IdxTy b4=out.id_count("doi",doi,0);
MM_ERR(MMPR2(doi,b4))
if (b4) continue;
//auto jj=seen.find(doi);
//if (jj!=seen.end()) { ++seen[doi]; continue; } 
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
//MM_ERR(" doi = "<<MMPR(StrTy(*ii)))
MM_ERR(" doi = "<<MMPR(doi))
//if (StrTy(*ii).length()<5) continue;
if (doi.length()<5) continue;

//IdxTy rcf=get_xref(out,in,StrTy(*ii),nm,0);
IdxTy rcf=get_xref(out,in,doi,nm+" (binscrape)" ,0);

MM_ERR(" breaking loop as later ones usually are refs or spurs")
//if (out.found()) break;
{
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nstart!=nnow) { break; } 
}

} // ii 

const IdxTy nnow=out.found();
if ((nnow==nstart)  || in.collect_all())
//if (!out.found()) 
{
// convert to text and try handledoi... 
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
Scrapes v; IdxTy fsz=0;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) { MM_ERR(" pdftotext did not generate output ") }
else {
// there is a potential problem here 
Src srcp=in.scraper().scrape_bin_file(v,dest,"txt",1);
fsz=v.size();
Src srct=in.scraper().scrape_file(v,dest,1);
MM_ERR(" found "<<MMPR(v.size()))
}
// need a local as the glob was disabled for now... 
// just use the same one. 
std::map<StrTy, IdxTy> seen;
//MM_LOOP(ii,v)
MM_SZ_LOOP(i,v,vsz)
{
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
//const StrTy doi=StrTy(*ii);
const StrTy doi=StrTy(v[i]);
// this is  confusing lol. 
//MM_ERR(" doi = "<<MMPR(doi))
MM_ERR(MMPR(doi))
if (doi.length()<5) continue;
//const IdxTy b4=out.id_count("doi",doi,1);
const IdxTy b4=out.id_count("doi",doi,0);
MM_ERR(" doi count  "<<MMPR2(doi,b4))
if (b4) continue;
if (seen.find(doi)!=seen.end()) continue;
++seen[doi];
const StrTy sfx=" (pdftotext "+((i>=fsz)?StrTy("txt"):StrTy("bin"))+")";
//IdxTy rcx=get_xref(out,in,doi,nm+" (pdftotext)",0);
// now including arXiv id, 
if (strncmp(doi.c_str(),"arXiv:",6)==0)
{
IdxTy rt=Recurse(in,out,doi);
//MM_ERR(" need to recurse to arxiv 4 here ")
}
else {
IdxTy rcx=get_xref(out,in,doi,nm+sfx,0);
} // normal DOI 

const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nnow!=nstart)  break;

//if (out.found())  break;
} // ii 
} // not found

out.exit(nm); 
return 0;
 } // handlepdf

//////////////////////////////////////////////////////////////////////////////////
 static  IdxTy handlepdfstrings(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlepdfstrings"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
// move this to the scraper... 
Scrapes v;
Src src=in.scraper().scrape_bin_file(v,fn,type,1);
MM_ERR(" found "<<MMPR(v.size()))
//std::map<StrTy,IdxTy> seen;
MM_LOOP(ii,v)
{
const StrTy doi=StrTy(*ii);
const IdxTy b4=out.id_count("doi",doi,1);
if (b4) continue;
//auto jj=seen.find(doi);
//if (jj!=seen.end()) { ++seen[doi]; continue; } 
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
//MM_ERR(" doi = "<<MMPR(StrTy(*ii)))
MM_ERR(" doi = "<<MMPR(doi))
//if (StrTy(*ii).length()<5) continue;
if (doi.length()<5) continue;

//IdxTy rcf=get_xref(out,in,StrTy(*ii),nm,0);
IdxTy rcf=get_xref(out,in,doi,nm,0);

MM_ERR(" breaking loop as later ones usually are refs or spurs")
//if (out.found()) break;
{
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nstart!=nnow) { break; } 
}

} // ii 
#if 0
const IdxTy nnow=out.found();
if ((nnow==nstart)  || in.collect_all())
//if (!out.found()) 
{
// convert to text and try handledoi... 
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
Scrapes v;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) { MM_ERR(" pdftotext did not generate output ") }
else {
Src srcp=in.scraper().scrape_bin_file(v,dest,"txt",1);
Src srct=in.scraper().scrape_file(v,dest,1);
MM_ERR(" found "<<MMPR(v.size()))
}
// just use the same one. 
//std::map<StrTy, IdxTy> seen;
MM_LOOP(ii,v)
{
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
const StrTy doi=StrTy(*ii);
MM_ERR(" doi = "<<MMPR(doi))
if (doi.length()<5) continue;
const IdxTy b4=out.id_count("doi",doi,1);
if (b4) continue;
//if (seen.find(doi)!=seen.end()) continue;
//++seen[doi];
IdxTy rcx=get_xref(out,in,doi,nm+" (pdftotext)",0);
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nnow!=nstart)  break;

//if (out.found())  break;
} // ii 
} // not found

#endif

out.exit(nm); 
return 0;
 } // handlepdfstrings

 static  IdxTy handlepdftext(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlepdftext"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
// move this to the scraper... 
Scrapes v;
#if 0
Src src=in.scraper().scrape_bin_file(v,fn,type,1);
MM_ERR(" found "<<MMPR(v.size()))
//std::map<StrTy,IdxTy> seen;
MM_LOOP(ii,v)
{
const StrTy doi=StrTy(*ii);
const IdxTy b4=out.id_count("doi",doi,1);
if (b4) continue;
//auto jj=seen.find(doi);
//if (jj!=seen.end()) { ++seen[doi]; continue; } 
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
//MM_ERR(" doi = "<<MMPR(StrTy(*ii)))
MM_ERR(" doi = "<<MMPR(doi))
//if (StrTy(*ii).length()<5) continue;
if (doi.length()<5) continue;

//IdxTy rcf=get_xref(out,in,StrTy(*ii),nm,0);
IdxTy rcf=get_xref(out,in,doi,nm,0);

MM_ERR(" breaking loop as later ones usually are refs or spurs")
//if (out.found()) break;
{
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nstart!=nnow) { break; } 
}

} // ii 
#endif

const IdxTy nnow=out.found();
if ((nnow==nstart)  || in.collect_all())
//if (!out.found()) 
{
// convert to text and try handledoi... 
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
Scrapes v;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) { MM_ERR(" pdftotext did not generate output ") }
else {
Src srcp=in.scraper().scrape_bin_file(v,dest,"txt",1);
Src srct=in.scraper().scrape_file(v,dest,1);
MM_ERR(" found "<<MMPR(v.size()))
}
// just use the same one. 
//std::map<StrTy, IdxTy> seen;
MM_LOOP(ii,v)
{
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
const StrTy doi=StrTy(*ii);
MM_ERR(" doi = "<<MMPR(doi))
if (doi.length()<5) continue;
const IdxTy b4=out.id_count("doi",doi,1);
MM_ERR(MMPR2(doi,b4))
if (b4) continue;
//if (seen.find(doi)!=seen.end()) continue;
//++seen[doi];
IdxTy rcx=get_xref(out,in,doi,nm+" (pdftotext)",0);
const IdxTy nnow=out.found();
if (!in.collect_all()) 
if (nnow!=nstart)  break;

//if (out.found())  break;
} // ii 
} // not found

out.exit(nm); 
return 0;
 } // handlepdftext


/////////////////////////////////////////////////////

// scrape links from common disseration formats
 static  IdxTy handlediss(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlediss"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
if (type!="pdf") { MM_ERR(" skipping "<<MMPR3(nm,type,fn))  out.exit(nm); return 0;  } 
// move this to the scraper... 

const IdxTy nnow=out.found();
if ((nnow==nstart)  || in.collect_all())
// convert to text and try handledoi... 
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
Scrapes v;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) { 
MM_ERR(" pdftotext did not generate output ") 
out.exit(nm); 
return 0;
}
std::vector<StrTy> vv;
Blob b;
b.load(dest);
 StrTy s=MutateOnly(StrTy(b), "grep \"This paper is posted\" | sed -e 's/ /\\n/g' | grep \"^http\"  " ,out);
MM_ERR(MMPR2(nm,s))
if (s.length()) vv.push_back(s);

// https://dc.etsu.edu/cgi/viewcontent.cgi?article=3352&context=etd
/*
pdftotext -layout junk/xxx.pdf - | head -n 100
                                                                                                 East Tennessee State University
                                                                                         Digital Commons @ East

Follow this and additional works at: https://dc.etsu.edu/etd
     Part of the Analytical Chemistry Commons

Recommended Citation
Ghann, William Emmanuel, "Studies of Surfactants Effect on Riboflavin Fluorescence and Its Determination in Commercial Food
Products and Vitamin Tablets." (2008). Electronic Theses and Dissertations. Paper 2000. https://dc.etsu.edu/etd/2000

*/
s=MutateOnly(StrTy(b), "grep -A 8 \"Recommended Citation\" | sed -e 's/ /\\n/g' | grep \"^http\"  " ,out);
if (s.length()) vv.push_back(s);
MM_ERR(MMPR3(nm,s,vv.size()))



//else {
//Src srcp=in.scraper().scrape_bin_file(v,dest,"txt",1);
//Src srct=in.scraper().scrape_file(v,dest,1);
//MM_ERR(" found "<<MMPR(v.size()))
//}
// just use the same one. 
//std::map<StrTy, IdxTy> seen;
// doi has leading and probably trailing junk tpoo.. 
// plus a spurious one 
MM_LOOP(ii,vv)
{
const StrTy surl=StrTy(*ii);
MM_ERR(" surl = "<<MMPR2(surl,nm))
if (surl.length()<5) continue;
const IdxTy b4=out.id_count("url",surl,1);
MM_ERR(MMPR2(surl,b4))
if (b4) continue;
IdxTy rt=Recurse(in,out,surl);
//Blob b; b=be.format(); MM_ERR(" fixing kv map "<<MMPR(StrTy(b))) b.save(fnbib);
//IdxTy rcfge= out.good_enough(b,fnbib,in,uin,nm,0);
const IdxTy nnow=out.found();
if (!in.collect_all()) if (nnow!=nstart)  break;
} // ii 

out.exit(nm); 
return 0;
 } // handlediss





////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////

 static  IdxTy handlebomtex(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlebomtex"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 

const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
if (type=="pdf"){ out.exit(nm);  return ~0; } 
typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
////////////////////////////////////////////
// x.load(li.words(),1); 
Ragged r;
const IdxTy fnty=1; // html and json parse 
IdxTy rc=bt.load_ragged(r,fn,fnty);
MM_ERR(" loaded "<<MMPR(r.size()))
std::map<StrTy,StrTy> m;
IdxTy rca=bt.assemble(m,r,0);
StrTy s;
IdxTy rcf=bt.format(s,m,0);
MM_ERR(MMPR4(nm,fn,rca,rcf)<< MMPR4(s,m.size(),r.size(),rc))
if ((r.size()<3)||(m.size()<3))
{
// move this to the scraper... 
MM_ERR(" r or m too small ")
out.exit(nm); 
return 0;
}
const StrTy fnbib=out.fn("bibtemp",0);
BibEntry be;
be.add(m);
be.type("article");
be.name("missing");
Blob b;
b=be.format();
MM_ERR(" fixing kv map "<<MMPR(StrTy(b)))
b.save(fnbib);
IdxTy rcfge= out.good_enough(b,fnbib,in,uin,nm,0);
/////////////////////////////////////////////
// move this to the scraper... 
out.exit(nm); 
return 0;
 } // handlebomtex
static IdxTy handlexmlformats(const InpTy & in , OutTy & out , const IdxTy xflags=0)
{
// These blocked the download but can mutate the 
// url at least, Zotero works 
const StrTy nm="handlexmlformats";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const bool fiucking_sht = false;

//if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; }

//const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
//const StrTy & uin=in.uin();
// should use headless... 
//Grc grc=in.getter().normalget(fn,uin,16);

// needs "1" to dump the dom doh although could also work with
// getting the pdf which seems to occur by accident now... 
//Grc grc=in.getter().headlessget(fn,in.uin(),1+(show_trial?32:0));
Grc grc=in.getter().headlessget(fn,in.uin(),1+16);

// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
MM_ERR(MMPR3(nm,type,uin))
if (type=="pdf"){ out.exit(nm); return ~0; }


if (type!="xml"){  // out.exit(nm); return ~0; }
Blob b;
// thi uses the command line utility now... 
Xrc xrc=in.xformer().html_to_parsed_hier(b,fn,0);
Ragged r;
r.kvp_splitter(0);
Ss ss;
//  IS SHOULD CTROR WITH STRIGN
ss<<(StrTy(b));
r.load(ss,1); // 
typedef std::vector<StrTy> LinkVec;
LinkVec lvmets,lvmarc;
/*
1 html 102 body 105 div 108 div 351 div 506 div 654 div 659 div 662 div 665 div 668 div 671 div 2652 div 2655 div 2658 div 2671 div 2674 div 2677 p 2685 text = .
1 html 102 body 105 div 108 div 351 div 506 div 654 div 659 div 662 div 665 div 668 div 671 div 2652 div 2655 div 2658 div 2671 div 2674 div 2686 a 2687 href 2688 text = https://ufdcimages.uflib.ufl.edu/AA/00/00/49/56/00001/marc.xml
1 html 102 body 105 div 108 div 351 div 506 div 654 div 659 div 662 

marchywka@happy:/home/documents/cpp/proj/toobib$ cat xxx | grep "xml$"
1 html 102 body 105 div 108 div 351 div 506 div 654 div 659 div 662 div 665 div 668 div 671 div 2585 div 2588 div 2591 div 2604 div 2607 div 2642 a 2643 href 2644 text = https://ufdcimages.uflib.ufl.edu/AA/00/00/49/56/00001/AA00004956_00001.mets.xml
1 html 102 body 105 div 108 div 351 div 506 div 654 div 659 div 662 div 665 div 668 div 671 div 2652 div 2655 div 2658 div 2671 div 2674 div 2686 a 2687 href 2688 text = https://ufdcimages.uflib.ufl.edu/AA/00/00/49/56/00001/marc.xml



*/


MM_LOOP(ii,r)
{
const IdxTy sz=(*ii).size();
if (sz==0) continue;
// should look for href etc... 
const StrTy & s=(*ii)[sz-1];
const char * p=s.c_str();
const IdxTy len=s.length();
//ZZ
const bool mets= StrUtil::ends_in( "mets.xml",8,p,len);
const bool marc=StrUtil::ends_in( "marc.xml",8,p,len);
if (mets) lvmets.push_back(s);
if (marc) lvmarc.push_back(s);

} // ii 

MM_LOOP(ii,lvmets)
{

const StrTy fntemp=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
// should be all with wget
//Grc grc1=in.getter().normalget(fntemp,url,1|(1<<10));
// select UA but use wget this should not be adversarial...
Grc grc1=in.getter().normalget(fntemp,(*ii),1+32);
MM_ERR(MMPR3(nm,fntemp,(*ii)))
MetsX mx;
StrTy bib=mx.mets_xform(fntemp,0);
Blob b;
b=bib;
b.save(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,(*ii),nm+" ( mets )" ,0);

const IdxTy nnow=out.found();
MM_ERR(MMPR3(nm,nstart,nnow))
if (!all) if ( nstart<nnow) { out.exit(nm); return 0; }
//IdxTy rcr=Recurse(in,out,url);

} // ii 
} // xml
else
{

const StrTy fnbib=out.fn("bibtex",1); //
MetsX mx;
StrTy bib=mx.mets_xform(fn,0);
Blob b;
b=bib;
b.save(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,uin,nm+" ( xml mets file )",0);

const IdxTy nnow=out.found();
MM_ERR(MMPR4(nm,all,nstart,nnow))
if (!all) if ( nstart<nnow) { out.exit(nm); return 0; }

} // xml 

// this sort of works, not a bad garbage bin... 
if (type=="xml" ) if (!fiucking_sht)
{
const StrTy url=MutateOnly(uin, "sed -e 's/\\/doi\\/[^/][^/]*\\/10/\\/doi\\/xml\\/10/'  " ,out );
MM_ERR(MMPR3(nm,url,uin))
// if the url==uin, the isbib thing should pick it up... 
//if ((url!=uin) &&(url.length()!=0))
if ((url.length()!=0))
{ // try friendly ua
Grc grc1=in.getter().normalget(fntemp,url,1|(1<<10));
JatsX jx;
StrTy bib=jx.jats_xform(fntemp,0);
Blob b;
b=bib;
b.save(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm+" ( xml jats  )",0);

const IdxTy nnow=out.found();
if (!all) if ( nstart<nnow) { out.exit(nm); return 0; }
//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 


} // !fiucking_sht

out.exit(nm);
return 0;
} // handlexmlformats 


////////////////////////////////////////////////////////////////////

// grep oup.com ref/jn1260112746.txt  | sort | uniq -c | sort -g -r | head -n 1
// 6 Downloaded from https://academic.oup.com/jn/article/126/11/2746/4724647 by guest on 06 December 2021


 static  IdxTy handlesilverchair(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="handlesilverchair"; // +out.partial();
out.enter(nm); 
const StrTy fn=out.fn(); //
MM_ERR(MMPR(nm))
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 
MM_ERR(MMPR(nm))
const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// identify the thing,
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
MM_ERR(MMPR3(nm,fn,type))

//const IdxTy nnow=out.found();
//if ((nnow==nstart)  || in.collect_all())
//if (!out.found()) 
//{
// convert to text and try handledoi... 
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
//Scrapes v;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
MM_ERR(MMPR4(nm,type,dest,fn))
if (!in.getter().exists(dest)) { MM_ERR(" pdftotext did not generate output ") }
else {
StrTy dlink;
const StrTy cmd="grep \"Downloaded from \" | grep \" on \"  | sort | uniq -c | sort -g -r | head -n 1 | awk '{print $4}' " ;
// silverchair download link 
MM_ERR(MMPR4(nm,type,dest,cmd))
Xrc xrcc= in.xformer().text_to_link( dlink, dest, cmd, 0);
//{ return  FileToLine( d,  fn,  cmd ,  flags); }
MM_ERR(MMPR4(nm,type,dlink,cmd))
if (dlink.length()!=0)
{
IdxTy rt=Recurse(in,out,dlink);
} // dlink 
} // text worked 
//Src srcp=in.scraper().scrape_bin_file(v,dest,"txt",1);
//Src srct=in.scraper().scrape_file(v,dest,1);
//MM_ERR(" found "<<MMPR(v.size()))
// just use the same one. 
out.exit(nm); 
return 0;
 } // handlesilverchair







////////////////////////////////////////////////////////////////////////////////////////////
#if 0
// scan the binary
const StrTy sfx="| sed -e 's/[.;)]$//' ";
const StrTy c1_doix="`grep -a -i  \"doi\\|10\\..*/\" \""+fn+"\" | sed -e 's/.this version posted.*//' | grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;

const StrTy c2_doix="grep -a -i  \"first published \" \""+fn+"\" |grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\\/\\-]10/10/' "+sfx;

const StrTy c3_doix=" sed -e 's/  */\\\\n/g' | grep 10| head -n 1 | sed -e 's/http.*org.//' "+sfx;

#endif

// and also convert to text 
// strings will not likely work unless this is not a pdf file... 

#if 0
nm="handlepdf"
ispdf=`strings "$fn" | head -n 4 | grep -i "%PDF"`
pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if [ -z "$temp_pdf_fick" ]
then
echo  pdftotext failed $temp_pdf_fick `cat $temp_pdf_fick` | smsg
fi
# 2018-03-28 added -i for microbiome file 
# 2020-05-25 preprint.org adds garbage 
# 2020-07-22  sciencemag unlabelled  
doix=`grep -a -i  "doi\|10\..*/" "$fnpdf" | sed -e 's/.this version posted.*//' | grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\/\-]10/10/' `
if [ "$doix" == "" ]
then
doix=`grep -a -i  "first published " "$fnpdf" |grep 10 |  sed -e 's/.*[^.0-9a-zA-Z\/\-]10/10/' `
fi
# this could be multiple words, get complicated 2017-12-18
# the url may be better but for now instead of doi.org just use xref 
doix=`echo $doix | sed -e 's/  */\\n/g' | grep 10| head -n 1 | sed -e 's/http.*org.//' `
# 2018-01-01 trailing period from physiology.org pdf file 
doix=`echo $doix| sed -e 's/\.$//' `
doix=`echo $doix| sed -e 's/[;)]$//' `
echo DOIXAAAAAA extract $doix |  smsg
getcrossref "$nm" "$doix" "$fnout" "$uin"
#endif

// AccessTy html_to_parsed(Blob & d, const StrTy & fn, const IdxTy flags)
//{ return  FileToBlob( d,  fn,  "html_to_parsed" ,  flags); }

#if 0 
../testHTML --mjmflags 8 --mjm xxx  2>zzzz | grep "^ 0 html  1 head  2 meta"  | grep -A 3 citation_ | grep text | more
 0 html  1 head  2 meta  3 name 4 text citation_title
 0 html  1 head  2 meta  3 content 4 text Targeting Mitochondria as Therapeutic Stra
tegy for Metabolic Disorders
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Sorriento, Daniela
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Pascale, Antonietta Valeria
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Finelli, Rosa
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Carillo, Anna Lisa
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Annunziata, Roberto
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Trimarco, Bruno
 0 html  1 head  2 meta  3 name 4 text citation_author
 0 html  1 head  2 meta  3 content 4 text Iaccarino, Guido
 0 html  1 head  2 meta  3 name 4 text citation_journal_title
 0 html  1 head  2 meta  3 content 4 text The Scientific World Journal
 0 html  1 head  2 meta  3 name 4 text citation_publisher

#endif

// https://www.bibsonomy.org/bibtex/11b0628b13ff45ef42f959d678ce6d6c9/huiyangsfsu
/*
 2468  wget -O xxx -S -v "https://www.bibsonomy.org/bibtex/11b0628b13ff45ef42f959d678ce6d6c9/huiyangsfsu"
 2469  cat xxx | sed -e 's/@/\n@/g'  | sed -n '/^@.*{/,/^}$/p'
 2470  cp xxx xxx.html 
 2471  lynx -dump xxx.html  | sed -e 's/@/\n@/g'  | sed -n '/^@.*{/,/^}$/p'
 2472  lynx -dump xxx.html 
 2473  lynx -dump xxx.html  | more
 2474  lynx -dump xxx.html  | sed -e 's/@/\n@/g'  | sed -n '/^@.*{/,/^}_*$/p'
 2475  lynx -dump xxx.html  | more
 2476  lynx -dump xxx.html  | sed -e 's/@/\n@/g'  | sed -n '/^@.*{/,/^}_*$/p'
 2477  lynx -dump xxx.html  | sed -e 's/@/\n@/g'  | sed -n '/^@.*{/,/^}_*$/p'| grep 


*/

// https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8036
// https://bonndoc.ulb.uni-bonn.de/xmlui/bitstream/handle/20.500.11811/8036/5508.pdf;jsessionid=8EE2D4C371DD10838C62271A1D5CEC0E?sequence=1
 static  IdxTy handleisbib(const InpTy & in , OutTy & out, const IdxTy xflags=0)  
{


const StrTy nm="handleisbib"; // +out.partial();
MM_ERR(" enter "<<MMPR(nm))
out.enter(nm); 
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
if (out.only(nm,in,"normalget")) { out.exit(nm);  return ~0; } 
const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
// https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8036
StrTy cmd= "sed -e 's/@/\\n@/g'  | sed -n '/^@.*{/,/^}$/p'";
//StrTy cmd= "sed -e 's/@/\\n@/g'  | sed -n '/^@.*{/,/^}$/p'| grep -v \"^[^0-9a-zA-Z=,_{}]*urn:\"";
IdxTy rc=filter_small_file(fntemp, fn, cmd, out, 0);
//cmd="lynx --force_html -dump -width=1000 "+fntemp;
cmd="lynx --force_html -dump -width=1000 "+fntemp+" | grep -v \"^ *urn:\"";
IdxTy rc2=filter_small_file(fnbib, "", cmd, out, 0);
inserturl(fnbib,uin,out,0);
StrTy url=uin;
Blob b;
b.load(fnbib);

// 2021-09-13 returning junk try to find some fields 
const StrTy maybebib=MutateOnly(StrTy(b), "grep -i \"author\" " ,out);
const StrTy maybebib2=MutateOnly(StrTy(b), "grep -i \"title\" " ,out);
MM_ERR(MMPR2(maybebib,maybebib2))
if ( maybebib.length()) if (maybebib2.length())
{
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
} // lengths
MM_ERR(" exit "<<MMPR(nm))
// identify the thing,
out.exit(nm);
return 0; 
} // handleisbib
// in google books now... 
 static  IdxTy handlebibtexlink(const InpTy & in , OutTy & out, const IdxTy xflags=0)  
{
const StrTy nm="handlebibtexlink"; // +out.partial();
out.enter(nm); 
MM_ERR(" enter "<<MMPR(nm))
const bool all=in.collect_all();
const StrTy fn=out.fn(); //
const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
if (out.only(nm,in,"normalget")) { 
MM_ERR(" SHOULD BE aborting since file does not exist DANGER WILL ROBINSON "<<MMPR2(__FUNCTION__,in.uin()))
//out.exit(nm);  return ~0; 
} 
// only gets if missing... doh 
// this will fail as it needs the headless get to get the right
// parsing crap otherwise the link is buried in js coded crap 
//Grc grc=in.getter().normalget(fn,uin,16);
// dom is 1 doh  this works, NEEDS headless get...  
Grc grc=in.getter().headlessget(fn,uin,16+1);
BomTex bt;
////////////////////////////////////////////
// x.load(li.words(),1); 
Ragged r,h,j;
const IdxTy fnty=0; // 16*0+3; // html only  
// skip the json 
IdxTy rc=bt. parse_both(r,  j,  fn, 0);
//IdxTy rc=bt.load_ragged(r,fn,fnty);
MM_ERR(" loaded "<<MMPR2(fn,r.size()))
MM_ERR(" XXXXXXX "<<MMPR2(fn,r.dump()))
std::map<StrTy,StrTy> m;
IdxTy rca=0; // bt.assemble(m,r,0);
MM_LOOP(ii,r)
{
const auto & l=(*ii);
const IdxTy lsz=l.size();
if (lsz<2) continue;
// this MAY need moficiation of there is an "=" between k and v 
const IdxTy poshttp=StrUtil::indexOf(l[lsz-1].c_str(),"http");
const IdxTy posbib=StrUtil::indexOf(l[lsz-1].c_str(),"bibtex");
if (poshttp>10) continue;
if (posbib==~0) continue;
const StrTy nurl=l[lsz-1];
MM_ERR(" trying "<<MMPR(nurl))
InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);

const IdxTy nnow=out.found();
if(!all) if (nnow!=nstart) { out.exit(nm); return 0; }
//if (lsz>=7) if(l[lsz-7]=="page") continue;
} // ii


MM_ERR(" exit "<<MMPR(nm))
out.exit(nm);
return 0; 
} // handlebibtexlink 

 static  IdxTy handlehasbib(const InpTy & in , OutTy & out, const IdxTy xflags=0)  
{
// need to try 1) plain, 2) pdftotext, 3) html render
// also email and based64 decode url decode etc 
// could try to render but let headless do that on recurse
const StrTy nm="handlehasbib"; // +out.partial();
out.enter(nm); 
MM_ERR(" enter "<<MMPR(nm))
const bool all=in.collect_all();
const StrTy fn=out.fn(); //
std::vector<StrTy> fns;
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
if (out.only(nm,in,"normalget")) { 
MM_ERR(" SHOULD BE aborting since file does not exist DANGER WILL ROBINSON "<<MMPR2(__FUNCTION__,in.uin()))
//out.exit(nm);  return ~0; 

} 
const IdxTy nstart=out.found();
//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);

fns.push_back(fn);
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
MM_ERR(MMPR(type))
if (type=="pdf")
{
fns.push_back(out.fn("temp",1));
Xrc xrc=in.xformer().pdf_to_text(fns.back(),fn,0);
if (!in.getter().exists(fns.back())) { fns.pop_back();  MM_ERR(" pdftotext did not generate output ") }
} // pdf
else {
// let headless render for now just add html render 
MM_ERR(MMPR(fns.size()))
fns.push_back(out.fn("temp",1));
MM_ERR(MMPR(fns.size()))
Xrc xrc=in.xformer().html_to_rendered(fns.back(),fn,0);
MM_ERR(MMPR(fns.size()))
MM_ERR(MMPR(fns.back()))
if (!in.getter().exists(fns.back())) { fns.pop_back();  MM_ERR(" html renderer did not generate output ") }

} // output


HasBib hb;
typename HasBib::scan_result sr;
MM_ERR(MMPR(fns.size()))
MM_LOOP(ii,fns) { 
hb.scan(sr, (*ii),0); 
MM_ERR(MMPR3(nm,(*ii),sr.m_v.size())) 

 } 

std::vector<StrTy> fnb;
MM_ERR(MMPR(sr.m_v.size()))
MM_LOOP(ii,sr.m_v)
{
Blob b;
b.load(*ii);
const StrTy & url=uin;
IdxTy rcf= out.good_enough(b,(*ii),in,url,nm,0);
const IdxTy nnow=out.found();
MM_ERR(MMPR3(nstart,nnow,(*ii)))
if (!all) if (nnow!=nstart) break;
} // m_v

// https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8036
//StrTy cmd= "sed -e 's/@/\\n@/g'  | sed -n '/^@.*{/,/^}$/p'";
//StrTy cmd= "sed -e 's/@/\\n@/g'  | sed -n '/^@.*{/,/^}$/p'| grep -v \"^[^0-9a-zA-Z=,_{}]*urn:\"";
//IdxTy rc=filter_small_file(fntemp, fn, cmd, out, 0);
//cmd="lynx --force_html -dump -width=1000 "+fntemp;
//cmd="lynx --force_html -dump -width=1000 "+fntemp+" | grep -v \"^ *urn:\"";
//IdxTy rc2=filter_small_file(fnbib, "", cmd, out, 0);
//inserturl(fnbib,uin,out,0);
//StrTy url=uin;
//Blob b;
//b.load(fnbib);

// 2021-09-13 returning junk try to find some fields 
//const StrTy maybebib=MutateOnly(StrTy(b), "grep -i \"author\" " ,out);
//const StrTy maybebib2=MutateOnly(StrTy(b), "grep -i \"title\" " ,out);
//if ( maybebib.length()) if (maybebib2.length())
//{
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
//} // lengths

// identify the thing,
MM_ERR(" exit "<<MMPR(nm))
out.exit(nm);
return 0; 
} // handlehasbib



 static  IdxTy handlegsmeta(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

// skip if pdf... 
const IdxTy nstart=out.found();
IdxTy rc=handlegsmeta(in,out,xflags,2);
IdxTy nnow=out.found();
//if (out.found()) return rc;
if (!in.collect_all()) if (nnow!=nstart) return rc;
rc=handlegsmeta(in,out,xflags,0);
nnow=out.found();
//if (out.found()) return rc;
if (!in.collect_all()) if (nnow!=nstart) return rc;
return handlegsmeta(in,out,xflags,1);
}


 static  IdxTy handlegsmeta(const InpTy & in , OutTy & out, const IdxTy xflags, const IdxTy flags )  
{ 
const bool scraping_parser=((flags&3)==0);
const bool html_parser=((flags&3)==1);
const bool real_parser=((flags&3)==2);
MM_ERR(MMPR4(flags,scraping_parser,html_parser,real_parser))
const StrTy parser=(scraping_parser?"scraper":(real_parser?"real":"html"));
const StrTy nm=StrTy("handlegsmeta(")+parser+")"; // +out.partial();;
out.enter(nm);
const StrTy fn=out.fn(); //
if (out.only(nm,in,"normalget")){ out.exit(nm);  return ~0; } 

//const StrTy fn=out.fn("normalget"); // in.ZZ
const StrTy & uin=in.uin();
Grc grc=in.getter().normalget(fn,uin,16);
StrTy href;
if (real_parser)
{
 //static  IdxTy handlegsmetareal(StrTy & href,const StrTy & nm, const InpTy & in , OutTy & out, const StrTy & fn,  const IdxTy flags )  
IdxTy rcr=handlegsmetareal(href,nm,in,out,fn,flags);

}
if (scraping_parser)
{
const StrTy awk="awk '{ x=\"\"; if ( $2 != \"content=\") print \"% unusual syntax \" $0; if ( $NF != \"/>\" ) if ( $NF != \">\" )  print \"% not terminated as expected \"$0; i=3; while ( i < NF) {x=x\" \"$i; i=i+1; } print $1\" \"x};' ";

const StrTy nocrap="sed -e 's/^[ \\t][ \\t]*//g' | sed -e 's/  */ /g'| htmldecode ";
//const StrTy isolate=" tr \"[\\r\\n]\" \" \" | sed -e 's/<meta /\\n<meta /g' | sed -e 's/>/>\\n/g' | grep \"^<meta \" | sed -e 's/xmlns=[^ ]* / /g' | grep \"citation_\" | sed -e 's/^.meta[^\"]*\"//' | sed -e 's/\"/ /g'  | sed -e 's/  */ /g'";
const StrTy isolate=" tr \"[\\r\\n]\" \" \" | sed -e 's/<meta /\\n<meta /g' | sed -e 's/>/>\\n/g' | grep \"^<meta \" | sed -e 's/xmlns=[^ ]* / /g' |grep -v \"citation_reference\" |  grep \"citation_\" | sed -e 's/^.meta[^\"]*\"//' | sed -e 's/\"/ /g'  | sed -e 's/  */ /g'";

// 2021-06-04 failing to remove bepress 
//const StrTy cmd1=isolate+"|"+ awk+"|" +nocrap+" | sed -e 's/citation_//' " ;
const StrTy cmd1=isolate+"|"+ awk+"|" +nocrap+" | sed -e 's/.*citation_//' " ;
StrTy cmd=cmd1;
cmd_exec(href,fn,cmd,out,2);
} // scraping parser
if (html_parser)
{
IdxTy rch=handlehtmlparse(href,nm,in,out,fn,flags);
#if 0 
StrTy type;
IdxTy rci= identify(type, fn,"",  out, 0);
//MM_ERR( MMPR(nm) << "see skip  "<<MMPR2(fn, type)<<" is pdf ")
if (type==StrTy("pdf"))
{
MM_ERR( MMPR(nm) << "skipping as "<<MMPR2(fn, type)<<" is pdf ")
out.exit(nm); 
 return 1;
}
Blob b;
MM_ERR(" try the parsing gsmeta thing")
Xrc xrc=in.xformer().html_to_parsed(b,fn,0);
const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/citation_//'  ";
StrTy s=b;
MM_ERR(" gsmeta fick")
cmd_exec(href,s,cmd,out,0);
MM_ERR(" gsmeta results "<<MMPR2(nm,href))
//../testHTML --mjmflags 8 --mjm xxx  2>zzzz | grep "^ 0 html  1 head  2 meta"  | grep -A 3 citation_ | grep text | more
 //0 html  1 head  2 meta  3 content 4 text The Scientific World Journal
 //0 html  1 head  2 meta  3 name 4 text citation_publisher

#endif

} // html_parser


if (href.length()==0) { out.exit(nm); return 1; } 

IdxTy rcfm=FixMeta(href, in,out,nm,uin);

out.exit(nm);
return 0; 


}
 static  IdxTy handlehtmlparse(StrTy & href,const StrTy & nm, const InpTy & in , OutTy & out, const StrTy & fn,  const IdxTy flags )  
{
StrTy type;
IdxTy rci= identify(type, fn,"",  out, 0);
//MM_ERR( MMPR(nm) << "see skip  "<<MMPR2(fn, type)<<" is pdf ")
if (type==StrTy("pdf"))
{
MM_ERR( MMPR(nm) << "skipping as "<<MMPR2(fn, type)<<" is pdf ")
//out.exit(nm); 
 return 1;
}
Blob b;
MM_ERR(" try the parsing gsmeta thing")
// this is failing now on clever html,
// https://ecommons.cornell.edu/handle/1813/41073
// but this is generally obserolete as the more general code 
// should get first shot at it... 
Xrc xrc=in.xformer().html_to_parsed(b,fn,0);
// although it does come up, 
// should be unambiguous now, 
//Xrc xrc=in.xformer().html_to_parsed_hier(b,fn,0);
/*

1 html 8 head 297 meta 300 name 301 text = citation_title
1 html 8 head 302 meta 303 content 304 text = Georgia Institute of Technology
1 html 8 head 302 meta 305 name 306 text = citation_publisher
1 html 8 head 307 meta 308 content 309 text = en_US
1 html 8 head 307 meta 310 name 311 text = citation_language
1 html 8 head 312 meta 313 content 314 text = Tillotson, John P.
1 html 8 head 312 meta 315 name 316 text = citation_author
1 html 8 head 317 meta 318 content 319 text = https://smartech.gatech.edu/bitstream/1853/60808/1/TILLOTSON-DISSERTATION-2018.pdf
1 html 8 head 317 meta 320 name 321 text = citation_pdf_url
1 html 8 head 322 meta 323 content 324 text = 2018-11-15
1 html 8 head 322 meta 325 name 326 text = citation_date
1 html 8 head 327 meta 328 content 329 text = https://smartech.gatech.edu/handle/1853/60808
1 html 8 head 327 meta 330 name 331 text = citation_abstract_html_url
1 html 332 body 333 text =


*/


// failing on bepress citations 2021-06-04
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/citation_//'  ";
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/.*citation_//'  ";
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/.*citation_//'  ";
const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee med2bib_test_file_sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | grep -v \"citation_reference\" | sed  -e 's/.*citation_//'  ";
StrTy s=b;
MM_ERR(" gsmeta huge problem ")
cmd_exec(href,s,cmd,out,0);
MM_ERR(" gsmeta results "<<MMPR2(nm,href))
//../testHTML --mjmflags 8 --mjm xxx  2>zzzz | grep "^ 0 html  1 head  2 meta"  | grep -A 3 citation_ | grep text | more
 //0 html  1 head  2 meta  3 content 4 text The Scientific World Journal
 //0 html  1 head  2 meta  3 name 4 text citation_publisher
return 0; 
} // handlehtmlparse

 static  IdxTy handlegsmetareal(StrTy & href,const StrTy & nm, const InpTy & in , OutTy & out, const StrTy & fn,  const IdxTy flags )  
{
StrTy type;
IdxTy rci= identify(type, fn,"",  out, 0);
//MM_ERR( MMPR(nm) << "see skip  "<<MMPR2(fn, type)<<" is pdf ")
if (type==StrTy("pdf"))
{
MM_ERR( MMPR(nm) << "skipping as "<<MMPR2(fn, type)<<" is pdf ")
//out.exit(nm); 
 return 1;
}
Blob b;
MM_ERR(" try the parsing gsmeta new real parser thing")
// but this is generally obserolete as the more general code 
//Xrc xrc=in.xformer().html_to_parsed(b,fn,0);
// should be unambiguous now, 
Xrc xrc=in.xformer().html_to_parsed_hier(b,fn,0);
HierUtil hu;
//typedef std::map<StrTy, std::map<StrTy, StrTy> > Hm;
typedef typename HierUtil::hier_map Hm;
Hm hm;
hu.organize(hm,StrTy(b),"meta",0);
/*
Ragged r;
r.kvp_splitter(0);
Ss ss; ss<<StrTy(b);
r.load(ss,false);
MM_ERR(MMPR2(nm,r.size()))
MM_SZ_LOOP(i,r,szr)
{
const Ragged::Line & l=r[i];
const IdxTy len=l.size();
if(len<5) {continue; }
for (IdxTy j=3; j<(len-3); ++j)
{
//MM_ERR(MMPR2(nm,len))
if ( l[j]=="meta")
{
//MM_ERR(MMPR2(nm,l[j]))
hm[l[j-1]][l[len-4]]=l[len-1];
MM_ERR(MMPR4(nm,l[j],l[j-1],l[len-4])<<MMPR(l[len-1]))

} // meta 

} // j 
} // i
*/

Ss rr;
MM_LOOP(ii,hm)
{
//const auto jj=(*ii);
const auto & key=(*ii).first;
StrTy k="";
StrTy v="";
MM_LOOP(kk,(*ii).second)
{
const StrTy & i1=(*kk).first;
const StrTy & i2=(*kk).second;
if (i1=="name") k+=i2;
else if (i1=="content") v+=i2;
else { MM_ERR(" unused meta "<<MMPR3(key,i1,i2))}
} // kk 
if (k!="") rr<<k<<" "<<v<<CRLF;
else {  MM_ERR(" no meta key "<<MMPR2(key,v)) } 
} // hm
href=rr.str();


/*
 equals is deleted in parse,

1 html 8 head 297 meta 300 name 301 text = citation_title
1 html 8 head 302 meta 303 content 304 text = Georgia Institute of Technology
1 html 8 head 302 meta 305 name 306 text = citation_publisher
1 html 8 head 307 meta 308 content 309 text = en_US
1 html 8 head 307 meta 310 name 311 text = citation_language
1 html 8 head 312 meta 313 content 314 text = Tillotson, John P.
1 html 8 head 312 meta 315 name 316 text = citation_author
1 html 8 head 317 meta 318 content 319 text = https://smartech.gatech.edu/bitstream/1853/60808/1/TILLOTSON-DISSERTATION-2018.pdf
1 html 8 head 317 meta 320 name 321 text = citation_pdf_url
1 html 8 head 322 meta 323 content 324 text = 2018-11-15
1 html 8 head 322 meta 325 name 326 text = citation_date
1 html 8 head 327 meta 328 content 329 text = https://smartech.gatech.edu/handle/1853/60808
1 html 8 head 327 meta 330 name 331 text = citation_abstract_html_url
1 html 332 body 333 text =


*/


// failing on bepress citations 2021-06-04
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/citation_//'  ";
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/.*citation_//'  ";
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | sed  -e 's/.*citation_//'  ";
//const StrTy cmd="grep \" 0 html  1 head\" | grep meta | grep -A 5 citation_ | grep text| tee sht1  | awk '{ i=1; while (i<8) { $i=\"\"; ++i; } if ($8==\"content\") { $8=\"\"; $9=\"\"; $10=\"\";  print name\" \"$0 ; name=\"\" } else if ($8==\"name\") {$8=\"\"; $9=\"\"; $10=\"\"; name=$0}}' | sed -e 's/^  *//' | sed -e 's/  */ /g'| htmldecode | grep -v \"citation_reference\" | sed  -e 's/.*citation_//'  ";
//StrTy s=b;
MM_ERR(" gsmeta real  huge problem ")
//cmd_exec(href,s,cmd,out,0);
MM_ERR(" gsmeta results "<<MMPR2(nm,href))
//../testHTML --mjmflags 8 --mjm xxx  2>zzzz | grep "^ 0 html  1 head  2 meta"  | grep -A 3 citation_ | grep text | more
 //0 html  1 head  2 meta  3 content 4 text The Scientific World Journal
 //0 html  1 head  2 meta  3 name 4 text citation_publisher
return 0; 
} // handlehtmlparse





// modify entries and save the bibtex if any... 
static IdxTy FixMeta( const StrTy & href,const InpTy & in, OutTy & out, const StrTy & nm, const StrTy & uin )
{
//typedef std::map<StrTy, StrTy> KV;
KV kv;
//static IdxTy build_kv_map(KV & kv, const StrTy & href, const IdxTy flags)
 IdxTy brc= build_kv_map(kv,  href, 0);
#if 0 
// buildentry
Ragged r;
Ss ss; ss<<href;
r.load(ss,false);
const IdxTy sz=r.size();
for(IdxTy i=0; i<sz; ++i)
{
const Ragged::Line l=r[i];
const IdxTy len=l.size();
if (len<1) continue;
StrTy key=kv[l[0]];
if (lut.find(key)!=lut.end()) key=lut[key];
StrTy& val=kv[l[0]];
const StrTy sep=" ";
for(IdxTy k=1; k<len; ++k) { if (val.length()) val=val+sep; val=val+l[k]; }

} // i 
#endif

const StrTy fp= kv["firstpage"];
const StrTy lp= kv["lastpage"];
if (fp!="") if (lp!="")
{
kv["pages"]=fp+"-"+lp;
kv.erase(kv.find("firstpage"));
kv.erase(kv.find("lastpage"));
}
remove_kv(kv,"reference"); // huge list of crap 
//auto ii=kv.find("publication_date");
auto ii=kv.find("date");
if (ii!=kv.end())
{
const StrTy  v=(*ii).second;
if (v.length()==10)
{
if ((v.c_str()[4]=='/') ||(v.c_str()[4]=='-'))
{
kv["year"]=StrTy(v.c_str(),0,4);
kv["month"]=StrTy(v.c_str(),5,2);
kv["day"]=StrTy(v.c_str(),8,2);
//remove_kv(kv,"publication_date");
}
} //10 
} // pub_Date

// kv map needs to be translated, need to add a name and type
// and comments
BibEntry be;
be.add(kv);
be.make_name();
if (be.name()=="") be.name("missing");
be.type("article");
const StrTy fnbib=out.fn("bibtexf",1); // in.ZZ
Blob b;
b=be.format();
MM_ERR(" fixing kv map "<<MMPR(StrTy(b)))
b.save(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,uin,nm,0);
return 0;
}
#if 0
awk_cits()
{
# this needs to send comments elsewhere it messes p otherwise ok line 
awk '{ x=""; if ( $2 != "content=") print "% unusual syntax " $0; if ( $NF != "/>" ) if ( $NF != ">" )  print "% not terminated as expected "$0; i=3; while ( i < NF) {x=x" "$i; i=i+1; } print $1" "x};'
}
nocrap()
{
sed -e 's/^[ \t][ \t]*//g' | sed -e 's/  */ /g'
}
isolate_cits()
{
# 2020-12-18  need  xml parser but just remove crlf 
# for now 
tr "[\r\n]" " " | sed -e 's/<meta /\n<meta /g' | sed -e 's/>/>\n/g' | grep "^<meta " | sed -e 's/xmlns=[^ ]* / /g' | grep "citation_" | sed -e 's/^.meta[^"]*"//' | sed -e 's/"/ /g'  | sed -e 's/  */ /g'




 more junk/bi_parse  | grep  -A 4 application/ld+ | grep "{" | sed -e 's/[^{]*//' | sed -e 's/\\\"/%27/g' | urldecode | while read ; do echo $REPLY | mjsonu.out -cmd "sax-json" ; done
#endif
// this is now obsolete. see json2
 static  IdxTy handleldjson(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="handleldjson";
out.enter(nm); 
///////////////////////////////////////////
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntemph=out.fn("temp",1); //
StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// from cmd line, the profiles , adding 8 to flags, seem to be a problem  
// adding 16 to only get once. 
if (false) //  (headlessunsafe(in,out,xflags)) 
{
MM_ERR(" headless not safe using the wget version "<<MMPR2(fnme,fn))
fnme=fn; 
}
else
{
// this can hang. known issue need puppeteer and newer nodejs  
Grc grc=in.getter().meget(fnme,uin,1+8+16);
}
MM_ERR(" the stupid  things shouldbe in messed upi  "<<MMPR(fnme))
//marchywka@happy:/home/documents/cpp/proj/toobib$ cat zzz | grep null | grep context | sed -e 's/[^{]*//' | mjsonu.out -cmd "sax-json"
// application/ld+json
MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))
const StrTy s1="testHTML --mjmflags 0 -mjm \""+fnme+"\"";
// isolate the json pieces 
//const StrTy s2="cat \""+fntemph+"\" | grep -A 4  \"application/ld+json\" | tail -n +1 |  grep null | grep \"{\" | grep \"@type\" |sed -e 's/[^{]*//' ";
// this runs a risk of getting the real one in the wrong block 
//const StrTy s2="cat \""+fntemph+"\" | grep -A 8  \"application/ld+json\" | tail -n +1 |  grep null | grep \"{\" | grep \"@type\" |sed -e 's/[^{]*//' ";
// this fails on the guardian for examled
//const StrTy s2="cat \""+fntemph+"\" | grep -A 8  \"application/ld+json\" | tail -n +1 |  grep null | grep \"null) {\" | grep \"@type\" |sed -e 's/[^{]*//' ";
const StrTy s2="cat \""+fntemph+"\" | grep -A 8  \"application/ld+json\" | tail -n +1 |  grep null | grep \"null) .*{\" | grep \"@type\" |sed -e 's/[^{]*//' ";


const StrTy s3="i=0; cat \""+fntempg+"\" | while read ; do echo $REPLY | mjsonu.out -cmd \"sax-json\" | awk -vf=\"$i\" '{print f\" \"$0 ;}'  ;i=$[ $i + 1 ] ;  done  ";
//const StrTy sc=s1+" | "+s2+" | "+s3;
//const StrTy sc=s2+" | "+s3;
// blank file names a problem... 
MM_ERR(" start ")
// parse html into fntemph 
IdxTy rch=cmd_exec(fntemph,StrTy(),s1,out,8);
//MM_MSG(" should have parsed html in "<<MMPR(fntemph))
IdxTy rcg=cmd_exec(fntempg,StrTy(),s2,out,8);
//MM_MSG(" should have json input  in "<<MMPR(fntempg))
Blob bio;
bio.load(fntempg);
//MM_MSG(MMPR(StrTy(bio)))

// need to parse the lines separately...
Ragged r;
r.load_lines(fntempg);
const IdxTy jsz=r.size();
MM_ERR(MMPR(jsz))
// NB in theory,any external size needs to be checked for
// limits of stack or even heap although stack can crash
// heap should return null or something 
BibEntry be[jsz];
for(IdxTy i=0; i<r.size(); ++i)
{
const Ragged::Line & ll=r[i];
if (ll.size()!=1)
{ MM_ERR(" bad line "<<MMPR2(i,ll.size()))
continue;
}
if (ll[1].size()==0)
{
MM_ERR(" there is no json probably writing js junk  out ")
continue; 

}
// this may not be a bug depending on validity of backslash quote
// the problem was not checking return code now it plows through
// but may still treat quote as invalid. 
// there is a bug in rapidjson with backslsh quote that infinite loops
// this is a kluge but may mess up other things.. 
 //const StrTy s4=" sed -e 's/\\\"/%27/g' | urldecode | mjsonu.out -cmd "sax-json" "; 
 //const StrTy s4=" sed -e 's/\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json\"  | sed -e 's/mjm_bs_bug/\\\"/'  "; 
 const StrTy s4=" sed -e 's/\\\\\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json\"  | sed -e 's/mjm_bs_bug/\\\\\\\"/g'  "; 
// use the easier to partse outut 
const bool new_format=true;
 const StrTy s4new=" sed -e 's/\\\\\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json - 1\"  | sed -e 's/mjm_bs_bug/\\\\\\\"/g'  "; 


//MM_MSG(" parsing "<<  MMPR(ll[0]))


IdxTy rcj=cmd_exec(fntempj,ll[0],(new_format?s4new:s4),out,8);
bio.clear();
bio.load(fntempj);
//MM_MSG(MMPR(StrTy(bio)))

const IdxTy bflags=0;
const StrTy fnbibnew=out.fn("bibtex",1); //
IdxTy rcs=new_format?sep_app_json_new(fnbibnew,nm,in,fntempj,out,bflags,xflags)
 : sep_app_json(fnbibnew,nm,in,fntempj,out,bflags,xflags);
//fntempj should now have the parsed json... 


// this parses the FIL not the blob 
//if (clean) in.finisher().light_clean(fnbib,fn("temp"),0);
Blob b;
// These will not parse because of special cars and maybe blank
// name. But files are ok.  
Frc frc=in.finisher().parse(b,fnbibnew,"",0);
//Frc frcb=0; // in.finisher().validate(b,fnbib,"",0);
BibEntry & ber=frc.be();
be[i]=ber;
//MM_MSG(MMPR3(i,fnbibnew,be[i].format()))

#if 0
 i=10; more junk/bi_parse  | grep  -A 4 application/ld+ | grep "{" | sed -e 's/[^{]*//' | sed -e 's/\\\"/%27/g' | urldecode | while read ; do echo $REPLY | mjsonu.out -cmd "sax-json"| awk -vf=$i '{print f" "$0}' ; i=$[ $i+ 1 ]  ; done

#endif


} // i 

///////////////////////////////////////////
out.exit(nm); 
 return 0;

} // handleldjson

///////////////////////////////////////////////////////////
// thi8s obsoletes the json version 
typedef typename HierKvp::candidate_nodes_type KvNodes;
typedef  std::map<StrTy, std::vector<StrTy>> Kvvp;
typedef std::vector<Kvp>  KvpVec; //  mj,mej;
typedef std::vector<Ragged>  RagVec; //  mj,mej;
#if 0 
 static  IdxTy handleldjson2Old(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="handleldjson2Old";
out.enter(nm); 
///////////////////////////////////////////
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntemph=out.fn("temp",1); //
StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const bool combine=true;
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
KvpVec kvps;
RagVec rags;
// get the dom 
// from cmd line, the profiles , adding 8 to flags, seem to be a problem  
// adding 16 to only get once. 
const StrTy maybepdf=MutateOnly(uin, "grep -i \"\\.pdf$\"" ,out);
// should try to get with wget and if zero len comes back, wait try with chrome 
const IdxTy insz=in.getter().size(fn);
StrTy ispdf;
//static IdxTy identify(StrTy & dest, const StrTy & fn, const StrTy &_cmd, OutTy & out, const IdxTy flags)
IdxTy idtc=identify(ispdf,fn,"",out,0);
//if (false) //  (headlessunsafe(in,out,xflags)) 
if ((maybepdf.length()!=0) ||( ispdf=="pdf" ))
{
MM_ERR(" headless not safe using the wget version "<<MMPR2(fnme,fn))
fnme=fn; 
}
else
{
// this can hang. known issue need puppeteer and newer nodejs  
// TODO use Chromate doh... 
// may already have it... 
MM_ERR(" NOT using MeGt bad idaea "<<MMPR(uin))
//Grc grc=in.getter().meget(fnme,uin,1+8+16);
}
MM_ERR(" the stupid  things shouldbe in messed up  "<<MMPR(fnme))
//marchywka@happy:/home/documents/cpp/proj/toobib$ cat zzz | grep null | grep context | sed -e 's/[^{]*//' | mjsonu.out -cmd "sax-json"
// application/ld+json
MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))
//////////////////////////////////////////////
// may want to use me...
StrTy fnin=fnme;
// only if missing... 
Grc grc=in.getter().normalget(fnin,uin,16+32);
const StrTy s1="testHTML --mjmflags hier -mjm \""+fnin+"\"";
MM_ERR(MMPR(s1))
IdxTy rch=cmd_exec(fntemph,StrTy(),s1,out,8);
Ragged r;
r.kvp_splitter(0);
r.load(fntemph); // this should now be parsed with value in the end of line field 

/////////////////////////////////////////////
// better but has change of "="
if (false) 
{
typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
Ragged r;
const IdxTy fnty=1; // html and json parse 
IdxTy rc=bt.load_ragged(r,fn,fnty);
MM_ERR(" loaded "<<MMPR(r.size()))
} 
//////////////////////////////////////////////



MM_ERR(" loaded "<<MMPR2(fntemph,r.size()))
HierKvp hkvp;
Kvp m,me;
IdxTy lhrc=load_hier_obj(hkvp,&m,r,in,out,0,xflags);
Blob b(hkvp.dump());
b.save(fntemp2);
//MM_MSG(" parsed outpint in "<<MMPR(fntemp2))
//typename HierKvp::name_list_type words;
//MM_MSG(" parsed outpint in "<<MMPR(fntemp2))
//typename HierKvp::name_list_type words;
//1 html 5 head 206 link 209 href 210 text = https://static-exp1.licdn.com/sc/h/7a4ujxvsdnd4s3x8tftdnnmt2
//1 html 5 head 211 script 212 type 213 text = application/ld+json
//1 html 5 head 211 script 214 (null) =       {"@context":"http://schema.org","@graph":[{"@type":"Person","address":{"@ty


StrTy pagetitle="";
HierWords words;
words.push_back("head");
words.push_back("title");
{
KvpLut lt=hkvp.lookup(words,7);
KvNodes no=hkvp.find(words,7);
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { pagetitle+=(*jj).second; } 
} /// ii 
if (pagetitle.length()==0)
{
MM_LOOP(ii,no) {
//me["location"]+=(*ii)->value(); 
StrTy v=(*ii)->value();
pagetitle+=v;
//IdxTy j=0;
//MM_LOOP(jj,cfv) { (*jj).m()["title"]+=v; ++j; } // jj 
}
}
}
//MM_MSG(MMPR(pagetitle))
words.clear();
#if 0 
/*
// they put it in the body, i
 0 html  1 body  2 div  3 div  4 main  5 div  6 div  7 div  8 div  9 article  10 script  11 type 12 text application/ld+json
 0 html  1 body  2 div  3 div  4 main  5 div  6 div  7 div  8 div  9 article  10 script  11 (null)     {"@context":"http:\/\/
 https://www.tctmd.com/news/anticoagulation-use-quickly-rose-meet-covid-19-need
*/
#endif

//words.push_back("head");
words.push_back("script");
//m.clear();
KvpLut lt=hkvp.lookup(words,1);
KvNodes no=hkvp.find(words,1);
//MM_MSG(MMPR(lt.size()))
KvNodes noa;
MM_LOOP(ii,no) { if ((*ii)->matching_value("application/ld+json",1)) noa.push_back((*ii)); }
no=noa;
//lt=lta;

lt.clear();
MM_LOOP(ii,no) { HierMap m; (*ii)->kvp_map_get(m); lt.push_back(m); }
//MM_MSG(MMPR(lt.size()))
IdxTy l=0; 
const bool dump_nodes=0;
MM_LOOP(ii,lt)
{
//if (dump_nodes) { MM_LOOP(jj,(*ii)) { MM_MSG(MMPR3(l,(*jj).first,(*jj).second)) } } // dump
if (dump_nodes) { MM_LOOP(jj,(*ii)) { MM_ERR(MMPR2(l,(*jj).first)) } } // dump
// these  keyts should be lower case earlier... 
//const typename HierKvp::kvp_map_type & kvp=(*ii);
StrTy nkey,nval;
HierMap & kvp=(*ii);
///IdxTy fkvrc=find_key_val(nkey,nval,kvp,l,0);
const StrTy json=kvp["(null)"];
//MM_MSG(" have a json "<<MMPR(kvp["(null)"] ))
IdxTy j=0; 
//MM_LOOP(jj,cfv) { (*jj).f()((*jj).m(),nkey,nval,kvp,0,0); ++j; } // jj 
const StrTy s4new=" sed -e 's/\\\\\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json - 1\"  | sed -e 's/mjm_bs_bug/\\\\\\\"/g'  "; 
IdxTy rcj=cmd_exec(fntempj,json,s4new,out,8);
Blob bio;
bio.clear();
bio.load(fntempj);
//MM_MSG(MMPR(StrTy(bio)))
///////////////////////////////
Ragged rj;
rj.kvp_splitter(0);
rj.load(fntempj); // this should now be parsed with value in the end of line field 
/// look for the "Person" graph on linked in and use that for bio biblio
if (combine) rags.push_back(rj);
HierKvp hkvpj;
Kvp mj,mej;
Kvvp kvvp;
IdxTy lhrc=load_hier_obj(hkvpj,&kvvp,&mj,rj,in,out,0,xflags);
Blob b(hkvpj.dump());
b.save(fntemp2);
bool li_person=false;
////////////////////////////////////////////////
// multiple graphs and ld+json entries need to eb accomodated...

words.clear();
words.push_back("@type");
KvNodes no=hkvpj.find(words,1);
//MM_MSG(MMPR(lt.size()))
KvNodes noa;
IdxTy minlvl=~0; 
IdxTy minplvl=~0; 

MM_ERR(MMPR(no.size()))
MM_LOOP(ii,no) 
{ 
const IdxTy lii=(*ii)->level();
if (lii<minlvl) minlvl=lii;
// this is not unique to LI wth wtf 
if ((*ii)->matching_value("Person",1))
{
if (lii<minplvl) minplvl=lii;
 noa.push_back((*ii)); 
}
MM_ERR(MMPR2((*ii)->info(),noa.size()))
}
if (minlvl!= ~0) if (minplvl==minlvl) li_person=true; 
MM_ERR(MMPR3(minlvl,minplvl, li_person));



//////////////////////////////////////////////
if (li_person) { 
words.clear();
words.push_back("addressLocality");
auto p=hkvpj.root()->matching_value("addressLocality",1); 
if (p) mj["location"]=p->value();
else
{ 
// this is probably wrong.. 
// li_person=false;
MM_ERR(" did not find addressLocality so reset ")
}  // p 
} // li_perseon
//////////////////////////////
//IdxTy rcs=sep_app_json_new(fnbibnew,nm,in,fntempj,out,bflags,xflags);

StrTy ename="";
const IdxTy bflags=0;
const StrTy fnbibnew=out.fn("bibtex",1); //
// this adds it... 
//IdxTy rcs=sep_app_json_new(fnbibnew,nm,in,fntempj,out,bflags,xflags);
//but this is the real logic,
// this is multiple words 
StrTy type=m["X_type"];
type="article";
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep=true||(type=="NewsArticle");
const bool include_all=true;
if (include_all) mej=mj;
// callinhg both does things like duplicate names etc...
StrTy url=uin;

// the precedence here is awful, it puts in blanks etc... 

IdxTy na=NewsArticleKvpMap(type,mej,mj,out,0);
if (m["X_type"]=="NewsArticle") { type="article"; me["subtype"]="news"; }
//IdxTy lip=LiPersonKvpMap(type,mej,mj,rj,out,0);
IdxTy lip=li_person?LiPersonKvpMapNewer(type,&kvvp,mej,mj,hkvpj,out,0):0;
if (mej.find("title")==mej.end()) mej["title"]=pagetitle;
mej["pagetitle"]=pagetitle;
MiscJson2(mej,out);
const IdxTy nstart=out.found();
const bool all=in.collect_all();

// X_mainEntity_author_name
if (combine) kvps.push_back(mej);
else
{
IdxTy rc=synthesize(in, out, url, nm, fnbibnew, ename, type, mej, ok_to_keep?0:1);
const IdxTy nfound=out.found();
MM_ERR("should have saved to "<<MMPR4(fnbibnew,type,nstart,nfound))

}
////////////////////////////////////////////////////////////////////////
++l; 
} // ii
if (combine)
{
Kvp mej;
StrTy url=uin;
StrTy  // type=m["X_type"];
type="article";
const bool ok_to_keep=true;
StrTy ename="";
IdxTy rcc=Combine(  mej,  kvps,  rags, 0);
//if (mej.find("real-type")!=mej.end()) type=mej["real-type"];
IdxTy rc=synthesize(in, out, url, nm, fnbib, ename, type, mej, ok_to_keep?0:1);
const IdxTy nfound=out.found();
MM_ERR("should have saved to "<<MMPR3(fnbib,type,nfound))
} // combine 
out.exit(nm); 
 return 0;

} // handleldjson2Old
#endif
////////////////////////////////////////////////////////////////////

 static  IdxTy handleldjson2(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const bool li=Bit(xflags,0);
const StrTy nm=StrTy("handleldjson2");
const bool all=in.collect_all();
const IdxTy nstart=out.found();
IdxTy rx=handleldjson2var(in,out,(li?0:1),xflags);
// only do this with "all" not just as a second attempt with good bib 
if (all ) handleldjson2var(in,out,0,xflags);
MM_ERR(MMPR4(nm,all,nstart,out.found()))
return 0;
}
 static  IdxTy handleldjson2var(const InpTy & in , OutTy & out ,const IdxTy opt,  const IdxTy xflags=0)  
{
const bool only_json=Bit(opt,0);
const bool only_good_json=Bit(opt,1);
const StrTy nm=StrTy("handleldjson2")+(only_json?"(only ldjson)":"(all)");
out.enter(nm); 
///////////////////////////////////////////
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
//StrTy fntemph=out.fn("temp",1); //
//StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const bool combine=true;
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
KvpVec kvps;
RagVec rags;
// should try to get with wget and if zero len comes back, wait try with chrome 
const IdxTy insz=in.getter().size(fn);
StrTy ispdf;
IdxTy idtc=identify(ispdf,fn,"",out,0);
//if (false) //  (headlessunsafe(in,out,xflags)) 
if (( ispdf=="pdf" ))
{
MM_ERR(" pointless to try pdf here "<<MMPR3(nm,fn,fn))
out.exit(nm); 
 return 0;
}
MM_ERR(" not getting as me seems ok   "<<MMPR(fn))
MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))
//////////////////////////////////////////////
// may want to use me...
/////////////////////////////////////////////
Ragged h,j;
// better but has change of "="
if (true) 
{
typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const IdxTy fnty=256+(only_good_json?8:0)+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rc=bt.parse_both(h,j,fn,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))
} 
//////////////////////////////////////////////

HierKvp hkvp;
Kvp m,me;
StrTy pagetitle="";
////////////////////////////////////////////////
IdxTy ptrc=GetPageTitle(pagetitle,in,out,h,m,hkvp,xflags);



#if 0 
/*
// they put it in the body, i
 0 html  1 body  2 div  3 div  4 main  5 div  6 div  7 div  8 div  9 article  10 script  11 type 12 text application/ld+json
 0 html  1 body  2 div  3 div  4 main  5 div  6 div  7 div  8 div  9 article  10 script  11 (null)     {"@context":"http:\/\/
 https://www.tctmd.com/news/anticoagulation-use-quickly-rose-meet-covid-19-need
*/
#endif

//words.push_back("head");
const bool dump_nodes=0;
/// look for the "Person" graph on linked in and use that for bio biblio
HierKvp hkvpj;
Kvp mj,mej;
Kvvp kvvp;
IdxTy lhrc=load_hier_obj(hkvpj,&kvvp,&mj,j,in,out,0,xflags);
{ Blob b(hkvpj.dump()); b.save(fntemp2); } 
////////////////////////////////////////////////
// multiple graphs and ld+json entries need to eb accomodated...
IdxTy tyty=0;
IdxTy rctype=FigureOutType( tyty, hkvpj,0);

bool li_person=(tyty&1);

//HierKvp hkvp;
//Kvp m,me;
IdxTy rcpe=PersonEtc(mj,hkvpj,tyty,0);

//StrTy ename="";
const IdxTy bflags=0;
const StrTy fnbibnew=out.fn("bibtex",1); //
//but this is the real logic,
// this is multiple words 
StrTy type=m["X_type"];
type="article";
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
//const bool ok_to_keep=true||(type=="NewsArticle");
const bool include_all=true;
if (include_all) mej=mj;
// callinhg both does things like duplicate names etc...
// this needs to be blank for local files with url info 
//StrTy url=uin;

// the precedence here is awful, it puts in blanks etc... 

IdxTy na=NewsArticleKvpMap(type,mej,mj,out,0);
if (m["X_type"]=="NewsArticle") { type="article"; me["subtype"]="news"; }
//IdxTy lip=LiPersonKvpMap(type,mej,mj,rj,out,0);
IdxTy lip=li_person?LiPersonKvpMapNewer(type,&kvvp,mej,mj,hkvpj,out,0):0;
if (mej.find("title")==mej.end()) mej["title"]=pagetitle;
if (mej.find("headline")!=mej.end()) mej["title"]=mej["headline"];
if (mej.find("X_headline")!=mej.end()) mej["title"]=mej["X_headline"];
mej["pagetitle"]=pagetitle;
MiscJson2(mej,out);
// X_mainEntity_author_name
 kvps.push_back(mej);
rags.push_back(j);
const IdxTy nfound=out.found();
MM_ERR("should have saved to "<<MMPR4(fnbibnew,type,nstart,nfound))

////////////////////////////////////////////////////////////////////////
//++l; 
// but have a crapped up one doh
//Kvp mej;
// not appropriate for local files
StrTy url=uin;
//StrTy  // type=m["X_type"];
type="article";
//const bool ok_to_keep=true;
StrTy ename="";
IdxTy rcc=Combine(  mej,  kvps,  rags, 0);
//if (mej.find("real-type")!=mej.end()) type=mej["real-type"];
IdxTy rc=synthesize(in, out, url, nm, fnbib, ename, type, mej, 0);
{const IdxTy nfound=out.found();
MM_ERR("should have saved to "<<MMPR3(fnbib,type,nfound))
}
out.exit(nm); 
 return 0;

} // handleldjson2

static IdxTy PersonEtc(Kvp & mj, HierKvp & hkvpj, const IdxTy tyty, const IdxTy flags)
{
const bool li_person=(tyty&1);
if (li_person) { 
HierWords words;
//words.clear();
words.push_back("addressLocality");
auto p=hkvpj.root()->matching_value("addressLocality",1); 
if (p) mj["location"]=p->value();
else
{ 
// this is probably wrong.. 
// li_person=false;
MM_ERR(" did not find addressLocality so reset ")
}  // p 
} // li_perseon
return 0; 
} // PersonEtc

static IdxTy FigureOutType(IdxTy&  tyty, HierKvp & hkvpj,const IdxTy flags)
{
HierWords words;
//words.clear();
words.push_back("@type");
KvNodes no=hkvpj.find(words,1);
//MM_MSG(MMPR(lt.size()))
KvNodes noa;
IdxTy minlvl=~0; 
IdxTy minplvl=~0; 

MM_ERR(MMPR(no.size()))
MM_LOOP(ii,no) 
{ 
const IdxTy lii=(*ii)->level();
if (lii<minlvl) minlvl=lii;
// this is not unique to LI wth wtf 
if ((*ii)->matching_value("Person",1))
{
if (lii<minplvl) minplvl=lii;
 noa.push_back((*ii)); 
}
MM_ERR(MMPR2((*ii)->info(),noa.size()))
} // ii 
if (minlvl!= ~0) if (minplvl==minlvl) tyty|=1; 
MM_ERR(MMPR3(minlvl,minplvl, tyty));
return 0; 
} // FigureOutType 

static IdxTy GetPageTitle(StrTy&pagetitle,const InpTy & in,OutTy & out,Ragged & h,Kvp & m,HierKvp & hkvp,const IdxTy xflags)
{
const bool remove_crap=!Bit(xflags,0);
StrTy fntemp=out.fn("temp",1); //
IdxTy lhrc=load_hier_obj(hkvp,&m,h,in,out,0,xflags);
Blob b(hkvp.dump()); b.save(fntemp);
MM_ERR(" saved hkvp blod "<<MMPR2(fntemp,__FUNCTION__))
//1 html 5 head 206 link 209 href 210 text = https://static-exp1.licdn.com/sc/h/7a4ujxvsdnd4s3x8tftdnnmt2
//1 html 5 head 211 script 212 type 213 text = application/ld+json
//1 html 5 head 211 script 214 (null) =       {"@context":"http://schema.org","@graph":[{"@type":"Person","address":{"@ty
HierWords words;
words.push_back("head");
words.push_back("title");
KvpLut lt=hkvp.lookup(words,7);
KvNodes no=hkvp.find(words,7);
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { pagetitle+=(*jj).second; } 
} /// ii 
if (pagetitle.length()==0) 
{ MM_LOOP(ii,no) { StrTy v=(*ii)->value(); pagetitle+=v;  } }
if (remove_crap)
{
const char *p=pagetitle.c_str();
const IdxTy len=pagetitle.length();
char c[len*2+2];
IdxTy pc=0;
for(IdxTy i=0; i<len; ++i)
{
 char cp=p[i];
if (cp=='|') cp=' ';
else if (cp=='#') cp=' ';
c[pc]=cp;
++pc;
} // i 
c[pc]=0;
pagetitle=StrTy(c);
} // remove_crap



//MM_MSG(MMPR(pagetitle))
///////////////////////////////////////////////////////
return 0;
} 

static IdxTy GetPublishedDate(StrTy&pagetitle,const StrTy & xxx,const InpTy & in,OutTy & out,Ragged & h,Kvp & m,HierKvp & hkvp,const IdxTy xflags)
{
const bool remove_crap=!Bit(xflags,0);
StrTy fntemp=out.fn("temp",1); //
IdxTy lhrc=load_hier_obj(hkvp,&m,h,in,out,0,xflags);
Blob b(hkvp.dump()); b.save(fntemp);
MM_ERR(" saved hkvp blod "<<MMPR2(fntemp,__FUNCTION__))
MM_ERR(" saved hkvp blod "<<MMPR2(hkvp.size(),__FUNCTION__))
//1 html 5 head 206 link 209 href 210 text = https://static-exp1.licdn.com/sc/h/7a4ujxvsdnd4s3x8tftdnnmt2
//1 html 5 head 211 script 212 type 213 text = application/ld+json
//1 html 5 head 211 script 214 (null) =       {"@context":"http://schema.org","@graph":[{"@type":"Person","address":{"@ty
HierWords words;
words.push_back("meta");
words.push_back(xxx);
KvpLut lt=hkvp.lookup(words,7);
KvNodes no=hkvp.find(words,7);
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { pagetitle+=(*jj).second; } 
} /// ii 
if (pagetitle.length()==0) 
{ MM_LOOP(ii,no) { StrTy v=(*ii)->value(); pagetitle+=v;  } }
if (remove_crap)
{
const char *p=pagetitle.c_str();
const IdxTy len=pagetitle.length();
char c[len*2+2];
IdxTy pc=0;
for(IdxTy i=0; i<len; ++i)
{
 char cp=p[i];
if (cp=='|') cp=' ';
else if (cp=='#') cp=' ';
c[pc]=cp;
++pc;
} // i 
c[pc]=0;
pagetitle=StrTy(c);
} // remove_crap



//MM_MSG(MMPR(pagetitle))
///////////////////////////////////////////////////////
return 0;
} 


////////////////////////////////////////////////////////////////////
static void AddAnd(Kvp & m, const StrTy & k, const StrTy & v, const IdxTy flags=0)
{
const bool space=Bit(flags,0);
if (v.length()==0) return;
const auto& ii=m.find(k);
if (ii==m.end()) {m[k]=v; return ; } 
const StrTy & vo=(*ii).second;
if (v==vo) return; 
if (vo.length()==0) { m[k]=v; return; } 
if (space) (*ii).second+=" "+v;
else (*ii).second+=" and "+v;
}
// https://www.sra-solder.com/bar-solder-96-4-lead-free
typedef std::map<StrTy, StrTy> LitMap;
static LitMap& Lits() 
{
static LitMap m;
if (m.size()) return m;
m["sku"]="";
m["email"]="";
m["telephone"]="";

return m;
} // Lits
static IdxTy Combine( Kvp & m, KvpVec & kvps, RagVec & rags, const IdxTy flags)
{
IdxTy rc=0;
MM_LOOP(ii,kvps)
{
MM_LOOP(jj,(*ii))
{
AddAnd(m,(*jj).first,(*jj).second,0); 
} // jj 

} // ii
LitMap lits=Lits();

MM_LOOP(ii,rags)
{
const Ragged & r=(*ii);
const IdxTy sz=r.size();
StrTy rtype="";
MM_LOOP(jj,r)
{
const Ragged::Line & l=(*jj);
const IdxTy len=l.size();
if (len<3) continue;
const IdxTy off=(l[len-2]=="=")?1:0; // equals should not really be there doh
const StrTy & k1=l[len-2-off];
 StrTy  v=l[len-1];
char fick[v.length()+10];
IdxTy pc=0;
for(IdxTy i=0; i<v.length(); ++i)
{
const char c=v.c_str()[i];
if (c>31) if (c<128) {fick[pc]=c; ++pc;} 
}
fick[pc]=0;
v=StrTy(fick);
// if here is at least one Product entry that is the real type
if (k1=="@type")
	{ if (v=="Product"){ m["type"]="product";  }rtype=v; continue;  }
if (lits.find(k1)!=lits.end()){  m[k1]=v; continue; } 

if (k1=="name")
{
if (rtype=="Brand") { m["publisher"]=v; continue; }
if (rtype=="Product") { m["title"]=v; continue; }
// https://cen.acs.org/physical-chemistry/Magnet-doubles-hydrogen-yield-water/97/web/2019/06
if (rtype=="Organization") {if (m.find("author")==m.end())  m["author"]=v; continue; }
} // name

if (k1=="streetAddress") { m["location"]=v; continue; }
if (k1=="addressLocality") { AddAnd(m, "location",v,1) ; continue; }
if (k1=="addressRegion") { AddAnd(m, "location",v,1) ; continue; }
if (k1=="postalCode") { AddAnd(m, "location",v,1) ; continue; }
if (k1=="addressCountry") { AddAnd(m, "location",v,1) ; m["country"]=v;  continue; }
if (k1=="sameAs") { AddAnd(m, "aliases",v) ; continue; }

if (k1=="price") { AddAnd(m, "cost",v,1) ; continue; }
if (k1=="priceCurrency") { AddAnd(m, "cost",v,1) ; continue; }
if (k1=="priceValidUntil") { AddAnd(m, "cost",StrTy("Valid :")+v,1) ; continue; }
if ((k1=="doi") ||(k1=="DOI"))
{
if (m.find("doi")==m.end()) m["doi"]=v; 
continue;
} // doi

} // jj

} // ii 
MM_LOOP(ii,m)
{
MM_ERR( "SAF "<<MMPR2((*ii).first,(*ii).second))

} // ii 
return rc;
} // combine
static IdxTy AddOneX(Kvp & _mej, OutTy & out )
{
Kvp m;
MM_LOOP(ii,_mej)
{
const StrTy k=(*ii).first;
const StrTy & v=(*ii).second;
const IdxTy len=k.length();
const char * p=k.c_str();
if (*p!='X') continue;
IdxTy n=0;
IdxTy i=0;
IdxTy off=0;
while (p[i])
{
if (p[i]=='_'){off=i+1;  ++n; if (n>1) break ; } 
++i;
} // p[i]
if (n!=1) continue; 
if (off==0) continue;
StrTy k2=StrTy(p+off);
if (_mej.find(k2)!=_mej.end()) if (_mej[k2].length()) continue ;  
m[k2]=v;

} // ii 
MM_LOOP(jj,m) { _mej[(*jj).first]=(*jj).second; }
return 0;
} // AddOneX
//typedef typename BibEntry::kvp_type Kvp;
static IdxTy MiscJson2(Kvp & mej, OutTy & out )
{

// this will generate spurious years for a year-oly date > 1999
MutateMissingKvp(mej,"date",mej,"X_mainEntity_datePublished","cat -",out);
MutateMissingKvp(mej,"day",mej,"date","date +%d -f-",out);
MutateMissingKvp(mej,"month",mej,"date","date +%m -f-",out);
MutateMissingKvp(mej,"year",mej,"date","date +%Y -f-",out);
MutateMissingKvp(mej,"abstract",mej,"X_mainEntity_description","cat -",out);
MutateMissingKvp(mej,"abstract",mej,"description","cat -",out);
MutateMissingKvp(mej,"author",mej,"X_mainEntity_author_name","cat -",out);
MutateMissingKvp(mej,"author",mej,"X_author_name","cat -",out);
MutateMissingKvp(mej,"author",mej,"creator","cat -",out);
MutateMissingKvp(mej,"author",mej,"X_creator_name","cat -",out);
MutateMissingKvp(mej,"issn",mej,"X_issn","cat -",out);
MutateMissingKvp(mej,"url",mej,"X_mainEntityOfPage_1_id","cat -",out);
MutateMissingKvp(mej,"affiliation",mej,"X_creator_affiliation","cat -",out);
MutateMissingKvp(mej,"journal",mej,"journal_title","cat -",out);
// x_inlanguage_name = {English}, // lc fick 
MutateMissingKvp(mej,"language",mej,"X_inLanguage_name","cat -",out);
MutateMissingKvp(mej,"keywords",mej,"X_keywords_keywords","cat -",out);
AddOneX(mej,out);
if ((mej.find("identifier")!=mej.end())
&&(mej.find("doi")==mej.end()))
{
const StrTy & x=mej["identifier"];
const IdxTy pos=mjm_strings::indexOf(x.c_str(),"doi");
if (pos!=~0) mej["doi"]=x;

} // identifier
/*
zenodo,
x_context = {https://schema.org/},
x_creator_1_affiliation = {self},
x_creator_1_name = {MIKE MARCHYWKA},
x_creator_1_type = {Person},
x_creator_affiliation = {self},
x_creator_name = {MIKE MARCHYWKA},
x_creator_type = {Person},
x_datepublished = {2022-03-27},
x_description = {<p>Informal description of dog with prior heartworm recovering strength from riboflavin.\&nbsp;</p>\x0a\x0a<p>Possible relationship to ongoing FDA investigation into mysterious diet related DCM.</p>},
x_headline = {Happy Again : Possible Canine Riboflavin Deficiency},
x_id = {https://doi.org/10.5281/zenodo.6388035},
x_identifier = {https://doi.org/10.5281/zenodo.6388035},
x_inlanguage_name = {English},
x_keywords_keywords = {riboflavin fda canine nutirition cardiomyopathy},
x_license = {https://creativecommons.org/licenses/by/2.0/legalcode},
x_type = {ScholarlyArticle},
x_url = {https://zenodo.org/record/6388035},
x_version = {0.50},





x_context = {https://schema.org/},
x_creator_1_affiliation = {self},
x_creator_1_name = {MIKE MARCHYWKA},
x_creator_1_type = {Person},
x_creator_affiliation = {self},
x_creator_name = {MIKE MARCHYWKA},
x_creator_type = {Person},

*/
return 0; 

} // MiscJson2



#if 0
//////////////////////////////////////////////
MM_ERR(" start ")
// the problem was not checking return code now it plows through
// but may still treat quote as invalid. 
// use the easier to partse outut 
 const StrTy s4new=" sed -e 's/\\\\\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json - 1\"  | sed -e 's/mjm_bs_bug/\\\\\\\"/g'  "; 

//MM_MSG(" parsing "<<  MMPR(ll[0]))
IdxTy rcj=cmd_exec(fntempj,ll[0],s4new,out,8);
bio.clear();
bio.load(fntempj);
//MM_MSG(MMPR(StrTy(bio)))
const IdxTy bflags=0;
const StrTy fnbibnew=out.fn("bibtex",1); //
// this adds it... 
IdxTy rcs=sep_app_json_new(fnbibnew,nm,in,fntempj,out,bflags,xflags);
//fntempj should now have the parsed json... 
// this parses the FIL not the blob 
//if (clean) in.finisher().light_clean(fnbib,fn("temp"),0);
Blob b;
// These will not parse because of special cars and maybe blank
// name. But files are ok.  
Frc frc=in.finisher().parse(b,fnbibnew,"",0);
//Frc frcb=0; // in.finisher().validate(b,fnbib,"",0);
BibEntry & ber=frc.be();
be[i]=ber;
//MM_MSG(MMPR3(i,fnbibnew,be[i].format()))

//} // i 
#endif

///////////////////////////////////////////












////////////////////////////////////////////////////////////////

/*
toobib -json1 junk/marco.phtml  | grep itle
./mjsonu.h581  cip.cmd()=sax-json good=1 eof=0 use_stdin=0 select= flags=1 cmd1= cmd2=
./mjsonu.h280 GLOBAL 368  1  573  2 @graph 1  6 jobTitle 1 jobTitle = Freelance Technical Writer
./mjsonu.h622   json dies before parse fick  ok=1 eof=1 good=0 n=2459705 cc=r



 1 html  5 head  17 script  18 id  19 text = litespeed-webfont-lib
 1 html  5 head  17 script  20 src  21 text = https://www.newfoodmagazine.com/wp-content/plugins/litespeed-cache/assets/js/webfontloader.min.js
 1 html  5 head  23 meta  24 name  25 text = viewport
 1 html  5 head  23 meta  26 content  27 text = width=device-width, initial-scale=1.0
 1 html  5 head  28 link  29 href  30 text = /favicon.ico?v=2
 1 html  5 head  28 link  31 type  32 text = image/x-icon
 1 html  5 head  28 link  33 rel  34 text = shortcut icon
 1 html  5 head  35 title  36 text = Lack of fruit and veg responsible for 18,000 early deaths says report


obsolete format, 
 html 1  head 5  link 9  href 12  text 13 = https://fonts.gstatic.com/




*/

static  IdxTy 
    load_hier_obj(HierKvp & hkvp ,  Kvp* mp,const Ragged & r,const InpTy & in ,  OutTy & out ,const IdxTy bflags,  const IdxTy xflags=0)  
{
Kvvp*  kvvp=0;
return load_hier_obj(hkvp,kvvp,mp,r,in,out,bflags,xflags); 
}

template <class Tkv > static  IdxTy 
    load_hier_obj(HierKvp & hkvp , Tkv * kvvp,  Kvp* mp,const Ragged & r,const InpTy & in ,  OutTy & out ,const IdxTy bflags,  const IdxTy xflags=0)  

{
MM_ERR(MMPR2(__FUNCTION__,r.size()))
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l=r[i];
//MM_MSG(r.debug_line(i))
const IdxTy len=l.size();
if (len<4)
{
MM_ERR(MMPR3(__FUNCTION__,i,r.debug_line(i)))
 continue;

}
//MM_ERR(" ok " << MMPR3(__FUNCTION__,i,r.debug_line(i)))
typename HierKvp::V vv;
const StrTy & val=l[len-1];
StrTy key="";
StrTy key_lvl="";
IdxTy j=len-2;
while (true)
{
const StrTy p=l[j];
if (p=="GLOBAL")
{
MM_ERR(MMPR4(__FUNCTION__,i,len,hkvp.size()))
 break;
}
//MM_MSG(MMPR4(j,l[j],key,val))
//if (p.length()==0) break;
if (p.length()!=0) {  
const char * c=p.c_str();
IdxTy ini=0;
if (c[0]=='@') ini=1;
if(0!=key.length())
{
key=StrTy(c+ini)+"_"+key;
Ss ss; ss<<l[j+1];
key_lvl=StrTy(c+ini)+"_"+ss.str()+"_"+key_lvl;
}
else{  key=StrTy(c+ini); key_lvl=key; } 
} // length
if (j<1)
{
MM_ERR(" non json input limit "<<MMPR4(__FUNCTION__,i,j,len))
 break;
}
vv.push_back(l[j-1]);
vv.push_back(p);
if (j<2) break;
j-=2;
//j-=1;
} // j
hkvp.add(vv,val);
if (mp)
{
Kvp & m=*mp; 
key=StrTy("X_")+key;
key_lvl=StrTy("X_")+key_lvl;
if (kvvp) (*kvvp)[key].push_back(val);
// use db or something here. 
auto ii=m.find(key);
if (ii==m.end()) m[key]=val;
else  (*ii).second+=StrTy(" and ")+val;
if (key!=key_lvl) { 
ii=m.find(key_lvl);
if (ii==m.end()) m[key_lvl]=val;
else  (*ii).second+=StrTy(" and ")+val;
} // X_headline on NPR
//MM_MSG(MMPR4(key,val,m[key],i))
} // mp
} // i
//MM_MSG(MMPR(hkvp.dump()))

return 0; 
}  // load_hier_obj



static  IdxTy 
    sep_app_json_new(const StrTy & fnbibnew , const StrTy & nm,const InpTy & in , const StrTy & jfn, OutTy & out ,const IdxTy bflags,  const IdxTy xflags=0)  
{
Ragged r;
r.kvp_splitter(0);
r.load(jfn); // this should now be parsed with value in the end of line field 
StrTy ename="";
StrTy type="";
#if 0
./mjsonu.h280 GLOBAL 10  8 datePublished = 2021-06-06T22:41:23+00:00
./mjsonu.h280 GLOBAL 10  9 publisher 1 @type = Organization
./mjsonu.h280 GLOBAL 10  9 publisher 2 name = Slate
./mjsonu.h280 GLOBAL 10  9 publisher 3 logo 1 @type = ImageObject
#endif
HierKvp hkvp;
Kvp m,me;
IdxTy lhrc=load_hier_obj(hkvp,&m,r,in,out,bflags,xflags);
const StrTy url=in.uin();
#if 0 
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l=r[i];
//MM_MSG(r.debug_line(i))
const IdxTy len=l.size();
if (len<2) continue;
typename HierKvp::V vv;
const StrTy & val=l[len-1];
StrTy key="";
IdxTy j=len-2;
while (true)
{
const StrTy p=l[j];
if (p=="GLOBAL") break;
//MM_MSG(MMPR4(j,l[j],key,val))
//if (p.length()==0) break;
if (p.length()!=0) {  
const char * c=p.c_str();
IdxTy ini=0;
if (c[0]=='@') ini=1;
if(0!=key.length())
key=StrTy(c+ini)+"_"+key;
else key=StrTy(c+ini);
} // length
vv.push_back(l[j-1]);
vv.push_back(p);
if (j<2) break;
j-=2;
//j-=1;
} // j
hkvp.add(vv,val);
key=StrTy("X_")+key;
// use db or something here. 
auto ii=m.find(key);
if (ii==m.end()) m[key]=val;
else  (*ii).second+=StrTy(" , ")+val;
//MM_MSG(MMPR4(key,val,m[key],i))
} // i
MM_MSG(MMPR(hkvp.dump()))
#endif

// this is multiple words 
type=m["X_type"];
type="article"; 
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep=(type=="NewsArticle");
const bool include_all=true;
if (include_all) me=m;
// callinhg both does things like duplicate names etc...

IdxTy na=NewsArticleKvpMap(type,me,m,out,0);
if (m["X_type"]=="NewsArticle") { type="article"; me["subtype"]="news"; }
IdxTy lip=LiPersonKvpMap(type,me,m,r,out,0);
IdxTy rc=synthesize(in, out, url, nm, fnbibnew, ename, type, me, ok_to_keep?0:1);
MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))
return rc; 
} //sep_app_json_new


static IdxTy
    sep_app_json(const StrTy & fnbibnew , const StrTy & nm,const InpTy & in , const StrTy & jfn, OutTy & out ,const IdxTy bflags,  const IdxTy xflags=0)  
{
//const StrTy nm="guesswikipedia";
//out.enter(nm);
//const StrTy & uin=in.uin();

// IdxTy rcs=sep_app_json(be[i],fntempj,in,out,bflags,xflags);
//IdxTy rc=0;
Blob b;
b.load(jfn);
const StrTy bs=StrTy(b);
//const StrTy datec=" date";
//const StrTy date=MutateOnly(bs,datec,out);
//const StrTy yearc=" date +%Y";
//const StrTy year=MutateOnly(bs,yearc,out);
// date +%Y-%m+%d  -d"2021-06-03T16:42:44.751Z"

Kvp m,me;
const StrTy url=in.uin();
//const StrTy name="";
//const StrTy type="person";
const StrTy typec="grep \"GLOBAL[^a-zA-Z]*@type=\" | sed -e 's/.*GLOBAL[^a-zA-Z]*@type=//' "; 
StrTy type=MutateOnly(bs,typec,out);
//const StrTy nm="";
StrTy ename="";

//12 ./mjsonu.h205 ;GLOBAL[1];[15];author[1];[1]:@context=http://schema.org/
const StrTy sragold="sed -e 's/[^;]*;//' | sed -e ':a;s/\\([a-zA-Z]\\)[^_a-zA-Z:][^_a-zA-Z]*\\([_a-zA-Z].*:\\)/\\1_\\2/g;ta' | sed -e ':a;s/\\(^[^=]*\\)\\[[0-9]*\\][:;]/\\1_/;ta' | sed -e 's/=/ /' | awk '{ $1=gensub(/[@_][@_]*/,\"_\",\"g\",$1); $1=gensub(/GLOBAL/,\"\",1,$1);  print $0; } '| sed -e 's/^__*//'   ";

// change the equals into a space and then get first for awk, 
//const StrTy srag="sed -e 's/^[^=]* ;GLOBAL/GLOBAL/' | sed -e 's/=/ /' | awk '{ $1=gensub(/\\[[0-9]*\\]/,\"_\",\"g\",$1); $1=gensub(/[@_][@_]*/,\"_\",\"g\",$1); $1=gensub(/[;:@]/,\"_\",\"g\",$1);$1=gensub(/__*/,\"_\",\"g\",$1);   print $0; }' | sed -e 's/^__*//' ";  
const StrTy srag="sed -e 's/^[^=]* ;GLOBAL/X_/' | sed -e 's/=/ /' | awk '{ $1=gensub(/\\[[0-9]*\\]/,\"_\",\"g\",$1); $1=gensub(/[@_][@_]*/,\"_\",\"g\",$1); $1=gensub(/[;:@]/,\"_\",\"g\",$1);$1=gensub(/__*/,\"_\",\"g\",$1);   print $0; }' | sed -e 's/^__*//' ";  

StrTy fnrag=out.fn("temp",1); //
IdxTy rcj=cmd_exec(fnrag,bs,srag,out,8);
MM_ERR(" ragged "<<MMPR(fnrag))
Ragged r;
r.code_splitter(1);
r.load(fnrag);
MM_SZ_LOOP(i,r,rsz)
{
const Ragged::Line & l=r[i];
const IdxTy sz=l.size();
if (sz!=2) { MM_ERR(" bad size "<<MMPR2(sz,i)) ;  continue; } 
auto ii=m.find(l[0]);
if (ii==m.end()) m[l[0]]=l[1];
else  (*ii).second+=StrTy(" , ")+l[1];
//m[l[0]]+=l[1];
} // i 

// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep=(type=="NewsArticle");
const bool include_all=true;
if (include_all) me=m;
// callinhg both does things like duplicate names etc...

IdxTy na=NewsArticleKvpMap(type,me,m,out,0);
if (m["X_type"]=="NewsArticle") { type="article"; me["subtype"]="news"; }



IdxTy lip=LiPersonKvpMap(type,me,m,bs,out,0);



IdxTy rc=synthesize(in, out, url, nm, fnbibnew, ename, type, me, ok_to_keep?0:1);
MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))

return rc; 

} // sep_add_json
#if 0
@missing{default,
X_context = {http://schema.org},
X_graph_address_addressLocality = {Greater Melbourne Area},
X_graph_address_type = {PostalAddress},
X_graph_alumniOf_name = {Box Hill Institute , Fortress Learning , University of Southern Mississippi , Universidade do Porto , Instituto de CiÃªncias BiomÃ©dicas Abel Salazar},
X_graph_alumniOf_type = {EducationalOrganization , EducationalOrganization , EducationalOrganization , EducationalOrganization , EducationalOrganization},
X_graph_alumniOf_url = {https://www.linkedin.com/school/box-hill-institute/ , https://www.linkedin.com/school/university-of-southern-mississippi/ , https://www.linkedin.com/school/universidadedoporto/},
X_graph_image_contentUrl = {https://media-exp3.licdn.com/dms/image/C5103AQGKOL38Bki0rA/profile-displayphoto-shrink_200_200/0/1518952072007?e=1628121600&v=beta&t=Z6MfIfR0Zn_H9ZX1SCCQbCCwbftan_MD2Mt2bAwOw2k},
X_graph_image_type = {ImageObject},
X_graph_name = {JosÃ© Gomes},
X_graph_reviewedBy_name = {JosÃ© Gomes},
X_graph_reviewedBy_type = {Person},
X_graph_sameAs = {https://au.linkedin.com/in/josegomestmp},
X_graph_type = {Person , WebPage},
X_graph_url = {https://au.linkedin.com/in/josegomestmp , https://au.linkedin.com/in/josegomestmp}
}

#endif

static IdxTy LiPersonKvpMapNewer(StrTy & type, Kvvp* kvvp,Kvp &me, const Kvp & m, HierKvp  & hkvpj,OutTy & out, const IdxTy flags)
{
me["real_type"]="person";
me["year"]=MutateOnly("","date +%Y",out);
me["month"]=MutateOnly("","date +%m",out);
me["day"]=MutateOnly("","date +%d",out);
// may have been picked up by other thing 
CopyKvp(me,"author",m,"X_graph_name",1);
//CopyKvp(me,"location",m,"X_graph_address_addressLocality");
me["name"]=me["author"];

//./mjsonu.h280 GLOBAL 368  1  573  2 @graph 1  6 jobTitle 1 jobTitle = Freelance Technical Writer

HierWords words;
words.push_back("addressLocality");
KvpLut lt=hkvpj.lookup(words,1);
KvNodes no=hkvpj.find(words,1);
//MM_MSG(MMPR2(lt.size(),no.size()))
MM_LOOP(ii,no) { me["location"]+=(*ii)->value(); } 
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { me[(*jj).first]+=(*jj).second; }
}
lt.clear();
MM_LOOP(ii,no) { HierMap mh; (*ii)->kvp_map_get(mh); lt.push_back(mh); //  }
me["location2"]+=mh["addressLocality"]; } 
words.clear();
words.push_back("alumniOf");
words.push_back("name");
no=hkvpj.find(words,1);
// there is some method for this 
MM_LOOP(ii,no) { me["alumni"]+=(*ii)->value(); } 
if (kvvp)
{
Kvvp& kvv=*kvvp;
GetIf(me,kvv,"title", "X_graph_jobTitle_jobTitle",1);
GetIf(me,kvv,"title", "X_graph_worksFor_name",1);
GetIf(me,kvv,"url", "X_graph_url",0);

} // kvvp 
CopyKvp(me,"title",m,"X_graph_alumniOf_name");
me["title"]=me["location"]+StrTy(" ")+me["title"];


return 0; 
} // LiPersonKvpMap

static IdxTy GetIf(Kvp & me, Kvvp& kvv , const StrTy&d , const StrTy & s,const IdxTy flags ) 
{
const bool append_ok=Bit(flags,0);
auto ii=me.find(d);
auto jj=kvv.find(s);
const bool have_alt=(jj!=kvv.end())?((*jj).second.size()!=0):false;
const bool have_already=(ii!=me.end())?((*ii).second!=""):false;
if (have_already&&!append_ok) return 0;
if (have_alt){ 
if (have_already) { me[d]=me[d]+" "+((*jj).second)[0]; }  
else me[d]=((*jj).second)[0]; return 1; } 
return ~0; 
} // GetOI


/*

HierWords words;
words.push_back("head");
words.push_back("script");
//m.clear();
KvpLut lt=hkvp.lookup(words,1);
typedef typename HierKvp::candidate_nodes_type KvNodes;
KvNodes no=hkvp.find(words,1);
MM_MSG(MMPR(lt.size()))
KvNodes noa;
MM_LOOP(ii,no) { if ((*ii)->matching_value("application/ld+json",1)) noa.push_back((*ii)); }
no=noa;
*/
//me["jt"]=jt;
//me["title"]+=jt;


static IdxTy LiPersonKvpMap(StrTy & type, Kvp &me, const Kvp & m, const Ragged  & r,OutTy & out, const IdxTy flags)
{
me["year"]=MutateOnly("","date +%Y",out);
me["month"]=MutateOnly("","date +%m",out);
me["day"]=MutateOnly("","date +%d",out);
// may have been picked up by other thing 
CopyKvp(me,"author",m,"X_graph_name",1);
CopyKvp(me,"location",m,"X_graph_address_addressLocality");
CopyKvp(me,"title",m,"X_graph_alumniOf_name");
me["name"]=me["author"];

//const StrTy title=MutateOnly(bs,titlec,out);
StrTy sep="";
IdxTy job=1;
StrTy jt;

MM_SZ_LOOP(i,r,rsz)
//while (true)
{
const Ragged::Line & l = r[i];
const IdxTy len=l.size();
// need the value plus name plus another field 
if (len<8) continue;
const StrTy & v=l[len-1];

//Ss ss; ss<<job; const StrTy si=ss.str();
// "works for comes first... "
//const StrTy tijc=
//" grep  -i  \"jobtitle\\["+si+"\\].*jobtitle=\\|worksfor\\["+si+"\\].*name=\" | sed -e 's/.*name=//' | sed -e 's/.*Title=//g'  | awk '{if (x==\"\") x=$0; else { x=$0\"(\"x\")\";  print x ;  x=\"\";  }     }' ";
//const StrTy tc=""; // MutateOnly(bs,tijc,out);
//if (tc=="") break;
//jt=jt+sep+tc;
//sep=",";
//++job;
//if (job>5) break;

} // i  while

me["jt"]=jt;
me["title"]+=jt;

return 0;
} // LiPersonKvpMap





static IdxTy LiPersonKvpMap(StrTy & type, Kvp &me, const Kvp & m, const StrTy & bs,OutTy & out, const IdxTy flags)
{
me["year"]=MutateOnly("","date +%Y",out);
me["month"]=MutateOnly("","date +%m",out);
me["day"]=MutateOnly("","date +%d",out);
// may have been picked up by other thing 
CopyKvp(me,"author",m,"X_graph_name",1);
CopyKvp(me,"location",m,"X_graph_address_addressLocality");
CopyKvp(me,"title",m,"X_graph_alumniOf_name");
me["name"]=me["author"];

//const StrTy title=MutateOnly(bs,titlec,out);
StrTy sep="";
IdxTy job=1;
StrTy jt;

while (true)
{
Ss ss;
ss<<job;
const StrTy si=ss.str();
// "works for comes first... "
const StrTy tijc=
" grep  -i  \"jobtitle\\["+si+"\\].*jobtitle=\\|worksfor\\["+si+"\\].*name=\" | sed -e 's/.*name=//' | sed -e 's/.*Title=//g'  | awk '{if (x==\"\") x=$0; else { x=$0\"(\"x\")\";  print x ;  x=\"\";  }     }' ";
const StrTy tc=MutateOnly(bs,tijc,out);
if (tc=="") break;
jt=jt+sep+tc;
sep=",";
++job;
if (job>5) break;

} // while

me["jt"]=jt;
me["title"]+=jt;

return 0;
} // LiPersonKvpMap

static IdxTy NewsArticleKvpMap(StrTy & type, Kvp &me, const Kvp & m, OutTy & out, const IdxTy flags)
{

// this flat approach needs to be fixed to hierarchy 
// some sites have the wrong escape level too 
CopyKvp(me,"author",m,"X_author_author",1);
CopyKvp(me,"author",m,"X_author_name",1);
CopyKvp(me,"author",m,"X_author_name_name",1);
CopyKvp(me,"author",m,"X_graph_author_author",1);
CopyKvp(me,"author",m,"X_graph_author_name",1);
CopyKvp(me,"author",m,"X_graph_name",1);

CopyKvp(me,"publisher",m,"X_publisher_name");
CopyKvp(me,"publisher",m,"X_graph_publisher_name");

CopyKvp(me,"journal",m,"X_publisher_name",1);
CopyKvp(me,"journal",m,"X_graph_publisher_name",1);
CopyKvp(me,"journal",m,"X_articleSection",1);
CopyKvp(me,"journal",m,"X_graph_articleSection",1);
CopyKvp(me,"title",m,"X_headline",1);
CopyKvp(me,"title",m,"X_graph_headline",1);
CopyKvp(me,"abstract",m,"X_description");
CopyKvp(me,"abstract",m,"X_graph_description");
const StrTy isolate_date="sed -e 's/.*=.//' | sed -e 's/[\">].*//'  | sed -e 's/ and .*//'  |date +%Y-%m-%d -f - | awk '{if ( $2 == \",\") { print $1} else { print $0} }' ";
//MutateKvp(me,"date",m,"X_datePublished","sed -e 's/.*=.//' | sed -e 's/[\">].*//'  |date +%Y-%m-%d -f -",out);
MutateKvp(me,"date",m,"X_datePublished",isolate_date+"  |date +%Y-%m-%d -f -",out);
//MutateKvp(me,"date",m,"X_graph_datePublished","sed -e 's/.*=.//' | sed -e 's/[\">].*//'  |date +%Y-%m-%d -f -",out);
MutateKvp(me,"date",m,"X_graph_datePublished",isolate_date+" |date +%Y-%m-%d -f -",out);
MutateKvp(me,"date",m,"X_uploadDate",isolate_date+" |date +%Y-%m-%d -f -",out);
MutateKvp(me,"date",m,"uploadDate",isolate_date+" |date +%Y-%m-%d -f -",out);
// final assembly should pick this up 
MutateKvp(me,"year",m,"X_datePublished","date +%Y -f -",out);
MutateKvp(me,"year",m,"X_graph_datePublished","date +%Y -f -",out);
MutateKvp(me,"month",m,"X_datePublished","date +%m -f -",out);
MutateKvp(me,"month",m,"X_graph_datePublished","date +%m -f -",out);
MutateKvp(me,"day",m,"X_datePublished","date +%d -f-",out);
MutateKvp(me,"day",m,"X_graph_datePublished","date +%d -f -",out);
return 0; 

}  // NewsArticleKvpMap

template <class Tm > 
static bool CopyKvp(Tm & d, const StrTy & ds, const Tm & s, const StrTy & ss, const IdxTy flags=0 ) 
{
const bool if_missing=Bit(flags,0);
const auto ii=s.find(ss);
if (ii==s.end()) return false;
auto jj=d.find(ds);
const bool missing=(jj==d.end());
if (if_missing&&!missing) return false;
if (missing) { d[ds]=(*ii).second; return true; } 
//(*jj).second+=StrTy(" , ")+(*ii).second;
(*jj).second+=StrTy(" and ")+(*ii).second;
return true;
}

template <class Tm > 
static bool MutateKvp(Tm & d, const StrTy & ds, const Tm & s, const StrTy & ss, const StrTy & cmd, OutTy & out) 
{
const auto ii=s.find(ss);
if (ii==s.end()) return false;
const StrTy cmdf=cmd; // "("+cmd+" ) ";
const StrTy news=MutateOnly((*ii).second,cmdf,out);
MM_ERR(" MutateOnly F "<<MMPR4((*ii).first,(*ii).second,cmdf,news))
if (news.length()==0) return false;
auto jj=d.find(ds);
if (jj==d.end())
{d[ds]=news; } // MutateOnly((*ii).second,cmd,out);}
else { d[ds]+=StrTy(" , ")+news; } // MutateOnly((*ii).second,cmd,out); } 
return true;
}

// MutateMissingKvp(mej,"year",mej,"date","date +%Y -f-",out);
//  echo 1999 | date +%Y -f-
//date: invalid date ‘1999’
//marchywka@happy:/home/documents/cpp/proj/toobib$ echo 2000 | date +%Y -f-
//2021
//marchywka@happy:/home/documents/cpp/proj/toobib$ 



template <class Tm > 
static bool MutateMissingKvp(Tm & d, const StrTy & ds, const Tm & s, const StrTy & ss, const StrTy & cmd, OutTy & out) 
{
const auto ii=s.find(ss);
if (ii==s.end()) return false;
auto jj=d.find(ds);
if (jj!=d.end()) if((*jj).second!="") return false;
{d[ds]=MutateOnly((*ii).second,cmd,out);
Ss ssx; ssx<<MMPR3(ds,(*ii).second,ss); 
d["X_toobib_Mutatemissing"]+=ssx.str();
}
//else { d[ds]+=StrTy(" , ")+MutateOnly((*ii).second,cmd,out); } 
return true;
}


/////////////////////////////////

typedef typename HierKvp::kvp_lookup_type  KvpLut;
typedef  typename HierKvp::kvp_map_type HierMap;
typedef typename HierKvp::name_list_type HierWords;
typedef  IdxTy  ( *CitFuncp )  (Kvp & , StrTy & , StrTy &,const HierMap &,const IdxTy, const IdxTy );
class kvp_extract {
public:
kvp_extract(const StrTy & s, const StrTy & nm, const CitFuncp p)
: m_string(s),m_name(nm),m_f(p) {}
CitFuncp f() { return m_f;}
Kvp & m() { return m_m; }
const StrTy & name() { return m_name; } 
private:
StrTy m_string,m_name;
Kvp m_m;
CitFuncp m_f;


}; // kvp_extract
static IdxTy find_key_val(StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
const bool itemprop_ok=Bit(flags,0);
bool is_item=false;
IdxTy rc=0;
auto jj=kvp.find("name");
if (jj==kvp.end()) jj=kvp.find("property");
if (jj==kvp.end()) 
{
// see also 
// https://food.ndtv.com/how-to/craving-pizza-try-these-delicious-but-healthy-recipes-this-weekend-2499326
 if( itemprop_ok) jj=kvp.find("itemprop"); // LEDSupply lol
if (jj!=kvp.end()) is_item=true;
}
if (jj==kvp.end()) 
{
// see also 
// https://food.ndtv.com/how-to/craving-pizza-try-these-delicious-but-healthy-recipes-this-weekend-2499326
 if( itemprop_ok) jj=kvp.find("http-equiv"); // LEDSupply lol
if (jj!=kvp.end()) is_item=true;
}

if (jj==kvp.end()) {
rc|=1;
{ 
MM_LOOP(jj,kvp) { MM_ERR("unusual node "<< MMPR3(l,(*jj).first,(*jj).second)) } 
} // dump after final check
} // nothing worked for key..  
else { nkey=(*jj).second; }

jj=kvp.find("content");
if (jj==kvp.end())
{
rc|=2;
{ MM_LOOP(jj,kvp) { MM_ERR("missing value  "<< MMPR3(l,(*jj).first,(*jj).second)) } } // dump
} 
else  // normal kvp
{ 
nval=(*jj).second;
}
if (is_item) MM_ERR(" itemprop "<<MMPR2(nkey,nval))
return rc;
} // find_key_val

static IdxTy munge(Kvp & m, StrTy & nkey, const StrTy & nval,const StrTy & str,const IdxTy l, const IdxTy flags)
{
// TODO this should have a case insensitive and fix common probems.. 
const bool comma=Bit(flags,0);
const bool fix_commons=Bit(flags,8);
//const bool omma=Bit(flags,0);
const char * s= nkey.c_str();
const char * cs=str.c_str(); // "citation_";
const bool zed=(cs[0]==0);
const IdxTy eo= StrUtil::endOf(s,cs);
const StrTy key=StrTy(s+eo);

//MM_MSG(MMPR4(nval,nkey,eo,key))
if (eo||zed)
{
if (m[key]=="") m[key]=nval;
else m[key]+=" and "+nval;
if (fix_commons)
{ // TODO DC.identifier doi code 


} // fix_commons


} // eo
return 0;
} // munge

static IdxTy get_all_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
const StrTy x=nkey;
// TODO FIXME this code causes intermittern crashes with line to long possible. ... 
//if (x==StrTy("citation_reference"))
//{} 
// { m[StrTy("reference")]=StrTy("available deleted for space") ; } 
//else
munge(m,nkey,nval,"",l,0);

return 0;
} // get_all_kvp


static IdxTy get_citation_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
const StrTy x=nkey;
// TODO FIXME this code causes intermittern crashes with line to long possible. ... 
if (x==StrTy("citation_reference"))
//{} 
 { m[StrTy("reference")]=StrTy("available deleted for space") ; } 
else
munge(m,nkey,nval,"citation_",l,0);

return 0;
} // get_citation_kvp

static IdxTy get_cse_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
const StrTy x=nkey;
// TODO FIXME this code causes intermittern crashes with line to long possible. ... 
if (x==StrTy("citation_reference")) {}
//{} 
// { m[StrTy("reference")]=StrTy("available deleted for space") ; } 
//else
munge(m,nkey,nval,"cse_",l,0);

return 0;
} // get_cse_kvp




// https://journals.lww.com/soilsci/Abstract/1974/02000/Effect_of_Pore_Size_on_Diffusion_Coefficients_in.2.aspx
static IdxTy get_wkhealth_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
const StrTy x=nkey;
// TODO FIXME this code causes intermittern crashes with line to long possible. ... 
if (x==StrTy("wkhealth_authors"))
{
const IdxTy len=nval.length();
const char * s=nval.c_str();
char * p= new char [1+(len<<3)];
IdxTy pc=0;
for(IdxTy i=0; i<len; ++i)
{
if (s[i]==';') { 
p[pc]=' '; ++pc;  
p[pc]='a'; ++pc;  
p[pc]='n'; ++pc;  
p[pc]='d'; ++pc;  
p[pc]=' '; ++pc;  
}
else { p[pc]=s[i]; ++pc; }  
} // i 
p[pc]=0;
m["author"]+=StrTy(p);

delete p;
}
else if (x==StrTy("wkhealth_date"))
{
const IdxTy len=nval.length();
const char * s=nval.c_str();
char * p= new char [1+(len<<3)];
IdxTy pc=0;
bool not_good=true;
for(IdxTy i=0; i<len; ++i)
{
if (s[i]==' ') 
{
m["year"]=StrTy(s+i+1);
p[pc]=0;
not_good=false;
const StrTy mo=StrTy(p);
if (mo=="January") m["month"]="1";
else if (mo=="February") m["month"]="2";
else if (mo=="March") m["month"]="3";
else if (mo=="April") m["month"]="4";
else if (mo=="May") m["month"]="5";
else if (mo=="June") m["month"]="6";
else if (mo=="July") m["month"]="7";
else if (mo=="August") m["month"]="8";
else if (mo=="September") m["month"]="9";
else if (mo=="October") m["month"]="10";
else if (mo=="November") m["month"]="11";
else if (mo=="December") m["month"]="12";
m["date"]=nval;
}
else { p[pc]=s[i]; ++pc; }
}
delete p;
if (not_good) munge(m,nkey,nval,"wkhealth_",l,0);
}
else if (x==StrTy("wkhealth_reference"))
//{} 
 { m[StrTy("reference")]=StrTy("available deleted for space") ; } 
else
munge(m,nkey,nval,"wkhealth_",l,0);

return 0;
} // get_wkhealth_kvp



static IdxTy get_DC_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
// need case insentive... 
munge(m,nkey,nval,"DC.",l,0);
munge(m,nkey,nval,"dc.",l,0);
munge(m,nkey,nval,"DCTERMS.",l,0);
munge(m,nkey,nval,"dcterms.",l,0);
const StrTy x=nkey;
// TODO lower case doh!
if ((x==StrTy("dc.identifier")) ||(x==StrTy("DC.identifier")))
{
const char * p=x.c_str();
if (strncmp(p,"DOI:",4)==0)
{
const auto ii=m.find("doi");
// if the compare succeedes this points to string or zero so ok.. 
// but don't take if there is already a doi there... 
if (ii==m.end()) m["doi"]=StrTy(p+4);

}
} // dc.identifier;
return 0;
} // get_citation_kvp


static IdxTy get_og_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
munge(m,nkey,nval,"og:",l,0);

return 0;
} // get_citation_kvp

// https://www.cdc.gov/coronavirus/2019-ncov/covid-data/covidview/index.html
static IdxTy get_cdc_kvp (Kvp & m, StrTy & nkey, StrTy & nval,const HierMap & kvp,const IdxTy l, const IdxTy flags)
{
munge(m,nkey,nval,"cdc:",l,0);
if (nkey=="cdc:last_updated" ) m["date"]=nval;
if (nkey=="cdc:content_source" ) m["author"]=nval;

return 0;
} // get_citation_kvp




///////////////////////////////////////////////////////////////////////


typedef std::vector<kvp_extract> ExtractList;
static IdxTy  ExtractKvp( ExtractList & cfv, const KvpLut & lt,const IdxTy findflags) 
{

IdxTy j=0; 
//MM_MSG(MMPR(lt.size()))
IdxTy l=0; 
const bool dump_nodes=0;
MM_LOOP(ii,lt)
{
if (dump_nodes) { MM_LOOP(jj,(*ii)) { MM_ERR(MMPR3(l,(*jj).first,(*jj).second)) } } // dump
// these  keyts should be lower case earlier... 
//const typename HierKvp::kvp_map_type & kvp=(*ii);
StrTy nkey,nval;
nkey="";
nval="";
const HierMap & kvp=(*ii);
IdxTy fkvrc=find_key_val(nkey,nval,kvp,l,findflags);
if (nkey=="")
{
MM_ERR(" nkey is blank "<<MMPR3(nkey,nval,l))
}
else
{
//IdxTy j=0; 
MM_LOOP(jj,cfv) { (*jj).f()((*jj).m(),nkey,nval,kvp,0,0); ++j; } // jj 
} // nkey not blank 

++l; 
} // ii
j=0;
MM_LOOP(jj,cfv) { j+=(*jj).m().size(); } // jj 

return j; 

} // ExtractKvp 
#if 0 
 static  IdxTy handleadhochtmlOld(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="handleadhochtmlOld";
out.enter(nm); 
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntemph=out.fn("temp",1); //
StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();


//kvp_extract(const StrTy & s, const StrTy & nm, const CitFuncp p)
//std::vector<CitFuncp> cfv;
// std::vector<kvp_extract> cfv;
ExtractList cfv;
// this needs a script version to be fully customizable 
cfv.push_back(kvp_extract("citation_", "citation", get_citation_kvp));
cfv.push_back(kvp_extract("DC.", "DC", get_DC_kvp));
cfv.push_back(kvp_extract("og:", "og", get_og_kvp));
cfv.push_back(kvp_extract("cdc:", "cdc", get_cdc_kvp));
cfv.push_back(kvp_extract("wkhealth_", "wkhealth", get_wkhealth_kvp));
// south china morn post
cfv.push_back(kvp_extract("cse_", "cse", get_cse_kvp));

// load entire document into m.... wow 
//is probably ok for small link data not entire web page lol. 


// get the dom 
// only if missing and debug  
Grc grc=in.getter().normalget(fn,uin,16+32);
const StrTy s1="testHTML --mjmflags hier -mjm \""+fn+"\"";
IdxTy rch=cmd_exec(fntemph,StrTy(),s1,out,8);
// use the more methodical html parser and hierarchial thing...
Ragged r;
r.kvp_splitter(0);
r.load(fntemph); // this should now be parsed with value in the end of line field 
MM_ERR(" loaded "<<MMPR2(fntemph,r.size()))
HierKvp hkvp;
Kvp m,me;
IdxTy lhrc=load_hier_obj(hkvp,&m,r,in,out,0,xflags);
Blob b(hkvp.dump());
b.save(fntemp2);
//MM_MSG(" parsed outpint in "<<MMPR(fntemp2))
//typename HierKvp::name_list_type words;
HierWords words;
//words.push_back("head");
words.push_back("meta");
//m.clear();
KvpLut lt=hkvp.lookup(words,1);
//IdxTy  ekrc=ExtractKvp(  cfv,  lt,0) ;
IdxTy  ekrc=ExtractKvp(  cfv,  lt,1) ;
// TODO this should not get here on file url's
MM_ERR(MMPR2(ekrc,uin))
//if (ekrc==0)
{
ExtractList cfvb;
cfvb.push_back(kvp_extract("all", "all", get_all_kvp));

words.clear();
//words.push_back("body");
words.push_back("meta");
//m.clear();
// doh this is the recurse flag not adj 
KvpLut lt=hkvp.lookup(words,1);
IdxTy  ekrc=ExtractKvp(  cfvb,  lt,1) ;
MM_LOOP(ii,cfvb) { cfv.push_back((*ii)); }
MM_ERR("full body "<< MMPR2(ekrc,uin))
 
} // ekrc

/*
// body ... form ... doh 
// https://www.indianjournals.com/ijor.aspx?target=ijor:ipo&volume=14&issue=2&article=092
85 html 86 body 115 form 841 table 856 tr 857 td 859 table 3532 tr 3533 td 3535 table 3546 tr 3898 td 3904 table 3917 tr 3918 td 3922 div 3926 div 4004 meta 4005 name 4006 text = citation_abstract_html_url
85 html 86 body 115 form 841 table 856 tr 857 td 859 table 3532 tr 3533 td 3535 table 3546 tr 3898 td 3904 table 3917 tr 3918 td 3922 div 3926 div 4009 meta 4010 name 4011 text = citation_pmid
*/
#if 0 
if (false)
{
ExtractList cfvb;
cfvb.push_back(kvp_extract("all", "all", get_all_kvp));

words.clear();
words.push_back("body");
words.push_back("form");
words.push_back("meta");
//m.clear();
// doh this is the recurse flag not adj 
KvpLut lt=hkvp.lookup(words,1);
IdxTy  ekrc=ExtractKvp(  cfvb,  lt,1) ;
MM_LOOP(ii,cfvb) { cfv.push_back((*ii)); }
MM_MSG("full body "<< MMPR2(ekrc,uin))
 
} // ekrc

#endif


// getting DC and citation in wrong order gicing DC as the default or best, 
// http://su.diva-portal.org/smash/record.jsf?dswid=-8098&pid=diva2%3A189901&c=1&searchType=SIMPLE&language=en&query=189901


words.clear();
//MM_MSG(" page title ")
words.push_back("head");
words.push_back("title");
////////////////////////////////////
StrTy pagetitle="";
{
KvpLut lt=hkvp.lookup(words,7);
KvNodes no=hkvp.find(words,7);
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { pagetitle+=(*jj).second; }
} /// ii 
if (pagetitle.length()==0)
{
MM_LOOP(ii,no) {
//me["location"]+=(*ii)->value(); 
StrTy v=(*ii)->value();
pagetitle+=v;
//IdxTy j=0;
//MM_LOOP(jj,cfv) { (*jj).m()["title"]+=v; ++j; } // jj 
}
}
}


//////////////////////////////////


const StrTy url=in.uin();
StrTy ename="";
StrTy type="";

//type=m["X_type"];
type="article"; 
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep= 1; // (type=="NewsArticle");
const bool include_all=true; // all;
if (include_all) me=m;
// callinhg both does things like duplicate names etc...
//IdxTy j=0; 
//IdxTy max=0,loc=0;
// the last entry should be "all" and only want that
// as the best if others lack important fields 
typedef std::map<int,std::vector<IdxTy> > Omap;
Omap szo;
IdxTy jji=0;
// this needs to take the "all" or last one only if 
// none of the othes are good enough... 
// the biggest is at the end... 
MM_LOOP(jj,cfv) { 
// put the "all" at the end... 
int idx=-(int)(*jj).m().size();
MM_ERR(MMPR4((*jj).name(),idx,(*jj).m().size(),jji))
if ( jji==(cfv.size()-1))  idx=0;
szo[idx].push_back(jji); ++jji; 
} 

MM_LOOP(jj,szo)
{
// this no longer works due to "all" set to zero 
const IdxTy sz=-(*jj).first;
const IdxTy n=(*jj).second.size();
MM_LOOP(kk,(*jj).second)
{
auto & ci=cfv[(*kk)];
auto  me=ci.m();
MM_ERR(MMPR4(ci.name(),me.size(),sz,n))
const IdxTy rsz=me.size();
//const IdxTy szii=meii.size();
//MM_MSG(MMPR4(j,(*ii).m().size(),loc,max))
if ( rsz<3) continue; // 
if (sz==0) if (rsz<5) continue; // default should ahve stuff 
// start to fix it up... need page title  
if (me.find("title")==me.end()) me["title"]=pagetitle;
if (me["title"]=="") me["title"]=pagetitle;
me["pagetitle"]=pagetitle;
IdxTy rc= FixAdHoc(in, out, url, nm, ci.name(),ename, type, me, ok_to_keep);

const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 
} // kk 

const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 

} // jj szp 


////////////////////////////////////////////////////////////////////
#if 0 
MM_LOOP(jj,cfv)
{
const IdxTy sz=(*jj).m().size();
if (sz>max) {  loc=j; max=sz; } 
++j; 
} // jj 
j=0;
// now only want "all" to lead if none of the others have several entries
// like author... 
MM_LOOP(ii,cfv)
{
const auto & meii=(*ii).m();
const IdxTy szii=meii.size();
//MM_MSG(MMPR4(j,(*ii).m().size(),loc,max))
if (!all) if (j!=loc) { ++j;  continue;  } 
// this includes fake entries, should fix... 
if (szii==0)  { ++j;  continue; }  
if (szii==1) if ( (*meii.begin()).second=="") if (j!=loc)  { ++j;  continue; }  
me=(*ii).m(); // ma[pi];
if (me.find("title")==me.end()) me["title"]=pagetitle;
if (me["title"]=="") me["title"]=pagetitle;
me["pagetitle"]=pagetitle;
// add fake entreis here...
auto & mej=me; // (*jj).m();
// this is done later 
const auto iid=mej.find("date");
if (iid!=mej.end())
{
if ( ((*iid).second.length())>4)
{
Ss sss;
sss<<MMPR2((*iid).second,(*iid).second.length());
mej["X_toobib_crap"]+=sss.str();
MutateMissingKvp(mej,"day",mej,"date","date +%d -f-",out);
MutateMissingKvp(mej,"month",mej,"date","date +%m -f-",out);
MutateMissingKvp(mej,"year",mej,"date","date +%Y -f-",out);
} // > 4

} // iid != end

MutateMissingKvp(mej,"abstract",mej,"description","cat -",out);
MutateMissingKvp(mej,"author",mej,"creator","cat -",out);
MutateMissingKvp(mej,"journal",mej,"journal_title","cat -",out);

const StrTy fnbibnew=out.fn("bibtex",1); //
out.enter((*ii).name());
IdxTy rc= synthesize(in, out, url, nm, fnbibnew, ename, type, me, ok_to_keep?0:1);
out.exit((*ii).name());
MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))
const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 

++j;
}
/////////////////////////////////////////////////////////////
#endif
//return rc; jjj
// note that headless may fail where wget/ua=Firefox works :)

out.exit(nm); 
return 0;

} // handleadhochtmlOld
#endif

/////////////////////////////////////////////////


 static  IdxTy handleadhochtml(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="handleadhochtml";
out.enter(nm); 
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const bool skip_dc=Bit(xflags,0);

//kvp_extract(const StrTy & s, const StrTy & nm, const CitFuncp p)
//std::vector<CitFuncp> cfv;
// std::vector<kvp_extract> cfv;
ExtractList cfv;
// this needs a script version to be fully customizable 
cfv.push_back(kvp_extract("citation_", "citation", get_citation_kvp));
if (!skip_dc) { cfv.push_back(kvp_extract("DC.", "DC", get_DC_kvp)); }
cfv.push_back(kvp_extract("og:", "og", get_og_kvp));
cfv.push_back(kvp_extract("cdc:", "cdc", get_cdc_kvp));
cfv.push_back(kvp_extract("wkhealth_", "wkhealth", get_wkhealth_kvp));
// south china morn post
cfv.push_back(kvp_extract("cse_", "cse", get_cse_kvp));

// load entire document into m.... wow 
//is probably ok for small link data not entire web page lol. 


// get the dom 
// only if missing and debug  
Grc grc=in.getter().normalget(fn,uin,16+32);
Ragged r;
if (true)
{
typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
// bit 4 skip json... 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const IdxTy fnty=256+16+8+1; //  not include eq, deb  , only good json 
Ragged & h=r;
Ragged j; // should have dummy 
IdxTy rc=bt.parse_both(h,j,fn,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))
} // true
if (false)
{
StrTy fntemph=out.fn("temp",1); //
const StrTy s1="testHTML --mjmflags hier -mjm \""+fn+"\"";
IdxTy rch=cmd_exec(fntemph,StrTy(),s1,out,8);
// use the more methodical html parser and hierarchial thing...
r.kvp_splitter(0);
r.load(fntemph); // this should now be parsed with value in the end of line field 
MM_ERR(" loaded "<<MMPR2(fntemph,r.size()))

} // false




HierKvp hkvp;
Kvp m,me,m2;
IdxTy lhrc=load_hier_obj(hkvp,&m,r,in,out,0,xflags);
Blob b(hkvp.dump());
b.save(fntemp2);
//MM_MSG(" parsed outpint in "<<MMPR(fntemp2))
//typename HierKvp::name_list_type words;
HierWords words;
//words.push_back("head");
words.push_back("meta");
//m.clear();
KvpLut lt=hkvp.lookup(words,1);
//IdxTy  ekrc=ExtractKvp(  cfv,  lt,0) ;
IdxTy  ekrc=ExtractKvp(  cfv,  lt,1) ;
// TODO this should not get here on file url's
MM_ERR(MMPR2(ekrc,uin))
//if (ekrc==0)
{
ExtractList cfvb;
cfvb.push_back(kvp_extract("all", "all", get_all_kvp));

words.clear();
//words.push_back("body");
words.push_back("meta");
//m.clear();
// doh this is the recurse flag not adj 
// this is a verctor of maps... 
KvpLut lt=hkvp.lookup(words,1);
IdxTy  ekrc=ExtractKvp(  cfvb,  lt,1) ;
MM_LOOP(ii,cfvb) { cfv.push_back((*ii)); }
MM_ERR("ALL full body "<< MMPR3(cfvb.size(),ekrc,uin))
 
} // ekrc

/*
// body ... form ... doh 
// https://www.indianjournals.com/ijor.aspx?target=ijor:ipo&volume=14&issue=2&article=092
85 html 86 body 115 form 841 table 856 tr 857 td 859 table 3532 tr 3533 td 3535 table 3546 tr 3898 td 3904 table 3917 tr 3918 td 3922 div 3926 div 4004 meta 4005 name 4006 text = citation_abstract_html_url
85 html 86 body 115 form 841 table 856 tr 857 td 859 table 3532 tr 3533 td 3535 table 3546 tr 3898 td 3904 table 3917 tr 3918 td 3922 div 3926 div 4009 meta 4010 name 4011 text = citation_pmid
*/
#if 0 
if (false)
{
ExtractList cfvb;
cfvb.push_back(kvp_extract("all", "all", get_all_kvp));

words.clear();
words.push_back("body");
words.push_back("form");
words.push_back("meta");
//m.clear();
// doh this is the recurse flag not adj 
KvpLut lt=hkvp.lookup(words,1);
IdxTy  ekrc=ExtractKvp(  cfvb,  lt,1) ;
MM_LOOP(ii,cfvb) { cfv.push_back((*ii)); }
MM_MSG("full body "<< MMPR2(ekrc,uin))
 
} // ekrc

#endif


// getting DC and citation in wrong order gicing DC as the default or best, 
// http://su.diva-portal.org/smash/record.jsf?dswid=-8098&pid=diva2%3A189901&c=1&searchType=SIMPLE&language=en&query=189901
// these do not work... 
StrTy pagetitle="";
StrTy published_date="";
IdxTy ptrc=GetPageTitle(pagetitle,in,out,r,m2,hkvp,xflags);
IdxTy ptrcpd=GetPublishedDate(published_date,"published_date",in,out,r,m2,hkvp,xflags);
if (pagetitle=="") GetPublishedDate(pagetitle,"title",in,out,r,m2,hkvp,xflags);
MM_ERR(MMPR2(pagetitle,published_date))


if (false)
{
words.clear();
//MM_MSG(" page title ")
words.push_back("head");
words.push_back("title");
////////////////////////////////////
{
KvpLut lt=hkvp.lookup(words,7);
KvNodes no=hkvp.find(words,7);
MM_LOOP(ii,lt)
{
MM_LOOP(jj,(*ii)) { pagetitle+=(*jj).second; }
} /// ii 
if (pagetitle.length()==0)
{
MM_LOOP(ii,no) {
//me["location"]+=(*ii)->value(); 
StrTy v=(*ii)->value();
pagetitle+=v;
//IdxTy j=0;
//MM_LOOP(jj,cfv) { (*jj).m()["title"]+=v; ++j; } // jj 
}
}
}
} // scope

//////////////////////////////////
if(published_date=="") 
{
MM_LOOP(jj,cfv) { 

const auto & m=(*jj).m();
const auto ii= m.find("published_date");
if (ii!=m.end()) {  published_date=(*ii).second; } 

}// jj
}
MM_ERR(MMPR2(pagetitle,published_date))

const StrTy url=in.uin();
StrTy ename="";
StrTy type="";

//type=m["X_type"];
type="article"; 
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep= 1; // (type=="NewsArticle");
const bool include_all=true; // all;
if (include_all) me=m;
// callinhg both does things like duplicate names etc...
//IdxTy j=0; 
//IdxTy max=0,loc=0;
// the last entry should be "all" and only want that
// as the best if others lack important fields 
typedef std::map<int,std::vector<IdxTy> > Omap;
Omap szo;
IdxTy jji=0;
// this needs to take the "all" or last one only if 
// none of the othes are good enough... 
// the biggest is at the end... 
MM_LOOP(jj,cfv) { 
// put the "all" at the end... 
const auto & foo=(*jj).m();
int idx=-(int)(*jj).m().size();
// it is picking stupid things, try to make author and title work more
auto ii=foo.find("author");
if (ii!=foo.end()) if ((*ii).second.length()>2) idx-=2;
ii=foo.find("date");
if (ii!=foo.end()) if ((*ii).second.length()>2) idx-=2;
Ss ss; MM_LOOP(ii,foo ) { ss<<(*ii).first<<" "; }
MM_ERR(MMPR4((*jj).name(),idx,(*jj).m().size(),jji)<<MMPR(ss.str()))
if ( jji==(cfv.size()-1))  idx=0;
szo[idx].push_back(jji); ++jji; 
} 

MM_LOOP(jj,szo)
{
//auto jjend=cfv.rbegin();
//auto xxx=(*szo.rbegin()).second;
//MM_ERR(MMPR(xxx.size()));
//if (xxx.size()==0) xxx.push_back(0);
//auto& jjend=cfv[xxx[xxx.size()-1]];
//auto mxend=(jjend).m();
// this no longer works due to "all" set to zero 
const IdxTy sz=-(*jj).first;
const IdxTy n=(*jj).second.size();
MM_LOOP(kk,(*jj).second)
{
auto & ci=cfv[(*kk)];
auto  me=ci.m();
MM_ERR(MMPR4(ci.name(),me.size(),sz,n))
const IdxTy rsz=me.size();
//const IdxTy szii=meii.size();
//MM_MSG(MMPR4(j,(*ii).m().size(),loc,max))
if ( rsz<3) continue; // 
if (sz==0) if (rsz<5) continue; // default should ahve stuff 
// start to fix it up... need page title  
if (me.find("title")==me.end()) me["title"]=pagetitle;
if (me.find("date")==me.end()) me["date"]=published_date;
//if (me.find("date")==me.end()) me["date"]=mxend["published_date"];
if (me["title"]=="") me["title"]=pagetitle;
me["pagetitle"]=pagetitle;
IdxTy rc= FixAdHoc(in, out, url, nm, ci.name(),ename, type, me, ok_to_keep);

const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 
} // kk 

const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 

} // jj szp 


////////////////////////////////////////////////////////////////////
#if 0 
MM_LOOP(jj,cfv)
{
const IdxTy sz=(*jj).m().size();
if (sz>max) {  loc=j; max=sz; } 
++j; 
} // jj 
j=0;
// now only want "all" to lead if none of the others have several entries
// like author... 
MM_LOOP(ii,cfv)
{
const auto & meii=(*ii).m();
const IdxTy szii=meii.size();
//MM_MSG(MMPR4(j,(*ii).m().size(),loc,max))
if (!all) if (j!=loc) { ++j;  continue;  } 
// this includes fake entries, should fix... 
if (szii==0)  { ++j;  continue; }  
if (szii==1) if ( (*meii.begin()).second=="") if (j!=loc)  { ++j;  continue; }  
me=(*ii).m(); // ma[pi];
if (me.find("title")==me.end()) me["title"]=pagetitle;
if (me["title"]=="") me["title"]=pagetitle;
me["pagetitle"]=pagetitle;
// add fake entreis here...
auto & mej=me; // (*jj).m();
// this is done later 
const auto iid=mej.find("date");
if (iid!=mej.end())
{
if ( ((*iid).second.length())>4)
{
Ss sss;
sss<<MMPR2((*iid).second,(*iid).second.length());
mej["X_toobib_crap"]+=sss.str();
MutateMissingKvp(mej,"day",mej,"date","date +%d -f-",out);
MutateMissingKvp(mej,"month",mej,"date","date +%m -f-",out);
MutateMissingKvp(mej,"year",mej,"date","date +%Y -f-",out);
} // > 4

} // iid != end

MutateMissingKvp(mej,"abstract",mej,"description","cat -",out);
MutateMissingKvp(mej,"author",mej,"creator","cat -",out);
MutateMissingKvp(mej,"journal",mej,"journal_title","cat -",out);

const StrTy fnbibnew=out.fn("bibtex",1); //
out.enter((*ii).name());
IdxTy rc= synthesize(in, out, url, nm, fnbibnew, ename, type, me, ok_to_keep?0:1);
out.exit((*ii).name());
MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))
const IdxTy nnow=out.found();
if (!all) if (nnow>nstart) break; 

++j;
}
/////////////////////////////////////////////////////////////
#endif
//return rc; jjj
// note that headless may fail where wget/ua=Firefox works :)

out.exit(nm); 
return 0;

} // handleadhochtml




////////////////////////////////////////////////////
static IdxTy FixAdHoc(const InpTy &in, OutTy & out, const IdioRL & url, const StrTy & nm, const StrTy & cn, const StrTy& ename, const StrTy & type, Kvp & mej, const bool ok_to_keep)
{
// this should be in one place, split with pagetitle in caller...
// This  is not working wtf... 
MutateMissingKvp(mej,"abstract",mej,"description","cat -",out);
MutateMissingKvp(mej,"author",mej,"creator","cat -",out);
MutateMissingKvp(mej,"author",mej,"creator.personalname","cat -",out);
MutateMissingKvp(mej,"doi",mej,"identifier.doi","cat -",out);
MutateMissingKvp(mej,"page",mej,"identifier.pagenumber","cat -",out);
MutateMissingKvp(mej,"journal",mej,"journal_title","cat -",out);

const StrTy fnbibnew=out.fn("bibtex",1); //
//out.enter(ci.name());
out.enter(cn);
IdxTy rc= synthesize(in, out, url, nm, fnbibnew, ename, type, mej, ok_to_keep?0:1);
out.exit(cn);

MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))
return 0; 
}// FixAdHoc



////////////////////////////////////////////////////////////////////////

/*
the date defaults assumptions are stupid, 
echo 1999 | date +%Y -f-
date: invalid date ‘1999’
marchywka@happy:/home/documents/cpp/proj/toobib$ echo 2000 | date +%Y -f-
2021
marchywka@happy:/home/documents/cpp/proj/toobib$ 
*/


////////////////////////////////////////////////////////////////////
#if 0 

 static  IdxTy handleadhochtml(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="handleadhochtml";
out.enter(nm); 
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fnme=out.fn("me"); //
StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntemph=out.fn("temp",1); //
StrTy fntempg=out.fn("temp",1); //
StrTy fntempj=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// only if missing and debug  
Grc grc=in.getter().normalget(fn,uin,16+32);
const StrTy s1="testHTML --mjmflags hier -mjm \""+fn+"\"";
IdxTy rch=cmd_exec(fntemph,StrTy(),s1,out,8);
// use the more methodical html parser and hierarchial thing...
Ragged r;
r.kvp_splitter(0);
r.load(fntemph); // this should now be parsed with value in the end of line field 
HierKvp hkvp;
Kvp m,me;
//kvp_extract(const StrTy & s, const StrTy & nm, const CitFuncp p)
//std::vector<CitFuncp> cfv;
std::vector<kvp_extract> cfv;
cfv.push_back(kvp_extract("citation_", "citation", get_citation_kvp));
cfv.push_back(kvp_extract("DC.", "DC", get_DC_kvp));
cfv.push_back(kvp_extract("og:", "og", get_og_kvp));
//cfv.push_back(get_DC_kvp);
//std::vector<StrTy> pfxx;
//pfxx.push_back("citation_");
//pfxx.push_back("DC.");
//pfxx.push_back("og:");
//Kvp ma[pfxx.size()];
MM_ERR(" loaded "<<MMPR2(fntemph,r.size()))
// load entire document into m.... wow 
//is probably ok for small link data not entire web page lol. 
IdxTy lhrc=load_hier_obj(hkvp,&m,r,in,out,0,xflags);
Blob b(hkvp.dump());
b.save(fntemp2);
MM_MSG(" parsed outpint in "<<MMPR(fntemp2))
//typename HierKvp::name_list_type words;
HierWords words;
words.push_back("head");
words.push_back("meta");
//m.clear();
KvpLut lt=hkvp.lookup(words,1);
MM_MSG(MMPR(lt.size()))
IdxTy l=0; 
const bool dump_nodes=0;
MM_LOOP(ii,lt)
{
if (dump_nodes) { MM_LOOP(jj,(*ii)) { MM_MSG(MMPR3(l,(*jj).first,(*jj).second)) } } // dump
// these  keyts should be lower case earlier... 
//const typename HierKvp::kvp_map_type & kvp=(*ii);
StrTy nkey,nval;
const HierMap & kvp=(*ii);
IdxTy fkvrc=find_key_val(nkey,nval,kvp,l,0);
IdxTy j=0; 
MM_LOOP(jj,cfv)
{
//(*jj)(ma[j],nkey,nval,kvp,0,0);
//out.enter((*jj).name());
(*jj).f()((*jj).m(),nkey,nval,kvp,0,0);
//out.exit((*jj).name());
++j; 
} // jj 
#if 0
auto jj=kvp.find("name");
if (jj==kvp.end()) jj=kvp.find("property");
if (jj==kvp.end())
{
{ MM_LOOP(jj,(*ii)) { MM_MSG("unusual node "<< MMPR3(l,(*jj).first,(*jj).second)) } } // dump
} 
else
{
const StrTy name=(*jj).second;
jj=kvp.find("content");
if (jj==kvp.end())
{
{ MM_LOOP(jj,(*ii)) { MM_MSG("missing value  "<< MMPR3(l,(*jj).first,(*jj).second)) } } // dump
} 
else  // normal kvp
{ 
const StrTy val=(*jj).second;
const char * s= name.c_str();

for(IdxTy pi=0; pi<pfxx.size(); ++pi )
{
const char * cs=pfxx[pi].c_str(); // "citation_";
const IdxTy eo= StrUtil::endOf(s,cs);
const StrTy key=StrTy(s+eo);
MM_MSG(MMPR4(val,name,eo,key))
if (eo)
{
auto & mm=ma[pi];
if (mm[key]=="") mm[key]=val;
else mm[key]+=" and "+val;
} // eo
} // pi
} // normal kvp  
} // kvp has name 
#endif

++l; 
} // ii

const StrTy url=in.uin();
StrTy ename="";
StrTy type="";

type=m["X_type"];
type="article"; 
// this will try to add it to the list of found bibs
// which is not good yet so set bit 0 of the flgs
// to stop adding ... 
const bool ok_to_keep= 1; // (type=="NewsArticle");
const bool include_all=true; // all;
if (include_all) me=m;
// callinhg both does things like duplicate names etc...
IdxTy j=0; 
IdxTy max=0,loc=0;
MM_LOOP(jj,cfv)
{
//(*jj)(ma[j],nkey,nval,kvp,0,0);
//out.enter((*jj).name());
const IdxTy sz=(*jj).m().size();
if (sz>max) {  loc=j; max=sz; } 
//out.exit((*jj).name());
++j; 
} // jj 

//IdxTy na=NewsArticleKvpMap(type,me,m,out,0);
//if (m["X_type"]=="NewsArticle") { type="article"; me["subtype"]="news"; }
///IdxTy lip=LiPersonKvpMap(type,me,m,r,out,0);

//for(IdxTy pi=0; pi<pfxx.size(); ++pi)
//for(IdxTy pi=0; pi<pfxx.size(); ++pi)
j=0;
MM_LOOP(ii,cfv)
{
if (!all) if (j!=loc) continue; 
if ((*ii).m().size()==0) if (j!=loc)  continue; 
me=(*ii).m(); // ma[pi];
const StrTy fnbibnew=out.fn("bibtex",1); //
out.enter((*ii).name());
IdxTy rc= synthesize(in, out, url, nm, fnbibnew, ename, type, me, ok_to_keep?0:1);
out.exit((*ii).name());
MM_ERR("should have saved to "<<MMPR2(fnbibnew,type))
const IdxTy nnow=out.found();
++j;
}

//return rc; 
// note that headless may fail where wget/ua=Firefox works :)

out.exit(nm); 
return 0;

} // handleadhochtml
////////////////////////////////////////////////////////////////////////


#endif






////////////////////////////
 static  IdxTy headlessunsafe(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
StrTy type;
const StrTy fn=out.fn(); // in.ZZ
// only if not already have 
Grc grc=in.getter().normalget(fn,in.uin(),16+32);
IdxTy rci= identify(type, fn,"",  out, 0);
if (type==StrTy("pdf")||(rci!=0))
{
MM_ERR( MMPR2(rci,out.nm()) << "skipping as "<<MMPR2(fn, type)<<" is pdf ")
 return 1;
}
return 0; 
}
 static  IdxTy handleheadless(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="handleheadless";
out.enter(nm); 
{
// this should be ok now with the chromate getter 
if (false)if (headlessunsafe(in,out,xflags)) 
{
//MM_ERR( MMPR2(rci,nm) << "skipping as "<<MMPR2(fn, type)<<" is pdf ")
out.exit(nm); 
 return 1;
}

}
const StrTy & uin=in.uin();
if (!in.uin().http())
{
MM_ERR(" skipp  because uin is not http  "<<MMPR2(nm,uin))
out.exit(nm); 
return 0; 
}
out.push_source("headless");
const StrTy fn=out.fn(); // in.ZZ
StrTy x;
// overwrite the normal get file...
// this should now be print-to-pdf with bit zero of the flags clear
Grc grc=in.getter().headlessget(fn,uin,0);

// TODO FIXME This can now recurse with the file as the url... 
IdxTy rt=Recurse(in,out,fn);

out.pop_source();
out.exit(nm); 
return rt;
}
