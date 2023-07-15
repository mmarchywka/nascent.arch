 



/////////////////////////////////////////////////////
// move to guessesh

 static  IdxTy FetchAsRis(const InpTy & in , OutTy & out, const StrTy & url, const StrTy & nm, const IdxTy flags )  { 
const StrTy fnbib=out.fn(nm); //
const StrTy fntemp=out.fn("temp"); //
Grc grc=in.getter().normalget(fntemp,url,16);
// this should now contain ris, need to xform . 
MM_ERR(MMPR2(fnbib,fntemp))
 
BibEntry bib;
Xrc xrc=in.xformer().ris_to_bib(bib,fntemp,0);
if (bib.name()=="") bib.make_name();
MM_ERR(MMPR(bib.format()))
const StrTy bibtext=(bib.format());
/////////////////////////////////////////
{
//Blob b;
//static IdxTy CleanUpBadBib(Blob & b, const StrTy & fntemp ,  const StrTy & fnbib, const IdxTy flags) 
//const IdxTy cubbrc=CleanUpBadBib( b, fntemp , fnbib, 0) ;
//BibEntry bib2;
//bib2.parse(
}


/////////////////////////////////
if (!false) { 
Blob b(bibtext);
b.save(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
#if 0 
Frc frc=in.finisher().parse(b,fnbib,"",0);
Frc frcb=0; // in.finisher().validate(b,fnbib,"",0);
MM_ERR(MMPR(IdxTy(frc)))
// keep anyway for edit
//if (frc==0)
{
// if b and fnbib both have contents both written doh... 
out.result( b, fnbib, url, nm,  0);
} // frc==0
#endif
}
return 0; 
}

 static  IdxTy FetchAsEnw(const InpTy & in , OutTy & out, const StrTy & url, const StrTy & nm, const IdxTy flags )  { 
const StrTy fnbib=out.fn(nm); //
Grc grc=in.getter().normalget(fnbib,url,16);
// this should now contain ris, need to xform . 
MM_ERR(MMPR(fnbib))
 
BibEntry bib;
Xrc xrc=in.xformer().enw_to_bib(bib,fnbib,0);
MM_ERR(MMPR(bib.format()))
const StrTy bibtext=(bib.format());
/////////////////////////////////
if (!false) { 
Blob b(bibtext);
b.save(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
}
return 0; 
} // FetchAsEnw

// this should fix bad line break and dos line endings. 
// fntemp is read-only as-received for debug 
static IdxTy CleanUpBadBib(Blob & b, const StrTy & fntemp ,  const StrTy & fnbib, const IdxTy flags) 
{
BibEntryMap bem;
b.load(fntemp);
b.nocrlf();
b.save(fnbib);
bem.parse(fnbib);
if (bem.size()!=1)
{
//MM_ERR(MMPR(nm)<< " fails "<<MMPR(bem.size())<<MMPR2(fn,fnbib))
MM_ERR(" parsing fails "<<MMPR(bem.size())<<MMPR2(fntemp,fnbib))
//out.exit(nm);
return 2;
}
BibEntry  be=*(bem.begin()); // ((*(bem.begin())).second)[0];
//be.make_name();
//Blob b;
if (be.name()=="") be.make_name();
b=be.format();
b.save(fnbib);
return 0;
} // CleanUpBadBib

// this could also use the partse... 
// but there is other code to insert missing srcurl and citeurl etc 
// use insert_missing in out... 
static IdxTy inserturl(const StrTy & fnbib, const StrTy & url, OutTy & out, const IdxTy flags)
{
// std::map<StrTy, StrTy> m;
// m["url"]=url;
// out.insert_missing(fnbib,m,0);
StrTy iucmd="awk -v v=\""+url+"\"  'BEGIN{have=0; lvl=0; done=0}{if ($1==\"url\") have=1;   split($0, chars, \"\"); out=\"\";  ln=length($0); for (i=1; i <= ln; i++) { c=chars[i]; if (c==\"{\") { lvl=lvl+1;if (lvl>max) max=lvl; }  if ( c==\"}\") { lvl=lvl-1; if (lvl==0) if (max>0){if (have==0)  {out=out\" , url={\"v\"}}\";}else { out=out\"}\";  } if ( done==0) { print out ; }  done=1; } }  if ( done==0) {out=out\"\"c; }   } if ( done==0) { if ( ( max==0) ||(lvl>0)) print out;} }'";
const StrTy fntemp=out.fn("temp");
const IdxTy rc=filter_small_file(fntemp,fnbib,iucmd,out,0);
StrTy x;
const StrTy cmdmo="mv \""+fntemp+"\" \""+fnbib+"\"";
cmd_exec(x,"",cmdmo,out,0);
return 0;
} // inserturl 

 static  const StrTy &  GuessConcat(const InpTy & in , OutTy & out )  
{

static const StrTy guessconcat="awk 'BEGIN{c=0}{ if (c==0) { if (($3!=\"{\")||( $4 != \"\"))  {print $0;} else {c=1; line=$0; }} else {if ( substr($0,0,1)==\"}\") {print line$0; line=\"\";c=0;} else {line=line $0}}}END { if (line!=\"\") print line; } '";

return guessconcat;
}

/////////////////////////////////////////////////// 
// redundant with get_xref in handlers header... 
 static  IdxTy guessdoi(const InpTy & in , OutTy & out , const IdxTy xflags =0  )  { 

const StrTy nm="guesssdoi";
out.enter(nm);
const StrTy fnbib=out.fn("bibtex"); // in.ZZ
const StrTy & _uin=in.uin();
/// may be a doi link although should already be handled... 
//  "http://doi.org/10.1523/jneurosci.0519-17.2017" 
// "https://doi.org/10.5438/0012",
// dx.doi.org/10.22541/au.159103680.00306295
const char * p =_uin.c_str();
const IdxTy pos=mjm_strings::indexOf(p,"10.");
if (pos==(~0))
{
MM_ERR(MMPR4(nm,__FUNCTION__, pos,_uin))
return 0; 
}
const StrTy uin=StrTy(p+pos);
MM_ERRF(MMPR4(nm,uin, pos,_uin))
Grc grc=in.getter().crossref(fnbib,uin,0);
const StrTy src=grc.name();
const StrTy srcx=(src.length())?(StrTy("(")+src+")"):StrTy();
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grc.url(),nm+srcx,0);
out.exit(nm);
return 0;
}
// guessdoi
 static  IdxTy guesscell(const InpTy & in , OutTy & out , const IdxTy xflags =0  )  { 
const StrTy nm="guessscell";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn("normalget"); // in.ZZ
Grc grc=in.getter().normalget(fn,uin,16);
StrTy url;
const StrTy cmdin="cat \""+fn+"\"| grep downloadFile ";
cmd_exec(url,"",cmdin,out,0);
const StrTy fnbib=out.fn(nm); //
MM_ERR(" "<<MMPR3(nm,uin,url))
Grc grc2=in.getter().normalget(fnbib,url,16);
Blob b;
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0; 
#if 0
guesscell()
uin="$1"
#dff=`lynx -dump -source "$uin" | grep downloadFile `
dff=`$WGET --no-check-certificate -O - -q  "$uin" | grep downloadFile `
echo guesscell think we want $dff | smsg
#This worked, but it wanted a post 
# http://www.cell.com/cell/abstract/S0092-8674(16)30213-6
# $WGET -O xxx -S -v "https://secure.jbs.elsevierhealth.com/action/downloadCitationSecure?code=cell-site&direct=true&downloadFile=marlin_cell165_579&objectUri=pii:S0092867416302136"
url="https://secure.jbs.elsevierhealth.com/action/downloadCitationSecure?code=cell-site&direct=true&downloadFile=marlin_cell165_579&objectUri=pii:S0092867416302136"
url="$dff"
#endif
} // guesscell


 static  IdxTy guessaip(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guessaip";
out.enter(nm);
const StrTy & uin=in.uin();
//const StrTy fn=out.fn("normalget"); // in.ZZ
//Grc grc=in.getter().normalget(fn,uin,16);
StrTy url=uin+"/cite/bibtex";
const StrTy fnbib=out.fn(nm); //
MM_ERR(" "<<MMPR3(nm,uin,url))
Grc grc=in.getter().normalget(fnbib,url,16);
Blob b;
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0; 
#if 0
guessaip()
uin="$1"
url="$uin/cite/bibtex"
guessget "$url" "$uin"

#endif
} // guessaip

#if 0
guesshind enw works but the things has special leading chars
and they do support ris , 
%0 Journal Article
%T Targeting Mitochondria as Therapeutic Strategy for Metabolic Disorders
%J The Scientific World Journal
%V 2014
%A Sorriento, Daniela
%A Pascale, Antonietta Valeria
%A Finelli, Rosa
%A Carillo, Anna Lisa
%A Annunziata, Roberto
%A Trimarco, Bruno
%A Iaccarino, Guido
%R 10.1155/2014/604685
%D 2014
%U http://dx.doi.org/10.1155/2014/604685
%] 604685
%P 9
marchywka@happy:/home/documents/cpp/proj/toobib/junk$ cat xxx | od -ax | more
0000000   o   ;   ?   %   0  sp   J   o   u   r   n   a   l  sp   A   r
           bbef    25bf    2030    6f4a    7275    616e    206c    7241
0000020   t   i   c   l   e  nl   %   T  sp   T   a   r   g   e   t   i
           6974    6c63    0a65    5425    5420    7261    6567    6974



#endif

 static  IdxTy guesshind(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
// this uses enw, may be better to use doi...
const StrTy nm="guesshind";
out.enter(nm);
const StrTy & uin=in.uin();
//const StrTy fn=out.fn("normalget"); // in.ZZ
//Grc grc=in.getter().normalget(fn,uin,16);
StrTy url;
IdxTy rcfar=0;
if (false)
{
const StrTy cmdin="sed -e 's/www.hindawi/files.hindawi/' | sed -e 's/\\/[a-z]*\\/$//'| sed -e 's;/$;;'"; // .enw
cmd_exec(url,uin,cmdin,out,0);
url=url+".enw";
// the doi crossref result is usally better.. 
rcfar= FetchAsEnw( in , out,  url,  nm, 0) ;
}
else
{
const StrTy cmdin=" sed -e 's/\\/[a-z]*\\/$//'| sed -e 's;/$;;'"; // .enw
cmd_exec(url,uin,cmdin,out,0);
url=url+".ris";
rcfar= FetchAsRis( in , out,  url,  nm, 0) ;
} 

/*const StrTy fnbib=out.fn(nm); //
MM_ERR(" "<<MMPR3(nm,uin,url))
Grc grc=in.getter().normalget(fnbib,url,16);
Blob b;
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

*/
out.exit(nm);
//MM_ERR(" not impl")  
return 0; 
} // guesshind


#if 0
nm="guesshind"
uin="$1"
fn="$temp1"
fndest="$temp5"
#http://www.hindawi.com/journals/tswj/2014/604685/cta/
#http://files.hindawi.com/journals/tswj/2014/604685.enw
#url=`echo $uin | sed -e 's/www.hindawi/files.hindawi/' | sed -e 's/\/cta\//.enw/'`
url=`echo $uin | sed -e 's/www.hindawi/files.hindawi/' | sed -e 's/\/[a-z]*\/$//'| sed -e 's;/$;;'`.enw
#guessget "$url" "$uin"
# this is worse than the doi version, this should just be skipped
if [ 1 == 1 ]
then 
echo $nm rose than doi skipping | smsg
else
normalget "$fn" "$url"
cat "$fn" | $ENW > "$fndest"




#endif

// TODO euroomc is not getting here  - guessEpmc
 static  IdxTy guesspmc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesspmc";
out.enter(nm);
const StrTy & uin=in.uin();
//num=`echo "$uin" | sed -e 's;/;\n;g'| grep "[0-9]" | grep PMC `
//citex="\\cite{$num}"
const StrTy cmd="sed -e 's;/;\\n;g'| grep \"[0-9]\" | grep PMC ";
StrTy pmcid;
cmd_exec(pmcid,uin,cmd,out,0);
// TODO WTHAT THE F this needs to check punt firs tfucj 
const bool get_anyway= !in.punt_pmc(); // !false;
MM_ERR(__FUNCTION__<<MMPR2(pmcid,get_anyway))
if (!get_anyway) if ( pmcid.length()==0) { out.exit(nm); return 0;  }

if (get_anyway)
{
const StrTy db="pmc";
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
//AccessTy pubmed(const StrTy & fn, const StrTy& id, const StrTy& db, const IdxTy flags)
Grc grc=in.getter().pubmed(fntemp,pmcid,db,0);
BibEntry be;
Xrc xrc=in.xformer().pmc_to_bib(be,fntemp,0);
be.set("url",uin);
Blob b;
b=be.format();
b.save(fnbib);
MM_ERR(__FUNCTION__<<MMPR3(pmcid,fnbib,fntemp))
StrTy url=grc.url();
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0;
} // get_anyway
out.cite(pmcid); 
 out.punt(); 
 MM_ERR("  punting "<<MMPR2(nm,pmcid))  

return 0; 
}

 //static  IdxTy guesspmcafick(const InpTy & in , OutTy & out ) 
 static  IdxTy guesspmceuro(const InpTy & in , OutTy & out , const IdxTy xflags=0) 

{ 
const StrTy nm="guesseuropmc";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy cmd=" sed -e 's;.*PMC/;PMC;'|sed -e 's;/.*;;' |  grep \"[0-9]\" | grep PMC ";
StrTy pmcid;
// TODO this may include a crlf 
cmd_exec(pmcid,uin,cmd,out,0);
if ( pmcid.length()!=0) {  out.punt(); out.cite(pmcid);  }
 MM_ERR("  punting "<<MMPR2(nm,pmcid))  
#if 0 
# https://europepmc.org/article/PMC/6174679
# these  pompous  hide the data in 
#  javascript or other  
guesspmcafick()
uin="$1"
url="$1"
# this is picking up fragments etc 
#num=`echo "$uin" | sed -e 's;/;\n;g'| grep "[0-9]"`
#num=`echo "$uin" | sed -e 's;/;\n;g'| grep "[0-9]" | grep PMC `
num=`echo "$uin" | sed -e 's;.*PMC/;PMC;'|sed -e 's;/.*;;' |  grep "[0-9]" | grep PMC `
citex="\\cite{$num}"
# note that this COULD be PMC not PMID 
#echo "\\\\cite{PMID$num}" | my_paste
echo "$citex" | my_paste
#endif
//MM_ERR(" not impl")  
out.exit(nm);
return 0; 
}

// What the is this 
 static  IdxTy guessepmc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guessepmc";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url=uin+"?javascript_support=no";
const StrTy temp=out.fn("temp");
Grc grc=in.getter().normalget(temp,url,0);
Blob b;
Xrc xrc=in.xformer().html_to_rendered(b,temp,0);
StrTy s=StrTy(b);
StrTy cmdin="grep -i doi | grep -i http | tail -n 1 |  sed -e 's/.*doi.org.//' | awk '{print $NF}' ";
StrTy doi;
// leaving  crlf 
//cmd_exec(doi,s,cmdin,out,0);
cmd_exec(doi,s,cmdin,out,4); // crlf 
MM_ERR(MMPR2(nm,doi))
IdxTy rcx=get_xref(out,in,doi,nm,0);
#if 0
nm=guessepmc
uin="$1"
fn="$temp1"
normalget "$fn" "$uin?javascript_support=no"
#doi=`lynx -dump --force-html "$fn" | grep -i doi | grep -i Cite | head -n 1 |  sed -e 's/.*]//' | awk '{print $NF}' `
doi=`lynxhtml "$fn" | grep -i doi | grep -i http | tail -n 1 |  sed -e 's/.*doi.org.//' | awk '{print $NF}' `
# dest doi
getcrossref "$nm" "$doi" "$fn" "$uin"
#endif
out.exit(nm);
return rcx; 
}

//1 html 243 body 340 div 747 main 748 div 754 div 758 div 785 div 789 div 858 div 862 table 865 tbody 1003 tr 1014 td 1017 a 1026 text = https://doi.org/10.1103/PhysRevLett.98.253005
 static  IdxTy doiorg(const InpTy & in , OutTy & out , const StrTy & nm,const StrTy doiurl,const IdxTy xflags=0)  { 
const char * p =doiurl.c_str();
const IdxTy pos=mjm_strings::indexOf(p,"10.");
if (pos==(~0))
{
MM_ERR(MMPR4(nm,__FUNCTION__, pos,doiurl))
return 0; 
}
const StrTy doi=StrTy(p+pos);
MM_ERR(MMPR4(nm,__FUNCTION__, pos,doi))
IdxTy rcx=get_xref(out,in,doi,nm,0);
return 0;
} // doiorg

 static  IdxTy guessasm(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
// apparently uses highwire or the gca thing  citmgr 
const StrTy nm="guessasm"; // +out.partial();
out.enter(nm);
//const StrTy & uin=in.uin();
StrTy url="";
//out.push(nm);
IdxTy rchw=handlehighwire(in,out);
IdxTy rcc;
if (!out.found())
{
rcc=handlecitmgr(in,out);
}
//out.pop();

out.exit(nm);
return 0; 

} // guessasm


 static  IdxTy guessnature(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
#if 0 
TY  - JOUR
AU  - Bongiovanni, Dario
AU  - Klug, Melissa
AU  - Lazareva, Olga
SN  - 2041-4889
UR  - https://doi.org/10.1038/s41419-020-03333-9
DO  - 10.1038/s41419-020-03333-9
ID  - Bongiovanni2021
ER  - 
url=`echo $uin | sed -e 's/[?#&].*//g' | sed -e 's;/abs/;/ris/;g'| sed -e 's;/full/;/ris/;g' |sed -e 's/html$/ris/'| sed -e 's/pdf$/ris/' `
xris=`echo $url | grep ris$`
if [ "$xris" == "" ]
then
#  
url="$url.ris"
#echo appending ris to nature url $url
fi
#echo getting $url
getriscite "$url" "$uin" "$fn" "$fndest"


#endif
const StrTy nm="guessnature";
out.enter(nm);
const IdxTy nstart=out.found();
const StrTy & uin=in.uin();
StrTy url="";
//failing ,  https://www.nature.com/articles/1600670.pdf?origin=ppub
// but conversion to this should work  https://www.nature.com/articles/1600670
// although ading ris does not work 
/*
From hanldemutate,
Mutations v;
v.push_back(" sed -e 's/[.\\/]pdf$//'");
v.push_back(" sed -e 's/\\/pdf$/\\/abs/'");
v.push_back(" sed -e 's/\\/pdf$/\\/full/'");
v.push_back(" sed -e 's;/full/;/abs/;'");
//Mutate(in,out," sed -e 's/[.\\/]pdf$//'");
//if (out.found()==0) { Mutate(in,out," sed -e 's/\\/pdf$/\\/abs/'"); }
//if (out.found()==0) { Mutate(in,out," sed -e 's/\\/pdf$/\\/full/'"); }
MM_LOOP(ii,v)
{
Mutate(in,out,(*ii));
const IdxTy nnow=out.found();
 

*/


//StrTy cmdin="sed -e 's/[?#&].*//g' | sed -e 's;/abs/;/ris/;g'| sed -e 's;/full/;/ris/;g' |sed -e 's/html$/ris/'| sed -e 's/pdf$/ris/'";
// epdf?
// https://www.nature.com/articles/s41598-021-92388-5.epdf?sharing_token=q7U6qoweacU-IeJZEjn-cdRgN0jAjWel9jnR3ZoTv0OyQvQ64BMbUR_MieBX_t2DxyK1tKdB6VoyjL4cSf-WenIByRAaYPixPyQgNE_1YroE4mSqV5PFpbFLruBo1rtFFibqK6Qyz-omWbJls5_OIpG6ghRwLoMkUv-XRiJJKBY%3D
//StrTy cmdin="sed -e 's/[?#&].*//g' | sed -e 's;/abs/;/ris/;g'| sed -e 's;/full/;/ris/;g' |sed -e 's/html$/ris/'| sed -e 's/pdf$/ris/'| sed -e 's/[.]epdf[?].*//' ";
// it has to do epdf first or else it leaves e when doing pdf doh 
StrTy cmdin="sed -e 's/[?#&].*//g' | sed -e 's;/abs/;/ris/;g'| sed -e 's;/full/;/ris/;g' |sed -e 's/html$/ris/'|sed -e 's/epdf$/ris/' | sed -e 's/pdf$/ris/'| sed -e 's/[.]epdf[?].*//' ";
{ 
//const StrTy cmd=cmdin; 
//static IdxTy cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, const IdxTy flags)
cmd_exec(url,uin,cmdin,out,0);
//Blob d,err,cout;
//d=uin;
//IdxTy rc=out.hand().fileio(cout,err,d,cmd,3);
//MM_ERR(MMPR4(cmd,StrTy(cout),StrTy(err),StrTy(d)))
//url=StrTy(cout);

if (url.length()<4) url=uin+".ris";
const IdxTy urllen=url.length();
const char * urlp=url.c_str()+urllen-4;
if (strcmp(".ris",urlp)!=0)  url=url+".ris";
//MM_ERR(MMPR4(StrTy(d),StrTy(err),StrTy(cout),rc)<<MMPR(url))
}

 IdxTy rcfar= FetchAsRis( in , out,  url,  nm, 0) ;
{const IdxTy nnow=out.found();
if (nnow!=nstart) if (!in.collect_all()) { out.exit(nm); return 0; } 
}
/////////////////////////////////

// the ris no longer works, this is the best it can do although json2
// seems to be a good guess.. 
const StrTy planb=" sed -e 's/\\([0-9]\\)\\..*/\\1/g'";
Mutate(in,out,planb);
{const IdxTy nnow=out.found();
if (nnow!=nstart) if (!in.collect_all()) { out.exit(nm); return 0; } 
}
IdxTy rcadhoc=handleadhochtml(in,out,xflags); // in.mom()->Guess(in2,out2,in.rflags());

out.exit(nm);
return 0;


 } // guessnature


//  static  IdxTy handledoilink(const InpTy & in , OutTy & out , const IdxTy xflags=0)  {


 static  IdxTy guesscdnscience(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guesscdnscience";
out.enter(nm);
const IdxTy nstart=out.found();
const StrTy & uin=in.uin();
//StrTy url="";

IdxTy rcadhoc=handledoilink(in,out,xflags); //

out.exit(nm);
return 0;


 } // guessnature




///////////////////////////////////////////////////////

 static  IdxTy guessjes(const InpTy & in , OutTy & out , const IdxTy xflags=0 )  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessacs(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guessacs";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
// sites uses post but accepts get with query 
//wget -O xxx -S -v "https://pubs.acs.org/action/downloadCitation?doi=10.1021%2Facssynbio.0c00022&downloadFileName=achs_asbcd69_1693&direct=true&format=bibtex&include=cit&submit=Download+Citation%28s%29"

StrTy urla="http://pubs.acs.org/action/downloadCitation?direct=true&doi=";
//+doi+
StrTy urlb="&downloadFileName=irrelevant_placeholder&format=bibtex&include=cit&submit=Download+Citation%28s%29";
const StrTy doi=MutateOnly(uin, " sed -e 's/.*doi\\///' | sed -e 's/.*pdf\\///' |  sed -e 's/.*abs\\///'|sed -e 's/.*full\\///' | sed -e 's/[?#&].*//'  | urlencode " ,out );
StrTy url=urla+doi+urlb;
//MM_MSG(url);
Grc grc1=in.getter().normalget(fntemp,url,0);
// const StrTy guessconcat="awk 'BEGIN{c=0}{ if (c==0) { if (($3!=\"{\")||( $4 != \"\"))  {print $0;} else {c=1; line=$0; }} else {if ( substr($0,0,1)==\"}\") {print line$0; line=\"\";c=0;} else {line=line $0}}}END { if (line!=\"\") print line; } '";

 const StrTy & guessconcat=  GuessConcat(in , out )  ;

StrTy cmd= " grep \"[^ \\t]\" |"+guessconcat+" | sed -e 's/  */ /g' ";
IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
inserturl(fnbib,uin,out,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
//const IdxTy nnow=out.found();
#if 0
guessacs()
{
uin="$1"
doi=`echo $uin | sed -e 's/.*abs\///'|sed -e 's/.*full\///' | sed -e 's/[?#&].*//'`
url="http://pubs.acs.org/action/downloadCitation?direct=true&doi=$doi&downloadFileName=irrelevant_placeholder&format=bibtex&include=abs"
guessget "$url" "$uin" | grep "[^ \t]" |guessconcat  | sed -e 's/  */ /g' | tee "$temp1"
citepaste "$temp1" "$uin" "$url"
guessconcat ()
{
awk 'BEGIN{c=0}{ if (c==0) { if (($3!="{")||( $4 != ""))  {print $0;} else {c=1; line=$0; }} else {if ( substr($0,0,1)=="}") {print line$0; line="";c=0;} else {line=line $0}}}END { if (line!="") print line; } '
}
#endif

out.exit(nm);
return 0; 

}
///////////////////////////////////

// https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8036
// https://bonndoc.ulb.uni-bonn.de/xmlui/bitstream/handle/20.500.11811/8036/5508.pdf;jsessionid=8EE2D4C371DD10838C62271A1D5CEC0E?sequence=1
//AddPair("guessbonndoc",&Myt::guessbonndoc,"/bonndoc\..*uni-bonn\.de/");

 static  IdxTy guessbonndoc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessbonndoc";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guessbonndoc



// https://iris.univpm.it/retrieve/handle/11566/263460/77778/Tesi_Cirilli.pdf
// https://iris.univpm.it/retrieve/handle/11566/263460/
// http://hdl.handle.net/11566/263460 
//AddPair("guesshandle",&Myt::guesshandle,"/handle/");
 static  IdxTy guesshandle(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guesshandle";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
StrTy nurl=MutateOnly(uin, "sed -e 's;.*handle\\(/[0-9][0-9]*/[0-9][0-9]*\\).*;\\1;'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
nurl="http://hdl.handle.net"+nurl;
MM_ERR(MMPR3(nm,uin,nurl))
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
IdxTy rcr=Recurse(in,out,nurl);
//IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guesshandle

// change handle to bistream and drop /1/....
// https://sapientia.ualg.pt/bitstream/10400.1/14504/1/Vincent%20Laiz%C3%A9%20-%20FERNANDEZ_2019_BBAGEN_Author%20version.pdf


 static  IdxTy guessbitstream(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessbitstream";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
StrTy nurl=MutateOnly(uin, "sed -e 's;\\(.*\\)/bitstream\\(/[0-9.][0-9.]*/[0-9.].[0-9]*\\).*;\\1/handle\\2;'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
//nurl="http://hdl.handle.net"+nurl;
MM_ERR(MMPR3(nm,uin,nurl))
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
IdxTy rcr=Recurse(in,out,nurl);
//IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guessbitstream


// 2011  #./chromate.out "mode 0" launch "wget xxx.pdf https://authors.library.caltech.edu/32187/1/BROsiamjssc85.pdf"
// 2012  # https://authors.library.caltech.edu/32187/

//AddPair("guesscaltech",&Myt::guesscaltech,"/authors\.library\.caltech"); 
// just remove filename...  should be a mutate gess
 static  IdxTy guesscaltech(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guesscaltech";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
StrTy nurl=MutateOnly(uin, "sed -e 's;\\(caltech.edu/[0-9][0-9]*/\\).*;\\1;'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
//nurl="http://hdl.handle.net"+nurl;
MM_ERR(MMPR3(nm,uin,nurl))
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
IdxTy rcr=Recurse(in,out,nurl);
//IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guesscaltech

//https://apps.dtic.mil/sti/pdfs/ADA121307.pdf
//https://apps.dtic.mil/sti/citations/ADA121307
// 2023-02-25
//AddPair("guessdticmil",&Myt::guesscticmil,"\\.dtic\\.mil/");

 static  IdxTy guessdticmil(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessdticmil";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
StrTy nurl=MutateOnly(uin, "sed -e 's;/pdfs/\\(.*\\).pdf$;/citations/\\1;'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
//nurl="http://hdl.handle.net"+nurl;
MM_ERR(MMPR3(nm,uin,nurl))
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
// this really needs to combine a few as the htlml ius fcked up
IdxTy rcr=Recurse(in,out,nurl);
//IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guessdticmil






       



//////////////////////////////////////

// http://su.diva-portal.org/smash/record.jsf?dswid=-8098&pid=diva2%3A189901&c=1&searchType=SIMPLE&language=en&query=189901&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all
// http://su.diva-portal.org/smash/get/diva2:189901/FULLTEXT01

 static  IdxTy guessdiva(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessdiva";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const StrTy nurl=MutateOnly(uin, "sed -e 's/xmlui\\/bitstream/xmlui/'| sed -e 's/\\/[0-9]*\\.pdf.*//'" ,out );
MM_ERR(MMPR3(nm,uin,nurl))
if ((nurl!=uin)&&( nurl.length()!=0))
{
InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rc=handleisbib(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);
}
// otherwise isbib should still work later... but check for others 
out.exit(nm);
return 0;

} // guessbonndoc



///////////////////////////////////////

#if 0 
guessrsc()
{
#http://pubs.rsc.org/en/content/articlelanding/2004/cp/b405084c#!divAbstract
base="http://pubs.rsc.org/en/content/getformatedresult/"
uin="$1"
fn="$temp1"
nm=`echo "$uin" | sed -e 's/#.*//' | sed -e 's/.*\///' `
url="$base$nm?downloadtype=article"
#$WGET --post-data="ResultAbstractFormat=BibTex" -O xxxx -S -v "http://pubs.rsc.org/en/content/getformatedresult/b405084c?downloadtype=article"
#$WGET --post-data="ResultAbstractFormat=BibTex" -O "$fn"  -S -v "$url"
$WGET -q --post-data="ResultAbstractFormat=BibTex" -O "$fn"  "$url"
cat "$fn"
citepaste "$fn" "$uin" "$url"

}

Note that this is also a doi, i
https://pubs.rsc.org/en/content/articlelanding/1994/FT/FT9949001581
Article information
DOI	https://doi.org/10.1039/FT9949001581

#endif
// now this works as a GET,
// https://pubs.rsc.org/en/content/formatedresult?markedids=ft9949001581&downloadtype=article&managertype=bibtex"

// $WGET -q --post-data="ResultAbstractFormat=BibTex" -O "$fn"  "$url"
// /tmp/fileIwPpyI.outguesses
 // this returns output in wrong format, ignoring post data?
// https://pubs.rsc.org/en/content/getformatedresult/ft9949001581?downloadtype=article
// ResultAbstractFormat=BibTex&go=Go
// https://pubs.rsc.org/en/content/formatedresult?markedids=ft9949001581&downloadtype=article&managertype=bibtex
// markedids=ft9949001581&downloadtype=article&managertype=bibtex

/*
--post-data="markedids=FT9949001581&downloadtype=article&managertype=bibtex"  -O "/tmp/file55EpdH.outguesses" "http://pubs.rsc.org/en/content/getformatedresult?markedids=FT9949001581&downloadtype=article&managertype=bibtex" < /tmp/fileL2mPZn.pawnoff >/tmp/fileRcdfM4.pawnoff 2>/tmp/fileL4eFyL.pawnoff ; echo $? >> /tmp/fileX3b5ks.pawnoff


*/
//const StrTy pd="markedids="+_uin+"&downloadtype=article&managertype=bibtex";


// after all of that the publisher's result is not as good as the doi...

 static  IdxTy guessrsc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guessrsc";
out.enter(nm);
const bool all=in.collect_all();
const IdxTy nstart=out.found();
const StrTy & uin=in.uin();
// old way 
const StrTy basep="http://pubs.rsc.org/en/content/getformatedresult/";
const StrTy baseg="https://pubs.rsc.org/en/content/formatedresult?markedids=";
// should do to lower... 
const StrTy _uin=MutateOnly(uin,"sed -e 's/#.*//' |sed -e 's/.*\\///'",out );
MM_ERR(MMPR3(nm,uin,_uin))
if ( _uin=="" ) { out.exit(nm); return 0; } 
const StrTy urlpost=basep+_uin+"?downloadtype=article";
const StrTy urlget=baseg+_uin+"&downloadtype=article&managertype=bibtex";

const StrTy fnbib=out.fn("bibtex"); //
Grc grc1=in.getter().normalget(fnbib,urlget,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,urlget,nm+" (get) ",0);


IdxTy nfound=out.found();
if ( (nstart!=nfound) && ! all )
{ out.exit(nm); return 0; } // found

const StrTy fnbibp=out.fn("bibtex"); //
const StrTy pd="ResultAbstractFormat=BibTex&go=Go";
Grc grc2=in.getter().normalpost(fnbibp,urlpost,pd,0);
IdxTy rcfp= out.good_enough(b,fnbibp,in,urlpost,nm+" (post) ",0);

nfound=out.found();
if ( (nstart!=nfound) && ! all )
{ out.exit(nm); return 0; } // found

const StrTy doi="10.1039/"+_uin;
MM_ERRF(MMPR4(nm,doi,nfound,nstart))
IdxTy rcx=get_xref(out,in,doi,nm+" (doi) " ,0);
out.exit(nm);
return 0; 

} // guessrsc

 static  IdxTy guessscidirect(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
//pii=`echo $uin | sed -e 's/.*\///g' | sed -e 's/[^0-9A-Za-z]//g'`
//debugstrat  guessscidirect uin  $uin pii  $pii
//url201803="https://www.sciencedirect.com/sdfe/arp/cite?pii=$pii&format=text/x-bibtex&withabstract=false"
const StrTy nm="guessscidirect";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
const StrTy urli="https://www.sciencedirect.com/sdfe/arp/cite?pii=";
const StrTy urlf="&format=text/x-bibtex&withabstract=false";
//const StrTy cmdin="( sed -e 's/.*\\///g' | sed -e 's/[^0-9A-Za-z]//g' ) ";
// now insert parens at the call 
const StrTy cmdin=" sed -e 's/.*\\///g' | sed -e 's/[^0-9A-Za-z]//g'  ";
cmd_exec(url,uin,cmdin,out,0);
url=urli+url+urlf;
const StrTy fnbib=out.fn("bibtex"); //
Grc grc=in.getter().normalget(fnbib,url,16);
Blob b;
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0; 
}// guessscidirect
// eutilsnew -out xxx -db pubmed -report medline 12633389
// DOPMC="$EU -v -out $eu_output -db pmc -report medline"

 static  IdxTy guesspubmed(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guesspubmed";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy cmd="sed -e 's/\\/$//'  | sed -e 's/.*\\///g' | sed -e 's/[?#].*//'";
const StrTy db="pubmed";
StrTy pmcid;
cmd_exec(pmcid,uin,cmd,out,0);
if ( pmcid.length()==0) {  out.exit(nm); return 0;   }
const bool get_anyway=!in.punt_pmc(); // !false;
if (get_anyway)
{
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const StrTy fnpmc=out.fn("pmc"); //
//AccessTy pubmed(const StrTy & fn, const StrTy& id, const StrTy& db, const IdxTy flags)
Grc grc=in.getter().pubmed(fnpmc,pmcid,db,0);
BibEntry be;
Xrc xrc=in.xformer().pmc_to_bib(be,fnpmc,0);
Blob b;
b=be.format();
b.save(fnbib);
//StrTy url=uin;
StrTy url=grc.url();
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0;
} // get_anyway
// we may have tried but punting is still ok if failure... 
out.punt(); 
out.cite(StrTy("PMID"+pmcid)); 
MM_ERR("  punting "<<MMPR2(nm,pmcid))  


#if 0 
# there may be a "#" or "?" etc  something stupid at the end 
num=`echo "$uin" | sed -e 's/\/$//'  | sed -e 's/.*\///g' | sed -e 's/[?#].*//'`
#inum=`echo $num | grep -v "[0-9 ]" `
##if [ "$inum" != "" ]
if [ "$num" == "" ]
then
echo wrong format try alt  | smsg
num=`echo "$uin" | sed -e 's/[^0-9]//g'`
fi
citex="\\cite{PMID$num}"

#endif

//MM_ERR(" not impl")  
out.exit(nm);
return 0; 

}


// https://www.ijvsbt.org/index.php/journal/article/view/1386
// https://www.ijvsbt.org/index.php/journal/article/download/1386/1058
//AddPair("guessijvsbt",&Myt::guessijvsbt,"\\.ijvsbt\\.org/");

//IdxTy rcr=Recurse(in,out,newlink);


 static  IdxTy guessijvsbt(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessijvsbt";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy _fn=out.fn(); //
// at this point mutating the pdf url would be good but then
// the file name can't pollute the cache. 
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fn=mut?fntemp:_fn;
const StrTy _uin=MutateOnly(uin, "(x=`cat -`; echo $x| sed -e 's;download/\\([0-9]*\\)/\\([0-9]*\\);view/\\1;' ) " ,out );
MM_ERR(MMPR3(nm,uin,_uin))
if ( _uin=="" ) { out.exit(nm); return 0; } 
if (_uin==uin)
{
IdxTy rc=handleadhochtml(in,out,xflags);
}
else
{
IdxTy rcr=Recurse(in,out,_uin);
}
out.exit(nm);
return 0; 

}
////////////////////////////////////////////////////////

 static  IdxTy guessmdpi(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesssmdpi";
out.enter(nm);
const bool all=in.collect_all();
const IdxTy nstart=out.found();
const StrTy & uin=in.uin();
const StrTy _fn=out.fn(); //
// at this point mutating the pdf url would be good but then
// the file name can't pollute the cache. 
const StrTy fntemp=out.fn("temp"); //
const bool mut=(strncmp(uin.c_str()+uin.length()-3,"pdf",3)==0);
const StrTy fn=mut?fntemp:_fn;
const StrTy _uin=mut?MutateOnly(uin, " sed -e 's/[.\\/]pdf$//'" ,out ):uin;
Grc grc=in.getter().normalget(fn,_uin,16);
StrTy uinclean,url,aid;
StrTy urlb="http://www.mdpi.com/export?articles_ids%5B%5D=";
StrTy urlc="&export_format_top=bibtex&export_submit_top=";
{
const StrTy cmdin="sed -e 's/\\/[a-z][a-z]*$//'";
cmd_exec(uinclean,uin,cmdin,out,0);
}
Grc grc1=in.getter().normalget(fntemp,uinclean,0);
{
const StrTy cmdin=" grep articles_ids | sed -e 's/[^0-9]//g' | sort | uniq";
//cmd_exec(aid,fntemp,cmdin,out,2);
cmd_exec(aid,fntemp,cmdin,out,2+4); // 4 is no trailing cr lf doh 
}
url=urlb+aid+urlc;
const StrTy pd=StrTy("articles_ids[]=")
	+aid+"&export_format_top=bibtex&export_submit_top=";

//const StrTy fnbib=out.fn(nm); //
const StrTy fnbib=out.fn("bibtex"); //
url="http://www.mdpi.com/export";
Grc grc2=in.getter().normalpost(fnbib,url,pd,0);
//fixcitename "$fn" "$uin" "$url"| inserturl "$uin"  > "$temp5"
BibEntryMap bem;
bem.parse(fnbib);
if (bem.size()!=1)
{
MM_ERR(MMPR(nm)<< " fails "<<MMPR(bem.size())<<MMPR2(fn,fnbib))
out.exit(nm);
return 2;
}
BibEntry  be=*(bem.begin()); // ((*(bem.begin())).second)[0];
// the mdpi result is not as good as crossref now.
// if it has a doi go get xref first... 
// avoiding canonicl here but may be a good idea.. 
//////////// similar to issue with rg
// TODO in theory there could be an entry for doi with zero
// length but not sure it is possible... 
{ // scoping
StrTy doi=be.get("DOI",0); // be.uin();
// should just canonicalize doh 
if (doi.length()==0) doi=be.get("doi",0); // be.uin();
if (doi!="")
{
const StrTy fntemp2=out.fn("temp"); //
Grc grc=in.getter().crossref(fntemp2,doi,0);
const bool doiok=(grc.rc()=="200");
//added_doi=doiok;
//const StrTy fnbest=doiok?fntemp2:fnbib;
const StrTy fnbest=fntemp2;
MM_ERR(MMPR3(doiok,grc.rc(),fnbest))
Blob b;
b.load(fnbest);
IdxTy rcf=doiok?  out.good_enough(b,fnbest,in,grc.url(),nm+" (xref) ",0)
: 0 ; // need to pick one place or the other fix logi here, want doi first if available. 
// : out.good_enough(b,fnbest,in,grc2.url(),nm+" (rg) ",0);
const IdxTy nfound=out.found();
if ( (nstart!=nfound) && ! all )
{ out.exit(nm); return 0; } // found
} // doi
}// scoping  


//////////////////////////////////////
be.make_name();
Blob b;
b=be.format();
b.save(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url+" "+pd,nm,0);
#if 0
uinclean=`echo $uin | sed -e 's/\/[a-z][a-z]*$//'`
aid=`lynx -dump -source "$uinclean" | grep articles_ids | sed -e 's/[^0-9]//g' | sort | uniq`
url="http://www.mdpi.com/export?articles_ids%5B%5D=$aid&export_format_top=bibtex&export_submit_top="
pd="articles_ids[]=$aid&export_format_top=bibtex&export_submit_top="
urlb="http://www.mdpi.com/export"
# this now requires a post not  get 
# wget -O xxxx -UMozilla -S -v "http://www.mdpi.com/export" --post-data="articles_ids[]=222355&export_format_top=bibtex&export_submit_top="
#normalget "$fn" "$url"
normalpost "$fn" "$urlb" "$pd"
#fixcitename "$fn" "$uin" "$url" > "$temp5"
fixcitename "$fn" "$uin" "$url"| inserturl "$uin"  > "$temp5"

#endif


out.exit(nm);
return 0; 

}


 static  IdxTy guessspiedigital(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessoxford(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessoxfordpcp(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesstandf(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessspringer(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
// had a problem with a book appendix,
// \url{https://link.springer.com/content/pdf/bbm%3A978-1-4615-2131-0%2F1.pdf }
// this link worked, but because it could find the doi, 
// https://link.springer.com/book/10.1007/978-1-4615-2131-0
// this returned something but not bibtex, 
// https://citations.springer.com/book?doi=10.1007/978-1-4615-2131-0&format=bibtex
const StrTy nm="guesssspringer";
out.enter(nm);
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //

const StrTy base="http://citation-needed.services.springer.com/v2/references/";

const StrTy sfx="?format=bibtex&flavour=citation";
// name=`echo $uin | sed -e 's;.*com/[a-z][a-z]*/;;' | sed -e 's/[?#].*//' | sed -e 's;%2F;/;g' `

StrTy name; 

{const StrTy cmdin="sed -e 's;.*com/[a-z][a-z]*/;;' | sed -e 's/[?#].*//' | sed -e 's;%2F.*;/;gi' | sed -e 's/^pdf.//' ";
cmd_exec(name,uin,cmdin,out,0);
MM_ERR(" url name .  "<<MMPR(name))
}
if (name.length()!=0) //  { out.exit(nm); return 0; }
{
StrTy url=base+name+sfx;
Grc grc1=in.getter().normalget(fntemp,url,0);
//cat "$fn"  |  grep "[^ \t]" |guessconcat  | sed -e 's/  */ /g' >  "$fndest"
//StrTy bib,b2;
Blob b;
//b.load(fntemp);
//bib=StrTy(b);
//b.save(fnbib);
//static IdxTy CleanUpBadBib( const StrTy & fnbib, const IdxTy flags) 
IdxTy cubbrc=CleanUpBadBib(b, fntemp,fnbib, 0) ;
//b.load(fnbib);
if (cubbrc!=0)
{
MM_ERR(" failed to convert bib "<<MMPR3(fnbib,fntemp,cubbrc))
//out.exit(nm);
//return cubbrc; 
}
else
{
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
}
}// name

const IdxTy nfound=out.found();
if ( (nstart!=nfound) &&! all ) { out.exit(nm); return 0; } 

// https://api.crossref.org/works/10.1007/978-1-4615-2131-0/1.pdf
// \url{https://link.springer.com/content/pdf/bbm%3A978-1-4615-2131-0%2F1.pdf }
// https://link.springer.com/book/10.1007/978-1-4615-2131-0
// assume always url encoded... 
// ./mjm_med2bib_funcs.h896  url name .   name=pdf/bbm%3A978-1-4615-2131-0/1.pdf
{const StrTy cmdin="sed -e 's/.*%3A//'  | sed -e 's/[?#].*//' | sed -e 's;%2F.*;/;g' ";
cmd_exec(name,uin,cmdin,out,0);
MM_ERR(" url name should be doi  .  "<<MMPR(name))
}
if (name.length()) {  // out.exit(nm); return 0; }
StrTy doi="10.1007/"+name;
if (doi!="")
{
const StrTy fntemp2=out.fn("temp"); //
Grc grc=in.getter().crossref(fntemp2,doi,0);
const bool doiok=true; // (grc.rc()=="200");
//added_doi=doiok;
const StrTy fnbest=fntemp2; // doiok?fntemp2:fnbib;
MM_ERR(MMPR3(doiok,grc.rc(),fnbest))
Blob b;
b.load(fnbest);
IdxTy rcf=doiok?  out.good_enough(b,fnbest,in,grc.url(),nm+" (xref) ",0)
: 0 ; // need to pick one place or the other fix logi here, want doi first if available. 
// : out.good_enough(b,fnbest,in,grc2.url(),nm+" (rg) ",0);
MM_ERR(MMPR4(nm,grc.rc(),out.found(),fnbest))
} // doi

 } // name 
out.exit(nm);
return 0; 
} // guessspringer 



// this can't possibly work 
// guessconcat ()
//{
//awk 'BEGIN{c=0}{ if (c==0) { if (($3!="{")||( $4 != ""))  {print $0;} else {c=1; line=$0; }} else {if ( substr($0,0,1)=="}") {print line$0; line="";c=0;} else {line=line $0}}}END { if (line!="") print line; } '
//}

#if 0 
{

//const StrTy awk="awk 'BEGIN{c=0}{ if (c==0) { if (($3!=\"{\")||( $4 != \"\"))  {print $0;} else {c=1; line=$0; }} else {if ( substr($0,0,1)==\"}\") {print line$0; line=\"\";c=0;} else {line=line $0}}}END { if (line!=\"\") print line; } '";
//const StrTy cmdin="grep \"[^ \t]\" |guessconcat  | sed -e 's/  */ /g'";
//const StrTy cmdin="grep \"[^ \\t]\" |"+awk+"| sed -e 's/  */ /g'";
const StrTy cmdin="grep \"[^ \\t]\" | sed -e 's/  */ /g'";
cmd_exec(b2,bib,cmdin,out,0);
MM_ERR(" g2 .  "<<MMPR(b2))
}
const StrTy fnbib=out.fn("bibtex"); //
b=b2;
b.save(fnbib);
Bibmap bm;
bm.parse(fnbib);
if (bm.size()!=0)
{


}
#endif




#if 0

guessspringer()
nm="guessspringer"
base="http://citation-needed.services.springer.com/v2/references/"
sfx="?format=bibtex&flavour=citation"
ty=`echo $uin | grep chapter`
name=`echo $uin | sed -e 's;.*com/[a-z][a-z]*/;;' | sed -e 's/[?#].*//' | sed -e 's;%2F;/;g' `
url="$base$name$sfx"
# http://link.springer.com/article/10.1007%2FBF02369136?LI=true
fn="$temp1"
fndest="$temp5"
normalget "$fn" "$url"
cat "$fn"  |  grep "[^ \t]" |guessconcat  | sed -e 's/  */ /g' >  "$fndest"
citepaste "$fndest" "$uin" "$url"




#endif


 static  IdxTy guessasce(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesswiley(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guesswiley";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
const StrTy name=
	MutateOnly(uin, "sed -e 's/.*\\/10\\./10./'| sed -e 's/%2[fF]/\\//g'" ,out);
if (name.length()==0) { out.exit(nm);  return 0; } 
const StrTy base="https://onlinelibrary.wiley.com/action/downloadCitation?doi=";
const StrTy sfx="&citation-type=bibtex&format=bibtex";
const StrTy url=base+name+sfx;
Grc grc1=in.getter().normalget(fntemp,url,0);
//IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
//Blob b; b.load(fnbib);
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
 const StrTy & guessconcat=  GuessConcat(in , out )  ;
//cat "$temp5"  |  grep "[^ \t]" |guessconcat  | sed -e 's/  */ /g' | sed -e 's/,,/,/g'  > "$fn"
IdxTy rc=filter_small_file(fnbib, fntemp,
"grep \"[^ \\t]\" |"+guessconcat+" | sed -e 's/  */ /g' | sed -e 's/,,/,/g'" 
 , out, 0);

Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

out.exit(nm);


#if 0

guesswiley()
uin="$1"
base="http://onlinelibrary.wiley.com/enhanced/getCitation/doi"
sfx="?citation-type=bibtex"
base="https://onlinelibrary.wiley.com/action/downloadCitation?doi="
sfx="&citation-type=bibtex&format=bibtex"

name=`echo "$uin" |  sed -e 's/.*\/10\./10./'| sed -e 's/%2[fF]/\//g'`
url="$base$name$sfx"
fn="$temp1"
normalget "$temp5"  "$url"
echo wiley double comma fick  | smsg
cat "$temp5"  |  grep "[^ \t]" |guessconcat  | sed -e 's/  */ /g' | sed -e 's/,,/,/g'  > "$fn"

#endif



return 0; 

} // guesswiley 


 static  IdxTy guesssci(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesssci";
out.enter(nm);
//const StrTy & uin=in.uin();
//const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// two approaches- get headless and find doi OR use highwire on .citation-tools
//const StrTy nurl=MutateOnly(uin, "sed -e 's/\\([0-9]\\)[.?#;%].*/\\1/g' | awk '{print $0\".citation-tools\"}'" ,out );
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
IdxTy rc=handlehighwire(in,out,xflags); // in.mom()->Guess(in2,out2,in.rflags());
//MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
IdxTy nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rc; }
 out.exit(nm); return rc; 

} // guesssci


 static  IdxTy guessrupress(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesscopyclip(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessfuturesci(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesspnas(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

MM_ERR(" obsolete, highwire and doi work now implement for speed  ")  
//return 0; 
const StrTy nm="guesspnas";
out.enter(nm);
const StrTy & uin=in.uin();
IdxTy rc= handledoilink(in,out,xflags); 
out.exit(nm);
return rc; 

} // guesspnas


 static  IdxTy guessscielo(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessjama(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesskarger(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guesskarger";
out.enter(nm);
const StrTy & uin=in.uin();
StrTy url="";
StrTy urlbase="http://www.karger.com/Article/CitationEndNote/";
//numb=`echo $uin | sed -e 's/[?#&].*//g' | sed -e 's/.*\///g'| sed -e 's/[^0-9]//g'  `
StrTy cmdin=" sed -e 's/[?#&].*//g' | sed -e 's/.*\\///g'| sed -e 's/[^0-9]//g' ";
StrTy numb;
cmd_exec(numb,uin,cmdin,out,0);
if (numb.length()==0) 
{ MM_ERR(" blank number "<<MMPR3(nm,uin,cmdin)) }
url=urlbase+numb;
// returning without a name.. 
 IdxTy rcf= FetchAsRis( in , out,  url,  nm, 0) ;
out.exit(nm);
return rcf; 
} // guesskarger

 static  IdxTy guesssage(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesssage"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
const IdxTy nstart=out.found();
const StrTy cmd=" sed -e 's/.*\\(10\\.[0-9][0-9][0-9][0-9a-zA-Z.]*\\/[-0-9a-zA-Z.]*\\).*/\\1/'"; 
const StrTy doi=MutateOnly(uin, cmd,out);
MM_ERR(MMPR3(uin,doi,cmd))
if (doi.length()>4 ) {
const StrTy fnbib=out.fn("bibtex"); //
// use sage for now 
//IdxTy rcx=get_xref(out,in,doi,nm,0);
const StrTy url="http://journals.sagepub.com/action/downloadCitation";
const StrTy pd="doi="+doi+"&include=abs&format=bibtex&direct=true&submit=Download+Citation";
Grc grc2=in.getter().normalpost(fnbib,url,pd,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url+" "+pd,nm,0);
IdxTy nnow=out.found();
if (!in.collect_all()) if (nstart!=nnow){ out.exit(nm); return rcf; }
} // doi 

out.exit(nm);
return 0; 
//MM_ERR(" not impl")  return 0; 
#if 0
#pd="articles=$aid&ArticleAction=export_bibtex"
url="http://journals.sagepub.com/action/downloadCitation"
#doi=`echo $uin| sed -e 's/.*\(10\.117[0-9a-zA-Z.]*\/[0-9a-zA-Z.]*\).*/\1/' `
# AFFCUK 
doi=`echo $uin| sed -e 's/.*\(10\.[0-9][0-9][0-9][0-9a-zA-Z.]*\/[-0-9a-zA-Z.]*\).*/\1/' `
pd="doi=$doi&include=abs&format=bibtex&direct=true&submit=Download+Citation"
echo blank response try post doi=$doi | smsg
fn="$temp1"
fn5="$temp5"
rm -f "$fn"
normalpost "$fn" "$url" "$pd"


#endif


} // guesssage


 static  IdxTy guessingenta(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

#if 0
guessciteseer()
{
nm="guessciteseer"
echo try $nm | smsg
#lynx -dump "http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.462.9736" | sed -n -e '/BibTeX/,/Share/p'
uin="$1"
base="http://citeseerx.ist.psu.edu/viewdoc/summary?doi="
doi=`echo $uin  | sed -e 's/.*doi=//' | sed -e 's/[?#&].*//'`
url="$base$doi"
fn="$temp1"
fnss="$temp5"
rm -f "$fn"
# this may put in CR 
lynx -width 10000 -dump "$url" | sed -n -e '/BibTeX/,/Share/p' | sed -e 's/^   //' | grep -v "^BibTeX$" | grep -v "^Share$" > "$fn"
citepastecat  "$fn" "$uin" "$url" "$nm" > "$fnss"
if [ "$handled" != "" ]
then
cat "$fnss"
else
echo $nm did not work `ls -al "$fnss"` | smsg
fi

}

#endif

 static  IdxTy guessciteseer(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guessciteseer";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fntemp3=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// from cmd line, the profiles , adding 16 to flags, seem to be a problem  
//Grc grc=in.getter().meget(fntemp,uin,1+8);
//MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))

const StrTy base="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=";
const StrTy cmd="sed -e 's/.*doi=//' | sed -e 's/[?#&].*//'";
const StrTy doi=MutateOnly(uin, cmd,out);
MM_ERR(MMPR3(uin,doi,cmd))
if (doi.length()>4 ) {
const StrTy url=base+doi;
// this should use the getter or renderer....
const StrTy bcmd="lynx -width 10000 -dump \""+url+"\" | sed -n -e '/BibTeX/,/Share/p' | sed -e 's/^   //' | grep -v \"^BibTeX$\" | grep -v \"^Share$\"";
/////////////////////////////////
StrTy txt;
IdxTy rcc=cmd_exec(txt,"",bcmd,out,0);

//IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
//const StrTy url=base+sfx+ssfx;
// should do dos2unix and concat... 
//Grc grc2=in.getter().normalget(fntemp,url,0);
Blob b;
b=txt;
b.save(fnbib);
//b.load(fnbib);
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
//StrTy cmd2= " dos2unix |"+guessconcat;
//IdxTy rc=filter_small_file(fnbib, fntemp, cmd2, out, 0);
inserturl(fnbib,uin,out,0);
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

////////////////////////////////
const IdxTy nnow=out.found();
if (!all) if (nstart!=nnow) { out.exit(nm); return 0; } 
// this can not work, wrong doi doh
IdxTy rcx=get_xref(out,in,doi,nm,0);


}// doi

out.exit(nm);
return 0;



} // guessciteseer



 static  IdxTy guessosa(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessscientificnet(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessiop(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
// https://iopscience.iop.org/article/10.1088/0034-4885/79/8/084502
const StrTy nm="guessiop"; // +out.partial();;
out.enter(nm);
const StrTy & uin=in.uin();
const IdxTy nstart=out.found();
if (true) 
{
const StrTy cmd=("sed -e 's/.*article\\/10/10/' | grep \"^10\" | grep \"/\" "); // doi=, plos
const StrTy doi=MutateOnly(uin, cmd,out);
MM_ERR(MMPR3(uin,doi,cmd))
if (doi.length()>4 ) {
IdxTy rcx=get_xref(out,in,doi,nm,0);
IdxTy nnow=out.found();
if (!in.collect_all()) if (nstart!=nnow){ out.exit(nm); return rcx; }  
}
}

const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //

//Grc grc=in.getter().normalget(fn,uin,17+128);
Grc grc=in.getter().headlessget(fntemp,uin,1);
//StrTy isolatedomain=" sed -e 's;\\([^/]\\)/[^/].*;\\1;'"; // isolatedomain
const StrTy base=MutateOnly(uin,isolatedomain(),out);
Blob b;
//Xrc xrc=in.xformer().html_to_rendered(b,fn,0);
b.load(fntemp);
const StrTy isolatehref="sed -e 's/.*href=[\"]//' | sed -e 's/[\"].*//'";
const StrTy grep_biblink=MutateOnly(StrTy(b), " grep -i export_bib ",out );
MM_ERR(MMPR(grep_biblink))
const StrTy link=MutateOnly(StrTy(b), " grep -i export_bib |"+isolatehref,out );
const StrTy url=base+StrTy("/")+link;
MM_ERR(MMPR3(url,base,link))
Grc grc1=in.getter().normalget(fnbib,url,0+128);
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
const IdxTy nnow=out.found();
if (!in.collect_all()) if( nnow==nstart)
{
Xrc xrc= in.xformer().html_to_parsed(b,fn,0);
const StrTy linkp=
	MutateOnly(StrTy(b), " grep -i export_bib |awk '{print $NF}'",out );
const StrTy url=base+StrTy("/")+linkp;
MM_ERR(MMPR3(url,base,linkp))
Grc grc1=in.getter().normalget(fnbib,url,0+128);
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

} // parser 
out.exit(nm);
  return 0; 
#if 0
guessiop()
# the export id does not match the  in the url  
# lynx -dump "http://iopscience.iop.org/export?articleId=0953-8984/15/1/302&exportFormat=iopexport_bib&exportType=abs&navsubmit=Export%2Babstract"
# lynx -dump "http://iopscience.iop.org/article/10.1088/0953-8984/15/1/302/meta" | grep -i export_bib
#  27. http://iopscience.iop.org/export?articleId=0953-8984/15/1/302&exportFormat=iopexport_bib&exportType=abs&navsubmit=Export%2Babstract
uin="$1"
# in theory we can just derive it or fetch it 
# the html would work too but needs differen fromat 
#url=`lynx -dump "$uin" | grep -i export_bib | awk '{print $2}' `
fck=`lynx -dump -source "$uin" | grep -i export_bib | isolatehref `
#base="http://iopscience.iop.org/export?articleId="
base=`echo $uin | isolatedomain`
# the name in the uin does not match the  xport tag
#sfx="&exportFormat=iopexport_bib&exportType=abs&navsubmit=Export%2Babstract"
#name=`echo  "$uin" | sed -e 's/.*iop.org\///g' | sed -e 's;/; ;g' | awk '{print $3"/"$4"/"$5"/"$6}'`
#url="$base$name$sfx"
url="$base/$fck"
citetempone "$uin" "$url"
#endif


} // guessiop
/*
uin="$1"
# 2453  myclip -degoogle
# 2454  lynx -dump "http://www.jci.org/articles/view/101223/pdf/render
# 2455  lynx -dump "http://www.jci.org/articles/view/101223"
# 2456  lynx -dump "http://www.jci.org/articles/view/101223/cite"
# 2457  lynx -dump "http://www.jci.org/articles/view/101223/cite/bibtex"

uintrim=`echo $uin | sed -e 's;/[^0-9][^0-9]*$;;'`
url="$uintrim/cite/bibtex"
citetempone   "$uin" "$url" "guessjci"
}

http://content-assets.jci.org/manuscripts/106000/106658/JCI71106658.pdf
https://www.jci.org/articles/view/106658
https://doi.org/10.1172/JCI106658

*/

 static  IdxTy guessjci(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
//MM_ERR(" not impl")  return 0; 

const StrTy nm="guessjci";
out.enter(nm);
const StrTy & uin=in.uin();
//const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const StrTy cmd="sed -e 's;/;\\n;g' |grep -v \"[^0-9]\"  | tail -n 1";
const StrTy id=MutateOnly(uin, cmd ,out );

//const StrTy base="http://www.jci.org/articles/view/101223/cite/bibtex"
const StrTy base="http://www.jci.org/articles/view/"; //101223/cite/bibtex"
const StrTy sfx="/cite/bibtex";
const StrTy url=base+id+sfx;
MM_ERR(MMPR4(nm,uin,id,url)<<MMPR(cmd))
//const StrTy url=base+name;
Grc grc1=in.getter().normalget(fnbib,url,0);
 
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

out.exit(nm);
  return 0; 


} // guessjci

 static  IdxTy guessgoopat(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }


#if 0
 2186  wget -O xxx.html -S -v "https://ieeexplore.ieee.org/abstract/document/9406427"
 2187  grep -i doi xxx.html
 2188  grep -i doi xxx.html | more
 2189  lynx -dump xxx.html | grep -i doi
 2190  toobib -hhtml xxx.html | grep doi
 2191  toobib -hhtml xxx.html | grep -i  doi
 2192  toobib -hhtml xxx.html | more
 2193  toobib -hhtml xxx.html | more | grep xplGlob
 2194  toobib -hhtml xxx.html | more | grep xplGlobal.document.metadata
 2195  toobib -hhtml xxx.html | more | grep xplGlobal.document.metadata | sed -e 's/.*xplGlobal//'
 2196  toobib -hhtml xxx.html | more | grep xplGlobal.document.metadata | sed -e 's/.*xplGlobal.document.metadata//'
 2197  toobib -hhtml xxx.html | more | grep xplGlobal.document.metadata | sed -e 's/.*xplGlobal.document.metadata//' | more
 2198  toobib -hhtml xxx.html | more | grep xplGlobal.document.metadata | sed -e 's/.*xplGlobal.document.metadata=//' > xxx.json
 2199  toobib -json xxx.json
 2200  viw toobib
 2201  toobib -json1 xxx.json
 2202  toobib -json1 xxx.json | more

# lynx -width 10000 -dump "http://ieeexplore.ieee.org/xpl/downloadCitations?download-format=download-bibtex&citations-format=citation-only&recordIds=293357" 
#http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&arnumber=293357
base="http://ieeexplore.ieee.org/xpl/downloadCitations?download-format=download-bibtex&citations-format=citation-only&recordIds="
# the doi is not included with citation only and extra comma doh 
base="http://ieeexplore.ieee.org/xpl/downloadCitations?download-format=download-bibtex&citations-format=abstract&recordIds="
arn=`echo "$uin" | sed -e 's/.*arnumber=//' | sed -e 's/[^0-9]//g'`
url="$base$arn"
normalget "$fn" "$url"
echo  >> "$fn"
fn2sht="$temp2"
#asfick="afick"
#rm -f "$asfick"
cat "$fn" | sed -e 's/<br>//g' | sed -e 's/,\}/\n\}/' > asfick


#endif

// https://ieeexplore.ieee.org/abstract/document/9406427

 static  IdxTy guessieee(const InpTy & in , OutTy & out , const IdxTy xflags=0)

  { 


const StrTy nm="guessieee";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fntemp=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// from cmd line, the profiles , adding 16 to flags, seem to be a problem  
//Grc grc=in.getter().meget(fntemp,uin,1+8);
//MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))

const StrTy base="http://ieeexplore.ieee.org/xpl/downloadCitations?download-format=download-bibtex&citations-format=abstract&recordIds=";
const StrTy cmd="sed -e 's/.*arnumber=//' | sed -e 's/[^0-9]//g'";
const StrTy arn=MutateOnly(uin, cmd,out);
MM_ERR(MMPR3(uin,arn,cmd))
if (arn.length()>4 ) {
const StrTy url=base+arn;
// this should use the getter or renderer....

//IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
//const StrTy url=base+sfx+ssfx;
// should do dos2unix and concat... 
Grc grc2=in.getter().normalget(fntemp,url,0);
const StrTy cmd=" sed -e 's/<br>//g' | sed -e 's/,\\}/\\n\\}/'";
IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
Blob b;
b.load(fnbib);
//b.load(fnbib);
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
//StrTy cmd2= " dos2unix |"+guessconcat;
//IdxTy rc=filter_small_file(fnbib, fntemp, cmd2, out, 0);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
}//arn 

////////////////////////////////
const IdxTy nnow=out.found();
if (!all) if (nstart!=nnow) { out.exit(nm); return 0; } 
//IdxTy rcx=get_xref(out,in,doi,nm,0);
// try to decode the json... 


out.exit(nm);
return 0;

} // guessieee




 static  IdxTy guessbegel(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessaps(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 


const StrTy nm="guessaps";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const StrTy base=MutateOnly(uin, isolatedomain() ,out );

const StrTy ssfx="?type=bibtex&download=true";
//const StrTy url=base+name;
Grc grc1=in.getter().normalget(fn,uin,16);
StrTy sfx;
 
//static IdxTy cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, OutTy & out, const IdxTy flags)
const StrTy cmd=" sed -e 's/\"/\\n/g' | grep xport  | grep -v \"[<>]\" | grep \"/\"| head -n 1";
IdxTy rcc=cmd_exec(sfx,fn,cmd,out,6);

//IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
const StrTy url=base+sfx+ssfx;
// should do dos2unix and concat... 
Grc grc2=in.getter().normalget(fntemp,url,0);
//Blob b;
//b.load(fnbib);
//dxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

 const StrTy & guessconcat=  GuessConcat(in , out )  ;

StrTy cmd2= " dos2unix |"+guessconcat;
IdxTy rc=filter_small_file(fnbib, fntemp, cmd2, out, 0);
inserturl(fnbib,uin,out,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

out.exit(nm);
  return 0; 

#if 0
}
guessaps()
{
uin="$1"
#http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.73.975
# http://journals.aps.org/prl/export/10.1103/PhysRevLett.73.975?type=bibtex&download=true
ssfx="?type=bibtex&download=true"
fn="$temp1"
$WGET -q -O "$fn" "$uin"
sfx=` cat "$fn" | sed -e 's/"/\n/g' | grep xport  | grep -v "[<>]" | grep "/"| head -n 1`
base=`echo "$uin" | isolatedomain`
url="$base$sfx$ssfx"
fn="$temp1"
rm -f "$fn"
citepastelynk  "$fn" "$uin" "$url"
#endif

} // guessaps


 static  IdxTy guessfuture(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

/*


base=`echo $uin | isolatedomain`
fn="$temp1"
rm -f "$fn"
$WGET --no-check-certificate -q -O "$fn"  "$uin"

# isolatehref  now returns entire utl  
# https://www.jstage.jst.go.jphttps://www.jstage.jst.go.jp/AF06S010Sho shJkuDld?sryCd=jcbn1986&noVol=15&noIssue=1&kijiCd=15_1_65&kijiLangKrke=ja&kijiTo olIdHkwtsh=AT0073&request_locale=JA

sfx=`cat "$fn"  | grep -i "bibtex\|BIB TEX"  | isolatehref| htmldecode`
url="$base$sfx"
url="$sfx"
$WGET --no-check-certificate -q -O "$fn"  "$url"

x=`grep -i "@article" "$fn" `
if [ "$x" != "" ]



https://www.jstage.jst.go.jp/article/cpb1958/33/9/33_9_3887/_pdf

v 1720 div 1975 div 1979 div 1999 div 2025 text =         
1 html 663 body 1663 div 1720 div 1975 div 1979 div 1999 div 2026 span 2027 class 2028 text = fa fa-file-o
1 html 663 body 1663 div 1720 div 1975 div 1979 div 1999 div 2029 a 2030 href 2031 text = https://www.jstage.jst.go.jp/AF06S010ShoshJkuDld?sryCd=cpb1958&noVol=33&noIssue=9&kijiCd=33_9_3887&kijiLangKrke=en&kijiToolIdHkwtsh=AT0073&request_locale=EN
1 html 663 body 1663 div 1720 div 1975 div 1979 div 1999 div 2029 a 2032 text = BIB TEX



 2080  wget -O xxx.html -S -v "https://www.jstage.jst.go.jp/article/cpb1958/33/9/33_9_3887/_article/-char/en/"
 2081  toobib -hhtml xxx.html | grep -i bib
 2082  toobib -hhtml xxx.html > xxx
 2083  vi xxx
 2084  toobib -hhtml xxx.html | grep -i bib
 2085  vi xxx
 2086  wget -O xxx.bib -S -v "https://www.jstage.jst.go.jp/AF06S010ShoshJkuDld?sryCd=cpb1958&noVol=33&noIssue=9&kijiCd=33_9_3887&kijiLangKrke=en&kijiToolIdHkwtsh=AT0073&request_locale=EN"

*/


 static  IdxTy guessjstage(const InpTy & in , OutTy & out , const IdxTy xflags=0)
{ 
const StrTy nm="guessjstage";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp",1); //
const StrTy fntemp2=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
StrTy url=MutateOnly(uin, "sed -e 's;/_[a-z].*;/_article/-char/en/;'" ,out);
// should see if uin=url
MM_ERR(MMPR3(nm,uin,url))
Grc grc1=in.getter().normalget(fntemp,url,0);
//IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
HierKvp hkvp;
Kvp m;
IdxTy lhrcp= parse_hier_html( hkvp, & m,  fntemp2,  fntemp,in,out,  0);
std::vector<StrTy> words;
//words.push_back("head");
words.push_back("a");
KvpLut lt=hkvp.lookup(words,1);
KvNodes no=hkvp.find(words,1);
MM_ERR(MMPR(lt.size()))
KvNodes noa;
MM_LOOP(ii,no) { 
if ((*ii)->matching_value("bibtex",1)) noa.push_back((*ii)); 
if ((*ii)->matching_value("BIB TEX",1)) noa.push_back((*ii)); 
}
no=noa;
MM_ERR(MMPR2(nm,noa.size()))
MM_LOOP(ii,noa)
{
Ss ss;
(*ii)->dump(ss,0);
MM_ERR(MMPR2(nm,ss.str()))
auto  * no=(*ii)->node("href",0);
if ( no)
{
StrTy ubib=no->value();
MM_ERR(MMPR3(nm,(*no).info(0),ubib))
// this worksb but it is merde 
auto  * noet=(*ii); // ->node("text",0);

if (noet) {
Ss ss;
(noet)->dump(ss,0);
MM_ERR(MMPR3(nm,(*noet).info(0),ss.str()))
MM_ERR(" looking for values ")
auto  * nov=noet->any_value(1); // node("text",0);
MM_ERR(" looking for values "<<MMPR(nov))
if ( nov)
{
ubib=nov->value();
MM_ERR(MMPR3(nm,(*nov).info(0),ubib))
if (ubib!="")
{
Grc grc1=in.getter().normalget(fnbib,ubib,0);
url=ubib;
break;
}}
}
}
else
MM_ERR(" no href wtf ")
} // ii noa 

Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);

return 0; 

//MM_ERR(" not impl")  return 0; 


} // guessjstage


 static  IdxTy guessannrev(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesstwente(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessplos(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guessplos";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
const StrTy url=MutateOnly(uin, "sed -e 's/article?id=/article\\/citation\\/bibtex?id=/'",out);
Grc grc1=in.getter().normalget(fnbib,url,0);
//IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);

return 0; 

#if 0
guessplos()
{
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119177
# http://journals.plos.org/plosone/article/citation/bibtex?id=10.1371/journal.pone.0119177
uin="$1"
url=`echo $uin | sed -e 's/article?id=/article\/citation\/bibtex?id=/' `
fn="$temp1"
citepastelynk  "$fn" "$uin" "$url" "guessplos"
handled=1
} # guessplos

#endif


} // guessplos

// https://www.worldcat.org/title/introduction-to-green-chemistry/oclc/865160441

// https://www.worldcat.org/oclc/865160441?page=easybib&client=worldcat.org-detailed_record

 static  IdxTy guessworldcat(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessworldcat";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fntemp2=out.fn("temp"); //
// TODO FIXME this can be erased when all  are rquested 
// See kluge in out.Fn( 
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
// if this is pointing to a pdf, get the summary page...
const StrTy fnsum=out.fn("temp");
StrTy u1=MutateOnly(uin, "sed -e 's/worldcat\\.org.*\\/oclc\\//worldcat.org\\/oclc\\//' ",out);
//if (u1==url) ||(u1=="")
MM_ERR(MMPR3(__FUNCTION__,uin,u1))
if (u1=="")
{

out.exit(nm);
return 0;

}
u1=u1+"?page=easybib&client=worldcat.org-detailed_record";
typename Getter::HdrMap hm;
hm["Referer"]=uin;
//Grc grc1=in.getter().normalget(fn,uin,2+128+8+32);
Grc grc1=in.getter().normalget(fnsum,u1,2+128+8+32);
BomTex bt;
Ragged j,h;
Ragged r;

//const bool xh=Bit(flags,0);
//const bool xj=!Bit(flags,1); // earlier default but bad pattern 
//const bool punt=Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); /


const IdxTy fnty=1+256; // html and json parse  only good ones 
IdxTy rc=bt.load_ragged(r,fnsum,fnty);
MM_ERR(" loaded "<<MMPR(r.size()))
std::map<StrTy,StrTy> m;
IdxTy rca=0; // bt.assemble(m,r,0);
StrTy type="article";
StrTy name="xxx"; // make from trial year or author etc. 
bt.assemble_easy_bib(m,r,0);
if (m.find("type-name")!=m.end()) type=m["type-name"];
IdxTy rcfick=synthesize(in, out, u1, nm, fnbib, name, type, m, 0);


out.exit(nm);
return 0; 

} // guessworldcat



// now they seem to have good and f-ed up doi's... 
// now also needs to accept "arXiv:" forms... 

// https://ui.adsabs.harvard.edu/abs/2007PhRvL..98y3005L/abstract
// https://arxiv.org/abs/quant-ph/0701208v2
 static  IdxTy guessarxivfour(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
//if (true) return 0; 
const StrTy nm="guessarxivfour";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fntemp2=out.fn("temp"); //
// TODO FIXME this can be erased when all  are rquested 
// See kluge in out.Fn( 
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
// if this is pointing to a pdf, get the summary page...
const StrTy fnsum=out.fn("temp");
StrTy u1=MutateOnly(uin, "sed -e 's/.pdf$//' | sed -e 's;/pdf/;/abs/;' ",out);
if (strncmp(uin.c_str(),"arXiv:",6)==0)
{
u1="https://arxiv.org/abs/"+StrTy(uin.c_str()+6);
}

MM_ERR(MMPR3(nm,uin,u1))
typename Getter::HdrMap hm;
hm["Referer"]=uin;
//Grc grc1=in.getter().normalget(fn,uin,2+128+8+32);
Grc grc1=in.getter().normalget(fnsum,u1,2+128+8+32);
BomTex bt;
Ragged j,h;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const bool only_json=true;
const IdxTy fnty=256+8+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rc=bt.parse_both(h,j,fnsum,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fnsum,h.size(),j.size()))

MM_LOOP(ii,h)
{
const auto & l=(*ii);
const IdxTy len=l.size();
if (len==0) continue;
const StrTy & x=l[len-1];
const IdxTy pos=mjm_strings::indexOf(x.c_str(),"doi");
const IdxTy posa=mjm_strings::indexOf(x.c_str(),"arXiv.");
// arXiv doi bad
//marchywka@happy:/home/documents/cpp/proj/toobib$ grep http junk/fick | grep doi1 html 243 body 340 div 747 main 748 div 754 div 758 div 785 div 789 div 858 div 862 table 865 tbody 930 tr 936 td 940 a 941 href 942 text = https://doi.org/10.48550/arXiv.quant-ph/0701208
//1 html 243 body 340 div 747 main 748 div 754 div 758 div 785 div 789 div 858 div 862 table 865 tbody 930 tr 936 td 940 a 943 text = https://doi.org/10.48550/arXiv.quant-ph/0701208
//1 html 243 body 340 div 747 main 748 div 754 div 758 div 785 div 789 div 858 div 862 table 865 tbody 1003 tr 1014 td 1017 a 1022 href 1023 text = https://doi.org/10.1103/PhysRevLett.98.253005
//1 html 243 body 340 div 747 main 748 div 754 div 758 div 785 div 789 div 858 div 862 table 865 tbody 1003 tr 1014 td 1017 a 1026 text = https://doi.org/10.1103/PhysRevLett.98.253005
if (pos==(~0)) continue;
if (posa!=(~0)) continue;
IdxTy doirc= doiorg(in, out, nm ,x,0)   ;
const IdxTy nfound=out.found();
if (!in.collect_all()) if (nfound!=nstart) break;
} // ii 
//Blob b;
//b.load(fnbib);
//IdxTy rcf= out.good_enough(b,fnbib,in,u1,nm,0);


out.exit(nm);
return 0; 

} // guessarxivthreesom


// this seems to work in one test but they also have a bibsonomy method 
// now 
// https://arxiv.org/ct?url=http://www.bibsonomy.org/BibtexHandler?requTask=upload&url=https://arxiv.org/abs/2102.11521&description=Whole-device+entanglement+in+a+65-qubit+superconducting+quantum+computer&v=b771322e
// https://arxiv.org/abs/2102.11521
// may be a general handler soon 
// this is wrong, the bibsonomy thing is an uploader and the googlsecholar thing needs more
// links 
 static  IdxTy guessarxivthreesom(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
//if (true) return 0; 
const StrTy nm="guessarxivthreepm";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fntemp2=out.fn("temp"); //
// TODO FIXME this can be erased when all  are rquested 
// See kluge in out.Fn( 
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();

typename Getter::HdrMap hm;
hm["Referer"]=uin;
Grc grc1=in.getter().normalget(fn,uin,2+128+8+32);
Blob b;
Xrc xrc= in.xformer().html_to_parsed(b,fn,0);
const StrTy linkp=
        //MutateOnly(StrTy(b), " sed -e 's/ /\\n/g' |  grep -i bibsonomy.org ",out );
        MutateOnly(StrTy(b), " sed -e 's/ /\\n/g' |  grep -i scholar.google ",out );
MM_ERR(MMPR2(fnbib,linkp))
// render and get the bibsonomy link e
Grc grc1x=in.getter().normalget(fnbib,linkp,2+128+32,hm);
//::sleep(1);
//hm["Referer"]=ucrap;
//Grc grc2=in.getter().normalget(fntemp,url,2+128,hm);
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,linkp,nm,0);
out.exit(nm);
return 0; 

} // guessarxivthreesom

 static  IdxTy guessarxivthree(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 


const StrTy nm="guessarxivthree";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fntemp2=out.fn("temp"); //
// See kluge in out.Fn( 
//const StrTy fnbib=out.fn("bibtexasdf"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
const StrTy axid=MutateOnly(uin, "sed -e 's;.*/;;g' | sed -e 's/\\.[a-z]*$//'",out);
const StrTy urla="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=arXiv:";
//const StrTy urla="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=arXiv:";
const StrTy urlb="&data_type=BIBTEX&db_key=PRE&nocookieset=1";
const StrTy url=urla+axid+urlb;
// first you need to get this , 
Grc grc1=in.getter().normalget(fntemp,url,128);
// then  lok for a  bib code, 
 // marchywka@happy:/home/documents/cpp/proj/toobib/junk$ lynx -dump -force_html xxx | grep -A 3 -i bibcode | grep "]" | sed -e 's/.*]//'
//2012arXiv1208.3684S
Blob bb;
Xrc xrc=in.xformer().html_to_rendered(bb,fntemp,0);
const StrTy bibcode=MutateOnly(StrTy(bb), " grep -A 3 -i bibcode | grep \"]\" | sed -e 's/.*]//' ",out);

const StrTy base1="http://adsabs.harvard.edu/abs/";
const StrTy base2="/exportcitation";
const StrTy url2=base1+bibcode+base2;
Grc grc2=in.getter().normalget(fntemp2,url2,128);

// IdxTy rc=filter_small_file(fnbib, fntemp2, " sed -n '/@[a-zA-Z]*{.*,/,/<\\/textarea>/p'| sed -e 's/textarea.*>@/textarea>\\n@/'| grep -v \"textarea>$\" ", out, 0);
//IdxTy rc=filter_small_file(fnbib, fntemp2, " sed -n '/@[a-zA-Z]*{.*,/,/<\\/textarea>/p'| sed -e 's/textarea.*>@/textarea>\\n@/'| grep -v \"textarea>$\" | htmldecode  ", out, 0);
IdxTy rc=filter_small_file(fnbib, fntemp2, " sed -e 's/#34//g' |  sed -n '/@[a-zA-Z]*{.*,/,/<\\/textarea>/p'| sed -e 's/textarea.*>@/textarea>\\n@/'| grep -v \"textarea>$\" | htmldecode  ", out, 0);
// wth,
//  title = #34{Magnetically Controlling the Explosion of Dirac Fermions during the Oxidation of Graphene}#34,

//IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url2,nm,0);
out.exit(nm);
return 0; 
#if 0

 wget -O xxx -S -v "http://adsabs.harvard.edu/abs/2012arXiv1208.3684S/exportcitation"
 2615  lynx -dump -force_html xxx | more
 2628  cat xxx | sed -n '/@[a-zA-Z]*{.*,/,/<\/textarea>/p'



guessarxivthree()
{

echo try guessarxivthree | smsg
uin="$1"
fn="$temp1"
axid=`echo $uin| sed -e 's;.*/;;g' | sed -e 's/\.[a-z]*$//'`
bibcode="$axid"
echo bibcode $bibcode | smsg
url="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=arXiv:$bibcode&data_type=BIBTEX&db_key=PRE&nocookieset=1"
normalget "$fn" "$url"
fn2="$temp5"
cat "$fn" | sed -n '/^@/,//p' > "$fn2"
citepaste "$fn" "$uin" "$url"
cat "$fn" | inserturl "$uin"

#endif



} // guessarxivthree

#if 0

#https://academic.oup.com/jcem/article-pdf/85/9/3121/9059437/jcem3121.pdf

        <form action="/Citation/Download" method="get" id="citationModal">^M
            <input type="hidden" name="resourceId" value="2660502" />^M
            <input type="hidden" name="resourceType" value="3" />^M
            <label for="selectFormat" class="hide">Select Format</label>^M
            <select required name="citationFormat" class="citation-download-format" id="selectFormat">^M
                <option selected disabled >Select format</option>^M
                        <option value="0" >.ris (Mendeley, Papers, Zotero)</option>^M
                        <option value="1" >.enw (EndNote)</option>^M
                        <option value="2" >.bibtex (BibTex)</option>^M
                        <option value="3" >.txt (Medlars, RefWorks)</option>^M
^M

#endif

// the pdf has a downloaded from link, 
//  more xxx.txt| grep http
//Downloaded from https://academic.oup.com/jcem/article/85/9/3121/2660502 by guest on 11 July 2021

 static  IdxTy guessoup(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 

const StrTy nm="guessoup";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fntemp2=out.fn("temp2"); //
const StrTy fnbib=out.fn("bibtex"); //
const StrTy fnbib2=out.fn("bibtex2"); //
const bool all=in.collect_all();
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
const StrTy numb=MutateOnly(uin, " sed -e 's;[^0-9]*$;;' | sed -e 's;.*/;;' ",out);
const StrTy baseurl="https://academic.oup.com/Citation/Download?resourceId=";
const StrTy endurl="&resourceType=3&citationFormat=2";
const StrTy url=baseurl+numb+endurl;
Grc grc1=in.getter().normalget(fntemp,url,2);
MM_ERR(MMPR4(fntemp,uin,numb,url))
IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
const IdxTy nfound=out.found();
if ( (nstart==nfound) || all )
{
//#https://academic.oup.com/jcem/article-pdf/85/9/3121/9059437/jcem3121.pdf
const StrTy fck=MutateOnly(uin, " sed -e 's;[-.]pdf;;g' ",out);
// should use cache with 16 but all messed up  
Grc grc1=in.getter().normalget(fn,fck,2+0*128+ 0*16);
            // <input type="hidden" name="resourceId" value="2660502" />^M
Blob blob;
blob.load(fn);
const StrTy rid=MutateOnly(blob, "grep \"resourceId\" | sed -e 's/.*value=//' sed -e 's/[^0-9]//g'  ",out);
MM_ERR(MMPR2(nm,rid))
const StrTy url2=baseurl+rid+endurl;
Grc grc2=in.getter().normalget(fntemp2,url2,2+0);
IdxTy rc=filter_small_file(fnbib2, fntemp2, "sed -n '/^@/,//p'", out, 0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url2,nm,0);

}


out.exit(nm);
return 0; 
#if 0
uin="$1"
nm="guessoup"
numb=`echo $uin| sed -e 's;/[^0-9]*$;;' | sed -e 's;.*/;;' `
#url="https://academic.oup.com/hmg/downloadcitation/$numb?format=bibtex"
url="https://academic.oup.com/Citation/Download?resourceId=$numb&resourceType=3&citationFormat=2"
normalget "$fn" "$url"
cat "$fn" | sed -n '/^@/,//p' > "$fn2"
cp "$fn2" "$fn"
# this does not get it again 
citepaste "$fn" "$uin" "$url"

#endif


} // guessoup

 static  IdxTy guessnejm(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessphysiology(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessbiomedcentralcom(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessnapedu(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessembo(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessembo";
out.enter(nm);

MM_ERR(" should use downcit... "<<MMPR(nm)); //   return 0; 
IdxTy rc= handledowncit(in,out,xflags); 
out.exit(nm);
return rc;

}

// the pdf eventually found the doi but that did not have the jounral info,

// https://www.cambridge.org/core/journals/epidemiology-and-infection/article/influence-of-age-severity-of-infection-and-coinfection-on-the-duration-of-respiratory-syncytial-virus-rsv-shedding/34472B89F9319BA11B2C905EA38474BC

// https://www.cambridge.org/core/services/aop-cambridge-core/content/view/34472B89F9319BA11B2C905EA38474BC/S0950268814001393a.pdf/influence-of-age-severity-of-infection-and-co-infection-on-the-duration-of-respiratory-syncytial-virus-rsv-shedding.pdf

 static  IdxTy guesscambridge(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guesscambridge";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
//id=`echo $uin | sed -e 's/.*\///' | sed -e 's/[#?.].*//g'`
//const StrTy muname="sed -e 's;/[^/A-Z0-9]*$;;g' | sed -e 's/.*\\///' | sed -e 's/[#?.].*//g'";
// remove extra junk for a pdf linke 
//const StrTy muname="sed -e 's;/[^/A-Z0-9]*$;;g' | sed -e 's/\\/[^/]*[a-z][^/]*//' | sed -e 's/.*\\///' | sed -e 's/[#?.].*//g'";
const StrTy muname="sed -e 's/\\/[^/]*pdf\\/[^/]*$//' | sed -e 's;/[^/A-Z0-9]*$;;g' |  sed -e 's/.*\\///' | sed -e 's/[#?.].*//g'";
const StrTy id=MutateOnly(uin, muname ,out );
const StrTy url1="https://www.cambridge.org/core/services/aop-easybib/export/?exportType=bibtex&productIds=";
const StrTy url2="&citationStyle=bibtex";
const StrTy url=url1+id+url2;
MM_ERR(MMPR4(nm,id,url,uin))
Grc grc1=in.getter().normalget(fntemp,url,0);
//const StrTy cmd=" sed -e 's/[<>]/\\n/g' | grep -A 1000 \"@article{\\|@inproceedings\" | inserturl \""+uin+"\"" ;
const StrTy cmd=" sed -e 's/ title=/\\ntitle=/'";
IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grc1.url(),nm,0);
out.exit(nm);
return 0; 

#if 0

guesscambridge()
{
# https://www.cambridge.org/core/services/aop-easybib/export/?exportType=bibtex&productIds=0323FEA9C519F454CF68D31CC14658C9&citationStyle=bibtex
# https://www.cambridge.org/core/journals/canadian-journal-of-neurological-sciences/article/coexisting-parkinsons-and-wilsons-disease-chance-or-connection/0323FEA9C519F454CF68D31CC14658C9

uin="$1"
fn="$temp1"
fn2="$temp5"
rm "$fn"
rm "$fn2"
# 2017-12-30 this fails with some suffix thing like core-reader
# hopefully the case thing works lol 
id=`echo $uin | sed -e 's/.*\///' | sed -e 's/[#?.].*//g'`
id=`echo $uin | sed -e 's;/[^/A-Z0-9]*$;;g' | sed -e 's/.*\///' | sed -e 's/[#?.].*//g'`

url="https://www.cambridge.org/core/services/aop-easybib/export/?exportType=bibtex&productIds=$id&citationStyle=bibtex"
# this should call the one that gets it 
# TODO see what bad with -Q as wrong option, where did stderr go
$WGET -O "$fn2" -q --no-check-certificate "$url"
echo changed guesscambridge to insert ctlf etc | smsg
#$WGET -O "$fn" -S -v  --no-check-certificate "$url"
cat "$fn2" | sed -e 's/ title=/\ntitle=/' > "$fn"
citepastecat  "$fn" "$uin" "$url" "guesscambridge"



#endif

//MM_ERR(" not impl")  return 0; 

} // guesscambridge


 static  IdxTy guessabo(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessmcmaster(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
// this site is an AFFCUFK 
// https://www.researchgate.net/profile/Jody-Ouweland/publication/343917910_Reduced_vitamin_K_status_as_a_potentially_modifiable_risk_factor_of_severe_COVID-19/links/5f48bab6299bf13c504629dd/Reduced-vitamin-K-status-as-a-potentially-modifiable-risk-factor-of-severe-COVID-19.pdf
 static  IdxTy guessresearchgate2(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 

const StrTy nm="guessresearchgate2";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const StrTy numb=MutateOnly(uin,
	"sed -e 's/.*publication\\/\\([0-9][0-9]*\\)[^0-9].*/\\1/'",out);
MM_ERR(MMPR2(nm,numb))
if ((numb.length()==0) ||(numb==uin))
{
MM_ERR(" numb not useful return  ")
out.exit(nm);
return 0;
}
//uid=`echo $uin| sed -e 's/.*publication\/\([0-9][0-9]*\)[^0-9].*/\1/'`

const StrTy url="https://www.researchgate.net/lite.publication.PublicationDownloadCitationModal.downloadCitation.html?fileType=BibTeX&citation=citationAndAbstract&publicationUid="+numb;
typename Getter::HdrMap hm;
hm["Referer"]=uin;
// try going through temp 
// https://www.researchgate.net/publication/223988781_Novel_pathway_of_iron-induced_blood_coagulation_Implications_for_diabetes_mellitus_and_its_complications
// https://www.researchgate.net/publication/223988781_Novel_pathway_of_iron-induced_blood_coagulation_Implications_for_diabetes_mellitus_and_its_complications/citation/download
const StrTy ucrap=uin+"/citation/download";
// why not use the file fn? 
Grc grc1=in.getter().normalget(fn,uin,2+128+8);
::sleep(1);
Grc grc1x=in.getter().normalget(fntemp,ucrap,2+128,hm);
::sleep(1);
hm["Referer"]=ucrap;
Grc grc2=in.getter().normalget(fntemp,url,2+128,hm);
//BibEntry be;
Blob b;
//static IdxTy CleanUpBadBib(Blob & b, const StrTy & fntemp ,  const StrTy & fnbib, const IdxTy flags) 
const IdxTy cubbrc=CleanUpBadBib( b, fntemp , fnbib, 0) ;
//b.load(fnbib);
if (cubbrc!=0)
{
MM_ERR(" failed to convert bib "<<MMPR3(fnbib,fntemp,cubbrc))
// doh 2022-06-20 bad exit lateer 
//out.exit(nm);
// 2022-02-05 this needs to keep trying then doh 
//return cubbrc;
}
// the researchgate bibs may not have a journal but the 
// doi will work with cross ref... 
const bool use_rg_doi=true;
bool added_doi=false;
if ( use_rg_doi)
{
BibEntryMap bm;
bm.parse(fnbib);
if (bm.size()!=0)
{
MM_ERR(" code not implemented yet to get doi of rg bibtex")
BibEntry  be=*(bm.begin()); // ((*(bem.begin())).second)[0];
if (be.name()=="article")
{
be.make_name();
be.save(fnbib,3);
}
// TODO in theory there could be an entry for doi with zero
// length but not sure it is possible... 
const StrTy doi=be.get("doi",0); // be.uin();
if (doi!="")
{
const StrTy fntemp2=out.fn("temp"); //
Grc grc=in.getter().crossref(fntemp2,doi,0);
const bool doiok=(grc.rc()=="200");
added_doi=doiok;
const StrTy fnbest=doiok?fntemp2:fnbib;
MM_ERR(MMPR3(doiok,grc.rc(),fnbest))
Blob b;
b.load(fnbest);
IdxTy rcf=doiok?  out.good_enough(b,fnbest,in,grc.url(),nm+" (xref) ",0)
: 0 ; // need to pick one place or the other fix logi here, want doi first if available. 
// : out.good_enough(b,fnbest,in,grc2.url(),nm+" (rg) ",0);
} // doi
else MM_ERR(" no doi for rg entry "<<MMPR2(uin,be.name()))
}// bm.size()

IdxTy nnow=out.found();
const bool all=in.collect_all();
MM_ERR(MMPR4(nm,nstart,nnow,all)<<MMPR(added_doi))
if (all||!added_doi)
{
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,grc2.url(),nm+" (rg) " ,0);
nnow=out.found();
}

if (nnow!=nstart) if(!in.collect_all()) { out.exit(nm); return 0;  } 

} // use_rg_doi
else { 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
}
//MM_ERR(" not impl")  
out.exit(nm);
return 0; 


} // guessresearchgate2
#if 0
guessresearchgate2()
{
h1="Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"
h2="--user-agent=\"Mozilla\" "
#h3="Referer: https://www.researchgate.net/publication/316906279_Structure_and_vibrational_spectroscopic_study_of_alpha-tocopherol" 
h3="Referer: $uin"
cook=" --load-cookies=crapcookie.txt --save-cookies=crapcookie.txt --keep-session-cookies"
uid=`echo $uin| sed -e 's/.*publication\/\([0-9][0-9]*\)[^0-9].*/\1/'`
echo researhgate2 uid=$uid | smsg
url="https://www.researchgate.net/lite.publication.PublicationDownloadCitationModal.downloadCitation.html?fileType=BibTeX&citation=citationAndAbstract&publicationUid=$uid" # 325907548
rm -f "$fn"
#echo using $WGET | smsg
# get the  cookie first 
$WGET -O "$fn" --header="Accept: $acchey" --user-agent="$uahey"   --no-use-server-timestamps $cook --no-check-certificate    "$uin"  | debugecho "getting cookie"

$WGET -O "$fn" --header="Accept: $acchey" --user-agent="$uahey"   --no-use-server-timestamps $cook --no-check-certificate  --header="$h3"    "$url"
#cat "$fn"
# this is returning a citation name "article"
rm "$temp5"
#fixcitename "$fn" "$uin" "$url" > "$temp5"
fixcitename "$fn" "$uin" "$url"| inserturl "$uin"  > "$temp5"
cp "$temp5" "$fn"
citepaste "$fn" "$uin" "$url"


#endif



 static  IdxTy guessagrisfao(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessrugnlportal(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesslibert(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessliebert";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
//const IdxTy nstart=out.found();
// cross ref seems to work with the entire url but still lol, 
// https://www.liebertpub.com/doi/pdf/10.1089/ars.2019.7803
const StrTy doi=MutateOnly(uin, 
//s"ed -e 's;.*/10\\.;10.;'| sed -e 's/[#].*//g'` ",out);
// this worked with the above link but not the right way, 
//"sed -e 's;.*doi\\.;;'| sed -e 's/[#].*//g' ",out);
// this can fail if the doi has "doi"  in it...  need non-greedy one 
"sed -e 's;.*doi[./];;'| sed -e 's/[#?].*//g' | sed -e 's;^pdf/10;10;' ",out);
IdxTy rcx=get_xref(out,in,doi,nm,0);
//Grc grc1=in.getter().normalget(fnbib,url,0);
//IdxTy rc=filter_small_file(fnbib, fntemp, "sed -n '/^@/,//p'", out, 0);
//Blob b;
//b.load(fnbib);
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0; 

#if 0
guesslibert()
{
#http://online.liebertpub.com/doi/abs/10.1089/jop.2011.0155
canhandle=`echo $uin | sed -e 's;.*/10\.;10.;'| sed -e 's/[#].*//g'`
url="http://api.crossref.org/works/$canhandle/transform/application/x-bibtex"
normalget "$fn" "$url"
something=`cat "$fn"`

#endif

MM_ERR(" not impl")  return 0; 


} // guesslibert

#if 0
https://www.ahajournals.org/doi/abs/10.1161/circ.144.suppl_1.10712
aha also has json,
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 1 publicationDate = Nov 16, 2021, 12:00:00 AM
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 2 authors = Gundry, Steven R
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 1  1 label = Acute coronary syndromes
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 1  2 searchUrl = /keyword/Acute+coronary+syndromes
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 2  1 label = Antibodies
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 2  2 searchUrl = /keyword/Antibodies
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 3  1 label = COVID-19
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 3  2 searchUrl = /keyword/COVID-19
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 4  1 label = Inflammation and inflammatory markers
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 4  2 searchUrl = /keyword/Inflammation+and+inflammatory+markers
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 5  1 label = Myocarditis
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 3 articleKeywords 5  2 searchUrl = /keyword/Myocarditis
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 4 format 1 format = ABSTRACT
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 5 assetTitle = Abstract 10712: Mrna COVID Vaccines Dramatically Increase Endothelial Inflammatory Markers and ACS Risk as Measured by the PULS Cardiac Test: a Warning
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 6 doi = 10.1161/circ.144.suppl_1.10712
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 7 issue = Suppl_1
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 8 volume = 144
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 9 sku 1 value = 0009-7322
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 9 sku 2 type = ISSN
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 9 sku 3 subtype = PPUB_IDSUBTYPE
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 9 sku 4 source = CONTENT
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 10 parentTitle = Circulation
./mjsonu.h622   json dies before parse   ok=1 eof=1 good=0 n=908484 cc=s

./mjsonu.h496  FICK  i=3 argc=4 args[i]=-quit s=-quit
./mjsonu.h280 GLOBAL 220  240  4  1  2  1  20  3 content 1 productInfo 11 viewType = Abstract
marchywka@happy:/home/documents/cpp/proj/toobib/junk$ toobib -hhtml aha.html  | toobib -json1 - | grep -v " GLOBAL = " | grep -i productInfo



#endif



 static  IdxTy guessahajournal(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
// https://www.ahajournals.org/doi/10.1161/circ.140.suppl_1.16393
//  wget -O xxx -S -v "https://www.ahajournals.org/action/downloadCitation?doi=10.1161%2Fcirc.140.suppl_1.16393&downloadFileName=aha_circ140_A10613&include=abs&format=bibtex&direct=&submit=Download"
// https://www.ahajournals.org/action/showCitFormats?doi=10.1161%2Fcirc.140.suppl_1.16393
const StrTy nm="guessahajournal";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy url="http://ahajournals.org/action/downloadCitation";
//StrTy cmd=" urldecode | sed -e 's/.*\\(10\\.[0-9][0-9][0-9][0-9][0-9a-zA-Z.]*\\/[-0-9a-zA-Z.]*\\).*/\\1/' ";
//StrTy cmd=" urldecode | sed -e 's;.*doi[=/];;' ";
//https://www.ahajournals.org/doi/abs/10.1161/circ.144.suppl_1.10712
StrTy cmd=" urldecode | sed -e 's;.*doi[=/];;'| sed -e 's/^[^1][^1]*//' ";
const StrTy name=MutateOnly(uin, cmd ,out );
// this seems to work, 
// https://www.ahajournals.org/action/downloadCitation?doi=10.1161%2Fcirc.140.suppl_1.16393&downloadFileName=aha_circ140_A10613&include=abs&format=bibtex&direct=&submit=Download

//const StrTy url=base+name;
//Grc grc1=in.getter().normalget(fntemp,url,0);
const StrTy pd="doi="+name+"&include=abs&format=bibtex&direct=&submit=Download";
MM_ERR(MMPR3(name,pd,cmd))
//Grc grc1=in.getter().normalpost(fntemp,url,pd,0);
Grc grc1=in.getter().normalget(fntemp,url+"?"+pd,0);

cmd=" sed -e 's/^@\\(.*\\)\\\\/@\\1/g'" ;
StrTy cmd2=" sed -e 's/^doi =\\(.*\\)\\\\/doi =\\1/g'" ;
// I think this is a bug in their code, no comma after doi entry 
StrTy cmd3=" sed -e 's/^doi =\\(.*\\)}$/doi =\\1},/'" ;
cmd=cmd+"|"+cmd2+"|"+cmd3;
IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
//inserturl(fnbib,uin,out,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
out.exit(nm);
return 0; 
} // guessahajournals
#if 0
guessahajournal()
{
echo guessahajournal | smsg
nm="guessahajournal"
uin="$1"
# apparently gets work well too 
#http://ahajournals.org/action/downloadCitation?doi=10.1161/ATVBAHA.108.172072&include=abs&format=bibtex&direct=&submit=Download

url="http://ahajournals.org/action/downloadCitation"
doi=`echo $uin| urldecode | sed -e 's/.*\(10\.[0-9][0-9][0-9][0-9][0-9a-zA-Z.]*\/[-0-9a-zA-Z.]*\).*/\1/' `
pd="doi=$doi&include=abs&format=bibtex&direct=&submit=Download"
url="$url?$pd"
fn="$temp1"
rm -f "$fn"
fn2="$temp5"
#$WGET -q  -O "$fn" "$url" --post-data="$pd"
#$WGET -S -v  --no-check-certificate  -O "$fn" "$url" --post-data="$pd"
#normalpost "$fn" "$url" "$pd"
normalget "$fn" "$url"

#endif 

 static  IdxTy guessatsjournal(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessarchpath(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }


/*
Toobib fails because it won't let wget download,
// https://www.jbc.org/article/S0021-9258(17)43386-2/pdf

 strings chromate_downloads/PIIS0021925817433862.pdf  | grep -i doi
         <dc:identifier>doi:10.1016/S0021-9258(17)43386-2</dc:identifier>
               <rdf:li xml:lang="x-default">Journal of Biological Chemistry, 259 (1984) 2537-2540. doi:10.1016/S0021-9258(17)43386-2</rdf:li>
         <prism:doi>10.1016/S0021-9258(17)43386-2</prism:doi>
         <prism:url>http://dx.doi.org/10.1016/S0021-9258(17)43386-2</prism:url>
         <crossmark:DOI>10.1016/S0021-9258(17)43386-2</crossmark:DOI>
         <pdfx:doi>10.1016/S0021-9258(17)43386-2</pdfx:doi>
// no need to scrape, it is in the exif... 


*/
 static  IdxTy guessjbc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
// https://www.jbc.org/article/S0021-9258(17)43386-2/pdf
const StrTy nm="guessjbc";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy base="10.1016/";;
const StrTy muname="sed -e 's/.pdf$//' | sed -e 's/.*\\///'";
const StrTy name=MutateOnly(uin, muname ,out );
MM_ERR(MMPR4(nm,uin,muname,name))
const StrTy doi=base+name;
IdxTy rcx=get_xref(out,in,doi,nm,0);
out.exit(nm);
return 0; 

} // guessjbc
/*
 2029  wget -O xxx -S -v "https://cds.cern.ch/record/1108310/files/p75.pdf"
 2034  wget -O zzz -S -v "https://cds.cern.ch/record/1108310"
*/

 static  IdxTy guesscern(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
// https://www.jbc.org/article/S0021-9258(17)43386-2/pdf
const StrTy nm="guesscern";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy base="10.1016/";;
const StrTy muname="sed -e 's/\\/files.*//'";
const StrTy name=MutateOnly(uin, muname ,out );
MM_ERR(MMPR4(nm,uin,muname,name))
if (name!=uin)
{
IdxTy rcr=Recurse(in,out,name);
MM_ERR(MMPR3(rcr,nm,out.found()))
} 

out.exit(nm);
return 0; 

} // guessjbc


 static  IdxTy guessmicrores(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessgenetics(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesssemantic(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guesssemantic";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy base="https://semanticscholar.org/paper/";
const StrTy muname="sed -e 's/.*\\.org\\///' | sed -e 's/\\..*//' | sed -e 's;/;;g'";
const StrTy name=MutateOnly(uin, muname ,out );
const StrTy url=base+name;
Grc grc1=in.getter().normalget(fntemp,url,0);
//const StrTy cmd=" sed -e 's/[<>]/\\n/g' | grep -A 1000 \"@article{\\|@inproceedings\" | inserturl \""+uin+"\"" ;
const StrTy cmd=" sed -e 's/[<>]/\\n/g' | grep -A 1000 \"@article{\\|@inproceedings\" " ;
IdxTy rc=filter_small_file(fnbib, fntemp, cmd, out, 0);
inserturl(fnbib,uin,out,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
const IdxTy nnow=out.found();
if (!all)  if (nnow!=nstart) { out.exit(nm); return 0; }
{
Grc grc=in.getter().normalget(fn,uin,16);
//const StrTy cmd2= "sed -e 's/[<>]/\\n/g' |sed -e 's/@/\\n@/g' |  grep -A 1000 "@article{\|@inproceedings" | inserturl \""+uin+"\""; 
const StrTy cmd2= "sed -e 's/[<>]/\\n/g' |sed -e 's/@/\\n@/g' |  grep -A 1000 \"@article{\\|@inproceedings\" "; 
IdxTy rc2=filter_small_file(fnbib, fn, cmd, out,  0);
inserturl(fnbib,uin,out,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
}
//cat "$fn" | sed -e 's/[<>]/\n/g' |sed -e 's/@/\n@/g' |  grep -A 1000 "@article{\|@inproceedings" | inserturl "$uin"  > "$fn2"
out.exit(nm);
return 0; 
#if 0
guesssemanticnu()
url="$uin"
normalget "$fn" "$url"

cat "$fn" | sed -e 's/[<>]/\n/g' |sed -e 's/@/\n@/g' |  grep -A 1000 "@article{\|@inproceedings" | inserturl "$uin"  > "$fn2"
echo twy2 handled= $handled check "$fn" if messed up | smsg

guesssemantic()
{
uin="$1"
base="https://semanticscholar.org/paper/"
name=`echo $uin | sed -e 's/.*\.org\///' | sed -e 's/\..*//' | sed -e 's;/;;g'`
url=$base$name
normalget "$fn" "$url"
# this is supposed to insert the url into the bibtex but fails.
#cat "$fn" | sed -e 's/[<>]/\n/g' | grep -A 1000 "@article{\|@inproceedings" | awk -v v="$uin" 'BEGIN{p=1}{ if ( $1=="@article") p=1 ; if ( $1=="}") print ",url={ "v"}"; if ( p!=0) print $0; if ($1=="}") p=0; }' > "$fn2"
cat "$fn" | sed -e 's/[<>]/\n/g' | grep -A 1000 "@article{\|@inproceedings" | inserturl "$uin"  > "$fn2"


#endif



} // guesssemantic

 static  IdxTy guessdiabetesjournal(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guessfasebj(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

#if 0 
this seems to still work but did not verify mayube  bad name, 
# 648  wget -O xxx -S -v --no-check-certificate --post-file=p2.txt --header="X-Requested-With: XMLHttpRequest" --user-agent="Mozilla" --header="Referer: https://www.cabdirect.org/cabdirect/abstract/19360800757" --header="Upgrade-Insecure-Requests: 1"  "https://www.cabdirect.org/cabdirect/utility/cd4mycabiajaxhandler/"
#  650  wget -O xxx -S -v --no-check-certificate --post-file=p2.txt  "https://www.cabdirect.org/cabdirect/utility/cd4mycabiajaxhandler/"
cabdata()
{
cdnum="$1"
# cat p2.txt  | urldecode
cat - << ----XXX---
methodData={"method":"downloadRecords","items":[{"type":"AbstractMarkLogic","itemUrls":["/cabdirect/abstract/$cdnum"]}],"recordSource":"SelectedRecords","pageSource":"unknown","numberRange":"","pageUrl":"https://www.cabdirect.org/cabdirect/abstract/$cdnum"}&userPrefs={"format":"BibTeX","downloadTarget":"DownloadFile","portions":"CitationOnly","destination":"EmailBodyText","exportEmail":"","SavedRecordsPageSize":50}
----XXX---
}

guesscab()
{
uin="$1"
url="https://www.cabdirect.org/cabdirect/utility/cd4mycabiajaxhandler/"
cdnum=`echo $uin| sed -e 's/.*[^0-9]\([0-9][0-9]\{6,20\}\)/\1/' `
echo guesscab  $cdnum | smsg
pd=`cabdata "$cdnum"|urlencode`
fn="$temp1"
fn2="$temp5"
rm -f "$fn"
rm -f "$fn2"
cabdata "$cdnum"|urlencode| sed -e 's/%3D/=/g'| sed -e 's/%26/\&/g' > "$fn2"
normalpostf "$fn" "$url" "$fn2"
#citepastecat  "$fn" "$uin" "$url " "guesscab $pd"
citepastecat  "$fn" "$uin" "$url " "guesscab"
echo handled $handled | smsg
#handled=1
#echo force handled for now | smsg

} # guesscab 


#endif



 static  IdxTy guesscab(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }




 static  IdxTy guessblood(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesshathi(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesscelt(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesspubag(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

#if 0



#https://www.jstor.org/stable/3282360?seq=1#page_scan_tab_contents
guessjstor()
{
nm=guessjstor
uin="$1"
echo $nm | smsg
fn="$temp1"
fndest="$temp2"
#getriscite "$url" "$uin" "$fn" "$fndest"
#fixcitename "$fndest" "$uin" "$url"| inserturl "$uin"  > "$temp5"
#cp "$temp5" "$fn"
normalget "$fn" "$uin"
doid=`cat "$fn" | grep /citation/info/ | head -n 1 | sed -e 's/.*citation\/info\///' | sed -e 's/".*//g' `
url="https://www.jstor.org/citation/text/$doid"
echo doid = $doid, url=$url | smsg
normalget "$fn" "$url"
citepaste "$fn" "$uin" "$url"



https://www.jstor.org/citation/text/10.2307/20170301
https://www.jstor.org/stable/20170301
#endif


static IdxTy guessjstor(const InpTy & in , OutTy & out , const IdxTy xflags=0) 
{ 
// These  blocked the download but can mutate the 
// url at least, Zotero works  
const StrTy nm="guessjstor";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const bool fiucking_sht = false; 
if (!fiucking_sht)
{
const StrTy url=MutateOnly(uin, "sed -e 's/[?#].*//' | sed -e 's;stable;citation/text/10.2307;'| grep citation " ,out );
MM_ERR(MMPR3(nm,url,uin))
// if the url==uin, the isbib thing should pick it up... 
//if ((url!=uin) &&(url.length()!=0))
if ((url.length()!=0))
{ // try friendly ua
Grc grc1=in.getter().normalget(fnbib,url,1|(1<<10));
Blob b;
b.load(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 


} // !fiucking_sht


if (fiucking_sht)
{
Grc grc1=in.getter().normalget(fn,uin,16);
Blob b;
b.load(fn);
StrTy page=StrTy(b);
const StrTy url=MutateOnly(page, "grep /citation/info | head -n 1 | sed -e 's/.*citation\\/info\\///' | sed -e 's/\".*//g'" ,out );
MM_ERR(MMPR3(nm,url,uin))
if ((url!=uin) &&(url.length()!=0))
{
Grc grc1=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 

} // fiucking_sht;

out.exit(nm);
return 0;
} // guessjstor

/*

Can either try to get xml link or make it up, 
 2659  wget -O jaxk.html -S -v "https://www.authorea.com/doi/full/10.22541/au.159103680.00306295"
 2663  wget -O jaxk.xxx -S -v "https://www.authorea.com/doi/xml/10.22541/au.159103680.00306295/v1"
 2882  history | grep wget
 2883  wget -O xxx -S -v "https://www.authorea.com/doi/xml/10.22541/au.159103680.00306295"
 2884  wget -O xxx -S -v "https://www.authorea.com/doi/abs/10.22541/au.159103680.00306295"
 2885  wget -O xxx -S -v "https://www.authorea.com/doi/full/10.22541/au.159103680.00306295"
 2891  history | grep wget



*/

// https://ufdc.ufl.edu/AA00004956/00001
//http://ufdcimages.uflib.ufl.edu/AA/00/00/49/56/00001/AA00004956.pdf
static IdxTy guessufl(const InpTy & in , OutTy & out , const IdxTy xflags=0) 
{
const StrTy nm="guessufl";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
//const StrTy url=MutateOnly(uin, "sed -e 's/\\/doi\\/[^/][^/]*\\/10/\\/doi\\/xml\\/10/'  " ,out );
const StrTy url=MutateOnly(uin, "grep ufdcimages | sed -e 's;[/.]; ;g;' |awk '{print \"https://ufdc.ufl.edu/\"$(NF-1)\"/\"$(NF-2)}' " ,out );

//static IdxTy handlexmlformats(const InpTy & in , OutTy & out , const IdxTy xflags=0)
const bool url_is_ok=(url.length()&&(url!=uin));
MM_ERR(MMPR(url_is_ok)<<MMPR4(nm,fn,uin,url))
if (!url_is_ok)
{
// may not fetch? 
IdxTy rch=handlexmlformats(in,out,xflags);
}
else
{
const StrTy nurl=url;
InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rch=handlexmlformats(in2,out2,xflags);
MM_ERR(" return kluge ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);
} // url_is_ok

out.exit(nm);
return 0;

} // guessufl

static IdxTy guessauthorea(const InpTy & in , OutTy & out , const IdxTy xflags=0) 
{ 
// These  blocked the download but can mutate the
// url at least, Zotero works 
const StrTy nm="guessaurhoterea";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const bool fiucking_sht = false; 
if (!fiucking_sht)
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
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

const IdxTy nnow=out.found();
if (!all) if ( nstart<nnow) { out.exit(nm); return 0; } 
//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 


} // !fiucking_sht


if (fiucking_sht)
{
Grc grc1=in.getter().normalget(fn,uin,16);
Blob b;
b.load(fn);
StrTy page=StrTy(b);
const StrTy url=MutateOnly(page, "grep /citation/info | head -n 1 | sed -e 's/.*citation\\/info\\///' | sed -e 's/\".*//g'" ,out );
MM_ERR(MMPR3(nm,url,uin))
if ((url!=uin) &&(url.length()!=0))
{
Grc grc1=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 

} // fiucking_sht;

out.exit(nm);
return 0;
} // guessauthorea



// https://digitalcommons.odu.edu/biomedicalsciences_etds/63/?utm_source=digitalcommons.odu.edu%2Fbiomedicalsciences_etds%2F63&utm_medium=PDF&utm_campaign=PDFCoverPages
// https://digitalcommons.odu.edu/cgi/viewcontent.cgi?article=1063&context=biomedicalsciences_etds
//pdftotext -layout  xxx.pdf xxx.txt
//marchywka@happy:/home/documents/cpp/proj/toobib/junk$ grep http xxx.txtFollow this and additional works at: https://digitalcommons.odu.edu/biomedicalsciences_etds
//https://digitalcommons.odu.edu/biomedicalsciences_etds/63


 static  IdxTy guessdigitalcommons(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessdigitalcommons";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex",1); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const bool fiucking_sht = false; 

//const StrTy cmd= "grep /citation/info | head -n 1 | sed -e 's/.*citation\\/info\\///' | sed -e 's/\".*//g'" ;
//const StrTy url=MutateOnly(url,cmd,out );

Grc grc=in.getter().normalget(fn,uin,16);
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",1);
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) 
{ MM_ERR(" pdftotext did not generate output ") 
out.exit(nm);
return 1;
}
const StrTy cmd="awk '{print $1}' | grep http | grep \"digitalcommons\\..*\\.edu\" | head -n 1 ";
StrTy newlink=ExtractOnly(dest,cmd,out);
MM_ERR(MMPR4(nm,newlink,cmd,dest))
if (newlink=="") { out.exit(nm); return 2; } // newlink 

IdxTy rcr=Recurse(in,out,newlink);
//nnow=out.found();
// this is not working wtf??? 
if (!false) { 
// 
const IdxTy nnow=out.found();
for (IdxTy ib=nstart; ib<nnow; ++ib)
{ // should check type for phdthesis but ok for now. 
const StrTy _fnbib=out.result_file(ib);
// links may be missing 
 std::map<StrTy, StrTy> m;
 m["school"]="dissertation_institution";
 m["school"]="publisher";
// m["type"]="phdthesis";
 out.insert_missing(_fnbib,m,1);
} // nnow!=nstart
//} // newlink 
} 


out.exit(nm);
return 0;

} // guessdigitalcommons




 static  IdxTy guessuridc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }



 static  IdxTy guesstracetenedu(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesskidint(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesselife(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { 
const StrTy nm="guesselife";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
const StrTy url=MutateOnly(uin," sed -e 's/\\.[a-a][a-z]*$//'| awk '{print $0\".bib\"}'",out);
Grc grc1=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
const IdxTy nnow=out.found();
if( !in.collect_all())  if (nnow!=nstart) { out.exit(nm); return 0; } 

Grc grc2=in.getter().normalget(fn,uin,16);
//AccessTy html_to_parsed(Blob & d, const StrTy & fn, const IdxTy flags)
b.clear();
Xrc xrc=in.xformer().html_to_rendered(b,fn,0);
MM_ERR(MMPR(StrTy(b)))
const StrTy doi=MutateOnly(StrTy(b), " grep -i doi | grep -i Cite | sed -e 's/.*]//' | grep \"[0-9]\" |  head -n 1 | awk '{print $NF}' ",out );
IdxTy rcx=get_xref(out,in,doi,nm,0);
#if 0
probably obserolete, just add bib afacit 
guesselife()
{
nm=guesselife
uin="$1"
fn="$temp1"
normalget "$fn" "$uin"
doi=`lynxhtml "$fn" | grep -i doi | grep -i Cite | head -n 1 |  sed -e 's/.*]//' | awk '{print $NF}' `
#doi=`lynx -dump --force-html "$fn" | grep -i doi | grep -i Cite | head -n 1 |  sed -e 's/.*]//' | awk '{print $NF}' `
# dest doi
getcrossref "$nm" "$doi" "$fn" "$uin"


#endif


out.exit(nm); // don't exit if we are mutating... 
return 0; 

}
 static  IdxTy guessjlr(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

/*
mv "$fn" "$fn.gz"
gunzip "$fn.gz"
doi=`cat "$fn" | doisfromtext `
echo $nm trying doi  $doi  | smsg
getcrossref "$nm" "$doi" "$fn" "$uin"
#normalget "$fn" "$url" 
#citepastecat  "$fn" "$uin" "$url" "$nm"

} # guessfcklww()


*/

// https://journals.lww.com/soilsci/Abstract/1974/02000/Effect_of_Pore_Size_on_Diffusion_Coefficients_in.2.aspx
 static  IdxTy guessfcklww(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessfcklww" ; // "guesssmdpi";
out.enter(nm);
const bool all=in.collect_all();
const IdxTy nstart=out.found();
const StrTy & uin=in.uin();
const StrTy _fn=out.fn(); //
// at this point mutating the pdf url would be good but then
// the file name can't pollute the cache. 
const StrTy fntemp=out.fn("temp"); //
const bool mut=(strncmp(uin.c_str()+uin.length()-3,"pdf",3)==0);
const StrTy fn=mut?fntemp:_fn;
const StrTy _uin=mut?MutateOnly(uin, " sed -e 's/[.\\/]pdf$//'" ,out ):uin;
Grc grc=in.getter().normalget(fn,_uin,16);

//MM_ERR(" not impl")  return 0; 

out.exit(nm);
return 0;

} // guessfcklww

 static  IdxTy guessthoracicsurgery(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }
 static  IdxTy guesscmajca(const InpTy & in , OutTy & out , const IdxTy xflags=0)  { MM_ERR(" not impl")  return 0; }

 static  IdxTy guessdora(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessdora";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const StrTy url=MutateOnly(uin, "sed -e 's/\\/datastream\\/.*//' " ,out );
if ((in.depth()<5)&&(url!=uin) &&(url.length()!=0))
{
IdxTy rcr=Recurse(in,out,url);
}
out.exit(nm);
return 0;
} // guessdora


 static  IdxTy guessfrontiersin(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessfrontiersin";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const StrTy url=MutateOnly(uin, "sed -e 's/\\([0-9]\\)\\/[a-z]*$/\\1\\/bibtex/'" ,out );
if ((in.depth()<5)&&(url!=uin) &&(url.length()!=0))
{
Grc grc1=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 

out.exit(nm);
return 0;

} // guessfrontiersin


// http://cogprints.org/8948/1/2012-4-19.pdf
//http://cogprints.org/cgi/export/eprint/9825/BibTeX/cogprints-eprint-9825.bib

 static  IdxTy guesscogprints(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesscogprints";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
//const bool all=in.collect_all();
const StrTy code=MutateOnly(uin, "sed -e 's/.*org\\///' | sed -e 's;/.*;;g'" ,out );
const StrTy  url="http://cogprints.org/cgi/export/eprint/"+code+"/BibTeX/cogprints-eprint-"+code+".bib";
if ((url!=uin) &&(url.length()!=0))
{
Grc grc1=in.getter().normalget(fnbib,url,0);
Blob b;
b.load(fnbib);
// this will fix dos2unix and remove blank lines 
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
} // mutate worked 

out.exit(nm);
return 0;
} // guessdora


// grep http xxx.txt | more
//Follow this and additional works at: https://lib.dr.iastate.edu/rtd
//https://lib.dr.iastate.edu/rtd/8306
//marchywka@happy:/home/documents/latex/proj/paradox/junk$ wget -O xxx -S -v "https://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=9305&context=rtd"
// for the rtd link, render to text and find the first matching link 

 static  IdxTy guessiastatertd(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessiastatertd";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const IdxTy nstart=out.found();
Grc grc=in.getter().normalget(fn,uin,16);
StrTy type;
IdxTy tci=identify(type,fn,"",out,0);
//pdftotext -q  "$fn" "$fnpdf" 2>&1 > "$temp_pdf_fick"
if (type!="pdf")
{ MM_ERR(" may not be pdf "<<MMPR(type)) }
const StrTy dest=out.fn("pdftotext",0);
//Scrapes v;
Xrc xrc=in.xformer().pdf_to_text(dest,fn,0);
if (!in.getter().exists(dest)) 
{ MM_ERR(" pdftotext did not generate output ") 
out.exit(nm);
return 1;

}
//else {
//IdxTy rc2=filter_small_file(fnbib, fn, cmd, out,  0);
// grep http xxx.txt | more
//Follow this and additional works at: https://lib.dr.iastate.edu/rtd
//https://lib.dr.iastate.edu/rtd/8306
const StrTy cmd="grep http | grep \"lib.dr.iastate.edu\" | grep \"/rtd/[0-9][0-9]*\" | head -n 1 ";
StrTy newlink=ExtractOnly(dest,cmd,out);
MM_ERR(MMPR3(newlink,cmd,dest))
if (newlink=="") { out.exit(nm); return 2; } // newlink 
// const IdxTy rc=handlegsmeta(in,out,xflags,0);

/*
The newllink has better gsmeta info, 
 0 html  1 head  2 meta  3 name 4 text bepress_citation_dissertation_institution
 0 html  1 head  2 meta  3 name 4 text bepress_citation_dissertation_name
 0 html  1 head  2 meta  3 name 4 text bepress_citation_author
 0 html  1 head  2 meta  3 name 4 text bepress_citation_author_institution
 0 html  1 head  2 meta  3 name 4 text bepress_citation_title
 0 html  1 head  2 meta  3 name 4 text bepress_citation_date
 0 html  1 head  2 meta  3 name 4 text bepress_citation_pdf_url
 0 html  1 head  2 meta  3 name 4 text bepress_citation_abstract_html_url
 0 html  1 head  2 meta  3 name 4 text bepress_citation_publisher
 0 html  1 head  2 meta  3 name 4 text bepress_citation_doi
 0 html  1 head  2 meta  3 name 4 text bepress_citation_online_date

IdxTy rch=handlehtmlparse(href,nm,in,out,fn,flags);
if (href.length()==0) { out.exit(nm); return 1; }
IdxTy rcfm=FixMeta(href, in,out,nm,uin);
out.exit(nm);
return 0;

*/

// should find the doi or meta
IdxTy nnow=out.found();
IdxTy nstartx=nnow;
// this needs to make the type thesis and change the bepress prefix..  
if (!true)
{
StrTy href;
// TODO  this creates a problem for recurse it has to get this twice... 
// making a url cache may be better but still issue with cookie etc
// could just ask for flush. 
const StrTy fnnew=out.fn("temp"); //
Grc grc=in.getter().normalget(fnnew,newlink,0);
IdxTy rch=handlehtmlparse(href,nm,in,out,fnnew,0);
//if (href.length()==0) { out.exit(nm); return 1; }
if (href.length()!=0)
{
IdxTy rcfm=FixMeta(href, in,out,nm+" ( meta ) ",uin);
} // href 
nnow=out.found();
nstartx=nnow;
} // true 
if (nnow!=nstart) if (!in.collect_all()) { out.exit(nm); return 0; }
IdxTy rcr=Recurse(in,out,newlink);
nnow=out.found();
// 
for (IdxTy ib=nstartx; ib<nnow; ++ib)
{ // should check type for phdthesis but ok for now. 
const StrTy _fnbib=out.result_file(ib);
// links may be missing 
 std::map<StrTy, StrTy> m;
 m["school"]="publisher";
 out.insert_missing(_fnbib,m,1);
} // nnow!=nstart
//} // newlink 


//} // have the text output 

//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
//const IdxTy nstart=out.found();
//const bool all=in.collect_all();
//const StrTy url=MutateOnly(uin, "sed -e 's/\\([0-9]\\)\\/[a-z]*$/\\1\\/bibtex/'" ,out );
//if ((in.depth()<5)&&(url!=uin) &&(url.length()!=0))
//{
//Grc grc1=in.getter().normalget(fnbib,url,0);
//Blob b;
//b.load(fnbib);
// this will fix dos2unix and remove blank lines 
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
//} // mutate worked 

out.exit(nm);
return 0;
} // guessiastatertd

// https://www.nber.org/papers/w28587
// https://www.nber.org/system/files/working_papers/w28587/w28587.pdf
 static  IdxTy guessnber(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessnber";
out.enter(nm);
const StrTy & uin=in.uin();
//const StrTy fn=out.fn(); //
//const StrTy fntemp=out.fn("temp"); //
//const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// two approaches- get headless and find doi OR use highwire on .citation-tools
// this is messing up bmj.n16135.short wait no it isnt' doh 
const StrTy numb=MutateOnly(uin, "sed -e 's/[^a-zA-Z0-9]/\\n/g' | grep \"^w[0-9][0-9]*$\" |   head -n 1 " ,out );

const StrTy url="https://www.nber.org/papers/"+numb;
MM_ERR(MMPR4(nm,uin,numb,url))
if (numb.length()==0) return 0;
if (url==uin)
{
MM_ERR(" have nirvana let the handlers find it doh ")
return 0; 
}
IdxTy rc=Recurse(in,out,url);
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
//IdxTy rc=handlehighwire(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from "<<nm<<" recurions adopting now ... "<<MMPR(url))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
IdxTy nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rc; } 
//Blob b;
//b.load(fnbib);
// this will fix dos2unix and remove blank lines 
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
//} // mutate worked 

out.exit(nm);
return 0;
} // guessnber

// this link needs to be fetched, 
// https://www.jocms.org/index.php/jcms/article/view/1100/499
// to scrape the last bit of this, 
//	<a href="https://www.jocms.org/index.php/jcms/article/download/1100/499/2212" class="download" download="">

/*
the code works but wtf this is easier just remove the last entry lol, 

 2371  wget -O xxx -S -v "https://www.jocms.org/index.php/jcms/article/view/1100/499"
 2372  more xxx
 2373  more xxx | grep -i doi
 2374  wget -O xxx -S -v "https://www.jocms.org/index.php/jcms/article/view/1209"
 2375  history
 2376  wget -O xxx -S -v "https://www.jocms.org/index.php/jcms/article/view/1100"
 2377  more xxx


*/



static IdxTy guessjocms(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guessjocms";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
const StrTy turl=MutateOnly(uin, "grep download" ,out );
if (turl!="") { 
MM_ERR(" exitind "<<MMPR(turl))
out.exit(nm); return 0; } 
// first try deleting the second termina number lol
const StrTy snurl=MutateOnly(uin, "sed -e 's;/\\([0-9][0-9]*\\)/[0-9][0-9]*$;/\\1;'" ,out );
MM_ERR(MMPR3(__FUNCTION__,uin,snurl))
// no longer recursing, url is already ok. 
//if (snurl==uin) { out.exit(nm); return 0; } // recursing and download wont wor
if( snurl!="") //  if (snurl!=uin)
{
// their doi entry is messed  up missing author but adhoc html is ok 
//IdxTy rcsn=Recurse(in,out,snurl);
InpTy in2(in,snurl);
OutTy out2(out,snurl);
// skip DC
IdxTy rc=handleadhochtml(in2,out2,1); // in.mom()->Guess(in2,out2,in.rflags());
//MM_ERR(" return from highwire recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);

//handleadhochtml
} 
if (!all) 
{const IdxTy nnow=out.found();
if (nstart!=nnow)  { out.exit(nm); return 0; } 
} // all 
if (true)  { out.exit(nm); return 0; } 
// needs "1" to dump the dom doh although could also work with
// getting the pdf which seems to occur by accident now... 
//Grc grc=in.getter().headlessget(fn,in.uin(),1+(show_trial?32:0));
Grc grc=in.getter().headlessget(fntemp,in.uin(),33);
// does not really need json but still 
BomTex bt;
Ragged j,h;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const bool only_json=!true;
const IdxTy fnty=256+8+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rcp=bt.parse_both(h,j,fntemp,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))
StrTy nurl="";
MM_ERR(MMPR(h.dump()))
MM_LOOP(ii,h)
{
const IdxTy len=(*ii).size();
if (len==0) continue;
const StrTy & str=(*ii)[len-1];
const IdxTy pos=mjm_strings::indexOf(str.c_str(),"download");
const IdxTy pos2=mjm_strings::indexOf(str.c_str(),"jocms.org");
if (pos!=(~0))if(pos2!=(~0)) { nurl=str;  break; }

} // ii 
MM_ERR(MMPR(nurl))

if (nurl=="") { out.exit(nm); return 0; }

// 1 html 85 body 89 header 111 a 112 href 113 text = https://www.jocms.org/index.php/jcms/article/download/1100/499/2212

// this should return a pdf but need to go get it, for now just let
// normal handlers work. 

//MM_ERRF(" trying new chromate download  "<<MMPR2(fn,StrTy(in.uin())))
//in.getter().clear(fn,0);
//Grc grc=in.getter().headlessdownload(fn,in.uin(),16+(show_trial?32:0));
//const IdxTy fnsz=in.getter().size(fn,0);
//headless_failed=(fnsz<1)||(IdxTy(grc)!=0);

IdxTy rcr=Recurse(in,out,nurl);


out.exit(nm);
return 0;
} // guessjocms

// 2022-10-29 started
// search this, 
// http://apps.who.int/iris/handle/10665/43407?search-result=true&query=0300-0923&scope=%2F&rpp=10&sort_by=score&order=desc
// using scraped ISSN from here, 
// https://apps.who.int/iris/rest/bitstreams/1425664/retrieve

static IdxTy guesswhoiris(const InpTy & in, OutTy & out , const IdxTy xflags=0)  
{
const StrTy nm="guesswhoiris";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
//const StrTy turl=MutateOnly(uin, "grep download" ,out );
//if (turl!="") { 
//MM_ERR(" exitind "<<MMPR(turl))
//out.exit(nm); return 0; } 
// first try deleting the second termina number lol


out.exit(nm);
return 0;
} // guesswhoiris



 static  IdxTy guessbmj(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessbmj";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// two approaches- get headless and find doi OR use highwire on .citation-tools
// this is messing up bmj.n16135.short wait no it isnt' doh 
const StrTy nurl=MutateOnly(uin, "sed -e 's/\\([0-9]\\)[.?#;%].*/\\1/g' | awk '{print $0\".citation-tools\"}'" ,out );
InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rc=handlehighwire(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from highwire recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);
IdxTy nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rc; } 
// try the headlessthing
//const StrTy code=MutateOnly(uin, "sed -e 's/.*org\\///' | sed -e 's;/.*;;g'" ,out );
//const StrTy  url="http://cogprints.org/cgi/export/eprint/"+code+"/BibTeX/cogprints-eprint-"+code+".bib";
//if ((url!=uin) &&(url.length()!=0))
//{
Grc grc1=in.getter().headlessget(fn,uin,7);
// this renders nicely although it may be scrapable as-is./
Blob b;
Xrc xrc=in.xformer().html_to_rendered(b,fn,0);
StrTy s=StrTy(b);
// picking up crap at end of line but crossref still works... 
const StrTy doiold=MutateOnly(s, "grep doi\\.org | sed -e 's/.*org.10/10/' | head -n 1  " ,out );
// 2021-07-25
const StrTy doi=MutateOnly(s, "grep doi\\.org | sed -e 's/.*org.10/10/' | sed -e 's/ .*//'  | head -n 1  " ,out );
MM_ERR(" changed to eliminate spaces "<< MMPR3(s,doi,doiold))
IdxTy rcx=get_xref(out,in,doi,nm+" ( doi )" ,0);
 nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rc; } 
StrTy href;
 IdxTy  hrcp=handlehtmlparse( href, nm+" ( parser )" ,  in ,  out,  fn,  xflags );

if ( href.length()==0)  { out.exit(nm); return rc; } 
IdxTy rcfm=FixMeta(href, in,out,nm+" ( parser )",uin);


// also has meta
// but it also seems to have highwire. 
// http://dx.doi.org/10.1136/hrt.52.1.82
//Blob b;
//b.load(fnbib);
// this will fix dos2unix and remove blank lines 
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
//} // mutate worked 

out.exit(nm);
return 0;
} // guessbmj
// https://discovery.ucl.ac.uk/id/eprint/1538561/1/Laycock_Thesis%20Joanne%20Laycock%20(No%20copyright).pdf

 static  IdxTy guessdiscoverac(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessdiscoverac";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// TODO bug in testHTML causes hang but can be avoided
// don't parse pdf, 
// https://core.ac.uk/download/pdf/82352096.pdf
// https://core.ac.uk/reader/82352096
const StrTy ufic="sed -e 's;\\(eprint/[^/]*\\).*;\\1;g' ";
 //static  IdxTy Mutate(const InpTy & in , OutTy & out, const StrTy & cmd ) 
IdxTy murc=Mutate(in,out,ufic);
//const StrTy ubase=MutateOnly(uin, ufic,out);
//const bool eq=(ubase==uin);
//MM_ERR(MMPR3(nm,eq,uin))
// let the handlers sort it out... 
//if (eq){ out.exit(nm); return 0; } 


out.exit(nm); return 0;
} // guessdiscoverac

#if 0 
// has meta but empty, can probably scrape doi but cool info in json within html
// https://core.ac.uk/reader/76958504?utm_source=linkout
 2041  testHTML --mjmflags 0 -mjm jsonh.html  | grep  "{" | tail -n 1 | sed -e 's/[^{]*//' | ./mjsonu.out -cmd "sax-json" 
 2042  testHTML --mjmflags 0 -mjm jsonh.html  | grep -B 3 -A 3   "{" | cut -c1-80
 2043  testHTML --mjmflags 0 -mjm jsonh.html  | more
 2044  history | grep http
 2045  testHTML --mjmflags 0 -mjm jsonh.html  | grep -B 3 -A 3   "{" | cut -c1-80

 testHTML --mjmflags 0 -mjm jsonh.html  | grep  "{" | tail -n 1 | sed -e 's/[^{]*//' | ./mjsonu.out -cmd "sax-json"  | grep doi
./mjsonu.h205 ;GLOBAL[1];[1];props[1];pageProps[13]:doi=10.1038/nprot.2016.036

#endif

 static  IdxTy guesscoreac(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesscoreac";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// TODO bug in testHTML causes hang but can be avoided
// don't parse pdf, 
// https://core.ac.uk/download/pdf/82352096.pdf
// https://core.ac.uk/reader/82352096
const StrTy ufic="sed -e 's/.*\\///g' | sed -e 's/[^0-9]//g' ";
const StrTy num=MutateOnly(uin, ufic,out);
const StrTy fixed_uin=StrTy("https://core.ac.uk/reader/")+num;
const bool eq=(uin==fixed_uin);
MM_ERR(MMPR3(eq,uin,fixed_uin))
if ( num=="" )  { MM_ERR(" skipping no number in uin ")  out.exit(nm); return 1; } 
const StrTy fnq=(eq)?fn:fntemp;
Grc grc=eq?in.getter().normalget(fnq,fixed_uin,16):
in.getter().normalget(fnq,fixed_uin,0);
// use the html parser, then get the json files
// wth does the xformer do???
Blob b;
Xrc xrc= in.xformer().html_to_parsed(b,fnq,0);
// TODO when run on a pdf file, this hangs maybe the testHTML is not terminating 
/*
This cmd , 
 testHTML --mjmflags 0 --mjm  /tmp/fileKdqK0U

Does put junk on the next prompt with no cr lol, 

/home/documents/latex/proj/paradox/junk$ 62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c62;c


*/

#if 0 
 testHTML --mjmflags 0 -mjm xxx | grep "{" | tail -n 1 | sed -e 's/[^{]*//' |  mjsonu.out -cmd "sax-json"  | grep doi 
./mjsonu.h205 ;GLOBAL[1];[1];props[1];pageProps[8];identifiers[1]:identifiers=doi:10.1016/S0014-5793(99)00131-3
./mjsonu.h205 ;GLOBAL[1];[1];props[1];pageProps[13]:doi=10.1016/S0014-5793(99)00131-3
marchywka@happy:/home/documents/latex/proj/paradox/junk$ testHTML --mjmflags 0 -mjm xxx | grep "{" | tail -n 1 | sed -e 's/[^{]*//' |  mjsonu.out -cmd "sax-json"  | grep doi |  sed -e 's/.*doi=//'
./mjsonu.h205 ;GLOBAL[1];[1];props[1];pageProps[8];identifiers[1]:identifiers=doi:10.1016/S0014-5793(99)00131-3
10.1016/S0014-5793(99)00131-3



#endif

const StrTy doi=
        MutateOnly(StrTy(b), 
 //StrTy("testHTML --mjmflags 0 -mjm \"")+fnq+"\"  | grep  \"{\" | tail -n 1 | sed -e 's/[^{]*//' | mjsonu.out -cmd \"sax-json\"  | grep doi | sed -e 's/.*doi=//'",out);
 StrTy("testHTML --mjmflags 0 -mjm \"")+fnq+"\"  | grep  \"{\" | tail -n 1 | sed -e 's/[^{]*//' | mjsonu.out -cmd \"sax-json\"  | grep \"identifiers=doi\" | head -n 1  | sed -e 's/.*doi://'",out);

MM_ERR(MMPR(doi))
IdxTy rcx=get_xref(out,in,doi,nm,0);

// two approaches- get headless and find doi OR use highwire on .citation-tools
//const StrTy nurl=MutateOnly(uin, "sed -e 's/\\([0-9]\\)[.?#;%].*/\\1/g' | awk '{print $0\".citation-tools\"}'" ,out );
//InpTy in2(in,nurl);
//OutTy out2(out,nurl);
//IdxTy rc=handlehighwire(in2,out2,xflags); // in.mom()->Guess(in2,out2,in.rflags());
//MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
//out.adopt(out2);
IdxTy nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rcx; } 
#if 0
// try the headlessthing
//const StrTy code=MutateOnly(uin, "sed -e 's/.*org\\///' | sed -e 's;/.*;;g'" ,out );
//const StrTy  url="http://cogprints.org/cgi/export/eprint/"+code+"/BibTeX/cogprints-eprint-"+code+".bib";
//if ((url!=uin) &&(url.length()!=0))
//{
Grc grc1=in.getter().headlessget(fn,uin,7);
// this renders nicely although it may be scrapable as-is./
Blob b;
Xrc xrc=in.xformer().html_to_rendered(b,fn,0);
StrTy s=StrTy(b);
const StrTy doi=MutateOnly(s, "grep doi\\.org | sed -e 's/.*org.10/10/' | head -n 1  " ,out );
MM_ERR(MMPR2(s,doi))
IdxTy rcx=get_xref(out,in,doi,nm+" ( doi )" ,0);
 nnow=out.found();
if ( nnow!=nstart) if (!all) { out.exit(nm); return rc; } 
StrTy href;
 IdxTy  hrcp=handlehtmlparse( href, nm+" ( parser )" ,  in ,  out,  fn,  xflags );

if ( href.length()==0)  { out.exit(nm); return rc; } 
IdxTy rcfm=FixMeta(href, in,out,nm+" ( parser )",uin);
#endif


// also has meta
// but it also seems to have highwire. 
// http://dx.doi.org/10.1136/hrt.52.1.82
//Blob b;
//b.load(fnbib);
// this will fix dos2unix and remove blank lines 
//IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);

//IdxTy rcr=Recurse(in,out,url);
//} // mutate worked 

out.exit(nm);
return 0;
} // guesscoreac


// should be able to get nm from stack
template <class Tm > 
static  IdxTy synthesize(const InpTy & in , OutTy & out , const StrTy & url, const StrTy & nm, const StrTy & fnbib, const StrTy & name, const StrTy & type, const Tm & kvp, const IdxTy flags=0)  
{
const bool skip_adding=Bit(flags,0);
const bool fix_keys=!Bit(flags,1);
BibEntry bib;
bib.name(name);
bib.type(type);
bib.add(kvp);
if (fix_keys||true) bib.fix_keys();

// this is also done at final assembly but only
// if this name is bad. TODO fix this as dates may need to be added etc. 
//if (bib.name()=="") bib.make_name();
if (bib.name()=="") bib.name("DUM");
//if (bib.name()=="") bib.name("default");
if (bib.name()=="") bib.name(bib.missing_name());
MM_ERR(MMPR2(fnbib,bib.format()))
const StrTy bibtext=(bib.format());
Blob b(bibtext);
b.save(fnbib);
if (skip_adding) return 0; 
MM_ERR(" deciding to add "<<MMPR3(fnbib,url,nm))
IdxTy rcf= out.good_enough(b,fnbib,in,url,nm,0);
return rcf;

} // synthesize

// https://en.wikipedia.org/wiki/Wikipedia:Citing_Wikipedia
#if 0
pem .
marchywka@happy:/home/documents/cpp/proj/toobib/junk$ testHTML -mjm -mjmflags=0 xxx  2> yyy | grep date| tail -n 1 | mjm zed 8 > j
marchywka@happy:/home/documents/cpp/proj/toobib/junk$ echo sax-json j | mjsonu.out 
mjm>sax-json j
./mjsonu.h205 ;GLOBAL[1];[1]:@context=https://schema.org
./mjsonu.h205 ;GLOBAL[1];[2]:@type=Article
./mjsonu.h205 ;GLOBAL[1];[3]:name=SGK1
./mjsonu.h205 ;GLOBAL[1];[4]:url=https://en.wikipedia.org/wiki/SGK1
./mjsonu.h205 ;GLOBAL[1];[5]:sameAs=http://www.wikidata.org/entity/Q18031469
./mjsonu.h205 ;GLOBAL[1];[6]:mainEntity=http://www.wikidata.org/entity/Q18031469
./mjsonu.h205 ;GLOBAL[1];[7];author[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[7];author[2]:name=Contributors to Wikimedia projects
./mjsonu.h205 ;GLOBAL[1];[8];publisher[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[8];publisher[2]:name=Wikimedia Foundation, Inc.
./mjsonu.h205 ;GLOBAL[1];[8];publisher[3];logo[1]:@type=ImageObject
./mjsonu.h205 ;GLOBAL[1];[8];publisher[3];logo[2]:url=https://www.wikimedia.org/static/images/wmf-hor-googpub.png
./mjsonu.h205 ;GLOBAL[1];[9]:datePublished=2013-09-24T19:01:18Z
./mjsonu.h205 ;GLOBAL[1];[10]:dateModified=2020-10-11T19:44:47Z
./mjsonu.h205 ;GLOBAL[1];[11]:image=https://upload.wikimedia.org/wikipedia/commons/e/e3/Protein_SGK1_PDB_2R5T.png
./mjsonu.h205 ;GLOBAL[1];[12]:headline=protein-coding gene in the species Homo sapiens


#endif

typedef typename BibEntry::kvp_type Kvp;
//typedef typename mjm_wovdb<Tr,StrTy> CrossLinks;
//typedef  mjm_wovdb<Tr,StrTy> CrossLinks;
typedef  std::map<IdxTy,StrTy> CrossLinks;

 static  IdxTy guessclinicaltrialsgov(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessclinicaltrialsgov";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
Grc grc=in.getter().normalget(fn,uin,16);
Blob b;
Xrc xrc=in.xformer().html_to_rendered(b,fn,0);
StrTy s=StrTy(b);
Ss ss(s);
// this needs a bash version for customization but hard code for now
Ragged r;
const bool rdebug=false;
// maybe make the colon the separater... 
r.sep(":");
r.load(ss,rdebug);
MM_ERR(MMPR2(r.size(),s.length()))
// now should have kvp's... 
Kvp m;
// want to extract from the text in case it is a local copy 
// needs to get link from clinical trials identifier... 
CrossLinks links;
ExtractKvp(m,links, r, 0);

const StrTy url=uin;
//const StrTy name="";
const StrTy type="article";
m["quality"]="usgov self synthesized";
m["real_type"]="clinicaltrial";
std::vector<StrTy> v;
v.push_back("Official_Title");
FirstOfKeys(m,"title",v,0);
v.clear();
v.push_back("First_Posted");
v.push_back("Last_Update_Posted");
FirstOfKeys(m,"date",v,0);
v.clear();
v.push_back("Please_refer_to_this_study_by_its_ClinicalTrials_gov_identifier__NCT_number_");
v.push_back("ClinicalTrials_gov_Identifier");
FirstOfKeys(m,"issue",v,0);

v.clear();
v.push_back("Please_refer_to_this_study_by_its_ClinicalTrials_gov_identifier__NCT_number__link");
v.push_back("ClinicalTrials_gov_Identifier_link");
v.push_back("ClinicalTrials_gov_Identifier_2_link");
FirstOfKeys(m,"url",v,0);
if (m.find("url")!=m.end())
{
const IdxTy lin=atoi(m["url"].c_str());
if ( links.find(lin)!=links.end()) m["url"]=links[lin];
else MM_ERR("missing link url "<<MMPR2(lin,m["url"]))
}

v.clear();
v.push_back("Contact");
v.push_back("Responsible_Party");
v.push_back("Principle_Investigator");
FirstOfKeys(m,"author",v,1);

//StrTy name="dummyname"; // make from trial year or author etc. 
StrTy name="xxx"; // make from trial year or author etc. 
auto ii=m.find("issue");
if (ii!=m.end()) name=(*ii).second;
// this needs to fix keys starting with a number....
std::map<char,IdxTy> bad; 
bad['0']=1; bad['1']=1; bad['2']=1; bad['3']=1; bad['4']=1; bad['5']=1;
bad['6']=1; bad['7']=1; bad['9']=1; bad['9']=1;

Kvp mnew;
MM_LOOP(ii,m)
{
StrTy  k=(*ii).first;
if (bad.find(k.c_str()[0])!=bad.end()) k=StrTy("X_")+k;
mnew[k]=(*ii).second;

} // ii 
m=mnew;

IdxTy rc=synthesize(in, out, url, nm, fnbib, name, type, m, 0);
out.exit(nm);
return rc;
} // guessclinicaltrialsgov
template <class Tv> 
static IdxTy FirstOfKeys(Kvp & m, const StrTy & ko, const Tv & v, const IdxTy flags)
{
const bool trim_to_comma=Bit(flags,0);
const StrTy* kp=&ko;
IdxTy i=0;
if ( m.find(*kp)!=m.end()) if ( m[*kp].length()) return i; 
MM_LOOP(ii,v)
{
++i;
kp=&(*ii);
if ( m.find(*kp)!=m.end()) if ( m[*kp].length()) 
{ 
StrTy newone=m[(*ii)];
if (trim_to_comma)
{
IdxTy len=newone.length();
char c[len+1];
IdxTy j=0;
for(;j<len; ++j)
{
char cp=newone.c_str()[j];
if (cp==',') break;
c[j]=cp; 
}
c[j]=0;
newone=StrTy(c);
} // trim_to_comma

m[ko]=newone; // (*ii); 

return i;  
} 
} // ii 

return i;
} // FirstOfKeys 

static IdxTy ExtractKvp(Kvp & m,CrossLinks & links, const Ragged & r, IdxTy eflags )
{
IdxTy refline=~0;
const IdxTy sz=r.size();
if (sz==0) return 0; 
IdxTy ix=sz;
for(IdxTy i=0; i<sz; ++i)
{
--ix;
const auto &  l=r[ix]; 
if ( l.size()!=1) continue;
if (l[0]=="References") { refline=ix; break; } 
} //i  < sz 
if ( refline==~0)
{
MM_ERR(" required lynx or similar references at end, ignoring ")
refline=sz-1;

} // refline 
for(IdxTy i=0; i<refline; ++i)
{
// remove lead/trail spaces and index refs
const auto & l=r[i];
const IdxTy len=l.size();
if ( len<2) continue;
// may have a colon at end meaning take v from next line,
// TODO this is making up bad keys
//
#if 0
/*
lynx -dump "https://www.clinicaltrials.gov/ct2/show/NCT04770740" | grep Plasma
    1. Plasma desmosine levels [ Time Frame: Day 1 until day 28 or until
       Plasma desmosine levels before and during vitamin K supplementation
    2. Plasma dp-ucMGP levels [ Time Frame: Day 1 until day 28 or until


@article{Margot_Visser_Phase_Double_Blind_Randomized,
1__Plasma_desmosine_levels = {Day 1 until day 28 or until discharge if this is earlier.},
1__Safety = {Day 1 until day 28},
2__Plasma_dp_ucMGP_levels = {Day 1 until day 28 or until discharge if this is earlier.},
2__Serum_PIVKA_II_levels = {Day 1 until day 28 or until discharge if this is earlier.},



*/
#endif
//
StrTy k=StrUtil::extraneous_ws(l[0]);
StrTy v=StrUtil::extraneous_ws(l[1]);
const StrTy keyref=StrUtil::first_brace(k);
k=StrUtil::brace_remove(k);
const StrTy oneref=StrUtil::first_brace(v);
v=StrUtil::brace_remove(v);
k=FixKey(StrUtil::extraneous_ws(k));
v=StrUtil::extraneous_ws(v);
// need to remove braces and assume inclusion is the link number


// need to make key parasable but not sure if that should
// be here or in bibtex entry. 

StrTy kn=k;
if (m.find(k)==m.end()) m[kn]=v;
else
{
IdxTy n=2;
while (true)
{
Ss ss;
ss<<k<<"_"<<n;
if (m.find(ss.str())==m.end()) { kn=ss.str(); break; } 
++n;
}
 //m[k]=m[k]+StrTy(" ")+v;
 m[kn]=v;
} // else

if ( oneref.length()) { m[kn+StrTy("_link")]=oneref; }
} // i <refline

// lines were parsed on colon, these should all be screwed up lol 
for(IdxTy i=(refline+1); i<sz; ++i)
{
const auto &  l=r[i]; 
//if ( l.size()!=1) continue;
//if (l[0]=="References") { refline=ix; break; } 
StrTy line;
MM_LOOP(ii,l) {if (line.length()) line=line+":";  line=line+(*ii); } 
const char * p=line.c_str();
const IdxTy no=atoi(p);
//IdxTy pc=0;
const char * ppc=p;
// there is a  leading space  
//if (*ppc) ++ppc;
while (*ppc) { if (*ppc!=' ') break; ++ppc; } 
while (*ppc) { if (*ppc==' ') break; ++ppc; } 
while (*ppc) { if (*ppc!=' ') break; ++ppc; } 
const StrTy link=StrTy(ppc);
links[no]=link;
} // i 

return 0; 
} // ExtractKvp

static StrTy FixKey(const StrTy & kin) 
{
const IdxTy slen=kin.length();
if (slen==0) return StrTy("_");
const char * p=kin.c_str();
ParseTable pt;
const IdxTy sz=slen+2;
char * c = new char[sz];
IdxTy i=~0;
IdxTy j=0;
bool waswhite=false;
while (true)
{
++i; // need a finally for continue lol 
const char cx=p[i];
if (cx==0 ) { c[j]=0; ++j; break; } 
//Iv v=pt.lut(cx);
const bool pr=pt.printable(cx);
const bool di=pt.digit(cx);
const bool an=pt.alphanum(cx);
const bool wh=pt.white(cx);
if (!pr) { continue; } 
if (wh&&!waswhite) { waswhite=true; c[j]='_'; ++j; continue; } 
if (wh) { continue; }
waswhite=false;
if (an) { c[j]=cx; ++j; continue; }
c[j]='_'; ++j;
} // true
StrTy k(c);
if (c[0]==0) return k=StrTy("_");
delete[] c;
return k;
}
 static  IdxTy guesswikipedia(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesswikipedia";
out.enter(nm);
const bool try_new_way=true;
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
Grc grc=in.getter().normalget(fn,uin,16);
StrTy pagetitle="";
// needs t use the parse_both thing
if (try_new_way)
{

const StrTy url=uin;
const StrTy name="";
const StrTy type="article";

BomTex bt;
Ragged j,h;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const bool only_json=!true;
const IdxTy fnty=256+8+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rcp=bt.parse_both(h,j,fn,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))

Kvp m,me;
Blob b;
HierKvp hkvp;
// me will contain all kinds of stuff for sorting later... 
IdxTy ptrc=GetPageTitle(pagetitle,in,out,h,me,hkvp,xflags);
MM_ERR(MMPR(pagetitle))
// they changed all the fields... 
MM_SZ_LOOP(i,j,szj)
{
const Ragged::Line & l =j[i];
const IdxTy len=l.size();
if (len<2) continue;
const bool eq=l[len-2]=="=";
const int  off=eq?1:0;
if (eq)
{
MM_ERR( " should have no eq wtf ")
Ss ss; MM_SZ_LOOP(j,l,szl) { ss<<"|"<<j<<"|"<<l[j]; } MM_ERR(ss.str())

}
StrTy  k=l[len-2];
StrTy  v=l[len-1];
/// "Title" confuses things

k=StrUtil::fancy_to_lower(k);



Ss ss; MM_SZ_LOOP(j,l,szl) { ss<<"|"<<j<<"|"<<l[j]; } MM_ERR(ss.str())
//if (k.c_str()[0]=='@') continue; // bad result... 
//  values are bad too but mayube escape later in assembly ....
const char * fuy=v.c_str();
while (true)
{
char c=(*fuy);
if (c==0) break;
if (c=='\\') break;
++fuy;
}
if (*fuy) continue;

/*
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|mainEntity|5|http://www.wikidata.org/entity/Q175621
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|author|5|365468|6|@type|7|Organization
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|author|5|365468|6|name|7|Contributors to Wikimedia projects
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|publisher|5|365469|6|@type|7|Organization
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|publisher|5|365469|6|name|7|Wikimedia Foundation, Inc.
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|publisher|5|365469|6|logo|7|365470|8|@type|9|ImageObject
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|publisher|5|365469|6|logo|7|365470|8|url|9|https://www.wikimedia.org/static/images/wmf-hor-googpub.png
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|datePublished|5|2002-05-29T18:47:46Z
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|dateModified|5|2022-07-02T21:54:00Z
./mjm_med2bib_funcs.h3926 |0|GLOBAL|1|365466|2||3|365467|4|image|5|https://upload.wikimedia.org/wikipedia/commons/1/15/Cholecalciferol2.svg


*/
MM_ERR(MMPR3(k,v,len))
if (len==6)
{
if (k.c_str()[0]=='@') continue; // bad result... 
const char * p=k.c_str();
if (p[0]=='w') if (p[1]=='g') k=StrTy(p+2); //   continue; 
StrTy & vx=m[k];
if (vx.length()==0) { vx=v; }
else if (vx!=v) {  vx=vx+" and "+v; } 
continue;
} // 5

if (len==8)
{
StrTy k2=l[len-4];
k2=StrUtil::fancy_to_lower(k2);
if (k=="@type") {k2=k2+"_type"; }
//else if (k!="name") continue;
else if (k!="name") k2=k2+"_"+k;
const char * p=k.c_str();
if (p[0]=='w') if (p[1]=='g') k=StrTy(p+2); //   continue; 
//if (p[0]=='w') if (p[1]=='g')  continue; 
StrTy & vx=m[k2];
if (vx.length()==0) { vx=v; }
// duplicates happen to be only anbd adjacent otherwise  ... 
else if (vx!=v) {  vx=vx+" and "+v; } 
continue;
} // 5



} // i
MM_ERR(MMPR2(nm,m.size()))
if (m.size()>2)
{
if (m["title"]=="") m["title"]=pagetitle;
m["quality"]="community";
m["journal"]="Wikipedia";
IdxTy rc=synthesize(in, out, url, nm+"(new)" , fnbib, name, type, m, xflags);
const IdxTy nhave=out.found();
MM_ERR(MMPR4(nm,nstart,nhave,all))
if (!all) if (nhave!=nstart) { out.exit(nm); return rc; } 
} // m ok 


} // try_new_way



const StrTy s1="testHTML --mjmflags 0 -mjm \""+fn+"\"";
const StrTy s2="grep date| tail -n 1 | mjm zed 8";
const StrTy s3="mjsonu.out -cmd \"sax-json\" "; 
const StrTy sc=s1+" | "+s2+" | "+s3;
const StrTy sauth = sc+" | grep -i author | grep -i name= | sed -e 's/.*name=//g' ";
const StrTy spub = sc+" | grep -i publisher | grep -i name= | sed -e 's/.*name=//g' ";
//const StrTy stitle = sc+" |  grep -i headline= | sed -e 's/.*headline=//g' ";
const StrTy stitle = sc+" |  grep \"name = \" | sed -e 's/.*name = //g' ";
const StrTy syear = sc+" |  grep -i datemodified= | sed -e 's/.*dateModified=//g' | cut -c1-4 ";
const StrTy smo = sc+" |  grep -i datemodified= | sed -e 's/.*dateModified=//g' | cut -c6-7 ";
const StrTy sday = sc+" |  grep -i datemodified= | sed -e 's/.*dateModified=//g' | cut -c9-10 ";
const StrTy sname = sc+" |  grep -i \":name=\" | sed -e 's/.*.name=//g' | head -n 1 | sed -e 's/[^0-9a-zA-Z]/_/g' ";


// use the html parser, then get the json files
// wth does the xformer do???
Blob b;
// wasteful but probably faster than writing temp file...
const StrTy auth=MutateOnly(StrTy(b),sauth,out);
const StrTy pub=MutateOnly(StrTy(b),spub,out);
const StrTy title=MutateOnly(StrTy(b),stitle,out);
const StrTy year=MutateOnly(StrTy(b),syear,out);
const StrTy mo=MutateOnly(StrTy(b),smo,out);
const StrTy day=MutateOnly(StrTy(b),sday,out);
const StrTy name=StrTy("Wikipedia_")+MutateOnly(StrTy(b),sname,out);
 
//Xrc xrc= in.xformer().html_to_parsed(b,fn,0);
//const StrTy doi=
//        MutateOnly(StrTy(b),
// StrTy("testHTML --mjmflags 0 -mjm \"")+fn+"\"  | grep  \"{\" | tail -n 1 | sed -e 's/[^{]*//' | mjsonu.out -cmd \"sax-json\"  | grep doi | sed -e 's/.*doi=//'",out);
//


Kvp m;
const StrTy url=uin;
//const StrTy name="";
const StrTy type="article";
m["quality"]="community";
m["author"]=auth;
m["publisher"]=pub;
if (pagetitle.length()>title.length()) m["title"]=pagetitle; else m["title"]=title;
m["year"]=year;
m["month"]=mo;
m["day"]=day;
m["url"]=uin;
m["journal"]="Wikipedia";

IdxTy rc=synthesize(in, out, url, nm+"( old )" , fnbib, name, type, m, xflags);
out.exit(nm);
return rc;
} // guesswikipedia

 static  IdxTy guesspubchem(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesspubchem";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
Grc grc=in.getter().normalget(fn,uin,16);
// needs t use the parse_both thing
BomTex bt;
Ragged j,h;
// 4 is for json
//const bool insure_neq=Bit(flags,0); // earlier default but bad pattern 
//const bool insure_eq=Bit(flags,1); // earlier default but bad pattern 
//const bool only_ldjson= Bit(flags,2); // earlier default but bad pattern 
//const bool only_good_json=Bit(flags,3); // earlier default but bad pattern 
//const bool debug=Bit(flags,8); // earlier default but bad pattern 
const bool only_json=true;
const IdxTy fnty=256+8+(only_json?4:0)+1; //  not include eq, deb  , only good json 
IdxTy rcp=bt.parse_both(h,j,fn,fnty);
MM_ERR(" loaded "<<MMPR4(nm,fn,h.size(),j.size()))


Kvp m,me;

Blob b;

HierKvp hkvp;
StrTy pagetitle="";
// me will contain all kinds of stuff for sorting later... 
IdxTy ptrc=GetPageTitle(pagetitle,in,out,h,me,hkvp,xflags);
// need to remove bar etc from pagetitle...
// should be in cleanup... 



// wasteful but probably faster than writing temp file...
//const StrTy auth=MutateOnly(StrTy(b),sauth,out);
const StrTy pub="National Library of Medicine (US), National Center for Biotechnology Information";
const StrTy title=pagetitle; // MutateOnly(StrTy(b),stitle,out);
//const StrTy year=MutateOnly(StrTy(b),syear,out);
//const StrTy mo=MutateOnly(StrTy(b),smo,out);
//const StrTy day=MutateOnly(StrTy(b),sday,out);
const StrTy name="xxx"; // StrTy("Wikipedia_")+MutateOnly(StrTy(b),sname,out);
 

const StrTy url=uin;
//const StrTy name="";
const StrTy type="article";
//m["quality"]="community";
m["author"]="PubChem";
m["publisher"]=pub;
m["city"]="Bethesda (Maryland USA)";
m["publisher-location"]="Bethesda (Maryland USA)";
m["title"]=title;
//m["year"]=year;
//m["month"]=mo;
//m["day"]=day;
m["url"]=uin;
m["journal"]="PubChem Compound Summary ";

IdxTy rc=synthesize(in, out, url, nm, fnbib, name, type, m, xflags);
out.exit(nm);
return rc;
} // guesspubchem


#if 0

marchywka@happy:/home/documents/cpp/proj/toobib$ cat zzz | grep null | grep context | sed -e 's/[^{]*//' | mjsonu.out -cmd "sax-json"


 0 html  1 head  2 script  3 type 4 text application/ld+json
 0 html  1 head  2 script  3 (null)       {"@context":"http://schema.org","@graph":[{"@type":"Person","address":{"@type":"PostalAddress","addressLocality":"Raleigh-Durham, North Carolina Area"},"alumniOf":[{"@type":"Organization","name":"Genentech","url":"https://www.linkedin.com/company/genentech?trk=ppro_cprof"},{"@type":"Organization","name":"Pfizer","url":"https://www.linkedin.com/company/pfizer?trk=ppro_cprof"},{"@type":"EducationalOrganization","name":"Villanova University","url":"https://www.linkedin.com/school/villanova-university/"}],"awards":[],"image":{"@type":"ImageObject","contentUrl":"https://media-exp1.licdn.com/dms/image/C4E03AQEesNIwmvEP5g/profile-displayphoto-shrink_200_200/0/1537279800226?e=1626307200&v=beta&t=jnJDQWveq3pR_Zyk75AmnHv0KZ0M0gb-nUR1v5zOnqY"},"jobTitle":["Principal and Founder","Ambassador","Scientific Advisor"],"name":"Nikki Schultek","sameAs":"https://www.linkedin.com/in/nikki-schultek-03634887","url":"https://www.linkedin.com/in/nikki-schultek-03634887","worksFor":[{"@type":"Organization","name":"Intracell Research Group","url":"https://www.linkedin.com/company/intracell-research-group?trk=ppro_cprof"},{"@type":"Organization","name":"BrainFit4Life"},{"@type":"Organization","name":"Atopic Medical, LLC"}]},{"@type":"WebPage","url":"https://www.linkedin.com/in/nikki-schultek-03634887","reviewedBy":{"@type":"Person","name":"Nikki Schultek"}}]}

./mjsonu.h205 ;GLOBAL[1];[1]:@context=http://schema.org
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[1]:@type=Person
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[2];address[1]:@type=PostalAddress
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[2];address[2]:addressLocality=Raleigh-Durham, North Carolina Area
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[1];[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[1];[2]:name=Genentech
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[1];[3]:url=https://www.linkedin.com/company/genentech?trk=ppro_cprof
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[2];[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[2];[2]:name=Pfizer
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[2];[3]:url=https://www.linkedin.com/company/pfizer?trk=ppro_cprof
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[3];[1]:@type=EducationalOrganization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[3];[2]:name=Villanova University
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[3];alumniOf[3];[3]:url=https://www.linkedin.com/school/villanova-university/
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[5];image[1]:@type=ImageObject
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[5];image[2]:contentUrl=https://media-exp1.licdn.com/dms/image/C4E03AQEesNIwmvEP5g/profile-displayphoto-shrink_200_200/0/1537279800226?e=1626307200&v=beta&t=jnJDQWveq3pR_Zyk75AmnHv0KZ0M0gb-nUR1v5zOnqY
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[6];jobTitle[1]:jobTitle=Principal and Founder
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[6];jobTitle[2]:jobTitle=Ambassador
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[6];jobTitle[3]:jobTitle=Scientific Advisor
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[7]:name=Nikki Schultek
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[8]:sameAs=https://www.linkedin.com/in/nikki-schultek-03634887
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[9]:url=https://www.linkedin.com/in/nikki-schultek-03634887
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[1];[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[1];[2]:name=Intracell Research Group
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[1];[3]:url=https://www.linkedin.com/company/intracell-research-group?trk=ppro_cprof
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[2];[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[2];[2]:name=BrainFit4Life
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[3];[1]:@type=Organization
./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[3];[2]:name=Atopic Medical, LLC
./mjsonu.h205 ;GLOBAL[1];[2];@graph[2];[1]:@type=WebPage
./mjsonu.h205 ;GLOBAL[1];[2];@graph[2];[2]:url=https://www.linkedin.com/in/nikki-schultek-03634887
./mjsonu.h205 ;GLOBAL[1];[2];@graph[2];[3];reviewedBy[1]:@type=Person
./mjsonu.h205 ;GLOBAL[1];[2];@graph[2];[3];reviewedBy[2]:name=Nikki Schultek

#endif
 static  IdxTy guesslinkedin(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesslinkedin";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
StrTy fntemp3=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// from cmd line, the profiles , adding 16 to flags, seem to be a problem  
// 2022-01-01 need to sign in... 
//Grc grc=in.getter().meget(fntemp,uin,1+8);
// setting bit 0 will force it to use the scriupt "me -bro" 
// otherwise do internally and there wa s a bug in quote.. 
// but it still needs "dump dom" or it resturns pdf of the screen
Grc grc=in.getter().meget(fn,uin,1);
//MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))

const IdxTy ldjrc=handleldjson2(in,out,1);
// missing json means blog post probably 


out.exit(nm);
return 0;

} // guesslinkedin
// profile code is fragile code, obsolete//// 
// need to salvage post citations however... 
 static  IdxTy oldguesslinkedin(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesslinkedin";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp",1); //
StrTy fntemp2=out.fn("temp",1); //
StrTy fntemp3=out.fn("temp",1); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
// get the dom 
// from cmd line, the profiles , adding 16 to flags, seem to be a problem  
Grc grc=in.getter().meget(fntemp,uin,1+8);
MM_ERR(" the stupd  things shouldbe in messed up  "<<MMPR(fntemp))
//marchywka@happy:/home/documents/cpp/proj/toobib$ cat zzz | grep null | grep context | sed -e 's/[^{]*//' | mjsonu.out -cmd "sax-json"
// application/ld+json
MM_ERR(MMPR4(fntemp,fntemp3,fntemp2,fnbib))
const StrTy s1="testHTML --mjmflags 0 -mjm \""+fntemp+"\"";
const StrTy s2="grep -A 4  \"application/ld+json\" | grep null | grep context  |sed -e 's/[^{]*//' ";
const StrTy s3="mjsonu.out -cmd \"sax-json\" "; 
//const StrTy sc=s1+" | "+s2+" | "+s3;
const StrTy sc=s2+" | "+s3;
// blank file names a problem... 
MM_ERR(" start ")
IdxTy rch=cmd_exec(fntemp3,StrTy(),s1,out,8);

MM_ERR(" partial "<< MMPR4(fntemp,fntemp3,fntemp2,fnbib))
MM_ERR(" crud   "<< MMPR(sc))
Blob b3;
b3.load(fntemp3);
const StrTy f3=b3;
// cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, OutTy & out, const IdxTy flags)
//const bool always=true; // Bit(flags,0);
//const bool file= Bit(flags,1);
//const bool notrailingcrlf= Bit(flags,2); // 4 
//const bool dfile= Bit(flags,3); /// 8
//IdxTy rcp=cmd_exec(fntemp2,StrTy(),sc,out,8); 
//IdxTy rcp=cmd_exec(fntemp2,f3,sc,out,8); 
//StrTy fick=MutateOnly(b3,sc,out); 
//StrTy fick=MutateOnly(b3,"echo f ",out); 
StrTy fck=MutateOnly(b3,sc,out); 
b3=fck;
b3.save(fntemp2);
Blob b;
b.load(fntemp2);
const StrTy bs=StrTy(b);
const StrTy datec=" date";
const StrTy date=MutateOnly(bs,datec,out);
const StrTy yearc=" date +%Y";
const StrTy year=MutateOnly(bs,yearc,out);
const StrTy concat=" awk 'BEGIN{x=\"\"; }{if (x!=\"\") x=x\",\"; x=x$0}END{print x;}'";
const StrTy calum=" grep  \"alumniOf.*name=\" | sed -e 's/.*name=//' |"+concat;
MM_ERR(" looking ... ")
const StrTy alume=MutateOnly(bs,calum,out);
MM_ERR(" return  ... ")
const StrTy pnamec=" grep  \"@graph[^a-z]*name=\" | sed -e 's/.*name=//' | sed -e 's/,/ /g' | sed -e 's/  */ /g'  "+concat;
const StrTy pname=MutateOnly(bs,pnamec,out);
//const StrTy enamec="  sed -e 's/[^a-zA-Z]/_/' |"+concat;
const StrTy enamec=" sed -e 's/\\/$//' |  sed -e 's/.*\\///g' | sed -e 's/[^0-9a-zA-Z]/_/g'   "; // +concat;
//const StrTy ename=StrTy("LIPRO_")+MutateOnly(pname,enamec,out);
// the entry needs to be unambiguous doh 
const StrTy ename=StrTy("LIPRO_")+MutateOnly(uin,enamec,out);

const StrTy locc=" grep  \"address.*addressLocality\" | sed -e 's/.*addressLocality=//' |"+concat;
const StrTy loc=MutateOnly(bs,locc,out);
const StrTy workc=" grep  -i  \"worksFor.*name=\" | sed -e 's/.*name=//' |"+concat;
const StrTy work=MutateOnly(bs,workc,out);
const StrTy titlec=" grep  -i  \"jobtitle.*jobtitle=\" | sed -e 's/.*Title=//' |"+concat;
const StrTy title=MutateOnly(bs,titlec,out);
IdxTy job=1;
StrTy jt;
StrTy sep="";
// ./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[6];jobTitle[3]:jobTitle=Scientific Advisor
// ./mjsonu.h205 ;GLOBAL[1];[2];@graph[1];[10];worksFor[1];[1]:@type=Organization
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


//typedef typename BibEntry::kvp_type Kvp;
Kvp m;
const StrTy url=uin;
//const StrTy name="";
const StrTy type="person";
m["quality"]="unverfiied";
m["alumni"]=alume;
m["name"]=pname;
m["author"]=pname;
m["location"]=loc;
m["address"]=loc;
m["jobs"]=work;
m["year"]=year;
m["date"]=date;
m["title"]=jt+" "+alume; // title;
m["jt"]=title;
//m["date"]=date;
//m["day"]=day;
m["url"]=uin;
//m["journal"]="Wikipedia";

//IdxTy rc=0;
IdxTy rc=synthesize(in, out, url, nm, fnbib, ename, type, m, xflags);
out.exit(nm);
return rc;
} // oldguesslinkedin

#if 0
marchywka@happy:/home/documents/cpp/proj/toobib$ me -bro "https://www.npr.org/sections/coronavirus-live-updates/2021/05/29/1001590855/vietnam-detects-new-highly-transmissible-coronavirus-variant" > npr.html

 0 html  1 head  2 script  3 type 4 text application/ld+json
 0 html  1 head  2 script  3 (null)    {"@type":"NewsArticle","publisher":{"@type":"Organization","name":"NPR","logo":{"@type":"ImageObject","url":"https:\/\/media.npr.org\/chrome\/npr-logo.jpg"}},"headline":"Vietnam Detects New Highly Transmissible Coronavirus Variant","mainEntityOfPage":{"@type":"WebPage","@id":"https:\/\/www.npr.org\/sections\/coronavirus-live-updates\/2021\/05\/29\/1001590855\/vietnam-detects-new-highly-transmissible-coronavirus-variant"},"datePublished":"2021-05-29T17:05:40-04:00","dateModified":"2021-05-29T17:05:40-04:00","author":{"@type":"Person","name":["Wynne Davis"]},"description":"Vietnam's health ministry announced the discovery of the new variant on Saturday that has characteristics of two other strains. The country is currently dealing with a recent spike in infections.","image":{"@type":"ImageObject","url":"https:\/\/media.npr.org\/assets\/img\/2021\/05\/29\/gettyimages-1232802621-61c88b4e62098fbb10d3c07a5f35903bb2fb74b1.jpg"},"@context":"http:\/\/schema.org"}
 0 html  1 head  2 script


#endif

 static  IdxTy guessnpr(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessnpr";
out.enter(nm);
const StrTy & uin=in.uin();
IdxTy rc=0; 
const IdxTy ldjrc=handleldjson(in,out,xflags);
out.exit(nm);
return rc;
} // guessnpr




 static  IdxTy guesswebofproc(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guesswebofproc";
out.enter(nm);
const StrTy & uin=in.uin();
IdxTy rc=0; 

out.exit(nm);
return rc;
}
// there is no easy way to relate one to the other yet, 
// https://etd.ohiolink.edu/apexprod/rws_olink/r/1501/10?p10_etd_subid=125009&clear=10
// https://etd.ohiolink.edu/apexprod/rws_etd/send_file/send?accession=osu1486561555617168&disposition=inline
 static  IdxTy guessohiolink(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
// not working ye.. 
const StrTy nm="guessohiolink";
//if (in.depth()>5) return 0;
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
//Grc grc=in.getter().normalget(fn,uin,16);
//const StrTy & uin=in.uin();
const StrTy nurl=uin;
if (nurl!=uin) 
{
Mutate(in,out,nurl); 
const IdxTy nnow=out.found();
if (!in.collect_all()) if (nnow!=nstart) { out.exit(nm); return 0; }  
// insert code here if possible to relate two tings. 


} // ii 
out.exit(nm); // don't exit if we are mutating... 
return 0; 
} // handlemutate



 static  IdxTy guessxxx(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessxxx";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtex"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
Grc grc=in.getter().normalget(fn,uin,16);
// use the html parser, then get the json files


out.exit(nm);
return 0;
} // guessxxx
// new format has a "bibtex" link lol
 static  IdxTy guessgooglebooks(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="guessgooblebooks";
out.enter(nm);
const StrTy & uin=in.uin();
const StrTy fn=out.fn(); //
const StrTy fntemp=out.fn("temp"); //
const StrTy fnbib=out.fn("bibtexx"); //
const IdxTy nstart=out.found();
const bool all=in.collect_all();
//Grc grc=in.getter().normalget(fn,uin,16);
Grc grc=in.getter().headlessget(fn,uin,1);
out.enter("new");
handlebibtexlink(in,out,xflags);
out.exit("new");
{const IdxTy nnow=out.found();
if (!all)  if( nnow!=nstart) { out.exit(nm); return 0; } 
}
// use the html parser, then get the json files
// from handlebomtex
//typedef mjm_bomtex_json<Tr> BomTex;
BomTex bt;
////////////////////////////////////////////
// x.load(li.words(),1); 
Ragged r;
const IdxTy fnty=1+8; // html and json parse  only good ones 
IdxTy rc=bt.load_ragged(r,fn,fnty);
MM_ERR(" loaded "<<MMPR(r.size()))
std::map<StrTy,StrTy> m;
IdxTy rca=0; // bt.assemble(m,r,0);

// picking wrong title

//21|GLOBAL|12||13|child_node|14||17|title|=|Books on Google Play
//643|GLOBAL|328||329|title|=|Tin-Based Antitumour Drugs
//644|GLOBAL|328||329|subtitle|=|
//766|GLOBAL|328||329|last_page|359|title|=|228



MM_LOOP(ii,r)
{
const auto & l=(*ii);
const IdxTy lsz=l.size();
if (lsz<5) continue;
// this MAY need moficiation of there is an "=" between k and v 
if (lsz>=7) if(l[lsz-7]=="page") continue;
if (lsz>=7) if(l[lsz-7]=="child_node") continue;
if (lsz>=11) if(l[lsz-11]=="additional_info") continue;
bool skip=false;
for(IdxTy i=0; i<(lsz-1); ++i)
{
const StrTy & x=l[i];
if (x=="last_page") { skip=true; break; } 

}
if (skip) continue;

const StrTy k=l[lsz-3];
const StrTy v=l[lsz-1];
if (m[k].length()) m[k]+=" and ";
m[k]+=v;

} // ii 
StrTy s;
if (m["author"]=="" ) m["author"]=m["attribution"];
if (m["date"]=="" )
{
const char * p=m["publication_date"].c_str();
const IdxTy len=strlen(p);
if (len==10)
{
char cp[11];
memcpy(cp,p,11);
cp[4]='/';
cp[7]='/';

cp[10]=0;
m["publication_date"]=StrTy(cp);
}

}  // date
IdxTy rcf=bt.format(s,m,0);
MM_ERR(MMPR4(nm,fn,rca,rcf)<< MMPR4(s,m.size(),r.size(),rc))
if ((r.size()<3)||(m.size()<3))
{
// move this to the scraper... 
MM_ERR(" r or m too small ")
out.exit(nm);
return 0;
}
//const StrTy fnbib=out.fn("bibtemp",0);
BibEntry be;
be.add(m);
be.type("article");
be.name("missing");
Blob b;
b=be.format();
MM_ERR(" fixing kv map "<<MMPR(StrTy(b)))
b.save(fnbib);
IdxTy rcfge= out.good_enough(b,fnbib,in,uin,nm,0);




out.exit(nm);
return 0;
} // guessxxx





//////////////////////////////////////////////////////////////
//
// New to add
/*

# https://www.spandidos-publications.com/downloadRisFile/10.3892/mmr.2018.8940
marchywka@happy:/home/documents/latex/proj/drafts/etcycle$ 
# https://www.spandidos-publications.com/10.3892/mmr.2018.8940


*/
/*
osti.gov did this automatically lol,
wget -O xxx.pdf -S -v "https://www.osti.gov/servlets/purl/4481349"
--2022-07-04 20:47:54--  https://www.osti.gov/servlets/purl/4481349
Resolving www.osti.gov (www.osti.gov)... 192.107.175.222
Connecting to www.osti.gov (www.osti.gov)|192.107.175.222|:443... connected.
HTTP request sent, awaiting response... 
  HTTP/1.1 302 
  Date: Tue, 05 Jul 2022 00:47:53 GMT
  Server: Apache
  Strict-Transport-Security: max-age=31536000; includeSubDomains; preload
  Location: https://www.osti.gov/biblio/4481349
  Content-Length: 0
  Set-Cookie: JSESSIONID=03C1AEB954E9760FA845DBBA78533D2B; Path=/; HttpOnly
  Set-Cookie: BIGipServerlbapp_tc3=3892314634.49300.0000; path=/; Httponly
  Keep-Alive: timeout=5, max=95
  Connection: Keep-Alive
  Set-Cookie: BIGipServerwww.osti.gov_pool=1132494278.20480.0000; path=/; Httponly; Secure
  X-Frame-Options: SAMEORIGIN
Location: https://www.osti.gov/biblio/4481349 [following]
--2022-07-04 20:47:54--  https://www.osti.gov/biblio/4481349
Reusing existing connection to www.osti.gov:443.
HTTP request sent, awaiting response... 
  HTTP/1.1 200 
  Date: Tue, 05 Jul 2022 00:47:53 GMT
  Server: Apache
  Strict-Transport-Security: max-age=31536000; includeSubDomains; preload
  Cache-Control: no-cache, no-store, must-revalidate
  Pragma: no-cache
  Expires: Thu, 01 Jan 1970 00:00:00 GMT
  Content-Type: text/html;charset=UTF-8
  Vary: Accept-Encoding
  Keep-Alive: timeout=5, max=94
  Connection: Keep-Alive
  Transfer-Encoding: chunked
  X-Frame-Options: SAMEORIGIN
Length: unspecified [text/html]
Saving to: ‘xxx.pdf’

xxx.pdf                      [  <=>                           ] 128.85K   486KB/s    in 0.3s    

2022-07-04 20:47:54 (486 KB/s) - ‘xxx.pdf’ saved [131946]



*/

///////////////////////////////////////////////////////////


// should be ni handlers lol wtf

 static  IdxTy handlemutate(const InpTy & in , OutTy & out , const IdxTy xflags=0)  
{ 
const StrTy nm="handlemutate";
if (in.depth()>5) return 0;
const StrTy & uin=in.uin();
if (!in.uin().http())
{
MM_ERR(" skipp  because uin is not http  "<<MMPR2(nm,uin))
return 0; 
}

out.enter(nm);

const IdxTy nstart=out.found();
//const StrTy & uin=in.uin();
typedef std::vector<StrTy> Mutations;
Mutations v;
// TODO problem with 200 response an 404 content doh 
// lower case seems important 
// https://www.aaai.org/Library/AAAI/2006/aaai06-167.php
//https://www.aaai.org/Papers/AAAI/2006/AAAI06-167.pdf
v.push_back(" sed -e 's/[.\\/]pdf$/.php/'");
v.push_back(" sed -e 's/[.\\/]pdf$//'");
v.push_back(" sed -e 's/[.\\/]epdf[?].*//'");
v.push_back(" sed -e 's/\\/pdf$/\\/abs/'"); 
// arxiv puts these in the middle doh 
v.push_back(" sed -e 's/\\/pdf\\//\\/abs\\//'"); 
v.push_back(" sed -e 's/\\/pdf$/\\/full/'"); 
v.push_back(" sed -e 's;/full/;/abs/;'"); 
//https://hal-cea.archives-ouvertes.fr/cea-01803833/document
// https://hal-cea.archives-ouvertes.fr/cea-01803833
// compiler esc? wth was this? 
//v.push_back(" sed -e 's;/\/document;;'"); 
v.push_back(" sed -e 's;/\\/document;;'"); 
// https://www.recercat.cat/bitstream/handle/2072/368596 wrong... 
//                          ^^^^^^^^^
// https://www.recercat.cat/handle/2072/368596
// https://www.recercat.cat/bitstream/handle/2072/368596/10.1038%40s41560-019-0404-4.pdf?sequence=1

v.push_back(" sed -e 's;/bitstream\\(/handle/[^/][^/]*/[^/][^/]*\\).*;\\1;'"); 

//wget -O xxx.pdf -S -v "https://www.osti.gov/servlets/purl/4481349"
//Location: https://www.osti.gov/biblio/4481349 [following]
v.push_back(" sed -e 's;/servlets/purl/;/biblio/;'"); 

// guessdiscoverac 
// https://discovery.ucl.ac.uk/id/eprint/1538561/1/Laycock_Thesis%20Joanne%20Laycock%20(No%20copyright).pdf
v.push_back("sed -e 's;\\(eprint/[^/]*\\).*;\\1;g' ");
// https://jscholarship.library.jhu.edu/handle/1774.2/59133
// https://jscholarship.library.jhu.edu/bitstream/handle/1774.2/59133/ZIEGLER-DISSERTATION-2018.pdf?sequence=1
v.push_back("sed -e 's;bitstream/\\(handle/[^/]*/[^/]*/[0-9]*\\).*;\\1;g' ");


// https://iris.univpm.it/retrieve/handle/11566/263460/77778/Tesi_Cirilli.pdf
// https://iris.univpm.it/retrieve/handle/11566/263460/
// http://hdl.handle.net/11566/263460 
//AddPair("guesshandle",&Myt::guesshandle,"/handle/");




//Mutate(in,out," sed -e 's/[.\\/]pdf$//'");
//if (out.found()==0) { Mutate(in,out," sed -e 's/\\/pdf$/\\/abs/'"); }
//if (out.found()==0) { Mutate(in,out," sed -e 's/\\/pdf$/\\/full/'"); }
MM_LOOP(ii,v)
{
MM_ERR(" trying "<<MMPR2((*ii),uin))
Mutate(in,out,(*ii)); 
const IdxTy nnow=out.found();
MM_ERR(" done trying "<<MMPR3((*ii),nnow,uin))
if (!in.collect_all()) if (nnow!=nstart) break; 
} // ii 
out.exit(nm); // don't exit if we are mutating... 
return 0; 
} // handlemutate

static StrTy MutateOnly(const StrTy & in, const StrTy & cmd,OutTy & out)
{
StrTy  nurl; // =out.nurl();
cmd_exec(nurl,in,cmd,out,4);
return nurl;
}
// in is a file name 
static StrTy ExtractOnly(const StrTy & in, const StrTy & cmd,OutTy & out)
{
StrTy  nurl; // =out.nurl();
cmd_exec(nurl,in,cmd,out,4+2);
return nurl;
}


 static  IdxTy Mutate(const InpTy & in , OutTy & out, const StrTy & cmd ) 
{
const StrTy & uin=in.uin();
StrTy  nurl; // =out.nurl();
cmd_exec(nurl,uin,cmd,out,4);
MM_ERR(" recurins "<<MMPR3(uin,nurl,nurl.length()))
//MM_MSG(" recurins "<<MMPR3(uin,nurl,nurl.length()))
if ((nurl!=uin) &&(nurl.length()!=0))
{
IdxTy rcr=Recurse(in,out,nurl);
}
else MM_ERR(" ignoing mutate "<<MMPR3(nurl,uin,cmd))
return 0;

}
// TODO someone should check the call stack for loops, put is in a map 
 static  IdxTy Recurse(const InpTy & in , OutTy & out, const StrTy & nurl ) 
{
//MM_MSG(" recursion "<<MMPR2(in.depth(),in.uin()))
if (in.depth()>3)
{
MM_ERR(" recursion too deep"<<MMPR(in.depth()))
return 0; // break;
}

InpTy in2(in,nurl);
OutTy out2(out,nurl);
IdxTy rc=in.mom()->Guess(in2,out2,in.rflags());
MM_ERR(" return from recurions adopting now ... "<<MMPR(nurl))
//MM_MSG(MMPR(out2.found()))
out.adopt(out2);
return 0;

} // Recurse  

#if 0
handlemutate()
{
nmxx="handlemutate"
nm="$nmxx"
uin="$1"
echo guess $nm  $uin | smsg
#uinx=`echo $uin| sed -e 's/\/pdf$//'`
uinx=`echo $uin| sed -e 's/[.\/]pdf$//'`
echo $nm  $uinx | smsg
fn="$temp1"
mutateduin "$uinx" "$nm" "$1"
if [ "$handled" == "" ]
then
# the other  code uses this global 
uin="$1"
nm="$nmxx"
uinx=`echo $uin| sed -e 's/\/pdf$/\/abs/'`
echo $nm  $uinx | smsg
mutateduin "$uinx" "$nm" "$1"
fi

} #handlemutate 

ZZ

mutateduin()
{
nm="$2(mutatedin) "
uinx="$1"
uinorig="$3"
if [ "$uinorig" == "" ]
then
uinorig="$uinx"
fi
fn="$temp1"
normalget "$fn" "$uinx" 
doi=`cat "$fn" | doisfrommeta `
if [ "$doi" == "" ]
then
doi=`cat "$fn" | doisfromtext `
fi

echo $nm trying $uinx  doi  $doi  | smsg
#getcrossref "$nm" "$doi" "$fn" "$uinx"
getcrossref "$nm" "$doi" "$fn" "$uinorig"

} # mutatedin

#endif



