#ifndef MJM_RULE_LIST_H__
#define MJM_RULE_LIST_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

#include "mjm_thread_util.h"
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

#define MJM_HAVE_REGEX
#ifdef MJM_HAVE_REGEX
#include <regex>
#endif



// Thu May 23 06:28:39 EDT 2019
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_rule_list   
// g++ -std=gnu++11 -DTEST_MJM_RULE_LIST -I. -I../../mjm/hlib -gdwarf-3 -O0  -x c++ mjm_rule_list.h  -lpthread

template <class Tr>
class mjm_rule_list 
{
 typedef mjm_rule_list Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef std::vector<StrTy> Line;
typedef std::map<StrTy,IdxTy> HitMissCache;
// typedef typename Tr::MyBlock  MyBlock;
// TODO use the tokenizer.... 

//class rule_value : public std::vector<StrTy>
class rule_value : public std::vector<Line>
{
public:
//template <class Ty> 
bool vset(StrTy & x) const  { return Vset(x); }
void add(const Line & l, const IdxTy n) { Add(l,n); } 
void add(const StrTy & s) { Add(s); } 
StrTy dump(const Line & l, const IdxTy flags=0) { return Dump(l,flags); } 

private:
void Add(const StrTy & s) {
Line l(1);
l[0]=s;
Add(l,0);
 } 
void Add(const Line & l, const IdxTy n)  
{
const IdxTy  sz=l.size();
if (n>=sz) { push_back(Line()); return; } 
Line x(sz-n);
for(IdxTy i=n; i<sz; ++i) x[i-n]=l[i];
push_back(x);

}


bool Vset(StrTy & x) const  { 
if (size()==0) return false; 
const Line & xx=back(); 
if (xx.size()==0) return false; 
x=xx.back();
return true; 

 }
StrTy  Dump(const Line & l, const IdxTy flags) { 
Ss ss;
MM_LOOP(ii,l) { ss<<" "<<(*ii); } 
return ss.str();

 } 

}; // rule_value


typedef rule_value RuleValue;

class rule
{
public:
rule() {}
//rule(const StrTy& exp, const Line & w, const IdxTy s): m_expression(exp),m_regex(m_expression){ Add(w,s); }
rule(const StrTy& exp, const Line & w, const IdxTy s): m_expression(exp){MakeRegex();  Add(w,s); }
//rule(const StrTy& exp,  const StrTy & w): m_expression(exp),m_regex(m_expression){ Add(w); }
rule(const StrTy& exp,  const StrTy & w): m_expression(exp){MakeRegex(); Add(w); }
bool matches(const StrTy & exp) { return Matches(exp); } 
const RuleValue & value() const { return m_value; } 
//const Line & value() const { return m_value.back(); } 
StrTy dump() { return Dump(); } 
private:
typedef std::map<StrTy,IdxTy> HitMissCache;
typedef std::regex Regex;
StrTy m_expression;
RuleValue m_value;
mutable HitMissCache m_cache;
Regex m_regex;
// this is great but it does not distinguish enries 
void Add( const Line & w, const IdxTy s)
{
m_value.add(w,s);
//const IdxTy sz=w.size();
//for(IdxTy i=s; i<sz; ++i) m_value.push_back(w[i]); 
//MM_MSG(" adding "<<MMPR2(s,sz))

}
void Add(  const StrTy & w  ) {
//MM_MSG(" adding "<<MMPR(w))

m_value.add(w);

//m_value.push_back(w);
  }
void MakeRegex()
{
try { 
if (m_expression=="*")
{
m_expression=".*";
MM_ERR(" m_expression changed from * to .* for "<<MMPR(m_expression)) 

}
m_regex=Regex(m_expression); 
} catch (std::regex_error&  e) // was catching by value.. 
{
MM_ERR( "regex error "<<MMPR(m_expression))

}

}
bool Matches(const StrTy & exp) 
{ 
// TODO FIXME incredibly slow should maintain mutable list of hits and misses 
//Regex rgx(m_expression);
auto ii=m_cache.find(exp);
if (ii!=m_cache.end()) return ((*ii).second==1);
const bool hit= ( std::regex_search(exp.c_str(),m_regex));
if (m_cache.size()>100) m_cache.clear();
m_cache[exp]=hit?1:0;
return hit;

}
StrTy Dump()
{
Ss ss;
ss<<m_expression;
//MM_LOOP(ii,m_value) {ss<<" "<<(*ii); } 
MM_LOOP(ii,m_value) {ss<<" "<<m_value.dump(*ii); } 
return ss.str();
} // Dump
 

}; // rule



typedef rule Rule;

class rule_list : public std::vector<Rule>
{

public:
void add(const StrTy& exp, const Line & w, const IdxTy s) { Add(exp,w,s); }
void add(const StrTy& exp, const StrTy& w) { Add(exp,w); }
Rule * first_match( const StrTy & exp,const IdxTy skip) { return FirstMatch(exp,skip); } 
StrTy dump() { return Dump(); } 
private:
void Add(const StrTy& exp, const Line & w, const IdxTy s) { this->push_back(Rule(exp,w,s)); } 
void Add(const StrTy& exp, const StrTy & w) { this->push_back(Rule(exp,w)); } 

// there is no point returning all of them even though it is
// slow to sort through re-matching should not be common enough 
// to change order. 
Rule * FirstMatch( const StrTy & exp, const IdxTy skip=0)  
{
IdxTy nfound=0; 
//MM_LOOP(ii,(*this)) { if ( (*ii).matches(exp)) return &(*ii); } // ii 
const IdxTy sz=(*this).size();
for(IdxTy i=0; i<sz; ++i) 
{ Rule & x=(*this)[sz-1-i];  if (x.matches(exp)){++nfound; if (nfound>skip)   return &x; }  } // ii 
return NULL;
} // first_match
StrTy Dump()
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).dump()<<CRLF;
} // ii 
return ss.str();
}// Dump 
}; // rule_list


typedef rule_list RuleList;
typedef std::map<StrTy, RuleList> RuleMap;
typedef typename RuleMap::iterator MapItor;

public:
mjm_rule_list():m_debug(0) {}
mjm_rule_list(const IdxTy deb):m_debug(deb) {}

template <class Ty>
bool apply(Ty & dest, const StrTy & rule_name, const StrTy & obj_name ) {return  Apply(dest,rule_name,obj_name); } 
void add( const Line & line  ) {Add(line); } 
template <class Ty, class Tf>  void add_lines(const Ty & ii, const Tf & ei) 
{Ty fck=ii;  while (fck!=ei) { add((*fck)); ++fck; } }
StrTy dump () { return Dump(); } 
private:
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};



mutable MutexVector m_mutex_vector;


void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }

void Init()
{
m_mutex_vector = MutexVector(MU_SZ);
}
Rule *  Applicable( const StrTy & rule_name, const StrTy & obj_name, const IdxTy skip=0 ) 
{
MapItor mi= m_map.find(rule_name);
if (mi==m_map.end()) return 0;
//MM_ERR(" found "<<MMPR(rule_name))
RuleList& rl = (*mi).second;
Rule * rvp =rl.first_match(obj_name,skip);
return rvp;
}

template <class Ty>
bool Apply(Ty & dest, const StrTy & rule_name, const StrTy & obj_name ) 
{ 
IdxTy skip=0; 
while (true) { 
Rule * rvp= Applicable(rule_name,obj_name,skip);
if (rvp==NULL) return false;
//MM_ERR(" found "<<MMPR((long int) rvp))
// return true even if it does not set it. 
StrTy ASUCK="";
bool success=rvp->value().vset(ASUCK); 
// this returns false if hit not found or a blank entry in latest hit. 
// should probably just return 
if ( success ){  Convert(dest,ASUCK); return true; }
++skip;
}
return !true; 
}
void  Convert(StrTy & y,  const StrTy & x) { y=x; }
void  Convert( D & y,  const StrTy & x) { y= atof(x.c_str()); }
void Convert( IdxTy & y, const StrTy & x) { y= myatoi(x.c_str()); }
void  Convert(int & y, const StrTy & x) { y=myatoi(x.c_str()); }
void  Convert(bool & y, const StrTy & x) { const char  p = x.c_str()[0]; y= (p=='1')||(p=='t')||(p=='T')||(p=='+'); }


template <class Ty>
//bool Apply(std::vector<StrTy>  & dest, const StrTy & rule_name, const StrTy & obj_name ) 
bool Apply(std::vector<Ty>  & dest, const StrTy & rule_name, const StrTy & obj_name ) 
{
IdxTy skip=0;
//MM_ERR(" apply "<<MMPR2(rule_name, obj_name)) 
Rule * rvp= Applicable(rule_name,obj_name,skip);
if (rvp==NULL) return false;
//MM_MSG(" found "<<MMPR((long int) rvp))
// return true even if it does not set it. 
StrTy ASUCK;
auto x=rvp->value();
// probably should append
//dest=x.back();
MM_LOOP(ii,x.back()){ // MM_ERR(" rule idx "<<(*ii))  
 dest.push_back((*ii));  } 
//MM_LOOP(ii,x){  dest.push_back((*ii));  } 
//.vset(ASUCK); 
//Convert(dest,ASUCK);

return true; 
}


//template <class Ty>
//bool Apply(std::vector<StrTy>  & dest, const StrTy & rule_name, const StrTy & obj_name ) 
bool Apply(std::vector<IdxTy>  & dest, const StrTy & rule_name, const StrTy & obj_name ) 
{
IdxTy skip=0;
MM_ERR(" apply "<<MMPR2(rule_name, obj_name)) 
Rule * rvp= Applicable(rule_name,obj_name,skip);
if (rvp==NULL) return false;
//MM_MSG(" found "<<MMPR((long int) rvp))
// return true even if it does not set it. 
//StrTy ASUCK;
auto x=rvp->value();
// probably should append
//dest=x.back();
MM_LOOP(ii,x.back()){ MM_ERR(" rule idx "<<(*ii))   dest.push_back(myatoi(*ii));  } 
//.vset(ASUCK); 
//Convert(dest,ASUCK);

return true; 
}





void Add( const Line & line  ) {
const IdxTy sz=line.size();
if (sz<3) return; // right now no blank values 
const StrTy & rule=line[0];
const StrTy & exp=line[1];
m_map[rule].add(exp,line,2);
 } 
StrTy Dump() 
{
Ss ss;
MM_LOOP(ii,m_map)
{
const StrTy k=(*ii).first;
ss<<" rule "<<MMPR(k)<<CRLF;
ss<<(*ii).second.dump();

} // ii 
return ss.str();
}
static int myatoi(const StrTy & s )  { return myatoi(s.c_str()); }
static int myatoi(const char * c)  { return ::strtol(c,0,0); }

RuleMap m_map;
IdxTy m_debug;

}; // mjm_rule_list
#ifdef  TEST_MJM_RULE_LIST
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::string  StrTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 

int main(int argc,char **args)
{
typedef mjm_rule_list<Tr>  Myt;
typedef Tr::StrTy StrTy;
//Myt x(argc,args);
Myt x;
//std::ifstream ifs(args[1]);
//std::istream & ifs=std::cin;
//while (!ifs.eof())
MM_MSG(" enter type rule exp val ")
MM_MSG(" type=add,apply  ")

char* ty = new char[100]; //  ru[100],ex[100],val[100];
char* ru = new char[100]; //  ru[100],ex[100],val[100];
char* ex = new char[100]; //  ru[100],ex[100],val[100];
char* val = new char[100]; //  ru[100],ex[100],val[100];
while (true) { 
scanf("%s %s %s %s",ty,ru,ex,val);
std::vector<StrTy> v;
v.push_back(ru);
v.push_back(ex);
v.push_back(val);
StrTy tyASUCK;
if( strcmp(ty,"add")==0) x.add( v  ); 
else x.apply( tyASUCK,  ru,  ex) ; 
MM_MSG(MMPR3(tyASUCK,ru,ex))
MM_MSG(x.dump())
}
//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_RULE_LIST_H__ 
