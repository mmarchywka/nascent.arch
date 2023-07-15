
class nomial_group
{
typedef nomial_group Myt;
protected:
//typedef multinomial_thing Multi;
typedef std::map<IdxTy, Nomials> CollMap;
typedef std::map<StrTy,Coef> NomialCharVars;
typedef std::map<IdxTy,Coef> NomialTokVars;
//typedef std::map<IdxTy, int > NomialVars;
typedef typename Term::NomialMap  NomialVars;

public:
//typedef std::map<StrTy, IdxTy> term_map;
typedef std::map<StrTy, int > term_map;
Nomials & multinomial(const StrTy & nm) { return Multinomial(nm.c_str()); } 
Nomials & multinomial(const char *  nm) { return Multinomial(nm); } 
void clear() { Clear(); } 
// only the isolate ones, leave the families and tables from
// derivers 
void clear_scratch() { m_map.clear(); } 
// take up a polynomial from different heap and string tokenizer
void adopt(const char * nm, const char * from, Myt & src)
{ Adopt(nm, from, src); }


template <class Ty> 
void add_polynomial(const char *nm, const char * var, const Ty & v,const IdxTy flags)
{ AddPolynomial(nm, var,  v,flags); }

void add_term(const char * nm, const term_map & tm, const Coef & c)
{ AddTerm(nm,  tm,c); }
bool remove_term(const char * nm, const term_map & tm)
{ return RemoveTerm(nm,  tm); }

void  add(const char * d, const char * s1, const char * s2 )
{ multinomial(d)=multinomial(s1)+multinomial(s2); } 
void  subtract(const char * d, const char * s1, const char * s2 )
{ multinomial(d)=multinomial(s1)-multinomial(s2); } 
void  multiply(const char * d, const char * s1, const char * s2 )
{ multinomial(d)=multinomial(s1)*multinomial(s2); } 
void  multiply(const char * d, const char * s1, const Coef & c )
{ multinomial(d)=multinomial(s1)*c; } 


// replace variable var in s with sv 
void  substitute(const char * d, const char * s, const char * var, const char * sv, const IdxTy flags )
{ Substitute(d, s,var, sv,flags ); } 

Coef evaluate(const char * nm, const NomialCharVars & vars)
{ return  Evaluate(nm,vars); }

void  evaluate_to(const char * d, const char * s, const NomialCharVars & v ) 
{ return EvaluateTo(d,s,v); }

void diff_int(const char *d , const char *s, const char *var, const IdxTy flags)
{
DiffInt(d , s, var,  flags);
}
// integrate a polynomial as if it were a trig polynomial
// and evaluate. Currently there is no way to handle functions and
// angles that change during an operation so this is the best I can do 
// although making variable into functions could make this work
// and return something with different arguments. 
void trig_diff_int(const char *d , const char *s, const char *varc, const char * vars, const Coef & a, const Coef & b, const IdxTy flags)
{ TrigDiffInt(d , s,varc, vars, a,  b,  flags); } 
//void TrigDiffInt(const char *d , const char *s, const char *varc, const char * vars, const Coef & a, const Coef & b, const IdxTy flags)

void  ad_hoc(const char * d, const char * s,    SpecParams * sp )
{  
//int var(const IdxTy i) const { return m_vars[i]; } 
//const StrTy  var_name(const IdxTy i) const { return m_var_names[i]; } 
//Myt & push_var(const IdxTy i ) { m_vars.push_back(i); return *this; } 
//Myt & push_var(const StrTy & i ) { m_var_names.push_back(i); return *this; } 
const IdxTy sz=sp->names_size();
for(IdxTy i=0; i<sz; ++i) { (*sp).push_var(m_st((*sp).var_name(i))); }
if (sp->create() ) { AdHocCreate(d,s,sp); return; } 
//MM_ERR(" normal ad hoc branck")
multinomial(d)=multinomial(s).ad_hoc(sp); 
} 

void  ad_hoc_general(const char * d, const char * s,    SpecParams * sp )
{  
const IdxTy sz=sp->names_size();
for(IdxTy i=0; i<sz; ++i) { (*sp).push_var(m_st((*sp).var_name(i))); }
//if (sp->create() ) { AdHocCreate(d,s,sp); return; } 
//MM_ERR(" normal ad hoc branck")
sp->st(&m_st);
multinomial(d)=multinomial(s).ad_hoc_general(sp); 
} 




StrTy dump(const IdxTy flags) { return Dump(flags); } 
StrTy dump(const char * nm, const IdxTy flags) { return Dump(nm,flags); } 
void count_refs(RefCounts & rc) { CountRefs( rc) ; } 
StrTy dump_refs(const IdxTy flags) { return  DumpRefs( flags) ; } 
private:


void Adopt(const char * nm, const char * from, Myt & src)
{

const Nomials & x= src.multinomial(from);
Nomials d;
const typename Nomials::terms_map & tm=x.terms();
// export each term in x and then import into the new one. 
MM_LOOP(ii,tm)
{
NomialCharVars cv;
const Term & tin = m_th.find((*ii).first);
//NomialTokVars tv;
//MM_LOOP(ii,vars) { tv[src.m_st((*ii).first)]=(*ii).second; } 
// the term does not have an m_st can't do it there. 
//MM_LOOP(jj,(*ii).fi) { tv[src.m_st((*ii).first)]=(*ii).second; } 

//NomialVars nv;
//MM_LOOP(ii,cv) { nv[m_st((*ii).first)]=(*ii).second; } 
Term t(tin,&m_st,&src.m_st);
IdxTy tidx=m_th.find(t);
d.add_term(tidx,(*ii).second); 
} // ii
multinomial(nm)=d;
} // Adopt



void  AdHocCreate(const char * d, const char * s,    SpecParams * sp )
{
//MM_ERR(" in adhoc create "<<MMPR2(d,s))
switch (sp->code())
{

case SpecParams::CREATE_MSIN_POLY: { AdHocMsin(d,  s, sp ); break; } 
//case SpecParams::CREATE_MSIN_POLY: { AdHocMsinAlt(d,  s, sp ); break; } 

default : MM_ERR(" bad create code "<<MMPR((sp->code())))

}  // switch 



} // AdHocCreate

// http://mathworld.wolfram.com/Multiple-AngleFormulas.html
// not working yet 
void  AdHocMsinAlt(const char * d, const char * s,    SpecParams * sp )
{
const IdxTy n = sp->count(0);
const IdxTy cidx=sp->var(0);
const IdxTy sidx=sp->var(1);
Coef bin=1;
NomialVars nv;
if (n==0)
{
Term t(nv);
IdxTy tidx=m_th.find(t);
multinomial(d).add_term(tidx,1); 
multinomial(s).add_term(tidx,0); 
return; 
}
nv[cidx]=n;
Term t(nv);
IdxTy tidx=m_th.find(t);
// the cos term 
multinomial(d).add_term(tidx,1<<(n-1)); 

Coef binc=1.0/Coef(n); 
const IdxTy n2=n>>1;
const IdxTy n2s=(n-1)>>1;
for(IdxTy k=0; k<=n2; ++k)
{

NomialVars nv;
nv[cidx]=n-(k<<1);
//  (n-k-1) 1
//  (   ) ---
//   (k-1)  k 
Coef cc=n*(((k&1)!=0)?(-1):1)*(1<<(n-(k<<1)-1))*binc;
Term t(nv);
IdxTy tidx=m_th.find(t);
multinomial(d).add_term(tidx,cc); 
Coef binc2=binc*(n-k-2)*(n-2*k)*Coef(n-2*k-1)/Coef(k+1);
if (k>n2s){ binc=binc2;  continue;}

//  (n-k-1) 
//  (   ) 
//    (k)  

Coef bins= binc*(n-2*k);

nv[cidx]=n-(k<<1)-1;
nv[sidx]=1;
Term t2(nv);
tidx=m_th.find(t2);
Coef sc=n*(((k&1)!=0)?(-1):1)*(1<<(n-(k<<1)-1))*bins;
multinomial(s).add_term(tidx,sc); 
binc=binc2;
}


} // AdHocsMsin

// these expressions are not ideal, see above 
// http://mathworld.wolfram.com/Multiple-AngleFormulas.html
void  AdHocMsin(const char * d, const char * s,    SpecParams * sp )
{
// doh these do not clear the polynomial first... 
 multinomial(s).clear();  
 multinomial(d).clear();  
//NomialVars nv;
//MM_LOOP(ii,tm) { nv[m_st((*ii).first)]=(*ii).second; } 
//Term t(nv);
//IdxTy tidx=m_th.find(t);
//multinomial(nm).add_term(tidx,c);
const IdxTy n = sp->count(0);
const IdxTy cidx=sp->var(0);
const IdxTy sidx=sp->var(1);
Coef bin=1;
//MM_ERR(MMPR3(n,cidx,sidx))
for(IdxTy k=0; k<=n; ++k)
{
NomialVars nv;
const int dir=n-k;
if (k!=0) nv[cidx]=k;
if ((dir)!=0) nv[sidx]=dir;
Term t(nv);
IdxTy tidx=m_th.find(t);
const int sw=dir&3;
// just to avoid counting properly for debug lol... 
//Coef sc=bin*sin((n-k)*.5*M_PI);
//Coef cc=bin*cos((n-k)*.5*M_PI);
// multinomial(d).add_term(tidx,cc); multinomial(s).add_term(tidx,sc);
switch (sw)
{
case 0 : { multinomial(d).add_term(tidx,bin); break; }
case 1 : { multinomial(s).add_term(tidx,bin); break; }
case 2 : { multinomial(d).add_term(tidx,-bin); break; }
case 3 : { multinomial(s).add_term(tidx,-bin); break; }

} // switch 

bin=bin*(dir)/(k+1);

}
if (n==0)  {
NomialVars nv;
Term t(nv); 
IdxTy tidx=m_th.find(t);
 multinomial(s).add_term(tidx,0);  }
} // AdHocMsin


// add a polynomial in variable v from a vector or similar object of coefficients
template <class Ty> void AddPolynomial(const char *nm, const char * var, const Ty & v,const IdxTy flags )
{
const IdxTy vidx=m_st(var);
multinomial(nm).add_polynomial(vidx,v,flags);

}


// replace variable var in s with sv 
void  Substitute(const char * d, const char * s, const char * var, const char * sv,const IdxTy flags )
{
const IdxTy vidx=m_st(var);
multinomial(d)=multinomial(s).substitute(vidx,multinomial(sv),flags);

}

void TrigDiffInt(const char *d , const char *s, const char *varc, const char * vars, const Coef & a, const Coef & b, const IdxTy flags)
{
multinomial(d)=multinomial(s).trig_diff_int(m_st(varc),m_st(vars),a,b,flags);
}
void DiffInt(const char *d , const char *s, const char *var, const IdxTy flags)
{
multinomial(d)=multinomial(s).diff_int(m_st(var),flags);
}
void  EvaluateTo(const char * d, const char * s, const NomialCharVars & vars ) 
{
NomialTokVars tv;
MM_LOOP(ii,vars) { tv[m_st((*ii).first)]=(*ii).second; } 
multinomial(d)=multinomial(s).evaluated(tv);
}

Coef Evaluate(const char * nm, const NomialCharVars & vars)
{
NomialTokVars tv;
MM_LOOP(ii,vars) { tv[m_st((*ii).first)]=(*ii).second; } 
return multinomial(nm).evaluate(tv);

}
void AddTerm(const char * nm, const term_map & tm, const Coef & c)
{
NomialVars nv;
MM_LOOP(ii,tm) { nv[m_st((*ii).first)]=(*ii).second; } 
Term t(nv);
IdxTy tidx=m_th.find(t);
multinomial(nm).add_term(tidx,c);
}
bool  RemoveTerm(const char * nm, const term_map & tm)
{
NomialVars nv;
MM_LOOP(ii,tm) { nv[m_st((*ii).first)]=(*ii).second; } 
Term t(nv);
IdxTy tidx=m_th.find(t);
return multinomial(nm).remove_term(tidx);
}



void CountRefs(RefCounts & rc)  
{
MM_LOOP(ii,m_map) {(*ii).second.count_refs(rc); } 
}
StrTy  DumpRefs(const IdxTy flags) {
Ss ss;
RefCounts rc;
CountRefs(rc);
MM_LOOP(ii,rc)
{
ss<<(*ii).first<<" "<<(m_th)((*ii).first).dump(&m_st,flags)
<<" "<<(*ii).second<<CRLF;
}
return ss.str();
} 

bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
StrTy Dump(const char * nm  , const IdxTy flags) {
Ss ss;
const auto ii=m_map.find(m_st(nm));
if (ii!=m_map.end())
//{ ss<<m_st((*ii).first)<<" "<<(*ii).second.dump(0)<<CRLF;  }
{ ss<<m_st((*ii).first)<<" "<<(*ii).second.dump(flags)<<CRLF;  }
else { ss<<nm<<" not found"<<CRLF; } 
return ss.str();
}
StrTy Dump(const IdxTy flags) {
Ss ss;
const bool dump_heap_index=!Bit(flags,0);
// return Dump(flags);
MM_LOOP(ii,m_map)
{
//ss<<m_st((*ii).first)<<" "<<(*ii).second.dump(0)<<CRLF; 
// needs to use flgs nice with others 
ss<<m_st((*ii).first)<<" "<<(*ii).second.dump(flags)<<CRLF; 

}
if (dump_heap_index)
{
ss<<CRLF<<" heap index "<<CRLF;
ss<<(m_th).dump(flags);
}

return ss.str();
 } 
Nomials & Multinomial(const char *  nm) { 
const IdxTy idx=m_st(nm);
if (m_map.find(idx)==m_map.end()) m_map[idx]= Nomials(&m_st, &m_th);
return m_map[idx]; 
} 
void Clear()
{
m_map.clear();
m_th= TermHeap();
m_st=Tokenizer();

}
protected:
//  the fudding shot compiler won't fudding let protected wotk
// FUDD SHOT 
//public:
// could get really clever and use different ones.. 
Tokenizer m_st;
TermHeap m_th;
CollMap m_map;

}; // nomial_group



class nomial_families : public nomial_group
{
typedef nomial_families Myt;
typedef nomial_group Super;
//typedef std::map<IdxTy, Nomials> CollMap;
typedef typename Super::CollMap CollMap;
typedef std::map<int , Nomials>  EntryMap;
typedef std::map<int , EntryMap> FirstMap;
typedef std::map<int , FirstMap> FamiliesMap;
typedef std::map<int , FamiliesMap> ExtFamiliesMap;

public:
nomial_families():Super() {}
void clear() { m_exfmap.clear();  m_fmap.clear(); Super::clear(); } 
// whty the ASSFUDD is this shot needed- it can't fuind tha base ass function?
// FUDD SHOT
Nomials & multinomial(const char *  nm) { 
return Super::multinomial(nm);
}
Nomials & multinomial(const char *  nm, const IdxTy m, const IdxTy n) { 
return  Multinomial( nm, m,  n)  ;} 
Nomials & multinomial(const char *  nm, const IdxTy m, const IdxTy n,const IdxTy p) { 
return  Multinomial( nm, m,  n,p)  ;} 


private: 
Nomials & Multinomial(const char *  nm, const int  m, const int  n,const int  p ) { 
const IdxTy idx=Super::m_st(nm);
FamiliesMap & efm = m_exfmap[idx];
FirstMap & fm = efm[m];
EntryMap & em= fm[n];
if ( em.find(p)==em.end()) em[n]=Nomials(&this->Super::m_st,&this->Super::m_th);
//if (m_fmap.find(idx)==m_fmap.end()) m_map[idx]= Nomials(&m_st, &m_th);
return em[n]; 


}
Nomials & Multinomial(const char *  nm, const int m, const int  n) { 
const IdxTy idx=Super::m_st(nm);
FirstMap & fm = m_fmap[idx];
EntryMap & em= fm[m];
auto assfudd =&this->Super::m_st; // &Super::m_st;
//if ( em.find(n)==em.end()) em[n]=Nomials(&Super::m_st,&Super::m_th);
//you NEVER used to have to derefernce this shot wthat that assfudd is
// the fudding point of hieratchy with this cumbersom shot ass notation
// FUCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
if ( em.find(n)==em.end()) em[n]=Nomials(assfudd,&this->Super::m_th);
//if (m_fmap.find(idx)==m_fmap.end()) m_map[idx]= Nomials(&m_st, &m_th);
return em[n]; 
} 


FamiliesMap m_fmap;
ExtFamiliesMap m_exfmap;


}; // nomial_families
