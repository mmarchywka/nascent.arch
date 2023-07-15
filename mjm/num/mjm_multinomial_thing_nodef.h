
// now needs mjm_format_text.h

class multinomial_thing
{
typedef multinomial_thing Myt;
typedef nomial_term Term;
typedef std::vector<IdxTy> Terms; 
typedef std::vector<Coef> Coefs; 
typedef std::map<IdxTy , Coef> TermsMap;

typedef OsTy Os;
// for analysis 
typedef std::vector<IdxTy> Vloc;
typedef std::map<Coef, Vloc > Cmap;
typedef std::map<IdxTy, Vloc > Fmap;
typedef typename Term::NomialMap Tmap;

//typedef std::map<IdxTy, int > NomialVars;
typedef Tmap NomialVars;

public:
multinomial_thing( ): m_st(0),m_th(0) {}
multinomial_thing(Tokenizer * st, TermHeap * th ): m_st(st),m_th(th) 
{

}
typedef TermsMap terms_map;
void clear() { m_map.clear(); } 
const terms_map & terms() const { return m_map;}
const IdxTy fudd() const { return m_map.size(); } 
Coef evaluate(const NomialValues & v ) { return Evaluate(v); }
// partial evaluation generating a new Myt
Myt evaluated(const NomialValues & v ) { return Evaluated(v); }
void add_term(const IdxTy tidx, const Coef & c) {AddTerm( tidx, c); } 
bool remove_term(const IdxTy tidx) {return RemoveTerm( tidx); } 

Myt diff_int( const IdxTy var, const IdxTy flags) { return  DiffInt(  var, flags); } 

// return a polynomial that has been integrated between points a and b with evaluation
// needed to overcome lack of "function" vars and multiple angles changing the 
// arguments 
Myt trig_diff_int( const IdxTy varc, const IdxTy vars, const Coef & a, const Coef & b, const IdxTy flags) 
{ 
return  TrigDiffInt(  varc,vars, flags); 

} 

typedef mjm_tab_data_spec<Tr> OutSpec;
void output(Os & os, OutSpec & outs) const { Output(os,outs); } 

StrTy dump(const IdxTy flags) const  { return Dump(flags); } 

void count_refs(RefCounts & rc) {  CountRefs(rc); } 
Myt operator+(const Myt & that) const  { return Add(that,0); }
Myt operator-(const Myt & that) const  { return Add(that,1); }
Myt operator*(const Myt & that) const  { return Times(that,0); }
Myt operator*(const Coef & c) const  { return Times(c,0); }
Myt operator=(const Myt & that) { return Assign(that,0); }
Myt substitute(const IdxTy var, const Myt & s,const IdxTy flags) 
{ return Substitute(var,s,flags); }
template <class Ty> 
void add_polynomial(const IdxTy vidx, const Ty & v,const IdxTy flags)
{ AddPolynomial(vidx, v,flags); }
template<class Td, class Tc> void set(const Td & d, const Tc & c, IdxTy flags=0){ Set(d,c,flags); }
Myt ad_hoc(  SpecParams * sp ) { return AdHoc(sp ); } 
Myt ad_hoc_general(  SpecParams * sp ) { return AdHocGeneral(sp ); } 

private:
template<class Td, class Tc> void  Set(const Td & d, const Tc & c, IdxTy flags){
clear();
MM_SZ_LOOP(i,d,sz)
{
const auto & pm =d[i];
NomialVars nv;
MM_LOOP(jj,pm){ //   MM_ERR(" setting "<<(*jj).first<<" "<<(*jj).second);  
nv[(*m_st)((*jj).first)]=(*jj).second; } 
Term t(nv);
const Coef &  cnew=(c[i]);
//MM_ERR(MMPR(cnew))
//Simplify(cnew);
add_term((*m_th).find(t),cnew);

} //i 

 }
Myt AdHocGeneral(   SpecParams * sp )
{
Myt x(m_st,m_th); // return value 

MM_LOOP(ii,m_map) 
{
//Myt y(m_st,m_th); // created by operation on each term 
free_term_collection y;
y.clear();
const IdxTy tidx=(*ii).first;
const Coef & ci=(*ii).second;
//MM_ERR(" adhocgeneral  "<<MMPR2(ci,((*m_th)(tidx).dump(m_st,0))))
if (ci==Coef(0)) continue;
//Coef cnew;
//Term nt=(*m_th)(tidx).ad_hoc_general(y,sp);
//(*m_th)(tidx).ad_hoc_general(y,sp,m_th);
(*m_th)(tidx).ad_hoc_general(&y,sp);
MM_SZ_LOOP(i,y,szy)
{
const Coef & cnew= y.c(i);
if (cnew==0) continue;
const Term & t=y.t(i);
 x.add_term((*m_th).find(t),cnew*ci);
//if (cnew !=Coef(0)) x.add_term((*m_th).find(nt),cnew*ci);
//x=x+y*ci;
}
//if (cnew !=Coef(0)) x.add_term((*m_th).find(nt),cnew*ci);

} // ii 



return x;
}


Myt AdHoc(   SpecParams * sp )
{
Myt x(m_st,m_th);
MM_LOOP(ii,m_map) 
{
const IdxTy tidx=(*ii).first;
const Coef & ci=(*ii).second;
//MM_ERR(" adhoc "<<MMPR2(ci,((*m_th)(tidx).dump(m_st,0))))
if (ci==Coef(0)) continue;
Coef cnew;
Term nt=(*m_th)(tidx).ad_hoc(cnew,sp);
if (cnew !=Coef(0)) x.add_term((*m_th).find(nt),cnew*ci);
} // ii 

return x; 
}

template <class Ty> void AddPolynomial(const IdxTy vidx, const Ty & v,const IdxTy flags )
{

MM_SZ_LOOP(i,v,sz)
{
NomialVars nv;
nv[vidx]=i;
Term t(nv);
//const Coef cnew=Coef(v[i]);
Coef cnew=Coef(v[i]);
Simplify(cnew);
add_term((*m_th).find(t),cnew);
}

} // AddPolynomial

bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// this attempts to hold references while adding
// things to the heap which may cause errors so the heap
// needs to be locked or refs re-evaluated after a find
Myt Substitute(const IdxTy var, const Myt & that, const IdxTy flags) { 
// the new polynomial 
Myt x(m_st,m_th);
MM_LOOP(ii,m_map) 
{ 
const IdxTy tidx=(*ii).first;
//const Coef & ci=(*ii).second;
 Coef & ci=(*ii).second;
if (ci==Coef(0)) continue;
Simplify(ci);
//x.m_map[(*ii).first]+=(*ii).second; 
// if the term contains var to a negative power, this fails.
// if it is positive we need to have computed powers of the new polynomial...
// this reference is not valid after a term find unless
// we mark the heap as fixed, 
const Term & t=(*m_th)(tidx);
const int p=t.power(var);
if (p<0) 
{ MM_ERR(" negative power "<<MMPR2(((*m_st)(var)),p)<<" "<<t.dump(m_st,0)) } 
// we should not have a zero-power term but remove it anyway 
if (p==0)
{
// remove the dead power
const Coef cold=x.m_map[tidx];
x.m_map.erase(x.m_map.find(tidx));
//x.m_map[tidx]+=ci;
Term nt=t;
nt.zero_term(var);
x.add_term((*m_th).find(nt),cold+ci);
continue;
}
//MM_ERR(MMPR2(tidx, t.dump(m_st,0)))
// for positive powers, this generates a lot of terms... 
Myt xtemp(m_st,m_th);
xtemp=that;
// this should probably use a temp hash manager 
for(IdxTy i=1; i<p; ++i) {
//MM_ERR(MMPR4(i,p,xtemp.dump(0),that.dump(0)))
// this operation destroys existing references unless the 
// heap move was disabled... 
xtemp=xtemp*that;
 } 
//MM_ERR("t after power " << MMPR2(tidx,t.dump(m_st,0)))
MM_LOOP(jj,xtemp.m_map)
{
//const Coef cnew= ci*(*jj).second;
 Coef cnew= ci*(*jj).second;
if (cnew==Coef(0)) continue;
Simplify(cnew);
Term nt;
// the multiply above and the find below have invalidated
// prior t refs if the heap moved 
const Term & t=(*m_th)(tidx);
const Term & told=(*m_th)((*jj).first);
//MM_ERR(MMPR2((*jj).first,told.dump(m_st,0)))
//MM_ERR("t again " << MMPR2(tidx,t.dump(m_st,0)))
//t.expand(nt,var,(*m_th)((*jj).first));
t.expand(nt,var,told);
// this find invalidates the old refeences 
x.add_term((*m_th).find(nt),cnew);
} // jj
} // ii  

return x; 

}


Myt Assign(const Myt & x,const IdxTy flags )
{
m_st=x.m_st;
m_th=x.m_th;
m_map=x.m_map;
return *this;
}
Myt Add(const Myt & that, const IdxTy flags)const { 
const bool sub=Bit(flags,0);
Myt x(m_st,m_th);
MM_LOOP(ii,m_map) { x.m_map[(*ii).first]+=(*ii).second; } 
if (sub) 
{MM_LOOP(ii,that.m_map) { x.m_map[(*ii).first]-=(*ii).second; }  } 
else { MM_LOOP(ii,that.m_map) { x.m_map[(*ii).first]+=(*ii).second; }} 
return x; 
}
Myt Times(const Myt & that, const IdxTy flags)const  { 
//const bool sub=Bit(flags,0);
MM_ONCE(" changed powers to powes_one no test ",)
Myt x(m_st,m_th);
MM_LOOP(ii,m_map) { 
//x.m_map[(*ii).first]+=(*ii).second; 
const IdxTy tidx=(*ii).first;
const Coef & ci=(*ii).second;
if (ci==Coef(0)) continue;
MM_LOOP(jj,that.m_map)
{
const IdxTy tidxj=(*jj).first;
const Coef & cj=(*jj).second;
if (cj==Coef(0)) continue;
Term t;


//(*m_th)(tidx).powers(t,(*m_th)(tidxj),1,1);
(*m_th)(tidx).powers_one(t,(*m_th)(tidxj));


// this "find" invalidates the prior references if the heap
// is moved, however we do not appear to have any 
const IdxTy fadd=(*m_th).find(t);
//MM_ERR(" adding T "<<MMPR4(tidx,tidxj,fadd,t.dump(m_st,0)))
//x.add_term((*m_th).find(t),cj*ci);
x.add_term(fadd,cj*ci);
} // jj 
} // ii

//if (sub) 
//{MM_LOOP(ii,that.m_map) { x.m_map[(*ii).first]-=(*ii).second; }  } 
//else { MM_LOOP(ii,that.m_map) { x.m_map[(*ii).first]+=(*ii).second; }} 

return x; 
}
Myt Times(const Coef & c, const IdxTy flags)const  { 
//const bool sub=Bit(flags,0);
Myt x(m_st,m_th);
MM_LOOP(ii,m_map) { 
x.m_map[(*ii).first]=(*ii).second*c; 
} // ii
return x; 
}

Myt TrigDiffInt( const IdxTy varc, const IdxTy vars, const Coef & a, const Coef & b, const IdxTy flags) 
{
Myt x(m_st,m_th);
MM_LOOP(ii,m_map)
{
const Coef & ci=(*ii).second;
if (ci==Coef(0)) continue; 
Term t;
Coef cnew=((*m_th)((*ii).first)).trig_diff_int(t,varc,vars,a,b,flags);
Simplify(cnew);
if(cnew!=Coef(0)) x.add_term((*m_th).find(t),cnew*ci);
}


return x;
}

Myt  DiffInt( const IdxTy var, const IdxTy flags)
{
Myt x(m_st,m_th);
MM_LOOP(ii,m_map)
{
const Coef & ci=(*ii).second;
if (ci==Coef(0)) continue; 
Term t;
Coef cnew=((*m_th)((*ii).first)).diff_int(t,var,flags);
if(cnew!=Coef(0)) x.add_term((*m_th).find(t),cnew*ci);
}
/*
MM_SZ_LOOP(i,m_terms,sz)
{
const Coef & ci=m_coefs[i];
if (ci==Coef(0)) continue; 
// this returns a term but that is not stored yet. 
Term t;
Coef cnew=((*m_th)(m_terms[i])).diff_int(t,var,flags);
//if(cnew!=Coef(0)) x.add_term(cnew*ci,(*m_th)(m_terms[i]).evaluated(v));
if(cnew!=Coef(0)) x.add_term((*m_th).find(t),cnew*ci);
//MM_ERR(MMPR3(i,ci,cnew));
}
*/

return x;
} // DiffInt

void AddTerm(const IdxTy tidx, const Coef & c)
{
m_map[tidx]+=c;
//m_terms.push_back(tidx); m_coefs.push_back(c);
}

bool RemoveTerm(const IdxTy tidx)
{
auto ii=m_map.find(tidx);
if (ii==m_map.end()) return false; 
m_map.erase(ii);
}

void CountRefs(RefCounts & rc)  
{
//MM_LOOP(ii,m_terms) { ++ rc[(*ii)]; } 
MM_LOOP(ii,m_map) { ++ rc[(*ii).first]; } 
}

Myt Evaluated(const NomialValues & v )
{
Myt x(m_st,m_th);
//Coef c=Coef(0);
MM_LOOP(ii,m_map)
{
const Coef & ci=(*ii).second;
if (ci==Coef(0)) continue; 
Term t;
Coef cnew=((*m_th)((*ii).first)).evaluated(t,v);
if(cnew!=Coef(0)) x.add_term((*m_th).find(t),cnew*ci);
}
/*
MM_SZ_LOOP(i,m_terms,sz)
{
const Coef & ci=m_coefs[i];
if (ci==Coef(0)) continue; 
// this returns a term but that is not stored yet. 
Term t;
Coef cnew=((*m_th)(m_terms[i])).evaluated(t,v);
//if(cnew!=Coef(0)) x.add_term(cnew*ci,(*m_th)(m_terms[i]).evaluated(v));
if(cnew!=Coef(0)) x.add_term((*m_th).find(t),cnew*ci);
//MM_ERR(MMPR3(i,ci,cnew));
}
*/

return x;
}

Coef Evaluate(const NomialValues & v )
{
Coef c=Coef(0);
MM_LOOP(ii,m_map)
{
const Coef & ci=(*ii).second;
if (ci==Coef(0)) continue; 
c+=ci*(*m_th)((*ii).first).evaluate(v);

}
/*
MM_SZ_LOOP(i,m_terms,sz)
{
const Coef & ci=m_coefs[i];
if (ci==Coef(0)) continue; 
c+=ci*(*m_th)(m_terms[i]).evaluate(v);
//MM_ERR(MMPR3(i,ci,c));
}
*/
return c;
}

void MakeCmap(Tmap & vars, Cmap & cmap) const 
{

MM_LOOP(ii,m_map) {
//const Coef & c= (*ii).second;
const IdxTy tidx=(*ii).first;
(*m_th)(tidx).vars(vars); 
Coef  c= (*ii).second;
if (c<0) c=-c;
//if (c==0) {MM_ONCE(" zero coef included need to check ",) } 
if (c==0) {MM_ONCE(" zero coef excluded need to check ",)  continue;  } 
//const IdxTy tidx=(*ii).first;
cmap[c].push_back(tidx);

} // ii 

} // MakeCmap

void MakeFmap(Fmap & fmap, Ss & ss, Cmap & cmap) const 
{
MM_LOOP(ii,cmap)
{
Coef  ci= (*ii).first;
Factors(fmap,ss,ci,1,cmap);

} //ii 
}

StrTy ReduceAndFactor(const IdxTy flags) const
{
Ss ss;
ss<<"analysing";
Tmap vars;
Cmap  cmap;
Fmap fmap;
MakeCmap(vars,cmap);
const bool var_online=Bit(flags,8);
/*
MM_LOOP(ii,m_map) {
//const Coef & c= (*ii).second;
const IdxTy tidx=(*ii).first;
(*m_th)(tidx).vars(vars); 
Coef  c= (*ii).second;
if (c<0) c=-c;
//const IdxTy tidx=(*ii).first;
cmap[c].push_back(tidx);
} // ii 
*/
MakeFmap(fmap,ss,cmap);

ss<<" 2pi ";
const Coef cpi=Coef(2*M_PI);
Factors(fmap,ss,cpi,8,cmap);
ss<<" pi ";
const Coef cpi1=Coef(M_PI);
Factors(fmap,ss,cpi1,8,cmap);

ss<<CRLF<<"analyze_vars ";
MM_LOOP(ii,vars) { ss<<(*m_st)((*ii).first)<<" "; } 
ss<<CRLF;
const Coef  scale=(cmap.size()!=0)?(*(cmap.begin())).first:Coef(1);
ss<<" analyze_low_scale "<<MMPR2(scale,double(scale))<<CRLF;
MM_LOOP(ii,m_map) {
const IdxTy tidx=(*ii).first;
// not valid after heap update NB
const Term & t= (*m_th)(tidx); // .vars(vars); 
const Coef & v=(*ii).second; // <<" ";
//ss<<(*ii).second<<" ";
Coef vs=v/scale;
// this fails for doubles, 
//vs.simplify();
Simplify(vs);
ss<<"analyze_terms "<<v<<" "<<vs<<" "<<double(vs)<<" ";
if (var_online)
MM_LOOP(jj,vars) { ss<<((*m_st)((*jj).first))<<" "<<t.power((*jj).first)<<" "; }
else MM_LOOP(jj,vars) { ss<<t.power((*jj).first)<<" "; }
ss<<CRLF;

} // ii 
return ss.str();
} // ReduceAndFactor

template <class Ts > void Simplify(Ts & c) const
{
c.simplify();
}

void Simplify(double  & c) const  {}

template <class Tfuc > 
IdxTy Factors(Fmap & fmap, Ss & ss, const Coef & ci,  const IdxTy df, const Tfuc & cmap) const
{
MM_LOOP(jj,cmap)
{
Coef  cj= (*jj).first;
if ( cj<=ci) continue;
// TODO need an int operator for rational
//IdxTy  f= IdxTy(double(cj/ci));
//Coef rctor= ci*f;
//if (rctor==cj) ss<<MMPR3(ci,cj,f)<<CRLF;
const IdxTy f= Factor(cj,ci,df);
if (f!=0) ss<<MMPR3(ci,cj,Coef(f)/df)<<CRLF;

} // jj 
return 0;
}
// ASSFUDD SHJIASDFASDFSDFS
IdxTy Factor(const Coef & cj, const Coef & ci, const IdxTy df=1) const
{
// TODO need an int operator for rational
IdxTy  f= IdxTy(double(cj*df/ci));
Coef rctor= ci*f/df;
if (rctor==cj) return f; //  ss<<MMPR3(ci,cj,f)<<CRLF;
return 0;
} 
void OutputLUT(Os & os, OutSpec & outs) const { 
Ss ss;
Tmap vars;
Cmap  cmap;
Fmap fmap;
MakeCmap(vars,cmap);
const bool include_var_online=outs.var_online();
ss<<CRLF<<"analyze_vars "<<outs.prefix()<<" coef approx divpi scaled approxs";
MM_LOOP(ii,vars) { ss<<" "<<(*m_st)((*ii).first); } 
ss<<CRLF;
 Coef  scale=(cmap.size()!=0)?(*(cmap.begin())).first:Coef(1);
if (scale==Coef(0)) scale=Coef(1);
//ss<<" analyze_low_scale "<<MMPR2(scale,double(scale))<<CRLF;
MM_LOOP(ii,m_map) {
const IdxTy tidx=(*ii).first;
// not valid after heap update NB
const Term & t= (*m_th)(tidx); // .vars(vars); 
const Coef & v=(*ii).second; // <<" ";
if (v==Coef(0))
{
MM_ONCE(" again exlude zero but whi is it here ?",) continue ;  
}
//ss<<(*ii).second<<" ";
Coef vs=v/scale;
// this fails for doubles, 
//vs.simplify();
Simplify(vs);
ss<<std::setprecision(outs.prec())<<"analyze_terms "<<outs.prefix()<<" "<<v<<" "<<double(v)<<" "<<(2*double(v/(2*M_PI)))<<" "<<vs<<" "<<double(vs);
if ( include_var_online ) 
{MM_LOOP(jj,vars) { ss<<" "<<(*m_st)((*jj).first)<<" "<<t.power((*jj).first); } } 
else { MM_LOOP(jj,vars) { ss<<" "<<t.power((*jj).first); } } 
// should make sure nothing is in t that we missed... 

ss<<CRLF;

} // ii 

os<<ss.str();
} // OutputLUT


void Output(Os & os, OutSpec & outs) const { 
OutputLUT(os,outs);


//Output(o,outs);


 } 

StrTy Dump(const IdxTy flags) const 
{
const bool coef_d=Bit(flags,1);
const bool analyze=Bit(flags,16);
//MM_ERR(MMPR3(coef_d,analyze,flags))
Ss ss;
StrTy sep="+";
IdxTy j=0;
MM_LOOP(ii,m_map)
{
const Coef & c= (*ii).second;
if (c!=0) 
{
const IdxTy tidx=(*ii).first;
if (j!=0) ss<<sep;  
if(coef_d) ss<<double(c);
else ss<<c;
ss<<(*m_th)(tidx).dump(m_st,flags); ++j; } 

}
if (analyze) 
{
ss<<CRLF<<ReduceAndFactor(flags)<<CRLF; 

}
/*
MM_SZ_LOOP(i,m_terms,sz)
{
const Coef & c= m_coefs[i];
if (c!=0) {if (j!=0) ss<<sep;  ss<<c<<(*m_th)(m_terms[i]).dump(m_st,flags); ++j; } 

}
*/
return ss.str();
}

Tokenizer * m_st;
TermHeap * m_th;
//Terms m_terms;
//Coefs m_coefs;
TermsMap m_map;

}; // multinomial_thing
