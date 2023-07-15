class nomial_term
{

typedef nomial_term Myt; 
typedef  mjm_numerical_expansions<Tr, Coef> NumEx; 
public:
// need to fix ctor for this 
// IIRC this was faster than unordered but need to make own ds here. 
// second try same result, try custom structure.
typedef std::map<IdxTy,int> NomialMap;
//typedef std::unordered_map<IdxTy,int> NomialMap;
// this really needs ONLY a map ctor but the stupid collections
// IIRC need defaults...and so does new [] wtf
nomial_term(): m_index_hash(0) {  } 
nomial_term(const NomialMap & nm): m_index_hash(0)  {  m_terms=nm; } 
nomial_term(const NomialMap & nm, const IdxTy h ): m_index_hash(h)  {  m_terms=nm; } 
// move from one tokenizer to another 
template <class St> nomial_term(const Myt & s, St * st, St * oldst)
{
// var, power
MM_LOOP(ii,s.m_terms)
{ m_terms[(*st)((*oldst)((*ii).first))]=(*ii).second; }
} // ctor


IdxTy hash() const { return m_index_hash; }
void zero_term( const IdxTy var ) 
 { ZeroTerm(*this,  var ); }  
// there was some reason not to do this- ??? 
//const NomialMap & vars() const { return m_map; } 

void vars(NomialMap & x) const { MM_LOOP(ii,m_terms) { x[(*ii).first]+=1; }} 
int power(const IdxTy var) const
{
const auto & f=m_terms.find(var);
if (f==m_terms.end()) return 0;
return (*f).second;
}
// this too is stupid as the term map should be set once
// and const forever after... 
void rehash() const { Rehash(); }  
// for now just let us calculate it 
//IdxTy hash(const IdxTy h )  { m_index_hash=h;  return m_index_hash; }
bool operator==(const Myt & that) const  { return Equals(that); } 
// need alt with powers stored... 
Coef evaluate(const NomialValues & v) const 
{ return Evaluate(v); } 
StrTy dump(Tokenizer * st, const IdxTy flags) const 
{ return Dump(st,flags); } 
// this should be the other way around, returning a 
// properly constructed const Myt  and updating a coef... 
Coef evaluated(Myt & n, const NomialValues & v) const
{ return Evaluated(n,v); } 
void powers(Myt & n, const Myt & that, const int pthis, const int pthat) const
{ Powers(n, that, pthis, pthat) ; }  
void powers_one(Myt & n, const Myt & that) const { PowersOne(n,  that); } 

// using a map of vars would eliminate temp terms of
// no interest... 
Coef diff_int(Myt & n, const IdxTy var, const IdxTy flags) const
{ return  DiffInt( n, var,  flags); } 
// this is at the term level now, the integration can create
// many terms but with differing arguments so need to evaluate as
// there is no function concept just vars
Coef trig_diff_int(Myt & n, const IdxTy varc, const IdxTy vars, const Coef & a, const Coef & b, const IdxTy flags)
{ return TrigDiffInt(n,varc, vars, a,  b,  flags); } 
// if possible, return single coef but now the single term 
// needs to be in a hokey array...  
//IdxTy ad_hoc(Myt ** n_p, Coef ** coefs_p,  const SpecParams & sp )
// use this sig when possible 
Myt ad_hoc( Coef & coef,   SpecParams * sp ) const
{return  AdHoc(  coef, sp ) ; } 

template <class Tnomial> IdxTy ad_hoc_general(Tnomial * resp, SpecParams * sp ) const
{ return AdHocGeneral(resp,sp); }

//t.expand(nt,var,(*jj).first);
void expand(Myt & n, const IdxTy var, const Myt & ti) const
{

Expand( n, var,  ti); //  const
}

private:
template <class Tnomial> IdxTy AdHocGeneral(Tnomial * resp, SpecParams * sp ) const
{
switch (sp->code())
{
case SpecParams::TRIG_POLY_ZERO_PI: { return ZeroToPIGeneral(resp,sp,1); } 
case SpecParams::GTRIG_POLY_ZERO_PI: { return ZeroToPIGeneral(resp,sp,1); } 
case SpecParams::TRIG_POLY_ZERO_2PI: { return ZeroToPIGeneral(resp,sp,2); } 
case SpecParams::GTRIG_POLY_ZERO_2PI: { return ZeroToPIGeneral(resp,sp,2); } 
default: MM_ERR(" bad code in AdHocGeneral "<<MMPR2(dump(0,0),((*sp).dump())))
} // switch 
return ~0; 

} // AdHocGeneral 


Myt AdHoc( Coef & coef,   SpecParams * sp ) const
{
switch (sp->code())
{

case SpecParams::TRIG_POLY_ZERO_PI: { return ZeroToPI(coef,sp,1); } 
case SpecParams::TRIG_POLY_ZERO_2PI: { return ZeroToPI(coef,sp,2); } 
default: MM_ERR(" bad code in AdHoc "<<MMPR2(dump(0,0),((*sp).dump())))
} // switch 
coef=Coef(1);
return *this; 

} //AdHoc

// resp needs to be able to add terms with new named variables so needs
// a string tokenizer etc. Not sure if it should add terms to heap
// had to see who has references etc that are not safe unless heap
// is locked leading to possible space exhaustion.  
template <class Tnomial> IdxTy ZeroToPIGeneral(Tnomial * resp, SpecParams * sp , const int npi  ) const
{
const int varc=sp->var(0);
const int vars=sp->var(1);
// this could be defered until it is determined which are zero 
Myt ndc=*this;
ZeroTerm(ndc,varc);
ZeroTerm(ndc,vars);
Myt nomega=ndc;

const int pc = VarPo(varc);
const int ps = VarPo(vars);
NumEx numex;
Coef norm=numex.cmsn_norm(pc,ps); 
const int lim=pc+ps;
// create a new term for the DC component adding a power
// of \pi
Coef cdc=numex.cmsn_w(0,pc,ps);
if (cdc!=0) { 
//if (npi==1) cdc=cdc*M_PI;
//else if (npi==2) cdc=cdc*(2.0*M_PI);
if (npi!=1) cdc=cdc*npi;
ndc.m_terms[(*(sp->st()))("\\pi")]=1;
cdc=cdc*norm;
resp->add(cdc,ndc);

 } // ==0

Coef c=Coef(0);
int omega=0;
if (npi==1)
{
while (omega<lim)
{
// need to get all the odd harmonics that have an
// odd piece left over ... 
++omega;
// remember this needs to integrate or divide by power ... 
Coef plus=numex.cmsn_w(omega,pc,ps);
Coef neg=numex.cmsn_w(-omega,pc,ps);
//MM_ERR(MMPR3(omega,plus,neg)<<MMPR2(pc,ps))
if ((neg*plus)<0) 
{
c+=Coef(Coef(-2)*plus/Coef(omega)); // normalization picks uip the sin denominator.
} // neg == -plus 
else if (neg!=plus) {
MM_ONCE(" plus minus not right "<<MMPR4(plus,neg,omega,cdc)<<MMPR3(pc,ps,npi),) 
//{MM_ERR(" plus minus not right "<<MMPR4(plus,neg,omega,cdc)<<MMPR3(pc,ps,npi))  }
}// - != + 
++omega; // omega has to be odd,,,, 
} // while
} // if npi==1
if (npi>2)
{ MM_ERR(" can't handle npi too big "<<MMPR(npi))}
//if (c!=0) c=c*norm;
// now have two terms, DC and omega differing in \pi factor
if( c!=0)  resp->add(c*norm,nomega);

return 0;
} // ZeroToPIGeneral

Myt ZeroToPI( Coef & c,  const SpecParams * sp, const IdxTy npi ) const 
{
Myt n=*this;
c=Coef(0);
const int varc=sp->var(0);
const int vars=sp->var(1);

ZeroTerm(n,varc);
ZeroTerm(n,vars);
const int pc = VarPo(varc);
const int ps = VarPo(vars);
NumEx numex;
// cos,sin
Coef norm=numex.cmsn_norm(pc,ps); 
const int lim=pc+ps;
//cmsn_w(const int twow, const IdxTy m, const IdxTy n )
Coef cdc=numex.cmsn_w(0,pc,ps);
// no longer supports rats... 
// just insert a new var for pi ... 
//n.m_terms[var]=0;
if (npi==1) c=cdc*M_PI;
//else if (npi==2) c=cdc*(M_PI+M_PI);
else if (npi==2) c=cdc*(2.0*M_PI);
int omega=0;
if (npi==1)
{
while (omega<lim)
{
// need to get all the odd harmonics that have an
// odd piece left over ... 
++omega;
// remember this needs to integrate or divide by power ... 
Coef plus=numex.cmsn_w(omega,pc,ps);
Coef neg=numex.cmsn_w(-omega,pc,ps);
//MM_ERR(MMPR3(omega,plus,neg)<<MMPR2(pc,ps))
//if (neg==(-plus)) 
if ((neg*plus)<0) 
{
c+=Coef(Coef(-2)*plus/Coef(omega)); // normalization picks uip the sin denominator.
} // neg == -plus 
else if (neg!=plus) {

MM_ONCE(" plus minus not right "<<MMPR4(plus,neg,omega,cdc)<<MMPR3(pc,ps,npi),) 

//{MM_ERR(" plus minus not right "<<MMPR4(plus,neg,omega,cdc)<<MMPR3(pc,ps,npi))  }
}// - != + 
++omega; // omega has to be odd,,,, 
} // while
} // if npi==1
if (npi>2)
{ MM_ERR(" can't handle npi too big "<<MMPR(npi))}

if (c!=0) c=c*norm;
//MM_ERR(MMPR3(c,pc,ps))
return n;

} // ZeroToPI


const bool Bit(const IdxTy flags, const IdxTy b) const
{ return ((flags&(1<<b))!=0); } 
void Expand(Myt & n, const IdxTy var, const Myt & ti) const
{
//MM_ERR("xxxadf")
n.m_terms=m_terms;
//MM_ERR("fooxxxadf")
n.m_terms[var]=0;
//MM_ERR("adsfasxxxadf")
n.m_terms.erase(n.m_terms.find(var));
//MM_ERR("adcasdcxxxadf")
Power(n.m_terms,ti.m_terms,1);
//MM_ERR("adf")
//void Power(NomialMap & nn, const NomialMap & nx, const int p) const

}
void ZeroTerm(Myt & n, const IdxTy var ) const
{
n.m_terms[var]=0;
n.m_terms.erase(n.m_terms.find(var));
}

const int VarPo(const IdxTy var) const
{
const auto ii=m_terms.find(var);
const bool have= (ii!=m_terms.end());
const int px=have?(*ii).second:0;
return px;
}
Coef TrigDiffInt(Myt & n, const IdxTy varc, const IdxTy vars, const Coef & a, const Coef & b, const IdxTy flags)
{
const int pc = VarPo(varc);
const int ps = VarPo(vars);
n.m_terms=m_terms;
n.m_terms[varc]=0;
n.m_terms[vars]=0;
n.m_terms.erase(n.m_terms.find(varc));
n.m_terms.erase(n.m_terms.find(vars));
Coef c=0;

return c;
}

Coef DiffInt(Myt & n, const IdxTy var, const IdxTy flags) const
{
//const auto ii=m_terms.find(var);
//const bool have= (ii!=m_terms.end());
//const int px=have?(*ii).second:0;
const int px = VarPo(var);
//const bool diff=Bit(flags,1);
const IdxTy br=flags&255;
switch (br)
{
case 0:  // differentiate 
{
if (px==0) return Coef(0);
n.m_terms=m_terms;
n.m_terms[var]=px-1;
return Coef(px); 
}
case 1:  // integrate 
{
//if (px==0) return Coef(0);
if (px==(-1)) MM_ERR(" unable to integrate -1 "<<MMPR(flags)<<" "<<Dump(0,0))
n.m_terms=m_terms;
n.m_terms[var]=px+1;
// should return nan or something similar... 
return Coef(1)/Coef(px+1); 
}
default : MM_ERR(" unrecognized op "<<MMPR(flags)<<" "<<Dump(0,0))
}; // switch 

return Coef(0); 
} // DiffInt

void Powers(Myt & n, const Myt & that, const int pthis, const int pthat) const
{
Power(n,*this,pthis);
Power(n,that,pthat);

}
void PowersOne(Myt & _n, const Myt & that) const
{
NomialMap & n=_n.m_terms;
MM_LOOP(ii,m_terms) { n[(*ii).first]+=(*ii).second; } 
MM_LOOP(ii,that.m_terms) { n[(*ii).first]+=(*ii).second; } 


}
void Power(Myt & n, const Myt & x, const int p) const
{
Power(n.m_terms,x.m_terms,p);
//MM_LOOP(ii,x.m_terms) { 
//const IdxTy i=(*ii).first; 
//const int px=(*ii).second; 
//n.m_terms[i]+=p*px;
//}
}

//NomialMap m_terms;
void Power(NomialMap & nn, const NomialMap & nx, const int p) const
{
MM_LOOP(ii,nx) { 
const IdxTy i=(*ii).first; 
const int px=(*ii).second; 
nn[i]+=p*px;
}
}

Coef Evaluated(Myt & n, const NomialValues & v) const
{
//MM_LOOP(ii,m_terms) { MM_ERR("terms "<< MMPR2((*ii).first,(*ii).second)) }
//MM_LOOP(ii,v) { MM_ERR("vals "<< MMPR2((*ii).first,(*ii).second)) }
Coef c=Coef(1);
MM_LOOP(ii,m_terms)
{
const IdxTy var=(*ii).first;
int poww=(*ii).second;
if (poww==0) continue;
const auto jj=v.find(var);
if (jj==v.end()) {n.m_terms[(*ii).first]=(*ii).second;  continue; } 
const Coef & val=(*jj).second;
PowEval( c,  val, poww);
} // ii 
return c;
} // Evaluated

Coef Evaluate(const NomialValues & v) const 
{
Coef c=Coef(1);
// a constant is ok doh 
//if (m_terms.size()==0) return 0; 
//MM_LOOP(ii,m_terms) { MM_ERR("terms "<< MMPR2((*ii).first,(*ii).second)) }
//MM_LOOP(ii,v) { MM_ERR("vals "<< MMPR2((*ii).first,(*ii).second)) }
MM_LOOP(ii,m_terms)
{
const IdxTy var=(*ii).first;
// should zero be pruned? 
int poww=(*ii).second;
//MM_ERR(MMPR2(var,poww))
if (poww==0) continue;
const auto jj=v.find(var);
if (jj==v.end()) return Coef(0);
const Coef & val=(*jj).second;
PowEval(c,val,poww);
if (c==Coef(0)) return c; // FIXME added recently... 
}
//MM_ERR(MMPR(c))
return c;
}
void PowEval(Coef & c, const Coef & val, int poww) const
{
// arrghhhh
if (poww>0) { while (poww>0) { c=c*val; --poww; } }
else {
const Coef vi=Coef(1)/val;
while (poww<0) { c=c*vi; ++poww; }
//MM_ERR(MMPR(c))
//while (poww<0) { c=c/val; --poww; }
}



} // PowEval


StrTy Dump(Tokenizer * st, const IdxTy flags) const 
{
Ss ss;
if (st!=0) 
{ MM_LOOP(ii,m_terms) { ss<<(*st)((*ii).first)<<"^"<<(*ii).second; } } 
else
{
{ MM_LOOP(ii,m_terms) { ss<<" .no st. "<<((*ii).first)<<"^"<<(*ii).second; } } 

}
return ss.str();
}
bool Equals(const Myt & that) const 
{
return ( m_terms==that.m_terms);
}
// this should never be needed 
void Rehash() const
{
m_index_hash=0;
MM_LOOP(ii,m_terms) { m_index_hash+=(*ii).second<<((*ii).first&15); }  
}
// this really needs to be const.... 
NomialMap m_terms;
// doh 
mutable IdxTy m_index_hash;
}; // nomial_term

// normally terms exist in a term heap manager for fast usage and managed
// storagei along with restrictions on access.  These terms have
// been tokenizer byuser provided string_tokenizer but not added
// to a term heap manager yet. These may result from manipulations
// of a single term creating more terms that need to added to a 
// multinomial 
// force access through the interface so it can be fixed later lol
class free_term_collection
{
typedef nomial_term Term;
typedef free_term_collection Myt;
typedef std::vector<Coef> C;
typedef std::vector<Term> T;
public:
void size(const IdxTy n ) { m_terms.resize(n); m_coefs.resize(n); }
void clear( ) { m_terms.clear(); m_coefs.clear(); }
IdxTy size() const { return m_terms.size(); } 
void add(const Coef & c,const Term & t) 
{ 
	m_coefs.push_back(c);
	m_terms.push_back(t); 
} 
const Coef & c( const IdxTy i ) { return m_coefs[i]; } 
const Term & t( const IdxTy i ) { return m_terms[i]; } 
private:
T m_terms;
C m_coefs;
}; // free_term_collection

