#ifndef MJM_SPARSE_2NOMIAL_H__
#define MJM_SPARSE_2NOMIAL_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
//#include "mjm_rational.h"
//#include "mjm_generic_iterators.h"
//#include "mjm_block_matrix.h"
//#include "mjm_instruments.h"

/*
2017-11-25
#ifdef TEST_SPARSE_2NOMIAL
 2408  g++ -DTEST_SPARSE_2NOMIAL -Wall -Wno-unused-function  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -I.. -gdwarf-3 -x c++ mjm_sparse_2nomial.h 
*/





// a polynomial of x and y up to max degreee n-1. 
// use block matrix for later extension to multinomial 
template <class Tv> class sparse_2nomial 
{
typedef sparse_2nomial<Tv> Myt;
class Tr{
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
//typedef LocTy location_type;
//typedef mjm_rational RatTy;
}; // Tr

typedef typename Tr::D D;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::StrTy StrTy;

typedef Tv value_type;
//typedef mjm_block_matrix<D> MyBlock;
//typedef mjm_block_matrix<IdxTy> MyBlockInt;
//typedef mjm_block_matrix<std::vector<D> > MyBlockVec;
//typedef mjm_generic_itor<2> LoopItor2;
//typedef mjm_block_matrix<Tv> Super;
//typedef Tv value_type;

template <class Tvel, class Ti, IdxTy n > class Element
{
typedef Element Myt;
public:
typedef std::vector<Ti> TermIdx;
typedef Tvel TermVal;
enum {INDEX_SIZE=n};
// really just need to make it invalid, that is all 
Element(): m_valid(false),m_v(0),m_idx(INDEX_SIZE) {}
// TODO the idx size could be violated here 
Element(const TermIdx & idx, const TermVal & v): m_valid(false),m_v(v),m_idx(idx) {Valid(); }
Element(const IdxTy & i, const TermVal & v): m_valid(false),m_v(v),m_idx(Idx(i)) {Valid(); }
Element(const IdxTy & i, const IdxTy & j, const TermVal & v): m_valid(false),m_v(v),m_idx(Idx(i,j)) {Valid(); }
const TermVal & value() const { return m_v; } 
const TermIdx & index() const { return m_idx; } 
const IdxTy & index(const IdxTy & i ) const { return m_idx[i]; } 
const bool valid() const { return m_valid; } 

bool operator==(const Myt & x ) const { return SameTerm(x)&&SameV(x); } 
bool operator!=(const Myt & x ) const { return !(SameTerm(x)&&SameV(x)); } 
Myt operator*(const Myt & that) const { Myt x=(*this); x.AddIdx(that); x.m_v=m_v*that.m_v; x.Valid(); return x; } 
Myt &operator*=(const Myt & that) { AddIdx(that); m_v=m_v*that.m_v; Valid(); return (*this); } 
Myt operator+(const Myt & that) const { if (!SameTerm(that)) return Myt(); Myt x=(*this);  x.m_v=m_v+that.m_v; x.Valid(); return x; } 
Myt operator-(const Myt & that) const { if (!SameTerm(that)) return Myt(); Myt x=(*this);  x.m_v=m_v-that.m_v; x.Valid(); return x; } 
Myt& operator+=(const Myt & that)  { if (SameTerm(that)) { m_v=m_v+that.m_v; Valid();} else {InValid(); }  return (*this); } 
Myt& operator-=(const Myt & that)  { if (SameTerm(that)) { m_v=m_v-that.m_v; Valid();} else {InValid(); }  return (*this); } 

private:
Myt & Differentiate(const IdxTy&  i)
{
Ti & p=m_idx[i];
if (p==0) InValid();
else { m_v=m_v*p; --p; } 
return *this; 
}
const TermIdx  Idx(const IdxTy&  i) const { TermIdx idx= TermIdx(n); idx[0]=i;  return idx; } 
const TermIdx  Idx(const IdxTy&  i, const IdxTy&  j) const { TermIdx idx= TermIdx(n); idx[0]=i; idx[1]=j; return idx; } 
bool SameTerm(const Myt & that) const { return  m_idx==that.m_idx; } 
bool SameV(const Myt & that) const { return  m_v==that.m_v; } 
void AddIdx(const Myt & that)
{ for (IdxTy i=0; i<n; ++i) m_idx[i]+=that.m_idx[i]; }
void InValid() { m_v=0; m_valid=false; } 
void Valid()
{
m_valid=(m_v!=0);

}

void Assign(Myt & d) const
{
d.m_valid=m_valid;
d.m_v=m_v;
d.m_idx=m_idx;

}


bool m_valid;
TermVal m_v;
TermIdx m_idx;
//MaxIdx m_max;
//Ti m_idx[n];

}; // Element

typedef Element<D,IdxTy,2> ElemTy;

template <class Tiv>
class IfVec : public std::vector<Tiv>
{
typedef IfVec Myt;
typedef  Tiv TermTy;
typedef std::vector<Tiv> Super;
typedef typename Tiv::TermIdx TermIdx;
typedef typename Tiv::TermIdx TermVal;
typedef std::map<TermIdx,TermVal> TermMap;
typedef std::vector<IdxTy> MaxIdx;
enum{INDEX_SIZE=ElemTy::INDEX_SIZE};
public:
IfVec(): Super(),m_max(INDEX_SIZE) {}
IfVec(const IdxTy n ): Super(n ),m_max(INDEX_SIZE) {}
/*
ElemTy & operator()(const IdxTy i) 
{
ElemTyj
}
*/
void push_back(const TermTy & v )
{
if (v.valid()) Super::push_back(v); 
}
const IdxTy & max(const IdxTy & i ) const { return m_max[i]; }
void index(TermMap & m ) 
{
//TermMap m=TermMap();
for(auto i=m_terms.begin(); i!=m_terms.end(); ++i)
{ m[(*i).index()]+=(*i).value(); }

}

void load(Myt & d, const TermMap & m)
{
for(auto i=m.begin(); i!=m.end(); ++i) d.push_back(TermTy((*i).first),(*i).second); 

}
void compile()
{
TermMap m=TermMap();
index(m);
clear();
load(*this,m);
}

MaxIdx m_max;

}; // IfVec

typedef IfVec<ElemTy> ElStack;


public:
// this is the sie of the matrix, degree 0 to n-1
sparse_2nomial() {}
Myt&  operator=(const Myt & that)
{
Myt&  x= *this; // (that.m_n);
return x;
}

//const IdxTy get_size() const { return m_n; } 
//Myt&  set_size(const IdxTy n )
//{
//Myt&  x= *this; // (that.m_n);
//return x;
//}

//dense_2_polynomial(const IdxTy n): Super(n,n),m_n(n) {}
template <class Tpoly, class Tval > sparse_2nomial(const Tpoly & px, const Tval  & cxp, const  Tval & cy , const Tval & c)
{
(*this)=  from_convolution( px, cxp, cy , c);
}
 sparse_2nomial(const Tv  & cx, const  Tv & cy , const Tv & c)
//:Super(3,3),m_n(3)
{
clear();
add_term(1,0,cx);
add_term(0,1,cy);
add_term(0,0,c);
compile();
//(*this)(1,0)=cx; (*this)(0,1)=cy; (*this)(0,0)=c;
}


template <class Tpoly > sparse_2nomial(const Tpoly & px, const Tpoly & py)
{
(*this)=  from_outer_product( px,py);
}


template<class Td> static 
	void get_vector( Td & dest, const Myt & x, const IdxTy k, const IdxTy dir)
{
dest=Td(x.m_terms.max(dir)+1);
for (auto ii=x.m_terms.begin(); ii!=x.m_terms.end(); ++ii)
{
const auto & xx=(*ii);
// TODO FIXME this convention is stufpu 
//MM_ONCE(" convention differes ",)
if (xx.index(1-dir)!=k) continue; 
dest[xx.index(dir)]+=xx.value();
//if (dir==0) dest[i]=x(i,k);
//else if (dir==1) dest[i]=x(k,i);
}

}
template<class Td> static 
	void set_vector( Myt & x, const Td & dest,  const IdxTy k, const IdxTy dir)
{
const IdxTy sz=dest.size();
//dest=Td(x.m_n);
for (IdxTy i=0; i<sz; ++i)
{
// TODO FIXME this is not right as the old terms will not be eliminated 
if (dir==0) x.add_term(i,k,dest[i]);
else  x.add_term(k,i,dest[i]);
//SP_ITOR_I
//if (i<sz)  {  
//if (dir==0) x(i,k)=dest[i];
//else if (dir==1) x(k,i)=dest[i];
//} else 
//{
//if (dir==0) x(i,k)=0;
//else if (dir==1) x(k,i)=0;
//}
// SP_ITOR_END
}
}
template<class Td> static 
	void add_vector( Myt & x, const Td & dest,  const IdxTy k, const IdxTy dir, const value_type & a)
{
//dest=Td(x.m_n);
const IdxTy sz=dest.size();
SP_ITOR_I
if (i<sz)  {  
if (dir==0) x(i,k)+=dest[i]*a;
else if (dir==1) x(k,i)+=dest[i]*a;
}
SP_ITOR_END

}



// this only works for same size, AFAICT but should
// check the block_matrix operator 
static void add( Myt & d, const Myt & x, const Myt & y)
{
const IdxTy szx=x.m_n;
const IdxTy szy=y.m_n;

const bool xbigger=(szx>szy);
const IdxTy szmax=(xbigger)?szx:szy;
const IdxTy szmin=(!xbigger)?szx:szy;
d= Myt(szmax);
//d=x+y;
if (xbigger)
{
SP_ITOR_IJ
d(i,j)=x(i,j);
if (i<szmin) if (j<szmin) d(i,j)+=y(i,j);
SP_ITOR_END
}
else
{
SP_ITOR_IJP
d(ip,jp)=y(ip,jp); // +y(ip,jp);
if (ip<szmin) if (jp<szmin) d(ip,jp)+=x(ip,jp);
SP_ITOR_END
}
}

 void accumulate(  const Myt & y)
{
Myt & x =(*this);
const IdxTy szx=x.m_n;
const IdxTy szy=y.m_n;

const bool xbigger=(szx>szy);
const IdxTy szmax=(xbigger)?szx:szy;
//d=x+y;
if (!xbigger)
{
SP_ITOR_IJ
x(i,j)+=y(i,j);
SP_ITOR_END
}
else
{
SP_ITOR_IJP
x(ip,jp)+=y(ip,jp);
SP_ITOR_END
}
}
template <class Tp> 
static void multiply( Myt & d, const Myt & x, const Tp & y, const bool dir, const value_type & scale=1)
{
const IdxTy szy=y.size();
if ((x.m_n==0) || (szy==0)) return; 
const IdxTy szd=1+(x.m_n-1)+(szy-1);
d= Myt(szd);
const bool dirx=dir;
//SP_ITOR_IJP
//for (IdxTy ip=0; ip<m_n; ++ip) for (IdxTy jp=0; jp<m_n; ++jp)
SP_ITOR_IJ
// TODO FIXME only do this if x(i,j) not zero 
for (IdxTy ip=0; ip<szy; ++ip)
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dirx) d(i+ip,j)+=x(i,j)*y[ip]*scale;
else d(i,j+ip)+=x(i,j)*y[ip]*scale;

//SP_ITOR_END
SP_ITOR_END

}




static void multiply( Myt & d, const Myt & x, const Myt & y, const value_type & scale=1)
{
if ((x.m_n==0) || (y.m_n==0)) return; 
const IdxTy szd=1+(x.m_n-1)+(y.m_n-1);
d= Myt(szd);
SP_ITOR_IJP
//for (IdxTy ip=0; ip<m_n; ++ip) for (IdxTy jp=0; jp<m_n; ++jp)
// TODO FIXME only do this if x(ip,jp) not zero 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
d(i+ip,j+jp)+=x(i,j)*y(ip,jp)*scale;

SP_ITOR_END
SP_ITOR_END

}




// considering x as a polynomial in variables x and y, replace the
// powers of "y" in x with the definition in y as a polynomial in 2 variables
// on or which is the x variable in x lol. 
template <class Txx> static Txx cpo(const Txx & v, const IdxTy n)
{
Txx r=1;
IdxTy i=n;
while (i>0){  r=v*r; --i; }  
return r;
}
static void composite( Myt & d, const Myt & x, const Myt & y, const IdxTy flags=1)
{
if ((x.m_n==0) || (y.m_n==0)) return; 
// may need to trim later but the biggest dimension is 
// x(m_n-1,m_n-1) and y(m_n-1,m_n-1)
const IdxTy szd=1+(x.m_n-1)*(y.m_n-1); // 1+(x.m_n-1)+(y.m_n-1);
const IdxTy imax=x.m_n;
const IdxTy jmax=x.m_n;
d= Myt(szd);
Myt yn=Myt(1);
yn(0,0)=1;
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
// the first index of "x" is for the variable being replcaed by y
const IdxTy xyz=flags;
switch (xyz)
{
case 0://{ SP_ITOR_IJP SP_ITOR_IJ d(i*ip,j+jp*i)+=x(i,j)*cpo(y(ip,jp),i); SP_ITOR_END SP_ITOR_END break; } 
{
//for (IdxTy j=0; j<jmax; ++j) {d(0,j)=x(0,j); }
for(IdxTy i=0; i<imax; ++i)
{
const IdxTy ipmax=yn.m_n;
const IdxTy jpmax=yn.m_n;

for (IdxTy j=0; j<jmax; ++j)
{
// TODO FIXME only for yn!=0 or just make these all SPARSE 
for (IdxTy ip=0; ip<ipmax; ++ip)
for (IdxTy jp=0; jp<jpmax; ++jp)

// for each row just do multiply times the thing to n power
// not tested, the dest power is just the yn power and NOT the sum as it is replcaced
d(ip,j+jp)+=x(i,j)*yn(ip,jp);

} // j 
//if ( i!=0) 
{ Myt ynext; multiply(ynext,yn,y); yn=ynext;}

} // i 
break; 
}
case 1://{SP_ITOR_IJP SP_ITOR_IJ d(i+j*ip,j*jp)+=x(i,j)*cpo(y(ip,jp),j);SP_ITOR_END SP_ITOR_END break;} 
{
//for (IdxTy j=0; j<jmax; ++j) {d(j,0)=x(j,0); }
for(IdxTy i=0; i<imax; ++i)
{
const IdxTy ipmax=yn.m_n;
const IdxTy jpmax=yn.m_n;

for (IdxTy j=0; j<jmax; ++j)
{
for (IdxTy ip=0; ip<ipmax; ++ip)
for (IdxTy jp=0; jp<jpmax; ++jp)
// this may work, the dest is NOT the sum of sources but is replcaed by the y power
// for each row just do multiply times the thing to n power
d(j+jp,ip)+=x(j,i)*yn(jp,ip);

} // j 
//if ( i!=0) 
{ Myt ynext; multiply(ynext,yn,y); yn=ynext;}

} // i 
break; 
}

default : MM_ERR(" bad flag in composite "<<MMPR(xyz))
}; // switch 

}

static void integrate( Myt & d, const Myt & x, const IdxTy dir)
{
// this is stupid to have SQUARE matrix when you do this crap 
const IdxTy szd=1+(x.m_n);
d= Myt(szd);
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dir==0) d(i+1,j)=x(i,j)/(i+1);
else if (dir==1) d(i,j+1)=x(i,j)/(j+1);
SP_ITOR_END
}
static void differentiate( Myt & y, const Myt & x, const IdxTy dir)
{
// this is stupid to have SQUARE matrix when you do this crap 
if (x.m_n==0) { y=Myt(0); return; } 
const IdxTy szd=(x.m_n)-0;
const IdxTy szmax=(x.m_n)-1;
y= Myt(szd);
SP_ITOR_IJP
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
if (dir==0) {if (ip<szmax) { y(ip,jp)=(ip+1.0)*x(ip+1,jp); }}
else {if (dir==1) {if (jp<szmax) y(ip,jp)=(jp+1.0)*x(ip,jp+1);}}
SP_ITOR_END
}




// this probably does not work right, use _new
static void integrate_e( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
MM_ONCE(" obsolete make work with _e_new", )
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(sz);

// this is stupid to have SQUARE matrix when you do this crap 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
	const IdxTy loc=(dir==0)?i:j;
	const IdxTy ord=(dir==0)?j:i;
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	// note that these could be cached but pretty easy to compute 
	//D res=mi.x_n_exp_integral(p, n, b);
	mi.x_n_erf_integral (perf, p,  ord, b,false);
	//derf.add_vector(derf,perf,loc,(dir==1)?1:0,x(i,j));
	derf.add_vector(derf,perf,ord,(dir==1)?1:0,x(i,j));
	//derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,dir,x(i,j));
//	dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	dexp.add_vector(dexp,p,ord,(dir==1)?1:0,x(i,j));
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	for (IdxTy ii=0; ii<p.size(); ++ii) p[ii]=0;
	D res=mi.x_n_exp_integral(p, ord, b,false);
	perf[0]=res*x(i,j);
//	mi.x_n_erf_integral (perf, p,  n, b);
//	derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	// dir==0 means add up and down  vary fiest index  i order is i  
	dexp.add_vector(dexp,p,ord,(dir==1)?1:0,x(i,j));
	//dexp.add_vector(dexp,p,loc,(dir==1)?1:0,x(i,j));
	if (dir==1) derf(i,0)+=perf[0];
	else derf(0,j)+=perf[0];

} // dir==1 
SP_ITOR_END
}
///////////////////////////////////////////////////////////

// this appears to work in the one case I tested
static void integrate_e_new( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(erfsz);
const IdxTy dirv=dir;
const bool rev=(dirv==0);
for (IdxTy i=0; i<sz; ++i)
{
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	mi.x_n_erf_integral (perf, p, i, b,false);
	MM_MSG(MMPR4(perf.size(),p.size(),derf.m_n,dexp.m_n))
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		derf.add_vector(derf,perf,j,dirv,v);
		dexp.add_vector(dexp,p,j,dirv,v);
	} // j 
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	D res=mi.x_n_exp_integral(p,i, b,false);
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		perf[0]=res*v;
		dexp.add_vector(dexp,p,j,dirv,v);
		if (dir==1) derf(j,0)+=perf[0];
		else derf(0,j)+=perf[0];
	}

} // dir==1 

} // i

}
/////////////////////////////////////////////////////////////////////
// user the rationzlied d "rerf" function which does integration ignoring the
// conversion of 2sqrt(b)/sqrt(pi) but user needs to add that back before calling erf 
static void integrate_e_rnew( Myt & dexp,Myt & derf, const Myt & x, const value_type & b
, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
typedef std::vector<value_type> Tpoly;
const IdxTy sz=x.m_n; // this means that the max power is (m_n-1) 
// for erf the power can go up one in erf but not in exp as int by parts
// puts (n+1)  pwer in integrand going back to n after integration 
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(erfsz);
const IdxTy dirv=dir;
const bool rev=(dirv==0);
for (IdxTy i=0; i<sz; ++i)
{
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	mi.x_n_rerf_integral (perf, p, i, b,false);
	// TODO this can not be trimmed right away as it could be accumulating although
	// should fix it 
	if (false) MM_MSG(MMPR4(perf.size(),p.size(),derf.m_n,dexp.m_n))
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		derf.add_vector(derf,perf,j,dirv,v);
		dexp.add_vector(dexp,p,j,dirv,v);
	} // j 
	}

	else if (factor==1)
{ // integrating the exp only produces a single rerf term 
//	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	D res=mi.x_n_rexp_integral(p,i, b,false);
	for (IdxTy j=0; j<sz; ++j)
	{
		const D v=rev?x(i,j):x(j,i);
		const D perf0=res*v; // 	perf[0]=res*v;
		dexp.add_vector(dexp,p,j,dirv,v);
		if (dir==1) derf(j,0)+=perf0; // perf[0];
		else derf(0,j)+=perf0; // perf[0];
	}

} // dir==1 

} // i

}
/////////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////

static void integrate_e_old( Myt & dexp,Myt & derf, const Myt & x, const value_type & b, const IdxTy factor, const IdxTy dir, const bool init_d=true)
{
 mjm_integrals mi;
MM_ONCE(" obsolete make work with _e_new", )
typedef std::vector<value_type> Tpoly;
//    typedef shape_function_integrals Pol;
 //   typedef Pol::simple_polynomial<D> Po;
//P dest,desterfexp,desterferf;
//D res=mi.x_n_exp_integral(dest, n, b);
//mi.x_n_erf_integral (desterferf, desterfexp,  n, b);
const IdxTy sz=x.m_n;
if (init_d) dexp=Myt(sz);
IdxTy erfsz=(factor==0)?(sz+1):sz;
if (init_d) derf=Myt(sz);

// this is stupid to have SQUARE matrix when you do this crap 
SP_ITOR_IJ
//for (IdxTy i=0; i<m_n; ++i) { for (IdxTy j=0; j<m_n; ++j) {
	const IdxTy loc=(dir==0)?i:j;
	const IdxTy ord=(dir==0)?j:i;
	Tpoly p(sz),perf(erfsz);
	if (factor==0)
	{
	// note that these could be cached but pretty easy to compute 
	//D res=mi.x_n_exp_integral(p, n, b);
	mi.x_n_erf_integral (perf, p,  ord, b,false);
	derf.add_vector(derf,perf,loc,(dir==0)?1:0,x(i,j));
	//derf.add_vector(derf,perf,loc,dir,x(i,j));
	//dexp.add_vector(dexp,p,loc,dir,x(i,j));
	dexp.add_vector(dexp,p,loc,dir,x(i,j));
	}

	else if (factor==1)
{
	perf=Tpoly(sz+2); // TODO FIXME this only needs to be 1 ???
	for (IdxTy ii=0; ii<p.size(); ++ii) p[ii]=0;
	D res=mi.x_n_exp_integral(p, ord, b,false);
	perf[0]=res*x(i,j);
//	mi.x_n_erf_integral (perf, p,  n, b);
//	derf.add_vector(derf,perf,loc,dir,x(i,j));
	dexp.add_vector(dexp,p,loc,dir,x(i,j));
	if (dir==0) derf(i,0)+=perf[0];
	else derf(0,j)+=perf[0];

} // dir==1 
SP_ITOR_END
}






// evaluate the 2 nomial at value v for one of the variable and 
// return a polynomial in dest
template<class Td> static void evaluate( Td & dest, const Myt & x, const value_type & v, const IdxTy dir)
{
dest=Td(x.m_n);
value_type vp[x.m_n];
vp[0]=1;
SP_ITOR_IJ
if (dir==0) dest[i]+=x(i,j)*vp[i];
else if (dir==1) dest[j]+=x(i,j)*vp[j];
SP_ITOR_END

}

value_type evaluate(  const value_type & a, const value_type & b ) const
{
const Myt & x=(*this);
const IdxTy sz=x.m_n;
if (sz==0) return 0; 
value_type vx[sz],vy[sz];
vx[0]=1; vy[0]=1;
for (IdxTy i=1; i<sz; ++i) { vx[i]=vx[i-1]*a; vy[i]=vy[i-1]*b; }
// TODO  consider making the power tables exportable 
value_type res=0;
SP_ITOR_IJ
res+=x(i,j)*vx[i]*vy[j];
SP_ITOR_END
return res;


}
// these are actually backwards but ok 
value_type sum_squares(  const value_type & x0, const value_type & x1
,  const value_type & xp0, const value_type & xp1,const value_type & d
) const 
{
value_type res=0;
value_type v=evaluate(x0-d,xp0-d); res+=v*v;
v=evaluate(x1-d,xp0-d); res+=v*v;
v=evaluate(x1-d,xp1-d); res+=v*v;
v=evaluate(x0-d,xp1-d); res+=v*v;
return res;

}
value_type evaluate_four(  const value_type & x0, const value_type & x1
,  const value_type & xp0, const value_type & xp1,
  const value_type & e00, const value_type & e11,
  const value_type & e01, const value_type & e10
, IdxTy * pflags=0
) const 
{
//value_type res=0;
//value_type v=evaluate(x0-d,xp0-d); res+=v*v;
//v=evaluate(x1-d,xp0-d); res+=v*v;
//v=evaluate(x1-d,xp1-d); res+=v*v;
//v=evaluate(x0-d,xp1-d); res+=v*v;
const Myt & x=(*this);
const IdxTy sz=x.m_n;
if (sz==0) return 0; 
bool no_x=true;
value_type vx[sz],vy[sz],vx2[sz],vy2[sz];
vx[0]=1; vy[0]=1;
vx2[0]=1; vy2[0]=1;
for (IdxTy i=1; i<sz; ++i) { vx[i]=vx[i-1]*x0; vy[i]=vy[i-1]*xp0; }
for (IdxTy i=1; i<sz; ++i) { vx2[i]=vx2[i-1]*x1; vy2[i]=vy2[i-1]*xp1; }
// TODO  consider making the power tables exportable 
value_type res=0;
SP_ITOR_IJ
if (x(i,j)==0) continue;
if (i!=0) if (j!=0)  no_x=false;
res+=x(i,j)*(vx[i]*vy[j]*e00+vx2[i]*vy2[j]*e11-vx[i]*vy2[j]*e10-vx2[i]*vy[j]*e01);
SP_ITOR_END
if (pflags!=0) *pflags=(no_x)?1:0;
return res;

}



// evaluate with a polynomial in the other variable, 
// the resulting polynomial dest is a higher order
template<class Td> static void evaluate( Td & dest, const Myt & x, const Td & vpoly, const IdxTy dir)
{
const IdxTy szx=x.m_n;
const IdxTy szp=vpoly.size();
if (szx==0) return;
if (szp==0) return;
const IdxTy sz=1+(szx-1)*(szp-1);
dest= Td(sz);
if (sz<1) return;
mjm_integrals mi;
const IdxTy vpoly_size=vpoly.size();
dest=Td(sz);
//dest[0]=x(0,0);
//value_type vp[x.m_n];
//vp[0]=1;
// composite
for (IdxTy i=0; i<szx; ++i) { 
Td fn; // =vpoly;
fn.push_back(1);
for (IdxTy k=0; k<fn.size(); ++k) {
for (IdxTy j=0; j<szx; ++j) {
if (dir==0) dest[i+k]+=x(i,j)*fn[i];
else if (dir==1) dest[i+k]+=x(i,j)*fn[j];
}
}
Td next;
mi.multiply_polynomials(next,vpoly,fn);
fn=next;

}
}

template <class Tpoly > static Myt from_outer_product(const Tpoly & px, const Tpoly & py)
{
const IdxTy szx=px.size();
const IdxTy szy=py.size();
const IdxTy sz=(szx>szy)?szx:szy;
//Tval cp[szx]; cp[0]=1;
Myt mt(sz);
Myt & x=mt;
SP_ITOR_IJ
if ((i>=szx) ||(j>=szy)) continue;
mt(i,j)=px[i]*py[j];
SP_ITOR_END
return mt;
}


#if 0
// make a _2_poly from a polynomail p(x) and composite with 
// x= cxp*x' + cy*y+ c 
template <class Tpoly, class Tval > static Myt from_convolution(const Tpoly & px, const Tval  & cx, const  Tval & cy , const Tval & c)
{
const IdxTy szx=px.size();
//Tval cp[szx]; cp[0]=1;
Tval cxp[szx], cyp[szx],cp[szx];
Myt mt(szx);
cxp[0]=1;  cyp[0]=1; cp[0]=1;
for (IdxTy i=1; i<szx; ++i) 
	{cxp[i]=cxp[i-1]*cx; cyp[i]=cyp[i-1]*cy; cp[i]=cp[i-1]*c; }
Myt & x=mt;
SP_ITOR_IJ
for (IdxTy ix=(i+j); ix<szx; ++ix) 
{
mt(i,j)+=px[ix]*cxp[i]*cyp[j]*cp[ix-i-j];

} // ix
SP_ITOR_END

return mt;
}

#endif



#if 0 

static void coef_table( Super & table, const IdxTy szx, const value_type & cx, 
	const value_type & cy, const value_type & c)
{
Super m(szx,szx,szx);
//typedef value_type Tval;
trinomial_coef(m,szx);
//Tval cx[szx], cyp[szx],cp[szx];
// any of these could be zero 
//cxp[0]=1;  cyp[0]=1; cp[0]=1;
//for (IdxTy i=1; i<szx; ++i) 
//	{cxp[i]=cxp[i-1]*cx; cyp[i]=cyp[i-1]*cy; cp[i]=cp[i-1]*c; }


table=m;
}
static IdxTy bffac(IdxTy n) { IdxTy i=1; while (n>0) { i=i*n; --n; }return i; }
static void trinomial_coefs(Super & t, const IdxTy sz)
{
const IdxTy szi=sz;
const IdxTy szj=sz;
const IdxTy szk=sz;
value_type fac=1;
for (IdxTy i=0; i<szi; ++i) { 
//fac=1;
for (IdxTy j=0; i<szj; ++j) { 
fac=1;
for (IdxTy k=0; i<szk; ++k) { 
t(i,j,k)=fac;
IdxTy bf=bffac(i)/(bffac(j)*bffac(k)*bffac(i-j-k));
if (t(i,j,k)!=bf) MM_ERR(MMPR4(t(i,j,k),bf,i,j)<<MMPR(k))
fac=fac*(i-j-k)/(k+1);
}
}
}


}

#endif


template<class Tx, class Ty, class Tn > bool  add_term( Tx & terms, Ty & ss, const Tn & t, const bool force=!true) const 
{
	if ((t!=0)||(force)) { if (terms!=0) { if ( t>=0) ss<<"+"; } ss<<t; ++terms; return true;  }
	return false; 
}
//StrTy to_string( const IdxTy flags=0)
StrTy to_string(const IdxTy flags=0,const StrTy & var=StrTy("x"),
	const StrTy & var2=StrTy("y")) const
{
std::stringstream ss;
//const int x2s=x2.size();
const Myt & x=*this;
int terms=0;
const IdxTy sz=x.m_n;
SP_ITOR_IJ
const IdxTy ii=sz-i-1;
const IdxTy jj=sz-j-1;
const value_type & t=x(ii,jj);
 if ( add_term(terms,ss,x(ii,jj)))   
{
//if (jj!=0) {if (jj==1) {ss<<"x";} else   {ss<<"x^"<<jj;} } 
if (ii!=0) {ss<<var; if (ii!=1) {ss<<"^"<<ii;}  } 
if (jj!=0) {ss<<var2; if (jj!=1) {ss<<"^"<<jj;}  } 
//if (ii!=0) {if (ii==1) {ss<<"y";} else { ss<<"y^"<<ii;} } 

} 
//if (x2s>0) {if (x2[0]!=0) {  if (terms!=0) ss<<"+"; ss<<x2[0];} }
//if (x2s>0) {add_term(terms,ss,x2[0]); }
SP_ITOR_END
return ss.str();
}



void clear() { m_terms.clear(); }
void compile()
{m_terms.compile();
}

void add_term( const IdxTy & i, const IdxTy& j, const Tv & v)
{
//m_terms.push_back(ElemTy(1,0,cx));
m_terms.push_back(ElemTy(i,j,v));
}

private:
ElStack m_terms;
//IdxTy m_n;

}; // sparse_2nomial

#ifdef TEST_SPARSE_2NOMIAL
//////////////////////////////////////////////////////////
int main(int argc,char **args)
{

typedef sparse_2nomial Myt;
Myt mi;
typedef unsigned int IdxTy;
typedef double D ;
//const IdxTy nu=1;
const IdxTy flags=0;
typedef std::vector<D> Dest;
Dest dest,desthc;
MM_MSG(" args are v,mu,nu,kmax");
int  pos=1;
const D v=(argc>pos)?atof(args[pos]):0; ++pos;
const IdxTy  mu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  nu=(argc>pos)?atoi(args[pos]):0; ++pos;
const IdxTy  kmax=(argc>pos)?atoi(args[pos]):0; ++pos;
const D  z0=(argc>pos)?atof(args[pos]):0; ++pos;
return 0;
} // main

#endif


#endif // guard 

