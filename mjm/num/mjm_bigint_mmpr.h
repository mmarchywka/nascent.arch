#ifndef MJM_BIGINT_MMPR_H__
#define MJM_BIGINT_MMPR_H__

#include "mjm_globals.h"
#include "mjm_generic_iterators.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>
#include "gmp.h"
#include "gmp-impl.h"

/*

Copied from mjm_block_matrix.h because the libmesh SparseMatrix is 
confusing and can't get sparsity pattern easily...

 g++ -DTEST_BIGINT_MMPR_MAIN__ -Wall -I.. -I gmp/gmp-6.1.2   -std=gnu++11  -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -Lgmp/gmp-6.1.2/.libs -lgmp -x c++ mjm_bigint_mmpr.h 2>fudd


*/
class mjm_bigint_mmpr ;

namespace mjm_int_ops {

double ratio(const mjm_bigint_mmpr & n, const mjm_bigint_mmpr & d);

};

class mjm_bigint_mmpr 
{
typedef mjm_bigint_mmpr Myt;
typedef mpz_t BigIntTy;
public:
typedef double D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
//typedef std::vector<IdxTy> LocTy;
//typedef mjm_generic_itor<2> LoopItor2;
//typedef LocTy location_type;
mjm_bigint_mmpr() {Init();}
mjm_bigint_mmpr(const int n) {Init(); Set(n); }
mjm_bigint_mmpr(const unsigned int n) {Init(); Set(n); }
// long maybe superfluous here, and then ref varies 
mjm_bigint_mmpr(const long unsigned int & n) {Init(); Set(n); }
mjm_bigint_mmpr(const double & n) {Init(); Set(n); }
// this needs a fudding copy ctor or else it copies a fudding pointer
mjm_bigint_mmpr(const Myt & that) {Init(); Set(that.m_n); }


~mjm_bigint_mmpr() {Dtor();}

StrTy string(const IdxTy base=10) const
{
return String(base);
}
// idiosyncratics 
double ratio( const Myt & d ) const { return Ratio(d); } 
// conversions
explicit operator double() const { return Double(); } 
explicit operator long unsigned int () const { return Lui(); } 
explicit operator unsigned int () const { return Lui(); } 
explicit operator int () const { return Li(); } 
explicit operator signed long int () const { return Li(); } 
Myt&  operator=(const Myt&  rhs) { Set(rhs.m_n); return *this; }
Myt&  operator=(const int rhs) { Set(rhs); return *this; }
Myt&  operator=(const unsigned int rhs) { Set(rhs); return *this; }
Myt&  operator=(const double & rhs) { Set(rhs); return *this; }
Myt&  operator--() { (*this)-=Myt(1);  return *this; }
Myt&  operator++() { (*this)+=Myt(1);  return *this; }
// f is for fudd 
//Myt&  operator=(const float & rhs) { Set(rhs); return *this; }

bool operator>(const Myt & that ) const { return Gt(that); } 
bool operator<(const Myt & that ) const { return Lt(that); } 
bool operator<=(const Myt & that ) const { return Le(that); } 
bool operator==(const Myt & that ) const { return Eq(that); } 
bool operator!=(const Myt & that ) const { return !Eq(that); } 

Myt  operator<<(const unsigned int rhs) const
{ Myt d; Lshift(d,*this,rhs); return d; }
Myt  operator<<(const Myt & rhs) const
{ Myt d; Lshift(d,*this,(unsigned int)(rhs)); return d; }

Myt&  operator<<=(const unsigned int rhs) 
{ Lshift(*this,*this,rhs); return *this; }
Myt&  operator<<=(const Myt & rhs) 
{  Lshift(*this,*this,(unsigned int)(rhs)); return *this; }





Myt  operator>>(const unsigned int rhs) const
{ Myt d; Rshift(d,*this,rhs); return d; }

Myt&  operator>>=(const unsigned int rhs) 
{ Rshift(*this,*this,rhs); return *this; }
Myt&  operator>>=(const Myt & rhs) 
{  Rshift(*this,*this,(unsigned int)(rhs)); return *this; }

Myt  operator+(const Myt & rhs) const
{ Myt d; Add(d,*this,rhs); return d; }
Myt&  operator+=(const Myt & rhs) 
{ Add(*this,*this,rhs); return *this; }
Myt  operator-(const Myt & rhs) const
{ Myt d; Sub(d,*this,rhs); return d; }
Myt&  operator-=(const Myt & rhs) 
{ Sub(*this,*this,rhs); return *this; }

Myt  operator*(const Myt & rhs) const
{ Myt d; Mult(d,*this,rhs); return d; }
Myt&  operator*=(const Myt & rhs) 
{ Mult(*this,*this,rhs); return *this; }


Myt  operator/(const Myt & rhs) const
{ Myt d; Div(d,*this,rhs); return d; }
Myt&  operator/=(const Myt & rhs) 
{ Div(*this,*this,rhs); return *this; }

Myt  operator%(const Myt & rhs) const
{ Myt d; Mod(d,*this,rhs); return d; }


Myt  operator&(const Myt & rhs) const
{ Myt d; And(d,*this,rhs); return d; }
Myt&  operator&=(const Myt & rhs) 
{ And(*this,*this,rhs); return *this; }

// results differ from Myt  operators
double  operator*(const double & rhs) const {return double(*this)*rhs; }
Myt  operator*(const int & rhs) const {return (*this)*Myt(rhs); }
Myt  operator*(const unsigned int & rhs) const {return (*this)*Myt(rhs); }
double  operator+(const double & rhs) const {return double(*this)+rhs; }
Myt  operator+(const int & rhs) const {return (*this)+Myt(rhs); }
Myt  operator+(const unsigned int & rhs) const {return (*this)+Myt(rhs); }
double  operator-(const double & rhs) const {return double(*this)-rhs; }
Myt  operator-(const int & rhs) const {return (*this)-Myt(rhs); }
// wtf unary minus??? 
Myt  operator-() const {return (*this)*(-1); }
double  operator/(const double & rhs) const {return double(*this)/rhs; }
Myt  operator/(const int & rhs) const {return (*this)/Myt(rhs); }
Myt  operator/(const unsigned int & rhs) const {return (*this)/Myt(rhs); }
bool  operator>(const double & rhs) const {return double(*this)>rhs; }
bool  operator>(const int & rhs) const {return (*this)>Myt(rhs); }
bool  operator>(const unsigned int & rhs) const {return (*this)>Myt(rhs); }
bool  operator<(const double & rhs) const {return double(*this)<rhs; }
bool  operator<(const int & rhs) const {return (*this)<Myt(rhs); }
bool  operator<(const unsigned int & rhs) const {return (*this)<Myt(rhs); }
bool  operator==(const double & rhs) const {return double(*this)==rhs; }
bool  operator==(const unsigned int & rhs) const {return (*this)==Myt(rhs); }
bool  operator==(const int & rhs) const {return (*this)==Myt(rhs); }
//double  operator*=(const double & rhs) { Mult(*this,*this,rhs); return *this; }



void SetString( const char * s,const IdxTy base=10)
{
mpz_set_str(m_n,s,base);
}
// fudding friend shot is fudded 
private:
StrTy String(const IdxTy base) const 
{
	char * str = mpz_get_str ((char *) 0, base, m_n);
	StrTy s(str);
	(*__gmp_free_func) (str, strlen (str) + 1);
	return s;
}

void Init()
{
mpz_init(m_n);

}
void Dtor()
{
mpz_clear(m_n);
}
void Set(const BigIntTy & n) { mpz_set(m_n,n); }
void Set(const int n) { mpz_set_si(m_n,n); }
void Set(const unsigned int n) { mpz_set_ui(m_n,n); }
void Set(const long unsigned int & n) { mpz_set_ui(m_n,n); }
void Set(const double & n) { mpz_set_d(m_n,n); }
// WTF is f????
//void Set(const float & n) { mpz_set_f(m_n,n); }
double Ratio(const Myt & d) const
{
signed long int en,ed;
const double mn=mpz_get_d_2exp(&en,m_n);
const double md=mpz_get_d_2exp(&ed,d.m_n);
signed long int de=en-ed;
double x= mn/md;
while (de>0) {x=x*2.0; --de; } 
while (de<0) {x=x/2.0; ++de; } 
return x; 
}

double Double() const
{
double x= mpz_get_d(m_n);
return x; 
}


long unsigned int Lui() const
{
return mpz_get_ui(m_n);
}
long  int Li() const
{
return mpz_get_si(m_n);
}

void Lshift(Myt & d, const Myt & s, const unsigned int rhs) const
{ mpz_mul_2exp(d.m_n,s.m_n,rhs); }

void Rshift(Myt & d, const Myt & s, const unsigned int rhs) const
{ mpz_fdiv_q_2exp(d.m_n,s.m_n,rhs); }

bool Gt(const Myt & that ) const 
{ int x=mpz_cmp(this->m_n,that.m_n); return (x>0); }

bool Lt(const Myt & that ) const 
{ int x=mpz_cmp(this->m_n,that.m_n); return (x<0); }
bool Le(const Myt & that ) const 
{ int x=mpz_cmp(this->m_n,that.m_n); return !(x>0); }


bool Eq(const Myt & that ) const 
{ int x=mpz_cmp(this->m_n,that.m_n); return (x==0); }




// mpz_tdiv_q (quotient2, dividend, divisor);
//      mpz_tdiv_r (remainder2, dividend, divisor);
void Add(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_add(d.m_n,r1.m_n,r2.m_n); }
void Sub(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_sub(d.m_n,r1.m_n,r2.m_n); }


void Div(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_tdiv_q(d.m_n,r1.m_n,r2.m_n); }
void Mod(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_tdiv_r(d.m_n,r1.m_n,r2.m_n); }
void Mult(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_mul(d.m_n,r1.m_n,r2.m_n); }

void And(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_and(d.m_n,r1.m_n,r2.m_n); }


void Or(Myt & d, const Myt & r1, const Myt & r2) const
{ mpz_ior(d.m_n,r1.m_n,r2.m_n); }





BigIntTy m_n;
// this is already heavy, just put more state crap here lol


friend   
    std::ostream&  operator<<(std::ostream& os, const mjm_bigint_mmpr  & r);
// this is a huge fudd nothing FUDDNG works FUDD 
//friend   
    std::istream&  operator>>(  mjm_bigint_mmpr  & r);
friend   

double mjm_int_ops::ratio(const mjm_bigint_mmpr & n, const mjm_bigint_mmpr & d);
//double ratio(const mjm_bigint_mmpr & n, const mjm_bigint_mmpr & d);

}; // mjm_bigint_mmpr

std::ostream&  operator<<(std::ostream& os,  const  mjm_bigint_mmpr & r)
{ os<<r.string(); return os; }

// TODO this should use the base indicted for the stream 
std::istream&  operator>>( std::istream& fudd,  mjm_bigint_mmpr & r)
//{ StrTy s; (*this)>>s; r.SetString(s.c_str(),10); return (*this); }
{ StrTy s; fudd>>s; r.SetString(s.c_str(),10); return (fudd); }

namespace mjm_int_ops
{

double ratio(const mjm_bigint_mmpr & n, const mjm_bigint_mmpr & d) 
{
signed long int en,ed;
const double mn=mpz_get_d_2exp(&en,n.m_n);
const double md=mpz_get_d_2exp(&ed,d.m_n);
signed long int de=en-ed;
double x= mn/md;
while (de>0) {x=x*2.0; --de; } 
while (de<0) {x=x/2.0; ++de; } 
return x; 
}
}; // mjm_int_ops

#ifdef TEST_BIGINT_MMPR_MAIN__
#define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
#include "mjm_rational.h"
int main(int argc, char ** argv)
{
typedef mjm_bigint_mmpr T;

T x,y; // =1;
x=1;
y=2;
T z=x+y;
//x=99+x;
MM_MSG(MMPR3(x,y,z))

return 0;
}

#endif // test sparse

#endif // guard

