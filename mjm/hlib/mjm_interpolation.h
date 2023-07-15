#ifndef MJM_INTERPOLANTS_H__
#define MJM_INTERPOLANTS_H__

#include "mjm_globals.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_sparse_matrix.h"

// eventually want to solve complicated models for fits
// #define USING_PETSC
#ifdef USING_PETSC_INTERPOLATION
// deviate from shooting to matrix 
#include "mjm_petsc_util.h"
#endif

/*
g++ -DTEST_INTERPOLANTS__ -Wall  -std=gnu++11 -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -gdwarf-3 -Wno-unused-function -I.. -x c++ mjm_interpolation.h

*/
class mjm_interpolants 
{
public:

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // Tr

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;

typedef std::vector<D> VecTy;

typedef mjm_interpolants Myt;

// finds a polynaminl \sigma {c_i*x_i} to estimate
// x , x!= x_i for any i. c_i= A*\pi{ (x-x_i) } 
// I thnk this is lagrange interpolation 
// ideally suited for rational number implementation 
mjm_interpolants() : m_points(0), m_A(1)
{}
// location value
template <class Ty>
void samples(const Ty & pairs)
{


}

template <class Tx,class Ty>
void samples(const Tx & x, const Ty & y)
{
m_locations=x;
m_values=y;
}
template <class Tx,class Ty>
void sample(const Tx & x, const Ty & y)
{
m_locations.push_back(x);
m_values.push_back(y);
}



// find coefficient i at position p 
D coef(const D & pos, const IdxTy i) const
{
myassert(i<m_points, " index out of range ",__LINE__);
const IdxTy sz=points();
D x=1;
for (IdxTy j=0; j<sz; ++j)
{
if (i==j) continue;
// location warping could be done once in compile, 
// and once for the pos on each call 
x=x*(pos-m_locations[j]);
} // j
// need rat numbers here lol 
return x*m_den[i]; 
}
// normally the calculation is expensive who cares about vft lu
// while slower, fwd and rev virtuals for each value lu work too 
virtual D value(const D & pos)
{
D x=0;
const IdxTy sz=points();
for (IdxTy i=0; i<sz; ++i) x+=coef(pos,i)*m_values[i];

return x/m_A;
}

template <class Tx, class Ty>
bool  myassert(const bool cond, const Tx & msg, const Ty & line) const
{
if (cond) return true;
MM_ERR(" assertion failed for "<<msg<<" at line "<<line)
return false;

}
// O(0) although probably just warping functionf from locations
// would work too and avoid re-writing these loops.  

virtual void compile()
{
m_den.clear();
m_points=m_values.size();
const IdxTy sz=points();
for (IdxTy i=0; i<sz; ++i)
{
D x=1;

for (IdxTy j=0; j<sz; ++j)
{
if (i==j) continue;
x=x*(m_locations[i]-m_locations[j]);
} // j
// need rat numbers here lol 
m_den.push_back(1.0/x);
} // i 
verify_at_points();
}

D  verify_at_points()
{
D x=0; 
const IdxTy sz=points();
for (IdxTy i=0; i<sz; ++i)
{
	const D pos=m_locations[i];
	D v=value(pos);
	const D pv=m_values[i];
	D l2=(v-pv)*(v-pv);
	 if (l2!=0)
	{
		MM_ERR(" interpolation in exact "<<i<<" pos "<<pos<<" "<<v<<" vs "<<pv)
		x+=l2;
	}
}
if (x!=0) { MM_ERR(" interpolation error l2 "<<x)}

return x;
}


const IdxTy points() const { return m_points; }

IdxTy m_points;
VecTy m_values;
VecTy m_locations;
VecTy m_den;
D m_A;

}; // mjm_interpolants


class mjm_exp_interp : public mjm_interpolants
{
typedef mjm_interpolants Super;
typedef mjm_exp_interp Myt;
typedef Super Tr;
public:

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;
//typedef std::vector<D> VecTy;
//typedef mjm_interpolants Myt;


D value(const D & pos) const
{
D x=0;
const IdxTy sz=points();
for (IdxTy i=0; i<sz; ++i) x+=coef(pos,i)*fwd(m_values[i]);

return rev(x/m_A);
}

D fwd( const D & x ) const { return ::log(x); } 
D rev( const D & x ) const { return ::exp(x); } 

}; // mjm_exp_interp



#ifdef TEST_INTERPOLANTS__

int main(int argc, char ** argv)
{
typedef mjm_interpolants Myt;
typedef double D;
typedef unsigned int IdxTy;
Myt x;
typedef std::vector<double> V;
V pos,val;
pos.push_back(0); val.push_back(1);
pos.push_back(2); val.push_back(3);
pos.push_back(10); val.push_back(10);
pos.push_back(4); val.push_back(2);
pos.push_back(15); val.push_back(0);

x.samples(pos,val);
x.compile();

for (IdxTy i=0; i<10; ++i)
{
const D loc=1.0*i/5;
D valx=x.value(loc);
MM_MSG(" location "<<loc<<" val "<<valx)
} 


Myt y;
V pos2,val2;
pos2.push_back(1); val2.push_back(2);
pos2.push_back(0); val2.push_back(1);
y.samples(pos2,val2);
y.compile();
const D loc=.25;
D valx=y.value(loc);
MM_MSG(" 2 point location "<<loc<<" val "<<valx)



return 0;
}

#endif // TEST_INTERPOLANTS__



#endif // MJM_INTERPOLANTS


