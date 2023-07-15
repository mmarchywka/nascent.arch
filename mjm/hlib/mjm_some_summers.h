#ifndef MJM_SOME_SUMMERS__H__
#define MJM_SOME_SUMMERS__H__
#include "common_includes.h"
#include "mjm_globals.h"

/*
Find various sums of coloumb like points and their derivates etc.
For jdftx ewald alternatives. 

*/

// mjm_sum_list and matrix are really just interfaces for the input
// and output things. j
class mjm_sum_traits
{


};
// can pickuip from tr...
class mjm_sum_typedefs
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::IdxTy IdxTy;
typedef double D;
typedef std::stringstream SS;

static OsTy & info() { return std::cout; }


};


// the sums operate on CARTESIAN coords so notion of vector really meaningless
class mjm_sum_list : public mjm_sum_typedefs
{
typedef mjm_sum_list Myt;
typedef std::vector<D>  Item;
typedef std::vector<Item> List;

List m_items;
IdxTy m_coords; // typically 3
public:
enum { X=0, Y=1, Z=2,ZEFF=3,MASS=4};
mjm_sum_list(): m_coords(3) {}
IdxTy coords() const { return size()*m_coords; } 
IdxTy size() const  { return m_items.size(); }
const Item & operator[](const IdxTy n) const { return m_items[n];}
// get attribute for item n 
const D &  get(const IdxTy n, const IdxTy attr) const { return m_items[n][attr];}
// get coord m in order 
const D & get(const IdxTy m) const { return m_items[m/m_coords][m%m_coords]; } 
void add(const Item & a) { m_items.push_back(a); } 

}; // mjm_sum_list
class mjm_sum_matrix :public mjm_sum_typedefs
{
typedef mjm_sum_matrix Myt;
typedef std::vector<D>  Row;
typedef std::vector<Row> Grid;
Grid m_matrix;

template <class Ty> void sz(Ty & x, const IdxTy sz) 
{
while ( x.size()<sz) x.push_back(typename Ty::value_type()); 

}
public:
void size(const IdxTy m, const IdxTy n)
{
sz(m_matrix,m);
for ( IdxTy i=0; i<m_matrix.size(); ++i) sz(m_matrix[i],n); 

}


}; // mjm_sum_matrix


class mjm_summer :public mjm_sum_typedefs
{

public:
template <class Ti > D d( const Ti & i1, const Ti & i2)
{
enum { X=Ti::X, Y=Ti::Y, Z=Ti::Z};
const D dx=i1[X]-i2[X];
const D dy=i1[Y]-i2[Y];
const D dz=i1[Z]-i2[Z];
return dx*dx+dy*dy+dz*dz;

}
template <class Ti > D inter( const Ti & i1, const Ti & i2, const IdxTy c1, const IdxTy c2)
{

return 0;
}
// change in energy of coord i due to change of location of coord j 
template <class Tlist, class Tmat> void coulomb_force(Tmat & d, const Tlist & list)
{



}



}; // mjm_summer

#endif

