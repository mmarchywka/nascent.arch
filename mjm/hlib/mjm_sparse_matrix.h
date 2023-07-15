#ifndef MJM_SPARSEMATRIX_H__
#define MJM_SPARSEMATRIX_H__

#include "mjm_globals.h"
#include "mjm_generic_iterators.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*

Copied from mjm_block_matrix.h because the libmesh SparseMatrix is 
confusing and can't get sparsity pattern easily...

g++ -DTEST_SPARSE_MATRIX_MAIN__ -Wall -I.. -std=gnu++11  -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -x c++ mjm_sparse_matrix.h 

g++ -DTEST_OUTER_MATRIX_MAIN__ -Wall -I.. -std=gnu++11  -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -x c++ mjm_block_matrix.h 


 g++ -DTEST_OUTER_MATRIX_MAIN__ -Wall -I.. -std=gnu++11  -x c++ mjm_block_matrix.h


 2330  g++ -DTEST_BLOCK_MATRIX_MAIN__ -Wall -I.. -x c++ mjm_block_matrix.h 
 2331  ./a.out

*/
template <class Ty> class mjm_sparse_matrix 
{
typedef mjm_sparse_matrix Myt;
public:
typedef Ty D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
typedef mjm_generic_itor<2> LoopItor2;
typedef LocTy location_type;


class Row
{
typedef mjm_sparse_matrix Tr;
typedef Row Myt;
typedef Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::SsTy SsTy;
typedef std::vector<D> RowData;
typedef std::vector<IdxTy> IdxData;
enum { BAD=~0};

public:
Row(): m_dummy(0) {}
D & operator() (const int i)  { return (*this)((IdxTy) i); }
D & operator() (const IdxTy i) 
{
	IdxTy j=find(i,0,true);
	if (j==IdxTy(BAD)) { 
		MM_ERR(" find fails to insert "<<j<<" for "<<i<<" at size "<<size())
		j=insert(i); }
	if (m_idx[j]!=i)
		MM_ERR(" bad index lookup returns "<<j<<" which is "<<m_idx[j]<<" != "<<i)
	return m_data[j];
}
bool exists(const IdxTy i) { return (IdxTy(BAD)!=find(i,0,false)); }

bool check()
{
	bool ok=true;
	if (m_idx.size()!=m_data.size())
		 { ok=false; MM_ERR(" size error "<<m_idx.size()<<" vs "<<m_data.size()) } 
	if (m_idx.size()!=0)
	{
		if (m_idx[0]>=m_idx[m_idx.size()-1])
			{ ok=false; MM_ERR(" index errors "<<m_idx[0]<<" vs "<<m_idx[m_idx.size()-1]) }

	}
	return ok;
}

IdxTy insert(const IdxTy i)
{


return BAD;
}
// insert an entry for idx at element location i
// user refs idx for our entry i
void insert_at(const IdxTy idx, const IdxTy i)
{

	const IdxTy sz=size(); 
	m_idx.push_back(0);
	m_data.push_back(0);
	for(IdxTy j=sz; j>i; --j) 
	{
		m_idx[j]=m_idx[j-1];
		m_data[j]=m_data[j-1];
	}
	if (i>=m_idx.size()) 
		{MM_ERR(" insert failed wtf "<<i<<" vs "<<idx) }
	m_idx[i]=idx;
	m_data[i]=0; // wtf the += operator was inserting and then adding to the old data lol fudd. 
}

void insert_block_at(const IdxTy idx, const IdxTy i, const IdxTy len)
{
// could do a lot better.. 
	const IdxTy sz=size(); 
	for (IdxTy j=0; j<len; ++j)
	{	m_idx.push_back(0); m_data.push_back(0);}
	const IdxTy szeff=sz+len-1;
	const IdxTy steff=i+len-1;
	for(IdxTy j=szeff; j>steff; --j) 
	{
		m_idx[j]=m_idx[j-len];
		m_data[j]=m_data[j-len];
	}
	if (i>=m_idx.size()) 
		{MM_ERR(" insert failed wtf "<<i<<" vs "<<idx) }
	for(IdxTy j=0; j<len; ++j) m_idx[i+j]=idx+j;
}



IdxTy find(const IdxTy idx, const IdxTy start=0, const bool or_insert=true)
{
const IdxTy sz=size(); 
for (IdxTy i=start; i<sz; ++i)
{
if (m_idx[i]==idx) {return i; }
if (m_idx[i]>idx)
{
if (or_insert) {insert_at(idx,i); return i; }
return BAD;
}

}
if (!or_insert)return BAD;
insert_at(idx,sz);
return sz; 


}
const IdxTy size() const { return m_idx.size(); }

StrTy string() const 
{
	SsTy s;
	dump(s,3);
	return s.str();
}
void multiply( const D & v)
{
	const IdxTy sz=size();
	for (IdxTy i=0; i<sz; ++i)  m_data[i]*=v;

}
// this is for the fudding row, now column 
template < class Os> 
	void dump(Os & os, const IdxTy flags, const IdxTy idx=0) const 
{
const StrTy sep=" ";
const IdxTy sz=size();
switch (flags)
{
case 0:	{  // serial col value
			for (IdxTy i=0; i<sz; ++i) 
				os<<i<<sep<<m_idx[i]<<sep<<m_data[i]<<CRLF;
			break; 
		}
case 1:	{  // x(i)=value; cr every few lines
			for (IdxTy i=0; i<sz; ++i) 
			{
				os<<"x("<<m_idx[i]<<")="<<m_data[i]<<"; ";
				if (((i&3)==3) || ( i==(sz-1))) os<<CRLF;
			}
			break; 
		}

case 2:	{  // x(label,i)=value;
			for (IdxTy i=0; i<sz; ++i) 
			{
				os<<"x("<<idx<<","<<m_idx[i]<<")="<<m_data[i]<<"; ";
				if (((i&3)==3) || ( i==(sz-1))) os<<CRLF;
			}
			break; 
		}

case 3:	{  // x(i)=value all on one line 
			for (IdxTy i=0; i<sz; ++i) 
			{
				os<<"x("<<m_idx[i]<<")="<<m_data[i]<<"; ";
			}
			break; 
		}
/*
case 4:	{  // log - unfortunately this requests that the data tpe
			// be convertable into a doubl3... 
			SsTy ss; 
			const double scale=1.0/::log(10);
			for (IdxTy i=0; i<sz; ++i) 
			{
				const D  & v=m_data[i];
				double vd=::fabs(double(v));
				if (vd<1e-100) vd=1e-100;	
				int x= int(scale(*::log(vd)));
				ss<<"x("<<m_idx[i]<<")="<<x<<"; ";
			}
			os<<ss.str();
			break; 
		}
*/
case 5:	{  // x(i)=value all on one line 
			D mmin=0;
			D mmax=0;
			if (sz!=0) { mmin=m_data[0]; mmax=mmin; } 
			IdxTy lmin=0, lmax=0; 
			for (IdxTy i=0; i<sz; ++i) 
			{
				const D & v=m_data[i];
				if (v>mmax) { lmax=m_idx[i]; mmax=v; }
				if (v<mmin) { lmin=m_idx[i]; mmin=v; }
			}

				os<<"min "<<mmin<<"@"<<lmin<<" and max "<<mmax<<"@"<<lmax<<"; "<<CRLF;
			break; 
		}



default :	{ 
	MM_ERR(" outputing default format for sparse matric "<<flags) 
			for (IdxTy i=0; i<sz; ++i) 
				os<<i<<sep<<m_idx[i]<<sep<<m_data[i]<<CRLF;
			break; 
		}

} // switch 

}

bool abs_gt( const D & x, const D & y) 
{
return (x*x)>(y*y); 
}
const D & max_abs()
{
	IdxTy loc=0;
	const IdxTy sz=size(); 
	if (sz==0) return m_dummy;
	for (IdxTy i=0; i<sz; ++i)
	{
		if (abs_gt(m_data[i],m_data[loc])) loc=i;
	}
	return m_data[loc];
}
// read only porbably a good idea since this is a temp
// impl 
const D & data(const IdxTy i) const { return m_data[i];}
const IdxTy & index(const IdxTy i) const { return m_idx[i];}

protected:
//IdxTy m_first, m_size;
IdxData m_idx;
RowData m_data;
D m_dummy;

friend // template <class Td> 
	std::ostream&  operator<<(std::ostream& os, const mjm_sparse_matrix<D>::Row  & r);

}; // Row

typedef std::map<IdxTy, Row> RowMap;

/*
In reality, a better one would use blocks but then you 
need to deal with block diagonal plus some ad hocs or conservation
equantions that may overlap. 
*/

mjm_sparse_matrix(  ): m_ok_to_add_elements(true)
,m_dummy(0)
//: m_dims(dimwit(1)),m_sizes(elements(m_dims))
//,m_top(m_dims.size()-1),m_ptr(valloc())
{
}

~mjm_sparse_matrix() {} // { delete [] m_ptr; m_ptr=0; }

D & operator()(const IdxTy i, const IdxTy j)
{
if (!m_ok_to_add_elements) if(!exists(i,j))
{
MM_ERR(" access to "<<i<<" "<<j<<" not exist")
return m_dummy;
}
 
return m_rows[i](j);
}

Row & operator()(const IdxTy i)
{
if (!m_ok_to_add_elements) if(!exists(i))
{
MM_ERR(" access to row "<<i<<" not exist")
return m_dummy_row;
}
return m_rows[i];
}
const Row & operator()(const IdxTy i ) const
{
if(!exists(i))
{
MM_ERR(" access to row "<<i<<" not exist")
return m_dummy_row;
}
// const fudd 
//return m_rows[i];
return (*(m_rows.find(i))).second;
}

const IdxTy size() const { return m_rows.size(); } 


const bool exists(const IdxTy i ) const 
{
// should use the find for the subsequent look up wtf
return (m_rows.find(i)!=m_rows.end()) ;
}

const bool exists(const IdxTy i, const IdxTy j)
{
// should use the find for the subsequent look up wtf
if (m_rows.find(i)==m_rows.end()) return false;
return  (m_rows[i].exists(j)) ;

}

#define MAT_ITOR for (IdxTy i=0; i<rows; ++i) for (IdxTy j=0; j<cols; ++j)
// this uses brackets for the dof index syntax so it will not
// work with libmesh but should work with std::vector
template <class Tm, class Td>
bool update(const Tm & mat, const IdxTy rows, const IdxTy cols, const Td & dofs, const IdxTy op)
{
if (!m_ok_to_add_elements)
{
	MM_ERR(" updating elements to closed matrix is not checked "<<op<<" "<<rows<<" "<<cols)
	return false;
} 
switch (op)
{
case 0: { MAT_ITOR m_rows[dofs[i]](dofs[j])=mat(i,j); break; }
case 1: { MAT_ITOR m_rows[dofs[i]](dofs[j])+=mat(i,j); break; }
default: { MM_ERR(" sparse mat update op invalid "<<op)  }
} // switch 
return true;
}
#undef MAT_ITOR

const D & max_abs_row(const IdxTy row)
{ return m_rows[row].max_abs();}

template < class Os> void dump(Os & os, const IdxTy flags, const IdxTy idx=0)
{
for( auto ii=m_rows.begin(); ii!=m_rows.end(); ++ii)
	(*ii).second.dump(os,2,(*ii).first);

}
void zero_row(const IdxTy dof) { m_rows.erase(dof); }
// create diagonal entries for each given dof in the vector
template <class Tv> void zero_rows(const Tv & dofs, const D & diag)
{
for (auto i=dofs.begin(); i!=dofs.end(); ++i)
{
// easiest way to zero it
m_rows[*i]=Row();
// set diag element 
m_rows[*i](*i)=diag;


}

}
void clear() { m_rows.clear(); } 

// members 
bool m_ok_to_add_elements;
RowMap m_rows;
D m_dummy;
Row m_dummy_row;

class sequential_iterator;
friend class Myt::sequential_iterator;

// friend classes and methods that use them
// currently no const version
class sequential_iterator
{

typedef mjm_sparse_matrix Tr;
typedef mjm_sparse_matrix Tgt;
typedef  sequential_iterator Myt;
typedef Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::SsTy SsTy;
typedef std::vector<D> RowData;
typedef std::vector<IdxTy> IdxData;
typedef typename Tgt::RowMap::iterator RowItor;
enum { BAD=~0};

public:

typedef  IdxData LocTy;
sequential_iterator( Tgt & mat) {Setup(mat);}
 operator bool() const { return m_ok; } 
 bool ok() const { return m_ok; } 
// this is stupid, it should just look it up by serial not user dense index 
D &  operator *() const { return (*m_row_itor).second(m_j); } 
//private:


Myt & inc()
{
	++m_col;
	while ((!col_ok())&&(m_ok)) inc_row();
	if (m_ok)
	{
		const auto ri=(*m_row_itor).second;
		m_j=ri.index(m_col);
	}
//MM_MSG(" inc "<<m_col<<" of "<<m_cols<<" and "<<m_j)
return *this; 
}
Myt & operator++( )  {inc(); return *this;}

void Setup( Tgt & mat)
{

	m_row_itor=mat.m_rows.begin();
	m_row_end=mat.m_rows.end();
	m_rows=mat.m_rows.size();
	m_row=0;
	m_col=0;
	if (!row_ok()) { m_ok=false; }	
	else
	{
		m_row=0;
		m_col=0;
		// need to check for zero len although not really possible 
		m_cols=(*m_row_itor).second.size();
		m_ok=true;
		while ((m_cols==0)&&(m_ok)) inc_row();
		if (ok())
		{
			const auto ri=(*m_row_itor).second;
			m_cols=ri.size();
			m_i=(*m_row_itor).first;
			m_j=ri.index(m_col);
		}
	}	
//	MM_ONCE(" delete this msg  setup itor "<<m_rows<<" "<<m_cols<<" "<<m_row<<" "<<m_col<<" "<<ok(),)
//	return *this;
}
// these are the mat indicies
const IdxTy & i() const { return m_i; }
const IdxTy & j() const { return m_j; }

IdxData loc() const 
{
IdxData x(2);
x[0]=i(); x[1]=j();
return x;
}
// these are SERIAL numbers 
const IdxTy & row() const { return m_row; }
const IdxTy & col() const { return m_col; }

bool row_ok() const { return m_row_itor!=m_row_end; }
bool col_ok() const { return m_col<m_cols; }

void inc_row()
{
	++m_row_itor; ++m_row; m_col=0;
	if (!row_ok()) { m_ok=false; return ; }
	const auto ri=(*m_row_itor).second;
	m_cols=ri.size();
	m_i=(*m_row_itor).first;
	m_j=ri.index(m_col);
}
StrTy to_string(const IdxTy flags=0) const
{
SsTy ss;
switch (flags)
{
case 1: {
if (m_ok) {ss<<i()<<" "<<j(); }
else {ss<<"-1 -1";}
break;
}
default:{if (m_ok) { ss<<"("<<i()<<","<<j()<<")"; }
else {ss<<" bad iterator "; }
}
}
return ss.str(); 
}
// these are the indexes into the sparse array
IdxTy m_row,m_col, m_rows,m_cols;
RowItor m_row_itor,m_row_end;
// these are the current dense values 
IdxTy m_i,m_j;
bool m_ok;

}; // sequential_itertor
typedef sequential_iterator SItor;

// meant to transfer this thing to the jacobian for libmesh 
template < class Tlm> void to_libmesh(Tlm & j)
{
	SItor ii(*this);
	while (ii) { j.set(ii.i(),ii.j(),*ii); ++ii; } // ii 
}

template < class Tlm> void to_paren(Tlm & j)
{
	SItor ii(*this);
	while (ii) { j(ii.i(),ii.j())=*ii; ++ii; } // ii 
}

//template < class Tlm> void to_pointer(Tlm * j, const IdxTy rows, const IdxTy order) const
// there is not fudding const itor fudd 
template < class Tlm> void to_pointer(Tlm * j, const IdxTy rows, const IdxTy order) 
{
	SItor ii(*this);
		if ( order==0){ // should work for petsc
	while (ii) { j[ii.i()+rows*ii.j()]=*ii; ++ii; } // ii 
	}
}


// user needs to clear result and insure sizes...
template < class Tv > void multiply( Tv & b, const Tv & x)
{
	SItor ii(*this);
	while (ii) 
	{
		b(ii.i())+=*ii*x(ii.j());
		++ii;
	}

}
// user needs to clear result and insure sizes...
// multiply by vector of all "1"'s
template < class Tv > void collapse( Tv & b)
{
	SItor ii(*this);
	while (ii) { b(ii.i())+=*ii; ++ii; }
}


void catalog()
{
const Myt & x=*this;
//std::map<StrTy, std::vector<SItor::LocTy> > values;
std::map<StrTy, std::vector<std::vector<IdxTy> > > values;
std::map<StrTy,StrTy> names;
const StrTy base="v";
IdxTy cnt=0;
	SItor ii(*this);
while (ii)
{
SsTy ss;
ss<<*ii;
if (values.find(ss.str())==values.end())
{
SsTy sn;
sn<<base<<cnt;
++cnt;
names[ss.str()]=sn.str();
}
values[ss.str()].push_back(ii.loc());
    ++ii;
} // itor

// std::cout<<generate(values,names)<<CRLF;


} //catalog


/*
template <class Tn, class Tv>
StrTy generate(const Tv & values,Tn & names)
{
const IdxTy sz=16; // m_dims[m_top];
const IdxTy szdim=16; // m_top;
const Myt & x=*this;
const StrTy vtype="ValTy";
const StrTy casttype="ValTy";
std::stringstream ss;
for (auto ii=names.begin(); ii!=names.end(); ++ii)
{
auto fudd_const = values.find((*ii).first);
const IdxTy fudd=(*fudd_const).second.size();
ss<<" const "<<vtype<<" "<<(*ii).second<<" = "
<<casttype<<"("<<(*ii).first<<"); // n="<<fudd<<CRLF; }
MM_ERR("// dims "<<x.m_dims.size()<<" "<<x.m_top<<" "<<x.m_dims[0]);
LoopItor2 itor(x.m_dims);
SItor ii(*this);
while (ii)
{
for (IdxTy j=0; j<sz; ++j)
{
ss<<"mat(";
for (IdxTy i=0; i<szdim; ++i) { ss<<itor[i];  ss<<",";  }
SsTy sn;
sn<<x(itor.cursor());
const StrTy varname=names[sn.str()];
ss<<itor[x.m_top]<<")="<<varname<<"; ";
++itor;
} // j
ss<<CRLF;
} //itor 
return ss.str();
} // generate
*/

template < class Os,class Ts > void dump(Os & os, const IdxTy flags, const Ts & label)
{
	const StrTy sep=" ";
	IdxTy serial=0;
	SItor ii(*this);
	IdxTy ilast=~0;
	D summ=D();
	const IdxTy inone=~0;
	SsTy ss; // this makes ctrl-C easier to at least get whole lines
	// but here just buffer whole array, dumb but easy to code	
	// this avoid dealing with real io but could use memory 

	switch (flags)
{
case 6:
case 7:
case 1:{
	while (ii) { 
		const IdxTy j=ii.j();	
		const IdxTy i=ii.i();	

		if (i!=ilast) 
		{ if (ilast!=inone) {  if (flags==7) {ss<<" summ="<<summ; summ=D(); } ss<<CRLF;  } 
				ss<<label<<sep<<serial<<sep<<" i="<<i; ilast=i;
		}
		const D & val=(*ii); 
		if (flags==1) {ss<<sep<<"(,"<<j<<")="<<*ii;}
		else if (flags==6) {ss<<sep<<"(,"<<(int(j)-int(i))<<")="<<*ii;}
		else if ((flags==7)&&(val!=0))  {summ+=(*ii); ss<<sep<<"(,"<<(int(j)-int(i))<<")="<<*ii;}

		++serial; ++ii; 
	}  
	if (flags==1) ss<<CRLF;
	if (flags==6) ss<<CRLF;
	if (flags==7) {ss<<" summ="<<summ;  ss<<CRLF;  } 

		break; }
// min and max only 
case 5:{
	D mmax=0;
	D mmin=0;
	IdxTy imin=~0; IdxTy imax=~0;
	IdxTy jmin=~0; IdxTy jmax=~0;

	while (ii) { 
		const IdxTy j=ii.j();	
		const IdxTy i=ii.i();	
	 	const D & v=*ii;
		if ((v>mmax)||(imax==(~0))) { mmax=v; imax=i; jmax=j; }	
		if ((v<mmin)||(imin==(~0))) { mmin=v; imin=i; jmin=j; }	
		++serial; ++ii; 
	}  
		ss<<label<<sep<<" min "<<mmin<<"("<<imin<<","<<jmin<<")";
		ss<<sep<<" max "<<mmax<<"("<<imax<<"."<<jmax<<")"<<CRLF;

		break; }





default :
	while (ii) { 

 ss<<label<<sep<<ii.to_string(1)<<" "<<*ii <<CRLF; 

		++serial; ++ii; 
	}  


} // switch 


	os<<ss.str();
}

template < class Os,class Tv > void dump_terms(Os & os, const IdxTy flags, const Tv & vector, const StrTy & label)
{
	const StrTy sep=" ";
	IdxTy serial=0;
	const IdxTy sz=vector.size();
	SItor ii(*this);
	while (ii) 
	{ 
		os<<label<<sep<<serial<<sep;
		const IdxTy j=ii.j();	
		if (j<sz) 
			os<<ii.i()<<sep<<j<<sep<<*ii<<sep<<vector(j)<<" x "<<(*ii*vector(j))<<" dof "<<j<<CRLF;
		else 
			os<<ii.i()<<sep<<j<<sep<<*ii<<CRLF;
		++serial; ++ii; 
	}  

}

//Myt & operator=(const Myt & that) { this->from(that); return * this; } 
// copy from 
/*void from( const Myt& that)
{
 	m_dims=that.m_dims;
	m_sizes=elements(m_dims);
	m_top=(m_dims.size()-1);
	realloc(that.m_ptr);
}
*/
/*
template <class Td> void copy_to(Td & d)
{
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) d(i)=m_ptr[i];
}
 LocTy m_sizes;
IdxTy m_top;
D * m_ptr;

typedef std::vector<LocTy> FormStack;
FormStack m_form_stack;
*/


}; // mjm_rational
//template <class Td> 
std::ostream&  operator<<(std::ostream& os,  typename mjm_sparse_matrix<double>::Row const & r)
{ os<<r.string(); return os; }



/*
#ifdef TEST_OUTER_MATRIX_MAIN__

int main(int argc, char ** argv)
{
typedef mjm_block_matrix<double> T;
//typedef T::location_type L;
const T::D  d1[]={1,2,3,4};

const T::D  d2[]={5,6,7,8};
T m1(d1,4), m2(d2,4);
MM_ERR(m1.to_string())
MM_ERR(m2.to_string())

MM_ERR("location")
MM_ERR((m1.outer_product(m2)).to_string())
return 0;
}  //main

#endif // TEST_OUTER_MATRIX_MAIN__


#ifdef TEST_BLOCK_MATRIX_MAIN__


int main(int argc, char ** argv)
{
typedef mjm_block_matrix<double> T;
typedef T::location_type L;
const T::IdxTy  szi[]={4,5,6,7};
const T::IdxTy  els=sizeof(szi)/sizeof(T::IdxTy);
MM_ERR("location")
L s=T::loc(szi,els);
MM_ERR("ctor")
T::IdxTy i=0;
T x(s);
MM_ERR("dtor")
L cursor(els);
do { x(cursor)=i; ++i; } while (x.inc(cursor,1)==0);
i=0;
L cursor2(els);
do { MM_ERR(i<<" "<<x(cursor2))  ++i; } while (x.inc(cursor2,1)==0);
x(1,2,3,4)=1234;
MM_ERR(" test "<<x(0,0,0,0)<<" "<<x(3,4,5,6)<<" "<<x(1,2,3,4))


return 0;
}


#endif // main

*/

#ifdef TEST_SPARSE_MATRIX_MAIN__
int main(int argc, char ** argv)
{
typedef mjm_sparse_matrix<double> T;
T::Row x;
x(99)=2;
x(10)=3;
x(100)=4;
x(30)=5;
x(10)=6;
MM_MSG(" output should be 10->6, 30->5, 99->2, 100->4")
x.dump(std::cout,1);


T y;
y(1,2)=3;
y(4,5)=6;
y(4,10)=610;
y(4,100)=6100;
y(4,50)=650;
y(7,8)=y(1,2)+y(4,5);
y(99999,1234)=1.234;
MM_MSG(" dumping y")
y.dump(std::cout,1);
T::sequential_iterator ii(y);

while (ii)
{
MM_MSG(" iterator r "<<ii.row()<<" c "<<ii.col()<<" i "<<ii.i()<<" j "<<ii.j())
MM_MSG(ii.ok())
MM_MSG("value ptr "<<(*ii))

++ii;
}

return 0;
}

#endif // test sparse

#endif // guard

