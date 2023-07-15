#ifndef MJM_TABLE_MANAGERS_H__
#define MJM_TABLE_MANAGERS_H__

#include "mjm_globals.h"
#include "mjm_instruments.h"
//#include "mjm_generic_iterators.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*

g++ -DTEST_ITERATOR_TRI_MAIN__ -O0 -gdwarf-3 -std=c++11 -Wall -I.. -x c++ mjm_iterator_base.h 


g++ -DTEST_ITERATOR_BASE_MAIN__ -O0 -gdwarf-3 -std=c++11 -Wall -I.. -x c++ mjm_iterator_base.h  


 g++ -DTEST_PART_ITERATORS_MAIN__ -std=c++11 -Wall -I.. -x c++ mjm_part_iterators.h


*/

// for now go with simple and virtuals 
//template <class Timp> class mjm_iterator_base  : public Timp 
template <class Ty>  class mjm_table_base  // : public Timp 
{
typedef mjm_table_base Myt;
//typedef Timp Super ;
public:
class Tr {
public:
//typedef Ty D;
typedef unsigned int IdxTy;
typedef unsigned int MUInt;
typedef double  D;
typedef double  MUD;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
//typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocTy;
}; // Tr

typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::MUInt MUIntTy;
typedef typename Tr::MUD MUDTy;
typedef typename Tr::D D;
typedef typename Tr::Sint Sint;
typedef typename Tr::StrTy StrTy;
typedef typename Tr::SsTy SsTy;
typedef typename Tr::SsTy Ss;
typedef typename Tr::LocTy LocTy;

typedef LocTy location_type;

typedef Ty data_type;
typedef std::vector<const data_type *> TableVector;

class TableInfo
{
public:
IdxTy sz;
}; // TableInfo
typedef TableInfo table_info_type;
typedef std::vector<table_info_type> TableInfoVector;
//typedef std::vector<LocTy> PartList;
//typedef std::vector<IdxTy> ExcludeList;
//mjm_fixed_sum_itor(const IdxTy n, const MUIntTy & sum )
//: m_sum(sum),m_n(n) {Init(); } 
mjm_table_base( )
//: m_k(0),m_n(0),m_maxmul(0) 
{Init(); } 
~mjm_table_base( ) { Free(); } 

data_type * add_table(const IdxTy i, const IdxTy sz )  
{
return AddTable(i,sz);
}
void adopt_table(const data_type * p , const IdxTy i, const IdxTy sz )  
{
return AdoptTable(p,i,sz);
}

const bool have(const IdxTy i) const 
{ 
if (m_tables.size()<=i) return false; 
return m_tables[i]!=0; }

const data_type * operator[]( const IdxTy i ) const { return table(i); } 
const data_type * table(const IdxTy i ) const {
if ( m_tables[i]==0) { MM_ERR(" no table "<<MMPR2(i,to_string())) }
return m_tables[i];

}

StrTy to_string(const IdxTy flags=0) const 
{
Ss ss;
ss<<" tables "; 
IdxTy i=0; 
MM_LOOP(ii,m_tables)  
{ if ((*ii)!=0) 
	{ss<<i<<",sz="<<m_infos[i].sz<<" "; } ++i; 
} 
return ss.str();
}


/*
const bool done() const { return m_done; }
const bool ok() const { return !m_done; }
operator bool() const { return !m_done; }
const IdxTy & operator[]( const IdxTy i) const {return m_cursor[i];} 
IdxTy & operator[]( const IdxTy i)  {return m_cursor[i];} 
// this would be better ... 
IdxTy & twiddle( const IdxTy i)  {return m_cursor[i];} 

Myt & operator++( )  {Inc(); return *this;} 
*/

private:

protected:

//virtual void Init() //  =0; 
 void Init() //  =0; 
{
}
 void Free() //  =0; 
{
MM_LOOP(ii, m_tables) { delete [] (*ii); }
m_tables.clear();
m_infos.clear();
}
data_type * AddTable(const IdxTy i, const IdxTy sz )  
{
ResizeFor(i);
delete[] m_tables[i];
m_tables[i]= new data_type[sz];
m_infos[i].sz=sz;
}
void AdoptTable(const data_type * p, const IdxTy i, const IdxTy sz )  
{
ResizeFor(i);
delete[] m_tables[i];
m_tables[i]= p; // new data_type[sz];
m_infos[i].sz=sz;
}

void ResizeFor(const IdxTy i)
{
if (m_tables.size()>i) return;
m_tables.resize(i+1);
m_infos.resize(i+1);

}

private:

TableVector m_tables;
TableInfoVector m_infos;

}; // mjm_table_base



#ifdef TEST_COUNT_RAG_MAIN__
int main(int argc, char ** argv)
{
typedef unsigned int IdxTy;
typedef std::vector<IdxTy> LocTy;
LocTy c;
for(int i=1; i<argc; ++i) c.push_back(atoi(argv[i]));
typedef mjm_ragged_card_itor Mrc;
Mrc mrc(c);
//std::vector<IdxTy> c;
//c.push_back(
std::map<IdxTy, IdxTy> cnts;
for(mrc.reset(); mrc.ok(); ++mrc)
{
//MM_MSG(mti.to_string());
++cnts[mrc.sum()];
}
MM_LOOP(ii,cnts) { MM_MSG(MMPR2((*ii).first,(*ii).second)) } 
MM_LOOP(ii,c) { MM_MSG("C "<<MMPR((*ii))) } 
return 0; 
}
#endif

#endif




