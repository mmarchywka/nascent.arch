#ifndef MJM_INDEXING_H__
#define MJM_INDEXING_H__

namespace indexing {
class index_traits
{
public:
typedef unsigned int IdxTy;


}; 

template < class Tr=index_traits> class basic_iterator
{
typedef basic_iterator Myt;
typedef Tr::IdxTy IdxTy;

IdxTy m_atom,m_coord,m_terminal_state,m_idx;
public:
basic_iterator(const IdxTy sz)
:m_atom(0),m_coord(0), m_termina_state(sz),m_idx(0){}

void clear() { m_atom=0; m_coord=0; m_idx=0; }
Myt & operator++() {Inc(); return *this; }
operator bool() const { return m_idx<m_terminal_state; }
IdxTy operator *() const { return m_idx; }
operator IdxTy() const { return m_idx; }
IdxTy atom() const { return m_atom; }
IdxTy coord() const { return m_coord; }
void atom_coord(IdxTy &a , IdxTy & c) const { a=m_atom; c=m_coord;}
IdxTy idx() const {return m_idx;}
private:
void Inc()
{
++m_coord; if ( m_coord==3) { m_coord=0; ++m_atom;}
++m_idx;

}
};


}; // ns

#endif

