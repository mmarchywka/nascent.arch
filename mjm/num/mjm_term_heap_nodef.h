
template <class Ty>
class term_heap
{
typedef Ty Term;
typedef std::vector<IdxTy> Hvec;
//typedef std::map<IdxTy, Hvec> Hmap;
// user neds to include this header shot 
typedef std::unordered_map<IdxTy, Hvec> Hmap;
public:
IdxTy bad() { return ~0; } 
term_heap(){ Init(); }
~term_heap(){ Release(); }
// I guess this works ok 
const Term & operator()(const IdxTy n) { return m_heap[n]; }
// the look up returns references, these are not
// valid after a heap move doh,
bool refs_safe() const { return !m_safe_to_move; } 
void safe_to_move(const bool tf)  { m_safe_to_move=tf; } 
IdxTy chunk() const { return m_chunk;}
void chunk(const IdxTy c )  {  m_chunk=c;}

IdxTy find(const Term & t) { return Find(t); } 
StrTy dump(const IdxTy flags) { return  Dump( flags) ; } 
private:
IdxTy Find(const Term & t)
{
// index later
//for(IdxTy i=0; i<m_used; ++i) if (m_heap[i]==t) return i; 
IdxTy loc=FindOnly(t);
if (loc!=bad()) return loc;
return Alloc(t);
}
IdxTy Alloc(const Term & t)
{
// danger will robinson - the () look up returns references
// which are addresses, all the expands screw that up ... 
// the hysterisis was turned off for testing but now appears to work
if (m_used>(m_sz>>1) ) 
Expand(m_sz+m_chunk);
if (m_used==m_sz)
{
MM_ERR(" term heap filled up will bomb now "<<MMPR3(m_used,m_sz,t.dump(0,0)))
}
m_heap[m_used]=t;
const IdxTy u=m_used;
auto & vv =m_index[t.hash()];
vv.push_back(u);
++m_used;
return u; 
}



IdxTy FindOnly(const Term & t)
{
// index later
t.rehash();
const auto & vv =m_index[t.hash()];
//for(IdxTy i=0; i<m_used; ++i) if (m_heap[i]==t) return i; 
for(IdxTy i=0; i<vv.size(); ++i) if (m_heap[vv[i]]==t) return vv[i]; 
return bad();
}

void Release()
{
delete [] m_heap;
m_sz=0;
m_used=0;

}
void New(const IdxTy n)
{
Release();
m_heap=new Term[n];
m_sz=n;
m_used=0;
m_index.clear();
}
void Expand(const IdxTy n)
{
if (!m_safe_to_move) return; 
if (n<m_sz) return;
Term * heap = new Term[n];
// in theory this should work for m_heap==0 
for(IdxTy i=0; i<m_used; ++i ) heap[i]=m_heap[i];
m_sz=n;
delete[] m_heap;
m_heap=heap;


}


void Init() 
{
m_heap=0;
m_sz=0;
m_used=0;
m_chunk=2048;
m_safe_to_move=true;
New(m_chunk);
}
StrTy Dump(const IdxTy flags)
{
Ss ss;
MM_LOOP(ii,m_index)
{
ss<<(*ii).first;
const auto & v=(*ii).second;
MM_LOOP(jj,v) ss<<" "<<(*jj);
ss<<CRLF;
}
return ss.str();
}
Term * m_heap;
IdxTy m_sz;
IdxTy m_used;
IdxTy m_chunk;
Hmap m_index;
bool m_safe_to_move;

}; // term_heap
