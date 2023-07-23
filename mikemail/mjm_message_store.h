#ifndef MJM_MESSAGE_STORE_H__
#define MJM_MESSAGE_STORE_H__




//# forgot I had this lol 
#if 0
template <class Tr, class Tm >
class mjm_message_junk
{
 typedef mjm_message_junk Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
typedef typename std::map<StrTy,StrTy> Map;

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;
enum { MAP_MU=0 , MU_SZ};
public:

private:

StrTy m_verbatim, m_headers, m_body;
StrTy m_flags;
// headers reconstructed from split line versions
StrTy m_concat_headers;
Map m_att;

MutexVector m_mutex_vector;

void enter_serial(const IdxTy i) {  m_mutex_vector.enter_serial(i ); }
void exit_serial(const IdxTy i) {  m_mutex_vector.exit_serial(i ); }
void Init()
{
m_mutex_vector = MutexVector(MU_SZ);

}

}; // message_junk
#endif



template <class Tr, class Tm >
class mjm_message_store 
{
 typedef mjm_message_store Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
typedef typename std::map<IdxTy, Tm> Map;

//typedef pthread_mutex_t Mutex;
//typedef typename Tr::Mutex Mutex;
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};


public:
mjm_message_store() { Init(); }

~mjm_message_store()
{
//pthread_mutex_destroy(mutex1());
//pthread_mutex_destroy(mutex2());

}
void clear() { m_map.clear(); } 
void append(const Tm & m)
{
enter_serial(MAP_MU);
IdxTy sz=m_map.size();
m_map[sz]=m;
exit_serial(MAP_MU);
}
Tm & operator[](const IdxTy n) 
{
enter_serial(MAP_MU);
// should be stable during map updates
Tm & x=m_map[n];
exit_serial(MAP_MU);
return x; 
}
const IdxTy size()
{
enter_serial(MAP_MU);
IdxTy sz=m_map.size();
exit_serial(MAP_MU);
return sz; 
}
private:

//volatile 
Map m_map;


MutexVector m_mutex_vector;
void enter_serial(const IdxTy i) {  m_mutex_vector.enter_serial(i ); }
void exit_serial(const IdxTy i) {  m_mutex_vector.exit_serial(i ); }

void Init()
{
m_mutex_vector = MutexVector(MU_SZ);

}

}; // mjm_message_store

#endif // MJM_MESSAGE_STORE_H__ 
