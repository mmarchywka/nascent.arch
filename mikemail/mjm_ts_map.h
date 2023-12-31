#ifndef MJM_TS_MAP_H__
#define MJM_TS_MAP_H__



#include "mjm_globals.h"
#include <mjm_thread_util.h>
#include <map>


// Sun May  5 13:27:25 EDT 2019
// generated by -classhdr mjm_ts_map  

/*

A grossly thread safe map with limited interface and limited
std::map compatibility. 


*/


template <class Tk, class Tv, class Tr>
class mjm_ts_map 
{
 typedef mjm_ts_map Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;
typedef typename std::map<Tk,Tv> Map;
typedef typename Map::iterator iterator;
typedef typename Map::const_iterator const_iterator;

typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};

public:
mjm_ts_map() { Init(); } 

const IdxTy size() const  
{
EnterSerial(MAP_MU);
IdxTy sz=m_map.size();
ExitSerial(MAP_MU);
return sz;
}
const_iterator  find(const Tk & k ) const
{
EnterSerial(MAP_MU);
const_iterator   x=m_map.find(k);
ExitSerial(MAP_MU);
return x ;
}
iterator  find(const Tk & k ) 
{
EnterSerial(MAP_MU);
iterator   x=m_map.find(k);
ExitSerial(MAP_MU);
return x ;
}


const_iterator  begin() const
{
EnterSerial(MAP_MU);
const_iterator   x=m_map.begin();
ExitSerial(MAP_MU);
return x ;
}
iterator  begin( ) 
{
EnterSerial(MAP_MU);
iterator   x=m_map.begin();
ExitSerial(MAP_MU);
return x ;
}


const_iterator  end() const
{
EnterSerial(MAP_MU);
const_iterator   x=m_map.end();
ExitSerial(MAP_MU);
return x ;
}

iterator  end( ) 
{
EnterSerial(MAP_MU);
iterator   x=m_map.end();
ExitSerial(MAP_MU);
return x ;
}


void clear( ) 
{
EnterSerial(MAP_MU);
m_map.clear();
ExitSerial(MAP_MU);
}



/*
const Tv & operator[](const Tk & k ) 
{
EnterSerial(MAP_MU);
Tv & x=m_map[k];
ExitSerial(MAP_MU);
return x; 
}
*/

Tv & operator[](const Tk & k ) 
{
EnterSerial(MAP_MU);
Tv & x=m_map[k];
ExitSerial(MAP_MU);
return x; 
}

private:

//volatile
 Map m_map;

mutable MutexVector m_mutex_vector;


void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }

void Init()
{
m_mutex_vector = MutexVector(MU_SZ);
}

}; // mjm_ts_map

#endif // MJM_TS_MAP_H__ 
