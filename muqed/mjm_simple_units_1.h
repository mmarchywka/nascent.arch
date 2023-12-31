#ifndef MJM_SIMPLE_UNITS_H__
#define MJM_SIMPLE_UNITS_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

// should just be ragged table but allows more saving. 
//#include "mjm_ragged_forms.h"
#include "mjm_collections.h"
// don't use this yet 
//#include "mjm_string_tokenizer.h"
#include "mjm_read_buffer.h"
#include "mjm_indexed_map.h"
#include "mjm_unit_crap.h"
#include "../toobib/mjm_vv_map.h"
#include "../toobib/mjm_wovdb.h"
//#include "mjm_constrained_file.h"

#include <map>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <stdexcept>


// Sun Sep  5 09:38:58 EDT 2021
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_simple_units   
// g++  -Wall -std=gnu++11 -DTEST_MJM_SIMPLE_UNITS -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_simple_units.h  -lpthread -lreadline


template <class Tr>
class mjm_cluster_map 
{
 typedef mjm_cluster_map Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef std::map<StrTy, IdxTy> Cmap;
typedef std::vector<StrTy> Cvec;
typedef std::map<IdxTy, Cvec> Clist;
public:
enum { BAD=~0};
mjm_cluster_map(): m_next(1) {}
IdxTy add(const StrTy& x, const StrTy & y) { return Add(x,y); }
IdxTy lookup(const StrTy & name) const { return Lookup(name); }
//  the cluster needs to exist doh...  
const Cvec & cluster(const IdxTy i) const { return (*(m_list.find(i))).second; }  
IdxTy bad() const { return BAD; } 
StrTy dump(const IdxTy flags=0) const  { return Dump(flags); } 
private:
StrTy Dump(const IdxTy flags=0) const 
{
Ss ss;
ss<<MMPR(m_next)<<CRLF;
//Cmap m_map;
MM_LOOP(ii,m_map) 
{ ss<<" cluster m_map "<<(*ii).first<<" "<<(*ii).second<<CRLF; }
//Clist m_list;
MM_LOOP(ii,m_list) 
{
ss<<(*ii).first;
const auto & xx=(*ii).second;
ss<<" cluster list ";
MM_LOOP(jj,xx) { ss<<" "<<(*jj); } 
ss<<CRLF;
} // ii 
//IdxTy m_next;

return ss.str();
} // Dump
IdxTy Add(const StrTy& x, const StrTy & y) { 
const IdxTy lx=Lookup(x);
const IdxTy ly=Lookup(y);
const bool xbad=(lx==IdxTy(BAD));
const bool ybad=(ly==IdxTy(BAD));
const IdxTy rc=(xbad?2:0)|(ybad?1:0); 
MM_ERR(MMPR4(lx,ly,xbad,ybad))
if (!xbad&&!ybad)
{
 if (lx==ly) return 8; 
Merge(lx,ly);
return rc;
}
if (xbad&&ybad)
{
const IdxTy newx=Seed(x);
Add(y,newx);
return rc;
}
if (ybad) { Add(y,lx); return rc; } 
Add(x,ly); 
return rc; 
} // Add
// move y into x 
IdxTy Merge(const IdxTy lx,const IdxTy ly)
{
// refs ok until merge lol. 
auto fx=m_list[lx];
const auto fy=m_list[ly];
MM_LOOP(ii,fy)
{
m_map[(*ii)]=lx;
fx.push_back((*ii));
} // ii 

m_list[lx]=fx;
// could delete from map, but 
//m_list[ly].clear();
m_list.erase(m_list.find(ly));
return 0;
} // Merge

IdxTy Add(const StrTy & x , const IdxTy xc) 
{
m_map[x]=xc;
m_list[xc].push_back(x);
return 0; 
} // Add
IdxTy Seed(const StrTy & x)
{
const IdxTy xc=m_next;
++m_next;
m_map[x]=xc;
m_list[xc].push_back(x);
return xc;
}
IdxTy Lookup(const StrTy & name) const { 
auto ii=m_map.find(name);
if (ii==m_map.end()) return BAD;
return (*ii).second; 
} // Lookup

Cmap m_map;
Clist m_list;
IdxTy m_next;
}; // mjm_cluster_map

template <class Tr>
class mjm_simple_units 
{
 typedef mjm_simple_units Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef string_tokenizer St;
typedef St::coded_type Co;
typedef mjm_read_buffer<Tr>  Rb;

typedef D F;
typedef mjm_unit_crap<Tr> UnitCrap;
typedef typename UnitCrap::unit_desc UnitDesc;
class _value_class
{

public:
// should be one... 
_value_class(): m_value(0) {}
_value_class( const StrTy &  unit, const StrTy & thing,  const F & f)
:m_unit(unit), m_thing(thing), m_value(f) {}

const StrTy & u() const { return m_unit; } 
const StrTy & t() const { return m_thing; } 
const F& v() const { return m_value; } 

StrTy dump() const
{
Ss ss;
ss<<MMPR3(m_unit,m_thing,m_value);

return ss.str(); 
} // dump

private:
StrTy m_unit,m_thing;
F m_value;

}; // _value_class

typedef _value_class ValueClass;

class _conversion_entry
{
typedef _conversion_entry Myt;
public:
_conversion_entry() : m_factor(0) {}
_conversion_entry( const StrTy &  from_unit, const StrTy & from_thing,
const StrTy & to_unit, const StrTy & to_thing, const F & f)
:m_from_unit(from_unit), m_from_thing(from_thing),
m_to_unit(to_unit), m_to_thing(to_thing),m_factor(f) {}
_conversion_entry( const StrTy &  from_unit,
const StrTy & to_unit,  const F & f)
:m_from_unit(from_unit),
m_to_unit(to_unit), m_factor(f) {}

const StrTy& from_unit() const { return m_from_unit;}
const StrTy& to_unit() const { return m_to_unit;}
const StrTy& to_thing() const { return m_to_thing;}
const StrTy& from_thing() const { return m_from_thing;}
ValueClass operator*(const ValueClass & that) const { return Mult(that); } 
Myt operator*(const Myt & that)const { return Mult(that); } 
Myt operator/(const Myt & that)const  { return Div(that); } 
Myt flip() const { return Flip(); } 
const bool generic() const 
{ return (m_from_thing.c_str()[0]|m_to_thing.c_str()[0])==0; }
StrTy dump() const
{
Ss ss;
ss<<MMPR4(m_to_unit,m_to_thing,m_from_unit,m_from_thing);
ss<<MMPR(m_factor); 
return ss.str();
} 
private:


// if the conversion factor lacks a thing, any is ok...  
bool Ok(const StrTy & tv, const StrTy &tc) const
{
if (tc=="") return true;
if (tc!=tv) return false;
return true; 
}
bool OkBi(const StrTy & tv, const StrTy &tc) const
{
if (tc=="") return true;
if (tv=="") return true;
if (tc!=tv) return false;
return true; 
}


// the value has to cancel the denominator or "from"
// with identical units and thing must be either lacking
// in factor or identical  
bool  Can_mult(const ValueClass & vc) const
{
const bool can_mult=Ok(vc.t(),m_from_thing);
return can_mult&&(vc.u()==m_from_unit);
} 
// one or the other set needs to cancel...
bool  Can_mult1(const Myt & that) const
{
//const bool can_mult=Ok(vc.t(),m_from_thing);
//return can_mult&&(vc.u()==m_from_unit);
bool can_mult1=OkBi(that.m_to_thing,m_from_thing);
 can_mult1&=(that.m_to_unit==m_from_unit);
return can_mult1;
} 
bool  Can_div1(const Myt & that) const
{
//const bool can_mult=Ok(vc.t(),m_from_thing);
//return can_mult&&(vc.u()==m_from_unit);
bool can_div1=OkBi(that.m_from_thing,m_from_thing);
 can_div1&=(that.m_from_unit==m_from_unit);
return can_div1;
} 


bool  Can_div(const Myt & that) const
{ return Can_div1(that)||that.Can_div1(*this); } 


bool  Can_mult(const Myt & that) const
{ return Can_mult1(that)||that.Can_mult1(*this); } 


bool  Can_div(const ValueClass & vc) const
{
const bool can_div=Ok(vc.t(),m_to_thing);
return can_div&&(vc.u==m_to_unit);
} 


IdxTy can(const ValueClass & vc) const
{
// if our "thing" is blank it goes with all 
// otherwise that needs to match... 
const bool can_mult=Can_mult(vc); // Ok(vc.t(),m_from_thing);
const bool can_div=Ok(vc.t(),m_to_thing);
// if the thing part is right, see if the units 
// are identical... 

return 0; 
}
Myt Flip( ) const
{
Myt x;
x.m_to_thing=m_from_thing;
x.m_from_thing=m_to_thing;
x.m_to_unit=m_from_unit;
x.m_from_unit=m_to_unit;
x.m_factor=1.0/m_factor;
return x; 
} // Flip 
Myt Div( const Myt & that) const
{
Myt x=that.Flip();
return Mult(x);

}
Myt Mult( const Myt & that) const
{
Myt x;
// can_mult1&=(that.m_to_unit==m_from_unit);
if (Can_mult1(that))
{
x.m_factor=m_factor*that.m_factor;
x.m_to_unit=m_to_unit;
x.m_from_unit=that.m_from_unit;
x.m_to_thing=(m_to_thing.length()==0)?that.m_to_thing:m_to_thing;
x.m_from_thing=(m_from_thing.length()==0)?that.m_from_thing:m_from_thing;
}
else if (that.Can_mult1(*this))
{
return that.Mult(*this);
}
else  throw std::logic_error(Dump(that));
return x;
} // operator*

ValueClass Mult(const ValueClass & x) const
{
if (!Can_mult(x)) throw std::logic_error(Dump(x));

ValueClass xr(m_to_unit,m_to_thing,m_factor*x.v());

return xr; 

}
StrTy Dump(const ValueClass & x) const
{
Ss ss;
ss<<" incompdatible "<<MMPR2(Dump(), x.dump());
return ss.str();
}
StrTy Dump(const Myt & x) const
{
Ss ss;
ss<<" incompdatible "<<MMPR2(Dump(), x.Dump());
return ss.str();
}


StrTy Dump() const
{
Ss ss;
ss<<MMPR4( m_from_unit, m_from_thing,  m_to_unit, m_to_thing)
<<MMPR( m_factor);
return ss.str();
}
StrTy m_from_unit, m_from_thing;
StrTy m_to_unit, m_to_thing;
F m_factor;

}; // _conversion_entry

typedef _conversion_entry ConvEntry;

typedef mjm_wovdb<Tr,ConvEntry> ConvDb;
typedef typename ConvDb::vector_type DbHits;
typedef std::vector<ConvEntry> ConvVec;
typedef std::map<StrTy,ConvVec> ConvMap;
typedef std::map<StrTy,StrTy> TypeMap;
typedef std::map<StrTy,ConvEntry> BaseEntries;

typedef mjm_cluster_map<Tr> Clusters;

public:
mjm_simple_units() {Init(); }
~mjm_simple_units() {}
IdxTy load(const StrTy & fn) { m_defs.load(fn); return size();}
IdxTy size() const { return m_defs.size(); }
template <class Tv> 
IdxTy convert(  F & vout, const F & vin, const StrTy & to_unit, 
const StrTy & to_thing,
const StrTy & from_unit, const StrTy & from_thing, 
const Tv & vspec,
const IdxTy flags)
{
return Convert(vout,vin,to_unit,to_thing,from_unit,from_thing,flags);
} 
IdxTy convert(  F & vout, const F & vin, const StrTy & to_unit, 
const StrTy & to_thing,
const StrTy & from_unit, const StrTy & from_thing, 
const IdxTy flags)
{
return Convert(vout,vin,to_unit,to_thing,from_unit,from_thing,ConvVec(),flags);
} 
IdxTy convert(  F & vout, const F & vin, const StrTy & to_unit, 
const StrTy & from_unit, const IdxTy flags)
{
return convert(vout,vin,to_unit,StrTy(),from_unit,StrTy(),flags);
} 


IdxTy parse(const IdxTy flags) {return  Parse(flags); } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  

//Ragged m_defs; IdxTy m_istart;

//ConvDb m_db;
for(IdxTy i=0; i<m_db.size(); ++i)
{
const auto & x=m_db[i];
ss<<" db entry "<<MMPR(i)<<x.dump()<<CRLF;
} 
//Clusters m_clust;
ss<<m_clust.dump();
//BaseEntries m_bases;
MM_LOOP(ii,m_bases)
{
//MM_ERR(MMPR2((*ii).first,(*ii).second.dump()))
ss<<" base facotr "<< MMPR2((*ii).first,(*ii).second.dump())<<CRLF;
} // ii 

return ss.str(); 
}
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

IdxTy Parse(const IdxTy flags) {
// the unit crap contains a tokenizer crap 
const bool quad_parse=Bit(flags,0); 
St st;
Ragged & r= m_defs;
IdxTy lno=0;
const IdxTy rsz=r.size(); // MM_LOOP(ii,r)
for(IdxTy i=m_istart; i<rsz; ++i)
{
const Line & l =r[i]; // (*ii);
if (quad_parse)
{

++lno;
continue; 
}  // quad_parse
if (l.size()==0)  { ++lno; continue; }
if (l[0].c_str()[0]=='#') { ++lno; continue; }
//MM_ERR(" parsing line "<<MMPR2(lno,l.size()))
UnitDesc x(l,st);
if (x.error())
{
// if length is zero no error possible
if (l.size()>0)
MM_ERR(" bad units abbrv "<<MMPR3(lno,l[0],x.error_string()))
else
MM_ERR(" error with zero length line wtf")
}
//m_map[m_st(x.name())]=x;
//m_map.add_multi(m_st(x.name()),x.labels(),x);
// find the type ( length mass etc ) and get all the aliases too
// 
//typedef _conversion_entry ConvEntry;
//ConvEntry ce(from,to,v);
// get name too, all using st() 
//const Aliases & aka() const { return m_aliases; }
// use unity conversion for aka 
//const ConvMap & conv_factors() const { return m_conv; }
const StrTy nm=st(x.name());
const auto & a=x.aka();
MM_LOOP(jj,a) { ConvEntry ce(nm,((*jj).first),1); Add(ce,0); } // jj 
const auto & b=x.conv_factors();
MM_LOOP(jj,b) 
{ ConvEntry ce(nm,st((*jj).first),(*jj).second); Add(ce,0); } // jj 

++lno;
} // ii 
m_istart=rsz;
return 0;

} // Parse
template <class Tv>
IdxTy Convert(  F & vout, const F & vin, 
const StrTy & to_unit, const StrTy & to_thing,
const StrTy & from_unit, const StrTy & from_thing,
const Tv & vspec
, const IdxTy flags)
{
IdxTy rc=0; 
MM_ERR(MMPR4(vin,to_unit,to_thing,from_unit)<<MMPR(from_thing))
if (vin==0) { vout=0; return rc; }
const bool same_thing=(to_thing==from_thing);
//const bool no_thing=same_thing&&(to_thing.length()==0);
const bool same_unit=(to_unit==from_unit); 
if (same_unit) if (same_thing) { vout=vin; return 0; } 
const IdxTy cto=m_clust.lookup(to_unit);
const IdxTy cfrom=m_clust.lookup(from_unit);
const bool good=(cto!=m_clust.bad())&&(cfrom!=m_clust.bad());
const bool have_both=good&&(cto==cfrom);
// see if these things are in the DB and can be reconciled
//const DbHits * inf=m_db.find("names",from_unit);
//const DbHits * outf=m_db.find("names",to_unit);
// if any contain both, have a direct conversion else
// need to find intermediates or customs
//const bool have_inf=(inf!=0)&&(inf->size());
//const bool have_outf=(outf!=0)&&(outf->size());
// TODO FIXME in theory custom factors could override fixed but
// ignore that for now  
// go back to base units.
MM_ERR(MMPR4(cto,cfrom,good,have_both))
if (have_both&&same_thing) //if simple 
{
auto ii=m_bases.find(to_unit);
if (ii==m_bases.end()) throw std::logic_error(to_unit); 
const ConvEntry & b1=(*ii).second; // this converts to to base
ii=m_bases.find(from_unit);
if (ii==m_bases.end()) throw std::logic_error(from_unit); 
const ConvEntry & b2=(*ii).second; // this converts from to base

ValueClass in(from_unit,"",vin);
ConvEntry b1flip=b1.flip();
MM_ERR(MMPR3(b1.dump(),b2.dump(),b1flip.dump()))
ValueClass out=b1flip*(b2*in);
vout=out.v();
if (out.u()!=to_unit)
{
Ss ss; 
ss<<" result is wrong "<<MMPR(vin)<< MMPR4(vout,from_unit,to_unit,out.u());
MM_ERR(ss.str())
} // !=
return 0;
} // simple 
///MM_ERR(" copmildated "<<MMPR(kkkkkkkkkkkkk 
// If not check the custom factors.... 

return 0;
} // Convert
IdxTy Add(const ConvEntry & ce, const IdxTy flags)
{
MM_ERR(" ADDING "<<MMPR(ce.dump()))
if (!ce.generic())
{
MM_ERR(" not adding specific "<<MMPR(ce.dump()))
return ~0;
}
const IdxTy loc=m_db.add(ce);
m_db.index(loc,"names",ce.from_unit());
m_db.index(loc,"names",ce.to_unit());

// cluster needs to add a made up base unit conversion...
// for a cluster of just two, no new one is needed
// otherwise, need to convert the new unit to the first one
// and add it 
// these are in the same cluster 
// either they alaready were, or a new one was added, or a merge
const IdxTy cadd=m_clust.add(ce.from_unit(),ce.to_unit());
// find base entry of cluster 
// both in same cluster, 
if (8==cadd) return 0;  // non zero means it had to cluster
//{
// find the members of the cluster.. 
const IdxTy cn=m_clust.lookup(ce.from_unit());
const IdxTy cntoo=m_clust.lookup(ce.to_unit());
if (cn!=cntoo){
Ss ss; ss<<MM_MARK<<" bad clusters after add "<<MMPR4(cadd,cn,cntoo,ce.dump()); 
 throw std::logic_error(ss.str()); }
const auto & v=m_clust.cluster(cn);
const IdxTy sz=v.size();
MM_ERR(MMPR3(sz,cn,cntoo))
const StrTy & cebase=v[0]; // should identity
IdxTy i=sz;
if (sz==2)
{ // new cluster just seed it... 
// should added from first ... 
MM_ERR(MMPR2(ce.dump(),cebase))
if (cebase==ce.from_unit())
{
ConvEntry bc(ce.from_unit(),ce.from_unit(),1);
m_bases[ce.from_unit()]=bc;
m_bases[ce.to_unit()]=ce.flip();
}
else
if (cebase==ce.to_unit())
{
ConvEntry bc(ce.to_unit(),ce.to_unit(),1);
m_bases[ce.to_unit()]=bc;
m_bases[ce.from_unit()]=ce; // .flip(); // .flip();
}
else throw std::logic_error(" wtf ");
i=0;
return 0;
} // sz==2
// either a merge occured between existin units
// or only one of them was new and it was added to the end... 
auto jj=m_bases.find(ce.from_unit());
auto kk=m_bases.find(ce.to_unit());
const bool from_missing=(jj==m_bases.end());
const bool to_missing=(kk==m_bases.end());
if (!from_missing&&!to_missing)
{ // merged, need to fix old ones.. 
// the old ones should be at the end ( need to add this
// to clustal class for consistency )
//const ConvEntry  clustbase
const ConvEntry b1=(*jj).second;
const ConvEntry b2=(*kk).second;
// vould both come back wrong? 
const bool b1isok=(b1.to_unit()==cebase); //from unit convert ok
const bool b2isok=(b2.to_unit()==cebase);
MM_ERR(MMPR4(cadd,b1isok,b2isok,ce.dump()))
MM_ERR(MMPR2(b1.dump(),b2.dump()))
if (!b1isok&&!b2isok)
{
MM_ERR(" danger will robinson neither is ok ")

}
// from unitis ok, conv base/from but to is oldbase/to
const ConvEntry old_to_new=  b1isok?   //  ce*b1/b2;
			b1/(ce*b2): // (base/from)/((to/from)*(oldbase/to) ) 
			b2/(ce*b1);
MM_ERR(MMPR(old_to_new.dump()))
while (i) // not a new but added or merged... 
{
--i;
// the last one needs a new entry 
const StrTy & cei=v[i];
// this should convert to cebase
// if there was no entyr, its a new factor and the only
// conversion is the current one which should already have
// a conversion to the base... 
// it had to already be there or this was not a merge... 
auto ii=m_bases.find(cei);
if (ii==m_bases.end())
{
Ss ss; ss<<MM_MARK<<" cluster wrong "<<MMPR3(cadd,cei,ce.dump());
throw std::logic_error(ss.str()); 
} // end 
// the entry exists, see if it needs to be modified.. 
ConvEntry old=(*ii).second; // m_bases[cei];
const StrTy ou=old.to_unit(); 
if (ou==cebase) break; // that should be the last one..  
// the old to_unit is wrong and this is part of a merged
// cluster.  
//ConvEntry bc(v,to,from);
MM_ERR(MMPR3(i,old.dump(),old_to_new.dump()))
ConvEntry bc=old*old_to_new ;   // (cebase,cei,v);
m_bases[cei]=bc;
} // while i 


return 0; 
} // neither missing, merfed
if (from_missing&&to_missing)
{
Ss ss;
ss<<" cluster failure "<<MM_MARK<<MMPR(ce.dump());
throw std::logic_error(ss.str()); 
}
if (from_missing)
{ // make conversion to base using current factor and other base
const ConvEntry x1=m_bases[ce.to_unit()];
m_bases[ce.from_unit()]=x1*ce;
}
else if (to_missing)
{ // make conversion to base using current factor and other base
const ConvEntry x1=m_bases[ce.from_unit()];
m_bases[ce.to_unit()]=x1/ce;
}
// 
//} // cadd

return 0; 
} // Add

void Init()
{
m_istart=0; 
} // Init 

// MEMBERS
Ragged m_defs;
IdxTy m_istart;
//ConvVec m_cons;
//ConvMap m_map;
//TypeMap m_types;
ConvDb m_db;
Clusters m_clust;
BaseEntries m_bases;


}; // mjm_simple_units

//////////////////////////////////////////////

template <class Tr>
class mjm_simple_units_map : public std::map<typename Tr::StrTy, mjm_simple_units< Tr > >  
{
 typedef mjm_simple_units_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_simple_units< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_simple_units_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }

private:

void Init()
{
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;


}
return ss.str();
// return Dump(flags); 

}




private:

}; // mjm_simple_units_map




////////////////////////////////////////////
#ifdef  TEST_MJM_SIMPLE_UNITS
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_simple_units <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_SIMPLE_UNITS "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_simple_units<Tr>  Myt;
//Myt x(argc,args);
Myt x;

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="convert") {
typedef double D;
D vin,vout;
vin=atof(cip.p1.c_str());
StrTy to=cip.p2;
StrTy from=cip.wif(3);
const IdxTy rc=x.convert(vout,vin,to,from,0);
MM_ERR(MMPR4(vout,vin,to,from)<<MMPR(rc))
} // convert
// 
//IdxTy convert(  F & vout, const F & vin, const StrTy & to_unit, 
//const StrTy & from_unit, const IdxTy flags)
else if (cmd=="load") { x.load(cip.p1); }
else if (cmd=="parse") { x.parse(0); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_SIMPLE_UNITS_H__ 
