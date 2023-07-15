#ifndef PARAMLIST_COMPILER_CPP_H__MJM
#define PARAMLIST_COMPILER_CPP_H__MJM

 
/*

2014-04-10 try to remove junk from proof of concept version
and reorganize a bit, remove most profanity.  Remove the R code.
In c++ with STL this could be quite nice way to interface libraries
with API's likely to change


*/

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "xc.h"

// all the c code shoudl have been removed
//#include "mjm_paramlist_compiler_libxc.h"

#ifdef __cplusplus 
#include <vector>
#include <map>
#include <string>
#include <iostream>



namespace nv_api
{
enum { BAD=~0};
typedef unsigned int IdxTy;
typedef std::ostream OsTy;
typedef  int return_code;
typedef IdxTy ParamListHandle;
// yes this is still here but user need not care for now...
typedef void * default_value_type;
// the basics
typedef std::string name_type;
typedef std::string StrTy;
typedef IdxTy parameter_id_type;


typedef std::vector<name_type> BaseParamNameList;
class vector_with_c_loading;
typedef vector_with_c_loading ParamNameList;
typedef std::vector<default_value_type> BaseParamValuesArray;
class vector_with_name_checking;
typedef vector_with_name_checking ParamValuesArray;
typedef std::map<name_type,default_value_type>  BaseParamNameValuePairList;  
class map_with_voids;  
typedef map_with_voids ParamNameValuePairList;  
typedef std::map<name_type,parameter_id_type>  BaseParamIDMap;  
class  map_for_id ;
typedef  map_for_id  ParamIDMap;  


// list of things that have been done in a class or instance
typedef std::map<name_type,StrTy>  CheckListMap;  




// I suppose this is kind of dumb as it allows the same
// name to be use many times, should just make this a MAP
// and implement vector like things as needed. 
class vector_with_c_loading : public  BaseParamNameList
{

public:
void load(const char ** pp)
{ const char ** p=pp; while (*p!=0) {push_back(*p); ++p; } }
//may want associative index?
// this will aactually be important very quickly , should justindex...
IdxTy get_index(const name_type & nm) const
{
// wanted a vector here, could index but this should not be a big deal
const IdxTy sz=size();
for (IdxTy i=0; i<sz; ++i)  if ( (*this)[i]==nm) return i; 
return BAD;
}

}; 


class vector_with_name_checking : public BaseParamValuesArray
{
public:
vector_with_name_checking(): m_check(0) {}
// really this is dumb if pnl goes out of scope etc we should just
// copy it bit then risk inconsistent copies etc. 
vector_with_name_checking(const ParamNameList & pnl):m_check(&pnl) {}

return_code  set(const name_type & nm, default_value_type  t)
{
if ( m_check==0) { push_back(t); return 0 ; } 
IdxTy loc=(*m_check).get_index(nm);
if ( loc==(IdxTy)BAD) return BAD;
if ( size()<=loc) resize(loc+1);
(*this)[loc]=t;
return 0;
} // set

const ParamNameList * m_check;
};



class map_with_voids : public BaseParamNameValuePairList
{
public:
template <class To> void set(const name_type & nm, To & v)
{ (*this)[nm]=(void*)&v;} 
};


class  map_for_id : public BaseParamIDMap
{

public:
void load(const char ** nm, const int * param)
{
const char ** p=nm;
int i=0;
while (*p) {(*this)[name_type(*p)]=param[i]; ++i; ++p; }

}

};

// this was a huge debat over templating or not, for now just one base class.
class ParamList  {

typedef ParamList Myt;
// this needs to be a pointer for polymorphism
typedef std::map<ParamListHandle, Myt  * >  HandleMap;
public:
// generally if a ostream exists, print to it else look at flags

enum {NO_COMPUTATIONS=1,BIBTEX_CITATIONS=2,
	HUMAN_CITATIONS=4,COMMENTARY=8, FAIL_ON_PARMS=16,
CHECK_PARAM_AND_FUNC=32
  };

typedef unsigned int IdxTy; // for subclassing 
typedef nv_api::ParamNameValuePairList::const_iterator Ii;
typedef  nv_api::name_type StrTy;
typedef  nv_api::default_value_type ValTy;
typedef  void ( * WorkerFunc )( const IdxTy sz,const IdxTy pcount,  Myt  & pl, ParamValuesArray & ) ; 

ParamList()
	:m_serial(serial()), m_debug(0),m_citations(0),m_err(0),
m_params(0), m_id_map(), m_map(0),m_mode(0),m_worker(0) {}
ParamList(const IdxTy & params, const ParamIDMap & names)
	:m_serial(serial()),m_debug(0),m_citations(0),m_err(0),
	m_params(params), m_id_map(names), m_map(new IdxTy[m_params]),m_mode(0),m_worker(0) {}

~ParamList() { clear();} 

ParamList(const Myt& that) { this->copy(*this,that); }

virtual ParamList & operator=(const Myt& that)
{ Clear(); this->copy(*this,that); return *this; }

void clear() { Clear(); }
// leaving all the implementations in right now, fwiw
// this is a member so it can be changed etc
virtual void copy(Myt& d, const Myt& s) { Copy(d,s); }
// not thread safe and note it is a reference non const
static IdxTy & serial() 
{ static IdxTy count=~0;  ++count; return count; } 
static HandleMap & handle_map() 
{static HandleMap hm; return hm;}
static Myt * get_pl(const ParamListHandle h) { return handle_map()[h]; }
static ParamListHandle & get_next_handle() 
{ static IdxTy next=0; return next; } 
static CheckListMap & check_list() { static CheckListMap clm; return clm; } 
static bool have_done(const StrTy & nm) { return check_list()[nm]!=""; }
static const StrTy & did(const StrTy & nm,const StrTy & v) { return check_list()[nm]+=v; return check_list()[nm]; }

// maybe should be pure virtual but allow this as default
virtual int parse_pair(const StrTy &nm, const ValTy & v) 
{ return Parse_pair(nm,v); }

// probably need not be virtual, 
virtual return_code parse( const ParamNameValuePairList & pnl)
{ return Parse(pnl); } 
// subclass nedds to put something here, virtual =0?
virtual return_code find_a_worker(nv_api::ParamNameList & pnl)
{
MM_MSG(" in default worker find code prbably not intended")
return 0;
}

// this is called by compile, but user can override and then call match1
virtual return_code match( const nv_api::ParamNameList & pnl)
{ return match1(m_id_map,pnl); }

virtual return_code match1( const nv_api::ParamIDMap & idm, const nv_api::ParamNameList & pnl)
{ return Match1(idm,pnl);  }

// not virtual  - try to make the recurring stuff faster
return_code  dispatch(const ParamListHandle plh, const IdxTy sz, ParamValuesArray&   pva)
{
(* m_worker)(sz,0,*this,pva);
return 0; 
}

bool debug() const { return m_debug!=0; }
bool citations() const { return m_citations!=0; }
bool errs() const { return m_err!=0; }

IdxTy m_serial;

OsTy * m_debug, * m_citations, * m_err;

IdxTy m_params;
// map of allowed parameter names to param id code
ParamIDMap m_id_map;

// without template parameter, this is now a copyuing mess..
//IdxTy  m_map[MAPSZ]; // maps input values to existing signature
IdxTy  *m_map; // maps input values to existing signature
IdxTy m_mode; // debug, information/citations only, status etc

WorkerFunc  m_worker;
// could have lots of things to add here in basic or subclasses

// ParamList Impl. wanted to avoid a cpp file....  
private:

void Clear() { if ( m_map!=0) delete [] m_map; }
void Copy(Myt& d, const Myt& s)
{
d.m_params=s.m_params;
d.m_id_map=s.m_id_map;
d.m_mode=s.m_mode;
d.m_map= new IdxTy[d.m_params];
memcpy(d.m_map,s.m_map,d.m_params*sizeof(IdxTy));
}

return_code Parse_pair(const StrTy &nm, const ValTy & v) 
{ return 0; }


return_code Parse_pair_base(const StrTy &nm, const ValTy & v) 
{ 
// should make a map and switch but again this is init code,
// can fix it lateer....
if (nm=="debug-stream") m_debug=(OsTy*)v;
else if (nm=="citation-stream") m_citations=(OsTy*)v;
else if (nm=="cerr") m_err=(OsTy*)v;
return 0; 

}

return_code Parse( const ParamNameValuePairList & pnl)
{
if(debug()) MM_MSG("parsing the libx thing")
//const IdxTy sz=pnl.size();
Ii ii=pnl.begin();
while (ii!=pnl.end())
{
const StrTy &  nm=(*ii).first;
const ValTy &  v=(*ii).second;
if(debug()) MM_MSG(" parsing "<<nm)
// wtf compiler can not find the over xxx?
Parse_pair_base(nm,v);
parse_pair(nm,v);
++ii;
}
return 0; 
}


 return_code Match1( const nv_api::ParamIDMap & idm, const nv_api::ParamNameList & pnl)
{
if(debug()) MM_MSG("matching now")
// these create problems with char**
typedef ParamNameList::value_type St;
typedef ParamNameList::iterator It;
typedef ParamIDMap::const_iterator Im;
clear();
// must be sequentially numbered...
m_map = new IdxTy[idm.size()];
const IdxTy pcount=pnl.size(); 
for ( IdxTy i=0; i<pcount; ++i)
{
const St & nm=pnl[i];
Im loc=idm.find(nm);
//int found=0;
if ( loc!=idm.end()) { // found=1; 
if(debug()) MM_MSG(" matching "<<nm<<" at "<<i<<" to "<<(*loc).second)
m_map[(*loc).second]=i;

 }
else 
{
// this is an error, we also need to see if the required
// params have all been supplied 
MM_MSG(" no parameter for "<<nm<< " at position "<<i)
}


}
if (debug()) MM_MSG("done matcyhing")

return 0;
}


}; //   ParamList;





// this needs to know what kijnd of ParamList tomake 
// I suppose it could be passed, but we wanted to hide all of that from user.
// not thread safe
template<class PL> ParamListHandle compile
( const nv_api::ParamNameList & pnl, const ParamNameValuePairList & nvpl)
{
// not thread safe.
ParamListHandle plh=PL::get_next_handle(); // note this is NOT thread safe.
PL&  pl= * new PL() ; // do not want it to go out of scope, need to delete etc 
MM_MSG(" have a handle of "<<plh)
// this needs to print banners etc.
pl.parse(nvpl);
// distinct from parse for error reporting and params and subclassing. 
pl.find_a_worker(pnl);

MM_MSG(" calling name match now")
// I suppose we could just make one call to pl to do it all, there
// was some point to this. 
pl.match(pnl);

// sore a ptr to this so handle look up later can find it. 
PL::handle_map()[plh]=&pl;
++PL::get_next_handle(); // note this is NOT thread safe.
return plh; 
} // compile
///////////////////////////
// sz should be compiled... 
int  dispatch(const ParamListHandle plh, const IdxTy sz, ParamValuesArray & pva)
{
// should use reference or pointer or something 
MM_MSG("dispatching handle "<<plh);

ParamList& pl=* ParamList::get_pl(plh);
// no idea why we finally did this?
return pl.dispatch(plh,sz,pva);
//(*(WorkerFunc ) pl.worker)(sz,0,&pl,&pva);

}


// basic dumping utility for debug, could also have column headers etc. 
void print_comp_values_hdr(OsTy & os, const ParamNameList all_names,const ParamValuesArray values)
{
const IdxTy allcount=values.size();
for (IdxTy  j=0; j<allcount; ++j)
{
// right now we exclude those with no values
if ( values[j]!=0) { os<<all_names[j]<<" "; } 
}

}

void print_comp_values(const ParamNameList & all_names,const ParamValuesArray & values,const unsigned int sz)
{
OsTy & os=std::cout;
const bool print_hdr=false;
const bool print_names=true;
if ( print_hdr) print_comp_values_hdr(os,all_names,values);
const IdxTy allcount=values.size();
for (IdxTy i=0; i<sz; ++i)
{
os<<i<<" ";
// this is where we need types and should cache array[i] lol. 
for (IdxTy  j=0; j<allcount; ++j)
{
// nulls are ok
if ( values[j]!=0)
{
if ( print_names) os<<all_names[j]<<" ";
os<<(((double*)values[j])[i])<<" ";
} 
} // jh
os<<"\n";// copied from the c code, should fix all of this

} // i 

} //print_comp_values




}; //nv_api ns

#else
#warning compiler thinks the nv_api is included in non-c++ code so ignored.

#endif // c++ version




#endif // guard

