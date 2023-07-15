#ifndef MJM_VIRTUAL_DATA_FRAME_H__
#define MJM_VIRTUAL_DATA_FRAME_H__

 

//#include <tcl.h>
#include <mjm_globals.h>
#include <mjm_templates.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <complex>

// want to move other crap that uses this from tcl classes.
#include <gsl_sf.h>
// this was included somewhere else 
// #define HAVE_FFTW

#include "mjm_fft_container.h"


// this was in by default, no ied
#ifndef LINKING_JDFTX
#define EXCLUDE_JDFTX_STUFF
#warning excluding jdfx from dump code
#else
#warning including  jdfx from dump code
 // needed for the dumpers code
#define INCLUDE_JDFTX_CODE
#include <electronic/Everything.h>
#endif


#include "mjm_dumpers.h"

/*
General container based on JDFTX code but for simpler case of research
on simple operations and memory management. On the one hand, it is great
to keep expensive results but if may be faster to re-create than
swap back in or worse thrash. Pre-loading and writing can be much better than
reliance on VM guessing. Note that threads for disk IO will normaly spend a lot
of time waiting and then do some block moves, hopefully like DMA or something.
This is different from computing threads. i

Want to look at things like Intel site but for now stuff like this would work,

http://stackoverflow.com/questions/11563963/writing-a-binary-file-in-c-very-fast   

*/

namespace mjm_virtual_df_functions
{
template <class Os,class Tdata  >  void foo(Os & os, const Tdata &  x) { os<<x; }
template <class Os, class Tdata >  void foo(std::ostream  & os,const std::complex<double> &  x) { os<<std::abs(x); }
// calling only foo without foo<> does not find the specialization
template <class Os, class Tdata,int n >  void foop(Os & os,const Tdata  &  x)
{ os.precision(n); foo<Os,Tdata>(os,x);  }



}; // 

class mjm_virtual_data_frame 
{
typedef mjm_virtual_data_frame Myt;
friend class mjm_df_iterator<Myt>;
//friend class mjm_df_iterator<>;

typedef mjm_grid_typedefs Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef Tr::D D;
typedef Tr::I I;
typedef Tr::Cmplx Cmplx;

typedef StrTy NameTy;
// really want doble too. 
typedef mjm_complex_gridded_data ContainerTy;
//typedef  mjm_virtual_df_functions funcs;
// something we can swap in and out 
// or push up and down?
// finally concede raw pointer access may be bad for predictive 
// caching vs speculative. 
public:
typedef ContainerTy container_type;
typedef Cmplx data_type; 
mjm_virtual_data_frame ():m_ptr(0),m_actives(0),m_active(0) {}


// add these columns and sequence the input to them. 
void have(const char **  p, const mjm_grid * gr ) {

const char ** px=p;
while ( *px)
{
StrTy nm=StrTy(*px);
//osx()<<MM_MARK<<" making "<<nm<<CRLF;
// these should be using the swapping allocator etc. 
if( m_locations.find(nm)!=m_locations.end())
{
//std::cerr<<MM_MARK<<" adding existing name to df, name = "<<nm<<CRLF;
//osx()<<MM_MARK<<" adding existing name to df, name = "<<nm<<CRLF;
}

ContainerTy * ct= new ContainerTy(*gr);
add(nm,ct);
++px;
} // px

}
void add(const StrTy & nm, ContainerTy * ct)
{
// duplicated.
if( m_locations.find(nm)!=m_locations.end())
{
//std::cerr<<MM_MARK<<" adding existing name to df, name = "<<nm<<CRLF;
//osx()<<MM_MARK<<" adding existing name to df, name = "<<nm<<CRLF;
}
m_locations[nm]=m_containers.size(); // if there was one there, it just caused a memory leak... 
m_names.push_back(nm);
m_containers.push_back(ct);// ref counting would benice here but that onlyh helps
++m_actives;
// if something else is already messed up as a debug tool.
}


void reset() { m_ptr=0; m_active=0; } 
data_type * raw(const StrTy & nm ) { return  m_containers[m_locations[nm]]->raw(); } 
container_type * column(const StrTy & nm ) { 
if ( m_locations.find(nm)==m_locations.end())
{
std::cerr<<MM_MARK<<" request for missing coluimn "<<nm<<" not handled well LOL"<<CRLF;
osx()<<MM_MARK<<" request for missing coluimn "<<nm<<" not handled well LOL"<<CRLF;
return 0; 

}

return  m_containers[m_locations[nm]]; } 
// these can accept inputs and store on bg thread etc. 
Myt & operator<<( const double &  a) { return * this; }
Myt & operator<<( const int &  a) { return * this; }
// this does not work for zero elements... 
// really needs a streaming buffer, going from row-to-column .
// this shold just mkae a local copy and distriute later as we are in write-only
// mode now. 
Myt & operator<<( const Cmplx &  a ) {
//osx()<<MM_MARK<<" "<<a<<" active="<<m_active<<" ptr="<<m_ptr<<CRLF;

 m_containers[m_active]->raw()[m_ptr]=a;Inx(); return * this; }
Myt & operator<<( const StrTy &  a ) { return * this; }
// this need not delete all the containers...

template <class Ti= void > void to_xplor(const StrTy & dest,  container_type * s, Ti * iter=0,const D & scale=.529)
{
data_type * src=s->raw();
to_xplor(dest,src,iter,scale);
}


template <class Ti= void > void to_xplor(const StrTy & dest, const StrTy & name, Ti * iter=0,const D & scale=.529)
{
IdxTy loc=m_locations[name];
data_type * src=m_containers[loc]->raw();
to_xplor(dest,src,iter,scale);
}



template <class Ti= void > void to_xplor(const StrTy & dest, data_type * src, Ti * iter=0,const D & scale=.529)
{typedef mjm_xplor_dumper Du;

Du du;
du.destination(dest.c_str());
du.comment(" test ");
du.xform((*iter).span(0)*scale);
du.xform((*iter).span(1)*scale);
du.xform((*iter).span(2)*scale);
du.dump(src, (*iter).size(0),(*iter).size(1),(*iter).size(2)); 

}
//template <class Os,class Ti= void > void dump(Os & os, const StrTy & label,const Ty & names, Ti * iter=0)
template <class Os,class Ti= void > void dump(Os & os, const StrTy & label, Ti * iter=0)
{ 
const StrTy sep=" ";
std::vector<Cmplx*> ptrs;
osx()<<MM_MARK<<CRLF; osx().flush();
//if ( names.size()==0)
{
for ( IdxTy i=0; i<m_containers.size(); ++i) ptrs.push_back(m_containers[i]->raw()); // arrgh raw again..  

}
//osx()<<MM_MARK<<CRLF; osx().flush();
os<<label;
for ( IdxTy i=0; i<(*iter).rank(); ++i) os<<sep<<(*iter).name(i);
for ( IdxTy i=0; i<m_names.size(); ++i) os<<sep<<m_names[i];
//osx()<<MM_MARK<<CRLF; osx().flush();
os<<CRLF;
//const Cmplx*  sp[]= {&tdf.s[0],&tdf.ztpss[0],&tdf.ztpssr[0],&tdf.ztpssnum[0],grho.raw(),&tdf.rho[0],gradpsi.raw(),&tdf.tau[0],gradrho.raw(),&tdf.sigma[0],0};
//osx()<<lbl<<" "<<"analyticpsi"<<" "<<"ztpssanal"<<" "<<"ztpssanalr"<<" "<<"ztpssnum"<<" "<<"grho"<<" "<<"rho"<<" "<<"gradpsi"<<" "<<"tau"<<" "<<"gradrho"<<" "<<"sigma"<<CRLF;
// this needs to just find the zero....
//if ( iter==NULL) DumpArrays(os, label, sep, sp, 10,int(npts),foop<OsTy,Cmplx,8>);
 DumpArrays(os, label, sep, ptrs, *iter,foop<OsTy,Cmplx,8>);

//osx()<<MM_MARK<<CRLF; osx().flush();
//++i; 

}




~mjm_virtual_data_frame () {}
private:
// generally order will need to be fixed. 
// i\logical groups depend on usage or name, usually a 
// point is a row, a type is a column. 
 std::map<NameTy, ContainerTy> m_items;

std::map<NameTy, IdxTy > m_locations; 
std::vector<ContainerTy* > m_containers;
std::vector<StrTy > m_names;
IdxTy m_ptr,m_actives,m_active;
// force all access through iterators that can call the access junk

void Inx()
{
//++m_ptr; // move to other place? 
++m_active; if (m_active==m_actives){ ++m_ptr;  m_active=0;}
}

 //DumpArrays(os, label, sep, ptrs, *iter,foop<OsTy,Cmplx,8>);
template <class Os,class Ts, class Ti, class Tdata, typename Ttri  > 
static void  DumpArrays(Os & os, const Ts & lbl, const Ts & sep
, const std::vector<Tdata*>  ptrs , Ti & iter,  Ttri *  fptr= mjm_virtual_df_functions::foo<Os,Tdata> )
{
const int narrays=ptrs.size();
const int rank=iter.rank();
int i=0;
//osx()<<MM_MARK<<CRLF; osx().flush();
while ( iter)
{
os<<lbl;
for ( int j=0; j<rank; ++j)
if ( fptr==0) os<<sep<<iter[j];
//else { os<<sep; (*fptr)(os, iter[j]); }
else { os<<sep; os<<iter[j]; }

//osx()<<MM_MARK<<CRLF; osx().flush();
for ( int j=0; j<narrays; ++j)
{
//osx()<<MM_MARK<<" "<<j<<" "<<i<<CRLF; osx().flush();
//os<<MM_MARK<<CRLF; os.flush();
const Tdata & p=(ptrs[j])[i];
if ( fptr==0) os<<sep<<p;
else { os<<sep; (*fptr)(os, p); }
}
os<<CRLF;
++iter;
++i;
}

}








// this should be a manager API
void prefetch() {}
void evict() {}
// see about std allocators. 
void alloc() {}
// maintaining equally size things, keep the block for quick re=iese
void free() {}
void will_need(){}
void done_for_now(){}
// often free just means we are done with this but will need an equaly sized large
// array later, this mean relase the memory to OS. 
void free_and_release(){}
// forcing a ro mode may be helpful to do manually
// or only have const ptr/itor active. 
void read_only() {}
void rmw() {}
void write_policy() {}


}; // mjm_virtual_data_frame


class mjm_fptr_in
{


};
class mjm_fptr_param
{


};
class mjm_fptr_out
{


};

#endif // include guard  


