#ifndef MJM_DUMPERS_TWO_BASE_H__
#define MJM_DUMPERS_TWO_BASE_H__
 

#include <stdlib.h>
#include <math.h>
#include <cmath>
//ZZ#include <mjm_globals.h>
//#include <mjm_templates.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <vector>
//#include "mjm_io.h"
/*
2013-06 marchywka grid dumpers 
Dumpers for outputing real or complex data from grids into common
formats. Implements xplor for now, 

 http://hincklab.uthscsa.edu/html/soft_packs/msi_docs/insight980/xplor/formats.html
http://en.wikipedia.org/wiki/Bravais_lattice

This is all implemented in header files for ease of integration but
these can be included in a single cpp file and implementations
moved if anyeone cares.

From JDFTX for example, this appears to work,

template < class StrTy> void mjm_dump(const StrTy & dest, const StrTy & comment, const GridInfo & g, const DataRptr & r)
{
typedef mjm_xplor_dumper Du;
Du du;
du.destination(dest);
du.comment(comment);
du.dump(r->dataPref(),g.S[0],g.S[1],g.S[2]);
}


*/



class mjm_local_generic_traits
{

public:

typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef unsigned int IdxTy;
typedef std::runtime_error ErTy;

};


class mjm_dumpers_typedefs
{
public:
typedef mjm_local_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef std::stringstream SsTy;
typedef double D;
//typedef std::complex CoTy;

static OsTy & info() { return std::cout; }


};

class mjm_dumpers 
{
typedef mjm_dumpers Myt;

public:
typedef mjm_dumpers_typedefs Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef Tr::SsTy SsTy;
typedef double D;
typedef std::vector<D> V;
typedef std::vector<V> BV;
typedef std::vector<D> Xform; // initial, final x,y,z
typedef std::vector<StrTy> Text; // initial, final x,y,z
static OsTy & info() { return Tr::info(); }

// where to send the stuff- could be url but now just a file 
void destination( const StrTy & dest) { m_dest=dest;}
// single word name for this dumper
void name( const StrTy & name) { m_name=name;}
const StrTy &  name( ) const { return m_name;}
// some attempt to make flexible interface for coord xforms from
// internal to external, still maing this up. 
// coordinate initialx,finalx, init y etc. 
void range(const Xform & r ) { m_range=r;}
void range(const int x  ) { m_range.push_back(x);}
void xform(const Xform & r ) { m_xform=r;}
void xform(const D & x  ) { m_xform.push_back(x);}

int size_limit() const { return m_size_limit; }
int size_limit(const int x )  { m_size_limit=x; return m_size_limit; }
static const int no_limit()  { return ~0; } 
// this code was never tested, should use instead of above 
// these are supposed to take lattuce vectors in obhects that support either [][] or
// (x,y) accessors. In theory they can then derive transforms to external units
// but this is now missing a scale factor for au to angstroms etc. 
// we also need to specifiy something like a view or iterator ordering
// 
template <class Ty> void vectors(const Ty & v) { Vectors(v); } 
template <class Ty> void commas(const Ty & v) { Commas(v); } 

// for formats that support them. 
void comment(const StrTy & c) { m_comments.push_back(c);} 
mjm_dumpers(): m_size_limit(no_limit()) {}
protected:
// ideally the destination can be handled abstractly depending on what 
// support exists for files, url, pipes, etc. 
// needs a base class... 
OsTy *  GetOs() const
{
StrTy def="xplor.xplor";
if  (m_dest.length()==0) 
{info()<<MM_MARK<<" setting output for "<<m_name<<" to " <<def<<CRLF;}
else def=m_dest;
OsTy * os= new std::ofstream(def.c_str(),std::ios::binary);
// should make sure it openeed etc. 
return os;
}
void CloseOs(OsTy & os) const
{
os.flush();
delete &os; // don't want to leave this here, also support sockets etc. 
}
// This allows imples to manipulate and write to a string stream and periodically
// publish it to a non-manipulated target stream. 
// probably should be Flush 
void flush(OsTy & tos, SsTy & ss) const 
{ ss.flush(); tos<<ss.str(); tos.flush(); ss.str(""); }

// See if the transforms have been specified or provied defaults.
// this really is specific to xplor
template <class Ty > void DumpCheck( const Ty * s, const IdxTy & x, const IdxTy & y, const IdxTy & z)  
{
IdxTy sz=m_comments.size();
if ( sz==0) comment(" comments left blank. Generated from JDFTX"); 
sz=m_range.size();
// good use for switch without breaks LOL
if (sz<2) { range(0); range(x-1); } 
if (sz<4) { range(0); range(y-1); } 
if (sz<6) { range(0); range(z-1); } 
sz=m_xform.size();

if ( sz==0)  { xform(10.0); xform(20.0); xform(10.0); }
if ( sz<6) { xform(90.0); xform(90.0); xform(90.0); } 

}
// derive things from lattic vectors  from object with [][] accessorts
template <class Ty> void Vectors(const Ty & v)
{
// should have a msg
if ( v.size()<3) return; 
for ( IdxTy i=0; i<3; ++i)
{
V x; 
x.push_back(v[i][0]);
x.push_back(v[i][1]);
x.push_back(v[i][2]);
m_bv.push_back(x);
}

}
template <class Ty> void Commas(const Ty & v)
{
// should have a msg
//if ( v.size()<3) return; 
for ( IdxTy i=0; i<3; ++i)
{
V x; 
x.push_back(v(i,0));
x.push_back(v(i,1));
x.push_back(v(i,2));
m_bv.push_back(x);
}

}




template <class Ty> D Len(const Ty & x, const Ty & y) const
{
D d1=x[0]-y[0];
D d2=x[1]-y[1];
D d3=x[2]-y[2];
return ::sqrt(d1*d1+d2*d2+d3*d3);
}
template <class Ty> D Dot(const Ty & x, const Ty & y) const
{ return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];}
template <class Ty> D Phi(const Ty & x, const Ty & y) const
{
static const D pi=4.0*::atan(1.0);
return  180.0/pi*::acos(Dot(x,y)/::sqrt(Dot(x,x)*Dot(y,y)));
}
bool  MakeLengthsAndAngles()
{
if ( m_bv.size()<3) return false;
// in the normal case orthorhombic this is waste...  
xform(::sqrt(Dot(m_bv[0],m_bv[0])));
xform(::sqrt(Dot(m_bv[1],m_bv[1])));
xform(::sqrt(Dot(m_bv[2],m_bv[2])));
xform(Phi(m_bv[0],m_bv[2]));
xform(Phi(m_bv[1],m_bv[2]));
xform(Phi(m_bv[0],m_bv[1]));
return true; 
}


//template <class CB> void DumpLine(const IdxTy i, const IdxTy j, const CB & cb, const complex * dc,const IdxTy states,const IdxTy * seq,OsTy & os)
template <class CB, class Tt > void DumpLine(const IdxTy i, const IdxTy j, const CB & cb, const Tt * dc,const IdxTy states,const IdxTy * seq,OsTy & os)
{
for ( IdxTy state=0; state<states; ++state)
{
switch ( seq[state])
{
case 0: { os<<(*dc).real(); break; } 
case 1: { os<<(*dc).imag(); break; } 
case 2: { os<<i; break; } 
case 3: { os<<j; break; } 
case 4: { os<<" "; break; } 
case 5: { os<<CRLF; break; } 
case 6: { 
double xk=(cb.basis->iGarr[j])[0];   
double yk=(cb.basis->iGarr[j])[1];   
double zk=(cb.basis->iGarr[j])[2];   
os<<(::sqrt(xk*xk+yk*yk+zk*zk));  break; } 

case 7: { os<<(*dc).abs(); break; } 
}; // state
} // state
}




// string specification for destination filename or url etc
StrTy m_dest;
// name for this dumper
StrTy m_name;
mutable Xform m_range;
mutable Xform m_xform;
mutable BV m_bv;
int m_size_limit;
Text m_comments;
} ;

// Dump as ascii test in various line oriented formats.

class mjm_ssv_dumper : public mjm_dumpers
{

public:
template <class Ty> void dump( const Ty * s, const IdxTy & x, const IdxTy & y, const IdxTy & z) 
{
DumpCheck( s, x,  y, z); Dump( s, x,  y, z); 
}

private:

template <class Ty> void Dump( const Ty * s, const IdxTy & x, const IdxTy & y, const IdxTy & z)
const 
{
const IdxTy FLUSH_MASK=(1<<10)-1;
const StrTy sep=" ";
OsTy & tos = * GetOs(); 
SsTy os;
flush(tos,os);
IdxTy bsz=x*y*z;
const IdxTy szlim=size_limit();
//const bool  limit_size=((szlim)!=(no_limit()))&&(szlim<bsz);
const bool  limit_size=(szlim^no_limit())&&(szlim<bsz);
const IdxTy sz=(limit_size)?szlim:bsz;
//IdxTy start=0;
IdxTy j=0;
IdxTy inc=1;
IdxTy n[3],m[3];
n[0]=0; n[1]=0; n[2]=0; 
m[0]=x; m[1]=y;m[2]=z;
for (IdxTy i=0; i<sz; ++i)
{
Fix(os)<<n[0]; os<<sep;
Fix(os)<<n[1]; os<<sep;
Fix(os)<<n[2]; os<<sep;
Float(os)<<s[j];
Eol(os);
if ( 0==(i&FLUSH_MASK)) flush(tos,os);
j+=inc;
++n[0]; 
if (n[0]==m[0]) { n[0]=0; ++n[1]; if ( n[1]==m[1]) { n[1]=0;++n[2];}}
} //plane
flush(tos,os);
CloseOs(tos);
}


enum { FP_DIGITS=18, FP_WIDTH=24,INT_DIGITS=8,INT_WIDTH=8};
template <class Os> Os & Eol(Os & os) const {os<<CRLF; return os;  }
template <class Os> Os & Float(Os & os)  const
{
os<<std::scientific<<std::showpos; 
os.precision(FP_DIGITS); os.width(FP_WIDTH); 
return os;
}
template <class Os> Os & Fix(Os & os) const 
{
os<<std::fixed<<std::right; 
os.precision(5); os.width(8); 
os<<std::fixed<<std::right; 

return os; 
}


};

#ifndef EXCLUDE_JDFTX_STUFF

template <class CB> class mjm_columnbundle_dumper : public mjm_dumpers
{

public:
//typedef ColumnBundle CB;
typedef complex DataTy;
void dump( const CB & cb ) { Dump(cb); } 
private:

void Dump (const CB& cb) 
{
OsTy & tos = * GetOs(); // std::cout;
SsTy os;
const DataTy * d=cb.data();
IdxTy cols=cb.nCols();
IdxTy rows=cb.colLength();
//const IdxTy seq[]={2,4,3,4,0,4,1,5}; const IdxTy states=8;
const IdxTy seq[]={6,4,2,4,3,4,7,4,5}; const IdxTy states=9;
for ( IdxTy i=0; i<cols; ++i)
{
// this assumes the order butshould be ok 
const DataTy * dc=d+cb.index(i,0);
for ( IdxTy j=0; j<rows; ++j)
{
DumpLine(i,j,cb,dc,states,seq,os);
++dc;
}
flush(tos,os);
}
CloseOs(tos);
}







}; // column bundle dumper


template <class CB> class mjm_data_dumper : public mjm_dumpers
{

public:
//typedef ColumnBundle CB;
typedef complex DataTy;
void dump( const CB & cb ) { Dump(cb); } 
private:

void Dump (const CB& cb) 
{
OsTy & tos = * GetOs(); // std::cout;
SsTy os;
const DataTy * d=cb->dataPref();
IdxTy cols=1; // cb.nCols();
IdxTy rows=cb.nElem;
//const IdxTy seq[]={2,4,3,4,0,4,1,5}; const IdxTy states=8;
const IdxTy seq[]={6,4,2,4,3,4,7,4,5}; const IdxTy states=9;
for ( IdxTy i=0; i<cols; ++i)
{
// this assumes the order butshould be ok 
const DataTy * dc=d; // +cb.index(i,0);
for ( IdxTy j=0; j<rows; ++j)
{
for ( IdxTy state=0; state<states; ++state)
{
switch ( seq[state])
{
case 0: { os<<(*dc); break; } 
case 1: { os<<(*dc); break; } 
case 2: { os<<i; break; } 
case 3: { os<<j; break; } 
case 4: { os<<" "; break; } 
case 5: { os<<CRLF; break; } 
case 6: { 
// this needs a coord xform api...
double xk=0; // (cb.basis->iGarr[j])[0];   
double yk=0; // (cb.basis->iGarr[j])[1];   
double zk=0; // (cb.basis->iGarr[j])[2];   
os<<(::sqrt(xk*xk+yk*yk+zk*zk));  break; } 

case 7: { os<<(*dc).abs(); break; } 
}; // state


} // state


++dc;
}
flush(tos,os);

}


CloseOs(tos);
}



}; // data g/r dumper




#endif //  EXCLUDE_JDFTX_STUFF













// Dump in xplor format 
class mjm_xplor_dumper : public mjm_dumpers
{
enum { FP_DIGITS=5, FP_WIDTH=12,INT_DIGITS=8,INT_WIDTH=8};
public:
template <class Ty> void dump( const Ty * s, const IdxTy & x, const IdxTy & y, const IdxTy & z) 
{ DumpCheck( s, x,  y, z); Dump( s, x,  y, z); }

private:

// Use a string stream for formatting and dump to target stream periodicallt.
template <class Ty> void Dump( const Ty * s, const IdxTy & x, const IdxTy & y, const IdxTy & z)
const 
{
OsTy & tos = * GetOs(); // std::cout;
SsTy os;
WriteComments(os);
// need to set FP precision to const size, set int to pad fields
WriteScale(os);
// this is just to avoid manipulating the target stream that may be stdout
flush(tos,os);
// really should use an iterator to hide all of this junk 
// and pick directions and offsets for unit cell origin. 
const IdxTy plane_size=x*y;
IdxTy start=0;
for (IdxTy plane=0; plane<z; ++plane)
{
WritePlane(os,s,plane,start,plane_size);
flush(tos,os);
start+=plane_size;
} //plane
WritePostlog(os);
flush(tos,os);
CloseOs(tos);
}

private:
template <class Os> Os & Eol(Os & os) const {os<<CRLF; return os;  }
template <class Os> Os & Float(Os & os)  const
{
os<<std::scientific<<std::showpos; 
os.precision(FP_DIGITS); os.width(FP_WIDTH); 
return os;
}
template <class Os> Os & Fix(Os & os) const 
{
// move tghe consts to enum 
os<<std::fixed<<std::right; os.precision(5); 
os.width(8); 
os<<std::fixed<<std::right; 

return os; 
}

template < class Os> void WriteComments(Os & os ) const 
{
const IdxTy headers=4;
IdxTy sz=m_comments.size();
os<<""; Eol(os);
//os<<"      "<<sz<<" !NTITLE"; Eol(os);
if ( sz>headers) sz=headers;
for ( IdxTy i=0; i<sz; ++i) { os<<m_comments[i]; Eol(os); }
while ( sz<headers) { os<<""; Eol(os); ++sz;}

}
const D output(const std::complex<double> & c) const { return abs(c); }
const D output(const D & c) const { return c; }

template < class Os> void WriteScale(Os & os ) const 
{
Fix(os);
// size must be 6
// NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX. 
IdxTy sz=m_range.size();
if ( sz!=6) info()<<MM_MARK<<" need 6 values for indixeces"<<CRLF;
// should be able to infer this or deduce during dump
while ( sz<6) {m_range.push_back(0); ++sz;}
for ( IdxTy i=1;  i<6; i+=2)
{
int  cnt=(int)(m_range[i]-m_range[i-1]+1); 
Fix(os)<<cnt;
Fix(os)<<(int)m_range[i-1];
Fix(os)<<(int)m_range[i];
} 

Eol(os); Float(os); 
// should be size 9 
sz=m_xform.size();
if ( sz!=6) info()<<MM_MARK<<" need 6 values for xform "<<CRLF;
for ( IdxTy i=0; i<sz; ++i) { os<<m_xform[i]; }Eol(os); 
// this seems to be critical or it dies in pymol load LOL. 
os<<"ZYX"; Eol(os);
}
template <class Os,class Ty> void WritePlane(Os & os,const Ty * s,
const IdxTy plane,const IdxTy start,const IdxTy plane_size) const 
{
Fix(os)<<plane; Eol(os);

Float(os);
const Ty * sp=s+start;
const IdxTy line_size=6;
IdxTy pos=0;
for ( IdxTy i=0; i<plane_size; ++i)
{
//os<<" ";
Float(os)<<output(sp[i]);
++pos;
if ( pos==line_size) { pos=0; Eol(os);} 

}
if ( pos) { Eol(os);} 


}

template < class Os> void WritePostlog(Os & os, const D x1=0, const D x2=0 ) const 
{
Fix(os)<<(-9999); Eol(os);
Float(os)<<x1<<" "<<x2; Eol(os);

}
};

#ifdef INCLUDE_JDFTX_CODE
template < class StrTy> void mjm_dump(const StrTy & dest, const StrTy & comment, const ColumnBundle & cb)
{

typedef mjm_columnbundle_dumper<ColumnBundle> Du;
Du du;
du.destination(dest);
du.comment(comment);
du.dump(cb);

}
template < class StrTy> void mjm_dump(const StrTy & dest, const StrTy & comment, const GridInfo & g, const DataGptr & r)
{
typedef mjm_xplor_dumper Du;
Du du;
du.destination(dest);
du.comment(comment);
// this needs an inverse or something. 
du.dump(r->dataPref(),g.S[0],g.S[1],g.S[2]);
}
template < class StrTy> void mjm_dump(const StrTy & dest, const StrTy & comment, const GridInfo & g, const DataRptr & r)
{
typedef mjm_xplor_dumper Du;
Du du;
du.destination(dest);
du.comment(comment);
// right now this is unit cell lengths in Angstroms, only work with cubic
const double fac=.5291772;
int i=0, j=0;
du.xform(g.R(i++,j++)*fac);
du.xform(g.R(i++,j++)*fac);
du.xform(g.R(i,j)*fac);
// could put 3 angles here, above should be unit cell lattice vector lengths... 
du.commas(g.R); // this fails to convert au to ang.
// this is not right???? 
du.dump(r->dataPref(),g.S[0],g.S[1],g.S[2]);
}
// this did not call the ssv code lol.
template < class StrTy> void mjm_dump_ssv(const StrTy & dest, 
const StrTy & comment, const GridInfo & g, const DataRptr & r,
 const int size_limit=mjm_dumpers::no_limit())
{
typedef mjm_ssv_dumper Du;
Du du;
du.destination(dest);
du.comment(comment);
du.size_limit(size_limit);
// right now this is unit cell lengths in Angstroms, only work with cubic
const double fac=.5291772;
int i=0, j=0;
du.xform(g.R(i++,j++)*fac);
du.xform(g.R(i++,j++)*fac);
du.xform(g.R(i,j)*fac);
// could put 3 angles here, above should be unit cell lattice vector lengths... 
du.commas(g.R); // this fails to convert au to ang.
du.dump(r->dataPref(),g.S[0],g.S[1],g.S[2]);
}




#endif //  INCLUDE_JDFTX_CODE

#endif

