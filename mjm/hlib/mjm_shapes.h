#ifndef MJM_SHAPES_H__
#define MJM_SHAPES_H__

#include "mjm_globals.h"

#include <vector>
#include <sstream>
#include <string>

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <math.h>




class mjm_shape_generator
{
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
//typedef mjm_block_matrix<D> MyBlock;
//typedef  data_model_error_log Dmel;
//typedef mjm_string_ordering Ordering;
public:
static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); }
static int myatoi(const char * c) { return ::strtol(c,0,0); }
template <class Tv> 
static void ngon( Tv & x, Tv & y, const D & r, const IdxTy & n
, const D & xz=0, const D & yz=0, const D & phiz=0, const IdxTy & m=1)
{
static const D px=M_PI/180.0;
for (IdxTy i=0; i<n; ++i)
{
const D f= 360.0*((i*m)%n)/n;
D xp=r*cos(px*(phiz+f));
D yp=r*sin(px*(phiz+f));
x.push_back(xp);
y.push_back(yp);

} // i  



}


}; // mjm_shape_generator

// put crap in one place for all the ad hoc stuff 
// this needs to work with rational types for D
template <class D> class SimpleRectClass 
{
public:
typedef SimpleRectClass Myt;
//typedef double D;
typedef unsigned int IdxTy;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
enum {T=0,B=1,L=2,R=3};
// this creates a rectangle at zero, but we may want an invalid
//SimpleRect():m_vec(4) {}
SimpleRectClass(const bool valid=true):m_vec(4),m_valid(valid) {}
SimpleRectClass(const D & t, const D & b, const D & l, const D & r):m_vec(4) {
//m_vec.push_back(t); m_vec.push_back(b); m_vec.push_back(l); m_vec.push_back(r);
m_vec[T]=t; m_vec[B]=b; m_vec[L]=l; m_vec[R]=r; 
m_valid=true;
}
SimpleRectClass(const D & x, const D & y):m_vec(4) {
m_vec[T]=y; m_vec[B]=y; m_vec[L]=x; m_vec[R]=x; 
m_valid=true;
}

void set (const D & x, const D & y) {
m_vec[T]=y; m_vec[B]=y; m_vec[L]=x; m_vec[R]=x; 
}
void set(const D & t, const D & b, const D & l, const D & r) {
m_vec[T]=t; m_vec[B]=b; m_vec[L]=l; m_vec[R]=r; 
}

static bool maxx(D & x, const D & val)
{ if (x<val) {x=val; return true; } return false; }
static bool minn(D & x, const D & val)
{ if (x>val) {x=val; return true; } return false; }


bool enclose(const D & x, const D & y)
{
bool expand=false;
if (!m_valid) { set(x,y); m_valid=true; return true; }
expand|=maxx(m_vec[R],x);
expand|=minn(m_vec[L],x);
expand|=maxx(m_vec[T],y);
expand|=minn(m_vec[B],y);
return expand;

}
void shrink(const D & dx, const D & dy)
{
m_vec[R]-=dx;
m_vec[L]+=dx;
m_vec[T]-=dy;
m_vec[B]+=dy;

}
static bool tween( const D & x, const D & xmin, const D & xmax, const bool bottom)
//{if (bottom) return (( x>=xmin)&&(x<xmax));
{if (bottom) return ((( x>xmin))&&(!(x>xmax)));
 //return (( x>=xmin)&&(x<xmax));
 return ((!( x<xmin))&&(x<xmax));
}
static void bbest(  D & d, const D & x,const D & xmin, const D & xmax)
{
 //if (( x-xmin)<=(xmax-x)) d=xmin; else d=xmax;
 if (!(( x-xmin)>(xmax-x))) d=xmin; else d=xmax;
}
static bool in(const PtTy & r, const D & x, const D & y)
{
if (r[T]>y) if (r[B]<y) if (r[L]<x) if (r[R]>x) return true;
return false; 
}
static bool in_or_on(const PtTy & r, const D & x, const D & y)
{
//if (r[T]>=y) if (r[B]<=y) if (r[L]<=x) if (r[R]>=x) return true;
// for the rational types, this was easier than making more operators 
if (!(r[T]<y)) if (!(r[B]>y)) if (!(r[L]>x)) if (!(r[R]<x)) return true;
return false; 
}
static bool near(const PtTy & r, const D & x, const D & y, const D & d)
{
if ((r[T]+d)>=y) if ((r[B]-d)<=y) if ((r[L]-d)<=x) if ((r[R]+d)>=x) return true;
return false; 
}


bool in_or_on(const D & x, const D & y) const { return in_or_on(m_vec,x,y); } 
bool in(const D & x, const D & y) const { return in(m_vec,x,y); } 

bool contains(const D & x, const D & y) const
{ return in_or_on(m_vec,x,y)||near(m_vec,x,y,5e-2); } 
// if you had to pick, return an edge closest to thie rect. 
// there is an edge between these, snapped to best of the 4
//bool snapped(D & xbest, D & ybest, const D & x, const D & y, const D & dx, const D & dy)
template <class Td>
bool bestbounds(Td & xbest, Td & ybest, const Td & x, const Td & y, const Td & dx, const Td & dy)
{
const Td xmax=x+dx;
const Td xmin=x-dx;
const Td ymax=y+dy;
const Td ymin=y-dy;
//const D xmid=x+.5*dx;
//warn const D ymid=y+.5*dy;
PtTy best(4);
if (tween(m_vec[T],y,ymax,false)) bbest(best[T],m_vec[T],y,ymax); else best[T]=m_vec[T];
//if (tween(m_vec[B],ymin,y,true)) bbest(best[B],m_vec[B],ymin,y); else best[B]=m_vec[B];
if (tween(m_vec[B],ymin,y,true)) bbest(best[B],m_vec[B],ymin,y); else best[B]=m_vec[B];
//if (tween(m_vec[B],y,ymax,false)) bbest(best[B],m_vec[B],y,ymax); else best[B]=m_vec[B];
if (tween(m_vec[R],x,xmax,false)) bbest(best[R],m_vec[R],x,xmax); else best[R]=m_vec[R];
//if (tween(m_vec[L],x,xmax,false)) bbest(best[L],m_vec[L],x,xmax); else best[L]=m_vec[L];
if (tween(m_vec[L],xmin,x,true)) bbest(best[L],m_vec[L],xmin,x); else best[L]=m_vec[L];
//if (tween(m_vec[L],x,xmax,true)) bbest(best[L],m_vec[L],x,xmax); else best[L]=m_vec[L];

return in(best,x,y);
}

template <class Td>
bool assfudd(const Td & x, const Td & y, const Td & dx, const Td & dy)
{
PtTy best(4);
best[T]=m_vec[T]+dy;
best[B]=m_vec[B]-dy;
best[R]=m_vec[R]+dx;
best[L]=m_vec[L]-dx;

return in(best,x,y);

}




D & operator[](const IdxTy i)  { return m_vec[i];}
const D & operator[](const IdxTy i) const { return m_vec[i];}
//enum {T=0,B=1,L=2,R=3};
StrTy string() const
{ SsTy ss;
ss<< " T= "<<m_vec[T]<<" ";
ss<< " B= "<<m_vec[B]<<" ";
ss<< " R= "<<m_vec[R]<<" ";
ss<< " L= "<<m_vec[L]<<" ";

return ss.str();
}
PtTy m_vec;
bool m_valid;

}; //SimpleRect
typedef SimpleRectClass<double> SimpleRect;


class SimpleSegment
{



}; // SimpleSegment

class SimplePolygon
{


}; // SimplePolygon



#endif

