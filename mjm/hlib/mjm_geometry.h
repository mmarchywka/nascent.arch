#ifndef MJM_GEOMETRY_H__
#define MJM_GEOMETRY_H__

#include "mjm_globals.h"

#include <vector>
#include <sstream>
#include <string>


class Perimeter1D
{
typedef Perimeter1D Myt;
typedef double D;
typedef unsigned int IdxTy;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
// mostly information about the pieces of the periphery with the addition
// of a value or density at that point. 
// Integral contribution should be dl*value;
// this may be useful for looking at J across an electrode
// to see where the corrosion is. 
class segment
{
public:
D dl; // incremental distance covered
D l; // single monotonic coordinate along the path
PtTy point; // center or similar 
PtTy normal; // outward on closed paths.
D value;
}; // segment
typedef segment Seg;
typedef std::vector<Seg> SegVec;
public:
Perimeter1D():m_name("def"), digits(32),area(0),integral(0) {}
Perimeter1D(const StrTy & nm):m_name(nm), digits(32),area(0),integral(0) {}

void clear() { m_segs.clear(); } 
void add(const PtTy& pt, const PtTy &n, const D & value, const D l, const D dl, const D di)
{
Seg seg;
seg.point=pt;
seg.normal=n;
seg.value=value;
seg.l=l;
seg.dl=dl;
area+=dl;
integral+=di;
m_segs.push_back(seg);
} 

template <class Ty> void ssv(SsTy & ss, const Ty & v, const StrTy &base,const StrTy & sep) const
{
	IdxTy sv=v.size();
	if (sv>0) ss<<base<<"x"<<sep;
	if (sv>1) ss<<base<<"y"<<sep;
	if (sv>2) ss<<base<<"z"<<sep;
	for (IdxTy i=3; i<sv; ++i) ss<<base<<i<<sep;

}
//		ss<<format_vector(s.point)<<" "<<format_vector(s.normal)<<" "<<s.l<<" "<<s.dl<<" "<<s.value;
StrTy header(const StrTy & label="label",const StrTy & sep=" ") const
{	
	SsTy ss;
	const IdxTy sz=m_segs.size();
	if (sz==0) return ss.str(); // no header?
	const Seg & s=m_segs[0];
	//ss<<"name"<<sep<<"label"<<sep;
	ss<<"name"<<sep<<label<<sep;
	ssv(ss,s.point,"point",sep);
	ssv(ss,s.normal,"norm",sep);
	ss<<"l"<<sep<<"dl"<<sep<<"value";
	return ss.str(); 
}

template < class Os> Os & dump ( Os & os, const StrTy & nm,const bool head=false) const
{
	const IdxTy sz=m_segs.size();
	if (sz==0) return os; // no header?
	if (head) { os<<header(nm," "); os<<std::endl;}
	for (IdxTy i=0; i<sz; ++i)
	{
		os<<m_name<<" "<<nm<<" "<<format_element(m_segs,i);
		os<<std::endl;
	}

	return os;
}
template <class Tv> StrTy format_vector(const Tv & v) const
{
	SsTy ss;
	ss.precision(digits);
	const IdxTy sz=v.size();
	if (sz>0)	ss<<v[0];
	for (IdxTy i=1; i<sz; ++i) { ss<<" "<<v[i]; }	
	return ss.str();

}
template <class Tv> StrTy format_element(const Tv & v, const IdxTy i) const
{
	SsTy ss;
	ss.precision(digits);
	const auto & s=v[i];
	ss<<format_vector(s.point)<<" "<<format_vector(s.normal)<<" "<<s.l<<" "<<s.dl<<" "<<s.value;
	return ss.str();
}


const StrTy m_name;
IdxTy digits;
D area,integral;
SegVec m_segs;



}; // Perimeter1D

// Real etc was defined in libmesh, moved out there for now 
// zero on edge, 1 strictly inside
template <class Real> int Classify1D(const Real & x, const Real &xmin, const Real & xmax)
{
if ( x==xmin) return 0;
if ( x==xmax) return 0;
if ( x>xmax) return -1;
if ( x<xmin) return -1;
return 1;
}
template <class Tx, class Ty> int ClassifyPoint(const Tx & qpoint,const Ty & rect)
{ // return 1 if inside 0 if on 
//const Real xp=qpoint(0);
const auto xp=qpoint(0);
//const Real yp=qpoint(1);
const auto yp=qpoint(1);
int vx=Classify1D(xp,rect[2],rect[3]);
if ( vx<0) return -1;
int vy=Classify1D(yp,rect[1],rect[0]);
if ( vy<0) return -1;

return vx*vy;
}
template <class Tx, class Ty> int ClassifyPointBrace(const Tx & qpoint,const Ty & rect)
{ // return 1 if inside 0 if on 
//const Real xp=qpoint(0);
const auto xp=qpoint[0];
//const Real yp=qpoint(1);
const auto yp=qpoint[1];
int vx=Classify1D(xp,rect[2],rect[3]);
if ( vx<0) return -1;
int vy=Classify1D(yp,rect[1],rect[0]);
if ( vy<0) return -1;
// neither is -1, 0 on edge
return vx*vy;

}
// this needs a rect or polygon class to keep order etc straight. 
template <class Ty, class Tr> void push_rect( Ty & rect,const Tr & top
, const Tr & bottom, const Tr & left, const Tr & right)
{
rect.push_back(top);
rect.push_back(bottom);
rect.push_back(left);
rect.push_back(right);
}


template <class Ty, class Tr> void push_point( Ty & rect,const Tr &x 
, const Tr & y)
{ rect.push_back(x); rect.push_back(y); }

// TBLR
template < class Elem, class Tr> int ClassifyElementRect(const Elem * el , const Tr & rect)
{
typedef unsigned int IdxTy;
//typedef double Real;
const IdxTy sz=(*el).n_nodes();
if ( sz<1) return -1;
//const IdxTy coord=_coord;
//const Node * ni=(*el).get_node(0);
//const auto * ni=(*el).get_node(0);
//Real xmax=(*ni)(coord);
//Real xmin=xmax;
bool inside=false;
bool outside=false;
for ( IdxTy i=0; i<sz; ++i)
{
const auto * ni=(*el).get_node(i);
// just want those that cross of are on...
// zed on edge, 1 inside 
int v=ClassifyPoint(*ni,rect);
if (v==0) return 0;
if (v==1) return 1; // inside=true;
if (v==-1) outside=true;
}
if (inside&& outside) return 1; 

return -1;

}


template < class Elem, class Real> int ClassifyElement(const Elem * el , const Real xzed, const unsigned int  _coord=0)
{
typedef unsigned int IdxTy;
const IdxTy sz=(*el).n_nodes();
if ( sz<1) return -1;
const IdxTy coord=_coord;
//const Node * ni=(*el).get_node(0);
const auto * ni=(*el).get_node(0);
Real xmax=(*ni)(coord);
Real xmin=xmax;

for ( IdxTy i=1; i<sz; ++i)
{
//const Node * ni=(*el).get_node(i);
const auto * ni=(*el).get_node(i);
const Real xnode=(*ni)(coord);
if ( xnode>xmax) xmax=xnode;
if (xnode<xmin) xmin=xnode;

}
//MM_MSG(" fuddkk "<<sz<<" "<<xzed<<" min="<<xmin<<" max="<<xmax)
if ( xmax==xzed) return 1;
//if ( xmin==xzed) return 1;
if ((xmax>xzed)&&(xmin<xzed)) return 2;
return 0;

}







#endif

