#ifndef MJM_ELECTRODE_GEOMETRY_H__
#define MJM_ELECTRODE_GEOMETRY_H__

#include <vector>
#include "mjm_libmesh/mjm_geometry.h" 
#include "mjm_libmesh/mjm_shapes.h" 

class _perimeter_vector: public std::vector<double>
{
public:
typedef double D;
_perimeter_vector():std::vector<D>(2) {}
_perimeter_vector(const D & x, const D& y):std::vector<D>(0) 
{
push_back(x);
push_back(y);
}

};
typedef _perimeter_vector  perimeter_vector;


class mjm_electrode
{
typedef double D;
    //std::vector< D> rect;
    //std::vector< D> qpoint;
    //push_point(qpoint,x,y);
    //push_rect(rect, top,bottom,left,right);
public:

enum { T=0,B=1,L=2,R=3};
mjm_electrode( const  D &  top,const D & bottom,const D & left,const D& right)
: rectangle(top,bottom,left,right)
, voltage(0)
{ push_rect(rect, top,bottom,left,right);mid(); snapped_rectangle=rectangle;}
mjm_electrode(const D & v,  const  D &  top,const D & bottom,const D & left,const D& right)
: rectangle(top,bottom,left,right)
, voltage(v)
//{ push_rect(rect, top,bottom,left,right); mid();}
{ push_rect(rect, top,bottom,left,right);mid(); snapped_rectangle=rectangle;}

const D & volts() const { return voltage; } 
template <class Ty> int classify(const Ty & qpoint) const
	{
    	int v=ClassifyPointBrace(qpoint,rect);
		return v; 
	}
template <class Ty> int classifyElement(const Ty & el) const
	{
    	int v=ClassifyElementRect(el,rect);
		return v; 
	}
const IdxTy perimeter() const { return  (dx+dy+dy+dx); } 

// return location of point i on perimeter for nominal grid spacings delx and y/
// retrn true if this is still within the range need to cover eleteode. 
// his should just reutrn vecotrs and avoid doing these divisions for each point... 
// a generator of sorts
bool  grid(D & px, D & py, D & nx, D & ny,const IdxTy i,const D & delx, const D & dely) const
{
	// this is not quite right 
	const IdxTy nptsx=IdxTy(dx/delx+.5);
	const IdxTy nptsy=IdxTy(dy/dely+.5);
	const IdxTy n=2*( nptsx+nptsy); 
	
	if (i>=n ) return false;
	if (i<nptsy)
	{
		py=rect[T]-i*dely; px=rect[L]; nx=-1; ny=0;return true;
	}
	if (i<(nptsx+nptsy))
	{
		py=rect[B]; px=rect[L]+(i-nptsy)*delx; nx=0; ny=-1;return true;
		return true;
	}
	if (i<(nptsx+2*nptsy))
	{
		py=rect[B]+(i-nptsx-nptsy)*dely; px=rect[R]; nx=1; ny=0;return true;

		return true;
	}

	py=rect[T]; px=rect[R]-(i-nptsy*2-nptsx)*delx; nx=0; ny=-1;return true;



	return true;
}


// this just attempts to project a relevant point onto the perimeter.
// does not need to be specific to electrode or rect really. 
// needs to be able to compute a normal and total distance from
// origin of surface. 
template <class Ty> int get(const Ty & qpoint, Ty & pt,Ty & n,D & l) const
{
//const D qx=qpoint[0]-midx;
//const D qy=qpoint[1]-midy;

// this should be on the surface but may not ce exact
D& x=pt[0];
x=qpoint[0];
D& y=pt[1];
y=qpoint[1];

// rules out left, then bottom leaving top and right. 
if (( x>midx) && (y>midy)  )
{
D mx=::fabs(x-rect[R]); D my=::fabs(y-rect[T]);
if ( mx<my) { x=rect[R]; l=dy+dx+y-rect[B];n[0]=1; n[1]=0;} 
else { y=rect[T];l=dy+dy+dx+rect[R]-x; n[0]=0; n[1]=1;}
}
if (( x>midx) && (y<midy)  )
{
D mx=::fabs(x-rect[R]); D my=::fabs(y-rect[B]);
if ( mx<my) { x=rect[R]; l=dx+dy+y-rect[B]; n[0]=1; n[1]=0; } 
else { y=rect[B]; l=dy+x-rect[L];n[0]=0; n[1]=-1; }
}
if (( x<midx) && (y<midy)  )
{
D mx=::fabs(x-rect[L]); D my=::fabs(y-rect[B]);
if ( mx<my) { x=rect[L];l=rect[T]-y; n[0]=-1; n[1]=0; } 
else { y=rect[B];l=dy+x-rect[L]; n[0]=0; n[1]=-1; }
}
if (( x<midx) && (y>midy)  )
{
D mx=::fabs(x-rect[L]); D my=::fabs(y-rect[T]);
if ( mx<my) { x=rect[L];l=rect[T]-y; n[0]=-1; n[1]=0; } 
else { y=rect[T]; l=dy+dy+dx+rect[R]-x; n[0]=0; n[1]=1; }
}




return 0; // compiler flagged this as missing, no idea what is does

}
void mid() //centroid?
{
midx=(rect[R]+rect[L])*.5;
midy=(rect[T]+rect[B])*.5;
dx=rect[R]-rect[L];
dy=rect[T]-rect[B];
}
bool contains(const D & x, const D & y) const 
{
	return snapped_rectangle.contains(x,y);
}
void snapped(const SimpleRect & s) const { snapped_rectangle=s;}
// this needs a polygon or shape library
    std::vector< D> rect;
	// right now this is not used just for migration. 
	SimpleRect rectangle;
	mutable SimpleRect snapped_rectangle;
	D midx,midy,dx,dy;
	D voltage;


};

class electrode_geometry
{
public:
typedef double D;
typedef mjm_electrode Elec;
typedef std::vector<Elec> ElecVec;
typedef _perimeter_vector  perimeter_vector;
//electrode_geometry():m_bc_only(true), m_v_plus(.1), m_v_minus(-.1)
electrode_geometry():m_bc_only(true), m_v_plus(.8), m_v_minus(-.8)
,m_y_max(::sqrt(.3)),m_y_min(-m_y_max), m_y_taper(::sqrt(.4))
,m_taper_scale(1.0/(m_y_taper-m_y_max))
{
//    const Real top=.7; const Real bottom=-.7;
    const Real top=.4; const Real bottom=-.405;
    const Real left=.45; const Real right=.675;
electrodes.push_back(Elec(m_v_plus,top,bottom,left,right)); //0 
electrodes.push_back(Elec(m_v_minus,top,bottom,-right,-left)); //1


}
const bool bc_only() const { return m_bc_only;}
const D  vplus() const { return m_v_plus;}
const D  vminus() const { return m_v_minus;}
const D ymax() const { return  m_y_max;}
const D ymin() const { return  m_y_min;}

const D guess_init_value(const D & x, const D & y) const
{
// this should use some estimate of the fields from 
// 2D electrostatics, it may help but hard to say if worht
// effort. 
const D le=electrodes[1].snapped_rectangle[mjm_electrode::R];
const D ri=electrodes[0].snapped_rectangle[mjm_electrode::L];
const D riorig=electrodes[0].rectangle[mjm_electrode::L];
static bool once=false;
if (ri==riorig)
{
if (!once) { MM_ERR(" ri has not been sanpped "<<ri)}
once=true; 

}

D eps=(x-le)/(ri-le);
if (eps>1) eps=1;
if (eps<0) eps=0; 
//return 0; 
return (m_v_minus+eps*(m_v_plus-m_v_minus));
}

// eventually determine if point is inside some electrode
// if true then v was updated with voltage
template <class Ty> const bool epoint(D & v, const Ty & p ) const
{ return false; }


const IdxTy n_electrodes() const { return  electrodes.size(); }
const IdxTy perimeter(const IdxTy e ) const { return  electrodes[e].perimeter(); }


 const bool epoint(D & v, const bool  is_edge, const D & x, const D & y, const D & z=0 ) const
{
if ( !is_edge) return false;
// duplicate the old grap
//const D ys=m_y_max; //::sqrt(ye2);
const D ye=m_y_taper; //::sqrt(ye2taper);
// xe>xs
const D yf=(y>0)?y:(-y);
// this is the distnace from edge to taper to zero.
// assumes symmetry around y=0
//Real w=(ye-yf)/(ye-ys);
Real w=(ye-yf)*m_taper_scale;
//if ( w<0) w=0;
if (w<=0) return false;
if ( w>1.0)  w=1.0;
if ( x>0) { v=w*m_v_plus; return true;}
// fudd this is X not Y...
// lol this left x=====0 gap fudd 
//if ( x<0) 
{ v=w*m_v_minus; return true;}

//return false; 

}





// this includes the edge
 const bool epoint2xxx(D & value, const bool  is_edge, const D & x, const D & y, const D & z=0 ) const
{
// try to cut out grid
//    if ( is_edge) return false;
    const Real top=.7; const Real bottom=-.7;
    const Real left=.5; const Real right=.625;
    std::vector< D> rect;
    std::vector< D> qpoint;
    push_point(qpoint,x,y);
    push_rect(rect, top,bottom,left,right);

    //Real value=m_v_plus;
    int v=ClassifyPointBrace(qpoint,rect);
    if (v>=0) {value=m_v_plus; return true; }
    rect.clear();
    push_rect(rect, top,bottom,-right,-left);
    v=ClassifyPointBrace(qpoint,rect);
    if (v>=0) { value=m_v_minus; return true; }
    return false;
}
// strictly inside
 const bool epoint2exxx(D & value, const bool  is_edge, const D & x, const D & y, const D & z=0 ) const
{
// try to cut out grid
//    if ( is_edge) return false;
    const Real top=.7; const Real bottom=-.7;
    const Real left=.5; const Real right=.625;
    std::vector< D> rect;
    std::vector< D> qpoint;
    push_point(qpoint,x,y);
    push_rect(rect, top,bottom,left,right);

    //Real value=m_v_plus;
    int v=ClassifyPointBrace(qpoint,rect);
    if (v>0) {value=m_v_plus; return true; }
    rect.clear();
    push_rect(rect, top,bottom,-right,-left);
    v=ClassifyPointBrace(qpoint,rect);
    if (v>0) { value=m_v_minus; return true; }
    return false;
}
// this should obsolete the epoint methods
 const bool contains(D & value, const bool  is_edge, const D & x, const D & y, const D & z=0 ) const
{
if(!is_edge) return false; 
const IdxTy sz=electrodes.size();
for (IdxTy i=0; i<sz; ++i)
{
	const Elec & e=electrodes[i];
	if (e.contains(x,y)) 
	{
		if (i==1) value=m_v_minus;
		if (i==0) value=m_v_plus;
		// MM_ERR(" imposing "<<x<<" "<<y<<" "<<value)
		return true;
	}
} // electrode
return false;
}

 const D &  voltage( const IdxTy electrode ) const { return electrodes[electrode].volts(); } 

const bool voltage(D & value, const IdxTy electrode, const D & x, const D & y, const D & z=0 ) const
{

		if (electrode==1) value=m_v_minus;
		else if (electrode==0) value=m_v_plus;
		else return false;

return true; 
}

//if( geo.elec_surface(el,electrode,pt,n,l))
// if the element is relevant to the electrode, return the electrode number, a
// middle point on surface, outward normal, and distance from arbutrary srface origin.
template <class Te, class Tf,class Tp>
const bool elec_surface( const Te & el, const Tf & _qpoint,IdxTy & electrode, Tp & pt, Tp & n, D & l) const
{
//perimeter_vector qpoint;
const Tp  qpoint(_qpoint(0),_qpoint(1));
const IdxTy sz=electrodes.size();
for (IdxTy i=0; i<sz; ++i)
{
	const Elec & e=electrodes[i];
	// thjis should work for any shape. 
	// right now it just calls ClassifyElement on el and rect for the electrode 
	const int cl=e.classifyElement(el);
	if (cl<0) continue;
	// on or inside 
	electrode=i;
	e.get(qpoint,pt,n,l);
	
	return true;
} 
return false; 

}

const D j(const D & x, const D & y) const
{
if ( x>.99) if (yok(y))  return -1000e25;
if ( x<-.99) if (yok(y))  return 1000e25;
return 0;
}

const bool yok(const D & y) const
{return  ((y<m_y_max)&&(y>m_y_min)) ; }

const bool m_bc_only;

const D m_v_plus;
const D m_v_minus;


const D m_y_max;
const D m_y_min;
const D m_y_taper;
const D m_taper_scale;

ElecVec electrodes;

}; // electrode_geometry




#endif


