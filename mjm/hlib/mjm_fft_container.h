#ifndef MJM_FFT_CONTAINER_H__
#define MJM_FFT_CONTAINER_H__

 

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
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

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


// a mess, baiscally not use except to contain fftw.

class mjm_fft_typedefs
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;
typedef Tr::IdxTy I;
typedef std::complex<double> Cmplx;
typedef std::vector<D> V;
// should be a machine sized thing
typedef unsigned long long MachTy;

};

#ifdef HAVE_FFTW
class fft_container
{
typedef fft_container Myt;

typedef mjm_fft_typedefs Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef Tr::D D;
typedef Tr::I I;
typedef Tr::Cmplx Cmplx;

typedef  fftw_complex FftC;
enum { FWD=FFTW_FORWARD, BACK=FFTW_BACKWARD, ZEDPLAN=FFTW_ESTIMATE};
public:
enum {IDENTITY=0,LPF=1,LPFE=2,HPF=3,XFORM=4};
	fft_container(  const I x, const I y, const I z)
:  m_d(3),m_x(x),m_y(y),m_z(z),m_sz(m_x*m_y*m_z),m_data(new Cmplx[m_sz])
{
Plan3D();
}
	fft_container(  const I x, const I y)
:  m_d(1),m_x(x),m_y(y),m_z(1),m_sz(m_x*m_y*m_z),m_data(new Cmplx[m_sz])
{ Plan2D(); }
	fft_container(  const I x)
:  m_d(1),m_x(x),m_y(1),m_z(1),m_sz(m_x*m_y*m_z),m_data(new Cmplx[m_sz])
{ Plan1D(); }

// npts should match, this is really kind of dump
void xform( Cmplx * out, Cmplx * in, const I npts,const bool fwd)
{
//fftw_plan fp = fftw_plan_dft_1d(m_x,( FftC*)in,(FftC*)out,0,FFTW_DESTROY_INPUT);
// the fftw malloc aligns on sse boundaries... 
//fftw_plan fp = fftw_plan_dft_1d(npts,( FftC*)in,(FftC*)out,fwd?FFTW_FORWARD:FFTW_BACKWARD,FFTW_ESTIMATE);
fftw_plan fp = fftw_plan_dft_1d(npts,( FftC*)in,(FftC*)out,bf(fwd),ZEDPLAN);
fftw_execute(fp);
fftw_destroy_plan(fp);
}
// this appears to work, origuinally the test function was discontinuous and it
// looked odd. A windowed psi(R), samples at +/- R works fine. 
// input points
void upsample( Cmplx * out, Cmplx * in, const I npts,const IdxTy  scale, const double fac)
{
const bool dumb_way=!true;
int opts=npts*scale;
Cmplx temp[opts];
xform(&temp[0],in,npts,true);
int npts2=npts>>1;
const Cmplx f=Cmplx(fac*1.0/npts,0); // this is not right... 
if ( dumb_way)
{ // this code is dumb. 
int is=npts-1;
for ( int i=opts-1; i>=0; --i)
{
for ( int j=scale; j>1; --j) { if ( i>=0) temp[i]=Cmplx(0,0); --i; } 
if ( i>=0) temp[i]=temp[is]*f;  --is; // this shoudl use shift or ++ 
} 

} else
{// move top part up, 
int topart=opts-npts;
// just so overlap works, do this BACKWARDS
for ( int i=npts-1; i>=npts2; --i) temp[i+topart]=f*temp[i];
// stupid factor
for ( int i=0; i<npts2; ++i) temp[i]=f*temp[i];
for ( int i=npts2; i<opts-npts2; ++i) temp[i]=Cmplx(0,0);
//std::cout<<MM_MARK<<"upsample coesfs "<<CRLF;
//for ( int i=0; i<opts; ++i) std::cout<<MM_MARK<<" upsampler  "<<i<<" "<<temp[i]<<CRLF;
} // less dumb way
xform(out,&temp[0],opts,false);

}
// this seems ok, again the original test function was a sampled orbital from 0 to larg r.
// making it go to zero at both ends works better. 
// inut points. 
// FORGOT TO RETEST WITH BW LIMIT IN PLACE !!!!!!!!!!
void downsample( Cmplx * out, Cmplx * in, const I npts,const IdxTy  scale, const double fac)
{
Cmplx temp[npts];
const IdxTy opts=npts/scale;
xform(&temp[0],in,npts,true);
const IdxTy opts2=opts>>1;
const IdxTy topside=npts-opts;
//const Cmplx f=Cmplx(fac*1.0/npts*::sqrt(1.0*scale),0); // this is not right... 
const Cmplx f=Cmplx(fac*1.0/npts,0); // this is not right... 
const bool fuscker=false; 
if ( fuscker)
{
// dumb way 
for ( IdxTy i=opts2; i<opts; ++i) temp[i]=f*temp[i+topside]; 
// need to scale in place stuff too lol 
for ( IdxTy i=0; i<opts2; ++i) temp[i]=f*temp[i]; 

xform(out,&temp[0],opts,false);
} else
{
//std::cout<<MM_MARK<<" DOWNSAMPLE FUCL "<<opts<<" "<<npts<<" "<<opts2<<" f "<<(npts-opts2)<<CRLF;
// whoops forgot to re include this... 
// AFAICT this is ok but need to check... 
for ( IdxTy i=opts2; i<(npts-opts2); ++i) temp[i]=Cmplx(0,0); 
//for( IdxTy i=0; i<npts; ++i) std::cout<<MM_MARK<<" downsample coeg "<<i<<" "<<temp[i]<<CRLF;
xform(&temp[0],&temp[0],npts,false);
for ( IdxTy i=0; i<opts; ++i) out[i]=f*temp[i*scale]; 
}

}
// power is nov the integration constant... 
void differentiate( Cmplx * out, Cmplx * in, const I npts,const double scale,
	 const int mode, const double power=0)
{
enum { DIFF=0, INTEGRATE=1, FRACTIONAL=2 };
Cmplx temp[npts];
// use xform?
xform(&temp[0],in,npts,true);
//fftw_plan fp = fftw_plan_dft_1d(npts,( FftC*)in,(FftC*)&temp[0],FWD,ZEDPLAN);
//fftw_execute(fp);
//fftw_destroy_plan(fp);
Cmplx fi=Cmplx(0,0);
//double df=1;
//double ff=1;
// this is der with respect to index, need dx/dn scale factor too
// scale is then points/distance of 1/distance between points. 
// adding would be faster but round off etc. 
const int  ny=npts>>1;
const D xfn=sqrt(npts*1.0);
switch (mode)
{
case DIFF : { //const double  c=-2.0*M_PI/npts*scale;
osx()<<MM_MARK<<" diff  "<<mode<<CRLF;
const D c=-2.0*M_PI*scale/npts;// one for fft, one for scale? 
//const D c=-2.0*M_PI*scale;// one for fft, one for scale? 
for ( int i=0; i<(int)npts; ++i) {int j=i; // =((i+ny)%((int)npts))-ny;  
if ( i>=ny) j-=(int)npts; 

fi=Cmplx(0,c*j) ;temp[i]*=fi; 
// there is some stpud phase error herr FUDD 
// unit impulse gives npts spike ok but -3.14159 nyquiest component too. 
//if ( i==ny){ temp[0]+=temp[i]; temp[i]=Cmplx(0,0); } 



  } break; }
// this seems to work, the code was hacked since the abs() was maing it look wrong LOL
// complex is also cofusing. 
case INTEGRATE : {  
const D c=2.0*M_PI*scale*npts;
//const D c=2.0*M_PI*scale; // *npts;
temp[0]= power /scale;
for ( int i=1; i<(int)npts; ++i) {
int j=i; if ( i>=ny) j-=(int)npts; 
 fi=Cmplx(0,1.0/(c*j));  temp[i]*=fi;  
 } break; }

// do not know how to handle DC  term LOL
// the prefeactor is still wrong.. 
case FRACTIONAL : { 
// check and fix this crap
const D c=-2.0*M_PI*scale/npts/xfn/::pow(M_PI*M_PI*npts*xfn,power);
temp[0]=0; // ?????
for ( int i=1; i<(int)npts; ++i) {
// could make 2 loops...
int j=i;   if ( i>=ny) j-=(int)npts; 
D cj=c*j; D sgn=(cj<0)?(-1):1;
D radius=::pow(::fabs(cj),power); D angle=sgn*M_PI/2.0*power; 
fi=Cmplx(radius*::cos(angle), radius*::sin(angle)); 
 temp[i]*=fi;   }

break; }





default:
{
osx()<<MM_MARK<<" bad mode "<<mode<<CRLF;

}


}; //mode
// could just call xform
xform(out,&temp[0],npts,false);
//fp = fftw_plan_dft_1d(npts,( FftC*)&temp[0],(FftC*)out,BACK,ZEDPLAN);
//fftw_execute(fp);
//fftw_destroy_plan(fp);

}


// n dim fftjjh
static void xform( Cmplx * out, Cmplx * in, const I npts,const bool fwd, const I rank)
{
int dims[rank];
for ( I i=0; i<rank; ++i) dims[i]=npts;
fftw_plan fp = fftw_plan_dft(rank,dims, ( FftC*)in,(FftC*)out,fwd?FFTW_FORWARD:FFTW_BACKWARD,FFTW_DESTROY_INPUT);
fftw_execute(fp);
fftw_destroy_plan(fp);
}
// need to verify orders erc
template < class GridTy > static void xform(Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid, const bool fwd)
{
int rank=sgrid.rank();
int dims[rank],dumdims[rank];
sgrid.get_dims(dumdims);
for ( int rr=0; rr<rank; ++rr) dims[rr]=dumdims[rank-rr-1];
for ( int rr=0; rr<rank; ++rr) std::cout<<MM_MARK<<" fft rank is "<<dims[rr]<<CRLF;
const IdxTy sz_dest=dgrid.size();
const IdxTy sz_src=sgrid.size();
if ( sz_dest<sz_src) osx()<<MM_MARK<<" grid sizes will not work "<<sz_dest<<" too small for " <<sz_src<<CRLF;
osx()<<MM_MARK<<" grid sizes fft k "<<sz_dest<<" for " <<sz_src<<CRLF;

fftw_plan fp = fftw_plan_dft(rank,dims, ( FftC*)in,(FftC*)out,bf(fwd),ZEDPLAN);
fftw_execute(fp);
fftw_destroy_plan(fp);

}
// xform in only one diretion 

template < class GridTy, class Ti > static void xform_dir(Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid, const bool fwd,
const IdxTy dir)
{
int rank=sgrid.rank();
int dims[rank],dumdims[rank];
sgrid.get_dims(dumdims);
// arrghhhh..... 
for ( int rr=0; rr<rank; ++rr) dims[rr]=dumdims[rank-rr-1];
//for ( int rr=0; rr<rank; ++rr) std::cout<<MM_MARK<<" fft rank is "<<dims[rr]<<CRLF;
// multiplied by all those below 
IdxTy fstride=1;
for ( IdxTy  rr=0; rr<dir; ++rr) fstride*=sgrid.size(rr); 
// it looks like we need to pick one plane and then loop over the others.
//so, for x this means each "next is the y stide. 
// doing x/y plane so we kloop over z,w,u,v etc 
const IdxTy dir2=(dir+1)%rank;
IdxTy fdist=sgrid.stride(dir2);
//IdxTy remaining=sgrid.size()/sgrid.size(dir)/sgrid.size(dir2);
Ti  residual(sgrid);
residual.limit(dir,0);
residual.limit(dir2,0);

int fdim[1];
fdim[0]=dims[dir];
const IdxTy sz_dest=dgrid.size();
const IdxTy sz_src=sgrid.size();
// this is for all in the loop, we only do one thing a  atime.
//const IdxTy fmany=sz_src/dumdims[dir]; // total points div those per seg of interest
const IdxTy fmany=dumdims[dir2]; // total points div those per seg of interest
int  * fbed=0; // no idea how to use this 
if ( sz_dest<sz_src) osx()<<MM_MARK<<" grid sizes will not work "<<sz_dest<<" too small for " <<sz_src<<CRLF;
//osx()<<MM_MARK<<" grid sizes fft k "<<sz_dest<<" for " <<sz_src<<CRLF;
MM_MSG("dir="<<dir<<" dir2="<<dir2<<" dim="<<fdim[0]<<" stride="<<fstride<<" dist="<<fdist)
// effective rank is 1
// this ony seems to work for one plane, a 3 D requires another stride etc. 
while (residual)
{
IdxTy off=residual.ptr();
//std::cout<<MM_MARK<<" "<<residual(0)<<" "<<residual(1)<<" "<<residual(2)<<" "<<off<<CRLF; std::cout.flush();

fftw_plan fp = fftw_plan_many_dft(1,fdim,fmany, ( FftC*)in+off,fbed, fstride,fdist,(FftC*)out+off,fbed,fstride,fdist,bf(fwd),ZEDPLAN);
fftw_execute(fp);
fftw_destroy_plan(fp);

++residual;
}

}




// may as well include up/down sample etc. 
// te ft domain needs own grid, but for now just ignore. 
template <class GridTy,class Ti> static void filter( Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid,
const IdxTy mode, const double * dir , const double param=0)
{
xform(out,in,dgrid,sgrid,true);
const bool need_loop=true; // (mode!=0);
const bool need_inverse=(mode!=XFORM); // (mode!=0);
IdxTy bits=0;
if ( dir==0) bits=IDENTITY;
else bits=LPF;
const IdxTy action=bits;


if ( need_loop)
{
Ti fitor(dgrid);
double fmag=0;
double fscale=dgrid.size()?(1.0/dgrid.size()):1;
const Cmplx cscale=Cmplx(fscale,0);
const Cmplx cscale2=Cmplx(::sqrt(fscale),0);
IdxTy i=0;
IdxTy dims=dgrid.rank(); 
if ( dir) { for ( IdxTy j=0; j<dims; ++j ) fmag+=dir[j]*dir[j]; } 
while (fitor)
{
double f=0;
//bool post=false;
if ( dir ) { for ( IdxTy j=0; j<dims; ++j ) f+=fitor[j]*dir[j]; f = ::fabs(f); } 
else { for ( IdxTy j=0; j<dims; ++j ) f+=fitor[j]*fitor[j];  } 
// rely on compiler to move this out of loop LOL
switch (action)
{
case IDENTITY: { out[i]*=cscale; break;  }
// these need the identity scale too
// want to do this post processing
case LPF: {if ( f>fmag) out[i]=Cmplx(0,0); else out[i]*=cscale; break;  }
case LPFE: {if ( f>=fmag) out[i]=Cmplx(0,0); break;  }
case HPF: {if ( f<fmag) out[i]=Cmplx(0,0); break;  }
case XFORM: { out[i]*=cscale2; break;  }

} // switch 

++i; // partof better itor? 
++fitor;

}
} // loop 
// later the outputs will not be compatible but nice for now. 
if ( need_inverse) xform(out,out,dgrid,dgrid,false);
}


template <class GridTy,class Ti> static void iter_filter( Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid,
//const IdxTy mode, const double * dir )
Ti & iter )
{
xform(out,in,dgrid,sgrid,true);
// the itor should know what to do ...
IdxTy i=0;
while (iter) { iter.sample(out[i]);   ++iter; ++i; }

xform(out,out,dgrid,dgrid,false);
if ( ! iter.needs_post()) return; 
std::cout<<MM_MARK<<" post processing "<<CRLF;
iter.reset(); i=0;
while (iter) { iter.post(out[i]);   ++iter; ++i; }

}

template <class GridTy,class Ti,class Tcrap> static void iter_filter_dir( Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid,
//const IdxTy mode, const double * dir )
Ti & iter, const IdxTy dir )
{
std::cout<<MM_MARK<<" "<<CRLF; std::cout.flush();
xform_dir<GridTy,Tcrap>(out,in,dgrid,sgrid,true,dir);
// the itor should know what to do ...
IdxTy i=0;
MM_MSG("calling SAMPLE") // std::cout<<MM_MARK<<" "<<CRLF; std::cout.flush();
while (iter) { iter.sample(out[i]);   ++iter; ++i; }
MM_MSG("returned from  SAMPLE") // 
//std::cout<<MM_MARK<<" "<<CRLF; std::cout.flush();

xform_dir<GridTy,Tcrap>(out,out,dgrid,dgrid,false,dir);
MM_MSG("dir xform inverted now ") // 
//std::cout<<MM_MARK<<" "<<CRLF; std::cout.flush();
if ( ! iter.needs_post()) return;
//else if (  iter.needs_post()) return;
 
std::cout<<MM_MARK<<" post processing "<<CRLF;
iter.reset(); i=0;
while (iter) { iter.post(out[i]);   ++iter; ++i; }


}


// I guess actually using a template solvesa prbolem. 
template <class GridTy,class Ti> static void differentiate( Cmplx * out, Cmplx * in, const GridTy & dgrid, const GridTy & sgrid,
const double * dir, Ti & fitor )
{
xform(out,in,dgrid,sgrid,true);
//mjm_fft_iterator fitor(dgrid);
// this seems to work now, could just make a integration itor. 
//Ti fitor(dgrid);
IdxTy i=0;
const IdxTy dims=dgrid.rank(); 
while (fitor)
{
double f=0;
for ( IdxTy j=0; j<dims; ++j ) f+=fitor[j]*dir[(j)];
//out[i]*=Cmplx(1.0/(128*128*128),0); // Cmplx(0,f);
out[i]= Cmplx(0,f)*out[i];
// std::cout<<MM_MARK<<" der f="<<f<<" "<<i<<" "<<dims<<" out="<<out[i]<<" k^="<<fitor[0]<<" "<<fitor[1]<<" "<<fitor[2]<<CRLF;

++i; // partof better itor? 
++fitor;
}
xform(out,out,dgrid,dgrid,false);
} 




// not sure m_fp copies etc. 
~fft_container() { fftw_destroy_plan(m_fp); delete [] m_data; } 
private:
const I m_d,m_x, m_y,m_z,m_sz;
const Cmplx * m_data;
fftw_plan m_fp;
//fft_container() {}
//fft_container(const Myt & ) {}
Myt& operator=(const Myt & ) {return *this; }

static int bf(const bool dir)  { return dir?FWD:BACK; } 

FftC * In() { return 0; } 
FftC * Out() { return 0; } 
// want 4D at ome point
void Plan3D(){m_fp=0;}//{ m_fp = fftw_plan_dft_3d(m_x,m_y,m_z,In(),Out(),0,FFTW_DESTROY_INPUT);}
void Plan2D(){m_fp=0;}//{ m_fp = fftw_plan_dft_2d(m_x,m_y,In(),Out(),0,FFTW_DESTROY_INPUT);}
void Plan1D(){m_fp=0;}//{ m_fp = fftw_plan_dft_1d(m_x,In(),Out(),0,FFTW_DESTROY_INPUT);}
/*
Myt & fft()
{
const IdxTy n1=m_dims[0];
const IdxTy n2=m_dims[1];
const IdxTy n3=m_dims[2];
const IdxTy sz=n1*n2*n3;
fftw_complex * xd=(fftw_complex*)m_data;
fftw_complex * out = new fftw_complex[sz];
fftw_complex * in = xd;
// check the last 2 params... 
fftw_plan fp = fftw_plan_dft_3d(n1,n2,n3,in,out,0,FFTW_DESTROY_INPUT);
fftw_execute(fp);
fftw_destroy_plan(fp);
m_data=(char*)out;
delete[] xd;
return *this;
}
*/
}; // fft_continer

#else // HAVE_FFTW


#endif // HAVE_FFTW





#endif // include guard  


