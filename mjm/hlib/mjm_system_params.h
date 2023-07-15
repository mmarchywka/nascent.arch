#ifndef MJM__SYSTEM_PARAMS_H__
#define MJM__SYSTEM_PARAMS_H__


// the libmesh variable numbers should be moved out to the mla template.
// etc. Just physics here. 

class system_params
{
public:
// generally cgs and ev
system_params():
 uv(0),mzedv(1),mplusv(2),eminusv(3)

/*
The forward reaction seems to be recombination, 
M+ + e- -> Mo and 
dM+
--- = [Mo]k_r - k_f[M+][e-]
dt 

For a normal reaction, I guess the concentration unit is moles/liter.
And k_r then is in 1/sec and k_f is in (1/s)(1/Molar). 
So k_f needs to be divided by 6e23 atoms/mole and then convert cc to liters
*/

// too low a concentration bombs
//,T(300),kT(.0259),mtotaleq(conc(1e-9)),kf(1e6),kr(1e8)
// the one below was for original testing 
//,T(300),kT(.0259),mtotaleq(conc(1e-29)),kf(1e8),kr(1e10)
// better
//,T(300),kT(.0259),mtotaleq(conc(1e-5)),kf(1e8),kr(1e10)
// wow looks like zero except the neutrals are way out but cool 
//,T(300),kT(.0259),mtotaleq(conc(1e-19)),kf(1e8),kr(1e16)
// BC's still ok but huge values 
//,T(300),kT(.0259),mtotaleq(conc(1e-9)),kf(1e8),kr(1e16)

// seems to start working in NL case 
//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(1e6),kr(1e14)
// ok too 
//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(1e13),kr(1e21)
// this could be a units problems, MOLES vs n ? 6e23 lol?
// these values seem to shield everything with zero field 
// and neutral bulk. 
//,T(300),kT(.0259),mtotaleq(conc(1e-9)),kf(1e23),kr(1e31)
//,T(300),kT(.0259),mtotaleq(conc(1e-5)),kf(1e28),kr(1e36)

// ok so far
//,T(300),kT(.0259),mtotaleq(conc(1e-3)),kf(invconc(1e3)),kr(1e3)
// try to get some current lol 
// no improvmenets 
//,T(300),kT(.0259),mtotaleq(conc(1)),kf(invconc(1e3)),kr(1e3)
// this solved ok but the charges were similar sign at both ends
//,T(300),kT(.0259),mtotaleq(conc(1e2)),kf(invconc(1e6)),kr(1e3)
// no improvment
//,T(300),kT(.0259),mtotaleq(conc(1e-3)),kf(invconc(1e8)),kr(1e3)
// ok now without the added current eqns
//,T(300),kT(.0259),mtotaleq(conc(1e-2)),kf(invconc(1e8)),kr(1e3)
// no imprvoment
//,T(300),kT(.0259),mtotaleq(conc(1e-5)),kf(invconc(1e8)),kr(1e3)
// still none after 5 iteration
//,T(300),kT(.0259),mtotaleq(conc(1e-3)),kf(invconc(1e8)),kr(1e3)
// looks good after 5 iterations, i=3.5414e-17
//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(invconc(1e8)),kr(1e3)
// current is decreasing as [] goes up lol, probably side artifact
//,T(300),kT(.0259),mtotaleq(conc(1e-2)),kf(invconc(1e8)),kr(1e3)
// made delta in numerical ders smaller see if that helps 
//,T(300),kT(.0259),mtotaleq(conc(1e-4)),kf(invconc(1e8)),kr(1e3)
// now this fails 
//,T(300),kT(.0259),mtotaleq(conc(1e-3)),kf(invconc(1e8)),kr(1e3)
// ok 
//,T(300),kT(.0259),mtotaleq(conc(1e6)),kf(invconc(1e8)),kr(1e3)
//,T(300),kT(.0259),mtotaleq(conc(1e1)),kf(invconc(1e8)),kr(1e3)
//,T(300),kT(1.38e-16*T/1.6e-19),mtotaleq(conc(1e-3)),kf(invconc(1e8)),kr(1e3)
//,T(300),kT(1.38e-16*T/1.6e-19),mtotaleq(conc(1e-3)),kf(invconc(1e8)),kr(1e3)
,T(300),kT(1.38e-16*T/1.6e-19),mtotaleq(conc(2e-3)),kf(invconc(1e8)),kr(1e3)
//,T(300),kT(1.38e-16*T/1.6e-19),mtotaleq(conc(2e-3)),kf(invconc(1e8)),kr(1e3)
// the current looks good here as it only did one updatelol 
//,T(300),kT(.0259),mtotaleq(conc(1e38)),kf(invconc(1e8)),kr(1e3)
// no updates,
//,T(300),kT(.0259),mtotaleq(conc(1e42)),kf(invconc(1e8)),kr(1e3)
// voltages 1e5 and carrier densities STILL exactly zero 
//,T(300),kT(.0259),mtotaleq(conc(1e52)),kf(invconc(1e8)),kr(1e3)


// try to expand depletion region 
// this ran all night but I upped the iter limit and killed 
//,T(300),kT(.0259),mtotaleq(conc(1e-5)),kf(invconc(1e3)),kr(1e3)


//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(1e2),kr(1e8)
// init conds unchanged
//,T(300),kT(.0259),mtotaleq(conc(1e-12)),kf(1e25),kr(1e33)


//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(1e23),kr(1e31)


//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(1e3),kr(1e4)
//,T(300),kT(.0259),mtotaleq(conc(1e-1)),kf(0),kr(0)
// never anything useful 
//,T(300),kT(.0259),mtotaleq(conc(1e-30)),kf(1e18),kr(1e20)
// the subtraction may fail, really need to check each to get
// huge range
// this does look like one problem as well as corners 
// using meqd here makes a huge mess even with FooLog
// was working 
// this is close probably 
//,mpluseq(meq(mtotaleq,kf,kr)),mzedeq(mtotaleq-mpluseq),eeq(mpluseq)
// this may be closer
,mpluseq(meq(mtotaleq,kf,kr)),mzedeq(kf/kr*mpluseq*mpluseq),eeq(mpluseq)

//,mpluseq(meqd(mtotaleq,kf,kr)),mzedeq(mtotaleq-mpluseq),eeq(mpluseq)


,q(1.6e-19),q_over_e(1.6e-19/8.854e-14),e_over_q(1.0/q_over_e)
,mumplus(1),ummplus(1.0/mumplus)
,mueminus(100),umeminus(1.0/mueminus)
,mueinstein(1000),umeinstein(1.0/mueinstein) // this is really D for einstein relation...
{
    add_variable("u",vars);
    add_variable("mzed", vars);
    add_variable("mplus", vars);
    add_variable("eminus", vars);
    this->dump(std::cout);
}
// many of these need to be reciprocals to stop divsisions lataer
// where no consts are easy to define.. 
const IdxTy uv,mzedv,mplusv,eminusv;
const Real T;
const Real kT;
const Real mtotaleq;
const Real kf;
const Real kr;
const Real mpluseq;


const Real mzedeq;
const Real eeq;
const Real q,q_over_e,e_over_q;
const Real mumplus,ummplus;
const Real mueminus,umeminus;
// mobility for neutral particle to get D easily... 
const Real mueinstein,umeinstein;
StrVec vars;
// convert moles/liter to n/cc
static Real conc(const Real & molar) { return molar*6.022e23*1e-3; }
static Real invconc(const Real & molar) { return molar/(6.022e23*1e-3); }
// this is likely going to encounver precisions problems rounding to zed
static Real meq(const Real & m, const Real & kf, const Real & kr)
{ return -.5*(kr/kf)*(1.0-::sqrt(1+4*kf/kr*m));}
static Real meqd(const Real & m, const Real & kf, const Real & kr)
{ const Real k=kr/kf; return m+.5*k-::sqrt(k*k*.25+m*k);}

static void add_variable(const StrTy & nm, StrVec & v) { v.push_back(nm); }

template <class Os> Os & dump( Os & os) const
{ 
os<<"vars u="<<uv<<" zed="<<mzedv<<" h+="<<mplusv<<" e-="<<eminusv<<CRLF;
os<<" T="<<T<< " kT="<<kT<<CRLF;
os<<" mtotaleq(n/cc)="<<mtotaleq<<" kf(1/sec)(cc/n)="<<kf<<" kr(1/sec)="<<kr<<CRLF;
os<<" mpluseq(n/cc)="<<mpluseq <<" mzedeq="<<mzedeq<<" eeq="<<eeq<<CRLF;
os<<" q_over_e(F/cm)="<<q_over_e<<" e_over_q="<<e_over_q<<CRLF;
os<<" mumplus(h mobility)="<<mumplus<<" ummplus="<<ummplus<<CRLF;
os<<" mueminus(e mobility)="<<mueminus<<" umeminus="<<umeminus<<CRLF;
os<<" mueinstein( einstein mobilty of neutral)="<<mueinstein<<" umeinstein="<<umeinstein<<CRLF;
os<< "units should be cgs with kT in ev. "<<CRLF;
//And k_r then is in 1/sec and k_f is in (1/s)(1/Molar). 
return os;
 }


    template<class Tx, class Ty, class Tp, class Tz>
    void rhs(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & values,
    const Tz & grads,const Tp &qpoint) const
{
    const Real gr_plus=
        (kf*values[mzedv]-kr*values[mplusv]*values[eminusv]);
    const Real gr_zed=-gr_plus;
                                                                                        


   const Real gr_e=gr_plus;
// all the geo crap has been moved
    const Real jelectrode=0; // g_geo.j(qpoint(0),qpoint(1));
// These need to be moved to the LHS were it makes sense... 
    f_vec(uv) = d2u; //  -cond.e_over_q*(values[cond.mplusv]/bf-values[cond.eminusv]*bf);
    f_vec(mzedv) =  umeinstein*gr_zed;
//    f_vec(cond.mplusv) = cond.ummplus*gr_plus+values[cond.mplusv]*(-f_vec(cond.uv))-grads[cond.mplusv].contract(grads[cond.uv]);
 //   f_vec(cond.eminusv) = cond.umeminus*gr_e-values[cond.eminusv]*f_vec(cond.uv)+grads[cond.eminusv].contract(grads[cond.uv]);
    f_vec(mplusv) = ummplus*gr_plus;
    f_vec(eminusv) = umeminus*gr_e+jelectrode;
MM_MSG(" should not get here ")

}



}; // system_params
////////////////////////////////////////////////////////////////////////////////////



class ds_system_params
{
public:
// generally cgs and ev
ds_system_params():
 uv(0),mzedv(1),sumv(2),diffv(3)
,T(300),kT(.0259),mtotaleq(conc(1e-9)),kf(1e6),kr(1e8)
// the subtraction may fail, really need to check each to get
// huge range
// this does look like one problem as well as corners 
// meq returns inf while meqd is sane for same conditionats
//,sumeq(2.0*meq(mtotaleq,kf,kr)),mzedeq(mtotaleq-.5*sumeq),diffeq(0)
,sumeq(2.0*meqd(mtotaleq,kf,kr)),mzedeq(mtotaleq-.5*sumeq),diffeq(0)
,     q(1.6e-19),q_over_e(1.6e-19/8.854e-14),e_over_q(1.0/q_over_e)
,mumplus(1),ummplus(1.0/mumplus)
,mueminus(100),umeminus(1.0/mueminus)
,mueinstein(1000),umeinstein(1.0/mueinstein) // this is really D for einstein relation...
{
    add_variable("u",vars);
    add_variable("mzed", vars);
    add_variable("sum", vars);
    add_variable("diff", vars);
    this->dump(std::cout);
}
// many of these need to be reciprocals to stop divsisions lataer
// where no consts are easy to define.. 
const IdxTy uv,mzedv,sumv,diffv;
const Real T;
const Real kT;
const Real mtotaleq;
const Real kf;
const Real kr;
const Real sumeq;
const Real mzedeq;
const Real diffeq;
const Real q,q_over_e,e_over_q;
const Real mumplus,ummplus;
const Real mueminus,umeminus;
// mobility for neutral particle to get D easily... 
const Real mueinstein,umeinstein;
StrVec vars;
// convert moles/liter to n/cc
static Real conc(const Real & molar) { return molar*6.022e23*1e-3; }
// this is likely going to encounver precisions problems rounding to zed
static Real meq(const Real & m, const Real & kf, const Real & kr)
{ return -.5*(kr/kf)*(1.0-::sqrt(1+4*kf/kr*m));}
// return the amount dissociated from the total 
static Real meqd(const Real & m, const Real & kf, const Real & kr)
{ const Real k=kr/kf; return m+.5*k-::sqrt(k*k*.25+m*k);}

static void add_variable(const StrTy & nm, StrVec & v) { v.push_back(nm); }

template <class Os> Os & dump( Os & os) const
{ return os<<"vars "<<uv<<" "<<mzedv<<" "<<sumv<<" "<<diffv<<" T="<<T<<
" kT="<<kT<<" mtotaleq="<<mtotaleq<<" kf="<<kf<<" kr="<<kr<<" sumeq="<<sumeq
<<" mzedeq="<<mzedeq<<" diffeq="<<diffeq<<" q_over_e="<<q_over_e<<" e_over_q="<<e_over_q
<<" mumplus="<<mumplus<<" ummplus="<<ummplus<<" mueminus="<<mueminus<<" umeminus="<<umeminus
<<" mueinstein="<<mueinstein<<" umeinstein="<<umeinstein; }

    template<class Tx, class Ty, class Tp, class Tz>
    void rhs(Ty & f_vec,const IdxTy nvar,const Real &d2u,const Tx & values,
    const Tz & grads,const Tp &qpoint) const
{
    const Real gr_plus=1.6e-19*(
        kf*values[mzedv]-kr*.25*
            ( values[sumv]-values[diffv])*(values[sumv]+values[diffv]));
    const Real gr_zed=-gr_plus;
    //const Point & qpoint=fe.q_point[qp];
// all the geo crp should have been moved 
    const Real jelectrode=0; // g_geo.j(qpoint(0),qpoint(1));
    const Real gr_e=gr_plus +jelectrode;
// These need to be moved to the LHS were it makes sense... 
    f_vec(uv) = d2u; //  -cond.e_over_q*(values[cond.mplusv]/bf-values[cond.eminusv]*bf);
    f_vec(mzedv) =  umeinstein*gr_zed; // /cond.q;
if (false)
{
//    f_vec(cond.mplusv) = cond.ummplus*gr_plus+values[cond.mplusv]*(-f_vec(cond.uv))-grads[cond.mplusv].contract(grads[cond.uv]);
 //   f_vec(cond.eminusv) = cond.umeminus*gr_e-values[cond.eminusv]*f_vec(cond.uv)+grads[cond.eminusv].contract(grads[cond.uv]);
}else {

  f_vec(sumv) = (ummplus+umeminus)*gr_plus+jelectrode; // /cond.q; 
    f_vec(diffv) = (ummplus-umeminus)*gr_e-jelectrode; // /cond.q;  
}


} // rhs


}; // system_params

#endif

