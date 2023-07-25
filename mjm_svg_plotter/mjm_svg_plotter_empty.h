#ifndef MJM_PLOTTER_EMPTY_H__
#define MJM_PLOTTER_EMPTY_H__

#include "mjm_globals.h"
#include "mjm_svg_plotter_basics.h"

mjm_global_credits::credit __credit__mjm_svg_plotter_empty("mjm_svg_plotter_basics"
, "  ");
///////////////////////////////////////////////////////////////////////////
template <class Tr>
class mjm_svg_plotter_empty : public mjm_svg_plotter_base<Tr>
{
 typedef mjm_svg_plotter_empty Myt;
 typedef mjm_svg_plotter_base<Tr> Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
typedef typename Super::Ragged Ragged;
typedef typename Ragged::Line  Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_canned_methods Canned;
typedef mjm_svg_writer Sw;
enum  {DEFS=Super::DEFS,BODY=Super::BODY,POST=Super::POST};
public:
// floating point equality is not assured 
//typedef D Zidx;
typedef IdxTy Zidx;
typedef std::map<D,IdxTy> Zmap;
// a component can have multiple Z-levels, traversal is
// potential mess ... 
typedef typename Super::_Point Point;
typedef typename Super::Points Points;
typedef typename Super::Values Values;
typedef typename Super::Zord Zord;
typedef typename Super::params_type params_type;
public:
// I guess these could be virtual, don't emember ctor order on hierarchy.. 
mjm_svg_plotter_empty(): Super() {
MM_ERR(" ctor adscascas ")

Init();}
virtual ~mjm_svg_plotter_empty() {Dtor();}
IdxTy ctor(const Ragged & r, const StrTy & type, const StrTy & nm, params_type & pt, const IdxTy flags) 
{
// no call super? 
MM_ERR("  ctor fasdfasdf")
return 0;
} // ctor 

// at least override this... 

//virtual IdxTy write_self_svg(OsTy & os, Sw & sw, Myt * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 

virtual IdxTy write_self_svg(OsTy & os, Sw & sw, Super * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{
// this confuses the  read
if (false) {
Ss ss; 
ss<<"first component here  ... "<< MM_MARKF<<" ";
ss<<MMPR3(Super::m_rtti,Super::m_name,Super::m_type)<<MMPR2(zord.dump(),pass);
 os<<sw.comment_text(ss);
os<<CRLF;
}
if (pass==DEFS)
{
//os<<sw.start_text(" test foo",m_xs,m_ys);
//os<<sw.frame_text(m_bg_color,m_xs,m_ys);
//os<<CRLF;
}
else
{
//os<<sw.r_start_text(" test foo",m_xs,m_ys,m_u);
//os<<CRLF;
}
//os<<"<defs>"<<CRLF; 

if (pass==BODY)
{
// this now takes strings... 
//os<<sw.rect_text(0,0,m_xs,m_ys,m_bg_color);
//os<<CRLF;

} // BODY 

return 0; 

} // write_self_svg
virtual IdxTy write_self_end_svg(OsTy & os, Sw & sw, Super * root, Zord & zord, const IdxTy pass,const IdxTy flags) // const 
{

//Ss ss;
//ss<<""<< MM_MARKF<<" ";
// "Super" is wtf
//ss<<MMPR3(Super::m_rtti,Super::m_name,Super::m_type)<<MMPR2(zord.dump(),pass);
// os<<sw.comment_text(ss);
//os<<CRLF;

//if (pass==DEFS){  os<<"</defs>"<<CRLF;  }
//if (pass==POST){  os<<"</svg>"<<CRLF;  }

return 0;

} // write_self_end_svg




#if 0 
template <class Ty,class Tv > void setup(Ty & d, const Tv & def, const StrTy & nm, ReadWriteMap&  rwm,const IdxTy flags=0)
{

//m_title="Set a title" ; rwm.get("title",m_title);
//m_title="Set a title" ; rwm.get("title",m_title);
d=def; rwm.get(nm,d);
return 0; 
} //setup 

#endif

StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
protected:
//bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
//typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

//enum { MAP_MU=0 , MU_SZ};
//mutable MutexVector m_mutex_vector;
//void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
//void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
void Init()
{
Super::m_rtti="frame";
m_xs=2000;
m_ys=2000;
m_u="px";
m_use_old_header=false;
m_bg_color="#ffffff";
MM_ERR("  vcrit adlkj  ")
} // Init
void Dtor()
{
Super::Dtor();
//ClearMap();

} // Dtor 

// MEMBERS
//Zmap m_zmap;
IdxTy m_xs,m_ys;
StrTy m_u,m_bg_color;
bool m_use_old_header;

}; // mjm_svg_plotter_empty


#endif // MJM_PLOTTER_EMPTY_H__ 
