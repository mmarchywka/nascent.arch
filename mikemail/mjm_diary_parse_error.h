#ifndef MJM_DAIRTY_PARSE_ERROR_H__
#define MJM_DAIRTY_PARSE_ERROR_H__

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


template <class Tr> 
class mjm_diary_parse_error
{
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;


public:
mjm_diary_parse_error() :m_w(),m_pos(0),m_msg() {}
mjm_diary_parse_error(const StrTy & w, const IdxTy pos, const char * msg)
:m_w(w),m_pos(pos),m_msg(msg) {}
mjm_diary_parse_error(const StrTy & w, const IdxTy pos, Ss & ss)
:m_w(w),m_pos(pos),m_msg(ss.str()) {}

mjm_diary_parse_error & scope(const StrTy & w) {
Ss ss; ss<<w;  ss<<"::"; ss<<m_w; m_w=ss.str(); return *this;
}
StrTy dump() const { Ss ss; return dump(ss).str(); }
Ss &   dump(Ss & ss) const {
ss<<MMPR3(m_w,m_pos,m_msg); // <<CRLF; 
return ss; }
//StrTy m_w; IdxTy m_pos; const char * m_msg;
StrTy m_w; IdxTy m_pos; StrTy m_msg;
}; // mjm_diary_parse_error



template <class Tr> 
class parse_settings
{
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;


static bool Bit(const IdxTy f, const IdxTy b)   { return  ((f>>b)&1)!=0; }
public:
enum {MARKUP, OUTCLEAN, BLANK, RECIPES, EVAL, PRESERVE };
parse_settings() : m_state(0) {}
parse_settings(const IdxTy s) : m_state(s) {}
void date_range(const StrTy & s, const StrTy & e) { m_start=s; m_end=e; }
void dates(const StrTy & s, const StrTy & e) { m_start=s; m_end=e; }
// output a new Ragged file with the original annotated
// with comment lines for errors 
bool eval() const { return Bit(m_state,EVAL); }
bool markup() const { return Bit(m_state,MARKUP); }
// output a cleaned up version which may include canonical dates and nouns
// with multiple spaced, parens, and "no" nouns removed.
bool outclean() const { return Bit(m_state,OUTCLEAN); }
// put everything in parens for use the next day. 
bool blank() const { return Bit(m_state,BLANK); }
bool preserve() const { return Bit(m_state,PRESERVE); }
bool dont_expand_recipes() const { return Bit(m_state,RECIPES); }
bool date_ok(const StrTy & d ) const { return DateOk(d); }
private:
bool DateOk(const StrTy & d ) const {
bool ok=true;
if ( m_start.c_str()[0]!=0) if (d<m_start) return false;
if ( m_end.c_str()[0]!=0) if (d>m_end) return false;

return ok;
}
IdxTy m_state;
StrTy m_start,m_end;

}; // parse_settings





#endif

