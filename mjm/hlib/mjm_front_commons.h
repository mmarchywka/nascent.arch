#ifndef MJM_FRONT_COMMONS_H__
#define MJM_FRONT_COMMONS_H__
#include "common_includes.h"

#include "mjm_globals.h"

class globals
{

typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;


typedef std::vector<StrTy> AboutInfo;
typedef StrTy AboutText;
typedef StrTy HelpText;



public:
globals(const StrTy & nm) { 

static  const AboutText &  about() { static StrTy m_about; return m_about; }
static  const HelpText &  help() { static StrTy m_help; return m_help; }
static  const StrTy &  sep() {  return CRLF; }
protected:
static  AboutText &  About() { static StrTy m_about; return m_about; }
static  HelpText &  Help() { static StrTy m_help; return m_help; }

static  AboutText &  About(const StrTy & x) { static StrTy m_about; m_about+=x+sep();return m_about; }
static  HelpText &  Help(const StrTy & x) { static StrTy m_help; m_help+=x+sep();return m_help; }





};
typedef globals global_junk_type;






#endif

