#ifndef MJM_LOAD_SAVE_IO_H__
#define  MJM_LOAD_SAVE_IO_H__

class mjm_load_save_io 
{

class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::ostream OsTy;
typedef std::istream IsTy;
typedef std::fistream FisTy;
typedef std::fostream FosTy;

}; // 
public:
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::StrTy StrTy;
typedef Tr::OsTy OsTy;
typedef Tr::IsTy IsTy;
typedef Tr::FisTy FisTy;
typedef Tr::FosTy FosTy;
mjm_load_save_io() : m_os(0),m_is(0) {}
~mjm_load_save_io() {}
void set_is(IsTy * is) { m_is=is; }
void set_os(OsTy * os) { m_os=os; }
void infile(const StrTy & fn) { m_is= new FisTy(fn.c_str()); }
void outfile(const StrTy & fn) { m_os= new FosTy(fn.c_str()); }

void close() // not called from dtor although user could null out streams to preserve
{
//if (m_os) (*m_os).close();
//if (m_is) (*m_is).close();
delete m_os;
delete m_is; 
}
private:

OsTy & Out(const void * blob, const IdxTy len) { (*m_os).write(blob,len); return *m_os; }
OsTy & In( void * blob, const IdxTy len) { (*m_is).read(blob,len); returm *m_is; }
void OutLen( const IdxTy len) { Out(&len,sizeof(IdxTy)); }
IdxTy InLen( )
{
	IdxTy len=0;
	In(&len,sizeof(IdxTy));
	return len;
}

void OutLabel(const StrTy & label )
{
	const char * s=label.c_str();
	const IdxTy len=strlen(s);
	// probably faster to concat locally and make one write call 
	OutLen(len);
	Out(s,len); // NOT len+1 although zed term good idea 
}
StrTy InLabel()
{
	StrTy x;
	const IdxTy len=InLen();
	char c[len+1];
	In(c,len);
	c[len]=0;
	x=StrTy(c);
	return x;
}
template <class Ty> void OutGeneric(const StrTy & label, const Ty & x)
{
OutLabel(label);
const IdxTy len=sizeof(x);
OutLen(len);
Out(&x, len);
}
template <class Ty> void InGeneric(const StrTy & label,  Ty & x)
{
	StrTy lreal= InLabel();
	CheckRead("label ",label,lreal);
	const IdxTy lenreal=InLen();
	const IdxTy len=sizeof(x);
	CheckRead(" len ", len, lenreal); 
	In(&x, len);
}
template <class Ty> void OutGenericPtr(const StrTy & label, const IdxTy n, const Ty * x)
{
OutLabel(label);
const IdxTy len=n*sizeof(Ty);
OutLen(len);
Out(x, len);
}
template <class Ty> void InGenericPtr(const StrTy & label,  Ty * & x)
{
	StrTy lreal= InLabel();
	CheckRead("label ",label,lreal);
	const IdxTy lenreal=InLen();
	// this should check mod etc 
	x = new Ty[lenreal/sizeof(Ty)];
	In(x, lenreal);
}



template <class Ty> 
bool CheckRead( const StrTy & label, const Ty & want, const Ty & got)
{
if (want!=got) 
{
MM_MSG(" rrror reading "<<label<<" wanted "<<want<<" got instead "<<got);
MM_ERR(" rrror reading "<<label<<" wanted "<<want<<" got instead "<<got);
return false; 
}

return true;
}



OsTy * m_os;
IsTy * m_is;

}; //mjm_load_save_io

#endif // guard

