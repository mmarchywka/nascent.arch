#ifndef MJM_PARANOID_OFSTREAM_H__
#define MJM_PARANOID_OFSTREAM_H__


#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <fcntl.h>

/*

This was supposed to do things like write to temp and atomically move
and allow fast async or sure-thing sync calls but now is
just calls to ostream and some locks lol. 


*/

template <class Tr>
class mjm_paranoid_ofstream 
{
 typedef mjm_paranoid_ofstream Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
 typedef typename Tr::MyBlock  MyBlock;

typedef typename  mjm_thread_util<Tr>::mutex_vector MutexVector;

typedef std::ios_base::openmode OpenMode;

enum { MU_SERIAL, MU_WRITE,MU_SZ};

static MutexVector & mv()
{
static MutexVector m = MutexVector(MU_SZ);  
return m;
}
static void SEnterSerial(const IdxTy i) { mv().enter_serial(i); }
static void SExitSerial(const IdxTy i)  { mv().exit_serial(i); }


static IdxTy  serial( const IdxTy d)
{
SEnterSerial(MU_SERIAL);
static IdxTy s=0;
IdxTy x=s;
s+=d;
SExitSerial(MU_SERIAL);
return x;

}

void Init()
{
m_os=0;
m_mutex_vector = MutexVector(MU_SZ);
m_bad=false;
//m_serial=0;
//std::ofstream ofs(dest,std::ios_base::app);
NewOs(m_fn);
}
void NewOs(const StrTy & fn)
{
delete m_os;
m_os=0;
if (fn.c_str()[0]!=0) m_os= new std::ofstream(m_fn,std::ios_base::app);
}
// do not use 
Myt& operator=(const Myt & x) {return *this; } 
mjm_paranoid_ofstream(const Myt & x) {  } 

public:
mjm_paranoid_ofstream() { Init(); } 
mjm_paranoid_ofstream(const StrTy & d):m_fn(d) { Init(); } 
mjm_paranoid_ofstream(const StrTy & d, const OpenMode & om ):m_fn(d) { Init(); } 
~mjm_paranoid_ofstream() { delete m_os; m_os=0;   } 
const StrTy &  fn( ) { return m_fn; } 
const StrTy &  set_file( const StrTy & fn)
{
m_fn=fn;
NewOs(m_fn);
return m_fn;
}


const StrTy &  set_file( const StrTy & dir, const StrTy & base, const StrTy & ext)
{
Ss ss;
ss<<dir;
ss<<base;
//EnterSerial(MU_SERIAL);
ss<<serial(1); // m_serial;
//++m_serial;
//ExitSerial(MU_SERIAL);
ss<<ext;
m_fn=ss.str();
NewOs(m_fn);
return m_fn;
}

StrTy friendly_file( const StrTy & fn)
{
//Ss ss;
const IdxTy sz=fn.length();
const char * p=fn.c_str();
char c[sz+1];
for(IdxTy i=0; i<sz; ++i)
{
if ( p[i]!=' ') c[i]=p[i];
else c[i]='_';
}
c[sz]=0;
return StrTy(&c[0]);
//return ss.str();
}
/*
StrTy friendly_unique_file( const StrTy & fn)
{
StrTy f1=friendly_file(fn);
return unique 


}
*/

StrTy unique( const StrTy & dir, const StrTy & _fnb, const StrTy & ext)
{
const StrTy fnb=friendly_file(_fnb);
// if dir does not exist or combo is invalid nothing 
// can make it exist lol 
StrTy fn=dir+"/"+fnb+"."+ext;
while (true)
{
//MM_ERR(" testing "<<MMPR(fn))
if ( ! Exists(fn))  break; 

Ss ss;
ss<< dir<<"/"<<fnb<<"_"<<serial(1)<<"._"<<ext;
fn=ss.str();


} // true;
return fn;
}

Myt & operator<<(const StrTy & s) {
// need to use fil lock this just for test
// TODO FIXME doh
SEnterSerial(MU_WRITE);
// (*m_os) <<s;
write_nothrow(s);
SExitSerial(MU_WRITE);
 return *this;  } 


// TODO FIXME there is some reason I did NOT want this...
Myt & operator<<(const IdxTy & s) {
// need to use fil lock this just for test
// TODO FIXME doh
SEnterSerial(MU_WRITE);
// (*m_os) <<s;
write_nothrow(s);
SExitSerial(MU_WRITE);
 return *this;  } 




bool sgood() const { return (*m_os).good(); } 
bool good() const { return sgood()&&!m_bad; } 
// this is supposed to insure that all output has completed not
// just write it to disk but also close or 
void finish()
{
(*m_os).flush(); 

}
template <class Ty> 
//void write_nothrow(const StrTy & v)
void write_nothrow(const Ty & v)
{
try{  
(*m_os)<<v;
if (!sgood()) 
{
m_bad=true;
}
} catch ( ...) { m_bad=true; CatchDot3(); } 
}


void write_nothrow(const char * p, const IdxTy len )
{
try{  
(*m_os).write(p,len);
if (!sgood()) 
{
m_bad=true;
}
} catch ( ...) { m_bad=true; CatchDot3(); } 
}


private:

StrTy m_fn;
OsTy* m_os;
MutexVector m_mutex_vector;
bool m_bad;
struct stat64 m_stat;
//IdxTy m_serial;
//SafeOsTy sos(dest);
//const StrTy v=m_formatter.convert_to_mbox(msg);
//sos<<v;
//sos.flush();
//if (!sos.good())


// https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool Exists(const StrTy & name)
{
// struct stat buffer;   
  return (stat64(name.c_str(), &m_stat) == 0); 

}

bool Mkdir(const StrTy & dir)
{
if (Exists(dir)) return true;
IdxTy e=mkdir(dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
if (Exists(dir)) return true;
return (e==0); 
}

void CatchDot3()
{
   // https://stackoverflow.com/questions/315948/c-catching-all-exceptions
  std::exception_ptr pe = std::current_exception();
  MM_ERR((pe ? pe.__cxa_exception_type()->name() : "null"))
}

#if 0 

void write_to_file( IdxTy & qc, const StrTy & dest,const StrTy & _v, const StrTy & sfs=StrTy() ) // nothrow
{
//MM_ERR(" begin save "<<MMPR(dest))
//MM_ERR(" begin save "<<MMPR(v.length()))
// stupidly we already have this info parsed but need it again... doh
StrTy v;//  =_v; 
const bool save_as_mbox=true;
if ( save_as_mbox) convert_to_mbox(v,_v,sfs);
else v=_v;

enter_serial2(); // TODO file lock, and rememebr not locks during IO WTF...  
try{ // scoping 
std::ofstream ofs(dest,std::ios_base::app);
ofs<<v;
if (!ofs.good()) 
{
//MM_ERR(" no good write e "<<MMPR2(dest,a))
qc=2;
//++fails; 
}
} catch ( ...) {exit_serial2();   /* throw; */ }exit_serial2();  /// scoping
//MM_ERR(" end save "<<MMPR(dest))

} // write_to_file

#endif


}; // mjm_paranoid_ofstream

#endif // MJM_PARANOID_OFSTREAM_H__ 
