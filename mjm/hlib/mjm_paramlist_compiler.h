#ifndef PARAMLIST_COMPILERH__MJM
#define PARAMLIST_COMPILERH__MJM

 
/*

 2014-04-08 skeleeton of usable interface based on r_interface.h see below
This is intended to be c-like avoiding c++ niceties but may require
c++ for libxc impl LOL.
i R-code included for illustration
 

R interface for libxc testing. For illustation only, 
Skeleton code suggestive of a way to use existing libxc test code to interface
with some R methods. Really could be integrated better, I'm just using this
with my own tcl script and dont care about speed at this interface but
it could be a big issue in real life. 

 Marchywka 3013-11-2x


*/


#if 0
http://stackoverflow.com/questions/6658168/passing-a-data-frame-from-to-r-and-c-using-call
http://r.789695.n4.nabble.com/Call-and-data-frames-td911520.html
http://r.789695.n4.nabble.com/Call-and-data-frames-td911520.html
https://stat.ethz.ch/pipermail/r-help/2003-May/034044.html

#endif


#include "xc.h"


#ifdef __cplusplus 
#warning thinks this is cplusplus 
extern"C"{
#else
#warning thinks this is c only NOT cplusplus 
#endif
typedef unsigned int IdxTy;
typedef IdxTy ParamListHandle;
typedef void ** ParamValuesArray;
typedef const char ** ParamNameList;
// why no easy way to take this from global?????j
     // xc_mgga(&func, 1, xc.rho, xc.sigma, xc.lapl, xc.tau, pzk, pvrho, pvsigma, pvlapl, pvtau, pv2rho2, pv2sigma2, pv2lapl2, pv2tau2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigmalapl, pv2sigmatau, pv2lapltau);
enum { LXC_RHO=0,LXC_PZK,LXC_PVRHO,LXC_PV2RHO2,LXC_PV3RHO3,
LXC_SIGMA, LXC_PVSIGMA, LXC_PV2RHOSIGMA, LXC_PV2SIGMA2, LXC_LAPL,LXC_TAU,
LXC_PVLAPL,LXC_PVTAU,LXC_PV2LAPL2,LXC_PV2TAU2,LXC_PV2RHOLAPL,
LXC_PV2RHOTAU, LXC_PV2SIGMALAPL, LXC_PV2SIGMATAU, LXC_PV2LAPLTAU,
 LXC_FUNC,LXC_NSPIN,LXC_LAST_VALUE } ;

enum { LXC_BAD=~0};

const char * libxc_names[] ={"rho","pzk","pvrho","pv2rho2","pv3rho3","sigma",
"pvsigma","pv2rhosigma","pv2sigma2","lapl","tau","pvlapl","pvtau","pv2lapl2",
"px2tau2","pv2rholapl","pv2rhotau","pv2sigmalapl","pv2sigmatau","pv2lapltau",
// "pv2rhotau","pc2sigmalapl","pv2sigmatau","pv2lapltau",
"lap","e","drho","dlap","dtau",0};
const int  libxc_params[]={ LXC_RHO,LXC_PZK,LXC_PVRHO,LXC_PV2RHO2,LXC_PV3RHO3,
LXC_SIGMA, LXC_PVSIGMA, LXC_PV2RHOSIGMA, LXC_PV2SIGMA2, LXC_LAPL,LXC_TAU,
LXC_PVLAPL,LXC_PVTAU,LXC_PV2LAPL2,LXC_PV2TAU2,LXC_PV2RHOLAPL,
LXC_PV2RHOTAU, LXC_PV2SIGMALAPL, LXC_PV2SIGMATAU, LXC_PV2LAPLTAU,
 LXC_LAPL,LXC_PZK,LXC_PVRHO,LXC_PVLAPL,LXC_PVTAU } ;

// make all the worker functions just take arrays of pointers along with required
// parameter count  pcount and array size sz 



typedef void (* stupid) (); 

typedef struct  {
//class ParamList  {
// stuct not support enum?????????
// does c support a dtor as we own this 
IdxTy  m_map[LXC_LAST_VALUE]; // maps input values to existing signature
IdxTy m_function_idx; // right now envisionzed as a huge switch statement
// should also store names and other const info maybe functional ans spin
// etc
  xc_func_type func;
IdxTy nspins,functional;

// fudding void start not work with memory fudd 
//WorkerFunc * worker;
//void * worker;
stupid  worker;

}  ParamList;
typedef  void ( * WorkerFunc )( const IdxTy sz,const IdxTy pcount,  ParamList *pl, ParamValuesArray * ) ; 
typedef struct
{
const char ** names;
const void** values;

} ParamNameValuePairList;
//}  ;

// c++ STL would be a bit better here lol... 
IdxTy mjm_global_param_count=0;
IdxTy mjm_global_param_count_limit=100;
// this does not like the variable ??? global scope issue 
ParamList mjm_global_paramlist[100];

int parse(ParamList *pl, const ParamNameValuePairList * pnl)
{
int i=0;
const char * p=(*pnl).names[i];
while (p)
{
if ( strcmp(p,"functional")==0) (*pl).functional=*(int*)(*pnl).values[i];
else if ( strcmp(p,"nspin")==0) (*pl).nspins=*(int*)(*pnl).values[i];
else printf ( " no match for %s\n", p);
 p=(*pnl).names[++i];
} //p 

return 0;
}
double * cx(const ParamList* pl, ParamValuesArray  * pvap,const int point, const int idx);


void lda_worker( const IdxTy sz,const IdxTy pcount,  ParamList * pl, ParamValuesArray *pva )  
{
printf("lda worker\n"); fflush(stdout);
//if ( 1 == 1 ) return; 
int i;
for (i=0; i<sz; ++i)
{
      //xc_lda(&func, 1, xc.rho, pzk, pvrho, pv2rho2, pv3rho3);
      xc_lda(&(*pl).func, 1, cx(pl,pva,i,LXC_RHO),
      					 cx(pl,pva,i,LXC_PZK),
      					 cx(pl,pva,i,LXC_PVRHO),
      					 cx(pl,pva,i,LXC_PV2RHO2),
      					 cx(pl,pva,i,LXC_PV3RHO3)
);
}
} // lda_worker

void gga_worker( const IdxTy sz,const IdxTy pcount,  ParamList * pl, ParamValuesArray *pva )  
{
printf("gga worker\n"); fflush(stdout);
int i;
for (i=0; i<sz; ++i)
{
    //  xc_gga(&func, 1, xc.rho, xc.sigma, pzk, pvrho, pvsigma, pv2rho2, pv2rhosigma, pv2sigma2);

 xc_gga(&(*pl).func, 1, cx(pl,pva,i,LXC_RHO),
      					 cx(pl,pva,i,LXC_SIGMA),
      					 cx(pl,pva,i,LXC_PZK),
      					 cx(pl,pva,i,LXC_PVRHO),
      					 cx(pl,pva,i,LXC_PVSIGMA),
      					 cx(pl,pva,i,LXC_PV2RHO2),
      					 cx(pl,pva,i,LXC_PV2RHOSIGMA),
      					 cx(pl,pva,i,LXC_PV2SIGMA2)
);
}
}

void mgga_worker( const IdxTy sz,const IdxTy pcount,  ParamList * pl, ParamValuesArray *pva )  
{
printf("mgga worker\n"); fflush(stdout);
printf("mgga worker\n"); fflush(stdout);
printf("mgga worker\n"); fflush(stdout);
printf("mgga worker\n"); fflush(stdout);
int i;
for (i=0; i<sz; ++i)
{
     // xc_mgga(&func, 1, xc.rho, xc.sigma, xc.lapl, xc.tau, pzk, pvrho, pvsigma, pvlapl, pvtau, pv2rho2, pv2sigma2, pv2lapl2, pv2tau2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigmalapl, pv2sigmatau, pv2lapltau);
     
 xc_mgga(&(*pl).func, 1, cx(pl,pva,i,LXC_RHO),
      					 cx(pl,pva,i,LXC_SIGMA),
      					 cx(pl,pva,i,LXC_LAPL),
      					 cx(pl,pva,i,LXC_TAU),
      					 cx(pl,pva,i,LXC_PZK),
      					 cx(pl,pva,i,LXC_PVRHO),
      					 cx(pl,pva,i,LXC_PVSIGMA),
      					 cx(pl,pva,i,LXC_PVLAPL),
      					 cx(pl,pva,i,LXC_PVTAU),
      					 cx(pl,pva,i,LXC_PV2RHO2),
      					 cx(pl,pva,i,LXC_PV2SIGMA2),
      					 cx(pl,pva,i,LXC_PV2LAPL2),
      					 cx(pl,pva,i,LXC_PV2TAU2),
      					 cx(pl,pva,i,LXC_PV2RHOSIGMA),
      					 cx(pl,pva,i,LXC_PV2RHOLAPL),
      					 cx(pl,pva,i,LXC_PV2RHOTAU),
      					 cx(pl,pva,i,LXC_PV2SIGMALAPL),
      					 cx(pl,pva,i,LXC_PV2SIGMATAU),
      					 cx(pl,pva,i,LXC_PV2LAPLTAU)
);
}
}


// fuddkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
int find_a_worker(ParamList * pl,ParamNameList pnl)
{
  switch((*pl).func.info->family)
{
    case XC_FAMILY_LDA: {  (*pl).worker=lda_worker; break;}
    case XC_FAMILY_GGA: {  (*pl).worker=gga_worker; break;}
    case XC_FAMILY_MGGA: {  (*pl).worker=mgga_worker; break;}
default :{  printf(" bad family for functional\n"); return LXC_BAD; }
}

return 0;
}
void match_names(ParamList * pl, const ParamNameList pnl,const IdxTy pcount)
{
int i;
// the pnl should also be null terminated. 
printf(" trying to match %d",pcount);
fflush(stdout);
for ( i=0; i<pcount; ++i)
{
const char * nm=pnl[i];
if (nm==0) return; // should have accuraet count however.. 
printf(" looking up %d %s\n",i,nm);
fflush(stdout);
int j=0;
int found=0; 
const char * lnm=libxc_names[j];
while (lnm!=0)
{
if ( strcmp(lnm,nm)==0) {
int fu=libxc_params[j];
printf(" found match at %d\n",fu);
 (*pl).m_map[libxc_params[j]]=i;++found;
printf(" setting m_map %d to %d\n",libxc_params[j],i);

}
++j;
lnm=libxc_names[j];
} // j 
if ( found!=1) printf(" bad matches for %s of %d at col %d\n",nm,found,i);
} // i 

}
 
// not thread safe
// this really should also require a size and maybe a functional name
// perhaps also have a const parameter list for bind or somethign
ParamListHandle compile(const IdxTy pcount, const ParamNameList pnl,
const ParamNameValuePairList nvpl
)
{
ParamListHandle plh=mjm_global_param_count; // note this is NOT thread safe.
ParamList pl;
printf("parsing\n");
parse(&pl,&nvpl);
printf("finding a worker\n");
fflush(stdout);
printf("try to find funcational %d and %d \n ",pl.functional,pl.nspins);
fflush(stdout);
// this needs to be done first to determine the parameter line ups. 
 if(xc_func_init(&pl.func, pl.functional, pl.nspins) != 0){ return LXC_BAD; } 
int x=find_a_worker(&pl,pnl);
fflush(stdout);
printf(" return code from finding worker is %d\n",x);
fflush(stdout);
// leaking funcs?


printf("match names %d\n",pcount);
match_names(&pl,pnl,pcount);
printf("done match names %d\n",pcount);
fflush(stdout);

mjm_global_paramlist[plh]=pl;
++mjm_global_param_count;  // not thread safe. 
return plh;


}
// need an object lol 
double * cx(const ParamList* pl, ParamValuesArray  * pvap,const int point, const int idx)
{
printf(" looking to cx for point %d and idx %d\n",point,idx);fflush(stdout);
printf(" looking to cx for which is %d and 0,0 %g \n",(*pl).m_map[idx],
((double**)*pvap)[(*pl).m_map[0]][0]

);fflush(stdout);

return  & ((double**)*pvap)[(*pl).m_map[idx]][point];
}





int  dispatch(const ParamListHandle plh, const IdxTy sz, ParamValuesArray  pva)
{
// should use reference or pointer or something 
printf("dispatching handle %d\n",plh);
fflush(stdout);

ParamList  pl=mjm_global_paramlist[plh];
//int i;
/// these could be const for arry making this much faster
//int * functional_list=(int*)pva[pl.m_map[LXC_FUNC]];
//int * nspins=(int*)pva[pl.m_map[LXC_NSPIN]];
// printf("dispatching 2 with info %d %d %u %u  %u\n", pl.functional,pl.nspins,pl.worker, &lda_worker,&mgga_worker);
// fflush(stdout);

// QWHTAT THE FKk

//for (i=0; i<sz; ++i)
//{
// FKKKKKKKK
//void lda_worker( const IdxTy sz,const IdxTy pcount,  ParamList * pl, ParamValuesArray *pva )  
//(*(WorkerFunc *) pl.worker)(sz,0,&pl,0);
printf(" fudd this shot MMMMMMMMMMMMMMMMMMMMMMMM \n"); fflush(stdout);
//lda_worker(sz,0,&pl,0);
//(*(WorkerFunc ) &lda_worker)(sz,0,&pl,0);
(*(WorkerFunc ) pl.worker)(sz,0,&pl,&pva);
printf(" fudd this shot \n"); fflush(stdout);
// newer
// these could be const for this sig 
//  xc_values_type xc;
/*
  switch(pl.func.info->family)
    {
    case XC_FAMILY_LDA:
      //xc_lda(&func, 1, xc.rho, pzk, pvrho, pv2rho2, pv3rho3);
      xc_lda(&pl.func, 1, cx(&pl,&pva,i,LXC_RHO),
      					 cx(&pl,&pva,i,LXC_PZK),
      					 cx(&pl,&pva,i,LXC_PVRHO),
      					 cx(&pl,&pva,i,LXC_PV2RHO2),
      					 cx(&pl,&pva,i,LXC_PV3RHO3)
);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga(&func, 1, xc.rho, xc.sigma,
	     pzk, pvrho, pvsigma, pv2rho2, pv2rhosigma, pv2sigma2);
      break;
    case XC_FAMILY_MGGA:
      xc_mgga(&func, 1, xc.rho, xc.sigma, xc.lapl, xc.tau,
	      pzk, pvrho, pvsigma, pvlapl, pvtau,
	      pv2rho2, pv2sigma2, pv2lapl2, pv2tau2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigmalapl, pv2sigmatau, pv2lapltau);
      break;
 */
//   }


 // xc_func_end(&func);

//} // i 
return 0;
}



#ifdef __cplusplus 
}
#endif

// this is really legacy code
#ifdef INCLUDE_R_INTERFACE_EXAMPLE

#include <Rdefines.h>


// mapping names between input and output arrays. Only a few things suppoorted right now. 
typedef struct  {

const char ** names;
const int * in_func_map;
const char ** results;
const int *out_func_map;
int NFUNC;
int NOUT;

}   ParamLocator;

// just in case it will link LOL. 
#ifdef __cplusplus 
extern"C"{
#endif

// Main code, call "foo" for each row in the dataframe aafter finding 
// names to match as specified in the arrays. See the sample
// calling program for details.
// The same data frame is used for input and output but 
// columns are specified by name and may be used for both.
// 
// df_in : the R data frame containing input and output columns
// foo : the function to be called on each row
// names : names of input columns that may be df_in
//    must have a null string at the end
// in_func_map: a list of indicies corredsponding to names
// indicaing which subscript is used to send name to foo.
// results, out_func_map: as with input indicating output field
//    names and locations. 
// sorry about mixed conventions, should fix names.. 
SEXP find_xc(SEXP df_in, int ( *foo)(int,double *, int,int, double *  )
,const char * names[], const int in_func_map[],const char * results[],
const int out_func_map[]
,const int NFUNC, const int NOUT
);

// alternative with params in struct
SEXP find_xc_param(SEXP df_in, int (*foo)(int,double*,int,int,double*) , const ParamLocator * pl)
{
return 
 find_xc( df_in,foo ,pl->names,pl->in_func_map,pl->results, pl->out_func_map ,pl->NFUNC, pl->NOUT);
}

///////////////////////////////////////////////
// 
// I don't remember what c has for associative containers, no c++ here...
// This should find a matching name to cname from column i in data frame
// and assign column number to results [ program parameter number ].
void name_map(const char * cname,const char **  results,int * osrc,const int* map,int i)
{
int idx=0;
while ( results[idx])
{
if (strcmp(cname,results[idx])==0)
{
// the value for src[x] is found in column i of the data frame
osrc[map[idx]]=i;
//printf(" %s mapping osrc[ %d ] = %d\n",cname,idx,osrc[idx]);
return;
}

++idx;
}

}


// the actual code

SEXP find_xc(SEXP df_in, int ( *foo)(int,double *, int,int, double *  )
,const char * names[], const int in_func_map[],const char * results[],
const int out_func_map[]
,const int NFUNC, const int NOUT
)
{
int i;
    SEXP df_out;

// this seems to give number of columsn ok 
const int size=Rf_length(df_in);
const int ncols=size;
// the columns are assumed to be double no ints for now. 
// really nspin and function should be moved out anyway. 
double * colptrs[ncols];
//printf("df size %d \n",size);
SEXP colnames = getAttrib(df_in, R_NamesSymbol) ;
//    PROTECT( nms = GET_COLNAMES( GET_DIMNAMES( df_in ) ) );
int row=0;
int rows=0;
int isrc[NFUNC],osrc[NOUT];
// put result y into dfin colum osrc[y]
// take column isrc[x] and use that index from df_in for param x to func
for ( row=0; row<NFUNC; ++row) isrc[row]=-1;
for ( row=0; row<NOUT; ++row) osrc[row]=-1;
for ( i=0; i<ncols; ++i)
{
	const char * cname  = (CHAR(STRING_ELT(colnames,i))) ;
	// this is really dumb but we have no assoc container AFAIC
	name_map(cname,names,isrc,in_func_map,i);
	name_map(cname,results,osrc,out_func_map,i);

	SEXP coldata = VECTOR_ELT(df_in,i); // (data for i-th column) 
	// just assume for now
	//if(isReal(coldata)) 
	{ 
		rows=Rf_length(coldata);
		colptrs[i]=REAL(coldata);
	}
}
PROTECT( df_out = NEW_NUMERIC( rows ) );
//printf("df dims are '%d' and '%d' rows '%d' \n\n\n\n",n,p,rows);
double inputs[NFUNC],outputs[NOUT];

for ( row=0; row<rows; ++row)
{
// this is really inefficient if the columns are long as memory thrashing is a mess. 
	int ii;
	for ( ii=0; ii<NFUNC; ++ii)
	{
		// the lut index is the column 
		const int idx=isrc[ii];
		if ( idx<0) inputs[ii]=0.0; else inputs[ii]=colptrs[idx][row];
//printf(" setting %d to %f from col %d\n",ii,inputs[ii],idx);	

	} // ii 
	// no real error checking etc, needs struct or something. 
	if ( foo) (*foo)(NFUNC,inputs,0,NOUT,outputs);
	// same as with input 
	for (  ii=0; ii<NOUT; ++ii)
	{
		const int idx=osrc[ii];
		if ( idx>=0) colptrs[idx][row]=outputs[ii];
		//printf("outputting %f to %d from %f \n",colptrs[idx][row],idx,outputs[ii]);
	}
	if ( NOUT) REAL(df_out)[row]=outputs[0];
} // row

    UNPROTECT( 1 );
return df_out;

} //find_xc

// not usable now. 
#if 0  

int main ( int argc, char ** argv)
{ 

SEXP x; 
find_xc(x,0);



return 0; }

#endif

#ifdef __cplusplus 
}
#endif


#endif //  INCLUDE_R_INTERFACE

#endif // guard

