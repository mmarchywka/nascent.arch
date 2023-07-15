#ifndef MJM_TIMELINE_H__
#define MJM_TIMELINE_H__

#include "mjm_globals.h"
#include "mjm_timeline_logic.h"
// some of these pickup the math libs... 
#include "mjm_rational.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_interpolation.h"
#include "mjm_instruments.h"

#include <algorithm>
#include <map>
#include <string>
#include <fstream>
#include <signal.h>
/*


*/

/*
mjm_timeline.h : 
Create ASCII art, latex, or maybe svg for timelines derived from
line orienrted date event data. Needs linux date to turn lexi
date into day number

Charts need to be tractable so batch vs stream is ok.
Accumulate points as maps and convert to dense arrays.
*/

/*


g++ -DTEST_TIMELINE__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_timeline.h



2579  ./a.out -set-param "ordinal_column 0" -set-param "time_column 1"  -init -add biotin  -add-event biotin biotin -add B-6 -add-event B-6 B-6 -add B-100 -add-event B-100 B-100 -add seizure -add-event seizure seizure -harvest kelseyd -make -hlatex -unused -quit 2> fudd | grep -v "^  *$"
 2580  . srecord
 2581  history >> hist_doing
 2582  g++ -DTEST_TIMELINE__ -Wall  -std=gnu++11 -gdwarf-3 -I.. -Wno-unused-variable -Wno-unused-function -Wno-sign-compare -Wno-non-template-friend -x c++ mjm_timeline.h 
 2583  vi mjm_timeline.h
 2584  svn commit -m " it generates correct latex but looks bad "
 2585  history | grep event
 2586  ./a.out -set-param "ordinal_column 0" -set-param "time_column 1"  -init -add foo -add-event foo syrup -add Cu -add-event Cu Cu -add honk -add-event honk honking -harvest kelseyd -make -htext -banner -quit 2> fudd | grep -v "^  *$"
 2587  ./a.out -set-param "ordinal_column 0" -set-param "time_column 1"  -init -add foo -add-event foo syrup -add Cu -add-event Cu Cu -add honk -add-event honk honking -harvest events -make -hlatex -banner -quit 2> fudd | grep -v "^  *$"


*/

////////////////////////////////////////////////////////////////



class mjm_timeline 
{
/*




*/


class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef mjm_block_matrix<D> MyBlock;
typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef mjm_timeline Myt;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::MyBlock  MyBlock;
typedef Tr::MySparse MySparse;

typedef std::ifstream IfTy;

typedef TimelineBranches BraTy;
typedef TimelineParams FlpTy;


typedef mjm_block_matrix<char> MyPageBase;
class my_page_image : public MyPageBase
{

typedef my_page_image Myt;
typedef MyPageBase Super;
public:
my_page_image() : Super(){init(); }
my_page_image(const IdxTy n, const IdxTy m) : Super(n,m) {init(); }

void init()
{
blank();
}
void blank() { for (IdxTy i=0; i<size(); ++i) m_ptr[i]=' '; }
StrTy string_image() const 
{
Ss ss;
IdxTy cnt=0;
LoopItor2 itor(m_dims);
while (itor)
{
ss<<(*this)(itor.cursor());
//if (' '!=(*this)(itor.cursor()))
//MM_MSG(" fudd "<<MMPR(itor.cursor()[0])<<MMPR(itor.cursor()[1]))

++cnt;
++itor;
const bool eol=((!itor)|| (itor[m_top]==0)); // || (cnt==n));
if (eol) { ss<<CRLF; cnt=0; }
}

return ss.str();
}

}; // my_page_image

// really shoulc be a char array but may have more
// than that although it could be nested array 
// this is really a almost  page image
typedef mjm_block_matrix<StrTy> MyChartBase;
class my_chart : public MyChartBase
{
typedef my_chart Myt;
typedef MyChartBase Super;
public:
my_chart() : Super(), m_points(0), m_tl(0) {}
my_chart(const IdxTy n, const IdxTy m) : Super(n,m),m_points(m),m_tl(n) {}
// the order is backwards here 
void set_size(const IdxTy ntl, const IdxTy points)
{
m_points=points;
m_tl=ntl;
resize(ntl,points);
names.clear();
}
void add_name(const StrTy nm ) { names.push_back(nm); }
void add_time(const StrTy t ) { times.push_back(t); }
void add_time(const IdxTy t ) { ord_times.push_back(t); }

StrTy to_ssv(const bool hor)
{
Ss ss;
const bool have_names= (names.size()==m_tl);
const bool have_times= (times.size()==m_points);
const bool have_ord_times= (ord_times.size()==m_points);
// the headers don't work since times are fudded up
// for entries with no events fudd 
const bool top_header=!true;
const bool row_label=!true;
IdxTy row_labels=0;
const StrTy sep=" ";
if( top_header)
{
if ( !hor)
{
if (row_label) { row_labels=2; ss<<"time ord "; } 
if (names.size()<m_tl) 
	MM_ERR(" sizes wrong for names "<<m_tl<<" v "<<names.size())
for (IdxTy i=0; i<m_tl; ++i) {if (i!=0) {ss<<sep;}  ss<<names[i]; }
ss<<CRLF;
}else
{
if (times.size()<m_points) 
	MM_ERR(" sizes wrong for times "<<m_points<<" v "<<times.size())
if (ord_times.size()<m_points) 
	MM_ERR(" sizes wrong for ord "<<m_points<<" v "<<ord_times.size())

if (row_label) { row_labels=1; ss<<"timeline "; } 
for (IdxTy i=0; i<m_points; ++i) {if (i!=0) {ss<<sep;}  ss<<times[i]; }
ss<<CRLF;
for (IdxTy i=0; i<m_points; ++i) {if (i!=0) {ss<<sep;}  ss<<ord_times[i]; }
ss<<CRLF;
}
}

const IdxTy row_limit=hor?m_tl:m_points;
const IdxTy line_limit=hor?m_points:m_tl;
// good case for the iterator lol 
for (IdxTy i=0; i<row_limit; ++i)
{ 
if (row_label)
{
if (hor) { ss<<names[i]<<sep; }
else { ss<<times[i]<<sep; if (row_labels>1) {ss<<ord_times[i]<<sep; } }
} // row_label
for (IdxTy j=0; j<line_limit; ++j)
{
if (j!=0) ss<<sep;
// subs are tl, point
const StrTy et=hor?(*this)(i,j):(*this)(j,i);
if (et.length()!=0) ss<<et;
}
ss<<CRLF;
}

return ss.str();
}

StrTy to_ascii_h() const
{
Ss ss;
const IdxTy page_height=50;
const IdxTy page_width=80;
my_page_image pi(page_height,page_width);
pi.blank();
// layout the timelines
for (IdxTy i=0; i<m_tl; ++i)
{

for (IdxTy j=0; j<m_points; ++j)
{
const StrTy et=(*this)(i,j);
if (et.length()!=0)
{
//MM_MSG(" shot fudd place "<<MMPR(i)<<MMPR(j)<<MMPR(et))
IdxTy r=i;
IdxTy c=j;

if (c>page_width) 
{
MM_ONCE(" wrapping "<<MMPR(r)<<MMPR(page_height),)
IdxTy nadv=c/page_width;
c=c % page_width; 
r+=nadv*(m_tl+3);
}
if (r>page_height) 
{
MM_ONCE(" wrapping "<<MMPR(r)<<MMPR(page_height),)
r=r % page_height; 
}



 pi(2*r,c)=et.c_str()[0];
}

} 
}
ss<<pi.string_image();
return ss.str();
}

StrTy to_ascii_v() const
{
Ss ss;
return ss.str();
}


StrTy to_latex_h() const
{
if (m_points<1) return "";
const IdxTy ppseg=40;
Ss ss;
ss<<latex_start();
//ss<<latex_tab(m_points+2,7,6);
const IdxTy segs=((m_points-1)/ppseg)+1;
if ( segs!=1) { ss<<latex_tab(ppseg+2,7,1);}
else { ss<<latex_tab(7*((m_points+6)/7)+2,7,1); }
for (IdxTy k=0; k<segs; ++k )
{
ss<<" \\hline "<<CRLF;
const IdxTy segstart=k*ppseg;
IdxTy segend=segstart+ppseg;
ss<<latex_weeks(segstart/7,segend-segstart);
if (segend>m_points) segend=m_points;
for (IdxTy i=0; i<m_tl; ++i) {
if (i<names.size()) ss<<names[i];
ss<<"&";
for (IdxTy j=segstart; j<segend; ++j) {
const StrTy et=(*this)(i,j);
if (et.length()!=0) ss<<et.c_str()[0];

ss<<"&";
}
ss<<"\\\\";
}
} // segs
ss<<" "<<latex_end();

return ss.str();
}
StrTy to_latex_v() const
{
Ss ss;
return ss.str();
}

private:

StrTy latex_weeks(const IdxTy week0,const IdxTy len) const
{
Ss ss;
IdxTy weekend=week0+len/7; // this could cause problems with partialweeks
for (IdxTy i=week0; i<weekend; ++i){ 
ss<<"&";
ss<<"\\multicolumn{7}{|c|}{Week "<<i<<"}";  }
ss<<"\\\\";
return ss.str();
}

StrTy latex_start() const
{
Ss ss;
ss<<"\\begin{table}\\tiny\\setlength\\tabcolsep{0.5pt}\\begin{tabular}" ; 
return ss.str();
}
//text & text & text
StrTy latex_tab(const IdxTy times, const IdxTy delline,const IdxTy off) const
{
Ss ss;
ss<<"{";
for (IdxTy i=0; i<times; ++i)
{
if ((i%delline)==off) ss<<"|";
ss<<"c";
}
ss<<"}";
return ss.str();
}
StrTy latex_end() const
{Ss ss;
ss<<"\\end{tabular} \\end{table}";
return ss.str();
}



IdxTy m_points, m_tl;
std::vector<StrTy> names;
std::vector<StrTy> times;
std::vector<IdxTy> ord_times;


}; // my_chart
typedef my_chart MyChart;

// translate input raw time designations into uniform ones
// such as lexi date into day number
class time_map
{
public:
typedef StrTy DateTy;
typedef IdxTy OrdTy;
typedef std::map<DateTy, OrdTy> TMap;
time_map(): m_warn(true) {}

void add_par(const DateTy & date, const OrdTy ord)
{
if (m_warn) warn_differ(date,ord);
m_map[date]=ord;
}
bool warn_differ(const DateTy & date, const OrdTy ord)
{
auto ii=m_map.find(date);
if (ii!=m_map.end())
{
const OrdTy & old=(*ii).second;
if (old!=ord)
	MM_ERR( "rededfining date "<<date<< " from "<<old<<" to "<<ord)
return old==ord;
}
return true;
}
OrdTy operator[](const DateTy & d) { return m_map[d]; } 
private:
bool m_warn;
TMap m_map;

}; // time_map

class time_line
{
public:
typedef StrTy TimeTy;
typedef StrTy EventTy;
//typedef std::map<TimeTy, EventTy> TMap;
typedef std::map<StrTy, EventTy> EventMap;
typedef std::map<EventTy,IdxTy> EventCounts;
typedef std::map<StrTy,IdxTy> SymbolCounts;
typedef std::map<EventTy,StrTy> EventSymbolMap;
typedef std::map<StrTy, EventCounts> TMap;
time_line(const StrTy & nm) : m_name(nm),m_total_events(0),m_sum_numbers(false)  {}
time_line() : m_name(""),m_total_events(0),m_sum_numbers(false)  {}

void sum_numbers(const bool x=true) { m_sum_numbers=x; }
//void sum_numbers() { m_sum_numbers=true; }
bool  add(const TimeTy & t, const EventTy & v, const IdxTy col=~0)
{ 
// if this is non-zero then add to accumulator for thistime 
if (m_sum_numbers)
{
// this is only a partial sum based on operand, need
// to sum the map later
IdxTy val=atoi(v.c_str());
if (val!=0)
{
	// need to sum the map later doh 
	m_map[t][v]+=val;
return true; 
}

}
if (m_events.find(v)==m_events.end()) return false; 
	++(m_map[t][v]);
++m_total_events;
return true;
 }
// retrain however the old ones 
void clear_events()  { m_events.clear(); m_sum_numbers=false; }
void clear_record()  { m_map.clear();}

bool add_event(const StrTy & event) 
{
if (m_events.find(event)!=m_events.end()) return false; 
m_events[event]="x"; // nothing for now a set would do 
return true;
}

// now this needs a way to output the events...
const EventTy get(const TimeTy & t) const
{
auto ii=m_map.find(t);
if (ii== m_map.end()) return EventTy();
// map is non empty
return  (*(*ii).second.begin()).first;

}
const StrTy get_symbol(const TimeTy & t) const
{
if (m_sum_numbers)
{
Ss ss;
EventCounts ec=get_all(t);
IdxTy sum=0;
for (auto ii=ec.begin(); ii!=ec.end(); ++ii) 
{
	const IdxTy fac=atoi((*ii).first.c_str());
	sum+=(*ii).second*fac;
}
if (sum<10) { ss<<sum; } else if (sum<20) {ss<<"+"; } else {ss<<"!"; }
return ss.str();
}

const EventTy & e= get(t);
return event_lut(e);
}

IdxTy n_events() const { return m_total_events; }
// this only makes a copy due to f-ing default value 
const EventCounts  get_all(const TimeTy & t) const
{
auto ii=m_map.find(t);
if (ii== m_map.end()) return EventCounts();
// map is non empty
return  ((*ii).second);
}

const SymbolCounts  get_all_symbols(const TimeTy & t) const
{
const EventCounts & ec=get_all(t);
SymbolCounts sc;
for (auto ii=ec.begin(); ii!=ec.end(); ++ii)
{
EventTy key=(*ii).first;
IdxTy count=(*ii).second;
sc[event_lut(key)]=count;
}

return sc;
}
StrTy event_lut(const EventTy & e) const 
{
auto ii=m_sym.find(e);
if (ii==m_sym.end()) return StrTy(e);
return (*ii).second;


}



const IdxTy  n_types(const TimeTy & t) const
{
auto ii=m_map.find(t);
if (ii== m_map.end()) return 0;
// map is non empty
return  ((*ii).second).size();

}



	//	(*ii).second.get_all_times(all_times);	
template <class Ty> void get_all_times(Ty & all_times)
{

for ( auto ii=m_map.begin(); ii!=m_map.end(); ++ii)
	++all_times[(*ii).first]; 

}
const StrTy & name() const { return m_name; } 
StrTy to_string() const
{
Ss ss;
ss<<MMPR(m_name)<<MMPR(m_map.size())<<MMPR(m_events.size());

return ss.str();

}
StrTy m_name;
IdxTy m_total_events;
TMap m_map;
EventMap m_events;
EventSymbolMap m_sym;
bool m_sum_numbers;
}; // time_line
///////////////////////////////////////////////////////////////////




typedef std::map<StrTy,time_line> TimeLines;
typedef std::map<StrTy,IdxTy> KeyCounts;

public :
//ficklanc():m_size(1),m_save_sol(false) {Init();}
mjm_timeline() {Init();}
mjm_timeline(int argc,char **_args) // : m_save_sol(false)
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<argc; ++i) args[i]=_args[i];
for (IdxTy i=argc; i<ikluge; ++i) args[i]=&dummy[0];

//m_size=1;
//m_points=2000; // this really should not have a default.. 
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
m_tree.config("-tree",i,argc,args);
m_flp.config("-params",i,argc,args);
//configi(m_points,"-points",i,argc,args);
m_flp.config_set("-set-param",  i,  argc, args);
m_tree.config_set("-set-branch",  i,  argc, args);
cmdlcmd( i, argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
}
////////////////////////////////////////////////////////
// command block

// this should be in the parameters map, nothing special here... 
 void configi(IdxTy & dest, const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
//const StrTy nm=StrTy(args[i]);
dest=::atoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==StrTy("-harvest")) { ++i;
//dest=::atoi(args[i]);
const StrTy fn=StrTy(args[i]);
process_file(fn);
if (confirm) MM_ERR(" done harvesting "<<fn)
++i; // consume param and cmd
return;
}
if ( s==StrTy("-htext")) { ++i;
if(confirm) MM_ERR(" output as htext ") // <<s<<" to "<<dest)
StrTy s=m_chart.to_ascii_h();
std::cout<<s<<CRLF;
//++i; // consume param and cmd
return;
}


if ( s==StrTy("-vtext")) { ++i;
//if(confirm) MM_ERR(" setting "<<s<<" to "<<dest)
if(confirm) MM_ERR(" output as vtext ") // <<s<<" to "<<dest)
StrTy s=m_chart.to_ascii_v();
std::cout<<s<<CRLF;
//++i; // consume param and cmd
return;
}

if ( s==StrTy("-cmd")) { ++i;
const StrTy cmd=StrTy(args[i]);
if(confirm) MM_ERR(" command mode exec of  "<<cmd) // <<s<<" to "<<dest)
command_mode(cmd); 
++i; // consume param and cmd
return;
}

// tl 
if ( s==StrTy("-hlatex")) { arg_cmd(i,args,0,"hlatex",confirm); return; }
if ( s==StrTy("-vssv")) { arg_cmd(i,args,0,"vssv",confirm); return; }
if ( s==StrTy("-hssv")) { arg_cmd(i,args,0,"hssv",confirm); return; }
if ( s==StrTy("-init")) { arg_cmd(i,args,0,"init",confirm); return; }
if ( s==StrTy("-quit")) { arg_cmd(i,args,0,"quit",confirm); return; }
if ( s==StrTy("-add")) { arg_cmd(i,args,1,"add",confirm); return; }
// tl and event
if ( s==StrTy("-add-event")) { arg_cmd(i,args,2,"add-event",confirm); return; }
if ( s==StrTy("-add-sum")) { arg_cmd(i,args,1,"add-sum",confirm); return; }
if ( s==StrTy("-set-sum")) { arg_cmd(i,args,1,"set-sum",confirm); return; }
if ( s==StrTy("-reset-sum")) { arg_cmd(i,args,1,"reset-sum",confirm); return; }
if ( s==StrTy("-make")) { arg_cmd(i,args,0,"make",confirm); return; }
if ( s==StrTy("-unused")) { arg_cmd(i,args,0,"unused",confirm); return; }
if ( s==StrTy("-banner")) { arg_cmd(i,args,0,"banner",confirm); return; }

} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 


void dump_unused()
{
Ss ss;
for (auto ii=m_unused.begin(); ii!=m_unused.end(); ++ii)
{
ss<<(*ii).first<<" "<<(*ii).second<<CRLF;
}
MM_MSG("unused:"<<CRLF<<ss.str())

}


// take current time lines and compile the chart
void make_chart()
{
	const bool  dump_timeline_info=false;
// for now just go through each time line and
// put symbols on the chart, need to create
// table of symbol meanings. 
	//typedef StrTy Event;
	//const IdxTy points=0;
	// now make a list of times to plot and get
	// eac for each
	auto ix=m_timelines.begin();
	// create a map of times
	std::map<StrTy, IdxTy> all_times;
	for (auto ii=ix ; ii!=m_timelines.end(); ++ii)
	{
		(*ii).second.get_all_times(all_times);	
	}
	std::vector<StrTy> time_vector;
	for (auto ii=all_times.begin() ; ii!=all_times.end(); ++ii)
	{
		time_vector.push_back((*ii).first); 
	}
	std::sort(time_vector.begin(), time_vector.end());
	const IdxTy times=time_vector.size();
	if (times==0) return;
	const IdxTy tzero=m_time_map[time_vector[0]];
	const IdxTy tend=m_time_map[time_vector[times-1]];
// there is no ordinal time if no event for that date... 
// fudd 
	m_chart.set_size(m_timelines.size(), times);
	m_chart.set_size(m_timelines.size(), tend-tzero+1);
	// place string at the right places, leaving it to the formatter
	// to use table etc. 
	for (IdxTy i=0; i<times; ++i)
	{
			const StrTy & t=time_vector[i];
	const IdxTy tn=(m_time_map[time_vector[i]]);
m_chart.add_time(t);
m_chart.add_time(tn);

	}
	IdxTy j=0;	
	for (auto ii=ix ; ii!=m_timelines.end(); ++ii)
	{
		m_chart.add_name((*ii).second.name());
		for (IdxTy i=0; i<times; ++i)
		{
			const StrTy & t=time_vector[i];
			//const Event e=(*ii).second.get(t);
			const StrTy e=(*ii).second.get_symbol(t);
	const IdxTy tn=(m_time_map[time_vector[i]]);
			if (e.length()!=0)
				{ m_chart(j,tn-tzero)=e;

//MM_MSG(" adding fudder "<<MMPR(t)<<MMPR(e)<<MMPR(i)<<MMPR(j))
}
			// need to convert into location to mark 
		}
		++j;
	} // ii 
	if ( dump_timeline_info)
	{
	IdxTy iii=0;
	for (auto ii=ix ; ii!=m_timelines.end(); ++ii)
	{
		MM_MSG(" time line i="<<iii<<" "<<(*ii).second.to_string())
		++iii;	
	}
	}
}

void add_event(const StrTy & nm, const StrTy & e)
{
m_timelines[nm].add_event(e);
}
void set_sum(const StrTy & nm) { m_timelines[nm].sum_numbers(); }
void reset_sum(const StrTy & nm) { m_timelines[nm].sum_numbers(); }

void all_timelines( const IdxTy op)
{

for (auto ii=m_timelines.begin(); ii!=m_timelines.end(); ++ii)
{
switch (op)
{
case 0: { (*ii).second.clear_events(); break; } 
case 1: { (*ii).second.clear_record(); break; } 

default : MM_ERR(" bad timline op "<<op); 
} 

}


}

void add_timeline(const StrTy & nm)
{
m_timelines[nm]=time_line(nm);
}

void process_file( const StrTy & nm)
{
IfTy ifs(nm.c_str());
process_stream(ifs);

}
bool have_ord_col(const IdxTy sz) const 
{
if ((m_ord_col>=0) &&(m_ord_col<sz)) return true;
return false; 
}
bool ok_to_process(const IdxTy sz) const 
{
return sz>m_t_col; 
}
bool event_col(const IdxTy i, const IdxTy cola, const IdxTy sz) const 
{
return (i!=m_t_col)&&(i!=m_ord_col); 
}

template <class Is> void process_stream(Is & is)
{
	const bool warn_skip=true;
	const IdxTy n_tl= m_timelines.size();
	LineIterator li(&is);
	while (li.nextok())
	{
		const IdxTy lineno=li.lineno();
		const IdxTy sz=li.size(); // number of words on the line
//		const StrTy time=li.word(m_t_col); // time index
		// now harvest the events indicated
		if (!ok_to_process(sz)) 
		{
			MM_INC_MSG(m_cm,"skip_input_line" )
			if (warn_skip) {MM_ERR(" skipping line "<<li.lineno()<<" : "<<li.line()) } 
			continue;
		}
		//MM_MSG(" processing line "<<li.line()) 
		const StrTy &  t=li.word(m_t_col);
		if (have_ord_col(sz)) 
		{
			m_time_map.add_par(t, atoi(li.word(m_ord_col).c_str()));
		}
		// for the rest of the line, try to extract events
		// although this may be a mess lola
		IdxTy cola=0; // ampersands can be expanded internally  
		for (IdxTy i=0; i<sz; ++i)
		{
			// check for ampersands, then split as needed
			const StrTy word=li.word(i);
			if ( event_col(i,cola,sz))
			{
				// need rules for assigning events to timelines
				// and tracking the debris
					// it is up to the timeline to reject it 
				bool used=false;
				auto ii=m_timelines.begin();
				for(IdxTy k=0; k<n_tl; ++k)
				{ used|=(*ii).second.add(t,word,i); ++ii; }
				if (!used) { ++m_unused[word]; }


			}
			//++cola;
		} // i<sz word
	} // nextok

}
void command_mode()
{
//LineIterator 
CommandInterpretter li(&std::cin);
command_mode(li);

}
void command_mode(const StrTy & cmd)
{

CommandInterpretter li;
li.set(cmd,1);
//li.set(cmd,0);
command_mode(li);
}


void command_mode(CommandInterpretter & li)
{
StrTy local_label="fick";
while (li.nextok())
{
const IdxTy sz=li.size();
MM_ERR(" processing "<<li.dump())
if (sz<1) continue;
const StrTy cmd=li.word(0);
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
if (cmd=="get-tree") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_tree.get(li.word(1))<<CRLF;  continue; } 

if (cmd=="set-tree") { if (li.cmd_ok(3))  m_tree.set(li.cmd_set());  continue; } 
if (cmd=="set-label") { if (li.cmd_ok(2))  local_label=(li.word(1));  continue; } 
//if (cmd=="points") { if (li.cmd_ok(2))  m_points=::atoi((li.word(1).c_str()));  continue; } 
//if (cmd=="steps") { if (li.cmd_ok(2))  steps(li.n(1));  continue; } 
//if (cmd=="dump") { dump_state_etc(std::cout,local_label,1);  continue; } 
//if (cmd=="dump_to") { if (li.cmd_ok(2)) dump_to_file(li.word(1),local_label,1);  continue; } 
if (cmd=="process") { if (li.cmd_ok(2)) process_file(li.word(1));  continue; } 
if (cmd=="add") { if (li.cmd_ok(2)) add_timeline(li.word(1));  continue; } 
if (cmd=="init") { Init();  continue; } 
if (cmd=="hlatex") { std::cout<<m_chart.to_latex_h();  continue; } 
if (cmd=="hssv") { std::cout<<m_chart.to_ssv(true);  continue; } 
if (cmd=="vssv") { std::cout<<m_chart.to_ssv(!true);  continue; } 
if (cmd=="add-event") { if (li.cmd_ok(3)) add_event(li.word(1),li.word(2));  continue; } 
if (cmd=="add-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
if (cmd=="set-sum") { if (li.cmd_ok(2)) set_sum(li.word(1));  continue; } 
if (cmd=="reset-sum") { if (li.cmd_ok(2)) reset_sum(li.word(1));  continue; } 
if (cmd=="make") { make_chart();  continue; } 
if (cmd=="unused") { dump_unused();  continue; } 
if (cmd=="clear-unused") { m_unused.clear();  continue; } 
//if (cmd=="legend") { dump_legend(std::cout,local_label,1);  continue; } 
//if (cmd=="dump-human") { dump_state_etc(std::cout,local_label,0);  continue; } 
if (cmd=="init") { Init();   continue; } 
//if (cmd=="save-solution") { m_save_sol=true;  continue; } 
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
/*
if (cmd=="test") { test(li);  continue; } 
//if (cmd=="refine") { if (cmd_ok(li,sz,2))  refine(::atof((li.word(1).c_str())));  continue; } 
if (cmd=="refine") { if (li.cmd_ok(2))  refine(::atof((li.word(1).c_str())));  continue; } 
*/
if (cmd=="quit") { clean_up(); return; } 

MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())

}




} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. fudd 
li.push(); // shift commands and remove invoking one, 
const StrTy & cmd=li.word(0);
/*
if (cmd=="flow") { test_flow(li); li.pop(); return ; } 
if (cmd=="symm") { test_flow_symmetry(li,true); li.pop(); return ; } 
if (cmd=="symmplot") { test_flow_symmetry(li,!false,1.01); li.pop(); return ; } 
if (cmd=="driftsf") { test_driftsf(li); li.pop(); return ; } 
if (cmd=="diffsf") { test_diffsf(li); li.pop(); return ; } 
if (cmd=="diffsfmat") { test_diffsf_mat(li); li.pop(); return ; } 
if (cmd=="diffsfolap") { test_diffsf_olap(li); li.pop(); return ; } 

//void test_diffsf_conserve( CommandInterpretter & li )
if (cmd=="diffsfconserve") { test_diffsf_conserve(li); li.pop(); return ; } 
if (cmd=="diffsfbc") { test_diffsf_bc(li); li.pop(); return ; } 
if (cmd=="diffimpulse") { test_diff_impulse(li); li.pop(); return ; } 
if (cmd=="gr") { test_gr_integrator(li); li.pop(); return ; } 
*/
//void test_gr_integrator( CommandInterpretter & li )
MM_ERR(" unrecignized TEST command "<<li.dump())

li.pop();
} // test

void dump_cm()
{
 m_cm.dump("solve_step"  , std::cout);
 MM_MSG("solve_times"<<CRLF<<m_cm.time_table("solver_times"))
}

void config_banner()
{
MM_INC_MSG(m_cm,"test test" )
MM_MSG(" configuration "<<m_flp.to_string())
MM_MSG(" logic "<<m_tree.to_string())
MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}

bool done() const  { return m_done; } 

private:
void Init()
{
m_ord_col=m_flp.ordinal_column();
m_t_col=m_flp.time_column();
m_done=false;

}
bool m_done;
BraTy m_tree;
FlpTy m_flp; // right now this is being SET in Init not being used to set m_h

TimeLines m_timelines;
KeyCounts m_unused;

MyChart m_chart;
time_map m_time_map;
IdxTy m_ord_col; // input column containing a time ordinal (such as day number)
IdxTy m_t_col; // input col with absolute time string such as date.


CounterMap m_cm;

}; //mjm_timeline 



/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////

#ifdef  TEST_TIMELINE__

int main(int argc,char **args)
{
typedef mjm_timeline Myt;

Myt x(argc,args);

if (!x.done()) x.command_mode();
return 0;
}

#endif // TEST_FICK__

#endif

