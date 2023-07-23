#ifndef MJM_STRING_SEQ_H__
#define MJM_STRING_SEQ_H__
 
#include "mjm_globals.h"
// for the presence absence vector 
#include "mjm_char_mat.h"
#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_canned_methods.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
//#include "../mjm_fasta_ii.h"
#include "mjm_fasta_ii.h"

// for delegating processing to s script... 
#include "mjm_pawnoff.h" 
#include "mjm_wscat_bot.h" 


#include "mjm_collections.h"
#include "mjm_svg_writer.h"

#ifdef PYTHON_BOOST_BUILD
#include <boost/python.hpp>
#endif


//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
//#include "mjm_biom_hdf5.h"
#include <algorithm>
#include <random>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>
#include <stdint.h>

#include <pthread.h>

#include <sys/time.h>

int my_random_seed() 
{
struct timeval tp;
gettimeofday(&tp, NULL);
typedef long long Tl;
Tl  t = ((Tl) tp.tv_sec)  * 1000000 + tp.tv_usec; //  / 1000;
if (tp.tv_usec==0) t=tp.tv_sec;
return int(t);
}

/*



2073  ./mjm_string_seq.out -cmd "set-param mflags 2"  -cmd "read-fasta knowns limis_k2.fasta" -cmd "index-fasta knowns 8"  -cmd "set-param maxscores 15" -cmd "set-param save_aligns 1"  -cmd "mt-explore-streaming-hits knowns ku.fasta" -cmd "cmd-align default \"set 45 red\"" -cmd "write-align xxx.svg default"  -quit # | tee linear_combs_cross._txt




This is probably the code used to generate the better histograms but there is no assurance that xxx5
is the same. Need to get teh svn version for that date too,

Today is April 30 but this seems like forever ago lol, 
ls -al ~/d/latex/keep/ *exrandomrat*
-rw-rw-r-- 1 marchywka marchywka 44682 2018-04-08 16:26 /home/marchywka/d/latex/keep/mjmexrandomrat.pdf

hist_doing7- 2740  ./mjm_string_seq.out -cmd "set-param mflags 2" -source xxx5 -quit  | grep bests | tee some_bests2r_random
hist_doing7- 2741  cat  some_bests2r_random | mjm eq | awk '{print $4" "$6" "$10" " $15" "$18" "$19}' > xxx
hist_doing7- 2742  cat xxx | awk '{print $5}' | sort | uniq -c | sort -g  | tail
hist_doing7- 2743  cat  some_bests2r_random |more 
hist_doing7- 2744  more xxx
hist_doing7- 2745  cat  some_bests2r_random | mjm eq | awk '{print $4" "$6" "$10" " $14$15$16" "$20" "$21}' > xxx
hist_doing7- 2746  more xxx
hist_doing7- 2747  cat  some_bests2r_random | mjm eq | awk '{print $4" "$6" "$10" " $14$15$16" "$19" "$20}' > xxx
hist_doing7- 2748  more xxx
hist_doing7- 2749  cat xxx | awk '{print $5}' | sort | uniq -c | sort -g  | tail
hist_doing7- 2750  cat xxx | awk '{print $5" "$6}' | sort | uniq -c | sort -g  | tail
hist_doing7- 2751  history
hist_doing7- 2752  cat xxx | awk '{print $5" "$6}' | sort | uniq -c | sort -g  | tail
hist_doing7- 2753  R
hist_doing7- 2754  mv mjmexrandomhisto.pdf ~/d/latex/keep 
hist_doing7: 2755  mv mjmexrandomrat.pdf ~/d/latex/keep 

 more xxx5
read-fasta nohits /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
read-fasta repset /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/otus/rep_set.fna
read-fasta nohitdemo /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
#make-fasta random random 100 1000
make-fasta random random 300 500
read-fasta knowns ncbitax_16S.fasta 
read-fasta misc misc.fasta 
#read-fasta q mjm_commons16s.fasta 
read-fasta q mjm_10_16s.fasta 
index-fasta knowns 8
#read-fasta q best_discriminators.fasta 
#read-fasta q best_16s_signatures.fasta
#read-fasta q sigs10.fasta 
#fasta-commons q samples 
#lowp-fasta-hits q samples 10
#explore-hits knowns samples 10
#explore-hits knowns nohits 10
#explore-hits knowns nohits
#explore-hits knowns repset 
#explore-hits knowns nohits 
#explore-hits knowns random 
#explore-hits knowns nohitdemo
explore-hits knowns misc
#explore-hits random random 10
#explore-hits samples samples 10
#fasta-summary q samples 






*/


/*

 3422  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus AAAAAAAA" -quit | tee fick
 3423  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$NF" "$3; print $2}' > seed8_pairs1.fasta
 3424  cp fick seed8_pairs1._txt
 3434  #cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$NF" "$3; print $2}' > seed8_pairs1.fasta
 3435  tail fick
 3436  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs1.fasta
 3443  tail fick
 3444  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus AAGCCGGT" -quit | tee fick
 3445  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs2.fasta
 3446  cp fick seed8_pairs2._txt
 3466  tail fick
 3467  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus CCGTGAGG" -quit | tee fick
 3468  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs3.fasta
 3469  cp fick seed8_pairs3._txt
 3496  tail fick
 3497  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus CTGATGTG" -quit | tee fick
 3498  cp fick seed8_pairs4._txt
 3499  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs4.fasta
 3524  tail fick
 3525  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus GGTAAGGT" -quit | tee fick
 3526  cp fick seed8_pairs5._txt
 3528  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs5.fasta
 3545  tail fick
 3546  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-digenus TCCGCCTG" -quit | tee fick
 3547  cp fick seed8_pairs6._txt
 3550  cat fick | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$7" "$8" "$3; print $2}' > seed8_pairs6.fasta


*/
/*
901  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-ngenus AGGATTAG 0 .8 10" -quit |tee fick
 3902  cp fick seed8_ngenus8_10_part2._txt
 3903  bzip2 seed8_ngenus8_10_part2._txt
 3904  cp seed8_ngenus8_10_part2.fasta.bz2  `backupdir`
 3905  cp seed8_ngenus8_10_part2.bz2  `backupdir`
 3906  cp seed8_ngenus8_10_part2._txt.bz2  `backupdir`
 3907  cat fick | sed -e 's/ len=/ /' | awk 'BEGIN{cnt=0;}{cnt=cnt+1; i=7; x=""; while (i<=NF) { x=x" "$i; i=i+1;} print ">genuseek"cnt" "$3" "x; print $2}' > seed8_ngenus8_10_part2.fasta
 3908  head seed8_ngenus8_10_part2.fasta 
 3909  bunzip2 -c seed8_ngenus8_10_part*._txt.bz2 | sed -e 's/ len=/ /' | awk 'BEGIN{cnt=0;}{cnt=cnt+1; i=7; x=""; while (i<=NF) { x=x" "$i; i=i+1;} print ">genuseek"cnt" "$3" "x; print $2}' > seed8_ngenus8_10c.fasta
 3910  more seed8_ngenus8_10c.fasta 
 3911  ls -al *.fasta
 3912  df
 3913  history
marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ bunzip2 -c seed8_ngenus8_10_part*._txt.bz2 | sort | uniq | sed -e 's/ len=/ /' | awk 'BEGIN{cnt=0;}{cnt=cnt+1; i=7; x=""; while (i<=NF) { x=x" "$i; i=i+1;} print ">genuseek"cnt" "$3""x; print $2}' > seed8_ngenus8_10c.fasta
*/


/*
more xxx | awk 'BEGIN{cnt=0;}{cnt=cnt+1; print ">genuseek"cnt" "$NF" "$3; print $2}' > seed9.fasta

3305  ./mjm_string_seq.out  -cmd "read-fasta x seed9.fasta"  -cmd "read-fasta fick /home/marchywka/d/zymo/run1/zr2097.180309.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta" -cmd "index-fasta fick 8" -cmd "assfick x " -quit | awk '{s=$(NF-1); n=$(NF-3)" "$(NF-2); if (s==l) {t=t" "n; }else {print s" "t; t=n;}l=s}END{print s" "t;}' | grep seq686
 3306  history >> hist_doing7
 3307  . srecord
 3308  ./mjm_string_seq.out  -cmd "read-fasta x ncbitax_16S.fasta" -cmd "index-fasta x 8" -cmd "grow-genus AAAAAAAA" -quit | tee xxx
 3309  cp xxx seed8_part1.txt



cat z20.tex 
read-fasta knowns ncbitax_16S.fasta
read-fasta commons mjm_commons16s.fasta
index-fasta knowns 4
seq-hits AGCGTTGCTGAA
seq-hits AAAAAACCAATCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCG
distro-commons commons



cat xxx3
#read-fasta samples /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
#read-fasta samples /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta

#read-fasta samples ncbitax_16S.fasta 
#read-fasta q mjm_commons16s.fasta 
#read-fasta q sigs10.fasta 
read-fasta q best_16s_signatures.fasta 
enumerate-fasta dfasta q 
write-fasta foofooenu.fasta dfasta 
#read-fasta q best_16s_signatures.fasta
#read-fasta q sigs10.fasta 
#fasta-commons q samples 
#fasta-hits q samples 
#fasta-summary q samples 
 cat xxx4
read-fasta samples /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
#read-fasta samples /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/otus/rep_set.fna
#read-fasta samples /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta

read-fasta knowns ncbitax_16S.fasta 
#read-fasta q mjm_commons16s.fasta 
read-fasta q mjm_10_16s.fasta 
#read-fasta q best_discriminators.fasta 
#read-fasta q best_16s_signatures.fasta
#read-fasta q sigs10.fasta 
#fasta-commons q samples 
lowp-fasta-hits q samples 10
#fasta-summary q samples 

marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ cat xxx2
#read-fasta samples /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
read-fasta samples /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta

#read-fasta samples ncbitax_16S.fasta 
read-fasta q mjm_commons16s.fasta 
#read-fasta q best_discriminators.fasta 
#read-fasta q best_16s_signatures.fasta
#read-fasta q sigs10.fasta 
#fasta-commons q samples 
#fasta-hits q samples 
fasta-summary q samples 

marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ cat xxx1
#read-fasta ncbi /home/marchywka/d/latex/all16s.fasta
read-fasta sigs sigs10.fasta 
read-fasta q /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
 
#vector-search q zymo 
#hit-or-miss q ncbi 
###load-tax tax-info
err checking coverage of the tax tree with the fasta file 
###coverage lut ncbi
err check queries now  
#query-coverage q ncbi lut 1 5
#distinctives q ncbi
fasta_hits sigs q
#descend fusobacteria ncbi 
#coverage ncbi
marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ cat z13.txt 
read-fasta ku cclus.fasta
#read-fasta ku asdf.fasta
align-fasta ku
write-fasta kua.aln ku
write-interleaved-fasta ku kuai.aln ku
marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ 







TODO FIXME
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 



////////////////////////////////////////////////////////////////

class string_seq_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
string_seq_params( const StrTy & nm) : Super(nm) {}
string_seq_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//IdxTy mesh_nx() const { return m_map.get_uint("mesh_nx",10); } // // 100;
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy start_date() const { return m_map.get_string("start_date","2017-04-22"); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy otu_format() const { return m_map.get_uint("otu_format",0); } // // 100;
IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
IdxTy maxdepth() const { return m_map.get_uint("maxdepth",3); } // // 100;
bool print_counts() const { return m_map.get_bool("print_counts",!true); }
bool print_haves() const { return m_map.get_bool("print_haves",!true); }
bool print_havenots() const { return m_map.get_bool("print_havenots",!true); }
bool print_if_have() const { return m_map.get_bool("print_if_have",!true); }
bool suppress_vector() const { return m_map.get_bool("suppress_vector",!true); }
bool add_level() const { return m_map.get_bool("add_level",true); }
bool print_hit() const { return m_map.get_bool("print_hit",true); }
StrTy knowns_fasta() const { return m_map.get_string("knowns_fasta","knowns"); }
StrTy fasta_dir() const { return m_map.get_string("fasta_dir",""); }
StrTy hit_file() const { return m_map.get_string("hit_file",""); }
IdxTy mflags() const { return m_map.get_uint("mflags",3); } // // 100;
IdxTy maxscores() const { return m_map.get_uint("maxscores",8); } // // 100;
bool write_align_marks() const { return m_map.get_bool("write_align_marks",!true); }
bool save_aligns() const { return m_map.get_bool("save_aligns",!true); }
StrTy align_name() const { return m_map.get_string("align_name","default"); }
IdxTy random_interval() const { return m_map.get_uint("random_interval",0); } // // 100;
IdxTy n_threads() const { return m_map.get_uint("n_threads",3); } // // 100;
IdxTy refpos() const { return m_map.get_uint("refpos",0); } // // 100;
//StrTy catagory_map() const { return m_map.get_string("catagory_map","."); }
//StrTy sample_filter() const { return m_map.get_string("sample_filter","."); }
//bool  use_sample_filter() const { return m_map.get_bool("use_sample_filter",false); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;

ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"accmode="<<accmode()<<sep;
ss<<"maxdepth="<<maxdepth()<<sep;
ss<<"print_counts="<<print_counts()<<sep;
ss<<"print_haves="<<print_haves()<<sep;
ss<<"print_havenots="<<print_havenots()<<sep;
ss<<"print_if_have="<<print_if_have()<<sep;
ss<<"suppress_vector="<<suppress_vector()<<sep;
ss<<"add_level="<<add_level()<<sep;
ss<<"print_hit="<<print_hit()<<sep;
ss<<"knowns_fasta="<<knowns_fasta()<<sep;
ss<<"fasta_dir="<<fasta_dir()<<sep;
ss<<"hit_file="<<hit_file()<<sep;
ss<<"mflags="<<mflags()<<sep;
ss<<"maxscores="<<maxscores()<<sep;
ss<<"write_align_marks="<<write_align_marks()<<sep;
ss<<"save_aligns="<<save_aligns()<<sep;
ss<<"align_name="<<align_name()<<sep;
ss<<"random_interval="<<random_interval()<<sep;
ss<<"n_threads"<<n_threads()<<sep;
ss<<"refpos"<<refpos()<<sep;
//ss<<"otu_format="<<otu_format()<<sep;
//ss<<"catagory_map="<<catagory_map()<<sep;
//ss<<"sample_filter="<<sample_filter()<<sep;
//ss<<"use_sample_filter="<<use_sample_filter()<<sep;
return ss.str();
}


}; // trees_and_tables_params


// reserved words or specific things with meanings different
// from canonical names .



namespace string_seq_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
class fake_int
{
typedef fake_int Myt;
typedef unsigned char Dr;
enum {SZ=32, SBM=1<<((sizeof(Dr)<<3)-1) };
public:
fake_int() { Init(); }
fake_int(const IdxTy x ) { Init(); Set(x); }
fake_int(const int x ) { Init(); Set(x); }
void Init() { memset(&m_a[0],0,sizeof(m_a));}
operator unsigned char () const { return m_a[0]; }
// this compiles but all the shift operators are wrong 
void left_shift(Myt & x) const 
{
for(IdxTy i=SZ-1; i>0; --i) 
{ 
Dr & t=x.m_a[i];
Dr & b=x.m_a[i-1];
t=t<<1;
if (b&SBM) t+=1;
} // i 
x.m_a[0]<<=1;
} // left_shift
void right_shift(Myt & x) const 
{
for(IdxTy i=0; i<(SZ-1); ++i) 
{ 
Dr & t=x.m_a[i];
Dr & b=x.m_a[i+1];
t=t>>1;
if (b&1) t+=SBM;
} // i 
x.m_a[SZ-1]>>=1;
} // left_shift



//Myt operator>>(const IdxTy n) const { Myt x=(*this); for(IdxTy i=1; i<SZ; ++i) { x.m_a[SZ-i-1]=x.m_a[SZ-i]; } x.m_a[SZ-1]=0; return x;  }// >>
Myt operator>>(const IdxTy n) const { Myt x=(*this); for(IdxTy i=0; i<n; ++i) right_shift(x); return x;  }// >>
Myt operator<<(const IdxTy n) const { Myt x=(*this); for(IdxTy i=0; i<n; ++i) left_shift(x); return x;  }// >>
Myt operator<<(const int n) const { Myt x=(*this); for(IdxTy i=0; i<n; ++i) left_shift(x); return x;  }// >>
//Myt operator<<(const IdxTy n) const { Myt x=(*this); for(IdxTy i=SZ-1; i<SZ; --i) { x.m_a[i]=x.m_a[i-1]; } x.m_a[0]=0; return x;  }// >>
//Myt operator<<(const int  n) const { Myt x=(*this); for(IdxTy i=SZ-1; i<SZ; --i) { x.m_a[i]=x.m_a[i-1]; } x.m_a[0]=0; return x;  }// >>
// check order.. 
Myt operator|(const Myt& that ) const { Myt x=(*this); for(IdxTy i=0; i<SZ; ++i) { x.m_a[i]|=that.m_a[i]; } return x;  }// >>
Myt operator&(const Myt& that ) const { Myt x=(*this); for(IdxTy i=0; i<SZ; ++i) { x.m_a[i]&=that.m_a[i]; } return x;  }// >>
Myt operator=(const Myt& that ) {  for(IdxTy i=0; i<SZ; ++i) { m_a[i]=that.m_a[i]; } return *this;  }// >>
bool operator>(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) 
{ int x=(int(m_a[i])&255)-(int(that.m_a[i])&255) ; if (x>0) return !false; if (x<0) return !true; }  
return false; 
}
bool operator>=(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) 
{ int x=(int(m_a[i])&255)-(int(that.m_a[i])&255) ; if (x>0) return !false; if (x<0) return !true; }  
return !false; 
}



bool operator<(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) 
{ int x=(int(m_a[i])&255)-(int(that.m_a[i])&255) ; if (x>0) return false; if (x<0) return true; }  
return false; 
}

bool operator<=(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) 
{ int x=(int(m_a[i])&255)-(int(that.m_a[i])&255) ; if (x>0) return false; if (x<0) return true; }  
return true; 
}
bool operator==(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) { if (m_a[i]!=that.m_a[i]) return false; }  
return true; 
}
bool operator!=(const Myt & that ) const
{
for(IdxTy i=(SZ-1); i<SZ; --i) { if (m_a[i]!=that.m_a[i]) return !false; }  
return !true; 
}




void Set(const IdxTy x) { for (IdxTy i=0; i<sizeof(IdxTy); ++i) m_a[i]=(x>>(i<<3))&255; } 
StrTy string() const
{
Ss ss;
ss<<std::hex<<std::setprecision(2);
for(IdxTy i=0; i<SZ; ++i ) ss<<m_a[i];
return ss.str();
}
Dr m_a[SZ];

}; // fake_int
//typedef unsigned int KeyCode;
//typedef unsigned __int128 KeyCode;
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
///typedef uint64_t KeyCode;
typedef fake_int  KeyCode;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
typedef std::vector<IdxTy > Locations;



//typedef mjm_sparse_matrix<D> MySparse;
}; // 

std::ostream & operator<<(std::ostream & os, const Tr::fake_int & f)
{
os<<f.string();
return os; 
} 



}; // trees_and_tables_traits
class mjm_string_query_group
{
typedef  string_seq_traits::Tr  Tr;
typedef mjm_string_query_group Myt;

typedef mjm_cli_ui<Myt> CliTy;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef std::vector<StrTy> Strings;
public:

void push_back(const StrTy & x) { m_s.push_back(x);}
const StrTy & operator[](const IdxTy i )const  { return m_s[i]; } 
StrTy & operator[](const IdxTy i ) { return m_s[i]; } 
IdxTy size() const { return m_s.size(); } 

Strings m_s;


}; // mjm_string_query_group


class mjm_index_bin_search
{
typedef  string_seq_traits::Tr  Tr;
typedef mjm_index_bin_search  Myt;

typedef mjm_cli_ui<Myt> CliTy;

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef std::vector<StrTy> Strings;
public:


const IdxTy& bad() const { static IdxTy i=~0; return i; } 
//template <class Ti > IdxTy find(const IdxTy & t, Ti & ii, const Ti & ee)
template <class Ti, class Tv > IdxTy find(const Ti & t, Tv & v, const IdxTy & sz) const
{
const IdxTy top=sz-1;
IdxTy ll=0;
IdxTy ul=top;
while (true) {
// TODO FIXME also provide linear data based interpolation etc 
IdxTy guess=(ll+ul)>>1;
Ti res=v[guess];
if (res==t) return guess;
if (ul<=ll) return bad();
if ( res>t) { 
if (guess==0) return bad(); 
ul=guess-1; 
} 
else ll=guess+1;
}
return bad();
}

template <class Ti>  
IdxTy divb(const IdxTy d, const Ti & t, const Ti & tlo, const Ti & thi) const
{
IdxTy bits=0;
IdxTy bit=1;
Ti p=t-tlo;
Ti q=thi-tlo;
if (q<10) { return (d>>1) ; } 
//if (true) { return (d>>1) ; } 
while (q>1)
{
q=q>>1;
if (p>=q) {bits|=(1<<bit); p=p-(p>>1); }
++bit; 
if (bit>8) break;
}  // true
IdxTy r=0;
for (IdxTy i=1; i<bit; ++i) { if ((bits&(1<<i))!=0) r=r+(d>>i); } 
//if (r==0) 
r+=(d>>bit);
//MM_ERR(MMPR3(r,bit,bits)<<MMPR4(d,t,tlo,thi))
return r;
} // divb

template <class Ti> IdxTy 
interp(const IdxTy d, const Ti & t, const Ti & tlo, const Ti & thi) const
{
if (thi<=tlo) return 0; 
if (t<=tlo) return 0; 
if (t>=thi) return d; 
if (true) return  divb( d, t, tlo,  thi);
// this could overflow 
const D uck = (D(t-tlo)/D(thi-tlo));
// D  fick= ((D(t)-D(tlo))/D(thi-tlo));
//if (fick>1) return d;
//if (fick<0) return 0; 
return IdxTy( d*uck); 
}
template <class Ti, class Tv > IdxTy find2(const Ti & t, const IdxTy & llast, Tv & v, const IdxTy & sz) const
{
const IdxTy top=sz-1;
IdxTy ll= llast;
IdxTy ul=top;
if (sz<1)
{
MM_ONCE(" size is zero wtf "<<MMPR3(t,llast,sz),)
return bad(); 
}
// only corner cases but otherwise requires check per iter... 
Ti  tlo=v[ll];
Ti thi=v[ul];
if (t>thi) return bad();
if (t<tlo) return bad();
while (true) {
// TODO FIXME also provide linear data based interpolation etc 
IdxTy guess=(ll+ul)>>1;
if (ul<ll) return bad();
//IdxTy guess=ll+interp(ul-ll,t,tlo,thi);
Ti res=v[guess];
if (res==t) return guess;
if (ul<=ll) return bad();
if ( res>t) { 
if (guess==0) return bad(); 
ul=guess-1; 
} 
else ll=guess+1;

tlo=v[ll];
thi=v[ul];
}
return bad();
}





/*
void push_back(const StrTy & x) { m_s.push_back(x);}
const StrTy & operator[](const IdxTy i )const  { return m_s[i]; } 
StrTy & operator[](const IdxTy i ) { return m_s[i]; } 
IdxTy size() const { return m_s.size(); } 

Strings m_s;
*/

}; // mjm_index_bin_search

class mjm_compact_string_index
{
typedef string_seq_traits::Tr  Tr;
typedef mjm_compact_string_index Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef Tr::KeyCode KeyCode;
typedef Tr::Locations Locations;

typedef  mjm_index_bin_search Bs;
public:
mjm_compact_string_index() 
: m_keys(0),m_n_keys(0),m_locations(0),m_n_locations(0)
,m_data(0),m_n_data(0),m_last_keycode(0),m_last_idx(0) {}
mjm_compact_string_index(const Myt & x) 
: m_keys(0),m_n_keys(0),m_locations(0),m_n_locations(0)
,m_data(0),m_n_data(0),m_last_keycode(0),m_last_idx(0)  {Copy(x);}

Myt& operator=(const Myt & that ) 
{
Copy(that);
return *this;
}

~mjm_compact_string_index() { Release(); }

const IdxTy&  bad() const {static  IdxTy x=~0;  return x; } 

const IdxTy & size() const { return m_n_keys;}
const KeyCode * keys() const { return m_keys;}
const KeyCode & key(const IdxTy& i ) const { return m_keys[i];}

//template <class Tv > IdxTy find(const IdxTy & t, Tv & v, const IdxTy & sz)
IdxTy count( const KeyCode & key)
{
IdxTy loc=m_bs.find(key,m_keys,m_n_keys);
if (loc==m_bs.bad())  return bad(); // TODO FIXME zero???? 
return m_data[m_locations[loc]];
}
class fake_vector
{
public:
fake_vector(const IdxTy * p, const IdxTy & n ) : m_size(n),m_ptr(p) {

//MM_ERR(MMPR2(m_size,m_ptr))
}
const IdxTy size() const { return m_size; }
const IdxTy& operator[](const IdxTy i ) const { return m_ptr[i]; }
// IdxTy& operator[](const IdxTy i )  { return m_ptr[i]; }

const IdxTy m_size;
const IdxTy * m_ptr;

}; // fake_vector

fake_vector find(const KeyCode & key) const
{
if (key<m_last_keycode) { m_last_idx=0; m_last_keycode=0; } 
// test kluge FIXME 2022
//if (true) { m_last_idx=0; m_last_keycode=0; } 
IdxTy loc=m_bs.find2(key,m_last_idx,m_keys,m_n_keys);
Ss sshex; sshex<<std::hex<<MMPR3(key,m_last_keycode,m_keys[0]);
//MM_ERR(MMPR(__FUNCTION__)<<MMPR3(sshex.str(),loc,m_last_idx))
//IdxTy loc=m_bs.find(key,m_keys,m_n_keys);
if (loc==m_bs.bad())  return fake_vector(0,0); // TODO FIXME zero???? 
m_last_keycode=key;
m_last_idx=loc;
const IdxTy base=m_locations[loc];
return fake_vector(m_data+base+1,m_data[base]); 
}
fake_vector find_index(const IdxTy & loc) const 
{
//IdxTy loc=m_bs.find(key,m_keys,m_n_keys);
if (loc==m_bs.bad())  return fake_vector(0,0); // TODO FIXME zero???? 
const IdxTy base=m_locations[loc];
//MM_ERR(MMPR4(base,loc,m_n_locations,m_n_data))
return fake_vector(m_data+base+1,m_data[base]); 
}



// could use DMEL lol 
template <class Tm> void load(const Tm & m)
{
Release();
KeyCode last=0;
const IdxTy sz=m.size();
m_keys= new KeyCode[sz]; m_n_keys=sz;
m_locations= new IdxTy[sz]; m_n_locations=sz;
IdxTy ds=0;
MM_LOOP(ii,m) { ds+=1+(*ii).second.size(); }
m_data=new IdxTy[ds]; m_n_data=ds; 
IdxTy i=0;
IdxTy data_ptr=0;
MM_LOOP(ii,m)
{
const KeyCode & k=(*ii).first;
//if (k<=last) { MM_ERR(" out of order "<<MMPR2(k,last)) } 
if (k<=last) { MM_ERR(" out of order ") } 
last=k;
m_keys[i]=k;
m_locations[i]=data_ptr;
const Locations & loc=(*ii).second;
m_data[data_ptr]=loc.size(); ++data_ptr;
MM_LOOP(jj,loc)
{
m_data[data_ptr]=(*jj);
++data_ptr;
} // jj 

++i;
} // ii 
if (m_n_data!=data_ptr) 
{ MM_ERR("DANGER WILL ROBINSON FACK "<<MMPR4(m_n_data,data_ptr,m_n_keys,m_n_locations)) } 

} // load

private:
void Release() {
delete [] m_keys; m_keys=0;
delete [] m_locations; m_locations=0;
delete [] m_data; m_data=0;
}
void Copy(const Myt & x) 
{
Release();
m_n_keys=x.m_n_keys;
m_keys=new KeyCode[m_n_keys];
m_n_locations=x.m_n_locations;
m_locations= new IdxTy[m_n_locations];
m_n_data=x.m_n_data;
m_data= new IdxTy[m_n_data];

::memcpy(m_keys,x.m_keys, m_n_keys*sizeof(KeyCode));
::memcpy(m_locations,x.m_locations, m_n_locations*sizeof(IdxTy));
::memcpy(m_data,x.m_data, m_n_data*sizeof(IdxTy));
}

// this is an array of KEYS in order
// the location of the kye points to index in m_locations
// which contains count folloed by n-values 
KeyCode * m_keys;
IdxTy m_n_keys;
IdxTy * m_locations;
IdxTy m_n_locations;
IdxTy * m_data;
IdxTy m_n_data;
Bs m_bs;
mutable KeyCode m_last_keycode;
mutable IdxTy m_last_idx;

}; // mjm_compact_string_index



class string_indexer 
{


typedef string_seq_traits::Tr  Tr;
typedef string_indexer Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef Tr::KeyCode KeyCode;
typedef Tr::Locations Locations;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
//typedef std::vector<IdxTy > Locations;
//typedef unsigned __int128 KeyCode;
//typedef unsigned int KeyCode;
typedef std::map<KeyCode, Locations> Sindex;
typedef  mjm_compact_string_index Cindex;
//typedef std::vector<StrTy> Strings;
public:

string_indexer(const char * c) :m_isize(4), m_s(c),m_compact(false) {Init();}
string_indexer(const char * c, const IdxTy & n ) :m_isize(n), m_s(c),m_compact(false)  {Init();}
void compact()
{
m_ci.load(m_index);
m_index.clear();
m_compact=true;
}

IdxTy bad() const { return ~0; } 
IdxTy find_first(const char * p ) const
{
if (m_compact) return find_first2(p);
return find_first1(p);
}
IdxTy find_first_agct(const char * p ) const
{
if (m_compact) return find_first_agct2(p);
return find_first1(p);
}




IdxTy find_first1(const char * p ) const
{
IdxTy i=0; 
const KeyCode mask=Mask();
typedef  mjm_index_bin_search Bs;
Bs bs; 
KeyCode key=Key(p+i);
//if (key==0) break; // return bad ???
auto li=m_index.find(key); 
// TODO should handle <4 better than this.. 
if ( li==m_index.end())  return bad(); 
const Locations & loc= (*li).second; // m_index[key];
MM_LOOP(ii,loc)
{
i=0;
const IdxTy locii=(*ii); // .first;
while (true){
i+=m_isize;
key=Key(p+i);
if (key==KeyCode(0)) return locii;
if (key<mask) // (1<<24))
{
bool partial_ok=true;
for (IdxTy j=0; j<m_isize; ++j)
{
if ( p[j+i]== 0 ) return locii;
if ( p[j+i]!=m_s[locii+j+i] ) { partial_ok=false; break;  }

} // j 
if (partial_ok) return locii; else break;
} // key 
li=m_index.find(key); 
if ( li==m_index.end()) break; //  return bad(); 
// jjjj
MM_ONCE(" fix or check  the index follow on lol",)
const IdxTy loc2=bs.find(locii+i,(*li).second,(*li).second.size());
if (loc2==bs.bad()) break; 
// instead of this, 
// if ( ((*li).second).find(locii+i)==((*li).second).end()) break;
// the vectors are ordered as per discovery algorightm, the
// matching 


} // true; 
//const Locations & locnext= m_index[key];
//if (locanext.find(locii+4)==locnext.end()) return bad(); } 
 
} // ii 
return bad(); 
}
////////////////////////////////////////////

IdxTy find_first2(const char * p ) const
{
MM_ONCE(" untested lol", ) 
IdxTy i=0; 
//const KeyCode mask=Mask();
//typedef  mjm_index_bin_search Bs;
//Bs bs; 
//MM_MSG(" ASSFACKL "<<MMPR3(strlen(p),i,p))
KeyCode key=Key(p+i);
//if (key==0) break; // return bad ???
//MM_MSG(" ASSFACKL "<<MMPR4(strlen(p),key,i,p))
if  (key<Mask())
{
MM_ONCE(" the  fails for key less then min this uck si ucked in bs uck "<<Mask(),)
return bad();
}
Cindex::fake_vector  li=m_ci.find(key); 
//MM_MSG(" ASSFACKL "<<MMPR(li.size()))
// TODO should handle <4 better than this.. 
if ( li.size()==0)  return bad(); 
if ( li.size()==bad())  return bad(); 
//const Locations & loc= (*li).second; // m_index[key];
MM_SZ_LOOP(jj,li,sz)
{
i=0;
const IdxTy locii=li[jj]; // (*ii); // .first;
// this COULD look at follow on keys but too slow 
i+=m_isize;
IdxTy j=0; // for (IdxTy j=0; j<m_isize; ++j)
while (true)
{ // only complete keys should be indexed 
const IdxTy k=j+i;
if ( p[j+i]== 0 ) return locii;
const char cs=m_s[locii+k];
if (cs==0 ) break; 
if ( p[k]!=cs ) {  break;  }
++j;
} // j 
} // ii 
return bad(); 
}

static char * goods() 
{
static char t[256];
static bool init=false;
if (!init)
{
::memset(&t[0],0,256);
t['A']=1;
t['G']=1;
t['C']=1;
t['T']=1;

init=true;
}

return &t[0];
}
IdxTy find_first_agct2(const char * p ) const
{
MM_ONCE(" untested lol", ) 
// TODO the first key chars must be exact 
KeyCode key=Key(p);
const char * g=goods();
Cindex::fake_vector  li=m_ci.find(key); 
if ( li.size()==0)  return bad(); 
if ( li.size()==bad())  return bad(); 
MM_SZ_LOOP(jj,li,sz)
{
const IdxTy locii=li[jj]; // (*ii); // .first;
IdxTy j=0; // for (IdxTy j=0; j<m_isize; ++j)
while (true)
{ // only complete keys should be indexed 
const IdxTy k=j+m_isize;
if ( p[j]== 0 ) return locii;
const char cs=m_s[locii+k];
// readoff is now ok 
if (cs==0 ) return locii; //if (cs==0 ) break; 
// try harder. 
if ( p[k]!=cs ) { if (g[IdxTy(cs)]!=0) if (g[IdxTy(p[k])]!=0)  break;  }
++j;
} // j 
} // ii 
return bad(); 
}

IdxTy find_best_agct2(const char * p ) const
{
MM_ONCE(" untested lol", ) 
//MM_ERR(MMPR2(p,strlen(p)))
// TODO the first key chars must be exact 
KeyCode key=Key(p);
IdxTy best=bad();

const char * g=goods();
Cindex::fake_vector  li=m_ci.find(key); 
if ( li.size()==0)  return bad(); 
if ( li.size()==bad())  return bad(); 
MM_SZ_LOOP(jj,li,sz)
{
IdxTy errors=0;
const IdxTy locii=li[jj]; // (*ii); // .first;
IdxTy j=0; // for (IdxTy j=0; j<m_isize; ++j)
while (true)
{ // only complete keys should be indexed 
const IdxTy k=j+m_isize;
if ( p[k]== 0 ) break; // return errors;
if ( p[j]== 0 ) break; // return errors;
const char cs=m_s[locii+k];
// readoff was  now ok 
if (cs==0 ) { errors+=strlen(p)-k;  break;  }
// try harder. 
if ( p[k]!=cs ) { if (g[IdxTy(cs)]!=0) if (g[IdxTy(p[k])]!=0)  ++errors;  }
++j;
} // j 
if (errors<best) best=errors; 
} // ii 
return best; // bad(); 
}




////////////////////////////////////////////



//////////////////////////////////////////
//void push_back(const StrTy & x) { m_s.push_back(x);}
//const StrTy & operator[](const IdxTy i )const  { return m_s[i]; } 
//StrTy & operator[](const IdxTy i ) { return m_s[i]; } 
//IdxTy size() const { return m_s.size(); } 


void matches( std::vector<StrTy> & hits, const IdxTy & minrl, const Myt & that,const IdxTy minex=12) const
{
if (m_compact)   Matches2( hits, minrl,that,that.m_ci,m_ci,minex) ;
else   Matches( hits, minrl, that,that.m_index, m_index) ;
}

template <class AlTy  > 
void alignment_points( std::vector<AlTy> & hits, const IdxTy & minrl, const Myt & that,const IdxTy minex=0) const
//j	void align_points( std::vector<AlTy> & hits, const IdxTy & minrl
//,	const Myt & that, const Ti & that_index, const Tx & _index, const IdxTy minex=12) const

{
if (m_compact)   AlignPts( hits, minrl,that,that.m_ci,m_ci,minex) ;
else MM_ONCE(" need to compact to align ",) 
//	 AlignPts(  hits, minrl , that,  that_index, x & _index, const IdxTy minex=12) const


}
typedef KeyCode key_code; // exposing stats means it is just easier to expose key too 
typedef  std::map<KeyCode, IdxTy> key_code_map; 
typedef std::vector<IdxTy> IdxVec;
typedef IdxVec index_vector;
typedef  std::map<KeyCode, IdxVec> key_code_index; 
key_code encode(const char * c, const IdxTy flags)
{
switch (flags)
{
case 0 : { return Key(c); } 
// TODO the strlen can be slow for long streings 
case 1 : { const IdxTy len=strlen(c); return Key(c+len-m_isize); } 

default: { MM_ERR(" bad encode flags  "<<MMPR(flags))  } 
} // flags

return 0; 
} // encode


StrTy yek(const key_code & kc ) const { return Yek(kc); } 
 void stats( key_code_map & cnts, const char * c) const
{
Stats(cnts,c);
}

private:
void Init()
{
Index();

}

template <class Ti, class Tx  > 
	void Matches( std::vector<StrTy> & hits, const IdxTy & minrl
,	const Myt & that, const Ti & that_index, const Tx & _index) const
{

MM_LOOP(ii,_index)
{
const KeyCode&  key=(*ii).first;
auto jj=that_index.find(key);
if (jj==that_index.end()) continue;
auto & thismap=(*ii).second;
auto & thatmap=(*jj).second;
//const IdxTy thissz=thismap.size();
//const IdxTy thatsz=thatmap.size();
MM_LOOP(kk,thismap)
{
MM_LOOP(ll,thatmap)
{
const IdxTy i1=(*kk); // .first;
const IdxTy i2=(*ll); // .first;
//MM_ERR(" ASSFACK "<<MMPR4(i1,i2,m_s+i1,that.m_s+i2))
IdxTy n=Neq(m_s+i1+m_isize, that.m_s+i2+m_isize);
//MM_ERR(" wtf "<<MMPR(n) << MMPR4(i1,i2,Seg(m_s+i1,12) ,Seg(that.m_s+i2,12)))
if (!false) if (strncmp(m_s+i1,that.m_s+i2,m_isize)!=0)
{
MM_ERR(" wtf "<<MMPR4(i1,i2,m_s+i1,that.m_s+i2))
}
//MM_ERR(MMPR(n))
if (n>4)
{
//MM_ERR(" wtf "<<MMPR(n) << MMPR4(i1,i2,Seg(m_s+i1,m_isize+n) ,Seg(that.m_s+i2,m_isize+n)))
const IdxTy uck=n+m_isize;
char cc[uck+1];
memcpy(cc,m_s+i1,uck);
cc[uck]=0;
hits.push_back(StrTy(cc));
}
} 
} // kk
}
} // Matches


template <class Ti, class Tx  > 
	void Matches2( std::vector<StrTy> & hits, const IdxTy & minrl
,	const Myt & that, const Ti & that_index, const Tx & _index, const IdxTy minex=12) const
{

MM_SZ_LOOP(i,_index,szkeys)
{
const KeyCode&  key=_index.key(i); // (*ii).first;
auto j=that_index.find(key);
if (j.size()==bad()) continue;
if (j.size()==0) continue;
auto  thismap=_index.find_index(i); // (*ii).second;
auto  thatmap=j; //that_index.find_index(j); // (*jj).second;
//MM_ERR(MMPR2(thismap.size(),thatmap.size()))
//const IdxTy thissz=thismap.size();
//const IdxTy thatsz=thatmap.size();
MM_SZ_LOOP(kk,thismap,szk)
{
MM_SZ_LOOP(ll,thatmap,szl)
{
const IdxTy i1=thismap[kk]; // (*kk); // .first;
const IdxTy i2=thatmap[ll]; // (*ll); // .first;
//MM_ERR(" ASSFACK "<<MMPR4(i1,i2,m_s+i1,that.m_s+i2))
const char * p1=m_s+i1;
const char * p2=that.m_s+i2;
bool start=(i1==0)||(i2==0);
if (!start) start=(*(p1-1)!=*(p2-1));
if (!start) continue;
IdxTy n=Neq(p1+m_isize, p2+m_isize);
//MM_ERR(" wtf "<<MMPR(n) << MMPR4(i1,i2,Seg(m_s+i1,12) ,Seg(that.m_s+i2,12)))
if (!false) if (strncmp(p1,p2,m_isize)!=0)
{
MM_ERR(" wtf "<<MMPR4(i1,i2,m_s+i1,that.m_s+i2))
}
//MM_ERR(MMPR(n))
if (n>minex)
{
//MM_ERR(" wtf "<<MMPR(n) << MMPR4(i1,i2,Seg(m_s+i1,m_isize+n) ,Seg(that.m_s+i2,m_isize+n)))
const IdxTy uck=n+m_isize;
char cc[uck+1];
memcpy(cc,m_s+i1,uck);
cc[uck]=0;
hits.push_back(StrTy(cc));
}
} 
} // kk
}
} // Matches

/////////////////////////////////////////////////////////
template <class Ti, class Tx,class AlTy  > 
	void AlignPts( std::vector<AlTy> & hits, const IdxTy & minrl
,	const Myt & that, const Ti & that_index, const Tx & _index, const IdxTy minex=12) const
{

MM_ERR(MMPR(m_isize)<<MMPR4(_index.size(),that_index.size(),minex,__FUNCTION__))
MM_SZ_LOOP(i,_index,szkeys)
{
const KeyCode&  key=_index.key(i); // (*ii).first;
auto j=that_index.find(key);
//MM_ERR(MMPR3(Yek(key),_index.find_index(i).size(),j.size()))
if (j.size()==bad()) continue;
if (j.size()==0) continue;
auto  thismap=_index.find_index(i); // (*ii).second;
auto  thatmap=j; //that_index.find_index(j); // (*jj).second;
MM_SZ_LOOP(kk,thismap,szk)
{
MM_SZ_LOOP(ll,thatmap,szl)
{
const IdxTy i1=thismap[kk]; // (*kk); // .first;
const IdxTy i2=thatmap[ll]; // (*ll); // .first;
const char * p1=m_s+i1;
const char * p2=that.m_s+i2;
bool start=(i1==0)||(i2==0);
if (!start) start=(*(p1-1)!=*(p2-1));
if (!start) continue;
IdxTy n=Neq(p1+m_isize, p2+m_isize);
if (!false) if (strncmp(p1,p2,m_isize)!=0)
{ MM_ERR(" wtf "<<MMPR4(i1,i2,m_s+i1,that.m_s+i2)) }
// fing equality probablty messed up all lol TODO FIXME 2018-04-11
if (n>=minex)
{
//MM_ERR(" wtf "<<MMPR(n) << MMPR4(i1,i2,Seg(m_s+i1,m_isize+n) ,Seg(that.m_s+i2,m_isize+n)))
const IdxTy uck=n+m_isize;
char cc[uck+1];
memcpy(cc,m_s+i1,uck);
cc[uck]=0;
hits.push_back(AlTy(i1,i2,StrTy(cc)));
}
} 
} // kk
}
} // Matches

/////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////





StrTy Seg(const char * c, const IdxTy len) const 
{
char cc[len+1];
IdxTy i=0;
while (i<len) {if (c[i]==0) break; cc[i]=c[i]; ++i; }
cc[i]=0;
return StrTy(cc);
}
// isize is in bytes or chars not bts?
// isize is in bytes, with ==8, 8<<3 -> 64 or zero lol.. .... 
// doh it is minus one so working as expected doh 
KeyCode Mask() const { return (KeyCode(1)<<((m_isize-1)<<3)); }
//KeyCode Mask() const { return 1; }

KeyCode Key(const char * p) const
{
KeyCode key=0;
IdxTy j=0;
const bool ucked=false; // ((m_isize&7)!=0)||(m_isize==0);
if (ucked) { MM_ERR(" FUADFASD "<<MMPR3(strlen(p),j,m_isize)) } 
// 2022 is this backwards? It probably does not matter but may be confusing at the end... 
while ((p[j]!=0)&&(j<m_isize)) { key=(key<<8) | KeyCode((unsigned char)(p[j])); ++j; } 
if (ucked) {MM_ERR(" FUADFASD "<<MMPR3(strlen(p),j,m_isize)) }
//MM_ERR(MMPR4(Yek(key),m_isize,j,p))
return key; 
}
StrTy Yek(const KeyCode & kc) const
{
const IdxTy sze=m_isize+1;
char c[sze];
c[m_isize]=0;
for(IdxTy i=0; i<m_isize; ++i ) { c[m_isize-1-i]=(kc>>(i<<3))&KeyCode(255);}
return StrTy(c);
}


void Index()
{
if (m_isize>sizeof(KeyCode))
{
// need dmel here TODO FIXME 
MM_ERR(" dnager will robinson "<<MMPR2(m_isize,sizeof(KeyCode))) 
}
MM_ERR(MMPR(m_isize))
m_index.clear();
// this is not a bit mask but a minimum key amount??????????? 
// fortunatious ly  if big enough mask becomes zero LOL 
// now set to zero
const KeyCode mask=Mask(); // (KeyCode(1)<<((m_isize-1)<<3));
const char * c=m_s;
IdxTy i=0;
while (c[i]!=0)
{
KeyCode  key=Key(c+i);
//MM_ERR(MMPR4(Yek(key),m_isize,i,(c+i)))
//if (key>mask) { m_index[key][i]=1; }  // this should be a set but will use values some how. 
// if (key>mask) { m_index[key].insert(i); }  // this should be a set but will use values some how. 
if (key>=mask) { m_index[key].push_back(i); }  // this should be a set but will use values some how. 
++i;
} // while 
//Ss ssmask; ssmask<<std::hex<<mask;
//MM_ERR(MMPR4(__FUNCTION__,m_index.size(),ssmask.str(),sizeof(KeyCode)))
} // Index 
//void Stats( std::map<char,IdxTy> & card, std::map<KeyCode, IdxTy> & cnts)
void Stats( key_code_map  & cnts, const char * c) const 
{
const KeyCode mask=Mask(); // (KeyCode(1)<<((m_isize-1)<<3));
//const char * c=m_s;
IdxTy i=0;
while (c[i]!=0)
{
KeyCode  key=Key(c+i);
if (key>=mask) { ++cnts[key]; }  
++i;
} // while 

} // Stats



IdxTy Neq(const char * c1, const char * c2) const
{
IdxTy n=0;
const char * c1i=c1;
const char * c2i=c2;
while (true)
{
if (*c1i==0 ) return n;
if (*c1i!=*c2i) return n;
++c1i; ++c2i;
++n; // ASSFACK

}

return n;
} // neq


IdxTy m_isize;
const char * m_s;
IdxTy m_len;
Sindex m_index;
Cindex m_ci;
bool m_compact;

}; // mjm_string_query_group



class mjm_string_index_collection 
{
typedef  string_seq_traits::Tr  Tr;
typedef mjm_string_index_collection Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;
typedef string_indexer Si;
typedef std::vector< Si> Fid;
typedef Fid::iterator Fitor;
typedef Fid::const_iterator Fcitor;
public:

 mjm_string_index_collection() {Init(); } 

 ~mjm_string_index_collection() {} 
void Init()
{
m_key_size=0; 
}
IdxTy key_size() const { return m_key_size; } 
//note that the source string needs to persist although normally no a problem
template <class Tx> void index(Tx & seqs,const bool comp=true, const IdxTy isz=8 )
{
m_key_size=isz;
MM_SZ_LOOP(i,seqs,sz)
{
const char * c=seqs.seq(i).c_str();
m_fid.push_back(Si(c,isz));
if (comp) m_fid.back().compact();
} //  seqs
} // index

typedef  Si::key_code key_code;
typedef  Si::key_code_map key_code_map;
typedef  Si::key_code_index key_code_index;

/// the index is a map of vectors with key being encode key and list of locations
// in source fasta. 
template<class Ts>
void index_terminal_table( key_code_index & kci, Ts & seqs, const IdxTy isz=8)
{
if (m_fid.size()==0) if (seqs.size()!=0) m_fid.push_back(Si(seqs.seq(0).c_str(),isz));
MM_SZ_LOOP(i,seqs,sz)
{
const char * c=seqs.seq(i).c_str();
const key_code key=m_fid.back().encode(c,1);
//m_fid.back().stats(m,c);
//m_fid.push_back(Si(c,isz));
kci[key].push_back(i);
}
} // index_table

template<class Ts>
void index_table( key_code_index & kci, Ts & seqs, const IdxTy isz=8)
{
if (m_fid.size()==0) if (seqs.size()!=0) m_fid.push_back(Si(seqs.seq(0).c_str(),isz));
MM_SZ_LOOP(i,seqs,sz)
{
const char * c=seqs.seq(i).c_str();
std::map<key_code,IdxTy> m;
const IdxTy len=strlen(c);
if (len<8) continue;
const IdxTy leneff=len-9;
for (IdxTy j=0; j<leneff; ++j)
{
bool ok=true;
for (IdxTy k=0; k<8; ++k) { const char cx=*(c+j+k); 
const bool cxa=(cx=='A');
const bool cxc=(cx=='C');
const bool cxg=(cx=='G');
const bool cxt=(cx=='T');
if (cxa||cxc||cxg||cxt) continue;
ok=false;
break;
}
if (!ok) continue; 
const key_code key=m_fid.back().encode(c+j,0);
++m[key];
}
//m_fid.back().stats(m,c);
//m_fid.push_back(Si(c,isz));
MM_LOOP(ii,m) { kci[(*ii).first].push_back(i); } 
}
} // index_table





template<class Ts>
void flag_left_subsets(std::vector<IdxTy> & red, Ts & seqs, const IdxTy isz=8)
{
key_code_index kci;
index_terminal_table(kci,seqs,isz);
red.clear();
red.resize(seqs.size());
// just iterate over map instead? 
MM_LOOP(ii,kci)
{
const auto & v=(*ii).second;
// some of these can be subsets of each other, maybe order in size? 
MM_LOOP(jj,v)
{

} // jj 

} // ii 

}  // flag_left_subsets


// TODO FIXME need to avoid a dummy instance 
StrTy yek( const key_code & kc ) const { return m_fid.back().yek(kc); } 

template <class Tm, class Tx> void stats(Tm & m, Tx & seqs, const IdxTy isz=8 )
{
// don't save the sequences, insetad just make a dummy instance
if (m_fid.size()==0) if (seqs.size()!=0) m_fid.push_back(Si(seqs.seq(0).c_str(),isz));
MM_SZ_LOOP(i,seqs,sz)
{
const char * c=seqs.seq(i).c_str();
m_fid.back().stats(m,c);
//m_fid.push_back(Si(c,isz));
//if (comp) m_fid.back().compact();
} //  seqs
} // index




IdxTy bad() const { return ~0U; } 

Fitor begin() { return m_fid.begin(); } 
Fitor end() { return m_fid.end(); } 
Fcitor begin()const  { return m_fid.begin(); } 
Fcitor end()const  { return m_fid.end(); } 
IdxTy size() const { return m_fid.size(); } 
const Si & operator[](const IdxTy& i )  const { return m_fid[i]; } 
Si & operator[](const IdxTy& i ) { return m_fid[i]; } 

private:

Fid m_fid;
IdxTy m_key_size;

}; //  mjm_string_index_collection 




class mjm_string_seq 
{
typedef  string_seq_traits::Tr  Tr;
typedef mjm_string_seq Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef string_seq_params Logic;
typedef mjm_logic_base VariableStore;

typedef  mjm_string_query_group TestStrings;
typedef std::map<StrTy, TestStrings>  TestStringMap;

typedef mjm_fasta::fasta_file Fasta;
typedef Fasta::Ssb Ssb;
typedef std::map<StrTy, Fasta> FastaMap;

typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_canned_methods Canned;


//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;
typedef  mjm_string_index_collection  Fidc;
typedef  mjm_string_index Msi;

typedef Msi::align_type align_type;
typedef Msi::space_type space_type;
//typedef std::vector<align_type> HiVec;
typedef Msi::HiVec HiVec;


typedef mjm_align_collection Agc;
typedef std::map<StrTy, Agc> AgcMap;



typedef mjm_pawnoff<Tr> Hand;


//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;

typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

public:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public :
mjm_string_seq():m_dmel(new Dmel()),m_fidc(0) {Init();}
mjm_string_seq(int argc,char **_args) : m_dmel(new Dmel()),m_fidc(0)
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<IdxTy(argc); ++i) args[i]=_args[i];
for (IdxTy i=argc; i<ikluge; ++i) args[i]=&dummy[0];
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
//m_tree.config("-tree",i,argc,args);
//m_flp.config("-params",i,argc,args);
//configi(m_points,"-points",i,argc,args);
//m_flp.config_set("-set-param",  i,  argc, args);
//m_tree.config_set("-set-branch",  i,  argc, args);
cmdlcmd( i, argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
}
~mjm_string_seq()
{
//clear_handlers();
delete m_dmel;
delete m_fidc;
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
dest=myatoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-source") { ++i; command_modef(args[i]); ++i; }
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { ++i; clean_up(); }
if (s=="-about") { ++i; about(); }
if (confirm) {}
} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 
void start_date(const StrTy & d) { m_flp.set("start_date",d); }
void end_date(const StrTy & d) { m_flp.set("end_date",d); }

void command_modef(const char * fn)
{ std::ifstream fin(fn); CommandInterpretter li(&fin); command_mode(li); }
void FACKBOOST() { command_mode(); } 
void command_mode() { CommandInterpretter li(&std::cin); command_mode(li); }
void command_mode(const StrTy & cmd) 
{ CommandInterpretter li; li.set(cmd,1); command_mode(li); }


typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

CmdMap m_cmd_map;
CompMap m_comp_map;

 void cli_cmd( CliTy::list_type & choices,  const char * frag)
{
//MM_ERR("cli_cmd"<<MMPR(frag))
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);

}

}
 void cli_param( CliTy::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
const StrTy cmd=CliTy::word(StrTy(_cmd),0);
auto ii=m_comp_map.find(cmd);
if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag); 
}

// this needs to ue a fifo, use stupid files for now... 
IdxTy PerformScript(Ragged & drag,const Fasta & sfasta,const StrTy & script,const IdxTy flags)
{
// make stupid temp input file ( use a fifo 
MM_SZ_LOOP(i, sfasta,sz)
{
const StrTy & sn=sfasta.name(i);
const IdxTy  len=sfasta.seq(i).length();
//MM_MSG(MMPR4(i,name,len,sn))
} //i 
// do the script

// now get the stupid output file into the ragged...

// dun... doh !

return 0; 
} // PerformScript

void _quiet() { mjm_global_flags::mm_err_enable=!true; }
void _verbose() { mjm_global_flags::mm_err_enable=true; }
void push_loudness() {m_loudness.push_back( mjm_global_flags::mm_err_enable); } 
void restore_loudness() 
{
mjm_global_flags::mm_err_enable=m_loudness.back();  
m_loudness.pop_back(); 
}

void cmd_quiet(Cip & cip , LocalVar & lv )  { _quiet(); }
void cmd_verbose(Cip & cip , LocalVar & lv )  { _verbose(); }

IdxTy StreamFasta(OsTy * os,Ssb  & ssb,const IdxTy flags)
{
const bool u_to_t=Bit(flags,0);
const bool rev=Bit(flags,1);
const bool cmpl=Bit(flags,2);
const bool modify=u_to_t||rev||cmpl;
MM_ERR(MMPR4(__FUNCTION__,u_to_t,rev,cmpl)<<MMPR(flags))
Fasta f;
push_loudness();
_quiet();
IdxTy line=0;
while (true)
{
// TODO this is dumb, the "f" has no relationship here and ssb.next() would work.. 
f.next(ssb);
if (!ssb.ok()) break; 
const StrTy & n=ssb.name;
const StrTy & s=ssb.seq;
(*os)<<">"<<n<<CRLF;
if (!modify) { (*os)<<s<<CRLF; }
else
{
const IdxTy fini=s.length();
int i=(rev)?(fini-1):0;
int inc=(rev)?-1:1;
while (rev?(i>=0):(i<fini) )
{
char c=s.c_str()[i];
if (u_to_t) if(c=='U') c='T';
if (cmpl) 
{
if (c=='G') { c='C'; }
else if (c=='C') { c='G'; }
else if (c=='T') { c='A'; }
else if (c=='A') { c='T'; }
} // cmpl
(*os)<<c;
i+=inc;
} // while
(*os)<<CRLF;
} // modify
MM_ERR(MMPR3(line,n,s))
// in theory this could now block doh.. 
++line;
} //true 
restore_loudness();
return 0; 
} // StreamFasta

IdxTy PerformScript(Ragged & drag,Ssb  & ssb,const StrTy & script,const IdxTy flags)
{
typedef mjm_wscat_bot<Tr> Piper;
Piper pp;
pp.set_bro(script);
pp.launch("",0);
Fasta f;
IdxTy line=0;
StrTy out;
// make stupid temp input file ( use a fifo 
_quiet();
while (true)
{
f.next(ssb);
if (!ssb.ok()) break; 
const StrTy & n=ssb.name;
const StrTy & s=ssb.seq;
MM_ERR(MMPR3(line,n,s))
// in theory this could now block doh.. 

//const StrTy r=pp.send_recv(s,0);
StrTy r=pp.send(s);
StrTy r2=pp.read(0);
Ragged::Line l;
l.push_back(StrTy("0"));
l.push_back(n);
//l.push_back(r);
//l.push_back(r2);
drag.add(l);
out+=r;
out+=r2;
MM_ERR(MMPR3(line,r,r2))
++line;
//MM_MSG(MMPR4(i,name,len,sn))
} //true 
pp.done_output();
out+=pp.to_eof(0);
sleep(1);
out+=pp.to_eof(0);

MM_ERR(MMPR(out))
Ragged dd;
std::stringstream iss(out);
dd.load(iss);
// originally the script generated just the line but now can output each pattern and name... 
MM_LOOP(ii,dd)
{
const IdxTy i=myatoi((*ii)[0]);
StrTy & x=drag[i-1][0];
// could maintain a int count map for speed but this is ok for now... 
const IdxTy ized=myatoi(x);
Ss ss;
ss<<(ized+1);
//drag[i-1][0]="1";
drag[i-1][0]=ss.str();

}
_verbose();
return 0; 
} // PerformScript





void cmd_modify_lines(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy action=cip.p1;
//const IdxTy flags=myatoi(cip.p1);
std::ostream & os=std::cout;
 IdxTy lines=0;
 IdxTy skip=0;
const  IdxTy interval=10000;
const IdxTy w=2;
std::ifstream is(fn);
    CommandInterpretter li(&is);
li.readline_ok(false);
    while (li.nextok())
    {
		const IdxTy sz=li.size();
		if (sz<=w){ ++skip;  continue; } 
		std::vector<StrTy>  owords=li.words();
		const IdxTy len=owords[1].length();
		const char * s=owords[1].c_str();
		char c[len+1];
		for (IdxTy i=0; i<len; ++i)
		{
			c[i]=s[len-1-i];
		}
		c[sz]=0;
		owords[1]=StrTy(c);
		os<<owords[0]; // length checked already 
		for (IdxTy i=1; i<sz; ++i ) os<<" "<<owords[i]; 
		os<<CRLF;
		++lines;
		if ((lines%interval)==0)
			{MM_STATUS(MMPR2(lines,li.word(w))<<MMPR(skip))}
//		if (li.line().size()!=0) ts.push_back(li.word(1)); 
	} // nextok()


} // cmd_modify_lines

void cmd_count_uniq(Cip & cip , LocalVar & lv ) 
{
//const StrTy name=cip.p1;
//const StrTy fn=cip.p2;
//TestStrings & ts= m_queries[name];
//std::ifstream is(fn); 
std::map<StrTy,IdxTy> m;
std::istream & is =  std::cin;
 IdxTy lines=0;
 IdxTy skip=0;
const  IdxTy interval=100000;
const IdxTy w=2;
    CommandInterpretter li(&is);
li.readline_ok(false);
//    li.set_split(1,'|');
    while (li.nextok())
    {
		const IdxTy sz=li.size();
		if (sz<=w){ ++skip;  continue; } 
		++m[li.word(w)];
		++lines;
		if ((lines%interval)==0)
			{MM_STATUS(MMPR4(m.size(),lines,li.word(w),m[li.word(w)])<<MMPR(skip))}
//		if (li.line().size()!=0) ts.push_back(li.word(1)); 
	} // nextok()
//MM_ERR(MMPR3(name,fn,ts.size()))
MM_LOOP(ii,m)
{
const StrTy & str=(*ii).first;
const IdxTy & cnt=(*ii).second;
MM_MSG(MMPR2(str,cnt))

}

}

void cmd_cmp_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy f1=cip.p1;
const StrTy f2=cip.p2;
const IdxTy flags=0;
Fasta & fasta1=m_fasta_map[f1];
Fasta & fasta2=m_fasta_map[f2];
MM_ERR("reading fasta fn into name "<<MMPR4(f1,fasta1.size(),f2,fasta2.size()))
//typedef mjm_string_index Msi;
typedef Msi::summary_type Mst;
Msi msi;
msi.base_map();
//IdxTy eq,diff,ncmp;
//typedef Summary summary_type;
//summary_type comp_base_sequences(const char * s1, const char * s2 , const IdxTy & flags)
MM_SZ_LOOP(i,fasta1,sz1)
{
const char * nm1=fasta1.name(i).c_str();
MM_SZ_LOOP(j,fasta2,sz2)
{
const char * nm2=fasta2.name(j).c_str();
Mst x= msi.comp_base_sequences(fasta1.seq(i).c_str(),fasta2.seq(j).c_str(),flags);
MM_MSG(" compare "<<MMPR4(i,j,x.eq,x.diff)<<MMPR4(x.ncmp, nm1,nm2,flags))
} // j 
} // i 
}


//void start() { initscr(); }
//void finish() { endwin(); }
//char show() { refresh(); return getch(); }
//void branch(const StrTy & base, const std::vector<StrTy> & kids)

void cmd_read_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_fasta_map[name].load(fn);
MM_ERR("reading fasta fn into name "<<MMPR3(fn,m_fasta_map[name].size(),name))
}
void cmd_add_to_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy seq=cip.p2;
Fasta & f=m_fasta_map[name];
StrTy sn=cip.wif(3);
if (sn.length()==0) 
{ Ss ss ; ss<<"seq"<<(f.size()+1); sn=ss.str(); } 
sn=StrTy(">") + sn;
f.add(sn,seq); // load(fn);
MM_ERR("adding to  "<<MMPR4(name,m_fasta_map[name].size(),sn,seq))
}


void cmd_copy_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
const StrTy key=cip.wif(3);
const IdxTy nmax=key.length();
Fasta & dfasta=m_fasta_map[dname];
const Fasta & sfasta=m_fasta_map[sname];
MM_ERR("copy fasta fn into name "<<MMPR(key)<<MMPR4(dname,dfasta.size(),sname,sfasta.size()))
MM_SZ_LOOP(i, sfasta,sz)
{
const StrTy & sn=sfasta.name(i);
for(IdxTy off=0; off<sn.length(); ++off)
{
if (strncmp(sn.c_str()+off,key.c_str(),nmax)==0)
{
MM_ERR(" include entry "<<MMPR2(dfasta.size(),sn))
const StrTy & sseq=sfasta.seq(i);
dfasta.add( StrTy(sn) ,sseq);
break;
}
}

} //i 
MM_ERR("exit copy fasta fn into name "<<MMPR(key)<<MMPR4(dname,dfasta.size(),sname,sfasta.size()))
}
//AgcMap m_aligns_map;
void cmd_list_align(Cip & cip , LocalVar & lv ) 
{
//const StrTy name=cip.p1;
//const StrTy sname=cip.p2;
//const StrTy key=cip.wif(3);
//const IdxTy nmax=key.length();
//Fasta & dfasta=m_fasta_map[dname];
//const Agc & salgn=m_aligns_map[name];
//MM_ERR("copy fasta fn into name "<<MMPR(key)<<MMPR4(dname,dfasta.size(),sname,sfasta.size()))
IdxTy i=0; 
MM_LOOP(ii, m_aligns_map)
{
const StrTy & sn=(*ii).first; // sfasta.name(i);
const IdxTy  len=(*ii).second.size();
MM_MSG(MMPR3(i,sn,len))
++i;
} //i 
}


void cmd_write_drift(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy aname=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_aligns_map[aname].write_drift(fn,flags);
MM_ERR(" cmd_write_drift "<<MMPR4(m_aligns_map[aname].size(),fn,aname,flags))
}


void cmd_write_align(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy aname=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_aligns_map[aname].write_svg(fn,flags);
MM_ERR(" cmd_write_align "<<MMPR4(m_aligns_map[aname].size(),fn,aname,flags))
}


void cmd_set_annotation(Cip & cip , LocalVar & lv ) 
{
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_aligns_map[dname].set_annotation( m_ragged_map[sname]);
MM_ERR(" cmd_set_annotation "<<MMPR4(m_aligns_map[dname].size(),m_ragged_map[sname].size(),dname,flags))
}





void cmd_cmd_align(Cip & cip , LocalVar & lv ) 
{
const StrTy aname=cip.p1;
const StrTy cmd=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
m_aligns_map[aname].cmd(cmd,flags);
MM_ERR(" cmd_cmd_align "<<MMPR4(m_aligns_map[aname].size(),cmd,aname,flags))
}





void cmd_list_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
//const StrTy sname=cip.p2;
//const StrTy key=cip.wif(3);
//const IdxTy nmax=key.length();
//Fasta & dfasta=m_fasta_map[dname];
const Fasta & sfasta=m_fasta_map[name];
//MM_ERR("copy fasta fn into name "<<MMPR(key)<<MMPR4(dname,dfasta.size(),sname,sfasta.size()))
MM_SZ_LOOP(i, sfasta,sz)
{
const StrTy & sn=sfasta.name(i);
const IdxTy  len=sfasta.seq(i).length();
MM_MSG(MMPR4(i,name,len,sn))
} //i 
//MM_ERR("exit copy fasta fn into name "<<MMPR(key)<<MMPR4(dname,dfasta.size(),sname,sfasta.size()))
}

template <class Tm>
void list_a_map(const StrTy & nm, const Tm & m, const IdxTy flags)
{
IdxTy i=0;
MM_MSG(" listing "<<MMPR(nm))
MM_LOOP(ii, m)
{
MM_MSG(MMPR3(i,(*ii).first,(*ii).second.size()))
++i;
} //ii 

} // list_a_map

void cmd_list(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
//const StrTy sname=cip.p2;
//const StrTy key=cip.wif(3);
//const IdxTy nmax=key.length();
//Fasta & dfasta=m_fasta_map[dname];

//AgcMap m_aligns_map;
list_a_map("fasta_map",m_fasta_map,0);
list_a_map("ragged_map",m_ragged_map,0);
MM_MSG(m_flp.dump())
MM_MSG(m_flp.to_string())
MM_MSG(MMPR( mjm_global_flags::mm_err_enable) <<" may be disabled at compile time however "  )
}// cmd_list


void cmd_script(Cip & cip , LocalVar & lv ) 
{
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
const StrTy script=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(3));
//const IdxTy nmax=key.length();
const Fasta & sfasta=m_fasta_map[sname];
Ragged & drag=m_ragged_map[dname];
MM_ERR(MMPR4(dname,sname,script,flags)<<MMPR2(sfasta.size(),drag.size()))
PerformScript(drag,sfasta,script,flags);

//AgcMap m_aligns_map;
}// cmd_script

void cmd_stream_script(Cip & cip , LocalVar & lv ) 
{
const StrTy dname=cip.p1;
const StrTy fname=cip.p2;
const StrTy script=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(3));
//const IdxTy nmax=key.length();
//const Fasta & sfasta=m_fasta_map[sname];

if (fname=="-")
{
Ssb ssb(&std::cin);
Ragged & drag=m_ragged_map[dname];
MM_ERR(MMPR4(dname,fname,script,flags)<<MMPR(drag.size()))
PerformScript(drag,ssb,script,flags);

}
else
{
std::ifstream is(fname);
Ssb ssb(&is);
Ragged & drag=m_ragged_map[dname];
MM_ERR(MMPR4(dname,fname,script,flags)<<MMPR3(is.good(),is.eof(),drag.size()))
PerformScript(drag,ssb,script,flags);
}



//AgcMap m_aligns_map;
}// cmd_stream_script


void cmd_stream_fasta(Cip & cip , LocalVar & lv ) 
{
if (cip.params()==0)
{
Ss ss;
ss<<cip.cmd()<<" cmd_stream_fasta dname fname flags ";

MM_ERR(ss.str())
return;
} // help 
const StrTy dname=cip.p1;
const StrTy fname=cip.p2;
//const StrTy script=cip.wif(3);
const IdxTy flags=myatoi(cip.wif(3));
//const IdxTy nmax=key.length();
//const Fasta & sfasta=m_fasta_map[sname];
if (fname!="-")
{
std::ifstream is(fname);
Ssb ssb(&is);
MM_ERR(MMPR3(dname,fname,flags)<<MMPR2(is.good(),is.eof()))

if (dname=="-")
{
StreamFasta(&std::cout,ssb,flags);
return;
}
std::ofstream fos(dname);
StreamFasta(&fos,ssb,flags);
}
else // fname 
{
Ssb ssb(&std::cin);
MM_ERR(MMPR3(dname,fname,flags))
if (dname=="-")
{
StreamFasta(&std::cout,ssb,flags);
return;
}
std::ofstream fos(dname);
StreamFasta(&fos,ssb,flags);

} // fname 
//AgcMap m_aligns_map;
}// cmd_stream_fasta







void cmd_read_interleaved_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_fasta_map[name].load_interleaved(fn.c_str());
MM_ERR("reading interleaved  fasta fn into name "<<MMPR3(fn,m_fasta_map[name].size(),name))
}
void cmd_index_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const IdxTy  n=myatoi(cip.p2);
const Fasta & fasta=m_fasta_map[name];
if (n==0)
{
MM_ERR(" Danger will robinson indexing zero length keys  "<<MMPR3(name,n,fasta.size()))
}
MM_ERR(" indexing "<<MMPR3(name,n,fasta.size()))
if ((n<4)||(n>16)) { MM_ERR( " danger will robinson need indexin size to be ok "<<MMPR(n))  } 
delete m_fidc;
m_fidc= new Fidc();
(*m_fidc).index(fasta,true,n);
m_indexed_fasta=name;
MM_ERR(" done indexing "<<MMPR4(name,n,fasta.size(),m_fidc->size()))
//MM_ERR("reading interleaved  fasta fn into name "<<MMPR3(fn,m_fasta_map[name].size(),name))
}



void cmd_write_interleaved_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_fasta_map[name].write_interleaved(fn.c_str());
MM_ERR("writing  interleaved  fasta fn  name "<<MMPR3(fn,m_fasta_map[name].size(),name))
}



void cmd_write_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const bool skip_gt=Bit(flags,0);

MM_ERR(" writting to fn "<<MMPR4(fn,m_fasta_map[name].size(),name,flags))
//m_fasta_map[name].load(fn);
std::ofstream ofs(fn);
Fasta & f=m_fasta_map[name];
MM_SZ_LOOP(i,f,sz)
{
const char * nm=f.name(i).c_str();
const char * s=f.seq(i).c_str();
if (!skip_gt) { ofs<<">"; }
ofs<<nm<<CRLF;
ofs<<s<<CRLF;

}

} // cmd_write_fasta




void cmd_add_ragged(Cip & cip , LocalVar & lv )
{
Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ;
}


void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged(cip,lv,m_ragged_map); }
void cmd_dump_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_dump_ragged(std::cout, cip,lv,m_ragged_map); }



void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m_ragged_map[name].load(fn);
MM_ERR(MMPR2(m_ragged_map[name].size(),name))
}

void cmd_modify_fasta_names(Cip & cip , LocalVar & lv ) 
{
const StrTy fname=cip.p1;
const StrTy rname=cip.p2;
Fasta & f=m_fasta_map[fname];
Ragged & r=m_ragged_map[rname];
OsTy & os=std::cout;
MM_ERR(MMPR4(r.size(),rname,f.size(),fname))
MM_SZ_LOOP(i,f,sz)
{
const char * nm=f.name(i).c_str();
const char * s=f.seq(i).c_str();
// TODO FIXME dnager will robinson the scanner returns some spurious
// values need to verify 
std::vector<StrTy>  w=r.line(i);
IdxTy iinfo=2;
StrTy ns="tooshort";
if (w.size()>1) ns=w[1];
if ( w.size()<4)
{
Ss ss; for(IdxTy ii=0; ii<w.size(); ++ii) { ss<<w[ii]<<" ";}
MM_ERR(" bad lut at "<<MMPR4(ns,i,w.size(),ss.str()))
while (w.size()<3) { w.push_back("wtf"); } 
} else
{
IdxTy ized=3;
while (ized<w.size()) 
{ if (w[ized]==StrTy("0")) { iinfo=ized-1;  break; } ++ized;}
 if (w[ized]!=StrTy("0")) 
{ 
Ss ss; for(IdxTy ii=0; ii<w.size(); ++ii) { ss<<w[ii]<<" ";}
MM_ERR(" bad lut at "<<MMPR2(ns,i)<<MMPR4(w.size(),iinfo,ized,ss.str()))
  }

}
//os<<">"<<nm<<" "<<MMPR2(w[1],w[2])<<CRLF;
os<<nm<<" "<<MMPR2(w[1],w[iinfo])<<CRLF;
os<<s<<CRLF;


} // fasta i loop 



}





void cmd_read_queries(Cip & cip , LocalVar & lv ) 
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
TestStrings & ts= m_queries[name];
std::ifstream is(fn); 
    CommandInterpretter li(&is);
//    li.set_split(1,'|');
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if (sz<2) continue;
		if (li.line().size()!=0) ts.push_back(li.word(1)); 
	} // nextok()
MM_ERR(MMPR3(name,fn,ts.size()))
}


void cmd_index_analysis(Cip & cip , LocalVar & lv ) 
{
const bool show_key_in_hex=false;
const bool total_char_cnt=!false;
const StrTy fasta=cip.p1;
Fasta & seqs=m_fasta_map[fasta];
MM_ERR(" index-analysis "<<MMPR2(fasta,seqs.size()))
if (false)
{
//typedef  mjm_string_index_collection  Fidc;
std::map<StrTy,IdxTy> m;
Fidc fidc;
fidc.index(seqs);
MM_SZ_LOOP(i,fidc,fidcsz)
{


} //i 
} // false
typedef  Fidc::key_code_map key_code_map;
key_code_map kcm;
Fidc fidc;
fidc.stats(kcm,seqs);
MM_ERR(" done making map  "<<MMPR(kcm.size()))
MM_LOOP(ii,kcm)
{
const StrTy key=fidc.yek((*ii).first);
std::map<char,IdxTy> cc;
const IdxTy len=key.length();
for (IdxTy i=0; i<len; ++i ) ++cc[key.c_str()[i]];
// TODO the ucking line number is effected by hex fick
Ss ss;
ss<<" "<<MMPR(key);
const auto & khex=(*ii).first;
const IdxTy & cnt=(*ii).second;
if ( show_key_in_hex) {ss<<MMPR(std::hex<<khex); ss<<std::dec; } 
ss<<MMPR(cnt);
const IdxTy A=cc['A'];
const IdxTy C=cc['C'];
const IdxTy G=cc['G'];
const IdxTy T=cc['T'];
const int  D=len-(A+C+G+T);
ss<<MMPR4(A,C,G,T)<<MMPR(D);
//MM_LOOP(jj,cc) { ss<<" "<<((*jj).first)<<"="<<(*jj).second; } 
MM_MSG(ss.str())
//MM_MSG(MMPR3(key,std::hex<<(*ii).first,(*ii).second))
} // ii 
if ( total_char_cnt)
{
MM_ERR(" making total char counts")
std::map<char,IdxTy> cc;
IdxTy len=0;
MM_SZ_LOOP(i,seqs,szs)
{
const char * s=seqs.seq(i).c_str();
const char * sin=s;
while (*s!=0) { ++cc[*s]; ++s; }
len+=(s-sin);
} // i  
if (len!=0)
{
const IdxTy A=cc['A'];
const IdxTy C=cc['C'];
const IdxTy G=cc['G'];
const IdxTy T=cc['T'];
const int  S=(A+C+G+T);
const int  D=len-(S);
Ss ss;
ss<<MMPR4(A,C,G,T)<<MMPR3(D,len,S);
ss<<MMPR4((1.0*A/S),(1.0*C/S),(1.0*G/S),(1.0*T/S))<<MMPR3(D,len,S);
MM_MSG(ss.str())
MM_ERR(ss.str())
//MM_LOOP(jj,cc) { ss<<" "<<((*jj).first)<<"="<<(*jj).second; } 
} // len 
} // total_char_cnt

} // index_analysis 
/*
class align_type {
public:
align_type ( const IdxTy ii, const IdxTy jj, const StrTy & ss): i(ii),j(jj),s(ss) {}
const IdxTy i,j;
const StrTy s;
};
class space_type {
public:
space_type ( const IdxTy ii, const IdxTy jj): x(ii),len(jj) {}
const IdxTy x,len;
//const StrTy s;
};
*/

void cmd_align_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy fasta=cip.p1;
const StrTy f1=cip.p2;
const StrTy f2=cip.wif(3);
Fasta & seqs=m_fasta_map[fasta];
Fasta & dest=m_fasta_map[fasta];

const IdxTy i=(f1.length()==0)?0:myatoi(f1);
const IdxTy j=(f2.length()==0)?0:myatoi(f2);
//const IdxTy j=1;
MM_ERR(" fasta align "<<MMPR3(fasta,i,j))
std::map<StrTy,IdxTy> m;
Msi msi;
Fidc fidc;
fidc.index(seqs,true,4);
HiVec hitsf,hitsr;
const IdxTy guesssz=15000;
hitsf.reserve(guesssz);
hitsr.reserve(guesssz);
fidc[i].alignment_points(hitsf,0,fidc[j]);
fidc[j].alignment_points(hitsr,0,fidc[i]);
StrTy d1,d2;
MM_ERR( " done finding markers now align "<<MMPR4(i,j,hitsf.size(),hitsr.size()))
msi.make_alignment(d1,d2,seqs.seq(i),seqs.seq(j),hitsf,hitsr);
dest.add( StrTy("aligned 1") ,d1);
dest.add( StrTy("aligned 2") ,d2);


} // cmd_align_strings


#if 0 
void make_alignment( StrTy & d1, StrTy&  d2, const StrTy s1, const StrTy & s2, const HiVec &  vf, const HiVec & vr) 
{
IdxTy sz1=s1.length();
IdxTy sz2=s2.length();
IdxTy sumzed=sz1+sz2+1;
typedef std::vector<space_type> SpVec;
SpVec s1s,s2s;
int offset=0;
int sp1=0;
int sp2=0;
MM_LOOP(ii,vf)
{
IdxTy loc1=(*ii).i;
IdxTy loc2=(*ii).j;
MM_ERR(" hit "<<MMPR3(loc1,loc2,(*ii).s))
int offnew=loc2-loc1;
int del=offnew-offset;
if (del>0)
{
if (loc1<=sp1) continue;
sp1=loc1;
s1s.push_back(space_type(sp1,del));
sp1+=del;
}
else if (del<0)
{
if (loc2<=sp2) continue;
sp2=loc2;
s2s.push_back(space_type(sp2,-del));
sp2-=del;

}
offset=offnew;
} // vf 

const char space='.';
char c1[10*sumzed], c2[10*sumzed];
const char * p1=s1.c_str();
const char * p2=s2.c_str();

IdxTy ps1=0;
IdxTy ps2=0;
IdxTy pd1=0;
IdxTy pd2=0;
IdxTy sz=0;
MM_LOOP(ii,s1s)
{
sz=(*ii).x-ps1;
::memcpy(c1+pd1,p1+ps1,sz); pd1+=sz; ps1+=sz;
sz=(*ii).len;
::memset(c1+pd1,space,sz); pd1+=sz;
}
sz=sz1-ps1;
//MM_ERR(" fick "<<MMPR4(pd1,ps1,sz,sz1))
::memcpy(c1+pd1,p1+ps1,sz); pd1+=sz; ps1+=sz;
MM_LOOP(ii,s2s)
{
sz=(*ii).x-ps2;
::memset(c2+pd2,space,sz); pd2+=sz;
sz=(*ii).len;
::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz;
}
sz=sz2-ps2;
//MM_ERR(" fick "<<MMPR4(pd2,ps2,sz,sz2))
::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz;
//while (true)
//{
//::memset(c2+pd2,space,sz); pd2+=sz;
//::memcpy(c2+pd2,p2+ps2,sz); pd2+=sz; ps2+=sz;
//break; 
//} // while 
c1[pd1]=0;
c2[pd2]=0;
d1=StrTy(c1);
d2=StrTy(c2);
//MM_ERR(" str "<<MMPR2(d1.length(),d2.length()))
}



#endif

void cmd_find_strings(Cip & cip , LocalVar & lv ) 
{
const StrTy fasta=cip.p1;
Fasta & seqs=m_fasta_map[fasta];

//typedef  mjm_string_index_collection  Fidc;
//typedef std::vector< string_indexer> Fid;
std::map<StrTy,IdxTy> m;
Fidc fidc;
fidc.index(seqs);
//Fid fid;
//MM_SZ_LOOP(i,seqs,sz)
//{
//const char * c=seqs.seq(i).c_str();
//fid.push_back(string_indexer(c,16));
//fid.push_back(string_indexer(c,8));
//} //  seqs


const IdxTy guesssz=15000;
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
//const IdxTy fidsz=10; // fidc.size();
const IdxTy fidsz= fidc.size();
for(IdxTy i=0; i<fidsz; ++i)
{
IdxTy hitss=0;
for(IdxTy j=(i+1); j<fidsz; ++j){  
std::vector<StrTy> hits;
hits.reserve(guesssz);
fidc[i].matches(hits,0,fidc[j]);
MM_LOOP(jj,hits) {++hitss;  ++m[(*jj)]; 
//MM_ERR(" at "<<MMPR3(i,m.size(),hits.size()))

}
//MM_ERR(MMPR((*jj)) )
}
//MM_ERR(" at "<<MMPR3(i,m.size(),hits.size()))
if (!false) { MM_STATUS(" at "<<MMPR3(i,m.size(),hitss)) } 

} // i 


MM_LOOP(ii,m)
{
const StrTy & s=(*ii).first;
const IdxTy & cnt=(*ii).second;
const D & h=entropy(s.c_str());
MM_MSG(MMPR4(s,cnt,s.length(),h))

} // ii 
} // cmd_find_strings

D entropy(const char * s ) const
{
D e=0;
const D fac=1.0/log(2.0);
const char * p=s;
IdxTy table[255];
for (IdxTy i=0; i<255; ++i ) table[i]=0;
IdxTy n=0; 
while (*p!=0) { ++table[IdxTy(*p)];++n;  ++p;}
if (n==0) return 0; 
for (IdxTy i=0; i<255; ++i ) 
{D q=table[i]*1.0/n; if (q!=0)  e-=log(q)*q*fac ; } 
return e;
} // entorpy


void cmd_vector_search(Cip & cip , LocalVar & lv ) 
{
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
std::map<StrTy, IdxTy> counts;
TestStrings & ts= m_queries[s];
Fasta & seqs=m_fasta_map[fasta];
MM_ERR(MMPR4(s,ts.size(),fasta,seqs.size()))
MM_SZ_LOOP(i,seqs,sz)
{
const char * c=seqs.seq(i).c_str();
string_indexer si(c);
//MM_ERR(" done indexing ")
StrTy vector="";
MM_SZ_LOOP(j,ts,szs)
{
const char * cq=ts[j].c_str();
//MM_ERR(" find "<<MMPR(cq))
const IdxTy hitloc=si.find_first(cq);
bool hit=( hitloc!=si.bad());
if (hit)
{
IdxTy test=strncmp(cq,c+hitloc,strlen(cq));
if (test!=0)
{
MM_ERR( " bad hit "<<MMPR3(hitloc,cq,c+hitloc))
}

}
//const IdxTy qlen=strlen(cq);
//IdxTy ip=0;
//IdxTy jp=0; 
//char cqp=cq[j+jp]; char cp=c[i+ip];
//while ((cqp!=0)&&(cp!=0))
//while ((cp!=0))
//{
//++ip; ++jp;
// cqp=cq[j+jp];  cp=c[i+ip];
//}
if (hit){  vector+=StrTy("1"); ++counts[StrTy(cq)];  }  else vector+=StrTy("0"); 

} // j 
MM_MSG(MMPR3(i,seqs.name(i),vector))
} // i 
MM_SZ_LOOP(j,ts,szs)
{
MM_ERR(MMPR3(j,ts[j],counts[ts[j]]))
}
//MM_LOOP(ii,counts) { MM_ERR(MMPR2((*ii).first,(*ii).second)) } 

} // cmd_vector_search

void cmd_write_svg(Cip & cip , LocalVar & lv ) 
{
m_char_mat.cmd_write_svg(cip,lv);
}

void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))

} 

}

// expand all combinations of wilcard bases can get huge quickly
void cmd_enumerate_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy dfasta=cip.p1;
const StrTy sfasta=cip.p2;
Fasta & seqs=m_fasta_map[sfasta];
Fasta & dseqs=m_fasta_map[dfasta];
mjm_string_index enw;
const IdxTy maxcard=128;
IdxTy ntoobig=0;
IdxTy nexp=0; 
MM_ERR(" enumerate non-CGAT "<<MMPR4(dfasta,dseqs.size(),sfasta,seqs.size()))
MM_SZ_LOOP(i,seqs,sz)
{
const char * n=seqs.name(i).c_str();
const char * s=seqs.seq(i).c_str();
//ofs<<">"<<nm<<CRLF;
std::vector<StrTy> enumerated;
enw.enumerate_wilds(enumerated,s,maxcard);
const bool toobig=(enumerated.size()==0);
const bool expanded=(enumerated.size()>1);
if (expanded) ++nexp;
if (toobig) ++ntoobig;
MM_SZ_LOOP(j,enumerated,ensz)
{
	Ss ss;
	ss<<n;
	if (expanded) { ss<<" :expansion "; ss<< j << " OF "<<ensz; }
	dseqs.add(ss.str(),enumerated[j]);
} // j 
MM_STATUS(MMPR4(i,dseqs.size(),nexp,ntoobig))
} // i 
MM_ERR(" exit enumerate non-CGAT "<<MMPR4(dfasta,dseqs.size(),sfasta,seqs.size()) <<MMPR3(dseqs.size(),nexp,ntoobig))
} //cmd_enumerate_fasta 

//////////////////////////////////
void cmd_check_fasta_taxon(Cip & cip , LocalVar & lv ) 
{

const StrTy fasta=cip.p1;
Fasta & seqs=m_fasta_map[fasta];
MM_ERR(" check_fasta)taxon "<<MMPR2(fasta,seqs.size()))
TaxTree & tt=m_tax_tree;
MM_SZ_LOOP(i,seqs,sz)
{
const char * n=seqs.name(i).c_str();
std::vector<StrTy> lineage;
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
const IdxTy taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
tt.lineage(lineage,taxon);
MM_MSG(" check " <<MMPR4(i,taxon,n,tt.l2s(lineage)))

} // i 

} // check_fasta_taxon



IdxTy sum_one(const StrTy & s ) const
//IdxTy sum_one(const char *  sc ) const
{
IdxTy cnt=0;
const char * c= s.c_str();
while (*c!=0) { if (*c!='0')  ++cnt; ++c; }
return cnt;
}
////////////////////////////////////////////////////////////////////////

void cmd_coverage(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy lut=cip.p1;
const StrTy fasta=cip.p2;
const bool lut_is_desc=(cmd=="descend");
const bool print_result=false;
Fasta & seqs=m_fasta_map[fasta];
TaxTree & tt=m_tax_tree;
//std::map<IdxTy, IdxTy>
// TODO FIXME this creaes spurious entry for decendent cmd 
Mii &  m=m_luts[lut];
MM_ERR(" checking fasta node ids for listed nodes "<<MMPR(lut_is_desc)<<MMPR4(fasta,seqs.size(),lut,m.size()))
MM_SZ_LOOP(i,seqs,sz)
{
	const char * n=seqs.name(i).c_str();
	IdxTy taxon=~0U;
	std::vector<IdxTy> lineage;
	std::vector<StrTy> name_words;
	cip.m_li.parse_full(name_words, n,' ' );
	const IdxTy nwsz=name_words.size();
	taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
	tt.lineage(lineage,taxon);
	if (!lut_is_desc) { MM_LOOP(ii,lineage) { ++m[*ii]; } } 
	else { MM_LOOP(ii,lineage) {

//MM_MSG(tt.node_name((*ii)))
 if (tt.node_name((*ii))==lut) {MM_MSG(n); break; }  } } 

} // i 
if (print_result) {  print_node_lut(m,tt);
/*
MM_LOOP(ii,m)
{
const IdxTy node=(*ii).first;
const IdxTy coverage =(*ii).second;
MM_MSG(MMPR4(node,coverage, tt.depth(node), tt.node_name(node)))

} // ii 
*/
} // print result

} // coverage
void cmd_print_node_lut(Cip & cip , LocalVar & lv ) 
{
TaxTree & tt=m_tax_tree;
const StrTy lut=cip.p1;
print_node_lut(m_luts[lut],tt);
}
void print_node_lut(const Mii & m, const TaxTree & tt)
{
//Mii &  m=m_luts[lut];
MM_LOOP(ii,m)
{
const IdxTy node=(*ii).first;
const IdxTy coverage =(*ii).second;
const StrTy & name=tt.node_name(node);
MM_MSG(MMPR4(node,coverage, tt.depth(node), name))

} // ii 

} // print_node_lut 

void cmd_hit_or_miss_fasta(Cip & cip , LocalVar & lv ) 
{ 
typedef std::map<StrTy, IdxTy> HaveMap;
//typedef  mjm_string_index_collection  Fidc;
//typedef std::vector< string_indexer> Fid;
const StrTy deflin=StrTy("nolineage");
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const bool print_haves=m_flp.print_haves();
const bool print_havenots=m_flp.print_havenots(); // false;
const bool print_counts=m_flp.print_counts(); // !false;
const bool print_if_have=m_flp.print_if_have(); // !false;
const bool suppress_vector=m_flp.suppress_vector(); // !false;
const bool mode=m_flp.add_level(); // !false;
const bool skip_counts=!(print_haves||print_havenots||print_counts);
const bool make_vector=skip_counts;
const bool print_hit=m_flp.print_hit();
const IdxTy accmode=skip_counts?4:m_flp.accmode();
const IdxTy maxdepth=m_flp.maxdepth();

MM_ERR(" hit_or_miss "<<MMPR4(s,fasta,accmode,print_counts)<<MMPR4(print_haves,print_havenots,print_if_have,mode)<<MMPR2(suppress_vector,print_hit))
//typedef std::vector<StrTy>  Vec;

//Vec vec;
CharMat &  vec=m_char_mat;
vec.clear();

//TestStrings & ts= m_queries[s];
Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
//if (make_vector) vec.reserve(seqs.size());
//if (make_vector) vec=Vec(seqs.size());
if (make_vector) vec.size(seqs.size(),fts.size());
std::map<StrTy,IdxTy> m;
Fidc fidc;
fidc.index(seqs,true);

MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
#if 0 
Ragged lineages;
MM_SZ_LOOP(i,seqs,sz)
{
const char * n=seqs.name(i).c_str();
std::vector<StrTy> name_words;
std::vector<StrTy> lineage;
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
const IdxTy taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
tt.lineage(lineage,taxon);
lineage.push_back(lineage);
} // i 
MM_ERR(" done finding lineages  "<<MMPR(lineages.size()))
#endif

TaxTree & tt=m_tax_tree;
MM_ERR(MMPR4(s,fts.size(),fasta,seqs.size()))
MM_ERR(MMPR(tt.size()))
MM_SZ_LOOP(j,fts,szs)
{
HaveMap have,havenot;
std::map<StrTy, IdxTy> counts;
const char * cq=fts.seq(j).c_str(); // ts[j].c_str();
MM_SZ_LOOP(i,seqs,sz)
{
//const char * c=seqs.seq(i).c_str();
const char * n=seqs.name(i).c_str();
IdxTy taxon=~0U;
// TODO this is slow, do just once  but then memory and tokenizing etc 
#if 1 
std::vector<StrTy> lineage;
if (!skip_counts) { 
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
tt.lineage(lineage,taxon);
}
#else
const std::vector<StrTy>&  lineage=lineages[i];

#endif
if (lineage.size()==0) lineage.push_back(deflin); 
//const IdxTy hitloc=fidc[i].find_first(cq);
const IdxTy hitloc=fidc[i].find_first_agct(cq);
bool hit=( hitloc!=fidc.bad());
HaveMap & mr =(hit)?have:havenot;
//IdxTy depth=0;
//TODO FIXME the lineage entries need to have whitespace removed 
switch (accmode)
{
case 0:{MM_SZ_LOOP(ii,lineage,szlin)  // all levels
		//{ ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{ucnts(mr,counts,lineage,ii,mode); } break; } //  ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
case 1:{MM_SZ_LOOP(ii,lineage,szlin)  // just top few 
		//{if ((ii+maxdepth)<szlin) continue; ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{if ((ii+maxdepth)<szlin) continue; ucnts(mr,counts,lineage,ii,mode);}break;} 
case 2:{MM_SZ_LOOP(ii,lineage,szlin)  // just phlya and equivalent s
		//{if ((ii+3)<szlin) continue; ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{if ((ii+3)<szlin) continue; ucnts(mr,counts,lineage,ii,mode);}break;} 
case 3:{MM_SZ_LOOP(ii,lineage,szlin)  // terminal name only 
		//{ ++mr[lineage[ii]];++counts[lineage[ii]];break; }break;} 
		{ucnts(mr,counts,lineage,ii,mode); ;break; }break;} 

case 4: { break; } // null mode 
default: MM_ERR(" unknown lineage counting mode "<<MMPR(accmode)) 
};
// taking more time than search? 
if (false) { MM_STATUS(MMPR4(i,have.size(),havenot.size(),hit)<<MMPR2(taxon,lineage.size())<<"  .. ") } 
if (hit){
if (print_hit) { MM_MSG(" hit "<<MMPR4(i,j,seqs.name(i),fts.name(j))) }
 } else { }

//MM_MSG(MMPR4(cq,c,n,hit))
//if (make_vector){  vec[i]+=hit?"1":"0";  } 
if (make_vector){  vec(i,j)=hit?'1':'0';  } 
} // i 
if (print_haves) { MM_LOOP(ii,have) { MM_MSG(" have "<<MMPR2(j,cq)<<MMPR2((*ii).first,(*ii).second)) } } 
if (print_havenots) { MM_LOOP(ii,havenot) { MM_MSG(" havenot "<<MMPR2(j,cq)<<MMPR2((*ii).first,(*ii).second)) } }
if (print_counts)
{
MM_LOOP(ii,counts)
{
const StrTy& name=(*ii).first; 
const D total=(*ii).second;
const IdxTy h=have[name];
const IdxTy hn=havenot[name];
const bool pr=( (print_if_have&&(h!=0) )|| !print_if_have);
if (pr) { const D del=D(h)-D(hn); MM_MSG(" counts "<<MMPR2(j,cq)<<MMPR4(name,total,h,hn)<<MMPR((del/total))) }
}
} // print_counts
if (!false) { MM_STATUS(MMPR3(j,have.size(),havenot.size())<<"  .. ") } 

} // j 
if (!suppress_vector)
{
//MM_SZ_LOOP(i,vec,sz)
for (IdxTy i=0; i<vec.rows(); ++i)
{
const char * n=seqs.name(i).c_str();
const IdxTy sumone=sum_one(vec[i]);
MM_MSG(MMPR4(i,n,vec[i],sumone))
}
} // suppressvector
} // cmd_hit_or_miss
/*
Look for exact string matches in [AGCT] ignoring others. The first parameter is the known library
of itnerstrint sequences and the second is the unknown samples. 
Output is organized differently depending on command and parrameteres.

3440  ./mjm_string_seq.out -source xxx2 -quit 
marchywka@marchywka-Latitude-E6510:~/d/jdft/jdftx/svn/cpp/mjm_libmesh$ cat xxx2
read-fasta samples /home/marchywka/d/zymo/run1/zr2097.180216.zymo/WithSoil.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta
#read-fasta samples /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/SVs_poor_annotated/SV_no_hits.fasta

#read-fasta samples ncbitax_16S.fasta 
read-fasta q mjm_commons16s.fasta 
#read-fasta q best_discriminators.fasta 
#read-fasta q best_16s_signatures.fasta
#read-fasta q sigs10.fasta 
#fasta-commons q samples 
fasta-hits q samples 
#fasta-summary q samples 



 ./mjm_string_seq.out -source xxx2 -quit  | head
mjm_string_seq.h1153  m_fasta_map[name].size()=137 name=samples
mjm_string_seq.h1153  m_fasta_map[name].size()=1133924 name=q
mjm_string_seq.h1698  fasta_hits  fasta=samples seqs.size()=137 s=q fts.size()=1133924 cmd=fasta-hits old=1 print_hits=1 print_tally=0
mjm_string_seq.h1701  done indexing  fasta=samples seqs.size()=137 fidc.size()=137
mjm_string_seq.h560 ONCE  untested lol
mjm_string_seq.h1709  Unknown sequence  j=0>seq975	426
mjm_string_seq.h1722 contains  j=0 i=38195 hitloc=153 iname=>MM_16SCOM_38196 0  cnt 5814 len 21 h 1.81903
mjm_string_seq.h1722 contains  j=0 i=38196 hitloc=153 iname=>MM_16SCOM_38197 0  cnt 280 len 22 h 1.81485
mjm_string_seq.h1722 contains  j=0 i=38197 hitloc=153 iname=>MM_16SCOM_38198 0  cnt 504 len 23 h 1.87418
mjm_string_seq.h1722 contains  j=0 i=38198 hitloc=153 iname=>MM_16SCOM_38199 0  cnt 56 len 24 h 1.89557
mjm_string_seq.h1722 contains  j=0 i=38199 hitloc=153 iname=>MM_16SCOM_38200 0  cnt 21396 len 25 h 1.90604
mjm_string_seq.h1722 contains  j=0 i=38200 hitloc=153 iname=>MM_16SCOM_38201 0  cnt 57 len 30 h 1.87391
mjm_string_seq.h1722 contains  j=0 i=55677 hitloc=142 iname=>MM_16SCOM_55678 0  cnt 14 len 22 h 1.90766
mjm_string_seq.h1722 contains  j=0 i=55678 hitloc=142 iname=>MM_16SCOM_55679 0  cnt 42 len 23 h 1.9258
mjm_string_seq.h1722 contains  j=0 i=55679 hitloc=142 iname=>MM_16SCOM_55680 0  cnt 14 len 28 h 1.90682


*/

// command line hits of sequence versus indexed fasta file
void cmd_seq_hits(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const bool list_hits=(cmd=="list-hits");
const StrTy seq=cip.p1;
const bool print_full=list_hits?true:false;
const bool genus_species=list_hits?false:true; 
std::map<StrTy, IdxTy> m;
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
const StrTy fasta=m_indexed_fasta;
Fasta & fts=m_fasta_map[fasta];
MM_ERR(" searchsing "<<MMPR4(cmd, fts.size(),fasta,seq))
Fidc&  fidc = *m_fidc;
const char * cq=seq.c_str(); 
MM_SZ_LOOP(i,fts,szs)
{
//const char * n=seqs.name(i).c_str();
const IdxTy hitloc=fidc[i].find_first(cq);
if (hitloc!=bad())
{
const char * qjname=fts.name(i).c_str(); // ts[j].c_str();
if (genus_species)
{
 add_ncbi_gs_to_map( m, qjname, cip, true );
/*
std::vector<StrTy> name_words;
const char * n = qjname; // seqs.name(all[i]).c_str();
cip.m_li.parse_full(name_words, n,' ' );
StrTy sep=""; 
StrTy key="";
MM_SZ_LOOP(k,name_words,nwsz) { if (k==0) continue; if (k>2) break; key=key+sep+name_words[k];sep=" "; } 
++m[key];
*/

}

if (print_full ) { MM_MSG(" hit "<<MMPR4(i,hitloc,seq,qjname)) } 

} // hitloc

} // i 
if (genus_species)
{
MM_LOOP(ii,m) { 
const StrTy & name=(*ii).first;
const IdxTy & n=(*ii).second;
MM_MSG(MMPR3(name,n, seq)) }
}
} // cmd_seq_hits

template <class Tm >
void add_ncbi_gs_to_map(Tm & m, const char * nm,  Cip & cip, const bool gs )
{
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, nm,' ' );
StrTy sep=""; 
StrTy key="";
const IdxTy kend=(gs?2:1);
MM_SZ_LOOP(k,name_words,nwsz) 
{ if (k==0) continue; if (k>kend) break; key=key+sep+name_words[k];sep=" "; } 
++m[key];
} // add_ncbi_gs_to_map
//////////////////////////////////////////////////////////////
template <class Tm >
void ncbi_normalize_genus_count(Tm & m, const StrTy &ncbifasta,  Cip & cip, const bool gs )
{

const Fasta & cseq=m_fasta_map[ncbifasta];
MM_ERR(" normalizing couts to  "<<MMPR2(ncbifasta, cseq.size()))
MM_SZ_LOOP(j,cseq,csz)
{
//const char * cq=cseq.seq(j).c_str(); 
const StrTy & nm=cseq.name(j);
add_ncbi_gs_to_map( m,  nm.c_str(),  cip,  gs );
} // j

} //ncbi_genus_count



/////////////////////////////////////////////////////////////
void SpecSplit(StrTy & base,IdxTy & n,IdxTy & m,const StrTy & key)
{
const IdxTy sz=key.length();
char c[sz+1];
const char * s=key.c_str();
::memcpy(c,s,sz+1);
IdxTy sp=0;
for(IdxTy i=0; i<sz; ++i)
{
if (c[i]==':') { c[i]=0; base=StrTy(c); sp=i+1; } 
if (c[i]=='/') { c[i]=0; n=myatoi(c+sp); sp=i+1; } 

}
m=myatoi(c+sp); 
if ((n==0)||(m==0)) { MM_ERR(" bad genus parse "<<MMPR4(key,base,n,m)) }


}


// analyze names created in distro, 
//>MM_16SCOM_390235 0  cnt 339 len 22 h 1.83713 Brevifollis:1/59 Deinococcus:49/59 Methylacidimicrobium:3/59 Prosthecobacter:4/59 Roseimicrobium:1/59 Terrimicrobium:1/59
/*
 ficking fields  i=0 imin=8 name_words[i]=>MM_16SCOM_575342
mjm_string_seq.h2158  ficking fields  i=1 imin=8 name_words[i]=0
mjm_string_seq.h2158  ficking fields  i=2 imin=8 name_words[i]=
mjm_string_seq.h2158  ficking fields  i=3 imin=8 name_words[i]=cnt
mjm_string_seq.h2158  ficking fields  i=4 imin=8 name_words[i]=24
mjm_string_seq.h2158  ficking fields  i=5 imin=8 name_words[i]=len
mjm_string_seq.h2158  ficking fields  i=6 imin=8 name_words[i]=30
mjm_string_seq.h2158  ficking fields  i=7 imin=8 name_words[i]=h
mjm_string_seq.h2158  ficking fields  i=8 imin=8 name_words[i]=1.91672
mjm_string_seq.h2158  ficking fields  i=9 imin=8 name_words[i]=Aquaspirillum:1/69

*/


// TODO FIXME add multiplicity corrections with mnornmalize 
template <class Tm >
void add_mjm_gs_to_map(Tm & m, const char * nm,  Cip & cip, const IdxTy flags, const D scale=1.0 )
{
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, nm,' ' );
StrTy sep=""; 
StrTy key="";
// there is a fcuking double blank in the name 
const IdxTy imin=9; // should be 8... 
const IdxTy lenpos=6; // should be 8... 
const IdxTy nwsz=name_words.size();
D lenfac=1;
if (nwsz>lenpos) lenfac=myatoi(name_words[lenpos]);
//for (IdxTy i=0; i<nwsz; ++i) { MM_ERR(" ucking fields "<<MMPR3(i,imin,name_words[i])) }
for (IdxTy i=imin; i<nwsz; ++i)
{ const StrTy & key=name_words[i]; 
// need to remove :n/m
IdxTy n,mm;
StrTy base;
SpecSplit(base,n,mm,key);
if (base=="") { MM_ERR(" bad key "<<MMPR4(base,key,n,mm))  continue; } 
switch (flags)
{
case 1 : { m[base]+=n*scale; break; } 
case 2 : { m[base]+=D(n)/D(mm)*scale; break; } 
case 3 : { m[base]+=D(n)/D(mm)*lenfac*scale; break; } 
default :  ++m[base];
} // switch

 }

} // add_mjm_gs_to_map


template <class Tm, class Tn >
void add_mjm_gs_to_map(Tm & m, const char * nm,  Tn & mnorma,  Cip & cip,  const IdxTy flags,const D scale=1.0 )
{
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, nm,' ' );
StrTy sep=""; 
StrTy key="";
// there is a fcuking double blank in the name 
const IdxTy imin=9; // should be 8... 
const IdxTy lenpos=6; // should be 8... 
const IdxTy nwsz=name_words.size();
D lenfac=1;
if (nwsz>lenpos) lenfac=myatoi(name_words[lenpos]);
std::vector<D> ni,mi;
std::vector<StrTy> gi;
D total=0;
//for (IdxTy i=0; i<nwsz; ++i) { MM_ERR(" ucking fields "<<MMPR3(i,imin,name_words[i])) }
for (IdxTy i=imin; i<nwsz; ++i)
{ const StrTy & key=name_words[i]; 
// need to remove :n/m
IdxTy n=0,mm=0;
StrTy base;
SpecSplit(base,n,mm,key);
if (base=="") { MM_ERR(" bad key "<<MMPR4(base,key,n,mm))  continue; } 
 D multi=mnorma[base];
if (multi==0 ) multi=1;
ni.push_back(n);
//mi.push_back(m);
mi.push_back(multi);
gi.push_back(base);
//total+=multi; // mnorma[base];
total+=n; // mnorma[base];
}

for (IdxTy i=0; i<gi.size(); ++i)
{
D n=ni[i];
D mm=mi[i];
StrTy & base=gi[i];

switch (flags)
{
case 1 : { m[base]+=n*scale; break; } 
case 2 : {
MM_ONCE(" using lenfac "<<MMPR4(lenfac,flags,total,mm),)

 m[base]+=(n)/total/mm*scale; break; } 
case 3 : {
MM_ONCE(" using lenfac "<<MMPR4(lenfac,flags,total,mm),)

 m[base]+=(n)/total/mm*lenfac*scale; break; } 
default :  ++m[base];
} // switch

 }

} // add_mjm_gs_to_map







IdxTy commons_count(const StrTy & nm,Cip & cip)
{

std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, nm.c_str(),' ' );
if (name_words.size()<5) return 0; 
return myatoi(name_words[4]);

}


// find genus distributions of common string selected for freq
// TODO FIXME this needs to normalize by number of genus recs in knowsn db lol 
void cmd_distro_commons(Cip & cip , LocalVar & lv ) 
{

const StrTy cmd=cip.cmd();
const StrTy cfasta=cip.p1;
const StrTy dfasta=cip.p2;
std::ofstream ofs(dfasta);
OsTy & os =ofs; // std::cout;
IdxTy tries=0;
//const bool print_full=false;
//const bool genus_species=true; 
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
const StrTy fasta=m_indexed_fasta;
Fasta & fts=m_fasta_map[fasta];
Fasta & cseq=m_fasta_map[cfasta];
MM_ERR(" searchsing "<<MMPR4(cmd, fts.size(),fasta,cfasta)<<MMPR2(cseq.size(),dfasta))
Fidc&  fidc = *m_fidc;
MM_SZ_LOOP(j,cseq,csz)
{
const char * cq=cseq.seq(j).c_str(); 
const StrTy & nm=cseq.name(j);
const IdxTy cnt=commons_count(nm,cip);
if (cnt==0) { MM_ERR(" dnager will roginson "<<MMPR(cnt)) } 
if ((cnt>500) || ( cnt<20)) continue;
//MM_ERR(" trying "<<MMPR2(cnt,nm))
if ((tries&1023)==0) MM_STATUS(MMPR3(tries,j,cnt)<<"     ")
++tries;
std::map<StrTy, IdxTy> m;
MM_SZ_LOOP(i,fts,szs)
{
//const char * n=seqs.name(i).c_str();
const IdxTy hitloc=fidc[i].find_first(cq);
if (hitloc!=bad()){
const char * qjname=fts.name(i).c_str(); // ts[j].c_str();
 add_ncbi_gs_to_map( m, qjname, cip, !true );
} // hitloc


} // ii
Ss ss;
IdxTy t=0;
MM_LOOP(ii,m) { t=t+(*ii).second; }
MM_LOOP(ii,m) { ss<<(*ii).first<<":"<<(*ii).second<<"/"<<t<<" "; } 
os<<cseq.name(j)<<" "<<ss.str()<<CRLF;
os<<cq<<CRLF;

} // j

}// cmd_distro_commons




/*
handles seveal commands for operating on fasta files

*/

void cmd_fasta_hits(Cip & cip , LocalVar & lv ) 
{

const StrTy cmd=cip.cmd();
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const StrTy ncbifasta=cip.wif(3);
const IdxTy scoring=3;
const bool old=(cmd!=StrTy("fasta-commons")); 
const bool summary=(cmd==StrTy("fasta-summary")); 
const bool print_hits=old&&!summary;
const bool tally_hits=!old;
const bool print_tally=!old;
const bool parse_gs_distro=true;
Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
std::map<StrTy,IdxTy> mnormalize;
MM_ERR(" making noramlize map")
ncbi_normalize_genus_count(mnormalize, ncbifasta, cip, false );
MM_ERR(" normalize map "<<MMPR(mnormalize.size()))
MM_ERR(" fasta_hits "<<MMPR4(fasta,seqs.size(),s,fts.size())<<MMPR4(cmd,old,print_hits,print_tally))
Fidc fidc;
fidc.index(seqs,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
//TaxTree & tt=m_tax_tree;
std::vector<IdxTy> cnt;
IdxTy tot=0;
if ( tally_hits) cnt.resize(fts.size());
MM_SZ_LOOP(j,seqs,sz)
{
std::map<StrTy,D> mdistro;
const char * qjname=seqs.name(j).c_str(); // ts[j].c_str();
if (print_hits) {  MM_MSG(" Unknown sequence "<<MMPR(j) <<qjname)  }
IdxTy jhits=0;
MM_SZ_LOOP(i,fts,szs)
{
const char * cq=fts.seq(i).c_str(); // ts[j].c_str();
//const char * n=seqs.name(i).c_str();
const IdxTy hitloc=fidc[j].find_first(cq);
if (hitloc!=bad())
{
if (tally_hits){ ++tot;  ++cnt[i]; }
//void add_mjm_gs_to_map(Tm & m, const char * nm,  Cip & cip, const IdxTy flags )
if (parse_gs_distro)
{
//add_mjm_gs_to_map(mdistro, fts.name(i).c_str(),  cip, scoring );
 add_mjm_gs_to_map(mdistro, fts.name(i).c_str(),mnormalize,  cip, scoring );
}
if (summary) { ++jhits;} 
if (print_hits) { 
const char * iname=fts.name(i).c_str(); // ts[j].c_str();
MM_MSG( "contains "<<MMPR4(j,i,hitloc,iname))
} // print_hits
} // hitloc 
else
{
if (parse_gs_distro)
{
//add_mjm_gs_to_map(mdistro, fts.name(i).c_str(),  cip, scoring,-1.0 );
// add_mjm_gs_to_map(mdistro, fts.name(i).c_str(),mnormalize,  cip, scoring,-1.0 );
}


} // not hitloc



} // i
if (parse_gs_distro)
{
D best=-1e100; StrTy bestk="none";
D pbest=-1e100; StrTy pbestk="none";
MM_LOOP(ii,mdistro)
{
const StrTy & genus=(*ii).first;
// TODO moved into the scorint system 
IdxTy normalize=1; // mnormalize[genus];
if (normalize==0) normalize=1;
const D &  score=(*ii).second/D(normalize);
MM_MSG(" genus score "<<MMPR4(j, qjname,genus,score))
if ( score>best) {pbest=best; pbestk=bestk;  best=score; bestk=genus; } 
else if ( score>pbest) {pbest=score; pbestk=genus; } 
} // mdistro
MM_MSG(" assign "<<MMPR4(j,qjname, bestk, best)<<" vs "<<MMPR2(pbestk,pbest))

} // parse_gs_distro

if (summary)
{
  MM_MSG(" sequence "<<MMPR3(j,jhits,szs) <<qjname) 
}
if (!print_hits&&!summary) { 
	if ( ( j%10)==0) 
		{ MM_STATUS(" cmd_fasta_hits "<<MMPR2(j,tot)<<"     ") }  
}
} // j
if (print_tally)
{
MM_SZ_LOOP(i,cnt,szq)
{
const char * cq=fts.seq(i).c_str(); // ts[j].c_str();
MM_MSG( MMPR4(i,cnt[i],tot,cq))
} // i 

} // print_tally 

} // fasta_hits






/*
similar to fasta_hits except interpret name of the query and determine knowns to which it belongs

*/

// TODO FIXME cache these as this is slow... 
StrTy name_2_gs(const StrTy & name)
{
Ss ss;
ss<<name;
StrTy x,g,s;
ss>> x;
ss>>g;
ss>>s;
return g+StrTy(" ")+s;
}
template <class Tm,class Tv > 
void freqtab(Tm & kfreq,Fasta & knowns, Tv * v, const bool cnt)
{
MM_SZ_LOOP(i,knowns,sz) {
 //++kfreq[name_2_gs(knowns.name(i))];
 const StrTy gs=name_2_gs(knowns.name(i));
if (v!=0) (*v)[i]=gs;
if (cnt) {  ++kfreq[gs]; }
else { kfreq[gs]+=knowns.seq(i).length(); }

 }
}

void modified_fasta( Fasta & d , const Fasta & s, const IdxTy flags)
{
const bool invert=((flags&1)!=0);
const bool rev=((flags&2)!=0);
MM_ONCE(" modifed flags "<<MMPR3(invert,rev,flags),)
MM_SZ_LOOP(i,s,sz)
{
const StrTy & seq=s.seq(i);
const StrTy & n=s.name(i);

const char * si=seq.c_str();
const IdxTy len=seq.length();
StrTy s2=seq;

char c[len+1];
for (IdxTy i=0; i<len; ++i)
{
if (rev)	 c[i]= si[len-1-i];
else 	 c[i]= si[i];
if (invert) {
	if (c[i]=='A') c[i]='T';
	else if (c[i]=='T') c[i]='A';
	else if (c[i]=='C') c[i]='G';
	else if (c[i]=='G') c[i]='C';
} // invert 
}

c[len]=0;
s2=StrTy(c);

d.add(n,s2);

}

} // modified_fasta

StrTy  random_string(const IdxTy len)
{
char ct[4];ct[0]='A'; ct[1]='T'; ct[2]='G'; ct[3]='C';
//static std::mt19937_64 mt;
static std::mt19937_64 mt(my_random_seed());
char c[len+1];
c[len]=0;
for (IdxTy j=0; j<len; ++j) 
{
c[j]=ct[mt()&3];
} // j 
return StrTy(c);
}

void cmd_make_fasta(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fastan=cip.p1;
const StrTy algol=cip.p2;
const IdxTy sz=myatoi(cip.wif(3));
const IdxTy len=myatoi(cip.wif(4));
MM_ERR(" making fasta "<<MMPR4(fastan,algol,sz,len))
Fasta & fasta=m_fasta_map[fastan];

if (algol=="random")
{
char ct[4];ct[0]='A'; ct[1]='T'; ct[2]='G'; ct[3]='C';
static std::mt19937_64 mt;
char c[len+1];
c[len]=0;
for (IdxTy i=0 ; i<sz; ++i)
{
for (IdxTy j=0; j<len; ++j) 
{
c[j]=ct[mt()&3];
} // j 

Ss ss;
ss<<"Random sequence "<<i; 
StrTy s2(c);
fasta.add(ss.str(),s2);
} // i 
MM_ERR(" made random "<<MMPR2(fastan,fasta.size()))
} //random 


}

IdxTy char_hits(HiVec & v, const bool i_vs_j)
{
std::map<IdxTy, IdxTy> hits;
MM_LOOP(ii,v) { IdxTy i=i_vs_j?(*ii).i:(*ii).j;
for (IdxTy k=0; k<(*ii).len; ++k) ++hits[i+k];
}
return hits.size();
}

void genus_list(std::vector<StrTy> & g, const StrTy & fname) //const
{
CommandInterpretter li;
Fasta & seqs=m_fasta_map[fname];
MM_SZ_LOOP(i,seqs,ssz)
{
std::vector<StrTy> w;
const StrTy & s=seqs.name(i);
li.parse_full(w, s.c_str(),' ' );
if ( w.size()>1) g.push_back(w[1]); else g.push_back(StrTy("none")); 

} // i

} // genus_map


void incseq(StrTy & seq)
{
const IdxTy sz=seq.length();
char c[sz+2];
c[sz+1]=0;
bool cin=true;
for (IdxTy i=0; i<seq.length(); ++i)
{
const char ci=seq.c_str()[i];
char cn=ci;
if (cin){
cin=false;
if (ci=='A') cn='C';
else if (ci=='C') cn='G';
else if (ci=='G') cn='T';
else if (ci=='T'){  cn='A'; cin=true; } 
} // cin

c[i+1]=cn;
}
if (cin) 
{
c[0]='A';
seq= StrTy(c);
return;
}
seq= StrTy(c+1); 
}

/////////////////////////////////////////////////////////////////////

void cmd_genus_disc(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy seq=cip.p1;
const IdxTy mflags=myatoi(cip.p2);
const bool min_size_ok=false; // ((mflags&1)!=0);
//const StrTy seqno=cip.wif(3);
//const StrTy knono=cip.wif(4);

///f.add(sn,seq); // load(fn);
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
MM_ERR(" running against  "<<MMPR4(seq,mflags,m_indexed_fasta,knowns.size()))
Fidc&  kfidc = *m_fidc;
std::vector<StrTy> gvec;
genus_list(gvec,m_indexed_fasta);
std::map<StrTy,IdxTy> gnorm;
MM_LOOP(ii,gvec) ++gnorm[*ii]; 
while (true)
{
IdxTy nhits=0;
std::map<StrTy,IdxTy> gmap;
Fasta test= Fasta();
const StrTy sn="userinput";
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);
std::map<IdxTy,IdxTy> hs;
//MM_SZ_LOOP(j,seqs,szs)
MM_SZ_LOOP(krep,knowns,szk)
{
HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
if (hitsf.size()!=0){ ++nhits;   ++gmap[gvec[krep]]; }
//++hs[hitsf.size()];
//const StrTy & n1=seqs.name(0);
//const StrTy & n2=knowns.name(krep);
//const bool i_vs_j=true;
//IdxTy lenf=char_hits(hitsf,i_vs_j);
//MM_SZ_LOOP(i,hitsf,szf)
{
//const StrTy&  s=hitsf[i].s;
//const int ii=hitsf[i].i;
//const int jj=hitsf[i].j;
//const IdxTy len=hitsf[i].len;
//MM_MSG(MMPR3(i,s,len)<<MMPR3(ii,jj,(ii-jj))<<MMPR3(szf,krep,n2))
} // hitsf
} // krep
//MM_LOOP(ii,hs) { const IdxTy hitsz=(*ii).first; const IdxTy cnt=(*ii).second; MM_ERR(MMPR2(hitsz,cnt)) } 
IdxTy tot=0;
Ss ss; //ss<<seq<<" "; 
// gnorm[nm] could be zero but that is bug
MM_LOOP(ii,gmap) 
{ 
	const StrTy & nm=(*ii).first; 
	const IdxTy cnt=(*ii).second; 
	D f=D(cnt)/D(gnorm[nm]); 
	tot+=cnt; ss<<" "<<nm<<"="<<cnt<<","<<f; 
}
MM_MSG(seq<<" "<<tot<<" "<<nhits<<" "<<gmap.size()<<" "<<ss.str())
incseq(seq);
} // seqloop 
} //genus_disc 

///////////////////////////////////////////////////////////////////////

void cmd_seg_grow(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy sseq=cip.p1;
const IdxTy mflags=myatoi(cip.p2);
const bool min_size_ok=false; // ((mflags&1)!=0);
//const StrTy seqno=cip.wif(3);
//const StrTy knono=cip.wif(4);
tree_string_iterator tsi(sseq);
///f.add(sn,seq); // load(fn);
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
MM_ERR(" running against  "<<MMPR4(sseq,mflags,m_indexed_fasta,knowns.size()))
Fidc&  kfidc = *m_fidc;
std::vector<StrTy> gvec;
genus_list(gvec,m_indexed_fasta);
std::map<StrTy,IdxTy> gnorm;
MM_LOOP(ii,gvec) ++gnorm[*ii]; 
while (tsi)
{
IdxTy nhits=0;
std::map<StrTy,IdxTy> gmap;
Fasta test= Fasta();
const StrTy sn="userinput";
const StrTy seq=(*tsi);
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);
std::map<IdxTy,IdxTy> hs;
std::vector<IdxTy> & cans=tsi.cans();
std::vector<IdxTy> nch;
//MM_SZ_LOOP(j,seqs,szs)
MM_SZ_LOOP(_krep,knowns,szk)
{

IdxTy krep=_krep;
if (cans.size()>0)
{
if (_krep>=cans.size()) break;
krep=cans[_krep];

}
HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
if (hitsf.size()!=0){ nch.push_back(krep); ++nhits;   ++gmap[gvec[krep]]; }
//++hs[hitsf.size()];
//const StrTy & n1=seqs.name(0);
//const StrTy & n2=knowns.name(krep);
//const bool i_vs_j=true;
//IdxTy lenf=char_hits(hitsf,i_vs_j);
//MM_SZ_LOOP(i,hitsf,szf)
{
//const StrTy&  s=hitsf[i].s;
//const int ii=hitsf[i].i;
//const int jj=hitsf[i].j;
//const IdxTy len=hitsf[i].len;
//MM_MSG(MMPR3(i,s,len)<<MMPR3(ii,jj,(ii-jj))<<MMPR3(szf,krep,n2))
} // hitsf
} // krep
//MM_LOOP(ii,hs) { const IdxTy hitsz=(*ii).first; const IdxTy cnt=(*ii).second; MM_ERR(MMPR2(hitsz,cnt)) } 
IdxTy tot=0;
Ss ss; //ss<<seq<<" "; 
// gnorm[nm] could be zero but that is bug
MM_LOOP(ii,gmap) 
{ 
	const StrTy & nm=(*ii).first; 
	const IdxTy cnt=(*ii).second; 
	D f=D(cnt)/D(gnorm[nm]); 
	tot+=cnt; ss<<" "<<nm<<"="<<cnt<<","<<f; 
}
if ( nch.size()!=0) 
{
MM_MSG(seq<<" len="<<seq.length()<<" tsi-size="<<tsi.size()<<" "<<tot<<" "<<nhits<<" "<<gmap.size()<<" genera "<<ss.str())
}
//incseq(seq);
//tsi.inc(tot!=0);
//tsi.inc(nhits>1);
//MM_ERR(MMPR2(tsi.cans().size(),nch.size()))
tsi.cans()=nch;
tsi.inc(gmap.size()>1);
//tsi.cans()=nch;
} // seqloop 
} //genus_disc 
/////////////////////////////////////////////////////////////////////////////////

#if 1  
void cmd_shrink_group(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy sseq=cip.p1;
const IdxTy mflags=myatoi(cip.p2);
const IdxTy cutoff=myatoi(cip.wif(3));
const IdxTy maxlen=50; // myatoi(cip.wif(2));
D _fmaxmin=.9;
IdxTy _gmszmin=1;

const D fmaxmin=_fmaxmin; // pairs?.5:.9;
const IdxTy gmszmin=_gmszmin; // pairs?2:1;

const StrTy sseqi=sseq; // =cip.p1;
const bool min_size_ok=false; // ((mflags&1)!=0);
//const bool min_size_ok=false; // ((mflags&1)!=0);
//const bool  doing_sequental_search=false; // ) incseq(sseq);
const IdxTy freq=100;
const IdxTy statoff=10;
IdxTy kci_iter=0;
MM_ERR(" shrink_group " <<MMPR4(cmd,fmaxmin,gmszmin,cutoff)<< MMPR4(maxlen,sseq,cmd,mflags))

if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
IdxTy iter=0;
Fidc&  kfidc = *m_fidc;
Fidc::key_code_index kci;
MM_ERR(" start  indexing kci")
kfidc.index_table(kci,knowns);
MM_ERR(" done indexing kci size "<< MMPR2(kci.size(),(1<<(8*2))))
MM_LOOP(kk,kci)
{

++kci_iter;
sseq=kfidc.yek((*kk).first);
{ if ( sseq<sseqi ) continue;  } // 
tree_string_iterator tsi(sseq);
while (tsi)
{
IdxTy nhits=0;
//std::map<StrTy,IdxTy> gmap;
Fasta test= Fasta();
const StrTy sn="userinput";
const StrTy seq=(*tsi);
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);
std::map<IdxTy,IdxTy> hs;
std::vector<IdxTy> & cans=tsi.cans();
std::vector<IdxTy> nch;

MM_SZ_LOOP(_krep,knowns,szk)
{

IdxTy krep=_krep;
if (cans.size()>0) { if (_krep>=cans.size()) break; krep=cans[_krep]; }
else
{
const auto & locs=(*kk).second; 
if (_krep>=locs.size()) break;
krep=locs[_krep]; 
}
HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
if (hitsf.size()!=0){ nch.push_back(krep); ++nhits;   //++gmap[gvec[krep]]; 

}
} // krep
//IdxTy tot=0;
//const bool fmaxok=true; // (fmax>fmaxmin);
const IdxTy gmsz=nch.size() ; // gmap.size();
const IdxTy seqlen=seq.length() ; // gmap.size();
// there is no ucking point continuning with a just one genus... 
std::map<StrTy,IdxTy> nmmap;
MM_LOOP(ii,nch) { ++nmmap[mangle((*ii),knowns.name(*ii))]; }
const IdxTy ngenus=nmmap.size();

//const bool success=(gmsz<20)&&(gmsz>0); //  ((gmsz!=0)&&(gmsz<=gmszmin)&&fmaxok); // (fmax>.9));
//const bool tryagain= (gmsz>19)&&(seqlen<50); // ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));

// this is the current higlhly selctive one that leaves open some genuses. 
const bool tryagain= (gmsz>cutoff)&&(seqlen<maxlen); // ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));
// the above worked well 
//const bool tryagain= (ngenus>cutoff)&&(seqlen<maxlen); // ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));
//const bool success=(gmsz<20)&&(gmsz>0); //  ((gmsz!=0)&&(gmsz<=gmszmin)&&fmaxok); // (fmax>.9));
const bool success=(!tryagain)&&(gmsz>0); //  ((gmsz!=0)&&(gmsz<=gmszmin)&&fmaxok); // (fmax>.9));
if (success) // ((gmap.size()==1)&&(fmax>.9999))
{
Ss ss;
//std::map<StrTy,IdxTy> nmmap;
//MM_LOOP(ii,nch) { ss<<" "<< mangle((*ii),knowns.name(*ii)); }
//MM_LOOP(ii,nch) { ++nmmap[mangle((*ii),knowns.name(*ii))]; }

MM_LOOP(ii,nmmap) { ss<<" "<< (*ii).first<<":"<<(*ii).second; }
//MM_MSG(seq<<" len="<<seq.length()<<" "<<tot<<" "<<nhits<<" "<<gmsz<<ss.str())
MM_MSG(seq<<" len="<<seqlen<<" "<<ngenus<<" "<<gmsz<<ss.str())
}
tsi.cans()=nch;
tsi.inc(tryagain);
if ((iter%freq)==statoff)
{
const IdxTy kcisz=kci.size();
const IdxTy tsisz=tsi.size();
const IdxTy gmpsz=0; // gmap.size();
MM_STATUS(MMPR2(kci_iter,kcisz)<<MMPR2(iter,tsisz)<<MMPR4(seq.size(),sseq,gmpsz,nhits)<<"         ")
}
++iter;
} // tsi
} // kk


} // shrink_group

class mt_shrink_group_param
{
public:
// https://www.cs.cmu.edu/afs/cs/academic/class/15492-f07/www/pthreads.html
static pthread_mutex_t*  mutex1()//  { = PTHREAD_MUTEX_INITIALIZER;
{
static pthread_mutex_t _mutex1 = PTHREAD_MUTEX_INITIALIZER;
return & _mutex1;
}
void enter_serial() { pthread_mutex_lock( mutex1() ); }
void exit_serial() { pthread_mutex_unlock( mutex1() ); }
void launch(const IdxTy nthreads, void *(*thread_function)(void* ) ,mt_shrink_group_param *msgp)
{
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, thread_function, msgp+i ); }
   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}

StrTy sseqi, sseqf;
Fidc::key_code_index*  pkci; // =(*msgp).pkci;
bool min_size_ok;
Myt * p;
IdxTy maxlen,cutoff,freq,statoff,nt,i;
}; // mt_shrink_group_param



template <class Ti, class Te>
//for  ( ; kk!=kci.end(); ioff(kk,(*msgp).nt,kci.end()))
static void ioff(Ti & kk, const IdxTy n, const Te & ke)
{
for(IdxTy i=0; i<n; ++i) 
{
if (kk==ke) return;
++kk;
}
} 


//static void mt_shrink_group(mt_shrink_group_param * msgp)
static void*  mt_shrink_group(void * _msgp)
{
 mt_shrink_group_param * msgp= (mt_shrink_group_param*)_msgp;

const StrTy sseqi=(*msgp).sseqi;
const StrTy sseqf=(*msgp).sseqf;
const bool min_size_ok=(*msgp).min_size_ok;
const IdxTy maxlen=(*msgp).maxlen;
const IdxTy cutoff=(*msgp).cutoff;
const IdxTy freq=(*msgp).freq;
const IdxTy statoff=(*msgp).statoff;
// this should be derived from msgp for consistency 
Myt * p=(*msgp).p;
Fasta & knowns=(*p).m_fasta_map[(*p).m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
Fidc&  kfidc = *(*p).m_fidc;
Fidc::key_code_index&  kci=*(*msgp).pkci;
IdxTy kci_iter=0;
IdxTy iter=0;
//MM_LOOP(kk,kci)
auto kk=kci.begin();
IdxTy offset=(*msgp).i;
while ((kk!=kci.end())&&(offset!=0)) {  ++kk; --offset; } 
//for  ( ; kk!=kci.end(); kk+=(*msgp).nt)
for  ( ; kk!=kci.end(); ioff(kk,(*msgp).nt,kci.end()))
//while  (  kk!=kci.end() )
{

++kci_iter;
const StrTy &  sseq=kfidc.yek((*kk).first);
if ( sseq<sseqi ) continue; 
if ( sseq>sseqf ) break; 
tree_string_iterator tsi(sseq);
while (tsi)
{
IdxTy nhits=0;
//std::map<StrTy,IdxTy> gmap;
Fasta test= Fasta();
const StrTy sn="userinput";
const StrTy seq=(*tsi);
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);
std::map<IdxTy,IdxTy> hs;
std::vector<IdxTy> & cans=tsi.cans();
std::vector<IdxTy> nch;

MM_SZ_LOOP(_krep,knowns,szk)
{

IdxTy krep=_krep;
if (cans.size()>0) { if (_krep>=cans.size()) break; krep=cans[_krep]; }
else
{
const auto & locs=(*kk).second; 
if (_krep>=locs.size()) break;
krep=locs[_krep]; 
}
HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
if (hitsf.size()!=0){ nch.push_back(krep); ++nhits;   //++gmap[gvec[krep]]; 

}
} // krep
const IdxTy gmsz=nch.size() ; // gmap.size();
const IdxTy seqlen=seq.length() ; // gmap.size();
// there is no ucking point continuning with a just one genus... 
std::map<StrTy,IdxTy> nmmap;
MM_LOOP(ii,nch) { ++nmmap[(*p).mangle((*ii),knowns.name(*ii))]; }
const IdxTy ngenus=nmmap.size();
const bool tryagain= (gmsz>cutoff)&&(seqlen<maxlen); // ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));
//const bool tryagain= (ngenus>cutoff)&&(seqlen<maxlen); // ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));
const bool success=(!tryagain)&&(gmsz>0); //  ((gmsz!=0)&&(gmsz<=gmszmin)&&fmaxok); // (fmax>.9));
if (success) // ((gmap.size()==1)&&(fmax>.9999))
{
Ss ss;
MM_LOOP(ii,nmmap) { ss<<" "<< (*ii).first<<":"<<(*ii).second; }
//MM_MSG(seq<<" len="<<seq.length()<<" "<<tot<<" "<<nhits<<" "<<gmsz<<ss.str())
(*msgp).enter_serial();
MM_MSG(seq<<" len="<<seqlen<<" "<<ngenus<<" "<<gmsz<<ss.str())
(*msgp).exit_serial();

}
tsi.cans()=nch;
tsi.inc(tryagain);
if ((iter%freq)==statoff)
{
const IdxTy kcisz=kci.size();
const IdxTy tsisz=tsi.size();
const IdxTy gmpsz=0; // gmap.size();
(*msgp).enter_serial();
MM_STATUS(MMPR2(kci_iter,kcisz)<<MMPR2(iter,tsisz)<<MMPR4(seq.size(),sseq,gmpsz,nhits)<<"         ")
(*msgp).exit_serial();

}
++iter;
} // tsi
} // kk
return 0;
} // mt_shrink_grop


void cmd_mt_shrink_group(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
StrTy sseq=cip.p1;
//const IdxTy mflags=myatoi(cip.p2);
const IdxTy cutoff=myatoi(cip.wif(3));
const IdxTy maxlen=50; // myatoi(cip.wif(2));
//D _fmaxmin=.9;
//IdxTy _gmszmin=1;

//const D fmaxmin=_fmaxmin; // pairs?.5:.9;
//const IdxTy gmszmin=_gmszmin; // pairs?2:1;

const StrTy sseqi=sseq; // =cip.p1;
const bool min_size_ok=false; // ((mflags&1)!=0);
//const bool min_size_ok=false; // ((mflags&1)!=0);
//const bool  doing_sequental_search=false; // ) incseq(sseq);
const IdxTy freq=100;
const IdxTy statoff=10;
//IdxTy kci_iter=0;

if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
//IdxTy iter=0;
Fidc&  kfidc = *m_fidc;
Fidc::key_code_index kci;
MM_ERR(" start  indexing kci")
kfidc.index_table(kci,knowns);
MM_ERR(" done indexing kci size "<< MMPR2(kci.size(),(1<<(8*2))))
IdxTy nt=3;
// TODO FIXME these need to be interleaved NOT segmented for memory thrashing.
// the tsi iterator needs to incrment by nt... 
const StrTy seqf="TTTTTTTT";
const IdxTy klen=seqf.length();
IdxTy si=0;
IdxTy sf=0;
agct_int( si, sseqi);
agct_int( sf, seqf);
IdxTy ni=si;
IdxTy dnf=(sf-si)/nt;
// StrTy sii= agct_int(  n, klen);

auto pfunc= Myt:: mt_shrink_group; // (mt_shrink_group_param * msgp)
mt_shrink_group_param mgsgp[nt];
for (IdxTy i=0; i<nt; ++i)
{
mt_shrink_group_param * p= mgsgp+i;
// these are indclusive 
if (false) { 
(*p).sseqi=agct_int(ni,klen);
IdxTy nf=ni+dnf;
ni=nf+1;
(*p).sseqf=agct_int(nf,klen);
if (i==(nt-1)) (*p).sseqf=seqf;
else (*p).sseqf=agct_int(nf,klen);
}
else { (*p).sseqi=sseqi; (*p).sseqf=seqf; }

(*p).nt=nt;
(*p).i=i;
(*p).pkci= &kci; // =(*msgp).pkci;
(*p).min_size_ok=min_size_ok;
 (*p).p=this;
(*p).maxlen=maxlen;
(*p).cutoff=cutoff;
(*p).freq=freq;
(*p).statoff=statoff;

}
mgsgp[0].launch(nt,pfunc,&mgsgp[0]);


} // cmd_mt_shrink_group
template <class Ti > void agct_int( Ti & n, const StrTy & s)
{
const IdxTy sz=s.length();
const char * c=s.c_str();
for (IdxTy i=0; i<sz; ++i)
{
const char ci=c[i];
if (ci=='A') { n<<=2; n|=0; }
else if (ci=='C') { n<<=2; n|=1; }
else if (ci=='G') { n<<=2; n|=2; }
else if (ci=='T') { n<<=2; n|=3; }
else MM_ERR( " agct_int bad char "<<MMPR3(i,s,n)) 
}
}
template <class Ti > StrTy agct_int( const Ti & n, const IdxTy len)
{
Ti ni=n;
char  c[len+1];;
for (IdxTy i=0; i<len; ++i)
{
const IdxTy ci=ni&3;
ni>>=2;
const IdxTy j=len-i-1;
switch (ci)
{
case 0: { c[j]='A'; break;}
case 1: { c[j]='C'; break;}
case 2: { c[j]='G'; break;}
case 3: { c[j]='T'; break;}
}
}
c[len]=0;
return StrTy(c);
}



StrTy mangle(const IdxTy & i, const StrTy & s) const
{
IdxTy sz=s.length();
char c[sz+1];
const char * si=s.c_str();
IdxTy p=0;
IdxTy pst=0;
IdxTy word=0;
while (true)
{
c[p]=si[p];
if (c[p]==0) break;
if (c[p]==' ') {++word; if (word>1) { c[p]=0; break; }else {pst=p+1; c[p]=';';}} 
++p;
}
return StrTy ( c+pst);
}
#endif


void cmd_grow_genus(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const bool pairs=(cmd=="grow-digenus");
const bool triples=(cmd=="grow-trigenus");
const bool ntuples=(cmd=="grow-ngenus");
//D _fmaxmin=pairs?.5:.9;
//IdxTy _gmszmin=pairs?2:1;

D _fmaxmin=.9;
IdxTy _gmszmin=1;
if (pairs) { _fmaxmin=.5; _gmszmin=2; }
if (triples) { _fmaxmin=.5; _gmszmin=3; }
if (ntuples) { _fmaxmin=atof(cip.wif(3).c_str()); _gmszmin=myatoi(cip.wif(4)); }

const D fmaxmin=_fmaxmin; // pairs?.5:.9;
const IdxTy gmszmin=_gmszmin; // pairs?2:1;

StrTy sseq=cip.p1;
const StrTy sseqi=sseq; // =cip.p1;
const IdxTy mflags=myatoi(cip.p2);
const bool min_size_ok=false; // ((mflags&1)!=0);
const bool  doing_sequental_search=false; // ) incseq(sseq);
const IdxTy freq=100;
const IdxTy statoff=10;
IdxTy kci_iter=0;
MM_ERR(" grow_genus " <<MMPR4(cmd,triples,fmaxmin,gmszmin)<< MMPR4(pairs,sseq,cmd,mflags))
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
IdxTy iter=0;
std::vector<StrTy> gvec;
std::map<StrTy,IdxTy> gnorm;
genus_list(gvec,m_indexed_fasta);
MM_LOOP(ii,gvec) ++gnorm[*ii]; 
Fidc&  kfidc = *m_fidc;
Fidc::key_code_index kci;
if (!doing_sequental_search)
{
MM_ERR(" start  indexing kci")
kfidc.index_table(kci,knowns);
MM_ERR(" done indexing kci size "<< MMPR2(kci.size(),(1<<(8*2))))
}
// 
if (doing_sequental_search) { MM_ERR(" forgot to fix compile time code changes lol ") } 
MM_LOOP(kk,kci)
//while (true)
{
++kci_iter;
if (!doing_sequental_search)  sseq=kfidc.yek((*kk).first);
if ( !doing_sequental_search) { if ( sseq<sseqi ) continue;  } // 
tree_string_iterator tsi(sseq);
if ( false) MM_ERR(" running against  "<<MMPR4(sseq,mflags,m_indexed_fasta,knowns.size()))
while (tsi)
{
IdxTy nhits=0;
std::map<StrTy,IdxTy> gmap;
Fasta test= Fasta();
const StrTy sn="userinput";
const StrTy seq=(*tsi);
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);
std::map<IdxTy,IdxTy> hs;
std::vector<IdxTy> & cans=tsi.cans();
std::vector<IdxTy> nch;
//MM_SZ_LOOP(j,seqs,szs)
MM_SZ_LOOP(_krep,knowns,szk)
{

IdxTy krep=_krep;
if (cans.size()>0)
{
if (_krep>=cans.size()) break;
krep=cans[_krep];
}
else
{
if (!doing_sequental_search)
{
const auto & locs=(*kk).second;
if (_krep>=locs.size()) break;
krep=locs[_krep];

}

}


HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
if (hitsf.size()!=0){ nch.push_back(krep); ++nhits;   ++gmap[gvec[krep]]; }
} // krep
//MM_LOOP(ii,hs) { const IdxTy hitsz=(*ii).first; const IdxTy cnt=(*ii).second; MM_ERR(MMPR2(hitsz,cnt)) } 
IdxTy tot=0;
D fmax=0;

Ss ss; //ss<<seq<<" "; 
// gnorm[nm] could be zero but that is bug
MM_LOOP(ii,gmap) 
{ 
	const StrTy & nm=(*ii).first; 
	const IdxTy cnt=(*ii).second; 
	D f=D(cnt)/D(gnorm[nm]); 
	if (f>fmax) fmax=f;
	tot+=cnt; ss<<" "<<nm<<"="<<cnt<<","<<f; 
}
//if ( nch.size()!=0) 
const bool fmaxok=(fmax>fmaxmin);
const IdxTy gmsz=gmap.size();
//const bool success= ((gmap.size()==1)&&fmaxok); // (fmax>.9));
const bool success= ((gmsz!=0)&&(gmsz<=gmszmin)&&fmaxok); // (fmax>.9));
//const bool tryagain= ((gmap.size()>1)&&fmaxok); // (fmax>.9));
const bool tryagain= ((gmsz>gmszmin)&&fmaxok); // (fmax>.9));
if (success) // ((gmap.size()==1)&&(fmax>.9999))
{
//MM_MSG(seq<<" len="<<seq.length()<<" tsi-size="<<tsi.size()<<" "<<tot<<" "<<nhits<<" "<<gmap.size()<<" genera "<<ss.str())
//MM_MSG(seq<<" len="<<seq.length()<<" "<<tot<<" "<<nhits<<" "<<gmap.size()<<" "<<ss.str())
MM_MSG(seq<<" len="<<seq.length()<<" "<<tot<<" "<<nhits<<" "<<gmsz<<" "<<ss.str())
}
//tsi.inc(tot!=0);
//tsi.inc(nhits>1);
//MM_ERR(MMPR2(tsi.cans().size(),nch.size()))
tsi.cans()=nch;
//tsi.inc(gmap.size()>1);
tsi.inc(tryagain);
//tsi.cans()=nch;
if ((iter%freq)==statoff)
{
const IdxTy kcisz=kci.size();
const IdxTy tsisz=tsi.size();
const IdxTy gmpsz=gmap.size();
MM_STATUS(MMPR2(kci_iter,kcisz)<<MMPR2(iter,tsisz)<<MMPR4(seq.size(),sseq,gmpsz,nhits)<<"         ")
}
++iter;
} // seqloop 
if ( doing_sequental_search) incseq(sseq);
} // true 
} //genus_disc 
/////////////////////////////////////////////////////////////////////////////////



















/////////////////////////////////////////////////////////////////////////
void cmd_seg_hits(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy seq=cip.p1;
const IdxTy mflags=myatoi(cip.p2);
const bool min_size_ok=((mflags&1)!=0);
//const StrTy seqno=cip.wif(3);
//const StrTy knono=cip.wif(4);
Fasta test= Fasta();
const StrTy sn="userinput";
test.add(sn,seq);
Fidc fidc;
fidc.index(test,true);

///f.add(sn,seq); // load(fn);
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
//if (kfasta!=m_indexed_fasta)
//{
Fasta & knowns=m_fasta_map[m_indexed_fasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
MM_ERR(" running against  "<<MMPR4(seq,mflags,m_indexed_fasta,knowns.size()))

//return; 
//}
//const StrTy fasta=m_indexed_fasta;
//Fasta & fts=m_fasta_map[fasta];
//MM_ERR(" searchsing "<<MMPR4(cmd, fts.size(),fasta,seq))
Fidc&  kfidc = *m_fidc;
std::map<IdxTy,IdxTy> hs;
//MM_SZ_LOOP(j,seqs,szs)
MM_SZ_LOOP(krep,knowns,szk)
{
HiVec hitsf;
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
++hs[hitsf.size()];
//const StrTy & n1=seqs.name(0);
const StrTy & n2=knowns.name(krep);
const bool i_vs_j=true;
IdxTy lenf=char_hits(hitsf,i_vs_j);
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
MM_MSG(MMPR3(i,s,len)<<MMPR3(ii,jj,(ii-jj))<<MMPR3(szf,krep,n2))
} // hitsf
} // krep
MM_LOOP(ii,hs) { const IdxTy hitsz=(*ii).first; const IdxTy cnt=(*ii).second; MM_ERR(MMPR2(hitsz,cnt)) } 
} // seg_hits

////////////////////////////////////////////////
// TODO no initializer, not even memset lol. 
class mt_explore_hits_param
{
public:
// 2022 add for refposp
mt_explore_hits_param() : refposp(~0) {}
// https://www.cs.cmu.edu/afs/cs/academic/class/15492-f07/www/pthreads.html
static pthread_mutex_t*  mutex1()//  { = PTHREAD_MUTEX_INITIALIZER;
{
static pthread_mutex_t _mutex1 = PTHREAD_MUTEX_INITIALIZER;
return & _mutex1;
}

void enter_serial() { pthread_mutex_lock( mutex1() ); }
void exit_serial() { pthread_mutex_unlock( mutex1() ); }

static pthread_mutex_t*  mutex2()//  { = PTHREAD_MUTEX_INITIALIZER;
{
static pthread_mutex_t _mutex2 = PTHREAD_MUTEX_INITIALIZER;
return & _mutex2;
}
void enter_serial2() { pthread_mutex_lock( mutex2() ); }
void exit_serial2() { pthread_mutex_unlock( mutex2() ); }



static IdxTy & idx()
{
static IdxTy i=0;
return i;
}
static CommandInterpretter *& source()
{
static CommandInterpretter * p=0;
return p;
}
static void next(Fasta & seqs, IdxTy & ij)
{
StrTy sn="";
StrTy seq="";
IdxTy n=0;
//MM_ONCE(" danger the j will not be consistent here doh ",)
{ pthread_mutex_lock( mutex1() ); }
IdxTy & j=idx();
IdxTy jn=j;
++j;
ij=jn;
CommandInterpretter & li= *source();
while (li.nextoks())
{
const StrTy & line=li.line();
const char * c=line.c_str();
if (*c==0) { ++n; continue; }
const bool is_name=(*c=='>');
if (is_name) 
{
if (n==0)  sn=line;
else { li.save_last_line(true); break; }
}
else 
{
if (n==0) { MM_ERR(" should have a name here "<<line) }
else seq+=line;
}
++n;
}
{ pthread_mutex_unlock( mutex1() ); }
seqs.add(sn,seq); // load(fn);
}
static IdxTy next() 
{
{ pthread_mutex_lock( mutex1() ); }
//static IdxTy j=0;
IdxTy & j=idx();
IdxTy jn=j;
++j;
{ pthread_mutex_unlock( mutex1() ); }
return jn;
}

void launch(const IdxTy nthreads, void *(*thread_function)(void* ) ,mt_explore_hits_param *msgp)
{
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, thread_function, msgp+i ); }
   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}
Myt * p;
IdxTy nt,i;

StrTy fasta;
IdxTy seqnon;
IdxTy knonon;
IdxTy mflags; // =m_flp.mflags(); // "knowns";
IdxTy maxscores; // =m_flp.mflags(); // "knowns";
bool only_print_muts; // =((mflags&4)!=0);;
Fasta *seqs; // =m_fasta_map[fasta];
Fasta* dummy;
Fasta *rev;
Fasta * knowns;// =m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
Fidc *fidc;
Fidc *rfidc;
Fidc *kfidc;
StrTy fasta_dir;
bool write_align_marks;
bool save_aligns;
StrTy align_name;
std::ostream * pos;
IdxTy first_k;
IdxTy refposp;
// StrTy fasta; // =cip.p1;
// bool print_each; //=false;
// bool group_hits;// =true;
// bool skip_pathologicals; // =true;

}; // mt_explore_hits_param

static void* mt_explore_hits(void* param ) 
{
mt_explore_hits_param * mafp=(mt_explore_hits_param*)param;
Myt * p = (*mafp).p;
const StrTy fasta=(*mafp).fasta; // cip.p1;
//const bool print_each=(*mafp).print_each; // false;
//const bool group_hits=(*mafp).group_hits; // true;
//const bool skip_pathologicals=(*mafp).skip_pathologicals; // true;
//const StrTy kfasta=(*p).m_indexed_fasta;
//Fasta & seqs=(*p).m_fasta_map[kfasta];

const IdxTy seqnon =(*mafp).seqnon ;
const IdxTy knonon =(*mafp).knonon ;
const IdxTy mflags =(*mafp).mflags ;
bool only_print_muts=(*mafp).only_print_muts ;
const bool include_control_align=false;
const Fasta &seqs =*(*mafp).seqs ;
//const Fasta& dummy =*(*mafp).dummy ;
//const Fasta & rev =*(*mafp).rev ;
const Fasta & knowns =*(*mafp).knowns ;
const Fidc &fidc =*(*mafp).fidc ;
const Fidc &rfidc =*(*mafp).rfidc ;
const Fidc &kfidc =*(*mafp).kfidc ;


// StrTy fasta; // =cip.p1;
/// this should get the next j from the param for better balance and ordering
// however... 
const IdxTy szs=seqs.size();
while (true)
//MM_SZ_LOOP(j,seqs,szs)
{
IdxTy j=(*mafp).next();
if (j==bad()) break;
if (j>=szs) break;
if (seqnon!=bad()) if (seqnon!=j) continue;
const IdxTy slen=seqs.seq(j).length();
if (slen==0) continue;
IdxTy kbest=bad();
IdxTy krbest=bad();
IdxTy lfbest=0;
IdxTy lrbest=0;

MM_SZ_LOOP(krep,knowns,szk)
{
if (knonon!=bad()) if (knonon!=krep) continue;
HiVec hitsf,hitsr;
kfidc[krep].alignment_points(hitsf,0,fidc[j]);
if (include_control_align) { kfidc[krep].alignment_points(hitsr,0,rfidc[j]); } 
const StrTy & n1=seqs.name(j);
const StrTy & n2=knowns.name(krep);
const bool i_vs_j=!true;
IdxTy lenf=(*p).char_hits(hitsf,i_vs_j);
IdxTy lenr=(*p).char_hits(hitsr,i_vs_j);
// this failes due to overlap lol 
//MM_LOOP(ii,hitsf) lenf+=(*ii).len;
//MM_LOOP(ii,hitsr) lenr+=(*ii).len;
if (lenf>lfbest) { lfbest=lenf; kbest=krep; } 
if (lenr>lrbest) { lrbest=lenr; krbest=krep; } 

if (only_print_muts) if (lenf>lenr) continue;
(*mafp).enter_serial();
MM_MSG(MMPR2(j,krep)<<MMPR4(n1,n2,hitsf.size(),hitsr.size())<<MMPR2(lenf,lenr))
(*mafp).exit_serial();
//fidc[j].alignment_points(hitsr,0,kfidc[krep]);
if (false)
{
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
(*mafp).enter_serial();
MM_MSG("fwd "<<MMPR4(mflags, i,s,len)<<MMPR3(ii,jj,(ii-jj)))
(*mafp).exit_serial();
}
MM_SZ_LOOP(i,hitsr,szr)
{
const StrTy&  s=hitsr[i].s;
const int  ii=hitsr[i].i;
const int  jj=hitsr[i].j;
const IdxTy len=hitsr[i].len;
(*mafp).enter_serial();
MM_MSG("rev "<<MMPR4(mflags,i,s,len)<<MMPR3(ii,jj,(ii-jj)))
(*mafp).exit_serial();
}
} // knowns 
} // false
const StrTy & n1=seqs.name(j);
const StrTy & n2=(kbest!=bad())?knowns.name(kbest):StrTy("");
const StrTy & n2r=(krbest!=bad())?knowns.name(krbest):StrTy("");
const D score=1.0*lfbest/slen;
const bool only_print_fwd=true;

(*mafp).enter_serial();
if (only_print_fwd) 
{MM_MSG(" bests "<<MMPR(j)<<MMPR3(lfbest,kbest,score)<<MMPR2(n1,n2))}
else
{MM_MSG(" bests "<<MMPR(j)<<MMPR4(lfbest,kbest,lrbest,krbest)<<MMPR3(n1,n2,n2r)) }

(*mafp).exit_serial();

} //j  


return 0; 
} // mt_explore_hits
//////////////////////////////////////////////////////////////////


void cmd_random_scores(Cip & cip , LocalVar & lv ) 
{
//mt_explore_hits_param * mafp=(mt_explore_hits_param*)param;
//Myt * p = (*mafp).p;
IdxTy len=1200;
//const StrTy fasta=(*mafp).fasta; // cip.p1;
IdxTy sjlen=400;
IdxTy lcnt=1;
const IdxTy l1=myatoi(cip.p1);
const IdxTy l2=myatoi(cip.p2);
const IdxTy l3=myatoi(cip.wif(3));
const IdxTy flags=myatoi(cip.wif(4));
const IdxTy jmax=myatoi(cip.wif(5));
if (l1!=0) len=l1;
if (l2!=0) sjlen=l2;
if (l3!=0) lcnt=l3;
IdxTy interv=10000/lcnt;
if (interv<10) interv=10;
const IdxTy flpint=m_flp.random_interval();
if (flpint!=0) interv=flpint;
//const bool print_each=(*mafp).print_each; // false;
//const bool group_hits=(*mafp).group_hits; // true;
//const bool skip_pathologicals=(*mafp).skip_pathologicals; // true;
//const StrTy kfasta=(*p).m_indexed_fasta;
//Fasta & seqs=(*p).m_fasta_map[kfasta];
//typedef std::map< D, std::vector < IdxTy > > ScoreMap;
typedef std::vector<IdxTy> ScoreVec;
//const IdxTy seqnon =(*mafp).seqnon ;
//const IdxTy knonon =(*mafp).knonon ;
//const IdxTy mflags =(*mafp).mflags ;
//const StrTy & fasta_dir=  ((*mafp).fasta_dir); // !true;
const bool write_histo=((flags&1)==0); // (*mafp).write_align_marks;
const bool write_align_marks=((flags&2)!=0); // (*mafp).write_align_marks;
const bool write_suspicious_marks=((flags&4)!=0); // (*mafp).write_align_marks;
const bool histo_abs_count=((flags&8)!=0); // (*mafp).write_align_marks;
const bool terminal_only=((flags&16)!=0); // (*mafp).write_align_marks;
MM_ERR("random_scores "<<MMPR4(l1,l2,l3,flags)<<MMPR4(jmax,write_histo,write_align_marks,write_suspicious_marks)<<MMPR2(histo_abs_count,terminal_only))
//const IdxTy max_score_size=10;
//bool only_print_muts=(*mafp).only_print_muts ;
//const bool include_control_align=false;
Fasta  knowns;// =*(*mafp).knowns ;
Fidc kfidc ;// =*(*mafp).kfidc ;
{Ss ss;
ss<<"Random sequence probe"; // <<i; 
for(IdxTy junk=0; junk<lcnt; ++junk){
const StrTy s2=random_string(len);
knowns.add(ss.str(),s2);
}
kfidc.index(knowns,true);
}
ScoreVec sv(histo_abs_count?(len+2):101);
IdxTy j=0; // (*mafp).next();
while (j<jmax)
{
//if (seqnon!=bad()) if (seqnon!=j) continue;
//const IdxTy jdummy=0;
Fasta seqs;
const StrTy s2=random_string(sjlen);
seqs.add(StrTy("seqj"),s2);

const IdxTy slen=seqs.seq(0).length();
if (slen==0) break; // continue;

Fidc fidc;
fidc.index(seqs,true);
D best=0;
MM_SZ_LOOP(krep,knowns,szk)
{
//if (knonon!=bad()) if (knonon!=krep) continue;
HiVec hitsf; // ,hitsr;
kfidc[krep].alignment_points(hitsf,0,fidc[0]);
//if (include_control_align) { kfidc[krep].alignment_points(hitsr,0,rfidc[jdummy]); } 
//const StrTy & n1=seqs.name(0);
//const StrTy & n2=knowns.name(krep);
const bool i_vs_j=!true;
IdxTy lenf=char_hits(hitsf,i_vs_j);
D score=histo_abs_count?(D(lenf)):(D(lenf)/sjlen); 
if (score<0) score=0; 
if (score>best) best=score;


//fidc[j].alignment_points(hitsr,0,kfidc[krep]);
if (write_align_marks)
{
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
//(*mafp).enter_serial();
MM_MSG("fwd "<<MMPR2(j,krep)<<MMPR2(s,len)<<MMPR3(ii,jj,(ii-jj)))
//(*mafp).exit_serial();
}
} // false


if ( write_suspicious_marks ) // =((flags&4)!=0); // (*mafp).write_align_marks;
{
const IdxTy slen=seqs.seq(0).length();
const IdxTy klen=knowns.seq(krep).length();
const char * sp=seqs.seq(0).c_str();
const char * kp=knowns.seq(krep).c_str();

const bool too_short=(lenf<8);
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
const bool fail_jj=(strncmp(sp+jj,s.c_str(),len)!=0);
const bool fail_ii=(strncmp(kp+ii,s.c_str(),len)!=0);

const bool bad_len=(len<8)||(len!=s.length());
// this seems backwards???
const bool oori=(ii+len)>klen;
const bool oorj=(jj+len)>slen;
const bool badal=too_short||bad_len||oori||oorj;

//(*mafp).enter_serial();
if (badal) { 
Ss ss;
ss<<MMPR4(i,szf,fail_ii,fail_jj);
ss<<MMPR4(too_short,bad_len,oori,oorj);
ss<<MMPR4(ii,jj,slen,klen);
ss<<MMPR3(s,s.length(),len);
MM_MSG("baddd "<<ss.str())
//MM_MSG("baddd "<<MMPR2(j,krep)<<MMPR2(s,len)<<MMPR3(ii,jj,(ii-jj)))
}
//(*mafp).exit_serial();
}

} // suspicious 



} // knowns 
if (histo_abs_count)
{
IdxTy loc=best;
if (loc>=sv.size()) { MM_ERR(" abs meessed pip "<<MMPR2(loc,sv.size()))}
++sv[loc];
}
else
{
IdxTy bin=IdxTy(best*100.0+0);
if (bin>100) bin=100;
//MM_ERR(MMPR2(seqs.seq(0),knowns.seq(0)))
++sv[bin];
}


if (write_histo) if (((j%interv)==(interv-1)&&!terminal_only) ||(j==(jmax-1)))
{
Ss ss; 
MM_SZ_LOOP(i,sv,svz)
{
if (histo_abs_count){ ss<<j<<" "<<(i)<<" "<<sv[i]<<CRLF; } 
else { ss<<j<<" "<<(1.0*i/100)<<" "<<sv[i]<<CRLF; }
}
//{MM_MSG(ss.str())}
std::cout<<ss.str();
}

++j;

} //j  
//return 0; 
} // mt_explore_hits


//////////////////////////////////////////////////////////////////

void make_fasta_from(Fasta & d, const Fasta & s, const IdxTy & flags)
{
//os<<seqs.name(j)<<CRLF;
//os<<seqs.seq(j)<<CRLF;
//const IdxTy slen=seqs.seq(j).length();
if (flags==0) { d.add(s.name(0),s.seq(0));  }
}

// 2022-09-04 sarcina


void make_fasta_from_pos(Fasta & d, const Fasta & s, const IdxTy & pos)
{
//os<<seqs.name(j)<<CRLF;
//os<<seqs.seq(j)<<CRLF;
//const IdxTy slen=seqs.seq(j).length();
{ d.add(s.name(pos),s.seq(pos));  }
}


static void* mt_explore_streaming_hits(void* param ) 
{
mt_explore_hits_param * mafp=(mt_explore_hits_param*)param;
Myt * p = (*mafp).p;
const StrTy fasta=(*mafp).fasta; // cip.p1;
const IdxTy mflags=(*mafp).mflags; // cip.p1;
const bool min_size_ok=((mflags&1)==0);
const IdxTy k0=(*mafp).first_k; // cip.p1;
const IdxTy refposp=(*mafp).refposp; // cip.p1;
std::ostream * ppos= (*mafp).pos; // cip.p1;
//const bool print_each=(*mafp).print_each; // false;
//const bool group_hits=(*mafp).group_hits; // true;
//const bool skip_pathologicals=(*mafp).skip_pathologicals; // true;
//const StrTy kfasta=(*p).m_indexed_fasta;
//Fasta & seqs=(*p).m_fasta_map[kfasta];
typedef std::map< D, std::vector < IdxTy > > ScoreMap;

const IdxTy seqnon =(*mafp).seqnon ;
const IdxTy knonon =(*mafp).knonon ;

(*mafp).enter_serial2();
MM_ONCE(" check output as mutex split needs to be verified lol ", ) 
MM_ONCE(" SPIROPIS MESSAGES MAY BE WRITTEN WITOUT SYNCE FOR FACKED DEBGU ", ) 

(*mafp).exit_serial2();



//const IdxTy mflags =(*mafp).mflags ;
//const bool min_size_ok=((mflags&1)!=0)
if (!min_size_ok)
{
MM_ONCE(" mflags bit 0 requires exact match should use seg_hit instead scores will be meaningless but alignmarks may be helpful. Also note polarity of bit is changed here due to default lol  ",)
}
const StrTy & fasta_dir=  ((*mafp).fasta_dir); // !true;
const bool write_align_marks=(*mafp).write_align_marks;
const bool write_bests=true;
const bool save_aligns=(*mafp).save_aligns;
const bool write_fasta= (fasta_dir.length()!=0); // ((*mafp).fasta_dir.length()!=0); // !true;
const StrTy align_name=(*mafp).align_name;
IdxTy maxscores=(*mafp).maxscores; // 7;
const IdxTy max_score_size=maxscores+20;
//bool only_print_muts=(*mafp).only_print_muts ;
const bool include_control_align=false;
//const Fasta &seqs =*(*mafp).seqs ;
//const Fasta& dummy =*(*mafp).dummy ;
//const Fasta & rev =*(*mafp).rev ;
const Fasta & knowns =*(*mafp).knowns ;
//const Fidc &fidc =*(*mafp).fidc ;
//const Fidc &rfidc =*(*mafp).rfidc ;
const Fidc &kfidc =*(*mafp).kfidc ;


// StrTy fasta; // =cip.p1;
/// this should get the next j from the param for better balance and ordering
// however... 
//const IdxTy szs=seqs.size();
while (true)
//MM_SZ_LOOP(j,seqs,szs)
{
IdxTy j=0; // (*mafp).next();
//if (j==bad()) break;
//if (j>=szs) break;
const IdxTy jdummy=0;
Fasta seqs;
(*mafp).next(seqs,j);
if (seqnon!=bad()) if (seqnon!=j) continue;
const IdxTy slen=seqs.seq(jdummy).length();
if (slen==0) break; // continue;

// get a sequence from source 
//Fasta rev;
//(*p).modified_fasta(rev ,seqs,mflags);
Fidc fidc;

// this needs the same index size as the knowns .. doh 
//fidc.index(seqs,true);
//fidc.index(seqs,true,6);
fidc.index(seqs,true,kfidc.key_size());


//Fidc rfidc;
//rfidc.index(rev,true);

//IdxTy kbest=bad();
//IdxTy krbest=bad();
//IdxTy lfbest=0;
//IdxTy lrbest=0;
ScoreMap scores;

//MM_SZ_LOOP(krep,knowns,szk)
const IdxTy szk=knowns.size();
for(IdxTy krep=k0; krep<szk; ++krep)
{
// 2022 sarcina
if (krep==refposp) continue; 
//MM_ERR(MMPR3(krep,k0,szk))
if (knonon!=bad()) if (knonon!=krep) continue;
//MM_ERR(MMPR3(krep,k0,szk))
HiVec hitsf; // ,hitsr;
// match entire sequence similar to seg_hit
if (!min_size_ok)
kfidc[krep].alignment_points(hitsf,0,fidc[jdummy],knowns.seq(krep).length()-8);
else
// original align
kfidc[krep].alignment_points(hitsf,0,fidc[jdummy]);

/*
code from seg_hits, 
const bool min_size_ok=((mflags&1)!=0);
if (!min_size_ok)
{ kfidc[krep].alignment_points(hitsf,0,fidc[0],seq.length()-8);}
else { kfidc[krep].alignment_points(hitsf,0,fidc[0]); } 
*/

MM_ERR(MMPR2(krep,hitsf.size()))

//if (include_control_align) { kfidc[krep].alignment_points(hitsr,0,rfidc[jdummy]); } 
const StrTy & n1=seqs.name(jdummy);
const StrTy & n2=knowns.name(krep);
const bool i_vs_j=!true;
IdxTy lenf=(*p).char_hits(hitsf,i_vs_j);
//IdxTy lenr=(*p).char_hits(hitsr,i_vs_j);
// this failes due to overlap lol 
//MM_LOOP(ii,hitsf) lenf+=(*ii).len;
//MM_LOOP(ii,hitsr) lenr+=(*ii).len;
//if (lenf>lfbest) { lfbest=lenf; kbest=krep; } 
//if (lenr>lrbest) { lrbest=lenr; krbest=krep; } 
scores[-D(lenf)/slen].push_back(krep); 
// fick
//if (scores.size()>10) {scores.erase(scores.end()-1); } 
// TODO WTF 
//while (scores.size()>max_score_size) 
//{ auto ie= scores.end(); --ie; scores.erase(ie); } 

//if (only_print_muts) if (lenf>lenr) continue;
if (include_control_align) { 
(*mafp).enter_serial2();
//MM_MSG(MMPR2(j,krep)<<MMPR4(n1,n2,hitsf.size(),hitsr.size())<<MMPR2(lenf,lenr))
MM_MSG(MMPR2(j,krep)<<MMPR3(n1,n2,hitsf.size())<<MMPR(lenf))
(*mafp).exit_serial2();
}
if (save_aligns) { 
(*mafp).enter_serial();
//MM_MSG(MMPR2(j,krep)<<MMPR4(n1,n2,hitsf.size(),hitsr.size())<<MMPR2(lenf,lenr))
//(*p).m_aligns_map[align_name].add(j,krep,hitsf);

const bool save_target_seq=true;
if (save_target_seq)
{
const StrTy & seq1=seqs.seq(jdummy);
const StrTy & seq2=knowns.seq(krep);
 (*p).m_aligns_map[align_name].add(j,krep,hitsf,n1,n2,seq1,seq2); 
}
else { (*p).m_aligns_map[align_name].add(j,krep,hitsf,n1,n2); }


(*mafp).exit_serial();
}




//fidc[j].alignment_points(hitsr,0,kfidc[krep]);
if (write_align_marks)
{
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
(*mafp).enter_serial2();
MM_MSG("fwd "<<MMPR2(j,krep)<<MMPR4(mflags, i,s,len)<<MMPR3(ii,jj,(ii-jj)))
(*mafp).exit_serial2();
}
/*
MM_SZ_LOOP(i,hitsr,szr)
{
const StrTy&  s=hitsr[i].s;
const int  ii=hitsr[i].i;
const int  jj=hitsr[i].j;
const IdxTy len=hitsr[i].len;
(*mafp).enter_serial();
MM_MSG("rev "<<MMPR4(mflags,i,s,len)<<MMPR3(ii,jj,(ii-jj)))
(*mafp).exit_serial();
}
*/
} // knowns 
} // false

if( write_fasta)
{
const StrTy dir=fasta_dir; // "./temp/";
Ss ss;
ss<<dir<<"comp_fasta"<<j<<".fastaa";
const StrTy fn=ss.str();
MM_ERR(" wite_fasta "<<MMPR2(scores.size(),maxscores))
(*p).write_diff_fasta(fn,jdummy,seqs,knowns,scores,maxscores);


}

if (write_bests)
{
Ss ss;
ss<<" bests "<<MMPR(j);
IdxTy nscores=0;
const StrTy & n1=seqs.name(jdummy);
ss<<" n1="<<(*p).shorten(n1,0)<< " bests| ";
MM_LOOP(ii,scores)
{
if (nscores>maxscores) break;
const D score = -(*ii).first; // /slen;
const std::vector< IdxTy> v=(*ii).second;
MM_LOOP(jj,v)
{
++nscores;
if (nscores>maxscores) break;
const StrTy & n2=knowns.name(*jj);
ss<<" "<<(*p).shorten(n2,0,score);
//const StrTy & n2=(kbest!=bad())?knowns.name(kbest):StrTy("");
//const StrTy & n2r=(krbest!=bad())?knowns.name(krbest):StrTy("");
//const D score=1.0*lfbest/slen;
//const bool only_print_fwd=true;
} // jj 

} // ii 
//ss<<CRLF;
(*mafp).enter_serial2();
MM_ERR(" bests "<<MMPR4(scores.size(),max_score_size,nscores,maxscores))
//if (only_print_fwd) 
//{MM_MSG(" bests "<<MMPR(j)<<MMPR3(lfbest,kbest,score)<<MMPR2(n1,n2))}
{MM_MSG(ss.str())}
if (ppos!=0) { (*ppos)<<ss.str()<<CRLF; }
//else
//{MM_MSG(" bests "<<MMPR(j)<<MMPR4(lfbest,kbest,lrbest,krbest)<<MMPR3(n1,n2,n2r)) }

(*mafp).exit_serial2();
} // write+bests


} //j  


return 0; 
} // mt_explore_hits










template <class Tm > 
void write_diff_fasta(const StrTy & fn,const IdxTy j, const Fasta & seqs,const Fasta& knowns
	,const Tm & scores, const IdxTy & maxscores)
{
std::ofstream os(fn);
os<<seqs.name(j)<<CRLF;
os<<seqs.seq(j)<<CRLF;
//const IdxTy slen=seqs.seq(j).length();
MM_ERR(" diff_fsta "<<MMPR2(scores.size(),maxscores))
std::map<IdxTy, D> locs;
const bool  try_to_find_noted=!true;
if (try_to_find_noted)
{
// TODO FIXME faster to index... 
CommandInterpretter li;
std::vector<StrTy> wg;
li.parse_full(wg, seqs.name(j).c_str(),' ' );
StrTy seqgenus="";
if (wg.size()>2) seqgenus=wg[2];
if (seqgenus.length()!=0)
{
MM_SZ_LOOP(i,knowns,ksz)
{
std::vector<StrTy> w;
const StrTy & s=knowns.name(i);
li.parse_full(w, s.c_str(),' ' );
if ( w.size()>1) { // g.push_back(w[1]); else g.push_back(StrTy("none"));
if (w[1]==seqgenus) locs[i]=0;
}}

}
if (locs.size()==0) { MM_ERR(" no knowns for genus "<<MMPR2(seqgenus,seqs.name(j))) } 
} // try_to_find_noted 

//IdxTy maxscores=3;
IdxTy nscores=0;
MM_LOOP(ii,scores)
{
if (nscores>maxscores) break;
MM_LOOP(jj,(*ii).second)
{
++nscores;
if (nscores>maxscores) break;
locs[(*jj)]=(*ii).first;
}
} // ii
// should get the known one too.. 

MM_LOOP(ii,locs)
{
//os<<knowns.name((*ii).first)<<" score "<<((-(*ii).second)/slen)<<CRLF;
os<<knowns.name((*ii).first)<<" score "<<((-(*ii).second))<<CRLF;
os<<knowns.seq((*ii).first)<<CRLF;

}

}
//////////////////////////////////////////////////////////////////
StrTy shorten(const StrTy & s, const IdxTy flags, const D & score)
{
Ss ss;
ss<<score<<",";
ss<< shorten(s,flags); 
return ss.str();
}
StrTy shorten(const StrTy & s, const IdxTy flags)
{
const IdxTy sz=s.length();
const char * c=s.c_str();
char d[sz+1];
IdxTy i=0;
IdxTy w=0;
IdxTy p=0;
for (; i<sz; ++i)
{
if (c[i]!=' ') { d[p]=c[i]; ++p; } 
else
{
if (w==2) break;
++w;
d[p]=':';
++p;

}
} 
d[p]=0;
return StrTy(d);
} // shortern


/////////////////////////////////////////////////////////////////

void cmd_mt_explore_hits(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy kfasta=cip.p1;
const StrTy fasta=cip.p2;
const StrTy seqno=cip.wif(3);
const StrTy knono=cip.wif(4);
const IdxTy nt=3; // number of threads 
IdxTy seqnon=bad();
IdxTy knonon=bad();
IdxTy maxscores=7; // =bad();
if (seqno.length()!=0) seqnon=myatoi(seqno);
if (knono.length()!=0) knonon=myatoi(knono);
std::ostream * hit_ostream=0;
//const StrTy kfasta=m_flp.knowns_fasta(); // "knowns";
const IdxTy mflags=m_flp.mflags(); // "knowns";
//const bool use_knowns=true;
const bool only_print_muts=((mflags&4)!=0);;
MM_ERR(" enter cmd_explore_hits "<<MMPR3(cmd,fasta,kfasta)<<MMPR3(mflags,int(seqnon),int(knonon)))
//Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
//Fasta dummy;
Fasta rev;
Fasta & knowns=m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
modified_fasta(rev ,seqs,mflags);
Fidc fidc;
fidc.index(seqs,true);
Fidc rfidc;
rfidc.index(rev,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
if (kfasta!=m_indexed_fasta)
{
MM_ERR(" wrong fasta infdex want "<<kfasta<<" have "<<m_indexed_fasta)
return; 
}
Fidc&  kfidc = *m_fidc;
MM_ERR(" previously  indexed knowns  "<<MMPR3(kfasta,knowns.size(),kfidc.size()))
auto pfunc= Myt:: mt_explore_hits; // (mt_shrink_group_param * msgp)
mt_explore_hits_param mehp[nt];
for (IdxTy i=0; i<nt; ++i)
{
mt_explore_hits_param * p= mehp+i;
(*p).nt=nt; (*p).i=i; (*p).p=this; 
(*p).fasta=fasta; // cip.p1;
//(*p).print_each=print_each; // false;
//(*p).group_hits=group_hits; // true;
//(*p).skip_pathologicals=skip_pathologicals; // true;
(*p).seqnon=seqnon;
(*p).knonon=knonon;
(*p).maxscores=maxscores; // knonon;
(*p).mflags=mflags; // =m_flp.mflags(); // "knowns";
(*p).only_print_muts=only_print_muts; // =((mflags&4)!=0);;
(*p).seqs=&seqs; // =m_fasta_map[fasta];
//(*p).dummy=&dummy;
(*p).rev=&rev;
(*p).knowns=&knowns;// =m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
(*p).fidc=&fidc;
(*p).rfidc=&rfidc;
(*p).kfidc=&kfidc;
(*p).pos=hit_ostream;
// StrTy fasta; // =cip.p1;
// bool print_each; //=false;
}
mehp[0].idx()=0;
mehp[0].launch(nt,pfunc,&mehp[0]);
if (hit_ostream!=0) delete hit_ostream;

} // cmd_mt_explore_hits
////////////////////////////////////////////

void cmd_mt_explore_streaming_hits(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy kfasta=cip.p1;
const StrTy fasta=cip.p2;
const StrTy seqno=cip.wif(3);
const StrTy knono=cip.wif(4);
const StrTy fasta_dir=m_flp.fasta_dir();
const StrTy hit_file=m_flp.hit_file();
const bool write_align_marks=m_flp.write_align_marks();
const bool save_aligns=m_flp.save_aligns();
const StrTy align_name=m_flp.align_name();
const IdxTy nt=m_flp.n_threads(); // 3;
const IdxTy refpos=m_flp.refpos(); //2022-09-04 for sarcina strains  
std::ostream * ppos=(hit_file.length()==0)?0:(new std::ofstream(hit_file));
IdxTy seqnon=bad();
IdxTy knonon=bad();
IdxTy maxscores=m_flp.maxscores();
if (seqno.length()!=0) seqnon=myatoi(seqno);
if (knono.length()!=0) knonon=myatoi(knono);

//const StrTy kfasta=m_flp.knowns_fasta(); // "knowns";
const IdxTy mflags=m_flp.mflags(); // "knowns";
//const bool use_knowns=true;
const bool only_print_muts=((mflags&4)!=0);
const bool make_unknowns_from_first_known=((mflags&8)!=0);
MM_ERR(" enter cmd_explore_hits "<<MMPR4(cmd,fasta,kfasta,make_unknowns_from_first_known)<<MMPR(maxscores)<<MMPR4(mflags,int(seqnon),int(knonon),fasta_dir))
//Fasta & fts=m_fasta_map[s];
//Fasta & seqs=m_fasta_map[fasta];
//Fasta dummy;
//Fasta rev;
Fasta & knowns=m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];

// note this does not CLEAR the dest fasta 
if (make_unknowns_from_first_known)
{
	Fasta f;
	if (knowns.size()==0)
	{
		MM_ERR(" danger will robinson knowns size is zero "<<MMPR(kfasta))
	}
	//make_fasta_from(m_fasta_map[fasta] ,knowns, 0);
	//make_fasta_from(f ,knowns, 0);
	make_fasta_from_pos(f ,knowns, refpos);
	//const StrTy & cuck=((const Fasta)f).name(0U);
	const StrTy & cuck=f.name(0U);
	//MM_ERR(" writing fasta "<<MMPR3(fasta,f.size(),StrTy(f.name(0U))))
	MM_ERR(" writing fasta "<<MMPR3(fasta,f.size(),cuck))
//MM_ERR(f.seq(0U))
// FACK this looks like a compiler oerload fick 
	//f.write(fasta);
	f.write(fasta.c_str());
	MM_ERR(" wrote "<<MMPR(fasta))
}
std::ifstream is(fasta);
CommandInterpretter li(&is);
li.readline_ok(false);



//modified_fasta(rev ,seqs,mflags);
//Fidc fidc;
//fidc.index(seqs,true);
//Fidc rfidc;
//rfidc.index(rev,true);
//MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
if (kfasta!=m_indexed_fasta)
{
MM_ERR(" wrong fasta infdex want "<<kfasta<<" have "<<m_indexed_fasta)
return; 
}
Fidc&  kfidc = *m_fidc;
MM_ERR(" previously  indexed knowns  "<<MMPR3(kfasta,knowns.size(),kfidc.size()))
auto pfunc= Myt:: mt_explore_streaming_hits; // (mt_shrink_group_param * msgp)
// TODO 2022-09 sarcina allowed picking a sequence other than the first
// but now it gets included and the first is skipped wtf 
const IdxTy k0= ((refpos==0)&&make_unknowns_from_first_known) ?1:0;
mt_explore_hits_param mehp[nt];

for (IdxTy i=0; i<nt; ++i)
{
mt_explore_hits_param * p= mehp+i;
(*p).nt=nt; (*p).i=i; (*p).p=this; 
(*p).fasta=fasta; // cip.p1;
//(*p).print_each=print_each; // false;
//(*p).group_hits=group_hits; // true;
//(*p).skip_pathologicals=skip_pathologicals; // true;
(*p).seqnon=seqnon;
(*p).knonon=knonon;
(*p).mflags=mflags; // =m_flp.mflags(); // "knowns";
(*p).maxscores=maxscores; // =m_flp.mflags(); // "knowns";
(*p).only_print_muts=only_print_muts; // =((mflags&4)!=0);;
//(*p).seqs=&seqs; // =m_fasta_map[fasta];
//(*p).dummy=&dummy;
//(*p).rev=&rev;
(*p).knowns=&knowns;// =m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
//(*p).fidc=&fidc;
//(*p).rfidc=&rfidc;
(*p).kfidc=&kfidc;
(*p).fasta_dir=fasta_dir;
(*p).write_align_marks=write_align_marks;
(*p).save_aligns=save_aligns;
(*p).align_name=align_name;
(*p).pos=ppos;
(*p).first_k=k0;
// only do this if it is one of the knowns otherwise it is inited to bad valuew
if (make_unknowns_from_first_known) (*p).refposp=refpos;
// StrTy fasta; // =cip.p1;
// bool print_each; //=false;
}
mehp[0].idx()=0;
mehp[0].source()=&li;
// this takes the one andonly unknown not the first known doh 
//if (make_unknowns_from_first_known) { Fasta foo; IdxTy junk;
//mehp[0].next(foo,junk); // eliminate the references as it messes up 
//}


mehp[0].launch(nt,pfunc,&mehp[0]);


} // cmd_mt_explore_hits
////////////////////////////////////////////








//////////////////////////////////////////


void cmd_explore_hits(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy kfasta=cip.p1;
const StrTy fasta=cip.p2;
const StrTy seqno=cip.wif(3);
const StrTy knono=cip.wif(4);
IdxTy seqnon=bad();
IdxTy knonon=bad();
if (seqno.length()!=0) seqnon=myatoi(seqno);
if (knono.length()!=0) knonon=myatoi(knono);

//const StrTy kfasta=m_flp.knowns_fasta(); // "knowns";
const IdxTy mflags=m_flp.mflags(); // "knowns";
//const bool use_knowns=true;
const bool only_print_muts=((mflags&4)!=0);;
MM_ERR(" enter cmd_explore_hits "<<MMPR3(cmd,fasta,kfasta)<<MMPR3(mflags,int(seqnon),int(knonon)))
//Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
Fasta dummy;
Fasta rev;
Fasta & knowns=m_fasta_map[kfasta]; // (!use_knowns)?dummy:m_fasta_map[kfasta];
modified_fasta(rev ,seqs,mflags);
Fidc fidc;
fidc.index(seqs,true);
Fidc rfidc;
rfidc.index(rev,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))

//Fidc kfidc;
//kfidc.index(knowns,true);

if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
if (kfasta!=m_indexed_fasta)
{
MM_ERR(" wrong fasta infdex want "<<kfasta<<" have "<<m_indexed_fasta)
return; 
}
//const StrTy fasta=m_indexed_fasta;
//Fasta & fts=m_fasta_map[fasta];
//MM_ERR(" searchsing "<<MMPR4(cmd, fts.size(),fasta,seq))
Fidc&  kfidc = *m_fidc;

MM_ERR(" previously  indexed knowns  "<<MMPR3(kfasta,knowns.size(),kfidc.size()))



MM_SZ_LOOP(j,seqs,szs)
{
if (seqnon!=bad()) if (seqnon!=j) continue;
IdxTy kbest=bad();
IdxTy krbest=bad();
IdxTy lfbest=0;
IdxTy lrbest=0;

MM_SZ_LOOP(krep,knowns,szk)
{
if (knonon!=bad()) if (knonon!=krep) continue;
HiVec hitsf,hitsr;
kfidc[krep].alignment_points(hitsf,0,fidc[j]);
kfidc[krep].alignment_points(hitsr,0,rfidc[j]);
const StrTy & n1=seqs.name(j);
const StrTy & n2=knowns.name(krep);
const bool i_vs_j=!true;
IdxTy lenf=char_hits(hitsf,i_vs_j);
IdxTy lenr=char_hits(hitsr,i_vs_j);
// this failes due to overlap lol 
//MM_LOOP(ii,hitsf) lenf+=(*ii).len;
//MM_LOOP(ii,hitsr) lenr+=(*ii).len;
if (lenf>lfbest) { lfbest=lenf; kbest=krep; } 
if (lenr>lrbest) { lrbest=lenr; krbest=krep; } 

if (only_print_muts) if (lenf>lenr) continue;
MM_MSG(MMPR2(j,krep)<<MMPR4(n1,n2,hitsf.size(),hitsr.size())<<MMPR2(lenf,lenr))
//fidc[j].alignment_points(hitsr,0,kfidc[krep]);
MM_SZ_LOOP(i,hitsf,szf)
{
const StrTy&  s=hitsf[i].s;
const int ii=hitsf[i].i;
const int jj=hitsf[i].j;
const IdxTy len=hitsf[i].len;
MM_MSG("fwd "<<MMPR4(mflags, i,s,len)<<MMPR3(ii,jj,(ii-jj)))
}
MM_SZ_LOOP(i,hitsr,szr)
{
const StrTy&  s=hitsr[i].s;
const int  ii=hitsr[i].i;
const int  jj=hitsr[i].j;
const IdxTy len=hitsr[i].len;
MM_MSG("rev "<<MMPR4(mflags,i,s,len)<<MMPR3(ii,jj,(ii-jj)))
}
} // knowns 
const StrTy & n1=seqs.name(j);
const StrTy & n2=knowns.name(kbest);
const StrTy & n2r=knowns.name(krbest);
MM_MSG(" bests "<<MMPR(j)<<MMPR4(lfbest,kbest,lrbest,krbest)<<MMPR3(n1,n2,n2r))


} //j  


//StrTy d1,d2;
//Msi msi;



} // cmd_explore_hits


void cmd_lowp_fasta_hits(Cip & cip , LocalVar & lv ) 
{
const bool debug_result=false; 
const StrTy cmd=cip.cmd();
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
IdxTy plim=myatoi(cip.wif(3));
if (plim==0) plim=4;
const bool print_status=true;
const bool use_knowns=true;
const bool total_genus_species=true;
const bool sort_gs=true;
const bool cnt_only=!true;
const bool print_contains=!total_genus_species;
const StrTy kfasta="knowns";
Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
Fasta dummy;
Fasta & knowns=(!use_knowns)?dummy:m_fasta_map[kfasta];
MM_ERR(" lowp_fasta_hits "<<MMPR4(fasta,seqs.size(),s,fts.size())<<MMPR2(cmd,plim))
Fidc fidc;
fidc.index(seqs,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
Fidc kfidc;
kfidc.index(knowns,true);
MM_ERR(" done indexing knowns  "<<MMPR3(kfasta,knowns.size(),kfidc.size()))
typedef std::map<StrTy, IdxTy> TabMap;
TabMap kfreq,ufreq;
std::vector<StrTy> gsv(knowns.size()); 
StrTy gsname="";
if (total_genus_species) { freqtab(kfreq,knowns,&gsv,cnt_only); } 

//TaxTree & tt=m_tax_tree;
//std::vector<IdxTy> cnt;
IdxTy tot=0;
//if ( tally_hits) cnt.resize(fts.size());
MM_SZ_LOOP(j,seqs,sz)
{
const char * qjname=seqs.name(j).c_str(); // ts[j].c_str();
MM_MSG(" Unknown sequence "<<MMPR(j) <<qjname) 
if (total_genus_species) { ufreq.clear();  } 
IdxTy jhits=0;
MM_SZ_LOOP(i,fts,szs)
{
const char * cq=fts.seq(i).c_str(); // ts[j].c_str();
//const char * n=seqs.name(i).c_str();
const IdxTy hitloc=fidc[j].find_first(cq);
if (hitloc!=bad())
{
++tot; ++jhits;
const char * iname=fts.name(i).c_str(); // ts[j].c_str();
if (print_contains) { MM_MSG( "contains "<<MMPR4(j,i,hitloc,iname))  } 
if (use_knowns)
{
const IdxTy lowp=namefreq(iname);
if (lowp==bad()) continue;
if (lowp==0) continue;
if (lowp>plim) continue;
MM_SZ_LOOP(k,knowns,ksz)
{
const IdxTy khitloc=kfidc[k].find_first(cq);
if (khitloc!=bad())
{
const char * kname=knowns.name(k).c_str(); // ts[j].c_str();
if (total_genus_species) { 
 if (cnt_only) { ++ufreq[gsv[k]];  }
else { ufreq[gsv[k]]+=fts.seq(i).length(); } 

 } 
if (print_contains ) { MM_MSG( "containsknown  "<<MMPR(k)<<MMPR4(j,i,hitloc,kname)) } 

}

} // k 
} // use_knowns 
} // hitloc 
} // i
if (total_genus_species&&!sort_gs)
{
MM_LOOP(ii,ufreq)
{
const StrTy k=(*ii).first;
const D rat=1.0*ufreq[k]/kfreq[k];
MM_MSG(" gssummary "<<MMPR4(j,k, ufreq[k],kfreq[k])<<MMPR2(rat,qjname))
}
}
if (total_genus_species&&sort_gs)
{
std::map<D, StrTy> xxx;
MM_LOOP(ii,ufreq)
{
const StrTy k=(*ii).first;
const D rat=1.0*ufreq[k]/kfreq[k];
xxx[-rat]=k;
//MM_MSG(" gssummary "<<MMPR4(j,k, ufreq[k],kfreq[k])<<MMPR2(rat,qjname))
}
IdxTy slim=10;
MM_LOOP(ii,xxx)
{
const StrTy k=(*ii).second;
const D rat=-(*ii).first;
// k is the genus-species string qjname is the unknown seq name  from the fasta 
MM_MSG(" gssummary "<<MMPR4(j,k, ufreq[k],kfreq[k])<<MMPR2(rat,qjname))
const bool find_to_align=true;
if (find_to_align)
{
IdxTy krep=0;
while (krep<knowns.size() )
{
if (gsv[krep]==k) break;
++krep;
}
MM_ERR(MMPR3(k,krep,knowns.size()))
HiVec hitsf,hitsr;
//const IdxTy guesssz=15000;
//hitsf.reserve(guesssz);
//hitsr.reserve(guesssz);
kfidc[krep].alignment_points(hitsf,0,fidc[j]);
fidc[j].alignment_points(hitsr,0,kfidc[krep]);
StrTy d1,d2;
Msi msi;
MM_ERR( " done finding markers now align "<<MMPR4(krep,j,hitsf.size(),hitsr.size()))
//msi.make_alignment(d1,d2,seqs.seq(j),fts.seq(krep),hitsf,hitsr);
msi.make_alignment(d1,d2,seqs.seq(j),knowns.seq(krep),hitsr,hitsf);
if (debug_result) { MM_ERR(d1);}
if (debug_result) { MM_ERR(d2);}
//const char * nm2=fasta2.name(j).c_str();
typedef Msi::summary_type Mst;
msi.base_map();
Mst x= msi.comp_base_sequences(d1.c_str(),d2.c_str(),0);
MM_MSG(" compare "<<MMPR2(x.eq,x.diff)<<MMPR(x.ncmp))
//dest.add( StrTy("aligned 1") ,d1);
//dest.add( StrTy("aligned 2") ,d2);


} // find_to_align 


--slim;
if (slim==0) break; 
}


}





if (print_status)
 { 
	if ( ( j%10)==0) 
		{ MM_STATUS(" cmd_fasta_hits "<<MMPR2(j,tot)<<"     ") }  
}
} // j

} // lowp_fasta_hits

IdxTy namefreq(const char * iname)
{
// h1722 contains  j=0 i=55679 hitloc=142 iname=>MM_16SCOM_55680 0  cnt 14 len 28 h 1.90682
Ss ss;
ss<<iname;
while (true) 
{ 
	StrTy x; ss>>x; 
	if (x.length()==0) { if( ss.eof()) break;  continue; } 
	if (x=="cnt") { ss>>x;  return myatoi(x); } 
} // while 
return bad();
} // namefreq




static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 
void cmd_blast(Cip & cip , LocalVar & lv ) 
{
//typedef std::map<StrTy, IdxTy> HaveMap;
//typedef  mjm_string_index_collection  Fidc;
//typedef std::vector< string_indexer> Fid;
const StrTy deflin=StrTy("nolineage");
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const bool print_haves=m_flp.print_haves();
const bool print_havenots=m_flp.print_havenots(); // false;
const bool print_counts=m_flp.print_counts(); // !false;
const bool print_if_have=m_flp.print_if_have(); // !false;
const bool suppress_vector=m_flp.suppress_vector(); // !false;
const bool mode=m_flp.add_level(); // !false;
const bool skip_counts=!(print_haves||print_havenots||print_counts);
//const bool make_vector=skip_counts;
const bool print_hit=m_flp.print_hit();
const IdxTy accmode=skip_counts?4:m_flp.accmode();
//const IdxTy maxdepth=m_flp.maxdepth();

MM_ERR(" mmblast "<<MMPR4(s,fasta,accmode,print_counts)<<MMPR4(print_haves,print_havenots,print_if_have,mode)<<MMPR2(suppress_vector,print_hit))

Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];


//std::map<StrTy,IdxTy> m;
Fidc fidc;
fidc.index(seqs,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))

TaxTree & tt=m_tax_tree;
MM_ERR(MMPR4(s,fts.size(),fasta,seqs.size()))
MM_ERR(MMPR(tt.size()))
MM_SZ_LOOP(j,fts,szs)
{
const char * cq=fts.seq(j).c_str(); // ts[j].c_str();
const char * qjname=fts.name(j).c_str(); // ts[j].c_str();
IdxTy best=~0U;
IdxTy ibest=0;
std::vector<IdxTy> errors(seqs.size()); 
MM_SZ_LOOP(i,seqs,sz)
{
//const char * n=seqs.name(i).c_str();
const IdxTy error=fidc[i].find_best_agct2(cq);
errors[i]=error;
if (error<best) { best=error; ibest=i; }
// taking more time than search? 
}
std::vector<IdxTy> all;
all.push_back(ibest);
for(IdxTy i=ibest; i<errors.size(); ++i) if (errors[i]==ibest) all.push_back(i);  

//MM_MSG( " best match "<< MMPR4(j,ibest,best,qjname)<<MMPR(seqs.name(ibest)))
Ss ss; 
MM_SZ_LOOP(i,all,allsz)
{
std::vector<StrTy> name_words;
const char * n = seqs.name(all[i]).c_str();
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
IdxTy taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
//tt.lineage(lineage,taxon);
ss<<MMPR2(taxon ,tt.lineage(taxon));

} // i 
MM_MSG(MMPR4(best,qjname,allsz,ss.str()))
} // j 
} // cmd_blast


/*
Current name format is id len genera... 
head mm_gen_pairs.fasta 
>MMGEN_1 9 Peptoanaerobacter=1,1 Thermoanaerobacter=2,0.0689655
AAAAAAACG
>MMGEN_2 9 Defluviitoga=1,1 Thiomicrospira=1,0.111111
AAAAAAACT
>MMGEN_3 10 Acetobacterium=10,0.909091 Thermoanaerobacterium=1,0.1
AAAAAAAGAC
>MMGEN_4 11 Desulfurella=1,0.166667 Ornatilinea=1,1
AAAAAAAGGAG
>MMGEN_5 11 Cylindrospermum=3,0.333333 Scytonema=1,1
*/
class genera_hit_mark
{
public:
class genera_share
{
public:
genera_share(const StrTy & g, const D & x, const IdxTy & n):m_genus(g),m_frac(x),m_n(n) {}
const StrTy m_genus;
const D m_frac;
const IdxTy m_n;
}; // genera_share
public:
int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

genera_hit_mark(const char * name, const IdxTy & loc, const StrTy & seq)
{
m_string=seq;
m_loc=loc;
// thisis slow and shouldbe cahced 
const IdxTy slen=strlen(name);
char c[slen+1];
c[slen]=0;
memcpy(c,name,slen+1);
IdxTy word=0;
IdxTy sp=0;
StrTy gen="";
D frac=0;
IdxTy form=~0;
IdxTy  n=0;
for(IdxTy i=0; i<slen; ++i)
{
if (word==0) if (c[i]==' ')
{
c[i]=0;
m_name=c;
sp=i+1;
++word;
continue;
}
if (word==1) if (c[i]==' ')
{
//c[i]=0;
while (c[i]==' ') { c[i]=0; ++i; }
m_len=myatoi(c+sp);
if (m_string.length()!=m_len)
{
MM_ERR(" danger will robinson name and sequence len differ  "<<MMPR2(m_len,m_string.length()))

}
sp=i; // +1;
++word;
continue;
}
//>MMGEN_5 11 Cylindrospermum=3,0.333333 Scytonema=1,1
if (word>1) if (c[i]=='=')
{
c[i]=0;
gen=(c+sp);
sp=i+1;
++word;
form=0;
continue;
}
// alt format with genera:hits without the =n,f stuff
if (word>1) if (c[i]==':')
{
c[i]=0;
gen=(c+sp);
sp=i+1;
++word;
form=1;
continue;
}
if (word>1) if (c[i]==',')
{
c[i]=0;
// TODO FIXME NOTE This is NOT the total number of entries just the 
// TOTAL NUMBER OS HITS, n/frac should give total hits... 
n=myatoi(c+sp);
sp=i+1;
++word;
continue;
}
const bool eol=(c[i+1]==0);
if (word>1) if ((c[i]==' ') ||(eol))
{
if ((i>sp)||eol) {
if (c[i+1]!=0) c[i]=0;
if ( form==0) { frac=atof(c+sp);}
else if ( form==1) { n=myatoi(c+sp); frac=1.0; }
else { MM_ERR(" bad form "<<MMPR3(int(form),sp,name)) } 
}
sp=i+1;
++word;
if ((frac>0) &&(n>0))
{ m_g.push_back(genera_share( gen, frac,  n)); }
else{
MM_ERR( " zero uck "<<MMPR4(name,loc,seq,form)<<MMPR2(strlen(name),sp))
}
gen="dup";
frac=0;
n=~0;
continue;
}

} //  i  

} // genera_hit_mark ctor
StrTy to_string() const
{
Ss ss;
ss<<MMPR4(m_loc,m_len,m_name,m_string);
MM_LOOP(ii,m_g) { ss<<MMPR((*ii).m_genus); } 
return ss.str();
}

IdxTy m_loc, m_len;
StrTy m_name, m_string;
std::vector<genera_share> m_g;

};  // genera_hit_mark

class ghm_vector : public  std::vector<genera_hit_mark> 
{

public:
//Ss ss;
void identify(Ss & ss, const char * uname, const StrTy & useq, const bool print=true)
{
const IdxTy ulen=useq.length();
typedef std::vector<D> GenCo;
typedef std::map<StrTy,GenCo> GenCoMap;
GenCoMap gcm;
MM_LOOP(ii,(*this))
{
genera_hit_mark & ghm= (*ii); // .second;
MM_LOOP(jj,ghm.m_g)
{
const D f=(*jj).m_n*(*jj).m_frac;
if (f==0)
{
MM_ERR(" ASSFACK "<<MMPR3(f,(*jj).m_n,(*jj).m_frac))
}
//const D f=(*jj).m_frac;
const StrTy & genus=(*jj).m_genus;
if (gcm[genus].size()==0) gcm[genus]=GenCo(ulen);
GenCo & gco=gcm[genus];
for(IdxTy i=ghm.m_loc; i<(ghm.m_len+ghm.m_loc); ++i)
{
gco[i]+=f;
} // i 
}
} // ii 
typedef std::map<StrTy , D> Scores;
Scores scores;
MM_LOOP(ii,gcm)
{
const GenCo & gc=(*ii).second;
const StrTy & genus=(*ii).first;
D tot=0;
MM_LOOP(jj,gc) {
if ((*jj)!=0) tot+=1;
// tot+=(*jj);

 } 
scores[genus]-=tot;

}

typedef  std::pair<D,StrTy > Pt;
std::vector< std::pair<D,StrTy > > scoress;
MM_LOOP(ii,scores) { scoress.push_back(Pt((*ii).second,(*ii).first)); }
std::sort(scoress.begin(),scoress.end());

IdxTy nlim=0;
//Ss ss;
IdxTy ntoshow=3;
MM_LOOP(ii,scoress)
{
const D sc=(*ii).first;
if (sc==0)
{
MM_ERR(" zeor score wtf "<<MMPR3((*this).size(),(*ii).first,(*ii).second))
IdxTy ie=0; MM_LOOP(ff,(*this)) {MM_ERR(MMPR(ie)<<(*ff).to_string())  ++ie; }
}
ss<<" "<<((*ii).second)<<":"<<(-(*ii).first/ulen);
if ( nlim>=ntoshow) break;
++nlim;
}
while (nlim<ntoshow) { ss<<" xxx:0" ;  ++nlim; } 
if (print ) { MM_MSG(MMPR2(uname,ss.str())) } 

} // identify

}; // ghm_vector
bool pathological(const StrTy & s) const
{
const IdxTy sz=s.length();
const char * c=s.c_str();
for(IdxTy i=0; i<sz; ++i)
{
const char ci=c[i];
if (ci=='A') continue;
if (ci=='C') continue;
if (ci=='G') continue;
if (ci=='T') continue;
return true;
}

return false; 
}
void cmd_massmuck(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fasta=cip.p1;
const bool print_each=false;
const bool group_hits=true;
const bool skip_pathologicals=true;
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
const StrTy kfasta=m_indexed_fasta;
Fasta & seqs=m_fasta_map[kfasta];
Fasta & fts=m_fasta_map[fasta];
MM_ERR(" searchsing "<<MMPR4(cmd, fasta,fts.size(),seqs.size()))
Fidc&  fidc = *m_fidc;

MM_SZ_LOOP(j,seqs,ucksz)
{
const StrTy & seq=seqs.seq(j);
if (skip_pathologicals) if (pathological(seq)) continue;
//const char * cq=seq.c_str();
const char * n=seqs.name(j).c_str();
typedef ghm_vector Ghmv;
Ghmv ghmv;
//MM_MSG(" shat "<<MMPR2(j,n))
MM_SZ_LOOP(i,fts,szs)
{
const StrTy & seqi=fts.seq(i);
const char * ci=seqi.c_str();
//MM_ERR(" trying to find  uck "<<MMPR(cq))
const IdxTy hitloc=fidc[j].find_first(ci);
if (hitloc!=bad())
{
//const char * qjname=seqs.name(j).c_str(); // ts[j].c_str();
const char * ijname=fts.name(i).c_str(); // ts[j].c_str();
//MM_MSG(" uck "<<MMPR3(hitloc,n,qjname))
if ( print_each) { MM_MSG(" uck "<<MMPR3(hitloc,ijname,n)) } 
if (group_hits) { 
genera_hit_mark ghm(ijname,hitloc,seqi);
ghmv.push_back(ghm);
}
} // hit
} // i 
if (group_hits) {Ss ss;  ghmv.identify(ss,n,seq,true); }
} // j 
} // massmuck 


class mt_massmuck_param
{
public:
// https://www.cs.cmu.edu/afs/cs/academic/class/15492-f07/www/pthreads.html
static pthread_mutex_t*  mutex1()//  { = PTHREAD_MUTEX_INITIALIZER;
{
static pthread_mutex_t _mutex1 = PTHREAD_MUTEX_INITIALIZER;
return & _mutex1;
}
void enter_serial() { pthread_mutex_lock( mutex1() ); }
void exit_serial() { pthread_mutex_unlock( mutex1() ); }
void launch(const IdxTy nthreads, void *(*thread_function)(void* ) ,mt_massmuck_param *msgp)
{
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, thread_function, msgp+i ); }
   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}
Myt * p;
IdxTy nt,i;

 StrTy fasta; // =cip.p1;
 bool print_each; //=false;
 bool group_hits;// =true;
 bool skip_pathologicals; // =true;

}; // mt_massmuck_param





static void* mt_massmuck(void* param ) 
{
mt_massmuck_param * mafp=(mt_massmuck_param*)param;
Myt * p = (*mafp).p;
const StrTy fasta=(*mafp).fasta; // cip.p1;
const bool print_each=(*mafp).print_each; // false;
const bool group_hits=(*mafp).group_hits; // true;
const bool skip_pathologicals=(*mafp).skip_pathologicals; // true;
const StrTy kfasta=(*p).m_indexed_fasta;
Fasta & seqs=(*p).m_fasta_map[kfasta];
Fasta & fts=(*p).m_fasta_map[fasta];
//MM_ERR(" searchsing "<<MMPR4(cmd, fasta,fts.size(),seqs.size()))
Fidc&  fidc = *(*p).m_fidc;
const IdxTy ucksz=seqs.size();
//MM_SZ_LOOP(j,seqs,ucksz)
IdxTy j=(*mafp).i;
const IdxTy del=(*mafp).nt;
for(; j<ucksz; j+=del)
{
const StrTy & seq=seqs.seq(j);
if (skip_pathologicals) if ((*p).pathological(seq)) continue;
//const char * cq=seq.c_str();
const char * n=seqs.name(j).c_str();
typedef ghm_vector Ghmv;
Ghmv ghmv;
//MM_MSG(" shat "<<MMPR2(j,n))
MM_SZ_LOOP(i,fts,szs)
{
const StrTy & seqi=fts.seq(i);
const char * ci=seqi.c_str();
//MM_ERR(" trying to find uck "<<MMPR(cq))
const IdxTy hitloc=fidc[j].find_first(ci);
if (hitloc!=bad())
{
//const char * qjname=seqs.name(j).c_str(); // ts[j].c_str();
const char * ijname=fts.name(i).c_str(); // ts[j].c_str();
//MM_MSG(" uck "<<MMPR3(hitloc,n,qjname))
if ( print_each) { MM_MSG(" uck "<<MMPR3(hitloc,ijname,n)) } 
if (group_hits) { 
genera_hit_mark ghm(ijname,hitloc,seqi);
ghmv.push_back(ghm);
}
} // hit
} // i 
if (group_hits) {Ss ss;  ghmv.identify(ss,n,seq,false);
const char * uname=n;
(*mafp).enter_serial();
 { MM_MSG(MMPR2(uname,ss.str())) } 
(*mafp).exit_serial();


 }
} // j uck 

return 0;
}

void cmd_mt_massmuck(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fasta=cip.p1;
const bool print_each=false;
const bool group_hits=true;
const bool skip_pathologicals=true;
if (m_fidc==0 )  { MM_ERR(" use index-fasta first ") return; }
const IdxTy nt=3;

auto pfunc= Myt:: mt_massmuck; // (mt_shrink_group_param * msgp)
mt_massmuck_param mgsgp[nt];
for (IdxTy i=0; i<nt; ++i)
{
mt_massmuck_param * p= mgsgp+i;
(*p).nt=nt; (*p).i=i; (*p).p=this; 
(*p).fasta=fasta; // cip.p1;
(*p).print_each=print_each; // false;
(*p).group_hits=group_hits; // true;
(*p).skip_pathologicals=skip_pathologicals; // true;

}
mgsgp[0].launch(nt,pfunc,&mgsgp[0]);


/*
const StrTy kfasta=m_indexed_fasta;
Fasta & seqs=m_fasta_map[kfasta];
Fasta & fts=m_fasta_map[fasta];
MM_ERR(" searchsing "<<MMPR4(cmd, fasta,fts.size(),seqs.size()))
Fidc&  fidc = *m_fidc;
const IdxTy ucksz=seqs.size();
//MM_SZ_LOOP(j,seqs,ucksz)
IdxTy j=0;
const IdxTy del=1;
for(; j<ucksz; j+=del)
{
const StrTy & seq=seqs.seq(j);
if (skip_pathologicals) if (pathological(seq)) continue;
//const char * cq=seq.c_str();
const char * n=seqs.name(j).c_str();
typedef ghm_vector Ghmv;
Ghmv ghmv;
//MM_MSG(" shat "<<MMPR2(j,n))
MM_SZ_LOOP(i,fts,szs)
{
const StrTy & seqi=fts.seq(i);
const char * ci=seqi.c_str();
//MM_ERR(" trying to find uck "<<MMPR(cq))
const IdxTy hitloc=fidc[j].find_first(ci);
if (hitloc!=bad())
{
//const char * qjname=seqs.name(j).c_str(); // ts[j].c_str();
const char * ijname=fts.name(i).c_str(); // ts[j].c_str();
//MM_MSG(" uck "<<MMPR3(hitloc,n,qjname))
if ( print_each) { MM_MSG(" uck "<<MMPR3(hitloc,ijname,n)) } 
if (group_hits) { 
genera_hit_mark ghm(ijname,hitloc,seqi);
ghmv.push_back(ghm);
}
} // hit
} // i 
if (group_hits) { ghmv.identify(n,seq); }
} // j uck 

*/
} // mt_massmuck 






////////////////////////////////////////////////////////////////


void cmd_distinctives(Cip & cip , LocalVar & lv ) 
{
//typedef  mjm_string_index_collection  Fidc;
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const IdxTy minlen=10;
TestStrings & ts= m_queries[s];
Fasta & seqs=m_fasta_map[fasta];
MM_ERR(" distinctives "<<MMPR4(fasta,seqs.size(),s,ts.size())<<MMPR(flags))
Fidc fidc;
fidc.index(seqs,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
TaxTree & tt=m_tax_tree;
MM_ERR(MMPR4(s,ts.size(),fasta,seqs.size()))
MM_ERR(MMPR(tt.size()))
MM_SZ_LOOP(i,seqs,sz)
{ // only do once doh 
	const char * n=seqs.name(i).c_str();
	const char * seq=seqs.seq(i).c_str();
	const IdxTy slen=strlen(seq);
	char hits[slen+1];
	for (IdxTy i=0; i<=slen; ++i) hits[i]=0;
	MM_SZ_LOOP(j,ts,szs)
	{
		const char * cq=ts[j].c_str();
		const IdxTy qlen=strlen(cq);
		IdxTy last=0;
		while (true) { 
		const IdxTy hitloc=fidc[i].find_first(cq+last);
		if (hitloc==fidc.bad()) break;
		const IdxTy endh=hitloc+last+qlen;
		const IdxTy lend=(endh<slen)?endh:slen;
		for( IdxTy k=(last+hitloc); k<lend; ++k)
		{
			hits[k]='1';
		}
		// TODO FIXME should worry about overlap but not now 
		last=endh; // hitloc+qlen;
		if (last>=slen) break; 
		} // while (hitloc!=(~0U));
	} // j query
	std::vector<StrTy>  remnants;
	std::vector<IdxTy>  remloc;
	StrTy x="";
	for(IdxTy k=0; k<slen; ++k) { if (hits[k]!=0) {if (x.length()>=minlen) {remloc.push_back(k-x.length()) ; remnants.push_back(x); x=""; }}
else x=x.append(1,seq[k]); }
if (x.length()>=minlen){remloc.push_back(slen-x.length());  remnants.push_back(x);  } 
MM_SZ_LOOP(ir,remnants,rsz) 
{const IdxTy & pos=remloc[ir]; 
const IdxTy & len=remnants[ir].length(); 
const StrTy & s=remnants[ir]; 
MM_MSG(MMPR3(i,ir,len)<<MMPR3(pos,s,n))}
} // i

} // cmd_distinctives






///////////////////////////////////////////////////////////////


void cmd_query_coverage(Cip & cip , LocalVar & lv ) 
{
//typedef  mjm_string_index_collection  Fidc;
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const StrTy ref=(cip.wif(3));
const IdxTy flags=myatoi(cip.wif(4));
const IdxTy level_depth=myatoi(cip.wif(5));
// dmel? TODO FIXME 
if ((level_depth==0) ||(level_depth>10)) 
	{MM_ERR( " level depth seems odd at "<<MMPR(level_depth))}
//const IdxTy level_depth=4;
const D fracmin=.5;
const bool only_print_singles=((flags&1)!=0);
const bool print_all=!only_print_singles;
TestStrings & ts= m_queries[s];
Fasta & seqs=m_fasta_map[fasta];
//std::map<StrTy,IdxTy> 
// TODO FIXME const Mii & mref=m_luts[ref] ;
Mii & mref=m_luts[ref] ;
MM_ERR(" query_coverage "<<MMPR2(level_depth,fracmin)<<MMPR4(fasta,seqs.size(),ref,mref.size())<<MMPR3(flags,only_print_singles,print_all))
Fidc fidc;
fidc.index(seqs,true);
MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
TaxTree & tt=m_tax_tree;
MM_ERR(MMPR4(s,ts.size(),fasta,seqs.size()))
MM_ERR(MMPR(tt.size()))
MM_SZ_LOOP(j,ts,szs)
{
	const char * cq=ts[j].c_str();
	Mii  m;
	m.clear(); // this wil eventually be a m_luts eleement/
	IdxTy hitnode=~0U;
	MM_SZ_LOOP(i,seqs,sz)
	{ // only do once doh 
		const char * n=seqs.name(i).c_str();
		const IdxTy hitloc=fidc[i].find_first(cq);
		bool hit=( hitloc!=fidc.bad());
		if (!hit) continue;
		IdxTy taxon=~0U;
		std::vector<IdxTy> lineage;
		std::vector<StrTy> name_words;
		cip.m_li.parse_full(name_words, n,' ' );
		const IdxTy nwsz=name_words.size();
		taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
		tt.lineage(lineage,taxon);
		MM_SZ_LOOP(k,lineage,lnsz) {const IdxTy & lk=lineage[k];  ++m[lk];}
		// TODO FIXME this seems stupib but this can abort early if crosses taxon
		if (lnsz>level_depth) 
		{ const IdxTy & rnode=lineage[lnsz-level_depth];   if (hitnode==~0U) hitnode=rnode; 
			else if ( rnode !=hitnode) { m.clear(); break; } 
		}	
	} // i seq 
	IdxTy hits=0;
	IdxTy hnode=0;
	IdxTy cnt=0;
	MM_LOOP(ii,m)
	{
		const IdxTy node=(*ii).first;
		const IdxTy coverage =(*ii).second;
		const IdxTy depth=tt.depth(node);
		if (only_print_singles) { if (depth==level_depth){  ++hits; hnode=node;cnt=coverage;  if (hits>1) break; }  continue; } 
		if (print_all) { MM_MSG(MMPR(j)<<MMPR4(node,coverage, depth, tt.node_name(node)))}
	} // ii 
	if (only_print_singles) { if (hits==1)
{
	const StrTy nname=tt.node_name(hnode);
	// TODO make this const or something to stop crapping up memory 
	const IdxTy & total=mref[hnode];
	const D frac=(total!=0)?(1.0*cnt/total):1;
 	if (frac>fracmin) { MM_MSG(MMPR4(j,hnode,cnt,nname)<<MMPR2(frac,total))  }

} 
}

	MM_STATUS(" query "<< MMPR2(j,m.size())<<"     " ) 


} // j query

} // cmd_query_coverage






void cmd_hit_or_miss(Cip & cip , LocalVar & lv ) 
{
typedef std::map<StrTy, IdxTy> HaveMap;
//typedef  mjm_string_index_collection  Fidc;
//typedef std::vector< string_indexer> Fid;
const StrTy deflin=StrTy("nolineage");
const StrTy s=cip.p1;
const StrTy fasta=cip.p2;
const bool print_haves=m_flp.print_haves();
const bool print_havenots=m_flp.print_havenots(); // false;
const bool print_counts=m_flp.print_counts(); // !false;
const bool print_if_have=m_flp.print_if_have(); // !false;
const bool suppress_vector=m_flp.suppress_vector(); // !false;
const bool mode=m_flp.add_level(); // !false;
const bool skip_counts=!(print_haves||print_havenots||print_counts);
const bool make_vector=skip_counts;
const IdxTy accmode=skip_counts?4:m_flp.accmode();
const IdxTy maxdepth=m_flp.maxdepth();

MM_ERR(" hit_or_miss "<<MMPR4(s,fasta,accmode,print_counts)<<MMPR4(print_haves,print_havenots,print_if_have,mode)<<MMPR(suppress_vector))
//typedef std::vector<StrTy>  Vec;

//Vec vec;
CharMat &  vec=m_char_mat;
vec.clear();

TestStrings & ts= m_queries[s];
//Fasta & fts=m_fasta_map[s];
Fasta & seqs=m_fasta_map[fasta];
//if (make_vector) vec.reserve(seqs.size());
//if (make_vector) vec=Vec(seqs.size());
if (make_vector) vec.size(seqs.size(),ts.size());
std::map<StrTy,IdxTy> m;
Fidc fidc;
fidc.index(seqs,true);

MM_ERR(" done indexing "<<MMPR3(fasta,seqs.size(),fidc.size()))
#if 0 
Ragged lineages;
MM_SZ_LOOP(i,seqs,sz)
{
const char * n=seqs.name(i).c_str();
std::vector<StrTy> name_words;
std::vector<StrTy> lineage;
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
const IdxTy taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
tt.lineage(lineage,taxon);
lineage.push_back(lineage);
} // i 
MM_ERR(" done finding lineages  "<<MMPR(lineages.size()))
#endif

TaxTree & tt=m_tax_tree;
MM_ERR(MMPR4(s,ts.size(),fasta,seqs.size()))
MM_ERR(MMPR(tt.size()))
MM_SZ_LOOP(j,ts,szs)
{
HaveMap have,havenot;
std::map<StrTy, IdxTy> counts;
const char * cq=ts[j].c_str();
MM_SZ_LOOP(i,seqs,sz)
{
//const char * c=seqs.seq(i).c_str();
const char * n=seqs.name(i).c_str();
IdxTy taxon=~0U;
// TODO this is slow, do just once  but then memory and tokenizing etc 
#if 1 
std::vector<StrTy> lineage;
if (!skip_counts) { 
std::vector<StrTy> name_words;
cip.m_li.parse_full(name_words, n,' ' );
const IdxTy nwsz=name_words.size();
taxon=(nwsz>0)?myatoi(name_words[nwsz-1]):(~0U);
tt.lineage(lineage,taxon);
}
#else
const std::vector<StrTy>&  lineage=lineages[i];

#endif
if (lineage.size()==0) lineage.push_back(deflin); 
const IdxTy hitloc=fidc[i].find_first(cq);
bool hit=( hitloc!=fidc.bad());
HaveMap & mr =(hit)?have:havenot;
//IdxTy depth=0;
//TODO FIXME the lineage entries need to have whitespace removed 
switch (accmode)
{
case 0:{MM_SZ_LOOP(ii,lineage,szlin)  // all levels
		//{ ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{ucnts(mr,counts,lineage,ii,mode); } break; } //  ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
case 1:{MM_SZ_LOOP(ii,lineage,szlin)  // just top few 
		//{if ((ii+maxdepth)<szlin) continue; ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{if ((ii+maxdepth)<szlin) continue; ucnts(mr,counts,lineage,ii,mode);}break;} 
case 2:{MM_SZ_LOOP(ii,lineage,szlin)  // just phlya and equivalent s
		//{if ((ii+3)<szlin) continue; ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
		{if ((ii+3)<szlin) continue; ucnts(mr,counts,lineage,ii,mode);}break;} 
case 3:{MM_SZ_LOOP(ii,lineage,szlin)  // terminal name only 
		//{ ++mr[lineage[ii]];++counts[lineage[ii]];break; }break;} 
		{ucnts(mr,counts,lineage,ii,mode); ;break; }break;} 

case 4: { break; } // null mode 
default: MM_ERR(" unknown lineage counting mode "<<MMPR(accmode)) 
};
// taking more time than search? 
if (false) { MM_STATUS(MMPR4(i,have.size(),havenot.size(),hit)<<MMPR2(taxon,lineage.size())<<"  .. ") } 
if (hit){ } else { }
//MM_MSG(MMPR4(cq,c,n,hit))
//if (make_vector){  vec[i]+=hit?"1":"0";  } 
if (make_vector){  vec(i,j)=hit?'1':'0';  } 
} // i 
if (print_haves) { MM_LOOP(ii,have) { MM_MSG(" have "<<MMPR2(j,cq)<<MMPR2((*ii).first,(*ii).second)) } } 
if (print_havenots) { MM_LOOP(ii,havenot) { MM_MSG(" havenot "<<MMPR2(j,cq)<<MMPR2((*ii).first,(*ii).second)) } }
if (print_counts)
{
MM_LOOP(ii,counts)
{
const StrTy& name=(*ii).first; 
const D total=(*ii).second;
const IdxTy h=have[name];
const IdxTy hn=havenot[name];
const bool pr=( (print_if_have&&(h!=0) )|| !print_if_have);
if (pr) { const D del=D(h)-D(hn); MM_MSG(" counts "<<MMPR2(j,cq)<<MMPR4(name,total,h,hn)<<MMPR((del/total))) }
}
} // print_counts
if (!false) { MM_STATUS(MMPR3(j,have.size(),havenot.size())<<"  .. ") } 

} // j 
if (!suppress_vector)
{
//MM_SZ_LOOP(i,vec,sz)
for (IdxTy i=0; i<vec.rows(); ++i)
{
const char * n=seqs.name(i).c_str();
const IdxTy sumone=sum_one(vec[i]);
MM_MSG(MMPR4(i,n,vec[i],sumone))
}
} // suppressvector
} // cmd_hit_or_miss




template <class Tx, class Ty,class Tz>
void ucnts(Tx & mr, Ty & counts, Tz& lineage, const IdxTy ii, const bool lbl) const
{
//{ ++mr[lineage[ii]];++counts[lineage[ii]];}break;} 
 StrTy x=lineage[ii];
if (lbl) { Ss ss; ss<<ii<<" "<<x; x=ss.str(); } 
 	++mr[x];
	++counts[x]; 
}

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("quiet")]=&Myt::cmd_quiet;
m_cmd_map[StrTy("verbose")]=&Myt::cmd_verbose;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;
m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
m_cmd_map[StrTy("copy-fasta")]=&Myt::cmd_copy_fasta;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("script")]=&Myt::cmd_script;
m_cmd_map[StrTy("stream-script")]=&Myt::cmd_stream_script;
m_cmd_map[StrTy("stream-fasta")]=&Myt::cmd_stream_fasta;
m_cmd_map[StrTy("list-fasta")]=&Myt::cmd_list_fasta;
m_cmd_map[StrTy("list-align")]=&Myt::cmd_list_align;
//void cmd_list_algn(Cip & cip , LocalVar & lv ) 
//void cmd_write_align(Cip & cip , LocalVar & lv ) 
m_cmd_map[StrTy("read-interleaved-fasta")]=&Myt::cmd_read_interleaved_fasta;
m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
m_cmd_map[StrTy("write-align")]=&Myt::cmd_write_align;
m_cmd_map[StrTy("set-annotation")]=&Myt::cmd_set_annotation;
m_cmd_map[StrTy("write-drift")]=&Myt::cmd_write_drift;
m_cmd_map[StrTy("cmd-align")]=&Myt::cmd_cmd_align;
m_cmd_map[StrTy("index-fasta")]=&Myt::cmd_index_fasta;
m_cmd_map[StrTy("write-interleaved-fasta")]=&Myt::cmd_write_interleaved_fasta;
m_cmd_map[StrTy("compare-fasta")]=&Myt::cmd_cmp_fasta;
m_cmd_map[StrTy("align-fasta")]=&Myt::cmd_align_fasta;
m_cmd_map[StrTy("enumerate-fasta")]=&Myt::cmd_enumerate_fasta;
m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("dump-ragged")]=&Myt::cmd_dump_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;
m_cmd_map[StrTy("modify-fasta-names")]=&Myt::cmd_modify_fasta_names;
m_cmd_map[StrTy("read-queries")]=&Myt::cmd_read_queries;
m_cmd_map[StrTy("vector-search")]=&Myt::cmd_vector_search;
m_cmd_map[StrTy("count-uniq")]=&Myt::cmd_count_uniq;
m_cmd_map[StrTy("find-strings")]=&Myt::cmd_find_strings;
m_cmd_map[StrTy("hit-or-miss")]=&Myt::cmd_hit_or_miss;
m_cmd_map[StrTy("hit-or-miss-fasta")]=&Myt::cmd_hit_or_miss_fasta;
m_cmd_map[StrTy("check-fasta-taxon")]=&Myt::cmd_check_fasta_taxon;
m_cmd_map[StrTy("seq-hits")]=&Myt::cmd_seq_hits;
m_cmd_map[StrTy("list-hits")]=&Myt::cmd_seq_hits;
m_cmd_map[StrTy("write-svg")]=&Myt::cmd_write_svg;
m_cmd_map[StrTy("index-analysis")]=&Myt::cmd_index_analysis;
m_cmd_map[StrTy("blast")]=&Myt::cmd_blast;
m_cmd_map[StrTy("coverage")]=&Myt::cmd_coverage;
m_cmd_map[StrTy("descend")]=&Myt::cmd_coverage;
m_cmd_map[StrTy("query-coverage")]=&Myt::cmd_query_coverage;
m_cmd_map[StrTy("print-node-lut")]=&Myt::cmd_print_node_lut;
m_cmd_map[StrTy("distinctives")]=&Myt::cmd_distinctives;
m_cmd_map[StrTy("fasta-hits")]=&Myt::cmd_fasta_hits;
m_cmd_map[StrTy("lowp-fasta-hits")]=&Myt::cmd_lowp_fasta_hits;
m_cmd_map[StrTy("explore-hits")]=&Myt::cmd_explore_hits;
m_cmd_map[StrTy("mt-explore-hits")]=&Myt::cmd_mt_explore_hits;
m_cmd_map[StrTy("mt-explore-streaming-hits")]=&Myt::cmd_mt_explore_streaming_hits;
m_cmd_map[StrTy("make-fasta")]=&Myt::cmd_make_fasta;
m_cmd_map[StrTy("fasta-commons")]=&Myt::cmd_fasta_hits;
m_cmd_map[StrTy("fasta-summary")]=&Myt::cmd_fasta_hits;
//void cmd_seg_hits(Cip & cip , LocalVar & lv ) 
m_cmd_map[StrTy("seg-hits")]=&Myt::cmd_seg_hits;
m_cmd_map[StrTy("distro-commons")]=&Myt::cmd_distro_commons;
//void cmd_genus_disc(Cip & cip , LocalVar & lv ) 
m_cmd_map[StrTy("genus-disc")]=&Myt::cmd_genus_disc;
//void cmd_seg_grow(Cip & cip , LocalVar & lv ) 
m_cmd_map[StrTy("seg-grow")]=&Myt::cmd_seg_grow;
m_cmd_map[StrTy("grow-genus")]=&Myt::cmd_grow_genus;
m_cmd_map[StrTy("grow-digenus")]=&Myt::cmd_grow_genus;
m_cmd_map[StrTy("grow-trigenus")]=&Myt::cmd_grow_genus;
m_cmd_map[StrTy("grow-ngenus")]=&Myt::cmd_grow_genus;
m_cmd_map[StrTy("shrink-group")]=&Myt::cmd_shrink_group;
m_cmd_map[StrTy("modify-lines")]=&Myt::cmd_modify_lines;
m_cmd_map[StrTy("massmuck")]=&Myt::cmd_massmuck;
m_cmd_map[StrTy("mt-massmuck")]=&Myt::cmd_mt_massmuck;
m_cmd_map[StrTy("mt-shrink-group")]=&Myt::cmd_mt_shrink_group;
m_cmd_map[StrTy("random-scores")]=&Myt::cmd_random_scores;
//void cmd_random_scores(Cip & cip , LocalVar & lv ) 
//void cmd_mt_shrink_group(Cip & cip , LocalVar & lv ) 
//void cmd_grow_genus(Cip & cip , LocalVar & lv ) 
//void cmd_make_fasta(Cip & cip , LocalVar & lv ) 
//void cmd_distro_commons(Cip & cip , LocalVar & lv ) 

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
TaxTree & tt = m_tax_tree;
StrTy local_label="tat";
//typedef void (Tsrc::* TargCmd)( ListTy & choices,  const char * frag);
//typedef void (Tsrc::* TargParam)( ListTy & choices, const char *  cmd, const char * frag);
m_cli.set_target(*this);
//void set_command_handler(TargCmd * p ) { m_targ_cmd=p; }
//void set_param_handler(TargParam * p ) { m_targ_param=p; }
m_cli.set_command_handler(&Myt::cli_cmd);
m_cli.set_param_handler(&Myt::cli_param);
//std::vector<StrTy> some;
//some.push_back(StrTy("load-tax-nodes"));
//m_cli.set_list(some);
m_cli.activate();
LocalVar mloc;
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
if (m_flp.log_commands()) 
		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd.c_str()[0]=='#' ) continue; 

if (m_cmd_map.find(cmd)!=m_cmd_map.end())
{
 CommandInterpretterParam  cip(li); 
(this->*m_cmd_map[cmd])(cip,mloc);
continue;

}
const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
if (cmd=="about") { about();  continue; } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="print") { MM_MSG(li.line())  continue; } 
if (cmd=="err") { MM_ERR(li.line())  continue; } 
if (cmd=="status") { MM_STATUS(li.line())  continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 

//m_tax_tree.standard_commnds(cmd,p1,p2,li);


if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }

/*
//if (cmd=="reset-solution") { m_save_sol=!true;  continue; } 
*/

if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");
if (m_flp.exit_on_err())
{
MM_ERR(" quiting "<<MMPR(m_flp.exit_on_err()))
clean_up();
return; 

}
}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
void about()
{
Ss ss;
ss<<" mjm_trees_and_tables "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com 2018-01-09 "<<CRLF;
ss<<" Code to read various files related to 16S rRNA analyses "<<CRLF;
ss<<"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837139/#SM10"<<CRLF;
ss<<"database link provided by Bradley.stevenson@ou.edu"<<CRLF;
ss<<"http://www.earthmicrobiome.org/data-and-code/ "<<CRLF;
ss<<"Sample data provided by echen@zymoresearch.com "<<CRLF;
ss<<"http://www.isppweb.org/smc_files/bull%20et%20al.%202010%20jpp%20list.pdf"<<CRLF;
ss<<"http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.690.647&rep=rep1&type=pdf"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/Pathogenic_bacteria"<<CRLF;
ss<<"https://en.wikipedia.org/wiki/List_of_bacteria_genera"<<CRLF;

std::ostream & os=std::cout;
//os<<ss;
os<<ss.str();

}

// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code.. uck 
li.push(); // shift commands and remove invoking one, 
// use CommandInterpretterParam 
//const StrTy & cmd=li.word(0);

//if (cmd=="int") { test_int(li);} //   continue; } 
//else if (cmd=="diff") { test_diff(li);} //   continue; } 
//else if (cmd=="fl") { test_fl(li);} //   continue; } 
//void test_gr_integrator( CommandInterpretter & li )
if (false) {}else { MM_ERR(" unrecignized TEST command "<<li.dump()) } 

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
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}
bool done() const  { return m_done; } 

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general
/*

StrTy tod_from_string(const StrTy & w)
{
if (w.length()!=6) return StrTy("BAD")+w;
return w.substr(0,4);
}

*/



void parse_dmel( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
// do nothing with this for now 
		if (m_dmel!=0) 
		{ (*m_dmel).event("userdmelentry" ,name, d,w); } 

}




/////////////////////////////////////////////////////////////////////
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
void Init()
{
//delete m_dmel;
//m_dmel= new Dmel();
m_done=false;
//load_literals();
//load_shares();
//load_handlers();
}

/*
void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<" "<<code<<CRLF;
    ss.str(StrTy(""));
}

*/

void DMel(const StrTy & e)
{
MM_ERR(e)
if (m_dmel!=0) {m_dmel->event("wtf",e); }
}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
	if (m_dmel!=0) {m_dmel->event("wtf",ss.str()); }
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; }
	if (m_dmel!=0) {m_dmel->event(code,ss.str()); }
    ss.str(StrTy(""));
}
void DumpDMel() // (OsTy & os)
{
if (m_dmel!=0) { MM_ERR(MMPR((m_dmel->string(1)))) } 
else { MM_ERR(" m_dmel is null ") }

}


// MEMBERS

bool m_done;

Logic m_flp;
VariableStore m_variables;
Dmel * m_dmel;
Ss m_ss;
FastaMap m_fasta_map;
RaggedMap m_ragged_map;
AgcMap m_aligns_map;
TestStringMap m_queries;
CharMat m_char_mat;
MiiMap m_luts;
TaxTree m_tax_tree; // now have sequences ID's to taxon 
CounterMap m_cm;
CliTy m_cli;
std::map<StrTy,Fidc> m_fidc_map;
Fidc * m_fidc;
StrTy m_indexed_fasta;
std::vector<bool> m_loudness;
}; //mjm_trees_and_tables


#ifdef PYTHON_BOOST_BUILD
char const* greet() { return "hello, world"; }
using namespace boost::python;
static mjm_string_seq * fick() { 
static mjm_string_seq * p = new mjm_string_seq();
return p;
}

class fick_python
{
public:
//fick_python() {}
//fick_python() {m_fick=new mjm_string_seq(); }
fick_python() {m_fick=fick(); }
//~fick_python() {delete m_fick;}
~fick_python() {}
void FACKBOOST() { m_fick->FACKBOOST(); } 
static void FACKBOOSTs() { fick()->FACKBOOST(); } 
static void cmd_str(const StrTy & s ) { fick()->command_mode(s); } 
mjm_string_seq *  m_fick; 
}; // fick_python
//    class_<fick_python>("mjm_string_seq")
BOOST_PYTHON_MODULE(mjm_string_seq)
{
        //.def("command_mode", &fick_python::FACKBOOSTs)
      //  boost::python::def("command_mode", &fick_python::FACKBOOSTs) ;
        boost::python::def("cmd", &fick_python::cmd_str) ;
}


#endif

/////////////////////////////////////////////////////////

#ifdef  TEST_STRING_SEQ__
int main(int argc,char **args)
{
typedef mjm_string_seq  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

