
This appears to work but in some cases intermediate names in
Sarcina changed requiring summing over g and s only with "add_identicals",

 2335  ./mjm_zymo.out  -cmd "bioms-to-ssv assdf r 1 1 data/2018-03-09/otu_table.biom data/2022-08-19/ASV_Table.biom" -cmd "set-param crp_n1 1" -cmd "set-param crp_n2 4" -cmd "set-param crp_ngroup 6" -cmd "count-key-select xx6 r 2" -cmd "quit" 2>&1 | tee fck


Reconciling sequences is hindered by trimming differences.
Here the names need to match in bioms-to-ssv in order touse the
sequences rather than pipeline dependent assignments. If flag 1<<3
is set, loaded assignment files are used. Apparently the loaded name
needs to match the count biom name exactly.   

2002  ./mjm_zymo.out -cmd "load-recon br data/2022-08-19/sv.seqs.fna 1" -cmd "load-recon y data/2022-04-19/zr6409_sv.seqs.fna 1" -cmd "load-recon z data/2022-07-11/zr6948_sv.seqs.fna 1"  -cmd "bioms-to-ssv assdf r 0 9 br y z" -cmd "set-param crp_n1 0" -cmd "set-param crp_n2 2" -cmd "set-param crp_ngroup 8" -cmd "count-key-select bb7 r 2" -cmd "quit" 2>&1 | tee fck


./mjm_zymo.out  -cmd "load-recon 2018 data/2018-03-09/zr2097_sv.seqs.fna 1"  -cmd "load-recon 202208 data/2022-08-19/sv.seqs.fna 1"  -cmd "bioms-to-ssv assdf r 0 9 2018 202208" -cmd "set-param crp_n1 1" -cmd "set-param crp_n2 4" -cmd "set-param crp_ngroup 10" -cmd "count-key-select xxr7 r 2" -cmd "quit" 2>&1 | tee fck

This left out Brownie, 
./mjm_sequence_reconcile.out  "load x data/2022-02-18/zr5958_sv.seq.fna 1 " "load-ass x data/2022-02-18/zr5958_tax_assignments.txt 0"  "load y data/2022-04-19/zr6409_sv.seqs.fna 1" "load-ass y data/2022-04-19/zr6409_tax_assignments.txt 0" "load a data/2018-02-16/rep_set.seqs.fna 1" "load-ass a data/2018-02-16/tax_assignments.txt 0" "load b data/2018-03-09/zr2097_sv.seqs.fna 1" "load-ass b data/2018-03-09/tax_assignments.txt" "load z data/2022-07-11/zr6948_sv.seqs.fna 1" "load-ass z data/2022-07-11/zr6948_tax_assignments.txt 0"  "dump-ass-fasta allanno.fasta 1" quit 2>&1 | tee fck 

 Include Brownie but leave out ITS ,

 2406  ./mjm_sequence_reconcile.out  "load x data/2022-02-18/zr5958_sv.seq.fna 1 " "load-ass x data/2022-02-18/zr5958_tax_assignments.txt 0"  "load y data/2022-04-19/zr6409_sv.seqs.fna 1" "load-ass y data/2022-04-19/zr6409_tax_assignments.txt 0" "load a data/2018-02-16/rep_set.seqs.fna 1" "load-ass a data/2018-02-16/tax_assignments.txt 0" "load b data/2018-03-09/zr2097_sv.seqs.fna 1" "load-ass b data/2018-03-09/tax_assignments.txt" "load z data/2022-07-11/zr6948_sv.seqs.fna 1" "load-ass z data/2022-07-11/zr6948_tax_assignments.txt 0" "load br data/2022-08-19/sv.seqs.fna 1" "load-ass br data/2022-08-19/tax_assignments.txt 0"  "dump-ass-fasta allanno2.fasta 1" quit 2>&1 | tee fck 

