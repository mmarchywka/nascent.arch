
 3331  history | grep stream-srcipt >> README.txt 
 3315  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep "W..W..W..C\|>"| grep -v -B 1  ">" | highlight -red W..W..W..C |  more
 3316  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep "W..W..W..C\|>"| grep -v -B 1  ">" | grep ">" | wc
 3318  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep "W..W..W..C\|>"| grep -v -B 1  ">" | highlight -red W..W..W..C |  more
 3320  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>"| grep -v -B 1  ">" | highlight -red W..W..W..C |  more
 3321  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  more
 3322  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  grep "6895[56]"
 3323  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  grep "6895[56]" | wc
 3324  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  grep "6895[56]" 
 3325  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  grep "6895[56]"  | sort | uniq -c
 3328  ./mjm_string_seq.out -cmd "stream-fasta - yyy" -cmd quit  | grep -o -n "W..W..W..C\|>.*"| grep -v -B 1  ">" | highlight -red W..W..W..C |  grep "6895[56]"  | sort | uniq -c
 3332  history | grep stream-fasta >> README.txt 


 3269  ./mjm_string_seq.out -cmd "stream-script rag xxx \"grep -n W..W \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged ragg rag" -cmd quit
 3284  ./mjm_string_seq.out -cmd "stream-script rag xxx \"grep -n W..W \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged ragg rag" -cmd quit
 3288  ./mjm_string_seq.out -cmd "stream-script rag xxx \"grep -n W..W..W..C \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged ragg rag" -cmd quit
 3298  ./mjm_string_seq.out -cmd "stream-script rag yyy \"grep -n W..W..W..C \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged raggu rag" -cmd quit
 3329  #./mjm_string_seq.out -cmd "stream-script rag yyy \"grep -n W..W..W..C \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged raggu rag" -cmd quit
 3334  history | grep stream-script >> README.txt 

# this seems to work but should be "off by one" including the ">" but that is only on the stream-fasta version... 
 3339  ./run_string_seq_new -compile
 3340  ./mjm_string_seq.out -cmd "stream-script rag yyy \"grep -o -n W..W..W..C \\| sed -e 's/:/ /g' \\|awk '{print \$1}' \"" -cmd "write-ragged raggu rag" -cmd quit
 3341  cat raggu | grep -v "|0|"

