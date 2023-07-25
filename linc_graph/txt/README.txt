
Current usage,

1) Edit the data file with daily.txt

2) This creates the updated data file used by 

 2138  viw mudaily.txt
 2139  mudaily.txt -reports 2000-01-01 2030-01-01
 2140  ls -al  /home/scripts/script_data/cases/out/dog_daily.txt
 2141  more /home/scripts/script_data/cases/out/dog_daily.txt  |  head
 2142  more /home/scripts/script_data/cases/out/dog_daily.txt  |  grep Andy

3) Generate report files, 

 2143  mudaily.txt -reports 2019-01-01 2030-01-01

4) Look for errors, 
grep -B 1 "#"   /tmp/cases_markup.txt

edit nouns if needed or cases file, 
 2164  mudaily.txt -nouns

5) then should be read tu run 

6) Most recent things, individual diet graphs and diet
tables, 

eog spicey_each/diettablex_ctbrothbs\(-\).svg 
 2092  cat /home/scripts/script_data/cases/out/dog_glob.txt  | grep Spicey | grep lysinehcl  | grep 2021-03 
 2093  rm spicey_each/*.svg
 2094  ls spicey_each/
 2095  ./backup
 2096  ./run_linc_graph -adp pin_spicey.txt 
 2097  history | grep ey2
 2098  ./run_linc_graph -dt-mo spicey2.txt 
 2099  history | grep texfr
 2100  texfrag -include xxxtable



