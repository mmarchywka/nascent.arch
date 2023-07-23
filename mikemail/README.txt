
can use the repo libetpan although the dev version may work
too. 
2023-05-30 new focal insteall had to copy glib-2.0 and can't
find libiculx wtf complains neomutt needs to be recomliated
for readline.  
 2001  apt-cache search sasl2
 2002  sudo apt-get install libsasl2-dev
 2003  apt-cache search tidy
 2004  apt-cache search libtidy
 2005  sudo apt-get install libtidy-dev
 2006  apt-cache search libctemplate
 2007  sudo apt-get install libctemplate-dev
 2008  apt-cache search libuuid
 2009  sudo apt-get install uuid-dev
 2010  apt-cache search libiculx
 2011  apt-cache search iculx
 2012  apt-cache search glib-2.0
 2013  find /usr/lib | grep glib
 2014  find /usr/lib | grep iculx
 2015  find /lib | grep iculx
 2016  find /lib | grep glib
 2017  ls -al /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0
 2018  history


2023-01-04 failed to get images from mms.att message on Andy from
Jan 3 2023. neomutt can display them but don't show in browser either.

2022-07-01 problems on non-inbox download. Reducing thread count ( last param)
from 25 to 2 seems to work not sure of problem... 
rm oldcrap
rm mboxx
rm attach/*

 ./mikemail.out -cmd "read-server foo hotmail.tex" -cmd "verbosity 10000" -cmd "download-folder foo church oldcrap mboxx 2"


 2148  ./mikemail.out -cmd "read-server foo hotmailm.tex" -cmd "start-session foo INBOX" -cmd "list-folders foo la" -cmd list

