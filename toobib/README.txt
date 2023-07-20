
TooBib by Mike Marchywka. See toobib.h for more complete 
information. This is being made public mostly for information
and is intended for mature developer audiences only.

The name "Toobib"  unforuntately has been used by other applications 
but this software is a collection of hacks to allow downloading
at least one bibtex entry for most webpages or pdf files that
point to a particular work. This may be pubmed, various publichser
site pages giving abstracts or full texts of  scientific works, 
pdf files with sufficient information for deduction.
It is perhaps similar to Zotero but appears to be more comprehensive
in its search finding additional information not available on
 from  https://zbib.org/  . Users are welcome to compare to products
such as this and report results or commeents back to marchywka@hotmail.com.

Toobib is only designed as a command line utility although its
console interface could be modified to function in server or demon
mode. As it does its best to find any or all bibtex entries the
quality of the results can not be guaranteed but as sites change
canned products like Zotero may fail or miss things like retraction
notices that are not always provided in "nice" places like crossref.

Toobib is amazingly simple and currently is build from a single 
command line to compile and link. It is entirely implemented
as c++ headers although a lot of work is delegated by std::system
calls. Some of this would be better integrated if there is sufficnet
interest. 

Toobib can use wget for web access but due to site issues generally
uses headless chrome and nodejs websockets as an interface. 

Among the issues include unicode and various field clean ups
although it does pretty well on the most common sources I use
in the English medical and scientific literature. 
The remote website initial gusses come from regex match
to url. Sites change and quick updates with dynamically loaded
shared objects may make some sense. However, as the build is
simple right now rebuild is just as easy. If the website's default
dedicated hnaler fails, the common handlers will try to find something
sane and usually end up working. While this can be slower, 
it is perfactly fine for most situations. 

A lot of information is hardcoded and a lot of temporary
files are generated in /tmp and the current directory.
Toobib deletes most temp files as they accumulate after a while.

For access to medline, Toobib uses my ancient "eutilsnew" which
I am still cleaning up for publication. 

requires
headless chrome
dos2unix
lynx
mjsonu
pdftotext
exiftool

eutilsnew 

