 AddPair("guessdoi",&Myt::guessdoi,"^10\\."); 
 AddPair("guessdoilink",&Myt::guessdoi,"http.*doi\\.org/"); 
AddPair("guesscell",&Myt::guesscell,"=www.cell.com"); 
 AddPair("guessaip",&Myt::guessaip,"aip.org"); 
 AddPair("guesshind",&Myt::guesshind,"hindawi.com"); 
 AddPair("guesspmc",&Myt::guesspmc,"ncbi.nlm.nih.*PMC[0-9]"); 
 AddPair("guesspmceuro",&Myt::guesspmceuro,"europepmc.org.*PMC/[0-9]"); 
 AddPair("guessepmc",&Myt::guessepmc,"europepmc.org"); 
 AddPair("guessasm",&Myt::guessasm,"asm.org"); 
 AddPair("guessnature",&Myt::guessnature,"nature.com"); 
 AddPair("guessjes",&Myt::guessjes,"jes.ecsdl.org"); 
 AddPair("guessacs",&Myt::guessacs,"pubs.acs.org"); 
 AddPair("guessrsc",&Myt::guessrsc,"pubs.rsc.org"); 
 AddPair("guessscidirect",&Myt::guessscidirect,"sciencedirect.com"); 
 AddPair("guesspubmed",&Myt::guesspubmed,"gov/pubmed/"); 
 AddPair("guesspubmed",&Myt::guesspubmed,"pubmed.ncbi.nlm.nih.gov/"); 
 AddPair("guessmdpi",&Myt::guessmdpi,"www.mdpi.com/"); 
 AddPair("guessspiedigital",&Myt::guessspiedigital,"spiedigitallibrary.org"); 
 AddPair("guessoxford",&Myt::guessoxford,"femsle.oxfordjournals.org/"); 
 AddPair("guessoxfordpcp",&Myt::guessoxfordpcp,"pcp.oxfordjournals.org"); 
 AddPair("guesstandf",&Myt::guesstandf,"tandfonline.com"); 
 AddPair("guessspringer",&Myt::guessspringer,"springer.com/"); 
 AddPair("guessasce",&Myt::guessasce,"ascelibrary.org"); 
 AddPair("guesswiley",&Myt::guesswiley,"wiley.com"); 
 AddPair("guesssci",&Myt::guesssci,"sciencemag.org"); 
 AddPair("guessrupress",&Myt::guessrupress,"rupress.org/"); 
 AddPair("guesscopyclip",&Myt::guesscopyclip,"@[a-zA-Z]*\\{"); 
 AddPair("guessfuturesci",&Myt::guessfuturesci,"future-science.com"); 
 AddPair("guesspnas",&Myt::guesspnas,"www.pnas.org"); 
 AddPair("guessscielo",&Myt::guessscielo,"www.scielo.br"); 
 AddPair("guessjama",&Myt::guessjama,"jamanetwork.com"); 
 AddPair("guesskarger",&Myt::guesskarger,"karger.com"); 
 AddPair("guesssage",&Myt::guesssage,"sagepub.com"); 
 AddPair("guessingenta",&Myt::guessingenta,"ingentaconnect.com"); 
 AddPair("guessciteseer",&Myt::guessciteseer,"citeseerx.ist.psu.edu"); 
 AddPair("guessosa",&Myt::guessosa,"osapublishing.org"); 
 AddPair("guessscientificnet",&Myt::guessscientificnet,"www.scientific.net"); 
 AddPair("guessiop",&Myt::guessiop,"iopscience.iop.org"); 
// needs to get other hosts or no www ok to fail on sites
// ending in jci.org need to check "." lol
// TODO escape the literal periods doh 
 //AddPair("guessjci",&Myt::guessjci,"www.jci.org"); 
 AddPair("guessjci",&Myt::guessjci,"[.]jci.org"); 
 AddPair("guessgoopat",&Myt::guessgoopat,"google.com/patents"); 
 AddPair("guessieee",&Myt::guessieee,"ieeexplore.ieee.org"); 
 AddPair("guessbegel",&Myt::guessbegel,"dl.begellhouse.com"); 
 AddPair("guessaps",&Myt::guessaps,"journals.aps.org"); 
 AddPair("guessfuture",&Myt::guessfuture,"futuremedicine.com"); 
 AddPair("guessjstage",&Myt::guessjstage,"jstage.jst.go.jp"); 
 AddPair("guessannrev",&Myt::guessannrev,"annualreviews.org"); 
 AddPair("guesstwente",&Myt::guesstwente,"utwente.nl"); 
 AddPair("guessplos",&Myt::guessplos,"plos.org"); 
 AddPair("guessarxivfour",&Myt::guessarxivfour,"arXiv:"); 
 AddPair("guessarxivfour",&Myt::guessarxivfour,"arxiv.org"); 
 AddPair("guessarxivthreesom",&Myt::guessarxivthreesom,"arxiv.org"); 
 AddPair("guessarxivthree",&Myt::guessarxivthree,"arxiv.org"); 
 AddPair("guessoup",&Myt::guessoup,"academic.oup.com"); 
 AddPair("guessnejm",&Myt::guessnejm,"nejm.org"); 
 AddPair("guessphysiology",&Myt::guessphysiology,"physiology.org"); 
 AddPair("guessbiomedcentralcom",&Myt::guessbiomedcentralcom,"biomedcentral.com"); 
 AddPair("guessnapedu",&Myt::guessnapedu,"nap.edu"); 
 AddPair("guessembo",&Myt::guessembo,"embopress.org"); 
 AddPair("guesscambridge",&Myt::guesscambridge,"cambridge.org"); 
 AddPair("guessabo",&Myt::guessabo,"ojs.abo.fi"); 
 AddPair("guessmcmaster",&Myt::guessmcmaster,"mcmaster.ca"); 
 AddPair("guessresearchgate2",&Myt::guessresearchgate2,"researchgate.net"); 
 AddPair("guessagrisfao",&Myt::guessagrisfao,"agris.fao.org"); 
 AddPair("guessrugnlportal",&Myt::guessrugnlportal,"rug.nl/research/portal"); 
 AddPair("guesslibert",&Myt::guesslibert,"liebertpub.com"); 
 AddPair("guessahajournal",&Myt::guessahajournal,"ahajournals.org"); 
 AddPair("guessatsjournal",&Myt::guessatsjournal,"atsjournals.org"); 
 AddPair("guessarchpath",&Myt::guessarchpath,"archivesofpathology.org"); 
 AddPair("guessjbc",&Myt::guessjbc,".jbc.org"); 
 AddPair("guesscern",&Myt::guesscern,"cds.cern.ch"); 
 AddPair("guessmicrores",&Myt::guessmicrores,".microbiologyresearch.org"); 
 AddPair("guessgenetics",&Myt::guessgenetics,"www.genetics.org"); 
 AddPair("guesssemantic",&Myt::guesssemantic,"semanticscholar.org"); 
 AddPair("guessdiabetesjournal",&Myt::guessdiabetesjournal,"diabetesjournals.org"); 
 AddPair("guessfasebj",&Myt::guessfasebj,"fasebj.org"); 
 AddPair("guesscab",&Myt::guesscab,"cabdirect.org"); 
 AddPair("guessblood",&Myt::guessblood,"bloodjournal.org"); 
 AddPair("guesshathi",&Myt::guesshathi,"hathitrust.org"); 
 AddPair("guesscelt",&Myt::guesscelt,"research.ucc.ie"); 
 AddPair("guesspubag",&Myt::guesspubag,"pubag.nal.usda.gov"); 
 AddPair("guessjstor",&Myt::guessjstor,"jstor.org"); 
 AddPair("guessauthorea",&Myt::guessauthorea,"authorea.com"); 

 AddPair("guessdigitalcommons",&Myt::guessdigitalcommons,"digitalcommons\\..*\\.edu"); 
 AddPair("guessuridc",&Myt::guessuridc,"digitalcommons.uri.edu"); 
 AddPair("guesstracetenedu",&Myt::guesstracetenedu,"trace.tennessee.edu"); 
 AddPair("guesskidint",&Myt::guesskidint,"kidney-international.org"); 
 AddPair("guesselife",&Myt::guesselife,"/elifesciences.org"); 
 AddPair("guessjlr",&Myt::guessjlr,"www.jlr.org"); 
 AddPair("guessfcklww",&Myt::guessfcklww,"\\.lww\\.com"); 
 AddPair("guessthoracicsurgery",&Myt::guessthoracicsurgery,"annalsthoracicsurgery.org"); 
 AddPair("guesscmajca",&Myt::guesscmajca,"cmaj.ca"); 
// moved tohandler for all 
 //AddPair("handlemutate",&Myt::handlemutate,"[/.]pdf$"); 
 AddPair("guessdora",&Myt::guessdora,"dora.lib4ri.ch.*datastream"); 
 AddPair("guessfrontiersin",&Myt::guessfrontiersin,"frontiersin.org"); 
// https://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=9305&context=rtd
// will use recursion so get specifici
 AddPair("guessiastatertd",&Myt::guessiastatertd,"lib.dr.iastate.edu.*viewcontent.*context=rtd"); 
// http://cogprints.org/8948/1/2012-4-19.pdf
//http://cogprints.org/cgi/export/eprint/9825/BibTeX/cogprints-eprint-9825.bib
 AddPair("guesscogprints",&Myt::guesscogprints,"/cogprints.org/"); 
// https://heart.bmj.com/content/52/1/82 
AddPair("guessbmj",&Myt::guessbmj,"\\.bmj\\.com/"); 
AddPair("guessnber",&Myt::guessnber,"\\.nber\\.org/"); 
AddPair("guessgooglebooks",&Myt::guessgooglebooks,"[./]books\\.google\\.com/"); 
AddPair("guessgooglebooks",&Myt::guessgooglebooks,"\\.google\\.com/books/"); 
AddPair("guessworldcat",&Myt::guessworldcat,"[\\/.]worldcat\\.org/"); 

// https://ufdc.ufl.edu/AA00004956/00001
//http://ufdcimages.uflib.ufl.edu/AA/00/00/49/56/00001/AA00004956.pdf
AddPair("guessufl",&Myt::guessufl,"\\.uflib\\.ufl\\.edu/"); 
AddPair("guessufl",&Myt::guessufl,"ufdc\\.ufl\\.edu/"); 

// https://www.ijvsbt.org/index.php/journal/article/view/1386
// https://www.ijvsbt.org/index.php/journal/article/download/1386/1058
AddPair("guessijvsbt",&Myt::guessijvsbt,"\\.ijvsbt\\.org/"); 

AddPair("guessohiolink",&Myt::guessohiolink,"\\.ohiolink\\.edu/"); 
// https://discovery.ucl.ac.uk/id/eprint/1538561/1/Laycock_Thesis%20Joanne%20Laycock%20(No%20copyright).pdf
AddPair("guessdiscoverac",&Myt::guessdiscoverac,"//discover.*\\.ac\\.uk/"); 
// https://core.ac.uk/reader/76958504?utm_source=linkout
AddPair("guesscoreac",&Myt::guesscoreac,"//core\\.ac\\.uk/"); 

// https://en.wikipedia.org/wiki/Wikipedia:Citing_Wikipedia
AddPair("guesswikipedia",&Myt::guesswikipedia,"wikipedia\\.org/"); 

AddPair("guesspubchem",&Myt::guesspubchem,"pubchem\\.ncbi\\.nlm\\.nih\\.gov/"); 

AddPair("guessclinicaltrialsgov",&Myt::guessclinicaltrialsgov,"[/.]clinicaltrials\\.gov/"); 
AddPair("guesslinkedin",&Myt::guesslinkedin,"\\.linkedin\\.com/"); 
AddPair("guessnpr",&Myt::guessnpr,"www\\.npr\\.org/"); 
AddPair("guessnpr",&Myt::guessnpr,"www\\.nbcnews\\.com/"); 
AddPair("guessnpr",&Myt::guessnpr,"www\\..*news.*\\.com/"); 

AddPair("guessdiva",&Myt::guessdiva,"\\.diva-portal\\.org/"); 
// huge filt but the doi in link works, 
// https://cdnsciencepub.com/doi/abs/10.1139/o59-099
AddPair("guesscdnscience",&Myt::guesscdnscience,"cdnsciencepub\\.com/"); 

// https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8036
// https://bonndoc.ulb.uni-bonn.de/xmlui/bitstream/handle/20.500.11811/8036/5508.pdf;jsessionid=8EE2D4C371DD10838C62271A1D5CEC0E?sequence=1
AddPair("guessbonndoc",&Myt::guessbonndoc,"/bonndoc\\..*uni-bonn\\.de/"); 



// this has a DOI in non-standard form, DOI: ... ZZ
// https://webofproceedings.org/proceedings_series/ALSMB/MedSBE%202020/ZZEM740.pdf
// https://webofproceedings.org/proceedings_series/article/artId/18734.html
// DOI: 10.25236/medsbe.2020.042
// the fcking doi is not registered with crossref fck 
AddPair("guesswegofproc",&Myt::guesswebofproc,"\\.webproceedings\\.org/"); 
// https://www.jocms.org/index.php/jcms/article/view/1100/499
AddPair("guessjocms",&Myt::guessjocms,"\\.jocms\\.org/"); 

// 2022-10-29 started
// search this, 
// http://apps.who.int/iris/handle/10665/43407?search-result=true&query=0300-0923&scope=%2F&rpp=10&sort_by=score&order=desc
// using scraped ISSN from here, 
// https://apps.who.int/iris/rest/bitstreams/1425664/retrieve
AddPair("guesswhoiris",&Myt::guesswhoiris,"\\.who\\.int/iris/"); 

// https://iris.univpm.it/retrieve/handle/11566/263460/77778/Tesi_Cirilli.pdf
// https://iris.univpm.it/retrieve/handle/11566/263460/
// http://hdl.handle.net/11566/263460 
AddPair("guesshandle",&Myt::guesshandle,"/handle/"); 
AddPair("guessbitstream",&Myt::guessbitstream,"/bitstream/"); 

// 2011  #./chromate.out "mode 0" launch "wget xxx.pdf https://authors.library.caltech.edu/32187/1/BROsiamjssc85.pdf"
// 2012  # https://authors.library.caltech.edu/32187/
// 2023-02-25
AddPair("guesscaltech",&Myt::guesscaltech,"/authors\\.library\\.caltech"); 

//https://apps.dtic.mil/sti/pdfs/ADA121307.pdf
//https://apps.dtic.mil/sti/citations/ADA121307
// 2023-02-25
AddPair("guessdticmil",&Myt::guessdticmil,"\\.dtic\\.mil/"); 



// https://api.crossref.org/works/10.3987/com-06-s(o)36/transform/application/x-bibtex
AddPair("handledoilink10",&Myt::handledoilink,"10\\.[0-9][0-9]"); 
AddPair("handledoilinkdoi",&Myt::handledoilink,"doi"); 
