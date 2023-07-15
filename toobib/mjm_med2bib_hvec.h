 AddHandlerPair("handlecitmgr",&Myt::handlecitmgr,"handlecitmgr"); 
 AddHandlerPair("handlespring",&Myt::handlespring,"handlespring"); 
 AddHandlerPair("handledowncit",&Myt::handledowncit,"handledowncit"); 
 AddHandlerPair("handlehighwire",&Myt::handlehighwire,"handlehighwire"); 
 AddHandlerPair("handlebiblink",&Myt::handlebiblink,"handlebiblink"); 
 AddHandlerPair("handledoilink",&Myt::handledoilink,"handledoilink"); 
// this tries to scrape the pdf first before strings failing on
// out of roder text, use pdftotext -layour first. 
 AddHandlerPair("handlepdftext",&Myt::handlepdftext,"handlepdftext"); 
// see comments above on why above was moved over this.  
AddHandlerPair("handledoi",&Myt::handledoi,"handledoi"); 
 AddHandlerPair("handleisbib",&Myt::handleisbib,"handleisbib"); 
 AddHandlerPair("handlehasbib",&Myt::handlehasbib,"handlehasbib"); 
 AddHandlerPair("handlebibtexlink",&Myt::handlebibtexlink,"handlehibtexlink"); 
 AddHandlerPair("handlebibxml",&Myt::handlebibxml,"handlebibxml"); 
 AddHandlerPair("handledoixml",&Myt::handledoixml,"handledoixml"); 
 AddHandlerPair("handlepdfexif",&Myt::handlepdfexif,"handlepdfexif"); 
 //AddHandlerPair("handlepdf",&Myt::handlepdf,"handlepdf"); 
 AddHandlerPair("handlepdf",&Myt::handlepdf,"handlepdf"); 
 AddHandlerPair("handlediss",&Myt::handlediss,"handlediss"); 
 AddHandlerPair("handlexmlformats",&Myt::handlexmlformats,"handlexmlformats"); 
 AddHandlerPair("handlesilverchair",&Myt::handlesilverchair,"handlesilverchair"); 
 AddHandlerPair("handlepdfstrings",&Myt::handlepdfstrings,"handlepdfstrings"); 
 AddHandlerPair("handlemutate",&Myt::handlemutate,"handlemutate"); 
// this is a real kluge even by the hacking standards here
// and is now failing on cleveer html but should be obsoleted by adhochtml
// AddHandlerPair("handlegsmeta",&Myt::handlegsmeta,"handlegsmeta"); 
 //AddHandlerPair("handleldjson",&Myt::handleldjson2,"handleldjson"); 
 AddHandlerPair("handleldjson2",&Myt::handleldjson2,"handleldjson2"); 
 AddHandlerPair("handleadhochtml",&Myt::handleadhochtml,"handleadhochtml"); 

// this is a real kluge even by the hacking standards here
// and is now failing on cleveer html but should be obsoleted by adhochtml
 AddHandlerPair("handlegsmeta",&Myt::handlegsmeta,"handlegsmeta"); 

 AddHandlerPair("handlebomtex",&Myt::handlebomtex,"handlebomtex"); 

 AddHandlerPair("handleheadless",&Myt::handleheadless,"handleheadless"); 
