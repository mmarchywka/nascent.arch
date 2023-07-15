#ifndef MJM_BIBTEX_DOMAINS_H__
#define MJM_BIBTEX_DOMAINS_H__



/*
//typedef IdxTy (*Guess)(Myt * mytp, TempProg & p, const IdxTy flags); 
m_gdv.push_back(GuessDomain("",&Myt::NullGuess));
LoadGdV();


}
#include "mjm_bibtex_domains.h"
//IdxTy NullGuess(Myt * mytp, TempProg & p, const IdxTy flags )
HANDLER(NullGuess) { return 0; }




*/

void LoadGdV()
{

//m_gdv.push_back(GuessDomain("",&Myt::NullGuess));
//m_gdv.push_back(GuessDomain("",&Myt::NullGuess));


//}

//HANDLER(NullGuess) { return 0; }



m_gdv.push_back(GuessDomain("www.cell.com",&Myt::guesscell));
m_gdv.push_back(GuessDomain("aip.org",&Myt::guessaip));
m_gdv.push_back(GuessDomain("hindawi.com",&Myt::guesshind));
m_gdv.push_back(GuessDomain("ncbi.nlm.nih.*PMC[0-9]",&Myt::guesspmc));
m_gdv.push_back(GuessDomain("asm.org",&Myt::guessasm));
m_gdv.push_back(GuessDomain("nature.com",&Myt::guessnature));
m_gdv.push_back(GuessDomain("jes.ecsdl.org",&Myt::guessjes));
m_gdv.push_back(GuessDomain("pubs.acs.org",&Myt::guessacs));
m_gdv.push_back(GuessDomain("pubs.rsc.org",&Myt::guessrsc));
m_gdv.push_back(GuessDomain("sciencedirect.com",&Myt::guessscidirect));
m_gdv.push_back(GuessDomain("gov/pubmed/",&Myt::guesspubmed));
m_gdv.push_back(GuessDomain("www.mdpi.com/",&Myt::guessmdpi));
m_gdv.push_back(GuessDomain("spiedigitallibrary.org",&Myt::guessspiedigital));
m_gdv.push_back(GuessDomain("femsle.oxfordjournals.org/",&Myt::guessoxford));
m_gdv.push_back(GuessDomain("pcp.oxfordjournals.org",&Myt::guessoxfordpcp));
m_gdv.push_back(GuessDomain("tandfonline.com",&Myt::guesstandf));
m_gdv.push_back(GuessDomain("springer.com/",&Myt::guessspringer));
m_gdv.push_back(GuessDomain("ascelibrary.org",&Myt::guessasce));
m_gdv.push_back(GuessDomain("wiley.com",&Myt::guesswiley));
m_gdv.push_back(GuessDomain("sciencemag.org",&Myt::guesssci));
m_gdv.push_back(GuessDomain("rupress.org/",&Myt::guessrupress));
m_gdv.push_back(GuessDomain("@[a-zA-Z]*{",&Myt::guesscopyclip));
m_gdv.push_back(GuessDomain("future-science.com",&Myt::guessfuturesci));
m_gdv.push_back(GuessDomain("www.pnas.org",&Myt::guesspnas));
m_gdv.push_back(GuessDomain("www.scielo.br",&Myt::guessscielo));
m_gdv.push_back(GuessDomain("jamanetwork.com",&Myt::guessjama));
m_gdv.push_back(GuessDomain("karger.com",&Myt::guesskarger));
m_gdv.push_back(GuessDomain("sagepub.com",&Myt::guesssage));
m_gdv.push_back(GuessDomain("ingentaconnect.com",&Myt::guessingenta));
m_gdv.push_back(GuessDomain("citeseerx.ist.psu.edu",&Myt::guessciteseer));
m_gdv.push_back(GuessDomain("osapublishing.org",&Myt::guessosa));
m_gdv.push_back(GuessDomain("www.scientific.net",&Myt::guessscientificnet));
m_gdv.push_back(GuessDomain("iopscience.iop.org",&Myt::guessiop));
m_gdv.push_back(GuessDomain("www.jci.org",&Myt::guessjci));
m_gdv.push_back(GuessDomain("google.com/patents",&Myt::guessgoopat));
m_gdv.push_back(GuessDomain("ieeexplore.ieee.org",&Myt::guessieee));
m_gdv.push_back(GuessDomain("dl.begellhouse.com",&Myt::guessbegel));
m_gdv.push_back(GuessDomain("journals.aps.org",&Myt::guessaps));
m_gdv.push_back(GuessDomain("futuremedicine.com",&Myt::guessfuture));
m_gdv.push_back(GuessDomain("jstage.jst.go.jp",&Myt::guessjstage));
m_gdv.push_back(GuessDomain("annualreviews.org",&Myt::guessannrev));
m_gdv.push_back(GuessDomain("utwente.nl",&Myt::guesstwente));
m_gdv.push_back(GuessDomain("plos.org",&Myt::guessplos));
m_gdv.push_back(GuessDomain("arxiv.org",&Myt::guessarxivthree));
m_gdv.push_back(GuessDomain("academic.oup.com",&Myt::guessoup));
m_gdv.push_back(GuessDomain("nejm.org",&Myt::guessnejm));
m_gdv.push_back(GuessDomain("physiology.org",&Myt::guessphysiology));
m_gdv.push_back(GuessDomain("biomedcentral.com",&Myt::guessbiomedcentralcom));
m_gdv.push_back(GuessDomain("nap.edu",&Myt::guessnapedu));
m_gdv.push_back(GuessDomain("embopress.org",&Myt::guessembo));
m_gdv.push_back(GuessDomain("cambridge.org",&Myt::guesscambridge));
m_gdv.push_back(GuessDomain("ojs.abo.fi",&Myt::guessabo));
m_gdv.push_back(GuessDomain("mcmaster.ca",&Myt::guessmcmaster));
m_gdv.push_back(GuessDomain("researchgate.net",&Myt::guessresearchgate2));
m_gdv.push_back(GuessDomain("agris.fao.org",&Myt::guessagrisfao));
m_gdv.push_back(GuessDomain("rug.nl/research/portal",&Myt::guessrugnlportal));
m_gdv.push_back(GuessDomain("liebertpub.com",&Myt::guesslibert));
m_gdv.push_back(GuessDomain("ahajournals.org",&Myt::guessahajournal));
m_gdv.push_back(GuessDomain("atsjournals.org",&Myt::guessatsjournal));
m_gdv.push_back(GuessDomain("archivesofpathology.org",&Myt::guessarchpath));
m_gdv.push_back(GuessDomain(".jbc.org",&Myt::guessjbc));
m_gdv.push_back(GuessDomain(".microbiologyresearch.org",&Myt::guessmicrores));
m_gdv.push_back(GuessDomain("www.genetics.org",&Myt::guessgenetics));
m_gdv.push_back(GuessDomain("semanticscholar.org",&Myt::guesssemantic));
m_gdv.push_back(GuessDomain("diabetesjournals.org",&Myt::guessdiabetesjournal));
m_gdv.push_back(GuessDomain("fasebj.org",&Myt::guessfasebj));
m_gdv.push_back(GuessDomain("cabdirect.org",&Myt::guesscab));
m_gdv.push_back(GuessDomain("bloodjournal.org",&Myt::guessblood));
m_gdv.push_back(GuessDomain("hathitrust.org",&Myt::guesshathi));
m_gdv.push_back(GuessDomain("research.ucc.ie",&Myt::guesscelt));
m_gdv.push_back(GuessDomain("pubag.nal.usda.gov",&Myt::guesspubag));
m_gdv.push_back(GuessDomain("jstor.org",&Myt::guessjstor));
m_gdv.push_back(GuessDomain("digitalcommons.uri.edu",&Myt::guessuridc));
m_gdv.push_back(GuessDomain("trace.tennessee.edu",&Myt::guesstracetenedu));
m_gdv.push_back(GuessDomain("kidney-international.org",&Myt::guesskidint));
}

 HANDLER(guesscell) {return 0; }
 HANDLER(guessaip) {return 0; }
 HANDLER(guesshind) {return 0; }
 HANDLER(guesspmc) {return 0; }
 HANDLER(guessasm) {return 0; }
 HANDLER(guessnature) {return 0; }
 HANDLER(guessjes) {return 0; }
 HANDLER(guessacs) {return 0; }
 HANDLER(guessrsc) {return 0; }
 HANDLER(guessscidirect) {return 0; }
 HANDLER(guesspubmed) {return 0; }
 HANDLER(guessmdpi) {return 0; }
 HANDLER(guessspiedigital) {return 0; }
 HANDLER(guessoxford) {return 0; }
 HANDLER(guessoxfordpcp) {return 0; }
 HANDLER(guesstandf) {return 0; }
 HANDLER(guessspringer) {return 0; }
 HANDLER(guessasce) {return 0; }
 HANDLER(guesswiley) {return 0; }
 HANDLER(guesssci) {return 0; }
 HANDLER(guessrupress) {return 0; }
 HANDLER(guesscopyclip) {return 0; }
 HANDLER(guessfuturesci) {return 0; }
 HANDLER(guesspnas) {return 0; }
 HANDLER(guessscielo) {return 0; }
 HANDLER(guessjama) {return 0; }
 HANDLER(guesskarger) {return 0; }
 HANDLER(guesssage) {return 0; }
 HANDLER(guessingenta) {return 0; }
 HANDLER(guessciteseer) {return 0; }
 HANDLER(guessosa) {return 0; }
 HANDLER(guessscientificnet) {return 0; }
 HANDLER(guessiop) {return 0; }
 HANDLER(guessjci) {return 0; }
 HANDLER(guessgoopat) {return 0; }
 HANDLER(guessieee) {return 0; }
 HANDLER(guessbegel) {return 0; }
 HANDLER(guessaps) {return 0; }
 HANDLER(guessfuture) {return 0; }
 HANDLER(guessjstage) {return 0; }
 HANDLER(guessannrev) {return 0; }
 HANDLER(guesstwente) {return 0; }
 HANDLER(guessplos) {return 0; }
 HANDLER(guessarxivthree) {return 0; }
 HANDLER(guessoup) {return 0; }
 HANDLER(guessnejm) {return 0; }
 HANDLER(guessphysiology) {return 0; }
 HANDLER(guessbiomedcentralcom) {return 0; }
 HANDLER(guessnapedu) {return 0; }
 HANDLER(guessembo) {return 0; }
 HANDLER(guesscambridge) {return 0; }
 HANDLER(guessabo) {return 0; }
 HANDLER(guessmcmaster) {return 0; }
 HANDLER(guessresearchgate2) {return 0; }
 HANDLER(guessagrisfao) {return 0; }
 HANDLER(guessrugnlportal) {return 0; }
 HANDLER(guesslibert) {return 0; }
 HANDLER(guessahajournal) {return 0; }
 HANDLER(guessatsjournal) {return 0; }
 HANDLER(guessarchpath) {return 0; }
 HANDLER(guessjbc) {return 0; }
 HANDLER(guessmicrores) {return 0; }
 HANDLER(guessgenetics) {return 0; }
 HANDLER(guesssemantic) {return 0; }
 HANDLER(guessdiabetesjournal) {return 0; }
 HANDLER(guessfasebj) {return 0; }
 HANDLER(guesscab) {return 0; }
 HANDLER(guessblood) {return 0; }
 HANDLER(guesshathi) {return 0; }
 HANDLER(guesscelt) {return 0; }
 HANDLER(guesspubag) {return 0; }
 HANDLER(guessjstor) {return 0; }
 HANDLER(guessuridc) {return 0; }
 HANDLER(guesstracetenedu) {return 0; }
 HANDLER(guesskidint) {return 0; }
#endif

