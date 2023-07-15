#ifndef MJM_XML_LUIT_H__
#define MJM_XML_LUIT_H__
void xml_lut(int x, const char ** c){
if (x==1 ){ *c="XML_ELEMENT_NODE"; return ; }
if (x==2 ){ *c="XML_ATTRIBUTE_NODE"; return ; }
if (x==3 ){ *c="XML_TEXT_NODE"; return ; }
if (x==4 ){ *c="XML_CDATA_SECTION_NODE"; return ; }
if (x==5 ){ *c="XML_ENTITY_REF_NODE"; return ; }
if (x==6 ){ *c="XML_ENTITY_NODE"; return ; }
if (x==7 ){ *c="XML_PI_NODE"; return ; }
if (x==8 ){ *c="XML_COMMENT_NODE"; return ; }
if (x==9 ){ *c="XML_DOCUMENT_NODE"; return ; }
if (x==10 ){ *c="XML_DOCUMENT_TYPE_NODE"; return ; }
if (x==11 ){ *c="XML_DOCUMENT_FRAG_NODE"; return ; }
if (x==12 ){ *c="XML_NOTATION_NODE"; return ; }
if (x==13 ){ *c="XML_HTML_DOCUMENT_NODE"; return ; }
if (x==14 ){ *c="XML_DTD_NODE"; return ; }
if (x==15 ){ *c="XML_ELEMENT_DECL"; return ; }
if (x==16 ){ *c="XML_ATTRIBUTE_DECL"; return ; }
if (x==17 ){ *c="XML_ENTITY_DECL"; return ; }
if (x==18 ){ *c="XML_NAMESPACE_DECL"; return ; }
if (x==19 ){ *c="XML_XINCLUDE_START"; return ; }
if (x==20 ){ *c="XML_XINCLUDE_END"; return ; }
if (x==21 ){ *c="XML_DOCB_DOCUMENT_NODE"; return ; }
}
#endif

