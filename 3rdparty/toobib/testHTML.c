// downloaded from https://dev.w3.org/XML/testHTML.c
/*
 * testHTML.c : a small tester program for HTML input.
 *
 * See Copyright for the status of this software.
 *
 * Daniel.Veillard@w3.org
 */


// mjm added other sample code use for general utility
// gcc -I/usr/include/libxml2 testHTML.c -o testHTML  -lxml2
// or with the char conversiotn, 
//g++ -std=c++11 -I/usr/include/libxml2 -I../../mjm/hlib  -x c++ testHTML.c -o testHTML  -lxml2


#ifdef WIN32
#include "win32config.h"
#undef LIBXML_DLL_IMPORT
#else
//#include "config.h"
#endif

#include <libxml/xmlversion.h>
#ifdef LIBXML_HTML_ENABLED

#include <stdio.h>
#include <string.h>
#include <stdarg.h>


#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <libxml/xmlmemory.h>
#include <libxml/HTMLparser.h>
#include <libxml/HTMLtree.h>
#include <libxml/debugXML.h>
#include <libxml/xmlerror.h>
#include "xml_lut.h" 


#define HAVE_MJM_STUFF 1
#ifdef HAVE_MJM_STUFF
#include "mjm_hierarchial_kvp.h"
#include "mjm_misc_parse.h" 
typedef std::map<void*,int> NodeMap;
NodeMap nm;
#endif


#ifdef LIBXML_DEBUG_ENABLED
static int debug = 0;
#endif
static int copy = 0;
static int sax = 0;
static int repeat = 0;
static int noout = 0;
static int push = 0;
static char *encoding = NULL;

///////////////////////////////////////////////////////////
 // from test1.c

static void print_attributes ( xmlAttributePtr p );
static void print_properties ( xmlAttr*  p, int lvl );
mjm_var_parse<mjm_generic_traits> var_parse;
typedef mjm_misc_parse<mjm_generic_traits> LutTy;

struct form_spec {
int stack=1;
int print_stacked=1;
int print_curr=0;
int remove_crlf=1;
int include_node_numbers=0;
int suppress_levels=0;
} ;

int mjm_format(char * d, const unsigned char * p, form_spec * fsp)
{
int i=0;
int j=0;
int n=0;
if (p==0) { d[0]=0; return n  ; } 
const bool check_eol=(fsp==0)||(fsp->remove_crlf!=0);
const unsigned int * lut=var_parse.lut();
while (p[i]!=0)
{

const unsigned int code= lut[p[i]];
//const unsigned int code= 0;
if (check_eol&&((code&LutTy::EOL)!=0)) { ++n; } 
else { d[j]=p[i]; ++j; } 
++i;
}
d[j]=0;
return n;
}


static void print_properties_stack (std::string pfx,  xmlAttr*  p, form_spec * fsp,  int lvl );

static void print_element_names_stack(std::string  pfx, xmlNode * a_node, form_spec * fsp, int lvl)
{
//printf("  xxx loc a\n");
//char cb[1<<24];
//char cb[1<<16];
//MM_MSG("asdf test")
char * cb= new char[1<<24];
//printf("  xxx loc b\n");
    xmlNode *cur_node = a_node;
	if (nm.find(cur_node)==nm.end()) nm[cur_node]=nm.size(); 
    for (; cur_node; cur_node = cur_node->next) 
	{
	if (nm.find(cur_node)==nm.end()) nm[cur_node]=nm.size(); 
int noden=nm[cur_node];
	cb[0]=0;
	const char * c="UNKNOWN";
	xml_lut(cur_node->type,&c);
	const unsigned char * p=cur_node->content;
	const int msz=(p==0)?10:(strlen((const char *)p)*3+3);
	char px[msz];
	int dn=mjm_format(px,p,fsp);
       //printf("_NT_: %s _lvl_ %d _name_ %s _content_ %s \n",c, lvl,cur_node->name,px);  
       if (fsp->print_curr) sprintf(cb,"%s _NT_: %s _lvl_ %d _name_ %s _content_ %s",pfx.c_str(),c, lvl,cur_node->name,px);  
int ok_to_print=1;
       //sprintf(cb,"%s %s %d %s %s",pfx.c_str(),c, lvl,cur_node->name,px);  
      if (fsp->print_stacked){
if (fsp->include_node_numbers!=0)
{
if (fsp->suppress_levels!=0)
  { 
if (px[0]||(!cur_node->children&&!cur_node->properties)) 
{ } else { ok_to_print=0; }    
//if (ok_to_print) sprintf(cb,"%s %s %d = %s",pfx.c_str(), cur_node->name,noden,px);  
if (ok_to_print) sprintf(cb,"%s%d %s = %s",pfx.c_str(), noden,cur_node->name,px);  
//else sprintf(cb,"%s %s %d ",pfx.c_str(), cur_node->name,noden);  
else sprintf(cb,"%s%d %s ",pfx.c_str(), noden,cur_node->name);  
//else ok_to_print=0;
}

  //sprintf(cb,"%s %s [ anode= %d %d  ]  %s",pfx.c_str(), cur_node->name,nm[a_node],noden,px);  
else   sprintf(cb,"%s %d %s [ %d   %d ]  %s",pfx.c_str(), lvl,cur_node->name,nm[a_node],noden,px);  
}
  else 
{
if (fsp->suppress_levels!=0)
sprintf(cb,"%s %s %s",pfx.c_str(), cur_node->name,px);  
else sprintf(cb,"%s %d %s %s",pfx.c_str(), lvl,cur_node->name,px);  
} // else numb

} // print_stacked

if (ok_to_print) 	printf("%s\n",cb);
	std::string pfx2=(fsp->stack)?std::string(cb):pfx;
       	print_properties_stack(pfx2,cur_node->properties,fsp,lvl+1); 
        print_element_names_stack(pfx2,cur_node->children,fsp, lvl+1);
	} // for
delete []cb;
}
static void print_properties_stack (std::string pfx,  xmlAttr*  p, form_spec * fsp,  int lvl )
{
xmlAttr*  c=NULL;
if (nm.find(p)==nm.end()) nm[p]=nm.size();
//char cb[1<<16];
char* cb = new  char[1<<16];
for ( c=p; c; c=c->next)
{
if (nm.find(c)==nm.end()) nm[c]=nm.size();
	cb[0]=0;
//printf("_property_ %d _name_ %s\n",lvl,c->name);
if (fsp->print_curr) sprintf(cb,"%s _property_ %d _name_ %s",pfx.c_str(),lvl,c->name);
if (fsp->print_stacked)
{
if (fsp->suppress_levels!=0)
{

if (fsp->include_node_numbers)
 sprintf(cb,"%s%d %s ",pfx.c_str(),nm[c],c->name);

else sprintf(cb,"%s %s ",pfx.c_str(),c->name);
}

else sprintf(cb,"%s %d %s",pfx.c_str(),lvl,c->name);

} // print_stacked 
// for compatibility, this is just a debug line doh 
if (fsp->include_node_numbers==0) printf("%s\n",cb);
//std::string pfx2=std::string(cb);
std::string pfx2=(fsp->stack)?std::string(cb):pfx;
if (fsp->include_node_numbers)
{
// the first one has no additional info 
//sprintf(cb," [ %d %d ] ",nm[p],nm[c]);
//moved up 
// sprintf(cb,"%d ",nm[c]);
//pfx2+=std::string(cb);
}
print_element_names_stack(pfx2,c->children,fsp, lvl+1);
}

delete [] cb;

}


static void
//print_element_names(xmlNode * a_node, int lvl=0)
print_element_names(xmlNode * a_node, int lvl)
{
    xmlNode *cur_node = NULL;

    for (cur_node = a_node; cur_node; cur_node = cur_node->next) {
	const char * c="UNKNOWN";
	xml_lut(cur_node->type,&c);
	const unsigned char * p=cur_node->content;
	const int msz=(p==0)?10:(strlen((const char *)p)*3+3);
	char px[msz];
	mjm_format(px,p,0);
//        if (cur_node->type == XML_ELEMENT_NODE) {
//            printf("node type: Element, %d  name: %s content %s \n", lvl,cur_node->name,cur_node->content); } 
//            else { printf("node type: Other, %d  name: %s content %s \n", lvl,cur_node->name,cur_node->content); }
//            printf("node type:  %s , %d  name: %s content %s \n",c, lvl,cur_node->name,cur_node->content);  
            //printf("node type:  %s , %d  name: %s content %s \n",c, lvl,cur_node->name,px);  
            printf("_NT_: %s _lvl_ %d _name_ %s _content_ %s \n",c, lvl,cur_node->name,px);  
       	print_properties(cur_node->properties,lvl+1); 

        print_element_names(cur_node->children,lvl+1);
    }
}
/*
static void print_attributes ( xmlAttributePtr p )
{
xmlAttributePtr c=NULL;
for ( c=p; c; c=c->next)
{
print("     attrbutbte  name= %s\n",c->name);
}
} // print_attributes 
*/

static void print_properties ( xmlAttr*  p, int lvl )
{
xmlAttr*  c=NULL;
for ( c=p; c; c=c->next)
{
//printf("     property  name= %s value= %s content = %s \n",c->name,c->children->name,c->children->content);
//printf("     property %d   name= %s   \n",lvl,c->name);
printf("_property_ %d _name_ %s\n",lvl,c->name);
print_element_names(c->children,lvl+1);
}

} // print_attributes 


int mjm_print=0;
int mjm_flags=0;
char * mjm_sflags=0; 
    //htmlDocPtr doc = NULL, tmp;
void print_doc(htmlDocPtr doc, form_spec * fsp )
{
  xmlNode *root_element = xmlDocGetRootElement(doc);
 //print_element_names(root_element,0); 
std::string pfx;

//printf(" xxx loc c  %d %d  %d  \n",root_element,fsp,pfx.length());
pfx="";
 print_element_names_stack(pfx,root_element,fsp,0); 
//printf("  xxx loc d\n");
}

/////////////////////////////////////////////////////
#if 1
static void mjm_xmlStructuredErrorFunc(void * ud, xmlErrorPtr err)
{ 

//MM_ERR( MMPR3(err->file,err->line, err->message)) 

}
static void mjm_xmlGenericErrorFunc(void * ctx, const char * msg, ...) //  ... ...  )
{ 

MM_ERR(MMPR(msg)) 


}


void (*foo)(void*, const char*, ...);

void quiet_mjm_xml()
{ 
xmlSetGenericErrorFunc(0,&mjm_xmlGenericErrorFunc);
foo=mjm_xmlGenericErrorFunc;
initGenericErrorDefaultFunc(&foo);
xmlSetStructuredErrorFunc(0,&mjm_xmlStructuredErrorFunc);
}


#endif

//////////////////////////////////////////////////////////////




/*


struct _xmlNode {
    void *	_private	: application data
    xmlElementType	type	: type number, must be second !
    const xmlChar *	name	: the name of the node, or the entity
    struct _xmlNode *	children	: parent->childs link
    struct _xmlNode *	last	: last child link
    struct _xmlNode *	parent	: child->parent link
    struct _xmlNode *	next	: next sibling link
    struct _xmlNode *	prev	: previous sibling link
    struct _xmlDoc *	doc	: the containing document End of common p
    xmlNs *	ns	: pointer to the associated namespace
    xmlChar *	content	: the content
    struct _xmlAttr *	properties	: properties list
    xmlNs *	nsDef	: namespace definitions on this node
    void *	psvi	: for type/PSVI informations
    unsigned short	line	: line number
    unsigned short	extra	: extra data for XPath/XSLT
}
struct _xmlAttribute {
    void *	_private	: application data
    xmlElementType	type	: XML_ATTRIBUTE_DECL, must be second !
    const xmlChar *	name	: Attribute name
    struct _xmlNode *	children	: NULL
    struct _xmlNode *	last	: NULL
    struct _xmlDtd *	parent	: -> DTD
    struct _xmlNode *	next	: next sibling link
    struct _xmlNode *	prev	: previous sibling link
    struct _xmlDoc *	doc	: the containing document
    struct _xmlAttribute *	nexth	: next in hash table
    xmlAttributeType	atype	: The attribute type
    xmlAttributeDefault	def	: the default
    const xmlChar *	defaultValue	: or the default value
    xmlEnumerationPtr	tree	: or the enumeration tree if any
    const xmlChar *	prefix	: the namespace prefix if any
    const xmlChar *	elem	: Element holding the attribute
}
struct _xmlAttr {
    void *	_private	: application data
    xmlElementType	type	: XML_ATTRIBUTE_NODE, must be second !
    const xmlChar *	name	: the name of the property
    struct _xmlNode *	children	: the value of the property
    struct _xmlNode *	last	: NULL
    struct _xmlNode *	parent	: child->parent link
    struct _xmlAttr *	next	: next sibling link
    struct _xmlAttr *	prev	: previous sibling link
    struct _xmlDoc *	doc	: the containing document
    xmlNs *	ns	: pointer to the associated namespace
    xmlAttributeType	atype	: the attribute type if validating
    void *	psvi	: for type/PSVI informations


struct _xmlElement {
    void *	_private	: application data
    xmlElementType	type	: XML_ELEMENT_DECL, must be second !
    const xmlChar *	name	: Element name
    struct _xmlNode *	children	: NULL
    struct _xmlNode *	last	: NULL
    struct _xmlDtd *	parent	: -> DTD
    struct _xmlNode *	next	: next sibling link
    struct _xmlNode *	prev	: previous sibling link
    struct _xmlDoc *	doc	: the containing document
    xmlElementTypeVal	etype	: The type
    xmlElementContentPtr	content	: the allowed element content
    xmlAttributePtr	attributes	: List of the declared attributes
    const xmlChar *	prefix	: the namespace prefix if any
    xmlRegexpPtr	contModel	: the validating regexp
    void *	contModel
}



*/


//    htmlEncodeEntities(output, &outlen, att, &attlen, '\'');
void    xxhtmlEncodeEntities(unsigned char * output, int * outlen, const unsigned char * att, int * attlen, char s)
{

    //htmlEncodeEntities(output, outlen, att, attlen, s);
*outlen=*attlen;
memcpy(output, att, *attlen);
output[*attlen]=0; 
printf(" xxhtmlEncode %s\n",att);
}


///////////////////////////////////////////////////////

xmlSAXHandler emptySAXHandlerStruct = {
    NULL, /* internalSubset */
    NULL, /* isStandalone */
    NULL, /* hasInternalSubset */
    NULL, /* hasExternalSubset */
    NULL, /* resolveEntity */
    NULL, /* getEntity */
    NULL, /* entityDecl */
    NULL, /* notationDecl */
    NULL, /* attributeDecl */
    NULL, /* elementDecl */
    NULL, /* unparsedEntityDecl */
    NULL, /* setDocumentLocator */
    NULL, /* startDocument */
    NULL, /* endDocument */
    NULL, /* startElement */
    NULL, /* endElement */
    NULL, /* reference */
    NULL, /* characters */
    NULL, /* ignorableWhitespace */
    NULL, /* processingInstruction */
    NULL, /* comment */
    NULL, /* xmlParserWarning */
    NULL, /* xmlParserError */
    NULL, /* xmlParserError */
    NULL, /* getParameterEntity */
};

xmlSAXHandlerPtr emptySAXHandler = &emptySAXHandlerStruct;
extern xmlSAXHandlerPtr debugSAXHandler;

/************************************************************************
 *									*
 *				Debug Handlers				*
 *									*
 ************************************************************************/

/**
 * isStandaloneDebug:
 * @ctxt:  An XML parser context
 *
 * Is this document tagged standalone ?
 *
 * Returns 1 if true
 */
int
isStandaloneDebug(void *ctx)
{
    fprintf(stdout, "SAX.isStandalone()\n");
    return(0);
}

/**
 * hasInternalSubsetDebug:
 * @ctxt:  An XML parser context
 *
 * Does this document has an internal subset
 *
 * Returns 1 if true
 */
int
hasInternalSubsetDebug(void *ctx)
{
    fprintf(stdout, "SAX.hasInternalSubset()\n");
    return(0);
}

/**
 * hasExternalSubsetDebug:
 * @ctxt:  An XML parser context
 *
 * Does this document has an external subset
 *
 * Returns 1 if true
 */
int
hasExternalSubsetDebug(void *ctx)
{
    fprintf(stdout, "SAX.hasExternalSubset()\n");
    return(0);
}

/**
 * hasInternalSubsetDebug:
 * @ctxt:  An XML parser context
 *
 * Does this document has an internal subset
 */
void
internalSubsetDebug(void *ctx, const xmlChar *name,
	       const xmlChar *ExternalID, const xmlChar *SystemID)
{
    fprintf(stdout, "SAX.internalSubset(%s,", name);
    if (ExternalID == NULL)
	fprintf(stdout, " ,");
    else
	fprintf(stdout, " %s,", ExternalID);
    if (SystemID == NULL)
	fprintf(stdout, " )\n");
    else
	fprintf(stdout, " %s)\n", SystemID);
}

/**
 * resolveEntityDebug:
 * @ctxt:  An XML parser context
 * @publicId: The public ID of the entity
 * @systemId: The system ID of the entity
 *
 * Special entity resolver, better left to the parser, it has
 * more context than the application layer.
 * The default behaviour is to NOT resolve the entities, in that case
 * the ENTITY_REF nodes are built in the structure (and the parameter
 * values).
 *
 * Returns the xmlParserInputPtr if inlined or NULL for DOM behaviour.
 */
xmlParserInputPtr
resolveEntityDebug(void *ctx, const xmlChar *publicId, const xmlChar *systemId)
{
    /* xmlParserCtxtPtr ctxt = (xmlParserCtxtPtr) ctx; */

    
    fprintf(stdout, "SAX.resolveEntity(");
    if (publicId != NULL)
	fprintf(stdout, "%s", (char *)publicId);
    else
	fprintf(stdout, " ");
    if (systemId != NULL)
	fprintf(stdout, ", %s)\n", (char *)systemId);
    else
	fprintf(stdout, ", )\n");
/*********
    if (systemId != NULL) {
        return(xmlNewInputFromFile(ctxt, (char *) systemId));
    }
 *********/
    return(NULL);
}

/**
 * getEntityDebug:
 * @ctxt:  An XML parser context
 * @name: The entity name
 *
 * Get an entity by name
 *
 * Returns the xmlParserInputPtr if inlined or NULL for DOM behaviour.
 */
xmlEntityPtr
getEntityDebug(void *ctx, const xmlChar *name)
{
    fprintf(stdout, "SAX.getEntity(%s)\n", name);
    return(NULL);
}

/**
 * getParameterEntityDebug:
 * @ctxt:  An XML parser context
 * @name: The entity name
 *
 * Get a parameter entity by name
 *
 * Returns the xmlParserInputPtr
 */
xmlEntityPtr
getParameterEntityDebug(void *ctx, const xmlChar *name)
{
    fprintf(stdout, "SAX.getParameterEntity(%s)\n", name);
    return(NULL);
}


/**
 * entityDeclDebug:
 * @ctxt:  An XML parser context
 * @name:  the entity name 
 * @type:  the entity type 
 * @publicId: The public ID of the entity
 * @systemId: The system ID of the entity
 * @content: the entity value (without processing).
 *
 * An entity definition has been parsed
 */
void
entityDeclDebug(void *ctx, const xmlChar *name, int type,
          const xmlChar *publicId, const xmlChar *systemId, xmlChar *content)
{
    fprintf(stdout, "SAX.entityDecl(%s, %d, %s, %s, %s)\n",
            name, type, publicId, systemId, content);

    fprintf(stdout, "wtf is this %s \n", content);

}

/**
 * attributeDeclDebug:
 * @ctxt:  An XML parser context
 * @name:  the attribute name 
 * @type:  the attribute type 
 *
 * An attribute definition has been parsed
 */
void
attributeDeclDebug(void *ctx, const xmlChar *elem, const xmlChar *name,
              int type, int def, const xmlChar *defaultValue,
	      xmlEnumerationPtr tree)
{
    fprintf(stdout, "SAX.attributeDecl(%s, %s, %d, %d, %s, ...)\n",
            elem, name, type, def, defaultValue);
}

/**
 * elementDeclDebug:
 * @ctxt:  An XML parser context
 * @name:  the element name 
 * @type:  the element type 
 * @content: the element value (without processing).
 *
 * An element definition has been parsed
 */
void
elementDeclDebug(void *ctx, const xmlChar *name, int type,
	    xmlElementContentPtr content)
{
    fprintf(stdout, "SAX.elementDecl(%s, %d, ...)\n",
            name, type);
//  this is toin to take foruck forever for this stht  
    fprintf(stdout, " decldebug wtf is this %s \n", content->name);
}

/**
 * notationDeclDebug:
 * @ctxt:  An XML parser context
 * @name: The name of the notation
 * @publicId: The public ID of the entity
 * @systemId: The system ID of the entity
 *
 * What to do when a notation declaration has been parsed.
 */
void
notationDeclDebug(void *ctx, const xmlChar *name,
	     const xmlChar *publicId, const xmlChar *systemId)
{
    fprintf(stdout, "SAX.notationDecl(%s, %s, %s)\n",
            (char *) name, (char *) publicId, (char *) systemId);
}

/**
 * unparsedEntityDeclDebug:
 * @ctxt:  An XML parser context
 * @name: The name of the entity
 * @publicId: The public ID of the entity
 * @systemId: The system ID of the entity
 * @notationName: the name of the notation
 *
 * What to do when an unparsed entity declaration is parsed
 */
void
unparsedEntityDeclDebug(void *ctx, const xmlChar *name,
		   const xmlChar *publicId, const xmlChar *systemId,
		   const xmlChar *notationName)
{
    fprintf(stdout, "SAX.unparsedEntityDecl(%s, %s, %s, %s)\n",
            (char *) name, (char *) publicId, (char *) systemId,
	    (char *) notationName);
}

/**
 * setDocumentLocatorDebug:
 * @ctxt:  An XML parser context
 * @loc: A SAX Locator
 *
 * Receive the document locator at startup, actually xmlDefaultSAXLocator
 * Everything is available on the context, so this is useless in our case.
 */
void
setDocumentLocatorDebug(void *ctx, xmlSAXLocatorPtr loc)
{
    fprintf(stdout, "SAX.setDocumentLocator()\n");
}

/**
 * startDocumentDebug:
 * @ctxt:  An XML parser context
 *
 * called when the document start being processed.
 */
void
startDocumentDebug(void *ctx)
{
    fprintf(stdout, "SAX.startDocument()\n");
}

/**
 * endDocumentDebug:
 * @ctxt:  An XML parser context
 *
 * called when the document end has been detected.
 */
void
endDocumentDebug(void *ctx)
{
    fprintf(stdout, "SAX.endDocument()\n");
}

/**
 * startElementDebug:
 * @ctxt:  An XML parser context
 * @name:  The element name
 *
 * called when an opening tag has been processed.
 */
void
startElementDebug(void *ctx, const xmlChar *name, const xmlChar **atts)
{
    int i;

    fprintf(stdout, "SAX.startElement(%s", (char *) name);
    if (atts != NULL) {
        for (i = 0;(atts[i] != NULL);i++) {
	    fprintf(stdout, ", %s", atts[i++]);
	    if (atts[i] != NULL) {
		unsigned char output[400000];
		const unsigned char *att = atts[i];
		int outlen, attlen;
	        fprintf(stdout, "='");
		while ((attlen = strlen((char*)att)) > 0) {
		    outlen = sizeof output - 1;
		    xxhtmlEncodeEntities(output, &outlen, att, &attlen, '\'');
		    //fprintf(stdout, "%.*s", outlen, output);
		    fprintf(stdout, "%.*s", outlen, output);
		    fprintf(stdout, "att char %s\n", att);
		    att += attlen;
		}
		fprintf(stdout, "'");
	    }
	}
    }
    fprintf(stdout, ")\n");
}

/**
 * endElementDebug:
 * @ctxt:  An XML parser context
 * @name:  The element name
 *
 * called when the end of an element has been detected.
 */
void
endElementDebug(void *ctx, const xmlChar *name)
{
    fprintf(stdout, "SAX.endElement(%s)\n", (char *) name);
}

/**
 * charactersDebug:
 * @ctxt:  An XML parser context
 * @ch:  a xmlChar string
 * @len: the number of xmlChar
 *
 * receiving some chars from the parser.
 * Question: how much at a time ???
 */
void
charactersDebug(void *ctx, const xmlChar *ch, int len)
{
    unsigned char output[400000];
    //int inlen = len, outlen = 30;
    int inlen = len, outlen = sizeof output - 1 ;

    xxhtmlEncodeEntities(output, &outlen, ch, &inlen, 0);
    output[outlen] = 0;

    fprintf(stdout, "SAX.characters(%s, %d)\n", output, len);
}

/**
 * cdataDebug:
 * @ctxt:  An XML parser context
 * @ch:  a xmlChar string
 * @len: the number of xmlChar
 *
 * receiving some cdata chars from the parser.
 * Question: how much at a time ???
 */
void
cdataDebug(void *ctx, const xmlChar *ch, int len)
{
    unsigned char output[400000];
    //int inlen = len, outlen = 30;
    int inlen = len, outlen = sizeof output - 1;

    xxhtmlEncodeEntities(output, &outlen, ch, &inlen, 0);
    output[outlen] = 0;

    fprintf(stdout, "SAX.cdata(%s, %d)\n", output, len);
}

/**
 * referenceDebug:
 * @ctxt:  An XML parser context
 * @name:  The entity name
 *
 * called when an entity reference is detected. 
 */
void
referenceDebug(void *ctx, const xmlChar *name)
{
    fprintf(stdout, "SAX.reference(%s)\n", name);
}

/**
 * ignorableWhitespaceDebug:
 * @ctxt:  An XML parser context
 * @ch:  a xmlChar string
 * @start: the first char in the string
 * @len: the number of xmlChar
 *
 * receiving some ignorable whitespaces from the parser.
 * Question: how much at a time ???
 */
void
ignorableWhitespaceDebug(void *ctx, const xmlChar *ch, int len)
{
    char output[400000];
    int i;

    for (i = 0;(i<len) && (i < 30);i++)
	output[i] = ch[i];
    output[i] = 0;

    fprintf(stdout, "SAX.ignorableWhitespace(%s, %d)\n", output, len);
}

/**
 * processingInstructionDebug:
 * @ctxt:  An XML parser context
 * @target:  the target name
 * @data: the PI data's
 * @len: the number of xmlChar
 *
 * A processing instruction has been parsed.
 */
void
processingInstructionDebug(void *ctx, const xmlChar *target,
                      const xmlChar *data)
{
    fprintf(stdout, "SAX.processingInstruction(%s, %s)\n",
            (char *) target, (char *) data);
}

/**
 * commentDebug:
 * @ctxt:  An XML parser context
 * @value:  the comment content
 *
 * A comment has been parsed.
 */
void
commentDebug(void *ctx, const xmlChar *value)
{
    fprintf(stdout, "SAX.comment(%s)\n", value);
}

/**
 * warningDebug:
 * @ctxt:  An XML parser context
 * @msg:  the message to display/transmit
 * @...:  extra parameters for the message display
 *
 * Display and format a warning messages, gives file, line, position and
 * extra parameters.
 */
void
warningDebug(void *ctx, const char *msg, ...)
{
    va_list args;

    va_start(args, msg);
    fprintf(stdout, "SAX.warning: ");
    vfprintf(stdout, msg, args);
    va_end(args);
}

/**
 * errorDebug:
 * @ctxt:  An XML parser context
 * @msg:  the message to display/transmit
 * @...:  extra parameters for the message display
 *
 * Display and format a error messages, gives file, line, position and
 * extra parameters.
 */
void
errorDebug(void *ctx, const char *msg, ...)
{
    va_list args;

    va_start(args, msg);
    fprintf(stdout, "SAX.error: ");
    vfprintf(stdout, msg, args);
    va_end(args);
}

/**
 * fatalErrorDebug:
 * @ctxt:  An XML parser context
 * @msg:  the message to display/transmit
 * @...:  extra parameters for the message display
 *
 * Display and format a fatalError messages, gives file, line, position and
 * extra parameters.
 */
void
fatalErrorDebug(void *ctx, const char *msg, ...)
{
    va_list args;

    va_start(args, msg);
    fprintf(stdout, "SAX.fatalError: ");
    vfprintf(stdout, msg, args);
    va_end(args);
}

xmlSAXHandler debugSAXHandlerStruct = {
    internalSubsetDebug,
    isStandaloneDebug,
    hasInternalSubsetDebug,
    hasExternalSubsetDebug,
    resolveEntityDebug,
    getEntityDebug,
    entityDeclDebug,
    notationDeclDebug,
    attributeDeclDebug,
    elementDeclDebug,
    unparsedEntityDeclDebug,
    setDocumentLocatorDebug,
    startDocumentDebug,
    endDocumentDebug,
    startElementDebug,
    endElementDebug,
    referenceDebug,
    charactersDebug,
    ignorableWhitespaceDebug,
    processingInstructionDebug,
    commentDebug,
    warningDebug,
    errorDebug,
    fatalErrorDebug,
    getParameterEntityDebug,
    cdataDebug,
    NULL
};

xmlSAXHandlerPtr debugSAXHandler = &debugSAXHandlerStruct;
/************************************************************************
 *									*
 *				Debug					*
 *									*
 ************************************************************************/

void parseSAXFile(char *filename) {
    htmlDocPtr doc = NULL;

    /*
     * Empty callbacks for checking
     */
    if (push) {
	FILE *f;

	f = fopen(filename, "r");
	if (f != NULL) {
	    int res, size = 3;
	    char chars[4096];
	    htmlParserCtxtPtr ctxt;

	    /* if (repeat) */
		size = 4096;
	    res = fread(chars, 1, 4, f);
	    if (res > 0) {
		ctxt = htmlCreatePushParserCtxt(emptySAXHandler, NULL,
			    chars, res, filename, XML_CHAR_ENCODING_NONE);
		while ((res = fread(chars, 1, size, f)) > 0) {
		    htmlParseChunk(ctxt, chars, res, 0);
		}
		htmlParseChunk(ctxt, chars, 0, 1);
		doc = ctxt->myDoc;
		htmlFreeParserCtxt(ctxt);
	    }
	    if (doc != NULL) {
		fprintf(stdout, "htmlSAXParseFile returned non-NULL\n");
		xmlFreeDoc(doc);
	    }
	    fclose(f);
	}
	if (!noout) {
	    f = fopen(filename, "r");
	    if (f != NULL) {
		int res, size = 3;
		char chars[4096];
		htmlParserCtxtPtr ctxt;

		/* if (repeat) */
		    size = 4096;
		res = fread(chars, 1, 4, f);
		if (res > 0) {
		    ctxt = htmlCreatePushParserCtxt(debugSAXHandler, NULL,
				chars, res, filename, XML_CHAR_ENCODING_NONE);
		    while ((res = fread(chars, 1, size, f)) > 0) {
			htmlParseChunk(ctxt, chars, res, 0);
		    }
		    htmlParseChunk(ctxt, chars, 0, 1);
		    doc = ctxt->myDoc;
		    htmlFreeParserCtxt(ctxt);
		}
		if (doc != NULL) {
		    fprintf(stdout, "htmlSAXParseFile returned non-NULL\n");
		    xmlFreeDoc(doc);
		}
		fclose(f);
	    }
	}
    } else {	
	doc = htmlSAXParseFile(filename, NULL, emptySAXHandler, NULL);
	if (doc != NULL) {
	    fprintf(stdout, "htmlSAXParseFile returned non-NULL\n");
	    xmlFreeDoc(doc);
	}

	if (!noout) {
	    /*
	     * Debug callback
	     */
	    doc = htmlSAXParseFile(filename, NULL, debugSAXHandler, NULL);
	    if (doc != NULL) {
		fprintf(stdout, "htmlSAXParseFile returned non-NULL\n");
		xmlFreeDoc(doc);
	    }
	}
    }
}

#define Bit(f,i) (((f>>i)&1)!=0)
 
void parseAndPrintFile(char *filename) {
    htmlDocPtr doc = NULL, tmp;

    /*
     * build an HTML tree from a string;
     */
    if (push) {
	FILE *f;

	f = fopen(filename, "r");
	if (f != NULL) {
	    int res, size = 3;
	    char chars[4096];
	    htmlParserCtxtPtr ctxt;

	    /* if (repeat) */
		size = 4096;
	    res = fread(chars, 1, 4, f);
	    if (res > 0) {
		ctxt = htmlCreatePushParserCtxt(NULL, NULL,
			    chars, res, filename, XML_CHAR_ENCODING_NONE);
		while ((res = fread(chars, 1, size, f)) > 0) {
		    htmlParseChunk(ctxt, chars, res, 0);
		}
		htmlParseChunk(ctxt, chars, 0, 1);
		doc = ctxt->myDoc;
		htmlFreeParserCtxt(ctxt);
	    }
	    fclose(f);
	}
    } else {	
	doc = htmlParseFile(filename, NULL);
    }
    if (doc == NULL) {
        xmlGenericError(xmlGenericErrorContext,
		"Could not parse %s\n", filename);
    }

    /*
     * test intermediate copy if needed.
     */
    if (copy) {
        tmp = doc;
	doc = xmlCopyDoc(doc, 1);
	xmlFreeDoc(tmp);
    }

    /*
     * print it.
     */
if (mjm_print)
{
form_spec fs;
form_spec * fsp= &fs;
if (mjm_flags!=0)
{

fsp->stack=Bit(mjm_flags,0)?0:1;
fsp->print_stacked=Bit(mjm_flags,1)?0:1;
fsp->print_curr=Bit(mjm_flags,2)?1:0;
fsp->remove_crlf=Bit(mjm_flags,3)?0:1;
fsp->include_node_numbers=Bit(mjm_flags,4)?1:0;
fsp->suppress_levels=Bit(mjm_flags,5)?1:0;
fprintf(stderr," fsp : stack %d print_stacked %d print_curr %d remove_Crlf %d include_node_numbers %d suppress_levels %d\n",fsp->stack, 
fsp->print_stacked, fsp->print_curr, fsp->remove_crlf
,fsp->include_node_numbers, fsp->suppress_levels);
}
if (mjm_sflags!=0)
{
if (strcmp(mjm_sflags,"hier")==0) 
{

fsp->stack=1;
fsp->print_stacked=1;
fsp->print_curr=0;
fsp->remove_crlf=1;
fsp->include_node_numbers=1;
fsp->suppress_levels=1;
fprintf(stderr," fsp : stack %d print_stacked %d print_curr %d remove_Crlf %d include_node_numbers %d suppress_levels %d\n",fsp->stack, 
fsp->print_stacked, fsp->print_curr, fsp->remove_crlf
,fsp->include_node_numbers, fsp->suppress_levels);

} // hier

} // mjm_sflags



#undef Bit
print_doc(doc,fsp);

} else
    if (!noout) { 
#ifdef LIBXML_DEBUG_ENABLED
	if (!debug) {
	    if (encoding)
		htmlSaveFileEnc("-", doc, encoding);
	    else
		htmlDocDump(stdout, doc);
	} else
	    xmlDebugDumpDocument(stdout, doc);
#else
	if (encoding)
	    htmlSaveFileEnc("-", doc, encoding);
	else
	    htmlDocDump(stdout, doc);
#endif
    }	

    /*
     * free it.
     */
    xmlFreeDoc(doc);
}
#ifndef TESTHTML_HEADER

int main(int argc, char **argv) {
    int i, count;
    int files = 0;
if (argc==0)
{
fprintf(stderr,"mjm  modified testHTML.c mjm %s %s\n",__FILE__,__DATE__);
fprintf(stderr," downloaded from https://dev.w3.org/XML/testHTML.c\n"); 
}

    for (i = 1; i < argc ; i++) {
#ifdef LIBXML_DEBUG_ENABLED
	if ((!strcmp(argv[i], "-debug")) || (!strcmp(argv[i], "--debug")))
	    debug++;
	else
#endif
	    if ((!strcmp(argv[i], "-copy")) || (!strcmp(argv[i], "--copy"))) copy++;
	else if ((!strcmp(argv[i], "-push")) || (!strcmp(argv[i], "--push"))) push++;
	else if ((!strcmp(argv[i], "-quiet")) || (!strcmp(argv[i], "--quiet"))) //push++;
 { quiet_mjm_xml(); } 
	else if ((!strcmp(argv[i], "-mjm")) || (!strcmp(argv[i], "--mjm"))){ mjm_print++;
//quiet_mjm_xml();
}

	else if ((!strcmp(argv[i], "-mjmflags")) || (!strcmp(argv[i], "--mjmflags"))) 
	{ ++i; if (i<argc) 
{
		 mjm_flags=atoi(argv[i]); 
		int len=strlen(argv[i])+2;
		mjm_sflags= new char[len];
	strcpy(mjm_sflags,argv[i]); 

} /// i<argc
} // --mjmflags 
	else if ((!strcmp(argv[i], "-sax")) || (!strcmp(argv[i], "--sax")))
	    sax++;
	else if ((!strcmp(argv[i], "-noout")) || (!strcmp(argv[i], "--noout")))
	    noout++;
	else if ((!strcmp(argv[i], "-repeat")) ||
	         (!strcmp(argv[i], "--repeat")))
	    repeat++;
	else if ((!strcmp(argv[i], "-encode")) ||
	         (!strcmp(argv[i], "--encode"))) {
	    i++;
	    encoding = argv[i];
        }
    }
    for (i = 1; i < argc ; i++) {
	if ((!strcmp(argv[i], "-encode")) ||
	         (!strcmp(argv[i], "--encode"))) {
	    i++;
	    continue;
        }
	if ((!strcmp(argv[i], "-mjmflags")) || (!strcmp(argv[i], "--mjmflags"))) { ++i; continue; }  
	if ((argv[i][0] != '-')||(argv[i][0]=='-')&&(strlen(argv[i])==1)) {
	    if (repeat) {
		for (count = 0;count < 100 * repeat;count++) {
		    if (sax)
			parseSAXFile(argv[i]);
		    else   
			parseAndPrintFile(argv[i]);
		}    
	    } else {
		if (sax)
		    parseSAXFile(argv[i]);
		else   
		    parseAndPrintFile(argv[i]);
	    }
	    files ++;
	}
    }
    if (files == 0) {
	printf("Usage : %s [--debug] [--copy] [--copy] HTMLfiles ...\n",
	       argv[0]);
	printf("\tParse the HTML files and output the result of the parsing\n");
#ifdef LIBXML_DEBUG_ENABLED
	printf("\t--debug : dump a debug tree of the in-memory document\n");
#endif
	printf("\t--copy : used to test the internal copy implementation\n");
	printf("\t--sax : debug the sequence of SAX callbacks\n");
	printf("\t--repeat : parse the file 100 times, for timing\n");
	printf("\t--noout : do not print the result\n");
	printf("\t--push : use the push mode parser\n");
	printf("\t--encode encoding : output in the given encoding\n");
    }
    xmlCleanupParser();
    xmlMemoryDump();
	delete [] mjm_sflags;
    return(0);
}
#endif //  TESTHTML_HEADER
#else /* !LIBXML_HTML_ENABLED */
#include <stdio.h>
int main(int argc, char **argv) {
    printf("%s : HTML support not compiled in\n", argv[0]);
    return(0);
}
#endif
