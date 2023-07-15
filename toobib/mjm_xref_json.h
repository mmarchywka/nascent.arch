#ifndef MJM_XREF_JSON_H__
#define MJM_XREF_JSON_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_hier_two.h"

#include <../mjsonu/mjsonu.h>



#include "mjm_pawnoff.h"
#include "mjm_collections.h"
#include "mjm_wovdb.h"
#include "mjm_strings.h"

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>

/*
Add datacite look up,
i

# GET /dois
$ curl https://api.test.datacite.org/dois/10.5438/0012
$ curl https://api.datacite.org/dois/10.5438/0012
https://support.datacite.org/docs/api-get-doi


{
  "data": {
    "id": "10.5438/0012",
    "type": "dois",
    "attributes": {
      "doi": "10.5438/0012",
      "prefix": "10.5438",
      "suffix": "0012",
      "identifiers": [
        {
          "identifier": "https://doi.org/10.5438/0012",
          "identifierType": "DOI"
        }
      ],
      "creators": [
        {
          "name": "DataCite Metadata Working Group"
        }
      ],
      "titles": [
        {
          "title": "DataCite Metadata Schema Documentation for the Publication and Citation of Research Data v4.0"
        }
      ],
      "publisher": "DataCite e.V.",
      "container": {},
      "publicationYear": 2016,
      "subjects": [],
      "contributors": [
        {
          "name": "Starr, Joan",
          "nameType": "Personal",
          "givenName": "Joan",
          "familyName": "Starr",
          "affiliation": "California Digital Library",
          "contributorType": "ProjectLeader",
          "nameIdentifiers": [
            {
              "nameIdentifier": "https://orcid.org/0000-0002-7285-027X",
              "nameIdentifierScheme": "ORCID"
            }
          ]
        },
        {
          "name": "Smaele, Madeleine de",
          "nameType": "Personal",
          "givenName": "Madeleine de",
          "familyName": "Smaele",
          "affiliation": "TU Delft",
          "contributorType": "ProjectLeader"
        },
        {
          "name": "Ashton, Jan",
          "nameType": "Personal",
          "givenName": "Jan",
          "familyName": "Ashton",
          "affiliation": "British Library",
          "contributorType": "Editor"
        },
        {
          "name": "Barton, Amy",
          "nameType": "Personal",
          "givenName": "Amy",
          "familyName": "Barton",
          "affiliation": "Purdue University Library",
          "contributorType": "Editor"
        },
        {
          "name": "Bradford, Tina",
          "nameType": "Personal",
          "givenName": "Tina",
          "familyName": "Bradford",
          "affiliation": "NRC/CISTI",
          "contributorType": "Editor"
        },
        {
          "name": "Ciolek-Figiel, Anne",
          "nameType": "Personal",
          "givenName": "Anne",
          "familyName": "Ciolek-Figiel",
          "affiliation": "Inist‐CNRS",
          "contributorType": "Editor"
        },
        {
          "name": "Dietiker, Stefanie",
          "nameType": "Personal",
          "givenName": "Stefanie",
          "familyName": "Dietiker",
          "affiliation": "ETH Zürich",
          "contributorType": "Editor"
        },
        {
          "name": "Elliot, Jannean",
          "nameType": "Personal",
          "givenName": "Jannean",
          "familyName": "Elliot",
          "affiliation": "DOE/OSTI",
          "contributorType": "Editor"
        },
        {
          "name": "Genat, Berrit",
          "nameType": "Personal",
          "givenName": "Berrit",
          "familyName": "Genat",
          "affiliation": "TIB",
          "contributorType": "Editor"
        },
        {
          "name": "Harzenetter, Karoline",
          "nameType": "Personal",
          "givenName": "Karoline",
          "familyName": "Harzenetter",
          "affiliation": "GESIS",
          "contributorType": "Editor"
        },
        {
          "name": "Hirschmann, Barbara",
          "nameType": "Personal",
          "givenName": "Barbara",
          "familyName": "Hirschmann",
          "affiliation": "ETH Zürich",
          "contributorType": "Editor",
          "nameIdentifiers": [
            {
              "nameIdentifier": "https://orcid.org/0000-0003-0289-0345",
              "nameIdentifierScheme": "ORCID"
            }
          ]
        },
        {
          "name": "Jakobsson, Stefan",
          "nameType": "Personal",
          "givenName": "Stefan",
          "familyName": "Jakobsson",
          "affiliation": "SND",
          "contributorType": "Editor"
        },
        {
          "name": "Mailloux, Jean-Yves",
          "nameType": "Personal",
          "givenName": "Jean-Yves",
          "familyName": "Mailloux",
          "affiliation": "NRC/CISTI",
          "contributorType": "Editor"
        },
        {
          "name": "Newbold, Elizabeth",
          "nameType": "Personal",
          "givenName": "Elizabeth",
          "familyName": "Newbold",
          "affiliation": "British Library",
          "contributorType": "Editor",
          "nameIdentifiers": [
            {
              "nameIdentifier": "https://orcid.org/0000-0002-8255-9013",
              "nameIdentifierScheme": "ORCID"
            }
          ]
        },
        {
          "name": "Nielsen, Lars Holm",
          "nameType": "Personal",
          "givenName": "Lars Holm",
          "familyName": "Nielsen",
          "affiliation": "CERN",
          "contributorType": "Editor",
          "nameIdentifiers": [
            {
              "nameIdentifier": "https://orcid.org/0000-0001-8135-3489",
              "nameIdentifierScheme": "ORCID"
            }
          ]
        },
        {
          "name": "Yahia, Mohamed",
          "nameType": "Personal",
          "givenName": "Mohamed",
          "familyName": "Yahia",
          "affiliation": "Inist-CNRS",
          "contributorType": "Editor"
        },
        {
          "name": "Ziedorn, Frauke",
          "nameType": "Personal",
          "givenName": "Frauke",
          "familyName": "Ziedorn",
          "affiliation": "TIB",
          "contributorType": "Supervisor",
          "nameIdentifiers": [
            {
              "nameIdentifier": "https://orcid.org/0000-0002-1143-781X",
              "nameIdentifierScheme": "ORCID"
            }
          ]
        }
      ],
      "dates": [
        {
          "date": "2016",
          "dateType": "Issued"
        }
      ],
      "language": "eng",
      "types": {
        "ris": "RPRT",
        "bibtex": "article",
        "citeproc": "article-journal",
        "schemaOrg": "ScholarlyArticle",
        "resourceType": "Documentation",
        "resourceTypeGeneral": "Text"
      },
      "relatedIdentifiers": [
        {
          "relationType": "Documents",
          "relatedIdentifier": "10.5438/0013",
          "relatedIdentifierType": "DOI"
        },
        {
          "relationType": "IsNewVersionOf",
          "relatedIdentifier": "10.5438/0010",
          "relatedIdentifierType": "DOI"
        }
      ],
      "sizes": [
        "45 pages"
      ],
      "formats": [
        "application/pdf"
      ],
      "version": null,
      "rightsList": [],
      "descriptions": [
        {
          "description": "1 Introduction\n1.1 The DataCite Consortium\n1.2 DataCite Community Participation\n1.3 The Metadata Schema\n1.4 Version 4.0 Update\n2 DataCite Metadata Properties\n2.1 Overview\n2.2 Citation\n2.3 DataCite Properties\n3 XML Example\n4 XML Schema\n5 Other DataCite Services\nAppendices\nAppendix 1: Controlled List Definitions\nAppendix 2: Earlier Version Update Notes",
          "descriptionType": "TableOfContents"
        }
      ],
      "geoLocations": [],
      "fundingReferences": [],
      "xml": "PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4NCjxyZXNvdXJjZSB4bWxucz0iaHR0cDovL2RhdGFjaXRlLm9yZy9zY2hlbWEva2VybmVsLTQiIHhtbG5zOnhzaT0iaHR0cDovL3d3dy53My5vcmcvMjAwMS9YTUxTY2hlbWEtaW5zdGFuY2UiIHhzaTpzY2hlbWFMb2NhdGlvbj0iaHR0cDovL2RhdGFjaXRlLm9yZy9zY2hlbWEva2VybmVsLTQgaHR0cDovL3NjaGVtYS5kYXRhY2l0ZS5vcmcvbWV0YS9rZXJuZWwtNC9tZXRhZGF0YS54c2QiPg0KCTxpZGVudGlmaWVyIGlkZW50aWZpZXJUeXBlPSJET0kiPjEwLjU0MzgvMDAxMjwvaWRlbnRpZmllcj4NCgk8Y3JlYXRvcnM+DQoJCTxjcmVhdG9yPg0KCQkJPGNyZWF0b3JOYW1lPkRhdGFDaXRlIE1ldGFkYXRhIFdvcmtpbmcgR3JvdXA8L2NyZWF0b3JOYW1lPg0KCQk8L2NyZWF0b3I+DQoJPC9jcmVhdG9ycz4NCgk8dGl0bGVzPg0KCQk8dGl0bGU+RGF0YUNpdGUgTWV0YWRhdGEgU2NoZW1hIERvY3VtZW50YXRpb24gZm9yIHRoZSBQdWJsaWNhdGlvbiBhbmQgQ2l0YXRpb24gb2YgUmVzZWFyY2ggRGF0YSB2NC4wPC90aXRsZT4NCgk8L3RpdGxlcz4NCgk8cHVibGlzaGVyPkRhdGFDaXRlIGUuVi48L3B1Ymxpc2hlcj4NCgk8cHVibGljYXRpb25ZZWFyPjIwMTY8L3B1YmxpY2F0aW9uWWVhcj4NCgk8Y29udHJpYnV0b3JzPg0KCQk8Y29udHJpYnV0b3IgY29udHJpYnV0b3JUeXBlPSJQcm9qZWN0TGVhZGVyIj4NCgkJCTxjb250cmlidXRvck5hbWU+U3RhcnIsIEpvYW48L2NvbnRyaWJ1dG9yTmFtZT4NCgkJCTxnaXZlbk5hbWU+Sm9hbjwvZ2l2ZW5OYW1lPg0KCQkJPGZhbWlseU5hbWU+U3RhcnI8L2ZhbWlseU5hbWU+DQoJCQk8bmFtZUlkZW50aWZpZXIgbmFtZUlkZW50aWZpZXJTY2hlbWU9Ik9SQ0lEIiBzY2hlbWVVUkk9Imh0dHA6Ly9vcmNpZC5vcmciPjAwMDAtMDAwMi03Mjg1LTAyN1g8L25hbWVJZGVudGlmaWVyPg0KCQkJPGFmZmlsaWF0aW9uPkNhbGlmb3JuaWEgRGlnaXRhbCBMaWJyYXJ5PC9hZmZpbGlhdGlvbj4NCgkJPC9jb250cmlidXRvcj4NCgkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iUHJvamVjdExlYWRlciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPlNtYWVsZSwgTWFkZWxlaW5lIGRlPC9jb250cmlidXRvck5hbWU+DQoJCQk8Z2l2ZW5OYW1lPk1hZGVsZWluZSBkZTwvZ2l2ZW5OYW1lPg0KCQkJPGZhbWlseU5hbWU+U21hZWxlPC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPlRVIERlbGZ0PC9hZmZpbGlhdGlvbj4NCgkJPC9jb250cmlidXRvcj4NCgkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iRWRpdG9yIj4NCgkJCTxjb250cmlidXRvck5hbWU+QXNodG9uLCBKYW48L2NvbnRyaWJ1dG9yTmFtZT4NCgkJCTxnaXZlbk5hbWU+SmFuPC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5Bc2h0b248L2ZhbWlseU5hbWU+DQoJCQk8YWZmaWxpYXRpb24+QnJpdGlzaCBMaWJyYXJ5PC9hZmZpbGlhdGlvbj4NCgkJPC9jb250cmlidXRvcj4NCgkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iRWRpdG9yIj4NCgkJCTxjb250cmlidXRvck5hbWU+QmFydG9uLCBBbXk8L2NvbnRyaWJ1dG9yTmFtZT4NCgkJCTxnaXZlbk5hbWU+QW15PC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5CYXJ0b248L2ZhbWlseU5hbWU+DQoJCQk8YWZmaWxpYXRpb24+UHVyZHVlIFVuaXZlcnNpdHkgTGlicmFyeTwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPkJyYWRmb3JkLCBUaW5hPC9jb250cmlidXRvck5hbWU+DQoJCQk8Z2l2ZW5OYW1lPlRpbmE8L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPkJyYWRmb3JkPC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPk5SQy9DSVNUSTwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPkNpb2xla+KAkEZpZ2llbCwgQW5uZTwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5Bbm5lPC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5DaW9sZWstRmlnaWVsPC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPkluaXN04oCQQ05SUzwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPkRpZXRpa2VyLCBTdGVmYW5pZTwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5TdGVmYW5pZTwvZ2l2ZW5OYW1lPg0KCQkJPGZhbWlseU5hbWU+RGlldGlrZXI8L2ZhbWlseU5hbWU+DQoJCQk8YWZmaWxpYXRpb24+RVRIIFrDvHJpY2g8L2FmZmlsaWF0aW9uPg0KCQk8L2NvbnRyaWJ1dG9yPg0KCQk8Y29udHJpYnV0b3IgY29udHJpYnV0b3JUeXBlPSJFZGl0b3IiPg0KCQkJPGNvbnRyaWJ1dG9yTmFtZT5FbGxpb3R0LCBKYW5uZWFuPC9jb250cmlidXRvck5hbWU+DQoJCQk8Z2l2ZW5OYW1lPkphbm5lYW48L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPkVsbGlvdDwvZmFtaWx5TmFtZT4NCgkJCTxhZmZpbGlhdGlvbj5ET0UvT1NUSTwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPkdlbmF0LCBCZXJyaXQ8L2NvbnRyaWJ1dG9yTmFtZT4NCgkJCTxnaXZlbk5hbWU+QmVycml0PC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5HZW5hdDwvZmFtaWx5TmFtZT4NCgkJCTxhZmZpbGlhdGlvbj5USUI8L2FmZmlsaWF0aW9uPg0KCQk8L2NvbnRyaWJ1dG9yPg0KCQk8Y29udHJpYnV0b3IgY29udHJpYnV0b3JUeXBlPSJFZGl0b3IiPg0KCQkJPGNvbnRyaWJ1dG9yTmFtZT5IYXJ6ZW5ldHRlciwgS2Fyb2xpbmU8L2NvbnRyaWJ1dG9yTmFtZT4NCgkJCTxnaXZlbk5hbWU+S2Fyb2xpbmU8L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPkhhcnplbmV0dGVyPC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPkdFU0lTPC9hZmZpbGlhdGlvbj4NCgkJPC9jb250cmlidXRvcj4NCgkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iRWRpdG9yIj4NCgkJCTxjb250cmlidXRvck5hbWU+SGlyc2NobWFubiwgQmFyYmFyYTwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5CYXJiYXJhPC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5IaXJzY2htYW5uPC9mYW1pbHlOYW1lPg0KCQkJPG5hbWVJZGVudGlmaWVyIG5hbWVJZGVudGlmaWVyU2NoZW1lPSJPUkNJRCIgc2NoZW1lVVJJPSJodHRwOi8vb3JjaWQub3JnIj4wMDAwLTAwMDMtMDI4OS0wMzQ1PC9uYW1lSWRlbnRpZmllcj4NCgkJCTxhZmZpbGlhdGlvbj5FVEggWsO8cmljaDwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPkpha29ic3NvbiwgU3RlZmFuPC9jb250cmlidXRvck5hbWU+DQoJCQk8Z2l2ZW5OYW1lPlN0ZWZhbjwvZ2l2ZW5OYW1lPg0KCQkJPGZhbWlseU5hbWU+SmFrb2Jzc29uPC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPlNORDwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPk1haWxsb3V4LCBKZWFu4oCQWXZlczwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5KZWFuLVl2ZXM8L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPk1haWxsb3V4PC9mYW1pbHlOYW1lPg0KCQkJPGFmZmlsaWF0aW9uPk5SQy9DSVNUSTwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCTxjb250cmlidXRvciBjb250cmlidXRvclR5cGU9IkVkaXRvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPk5ld2JvbGQsIEVsaXphYmV0aDwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5FbGl6YWJldGg8L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPk5ld2JvbGQ8L2ZhbWlseU5hbWU+DQoJCQk8bmFtZUlkZW50aWZpZXIgbmFtZUlkZW50aWZpZXJTY2hlbWU9Ik9SQ0lEIiBzY2hlbWVVUkk9Imh0dHA6Ly9vcmNpZC5vcmciPjAwMDAtMDAwMi04MjU1LTkwMTM8L25hbWVJZGVudGlmaWVyPg0KCQkJPGFmZmlsaWF0aW9uPkJyaXRpc2ggTGlicmFyeTwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJCQkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iRWRpdG9yIj4NCgkJCTxjb250cmlidXRvck5hbWU+TmllbHNlbiwgTGFycyBIb2xtIDwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5MYXJzIEhvbG08L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPk5pZWxzZW48L2ZhbWlseU5hbWU+DQoJCQk8bmFtZUlkZW50aWZpZXIgbmFtZUlkZW50aWZpZXJTY2hlbWU9Ik9SQ0lEIiBzY2hlbWVVUkk9Imh0dHA6Ly9vcmNpZC5vcmciPjAwMDAtMDAwMS04MTM1LTM0ODk8L25hbWVJZGVudGlmaWVyPg0KCQkJPGFmZmlsaWF0aW9uPkNFUk48L2FmZmlsaWF0aW9uPg0KCQk8L2NvbnRyaWJ1dG9yPg0KCQk8Y29udHJpYnV0b3IgY29udHJpYnV0b3JUeXBlPSJFZGl0b3IiPg0KCQkJPGNvbnRyaWJ1dG9yTmFtZT5ZYWhpYSwgTW9oYW1lZDwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5Nb2hhbWVkPC9naXZlbk5hbWU+DQoJCQk8ZmFtaWx5TmFtZT5ZYWhpYTwvZmFtaWx5TmFtZT4NCgkJCTxhZmZpbGlhdGlvbj5JbmlzdC1DTlJTPC9hZmZpbGlhdGlvbj4NCgkJPC9jb250cmlidXRvcj4NCgkJPGNvbnRyaWJ1dG9yIGNvbnRyaWJ1dG9yVHlwZT0iU3VwZXJ2aXNvciI+DQoJCQk8Y29udHJpYnV0b3JOYW1lPlppZWRvcm4sIEZyYXVrZTwvY29udHJpYnV0b3JOYW1lPg0KCQkJPGdpdmVuTmFtZT5GcmF1a2U8L2dpdmVuTmFtZT4NCgkJCTxmYW1pbHlOYW1lPlppZWRvcm48L2ZhbWlseU5hbWU+DQoJCQk8bmFtZUlkZW50aWZpZXIgbmFtZUlkZW50aWZpZXJTY2hlbWU9Ik9SQ0lEIiBzY2hlbWVVUkk9Imh0dHA6Ly9vcmNpZC5vcmciPjAwMDAtMDAwMi0xMTQzLTc4MVg8L25hbWVJZGVudGlmaWVyPg0KCQkJPGFmZmlsaWF0aW9uPlRJQjwvYWZmaWxpYXRpb24+DQoJCTwvY29udHJpYnV0b3I+DQoJPC9jb250cmlidXRvcnM+DQoJPGxhbmd1YWdlPmVuZzwvbGFuZ3VhZ2U+DQoJPHJlc291cmNlVHlwZSByZXNvdXJjZVR5cGVHZW5lcmFsPSJUZXh0Ij5Eb2N1bWVudGF0aW9uPC9yZXNvdXJjZVR5cGU+DQoJPHJlbGF0ZWRJZGVudGlmaWVycz4NCgkJPHJlbGF0ZWRJZGVudGlmaWVyIHJlbGF0ZWRJZGVudGlmaWVyVHlwZT0iRE9JIiByZWxhdGlvblR5cGU9IkRvY3VtZW50cyI+MTAuNTQzOC8wMDEzPC9yZWxhdGVkSWRlbnRpZmllcj4NCgkJPHJlbGF0ZWRJZGVudGlmaWVyIHJlbGF0ZWRJZGVudGlmaWVyVHlwZT0iRE9JIiByZWxhdGlvblR5cGU9IklzTmV3VmVyc2lvbk9mIj4xMC41NDM4LzAwMTA8L3JlbGF0ZWRJZGVudGlmaWVyPg0KCTwvcmVsYXRlZElkZW50aWZpZXJzPg0KCTxzaXplcz4NCgkJPHNpemU+NDUgcGFnZXM8L3NpemU+DQoJPC9zaXplcz4NCgk8Zm9ybWF0cz4NCgkJPGZvcm1hdD5hcHBsaWNhdGlvbi9wZGY8L2Zvcm1hdD4NCgk8L2Zvcm1hdHM+DQoJPHZlcnNpb24+NC4wPC92ZXJzaW9uPg0KCTxkZXNjcmlwdGlvbnM+DQoJCTxkZXNjcmlwdGlvbiBkZXNjcmlwdGlvblR5cGU9IlRhYmxlT2ZDb250ZW50cyI+MSBJbnRyb2R1Y3Rpb248YnIvPg0KMS4xIFRoZSBEYXRhQ2l0ZSBDb25zb3J0aXVtPGJyLz4NCjEuMiBEYXRhQ2l0ZSBDb21tdW5pdHkgUGFydGljaXBhdGlvbjxici8+DQoxLjMgVGhlIE1ldGFkYXRhIFNjaGVtYTxici8+DQoxLjQgVmVyc2lvbiA0LjAgVXBkYXRlPGJyLz4NCjIgRGF0YUNpdGUgTWV0YWRhdGEgUHJvcGVydGllczxici8+DQoyLjEgT3ZlcnZpZXc8YnIvPg0KMi4yIENpdGF0aW9uPGJyLz4NCjIuMyBEYXRhQ2l0ZSBQcm9wZXJ0aWVzPGJyLz4NCjMgWE1MIEV4YW1wbGU8YnIvPg0KNCBYTUwgU2NoZW1hPGJyLz4NCjUgT3RoZXIgRGF0YUNpdGUgU2VydmljZXM8YnIvPg0KQXBwZW5kaWNlczxici8+DQpBcHBlbmRpeCAxOiBDb250cm9sbGVkIExpc3QgRGVmaW5pdGlvbnM8YnIvPg0KQXBwZW5kaXggMjogRWFybGllciBWZXJzaW9uIFVwZGF0ZSBOb3RlczwvZGVzY3JpcHRpb24+DQoJPC9kZXNjcmlwdGlvbnM+DQo8L3Jlc291cmNlPg0K",
      "url": "https://schema.datacite.org/meta/kernel-4.0/index.html",
      "contentUrl": null,
      "metadataVersion": 1,
      "schemaVersion": "http://datacite.org/schema/kernel-4",
      "source": null,
      "isActive": true,
      "state": "findable",
      "reason": null,
      "created": "2016-09-19T21:53:56.000Z",
      "registered": "2016-09-19T21:53:56.000Z",
      "updated": "2018-12-05T04:01:33.000Z"
    },
    "relationships": {
      "client": {
        "data": {
          "id": "datacite.datacite",
          "type": "clients"
        }
      },
      "media": {
        "data": []
      }
    }
  },
  "included": [
    {
      "id": "datacite.datacite",
      "type": "clients",
      "attributes": {
        "name": "DataCite",
        "symbol": "DATACITE.DATACITE",
        "year": 2011,
        "contactName": "DataCite",
        "contactEmail": "support@datacite.org",
        "domains": "*",
        "url": null,
        "created": "2011-12-07T13:43:39.000Z",
        "updated": "2018-08-26T02:35:09.000Z",
        "isActive": true,
        "hasPassword": true
      },
      "relationships": {
        "provider": {
          "data": {
            "id": "datacite",
            "type": "providers"
          }
        },
        "prefixes": {
          "data": [
            {
              "id": "10.5438",
              "type": "prefixes"
            }
          ]
        }
      }
    }
  ]
}
*/



// Sun Oct  3 18:38:51 EDT 2021
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_xref_json   
// g++  -Wall -std=gnu++11 -DTEST_MJM_XREF_JSON -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_xref_json.h  -lpthread -lreadline

template <class Tr>
class mjm_xref_json 
{
 typedef mjm_xref_json Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;

typedef mjm_pawnoff<Tr> Hand;
typedef typename Hand::blob Blob;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef Ragged::Word Word;
typedef mjm_ragged_cursor RaggedCursor;



typedef mjm_wovdb<Tr,StrTy> Tdb;


typedef mjm_hier_two<Tr> HierUtil;


public:
mjm_xref_json() {}
~mjm_xref_json() {}
// Zenodo
StrTy datacite_xform( const StrTy & fn, const IdxTy  flags) const
{
Ragged r;
IdxTy rc=load_ragged(r,fn,0);
std::map<StrTy,StrTy> m;
IdxTy rca=AssembleDataCite(m,r,0);
StrTy s;
IdxTy rcf=format(s,m,0);
MM_ERR(MMPR(fn)<<MMPR4(s,m.size(),r.size(),rc))
return s; 
} // datacite_xform 
StrTy crossref_xform( const StrTy & fn, const IdxTy  flags) const
{
Ragged r;
IdxTy rc=load_ragged(r,fn,0);
std::map<StrTy,StrTy> m;
IdxTy rca=assemble(m,r,0);
StrTy s;
IdxTy rcf=format(s,m,0);
MM_ERR(MMPR(fn)<<MMPR4(s,m.size(),r.size(),rc))
return s; 
} // crossref_xform


IdxTy load_ragged_orig(Ragged & r, const StrTy & fn, const IdxTy  flags) const
{
const StrTy cmd=" sed -e 's/\\\\\\\"/mjm_bs_bug/g' | mjsonu.out -cmd \"sax-json - 1\"  | sed -e 's/mjm_bs_bug/\\\\\\\"/g'  ";
Blob out,err,data;
data.load(fn);
IdxTy c=m_hand.fileio(out,err,data,cmd,2);
//MM_ERR(MMPR(StrTy(err)));
//MM_ERR(MMPR(StrTy(out)));
///Ragged r;
r.kvp_splitter(0);
Ss ss;
ss<<StrTy(out);
r.load(ss); // this should now be parsed with value in the end of line field 

MM_ERR(r.dump())
return 0; 
}



IdxTy load_ragged(Ragged & r, const StrTy & fn, const IdxTy  flags) const
{
//MM_MSG(" using integrated json for xref may not work doh ")
typedef mjm_mjsonu JsonTy;
JsonTy json;
std::ifstream is(fn);
auto pr=json.dirty_parse(is,0);
//Ragged & rn=pr.res;
 r=pr.res;
MM_ERR(r.dump())
return 0; 
}



template <class Tm > IdxTy assemble(Tm & m, const Ragged & r,  const IdxTy  flags) const
{
return Assemble(m,r,flags); 
}
template <class Tm > IdxTy format(StrTy & s, Tm & m,  const IdxTy  flags) const
{
return Format(s,m,flags);
}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);
/*
mjm_xref_json.h74 0|./mjsonu.h280|GLOBAL|1||1|status|ok
1|./mjsonu.h280|GLOBAL|1||2|message-type|work
2|./mjsonu.h280|GLOBAL|1||3|message-version|1.0.0
3|./mjsonu.h280|GLOBAL|1||4|message|1|indexed|1|date-parts|1||1|date-parts|2021
4|./mjsonu.h280|GLOBAL|1||4|message|1|indexed|1|date-parts|1||2|date-parts|10
5|./mjsonu.h280|GLOBAL|1||4|message|1|indexed|1|date-parts|1||3|date-parts|3
6|./mjsonu.h280|GLOBAL|1||4|message|1|indexed|2|date-time|2021-10-03T16:45:09Z
7|./mjsonu.h280|GLOBAL|1||4|message|1|indexed|3|timestamp|1633279509242
8|./mjsonu.h280|GLOBAL|1||4|message|2|reference-count|6
9|./mjsonu.h280|GLOBAL|1||4|message|3|publisher|Canadian Science Publishing
10|./mjsonu.h280|GLOBAL|1||4|message|4|issue|8
11|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||1|start|1|date-parts|1||1|date-parts|1959
12|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||1|start|1|date-parts|1||2|date-parts|8
13|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||1|start|1|date-parts|1||3|date-parts|1
14|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||1|start|2|date-time|1959-08-01T00:00:00Z
15|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||1|start|3|timestamp|-328838400000
16|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||2|content-version|tdm
17|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||3|delay-in-days|0
18|./mjsonu.h280|GLOBAL|1||4|message|5|license|1||4|URL|http://www.nrcresearchpress.com/page/about/CorporateTextAndDataMining


*/

class _msg_state
{
public:
_msg_state():mpos(~0) {} 
StrTy dump() const
{
Ss ss;
ss<<MMPR4(mpos,type,auno,status);
return ss.str();}
IdxTy mpos;
StrTy type,auno,status;

}; 

template <class Tm > IdxTy AddAuthors(Tm & m, const Ragged & r, const IdxTy flags) const 
{ 
HierUtil hu;
//typedef std::map<StrTy, std::map<StrTy, StrTy> > Hm;
typedef typename HierUtil::hier_map Hm;
Hm hm;
//hu.organize(hm,StrTy(b),"author",0);
hu.organize(hm,r,"author",0);
StrTy au=m["author"]; 
StrTy af=m["affiliation"];
/*
30|GLOBAL|0||0|message|1|author|20||21|given|=|M rio Felipe
31|GLOBAL|0||0|message|1|author|20||21|family|=|Balaro
32|GLOBAL|0||0|message|1|author|20||21|sequence|=|first
33|GLOBAL|0||0|message|1|author|20||21|affiliation|22||23|name|=|Universidade Federal Fluminense
34|GLOBAL|0||0|message|1|author|20||24|given|=|Fernanda Gon


*/
MM_LOOP(ii,hm)
{
// make a copy and use non const  [] access
auto   mm=(*ii).second; // a map of kv-pairs
if ( mm["name"].length())
{
//MM_ERR(MMPR3(af,nm["name"],nm["name"].length()))
if (af.length()) af+= " and ";
//af+=mm["affiliation"];
af+=mm["name"];
} 
const StrTy fn=mm["family"];
const StrTy gn=mm["given"];
const IdxTy lena=fn.length()+gn.length();
MM_ERR(MMPR3(lena,fn.length(),gn.length()))
if (lena)
{
if (au.length()) au+= " and ";
au+=gn+" "+fn;
} // lena
Ss ss; 
MM_LOOP(jj,mm) { ss<<(*jj).first<<" "; }
MM_ERR(" keys "<<MMPR2((*ii).first,ss.str()))

} // ii 
m["affiliation"]=af;
m["author"]=au;
return 0; 
} // Add authors
///////////////////////////////////////////////////////
/*

Data Cite:

0|GLOBAL|0||0|data|1|id|=|10.13140/rg.2.2.24534.27207
1|GLOBAL|0||0|data|1|type|=|dois
2|GLOBAL|0||0|data|1|attributes|2|doi|=|10.13140/rg.2.2.24534.27207
3|GLOBAL|0||0|data|1|attributes|2|prefix|=|10.13140
4|GLOBAL|0||0|data|1|attributes|2|suffix|=|rg.2.2.24534.27207
5|GLOBAL|0||0|data|1|attributes|2|creators|5||6|name|=|Lundeen, James
6|GLOBAL|0||0|data|1|attributes|2|creators|5||6|nameType|=|Personal
7|GLOBAL|0||0|data|1|attributes|2|creators|5||6|givenName|=|James
8|GLOBAL|0||0|data|1|attributes|2|creators|5||6|familyName|=|Lundeen
9|GLOBAL|0||0|data|1|attributes|2|titles|9||10|title|=|The Spectrum of COVID-19 Disease Manifestations and Rational Treatment Approaches
10|GLOBAL|0||0|data|1|attributes|2|publisher|=|Unpublished
11|GLOBAL|0||0|data|1|attributes|2|publicationYear|=|2020
12|GLOBAL|0||0|data|1|attributes|2|dates|14||15|date|=|2020
13|GLOBAL|0||0|data|1|attributes|2|dates|14||15|dateType|=|Issued
14|GLOBAL|0||0|data|1|attributes|2|language|=|en
15|GLOBAL|0||0|data|1|attributes|2|types|16|ris|=|RPRT
16|GLOBAL|0||0|data|1|attributes|2|types|16|bibtex|=|article
17|GLOBAL|0||0|data|1|attributes|2|types|16|citeproc|=|article-journal
18|GLOBAL|0||0|data|1|attributes|2|types|16|schemaOrg|=|ScholarlyArticle
19|GLOBAL|0||0|data|1|attributes|2|types|16|resourceType|=|Preprint
20|GLOBAL|0||0|data|1|attributes|2|types|16|resourceTypeGeneral|=|Text
21|GLOBAL|0||0|data|1|attributes|2|version|=|
23|GLOBAL|0||0|data|1|attributes|2|url|=|http://rgdoi.net/10.13140/RG.2.2.24534.27207
24|GLOBAL|0||0|data|1|attributes|2|contentUrl|=|
25|GLOBAL|0||0|data|1|attributes|2|metadataVersion|=|0
26|GLOBAL|0||0|data|1|attributes|2|schemaVersion|=|http://datacite.org/schema/kernel-4
27|GLOBAL|0||0|data|1|attributes|2|source|=|mds
28|GLOBAL|0||0|data|1|attributes|2|isActive|=|true
29|GLOBAL|0||0|data|1|attributes|2|state|=|findable
30|GLOBAL|0||0|data|1|attributes|2|reason|=|
31|GLOBAL|0||0|data|1|attributes|2|viewCount|=|0
32|GLOBAL|0||0|data|1|attributes|2|downloadCount|=|0
33|GLOBAL|0||0|data|1|attributes|2|referenceCount|=|0
34|GLOBAL|0||0|data|1|attributes|2|citationCount|=|0
35|GLOBAL|0||0|data|1|attributes|2|partCount|=|0
36|GLOBAL|0||0|data|1|attributes|2|partOfCount|=|0
37|GLOBAL|0||0|data|1|attributes|2|versionCount|=|0
38|GLOBAL|0||0|data|1|attributes|2|versionOfCount|=|0
39|GLOBAL|0||0|data|1|attributes|2|created|=|2020-06-02T04:32:31.000Z
40|GLOBAL|0||0|data|1|attributes|2|registered|=|2020-06-02T04:32:33.000Z
41|GLOBAL|0||0|data|1|attributes|2|published|=|2020
42|GLOBAL|0||0|data|1|attributes|2|updated|=|2020-08-01T20:11:25.000Z
43|GLOBAL|0||0|data|1|relationships|27|client|28|data|29|id|=|rg.rg
44|GLOBAL|0||0|data|1|relationships|27|client|28|data|29|type|=|clients
45|GLOBAL|0||0|data|1|relationships|27|provider|30|data|31|id|=|rg
46|GLOBAL|0||0|data|1|relationships|27|provider|30|data|31|type|=|providers
47|GLOBAL|0||0|data|1|relationships|27|media|32|data|33|id|=|10.13140/rg.2.2.24534.27207
48|GLOBAL|0||0|data|1|relationships|27|media|32|data|33|type|=|media


*/
template <class Tm> void AddAnd(Tm& m,const StrTy & s, const StrTy & v) const
{
auto ii=m.find(s);
if (ii==m.end()){ m[s]=v; return ; } 
if ((*ii).second.length()) m[s]=m[s]+" and "+v;
else m[s]=v;
} // AddAnd
template <class Tm > IdxTy AssembleDataCite(Tm & m, const Ragged & _r,  const IdxTy  flags) const
{
IdxTy rc=0;
m["type-name"]="article";
m["name-name"]="xxx";

const Ragged&  r=_r;
IdxTy attrpos=~0;
MM_SZ_LOOP(i,r,szr)
{
const Line & l=r[i];
const IdxTy len=l.size();
if (len<2) continue; // no kvp 
//MM_ERR(MMPR(l[len-1]))
const IdxTy off=(l[len-2]=="=")?1:0;
if (len<(2+off)) continue; // no kvp 
const StrTy & v=l[len-1];
const StrTy & k1=l[len-2-off];
if (attrpos!=~0){ if ( len<=(attrpos+3) ) attrpos=~0; } 
if (attrpos!=~0){ if ( l[attrpos]!="attributes") attrpos=~0; } 
if (attrpos==~0)
{
// this may be dangeous though... 
MM_SZ_LOOP(j,l,llen) { if (l[j]=="attributes"){ attrpos=j; break; }}
} // attrpos
if (attrpos==~0) continue;
if (attrpos!=~0) if ( len<=(attrpos+3) ) { attrpos=~0; continue; }  
if (l[attrpos+2]=="creators") 
{
if (k1=="name") AddAnd(m,"authors",v); 
continue;
}// creators
AddAnd(m,k1,v);
//if (k1=="bibtex") m["type-name"]=v;
} // i 
return rc;
}// AssembleDataCite

///////////////////////////////////////////////////////

typedef _msg_state msg_state;
template <class Tm > IdxTy Assemble(Tm & m, const Ragged & _r,  const IdxTy  flags) const
{
// TODO FIX THIS KLUGE LOL
// now there is an equals doh 
// 1|GLOBAL|0||0|message-type|=|work
Ragged r;
MM_SZ_LOOP(i,_r,szr)
{
const Ragged::Line & line=_r[i];
const IdxTy len=line.size();
if (len>3) if ( line[len-2]=="=")
{
Ragged::Line rl;
for(IdxTy j=0; j<(len-2); ++j) rl.push_back(line[j]);
rl.push_back(line[len-1]);
r.add(rl);
continue;
}
r.add(_r[i]);

} // i 

Ragged ra;
msg_state ms;
const IdxTy sz=r.size();
for(IdxTy i=0; i<sz; ++i)
{
const Line & l=r[i];
const IdxTy len=l.size();
if (len<2) continue; // no kvp 
MM_ERR(MMPR(l[len-1]))
if (ms.mpos==~0U)
{
if (l[len-2]=="message-type")
{
ms.mpos=len-2;
ms.type=l[len-1];
continue;
} // msg-type
if (l[len-2]=="status")
{
ms.status=l[len-1];
continue;

} // status

} // mpos==-1
if (len<=ms.mpos) continue;
MM_ERR(MMPR(l[ms.mpos]))
if (l[ms.mpos]!="message") continue;
//9|./mjsonu.h280|GLOBAL|1||4|message|3|publisher|Canadian Science Publishing
if (len<=( ms.mpos+3)) continue;
const StrTy & k1=l[ms.mpos+2];
// some fckers put in cr and lf fck 
//const StrTy & v=l[len-1];
StrTy _v=l[len-1];

const IdxTy lv=_v.length();
const char *  cp=_v.c_str();
char * vc=new char[lv+2];
IdxTy p=0;
IdxTy pout=0;
for(; p<lv; ++p)
{
if (cp[p]=='\\')
{
if (cp[p+1]=='x') p+=3; continue;
}
{ vc[pout]=_v.c_str()[p]; ++pout; } 
} 
vc[pout]=0;

const StrTy v= StrTy(vc);
delete[] vc;
// TODO most of these are the same, can use a map... 
if (k1=="publisher"){  m[k1]=v; continue; } 
// problem parsing.. 
// https://efsa.onlinelibrary.wiley.com/doi/pdf/10.2903/sp.efsa.2016.EN-1027
if (k1=="reference"){ m[k1]="deleted";   continue; } 
if (k1=="references"){ m[k1]="deleted";  continue; } 
/*
97|GLOBAL|1||2|message|3|reference|69||70|key|=|10.2903/sp.efsa.2016.EN-1027_ref1
98|GLOBAL|1||2|message|3|reference|69||70|unstructured|=|AGES , 2013‐ Austrian Agency for Health & Food Safety, Institute for Plant Protection Products. «External scientific report. Case studies for the application of the Guidance of EFSA on Submission of scientific peer‐reviewed open literature for the approval of pesticide active substances under Regulation (EC) No 1107/2009, using substances for which dossiers are submitted under Regulation (EU) No 1141/20101." Report and appendixes available on the internet (http://www.efsa.europa.eu/en/supporting/doc/511e.pdf).
99|GLOBAL|1||2|message|3|reference|69||71|issue|=|3
100|GLOBAL|1||2|message|3|reference|69||71|key|=|10.2903/sp.efsa.2016.EN-1027_ref2
101|GLOBAL|1||2|message|3|reference|69||71|doi-asserted-by|=|crossref
102|GLOBAL|1||2|message|3|reference|69||71|first-page|=|185
103|GLOBAL|1||2|message|3|reference|69||71|DOI|=|10.1016/j.tox.2006.07.022
104|GLOBAL|1||2|message|3|reference|69||71|article-title|=|A dose‐response study following in utero and lactational exposure to di‐(2‐ethylhexyl)‐phthalate (DEHP): non‐monotonic dose‐response and low dose effects on rat brain aromatase activity
105|GLOBAL|1||2|message|3|reference|69||71|volume|=|227
106|GLOBAL|1||2|message|3|reference|69||71|author|=|Andrade
107|GLOBAL|1||2|message|3|reference|69||71|year|=|2006
108|GLOBAL|1||2|message|3|reference|69||71|journal-title|=|Toxicology
109|GLOBAL|1||2|message|3|reference|69||72|issue|=|1


*/
if (k1=="issue"){  m[k1]=v; continue; } 
if (k1=="volume"){  m[k1]=v; continue; } 
if (k1=="title"){  m[k1]=v; continue; } 
if (k1=="page"){  m[k1]=v; continue; } 
if (k1=="abstract"){  m[k1]=v; continue; } 
if (k1=="DOI"){  m["doi"]=v; continue; } 
if (k1=="source"){  m["bib-source"]=v; continue; } 
if (k1=="subject"){  m["subject"]=v; continue; } 
if (k1=="ISSN"){  m["issn"]=v; continue; } 
if (k1=="URL"){  m["url"]=v; continue; } 
if (k1=="short-container-title"){  m["abbrvjrnl"]=v; continue; } 
if (k1=="container-title"){  m["journal"]=v; continue; }
// TODO dates about the same too 
// this is messing up when there is a journal-issue date... 
if (k1=="published"){ 
//StrTy pkk="date"; // this confuses publication with journal date...
StrTy pkk="publication-date"; // this confuses publication with journal date...
if (m[pkk].length()!=0) m[pkk]+="-";  
if (v.length()==1) m[pkk]+=StrTy("0")+v; else m[pkk]+=v; 
continue; 
} 
if (k1=="type"){  m["type-name"]=v; 
  m["type"]=v; 
if (v=="journal-article") m["type-name"]="article"; 
if (v=="posted-content") m["type-name"]="article"; 
continue; }
// need to organize not parse 
if (!false) if (k1=="author"){  ra.add(l); continue;  }
#if 0
if (false) if (k1=="author"){  
/*
38|./mjsonu.h280|GLOBAL|1||4|message|19|author|1||1|given|E. G.
39|./mjsonu.h280|GLOBAL|1||4|message|19|author|1||2|family|Bligh
40|./mjsonu.h280|GLOBAL|1||4|message|19|author|1||3|sequence|first
41|./mjsonu.h280|GLOBAL|1||4|message|19|author|2||1|given|W. J.
42|./mjsonu.h280|GLOBAL|1||4|message|19|author|2||2|family|Dyer
43|./mjsonu.h280|GLOBAL|1||4|message|19|author|2||3|sequence|addi

onu.h280|GLOBAL|1||4|message|16|author|9||3|sequence|additional
75|./mjsonu.h280|GLOBAL|1||4|message|16|author|9||4|affiliation|1||1|name|Department of Functional Brain Imaging, National Institute of Radiological Sciences, National Institutes for Quantum and Radiological Sciences and Technology, Chiba 263-8555, Japan.
76|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||1|ORCID|http://orcid.org/0000-0003-2239-1754
77|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||2|authenticated-orcid|true
78|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||3|given|Maiko
79|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||4|family|Ono
80|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||5|sequence|additional
81|./mjsonu.h280|GLOBAL|1||4|message|16|author|10||6|affiliation|1||1|name|Department of Functional Brain Imaging, National Institute of Radiological Sciences, National Institutes for Quantum and Radiological Sciences and Technology, Chiba 263-8555, Japan.
82|./mjsonu.h280|GLOBAL|1||4|message|16|author|11||1|given|Jun
83|./mjsonu.h280|GLOBAL|1||4|message|16|author|11||2|family|Maeda

// refs in wrong pos so automatically dropped... 
cat aaaa| grep GLOBAL | grep author | grep refer | tail -n 2
546|./mjsonu.h280|GLOBAL|1||4|message|18|reference|50||7|author|Sahara N.
556|./mjsonu.h280|GLOBAL|1||4|message|18|reference|51||7|author|Maeda J.


const StrTy & k1=l[ms.mpos+2];

*/
//const StrTy & k1=l[ms.mpos+2];
const IdxTy pmod=ms.mpos+2+4;
bool aff=(pmod<len);

#if 0
/*
this needs to use the number before the autho to group the fking names 
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  27 given = Chien-Che
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  27 family = Chang
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  27 sequence = additional
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  27 affiliation 28  29 name = Depa
rtment and Graduate Institute of Clinical Laboratory Sciences and Medical Biotec
hnology, College of Medicine, National Taiwan University, Taipei, Taiwan, Republ
ic of China
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  30 given = Yu-Han
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  30 family = Yuan
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  30 sequence = additional
./mjsonu.h312 GLOBAL 0  0 message 1 author 23  30 affiliation 31  32 name = Depa
rtment and Graduate Institute of Clinical Laboratory Sciences and Medical Biotec
hnology, College of Medicine, National Taiwan University, Taipei, Taiwan, Republ
ic of China
*/


42|GLOBAL|0||0|message|1|author|20||30|given|=|Isabel
43|GLOBAL|0||0|message|1|author|20||30|family|=|Cosentino
44|GLOBAL|0||0|message|1|author|20||30|sequence|=|additional
45|GLOBAL|0||0|message|1|author|20||30|affiliation|31||32|name|=|Universidade Federal Fluminense
ZZ


#endif



const StrTy ty=l[len-2];
const StrTy no=l[len-5];
MM_ERR(MMPR3(ty,no,v))
//MM_ERR("  affiliation "<<MMPR4(ty,l[pmod],aff
StrTy newt;
const IdxTy la=m["author"].length();
if (ty=="sequence" ) continue;
if(  (ty=="affiliation" )
 ||((aff)&&(l[pmod]=="affiliation" ))
 ||(ty=="authenticated-orcid" )
 ||(ty=="ORCID" ))
{
StrTy tyn=ty;
 if (aff) if (l[pmod]=="affiliation" ) tyn=l[pmod];
Ss ss;
ss<<"("<<no<<")"<<v; 
m[tyn]+= ss.str();
 continue;
} // affiliation


if (la)  m["author"]+=" ";
if (la) if( ms.auno!=no) m["author"]+="and ";

if (la) if( ms.auno==no)if (ty=="given")  m["author"]+=", ";
// above just puts in the conjection punctuation, 
if ((ty=="given")||(ty=="family")) m["author"]+=v; // newt;
else m["author"]+=" unknown "+v;

ms.auno=no;
continue; }  // author, gated out 
#endif


if (k1=="indexed"){  
// ./mjsonu.h280 GLOBAL 1  4 message 1 indexed 2 date-time = 2021-10-03T16:45:09Z
const StrTy ty=l[len-2];
if (ty=="date-time") m["date-indexed"]=v; 
continue; 
} // indexed

if (k1=="created"){  
// ./mjsonu.h280 GLOBAL 1  4 message 1 indexed 2 date-time = 2021-10-03T16:45:09Z
const StrTy ty=l[len-2];
if (ty=="date-time") m["date-created"]=v; 
continue; 
} // indexed
//continue;
// } // mpos==0
// a left-over 

/*
74|GLOBAL|0||0|message|1|member|=|9829
75|GLOBAL|0||0|message|1|deposited|56|date-parts|57||58|date-parts|=|2020
76|GLOBAL|0||0|message|1|deposited|56|date-parts|57||58|date-parts|=|6
77|GLOBAL|0||0|message|1|deposited|56|date-parts|57||58|date-parts|=|1
78|GLOBAL|0||0|message|1|deposited|56|date-time|=|2020-06-01T18:40:05Z
79|GLOBAL|0||0|message|1|deposited|56|timestamp|=|1591036805000
80|GLOBAL|0||0|message|1|score|=|1
81|GLOBAL|0||0|message|1|issued|61|date-parts|62||63|date-parts|=|2020
82|GLOBAL|0||0|message|1|issued|61|date-parts|62||63|date-parts|=|6
83|GLOBAL|0||0|message|1|issued|61|date-parts|62||63|date-parts|=|1
84|GLOBAL|0||0|message|1|references-count|=|0
85|GLOBAL|0||0|message|1|URL|=|http://dx.doi.org/10.22541/au.159103680.00306295
86|GLOBAL|0||0|message|1|published|65|date-parts|66||67|date-parts|=|2020
87|GLOBAL|0||0|message|1|published|65|date-parts|66||67|date-parts|=|6
88|GLOBAL|0||0|message|1|published|65|date-parts|66||67|date-parts|=|1
89|GLOBAL|0||0|message|1|subtype|=|preprint


*/

// now we have a key and a type 
//if (k1=="created"){  
const StrTy ty=l[len-2];
if ( ty=="date-parts")
{
StrTy kf=StrTy("date-")+k1;
if (m[kf].length()) m[kf]+="-";
if (v.length()==1) m[kf]+="0";
m[kf]+=v;
continue;
}
if ( ty=="date-time")
{
StrTy kf=StrTy("date-")+k1;
m[kf]+=v;
continue;
}
if (m[k1].length()) m[k1]+=", ";
m[k1]+=v;
} // i 
AddAuthors(m,ra,0);
m["name-name"]="xxx";
return 0; 
} // Assemble 
template <class Tm > IdxTy Format(StrTy & s, Tm & m,  const IdxTy  flags) const
{
Ss ss;
ss<<"@"<<m["type-name"]<<"{"<<m["name-name"]; // <<","<<CRLF;
MM_LOOP(ii,m)
{
const StrTy & key=(*ii).first;
if (key=="type-name") continue;
if (key=="name-name") continue;

ss<<","<<CRLF;
const StrTy & v=(*ii).second;
ss<<key<<" = {"<<v<<"} ";

}
ss<<CRLF<<"}"<<CRLF;

s=ss.str();
return 0; 
}// Format 

// MEMBERS
mutable Hand m_hand;


}; // mjm_xref_json

//////////////////////////////////////////////

template <class Tr>
class mjm_xref_json_map : public std::map<typename Tr::StrTy, mjm_xref_json< Tr > >  
{
 typedef mjm_xref_json_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_xref_json< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_xref_json_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }

private:

void Init()
{
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;


}
return ss.str();
// return Dump(flags); 

}




private:

}; // mjm_xref_json_map




////////////////////////////////////////////
#ifdef  TEST_MJM_XREF_JSON
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_xref_json <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_XREF_JSON "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

typedef mjm_ragged_table Ragged;
int main(int argc,char **args)
{
about();
typedef mjm_xref_json<Tr>  Myt;
//Myt x(argc,args);
Myt x;

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
else if (cmd=="load") {
// x.load(li.words(),1); 
Ragged r; 
IdxTy rc=x.load_ragged(r,cip.p1,0);
std::map<StrTy,StrTy> m;
IdxTy rca=x.assemble(m,r,0);
StrTy s; 
IdxTy rcf=x.format(s,m,0);
MM_ERR(MMPR4(s,m.size(),r.size(),rc))

//IdxTy load_ragged(Ragged & r, const StrTy & fn, const IdxTy  flags) const
}
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_XREF_JSON_H__ 
