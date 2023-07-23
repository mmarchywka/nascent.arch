/*========================================================================*\

Copyright (c) 1990-2013  Paul Vojta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
PAUL VOJTA OR ANY OTHER AUTHOR OF THIS SOFTWARE BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

NOTE:
xdvi is based on prior work, as noted in the modification history
in xdvi.c.

\*========================================================================*/

/*
#include "xdvi-config.h"
#include "xdvi.h"
#include "dvi-init.h"
//#include "dvi-draw.h"
#include "util.h"
#include "x_util.h"
#include "exit-handlers.h"
#include "mime.h"
#include "pagesel.h"
#include "special.h"
#include "hypertex.h"
#include "kpathsea/c-fopen.h"
#include "kpathsea/c-stat.h"
#include "kpathsea/magstep.h"
#include "kpathsea/tex-glyph.h"
#include "dvi.h"
#include "string-utils.h"
#include "browser.h"
#include "sfSelFile.h"
#include "xm_toolbar.h"
#include "pagehist.h"
#include "message-window.h"
#include "search-internal.h"
#include "statusline.h"
#include "events.h"
#include "font-open.h"
*/
/////////////////////////////////////////////////////////
/*

#include "xdvi-config.h"
#include "xdvi.h"

#include <stdarg.h>

#include <stdlib.h>
#include <ctype.h>

#include <setjmp.h>

#define USE_HASH

#include <ctype.h>
#include "kpathsea/c-fopen.h"
#include "kpathsea/c-stat.h"
#include "kpathsea/magstep.h"
#include "kpathsea/tex-file.h"

#include <string.h>

#include "dvi.h"
#include "string-utils.h"
#include "util.h"
#include "x_util.h"
#include "events.h"
#include "dvi-init.h"
#include "font-open.h"
#include "statusline.h"
#include "hypertex.h"
#include "special.h"
#include "my-snprintf.h"
#include "kpathsea/tex-file.h"
#include "mag.h"
#include "message-window.h"

#include "dvi-draw.h"
#include "search-internal.h"
#include "encodings.h"
#include "pagesel.h"
#include "pagehist.h"

*/



///////////////////////////////////////////////////////////




void mjm_run_dvi_file(const char *filename, void *data);
void mjm_do_pages(void);
extern int mjm_made_it_to_mjm_do_pages;

extern int mjm_failed_to_init;

#include "xdvi-config.h"

#include "c-auto.h"
#include "dvi-init.h"
Boolean mjm_load_dvi_file(
#if !DELAYED_MKTEXPK
	      Boolean load_fonts,
#endif
	      dviErrFlagT *errflag);


