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
// this is c code doh 
//#include "mjm_globals.h"
#include "xdvi_code_pieces.h"
//#include "xdvi_dvi-init.h"
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



extern  int source_show_all;
extern int source_reverse_x, source_reverse_y;
extern FILE *m_dvi_fp;
extern struct stat fstatbuf;
extern home_proc home_action;
extern int file_history_get_page(void);
extern int check_goto_page(int pageno, Boolean insert_into_pagehist);
extern Position m_window_x, m_window_y;
extern Arg arg_xy[2]; 
#define get_xy() XtGetValues(globals.widgets.draw_widget, arg_xy, XtNumber(arg_xy))




/*
xdvi_code_pieces.c:178:19: error: ‘event_freq’ undeclared (first use in this function)
xdvi_code_pieces.c:193:34: error: ‘all_signals’ undeclared (first use in this function)
xdvi_code_pieces.c:207:14: error: ‘sig_flags’ undeclared (first use in this function)
xdvi_code_pieces.c:208:4: error: ‘flags_to_sigproc’ undeclared (first use in this function)
xdvi_code_pieces.c:262:7: error: ‘io_dirty’ undeclared (first use in this function)
xdvi_code_pieces.c:265:11: error: ‘num_fds’ undeclared (first use in this function)
xdvi_code_pieces.c:265:21: error: ‘max_fds’ undeclared (first use in this function)
xdvi_code_pieces.c:266:8: error: ‘fds’ undeclared (first use in this function)
xdvi_code_pieces.c:273:17: error: ‘iorecs’ undeclared (first use in this function)
xdvi_code_pieces.c:398:6: error: ‘resized’ undeclared (first use in this functio


*/


typedef	void	(*signalproc)(void); 
extern VOLATILE int event_freq;
extern sigset_t all_signals;
extern VOLATILE unsigned int sig_flags ;
extern const signalproc flags_to_sigproc[64] ;
    
extern struct xio      *iorecs ; /* head of xio list */

#if HAVE_POLL
extern struct pollfd   *fds ;
extern int             num_fds ;    /* current number of fds */
extern int             max_fds;    /* max allocated number of fds */
extern Boolean         io_dirty; /* need to recompute fds[] array */
#else
static  int             numfds  = 0;
static  fd_set          readfds;
static  fd_set          writefds;
#endif


extern Boolean resized;





/////////////////////////////////////
/*======================================================================*\

Copyright (c) 1990-2016  Paul Vojta and others

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
 *	Since redrawing the screen is (potentially) a slow task, xdvi checks
 *	for incoming events while this is occurring.  It does not register
 *	a work proc that draws and returns every so often, as the toolkit
 *	documentation suggests.  Instead, it checks for events periodically
 *	(or not, if SIGPOLL can be used instead) and processes them in
 *	a subroutine called by the page drawing routine.  This routine (below)
 *	checks to see if anything has happened and processes those events and
 *	signals.  (Or, if it is called when there is no redrawing that needs
 *	to be done, it blocks until something happens.)
 *
 *	Ultimately, the goal is to have this be the only place in xdvi where
 *	blocking occurs.
 *
 *	The argument to this function should be a mask of event types (EV_*)
 *	indicating which event types should cause read_events to return instead
 *	of waiting for more events.  This function will always process all
 *	pending events and signals before returning.
 *	The return value is the value of globals.ev.flags.
 */
///////////////////////////////////////////////////////////////



void
mjm_redraw_page(void)
{
#if COLOR
    const struct rgb *rgbp;
#endif
    TRACE_FILES((stderr, "Redraw page on %p", (void *)globals.dvi_file.bak_fp));

    if (globals.dvi_file.bak_fp == NULL)
        return;

    if (scanned_page < current_page) {
        TRACE_FILES((stderr, "redraw_page: scanned_page = %d, current_page = %d, prescanning %p\n",
                     scanned_page, current_page, (void *)globals.dvi_file.bak_fp));

        prescan(globals.dvi_file.bak_fp);

        if (globals.ev.flags & EV_GE_NEWPAGE) { /* if we need to re-prescan */
            return;
        }
    }

    TRACE_FILES((stderr, "redraw_page: current_page = %d", current_page));
    if (pageinfo_get_window_width(current_page) != globals.page.unshrunk_w
        || pageinfo_get_window_height(current_page) != globals.page.unshrunk_h) {
        TRACE_FILES((stderr, "NEW SIZE: %dx%d",
                     pageinfo_get_window_width(current_page), pageinfo_get_window_height(current_page)));
        init_page();
        reconfig();
    }

    /* We can't call home() without proper unshrunk_page_*, which requires
     * prescan(), which can't be done from within read_events() */

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUG ALERT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       There's some complicated interaction with Postscript specials
       here: if home_action check comes before the gs stuff, psp.drawfile
       might not get initialized correctly, resulting in empty PS figures
       (bounding box instead of figure). This is different in xdvi, due to
       different handling of the home_action stuff, but at the moment I can't
       remember the reason for this ...
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUG ALERT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
#if PS_GS
    if (gs_postpone_prescan) {
        if (!setjmp(globals.ev.canit)) {
            gs_resume_prescan();
        }
        else
            return;
    }
#endif
    if (home_action != NULL) {
        home_action(False);
        home_action = NULL;
        /* This discards the expose event generated by home()
           (1 for each page) */
printf(" need to change - want to call mjm_read_Events lol\n");
        if (read_events(EV_NOWAIT) & EV_GE_NEWPAGE) {
            return;
        }

        can_exposures(&mane);
    }

#if COLOR
    rgbp = &bg_initial;
    if (page_colors.stack != NULL) {
        ASSERT(current_page < (int)page_colors.size, "page_colors.size too small");
        rgbp = &page_colors.stack[current_page].bg;
    }

    /* Set background color */
    if (bg_current == NULL
        || rgbp->r != bg_current->color.r
        || rgbp->g != bg_current->color.g
        || rgbp->b != bg_current->color.b) {
        struct bgrec **bgpp;

        for (bgpp = &bg_head;;) {
            bg_current = *bgpp;
            if (bg_current == NULL) {   /* if bg is not in list */
                bg_current = *bgpp = xmalloc(sizeof *bg_current);
                bg_current->next = NULL;
                bg_current->color = *rgbp;
                bg_current->fg_head = NULL;
                bg_current->pixel_good = False;
                                                                          
               break;
            }
            if (bg_current->color.r == rgbp->r
                && bg_current->color.g == rgbp->g
                && bg_current->color.b == rgbp->b)
                break;
            bgpp = &bg_current->next;
        }
        fg_current = NULL;      /* force change of foreground color */
        /* globals.gc.high is only used in XDrawRectangle, so its background color
           doesn't need to be changed.  */
        if (globals.debug & DBG_DVI)
            printf("Changing background color to %5d %5d %5d\n",
                   bg_current->color.r, bg_current->color.g,
                   bg_current->color.b);

        if (!bg_current->pixel_good) {
            bg_current->pixel = alloc_color(&bg_current->color,
                                            color_data[1].pixel);
            bg_current->pixel_good = True;
        }
        XSetWindowBackground(DISP, mane.win, bg_current->pixel);
#if MOTIF && !FIXED_FLUSHING_PAGING
        XtVaSetValues(XtParent(globals.widgets.draw_widget), XtNbackground, bg_current->pixel, NULL);
#endif
        /*      XSetWindowBackground(DISP, mane.win, bg_current->pixel); */
        /*      XClearWindow(DISP, mane.win); */
#if 0 /* don't recolor the cursor - gives too low contrast on color backgrounds,
         and bad appearance when part of the background is white and cursor mask
         is colored */
        {
            XColor      bg_Color;
            bg_Color.pixel = bg_current->pixel;
            XQueryColor(DISP, G_colormap, &bg_Color);
            XRecolorCursor(DISP, globals.cursor.ready, &globals.cr_color, &bg_Color);
            XRecolorCursor(DISP, globals.cursor.wait, &globals.cr_color, &bg_Color);
        }
#endif
    }
#endif /* COLOR */

    if (!globals.pausing.flag) {
        XClearWindow(DISP, mane.win);
    }

    if (G_backing_store != NotUseful) {
        mane.min_x = mane.min_y = 0;
        mane.max_x = globals.page.w;
        mane.max_y = globals.page.h;
    }
    else {
        get_xy();
        mane.min_x = -m_window_x;
        mane.max_x = -m_window_x + mane.width;
        mane.min_y = -m_window_y;
        mane.max_y = -m_window_y + mane.height;
    }

    /*      update_TOC(); */
    redraw(&mane);
}

////////////////////////////////////////////////////////////////
unsigned int
mjm_read_events(unsigned int ret_mask)
{
static int cnt=0;
    XEvent event;

#if !HAVE_POLL
    if (numfds == 0)
	numfds = ConnectionNumber(DISP) + 1;
#endif

    if (globals.debug & DBG_EVENT)
	fprintf(stderr, "%s:%d: mjm_read_events %u\n", __FILE__, __LINE__, ret_mask);
    for (;;) {
	globals.ev.ctr = event_freq;
	/*
	 * The above line clears the flag indicating that an event is
	 * pending.  So if an event comes in right now, the flag will be
	 * set again needlessly, but we just end up making an extra call.
	 * Also, be careful about destroying the magnifying glass while
	 * drawing on it.
	 */

#if !FLAKY_SIGPOLL

	if (event_freq < 0) {	/* if SIGPOLL works */
	    if (!XtPending()) {
		sigset_t oldsig;

		(void) sigprocmask(SIG_BLOCK, &all_signals, &oldsig);
			// mjm 

// print oldsig
//$1 = {__val = {0, 0, 0, 0, 0, 0, 0, 0, 0, 140733461823488, 12, 140737314969650, 1, 0, 7634233, 7634223}}
//  $2 = {__val = {0, 0, 0, 0, 0, 0, 0, 0, 0, 140733461823488, 12, 140737314969650, 1, 0, 0, 0}}
		//	oldsig.__val[14]=0; oldsig.__val[15]=0; 
		for (;;) {
#ifdef SHOW_SIG_FLAGS
		    /* this gives HUGE output ... */
		    if (globals.debug & DBG_EVENT)
			fprintf(stderr, "%s:%d: sig_flags = %d\n",
				__FILE__, __LINE__, sig_flags);
#endif
		    while (sig_flags) {
			flags_to_sigproc[sig_flags]();
		    }

		    if (XtPending()) { // mjm 
//			fprintf(stderr, " XtPending true %d     \n",cnt); ++cnt;
			break; }
//			fprintf(stderr, " NO XtPending true %d     \n",cnt); ++cnt;

		    if (globals.ev.flags & ret_mask) {
			(void) sigprocmask(SIG_SETMASK, &oldsig, (sigset_t *) NULL);
			return globals.ev.flags;
		    }
		    (void) sigsuspend(&oldsig);
		} // for mjm 
		(void) sigprocmask(SIG_SETMASK, &oldsig, (sigset_t *) NULL);
	    }
	}
	else

#endif /* not FLAKY_SIGPOLL */

	{
	    for (;;) {
		struct xio	*ip;

		if (globals.debug & DBG_EVENT)
		    fprintf(stderr, "%s:%d: (flaky) sig_flags = %d\n",
			    __FILE__, __LINE__, sig_flags);
		while (sig_flags) {
		    sigset_t oldsig;

		    (void) sigprocmask(SIG_BLOCK, &all_signals, &oldsig);

		    while (sig_flags) {
			flags_to_sigproc[sig_flags]();
		    }

		    (void) sigprocmask(SIG_SETMASK, &oldsig,
				       (sigset_t *) NULL);
		}

		if (XtPending())
		    break;

		if (globals.ev.flags & ret_mask)
		    return globals.ev.flags;

		/* If a SIGUSR1 signal comes right now, then it will wait
		   until an X event or another SIGUSR1 signal arrives. */

#if HAVE_POLL
		if (globals.debug & DBG_EVENT)
		    fprintf(stderr, "%s:%d: have_poll!\n",
			    __FILE__, __LINE__);
		if (io_dirty) {
		    struct pollfd *fp;

		    if (num_fds > max_fds) {
			if (fds != NULL) free(fds);
			fds = xmalloc(num_fds * sizeof *fds);
			max_fds = num_fds;
			fds->fd = ConnectionNumber(DISP);
			fds->events = POLLIN;
		    }
		    fp = fds + 1;
		    for (ip = iorecs; ip != NULL; ip = ip->next) {
			fp->fd = ip->fd;
			fp->events = ip->xio_events;
			ip->pfd = fp;
			++fp;
		    }
		    io_dirty = False;
		}

		for (;;) {
		    if (poll(fds, num_fds, -1) >= 0) {
			for (ip = iorecs; ip != NULL; ip = ip->next) {
			    int revents = ip->pfd->revents;

			    if (revents & POLLIN && ip->read_proc != NULL)
				(ip->read_proc)(ip->fd, ip->data);
			    if (revents & POLLOUT && ip->write_proc != NULL)
				(ip->write_proc)(ip->fd, ip->data);
			}
			break;
		    }

		    if (errno == EINTR)
			break;

		    if (errno != EAGAIN) {
			perror("xdvi: poll");
			break;
		    }
		}
#else /* HAVE_POLL */
		if (globals.debug & DBG_EVENT)
		    fprintf(stderr, "%s:%d: NOT have_poll!\n",
			    __FILE__, __LINE__);
		FD_ZERO(&readfds);
		FD_ZERO(&writefds);
		FD_SET(ConnectionNumber(DISP), &readfds);
		for (ip = iorecs; ip != NULL; ip = ip->next) {
		    if (ip->xio_events & XIO_IN)
			FD_SET(ip->fd, &readfds);
		    if (ip->xio_events & XIO_OUT)
			FD_SET(ip->fd, &writefds);
		}

		for (;;) {
		    if (select(numfds, &readfds, &writefds, (fd_set *) NULL,
			       (struct timeval *) NULL) >= 0) {
			for (ip = iorecs; ip != NULL; ip = ip->next) {
			    if (FD_ISSET(ip->fd, &readfds) && ip->read_proc != NULL) {
				if (globals.debug & DBG_EVENT)
				    fprintf(stderr, "%s:%d: reading from %d\n",
					    __FILE__, __LINE__, ip->fd);
				(ip->read_proc)(ip->fd, ip->data);
			    }
			    if (FD_ISSET(ip->fd, &writefds) && ip->write_proc != NULL) {
				if (globals.debug & DBG_EVENT)
				    fprintf(stderr, "%s:%d: writing to %d\n",
					    __FILE__, __LINE__, ip->fd);
				(ip->write_proc)(ip->fd, ip->data);
			    }
			}
			break;
		    }

		    if (errno == EINTR)
			break;

		    if (errno != EAGAIN) {
			perror("xdvi: select");
			break;
		    }
		}
#endif /* HAVE_POLL */
	    }
	}

	XtAppNextEvent(globals.app, &event);

#ifdef MOTIF
	if ((resource.expert_mode & XPRT_SHOW_TOOLBAR) != 0)
	    TipAppHandle(globals.app, &event);
#endif

#if HAVE_XI21
	if (event.xany.type == GenericEvent
	  && event.xcookie.extension == xi2_opcode) {
	    if (!XGetEventData(DISP, &event.xcookie)) {
		TRACE_EVENTS((stderr,
		  "Received XI2 event, of type %d, with no cookie",
		  event.xcookie.evtype));
	    }
	    else {
		switch (event.xcookie.evtype) {

		    case XI_HierarchyChanged:
			xi2_ev_hierchange(event.xcookie.data);
			break;

		    case XI_DeviceChanged:
			xi2_ev_devchange(event.xcookie.data);
			break;

		    case XI_Enter:
			xi2_ev_enter(event.xcookie.data);
			break;

		    case XI_Motion:
			/* This needs to be filled in by the client */
			((XIDeviceEvent *) event.xcookie.data)->serial
			  = event.xany.serial;
			xi2_ev_motion(event.xcookie.data);
			break;

		    default:
			TRACE_EVENTS((stderr,
			  "Received XI2 event of unknown type %d",
			  event.xcookie.evtype));

		}
		XFreeEventData(DISP, &event.xcookie);
	    }
	    continue;
	}
#endif /* HAVE_XI21 */

	if (resized)
	    get_geom();

	if (event.xany.window == magnifier.win && event.type == Expose) {
	    handle_expose((Widget) NULL, (XtPointer)&magnifier, &event,
			  (Boolean *) NULL);
	    continue;
	}
	else if (globals.broken_motif_event_handling &&
		 (globals.cursor.flags & (CURSOR_RULER | CURSOR_TEXT))) {
	    /* In this case, Act_motion() and Act_release() are not called properly
	     * for updating the ruler/text selection (it works with the magnifier though),
	     * so we need to invoke them ourselves here: */
	    if (event.type == MotionNotify)
		Act_motion(NULL, &event, NULL, NULL);
	    else if (event.type == ButtonRelease)
		Act_release(NULL, &event, NULL, NULL);
	}

#ifdef MOTIF
	if (XtIsRealized(globals.widgets.top_level)
	    && event.xany.window == XtWindow(globals.widgets.clip_widget)
	    && event.type == KeyPress) { /* workaround for #610206 */
	    motif_translations_hack();
	}
#else
	if (resource.expert_mode & XPRT_SHOW_BUTTONS)
	    SubMenuHandleEvent(globals.app, &event);
#endif
	XtDispatchEvent(&event);
    }
}






















//////////////////////////////////////////////////

#if 1 


Boolean mjm_load_dvi_file(
#if !DELAYED_MKTEXPK
	      Boolean load_fonts,
#endif
	      dviErrFlagT *errflag)
{
    unsigned int old_page_w, old_page_h;
    static ino_t dvi_inode = 0;

    TRACE_FILES((stderr, "load_dvi_file: going to read %p", (void *)m_dvi_fp));
    
    if (resource.use_temp_fp && m_dvi_fp != NULL) {
	/* in this case, reload only if file has been written completely */
	*errflag = NO_ERROR;
	fseek(m_dvi_fp, 0L, SEEK_SET);
	if (!process_preamble(m_dvi_fp, errflag)
	    || !find_postamble(m_dvi_fp, errflag)
	    || !read_postamble(m_dvi_fp, errflag,
#if DELAYED_MKTEXPK
			       False, False
#else
			       True
#endif
			       )) {
	    TRACE_FILES((stderr, "reading of %p failed: %s!",
			 (void *)m_dvi_fp,
			 get_dvi_error(*errflag)));
	    return False;
	}
    }

    old_page_w = globals.page.w;
    old_page_h = globals.page.h;

    *errflag = NO_ERROR;

#if DELAYED_MKTEXPK
    /* use same trick with reading postamble twice to first find names of PK fonts
       that need to be created. */
    reset_missing_font_count();
    kpse_set_program_enabled(kpse_any_glyph_format, False, kpse_src_compile);
#endif
    
    if (!internal_open_dvi(globals.dvi_name, errflag,
#if DELAYED_MKTEXPK
			   True, False
#else
			   load_fonts
#endif
			   )) {
	XClearWindow(DISP, mane.win);
	xdvi_bell();
	statusline_info(STATUS_MEDIUM, "%s: %s%s", globals.dvi_name, get_dvi_error(*errflag),
			 resource.watch_file > 0.0 ? "; will try to reload ..." : " (click on window to reload)");
	close_old_filep();

	return False;
    }
#if DELAYED_MKTEXPK
    /* second time */
    kpse_set_program_enabled(kpse_any_glyph_format, True, kpse_src_compile);
    
    if (!internal_open_dvi(globals.dvi_name, errflag, True, True)) {
	XClearWindow(DISP, mane.win);
	xdvi_bell();
	statusline_info(STATUS_MEDIUM, "%s: %s%s", globals.dvi_name, get_dvi_error(*errflag),
			 resource.watch_file > 0.0 ? "; will try to reload ..." : " (click on window to reload)");
	close_old_filep();
	return False;
    }
#endif
    else { /* success */
	if (fstatbuf.st_ino != dvi_inode) {
	    dvi_inode = fstatbuf.st_ino;
	    form_dvi_property();
	    set_dvi_property();
	}
	if (globals.page.w != old_page_w || globals.page.h != old_page_h)
	    reconfig();

	htex_reinit();

	globals.cursor.flags &= ~CURSOR_CORRUPTED;
	return True;
    }
}
/*======================================================================*\

Copyright (c) 1990-2016  Paul Vojta and others

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


void mjm_do_pages(void)
{
mjm_made_it_to_mjm_do_pages=1;
    if (globals.debug & DBG_BATCH) {
	
	(void)mjm_read_events(EV_GT_IDLE);
	for (current_page = 0; current_page < total_pages; ++current_page) {
	    if (resource.keep_flag) {
		home_action = NULL;
	    }
	    else {
		home_action = home;
	    }
	    globals.warn_spec_now = resource.warn_spec;
#if PS_GS
	    for (;;) {
		mjm_redraw_page();
		(void) mjm_read_events(EV_NOWAIT);
		if (!(globals.ev.flags & (EV_NEWPAGE | EV_NEWDOC | EV_RELOAD)))
		    break;
		globals.ev.flags = EV_IDLE;
	    }
#else
	    mjm_redraw_page();
#endif
	}
	xdvi_exit(EXIT_SUCCESS);
    }
    else {
	for (;;) {	/* normal operation */
	    (void) mjm_read_events(EV_GT_IDLE);
// mjm 	    
if (globals.ev.flags & ( EV_RELOAD )) {
//MM_ERR(" reload")
		      		    fprintf(stderr, "  reload! %s %d \n", __FILE__,__LINE__); 

}
 	    TRACE_EVENTS((stderr, "globals.ev.flags: %d; ev_newpage: %d, ev_newdoc: %d, ev_reload: %d\n",
			  globals.ev.flags, EV_NEWPAGE, EV_NEWDOC, EV_RELOAD));
	    /* NOTE: reloading must be checked first! */
	    if (globals.ev.flags & (EV_NEWPAGE | EV_NEWDOC | EV_RELOAD | EV_PS_TOGGLE)) {
		TRACE_EVENTS((stderr, "EV_NEWPAGE | ..."));
		globals.ev.flags &= ~(EV_NEWPAGE | EV_EXPOSE | EV_PS_TOGGLE);
		if (globals.ev.flags & EV_RELOAD) {
		    dviErrFlagT errflag;

		    globals.ev.flags &= ~EV_RELOAD;
		    if (mjm_load_dvi_file(
#if !DELAYED_MKTEXPK
				      True,
#endif
				      &errflag)) {
#if PS_GS
			if (resource.gs_alpha) {
			    /* restart gs so that user has a method for fixing GS artifacts with gs_alpha
			       by using `reload' (see also GS_PIXMAP_CLEARING_HACK) */
			    ps_destroy();
			}
#endif
			statusline_info(STATUS_SHORT, "File reloaded.");
		    }
		    /* 		    else { */
		    /* 			statusline_info(STATUS_MEDIUM, "File corrupted, not reloading."); */
		    /* 		    } */
		}
		if (globals.ev.flags & EV_NEWDOC) {
		    dviErrFlagT errflag;
		    TRACE_EVENTS((stderr, "EV_NEWDOC!"));
		    /*  		    fprintf(stderr, "newdoc!\n"); */
		    TRACE_FILES((stderr, "current page: %d", current_page));
		    /*  		    file_history_set_page(current_page); */
		    globals.ev.flags &= ~EV_NEWDOC;
		    if (load_dvi_file(
#if !DELAYED_MKTEXPK
				      True,
#endif
				      &errflag)) {
			statusline_append(STATUS_SHORT, "Opened ", "Opened \"%s\"", globals.dvi_name);
			/* 			statusline_info(STATUS_SHORT, "Opened \"%s\"", globals.dvi_name); */
			TRACE_FILES((stderr, "Adding to history: |%s|\n", globals.dvi_name));
			if (file_history_push(globals.dvi_name)) { /* it's a new file, add to history */
			    TRACE_FILES((stderr, "New entry!"));
			    filehist_menu_add_entry(globals.dvi_name);
			}
			else { /* only need to move existing elements to new positions */
			    TRACE_FILES((stderr, "Existing entry!\n"));
			    filehist_menu_refresh();
			}
		    }
		}

		can_exposures(&mane);
		can_exposures(&magnifier);

#if PS && PS_GS && GS_PIXMAP_CLEARING_HACK
		if (had_ps_specials && !MAGNIFIER_ACTIVE) {
		    erasepage_gs();
		    had_ps_specials = False;
		}
#endif /* PS && PS_GS && GS_PIXMAP_CLEARING_HACK */

		if (globals.dvi_file.bak_fp != NULL) {
		    TRACE_EVENTS((stderr, "mjm_redraw_page()"));
		    mjm_redraw_page();
		}
		else {
		    TRACE_EVENTS((stderr, "dvi_file_changed()"));
		    (void)dvi_file_changed();
		}
	    }
	    else if (globals.ev.flags & EV_PAGEHIST_GOTO_PAGE) {
		int pageno;
		globals.ev.flags &= ~EV_PAGEHIST_GOTO_PAGE;
		pageno = check_goto_page(page_history_get_page(), False);
		goto_page(pageno, resource.keep_flag ? NULL : home, False);
		TRACE_FILES((stderr, "got page: %d", pageno));
	    }
	    else if (globals.ev.flags & EV_FILEHIST_GOTO_PAGE) {
		int pageno;
		globals.ev.flags &= ~EV_FILEHIST_GOTO_PAGE;
		pageno = check_goto_page(file_history_get_page(), True);
		goto_page(pageno, resource.keep_flag ? NULL : home, False);
		TRACE_FILES((stderr, "got page: %d", pageno));
	    }
	    else if (globals.ev.flags & EV_PAGEHIST_INSERT) {
		globals.ev.flags &= ~EV_PAGEHIST_INSERT;
		page_history_insert(current_page);
	    }
	    else if (globals.ev.flags & EV_FIND_CANCEL) {
		/* NOTE: This must be done before checking for expose() */
		globals.ev.flags &= ~EV_FIND_CANCEL;
	    }
	    else if (globals.ev.flags & EV_ANCHOR) {
		/*
		 * Similar to forward search: search for a htex anchor.
		 * This needs to come before the next case which does the redraw_page(),
		 * otherwise anchors for the current page might not be drawn at all:
		 * anchor_search() sets the info later used by htex_draw_anchormarkers(),
		 * which is invoked by redraw_page().
		 */
		
		/* switch off the link cursor */
		globals.cursor.flags &= ~CURSOR_LINK;

		if (dvi_file_changed())
		    continue;
		
		anchor_search(g_anchor_pos);
		
		/* added, otherwise anchors on same page may not be drawn */
		/* NOTE: This caused a crash when clicking on link "Langer" on p3 of diss.dvi
		   since the color stack wasn't allocated for the target pages. Removed for the
		   time being, since I can't reproduce the problem for which it was introduced.
		   CVS commit message was:
		   "fix anchor drawing for other window"
		*/
		/* redraw(&mane);  */
		
		globals.ev.flags &= ~EV_ANCHOR;
	    }
	    else if (globals.ev.flags & EV_SRC) {
		/*
		 * Source special operations are deferred to here because
		 * they call geom_scan(), which may call define_font(),
		 * which may call makefont(), which may call read_events()
		 * recursively.
		 */
		if (globals.src.fwd_string != NULL) {
		    const char *s = globals.src.fwd_string;

		    if (dvi_file_changed())
			continue;
		    
		    source_forward_search(s);
		    globals.ev.flags &= ~EV_SRC;
		    globals.src.fwd_string = NULL;

		    /* de-iconify window if needed, and raise it */
		    XMapRaised(XtDisplay(globals.widgets.top_level), XtWindow(globals.widgets.top_level));
		    raise_message_windows();
		}
		else if (source_reverse_x != -1) {
		    if (dvi_file_changed())
			continue;
		    
		    source_reverse_search(source_reverse_x, source_reverse_y, True);
		    globals.ev.flags &= ~EV_SRC;
		}
		else {
		    source_special_show(source_show_all);
		    globals.ev.flags &= ~EV_SRC;
		}
	    }
	    /* support for `-findstring' */
	    else if (globals.ev.flags & EV_FIND) {
		if (dvi_file_changed())
		    continue;

		if (resource.find_string != NULL) { /* not first call */
		    dvi_find_string(resource.find_string, False);
		    resource.find_string = NULL;
		}
		else { /* actually should never arrive here?? */
		    dvi_find_string(NULL, True);
		}
		globals.ev.flags &= ~EV_FIND;
	    }
	    else if (globals.ev.flags & EV_MAG_MOVE) {
		MYTRACE((stderr, "moving mag!\n"));
		move_magnifier();
	    }
	    else if (globals.ev.flags & EV_EXPOSE) {
		if (magnifier.min_x < MAXDIM) {
		    /*  		    fprintf(stderr, "magnifier < maxdim!\n"); */
		    if (mane.min_x >= MAXDIM) {
			/*  			fprintf(stderr, "mane >= maxdim!\n"); */
			globals.ev.flags &= ~EV_EXPOSE;
		    }
		    redraw(&magnifier);
		}
		else {
		    /* see comment in mag.c */
		    globals.ev.flags &= ~EV_EXPOSE;
		    if (mane.min_x < MAXDIM) {
			redraw(&mane);
		    }
		}
	    }
	    else if (globals.ev.flags & EV_CURSOR) {
		/*
		 * This code eliminates unnecessary calls to XDefineCursor,
		 * since this is a slow operation on some hardware (e.g., S3
		 * chips).
		 */
		XSync(DISP, False);
		if (!XtPending()) {
		    Cursor curr;
			
		    if (globals.cursor.flags & CURSOR_DRAG_V)
			curr = globals.cursor.drag_v;
		    else if (globals.cursor.flags & CURSOR_DRAG_H)
			curr = globals.cursor.drag_h;
		    else if (globals.cursor.flags & CURSOR_DRAG_A)
			curr = globals.cursor.drag_a;
		    else if (resource.mouse_mode == MOUSE_MODE3)
			curr = globals.cursor.mode3;
		    else if (resource.mouse_mode == MOUSE_MODE2)
			curr = globals.cursor.mode2;
		    else if (globals.cursor.flags & CURSOR_LINK)
			curr = globals.cursor.link;
		    else if (globals.cursor.flags & CURSOR_MAG)
			curr = globals.cursor.empty;
		    else if (globals.cursor.flags & CURSOR_CORRUPTED)
			curr = globals.cursor.corrupted;
		    else if (globals.pausing.flag)
			curr = globals.cursor.pause;
		    else
			curr = globals.cursor.mode1;
		    XDefineCursor(DISP, CURSORWIN, curr);
		    globals.ev.flags &= ~EV_CURSOR;
		}
	    }
	    XFlush(DISP);
	}
    }
}


#if 0 
////////////////////////////////////////////////////////////
    /* add_w has been initialized by create_menu_buttons() call above.
       Reset to 0 if we're in expert mode. */
    if (!(resource.expert_mode & XPRT_SHOW_BUTTONS)) {
	*add_w = 0;
    }
    if (resource.expert_mode & XPRT_SHOW_STATUSLINE) {
	/* FIXME: Unfortunately the statusline hasn't been created at this point for Xaw,
	   so the value is still the built-in default.
	*/
	*add_h += get_statusline_height();
    }
}
#endif



/*
 * Unfortunately this must be a callback, for the file selector ...
 * This is the second part of Main: Create all widgets, initialize the DVI file,
 * and enter the main event loop.
 */

void mjm_run_dvi_file(const char *filename, void *data)
{
    Boolean tried_dvi_ext = False;
    struct startup_info *cb = (struct startup_info *)data;
    
#ifdef MOTIF
    Widget tool_bar = 0;
    Widget form = 0;
#endif

    char *title_name = NULL;
    char *icon_name = NULL;
    dviErrFlagT errflag = NO_ERROR;
    
    int add_w = 0, add_h = 0;
    Dimension main_win_w, main_win_h;
    
    UNUSED(data);
    ASSERT(filename != NULL, "filename must have been initialized here!");
    
    globals.dvi_name = xstrdup(filename);
    file_history_push(globals.dvi_name);
    
    TRACE_FILES((stderr, "globals.dvi_name is: |%s| %p\n", globals.dvi_name, globals.dvi_name));

    globals.dvi_file.dirname = get_dir_component(globals.dvi_name);
    xdvi_assert(XDVI_VERSION_INFO, __FILE__, __LINE__,
		globals.dvi_file.dirname != NULL,
		"globals.dvi_name (%s) must contain a dir component",
		globals.dvi_name);
    globals.dvi_file.dirlen = strlen(globals.dvi_file.dirname);

    form_dvi_property();
    
    /*
      If `unique' is active, we may need to pass the file to a different instance of xdvi:
    */
    if (resource.unique) {
	Window w1 = 0, w2 = 0;
	if ((w1 = get_xdvi_window_id(True, NULL)) != 0 || (w2 = get_xdvi_window_id(False, NULL)) != 0) {
	    if (w1 != 0) { /* another xdvi instance, same file: reload */
		w2 = w1;
		set_string_property("", atom_reload(), w2);
	    }
	    else { /* another xdvi instance, different file: load new file */
		set_string_property(globals.dvi_name, atom_newdoc(), w2);
	    }
	    if (cb->page_arg != NULL) { /* switch to different page */
		if (strlen(cb->page_arg) == 0) { /* special case: treat `+' as last page */
		    set_string_property("+", atom_newpage(), w2);
		}
		else {
		    set_string_property(cb->page_arg, atom_newpage(), w2);
		}
	    }
	    else {
		/* if page_arg is empty, go to 1st page so that this page is
		   inserted into the page history (fix for #1044891) */
		set_string_property("1", atom_newpage(), w2);
	    }
	    /* in all cases, raise the window of the other instance */
	    set_string_property("", atom_raise(), w2);
	    xdvi_exit(EXIT_SUCCESS);
	}
	else if (resource.src_fork) {
	    do_fork();
	}
    }

    /*
      Similar for forward search or string search:
    */
    if (resource.src_pos != NULL || resource.find_string != NULL) {
	Window w;
	if ((w = get_xdvi_window_id(True, NULL)) != 0) {
	    /* another instance of xdvi running, displaying the same file */
	    TRACE_CLIENT((stderr, "Match; changing property of client and exiting ..."));
	    if (resource.src_pos != NULL)
		set_sourceposition_property(resource.src_pos, w);
	    else
		set_stringsearch_property(resource.find_string, w);
	    xdvi_exit(EXIT_SUCCESS);
	}
	else if (resource.src_fork) {
	    do_fork();
	}
    }

    /* Needed for source specials and for calling ghostscript. */
    xputenv("DISPLAY", XDisplayString(DISP));


    if (globals.debug) {
	fprintf(stderr, "%s %s, kpathsea: %s\n", XDVIK_PROGNAME, XDVI_VERSION_INFO, kpathsea_version_string);
	fprintf(stderr,
		"configured with: ppi=%d shrink=%d mfmode=%s alt_font=%s paper=%s\n",
		resource.pixels_per_inch,
		currwin.shrinkfactor,
		resource.mfmode ? resource.mfmode : "<NONE>",
		resource.alt_font,
		resource.paper);
    }
    
    kpse_set_program_enabled(kpse_any_glyph_format, resource.makepk, kpse_src_compile);
    /* janl 16/11/98: I have changed this. The above line used to
       say the settings in resource.makepk was supplied on the
       commandline, resulting in it overriding _all other_
       settings, derived from the environment or texmf.cnf, no
       matter what the value. The value in resource.makepk could
       be the compile-time default...

       Personaly I like the environment/texmf.cnf to override
       resources and thus changed the 'level' of this setting to
       kpse_src_compile so the environment/texmf.cnf will override
       the values derived by Xt.

       Previous comment here:

       ``Let true values as an X resource/command line override false
       values in texmf.cnf/envvar.''  */

    /*
     *		Step 2:  Settle colormap issues.  This should be done before
     *		other widgets are created, so that they get the right
     *		pixel values.  (The top-level widget won't have the right
     *		values, but I don't think that makes any difference.)
     */
    
#ifdef XSERVER_INFO
    print_xserver_info();
#endif

    create_colormaps();

    
#ifdef TESTING_OPEN_FILES
    fprintf(stderr, "open_max: %ld\n", OPEN_MAX);
    for (i = 0; i < OPEN_MAX - 10; i++) {
	FILE *fp;
	if ((fp = fopen("/tmp/foo", "r")) == NULL) {
	    perror("fopen");
	    xdvi_exit(EXIT_FAILURE);
	}
    }
    fprintf(stderr, "opened %d files.\n", i);
#endif
    

    /* toolbar code may open files, but we have no check close_a_file() in
       the toolbar code; so do this before prescan() possibly opens lots of files.
    */
#ifdef MOTIF
    globals.widgets.main_row = XmCreateMainWindow(globals.widgets.top_level, "main", NULL, 0);
    
    create_menu_buttons(globals.widgets.main_row, &globals.widgets.menu_bar);

    /* seems to be needed for enabling `XmNhighlightOnEnter' for the toolbar buttons
       - is this the correct place to do it? */
    XtVaSetValues(globals.widgets.top_level, XmNkeyboardFocusPolicy, (XtArgVal)XmPOINTER, NULL);

    form = XtVaCreateWidget("form", xmFormWidgetClass, globals.widgets.main_row,
			    XmNshadowThickness, 0,
			    NULL);

    if (resource.main_translations != NULL) {
	XtOverrideTranslations(form, XtParseTranslationTable(resource.main_translations));
    }

    
#if defined(USE_PANNER) && USE_XAW_PANNER
    panner = XtVaCreateWidget("panner", pannerWidgetClass, form,
			      XmNtopAttachment, XmATTACH_FORM,
			      XmNleftAttachment, XmATTACH_FORM,
			      XmNleftOffset, 20,
			      XmNrightAttachment, XmATTACH_OPPOSITE_FORM,
			      XmNrightOffset, -resource.pagelist_width + 20,
			      XmNtopOffset, 2,
			      XmNleftOffset, 2,
			      XtNheight, 60,
			      XtNwidth, resource.pagelist_width - 50,
			      XtNsliderX, 5,
			      XtNsliderY, 7,
			      XtNinternalSpace, 0,
			      XtNshadowThickness, 0,
			      NULL);
#endif
#endif

    /* use bounding box for highlighting if our visual isn't TrueColor
       (not worth the trouble ...) */
    if (G_visual->class != TrueColor) {
	resource.match_highlight_inverted = False;
    }
    
    /*
     *		Step 3:  Initialize the dvi file and set titles.
     */

#if FREETYPE
    /*
      At this point DISP, G_visual, G_depth and G_colormap must
      be defined. Also, init_t1_lookup() must go before internal_open_dvi(),
      since read_postamble will define some fonts and insert them into
      fontmaps_hash, but we need a clean fontmaps_hash for detecting
      duplicate entries in the map file.
    */

    if (resource.freetype) {
	if (!init_t1_lookup()) {
	    /* nag 'em with a popup so that they'll do something about this */
	    popup_message(globals.widgets.top_level,
		      MSG_ERR,
		      "Direct Type 1 font rendering via FreeType gives you "
		      "many benefits, such as:\n"
		      " - quicker startup time, since no bitmap fonts need "
		      "to be generated;\n"
		      " - saving disk space for storing the bitmap fonts.\n"
		      "To fix this error, check that the file `ps2pk.map' "
		      "is located somewhere in your XDVIINPUTS path.  "
		      "Have a look at the xdvi wrapper shell script "
		      "(type \"which xdvi\" to locate that shell script) "
		      "for the current setting of XDVIINPUTS.",
		      "Could not load any of the map files listed in xdvi.cfg "
		      "- disabling FreeType.");
	    resource.freetype = False;
	}
    }
#endif /* FREETYPE */

#if DELAYED_MKTEXPK
    /* Open and initialize the DVI file. First, disable creation of PK fonts
     * so that we can count the missing fonts that are to be generated. */
    kpse_set_program_enabled(kpse_any_glyph_format, False, kpse_src_compile);
#endif
    
    setup_signal_handlers(False);

#if !DELAYED_MKTEXPK
    /* Notify users that fonts are being created. This is just a hack
       and no replacement for true asynchronous font creation since it
       doesn't give details (is just invoked if startup takes somewhat
       longer) and freezes during font creation.
    */
    register_font_popup();
#endif
    
    /* open and initialize the DVI file, but don't read the fonts yet */
    if (!internal_open_dvi(globals.dvi_name, &errflag, True
#if DELAYED_MKTEXPK
			   , False /* read fonts, but don't initialize data structures */
#endif
			   )) {
	if (tried_dvi_ext) {
	    XDVI_FATAL((stderr, "Could not open %s: %s, and %s.dvi doesn't exist either - exiting.",
			globals.dvi_name, get_dvi_error(errflag), globals.dvi_name));
	}
	else {
	    XDVI_FATAL((stderr, "Could not open %s: %s.",
			globals.dvi_name, get_dvi_error(errflag)));
	}
    }

#if DELAYED_MKTEXPK
    fprintf(stderr, "after opening ...\n");
    /* Now re-enable PK creation and read the postamble for a second time.
     * FIXME: Actually we don't need this re-reading, could as well read the
     * entire thing in the first run, not quit early and correctly initialize
     * the fonts without creating them. */
    kpse_set_program_enabled(kpse_any_glyph_format, resource.makepk, kpse_src_compile);

    if (!internal_open_dvi(globals.dvi_name, &errflag, True, True)) {
	if (tried_dvi_ext) {
	    XDVI_FATAL((stderr, "Could not open %s: %s, and %s.dvi doesn't exist either - exiting.",
			globals.dvi_name, get_dvi_error(errflag), globals.dvi_name));
	}
	else {
	    XDVI_FATAL((stderr, "Could not open %s: %s.",
			globals.dvi_name, get_dvi_error(errflag)));
	}
    }
#else
    unregister_font_popup();
#endif
    
    if (cb->page_arg != NULL) {
	if (cb->page_arg[0] == '\0') { /* empty page_arg -> goto last page */
	    current_page = total_pages - 1;
	    page_history_insert(current_page);
	}
	else {
	    char *testptr;
	    current_page = strtoul(cb->page_arg, &testptr, 10) - 1;
	    if (*testptr != '\0') {
		XDVI_FATAL((stderr, "Invalid page number: `%s'.", cb->page_arg));
	    }
	    current_page = check_goto_page(current_page, True);
	}
    }
    else {
	page_history_insert(current_page);
    }
    file_history_set_page(current_page);
    
    ASSERT(globals.dvi_file.bak_fp != NULL, "Backup file pointer must have been initialized here");
    if (resource.prescan) {
	prescan(globals.dvi_file.bak_fp);
    }

    globals.page.unshrunk_w = pageinfo_get_page_width(current_page);
    globals.page.unshrunk_h = pageinfo_get_page_height(current_page);
    TRACE_GUI((stderr, "globals.page.unshrunk_w: %d, h: %d; window: %d, %d",
	       globals.page.unshrunk_w, globals.page.unshrunk_h,
	       pageinfo_get_window_width(current_page),
	       pageinfo_get_window_height(current_page)));
    
    init_page();

    /*
     *		Step 4:  Create widgets, and set initial window size.
     */

    /* currently these override expert mode - using this is deprecated
       in favour of `-expertmode'; inform user about this: */
    if (resource.statusline) {
	XDVI_WARNING((stderr, "The option/X resource `statusline' is obsolete; "
		      "use `-expertmode <flag>' instead, e.g. `-expertmode 1'\n"
		      "to switch on the status line, or `-expertmode 6'\n"
		      "to switch it off. See the xdvi man page for details."));
	resource.expert_mode |= XPRT_SHOW_STATUSLINE;
    }

    /*      XtRealizeWidget(globals.widgets.top_level); */

#ifdef MOTIF
    tool_bar = create_toolbar(globals.widgets.main_row, globals.widgets.menu_bar);
    if (resource.main_translations != NULL) {
	XtOverrideTranslations(tool_bar, XtParseTranslationTable(resource.main_translations));
    }
#endif

    create_widgets(
#ifdef MOTIF
		   tool_bar, form,
#endif
		   &add_w, &add_h);

    TRACE_GUI((stderr, "add_w = %d, add_h = %d\n", add_w, add_h));
    /*  fprintf(stderr, "geometry xdvirc: |%s|, orig: |%s|\n", resource.xdvirc_geometry, resource.geometry); */
    
    /*
     *	Set initial window size.
     *	This needs to be done before colors are assigned because if
     *	-s 0 is specified, we need to compute the shrink factor
     *	(which in turn affects whether init_pix is called).
     */
    set_windowsize(&main_win_w, &main_win_h, add_w, add_h, False);

    realize_widgets(main_win_w, main_win_h, add_w);

    /* this needs to be done after total_pages is known (via internal_open_dvi) */
    get_icon_and_title(globals.dvi_name, &icon_name, &title_name);
    add_icon(globals.widgets.top_level, title_name, icon_name);
    /* this needs to be done after the widgets have been created */
    set_icon_and_title(icon_name, title_name);
    free(icon_name);
    free(title_name);
    icon_name = title_name = NULL;

    
    G_image = XCreateImage(DISP, G_visual, 1, XYBitmap, 0,
			   (char *)NULL, 0, 0, BMBITS, 0);
    G_image->bitmap_unit = BMBITS;
#ifdef	WORDS_BIGENDIAN
    G_image->bitmap_bit_order = MSBFirst;
#else
    G_image->bitmap_bit_order = LSBFirst;
#endif
    {
	short endian = MSBFirst << 8 | LSBFirst;
	G_image->byte_order = *((char *)&endian);
    }

    /* Store window id for use by get_xdvi_window_id().  */
    {
	long data = XtWindow(globals.widgets.top_level);

	XChangeProperty(DISP, DefaultRootWindow(DISP),
			atom_xdvi_windows(), atom_xdvi_windows(), 32,
			PropModePrepend, (unsigned char *)&data, 1);
	set_dvi_property();
    }


#if HAVE_XI21
    xi2_init();	/* Set up hi-res (smooth) scrolling */
#endif

    /*
     *	Step 5:  Assign colors and GCs.
     *		 Because of the latter, this has to go after the widgets are realized.
     */

    create_gcs();

    create_cursors();

#ifdef MOTIF
#if defined(USE_PANNER) && USE_XAW_PANNER
    XtVaSetValues(panner, XtNsliderWidth, globals.page.w / 2,
		  XtNsliderHeight, globals.page.h / 2,
		  XtNcanvasWidth, globals.page.w,
		  XtNcanvasHeight, globals.page.h,
		  NULL);
    XtManageChild(panner);
    XtAddCallback(panner, XtNreportCallback, panner_cb, (XtPointer)NULL);
#endif
    create_pagelist();
#endif

    /* trigger forward search */
    do_forward_search(resource.src_pos);

    /* trigger string search */
    if (resource.find_string != NULL) {
	globals.ev.flags |= EV_FIND;
    }

    /* trigger anchor search */
    if (resource.anchor_pos != NULL) {
	g_anchor_pos = xstrdup(resource.anchor_pos);
	g_anchor_len = strlen(g_anchor_pos);
	globals.ev.flags |= EV_ANCHOR;
    }

#if defined(MOTIF) && HAVE_XPM
    tb_check_navigation_sensitivity(current_page);
#endif

#if CHECK_APP_FILEVERSION
    check_app_defaults_fileversion();
#endif

    /* can do this only after scrollbars have been realized */
    if (!BROKEN_RECONFIG && (resource.expert_mode & XPRT_SHOW_SCROLLBARS) == 0) {
	toggle_scrollbars();
    }

    /* initialize file watching */
    if (resource.watch_file > 0.0) {
#if XDVI_XT_TIMER_HACK
	watch_file_cb(NULL, NULL);
#else
	XDVI_WARNING((stderr, "Could not redefine XtAppAddTimeOut(); `watchfile' not available."));
#endif
    }

    /* raise `early' message windows */
    raise_message_windows();

    {    
	String args[1];
	char mmode[LENGTH_OF_INT];
	snprintf(mmode, LENGTH_OF_INT, "%d", resource.mouse_mode);
	args[0] = mmode;
	XtCallActionProc(globals.widgets.top_level, "switch-mode", NULL, args, 1);
    }

    /* go to home position on first page to honor -(side|top)margin flags */
    if (!resource.keep_flag)
    	home(False);

    mjm_do_pages();
}

#endif

