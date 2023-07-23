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

unsigned int
read_events(unsigned int ret_mask)
{
static int cnt=0;
    XEvent event;

#if !HAVE_POLL
    if (numfds == 0)
	numfds = ConnectionNumber(DISP) + 1;
#endif

    if (globals.debug & DBG_EVENT)
	fprintf(stderr, "%s:%d: read_events %u\n", __FILE__, __LINE__, ret_mask);
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

