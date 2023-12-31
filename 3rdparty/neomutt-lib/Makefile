###############################################################################
# substitutions by autosetup

PACKAGE=	neomutt
PACKAGE_VERSION=20200501
PACKAGE_DATE=	2020-05-01

# Build-time compiler and flags. These are used for building executables that
# are only used at build-time, e.g., docs/makedoc. These are different from CC /
# CFLAGS / LDFLAGS because the host and target systems might be different.
CC_FOR_BUILD=		 cc
CFLAGS_FOR_BUILD=	-std=c99 -fPIC  -D_ALL_SOURCE=1 -D_GNU_SOURCE=1 -D__EXTENSIONS__ $(EXTRA_CFLAGS_FOR_BUILD)
LDFLAGS_FOR_BUILD=	 $(EXTRA_LDFLAGS_FOR_BUILD)

PWD=		/home/ubuntu/dev/neomutt/neomutt-20200501
AR=		ar
CC=		 cc
CPP=		cc -E
CFLAGS=		 -g -O2 -std=c99 -fPIC  -D_ALL_SOURCE=1 -D_GNU_SOURCE=1 -D__EXTENSIONS__ -I/usr/include -DNCURSES_WIDECHAR -isystem /usr/include/mit-krb5 -I. -I. -Wall $(EXTRA_CFLAGS)
LDFLAGS=	-L/usr/lib -L/usr/lib/x86_64-linux-gnu/mit-krb5 -Wl,-Bsymbolic-functions -Wl,-z,relro -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err $(EXTRA_LDFLAGS)
EXEEXT=		
LIBS=		-ltokyocabinet -lreadline -lidn -lgnutls -lncursesw -ltinfo -lsasl2 -lgpgme -lanl  -L/usr/lib/x86_64-linux-gnu -lgpgme -lassuan -lgpg-error -L/usr/lib/x86_64-linux-gnu -lgpg-error
INTLLIBS=	 
RANLIB=		ranlib
SRCDIR=		.
INSTALL=	install -c
INSTALL_DATA=	install -m 644
MKDIR_P=	mkdir -p
RM=		rm -fr
DEPFILES=	$(ALLOBJS:.o=.Po)
BINFILES=	$(NEOMUTT)
LIBBINFILES+=	$(PGPEWRAP)
LIBBINFILES+=	$(SRCDIR)/contrib/smime_keys

# paths
bindir=		/usr/bin
datadir=	/usr/share/neomutt
docdir=		/usr/share/doc/neomutt
mandir=		/usr/share/man
libdir=		/usr/lib
libexecdir=	/usr/libexec
sysconfdir=	/etc
textdomaindir=	/usr/share/locale

# targets for specific subdirectories
ALL_TARGETS=		all-po all-docs all-contrib
CLEAN_TARGETS=		clean-po clean-docs clean-contrib
INSTALL_TARGETS=	install-po install-docs install-contrib
UNINSTALL_TARGETS=	uninstall-po uninstall-docs uninstall-contrib

VPATH=		$(SRCDIR):$(SRCDIR)/po:$(SRCDIR)/docs:$(SRCDIR)/contrib

ALL_FILES!=	(cd $(SRCDIR) && git ls-files 2>/dev/null) || true

###############################################################################
# neomutt
NEOMUTT=	neomutt$(EXEEXT)
NEOMUTTOBJS=	addrbook.o alias.o browser.o commands.o \
		complete.o compose.o conststrings.o context.o copy.o \
		edit.o editmsg.o enriched.o enter.o flags.o functions.o \
		git_ver.o handler.o hdrline.o help.o hook.o icommands.o index.o init.o \
		keymap.o mailcap.o main.o  menu.o mutt_account.o mutt_attach.o \
		mutt_body.o mutt_header.o mutt_history.o mutt_logging.o mutt_mailbox.o \
		mutt_parse.o mutt_signal.o mutt_socket.o mutt_thread.o \
		muttlib.o mx.o myvar.o pager.o pattern.o postpone.o progress.o query.o \
		recvattach.o recvcmd.o resize.o rfc3676.o score.o send.o sendlib.o \
		sidebar.o smtp.o sort.o state.o status.o system.o version.o \
		mutt_commands.o mutt_config.o command_parse.o

NEOMUTTOBJS+=	monitor.o
NEOMUTTOBJS+=	mjm_main_structs.o mjm_mutt_server.o

CLEANFILES+=	$(NEOMUTT) $(NEOMUTTOBJS)
ALLOBJS+=	$(NEOMUTTOBJS)

###############################################################################
# libautocrypt

###############################################################################
# libpop
LIBPOP=	libpop.a
LIBPOPOBJS=	pop/pop.o pop/pop_auth.o pop/pop_lib.o
CLEANFILES+=	$(LIBPOP) $(LIBPOPOBJS)
MUTTLIBS+=	$(LIBPOP)
ALLOBJS+=	$(LIBPOPOBJS)

###############################################################################
# libnntp
LIBNNTP=	libnntp.a
LIBNNTPOBJS=	nntp/browse.o nntp/complete.o nntp/newsrc.o nntp/nntp.o
CLEANFILES+=	$(LIBNNTP) $(LIBNNTPOBJS)
MUTTLIBS+=	$(LIBNNTP)
ALLOBJS+=	$(LIBNNTPOBJS)

###############################################################################
# libcompmbox
LIBCOMPMBOX=	libcompmbox.a
LIBCOMPMBOXOBJS=compmbox/compress.o
CLEANFILES+=	$(LIBCOMPMBOX) $(LIBCOMPMBOXOBJS)
MUTTLIBS+=	$(LIBCOMPMBOX)
ALLOBJS+=	$(LIBCOMPMBOXOBJS)

###############################################################################
# libstore
LIBSTOREOBJS+=	store/tc.o
LIBSTORE=	libstore.a
LIBSTOREOBJS+=	store/store.o
CLEANFILES+=	$(LIBSTORE) $(LIBSTOREOBJS)
MUTTLIBS+=	$(LIBSTORE)
ALLOBJS+=	$(LIBSTOREOBJS)

###############################################################################
# libgui
LIBGUI=	libgui.a
LIBGUIOBJS=	gui/color.o gui/curs_lib.o gui/mutt_curses.o gui/mutt_window.o gui/reflow.o gui/terminal.o
CLEANFILES+=	$(LIBGUI) $(LIBGUIOBJS)
MUTTLIBS+=	$(LIBGUI)
ALLOBJS+=	$(LIBGUIOBJS)

###############################################################################
# libdebug

###############################################################################
# libmbox
LIBMBOX=	libmbox.a
LIBMBOXOBJS=	mbox/mbox.o
CLEANFILES+=	$(LIBMBOX) $(LIBMBOXOBJS)
MUTTLIBS+=	$(LIBMBOX)
ALLOBJS+=	$(LIBMBOXOBJS)

###############################################################################
# libnotmuch

###############################################################################
# libmaildir
LIBMAILDIR=	libmaildir.a
LIBMAILDIROBJS=	maildir/maildir.o maildir/mh.o maildir/shared.o
CLEANFILES+=	$(LIBMAILDIR) $(LIBMAILDIROBJS)
MUTTLIBS+=	$(LIBMAILDIR)
ALLOBJS+=	$(LIBMAILDIROBJS)

###############################################################################
# libncrypt
LIBNCRYPT=	libncrypt.a
LIBNCRYPTOBJS=	ncrypt/crypt.o ncrypt/cryptglue.o ncrypt/crypt_mod.o
LIBNCRYPTOBJS+=	ncrypt/crypt_gpgme.o ncrypt/crypt_mod_pgp_gpgme.o \
		ncrypt/crypt_mod_smime_gpgme.o
LIBNCRYPTOBJS+=	ncrypt/crypt_mod_pgp_classic.o ncrypt/gnupgparse.o \
		ncrypt/pgp.o ncrypt/pgpinvoke.o ncrypt/pgpkey.o ncrypt/pgplib.o \
		ncrypt/pgpmicalg.o ncrypt/pgppacket.o
LIBNCRYPTOBJS+=	ncrypt/crypt_mod_smime_classic.o ncrypt/smime.o
CLEANFILES+=	$(LIBNCRYPT) $(LIBNCRYPTOBJS)
MUTTLIBS+=	$(LIBNCRYPT)
ALLOBJS+=	$(LIBNCRYPTOBJS)

###############################################################################
# libimap
LIBIMAP=	libimap.a
LIBIMAPOBJS=	imap/auth.o \
		imap/auth_login.o imap/auth_oauth.o imap/auth_plain.o imap/browse.o \
		imap/command.o imap/imap.o imap/message.o imap/utf7.o imap/util.o
LIBIMAPOBJS+=	imap/auth_gss.o
LIBIMAPOBJS+=	imap/auth_sasl.o
CLEANFILES+=	$(LIBIMAP) $(LIBIMAPOBJS)
MUTTLIBS+=	$(LIBIMAP)
ALLOBJS+=	$(LIBIMAPOBJS)

###############################################################################
# libconn
LIBCONN=	libconn.a
LIBCONNOBJS=	conn/connaccount.o conn/conn_globals.o conn/getdomain.o \
		conn/raw.o conn/sasl_plain.o conn/socket.o conn/tunnel.o
LIBCONNOBJS+=	conn/sasl.o
LIBCONNOBJS+=	conn/gui.o
LIBCONNOBJS+=	conn/gnutls.o
CLEANFILES+=	$(LIBCONN) $(LIBCONNOBJS)
MUTTLIBS+=	$(LIBCONN)
ALLOBJS+=	$(LIBCONNOBJS)

###############################################################################
# libhcache
LIBHCACHE=	libhcache.a
LIBHCACHEOBJS=	hcache/hcache.o hcache/serialize.o
CLEANFILES+=	$(LIBHCACHE) $(LIBHCACHEOBJS)
MUTTLIBS+=	$(LIBHCACHE)
ALLOBJS+=	$(LIBHCACHEOBJS)

###############################################################################
# libcompress

###############################################################################
# libbcache
LIBBCACHE=	libbcache.a
LIBBCACHEOBJS=	bcache/bcache.o

CLEANFILES+=	$(LIBBCACHE) $(LIBBCACHEOBJS)
MUTTLIBS+=	$(LIBBCACHE)
ALLOBJS+=	$(LIBBCACHEOBJS)

###############################################################################
# libhistory
LIBHISTORY=	libhistory.a
LIBHISTORYOBJS=	history/history.o

CLEANFILES+=	$(LIBHISTORY) $(LIBHISTORYOBJS)
MUTTLIBS+=	$(LIBHISTORY)
ALLOBJS+=	$(LIBHISTORYOBJS)

###############################################################################
# pgpewrap
PGPEWRAP=	pgpewrap$(EXEEXT)
PGPEWRAPOBJS=	pgpewrap.o
CLEANFILES+=	$(PGPEWRAP) $(PGPEWRAPOBJS)
ALLOBJS+=	$(PGPEWRAPOBJS)

###############################################################################
# libcore
LIBCORE=	libcore.a
LIBCOREOBJS=	core/account.o core/mailbox.o core/neomutt.o

CLEANFILES+=	$(LIBCORE) $(LIBCOREOBJS)
MUTTLIBS+=	$(LIBCORE)
ALLOBJS+=	$(LIBCOREOBJS)

###############################################################################
# libconfig
LIBCONFIG=	libconfig.a
LIBCONFIGOBJS=	config/address.o config/bool.o config/dump.o config/enum.o \
		config/long.o config/mbtable.o config/number.o config/path.o config/quad.o \
		config/regex.o config/set.o config/slist.o config/sort.o \
		config/string.o config/subset.o

CLEANFILES+=	$(LIBCONFIG) $(LIBCONFIGOBJS)
MUTTLIBS+=	$(LIBCONFIG)
ALLOBJS+=	$(LIBCONFIGOBJS)

###############################################################################
# libemail
LIBEMAIL=	libemail.a
LIBEMAILOBJS=	email/attach.o email/body.o email/email.o \
		email/email_globals.o email/envelope.o email/from.o email/mime.o \
		email/parameter.o email/parse.o email/rfc2047.o email/rfc2231.o \
		email/tags.o email/thread.o email/url.o
CLEANFILES+=	$(LIBEMAIL) $(LIBEMAILOBJS)
MUTTLIBS+=	$(LIBEMAIL)
ALLOBJS+=	$(LIBEMAILOBJS)

###############################################################################
# libaddress
LIBADDRESS=	libaddress.a
LIBADDRESSOBJS=	address/address.o address/group.o address/idna.o

CLEANFILES+=	$(LIBADDRESS) $(LIBADDRESSOBJS)
MUTTLIBS+=	$(LIBADDRESS)
ALLOBJS+=	$(LIBADDRESSOBJS)

###############################################################################
# libmutt
LIBMUTT=	libmutt.a
LIBMUTTOBJS=	mutt/base64.o mutt/buffer.o mutt/charset.o mutt/date.o \
		mutt/envlist.o mutt/exit.o mutt/file.o mutt/filter.o \
		mutt/hash.o mutt/list.o mutt/logging.o mutt/mapping.o \
		mutt/mbyte.o mutt/md5.o mutt/memory.o mutt/notify.o \
		mutt/path.o mutt/pool.o mutt/prex.o mutt/regex.o \
		mutt/signal.o mutt/slist.o mutt/string.o
CLEANFILES+=	$(LIBMUTT) $(LIBMUTTOBJS)
MUTTLIBS+=	$(LIBMUTT)
ALLOBJS+=	$(LIBMUTTOBJS)

###############################################################################
# generated
GENERATED=	git_ver.c hcache/hcversion.h
CLEANFILES+=	$(GENERATED)

##############################################################################
# targets
all: $(BINFILES) $(LIBBINFILES) $(ALL_TARGETS)

# compile + dependencies
.c.o:
	$(CC) $(CFLAGS) -MT $@ -MD -MP -MF $(@:.o=.Tpo) -c -o $@ $<
	@mv $(@:.o=.Tpo) $(@:.o=.Po)

$(ALLOBJS):

# mutt
$(NEOMUTT): $(GENERATED) $(NEOMUTTOBJS) $(MUTTLIBS)
	$(CC) -o $@ $(NEOMUTTOBJS) $(MUTTLIBS) $(LDFLAGS) $(LIBS)

# libmutt
$(LIBMUTT): $(PWD)/mutt $(LIBMUTTOBJS)
	$(AR) cr $@ $(LIBMUTTOBJS)
	$(RANLIB) $@
$(PWD)/mutt:
	$(MKDIR_P) $(PWD)/mutt

# libaddress
$(LIBADDRESS): $(PWD)/address $(LIBADDRESSOBJS)
	$(AR) cr $@ $(LIBADDRESSOBJS)
	$(RANLIB) $@
$(PWD)/address:
	$(MKDIR_P) $(PWD)/address

# libautocrypt
$(LIBAUTOCRYPT): $(PWD)/autocrypt $(LIBAUTOCRYPTOBJS)
	$(AR) cr $@ $(LIBAUTOCRYPTOBJS)
	$(RANLIB) $@
$(PWD)/autocrypt:
	$(MKDIR_P) $(PWD)/autocrypt

# libemail
$(LIBEMAIL): $(PWD)/email $(LIBEMAILOBJS)
	$(AR) cr $@ $(LIBEMAILOBJS)
	$(RANLIB) $@
$(PWD)/email:
	$(MKDIR_P) $(PWD)/email

# libpop
$(LIBPOP): $(PWD)/pop $(LIBPOPOBJS)
	$(AR) cr $@ $(LIBPOPOBJS)
	$(RANLIB) $@
$(PWD)/pop:
	$(MKDIR_P) $(PWD)/pop

# libnntp
$(LIBNNTP): $(PWD)/nntp $(LIBNNTPOBJS)
	$(AR) cr $@ $(LIBNNTPOBJS)
	$(RANLIB) $@
$(PWD)/nntp:
	$(MKDIR_P) $(PWD)/nntp

# libcompmbox
$(LIBCOMPMBOX): $(PWD)/compmbox $(LIBCOMPMBOXOBJS)
	$(AR) cr $@ $(LIBCOMPMBOXOBJS)
	$(RANLIB) $@
$(PWD)/compmbox:
	$(MKDIR_P) $(PWD)/compmbox

# libstore
$(LIBSTORE): $(PWD)/store $(LIBSTOREOBJS)
	$(AR) cr $@ $(LIBSTOREOBJS)
	$(RANLIB) $@
$(PWD)/store:
	$(MKDIR_P) $(PWD)/store

# libcompress
$(LIBCOMPRESS): $(PWD)/compress $(LIBCOMPRESSOBJS)
	$(AR) cr $@ $(LIBCOMPRESSOBJS)
	$(RANLIB) $@
$(PWD)/compress:
	$(MKDIR_P) $(PWD)/compress

# libgui
$(LIBGUI): $(PWD)/gui $(LIBGUIOBJS)
	$(AR) cr $@ $(LIBGUIOBJS)
	$(RANLIB) $@
$(PWD)/gui:
	$(MKDIR_P) $(PWD)/gui

# libdebug
$(LIBDEBUG): $(PWD)/debug $(LIBDEBUGOBJS)
	$(AR) cr $@ $(LIBDEBUGOBJS)
	$(RANLIB) $@
$(PWD)/debug:
	$(MKDIR_P) $(PWD)/debug

# libmbox
$(LIBMBOX): $(PWD)/mbox $(LIBMBOXOBJS)
	$(AR) cr $@ $(LIBMBOXOBJS)
	$(RANLIB) $@
$(PWD)/mbox:
	$(MKDIR_P) $(PWD)/mbox

# libnotmuch
$(LIBNOTMUCH): $(PWD)/notmuch $(LIBNOTMUCHOBJS)
	$(AR) cr $@ $(LIBNOTMUCHOBJS)
	$(RANLIB) $@
$(PWD)/notmuch:
	$(MKDIR_P) $(PWD)/notmuch

# libmaildir
$(LIBMAILDIR): $(PWD)/maildir $(LIBMAILDIROBJS)
	$(AR) cr $@ $(LIBMAILDIROBJS)
	$(RANLIB) $@
$(PWD)/maildir:
	$(MKDIR_P) $(PWD)/maildir

# libncrypt
$(LIBNCRYPT): $(PWD)/ncrypt $(LIBNCRYPTOBJS)
	$(AR) cr $@ $(LIBNCRYPTOBJS)
	$(RANLIB) $@
$(PWD)/ncrypt:
	$(MKDIR_P) $(PWD)/ncrypt

# libimap
$(LIBIMAP): $(PWD)/imap $(LIBIMAPOBJS)
	$(AR) cr $@ $(LIBIMAPOBJS)
	$(RANLIB) $@
$(PWD)/imap:
	$(MKDIR_P) $(PWD)/imap

# libconn
$(LIBCONN): $(PWD)/conn $(LIBCONNOBJS)
	$(AR) cr $@ $(LIBCONNOBJS)
	$(RANLIB) $@
$(PWD)/conn:
	$(MKDIR_P) $(PWD)/conn

# libcore
$(LIBCORE): $(PWD)/core $(LIBCOREOBJS)
	$(AR) cr $@ $(LIBCOREOBJS)
	$(RANLIB) $@
$(PWD)/core:
	$(MKDIR_P) $(PWD)/core

# libconfig
$(LIBCONFIG): $(PWD)/config $(LIBCONFIGOBJS)
	$(AR) cr $@ $(LIBCONFIGOBJS)
	$(RANLIB) $@
$(PWD)/config:
	$(MKDIR_P) $(PWD)/config

# libhcache
hcache/hcache.o:	hcache/hcversion.h
$(LIBHCACHE): $(PWD)/hcache $(LIBHCACHEOBJS)
	$(AR) cr $@ $(LIBHCACHEOBJS)
	$(RANLIB) $@
$(PWD)/hcache:
	$(MKDIR_P) $(PWD)/hcache

# libbcache
$(LIBBCACHE): $(PWD)/bcache $(LIBBCACHEOBJS)
	$(AR) cr $@ $(LIBBCACHEOBJS)
	$(RANLIB) $@
$(PWD)/bcache:
	$(MKDIR_P) $(PWD)/bcache

# libhistory
$(LIBHISTORY): $(PWD)/history $(LIBHISTORYOBJS)
	$(AR) cr $@ $(LIBHISTORYOBJS)
	$(RANLIB) $@
$(PWD)/history:
	$(MKDIR_P) $(PWD)/history

# pgpewrap
$(PGPEWRAP): $(PGPEWRAPOBJS)
	$(CC) $(LDFLAGS) -o $@ $(PGPEWRAPOBJS)

# generated
git_ver.c: $(ALL_FILES)
	version=`git describe --dirty --abbrev=6 --match "20[0-9][0-9][0-9][0-9][0-9][0-9]" 2> /dev/null | \
		sed -e 's/^[0-9]\{8\}//; s/-g\([a-z0-9]\{6\}\)/-\1/'`; \
	echo 'const char *GitVer = "'$$version'";' > $@.tmp; \
	cmp -s $@.tmp $@ || mv $@.tmp $@; \
	rm -f $@.tmp

hcache/hcversion.h:	$(SRCDIR)/address/address.h $(SRCDIR)/email/body.h \
			$(SRCDIR)/email/email.h $(SRCDIR)/email/envelope.h \
			$(SRCDIR)/email/parameter.h $(SRCDIR)/hcache/hcachever.sh \
			$(SRCDIR)/mutt/buffer.h $(SRCDIR)/mutt/list.h
	$(MKDIR_P) $(PWD)/hcache
	( echo '#include "config.h"'; \
	echo '#include "address/address.h"'; \
	echo '#include "email/body.h"'; \
	echo '#include "email/email.h"'; \
	echo '#include "email/envelope.h"'; \
	echo '#include "email/parameter.h"'; \
	echo '#include "mutt/buffer.h"'; \
	echo '#include "mutt/list.h"';) | $(CPP) $(CFLAGS) - | \
	sh $(SRCDIR)/hcache/hcachever.sh hcache/hcversion.h

# clean
clean: $(CLEAN_TARGETS)
	rm -f $(CLEANFILES)

# install
install: all $(INSTALL_TARGETS)
	# Install binaries
	$(MKDIR_P) $(DESTDIR)$(bindir)
	for f in $(BINFILES); do \
		$(INSTALL) $$f $(DESTDIR)$(bindir)/`basename $$f` || exit 1; \
	done
	$(MKDIR_P) $(DESTDIR)$(libexecdir)/$(PACKAGE)
	for f in $(LIBBINFILES); do \
		$(INSTALL) $$f $(DESTDIR)$(libexecdir)/$(PACKAGE)/`basename $$f` || exit 1; \
	done

# uninstall
uninstall: $(UNINSTALL_TARGETS)
	# Uninstall binaries
	for f in $(BINFILES); do \
		$(RM) $(DESTDIR)$(bindir)/`basename $$f`; \
	done
	for f in $(LIBBINFILES); do \
		$(RM) $(DESTDIR)$(libexecdir)/$(PACKAGE)/`basename $$f`; \
	done

# distclean
distclean: clean
	$(RM) $(DEPFILES) autosetup/jimsh0 config.h config.log conststrings.c \
		contrib/Makefile docs/Makefile docs/neomutt.1 html Makefile po/Makefile \
		test/Makefile .clang_complete
	$(RM) *.gc?? */*.gc?? test/*/*.gc??
	$(RM) lcov.info lcov


##############################################################################
# include generated dependency files
-include $(DEPFILES)

##############################################################################
# include special-purpose makefiles, each one of which MUST define the
# dedicated all-, clean-, install-, and uninstall- targets.

MSGID_BUGS_ADDRESS = neomutt-devel@neomutt.org
GMSGFMT		= msgfmt
localedir	= /usr/share/locale
XGETTEXT	= xgettext
XGETTEXT_OPTIONS= --keyword=_ --keyword=N_ --no-wrap --add-comments=L10N
MSGMERGE	= msgmerge
MSGMERGE_UPDATE = msgmerge --update --no-wrap

MOFILES		= po/bg.mo po/ca.mo po/cs.mo po/da.mo po/de.mo po/el.mo \
		po/en_GB.mo po/eo.mo po/es.mo po/et.mo po/eu.mo po/fi.mo po/fr.mo \
		po/ga.mo po/gl.mo po/hu.mo po/id.mo po/it.mo po/ja.mo po/ko.mo \
		po/lt.mo po/nl.mo po/pl.mo po/pt_BR.mo po/ru.mo po/sk.mo po/sv.mo \
		po/tr.mo po/uk.mo po/zh_CN.mo po/zh_TW.mo
POTFILE		= po/$(PACKAGE).pot

all-po:	po $(MOFILES)

.PHONY: po
po:
	$(MKDIR_P) po

.SUFFIXES: .mo .po
.po.mo:
	$(GMSGFMT) -c -o $@ $<

clean-po:
	$(RM) $(MOFILES)
	$(RM) $(POTFILE)
	$(RM) po/*~

install-po: all-po
	@catalogs='$(MOFILES)'; \
	for cat in $$catalogs; do \
	  lang=`echo $$cat | sed -e 's/\.mo$$//' -e 's|^po/||'`; \
	  dir=$(DESTDIR)$(localedir)/$$lang/LC_MESSAGES; \
	  $(MKDIR_P) $$dir; \
	  $(INSTALL_DATA) $$cat $$dir/$(PACKAGE).mo || exit 1; \
	  echo "installing $$dir/$(PACKAGE).mo"; \
	done

uninstall-po:
	@catalogs='$(MOFILES)'; \
	for cat in $$catalogs; do \
	  lang=`echo $$cat | sed -e 's/\.mo$$//' -e 's|^po/||'`; \
	  dir=$(DESTDIR)$(localedir)/$$lang/LC_MESSAGES; \
	  $(RM) $$dir/$(PACKAGE).mo; \
	  echo "removing $$dir/$(PACKAGE).mo"; \
	done

update-po: clean-po
	$(MAKE) $(PACKAGE).pot-update
	$(MAKE) all-po

$(PACKAGE).pot-update:
	$(XGETTEXT) --default-domain=$(PACKAGE) --directory=$(SRCDIR) \
	    $(XGETTEXT_OPTIONS) \
	    --files-from=po/POTFILES.in \
	    --copyright-holder='$(COPYRIGHT_HOLDER)' \
	    --package-name="$(PACKAGE)" \
	    --package-version="$(PACKAGE_VERSION)" \
	    --msgid-bugs-address="$(MSGID_BUGS_ADDRESS)" \
	    --output=$(POTFILE)
	@catalogs='$(MOFILES)'; \
	for cat in $$catalogs; do \
	  lang=`echo $$cat | sed -e 's/\.mo$$//'`; \
	  echo "$(MSGMERGE_UPDATE) $${lang}.po $(POTFILE)"; \
	  $(MSGMERGE_UPDATE) $${lang}.po $(POTFILE); \
	done


# vim: set ts=8 noexpandtab:
SAMPLES=	colors.default colors.linux gpg.rc mairix_filter.pl Mush.rc Pine.rc \
		sample.mailcap sample.neomuttrc sample.neomuttrc-starter \
		sample.neomuttrc-tlr smime.rc smime_keys_test.pl Tin.rc

CONTRIB_DIRS=	colorschemes hcache-bench keybase logo lua vim-keys

all-contrib:
clean-contrib:

install-contrib:
	$(INSTALL) -d -m 755 $(DESTDIR)$(docdir)/samples
	for f in $(SAMPLES); do \
		sed -e 's!/usr/libexec!$(libexecdir)!g' $(SRCDIR)/contrib/$$f > $(SRCDIR)/contrib/$$f.tmp; \
		$(INSTALL) -m 644 $(SRCDIR)/contrib/$$f.tmp $(DESTDIR)$(docdir)/samples/$$f || exit 1; \
		rm -f -- $(SRCDIR)/contrib/$$f.tmp; \
	done
	for d in $(CONTRIB_DIRS); do \
		echo "Creating directory $(DESTDIR)$(docdir)/$$d"; \
		$(INSTALL) -d -m 755 $(DESTDIR)$(docdir)/$$d || exit 1; \
		for f in $(SRCDIR)/contrib/$$d/*; do \
			echo "Installing $$f"; \
			$(INSTALL) -m 644 $$f $(DESTDIR)$(docdir)/$$d || exit 1; \
		done \
	done
	chmod +x $(DESTDIR)$(docdir)/keybase/*.sh

uninstall-contrib:
	for f in $(SAMPLES); do \
		$(RM) $(DESTDIR)$(docdir)/samples/$$f; \
	done
	for d in $(CONTRIB_DIRS); do \
		$(RM) $(DESTDIR)$(docdir)/$$d; \
	done
	-rmdir $(DESTDIR)$(docdir)/samples
	-rmdir $(DESTDIR)$(docdir)

# vim: set ts=8 noexpandtab:
MAKEDOC_CPP = $(CC_FOR_BUILD) $(CFLAGS_FOR_BUILD) -D_MAKEDOC -E -C -I.

docs/makedoc$(EXEEXT): $(SRCDIR)/docs/makedoc.c
	$(CC_FOR_BUILD) -I. $(CFLAGS_FOR_BUILD) $(LDFLAGS_FOR_BUILD) \
	    -o $@ $(SRCDIR)/docs/makedoc.c

docs/neomuttrc: docs/makedoc$(EXEEXT) $(SRCDIR)/docs/neomuttrc.head \
		$(SRCDIR)/mutt_config.c
	sed -e 's,@docdir@,$(docdir),' $(SRCDIR)/docs/neomuttrc.head \
	    > docs/neomuttrc
	$(MAKEDOC_CPP) $(SRCDIR)/mutt_config.c | docs/makedoc$(EXEEXT) -c \
	    >> docs/neomuttrc

.PHONY: docs
docs:
	$(MKDIR_P) docs

# Let's generate neomuttrc in all cases: it doesn't require any additional 3rd
# party dependencies and distributions tend to rely on having it.
all-docs: docs docs/neomuttrc

clean-docs:
	$(RM) docs/neomuttrc docs/makedoc$(EXEEXT)

install-docs: all-docs
	$(MKDIR_P) $(DESTDIR)$(sysconfdir)
	$(INSTALL) -m 644 docs/neomuttrc $(DESTDIR)$(sysconfdir)/neomuttrc

uninstall-docs:
	$(RM) $(DESTDIR)$(sysconfdir)/neomuttrc

# vim: set ts=8 noexpandtab:
ACCOUNT_OBJS	= test/account/account_free.o \
		  test/account/account_mailbox_add.o \
		  test/account/account_mailbox_remove.o \
		  test/account/account_new.o

ADDRESS_OBJS	= test/address/mutt_addr_cat.o \
		  test/address/mutt_addr_cmp.o \
		  test/address/mutt_addr_copy.o \
		  test/address/mutt_addr_create.o \
		  test/address/mutt_addr_for_display.o \
		  test/address/mutt_addr_free.o \
		  test/address/mutt_addr_new.o \
		  test/address/mutt_addr_to_intl.o \
		  test/address/mutt_addr_to_local.o \
		  test/address/mutt_addr_valid_msgid.o \
		  test/address/mutt_addr_write.o \
		  test/address/mutt_addrlist_append.o \
		  test/address/mutt_addrlist_clear.o \
		  test/address/mutt_addrlist_copy.o \
		  test/address/mutt_addrlist_count_recips.o \
		  test/address/mutt_addrlist_dedupe.o \
		  test/address/mutt_addrlist_equal.o \
		  test/address/mutt_addrlist_parse.o \
		  test/address/mutt_addrlist_parse2.o \
		  test/address/mutt_addrlist_prepend.o \
		  test/address/mutt_addrlist_qualify.o \
		  test/address/mutt_addrlist_remove.o \
		  test/address/mutt_addrlist_remove_xrefs.o \
		  test/address/mutt_addrlist_search.o \
		  test/address/mutt_addrlist_to_intl.o \
		  test/address/mutt_addrlist_to_local.o \
		  test/address/mutt_addrlist_write.o

ATTACH_OBJS	= test/attach/mutt_actx_add_attach.o \
		  test/attach/mutt_actx_add_body.o \
		  test/attach/mutt_actx_add_fp.o \
		  test/attach/mutt_actx_entries_free.o \
		  test/attach/mutt_actx_free.o \
		  test/attach/mutt_actx_new.o

BASE64_OBJS	= test/base64/mutt_b64_buffer_decode.o \
		  test/base64/mutt_b64_buffer_encode.o \
		  test/base64/mutt_b64_decode.o \
		  test/base64/mutt_b64_encode.o

BODY_OBJS	= test/body/mutt_body_cmp_strict.o \
		  test/body/mutt_body_free.o \
		  test/body/mutt_body_new.o

BUFFER_OBJS	= test/buffer/mutt_buffer_add_printf.o \
		  test/buffer/mutt_buffer_addch.o \
		  test/buffer/mutt_buffer_addstr.o \
		  test/buffer/mutt_buffer_addstr_n.o \
		  test/buffer/mutt_buffer_alloc.o \
		  test/buffer/mutt_buffer_concat_path.o \
		  test/buffer/mutt_buffer_concatn_path.o \
		  test/buffer/mutt_buffer_copy.o \
		  test/buffer/mutt_buffer_dealloc.o \
		  test/buffer/mutt_buffer_fix_dptr.o \
		  test/buffer/mutt_buffer_init.o \
		  test/buffer/mutt_buffer_is_empty.o \
		  test/buffer/mutt_buffer_len.o \
		  test/buffer/mutt_buffer_make.o \
		  test/buffer/mutt_buffer_pool_free.o \
		  test/buffer/mutt_buffer_pool_get.o \
		  test/buffer/mutt_buffer_pool_release.o \
		  test/buffer/mutt_buffer_printf.o \
		  test/buffer/mutt_buffer_reset.o \
		  test/buffer/mutt_buffer_strcpy.o \
		  test/buffer/mutt_buffer_strcpy_n.o \
		  test/buffer/mutt_buffer_strdup.o \
		  test/buffer/mutt_buffer_substrcpy.o

CHARSET_OBJS	= test/charset/mutt_ch_canonical_charset.o \
		  test/charset/mutt_ch_charset_lookup.o \
		  test/charset/mutt_ch_check.o \
		  test/charset/mutt_ch_check_charset.o \
		  test/charset/mutt_ch_choose.o \
		  test/charset/mutt_ch_chscmp.o \
		  test/charset/mutt_ch_convert_nonmime_string.o \
		  test/charset/mutt_ch_convert_string.o \
		  test/charset/mutt_ch_fgetconv.o \
		  test/charset/mutt_ch_fgetconv_close.o \
		  test/charset/mutt_ch_fgetconv_open.o \
		  test/charset/mutt_ch_fgetconvs.o \
		  test/charset/mutt_ch_get_default_charset.o \
		  test/charset/mutt_ch_get_langinfo_charset.o \
		  test/charset/mutt_ch_iconv.o \
		  test/charset/mutt_ch_iconv_lookup.o \
		  test/charset/mutt_ch_iconv_open.o \
		  test/charset/mutt_ch_lookup_add.o \
		  test/charset/mutt_ch_lookup_remove.o \
		  test/charset/mutt_ch_set_charset.o

CONFIG_OBJS	= test/config/account.o \
		  test/config/address.o \
		  test/config/bool.o \
		  test/config/common.o \
		  test/config/dump.o \
		  test/config/enum.o \
		  test/config/initial.o \
		  test/config/long.o \
		  test/config/mbtable.o \
		  test/config/number.o \
		  test/config/path.o \
		  test/config/quad.o \
		  test/config/regex.o \
		  test/config/set.o \
		  test/config/slist.o \
		  test/config/sort.o \
		  test/config/string.o \
		  test/config/subset.o \
		  test/config/synonym.o

DATE_OBJS	= test/date/mutt_date_add_timeout.o \
		  test/date/mutt_date_check_month.o \
		  test/date/mutt_date_epoch.o \
		  test/date/mutt_date_epoch_ms.o \
		  test/date/mutt_date_gmtime.o \
		  test/date/mutt_date_is_day_name.o \
		  test/date/mutt_date_local_tz.o \
		  test/date/mutt_date_localtime.o \
		  test/date/mutt_date_localtime_format.o \
		  test/date/mutt_date_make_date.o \
		  test/date/mutt_date_make_imap.o \
		  test/date/mutt_date_make_time.o \
		  test/date/mutt_date_make_tls.o \
		  test/date/mutt_date_normalize_time.o \
		  test/date/mutt_date_parse_date.o \
		  test/date/mutt_date_parse_imap.o \
		  test/date/mutt_date_sleep_ms.o

EMAIL_OBJS	= test/email/common.o \
		  test/email/email_cmp_strict.o \
		  test/email/email_free.o \
		  test/email/email_new.o \
		  test/email/email_size.o \
		  test/email/emaillist_add_email.o \
		  test/email/emaillist_clear.o \
		  test/email/mutt_autocrypthdr_new.o

ENVELOPE_OBJS	= test/envelope/mutt_autocrypthdr_free.o \
		  test/envelope/mutt_env_cmp_strict.o \
		  test/envelope/mutt_env_free.o \
		  test/envelope/mutt_env_merge.o \
		  test/envelope/mutt_env_new.o \
		  test/envelope/mutt_env_to_intl.o \
		  test/envelope/mutt_env_to_local.o

ENVLIST_OBJS	= test/envlist/mutt_envlist_free.o \
		  test/envlist/mutt_envlist_getlist.o \
		  test/envlist/mutt_envlist_init.o \
		  test/envlist/mutt_envlist_set.o \
		  test/envlist/mutt_envlist_unset.o

FILE_OBJS	= test/file/common.o \
		  test/file/mutt_buffer_file_expand_fmt_quote.o \
		  test/file/mutt_buffer_quote_filename.o \
		  test/file/mutt_file_check_empty.o \
		  test/file/mutt_file_chmod.o \
		  test/file/mutt_file_chmod_add.o \
		  test/file/mutt_file_chmod_add_stat.o \
		  test/file/mutt_file_chmod_rm.o \
		  test/file/mutt_file_chmod_rm_stat.o \
		  test/file/mutt_file_copy_bytes.o \
		  test/file/mutt_file_copy_stream.o \
		  test/file/mutt_file_decrease_mtime.o \
		  test/file/mutt_file_expand_fmt.o \
		  test/file/mutt_file_fclose.o \
		  test/file/mutt_file_fopen.o \
		  test/file/mutt_file_fsync_close.o \
		  test/file/mutt_file_get_size.o \
		  test/file/mutt_file_get_stat_timespec.o \
		  test/file/mutt_file_iter_line.o \
		  test/file/mutt_file_lock.o \
		  test/file/mutt_file_map_lines.o \
		  test/file/mutt_file_mkdir.o \
		  test/file/mutt_file_mkstemp_full.o \
		  test/file/mutt_file_open.o \
		  test/file/mutt_file_quote_filename.o \
		  test/file/mutt_file_read_keyword.o \
		  test/file/mutt_file_read_line.o \
		  test/file/mutt_file_rename.o \
		  test/file/mutt_file_resolve_symlink.o \
		  test/file/mutt_file_rmtree.o \
		  test/file/mutt_file_safe_rename.o \
		  test/file/mutt_file_sanitize_filename.o \
		  test/file/mutt_file_sanitize_regex.o \
		  test/file/mutt_file_set_mtime.o \
		  test/file/mutt_file_stat_compare.o \
		  test/file/mutt_file_stat_timespec_compare.o \
		  test/file/mutt_file_symlink.o \
		  test/file/mutt_file_timespec_compare.o \
		  test/file/mutt_file_touch_atime.o \
		  test/file/mutt_file_unlink.o \
		  test/file/mutt_file_unlink_empty.o \
		  test/file/mutt_file_unlock.o

FILTER_OBJS	= test/filter/filter_create.o \
		  test/filter/filter_create_fd.o \
		  test/filter/filter_wait.o

FROM_OBJS	= test/from/is_from.o

GROUP_OBJS	= test/group/mutt_group_match.o \
		  test/group/mutt_grouplist_add.o \
		  test/group/mutt_grouplist_add_addrlist.o \
		  test/group/mutt_grouplist_add_regex.o \
		  test/group/mutt_grouplist_clear.o \
		  test/group/mutt_grouplist_destroy.o \
		  test/group/mutt_grouplist_free.o \
		  test/group/mutt_grouplist_init.o \
		  test/group/mutt_grouplist_remove_addrlist.o \
		  test/group/mutt_grouplist_remove_regex.o \
		  test/group/mutt_pattern_group.o

GUI_OBJS	= test/gui/reflow.o

HASH_OBJS	= test/hash/mutt_hash_delete.o \
		  test/hash/mutt_hash_find.o \
		  test/hash/mutt_hash_find_bucket.o \
		  test/hash/mutt_hash_find_elem.o \
		  test/hash/mutt_hash_free.o \
		  test/hash/mutt_hash_insert.o \
		  test/hash/mutt_hash_int_delete.o \
		  test/hash/mutt_hash_int_find.o \
		  test/hash/mutt_hash_int_insert.o \
		  test/hash/mutt_hash_int_new.o \
		  test/hash/mutt_hash_new.o \
		  test/hash/mutt_hash_set_destructor.o \
		  test/hash/mutt_hash_typed_insert.o \
		  test/hash/mutt_hash_walk.o

HISTORY_OBJS	= test/history/mutt_hist_add.o \
		  test/history/mutt_hist_at_scratch.o \
		  test/history/mutt_hist_free.o \
		  test/history/mutt_hist_init.o \
		  test/history/mutt_hist_next.o \
		  test/history/mutt_hist_prev.o \
		  test/history/mutt_hist_read_file.o \
		  test/history/mutt_hist_reset_state.o \
		  test/history/mutt_hist_save_scratch.o \
		  test/history/mutt_hist_search.o

IDNA_OBJS	= test/idna/mutt_idna_intl_to_local.o \
		  test/idna/mutt_idna_local_to_intl.o \
		  test/idna/mutt_idna_print_version.o \
		  test/idna/mutt_idna_to_ascii_lz.o

LIST_OBJS	= test/list/common.o \
		  test/list/mutt_list_clear.o \
		  test/list/mutt_list_compare.o \
		  test/list/mutt_list_find.o \
		  test/list/mutt_list_free.o \
		  test/list/mutt_list_free_type.o \
		  test/list/mutt_list_insert_after.o \
		  test/list/mutt_list_insert_head.o \
		  test/list/mutt_list_insert_tail.o \
		  test/list/mutt_list_match.o \
		  test/list/mutt_list_str_split.o

LOGGING_OBJS	= test/logging/log_disp_file.o \
		  test/logging/log_disp_null.o \
		  test/logging/log_disp_queue.o \
		  test/logging/log_disp_terminal.o \
		  test/logging/log_file_close.o \
		  test/logging/log_file_open.o \
		  test/logging/log_file_running.o \
		  test/logging/log_file_set_filename.o \
		  test/logging/log_file_set_level.o \
		  test/logging/log_file_set_version.o \
		  test/logging/log_queue_add.o \
		  test/logging/log_queue_empty.o \
		  test/logging/log_queue_flush.o \
		  test/logging/log_queue_save.o \
		  test/logging/log_queue_set_max_size.o

MAILBOX_OBJS	= test/mailbox/mailbox_changed.o \
		  test/mailbox/mailbox_find.o \
		  test/mailbox/mailbox_find_name.o \
		  test/mailbox/mailbox_free.o \
		  test/mailbox/mailbox_new.o \
		  test/mailbox/mailbox_set_subset.o \
		  test/mailbox/mailbox_size_add.o \
		  test/mailbox/mailbox_size_sub.o \
		  test/mailbox/mailbox_update.o

MAPPING_OBJS	= test/mapping/mutt_map_get_name.o \
		  test/mapping/mutt_map_get_value.o

MBYTE_OBJS	= test/mbyte/mutt_mb_charlen.o \
		  test/mbyte/mutt_mb_filter_unprintable.o \
		  test/mbyte/mutt_mb_get_initials.o \
		  test/mbyte/mutt_mb_is_display_corrupting_utf8.o \
		  test/mbyte/mutt_mb_is_lower.o \
		  test/mbyte/mutt_mb_is_shell_char.o \
		  test/mbyte/mutt_mb_mbstowcs.o \
		  test/mbyte/mutt_mb_wcstombs.o \
		  test/mbyte/mutt_mb_wcswidth.o \
		  test/mbyte/mutt_mb_wcwidth.o \
		  test/mbyte/mutt_mb_width.o \
		  test/mbyte/mutt_mb_width_ceiling.o

MD5_OBJS	= test/md5/common.o \
		  test/md5/mutt_md5.o \
		  test/md5/mutt_md5_bytes.o \
		  test/md5/mutt_md5_finish_ctx.o \
		  test/md5/mutt_md5_init_ctx.o \
		  test/md5/mutt_md5_process.o \
		  test/md5/mutt_md5_process_bytes.o \
		  test/md5/mutt_md5_toascii.o

MEMORY_OBJS	= test/memory/mutt_mem_calloc.o \
		  test/memory/mutt_mem_free.o \
		  test/memory/mutt_mem_malloc.o \
		  test/memory/mutt_mem_realloc.o

NEOMUTT_OBJS	= test/neo/neomutt_account_add.o \
		  test/neo/neomutt_account_remove.o \
		  test/neo/neomutt_free.o \
		  test/neo/neomutt_mailboxlist_clear.o \
		  test/neo/neomutt_mailboxlist_get_all.o \
		  test/neo/neomutt_new.o

NOTIFY_OBJS	= test/notify/notify_free.o \
		  test/notify/notify_new.o \
		  test/notify/notify_observer_add.o \
		  test/notify/notify_observer_remove.o \
		  test/notify/notify_send.o \
		  test/notify/notify_set_parent.o

PARAMETER_OBJS	= test/parameter/mutt_param_cmp_strict.o \
		  test/parameter/mutt_param_delete.o \
		  test/parameter/mutt_param_free.o \
		  test/parameter/mutt_param_free_one.o \
		  test/parameter/mutt_param_get.o \
		  test/parameter/mutt_param_new.o \
		  test/parameter/mutt_param_set.o

PARSE_OBJS	= test/parse/mutt_auto_subscribe.o \
		  test/parse/mutt_check_encoding.o \
		  test/parse/mutt_check_mime_type.o \
		  test/parse/mutt_extract_message_id.o \
		  test/parse/mutt_is_message_type.o \
		  test/parse/mutt_matches_ignore.o \
		  test/parse/mutt_parse_content_type.o \
		  test/parse/mutt_parse_mailto.o \
		  test/parse/mutt_parse_multipart.o \
		  test/parse/mutt_parse_part.o \
		  test/parse/mutt_read_mime_header.o \
		  test/parse/mutt_rfc822_parse_line.o \
		  test/parse/mutt_rfc822_parse_message.o \
		  test/parse/mutt_rfc822_read_header.o \
		  test/parse/mutt_rfc822_read_line.o

PATH_OBJS	= test/path/mutt_path_abbr_folder.o \
		  test/path/mutt_path_basename.o \
		  test/path/mutt_path_canon.o \
		  test/path/mutt_path_concat.o \
		  test/path/mutt_path_dirname.o \
		  test/path/mutt_path_escape.o \
		  test/path/mutt_path_getcwd.o \
		  test/path/mutt_path_parent.o \
		  test/path/mutt_path_pretty.o \
		  test/path/mutt_path_realpath.o \
		  test/path/mutt_path_tidy.o \
		  test/path/mutt_path_tidy_dotdot.o \
		  test/path/mutt_path_tidy_slash.o \
		  test/path/mutt_path_tilde.o \
		  test/path/mutt_path_to_absolute.o

PATTERN_OBJS	= pattern.o \
		  test/pattern/comp.o \
		  test/pattern/dummy.o \
		  test/pattern/extract.o

REGEX_OBJS	= test/regex/mutt_regex_capture.o \
		  test/regex/mutt_regex_compile.o \
		  test/regex/mutt_regex_free.o \
		  test/regex/mutt_regex_match.o \
		  test/regex/mutt_regex_new.o \
		  test/regex/mutt_regexlist_add.o \
		  test/regex/mutt_regexlist_free.o \
		  test/regex/mutt_regexlist_match.o \
		  test/regex/mutt_regexlist_new.o \
		  test/regex/mutt_regexlist_remove.o \
		  test/regex/mutt_replacelist_add.o \
		  test/regex/mutt_replacelist_apply.o \
		  test/regex/mutt_replacelist_free.o \
		  test/regex/mutt_replacelist_match.o \
		  test/regex/mutt_replacelist_new.o \
		  test/regex/mutt_replacelist_remove.o

RFC2047_OBJS	= test/rfc2047/common.o \
		  test/rfc2047/rfc2047_decode.o \
		  test/rfc2047/rfc2047_decode_addrlist.o \
		  test/rfc2047/rfc2047_decode_envelope.o \
		  test/rfc2047/rfc2047_encode.o \
		  test/rfc2047/rfc2047_encode_addrlist.o \
		  test/rfc2047/rfc2047_encode_envelope.o

RFC2231_OBJS	= test/rfc2231/rfc2231_decode_parameters.o \
		  test/rfc2231/rfc2231_encode_string.o

SIGNAL_OBJS	= test/signal/mutt_sig_allow_interrupt.o \
		  test/signal/mutt_sig_block.o \
		  test/signal/mutt_sig_block_system.o \
		  test/signal/mutt_sig_empty_handler.o \
		  test/signal/mutt_sig_exit_handler.o \
		  test/signal/mutt_sig_init.o \
		  test/signal/mutt_sig_unblock.o \
		  test/signal/mutt_sig_unblock_system.o

SLIST_OBJS	= test/slist/slist_add_list.o \
		  test/slist/slist_add_string.o \
		  test/slist/slist_compare.o \
		  test/slist/slist_dup.o \
		  test/slist/slist_empty.o \
		  test/slist/slist_free.o \
		  test/slist/slist_is_member.o \
		  test/slist/slist_parse.o \
		  test/slist/slist_remove_string.o

STRING_OBJS	= test/string/mutt_str_adjust.o \
		  test/string/mutt_str_append_item.o \
		  test/string/mutt_str_asprintf.o \
		  test/string/mutt_str_atoi.o \
		  test/string/mutt_str_atol.o \
		  test/string/mutt_str_atos.o \
		  test/string/mutt_str_atoui.o \
		  test/string/mutt_str_atoul.o \
		  test/string/mutt_str_atoull.o \
		  test/string/mutt_str_dequote_comment.o \
		  test/string/mutt_str_find_word.o \
		  test/string/mutt_str_getenv.o \
		  test/string/mutt_str_inline_replace.o \
		  test/string/mutt_str_is_ascii.o \
		  test/string/mutt_str_is_email_wsp.o \
		  test/string/mutt_str_lws_len.o \
		  test/string/mutt_str_lws_rlen.o \
		  test/string/mutt_str_next_word.o \
		  test/string/mutt_str_remall_strcasestr.o \
		  test/string/mutt_str_remove_trailing_ws.o \
		  test/string/mutt_str_replace.o \
		  test/string/mutt_str_rstrnstr.o \
		  test/string/mutt_str_skip_email_wsp.o \
		  test/string/mutt_str_skip_whitespace.o \
		  test/string/mutt_str_startswith.o \
		  test/string/mutt_str_strcasecmp.o \
		  test/string/mutt_str_strcasestr.o \
		  test/string/mutt_str_strcat.o \
		  test/string/mutt_str_strchrnul.o \
		  test/string/mutt_str_strcmp.o \
		  test/string/mutt_str_strcoll.o \
		  test/string/mutt_str_strdup.o \
		  test/string/mutt_str_strfcpy.o \
		  test/string/mutt_str_stristr.o \
		  test/string/mutt_str_strlen.o \
		  test/string/mutt_str_strlower.o \
		  test/string/mutt_str_strncasecmp.o \
		  test/string/mutt_str_strncat.o \
		  test/string/mutt_str_strncmp.o \
		  test/string/mutt_str_strnfcpy.o \
		  test/string/mutt_str_strnlower.o \
		  test/string/mutt_str_substr_copy.o \
		  test/string/mutt_str_substr_dup.o \
		  test/string/mutt_str_sysexit.o \
		  test/string/mutt_str_word_casecmp.o

TAGS_OBJS	= test/tags/driver_tags_free.o \
		  test/tags/driver_tags_get.o \
		  test/tags/driver_tags_get_transformed.o \
		  test/tags/driver_tags_get_transformed_for.o \
		  test/tags/driver_tags_get_with_hidden.o \
		  test/tags/driver_tags_replace.o

THREAD_OBJS	= test/thread/clean_references.o \
		  test/thread/find_virtual.o \
		  test/thread/insert_message.o \
		  test/thread/is_descendant.o \
		  test/thread/mutt_break_thread.o \
		  test/thread/thread_hash_destructor.o \
		  test/thread/unlink_message.o

URL_OBJS	= test/url/url_check_scheme.o \
		  test/url/url_free.o \
		  test/url/url_parse.o \
		  test/url/url_pct_decode.o \
		  test/url/url_pct_encode.o \
		  test/url/url_tobuffer.o \
		  test/url/url_tostring.o

BUILD_DIRS	= $(PWD)/test/account $(PWD)/test/address $(PWD)/test/attach \
		  $(PWD)/test/base64 $(PWD)/test/body $(PWD)/test/buffer \
		  $(PWD)/test/charset $(PWD)/test/config $(PWD)/test/date \
		  $(PWD)/test/email $(PWD)/test/envelope $(PWD)/test/envlist \
		  $(PWD)/test/file $(PWD)/test/filter $(PWD)/test/from \
		  $(PWD)/test/group $(PWD)/test/gui $(PWD)/test/hash \
		  $(PWD)/test/history $(PWD)/test/idna $(PWD)/test/list \
		  $(PWD)/test/logging $(PWD)/test/mailbox $(PWD)/test/mapping \
		  $(PWD)/test/mbyte $(PWD)/test/md5 $(PWD)/test/memory \
		  $(PWD)/test/neo $(PWD)/test/notify $(PWD)/test/parameter \
		  $(PWD)/test/parse $(PWD)/test/path $(PWD)/test/pattern \
		  $(PWD)/test/regex $(PWD)/test/rfc2047 $(PWD)/test/rfc2231 \
		  $(PWD)/test/signal $(PWD)/test/slist $(PWD)/test/string \
		  $(PWD)/test/tags $(PWD)/test/thread $(PWD)/test/url

TEST_OBJS	= test/main.o test/common.o \
		  $(ACCOUNT_OBJS) \
		  $(ADDRESS_OBJS) \
		  $(ATTACH_OBJS) \
		  $(BASE64_OBJS) \
		  $(BODY_OBJS) \
		  $(BUFFER_OBJS) \
		  $(CHARSET_OBJS) \
		  $(CONFIG_OBJS) \
		  $(DATE_OBJS) \
		  $(EMAIL_OBJS) \
		  $(ENVELOPE_OBJS) \
		  $(ENVLIST_OBJS) \
		  $(FILE_OBJS) \
		  $(FILTER_OBJS) \
		  $(FROM_OBJS) \
		  $(GROUP_OBJS) \
		  $(GUI_OBJS) \
		  $(HASH_OBJS) \
		  $(HISTORY_OBJS) \
		  $(IDNA_OBJS) \
		  $(LIST_OBJS) \
		  $(LOGGING_OBJS) \
		  $(MAILBOX_OBJS) \
		  $(MAPPING_OBJS) \
		  $(MBYTE_OBJS) \
		  $(MD5_OBJS) \
		  $(MEMORY_OBJS) \
		  $(NEOMUTT_OBJS) \
		  $(NOTIFY_OBJS) \
		  $(PARAMETER_OBJS) \
		  $(PARSE_OBJS) \
		  $(PATH_OBJS) \
		  $(PATTERN_OBJS) \
		  $(REGEX_OBJS) \
		  $(RFC2047_OBJS) \
		  $(RFC2231_OBJS) \
		  $(SIGNAL_OBJS) \
		  $(SLIST_OBJS) \
		  $(STRING_OBJS) \
		  $(TAGS_OBJS) \
		  $(THREAD_OBJS) \
		  $(URL_OBJS)

CFLAGS	+= -I$(SRCDIR)/test

TEST_BINARY = test/neomutt-test$(EXEEXT)

.PHONY: test
test: $(TEST_BINARY)
	$(TEST_BINARY)

$(BUILD_DIRS):
	$(MKDIR_P) $@

$(TEST_BINARY): $(BUILD_DIRS) $(MUTTLIBS) $(TEST_OBJS)
	$(CC) -o $@ $(TEST_OBJS) $(MUTTLIBS) $(LDFLAGS) $(LIBS)

all-test: $(TEST_BINARY)

clean-test:
	$(RM) $(TEST_BINARY) $(TEST_OBJS) $(TEST_OBJS:.o=.Po)

install-test:
uninstall-test:

TEST_DEPFILES = $(TEST_OBJS:.o=.Po)
-include $(TEST_DEPFILES)

# vim: set ts=8 noexpandtab:

# vim: set ts=8 noexpandtab:
