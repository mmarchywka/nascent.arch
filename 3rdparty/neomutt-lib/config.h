#ifndef _CONFIG_H
#define _CONFIG_H
#define ALL_TARGETS "all-po all-docs all-contrib"
#define BINDIR "/usr/bin"
#define CLEAN_TARGETS "clean-po clean-docs clean-contrib"
#define CRYPT_BACKEND_CLASSIC_PGP 1
#define CRYPT_BACKEND_CLASSIC_SMIME 1
#define CRYPT_BACKEND_GPGME 1
#define ENABLE_NLS 1
#define HAVE_BIND_TEXTDOMAIN_CODESET 1
#define HAVE_BKGDSET 1
#define HAVE_CFLAG_STD_C99 1
#define HAVE_CLOCK_GETTIME 1
#define HAVE_COLOR 1
#define HAVE_CURS_SET 1
#define HAVE_DECL_GNUTLS_VERIFY_DISABLE_TIME_CHECKS 1
#define HAVE_DECL_SYS_SIGLIST 1
#define HAVE_FGETC_UNLOCKED 1
#define HAVE_FUTIMENS 1
#define HAVE_GETADDRINFO 1
#define HAVE_GETADDRINFO_A 1
#define HAVE_GETHOSTENT 1
#define HAVE_GETSID 1
#define HAVE_GNUTLS_CERTIFICATE_CREDENTIALS_T 1
#define HAVE_GNUTLS_CERTIFICATE_STATUS_T 1
#define HAVE_GNUTLS_CHECK_VERSION 1
#define HAVE_GNUTLS_DATUM_T 1
#define HAVE_GNUTLS_DIGEST_ALGORITHM_T 1
#define HAVE_GNUTLS_GNUTLS_H 1
#define HAVE_GNUTLS_PRIORITY_SET_DIRECT 1
#define HAVE_GNUTLS_SESSION_T 1
#define HAVE_GNUTLS_TRANSPORT_PTR_T 1
#define HAVE_GNUTLS_X509_CRT_T 1
#define HAVE_GNUTLS_X509_H 1
#define HAVE_GPGME 1
#define HAVE_GPGME_OP_EXPORT_KEYS 1
#define HAVE_ICONV 1
#define HAVE_ICONV_H 1
#define HAVE_IDNA_H 1
#define HAVE_IDNA_TO_ASCII_8Z 1
/* #undef HAVE_IDNA_TO_ASCII_FROM_LOCALE */
/* #undef HAVE_IDNA_TO_ASCII_FROM_UTF8 */
#define HAVE_IDNA_TO_ASCII_LZ 1
#define HAVE_IDNA_TO_UNICODE_8Z8Z 1
/* #undef HAVE_IDNA_TO_UNICODE_UTF8_FROM_UTF8 */
#define HAVE_INOTIFY_ADD_WATCH 1
#define HAVE_INOTIFY_INIT 1
#define HAVE_INOTIFY_INIT1 1
#define HAVE_INOTIFY_RM_WATCH 1
/* #undef HAVE_IOCTL_H */
#define HAVE_ISWBLANK 1
#define HAVE_LFS 1
#define HAVE_LIBIDN 1
#define HAVE_LITTLE_ENDIAN 1
#define HAVE_META 1
/* #undef HAVE_MINIX_CONFIG_H */
#define HAVE_MKDTEMP 1
#define HAVE_NCURSESW_NCURSES_H 1
#define HAVE_PGP 1
#define HAVE_SASL 1
#define HAVE_SASL_CLIENT_DONE 1
#define HAVE_SASL_ENCODE64 1
#define HAVE_SASL_SASL_H 1
#define HAVE_SETSOCKOPT 1
#define HAVE_SIGNAL_H 1
#define HAVE_SMIME 1
#define HAVE_START_COLOR 1
#define HAVE_STDLIB_H 1
#define HAVE_STRINGPREP_CHECK_VERSION 1
#define HAVE_STRINGPREP_H 1
#define HAVE_STRSEP 1
#define HAVE_STRUCT_STAT_ST_ATIM_TV_NSEC 1
#define HAVE_STRUCT_TIMESPEC 1
#define HAVE_SYSCALL_H 1
#define HAVE_SYSEXITS_H 1
#define HAVE_SYS_INOTIFY_H 1
#define HAVE_SYS_IOCTL_H 1
#define HAVE_SYS_PARAM_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_SYS_SYSCALL_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_TC 1
#define HAVE_TCBDBOPEN 1
#define HAVE_TCBDB_H 1
#define HAVE_TGETENT 1
#define HAVE_TIME_H 1
#define HAVE_TYPEAHEAD 1
#define HAVE_UNISTD_H 1
#define HAVE_USE_DEFAULT_COLORS 1
#define HAVE_USE_EXTENDED_NAMES 1
/* #undef HAVE_UTIMESNSAT */
#define HAVE_VASPRINTF 1
#define HAVE_WADDNWSTR 1
#define HAVE_WCSCASECMP 1
#define ICONV_CONST 
#define INSTALL_TARGETS "install-po install-docs install-contrib"
#define LOFF_T off_t
#define MAILPATH "/var/mail"
#define MUTTLOCALEDIR "/usr/share/locale"
#define OFF_T_FMT "%" PRId64
#define PACKAGE "neomutt"
#define PACKAGE_VERSION "20200501"
#define PKGDATADIR "/usr/share/neomutt"
#define PKGDOCDIR "/usr/share/doc/neomutt"
#define SENDMAIL "/usr/sbin/sendmail"
#define SIG_ATOMIC_VOLATILE_T volatile sig_atomic_t
#define SIZEOF_OFF_T 8
#define SUN_ATTACHMENT 1
#define SYSCONFDIR "/etc"
#define TEST_CASE_MAXSIZE 1024
#define UNINSTALL_TARGETS "uninstall-po uninstall-docs uninstall-contrib"
#define USE_COMP_MBOX 1
#define USE_FCNTL 1
/* #undef USE_FMEMOPEN */
#define USE_GSS 1
// mjm kluge 
//#define USE_HCACHE 1
#define USE_IMAP 1
#define USE_INOTIFY 1
#define USE_NNTP 1
#define USE_POP 1
#define USE_SASL 1
#define USE_SIDEBAR 1
#define USE_SMTP 1
#define USE_SOCKET 1
#define USE_SSL 1
#define USE_SSL_GNUTLS 1
#define VPATH "$(SRCDIR):$(SRCDIR)/po:$(SRCDIR)/docs:$(SRCDIR)/contrib"
#endif
