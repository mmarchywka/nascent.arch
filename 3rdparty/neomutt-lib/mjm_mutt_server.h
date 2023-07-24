#ifndef MJM_SERVER_H__
#define MJM_SERVER_H__


#include <stdbool.h>
#include <stdio.h>

/*

#include "mutt/observer.h" 
#include "mutt/notify_type.h" 
#include "mutt/string2.h" 
#include "config/subset.h" 

#include "globals.h"
#include "config.h"
#include <stdint.h>
#include "mutt/lib.h"
#include "mutt.h"
#include "mx.h"
#include "context.h"
#include "email/email_globals.h"
#include "core/mailbox.h"
#include "email/envelope.h"
#include "email/email.h"
#include "email/lib.h"
#include "imap/lib.h"
#include "send.h"

*/

#ifdef __cplusplus
extern "C" {
#endif
//#include "config.h"
#include "email/email.h"

//./email/email_globals.h:extern struct Regex *C_ReplyRegex;
extern struct Regex *C_ReplyRegex;
//extern struct NotifyCallback; //  *C_ReplyRegex;
extern struct Regex *C_ReplyRegex;
extern int  mjm_input_loop();
extern int  shutdown_mutt(); //   cls.rc = 0;
extern int  mutt_set_init_only(bool x); //   cls.rc = 0;
extern struct Email ** mutt_email_list(int flags);

extern int mjm_mutt_move_msg(struct Email * e, const char * dest, const int flags );

extern int set_read_flag(struct Email * e);
extern int reset_read_flag(struct Email * e);

/* message.c */

//int imap_copy_messages(struct Mailbox *m, struct EmailList *el, const char *dest, bool delete_original);

//  struct EmailList el = STAILQ_HEAD_INITIALIZER(el);
 // emaillist_add_email(&el, e);



extern  bool mjm_daemon_mode;
extern  bool mjm_init_only;
extern  bool mjm_debug_obs;
//int             mutt_sb_observer        (struct NotifyCallback *nc);

#define QFCKBOOST(x) #x
#ifndef QUOTE
#define QUOTE(x) QFCKBOOST(x)
#endif

#define MJM_CB_SIGS(x) int x  (struct NotifyCallback *nc)
#define MJM_CB_SIG(x) int x  (struct NotifyCallback *nc) { const char * nm=QUOTE(x);
#define MJM_CB_START int rc=0;    if ( mjm_debug_obs){printf("%s\n",nm);  mjm_note(nc);}

#define MJM_CB_END return  rc;


extern MJM_CB_SIGS(mjm_hist_observer) ;
extern MJM_CB_SIGS(mjm_log_observer);
extern MJM_CB_SIGS(mjm_menu_config_observer);
extern MJM_CB_SIGS(mjm_reply_observer);
extern MJM_CB_SIGS(mjm_abort_key_config_observer);
extern MJM_CB_SIGS(mjm_dlg_index_observer);


extern int mjm_send_via_mutt
//(char ** lto, char * subject, char * text, char ** attach, struct Email * e);
(const char ** lto, const char * subject, const char * text, const char ** attach, struct Email * e);

extern int mutt_main_server(const char * config, char * envp[]);
//extern char * message_text(struct Email * e);
extern void mutt_mjm_free(void * p); //  { free(p); }

extern char * mutt_mjm_message_text(struct Email * e); // 



#ifdef __cplusplus
}
#endif
/*
enum NotifyType
{
  NT_ACCOUNT, ///< Account has changed
  NT_COLOR,   ///< Colour has changed
  NT_COMMAND, ///< A Command has been executed
  NT_CONFIG,  ///< Config has changed
  NT_CONTEXT, ///< Context has changed
  NT_EMAIL,   ///< Email has changed
  NT_GLOBAL,  ///< Not object-related
  NT_MAILBOX, ///< Mailbox has changed
};


*/
#endif // MJM_SERVER_H__



