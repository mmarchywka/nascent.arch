#ifndef MJM_MAIN_STRUCTS_HH__
#define MJM_MAIN_STRUCTS_HH__

/**
 * @file
 * Command line processing
 *
 * @authors
 * Copyright (C) 1996-2007,2010,2013 Michael R. Elkins <me@mutt.org>
 * Copyright (C) 1999-2007 Thomas Roessler <roessler@does-not-exist.org>
 * Copyright (C) 2004 g10 Code GmbH
 * Copyright (C) 2019 Pietro Cerutti <gahr@gahr.ch>
 *
 * @copyright
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @page main Command line processing
 *
 * Command line processing
 */
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <locale.h>
#include <pwd.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stdbool.h>
#include <stdio.h>
#include "mutt/observer.h"
#include "mutt/notify_type.h"
#include "mutt/string2.h"
#include "config/subset.h"
#include "mutt_commands.h"

#include "globals.h"
#include "config.h"
#include <stdint.h>
#include "mutt/lib.h"
#include "mutt.h"
#include "context.h"
#include "email/email_globals.h"
#include "email/url.h"
#include "core/mailbox.h"
#include "email/envelope.h"
#include "email/email.h"
#include "email/lib.h"
#include "send.h"

#include "mutt/lib.h"
#include "address/lib.h"
#include "config/lib.h"
#include "email/lib.h"
#include "core/lib.h"
#include "conn/lib.h"
#include "gui/lib.h"
#include "debug/lib.h"
#include "alias.h"
#include "browser.h"
#include "commands.h"
#include "context.h"
#include "globals.h"
#include "hook.h"
#include "index.h"
#include "init.h"
#include "keymap.h"
#include "mutt_attach.h"
#include "mutt_history.h"
#include "mutt_logging.h"
#include "mutt_mailbox.h"
#include "mutt_menu.h"
#include "muttlib.h"
#include "mx.h"
#include "myvar.h"
#include "options.h"
#include "protos.h"
#include "send.h"
#include "sendlib.h"
#include "version.h"
#include "ncrypt/lib.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#endif


struct MjmMsgLists  {
// 
  struct ListHead attach; //  = STAILQ_HEAD_INITIALIZER(attach);
  struct ListHead commands; //  = STAILQ_HEAD_INITIALIZER(commands);
  struct ListHead queries; //  = STAILQ_HEAD_INITIALIZER(queries);
struct ListHead alias_queries; //  = STAILQ_HEAD_INITIALIZER(alias_queries);
 // struct ListHead alias_queries; //  = STAILQ_HEAD_INITIALIZER(alias_queries);
  struct ListHead cc_list; // = STAILQ_HEAD_INITIALIZER(cc_list);
  struct ListHead bcc_list;// = STAILQ_HEAD_INITIALIZER(bcc_list);

} ;


struct main_flags {
  bool explicit_folder; //  = false;
  bool dump_variables; //  = false;
  bool hide_sensitive; //  = false;
  bool batch_mode; //  = false;
  bool edit_infile; //  = false;
#ifdef USE_DEBUG_PARSE_TEST
  bool test_config;//  = false;
#endif
};


struct main_send_flags {

  char *subject; //  = NULL;
  char *include_file; //  = NULL;
  char *draft_file; //  = NULL;
  char *new_type;//  = NULL;
  char *dlevel; //  = NULL;
  char *dfile; //  = NULL;
#ifdef USE_NNTP
  char *cli_nntp; //  = NULL;
#endif


}; 

struct cmd_line_stuff
{ 
 int double_dash; //  = argc, nargc = 1;
 int  nargc; //  = 1;
  int rc; //  = 1;
  bool repeat_error; //  = false;
  struct Buffer folder; //  = mutt_buffer_make(0);
  struct Buffer expanded_infile; //  = mutt_buffer_make(0);
  struct Buffer tempfile; //  = mutt_buffer_make(0);
  struct ConfigSet *cs; //  = NULL;

}; //cmd_line_stuff


struct compose_files {   
 FILE *fp_in; //  = NULL;
    FILE *fp_out; //  = NULL;
    char *infile; //  = NULL;
    char *bodytext; //  = NULL;
    const char *bodyfile; //  = NULL;
    int rv; //  = 0;
};



void fkss(struct MjmMsgLists * pmml, struct ListHead*  pfuk);
void afck(struct MjmMsgLists*  pmml);
void main_flags_init(struct main_flags * pmf);
void cmd_line_stuff_init(struct cmd_line_stuff * clsp, int argc);


void zero_my_hero(void * p, int n); //  { memset(p,0,n); } ;

int formalities(int argc, char *argv[], char *envp[]);
void get_cc_lists(struct Email ** ep, struct MjmMsgLists mml);
int draft_proc(struct Email * e,struct compose_files cf,struct main_send_flags  msf,struct cmd_line_stuff cls
,            bool C_ResumeEditedDraftFiles);


//struct MjmMsgLists mml;
int unf_attach_stuff(struct Email * e, const char** nm) ;

int attach_stuff(struct Email * e,// struct MjmMsgLists mml,
struct MjmMsgLists * pmml) ;

int  prepapre_fp_in(struct Email * e,struct compose_files * pcf,
struct main_flags * pmf,struct cmd_line_stuff * pcls);

int add_mailto(struct Email * e,char ** bodytext,char * arg);
#endif

