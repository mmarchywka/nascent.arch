/**
 * @file
 * API for mailboxes
 *
 * @authors
 * Copyright (C) 1996-2002,2013 Michael R. Elkins <me@mutt.org>
 * Copyright (C) 1999-2002 Thomas Roessler <roessler@does-not-exist.org>
 * Copyright (C) 2017-2018 Richard Russon <rich@flatcap.org>
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

#ifndef MUTT_MX_H
#define MUTT_MX_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "config/lib.h"
#include "core/lib.h"

struct Email;
struct Context;
struct stat;

extern const struct MxOps *mx_ops[];

/* These Config Variables are only used in mx.c */
extern unsigned char C_CatchupNewsgroup;
extern bool          C_KeepFlagged;
extern unsigned char C_MboxType;
extern unsigned char C_Move;
extern char *        C_Trash;

extern struct EnumDef MboxTypeDef;

/* flags for mutt_open_mailbox() */
typedef uint8_t OpenMailboxFlags;   ///< Flags for mutt_open_mailbox(), e.g. #MUTT_NOSORT
#define MUTT_OPEN_NO_FLAGS       0  ///< No flags are set
#define MUTT_NOSORT        (1 << 0) ///< Do not sort the mailbox after opening it
#define MUTT_APPEND        (1 << 1) ///< Open mailbox for appending messages
#define MUTT_READONLY      (1 << 2) ///< Open in read-only mode
#define MUTT_QUIET         (1 << 3) ///< Do not print any messages
#define MUTT_NEWFOLDER     (1 << 4) ///< Create a new folder - same as #MUTT_APPEND,
                                    ///< but uses mutt_file_fopen() with mode "w" for mbox-style folders.
                                    ///< This will truncate an existing file.
#define MUTT_PEEK          (1 << 5) ///< Revert atime back after taking a look (if applicable)
#define MUTT_APPENDNEW     (1 << 6) ///< Set in mx_open_mailbox_append if the mailbox doesn't exist.
                                    ///< Used by maildir/mh to create the mailbox.

typedef uint8_t MsgOpenFlags;      ///< Flags for mx_msg_open_new(), e.g. #MUTT_ADD_FROM
#define MUTT_MSG_NO_FLAGS       0  ///< No flags are set
#define MUTT_ADD_FROM     (1 << 0) ///< add a From_ line
#define MUTT_SET_DRAFT    (1 << 1) ///< set the message draft flag

/**
 * enum MxCheckReturns - Return values from mx_mbox_check()
 */
enum MxCheckReturns
{
  MUTT_NEW_MAIL = 1, ///< New mail received in Mailbox
  MUTT_LOCKED,       ///< Couldn't lock the Mailbox
  MUTT_REOPENED,     ///< Mailbox was reopened
  MUTT_FLAGS,        ///< Nondestructive flags change (IMAP)
};

/**
 * struct Message - A local copy of an email
 */
struct Message
{
  FILE *fp;             ///< pointer to the message data
  char *path;           ///< path to temp file
  char *committed_path; ///< the final path generated by mx_msg_commit()
  bool write;           ///< nonzero if message is open for writing
  struct
  {
    bool read : 1;
    bool flagged : 1;
    bool replied : 1;
    bool draft : 1;
  } flags;
  time_t received; ///< the time at which this message was received
};

/**
 * struct MxOps - The Mailbox API
 *
 * Each backend provides a set of functions through which the Mailbox, messages,
 * tags and paths are manipulated.
 */
struct MxOps
{
  enum MailboxType type;  ///< Mailbox type, e.g. #MUTT_IMAP
  const char *name;       ///< Mailbox name, e.g. "imap"
  bool is_local;          ///< True, if Mailbox type has local files/dirs

  /**
   * ac_find - Find an Account that matches a Mailbox path
   * @param a    Account to search
   * @param path Path to search for
   * @retval  0 Success
   * @retval -1 Error
   */
  struct Account *(*ac_find)  (struct Account *a, const char *path);
  /**
   * ac_add - Add a Mailbox to an Account
   * @param a Account to add to
   * @param m Mailbox to add
   * @retval  0 Success
   * @retval -1 Error
   */
  int             (*ac_add)   (struct Account *a, struct Mailbox *m);
  /**
   * mbox_open - Open a Mailbox
   * @param m Mailbox to open
   * @retval  0 Success
   * @retval -1 Error
   * @retval -2 Aborted
   */
  int (*mbox_open)       (struct Mailbox *m);
  /**
   * mbox_open_append - Open a Mailbox for appending
   * @param m     Mailbox to open
   * @param flags Flags, see #OpenMailboxFlags
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*mbox_open_append)(struct Mailbox *m, OpenMailboxFlags flags);
  /**
   * mbox_check - Check for new mail
   * @param m          Mailbox
   * @param index_hint Remember our place in the index
   * @retval >0 Success, e.g. #MUTT_REOPENED
   * @retval -1 Error
   */
  int (*mbox_check)      (struct Mailbox *m, int *index_hint);
  /**
   * mbox_check_stats - Check the Mailbox statistics
   * @param m     Mailbox to check
   * @param flags Function flags
   * @retval  0 Success, no new mail
   * @retval >0 Success, number of new emails
   * @retval -1 Failure
   */
  int (*mbox_check_stats)(struct Mailbox *m, int flags);
  /**
   * mbox_sync - Save changes to the Mailbox
   * @param m          Mailbox to sync
   * @param index_hint Remember our place in the index
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*mbox_sync)       (struct Mailbox *m, int *index_hint);
  /**
   * mbox_close - Close a Mailbox
   * @param m Mailbox to close
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*mbox_close)      (struct Mailbox *m);
  /**
   * msg_open - Open an email message in a Mailbox
   * @param m     Mailbox
   * @param msg   Message to open
   * @param msgno Index of message to open
   * @retval  0 Success
   * @retval -1 Error
   */
  int (*msg_open)        (struct Mailbox *m, struct Message *msg, int msgno);
  /**
   * msg_open_new - Open a new message in a Mailbox
   * @param m   Mailbox
   * @param msg Message to open
   * @param e   Email
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*msg_open_new)    (struct Mailbox *m, struct Message *msg, struct Email *e);
  /**
   * msg_commit - Save changes to an email
   * @param m   Mailbox
   * @param msg Message to commit
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*msg_commit)      (struct Mailbox *m, struct Message *msg);
  /**
   * msg_close - Close an email
   * @param m   Mailbox
   * @param msg Message to close
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*msg_close)       (struct Mailbox *m, struct Message *msg);
  /**
   * msg_padding_size - Bytes of padding between messages
   * @param m Mailbox
   * @retval num Bytes of padding
   */
  int (*msg_padding_size)(struct Mailbox *m);
  /**
   * msg_save_hcache - Save message to the header cache
   * @param m Mailbox
   * @param e Email
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*msg_save_hcache) (struct Mailbox *m, struct Email *e);
  /**
   * tags_edit - Prompt and validate new messages tags
   * @param m      Mailbox
   * @param tags   Existing tags
   * @param buf    Buffer to store the tags
   * @param buflen Length of buffer
   * @retval -1 Error
   * @retval  0 No valid user input
   * @retval  1 Buf set
   */
  int (*tags_edit)       (struct Mailbox *m, const char *tags, char *buf, size_t buflen);
  /**
   * tags_commit - Save the tags to a message
   * @param m Mailbox
   * @param e Email
   * @param buf Buffer containing tags
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*tags_commit)     (struct Mailbox *m, struct Email *e, char *buf);
  /**
   * path_probe - Does this Mailbox type recognise this path?
   * @param path Path to examine
   * @param st   stat buffer (for local filesystems)
   * @retval num Type, e.g. #MUTT_IMAP
   */
  enum MailboxType (*path_probe)(const char *path, const struct stat *st);
  /**
   * path_canon - Canonicalise a Mailbox path
   * @param buf    Path to modify
   * @param buflen Length of buffer
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*path_canon)      (char *buf, size_t buflen);
  /**
   * path_pretty - Abbreviate a Mailbox path
   * @param buf    Path to modify
   * @param buflen Length of buffer
   * @param folder Base path for '=' substitution
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*path_pretty)     (char *buf, size_t buflen, const char *folder);
  /**
   * path_parent - Find the parent of a Mailbox path
   * @param buf    Path to modify
   * @param buflen Length of buffer
   * @retval  0 Success
   * @retval -1 Failure
   */
  int (*path_parent)     (char *buf, size_t buflen);
};

/* Wrappers for the Mailbox API, see MxOps */
int             mx_mbox_check      (struct Mailbox *m, int *index_hint);
int             mx_mbox_check_stats(struct Mailbox *m, int flags);
int             mx_mbox_close      (struct Context **ptr);
struct Context *mx_mbox_open       (struct Mailbox *m, OpenMailboxFlags flags);
int             mx_mbox_sync       (struct Mailbox *m, int *index_hint);
int             mx_msg_close       (struct Mailbox *m, struct Message **msg);
int             mx_msg_commit      (struct Mailbox *m, struct Message *msg);
struct Message *mx_msg_open_new    (struct Mailbox *m, struct Email *e, MsgOpenFlags flags);
struct Message *mx_msg_open        (struct Mailbox *m, int msgno);
int             mx_msg_padding_size(struct Mailbox *m);
int             mx_save_hcache     (struct Mailbox *m, struct Email *e);
int             mx_path_canon      (char *buf, size_t buflen, const char *folder, enum MailboxType *type);
int             mx_path_canon2     (struct Mailbox *m, const char *folder);
int             mx_path_parent     (char *buf, size_t buflen);
int             mx_path_pretty     (char *buf, size_t buflen, const char *folder);
enum MailboxType mx_path_probe     (const char *path);
struct Mailbox *mx_path_resolve    (const char *path);
int             mx_tags_commit     (struct Mailbox *m, struct Email *e, char *tags);
int             mx_tags_edit       (struct Mailbox *m, const char *tags, char *buf, size_t buflen);

struct Account *mx_ac_find   (struct Mailbox *m);
struct Mailbox *mx_mbox_find (struct Account *a, const char *path);
struct Mailbox *mx_mbox_find2(const char *path);
int             mx_ac_add    (struct Account *a, struct Mailbox *m);
int             mx_ac_remove (struct Mailbox *m);

int                 mx_access           (const char *path, int flags);
void                mx_alloc_memory     (struct Mailbox *m);
int                 mx_check_empty      (const char *path);
void                mx_fastclose_mailbox(struct Mailbox *m);
const struct MxOps *mx_get_ops          (enum MailboxType type);
bool                mx_tags_is_supported(struct Mailbox *m);

#endif /* MUTT_MX_H */