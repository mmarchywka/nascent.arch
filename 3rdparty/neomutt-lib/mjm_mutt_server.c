
#include "mutt/observer.h" 
#include "mutt/notify_type.h" 
#include "mutt/string2.h" 
#include "config/subset.h" 

#include "globals.h"
#include "config.h"
#include <stdint.h>
#include "mutt/lib.h"
#include "mutt.h"
#include "context.h"
#include "email/email_globals.h"
#include "core/mailbox.h"
#include "email/envelope.h"
#include "email/email.h"
#include "email/lib.h"
#include "send.h"


#include "mjm_mutt_server.h"
#include "mjm_readline.h"
#include "mjm_main_structs.h"
#include "mutt_commands.h"
//#include "debug/notify.c"
#include <string.h>

 MJM_CB_SIGS(mjm_h1);

 bool mjm_daemon_mode=false;
 bool mjm_init_only=false;
 bool mjm_debug_obs=false;
extern int  mutt_set_init_only(bool x) //   cls.rc = 0;
{ mjm_init_only=x; return 0; } 
int mjm_note(struct NotifyCallback*nc )
{
// notify_dump_command(nc);


 struct Command *cmd = nc->event_data;
if (cmd!=0) {
  if (cmd->data < 4096)
   printf("\tCommand: %s, data: %ld\n", cmd->name, cmd->data);
  else
    printf("\tCommand: %s, data: %p\n", cmd->name, (void *) cmd->data);
}


int et=nc->event_type;
  if (et==NT_ACCOUNT) printf( " Account has changed");
else if(et==  NT_COLOR) printf("   ///< Colour has changed");
 else if(et==   NT_COMMAND)printf(" ///< A Command has been executed");
 else if(et==   NT_CONFIG)printf("  ///< Config has changed");
 else if(et==   NT_CONTEXT)printf(" ///< Context has changed");
 else if(et==   NT_EMAIL) printf("   ///< Email has changed");
 else if(et==   NT_GLOBAL) printf("  ///< Not object-related");
 else if(et==   NT_MAILBOX) printf("  ///< Mailbox has changed");
 else  printf("  bbad event %d ",et );
printf("\n");

return 0; 

}

//struct Context * Context;

#define BOOLPRINT(x) printf("... %s  %s  ",#x, ((x)?"set":"clear") );

 MJM_CB_SIG(mjm_h1)
/*
printf("asdcasdc\n");
 if (!nc->event_data)
    return -1;
  if (nc->event_type != NT_CONFIG)
    return 0;
printf("xxxxxxxxxa\n");

  struct EventConfig *ec = nc->event_data;

  if (mutt_str_strcmp(ec->name, "reply_regex") != 0)
    return 0;
printf("eeeeeeeee\n");
*/
printf("context %d \n", Context );
if ( Context ) printf("mbox %d \n", Context->mailbox  );
  if (!Context || !Context->mailbox)
    return 0;

  regmatch_t pmatch[1];
  struct Mailbox *m = Context->mailbox;

printf(" msg count %d \n", m->msg_count );
  for (int i = 0; i < m->msg_count; i++)
  {
    struct Email *e = m->emails[i];
printf(" email  %d\n",e);
    if (!e)
      break;
    struct Envelope *env = e->env;
    if (!env || !env->subject)
      continue;
BOOLPRINT((e->read));
BOOLPRINT((e->old));
printf("   asdfasd  %s %s\n",env->subject, env->message_id);

// .. try to read here.
if ((i==0)||!(e->read))
{
struct Mailbox * m=Context->mailbox;
  mutt_parse_mime_message(m, e);
FILE * fp_out=stdout;
int  res = mutt_copy_message(fp_out, m, e, 0, 0, 80000);


} // read 


    if (mutt_regex_capture(C_ReplyRegex, env->subject, 1, pmatch))
    {
      env->real_subj = env->subject + pmatch[0].rm_eo;
      continue;
    }
// this code never executes.. 
    env->real_subj = env->subject;
printf("   %s %s\n",env->subject, env->message_id);
  }

//  OptResortInit = true; /* trigger a redraw of the index */
  return 0;
}

struct Email ** mutt_email_list(int flags)
{
bool uro=((flags&1)!=0);
//printf("context %d \n", Context );
//if ( Context ) printf("mbox %d \n", Context->mailbox  );
  if (!Context || !Context->mailbox)
{
printf(" Context or mailbox is null in mutt_email_list\n");
    return 0;
}
  struct Mailbox *m = Context->mailbox;
printf(" mutt_email_list  msg count %d \n", m->msg_count );
struct Email ** elist=malloc((m->msg_count+1)*sizeof(struct Email*));
int j=0;
  for (int i = 0; i < m->msg_count; i++)
  {
    struct Email *e = m->emails[i];
//printf(" email  %d\n",e);
    if (!e) break;
    struct Envelope *env = e->env;
    if (!env ) continue;
if (uro) if (e->read) continue;
if (e->deleted) continue;
elist[j]=e; ++j;
//BOOLPRINT((e->read)); BOOLPRINT((e->old));
// .. try to read here.

}
elist[j]=0;
return elist;
} // mutt_email_list


void mutt_mjm_free(void * p) { free(p); } 

char * mutt_mjm_message_text(struct Email * e)
{
struct Mailbox * m=Context->mailbox;
  mutt_parse_mime_message(m, e);
const char* fn=".mjm_temp_xxx.txt";
FILE * fp_out=fopen(fn,"w"); // stdout;
int  res = mutt_copy_message(fp_out, m, e, 0, 0, 80000);
fclose(fp_out);
FILE * f = fopen (fn, "rb");
char * buffer=0;
if (f)
{
  fseek (f, 0, SEEK_END);
int   length = ftell (f);
  fseek (f, 0, SEEK_SET);
   buffer = malloc (length+1);
  if (buffer) { fread (buffer, 1, length, f); }
buffer[length]=0;
  fclose (f);
}
return buffer;
}

int mjm_mutt_move_msg(struct Email * e, const char * dest, const int flags )
{
int myflags=1; // headless
struct EmailList el = STAILQ_HEAD_INITIALIZER(el);
emaillist_add_email(&el, e);
bool del_orig=((flags&1)!=0);
bool create_missing=((flags&2)!=0);
if (create_missing) myflags|=2; // create_missing 
int rc1= imap_copy_messages_mjm(Context->mailbox, &el, dest, del_orig,myflags);
if (rc1)  return rc1;
int rc=imap_sync_mailbox(Context->mailbox, false,false);
return rc;
/* message.c */

//int imap_copy_messages(struct Mailbox *m, struct EmailList *el, const char *dest, bool delete_original);

//  struct EmailList el = STAILQ_HEAD_INITIALIZER(el);
 // emaillist_add_email(&el, e);

}

//int set_read_flag(struct Email * e) { mutt_set_flag(m, e, MUTT_READ, true); }
int set_read_flag(struct Email * e) 
{ mutt_set_flag(Context->mailbox, e, MUTT_READ, true);
//Context->mailbox->append=true;
int rc=imap_sync_mailbox(Context->mailbox, false,false);
//int rc=mx_msg_commit(Context->mailbox, e);
printf(" seetinc commit %d\n",rc);
return rc; 
}
int reset_read_flag(struct Email * e) 
{ mutt_set_flag(Context->mailbox, e, MUTT_READ, !true);
//Context->mailbox->append=true;
int rc=imap_sync_mailbox(Context->mailbox, false,false);
//int rc=mx_msg_commit(Context->mailbox, e);
printf(" resseetinc commit %d\n",rc);
return rc;
 }

//int mx_msg_commit(struct Mailbox *m, struct Message *msg)


 MJM_CB_SIG(mjm_reply_observer)
 MJM_CB_START 
//mjm_h1(nc);

MJM_CB_END } // reply_observer


MJM_CB_SIG(mjm_hist_observer) 
 MJM_CB_START 

//mjm_h1(nc);
MJM_CB_END }
 MJM_CB_SIG(mjm_log_observer)
 MJM_CB_START 
//mjm_h1(nc);

MJM_CB_END }
 MJM_CB_SIG(mjm_menu_config_observer)
 
//mjm_h1(nc);

MJM_CB_START MJM_CB_END }



 MJM_CB_SIG(mjm_abort_key_config_observer)  MJM_CB_START

//mjm_h1(nc);

 MJM_CB_END }
 MJM_CB_SIG(mjm_dlg_index_observer)  MJM_CB_START
//mjm_h1(nc);

 MJM_CB_END }

////////////////////////////////////////////////////////////////

int old_mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach,
struct Email * e)
{
int rc=0;
char * to=(lto==0)?0:*lto;
int sendflags=0;
//if (e==0)       struct Email * e = email_new();
if (e==0)     {   e = email_new();} else {  sendflags=SEND_REPLY; } 
//  e = email_new(); // } else {  sendflags=SEND_REPLY; } 
char * body=0;
    if (e->env==0)   e->env = mutt_env_new();
while (true)
{
      if (0!=add_mailto(e,&body,to))
      {
        //mutt_error(_("Failed to parse mailto: link"));
        printf("Failed to parse  %s\n",to);
	rc=-1;
	break; // "goto" lol 
      }
      e->env->subject = mutt_str_strdup(subject);
 // if (0!= prepapre_fp_in(e,&cf,&mf,&cls))
//
const char * tempf="tempfile";
 struct MjmMsgLists mml;
afck(&mml);
FILE * fout= fopen(tempf, "w");
         if (body!=0) fputs(body, fout);
      if (body!=0) printf("mjm_send_via_mutt body %s  \n",body);
      if (text!=0) printf("mjm_send_via_mutt text %s  \n",text);
         fputs(text, fout);
//fclose(fout);

/*
if (attach!=0)
{
int  i=0; while (attach[i]!=0)
{ 
printf(" attaching %s\n",attach[i]);
mutt_list_insert_tail(&mml.attach, mutt_str_strdup(attach[i])); ++i; } // i 
} // attach==0

if (0!=attach_stuff(e,&mml))
{ printf(" forgot code here in mml check \n"); }
*/
unf_attach_stuff(e,attach);
//mutt_write_mime_body(e->content,fout);
fclose(fout);
if (attach!=0)
{ // mjm 2020-11-28 added SEND_TO_SENDER to avoid prompt in reply
// to list or spam email. 
  //sendflags=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER;
  sendflags=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER|SEND_TO_SENDER;
int  rc = mutt_send_message(sendflags, e,tempf, 0,0);
printf(" attached rc=%d\n",rc);
return rc;
}

//e->content=mutt_make_multipart(e->content);
//if (attach!=0) e->content->next=0;
printf(" done attaching check \n");
   //rc = mutt_send_message(sendflags, e, tempf, NULL, NULL);
   //rc = mutt_send_message(0, e, tempf, NULL, NULL);
    if ( TAILQ_EMPTY(&e->env->to) &&
        TAILQ_EMPTY(&e->env->cc))
    {
      printf("No recipients specified\n");
    }

      printf(" 1trying ... \n");

struct Buffer fck= mutt_buffer_make(0);
mutt_buffer_strcpy(&fck, tempf);
          mutt_buffer_expand_path(&fck);
          //(*pcf).fp_in = fopen(mutt_b2s(fck), "r");
   char *        carp=mutt_b2s(&fck);
printf(" file expa %s\n",carp);



 //rc = mutt_send_message(SEND_BATCH, e, tempf, Context, NULL);
// rc = mutt_send_message(SEND_BATCH, e, carp, Context, NULL);

  struct EmailList el = STAILQ_HEAD_INITIALIZER(el);
  //emaillist_add_email(&el, e_cur);
bool kluge=true;
//if (!kluge) 
 emaillist_add_email(&el, e);
/////////////////////////////////////////
//print(" attach to el doh\n");
//if (0!=attach_stuff(el,&mml))
//{ printf(" forgot code here in mml check \n"); }


///////////////////////////////////////

  //int rc = mutt_send_message(SEND_RESEND, e_new, NULL, ctx, &el);
  //rc = mutt_send_message(SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER, e, tempf, Context, &el);
  sendflags|=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER;
  if (sendflags==0) rc = mutt_send_message(sendflags, e, tempf, Context, &el);
// this did send a reply but also 2 dups 
// see if this does atach 
  else
{
if (!kluge)  rc = mutt_send_message(sendflags, 0,tempf, Context, &el);
 //else rc = mutt_send_message(sendflags, e,tempf, Context, &el);
 else
{

  //sendflags=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER;
 rc = mutt_send_message(sendflags, e,tempf, 0,0);
}

}
  emaillist_clear(&el);

#if 0
 mutt_buffer_mktemp(&cls.tempfile);

        cf.fp_out = mutt_file_fopen(mutt_b2s(&cls.tempfile), "w");
        if (!cf.fp_out)
        {
          mutt_file_fclose(&cf.fp_in);
          mutt_perror(mutt_b2s(&cls.tempfile));
          email_free(&e);
          goto main_curses; // TEST29: neomutt -H existing-file (where tmpdir=/path/to/FILE blocking tmpdir)
        }
        if (cf.fp_in)
        {
          mutt_file_copy_stream(cf.fp_in, cf.fp_out);
          if (cf.fp_in != stdin)
            mutt_file_fclose(&cf.fp_in);
        }
        else if (cf.bodytext)
          fputs(cf.bodytext, cf.fp_out);
        mutt_file_fclose(&cf.fp_out);

        cf.fp_in = fopen(mutt_b2s(&cls.tempfile), "r");
    FREE(&cf.bodytext);

if (0!=attach_stuff(e,&mml))

#endif



//    rc = mutt_send_message(sendflags, e, cf.bodyfile, NULL, NULL);

break;
} // loop 
printf(" meed to delete the entire elist when done doh\n");
//email_free(&e);
return rc; 
} // mjm_send_via_mutt

////////////////////////////////////////////////////////////////////////////

//int mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach,
int mjm_send_via_mutt(const char ** lto, const char * subject, const char * text, const char ** attach,
struct Email * e)
{
int rc=0;
const char * to=(lto==0)?0:*lto;
int sendflags=0;
bool wemade=false;
if (e==0)     { wemade=true;   e = email_new();} else {  sendflags=SEND_REPLY|SEND_TO_SENDER; } 
char * body=0;
if (e->env==0)   e->env = mutt_env_new();
// never loops, scoping and "finally" block 
while (true)
{
      if (0!=add_mailto(e,&body,to))
      {
        printf("Failed to parse  %s\n",to);
	rc=-1;
	break; // "goto" lol 
      }
      e->env->subject = mutt_str_strdup(subject);
//
const char * tempf="tempfile";
printf(" dangerous tempfile name need to fix \n");
// struct MjmMsgLists mml; afck(&mml);
FILE * fout= fopen(tempf, "w");
         if (body!=0) fputs(body, fout);
      if (body!=0) printf("mjm_send_via_mutt body %s  \n",body);
      if (text!=0) printf("mjm_send_via_mutt text %s  \n",text);
         fputs(text, fout);

unf_attach_stuff(e,attach);
fclose(fout);
if (attach!=0)
{
// see if this works in reply mode
  sendflags=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER|SEND_TO_SENDER;
//  sendflags|=SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER;
 rc = mutt_send_message(sendflags, e,tempf, 0,0);
printf(" attached rc=%d\n",rc);
break; // return rc;
}
printf(" done attaching check \n");
    if ( TAILQ_EMPTY(&e->env->to) &&
        TAILQ_EMPTY(&e->env->cc))
    { printf("No recipients specified\n"); }

      printf("2trying ... \n");

  struct EmailList el = STAILQ_HEAD_INITIALIZER(el);
 emaillist_add_email(&el, e);

  sendflags|=SEND_TO_SENDER|SEND_BATCH|SEND_MJM_ASIS|SEND_NO_FREE_HEADER;
  if (sendflags==0) rc = mutt_send_message(sendflags, e, tempf, Context, &el);
  else
{
rc = mutt_send_message(sendflags, 0,tempf, Context, &el);

}
  emaillist_clear(&el);

// never loops, just wanted a finally block but now no need 
break;
} // fake loop 
printf(" meed to delete the entire elist when done doh\n");
if (wemade) email_free(&e);
return rc; 
} // mjm_send_via_mutt



//////////////////////////////////////////////////////////////////////
void reply_to_a_msg(int n)
{
printf(" tried to find %d \n",n);
struct Email ** el= mutt_email_list(0);
if (el==0) return;
//int mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach, struct Email * e)
int i=0;
while ( el[i]!=0)
{
if (i==n)
{
//char ** lto = 0; // new char*[4];
char ** lto =malloc(sizeof( char*)*4);
lto[0]="marchywka@hotmail.com";
lto[1]=0;
char * subject="test reply";
char * text=" added to madg";

//int rc= mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach, struct Email * e)
struct Envelope *env = el[i]->env;
if (!env || !env->subject) continue;
BOOLPRINT((el[i]->read));
BOOLPRINT((el[i]->old));
printf(" try   %s %s\n",env->subject, env->message_id);
printf(" tri to send %d %d\n",i,n);
int rc= mjm_send_via_mutt( lto,  subject,  text, 0, el[i]);

printf(" tried to send %d %d %d\n",rc,i,n);
}
else email_free(&el[i]);
++i;
}


free(el);
} // reply_to_a_msg


void reply_to_a_msg_idx(int n)
{
printf(" tried to find %d \n",n);
struct Email ** el= mutt_email_list(0);
if (el==0) return;
//int mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach, struct Email * e)
int i=0;
while ( el[i]!=0)
{
if ((el[i]->index)==n)
{
//char ** lto = 0; // new char*[4];
const char ** lto =malloc(sizeof( char*)*4);
lto[0]="marchywka@hotmail.com";
lto[1]=0;
const char * subject="test reply";
const char * text=" added to madg";

//int rc= mjm_send_via_mutt(char ** lto, char * subject, char * text, char ** attach, struct Email * e)
struct Envelope *env = el[i]->env;
if (!env || !env->subject) continue;
BOOLPRINT((el[i]->read));
BOOLPRINT((el[i]->old));
printf(" try   %s %s\n",env->subject, env->message_id);
printf(" tri to send %d %d\n",i,n);
int rc= mjm_send_via_mutt( lto,  subject,  text, 0, el[i]);

printf(" tried to send %d %d %d\n",rc,i,n);
}
else email_free(&el[i]);
++i;
}


free(el);
} // reply_to_a_msg


int  cmd_line_send(const char * _s) 
{
const int sz=strlen(_s)+10;
char s[sz];
memcpy(s,_s,sz);
//FILE * rinput=rl_instream;
//rl_instream=~0; // disable for now 
const char * tok=strtok(s," ");
//while (tok!=0) { tok=strtok(0," "); } // while 
const char * to=strtok(0," ");
const char * subj=strtok(0," ");
const char * text=s; // +2+strlen(to)+strlen(subj);
printf(" to %s\n", to);
printf(" subh  %s\n", subj);
printf(" text %s\n", text);
if (text==0) return -1;
int rc= mjm_send_via_mutt(&to,  subj,  text,0,0);
printf(" rc = %d\n",rc);
//rl_instream=rinput; // enable for now 
return 0;
}


int  cmd_line_send_file(const char * _s) 
{
const int sz=strlen(_s)+10;
char s[sz];
memcpy(s,_s,sz);
const char * tok=strtok(s," ");
const char * to=strtok(0," ");
const char * subj=strtok(0," ");
const char * fn=strtok(0," ");
//const char * text=s; // +2+strlen(to)+strlen(subj);
// https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
char * buffer = 0;
long length;
printf(" try to open %s\n",fn);
FILE * f = fopen (fn, "rb");
if (f)
{
  fseek (f, 0, SEEK_END);
  length = ftell (f);
  fseek (f, 0, SEEK_SET);
  buffer = malloc (length);
  if (buffer) { fread (buffer, 1, length, f); }
  fclose (f);
}
else return -1;
printf(" to %s\n", to);
printf(" subh  %s\n", subj);
printf(" text %s\n", buffer);
if (buffer==0) return -1;
int rc= mjm_send_via_mutt(&to,  subj,  buffer,0,0);
printf(" rc = %d\n",rc);
free(buffer);
return 0;
}

int  cmd_check(const char * s) 
{
bool f=false;
if (s[5]!=0) if (s[6]=='1') f=true;
printf(" checking with flag %s \n",f?"T":"F");
int c1= imap_check_mailbox(Context->mailbox, f);
printf(" checking with flag %s gives %d  \n",f?"T":"F",c1);
return 0;
}


int  mjm_input_loop()
{
//printf(" just test without input\n");
//{ reply_to_a_msg(3); }
while (true)
{
rl_gets();
if (line_read==0) break;
if (*line_read==0) break;
if (strcmp(line_read,"dump")==0 ) { mjm_h1(0); } 
else if (strcmp(line_read,"verbose")==0 ) { mjm_debug_obs=true; } 
//int mjm_send_via_mutt(char * to, char * subject, char * text)
else if (strcmp(line_read,"quiet")==0 ) { mjm_debug_obs=!true; } 
else if (strncmp(line_read,"sendf",5)==0 ) { cmd_line_send_file(line_read); } 
else if (strncmp(line_read,"send",4)==0 ) { cmd_line_send(line_read); } 
else if (strncmp(line_read,"check",5)==0 ) { cmd_check(line_read); } 
else if (strncmp(line_read,"reply",5)==0 ) 
{ reply_to_a_msg(atoi(line_read+6)); }
else printf(" bad command %s\n",line_read);
 
}// true
return 0; 
} // mjm_input_loop


