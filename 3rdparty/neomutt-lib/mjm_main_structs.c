
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
#include "mjm_main_structs.h"
extern void init_locale();

void fkss(struct MjmMsgLists * pmml, struct ListHead*  pfuk)
{
  struct ListHead * sht =pfuk; 
struct ListHead afck = STAILQ_HEAD_INITIALIZER((*sht));
*pfuk=afck;

}
void afck(struct MjmMsgLists*  pmml)
{

//  struct ListHead * sht =&(*pmml).alias_queries; 
//struct ListHead afck = STAILQ_HEAD_INITIALIZER((*sht));
//(*pmml).alias_queries=afck;
fkss(pmml, &((*pmml).attach));
fkss(pmml, &((*pmml).commands));
fkss(pmml, &((*pmml).queries));
fkss(pmml, &((*pmml).alias_queries));
fkss(pmml, &((*pmml).cc_list));
fkss(pmml, &((*pmml).bcc_list));
}

void main_flags_init(struct main_flags * pmf)
{
  (*pmf).explicit_folder = false;
  (*pmf).dump_variables = false;
  (*pmf).hide_sensitive = false;
  (*pmf).batch_mode = false;
  (*pmf).edit_infile = false;
#ifdef USE_DEBUG_PARSE_TEST
  (*pmf).test_config = false;
#endif
}


void cmd_line_stuff_init(struct cmd_line_stuff * clsp, int argc)
{
 (*clsp).double_dash = argc;
 (*clsp).nargc  = 1;
 (*clsp).rc = 1;
 (*clsp).repeat_error = false;
 (*clsp).folder = mutt_buffer_make(0);
 (*clsp).expanded_infile = mutt_buffer_make(0);
 (*clsp).tempfile = mutt_buffer_make(0);
 (*clsp).cs  = NULL;


} 


void zero_my_hero(void * p, int n) { memset(p,0,n); } 

int formalities(int argc, char *argv[], char *envp[])
{


  /* sanity check against stupid administrators */
  if (getegid() != getgid())
  {
    mutt_error("%s: I don't want to run with privileges!", argv[0]);
   return ~0; //  goto main_exit; // TEST01: neomutt (as root, chgrp mail neomutt; chmod +s neomutt)
  }

  init_locale();
{
  int out = 0;
  if (mutt_randbuf(&out, sizeof(out)) < 0)
   return ~0; //  goto main_exit; // TEST02: neomutt (as root on non-Linux OS, rename /dev/urandom)
}
  umask(077);

  mutt_envlist_init(envp);
 

return 0; 

}
// copy mml but all pointers... 
void get_cc_lists(struct Email ** ep, struct MjmMsgLists mml)
{
  if (!STAILQ_EMPTY(&mml.cc_list) || !STAILQ_EMPTY(&mml.bcc_list))
  {
    (*ep) = email_new();
    (*ep)->env = mutt_env_new();

    struct ListNode *np = NULL;
    STAILQ_FOREACH(np, &mml.bcc_list, entries)
    {
      mutt_addrlist_parse(&(*ep)->env->bcc, np->data);
    }

    STAILQ_FOREACH(np, &mml.cc_list, entries)
    {
      mutt_addrlist_parse(&(*ep)->env->cc, np->data);
    }

    mutt_list_free(&mml.bcc_list);
    mutt_list_free(&mml.cc_list);
  }
}

int draft_proc(struct Email * e,struct compose_files cf,struct main_send_flags  msf,struct cmd_line_stuff cls
,            bool C_ResumeEditedDraftFiles)
{

        struct Envelope *opts_env = e->env;
        struct stat st;
// mmoved up 
    //    sendflags |= SEND_DRAFT_FILE;

        /* Set up a tmp Email with just enough information so that
         * mutt_prepare_template() can parse the message in fp_in.  */
        struct Email *e_tmp = email_new();
        e_tmp->offset = 0;
        e_tmp->content = mutt_body_new();
        if (fstat(fileno(cf.fp_in), &st) != 0)
        {
          mutt_perror(msf.draft_file);
          //email_free(&e_tmp);
          //email_free(&e);
          //goto main_curses; // TEST31: can't test
		return ~0; 
        }
        e_tmp->content->length = st.st_size;

        if (mutt_prepare_template(cf.fp_in, NULL, e, e_tmp, false) < 0)
        {
          mutt_error(_("Can't parse message template: %s"), msf.draft_file);
      //    email_free(&e);
          email_free(&e_tmp);
     //     goto main_curses;
		return ~0; 
        }

        /* Scan for neomutt header to set C_ResumeDraftFiles */
        struct ListNode *np = NULL, *tmp = NULL;
        STAILQ_FOREACH_SAFE(np, &e->env->userhdrs, entries, tmp)
        {
          if (mutt_str_startswith(np->data, "X-Mutt-Resume-Draft:", CASE_IGNORE))
          {
            if (C_ResumeEditedDraftFiles)
              cs_str_native_set(cls.cs, "resume_draft_files", true, NULL);

            STAILQ_REMOVE(&e->env->userhdrs, np, ListNode, entries);
            FREE(&np->data);
            FREE(&np);
          }
        }

        mutt_addrlist_copy(&e->env->to, &opts_env->to, false);
        mutt_addrlist_copy(&e->env->cc, &opts_env->cc, false);
        mutt_addrlist_copy(&e->env->bcc, &opts_env->bcc, false);
        if (opts_env->subject)
          mutt_str_replace(&e->env->subject, opts_env->subject);

        mutt_env_free(&opts_env);
        email_free(&e_tmp);


return 0;

}


//struct MjmMsgLists mml;
int attach_stuff(struct Email * e,// struct MjmMsgLists mml,
struct MjmMsgLists * pmml) 
{
    if (!STAILQ_EMPTY(&(*pmml).attach))
    {
      struct Body *b = e->content;

      while (b && b->next)
        b = b->next;

      struct ListNode *np = NULL;
      STAILQ_FOREACH(np, &(*pmml).attach, entries)
      {
	printf(" in mjm_main_structs.c  lis attaching %s\n",np->data);   
     if (b)
        {
          b->next = mutt_make_file_attach(np->data);
          b = b->next;
        }
        else
        {
          b = mutt_make_file_attach(np->data);
          e->content = b;
        }
        if (!b)
        {
          printf("%s: unable to attach file", np->data);
          mutt_error(_("%s: unable to attach file"), np->data);
          mutt_list_free(&(*pmml).attach);
//          email_free(&e);
  //        goto main_curses; // TEST32: neomutt john@example.com -a missing
		return ~0;
        }
      }
      mutt_list_free(&(*pmml).attach);
    }
struct Body * b=e->content;
while (b) { printf(" walk body %s\n",b->filename); b=b->next; } 
	return 0;
}

int unf_attach_stuff(struct Email * e,// struct MjmMsgLists mml,
const char** nm) 
{
	if (nm==0) return;
int i=0;
      struct Body *b = e->content;
      while (b && b->next) b = b->next;
while (nm[i]!=0)
{
printf(" unf %d %s \n",i,nm[i]);
     if (b)
        {
          b->next = mutt_make_file_attach(nm[i]);
          b = b->next;
        }
        else
        {
          b = mutt_make_file_attach(nm[i]);
          e->content = b;
        }
        if (!b)
        {
          printf("%s: unable to attach file", nm[i]);
//          mutt_list_free(&(*pmml).attach);
//          email_free(&e);
  //        goto main_curses; // TEST32: neomutt john@example.com -a missing
		return ~0;
        }
i=i+1;
      }
 //     mutt_list_free(&(*pmml).attach);
 b=e->content;
while (b) { printf(" walk body %s\n",b->filename); b=b->next; } 
	return 0;
}




int  prepapre_fp_in(struct Email * e,struct compose_files * pcf,
struct main_flags * pmf,struct cmd_line_stuff * pcls)
	{
      if ((*pcf).infile)
      {
        if (mutt_str_strcmp("-", (*pcf).infile) == 0)
        {
          if ((*pmf).edit_infile)
          {
            mutt_error(_("Can't use -E flag with stdin"));
		return ~0; 
//            email_free(&e);
//            goto main_curses; // TEST27: neomutt -E -H -
          }
          (*pcf).fp_in = stdin;
        }
        else
        {
          mutt_buffer_strcpy(&(*pcls).expanded_infile, (*pcf).infile);
          mutt_buffer_expand_path(&(*pcls).expanded_infile);
          (*pcf).fp_in = fopen(mutt_b2s(&(*pcls).expanded_infile), "r");
          if (!(*pcf).fp_in)
          {
            mutt_perror(mutt_b2s(&(*pcls).expanded_infile));
		return ~0; 
    //        email_free(&e);
    //        goto main_curses; // TEST28: neomutt -E -H missing
          }
        }
      }
return 0;
}

int add_mailto(struct Email * e,char ** bodytext,char * arg)
{
#if 1 
      if (url_check_scheme(arg) == U_MAILTO)
      {
        if (!mutt_parse_mailto(e->env, bodytext, arg))
        {
	
          mutt_error(_("Failed to parse mailto: link"));
	return ~0;
       //   email_free(&e);
       //   goto main_curses; // TEST25: neomutt mailto:?
        }
	else printf(" may have added ok %s\n",arg);
      }
      else
{
        mutt_addrlist_parse(&e->env->to, arg);
	 printf(" tacked on may have added ok %s\n",arg);
}
#endif
return 0;
}

