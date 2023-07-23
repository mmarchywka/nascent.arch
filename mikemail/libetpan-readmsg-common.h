

// copied from the libetpan examples 

#ifndef READMSG_COMMON_H

#define READMSG_COMMON_H

#include <libetpan/libetpan.h>

#define DEST_CHARSET "iso-8859-1"

enum {
/* SEB */
#ifndef NO_ERROR 
  NO_ERROR,
#endif
  ERROR_FILE = 1,
  ERROR_MEMORY = 2,
  ERROR_INVAL = 3,
  ERROR_FETCH = 4
};
/*
int etpan_mime_is_text(struct mailmime * build_info);

int show_part_info(FILE * f,
    struct mailmime_single_fields * mime_fields,
    struct mailmime_content * content);

int etpan_fetch_message(mailmessage * msg_info,
    struct mailmime * mime_part,
    struct mailmime_single_fields * fields,
    char ** result, size_t * result_len);

struct mailimf_fields * fetch_fields(mailmessage * msg_info,
    struct mailmime * mime);

int fields_write(FILE * f, int * col,
    struct mailimf_fields * fields);

*/

#endif

#if 0
#include <sys/stat.h>
#ifndef WIN32
#	include <sys/mman.h>
#	include <unistd.h>
#endif
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#endif

/* returns TRUE is given MIME part is a text part */

int etpan_mime_is_text(struct mailmime * build_info)
{
  if (build_info->mm_type == MAILMIME_SINGLE) {
    if (build_info->mm_content_type != NULL) {
      if (build_info->mm_content_type->ct_type->tp_type ==
          MAILMIME_TYPE_DISCRETE_TYPE) {
        if (build_info->mm_content_type->ct_type->tp_data.tp_discrete_type->dt_type ==
            MAILMIME_DISCRETE_TYPE_TEXT)
          return 1;
      }
    }
    else
      return 1;
  }

  return 0;
}


/* display content type */

int show_part_info(FILE * f,
    struct mailmime_single_fields * mime_fields,
    struct mailmime_content * content)
{
  char * description;
  char * filename;
  int col;
  int r;

  description = mime_fields->fld_description;
  filename = mime_fields->fld_disposition_filename;

  col = 0;

  r = fprintf(f, " [ Part ");
  if (r < 0)
    goto err;

  if (content != NULL) {
    r = mailmime_content_type_write(f, &col, content);
    if (r != MAILIMF_NO_ERROR)
      goto err;
  }

  if (filename != NULL) {
    r = fprintf(f, " (%s)", filename);
    if (r < 0)
      goto err;
  }

  if (description != NULL) {
    r = fprintf(f, " : %s", description);
    if (r < 0)
      goto err;
  }

  r = fprintf(f, " ]\n\n");
  if (r < 0)
    goto err;

  return NO_ERROR;
  
 err:
  return ERROR_FILE;
}

/*
  fetch the data of the mailmime_data structure whether it is a file
  or a string.

  data must be freed with mmap_string_unref()
*/

#if 0
static int fetch_data(struct mailmime_data * data,
    char ** result, size_t * result_len)
{
  int fd;
  int r;
  char * text;
  struct stat buf;
  int res;
  MMAPString * mmapstr;

  switch (data->dt_type) {
  case MAILMIME_DATA_TEXT:
    mmapstr = mmap_string_new_len(data->dt_data.dt_text.dt_data,
        data->dt_data.dt_text.dt_length);
    if (mmapstr == NULL) {
      res = ERROR_MEMORY;
      goto err;
    }

    * result = mmapstr->str;
    * result_len = mmapstr->len;

    return NO_ERROR;

  case MAILMIME_DATA_FILE:
    fd = open(data->dt_data.dt_filename, O_RDONLY);
    if (fd < 0) {
      res = ERROR_FILE;
      goto err;
    }

    r = fstat(fd, &buf);
    if (r < 0) {
      res = ERROR_FILE;
      goto close;
    }

    if (buf.st_size != 0) {
      text = mmap(NULL, buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
      if (text == (char *)MAP_FAILED) {
	res = ERROR_FILE;
	goto close;
      }

      mmapstr = mmap_string_new_len(text, buf.st_size);
      if (mmapstr == NULL) {
        res = r;
        goto unmap;
      }
      
      munmap(text, buf.st_size);
    }
    else {
      mmapstr = mmap_string_new("");
      if (mmapstr == NULL) {
        res = r;
        goto close;
      }
    }

    close(fd);

    * result = mmapstr->str;
    * result_len = mmapstr->len;

    return NO_ERROR;

  default:
    return ERROR_INVAL;
  }
  
 unmap:
  munmap(text, buf.st_size);
 close:
  close(fd);
 err:
  return res;
}
#endif

/* fetch message and decode if it is base64 or quoted-printable */

int etpan_fetch_message(mailmessage * msg_info,
    struct mailmime * mime_part,
    struct mailmime_single_fields * fields,
    char ** result, size_t * result_len)
{
  char * data;
  size_t len;
  int r;
  int encoding;
  char * decoded;
  size_t decoded_len;
  size_t cur_token;
  int res;
  int encoded;

  encoded = 0;

  r = mailmessage_fetch_section(msg_info,
      mime_part, &data, &len);
  if (r != MAIL_NO_ERROR) {
    res = ERROR_FETCH;
    goto err;
  }

  encoded = 1;

  /* decode message */

  if (encoded) {
    if (fields->fld_encoding != NULL)
      encoding = fields->fld_encoding->enc_type;
    else 
      encoding = MAILMIME_MECHANISM_8BIT;
  }
  else {
    encoding = MAILMIME_MECHANISM_8BIT;
  }

  cur_token = 0;
  r = mailmime_part_parse(data, len, &cur_token,
			  encoding, &decoded, &decoded_len);
  if (r != MAILIMF_NO_ERROR) {
    res = ERROR_FETCH;
    goto free; 
  }

  mailmessage_fetch_result_free(msg_info, data);
  
  * result = decoded;
  * result_len = decoded_len;
  
  return NO_ERROR;
  
 free:
  mailmessage_fetch_result_free(msg_info, data);
 err:
  return res;
}


/* fetch fields */

struct mailimf_fields * fetch_fields(mailmessage * msg_info,
    struct mailmime * mime)
{
  char * data;
  size_t len;
  int r;
  size_t cur_token;
  struct mailimf_fields * fields;

  r = mailmessage_fetch_section_header(msg_info, mime, &data, &len);
  if (r != MAIL_NO_ERROR)
    return NULL;

  cur_token = 0;
  r = mailimf_fields_parse(data, len, &cur_token, &fields);
  if (r != MAILIMF_NO_ERROR) {
    mailmessage_fetch_result_free(msg_info, data);
    return NULL;
  }

  mailmessage_fetch_result_free(msg_info, data);

  return fields;
}



#define MAX_MAIL_COL 72

/* write decoded mailbox */

 int
etpan_mailbox_write(FILE * f, int * col,
    struct mailimf_mailbox * mb)
{
  int r;

  if (* col > 1) {
    
    if (* col + strlen(mb->mb_addr_spec) >= MAX_MAIL_COL) {
      r = mailimf_string_write(f, col, "\r\n ", 3);
      if (r != MAILIMF_NO_ERROR)
	return ERROR_FILE;
      * col = 1;
    }
  }
  
  if (mb->mb_display_name) {
    char * decoded_from;
    size_t cur_token;

    cur_token = 0;
    r = mailmime_encoded_phrase_parse(DEST_CHARSET,
        mb->mb_display_name, strlen(mb->mb_display_name),
        &cur_token, DEST_CHARSET,
        &decoded_from);
    if (r != MAILIMF_NO_ERROR) {
      decoded_from = strdup(mb->mb_display_name);
      if (decoded_from == NULL)
        return ERROR_MEMORY;
    }

    r = mailimf_quoted_string_write(f, col, decoded_from,
        strlen(decoded_from));
    if (r != MAILIMF_NO_ERROR) {
      free(decoded_from);
      return ERROR_FILE;
    }

    if (* col > 1) {
      
      if (* col + strlen(decoded_from) + 3 >= MAX_MAIL_COL) {
	r = mailimf_string_write(f, col, "\r\n ", 3);
	if (r != MAILIMF_NO_ERROR) {
          free(decoded_from);
	  return r;
        }
	* col = 1;
      }
    }

    free(decoded_from);
    
    r = mailimf_string_write(f, col, " <", 2);
    if (r != MAILIMF_NO_ERROR)
      return ERROR_FILE;

    r = mailimf_string_write(f, col,
        mb->mb_addr_spec, strlen(mb->mb_addr_spec));
    if (r != MAILIMF_NO_ERROR)
      return ERROR_FILE;

    r = mailimf_string_write(f, col, ">", 1);
    if (r != MAILIMF_NO_ERROR)
      return ERROR_FILE;
  }
  else {
    r = mailimf_string_write(f, col,
        mb->mb_addr_spec, strlen(mb->mb_addr_spec));
    if (r != MAILIMF_NO_ERROR)
      return ERROR_FILE;
  }


  return NO_ERROR;

}

/* write decoded mailbox list */

int
etpan_mailbox_list_write(FILE * f, int * col,
    struct mailimf_mailbox_list * mb_list)
{
  clistiter * cur;
  int r;
  int first;

  first = 1;

  for(cur = clist_begin(mb_list->mb_list) ; cur != NULL ;
      cur = clist_next(cur)) {
    struct mailimf_mailbox * mb;

    mb = (mailimf_mailbox*) cur->data;

    if (!first) {
      r = mailimf_string_write(f, col, ", ", 2);
      if (r != MAILIMF_NO_ERROR)
	return ERROR_FILE;
    }
    else {
      first = 0;
    }

    r = etpan_mailbox_write(f, col, mb);
    if (r != NO_ERROR)
      return r;
  }

  return NO_ERROR;
}

/* write decoded group */

 int
etpan_group_write(FILE * f, int * col,
    struct mailimf_group * group)
{
  int r;

  r = mailimf_string_write(f, col, group->grp_display_name,
			   strlen(group->grp_display_name));
  if (r != MAILIMF_NO_ERROR)
    return ERROR_FILE;

  r = mailimf_string_write(f, col, ": ", 2);
  if (r != MAILIMF_NO_ERROR)
    return ERROR_FILE;
  
  if (group->grp_mb_list != NULL) {
    r = etpan_mailbox_list_write(f, col, group->grp_mb_list);
    if (r != NO_ERROR)
      return r;
  }

  r = mailimf_string_write(f, col, ";", 1);
  if (r != MAILIMF_NO_ERROR)
    return ERROR_FILE;

  return NO_ERROR;
}

/* write decoded address */

int
etpan_address_write(FILE * f, int * col,
    struct mailimf_address * addr)
{
  int r;

  switch(addr->ad_type) {
  case MAILIMF_ADDRESS_MAILBOX:
    r = etpan_mailbox_write(f, col, addr->ad_data.ad_mailbox);
    if (r != NO_ERROR)
      return r;

    break;

  case MAILIMF_ADDRESS_GROUP:
    r = etpan_group_write(f, col, addr->ad_data.ad_group);
    if (r != NO_ERROR)
      return r;
    
    break;
  }

  return MAILIMF_NO_ERROR;
}

/* write decoded address list */

int
etpan_address_list_write(FILE * f, int * col,
    struct mailimf_address_list * addr_list)
{
  clistiter * cur;
  int r;
  int first;

  first = 1;

  for(cur = clist_begin(addr_list->ad_list) ; cur != NULL ;
      cur = clist_next(cur)) {
    struct mailimf_address * addr;

    addr =( mailimf_address* )  clist_content(cur);

    if (!first) {
      r = mailimf_string_write(f, col, ", ", 2);
      if (r != MAILIMF_NO_ERROR)
	return ERROR_FILE;
    }
    else {
      first = 0;
    }

    r = etpan_address_write(f, col, addr);
    if (r != NO_ERROR)
      return r;
  }

  return NO_ERROR;
}

/* write decoded subject */

 int etpan_subject_write(FILE * f, int * col,
    char * subject)
{
  int r;
  char * decoded_subject;
  size_t cur_token;
  
  r = mailimf_string_write(f, col, "Subject: ", 9);
  if (r != MAILIMF_NO_ERROR) {
    return ERROR_FILE;
  }
  
  cur_token = 0;
  r = mailmime_encoded_phrase_parse(DEST_CHARSET,
      subject, strlen(subject),
      &cur_token, DEST_CHARSET,
      &decoded_subject);
  if (r != MAILIMF_NO_ERROR) {
    decoded_subject = strdup(subject);
    if (decoded_subject == NULL)
      return ERROR_MEMORY;
  }
  
  r = mailimf_string_write(f, col, decoded_subject, strlen(decoded_subject));
  if (r != MAILIMF_NO_ERROR) {
    free(decoded_subject);
    return ERROR_FILE;
  }

  free(decoded_subject);

  r = mailimf_string_write(f, col, "\r\n", 2);
  if (r != MAILIMF_NO_ERROR) {
    return ERROR_FILE;
  }
  * col = 0;

  return NO_ERROR;
}

/* write decoded fields */

int fields_write(FILE * f, int * col,
    struct mailimf_fields * fields)
{
  clistiter * cur;
  int r;
  
  for(cur = clist_begin(fields->fld_list) ; cur != NULL ;
      cur = clist_next(cur)) {
    struct mailimf_field * field;

    field =( mailimf_field* )  clist_content(cur);

    switch (field->fld_type) {
    case MAILIMF_FIELD_FROM:
      r = mailimf_string_write(f, col, "From: ", 6);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_mailbox_list_write(f, col,
          field->fld_data.fld_from->frm_mb_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_REPLY_TO:
      r = mailimf_string_write(f, col, "Reply-To: ", 10);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_address_list_write(f, col,
          field->fld_data.fld_reply_to->rt_addr_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_TO:
      r = mailimf_string_write(f, col, "To: ", 4);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_address_list_write(f, col,
          field->fld_data.fld_to->to_addr_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_CC:
      r = mailimf_string_write(f, col, "Cc: ", 4);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_address_list_write(f, col,
          field->fld_data.fld_cc->cc_addr_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_BCC:
      r = mailimf_string_write(f, col, "Bcc: ", 10);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      if (field->fld_data.fld_bcc->bcc_addr_list != NULL) {
        r = etpan_address_list_write(f, col,
            field->fld_data.fld_bcc->bcc_addr_list);
        if (r != NO_ERROR)
          goto err;
      }

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_SUBJECT:
      r = etpan_subject_write(f, col, field->fld_data.fld_subject->sbj_value);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      break;

    case MAILIMF_FIELD_RESENT_FROM:
      r = mailimf_string_write(f, col, "Resent-From: ", 13);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_mailbox_list_write(f, col,
          field->fld_data.fld_resent_from->frm_mb_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;
      break;

    case MAILIMF_FIELD_RESENT_TO:
      r = mailimf_string_write(f, col, "Resent-To: ", 11);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_address_list_write(f, col,
          field->fld_data.fld_resent_to->to_addr_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;
    case MAILIMF_FIELD_RESENT_CC:
      r = mailimf_string_write(f, col, "Resent-Cc: ", 11);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      r = etpan_address_list_write(f, col,
          field->fld_data.fld_resent_cc->cc_addr_list);
      if (r != NO_ERROR)
        goto err;

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;
    case MAILIMF_FIELD_RESENT_BCC:
      r = mailimf_string_write(f, col, "Resent-Bcc: ", 12);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      
      if (field->fld_data.fld_resent_bcc->bcc_addr_list != NULL) {
        r = etpan_address_list_write(f, col,
            field->fld_data.fld_resent_bcc->bcc_addr_list);
        if (r != NO_ERROR)
          goto err;
      }

      r = mailimf_string_write(f, col, "\r\n", 2);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      * col = 0;

      break;

    case MAILIMF_FIELD_ORIG_DATE:
    case MAILIMF_FIELD_RESENT_DATE:
      r = mailimf_field_write(f, col, field);
      if (r != MAILIMF_NO_ERROR)
        goto err;
      break;

    case MAILIMF_FIELD_OPTIONAL_FIELD:
      if ((strcasecmp(field->fld_data.fld_optional_field->fld_name,
               "X-Mailer") == 0)
          || (strncasecmp(field->fld_data.fld_optional_field->fld_name,
                  "Resent-", 7) == 0)
          || (strcasecmp(field->fld_data.fld_optional_field->fld_name,
                  "Newsgroups") == 0)
          || (strcasecmp(field->fld_data.fld_optional_field->fld_name,
                  "Followup-To") == 0)
          || (strcasecmp(field->fld_data.fld_optional_field->fld_name,
                  "User-Agent") == 0)) {
        r = mailimf_field_write(f, col, field);
        if (r != MAILIMF_NO_ERROR)
          goto err;
      }
      break;

    case MAILIMF_FIELD_MESSAGE_ID:
    case MAILIMF_FIELD_SENDER:
    case MAILIMF_FIELD_IN_REPLY_TO:
    case MAILIMF_FIELD_REFERENCES:
    default:
      break;
    }
  }

  return NO_ERROR;

 err:
  return ERROR_FILE;
}
enum {
  POP3_STORAGE = 0,
  IMAP_STORAGE,
  NNTP_STORAGE,
  MBOX_STORAGE,
  MH_STORAGE,
  MAILDIR_STORAGE,
  FEED_STORAGE
};


int init_storage(struct mailstorage * storage,
    int driver, const char * server, int port,
    int connection_type, const char * user, const char * password, int auth_type, bool xoauth2,
    const char * path, const char * cache_directory, const char * flags_directory)
{
  IdxTy &  r= m_r;
  int cached;

  cached = (cache_directory != NULL);

  switch (driver) {
  case POP3_STORAGE:
    r = pop3_mailstorage_init(storage, server, port, NULL, connection_type,
        auth_type, user, password, cached, cache_directory,
        flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing POP3 storage\n");
      goto err;
    }
    break;

  case IMAP_STORAGE:
    if (xoauth2) {
      r = imap_mailstorage_init_sasl(storage, server, port, NULL, connection_type,
          "xoauth2", NULL, NULL, NULL, NULL, user, password, NULL, cached, cache_directory);
    } else {
      r = imap_mailstorage_init(storage, server, port, NULL, connection_type,
          IMAP_AUTH_TYPE_PLAIN, user, password, cached, cache_directory);
    }
    if (r != MAIL_NO_ERROR) {
      printf("error initializing IMAP storage\n");
      goto err;
    }
    break;

  case NNTP_STORAGE:
    r = nntp_mailstorage_init(storage, server, port, NULL, connection_type,
        NNTP_AUTH_TYPE_PLAIN, user, password, cached, cache_directory,
      flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing NNTP storage\n");
      goto err;
    }
    break;

  case MBOX_STORAGE:
    r = mbox_mailstorage_init(storage, path, cached, cache_directory,
        flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing mbox storage\n");
      goto err;
    }
    break;

  case MH_STORAGE:
    r = mh_mailstorage_init(storage, path, cached, cache_directory,
        flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing MH storage\n");
      goto err;
    }
    break;
  case MAILDIR_STORAGE:
    r = maildir_mailstorage_init(storage, path, cached, cache_directory,
        flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing maildir storage\n");
      goto err;
    }
    break;
  case FEED_STORAGE:
    r = feed_mailstorage_init(storage, path, cached, cache_directory,
        flags_directory);
    if (r != MAIL_NO_ERROR) {
      printf("error initializing feed storage\n");
      goto err;
    }
    break;

 }

  return MAIL_NO_ERROR;

 err:
  return r;
}




 int etpan_render_mime(FILE * f, mailmessage * msg_info,
    struct mailmime * mime)
{
  IdxTy &  r=m_r;
  clistiter * cur;
  int col;
  int text;
  int show;
  struct mailmime_single_fields fields;
  int res;

  mailmime_single_fields_init(&fields, mime->mm_mime_fields,
      mime->mm_content_type);

  text = etpan_mime_is_text(mime);

  r = show_part_info(f, &fields, mime->mm_content_type);
  if (r != NO_ERROR) {
    res = r;
    goto err;
  }

  switch(mime->mm_type) {
  case MAILMIME_SINGLE:
    show = 0;
    if (text)
      show = 1;

    if (show) {
      char * data;
      size_t len;
      char * converted;
      size_t converted_len;
      const char * source_charset;
      size_t write_len;

      /* viewable part */

      r = etpan_fetch_message(msg_info, mime,
          &fields, &data, &len);
      if (r != NO_ERROR) {
        res = r;
        goto err;
      }
    source_charset = fields.fld_content_charset;
      if (source_charset == NULL)
        source_charset = DEST_CHARSET;

      r = charconv_buffer(source_charset, DEST_CHARSET,
          data, len, &converted, &converted_len);
      if (r != MAIL_CHARCONV_NO_ERROR) {

        r = fprintf(f, "[ error converting charset from %s to %s ]\n",
            source_charset, DEST_CHARSET);
          if (r < 0) {
            res = ERROR_FILE;
            goto err;
          }

          write_len = fwrite(data, 1, len, f);
          if (write_len != len) {
            mailmime_decoded_part_free(data);
            res = r;
            goto err;
          }
        }
        else {
          write_len = fwrite(converted, 1, converted_len, f);
          if (write_len != len) {
            charconv_buffer_free(converted);
            mailmime_decoded_part_free(data);
            res = r;
            goto err;
          }

          charconv_buffer_free(converted);
        }

        write_len = fwrite("\r\n\r\n", 1, 4, f);
        if (write_len < 4) {
          mailmime_decoded_part_free(data);
          res = ERROR_FILE;
          goto err;
        }

      mailmime_decoded_part_free(data);
    }
    else {
/* not viewable part */

      r = fprintf(f, "   (not shown)\n\n");
      if (r < 0) {
        res = ERROR_FILE;
        goto err;
      }
    }

    break;

  case MAILMIME_MULTIPLE:

    if (strcasecmp(mime->mm_content_type->ct_subtype, "alternative") == 0) {
      struct mailmime * prefered_body;
      int prefered_score;

      /* case of multiple/alternative */

      /*
        we choose the better part,
        alternative preference :

        text/plain => score 3
        text/xxx   => score 2
        other      => score 1
      */

      prefered_body = NULL;
      prefered_score = 0;

      for(cur = clist_begin(mime->mm_data.mm_multipart.mm_mp_list) ;
          cur != NULL ; cur = clist_next(cur)) {
        struct mailmime * submime;
        int score;

        score = 1;
        submime = (mailmime*) clist_content(cur);
        if (etpan_mime_is_text(submime))
          score = 2;

        if (submime->mm_content_type != NULL) {
          if (strcasecmp(submime->mm_content_type->ct_subtype, "plain") == 0)
            score = 3;
        }
  if (score > prefered_score) {
          prefered_score = score;
          prefered_body = submime;
        }
      }

      if (prefered_body != NULL) {
        r = etpan_render_mime(f, msg_info, prefered_body);
        if (r != NO_ERROR) {
          res = r;
          goto err;
        }
      }
    }
    else {
      for(cur = clist_begin(mime->mm_data.mm_multipart.mm_mp_list) ;
          cur != NULL ; cur = clist_next(cur)) {

        r = etpan_render_mime(f, msg_info, (mailmime*)clist_content(cur));
        if (r != NO_ERROR) {
          res = r;
          goto err;
        }
      }
    }

    break;

  case MAILMIME_MESSAGE:

    if (mime->mm_data.mm_message.mm_fields != NULL) {
      struct mailimf_fields * msg_fields;

      if (msg_info != NULL) {
        msg_fields = fetch_fields(msg_info, mime);
        if (msg_fields == NULL) {
          res = ERROR_FETCH;
          goto err;
        }

        col = 0;
        r = fields_write(f, &col, msg_fields);
        if (r != NO_ERROR) {
          mailimf_fields_free(msg_fields);
          res = r;
        goto err;
        }

        mailimf_fields_free(msg_fields);
      }
      else {
        col = 0;
        r = fields_write(f, &col, mime->mm_data.mm_message.mm_fields);
        if (r != NO_ERROR) {
          res = r;
          goto err;
        }
      }

      r = fprintf(f, "\r\n");
      if (r < 0) {
        res = ERROR_FILE;
        goto err;
      }
    }

    if (mime->mm_data.mm_message.mm_msg_mime != NULL) {
      r = etpan_render_mime(f, msg_info, mime->mm_data.mm_message.mm_msg_mime);
      if (r != NO_ERROR) {
        res = r;
        goto err;
      }
    }

    break;
  }

  return NO_ERROR;

 err:
  return res;
}

