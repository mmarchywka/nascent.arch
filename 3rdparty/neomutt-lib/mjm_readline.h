#ifndef MJM_READLINE_H__
#define MJM_READLINE_H__


#include <stdbool.h>


#include <stdlib.h>
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>

//http://web.mit.edu/gnu/doc/html/rlman_2.html

/* A static variable for holding the line. */
static char *line_read = (char *)NULL;


/* Read a string, and return a pointer to it.  Returns NULL on EOF. */
char *
rl_gets ()
{
const char * _rl_prompt="mjm>";
  /* If the buffer has already been allocated, return the memory
     to the free pool. */
  if (line_read)
    {
      free (line_read);
      line_read = (char *)NULL;
    }

  /* Get a line from the user. */
  line_read = readline (_rl_prompt);

  /* If the line has any text in it, save it on the history. */
  if (line_read && *line_read)
    add_history (line_read);

  return (line_read);
}


#endif // MJM_READLINE_H__



