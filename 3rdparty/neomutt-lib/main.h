/**
 * @file
 * Command line processing
 *
 * @authors
 * Copyright (C) 2018 Richard Russon <rich@flatcap.org>
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

#ifndef MUTT_MAIN_H
#define MUTT_MAIN_H

#include <stdbool.h>

/* These Config Variables are only used in main.c */
extern bool C_ResumeEditedDraftFiles;
extern  bool mjm_daemon_mode;
extern  bool mjm_init_only;

extern int shutdown_mutt();
#endif /* MUTT_MAIN_H */
