/* > cmdline.c
 * Wrapper for GNU getopt which deals with standard options
 * Copyright (C) 1998-2001 Olly Betts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "cmdline.h"
#include "filename.h"

#include "message.h"

/* It might be useful to be able to disable all long options on small
 * platforms like pre-386 DOS and perhaps some PDAs...
 */
#if 0
# define getopt_long(ARGC, ARGV, STR, OPTS, PTR) getopt(ARGC, ARGV, STR)
#endif

/*
 * bad command line give:
 * <problem>
 *
 * <short syntax>
 *
 * --help gives:
 * <version>
 *
 * <short syntax>
 *
 * <table>
 *
 * <blurb>
 *
 * --version gives:
 * <version>
 */

/*
 * want to cope with optional/required parameters on long options
 * and also parameters on short options
 */

static const char newline_tabs[] = "\n\t\t\t\t";

static int argc;
static char * const *argv;
static const char *shortopts;
static const struct option *longopts;
static int *longind;
static const struct help_msg *help;
static int min_args, max_args;
static const char *args_msg = NULL, *extra_msg = NULL;

static const char *argv0 = NULL;

void
cmdline_help(void)
{
   while (help->opt) {
      const char *longopt = 0;
      int opt = help->opt;
      const struct option *o = 0;

      if (HLP_ISLONG(opt)) {
	 o = longopts + HLP_DECODELONG(opt);
	 longopt = o->name;
	 opt = o->val;
      }

      if (isalnum(opt))
	 printf("  -%c%c", opt, longopt ? ',' : ' ');
      else
	 fputs("     ", stdout);

      if (longopt) {
	 int len = strlen(longopt);
	 printf(" --%s", longopt);
	 if (o && o->has_arg) {
	    const char *p;
	    len += len + 1;

	    if (o->has_arg == optional_argument) {
	       putchar('[');
	       len += 2;
	    }

	    putchar('=');

	    for (p = longopt; *p ; p++) putchar(toupper(*p));

	    if (o->has_arg == optional_argument) putchar(']');
	 }
	 len = (len >> 3) + 2;
	 if (len > 4) len = 0;
	 fputs(newline_tabs + len, stdout);
      } else {
	 fputs(newline_tabs + 1, stdout);
      }

      puts(help->msg);
      help++;
   }
   printf("      --help\t\t\t%s\n", msg(/*display this help and exit*/150));
   printf("      --version\t\t\t%s\n", msg(/*output version information and exit*/151));

   if (extra_msg) {
      putnl();
      puts(extra_msg);
   }

   exit(0);
}

void
cmdline_version(void)
{
   printf("%s - "PACKAGE" "VERSION"\n", argv0);
}

void
cmdline_syntax(void)
{
   printf("\n%s: %s", msg(/*Syntax*/49), argv0);
   if (help->opt) printf(" [%s]...", msg(/*OPTION*/153));
   if (args_msg) {
      putchar(' ');
      puts(args_msg);
      return;
   }
   if (min_args) {
      int i = min_args;
      while (i--) printf(" %s", msg(/*FILE*/124));
   }
   if (max_args == -1) {
      if (!min_args) printf(" [%s]", msg(/*FILE*/124));
      fputs("...", stdout);
   } else if (max_args > min_args) {
      int i = max_args - min_args;
      while (i--) printf(" [%s]", msg(/*FILE*/124));
   }
   putnl();
}

void
cmdline_set_syntax_message(const char *args, const char *extra)
{
   args_msg = args;
   extra_msg = extra;
}

int
cmdline_int_arg(void)
{
   long result;
   char *endptr;

   errno = 0;

   result = strtol(optarg, &endptr, 10);

   if (errno == ERANGE || result > INT_MAX || result < INT_MIN) {
      fprintf(stderr, "%s: ", argv0);
      fprintf(stderr, msg(/*numeric argument `%s' out of range*/185), optarg);
      fputnl(stderr);
      cmdline_syntax();
      exit(1);
   } else if (*optarg == '\0' || *endptr != '\0') {
      fprintf(stderr, "%s: ", argv0);
      fprintf(stderr, msg(/*argument `%s' not an integer*/186), optarg);
      fputnl(stderr);
      cmdline_syntax();
      exit(1);
   }

   return (int)result;
}

double
cmdline_double_arg(void)
{
   double result;
   char *endptr;

   errno = 0;

   result = strtod(optarg, &endptr);

   if (errno == ERANGE) {
      fprintf(stderr, "%s: ", argv0);
      fprintf(stderr, msg(/*numeric argument `%s' out of range*/185), optarg);
      fputnl(stderr);
      cmdline_syntax();
      exit(1);
   } else if (*optarg == '\0' || *endptr != '\0') {
      fprintf(stderr, "%s: ", argv0);
      fprintf(stderr, msg(/*argument `%s' not a number*/187), optarg);
      fputnl(stderr);
      cmdline_syntax();
      exit(1);
   }

   return result;
}

float
cmdline_float_arg(void)
{
   double result = cmdline_double_arg();
   if (fabs(result) > FLT_MAX) {
      fprintf(stderr, "%s: ", argv0);
      fprintf(stderr, msg(/*numeric argument `%s' out of range*/185), optarg);
      fputnl(stderr);
      cmdline_syntax();
      exit(1);
   }
   return (float)result;
}

void
cmdline_init(int argc_, char *const *argv_, const char *shortopts_,
	     const struct option *longopts_, int *longind_,
	     const struct help_msg *help_,
	     int min_args_, int max_args_)
{
   if (!argv0) argv0 = argv_[0];

   argc = argc_;
   argv = argv_;
   shortopts = shortopts_;
   longopts = longopts_;
   longind = longind_;
   help = help_;
   min_args = min_args_;
   max_args = max_args_;
}

int
cmdline_getopt(void)
{
   int opt = getopt_long(argc, argv, shortopts, longopts, longind);

   if (opt == EOF) {
      /* check valid # of args given - if not give syntax message */
      if (argc - optind < min_args) {
	 fprintf(stderr, "%s: %s\n", argv0, msg(/*too few arguments*/122));
	 opt = '?';
      } else if (max_args >= 0 && argc - optind > max_args) {
	 fprintf(stderr, "%s: %s\n", argv0, msg(/*too many arguments*/123));
	 opt = '?';
      }
   }

   switch (opt) {
    case ':': /* parameter missing */
    case '?': /* unknown opt, ambiguous match, or extraneous param */
      /* getopt displays a message for us */
      cmdline_syntax();
      fprintf(stderr, msg(/*Try `%s --help' for more information.\n*/157),
	      argv0);
      exit(1);
    case HLP_VERSION: /* --version */
      cmdline_version();
      exit(0);
    case HLP_HELP: /* --help */
      cmdline_version();
      cmdline_syntax();
      putchar('\n');
      cmdline_help();
      exit(0);
   }
   return opt;
}
