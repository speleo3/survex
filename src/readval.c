/* > readval.c
 * Routines to read a prefix or number from the current input file
 * Copyright (C) 1991-2001 Olly Betts
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
# include <config.h>
#endif

#include "cavern.h"
#include "debug.h"
#include "filename.h"
#include "message.h"
#include "readval.h"
#include "datain.h"
#include "netbits.h"
#include "osalloc.h"
#include "str.h"

#ifdef HAVE_SETJMP_H
# define LONGJMP(JB) longjmp((JB), 1)
#else
# define LONGJMP(JB) exit(1)
#endif

int root_depr_count = 0;

/* Dinky macro to handle any case forcing needed */
#define docase(X) (pcs->Case == OFF ? (X) :\
                   (pcs->Case == UPPER ? toupper(X) : tolower(X)))

/* if prefix is omitted: if fOmit return NULL, otherwise use longjmp */
static prefix *
read_prefix_(bool fOmit, bool fSurvey, bool fSuspectTypo, bool fAllowRoot)
{
   prefix *back_ptr, *ptr;
   char *name;
   size_t name_len = 32;
   size_t i;
   bool fNew;
   bool fImplicitPrefix = fTrue;
   int depth = -1;

   skipblanks();
   if (isRoot(ch)) {
      if (!fAllowRoot) {
	 compile_error(/*ROOT is deprecated*/25);
	 skipline();
	 LONGJMP(file.jbSkipLine);	 
      }
      if (root_depr_count < 5) {
	 compile_warning(/*ROOT is deprecated*/25);
	 if (++root_depr_count == 5)
	    compile_warning(/*No further uses of this deprecated feature will be reported*/95);
      }
      nextch();
      ptr = root;
      if (!isNames(ch)) {
	 if (!isSep(ch)) return ptr;
	 nextch();
      }
      fImplicitPrefix = fFalse;
   } else {
      ptr = pcs->Prefix;
   }

   i = 0;
   do {
      fNew = fFalse;
      if (isSep(ch)) fImplicitPrefix = fFalse;
      /* get a new name buffer for each level */
      name = osmalloc(name_len);
      /* i==0 iff this is the first pass */
      if (i) {
	 i = 0;
	 nextch();
      }
      while (isNames(ch)) {
	 if (i < pcs->Truncate) {
	    /* truncate name */
	    name[i++] = docase(ch);
	    if (i >= name_len) {
	       name_len = name_len + name_len;
	       name = osrealloc(name, name_len);
	    }
	 }
	 nextch();
      }
      if (i == 0) {
	 osfree(name);
	 if (!fOmit) {
	    if (isEol(ch)) {
	       if (fSurvey) {
		  compile_error(/*Expecting survey name*/89);
	       } else {
		  compile_error(/*Expecting station name*/28);
	       }
	    } else {
	       compile_error(/*Character `%c' not allowed in station name (use *SET NAMES to set allowed characters)*/7, ch);
	    }
	    skipline();
	    LONGJMP(file.jbSkipLine);
	 }
	 return (prefix *)NULL;
      }

      name[i++] = '\0';
      name = osrealloc(name, i);

      back_ptr = ptr;
      ptr = ptr->down;
      if (ptr == NULL) {
	 /* Special case first time around at each level */
	 ptr = osnew(prefix);
	 ptr->ident = name;
	 ptr->right = ptr->down = NULL;
	 ptr->pos = NULL;
	 ptr->shape = 0;
	 ptr->stn = NULL;
	 ptr->up = back_ptr;
	 ptr->filename = NULL;
	 ptr->min_export = ptr->max_export = 0;
	 ptr->sflags = BIT(SFLAGS_SURVEY);
	 if (fSuspectTypo && !fImplicitPrefix)
	    ptr->sflags |= BIT(SFLAGS_SUSPECTTYPO);
	 back_ptr->down = ptr;
	 fNew = fTrue;
#ifdef NEW3DFORMAT
	 if (fUseNewFormat) {
	    create_twig(ptr,file.filename);
	 }
#endif
      } else {
	 prefix *ptrPrev = NULL;
	 int cmp = 1; /* result of strcmp ( -ve for <, 0 for =, +ve for > ) */
	 while (ptr && (cmp = strcmp(ptr->ident, name))<0) {
	    ptrPrev = ptr;
	    ptr = ptr->right;
	 }
	 if (cmp) {
	    /* ie we got to one that was higher, or the end */
	    prefix *newptr = osnew(prefix);
	    newptr->ident = name;
	    if (ptrPrev == NULL)
	       back_ptr->down = newptr;
	    else
	       ptrPrev->right = newptr;
	    newptr->right = ptr;
	    newptr->down = NULL;
	    newptr->pos = NULL;
	    newptr->shape = 0;
	    newptr->stn = NULL;
	    newptr->up = back_ptr;
	    newptr->filename = NULL;
  	    newptr->min_export = newptr->max_export = 0;
	    newptr->sflags = BIT(SFLAGS_SURVEY);
	    if (fSuspectTypo && !fImplicitPrefix)
	       newptr->sflags |= BIT(SFLAGS_SUSPECTTYPO);
	    ptr = newptr;
	    fNew = fTrue;
#ifdef NEW3DFORMAT
	    if (fUseNewFormat) {
	       create_twig(ptr,file.filename);
	    }
#endif
	 }
      }
      depth++;
      fOmit = fFalse; /* disallow after first level */
   } while (isSep(ch));
   /* don't warn about a station that is refered to twice */
   if (!fNew) ptr->sflags &= ~BIT(SFLAGS_SUSPECTTYPO);
   
   if (fNew) {
      /* fNew means SFLAGS_SURVEY is currently set */
      ASSERT(TSTBIT(ptr->sflags, SFLAGS_SURVEY));
      if (!fSurvey) ptr->sflags &= ~BIT(SFLAGS_SURVEY);
   } else {
      /* check that the same name isn't being used for a survey and station */
      if (fSurvey ^ TSTBIT(ptr->sflags, SFLAGS_SURVEY)) {
	 compile_error(/*`%s' can't be both a station and a survey*/27,
		       sprint_prefix(ptr));
      }
   }

   /* check the export level */
#if 0
   printf("R min %d max %d depth %d pfx %s\n",
	  ptr->min_export, ptr->max_export, depth, sprint_prefix(ptr));
#endif
   if (ptr->min_export == 0) {
      if (depth > ptr->max_export) ptr->max_export = depth;
   } else if (ptr->max_export < depth) {
      const char *filename_store = file.filename;
      unsigned int line_store = file.line;
      prefix *survey = ptr;
      char *s;
      int level;
      for (level = ptr->max_export + 1; level; level--) {
	 survey = survey->up;
	 ASSERT(survey);
      }
      s = osstrdup(sprint_prefix(survey));
      if (survey->filename) {
	 file.filename = survey->filename;
	 file.line = survey->line;
      }
      compile_error(/*Station `%s' not exported from survey `%s'*/26,
		    sprint_prefix(ptr), s);
      if (survey->filename) {
	 file.filename = filename_store;
	 file.line = line_store;
      }
      osfree(s);
#if 0
      printf(" *** pfx %s warning not exported enough depth %d "
             "ptr->max_export %d\n", sprint_prefix(ptr),
	     depth, ptr->max_export);
#endif
   }
   return ptr;
}

/* if prefix is omitted: if fOmit return NULL, otherwise use longjmp */
extern prefix *
read_prefix_survey(bool fOmit, bool fAllowRoot)
{
   return read_prefix_(fOmit, fTrue, fFalse, fAllowRoot);
}

/* if prefix is omitted: if fOmit return NULL, otherwise use longjmp */
extern prefix *
read_prefix_stn(bool fOmit, bool fAllowRoot)
{
   return read_prefix_(fOmit, fFalse, fFalse, fAllowRoot);
}

/* if prefix is omitted: if fOmit return NULL, otherwise use longjmp */
/* Same as read_prefix_stn but implicit checks are made */
extern prefix *
read_prefix_stn_check_implicit(bool fOmit, bool fAllowRoot)
{
   return read_prefix_(fOmit, fFalse, fTrue, fAllowRoot);
}

/* if numeric expr is omitted: if fOmit return HUGE_REAL, else longjmp */
extern real
read_numeric(bool fOmit)
{
   bool fPositive, fDigits = fFalse;
   real n = (real)0.0;
   filepos fp;
   int ch_old;

   skipblanks();
   get_pos(&fp);
   ch_old = ch;
   fPositive = !isMinus(ch);
   if (isSign(ch)) nextch();

   while (isdigit(ch)) {
      n = n * 10.0f + (char)(ch - '0');
      nextch();
      fDigits = fTrue;
   }

   if (isDecimal(ch)) {
      real mult = (real)1.0;
      nextch();
      while (isdigit(ch)) {
	 mult *= (real).1;
	 n += (char)(ch - '0') * mult;
	 fDigits = fTrue;
	 nextch();
      }
   }

   /* !'fRead' => !fDigits so fDigits => 'fRead' */
   if (fDigits) return (fPositive ? n : -n);

   /* didn't read a valid number.  If it's optional, reset filepos & return */
   if (fOmit) {
      set_pos(&fp);
      return HUGE_REAL;
   }

   if (isOmit(ch_old)) {
      compile_error(/*Field may not be omitted*/8);
   } else {
      compile_error(/*Expecting numeric field*/9);
   }
   showandskipline(NULL, 1);
   LONGJMP(file.jbSkipLine);
   return 0.0; /* for brain-fried compilers */
}

/* read numeric expr or omit (return HUGE_REAL); else longjmp */
extern real
read_numeric_or_omit(void)
{
   real v = read_numeric(fTrue);
   if (v == HUGE_REAL) {
      if (!isOmit(ch)) {
	 compile_error(/*Expecting numeric field*/9);
	 showandskipline(NULL, 1);
	 LONGJMP(file.jbSkipLine);
	 return 0.0; /* for brain-fried compilers */
      }
      nextch();
   }
   return v;
}

extern unsigned int
read_uint(void)
{
   unsigned int n = 0;
   skipblanks();
   if (!isdigit(ch)) {
      compile_error(/*Expecting numeric field*/9);
      showandskipline(NULL, 1);
      LONGJMP(file.jbSkipLine);
   }
   while (isdigit(ch)) {
      n = n * 10 + (char)(ch - '0');
      nextch();
   }
   return n;
}

extern void
read_string(char **pstr, int *plen)
{
   s_zero(pstr);

   skipblanks();
   if (ch == '\"') {
      nextch();
      while (1) {
	 if (isEol(ch)) {
	    compile_error(/*Missing &quot;*/69);
	    skipline();
	    return;
	 }

	 if (ch == '\"') break;

	 s_catchar(pstr, plen, ch);
	 nextch();
      }
   } else {
      while (1) {
	 if (isEol(ch) || isComm(ch)) {
	    if (!*pstr || !(*pstr)[0]) {
	       compile_error(/*Expecting string field*/121);
	       skipline();
	    }
	    return;
	 }

	 if (isBlank(ch)) break;

	 s_catchar(pstr, plen, ch);
	 nextch();
      }
   }
   
   nextch();
}
