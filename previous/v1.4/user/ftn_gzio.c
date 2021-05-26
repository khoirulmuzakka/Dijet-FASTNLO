/*
*
* $Id: ftn_gzio.c,v 1.2 2004/05/04 07:41:58 hvogt Exp $
*
* $Log: ftn_gzio.c,v $
* Revision 1.2  2004/05/04 07:41:58  hvogt
* modified for x86_64 Architectures
*
* Revision 1.1.1.1  2002/06/05 16:48:52  hvogt
* Initial Installation
*
*
*/
/* ftn_gzio - usage of the zlib compression library by fortran calls
 * Copyright (C) 2001 Harald Vogt
 * derived from gzio
 * Copyright (C) 1995-1998 Jean-loup Gailly.
 */


#include <stdio.h>
#include <string.h>
#include "zlib.h"

extern void exit  OF((int));


/* ===========================================================================
 * fortran interface to open .gz files for read/write
 *
 *    SUBROTINE GZIOOP (FILDES, MODE, FILENAME, LGFNAME, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       MODE       string selecting IO mode (r, w)
 *       FILENAME   name of the file (including the null-terminator)
 *       LGFNAME    length of the filename string
 *      *ISTAT      status, =zero if success
 */

void gzioop_(gzFile *fildes,char *fmode,char *fname,int *lgfname,int *stat)
{
   char *pttext;
   gzFile file;

   *stat   = -1;

   if(*lgfname==0) {
        perror("GZIOOP: no filename");
        *stat = 107;
        return;
   }

   if(!( pttext=(char*)malloc(*lgfname+8))) {
         fprintf(stderr, "GZIOOP: no malloc possible\n");
         exit(1);
   }
   strncpy(pttext,fname,*lgfname);
   pttext[*lgfname]='\0';

   switch (*fmode) {
     case 'r':
           file = gzopen(pttext, "rb");
           if (file == NULL) {
               fprintf(stderr, "GZIOOP: input file open error\n");
               *stat = 108;
               return;
           }
           break;
     case 'w':
           file = gzopen(pttext, "wb");
           if (file == NULL) {
               fprintf(stderr, "GZIOOP: output file oopen error\n");
               *stat = 108;
               return;
           }
           break;
     default:  fprintf(stderr,"GZIOOP: unknown i/o mode %s\n",fmode);
               *stat = 108;
               return;
   }
/* myprint_(file);  */
   *fildes = file;
   *stat   = 0;
   return;
}

/*
 *
 *    SUBROTINE GZCLOS (FILDES, MODE, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       MODE       string selecting IO mode (r, w)
 *      *ISTAT      status, =zero if success
 */

void gzclos_(gzFile *fildes,char *fmode,int *stat)
{
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   switch (*fmode) {
     case 'r':
           gzclose(file);
           if (file == NULL) {
               fprintf(stderr, "GZCLOS: no file defined to close\n");
               exit(1);
           }
           break;
     case 'w':
/*           gzseek(file, 1L, SEEK_CUR);  */ /* add one zero byte */
           gzclose(file);
           if (file == NULL) {
               fprintf(stderr, "GZCLOS: no file defined to close\n");
               exit(1);
           }
           break;
     default:  fprintf(stderr,"GZCLOS: unknown i/o mode %s\n",fmode);
               exit(1);
   }
   *stat   = 0;
   return;
}

/*
 *
 *    SUBROTINE GZREWIND (FILDES)
 *
 *       FILDES     address of file descriptor area
 */

void gzrewind_(gzFile *fildes)
{
   gzFile file;

   file    = *fildes;
   gzrewind(file);
   return;
}

/*
 *
 *    SUBROTINE GZPUTS (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes to be written
 *       MBUF       vector to be written
 *      *ISTAT      status, =zero if success
 */

void gzputs_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   int  nbdn, nbdo;
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   nbdo   = *nbrec;

   if (file == NULL) {
       fprintf(stderr, "GZPUTS: no file defined to write\n");
       exit(1);
   }
   nbdn    = gzwrite(file,  mbuf, nbdo);
   *stat   = nbdn;
   return;
}

/*
 *
 *    SUBROTINE GZGETS (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes of array MBUF
 *       MBUF       array for storing unzipped info
 *      *ISTAT      status, =NBREC if ok., =zero if EOF encountered
 */

void gzgets_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   unsigned len;
   char *buf;
   gzFile file;

   *stat   = 0;
   file    = *fildes;

   len = *nbrec;
   buf = mbuf;
   if (file == NULL) {
       fprintf(stderr, "GZGETS: no file defined to write\n");
       exit(1);
   }
   if (buf == Z_NULL || len <= 0) *stat = 0;
   while (--len > 0 && gzread(file, buf, 1) == 1 && *buf++ != '\n') ;
   *buf = '\0';
   *stat   = buf - mbuf;

   return;
}

/*
 *
 *    SUBROTINE GZGETB (FILDES, NBREC, MBUF, ISTAT)
 *
 *       FILDES     address of file descriptor area
 *       NBREC      record size, number of bytes of array MBUF
 *       MBUF       array for storing unzipped info
 *      *ISTAT      status, =NBREC if ok., =zero if EOF encountered
 */

void gzgetb_(gzFile *fildes,int *nbrec, char *mbuf, int *stat)
{
   int  nbdn, nbdo;
   gzFile file;

   *stat   = -1;
   file    = *fildes;

   nbdo   = *nbrec;

   if (file == NULL) {
       fprintf(stderr, "GZGETB: Error - no input file defined\n");
       exit(1);
   }
   nbdn    = gzread(file,  mbuf, nbdo);
   *stat   = nbdn;
   return;
}

