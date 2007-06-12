*  
* $Id: gzio.inc,v 1.1 2002/11/12 11:50:02 hvogt Exp $
*  
* $Log: gzio.inc,v $
* Revision 1.1  2002/11/12 11:50:02  hvogt
* gzio modified for several I/O streams
*
*  
      INTEGER         LUNMAX, STDLUN
      PARAMETER       (LUNMAX=10, STDLUN=6)
      COMMON /GZIOFD/ LUNS(LUNMAX), LUNMODS(LUNMAX), 
     &                FILDES(8,LUNMAX)
      INTEGER         LUNS, LUNMODS, FILDES
*
