*  
* $Id: gzopen.F,v 1.2 2002/11/12 11:51:30 hvogt Exp $
*  
* $Log: gzopen.F,v $
* Revision 1.2  2002/11/12 11:51:30  hvogt
* gzio modified for several I/O streams
*
* Revision 1.1.1.1  2002/06/05 16:48:52  hvogt
* Initial Installation
*
*  
      SUBROUTINE GZOPEN (LUN, MODE, FILNAM, ISTAT)
*
************************************************************************
*
*  Opens a file for zipped IO
*
*  Arguments:
*     LUN        logical unit number used for this file
*     MODE       string selecting IO mode (r, w)
*     FILNAM     name of the file (including the null-terminator)
*     ISTAT      status, =zero if success, #0 otherwise
*
*  Input  file: Extension '.gz' is checked
*  Output file: Extension '.gz' is checked
*
************************************************************************
*
      IMPLICIT NONE
#include "gzio.inc"
*
      CHARACTER*1   MODE, MYMODE
      CHARACTER*(*) FILNAM
      CHARACTER*256 MYFNAM
      INTEGER       LUN, ILUN, GZFINDL, FLNLEN, ISTAT, LENOCC
*
      INTEGER       NLUNS
      DATA          NLUNS / 0/
ckr      SAVE          NLUNS
*
      ISTAT = 0
      ILUN = GZFINDL (LUN)
      IF (ILUN .NE. 0)  THEN
        WRITE (STDLUN,'(A,I4)')
     &         'GZOPEN: Error - LUN used already, LUN =', LUN
        ISTAT = 1
        GO  TO  90
      END IF
*
      IF (NLUNS+1 .GT. LUNMAX)  THEN
        WRITE (STDLUN,'(A,I4)')
     &         'GZOPEN: Error - LUN table overflow for LUN =', LUN
        ISTAT = 1
        GO  TO  90
      END IF
*
      FLNLEN = LENOCC (FILNAM)
      IF (MODE .EQ. 'r' .OR. MODE .EQ. 'R')  THEN
        MYMODE='r'
      ELSE IF (MODE .EQ. 'w' .OR. MODE .EQ. 'W')  THEN
        MYMODE='w'
      ELSE
        WRITE (STDLUN,'(A,I4)')
     &         'GZOPEN: Error - undefined mode for LUN =', LUN
      END IF
*
      IF (MYMODE .EQ. 'r')  THEN
        IF (FILNAM(FLNLEN-2:FLNLEN) .EQ. '.gz')  THEN
          NLUNS = NLUNS + 1
          LUNS(NLUNS) = LUN
          LUNMODS(NLUNS) = 1
          CALL GZIOOP (FILDES(1,NLUNS), 'r', FILNAM(1:FLNLEN), FLNLEN,
     &                 ISTAT)
          GO  TO  90
        ELSE
          WRITE (STDLUN,'(A,A,I4)')
     &          'GZOPEN: Error - Input file name has no .gz extension',
     &          ' LUN =', LUN
          ISTAT = 1
          GO  TO  90
        END IF
*
      ELSE IF (MYMODE .EQ. 'w')  THEN
        IF (FILNAM(FLNLEN-2:FLNLEN) .EQ. '.gz')  THEN
          NLUNS = NLUNS + 1
          LUNS(NLUNS) = LUN
          LUNMODS(NLUNS) = 2
          CALL GZIOOP (FILDES(1,NLUNS), 'w', FILNAM(1:FLNLEN), FLNLEN,
     &                 ISTAT)
          GO  TO  90
        ELSE
          WRITE (STDLUN,'(A,A,I4)')
     &          'GZOPEN: Error - Output file name has no .gz extension',
     &          ' LUN =', LUN
          ISTAT = 1
        END IF
      END IF
*
   90 RETURN
*
      END
*
*  
* $Id: gzclose.F,v 1.4 2003/01/21 14:51:39 hvogt Exp $
*  
* $Log: gzclose.F,v $
* Revision 1.4  2003/01/21 14:51:39  hvogt
* reset lun table if lun closed successfully
*
* Revision 1.3  2003/01/21 14:47:28  hvogt
* reset lun table if lun closed successfully
*
* Revision 1.2  2002/11/12 11:51:21  hvogt
* gzio modified for several I/O streams
*
* Revision 1.1.1.1  2002/06/05 16:48:52  hvogt
* Initial Installation
*
*  
      SUBROUTINE GZCLOSE (LUN, ISTAT)
*
*
************************************************************************
*
*  Closes the zipped file
*
*  Argument:
*     LUN        Logical Unit Number
*     ISTAT      =0 ok, =1 error occured
*
************************************************************************
*
      IMPLICIT NONE
#include "gzio.inc"
*
      INTEGER       LUN, ISTAT
      INTEGER       ILUN, GZFINDL
*
      ISTAT = 0
      ILUN = GZFINDL (LUN)
      IF (ILUN .EQ. 0)  THEN
        WRITE (STDLUN,'(A,I4)')
     &         'GZCLOSE: Error - undefined LUN =', LUN
        ISTAT = 1
        GO  TO  90
      END IF
*
      IF (LUNMODS(ILUN) .EQ. 1)  THEN
        CALL GZCLOS (FILDES(1,ILUN),'r',ISTAT)
        IF (ISTAT .NE. 0)  THEN
          WRITE (STDLUN,'(A,A,I4)')
     &          'GZCLOSE: Error - Input file could not be closed',
     &          'LUN =',LUN
        ELSE
          LUNS(ILUN) = 0
        END IF
      ELSE 
        CALL GZCLOS (FILDES(1,ILUN),'w',ISTAT)
        IF (ISTAT .NE. 0)  THEN
          WRITE (STDLUN,'(A,A,I4)')
     &          'GZCLOSE: Error - Output file could not be closed',
     &          'LUN =',LUN
        ELSE
          LUNS(ILUN) = 0
        END IF
      END IF
*
   90 RETURN
*
      END
*
*  
* $Id: gzread.F,v 1.3 2003/05/13 15:04:40 hvogt Exp $
*  
* $Log: gzread.F,v $
* Revision 1.3  2003/05/13 15:04:40  hvogt
* exit replaced by CERN/packlib's exitf
*
* Revision 1.2  2002/11/12 11:52:07  hvogt
* gzio modified for several I/O streams
*
* Revision 1.1.1.1  2002/06/05 16:48:52  hvogt
* Initial Installation
*
*  
      SUBROUTINE GZREAD (LUN, CHBUF, NBBUF, ISTAT)
************************************************************************
*
*  Reads a line from a zipped ascii file
*
*  Arguments:
*     LUN        Logical Unit Number for reading
*     CHBUF      array for storing unzipped info
*     NBBUF      record size, number of bytes of array MBUF
*     ISTAT      status, =NBREC if ok., =zero if EOF encountered
************************************************************************
*
      IMPLICIT NONE
#include "gzio.inc"
*
      CHARACTER*(*)   CHBUF
      INTEGER         LUN, NBBUF, ISTAT
      INTEGER         ILUN, GZFINDL
*
      ILUN = GZFINDL (LUN)
      IF (ILUN .EQ. 0)  THEN
        WRITE (STDLUN,'(A,I4)')
     &         'GZREAD: Error - undefined LUN =', LUN
        CALL EXITF (1)
      END IF
*
      IF (LUNMODS(ILUN) .EQ. 1)  THEN
        CALL GZGETS (FILDES(1,ILUN),NBBUF,CHBUF,ISTAT)
      ELSE
        WRITE (STDLUN,'(A,A,I4)')
     &         'GZREAD: Error - attempt to read from a file',
     &         'opened for writing, LUN =', LUN
        CALL EXITF (1)
      END IF
*
      END
*
*  
* $Id: gzfindl.F,v 1.2 2003/05/13 15:04:40 hvogt Exp $
*  
* $Log: gzfindl.F,v $
* Revision 1.2  2003/05/13 15:04:40  hvogt
* exit replaced by CERN/packlib's exitf
*
* Revision 1.1  2002/11/12 11:50:39  hvogt
* gzio modified for several I/O streams
*
*  
      INTEGER FUNCTION GZFINDL (LUN)
************************************************************************
*
*  check if LUN is defined and returns the correspondig pointer
*  as function value
*
*  Arguments:
*     LUN        Logical Unit Number for reading
*     ILUN       Pointer to the file descriptor area
************************************************************************
*
      IMPLICIT NONE
#include "gzio.inc"
*
      INTEGER         LUN, J
*
      IF (LUN .LT. 1 .OR. LUN .GT. 99)  THEN
        WRITE (STDLUN,'(A,I4)') 
     &         'GZFINDL: Error - LUN outside limits LUN =', LUN
        CALL EXITF (1)
      END IF
*
      DO  J = 1,LUNMAX
        IF (LUN .EQ. LUNS(J))  THEN
          GZFINDL = J
          RETURN
        END IF
      END DO
      GZFINDL = 0
*
      END
