###########################################################################################
#
#  simple Makefile for FastNLOReader
#
###########################################################################################

OBJS 	= ReadFNLOTable.o FastNLOBlockB.o FastNLOReader.o qcdnumdummy.o
CC 	= gcc
CFLAGS 	= -Wall -g
LDFLAGS	= -lm -lpthread `gtk-config --cflags` `gtk-config --libs` -lgthread

LHAPDFVERSION	= 5.8.5
LHAPDFDIR	= /afs/desy.de/group/alliance/mcg/public/MCGenerators/lhapdf
LHAPDFLIBPATH	= $(LHAPDFDIR)/$(LHAPDFVERSION)/$(SYSNAME)/lib/
LHAPDFLIB	= LHAPDF

ReadFNLOTable: $(OBJS)	
	$(CC) $(OBJS) \
	-L$(LHAPDFLIBPATH) -l$(LHAPDFLIB) \
	-o ReadFNLOTable

clean:
	rm -f *.o *.a *.so *Dict* 2>/dev/null

%.o: %.c

all: ReadFNLOTable
