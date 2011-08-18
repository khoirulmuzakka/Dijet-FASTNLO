###########################################################################################
#
#  simple Makefile for FastNLOReader
#
###########################################################################################

OBJS 	= ReadFNLOTable.o FastNLOBlockB.o FastNLOReader.o qcdnumdummy.o
CC 	= gcc
CFLAGS 	= -Wall -g
LDFLAGS	= -lm -lpthread `gtk-config --cflags` `gtk-config --libs` -lgthread

ReadFNLOTable: $(OBJS)	
	$(CC) $(OBJS) \
	-L/afs/desy.de/group/alliance/mcg/public/MCGenerators/lhapdf/5.8.4/i586_rhel50/lib/ -lLHAPDF \
	-o ReadFNLOTable

clean:
	rm -f *.o *.a *.so *Dict* 2>/dev/null

%.o: %.c

all: ReadFNLOTable
