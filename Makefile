CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=./bin

all:fqvendorfail

fqvendorfail:fqvendorfail.c seqtk-1.3/kseq.h
		$(CC) $(CFLAGS) fqvendorfail.c -o $@ -lz -lm

#install:all
#		install fqvendorfail $(BINDIR)

clean:
		rm -fr gmon.out *.o ext/*.o a.out fqvendorfail *~ *.a *.dSYM session*
