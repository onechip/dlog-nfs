CC		= gcc
CPP		= g++
LINK		= g++

NTLPREFIX	= /usr/local
NTLLIB		= -L$(NTLPREFIX)/lib -lntl -lgmp
NTLINCLUDE	= -I$(NTLPREFIX)/include

PREFFLAGS	=
CPPFLAGS	= $(PREFFLAGS) -O2 -Wall -Wno-deprecated $(NTLINCLUDE) -I.
LINKFLAGS	= $(NTLLIB)

ALL_PROGS	= dlog-test
COMMON_OBJS	= pair_ZZ_long.o ZZFactoring.o \
		  mat_long.o svector.o smatrix.o \
		  FactorBase.o AlgebraicFactorBase.o \
		  DiscreteLog.o IndexCalculus.o NumberFieldSieve.o
COMMON_HEADERS	=


.SUFFIXES: .cpp .o

.cpp.o:
	$(CPP) -c $(CPPFLAGS) $*.cpp


all:	$(ALL_PROGS)

dlog-test:	dlog-test.o $(COMMON_OBJS)
	$(LINK) -o $@ $^ $(LINKFLAGS)

clean:
	rm -f *% *~ *.o core a.out $(ALL_PROGS)

