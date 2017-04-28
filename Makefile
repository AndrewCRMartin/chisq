BINDIR = $(HOME)/bin
LIBSRC = $(HOME)/libsrc
LIB    = $(HOME)/lib
INC    = $(HOME)/include

GCC = /usr/bin/gcc -L$(LIB) -I$(INC) -Wall -pedantic -ansi
G++ = /usr/bin/g++ -L$(LIB) -I$(INC) -Wall -pedantic -ansi

EXE = chisq chisig chitab chisq3
OFILES = chisq.o chisig.o chitab.o chisq3.o

all : $(EXE)

chisq : chisq.o
	$(GCC) -o $@ $< -lgen -lm

chisq3 : chisq3.o
	$(GCC) -o $@ $< -lgen -lm

chisig : chisig.o
	$(G++) -o $@ $< -lnumerics -lm

chitab : chitab.o
	$(G++) -o $@ $< -lnumerics -lm

chisq.o : chisq.c
	$(GCC) -c -o $@ $<

.c.o :
	$(G++) -c -o $@ $<

install :
	cp $(EXE) $(BINDIR)

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)
