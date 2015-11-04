BINDIR = $(HOME)/bin
LIBSRC = $(HOME)/libsrc
LIB    = $(HOME)/lib
INC    = $(HOME)/include

GCC = gcc -L$(LIB) -I$(INC) -Wall -pedantic -ansi
G++ = g++ -L$(LIB) -I$(INC) -Wall -pedantic -ansi

EXE = chisq chisig chitab
OFILES = chisq.o chisig.o chitab.o

all : $(EXE)

chisq : chisq.o
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
