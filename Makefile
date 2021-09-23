MF=	Makefile

#COMPILER

# On ARCHER2, cc is a wrapper for whichever C compiler - GNU (gcc),
# AMD (aocc), or Cray (craycc) - has been chosen by loading the
# appropriate PrgEnv module prior to compilation. On other systems,
# you will probably have to use "mpicc".

CC=	cc
CFLAGS=	-g
LFLAGS=	-lm


EXE=	fractal

SRC= \
	arralloc.c \
	fractal.c \
	read_options.c \
	write_ppm.c

INC = \
	arralloc.h \
	read_options.h \
	write_ppm.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) $(INC)

clean:
	rm -f $(OBJ) $(EXE) output.ppm core
