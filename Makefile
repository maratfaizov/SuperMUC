####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = gcc
CFLAGS = -g -O0 -std=c99
LIBS = -lm

# Make sure you have loaded the papi module before uncommenting these and include papi.h in the sources
#CFLAGS += $(PAPI_INC)

#LIBS += $(PAPI_LIB)

CFLAGS += -I/usr/local/include
LIBS += /usr/local/lib/libpapi.a

SRCS = xread.c xwrite.c gccg.c vol2mesh.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: $(OBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf *.o gccg 
