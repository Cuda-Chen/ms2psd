CC = gcc

COMMON += -I../src/

LDFLAGS = -L/usr/local
LDLIBS = -lmseed -lm -lliquid -lfftw3 

SRCS := $(sort $(wildcard *.c))
BINS := $(SRCS:%.c=%)
OBJS := $(sort $(wildcard ../src/*.o))

.PHONY: all clean

all: $(BINS)

$(BINS): % : %.c
	@printf 'Building $<\n';
	$(CC) $(COMMON) $(CFLAGS) -o $@ $< $(OBJS) $(LDFLAGS) $(LDLIBS)

clean: 
	@rm -rf *.o $(BINS) *.dSYM
