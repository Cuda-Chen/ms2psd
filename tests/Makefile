DEBUG = 1
CC = gcc
CFLAGS = 
COMMON += -I../src/ -I../include/ -I/usr/local/include/liquid -I/usr/local/include/ -I/usr/local/

LDFLAGS = -L/usr/local
LDLIBS = -lmseed -lm -lliquid -lfftw3 

SRCS := $(sort $(wildcard *.c))
BINS := $(SRCS:%.c=%)
OBJS := $(sort $(wildcard ../src/*.o))

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
endif

OS := $(shell uname -s)
# Add library path for macOS/Darwin as libmseed stores differently
ifeq ($(OS), Darwin)
LDFLAGS += -L/usr/local/lib
endif

.PHONY: all clean

all: $(BINS)

$(BINS): % : %.c
	@printf 'Building $<\n';
	$(CC) $(COMMON) $(CFLAGS) -o $@ $< $(OBJS) $(LDFLAGS) $(LDLIBS)

clean: 
	@rm -rf *.o $(BINS) *.dSYM
