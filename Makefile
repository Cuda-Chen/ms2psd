DEBUG = 1

CC = gcc
EXEC = ms2psd
#COMMON = -I./libmseed/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid
COMMON = -I/usr/local/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid
CFLAGS =  -Wall
#LDFLAGS = -L./libmseed -Wl,-rpath,./libmseed
#LDLIBS = -Wl,-Bstatic -lmseed -Wl,-Bdynamic -lm -lliquid
LDFLAGS = -L/usr/local
LDLIBS = -lmseed -lm -lliquid -lfftw3

#OBJS = main.o src/parse_miniSEED.o src/bandpass_filter.o src/output2Octave.o src/autocorrelation.o src/spgram.o src/fft.o src/cosine_taper.o src/autocorr.o src/parse_sacpz.o

SRCS := $(sort $(wildcard src/*.c))
OBJS := $(SRCS:%.c=%.o)
OBJS += main.o

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
endif

.PHONY: all clean tests

all: $(EXEC)

$(EXEC): $(OBJS)
	#$(MAKE) -C libmseed/ static
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

tests: $(OBJS)
	@$(MAKE) -C tests

%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

clean:
	#$(MAKE) -C libmseed/ clean
	rm -rf $(OBJS) $(EXEC)
	@$(MAKE) -C tests clean
