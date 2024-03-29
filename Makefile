DEBUG = 0
SHOWEACHTRACE = 0

CC = gcc
EXEC = ms2psd
#COMMON = -I./libmseed/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid
COMMON = -I/usr/local/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid -I/usr/local/include/
CFLAGS =  -Wall
#LDFLAGS = -L./libmseed -Wl,-rpath,./libmseed
#LDLIBS = -Wl,-Bstatic -lmseed -Wl,-Bdynamic -lm -lliquid
LDFLAGS = -L/usr/local/
LDLIBS = -lmseed -lm -lliquid -lfftw3

#OBJS = main.o src/parse_miniSEED.o src/bandpass_filter.o src/output2Octave.o src/autocorrelation.o src/spgram.o src/fft.o src/cosine_taper.o src/autocorr.o src/parse_sacpz.o

SRCS := $(sort $(wildcard src/*.c))
OBJS := $(SRCS:%.c=%.o)
OBJS += main.o

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
else
CFLAGS += -O1
endif

OS := $(shell uname -s)
# Add library path for macOS/Darwin as libmseed stores differently
ifeq ($(OS), Darwin)
LDFLAGS += -L/usr/local/lib
endif

.PHONY: all clean tests format

all: $(EXEC)

$(EXEC): $(OBJS)
	#$(MAKE) -C libmseed/ static
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

tests: $(OBJS)
	@$(MAKE) -C tests

%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

format:
	clang-format -i main.c src/*.[hc] tests/*.[hc]

clean:
	#$(MAKE) -C libmseed/ clean
	rm -rf $(OBJS) $(EXEC)
	@$(MAKE) -C tests clean
