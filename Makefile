DEBUG = 1

CC = gcc
EXEC = ms2filter
#COMMON = -I./libmseed/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid
COMMON = -I/usr/local/ -Iinclude/ -Isrc/ -I/usr/local/include/liquid
CFLAGS =  -Wall
#LDFLAGS = -L./libmseed -Wl,-rpath,./libmseed
#LDLIBS = -Wl,-Bstatic -lmseed -Wl,-Bdynamic -lm -lliquid
LDFLAGS = -L/usr/local
LDLIBS = -lmseed -lm -lliquid

OBJS = main.o src/parse_miniSEED.o src/bandpass_filter.o src/write2script.o src/write2file.o src/write2mseed.o

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
endif

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	#$(MAKE) -C libmseed/ static
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

clean:
	#$(MAKE) -C libmseed/ clean
	rm -rf $(OBJS) $(EXEC)
