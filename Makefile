DEBUG = 0
DUMPDATA = 1

CC = gcc
EXEC = ms2fft
COMMON = -I./libmseed/ -I.
CFLAGS =  -Wall
LDFLAGS = -L./libmseed -Wl,-rpath,./libmseed
LDLIBS = -Wl,-Bstatic -lmseed -Wl,-Bdynamic -lm

OBJS = main.o

ifeq ($(DEBUG), 1)
CFLAGS += -O0 -g -DDEBUG=1
endif

ifeq ($(DUMPDATA), 1)
CFLAGS += -DDUMPDATA=1
endif

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(MAKE) -C libmseed/ static
	$(CC) $(COMMON) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(COMMON) $(CFLAGS) -c $< -o $@

clean:
	$(MAKE) -C libmseed/ clean
	rm -rf $(OBJS) $(EXEC)