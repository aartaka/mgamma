.PHONY: all indent clean

CC = gcc
CFLAGS += -g $(shell pkg-config --cflags openblas guile-3.0)
LDFLAGS += $(shell pkg-config --libs openblas guile-3.0)

libmgamma.so: extension.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -shared -fPIC -o $@

all: libmgamma.so

indent: extension.c
	indent -linux --procnames-start-lines $^

clean:
	-rm libmgamma.so
