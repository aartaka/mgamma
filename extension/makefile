.PHONY: all clean

CC = gcc
CFLAGS += -g $(shell pkg-config --cflags openblas guile-3.0)
LDFLAGS += $(shell pkg-config --libs openblas guile-3.0)

libmgamma.so: extension.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -shared -fPIC -o $@

all: libmgamma.so

clean:
	-rm libmgamma.so
