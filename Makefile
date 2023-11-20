.PHONY: all
all: bgex

ZSTD_DIR = ./zstd-1.5.2/lib

bgex: *.cpp
	gcc -Wall -pthread sqlite/sqlite3.c -c -o sqlite/sqlite3.o
	g++ -O2 -Wall -pthread *.cpp -c
	g++ -O2 -pthread -o bgex *.o sqlite/sqlite3.o -ldl -I./zlib-1.2.13 -I./zstd -lz -L./zstd/lib/ -lzstd

.PHONY: install
install:
	cp bgex /usr/local/bin/

.PHONY: clean
clean:
	rm -f bgex
