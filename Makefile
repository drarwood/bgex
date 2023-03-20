.PHONY: all
all: bgex

bgex: *.cpp
	g++ -O2 -Wall -pthread *.cpp -c
	g++ -O2 -pthread -o bgex *.o sqlite/sqlite3.o -ldl -I./zlib-1.2.13 -lz

.PHONY: install
install:
	cp bgex /usr/local/bin/

.PHONY: clean
clean:
	rm -f bgex
