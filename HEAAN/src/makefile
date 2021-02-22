all: clean main

main:
	g++ main.cpp -o main -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread

sanitize:
	g++ main.cpp -o main -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -g -I/../src -pthread -fsanitize=address,undefined

clean:
	rm -rf main
