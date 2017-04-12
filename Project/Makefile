main: main_linux.o cfd.o
	gcc -o main bin/main_linux.o bin/cfd.o

main_linux.o:
	gcc -o bin/main_linux.o -c src/main_linux.c -W -Wall -lm -std=gnu99

cfd.o:
	gcc -o bin/cfd.o -c src/cfd.c -W -Wall -lm -std=gnu99
