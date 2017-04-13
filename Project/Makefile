main: main_linux.o cfd_main.o cfd_numerical
	gcc -o main bin/main_linux.o bin/cfd_main.o bin/cfd_numerical.o

main_linux.o:
	gcc -o bin/main_linux.o -c src/main_linux.c -W -Wall -lm -std=gnu99

cfd_main.o:
	gcc -o bin/cfd_main.o -c src/cfd_main.c -W -Wall -lm -std=gnu99

cfd_numerical.o:
	gcc -o bin/cfd_numerical.o -c src/cfd_numerical.c -W -Wall -lm -std=gnu99
