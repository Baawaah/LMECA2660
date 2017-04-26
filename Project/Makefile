main: main_linux.o cfd_main.o cfd_numerical.o cfd_test.o cfd_integrator.o
	gcc -o main bin/main_linux.o bin/cfd_main.o bin/cfd_numerical.o bin/cfd_integrator.o bin/cfd_test.o -W -Wall -lm -std=gnu99

main_linux.o:
	gcc -o bin/main_linux.o -c src/main_linux.c -W -Wall -lm -std=gnu99

cfd_main.o:
	gcc -o bin/cfd_main.o -c src/cfd_main.c -W -Wall -lm -std=gnu99

cfd_numerical.o:
	gcc -o bin/cfd_numerical.o -c src/cfd_numerical.c -W -Wall -lm -std=gnu99

cfd_integrator.o:
	gcc -o bin/cfd_integrator.o -c src/cfd_integrator.c -W -Wall -lm -std=gnu99

cfd_test.o:
	gcc -o bin/cfd_test.o -c src/cfd_test.c -W -Wall -lm -std=gnu99
