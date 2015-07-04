#This should bild everything
CC = gcc
CFLAGS = -O3 -flto -Wall -fmessage-length=0  -std=gnu11
LIBS = -lm

all: green.o 
	$(CC) $(CFLAGS) green.o -o green.run $(LIBS)
green: green.c
	$(CC) $(CFLAGS) -c green.c
clean:
	rm -f *.o *.run