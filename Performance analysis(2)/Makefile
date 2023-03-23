INCDIR   := inc 
SRCDIR   := src

CC=gcc

HEADER = -I./inc -I./inc/miracl 
CFLAGS= -O3

ALL: #youngkyung
	$(CC) $(CFLAGS) cpucycles.c main.c core.a -o main $(HEADER)
		
clean:
	rm main