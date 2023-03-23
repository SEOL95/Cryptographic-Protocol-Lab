INCDIR   := inc/miracl
SRCDIR   := src

CC=gcc

SOURCE = main.c core.a

HEADER = -I./$(INCDIR) -I./$(INCDIR)/miracl
CFLAGS  = -march=native -mtune=native -O3 -fomit-frame-pointer

all:
   $(CC) $(SOURCE) $(HEADER) -o main