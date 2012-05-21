CC=gcc
CFLAGS=-O1 -g

all:	ushuffle fasta_ushuffle

ushuffle:	ushuffle.o	main.o

fasta_ushuffle:	ushuffle.o	fasta_ushuffle.o

clean:
	rm -f *.o ushuffle fasta_ushuffle
