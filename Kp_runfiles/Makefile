CC	= g++
LD	= g++

CCFLAGS	= -O0 -g -c `root-config --cflags`
LDFLAGS	= `root-config --libs`

fitkpalltheorysircentallwgtbgfixnewalppcs: fitkpalltheorysircentallwgtbgfixnewalppcs.o
	$(CC) $^ -o $@ $(LDFLAGS)

fitkpalltheorysircentallwgtbgfixnewalppcs.o: fitkpalltheorysircentallwgtbgfixnewalppcs.cxx
	$(CC) $^ -o $@ $(CCFLAGS)
