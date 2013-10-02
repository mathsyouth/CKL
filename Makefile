CC = g++

CFLAGS		= -O2 -I.

.SUFFIXES: .C .cpp

OBJECTS		= lib/CKL.o \
		  lib/BV.o \
		  lib/Build.o \
		  lib/TriDist.o

CLEAN		= $(OBJECTS) lib/libCKL.a include/*.h

library: $(OBJECTS)
	/bin/rm -f lib/libCKL.a
	ar ruv lib/libCKL.a $(OBJECTS)
	cp src/CKL.h include/
	cp src/CKL_Compile.h include/
	cp src/CKL_Internal.h include/
	cp src/BV.h include/
	cp src/Tri.h include/
	cp src/MatVec.h include/

lib/BV.o: src/BV.cpp
	$(CC) $(CFLAGS) -c src/BV.cpp -o lib/BV.o
lib/CKL.o: src/CKL.cpp
	$(CC) $(CFLAGS) -c src/CKL.cpp -o lib/CKL.o
lib/Build.o: src/Build.cpp
	$(CC) $(CFLAGS) -c src/Build.cpp -o lib/Build.o
lib/TriDist.o: src/TriDist.cpp
	$(CC) $(CFLAGS) -c src/TriDist.cpp -o lib/TriDist.o

clean:
	/bin/rm -f $(CLEAN)
