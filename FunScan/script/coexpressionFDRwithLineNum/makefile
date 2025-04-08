CC=g++
CXXC=g++
LIBS=-lm
INCLUDE=-I.
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
alglibsrc = alglibsrc


coExpressionFDR: coExpressionFDR.o snoSeeker_utils.o
	$(CC) $(CFLAGS) $(HG_DEFS) $(INCLUDE) -o coExpressionFDR coExpressionFDR.o snoSeeker_utils.o -I. -I./alglibsrc/ ./alglibsrc/alglib.a $(LIBS)


coExpressionFDR.o: coExpressionFDR.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c coExpressionFDR.cpp

	
snoSeeker_utils.o: snoSeeker_utils.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c snoSeeker_utils.cpp
			
clean:
	cd $(alglibsrc); make clean -f makefile
	rm -f *.o
