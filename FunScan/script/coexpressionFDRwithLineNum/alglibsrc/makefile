CC=g++
CXXC=g++
LIBS=-lm
INCLUDE=-I.
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE

O = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o

alglib.a: $(O) 
	ar rcus alglib.a $(O)
	
alglibinternal.o: alglibinternal.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c alglibinternal.cpp

alglibmisc.o: alglibmisc.cpp alglibmisc.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c alglibmisc.cpp

ap.o: ap.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c ap.cpp

dataanalysis.o: dataanalysis.cpp dataanalysis.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c dataanalysis.cpp

diffequations.o: diffequations.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c diffequations.cpp

fasttransforms.o: fasttransforms.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c fasttransforms.cpp

integration.o: integration.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c integration.cpp

interpolation.o: interpolation.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c interpolation.cpp

linalg.o: linalg.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c linalg.cpp
	
optimization.o: optimization.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c optimization.cpp

solvers.o: solvers.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c solvers.cpp
	
specialfunctions.o: specialfunctions.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c specialfunctions.cpp

statistics.o: statistics.cpp
	$(CC) $(CFLAGS) $(HG_DEFS) -c statistics.cpp

clean:
	rm -f *.o
