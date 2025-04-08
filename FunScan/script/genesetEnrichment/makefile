CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type

geneSetEnrichment: geneSetEnrichment.o geneSetEnrichmentMain.o bioUtils.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} -o geneSetEnrichment geneSetEnrichmentMain.o geneSetEnrichment.o bioUtils.o $(LIBS) 

geneSetEnrichment.o: geneSetEnrichment.cpp geneSetEnrichment.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} -c geneSetEnrichment.cpp
	
geneSetEnrichmentMain.o: geneSetEnrichmentMain.cpp geneSetEnrichment.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} -c geneSetEnrichmentMain.cpp

bioUtils.o: bioUtils.cpp bioUtils.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} -c bioUtils.cpp

clean:
	rm -f *.o
