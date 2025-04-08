/*
geneSetEnrichment $2014/12/09/$ @Jian-Hua Yang yangjh7@mail.sysu.edu.cn
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <getopt.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"

#include "geneSetEnrichment.h"

char version[] = "geneSetEnrichment version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *bgFile   = NULL; // back-ground file
  char *gsFile   = NULL; // gene set file
  char *goFile   = NULL; // gene ontology
  char *outFile  = NULL; // output file
  FILE *outfp    = NULL;
  FILE *bgfp     = NULL;
  FILE *gsfp     = NULL;
  FILE *gofp     = NULL;
  int showVersion = 0;
  int showHelp    = 0;
  int c           = 0;
  int i           = 0;
  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVo:b:g:s:p:q:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "output" , required_argument, NULL, 'o' },
    { "back-ground" , required_argument , NULL, 'b' },
    { "gene-ontology" , required_argument, NULL, 'g' },
    { "gene-set" , required_argument , NULL, 's' },
    { "pvalue" , required_argument, NULL, 'p' },
    { "qvalue" , required_argument, NULL, 'q' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose    = 0;
  paraInfo.pvalue     = 0.05;
  paraInfo.qvalue     = 0.05;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'o':
      outFile  = optarg;
      break;
    case 'b':
      bgFile  = optarg;
      break;
    case 'g':
      goFile = optarg;
      break;
    case 's':
      gsFile  = optarg;
      break;
    case 'p':
      paraInfo.pvalue = atof(optarg);
      break;
    case 'q':
      paraInfo.qvalue = atof(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    } // switch end
  }// while end

  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  if (bgFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: -b <background gene file>\n");
    usage();
  }
  else {
    bgfp = (FILE *) fopen(bgFile, "r");
    if (bgfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open background gene file: %s\n", bgFile);
      usage();
    }
  }

  if (goFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: -g <gene ontology/pathway file>\n");
    usage();
  }
  else {
    gofp = (FILE *) fopen(goFile, "r");
    if (gofp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open background gene file: %s\n", goFile);
      usage();
    }
  }

  if (gsFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: -s <gene set file>\n");
    usage();
  }
  else {
    gsfp = (FILE *) fopen(gsFile, "r");
    if (gsfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open background gene file: %s\n", gsFile);
      usage();
    }
  }

  if (outFile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outFile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outFile);
      usage();
    }
  }

  geneSetEnrichmentAnalysis(&paraInfo, outfp, bgfp, gofp, gsfp);

  fclose(bgfp);
  fclose(gofp);
  fclose(gsfp);
  fclose(outfp);
  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  geneSetEnrichment [options] -b <background genes> -g <genes of ontology/pathway> -s <gene set for the enrichment>\n\
geneSetEnrichment: for Gene set enrichment analysis\n\
[options]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : geneSetEnrichment version\n\
-h/--help                      : help informations\n\
-o/--outfile <string>          : output file\n\
-b/--back-ground <string>      : background gene file\n\
-g/--gene-ontology <string>    : genes of ontology/pathway\n\
-s/--gene-set <string>         : gene set for the enrichment\n\
-p/--pvalue <double>           : p-value for the enrichment[p<0.05]\n\
-q/--pvalue <double>           : q-value/FDR for the enrichment[q<0.05]\n\
");
  exit(1);
}
