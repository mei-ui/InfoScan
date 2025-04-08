/* API for bed format */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
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

void geneSetEnrichmentAnalysis(paraInfo *paraInfo, FILE *outfp, FILE *bgfp, FILE *gofp, FILE *gsfp)
{
  int minGeneNum = 1;
  paraInfo->pvalue = log10Val(paraInfo->pvalue);
  paraInfo->qvalue = log10Val(paraInfo->qvalue);
  geneMap allgeneMap;
  ontologyMap genesetMap;
  ontologyMap ontogeneMap;
  getBackgroundGenes(paraInfo, bgfp, allgeneMap);
  getInputGenes(paraInfo, gsfp, genesetMap);
  getOntologyGenes(paraInfo, gofp, ontogeneMap);

  fprintf(outfp, "queryName\ttermName\tlog10(pval)\tlog10(qval)\t");
  fprintf(outfp, "totalGeneNum\ttermGeneNum\tinputGeneNum\tcommonGeneNum\n");

  ontologyMap::iterator git;
  ontologyMap::iterator oit;
  int bgGeneNum = allgeneMap.size();
  for (git = genesetMap.begin(); git != genesetMap.end(); git++)
  {
    gseaVector geneEnrichVector;
    string queryName = git->first;
    geneMap inputMap = git->second;
    int genesetNum = getCommonGeneNum(inputMap, allgeneMap);
    for (oit = ontogeneMap.begin(); oit != ontogeneMap.end(); oit++)
    {
      string termName = oit->first;
      geneMap termMap = oit->second;
      int termGeneNum = termMap.size();
      int commonNum = getCommonGeneNum(inputMap, termMap);
      if (commonNum >= minGeneNum && termGeneNum >= minGeneNum)
      {
        gseaInfo *gsea = (gseaInfo *)safeMalloc(sizeof(gseaInfo));
        double hyperp  = hypergeometric(bgGeneNum, (double) termGeneNum / (double)bgGeneNum, genesetNum, commonNum);
        gsea->name     = strClone(const_cast<char*>(queryName.c_str()));
        gsea->termName = strClone(const_cast<char*>(termName.c_str()));
        gsea->pval     = log10Val(hyperp);
        gsea->qval     = 0;
        gsea->bgGeneNum   = bgGeneNum;
        gsea->termGeneNum = termGeneNum;
        gsea->genesetNum  = genesetNum;
        gsea->commonNum   = commonNum;
        geneEnrichVector.push_back(gsea);
      } // if
    } // for oit
    if (geneEnrichVector.size() > 0) // if loop
    {
      // output GSEA information
      getGseaFdr(geneEnrichVector);
      outputGseaInfo(paraInfo, outfp, geneEnrichVector);
      freeGseaVector(geneEnrichVector);
    }
  } // for git
}

void outputGseaInfo(paraInfo *paraInfo, FILE *outfp, gseaVector &geneEnrichVector)
{
  for (gseaVector::iterator vecItr = geneEnrichVector.begin(); vecItr != geneEnrichVector.end(); vecItr++)
  {
    gseaInfo *gsea  = *vecItr;
    if (gsea->pval < paraInfo->pvalue && gsea->qval < paraInfo->qvalue)
    {
      fprintf(outfp, "%s\t%s\t%.5f\t%.5f\t", gsea->name, gsea->termName, gsea->pval, gsea->qval);
      fprintf(outfp, "%d\t%d\t%d\t%d\n", gsea->bgGeneNum, gsea->termGeneNum, gsea->genesetNum, gsea->commonNum);
    } // if cutoff
  } // for GSEA
}

int cmpGsea(const gseaInfo *x, const gseaInfo *y)
{
  return x->pval < y->pval;
}

double log10Val(double pval)
{
  double returnVal = 0;
  if (pval == 0)
  {
    returnVal = -325.00;
  }
  else {
    returnVal = log(pval) / log(10);
  }
  return returnVal;
}


void getGseaFdr(gseaVector &BHhash)
{
  sort(BHhash.begin(), BHhash.end(), cmpGsea);
  int rank = 1;
  int increaseRank = 1;
  double lastVal = 1000;
  int arrayNum = BHhash.size();
  for (gseaVector::iterator curr = BHhash.begin(); curr != BHhash.end(); curr++) {
    gseaInfo *gsea = *curr;
    double pval = gsea->pval;
    if (lastVal != pval)
    {
      rank = increaseRank;
    }
    increaseRank++;
    lastVal = pval;
    double FDR = pval + log10Val(arrayNum) - log10Val(rank);
    if (FDR > 0) {
      FDR = 0;
    }
    gsea->qval = FDR;
  }
}

int getCommonGeneNum(geneMap &geneset1Map, geneMap &geneset2Map)
{
  int commonGeneNum = 0;
  geneMap::iterator it;
  for (it = geneset1Map.begin(); it != geneset1Map.end(); it++) {
    if (geneset2Map.find(it->first) != geneset2Map.end()) {
      commonGeneNum++;
    }
  }
  return commonGeneNum;
}

void getBackgroundGenes(paraInfo *paraInfo, FILE *bgfp, geneMap &allgeneMap)
{
  char *line      = NULL;
  char delims[]   = "\t";
  int fieldNum    = 0;
  char **goFields = NULL;
  // all gene analysis
  while (line = getLine(bgfp)) {
    if (feof(bgfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    goFields = splitString(line, delims, &fieldNum);
    if (fieldNum < 2) {
      continue;
    }
    string gene(goFields[0]);
    allgeneMap[gene] = atoi(goFields[1]);
    freeWords(goFields, fieldNum);
    safeFree(line);
  } // bgfile while loops
}

void getInputGenes(paraInfo *paraInfo, FILE *gsfp, ontologyMap &genesetMap)
{
  int i = 0;
  int startIdx = 1;
  char *line      = NULL;
  char delims[]   = "\t";
  int fieldNum    = 0;
  char **goFields = NULL;
  // all gene analysis
  while (line = getLine(gsfp)) {
    if (feof(gsfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    goFields = splitString(line, delims, &fieldNum);
    if (fieldNum < 2) {
      continue;
    }
    string genesetName(goFields[0]);
    for (i = startIdx; i < fieldNum; i++)
    {
      string gene(goFields[i]);
      genesetMap[genesetName][gene] = 1;
    }
    string gene(goFields[0]);
    freeWords(goFields, fieldNum);
    safeFree(line);
  } // gsfile while loops
}

void getOntologyGenes(paraInfo *paraInfo, FILE *gofp, ontologyMap &ontogeneMap)
{
  int i = 0;
  int startIdx = 2;
  char *line      = NULL;
  char delims[]   = "\t";
  int fieldNum    = 0;
  char **goFields = NULL;
  // all gene analysis
  while (line = getLine(gofp)) {
    if (feof(gofp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    goFields = splitString(line, delims, &fieldNum);
    if (fieldNum < 3) {
      continue;
    }
    string pathwayName(goFields[0]);
    for (i = startIdx; i < fieldNum; i++)
    {
      string gene(goFields[i]);
      ontogeneMap[pathwayName][gene] = 1;
    }
    freeWords(goFields, fieldNum);
    safeFree(line);
  } // gofile while loops
}

void freeGsea(gseaInfo *gsea)
{
  safeFree(gsea->name);
  safeFree(gsea->termName);
  safeFree(gsea);
}

void freeGseaVector(gseaVector &sList)
{
  for (gseaVector::iterator vecItr = sList.begin(); vecItr != sList.end(); vecItr++)
  {
    gseaInfo *gsea = *vecItr;
    freeGsea(gsea);
  }
  sList.clear();
}

// hypergeometric($refgeneNum,$idNum/$refgeneNum,$totalGeneNum,$geneNum);
long double hypergeometric(int n, long double p, int k, int r) {
  /*
   n = refgene number
   p = gene number in goterm/refgene-number
   k = gene number picked by methods
   r = gene number picked in goterm
   */
  long double q;
  int np;
  int nq;
  int top;
  int i;
  long double logNchooseK;
  long double lfoo;
  long double sum;
  q = (1 - p);
  np = floor((long double) n * p + 0.5);
  nq = floor((long double) n * q + 0.5);
  logNchooseK = lNchooseK(n, k);
  top = k;
  if (np < k) {
    top = np;
  }
  lfoo = lNchooseK(np, top) + lNchooseK(nq, k - top);
//fprintf(stderr, "%d %.5f %d %d\n", n, lfoo, k, r);
  sum = 0;

  for (i = top; i >= r; i--) {
    sum = sum + exp(lfoo - logNchooseK);
    if (i > r) {
      lfoo = lfoo + log((long double) i / (long double) (np - i + 1))
             + log(
               (long double) (nq - k + i)
               / (long double) (k - i + 1));
    }
  }
  return sum;
}

// ln factorial subroutine
long double lFactorial(int number) {
  long double returnValue = 0;
  for (int i = 2; i <= number; i++) {
    returnValue = returnValue + log(i);
  }
  return returnValue;
}

// ln N choose K subroutine
long double lNchooseK(int n, int k) {
  long double answer = 0;
  int i = 0;
  if (k > (n - k)) {
    k = (n - k);
  }
  for (i = n; i > (n - k); i--) {
    answer = answer + log(i);
  }
  answer = answer - lFactorial(k);
  return answer;
}
