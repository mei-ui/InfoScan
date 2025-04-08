#ifndef geneSetEnrichment_HEAD_H
#define geneSetEnrichment_HEAD_H

typedef struct parameterInfo
{
  int    verbose;
  double pvalue;
  double qvalue;
} paraInfo;

typedef struct gseaInfo
{
  char*  name;
  char*  termName;
  double pval;
  double qval;
  int    bgGeneNum;
  int    termGeneNum;
  int    genesetNum;
  int    commonNum;
} gseaInfo;

typedef vector<gseaInfo *> gseaVector;

typedef map<string, int> geneMap;

typedef map<string, geneMap> ontologyMap;

void geneSetEnrichmentAnalysis(paraInfo *paraInfo, FILE *outfp, FILE *bgfp, FILE *gofp, FILE *gsfp);

void outputGseaInfo(paraInfo *paraInfo, FILE *outfp, gseaVector &geneEnrichVector);

int cmpGsea(const gseaInfo *x, const gseaInfo *y);

double log10Val(double pval);

void getGseaFdr(gseaVector &BHhash);

int getCommonGeneNum(geneMap &geneset1Map, geneMap &geneset2Map);

void getBackgroundGenes(paraInfo *paraInfo, FILE *bgfp, geneMap &allgeneMap);

void getInputGenes(paraInfo *paraInfo, FILE *gsfp, ontologyMap &genesetMap);

void getOntologyGenes(paraInfo *paraInfo, FILE *gofp, ontologyMap &ontogeneMap);

void freeGsea(gseaInfo *gsea);

void freeGseaVector(gseaVector &sList);

long double hypergeometric(int n, long double p, int k, int r);

long double lFactorial(int number);

long double lNchooseK(int n, int k);

#endif /* End geneSetEnrichment_HEAD_H */
