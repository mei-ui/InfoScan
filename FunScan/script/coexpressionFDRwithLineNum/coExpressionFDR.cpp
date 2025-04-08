#include <math.h>
#include "./alglibsrc/stdafx.h"
#include "./alglibsrc/alglibmisc.h"
#include "./alglibsrc/alglibinternal.h"
#include "./alglibsrc/linalg.h"
#include "./alglibsrc/statistics.h"
#include "./alglibsrc/dataanalysis.h"
#include "./alglibsrc/specialfunctions.h"
#include "./alglibsrc/solvers.h"
#include "./alglibsrc/optimization.h"
#include "./alglibsrc/diffequations.h"
#include "./alglibsrc/fasttransforms.h"
#include "./alglibsrc/integration.h"
#include "./alglibsrc/interpolation.h"
#include "snoSeeker_utils.h"
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
using namespace alglib;

void usage(void);

template<typename T>
string NumberToString(T Number) {
	ostringstream ss;
	ss << Number;
	return ss.str();
}
typedef pair<string, double> PAIR;
int cmp(const PAIR& x, const PAIR& y) {
	return x.second < y.second;
}

char version[] = "coExpression version 0.1";

int main(int argc, char **argv) {
	int n = 9;
	int i = 0;
	int j = 0;
	int k = 0;
	int s = 0;
	FILE *fp;
	FILE *outfp;
	char *strLine;
	char *tmpLine;
	int fieldNum = 0;
	char **fields = NULL;
	char *outfile = NULL;
	int    method = 1;
	int showVersion = 0;
	int showHelp = 0;
	int genesetNum = 1;
	int geneNum = 0;
	int sampleNum = -1;
	double correctedPval = 0;
	double pValue = 0.05;
	double FDRcutoff = 0.05;
	double relationVal = 0.00000001;
	real_2d_array expArray;
	char **geneNames = NULL;
	char **geneSymbols = NULL;
	if (argc == 1) {
		usage();
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 'o':
				outfile = argv[++i];
				break;
			case 'm':
				method = atoi(argv[++i]);
				break;
			case 'p':
				pValue = atof(argv[++i]);
				break;
			case 'q':
				FDRcutoff = atof(argv[++i]);
				break;
			case 'n':
				genesetNum = atoi(argv[++i]);
				break;
			case 'r':
				relationVal = atof(argv[++i]);
				break;
			case 'v':
				showVersion = 1;
				break;
			case 'h':
				showHelp = 1;
				break;
			default:
				usage();
			}
		} else { /* doesn't start with '-' should be filename */
			if (i != argc - 1)
				usage();
			fp = (FILE *) fopen(argv[i], "r");
			if (fp == NULL) {
				fprintf(stderr, "ERROR: Can't open %s\n", argv[i]);
				usage();
			}
			break;
		}
	}

	if (outfile == NULL) {
		outfp = stdout;
	} else {
		outfp = (FILE *) fopen(outfile, "w");
		if (outfp == NULL) {
			fprintf(stderr, "ERROR: Can't open %s\n", outfile);
			usage();
		}
	}

	if (showVersion) {
		fprintf(stderr, "%s", version);
		exit(1);
	}

	if (showHelp) {
		usage();
		exit(1);
	}

	i = 0;
	while (strLine = getLine(fp)) {  //第一次读入文件，记录最小sample number及geneNum
		tmpLine = strLine;
		tmpLine = skipStartWhitespace(tmpLine);
		fieldNum = 0;
		fields = NULL;
		if (feof(fp) || tmpLine == NULL) {
			safeFree(strLine);
			break;
		}
		fields = splitWhitespace(tmpLine, &fieldNum);
		if (fieldNum >= 3) { // minimum number of expression values
			if (geneNum > 1 && sampleNum != fieldNum - 2) {
				fprintf(stderr, "sample number is not constant in %s\n",
						tmpLine);
				exit(1);
			}
			sampleNum = fieldNum - 2; //计算sampleNum
			geneNum += 1;             //计算geneNum
		}
		freeWords(fields, fieldNum);
		safeFree(strLine);
	}

	geneNames = (char **) safeMalloc(sizeof(char *) * (geneNum));  //初始化geneNames数组
	geneSymbols = (char **) safeMalloc(sizeof(char *) * (geneNum)); //初始化geneSymbols数组

	// set the 2 dimension array to save the expressed values
	expArray.setlength(geneNum, sampleNum);   //初始化expArray

	i = 0;
	j = 0;
	k = 0;
	s = 1;
	fseek(fp, 0, SEEK_SET); // Reposition stream position indicator
	while (strLine = getLine(fp)) { //第二次读取文件，给geneNames，geneSymbols，expArray赋值
		tmpLine = strLine;
		tmpLine = skipStartWhitespace(tmpLine);
		fieldNum = 0;
		fields = NULL;
		if (feof(fp) || tmpLine == NULL) {
			safeFree(strLine);
			break;
		}
		fields = splitWhitespace(tmpLine, &fieldNum);
		if (fieldNum >= 3) { // minimum number of expression values
			geneNames[i] = strClone(fields[0]);
			geneSymbols[i] = strClone(fields[1]);
			for (j = 2; j < fieldNum; j++) { //reading expression matrix into variable expArray
				expArray[i][j - 2] = atof(fields[j]);
			}
			i++;
		}
		freeWords(fields, fieldNum);
		safeFree(strLine);
	}

	fprintf(outfp,
			"geneOrder\tqueryGene\tquerySymbol\tcoexpGene\tcoexpSymbol\tsampleNum\tcoe\tpValue\tGeneNum\trank\tBFpvalue\tFDR\n"); // output the test
	// start coexpression analysis
		real_1d_array x;
		x.setlength(sampleNum);
		for (k = 0; k < sampleNum; k++) {
			x[k] = expArray[genesetNum][k];
		}
		vector<PAIR> BHhash;

//	for (i = 0; i < genesetNum; i++) {  //对于每一个待计算的基因（文件中的前n行），循环计算其他行与之的相关性
//		real_1d_array x;
//		x.setlength(sampleNum);
//		for (k = 0; k < sampleNum; k++) {
//			x[k] = expArray[i][k];
//		}
//		vector<PAIR> BHhash;
		
		for (j = 0; j < geneNum; j++) {

			if (strcmp(geneNames[genesetNum], geneNames[j]) == 0
					|| strcmp(geneSymbols[genesetNum], geneSymbols[j]) == 0) {
				continue;
			}
			real_1d_array y;
			y.setlength(sampleNum);
			for (k = 0; k < sampleNum; k++) {
				y[k] = expArray[j][k];
			}
			double bothtails;
			double lefttail;
			double righttail;
			double r;
			if (method == 1) {
			r = pearsoncorr2(x, y); // pearson correlation
			//fprintf(stderr, "%s %s %.5f\n", geneNames[i], geneNames[j], x[0]);
			pearsoncorrelationsignificance(r, sampleNum, bothtails, lefttail,
					righttail); // pearson correlation test
			}
			else if (method == 2) {
			r = spearmancorr2(x, y); // sperarman correlation
			spearmanrankcorrelationsignificance(r, sampleNum, bothtails, lefttail,
					righttail); // sperarman correlation test
			}
			else {
				fprintf(stderr, "error: please select correlation method use 1 or 2\n");
				usage();
			}
			if (bothtails<0) {
				fprintf(stderr, "error: correlation test r:%g\tpvalue: %g\tleft:%g\tright:%g\n", r, bothtails, lefttail, righttail);
				bothtails = 0;
			}
			string queryGeneName = geneNames[genesetNum];
			string queryGeneSymbol = geneSymbols[genesetNum];
			string targetGeneName = geneNames[j];
			string targetGeneSymbol = geneSymbols[j];
			string key = NumberToString(j+1) + "\t" + queryGeneName + "\t" + queryGeneSymbol + "\t"
					+ targetGeneName + "\t" + targetGeneSymbol + "\t"
					+ NumberToString(sampleNum) + "\t" + NumberToString(r)
					+ "\t" + NumberToString(bothtails);
			BHhash.push_back(make_pair(key, bothtails));
			//if (correctedPval<pValue) {
			//fprintf(outfp, "%d\t%s\t%s\t%s\t%s\t%d\t%.9f\t%g\n", s, geneNames[i], geneSymbols[i], geneNames[j], geneSymbols[j], sampleNum, r, bothtails); // output the test
			//s++;
			//fflush(outfp);			
			//}
			//fprintf(outfp, "%.15f,", bothtails);
		} /// for geneset Number end 
		  ///////////// FDR ///////////////
		  // FDR line
		sort(BHhash.begin(), BHhash.end(), cmp);
		int rank = 1;
		int increaseRank = 1;
		double lastVal = -100;
		for (vector<PAIR>::iterator curr = BHhash.begin(); curr != BHhash.end();
				curr++) {
			string key((*curr).first);
			double pval = (*curr).second;
			if (lastVal != pval) {
				rank = increaseRank;
			}
			increaseRank++;
			lastVal = pval;
			double FDR = pval * BHhash.size() / (double) rank;
			double correctedPvalue = pval * BHhash.size();
			if (correctedPvalue > 1) {
				correctedPvalue = 1.0;
			}
			if (FDR > 1) {
				FDR = 1.0;
			}
			if (FDR < FDRcutoff &&  pval<pValue) {
				fprintf(outfp, "%s\t%d\t%d\t%g\t%g\n", key.c_str(), BHhash.size(),
						rank, correctedPvalue, FDR);
				fflush(outfp);
			}
		}
		// fdr end
		//////////////FDR END/////////////
//	} // analysis end
	  //fprintf(outfp, "finished...\n");
	fflush(outfp);
	freeWords(geneNames, geneNum);
	freeWords(geneSymbols, geneNum);
	fclose(fp);
	fclose(outfp);
	return 0;
}

void usage(void) {

	fprintf(stderr, "%s",
			"Usage:  coExpression [options] <expression file>\n\
	[options]\n\
	-o<out file>    : output file\n\
	-p<double>      : p value cutoff[default<0.05]\n\
	-q<double>      : q(false discovery rate, FDR) value cutoff[default<0.05]\n\
	-n<int>         : gene set number(>=1)[default=1]\n\
	-r<double>      : absolute (relation Val)[>0.2]\n\
	-m<int>         : methods[pearson=1, spearman=2] default=1 .pearson\n\
  -v<version>     : coExpression version\n\
	-h<help>        : help informations\n\
	");

	exit(1);
}
