import pyranges as pr
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Test for argparse')
parser.add_argument("-a","--protein_coding", help='input protein_coding gtf file',required=True)
parser.add_argument("-b","--annotated_lncRNA", help='input annotated_lncRNA gtf file',required=True)
parser.add_argument("-c","--unannotated_lncRNA", help='input unannotated_lncRNA gtf file',required=True)
parser.add_argument("-e","--ercc", help='input ERCC gtf file',required=True)
parser.add_argument("-f","--pseudo_gene", help='input pseudogene gtf file',required=True)
parser.add_argument("-g","--other_rna", help='input other_rnas gtf file',required=True)
parser.add_argument("-o","--output", help='output')
parser.add_argument("-d","--unannotated_protein", help='input unannotated_protein gtf file',required=True)
args = parser.parse_args()

def merge_gtf(protein_coding,annotated_lncRNA,pseudo_gene,other_rna,unannotated_lncRNA,unannotated_protein,ercc,output):
	protein_coding=pr.read_gtf(protein_coding)
	protein_coding=protein_coding.df

	annotated_lncRNA=pr.read_gtf(annotated_lncRNA)
	annotated_lncRNA=annotated_lncRNA.df

	pseudo_gene=pr.read_gtf(pseudo_gene)
	pseudo_gene=pseudo_gene.df

	other_rna=pr.read_gtf(other_rna)
	other_rna=other_rna.df

	unannotated_lncRNA=pr.read_gtf(unannotated_lncRNA)
	unannotated_lncRNA=unannotated_lncRNA.df
	unannotated_lncRNA['gene_name']=unannotated_lncRNA['gene_id']
	unannotated_lncRNA['transcript_id'] = unannotated_lncRNA['transcript_id'].replace(";","",regex=True)

	unannotated_protein=pr.read_gtf(unannotated_protein)
	unannotated_protein=unannotated_protein.df
	unannotated_protein['gene_name']=unannotated_protein['gene_id']
	unannotated_protein['transcript_id'] = unannotated_protein['transcript_id'].replace(";","",regex=True)

	ercc=pr.read_gtf(ercc)
	ercc=ercc.df
	ercc['gene_name']=ercc['gene_id']

	protein_coding_1=protein_coding['gene_name']
	protein_coding_1=pd.DataFrame(protein_coding_1)
	protein_coding_1=np.unique(protein_coding_1)
	protein_coding_1=pd.DataFrame(protein_coding_1)
	protein_coding_1.columns=['gene_name']
	protein_coding_1['type']='protein_coding'

	annotated_lncRNA_1=annotated_lncRNA['gene_name']
	annotated_lncRNA_1=pd.DataFrame(annotated_lncRNA_1)
	annotated_lncRNA_1=np.unique(annotated_lncRNA_1)
	annotated_lncRNA_1=pd.DataFrame(annotated_lncRNA_1)
	annotated_lncRNA_1.columns=['gene_name']
	annotated_lncRNA_1['type']='annotated_lncRNA'

	pseudo_gene_1=pseudo_gene['gene_name']
	pseudo_gene_1=pd.DataFrame(pseudo_gene_1)
	pseudo_gene_1=np.unique(pseudo_gene_1)
	pseudo_gene_1=pd.DataFrame(pseudo_gene_1)
	pseudo_gene_1.columns=['gene_name']
	pseudo_gene_1['type']='pseudo_gene'

	other_rna_1=other_rna['gene_name']
	other_rna_1=pd.DataFrame(other_rna_1)
	other_rna_1=np.unique(other_rna_1)
	other_rna_1=pd.DataFrame(other_rna_1)
	other_rna_1.columns=['gene_name']
	other_rna_1['type']='other_rna'

	unannotated_lncRNA_1=unannotated_lncRNA['transcript_id']
	unannotated_lncRNA_1=pd.DataFrame(unannotated_lncRNA_1)
	unannotated_lncRNA_1=np.unique(unannotated_lncRNA_1)
	unannotated_lncRNA_1=pd.DataFrame(unannotated_lncRNA_1)
	unannotated_lncRNA_1.columns=['gene_name']
	unannotated_lncRNA_1['type']='unannotated_lncRNA'
	
	unannotated_protein_1=unannotated_protein['transcript_id']
	unannotated_protein_1=pd.DataFrame(unannotated_protein_1)
	unannotated_protein_1=np.unique(unannotated_protein_1)
	unannotated_protein_1=pd.DataFrame(unannotated_protein_1)
	unannotated_protein_1.columns=['gene_name']
	unannotated_protein_1['type']='unannotated_coding_transcript'

	un_pcl_1=ercc['gene_id']
	un_pcl_1=pd.DataFrame(un_pcl_1)
	un_pcl_1=np.unique(un_pcl_1)
	un_pcl_1=pd.DataFrame(un_pcl_1)
	un_pcl_1.columns=['gene_name']
	un_pcl_1['type']='ERCC'

	all_txt=pd.concat([protein_coding_1,annotated_lncRNA_1,unannotated_lncRNA_1,unannotated_protein_1,un_pcl_1,pseudo_gene_1,other_rna_1])
	all_txt.to_csv(output, sep='\t', index=False,header=False)

if __name__ == '__main__':
	try:
		merge_gtf(args.protein_coding,args.annotated_lncRNA,args.pseudo_gene,args.other_rna,args.unannotated_lncRNA,args.unannotated_protein,args.ercc,args.output)
	except Exception as e:
		print(e)