import pyranges as pr
import pandas as pd
import csv
import argparse
import re
from pyranges import PyRanges


parser = argparse.ArgumentParser(description='Test for argparse')
parser.add_argument('--gtffile', '-f', help='input gtf file',required=True)
parser.add_argument('--txtfile', '-i', help='inpu txt file',required=True)
parser.add_argument('--outfile', '-o', help='output gtf file',default="./output.gtf")
args = parser.parse_args()

def to_rows(anno):
    rowdicts = []
    try:
        l = anno.head(1)
        for l in l:
            l.replace('"', '').replace(";", "").split()
    except AttributeError:
        raise Exception("Invalid attribute string: {l}. If the file is in GFF3 format, use pr.read_gff3 instead.".format(l=l))
    for l in anno:
        rowdicts.append({k: v
                         for k, v in [kv.replace('""', '"NA"').replace('"', '').split(None, 1)
                                      for kv in l[:-1].split("; ")]})
    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def to_rows_keep_duplicates(anno):
    rowdicts = []
    for l in anno:
        rowdict = {}
        # l[:-1] removes final ";" cheaply
        for k, v in (kv.replace('"', '').split(None, 1) for kv in l[:-1].split("; ")):
            if k not in rowdict:
                rowdict[k] = v
            elif k in rowdict and isinstance(rowdict[k], list):
                rowdict[k].append(v)
            else:
                rowdict[k] = [rowdict[k], v]
        rowdicts.append({
            k: ','.join(v) if isinstance(v, list) else v
            for k, v in rowdict.items()
        })
    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)

def read_gtf_full(f, as_df=False, nrows=None, skiprows=0,duplicate_attr=False):
    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }
    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
    df_iter = pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=names,
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows,
        nrows=nrows,comment="#")
    _to_rows = to_rows_keep_duplicates if duplicate_attr else to_rows
    dfs = []
    for df in df_iter:
        extra = _to_rows(df.Attribute)
        df = df.drop("Attribute", axis=1)
        ndf = pd.concat([df, extra], axis=1, sort=False)
        dfs.append(ndf)
    df = pd.concat(dfs, sort=False)
    df.loc[:, "Start"] = df.Start
    if not as_df:
        return PyRanges(df)
    else:
        return df

def writeGTF(inGTF,file_path):
    """
    Write a GTF dataframe into a file
    :param inGTF: GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
    :param file_path: path/to/the/file.gtf
    :returns: nothing
    """
    cols=inGTF.columns.tolist()
    if len(cols) == 9:
        if 'attribute' in cols:
            df=inGTF
    else:
        df=inGTF[cols[:8]]
        df['attribute']=""
        for c in cols[8:]:
            if c == cols[len(cols)-1]:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'";'
            else:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'"; '
    df.to_csv(file_path, sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)

def main(gtffile,txtfile,outfile):
    #read gtf
    gtf_file = read_gtf_full(gtffile,as_df=True)
    # as DataFrame
    df = gtf_file
    #read txt 
    txt_file = pd.read_table(txtfile,header=None)
    rename=["id"]
    txt_file.columns=rename
    #inner_join
    gtf_file_2 = pd.merge(df,txt_file,how='inner',left_on="transcript_id",right_on="id",copy=False)
    gtf_file_2.drop('id',axis=1, inplace=True)
    writeGTF(gtf_file_2,outfile)

if __name__ == '__main__':
    try:
        main(args.gtffile, args.txtfile, args.outfile)
    except Exception as e:
        print(e)

