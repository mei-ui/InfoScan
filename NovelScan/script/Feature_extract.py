import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
from multiprocessing import Pool, cpu_count

def find_orf(sequence):
    sequence = str(sequence).upper()
    start_codons = [m.start() for m in re.finditer("ATG", sequence)]
    if not start_codons:
        return 0, 0, None
    
    stop_codons = sorted([m.start() for stop in ["TAA", "TAG", "TGA"] for m in re.finditer(stop, sequence)])
    seq_length = len(sequence)
    
    orf_lengths = []
    orf_starts = []
    orf_stops = []
    
    for start in start_codons:
        for stop in stop_codons:
            if stop > start and (stop - start) % 3 == 0:
                orf_lengths.append(stop - start + 3)
                orf_starts.append(start)
                orf_stops.append(stop)
                break
    
    if not orf_lengths:
        max_len = 0
        max_cov = 0
        orf_seq = None
    else:
        max_idx = np.argmax(orf_lengths)
        max_len = orf_lengths[max_idx]
        max_cov = max_len / seq_length
        orf_seq = sequence[orf_starts[max_idx]: orf_stops[max_idx] + 3]
    
    return max_len, max_cov, orf_seq

def count_base(sequence, base):
    return sequence.count(base)

def look_up_content_probability(value, base, content_prob, content_weight, content_para):
    if value < 0:
        return 0
    for i, threshold in enumerate(content_para):
        if value >= threshold:
            return content_prob[base][i] * content_weight[base]
    return 0

def look_up_position_probability(value, base, position_prob, position_weight, position_para):
    if value < 0:
        return 0
    for i, threshold in enumerate(position_para):
        if value >= threshold:
            return position_prob[base][i] * position_weight[base]
    return 0

def slice_sequence(sequence, start, step):
    return sequence[start::step]

def fickett_value(sequence):
    sequence = sequence.upper()
    if len(sequence) < 2:
        return 0
    
    total_base = len(sequence)
    base_contents = {base: count_base(sequence, base) / total_base for base in "ACGT"}
    
    phase_0, phase_1, phase_2 = slice_sequence(sequence, 0, 3), slice_sequence(sequence, 1, 3), slice_sequence(sequence, 2, 3)
    base_positions = {
        base: max(count_base(phase_0, base), count_base(phase_1, base), count_base(phase_2, base)) /
              (min(count_base(phase_0, base), count_base(phase_1, base), count_base(phase_2, base)) + 1.0)
        for base in "ACGT"
    }
    
    fickett_score = sum(
        look_up_content_probability(base_contents[base], base, content_prob, content_weight, content_para)
        + look_up_position_probability(base_positions[base], base, position_prob, position_weight, position_para)
        for base in "ACGT"
    )
    
    return fickett_score

# 数据映射
position_prob = {
    "A": [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
    "C": [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
    "G": [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
    "T": [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]
}
position_weight = {"A": 0.26, "C": 0.18, "G": 0.31, "T": 0.33}
position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]

content_prob = {
    "A": [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
    "C": [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
    "G": [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
    "T": [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]
}
content_weight = {"A": 0.11, "C": 0.12, "G": 0.15, "T": 0.14}
content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]


def calculate_physicochemical(sequence):
    sequence = str(sequence).upper().replace("N", "")
    if len(sequence) < 2:
        return 0, 0, 0, 0

    protein_seq = str(Seq(sequence).translate()).replace("*", "")  # 移除终止密码子

    if len(protein_seq) < 1:
        return 0, 0, 0, 0  # 避免空序列报错

    analysis = ProteinAnalysis(protein_seq)

    return (
        len(protein_seq),
        analysis.molecular_weight(),
        analysis.gravy(),
        analysis.isoelectric_point()
    )

def process_sequence(seq_record):
    seq_id = seq_record.id
    sequence = str(seq_record.seq)
    
    orf_len, orf_cov, _ = find_orf(sequence)
    fickett = fickett_value(sequence)
    physico = calculate_physicochemical(sequence)
    
    return seq_id, orf_len, orf_cov, fickett, *physico

def extract_features(fasta_file, conservation_file, output_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    with Pool(cpu_count()) as pool:
        results = pool.map(process_sequence, sequences)
    
    df = pd.DataFrame(results, columns=["gene", "ORF_Max_Len", "ORF_Max_Cov", "Fickett_Score", "Peptides_Length", "Molecular_Weight", "Hydrophobicity", "Isoelectric_Point"])
    
    #conservation_df = pd.read_csv(conservation_file, sep="\t", header=None, names=["gene", "transcript_length", "conservation"])
    conservation_df = pd.read_csv(conservation_file, sep='\t',header=None )
    selected_df = conservation_df.iloc[:, [0,1,5]].rename(
    columns={
        0: 'gene',
        1: 'transcript_length',
        5: 'conservation'
    })
    print("Features DataFrame Columns:", df.columns)
    print("Conservation DataFrame Columns:", selected_df.columns)

    df = df.merge(selected_df, on="gene", how="inner")
    
    df = df[["gene","conservation", "ORF_Max_Len", "ORF_Max_Cov", "Fickett_Score"]]
    df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequence features")
    parser.add_argument("-a", "--seq", required=True, help="Input sequence in FASTA format")
    parser.add_argument("-d", "--conservation", required=True, help="Input conservation annotation file")
    parser.add_argument("-o", "--feature", required=True, help="Output features file")
    
    args = parser.parse_args()
    extract_features(args.seq, args.conservation, args.feature)
