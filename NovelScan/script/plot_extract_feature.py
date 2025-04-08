import argparse
import pandas as pd
from plotnine import *
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

def extract_features(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    with Pool(cpu_count()) as pool:
        results = pool.map(process_sequence, sequences)
    
    df = pd.DataFrame(results, columns=["gene", "ORF_Max_Len", "ORF_Max_Cov", "Fickett_Score", "Peptides_Length", "Molecular_Weight", "Hydrophobicity", "Isoelectric_Point"])

    df = df[["gene", "ORF_Max_Len", "ORF_Max_Cov", "Fickett_Score"]]

    return df

def plot_ecdf_plotnine(df, feature, label_col, out_file, x_label=None, xlim=None):
    """
    使用 plotnine 绘制 ECDF 图并保存
    """
    x_label = x_label if x_label else feature
    
    p = (
        ggplot(df, aes(x=feature, color=label_col)) +
        stat_ecdf(alpha=1) +
        xlab(x_label) +
        scale_color_brewer(type='qual', palette='Set1') +
        theme_bw()
    )
    if xlim:
        p += scale_x_continuous(limits=xlim)
    
    p.save(out_file, width=7, height=5, verbose=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--mRNA_seq", required=True, help="Input mRNA FASTA")
    parser.add_argument("-b", "--lncRNA_seq", required=True, help="Input lncRNA FASTA")
    parser.add_argument("-c", "--unlncRNA_seq", required=True, help="Input unlncRNA FASTA")
    parser.add_argument("-d", "--unmRNA_seq", required=True, help="Input unmRNA FASTA")
    parser.add_argument("-e", "--Fickett_Score", required=True, help="Output Fickett_Score Feature")
    parser.add_argument("-f", "--MaxORFlength", required=True, help="Output MaxORFlength Feature")
    parser.add_argument("-g", "--MaxORFcoverage", required=True, help="Output MaxORFcoverage Feature")
    parser.add_argument("-k", "--Fickett_Score_png", required=True, help="Output Fickett_Score Feature")
    parser.add_argument("-i", "--MaxORFlength_png", required=True, help="Output MaxORFlength Feature")
    parser.add_argument("-j", "--MaxORFcoverage_png", required=True, help="Output MaxORFcoverage Feature")
    parser.add_argument("-o", "--feature", required=True, help="Output Feature.txt")
    args = parser.parse_args()

    # 读取序列并提取特征（关键修改点）
    features_mRNA_seq = extract_features(args.mRNA_seq)
    features_mRNA_seq['Label'] = "protein_coding"

    features_lncRNA_seq = extract_features(args.lncRNA_seq)
    features_lncRNA_seq['Label'] = "annotated_lncRNA"

    features_unlncRNA_seq = extract_features(args.unlncRNA_seq)
    features_unlncRNA_seq['Label'] = "unannotated_lncRNA"

    features_unmRNA_seq = extract_features(args.unmRNA_seq)
    features_unmRNA_seq['Label'] = "unannotated_coding_transcripts"

    final_df = pd.concat([features_mRNA_seq, features_lncRNA_seq,features_unlncRNA_seq,features_unmRNA_seq], axis=0, ignore_index=True)

    plot_ecdf_plotnine(final_df, 'Fickett_Score', 'Label', args.Fickett_Score)
    plot_ecdf_plotnine(final_df, 'ORF_Max_Len', 'Label', args.MaxORFlength, xlim=(0, 5000))
    plot_ecdf_plotnine(final_df, 'ORF_Max_Cov', 'Label', args.MaxORFcoverage)

    plot_ecdf_plotnine(final_df, 'Fickett_Score', 'Label', args.Fickett_Score_png)
    plot_ecdf_plotnine(final_df, 'ORF_Max_Len', 'Label', args.MaxORFlength_png, xlim=(0, 5000))
    plot_ecdf_plotnine(final_df, 'ORF_Max_Cov', 'Label', args.MaxORFcoverage_png)

    # 保存结果
    final_df = final_df[["gene","Label", "ORF_Max_Len", "ORF_Max_Cov", "Fickett_Score"]]
    final_df.to_csv(args.feature, sep="\t", index=False)

if __name__ == "__main__":
    main()