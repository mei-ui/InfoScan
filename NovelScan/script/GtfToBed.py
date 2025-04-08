#!/usr/bin/env python3
import re
import sys
from collections import defaultdict

def parse_attributes(attr_str):
    """解析GTF属性字段，处理带引号和不带引号的值"""
    attributes = {}
    # 匹配键值对，兼容 key "value"; 和 key value; 格式
    pattern = re.compile(r'(\w+)\s+("[^"]+"|\S+);')
    for match in re.finditer(pattern, attr_str):
        key = match.group(1)
        value = match.group(2).strip('"')
        attributes[key] = value
    return attributes

def gtf_to_bed(gtf_path, bed_path, feature_type='exon'):
    """转换GTF到BED格式，支持多exon合并"""
    transcripts = defaultdict(list)
    
    # 第一遍扫描：收集所有exon信息
    with open(gtf_path) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != feature_type:
                continue
                
            chrom = fields[0]
            start = int(fields[3]) - 1  # 转换为0-based
            end = int(fields[4])
            strand = fields[6]
            attr_str = fields[8]
            
            attrs = parse_attributes(attr_str)
            transcript_id = attrs.get('transcript_id', 'NA')
            
            # 存储exon信息用于后续合并
            transcripts[transcript_id].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': attrs.get('gene_id', 'NA'),
                'gene_name': attrs.get('gene_name', attrs.get('gene_id', 'NA')),
                'exon_number': int(attrs.get('exon_number', 0))
            })
    
    # 第二遍处理：生成BED条目
    with open(bed_path, 'w') as fout:
        for transcript_id, exons in transcripts.items():
            if not exons:
                continue
                
            # 按exon_number排序
            exons.sort(key=lambda x: x['exon_number'])
            
            # 提取基本信息
            chrom = exons[0]['chrom']
            strand = exons[0]['strand']
            gene_name = exons[0]['gene_name']
            
            # 计算转录本坐标
            tx_start = min(e['start'] for e in exons)
            tx_end = max(e['end'] for e in exons)
            
            # 构建BED12格式
            block_count = len(exons)
            block_sizes = []
            block_starts = []
            for exon in exons:
                block_sizes.append(str(exon['end'] - exon['start']))
                block_starts.append(str(exon['start'] - tx_start))
            
            bed_entry = [
                chrom,
                str(tx_start),
                str(tx_end),
                f"{transcript_id}",  # Name字段
                "0",                              # Score
                strand,                          # Strand
                str(tx_start),                   # Thick start
                str(tx_end),                     # Thick end
                "0,0,0",                         # RGB
                str(block_count),                # Block count
                ",".join(block_sizes),           # Block sizes
                ",".join(block_starts)           # Block starts
            ]
            
            fout.write('\t'.join(bed_entry) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.gtf output.bed")
        sys.exit(1)
    
    gtf_to_bed(sys.argv[1], sys.argv[2])