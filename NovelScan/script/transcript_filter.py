import re
import argparse

def parse_gtf(file_path):
    transcripts = {}
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            transcript_id_match = re.search('transcript_id \"(.*?)\"', attributes)
            class_code_match = re.search('class_code \"(.*?)\"', attributes)

            if not transcript_id_match:
                continue

            transcript_id = transcript_id_match.group(1)
            class_code = class_code_match.group(1) if class_code_match else "unknown"

            # 初始化转录本信息
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    "exons": [],
                    "class_code": class_code,
                    "lines": []
                }

            # 收集外显子信息
            if feature_type == "exon":
                start = int(fields[3])
                end = int(fields[4])
                transcripts[transcript_id]["exons"].append((start, end))
                # 添加调试信息
                #print(f"Adding exon for {transcript_id}: {start}-{end}")

            transcripts[transcript_id]["lines"].append(line)
            
            # 更新class_code（以防在transcript行之后出现）
            if class_code_match:
                transcripts[transcript_id]["class_code"] = class_code

    # 添加调试信息：打印没有外显子的转录本
    for transcript_id, data in transcripts.items():
        if not data["exons"]:
            print(f"Warning: No exons found for transcript {transcript_id}")
            
    return transcripts

def calculate_transcript_length(exons):
    if not exons:  # 如果外显子列表为空，返回长度为0
        return 0
    exons = sorted(exons, key=lambda x: x[0])
    merged_exons = []
    current_start, current_end = exons[0]

    for start, end in exons[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_exons.append((current_start, current_end))
            current_start, current_end = start, end
    merged_exons.append((current_start, current_end))

    return sum(end - start + 1 for start, end in merged_exons)

def filter_and_write(transcripts, output_path, length_threshold):
    filtered_count = 0
    zero_length_count = 0
    with open(output_path, "w") as file:
        for transcript_id, data in transcripts.items():
            transcript_length = calculate_transcript_length(data["exons"])
            if transcript_length == 0:
                zero_length_count += 1
                print(f"Warning: Zero length transcript found: {transcript_id}, Class code: {data['class_code']}")
            if transcript_length > length_threshold and data["class_code"] in {"u", "i", "j"}:
                file.writelines(data["lines"])
                filtered_count += 1

    print(f"\nSummary:")
    print(f"Total transcripts processed: {len(transcripts)}")
    print(f"Transcripts with zero length: {zero_length_count}")
    print(f"Transcripts passing filters: {filtered_count}")

def main():
    parser = argparse.ArgumentParser(description="Filter transcripts from a GTF file.")
    parser.add_argument("input_file", type=str, help="Path to the input GTF file.")
    parser.add_argument("output_file", type=str, help="Path to the output GTF file.")
    parser.add_argument("--length", type=int, default=200, help="Length threshold for filtering transcripts.")

    args = parser.parse_args()

    transcripts_data = parse_gtf(args.input_file)
    filter_and_write(transcripts_data, args.output_file, args.length)

if __name__ == "__main__":
    main()
