import os
import sys
import argparse
import datetime
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

def max_value(numbers):
    return max(numbers)

def min_value(numbers):
    return min(numbers)

def mean_value(numbers):
    return sum(numbers) / len(numbers)

def parse_arguments():
    parser = argparse.ArgumentParser(description='SMART-cluster equivalent in Python')
    parser.add_argument('-t', type=int, default=1, help='Number of threads (default: 1)')
    parser.add_argument('-o', type=str, default='result', help='Output directory (default: result)')
    parser.add_argument('-i', type=float, default=1.5, help='Inflation value for MCL algorithm (default: 1.5)')
    parser.add_argument('-e', type=float, default=0.05, help='Cutoff e-value after BLAST results (default: 0.05)')
    parser.add_argument('-p', type=str, default='primer.txt', help='The primer file')
    return parser.parse_args()

def read_primers(primer_file):
    if not os.path.exists(primer_file):
        print("WARNING: Primer file not found.")
        return '', '', False
    with open(primer_file, 'r') as f:
        lines = f.readlines()
    for_prim = lines[0].strip() if len(lines) > 0 else ''
    back_prim = lines[1].strip() if len(lines) > 1 else ''
    return for_prim, back_prim, True

def run_blast(output_dir, evalue, threads):
    print("***STEP: running blast***")
    makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=f"{output_dir}/uniq_aptamer.fasta", out=f"{output_dir}/goodsequence")
    makeblastdb_cline()
    
    blastn_cline = NcbiblastnCommandline(task="blastn-short", evalue=evalue, db=f"{output_dir}/goodsequence", query=f"{output_dir}/uniq_aptamer.fasta", outfmt=7, out=f"{output_dir}/goodsequence_blastn.out", dust="no", num_threads=threads)
    blastn_cline()

    with open(f"{output_dir}/goodsequence_blastn.out", 'r') as fin, open(f"{output_dir}/goodsequence_v1_blastn.out", 'w') as fout:
        for line in fin:
            if not line.startswith('#'):
                fout.write(line)

def generate_mcl_graph(output_dir):
    print("***STEP: generating mcl graph***")
    similarity = {}
    scores = {}
    
    with open(f"{output_dir}/goodsequence_v1_blastn.out", 'r') as fin:
        for line in fin:
            info = line.strip().split('\t')
            score = float(info[11])
            seq1, seq2 = info[0], info[1]
            
            if seq1 not in similarity:
                similarity[seq1] = {}
            if seq2 not in similarity:
                similarity[seq2] = {}
                
            similarity[seq1][seq2] = similarity[seq1].get(seq2, 0) + score
            similarity[seq2][seq1] = similarity[seq2].get(seq1, 0) + score
            scores[f"{seq1}-{seq2}"] = scores.get(f"{seq1}-{seq2}", 0) + score

    score_values = list(scores.values())
    max_val = max_value(score_values)
    min_val = min_value(score_values)
    mean_val = mean_value(score_values)

    mcl_map = {}
    r_mcl_map = {}
    with open(f"{output_dir}/uniq_aptamer.fasta", 'r') as fin:
        nu = 0
        for line in fin:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
                mcl_map[nu] = seq_id
                r_mcl_map[seq_id] = nu
                nu += 1

    nu -= 1
    with open(f"{output_dir}/graph.txt", 'w') as fout:
        fout.write(f"(mclheader\nmcltype matrix\ndimensions {nu+1}x{nu+1}\n)\n(mclmatrix\nbegin\n\n")
        for i in range(nu+1):
            fout.write(f"{i}")
            for j in similarity.get(mcl_map[i], {}):
                sim_val = similarity[mcl_map[i]][j] / max_val
                fout.write(f"\t{r_mcl_map[j]}:{sim_val}")
            fout.write("\t$\n")
        fout.write(")\n")
    
    return mcl_map, r_mcl_map

def run_mcl(output_dir, inflation, threads):
    print("***STEP: running mcl***")
    os.system(f"mcl {output_dir}/graph.txt -I {inflation} -o {output_dir}/mcl_cluster -te {threads} -V all")

# def parse_mcl_cluster(output_dir):
#     clusters = []
    
#     with open(f"{output_dir}/mcl_cluster", 'r') as file:
#         lines = file.readlines()
    
#     start_parsing = False
#     for line in lines:
#         line = line.strip()
#         if line == 'begin':
#             start_parsing = True
#             continue
#         if line == '$':
#             start_parsing = False
#             continue
#         if start_parsing:
#             parts = line.split()
#             if len(parts) > 1:  
#                 cluster_number = int(parts[0]) + 1  
#                 sequences = [seq for seq in parts[1:] if seq != '$'] 
#                 clusters.append((cluster_number, sequences))
    
#     result = []
#     for index, (cluster_number, sequences) in enumerate(clusters):
#         valid_sequences = [seq for seq in sequences if seq.isdigit()]
#         sequences_plus_one = [str(int(seq) + 1) for seq in valid_sequences]
#         count = len(sequences_plus_one)
#         sequences_str = ':'.join(f"Apt-{seq}" for seq in sequences_plus_one)
#         result.append(f"Clust-{index+1} {count} {sequences_str}")
    
#     with open(f"{output_dir}/aptamer_clusters", 'w') as file:
#         for line in result:  
#             file.write(line + '\n')


def parse_mcl_cluster(output_dir):
    with open(f"{output_dir}/mcl_cluster", 'r') as file:
        content = file.read()

    # Step 1: Find 'begin' and ')' to extract the content between them.
    begin_index = content.find('begin') + len('begin')
    end_index = content.find(')', begin_index)
    
    if begin_index == -1 or end_index == -1:
        raise ValueError("Invalid format: 'begin' or ')' not found.")

    # Extract the content between 'begin' and ')'
    matrix_content = content[begin_index:end_index].strip()

    # Step 2: Split the content by the '$' symbol to get individual groups
    groups = matrix_content.split('$')

    # Result to store the formatted output
    result = []

    for group in groups:
        lines = group.strip().split('\n')
        if not lines:
            continue

        # Step 3: Merge all lines in the current group into a single line
        merged_line = ' '.join(line.strip() for line in lines).strip()

        # Split the merged line into numbers
        numbers = merged_line.split()
        
        if not numbers:
            continue

        # The first number is the group identifier
        group_id = numbers[0]

        # Count frequencies of the rest of the numbers
        freq_count = {}
        for num in numbers[1:]:
            if num in freq_count:
                freq_count[num] += 1
            else:
                freq_count[num] = 1
        
        # Format the output for this group
        freq_string = ':'.join(f"Apt-{int(num)+1}" for num, count in freq_count.items())
        result.append(f"Clust-{int(group_id)+1} {len(freq_count)} {freq_string}")

    with open(f"{output_dir}/aptamer_clusters", 'w') as file:
        for line in result:  
            file.write(line + '\n')

def main():
    args = parse_arguments()

    output_dir = args.o.rstrip('/')
    os.makedirs(output_dir, exist_ok=True)

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"start analysis: {now}")

    for_prim, back_prim, primer_ok = read_primers(args.p)
    if not primer_ok:
        print("Primer file is missing or improperly formatted.")
    
    run_blast(output_dir, args.e, args.t)
    generate_mcl_graph(output_dir)
    run_mcl(output_dir, args.i, args.t)
    parse_mcl_cluster(output_dir)

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"end summary: {now}")

if __name__ == "__main__":
    main()
