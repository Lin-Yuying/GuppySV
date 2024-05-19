import sys

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_name = ''
        sequence_data = ''
        for line in file:
            if line.startswith('>'):
                if sequence_name:
                    sequences[sequence_name] = sequence_data
                sequence_name = line.strip()
                sequence_data = ''
            else:
                sequence_data += line.strip()
        if sequence_name:
            sequences[sequence_name] = sequence_data
    return sequences

def write_fasta(sequences, output_file_path):
    with open(output_file_path, 'w') as file:
        for sequence_name, sequence_data in sequences.items():
            file.write(sequence_name + '\n')
            for i in range(0, len(sequence_data), 60):
                file.write(sequence_data[i:i+60] + '\n')

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(sequence))

def invert_region(sequence, start, end):
    region = sequence[start:end]
    inverted_region = reverse_complement(region)
    return sequence[:start] + inverted_region + sequence[end:]

def read_ped(file_path):
    regions = []
    with open(file_path, 'r') as file:
        header = file.readline().strip().split('\t')
        for line in file:
            parts = line.strip().split('\t')
            chrom = parts[header.index("CHROM")]
            start = int(parts[header.index("START")])
            end = int(parts[header.index("END")])
            regions.append((chrom, start, end))
    return regions

def main(input_fasta, ped_file, output_fasta):
    sequences = read_fasta(input_fasta)
    regions = read_ped(ped_file)
    
    for chrom, start, end in regions:
        for sequence_name in sequences:
            if chrom in sequence_name:
                sequences[sequence_name] = invert_region(sequences[sequence_name], start, end)
                break
    
    write_fasta(sequences, output_fasta)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python invert_fasta.py <input_fasta> <ped_file> <output_fasta>")
    else:
        input_fasta = sys.argv[1]
        ped_file = sys.argv[2]
        output_fasta = sys.argv[3]
        main(input_fasta, ped_file, output_fasta)

