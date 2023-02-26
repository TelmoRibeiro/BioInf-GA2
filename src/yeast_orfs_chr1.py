def read_fasta_file(file_name, n):
    dna_seq = ""
    with open(file_name, 'r') as fn:    
        fn.readline()                       # 1st line does not contain base pairs
        dna_seq = fn.read()
        dna_seq = dna_seq.replace('\n', '')
        dna_seq = dna_seq[:n]
    return dna_seq

def base_pair_frequency(dna_seq, dna_seq_len):
    dictionary = {}
    for base_pair in dna_seq:
        if base_pair in dictionary:
            dictionary[base_pair] += 1
        else:
            dictionary[base_pair] =  1
    for base_pair in dictionary:
        dictionary[base_pair] /= (dna_seq_len)
    return dictionary

def gc_content(dictionary):
    return dictionary['G'] + dictionary['C']

def start_codon_number(dna_seq, dna_seq_len):
    count = 0
    for i in range(0, dna_seq_len - 2):
        if dna_seq[i] == 'T' and dna_seq[i + 1] == 'A' and dna_seq[i + 2] == 'C':
            count += 1
    return count

def stop_codon_number(dna_seq, dna_seq_len):
    count = 0
    for i in range(0, dna_seq_len - 2):
        if   dna_seq[i] == 'A' and dna_seq[i + 1] == 'T' and dna_seq[i + 2] == 'T':
            count += 1
        elif dna_seq[i] == 'A' and dna_seq[i + 1] == 'T' and dna_seq[i + 2] == 'C':
            count += 1
        elif dna_seq[i] == 'A' and dna_seq[i + 1] == 'C' and dna_seq[i + 2] == 'T':
            count += 1
    return count

def codon_frequency(dna_seq, dna_seq_len):
    dictionary = {}
    for i in range(0, dna_seq_len - 2):
        codon = dna_seq[i] + dna_seq[i + 1] + dna_seq[i + 2]
        if codon in dictionary:
            dictionary[codon] += 1
        else:
            dictionary[codon] =  1
    return dictionary           

def main():
    dna_seq = read_fasta_file("sequence_chr1.fasta", 30000)
    first   = len(dna_seq)
    second  = base_pair_frequency(dna_seq, first)
    third   = gc_content(second)
    fourth  = start_codon_number(dna_seq, first)
    fifth   = stop_codon_number(dna_seq, first)
    sixth   = codon_frequency(dna_seq, first)

main()