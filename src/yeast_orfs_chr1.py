#################### Authors: #################### 
# Diogo Ferreira ## 201805258 ##             MCC # 
# Sara Rescalli  ## 202210943 ##      Mobilidade #  
# Telmo Ribeiro  ## 201805124 ##             MCC #
##################################################



def read_file(file_name, nucleotides):
    dna_seq = ""
    with open(file_name, 'r') as fn:    
        fn.readline()                   # consume first line
        dna_seq = fn.read()             # read file
        fn.close()
    dna_seq = dna_seq.replace('\n', '') # remove all newlines
    dna_seq = dna_seq[:nucleotides]     # keep the first n nucleotides
    return dna_seq

def nucleotide_frequency(dna_seq):
    freq_dictionary = {}
    for nucleotide in dna_seq:
        if nucleotide in freq_dictionary: freq_dictionary[nucleotide] += 1 
        else:                             freq_dictionary[nucleotide]  = 1
    for nucleotide in freq_dictionary:
        freq_dictionary[nucleotide] /= len(dna_seq)
    return freq_dictionary

def gc_content(freq_dictionary):
    return freq_dictionary['G'] + freq_dictionary['C']

def start_codon_number(dna_seq):
    counter = 0
    for i in range(0, len(dna_seq) - 2):
        if dna_seq[i] == 'A' and dna_seq[i+1] == 'T' and dna_seq[i+2] == 'G': counter += 1
    return counter

def stop_codon_number(dna_seq):
    counter = 0
    for i in range(0, len(dna_seq) - 2):
        if   dna_seq[i] == 'T' and dna_seq[i+1] == 'A' and dna_seq[i+2] == 'A': counter += 1
        elif dna_seq[i] == 'T' and dna_seq[i+1] == 'A' and dna_seq[i+2] == 'G': counter += 1
        elif dna_seq[i] == 'T' and dna_seq[i+1] == 'G' and dna_seq[i+2] == 'A': counter += 1
    return counter

def codon_frequency(dna_seq):
    freq_dictionary = {}
    for i in range(1, len(dna_seq) - 2):
        codon = dna_seq[i] + dna_seq[i+1] + dna_seq[i+2]
        if codon in freq_dictionary: freq_dictionary[codon] += 1
        else:                        freq_dictionary[codon]  = 1
    return (min(freq_dictionary, key=freq_dictionary.get), max(freq_dictionary, key=freq_dictionary.get))
    
def translate_codon(codon):
    translation_table = {
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "TGT":"C", "TGC":"C",
        "GAT":"D", "GAC":"D",
        "GAA":"E", "GAG":"E",
        "TTT":"F", "TTC":"F",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        "CAT":"H", "CAC":"H",
        "ATA":"I", "ATT":"I", "ATC":"I",
        "AAA":"K", "AAG":"K",
        "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "ATG":"M", 
        "AAT":"N", "AAC":"N",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "TGG":"W",
        "TAT":"Y", "TAC":"Y",
        "TAA":"_", "TAG":"_", "TGA":"_"
    }
    if codon in translation_table: return translation_table[codon]
    return None

def translate_seq(dna_seq, initial_position):
    translated_seq = ""
    i = initial_position
    while i < len(dna_seq) - 2:
        codon = dna_seq[i] + dna_seq[i+1] + dna_seq[i+2]
        translated_seq += translate_codon(codon)
        i = i + 3
    return translated_seq

def nucleotide_pair(nucleotide):
    pair_table = {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"
    }
    if nucleotide in pair_table: return pair_table[nucleotide]
    return None

def reverse_complement(dna_seq):
    rc = ""
    i  = len(dna_seq) - 1
    while i >= 0:
        rc += nucleotide_pair(dna_seq[i])
        i   = i - 1 
    return rc

def reading_frames(dna_seq):
    rfs = []
    rc = reverse_complement(dna_seq)
    rfs.append(translate_seq(dna_seq, 0))
    rfs.append(translate_seq(dna_seq, 1))
    rfs.append(translate_seq(dna_seq, 2))
    rfs.append(translate_seq(rc,      0))
    rfs.append(translate_seq(rc,      1))
    rfs.append(translate_seq(rc,      2))
    return rfs

def all_proteins_rf(aa_seq):
    current_proteins = []
    proteins = []
    for aa in aa_seq:
        if aa == '_':
            if current_proteins:
                for p in current_proteins:
                    proteins.append(p)
                current_proteins = []
        else:
            if aa == 'M':
                current_proteins.append("")
            for i in range(len(current_proteins)):
                current_proteins[i] += aa
    return proteins

def insert_protein_ord(protein, proteins_list):
    i = 0
    while i < len(proteins_list) and len(protein) < len(proteins_list[i]):
        i = i + 1
    proteins_list.insert(i, protein)

def all_orfs_ord(dna_seq, minimum_size):
    proteins_list = []
    rfs = reading_frames(dna_seq)
    for aa_seq in rfs:
        proteins = all_proteins_rf(aa_seq)
        maximum_size = len(proteins[0])
        maximum_protein = proteins[0] 
        for p in proteins:
            if len(p) > maximum_size:
                maximum_size = len(p)
                maximum_protein = p 
        if maximum_size >= minimum_size: insert_protein_ord(maximum_protein, proteins_list)
    return proteins_list    

def main():
    dna_seq = read_file("sequence_chr1.fasta", 30000)
    print(f"Sequence Length: {len(dna_seq)}")                                       # 1st exercise
    freq_dictionary = nucleotide_frequency(dna_seq)
    for n in freq_dictionary:
        print(f"{n}: {freq_dictionary[n]}\t", end = " ")                            # 2nd exercise
    print()
    print(f"GC Content: {gc_content(freq_dictionary)}")                             # 3rd exercise
    print(f"Number of Start Codons: {start_codon_number(dna_seq)}")                 # 4th exercise
    print(f"Number of Stop  Codons: {stop_codon_number(dna_seq)}")                  # 5th exercise
    tuple = codon_frequency(dna_seq)
    print(f"Least Frequent Codons: {tuple[0]}\tMost Frequent Codons: {tuple[1]}")   # 6th exercise
    
    seven   = all_orfs_ord(dna_seq, 50)
    print(len(seven))
    with open("all_potential_proteins.txt", "w") as fn:
        for orf in seven:
            fn.write(orf + "\n")
        fn.close()

main()