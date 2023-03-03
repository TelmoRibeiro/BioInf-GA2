#################### Authors: #################### 
# Diogo Ferreira ## 201805258 ##             MCC # 
# Sara Rescalli  ## 202210943 ##      Mobilidade #  
# Telmo Ribeiro  ## 201805124 ##             MCC #
##################################################

import sys 

def read_file(file_name, nucleotides):
    dna_seq = ""
    with open(file_name, 'r') as fn:    
        fn.readline()                   # consume first line
        dna_seq = fn.read()             # read file
        fn.close()
    dna_seq = dna_seq.replace('\n', '') # remove all newlines
    dna_seq = dna_seq[:nucleotides]     # keep the first n nucleotides
    return dna_seq

def validate_dna (dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    if valid == len(seqm): return True
    else: return False

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

def translate_seq (dna_seq):
    """ Translates a DNA sequence into an aminoacid sequence. """
    
    seq_aa = ""

    # Generate aminoacid string
    for i in range(0, len(dna_seq) - 3 - 1, 3):
        seq_aa += translate_codon(dna_seq[i:i+3])
    return seq_aa

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

#####################################
# READING FRAMES AS CODONS
def reading_frame (dna_seq, ini_pos):
    """ Obtains a reading frame for a given DNA sequence. """
    # ini_pos = 0 > frame 1
    # ini_pos = 1 > frame 2
    # ....
    rf = ""

    for i in range(ini_pos, len(dna_seq), 3):
        rf += dna_seq[i:i+3]
    return rf

def all_reading_frames (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement)."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    
    dna_seq = dna_seq.upper()

    res = []
    res.append(reading_frame(dna_seq,0))
    res.append(reading_frame(dna_seq,1))
    res.append(reading_frame(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(reading_frame(rc,0))
    res.append(reading_frame(rc,1))
    res.append(reading_frame(rc,2))
    return res
#####################################

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

def all_orfs_ord(rfs, minimum_size):
    proteins_list = []
    
    aa_rfs = []

    for r in rfs: 
        aa_rfs.append(translate_seq(r))

    for aa_seq in aa_rfs:
        proteins = all_proteins_rf(aa_seq)

        maximum_size = 0
        maximum_protein = ""
        if(len(proteins) > 0): 
            maximum_size = len(proteins[0])
            maximum_protein = proteins[0] 

        for p in proteins:
            if len(p) > maximum_size:
                maximum_size = len(p)
                maximum_protein = p 
        if maximum_size >= minimum_size: insert_protein_ord(maximum_protein, proteins_list)
    return proteins_list    

##############################################
# EXERCISE 8
def convert_coords_rf(c, nrframe, seq_length):
    coord = 0
    if nrframe >= 1 and nrframe <= 3: # positive strand
        coord = c + nrframe
    else: # negative strand
        coord = seq_length - (c + (nrframe - 3 ) - 1) 
    return coord

def find_orfs_coords(dna_seq, minsize, rfs):
    """ Function to find coordinates of genomic sequence for all ORFs. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    
    # set start and stop codons
    start_codons = ['ATG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    # create list to store ORF coordinates
    orfs_coords = []
    
    nrframe = 0 # reading frame number 1-6

    for rf in rfs:
        start_index = 0
        first = False # to select longest ORF, if alternative start codons are found
        end_index = 0
        nrframe += 1  
        for i in range(0, len(rf), 3):
            if rf[i:i+3] in start_codons and not first:
                start_index = convert_coords_rf(i, nrframe, len(dna_seq))
                first = True

            elif rf[i:i+3] in stop_codons and first:
                end_index = convert_coords_rf(i+2, nrframe, len(dna_seq))

                if abs(end_index - start_index) >= minsize:
                    orfs_coords.append([start_index, end_index, nrframe])

                start_index = 0
                first = False
                end_index = 0
                
    rf_dict = {1:"+1", 2:"+2", 3:"+3", 4:"-1", 5:"-2", 6:"-3"}

    file = open('orf_coordinates.txt','w')
    i = 0
    for orf in orfs_coords: 
        orf_start = str(orf[0])
        orf_end = str(orf[1])
        length = str(abs(orf[1] - orf[0]) + 1)
        i += 1
        file.write(orf_start + ", " + orf_end + ", " + "ORF" + str(i) + "    " + rf_dict[orf[2]] + ", " + " length:" + length + "\n")
    file.close()

    # return list of ORFs & coordinates
    return orfs_coords
##############################################

##############################################
# EXERCISE 9 
def overlap(orfs_coords, gtf_file):
    """ Function to compare the results obtained with those from the annotation """
    # Consider only entries of the type Exon (3rd field)
    exons = []

    fh = open(gtf_file, "r")
    for line in fh:
        line = line.split()
        if(line[2] == "exon"):
            exons.append(line)
    #print(exons)
    fh.close()

    # Find Longest Overlap
    for e in exons:
        longest_overlap = 0
        gene_id = "" 
        a = [*range(int(e[3]), int(e[4]) + 1)] # Convert genomic region to list

        for o in orfs_coords:
            if(o[0] > o[1]):
                b = [*range(o[0], o[1] - 1, -1)]
            else:
                b = [*range(o[0], o[1] + 1)]

            overlap = len(list(set(a) & set(b))) / len(list(set(a)))

            if overlap >= longest_overlap:
                longest_overlap = overlap
                result = f'{longest_overlap:.1%}' # formatted percentage
        gene_id = e[9].split('"')[1]
        print(gene_id + "    ", result)
##############################################

def main():
    genome = sys.argv[1]    
    dna_seq = read_file(genome, 30000)
    rfs = all_reading_frames(dna_seq)

    for rf in rfs:
        print(f"Sequence Length: {len(rf)}")                                       # 1st exercise
        freq_dictionary = nucleotide_frequency(rf)
        for n in freq_dictionary:
            print(f"{n}: {freq_dictionary[n]}\t", end = " ")                       # 2nd exercise
        print()
        print(f"GC Content: {gc_content(freq_dictionary)}")                        # 3rd exercise
        print(f"Number of Start Codons: {start_codon_number(rf)}")                 # 4th exercise
        print(f"Number of Stop  Codons: {stop_codon_number(rf)}")                  # 5th exercise
        tuple = codon_frequency(rf)
        print(f"Least Frequent Codons: {tuple[0]}\tMost Frequent Codons: {tuple[1]}")   # 6th exercise
        print()
    
    seven = all_orfs_ord(rfs, 1)
    # print(len(seven))
    with open("all_potential_proteins.txt", "w") as fn:
        for orf in seven:
            fn.write(orf + "\n")
        fn.close()

    orf_coords = find_orfs_coords(dna_seq, 150, rfs)                     # 8th exercise
    annot = sys.argv[2]
    overlap(orf_coords, annot)                                           # 9th exercise

main()