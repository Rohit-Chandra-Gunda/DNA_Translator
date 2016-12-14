def translate(sequence):
    """
    Translates a given nucleotide sequence to a string or sequence of amino acids.
    returns a string of sequence of amino acids.
    :param sequence:
    :return:
    """
    valid = (len(sequence) % 3 == 0)
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W', }
    protein = ""
    if valid:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            try:
                protein += table[codon]
            except:
                valid = False
                break
    if valid:
        return protein
    else:
        print("The sequence of nucleotides is not valid.")
        return "INVALID SEQUENCE OF NUCLEOTIDES."

def read_seq(input_file):
    f = open(input_file, "r")
    sequence = f.read()
    f.close()
    sequence = sequence.replace("\n", "")
    sequence = sequence.replace("\r", "")
    return sequence

nucleotide_sequence = read_seq("dna.txt")
l1 = int(input("Enter the beginning position of the sequence."))
l2 = int(input("Enter the ending position of the sequence."))
nucleotide_sequence = nucleotide_sequence[l1-1:l2-3]
protein_translated = translate(nucleotide_sequence)
f = open("protein.txt", "w")
f.write(protein_translated)
f.close()
