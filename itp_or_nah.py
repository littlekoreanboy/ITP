def readFile(file) -> dict:
    file_dict = {}

    with open(file) as f:
        id = ''
        seq = ''
        for lines in f:
            stripped_lines = lines.strip()
            if stripped_lines.startswith('>'):
                if id and seq:
                    file_dict[id] = seq
                id = stripped_lines
                seq = ''
            else:
                seq += stripped_lines
        if id and seq:
            file_dict[id] = seq

    return file_dict # K = ID, V = Sequence

def countCYS(dict) -> dict:
    cys_dict = {}
    for id, seq in dict.items():
        count = 0
        for aminoAcid in seq:
            if aminoAcid == 'C':
                count += 1
        cys_dict[id] = {count : seq}

    return cys_dict # K = ID, V = {K = Count, V = Sequence}

def cysCoord(seq) -> list:
    cys_found = []
    for position, aminoAcid in enumerate(seq):
        if aminoAcid == 'C' or aminoAcid == 'c':
            cys_found.append(position)
    cys_found.append(len(seq))
    return cys_found

def conservedAA(cys_coord) -> list:
    conserved_aminoAcid = []
    n = 0
    while n < len(cys_coord) - 1:
        x1 = n + 1
        x0 = n
        conserved_aminoAcid.append(cys_coord[x1] - cys_coord[x0])
        n += 1
    return conserved_aminoAcid

def motif(cys_coord, seq) -> list:
    motif = []
    m = 0
    while m < len(cys_coord) - 1:
        pre_motif = ''
        y1 = m + 1
        y0 = m
        for i in range(cys_coord[y0] + 1, cys_coord[y1]):
            pre_motif += seq[i]
        motif.append(pre_motif)
        m += 1

    length_motif = []
    for i in motif:
        length_motif.append(len(i))
    #return motif
    return length_motif

def is_subsequence(pattern, sequence):
    pattern_length = len(pattern)
    for i in range(len(sequence) - pattern_length + 1):
        if sequence[i:i + pattern_length] == pattern:
            return True
    return False

def checkITP(dict) -> dict:
    # Define the conserved pattern of lengths between CYS residues (this can be customized)
    conserved_pattern = [15, 2, 12, 3, 8] 
    
    summary = {}
    for id, count_seq in dict.items():
        for count, seq in count_seq.items():
            if count >= 6:  # At least 6 CYS are required
                cys_found = cysCoord(seq)
                conserved_aminoAcid = conservedAA(cys_found)
                length_motif = motif(cys_found, seq)

                # Check if the conserved pattern is a subsequence within the length motif
                if is_subsequence(conserved_pattern, length_motif):
                    summary[id] = {"Yes" : f"Found {count} CYS residues. The conserved pattern of lengths between CYS residues matched {conserved_pattern}"}  # Sequence has ITP
                else:
                    summary[id] = {"No" : f"Found {count} CYS residues. However, the conserved amino acid length of {length_motif} did not meet the requirement of {conserved_pattern}"}   # Sequence doesn't match conserved pattern
            else:
                summary[id] = {"No" : f"Fewer than 6 CYS found. Found {count} CYS"}   # Fewer than 6 CYS
    return summary

def main():
    itp_BLAST = 'raw/raw_itp_BLAST_results.fa'
    eclosion_BLAST = 'raw/raw_eclosion_hormone_BLAST_results.fa'

    #ITP_dict = readFile(itp_BLAST)
    ITP_dict = readFile(eclosion_BLAST)

    itp_CYS_dict = countCYS(ITP_dict)

    itp_check = checkITP(itp_CYS_dict)

    for id, present in itp_check.items():
        for result, reason in present.items():
            print(f'{id} : {present}')

if __name__ == '__main__':
    main()