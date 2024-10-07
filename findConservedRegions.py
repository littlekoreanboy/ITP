import sys, subprocess

def muscle(input, output):
    muscleCommand = f'muscle -align {input} -output {output}'
    subprocess.run(muscleCommand, shell = True, check = True)

def table(input, output, sort_muscle):
    conserved_motifs = {
        'cys1' : 57,
        'cys2' : 73,
        'cys3' : 76,
        'cys4' : 89,
        'cys5' : 93,
        'cys6' : 102
    }
    order_sequence = {}
    with open(input, 'r') as file:
        order = ''
        sequence = ''
        for lines in file:
            stripped_lines = lines.strip()
            if stripped_lines.startswith('>'):
                if order and sequence:
                    order_sequence[order] = sequence
                order = stripped_lines.split('>')[1]
                sequence = ''
            else:
                sequence += stripped_lines
        if order and sequence:
            order_sequence[order] = sequence

    order_sequence_conserved = []
    for k, v in order_sequence.items():
        order = k.split('|')[0]
        species = k.split('|')[1]
        accessionID = k.split('|')[2]
        conserved_list = []
        for cys, coord in conserved_motifs.items():
            conserved_list.append(v[coord])
        order_sequence_conserved.append({'accession' : accessionID, 
                                         'order' : order, 
                                         'species' : species, 
                                         'sequence' : v, 
                                         'motifs' : conserved_list})

    order_sequence_conserved = sorted(order_sequence_conserved, key=lambda x: x['order'])
    for series in order_sequence_conserved:
        print(series)
    # Order, Species, Accession ID, Yes/No
    with open(output, 'w') as out:
        out.write('Order\tSpecies\tAccession ID\tITP Present\n')
        for blast in order_sequence_conserved:
            count_cys = 0
            for i in blast['motifs']:
                if i == 'C':
                    count_cys += 1
            if count_cys == 6:
                present = 'yes'
                out.write(f"{blast['order']}\t{blast['species']}\t{blast['accession']}\t{present}\n")
            else:
                present = 'no'
                out.write(f"{blast['order']}\t{blast['species']}\t{blast['accession']}\t{present}\n") 

    with open(sort_muscle, 'w') as muscleout:
        for blast in order_sequence_conserved:
            muscleout.write(f'>{blast['order']}|{blast['species']}|{blast['accession']}\n{blast['sequence']}\n')

def visualization(muscle_path, weblogo_out, pymsaviz_out):
    weblogo_command = f'weblogo -f {muscle_path} -D fasta -o {weblogo_out} -F pdf -n 80 --resolution 400'
    subprocess.run(weblogo_command, shell = True, check = True)
    pymsaviz_command = f'pymsaviz -i {muscle_path} -o {pymsaviz_out} --wrap_length 80 --color_scheme Taylor --show_consensus --show_count'
    subprocess.run(pymsaviz_command, shell = True, check = True)

def main():
    fastaFile = 'itp_BLAST_results.fa'
    muscleOut = 'test_muscled_itp_BLAST_results.fa'
    tableOut = 'Table.csv'
    sortmuscle = 'sorted_muscle_file.fa'

    muscleFile = 'sorted_muscled_itp_BLAST_results.fa'
    weblogo = 'weblogo.pdf'
    pymsaviz = 'pymsaviz.pdf'

    if sys.argv[1] == 'muscle':
        muscle(fastaFile, muscleOut)

    if sys.argv[1] == 'table':
        table(muscleOut, tableOut, sortmuscle)

    if sys.argv[1] == 'visualize':
        visualization(muscleFile, weblogo, pymsaviz)

if __name__ == '__main__':
    main()