from Bio.Blast import NCBIWWW, NCBIXML
import yaml, re, subprocess, os

def create_folders(rawFolder, outputFolder):
    raw_folder = rawFolder
    if not os.path.exists(raw_folder):
        os.makedirs(raw_folder)

    output_folder = outputFolder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  
    
    files_folder = os.path.join(output_folder, 'files')
    if not os.path.exists(files_folder):
        os.makedirs(files_folder)

    figures_folder = os.path.join(output_folder, 'figures')
    if not os.path.exists(figures_folder):
        os.makedirs(figures_folder)

def blastSearch(program, db, seq, order, hit):
    blastResults = {}
    uniqueSpecieSet = set()

    for name, taxID in order.items():
        queryBlast = NCBIWWW.qblast(program = program,
                                    database = db,
                                    sequence = seq,
                                    entrez_query = f'txid{taxID}[ORGN]',
                                    hitlist_size = hit)
        runBlast = NCBIXML.read(queryBlast)
        currentSpecieBlast = {}

        for alignment in runBlast.alignments:
            speciesTitle = alignment.title
            accession = alignment.accession
            speciesMatch = re.search(r'\[(.*?)\]', speciesTitle)

            if speciesMatch:
                species = speciesMatch.group(1)
                for hsp in alignment.hsps:
                    blastSequence = hsp.sbjct 
                    if species not in currentSpecieBlast and len(currentSpecieBlast) < hit:
                        currentSpecieBlast[species] = (accession, blastSequence)
                        uniqueSpecieSet.add(species)

        blastResults[name] = currentSpecieBlast

    return blastResults

def writeBlastResults(blastResults, output):
    with open(output, 'w') as file:
        for order, speciesResults in blastResults.items():
            for species, (accession, sequence) in speciesResults.items():
                file.write(f'>{order}|{species}|{accession}\n')
                file.write(f'{sequence}\n')

def muscle(fasta, output):
    muscleCommand = f'muscle -align {fasta} -output {output}'
    subprocess.run(muscleCommand, shell = True, check = True)

def itpTable(fasta, tableoutput, musclesortedoutput, coordinates):
    order_sequence = {}
    with open(fasta, 'r') as file:
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
        for cys, coord in coordinates.items():
            conserved_list.append(v[coord])
        order_sequence_conserved.append({'accession' : accessionID, 
                                         'order' : order, 
                                         'species' : species, 
                                         'sequence' : v, 
                                         'motifs' : conserved_list})
            
    order_sequence_conserved = sorted(order_sequence_conserved, key=lambda x: x['order'])

    # Order, Species, Accession ID, Yes/No
    with open(tableoutput, 'w') as out:
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

    # Alphabetically organize muscle file
    with open(musclesortedoutput, 'w') as muscleout:
        for blast in order_sequence_conserved:
            muscleout.write(f'>{blast['order']}|{blast['species']}|{blast['accession']}\n{blast['sequence']}\n')

def visualization(muscle_path, weblogo_out, pymsaviz_out):
    
    weblogo_command = f'weblogo -f {muscle_path} -D fasta -o {weblogo_out} -F pdf -n 80 --resolution 400'
    subprocess.run(weblogo_command, shell = True, check = True)
    pymsaviz_command = f'pymsaviz -i {muscle_path} -o {pymsaviz_out} --wrap_length 80 --color_scheme Taylor --show_consensus --show_count'
    subprocess.run(pymsaviz_command, shell = True, check = True)

def main():
    config_file = 'config.yaml'
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    raw = config['itp']['raw_folder']
    output = config['itp']['output_folder']
    create_folders(raw, output)

    # BLAST Variables
    blastProgram = config['itp']['program']
    blastDB = config['itp']['database']
    blastSeq = config['itp']['sequence']
    blastOrder = config['itp']['order']
    blastHit = config['itp']['hit_list']

    blastOutput = config['itp']['outfile']
    muscleOutput = config['itp']['muscleOutput']
    coordinates = config['itp']['itp_conserved_regions']
    sorted_muscleFile = config['itp']['sorted_muscleFile']
    table = config['itp']['tableFile']
    weblogo = config['itp']['weblogo']
    pymsaviz = config['itp']['pymsaviz']

    results = blastSearch(blastProgram, blastDB, blastSeq, blastOrder, blastHit)

    writeBlastResults(results, blastOutput)

    muscle(blastOutput, muscleOutput)

    itpTable(muscleOutput, table, sorted_muscleFile, coordinates)

    visualization(sorted_muscleFile, weblogo, pymsaviz)

if __name__ == '__main__':
    main()