import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]
min_associations = sys.argv[4]

barcode_oligo_txt_file = f'./{cells}/{library}{adaptor}/output/DNA_barcode_associations_{min_associations}/oligos_to_barcodes_comb_{library}{adaptor}_30_{min_associations}.txt'
fasta_file = f'./{cells}/{library}{adaptor}/output/DNA_barcode_associations_{min_associations}/oligos_to_barcodes_comb_{library}{adaptor}_30_{min_associations}.fasta'

with open(barcode_oligo_txt_file) as f:
    
    with open(fasta_file, 'w') as outfile:
        count = 0
        for line in f:
            count += 1
            if count == 1:  # skip the header line
                continue
            info = [x.strip() for x in line.split(',')]
            oligo = info[0]
            bc_seq = info[1]
            outfile.write('>'+oligo+'_'+bc_seq)
            outfile.write('\n')
            outfile.write(bc_seq)
            outfile.write('\n')
            
        
    

