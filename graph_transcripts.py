
# this script should be run to get the transcripts as a fasta file of a gene 
# geneId and path, to which the fasta file is saved, is being asked as an input 
# as needed files: - tsv file with all genes
#                  - majiq output
#                  - fasta file with exon sequences
#                  - fasta file with transcripts
# more details about the files below
# this script uses the scripts 'graph_onegene.py' and transcripts_one_gene.py' where no additional files are needed

import csv
from graph_onegene import graph
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



def transcripts():

    all_genes = []
    # mart_export.tsv is a tsv file with all genes from biomart 
    # the columns are: Gene stable ID, Gene stable ID version, Exon stable ID, Exon region start (bp), Exon region end (bp),Transcription start site (TSS), Strand, Gene start(bp), Gene end (bp)
    with open("/home/sarah/Downloads/Praktikum/mart_export.tsv", "r") as output:
        mart_output = csv.reader(output, delimiter="\t")
        for rec in mart_output:
            all_genes.append(rec)

    #print(all_genes[0])
    #print(all_genes[1])

    gene_id = input("Gene ID: ")
    path = input("Path: ")

    majiq_genes = [] # majiq junctions
    
    # Sclerotic_No_Category.deltapsi.tsv is the majiq ouput
    with open("/home/sarah/Downloads/Praktikum/Sclerotic_No_Category.deltapsi.tsv", "r") as output:
        majiq_output = csv.reader(output, delimiter="\t")
        for record in majiq_output:
            majiq_genes.append(record)

    # get all junctions of the gene
    gene_info = graph(gene_id,all_genes, majiq_genes)

    if len(gene_info) == 0:
        print("No majiq junctions corresponding to the gene")
        return
    if gene_info[0] == "not annotated":
        print("Not annotated exons required")
        return


    #print(len(gene_info))
    #print(gene_info)

    # generate transcripts
    from transcripts_one_gene import transcripts

    gene_trans = []  # = [gene_id, [[junction1], [junction2],..., [junctionN]]]
    gene_trans.append(gene_info[0][0])
    transcript = transcripts(gene_info) # = [[exonID1,exonID2,exonID3,...],[..],...,[..]]
    gene_trans.append(transcript)
    #print(gene_trans)


    # fasta file for the sequences
    records2 = []
    
    # martquery_0128151332_8.fasta = fasta file from biomart with exon sequences 
    # an entry looks f.ex. like follows: 
    #>ENSG00000271254|4612|29626|ENSE00003714086|10911|11083|-1
    #GCTCAGCAGGGAGCTGCTGGATGAGAAAGGGCCTGAAGTCTTGCAGGACTCACTGGATAG(...)
    # geneID, gene start, gene end, exonID, exon start, exon end, strand, sequence
    with open("/home/sarah/Downloads/Praktikum/martquery_0128151332_8.fasta", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records2.append(record)

    records = [] # all exons+sequences of the gene
    for i in range(len(records2)):
        exon = [] # = [geneID, exonID, sequence]
        inform = records2[i].id.split("|")
        exon.append(inform[0])
        exon.append(inform[3])
        exon.append(str(records2[i].seq))
        records.append(exon)

    #print(records[0])

    # for one gene: do not need the code as a function
    # -> is only there because the code before was for multiple genes
    # could change that
    def sequences(gene_id, transcripts):
        transcripts2 = []
        transcripts2.append(gene_id)

        ex_seq = {} #key=exonid : value=sequence
        # correponing exons of that gene
        for i in range(len(records)):
            if records[i][0] == gene_id:
                ex_seq[records[i][1]] = records[i][2]


        for j in range(len(transcripts)):  # [[trans1], [transcript2], ...[]]
            sequence = ''
            for l in range(len(transcripts[j])):  # [exonid, exonid2,....]
                sequence = sequence + ex_seq[transcripts[j][l]]

            transcripts2.append(sequence)

        return transcripts2

    # generate transcript sequences
    g_t = sequences(gene_trans[0], gene_trans[1])
    #print(len(g_t))


    # comparing to annotated transcripts and generating the output
    ann_transcripts = []
    
    # fasta file from biomart with transcripts 
    # an entry is: geneID, exonID, transcript
    with open("/home/sarah/Downloads/Praktikum/martquery_allTranscripts.fasta", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ann_transcripts.append(record)


    ann_or_new = []
    for i in range(1,len(g_t)):
        for k in range(len(ann_transcripts)):
            if ann_transcripts[k].id[0:15] == g_t[0]:  # same gene
                if str(ann_transcripts[k].seq) in g_t[i]:  # if substring
                    ann_or_new.append('annotated')
                    break
            if k == len(ann_transcripts) - 1:
                ann_or_new.append('not annotated')

   # print(len(ann_or_new))
   # print(ann_or_new[0])

    sequences = []
    # creating the sequence records
    for i in range(1,len(g_t)):
        record = SeqRecord(
                Seq(g_t[i]),
                id=g_t[0],  # geneid
                description=ann_or_new[i-1]
            )
        sequences.append(record)

    #print(len(sequences))
    #print(sequences[0])

    
    with open(path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")


transcripts()
