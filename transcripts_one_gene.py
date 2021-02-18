

# output of 'transcripts()': transcripts as lists of exonids


def transcripts(junctions):
    

    exon_ids_from = []
    exon_ids_to = []
    for i in range(len(junctions)):
        exon_ids_from.append(junctions[i][1])
        exon_ids_to.append(junctions[i][2])


    # starting exon when no junction with exonidto = exonid
    starting_exons = []
    for i in range(len(exon_ids_from)):
        if exon_ids_from[i] not in exon_ids_to:
            starting_exons.append(exon_ids_from[i])



    exonids = [] # transcripts as lists of exonids


    def create_transcripts(exons):
        corr_junctions = []

        # get all junctions corresponding to the last exon of 'exons'
        for i in range(len(junctions)):
            if exons[len(exons)-1] == junctions[i][1]:
                corr_junctions.append(junctions[i])

        if len(corr_junctions) == 0: # no junction with exonid = exonidfrom => return
            exonids.append(exons)
            return

        for j in range(len(corr_junctions)):
            exons.append(corr_junctions[j][2])
            create_transcripts(exons)


    # for all starting exons
    for i in range(len(starting_exons)):
        exons = []
        exons.append(starting_exons[i])
        create_transcripts(exons)


    return exonids


