
 # the function ' graph' creates the graph with the constitutive and majiq junctions 
 

import networkx as nx


def graph(gene_id,biomart_genes, majiq_genes):
    
    # first, get majiq junctions of the gene
    junctions = [] # majiq junctions corresponding to the given geneID
    gene = 'gene:' + gene_id
    for j in range(1, len(majiq_genes)):
        if majiq_genes[j] != []:
            if majiq_genes[j][0] == gene:
                juncs = majiq_genes[j][-2]  # all junctions as a string
                psi_values = majiq_genes[j][6].split(";")
                j = juncs.split(';') # split string of junctions
                for i in range(len(j)):
                    alszahl = j[i].split('-') # splitting junction in start & end
                    liste = []
                    liste.append(int(alszahl[0]))
                    liste.append(int(alszahl[1]))
                    liste.append(psi_values[i])
                    junctions.append(liste)


    if (len(junctions) == 0):
        return []


    graph = nx.DiGraph()


    # corresponding exons of the gene 
    exons = [] # exon-strings for the nodes of the graph
    exons_raw = [] # list of exons as [start,end]
    gene_end = 0
    gene_start = 0

    for i in range(1,len(biomart_genes)):
        if biomart_genes[i][0] == gene_id:
            gene_end = int(biomart_genes[i][-1])
            gene_start = int(biomart_genes[i][-2])
            exon_r = []
            exonid = []
            exon = 'start: ' + biomart_genes[i][-6] + ' end: ' + biomart_genes[i][-5] # exon start
            exon_r.append(int(biomart_genes[i][-6])) # exon start
            exon_r.append(int(biomart_genes[i][-5]))  # exon end
            exonid.append(str(biomart_genes[i][2])) # exonid
            exonid.append(int(biomart_genes[i][-6])) # exon start
            exonid.append(int(biomart_genes[i][-5])) # exon end
            exon_r.append(exonid)
            #exon_r.append(biomart_genes[i][-3]) # up or downstream
            if exon not in exons:
                exons.append(exon)
                exons_raw.append(exon_r) # appended: [start,end,[exonid, start, end]]
                # exon_r[2] only necessary for the result of this whole function

    if len(exons_raw) == 0:
        return []

    
    # Step 1: creating graph with constitutive junctions:

    # find 'first' exon (smallest start)
    first_exon_raw = exons_raw[0]
    first_exons_raw = []
    first_exon = exons[0]
    first_exons = []
    for i in range(1,len(exons_raw)):
        if exons_raw[i][0] < first_exon_raw[0]:
            first_exon_raw = exons_raw[i]
            first_exon = exons[i]

    # find overlapping exons of the first
    for i in range(len(exons_raw)):
        if first_exon_raw[0] == exons_raw[i][0] or first_exon_raw[1] == exons_raw[i][1]:
            first_exons_raw.append(exons_raw[i])
            first_exons.append(exons[i])
        elif exons_raw[i][0] in range(first_exon_raw[0], first_exon_raw[1]) and exons_raw[i][1] in range(first_exon_raw[0], first_exon_raw[1]):
            first_exons_raw.append(exons_raw[i])
            first_exons.append(exons[i])




    def makeGraph(exonraw, exonstring):
        nearest = []
        nearest_raw = []
        nearest_one_raw = [0, 0]
        distance = gene_end
        current_distance = 0

        # find one of the nearest exons to exonraw
        for i in range(len(exons_raw)):
            if exons_raw[i][0] < exonraw[1]:
                current_distance = gene_end
            else:
                current_distance = abs(exons_raw[i][0] - exonraw[1])
                if exons_raw[i] != exonraw and current_distance <= distance:  # kann auch < -> ist egal
                    # nimm hier das naeheste exon mit maximaler laenge
                    nearest_one_raw = exons_raw[i]
                    distance = current_distance

        if nearest_one_raw == [0, 0]:
            return

        # take all exons with the same start
        near = []
        for l in range(len(exons_raw)):
            if exons_raw[l][0] == nearest_one_raw[0]:
                near.append(exons_raw[l])

        # compare all exons to the nearest
        # thought behind it: could also take exon with the same start and the biggest length
        # doing it like that to be safe, could be changed
        # otherwise when f.ex. only the shortest: could be wrong results regarding the overlapping exons

        # get all overlapping exons
        for j in range(len(exons_raw)):
            for h in range(len(near)):
                if exons_raw[j][0] == near[h][0]:
                    nearest_raw.append(exons_raw[j])
                    nearest.append(exons[j])
                elif exons_raw[j][1] == near[h][1] and exons_raw[j][0] > exonraw[1]:
                    nearest_raw.append(exons_raw[j])
                    nearest.append(exons[j])

                elif exons_raw[j][0] >= near[h][0] and exons_raw[j][0] <= near[h][1]:
                    nearest_raw.append(exons_raw[j])
                    nearest.append(exons[j])

        # connect all nearest exons to the given
        # repeat for the nearest
        for k in range(len(nearest_raw)):
            graph.add_node(nearest[k], exonid = nearest_raw[k][2])
            graph.add_edge(exonstring, nearest[k])
            makeGraph(nearest_raw[k], nearest[k])


    # for all first exons
    for i in range(len(first_exons_raw)):
        graph.add_node(first_exons[i], exonid = first_exons_raw[i][2])
        makeGraph(first_exons_raw[i], first_exons[i])


        
    # Step 2: add all majiq junctions to the given graph

    
    # start exon = exon with no junction in majiq output with end_exon = start_junction
    start_exons = []
    for i in range(len(exons_raw)):
        for j in range(len(junctions)):
            if exons_raw[i][0] == junctions[j][1]:
                break;
            if j == len(junctions)-1:
                start_exons.append(exons_raw[i])
                start_exons.append(exons[i])



    def creategraph(exon,exonstring):
        junctions_e = []
        # get all junctions starting at that exon
        for i in range(len(junctions)):
            if exon[1] == junctions[i][0]:
                junctions_e.append(junctions[i])


        if len(junctions_e) == 0:
            return

        # get exons where the junctions end
        # call function for thos exons
        for j in range(len(junctions_e)):
            for k in range(len(exons_raw)):
                if junctions_e[j][1] == exons_raw[k][0]:
                    graph.add_edge(exonstring, exons[k], weight = junctions_e[j][2]) #( weight = psivalue)
                    junctions_used.append(junctions_e[j])
                    creategraph(exons_raw[k], exons[k])


    junctions_used = []


    # calling function for starting exons
    for i in range(0,len(start_exons), 2):
        creategraph(start_exons[i], start_exons[i+1])


    for i in range(len(junctions_used)):
        if junctions_used[i] in junctions:
            junctions.remove(junctions_used[i])


    if len(junctions) != 0:
        return ['not annotated']

    # creating the output as: ['GeneID', 'ExonID_From', 'ExonID_To', 'Pos_From', 'Pos_To', 'Psi_value']

    info = [] #output

    # exonid, start, end as a 'exonid' of the node
    # psi values as 'weight' of the edges
    edges = list(graph.edges.data())
    nodes = list(graph.nodes.data())
    for i in range(len(edges)):
        junction = []
        junction.append(gene_id)
        for j in range(len(nodes)):
            if nodes[j][0] == edges[i][0]: # nodefrom
                dict = nodes[j][1]
                node_info = dict['exonid']
                junction.insert(1,node_info[0]) # exonidfrom
                junction.insert(3, node_info[2]) # positionfrom = end of the exon

            if nodes[j][0] == edges[i][1]: #nodeTo
                dict2 = nodes[j][1]
                node_info2 = dict2['exonid']
                junction.insert(2, node_info2[0])  # exonidto
                junction.insert(4, node_info2[1])  # positionto = start of the exon
        attr = edges[i][2] # weight
        if len(attr) != 0:
            weight = attr['weight']
            junction.append(weight)
        else:
            junction.append(0) # = constitutive junction

        info.append(junction)

    return info







