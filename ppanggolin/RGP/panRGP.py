#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import argparse

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome, Region
from ppanggolin.formats import checkPangenomeInfo, writePangenome

class MatriceNode:
    def __init__(self, state, score, prev, gene):
        self.state = state  # state of the node. 1 for RGP and 0 for not RGP.
        self.score = score if score > 0 else 0  # current score of the node
        self.prev = prev  # previous matriceNode
        self.gene = gene  # gene this node corresponds to

    def changes(self, state, score):
        # state of the node. 1 for RGP and 0 for not RGP.
        self.state = 1 if score >= 0 else 0
        # current score of the node. If the given score is negative, set to 0.
        self.score = score if score >= 0 else 0

def extractRGP(contig, node, ID):
    """
        Extract the region from the given starting node
    """
    new_region = Region(contig.name + "#" + str(ID))
    while node.state:
        new_region.append(node.gene)
        node.state = 0
        node.score = 0
        node = node.prev
        if node is None:#it's the end of the contig and the end of the region.
            break
    return new_region

def rewriteMatrix(contig, matrix, index, persistent, continuity, multi):
    """
        ReWrite the matrice from the given index of the node that started a region.
    """
    prev = matrix[index]
    index += 1
    if index > len(matrix) and contig.is_circular:
        index = 0
    # else the node was the last one of the contig, and there is nothing to do
    if index < len(matrix):
        nextNode = matrix[index]
        nbPerc = 0
        while nextNode.state:  # while the old state is not 0, recompute the scores.
            if nextNode.gene.family.namedPartition == "persistent" and nextNode.gene.family not in multi:
                modif = -pow(persistent, nbPerc)
                nbPerc += 1
            else:
                modif = continuity
                nbPerc = 0

            curr_score = modif + prev.score
            curr_state = 1 if curr_score >= 0 else 0
            # scores can't be negative. If they are, they'll be set to 0.
            matrix[index].changes(curr_state, curr_score)
            index += 1
            if index >= len(matrix):
                if contig.is_circular:
                    index = 0
                else:
                    # else we're at the end of the contig, so there are no more computations. Get out of the loop
                    break

            prev = nextNode
            nextNode = matrix[index]

def initMatrices(contig, persistent_penalty, variable_gain, multi ):
    """initialize the vector of score/state nodes"""
    mat = []
    prev = None
    nbPerc = 0
    for gene in contig.genes:
        if gene.family.namedPartition == "persistent" and gene.family not in multi:
            modif = -pow(persistent_penalty, nbPerc)
            nbPerc += 1
        else:
            modif = variable_gain
            nbPerc = 0

        curr_score = modif + prev.score if prev is not None else modif
        curr_state = 1 if curr_score >= 0 else 0
        prev = MatriceNode(curr_state, curr_score, prev, gene)
        mat.append(prev)

    # if the contig is circular and we're in a rgp state, we need to continue from the "starting" gene until we leave rgp state.
    if contig.is_circular and curr_state:
        # the previous node of the first processed gene is the last node.
        mat[0].prev = prev
        lastNode = prev  # saving the last node that was inserted.
        curr_score = prev.score
        c = 0
        nbPerc = 0
        while curr_state:  # while state is rgp.
            matNode = mat[c]
            if matNode == lastNode:  # then we've parsed the entire contig twice. The whole sequence is a rgp so we're stopping the iteration now, otherwise we'll loop indefinitely
                break

            if matNode.gene.family.namedPartition == "persistent" and matNode.gene.family not in multi:
                modif = -pow(persistent_penalty, nbPerc)
                nbPerc += 1
            else:
                modif = variable_gain
                nbPerc = 0

            curr_score = modif + prev.score
            curr_state = 1 if curr_score >= 0 else 0
            matNode.changes(curr_state, curr_score)
            c += 1
    return mat

def _writeMat_debug(matrix, outname):
        """
            For dev, function to call between steps when you think there is a problem.
        """
        f = open(outname, "w")
        for matriceNode in matrix:
            f.write(matriceNode.gene.ID + "\t" + matriceNode.gene.family.namedPartition +
                    "\t" + str(matriceNode.score) + "\t" + str(matriceNode.state) + "\n")
        f.close()

def mkRegions(contig, matrix, min_length, min_score, persistent, continuity, multi):
    # processing matrix and 'emptying' it to get the regions.
    def maxIndexNode(lst):
        """gets the last node with the highest score from a list of matriceNode"""
        if isinstance(lst, list):
            # init with the first element of the list
            maxScore = lst[0].score
            maxIndex = 0
            for index, node in enumerate(lst):
                if node.score >= maxScore:
                    maxScore = node.score
                    maxIndex = index
            return (maxScore, maxIndex)
        else:
            raise TypeError("List of matriceNode is expected. The detected type was " + type(lst))

    contigRegions = set()
    val, index = maxIndexNode(matrix)
    c = 0
    cc = 0
    _writeMat_debug(matrix, 'tototo/' + contig.name + "#" + str(cc))
    cc+=1
    while val >= min_score:
        new_region = extractRGP(contig, matrix[index], len(contigRegions))
        new_region.score = val
        if (new_region[0].stop - new_region[-1].start) > min_length:
            contigRegions.add(new_region)
        rewriteMatrix(contig, matrix, index, persistent, continuity, multi)
        _writeMat_debug(matrix, 'tototo/' +contig.name + "#" + str(cc))
        cc+=1
        val, index = maxIndexNode(matrix)
        c += 1
    return contigRegions

def compute_org_rgp(organism, persistent_penalty, variable_gain, min_length, min_score, multigenics):
    orgRegions = set()
    for contig in organism.contigs:
        matrix = initMatrices(contig, persistent_penalty, variable_gain, multigenics)
        orgRegions |= mkRegions(contig, matrix, min_length, min_score, persistent_penalty, variable_gain, multigenics)

    return orgRegions

def get_multigenics(pangenome, dup_margin):
    """
        Returns the multigenic persistent families of the pangenome graph. A family will be considered multigenic if it is duplicated in more than 5% of the genomes where it is present.
    """
    multigenics = set()
    for fam in pangenome.geneFamilies:
        if fam.namedPartition == "persistent":
            dup=len([genes for org, genes in fam.getOrgDict().items() if len(genes) > 1])
            if (dup / len(fam.organisms)) >= dup_margin:#tot / nborgs >= 1.05
                multigenics.add(fam)
    logging.getLogger().info(f"{len(multigenics)} gene families are defined as being multigenic. (duplicated in more than {dup_margin} of the genomes)")
    return multigenics

def predictRGP(pangenome, persistent_penalty = 3, variable_gain = 1, min_length = 3000, min_score = 4, dup_margin = 0.05, cpu = 1):
    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=True, needPartitions = True)

    logging.getLogger().info("Detecting multigenic families...")
    multigenics = get_multigenics(pangenome, dup_margin)
    pangenomeRGP = set()
    logging.getLogger().info("Compute Regions of Genomic Plasticity ...")
    bar = tqdm(pangenome.organisms, unit = "genomes")
    for org in bar:
        pangenomeRGP |= compute_org_rgp(org, persistent_penalty, variable_gain, min_length, min_score, multigenics)
    logging.getLogger().info(f"Predicted {len(pangenomeRGP)} RGP")
    #save parameters and save status
    pangenome.parameters["RGP"] = {}
    pangenome.parameters["RGP"]["persistent_penalty"] = persistent_penalty
    pangenome.parameters["RGP"]["variable_gain"] = variable_gain
    pangenome.parameters["RGP"]["min_length"] = min_length
    pangenome.parameters["RGP"]["min_score"] = min_score
    pangenome.parameters["RGP"]["dup_margin"] = dup_margin
    pangenome.status['predictedRGP'] = "Computed"

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    predictRGP(pangenome, args.persistent_penalty, args.variable_gain, args.min_length, args.min_score, args.dup_margin, args.cpu)

def rgpSubparser(subparser):
    parser = subparser.add_parser("rgp", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--persistent_penalty', required=False, type=int, default=3, help="Penalty score to apply to persistent genes")
    optional.add_argument('--variable_gain', required=False, type=int, default=1, help="Gain score to apply to variable genes")
    optional.add_argument('--min_score', required=False, type=int, default=4, help="Minimal score wanted for considering a region as being a RGP")
    optional.add_argument('--min_length', required=False, type=int, default=3000, help="Minimum length (bp) of a region to be considered a RGP")
    optional.add_argument("--dup_margin", required = False, type=int, default=0.05, help="Minimum ratio of organisms where the family is present in which the family must have multiple genes for it to be considered 'duplicated'" )
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser