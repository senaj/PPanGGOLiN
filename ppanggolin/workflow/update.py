#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import time
import argparse
import logging

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.rarefaction import makeRarefactionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, writeFlatFiles
from ppanggolin.figures import drawTilePlot, drawUCurve
from ppanggolin.info import printInfo
from ppanggolin.RGP.genomicIsland import predictRGP
from ppanggolin.RGP.spot import predictHotspots
# a workflow to add genomes to an existing pangenome


def launch(args):
    pass

def updateSubparser(subparser):
    parser = subparser.add_parser("update", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone.")
    parser.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    #required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")
    parser.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")

    return parser
