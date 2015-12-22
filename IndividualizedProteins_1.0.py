#! /usr/bin/python

import sys
sys.path.append('/lustre_cfc/qbic/chris/Fred2')

import IO
from IO.MartsAdapter import MartsAdapter
from IO.RefSeqAdapter import RefSeqAdapter
from Core.Transcript import Transcript
from Prediction.PSSM import Syfpeithi
from Prediction.NetMHC import NetMHC
from IO.UniProtAdapter import UniProtDB
from Core.Base import AASequence

import time
import timeit
import pickle
import argparse
import datetime


VERSION = 'Version 1.0'

parser = argparse.ArgumentParser(description='Generate Individualized Proteins', prog='iProteins')

parser.add_argument(
        '--basedir', '-b',
        required=True,
        help = 'Base directory of files and for output.'
        )
parser.add_argument(
        '--somatic', '-s',
        required=True,
        help = 'File with somatic mutations.'
        )
parser.add_argument(
        '--germline', '-g',
        required=False,
        help = 'File with germline mutations.'
        )
parser.add_argument(
        '--sampleID', '-d',
        required=True,
        help = 'The ID of the processed sample.'
        )

args = parser.parse_args()

# Starting the actual IRMA functions
mart_db = MartsAdapter()
vl = IO.read_GSvar(args.somatic, args.sampleID + '_s')

if args.germline is None:
        number_of_variants = len(vl)

        for v in vl:
                tmp = v.find_zygosity()
                tmp = v.find_coding()
                tmp = v.find_gene(mart_db) # this could be repl. by get_all_variant_genes

        genes = set(filter(None, [v.gene for v in vl]))

        #~ filter only those with a gene and mutation syntax (coding)
        soon_transcripts = list()
        for g in genes:
                soon_transcripts.append([var for var in vl if var.gene == g and var.coding and 'variant_details' in var.metadata and 'nonsynonymous SNV' in var.metadata['variant_details']])
else:
        vl_normal = IO.read_GSvar(args.germline, args.sampleID + '_g')

        number_of_variants = len(vl) + len(vl_normal)

        for v in vl+vl_normal:
                tmp = v.find_zygosity()
                tmp = v.find_coding()
                tmp = v.find_gene(mart_db) # this could be repl. by get_all_variant_genes

        # What for ?
        for v in vl_normal:
                v.specific = True

        genes = set(filter(None, [v.gene for v in vl+vl_normal]))

        #~ filter only those with a gene and mutation syntax (coding)
        soon_transcripts = list()
        for g in genes:
                soon_transcripts.append([var for var in vl+vl_normal if var.gene == g and var.coding and 'variant_details' in var.metadata and 'nonsynonymous SNV' in var.metadata['variant_details']])

soon_transcripts = filter(None,soon_transcripts)
refseq_db = RefSeqAdapter('/lustre_cfc/qbic/reference_genomes/IRMA_references/RefSeq_human.protein.faa', 66, '/lustre_cfc/qbic/reference_genomes/IRMA_references/RefSeq_human.rna.fna', 66)

ids = mart_db.get_all_variant_ids(genes=list(genes))

gwt = set() # genes without found transcripts
transcripts = list()
for st in soon_transcripts:
        found_id = False
        for id in ids:
                pi, ps, ti, ts = [None]*4
                if 'RefSeq mRNA [e.g. NM_001195597]' in id and 'RefSeq Protein ID [e.g. NP_001005353]' in id and st[0].gene == id['UniProt Gene Name']:
                        ti = id['RefSeq mRNA [e.g. NM_001195597]']
                        pi = id['RefSeq Protein ID [e.g. NP_001005353]']
                        found_id = True
                else:
                        continue
                #TODO find & document which coding NM_s were not found
                for var in st:
                        if ti not in var.coding.keys() and pi not in var.coding.keys():
                                continue
                try:
                        ts = str(refseq_db.get_transcript_sequence(ti).seq)
                        ps = str(refseq_db.get_product_sequence(pi).seq)
                        if ts and ps:
                                tr = Transcript(ti, ts, pi, ps)
                        for var in st:
                                tr.add_variant(var)
                                transcripts.append(tr)
                except:
                        print 'unsuccessful retrieval',ti,refseq_db.get_transcript_sequence(ti),pi,refseq_db.get_product_sequence(pi)
        if not found_id:
                gwt.add(st[0].gene)
for u in gwt:
        print 'Gene had no associated NM/NP pair: ',u

c1 = 0
c2 = 0
peptides = list()
protein_list = list()
variantMapping = {}
with open('%s%s_proteins.fasta' % (args.basedir, args.sampleID), 'w') as indProts:
        for t in transcripts:
                if t.mrna and t.protein:
                        try:
                                protein_list = t.to_proteins()
                                c1 +=1
                        except:
                                c2+= 1
                                exc_type, exc_obj, exc_tb = sys.exc_info()
                        for p in protein_list:
                                indProts.write('>%s\n%s\n' % (p.id.split('|')[1], p.seq))
