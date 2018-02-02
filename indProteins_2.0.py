#! /usr/bin/python
__author__ = 'mohr'

import os
import sys
import logging
import csv
import re
import vcf
import subprocess
import argparse
import urllib2
import itertools
import pandas as pd
import numpy as np
import Fred2.Core.Generator as generator

from collections import defaultdict
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.ADBAdapter import EIdentifierTypes
from Fred2.IO.UniProtAdapter import UniProtDB
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.IO import FileReader

from datetime import datetime
from string import Template

ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
transcriptProteinMap = {}

VERSION = "2.0"


def get_fred2_annotation(vt, p, r, alt):
    if vt == VariationType.SNP:
        return p, r, alt
    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
        # more than one observed ?
        if alt != '-':
            alternative = '-'
            reference = r[len(alt):]
            position = p + len(alt)
        else:
            return p, r, alt
    elif vt == VariationType.INS or vt == VariationType.FSINS:
        if r != '-':
            position = p
            reference = '-'
            if alt != '-':
                alt_new = alt[len(r):]
                alternative = alt_new
            else:
                alternative = str(alt)
        else:
            return p, r, alt

    return position, reference, alternative


def check_min_req_GSvar(row):
    """
    checking the presence of mandatory columns
    :param row: dictionary of a GSvar row
    :return: boolean, True if min req met
    """
    if ("#chr" in row.keys() and "start" in row.keys() and "end" in row.keys() and "ref" in row.keys() and "obs" in row.keys() and ("coding_and_splicing_details" in row.keys() or "coding" in row.keys() or "coding_and_splicing" in row.keys())):
        return True
    return False


def read_GSvar(filename, pass_only=True):
    """
    reads GSvar and tsv files (tab sep files in context of genetic variants), omitting and warning about rows missing
    mandatory columns
    :param filename: /path/to/file
    :return: list FRED2 variants
    """
    global ID_SYSTEM_USED
    RE = re.compile("(\w+):([\w.]+):([&\w]+):\w*:exon(\d+)\D*\d*:(c.\D*([_\d]+)\D*):(p.\D*(\d+)\w*)")

    metadata_list = ["vardbid", "normal_dp", "tumor_dp", "tumor_af", "normal_af", "rna_tum_freq", "rna_tum_depth"]

    list_vars = list()
    lines = list()
    transcript_ids = []
    dict_vars = {}

    cases = 0

    with open(filename, 'rb') as tsvfile:
        tsvreader = csv.DictReader((row for row in tsvfile if not row.startswith('##')), delimiter='\t')
        for row in tsvreader:
            if not check_min_req_GSvar(row):
                logging.warning("read_GSvar: Omitted row! Mandatory columns not present in: \n"+str(row))
                continue
            lines.append(row)
    for mut_id, line in enumerate(lines):
        if "filter" in line and pass_only and line["filter"].strip():
            continue
        genome_start = int(line["start"]) - 1
        genome_stop = int(line["end"]) - 1
        chrom = line["#chr"]
        ref = line["ref"]
        alt = line["obs"]

        # metadata
        variation_dbid = line.get("dbSNP", '')
        norm_depth = line.get("normal_dp", '')
        tum_depth = line.get("tumor_dp", '')
        tum_af = line.get("tumor_af", '')
        normal_af = line.get("normal_af", '')
        rna_tum_freq = line.get("rna_tum_freq", '')
        rna_tum_dp = line.get("rna_tum_depth", '')

        gene = line.get("gene", '')

        isHomozygous = True if (('tumour_genotype' in line) and (line['tumour_genotype'].split('/')[0] == line['tumour_genotype'].split('/')[1])) else False
        # old GSvar version
        if "coding_and_splicing_details" in line:
            mut_type = line.get("variant_details", '')
            annots = RE.findall(line["coding_and_splicing_details"])
        elif "variant_type" in line:
            mut_type = line.get("variant_type", '')
            annots = RE.findall(line["coding_and_splicing"])
        else:
            mut_type = line.get("variant_details", '')
            annots = RE.findall(line["coding"])
        isyn = mut_type == "synonymous_variant"

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if mut_type == 'missense_variant' or 'missense_variant' in mut_type:
            vt = VariationType.SNP
        elif mut_type == 'frameshift_variant':
            if (ref == '-') or (len(ref) < len(alt)):
                vt = VariationType.FSINS
            else:
                vt = VariationType.FSDEL
        elif mut_type == "inframe_deletion":
            vt = VariationType.DEL
        elif mut_type == "inframe_insertion":
            vt = VariationType.INS

        coding = dict()

        for annot in annots:
            a_gene, nm_id, a_mut_type, exon, trans_coding, trans_pos, prot_coding, prot_start = annot
            if 'NM' in nm_id:
                ID_SYSTEM_USED = EIdentifierTypes.REFSEQ
            if "stop_gained" not in mut_type:
                if not gene:
                    gene = a_gene
                if not mut_type:
                    mut_type = a_mut_type
                nm_id = nm_id.split(".")[0]

                coding[nm_id] = MutationSyntax(nm_id, int(trans_pos.split('_')[0])-1, int(prot_start)-1, trans_coding, prot_coding)
                transcript_ids.append(nm_id)  
        if coding:
            var = Variant(mut_id, vt, chrom.strip('chr'), int(genome_start), ref.upper(), alt.upper(), coding, isHomozygous, isSynonymous=isyn)
            var.gene = gene
            var.log_metadata("vardbid", variation_dbid)
            var.log_metadata("normal_dp", norm_depth)
            var.log_metadata("tumor_dp", tum_depth)
            var.log_metadata("tumor_af", tum_af)
            var.log_metadata("normal_af", normal_af)
            var.log_metadata("rna_tum_freq", rna_tum_freq)
            var.log_metadata("rna_tum_depth", rna_tum_dp)
            dict_vars[var] = var
            list_vars.append(var)

    transToVar = {}

    # fix because of memory/timing issues due to combinatoric explosion
    for v in list_vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        if len(vs) > 10:
            cases += 1
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in metadata_list:
                    vs_new.log_metadata(m, v.get_metadata(m)[0])
                dict_vars[v] = vs_new
    return dict_vars.values(), transcript_ids, metadata_list


def read_vcf(filename, pass_only=True):
    """
    reads vcf files
    returns a list of FRED2 variants
    :param filename: /path/to/file
    :return: list of FRED2 variants
    """
    global ID_SYSTEM_USED

    vl = list()
    with open(filename, 'rb') as tsvfile:
        vcf_reader = vcf.Reader(tsvfile)
        vl = [r for r in vcf_reader]

    dict_vars = {}
    list_vars = []
    transcript_ids = []
    genotye_dict = {"het": False, "hom": True, "ref": True}

    for num, record in enumerate(vl):
        c = record.CHROM.strip('chr')
        p = record.POS - 1
        variation_dbid = record.ID
        r = str(record.REF)
        v_list = record.ALT
        f = record.FILTER
        if pass_only and f:
            continue

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if record.is_snp:
            vt = VariationType.SNP
        elif record.is_indel:
            if len(v_list) % 3 == 0:  # no frameshift
                if record.is_deletion:
                    vt = VariationType.DEL
                else:
                    vt = VariationType.INS
            else:  # frameshift
                if record.is_deletion:
                    vt = VariationType.FSDEL
                else:
                    vt = VariationType.FSINS
        gene = ''

        for alt in v_list:
            isHomozygous = False
            if 'HOM' in record.INFO:
                isHomozygous = record.INFO['HOM'] == 1
            elif 'SGT' in record.INFO:
                zygosity = record.INFO['SGT'].split("->")[1]
                if zygosity in genotye_dict:
                    isHomozygous = genotye_dict[zygosity]
                else:
                    if zygosity[0] == zygosity[1]:
                        isHomozygous = True
                    else:
                        isHomozygous = False
            else:
                for sample in record.samples:
                    if 'GT' in sample.data:
                        isHomozygous = sample.data['GT'] == '1/1'

            if record.INFO['ANN']:
                isSynonymous = False
                coding = dict()
                types = []
                for annraw in record.INFO['ANN']:  # for each ANN only add a new coding! see GSvar
                    annots = annraw.split('|')
                    obs, a_mut_type, impact, a_gene, a_gene_id, feature_type, transcript_id, exon, tot_exon, trans_coding, prot_coding, cdna, cds, aa, distance, warnings = annots
                    types.append(a_mut_type)

                    tpos = 0
                    ppos = 0
                    positions = ''

                    # get cds/protein positions and convert mutation syntax to FRED2 format
                    if trans_coding != '':
                        positions = re.findall(r'\d+', trans_coding)
                        ppos = int(positions[0]) - 1

                    if prot_coding != '':
                        positions = re.findall(r'\d+', prot_coding)
                        tpos = int(positions[0]) - 1

                    isSynonymous = (a_mut_type == "synonymous_variant")

                    gene = a_gene_id
                    # there are no isoforms in biomart
                    transcript_id = transcript_id.split(".")[0]

                    if 'NM' in transcript_id:
                        ID_SYSTEM_USED = EIdentifierTypes.REFSEQ

                    #take online coding variants into account, FRED2 cannot deal with stopgain variants right now
                    if not prot_coding or 'stop_gained' in a_mut_type:
                        continue

                    coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                    transcript_ids.append(transcript_id)

                if coding:
                    pos, reference, alternative = get_fred2_annotation(vt, p, r, str(alt))
                    var = Variant("line" + str(num), vt, c, pos, reference, alternative, coding, isHomozygous, isSynonymous)
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    dict_vars[var] = var
                    list_vars.append(var)

    transToVar = {}

    # fix because of memory/timing issues due to combinatoric explosion
    for v in list_vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        if len(vs) > 10:
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in ["vardbid"]:
                    vs_new.log_metadata(m, v.get_metadata(m))
                dict_vars[v] = vs_new

    return dict_vars.values(), transcript_ids


def __main__():
    parser = argparse.ArgumentParser(description="""Individualized Proteins 2.0 \n Script for generation of protein sequences which contain provided variants.""", version=VERSION)
    parser.add_argument('-s', "--somatic_mutations", help='Somatic variants')
    parser.add_argument('-g', "--germline_mutations", help="Germline variants")
    parser.add_argument('-i', "--identifier", help="<Required> Predictions will be written with this name prefix", required=True)
    parser.add_argument('-r', "--reference", help="Reference, retrieved information will be based on this ensembl version", required=False, default='GRCh37', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('-db', "--database", help="Proteome sequence reference database to be attached to individualized sequences", required=True)
    parser.add_argument('-o', "--output_dir", help="All files written will be put in this directory")

    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    logging.basicConfig(filename=os.path.join(args.output_dir,'{}_indproteinsDB.log'.format(args.identifier)), filemode='w+',
                        level=logging.DEBUG)
    logging.info("Starting generation of protein sequences at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    if args.output_dir is not None:
        try:
            os.chdir(args.output_dir)
            logging.info("Using provided data directory: {}".format(str(args.output_dir)))
        except:
            logging.info("No such directory, using current.")
    else:
        logging.info("Using current data directory.")

    '''start the actual IRMA functions'''
    metadata = []
    #references = {'GRCh37': 'http://grch37.ensembl.org', 'GRCh38': 'http://ensembl.org'}
    references = {'GRCh37': 'http://feb2014.archive.ensembl.org', 'GRCh38': 'http://dec2016.archive.ensembl.org'}
    global transcriptProteinMap

    '''read in variants'''
    if args.somatic_mutations.endswith('.GSvar') or args.somatic_mutations.endswith('.tsv'):
        vl, transcripts, metadata = read_GSvar(args.somatic_mutations)
    elif args.somatic_mutations.endswith('.vcf'):
        vl, transcripts = read_vcf(args.somatic_mutations)

    if args.germline_mutations is not None:
        if args.germline_mutations.endswith('.GSvar') or args.germline_mutations.endswith('.tsv'):
            vl_normal, transcripts_germline, metadata = read_GSvar(args.germline_mutations)
        elif args.germline_mutations.endswith('.vcf'):
            vl_normal, transcripts_germline = read_vcf(args.germline_mutations)

        # combine germline and somatic variants
        vl = vl + vl_normal
        transcripts = transcripts_germline + transcripts

    transcripts = list(set(transcripts))

    # initialize MartsAdapter, GRCh37 or GRCh38 based
    ma = MartsAdapter(biomart=references[args.reference])

    #generate transcripts containing variants, filter for unmutated sequences
    transcripts = [g for g in generator.generate_transcripts_from_variants(vl, ma, ID_SYSTEM_USED) if g.vars]

    #generate proteins from transcripts, table='Standard', stop_symbol='*', to_stop=True, cds=False
    proteins = generator.generate_proteins_from_transcripts(transcripts)

    diff_sequences = {}
    
    out_ref = args.database.split('/')[-1].replace('.fasta','_{}_individualized_protein_DB.fasta'.format(args.identifier))    

    cpRef = 'cp {f} {o}'.format(f=args.database,o=out_ref)
    subprocess.call(cpRef.split())

    with open(out_ref, 'a') as outfile:
        for p in proteins:

            variants = []
            for v in p.vars:
                variants = variants + p.vars[v]

            c = [x.coding.values() for x in variants]
            cf = list(itertools.chain.from_iterable(c))

            cds = ','.join([y.cdsMutationSyntax for y in set(cf)])
            aas = ','.join([y.aaMutationSyntax for y in set(cf)])

            outfile.write('>{}:{}\n'.format(p.transcript_id, aas))
            outfile.write('{}\n'.format(str(p)))

    logging.info("Finished generation of protein sequences at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


if __name__ == "__main__": __main__()
