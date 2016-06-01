#!/usr/bin/env python3
import re
import fileinput
import functools
import textwrap
import multiprocessing as mp
import argparse
import queue
import sys
import tempfile
import collections
from operator import itemgetter
from typing import List
from typing import Dict

__author__ = 'uvesten'

# definitions and regexp used globally (for now)
start_codons = ['ttg', 'ctg', 'att', 'atc', 'ata', 'atg', 'gtg']
stop_codons = ['taa', 'tag', 'tga']

start_regexp = re.compile('at.|[atcg]tg')
stop_regexp = re.compile('ta[ga]|tga')


shine_regexp = re.compile('a?g?gagg|ggag|gg.{1}gg')


gc_translation = str.maketrans("cg", "11", "atn")

# Data types
# named tuple to hold configuration data
GpConfig = collections.namedtuple(
    'GpConfig',
    'shortest_gene, allow_runoff, intra_gene_gap, shine_box_distance')

# named tuple for the results of gene finding
GeneData = collections.namedtuple(
    'GeneData',
    'current_gene_no, start_pos, end_pos, gc_content, shine_box, strand')

# named tuple for working on a scaffold
WorkUnit = collections.namedtuple(
    'WorkUnit',
    'scaffold_name, scaffold_sequence, tot_gc_pct, config')

GCFractions = collections.namedtuple('GCFractions', 'tot_gc, gcp1, gcp2, gcp3')

ScaffoldData = collections.namedtuple(
    'ScaffoldData', 'scaffold_sequence, genes_data')

complement = str.maketrans("actg", "tgac")

# utility functions


def reverse_complement(dna):
    """
    Generates the reverse complement of a dna sequence
    """

    dna_comp = dna.translate(complement)
    return dna_comp[::-1]


def get_gc_content(dna):
    """
    Gets the total gc content, and for positions 1-3
    of a DNA sequence
    """

    dna_length = len(dna)

    gc1 = len(dna[::3].translate(gc_translation))
    gc2 = len(dna[1::3].translate(gc_translation))
    gc3 = len(dna[2::3].translate(gc_translation))

    # calculate total gc content

    tot_gc = (gc1 + gc2 + gc3) / dna_length
    part_length = dna_length / 3

    gcp1 = gc1 / part_length
    gcp2 = gc2 / part_length
    gcp3 = gc3 / part_length

    return GCFractions(tot_gc, gcp1, gcp2, gcp3)


class GenePredictor:
    """
    Generic base class for a gene predictor.
    """

    def find_genes(self, work_unit: WorkUnit) -> List[GeneData]:
        raise NotImplementedError


class AdvancedPredictor(GenePredictor):
    """
    Advanced predictor, not implemented yet
    """
    pass


class NaivePredictor(GenePredictor):
    """
    Naive predictor, greedy, and takes the first match.
    Also slow.
    """

    def _find_start_codon(self, dna, start, stop):
        match = start_regexp.search(dna, start, stop)
        if match:
            return match.start(0)
        else:
            return -1

    def _find_stop_codon(self, dna, start, stop):
        # match = stop_regexp.match(dna, start, stop)

        if dna[start:stop] in stop_codons:
            return True
        else:
            return None

    def _find_shine_box(self, scaffold, start_pos,
                        end_pos, shine_box_distance):
        if start_pos - shine_box_distance < 0:
            return None

        match = shine_regexp.search(
            scaffold, start_pos - shine_box_distance, start_pos - 6)
        if match:
            return (match.group(0), start_pos - match.end(0))

    def _probable_gene(self, scaffold, start_pos, end_pos,
                       gc_content, total_gc, has_shine_box):
        allowedDiff = 0.12
        gc3 = total_gc - \
            allowedDiff <= gc_content[3] and total_gc + \
            allowedDiff >= gc_content[3]
        totGc = total_gc - \
            allowedDiff <= gc_content[0] and total_gc + \
            allowedDiff >= gc_content[0]

        if gc3 and totGc and has_shine_box:
            return True
        else:
            return False

    def _find_genes_on_strand(self, work_unit: WorkUnit,
                              strand: int) -> List[GeneData]:

        (scaffold_name, scaffold_sequence,
         total_gc, config) = work_unit
        if not (scaffold_name and scaffold_sequence):
            return

        if strand == -1:
            scaffold_sequence = reverse_complement(scaffold_sequence)

        currentPos = 0
        scaffold_length = len(scaffold_sequence)
        currentGeneNo = 1
        gene_data = []

        while scaffold_length - currentPos >= config.shortest_gene:
            currentPos = self._find_start_codon(
                scaffold_sequence, currentPos, scaffold_length)
            if currentPos != -1:
                posInGene = currentPos + 3
                stop_found = False
                possible_stops = stop_regexp.finditer(
                    scaffold_sequence[posInGene:])

                for stop_match in possible_stops:
                    if stop_match.start(0) % 3 == 0:
                        posInGene += stop_match.start(0)
                        stop_found = True
                        break

                if posInGene + 3 - currentPos >= config.shortest_gene:
                    if stop_found or config.allow_runoff:
                        gcContent = get_gc_content(
                            scaffold_sequence[currentPos:posInGene + 3])
                        shineBox = self._find_shine_box(scaffold_sequence, currentPos,
                                                        posInGene + 3, config.shine_box_distance)
                        hasShineBox = True if shineBox else False
                        if self._probable_gene(scaffold_sequence, currentPos,
                                               posInGene + 3, gcContent, total_gc,
                                               hasShineBox):
                            gene_data.append(GeneData(currentGeneNo,
                                                      currentPos, posInGene + 3, gcContent,
                                                      shineBox, strand))
                            currentGeneNo += 1
                            currentPos = posInGene + 3 + config.intra_gene_gap
                            continue
                        else:
                            currentPos += 1
                            continue

                currentPos += 1

            else:
                break

        return gene_data

    def find_genes(self, work_unit: WorkUnit) -> List[GeneData]:

        gene_data_list = self._find_genes_on_strand(work_unit, 1)
        gene_data_list.extend(self._find_genes_on_strand(work_unit, -1))

        return gene_data_list


class OutputFormatter:
    """
    Base class for an output formatter
    """

    def format_output(
            self, data: Dict[str, ScaffoldData], gc_percentage: int) -> List[str]:
        raise NotImplementedError


class ProdigalFormatter(OutputFormatter):
    """
    Formats output in the prodigal format, more or less
    """

    def _output_gene(self, gene_data: GeneData,
                     scaffold_name: str, sequence: str) -> List[str]:

        (current_gene_no, start_pos,
            end_pos, gc_content, shine_box, strand) = gene_data

        header_string = ""

        direction = "fwd" if strand == 1 else "rev"

        header_string = '>{0}_{6}_{1} # {2} # {3} # {4} # start_type={5};'.format(
            scaffold_name, current_gene_no, start_pos + 1, end_pos,
            strand, sequence[:3], direction)
        if shine_box:
            header_string += 'rbs_motif={0};rbs_spacer={1}bp;'.format(
                shine_box[0], shine_box[1])

        header_string += 'GC1={0:.2f};GC2={1:.2f};GC3={2:.2f};gc_cont={3:.2f}'.format(
            gc_content[1], gc_content[2], gc_content[3], gc_content[0])

        assert(len(sequence) % 3 == 0)

        output_lines = [header_string]

        output_lines.extend(textwrap.wrap(sequence))

        return output_lines

    def format_output(
            self, data: Dict[str, ScaffoldData], gc_percentage: int) -> List[str]:

        output = ["### Total GC content of genome is {0:.2f}% ###".format(gc_percentage)]

        for scaffold_name, contents in data.items():
            for gene_data in contents.genes_data:
                sequence = contents.scaffold_sequence[
                    gene_data.start_pos: gene_data.end_pos]
                output.extend(
                    self._output_gene(
                        gene_data,
                        scaffold_name,
                        sequence))

        return output


class GFF3Formatter(OutputFormatter):
    """
    Formats gff 3 output. 
    See http://www.sequenceontology.org/gff3.shtml
    """

    def _format_gene_gff3(self, gene_data: GeneData, scaffold_name: str, scaffold_length: int,
                          gff3_number: str) -> List[str]:

        (current_gene_no, start_pos,
            end_pos, gc_content, shine_box, strand) = gene_data

        formatted_output = []

        direction = "+"

        # set the correct coordinates for gff3 if we're on the reverse strand

        length = end_pos - start_pos

        assert(length % 3 == 0)
        if strand == -1:
            direction = "-"

            start_pos = scaffold_length - end_pos
            end_pos = start_pos + length

        start_pos += 1  # 1-based indexing

        gff3_line = scaffold_name + "\tGenePredictor\t{0}\t" + str(
            start_pos) + "\t" + str(end_pos) + "\t.\t" + direction + "\t.\t{1}"

        gene = gff3_line.format("gene", "ID=gene" + gff3_number)
        mrna = gff3_line.format(
            "mRNA",
            "ID=tran" +
            gff3_number +
            ";Parent=gene" +
            gff3_number)

        exon = gff3_line.format("exon", "Parent=tran" + gff3_number)
        cds = gff3_line.format("CDS", "Parent=tran" + gff3_number)

        formatted_output.extend([gene, mrna, exon, cds])

        return formatted_output

    def format_output(
            self, data: Dict[str, ScaffoldData], gc_percentage: int) -> List[str]:

        formatted_output = ["##gff-version 3.2.1", "# Total GC content of genome is {0:.2f}%".format(gc_percentage)]
        for scaffold_name, contents in data.items():
            sequence = contents.scaffold_sequence

            formatted_output.append(
                "##sequence-region\t{0}\t{1}\t{2}".format(scaffold_name, 1, len(sequence)))

            for index, gene in enumerate(contents.genes_data):
                gene_output = self._format_gene_gff3(
                    gene, scaffold_name, len(sequence), str(index + 1))
                formatted_output.extend(gene_output)

        return formatted_output


def handleInput(file_obj, config: GpConfig):

    tf = tempfile.SpooledTemporaryFile(mode='w+t')

    genome_length = 0
    genome_gc_count = 0
    for line in file_obj:
        tf.write(line)
        stripped = line.strip().upper()
        if (stripped[0] != '>'):
            for char in stripped:
                if char == 'C' or char == 'G':
                    genome_gc_count += 1

                genome_length += 1

    genome_gc_fraction = genome_gc_count / genome_length

    output_data = collections.OrderedDict()

    scaffold_sequence = ""
    scaffold_name = None

    predictor = NaivePredictor()

    tf.seek(0)
    line = tf.readline()
    while line:

        if line[0] == '>':
            scaffold_name = line[1:-1].strip()

            line = tf.readline()
            scaffold_sequence = ""

            while line and line[0] != '>':

                stripped = line.strip().lower()
                scaffold_sequence = scaffold_sequence + stripped
                line = tf.readline()

            output_data[scaffold_name] = ScaffoldData(scaffold_sequence, [])

            output_data[scaffold_name].genes_data.extend(predictor.find_genes(WorkUnit(scaffold_name,
                                                                                       scaffold_sequence, genome_gc_fraction, config)))

    return (genome_gc_fraction, output_data)

if __name__ == '__main__':
    # argparse

    usage = """
    Finds possible genes in prokaryotic genomes
    """

    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-s', '--shortest_gene',
        help='The shortest possible genes to consider',
        metavar='SHORTEST_GENE',
        dest='shortest_gene',
        type=int,
        default=50
    )

    parser.add_argument(
        '-g', '--intra_gene_gap',
        help='The minimum required gap between genes.',
        metavar='GENE_GAP',
        dest='intra_gene_gap',
        type=int,
        default=40
    )
    parser.add_argument(
        '-b', '--shine_box_distance',
        help='The maximum distance between the Shine-Delgarno sequence (box) and the predicted gene',
        metavar='SHINE_DISTANCE',
        dest='shine_box_distance',
        type=int,
        default=14
    )

    parser.add_argument(
        '--allow_runoff',
        dest='allow_runoff',
        action='store_true',
        help='Allow genes to run off a scaffold')

    parser.add_argument(
        '--prodigal_output',
        dest='prodigal_output',
        action='store_true',
        help='Produce output in semi-prodigal format. (If not set, GFF 3 is used)')

    parser.set_defaults(allow_runoff=False)
    parser.set_defaults(prodigal_output=False)

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help='Input file in FASTA format')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file')

    args = parser.parse_args()

    config = GpConfig(
        args.shortest_gene,
        args.allow_runoff,
        args.intra_gene_gap,
        args.shine_box_distance)

    # redirect stdout
    sys.stdout = args.outfile

    # calculate total genome GC content
    gc_percentage, data = handleInput(args.infile, config)

    output_formatter = None

    if args.prodigal_output:

        output_formatter = ProdigalFormatter()

    else: 

        output_formatter = GFF3Formatter()

    formatted_output = output_formatter.format_output(data, gc_percentage)

    for line in formatted_output:
        print(line)
