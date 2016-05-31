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

__author__ = 'uvesten'

start_codons = ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
stop_codons = ['TAA', 'TAG', 'TGA']

start_regexp = re.compile('AT.|[ATCG]TG')
stop_regexp = re.compile('TA[GA]|TGA')


gc_translation = str.maketrans("CG", "11", "ATN")

# named tuple to hold configuration data
GpConfig = collections.namedtuple(
    'GpConfig', [
        'shortest_gene', 'allow_runoff', 'intra_gene_gap', 'shine_box_distance'])


def get_gc_content(dna):

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

    return (tot_gc, gcp1, gcp2, gcp3)


def find_start_codon(dna, start, stop):
    match = start_regexp.search(dna, start, stop)
    if match:
        return match.start(0)
    else:
        return -1


def find_stop_codon(dna, start, stop):
    #match = stop_regexp.match(dna, start, stop)

    if dna[start:stop] in stop_codons:
        return True
    else:
        return None


def output_gene(name, count, scaffold, start_pos,
                end_pos, gc_content, shine_box, strand):
    print('>{0}_{1} # {2} # {3} # {4} # start_type={5};'.format(
        name, count, start_pos + 1, end_pos,
        strand, scaffold[start_pos:start_pos + 3]), end="")
    if shine_box:
        print('rbs_motif={0};rbs_spacer={1}bp;'.format(
            shine_box[0], shine_box[1]), end="")

    print('GC1={0:.2f};GC2={1:.2f};GC3={2:.2f};gc_cont={3:.2f}'.format(
        gc_content[1], gc_content[2], gc_content[3], gc_content[0]))

    assert(len(scaffold[start_pos:end_pos]) % 3 == 0)
    for line in textwrap.wrap(scaffold[start_pos:end_pos]):
        print(line)

shine_regexp = re.compile('A?G?GAGG|GGAG|GG.{1}GG')


def find_shine_box(scaffold, start_pos, end_pos, shine_box_distance):
    if start_pos - shine_box_distance < 0:
        return None

    match = shine_regexp.search(
        scaffold, start_pos - shine_box_distance, start_pos - 6)
    if match:
        return (match.group(0), start_pos - match.end(0))


def probable_gene(scaffold, start_pos, end_pos,
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


def find_genes(name, scaffold_no, scaffold,
               total_gc, strand, config: GpConfig):
    if not (name and scaffold):
        return

    currentPos = 0
    scaffoldLength = len(scaffold)
    currentGeneNo = 1
    gene_data = []

    while scaffoldLength - currentPos >= config.shortest_gene:
        currentPos = find_start_codon(scaffold, currentPos, scaffoldLength)
        if currentPos != -1:
            posInGene = currentPos + 3
            stop_found = False
            possible_stops = stop_regexp.finditer(scaffold[posInGene:])

            for stop_match in possible_stops:
                if stop_match.start(0) % 3 == 0:
                    posInGene += stop_match.start(0)
                    stop_found = True
                    break

            if posInGene + 3 - currentPos >= config.shortest_gene:
                if stop_found or config.allow_runoff:
                    gcContent = get_gc_content(
                        scaffold[currentPos:posInGene + 3])
                    shineBox = find_shine_box(scaffold, currentPos,
                                              posInGene + 3, config.shine_box_distance)
                    hasShineBox = True if shineBox else False
                    if probable_gene(scaffold, currentPos,
                                     posInGene + 3, gcContent, total_gc,
                                     hasShineBox):
                        gene_data.append((name, currentGeneNo, scaffold,
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


complement = str.maketrans("ACTG", "TGAC")


def reverse_complement(dna):

    dna_comp = dna.translate(complement)
    return dna_comp[::-1]


def worker(work_queue, results_queue):
    while True:
        try:
            (current_name, current_scaffold_number,
             current_scaffold, tot_gc_pct, strand, config) = work_queue.get_nowait()

            gene_data = find_genes(current_name, current_scaffold_number,
                                   current_scaffold, tot_gc_pct, strand, config)

            results_queue.put(gene_data)

            work_queue.task_done()

        except queue.Empty:
            break


def handleInput(file_obj, config, called_from_cmdline):

    tf = tempfile.SpooledTemporaryFile(mode='w+t')

    totCount = 0
    gcCount = 0
    for line in file_obj:
        tf.write(line)
        stripped = line.strip().upper()
        if (stripped[0] != '>'):
            for char in stripped:
                if char == 'C' or char == 'G':
                    gcCount += 1

                totCount += 1

    totGcPct = gcCount / totCount

    print("### Total GC content of genome is {0:.2f}% ###".format(totGcPct))

    # then set up queues for multiprocessing

    work_queue = mp.JoinableQueue()
    results_queue = mp.Queue()

    # start as many workers as we have cpu cores -1 (for the main process.)

    max_children = mp.cpu_count() - 1

    processes = [
        mp.Process(
            target=worker,
            args=(
                work_queue,
                results_queue)) for i in range(max_children)]

    current_scaffold = ""
    current_scaffold_number = 0
    current_name = None

    job_counter = 0

    tf.seek(0)
    for line in tf:
        if line[0] == '>':
            if current_scaffold != "":
                work_queue.put((current_name, current_scaffold_number,
                                current_scaffold, totGcPct, 1, config))

                work_queue.put((current_name, current_scaffold_number,
                                reverse_complement(current_scaffold), totGcPct, -1, config))

                job_counter += 2
            current_name = line[1:-1].strip()
            current_scaffold = ""
            current_scaffold_number += 1
        else:

            stripped = line.strip().upper()
            current_scaffold = current_scaffold + stripped

    work_queue.put((
        current_name,
        current_scaffold_number,
        current_scaffold,
        totGcPct,
        1,
        config))

    work_queue.put((current_name, current_scaffold_number,
                    reverse_complement(current_scaffold), totGcPct, -1, config))

    job_counter += 2

    for proc in processes:
        proc.start()

    web_data = []
    while job_counter > 0:
        gene_data = results_queue.get()
        job_counter -= 1
        if called_from_cmdline:
            for item in gene_data:
                output_gene(*item)
        else:
            web_data.extend(gene_data)

    if not called_from_cmdline:
        return (totGcPct, web_data)


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

    parser.set_defaults(allow_runoff=False)

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help='Input file in FASTA format')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file, in semi-Prodigal format')

    args = parser.parse_args()

    config = GpConfig(
        args.shortest_gene,
        args.allow_runoff,
        args.intra_gene_gap,
        args.shine_box_distance)

    # redirect stdout
    sys.stdout = args.outfile

    # calculate total genome GC content
    handleInput(args.infile, config, True)
