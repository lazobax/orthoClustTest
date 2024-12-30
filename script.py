#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import logging
import tempfile
import shutil
from collections import Counter
import datetime
from itertools import combinations

def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info('Script started.')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Protein clustering script using MMseqs2 and MCL (optionally DIAMOND).')
    parser.add_argument('input_directory', help='Directory containing FASTA proteome files.')
    parser.add_argument('prefix', help='Prefix for output files and directory.')

    # Threads
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of CPU threads to use for the aligner. Default is 1.'
    )

    # Aligner choice: 'm' = MMseqs2, 'd' = DIAMOND
    parser.add_argument(
        '--aligner', '-a',
        choices=['m', 'd'],
        default='m',
        help='Choose the aligner: mmseqs2 (m) or DIAMOND (d). Default is mmseqs2.'
    )

    # Debug mode
    parser.add_argument(
        '--debug',
        action='store_true',
        help='If set, keep temporary files upon error for debugging.'
    )

    return parser.parse_args()

def run_mmseqs2(input_dir, tmp_dir, threads):
    logging.info('Running MMseqs2 all-vs-all search.')
    db_dir = os.path.join(tmp_dir, 'db')
    result_dir = os.path.join(tmp_dir, 'result')

    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    # Concatenate all FASTA files
    concatenated_fasta = os.path.join(tmp_dir, 'all_proteins.fasta')
    with open(concatenated_fasta, 'w') as outfile:
        for filename in os.listdir(input_dir):
            if filename.endswith('.fasta') or filename.endswith('.fa'):
                with open(os.path.join(input_dir, filename), 'r') as infile:
                    shutil.copyfileobj(infile, outfile)

    # Create MMseqs2 database
    cmd = [
        'mmseqs', 'createdb',
        concatenated_fasta,
        os.path.join(db_dir, 'proteins')
    ]
    subprocess.run(cmd, check=True)

    # Run MMseqs2 all-vs-all search
    search_result = os.path.join(result_dir, 'search_result')
    cmd = [
        'mmseqs', 'search',
        os.path.join(db_dir, 'proteins'),
        os.path.join(db_dir, 'proteins'),
        search_result,
        tmp_dir,
        '-a',                  # Include alignment
        '--threads', str(threads)  # Number of threads
    ]
    subprocess.run(cmd, check=True)

    # Convert results to TSV (m8) format
    m8_result = os.path.join(result_dir, 'result.m8')
    cmd = [
        'mmseqs', 'convertalis',
        os.path.join(db_dir, 'proteins'),
        os.path.join(db_dir, 'proteins'),
        search_result,
        m8_result,
        '--format-output', 'query,target,pident',
        '-t', str(threads)
    ]
    subprocess.run(cmd, check=True)

    logging.info('MMseqs2 search completed.')
    return m8_result

def run_diamond(input_dir, tmp_dir, threads):
    logging.info('Running DIAMOND all-vs-all search.')
    result_dir = os.path.join(tmp_dir, 'result')
    os.makedirs(result_dir, exist_ok=True)

    # Concatenate all FASTA files
    concatenated_fasta = os.path.join(tmp_dir, 'all_proteins.faa')
    with open(concatenated_fasta, 'w') as outfile:
        for filename in os.listdir(input_dir):
            if filename.endswith('.fasta') or filename.endswith('.fa'):
                with open(os.path.join(input_dir, filename), 'r') as infile:
                    shutil.copyfileobj(infile, outfile)

    # Make DIAMOND DB
    diamond_db = os.path.join(tmp_dir, 'diamond_db')
    cmd = [
        'diamond', 'makedb',
        '--in', concatenated_fasta,
        '--db', diamond_db,
        '--ignore-warnings'
    ]
    subprocess.run(cmd, check=True)

    # Run DIAMOND all-vs-all
    m8_result = os.path.join(result_dir, 'result.m8')
    cmd = [
        'diamond', 'blastp',
        '--db', diamond_db,
        '--query', concatenated_fasta,
        '--out', m8_result,
        '--outfmt', '6 qseqid sseqid pident',  # Format 6
        '--threads', str(threads),
        '--evalue', '1e-5',
        '--ignore-warnings'
    ]
    subprocess.run(cmd, check=True)

    logging.info('DIAMOND search completed.')
    return m8_result

def build_graph(m8_file):
    logging.info('Building graph from alignment results.')
    edges = []
    nodes = set()
    with open(m8_file, 'r') as infile:
        for line in infile:
            query, target, pident = line.strip().split('\t')
            pident = float(pident)
            if pident > 0:
                edges.append((query, target, pident))
                nodes.update([query, target])
    logging.info(f'Graph has {len(nodes)} nodes and {len(edges)} edges.')
    return nodes, edges

def write_mcl_input(edges, mcl_input_file):
    logging.info('Writing MCL input file.')
    with open(mcl_input_file, 'w') as outfile:
        for source, target, weight in edges:
            outfile.write(f'{source}\t{target}\t{weight}\n')

def run_mcl(mcl_input_file, mcl_output_file):
    logging.info('Running MCL clustering.')
    cmd = ['mcl', mcl_input_file, '--abc', '-o', mcl_output_file]
    subprocess.run(cmd, check=True)
    logging.info('MCL clustering completed.')

def parse_fasta_headers(input_dir):
    logging.info('Parsing FASTA headers.')
    headers = {}
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            filepath = os.path.join(input_dir, filename)
            with open(filepath, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        header = line.strip()
                        # Attempt to parse by '|'
                        # Typical format: >db|protein_id|something_else ...
                        parts = header.split('|')
                        if len(parts) >= 3:
                            short_id = parts[1] + '|' + parts[2]
                            # Also store the original protein_id key (the full first token after '>')
                            full_id = header.split()[0][1:]
                            headers[full_id] = short_id
                        else:
                            full_id = header.split()[0][1:]
                            headers[full_id] = full_id
    return headers

def get_top_words(cluster_headers, stopwords):
    words = []
    for header in cluster_headers:
        annotation = ' '.join(header.split()[1:])
        words.extend(annotation.split())
    words = [word for word in words if word.lower() not in stopwords]
    word_counts = Counter(words)
    top_words = [word for word, count in word_counts.most_common(3)]
    return ', '.join(top_words)

def generate_csv(clusters, headers, input_dir, output_dir, prefix):
    logging.info('Generating output CSV files.')
    files = [f for f in os.listdir(input_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    files.sort()
    file_proteins = {f: set() for f in files}
    for f in files:
        filepath = os.path.join(input_dir, f)
        with open(filepath, 'r') as infile:
            for line in infile:
                if line.startswith('>'):
                    protein_id = line.strip().split()[0][1:]
                    file_proteins[f].add(protein_id)

    main_csv = os.path.join(output_dir, f'{prefix}_clusters.csv')
    unassigned_csv = os.path.join(output_dir, f'{prefix}_unassigned_orthogroups.csv')
    cluster_names_file = os.path.join(output_dir, f'{prefix}_cluster_names.txt')
    pairs_file = os.path.join(output_dir, f'{prefix}_ortholog_pairs.tsv')

    stopwords = ['protein', 'domain', 'those', 'hypothetical', 'conserved', 'putative']

    with open(main_csv, 'w') as main_outfile, \
         open(unassigned_csv, 'w') as unassigned_outfile, \
         open(cluster_names_file, 'w') as names_outfile, \
         open(pairs_file, 'w') as pairs_outfile:

        # Write CSV headers
        header_line = 'Orthogroup,' + ','.join(files) + '\n'
        main_outfile.write(header_line)
        unassigned_outfile.write(header_line)

        orthogroup_counter = 1

        for cluster in clusters:
            cluster_proteins = cluster
            cluster_headers = [headers.get(pid, '') for pid in cluster_proteins]

            # Determine which files the proteins come from
            species_present = set()
            proteins_in_species = {f: [] for f in files}
            for protein_id in cluster_proteins:
                for f in files:
                    if protein_id in file_proteins[f]:
                        species_present.add(f)
                        short_id = headers.get(protein_id, protein_id)
                        proteins_in_species[f].append(short_id)

            # Assign orthogroup name
            orthogroup_name = f'OG{orthogroup_counter:06d}'
            orthogroup_counter += 1

            # Get top words for cluster naming
            top_words = get_top_words(cluster_headers, stopwords)

            # Write to cluster names file
            names_outfile.write(f'{orthogroup_name}: {top_words}\n')

            # Build CSV row
            row = [orthogroup_name]
            for f in files:
                row.append(';'.join(proteins_in_species[f]))

            # Decide if assigned or unassigned
            if len(species_present) > 1:
                main_outfile.write(','.join(row) + '\n')
            else:
                unassigned_outfile.write(','.join(row) + '\n')

            # Write all pairwise combos
            short_cluster_proteins = [headers.get(pid, pid) for pid in cluster_proteins]
            for p1, p2 in combinations(short_cluster_proteins, 2):
                pairs_outfile.write(f'{p1}\t{p2}\n')

def main():
    args = parse_arguments()

    # Create output directory
    output_dir = args.prefix
    os.makedirs(output_dir, exist_ok=True)

    # Setup logging
    log_file = os.path.join(output_dir, f'{args.prefix}_protein_clustering.log')
    setup_logging(log_file)

    # Temporary directory
    tmp_dir = tempfile.mkdtemp(prefix='tmp_', dir=output_dir)
    success = False  # Track if the script ran successfully

    try:
        # Decide which aligner to use
        if args.aligner == 'm':
            m8_file = run_mmseqs2(args.input_directory, tmp_dir, args.threads)
        else:
            m8_file = run_diamond(args.input_directory, tmp_dir, args.threads)

        # Build graph
        nodes, edges = build_graph(m8_file)

        # Write MCL input
        mcl_input_file = os.path.join(tmp_dir, 'mcl_input.txt')
        write_mcl_input(edges, mcl_input_file)

        # Run MCL
        mcl_output_file = os.path.join(tmp_dir, 'mcl_output.txt')
        run_mcl(mcl_input_file, mcl_output_file)

        # Parse headers
        headers = parse_fasta_headers(args.input_directory)

        # Read clusters
        with open(mcl_output_file, 'r') as infile:
            clusters = [line.strip().split('\t') for line in infile]

        # Generate CSV
        generate_csv(clusters, headers, args.input_directory, output_dir, args.prefix)

        logging.info('Script completed successfully.')
        logging.info('Script ended.')
        success = True

    except Exception as e:
        logging.error('An error occurred:', exc_info=True)

    finally:
        # If the script did NOT fail or if --debug is NOT set, remove temp files
        # In other words, only keep temp files if there's an error AND --debug is True
        if success or not args.debug:
            shutil.rmtree(tmp_dir)

if __name__ == '__main__':
    main()
