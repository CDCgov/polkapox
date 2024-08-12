#!/usr/bin/env python3

### Written by S.Morrison, K. Knipe 20221018
### Refactored by Kyle O'Connell with help from GPT 4o and Gemini 1.5 Pro 20240730

import gzip
import argparse
import os
import sys
import fnmatch
import subprocess
import re
import gfapy
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import json
import shutil

def clean_graph_tags(input_file, output_file):
    # Define the regex pattern to match the LB and CL tags
    tag_pattern = re.compile(r'\tLB:z:[^\t]+\tCL:z:[^\t]+')

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Remove the tags from the line
            cleaned_line = re.sub(tag_pattern, '\n', line)
            # Remove trailing tabs
            cleaned_line = cleaned_line.rstrip('\t')
            # Write the cleaned line to the output file, ensuring the newline character is preserved
            outfile.write(cleaned_line)

def read_gfa_file(gfa_path):
    """Read GFA file into gfapy python structure."""
    try:
        gfa_graph = gfapy.Gfa.from_file(gfa_path)
        # print(gfa_graph.edges)
        if len(gfa_graph.edges) == 0:
            return None, "WARNING: GFA file only contains segment lines"
        return gfa_graph, "PASS"
    except Exception as e:
        
        return None, f"WARNING: Issue with GFA file : {str(e)}"

def remove_self_loops(links):
    """Remove self-loops from link information."""
    pruned_links = [link for link in links if link.from_name != link.to_name]

    return pruned_links, "PASS"

def identify_low_coverage_contigs(segments):
    """Identify contigs below 0.5x coverage."""
    low_cov = []
    for seg in segments:
        dp_value = seg.get('dp')  # This will return None if 'dp' is not a field in the segment
        if dp_value is not None and float(dp_value) < 0.5:
            low_cov.append(seg.get('name'))

    return low_cov, "PASS"

def filter_links_by_coverage(low_cov, links):
    """Filter out links involving low coverage contigs."""
    filtered_links = [link for link in links if link.from_name not in low_cov and link.to_name not in low_cov]
    
    return filtered_links, "PASS"

def find_longest_contig(segments, output_dir):
    """Find the longest contig based on length."""
    longest_contig = max(segments, key=lambda x: x.get('LN'))
    # write longest contig to file
    longest_contig_file = os.path.join(output_dir, "longest_contig.fasta")
    with open(longest_contig_file, "w") as gfaContigs:
        seqName = longest_contig.get('name')
        gfaContigs.write(">" + seqName + "\n")
        gfaContigs.write(longest_contig.get('sequence') + "\n")
    
    return longest_contig.get('name'), "PASS", longest_contig_file

def filter_segments_by_graph(segments, graph):
    """Filter segments to include only those that are part of the graph."""
    graph_nodes = set(graph.nodes())  # Get all nodes in the graph as a set for quick lookup
    filtered_segments = [seg for seg in segments if seg.get('name') in graph_nodes]
    
    return filtered_segments

def create_filtered_graph(links):
    """Create a filtered NetworkX graph from links and remove disconnected nodes."""
    graph = nx.Graph()  # Using a simple Graph instead of MultiGraph if you don't handle multiple edges between the same nodes
    for link in links:
        graph.add_edge(link.from_name, link.to_name)
    
    # Remove disconnected nodes
    disconnected_nodes = [node for node in graph.nodes if graph.degree(node) == 0]
    graph.remove_nodes_from(disconnected_nodes)
    
    if not nx.is_connected(graph):
        return None, "WARNING: Graph is not connected"
    
    return graph, "PASS"

def identify_itr( gfa_graph, segments ):
    """Identify all Inverted Terminal Repeats (ITRs) based on depth criteria and check if they are connected in the graph."""
    depth_data = {seg.get('name'): float(seg.get('dp')) for seg in segments if 'dp' in seg.tagnames}
    lower_bound = 1.5
    upper_bound = 20
    itrs = [contig for contig, depth in depth_data.items() if lower_bound < depth < upper_bound]
    
    itr_length = 0
    visited = set()

    for seg in segments:
        if seg.name in itrs and seg.name not in visited:
            visited.add(seg.name)  # Mark this segment as visited
            itr_length += len(seg.sequence)  # Accumulate the sequence length

    # print("Total ITR sequence length:", itr_length)
    return itrs, itr_length

def get_final_path(gfa_graph, filtered_graph, segments):
    """Find all longest paths in the graph starting from any ITR."""
    itrs, itr_length = identify_itr(gfa_graph, segments)
    if not itrs:
        return [], "WARNING: No suitable ITRs found based on depth criteria."
    
    longest_paths = []
    max_length = 0
    for itr in itrs:
        for contig in filtered_graph.nodes():
            if contig != itr:
                for path in nx.all_simple_paths(filtered_graph, source=itr, target=contig):
                    path_length = len(path)
                    if path_length > max_length:
                        longest_paths = [path]  # Start a new list with the new longest path
                        max_length = path_length
                    elif path_length == max_length:
                        longest_paths.append(path)  # Add path to the list of longest paths
    itr_order = []
    if longest_paths:
        final_path = longest_paths[0]
        itr_order.extend(c for c in final_path if c in itrs)
        final_path.extend(itr_order[::-1])
        return final_path, f"PASS: Found {len(longest_paths)} longest path(s) of length {max_length}"
    else:
        return [], "WARNING: No path found starting from any ITR."

def orient_longest_contig(query, reference, blast_db_dir):
    blastResults = query + "_blast.out"
    blast_db_path = os.path.join(blast_db_dir, query + "_DB")
    makeblastdb_command = ["makeblastdb", "-in", query, "-out", blast_db_path, "-dbtype", "nucl"]
    blastn_command = ["blastn", "-query", reference, "-db", blast_db_path, "-evalue", ".00001", "-outfmt", "6 qseqid sseqid pident length qcovs sstart send sseq mismatch gapopen qlen slen bitscore", "-out", blastResults]
    
    print(f"Running command: {' '.join(makeblastdb_command)}")
    subprocess.call(makeblastdb_command, shell=False)
    
    print(f"Running command: {' '.join(blastn_command)}")
    subprocess.call(blastn_command, shell=False)
    
    # Check if blastResults file is created and contains results
    if not os.path.exists(blastResults):
        print(f"BLAST results file not created: {blastResults}")
        return None, "FAIL"
    
    with open(blastResults, 'r') as e:
        for line in e:
            f13LInfo = line.split('\t')
            sstart = int(f13LInfo[5])
            send = int(f13LInfo[6])
    orientation = 1 if sstart < send else int(-1)
    return orientation, "PASS"

def get_final_orientation(final_path, lnks, longest_contig, longest_orient):
    contig_orientation = [0] * len(final_path)
    for i in range(final_path.index(longest_contig), len(final_path)):
        if i == final_path.index(longest_contig):
            contig_orientation[i] = longest_orient
        else:
            flag=0
            for link in lnks:
                if (link.from_name == final_path[i-1] and link.to_name == final_path[i]):
                    from_orient = 1 if link.from_orient == '+' else int(-1)
                    to_orient = 1 if link.to_orient == '+' else int(-1)
                    contig_orientation[i] = contig_orientation[i-1] * from_orient * to_orient
                    flag=1
            if flag == 0:
                for link in lnks:
                    if (link.to_name == final_path[i-1] and link.from_name == final_path[i]):
                        from_orient = 1 if link.from_orient == '+' else int(-1)
                        to_orient = 1 if link.to_orient == '+' else int(-1)
                        contig_orientation[i] = contig_orientation[i-1] * from_orient * to_orient
                        flag=1

    for i in range(final_path.index(longest_contig)-1, -1 ,-1):
        flag=0
        for link in lnks:
            if (link.to_name == final_path[i+1] and link.from_name == final_path[i]):
                from_orient = 1 if link.from_orient == '+' else int(-1)
                to_orient = 1 if link.to_orient == '+' else int(-1)
                contig_orientation[i] = contig_orientation[i+1] * from_orient * to_orient
                flag=1
        if flag == 0:
            for link in lnks:
                if (link.from_name == final_path[i+1] and link.to_name == final_path[i]):
                    from_orient = 1 if link.from_orient == '+' else int(-1)
                    to_orient = 1 if link.to_orient == '+' else int(-1)
                    contig_orientation[i] = contig_orientation[i+1] * from_orient * to_orient
                    flag=1
        check="PASS"
    # print(final_path)
    # print(contig_orientation)
    return contig_orientation, check

def get_final_sequence(contig_order, contig_orientation, segments):
    segment_info = {}
    for segment in segments:
        segment_count = contig_order.count(segment.get('name'))
        coverage = float(segment.get('dp')) if segment_count == 0 else float(segment.get('dp')) / contig_order.count(segment.get('name'))
        segment_info[segment.get('name')] = {
            'coverage': coverage,
            'sequence': Seq(segment.get('sequence')),
            'sequence_rc': Seq(segment.get('sequence')).reverse_complement()
        }
        
    final_sequence = ''
    final_order_orientation_copy_number = []
    for i in range(len(contig_order)):
        segment_name = contig_order[i].strip('+').strip('-')  # Strip orientation symbols
        # print(f"Processing segment: {segment_name}")
        sequence = segment_info[segment_name]['sequence'] * round(segment_info[segment_name]['coverage'])
        if contig_orientation[i] == -1:
            sequence = segment_info[segment_name]['sequence_rc'] * round(segment_info[segment_name]['coverage'])
        orientation = '+' if contig_orientation[i] == 1 else '-'
        order_orientation_copy_number = '%s%s' % (orientation, segment_name)
        if round(segment_info[segment_name]['coverage']) > 1:
            order_orientation_copy_number = '%sx%s' % (order_orientation_copy_number, round(segment_info[segment_name]['coverage']))
        final_sequence = final_sequence + sequence
        final_order_orientation_copy_number.append(order_orientation_copy_number)
    check = "PASS"
    # are all segments in final_sequence?
    cleaned_contig_order = [contig.strip('+-') for contig in contig_order]
    for segment in segment_info:
        if (segment_info[segment]['coverage'] > 0.5) and (segment not in cleaned_contig_order):
            check = 'WARNING: missing segments'
    return final_sequence, len(final_sequence), " ".join(final_order_orientation_copy_number), check

def write_oriented_fasta(final_path, segments, output_file, input_file):
    """Write the segments in the specified orientation to a single FASTA entry with the input file name as the header."""
    segment_dict = {seg.name: seg.sequence for seg in segments}  # Create a dictionary of segment names to sequences
    full_sequence = ""  # Initialize an empty string to accumulate full sequence

    for segment_name in final_path:
        seq = segment_dict[segment_name.strip('+').strip('-')]  # Remove orientation symbols for lookup
        # Check if the segment should be reversed
        if segment_name.endswith('-'):
            seq = str(Seq(seq).reverse_complement())  # Reverse complement the sequence if needed

        full_sequence += seq  # Concatenate sequence

    # Create a single SeqRecord with all concatenated sequences
    # Use the base name of the input file as the ID for the SeqRecord
    input_base_name = os.path.splitext(os.path.basename(input_file))[0]
    record = SeqRecord(Seq(full_sequence), id=input_base_name, description="All contigs concatenated based on the final path")
    
    # Write the single record to the output file
    with open(output_file, 'w') as f:
        SeqIO.write(record, f, "fasta")

def write_all_contigs(segments, output_file):
    """Write all contigs to a FASTA file without considering orientation."""
    with open(output_file, 'w') as f:
        for seg in segments:
            record = SeqRecord(Seq(seg.sequence), id=seg.name, description="")
            SeqIO.write(record, f, "fasta")

def write_log_and_exit(log, status):
    # print(log)
    log_file = os.path.join(log['00']['input']['output_dir'], log['00']['input']['sample_name'] + ".log")
    summary_file = os.path.join(log['00']['input']['output_dir'], log['00']['input']['sample_name'] + ".summary")

    headers = ['sample', 'status', 'contig_order', 'contig_orientation', 'contig_order_orientation_copy_number', 'assembly_length', 'itr_length']
    data = []
    if '00' in log:
        data.append(log['00']['input']['sample_name'])
        data.append(status)
    if '03' in log:
        data.append(",".join(str(i) for i in log['03']['output']['final_path']))
    if '06' in log:
        data.append(",".join(str(i) for i in log['06']['output']['final_orientation']))
    if '07' in log:
        data.append(str(log['07']['output']['final_order_orientation_copy_number']))
        data.append(str(log['07']['output']['final_sequence_length']))
    if '09' in log:
        data.append(str(log['09']['output']['final_itr_length']))

    with open(summary_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join(data) + '\n')

    with open(log_file, 'w') as f:
        json_object = json.dumps(log, indent=4)
        f.write(json_object)
    sys.exit(0)

def process_graph(gfa_graph, output_dir, input_file, reference):
    """Process the graph and write output."""
    log = {}
    input_with_ext = os.path.basename(input_file)
    input_base, _ = os.path.splitext(input_with_ext)
    sample_name = input_base

    # Initial log entry
    log['00'] = {'step_name': "initialization",
                 'step_description': "Initialize the logging process",
                 'status': "PASS",
                 'input': {
                     'gfa': input_file,
                     'sample_name': sample_name,
                     'output_dir': output_dir
                 },
                 'output': {}}

    # Always write all contigs to a FASTA file at the start
    write_all_contigs(gfa_graph.segments, os.path.join(output_dir, input_base + ".all_contigs.fasta"))

    # Start processing the graph
    filtered_edges, status = remove_self_loops(gfa_graph.edges)
    log['01'] = {'step_name': "remove_self_loops",
                 'step_description': "Remove self-loops from link information",
                 'status': status,
                 'output': {'filtered_edges': [str(edge) for edge in filtered_edges]}}
    if status != "PASS":
        write_log_and_exit(log, status)

    filtered_graph, status = create_filtered_graph(filtered_edges)
    log['02'] = {'step_name': "create_filtered_graph",
                 'step_description': "Create a filtered NetworkX graph from links and remove disconnected nodes",
                 'status': status,
                 'output': {'filtered_graph': str(filtered_graph)}}
    if status != "PASS":
        write_log_and_exit(log, status)

    filtered_segments = filter_segments_by_graph(gfa_graph.segments, filtered_graph)
    final_path, status = get_final_path(gfa_graph, filtered_graph, filtered_segments)
    log['03'] = {'step_name': "get_final_path",
                 'step_description': "Find all longest paths in the graph starting from any ITR",
                 'status': status,
                 'output': {'final_path': final_path}}
    if not status.startswith("PASS"):
        write_log_and_exit(log, status)

    # Find the longest contig
    longest_contig, status, longest_contig_file = find_longest_contig(filtered_segments, output_dir)
    log['04'] = {'step_name': "find_longest_contig",
                 'step_description': "Find the longest contig based on length",
                 'status': status,
                 'output': {'longest_contig': longest_contig}}
    if status != "PASS":
        write_log_and_exit(log, status)

    # Create blast database directory
    blast_db_dir = os.path.join(output_dir, 'blast_db')
    os.makedirs(blast_db_dir, exist_ok=True)  # Create if it doesn't exist

    # Determine the orientation of the longest contig
    longest_orient, status = orient_longest_contig(longest_contig_file, reference, blast_db_dir)
    log['05'] = {'step_name': "orient_longest_contig",
                 'step_description': "Determine the orientation of the longest contig",
                 'status': status,
                 'output': {'longest_orient': longest_orient}}
    if status != "PASS":
        write_log_and_exit(log, status)

    # Get the final orientation for all contigs
    contig_orientation, status = get_final_orientation(final_path, filtered_edges, longest_contig, longest_orient)
    log['06'] = {'step_name': "get_final_orientation",
                 'step_description': "Get the final orientation for all contigs",
                 'status': status,
                 'output': {'final_orientation': contig_orientation}}
    if status != "PASS":
        write_log_and_exit(log, status)

    # Re-orient the final path based on determined orientation
    oriented_final_path = [contig + ('-' if contig_orientation[i] == -1 else '+') for i, contig in enumerate(final_path)]

    # Get the final sequence and order
    final_sequence, seq_len, final_order_orientation_copy_number, status = get_final_sequence(oriented_final_path, contig_orientation, filtered_segments)
    log['07'] = {'step_name': "get_final_sequence",
                 'step_description': "Get the final sequence and order",
                 'status': status,
                 'output': {'final_sequence': str(final_sequence),  # Convert Seq object to string
                            'final_sequence_length': seq_len,
                            'final_order_orientation_copy_number': final_order_orientation_copy_number}}
    if status != "PASS":
        write_log_and_exit(log, status)

    # Write the oriented contigs to a FASTA file
    asm_name = input_base + ".assembly_asm.fasta"
    record = SeqRecord(Seq(final_sequence), id=input_base, description=final_order_orientation_copy_number)
    with open(os.path.join(output_dir, asm_name), 'w') as f:
        SeqIO.write(record, f, "fasta")

    log['08'] = {'step_name': "write_oriented_fasta",
                 'step_description': "Write the oriented contigs to a FASTA file",
                 'status': "PASS",
                 'output': {'assembly_file': asm_name}}

    # Remove the blast_db directory
    shutil.rmtree(blast_db_dir)  # Remove the directory and its contents

    # Identify ITRs and their length
    itrs, itr_length = identify_itr(gfa_graph, filtered_segments)
    log['09'] = {'step_name': "identify_itr",
                 'step_description': "Identify all Inverted Terminal Repeats (ITRs) based on depth criteria",
                 'status': "PASS",
                 'output': {'itrs': itrs,
                            'final_itr_length': itr_length}}

    # Write the final log and summary
    write_log_and_exit(log, "PASS")

def main(arguments):
    parser = argparse.ArgumentParser(description="GFA parser to construct assembly from Unicycler")
    parser.add_argument('-i', type=str, required=True, help="GFA file from Unicycler")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-r', type=str, required=True, help="Reference file for orientation")
    args = parser.parse_args(arguments)

    gfa_file = args.i
    output_dir = args.o
    reference_file = args.r

    # clean up gfa tags
    cleanedGfa = 'graph_cleaned.gfa'
    clean_graph_tags(gfa_file, cleanedGfa)

    gfa_graph, status = read_gfa_file(cleanedGfa)
    if status != "PASS":
        return

    if gfa_graph is None:
        print("Failed to read GFA file.")
        return

    process_graph(gfa_graph, output_dir, gfa_file, reference_file)

if __name__ == '__main__':
    main(sys.argv[1:])