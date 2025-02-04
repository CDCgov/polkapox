#!/usr/bin/env python3

### Written by S.Morrison, K. Knipe 20221018
### Refactored by Kyle O'Connell with help from GPT 4o, o1-preview, and Gemini 1.5 Pro 20240730

import os
import re
import sys
import json
import gfapy
import shutil
import argparse
import subprocess
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        if len(gfa_graph.edges) == 0:
            return None, "WARNING: GFA file only contains segment lines"
        return gfa_graph, "PASS"
    except Exception as e:
        print('Issue with gfa file')
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

def create_filtered_graph(links, segments):
    """Create a filtered NetworkX graph from links and remove disconnected nodes."""
    graph = nx.MultiGraph()

    for link in links:
        graph.add_edge(link.from_name, link.to_name)
    
    # Remove disconnected nodes
    disconnected_nodes = [node for node in graph.nodes if graph.degree(node) == 0]
    graph.remove_nodes_from(disconnected_nodes)
    if not nx.is_connected(graph):
        return None, "WARNING: Graph is not connected"
    
    # Check if the graph forms a cycle
    # get the number of edges for each node
    degrees = dict(graph.degree())
    # figure out how many nodes have more or less than two edges, count 
    count_not_two = sum(1 for degree in degrees.values() if degree < 2)
    all_degree_two = count_not_two <= 1

    num_edges = graph.number_of_edges()
    num_nodes = graph.number_of_nodes()

    # make sure the number of edges and nodes are equal (perfect loop) or one off (loop + ITR)

    if all_degree_two and num_edges == num_nodes or num_edges == num_nodes+1:
        print('The graph forms a complete cycle')
        return graph, "PASS"
    else:
        print("graph is not circular")
        return None, "WARNING: Graph is not circular"

def identify_itr(filtered_edges, segments):
    """
    Identify all Inverted Terminal Repeats (ITRs) based on the orientation of their connections.
    A terminal ITR is identified as a contig that has connections only on the same end (orientation),
    whereas other contigs are connected on both '+' and '-' ends.
    """
    # Create a dictionary to store the set of connected orientations for each contig
    depth_data = {seg.get('name'): float(seg.get('dp')) for seg in segments if 'dp' in seg.tagnames}
    lower_bound = 1.2
    filt_depth = 0.6
    potential_itrs = [contig for contig, depth in depth_data.items() if lower_bound < depth]
    filt_contigs = [contig for contig, depth in depth_data.items() if filt_depth >= depth]
    print("potential ITRs", potential_itrs)

    # Identify the terminal ITR
    contig_ends = []
    for itr in potential_itrs:
        print('itr in loop', itr)
        from_orients = []
        to_orients = []
        for edge in filtered_edges:
            if edge.from_name == itr or edge.to_name == itr:
                print(edge)
            if edge.from_name != edge.to_name:
                if edge.from_name == itr:
                    from_orients.append(edge.from_orient)
                elif edge.to_name == itr:
                    to_orients.append(edge.to_orient)
        # Create set of each list
        from_orients = set(from_orients)
        to_orients = set(to_orients)
        print(itr, from_orients, to_orients)

        # Determine if the contig is a terminal ITR based on orientations
        if (len(from_orients) == 1 and len(to_orients) == 0) or (len(from_orients) == 0 and len(to_orients) == 1):
            contig_ends.append(itr)
        elif len(from_orients) == 1 and len(to_orients) == 1:
            if list(from_orients)[0] != list(to_orients)[0]:
                contig_ends.append(itr)

    print("Terminal ITR contigs based on links:", contig_ends)

    # Create a subgraph with only potential ITRs to check for connected ITR contigs
    subgraph = nx.Graph()
    for edge in filtered_edges:
        if edge.from_name in potential_itrs and edge.to_name in potential_itrs:
            subgraph.add_edge(edge.from_name, edge.to_name)

    # Find connected components in the subgraph
    connected_components = list(nx.connected_components(subgraph))
    print('connected components', connected_components)

    # Select the connected component containing the terminal ITR
    if connected_components:
        for group in connected_components:
            if contig_ends and contig_ends[0] in group:
                itrs = max(connected_components, key=len)
                print('these are the itrs', itrs)
                break  # Exit the loop since we've found the group
        else:
            print('Warning, contig end not in connected components, assigning itr to terminal itr')
            itrs = contig_ends
    else:
        itrs = contig_ends
        print('Warning, no connected components, assigning itr to only terminal itr')

    itr_length = sum(len(seg.sequence) for seg in segments if seg.name in itrs)

    if itrs and len(contig_ends) < 2:
        status = 'PASS: itrs conform to expectation'
    else:
        status = 'FAIL: Issues with ITR identification'

    print("contigs in ITR", itrs, "with total length", itr_length)
    print(status)
    return list(itrs), itr_length, status

def get_final_paths(filtered_edges, filtered_graph, segments):
    """Find all longest paths in the graph starting from any ITR."""
    itrs, itr_length, itr_status = identify_itr(filtered_edges, segments)
    if not itrs or 'FAIL' in itr_status:
        status = itr_status
        return [], [], [], status

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
                        if path not in longest_paths:
                            longest_paths.append(path)  # Add path to the list of longest paths

    final_paths = []
    if longest_paths:
        print("longest path list =", longest_paths)

        # Ensure ITRs are not duplicated and are correctly oriented
        for path in longest_paths:
            left_itrs = [itr for itr in path if itr in itrs]
            right_itrs = [itr for itr in reversed(path) if itr in itrs]

            # Construct the final path with ITRs on both ends, ensuring no duplicates
            middle_path = [node for node in path if node not in itrs]
            final_path = left_itrs + middle_path + right_itrs
            final_paths.append(final_path)
        print("Found Longest Path")
        status = f"PASS: Found {len(longest_paths)} longest path(s) of length {max_length}"
        return final_paths, itrs, itr_length, status
    else:
        status = "FAIL: No path found starting from any ITR."
        return [], [], [], status
        
def orient_longest_contig(query, reference, blast_db_dir):
    blastResults = query + "_blast.out"
    blast_db_path = os.path.join(blast_db_dir, query + "_DB")
    makeblastdb_command = ["makeblastdb", "-in", query, "-out", blast_db_path, "-dbtype", "nucl"]
    blastn_command = ["blastn", "-query", reference, "-db", blast_db_path, "-evalue", "1e-10", "-word_size", "28", "-outfmt", "6 qseqid sseqid pident length qcovs sstart send mismatch gapopen qlen slen bitscore", "-out", blastResults]
    #print(f"Running command: {' '.join(makeblastdb_command)}")
    result = subprocess.run(makeblastdb_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #print(f"Running command: {' '.join(blastn_command)}")
    result2 = subprocess.run(blastn_command, shell=False)
    
    # Check if blastResults file is created and contains results
    if not os.path.exists(blastResults):
        print(f"BLAST results file not created: {blastResults}")
        return None, "FAIL"

    with open(blastResults, 'r') as e:
        for line in e:
            blastHits = line.split('\t')
           #print(blastHits)
            sstart = int(blastHits[5])
            send = int(blastHits[6])
    orientation = '1 +' if sstart < send else '1 -'
    print('orientation of longest contig',orientation)
    return orientation, "PASS"

def get_final_orientation(final_paths, lnks, longest_contig, longest_orient, itrs):
    # Build a dictionary for quick lookup of links
    links_dict = {}
    for link in lnks:
        print(link)
        
        # Map direct links
        key = (link.from_name, link.to_name)
        links_dict[key] = (link.from_orient, link.to_orient)
        # Map reverse links with flipped orientations
        reverse_key = (link.to_name, link.from_name)
        # Flip the orientations for the reverse mapping
        reversed_from_orient = '-' if link.to_orient == '+' else '+'
        reversed_to_orient = '-' if link.from_orient == '+' else '+'
        links_dict[reverse_key] = (reversed_from_orient, reversed_to_orient)

    final_oriented_path = []

    if len(final_paths) > 1:
        # Iterate over each path
        for path in final_paths:
            contig_orientations = []
            print('Processing path:', path)

            if not path:
                continue  # Skip empty paths

            # Assign orientation to the first contig
            first_contig = path[0]
            first_contig_orient = None
            # Find an initial orientation for the first contig
            for link in lnks:
                if first_contig == link.from_name:
                    first_contig_orient = link.from_orient
                    break
                elif first_contig == link.to_name:
                    # Flip the orientation since it's the 'to_name' in the link
                    first_contig_orient = '-' if link.to_orient == '+' else '+'
                    break
            if first_contig_orient is None:
                print(f"Could not determine orientation for contig {first_contig}")
                continue  # Skip this path if orientation is unknown

            contig_orientations.append(f"{first_contig} {first_contig_orient}")
            previous_contig = first_contig
            previous_orient = first_contig_orient

            # Process the rest of the contigs in the path
            for current_contig in path[1:]:
                key = (previous_contig, current_contig)
                if key in links_dict:
                    prev_orient_in_link, curr_orient_in_link = links_dict[key]
                    # Adjust current orientation based on the previous orientation
                    if previous_orient == prev_orient_in_link:
                        current_orient = curr_orient_in_link
                    else:
                        # If previous orientation doesn't match, flip the current orientation
                        current_orient = '-' if curr_orient_in_link == '+' else '+'
                    contig_orientations.append(f"{current_contig} {current_orient}")
                    previous_contig = current_contig
                    previous_orient = current_orient
                else:
                    print(f"No link found between {previous_contig} and {current_contig}")
                    break  # Cannot process further without a link

            print('Final contig orientations:', contig_orientations)
            # Check if this path contains the longest contig with the correct orientation
            for contig_info in contig_orientations:
                contig_name, contig_orient = contig_info.split(' ')
                if contig_name == longest_contig and contig_orient == longest_orient.split(' ')[1]:
                    final_oriented_path = contig_orientations
                    break

        # After all paths have been processed
        print('Final oriented path with correct longest contig orientation:', final_oriented_path)
        return path, final_oriented_path, "PASS"

    else:
        # Only one longest path; apply new logic
        return None, [], "FAIL: Unable to resolve final orientation, likely only one longest contig and ITRs"
        
def get_final_sequence(oriented_path, segments):
    """Generate the final sequence for a given path and orientation."""
    # Create a dictionary mapping segment names to sequences
    segment_dict = {seg.name: seg.sequence for seg in segments}
    final_sequence = ""  # Initialize an empty string to accumulate the final sequence

    for segment_name in oriented_path:
        #print(f"Processing segment: {segment_name}")
        # Split the segment name and orientation
        contig, orient = segment_name.strip().split(' ')
        #print(f"Contig: {contig}, Orientation: {orient}")
        # Get the sequence corresponding to the contig name
        seq = segment_dict.get(contig)
        if seq is None:
            print(f"Warning: Sequence for contig '{contig}' not found.")
            continue  # Skip if the contig is not found

        # Check if the segment should be reversed
        if orient == '-':
            seq = str(Seq(seq).reverse_complement())  # Reverse complement the sequence if needed

        final_sequence += seq  # Concatenate the sequence

    final_sequence_length = len(final_sequence)
    print('final_sequence_length = ', final_sequence_length)
    # Create a comma-separated string of the oriented path
    final_order_orientation_copy_number = ",".join(oriented_path)
    return final_sequence, final_sequence_length, final_order_orientation_copy_number, "PASS"

def write_oriented_fasta(final_path, segments, output_file, input_file):
    """Write the segments in the specified orientation to a single FASTA entry with the input file name as the header."""
    print('writing final oriented fasta to out file')
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
    log_file = os.path.join(log['00']['input']['output_dir'], log['00']['input']['sample_name'] + ".assembly.log")
    summary_file = os.path.join(log['00']['input']['output_dir'], log['00']['input']['sample_name'] + ".assembly.summary")

    headers = ['sample', 'status', 'longest_paths', 'final_contig_order_orientation', 'itr_length', 'assembly_length']
    data = []
    if '00' in log:
        data.append(log['00']['input']['sample_name'])  # 'sample'
        data.append(status)  # 'status'
    if '05' in log:
        data.append(",".join(str(i) for i in log['05']['output']['final_paths']))  # 'longest_paths'
    if '08' in log:
        data.append(str(log['08']['output']['final_path']))  # 'final_path'
        data.append(",".join(str(i) for i in log['08']['output']['final_orientation']))  # 'final_contig_order_orientation'    
    else:
        data.append('')  # Empty string if '08' not in log
    if '09' in log:
        data.append(str(log['09']['output']['final_sequence_length']))  # 'assembly_length'
    else:
        data.append('')  # Empty string if '09' not in log
    if '05' in log:
        data.append(str(log['05']['output']['itr_length']))  # 'itr_length'
        data.append(str(log['05']['output']['itrs']))  # 'itr_length'

    with open(summary_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join(data) + '\n')

    with open(log_file, 'w') as f:
        json_object = json.dumps(log, indent=4)
        f.write(json_object)
    sys.exit(0)

def process_graph(gfa_graph, output_dir, input_file, reference, log):
    """Process the graph and write output."""
    try:
        input_with_ext = os.path.basename(input_file)
        input_base, _ = os.path.splitext(input_with_ext)
        sample_name = input_base.replace('.bridges_applied', '')

        # Write all contigs to a FASTA file at the start
        write_all_contigs(gfa_graph.segments, os.path.join(output_dir, sample_name + ".contigs.fasta"))

        # Start processing the graph
        filtered_edges, status = remove_self_loops(gfa_graph.edges)
        log['01'] = {
            'step_name': "remove_self_loops",
            'step_description': "Remove self-loops from link information",
            'status': status,
            'output': {'filtered_edges': [str(edge) for edge in filtered_edges]}
        }
        if status != "PASS":
            write_log_and_exit(log, status)

        # Identify low-coverage contigs
        low_cov, status = identify_low_coverage_contigs(gfa_graph.segments)
        log['02'] = {
            'step_name': "identify_low_coverage_contigs",
            'step_description': "Identify contigs with coverage below 0.5x",
            'status': status,
            'output': {'low_coverage_contigs': low_cov}
        }
        if status != "PASS":
            write_log_and_exit(log, status)

        # Filter links to exclude those involving low-coverage contigs
        filtered_links, status = filter_links_by_coverage(low_cov, filtered_edges)
        log['03'] = {
            'step_name': "filter_links_by_coverage",
            'step_description': "Filter links to exclude those involving low-coverage contigs",
            'status': status,
            'output': {'filtered_links': [str(edge) for edge in filtered_links]}
        }
        if status != "PASS":
            write_log_and_exit(log, status)
        # Create filtered graph using filtered links
        filtered_graph, status = create_filtered_graph(filtered_links, gfa_graph.segments)
        log['04'] = {
            'step_name': "create_filtered_graph",
            'step_description': "Create a filtered NetworkX graph from links and remove disconnected nodes",
            'status': status,
            'output': {'filtered_graph': str(filtered_graph)}
        }
        if status != "PASS":
            write_log_and_exit(log, status)

        # Filter segments to include only those that are part of the filtered graph
        filtered_segments = filter_segments_by_graph(gfa_graph.segments, filtered_graph)

        # Find all longest paths        
        final_paths, itrs, itr_length, status = get_final_paths(filtered_links, filtered_graph, filtered_segments)
        log['05'] = {
            'step_name': "get_final_paths",
            'step_description': "Find all longest paths in the graph starting from any ITR",
            'status': status,
            'output': {'final_paths': final_paths, 'itr_length': itr_length, 'itrs': itrs}
        }
        
        if not status.startswith("PASS"):
            write_log_and_exit(log, status)

        # Check for problems with the ITRs
        if not itrs or len(itrs) == 0:
            status = "FAIL: No ITRs identified in the graph."
            log['05']['status'] = status
            write_log_and_exit(log, status)

        # Find the longest contig
        longest_contig, status, longest_contig_file = find_longest_contig(filtered_segments, output_dir)
        log['06'] = {
            'step_name': "find_longest_contig",
            'step_description': "Find the longest contig based on length",
            'status': status,
            'output': {'longest_contig': longest_contig}
        }
        if status != "PASS":
            write_log_and_exit(log, status)
        # Create blast database directory
        blast_db_dir = os.path.join(output_dir, 'blast_db')
        os.makedirs(blast_db_dir, exist_ok=True)  # Create if it doesn't exist
        # Determine the orientation of the longest contig
        longest_orient, status = orient_longest_contig(longest_contig_file, reference, blast_db_dir)
        log['07'] = {
            'step_name': "orient_longest_contig",
            'step_description': "Determine the orientation of the longest contig",
            'status': status,
            'output': {'longest_orient': longest_orient}
        }
        if status != "PASS":
            write_log_and_exit(log, status)
        # Get the final orientation for all contigs in each path using filtered links
        final_path, final_oriented_path, status = get_final_orientation(
            final_paths, filtered_links, longest_contig, longest_orient, itrs)
        
        log['08'] = {
            'step_name': "get_final_orientation",
            'step_description': "Get the final orientation for all contigs",
            'status': status,
            'output': {'final_path': final_path, 'final_orientation': final_oriented_path}
        }
        if status != "PASS" or final_path is None:
            write_log_and_exit(log, "FAIL: No valid path found with the longest contig in the correct orientation.")

        # Generate the final sequence for this path
        final_sequence, seq_len, final_order_orientation_copy_number, status = get_final_sequence(
            final_oriented_path, filtered_segments)
        if status != "PASS":
            write_log_and_exit(log, "FAIL: Sequence generation failed for the final path.")

        # Write the sequence to a FASTA file
        asm_name = sample_name + ".assembly_asm.fasta"
        assembly_file = os.path.join(output_dir, asm_name)
        record = SeqRecord(Seq(final_sequence), id=input_base, description=final_order_orientation_copy_number)
        with open(assembly_file, 'w') as f:
            SeqIO.write(record, f, "fasta")

        # Update log with the final selected assembly        
        log['09'] = {
            'step_name': "generate_final_sequence",
            'step_description': "Generate the final sequence based on the best path and orientation",
            'status': "PASS",
            'output': {
                'final_sequence_length': seq_len,
                'final_order_orientation_copy_number': final_order_orientation_copy_number,
                'assembly_file': asm_name,
                'final_path': final_path,
                'final_orientation': final_oriented_path
            }
        }

        # Remove the blast_db directory
        shutil.rmtree(blast_db_dir)  # Remove the directory and its contents

        # Write the final log and summary
        write_log_and_exit(log, "PASS")

    except Exception as e:
        log['Exception'] = {'error': str(e)}
        write_log_and_exit(log, "FAIL: Exception occurred during processing")
        return
    
def main(arguments):
    parser = argparse.ArgumentParser(description="GFA parser to construct assembly from Unicycler")
    parser.add_argument('-i', type=str, required=True, help="GFA file from Unicycler")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-r', type=str, required=True, help="Reference file for orientation")
    args = parser.parse_args(arguments)

    gfa_file = args.i
    output_dir = args.o
    reference_file = args.r

    input_with_ext = os.path.basename(gfa_file)
    input_base, _ = os.path.splitext(input_with_ext)
    sample_name = input_base.replace('.bridges_applied', '')

    # Initialize log
    log = {}
    log['00'] = {'step_name': "initialization",
                 'step_description': "Initialize the logging process",
                 'status': "PASS",
                 'input': {
                     'gfa': gfa_file,
                     'sample_name': sample_name,
                     'output_dir': output_dir
                 },
                 'output': {}}

    try:
        # clean up gfa tags
        cleanedGfa = 'graph_cleaned.gfa'
        clean_graph_tags(gfa_file, cleanedGfa)

        gfa_graph, status = read_gfa_file(cleanedGfa)
        if status != "PASS":
            write_log_and_exit(log, status)

        if gfa_graph is None:
            print("Failed to read GFA file.")
            write_log_and_exit(log, "FAIL: Failed to read GFA file")
            return

        process_graph(gfa_graph, output_dir, gfa_file, reference_file, log)
    except Exception as e:
        log['Exception'] = {'error': str(e)}
        write_log_and_exit(log, "FAIL: Exception occurred in main")

if __name__ == '__main__':
    main(sys.argv[1:])