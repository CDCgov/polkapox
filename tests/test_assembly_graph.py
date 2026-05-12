#!/usr/bin/env python3
"""
Unit tests for AssemblyGraph_gfaPy.py

Run with: pytest tests/test_assembly_graph.py -v

Requires: gfapy, networkx, biopython (available inside the staphb/polkapox container
or via: pip install gfapy networkx biopython)

Test GFAs are in tests/nf-test-data/ and cover:
  - Clade II standard (PASS)
  - Clade I standard (PASS x2)
  - Clade II single-ITR (PASS)
  - Clade II broken graph (FAIL)
  - Clade I broken graph (FAIL)
  - Clade II recombinant (FAIL)
"""

import os
import sys
import json
import pytest
import tempfile

# Make bin/ importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

pytest.importorskip("gfapy")
pytest.importorskip("networkx")
pytest.importorskip("Bio")

import networkx as nx
from AssemblyGraph_gfaPy import (
    clean_graph_tags,
    read_gfa_file,
    remove_self_loops_and_terminal_nodes,
    remove_self_loops,
    identify_low_coverage_contigs,
    filter_links_by_coverage,
    create_filtered_graph,
    find_longest_contig,
    filter_segments_by_graph,
    write_log_and_exit,
    Link,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DATA = os.path.join(os.path.dirname(__file__), 'nf-test-data')

PASS_GFAS = {
    'cladeII_standard':  os.path.join(DATA, 'cladeII_standard.bridges_applied.gfa'),
    'cladeII_one_itr':   os.path.join(DATA, 'cladeII_one_itr.bridges_applied.gfa'),
    'cladeI_standard':   os.path.join(DATA, 'cladeI_standard.bridges_applied.gfa'),
    'cladeI_standard2':  os.path.join(DATA, 'cladeI_standard2.bridges_applied.gfa'),
}

FAIL_GFAS = {
    'cladeII_broken':      os.path.join(DATA, 'cladeII_broken.bridges_applied.gfa'),
    'cladeI_broken':       os.path.join(DATA, 'cladeI_broken.bridges_applied.gfa'),
    'cladeII_recombinant': os.path.join(DATA, 'cladeII_recombinant.bridges_applied.gfa'),
}


def load_cleaned_gfa(gfa_path):
    """Clean tags and load a GFA file, returning (gfa_graph, status)."""
    with tempfile.NamedTemporaryFile(suffix='.gfa', delete=False, mode='w') as tmp:
        tmp_path = tmp.name
    try:
        clean_graph_tags(gfa_path, tmp_path)
        return read_gfa_file(tmp_path)
    finally:
        os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# clean_graph_tags
# ---------------------------------------------------------------------------

def test_clean_graph_tags_removes_lb_cl():
    """LB and CL tags should be stripped from segment lines."""
    raw = "S\t1\tACGT\tLN:i:4\tLB:z:somelib\tCL:z:somecolor\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False) as f_in:
        f_in.write(raw)
        in_path = f_in.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False) as f_out:
        out_path = f_out.name
    try:
        clean_graph_tags(in_path, out_path)
        content = open(out_path).read()
        assert 'LB:z:' not in content
        assert 'CL:z:' not in content
        assert 'ACGT' in content
    finally:
        os.unlink(in_path)
        os.unlink(out_path)


# ---------------------------------------------------------------------------
# read_gfa_file
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,path", list(PASS_GFAS.items()))
def test_read_gfa_pass(label, path):
    """All PASS GFAs should load without error."""
    gfa, status = load_cleaned_gfa(path)
    assert gfa is not None, f"{label}: expected successful load, got status={status}"
    assert status == "PASS"
    assert len(gfa.segments) > 0


def test_read_gfa_empty_edges(tmp_path):
    """A GFA with only segment lines (no edges) should return WARNING status."""
    gfa_content = "S\t1\tACGT\tLN:i:4\tdp:f:1.0\n"
    p = tmp_path / "no_edges.gfa"
    p.write_text(gfa_content)
    gfa, status = read_gfa_file(str(p))
    assert gfa is None
    assert "WARNING" in status


# ---------------------------------------------------------------------------
# remove_self_loops_and_terminal_nodes
# ---------------------------------------------------------------------------

def test_remove_self_loops():
    """Self-loop links (from == to) should be removed."""
    links = [Link('1', '1'), Link('1', '2'), Link('2', '3')]
    result, status = remove_self_loops(links)
    names = [(l.from_name, l.to_name) for l in result]
    assert ('1', '1') not in names
    assert ('1', '2') in names
    assert status == "PASS"


# ---------------------------------------------------------------------------
# create_filtered_graph — correctness of Graph vs MultiGraph
# ---------------------------------------------------------------------------

def test_create_filtered_graph_uses_simple_graph():
    """Graph should be nx.Graph (not MultiGraph) so parallel GFA edges don't double degrees."""
    # Two links between same pair (as Unicycler emits for both strand orientations)
    links = [Link('1', '2'), Link('2', '1'), Link('2', '3'), Link('3', '1')]
    segments = []  # not used for graph structure
    graph, status = create_filtered_graph(links, segments)
    assert graph is not None, f"Expected valid graph, got status={status}"
    # In a simple graph, node 1 should have degree 2 (connected to 2 and 3), not 3 or 4
    assert graph.degree('1') == 2


def test_create_filtered_graph_cycle_detection():
    """A 3-node cycle should be detected as circular."""
    links = [Link('A', 'B'), Link('B', 'C'), Link('C', 'A')]
    graph, status = create_filtered_graph(links, [])
    assert graph is not None
    assert status == "PASS"


def test_create_filtered_graph_non_circular_fails():
    """A path graph (no cycle) should return WARNING."""
    links = [Link('A', 'B'), Link('B', 'C')]
    graph, status = create_filtered_graph(links, [])
    assert graph is None
    assert "WARNING" in status


# ---------------------------------------------------------------------------
# find_longest_contig — LN fallback
# ---------------------------------------------------------------------------

class FakeSeg:
    def __init__(self, name, ln=None, sequence=''):
        self.name = name
        self.sequence = sequence
        self._ln = ln
    def get(self, tag):
        if tag == 'LN':
            return self._ln
        if tag == 'name':
            return self.name


def test_find_longest_contig_with_ln(tmp_path):
    segs = [FakeSeg('1', ln=100, sequence='A'*100),
            FakeSeg('2', ln=200, sequence='A'*200),
            FakeSeg('3', ln=50,  sequence='A'*50)]
    name, status, _ = find_longest_contig(segs, str(tmp_path))
    assert name == '2'
    assert status == "PASS"


def test_find_longest_contig_without_ln(tmp_path):
    """Should fall back to len(sequence) when LN tag is absent."""
    segs = [FakeSeg('1', ln=None, sequence='A'*100),
            FakeSeg('2', ln=None, sequence='A'*300),
            FakeSeg('3', ln=None, sequence='A'*50)]
    name, status, _ = find_longest_contig(segs, str(tmp_path))
    assert name == '2'
    assert status == "PASS"


# ---------------------------------------------------------------------------
# identify_low_coverage_contigs
# ---------------------------------------------------------------------------

class FakeSegWithDp:
    def __init__(self, name, dp):
        self.name = name
        self._dp = dp
        self.tagnames = ['dp'] if dp is not None else []
    def get(self, tag):
        if tag == 'name': return self.name
        if tag == 'dp': return self._dp
        return None


def test_identify_low_coverage_contigs():
    segs = [FakeSegWithDp('high', 2.5),
            FakeSegWithDp('low', 0.3),
            FakeSegWithDp('borderline', 0.6)]
    low_cov, status = identify_low_coverage_contigs(segs)
    assert 'low' in low_cov
    assert 'high' not in low_cov
    assert status == "PASS"


# ---------------------------------------------------------------------------
# write_log_and_exit — exit code correctness
# ---------------------------------------------------------------------------

def _make_minimal_log(tmp_path, sample='test_sample'):
    return {
        '00': {
            'step_name': 'initialization',
            'step_description': '',
            'status': 'PASS',
            'input': {'gfa': 'test.gfa', 'sample_name': sample, 'output_dir': str(tmp_path)},
            'output': {}
        }
    }


def test_write_log_and_exit_pass_exits_zero(tmp_path):
    log = _make_minimal_log(tmp_path)
    with pytest.raises(SystemExit) as exc_info:
        write_log_and_exit(log, "PASS")
    assert exc_info.value.code == 0


def test_write_log_and_exit_fail_exits_nonzero(tmp_path):
    log = _make_minimal_log(tmp_path)
    with pytest.raises(SystemExit) as exc_info:
        write_log_and_exit(log, "FAIL: something went wrong")
    assert exc_info.value.code != 0


def test_write_log_and_exit_warning_exits_nonzero(tmp_path):
    log = _make_minimal_log(tmp_path)
    with pytest.raises(SystemExit) as exc_info:
        write_log_and_exit(log, "WARNING: graph is not circular")
    assert exc_info.value.code != 0


def test_write_log_and_exit_writes_summary_columns(tmp_path):
    log = _make_minimal_log(tmp_path, sample='mysample')
    log['09'] = {'output': {'final_sequence_length': 197000}}
    log['08'] = {'output': {'final_orientation': ['1 +', '2 -']}}
    log['05'] = {'output': {'itr_length': 6500, 'final_paths': [], 'itrs': ['1']}}
    with pytest.raises(SystemExit):
        write_log_and_exit(log, "PASS")
    summary = open(tmp_path / 'mysample.assembly.summary').read().splitlines()
    headers = summary[0].split('\t')
    values  = summary[1].split('\t')
    assert 'sample' in headers
    assert 'assembly_length' in headers
    assert 'itr_length' in headers
    assert len(headers) == len(values), "Header and data row must have the same number of columns"
    idx_sample = headers.index('sample')
    assert values[idx_sample] == 'mysample'
    idx_asm = headers.index('assembly_length')
    assert values[idx_asm] == '197000'


def test_write_log_and_exit_writes_json_log(tmp_path):
    log = _make_minimal_log(tmp_path)
    with pytest.raises(SystemExit):
        write_log_and_exit(log, "PASS")
    log_data = json.loads(open(tmp_path / 'test_sample.assembly.log').read())
    assert '00' in log_data


# ---------------------------------------------------------------------------
# Integration: load real GFAs and run through graph processing steps
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,path", list(PASS_GFAS.items()))
def test_pass_gfa_graph_is_circular(label, path):
    """PASS GFAs should produce a connected circular graph after standard filtering."""
    gfa, status = load_cleaned_gfa(path)
    assert gfa is not None, f"{label}: GFA load failed: {status}"

    links_obj = [Link(e.from_name, e.to_name) for e in gfa.edges]
    filtered, _ = remove_self_loops_and_terminal_nodes(links_obj, gfa.segments)
    low_cov, _ = identify_low_coverage_contigs(gfa.segments)
    filtered, _ = filter_links_by_coverage(low_cov, filtered)
    graph, graph_status = create_filtered_graph(filtered, gfa.segments)

    assert graph is not None, f"{label}: expected circular graph, got: {graph_status}"
    assert graph_status == "PASS", f"{label}: {graph_status}"


@pytest.mark.parametrize("label,path", list(FAIL_GFAS.items()))
def test_fail_gfa_graph_is_not_circular(label, path):
    """FAIL GFAs should either fail to load or produce a non-circular graph."""
    gfa, load_status = load_cleaned_gfa(path)
    if gfa is None:
        # Acceptable — some broken GFAs fail at load time
        assert "WARNING" in load_status or "FAIL" in load_status
        return

    links_obj = [Link(e.from_name, e.to_name) for e in gfa.edges]
    filtered, _ = remove_self_loops_and_terminal_nodes(links_obj, gfa.segments)
    low_cov, _ = identify_low_coverage_contigs(gfa.segments)
    filtered, _ = filter_links_by_coverage(low_cov, filtered)
    graph, graph_status = create_filtered_graph(filtered, gfa.segments)

    assert graph is None or "WARNING" in graph_status, (
        f"{label}: expected non-circular graph but got PASS"
    )
