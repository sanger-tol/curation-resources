#!/usr/bin/env python3
"""
Rename and orient chromosomes based on reference alignment (PAF format).

This script renames chromosomes in a FASTA file and changes their orientation
based on alignment to a reference genome.
"""

__version__ = "1.1.1"

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

def _pearson_r(xs: List[float], ys: List[float]) -> float:
    """
    Pearson correlation coefficient between two lists of floats.
    Returns 0.0 if fewer than 2 points or zero variance in either variable.
    """
    n = len(xs)
    if n < 2:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx  = sum((x - mx) ** 2 for x in xs) ** 0.5
    sy  = sum((y - my) ** 2 for y in ys) ** 0.5
    if sx == 0 or sy == 0:
        return 0.0
    return num / (sx * sy)


def needs_reverse_complement_by_correlation(
    records: List,
    target_length: int,
) -> bool:
    """
    Decide whether a query chromosome needs reverse-complement using Pearson
    correlation between query midpoints and target midpoints.

    For '+' blocks: r = corr(query_mid, target_mid).
    For '-' blocks: r = corr(query_mid, target_length - target_mid).

    The strand whose r is higher wins.  This correctly ignores large structural
    inversions (their blocks are scattered → low r) while syntenic blocks
    form a clean diagonal → high r.

    Args:
        records:       PAFRecord objects for one query→best_target pair.
        target_length: Length of the target chromosome.

    Returns:
        True  if reverse complement is needed,
        False if the chromosome is already correctly oriented.
    """
    plus_recs  = [r for r in records if r.strand == '+']
    minus_recs = [r for r in records if r.strand == '-']

    r_plus = 0.0
    if len(plus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in plus_recs]
        ys = [(r.target_start + r.target_end) / 2 for r in plus_recs]
        r_plus = _pearson_r(xs, ys)

    r_minus = 0.0
    if len(minus_recs) >= 2:
        xs = [(r.query_start + r.query_end) / 2 for r in minus_recs]
        ys = [target_length - (r.target_start + r.target_end) / 2 for r in minus_recs]
        r_minus = _pearson_r(xs, ys)

    # If correlation is uninformative for both strands (< 2 blocks each),
    # fall back to naive length comparison.
    if r_plus == 0.0 and r_minus == 0.0:
        plus_len  = sum(r.query_end - r.query_start for r in plus_recs)
        minus_len = sum(r.query_end - r.query_start for r in minus_recs)
        return minus_len > plus_len

    return r_minus > r_plus


@dataclass
class PAFRecord:
    """Single PAF alignment record."""
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    alignment_length: int
    mapping_quality: int
    
    @property
    def alignment_block_length(self) -> int:
        """Length of the alignment block on query."""
        return self.query_end - self.query_start


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rename and orient chromosomes based on reference alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"%(prog)s {__version__}"
    )
    
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        type=Path,
        help="Input FASTA file (can be gzipped)"
    )
    
    # Alignment source: either a PAF file or a pre-built mapping table
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--paf", "-p",
        type=Path,
        default=None,
        help="PAF file with alignment to reference"
    )
    source_group.add_argument(
        "--mapping-table", "-mt",
        type=Path,
        default=None,
        dest="mapping_table",
        help="Pre-built mapping TSV (output of a previous run) to rename/orient a second "
             "haplotype without re-running alignment. Columns used: query, renamed_to, "
             "needs_reverse_complement."
    )
    
    parser.add_argument(
        "--output-dir", "-d",
        type=Path,
        default=Path("./rename_and_orient"),
        help="Output directory for generated files (will be created if it doesn't exist)"
    )
    
    parser.add_argument(
        "--output-prefix", "-o",
        type=str,
        help="Prefix for output file names (default: derived from input FASTA file name)"
    )
    
    parser.add_argument(
        "--min-coverage", "-c",
        type=float,
        default=0.5,
        help="Minimum coverage threshold for renaming (0.0-1.0)"
    )
    
    parser.add_argument(
        "--query-chromosome-prefix", "-q",
        type=str,
        default="SUPER_",
        help="Prefix for query chromosome names in input FASTA (e.g., SUPER_, scaffold_, contig_)"
    )
    
    parser.add_argument(
        "--output-chromosome-prefix", "-x",
        type=str,
        default="SUPER_",
        help="Prefix for output chromosome names (e.g., SUPER_, chr_, chr, or empty string for no prefix)"
    )
    
    parser.add_argument(
        "--reference-chromosome-prefix", "-r",
        type=str,
        default=None,
        help="Prefix for reference chromosome names in PAF target (e.g., chr_, chr, SUPER_, scaffold_). "
             "Auto-detected from PAF if not specified."
    )
    
    parser.add_argument(
        "--plot-alignments", "-P",
        action="store_true",
        default=False,
        help="Generate scatter plots of PAF alignment blocks for each chromosome "
             "(saved to <output-dir>/plots/). Requires matplotlib."
    )

    args = parser.parse_args()
    
    # Set default output prefix if not provided
    if args.output_prefix is None:
        args.output_prefix = args.fasta.stem
    
    return args


def read_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Read FASTA file (supports gzip compression).
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    # Determine if file is gzipped
    open_func = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'
    
    with open_func(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                
                # Start new sequence - take first word as name
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def parse_paf(paf_path: Path) -> List[PAFRecord]:
    """
    Parse PAF file.
    
    PAF format (tab-separated):
    1. Query sequence name
    2. Query sequence length
    3. Query start (0-based)
    4. Query end (0-based, open)
    5. Strand ('+' or '-')
    6. Target sequence name
    7. Target sequence length
    8. Target start (0-based)
    9. Target end (0-based, open)
    10. Number of matching bases
    11. Alignment block length
    12. Mapping quality (0-255)
    
    Args:
        paf_path: Path to PAF file
        
    Returns:
        List of PAFRecord objects
    """
    records = []
    
    with open(paf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 12:
                continue
            
            record = PAFRecord(
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                num_matches=int(fields[9]),
                alignment_length=int(fields[10]),
                mapping_quality=int(fields[11])
            )
            records.append(record)
    
    return records


def detect_reference_prefix(records: List[PAFRecord]) -> str:
    """
    Auto-detect reference chromosome prefix from PAF target names.

    Sums query alignment length per prefix and returns the winner, which
    correctly favours main chromosomes (chr_, SUPER_) over many short-hit
    unplaced contigs (NW_, NT_).  Logs the result to stderr.

    Raises ValueError if no alphabetic prefix is found — use
    --reference-chromosome-prefix to specify it explicitly.
    """
    prefix_aln_length: Dict[str, int] = defaultdict(int)

    for record in records:
        m = re.match(r'^([A-Za-z_\.]+)', record.target_name)
        if m:
            prefix_aln_length[m.group(1)] += record.query_end - record.query_start

    if not prefix_aln_length:
        raise ValueError(
            "Could not detect reference chromosome prefix from PAF file: "
            "no target names with a recognisable alphabetic prefix were found. "
            "Please specify --reference-chromosome-prefix explicitly."
        )

    detected = max(prefix_aln_length, key=lambda p: prefix_aln_length[p])
    print(f"[detect_reference_prefix] Auto-detected reference prefix: '{detected}' "
          f"(total aligned bases: {prefix_aln_length[detected]:,})", file=sys.stderr)
    return detected


def filter_paf_records(records: List[PAFRecord], 
                       query_chromosome_prefix: str = "SUPER_",
                       target_prefix: str = None) -> Tuple[List[PAFRecord], str]:
    """
    Filter PAF records to keep only those matching query and target prefixes.
    Auto-detects target prefix if not specified.
    
    Args:
        records: List of PAF records
        query_chromosome_prefix: Required prefix for query names
        target_prefix: Required prefix for target names (auto-detected if None)
        
    Returns:
        Tuple of (filtered records, detected target prefix)
    """
    if target_prefix is None:
        target_prefix = detect_reference_prefix(records)
    
    filtered = []
    for record in records:
        if record.query_name.startswith(query_chromosome_prefix) and \
           record.target_name.startswith(target_prefix):
            filtered.append(record)
    return filtered, target_prefix


def validate_paf_fasta_consistency(
    paf_records: List[PAFRecord],
    fasta_sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_"
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate consistency between PAF and FASTA files for chromosomes with given prefix.
    
    Checks that all chromosomes with query_chromosome_prefix in PAF exist in FASTA and vice versa.
    
    Args:
        paf_records: List of PAF records
        fasta_sequences: Dictionary of FASTA sequences
        query_chromosome_prefix: Prefix for chromosome names
        
    Returns:
        Tuple of:
        - is_valid: True if all chromosomes match
        - in_paf_not_fasta: List of chromosomes in PAF but not in FASTA
        - in_fasta_not_paf: List of chromosomes in FASTA but not in PAF
    """
    # Get unique chromosome names from PAF (excluding unloc)
    paf_chromosomes = set()
    for record in paf_records:
        if record.query_name.startswith(query_chromosome_prefix) and "_unloc_" not in record.query_name:
            paf_chromosomes.add(record.query_name)
    
    # Get chromosome names from FASTA (excluding unloc)
    fasta_chromosomes = set()
    for name in fasta_sequences.keys():
        if name.startswith(query_chromosome_prefix) and "_unloc_" not in name:
            fasta_chromosomes.add(name)
    
    # Find mismatches
    in_paf_not_fasta = sorted(paf_chromosomes - fasta_chromosomes)
    in_fasta_not_paf = sorted(fasta_chromosomes - paf_chromosomes)
    
    is_valid = len(in_paf_not_fasta) == 0 and len(in_fasta_not_paf) == 0
    
    return is_valid, in_paf_not_fasta, in_fasta_not_paf


@dataclass
class ChromosomeMapping:
    """Mapping information for a single chromosome."""
    query_name: str
    query_length: int
    target_name: str
    total_alignment_length: int
    coverage: float
    plus_strand_length: int
    minus_strand_length: int
    needs_reverse_complement: bool
    target_prefix: str = "chr_"  # Reference prefix (chr_ or chr)
    
    @property
    def target_suffix(self) -> str:
        """Extract suffix from target name (e.g., 'chr_5' -> '5', 'chrW' -> 'W')."""
        if self.target_name.startswith(self.target_prefix):
            return self.target_name[len(self.target_prefix):]
        return self.target_name


@dataclass
class UnlocMapping:
    """Mapping information for an unlocalized contig (SUPER_N_unloc_M)."""
    contig_name: str
    parent_chromosome: str  # e.g., "SUPER_5" for "SUPER_5_unloc_1"
    unloc_number: int       # e.g., 1 for "SUPER_5_unloc_1"
    needs_reverse_complement: bool  # Inherited from parent chromosome


@dataclass
class FinalChromosomeAssignment:
    """Final assignment for a chromosome after conflict resolution."""
    original_name: str           # e.g., "SUPER_1"
    new_name: str               # e.g., "SUPER_13" (after renaming)
    new_suffix: str             # e.g., "13" or "W"
    needs_reverse_complement: bool
    is_sex_chromosome: bool
    

# Sex chromosome patterns (letters or letters with numbers like Z1, Z2, X1X2, W1W2, B1B2)
SEX_CHROMOSOME_SUFFIXES = {'W', 'X', 'Y', 'Z', 'B'}


def is_sex_chromosome_suffix(suffix: str) -> bool:
    """
    Check if a suffix indicates a sex chromosome.
    
    Sex chromosomes have non-numeric suffixes like W, Z, X, Y, Z1, Z2, etc.
    
    Args:
        suffix: The suffix to check (e.g., "5", "W", "Z1")
        
    Returns:
        True if this is a sex chromosome suffix
    """
    if not suffix:
        return False
    
    # Check if the first character is a known sex chromosome letter
    first_char = suffix[0].upper()
    if first_char in SEX_CHROMOSOME_SUFFIXES:
        # Remaining characters should be digits or empty (W, Z, Z1, Z2, etc.)
        remaining = suffix[1:]
        return remaining == '' or remaining.isdigit()
    
    return False


def extract_chromosome_suffix(name: str, prefix: str = "SUPER_") -> str:
    """
    Extract chromosome suffix from name.
    
    Args:
        name: Chromosome name (e.g., "SUPER_5", "chr_W", "chrZ1", "scaffold_10")
        prefix: Prefix to strip
        
    Returns:
        Suffix string (e.g., "5", "W", "Z1")
    """
    if name.startswith(prefix):
        return name[len(prefix):]
    return name


def is_autosome_suffix(suffix: str) -> bool:
    """
    Check if suffix represents an autosome (purely numeric).
    
    Args:
        suffix: Chromosome suffix
        
    Returns:
        True if this is an autosome suffix (numeric)
    """
    return suffix.isdigit()
    

def is_unloc_contig(name: str) -> bool:
    """
    Check if sequence name is an unlocalized contig.
    
    Unloc contigs have format: SUPER_N_unloc_M
    
    Args:
        name: Sequence name
        
    Returns:
        True if name matches unloc pattern
    """
    return "_unloc_" in name


def parse_unloc_name(name: str) -> Tuple[str, int]:
    """
    Parse unlocalized contig name to extract parent chromosome and unloc number.
    
    Args:
        name: Contig name like "SUPER_5_unloc_1"
        
    Returns:
        Tuple of (parent_chromosome, unloc_number), e.g., ("SUPER_5", 1)
    """
    # SUPER_5_unloc_1 -> ["SUPER_5", "1"]
    parts = name.split("_unloc_")
    parent = parts[0]
    unloc_num = int(parts[1]) if len(parts) > 1 else 0
    return parent, unloc_num


def group_alignments_by_query(records: List[PAFRecord]) -> Dict[str, List[PAFRecord]]:
    """
    Group PAF records by query chromosome name.
    
    Args:
        records: List of PAF records
        
    Returns:
        Dictionary mapping query names to lists of PAF records
    """
    groups = defaultdict(list)
    for record in records:
        groups[record.query_name].append(record)
    return dict(groups)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping intervals into non-overlapping ones.
    
    Args:
        intervals: List of (start, end) tuples
        
    Returns:
        Sorted list of non-overlapping (start, end) tuples
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def calculate_target_alignments(records: List[PAFRecord]) -> Dict[str, Dict[str, int]]:
    """
    Calculate total alignment lengths for each query-target pair.
    
    Overlapping query intervals on the same target are merged before summing
    to avoid double-counting repeated alignments to the same region
    (e.g. when a repetitive reference chromosome attracts many duplicate hits).
    Strand totals use the same de-duplicated intervals per strand.
    
    Args:
        records: List of PAF records for a single query
        
    Returns:
        Dictionary mapping target names to dict with 'total', 'plus', 'minus' lengths
    """
    # Collect intervals per target, separated by strand
    target_intervals: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(
        lambda: {'+': [], '-': []}
    )

    for record in records:
        strand = record.strand
        interval = (record.query_start, record.query_end)
        target_intervals[record.target_name][strand].append(interval)

    target_stats: Dict[str, Dict[str, int]] = {}
    for target_name, strands in target_intervals.items():
        plus_len = sum(e - s for s, e in merge_intervals(strands['+'
                                                                    ]))
        minus_len = sum(e - s for s, e in merge_intervals(strands['-']))
        # For the total, merge all intervals together (regardless of strand)
        all_intervals = strands['+'] + strands['-']
        total_len = sum(e - s for s, e in merge_intervals(all_intervals))
        target_stats[target_name] = {
            'total': total_len,
            'plus': plus_len,
            'minus': minus_len,
        }

    return target_stats


def determine_best_target(target_stats: Dict[str, Dict[str, int]]) -> Tuple[str, Dict[str, int]]:
    """
    Determine the best target chromosome based on total alignment length.
    
    Args:
        target_stats: Dictionary from calculate_target_alignments
        
    Returns:
        Tuple of (best_target_name, stats_dict)
    """
    best_target = None
    best_stats = None
    max_length = 0
    
    for target_name, stats in target_stats.items():
        if stats['total'] > max_length:
            max_length = stats['total']
            best_target = target_name
            best_stats = stats
    
    return best_target, best_stats


def build_chromosome_mappings(
    records: List[PAFRecord],
    min_coverage: float = 0.5,
    target_prefix: str = "chr_"
) -> List[ChromosomeMapping]:
    """
    Build chromosome mappings from PAF records.
    
    For each query chromosome:
    1. Calculate alignment lengths to each target
    2. Select best target (maximum alignment length)
    3. Check coverage threshold
    4. Determine orientation based on strand statistics
    
    Args:
        records: Filtered PAF records
        min_coverage: Minimum coverage threshold (0.0-1.0)
        target_prefix: Prefix used in reference chromosome names (chr_ or chr)
        
    Returns:
        List of ChromosomeMapping objects
    """
    # Get query lengths from records
    query_lengths = {}
    for record in records:
        query_lengths[record.query_name] = record.query_length
    
    # Group alignments by query
    groups = group_alignments_by_query(records)
    
    mappings = []
    
    for query_name, query_records in groups.items():
        query_length = query_lengths[query_name]
        
        # Calculate alignment stats for each target
        target_stats = calculate_target_alignments(query_records)
        
        # Find best target
        best_target, best_stats = determine_best_target(target_stats)
        
        if best_target is None:
            continue
        
        # Calculate coverage
        coverage = best_stats['total'] / query_length
        
        # Skip if coverage below threshold
        if coverage < min_coverage:
            if not is_unloc_contig(query_name):
                print(f"  Warning: {query_name} -> {best_target} coverage {coverage:.2%} below threshold")
            continue
        
        # Determine orientation using Pearson correlation method.
        # Correlation between query_mid and target_mid identifies the dominant
        # syntenic strand while ignoring large structural inversions.
        best_target_records = [r for r in query_records if r.target_name == best_target]
        target_length = best_target_records[0].target_length if best_target_records else 0
        needs_rc = needs_reverse_complement_by_correlation(best_target_records, target_length)

        mapping = ChromosomeMapping(
            query_name=query_name,
            query_length=query_length,
            target_name=best_target,
            total_alignment_length=best_stats['total'],
            coverage=coverage,
            plus_strand_length=best_stats['plus'],
            minus_strand_length=best_stats['minus'],
            needs_reverse_complement=needs_rc,
            target_prefix=target_prefix
        )
        mappings.append(mapping)
    
    return mappings


def build_unloc_mappings(
    sequences: Dict[str, str],
    chromosome_mappings: List[ChromosomeMapping]
) -> List[UnlocMapping]:
    """
    Build mappings for unlocalized contigs based on parent chromosome mappings.
    
    Unloc contigs (SUPER_N_unloc_M) inherit orientation from their parent
    chromosome (SUPER_N). They don't need alignment-based mapping.
    
    Args:
        sequences: Dictionary of all sequences (name -> sequence)
        chromosome_mappings: List of chromosome mappings
        
    Returns:
        List of UnlocMapping objects
    """
    # Build lookup: parent chromosome name -> needs_reverse_complement
    parent_rc_lookup = {m.query_name: m.needs_reverse_complement for m in chromosome_mappings}
    
    unloc_mappings = []
    
    for seq_name in sequences.keys():
        if not is_unloc_contig(seq_name):
            continue
        
        parent_chr, unloc_num = parse_unloc_name(seq_name)
        
        # Get orientation from parent (default to False if parent not mapped)
        needs_rc = parent_rc_lookup.get(parent_chr, False)
        
        unloc_mapping = UnlocMapping(
            contig_name=seq_name,
            parent_chromosome=parent_chr,
            unloc_number=unloc_num,
            needs_reverse_complement=needs_rc
        )
        unloc_mappings.append(unloc_mapping)
    
    return unloc_mappings

def natural_sort_key(name: str, prefix: str = "SUPER_") -> Tuple:
    """
    Generate a sort key for natural sorting of chromosome names.
    
    Sorts autosomes numerically (1, 2, ... 10, 11, ...) 
    and sex chromosomes alphabetically after autosomes.
    
    Args:
        name: Chromosome name (e.g., "SUPER_1", "SUPER_W")
        prefix: Prefix to strip (default "SUPER_")
        
    Returns:
        Tuple for sorting: (is_sex_chr, numeric_or_alpha_key)
    """
    suffix = extract_chromosome_suffix(name, prefix)
    
    if is_sex_chromosome_suffix(suffix):
        # Sex chromosomes come after autosomes, sorted alphabetically
        return (1, suffix)
    else:
        # Autosomes sorted numerically
        try:
            return (0, int(suffix))
        except ValueError:
            return (0, float('inf'))  # Unknown format goes to end


def print_mapping_summary(mappings: List[ChromosomeMapping]) -> None:
    """
    Print summary of chromosome mappings.
    
    Args:
        mappings: List of ChromosomeMapping objects
    """
    print("\nChromosome mapping summary:")
    print("-" * 80)
    print(f"{'Query':<15} {'Target':<10} {'Coverage':>10} {'Plus':>12} {'Minus':>12} {'RC?':>5}")
    print("-" * 80)
    
    for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
        print(f"{m.query_name:<15} {m.target_name:<10} {m.coverage:>9.1%} "
              f"{m.plus_strand_length:>12,} {m.minus_strand_length:>12,} "
              f"{'Yes' if m.needs_reverse_complement else 'No':>5}")


def save_mapping_tsv(
    mappings: List[ChromosomeMapping],
    assignments: List[FinalChromosomeAssignment],
    output_path: Path
) -> None:
    """
    Save chromosome mapping summary to TSV file.
    
    Args:
        mappings: List of ChromosomeMapping objects
        assignments: List of FinalChromosomeAssignment for renamed_to column
        output_path: Path to output TSV file
    """
    # Build lookup: original_name -> new_name
    rename_lookup = {a.original_name: a.new_name for a in assignments}
    
    with open(output_path, 'w') as f:
        # Header
        f.write("query\ttarget\trenamed_to\tquery_length\talignment_length\tcoverage\t"
                "plus_strand\tminus_strand\tneeds_reverse_complement\n")
        
        for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
            renamed_to = rename_lookup.get(m.query_name, m.query_name)
            f.write(f"{m.query_name}\t{m.target_name}\t{renamed_to}\t{m.query_length}\t"
                    f"{m.total_alignment_length}\t{m.coverage:.4f}\t"
                    f"{m.plus_strand_length}\t{m.minus_strand_length}\t"
                    f"{'yes' if m.needs_reverse_complement else 'no'}\n")
    
    print(f"Mapping summary saved to: {output_path}")


def calculate_genome_length(fasta_path: Path) -> int:
    """
    Calculate total genome length from FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Total length of all sequences
    """
    sequences = read_fasta(fasta_path)
    return sum(len(seq) for seq in sequences.values())


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    
    # Handle non-standard bases by keeping them as-is
    result = []
    for base in reversed(seq):
        result.append(complement.get(base, base))
    
    return ''.join(result)


def resolve_chromosome_assignments(
    mappings: List[ChromosomeMapping],
    sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_"
) -> Tuple[List[FinalChromosomeAssignment], Dict[str, bool]]:
    """
    Resolve chromosome assignments with conflict handling and sex chromosome logic.
    
    This function handles:
    1. Sex chromosome detection (W, Z, X, Y, Z1, Z2, etc.)
    2. Conflict resolution when multiple chromosomes map to same target
    3. Autosome -> sex chromosome mapping (reassign to max_autosome + N)
    4. Sex chromosome -> autosome mapping (skip the autosome number)
    5. Unmapped chromosomes keep original names (may cause mapped chr to shift)
    
    Args:
        mappings: List of ChromosomeMapping objects from PAF analysis
        sequences: Dictionary of all sequences (for finding unmapped ones)
        query_chromosome_prefix: Input chromosome prefix (e.g., "SUPER_", "scaffold_")
        output_prefix: Output chromosome prefix (e.g., "SUPER_", "chr_", "chr", "")
        
    Returns:
        Tuple of:
        - List of FinalChromosomeAssignment for all chromosomes
        - Dictionary mapping original names to needs_reverse_complement flag
    """
    assignments = []
    rc_lookup = {}  # original_name -> needs_rc
    
    # First, identify all chromosomes (excluding unloc contigs)
    all_chr_names = {name for name in sequences.keys() 
                     if name.startswith(query_chromosome_prefix) and "_unloc_" not in name}
    mapped_names = {m.query_name for m in mappings if "_unloc_" not in m.query_name}
    unmapped_names = all_chr_names - mapped_names
    
    # Separate autosomes and sex chromosomes based on QUERY names
    autosome_mappings = []
    sex_mappings = []
    
    for m in mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_chromosome_prefix)
        if "_unloc_" in m.query_name:
            continue
            
        if is_sex_chromosome_suffix(query_suffix):
            sex_mappings.append(m)
        else:
            autosome_mappings.append(m)
    
    # Reserved numbers: numbers used by unmapped chromosomes
    # These chromosomes keep their original names, so we can't assign these numbers
    reserved_numbers = set()
    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_chromosome_prefix)
        if is_autosome_suffix(suffix):
            reserved_numbers.add(int(suffix))
            print(f"  Reserved: number {suffix} (unmapped {name} keeps original name)")
    
    # Track which target numbers are reserved by sex chromosomes in reference
    # (when sex -> autosome, we skip that autosome number)
    skipped_numbers = set()
    
    # Process sex chromosomes first - they keep their original suffix
    for m in sex_mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_chromosome_prefix)
        target_suffix = m.target_suffix
        
        # If sex chromosome maps to autosome target, mark that number as skipped
        if is_autosome_suffix(target_suffix):
            skipped_numbers.add(int(target_suffix))
            print(f"  Note: {m.query_name} (sex chr) -> {m.target_name} (autosome): "
                  f"number {target_suffix} will be skipped")
        
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{query_suffix}",  # Keep original suffix
            new_suffix=query_suffix,
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=True
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement
    
    # Build mapping: target number -> list of (mapping, alignment_length)
    target_to_autosomes = defaultdict(list)
    autosomes_to_sex_target = []  # Autosomes mapping to sex chromosome targets
    
    for m in autosome_mappings:
        target_suffix = m.target_suffix
        
        if is_sex_chromosome_suffix(target_suffix):
            # Autosome maps to sex chromosome target - will get new number
            autosomes_to_sex_target.append(m)
            print(f"  Note: {m.query_name} (autosome) -> {m.target_name} (sex chr): "
                  f"will be reassigned")
        elif is_autosome_suffix(target_suffix):
            target_num = int(target_suffix)
            target_to_autosomes[target_num].append((m, m.total_alignment_length))
    
    # Resolve conflicts: for each target, pick the one with most alignment
    assigned_numbers = set()
    autosome_assignments = []  # (mapping, assigned_number)
    deferred_autosomes = []  # Mappings that lost in conflict resolution
    
    for target_num, candidates in target_to_autosomes.items():
        # Skip this number if it's reserved by an unmapped chromosome
        if target_num in reserved_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            print(f"  Number {target_num} reserved - deferring: {[m.query_name for m, _ in candidates]}")
            continue
            
        # Skip if blocked by sex chromosome mapping
        if target_num in skipped_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            continue
            
        # Sort by alignment length, descending
        candidates.sort(key=lambda x: x[1], reverse=True)
        
        # Winner gets the number
        winner, _ = candidates[0]
        autosome_assignments.append((winner, target_num))
        assigned_numbers.add(target_num)
        
        # Losers are deferred
        for m, _ in candidates[1:]:
            deferred_autosomes.append(m)
            print(f"  Conflict: {m.query_name} lost to {winner.query_name} for number {target_num}")
    
    # Find max assigned autosome number (including reserved)
    all_used_numbers = assigned_numbers | reserved_numbers
    max_autosome = max(all_used_numbers) if all_used_numbers else 0
    
    # Assign deferred autosomes (from conflicts) to next available numbers
    next_available = max_autosome + 1
    unavailable = assigned_numbers | reserved_numbers | skipped_numbers
    
    for m in deferred_autosomes:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        print(f"  Reassigned: {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1
    
    # Assign autosomes that mapped to sex chromosome targets
    for m in autosomes_to_sex_target:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        print(f"  Reassigned (was sex target): {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1
    
    # Create final assignments for autosomes
    for m, num in autosome_assignments:
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{num}",
            new_suffix=str(num),
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=False
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement
    
    # Handle unmapped chromosomes - they keep original suffix with new prefix
    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_chromosome_prefix)
        assignment = FinalChromosomeAssignment(
            original_name=name,
            new_name=f"{output_prefix}{suffix}",
            new_suffix=suffix,
            needs_reverse_complement=False,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix)
        )
        assignments.append(assignment)
        rc_lookup[name] = False
        print(f"  Unmapped: {name} -> {output_prefix}{suffix} (keeping orientation)")
    
    return assignments, rc_lookup


def sort_assignments_for_output(
    assignments: List[FinalChromosomeAssignment]
) -> List[FinalChromosomeAssignment]:
    """
    Sort assignments for output: autosomes by number, then sex chromosomes alphabetically.
    
    Args:
        assignments: List of FinalChromosomeAssignment
        
    Returns:
        Sorted list
    """
    autosomes = [a for a in assignments if not a.is_sex_chromosome]
    sex_chrs = [a for a in assignments if a.is_sex_chromosome]
    
    # Sort autosomes by numeric suffix
    autosomes.sort(key=lambda x: int(x.new_suffix))
    
    # Sort sex chromosomes alphabetically
    sex_chrs.sort(key=lambda x: x.new_suffix)
    
    return autosomes + sex_chrs


def write_fasta(
    sequences: Dict[str, str],
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_",
    line_width: int = 60
) -> None:
    """
    Write FASTA file with renamed chromosomes and reverse complement where needed.
    
    Orientation is taken from a.needs_reverse_complement (chromosomes) and
    unloc.needs_reverse_complement (unlocalized contigs).

    Args:
        sequences: Original sequences dictionary
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output FASTA file
        output_prefix: Prefix for output chromosome names
        line_width: Line width for sequence output (default 60)
    """
    # Build unloc parent lookup: parent_chr_original -> parent_chr_new_suffix
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    
    # Track which sequences have been processed
    processed_sequences = set()
    
    # Build per-parent unloc lookup: original_parent_name -> sorted list of UnlocMapping
    unloc_by_parent: Dict[str, List[UnlocMapping]] = defaultdict(list)
    for unloc in unloc_mappings:
        unloc_by_parent[unloc.parent_chromosome].append(unloc)
    for parent in unloc_by_parent:
        unloc_by_parent[parent].sort(key=lambda x: x.unloc_number)

    with open(output_path, 'w') as f:
        # Write chromosomes in sorted order, each followed by its unloc contigs
        for a in assignments:
            seq = sequences.get(a.original_name, '')
            if not seq:
                print(f"  Warning: No sequence found for {a.original_name}")
                continue
            
            processed_sequences.add(a.original_name)
            
            # Apply reverse complement if needed
            if a.needs_reverse_complement:
                seq = reverse_complement(seq)
            
            # Write header and sequence
            f.write(f">{a.new_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
            
            # Write unloc contigs belonging to this parent immediately after
            for unloc in unloc_by_parent.get(a.original_name, []):
                unloc_seq = sequences.get(unloc.contig_name, '')
                if not unloc_seq:
                    print(f"  Warning: No sequence found for {unloc.contig_name}")
                    continue
                
                processed_sequences.add(unloc.contig_name)
                
                # Apply reverse complement if needed (inherited from parent)
                if unloc.needs_reverse_complement:
                    unloc_seq = reverse_complement(unloc_seq)
                
                # Generate new name: <output_prefix><new_parent_suffix>_unloc_<N>
                new_parent_suffix = parent_suffix_lookup.get(unloc.parent_chromosome,
                                                             extract_chromosome_suffix(unloc.parent_chromosome, output_prefix))
                new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"
                
                f.write(f">{new_name}\n")
                for i in range(0, len(unloc_seq), line_width):
                    f.write(unloc_seq[i:i+line_width] + '\n')
        
        # Write remaining sequences (non-SUPER_ contigs) with original names and orientation
        for seq_name, seq in sequences.items():
            if seq_name in processed_sequences:
                continue
            
            # Keep original name and orientation
            f.write(f">{seq_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
    
    print(f"FASTA written to: {output_path}")


def write_chromosome_list(
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_"
) -> None:
    """
    Write chromosome list CSV file.
    
    Format: name,suffix,yes/no (yes for main chr, no for unloc)
    
    Args:
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output CSV file
        output_prefix: Prefix for output chromosome names
    """
    # Build parent suffix lookup
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    
    # Build per-parent unloc lookup
    unloc_by_parent_csv: Dict[str, List[UnlocMapping]] = defaultdict(list)
    for unloc in unloc_mappings:
        unloc_by_parent_csv[unloc.parent_chromosome].append(unloc)
    for parent in unloc_by_parent_csv:
        unloc_by_parent_csv[parent].sort(key=lambda x: x.unloc_number)

    with open(output_path, 'w') as f:
        # Write main chromosomes each followed by their unloc contigs
        for a in assignments:
            f.write(f"{a.new_name},{a.new_suffix},yes\n")
            
            for unloc in unloc_by_parent_csv.get(a.original_name, []):
                new_parent_suffix = parent_suffix_lookup.get(unloc.parent_chromosome,
                                                             extract_chromosome_suffix(unloc.parent_chromosome, output_prefix))
                new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"
                f.write(f"{new_name},{new_parent_suffix},no\n")
    
    print(f"Chromosome list written to: {output_path}")


def plot_chromosome_alignments(
    records: List[PAFRecord],
    mappings: List[ChromosomeMapping],
    output_dir: Path,
) -> None:
    """
    Generate scatter plots of PAF alignment blocks for each mapped chromosome.

    Left panel:  raw alignment blocks (query_mid vs target_mid), coloured by strand.
    Right panel: same data after applying the orientation decision (target axis
                 flipped when RC is needed), so a correct decision shows a clean
                 diagonal from bottom-left to top-right.

    Args:
        records:    Filtered PAF records.
        mappings:   Chromosome mappings with orientation decisions.
        output_dir: Directory where PNG files are written.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "  Warning: matplotlib not installed — skipping plots. "
            "Install with: pip install matplotlib",
            file=sys.stderr,
        )
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Index records by query→target
    blocks_index: Dict[str, List[PAFRecord]] = defaultdict(list)
    for r in records:
        blocks_index[f"{r.query_name}->{r.target_name}"].append(r)

    for mapping in mappings:
        if is_unloc_contig(mapping.query_name):
            continue
        key = f"{mapping.query_name}->{mapping.target_name}"
        recs = blocks_index.get(key, [])
        if not recs:
            continue

        plus_recs  = [r for r in recs if r.strand == "+"]
        minus_recs = [r for r in recs if r.strand == "-"]
        target_len = max(r.target_end for r in recs)
        needs_rc   = mapping.needs_reverse_complement

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        # ---- Left: raw scatter ----
        ax = axes[0]
        if plus_recs:
            ax.scatter(
                [(r.target_start + r.target_end) / 2e6 for r in plus_recs],
                [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                s=1, alpha=0.4, color="steelblue", label="+ strand",
            )
        if minus_recs:
            ax.scatter(
                [(r.target_start + r.target_end) / 2e6 for r in minus_recs],
                [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                s=1, alpha=0.4, color="firebrick", label="- strand",
            )
        ax.set_xlabel(f"{mapping.target_name} position (Mb)")
        ax.set_ylabel(f"{mapping.query_name} position (Mb)")
        ax.set_title("All alignment blocks")
        ax.legend(markerscale=8, loc="upper left")

        # ---- Right: after applying orientation decision ----
        ax2 = axes[1]
        if needs_rc:
            # Flip target axis so the signal reads left→right
            if plus_recs:
                ax2.scatter(
                    [(target_len - (r.target_start + r.target_end) / 2) / 1e6 for r in plus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                    s=1, alpha=0.4, color="steelblue", label="+ strand",
                )
            if minus_recs:
                ax2.scatter(
                    [(target_len - (r.target_start + r.target_end) / 2) / 1e6 for r in minus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                    s=1, alpha=0.4, color="firebrick", label="- strand",
                )
            ax2.set_title("After RC (target axis flipped)")
        else:
            if plus_recs:
                ax2.scatter(
                    [(r.target_start + r.target_end) / 2e6 for r in plus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in plus_recs],
                    s=1, alpha=0.4, color="steelblue", label="+ strand",
                )
            if minus_recs:
                ax2.scatter(
                    [(r.target_start + r.target_end) / 2e6 for r in minus_recs],
                    [(r.query_start + r.query_end) / 2e6 for r in minus_recs],
                    s=1, alpha=0.4, color="firebrick", label="- strand",
                )
            ax2.set_title("As-is (no RC applied)")

        ax2.set_xlabel(f"{mapping.target_name} position (Mb)")
        ax2.set_ylabel(f"{mapping.query_name} position (Mb)")
        ax2.legend(markerscale=8, loc="upper left")

        rc_str  = "RC" if needs_rc else "no RC"
        cov_str = f"{mapping.coverage:.1%}"
        fig.suptitle(
            f"{mapping.query_name}  →  {mapping.target_name}   |   "
            f"decision: {rc_str}   coverage: {cov_str}   "
            f"+: {mapping.plus_strand_length:,}  −: {mapping.minus_strand_length:,}",
            fontsize=11, fontweight="bold",
        )
        fig.tight_layout()
        out_path = plots_dir / f"{mapping.query_name}_vs_{mapping.target_name}.png"
        fig.savefig(out_path, dpi=120)
        plt.close(fig)

    print(f"  Plots saved to: {plots_dir}")


def load_mapping_table_assignments(
    mapping_table_path: Path,
    sequences: Dict[str, str],
    query_chromosome_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_",
) -> Tuple[List[FinalChromosomeAssignment], List[UnlocMapping]]:
    """
    Build assignments and unloc mappings from a pre-built mapping TSV
    (produced by a previous run on haplotype 1). Columns used: query,
    renamed_to, needs_reverse_complement. Unlocs present only in this
    haplotype are mapped via their parent chromosome; they are never RC'd.
    """
    # Parse table: query -> {renamed_to, needs_rc}
    table_map: Dict[str, Dict] = {}
    with open(mapping_table_path) as fh:
        header = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                missing = {"query", "renamed_to", "needs_reverse_complement"} - set(header)
                if missing:
                    raise ValueError(f"Mapping table missing columns: {', '.join(sorted(missing))}")
                continue
            row = dict(zip(header, fields))
            table_map[row["query"]] = {
                "renamed_to": row["renamed_to"],
                "needs_rc": row["needs_reverse_complement"].strip().lower() == "yes",
            }
    if not table_map:
        raise ValueError(f"Mapping table {mapping_table_path} is empty.")
    print(f"  Loaded {len(table_map)} entries from mapping table")

    def _suffix_from_renamed(renamed_to: str) -> str:
        """Extract the bare chromosomal suffix (e.g. '1', 'X') from renamed_to
        regardless of which prefix was used in the original run."""
        m = re.search(r'([A-Z]\d*|\d+)$', renamed_to, re.IGNORECASE)
        return m.group(1) if m else renamed_to

    # parent query name -> new suffix (main chrs only, for unloc naming)
    parent_new_suffix = {
        q: _suffix_from_renamed(info["renamed_to"])
        for q, info in table_map.items() if not is_unloc_contig(q)
    }

    # Build assignments for main chromosomes in FASTA
    assignments: List[FinalChromosomeAssignment] = []
    for orig in (n for n in sequences if n.startswith(query_chromosome_prefix) and not is_unloc_contig(n)):
        if orig not in table_map:
            suffix = extract_chromosome_suffix(orig, query_chromosome_prefix)
            print(f"  Warning: {orig} not in mapping table — keeping original suffix")
        else:
            suffix = _suffix_from_renamed(table_map[orig]["renamed_to"])
        needs_rc = table_map[orig]["needs_rc"] if orig in table_map else False
        assignments.append(FinalChromosomeAssignment(
            original_name=orig,
            new_name=f"{output_prefix}{suffix}",
            new_suffix=suffix,
            needs_reverse_complement=needs_rc,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix),
        ))

    # Build unloc mappings (never RC'd)
    unloc_mappings: List[UnlocMapping] = []
    for contig in (n for n in sequences if n.startswith(query_chromosome_prefix) and is_unloc_contig(n)):
        parent, unloc_num = parse_unloc_name(contig)
        if parent not in parent_new_suffix:
            print(f"  Warning: parent '{parent}' for '{contig}' not in mapping table — skipping")
            continue
        unloc_mappings.append(UnlocMapping(
            contig_name=contig, parent_chromosome=parent,
            unloc_number=unloc_num, needs_reverse_complement=False,
        ))

    print(f"  Built {len(assignments)} assignments and {len(unloc_mappings)} unloc mappings")
    return assignments, unloc_mappings


def _write_and_validate(
    args: argparse.Namespace,
    sequences: Dict[str, str],
    sorted_assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_chromosome_prefix: str,
    mapping_tsv_path: Path = None,
) -> None:
    """Print assignment table, write FASTA + CSV, validate genome length."""
    print("\nFinal chromosome assignments:")
    print("-" * 70)
    print(f"{'Original':<20} {'New Name':<20} {'Suffix':<10} {'RC?':>5} {'Sex?':>5}")
    print("-" * 70)
    for a in sorted_assignments:
        print(f"{a.original_name:<20} {a.new_name:<20} {a.new_suffix:<10} "
              f"{'Yes' if a.needs_reverse_complement else 'No':>5} "
              f"{'Yes' if a.is_sex_chromosome else 'No':>5}")

    fasta_out = args.output_dir / f"{args.output_prefix}.fa"
    csv_out   = args.output_dir / f"{args.output_prefix}.chromosome.list.csv"

    print("\nWriting output files...")
    write_fasta(sequences, sorted_assignments, unloc_mappings, fasta_out, output_chromosome_prefix)
    write_chromosome_list(sorted_assignments, unloc_mappings, csv_out, output_chromosome_prefix)

    print("\nValidating genome length...")
    in_len  = sum(len(s) for s in sequences.values())
    out_len = calculate_genome_length(fasta_out)
    if in_len == out_len:
        print(f"  OK: Genome length matches ({in_len:,} bp)")
    else:
        print(f"  ERROR: Genome length mismatch! Input: {in_len:,} bp, Output: {out_len:,} bp", file=sys.stderr)
        sys.exit(1)

    print(f"\nDone! Output files:")
    print(f"  - FASTA: {fasta_out}")
    print(f"  - Chromosome list: {csv_out}")
    if mapping_tsv_path:
        print(f"  - Mapping summary: {mapping_tsv_path}")


def main():
    """Main entry point."""
    args = parse_args()

    qpfx = args.query_chromosome_prefix
    opfx = args.output_chromosome_prefix

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading FASTA file: {args.fasta}")
    sequences = read_fasta(args.fasta)
    print(f"  Found {len(sequences)} sequences")

    # --- Mode A: mapping table (second haplotype, no alignment needed) ---
    if args.mapping_table is not None:
        if not args.mapping_table.exists():
            print(f"Error: Mapping table not found: {args.mapping_table}", file=sys.stderr)
            sys.exit(1)
        # Warn about PAF-only options that have no effect in this mode
        if args.min_coverage != 0.5:
            print("  Warning: --min-coverage is ignored in --mapping-table mode", file=sys.stderr)
        if args.plot_alignments:
            print("  Warning: --plot-alignments is ignored in --mapping-table mode", file=sys.stderr)
        if args.reference_chromosome_prefix is not None:
            print("  Warning: --reference-chromosome-prefix is ignored in --mapping-table mode", file=sys.stderr)
        print(f"\nUsing pre-built mapping table: {args.mapping_table}")
        print(f"  Input prefix: '{qpfx}' -> Output prefix: '{opfx}'")
        assignments, unloc_mappings = load_mapping_table_assignments(
            args.mapping_table, sequences, qpfx, opfx
        )
        _write_and_validate(args, sequences, sort_assignments_for_output(assignments),
                            unloc_mappings, opfx)
        return

    # --- Mode B: full PAF-based alignment ---
    print(f"Parsing PAF file: {args.paf}")
    paf_records = parse_paf(args.paf)
    print(f"  Found {len(paf_records)} alignment records")

    ref_prefix_arg = args.reference_chromosome_prefix
    filtered_records, ref_prefix = filter_paf_records(paf_records, qpfx, ref_prefix_arg)
    label = "provided" if ref_prefix_arg else "auto-detected"
    print(f"  Reference prefix ({label}): '{ref_prefix}'")
    print(f"  After filtering ({qpfx}* -> {ref_prefix}*): {len(filtered_records)} records")

    print("\nValidating PAF/FASTA consistency...")
    is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
        paf_records, sequences, qpfx
    )
    if in_paf_not_fasta:
        print(f"  WARNING: In PAF but not FASTA: {', '.join(in_paf_not_fasta)}")
    if in_fasta_not_paf:
        print(f"  WARNING: In FASTA but not PAF: {', '.join(in_fasta_not_paf)} (will keep original suffix)")
    if is_valid:
        print(f"  OK: All {qpfx}* chromosomes match")

    print(f"\nBuilding chromosome mappings (min coverage: {args.min_coverage:.0%})...")
    mappings = build_chromosome_mappings(filtered_records, args.min_coverage, ref_prefix)
    print(f"  Successfully mapped {len(mappings)} chromosomes")

    if args.plot_alignments:
        print("\nGenerating alignment scatter plots...")
        plot_chromosome_alignments(filtered_records, mappings, args.output_dir)

    print_mapping_summary(mappings)

    print("\nResolving chromosome assignments...")
    print(f"  Input prefix: '{qpfx}' -> Output prefix: '{opfx}'")
    assignments, _ = resolve_chromosome_assignments(mappings, sequences, qpfx, opfx)

    unloc_mappings = build_unloc_mappings(sequences, mappings)
    if unloc_mappings:
        print(f"  Found {len(unloc_mappings)} unlocalized contigs")
    for unloc in unloc_mappings:
        unloc.needs_reverse_complement = False

    sorted_assignments = sort_assignments_for_output(assignments)
    mapping_tsv_path = args.output_dir / f"{args.output_prefix}.mapping.tsv"
    save_mapping_tsv(mappings, sorted_assignments, mapping_tsv_path)

    _write_and_validate(args, sequences, sorted_assignments, unloc_mappings, opfx, mapping_tsv_path)


if __name__ == "__main__":
    main()
