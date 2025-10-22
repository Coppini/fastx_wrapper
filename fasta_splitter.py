#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
from sys import stderr
from typing import Any, Generator, Iterable

import tqdm

from fastx_io_wrapper import (
    FASTA_EXTENSIONS,
    FASTQ_EXTENSIONS,
    GZIP_EXTENSIONS,
    FastaFields,
    FastaGenerator,
    FastxFormat,
    FastxIOWrapper,
    PairedFastxIOWrapper,
    QualitySimulator,
)

COMPLEMENT_DICT = {
    "A": "T",
    "C": "G",
    "T": "A",
    "G": "C",
    "Y": "R",
    "R": "Y",
    "W": "W",
    "S": "S",
    "K": "M",
    "M": "K",
    "D": "H",
    "V": "B",
    "H": "D",
    "B": "V",
    "N": "N",
    "X": "N",
    "-": "-",
}
for key, value in list(COMPLEMENT_DICT.items()):
    COMPLEMENT_DICT[key.lower()] = value.lower()


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT_DICT[base] for base in reversed(seq))


DEFAULT_WINDOW_SIZE = 150
DEFAULT_COVERAGE = 1.0
DEFAULT_QUALITY_SIMULATOR = QualitySimulator()

PairedFastaGenerator = Generator[
    tuple[FastaFields, FastaFields], None, None
]  # yields ((header, sequence), (header, sequence))


def split_sequence_by_coverage_into_single(
    sequence: str,
    seq_number: int = 0,
    header: str = "",
    window_size: int = 150,
    coverage: float = 1.0,
) -> FastaGenerator:
    seq_len = len(sequence)
    window_size = min(seq_len, window_size)
    n_windows = math.ceil((seq_len * coverage) / window_size)

    if n_windows * window_size >= seq_len:
        full_tiling = True
        step = (seq_len - window_size) / max(1, n_windows - 1)
    else:
        full_tiling = False
        step = (seq_len - (n_windows * window_size)) / (n_windows + 1)

    for i in range(n_windows):
        # start = int(round(((i + int(not full_tiling)) * step) + (window_size * int(not full_tiling) * i)))
        start = int(
            round((i * step) if full_tiling else ((i + 1) * step + (window_size * i)))
        )
        end = start + window_size
        yield f"{seq_number}_{i}_{header}:{start}-{end}", sequence[start:end]
        if end > seq_len:
            raise RuntimeError(
                f"Calculated end position {end} is greater than sequence length {seq_len} for sequence {header}"
            )
    if full_tiling and end < seq_len:
        raise RuntimeError(
            f"Did not reach the end of the sequence for sequence {header}, the final position was {end}/{seq_len}"
        )


def split_sequence_by_coverage_into_pairs(
    sequence: str,
    seq_number: int = 0,
    header: str = "",
    window_size: int = 150,
    coverage: float = 1.0,
    gap_between_pair: int | None = None,  # negative for overlapping paired-ends
) -> PairedFastaGenerator:
    seq_len = len(sequence)
    window_size = min(seq_len, window_size)
    pair_bases = window_size * 2
    n_pairs = math.ceil((seq_len * coverage) / pair_bases)

    if gap_between_pair is None:
        gap_between_pair = window_size
    max_gap = int(max(0, (seq_len / 2) - window_size - coverage))
    gap_between_pair = min(gap_between_pair, max_gap)
    insert_size = min(seq_len - gap_between_pair, pair_bases + gap_between_pair)
    gap_between_pair = insert_size - pair_bases  # adjust gap if near the end

    if n_pairs * pair_bases >= seq_len:
        full_tiling = True
        step = (seq_len - insert_size) / max(1, n_pairs - 1)
    else:
        full_tiling = False
        step = (seq_len - (n_pairs * insert_size)) / (n_pairs + 1)

    for i in range(n_pairs):
        # start = int(round(((i + int(not full_tiling)) * step) + (window_size * int(not full_tiling) * i)))
        paired_start = int(
            round((i * step) if full_tiling else ((i + 1) * step + (window_size * i)))
        )
        paired_end = paired_start + insert_size
        paired_header = f"{seq_number}_{i}_{header}:{paired_start}-{paired_end}"
        # paired_sequence = sequence[paired_start:paired_end]
        # read1_sequence = paired_sequence[:window_size]
        # read2_sequence = reverse_complement(paired_sequence[-window_size:])
        read1_fields = (
            f"{paired_header}/1",
            sequence[paired_start : paired_start + window_size],
        )
        read2_fields = (
            f"{paired_header}/2",
            reverse_complement(sequence[paired_end - window_size : paired_end]),
        )
        if paired_end > seq_len:
            raise RuntimeError(
                f"Calculated end position {paired_end} is greater than sequence length {seq_len} for sequence {header}"
            )
        yield read1_fields, read2_fields
    if full_tiling and paired_end < seq_len:
        raise RuntimeError(
            f"Did not reach the end of the sequence for sequence {header}, the final position was {paired_end}/{seq_len}"
        )


def positive_int_window_size(window_size: str) -> int:
    s_window_size = str(window_size).strip().lower()
    i_window_size = int(s_window_size.replace("bp", ""))
    if i_window_size <= 0:
        raise argparse.ArgumentTypeError(
            "Window-size must be an integer greater than 0."
        )
    return i_window_size


def positive_float_coverage(coverage: str) -> float:
    s_coverage = str(coverage).strip().lower().rstrip("x")
    f_coverage = (
        (float(s_coverage.rstrip("%")) / 100.0)
        if s_coverage.endswith("%")
        else float(s_coverage)
    )
    if f_coverage <= 0:
        raise argparse.ArgumentTypeError(
            "Coverage must be a float greater than 0."
            " If you want to skip bases, you can use a decimal coverage,"
            " with 'coverage < 0', but negative coverage does not exist."
        )
    return f_coverage


def tqdm_iter(items: Iterable, verbose: bool = False, **kwargs) -> Iterable[Any]:
    return tqdm.tqdm(items, **kwargs) if verbose else items


def split_fasta_into_reads(
    inp_io: FastxIOWrapper,
    window_size: int = 150,
    coverage: float = 1.0,
    paired: bool = False,
    gap_between_pair: int | None = None,  # negative for overlapping paired-ends
) -> PairedFastaGenerator:
    splitter_function = (
        split_sequence_by_coverage_into_pairs
        if paired
        else split_sequence_by_coverage_into_single
    )
    settings = {"window_size": window_size, "coverage": coverage}
    if paired:
        settings["gap_between_pair"] = gap_between_pair
    with inp_io as inp_fh:
        for i, (header, sequence) in enumerate(inp_fh.iter_sequences()):
            yield from splitter_function(
                sequence=sequence,
                seq_number=i,
                header=header,
                **settings,
            )


def split_from_fasta_and_write(
    inp_io: FastxIOWrapper,
    out_io: FastxIOWrapper | PairedFastxIOWrapper,
    window_size: int = 150,
    coverage: float = 1.0,
    paired: bool = False,
    quality_simulator: QualitySimulator | None = None,
    gap_between_pair: int | None = None,  # negative for overlapping paired-ends
    verbose: bool = True,
) -> tuple[int, int]:
    with out_io as out_fh:
        results = out_fh.write_fastx(
            sequences=tqdm_iter(
                split_fasta_into_reads(
                    inp_io=inp_io,
                    window_size=window_size,
                    coverage=coverage,
                    paired=paired,
                    gap_between_pair=gap_between_pair,
                ),
                verbose=verbose,
                desc="Splitting sequences",
            ),
            quality_simulator=quality_simulator,
        )
        out_fh.flush()
        return results


def determine_output_format(args: argparse.Namespace) -> FastxFormat | None:
    if args.fastq_mode:
        return FastxFormat.FASTQ
    if args.fasta_mode:
        return FastxFormat.FASTA
    if not (args.output or args.paired_output) or str(args.output) == "-":
        return FastxFormat.FASTA

    output_path = args.paired_output[0] if args.paired_output else args.output
    out_suffix = output_path.suffix.lower()
    if out_suffix in FASTQ_EXTENSIONS:
        return FastxFormat.FASTQ
    if out_suffix in GZIP_EXTENSIONS:
        out_suffix = output_path.with_suffix("").suffix.lower()
        if out_suffix in FASTQ_EXTENSIONS:
            return FastxFormat.FASTQ
        elif out_suffix in FASTA_EXTENSIONS:
            return FastxFormat.FASTA
    return None


if __name__ == "__main__":

    def _require_together(
        args, parser: argparse.ArgumentParser, required: str, dependent: str
    ) -> None:
        if getattr(args, dependent) and getattr(args, required) is None:
            parser.error(
                f"--{dependent.replace('_', '-')} requires --{required.replace('_', '-')}"
            )

    parser = argparse.ArgumentParser(
        description="Split FASTA sequences into smaller sequences."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input FASTA file",
    )

    parser.add_argument(
        "-w",
        "--window-size",
        "--length",
        "-l",
        default=DEFAULT_WINDOW_SIZE,
        type=positive_int_window_size,
        help=f"Desired window_size for split sequences (must be > 0) [default: {DEFAULT_WINDOW_SIZE}]",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        "--min-coverage",
        type=positive_float_coverage,
        default=DEFAULT_COVERAGE,
        help=f"Minimum coverage desired for the sequences [default: {DEFAULT_COVERAGE}]",
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "--verbose", action="store_true", help="Show progress bar"
    )
    verbosity_group.add_argument(
        "--quiet", action="store_true", help="Hide progress bar"
    )

    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output FASTA file (if not provided, will print to STDOUT)",
    )
    output_group.add_argument(
        "-p",
        "--paired-output",
        "--paired",
        type=Path,
        nargs=2,
        help="If provided, generates two files with paired sequences (e.g. for paired-end reads). Two filepaths must be provided.",
    )

    fast_group = parser.add_mutually_exclusive_group()
    fast_group.add_argument(
        "-f",
        "--fastq-mode",
        "--fastq",
        action="store_true",
        help="If provided, outputs in FASTQ format with dummy quality scores (by default, outputs FASTA unless output filename ends with .fq/.fastq).",
    )
    fast_group.add_argument(
        "-a",
        "--fasta-mode",
        "--fasta",
        action="store_true",
        help="If provided, outputs in FASTA format (by default, outputs FASTA unless output filename ends with .fq/.fastq).",
    )

    parser.add_argument(
        "-g",
        "--gap-between-pair",
        type=int,
        default=None,
        help="Gap between the two reads when using paired mode. If negative, reads will overlap. [default: same as window_size]",
    )

    quality_group = parser.add_argument_group(
        title="Quality score simulation options",
        description="Options for simulating quality scores when outputting FASTQ files.",
    )
    quality_group.add_argument(
        "--mean-q",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.mean_q,
        help=f"Mean quality score for simulated FASTQ [default: {DEFAULT_QUALITY_SIMULATOR.mean_q}]",
    )
    quality_group.add_argument(
        "--stdev-q",
        type=float,
        default=DEFAULT_QUALITY_SIMULATOR.stdev_q,
        help=f"Standard deviation of quality scores [default: {DEFAULT_QUALITY_SIMULATOR.stdev_q}]",
    )
    quality_group.add_argument(
        "--min-q",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.min_q,
        help=f"Minimum quality score [default: {DEFAULT_QUALITY_SIMULATOR.min_q}]",
    )
    quality_group.add_argument(
        "--max-q",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.max_q,
        help=f"Maximum quality score [default: {DEFAULT_QUALITY_SIMULATOR.max_q}]",
    )
    quality_group.add_argument(
        "--max-start-decay",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.max_start_decay,
        help=f"Maximum decay at start of read [default: {DEFAULT_QUALITY_SIMULATOR.max_start_decay}]",
    )
    quality_group.add_argument(
        "--max-end-decay",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.max_end_decay,
        help=f"Maximum decay at end of read [default: {DEFAULT_QUALITY_SIMULATOR.max_end_decay}]",
    )
    quality_group.add_argument(
        "--phred-offset",
        type=int,
        default=DEFAULT_QUALITY_SIMULATOR.phred_offset,
        help=f"PHRED offset for quality scores [default: {DEFAULT_QUALITY_SIMULATOR.phred_offset}]",
    )

    args = parser.parse_args()

    _require_together(
        args, parser, required="paired_output", dependent="gap_between_pair"
    )

    window_size = args.window_size
    coverage = args.coverage
    verbose = args.verbose or (not args.quiet and str(args.output) != "-")

    paired = bool(args.paired_output)
    if args.gap_between_pair is not None and not paired:
        print(
            "WARNING: --gap-between-pair is ignored when not using paired mode.",
            file=stderr,
        )

    inp_io = FastxIOWrapper(file=args.input, mode="rt")
    if inp_io.fastx_format is None:
        inp_io.fastx_format = FastxFormat.FASTA
    output_fastx_format = determine_output_format(args)
    out_io = (
        PairedFastxIOWrapper(
            read1_io=FastxIOWrapper(
                file=args.paired_output[0], mode="wt", fastx_format=output_fastx_format
            ),
            read2_io=FastxIOWrapper(
                file=args.paired_output[1], mode="wt", fastx_format=output_fastx_format
            ),
        )
        if paired
        else FastxIOWrapper(
            file=args.output, mode="wt", fastx_format=output_fastx_format
        )
    )
    quality_simulator = QualitySimulator(
        mean_q=args.mean_q,
        stdev_q=args.stdev_q,
        min_q=args.min_q,
        max_q=args.max_q,
        max_start_decay=args.max_start_decay,
        max_end_decay=args.max_end_decay,
        phred_offset=args.phred_offset,
    )

    split_from_fasta_and_write(
        inp_io=inp_io,
        out_io=out_io,
        window_size=window_size,
        coverage=coverage,
        paired=paired,
        quality_simulator=quality_simulator,
        gap_between_pair=args.gap_between_pair,
        verbose=verbose,
    )
