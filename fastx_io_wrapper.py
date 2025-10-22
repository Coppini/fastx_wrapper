#!/usr/bin/env python3

import gzip
from enum import Enum
from itertools import islice, tee, zip_longest
from pathlib import Path
from random import gauss
from sys import stderr, stdin, stdout
from typing import IO, Any, Callable, Generator, Iterable

from Bio import bgzf
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

BGZIP_EXTENSIONS = {".bgz", ".bgzip"}
GZIP_EXTENSIONS = {".gz", ".gzip"} | BGZIP_EXTENSIONS

FASTA_ALIASES = {"fasta", "fa", "faa", "fna", "fsa", "fas", "ffn"}
FASTQ_ALIASES = {"fastq", "fq"}
FASTX_ALIASES = FASTA_ALIASES | FASTQ_ALIASES

FASTA_EXTENSIONS = {f".{alias}" for alias in FASTA_ALIASES}
FASTQ_EXTENSIONS = {f".{alias}" for alias in FASTQ_ALIASES}
FASTX_EXTENSIONS = FASTA_EXTENSIONS | FASTQ_EXTENSIONS

DEFAULT_ILLUMINA_PHRED_OFFSET = 33
MAX_ILLUMINA_QUALITY_CHAR = chr(DEFAULT_ILLUMINA_PHRED_OFFSET + 40)  # 'I'


class FastxFormat(Enum):
    FASTA = "fasta"
    FASTQ = "fastq"


def parse_fastx_format(value: str) -> FastxFormat:
    v = value.lower().strip()
    if v.startswith("."):
        v = v[1:]
    if v in (FASTA_ALIASES | {"a"}):
        return FastxFormat.FASTA
    elif v in (FASTQ_ALIASES | {"q"}):
        return FastxFormat.FASTQ
    else:
        raise ValueError(
            f"Unknown format '{value}' is not recognized as neither FASTA nor FASTQ."
        )


class UnknownFastxFormatError(ValueError):
    pass


class IterableMismatchedCounts(ValueError):
    def __init__(
        self,
        message: str = "",
        missing_data: str = "<missing>",
        index: int = 0,
        function_name: str = "",
    ) -> None:
        self.function_name = function_name
        self.missing_data = missing_data  # "read1" or "read2"
        self.index = index  # 1-based index where mismatch occurred
        self.message = message or (
            f"{function_name} mismatched count of records in paired iterables - "
            f"{self.missing_data} missing data on record {self.index} (1-based)."
        )
        super().__init__(self.message)


def chunk(
    iterable: Iterable, size: int, raise_on_mismatch: bool = True
) -> Iterable[tuple]:
    """Yield tuples of `size` items from `iterator` until exhausted."""
    iterator = iter(iterable)
    index = 1
    while chunk := tuple(islice(iterator, size)):
        if raise_on_mismatch and any(item is None for item in chunk):
            if all(item is None for item in chunk):
                break
            raise IterableMismatchedCounts(
                index=index,
                function_name="chunk",
            )
        yield chunk
        index += 1


FastaFields = tuple[str, str]  # (header, sequence)
FastqFields = tuple[str, str, str]  # (header, sequence, quality)
FastxFields = FastaFields | FastqFields  # (header, sequence[, quality])


FastaGenerator = Generator[FastaFields, None, None]  # yields (header, sequence)
FastaIterator = Iterable[FastaFields]  # yields (header, sequence)

FastqGenerator = Generator[
    FastqFields, None, None
]  # yields (header, sequence, quality)
FastqIterator = Iterable[FastqFields]  # yields (header, sequence, quality)

FastxIterator = FastaIterator | FastqIterator  # yields (header, sequence[, quality])
FastxGenerator = FastaGenerator | FastqGenerator  # yields (header, sequence[, quality])

PairedFastxFields = tuple[
    FastxFields, FastxFields
]  # ((header, sequence[, quality]), (header, sequence[, quality]))
PairedFastqFields = tuple[
    FastqFields, FastqFields
]  # ((header, sequence, quality), (header, sequence, quality))
PairedFastaFields = tuple[
    FastaFields, FastaFields
]  # ((header, sequence), (header, sequence))

PairedFastqIterator = Iterable[PairedFastqFields]
PairedFastaIterator = Iterable[PairedFastaFields]
PairedFastxIterator = PairedFastaIterator | PairedFastqIterator

PairedFastaGenerator = Generator[PairedFastaFields, None, None]
PairedFastqGenerator = Generator[PairedFastqFields, None, None]
PairedFastxGenerator = PairedFastaGenerator | PairedFastqGenerator


class QualitySimulator:
    def __init__(
        self,
        mean_q: int = 38,
        stdev_q: float = 0.0,
        min_q: int = 38,
        max_q: int = 38,
        max_start_decay: int = 0,
        max_end_decay: int = 0,
        phred_offset: int | None = None,
    ):
        self.mean_q = mean_q
        self.stdev_q = stdev_q
        self.min_q = min_q
        self.max_q = max_q
        self.max_start_decay = max_start_decay
        self.max_end_decay = max_end_decay
        self.phred_offset = (
            phred_offset if phred_offset is not None else DEFAULT_ILLUMINA_PHRED_OFFSET
        )
        self.active = self.get_active()

        if self.min_q <= self.mean_q <= self.max_q:
            self.default_quality = self.mean_q
        elif self.min_q > self.mean_q:
            self.default_quality = self.min_q
            print(
                f"Warning: mean_q ({self.mean_q}) is lower than min_q ({self.min_q}). Default quality will be set to the min_q={self.min_q}",
                file=stderr,
            )
        elif self.max_q < self.mean_q:
            self.default_quality = self.max_q
            print(
                f"Warning: mean_q ({self.mean_q}) is greater than max_q ({self.max_q}). Default quality will be set to the max_q={self.max_q}",
                file=stderr,
            )
        else:  # if self.min_q > self.max_q:
            raise ValueError(
                f"max_q must be greater than or equal to min_q (max_q={self.max_q} < min_q={self.min_q})"
            )

        self.default_quality_char = chr(self.phred_offset + self.default_quality)

        self.simulate_scores = self.get_float_simulation_method()
        self.simulate_quality_string = self.get_string_simulation_method()

    def get_active(self):
        self.active = bool(
            (self.min_q != self.max_q)
            and (
                self.stdev_q
                or (
                    self.min_q < self.mean_q
                    and (self.max_end_decay or self.max_start_decay)
                )
            )
        )
        return self.active

    def get_float_simulation_method(
        self,
    ) -> Callable[[int], Generator[float, None, None]]:
        if not self.active:
            return self.mock_quality
        if self.stdev_q:
            if self.max_end_decay:
                if self.max_start_decay:
                    return self.apply_all
                return self.apply_std_dev_and_end_decay
            elif self.max_start_decay:
                return self.apply_std_dev_and_start_decay
            return self.apply_std_dev
        if self.max_end_decay:
            if self.max_start_decay:
                return self.apply_start_and_end_decay
            return self.apply_end_decay
        if self.max_start_decay:
            return self.apply_start_decay
        return self.mock_quality

    def get_string_simulation_method(self) -> Callable[[int], str]:
        if not self.active:
            return self.mock_quality_string
        return self.simulate_quality

    def apply_all(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string applying standard deviation and decay in both start and end."""
        yield from (
            gauss(self.mean_q, self.stdev_q)
            - (
                ((i / length) * self.max_end_decay)
                + ((length - i / length) * self.max_start_decay),
            )
            for i in range(length)
        )

    def apply_std_dev(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with standard deviation around the mean quality."""
        yield from (gauss(self.mean_q, self.stdev_q) for _ in range(length))

    def apply_std_dev_and_start_decay(
        self, length: int
    ) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with standard deviation and start decay."""
        yield from (
            gauss(self.mean_q, self.stdev_q)
            - ((length - i / length) * self.max_start_decay)
            for i in range(length)
        )

    def apply_std_dev_and_end_decay(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with standard deviation and end decay."""
        yield from (
            gauss(self.mean_q, self.stdev_q) - ((i / length) * self.max_end_decay)
            for i in range(length)
        )

    def apply_start_and_end_decay(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with dual decay from both ends."""
        yield from (
            self.mean_q
            - (
                ((i / length) * self.max_end_decay)
                + ((length - i / length) * self.max_start_decay),
            )
            for i in range(length)
        )

    def apply_start_decay(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with start decay."""
        yield from (
            gauss(self.mean_q, self.stdev_q)
            - ((length - i / length) * self.max_start_decay)
            for i in range(length)
        )

    def apply_end_decay(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with end decay."""
        yield from (
            gauss(self.mean_q, self.stdev_q) - ((i / length) * self.max_end_decay)
            for i in range(length)
        )

    def mock_quality(self, length: int) -> Generator[float, None, None]:
        """Generate a FASTQ quality string with constant quality."""
        yield from [self.default_quality] * length

    def mock_quality_string(self, length: int) -> str:
        """Generate a FASTQ quality string with constant quality."""
        return self.default_quality_char * length

    def simulate_quality(
        self,
        length: int,
    ) -> str:
        """Simulate a FASTQ quality string."""
        return "".join(
            chr(int(self.phred_offset + max(self.min_q, min(self.max_q, f_quality))))
            for f_quality in self.simulate_scores(length)
        )


def convert_sequence_to_fastq(
    fields: FastxFields, quality_simulator: QualitySimulator | None
) -> FastqFields:
    if len(fields) > 2:
        return fields[0], fields[1], fields[2]
    header, sequence = fields[0], fields[1]
    if not quality_simulator:
        return header, sequence, (len(sequence) * MAX_ILLUMINA_QUALITY_CHAR)
    return header, sequence, quality_simulator.simulate_quality_string(sequence)


def mock_qualities(
    sequences: Iterable[FastxFields], default_char: str
) -> FastqIterator:
    last_length = 0
    quality = ""
    for fields in sequences:
        if len(fields) > 2:
            yield fields[0], fields[1], fields[2]
        header, sequence = fields[0], fields[1]
        if (seq_len := len(sequence)) == last_length:
            yield header, sequence, quality
        yield header, sequence, (quality := (seq_len * default_char))


def convert_sequences_to_fastq(
    sequences: Iterable[FastxFields], quality_simulator: QualitySimulator | None
) -> FastqGenerator:
    if quality_simulator is None:
        yield from mock_qualities(sequences, MAX_ILLUMINA_QUALITY_CHAR)
    elif not quality_simulator.active:
        yield from mock_qualities(sequences, quality_simulator.default_quality_char)
    else:
        for fields in sequences:
            if len(fields) > 2:
                yield fields[0], fields[1], fields[2]
            header, sequence = fields[0], fields[1]
            yield (
                header,
                sequence,
                quality_simulator.simulate_quality_string(len(sequence)),
            )


def convert_sequence_pairs_to_fastq(
    sequence_pairs: PairedFastxIterator,
    quality_simulator: QualitySimulator | None,
) -> PairedFastqGenerator:
    r1_iter, r2_iter = tee(sequence_pairs)
    return (
        convert_sequences_to_fastq(
            (r1_seq for r1_seq, _ in r1_iter), quality_simulator=quality_simulator
        ),
        convert_sequences_to_fastq(
            (r2_seq for _, r2_seq in r2_iter), quality_simulator=quality_simulator
        ),
    )


def trim_sequence_to_length(fields: FastxFields, length: int) -> FastxFields:
    if len(fields) == 3:
        header, sequence, quality = fields
        return header, sequence[:length], quality[:length]
    elif len(fields) == 2:
        header, sequence = fields
        return header, sequence[:length]
    else:
        raise TypeError(
            f"Invalid number of fields in FastxFields: expected 2 or 3, got {len(fields)}"
        )


class FastxIOWrapper:
    def __init__(
        self,
        file: Path | str | None = None,
        handler: IO | Iterable | None = None,
        fastx_format: FastxFormat | None = None,
        mode: str = "r",
        mkdir: bool = False,
    ) -> None:
        self.file: Path | str | None = file
        self.handler: IO | Iterable | None = handler
        self.fastx_format: FastxFormat | None = fastx_format
        self.mode: str = mode
        self._prefetched_lines: list[str | bytes] = list()

        if self.file is not None and self.fastx_format is None:
            self.fastx_format = self.infer_format_from_filepath()

        if not isinstance(self.file, str) and str(self.file) == "-":
            self.file = "-"

        if isinstance(self.file, str) and not self.is_std_handle:
            self.file = Path(self.file)
        if (
            mkdir
            and isinstance(self.file, Path)
            and "w" in mode
            and not self.file.parent.exists()
        ):
            self.file.parent.mkdir(parents=True, exist_ok=True)

    def __repr__(self) -> str:
        return (
            f"FastxIOWrapper(file={self.file}, handler={self.handler}, "
            f"fastx_format={self.fastx_format}, mode='{self.mode}')"
        )

    @property
    def is_std_handle(self) -> bool:
        return (
            self.file == "-"
            or self.file is None
            or any(
                foh is value
                for foh in (self.file, self.handler)
                for value in (stdin, stdout, stderr, "-")
            )
            # or file is stdout
            # or file is stdin
            # or file is stderr
            # or handler is stdout
            # or handler is stdin
            # or handler is stderr
        )

    @property
    def source(self) -> str:
        if self.is_std_handle:
            return "stdout" if "w" in self.mode else "stdin"
        if self.file is not None:
            return str(self.file.name)
        if self.handler is not None:
            return str(self.handler)
        return "{unknown source}"

    def close(self) -> None:
        if self.is_std_handle:
            return
        return self.handler.close()

    @property
    def closed(self) -> bool:
        if self.is_std_handle:
            return False
        if self.handler is None:
            return True
        return self.handler.closed

    @property
    def readable(self) -> bool:
        if self.is_std_handle:
            return "r" in self.mode
        if self.handler is None:
            return False
        return hasattr(self.handler, "readable") and self.handler.readable()

    @property
    def writable(self) -> bool:
        if self.is_std_handle:
            return "w" in self.mode
        if self.handler is None:
            return False
        return hasattr(self.handler, "writable") and self.handler.writable()

    def flush(self) -> None:
        return self.handler.flush() if hasattr(self.handler, "flush") else None

    def __enter__(self):
        try:
            self.open()
        except Exception:
            self.close()
            raise
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> None:
        try:
            self.close()
        except Exception as close_exc:
            if exc_val is not None:
                raise close_exc from exc_val
            raise close_exc

    def open(self) -> "FastxIOWrapper":
        if self.handler is not None:
            if "w" in self.mode:
                if self.writable():
                    return self
                raise IOError(
                    f"File '{self.source}' is not writable and mode is '{self.mode}'"
                )
            if self.fastx_format is None:
                self.fastx_format = self.infer_format_from_handler()
            if self.readable():
                return self
            raise IOError(
                f"File '{self.source}' is not readable and mode is '{self.mode}'"
            )

        if self.handler is None and (str(self.file) == "-" or self.file is None):
            self.handler = stdout if "w" in self.mode else stdin
            self.fastx_format = self.infer_format_from_handler()
            return self

        path = Path(self.file)
        if not path.exists() and "w" not in self.mode:
            raise FileNotFoundError(
                f"File {self.file} does not exist for '{self.source}' with mode '{self.mode}'."
            )
        self.file = path
        if self.fastx_format is None:
            self.fastx_format = self.infer_format_from_filepath()
        suffix = self.file.suffix.lower()
        if suffix in GZIP_EXTENSIONS and "w" in self.mode:
            self.handler = bgzf.BgzfWriter(self.file, self.mode)
            return self
        opener = (
            bgzf.BgzfReader
            if suffix in BGZIP_EXTENSIONS
            else gzip.open
            if suffix in GZIP_EXTENSIONS
            else open
        )
        self.handler = opener(self.file, self.mode)
        if self.fastx_format is None:
            self.fastx_format = self.infer_format_from_handler()
        return self

    def infer_format_from_filepath(self) -> FastxFormat | None:
        if self.file is None:
            return
        suffix = self.file.suffix.lower()
        if suffix in GZIP_EXTENSIONS:
            suffix = self.file.with_suffix("").suffix.lower()
        if suffix in FASTQ_EXTENSIONS:
            return FastxFormat.FASTQ
        if suffix in FASTA_EXTENSIONS:
            return FastxFormat.FASTA

    def infer_format_from_handler(self) -> FastxFormat:
        if self.fastx_format is not None:
            return self.fastx_format
        if self.handler is None or not self.handler.readable() or self.handler.closed:
            raise UnknownFastxFormatError(
                f"Tried to infer fastX format for {self.file} from handler but could not read from handler {self.handler}."
            )

        first_line = self.handler.readline()
        first_line_first_40_chars = (
            first_line.replace("\r", "").replace("\n", "").strip()[:40]
        )
        if not first_line_first_40_chars:
            raise UnknownFastxFormatError(
                f"Could not determine fastX format from the handler for {self.source}."
                " Empty input stream."
            )
        if first_line_first_40_chars.startswith(">"):
            self.fastx_format = FastxFormat.FASTA
        if first_line_first_40_chars.startswith("@"):
            self.fastx_format = FastxFormat.FASTQ

        try:
            self.handler.seek(0)
        except (OSError, IOError):
            self._prefetched_lines.append(first_line)
        if self.fastx_format is None:
            raise UnknownFastxFormatError(
                f"Could not determine fastX format from the handler for {self.source}."
                f" First line does not start with '>' nor '@': {first_line_first_40_chars}"
            )
        return self.fastx_format

    def infer_format(self):
        if self.fastx_format is not None:
            return self.fastx_format
        if self.file is not None:
            self.fastx_format = self.infer_format_from_filepath()
            if self.fastx_format is not None:
                return self.fastx_format
        if self.handler is not None:
            self.fastx_format = self.infer_format_from_handler()
            if self.fastx_format is not None:
                return self.fastx_format
        raise UnknownFastxFormatError(
            f"Could not determine fastX format for {self.source}"
        )

    def __iter__(self) -> Generator[str | bytes, None, None]:
        if self.handler is None:
            raise ValueError(
                "Handler is not set. Did you forget to open the FastxIOWrapper?"
            )
        while self._prefetched_lines:
            yield self._prefetched_lines.pop(0)
        yield from self.handler

    def iter_sequences(self) -> FastxGenerator:
        if self.fastx_format == FastxFormat.FASTQ:
            fastq_fields: FastqGenerator = FastqGeneralIterator(self.__iter__())
            yield from fastq_fields
        elif self.fastx_format == FastxFormat.FASTA:
            fasta_fields: FastaGenerator = SimpleFastaParser(self.__iter__())
            yield from fasta_fields
        else:
            raise UnknownFastxFormatError(
                f"Unknown fastX format for '{self.source}' cannot be iterated over sequences."
            )

    def write_fastq_sequence(self, sequence_fields: FastqFields) -> tuple[int]:
        header, sequence, quality = sequence_fields
        self.handler.write(f"@{header}\n{sequence}\n+\n{quality}\n")
        return len(sequence)

    def write_fastq_sequences(self, sequences: FastqIterator) -> Iterable[int]:
        for sequence_fields in sequences:
            yield self.write_fastq_sequence(sequence_fields)

    def write_fastq(self, sequences: FastqIterator) -> tuple[int, int]:
        bases_written, sequences_written = (0, 0)
        for sequences_written, seq_len in enumerate(
            self.write_fastq_sequences(sequences), start=1
        ):
            bases_written += seq_len
        self.flush()
        return bases_written, sequences_written

    def write_fasta_sequence(
        self,
        sequence_fields: FastaFields,
        max_bases_per_line_for_multifasta: int | None = None,
    ) -> tuple[int]:
        header = sequence_fields[0]
        sequence = sequence_fields[1]
        seq_len = len(sequence)
        if (
            max_bases_per_line_for_multifasta is not None
            and seq_len > max_bases_per_line_for_multifasta > 0
        ):
            sequence = "\n".join(
                sequence[i : i + max_bases_per_line_for_multifasta]
                for i in range(0, seq_len, max_bases_per_line_for_multifasta)
            )
        self.handler.write(f">{header}\n{sequence}\n")
        return seq_len

    def write_fasta(
        self,
        sequences: FastaIterator,
        max_bases_per_line_for_multifasta: int | None = None,
    ) -> tuple[int, int]:
        bases_written, sequences_written = (0, 0)
        if not max_bases_per_line_for_multifasta:
            for sequences_written, (header, sequence) in enumerate(sequences, start=1):
                self.handler.write(f">{header}\n{sequence}\n")
                bases_written += len(sequence)
            return bases_written, sequences_written
        for sequences_written, (header, sequence) in enumerate(sequences, start=1):
            seq_len = len(sequence)
            if seq_len > max_bases_per_line_for_multifasta:
                sequence = "\n".join(
                    sequence[i : i + max_bases_per_line_for_multifasta]
                    for i in range(0, seq_len, max_bases_per_line_for_multifasta)
                )
            bases_written += seq_len
        self.flush()
        return bases_written, sequences_written

    def write_fastx_sequence(
        self,
        sequence_fields: FastxFields,
        max_bases_per_line_for_multifasta: int | None = None,
        quality_simulator: QualitySimulator | None = None,
    ) -> tuple[int]:
        if self.fastx_format == FastxFormat.FASTQ:
            fastq_fields = convert_sequence_to_fastq(
                sequence_fields, quality_simulator=quality_simulator
            )
            return self.write_fastq_sequence(fastq_fields)
        elif self.fastx_format == FastxFormat.FASTA:
            return self.write_fasta_sequence(
                sequence_fields,
                max_bases_per_line_for_multifasta=max_bases_per_line_for_multifasta,
            )
        raise UnknownFastxFormatError(
            f"Unknown fastX format for '{self.source}' - cannot write sequence."
        )

    def write_fastx(
        self,
        sequences: FastxIterator,
        max_bases_per_line_for_multifasta: int | None = None,
        quality_simulator: QualitySimulator | None = None,
    ) -> tuple[int, int]:
        if self.fastx_format == FastxFormat.FASTQ:
            fastq_sequences: FastqGenerator = convert_sequences_to_fastq(
                sequences, quality_simulator=quality_simulator
            )
            result = self.write_fastq(sequences=fastq_sequences)
        elif self.fastx_format == FastxFormat.FASTA:
            fasta_sequences: FastaGenerator = (
                (fields[0], fields[1]) for fields in sequences
            )
            result = self.write_fasta(
                sequences=fasta_sequences,
                max_bases_per_line_for_multifasta=max_bases_per_line_for_multifasta,
            )
        else:
            raise UnknownFastxFormatError(
                f"Unknown fastX format for '{self.source}' cannot be written."
            )
        self.flush()
        return result


class PairedFastxIOWrapper:
    def __init__(self, read1_io: FastxIOWrapper, read2_io: FastxIOWrapper):
        self.read1_io = read1_io
        self.read2_io = read2_io
        assert self.read1_io.fastx_format == self.read2_io.fastx_format, (
            f"PairedFastxIOWrapper requires both read1 and read2 to have the same fastx_format, got {self.read1_io.fastx_format} and {self.read2_io.fastx_format}"
        )
        assert self.read1_io.mode == self.read2_io.mode, (
            f"PairedFastxIOWrapper requires both read1 and read2 to have the same mode, got {self.read1_io.mode} and {self.read2_io.mode}"
        )
        self._toggle_write = True
        self._toggle_read = True
        self._toggle_general = True

    def __repr__(self) -> str:
        return (
            f"PairedFastxIOWrapper(read1={self.read1_io}, read2={self.read2_io})"
            f"fastx_format={self.fastx_format}, mode='{self.mode}')"
        )

    @property
    def fastx_format(self) -> FastxFormat:
        if self.read1_io.fastx_format != self.read2_io.fastx_format:
            raise ValueError(
                f"PairedFastxIOWrapper requires both read1 and read2 to have the same fastx_format, got {self.read1_io.fastx_format} and {self.read2_io.fastx_format}"
            )
        return self.read1_io.fastx_format

    @property
    def mode(self) -> str:
        if self.read1_io.mode != self.read2_io.mode:
            raise ValueError(
                f"PairedFastxIOWrapper requires both read1 and read2 to have the same mode, got {self.read1_io.mode} and {self.read2_io.mode}"
            )
        return self.read1_io.mode

    def get_file_io(self, index: int) -> FastxIOWrapper:
        if index == 1:
            return self.read1_io
        elif index == 2:
            return self.read2_io
        else:
            raise ValueError(f"Index must be 1 or 2, got {index}")

    def write(self, s: str) -> int:
        target = self.read1_io if self._toggle_write else self.read2_io
        self._toggle_write = not self._toggle_write
        return target.write(s)

    def read(self) -> str | bytes:
        target = self.read1_io if self._toggle_read else self.read2_io
        self._toggle_read = not self._toggle_read
        return target.read()

    def getattr_alternate(self, attr_name: str, default: Any = None) -> Any:
        target = self.read1_io if self._toggle_general else self.read2_io
        self._toggle_general = not self._toggle_general
        return getattr(target, attr_name, default)

    def call_alternate(self, func_name: str, *args, **kwargs) -> Any:
        return self.getattr_alternate(func_name)(*args, **kwargs)

    def getattr_both(self, attr_name: str) -> tuple[Any, Any]:
        exceptions = list()
        attrs = list()
        for file_io in (self.read1_io, self.read2_io):
            try:
                attrs.append(getattr(file_io, attr_name))
            except Exception as exc:
                exceptions.append(exc)
        if not exceptions:
            return tuple(attrs)
        raise ExceptionGroup(
            f"Multiple exceptions in PairedFastxIOWrapper.getattr_both({attr_name})",
            exceptions,
        )

    def call_both(self, func_name: str, *args, **kwargs) -> tuple[Any, Any]:
        exceptions = list()
        results = list()
        for func in self.getattr_both(func_name):
            try:
                results.append(func(*args, **kwargs))
            except Exception as exc:
                exceptions.append(exc)
        if not exceptions:
            return tuple(results)
        args_str = ", ".join(map(str, args))
        kwargs_str = ", ".join(f"{k}={v}" for k, v in kwargs.items())
        full_args_str = ", ".join(filter(None, [args_str, kwargs_str]))
        raise ExceptionGroup(
            f"Multiple exceptions in PairedFastxIOWrapper.call_both.{func_name}({full_args_str})",
            exceptions,
        )

    def call_both_as_generators(
        self, func_name: str, *args, **kwargs
    ) -> Generator[tuple[Any, Any], None, None]:
        r1_gen, r2_gen = self.call_both(func_name, *args, **kwargs)
        sentinel_value = object()
        if isinstance(r1_gen, Generator) and isinstance(r2_gen, Generator):
            yield from self.iter_on_iterable_of_paired_values(
                zip_longest(r1_gen, r2_gen, fillvalue=sentinel_value),
                function_name=f"call_both_as_generators.{func_name}",
                sentinel_value=sentinel_value,
            )
            return
        raise TypeError(
            f"PairedFastxIOWrapper.call_both.{func_name} expected generator return values from both, but got {type(r1_gen)} and {type(r2_gen)}."
        )

    def call_both_with_different_parameters_return_tuple(
        self,
        func_name: str,
        r1_args: tuple | None = None,
        r2_args: tuple | None = None,
        r1_kwargs: dict[str, Any] | None = None,
        r2_kwargs: dict[str, Any] | None = None,
    ) -> tuple[Any, Any]:
        exceptions = list()
        results = list()
        r1_func, r2_func = self.getattr_both(func_name)
        try:
            results.append(r1_func(*(r1_args or ()), **(r1_kwargs or {})))
        except Exception as exc:
            exceptions.append(exc)
        try:
            results.append(r2_func(*(r2_args or ()), **(r2_kwargs or {})))
        except Exception as exc:
            exceptions.append(exc)
        if not exceptions:
            return tuple(results)
        r1_args_str = ", ".join(map(str, r1_args or ()))
        r1_kwargs_str = ", ".join(f"{k}={v}" for k, v in (r1_kwargs or {}).items())
        r1_full_args_str = ", ".join(filter(None, [r1_args_str, r1_kwargs_str]))
        r2_args_str = ", ".join(map(str, r2_args or ()))
        r2_kwargs_str = ", ".join(f"{k}={v}" for k, v in (r2_kwargs or {}).items())
        r2_full_args_str = ", ".join(filter(None, [r2_args_str, r2_kwargs_str]))
        raise ExceptionGroup(
            f"Multiple exceptions in PairedFastxIOWrapper.call_both.{func_name}(r1=[{r1_full_args_str}], r2=[{r2_full_args_str}])",
            exceptions,
        )

    def call_both_with_different_parameters_return_generator(
        self,
        func_name: str,
        r1_args: tuple | None = None,
        r2_args: tuple | None = None,
        r1_kwargs: dict[str, Any] | None = None,
        r2_kwargs: dict[str, Any] | None = None,
    ) -> tuple[Any, Any] | Generator[tuple[Any, Any], None, None]:
        r1_gen, r2_gen = self.call_both_with_different_parameters_return_tuple(
            func_name=func_name,
            r1_args=r1_args,
            r2_args=r2_args,
            r1_kwargs=r1_kwargs,
            r2_kwargs=r2_kwargs,
        )
        if not (isinstance(r1_gen, Generator) and isinstance(r2_gen, Generator)):
            raise TypeError(
                f"PairedFastxIOWrapper.call_both.{func_name} expected generator return values from both, but got {type(r1_gen)} and {type(r2_gen)}."
            )
        sentinel_value = object()
        yield from PairedFastxIOWrapper.iter_on_iterable_of_paired_values(
            zip_longest(r1_gen, r2_gen, fillvalue=sentinel_value),
            function_name=f"call_both.{func_name}",
            sentinel_value=sentinel_value,
        )

    def call_both_with_different_parameters(
        self,
        func_name: str,
        r1_args: tuple | None = None,
        r2_args: tuple | None = None,
        r1_kwargs: dict[str, Any] | None = None,
        r2_kwargs: dict[str, Any] | None = None,
    ) -> tuple[Any, Any] | Generator[tuple[Any, Any], None, None]:
        values = self.call_both_with_different_parameters_return_tuple(
            func_name=func_name,
            r1_args=r1_args,
            r2_args=r2_args,
            r1_kwargs=r1_kwargs,
            r2_kwargs=r2_kwargs,
        )
        if (isinstance(value, Generator) for value in values):
            for value1, value2 in zip(*values):
                yield (value1, value2)
            return
        return values

    def flush(self) -> None:
        self.call_both("flush")

    def close(self) -> None:
        self.call_both("close")

    def open(self) -> "PairedFastxIOWrapper":
        self.read1_io, self.read2_io = self.call_both("open")
        return self

    def __enter__(self):
        self.call_both("__enter__")
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> None:
        try:
            self.call_both("__exit__", exc_type, exc_val, exc_tb)
        except Exception as exc_group:
            if exc_val is not None:
                raise exc_group from exc_val
            raise exc_group

    @property
    def readable(self) -> bool:
        return all(self.call_both("readable"))

    @property
    def writable(self) -> bool:
        return all(self.call_both("writable"))

    @property
    def closed(self) -> bool:
        return any(self.call_both("closed"))

    @staticmethod
    def iter_on_iterable_of_paired_values(
        iterable: Iterable[tuple[Any, Any]],
        function_name: str = "",
        sentinel_value: Any = None,
    ) -> Generator[tuple[Any, Any], None, None]:
        for i, (r1_value, r2_value) in enumerate(iterable, start=1):
            if r1_value is sentinel_value:
                raise IterableMismatchedCounts(
                    function_name=function_name, missing_data="read1", index=i
                )
            elif r2_value is sentinel_value:
                raise IterableMismatchedCounts(
                    function_name=function_name, missing_data="read2", index=i
                )
            yield (r1_value, r2_value)

    def __iter__(
        self,
    ) -> Generator[tuple[str | bytes, str | bytes], None, None]:
        yield from self.call_both_as_generators("__iter__")

    def iter_sequences(
        self,
    ) -> PairedFastxGenerator:
        if self.fastx_format not in FastxFormat:
            if self.read1_io.fastx_format != self.read2_io.fastx_format:
                raise ValueError(
                    "PairedFastxIOWrapper requires both read1 and read2 to have the same fastx_format,"
                    f" got {self.read1_io.fastx_format} != {self.read2_io.fastx_format}"
                    f" for '{self.read1_io.source}' and '{self.read2_io.source}'"
                )
            self.fastx_format = self.read1_io.fastx_format
        if self.fastx_format == FastxFormat.FASTQ:
            paired_fastq_generator: PairedFastqGenerator = self.call_both_as_generators(
                "iter_sequences"
            )
            yield from paired_fastq_generator
        elif self.fastx_format == FastxFormat.FASTA:
            paired_fasta_generator: PairedFastaGenerator = self.call_both_as_generators(
                "iter_sequences"
            )
            yield from paired_fasta_generator
        else:
            raise UnknownFastxFormatError(
                f"Unknown fastX format for '{self.read1_io.source}' and '{self.read2_io.source}' cannot be iterated over sequences."
            )

    def write_fastx_to_fasta(
        self,
        sequences: PairedFastxIterator,
        max_bases_per_line_for_multifasta: int | None = None,
    ):
        bases_written_r1 = 0
        bases_written_r2 = 0
        for sequences_written, (fields1, fields2) in enumerate(sequences, start=1):
            if fields1 is None and fields2 is None:
                break
            elif fields1 is None:
                raise IterableMismatchedCounts(
                    function_name="PairedFastxIOWrapper.write_fastx",
                    missing_data="read1",
                    index=sequences_written,
                )
            elif fields2 is None:
                raise IterableMismatchedCounts(
                    function_name="PairedFastxIOWrapper.write_fastx",
                    missing_data="read2",
                    index=sequences_written,
                )
            r1_seq_len = self.read1_io.write_fasta_sequence(
                fields1,
                max_bases_per_line_for_multifasta=max_bases_per_line_for_multifasta,
            )
            r2_seq_len = self.read2_io.write_fasta_sequence(
                fields2,
                max_bases_per_line_for_multifasta=max_bases_per_line_for_multifasta,
            )
            bases_written_r1 += r1_seq_len
            bases_written_r2 += r2_seq_len
            sequences_written += 1
        self.flush()
        return (bases_written_r1 + bases_written_r2), (sequences_written * 2)

    def write_fastx_to_fastq(
        self,
        sequences: PairedFastxIterator,
        quality_simulator: QualitySimulator | None = None,
    ):
        bases_written_r1 = 0
        bases_written_r2 = 0

        # convert_sequence_pairs_to_fastq yields pairs of FastqFields tuples:
        r1_fastq_sequences, r2_fastq_sequences = convert_sequence_pairs_to_fastq(
            sequences, quality_simulator=quality_simulator
        )
        # Write sequences and get bases written for each read
        for r1_bases, r2_bases in self.call_both_with_different_parameters(
            "write_fastq_sequences",
            r1_args=(r1_fastq_sequences,),
            r2_args=(r2_fastq_sequences,),
        ):
            bases_written_r1 += r1_bases
            bases_written_r2 += r2_bases
        self.flush()
        return bases_written_r1, bases_written_r2

    def write_fastx(
        self,
        sequences: PairedFastxIterator,
        max_bases_per_line_for_multifasta: int | None = None,
        quality_simulator: QualitySimulator | None = None,
    ) -> tuple[int, int]:
        if self.read1_io.fastx_format != self.read2_io.fastx_format:
            raise ValueError(
                "PairedFastxIOWrapper requires both read1 and read2 to have the same fastx_format,"
                f" got {self.read1_io.fastx_format} != {self.read2_io.fastx_format}"
                f" for '{self.read1_io.source}' and '{self.read2_io.source}'"
            )
        if self.fastx_format not in FastxFormat:
            self.fastx_format = self.read1_io.fastx_format
        if self.fastx_format == FastxFormat.FASTQ:
            result = self.write_fastx_to_fastq(
                sequences=sequences,
                quality_simulator=quality_simulator,
            )
            self.flush()
            return result
        if self.fastx_format == FastxFormat.FASTA:
            result = self.write_fastx_to_fasta(
                sequences=sequences,
                max_bases_per_line_for_multifasta=max_bases_per_line_for_multifasta,
            )
            return result
        raise UnknownFastxFormatError(
            f"Unknown fastX format for '{self.read1_io.source}' and '{self.read2_io.source}' cannot be written."
        )
