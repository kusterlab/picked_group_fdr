import csv
from typing import List
from contextlib import contextmanager

# csv.field_size_limit(sys.maxsize)
csv.field_size_limit(2147483647)


@contextmanager
def get_tsv_reader(filename: str, delimiter: str = "\t"):
    with open(filename, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter=delimiter)
        yield reader


@contextmanager
def get_tsv_writer(filename: str, delimiter: str = "\t"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        yield writer


def get_delimiter(filename: str):
    if filename.endswith(".csv"):
        return ","
    else:
        return "\t"


def get_column_index(
    headers: List[str], column_name: str, is_optional: bool = False
) -> int:
    if column_name not in headers:
        if is_optional:
            return -1
        else:
            raise ValueError(
                f"Column {column_name} is missing. Please check your input file."
            )
    return headers.index(column_name)


def get_header_col_func(headers: List[str]):
    def get_header_col(name, required=False):
        if required:
            return get_column_index(headers, name)
        else:
            if name in headers:
                return get_column_index(headers, name)
            return -1

    return get_header_col


def get_header_cols_starting_with_func(headers: List[str]):
    def get_header_cols_starting_with(name):
        cols = [idx for idx, h in enumerate(headers) if h.startswith(name)]
        return cols

    return get_header_cols_starting_with
