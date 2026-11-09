"""Kept for the prototype scripts: the reader now lives in tools/pkcreader/pkcreader.py."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pkcreader"))
from pkcreader import HEADER_SIZE, PKCFile, decode_block, parse_header, read_snapshot  # noqa: E402,F401
