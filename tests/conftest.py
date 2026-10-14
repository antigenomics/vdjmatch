"""Keep the in-process unit-test table pool within the smallest API budget."""
import os

os.environ.setdefault("POLARS_MAX_THREADS", "1")
