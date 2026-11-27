"""Plan spot workers and matrix tiles within a portable memory budget."""
from __future__ import annotations

import os
from pathlib import Path
import re


def memory_bytes(value):
    if value is None:
        return None
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([KMGT]?)(?:i?B)?", str(value).strip(), re.I)
    if not match or float(match[1]) <= 0:
        raise ValueError("Memory limit must be positive, for example 512M, 4G or 16GiB")
    return int(float(match[1]) * 1024 ** ("", "K", "M", "G", "T").index(match[2].upper()))


def available_memory():
    limits = []
    try:
        fields = Path("/proc/meminfo").read_text().splitlines()
        available = next(line for line in fields if line.startswith("MemAvailable:"))
        limits.append(int(available.split()[1]) * 1024)
    except (OSError, StopIteration, ValueError):
        try:
            limits.append(os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))
        except (ValueError, OSError, AttributeError):
            pass
    for limit_path, usage_path in (
        ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
        ("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes"),
    ):
        try:
            limit = int(Path(limit_path).read_text())
            usage = int(Path(usage_path).read_text())
            limits.append(max(0, limit - usage))
        except (OSError, ValueError):
            pass
    return min(limits, default=8 * 2**30)


def worker_count(value):
    if value not in (None, "auto"):
        return int(value)
    try:
        count = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        count = os.cpu_count() or 1
    return max(1, min(4, count - 1))


def plan_blocks(source, n_genes, n_spots, n_domains, *, threads=None, block_size=None,
                projection_size=None, memory_limit=None):
    """Return conservative working-set estimates, including output-table reserve.

    The budget guides allocation. Python/Arrow allocator overhead, the OS page
    cache, and plotting can cause process RSS to differ from the estimate.
    """
    automatic_threads = threads in (None, "auto")
    workers = worker_count(threads)
    block = min(block_size or 512, n_spots)
    projection = min(projection_size or 4096, n_spots)
    budget = memory_bytes(memory_limit)
    if budget is None:
        budget = min(4 * 2**30, available_memory() // 2)
    fixed = 128 * 2**20 + n_spots * 1024 + n_genes * n_domains * 24

    def io_bytes():
        if source.storage == "columnar":
            return len(source.genes) * projection * 16
        extra = len(source.genes) * 8192 if source.layout == "gene-major" else 0
        return len(source.genes) * block * 8 + extra

    def compute_bytes():
        return (workers + 1) * block * n_genes * 40

    def estimate():
        return fixed + compute_bytes() + io_bytes()

    # Preserve requested values. Split reductions between input projections and
    # compute blocks so small budgets do not force one projection per spot.
    while estimate() > budget:
        shrink_projection = (source.storage == "columnar" and projection_size is None
                             and projection > 1)
        shrink_block = block_size is None and block > 1
        if shrink_projection and (not shrink_block or io_bytes() >= compute_bytes()):
            projection = max(1, projection // 2)
        elif shrink_block:
            block = max(1, block // 2)
        elif automatic_threads and workers > 1:
            workers -= 1
        else:
            raise ValueError(
                f"Estimated analysis working set {estimate() / 2**30:.2f} GiB exceeds "
                f"--memory-limit {budget / 2**30:.2f} GiB. Increase the budget or reduce "
                "--threads, --block-size and --projection-size."
            )
    return {"threads": workers, "spot_block": block,
            "projection_size": projection if source.storage == "columnar" else None,
            "memory_budget_bytes": budget, "estimated_working_set_bytes": estimate(),
            "memory_policy": "estimated analysis working set; plotting and allocator overhead vary"}
