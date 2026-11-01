"""Measure one fresh CLI invocation, including its descendant processes.

Used by the example notebooks; install vdjmatch[bench,notebook]. RSS is sampled
at 50 ms and summed across live processes (shared pages may be counted twice).
CPU peaks use one-second windows; shorter windows amplify delayed OS accounting.
"""
from pathlib import Path
import json
import subprocess
import time
import sys
import ctypes
import errno
from contextlib import nullcontext

import psutil


if sys.platform == 'darwin':
    # sys/resource.h rusage_info_v2: footprint includes compressed/swapped pages.
    class _Rusage(ctypes.Structure):
        _fields_ = [('uuid', ctypes.c_uint8 * 16), ('values', ctypes.c_uint64 * 18)]
    _libproc = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)
    _libproc.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
    _libproc.proc_pid_rusage.restype = ctypes.c_int


def _accounted_memory(pid, rss):
    """Charge at least RSS, including macOS compressed/swapped physical footprint."""
    if sys.platform != 'darwin':
        return rss
    usage = _Rusage()
    if _libproc.proc_pid_rusage(pid, 2, ctypes.byref(usage)):
        error = ctypes.get_errno()
        if error == errno.ESRCH:
            raise psutil.NoSuchProcess(pid)
        raise OSError(error, 'proc_pid_rusage failed')
    return max(rss, usage.values[7])


def profile(command, log, *, timeout=300, max_rss_bytes=None, trace_path=None):
    """Measure a fresh command; optional RSS ceiling fails and kills its process tree."""
    start = time.perf_counter()
    cpu_by_pid = {}
    peak_rss = peak_accounted = peak_threads = 0
    peak_cores = None
    previous_cpu = 0.0
    previous_time = start
    completed = False
    failure_kind = None
    processes = []
    with Path(log).open('w') as stream, (Path(trace_path).open('w') if trace_path else nullcontext()) as trace:
        if trace:
            trace.write('seconds\ttree_rss_bytes\taccounted_bytes\tcpu_seconds\tlive_threads\n')
        child = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT)
        try:
            parent = psutil.Process(child.pid)
            processes = [parent]
            while child.poll() is None:
                now = time.perf_counter()
                if now - start > timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
                try:
                    processes = [parent, *parent.children(recursive=True)]
                except psutil.NoSuchProcess:
                    break
                rss = accounted = threads = 0
                for process in processes:
                    try:
                        with process.oneshot():
                            resident = process.memory_info().rss
                            rss += resident
                            accounted += _accounted_memory(process.pid, resident)
                            threads += process.num_threads()
                            usage = process.cpu_times()
                            cpu_by_pid[process.pid] = usage.user + usage.system
                    except psutil.NoSuchProcess:
                        continue
                cpu = sum(cpu_by_pid.values())
                observed = time.perf_counter()
                if trace:
                    trace.write(f'{observed-start:.6f}\t{rss}\t{accounted}\t{cpu:.6f}\t{threads}\n')
                    trace.flush()
                peak_rss = max(peak_rss, rss)
                peak_accounted = max(peak_accounted, accounted)
                peak_threads = max(peak_threads, threads)
                if max_rss_bytes is not None and accounted > max_rss_bytes:
                    raise MemoryError(f'process-tree accounted memory {accounted} (RSS {rss}) exceeds budget {max_rss_bytes} bytes')
                if previous_time == start:
                    previous_cpu, previous_time = cpu, observed
                elif observed - previous_time >= 1.0:
                    peak_cores = max(peak_cores or 0., (cpu - previous_cpu)/(observed - previous_time))
                    previous_cpu, previous_time = cpu, observed
                time.sleep(.05)
            if child.wait():
                raise subprocess.CalledProcessError(child.returncode, command)
            completed = True
        except BaseException as error:
            failure_kind = type(error).__name__
            raise
        finally:
            if child.poll() is None:
                for process in reversed(processes):
                    try:
                        process.kill()
                    except psutil.NoSuchProcess:
                        pass
                child.wait()
            elapsed = time.perf_counter() - start
            result = {'seconds': elapsed, 'peak_tree_rss_GB': peak_rss / 1e9,
                      'peak_tree_accounted_GB': peak_accounted / 1e9,
                      'memory_limit_basis': 'sum max(RSS, macOS physical footprint) per process',
                      'cpu_seconds': sum(cpu_by_pid.values()),
                      'mean_cpu_cores': sum(cpu_by_pid.values()) / elapsed,
                      'peak_sampled_cpu_cores': peak_cores,
                      'peak_live_threads': peak_threads, 'sample_interval_seconds': .05,
                      'cpu_peak_window_seconds': 1.0,
                      'completed': completed, 'exit_code': child.returncode,
                      'failure_kind': failure_kind}
            Path(log).with_suffix('.resources.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
