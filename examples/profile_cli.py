"""Measure one fresh CLI invocation, including its descendant processes.

Used by the example notebooks; install vdjmatch[bench,notebook]. RSS is sampled
at 50 ms and summed across live processes (shared pages may be counted twice).
CPU peaks use one-second windows; shorter windows amplify delayed OS accounting.
"""
from pathlib import Path
import json
import subprocess
import time

import psutil


def profile(command, log, *, timeout=300, max_rss_bytes=None):
    """Measure a fresh command; optional RSS ceiling fails and kills its process tree."""
    start = time.perf_counter()
    cpu_by_pid = {}
    peak_rss = peak_threads = 0
    peak_cores = None
    previous_cpu = 0.0
    previous_time = start
    completed = False
    failure_kind = None
    processes = []
    with Path(log).open('w') as stream:
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
                rss = threads = 0
                for process in processes:
                    try:
                        with process.oneshot():
                            rss += process.memory_info().rss
                            threads += process.num_threads()
                            usage = process.cpu_times()
                            cpu_by_pid[process.pid] = usage.user + usage.system
                    except psutil.NoSuchProcess:
                        continue
                cpu = sum(cpu_by_pid.values())
                observed = time.perf_counter()
                peak_rss = max(peak_rss, rss)
                peak_threads = max(peak_threads, threads)
                if max_rss_bytes is not None and rss > max_rss_bytes:
                    raise MemoryError(f'process-tree RSS {rss} exceeds budget {max_rss_bytes} bytes')
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
                      'cpu_seconds': sum(cpu_by_pid.values()),
                      'mean_cpu_cores': sum(cpu_by_pid.values()) / elapsed,
                      'peak_sampled_cpu_cores': peak_cores,
                      'peak_live_threads': peak_threads, 'sample_interval_seconds': .05,
                      'cpu_peak_window_seconds': 1.0,
                      'completed': completed, 'exit_code': child.returncode,
                      'failure_kind': failure_kind}
            Path(log).with_suffix('.resources.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
