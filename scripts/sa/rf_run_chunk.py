#!/usr/bin/env python3
"""Run a range of RF parameter samples (sample<N>.json) through testRun.

Each sample runs once with --seed N, writes into <out>/sample_<N>/, and its
Output_Biomarkers.csv is copied to <home>/output_<N>.csv. Finished samples
(complete CSV already in <home>) are skipped, so a chunk can be resubmitted.

Usage (from the VFB-ABM repo root):
  python3 scripts/sa/rf_run_chunk.py --samples sa/samples --start 1 --end 16 --jobs 16
"""
import argparse, os, shutil, subprocess
from concurrent.futures import ProcessPoolExecutor

def complete(path, numticks):
    return os.path.exists(path) and sum(1 for _ in open(path)) >= numticks + 2

def run_one(args):
    a, n = args
    home_csv = os.path.join(a.home, f"output_{n}.csv")
    if complete(home_csv, a.numticks):
        return n, "skipped"
    cfg = os.path.join(a.samples, f"sample{n}.json")
    if not os.path.exists(cfg):
        return n, "missing config"
    out = os.path.join(a.out, f"sample_{n}")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "stdout.log"), "w") as log:
        rc = subprocess.run([a.bin, "--config", cfg, "--seed", str(n),
                             "--numticks", str(a.numticks),
                             "--wxw", a.world, "--wyw", a.world, "--wzw", a.world,
                             "--output-dir", out],
                            stdout=log, stderr=subprocess.STDOUT).returncode
    csv_path = os.path.join(out, "Output_Biomarkers.csv")
    if rc != 0 or not complete(csv_path, a.numticks):
        return n, f"failed (exit {rc})"
    shutil.copyfile(csv_path, home_csv)
    return n, "ok"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--samples", default="sa/samples")
    p.add_argument("--out", default="sa/runs")
    p.add_argument("--home", default="sa/output_home")
    p.add_argument("--bin", default="build/bin/testRun")
    p.add_argument("--world", default="0.6")
    p.add_argument("--numticks", type=int, default=432)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)   # inclusive
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    a = p.parse_args()
    os.makedirs(a.home, exist_ok=True)
    os.environ["LD_LIBRARY_PATH"] = "build/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
    with ProcessPoolExecutor(a.jobs) as ex:
        for n, status in ex.map(run_one, [(a, n) for n in range(a.start, a.end + 1)]):
            print(f"sample {n}: {status}", flush=True)

if __name__ == "__main__":
    main()
