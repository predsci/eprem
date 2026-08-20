#!/usr/bin/env python3
"""
epremdigest.py - Build a small, text-only digest of EPREM netCDF output.

Instead of storing large binary reference netCDF files, this reads each
obs*.nc in a run directory and writes a compact text digest containing,
per variable:
  - name, shape, dtype
  - global stats: min, max, mean, sum, rms  (rounded to 6 sig figs)
  - SHA-256 of the (little-endian normalized) raw variable bytes

Optionally, per-slice aggregates (max, sum) along a chosen axis so a
localized error can't be hidden in the global mean.

Usage:
  eprem_digest.py <rundir> [-o digest.txt] [--granularity {global,time,node}]

  <rundir>           directory containing obs*.nc files
  -o, --output       output digest file (default: stdout)
  --granularity      global (default) | time (per first axis) | node (per 2nd axis)
  --tol              not used for digest creation; informational only
"""
import argparse
import glob
import hashlib
import os
import sys

import numpy as np
import netCDF4 as nc


def sigfig(x, n=6):
    """Round a float to n significant figures, returning a stable string."""
    if x is None:
        return "nan"
    x = float(x)
    if np.isnan(x):
        return "nan"
    if np.isinf(x):
        return "inf" if x > 0 else "-inf"
    if x == 0.0:
        return "0.0"
    from math import floor, log10
    d = n - int(floor(log10(abs(x)))) - 1
    val = round(x, d)
    return "%.6g" % val


def var_hash(arr):
    """SHA-256 of the variable's raw bytes, normalized to little-endian."""
    a = np.asarray(arr)
    if a.dtype.kind == "f":
        a = a.astype("<" + a.dtype.str[1:])  # normalize endianness
    return hashlib.sha256(a.tobytes()).hexdigest()


def global_stats(arr):
    a = np.asarray(arr, dtype="float64")
    if a.size == 0:
        return "min=nan max=nan mean=nan sum=0.0 rms=nan"
    with np.errstate(invalid="ignore", divide="ignore"):
        mn = float(np.nanmin(a))
        mx = float(np.nanmax(a))
        mean = float(np.nanmean(a))
        s = float(np.nansum(a))
        rms = float(np.sqrt(np.nanmean(a * a)))
    return "min=%s max=%s mean=%s sum=%s rms=%s" % (
        sigfig(mn), sigfig(mx), sigfig(mean), sigfig(s), sigfig(rms))


def slice_stats(arr, axis):
    """Per-slice max and sum along the given axis (one line per slice)."""
    a = np.asarray(arr, dtype="float64")
    if a.ndim <= axis:
        return []
    lines = []
    with np.errstate(invalid="ignore"):
        for i in range(a.shape[axis]):
            sl = np.take(a, i, axis=axis)
            lines.append("  slice[%d]=%d max=%s sum=%s" % (
                i, axis, sigfig(float(np.nanmax(sl))), sigfig(float(np.nansum(sl)))))
    return lines


def digest_variable(name, var, granularity):
    arr = np.asarray(var)
    out = []
    out.append("var %s" % name)
    out.append("  shape=%s dtype=%s" % (str(arr.shape), str(arr.dtype)))
    out.append("  " + global_stats(arr))
    out.append("  sha256=%s" % var_hash(arr))
    if granularity == "time" and arr.ndim >= 1:
        out.extend(slice_stats(arr, 0))
    elif granularity == "node" and arr.ndim >= 2:
        out.extend(slice_stats(arr, 1))
    return out


def digest_file(path, granularity):
    lines = []
    ds = nc.Dataset(path, "r")
    try:
        lines.append("file %s" % os.path.basename(path))
        for k in sorted(ds.variables.keys()):
            v = ds.variables[k]
            # Skip scalar attributes-only entries; only digest array vars.
            if len(v.shape) == 0:
                continue
            lines.extend(digest_variable(k, v, granularity))
    finally:
        ds.close()
    return lines


def main():
    ap = argparse.ArgumentParser(description="Build a text digest of EPREM netCDF output.")
    ap.add_argument("rundir", help="Directory containing obs*.nc files")
    ap.add_argument("-o", "--output", default=None, help="Output digest file (default stdout)")
    ap.add_argument("--granularity", choices=["global", "time", "node"], default="global")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.rundir, "obs*.nc")))
    if not files:
        sys.stderr.write("WARNING: no obs*.nc files found in %s\n" % args.rundir)

    all_lines = ["# eprem_digest granularity=%s" % args.granularity]
    for f in files:
        all_lines.extend(digest_file(f, args.granularity))
    all_lines.append("")

    text = "\n".join(all_lines)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(text)
        #sys.stderr.write("Wrote digest for %d file(s) to %s\n" % (len(files), args.output))
    else:
        print(text)


if __name__ == "__main__":
    main()