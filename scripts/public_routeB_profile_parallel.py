#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

def one(line: str, db: Path, threads: int, scratch: Path) -> tuple[str, str, str]:
    sample, fastq, outdir = line.split('\t')
    fastq, outdir = Path(fastq), Path(outdir)
    final = outdir / 'noGC_glm' / 'output.tsv'
    if final.exists() and final.stat().st_size > 0:
        return sample, 'skipped', str(final)
    task = scratch / sample
    if task.exists():
        shutil.rmtree(task)
    task.mkdir(parents=True)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd_count = [
        os.environ.get('SK2BGROW_BIN', 'sk2bgrow'), 'profile', str(fastq),
        '--db', str(db), '--output', str(task), '--mode', '2brad',
        '--max-mismatch', '0', '--threads', str(threads), '--no-stats', '--quiet'
    ]
    with (outdir / 'count.log').open('w') as h:
        subprocess.run(cmd_count, stdout=h, stderr=subprocess.STDOUT, check=True)
    counts = list(task.glob('*.counts.tsv'))
    if not counts:
        raise RuntimeError('no counts produced')
    envpy = os.environ.get('SK2BGROW_PY')
    cmd_stats = [
        envpy, '-m', 'sk2bgrow.cli', 'profile', str(counts[0]), '--db', str(db),
        '--output', str(outdir / 'noGC_glm'), '--no-gc-correct', '--method', 'glm'
    ]
    with (outdir / 'stats.log').open('w') as h:
        subprocess.run(cmd_stats, stdout=h, stderr=subprocess.STDOUT, check=True)
    if not (outdir / 'noGC_glm' / 'output.tsv').exists():
        raise RuntimeError('stats did not create output.tsv')
    shutil.rmtree(task, ignore_errors=True)
    return sample, 'complete', str(final)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True, type=Path)
    ap.add_argument('--db', required=True, type=Path)
    ap.add_argument('--scratch', required=True, type=Path)
    ap.add_argument('--threads-per-sample', type=int, default=2)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--status', required=True, type=Path)
    a = ap.parse_args()
    lines = [x.rstrip('\n') for x in a.manifest.open() if x.strip()]
    a.scratch.mkdir(parents=True, exist_ok=True)
    a.status.parent.mkdir(parents=True, exist_ok=True)
    failed = []
    with a.status.open('w') as status, ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one, line, a.db, a.threads_per_sample, a.scratch): line for line in lines}
        for n, fut in enumerate(as_completed(futs), 1):
            line = futs[fut]
            sample = line.split('\t')[0]
            try:
                sample, state, final = fut.result()
                status.write(f'{sample}\t{state}\t{final}\n')
                print(f'[{n}/{len(lines)}] {sample} {state}', flush=True)
            except Exception as e:
                failed.append(sample)
                status.write(f'{sample}\tfailed\t{e}\n')
                print(f'[{n}/{len(lines)}] {sample} FAILED: {e}', file=sys.stderr, flush=True)
    print(f'FAILED_COUNT={len(failed)}', file=sys.stderr)
    return 1 if failed else 0

if __name__ == '__main__':
    raise SystemExit(main())
