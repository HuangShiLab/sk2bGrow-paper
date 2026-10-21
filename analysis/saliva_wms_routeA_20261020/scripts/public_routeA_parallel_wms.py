#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

def run(cmd, log):
    with log.open('w') as h:
        subprocess.run(cmd, stdout=h, stderr=subprocess.STDOUT, check=True)

def one(line, db, threads, scratch):
    sample, fastq_field, outdir_name = line.split('\t')
    fastqs=[Path(x) for x in fastq_field.split(';') if x]
    outdir=Path(outdir_name)
    finals=[outdir/method/'output.tsv' for method in ('noGC_auto','noGC_glm')]
    if all(p.exists() and p.stat().st_size>0 for p in finals):
        return sample,'skipped',str(finals[0])
    missing=[str(p) for p in fastqs if not p.exists() or p.stat().st_size==0]
    if missing: raise RuntimeError('missing FASTQ(s): '+','.join(missing))
    task=scratch/sample
    if task.exists(): shutil.rmtree(task)
    task.mkdir(parents=True)
    outdir.mkdir(parents=True,exist_ok=True)
    binpath=os.environ.get('SK2BGROW_BIN','sk2bgrow')
    count_cmd=[binpath,'profile',*[str(p) for p in fastqs], '--db',str(db),
               '--output',str(task),'--mode','wms','--max-mismatch','2',
               '--windowing','anchors','--window-anchors','100',
               '--threads',str(threads),'--no-stats','--quiet']
    run(count_cmd,outdir/'count.log')
    counts=sorted(task.glob('*.counts.tsv'))
    if len(counts)!=1: raise RuntimeError(f'expected one count table, got {[x.name for x in counts]}')
    py=os.environ.get('SK2BGROW_PY')
    for method in ('auto','glm'):
        run([py,'-m','sk2bgrow.cli','profile',str(counts[0]),'--db',str(db),
             '--output',str(outdir/('noGC_'+method)),'--no-gc-correct','--method',method],
            outdir/f'stats_{method}.log')
    for p in finals:
        if not p.exists() or p.stat().st_size==0: raise RuntimeError(f'missing {p}')
    shutil.rmtree(task,ignore_errors=True)
    return sample,'complete',str(finals[0])

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest',required=True,type=Path)
    ap.add_argument('--db',required=True,type=Path)
    ap.add_argument('--scratch',required=True,type=Path)
    ap.add_argument('--threads-per-sample',type=int,default=16)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--status',required=True,type=Path)
    a=ap.parse_args()
    lines=[x.rstrip('\n') for x in a.manifest.open() if x.strip()]
    if lines and lines[0].split('\t')[0]=='sample_id': lines=lines[1:]
    a.scratch.mkdir(parents=True,exist_ok=True); a.status.parent.mkdir(parents=True,exist_ok=True)
    failed=[]
    with a.status.open('w') as status, ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs={ex.submit(one,line,a.db,a.threads_per_sample,a.scratch):line for line in lines}
        for n,fut in enumerate(as_completed(futs),1):
            line=futs[fut]; sample=line.split('\t')[0]
            try:
                sample,state,final=fut.result(); status.write(f'{sample}\t{state}\t{final}\n'); status.flush()
                print(f'[{n}/{len(lines)}] {sample} {state}',flush=True)
            except Exception as e:
                failed.append(sample); status.write(f'{sample}\tfailed\t{e}\n'); status.flush()
                print(f'[{n}/{len(lines)}] {sample} FAILED: {e}',file=sys.stderr,flush=True)
    print(f'FAILED_COUNT={len(failed)}',file=sys.stderr)
    return 1 if failed else 0
if __name__=='__main__': raise SystemExit(main())
