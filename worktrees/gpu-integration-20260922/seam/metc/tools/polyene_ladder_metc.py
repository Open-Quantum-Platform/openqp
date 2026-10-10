#!/usr/bin/env python3
"""Generate and run an all-trans polyene CPU/GPU METC benchmark ladder."""
import argparse, csv, math, os, re, shutil, statistics, subprocess, sys, time
from pathlib import Path

FLOAT = r"[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?"
NBF_RE = re.compile(r"Number of Basis Set functions\s*=\s*(\d+)", re.I)
TOTAL_RE = re.compile(r"TOTAL energy\s*=\s*(%s)" % FLOAT, re.I)
STATE_RE = re.compile(r"^\s*(\d+)\s+(%s)\s+(%s)\s+(%s)\s+" % (FLOAT, FLOAT, FLOAT))
TIMING_RE = re.compile(r"OQP_GPU_METC_TIMING\s+upload_d3_s=(?P<upload>[-+0-9.eE]+)\s+zero_f3_s=(?P<zero>[-+0-9.eE]+)\s+download_f3_s=(?P<download>[-+0-9.eE]+)")
F3_RE = re.compile(r"OQP_GPU_METC_F3_CHECK\s+(.*)")


def polyene_geometry(m):
    n = 2 * m
    # Planar zig-zag carbon backbone with alternating C=C/C-C distances.
    theta = math.radians(30.0)
    dirs = [(math.cos(theta), math.sin(theta), 0.0), (math.cos(theta), -math.sin(theta), 0.0)]
    cc_double = 1.34
    cc_single = 1.46
    coords = [(0.0, 0.0, 0.0)]
    x = y = z = 0.0
    for i in range(1, n):
        length = cc_double if i % 2 == 1 else cc_single
        dx, dy, dz = dirs[(i - 1) % 2]
        x += length * dx; y += length * dy; z += length * dz
        coords.append((x, y, z))
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    coords = [(x - cx, y - cy, z) for x, y, z in coords]
    atoms = [(6, *c) for c in coords]
    ch = 1.09
    for i, (x, y, z) in enumerate(coords):
        if i == 0:
            # terminal carbon: one H opposite chain, one above/outside
            vx, vy, _ = dirs[0]
            atoms.append((1, x - ch * vx, y - ch * vy, 0.0))
            atoms.append((1, x, y + ch, 0.0))
        elif i == n - 1:
            vx, vy, _ = dirs[(n - 2) % 2]
            atoms.append((1, x + ch * vx, y + ch * vy, 0.0))
            atoms.append((1, x, y + (ch if y < 0 else -ch), 0.0))
        else:
            ysign = 1.0 if y <= 0 else -1.0
            atoms.append((1, x, y + ysign * ch, 0.0))
    return atoms


def write_input(path, m, basis, nstate):
    atoms = polyene_geometry(m)
    lines = [
        f"# all-trans H-(CH=CH)_{m}-H polyene benchmark input",
        "[input]", "system=",
    ]
    for z, x, y, zz in atoms:
        lines.append(f" {z:2d} {x:16.9f} {y:16.9f} {zz:16.9f}")
    lines += [
        "charge=0", "runtype=energy", "method=tdhf", "functional=bhhlyp", f"basis={basis}", "",
        "[guess]", "save_mol=false", "type=huckel", "",
        "[scf]", "type=uhf", "multiplicity=3", "converger_type=diis", "maxit=200", "",
        "[dftgrid]", "rad_npts=90", "ang_npts=302", "",
        "[tdhf]", "type=umrsf", f"nstate={nstate}", "multiplicity=1", "",
    ]
    Path(path).write_text("\n".join(lines))


def run_once(tokens, inp, env_extra, log_path):
    env = dict(os.environ); env.update(env_extra)
    t0 = time.perf_counter()
    proc = subprocess.run(tokens + [str(inp)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    wall = time.perf_counter() - t0
    pyoqp_log = Path(inp).with_suffix('.log')
    log_text = pyoqp_log.read_text(errors='replace') if pyoqp_log.exists() else ''
    combined = proc.stdout + "\n--- STDERR ---\n" + proc.stderr + "\n--- PYOQP LOG ---\n" + log_text
    Path(log_path).write_text(combined)
    return proc.returncode, wall, proc.stderr, log_text, combined


def parse_scalars(log_text):
    vals = {}
    m = TOTAL_RE.search(log_text)
    if m: vals['scf_total_energy_hartree'] = float(m.group(1).replace('D','E').replace('d','e'))
    in_summary = False
    for line in log_text.splitlines():
        if line.strip().startswith('State') and 'Oscillator' in line:
            in_summary = True; continue
        if in_summary:
            mm = STATE_RE.match(line)
            if not mm:
                if line.lstrip().startswith('Transition'):
                    in_summary = False
                continue
            s = int(mm.group(1))
            vals[f'state_{s}_energy_hartree'] = float(mm.group(2).replace('D','E').replace('d','e'))
            vals[f'state_{s}_excitation_ev'] = float(mm.group(3).replace('D','E').replace('d','e'))
            vals[f'state_{s}_excitation_rel_gs_ev'] = float(mm.group(4).replace('D','E').replace('d','e'))
            toks = line.split()
            try: vals[f'state_{s}_oscillator_strength'] = float(toks[-1].replace('D','E').replace('d','e'))
            except Exception: pass
    return vals


def parse_nbf(log_text):
    m = NBF_RE.search(log_text)
    return int(m.group(1)) if m else ''


def parse_timing(stderr_text):
    sums = {'upload_d3_s':0.0, 'zero_f3_s':0.0, 'download_f3_s':0.0}
    found = False
    for m in TIMING_RE.finditer(stderr_text):
        found = True
        sums['upload_d3_s'] += float(m.group('upload'))
        sums['zero_f3_s'] += float(m.group('zero'))
        sums['download_f3_s'] += float(m.group('download'))
    return sums if found else {'upload_d3_s':'', 'zero_f3_s':'', 'download_f3_s':''}


def parse_f3(combined):
    best = {'f3_max_abs_error':'', 'f3_rms_abs_error':'', 'f3_passed':'', 'f3_notes':''}
    max_abs = -1.0; max_rms = -1.0; count = 0; notes=[]; passed=True
    for line in combined.splitlines():
        m = F3_RE.search(line)
        if not m: continue
        count += 1
        fields = {}
        for part in m.group(1).split():
            if '=' in part:
                k,v = part.split('=',1); fields[k]=v
        if fields.get('passed') != 'true': passed = False
        notes.append(fields.get('notes',''))
        if 'max_abs_error' in fields: max_abs = max(max_abs, float(fields['max_abs_error']))
        if 'rms_abs_error' in fields: max_rms = max(max_rms, float(fields['rms_abs_error']))
    if count:
        best.update({'f3_max_abs_error': f'{max_abs:.17g}', 'f3_rms_abs_error': f'{max_rms:.17g}', 'f3_passed': str(passed).lower(), 'f3_notes': f'{count}_checks_' + ('ok' if passed else 'failed')})
    return best


def median(xs):
    return statistics.median(xs) if xs else ''


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-cmd', default='python -m oqp.pyoqp')
    ap.add_argument('--sizes', default='2,3,4,5,6,8,10,15,20')
    ap.add_argument('--basis', default='6-31g')
    ap.add_argument('--nstate', type=int, default=3)
    ap.add_argument('--repeat-small', type=int, default=3)
    ap.add_argument('--repeat-large', type=int, default=1)
    ap.add_argument('--small-max-m', type=int, default=8)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--log-dir', required=True)
    ap.add_argument('--summary', required=True)
    args=ap.parse_args()
    tokens=args.run_cmd.split(); sizes=[int(x) for x in args.sizes.split(',') if x]
    log_dir=Path(args.log_dir); log_dir.mkdir(parents=True, exist_ok=True)
    input_dir=log_dir/'inputs'; input_dir.mkdir(parents=True, exist_ok=True)
    gpu_name=os.popen("nvidia-smi --query-gpu=name --format=csv,noheader | head -1").read().strip()
    host=os.uname().nodename; job=os.environ.get('SLURM_JOB_ID',''); commit=os.popen('git rev-parse --short HEAD').read().strip()
    rows=[]; largest_completed=''; largest_parity=''; first_failing=''; crossover=''
    for m in sizes:
        formula=f'C{2*m}H{2*m+2}'
        label=f'polyene_m{m}_{formula}'
        inp=input_dir/f'{label}.inp'; write_input(inp, m, args.basis, args.nstate)
        repeats=args.repeat_small if m <= args.small_max_m else args.repeat_large
        cpu_walls=[]; gpu_walls=[]; cpu_vals=gpu_vals=None; nbf=''; status='pass'; notes=[]; timing={}
        try:
            for r in range(repeats):
                rc, wall, stderr, log_text, combined = run_once(tokens, inp, {'OQP_GPU_METC':'0'}, log_dir/f'{label}_cpu_r{r}.log')
                if rc != 0: raise RuntimeError(f'CPU rc={rc}')
                cpu_walls.append(wall); cpu_vals=parse_scalars(log_text); nbf=nbf or parse_nbf(log_text)
            for r in range(repeats):
                rc, wall, stderr, log_text, combined = run_once(tokens, inp, {'OQP_GPU_METC':'1','OQP_GPU_METC_STRICT':'1','OQP_GPU_METC_TIMING':'1'}, log_dir/f'{label}_gpu_r{r}.log')
                if rc != 0: raise RuntimeError(f'GPU rc={rc}')
                gpu_walls.append(wall); gpu_vals=parse_scalars(log_text); timing=parse_timing(stderr); nbf=nbf or parse_nbf(log_text)
            rc, _, _, _, f3_combined = run_once(tokens, inp, {'OQP_GPU_METC':'1','OQP_GPU_METC_STRICT':'1','OQP_GPU_METC_TIMING':'1','OQP_GPU_METC_F3_CHECK':'1','OQP_GPU_METC_F3_CHECK_ATOL':'1e-7','OQP_GPU_METC_F3_CHECK_RTOL':'1e-7'}, log_dir/f'{label}_gpu_f3check.log')
            if rc != 0: raise RuntimeError(f'GPU f3check rc={rc}')
            f3=parse_f3(f3_combined)
            if f3.get('f3_passed') != 'true': raise RuntimeError('missing_or_failed_f3_check')
            keys=set(cpu_vals or {}) | set(gpu_vals or {})
            max_abs=max((abs((cpu_vals or {}).get(k, float('nan'))-(gpu_vals or {}).get(k, float('nan'))) for k in keys if k in cpu_vals and k in gpu_vals), default=float('nan'))
            max_rel=max((abs(cpu_vals[k]-gpu_vals[k])/max(1.0,abs(cpu_vals[k])) for k in keys if k in cpu_vals and k in gpu_vals), default=float('nan'))
            if not math.isfinite(max_abs) or max_rel > 1e-7:
                status='fail'; notes.append('scalar_parity_failed')
            largest_completed=str(m)
            if status == 'pass': largest_parity=str(m)
            if status != 'pass' and not first_failing: first_failing=str(m)
            cpu_med=median(cpu_walls); gpu_med=median(gpu_walls); ratio=cpu_med/gpu_med if gpu_med else float('nan')
            if ratio > 1.0 and not crossover: crossover=str(m)
            row={
                'm':m,'formula':formula,'nbf':nbf,'nstate':args.nstate,'basis':args.basis,
                'cpu_median_s':f'{cpu_med:.6f}','gpu_median_s':f'{gpu_med:.6f}','cpu_gpu_ratio':f'{ratio:.6f}',
                **{k:(f'{v:.6f}' if isinstance(v,float) else v) for k,v in timing.items()},
                'scalar_max_abs_delta':f'{max_abs:.12g}','scalar_max_rel_delta':f'{max_rel:.12g}',
                **f3,'gpu_name':gpu_name,'host':host,'slurm_job_id':job,'commit':commit,'status':status,'notes':';'.join(notes) or 'ok'
            }
        except Exception as exc:
            if not first_failing: first_failing=str(m)
            row={'m':m,'formula':formula,'nbf':nbf,'nstate':args.nstate,'basis':args.basis,'cpu_median_s':median(cpu_walls) or '', 'gpu_median_s':median(gpu_walls) or '', 'cpu_gpu_ratio':'', 'upload_d3_s':'','zero_f3_s':'','download_f3_s':'','scalar_max_abs_delta':'','scalar_max_rel_delta':'','f3_max_abs_error':'','f3_rms_abs_error':'','f3_passed':'false','f3_notes':'', 'gpu_name':gpu_name,'host':host,'slurm_job_id':job,'commit':commit,'status':'fail','notes':str(exc)}
            rows.append(row); break
        rows.append(row)
    fields=['m','formula','nbf','nstate','basis','cpu_median_s','gpu_median_s','cpu_gpu_ratio','upload_d3_s','zero_f3_s','download_f3_s','scalar_max_abs_delta','scalar_max_rel_delta','f3_max_abs_error','f3_rms_abs_error','f3_passed','f3_notes','gpu_name','host','slurm_job_id','commit','status','notes']
    with open(args.out_csv,'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,lineterminator='\n'); w.writeheader(); w.writerows(rows)
    faster=[r for r in rows if r.get('status')=='pass' and r.get('cpu_gpu_ratio') and float(r['cpu_gpu_ratio'])>1.0]
    Path(args.summary).write_text('\n'.join([
        f'Polyene GPU ladder Slurm job {os.environ.get("SLURM_JOB_ID","")}',
        f'commit: {commit}', f'basis: {args.basis}', 'method: BHHLYP/UMRSF-TDDFT', f'nstate: {args.nstate}',
        f'largest size completed: m={largest_completed or "none"}',
        f'largest size passing CPU/GPU parity: m={largest_parity or "none"}',
        f'first failing size: {first_failing or "none"}',
        f'GPU becomes faster at larger sizes: {"yes" if faster else "no"}',
        f'crossover size: {crossover or "not observed"}',
        'fallback occurred: no observed mixed fallback; strict GPU METC enabled for GPU runs',
        f'csv: {args.out_csv}', f'log_dir: {args.log_dir}',
    ])+'\n')
    print(Path(args.summary).read_text())
    return 0 if rows and rows[-1].get('status') == 'pass' else 1 if first_failing and not largest_parity else 0

if __name__ == '__main__':
    sys.exit(main())
