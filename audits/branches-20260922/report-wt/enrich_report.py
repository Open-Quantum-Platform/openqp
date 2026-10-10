import collections,concurrent.futures,json,re,subprocess
from pathlib import Path
R=Path(__file__).resolve().parent.parent;E=R/'evidence';G=['git','--git-dir='+str(R/'repo.git')]
d=json.loads((E/'tips-enriched.json').read_text());refs=json.loads((E/'refs.json').read_text());bysha={x['sha']:x for x in d}
main=set(subprocess.check_output(G+['rev-list','refs/remotes/gitlab/main'],text=True).splitlines())
rules=[('GPU',r'gpu|metc|routec|rot-seam'),('ERI',r'ispher|\brys\b|integral|int2|perf/int|spherical|libint|rot-rys'),('NMR',r'nmr|giao|shielding'),('REKS',r'reks'),('QMRSF',r'qmrsf|doublet.quartet|dk-|dk/|xqmrsf|s3r|fth-dc'),('HESSIAN',r'hess|curvature|frequency|s1-freq'),('NAC_NAMD',r'nac|namd|fssh|hopping|zhu.nakamura|state.propagation|rt-mrsf|uracil.*campaign'),('SOC_RELATIVITY',r'soc|spin.spin|ssc|x2c|zfs|umrsf'),('EKT_SPECTROSCOPY',r'ekt|dyson|pecd|photo|gelius|xas|spectro|spec-|ip-ea|ea-fock'),('WAVEFUNCTION',r'casscf|caspt2|nevpt|ccsd|fci|wf-|mp2|ncsf|ci.phase|ci.irrep|ci.spin|quantum.comput'),('SOLVATION',r'pcm|ddx|solvent|cosmo'),('DFTB_XTB_QMMM',r'dftb|dtcam|xtb|qmmm|droplet|espf'),('OPTIMIZATION',r'optim|meci|mecp|geometric|dlc|neb|coord.degeneracy'),('SCF_RESPONSE',r'scf|guess|trah|davidson|z.vector|zvector|response.phase|spin.label|mrsf.fock|higher.root|gradient'),('DFT_GRID',r'xc|dft.grid|sg1|sg.grids'),('SYMMETRY',r'symmetry|petite|irrep'),('BUILD_CI',r'build|ci/|/ci-|/ci\b|blas|ilp64|lp64|docker|macos|windows|intel|wheel|external|release|licen|enforcement|policy|test.on.CI|review|auto.request|concurrency|preflight'),('DOCS_API',r'docs/|readme|input|pythonic|api|log-|logging|molden|export|python.version|banner|author'),('ARCHIVE',r'gh-pages|karmachoi/README|backup|archive')]
def category(x):
 names=' '.join(a['name'].lower() for a in x['aliases'])
 for c,p in rules:
  if re.search(p,names,re.I):return c,'branch-name'
 text=' '.join(c['subject'] for c in x['commits'][:10])+' '+x['subject']
 for c,p in rules:
  if re.search(p,text,re.I):return c,'commit-subject'
 return 'OTHER','unclassified'
prs=[]
for fn,source in [('github-upstream-prs.json','upstream'),('github-personal-prs.json','personal'),('github-private-prs.json','private')]:
 for p in json.loads((E/fn).read_text()):
  prs.append(dict(source=source,number=p['number'],title=p['title'],url=p['html_url'],state=p['state'],merged_at=p['merged_at'],merge_sha=p.get('merge_commit_sha'),head_sha=p['head']['sha'],head_branch=p['head']['ref'],head_repo=p['head']['repo']['full_name'] if p['head']['repo'] else None,merge_in_main=p.get('merge_commit_sha') in main))
mrs=json.loads((E/'gitlab-metadata.json').read_text())['mrs']
for x in d:
 x['category'],x['category_basis']=category(x)
 x['prs_exact_tip']=[p for p in prs if p['head_sha']==x['sha']]
 names={a['name'].removeprefix('github-private/').removeprefix('github-karmachoi/') for a in x['aliases']}
 x['prs_same_name']=[p for p in prs if p['head_branch'] in names and p not in x['prs_exact_tip'] and p['head_repo'] in ['karmachoi/openqp','karmachoi/openqp-private']]
 x['mrs_exact_tip']=[dict(iid=p['iid'],title=p['title'],state=p['state'],url=p['web_url']) for p in mrs if p.get('sha')==x['sha']]
 if x['in_main']:x['disposition']='MAIN_ANCESTOR'
 elif any(p['merged_at'] and p['merge_in_main'] for p in x['prs_exact_tip']):x['disposition']='EXACT_TIP_MERGED_PR'
 elif x['patch_comparison']=='all-nonempty-patches-equivalent':x['disposition']='PATCH_EQUIVALENT_CHECK_MERGES'
 elif not x['merge_base']:x['disposition']='SEPARATE_HISTORY'
 else:x['disposition']='REQUIRES_REVIEW'
def anc(x):return x['sha'],set(subprocess.check_output(G+['rev-list',x['sha']],text=True).splitlines())
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:ancestry=dict(ex.map(anc,d))
for x in d:
 supers=[y for y in d if y['sha']!=x['sha'] and x['sha'] in ancestry[y['sha']] and y['category']==x['category']]
 supers.sort(key=lambda y:len(ancestry[y['sha']]-ancestry[x['sha']]))
 x['same_category_descendants']=[dict(sha=y['sha'],aliases=[a['source']+'/'+a['name'] for a in y['aliases']],distance=len(ancestry[y['sha']]-ancestry[x['sha']])) for y in supers[:4]]
groups=collections.defaultdict(list)
for x in d:
 if x.get('patch_distinct_main'):groups[x['patch_set_digest']].append(x['sha'])
for x in d:x['same_remaining_patch_set']=[s for s in groups.get(x.get('patch_set_digest'),[]) if s!=x['sha']]
(E/'review-index.json').write_text(json.dumps(d,ensure_ascii=False,indent=2));(E/'prs-normalized.json').write_text(json.dumps(prs,ensure_ascii=False,indent=2))
print('disposition',dict(collections.Counter(x['disposition'] for x in d)))
print('categories',dict(collections.Counter(x['category'] for x in d)))
amb=[x for x in d if x['category_basis']!='branch-name']
print('ambiguous categories',len(amb))
manifest=[]
for x in amb:
 cs=[c['subject'] for c in x['commits'] if c['sha'] in x.get('patch_distinct_main',[])][:5]
 manifest.append(dict(id=x['sha'],branch=x['aliases'][0]['name'],subjects=cs or [x['subject']],files=[f['path'] for f in x['feature_files'][:8]]))
(E/'qwen-input.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
