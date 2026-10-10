import collections,json,subprocess,os
from pathlib import Path
R=Path(__file__).resolve().parent.parent;E=R/'evidence';G=['git','--git-dir='+str(R/'repo.git')]
d=json.loads((E/'tips.json').read_text());refs=[x['sha'] for x in d if x['merge_base']]
env={**os.environ,'GIT_OPTIONAL_LOCKS':'0','GIT_TERMINAL_PROMPT':'0'}
with (E/'patch-errors.log').open('w') as err,(E/'patch-ids.txt').open('w') as out:
 a=subprocess.Popen(G+['log','--no-merges','--pretty=format:%H','-p','--no-ext-diff','--no-renames',*refs],stdout=subprocess.PIPE,stderr=err,env=env)
 b=subprocess.Popen(G+['patch-id','--stable'],stdin=a.stdout,stdout=out,stderr=err,env=env);a.stdout.close()
 br=b.wait();ar=a.wait()
 if ar or br:raise SystemExit(f'failed log={ar} patchid={br}')
pid={};byid=collections.defaultdict(list)
for l in (E/'patch-ids.txt').read_text().splitlines():
 p,c=l.split();pid[c]=p;byid[p].append(c)
main=set(subprocess.check_output(G+['rev-list','refs/remotes/gitlab/main'],text=True).splitlines())
mainids={pid[c] for c in main if c in pid}
for x in d:
 if not x['merge_base']:x['patch_comparison']='unrelated-history';continue
 if not x['ahead']:x['patch_comparison']='ancestor';continue
 x['patch_equivalent_main']=[c['sha'] for c in x['commits'] if pid.get(c['sha']) in mainids]
 x['patch_distinct_main']=[c['sha'] for c in x['commits'] if c['sha'] in pid and pid[c['sha']] not in mainids]
 x['merge_or_empty_commits']=[c['sha'] for c in x['commits'] if c['sha'] not in pid]
 x['patch_comparison']='distinct-patches' if x['patch_distinct_main'] else 'all-nonempty-patches-equivalent' if x['patch_equivalent_main'] else 'merge-or-empty-only'
 x['patch_set_digest']=__import__('hashlib').sha256('\n'.join(sorted({pid[c] for c in x['patch_distinct_main']})).encode()).hexdigest()
(E/'tips-enriched.json').write_text(json.dumps(d,ensure_ascii=False,indent=2))
(E/'patch-duplicates.json').write_text(json.dumps({p:cs for p,cs in byid.items() if len(cs)>1},indent=2))
print('patches',len(pid),'unique ids',len(byid),'status',dict(collections.Counter(x['patch_comparison'] for x in d)),flush=True)
