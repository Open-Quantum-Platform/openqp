#!/usr/bin/env python3
"""Read-only ref inventory. Writes evidence only under this audit directory."""
import collections, concurrent.futures, datetime, hashlib, json, os, re, subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
REPO=ROOT/'repo.git'
OUT=ROOT/'evidence'
MAIN='refs/remotes/gitlab/main'
ENV={**os.environ,'GIT_OPTIONAL_LOCKS':'0','GIT_TERMINAL_PROMPT':'0'}
def git(*args,check=True,timeout=180):
    p=subprocess.run(['git','--git-dir='+str(REPO),*args],text=True,capture_output=True,env=ENV,timeout=timeout)
    if check and p.returncode:raise RuntimeError(p.stderr[:1000])
    return p.stdout.strip() if check else p
def save(name,data):
    (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2))

refs=[]
for line in git('for-each-ref','--format=%(refname)\t%(objectname)','refs/remotes').splitlines():
    ref,sha=line.split('\t');source=ref.split('/')[2];name='/'.join(ref.split('/')[3:])
    refs.append(dict(ref=ref,source=source,name=name,sha=sha,archive=bool(re.search(r'(?:^|/)(?:github(?:-private)?/)?pull/\d+/',name))))
save('refs.json',refs)
main=git('rev-parse',MAIN)
bysha=collections.defaultdict(list)
for r in refs:bysha[r['sha']].append(r)
main_tree={}
for line in git('ls-tree','-r',main).splitlines():
    meta,path=line.split('\t',1);main_tree[path]=meta.split()[2]
main_ancestors=set(git('rev-list',main).splitlines())
gitlab_tips=[r['sha'] for r in refs if r['source']=='gitlab']
reachable_gitlab=set(git('rev-list',*sorted(set(gitlab_tips))).splitlines())
save('snapshot.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),main=main,sources=dict(collections.Counter(r['source'] for r in refs)),unique_tips=len(bysha)))

def inspect(item):
    sha,aliases=item
    cached=OUT/'tips'/f'{sha}.json'
    if cached.exists():return json.loads(cached.read_text())
    fields=git('show','-s','--format=%H%x00%T%x00%aI%x00%cI%x00%an%x00%s',sha).split('\0')
    d=dict(sha=sha,tree=fields[1],author_date=fields[2],commit_date=fields[3],author=fields[4],subject=fields[5],aliases=aliases,in_main=sha in main_ancestors,in_gitlab=sha in reachable_gitlab)
    base=git('merge-base',main,sha,check=False)
    d['merge_base']=base.stdout.strip() if base.returncode==0 else None
    d['ahead']=int(git('rev-list','--count',f'{main}..{sha}'))
    d['behind']=int(git('rev-list','--count',f'{sha}..{main}'))
    raw=git('log','--format=%H%x09%aI%x09%an%x09%s',f'{main}..{sha}')
    d['commits']=[dict(zip(['sha','date','author','subject'],l.split('\t',3))) for l in raw.splitlines()]
    if d['merge_base']:
        ns=git('diff','--no-renames','--name-status',d['merge_base'],sha)
        d['feature_files']=[dict(status=l.split('\t')[0],path=l.split('\t')[1]) for l in ns.splitlines()]
        branch_tree={}
        for line in git('ls-tree','-r',sha).splitlines():
            meta,path=line.split('\t',1);branch_tree[path]=meta.split()[2]
        d['feature_files_equal_main']=[x['path'] for x in d['feature_files'] if branch_tree.get(x['path'])==main_tree.get(x['path'])]
        d['feature_files_different_main']=[x['path'] for x in d['feature_files'] if branch_tree.get(x['path'])!=main_tree.get(x['path'])]
    else:
        d['feature_files']=[];d['feature_files_equal_main']=[];d['feature_files_different_main']=[]
    d['patch_comparison']='pending' if d['ahead'] else 'ancestor'
    cached.write_text(json.dumps(d,ensure_ascii=False,indent=2))
    return d
(OUT/'tips').mkdir(exist_ok=True)
tips=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
    for n,d in enumerate(ex.map(inspect,bysha.items()),1):
        tips.append(d)
        if n%40==0:print('inspected',n,'/',len(bysha),flush=True)
save('tips.json',tips)
remote=[r for r in refs if r['source'].startswith('github-')]
missing=[r for r in remote if r['sha'] not in reachable_gitlab]
local_missing=[r for r in refs if r['source'].startswith('local-') and r['sha'] not in reachable_gitlab]
save('missing-remote-refs.json',missing);save('missing-local-refs.json',local_missing)
print(json.dumps(dict(tips=len(tips),ancestor_tips=sum(d['in_main'] for d in tips),missing_remote=len(missing),missing_local=len(local_missing),without_common_base=sum(not d['merge_base'] for d in tips),ahead_buckets=dict(collections.Counter('0' if d['ahead']==0 else '1-100' if d['ahead']<=100 else '101-500' if d['ahead']<=500 else '>500' for d in tips))),indent=2),flush=True)
