import collections,csv,html,json,re
from pathlib import Path
from urllib.parse import quote

R=Path(__file__).resolve().parent;E=R.parent/'evidence'
d=json.loads((E/'review-index.json').read_text());refs=json.loads((E/'refs.json').read_text());snap=json.loads((E/'snapshot.json').read_text());worktrees=json.loads((E/'local-worktree-status.json').read_text());bysha={x['sha']:x for x in d}
labels={'MAIN_ANCESTOR':'main 이력에 포함','EXACT_TIP_MERGED_PR':'동일 끝점 PR 병합 확인','PATCH_EQUIVALENT_CHECK_MERGES':'non-merge patch 일치·merge 검토','REQUIRES_REVIEW':'고유 변경·반영 여부 검토','SEPARATE_HISTORY':'별도 이력'}
categories={'GPU':'GPU','ERI':'ERI / integrals','NMR':'NMR','REKS':'REKS','QMRSF':'QMRSF / higher spin','HESSIAN':'Hessian / frequency','NAC_NAMD':'NAC / NAMD','SOC_RELATIVITY':'SOC / SSC / X2C','EKT_SPECTROSCOPY':'EKT / PBC / spectroscopy','WAVEFUNCTION':'Wavefunction methods','SOLVATION':'Solvation / PCM','DFTB_XTB_QMMM':'DFTB / xTB / QM/MM','OPTIMIZATION':'Geometry optimization','SCF_RESPONSE':'SCF / response','DFT_GRID':'DFT / grid','SYMMETRY':'Symmetry','BUILD_CI':'Build / CI / release','DOCS_API':'Docs / API / input','ARCHIVE':'Historical documents','OTHER':'Other / mixed'}
q=[]
for p in sorted((E/'qwen').glob('[0-9][0-9][0-9].json')):q.extend(json.loads(p.read_text())['records'])
qby={x['id']:x for x in q}
for x in d:
 if x['sha'] in qby:
  x['category']=qby[x['sha']]['category'];x['category_basis']='Qwen metadata classification; topical index only'
 # Avoid substring errors such as SSC inside CASSCF. Explicit scientific branch names win.
 names=' '.join(a['name'].lower() for a in x['aliases'])
 if re.search(r'casscf|caspt2|nevpt2|ccsd',names):x['category']='WAVEFUNCTION'
 if re.search(r'analytic-ht-transition|mrsf-2pa',names):x['category']='EKT_SPECTROSCOPY'
 if re.search(r'rys-hessian',names):x['category']='HESSIAN'
 for pattern,category in [
  (r'v130-exclude-private-backends','BUILD_CI'),
  (r'xc-grid-gradient','DFT_GRID'),
  (r'nmr-integral-response-validation','NMR'),
  (r'claude/determined-galileo','WAVEFUNCTION'),
  (r'space-separated-method-basis|export-mo-frequency-formats|export-mrsf-analysis','DOCS_API'),
  (r'mrsf-h2-zero-closed','SCF_RESPONSE'),
  (r'odp-umbrella-native','NAC_NAMD'),
 ]:
  if re.search(pattern,names):x['category']=category
 x['display_name']=next((a['name'] for a in x['aliases'] if a['source']=='gitlab' and not a['archive'] and not a['name'].startswith('github-')),x['aliases'][0]['name'])

def esc(s):return str(s).replace('|','\\|').replace('\n',' ')
def url(a,sha=None):
 if a['source']=='gitlab':base='https://qchemlab.knu.ac.kr/open-quantum-platform/internal/openqp/-/tree/'
 elif a['source']=='github-personal':base='https://github.com/karmachoi/openqp/tree/'
 elif a['source']=='github-private':base='https://github.com/karmachoi/openqp-private/tree/'
 else:return None
 return base+(sha or quote(a['name'],safe=''))

(R/'branches').mkdir(exist_ok=True)
for x in d:
 s=[f"# {x['display_name']}",f"\nSHA: `{x['sha']}`  ",f"판정: **{labels[x['disposition']]}**  ",f"분야: {categories[x['category']]} (검색용 분류)  ",f"마지막 commit: {x['commit_date']} / {x['author']}  ",f"제목: {x['subject']}\n",'[검색 가능한 전체 목록](../index.html) · [조사 결과](../README.md)\n','## 같은 끝점을 가리키는 모든 원본\n']
 for a in x['aliases']:
  u=url(a,x['sha']);txt=f"{a['source']} / {a['name']}"
  s.append(f'- [{txt}]({u})' if u else f'- `{txt}`')
 s+=['\n## 현재 GitLab main과의 비교\n',f"- 기준 main: `{snap['main']}`",f"- 공통 조상: `{x['merge_base'] or '없음'}`",f"- 앞선 커밋 {x['ahead']} / 뒤처진 커밋 {x['behind']}",f"- non-merge patch: main과 일치 {len(x.get('patch_equivalent_main',[]))}, 다름 {len(x.get('patch_distinct_main',[]))}",f"- merge/empty 등 patch 비교 제외: {len(x.get('merge_or_empty_commits',[]))}",f"- GitLab 어느 브랜치에서도 끝점 도달 가능: {x['in_gitlab']}",f"- GitLab 전체에서 동일 patch를 못 찾은 커밋: {len(x['patches_absent_gitlab'])}",f"- 공통 조상 이후 변경 파일: {len(x['feature_files'])}; 그중 현재 main과 동일 {len(x['feature_files_equal_main'])}, 다름 {len(x['feature_files_different_main'])}", '\n앞선 커밋 수에는 inherited private 작업이 포함될 수 있다. 개수만으로 미반영 기능 수를 판단하지 않는다.\n']
 if x['prs_exact_tip'] or x['mrs_exact_tip']:
  s.append('## 이 끝점과 정확히 일치하는 PR/MR\n')
  for p in x['prs_exact_tip']:s.append(f"- [{p['source']} PR #{p['number']}: {p['title']}]({p['url']}) — state={p['state']}; merged={p['merged_at'] or '아니오'}; merge commit in main={p['merge_in_main']}")
  for p in x['mrs_exact_tip']:s.append(f"- [GitLab MR !{p['iid']}: {p['title']}]({p['url']}) — {p['state']}")
 if x['prs_same_name']:
  s.append('\n## 같은 이름의 PR — 끝점이 다르므로 전체 반영 증거가 아님\n')
  for p in x['prs_same_name']:s.append(f"- [{p['source']} #{p['number']}: {p['title']}]({p['url']}) — {p['state']}; PR head `{p['head_sha'][:12]}`; merged={p['merged_at'] or '아니오'}")
 if x['same_category_descendants']:
  s.append('\n## 이 끝점을 포함하는 같은 분야의 후속 브랜치 후보\n')
  for y in x['same_category_descendants']:s.append(f"- [{', '.join(y['aliases'][:5])}]({y['sha'][:12]}.md) — 추가 커밋 {y['distance']}")
  s.append('\n이 목록은 ancestry 관계이며 기능 대체나 최선의 통합 대상을 보증하지 않는다.')
 if x['same_remaining_patch_set']:
  s.append('\n## 남은 non-merge patch 집합이 같은 다른 끝점\n')
  for sha in x['same_remaining_patch_set']:s.append(f'- [{bysha[sha]["display_name"] if "display_name" in bysha[sha] else sha[:12]}]({sha[:12]}.md)')
 s+=['\n## 추가 커밋 전체\n','| SHA | 날짜 | main patch 판정 | GitLab 전체 patch | 제목 |','| --- | --- | --- | --- | --- |']
 cs=x['commits'] if x['merge_base'] else x['commits'][:20]
 for c in cs:
  cp='일치' if c['sha'] in x.get('patch_equivalent_main',[]) else '다름' if c['sha'] in x.get('patch_distinct_main',[]) else 'merge/empty/비교 제외'
  gp='미발견' if c['sha'] in x['patches_absent_gitlab'] else '다른 SHA의 동일 patch' if c['sha'] in x['patches_equivalent_gitlab'] else 'commit 보존 또는 비교 제외'
  s.append(f"| `{c['sha'][:12]}` | {c['date'][:10]} | {cp} | {gp} | {esc(c['subject'])} |")
 if not x['merge_base'] and len(x['commits'])>20:s.append(f'\n별도 문서 이력 {len(x["commits"])}개 중 최근 20개만 표시. 전체 목록은 evidence/tips.json에 보존했다.')
 s+=['\n## 변경 파일 전체\n','| 상태 | 경로 | 현재 main과 내용 |','| --- | --- | --- |']
 eq=set(x['feature_files_equal_main'])
 for f in x['feature_files']:s.append(f"| {f['status']} | `{esc(f['path'])}` | {'동일' if f['path'] in eq else '다름'} |")
 local=[w for w in worktrees if w.get('HEAD')==x['sha']]
 if local:
  s.append('\n## 이 끝점을 사용 중인 로컬 worktree\n')
  for w in local:s.append(f"- `{w['worktree']}` — exists={w['exists']}, status 항목 {len(w.get('status_entries',[]))}; {w.get('prunable','')}")
 (R/'branches'/f"{x['sha'][:12]}.md").write_text('\n'.join(s)+'\n')

with (R/'all-branches.csv').open('w',newline='') as f:
 fields=['source','name','sha','category','disposition','in_gitlab','ahead','behind','patches_different_main','patches_absent_gitlab','modified_files','last_commit_date','subject','detail']
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for a in refs:
  x=bysha[a['sha']];w.writerow(dict(source=a['source'],name=a['name'],sha=x['sha'],category=x['category'],disposition=x['disposition'],in_gitlab=x['in_gitlab'],ahead=x['ahead'],behind=x['behind'],patches_different_main=len(x.get('patch_distinct_main',[])),patches_absent_gitlab=len(x['patches_absent_gitlab']),modified_files=len(x['feature_files']),last_commit_date=x['commit_date'],subject=x['subject'],detail='branches/'+x['sha'][:12]+'.md'))

local=['# GitLab에 끝점이 없는 로컬 브랜치\n','같은 내용의 다른 patch, 후속 브랜치의 변형, merge commit일 수 있다. 미발견만으로 삭제·이동하지 않는다.\n','| 원본 | 브랜치 | SHA | GitLab 미발견 patch | 상세 |','| --- | --- | --- | --- | --- |']
for a in refs:
 if a['source'].startswith('local-') and not bysha[a['sha']]['in_gitlab']:
  x=bysha[a['sha']];local.append(f"| {a['source']} | `{a['name']}` | `{a['sha'][:12]}` | {len(x['patches_absent_gitlab'])} | [상세](branches/{a['sha'][:12]}.md) |")
local+=['\n## 로컬 원본 경로\n']
srcs={r['source']:r['common_dir'] for r in json.loads((E/'local-refs.json').read_text())}
for k,v in srcs.items():local.append(f'- `{k}`: `{v}`')
(R/'local-only.md').write_text('\n'.join(local)+'\n')
ws=['# 로컬 worktree 상태\n','tracked/untracked 경로만 확인했으며 파일을 변경하지 않았다. 상태는 조사 시점의 스냅샷이다.\n']
for w in worktrees:
 ws+=[f"## {w['worktree']}\n",f"- 원본 `{w['source']}`, branch `{w.get('branch','detached')}`, HEAD `{w.get('HEAD','')}`",f"- 존재: {w['exists']}; status 항목: {len(w.get('status_entries',[]))}"]
 if w.get('prunable'):ws.append('- Git 표시: '+w['prunable'])
 if w.get('status_entries'):ws+=['\n```text',*w['status_entries'],'```']
(R/'worktrees.md').write_text('\n'.join(ws)+'\n')

counts=collections.Counter(x['disposition'] for x in d)
summary=['## 전체 집계\n',f"스냅샷 UTC: `{snap['utc']}`  ",f"비교 main: `{snap['main']}`\n",'| 조사 대상 | 브랜치/ref 수 |','| --- | ---: |','| GitHub karmachoi/openqp | 244 |','| GitHub karmachoi/openqp-private | 152 |','| GitLab internal/openqp | 390 (개발·일반 379, PR 보존 11) |','| Ultra 로컬 16개 common directory | 326 |',f'| 합계 / 고유 SHA | {len(refs)} / {len(d)} |','\n| 고유 끝점 판정 | 개수 |','| --- | ---: |']
for k,n in counts.items():summary.append(f'| {labels[k]} | {n} |')
summary+=['\n## 상세 자료\n','- [검색·필터 가능한 전체 조사표](index.html)','- [원본별 전체 브랜치 CSV](all-branches.csv)','- [GitLab 미보존 로컬 브랜치](local-only.md)','- [215개 등록 worktree 상태](worktrees.md)','- [원시 스냅샷·patch·PR 증거](../evidence/)','\n## 분야별 브랜치 찾아보기\n']
for cat in categories:
 rows=sorted((x for x in d if x['category']==cat),key=lambda x:x['display_name'])
 if not rows:continue
 summary.append(f'### {categories[cat]} ({len(rows)})\n')
 for x in rows:summary.append(f"- [{x['display_name']} · {x['sha'][:8]}](branches/{x['sha'][:12]}.md) — {labels[x['disposition']]}")
summary.append(f'\n분야 분류는 검색을 돕기 위한 것이다. {len(q)}개 이름 불명확 항목에 로컬 Qwen의 구조 검증된 분류를 참고했으며, SHA·개수·병합 판정은 모두 Git/API에서 계산했다. 분류 자체는 과학적 검증이 아니다.')
text=(R/'findings.md').read_text();split=text.index('## 조사 방법')
(R/'README.md').write_text(text[:split]+'\n'+'\n'.join(summary[:summary.index('\n## 분야별 브랜치 찾아보기\n')])+'\n\n'+text[split:]+'\n'+'\n'.join(summary[summary.index('\n## 분야별 브랜치 찾아보기\n'):]))

data=[]
for x in d:data.append(dict(sha=x['sha'],name=x['display_name'],category=categories[x['category']],status=labels[x['disposition']],in_gitlab=x['in_gitlab'],ahead=x['ahead'],behind=x['behind'],missing=len(x['patches_absent_gitlab']),files=len(x['feature_files']),date=x['commit_date'][:10],subject=x['subject'],aliases=[a['source']+'/'+a['name'] for a in x['aliases']],paths=[f['path'] for f in x['feature_files']],detail='branches/'+x['sha'][:12]+'.md'))
payload=json.dumps(data,ensure_ascii=False).replace('<','\\u003c')
page='''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>OpenQP branch audit</title><style>
:root{color-scheme:dark}body{font:15px/1.55 system-ui,-apple-system,sans-serif;background:#10151e;color:#dce5f3;max-width:1600px;margin:auto;padding:28px}h1{font-size:28px;margin-bottom:6px}a{color:#87baff}small,.muted{color:#a6b4c9}.cards{display:flex;gap:16px;flex-wrap:wrap;margin:22px 0}.card{background:#1a2535;border:1px solid #34465d;border-radius:10px;padding:14px 22px}.card b{font-size:25px;display:block}input,select{font:inherit;color:inherit;background:#1b2738;border:1px solid #425675;border-radius:6px;padding:9px}input{min-width:340px}nav{display:flex;gap:12px;flex-wrap:wrap;margin:18px 0}table{width:100%;border-collapse:collapse}th{position:sticky;top:0;background:#24334a;text-align:left}td,th{padding:10px;border-bottom:1px solid #2e3c50;vertical-align:top}td:first-child{max-width:450px;overflow-wrap:anywhere}.sha{font-family:monospace;color:#9eafc7}details{font-size:12px;max-width:490px}summary{cursor:pointer;color:#98bfff}.warn{color:#ffd287}.note{background:#1a2535;border-left:4px solid #749edf;padding:12px 18px;margin:20px 0}button{padding:8px;background:#30496c;color:white;border:0;border-radius:5px;cursor:pointer}</style>
<h1>OpenQP 개발 브랜치 전수 조사</h1><div class="muted">2026-09-22 · GitHub personal/private + GitLab internal + Ultra local · 원본 변경 없음</div>
<nav><a href="README.md">세부 조사 보고서</a><a href="all-branches.csv">전체 CSV</a><a href="local-only.md">로컬 미보존 브랜치</a><a href="worktrees.md">worktree 상태</a></nav>
<div class="cards"><div class="card"><b>1,112</b>원본별 refs</div><div class="card"><b>533</b>고유 끝점</div><div class="card"><b>210</b>main 이력 / 동일 PR 병합</div><div class="card"><b>62</b>GitLab에 끝점 없는 SHA</div></div>
<div class="note">REKS 최신 끝점은 GitLab의 가장 앞선 보존본보다 14 commits 앞섭니다. 미발견 patch는 동일 변경을 못 찾았다는 뜻이며, 독립 기능 수가 아닙니다. 이 조사는 전체 branch의 Git/patch/PR 검토이며 전체 수학·분자 계산 검증은 아닙니다.</div>
<nav><input id="q" placeholder="branch, SHA, 파일 경로 검색" aria-label="검색"><select id="cat" aria-label="분야"><option value="">모든 분야</option></select><select id="status" aria-label="판정"><option value="">모든 판정</option></select><label><input type="checkbox" id="missing" style="min-width:0"> GitLab에 끝점 없음</label><button id="reset">초기화</button></nav><p id="count"></p>
<table><thead><tr><th>브랜치 · SHA · 원본</th><th>분야 / 판정</th><th>main<br>앞/뒤</th><th>변경<br>파일</th><th>GitLab 미발견<br>patch</th><th>최종 commit</th></tr></thead><tbody id="rows"></tbody></table>
<script>const DATA=__DATA__;
const $=id=>document.getElementById(id);function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n}
for(const [id,key] of [['cat','category'],['status','status']])for(const value of [...new Set(DATA.map(x=>x[key]))].sort()){const o=el('option',value);o.value=value;$(id).append(o)}
function draw(){const q=$('q').value.toLowerCase();const rows=DATA.filter(x=>(!$('cat').value||x.category===$('cat').value)&&(!$('status').value||x.status===$('status').value)&&(!$('missing').checked||!x.in_gitlab)&&(!q||[x.name,x.sha,x.subject,...x.aliases,...x.paths].join(' ').toLowerCase().includes(q)));$('count').textContent=`${rows.length} / ${DATA.length} 고유 끝점 표시 · 각 상세 문서에 전체 commit/file 목록 포함`;$('rows').replaceChildren();for(const x of rows){const tr=el('tr');const a=el('a',x.name);a.href=x.detail;const td=el('td');td.append(a,el('div',x.sha.slice(0,12),'sha'));const de=el('details');de.append(el('summary',`${x.aliases.length}개 원본 · ${x.subject}`));for(const s of x.aliases)de.append(el('div',s));td.append(de);tr.append(td);const c=el('td');c.append(el('div',x.category),el('small',x.status));tr.append(c,el('td',`${x.ahead} / ${x.behind}`),el('td',x.files),el('td',x.missing,x.missing?'warn':''),el('td',x.date));$('rows').append(tr)}}
for(const id of ['q','cat','status','missing'])$(id).addEventListener('input',draw);$('reset').onclick=()=>{$('q').value='';$('cat').value='';$('status').value='';$('missing').checked=false;draw()};draw();</script></html>'''
(R/'index.html').write_text(page.replace('__DATA__',payload))
(R/'report-data.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
print('rendered',len(d),'detail pages,',len(refs),'CSV rows; Qwen records',len(q))
