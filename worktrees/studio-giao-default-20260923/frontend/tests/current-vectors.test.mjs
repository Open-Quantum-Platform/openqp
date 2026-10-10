import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import ts from 'typescript';
import { parseHTML } from 'linkedom';

const source = readFileSync(new URL('../src/main.ts', import.meta.url), 'utf8');
const tree = ts.createSourceFile('main.ts', source, ts.ScriptTarget.Latest, true);
const names = new Set(['cubeNameFromUrl', 'currentComponentNames', 'currentVectorSettings',
  'showCurrentVectorControls', 'releaseCurrentComponents', 'loadCurrentComponents',
  'showDirectCube', 'showSurface', 'hideResultPanels']);
const selected = tree.statements.filter(s => ts.isFunctionDeclaration(s) && names.has(s.name?.text));
const functions = ts.transpileModule(selected.map(s => s.getText(tree)).join('\n'), {
  compilerOptions: { target: ts.ScriptTarget.ES2022 },
}).outputText;

function fixture() {
  // Parse the real page: controls must not inherit orbitalCard's hidden state.
  const { document } = parseHTML(readFileSync(new URL('../index.html', import.meta.url), 'utf8'));
  const element = id => { const e = document.getElementById(id); assert.ok(e, id); return e; };
  const select = (id, value) => {
    element(id).innerHTML = `<option value="${value}" selected>${value}</option>`;
  };
  select('surfacePrimary', 'water_acid.cube');
  select('surfaceSecondary', 'other.cube');
  select('surfaceOperation', 'display');
  select('surfaceSides', 'positive');
  select('moSides', 'both');
  element('currentVectorToggle').checked = true;
  element('surfaceIso').value = '0.002';
  const messages = [], fetched = [], revoked = [];
  let nextUrl = 0;
  const context = vm.createContext({
    selectedJob: 'job', volumetricRequestId: 0, activeCubeUrl: null, activeCurrentUrls: [],
    activeMapSource: null, activeDirectCube: null, activeSurfaceSelection: null,
    currentMolden: null, orbitalSources: new Map(),
    resultCubeFiles: ['water_acid.cube', 'water_jx.cube', 'water_jy.cube', 'water_jz.cube'],
    orbitalCard: element('orbitalCard'), modeCard: element('modeCard'),
    URL: { createObjectURL: () => `blob:${++nextUrl}`, revokeObjectURL: u => revoked.push(u) },
    URLSearchParams,
    $: element, pushOrbitalStyle() {}, pushArtSources() {},
    pushToResultViewer: message => messages.push(message),
    fetch: async url => {
      fetched.push(url);
      return url.includes('cube-geometry')
        ? Response.json({ xyz: '1\nwater\nH 0 0 0\n' }) : new Response('cube data');
    },
  });
  vm.runInContext(functions, context);
  return { context, document, element, select, messages, fetched, revoked };
}

function controlsVisible(f) {
  for (let e = f.element('currentVectorControls'); e; e = e.parentElement) {
    if (e.id === 'panel-analysis') break; // Tab selection is outside this controller test.
    if (e.style?.display === 'none') return false;
  }
  return true;
}

test('direct ACID opening exposes independent controls after hiding result panels', async () => {
  const f = fixture();
  f.context.hideResultPanels();
  await f.context.showDirectCube('job', '/api/jobs/job/files/water_acid.cube');
  assert.equal(f.element('orbitalCard').style.display, 'none');
  assert.equal(controlsVisible(f), true);
  assert.equal(Object.keys(f.messages.at(-1).components).length, 3);
});

test('surface display and changed isovalue retain current components and release old blobs', async () => {
  const f = fixture();
  await f.context.showSurface();
  const first = f.messages.at(-1);
  assert.equal(first.iso, 0.002);
  assert.equal(Object.keys(first.components).length, 3);
  assert.equal(controlsVisible(f), true);
  f.element('surfaceIso').value = '0.004';
  await f.context.showSurface();
  const second = f.messages.at(-1);
  assert.equal(second.iso, 0.004);
  assert.equal(Object.keys(second.components).length, 3);
  for (const url of [first.cube, ...Object.values(first.components)]) assert.ok(f.revoked.includes(url));
});

test('arithmetic surfaces and scalar-only cubes hide controls and discard components', async () => {
  const f = fixture();
  await f.context.showSurface();
  f.select('surfaceOperation', 'difference');
  await f.context.showSurface();
  assert.equal(f.messages.at(-1).components, undefined);
  assert.equal(controlsVisible(f), false);
  f.select('surfaceOperation', 'display');
  f.select('surfacePrimary', 'other.cube');
  await f.context.showSurface();
  assert.equal(f.messages.at(-1).components, undefined);
  assert.equal(controlsVisible(f), false);
});

test('toggling vectors off retains controls so they can be enabled again', async () => {
  const f = fixture();
  f.element('currentVectorToggle').checked = false;
  await f.context.showSurface();
  assert.equal(f.messages.at(-1).components, undefined);
  assert.equal(controlsVisible(f), true);
  assert.equal(f.fetched.filter(u => /_j[xyz]\.cube/.test(u)).length, 0);
  f.element('currentVectorToggle').checked = true;
  await f.context.showSurface();
  assert.equal(Object.keys(f.messages.at(-1).components).length, 3);
});

test('an older component failure cannot overwrite or revoke a newer scalar selection', async () => {
  const f = fixture();
  const originalFetch = f.context.fetch;
  let failOld;
  f.context.fetch = url => url.endsWith('_jx.cube')
    ? new Promise((_, reject) => { failOld = reject; }) : originalFetch(url);
  const old = f.context.showSurface();
  for (let i = 0; i < 100 && !failOld; i++) {
    await new Promise(resolve => setImmediate(resolve));
  }
  assert.equal(typeof failOld, "function", "surface display must request current components");
  f.select('surfacePrimary', 'other.cube');
  await f.context.showSurface();
  const latest = f.messages.at(-1);
  failOld(new Error('old download failed'));
  await old;
  assert.equal(f.messages.length, 1);
  assert.equal(f.messages.at(-1), latest);
  assert.ok(!f.revoked.includes(latest.cube));
  assert.equal(f.element('currentVectorStatus').textContent, '');
});
