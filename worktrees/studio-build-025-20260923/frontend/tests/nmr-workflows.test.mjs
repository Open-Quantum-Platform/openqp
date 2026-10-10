import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import ts from 'typescript';
import { parseHTML } from 'linkedom';

const source = readFileSync(new URL('../src/main.ts', import.meta.url), 'utf8');
const tree = ts.createSourceFile('main.ts', source, ts.ScriptTarget.Latest, true);
const names = new Set(['ALL_THEORIES', 'DFT_RESPONSE', 'LABELED_RESPONSE', 'MRSF_ONLY',
  'GRADIENT_THEORIES', 'WORKFLOWS', 'currentWf']);
const functions = new Set(['selectWorkflow', 'syncNmrControls', 'syncWorkflowOptions',
  'buildWorkflowSelect', 'generateInp', 'fieldValue']);
const selected = tree.statements.filter(s =>
  (ts.isFunctionDeclaration(s) && functions.has(s.name?.text)) ||
  (ts.isVariableStatement(s) && s.declarationList.declarations.some(d => names.has(d.name.getText(tree)))));
const code = ts.transpileModule(selected.map(s => s.getText(tree)).join('\n'), {
  compilerOptions: { target: ts.ScriptTarget.ES2022 },
}).outputText;

function fixture() {
  const { document, window } = parseHTML(readFileSync(new URL('../index.html', import.meta.url), 'utf8'));
  // Linkedom does not implement the browser's select.value setter.
  for (const select of document.querySelectorAll('select')) {
    Object.defineProperty(select, 'value', {
      configurable: true,
      get() { return this.querySelector('option[selected]')?.value ?? this.querySelector('option')?.value ?? ''; },
      set(value) {
        for (const option of this.querySelectorAll('option')) {
          option.toggleAttribute('selected', option.value === value);
        }
      },
    });
  }
  const $ = id => { const e = document.getElementById(id); assert.ok(e, id); return e; };
  $('nmrAcid').checked = false;
  const context = vm.createContext({
    document, $, theorySel: $('theory'), optionsCard: $('optionsCard'),
    hessTypeSel: $('hessType'), xyzArea: { value: 'H 0 0 0' },
    parseAtoms: () => [['H', 0, 0, 0]], currentBasis: () => 'sto-3g',
    currentFunctional: () => 'pbe', inputValidation: () => null,
    optionList: () => '', SCF_DEFAULTS: {}, pdbSource: null,
    syncPcmReference() {}, syncWorkflowDetails() {}, updateInpPreview() {},
  });
  vm.runInContext(code, context);
  context.syncFieldStates = context.syncNmrControls;
  const run = script => vm.runInContext(script, context);
  run('buildWorkflowSelect()');
  return { $, run, window, choose: key => run(`selectWorkflow(WORKFLOWS.find(w => w.key === ${JSON.stringify(key)}))`) };
}

test('ACID and AICD searches reveal a selectable dedicated calculation', () => {
  const { $, window, choose, run } = fixture();
  for (const query of ['ACID', 'AICD']) {
    $('workflowSearch').value = query;
    $('workflowSearch').dispatchEvent(new window.Event('input'));
    assert.ok($('workflowSel').querySelector('option[value="acid"]'));
  }
  choose('acid');
  assert.equal($('nmrAcid').checked, true);
  assert.equal($('nmrAcid').disabled, true);
  assert.equal($('nmrGauge').value, 'giao');
  assert.equal($('nmrGauge').disabled, true);
  assert.ok($('nmrAcid').closest('.wf-opt').classList.contains('on'));
  assert.equal($('acidGridRow').style.display, '');
  assert.match(run('generateInp()'), /\nnmr\(gauge=giao,acid=true,acid_spacing=0\.2,acid_padding=5\.0\)\n/);
});

test('NMR starts with GIAO and returning from ACID restores optional controls', () => {
  const { $, choose, run } = fixture();
  choose('nmr');
  assert.match(run('generateInp()'), /\nnmr\(gauge=giao\)\n/);
  choose('acid');
  choose('nmr');
  assert.equal($('nmrAcid').checked, false);
  assert.equal($('nmrAcid').disabled, false);
  assert.equal($('nmrGauge').disabled, false);
  assert.equal($('acidGridRow').style.display, 'none');
  $('nmrGauge').value = 'cgo';
  assert.match(run('generateInp()'), /\nnmr\(gauge=cgo\)\n/);
  assert.doesNotMatch(run('generateInp()'), /acid=true/);
});

test('restoring inconsistent ACID recipe controls cannot produce CGO or omit ACID', () => {
  const { $, choose, run } = fixture();
  choose('acid');
  $('nmrGauge').value = 'cgo';
  $('nmrAcid').checked = false;
  run('syncNmrControls()');
  assert.equal($('nmrGauge').value, 'giao');
  assert.equal($('nmrAcid').checked, true);
  assert.match(run('generateInp()'), /nmr\(gauge=giao,acid=true/);
});
