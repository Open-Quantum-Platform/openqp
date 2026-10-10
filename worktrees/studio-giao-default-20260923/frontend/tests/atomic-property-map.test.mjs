import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import ts from 'typescript';

const html = readFileSync(new URL('../public/builder3d.html', import.meta.url), 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
const tree = ts.createSourceFile('builder.js', script, ts.ScriptTarget.Latest, true);
const functions = [];
function visit(node) {
  if (ts.isFunctionDeclaration(node) && ['addAtomicPropertyMap', 'drawAtomicProperty'].includes(node.name?.text)) {
    functions.push(node.getText(tree));
  }
  ts.forEachChild(node, visit);
}
visit(tree);

function fixture() {
  const representations = [], models = [], cameras = [], labels = [];
  const plugin = {
    clear: async () => {},
    canvas3d: {
      camera: { getSnapshot: () => ({ position: [1, 2, 3] }) },
      requestCameraReset: value => cameras.push(value),
    },
    builders: {
      data: { rawData: async value => { models.push(value.data); return value; } },
      structure: {
        parseTrajectory: async data => data,
        createModel: async data => data,
        createStructure: async data => data,
        tryCreateComponentStatic: async data => data,
        representation: { addRepresentation: async (_, value) => representations.push(value) },
      },
    },
  };
  const context = vm.createContext({ viewer: { plugin },
    addAtomicValueLabels: async (rows, values) => labels.push({ rows, values }) });
  vm.runInContext(functions.join('\n'), context);
  return { context, representations, models, cameras, labels };
}

const source = { format: 'xyz', data: '2\nexample\nH 0 0 0\nH 0 0 0.74\n' };

test('atomic map draws a single model with thin neutral bonds and matched scale, without default labels', async () => {
  const { context, representations: reps, models, cameras } = fixture();
  await context.drawAtomicProperty(source, [-12.5, 30]);
  assert.equal(models.length, 1);
  assert.equal(reps.length, 2);
  assert.equal(reps[0].type, 'ball-and-stick');
  assert.ok(!reps[0].typeParams.visuals.includes('element-sphere'));
  assert.equal(reps[0].color, 'uniform');
  assert.equal(reps[1].type, 'spacefill');
  const css = readFileSync(new URL('../index.html', import.meta.url), 'utf8');
  // Mol* reverses this theme's palette; the page legend must run low to high.
  const colors = Array.from(reps[1].colorParams.list.colors).reverse().map(c => `#${c.toString(16)}`);
  assert.ok(css.includes(`linear-gradient(90deg, ${colors[0]}, ${colors[1]} 50%, ${colors[2]})`));
  assert.equal(Number(models[0].split('\n')[0].slice(60, 66)), 0);
  assert.equal(Number(models[0].split('\n')[1].slice(60, 66)), 100);
  assert.deepEqual(Array.from(cameras[0].snapshot.position), [1, 2, 3]);
});

test('values can be enabled explicitly and equal values use the neutral midpoint', async () => {
  const { context, representations: reps, models, labels } = fixture();
  await context.drawAtomicProperty(source, [20, 20], true);
  assert.equal(reps.length, 2);
  assert.deepEqual(Array.from(labels[0].values), [20, 20]);
  assert.equal(labels[0].rows.length, 2);
  for (const row of models[0].split('\n').slice(0, 2)) assert.equal(Number(row.slice(60, 66)), 50);
});
