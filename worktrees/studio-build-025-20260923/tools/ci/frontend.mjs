// Match GitHub's Node 24 install, audit and build checks on native OS runners.
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const root = fileURLToPath(new URL('../../', import.meta.url));
const actual = { linux: 'linux', darwin: 'macos', win32: 'windows' }[process.platform];
if (process.env.STUDIO_CI_PLATFORM !== actual) {
  throw new Error(`Mismatched runner: expected ${process.env.STUDIO_CI_PLATFORM}, found ${actual}`);
}
if (Number(process.versions.node.split('.')[0]) !== 24) {
  throw new Error('Studio frontend CI requires Node 24, matching GitHub CI');
}
console.log(`Frontend CI: ${actual}, Node ${process.version}`);
function run(command, args, cwd) {
  // npm is a .cmd launcher on Windows. All command arguments below are fixed.
  const result = spawnSync(command, args, {
    cwd, stdio: 'inherit', shell: process.platform === 'win32',
  });
  if (result.error) throw result.error;
  if (result.status !== 0) process.exit(result.status || 1);
}
const frontend = path.join(root, 'frontend');
run('npm', ['ci', '--no-audit', '--no-fund'], frontend);
run('node', ['--test', 'frontend/tests/responsiveness.test.mjs',
  'frontend/tests/current-vectors.test.mjs', 'frontend/tests/nmr-workflows.test.mjs',
  'frontend/tests/atomic-property-map.test.mjs'], root);
run('npm', ['audit', '--audit-level=moderate'], frontend);
run('npm', ['run', 'build'], frontend);
