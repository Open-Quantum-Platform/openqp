import assert from 'node:assert/strict';
import { test } from 'node:test';
import { pixelRatio, rayScale, rayStopReason, MAX_RENDER_PIXELS, MAX_RAY_PIXELS } from '../src/render-policy.ts';
import { withAbort } from '../src/api.ts';

test('4K/Retina rendering stays within pixel budgets', () => {
  for (const [width, height, dpi] of [[3840, 2160, 2], [7680, 4320, 3], [800, 600, 1]]) {
    const ratio = pixelRatio(width, height, dpi);
    const pixels = width * height * ratio ** 2;
    assert.ok(pixels <= MAX_RENDER_PIXELS + 1);
    const scale = rayScale(width * ratio, height * ratio, 1);
    assert.ok(pixels * scale ** 2 <= MAX_RAY_PIXELS + 1);
  }
});
test('ray tracing stops at time and sample limits and falls back on a slow frame', () => {
  assert.equal(rayStopReason(3, 100, 5), '');
  assert.match(rayStopReason(64, 100, 5), /complete/);
  assert.match(rayStopReason(3, 20000, 5), /20 second/);
  assert.match(rayStopReason(3, 100, 200), /interactive preview restored/);
});
test('IPC cancellation rejects promptly even while native work is still running', async () => {
  const controller = new AbortController();
  const waiting = withAbort(new Promise(() => {}), controller.signal);
  controller.abort();
  await assert.rejects(waiting, { name: 'AbortError' });
});
test('IPC normal replies and failures are preserved', async () => {
  assert.equal(await withAbort(Promise.resolve(42), new AbortController().signal), 42);
  await assert.rejects(withAbort(Promise.reject(new Error('sidecar stopped')),
    new AbortController().signal), /sidecar stopped/);
});
