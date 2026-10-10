// Conservative budgets for interactive desktop use, independent of display DPI.
export const MAX_RENDER_PIXELS = 1_000_000;
export const MAX_RAY_PIXELS = 260_000;
export const MAX_RAY_SAMPLES = 64;
export const MAX_RAY_MS = 20_000;
export const MAX_ART_ATOMS = 300;

export function pixelRatio(width: number, height: number, deviceRatio: number): number {
  return Math.min(Math.max(0.1, deviceRatio || 1), 1.25,
    Math.sqrt(MAX_RENDER_PIXELS / Math.max(1, width * height)));
}

export function rayScale(width: number, height: number, requested: number): number {
  return Math.min(Math.max(0.1, requested || 0.55), 0.75,
    Math.sqrt(MAX_RAY_PIXELS / Math.max(1, width * height)));
}

export function rayStopReason(samples: number, elapsedMs: number, frameMs: number): string {
  if (frameMs > 180) return "Ray tracing was too slow; interactive preview restored.";
  if (samples >= MAX_RAY_SAMPLES) return "Ray tracing complete (64 samples).";
  if (elapsedMs >= MAX_RAY_MS) return "Ray tracing paused at the 20 second limit.";
  return "";
}
