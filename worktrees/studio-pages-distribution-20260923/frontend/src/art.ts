import * as THREE from "three";
import { RectAreaLightUniformsLib } from "three/examples/jsm/lights/RectAreaLightUniformsLib.js";
import { GenerateMeshBVHWorker } from "three-mesh-bvh/src/workers/index.js";
import { pixelRatio, rayScale, rayStopReason, MAX_ART_ATOMS, MAX_RAY_SAMPLES } from "./render-policy";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { surfaceNet } from "three/examples/jsm/libs/surfaceNet.js";
import { PhysicalCamera, WebGLPathTracer } from "three-gpu-pathtracer";

type Atom = [string, number, number, number];
type OrbitalStyle = {
  positive?: number;
  negative?: number;
  alpha?: number;
  sides?: string;
};
type ArtScene = {
  atoms: Atom[];
  cube?: string;
  iso?: number;
  sides?: string;
  orbital?: OrbitalStyle;
  label?: string;
};
type ArtSource = { value: string; label: string; disabled?: boolean; reason?: string };
type ArtSettings = Record<string, string>;
type CubeGrid = {
  shape: [number, number, number];
  origin: THREE.Vector3;
  axes: [THREE.Vector3, THREE.Vector3, THREE.Vector3];
  values: Float32Array;
};

const BOHR_TO_ANGSTROM = 0.529177210903;
const MAX_SURFACE_CELLS = 60_000;
const ELEMENTS: Record<string, { color: number; covalent: number; display: number }> = {
  H: { color: 0xf2f2f2, covalent: 0.31, display: 0.34 },
  C: { color: 0x343a42, covalent: 0.76, display: 0.48 },
  N: { color: 0x315ed1, covalent: 0.71, display: 0.46 },
  O: { color: 0xe63226, covalent: 0.66, display: 0.45 },
  F: { color: 0x62c96b, covalent: 0.57, display: 0.43 },
  P: { color: 0xe28c28, covalent: 1.07, display: 0.56 },
  S: { color: 0xe5c83c, covalent: 1.05, display: 0.55 },
  Cl: { color: 0x44ad51, covalent: 1.02, display: 0.56 },
  Br: { color: 0x8d3328, covalent: 1.20, display: 0.61 },
  I: { color: 0x7650a8, covalent: 1.39, display: 0.66 },
};

const stage = document.getElementById("stage")!;
const progress = document.getElementById("progress")!;
const empty = document.getElementById("empty")!;
const contentSelect = document.getElementById("content") as HTMLSelectElement;
const surfaceMaterialSelect = document.getElementById("surfaceMaterial") as HTMLSelectElement;
const surfaceOpacity = document.getElementById("surfaceOpacity") as HTMLInputElement;
const surfacePhases = document.getElementById("surfacePhases") as HTMLSelectElement;
const backgroundSelect = document.getElementById("background") as HTMLSelectElement;
const artSettingIds = [
  "content", "quality", "exposure", "background", "dof",
  "surfaceMaterial", "surfaceOpacity", "surfacePhases",
] as const;

function currentArtSettings(): ArtSettings {
  return Object.fromEntries(artSettingIds.map((id) => [
    id, (document.getElementById(id) as HTMLInputElement | HTMLSelectElement).value,
  ]));
}

function reportArtSettings(): void {
  window.parent.postMessage(
    { type: "oqp-art-settings", settings: currentArtSettings() }, window.location.origin,
  );
}
const renderer = new THREE.WebGLRenderer({
  antialias: true, powerPreference: "default", preserveDrawingBuffer: true,
});

// Avoid the previous synchronous WebKit shader-compile workaround: a driver
// compile cannot be interrupted by a JavaScript timer, even with a Cancel UI.
const webkit = /AppleWebKit/.test(navigator.userAgent) &&
  !/(Chrome|Chromium|Edg)\//.test(navigator.userAgent);
const gl = renderer.getContext();
const debugRenderer = gl.getExtension("WEBGL_debug_renderer_info");
const gpuName = debugRenderer ? String(gl.getParameter(debugRenderer.UNMASKED_RENDERER_WEBGL)) : "";
const canTrace = Boolean(gpuName) && !webkit && (navigator.hardwareConcurrency || 2) > 2 &&
  !/(swiftshader|llvmpipe|software)/i.test(gpuName) &&
  renderer.extensions.has("KHR_parallel_shader_compile");
renderer.setPixelRatio(pixelRatio(stage.clientWidth, stage.clientHeight, devicePixelRatio));
renderer.setSize(Math.max(1, stage.clientWidth), Math.max(1, stage.clientHeight));
RectAreaLightUniformsLib.init();
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.05;
stage.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x111318);
const camera = new PhysicalCamera(34, stage.clientWidth / stage.clientHeight, 0.05, 500);
camera.position.set(5, 3.5, 7);
camera.fStop = 14;
camera.focusDistance = 8;
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

let tracer: WebGLPathTracer | null = null;
let rayGeneration = 0;
let buildingRay = false;
let cancelRayBuild: (() => void) | null = null;
let pendingArtwork: ArtScene | null = null;
let active = window.parent === window;
let dirty = true;
let notice = "Interactive preview · ray tracing is optional";
let rayStarted = 0;
let lastFrame = 0;
const mode = document.getElementById("renderMode") as HTMLSelectElement;
if (!canTrace) {
  mode.querySelector<HTMLOptionElement>('option[value="ray"]')!.disabled = true;
  notice = "Interactive preview · safe mode for this graphics device";
  mode.title = "Ray tracing is disabled on WebKit, software GPUs, or devices without asynchronous shader compilation.";
}

function preview(message = "Interactive preview"): void {
  rayGeneration += 1;
  cancelRayBuild?.();
  tracer?.dispose();
  tracer = null;
  mode.value = "preview";
  notice = message;
  paused = false;
  document.getElementById("pause")!.textContent = "Pause";
  dirty = true;
}

async function startRayTracing(): Promise<void> {
  if (!canTrace || buildingRay || loadingArtwork || !lastScene || sceneError || !active) {
    mode.value = "preview";
    return;
  }
  const generation = ++rayGeneration;
  buildingRay = true;
  notice = "Preparing ray tracing… choose Interactive to cancel";
  let candidate: WebGLPathTracer | null = null;
  let worker: GenerateMeshBVHWorker | null = null;
  const timeout = window.setTimeout(() => {
    if (generation === rayGeneration) preview("Ray tracing preparation timed out; interactive preview restored.");
  }, 15_000);
  try {
    // Paint the status and keep the preview available before preparing a BVH.
    await new Promise<void>((resolve) => window.setTimeout(resolve, 50));
    if (generation !== rayGeneration) return;
    candidate = new WebGLPathTracer(renderer);
    worker = new GenerateMeshBVHWorker();
    candidate.setBVHWorker(worker);
    const [scale, bounces] = (document.getElementById("quality") as HTMLSelectElement).value.split(",").map(Number);
    candidate.renderScale = rayScale(renderer.domElement.width, renderer.domElement.height, scale);
    candidate.bounces = Math.min(4, bounces);
    candidate.transmissiveBounces = 2;
    candidate.tiles.set(8, 8);
    candidate.dynamicLowRes = false;
    candidate.renderDelay = 200;
    candidate.textureSize.set(512, 512);
    const cancelled = new Promise<never>((_resolve, reject) => {
      cancelRayBuild = () => reject(new Error("Ray tracing cancelled"));
    });
    await Promise.race([candidate.setSceneAsync(scene, camera), cancelled]);
    if (generation !== rayGeneration) return;
    tracer = candidate;
    candidate = null;
    rayStarted = performance.now();
    paused = false;
    notice = "";
  } catch (error) {
    if (generation === rayGeneration) preview(`Ray tracing unavailable: ${String(error)}`);
  } finally {
    window.clearTimeout(timeout);
    candidate?.dispose();
    // Upstream's worker type omits WorkerBase.dispose(), present at runtime.
    (worker as (GenerateMeshBVHWorker & { dispose(): void }) | null)?.dispose();
    buildingRay = false;
    cancelRayBuild = null;
  }
}
mode.addEventListener("change", () => {
  if (mode.value === "ray") void startRayTracing();
  else preview();
});

let artwork: THREE.Group | null = null;
let surfaceMeshes: THREE.Mesh[] = [];
let paused = false;
let lastScene: ArtScene | null = null;
let renderRequest = 0;
let sceneError = "";
let loadingArtwork = false;

function setArtworkLoading(loading: boolean): void {
  loadingArtwork = loading;
  mode.disabled = loading;
  (document.getElementById("save") as HTMLButtonElement).disabled = loading;
}
let receivedParentScene = false;

function clearRenderedFrame(): void {
  renderer.setRenderTarget(null);
  renderer.setClearColor(scene.background instanceof THREE.Color ? scene.background : 0x111318, 1);
  renderer.clear(true, true, true);
}

function element(symbol: string) {
  return ELEMENTS[symbol] ?? { color: 0xb7bec9, covalent: 0.77, display: 0.49 };
}

function cylinderBetween(a: THREE.Vector3, b: THREE.Vector3, material: THREE.Material): THREE.Mesh {
  const delta = new THREE.Vector3().subVectors(b, a);
  const mesh = new THREE.Mesh(
    new THREE.CylinderGeometry(0.105, 0.105, delta.length(), 12), material,
  );
  mesh.position.copy(a).add(b).multiplyScalar(0.5);
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), delta.normalize());
  return mesh;
}

function disposeObject(object: THREE.Object3D): void {
  object.traverse((child) => {
    if (!(child instanceof THREE.Mesh)) return;
    child.geometry.dispose();
    const materials = Array.isArray(child.material) ? child.material : [child.material];
    materials.forEach((material) => material.dispose());
  });
}

function finiteNumbers(line: string): number[] {
  const trimmed = line.trim();
  if (!trimmed) return [];
  const values = trimmed.split(/\s+/).map((token) =>
    Number(token.replace(/[dD]/g, "E")));
  if (values.some((value) => !Number.isFinite(value))) throw new Error("cube has invalid numbers");
  return values;
}

function parseCube(text: string): CubeGrid {
  const lines = text.split(/\r?\n/);
  if (lines.length < 6) throw new Error("cube header is incomplete");
  const originRecord = finiteNumbers(lines[2]);
  const axisRecords = [3, 4, 5].map((index) => finiteNumbers(lines[index]));
  if (originRecord.length < 4 || axisRecords.some((record) => record.length !== 4)) {
    throw new Error("cube header is invalid");
  }
  const atomCount = Math.trunc(originRecord[0]);
  const counts = axisRecords.map((record) => Math.trunc(record[0]));
  if (!Number.isInteger(originRecord[0]) || counts.some((count, index) =>
    !Number.isInteger(axisRecords[index][0]) || count === 0)) {
    throw new Error("cube dimensions are invalid");
  }
  const signs = new Set(counts.map((count) => count < 0));
  if (signs.size !== 1) throw new Error("cube axes mix coordinate units");
  const factor = counts[0] < 0 ? 1 : BOHR_TO_ANGSTROM;
  const shape = counts.map(Math.abs) as [number, number, number];
  const pointCount = shape[0] * shape[1] * shape[2];
  if (!Number.isSafeInteger(pointCount) || pointCount < 8 || pointCount > 2_000_000) {
    throw new Error("cube grid is outside the supported size");
  }
  const origin = new THREE.Vector3(...originRecord.slice(1, 4)).multiplyScalar(factor);
  const axes = axisRecords.map((record) =>
    new THREE.Vector3(...record.slice(1, 4)).multiplyScalar(factor)) as CubeGrid["axes"];
  let cursor = 6 + Math.abs(atomCount);
  if (cursor > lines.length) throw new Error("cube atom header is incomplete");
  let datasets = Math.max(1, Math.trunc(originRecord[4] || 1));
  if (atomCount < 0) {
    const identifiers: number[] = [];
    let identifierCount: number | null = null;
    while (cursor < lines.length &&
           (identifierCount === null || identifiers.length < identifierCount)) {
      for (const value of finiteNumbers(lines[cursor++])) {
        if (identifierCount === null) identifierCount = Math.trunc(value);
        else identifiers.push(Math.trunc(value));
      }
    }
    if (!identifierCount || identifiers.length < identifierCount) {
      throw new Error("cube dataset identifiers are incomplete");
    }
    datasets = identifierCount;
  }
  const totalValueCount = pointCount * datasets;
  if (!Number.isSafeInteger(totalValueCount) || totalValueCount > 2_000_000) {
    throw new Error("cube datasets exceed the supported size");
  }
  const raw: number[] = [];
  for (; cursor < lines.length; cursor += 1) {
    if (!lines[cursor].trim()) continue;
    raw.push(...finiteNumbers(lines[cursor]));
    if (raw.length > totalValueCount) break;
  }
  if (raw.length !== totalValueCount) throw new Error("cube grid value count is invalid");
  const values = new Float32Array(pointCount);
  for (let index = 0; index < pointCount; index += 1) values[index] = raw[index * datasets];
  return { shape, origin, axes, values };
}

function sampledIndices(size: number, stride: number): number[] {
  const result: number[] = [];
  for (let index = 0; index < size; index += stride) result.push(index);
  if (result[result.length - 1] !== size - 1) result.push(size - 1);
  return result;
}

function surfaceGeometry(
  grid: CubeGrid, iso: number, phase: 1 | -1, center: THREE.Vector3,
): THREE.BufferGeometry | null {
  const cells = (grid.shape[0] - 1) * (grid.shape[1] - 1) * (grid.shape[2] - 1);
  const stride = Math.max(1, Math.ceil(Math.cbrt(cells / MAX_SURFACE_CELLS)));
  const indices = grid.shape.map((size) => sampledIndices(size, stride)) as
    [number[], number[], number[]];
  const dims = indices.map((axis) => axis.length) as [number, number, number];
  const sampled = new Float32Array(dims[0] * dims[1] * dims[2]);
  let offset = 0;
  for (const x of indices[0]) for (const y of indices[1]) for (const z of indices[2]) {
    sampled[offset++] = grid.values[(x * grid.shape[1] + y) * grid.shape[2] + z];
  }
  const valueAt = (x: number, y: number, z: number) => {
    const i = Math.max(0, Math.min(dims[0] - 1, Math.round(x)));
    const j = Math.max(0, Math.min(dims[1] - 1, Math.round(y)));
    const k = Math.max(0, Math.min(dims[2] - 1, Math.round(z)));
    return phase * sampled[(i * dims[1] + j) * dims[2] + k] - iso;
  };
  const net = surfaceNet(dims, valueAt, [[0, 0, 0], dims]);
  if (!net.positions.length || !net.cells.length) return null;

  const coordinate = (value: number, axis: number) => {
    const source = indices[axis];
    const low = Math.max(0, Math.min(source.length - 1, Math.floor(value)));
    const high = Math.min(source.length - 1, low + 1);
    return THREE.MathUtils.lerp(source[low], source[high], value - low);
  };
  const positions = new Float32Array(net.positions.length * 3);
  net.positions.forEach((position, index) => {
    const x = coordinate(position[0], 0);
    const y = coordinate(position[1], 1);
    const z = coordinate(position[2], 2);
    const point = grid.origin.clone()
      .addScaledVector(grid.axes[0], x)
      .addScaledVector(grid.axes[1], y)
      .addScaledVector(grid.axes[2], z)
      .sub(center);
    positions.set(point.toArray(), index * 3);
  });
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setIndex(net.cells.flat());
  geometry.computeVertexNormals();
  return geometry;
}

function orbitalStyle(): Required<OrbitalStyle> {
  const source = lastScene?.orbital ?? {};
  return {
    positive: Number(source.positive ?? 0x4f8fdd),
    negative: Number(source.negative ?? 0xdd6a4f),
    alpha: Number(source.alpha ?? 0.85),
    sides: String(lastScene?.sides ?? source.sides ?? "both"),
  };
}

function makeSurfaceMaterial(color: number): THREE.MeshPhysicalMaterial {
  const mode = surfaceMaterialSelect.value;
  const alpha = THREE.MathUtils.clamp(orbitalStyle().alpha * +surfaceOpacity.value, 0.05, 1);
  return new THREE.MeshPhysicalMaterial({
    color,
    roughness: mode === "matte" ? 0.5 : mode === "gloss" ? 0.14 : 0.08,
    metalness: mode === "gloss" ? 0.08 : 0,
    transmission: mode === "glass" ? 0.38 : 0,
    thickness: mode === "glass" ? 0.45 : 0,
    ior: 1.35,
    opacity: alpha,
    transparent: alpha < 0.999,
    side: THREE.DoubleSide,
  });
}

function phaseVisible(phase: "positive" | "negative"): boolean {
  const selected = surfacePhases.value === "inherit" ? orbitalStyle().sides : surfacePhases.value;
  return selected === "both" || selected === phase;
}

function updateSurfaceAppearance(_rebuild = true): void {
  if (loadingArtwork) return;
  preview("Surface changed · interactive preview");
  const style = orbitalStyle();
  for (const mesh of surfaceMeshes) {
    const phase = mesh.userData.phase as "positive" | "negative";
    mesh.visible = phaseVisible(phase);
    const materials = Array.isArray(mesh.material) ? mesh.material : [mesh.material];
    materials.forEach((material) => material.dispose());
    mesh.material = makeSurfaceMaterial(phase === "positive" ? style.positive : style.negative);
  }
  if (!surfaceMeshes.length) return;
  if (!surfaceMeshes.some((mesh) => mesh.visible)) {
    const iso = Math.max(1e-8, Math.abs(lastScene?.iso ?? 0.05));
    sceneError = `no requested-phase isosurface crosses ${iso}`;
    clearRenderedFrame();
    return;
  }
  sceneError = "";
  dirty = true;
}

function addMolecule(group: THREE.Group, atoms: Atom[], center: THREE.Vector3): void {
  const points = atoms.map(([, x, y, z]) => new THREE.Vector3(x, y, z).sub(center));
  atoms.forEach(([symbol], index) => {
    const spec = element(symbol);
    const material = new THREE.MeshStandardMaterial({
      color: spec.color, roughness: 0.26, metalness: 0.03,
    });
    const atom = new THREE.Mesh(new THREE.SphereGeometry(spec.display, 24, 16), material);
    atom.position.copy(points[index]);
    group.add(atom);
  });
  const bondMaterial = new THREE.MeshStandardMaterial({
    color: 0xaeb5bf, roughness: 0.38, metalness: 0.02,
  });
  for (let i = 0; i < atoms.length; i += 1) {
    for (let j = i + 1; j < atoms.length; j += 1) {
      const cutoff = 1.18 * (element(atoms[i][0]).covalent + element(atoms[j][0]).covalent);
      const distance = points[i].distanceTo(points[j]);
      if (distance > 0.35 && distance <= cutoff) {
        group.add(cylinderBetween(points[i], points[j], bondMaterial.clone()));
      }
    }
  }
  bondMaterial.dispose();
}

async function setArtwork(data: ArtScene): Promise<void> {
  if (!data.atoms.length && !data.cube) return;
  const request = ++renderRequest;
  preview();
  setArtworkLoading(true);
  clearRenderedFrame();
  try {
    await buildArtwork(data, request);
  } finally {
    // An older download must not enable controls for a newer pending scene.
    if (request === renderRequest) setArtworkLoading(false);
  }
}

async function buildArtwork(data: ArtScene, request: number): Promise<void> {
  if (data.atoms.length > MAX_ART_ATOMS) {
    if (artwork) {
      scene.remove(artwork);
      disposeObject(artwork);
      artwork = null;
    }
    surfaceMeshes = [];
    lastScene = null;
    sceneError = `Art supports up to ${MAX_ART_ATOMS} atoms; use Analysis for larger structures.`;
    progress.textContent = sceneError;
    clearRenderedFrame();
    return;
  }
  lastScene = data;
  sceneError = "";
  progress.textContent = data.cube ? "loading volumetric surface" : "building preview";
  let cube: CubeGrid | null = null;
  if (data.cube) {
    try {
      const response = await fetch(data.cube);
      if (!response.ok) throw new Error(`cube request failed (${response.status})`);
      const text = await response.text();
      if (text.length > 24_000_000) throw new Error("cube is too large for interactive Art; use Analysis");
      if (request !== renderRequest) return;
      cube = parseCube(text);
    } catch (error) {
      if (request !== renderRequest) return;
      sceneError = `surface unavailable: ${error instanceof Error ? error.message : String(error)}`;
      if (artwork) {
        scene.remove(artwork);
        disposeObject(artwork);
        artwork = null;
      }
      surfaceMeshes = [];
      lastScene = null;
      progress.textContent = sceneError;
      clearRenderedFrame();
      return;
    }
  }
  if (request !== renderRequest) return;
  if (!data.atoms.length && !cube) return;
  if (artwork) {
    scene.remove(artwork);
    disposeObject(artwork);
  }
  artwork = new THREE.Group();
  surfaceMeshes = [];
  const center = new THREE.Vector3();
  if (data.atoms.length) {
    data.atoms.forEach(([, x, y, z]) => center.add(new THREE.Vector3(x, y, z)));
    center.multiplyScalar(1 / data.atoms.length);
    addMolecule(artwork, data.atoms, center);
  } else if (cube) {
    center.copy(cube.origin)
      .addScaledVector(cube.axes[0], (cube.shape[0] - 1) / 2)
      .addScaledVector(cube.axes[1], (cube.shape[1] - 1) / 2)
      .addScaledVector(cube.axes[2], (cube.shape[2] - 1) / 2);
  }

  if (cube) {
    progress.textContent = "building positive and negative isosurfaces";
    const iso = Math.max(1e-8, Math.abs(data.iso ?? 0.05));
    for (const [phase, sign] of [["positive", 1], ["negative", -1]] as const) {
      const geometry = surfaceGeometry(cube, iso, sign, center);
      if (!geometry) continue;
      const style = orbitalStyle();
      const mesh = new THREE.Mesh(
        geometry, makeSurfaceMaterial(phase === "positive" ? style.positive : style.negative),
      );
      mesh.userData.phase = phase;
      mesh.visible = phaseVisible(phase);
      artwork.add(mesh);
      surfaceMeshes.push(mesh);
    }
  }
  if (cube && !surfaceMeshes.length) {
    sceneError = `no isosurface crosses ${Math.max(1e-8, Math.abs(data.iso ?? 0.05))}`;
    disposeObject(artwork);
    artwork = null;
    surfaceMeshes = [];
    lastScene = null;
    progress.textContent = sceneError;
    clearRenderedFrame();
    return;
  }
  if (cube && !surfaceMeshes.some((mesh) => mesh.visible)) {
    sceneError = `no requested-phase isosurface crosses ${Math.max(1e-8, Math.abs(data.iso ?? 0.05))}`;
    progress.textContent = sceneError;
    clearRenderedFrame();
  }
  scene.add(artwork);
  const box = new THREE.Box3().setFromObject(artwork);
  const sphere = box.getBoundingSphere(new THREE.Sphere());
  floor.position.y = box.min.y - 1.25;
  const distance = Math.max(4.5, sphere.radius / Math.tan(THREE.MathUtils.degToRad(camera.fov * 0.45)));
  camera.position.set(distance * 0.65, distance * 0.42, distance);
  camera.near = Math.max(0.02, distance / 100);
  camera.far = distance * 20;
  camera.focusDistance = camera.position.length();
  camera.updateProjectionMatrix();
  controls.target.set(0, 0, 0);
  controls.saveState();
  controls.update();
  empty.style.display = "none";
  progress.textContent = sceneError || (surfaceMeshes.length
    ? `building ${data.label ?? "volumetric"} preview`
    : "building molecular preview");
  dirty = true;
  notice = canTrace ? "Interactive preview · select Ray tracing when needed"
    : "Interactive preview · safe mode for this graphics device";
}

const floor = new THREE.Mesh(
  new THREE.PlaneGeometry(80, 80),
  new THREE.MeshStandardMaterial({ color: 0x242830, roughness: 0.68, metalness: 0 }),
);
floor.rotation.x = -Math.PI / 2;
floor.position.y = -2.2;
scene.add(floor);

const BACKGROUNDS: Record<string, { scene: number; floor: number }> = {
  studio: { scene: 0x111318, floor: 0x242830 },
  neutral: { scene: 0x30343a, floor: 0x565d66 },
  white: { scene: 0xf4f6f8, floor: 0xd8dde3 },
  black: { scene: 0x000000, floor: 0x101216 },
};

function setBackground(): void {
  const preset = BACKGROUNDS[backgroundSelect.value] ?? BACKGROUNDS.studio;
  scene.background = new THREE.Color(preset.scene);
  (floor.material as THREE.MeshStandardMaterial).color.setHex(preset.floor);
  preview("Background changed · interactive preview");
}
const key = new THREE.RectAreaLight(0xffffff, 28, 5, 5);
key.position.set(4, 6, 5);
key.lookAt(0, 0, 0);
scene.add(key);
const fill = new THREE.RectAreaLight(0x8cbfff, 18, 4, 4);
fill.position.set(-5, 2, 3);
fill.lookAt(0, 0, 0);
scene.add(fill);
const rim = new THREE.RectAreaLight(0xffd5aa, 16, 3, 3);
rim.position.set(1, 3, -5);
rim.lookAt(0, 0, 0);
scene.add(rim);

controls.addEventListener("change", () => {
  camera.focusDistance = camera.position.distanceTo(controls.target);
  tracer?.updateCamera();
  tracer?.reset();
  dirty = true;
});

document.getElementById("quality")!.addEventListener("change", () => {
  preview("Quality selected · start Ray tracing to apply");
});
document.getElementById("exposure")!.addEventListener("input", (event) => {
  renderer.toneMappingExposure = +(event.target as HTMLInputElement).value;
  tracer?.reset();
  dirty = true;
});
backgroundSelect.addEventListener("change", setBackground);
contentSelect.addEventListener("change", () => {
  progress.textContent = `loading ${contentSelect.selectedOptions[0]?.textContent ?? "content"}`;
  window.parent.postMessage(
    { type: "oqp-art-source-request", value: contentSelect.value }, window.location.origin,
  );
});
document.getElementById("dof")!.addEventListener("change", (event) => {
  const strength = +(event.target as HTMLSelectElement).value;
  camera.fStop = strength === 0 ? 1000 : strength === 1 ? 14 : 5.6;
  camera.focusDistance = camera.position.distanceTo(controls.target);
  tracer?.updateCamera();
  tracer?.reset();
  dirty = true;
});
surfaceMaterialSelect.addEventListener("change", () => updateSurfaceAppearance(false));
surfaceOpacity.addEventListener("input", () => updateSurfaceAppearance(false));
surfacePhases.addEventListener("change", () => updateSurfaceAppearance(true));
document.getElementById("pause")!.addEventListener("click", (event) => {
  if (buildingRay) { preview("Ray tracing cancelled · interactive preview"); return; }
  paused = !paused;
  if (!paused) {
    rayStarted = performance.now();
    if (tracer && tracer.samples >= MAX_RAY_SAMPLES) tracer.reset();
  }
  dirty = true;
  (event.target as HTMLButtonElement).textContent = paused ? "Resume" : "Pause";
});
document.getElementById("reset")!.addEventListener("click", () => {
  controls.reset();
  tracer?.updateCamera();
  tracer?.reset();
  dirty = true;
});
document.getElementById("save")!.addEventListener("click", () => {
  if (loadingArtwork || sceneError || !artwork) return;
  const link = document.createElement("a");
  link.download = "oqp-art.png";
  link.href = renderer.domElement.toDataURL("image/png");
  link.click();
});

window.addEventListener("message", (event) => {
  if (event.origin !== window.location.origin || event.source !== window.parent) return;
  if (event.data?.type === "oqp-art-active") {
    active = Boolean(event.data.active);
    if (!active) preview("Rendering stopped while Art is hidden");
    dirty = true;
    if (active) {
      resize();
      if (pendingArtwork) { const next = pendingArtwork; pendingArtwork = null; void setArtwork(next); }
    }
    return;
  }
  if (event.data?.type === "oqp-art-scene") {
    receivedParentScene = true;
    if (active) void setArtwork(event.data.scene as ArtScene);
    else pendingArtwork = event.data.scene as ArtScene;
  } else if (event.data?.type === "oqp-art-clear") {
    receivedParentScene = true;
    renderRequest += 1;
    setArtworkLoading(false);
    pendingArtwork = null;
    preview();
    if (artwork) {
      scene.remove(artwork);
      disposeObject(artwork);
      artwork = null;
    }
    surfaceMeshes = [];
    lastScene = null;
    sceneError = "";
    empty.style.display = "";
    progress.textContent = "waiting for a structure";
    clearRenderedFrame();
  } else if (event.data?.type === "oqp-art-sources") {
    const sources = event.data.sources as ArtSource[];
    contentSelect.replaceChildren(...sources.map((source) => {
      const option = document.createElement("option");
      option.value = source.value;
      option.textContent = source.label;
      option.disabled = Boolean(source.disabled);
      if (source.reason) option.title = source.reason;
      return option;
    }));
    contentSelect.value = String(event.data.selected ?? "molecule");
  } else if (event.data?.type === "oqp-art-style") {
    if (!lastScene) return;
    lastScene = { ...lastScene, orbital: event.data.orbital as OrbitalStyle };
    updateSurfaceAppearance(false);
  } else if (event.data?.type === "oqp-art-settings-apply") {
    const settings = event.data.settings as ArtSettings;
    for (const id of artSettingIds) {
      const control = document.getElementById(id) as HTMLInputElement | HTMLSelectElement;
      if (settings[id] === undefined) continue;
      const value = String(settings[id]);
      if (control instanceof HTMLSelectElement &&
          ![...control.options].some((option) => option.value === value)) continue;
      if (control.value === value) continue;
      control.value = value;
      control.dispatchEvent(new Event(control instanceof HTMLInputElement ? "input" : "change"));
    }
    reportArtSettings();
  }
});
window.setTimeout(() => {
  if (!receivedParentScene) {
    void setArtwork({ atoms: [
      ["O", 0.001, 0.398, 0],
      ["H", -0.764, -0.197, 0],
      ["H", 0.763, -0.201, 0],
    ] });
  }
}, 350);
function resize(): void {
  const width = Math.max(1, stage.clientWidth);
  const height = Math.max(1, stage.clientHeight);
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
  renderer.setPixelRatio(pixelRatio(width, height, devicePixelRatio));
  renderer.setSize(width, height);
  preview("View resized · interactive preview");
}
window.addEventListener("resize", resize);
document.addEventListener("visibilitychange", () => {
  if (document.hidden) preview("Rendering stopped while the window is hidden");
  dirty = true;
});
renderer.domElement.addEventListener("webglcontextlost", (event) => {
  event.preventDefault();
  preview("Graphics context lost; reopen Art to recover");
  sceneError = notice;
});

function animate(now: number): void {
  requestAnimationFrame(animate);
  if (!active || document.hidden || now - lastFrame < 50) return;
  lastFrame = now;
  controls.update();
  if (lastScene && !sceneError && !loadingArtwork && !paused) {
    try {
      if (tracer) {
        const started = performance.now();
        tracer.renderSample();
        const reason = rayStopReason(tracer.samples, now - rayStarted, performance.now() - started);
        if (reason.startsWith("Ray tracing was too slow")) preview(reason);
        else if (reason) {
          paused = true;
          notice = reason;
          document.getElementById("pause")!.textContent = "Resume";
        }
      } else if (dirty) {
        renderer.render(scene, camera);
        dirty = false;
      }
    } catch (error) {
      preview(`Interactive mode restored: ${String(error)}`);
    }
  }
  progress.classList.toggle("busy", loadingArtwork || buildingRay || Boolean(tracer && !paused));
  progress.textContent = loadingArtwork ? "Loading artwork…" : sceneError || (buildingRay ? notice : tracer
    ? `${notice || (paused ? "Paused" : "Tracing")} · ${tracer.samples.toFixed(0)} / 64 samples`
    : notice);
}
requestAnimationFrame(animate);

document.getElementById("toolbar")!.addEventListener("input", reportArtSettings);
document.getElementById("toolbar")!.addEventListener("change", reportArtSettings);
window.parent.postMessage(
  { type: "oqp-art-ready", settings: currentArtSettings() }, window.location.origin,
);
