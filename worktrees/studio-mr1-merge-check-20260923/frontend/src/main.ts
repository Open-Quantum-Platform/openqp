// OQP Studio frontend — Builder → Method → Run → Results workflow.
// Vanilla TS for now; Phase 1/2 swap the preview canvas for Mol* and add
// the Ketcher sketcher and keyword-schema forms.

import { installApiFetch } from "./api";

installApiFetch();

type Atom = [string, number, number, number];

// Geometries generated with RDKit ETKDGv3 + MMFF94 (seed 0xF00D).
const SAMPLES: Record<string, { label: string; atoms: Atom[] }> = {
  water: { label: "Water (H₂O)", atoms: [["O",0.001,0.398,0.000],["H",-0.764,-0.197,0.000],["H",0.763,-0.201,0.000]] },
  ammonia: { label: "Ammonia (NH₃)", atoms: [["N",0.021,-0.007,0.295],["H",0.916,-0.162,-0.168],["H",-0.612,-0.716,-0.073],["H",-0.325,0.886,-0.054]] },
  methane: { label: "Methane (CH₄)", atoms: [["C",-0.000,0.000,0.000],["H",-0.635,-0.826,-0.329],["H",-0.442,0.947,-0.316],["H",0.085,-0.017,1.089],["H",0.992,-0.105,-0.444]] },
  methanol: { label: "Methanol (CH₃OH)", atoms: [["C",-0.371,-0.013,0.018],["O",0.950,-0.078,0.525],["H",-1.074,-0.189,0.835],["H",-0.550,0.978,-0.407],["H",-0.506,-0.780,-0.748],["H",1.551,0.082,-0.222]] },
  ethanol: { label: "Ethanol (C₂H₅OH)", atoms: [["C",0.872,0.199,0.135],["C",-0.461,-0.516,0.045],["O",-1.326,0.182,-0.839],["H",0.739,1.219,0.511],["H",1.556,-0.335,0.800],["H",1.332,0.283,-0.855],["H",-0.939,-0.580,1.027],["H",-0.325,-1.530,-0.341],["H",-1.448,1.078,-0.483]] },
  formaldehyde: { label: "Formaldehyde (CH₂O)", atoms: [["C",-0.012,0.002,-0.001],["O",1.197,-0.164,0.098],["H",-0.718,-0.842,-0.050],["H",-0.467,1.004,-0.047]] },
  acetone: { label: "Acetone (C₃H₆O)", atoms: [["C",-1.269,-0.226,-0.089],["C",0.006,0.198,0.589],["O",0.017,0.591,1.754],["C",1.266,0.119,-0.230],["H",-1.190,-1.272,-0.395],["H",-1.450,0.410,-0.959],["H",-2.108,-0.123,0.605],["H",1.436,-0.915,-0.541],["H",1.176,0.767,-1.105],["H",2.117,0.452,0.371]] },
  aceticacid: { label: "Acetic acid (CH₃COOH)", atoms: [["C",-0.956,-0.066,0.101],["C",0.477,0.276,-0.141],["O",0.909,1.333,-0.567],["O",1.306,-0.738,0.169],["H",-1.109,-0.293,1.159],["H",-1.244,-0.916,-0.521],["H",-1.582,0.791,-0.166],["H",2.198,-0.386,-0.034]] },
  ethylene: { label: "Ethylene (C₂H₄)", atoms: [["C",0.570,0.134,-0.322],["C",-0.570,-0.134,0.322],["H",0.813,-0.359,-1.258],["H",1.281,0.852,0.075],["H",-1.281,-0.852,-0.075],["H",-0.813,0.359,1.258]] },
  acetylene: { label: "Acetylene (C₂H₂)", atoms: [["C",0.595,-0.007,-0.106],["C",-0.595,0.007,0.044],["H",1.653,-0.020,-0.239],["H",-1.653,0.020,0.177]] },
  butadiene: { label: "1,3-Butadiene (C₄H₆)", atoms: [["C",-1.607,0.846,-0.458],["C",-0.354,0.609,-0.053],["C",0.354,-0.609,-0.364],["C",1.607,-0.846,0.041],["H",-2.169,0.131,-1.049],["H",-2.101,1.778,-0.200],["H",0.163,1.362,0.539],["H",-0.163,-1.362,-0.955],["H",2.169,-0.131,0.633],["H",2.101,-1.778,-0.217]] },
  benzene: { label: "Benzene (C₆H₆)", atoms: [["C",-0.303,-1.361,-0.009],["C",-1.331,-0.418,0.015],["C",-1.027,0.943,0.024],["C",0.303,1.361,0.009],["C",1.331,0.418,-0.015],["C",1.027,-0.943,-0.024],["H",-0.539,-2.422,-0.016],["H",-2.367,-0.744,0.027],["H",-1.828,1.678,0.042],["H",0.539,2.422,0.016],["H",2.367,0.744,-0.027],["H",1.828,-1.678,-0.042]] },
  toluene: { label: "Toluene (C₇H₈)", atoms: [["C",2.221,0.034,-0.027],["C",0.721,0.010,-0.034],["C",0.031,-1.206,0.028],["C",-1.365,-1.227,0.054],["C",-2.081,-0.032,0.029],["C",-1.403,1.185,-0.018],["C",-0.007,1.206,-0.044],["H",2.590,0.071,1.002],["H",2.628,-0.856,-0.518],["H",2.600,0.906,-0.571],["H",0.577,-2.147,0.056],["H",-1.893,-2.176,0.097],["H",-3.168,-0.048,0.050],["H",-1.960,2.117,-0.032],["H",0.509,2.163,-0.073]] },
  phenol: { label: "Phenol (C₆H₅OH)", atoms: [["O",-2.475,0.279,0.115],["C",-1.120,0.140,0.058],["C",-0.341,1.202,0.503],["C",1.050,1.102,0.462],["C",1.651,-0.059,-0.023],["C",0.862,-1.121,-0.468],["C",-0.530,-1.024,-0.429],["H",-2.882,-0.534,-0.226],["H",-0.815,2.104,0.879],["H",1.663,1.930,0.809],["H",2.735,-0.137,-0.055],["H",1.333,-2.025,-0.846],["H",-1.131,-1.857,-0.778]] },
  aniline: { label: "Aniline (C₆H₅NH₂)", atoms: [["N",2.275,-0.456,0.086],["C",0.912,-0.153,-0.017],["C",0.023,-1.083,-0.562],["C",-1.353,-0.840,-0.564],["C",-1.856,0.320,0.020],["C",-0.984,1.226,0.620],["C",0.391,0.977,0.618],["H",2.867,0.362,0.190],["H",2.608,-1.084,-0.638],["H",0.396,-2.009,-0.992],["H",-2.031,-1.564,-1.007],["H",-2.926,0.508,0.025],["H",-1.373,2.118,1.102],["H",1.054,1.677,1.119]] },
  naphthalene: { label: "Naphthalene (C₁₀H₈)", atoms: [["C",2.427,0.669,0.269],["C",2.429,-0.710,0.077],["C",1.224,-1.391,-0.105],["C",0.001,-0.700,-0.098],["C",-1.220,-1.370,-0.279],["C",-2.427,-0.669,-0.269],["C",-2.429,0.710,-0.077],["C",-1.224,1.391,0.105],["C",-0.001,0.700,0.098],["C",1.220,1.370,0.280],["H",3.364,1.200,0.411],["H",3.368,-1.256,0.069],["H",1.246,-2.468,-0.253],["H",-1.238,-2.447,-0.431],["H",-3.364,-1.200,-0.412],["H",-3.368,1.256,-0.069],["H",-1.246,2.468,0.253],["H",1.238,2.447,0.431]] },
  anthracene: { label: "Anthracene (C₁₄H₁₀)", atoms: [["C",-3.660,0.660,-0.212],["C",-3.477,-0.147,-1.331],["C",-2.206,-0.636,-1.636],["C",-1.101,-0.324,-0.825],["C",0.184,-0.808,-1.120],["C",1.288,-0.495,-0.309],["C",2.573,-0.980,-0.602],["C",3.660,-0.660,0.212],["C",3.477,0.147,1.331],["C",2.206,0.636,1.636],["C",1.101,0.324,0.825],["C",-0.184,0.808,1.120],["C",-1.288,0.495,0.309],["C",-2.573,0.980,0.602],["H",-4.649,1.042,0.027],["H",-4.322,-0.397,-1.967],["H",-2.082,-1.266,-2.514],["H",0.328,-1.440,-1.994],["H",2.737,-1.612,-1.472],["H",4.649,-1.042,-0.027],["H",4.322,0.397,1.967],["H",2.082,1.266,2.514],["H",-0.328,1.440,1.994],["H",-2.737,1.612,1.472]] },
  pyridine: { label: "Pyridine (C₅H₅N)", atoms: [["C",-0.329,1.137,-0.014],["C",1.019,0.795,-0.050],["C",1.355,-0.549,-0.032],["N",0.448,-1.550,0.019],["C",-0.852,-1.187,0.053],["C",-1.286,0.129,0.039],["H",-0.631,2.180,-0.027],["H",1.789,1.557,-0.092],["H",2.393,-0.867,-0.059],["H",-1.561,-2.009,0.094],["H",-2.344,0.363,0.068]] },
  pyrrole: { label: "Pyrrole (C₄H₅N)", atoms: [["C",-1.069,0.511,-0.153],["C",-0.819,-0.869,0.034],["C",0.542,-1.025,0.170],["N",1.122,0.213,0.071],["C",0.147,1.156,-0.126],["H",-2.033,0.984,-0.293],["H",-1.553,-1.664,0.066],["H",1.146,-1.908,0.328],["H",2.114,0.401,0.133],["H",0.402,2.202,-0.230]] },
  furan: { label: "Furan (C₄H₄O)", atoms: [["C",0.567,-0.858,0.017],["C",-0.824,-0.611,-0.079],["C",-0.972,0.756,-0.050],["O",0.238,1.363,0.058],["C",1.170,0.376,0.098],["H",1.070,-1.814,0.026],["H",-1.620,-1.337,-0.160],["H",-1.825,1.419,-0.094],["H",2.196,0.706,0.185]] },
  imidazole: { label: "Imidazole (C₃H₄N₂)", atoms: [["C",-1.014,0.697,-0.096],["C",-0.759,-0.653,-0.071],["N",0.599,-0.755,0.057],["C",1.102,0.514,0.104],["N",0.150,1.413,0.014],["H",-1.973,1.192,-0.186],["H",-1.402,-1.519,-0.132],["H",1.138,-1.608,0.108],["H",2.160,0.718,0.204]] },
  thiophene: { label: "Thiophene (C₄H₄S)", atoms: [["C",0.698,-0.723,0.012],["C",-0.715,-0.686,-0.170],["C",-1.212,0.596,-0.094],["S",0.021,1.748,0.196],["C",1.225,0.532,0.220],["H",1.298,-1.625,-0.009],["H",-1.336,-1.556,-0.348],["H",-2.240,0.916,-0.194],["H",2.260,0.798,0.386]] },
  formamide: { label: "Formamide (CH₃NO)", atoms: [["C",0.644,0.120,-0.119],["O",1.512,-0.722,0.055],["N",-0.674,-0.109,0.119],["H",0.847,1.139,-0.484],["H",-1.384,0.594,-0.028],["H",-0.946,-1.023,0.457]] },
  urea: { label: "Urea (CH₄N₂O)", atoms: [["N",-1.109,-0.321,0.075],["C",-0.022,0.478,0.009],["O",-0.079,1.695,0.033],["N",1.134,-0.213,-0.085],["H",-0.993,-1.167,0.614],["H",-1.959,0.190,0.281],["H",1.932,0.382,-0.270],["H",1.096,-1.045,-0.657]] },
  glycine: { label: "Glycine (C₂H₅NO₂)", atoms: [["N",-1.399,0.591,0.239],["C",-0.631,-0.480,-0.429],["C",0.862,-0.396,-0.133],["O",1.700,-1.230,-0.439],["O",1.234,0.735,0.506],["H",-1.244,0.509,1.245],["H",-0.953,1.480,-0.002],["H",-0.773,-0.399,-1.510],["H",-1.001,-1.448,-0.080],["H",2.205,0.637,0.604]] },
  uracil: { label: "Uracil (C₄H₄N₂O₂)", atoms: [["O",2.582,-0.513,0.003],["C",1.378,-0.280,0.002],["C",0.874,1.108,-0.006],["C",-0.446,1.311,-0.008],["N",-1.316,0.258,-0.002],["C",-0.921,-1.054,0.006],["O",-1.730,-1.975,0.011],["N",0.429,-1.263,0.007],["H",1.602,1.908,-0.011],["H",-0.880,2.305,-0.013],["H",-2.318,0.414,-0.003],["H",0.747,-2.219,0.013]] },
  thymine: { label: "Thymine (C₅H₆N₂O₂)", atoms: [["C",-2.213,-0.190,0.151],["C",-0.725,-0.124,0.057],["C",0.089,-1.184,0.136],["N",1.446,-1.048,0.042],["C",2.077,0.152,-0.138],["O",3.298,0.238,-0.219],["N",1.259,1.241,-0.221],["C",-0.106,1.215,-0.139],["O",-0.783,2.234,-0.222],["H",-2.572,0.404,0.997],["H",-2.565,-1.218,0.294],["H",-2.673,0.194,-0.766],["H",-0.281,-2.194,0.279],["H",2.052,-1.858,0.104],["H",1.698,2.138,-0.354]] },
  cytosine: { label: "Cytosine (C₄H₅N₃O)", atoms: [["N",2.291,-0.133,-0.279],["C",0.932,-0.007,-0.190],["C",0.156,-0.997,0.595],["C",-1.160,-0.795,0.634],["N",-1.703,0.268,-0.021],["C",-0.966,1.174,-0.747],["O",-1.547,2.101,-1.304],["N",0.394,1.000,-0.808],["H",2.720,0.710,-0.651],["H",2.770,-0.490,0.538],["H",0.654,-1.823,1.079],["H",-1.841,-1.449,1.167],["H",-2.699,0.442,-0.012]] },
  adenine: { label: "Adenine (C₅H₅N₅)", atoms: [["N",-2.061,-1.393,-0.699],["C",-1.266,-0.336,-0.275],["N",-1.875,0.818,0.088],["C",-1.102,1.850,0.496],["N",0.241,1.898,0.614],["C",0.785,0.726,0.258],["N",2.110,0.401,0.256],["C",2.216,-0.897,-0.166],["N",1.039,-1.424,-0.436],["C",0.133,-0.412,-0.183],["H",-3.012,-1.314,-0.364],["H",-1.605,-2.295,-0.622],["H",-1.635,2.758,0.764],["H",2.862,1.018,0.523],["H",3.170,-1.400,-0.255]] },
  azobenzene: { label: "Azobenzene (C₁₂H₁₀N₂)", atoms: [["C",-4.447,0.576,-0.037],["C",-3.834,-0.592,-0.489],["C",-2.462,-0.775,-0.305],["C",-1.688,0.210,0.324],["C",-2.319,1.371,0.784],["C",-3.691,1.557,0.601],["N",-0.310,0.041,0.552],["N",0.305,-0.082,-0.530],["C",1.685,-0.239,-0.303],["C",2.198,-1.329,0.413],["C",3.575,-1.476,0.596],["C",4.452,-0.537,0.056],["C",3.955,0.545,-0.669],["C",2.577,0.689,-0.850],["H",-5.515,0.719,-0.180],["H",-4.422,-1.361,-0.984],["H",-1.988,-1.689,-0.653],["H",-1.732,2.136,1.286],["H",-4.168,2.466,0.959],["H",1.517,-2.067,0.830],["H",3.960,-2.323,1.158],["H",5.524,-0.650,0.198],["H",4.637,1.276,-1.094],["H",2.192,1.531,-1.419]] },
  stilbene: { label: "trans-Stilbene (C₁₄H₁₂)", atoms: [["C",-4.650,-0.513,-0.180],["C",-4.028,0.373,0.697],["C",-2.636,0.378,0.817],["C",-1.851,-0.496,0.053],["C",-2.490,-1.393,-0.813],["C",-3.882,-1.398,-0.933],["C",-0.387,-0.544,0.172],["C",0.387,0.544,0.318],["C",1.851,0.496,0.437],["C",2.636,-0.378,-0.327],["C",4.028,-0.373,-0.206],["C",4.650,0.513,0.670],["C",3.882,1.398,1.424],["C",2.490,1.393,1.303],["H",-5.733,-0.518,-0.273],["H",-4.626,1.057,1.293],["H",-2.174,1.064,1.522],["H",-1.907,-2.092,-1.409],["H",-4.365,-2.094,-1.614],["H",0.052,-1.538,0.127],["H",-0.052,1.538,0.363],["H",2.174,-1.064,-1.032],["H",4.626,-1.057,-0.803],["H",5.733,0.518,0.763],["H",4.365,2.094,2.105],["H",1.907,2.092,1.899]] },
  retinalPSB: { label: "PSB3 model (C₅H₈N⁺)", atoms: [["C",-0.579,0.222,0.149],["C",0.631,-0.354,0.187],["C",1.784,0.474,0.566],["N",2.978,0.016,0.631],["C",-1.773,-0.507,-0.211],["C",-2.984,0.061,-0.251],["H",-0.714,1.276,0.391],["H",0.767,-1.404,-0.054],["H",1.623,1.538,0.806],["H",3.763,0.614,0.897],["H",3.202,-0.957,0.422],["H",-1.688,-1.563,-0.461],["H",-3.153,1.108,-0.016],["H",-3.858,-0.524,-0.527]] },
  caffeine: { label: "Caffeine (C₈H₁₀N₄O₂)", atoms: [["C",3.294,0.253,0.440],["N",2.076,-0.507,0.312],["C",1.964,-1.871,0.367],["N",0.714,-2.262,0.215],["C",0.013,-1.105,0.059],["C",0.820,-0.009,0.113],["C",0.338,1.321,-0.022],["O",1.075,2.301,0.029],["N",-1.045,1.378,-0.215],["C",-1.660,2.682,-0.367],["C",-1.914,0.268,-0.277],["O",-3.126,0.421,-0.452],["N",-1.343,-0.997,-0.133],["C",-2.167,-2.191,-0.185],["H",3.200,0.918,1.302],["H",3.445,0.826,-0.478],["H",4.135,-0.429,0.591],["H",2.813,-2.526,0.518],["H",-2.167,2.725,-1.337],["H",-2.409,2.816,0.420],["H",-0.936,3.498,-0.310],["H",-1.826,-2.825,-1.009],["H",-2.070,-2.733,0.761],["H",-3.222,-1.951,-0.342]] },
};

const VIEWABLE = [".log", ".out", ".txt", ".json", ".molden", ".cube", ".cub",
                  ".xyz", ".trj", ".inp", ".oqp"];

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

function escapeHtml(value: unknown): string {
  return String(value).replace(/[&<>"']/g, (character) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  })[character]!);
}

function busyStatus(element: HTMLElement, label: string): () => void {
  const started = performance.now();
  element.classList.add("busy");
  element.setAttribute("aria-busy", "true");
  const update = () => { element.textContent = `${label} · ${Math.floor((performance.now() - started) / 1000)}s`; };
  update();
  const timer = window.setInterval(update, 1000);
  return () => {
    window.clearInterval(timer);
    element.classList.remove("busy");
    element.setAttribute("aria-busy", "false");
  };
}

// ---------- tabs ----------
const nav = $<HTMLElement>("nav");
document.querySelectorAll<HTMLElement>(".panel").forEach((panel) => {
  panel.setAttribute("aria-hidden", String(!panel.classList.contains("active")));
});
const artFrame = $<HTMLIFrameElement>("artFrame");
let artReady = false;
let artAtoms: Atom[] = [];
type ArtScene = {
  atoms: Atom[];
  cube?: string;
  iso?: number;
  sides?: string;
  orbital?: Record<string, unknown>;
  label?: string;
};
type ArtSource = { value: string; label: string; disabled?: boolean; reason?: string };
let artScene: ArtScene = { atoms: [] };
let artSettings: Record<string, string | number | boolean> = {
  content: "molecule", quality: "0.55,3", exposure: "1.05", background: "studio",
  dof: "0", surfaceMaterial: "gloss", surfaceOpacity: "1", surfacePhases: "inherit",
};

function pushArtSettings(): void {
  if (!artReady) return;
  artFrame.contentWindow?.postMessage(
    { type: "oqp-art-settings-apply", settings: artSettings }, window.location.origin,
  );
}

function cubeNameFromUrl(url: string): string {
  return decodeURIComponent(url.split("/").pop()?.split("?")[0] ?? "");
}

function setArtStructure(atoms: Atom[]): void {
  if (!atoms.length) return;
  artAtoms = atoms;
  artScene = { atoms };
  pushToArt();
}

function pushToArt(): void {
  if ((!artScene.atoms.length && !artScene.cube) || !artReady) return;
  artFrame.contentWindow?.postMessage(
    { type: "oqp-art-scene", scene: artScene }, window.location.origin,
  );
  pushArtSources();
}

function clearArtScene(): void {
  artAtoms = [];
  artScene = { atoms: [] };
  if (!artReady) return;
  artFrame.contentWindow?.postMessage({ type: "oqp-art-clear" }, window.location.origin);
  pushArtSources();
}

function artSources(): { sources: ArtSource[]; selected: string } {
  const sources: ArtSource[] = artScene.atoms.length
    ? [{ value: "molecule", label: "Molecule" }]
    : [{ value: "", label: "Choose content…", disabled: true }];
  if (currentMolden) {
    for (const option of [...orbitalSel.options]) {
      if (!option.dataset.kind || !orbitalSources.has(option.value)) continue;
      const prefix = option.dataset.kind === "dyson" ? "Dyson" : "SCF";
      sources.push({
        value: `orbital-entry:${encodeURIComponent(option.value)}`,
        label: `${prefix}: ${option.textContent ?? option.value}`,
      });
    }
    const hasScfOrbitals = [...orbitalSel.options]
      .some((option) => option.dataset.kind === "scf" && orbitalSources.has(option.value));
    if (hasScfOrbitals) {
      for (const option of [...mapKind.options]) {
        if (option.hidden || option.value === "mo" || option.value === "dyson") continue;
        sources.push({ value: `orbital:${option.value}`, label: option.textContent ?? option.value });
      }
    }
  }
  if (excitedAnalysis?.available) {
    for (const option of [...excitedMap.options]) {
      sources.push({ value: `excited:${option.value}`, label: option.textContent ?? option.value });
    }
  } else if (excitedAnalysis?.reason) {
    sources.push({
      value: "excited:unavailable",
      label: "NTO / excited-state densities unavailable",
      disabled: true,
      reason: excitedAnalysis.reason,
    });
  }
  for (const name of resultCubeFiles) {
    sources.push({ value: `surface:${name}`, label: `Cube: ${name}` });
  }
  if (activeSurfaceSelection && activeSurfaceSelection.operation !== "display") {
    const { operation, primary, secondary } = activeSurfaceSelection;
    const symbol = operation === "difference" ? "−" : "+";
    sources.push({
      value: "surface-combination",
      label: `${operation === "difference" ? "Difference" : "Sum"}: ${primary} ${symbol} ${secondary}`,
    });
  }
  let selected = sources[0]?.value ?? "";
  if (activeMapSource === "orbital" &&
      (mapKind.value === "mo" || mapKind.value === "dyson") && orbitalSel.value) {
    selected = `orbital-entry:${encodeURIComponent(orbitalSel.value)}`;
  } else if (activeMapSource === "orbital") selected = `orbital:${mapKind.value}`;
  else if (activeMapSource === "excited") selected = `excited:${excitedMap.value}`;
  else if (activeMapSource === "surface" && activeSurfaceSelection &&
           activeSurfaceSelection.operation !== "display") {
    selected = "surface-combination";
  } else if (activeMapSource === "surface" && activeSurfaceSelection?.primary) {
    selected = `surface:${activeSurfaceSelection.primary}`;
  } else if (activeMapSource === "direct" && activeDirectCube) {
    const directName = cubeNameFromUrl(activeDirectCube.url);
    if (resultCubeFiles.includes(directName)) selected = `surface:${directName}`;
  }
  return { sources, selected };
}

function pushArtSources(): void {
  if (!artReady) return;
  artFrame.contentWindow?.postMessage(
    { type: "oqp-art-sources", ...artSources() }, window.location.origin,
  );
  window.setTimeout(pushArtSettings, 0);
}

function pushArtOrbitalStyle(): void {
  const orbital = currentOrbitalStyle();
  artScene = { ...artScene, orbital };
  if (!artReady) return;
  artFrame.contentWindow?.postMessage(
    { type: "oqp-art-style", orbital }, window.location.origin,
  );
}

window.addEventListener("message", (event) => {
  if (event.origin !== window.location.origin || event.source !== artFrame.contentWindow) return;
  if (event.data?.type === "oqp-art-ready") {
    artReady = true;
    artSettings = { ...(event.data.settings ?? {}), ...artSettings };
    if (artScene.atoms.length || artScene.cube) pushToArt();
    else clearArtScene();
    pushArtSettings();
    artFrame.contentWindow?.postMessage({ type: "oqp-art-active",
      active: document.getElementById("panel-art")!.classList.contains("active") }, window.location.origin);
  } else if (event.data?.type === "oqp-art-source-request") {
    showArtSource(String(event.data.value ?? "molecule"));
  } else if (event.data?.type === "oqp-art-settings") {
    artSettings = { ...artSettings, ...(event.data.settings ?? {}) };
  }
});

function showTab(name: string): void {
  document.querySelectorAll<HTMLElement>(".panel").forEach((p) => {
    const active = p.id === `panel-${name}`;
    p.classList.toggle("active", active);
    p.setAttribute("aria-hidden", String(!active));
  });
  nav.querySelectorAll("button").forEach((b) => {
    b.classList.toggle("active", b.dataset.tab === name);
  });
  if (name === "analysis") {
    if (!resultFrame.src) {
      resultViewerReady = false;
      resultFrame.src = "/builder3d.html";
    }
    refreshJobs();
  }
  if (name === "art") {
    if (!selectedJob && !artScene.atoms.length && !artScene.cube) {
      const atoms = parseAtoms(xyzArea.value);
      if (atoms.length) {
        artAtoms = atoms;
        artScene = { atoms };
      }
    }
    if (!artFrame.src) {
      artReady = false;
      artFrame.src = "/art.html";
    } else {
      pushToArt();
    }
  }
  // Replace any deferred scene before activation so reopening builds it once.
  artFrame.contentWindow?.postMessage({ type: "oqp-art-active", active: name === "art" }, window.location.origin);
  if (name === "method") updateInpPreview();
}
nav.addEventListener("click", (e) => {
  const btn = (e.target as HTMLElement).closest("button");
  if (btn?.dataset.tab) showTab(btn.dataset.tab);
});

// ---------- builder ----------
const xyzArea = $<HTMLTextAreaElement>("xyz");
const builderStatus = $<HTMLSpanElement>("builderStatus");
const projectName = $<HTMLInputElement>("projectName");

// WebKit can briefly expose an iframe's default white surface while its WebGL
// backing store is resized. Cover it only for the active resize gesture.
const viewerResizeCovers = Array.from(document.querySelectorAll<HTMLElement>(".viewer-wrap"))
  .map((wrap) => {
    const cover = document.createElement("div");
    cover.className = "viewer-resize-cover";
    cover.setAttribute("aria-hidden", "true");
    wrap.appendChild(cover);
    return cover;
  });
let viewerResizeTimer = 0;
window.addEventListener("resize", () => {
  viewerResizeCovers.forEach((cover) => cover.classList.add("on"));
  clearTimeout(viewerResizeTimer);
  viewerResizeTimer = window.setTimeout(() => {
    viewerResizeCovers.forEach((cover) => cover.classList.remove("on"));
  }, 150);
});

function atomsToText(atoms: Atom[]): string {
  return atoms
    .map(([el, x, y, z]) =>
      `${el.padEnd(2)} ${x.toFixed(6).padStart(11)} ${y.toFixed(6).padStart(11)} ${z.toFixed(6).padStart(11)}`)
    .join("\n");
}

function parseAtoms(text: string): Atom[] {
  const atoms: Atom[] = [];
  for (const raw of text.split("\n")) {
    const parts = raw.trim().split(/\s+/);
    if (parts.length < 4) continue;
    const [el, x, y, z] = [parts[0], +parts[1], +parts[2], +parts[3]];
    if (/^[A-Za-z]{1,2}$/.test(el) && [x, y, z].every(Number.isFinite)) {
      atoms.push([el, x, y, z]);
    }
  }
  return atoms;
}

$<HTMLSelectElement>("sample").addEventListener("change", (e) => {
  const key = (e.target as HTMLSelectElement).value;
  if (SAMPLES[key]) {
    clearPdbSource();
    xyzArea.value = atomsToText(SAMPLES[key].atoms);
    updatePreview();
  }
});

// --- open an existing file: OpenQP inputs/outputs, trajectories, and the
// common interchange formats. Multi-geometry files become a frame slider. ---
let loadedFrames: { label: string; atoms: Atom[] }[] = [];
interface PdbSource {
  name: string;
  text: string;
}
let pdbSource: PdbSource | null = null;
const qmAtoms = new Map<number, string>();

function pdbFilename(name: string): string {
  const stem = name.trim().replace(/^.*[\\/]/, "").replace(/\.(pdb|ent)$/i, "") || "structure";
  return `${stem.replace(/[^A-Za-z0-9._-]/g, "_")}.pdb`;
}

function updateQmAtoms(): void {
  const selected = [...qmAtoms.entries()].sort(([a], [b]) => a - b);
  $<HTMLSpanElement>("qmSelectionStatus").textContent =
    `${selected.length} QM atom${selected.length === 1 ? "" : "s"} selected`;
  $<HTMLSpanElement>("qmAtomList").textContent = selected.length
    ? selected.slice(0, 12).map(([, label]) => label).join(", ") +
      (selected.length > 12 ? ` +${selected.length - 12}` : "")
    : "";
}

function pushPdbPreview(): void {
  if (!pdbSource || !viewerReady) return;
  const target = $<HTMLIFrameElement>("previewFrame").contentWindow;
  target?.postMessage(
    { type: "oqp-file", text: pdbSource.text, format: "pdb" }, window.location.origin);
  target?.postMessage({ type: "oqp-qmmm-mode", enabled: true }, window.location.origin);
}

function setPdbSource(name: string, text: string): void {
  pdbSource = { name: pdbFilename(name), text };
  qmAtoms.clear();
  $<HTMLDivElement>("qmmmCard").style.display = "";
  updateQmAtoms();
  const frame = $<HTMLIFrameElement>("previewFrame");
  frame.src ||= "/builder3d.html";
  if (viewerReady) pushPdbPreview();
}

function clearPdbSource(): void {
  if (!pdbSource) return;
  pdbSource = null;
  qmAtoms.clear();
  $<HTMLDivElement>("qmmmCard").style.display = "none";
  $<HTMLIFrameElement>("previewFrame").contentWindow?.postMessage(
    { type: "oqp-qmmm-mode", enabled: false }, window.location.origin);
}

function qmmmValidation(): string | null {
  if (!pdbSource) return null;
  if (!qmAtoms.size) return "Select at least one QM atom in the PDB viewer before running.";
  if (currentWf.key === "pcm") return "PCM is not available for a PDB QM/MM calculation.";
  return null;
}

function workflowValidation(inputText = ""): string | null {
  if (currentWf.enabled === false) {
    return currentWf.unavailableReason ?? `${currentWf.title} is unavailable.`;
  }
  if (currentWf.key === "scan") {
    if (pdbSource) return "Bond scans currently require an inline molecular geometry, not a PDB QM/MM structure.";
    const atoms = parseAtoms(xyzArea.value);
    const atomA = +fieldValue("scanAtomA");
    const atomB = +fieldValue("scanAtomB");
    const start = +fieldValue("scanStart");
    const end = +fieldValue("scanEnd");
    const points = +fieldValue("scanPoints");
    if (!Number.isInteger(atomA) || !Number.isInteger(atomB) || atomA < 1 || atomB < 1 ||
        atomA > atoms.length || atomB > atoms.length || atomA === atomB) {
      return `Choose two different atom numbers between 1 and ${atoms.length}.`;
    }
    if (!(start > 0) || !(end > 0) || start === end) {
      return "Scan start and end distances must be positive and different.";
    }
    if (!Number.isInteger(points) || points < 2 || points > 101) {
      return "A scan requires 2 to 101 points.";
    }
    return null;
  }
  if (currentWf.key !== "nacme") return null;
  const hasPreviousGeometry = Boolean(fieldValue("nacmeGeometry")) || /(?:^|\n)geom2\s*=/.test(inputText);
  return hasPreviousGeometry
    ? null
    : "NACME requires a previous geometry (.xyz) in addition to the current structure.";
}

function inputValidation(inputText = ""): string | null {
  if (/^\s*\[[^\]\n]+\]\s*$/m.test(inputText)) {
    return "Studio accepts OpenQP canonical .oqp input only; legacy sectioned input is not supported.";
  }
  return qmmmValidation() ?? workflowValidation(inputText);
}

function showFrame(index: number): void {
  const frame = loadedFrames[index - 1];
  if (!frame) return;
  xyzArea.value = atomsToText(frame.atoms);
  $<HTMLSpanElement>("frameLabel").textContent =
    loadedFrames.length > 1
      ? `${index} / ${loadedFrames.length}  ${frame.label}`
      : frame.label;
  updatePreview();
}

$<HTMLInputElement>("frameRange").addEventListener("input", (event) => {
  showFrame(+(event.target as HTMLInputElement).value);
});

$<HTMLInputElement>("fileInput").addEventListener("change", async (event) => {
  const file = (event.target as HTMLInputElement).files?.[0];
  if (!file) return;
  showBusy("Reading structure", file.name);
  try {
    const body = new FormData();
    body.append("file", file);
    const res = await fetch("/api/structure/open", { method: "POST", body });
    if (!res.ok) {
      builderStatus.textContent = `${file.name}: ${(await res.json()).detail ?? res.status}`;
      return;
    }
    const data = await res.json();
    if (/\.(pdb|ent)$/i.test(file.name)) {
      // Mol* reads PDB natively; going through plain coordinates would throw
      // away the residues and chains a cartoon needs.
      const text = await file.text();
      setPdbSource(file.name, text);
    } else {
      clearPdbSource();
    }
    loadedFrames = data.frames;
    const slider = $<HTMLInputElement>("frameRange");
    slider.max = String(loadedFrames.length);
    slider.value = String(loadedFrames.length);   // land on the final geometry
    $<HTMLDivElement>("frameRow").style.display =
      loadedFrames.length > 1 ? "flex" : "none";
    showFrame(loadedFrames.length);
    builderStatus.textContent =
      `${data.format} · ${loadedFrames[0].atoms.length} atoms` +
      (loadedFrames.length > 1 ? ` · ${loadedFrames.length} frames` : "");
  } finally {
    hideBusy();
  }
});

// Proteins and nucleic acids come from the RCSB Protein Data Bank, not from
// PubChem; the file is handed to the viewer whole so chains and residues
// survive for the cartoon representation.
$<HTMLButtonElement>("pdbFetch").addEventListener("click", async () => {
  const code = $<HTMLInputElement>("pdbCode").value.trim();
  if (!code) return;
  showBusy("Fetching from the PDB", `entry ${code.toUpperCase()}…`);
  try {
    const response = await fetch(`/api/pdb/${encodeURIComponent(code)}`);
    if (!response.ok) throw new Error((await response.json()).detail ?? response.statusText);
    const data = await response.json();
    setPdbSource(`${data.code}.pdb`, data.pdb);

    // Keep coordinates visible for inspection, but PDB runs use the original
    // topology and a deliberately selected QM region rather than all atoms.
    const body = new FormData();
    body.append("file", new File([data.pdb], `${data.code}.pdb`));
    const parsed = await fetch("/api/structure/open", { method: "POST", body });
    if (parsed.ok) {
      const structure = await parsed.json();
      loadedFrames = structure.frames;
      xyzArea.value = atomsToText(loadedFrames[0].atoms);
      await updatePreview();
      builderStatus.textContent =
        `${data.code} · ${loadedFrames[0].atoms.length} atoms from the PDB`;
    }
  } catch (err) {
    const message = (err as Error).message;
    builderStatus.textContent = `PDB: ${message}`;
    if (/certificate/i.test(message)) {
      hideBusy();
      await openNetworkSettings();
    }
  } finally {
    hideBusy();
  }
});

$<HTMLButtonElement>("pubchemFetch").addEventListener("click", async () => {
  const name = $<HTMLInputElement>("pubchemName").value.trim();
  if (!name) return;
  showBusy("Fetching from PubChem", `looking up “${name}”…`);
  try {
    const res = await fetch(`/api/pubchem/${encodeURIComponent(name)}`);
    if (!res.ok) throw new Error((await res.json()).detail ?? res.statusText);
    const data = await res.json();
    clearPdbSource();
    xyzArea.value = atomsToText(data.atoms);
    builderStatus.textContent = `${data.atoms.length} atoms loaded from PubChem`;
    await updatePreview();
  } catch (err) {
    const message = (err as Error).message;
    builderStatus.textContent = `PubChem: ${message}`;
    // A certificate failure has a fix, so take the user straight to it.
    if (/certificate/i.test(message)) {
      hideBusy();
      await openNetworkSettings();
    }
  } finally {
    hideBusy();
  }
});

// A spinner for the steps that take a noticeable moment: the first
// sketch-to-3D conversion loads RDKit, and PubChem is a network round trip.
const busy = $<HTMLDivElement>("busy");
let busyTimer: number | undefined;

function showBusy(message: string, sub = ""): void {
  $<HTMLDivElement>("busyMsg").textContent = message;
  const subEl = $<HTMLDivElement>("busySub");
  subEl.textContent = sub;
  busy.classList.add("on");
  const started = Date.now();
  clearInterval(busyTimer);
  busyTimer = window.setInterval(() => {
    const seconds = Math.round((Date.now() - started) / 1000);
    subEl.textContent = seconds >= 3 ? `${sub} ${seconds}s` : sub;
  }, 1000);
}

function hideBusy(): void {
  clearInterval(busyTimer);
  busy.classList.remove("on");
}

// --- 2D sketcher (Ketcher in a modal; RDKit turns the sketch into 3D) ---
const sketchOverlay = $<HTMLDivElement>("sketchOverlay");
$<HTMLButtonElement>("sketchBtn").addEventListener("click", () => {
  const frame = $<HTMLIFrameElement>("sketchFrame");
  if (!frame.src) frame.src = "/sketcher.html";
  sketchOverlay.style.display = "";
});
$<HTMLButtonElement>("sketchClose").addEventListener("click", () => {
  sketchOverlay.style.display = "none";
});
window.addEventListener("message", async (event) => {
  if (event.origin !== window.location.origin) return;
  if (event.data?.type !== "oqp-sketch") return;
  sketchOverlay.style.display = "none";
  showBusy("Building 3D structure", "embedding and optimizing with RDKit…");
  try {
    const res = await fetch("/api/structure3d", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ molfile: event.data.molfile }),
    });
    if (!res.ok) {
      builderStatus.textContent = `3D build failed: ${(await res.json()).detail ?? res.status}`;
      return;
    }
    const data = await res.json();
    xyzArea.value = atomsToText(data.atoms);
    builderStatus.textContent = `${data.atoms.length} atoms from sketch (RDKit ETKDG + MMFF)`;
    await updatePreview();
  } finally {
    hideBusy();
  }
});

// Display options are shared by both viewers so a molecule looks the same
// wherever it is shown. "auto" keeps Mol*'s own preset, which is what the app
// has always drawn — and what turns a PDB into a protein cartoon.
const displayStyle: Record<string, string> = {
  type: "auto",
  color: "element-symbol",
  labels: "none",
};

function pushStyle(): void {
  const message = { type: "oqp-style", style: displayStyle };
  for (const id of ["previewFrame", "resultFrame"]) {
    const frame = document.getElementById(id) as HTMLIFrameElement | null;
    if (frame?.src) frame.contentWindow?.postMessage(message, window.location.origin);
  }
}

const rendering: Record<string, boolean | number | string> = {
  occlusion: true,
  outline: true,
  shadow: false,
  depthCue: true,
  antialiasing: "smaa",
  occlusionQuality: "balanced",
  background: 0x1b1e24,
};

function pushToViewers(message: Record<string, unknown>): void {
  for (const id of ["previewFrame", "resultFrame"]) {
    const frame = document.getElementById(id) as HTMLIFrameElement | null;
    if (frame?.src) frame.contentWindow?.postMessage(message, window.location.origin);
  }
}

document.querySelectorAll<HTMLInputElement | HTMLSelectElement>(".render-input").forEach((input) => {
  input.addEventListener("change", () => {
    const key = input.dataset.render!;
    rendering[key] = input instanceof HTMLInputElement && input.type === "checkbox"
      ? input.checked
      : key === "background" ? parseInt(input.value, 16) : input.value;
    document.querySelectorAll<HTMLInputElement | HTMLSelectElement>(
      `.render-input[data-render="${key}"]`,
    ).forEach((twin) => {
      if (twin instanceof HTMLInputElement && twin.type === "checkbox") {
        twin.checked = input instanceof HTMLInputElement && input.checked;
      } else {
        twin.value = input.value;
      }
    });
    pushToViewers({ type: "oqp-rendering", rendering });
  });
});

// Orbital phase colours, opacity and surface kind, as MacMolPlt and IQmol offer.
function currentOrbitalStyle(): Record<string, unknown> {
  const [positive, negative] = $<HTMLSelectElement>("moColors").value.split(",");
  return {
    positive: parseInt(positive, 16),
    negative: parseInt(negative, 16),
    alpha: +$<HTMLInputElement>("moAlpha").value,
    visuals: $<HTMLSelectElement>("moVisual").value.split(","),
    sides: $<HTMLSelectElement>("moSides").value,
  };
}

function pushOrbitalStyle(): void {
  pushToViewers({ type: "oqp-orbital-style", orbital: currentOrbitalStyle() });
  pushArtOrbitalStyle();
}

// Ready-made orbital looks. "studio" is what the app has always drawn and
// stays the default; the MacMolPlt entries reproduce that program's opaque
// red/blue lobes, solid or as a mesh.
const MO_PRESETS: Record<string, { colors: string; visual: string; alpha: string }> = {
  studio: { colors: "4f8fdd,dd6a4f", visual: "solid", alpha: "0.85" },
  macmolplt: { colors: "d64a3b,3b6fd6", visual: "solid", alpha: "1" },
  "macmolplt-mesh": { colors: "d64a3b,3b6fd6", visual: "wireframe", alpha: "1" },
  iqmol: { colors: "3b6fd6,e0b93a", visual: "solid", alpha: "0.55" },
  print: { colors: "d9dde4,7b8494", visual: "solid,wireframe", alpha: "1" },
};

const moPreset = $<HTMLSelectElement>("moPreset");
const moColors = $<HTMLSelectElement>("moColors");
const moVisual = $<HTMLSelectElement>("moVisual");
const moAlpha = $<HTMLInputElement>("moAlpha");

moPreset.addEventListener("change", () => {
  const preset = MO_PRESETS[moPreset.value];
  if (!preset) return;                       // "custom" leaves the controls alone
  moColors.value = preset.colors;
  moVisual.value = preset.visual;
  moAlpha.value = preset.alpha;
  redrawActiveVolume();                      // colours are baked into the volume
});

for (const id of ["moColors", "moVisual", "moAlpha"]) {
  $<HTMLElement>(id).addEventListener("change", () => {
    moPreset.value = "custom";
    // Colours are baked in when the volume loads, so redraw for those.
    if (id === "moColors") redrawActiveVolume();
    else pushOrbitalStyle();
  });
}

document.querySelectorAll<HTMLSelectElement>(".style-input").forEach((select) => {
  select.addEventListener("change", () => {
    displayStyle[select.dataset.style!] = select.value;
    // Keep the duplicate control in the other tab in step.
    document.querySelectorAll<HTMLSelectElement>(
      `.style-input[data-style="${select.dataset.style}"]`,
    ).forEach((twin) => { twin.value = select.value; });
    pushStyle();
  });
});

// The analysis viewer reports what the current-density layer actually drew.
window.addEventListener("message", (event) => {
  if (event.origin !== window.location.origin) return;
  if (event.data?.type !== "oqp-current-vectors") return;
  const count = Number(event.data.count) || 0;
  const peak = Number(event.data.maxMagnitude);
  $<HTMLDivElement>("currentVectorStatus").textContent = count
    ? `${count} current vectors on the ACID isosurface`
      + (Number.isFinite(peak) ? `; max |J| = ${peak.toExponential(2)} a.u.` : "")
    : "no current vectors on this isosurface — lower the isovalue or raise the density";
});

// The preview iframe hosts one long-lived Mol* viewer; structures are pushed
// into it rather than reloading the page, which would re-initialize WebGL.
let viewerReady = false;
window.addEventListener("message", (event) => {
  if (event.origin === window.location.origin && event.data?.type === "oqp-viewer-ready" &&
      event.source === $<HTMLIFrameElement>("previewFrame").contentWindow) {
    viewerReady = true;
    if (pdbSource) pushPdbPreview();
    else updatePreview();
  }
});

async function updatePreview(): Promise<void> {
  invalidateSymmetryIfCoordinatesChanged();
  const atoms = parseAtoms(xyzArea.value);
  if (atoms.length && $<HTMLElement>("panel-builder").classList.contains("active")) {
    setArtStructure(atoms);
  }
  if (pdbSource) {
    pushPdbPreview();
    return;
  }
  if (!atoms.length) return;
  const xyz = `${atoms.length}\nOQP Studio preview\n${atomsToText(atoms)}\n`;
  const frame = $<HTMLIFrameElement>("previewFrame");
  if (!frame.src) {
    frame.src = "/builder3d.html";   // pushes arrive once it reports ready
    return;
  }
  if (!viewerReady) return;
  frame.contentWindow?.postMessage(
    { type: "oqp-structure", xyz },
    window.location.origin,
  );
}
$<HTMLButtonElement>("previewBtn").addEventListener("click", updatePreview);

type SymmetryResult = {
  point_group: string;
  tolerance_angstrom: number;
  max_deviation_angstrom: number;
  operations: string[];
  operation_count: number;
  equivalent_atoms: number[][];
  aligned_atoms: Atom[];
};
let symmetryResult: SymmetryResult | null = null;
let symmetrySource = "";
let symmetryPendingSource = "";
let symmetryTolerance: number | null = null;
let symmetryPendingTolerance: number | null = null;
let symmetryRequestId = 0;

function invalidateSymmetry(message: string): void {
  if (!symmetryResult && !symmetryPendingSource) return;
  symmetryRequestId += 1;
  symmetryPendingSource = "";
  symmetryPendingTolerance = null;
  symmetryResult = null;
  symmetrySource = "";
  symmetryTolerance = null;
  $<HTMLButtonElement>("symmetryAlign").disabled = true;
  $<HTMLDivElement>("symmetryResult").textContent = message;
}

function invalidateSymmetryIfCoordinatesChanged(): void {
  const current = xyzArea.value;
  if (symmetrySource === current || symmetryPendingSource === current) return;
  invalidateSymmetry("Coordinates changed; analyze symmetry again.");
}

xyzArea.addEventListener("input", invalidateSymmetryIfCoordinatesChanged);
$<HTMLInputElement>("symmetryTolerance").addEventListener("change", () => {
  const current = +$<HTMLInputElement>("symmetryTolerance").value;
  if (symmetryTolerance === current || symmetryPendingTolerance === current) return;
  invalidateSymmetry("Tolerance changed; analyze symmetry again.");
});

function pointGroupLabel(value: string): string {
  return value === "Dinfh" ? "D∞h" : value === "Cinfv" ? "C∞v" : value;
}

$<HTMLButtonElement>("symmetryAnalyze").addEventListener("click", async () => {
  const status = $<HTMLDivElement>("symmetryResult");
  const tolerance = +$<HTMLInputElement>("symmetryTolerance").value;
  const source = xyzArea.value;
  const requestId = ++symmetryRequestId;
  symmetryPendingSource = source;
  symmetryPendingTolerance = tolerance;
  symmetryResult = null;
  $<HTMLButtonElement>("symmetryAlign").disabled = true;
  status.textContent = "Analyzing coordinates…";
  let response: Response;
  try {
    response = await fetch("/api/symmetry", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ xyz: source, tolerance }),
    });
  } catch (error) {
    if (requestId !== symmetryRequestId || source !== xyzArea.value ||
        tolerance !== +$<HTMLInputElement>("symmetryTolerance").value) return;
    symmetryPendingSource = "";
    symmetryPendingTolerance = null;
    status.textContent = `symmetry analysis failed: ${error instanceof Error ? error.message : String(error)}`;
    return;
  }
  const payload = await response.json().catch(() => null);
  if (requestId !== symmetryRequestId || source !== xyzArea.value ||
      tolerance !== +$<HTMLInputElement>("symmetryTolerance").value) return;
  symmetryPendingSource = "";
  symmetryPendingTolerance = null;
  if (!response.ok) {
    status.textContent = payload?.detail ?? `symmetry analysis failed (${response.status})`;
    return;
  }
  if (!payload) {
    status.textContent = "symmetry analysis returned an invalid response";
    return;
  }
  symmetryResult = payload;
  symmetrySource = source;
  symmetryTolerance = tolerance;
  const equivalents = symmetryResult!.equivalent_atoms
    .filter((group) => group.length > 1)
    .map((group) => group.join(", "))
    .join("; ") || "none";
  status.textContent = `Likely ${pointGroupLabel(symmetryResult!.point_group)} at ` +
    `${symmetryResult!.tolerance_angstrom.toFixed(4)} Å tolerance; ` +
    `${symmetryResult!.operation_count} accepted operations, maximum residual ` +
    `${symmetryResult!.max_deviation_angstrom.toExponential(2)} Å. ` +
    `Equivalent atoms: ${equivalents}.`;
  $<HTMLButtonElement>("symmetryAlign").disabled = false;
});

$<HTMLButtonElement>("symmetryAlign").addEventListener("click", () => {
  invalidateSymmetryIfCoordinatesChanged();
  if (!symmetryResult) return;
  clearPdbSource();
  const aligned = atomsToText(symmetryResult.aligned_atoms);
  xyzArea.value = aligned;
  symmetrySource = aligned;
  builderStatus.textContent = "centered and aligned to the principal axes";
  void updatePreview();
});

function buildSampleList(): void {
  const select = $<HTMLSelectElement>("sample");
  for (const [key, sample] of Object.entries(SAMPLES)) {
    const option = document.createElement("option");
    option.value = key;
    option.textContent = sample.label;
    select.appendChild(option);
  }
}

// ---------- method: workflow catalog ----------
// First pick a workflow (like app.openqp.org); details open on demand.
interface Workflow {
  key: string;
  title: string;
  desc: string;
  detail: string;
  theories: string[];
  defaultTheory?: string;
  enabled?: boolean;
  unavailableReason?: string;
}

const ALL_THEORIES = [
  "hf", "dft", "mp2", "ccsd", "ccsd_t", "fci", "casci", "casscf", "sa-casscf",
  "caspt2", "ms-caspt2", "xms-caspt2", "nevpt2", "sc-nevpt2", "mrmp2", "mcqdpt2", "xmcqdpt2",
  "tddft", "tda", "sf", "mrsf", "umrsf",
];
const DFT_RESPONSE = ["tddft", "tda", "sf", "mrsf", "umrsf"];
const LABELED_RESPONSE = ["tddft", "tda", "mrsf"];
const MRSF_ONLY = ["mrsf"];
const GRADIENT_THEORIES = ["hf", "dft", "mp2", ...DFT_RESPONSE, "casscf", "sa-casscf", "caspt2", "ms-caspt2", "xms-caspt2", "nevpt2", "sc-nevpt2"];

const WORKFLOWS: Workflow[] = [
  { key: "energy", title: "Single-point energy", desc: "Electronic energy at the supplied geometry.", detail: "Keeps the nuclear coordinates fixed and reports the energy and electronic-state data required for subsequent analysis.", theories: ALL_THEORIES, defaultTheory: "dft" },
  { key: "grad", title: "Energy gradient", desc: "Nuclear energy gradient for the ground or selected electronic state.", detail: "Use the gradient to characterize forces at a fixed geometry or as the derivative input for structural optimization and path calculations.", theories: GRADIENT_THEORIES, defaultTheory: "dft" },
  { key: "opt", title: "Geometry optimization", desc: "Optimize the molecular structure on the selected potential-energy surface.", detail: "For response theories, choose S0 or a specified excited state. Convergence is assessed from energy changes, gradients, and nuclear displacements.", theories: GRADIENT_THEORIES, defaultTheory: "dft" },
  { key: "scan", title: "Bond-distance scan", desc: "Map the energy while changing one internuclear distance.", detail: "Choose two atoms, a distance interval, and the number of points. A rigid scan evaluates fixed structures; a relaxed scan optimizes every point while OpenQP holds the selected bond distance fixed.", theories: GRADIENT_THEORIES, defaultTheory: "dft" },
  { key: "hess", title: "Frequencies (Hessian)", desc: "Vibrational frequencies, normal modes, and thermochemical quantities.", detail: "Calculates the second derivative matrix at the supplied structure, then obtains IR intensities, zero-point energy, and temperature-dependent thermochemistry.", theories: ["hf", "dft", ...LABELED_RESPONSE], defaultTheory: "dft" },
  { key: "prop", title: "Molecular properties", desc: "Electronic populations, multipoles, and selected state properties.", detail: "Use this fixed-geometry calculation to inspect charge distribution and electrostatic or state-resolved quantities without changing the molecular structure.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "nmr", title: "NMR shielding", desc: "Ground-state nuclear magnetic shielding constants.", detail: "Computes the shielding tensor at the supplied geometry; isotropic shifts and tensor components can then be compared among chemically distinct nuclei.", theories: ["hf", "dft"], defaultTheory: "dft" },
  { key: "pcm", title: "PCM solvation", desc: "Reference-SCF energy in a polarizable continuum solvent.", detail: "Temporarily unavailable while ddX is disabled in distributed OpenQP engines.", theories: ["hf", "dft"], defaultTheory: "hf", enabled: false, unavailableReason: "PCM is temporarily unavailable because ddX is disabled in this release." },
  { key: "abs", title: "Vertical excited states", desc: "Excitation energies and oscillator strengths at the S0 geometry.", detail: "The nuclei remain at the reference structure, so the state energies and transition moments provide the stick spectrum used for absorption broadening.", theories: DFT_RESPONSE, defaultTheory: "mrsf" },
  { key: "exgrad", title: "Excited-state gradient", desc: "Nuclear gradient of a selected excited electronic state.", detail: "Evaluates the derivative on the chosen state at the present geometry; it is the starting point for excited-state structural relaxation or dynamics.", theories: LABELED_RESPONSE, defaultTheory: "mrsf" },
  { key: "exopt", title: "Excited-state optimization", desc: "Optimize the geometry on a selected excited-state surface.", detail: "Choose the target state explicitly. The resulting relaxed structure is appropriate for emission or excited-state absorption analysis when the state remains well characterized.", theories: LABELED_RESPONSE, defaultTheory: "mrsf" },
  { key: "meci", title: "MECI search", desc: "Locate a minimum-energy conical intersection between two electronic states.", detail: "Optimizes the mean energy while driving the selected state energy gap toward zero, yielding a crossing geometry relevant to internal conversion.", theories: ["mrsf", "casscf", "sa-casscf"], defaultTheory: "mrsf" },
  { key: "mecp", title: "MECP search", desc: "Locate a minimum-energy crossing point between states of different spin.", detail: "Optimizes the crossing geometry while enforcing equality of the selected state energies, for example between singlet and triplet surfaces.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "tci", title: "Three-state intersection", desc: "Search for a geometry where three selected electronic states approach degeneracy.", detail: "The calculation treats the three state energies together and is intended for multistate crossing regions rather than a single pairwise crossing.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "ts", title: "Transition-state search", desc: "Find a first-order saddle point on the selected potential-energy surface.", detail: "A successful structure has one unstable normal coordinate. Follow with a Hessian and an IRC calculation to verify the reaction connection.", theories: GRADIENT_THEORIES, defaultTheory: "dft" },
  { key: "irc", title: "Intrinsic reaction coordinate", desc: "Trace the steepest-descent path from a transition-state structure.", detail: "Integrates downhill in both directions to identify the connected reactant and product valleys on the same potential-energy surface.", theories: ["hf", "dft"], defaultTheory: "dft" },
  { key: "mep", title: "Minimum-energy path", desc: "Trace a reaction path on the selected electronic-state surface.", detail: "Uses state-specific gradients to follow a path through nuclear-coordinate space and records the energy profile along that path.", theories: GRADIENT_THEORIES, defaultTheory: "mrsf" },
  { key: "neb", title: "Nudged elastic band", desc: "Optimize a discretized path between reactant and product structures.", detail: "Uses intermediate images and spring forces to resolve a minimum-energy path when both endpoint geometries are available.", theories: ["hf", "dft"], defaultTheory: "dft" },
  { key: "nac", title: "Nonadiabatic coupling", desc: "Derivative coupling vector between two selected electronic states.", detail: "The coupling identifies nuclear motions that mix the states and is needed near crossings or for nonadiabatic nuclear dynamics.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "nacme", title: "NACME", desc: "Nonadiabatic coupling matrix element between two selected states.", detail: "Reports the state-to-state coupling quantity for the specified electronic-state pair at the current nuclear geometry.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "soc", title: "Spin-orbit coupling", desc: "Spin-orbit matrix elements between singlet and triplet states.", detail: "Use the coupling and state energy gaps to assess intersystem crossing pathways and spin-mixed electronic states.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "ekt", title: "Ionization / EA (EKT)", desc: "Ionization and electron-affinity states from the extended Koopmans theorem.", detail: "Reports IP/EA state energies and Dyson orbitals. Dyson strength gives an orbital occupation contribution through twice the reported strength.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
  { key: "namd", title: "Nonadiabatic dynamics", desc: "Surface-hopping nuclear dynamics across coupled electronic states.", detail: "Propagates nuclear trajectories while evaluating electronic populations and state transitions; choose the initial state, time step, and total number of steps.", theories: MRSF_ONLY, defaultTheory: "mrsf" },
];

let currentWf: Workflow = WORKFLOWS[0];

const theorySel = $<HTMLSelectElement>("theory");
const functionalSel = $<HTMLSelectElement>("functional");
const functionalCustom = $<HTMLInputElement>("functionalCustom");
const basisSel = $<HTMLSelectElement>("basis");
const basisCustom = $<HTMLInputElement>("basisCustom");
const optionsCard = $<HTMLDivElement>("optionsCard");
const hessTypeSel = $<HTMLSelectElement>("hessType");
const hessTypeHint = $<HTMLSpanElement>("hessTypeHint");

function syncWorkflowOptions(): void {
  document.querySelectorAll<HTMLElement>(".wf-opt").forEach((row) => {
    const workflows = (row.dataset.for ?? "").split(" ");
    const theories = (row.dataset.theories ?? "").split(" ").filter(Boolean);
    row.classList.toggle("on", (workflows.includes("*") || workflows.includes(currentWf.key)) &&
      (!theories.length || theories.includes(theorySel.value)));
  });
}

function selectWorkflow(wf: Workflow): void {
  const firstSelection = optionsCard.style.display === "none";
  currentWf = wf;
  $<HTMLSelectElement>("workflowSel").value = wf.key;
  $<HTMLDivElement>("workflowSummary").textContent = wf.desc;
  $<HTMLDivElement>("workflowDetail").textContent = wf.detail;
  for (const option of Array.from(theorySel.options)) {
    option.disabled = !wf.theories.includes(option.value);
  }
  if (firstSelection || !wf.theories.includes(theorySel.value)) {
    theorySel.value = wf.defaultTheory ?? wf.theories[0];
  }
  if (wf.key === "opt") $<HTMLInputElement>("targetState").value = "0";
  if (wf.key === "exopt" && fieldValue("targetState") === "0") {
    $<HTMLInputElement>("targetState").value = "1";
  }
  if (wf.key === "hess" && ["hf", "dft"].includes(theorySel.value)) {
    hessTypeSel.value = "analytical";
  }
  $<HTMLSpanElement>("wfLabel").textContent = `— ${wf.title}`;
  optionsCard.style.display = "";
  syncFieldStates();
  syncPcmReference();
  syncWorkflowOptions();
  syncWorkflowDetails();
  updateInpPreview();
}

function syncWorkflowDetails(): void {
  document.querySelectorAll<HTMLElement>(".wf-detail").forEach((section) => {
    const workflows = (section.dataset.for ?? "").split(" ");
    const theories = (section.dataset.theories ?? "").split(" ").filter(Boolean);
    section.classList.toggle("on", (workflows.includes("*") || workflows.includes(currentWf.key)) &&
      (!theories.length || theories.includes(theorySel.value)));
  });
}

function buildWorkflowSelect(): void {
  const select = $<HTMLSelectElement>("workflowSel");
  const render = (query = "") => {
    const normalized = query.trim().toLowerCase();
    const visible = WORKFLOWS.filter((wf) =>
      !normalized || `${wf.title} ${wf.desc} ${wf.detail} ${wf.key}`.toLowerCase().includes(normalized));
    select.replaceChildren(...visible.map((wf) => {
      const option = document.createElement("option");
      option.value = wf.key;
      option.textContent = wf.enabled === false ? `${wf.title} (unavailable)` : wf.title;
      option.disabled = wf.enabled === false;
      return option;
    }));
    if (visible.some((wf) => wf.key === currentWf.key)) select.value = currentWf.key;
  };
  render();
  $<HTMLInputElement>("workflowSearch").addEventListener("input", (event) => {
    render((event.target as HTMLInputElement).value);
  });
  select.addEventListener("change", () => {
    const workflow = WORKFLOWS.find((wf) => wf.key === select.value);
    if (workflow) selectWorkflow(workflow);
  });
}

type Recipe = {
  schema?: string;
  id: string;
  name: string;
  description: string;
  tags: string[];
  available: boolean;
  unavailable_reason: string;
  workflow: string;
  method: Record<string, string | number | boolean>;
  execution: Record<string, string | number | boolean>;
  analysis: Record<string, string | number | boolean>;
  art: Record<string, string | number | boolean>;
  postprocess: { enabled: boolean; code: string; timeout_seconds: number };
};
type RecipeRecord = { recipe: Recipe; builtin: boolean; trusted: boolean };

let recipeRecords: RecipeRecord[] = [];
let selectedRecipeId = "";
let activeRecipeId = "";
let pendingRecipeAnalysis: Record<string, string | number | boolean> = {};

function recipeId(name: string): string {
  const slug = name.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
  return (slug || `recipe-${Date.now()}`).slice(0, 64);
}

function methodSettings(): Record<string, string | boolean> {
  const values: Record<string, string | boolean> = {};
  optionsCard.querySelectorAll<HTMLInputElement | HTMLSelectElement | HTMLTextAreaElement>(
    "input[id], select[id], textarea[id]",
  ).forEach((control) => {
    values[control.id] = control instanceof HTMLInputElement && control.type === "checkbox"
      ? control.checked : control.value;
  });
  return values;
}

function setControls(values: Record<string, string | number | boolean>): void {
  for (const [id, value] of Object.entries(values)) {
    const control = document.getElementById(id) as HTMLInputElement | HTMLSelectElement | null;
    if (!control || !optionsCard.contains(control)) continue;
    if (control instanceof HTMLInputElement && control.type === "checkbox") {
      control.checked = Boolean(value);
    } else {
      control.value = String(value);
    }
    control.dispatchEvent(new Event("input", { bubbles: true }));
    control.dispatchEvent(new Event("change", { bubbles: true }));
  }
}

function renderRecipeSelect(query = ""): void {
  const select = $<HTMLSelectElement>("recipeSel");
  const normalized = query.trim().toLowerCase();
  const records = recipeRecords.filter(({ recipe }) =>
    !normalized || `${recipe.name} ${recipe.description} ${recipe.tags.join(" ")} ${recipe.workflow}`
      .toLowerCase().includes(normalized));
  const placeholder = document.createElement("option");
  placeholder.value = "";
  placeholder.textContent = records.length ? "Choose a recipe…" : "No matching recipes";
  select.replaceChildren(placeholder, ...records.map(({ recipe, builtin }) => {
    const option = document.createElement("option");
    option.value = recipe.id;
    option.textContent = `${recipe.name}${builtin ? "" : " · custom"}`;
    return option;
  }));
  if (records.some(({ recipe }) => recipe.id === selectedRecipeId)) select.value = selectedRecipeId;
}

function showRecipeRecord(): void {
  const id = $<HTMLSelectElement>("recipeSel").value;
  selectedRecipeId = id;
  const record = recipeRecords.find(({ recipe }) => recipe.id === id);
  const apply = $<HTMLButtonElement>("recipeApply");
  const exportButton = $<HTMLButtonElement>("recipeExport");
  const deleteButton = $<HTMLButtonElement>("recipeDelete");
  const pythonReview = $<HTMLDetailsElement>("recipePythonReview");
  apply.disabled = !record || !record.recipe.available;
  exportButton.disabled = !record;
  deleteButton.disabled = !record || record.builtin;
  pythonReview.style.display = record?.recipe.postprocess.code ? "" : "none";
  $<HTMLTextAreaElement>("recipePythonSource").value = record?.recipe.postprocess.code ?? "";
  $<HTMLButtonElement>("recipeTrust").disabled = !record || record.builtin || record.trusted;
  if (!record) {
    $<HTMLDivElement>("recipeSummary").textContent =
      "A recipe preserves calculation, execution, analysis, and Art settings.";
    $<HTMLDivElement>("recipeDetail").textContent = "";
    return;
  }
  $<HTMLDivElement>("recipeSummary").textContent = record.recipe.description;
  const detail = $<HTMLDivElement>("recipeDetail");
  detail.className = `detail${record.recipe.available ? "" : " recipe-unavailable"}`;
  detail.textContent = record.recipe.available
    ? `${record.recipe.tags.join(" · ")}${record.recipe.postprocess.enabled
      ? ` · Python ${record.trusted ? "trusted" : "disabled until trusted"}` : ""}`
    : record.recipe.unavailable_reason;
}

async function loadRecipes(): Promise<void> {
  const response = await fetch("/api/recipes");
  if (!response.ok) throw new Error(`recipe catalog failed (${response.status})`);
  recipeRecords = await response.json();
  renderRecipeSelect($<HTMLInputElement>("recipeSearch").value);
  showRecipeRecord();
}

function applyRecipe(recipe: Recipe): void {
  const workflow = WORKFLOWS.find((candidate) => candidate.key === recipe.workflow);
  if (!workflow || !recipe.available) return;
  activeRecipeId = recipe.id;
  const workflowSearch = $<HTMLInputElement>("workflowSearch");
  workflowSearch.value = "";
  workflowSearch.dispatchEvent(new Event("input"));
  selectWorkflow(workflow);
  setControls(recipe.method);
  if (recipe.execution.runner &&
      runnerSelect.querySelector(`option[value="${String(recipe.execution.runner)}"]`)) {
    runnerSelect.value = String(recipe.execution.runner);
    runnerSelect.dispatchEvent(new Event("change"));
  }
  if (recipe.execution.threads) threadsInput.value = String(recipe.execution.threads);
  pendingRecipeAnalysis = { ...recipe.analysis };
  artSettings = { ...artSettings, ...recipe.art };
  pushArtSettings();
  updateInpPreview();
  $<HTMLDivElement>("recipeStatus").textContent = `Applied ${recipe.name}`;
}

function applyPendingAnalysis(): void {
  const values = pendingRecipeAnalysis;
  const kind = String(values.spectrum_kind ?? "");
  if (kind && [...specKind.options].some((option) => option.value === kind)) specKind.value = kind;
  if (values.spectrum_shape) specShape.value = String(values.spectrum_shape);
  if (values.spectrum_width) specWidth.value = String(values.spectrum_width);
  if (values.spectrum_state) specState.value = String(values.spectrum_state);
  syncSpectrumControls();
}

function applyPendingPropertyMap(): void {
  const key = String(pendingRecipeAnalysis.map_kind ?? "");
  if (!key || !atomicProperties.has(key)) return;
  $<HTMLSelectElement>("atomicProperty").value = key;
  showAtomicProperty(key);
}

async function saveCurrentRecipe(): Promise<void> {
  const name = $<HTMLInputElement>("recipeName").value.trim();
  if (!name) {
    $<HTMLDivElement>("recipeStatus").textContent = "Give the recipe a name.";
    return;
  }
  const pythonEnabled = $<HTMLInputElement>("recipePythonEnabled").checked;
  const code = $<HTMLTextAreaElement>("recipePython").value;
  const recipe: Recipe = {
    schema: "oqp-studio-recipe/v1", id: recipeId(name), name,
    description: `End-to-end ${currentWf.title.toLowerCase()} recipe.`,
    tags: $<HTMLInputElement>("recipeTags").value.split(",").map((tag) => tag.trim()).filter(Boolean),
    available: true, unavailable_reason: "", workflow: currentWf.key,
    method: methodSettings(),
    execution: { runner: runnerSelect.value, threads: threadCount() },
    analysis: {
      spectrum_kind: specKind.value, spectrum_shape: specShape.value,
      spectrum_width: specWidth.value, spectrum_state: specState.value,
    },
    art: { ...artSettings },
    postprocess: { enabled: pythonEnabled, code, timeout_seconds: 60 },
  };
  const response = await fetch("/api/recipes", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      recipe, trust_python: pythonEnabled && $<HTMLInputElement>("recipePythonTrust").checked,
      replace: false,
    }),
  });
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    $<HTMLDivElement>("recipeStatus").textContent = detail?.detail ?? `save failed (${response.status})`;
    return;
  }
  selectedRecipeId = recipe.id;
  await loadRecipes();
  $<HTMLSelectElement>("recipeSel").value = recipe.id;
  showRecipeRecord();
  $<HTMLDivElement>("recipeStatus").textContent = `Saved ${name}`;
}

function installRecipeControls(): void {
  $<HTMLInputElement>("recipeSearch").addEventListener("input", (event) => {
    renderRecipeSelect((event.target as HTMLInputElement).value);
    showRecipeRecord();
  });
  $<HTMLSelectElement>("recipeSel").addEventListener("change", showRecipeRecord);
  $<HTMLButtonElement>("recipeApply").addEventListener("click", () => {
    const record = recipeRecords.find(({ recipe }) => recipe.id === selectedRecipeId);
    if (record) applyRecipe(record.recipe);
  });
  $<HTMLButtonElement>("recipeSave").addEventListener("click", () => { void saveCurrentRecipe(); });
  $<HTMLButtonElement>("recipeTrust").addEventListener("click", async () => {
    const record = recipeRecords.find(({ recipe }) => recipe.id === selectedRecipeId);
    if (!record || record.builtin || record.trusted || !record.recipe.postprocess.code) return;
    const response = await fetch("/api/recipes", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ recipe: record.recipe, trust_python: true, replace: true }),
    });
    if (!response.ok) {
      const detail = await response.json().catch(() => null);
      $<HTMLDivElement>("recipeStatus").textContent = detail?.detail ?? "Trust failed";
      return;
    }
    await loadRecipes();
    $<HTMLSelectElement>("recipeSel").value = selectedRecipeId;
    showRecipeRecord();
    $<HTMLDivElement>("recipeStatus").textContent = "Python trusted on this computer.";
  });
  $<HTMLButtonElement>("recipeExport").addEventListener("click", () => {
    const record = recipeRecords.find(({ recipe }) => recipe.id === selectedRecipeId);
    if (record) download(`${record.recipe.id}.oqp-recipe.json`, JSON.stringify(record.recipe, null, 2) + "\n");
  });
  $<HTMLButtonElement>("recipeDelete").addEventListener("click", async () => {
    const record = recipeRecords.find(({ recipe }) => recipe.id === selectedRecipeId);
    if (!record || record.builtin || !window.confirm(`Delete recipe “${record.recipe.name}”?`)) return;
    const response = await fetch(`/api/recipes/${encodeURIComponent(record.recipe.id)}`, { method: "DELETE" });
    $<HTMLDivElement>("recipeStatus").textContent = response.ok ? `Deleted ${record.recipe.name}` : "Delete failed";
    selectedRecipeId = "";
    await loadRecipes();
  });
  $<HTMLButtonElement>("recipeImport").addEventListener("click", () => $<HTMLInputElement>("recipeFile").click());
  $<HTMLInputElement>("recipeFile").addEventListener("change", async (event) => {
    const file = (event.target as HTMLInputElement).files?.[0];
    if (!file) return;
    try {
      const recipe = JSON.parse(await file.text());
      const response = await fetch("/api/recipes", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ recipe, trust_python: false, replace: false }),
      });
      if (!response.ok) throw new Error((await response.json()).detail ?? `import failed (${response.status})`);
      selectedRecipeId = recipe.id;
      await loadRecipes();
      $<HTMLSelectElement>("recipeSel").value = recipe.id;
      showRecipeRecord();
      $<HTMLDivElement>("recipeStatus").textContent =
        `Imported ${recipe.name}; embedded Python remains disabled until you review and trust it.`;
    } catch (error) {
      $<HTMLDivElement>("recipeStatus").textContent = `Import failed: ${error}`;
    }
  });
}

// The Basis Set Exchange, which the engine reads its basis sets from, knows
// the polarised Pople sets by their asterisk spelling: "6-31G*" is a name it
// has, "6-31G(d)" is not -- and asking for a name it does not have aborts the
// job at basis setup with a KeyError, after the settings dump, which reads
// like the calculation failed rather than like a typo.  Both spellings mean
// the same basis, so translate rather than refuse.
function normalizeBasis(name: string): string {
  return name.trim().replace(
    /^(\d-\d+\+*g)\((d|d,p)\)$/i,
    (_, stem: string, pol: string) => stem + (pol.toLowerCase() === "d" ? "*" : "**"),
  );
}

function currentBasis(): string {
  return normalizeBasis(
    basisSel.value === "__custom__" ? basisCustom.value : basisSel.value,
  );
}

function currentFunctional(): string {
  return functionalSel.value === "__custom__"
    ? functionalCustom.value.trim() || "bhhlyp"
    : functionalSel.value;
}

function fieldValue(id: string): string {
  return $<HTMLInputElement | HTMLSelectElement>(id).value.trim();
}

function optionList(entries: [string, string][], defaults: Record<string, string> = {}): string {
  return entries
    .filter(([key, value]) => value !== "" && value !== defaults[key])
    .map(([key, value]) => `${key}=${value}`)
    .join(",");
}

const SCF_DEFAULTS: Record<string, string> = {
  maxit: "30", conv: "1e-6", maxdiis: "7", diis_type: "cdiis", vshift: "0",
  converger_type: "diis", alternative_scf: "trah", escalation: "diis,soscf,trah",
  stability: "false", soscf_lvl_shift: "0", mom: "false", mom_switch: "0.003",
  pfon: "false", pfon_start_temp: "2000", pfon_cooling_rate: "50",
  trh_stab: "false", trh_ls: "false", trh_sub_solver: "davidson", trh_r0: "0.4",
  trh_nmic: "50", trh_nrtv: "1", trh_jd_start: "30", trh_gred: "0.001",
  trh_lred: "0.0001", trh_impl: "auto", init_scf: "no", init_it: "15",
  init_conv: "0.001", rstctmo: "false", incremental: "true", pscreen: "false",
  pscreen_k: "1e-2", pscreen_cap: "1e-8", pscreen_tight: "1e-4", verbose: "1",
};

const RESPONSE_DEFAULTS: Record<string, string> = {
  maxit: "50", conv: "1e-6", maxit_zv: "50", zvconv: "1e-6", nvdav: "50",
  z_solver: "0", gmres_dim: "50", resp_cutoff: "auto",
};

const HESS_DEFAULTS: Record<string, string> = {
  dx: "0.01", nproc: "1", temperature: "298.15", symmetry_unique: "false",
};

const NAC_DEFAULTS: Record<string, string> = {
  type: "numerical", dx: "1e-4", nproc: "1",
};

const OPTIMISATION_DEFAULTS: Record<string, string> = {
  maxit: "100", rmsd_grad: "0.003", max_grad: "0.004",
  rmsd_step: "0.004", max_step: "0.008",
};

const GEOMETRY_DEFAULTS: Record<string, string> = {
  coordsys: "dlc", trust: "0.1", trust_max: "0.3",
};

const CAS_DEFAULTS: Record<string, string> = {
  active_electrons: "0", active_orbitals: "0", frozen_core: "0",
};

const CI_DEFAULTS: Record<string, string> = {
  nroot: "1", solver: "auto", eig_tol: "1e-10", davidson_maxiter: "100",
};

const PCM_DEFAULTS: Record<string, string> = {
  model: "ddpcm", epsilon: "78.3553", radii: "uff",
};

const PCM_SOLVENT_DIELECTRICS: Record<string, string> = {
  water: "78.3553",
  acetonitrile: "35.688",
  methanol: "32.613",
  ethanol: "24.852",
  dimethylsulfoxide: "46.826",
  dimethylformamide: "36.71",
  acetone: "20.493",
  tetrahydrofuran: "7.4257",
  dichloromethane: "8.93",
  chloroform: "4.7113",
  "diethyl-ether": "4.24",
  toluene: "2.3741",
  benzene: "2.2706",
  dioxane: "2.2099",
  "ethyl-acetate": "6.02",
};

function withOptions(driver: string, options: string): string {
  if (!options) return driver;
  const opening = driver.indexOf("(");
  return opening >= 0 ? `${driver.slice(0, -1)},${options})` : `${driver}(${options})`;
}

function optimisationOptions(): string {
  // Concise .oqp geometry drivers select OpenQP's native optimizer automatically.
  // ``lib`` belongs only to the legacy [optimize] section and is rejected here.
  const options: [string, string][] = [
    ["maxit", fieldValue("geomMaxit")],
    ["rmsd_grad", fieldValue("rmsdGrad")],
    ["max_grad", fieldValue("maxGrad")],
    ["rmsd_step", fieldValue("rmsdStep")],
    ["max_step", fieldValue("maxStep")],
  ];
  return optionList(options, OPTIMISATION_DEFAULTS);
}

function geometryOptions(): string {
  // The native geometry engine is automatic.  Its controls belong to the
  // geometry driver itself; there is no separate user-facing optimizer choice.
  return optionList([
    ["coordsys", fieldValue("geomCoordsys") || "dlc"],
    ["trust", fieldValue("geomTrust")],
    ["trust_max", fieldValue("geomTrustMax")],
  ], GEOMETRY_DEFAULTS);
}

function hessianOptions(): string {
  return optionList([
    ["dx", fieldValue("hessDx")], ["nproc", fieldValue("hessNproc")],
    ["temperature", fieldValue("hessTemperature")], ["symmetry_unique", fieldValue("hessSymmetry")],
  ], HESS_DEFAULTS);
}

function nacOptions(): string {
  return optionList([
    ["type", fieldValue("nacType")], ["dx", fieldValue("nacDx")],
    ["nproc", fieldValue("nacNproc")],
  ], NAC_DEFAULTS);
}

function pcmModifier(): string {
  const solvent = fieldValue("pcmSolvent");
  const customSolvent = fieldValue("pcmSolventCustom") || "custom";
  const solventArgument = solvent === "__custom__"
    ? `solvent=${JSON.stringify(customSolvent)}`
    : solvent;
  const options = optionList([
    ["model", fieldValue("pcmModel")],
    ["epsilon", fieldValue("pcmEpsilon")],
    ["radii", fieldValue("pcmRadii")],
  ], PCM_DEFAULTS);
  return `pcm(${solventArgument}${options ? `,${options}` : ""})`;
}

function syncPcmSolvent(): void {
  const custom = fieldValue("pcmSolvent") === "__custom__";
  const customWrap = $<HTMLElement>("pcmCustomSolventWrap");
  customWrap.style.display = custom ? "" : "none";
  if (!custom) {
    $<HTMLInputElement>("pcmEpsilon").value = PCM_SOLVENT_DIELECTRICS[fieldValue("pcmSolvent")];
  }
}

function syncPcmReference(): void {
  const rohf = $<HTMLSelectElement>("pcmReference").querySelector<HTMLOptionElement>('option[value="rohf"]');
  const dftPcm = currentWf.key === "pcm" && theorySel.value === "dft";
  if (rohf) rohf.disabled = dftPcm;
  if (dftPcm) $<HTMLSelectElement>("pcmReference").value = "rhf";
}

// Generates OpenQP canonical .oqp input: route, one primary driver, globals, geometry.
function generateInp(): string {
  const atoms = parseAtoms(xyzArea.value);
  const theory = theorySel.value;
  const charge = +$<HTMLInputElement>("charge").value || 0;
  const mult = +$<HTMLInputElement>("mult").value || 1;
  const nstate = +$<HTMLInputElement>("nstate").value || 3;
  const target = +$<HTMLInputElement>("targetState").value || 0;
  const basis = currentBasis();
  const functional = currentFunctional();

  const inputError = inputValidation();
  if (inputError) return `# ${inputError}\n`;

  const routes: Record<string, string> = {
    hf: `hf/${basis}`,
    dft: `dft/${functional}/${basis}`,
    mp2: `mp2/${basis}`,
    tddft: `tddft(nstate=${nstate})/${functional}/${basis}`,
    tda: `tda(nstate=${nstate})/${functional}/${basis}`,
    sf: `sf(nstate=${nstate})/${functional}/${basis}`,
    mrsf: `mrsf(nstate=${nstate})/${functional}/${basis}`,
    umrsf: `umrsf(nstate=${nstate})/${functional}/${basis}`,
    ccsd: `ccsd/${basis}`,
    ccsd_t: `ccsd_t/${basis}`,
    fci: `fci/${basis}`,
    casci: `casci/${basis}`,
    casscf: `casscf/${basis}`,
    "sa-casscf": `sa-casscf/${basis}`,
    caspt2: `caspt2/${basis}`,
    "ms-caspt2": `ms-caspt2/${basis}`,
    "xms-caspt2": `xms-caspt2/${basis}`,
    nevpt2: `nevpt2/${basis}`,
    "sc-nevpt2": `sc-nevpt2/${basis}`,
    mrmp2: `mrmp2/${basis}`,
    mcqdpt2: `mcqdpt2/${basis}`,
    xmcqdpt2: `xmcqdpt2/${basis}`,
  };

  const usesStates = DFT_RESPONSE.includes(theory);
  const stateArg = usesStates && target > 0 ? `(S${target})` : "";
  let driver: string;
  switch (currentWf.key) {
    case "energy": driver = `energy${stateArg}`; break;
    case "grad": driver = `grad${stateArg}`; break;
    case "opt": driver = usesStates ? `opt(S${target})` : "opt"; break;
    case "scan": {
      const state = usesStates && target > 0 ? `S${target}` : "";
      driver = fieldValue("scanMode") === "relaxed"
        ? `opt(${[state, `freeze=distance(${fieldValue("scanAtomA")},${fieldValue("scanAtomB")})`].filter(Boolean).join(",")})`
        : `energy${state ? `(${state})` : ""}`;
      break;
    }
    case "hess": {
      const hessState = stateArg ? `${stateArg.slice(0, -1)},` : "(";
      driver = withOptions(`hess${hessState}type=${hessTypeSel.value})`, hessianOptions());
      break;
    }
    case "abs": driver = "energy"; break;
    case "exgrad": driver = `grad(S${target || 1})`; break;
    case "exopt": driver = `opt(S${target || 1})`; break;
    case "meci": {
      const a = +$<HTMLInputElement>("meciA").value || 0;
      const b = +$<HTMLInputElement>("meciB").value || 1;
      driver = `meci(S${Math.min(a, b)},S${Math.max(a, b)})`;
      break;
    }
    case "mecp": driver = `mecp(${fieldValue("mecpA") || "S0"},${fieldValue("mecpB") || "T0"})`; break;
    case "tci": driver = `tci(${fieldValue("tciA") || "S0"},${fieldValue("tciB") || "S1"},${fieldValue("tciC") || "S2"})`; break;
    case "ts": driver = `ts${stateArg}`; break;
    case "irc": driver = `irc${stateArg}`; break;
    case "mep": driver = `mep(S${target || 1})`; break;
    case "neb": {
      const product = fieldValue("nebProduct");
      driver = `neb(S${target || 0},${optionList([["maxit", fieldValue("geomMaxit")], ["product", product ? `\"${product}\"` : ""], ["nimage", fieldValue("nebImages")], ["spring", fieldValue("nebSpring")]])})`;
      break;
    }
    case "prop": driver = `prop${stateArg}`; break;
    case "nmr": {
      // ACID reads the GIAO magnetic response, so asking for it selects GIAO
      // rather than emitting a combination OpenQP will refuse.  The grid
      // keywords only travel with acid=true, which is the same rule the .oqp
      // surface enforces.
      const acid = $<HTMLInputElement>("nmrAcid").checked;
      const gauge = acid ? "giao" : fieldValue("nmrGauge");
      const options = [`gauge=${gauge}`];
      if (acid) {
        options.push("acid=true");
        const spacing = fieldValue("acidSpacing");
        const padding = fieldValue("acidPadding");
        if (spacing) options.push(`acid_spacing=${spacing}`);
        if (padding) options.push(`acid_padding=${padding}`);
      }
      driver = `nmr(${options.join(",")})`;
      break;
    }
    case "nac": driver = withOptions(`nac(${fieldValue("couplingA") || "S0"},${fieldValue("couplingB") || "S1"})`, nacOptions()); break;
    case "nacme": driver = `nacme(${fieldValue("couplingA") || "S0"},${fieldValue("couplingB") || "S1"})`; break;
    case "soc": driver = "soc"; break;
    case "ekt": {
      const opts = [];
      if ($<HTMLInputElement>("ektIp").checked) opts.push("ip=true");
      if ($<HTMLInputElement>("ektEa").checked) opts.push("ea=true");
      driver = `ekt(${opts.join(",") || "ip=true"})`;
      break;
    }
    case "namd": driver = `namd(S${fieldValue("namdState") || "1"},nstep=${fieldValue("namdSteps") || "100"},dt=${fieldValue("namdDt") || "0.5"},decoherence=${fieldValue("namdDecoherence") || "edc"})`; break;
    default: driver = "energy";
  }

  const route = currentWf.key === "pcm" && theory === "hf" && fieldValue("pcmReference") === "rohf"
    ? `rohf/${basis}`
    : routes[theory];
  const lines: string[] = [pdbSource ? `${route} qmmm_flag=true` : route, driver];
  const scfEntries: [string, string][] = [
    ["maxit", fieldValue("scfMaxit")], ["conv", fieldValue("scfConv")],
    ["maxdiis", fieldValue("scfMaxDiis")], ["diis_type", fieldValue("scfDiisType")],
    ["vshift", fieldValue("scfLevelShift")],
    ["converger_type", fieldValue("scfConverger")], ["alternative_scf", fieldValue("scfFallback")],
    ["escalation", fieldValue("scfEscalation")], ["stability", fieldValue("scfStability")],
    ["soscf_lvl_shift", fieldValue("soscfLevelShift")], ["mom", fieldValue("scfMom")],
    ["mom_switch", fieldValue("scfMomSwitch")], ["pfon", fieldValue("scfPfon")],
    ["pfon_start_temp", fieldValue("scfPfonTemp")], ["pfon_cooling_rate", fieldValue("scfPfonCooling")],
    ["trh_stab", fieldValue("trahStability")], ["trh_ls", fieldValue("trahLineSearch")],
    ["trh_sub_solver", fieldValue("trahSolver")], ["trh_r0", fieldValue("trahRadius")],
    ["trh_nmic", fieldValue("trahMicro")], ["trh_nrtv", fieldValue("trahTrialVectors")],
    ["trh_jd_start", fieldValue("trahJdStart")], ["trh_gred", fieldValue("trahGlobalReduction")],
    ["trh_lred", fieldValue("trahLocalReduction")], ["trh_impl", fieldValue("trahImplementation")],
    ["init_scf", fieldValue("scfInitial")], ["init_it", fieldValue("scfInitialIt")],
    ["init_conv", fieldValue("scfInitialConv")], ["rstctmo", fieldValue("scfRestart")],
    ["incremental", fieldValue("scfIncremental")], ["pscreen", fieldValue("scfPrescreen")],
    ["pscreen_k", fieldValue("scfPrescreenK")], ["pscreen_cap", fieldValue("scfPrescreenCap")],
    ["pscreen_tight", fieldValue("scfPrescreenTight")], ["verbose", fieldValue("scfVerbose")],
  ];
  const scfOptions = optionList(scfEntries, SCF_DEFAULTS);
  if (scfOptions) lines.push(`scf(${scfOptions})`);
  if (currentWf.key === "pcm") lines.push(pcmModifier());
  if (usesStates) {
    const responseOptions = optionList([
      ["maxit", fieldValue("responseMaxit")], ["conv", fieldValue("responseConv")],
      ["maxit_zv", fieldValue("zvectorMaxit")], ["zvconv", fieldValue("zvectorConv")],
      ["nvdav", fieldValue("davidsonSubspace")], ["z_solver", fieldValue("zvectorSolver")],
      ["gmres_dim", fieldValue("gmresDim")], ["resp_cutoff", fieldValue("responseCutoff")],
    ], RESPONSE_DEFAULTS);
    if (responseOptions) lines.push(`tdhf(${responseOptions})`);
  }
  if (theory === "mrsf") lines.push("guess(save_mol=true)");
  if (["fci", "casci", "casscf", "sa-casscf", "caspt2", "ms-caspt2", "xms-caspt2", "nevpt2", "sc-nevpt2", "mrmp2", "mcqdpt2", "xmcqdpt2"].includes(theory)) {
    const casOptions = optionList([
      ["active_electrons", fieldValue("activeElectrons")], ["active_orbitals", fieldValue("activeOrbitals")],
      ["frozen_core", fieldValue("frozenCore")],
    ], CAS_DEFAULTS);
    if (casOptions) lines.push(`cas(${casOptions})`);
    const ciOptions = optionList([
      ["nroot", fieldValue("ciRoots")], ["solver", fieldValue("ciSolver")],
      ["eig_tol", fieldValue("ciTolerance")], ["davidson_maxiter", fieldValue("ciMaxit")],
    ], CI_DEFAULTS);
    if (ciOptions) lines.push(`ci(${ciOptions})`);
  }
  if (["opt", "exopt", "ts"].includes(currentWf.key) ||
      (currentWf.key === "scan" && fieldValue("scanMode") === "relaxed")) {
    const opts = [optimisationOptions(), geometryOptions()].filter(Boolean).join(",");
    lines[1] = withOptions(driver, opts);
  }
  if (["meci", "mecp", "tci"].includes(currentWf.key)) {
    const crossing = optionList([
      ["maxit", fieldValue("geomMaxit")], ["energy_shift", fieldValue("energyShift")],
      ["energy_gap", fieldValue("energyGap")], ["rmsd_grad", fieldValue("rmsdGrad")],
      ["max_grad", fieldValue("maxGrad")], ["rmsd_step", fieldValue("rmsdStep")],
      ["max_step", fieldValue("maxStep")],
    ]);
    const algorithm = fieldValue("crossingAlgorithm");
    const algorithmOption = currentWf.key === "meci"
      ? `algorithm=${algorithm}`
      : currentWf.key === "mecp" ? `mecp_search=${algorithm}` : "";
    const driverOptions = [algorithmOption, crossing, geometryOptions()].filter(Boolean).join(",");
    lines[1] = withOptions(driver, driverOptions);
  }
  if (charge !== 0) lines.push(`charge=${charge}`);
  if (currentWf.key === "nacme") {
    const previousGeometry = fieldValue("nacmeGeometry");
    if (previousGeometry) lines.push(`geom2=${JSON.stringify(previousGeometry)}`);
  }
  // MRSF selects its high-spin working reference automatically — no mult.
  if (mult !== 1 && theory !== "mrsf") lines.push(`mult=${mult}`);
  if (pdbSource) {
    const selected = [...qmAtoms.keys()].sort((a, b) => a - b);
    const forcefield = fieldValue("qmmmForcefield") || "amber14-all.xml";
    lines.push(`qmmm(pdb_file="${pdbSource.name}",forcefield_files="${forcefield}",qm_atoms="${selected.join(" ")}")`);
    lines.push(`geom="${pdbSource.name} ${selected.join(" ")}"`);
  } else {
    lines.push('geom="""');
    for (const [el, x, y, z] of atoms) {
      lines.push(`${el.padEnd(2)} ${x.toFixed(6).padStart(11)} ${y.toFixed(6).padStart(11)} ${z.toFixed(6).padStart(11)}`);
    }
    lines.push('"""');
  }
  return lines.join("\n") + "\n";
}

const inpPreview = $<HTMLPreElement>("inpPreview");
const inpPreviewEditor = $<HTMLTextAreaElement>("inpPreviewEditor");
const editPreviewButton = $<HTMLButtonElement>("editPreview");
let previewIsEditing = false;

function previewInput(): string {
  return previewIsEditing ? inpPreviewEditor.value : generateInp();
}

function updateInpPreview(): void {
  if (!previewIsEditing) inpPreview.textContent = generateInp();
}

editPreviewButton.addEventListener("click", () => {
  previewIsEditing = !previewIsEditing;
  if (previewIsEditing) {
    inpPreviewEditor.value = inpPreview.textContent ?? generateInp();
    inpPreview.style.display = "none";
    inpPreviewEditor.style.display = "block";
    inpPreviewEditor.focus();
    editPreviewButton.textContent = "Done editing";
  } else {
    inpPreview.textContent = inpPreviewEditor.value;
    inpPreview.style.display = "";
    inpPreviewEditor.style.display = "none";
    editPreviewButton.textContent = "Edit input";
  }
});
function syncFieldStates(): void {
  const usesFunctional = ["dft", ...DFT_RESPONSE].includes(theorySel.value);
  functionalSel.disabled = !usesFunctional;
  functionalCustom.disabled = !usesFunctional || functionalSel.value !== "__custom__";
  basisCustom.disabled = basisSel.value !== "__custom__";

  const analyticHessian = !theorySel.disabled &&
    (theorySel.value === "hf" || theorySel.value === "dft");
  if (!analyticHessian) hessTypeSel.value = "numerical";
  hessTypeSel.disabled = !analyticHessian;
  hessTypeHint.textContent = analyticHessian
    ? "Analytical is available for the ground state; numerical remains selectable."
    : "Analytical Hessian is unavailable for this theory; numerical is required.";
}

for (const id of ["theory", "functional", "functionalCustom", "basis", "basisCustom", "charge", "mult",
                  "nstate", "targetState", "nmrGauge", "nmrAcid", "acidSpacing", "acidPadding", "meciA", "meciB", "mecpA", "mecpB", "tciA", "tciB", "tciC",
                  "couplingA", "couplingB", "nacmeGeometry", "nebProduct", "nebImages", "nebSpring", "pcmReference", "pcmSolvent", "pcmSolventCustom", "pcmModel",
                  "pcmEpsilon", "pcmRadii", "ektIp", "ektEa", "hessType", "scfMaxit", "scfConv", "responseMaxit",
                  "responseCutoff", "responseConv", "zvectorMaxit", "zvectorConv", "davidsonSubspace", "zvectorSolver", "gmresDim",
                  "scfMaxDiis", "scfDiisType", "scfLevelShift", "scfConverger", "scfFallback", "scfEscalation", "scfStability", "soscfLevelShift",
                  "scfMom", "scfMomSwitch", "scfPfon", "scfPfonTemp", "scfPfonCooling", "trahStability", "trahLineSearch", "trahSolver", "trahRadius",
                  "trahMicro", "trahTrialVectors", "trahJdStart", "trahGlobalReduction", "trahLocalReduction", "trahImplementation", "scfInitial", "scfInitialIt",
                  "scfInitialConv", "scfRestart", "scfIncremental", "scfPrescreen", "scfPrescreenK", "scfPrescreenCap", "scfPrescreenTight", "scfVerbose",
                  "geomCoordsys", "geomTrust", "geomTrustMax",
                  "hessDx", "hessNproc", "hessTemperature", "hessSymmetry", "nacType", "nacDx", "nacNproc",
                  "namdState", "namdSteps", "namdDt", "namdDecoherence", "activeElectrons", "activeOrbitals", "frozenCore", "ciRoots", "ciSolver", "ciTolerance", "ciMaxit",
                  "geomMaxit", "rmsdGrad", "maxGrad", "rmsdStep", "maxStep", "energyShift",
                  "energyGap", "trustRadius", "crossingAlgorithm", "scanAtomA", "scanAtomB",
                  "scanMode", "scanStart", "scanEnd", "scanPoints"]) {
  $<HTMLElement>(id).addEventListener("input", () => {
    syncFieldStates();
    if (id === "pcmSolvent") syncPcmSolvent();
    if (id === "theory") {
      syncPcmReference();
      syncWorkflowOptions();
    }
    syncWorkflowDetails();
    updateInpPreview();
  });
}

$<HTMLButtonElement>("generate").addEventListener("click", () => {
  const inputText = previewInput();
  const inputError = inputValidation(inputText);
  if (inputError) {
    builderStatus.textContent = inputError;
    showTab("builder");
    return;
  }
  $<HTMLTextAreaElement>("input").value = inputText;
  const name = projectName.value.trim() || "input";
  $<HTMLInputElement>("inputName").value = `${name}.oqp`;
  showTab("run");
});

// ---------- run ----------
const runStatus = $<HTMLSpanElement>("runStatus");
const runLog = $<HTMLPreElement>("runLog");
const runnerSelect = $<HTMLSelectElement>("runner");
const runButton = $<HTMLButtonElement>("runBtn");
const threadsInput = $<HTMLInputElement>("threads");
let admissionPermitted = true;
let submitting = false;
let admissionTimer: number | undefined;
let threadsInitialized = false;
const RUNNER_LABELS: Record<string, string> = {
  local: "OpenQP (local)",
  bundled: "OpenQP (bundled)",
  wsl: "OpenQP (WSL)",
};

function runnerLabel(name: string, version?: string): string {
  const label = RUNNER_LABELS[name] ?? name;
  return version ? `${label} v${version}` : label;
}

type HostStatus = {
  platform: string;
  physical_cores: number;
  logical_cores: number;
  memory_total_bytes: number | null;
  memory_available_bytes: number | null;
  estimated_memory_bytes?: number;
  permitted?: boolean;
  reason?: string | null;
};

function gib(bytes: number | null | undefined): string {
  return bytes == null ? "unavailable" : `${(bytes / 1024 ** 3).toFixed(1)} GiB`;
}

function threadCount(): number {
  return Math.max(1, Math.floor(+threadsInput.value || 1));
}

function renderHostStatus(host: HostStatus): void {
  const physical = host.physical_cores;
  threadsInput.max = String(host.logical_cores);
  if (!threadsInitialized) {
    threadsInput.value = String(physical);
    threadsInitialized = true;
  }
  const info = $<HTMLDivElement>("hostInfo");
  info.replaceChildren();
  const cores = document.createElement("div");
  cores.textContent = `${physical} physical cores · ${host.logical_cores} logical cores`;
  const memory = document.createElement("div");
  memory.textContent = `RAM: ${gib(host.memory_total_bytes)} installed · ${gib(host.memory_available_bytes)} available now`;
  info.append(cores, memory);
}

async function checkMemoryAdmission(): Promise<void> {
  const input = $<HTMLTextAreaElement>("input").value;
  if (!input.trim()) {
    $<HTMLDivElement>("memoryInfo").textContent = "Waiting for input…";
    admissionPermitted = true;
    return;
  }
  try {
    const response = await fetch("/api/host/admission", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ input_text: input, threads: threadCount() }),
    });
    const check: HostStatus = await response.json();
    admissionPermitted = check.permitted !== false;
    $<HTMLDivElement>("memoryInfo").textContent = admissionPermitted
      ? `Estimated ${gib(check.estimated_memory_bytes)}; admitted`
      : `RAM limit: estimated ${gib(check.estimated_memory_bytes)}, available ${gib(check.memory_available_bytes)}`;
    $<HTMLDivElement>("memoryInfo").style.color = admissionPermitted ? "" : "var(--err)";
  } catch {
    admissionPermitted = true;
    $<HTMLDivElement>("memoryInfo").textContent = "Memory preflight unavailable; the server will check before running.";
  }
}

function scheduleMemoryAdmission(): void {
  if (admissionTimer !== undefined) window.clearTimeout(admissionTimer);
  admissionTimer = window.setTimeout(() => { void checkMemoryAdmission(); }, 250);
}

async function refreshExecutionHost(): Promise<void> {
  try {
    const host: HostStatus = await (await fetch("/api/host")).json();
    renderHostStatus(host);
    await checkMemoryAdmission();
  } catch {
    $<HTMLDivElement>("hostInfo").textContent = "Hardware status unavailable.";
  }
}

async function loadRunners(): Promise<void> {
  const info = $<HTMLSpanElement>("runnersInfo");
  const stop = busyStatus(info, "Finding OpenQP engines");
  runnerSelect.disabled = true;
  runButton.disabled = true;
  try {
    const response = await fetch("/api/runners/detail");
    if (!response.ok) throw new Error(`engine discovery failed (${response.status})`);
    const detail = await response.json();
    const runners: Record<string, boolean> = detail.available;
    const versions: Record<string, string | undefined> = detail.versions ?? {};
    runnerSelect.innerHTML = "";
    const prompt = new Option("Choose OpenQP runner…", "", true, true);
    prompt.disabled = true;
    runnerSelect.appendChild(prompt);
    for (const [name, available] of Object.entries(runners)) {
      const label = runnerLabel(name, versions[name]);
      const option = new Option(available ? label : `${label} (unavailable)`, name);
      option.disabled = !available;
      runnerSelect.appendChild(option);
    }
    stop();
    info.textContent = "runners: " + Object.entries(runners)
      .map(([name, ok]) => `${runnerLabel(name, versions[name])} ${ok ? "✓" : "✗"}`).join(" · ");
    info.title = [detail.openqp && `local: ${detail.openqp}`,
      detail.bundled_openqp && `bundled: ${detail.bundled_openqp}`].filter(Boolean).join("\n");
    runStatus.textContent = "Choose an OpenQP runner";
    runnerSelect.disabled = false;
  } catch (error) {
    stop();
    info.textContent = `Engine discovery unavailable: ${String(error)}`;
    runStatus.textContent = "Restart Studio or check the engine installation.";
  }
}

runnerSelect.addEventListener("change", () => {
  runButton.disabled = submitting || !runnerSelect.value;
  runStatus.textContent = runnerSelect.value ? "" : "Choose an OpenQP runner";
});

threadsInput.addEventListener("input", scheduleMemoryAdmission);
$<HTMLTextAreaElement>("input").addEventListener("input", scheduleMemoryAdmission);
$<HTMLButtonElement>("hostRefresh").addEventListener("click", () => { void refreshExecutionHost(); });
$<HTMLButtonElement>("executionWorkspace").addEventListener("click", () => { void openWorkspaceSettings(); });

async function pollJob(jobId: string): Promise<void> {
  const infoResponse = await fetch(`/api/jobs/${jobId}`);
  if (!infoResponse.ok) {
    const detail = await infoResponse.json().then((data) => data?.detail).catch(() => null);
    runStatus.textContent = detail ?? `job lookup failed (${infoResponse.status})`;
    runLog.textContent = "(job output is unavailable)";
    return;
  }
  const info = await infoResponse.json();
  const tailResponse = await fetch(`/api/jobs/${jobId}/log`);
  const tail = tailResponse.ok ? await tailResponse.json() : { log: "" };
  runLog.textContent = tail.log || "(no output yet)";
  runStatus.textContent = info.status + (info.error ? ` — ${info.error}` : "") +
    (info.postprocess_status ? ` · Python ${info.postprocess_status}` : "") +
    (info.postprocess_error ? ` — ${info.postprocess_error}` : "");
  const busy = ["queued", "running", "cancelling"].includes(info.status);
  runStatus.classList.toggle("busy", busy);
  runStatus.setAttribute("aria-busy", String(busy));
  if (busy) {
    setTimeout(() => pollJob(jobId), 1000);
  } else {
    if (info.status === "not_converged") {
      window.alert(`Geometry optimization did not converge.\n\n${info.error ?? "The best geometry was retained."}\n\nUse Restart from best geometry in Analysis to continue.`);
    }
    selectJob(jobId);
  }
}

async function pollScan(groupId: string): Promise<void> {
  const response = await fetch(`/api/scans/${groupId}`);
  if (!response.ok) {
    runStatus.textContent = `scan lookup failed (${response.status})`;
    return;
  }
  const data: { points: ScanPoint[] } = await response.json();
  const running = data.points.find((point) => point.status === "running");
  const done = data.points.filter((point) => point.status === "done").length;
  const terminal = new Set(["done", "failed", "not_converged", "cancelled"]);
  runStatus.textContent = running
    ? `scan ${done}/${data.points.length} complete · ${running.value.toFixed(4)} Å running`
    : `scan ${done}/${data.points.length} complete`;
  if (running) {
    const logResponse = await fetch(`/api/jobs/${running.job_id}/log`);
    const tail = logResponse.ok ? await logResponse.json() : { log: "" };
    runLog.textContent = tail.log || "(no output yet)";
  }
  await refreshJobs();
  if (data.points.some((point) => !terminal.has(point.status))) {
    window.setTimeout(() => { void pollScan(groupId); }, 1000);
    return;
  }
  const completed = [...data.points].reverse().find((point) => point.status === "done");
  showTab("analysis");
  if (completed) await selectJob(completed.job_id);
  else runStatus.textContent += " · no point completed successfully";
}

runButton.addEventListener("click", async () => {
  const inputError = inputValidation($<HTMLTextAreaElement>("input").value);
  if (inputError) {
    runStatus.textContent = inputError;
    return;
  }
  if (!runnerSelect.value) {
    runStatus.textContent = "Choose an OpenQP runner";
    return;
  }
  if (submitting || runButton.disabled) return;
  submitting = true;
  runButton.disabled = true;
  const stop = busyStatus(runStatus, "Preparing OpenQP");
  try {
  await checkMemoryAdmission();
  if (!admissionPermitted) {
    runStatus.textContent = "RAM limit: free memory or reduce the calculation before running.";
    return;
  }
  runStatus.textContent = "submitting…";
  const isScan = currentWf.key === "scan";
  const res = await fetch(isScan ? "/api/scans" : "/api/jobs", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      input_text: $<HTMLTextAreaElement>("input").value,
      input_name: $<HTMLInputElement>("inputName").value,
      name: projectName.value.trim() || "job",
      runner: runnerSelect.value,
      threads: threadCount(),
      ...(activeRecipeId ? { recipe_id: activeRecipeId } : {}),
      ...(isScan ? {
        atom_a: +fieldValue("scanAtomA"), atom_b: +fieldValue("scanAtomB"),
        start: +fieldValue("scanStart"), end: +fieldValue("scanEnd"),
        points: +fieldValue("scanPoints"), relaxed: fieldValue("scanMode") === "relaxed",
      } : {}),
      ...(pdbSource ? { pdb_text: pdbSource.text, pdb_name: pdbSource.name } : {}),
    }),
  });
  if (!res.ok) {
    // The backend says why — usually a job directory it cannot write to.
    const detail = await res.json().then((d) => d?.detail).catch(() => null);
    runStatus.textContent = detail ?? `submit failed (${res.status})`;
    return;
  }
  const submitted = await res.json();
  stop();
  if (isScan) void pollScan(submitted.group_id);
  else void pollJob(submitted.id);
  } catch (error) {
    runStatus.textContent = `Could not start OpenQP: ${String(error)}`;
  } finally {
    stop();
    submitting = false;
    runButton.disabled = !runnerSelect.value;
  }
});

// ---------- results ----------
let selectedJob = "";
let resultMoldenFiles: string[] = [];
let resultCubeFiles: string[] = [];
type ListedJob = { id: string; name: string; runner: string; status: string; error?: string | null };
let listedJobs: ListedJob[] = [];

type ScanPoint = {
  job_id: string; name: string; status: string;
  value: number; unit: string; energy: number | null;
};

async function refreshJobs(): Promise<void> {
  const jobs: ListedJob[] = await (await fetch("/api/jobs")).json();
  listedJobs = jobs;
  const body = $<HTMLTableSectionElement>("jobsBody");
  body.innerHTML = "";
  for (const job of jobs) {
    const tr = document.createElement("tr");
    if (job.id === selectedJob) tr.className = "sel";
    const name = document.createElement("td");
    name.textContent = job.name;
    const runner = document.createElement("td");
    runner.textContent = job.runner;
    const status = document.createElement("td");
    const badge = document.createElement("span");
    badge.classList.add("badge", job.status);
    badge.textContent = job.status;
    status.appendChild(badge);
    const actions = document.createElement("td");
    const action = (label: string, className: string, title = "") => {
      const button = document.createElement("button");
      button.className = `ghost ${className}`;
      button.type = "button";
      button.dataset.jobId = job.id;
      button.textContent = label;
      button.title = title;
      actions.append(button, " ");
      return button;
    };
    if (job.status === "not_converged") {
      action("Restart", "job-restart", "Start from the retained best geometry");
    }
    if (job.status === "running" || job.status === "queued") {
      action("Cancel", "job-cancel");
    }
    const running = ["running", "queued", "cancelling"].includes(job.status);
    const deleteButton = action("Delete", "job-delete",
      running ? "Project is running" : "Delete project and its files");
    deleteButton.disabled = running;
    deleteButton.setAttribute("aria-label", `Delete ${job.name}`);
    tr.append(name, runner, status, actions);
    tr.addEventListener("click", () => selectJob(job.id));
    tr.querySelector<HTMLButtonElement>(".job-delete")!.addEventListener("click", (event) => {
      event.stopPropagation();
      void deleteJob(job.id, job.name);
    });
    tr.querySelector<HTMLButtonElement>(".job-restart")?.addEventListener("click", (event) => {
      event.stopPropagation();
      void restartJob(job.id, job.name);
    });
    tr.querySelector<HTMLButtonElement>(".job-cancel")?.addEventListener("click", (event) => {
      event.stopPropagation();
      void cancelJob(job.id, job.name);
    });
    body.appendChild(tr);
  }
  syncComparisonProjects();
}
$<HTMLButtonElement>("jobsRefresh").addEventListener("click", refreshJobs);

function syncComparisonProjects(): void {
  const select = $<HTMLSelectElement>("comparisonReference");
  const previous = select.value;
  const current = listedJobs.find((job) => job.id === selectedJob);
  const currentIsTerminal = !!current && ["done", "not_converged"].includes(current.status);
  const candidates = listedJobs.filter((job) =>
    job.id !== selectedJob && ["done", "not_converged"].includes(job.status));
  select.innerHTML = '<option value="">Choose a reference…</option>' + candidates.map((job) =>
    `<option value="${escapeMarkup(job.id)}">${escapeMarkup(job.name)}</option>`).join("");
  if (candidates.some((job) => job.id === previous)) select.value = previous;
  $<HTMLDivElement>("comparisonCard").style.display = currentIsTerminal && candidates.length ? "" : "none";
  $<HTMLButtonElement>("comparisonRun").disabled = !currentIsTerminal || !select.value;
}

let comparisonRequestId = 0;
$<HTMLSelectElement>("comparisonReference").addEventListener("change", (event) => {
  comparisonRequestId += 1;
  $<HTMLButtonElement>("comparisonRun").disabled = !(event.target as HTMLSelectElement).value;
});

$<HTMLButtonElement>("comparisonRun").addEventListener("click", async () => {
  const left = $<HTMLSelectElement>("comparisonReference").value;
  const right = selectedJob;
  if (!left || !right) return;
  const requestId = ++comparisonRequestId;
  const isCurrent = () => requestId === comparisonRequestId && right === selectedJob &&
    left === $<HTMLSelectElement>("comparisonReference").value;
  const body = $<HTMLDivElement>("comparisonBody");
  body.innerHTML = '<div class="hint">Comparing projects…</div>';
  let response: Response;
  try {
    response = await fetch(
      `/api/comparison?left=${encodeURIComponent(left)}&right=${encodeURIComponent(right)}`,
    );
  } catch (error) {
    if (isCurrent()) {
      body.innerHTML = `<div class="hint">${escapeMarkup(
        `comparison failed: ${error instanceof Error ? error.message : String(error)}`,
      )}</div>`;
    }
    return;
  }
  if (!isCurrent()) return;
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    if (!isCurrent()) return;
    body.innerHTML = `<div class="hint">${escapeMarkup(detail?.detail ?? `comparison failed (${response.status})`)}</div>`;
    return;
  }
  const data = await response.json();
  if (!isCurrent()) return;
  const values: [string, string][] = [];
  if (data.left.energy != null) values.push([escapeMarkup(`${data.left.name} energy (Ha)`), fixed(data.left.energy, 8)]);
  if (data.right.energy != null) values.push([escapeMarkup(`${data.right.name} energy (Ha)`), fixed(data.right.energy, 8)]);
  if (data.energy_delta_hartree != null) values.push(["ΔE, current − reference (Ha)", fixed(data.energy_delta_hartree, 8)]);
  if (data.energy_delta_kcal_mol != null) values.push(["ΔE, current − reference (kcal/mol)", fixed(data.energy_delta_kcal_mol, 3)]);
  if (data.geometry.available) values.push(["Aligned structure RMSD (Å)", fixed(data.geometry.rmsd_angstrom, 5)]);
  else values.push(["Structure RMSD", data.geometry.reason]);
  if (data.dipole_delta_debye != null) values.push(["Δ|μ| (Debye)", fixed(data.dipole_delta_debye, 4)]);
  const stateRows = data.states.length ? `<div class="sum-title">Matched excited states</div><table class="sum">
    <tr><th>State</th><th>Reference (eV)</th><th>Current (eV)</th><th>Δ (eV)</th></tr>${data.states.map(
      (state: { index: number; left_ev: number; right_ev: number; delta_ev: number }) =>
        `<tr><td class="k">S${state.index}</td><td class="v">${fixed(state.left_ev, 4)}</td><td class="v">${fixed(state.right_ev, 4)}</td><td class="v">${fixed(state.delta_ev, 4)}</td></tr>`,
    ).join("")}</table>` : "";
  body.innerHTML = `<div class="sum-title">Current relative to reference</div>${rows(values)}${stateRows}`;
});

async function cancelJob(jobId: string, name: string): Promise<void> {
  if (!window.confirm(`Cancel the running calculation “${name}”?`)) return;
  const response = await fetch(`/api/jobs/${jobId}/cancel`, { method: "POST" });
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    importStatus.textContent = `cancel failed: ${detail?.detail ?? response.status}`;
    return;
  }
  importStatus.textContent = `cancelling ${name}`;
  await refreshJobs();
}

async function restartJob(jobId: string, name: string): Promise<void> {
  if (!window.confirm(`Prepare a restart for “${name}” from its retained best geometry?`)) return;
  const response = await fetch(`/api/jobs/${jobId}/restart`, { method: "POST" });
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    importStatus.textContent = `restart failed: ${detail?.detail ?? response.status}`;
    return;
  }
  const restarted: { input_text: string; input_name: string; name: string; runner: string | null; threads: number } = await response.json();
  $<HTMLTextAreaElement>("input").value = restarted.input_text;
  $<HTMLInputElement>("inputName").value = restarted.input_name;
  projectName.value = restarted.name;
  threadsInput.value = String(restarted.threads);
  if (restarted.runner && runnerSelect.querySelector(`option[value="${restarted.runner}"]`)) {
    runnerSelect.value = restarted.runner;
  }
  runButton.disabled = submitting || !runnerSelect.value;
  runStatus.textContent = runnerSelect.value
    ? "Restart input prepared from the retained best geometry."
    : "Choose an OpenQP runner to start the restart.";
  importStatus.textContent = `restart input prepared for ${name}`;
  showTab("run");
  scheduleMemoryAdmission();
}

async function deleteJob(jobId: string, name: string): Promise<void> {
  if (!window.confirm(`Delete project “${name}” and all of its files? This cannot be undone.`)) return;
  const response = await fetch(`/api/jobs/${jobId}`, { method: "DELETE" });
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    importStatus.textContent = `delete failed: ${detail?.detail ?? response.status}`;
    return;
  }
  if (selectedJob === jobId) {
    selectedJob = "";
    resetAnalysisProject();
    $<HTMLUListElement>("resultFiles").innerHTML = '<li class="hint">select a job</li>';
  }
  importStatus.textContent = `deleted ${name}`;
  await refreshJobs();
}

// Results computed elsewhere — on a cluster, by the standalone command-line
// engine, or in an earlier session — are copied into a job directory of their
// own. Analysis is built on top of a job directory, so from that point on they
// behave exactly like a run this app started: summary, spectra, orbitals,
// normal modes and property maps all apply.
const importFiles = $<HTMLInputElement>("importFiles");
const importStatus = $<HTMLSpanElement>("importStatus");

async function openExistingResults(): Promise<void> {
  importFiles.value = "";
  importFiles.click();
}

importFiles.addEventListener("change", async () => {
  const chosen = Array.from(importFiles.files ?? []);
  if (!chosen.length) return;
  importStatus.textContent = `importing ${chosen.length} file(s)…`;
  const body = new FormData();
  for (const file of chosen) body.append("files", file, file.name);
  body.append("name", chosen[0].name.replace(/\.[^.]+$/, ""));
  try {
    const res = await fetch("/api/jobs/import", { method: "POST", body });
    if (!res.ok) {
      const detail = await res.json().catch(() => null);
      importStatus.textContent = `import failed: ${detail?.detail ?? res.status}`;
      return;
    }
    const info = await res.json();
    importStatus.textContent = `imported as ${info.name}`;
    showTab("analysis");
    await refreshJobs();
    await selectJob(info.id);
  } catch (err) {
    importStatus.textContent = `import failed: ${err}`;
  }
});

$<HTMLButtonElement>("importBtn").addEventListener("click", () => {
  void openExistingResults();
});

async function selectJob(jobId: string): Promise<void> {
  selectedJob = jobId;
  resetAnalysisProject();
  refreshJobs();
  let scanGroup = "";
  const infoResponse = await fetch(`/api/jobs/${jobId}`);
  if (infoResponse.ok && jobId === selectedJob) {
    const info: { group_id?: string | null; recipe_id?: string | null } = await infoResponse.json();
    scanGroup = info.group_id ?? "";
    if (info.recipe_id) {
      activeRecipeId = info.recipe_id;
      const record = recipeRecords.find(({ recipe }) => recipe.id === info.recipe_id);
      if (record) {
        selectedRecipeId = record.recipe.id;
        pendingRecipeAnalysis = { ...record.recipe.analysis };
        artSettings = { ...artSettings, ...record.recipe.art };
      }
    } else {
      activeRecipeId = "";
      pendingRecipeAnalysis = {};
    }
    if (scanGroup) void loadScanGroup(scanGroup);
  }
  const summaryLoading = loadSummary(jobId).catch(() => {});
  const optimization = loadOptimizationHistory(jobId).catch(() => false);
  const files: { name: string; size: number }[] =
    await (await fetch(`/api/jobs/${jobId}/files`)).json();
  const hasOptimization = await optimization;
  if (jobId !== selectedJob) return;
  resultMoldenFiles = files
    .map((file) => file.name)
    .filter((name) => name.toLowerCase().endsWith(".molden"));
  resultCubeFiles = files
    .map((file) => file.name)
    .filter((name) => /\.(?:cube|cub)$/i.test(name));
  syncSurfaceFiles();
  loadExcitedAnalysis(jobId).catch(() => {});
  const list = $<HTMLUListElement>("resultFiles");
  list.innerHTML = files.length ? "" : '<li class="hint">no files</li>';
  for (const f of files) {
    const url = `/api/jobs/${jobId}/files/${encodeURIComponent(f.name)}`;
    const li = document.createElement("li");
    const viewable = VIEWABLE.some((ext) => f.name.toLowerCase().endsWith(ext));
    li.append(document.createTextNode(`${f.name} `));
    const size = document.createElement("span");
    size.className = "hint";
    size.textContent = `(${f.size.toLocaleString()} B)`;
    li.appendChild(size);
    if (viewable) {
      const view = document.createElement("a");
      view.href = "#";
      view.dataset.url = url;
      view.dataset.name = f.name;
      view.textContent = "view";
      li.append(" ", view);
    }
    const download = document.createElement("a");
    download.href = url;
    download.download = "";
    download.textContent = "download";
    li.append(" ", download);
    list.appendChild(li);
  }
  list.querySelectorAll<HTMLAnchorElement>("a[data-url]").forEach((a) => {
    a.addEventListener("click", (e) => {
      e.preventDefault();
      viewResultFile(jobId, a.dataset.name!, a.dataset.url!);
    });
  });

  // Show something straight away. Picking a job, or importing one, is the
  // user saying "let me see this" -- making them hunt for a "view" link
  // first is a step that never had a reason to exist.
  const names = files.map((f) => f.name);
  const best = scanGroup
    ? names.find((name) => /\.(?:oqp|inp)$/i.test(name)) ?? bestFileToShow(names)
    : bestFileToShow(names);
  if (best && !(scanGroup && hasOptimization)) {
    await viewResultFile(jobId, best,
                         `/api/jobs/${jobId}/files/${encodeURIComponent(best)}`);
  }
  await summaryLoading;
  if (jobId === selectedJob) applyPendingPropertyMap();
}

function resetAnalysisProject(): void {
  summary = null;
  lastSpectrum = null;
  resultMoldenFiles = [];
  resultCubeFiles = [];
  resultFrames = [];
  optimizationSteps = [];
  selectedScanGroup = "";
  pendingResultMessage = null;
  if (activeCubeUrl) {
    URL.revokeObjectURL(activeCubeUrl);
    activeCubeUrl = null;
  }
  releaseCurrentComponents();
  hideResultPanels();
  $<HTMLDivElement>("summaryCard").style.display = "none";
  $<HTMLDivElement>("summaryBody").innerHTML = "";
  $<HTMLDivElement>("spectrumCard").style.display = "none";
  excitedAnalysis = null;
  $<HTMLDivElement>("excitedAnalysisCard").style.display = "none";
  excitedSource.innerHTML = "";
  excitedTarget.innerHTML = "";
  excitedPair.innerHTML = "";
  $<HTMLDivElement>("excitedMetrics").textContent = "";
  specKind.innerHTML = "";
  specShape.value = "lorentzian";
  specWidth.value = "20";
  specState.value = "1";
  $<HTMLDivElement>("specStateWrap").style.display = "none";
  $<HTMLDivElement>("specPlot").innerHTML = "";
  $<HTMLDivElement>("specNote").textContent = "";
  orbitalSel.innerHTML = '<option value="">Choose an orbital…</option>';
  mapKind.value = "mo";
  $<HTMLOptionElement>("dysonMapKind").hidden = true;
  $<HTMLOptionElement>("dysonMapKind").disabled = true;
  modeSel.innerHTML = '<option value="">Choose a normal mode…</option>';
  modePlay.disabled = true;
  modePause.disabled = true;
  modeArrows.disabled = true;
  modeReset.disabled = true;
  $<HTMLDivElement>("resultFrameCard").style.display = "none";
  $<HTMLElement>("optimizationStepControls").style.display = "none";
  $<HTMLSelectElement>("optimizationStepSelect").innerHTML = "";
  $<HTMLDivElement>("optimizationStepValues").textContent = "";
  $<HTMLHeadingElement>("resultFrameTitle").textContent = "Trajectory";
  $<HTMLDivElement>("pathCard").style.display = "none";
  $<HTMLDivElement>("pathPlot").innerHTML = "";
  $<HTMLDivElement>("pathNote").textContent = "";
  $<HTMLDivElement>("comparisonBody").innerHTML = "";
  $<HTMLDivElement>("atomicPropertyCard").style.display = "none";
  $<HTMLSelectElement>("atomicProperty").innerHTML = "";
  $<HTMLDivElement>("propertyLegend").style.display = "none";
  $<HTMLDivElement>("propertyStatus").textContent = "";
  displayedResultAtomCount = 0;
  atomicProperties.clear();
  $<HTMLDivElement>("surfaceCard").style.display = "none";
  $<HTMLSelectElement>("surfacePrimary").innerHTML = "";
  $<HTMLSelectElement>("surfaceSecondary").innerHTML = "";
  $<HTMLSelectElement>("surfaceOperation").value = "display";
  $<HTMLDivElement>("surfaceSecondaryWrap").style.display = "none";
  $<HTMLSelectElement>("surfaceSides").value = "both";
  $<HTMLInputElement>("surfaceIso").value = "0.05";
  $<HTMLDivElement>("surfaceStatus").textContent = "";
  volumetricRequestId += 1;
  activeMapSource = null;
  activeDirectCube = null;
  activeSurfaceSelection = null;
  clearArtScene();
}

// What to open of its own accord, best first: a molden file carries orbitals
// and normal modes, a log or an optimisation trajectory carries the geometry
// as it moved, and the rest are still better than an empty panel.
const SHOW_FIRST = [".molden", ".trj", ".log", ".out", ".xyz", ".json", ".oqp", ".inp"];

function bestFileToShow(names: string[]): string | undefined {
  for (const extension of SHOW_FIRST) {
    const match = names.find((name) => name.toLowerCase().endsWith(extension));
    if (match) return match;
  }
  return undefined;
}

// Energies, states, frequencies and properties, read back from whatever the
// job wrote — the JSON export when there is one, the log otherwise.
async function loadSummary(jobId: string): Promise<void> {
  const response = await fetch(`/api/jobs/${jobId}/summary?refresh=true`);
  if (!response.ok) return;
  const data: Summary = await response.json();
  if (jobId !== selectedJob) return;
  summary = data;
  renderSummary(data);
  buildSpectrumList(data);
  if ($<HTMLDivElement>("spectrumCard").style.display === "") await loadSpectrum();
}

interface ExcitedAnalysis {
  available: boolean;
  reason?: string;
  source_state: number;
  target_state: number;
  transition_ev: number;
  states: { index: number; label: string; relative_ev: number }[];
  nto_pairs: { index: number; singular_value: number; weight: number; fraction: number }[];
  nto_participation_ratio: number;
  transition_density_norm: number;
  n_promoted: number;
  n_attach: number;
  n_detach: number;
  molden_file: string;
}

const excitedSource = $<HTMLSelectElement>("excitedSource");
const excitedTarget = $<HTMLSelectElement>("excitedTarget");
const excitedMap = $<HTMLSelectElement>("excitedMap");
const excitedPair = $<HTMLSelectElement>("excitedPair");
let excitedAnalysis: ExcitedAnalysis | null = null;
let excitedAnalysisRequestId = 0;

function stateOptions(states: ExcitedAnalysis["states"]): string {
  return states.map((state) =>
    `<option value="${state.index}">${state.label} (${state.relative_ev.toFixed(3)} eV rel. S0)</option>`
  ).join("");
}

function renderExcitedAnalysis(data: ExcitedAnalysis): void {
  excitedAnalysis = data;
  const card = $<HTMLDivElement>("excitedAnalysisCard");
  card.style.display = data.available ? "" : "none";
  if (!data.available) {
    pushArtSources();
    return;
  }
  if (!excitedSource.options.length) {
    const options = stateOptions(data.states);
    excitedSource.innerHTML = options;
    excitedTarget.innerHTML = options;
  }
  excitedSource.value = String(data.source_state);
  excitedTarget.value = String(data.target_state);
  excitedPair.innerHTML = data.nto_pairs.map((pair) =>
    `<option value="${pair.index}">Pair ${pair.index + 1} · ${(100 * pair.fraction).toFixed(2)}%</option>`
  ).join("");
  $<HTMLDivElement>("excitedPairWrap").style.display =
    excitedMap.value.startsWith("nto_") ? "" : "none";
  $<HTMLDivElement>("excitedMetrics").textContent = [
    `${data.transition_ev.toFixed(4)} eV`,
    `NTO participation ${data.nto_participation_ratio.toFixed(3)}`,
    `promoted charge ${data.n_promoted.toFixed(4)} e`,
  ].join(" · ");
  pushArtSources();
}

async function loadExcitedAnalysis(jobId: string): Promise<void> {
  const currentRef = () => excitedSource.value || "0";
  const currentTarget = () => excitedTarget.value || "1";
  const ref = currentRef();
  const target = currentTarget();
  if (ref === target) return;
  const requestId = ++excitedAnalysisRequestId;
  const response = await fetch(
    `/api/jobs/${jobId}/excited-analysis?ref=${encodeURIComponent(ref)}&target=${encodeURIComponent(target)}`
  );
  if (!response.ok || requestId !== excitedAnalysisRequestId || jobId !== selectedJob ||
      ref !== currentRef() || target !== currentTarget()) return;
  const data = await response.json();
  if (requestId !== excitedAnalysisRequestId || jobId !== selectedJob ||
      ref !== currentRef() || target !== currentTarget()) return;
  renderExcitedAnalysis(data);
}

async function showExcitedMap(): Promise<void> {
  if (!selectedJob || !excitedAnalysis?.available) return;
  const jobId = selectedJob;
  const requestId = ++volumetricRequestId;
  const kind = excitedMap.value;
  const ref = excitedSource.value;
  const target = excitedTarget.value;
  const rank = excitedPair.value || "0";
  const moldenFile = excitedAnalysis.molden_file;
  const selectionIsCurrent = () => jobId === selectedJob && kind === excitedMap.value &&
    ref === excitedSource.value && target === excitedTarget.value &&
    rank === (excitedPair.value || "0") && moldenFile === excitedAnalysis?.molden_file;
  const reportFailure = (message: string) => {
    if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
    $<HTMLDivElement>("excitedMetrics").textContent = message;
    pushArtSources();
  };
  const query = new URLSearchParams({
    kind,
    ref,
    target,
    rank,
  });
  const base = `/api/jobs/${jobId}`;
  let geometryResponse: Response;
  let cubeResponse: Response;
  try {
    [geometryResponse, cubeResponse] = await Promise.all([
      fetch(`${base}/molden/${encodeURIComponent(moldenFile)}/geom.xyz`),
      fetch(`${base}/excited-analysis/cube?${query}`),
    ]);
  } catch (error) {
    reportFailure(`map generation failed: ${error instanceof Error ? error.message : String(error)}`);
    return;
  }
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  if (!geometryResponse.ok || !cubeResponse.ok) {
    reportFailure(`map generation failed (${geometryResponse.status}/${cubeResponse.status})`);
    return;
  }
  const [xyz, cube] = await Promise.all([geometryResponse.text(), cubeResponse.text()]);
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  if (activeCubeUrl) URL.revokeObjectURL(activeCubeUrl);
  releaseCurrentComponents();
  activeCubeUrl = URL.createObjectURL(new Blob([cube], { type: "text/plain" }));
  activeMapSource = "excited";
  pushOrbitalStyle();
  pushToResultViewer({
    type: "oqp-cube",
    xyz,
    cube: activeCubeUrl,
    iso: kind.startsWith("nto_") ? 0.05 : 0.002,
    sides: $<HTMLSelectElement>("moSides").value,
  });
}

for (const select of [excitedSource, excitedTarget]) {
  select.addEventListener("change", () => {
    volumetricRequestId += 1;
    if (excitedSource.value === excitedTarget.value) {
      const states = excitedAnalysis?.states ?? [];
      const current = Number(excitedSource.value);
      const alternative = states.find((state) => state.index !== current);
      if (!alternative) return;
      excitedTarget.value = String(alternative.index);
    }
    void loadExcitedAnalysis(selectedJob);
  });
}
excitedMap.addEventListener("change", () => {
  volumetricRequestId += 1;
  $<HTMLDivElement>("excitedPairWrap").style.display =
    excitedMap.value.startsWith("nto_") ? "" : "none";
});
excitedPair.addEventListener("change", () => { volumetricRequestId += 1; });
$<HTMLButtonElement>("excitedShow").addEventListener("click", () => void showExcitedMap());
$<HTMLButtonElement>("excitedReset").addEventListener("click", () => {
  if (!excitedAnalysis) return;
  const requestId = ++volumetricRequestId;
  activeMapSource = null;
  if (activeCubeUrl) {
    URL.revokeObjectURL(activeCubeUrl);
    activeCubeUrl = null;
  }
  releaseCurrentComponents();
  fetch(`/api/jobs/${selectedJob}/molden/${encodeURIComponent(excitedAnalysis.molden_file)}/geom.xyz`)
    .then((response) => response.text())
    .then((xyz) => {
      if (requestId === volumetricRequestId) {
        pushToResultViewer({ type: "oqp-structure", xyz });
      }
    });
});

const resultFrame = $<HTMLIFrameElement>("resultFrame");
const orbitalCard = $<HTMLDivElement>("orbitalCard");
const orbitalSel = $<HTMLSelectElement>("orbitalSel");
const orbitalReset = $<HTMLButtonElement>("orbitalReset");
const isoRange = $<HTMLInputElement>("isoRange");
const modeCard = $<HTMLDivElement>("modeCard");
const modeSel = $<HTMLSelectElement>("modeSel");
const amplitudeRange = $<HTMLInputElement>("ampRange");
const modePlay = $<HTMLButtonElement>("modePlay");
const modePause = $<HTMLButtonElement>("modePause");
const modeReset = $<HTMLButtonElement>("modeReset");
const modeArrows = $<HTMLInputElement>("modeArrows");
let currentMolden: { jobId: string; name: string } | null = null;
type OrbitalSource = { name: string; index: number };
const orbitalSources = new Map<string, OrbitalSource>();
let activeCubeUrl: string | null = null;
// Jx/Jy/Jz blobs that go with an ACID scalar cube. Held separately from
// activeCubeUrl so they are revoked together with it and never outlive the
// surface they belong to.
let activeCurrentUrls: string[] = [];

// The arrow controls only mean anything while a current is on screen, so they
// are revealed with one and hidden again when the surface has none.
function showCurrentVectorControls(visible: boolean): void {
  $<HTMLDivElement>("currentVectorControls").style.display = visible ? "" : "none";
  if (!visible) $<HTMLDivElement>("currentVectorStatus").textContent = "";
}

function currentVectorSettings(): { current_density: number; current_scale: number } {
  return {
    current_density: Number($<HTMLInputElement>("currentVectorDensity").value) || 24,
    current_scale: Number($<HTMLInputElement>("currentVectorScale").value) || 1,
  };
}

function releaseCurrentComponents(): void {
  activeCurrentUrls.forEach((url) => URL.revokeObjectURL(url));
  activeCurrentUrls = [];
  showCurrentVectorControls(false);
}

// An ACID set is <stem>_acid.cube with <stem>_jx/_jy/_jz.cube beside it. The
// viewer refuses a mismatched grid anyway, but naming one dataset is what stops
// two calculations of the same molecule being combined in the first place.
function currentComponentNames(scalarName: string): { x: string; y: string; z: string } | null {
  const match = /^(.*)_acid\.(cube|cub)$/i.exec(scalarName);
  if (!match) return null;
  const [, stem, ext] = match;
  const names = { x: `${stem}_jx.${ext}`, y: `${stem}_jy.${ext}`, z: `${stem}_jz.${ext}` };
  const have = new Set(resultCubeFiles);
  return (have.has(names.x) && have.has(names.y) && have.has(names.z)) ? names : null;
}
// Shared by direct result-file opening and the surface panel's display mode.
// Callers release the previous surface first and recheck currency before posting.
async function loadCurrentComponents(
  jobId: string, url: string, isCurrent: () => boolean,
): Promise<Record<string, string> | undefined> {
  const componentNames = currentComponentNames(cubeNameFromUrl(url));
  showCurrentVectorControls(Boolean(componentNames));
  if (!componentNames || !$<HTMLInputElement>("currentVectorToggle").checked) return;
  try {
    const base = url.slice(0, url.lastIndexOf("/") + 1);
    const axes = ["x", "y", "z"] as const;
    const blobs = await Promise.all(axes.map(async (axis) => {
      const reply = await fetch(base + encodeURIComponent(componentNames[axis]));
      if (!reply.ok) throw new Error(`${componentNames[axis]} (${reply.status})`);
      return reply.blob();
    }));
    if (!isCurrent() || jobId !== selectedJob) return;
    const urls = blobs.map((blob) => URL.createObjectURL(blob));
    activeCurrentUrls = urls;
    return { x: urls[0], y: urls[1], z: urls[2] };
  } catch (error) {
    // An old failed request must not revoke a newer selection's blob URLs.
    if (!isCurrent() || jobId !== selectedJob) return;
    $<HTMLDivElement>("currentVectorStatus").textContent =
      `current-density components unavailable: ${error instanceof Error ? error.message : String(error)}`;
  }
}

let resultFrames: { label: string; atoms: Atom[] }[] = [];
type OptimizationStep = {
  index: number; label: string; atoms: [string, number, number, number][];
  energy: number | null; energy_shift?: number;
  rmsd_step?: number; max_step?: number; rmsd_grad?: number; max_grad?: number;
  states: Summary["states"]; transitions: Summary["transitions"];
};
let optimizationSteps: OptimizationStep[] = [];
let selectedScanGroup = "";
let resultViewerReady = false;
let pendingResultMessage: Record<string, unknown> | null = null;
let displayedResultAtomCount = 0;
let volumetricRequestId = 0;
let viewerRenderGeneration = 0;
let activeMapSource: "direct" | "excited" | "orbital" | "surface" | null = null;
let activeDirectCube: { jobId: string; url: string } | null = null;
let activeSurfaceSelection: { operation: string; primary: string; secondary: string } | null = null;

async function showDirectCube(jobId: string, url: string): Promise<void> {
  if (jobId !== selectedJob) return;
  const requestId = ++volumetricRequestId;
  const name = cubeNameFromUrl(url);
  const reportFailure = (message: string) => {
    if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
    $<HTMLDivElement>("surfaceStatus").textContent = message;
    pushArtSources();
  };
  let response: Response;
  let geometryResponse: Response;
  try {
    [response, geometryResponse] = await Promise.all([
      fetch(url),
      fetch(`/api/jobs/${jobId}/cube-geometry?name=${encodeURIComponent(name)}`),
    ]);
  } catch (error) {
    reportFailure(`cube display failed: ${error instanceof Error ? error.message : String(error)}`);
    return;
  }
  if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
  if (!response.ok) {
    reportFailure(`cube display failed (${response.status})`);
    return;
  }
  const cube = await response.blob();
  if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
  const geometry = geometryResponse.ok ? await geometryResponse.json() : null;
  if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
  if (activeCubeUrl) URL.revokeObjectURL(activeCubeUrl);
  // Whatever surface replaces this one, the old components go with it: a
  // scalar-only cube never reaches the fetch below and would otherwise leave
  // three cube-sized blobs behind.
  releaseCurrentComponents();
  activeCubeUrl = URL.createObjectURL(cube);
  activeMapSource = "direct";
  activeDirectCube = { jobId, url };
  pushOrbitalStyle();
  const components = await loadCurrentComponents(jobId, url,
    () => requestId === volumetricRequestId && jobId === selectedJob);
  if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
  pushToResultViewer({
    type: "oqp-cube",
    ...(typeof geometry?.xyz === "string" ? { xyz: geometry.xyz } : {}),
    cube: activeCubeUrl, iso: 0.05,
    sides: $<HTMLSelectElement>("moSides").value,
    ...(components ? { components, ...currentVectorSettings() } : {}),
  });
}

function syncSurfaceFiles(): void {
  const options = resultCubeFiles.map((name) =>
    `<option value="${escapeMarkup(name)}">${escapeMarkup(name)}</option>`).join("");
  $<HTMLSelectElement>("surfacePrimary").innerHTML = options;
  $<HTMLSelectElement>("surfaceSecondary").innerHTML = options;
  if (resultCubeFiles.length > 1) $<HTMLSelectElement>("surfaceSecondary").selectedIndex = 1;
  $<HTMLDivElement>("surfaceCard").style.display = resultCubeFiles.length ? "" : "none";
  pushArtSources();
}

function syncSurfaceOperation(): void {
  const arithmetic = $<HTMLSelectElement>("surfaceOperation").value !== "display";
  $<HTMLDivElement>("surfaceSecondaryWrap").style.display = arithmetic ? "" : "none";
}

async function showSurface(): Promise<void> {
  if (!selectedJob) return;
  const jobId = selectedJob;
  const requestId = ++volumetricRequestId;
  const primary = $<HTMLSelectElement>("surfacePrimary").value;
  const secondary = $<HTMLSelectElement>("surfaceSecondary").value;
  const operation = $<HTMLSelectElement>("surfaceOperation").value;
  const selectionIsCurrent = () => jobId === selectedJob &&
    primary === $<HTMLSelectElement>("surfacePrimary").value &&
    secondary === $<HTMLSelectElement>("surfaceSecondary").value &&
    operation === $<HTMLSelectElement>("surfaceOperation").value;
  const status = $<HTMLDivElement>("surfaceStatus");
  const reportFailure = (message: string) => {
    if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
    status.textContent = message;
    pushArtSources();
  };
  if (!primary || (operation !== "display" && !secondary)) return;
  const sourceUrl = `/api/jobs/${jobId}/files/${encodeURIComponent(primary)}`;
  status.textContent = operation === "display" ? `Displaying ${primary}` : `Computing ${operation}…`;
  const query = new URLSearchParams({ left: primary, right: secondary, operation });
  let response: Response;
  let geometryResponse: Response;
  try {
    [response, geometryResponse] = await Promise.all([
      fetch(operation === "display" ? sourceUrl : `/api/jobs/${jobId}/cube-combine?${query}`),
      fetch(`/api/jobs/${jobId}/cube-geometry?name=${encodeURIComponent(primary)}`),
    ]);
  } catch (error) {
    reportFailure(`surface operation failed: ${error instanceof Error ? error.message : String(error)}`);
    return;
  }
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  if (!response.ok) {
    const detail = await response.json().catch(() => null);
    if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
    reportFailure(detail?.detail ?? `surface operation failed (${response.status})`);
    return;
  }
  const cube = await response.blob();
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  const geometry = geometryResponse.ok ? await geometryResponse.json() : null;
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  if (activeCubeUrl) URL.revokeObjectURL(activeCubeUrl);
  releaseCurrentComponents();
  activeCubeUrl = URL.createObjectURL(cube);
  const cubeUrl = activeCubeUrl;
  if (operation !== "display") {
    status.textContent = `${operation}: ${primary} and ${secondary}`;
  }
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  activeMapSource = "surface";
  activeSurfaceSelection = { operation, primary, secondary };
  const components = operation === "display"
    ? await loadCurrentComponents(jobId, sourceUrl,
      () => requestId === volumetricRequestId && selectionIsCurrent())
    : undefined;
  if (requestId !== volumetricRequestId || !selectionIsCurrent()) return;
  pushOrbitalStyle();
  pushToResultViewer({
    type: "oqp-cube",
    ...(typeof geometry?.xyz === "string" ? { xyz: geometry.xyz } : {}),
    cube: cubeUrl,
    iso: Math.max(0.0001, Math.abs(+$<HTMLInputElement>("surfaceIso").value) || 0.05),
    sides: $<HTMLSelectElement>("surfaceSides").value,
    ...(components ? { components, ...currentVectorSettings() } : {}),
  });
}

$<HTMLSelectElement>("surfaceOperation").addEventListener("change", () => {
  volumetricRequestId += 1;
  syncSurfaceOperation();
});
for (const id of ["surfacePrimary", "surfaceSecondary"]) {
  $<HTMLSelectElement>(id).addEventListener("change", () => { volumetricRequestId += 1; });
}
$<HTMLButtonElement>("surfaceShow").addEventListener("click", () => { void showSurface(); });
$<HTMLSelectElement>("surfaceSides").addEventListener("change", () => { void showSurface(); });

// Everything in Analysis renders through the same Mol* page the Builder uses,
// so results look like the rest of the app rather than a second program.
window.addEventListener("message", (event) => {
  if (event.origin === window.location.origin && event.data?.type === "oqp-viewer-ready" &&
      event.source === resultFrame.contentWindow) {
    resultViewerReady = true;
    if (pendingResultMessage) {
      resultFrame.contentWindow?.postMessage(pendingResultMessage, window.location.origin);
      pendingResultMessage = null;
    }
  }
});

function pushToResultViewer(message: Record<string, unknown>): void {
  if (["oqp-structure", "oqp-normal-mode", "oqp-normal-mode-reset", "oqp-atomic-property"]
      .includes(String(message.type))) {
    volumetricRequestId += 1;
    activeMapSource = null;
    activeDirectCube = null;
  }
  if (typeof message.xyz === "string") {
    const displayedAtoms = parseAtoms(message.xyz);
    if (displayedAtoms.length) displayedResultAtomCount = displayedAtoms.length;
  } else if (message.type === "oqp-cube") {
    displayedResultAtomCount = 0;
  }
  const type = String(message.type);
  const coordinateText = typeof message.xyz === "string" ? message.xyz
    : type === "oqp-file" && typeof message.text === "string" ? message.text : "";
  const sceneAtoms = coordinateText ? parseAtoms(coordinateText) : [];
  if (type === "oqp-cube" && typeof message.cube === "string") {
    const orbital = currentOrbitalStyle();
    if (sceneAtoms.length) artAtoms = sceneAtoms;
    artScene = {
      atoms: sceneAtoms,
      cube: message.cube,
      iso: typeof message.iso === "number" ? Math.abs(message.iso) : 0.05,
      sides: typeof message.sides === "string" ? message.sides : String(orbital.sides ?? "both"),
      orbital,
      label: activeMapSource ?? "volumetric surface",
    };
    pushToArt();
  } else if (["oqp-structure", "oqp-normal-mode", "oqp-normal-mode-reset", "oqp-file"]
      .includes(type) && sceneAtoms.length) {
    setArtStructure(sceneAtoms);
  } else if (type === "oqp-atomic-property" || type === "oqp-nmr-shielding") {
    setArtStructure(artAtoms);
  }
  const renderTypes = [
    "oqp-structure", "oqp-normal-mode", "oqp-normal-mode-reset", "oqp-file", "oqp-cube",
  ];
  const styledMessage = renderTypes.includes(String(message.type))
    ? { ...message, style: { ...displayStyle },
        ...(message.type === "oqp-cube" ? { orbital: currentOrbitalStyle() } : {}) }
    : message;
  const outgoing = renderTypes.includes(String(message.type))
    ? { ...styledMessage, renderGeneration: ++viewerRenderGeneration }
    : styledMessage;
  pendingResultMessage = outgoing;
  if (!resultFrame.src) {
    resultViewerReady = false;
    resultFrame.src = "/builder3d.html";
    return;
  }
  if (resultViewerReady) {
    resultFrame.contentWindow?.postMessage(outgoing, window.location.origin);
    pendingResultMessage = null;
  }
}

function showResultFrame(index: number): void {
  const frame = resultFrames[index - 1];
  if (!frame) return;
  const xyz = `${frame.atoms.length}\n${frame.label}\n${atomsToText(frame.atoms)}\n`;
  $<HTMLSpanElement>("resultFrameLabel").textContent =
    resultFrames.length > 1 ? `${index} / ${resultFrames.length}  ${frame.label}` : frame.label;
  pushToResultViewer({ type: "oqp-structure", xyz });
}

$<HTMLInputElement>("resultFrameRange").addEventListener("input", (event) => {
  const position = +(event.target as HTMLInputElement).value;
  if (optimizationSteps.length) selectOptimizationStep(position);
  else showResultFrame(position);
});

$<HTMLSelectElement>("optimizationStepSelect").addEventListener("change", (event) => {
  const index = +(event.target as HTMLSelectElement).value;
  const position = optimizationSteps.findIndex((step) => step.index === index) + 1;
  if (position > 0) selectOptimizationStep(position);
});

$<HTMLButtonElement>("editResultStructure").addEventListener("click", () => {
  const position = +$<HTMLInputElement>("resultFrameRange").value;
  const frame = resultFrames[position - 1];
  if (!frame) return;
  clearPdbSource();
  loadedFrames = [{ label: frame.label || "result structure", atoms: frame.atoms }];
  $<HTMLInputElement>("frameRange").value = "1";
  $<HTMLInputElement>("frameRange").max = "1";
  $<HTMLDivElement>("frameRow").style.display = "none";
  showTab("builder");
  showFrame(1);
  builderStatus.textContent = `${frame.label || "result structure"} loaded for editing`;
});

async function loadOptimizationHistory(jobId: string): Promise<boolean> {
  const response = await fetch(`/api/jobs/${jobId}/optimization`);
  if (!response.ok || jobId !== selectedJob) return false;
  const data: { steps: OptimizationStep[] } = await response.json();
  if (!data.steps.length || jobId !== selectedJob) return false;
  optimizationSteps = data.steps;
  resultFrames = data.steps.map((step) => ({ label: step.label, atoms: step.atoms }));
  const select = $<HTMLSelectElement>("optimizationStepSelect");
  select.replaceChildren(...data.steps.map((step) => {
    const option = document.createElement("option");
    option.value = String(step.index);
    option.textContent = `Step ${step.index}`;
    return option;
  }));
  const slider = $<HTMLInputElement>("resultFrameRange");
  slider.max = String(data.steps.length);
  $<HTMLDivElement>("resultFrameCard").style.display = "";
  $<HTMLHeadingElement>("resultFrameTitle").textContent = "Optimization path";
  $<HTMLElement>("optimizationStepControls").style.display = "";
  if (!selectedScanGroup) renderOptimizationPath();
  selectOptimizationStep(data.steps.length);
  return true;
}

function selectOptimizationStep(position: number): void {
  const step = optimizationSteps[position - 1];
  if (!step) return;
  $<HTMLInputElement>("resultFrameRange").value = String(position);
  $<HTMLSelectElement>("optimizationStepSelect").value = String(step.index);
  if (activeCubeUrl) {
    URL.revokeObjectURL(activeCubeUrl);
    activeCubeUrl = null;
  }
  releaseCurrentComponents();
  orbitalSel.value = "";
  showResultFrame(position);
  const value = (label: string, number: number | null | undefined, digits = 6) =>
    number == null ? "" : `${label} ${number.toFixed(digits)}`;
  $<HTMLDivElement>("optimizationStepValues").textContent = [
    value("Surface energy (Ha)", step.energy, 8),
    value("RMS gradient", step.rmsd_grad), value("Max gradient", step.max_grad),
    value("RMS step", step.rmsd_step), value("Max step", step.max_step),
  ].filter(Boolean).join("  |  ");
  if (!selectedScanGroup) renderOptimizationPath(position);
  void loadSpectrum();
}

type EnergyPathPoint = { x: number; energy: number; label: string; id: string };

function escapeMarkup(value: string): string {
  return value.replace(/[&<>"']/g, (character) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  })[character] ?? character);
}

function renderEnergyPath(
  points: EnergyPathPoint[], title: string, xLabel: string, selectedId: string,
  onSelect: (id: string) => void,
): void {
  const card = $<HTMLDivElement>("pathCard");
  const plot = $<HTMLDivElement>("pathPlot");
  if (!points.length) {
    card.style.display = "none";
    plot.innerHTML = "";
    return;
  }
  const width = 640;
  const height = 230;
  const left = 58;
  const right = 18;
  const top = 18;
  const bottom = 45;
  const energies = points.map((point) => point.energy);
  const minimum = Math.min(...energies);
  const relative = energies.map((energy) => (energy - minimum) * 627.509474);
  const xMin = Math.min(...points.map((point) => point.x));
  const xMax = Math.max(...points.map((point) => point.x));
  const yMax = Math.max(...relative, 0.1);
  const xAt = (x: number) => left + (x - xMin) / (xMax - xMin || 1) * (width - left - right);
  const yAt = (y: number) => top + (1 - y / yMax) * (height - top - bottom);
  const coords = points.map((point, index) => `${xAt(point.x).toFixed(1)},${yAt(relative[index]).toFixed(1)}`);
  const circles = points.map((point, index) =>
    `<circle class="path-point${point.id === selectedId ? " selected" : ""}" data-path-id="${escapeMarkup(point.id)}" cx="${xAt(point.x).toFixed(1)}" cy="${yAt(relative[index]).toFixed(1)}" r="5"><title>${escapeMarkup(point.label)}: ${relative[index].toFixed(3)} kcal/mol</title></circle>`,
  ).join("");
  plot.innerHTML = `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="${title}">
    <line class="path-axis" x1="${left}" y1="${height - bottom}" x2="${width - right}" y2="${height - bottom}" />
    <line class="path-axis" x1="${left}" y1="${top}" x2="${left}" y2="${height - bottom}" />
    <polyline class="path-line" points="${coords.join(" ")}" />${circles}
    <text class="path-label" x="${left}" y="${height - 18}" text-anchor="middle">${xMin.toFixed(2)}</text>
    <text class="path-label" x="${width - right}" y="${height - 18}" text-anchor="middle">${xMax.toFixed(2)}</text>
    <text class="path-label" x="${(left + width - right) / 2}" y="${height - 3}" text-anchor="middle">${xLabel}</text>
    <text class="path-label" x="${left - 8}" y="${height - bottom + 4}" text-anchor="end">0</text>
    <text class="path-label" x="${left - 8}" y="${top + 4}" text-anchor="end">${yMax.toFixed(2)}</text>
    <text class="path-label" transform="translate(13 ${(top + height - bottom) / 2}) rotate(-90)" text-anchor="middle">Relative energy (kcal/mol)</text>
  </svg>`;
  plot.querySelectorAll<SVGCircleElement>("[data-path-id]").forEach((circle) => {
    circle.addEventListener("click", () => onSelect(circle.dataset.pathId ?? ""));
  });
  $<HTMLHeadingElement>("pathTitle").textContent = title;
  $<HTMLDivElement>("pathNote").textContent = "Energy is referenced to the lowest displayed structure. Select a point to inspect it.";
  card.style.display = "";
}

function renderOptimizationPath(selectedPosition = optimizationSteps.length): void {
  const points = optimizationSteps.flatMap((step, position) => step.energy == null ? [] : [{
    x: step.index, energy: step.energy, label: `Step ${step.index}`, id: String(position + 1),
  }]);
  renderEnergyPath(points, "Reaction path", "Structure step", String(selectedPosition), (id) => {
    selectOptimizationStep(+id);
  });
}

async function loadScanGroup(groupId: string): Promise<void> {
  selectedScanGroup = groupId;
  const response = await fetch(`/api/scans/${groupId}`);
  if (!response.ok || selectedScanGroup !== groupId) return;
  const data: { points: ScanPoint[] } = await response.json();
  if (selectedScanGroup !== groupId) return;
  const points = data.points.flatMap((point) => point.energy == null ? [] : [{
    x: point.value, energy: point.energy, label: point.name, id: point.job_id,
  }]);
  renderEnergyPath(points, "Bond-distance scan", "Distance (Å)", selectedJob, (id) => {
    if (id !== selectedJob) void selectJob(id);
  });
}

function hideResultPanels(): void {
  showCurrentVectorControls(false);
  orbitalCard.style.display = "none";
  modeCard.style.display = "none";
  currentMolden = null;
  orbitalSources.clear();
}

async function viewResultFile(jobId: string, name: string, url: string): Promise<void> {
  if (jobId !== selectedJob) return;
  const requestId = ++volumetricRequestId;
  const lower = name.toLowerCase();
  hideResultPanels();

  if (lower.endsWith(".cube") || lower.endsWith(".cub")) {
    showDirectCube(jobId, url);
    return;
  }

  if (lower.endsWith(".molden")) {
    const base = `/api/jobs/${jobId}/molden/${encodeURIComponent(name)}`;
    currentMolden = { jobId, name };
    const moldenFiles = [...new Set([...resultMoldenFiles, name])];
    const [orbitalFiles, modes] = await Promise.all([
      Promise.all(moldenFiles.map(async (moldenName) => {
        const sourceBase = `/api/jobs/${jobId}/molden/${encodeURIComponent(moldenName)}`;
        const data = await fetch(`${sourceBase}/orbitals`)
          .then((r) => (r.ok ? r.json() : null)).catch(() => null);
        return { name: moldenName, data };
      })),
      fetch(`${base}/modes`).then((r) => (r.ok ? r.json() : null)).catch(() => null),
    ]);
    if (requestId !== volumetricRequestId || jobId !== selectedJob) return;
    const scfFile = orbitalFiles.find((file) =>
      file.data?.orbitals?.some((orbital: { kind?: string }) => orbital.kind === "scf"));
    if (scfFile) currentMolden = { jobId, name: scfFile.name };
    const allOrbitals = orbitalFiles.flatMap((file) =>
      (file.data?.orbitals ?? []).map((orbital: Record<string, unknown>) => ({
        ...orbital,
        sourceName: file.name,
      })),
    );
    if (allOrbitals.length) {
      orbitalSources.clear();
      orbitalSel.innerHTML = '<option value="">Choose an orbital…</option>';
      const scfGroup = document.createElement("optgroup");
      scfGroup.label = "SCF orbitals";
      scfGroup.dataset.kind = "scf";
      const dysonGroup = document.createElement("optgroup");
      dysonGroup.label = "Dyson orbitals";
      dysonGroup.dataset.kind = "dyson";
      for (const o of allOrbitals) {
        const opt = document.createElement("option");
        const sourceName = String(o.sourceName);
        const index = Number(o.index);
        const sourceId = `${sourceName}\u0000${index}`;
        orbitalSources.set(sourceId, { name: sourceName, index });
        opt.value = sourceId;
        if (o.kind === "dyson") {
          opt.dataset.kind = "dyson";
          const strength = o.strength == null ? "" : ` strength=${Number(o.strength).toFixed(6)}`;
          const occupation = o.occupation == null ? "" : ` occupation=${Number(o.occupation).toFixed(6)}`;
          const sourceState = typeof o.source_state === "string" ? ` from ${o.source_state}` : "";
          opt.textContent = `Dyson ${o.dyson_kind}${sourceState}, root ${o.state_index}${strength}${occupation}`;
          dysonGroup.appendChild(opt);
        } else {
          opt.dataset.kind = "scf";
          const occ = o.occupancy == null ? "" : ` occ=${Number(o.occupancy).toFixed(4)}`;
          opt.textContent = `MO ${index}  E=${Number(o.energy).toFixed(4)} Ha${occ} ${o.spin}`;
          scfGroup.appendChild(opt);
        }
      }
      if (scfGroup.children.length) orbitalSel.appendChild(scfGroup);
      if (dysonGroup.children.length) orbitalSel.appendChild(dysonGroup);
      // IP/EA output carries a dedicated Dyson Molden file. Use that durable
      // output contract as well as the parsed labels, so a partial orbital
      // metadata response cannot hide the Dyson display mode.
      const hasDyson = dysonGroup.children.length > 0 ||
        moldenFiles.some((moldenName) => /dyson/i.test(moldenName));
      const dysonMapOption = $<HTMLOptionElement>("dysonMapKind");
      dysonMapOption.hidden = !hasDyson;
      dysonMapOption.disabled = !hasDyson;
      if (!hasDyson && mapKind.value === "dyson") mapKind.value = "mo";
      orbitalCard.style.display = "";
      selectOrbitalClass(mapKind.value);
      pushArtSources();
    }
    if (modes?.modes?.length) {
      modeSel.innerHTML = '<option value="">Choose a normal mode…</option>';
      for (const m of modes.modes) {
        const opt = document.createElement("option");
        opt.value = String(m.index);
        const ir = m.intensity != null ? `  IR ${m.intensity.toFixed(2)}` : "";
        opt.textContent = `Mode ${m.index}  ${m.frequency.toFixed(1)} cm⁻¹${ir}`;
        modeSel.appendChild(opt);
      }
      modeCard.style.display = "";
      modePlay.disabled = true;
      modePause.disabled = true;
      modeArrows.disabled = true;
      modeReset.disabled = true;
    }
    if (!allOrbitals.length && !modes?.modes?.length) {
      await openAsStructure(name, url, requestId);
    } else {
      // Opening a Molden result must show its molecule before the user picks
      // an orbital. Reset map already did this, which is why it exposed the
      // missing initial structure message.
      showMoldenStructure();
    }
    return;
  }

  await openAsStructure(name, url, requestId);
}

// Reads the file through the same importer the Builder uses, so OpenQP logs,
// JSON, inputs and trajectories all render as geometry here.
async function openAsStructure(
  name: string, url: string, requestId = ++volumetricRequestId,
): Promise<void> {
  try {
    const response = await fetch(url);
    if (requestId !== volumetricRequestId) return;
    const blob = await response.blob();
    if (requestId !== volumetricRequestId) return;
    const body = new FormData();
    body.append("file", new File([blob], name));
    const res = await fetch("/api/structure/open", { method: "POST", body });
    if (requestId !== volumetricRequestId) return;
    if (!res.ok) return;
    const data = await res.json();
    if (requestId !== volumetricRequestId) return;
    resultFrames = data.frames;
    const slider = $<HTMLInputElement>("resultFrameRange");
    slider.max = String(resultFrames.length);
    slider.value = String(resultFrames.length);
    $<HTMLDivElement>("resultFrameCard").style.display =
      resultFrames.length > 1 ? "" : "none";
    $<HTMLElement>("optimizationStepControls").style.display = "none";
    $<HTMLHeadingElement>("resultFrameTitle").textContent = "Trajectory";
    showResultFrame(resultFrames.length);
  } catch {
    // Nothing renderable in this file; the download link still works.
  }
}

// Isovalues that suit each field: orbitals peak near 0.05, the density is
// conventionally drawn at 0.002, and the MEP is a potential in Hartree/e.
const MAP_ISO: Record<string, number> = {
  mo: 0.05, dyson: 0.05, density: 0.002, spin: 0.004, esp: 0.03,
};

const mapKind = $<HTMLSelectElement>("mapKind");

function selectOrbitalClass(map: string): void {
  const wanted = map === "dyson" ? "dyson" : "scf";
  orbitalSel.querySelectorAll<HTMLOptGroupElement>("optgroup[data-kind]").forEach((group) => {
    group.hidden = group.dataset.kind !== wanted;
  });
  const selected = orbitalSel.selectedOptions[0];
  if (selected?.dataset.kind === wanted) return;
  const first = [...orbitalSel.options].find((option) => option.dataset.kind === wanted);
  orbitalSel.value = first?.value ?? "";
}

function selectedOrbitalSource(): OrbitalSource | null {
  return orbitalSources.get(orbitalSel.value) ?? null;
}

function showMoldenStructure(): void {
  if (!currentMolden) return;
  const requestId = ++volumetricRequestId;
  if (activeCubeUrl) {
    URL.revokeObjectURL(activeCubeUrl);
    activeCubeUrl = null;
  }
  releaseCurrentComponents();
  const molden = currentMolden;
  const base = `/api/jobs/${molden.jobId}/molden/${encodeURIComponent(molden.name)}`;
  fetch(`${base}/geom.xyz`)
    .then((response) => response.text())
    .then((xyz) => {
      if (requestId !== volumetricRequestId) return;
      if (currentMolden?.jobId !== molden.jobId || currentMolden.name !== molden.name) return;
      pushToResultViewer({ type: "oqp-structure", xyz });
    });
}

function mapUrl(base: string, source: OrbitalSource | null): string {
  return (mapKind.value === "mo" || mapKind.value === "dyson")
    ? `/api/jobs/${currentMolden!.jobId}/molden/${encodeURIComponent(source!.name)}/cube?mo=${source!.index}`
    : `${base}/map?kind=${mapKind.value}`;
}

function showOrbital(): void {
  if (!currentMolden) return;
  const molden = currentMolden;
  const requestId = ++volumetricRequestId;
  pushOrbitalStyle();
  const isOrbital = mapKind.value === "mo" || mapKind.value === "dyson";
  $<HTMLDivElement>("orbitalPick").style.display = isOrbital ? "" : "none";
  const source = selectedOrbitalSource();
  if (isOrbital && !source) return;
  const base = `/api/jobs/${molden.jobId}/molden/${encodeURIComponent(molden.name)}`;
  const iso = +isoRange.value;
  // Mol* runs inside an iframe. In the standalone app only this top-level
  // window's fetch is routed through the sidecar, so pass its cube a local
  // blob URL instead of leaving the iframe to request /api itself.
  Promise.all([
    fetch(`${base}/geom.xyz`).then((response) => response.text()),
    fetch(mapUrl(base, source)).then((response) => response.text()),
  ])
    .then(([xyz, cube]) => {
      if (requestId !== volumetricRequestId) return;
      if (currentMolden?.jobId !== molden.jobId || currentMolden.name !== molden.name) return;
      if (activeCubeUrl) URL.revokeObjectURL(activeCubeUrl);
      releaseCurrentComponents();
      activeCubeUrl = URL.createObjectURL(new Blob([cube], { type: "text/plain" }));
      activeMapSource = "orbital";
      showCurrentVectorControls(false);
      pushToResultViewer({
        type: "oqp-cube",
        xyz,
        cube: activeCubeUrl,
        iso,
        sides: $<HTMLSelectElement>("moSides").value,
      });
    });
}
orbitalSel.addEventListener("change", showOrbital);
orbitalReset.addEventListener("click", () => {
  if (!currentMolden) return;
  orbitalSel.value = "";
  showMoldenStructure();
});
isoRange.addEventListener("change", showOrbital);
function redrawActiveVolume(): void {
  if (activeMapSource === "direct" && activeDirectCube) {
    showDirectCube(activeDirectCube.jobId, activeDirectCube.url);
  } else if (activeMapSource === "excited") void showExcitedMap();
  else if (activeMapSource === "surface") void showSurface();
  else if (activeMapSource === "orbital") showOrbital();
}

function showArtSource(value: string): void {
  if (value === "molecule") {
    const molecule = artScene.atoms;
    if (!molecule.length) return;
    volumetricRequestId += 1;
    activeMapSource = null;
    activeDirectCube = null;
    setArtStructure(molecule);
    return;
  }
  const separator = value.indexOf(":");
  const category = separator < 0 ? "" : value.slice(0, separator);
  const source = separator < 0 ? "" : value.slice(separator + 1);
  if (category === "orbital-entry" && currentMolden) {
    const orbitalValue = decodeURIComponent(source);
    const option = [...orbitalSel.options].find((candidate) => candidate.value === orbitalValue);
    if (!option?.dataset.kind || !orbitalSources.has(orbitalValue)) return;
    orbitalSel.value = orbitalValue;
    mapKind.value = option.dataset.kind === "dyson" ? "dyson" : "mo";
    selectOrbitalClass(mapKind.value);
    isoRange.value = String(MAP_ISO[mapKind.value]);
    showOrbital();
  } else if (category === "orbital" && currentMolden && [...mapKind.options]
      .some((option) => option.value === source && !option.hidden)) {
    mapKind.value = source;
    if (source === "mo" || source === "dyson") selectOrbitalClass(source);
    isoRange.value = String(MAP_ISO[source] ?? 0.05);
    showOrbital();
  } else if (category === "excited" && excitedAnalysis?.available && [...excitedMap.options]
      .some((option) => option.value === source)) {
    excitedMap.value = source;
    $<HTMLDivElement>("excitedPairWrap").style.display = source.startsWith("nto_") ? "" : "none";
    void showExcitedMap();
  } else if (value === "surface-combination") {
    if (!activeSurfaceSelection || activeSurfaceSelection.operation === "display") return;
    $<HTMLSelectElement>("surfaceOperation").value = activeSurfaceSelection.operation;
    $<HTMLSelectElement>("surfacePrimary").value = activeSurfaceSelection.primary;
    $<HTMLSelectElement>("surfaceSecondary").value = activeSurfaceSelection.secondary;
    syncSurfaceOperation();
    void showSurface();
  } else if (category === "surface" && resultCubeFiles.includes(source)) {
    void showDirectCube(selectedJob,
      `/api/jobs/${selectedJob}/files/${encodeURIComponent(source)}`);
  }
}

$<HTMLInputElement>("nmrAcid").addEventListener("change", () => {
  const on = $<HTMLInputElement>("nmrAcid").checked;
  $<HTMLDivElement>("acidGridRow").style.display = on ? "" : "none";
  // Offering CGO beside a ticked ACID box would propose a run OpenQP refuses.
  const gauge = $<HTMLSelectElement>("nmrGauge");
  if (on) gauge.value = "giao";
  gauge.disabled = on;
});

$<HTMLSelectElement>("moSides").addEventListener("change", redrawActiveVolume);
["currentVectorToggle", "currentVectorDensity", "currentVectorScale"].forEach((id) => {
  // change, not input: each redraw re-fetches and re-samples the components,
  // which is far too much work to do on every step of a drag.
  $<HTMLInputElement>(id).addEventListener("change", redrawActiveVolume);
});
mapKind.addEventListener("change", () => {
  // Each field has its own natural contour, so move the slider with it.
  isoRange.value = String(MAP_ISO[mapKind.value] ?? 0.05);
  if (mapKind.value === "mo" || mapKind.value === "dyson") selectOrbitalClass(mapKind.value);
  showOrbital();
});

function modeEquilibriumXyz(mode: { atoms: { element: string; position: number[] }[] }): string {
  return `${mode.atoms.length}\nnormal-mode equilibrium geometry\n` + mode.atoms
    .map((atom) => `${atom.element} ${atom.position.map((value) => value.toFixed(6)).join(" ")}`)
    .join("\n") + "\n";
}

async function showMode(): Promise<void> {
  modePlay.disabled = true;
  modePause.disabled = true;
  modeArrows.disabled = true;
  if (!currentMolden || !modeSel.value) return;
  const molden = currentMolden;
  const selectedMode = modeSel.value;
  const requestId = ++volumetricRequestId;
  const base = `/api/jobs/${molden.jobId}/molden/${encodeURIComponent(molden.name)}`;
  const vectors = await fetch(`${base}/mode?mode=${selectedMode}`);
  if (requestId !== volumetricRequestId || currentMolden !== molden || modeSel.value !== selectedMode) return;
  if (!vectors.ok) return;
  const mode = await vectors.json();
  if (requestId !== volumetricRequestId || currentMolden !== molden || modeSel.value !== selectedMode) return;
  pushToResultViewer({
    type: "oqp-normal-mode",
    xyz: modeEquilibriumXyz(mode),
    mode,
    amplitude: Number(amplitudeRange.value),
    arrows: modeArrows.checked,
    playing: false,
  });
  modePlay.disabled = false;
  modePause.disabled = false;
  modeArrows.disabled = false;
  modeReset.disabled = false;
}

modeSel.addEventListener("change", showMode);
amplitudeRange.addEventListener("change", showMode);
modeArrows.addEventListener("change", () =>
  pushToResultViewer({ type: "oqp-normal-mode-arrows", arrows: modeArrows.checked }));
modePlay.addEventListener("click", () =>
  pushToResultViewer({ type: "oqp-normal-mode-play" }));
modePause.addEventListener("click", () =>
  pushToResultViewer({ type: "oqp-normal-mode-pause" }));
modeReset.addEventListener("click", async () => {
  if (!currentMolden) return;
  const molden = currentMolden;
  const requestId = ++volumetricRequestId;
  const base = `/api/jobs/${molden.jobId}/molden/${encodeURIComponent(molden.name)}`;
  const response = await fetch(`${base}/geom.xyz`);
  if (requestId !== volumetricRequestId || currentMolden !== molden) return;
  if (!response.ok) return;
  const xyz = await response.text();
  if (requestId !== volumetricRequestId || currentMolden !== molden) return;
  pushToResultViewer({ type: "oqp-normal-mode-reset", xyz });
  modeSel.value = "";
  modePlay.disabled = true;
  modePause.disabled = true;
  modeArrows.disabled = true;
  modeReset.disabled = true;
});

// ---------- results summary ----------
type Summary = {
  energy: { total?: number; components?: Record<string, number>;
            final_states?: Record<string, number> };
  scf: { method?: string; energy?: number; iterations?: number; converged?: boolean };
  states: { index: number; total: number; excitation_ev: number;
            excitation_nm: number | null; oscillator: number | null }[];
  transitions: { from: number; to: number; excitation_ev: number;
                 oscillator: number | null }[];
  frequencies: { index: number; frequency: number; ir: number | null; raman: number | null }[];
  thermochemistry: Record<string, number>;
  charges: Record<string, number[]>;
  nmr: { atom: number; dia: number; para_uncoupled: number; para_coupled: number;
         total_uncoupled: number; total_coupled: number }[];
  dipole: { x: number; y: number; z: number; total_au: number; total_debye: number } | null;
  symmetry: { point_group: string; detected: string | null; enabled: boolean } | null;
  units: { ir: string; raman: string };
  has_frequencies: boolean;
  has_states: boolean;
  has_oscillators: boolean;
  ekt: { ip: { index: number; binding_ev: number; strength: number }[];
         ea: { index: number; binding_ev: number; strength: number }[] };
  has_ekt_ip: boolean;
  has_ekt_ea: boolean;
  excited_state_optimized: number | null;
  has_nmr: boolean;
};

let summary: Summary | null = null;
type AtomicProperty = { label: string; values: number[]; unit: string };
const atomicProperties = new Map<string, AtomicProperty>();

const ENERGY_LABELS: Record<string, string> = {
  total: "Total energy",
  one_electron: "One-electron energy",
  two_electron: "Two-electron energy",
  nuclear_repulsion: "Nuclear repulsion",
  potential: "Potential energy",
  kinetic: "Kinetic energy",
  virial_ratio: "Virial ratio (V/T)",
};

const THERMO_LABELS: Record<string, string> = {
  temperature: "Temperature (K)",
  pressure: "Pressure (atm)",
  zpe: "Zero-point energy",
  internal_energy: "Internal energy U",
  enthalpy: "Enthalpy H",
  gibbs_free_energy: "Gibbs free energy G",
};

function rows(pairs: [string, string][]): string {
  return pairs.length
    ? `<table class="sum">${pairs
        .map(([k, v]) => `<tr><td class="k">${escapeHtml(k)}</td><td class="v">${escapeHtml(v)}</td></tr>`)
        .join("")}</table>`
    : "";
}

function fixed(value: number, digits = 6): string {
  return Number.isFinite(value) ? value.toFixed(digits) : "—";
}

function renderSummary(data: Summary): void {
  atomicProperties.clear();
  const parts: string[] = [];
  const finalStates = Object.entries(data.energy.final_states ?? {})
    .map(([index, value]) => [Number(index), value] as const)
    .sort(([a], [b]) => a - b);

  const energy: [string, string][] = [];
  if (data.scf.energy !== undefined) {
    energy.push([`${data.scf.method ?? "SCF"} energy (Ha)`, fixed(data.scf.energy, 8)]);
  }
  if (data.energy.total !== undefined) {
    energy.push(["Total energy (Ha)", fixed(data.energy.total, 8)]);
  }
  if (data.scf.iterations !== undefined) {
    energy.push(["SCF iterations", String(data.scf.iterations)]);
  }
  if (data.scf.converged !== undefined) {
    energy.push(["Converged", data.scf.converged ? "yes" : "no"]);
  }
  for (const [key, label] of Object.entries(ENERGY_LABELS)) {
    const value = data.energy.components?.[key];
    if (value === undefined) continue;
    // The virial ratio is dimensionless; everything else here is in Hartree.
    energy.push([key === "virial_ratio" ? label : `${label} (Ha)`, fixed(value)]);
  }
  if (finalStates.length) {
    const stateRows = finalStates.map(([index, value]) =>
      `<tr><td class="k">State ${index}</td><td class="v">${fixed(value, 8)} Ha</td></tr>`,
    ).join("");
    parts.push(`<div class="sum-title">Final state energies</div><table class="sum">${stateRows}</table>`);
  }
  if (energy.length) {
    const title = finalStates.length ? "SCF reference" : "Energy";
    parts.push(`<div class="sum-title">${title}</div>${rows(energy)}`);
  }

  if (data.symmetry?.point_group) {
    const upper = (text: string) =>
      text.charAt(0).toUpperCase() + text.slice(1);
    parts.push('<div class="sum-title">Symmetry</div>' + rows([
      ["Point group", upper(data.symmetry.point_group)],
      ["Detected", upper(data.symmetry.detected ?? data.symmetry.point_group)],
      ["Used in the run", data.symmetry.enabled ? "yes" : "no"],
    ]));
  }

  if (data.dipole) {
    parts.push('<div class="sum-title">Dipole moment</div>' + rows([
      ["x, y, z (a.u.)",
       `${fixed(data.dipole.x, 4)}, ${fixed(data.dipole.y, 4)}, ${fixed(data.dipole.z, 4)}`],
      ["Magnitude (Debye)", fixed(data.dipole.total_debye, 4)],
    ]));
  }

  if (data.states.length > 1) {
    const s0 = data.states.find((state) => state.index === 0);
    const head = "<tr><th>State</th><th>E (Ha)</th><th>ΔE from S0 (eV)</th><th>λ (nm)</th><th>f</th></tr>";
    const body = data.states.map((state) =>
      `<tr><td class="k">S${state.index}</td>` +
      `<td class="v">${state.total.toFixed(8)}</td>` +
      `<td class="v">${state.excitation_ev.toFixed(3)}</td>` +
      `<td class="v">${state.excitation_nm ? state.excitation_nm.toFixed(1) : "—"}</td>` +
      `<td class="v">${state.oscillator != null ? state.oscillator.toFixed(4) : "—"}</td></tr>`,
    ).join("");
    parts.push(`<div class="sum-title">Excited states</div>` +
      `<table class="sum">${head}${body}</table>` +
      (s0 ? `<div class="hint">Reference energy: S0 = ${s0.total.toFixed(8)} Ha</div>` : "") +
      (data.has_oscillators ? "" :
        '<div class="hint">the output carries no oscillator strengths</div>'));
  }

  if (data.transitions.length) {
    const head = "<tr><th>Transition</th><th>ΔE (eV)</th><th>λ (nm)</th><th>f</th></tr>";
    const body = data.transitions.map((transition) =>
      `<tr><td class="k">S${transition.from} → S${transition.to}</td>` +
      `<td class="v">${transition.excitation_ev.toFixed(4)}</td>` +
      `<td class="v">${(1239.841984 / transition.excitation_ev).toFixed(1)}</td>` +
      `<td class="v">${transition.oscillator != null ? transition.oscillator.toFixed(4) : "—"}</td></tr>`,
    ).join("");
    parts.push(`<div class="sum-title">State-to-state transitions</div><table class="sum">${head}${body}</table>`);
  }

  for (const [kind, roots] of [["IP", data.ekt.ip], ["EA", data.ekt.ea]] as const) {
    if (!roots.length) continue;
    const head = "<tr><th>Dyson root</th><th>Binding energy (eV)</th><th>Strength</th></tr>";
    const body = roots.map((root) =>
      `<tr><td class="k">${root.index}</td><td class="v">${root.binding_ev.toFixed(4)}</td>` +
      `<td class="v">${root.strength.toFixed(6)}</td></tr>`,
    ).join("");
    parts.push(`<div class="sum-title">EKT ${kind} Dyson roots</div><table class="sum">${head}${body}</table>`);
  }

  if (data.frequencies.length) {
    const head = `<tr><th>Mode</th><th>cm⁻¹</th><th>IR (${escapeHtml(data.units.ir)})</th>` +
      `<th>Raman (${escapeHtml(data.units.raman)})</th></tr>`;
    const body = data.frequencies.map((mode) =>
      `<tr><td class="k">${mode.index}</td><td class="v">${mode.frequency.toFixed(1)}</td>` +
      `<td class="v">${mode.ir != null ? mode.ir.toFixed(3) : "—"}</td>` +
      `<td class="v">${mode.raman != null ? mode.raman.toFixed(2) : "—"}</td></tr>`,
    ).join("");
    const imaginary = data.frequencies.filter((m) => m.frequency < 0).length;
    parts.push(`<div class="sum-title">Vibrations</div>` +
      `<table class="sum">${head}${body}</table>` +
      (imaginary
        ? `<div class="hint">${imaginary} imaginary mode(s): this is a saddle point</div>`
        : ""));
  }

  if (data.nmr.length) {
    atomicProperties.set("nmr", {
      label: "NMR total coupled shielding", unit: "ppm",
      values: data.nmr.map((row) => row.total_coupled),
    });
    const head = "<tr><th>Atom</th><th>σdia</th><th>σpara (u)</th><th>σpara (c)</th><th>σtotal (u)</th><th>σtotal (c)</th></tr>";
    const body = data.nmr.map((row) =>
      `<tr><td class="k">${escapeHtml(row.atom)}</td><td class="v">${fixed(row.dia, 3)}</td>` +
      `<td class="v">${fixed(row.para_uncoupled, 3)}</td><td class="v">${fixed(row.para_coupled, 3)}</td>` +
      `<td class="v">${fixed(row.total_uncoupled, 3)}</td><td class="v">${fixed(row.total_coupled, 3)}</td></tr>`,
    ).join("");
    parts.push(`<div class="sum-title">NMR shielding (ppm) <button class="ghost" id="showNmrMap" type="button" style="float:right;padding:.15rem .7rem;font-size:12px">Show 3D map</button></div>` +
      `<table class="sum">${head}${body}</table><div class="hint">σtotal (c) includes the coupled magnetic response.</div>`);
  }

  const thermo = Object.entries(THERMO_LABELS)
    .filter(([key]) => data.thermochemistry[key] !== undefined)
    .map(([key, label]) => [label, fixed(data.thermochemistry[key])] as [string, string]);
  if (thermo.length) {
    parts.push(`<div class="sum-title">Thermochemistry</div>${rows(thermo)}`);
  }

  const charge = Object.entries(data.charges)[0];
  if (charge) {
    const [source, values] = charge;
    parts.push(`<div class="sum-title">Partial charges (${escapeHtml(source)})</div>` +
      rows(values.map((q, i) => [`Atom ${i + 1}`, fixed(q, 4)] as [string, string])));
  }
  for (const [source, values] of Object.entries(data.charges)) {
    if (values.length) atomicProperties.set(`charge:${source}`, {
      label: `${source.toUpperCase()} partial charge`, values, unit: "e",
    });
  }

  const body = $<HTMLDivElement>("summaryBody");
  body.innerHTML = parts.length
    ? parts.join("")
    : '<div class="hint">nothing summarisable in this job\u2019s output yet</div>';
  $<HTMLDivElement>("summaryCard").style.display = parts.length ? "" : "none";
  document.getElementById("showNmrMap")?.addEventListener("click", () => {
    showAtomicProperty("nmr");
  });
  const propertySelect = $<HTMLSelectElement>("atomicProperty");
  propertySelect.innerHTML = [...atomicProperties.entries()].map(([key, property]) =>
    `<option value="${escapeMarkup(key)}">${escapeMarkup(property.label)}</option>`).join("");
  $<HTMLDivElement>("atomicPropertyCard").style.display = atomicProperties.size ? "" : "none";
}

function showAtomicProperty(key = $<HTMLSelectElement>("atomicProperty").value): void {
  const property = atomicProperties.get(key);
  if (!property) return;
  const status = $<HTMLDivElement>("propertyStatus");
  if (!displayedResultAtomCount) {
    status.textContent = "Show a result structure before mapping this property.";
    $<HTMLDivElement>("propertyLegend").style.display = "none";
    return;
  }
  if (property.values.length !== displayedResultAtomCount) {
    status.textContent = `This property has ${property.values.length} values, but the displayed ` +
      `structure has ${displayedResultAtomCount} atoms.`;
    $<HTMLDivElement>("propertyLegend").style.display = "none";
    return;
  }
  status.textContent = "";
  pushToResultViewer({
    type: "oqp-atomic-property", values: property.values,
    label: property.label, unit: property.unit,
  });
  const minimum = Math.min(...property.values);
  const maximum = Math.max(...property.values);
  $<HTMLSpanElement>("propertyMinimum").textContent = `${fixed(minimum, 3)} ${property.unit}`;
  $<HTMLSpanElement>("propertyMaximum").textContent = `${fixed(maximum, 3)} ${property.unit}`;
  $<HTMLDivElement>("propertyLegend").style.display = "";
}

$<HTMLButtonElement>("showAtomicProperty").addEventListener("click", () => showAtomicProperty());

// ---------- spectra ----------
const SPECTRA: { value: string; label: string; needs: "freq" | "states" | "exopt" | "ip" | "ea" }[] = [
  { value: "ir", label: "IR absorption", needs: "freq" },
  { value: "raman", label: "Raman", needs: "freq" },
  { value: "absorption", label: "Absorption (S0 geometry)", needs: "states" },
  { value: "emission", label: "Emission (excited-state geometry)", needs: "exopt" },
  { value: "esa", label: "Excited-state absorption (excited-state geometry)", needs: "states" },
  { value: "photoelectron", label: "Photoelectron spectrum (IP)", needs: "ip" },
  { value: "inverse_photoelectron", label: "Inverse photoelectron spectrum (EA)", needs: "ea" },
];

const specKind = $<HTMLSelectElement>("specKind");
const specShape = $<HTMLSelectElement>("specShape");
const specWidth = $<HTMLInputElement>("specWidth");
const specState = $<HTMLInputElement>("specState");

function buildSpectrumList(data: Summary): void {
  const usable = SPECTRA.filter((entry) =>
    entry.needs === "freq" ? data.has_frequencies
      : entry.needs === "states" ? data.has_states
        : entry.needs === "exopt" ? data.has_states && data.excited_state_optimized !== null
        : entry.needs === "ip" ? data.has_ekt_ip : data.has_ekt_ea);
  specKind.innerHTML = "";
  for (const entry of usable) {
    const option = document.createElement("option");
    option.value = entry.value;
    option.textContent = entry.label;
    specKind.appendChild(option);
  }
  if (data.has_states) specKind.value = "absorption";
  else if (data.has_ekt_ip) specKind.value = "photoelectron";
  else if (data.has_ekt_ea) specKind.value = "inverse_photoelectron";
  applyPendingAnalysis();
  $<HTMLDivElement>("spectrumCard").style.display = usable.length ? "" : "none";
  if (usable.length) syncSpectrumControls();
}

// Vibrational widths are in cm-1, electronic ones in eV, so the slider has to
// change meaning with the spectrum.
function vibrationalKind(): boolean {
  return specKind.value === "ir" || specKind.value === "raman";
}

function syncSpectrumControls(): void {
  const vibrational = vibrationalKind();
  specWidth.min = vibrational ? "1" : "1";
  specWidth.max = vibrational ? "100" : "100";
  const width = widthValue();
  $<HTMLSpanElement>("specWidthLabel").textContent =
    vibrational ? `${width.toFixed(0)} cm⁻¹` : `${width.toFixed(2)} eV`;
  $<HTMLDivElement>("specStateWrap").style.display =
    specKind.value === "emission" || specKind.value === "esa" ? "" : "none";
}

function widthValue(): number {
  const raw = +specWidth.value;
  return vibrationalKind() ? raw : raw / 100;   // slider steps of 0.01 eV
}

async function loadSpectrum(): Promise<void> {
  const jobId = selectedJob;
  if (!jobId) return;
  syncSpectrumControls();
  const query = new URLSearchParams({
    kind: specKind.value,
    shape: specShape.value,
    fwhm: String(widthValue()),
    state: specState.value,
  });
  const selectedStep = optimizationSteps[+$<HTMLInputElement>("resultFrameRange").value - 1];
  if (selectedStep) query.set("step", String(selectedStep.index));
  const data = await (await fetch(`/api/jobs/${jobId}/spectrum?${query}`)).json();
  if (jobId !== selectedJob) return;
  drawSpectrum(data);
}

type SpectrumData = {
  available: boolean; reason?: string;
  x: number[]; y: number[]; x_nm?: number[];
  sticks: { position: number; intensity: number; position_nm?: number }[];
  x_label: string; y_label: string; reverse_x: boolean;
  estimated_intensities?: boolean; title?: string;
};

function drawSpectrum(data: SpectrumData): void {
  lastSpectrum = data;
  const plot = $<HTMLDivElement>("specPlot");
  const note = $<HTMLDivElement>("specNote");
  if (!data.available || !data.x.length) {
    plot.innerHTML = "";
    note.textContent = data.reason ?? "no spectrum for this job";
    return;
  }

  const w = 520, h = 240, pad = { l: 52, r: 12, t: 10, b: 34 };
  const xs = data.x, ys = data.y;
  const xMin = Math.min(...xs), xMax = Math.max(...xs);
  const yMax = Math.max(...ys, 1e-12);
  const flip = data.reverse_x;
  const sx = (x: number) => {
    const t = (x - xMin) / (xMax - xMin || 1);
    return pad.l + (flip ? 1 - t : t) * (w - pad.l - pad.r);
  };
  const sy = (y: number) => h - pad.b - (y / yMax) * (h - pad.t - pad.b);

  const path = xs.map((x, i) => `${i ? "L" : "M"}${sx(x).toFixed(1)},${sy(ys[i]).toFixed(1)}`)
    .join("");
  const stickMax = Math.max(...data.sticks.map((s) => s.intensity), 1e-12);
  const sticks = data.sticks.map((stick) => {
    const x = sx(stick.position).toFixed(1);
    const top = sy((stick.intensity / stickMax) * yMax * 0.92).toFixed(1);
    return `<line x1="${x}" y1="${h - pad.b}" x2="${x}" y2="${top}" ` +
      `stroke="var(--text-dim)" stroke-width="1" opacity="0.75"><title>` +
      `${stick.position.toFixed(2)}  I=${stick.intensity.toPrecision(3)}</title></line>`;
  }).join("");

  // Four ticks per axis is enough to read a spectrum without clutter.
  const ticks = [0, 1, 2, 3, 4].map((i) => {
    const value = xMin + ((xMax - xMin) * i) / 4;
    return `<text x="${sx(value).toFixed(1)}" y="${h - pad.b + 14}" fill="var(--text-dim)" ` +
      `font-size="10" text-anchor="middle">${value >= 100 ? value.toFixed(0) : value.toFixed(2)}</text>`;
  }).join("");
  const yTicks = [0, 0.5, 1].map((f) =>
    `<text x="${pad.l - 6}" y="${(sy(f * yMax) + 3).toFixed(1)}" fill="var(--text-dim)" ` +
    `font-size="10" text-anchor="end">${(f * yMax).toPrecision(2)}</text>`).join("");

  plot.innerHTML =
    `<svg viewBox="0 0 ${w} ${h}" width="100%" role="img" aria-label="${escapeHtml(data.x_label)}">` +
    `<rect x="${pad.l}" y="${pad.t}" width="${w - pad.l - pad.r}" height="${h - pad.t - pad.b}" ` +
    `fill="none" stroke="var(--border)"/>` +
    sticks +
    `<path d="${path}" fill="none" stroke="var(--accent)" stroke-width="1.6"/>` +
    ticks + yTicks +
    `<text x="${(pad.l + w - pad.r) / 2}" y="${h - 4}" fill="var(--text-dim)" font-size="11" ` +
    `text-anchor="middle">${escapeHtml(data.x_label)}</text>` +
    `<text x="12" y="${(h - pad.b + pad.t) / 2}" fill="var(--text-dim)" font-size="11" ` +
    `text-anchor="middle" transform="rotate(-90 12 ${(h - pad.b + pad.t) / 2})">${escapeHtml(data.y_label)}</text>` +
    `</svg>`;

  const shape = specShape.value === "lorentzian" ? "Lorentzian"
    : specShape.value === "gaussian" ? "Gaussian" : "pseudo-Voigt";
  note.textContent = `${data.title ?? ""} ${shape} broadening, ` +
    `${data.sticks.length} transition${data.sticks.length === 1 ? "" : "s"}` +
    (data.estimated_intensities
      ? " — no oscillator strengths in the output, so every line is drawn at equal intensity"
      : "");
}

for (const control of [specKind, specShape, specState]) {
  control.addEventListener("change", () => { void loadSpectrum(); });
}
specWidth.addEventListener("input", syncSpectrumControls);
specWidth.addEventListener("change", () => { void loadSpectrum(); });

// ---------- image and data export ----------
document.querySelectorAll<HTMLButtonElement>(".snapshot").forEach((button) => {
  button.addEventListener("click", () => {
    const frame = document.getElementById(button.dataset.frame!) as HTMLIFrameElement | null;
    frame?.contentWindow?.postMessage({ type: "oqp-snapshot" }, window.location.origin);
  });
});

window.addEventListener("message", (event) => {
  if (event.origin !== window.location.origin) return;
  if (event.data?.type !== "oqp-snapshot-data") return;
  const image = event.data.image as string | null;
  if (!image) return;
  const link = document.createElement("a");
  link.href = image;
  link.download = "oqp-studio-view.png";
  link.click();
});

let lastSpectrum: SpectrumData | null = null;

$<HTMLButtonElement>("specCsv").addEventListener("click", () => {
  if (!lastSpectrum?.available) return;
  const header = `# ${lastSpectrum.x_label}, ${lastSpectrum.y_label}\n`;
  const curve = lastSpectrum.x
    .map((x, i) => `${x.toPrecision(8)},${lastSpectrum!.y[i].toPrecision(8)}`)
    .join("\n");
  const sticks = lastSpectrum.sticks
    .map((stick) => `# stick,${stick.position.toPrecision(8)},${stick.intensity.toPrecision(6)}`)
    .join("\n");
  download(`${specKind.value}-spectrum.csv`, `${header}${sticks}\n${curve}\n`);
});

// ---------- measurements ----------
// The viewer reports what the user clicked; both tabs show the read-out.
window.addEventListener("message", (event) => {
  if (event.origin !== window.location.origin) return;
  if (event.data?.type === "oqp-qmmm-selection") {
    qmAtoms.clear();
    for (const atom of event.data.atoms as { index: number; label: string }[]) {
      if (Number.isInteger(atom.index) && atom.index >= 0) qmAtoms.set(atom.index, atom.label);
    }
    updateQmAtoms();
    updateInpPreview();
    return;
  }
  if (event.data?.type !== "oqp-measure") return;
  const { labels, kind, value, unit } = event.data as
    { labels: string[]; kind: string | null; value: number | null; unit: string | null };
  for (const id of ["measurePreview", "measureResult"]) {
    const box = document.getElementById(id);
    if (!box) continue;
    box.replaceChildren();
    if (labels.length) {
      const picked = document.createElement("span");
      picked.className = "picked";
      picked.textContent = labels.join(" – ");
      box.appendChild(picked);
      if (kind && value != null) {
        const result = document.createElement("span");
        result.className = "value";
        result.textContent = `${kind}: ${value.toFixed(3)} ${unit ?? ""}`;
        box.appendChild(result);
      }
      if (labels.length < 2) {
        const prompt = document.createElement("span");
        prompt.className = "picked";
        prompt.textContent = "pick another atom";
        box.appendChild(prompt);
      }
    }
    box.classList.toggle("on", Boolean(labels.length));
  }
});

$<HTMLButtonElement>("clearQmAtoms").addEventListener("click", () => {
  qmAtoms.clear();
  updateQmAtoms();
  $<HTMLIFrameElement>("previewFrame").contentWindow?.postMessage(
    { type: "oqp-qmmm-clear" }, window.location.origin);
  updateInpPreview();
});

// ---------- menu bar ----------
// A classic pull-down bar: click a title to open it, hover to walk across the
// open menus, click anywhere else (or Esc) to close.
const menuBar = document.querySelector<HTMLElement>(".mbar")!;
const menuNote = $<HTMLDivElement>("menuNote");

function closeMenus(): void {
  menuBar.querySelectorAll(".m.open").forEach((m) => m.classList.remove("open"));
}

menuBar.querySelectorAll<HTMLElement>(".m").forEach((group) => {
  group.querySelector(".mtitle")!.addEventListener("click", (event) => {
    event.stopPropagation();
    const wasOpen = group.classList.contains("open");
    closeMenus();
    if (!wasOpen) group.classList.add("open");
  });
  group.addEventListener("mouseenter", () => {
    if (menuBar.querySelector(".m.open")) {
      closeMenus();
      group.classList.add("open");
    }
  });
});
document.addEventListener("click", closeMenus);
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape") closeMenus();
});

// The Tauri webview blocks window.open, so the backend hands the URL to the
// system browser; the direct call stays as the fallback for a plain browser.
function openExternally(url: string): void {
  fetch("/api/open-external", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ url }),
  })
    .then((response) => {
      if (!response.ok) window.open(url, "_blank", "noopener");
    })
    .catch(() => { window.open(url, "_blank", "noopener"); });
}

function download(name: string, text: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  URL.revokeObjectURL(url);
}

const COMMANDS: Record<string, () => void> = {
  "open-file": () => $<HTMLInputElement>("fileInput").click(),
  "open-results": () => { void openExistingResults(); },
  workspace: () => { void openWorkspaceSettings(); },
  pubchem: () => { showTab("builder"); $<HTMLInputElement>("pubchemName").focus(); },
  "save-xyz": () => {
    const atoms = parseAtoms(xyzArea.value);
    download("molecule.xyz", `${atoms.length}\nOQP Studio\n${atomsToText(atoms)}\n`);
  },
  "save-oqp": () => download("job.oqp", $<HTMLTextAreaElement>("input").value),
  "clear-geom": () => {
    xyzArea.value = "";
    invalidateSymmetry("Geometry cleared; analyze symmetry again.");
    builderStatus.textContent = "geometry cleared";
  },
  sketch: () => $<HTMLButtonElement>("sketchBtn").click(),
  preview: () => { void updatePreview(); },
  generate: () => $<HTMLButtonElement>("generate").click(),
  run: () => $<HTMLButtonElement>("runBtn").click(),
  "refresh-jobs": () => { void refreshJobs(); },
  "toggle-nav": () => setNavVisible(document.body.classList.contains("nav-hidden")),
  "reset-layout": resetLayout,
  docs: () => openExternally("https://docs.openqp.org"),
  network: () => { void openNetworkSettings(); },
  engine: () => { void installEngine(); },
  "tab-builder": () => showTab("builder"),
  "tab-method": () => showTab("method"),
  "tab-run": () => showTab("run"),
  "tab-analysis": () => showTab("analysis"),
  "tab-art": () => showTab("art"),
};

menuBar.addEventListener("click", (event) => {
  const item = (event.target as HTMLElement).closest<HTMLElement>("[data-cmd]");
  if (!item) return;
  closeMenus();
  COMMANDS[item.dataset.cmd!]?.();
});

// Help entries keep the menu open so their result stays readable.
$<HTMLElement>("menuNote").addEventListener("click", (event) => event.stopPropagation());
for (const id of ["updateBtn", "releasesBtn", "aboutBtn"]) {
  $<HTMLElement>(id).addEventListener("click", (event) => event.stopPropagation());
}

$<HTMLButtonElement>("releasesBtn").addEventListener("click", () =>
  openExternally("https://github.com/Open-Quantum-Platform/oqp-studio/releases/latest"));

const COPYRIGHT = "© 2026 Open Quantum, Inc. All rights reserved.";

$<HTMLButtonElement>("aboutBtn").addEventListener("click", async () => {
  const health = await (await fetch("/api/health")).json();
  menuNote.innerHTML =
    `<div><strong>OQP Studio ${escapeHtml(health.version)}</strong></div>` +
    `<div>A graphical interface for the Open Quantum Platform.</div>` +
    `<div style="margin-top:.4rem">${COPYRIGHT}</div>`;
});

$<HTMLButtonElement>("updateBtn").addEventListener("click", async () => {
  menuNote.textContent = "checking…";
  try {
    const info = await (await fetch("/api/update-check")).json();
    if (info.available) {
      menuNote.innerHTML =
        `Version ${escapeHtml(info.latest)} is available (you have ${escapeHtml(info.current)}).<br />` +
        `<button class="ghost" id="getUpdate" style="margin-top:.4rem;padding:.3rem .8rem">` +
        `Download and install</button> ` +
        `<a href="#" id="openReleases" style="color:var(--accent)">open the release page</a>`;
      $<HTMLAnchorElement>("openReleases").addEventListener("click", (event) => {
        event.preventDefault();
        openExternally(info.url);
      });
      $<HTMLButtonElement>("getUpdate").addEventListener("click", () => {
        void installUpdate();
      });
    } else if (info.latest) {
      menuNote.textContent = `${info.current} is the latest version.`;
    } else {
      menuNote.textContent = info.detail ?? "could not check for updates";
    }
  } catch {
    menuNote.textContent = "could not check for updates";
  }
});

// ---------- network settings ----------
// A TLS-inspecting proxy is the usual reason a PubChem lookup fails, and the
// user needs a way to hand the app their network's root certificate.
const networkOverlay = $<HTMLDivElement>("networkOverlay");
const caBundle = $<HTMLInputElement>("caBundle");
const insecureSsl = $<HTMLInputElement>("insecureSsl");
const networkNote = $<HTMLSpanElement>("networkNote");

type NetworkStatus = {
  system_trust_store: boolean; ca_bundle: string; ca_bundle_found: boolean;
  insecure: boolean; settings_path: string;
};

function renderNetworkStatus(status: NetworkStatus): void {
  const lines: [string, string][] = [
    ["System trust store",
     status.system_trust_store ? "in use" : "unavailable — falling back to the bundled CA list"],
    ["Root certificate file",
     status.ca_bundle ? (status.ca_bundle_found ? status.ca_bundle : `${status.ca_bundle} (missing)`) : "none"],
    ["Certificate checks", status.insecure ? "skipped (you switched this on)" : "on"],
  ];
  $<HTMLDivElement>("networkStatus").innerHTML = rows(lines);
  caBundle.value = status.ca_bundle;
  insecureSsl.checked = status.insecure;
}

// ---------- results folder ----------
// Runs write their outputs here. That is the user's data, so the folder is
// the user's choice; the app only picks a sensible one until they say.
const workspaceOverlay = $<HTMLDivElement>("workspaceOverlay");
const workspaceDir = $<HTMLInputElement>("workspaceDir");
const workspaceNote = $<HTMLSpanElement>("workspaceNote");

interface Workspace {
  jobs_dir: string;
  active: string;
  default: string;
  overridden: boolean;
}

function renderWorkspace(status: Workspace): void {
  workspaceDir.value = status.jobs_dir;
  const chosen = status.jobs_dir
    ? "your choice"
    : status.overridden
      ? "set by OQP_STUDIO_JOBS"
      : "chosen by the app";
  const summary = $<HTMLDivElement>("workspaceStatus");
  summary.replaceChildren();
  const active = document.createElement("div");
  active.append("Writing to ");
  const path = document.createElement("strong");
  path.textContent = status.active;
  active.appendChild(path);
  const source = document.createElement("div");
  source.className = "hint";
  source.textContent = chosen;
  summary.append(active, source);
  $<HTMLSpanElement>("executionOutputDir").textContent = status.active;
}

async function openWorkspaceSettings(): Promise<void> {
  workspaceNote.textContent = "";
  workspaceOverlay.style.display = "block";
  try {
    renderWorkspace(await (await fetch("/api/workspace")).json());
  } catch {
    workspaceNote.textContent = "could not read the current folder";
  }
}

$<HTMLButtonElement>("workspaceClose").addEventListener("click", () => {
  workspaceOverlay.style.display = "none";
});

$<HTMLButtonElement>("workspaceSave").addEventListener("click", async () => {
  workspaceNote.textContent = "saving…";
  const response = await fetch("/api/workspace", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ jobs_dir: workspaceDir.value.trim() }),
  });
  const data = await response.json().catch(() => null);
  if (response.ok) {
    renderWorkspace(data);
    workspaceNote.textContent = "saved";
    await refreshJobs();
  } else {
    workspaceNote.textContent = data?.detail ?? "could not save";
  }
});

$<HTMLButtonElement>("workspaceReveal").addEventListener("click", async () => {
  const response = await fetch("/api/reveal-results", { method: "POST" });
  if (!response.ok) {
    workspaceNote.textContent =
      (await response.json().catch(() => null))?.detail ?? "could not open the folder";
  }
});

async function openNetworkSettings(): Promise<void> {
  networkNote.textContent = "";
  networkOverlay.style.display = "block";
  try {
    renderNetworkStatus(await (await fetch("/api/network")).json());
  } catch {
    networkNote.textContent = "could not read the current settings";
  }
}

$<HTMLButtonElement>("networkClose").addEventListener("click", () => {
  networkOverlay.style.display = "none";
});

// The file picker only yields a name, not a path, so the file is read and
// written next to the settings; typing a path stays available.
$<HTMLButtonElement>("caBrowse").addEventListener("click", () =>
  $<HTMLInputElement>("caFile").click());

$<HTMLInputElement>("caFile").addEventListener("change", async (event) => {
  const file = (event.target as HTMLInputElement).files?.[0];
  if (!file) return;
  const body = new FormData();
  body.append("file", file);
  const response = await fetch("/api/network/certificate", { method: "POST", body });
  if (response.ok) {
    const status = await response.json();
    renderNetworkStatus(status);
    networkNote.textContent = "certificate stored";
  } else {
    networkNote.textContent = "could not read that certificate";
  }
});

$<HTMLButtonElement>("networkSave").addEventListener("click", async () => {
  networkNote.textContent = "saving…";
  const response = await fetch("/api/network", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ ca_bundle: caBundle.value.trim(), insecure: insecureSsl.checked }),
  });
  if (response.ok) {
    renderNetworkStatus(await response.json());
    networkNote.textContent = "saved";
  } else {
    networkNote.textContent = (await response.json()).detail ?? "could not save";
  }
});

$<HTMLButtonElement>("networkTest").addEventListener("click", async () => {
  networkNote.textContent = "testing…";
  try {
    const result = await (await fetch("/api/network/test")).json();
    networkNote.textContent = result.detail;
  } catch {
    networkNote.textContent = "the test could not run";
  }
});

// One click installs: the backend downloads this platform's installer and,
// on macOS, stages a swap that runs as soon as the app quits.
async function installUpdate(): Promise<void> {
  menuNote.textContent = "starting…";
  try {
    const started = await fetch("/api/update/install", { method: "POST" });
    if (!started.ok) {
      menuNote.textContent = (await started.json()).detail ?? "could not start the update";
      return;
    }
  } catch {
    menuNote.textContent = "could not start the update";
    return;
  }

  // Poll until it finishes; the menu stays open so the progress is visible.
  for (;;) {
    await new Promise((resolve) => setTimeout(resolve, 700));
    let state: { status: string; percent: number; detail: string };
    try {
      state = await (await fetch("/api/update/status")).json();
    } catch {
      menuNote.textContent = "lost contact with the backend";
      return;
    }
    if (state.status === "downloading") {
      menuNote.textContent = `downloading… ${state.percent}%  (${state.detail})`;
    } else if (state.status === "installing") {
      menuNote.textContent = `installing… ${state.detail}`;
    } else if (state.status === "ready") {
      menuNote.textContent =
        "The new version is staged. Quit OQP Studio and it will replace itself " +
        "and reopen automatically.";
      return;
    } else if (state.status === "failed") {
      menuNote.textContent = `update failed: ${state.detail}`;
      return;
    }
  }
}

// The compute engine ships as its own release archive — unzip and run, no
// Python or BLAS needed — so the app can fetch the one that matches and run
// jobs without the user assembling a toolchain. The same archive works from a
// terminal for anyone who would rather not use the GUI.
async function installEngine(): Promise<void> {
  menuNote.textContent = "checking…";
  try {
    const current = await (await fetch("/api/engine")).json();
    if (current.installed) {
      // An all-in-one installer already put one here; say so, rather than
      // leaving the user wondering whether they still have to fetch it.
      const where = current.source ? ` (${current.source})` : "";
      menuNote.textContent = `the engine is already available${where}: ${current.path}`;
      return;
    }
    const started = await fetch("/api/engine/install", { method: "POST" });
    if (!started.ok) {
      menuNote.textContent = (await started.json()).detail ?? "could not start the download";
      return;
    }
  } catch {
    menuNote.textContent = "could not reach the backend";
    return;
  }

  for (;;) {
    await new Promise((resolve) => setTimeout(resolve, 700));
    let state: { status: string; percent: number; detail: string };
    try {
      state = await (await fetch("/api/engine/status")).json();
    } catch {
      menuNote.textContent = "lost contact with the backend";
      return;
    }
    if (state.status === "downloading") {
      menuNote.textContent = `downloading the engine… ${state.percent}% (${state.detail})`;
    } else if (state.status === "ready") {
      menuNote.textContent = "the engine is installed — OpenQP is ready";
      await loadRunners();
      return;
    } else if (state.status === "failed") {
      menuNote.textContent = `engine install failed: ${state.detail}`;
      return;
    }
  }
}

// ---------- layout: hideable sidebar and draggable splitters ----------
const NAV_DEFAULT = 210;
const SPLIT_DEFAULT = "44%";
const shield = $<HTMLDivElement>("dragShield");

function setNavVisible(visible: boolean): void {
  document.body.classList.toggle("nav-hidden", !visible);
  const tick = document.querySelector<HTMLElement>('[data-tick="nav"]');
  if (tick) tick.style.visibility = visible ? "visible" : "hidden";
  localStorage.setItem("oqp.nav.visible", String(visible));
}

function setNavWidth(px: number): void {
  const width = Math.min(360, Math.max(120, px));
  document.documentElement.style.setProperty("--nav-w", `${width}px`);
  localStorage.setItem("oqp.nav.width", String(width));
}

$<HTMLButtonElement>("navToggle").addEventListener("click", () =>
  setNavVisible(document.body.classList.contains("nav-hidden")));

// A drag needs a full-window shield: without it the pointer entering an
// iframe (the Mol* viewer) would stop delivering mousemove to this document.
function beginDrag(handle: HTMLElement, onMove: (x: number) => void): void {
  handle.classList.add("dragging");
  shield.classList.add("on");
  const move = (event: MouseEvent) => onMove(event.clientX);
  const up = () => {
    handle.classList.remove("dragging");
    shield.classList.remove("on");
    window.removeEventListener("mousemove", move);
    window.removeEventListener("mouseup", up);
  };
  window.addEventListener("mousemove", move);
  window.addEventListener("mouseup", up);
}

$<HTMLDivElement>("navGrip").addEventListener("mousedown", (event) => {
  event.preventDefault();
  beginDrag(event.currentTarget as HTMLElement, (x) => setNavWidth(x));
});

// Every two-column panel gets a splitter between its controls and its viewer.
document.querySelectorAll<HTMLElement>(".cols").forEach((cols, index) => {
  const splitter = document.createElement("div");
  splitter.className = "splitter";
  splitter.title = "Drag to resize";
  cols.insertBefore(splitter, cols.children[1]);

  const key = `oqp.split.${index}`;
  const saved = localStorage.getItem(key);
  if (saved) cols.style.setProperty("--split", saved);

  splitter.addEventListener("mousedown", (event) => {
    event.preventDefault();
    beginDrag(splitter, (x) => {
      const box = cols.getBoundingClientRect();
      const width = Math.min(box.width - 260, Math.max(260, x - box.left));
      const value = `${Math.round(width)}px`;
      cols.style.setProperty("--split", value);
      localStorage.setItem(key, value);
    });
  });

  // Double-click restores the default proportion, as most editors do.
  splitter.addEventListener("dblclick", () => {
    cols.style.setProperty("--split", SPLIT_DEFAULT);
    localStorage.removeItem(key);
  });
});

function resetLayout(): void {
  setNavVisible(true);
  setNavWidth(NAV_DEFAULT);
  document.querySelectorAll<HTMLElement>(".cols").forEach((cols, index) => {
    cols.style.setProperty("--split", SPLIT_DEFAULT);
    localStorage.removeItem(`oqp.split.${index}`);
  });
}

setNavVisible(localStorage.getItem("oqp.nav.visible") !== "false");
setNavWidth(+(localStorage.getItem("oqp.nav.width") ?? NAV_DEFAULT));

// ---------- boot ----------
xyzArea.value = atomsToText(SAMPLES.water.atoms);
buildSampleList();
buildWorkflowSelect();
installRecipeControls();
syncPcmSolvent();
syncPcmReference();
selectWorkflow(WORKFLOWS[0]);
$<HTMLSpanElement>("copyright").textContent = COPYRIGHT;
async function startBackend(): Promise<void> {
  const health = $<HTMLSpanElement>("health");
  const stop = busyStatus(health, "Starting OQP Studio");
  const abort = new AbortController();
  const timeout = window.setTimeout(() => abort.abort(), 180_000);
  try {
    const response = await fetch("/api/health", { signal: abort.signal });
    if (!response.ok) throw new Error(`startup failed (${response.status})`);
    const data = await response.json();
    stop();
    health.textContent = `OQP Studio v${data.version} · ready`;
    void loadRunners();
    void loadRecipes().catch((error) => { $<HTMLDivElement>("recipeStatus").textContent = String(error); });
    void refreshExecutionHost();
    void updatePreview().catch(() => {});
  } catch (error) {
    stop();
    health.textContent = `Studio could not start; restart the app. ${String(error)}`;
  } finally {
    window.clearTimeout(timeout);
  }
}
void startBackend();
