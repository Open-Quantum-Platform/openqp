declare module "three/examples/jsm/libs/surfaceNet.js" {
  export function surfaceNet(
    dimensions: [number, number, number],
    potential: (x: number, y: number, z: number) => number,
    bounds?: [[number, number, number], [number, number, number]],
  ): {
    positions: [number, number, number][];
    cells: [number, number, number][];
  };
}
