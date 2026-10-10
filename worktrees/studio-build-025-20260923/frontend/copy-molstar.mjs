// Copies the Mol* standalone viewer into public/ so the Builder 3D page can
// load it without bundling molstar into the main app chunk.
import { copyFileSync, mkdirSync } from "node:fs";

mkdirSync("public/molstar", { recursive: true });
for (const file of ["molstar.js", "molstar.css"]) {
  copyFileSync(`node_modules/molstar/build/viewer/${file}`, `public/molstar/${file}`);
}
console.log("molstar viewer copied to public/molstar/");
