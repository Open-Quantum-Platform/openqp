import { resolve } from "node:path";
import { defineConfig } from "vite";

export default defineConfig({
  define: {
    // ketcher packages read Node globals at runtime
    "process.env": {},
    global: "globalThis",
  },
  build: {
    // ketcher ships modules that mix require() with ESM
    commonjsOptions: { transformMixedEsModules: true },
    rollupOptions: {
      input: {
        main: resolve(import.meta.dirname, "index.html"),
        art: resolve(import.meta.dirname, "art.html"),
        sketcher: resolve(import.meta.dirname, "sketcher.html"),
      },
    },
  },
  server: {
    proxy: {
      "/api": process.env.OQP_STUDIO_BACKEND_URL ?? "http://127.0.0.1:8814",
    },
  },
});
