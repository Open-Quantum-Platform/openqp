// ChemDraw-style 2D sketcher page (Ketcher), embedded by the Builder tab.
// "Use structure" sends the molfile to the parent window, which asks the
// backend (RDKit) for 3D coordinates.

import { useRef, useState } from "react";
import { createRoot } from "react-dom/client";
import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import type { Ketcher } from "ketcher-core";
import "ketcher-react/dist/index.css";

const structServiceProvider = new StandaloneStructServiceProvider();

function SketcherApp() {
  const ketcherRef = useRef<Ketcher | null>(null);
  const [status, setStatus] = useState("");

  async function useStructure() {
    const ketcher = ketcherRef.current;
    if (!ketcher) return;
    const molfile = await ketcher.getMolfile();
    if (!molfile || !molfile.trim()) {
      setStatus("draw a structure first");
      return;
    }
    setStatus("sending…");
    window.parent.postMessage(
      { type: "oqp-sketch", molfile },
      window.location.origin,
    );
  }

  return (
    <>
      <div className="editor-host">
        <Editor
          staticResourcesUrl=""
          structServiceProvider={structServiceProvider}
          errorHandler={(message) => console.error("ketcher:", message)}
          onInit={(ketcher) => {
            ketcherRef.current = ketcher;
          }}
        />
      </div>
      <div className="bar">
        <span className="status">{status}</span>
        <button onClick={useStructure}>Use structure → 3D</button>
      </div>
    </>
  );
}

createRoot(document.getElementById("root")!).render(<SketcherApp />);
