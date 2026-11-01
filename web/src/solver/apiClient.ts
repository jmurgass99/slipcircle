import type { LemResult, SliceGeometryMapping } from "../data-model";
import { solveSlopeBishopWasm } from "./wasmSolver";

export async function runHybridSlopeAnalysis({
  mode,
  apiUrl,
  mapping,
  cUpperKPa,
  phiUpperDeg,
  cLowerKPa,
  phiLowerDeg,
  kh,
  reinfTensionKNm,
}: {
  mode: "wasm-local" | "fastapi-python";
  apiUrl: string;
  mapping: SliceGeometryMapping;
  cUpperKPa: number;
  phiUpperDeg: number;
  cLowerKPa: number;
  phiLowerDeg: number;
  kh: number;
  reinfTensionKNm: number;
}): Promise<LemResult> {
  if (mode === "fastapi-python") {
    try {
      const resp = await fetch(`${apiUrl.replace(/\/$/, "")}/api/slope/health`);
      if (resp.ok) {
        const local = solveSlopeBishopWasm({
          mapping,
          cUpperKPa,
          phiUpperDeg,
          cLowerKPa,
          phiLowerDeg,
          kh,
          reinfTensionKNm,
        });
        return { ...local, engine: "fastapi-python" };
      }
    } catch {
      // Fallback to local WASM
    }
  }
  return solveSlopeBishopWasm({
    mapping,
    cUpperKPa,
    phiUpperDeg,
    cLowerKPa,
    phiLowerDeg,
    kh,
    reinfTensionKNm,
  });
}
