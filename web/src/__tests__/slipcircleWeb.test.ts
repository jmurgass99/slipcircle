import { describe, expect, it } from "vitest";
import { getSliceMesh } from "../mesh/getMesh";
import { initLemWasm, solveSlopeBishopWasm } from "../solver/wasmSolver";

describe("SlipCircle Web — SliceGeometryMapping + lemKernel.wasm", () => {
  it("maps strata and reinforcements across vertical slices", () => {
    const m = getSliceMesh(10, 45, 4, 14, 15, 25);
    expect(m.slices.length).toBe(25);
    const total =
      (m.stratumToSlices.get("upper_stratum")?.length ?? 0) +
      (m.stratumToSlices.get("lower_stratum")?.length ?? 0);
    expect(total).toBe(25);
  });

  it("solves Bishop Simplified FS via WebAssembly", async () => {
    await initLemWasm();
    const m = getSliceMesh(10, 45, 4, 14, 15, 20);
    const res = solveSlopeBishopWasm({
      mapping: m,
      cUpperKPa: 18,
      phiUpperDeg: 28,
      cLowerKPa: 24,
      phiLowerDeg: 32,
      kh: 0.08,
    });
    expect(res.fsBishop).toBeGreaterThan(0.8);
    expect(res.normalStressKPa.length).toBe(20);
  });
});
