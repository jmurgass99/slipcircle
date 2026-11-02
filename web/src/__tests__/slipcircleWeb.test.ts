import assert from "node:assert/strict";
import { getSliceMesh } from "../mesh/getMesh.ts";
import { initLemWasm, solveSlopeBishopWasm } from "../solver/wasmSolver.ts";

// 1. Slice geometry mapping
const m = getSliceMesh(10, 45, 4, 14, 15, 25);
assert.equal(m.slices.length, 25);
const total =
  (m.stratumToSlices.get("upper_stratum")?.length ?? 0) +
  (m.stratumToSlices.get("lower_stratum")?.length ?? 0);
assert.equal(total, 25);

// 2. Bishop Simplified in WASM
await initLemWasm();
const m20 = getSliceMesh(10, 45, 4, 14, 15, 20);
const res = solveSlopeBishopWasm({
  mapping: m20,
  cUpperKPa: 18,
  phiUpperDeg: 28,
  cLowerKPa: 24,
  phiLowerDeg: 32,
  kh: 0.08,
});

assert.ok(res.fsBishop > 0.8);
assert.equal(res.normalStressKPa.length, 20);

console.log(
  `SlipCircle WASM test OK: FS_Bishop=${res.fsBishop.toFixed(3)}, n_slices=${res.normalStressKPa.length}`
);
