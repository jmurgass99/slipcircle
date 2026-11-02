import type { LemResult, SliceGeometryMapping } from "../data-model.ts";

const LEM_WASM_BASE64 =
  "AGFzbQEAAAABCQFgBHx8fHwBfAMDAgAAByMCDGJpc2hvcE1BbHBoYQAAEHNsaWNlRHJpdmluZ1Rlcm0AAQogAg0AIAAgASACoiADo6ALEAAgACABoiAAIAKiIAOioAs=";

export interface LemWasmExports {
  bishopMAlpha: (cosA: number, sinA: number, tanPhi: number, fs: number) => number;
  sliceDrivingTerm: (w: number, sinA: number, kh: number, arm: number) => number;
}

let cached: LemWasmExports | null = null;

function decode(b64: string): Uint8Array {
  if (typeof Buffer !== "undefined") return new Uint8Array(Buffer.from(b64, "base64"));
  const s = atob(b64);
  const out = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}

export async function initLemWasm(): Promise<LemWasmExports> {
  if (cached) return cached;
  const bytes = decode(LEM_WASM_BASE64);
  const { instance } = await WebAssembly.instantiate(bytes.buffer as ArrayBuffer);
  cached = instance.exports as unknown as LemWasmExports;
  return cached;
}

export function getLemWasmSync(): LemWasmExports {
  if (cached) return cached;
  const bytes = decode(LEM_WASM_BASE64);
  const mod = new WebAssembly.Module(bytes.buffer as ArrayBuffer);
  const instance = new WebAssembly.Instance(mod);
  cached = instance.exports as unknown as LemWasmExports;
  return cached;
}

export function solveSlopeBishopWasm({
  mapping,
  cUpperKPa,
  phiUpperDeg,
  cLowerKPa,
  phiLowerDeg,
  kh = 0.0,
  reinfTensionKNm = 25.0,
}: {
  mapping: SliceGeometryMapping;
  cUpperKPa: number;
  phiUpperDeg: number;
  cLowerKPa: number;
  phiLowerDeg: number;
  kh?: number;
  reinfTensionKNm?: number;
}): LemResult {
  const wasm = getLemWasmSync();
  const gammaSoil = 19.0; // kN/m^3
  const tanPhiU = Math.tan((phiUpperDeg * Math.PI) / 180);
  const tanPhiL = Math.tan((phiLowerDeg * Math.PI) / 180);

  let drivingSum = 0;
  let felleniusNum = 0;

  for (const s of mapping.slices) {
    const w = gammaSoil * s.widthM * s.heightM;
    const cosA = Math.cos(s.alphaRad);
    const sinA = Math.sin(s.alphaRad);
    const lBase = s.widthM / Math.max(cosA, 0.15);
    const c = s.stratum === "upper_stratum" ? cUpperKPa : cLowerKPa;
    const tanPhi = s.stratum === "upper_stratum" ? tanPhiU : tanPhiL;

    drivingSum += wasm.sliceDrivingTerm(w, sinA, kh, 0.45);
    felleniusNum += c * lBase + Math.max(0, w * cosA) * tanPhi;
  }

  let reinfResist = 0;
  for (const [, sliceList] of mapping.reinforcementToSlices.entries()) {
    if (sliceList.length > 0) reinfResist += reinfTensionKNm * 0.85;
  }

  const denom = Math.max(drivingSum, 1.0);
  const fsFellenius = (felleniusNum + reinfResist) / denom;
  let fs = Math.max(0.6, fsFellenius);

  for (let iter = 0; iter < 15; iter++) {
    let num = reinfResist;
    for (const s of mapping.slices) {
      const w = gammaSoil * s.widthM * s.heightM;
      const cosA = Math.cos(s.alphaRad);
      const sinA = Math.sin(s.alphaRad);
      const c = s.stratum === "upper_stratum" ? cUpperKPa : cLowerKPa;
      const tanPhi = s.stratum === "upper_stratum" ? tanPhiU : tanPhiL;
      const mAlpha = Math.max(0.2, wasm.bishopMAlpha(cosA, sinA, tanPhi, fs));
      num += (c * s.widthM + w * tanPhi) / mAlpha;
    }
    const nextFs = num / denom;
    if (Math.abs(nextFs - fs) < 1e-4) {
      fs = nextFs;
      break;
    }
    fs = nextFs;
  }

  const normalStressKPa: number[] = [];
  const mobilizedShearKPa: number[] = [];
  for (const s of mapping.slices) {
    const w = gammaSoil * s.widthM * s.heightM;
    const cosA = Math.cos(s.alphaRad);
    const sinA = Math.sin(s.alphaRad);
    const c = s.stratum === "upper_stratum" ? cUpperKPa : cLowerKPa;
    const tanPhi = s.stratum === "upper_stratum" ? tanPhiU : tanPhiL;
    const mAlpha = Math.max(0.2, wasm.bishopMAlpha(cosA, sinA, tanPhi, fs));
    const sigmaN = (w / s.widthM - (c * sinA) / fs) / mAlpha;
    const tauMob = (c + Math.max(0, sigmaN) * tanPhi) / fs;
    normalStressKPa.push(Number(Math.max(0, sigmaN).toFixed(2)));
    mobilizedShearKPa.push(Number(Math.max(0, tauMob).toFixed(2)));
  }

  return {
    engine: "wasm-local",
    fsBishop: Number(fs.toFixed(3)),
    fsFellenius: Number(fsFellenius.toFixed(3)),
    normalStressKPa,
    mobilizedShearKPa,
  };
}
