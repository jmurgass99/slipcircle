export interface SlopeSlice {
  index: number;
  xMidM: number;
  widthM: number;
  yTopM: number;
  yBaseM: number;
  heightM: number;
  alphaRad: number;
  stratum: "upper_stratum" | "lower_stratum";
}

export interface SliceGeometryMapping {
  stratumToSlices: Map<string, number[]>;
  reinforcementToSlices: Map<number, number[]>;
  slices: SlopeSlice[];
}

export interface LemResult {
  engine: "wasm-local" | "fastapi-python";
  fsBishop: number;
  fsFellenius: number;
  normalStressKPa: number[];
  mobilizedShearKPa: number[];
}
