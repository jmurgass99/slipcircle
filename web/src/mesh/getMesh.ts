import type { SliceGeometryMapping, SlopeSlice } from "../data-model";

export function getSliceMesh(
  heightM: number,
  slopeAngleDeg: number,
  cxM: number,
  cyM: number,
  radiusM: number,
  nSlices: number,
  layerBoundaryDepthM = 4.5,
): SliceGeometryMapping {
  const nS = Math.max(4, Math.floor(nSlices));
  const cotB = 1 / Math.tan((Math.max(15, Math.min(75, slopeAngleDeg)) * Math.PI) / 180);
  const xCrest = heightM * cotB;

  const xLeft = Math.max(0, cxM - Math.sqrt(Math.max(radiusM * radiusM - cyM * cyM, 1)));
  const xRight = cxM + Math.sqrt(Math.max(radiusM * radiusM - (cyM - heightM) * (cyM - heightM), 1));
  const dx = (xRight - xLeft) / nS;

  const stratumToSlices = new Map<string, number[]>([
    ["upper_stratum", []],
    ["lower_stratum", []],
  ]);
  const reinforcementToSlices = new Map<number, number[]>([
    [0, []],
    [1, []],
  ]);
  const reinfY = [0.3 * heightM, 0.6 * heightM];
  const slices: SlopeSlice[] = [];

  for (let i = 0; i < nS; i++) {
    const xm = xLeft + (i + 0.5) * dx;
    let yTop = heightM;
    if (xm <= 0) yTop = 0;
    else if (xm <= xCrest) yTop = (xm / Math.max(xCrest, 1e-6)) * heightM;

    const sq = Math.max(radiusM * radiusM - (xm - cxM) * (xm - cxM), 0.1);
    const yBase = cyM - Math.sqrt(sq);
    const hSlice = Math.max(0.15, yTop - yBase);
    const sinA = Math.max(-0.92, Math.min(0.92, (xm - cxM) / radiusM));
    const alphaRad = Math.asin(sinA);

    const localDepth = yTop - yBase;
    const stratum: "upper_stratum" | "lower_stratum" =
      localDepth <= layerBoundaryDepthM ? "upper_stratum" : "lower_stratum";
    stratumToSlices.get(stratum)!.push(i);

    reinfY.forEach((yr, rIdx) => {
      if (yBase <= yr && yr <= yTop) {
        reinforcementToSlices.get(rIdx)!.push(i);
      }
    });

    slices.push({
      index: i,
      xMidM: xm,
      widthM: dx,
      yTopM: yTop,
      yBaseM: yBase,
      heightM: hSlice,
      alphaRad,
      stratum,
    });
  }

  return { stratumToSlices, reinforcementToSlices, slices };
}
