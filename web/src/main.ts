import van from "vanjs-core";
import type { SliceGeometryMapping } from "./data-model";
import { getSliceMesh } from "./mesh/getMesh";
import { runHybridSlopeAnalysis } from "./solver/apiClient";
import { initLemWasm } from "./solver/wasmSolver";
import { getSliceCharts } from "./ui/getSliceCharts";
import { getViewer3D } from "./viewer/getViewer3D";

const { div, h2, h3, input, label, span } = van.tags;

await initLemWasm();

const heightM = van.state(10);
const slopeAngleDeg = van.state(45);
const nSlices = van.state(20);
const cUpperKPa = van.state(18);
const phiUpperDeg = van.state(28);
const kh = van.state(0.1);
const reinfTensionKNm = van.state(25);

const mapping = van.state<SliceGeometryMapping>(
  getSliceMesh(heightM.val, slopeAngleDeg.val, 4, 14, 15, nSlices.val),
);
const fsBishop = van.state(1.4);
const fsFellenius = van.state(1.25);
const normalStressKPa = van.state<number[]>([]);
const mobilizedShearKPa = van.state<number[]>([]);

van.derive(() => {
  mapping.val = getSliceMesh(heightM.val, slopeAngleDeg.val, 4, 14, 15, nSlices.val);
});

van.derive(() => {
  runHybridSlopeAnalysis({
    mode: "wasm-local",
    apiUrl: "http://127.0.0.1:8000",
    mapping: mapping.val,
    cUpperKPa: cUpperKPa.val,
    phiUpperDeg: phiUpperDeg.val,
    cLowerKPa: cUpperKPa.val * 1.35,
    phiLowerDeg: phiUpperDeg.val + 4,
    kh: kh.val,
    reinfTensionKNm: reinfTensionKNm.val,
  }).then((res) => {
    fsBishop.val = res.fsBishop;
    fsFellenius.val = res.fsFellenius;
    normalStressKPa.val = res.normalStressKPa;
    mobilizedShearKPa.val = res.mobilizedShearKPa;
  });
});

const left = div(
  { class: "sc-panel" },
  h2({ style: "font-size: 1.05rem;" }, "SlipCircle 3D — Dovelas & Sismo"),
  div(
    { class: "sc-card" },
    h3("Malla de Dovelas (stratumToSlices / reinforcementToSlices)"),
    div(
      { class: "sc-row" },
      label("Número de dovelas N:"),
      input({
        type: "number",
        min: "6",
        max: "60",
        value: () => String(nSlices.val),
        oninput: (e: Event) => (nSlices.val = Math.max(6, Number((e.target as HTMLInputElement).value) || 20)),
      }),
    ),
    div(
      { class: "sc-row" },
      label("Ángulo talud β [°]:"),
      input({
        type: "number",
        min: "20",
        max: "70",
        value: () => String(slopeAngleDeg.val),
        oninput: (e: Event) => (slopeAngleDeg.val = Math.max(20, Number((e.target as HTMLInputElement).value) || 45)),
      }),
    ),
    div(
      { class: "sc-row" },
      label("Cohesión c' [kPa]:"),
      input({
        type: "number",
        value: () => String(cUpperKPa.val),
        oninput: (e: Event) => (cUpperKPa.val = Math.max(0, Number((e.target as HTMLInputElement).value) || 18)),
      }),
    ),
    div(
      { class: "sc-row" },
      label("Fricción φ' [°]:"),
      input({
        type: "number",
        value: () => String(phiUpperDeg.val),
        oninput: (e: Event) => (phiUpperDeg.val = Math.max(5, Number((e.target as HTMLInputElement).value) || 28)),
      }),
    ),
    div(
      { class: "sc-row" },
      label("Coeficiente sísmico kh:"),
      input({
        type: "number",
        step: "0.02",
        value: () => String(kh.val),
        oninput: (e: Event) => (kh.val = Math.max(0, Number((e.target as HTMLInputElement).value) || 0)),
      }),
    ),
    div(
      { class: "sc-row" },
      label("FS Bishop / Fellenius:"),
      span(() => `${fsBishop.val.toFixed(3)} / ${fsFellenius.val.toFixed(3)}`),
    ),
  ),
);

const right = div(
  { class: "sc-panel right" },
  getSliceCharts({
    normalStressKPa,
    mobilizedShearKPa,
  }),
);

document.body.appendChild(
  div(
    { class: "sc-layout" },
    left,
    getViewer3D({
      mapping,
      mobilizedShearKPa,
    }),
    right,
  ),
);
