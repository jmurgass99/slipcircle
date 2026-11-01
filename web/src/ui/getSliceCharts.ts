import van, { type State } from "vanjs-core";

export function getSliceCharts({
  normalStressKPa,
  mobilizedShearKPa,
}: {
  normalStressKPa: State<number[]>;
  mobilizedShearKPa: State<number[]>;
}): HTMLElement {
  const wrap = document.createElement("div");
  wrap.style.display = "flex";
  wrap.style.flexDirection = "column";
  wrap.style.gap = "12px";

  const cSig = document.createElement("canvas");
  cSig.className = "sc-chart";
  const cTau = document.createElement("canvas");
  cTau.className = "sc-chart";

  const makeCard = (t: string, c: HTMLCanvasElement) => {
    const d = document.createElement("div");
    d.className = "sc-card";
    const h = document.createElement("h3");
    h.textContent = t;
    d.append(h, c);
    return d;
  };

  wrap.append(
    makeCard("Tensión Normal Efectiva en la Base σ'n [kPa]", cSig),
    makeCard("Cortante Movilizado por Dovela τmob [kPa]", cTau),
  );

  const drawBars = (canvas: HTMLCanvasElement, vals: number[], color: string) => {
    const ctx = canvas.getContext("2d");
    if (!ctx || vals.length === 0) return;
    const w = (canvas.width = canvas.clientWidth || 320);
    const h = (canvas.height = canvas.clientHeight || 165);
    ctx.clearRect(0, 0, w, h);
    const maxV = Math.max(...vals, 1) * 1.1;
    const bw = (w - 30) / vals.length - 2;
    vals.forEach((v, i) => {
      const x = 18 + i * (bw + 2);
      const bh = ((h - 28) * v) / maxV;
      ctx.fillStyle = color;
      ctx.fillRect(x, h - 18 - bh, bw, bh);
    });
  };

  van.derive(() => {
    drawBars(cSig, normalStressKPa.val, "#38bdf8");
    drawBars(cTau, mobilizedShearKPa.val, "#10b981");
  });

  return wrap;
}
