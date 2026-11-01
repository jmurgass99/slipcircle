import * as THREE from "three";
import van, { type State } from "vanjs-core";
import type { SliceGeometryMapping } from "../data-model";

export function getViewer3D({
  mapping,
  mobilizedShearKPa,
}: {
  mapping: State<SliceGeometryMapping>;
  mobilizedShearKPa: State<number[]>;
}): HTMLElement {
  const container = document.createElement("div");
  container.className = "sc-viewer";

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 500);
  camera.position.set(8, -26, 8);
  camera.up.set(0, 0, 1);
  camera.lookAt(8, 0, 4);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  container.appendChild(renderer.domElement);

  const group = new THREE.Group();
  scene.add(group);

  van.derive(() => {
    while (group.children.length > 0) group.remove(group.children[0]);
    const m = mapping.val;
    const tau = mobilizedShearKPa.val;
    const maxTau = Math.max(...tau, 1);

    for (const s of m.slices) {
      const ratio = (tau[s.index] ?? 0) / maxTau;
      const geom = new THREE.BoxGeometry(s.widthM * 0.92, 4.0, s.heightM);
      const col =
        s.stratum === "upper_stratum"
          ? new THREE.Color().setHSL(0.38 - ratio * 0.28, 0.82, 0.48)
          : new THREE.Color().setHSL(0.55 - ratio * 0.35, 0.85, 0.45);
      const box = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({ color: col }));
      box.position.set(s.xMidM, 0, s.yBaseM + 0.5 * s.heightM);
      group.add(box);
    }

    renderer.render(scene, camera);
  });

  const ro = new ResizeObserver(() => {
    const w = container.clientWidth || 700;
    const h = container.clientHeight || 500;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
    renderer.render(scene, camera);
  });
  ro.observe(container);

  return container;
}
