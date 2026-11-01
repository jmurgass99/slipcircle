from __future__ import annotations

from slipcircle.geometry_mapping import build_slope_slice_mesh


def test_build_slope_slice_mesh_geometry_mapping():
    m = build_slope_slice_mesh(n_slices=24)
    assert len(m.slices) == 24
    total_mapped = len(m.stratum_to_slices["upper_stratum"]) + len(m.stratum_to_slices["lower_stratum"])
    assert total_mapped == 24
    assert 0 in m.reinforcement_to_slices
