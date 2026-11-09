from __future__ import annotations

import pytest

from scripts.r16_projection_audit import reconstruct_projection


def test_monotone_source_q_is_not_pooled_and_pressure_is_reintegrated():
    mu = [2.0, 3.0, 4.0]
    n = [1.0, 2.25, 4.0]  # q = n/mu = [0.5, 0.75, 1.0]
    p = [1.0, 1.75, 3.0]
    eps = [0.1, 0.5, 1.0]

    out = reconstruct_projection(mu, n, p, eps)

    assert out["pav_block_count"] == 3
    assert out["pav_pooled_block_count"] == 0
    assert out["source_q_monotone"] is True
    assert out["projected_n"] == pytest.approx(n, rel=0.0, abs=1e-14)
    assert out["projected_p"] == pytest.approx([1.0, 31.0 / 12.0, 17.0 / 3.0], rel=0.0, abs=1e-14)
    assert out["projected_q"] == pytest.approx([0.5, 0.75, 1.0], rel=0.0, abs=1e-14)


def test_nonmonotone_source_q_is_pooled_but_projected_n_is_not_silently_preserved():
    mu = [1.0, 2.0, 3.0]
    n = [1.0, 0.4, 3.0]  # q = [1.0, 0.2, 1.0]
    p = [0.1, 0.5, 2.0]
    eps = [mu[i] * n[i] - p[i] for i in range(3)]

    out = reconstruct_projection(mu, n, p, eps)

    assert out["source_q_monotone"] is False
    assert out["source_q_violation_count"] == 1
    assert out["pav_block_count"] == 2
    assert out["pav_pooled_block_count"] == 1
    assert out["q_projected"][0] == pytest.approx(out["q_projected"][1])
    assert out["n_projected"][1] != pytest.approx(n[1])
    assert all(b >= a for a, b in zip(out["q_projected"], out["q_projected"][1:]))


def test_reconstruction_rejects_nonpositive_or_nonmonotone_source_knots():
    with pytest.raises(ValueError, match="strictly increasing"):
        reconstruct_projection([1.0, 1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [0.1, 0.2, 0.3])

    with pytest.raises(ValueError, match="positive"):
        reconstruct_projection([1.0, 2.0, 3.0], [1.0, 0.0, 3.0], [1.0, 2.0, 3.0], [0.1, 0.2, 0.3])

from scripts.r16_projection_audit import hermite_cell_diagnostics


def test_sparse_source_cell_can_fail_raw_hermite_convexity_while_projected_cell_passes():
    mu0, mu1 = 939.5626728207, 941.9599599398
    n0, n1 = 0.000269021, 0.0020000001
    p0_source, p1_source = 0.0004981910045986, 0.001960000017999996
    p0_projected, p1_projected = 0.0004979837658772507, 0.0032159718477015545

    raw = hermite_cell_diagnostics(mu0, p0_source, n0, mu1, p1_source, n1)
    projected = hermite_cell_diagnostics(mu0, p0_projected, n0, mu1, p1_projected, n1)

    assert raw["stable"] is False
    assert raw["min_d2p"] < 0.0
    assert projected["stable"] is True
    assert projected["causal"] is True
    assert projected["min_d2p"] > 0.0
    assert projected["min_causal_margin"] >= 0.0


def test_projection_exposes_cell_level_raw_and_projected_admissibility():
    out = reconstruct_projection(
        mu=[2.0, 3.0, 4.0],
        n_source=[1.0, 2.25, 4.0],
        p_source=[1.0, 1.75, 3.0],
        eps_source=[0.1, 0.5, 1.0],
    )

    assert len(out["source_cell_diagnostics"]) == 2
    assert len(out["projected_cell_diagnostics"]) == 2
    assert all(cell["stable"] and cell["causal"] for cell in out["projected_cell_diagnostics"])
    assert all(cell["min_d2p"] > 0.0 for cell in out["projected_cell_diagnostics"])
