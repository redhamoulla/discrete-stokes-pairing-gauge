"""Reproducible finite-dimensional Stokes checks for several trial spaces."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import qr


def matrix_rank(a: np.ndarray, rtol: float = 1e-10) -> int:
    if a.size == 0:
        return 0
    singular_values = np.linalg.svd(a, compute_uv=False)
    if not singular_values.size or not singular_values[0]:
        return 0
    return int(np.sum(singular_values > rtol * max(1.0, singular_values[0])))


def synthesize(
    name,
    basis,
    derivative,
    a: float = 0.0,
    b: float = 1.0,
    quadrature_order: int = 240,
):
    xg, wg = leggauss(quadrature_order)
    x = 0.5 * (b - a) * xg + 0.5 * (a + b)
    w = 0.5 * (b - a) * wg

    phi = np.column_stack([function(x) for function in basis])
    dphi = np.column_stack([function(x) for function in derivative])
    effort_dimension = phi.shape[1]

    # Select a numerically independent basis of d(span(phi)).
    _, r_matrix, pivots = qr(dphi, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(r_matrix))
    tolerance = 1e-11 * (diagonal[0] if diagonal.size else 1.0)
    flow_dimension = int(np.sum(diagonal > tolerance))
    selected = list(map(int, pivots[:flow_dimension]))
    psi = dphi[:, selected] if flow_dimension else np.zeros((len(x), 0))
    derivative_matrix = (
        np.linalg.lstsq(psi, dphi, rcond=None)[0]
        if flow_dimension
        else np.zeros((0, effort_dimension))
    )

    pairing = phi.T @ (w[:, None] * psi)
    trace_a = np.array([[function(np.array([a]))[0] for function in basis]])
    trace_b = np.array([[function(np.array([b]))[0] for function in basis]])
    trace = np.vstack([trace_a, trace_b])
    boundary = trace.T @ np.diag([-1.0, 1.0]) @ trace
    defect = pairing @ derivative_matrix + derivative_matrix.T @ pairing.T - boundary

    pairing_rank = matrix_rank(pairing)
    kernel_dimension = effort_dimension - pairing_rank

    # Find the smallest subset of endpoint traces completing row(pairing.T).
    completion = None
    for number_of_rows in range(3):
        for subset in itertools.combinations(range(2), number_of_rows):
            stack = np.vstack([pairing.T, trace[list(subset), :]])
            if matrix_rank(stack) == effort_dimension:
                completion = list(subset)
                break
        if completion is not None:
            break

    # Canonical balanced rank factorization pairing = Pe.T @ Pf.
    u, singular_values, vt = np.linalg.svd(pairing, full_matrices=False)
    threshold = (
        1e-11 * max(1.0, singular_values[0])
        if singular_values.size and singular_values[0]
        else 0.0
    )
    keep = singular_values > threshold
    square_roots = np.sqrt(singular_values[keep])
    effort_factor = square_roots[:, None] * u[:, keep].T
    flow_factor = square_roots[:, None] * vt[keep, :]
    factor_defect = effort_factor.T @ flow_factor - pairing

    return {
        "name": name,
        "effort_dim": effort_dimension,
        "derived_flow_dim": flow_dimension,
        "mass_rank": pairing_rank,
        "kernel_dim_MT": kernel_dimension,
        "trace_rank": matrix_rank(trace),
        "boundary_rank": matrix_rank(boundary),
        "minimal_trace_completion": completion,
        "stokes_relative_defect": float(
            np.linalg.norm(defect)
            / max(1.0, np.linalg.norm(boundary), np.linalg.norm(pairing @ derivative_matrix))
        ),
        "rank_factorization_relative_defect": float(
            np.linalg.norm(factor_defect) / max(1.0, np.linalg.norm(pairing))
        ),
        "mass_condition_nonzero": (
            float(singular_values[0] / singular_values[keep][-1])
            if np.any(keep)
            else None
        ),
    }


def constant(value: float = 1.0):
    return lambda x: np.full_like(x, value, dtype=float)


def build_cases():
    rng = np.random.default_rng(7)

    # A deliberately non-nodal, non-orthogonal polynomial basis.
    polynomial_coefficients = rng.normal(size=(5, 5))
    while abs(np.linalg.det(polynomial_coefficients)) < 1e-3:
        polynomial_coefficients = rng.normal(size=(5, 5))

    def polynomial_function(coefficients):
        return lambda x: sum(
            coefficients[j] * x**j for j in range(len(coefficients))
        )

    def polynomial_derivative(coefficients):
        return lambda x: sum(
            j * coefficients[j] * x ** (j - 1)
            for j in range(1, len(coefficients))
        )

    polynomial_basis = [
        polynomial_function(polynomial_coefficients[:, i]) for i in range(5)
    ]
    polynomial_derivatives = [
        polynomial_derivative(polynomial_coefficients[:, i]) for i in range(5)
    ]

    trigonometric_basis = [
        constant(),
        lambda x: np.sin(2 * np.pi * x),
        lambda x: np.cos(2 * np.pi * x),
        lambda x: np.sin(4 * np.pi * x),
        lambda x: np.cos(4 * np.pi * x),
    ]
    trigonometric_derivatives = [
        constant(0.0),
        lambda x: 2 * np.pi * np.cos(2 * np.pi * x),
        lambda x: -2 * np.pi * np.sin(2 * np.pi * x),
        lambda x: 4 * np.pi * np.cos(4 * np.pi * x),
        lambda x: -4 * np.pi * np.sin(4 * np.pi * x),
    ]

    exponential_rates = [0.0, -1.25, 0.7, 2.1, 3.4]
    exponential_basis = [
        lambda x, rate=rate: np.exp(rate * x) for rate in exponential_rates
    ]
    exponential_derivatives = [
        lambda x, rate=rate: rate * np.exp(rate * x) for rate in exponential_rates
    ]

    centers = np.array([0.08, 0.31, 0.57, 0.83])
    epsilon = 9.0
    rbf_basis = [constant()] + [
        lambda x, center=center: np.exp(-epsilon * (x - center) ** 2)
        for center in centers
    ]
    rbf_derivatives = [constant(0.0)] + [
        lambda x, center=center: -2
        * epsilon
        * (x - center)
        * np.exp(-epsilon * (x - center) ** 2)
        for center in centers
    ]

    tanh_a = np.array([0.9, 1.7, 2.8, 4.1])
    tanh_b = np.array([-0.4, 0.2, -1.1, 0.7])
    neural_basis = [constant()] + [
        lambda x, aa=aa, bb=bb: np.tanh(aa * x + bb)
        for aa, bb in zip(tanh_a, tanh_b)
    ]
    neural_derivatives = [constant(0.0)] + [
        lambda x, aa=aa, bb=bb: aa / np.cosh(aa * x + bb) ** 2
        for aa, bb in zip(tanh_a, tanh_b)
    ]

    closed_sine_basis = [lambda x: np.sin(np.pi * x)]
    closed_sine_derivatives = [lambda x: np.pi * np.cos(np.pi * x)]

    return [
        synthesize(
            "transformed_polynomials", polynomial_basis, polynomial_derivatives
        ),
        synthesize(
            "periodic_trigonometric",
            trigonometric_basis,
            trigonometric_derivatives,
        ),
        synthesize("exponentials", exponential_basis, exponential_derivatives),
        synthesize("gaussian_rbfs", rbf_basis, rbf_derivatives),
        synthesize("neural_tanh_features", neural_basis, neural_derivatives),
        synthesize(
            "boundary_vanishing_single_sine",
            closed_sine_basis,
            closed_sine_derivatives,
        ),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="Write the JSON results to this path instead of standard output.",
    )
    arguments = parser.parse_args()
    rendered = json.dumps(build_cases(), indent=2) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

