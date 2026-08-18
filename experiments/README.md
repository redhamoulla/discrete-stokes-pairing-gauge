# Numerical checks

`automatic_stokes_examples.py` tests the finite-dimensional construction on several non-nodal and non-orthogonal function families over `[0, 1]`.

For each trial space, the script:

1. selects an independent basis of the differentiated space;
2. constructs the mixed pairing matrix by high-order Gauss-Legendre integration;
3. forms the endpoint boundary matrix;
4. measures the residual of the discrete Stokes identity;
5. finds the smallest subset of physical endpoint traces that completes the rank;
6. verifies a balanced rank factorization of the mixed pairing.

The final single-sine example is intentionally boundary-vanishing. It records the obstruction that the Stokes identity may hold while physical endpoint traces fail to provide a nondegenerate open completion.

Run from the repository root:

```bash
make experiments
```

