# RiNNAL+

**RiNNAL+ — a MATLAB software for SDP-RLT and DNN relaxations of mixed-binary quadratic programming problems**

### [Di Hou](https://houdiopt.github.io), [Tianyun Tang](https://ttymath.github.io/tianyuntang.github.io/), [Kim-Chuan Toh](https://blog.nus.edu.sg/mattohkc/)

This software is designed to solve SDP-RLT and DNN relaxations of mixed-binary quadratic programming (MBQP) problems.

**Mixed-Binary Quadratic Programming (MBQP):**

```tex
minimize   x'Qx + 2c'x
subject to Ax = b                                (equality)
           Bx <= d                               (inequality)
           x_i ∈ {0,1}, for i ∈ B                (binary)
           x_i·x_j = 0, for (i, j) ∈ E           (complementarity)
           x ≥ 0                                 (nonnegative)
```

**SDP-RLT Relaxation:**

```tex
minimize   ⟨Q, X⟩ + 2c'x
subject to Ax = b
           A X = b x'
           diag(X)_B = x_B
           Bx <= d         
           B X <= d x'
           B X B' - B x d' - d x' B' + d d' >= 0
           X_{ij} = 0, for (i, j) ∈ E
           Y >= 0
           Y := [1, x'; x, X] ⪰ 0
```

---

## Citation

- Di Hou, Tianyun Tang, Kim-Chuan Toh, **RiNNAL+: a Riemannian ALM solver for SDP-RLT relaxations of mixed-binary quadratic programs,** [arXiv:2507.13776]([https://arxiv.org/abs/2507.13776])

```bibtex
@article{hou2025rinnal+,
  title={RiNNAL+: a Riemannian ALM solver for SDP-RLT relaxations of mixed-binary quadratic programs},
  author={Hou, Di and Tang, Tianyun and Toh, Kim-Chuan},
  journal={arXiv preprint arXiv:2507.13776},
  year={2025}
}
```

---

**Important note**

- The software is still under development. Thus it will invariably be buggy. We would appreciate your feedback and bug reports.
- This is research software and is not intended to be general-purpose at the moment.

## Setup

1. Set the path of RiNNAL+ in `setup_path.m`.
2. Set the path of SDPNAL+ in `setup_path.m` .
3. MEX files are necessary for SDPNAL+ but not for RiNNAL+.

## Examples

1. Run `runs_demo_BIQ.m` in MATLAB for a BIQ example.
2. Run `runs_demo_QKP.m` in MATLAB for a QKP example.

## Complete computational experiments

We use **demo scripts** for illustration examples and **runs_*.m** scripts for complete results.
To reproduce the full results reported in the paper, run the corresponding `runs_*.m` script listed below.
To generate the exact table TeX, see generators in `results/Summary/(specific problem class)`.

## Script ↔ Table/Figure mapping

The scripts below correspond to the main figures/tables in the paper.

- **Figure 1 (rank evolution)**: `runs_PGvsAdaptive.m`.
- **Figure 2 (performance profiles)**: run `plot_performance_profile.m`.
  By default it uses existing `.mat` files (`regenerate = false`).
  To recompute the `.mat` files, set `regenerate = true` at the top of the script.
- **Table 4 (tightness for BIQ)**: `runs_tightness_BIQ.m` (requires Gurobi).
- **Table 5 (BIQ)**: subset of Table 13 (use `runs_BIQ.m` for the full table).
- **Table 6 (Theta)**: subset of Table 14 (use `runs_theta.m` for the full table).
- **Table 7 (QKP)**: `runs_QKP.m`.
- **Table 8 (ccMSSC)**: data statistics only; no computation needed.
- **Table 9 (ccMSSC)**: subset of Table 15 (use `runs_ccMSSC.m` for the full table).
- **Table 10 (SStQP)**: `runs_SStQP.m`.
- **Table 11 (QMSTP)**: subset of Table 16 (use `runs_QMSTP.m` for the full table).
- **Table 12 (QMSTP)**: summary computed from Table 16 results (see `results/Summary/QMSTP/`).
- **Table 13 (BIQ)**: `runs_BIQ.m`.
- **Table 14 (Theta)**: `runs_theta.m`.
- **Table 15 (ccMSSC)**: `runs_ccMSSC.m`.
- **Table 16 (QMSTP)**: `runs_QMSTP.m`.
