# Random Matrices: Tracy-Widom Distribution of the Largest Eigenvalue

This project numerically verifies that the largest eigenvalue of large random symmetric
matrices follows the **Tracy-Widom distribution**, as predicted by random matrix theory,
rather than a Gaussian distribution.

## Method

Random symmetric matrices (off-diagonal ~ N(0,1), diagonal ~ N(0,2)) are generated in
parallel using an **MPI manager-worker model**: a manager distributes random seeds to
workers, each of which builds one matrix and computes its largest eigenvalue via **power
iteration**, falling back to a full **Schur decomposition** (LAPACK `DGEEV`) when power
iteration fails to converge. The resulting sample distribution is standardized and compared
against the theoretical Tracy-Widom GOE distribution using curve fitting and a
**Kolmogorov-Smirnov test**.

## Results

For n = 1500 matrices, the KS statistic against Tracy-Widom (≈0.009) is notably smaller
than against a matched Gaussian (≈0.015), confirming the largest eigenvalue follows
Tracy-Widom rather than Gaussian statistics, as theory predicts. Parallel scaling was also
studied: speedup is close to linear for large matrices, but MPI overhead limits efficiency
for small ones.

See `report.pdf` for full derivations, pseudocode, and plots.

## Contents

- `code/` — MPI source implementing the manager-worker eigenvalue computation
- `project3code.ipynb` — analysis notebook: distribution fitting, KS tests, plots
- `report.pdf` — full write-up

## Running

Open `project3code.ipynb` and run it to reproduce the plots and KS-test results. It expects
the eigenvalue output files (`eigs_*.txt`) produced by the MPI code in `code/`.
