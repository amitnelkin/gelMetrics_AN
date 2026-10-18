Rank-1 lambda / k0 (n = amplitudes fitted after the top-10% trim) per stage. * = k0 at a grid edge. Paper = MS.pdf, negative bursts only. SIM row: SimulationRunWithAnalysis.m, rng(1), true lambda = 2, k0 = 60.

| dataset           | sign   | paper           | baseline        | after_bugfixes   |
|:------------------|:-------|:----------------|:----------------|:-----------------|
| protein_RNA_1_1   | Neg    | 1 / 39 (n=91)   | 1 / 39 (n=82)   | 1 / 35 (n=82)    |
| protein_RNA_1_1   | Pos    |                 | 1 / 53 (n=64)   | 1 / 39 (n=64)    |
| protein_RNA_10_1  | Neg    | 1 / 63 (n=96)   | 1 / 62 (n=86)   | 1 / 62 (n=86)    |
| protein_RNA_10_1  | Pos    |                 | 1 / 139 (n=70)  | 1 / 88 (n=70)    |
| protein_RNA_100_1 | Neg    | 1 / 139 (n=103) | 1 / 139 (n=93)  | 1 / 110 (n=93)   |
| protein_RNA_100_1 | Pos    |                 | 1 / 145 (n=57)  | 1 / 121 (n=57)   |
| 4PP7              | Neg    | 1 / ? (n=63)    | 1 / 188 (n=57)  | 1 / 199 (n=57)   |
| 4PP7              | Pos    |                 | 1 / 379 (n=59)  | 1 / 212 (n=59)   |
| 4PP7-FUS          | Neg    | 2 / ? (n=210)   | 2 / 51 (n=189)  | 2 / 53 (n=189)   |
| 4PP7-FUS          | Pos    |                 | 1 / 73 (n=186)  | 1 / 79 (n=186)   |
| 30PP7             | Neg    | 1 / ? (n=112)   | 1 / 192 (n=101) | 1 / 204 (n=101)  |
| 30PP7             | Pos    |                 | 1 / 138 (n=78)  | 1 / 135 (n=78)   |
| 30PP7-FUS         | Neg    | 1 / ? (n=86)    | 1 / 97 (n=77)   | 1 / 118 (n=77)   |
| 30PP7-FUS         | Pos    |                 | 1 / 117 (n=77)  | 1 / 163 (n=77)   |
| 5Qb               | Neg    | 4 / 90 (n=134)  | 4 / 90 (n=121)  | 3 / 116 (n=121)  |
| 5Qb               | Pos    |                 | 1 / 270 (n=140) | 1 / 322 (n=140)  |
| 10Qb              | Neg    | 4 / 110 (n=181) | 4 / 112 (n=163) | 5 / 103 (n=163)  |
| 10Qb              | Pos    |                 | 3 / 181 (n=159) | 4 / 148 (n=159)  |
| SIM(lam=2,k0=60)  | Neg    | 4 / ?           | 2 / 63 (n=1845) | 4 / 41* (n=1663) |
| SIM(lam=2,k0=60)  | Pos    |                 | 2 / 72 (n=1688) | 1 / 107 (n=1425) |
