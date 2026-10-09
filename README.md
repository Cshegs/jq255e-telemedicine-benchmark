# jq255e Telemedicine Benchmark

Reference implementation, benchmark scripts, raw data, and statistical audit for the paper:

> **A Genus-1 Reformulation of an α–β Modified Edwards-Type Curve for Secure Telemedicine Data Transmission**
> P. Adeola, A. T. Akinwale, B. B. Olofin, A. E. Babalola
> Department of Computing, Anchor University, Lagos, Nigeria
> Manuscript 1619

---

## Repository contents

| File / folder | Description |
|---|---|
| `jq255e.py` | Pure-Python reference implementation of the jq255e Jacobi-quartic curve |
| `benchmark.py` | Benchmarks jq255e, X25519, and P-256 across three plaintext sizes |
| `statistical_audit.py` | Full statistical analysis: Shapiro-Wilk, Wilcoxon, effect sizes, CIs |
| `requirements.txt` | Python package dependencies |
| `LICENSE` | MIT License |
| `data/paired_benchmark_raw.csv` | Raw timing observations (50 trials × 3 curves × 3 operations × 3 msg sizes) |
| `data/statistical_audit_report.txt` | Full audit report generated from the raw CSV |

---

## Curve summary

| Parameter | Value |
|---|---|
| Equation | E² = 8u⁴ + 1 (Jacobi quartic) |
| Field | F\_q, q = 2²⁵⁵ − 18651 |
| Group order | #E(F\_q) = 2r |
| Prime subgroup order | r = 2²⁵⁴ − 131528281291764213006042413802501683931 |
| Cofactor | h = 2 |
| Generator H (prime-order) | H = 2G, where G = (E=3, u=1) |
| Security level | ≈ 128 bits (classical) |
| Published specification | https://doubleodd.group/params-do255e.html |

---

## Installation

Python 3.10 or later is required.

```bash
git clone https://github.com/<your-username>/jq255e-telemedicine-benchmark.git
cd jq255e-telemedicine-benchmark
pip install -r requirements.txt
```

---

## Reproducing the results

### Step 1 — Run the statistical audit on the provided raw data

This reproduces all statistical results in the paper (Table 4, Table 5, Section 5.3)
without re-running the full benchmark:

```bash
python statistical_audit.py
```

The full audit report is printed to the console and saved to
`data/statistical_audit_report.txt`.

### Step 2 — Re-run the benchmark (optional)

This re-generates the raw timing CSV. Note that execution times will differ from
those reported in the manuscript because:

- The jq255e construction is implemented in pure Python using affine coordinates,
  whereas X25519 and P-256 use optimised C-library routines.
- Absolute timings are hardware- and OS-dependent.
- The manuscript results were obtained on Windows 11, Intel Core i5, 8 GB RAM,
  Python 3.12.2.

```bash
python benchmark.py
```

The script runs a warm-up phase (5 iterations) followed by 50 paired trials for
each combination of curve, operation, and message size. Key-generation timings
are printed to the console. All other results are saved to
`data/paired_benchmark_raw.csv`.

---

## Raw data format

`data/paired_benchmark_raw.csv` contains one row per measurement:

| Column | Description |
|---|---|
| `message_size_bytes` | Plaintext size: 13, 55, or 1000 bytes |
| `operation` | `ecdh`, `encryption`, or `decryption` |
| `curve` | `jq255e`, `X25519`, or `P-256` |
| `trial` | Trial index (0–49) |
| `time_us` | Execution time in microseconds (µs) |

All three curves share the same trial index, enabling paired comparison.

---

## Important limitations

The jq255e implementation in `jq255e.py` is a **reference implementation**
intended for research validation only. It is:

- **Not constant-time.** It is vulnerable to timing side-channel attacks.
- **Not complete.** The affine Jacobi-quartic addition formulas have an
  exceptional case (denominator D = 0) that is not handled. A complete,
  exception-free law requires projective or extended coordinates.
- **Not production-ready.** It must not be deployed in clinical or security-
  critical systems without a reviewed, constant-time implementation and
  independent cryptographic audit.

---

## Dependencies

| Package | Purpose |
|---|---|
| `cryptography` | X25519, P-256 ECDH, AES-256-GCM |
| `pandas` | CSV loading and manipulation |
| `scipy` | Shapiro-Wilk, Wilcoxon signed-rank, t-test |
| `numpy` | Numerical utilities |

---

## Citation

If you use this code or data in your work, please cite:

```
P. Adeola, A. T. Akinwale, B. B. Olofin, and A. E. Babalola,
"A Genus-1 Reformulation of an alpha-beta Modified Edwards-Type Curve
for Secure Telemedicine Data Transmission," Manuscript 1619.
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
