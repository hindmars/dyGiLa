# dyGiLa
![dyGiLa Logo](https://github.com/dyGiLa/Images/blob/main/reposity-readme-logo-pic.png)

[![Home Page](https://img.shields.io/badge/Home-Website-blue?style=plastic&logo=Google-Chrome&logoColor=white)](https://dygila.github.io)
[![YouTube Channel](https://img.shields.io/badge/YouTube-Channel-red?style=plastic&logo=YouTube&logoColor=white)](https://www.youtube.com/@Superfluid-Field-Sim)
[![Releases](https://img.shields.io/github/v/release/dyGiLa/dyGiLa?style=plastic&logo=GitHub)](https://github.com/dyGiLa/dyGiLa/releases)
[![Contributors](https://img.shields.io/github/contributors/dyGiLa/dyGiLa?style=plastic&logo=GitHub)](https://github.com/dyGiLa/dyGiLa/graphs/contributors)


The dyGiLa project is an effort to apply algorithms of lattice field theory
and enabling GPU-CPU hybrid technologies with C++ fashion for efficient HPC simulation of
large-scale, 3-dimensional tensor-valued p-wave Ginzburg-Landau effective theory.

### Documentation

- [Dependencies and Compilation](https://dygila.github.io/install)

- [Benchmarks](https://dygila.github.io/benchmark)


### Phase marker and D<sub>B</sub> configuration

The parameter file is read in a fixed order, so these keys must appear exactly
where shown. **Configuration files written for v0.0.6 need these keys added;
without them dyGiLa stops with a missing-key error.** A complete example is
[`paras_conf/sim_config_hot-blob_DB.txt`](paras_conf/sim_config_hot-blob_DB.txt).

Directly after `pario_compute_phaseMarker`:

```
pario_compute_DB         yes        # D_B = sum_a (R_a - R_a^B)^2, distance from the B phase
invariants_gapA_min      0.2        # gap below this: phaseMarker 0 (normal), D_B masked
```

- `pario_compute_DB` (`yes`/`no`) computes D<sub>B</sub> = Σ<sub>a</sub>(R<sub>a</sub> − R<sub>a</sub><sup>B</sup>)²,
  the distance from the B phase in the space of the five invariants R<sub>1</sub>…R<sub>5</sub>.
- `invariants_gapA_min` is required if `pario_compute_phaseMarker` or `pario_compute_DB` is `yes`.
  The invariants are computed from the normalised order parameter, so they are just noise
  where the gap √tr(AA<sup>†</sup>) is small. Sites with a gap below this value
  (of the low-pass-filtered A if `useGaussianLP_filter` is `yes`) get phase marker 0 and a
  masked D<sub>B</sub>. Set it to 0 to reproduce the v0.0.6 phase marker. A value of 0.2
  separates the normal region cleanly in the hot-blob runs; 0.1 is within the noise there.

Directly before `remove_screen_annotations`:

```
do_DB_slice              yes        # needs pario_compute_DB yes
DBSlice_params_array     258.6, 258.6, 129.3, 0, 0, 1   # point x, y, z; normal x, y, z
DB_range_array           0.0, 2.34  # fixed colour range
DB_color_tables          "Jet"      # quoted, comma-separated; one image set per table
DB_mask_rgb              0.5, 0.5, 0.5   # colour of masked sites
```

The last four keys are read only if `do_DB_slice` is `yes`. Slice points are in
physical coordinates (lattice index × `dx`). Images are written to
`insitu/DB-slice/DB-slice_<table>_t-<cycle>.png` (non-alphanumeric characters in the
table name become `_`), so the run directory needs `insitu/DB-slice/`.

Phase-marker values and D<sub>B</sub> for the bulk phases:

| marker | phase | D<sub>B</sub> |
|---|---|---|
| 0 | normal (gap below `invariants_gapA_min`; drawn grey, not in the colour bar) | masked |
| 1 | unclassified | — |
| 2 | β | 2.0 |
| 3 | planar | 0.083 |
| 4 | γ | 1.667 |
| 5 | B | 0 |
| 6 | α | 1.0 |
| 7 | polar | 1.333 |
| 8 | bipolar | 1.083 |
| 9 | A | 2.0 |

`remove_screen_annotations` now controls the colour bar and title on every slice image,
including phase-marker slice 2. The bar sits across the top, at
`CBxMin`…`CBxMax`, `CByMin`…`CByMax`, so the slice may need a smaller `zoom1`/`zoom2`
to leave room for it.

### License

dyGiLa is licensed on under GPL-v2 license. 

### Contributing

We appreciate your willingness to contribute to the dyGiLa project. 