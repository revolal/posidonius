# Kaula tides: equations and implementation

`kaula_tides.pdf` derives, step by step, the tidal force and the spin torque of the Kaula model
(`TidalModel::Kaula`) and describes their implementation in Posidonius. It covers the tide raised in a
planet by its star and the tide raised in the star by each planet, at degree two, for an arbitrary
rheology given by a tabulated complex Love number.

The equation numbers quoted in the code (`src/effects/tides/kaula/kaula_tide.rs`, `kaula.rs`, `tides.rs`)
are those of this document.

## How to read it

| Goal | Chapters |
|---|---|
| Use the model | 8 (equations, mean rates, truncation, limits), 9 (tide in the star), 10 (inputs, messages, known limitations) |
| Check the code | 8.9 (algorithm) and the table of Section 10.1 (equation → Rust function) |
| Follow the derivation | 2 to 7, then 8.3 (from the general gradient to the implemented form) |
| Conventions and notation | Appendix A |

## Files

| File | Content |
|---|---|
| `kaula_tides.pdf` | the document |
| `G_exact_fractions.pdf` | exact rational coefficients of the eccentricity functions printed in Tables C.3–C.4 |
| `data/G_exact_l2.txt` | complete series of G_20q and G_21q, \|q\| ≤ 30, to e^40, as exact fractions (one line `l p q n c_n` per coefficient) |
| `data/G_series.json` | the same, machine-readable, also for degree 3 |
| `src/` | LaTeX sources |

## Building the PDFs

```
cd src
latexmk -pdf kaula_tides.tex
cd annex_G_exact
latexmk -pdf G_exact_fractions.tex
```

The packages needed are those of a standard TeX Live installation (`amsmath`, `mathtools`, `natbib`,
`longtable`, `booktabs`, `xurl`, `hyperref`). The tables `src/chapters/tables_F.tex`, `tables_G.tex` and
`src/annex_G_exact/G_exact_terms.tex` are generated from the exact series and should not be edited by hand.

## Verification

Every formula was checked numerically, independently of its derivation, against the exact linear
response of the deformed body computed without Kaula's expansion. The Rust implementation is checked by
the tests listed in Section 10.9 (`tests/test_kaula_*.rs`), against reference values computed by an
independent Python implementation (`tests/data/kaula_planet/`, `tests/data/kaula_stellar/`).
