# Rank-pair segregation: discussion draft

`main.tex` is the editable manuscript. Open this directory in the Codex LaTeX editor for a PDF preview. The standalone PDF is `output/pdf/main.pdf`.

The draft begins with segregation between income-rank slices across the full distribution and the need to separate sorting from monetary income inequality. It develops a common pairwise representation, derives Reardon's three measures as special cases, and introduces beta-product and power-gap weights as alternative evaluative specifications. Toy cases distinguish inequality from sorting and show how weights can reverse comparisons of residential arrangements. Empirical sections discuss interpolation, coherent density estimation, and an application strategy; the final section sets out unresolved mathematical and inferential questions. The motivation cites research on differential patterns across the income hierarchy, including existing work on middle-income groups.

No author or prospective collaborator has been assigned authorship in this draft.

The slice summaries include an exact decomposition into contributions from pairs of reporting bands, including comparisons within each band. Fixed rank weights give a corresponding decomposition of changes. A quartile example shows offsetting changes beneath an unchanged summary; the discussion distinguishes this accounting from the potentially signed contributions in a general pairwise representation of Reardon's indices.

## Reproduction

Build the PDF with an installed TeX distribution:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=output/pdf main.tex
```

Run the standalone numerical checks with base R:

```sh
Rscript check-identities.R
```

The checks compare threshold integration with triangular pair integration for all three variation functions and five continuous CDFs, verify normalizing constants, independently calculate the illustrative table, and check a population-mixture example. They also verify the slice-mixing formulas, translate three slice-weight specifications to pair contributions, and recover Reardon's indices from distribution-weighted slice mixing for smooth, contiguous, and interleaved neighborhood distributions. The latter cases exercise zero-density regions. Additional checks verify beta-product normalization, independently integrate the quarter-separation example, and verify the inequality-changing, rank-preserving income schedules. The checks support the identities; admissibility of arbitrary valuation functions and overlapping-neighborhood extensions remains a research question.

Generated compilation files are ignored by Git. The parent `papers` directory is excluded from R package builds.

Decomposition checks compare exact beta-weighted cell contributions with independent integration of the full surface, verify that refining reporting bands preserves contributions, and reproduce the quartile change table for all three index families.
