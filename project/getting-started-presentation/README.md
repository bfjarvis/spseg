# Getting-started presentation

Self-contained Quarto Beamer source using the standard Madrid theme (16:9).
The R chunks embed the vignette's simulation, analysis, aggregation, parallel
calculations, and ggplot figures. Rendering executes them in order. Selected
syntax is displayed; setup, longer preparations, and figure code stay in the
source with echo/include controls to keep slides visual.

```sh
cd project/getting-started-presentation
quarto render spseg-slides.qmd
```

Requires Quarto, a LaTeX distribution with Beamer, and R with the current spseg,
sf, dplyr, tidyr, purrr, ggplot2, igraph, and knitr. Install the current package first.
No external figure-generation script or saved input figures are required.
Quarto/knitr manages temporary figure files during rendering.

- `spseg-slides.qmd`: executable slides and presenter notes.
- `_quarto.yml`: Beamer configuration and execution defaults.
- `spseg-slides.pdf`: rendered presentation.

The seeded example uses 2,000 people and 80 Voronoi polygons. Single-pattern
analyses use radius bands; parallel pattern comparisons use population kNN.
Hidden assertions check count conservation and parallel/sequential agreement.
The measure table follows the measure-definitions vignette, with footnoted
citations and explicit labels for package dissimilarity supplements.

The status examples include both observed micro-level categories and local shares
of category 3. The network example samples 100 people (seed 2028), forms the
undirected union of four-nearest-neighbor nominations, and greedily retains edges
in increasing length if they cross no retained edge. Shared endpoints are allowed.
Final degrees can differ from four; no bridges are added. The embedded code keeps
component labels and verifies no proper crossings. Network calculations use group membership only. The slides show the final graph
and group-B shares at one, two, and three hops in a single row. Status maps remain
in the geographic microdata example.

Network figures share igraph's Fruchterman--Reingold layout (seed 2029).
The original spatial coordinates define the candidate/pruning graph; its retained
edges define hop distances. Display-layout coordinates are carried only for
mapping. An embedded check verifies identical indices and local environments
when coordinates are omitted from the network calculation.
