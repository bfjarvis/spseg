# Third-party notices

## Boost

Euclidean neighborhood searches use Boost.Geometry; network searches use
Boost.Graph breadth-first and Dijkstra traversal. Boost headers are supplied
by the R package BH through LinkingTo; no Boost source is vendored or modified.

- Project: https://github.com/boostorg/geometry
- Project: https://github.com/boostorg/graph
- R distribution: https://CRAN.R-project.org/package=BH
- License: Boost Software License 1.0 (BSL-1.0).
- License text: https://www.boost.org/LICENSE_1_0.txt

Upstream headers retain their copyright and license notices. Record the BH
version with `packageVersion("BH")` when reporting compiled-build provenance.
