# Project goals

Status: draft for user review. Owner: management.

## Mission

Make Reardon and colleagues' spatial, multigroup, ordinal, and income segregation
measures widely usable by sociologists through an efficient, portable, accessible,
open-source R package. Support varied populations and datasets, and enable carefully
specified extensions to neighborhoods defined by networks, including road networks.

## Success criteria

- **Scientific fidelity:** implemented measures map to original definitions, with
  documented assumptions, conventions, limitations, and independent reference cases.
  Package supplements and research proposals are clearly distinguished.
- **Useful workflows:** researchers can analyze individual records and aggregate
  counts, compare scales, obtain local environments, and interpret results using
  runnable examples. Binning, zoning, and surface assumptions remain visible.
- **Geographic and network support:** bandwidth and population-nearest-neighbor
  workflows have explicit distance/cost, eligibility, weighting, and tie semantics.
  Claims about network generalizations have supporting scientific evidence.
- **Practical performance:** representative workloads establish runtime, memory,
  and numerical agreement; optimization targets are set from those baselines.
- **Portability and maintenance:** supported environments are explicitly identified
  and checked; dependencies, installation, documentation, and licensing support reuse.

## Boundaries and priorities

Scientific correctness and clear interpretation take precedence over speed.
Favor a coherent public interface and reusable measurement components. Do not assume
that properties proven for geographic environments hold for arbitrary networks.

Proposed initial boundary: accept prepared networks and documented population
allocations; leave road topology repair and general network preparation to other tools.
Research manuscripts do not establish validated package functionality.

Original papers are in `project/references/`; research drafts are in `project/papers/`.
Workflow intent is documented in `README.md`
and `vignettes/`, especially the unfinished getting-started vignette. Those materials
guide the baseline audit but are not themselves proof of implementation correctness.
