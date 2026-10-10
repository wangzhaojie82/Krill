# Krill Source Code Repository

This repository contains the source code for **Krill**, a data-stream algorithm for detecting diffuse heavy flows under tight memory budgets.

A diffuse heavy flow is a high-frequency flow whose volume is substantially supported by many distinct elements with limited repetition. This property is useful in applications such as content access profiling, where the goal is to identify objects or flows with broad, distributed support. Krill detects such flows at the flow level. It uses a Diffuse Support Score (DSS) that rewards first occurrences of flow elements and penalizes repetitions. Instead of maintaining costly per-element statistics, Krill uses monotone hash-rank registers, called PairWitness, to obtain one-sided evidence that a flow-element pair is new. Since only a subset of new elements triggers register updates, Krill applies inverse-probability weighting to compensate for unobserved new-element contributions. The resulting DSS guides candidate admission and retention in FlowResidence, while frequency statistics determine whether a candidate is heavy.

## Repository Structure

```text
Krill/
|-- README.md
|-- MainExperiment/
`-- AblationStudy/
```

### MainExperiment

`MainExperiment/` contains the standalone implementation used for the main experimental evaluation.

It includes:

- Krill implementation
- Baselines used in the comparison
- Dataset loaders for the datasets
- Main experiment switches for baseline comparison, parameter sensitivity, stream-order robustness, and error breakdown

See `MainExperiment/README.md` for build instructions, dataset selection, and reproduction steps.

### AblationStudy

`AblationStudy/` contains the standalone implementation used for ablation experiments.

It studies the effects of Krill's main components, including:

- PairWitness and alternative pair-identification mechanisms
- FlowResidence and ideal residence settings
- Challenge strategy variants
- Compensation and non-compensation variants

See `AblationStudy/README.md` for build instructions and reproduction steps.
