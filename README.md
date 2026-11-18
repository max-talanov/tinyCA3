# tinyCA3

The NEURON simulation of synaptic tagging and capture (STC) during LTP in a tiny CA3 circuit, with short-term
plasticity on every synapse and replay through the CA3 → EC → DG → CA3 loop.

- Circuit: [docs/ca3_basic_circuit_architecture.md](docs/ca3_basic_circuit_architecture.md)
- Implementation plan and stages: [docs/implementation_plan.md](docs/implementation_plan.md)
- STC and replay algorithm: [docs/stc_replay_algorithm_tinyCA3.md](docs/stc_replay_algorithm_tinyCA3.md)
- Literature: [docs/synaptic_tagging_literature_review.md](docs/synaptic_tagging_literature_review.md)
- Results per stage: [docs/RESULTS_2026-11-17.md](docs/RESULTS_2026-11-17.md)

## Setup

NEURON runs from a Python 3.13 venv with pip's `neuron`. Keep `PYTHONPATH` unset: the NEURON app bundle in
`PYTHONPATH` shadows pip's and does not load under Homebrew's Python 3.14.

```bash
python3.13 -m venv .venv
```

```bash
env -u PYTHONPATH .venv/bin/pip install -r requirements.txt
```

```bash
env -u PYTHONPATH .venv/bin/python -m pytest -q tests/
```

The mechanisms in `neuron/` are compiled on first use by `nrn_setup.load_mechanisms`.

## Layout

| Path | Content |
|---|---|
| `neuron/` | NMODL: the Tsodyks–Markram AMPA/NMDA synapse, its STC version `TMExp2SynNMDA_STC.mod`, and the Bezaire et al. 2016 channels (`BEZAIRE_CHANNELS.md`) |
| `reference/` | Numpy references (`tm_numpy.py`) and the original Bezaire templates used to check the ports |
| `params.py` | Every value with source, status and recording conditions |
| `cells.py` | T1: Hodgkin–Huxley-type cells: `build_pyr`, `build_bc`, `build_rs`; DC measurements |
| `run_t1_cells.py` | T1 measurements → `results/t1_cells.json`, `results/t1_fI.png` |
| `circuit.py` | Pilot circuit: all projections P1–P8, patterns, online blocks |
| `stc.py` | PRP pools, DA/5-HT schedules, offline integration of tags and captured weights |
| `run_pilot.py`, `make_figures.py` | 8-h pilot (`results/pilot.npz`; the 81-MB `pilot_events.pkl` is not committed) and Figures 2–6 |
| `replay_stats.py` | Null readout (hypergeometric chance level) and 5-seed statistics → `results/replay_stats.json`, Figure 7; needs `run_pilot.py --seed N --light [--no-replay]` for seeds 1–5 (`results/seeds/`) |
| `da_stats.py` | DA during replay (VTA mode, `run_pilot.py --replay-da 10`) against replay and no-replay, 5 seeds → `results/da_stats.json`, Figure 8 |
| `tests/` | One test file per stage |
