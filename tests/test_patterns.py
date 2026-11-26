"""Pattern layout B (docs/architecture_update_plan_2026-11-26.md): disjoint basic patterns, same-kind combinations,
clean controls and the P1 synapse classes. Pure Python: no NEURON objects are built."""
import importlib
import os
import sys
from collections import Counter
from itertools import permutations
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="module")
def cB():
    old = os.environ.get("TINYCA3_LAYOUT")
    os.environ["TINYCA3_LAYOUT"] = "B"
    import circuit
    mod = importlib.reload(circuit)
    yield mod
    if old is None:
        os.environ.pop("TINYCA3_LAYOUT")
    else:
        os.environ["TINYCA3_LAYOUT"] = old
    importlib.reload(circuit)


def test_patterns_match_figure(cB):
    P = cB.PATTERNS
    assert cB.PATTERN_ORDER == ("X", "+", "=", "‖", "X∪+", "=∪‖")
    assert set(P["X"]) == {22, 24, 42, 44}
    assert set(P["+"]) == {23, 32, 33, 34, 43}
    assert set(P["="]) == {10 * r + c for r in (1, 5) for c in range(1, 6)}
    assert set(P["‖"]) == {10 * r + c for r in (2, 3, 4) for c in (1, 5)}


def test_basic_patterns_disjoint_and_cover_grid(cB):
    P = cB.PATTERNS
    assert cB.BASIC_DISJOINT
    cells = [c for k in cB.BASIC for c in P[k]]
    assert len(cells) == len(set(cells)) == 25
    assert set(cells) == set(cB.CELL_IDS)


def test_combinations_are_unions_of_two_same_kind(cB):
    P = cB.PATTERNS
    assert set(P["X∪+"]) == set(P["X"]) | set(P["+"])
    assert set(P["=∪‖"]) == set(P["="]) | set(P["‖"])
    for parts in cB.COMBOS.values():
        assert len(parts) == 2
        kinds = {k in cB.MOD_PATTERN.values() for k in parts}
        assert len(kinds) == 1, "a combination mixes a modulated and an unmodulated pattern"


def test_controls_never_get_a_modulator(cB):
    import stc
    stc = importlib.reload(stc)
    mod_cells = set().union(*stc.MOD_CELLS.values())
    assert stc.MOD_CELLS == {"DA": set(cB.PATTERNS["X"]), "5-HT": set(cB.PATTERNS["+"])}
    for k in ("=", "‖", "=∪‖"):
        assert not set(cB.PATTERNS[k]) & mod_cells
        assert k not in cB.MODULATOR_OF
    assert cB.MODULATOR_OF == {"X": "DA", "+": "5-HT", "X∪+": "DA+5-HT"}


def test_synapse_class_counts(cB):
    cls = Counter(cB.synapse_class(a, b) for a, b in permutations(cB.CELL_IDS, 2))
    assert cls == {"X": 12, "+": 20, "X→+": 20, "+→X": 20, "=": 90, "‖": 30, "=→‖": 60, "‖→=": 60, "none": 288}
    assert sum(cls.values()) == 600


def test_probes_cover_classes(cB):
    for name, (a, b) in cB.PROBES.items():
        assert a != b and a in cB.CELL_IDS and b in cB.CELL_IDS
    assert {cB.synapse_class(*ab) for ab in cB.PROBES.values()} == {"X", "+", "=", "‖", "X→+", "=→‖", "none"}
