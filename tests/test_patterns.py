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


def _layout(name):
    old = os.environ.get("TINYCA3_LAYOUT")
    os.environ["TINYCA3_LAYOUT"] = name
    import circuit
    mod = importlib.reload(circuit)
    yield mod
    if old is None:
        os.environ.pop("TINYCA3_LAYOUT")
    else:
        os.environ["TINYCA3_LAYOUT"] = old
    importlib.reload(circuit)
    import stc
    importlib.reload(stc)


@pytest.fixture
def cB():
    yield from _layout("B")


@pytest.fixture
def cB7():
    yield from _layout("B7")


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


# Layout B7: B plus the mixed combination X∪‖ (a modulated and an unmodulated basic pattern).

def test_b7_order_and_mixed_pattern(cB7):
    P = cB7.PATTERNS
    assert cB7.PATTERN_ORDER == ("X", "+", "=", "‖", "X∪+", "=∪‖", "X∪‖")
    assert set(P["X∪‖"]) == set(P["X"]) | set(P["‖"])
    assert cB7.PMIX == "X∪‖"
    assert cB7.BASIC_DISJOINT
    assert len(cB7.PCOL) == 7 and cB7.RESULTS.name == "layout_B7"
    for k in ("X", "+", "=", "‖", "X∪+", "=∪‖"):
        assert set(P[k]) == set(cB7.LAYOUTS["B"]["patterns"][k])


def test_b7_mixed_combination_kinds(cB7):
    mixed = [k for k, parts in cB7.COMBOS.items()
             if len({p in cB7.MOD_PATTERN.values() for p in parts}) == 2]
    assert mixed == ["X∪‖"]


def test_b7_controls_never_get_a_modulator(cB7):
    import stc
    stc = importlib.reload(stc)
    mod_cells = set().union(*stc.MOD_CELLS.values())
    assert stc.MOD_CELLS == {"DA": set(cB7.PATTERNS["X"]), "5-HT": set(cB7.PATTERNS["+"])}
    for k in ("=", "‖", "=∪‖"):
        assert not set(cB7.PATTERNS[k]) & mod_cells
        assert k not in cB7.MODULATOR_OF
    assert cB7.MODULATOR_OF == {"X": "DA", "+": "5-HT", "X∪+": "DA+5-HT", "X∪‖": "DA"}


def test_b7_synapse_class_counts(cB7):
    cls = Counter(cB7.synapse_class(a, b) for a, b in permutations(cB7.CELL_IDS, 2))
    assert cls == {"X": 12, "+": 20, "X→+": 20, "+→X": 20, "=": 90, "‖": 30, "=→‖": 60, "‖→=": 60,
                   "X→‖": 24, "‖→X": 24, "none": 240}
    assert sum(cls.values()) == 600


def test_b7_probes_cover_classes(cB7):
    for name, (a, b) in cB7.PROBES.items():
        assert a != b and a in cB7.CELL_IDS and b in cB7.CELL_IDS
        assert cB7.synapse_class(a, b) == name or name == "X→="
    assert {cB7.synapse_class(*ab) for ab in cB7.PROBES.values()} == \
        {"X", "+", "=", "‖", "X→+", "=→‖", "X→‖", "‖→X", "none"}
