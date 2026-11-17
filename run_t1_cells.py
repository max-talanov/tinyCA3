"""T1: measure the three cell types with DC current and save results/t1_cells.json and results/t1_fI.png.

Per type: resting potential, input resistance (-20 pA step), rheobase (500-ms DC, bisection to 0.1 pA) and the f-I
curve at 1, 1.5, 2, 3 and 5x rheobase (rate, first/last ISI, adaptation). Deterministic: no noise, fixed dt 0.025 ms.
Usage: .venv/bin/python run_t1_cells.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cells import (build_bc, build_pyr, build_rs, f_i_curve, h, input_resistance,  # noqa: E402
                   resting_potential, rheobase)

OUT = Path(__file__).resolve().parent / "results"
MULTS = (1.0, 1.5, 2.0, 3.0, 5.0)


def measure(name, builder):
    c = builder(name=name)
    s = c["soma"]
    rb = rheobase(s)
    res = dict(v_rest_mV=resting_potential(s), r_in_MOhm=input_resistance(s), rheobase_nA=rb,
               f_i=f_i_curve(s, [m * rb for m in MULTS]))
    del c
    for sec in list(h.allsec()):
        h.delete_section(sec=sec)
    return res


def main():
    OUT.mkdir(exist_ok=True)
    results = {"PYR": measure("pyr", build_pyr), "BC": measure("bc", build_bc), "RS": measure("rs", build_rs)}
    (OUT / "t1_cells.json").write_text(json.dumps(results, indent=1))
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for name, r in results.items():
        ax.plot([row["amp_nA"] / r["rheobase_nA"] for row in r["f_i"]], [row["rate_hz"] for row in r["f_i"]],
                "o-", label=f"{name} (rheobase {1000 * r['rheobase_nA']:.0f} pA)")
    ax.set_xlabel("DC current / rheobase")
    ax.set_ylabel("rate over 500 ms (Hz)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "t1_fI.png", dpi=150)
    for name, r in results.items():
        print(f"{name}: Vrest {r['v_rest_mV']:.1f} mV, Rin {r['r_in_MOhm']:.1f} MOhm, rheobase {1000 * r['rheobase_nA']:.1f} pA, "
              + ", ".join(f"{row['rate_hz']:.0f} Hz" for row in r["f_i"]))


if __name__ == "__main__":
    main()
