from pathlib import Path
import matplotlib.pyplot as plt
import demesdraw
import msprime

from msprime_from_macs_scenarios import parse_macs_args, build_demography

NREF = 10_000  # choose the same Nref you use for MaCS<->msprime scaling
OUT = Path("testData/out_phase1_useMacsMut_FALSE")
OUT.mkdir(parents=True, exist_ok=True)

scenarios = [
    ("single_const", "8 100000 -t 1e-3 -r 1e-4 -s "),
    ("single_eN", "8 100000 -t 1e-3 -r 1e-4 -eN 0.2 2.0 -eN 1.0 0.5 -s "),
    ("I2_migration", "8 100000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -eM 0.5 5e-3 -ej 1.0 2 1 -s "),
    ("I2_en_join", "8 100000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -en 0.2 2 0.5 -ej 1.0 2 1 -s "),
]

for name, args in scenarios:
    macs_arg = parse_macs_args(args)
    dem, samples, rec_rate, mut_rate = build_demography(macs_arg, NREF)

    # Save textual check of event interpretation
    (OUT / f"{name}.debug.txt").write_text(str(dem.debug()), encoding="utf-8")

    # Demography plot (not tree-sequence plot)
    graph = msprime.Demography.to_demes(dem)
    fig, ax = plt.subplots(figsize=(8, 5))
    demesdraw.tubes(graph, ax=ax, seed=1, log_time=True)
    ax.set_title(f"{name}\n{args}")
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.demography.png", dpi=200)
    plt.close(fig)

print(f"Saved plots to: {OUT}")
