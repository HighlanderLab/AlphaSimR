"""
Generate msprime tree-sequence replicates from MaCS-style scenario strings.

This script maps a supported subset of MaCS args to msprime demography:
  - sample_size sequence_length
  -t, -r
  -I (with optional global migration parameter)
  -eN, -en, -eM, -em, -ej
  -s (ignored placeholder)

Outputs:
  1) .trees files for each scenario/replicate/chromosome
  2) manifest CSV with paths/seeds/basic TS counts
"""



import argparse
import csv
from pathlib import Path

import msprime

def parse_macs_args(args):
    tok = [x for x in args.strip().split() if x]
    if len(tok) < 2:
        raise ValueError("MaCS args must start with '<sample_size> <sequence_length>'.")

    macs_arg = {'sample_size': int(tok[0]),
                'sequence_length': int(tok[1]),
                'num_pops': 1,
                'pop_samples':[],
                'events':[]}

    i = 2
    n = len(tok)
    while i < n:
        flag = tok[i]
        if flag == "-t":
            macs_arg['theta'] = float(tok[i + 1])
            i += 2
        elif flag == "-r":
            macs_arg['rec'] = float(tok[i + 1])
            i += 2
        elif flag == "-I":
            k = int(tok[i + 1])
            pop_samples = [int(tok[i + 2 + j]) for j in range(k)]
            j = i + 2 + k
            mig = 0.0
            if j < n and not tok[j].startswith("-"):
                mig = float(tok[j])
                j += 1
            macs_arg['num_pops'] = k
            macs_arg['pop_samples'] = pop_samples
            macs_arg['global_migration'] = mig
            i = j
        elif flag == "-eN":
            macs_arg['events'].append(
                {"type": "eN",
                 "t": float(tok[i + 1]),
                 "x": float(tok[i + 2])})
            i += 3
        elif flag == "-en":
            macs_arg['events'].append(
                {
                    "type": "en",
                    "t": float(tok[i + 1]),
                    "pop": int(tok[i + 2]),
                    "x": float(tok[i + 3]),
                }
            )
            i += 4
        elif flag == "-eM":
            macs_arg['events'].append(
                {"type": "eM",
                 "t": float(tok[i + 1]),
                 "M": float(tok[i + 2])})
            i += 3
        elif flag == "-em":
            macs_arg['events'].append(
                {
                    "type": "em",
                    "t": float(tok[i + 1]),
                    "src": int(tok[i + 2]),
                    "dst": int(tok[i + 3]),
                    "Mij": float(tok[i + 4]),
                }
            )
            i += 5
        elif flag == "-ej":
            macs_arg['events'].append(
                {
                    "type": "ej",
                    "t": float(tok[i + 1]),
                    "src": int(tok[i + 2]),
                    "dst": int(tok[i + 3]),
                }
            )
            i += 4
        elif flag == "-s":
            # In AlphaSimR wrappers this is a seed placeholder.
            # Ignore optional value if present.
            if i + 1 < n and not tok[i + 1].startswith("-"):
                i += 2
            else:
                i += 1
        else:
            raise ValueError(f"Unsupported MaCS token: {flag}")

    if macs_arg['num_pops'] == 1 and len(macs_arg['pop_samples'])==0:
        macs_arg['pop_samples'] = [macs_arg['sample_size']]

    if sum(macs_arg['pop_samples']) != macs_arg['sample_size']:
        raise ValueError(
            f"-I population samples sum to {sum(macs_arg['pop_samples'])}, "
            f"but sample_size is {macs_arg['sample_size']}."
        )
    return macs_arg


def directed_pairs(pop_names):
    return [(src, dst) for src in pop_names for dst in pop_names if src != dst]


def scaled_time_to_generations(t_coal, nref):
    return t_coal * 4.0 * nref


def scaled_rec_to_per_bp(r_arg, nref):
    return r_arg / (4.0 * nref)


def scaled_mut_to_per_bp(theta, nref):
    return theta / (4.0 * nref)


def scaled_global_M_to_pairwise_m(M, k_pops, nref):
    if k_pops <= 1:
        return 0.0
    return (M / (k_pops - 1.0)) / (4.0 * nref)


def scaled_pair_Mij_to_m(Mij, nref):
    return Mij / (4.0 * nref)


def build_demography(macs_arg, nref):
    dem = msprime.Demography()
    pop_names = [f"p{i + 1}" for i in range(macs_arg['num_pops'])]
    for name in pop_names:
        dem.add_population(name=str(name), initial_size=nref, initially_active=True)

    if macs_arg['num_pops'] > 1:
        m0 = scaled_global_M_to_pairwise_m(macs_arg['global_migration'], macs_arg['num_pops'], nref)
        if m0!=0:
            for src, dst in directed_pairs(pop_names):
                dem.set_migration_rate(str(src), str(dst), m0)

    events = sorted(macs_arg['events'], key=lambda e: e["t"])
    for ev in events:
        t_gen = scaled_time_to_generations(ev["t"], nref)
        typ = ev["type"]
        if typ == "eN":
            size = ev["x"] * nref
            for name in pop_names:
                dem.add_population_parameters_change(time=t_gen, population=str(name), initial_size=size)
        elif typ == "en":
            pop_name = f"p{ev['pop']}"
            size = ev["x"] * nref
            dem.add_population_parameters_change(time=t_gen, population=str(pop_name), initial_size=size)
        elif typ == "eM":
            m = scaled_global_M_to_pairwise_m(ev["M"], macs_arg['num_pops'], nref)
            for src, dst in directed_pairs(pop_names):
                dem.add_migration_rate_change(time=t_gen, source=str(src), dest=str(dst), rate=m)
        elif typ == "em":
            src = f"p{ev['src']}"
            dst = f"p{ev['dst']}"
            m = scaled_pair_Mij_to_m(ev["Mij"], nref)
            dem.add_migration_rate_change(time=t_gen, source=str(src), dest=str(dst), rate=m)
        elif typ == "ej":
            src = f"p{ev['src']}"
            dst = f"p{ev['dst']}"
            dem.add_population_split(time=t_gen, derived=[str(src)], ancestral=str(dst))
        else:
            raise ValueError(f"Unhandled event type: {typ}")

    rec_rate = scaled_rec_to_per_bp(macs_arg['rec'], nref)
    mut_rate = scaled_mut_to_per_bp(macs_arg['theta'], nref)
    pop_hap = [int(macs_arg["pop_samples"][i]) for i in range(macs_arg["num_pops"])]
    if any(n % 2 != 0 for n in pop_hap):
        raise ValueError("Cannot map MaCS haploid sample counts to ploidy=2 (odd count present).")
    samples = {f"p{i + 1}": pop_hap[i] // 2 for i in range(macs_arg["num_pops"])}
    # ploidy = 1
    #samples = {f"p{i + 1}": int(macs_arg['pop_samples'][i]) for i in range(macs_arg['num_pops'])}
    print(dem, samples, rec_rate)
    return dem, samples, rec_rate, mut_rate


def default_scenarios():
    return [
        {
            "id": 1,
            "name": "single_const",
            "args": "8 100000 -t 1e-3 -r 1e-4 -s ",
        },
        {
            "id": 2,
            "name": "single_eN",
            "args": "8 100000 -t 1e-3 -r 1e-4 -eN 0.2 2.0 -eN 1.0 0.5 -s ",
        },
        {
            "id": 3,
            "name": "I2_migration",
            "args": "8 100000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -eM 0.5 5e-3 -ej 1.0 2 1 -s ",
        },
        {
            "id": 4,
            "name": "I2_en_join",
            "args": "8 100000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -en 0.2 2 0.5 -ej 1.0 2 1 -s ",
        },
    ]


def make_seed(base_seed, scenario_id, rep_id, chr_id):
    return int(base_seed + scenario_id * 100000 + rep_id * 1000 + (chr_id - 1))


def run(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scenarios = default_scenarios()
    rows = []

    for sc in scenarios:
        macs_arg = parse_macs_args(sc["args"])
        dem, samples, rec_rate, mut_rate = build_demography(macs_arg, args.nref)
        for rep in range(1, args.nrep + 1):
            for chr_id in range(1, args.nchr + 1):
                seed = make_seed(args.base_seed, int(sc["id"]), rep, chr_id)
                ts = msprime.sim_ancestry(
                    samples=samples,
                    ploidy=2,
                    demography=dem,
                    sequence_length=macs_arg['sequence_length'],
                    recombination_rate=rec_rate,
                    model=args.model,
                    random_seed=seed,
                )
                mts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed+100)
                tree_path = out_dir / f"{sc['name']}_rep{rep:02d}_chr{chr_id:02d}.trees"
                mts.dump(tree_path)

                rows.append(
                    {
                        "scenario_id": sc["id"],
                        "scenario": sc["name"],
                        "rep": rep,
                        "chr": chr_id,
                        "args": sc["args"],
                        "nref": args.nref,
                        "seed_chr": seed,
                        "model": args.model,
                        "sequence_length": macs_arg['sequence_length'],
                        "rec_rate_bp": rec_rate,
                        "mut_rate_bp": mut_rate,
                        "num_trees": ts.num_trees,
                        "num_nodes": ts.num_nodes,
                        "num_edges": ts.num_edges,
                        "num_mutations": mts.num_mutations,
                        "max_root_time": float(ts.max_root_time),
                        "tree_path": str(tree_path),
                    }
                )
                if args.verbose:
                    print(
                        f"[ok] {sc['name']} rep={rep:02d} chr={chr_id:02d} "
                        f"trees={ts.num_trees} nodes={ts.num_nodes} edges={ts.num_edges}"
                    )

    manifest_path = out_dir / "msprime_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved manifest: {manifest_path}")
    print(f"Saved trees under: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nrep", type=int, default=50, help="Replicates per scenario.")
    p.add_argument("--nchr", type=int, default=1, help="Independent chromosomes per replicate.")
    p.add_argument("--nref", type=float, default=10000.0, help="Reference Ne for MaCS->msprime scaling.")
    p.add_argument("--base-seed", type=int, default=700000, help="Base seed for deterministic seed schedule.")
    p.add_argument("--model", type=str, default="smc_prime", help="msprime ancestry model.")
    p.add_argument(
        "--out-dir",
        type=str,
        default="testData/out_msprime_from_macs",
        help="Output directory for .trees and manifest CSV.",
    )
    p.add_argument("--verbose", action="store_true", help="Print progress lines.")
    return p


if __name__ == "__main__":
    parser = build_parser()
    ns = parser.parse_args()
    run(ns)
