from pathlib import Path
import csv
import statistics as st

def print_phase1_summary(
    macsts_manifest="testData/out_phase1_useMacsMut_FALSE/phase1_macsTS_useMacsMut_FALSE_manifest.csv",
    msprime_manifest="testData/out_msprime_from_macs/msprime_manifest.csv",
    try_load_macsts_trees=True,
    comparison_table_path="testData/method_comparison_long.csv",
):
    macsts_manifest = Path(macsts_manifest)
    msprime_manifest = Path(msprime_manifest)

    def read_csv_rows(path):
        with path.open(newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))

    def mean_sd(xs):
        xs = [float(x) for x in xs if x is not None]
        if len(xs) == 0:
            return float("nan"), float("nan")
        if len(xs) == 1:
            return xs[0], 0.0
        return st.mean(xs), st.stdev(xs)

    macsts_rows = read_csv_rows(macsts_manifest)
    msprime_rows = read_csv_rows(msprime_manifest)

    by_scn_macsts = {}
    for r in macsts_rows:
        by_scn_macsts.setdefault(r["scenario"], []).append(r)

    by_scn_msprime = {}
    for r in msprime_rows:
        by_scn_msprime.setdefault(r["scenario"], []).append(r)

    # Optional: load MaCSTS trees for edge/node/tree/root stats
    macsts_tree_stats = {}
    if try_load_macsts_trees:
        try:
            import tskit
            for scn, rows in by_scn_macsts.items():
                e, n, t, h, m = [], [], [], [], []
                for r in rows:
                    p = Path(r["tree_path"])
                    if not p.exists():
                        p = macsts_manifest.parent.parent / p  # handle relative path in CSV
                    ts = tskit.load(str(p))
                    e.append(ts.num_edges)
                    n.append(ts.num_nodes)
                    t.append(ts.num_trees)
                    h.append(float(ts.max_root_time))
                    m.append(ts.num_mutations)
                macsts_tree_stats[scn] = {
                    "num_edges": mean_sd(e),
                    "num_nodes": mean_sd(n),
                    "num_trees": mean_sd(t),
                    "max_root_time": mean_sd(h),
                    "num_mutations": mean_sd(m)
                }
        except Exception as ex:
            print(f"[note] skipped MaCSTS tree loading: {ex}")

    scenarios = sorted(set(by_scn_macsts) | set(by_scn_msprime))

    for scn in scenarios:
        print(f"\n=== {scn} ===")

        # Mutation summary from MaCS vs MaCSTS manifest
        if scn in by_scn_macsts:
            d = by_scn_macsts[scn]
            macs = [float(r["macs_num_mutations"]) for r in d]
            macsts = [float(r["macsts_num_mutations"]) for r in d]
            delta = [a - b for a, b in zip(macs, macsts)]

            macs_mean, macs_sd = mean_sd(macs)
            macsts_mean, macsts_sd = mean_sd(macsts)
            delta_mean, delta_sd = mean_sd(delta)
            rel_diff = (delta_mean / macs_mean) if macs_mean != 0 else float("nan")

            print(
                f"scenario={scn}, n={len(d)}, "
                f"macs_mut_mean={macs_mean:.4f}, macs_mut_sd={macs_sd:.4f}, "
                f"macsts_mut_mean={macsts_mean:.4f}, macsts_mut_sd={macsts_sd:.4f}, "
                f"mut_diff_mean={delta_mean:.4f}, mut_diff_sd={delta_sd:.4f}, "
                f"rel_diff={rel_diff:.6f}"
            )

        # msprime stats from manifest (already has edges/nodes/tree/root)
        if scn in by_scn_msprime:
            d = by_scn_msprime[scn]
            edges = [float(r["num_edges"]) for r in d]
            nodes = [float(r["num_nodes"]) for r in d]
            trees = [float(r["num_trees"]) for r in d]
            roots = [float(r["max_root_time"]) for r in d]
            muts = [float(r["num_mutations"]) for r in d] if "num_mutations" in d[0] else []

            e_m, e_sd = mean_sd(edges)
            n_m, n_sd = mean_sd(nodes)
            t_m, t_sd = mean_sd(trees)
            h_m, h_sd = mean_sd(roots)
            print(
                f"msprime: edges_mean={e_m:.4f}, edges_sd={e_sd:.4f}, "
                f"nodes_mean={n_m:.4f}, nodes_sd={n_sd:.4f}, "
                f"trees_mean={t_m:.4f}, trees_sd={t_sd:.4f}, "
                f"root_time_mean={h_m:.4f}, root_time_sd={h_sd:.4f}"
            )
            if muts:
                m_m, m_sd = mean_sd(muts)
                print(f"msprime: num_mut_mean={m_m:.4f}, num_mut_sd={m_sd:.4f}")

        # Optional MaCSTS tree stats
        if scn in macsts_tree_stats:
            s = macsts_tree_stats[scn]
            print(
                f"macsts-ts: edges_mean={s['num_edges'][0]:.4f}, edges_sd={s['num_edges'][1]:.4f}, "
                f"nodes_mean={s['num_nodes'][0]:.4f}, nodes_sd={s['num_nodes'][1]:.4f}, "
                f"trees_mean={s['num_trees'][0]:.4f}, trees_sd={s['num_trees'][1]:.4f}, "
                f"root_time_mean={s['max_root_time'][0]:.4f}, root_time_sd={s['max_root_time'][1]:.4f},"
                f"num_mut_mean={s['num_mutations'][0]:.4f}, mut_sd={s['num_mutations'][1]:.4f}"
            )

    # Long-format per-replicate table for downstream comparison
    long_rows = []
    macsts_tree_by_key = {}
    if try_load_macsts_trees:
        try:
            for r in macsts_rows:
                p = Path(r["tree_path"])
                if not p.exists():
                    p = macsts_manifest.parent.parent / p
                ts = tskit.load(str(p))
                macsts_tree_by_key[(r["scenario"], int(r["rep"]), int(r["chr"]))] = {
                    "num_trees": ts.num_trees,
                    "num_edges": ts.num_edges,
                    "num_nodes": ts.num_nodes,
                    "max_root_time": float(ts.max_root_time),
                }
        except Exception as ex:
            print(f"[note] could not enrich macsTS tree stats in long table: {ex}")

    for r in macsts_rows:
        key = (r["scenario"], int(r["rep"]), int(r["chr"]))
        ts_stats = macsts_tree_by_key.get(key, {})
        long_rows.append({
            "Scenarios": r["scenario"],
            "Methods": "macs",
            "rep_index": int(r["rep"]),
            "num_mut": int(r["macs_num_mutations"]),
            "num_trees": "NA",
            "num_edges": "NA",
            "num_nodes": "NA",
            "max_root_time": "NA",
        })
        long_rows.append({
            "Scenarios": r["scenario"],
            "Methods": "macsTS",
            "rep_index": int(r["rep"]),
            "num_mut": int(r["macsts_num_mutations"]),
            "num_trees": ts_stats.get("num_trees", "NA"),
            "num_edges": ts_stats.get("num_edges", "NA"),
            "num_nodes": ts_stats.get("num_nodes", "NA"),
            "max_root_time": ts_stats.get("max_root_time", "NA"),
        })

    for r in msprime_rows:
        long_rows.append({
            "Scenarios": r["scenario"],
            "Methods": "msprime",
            "rep_index": int(r["rep"]),
            "num_mut": int(float(r["num_mutations"])) if "num_mutations" in r else "NA",
            "num_trees": int(float(r["num_trees"])) if "num_trees" in r else "NA",
            "num_edges": int(float(r["num_edges"])) if "num_edges" in r else "NA",
            "num_nodes": int(float(r["num_nodes"])) if "num_nodes" in r else "NA",
            "max_root_time": float(r["max_root_time"]) if "max_root_time" in r else "NA",
        })

    comparison_table_path = Path(comparison_table_path)
    comparison_table_path.parent.mkdir(parents=True, exist_ok=True)
    with comparison_table_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "Scenarios",
                "Methods",
                "rep_index",
                "num_mut",
                "num_trees",
                "num_edges",
                "num_nodes",
                "max_root_time",
            ],
        )
        w.writeheader()
        w.writerows(long_rows)
    print(f"\nSaved long comparison table: {comparison_table_path}")

print_phase1_summary()
