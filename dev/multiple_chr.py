import io
import msprime
from pathlib import Path

output_path = "/Users/jliang2/Projects/test_TSK2ASR/data/simulations/normal"
Path(output_path).mkdir(parents=True, exist_ok=True)

mut_rate = 1.25e-8
mut_random_seed = 5678

ped_txt = """\
# id parent0 parent1 time is_sample
0   2   3   0.0 1
1   4   5   0.0 1
2   6   7   1.0 0
3   8   9   1.0 0
4   6   7   1.0 0
5   10  11  1.0 0
6   .   .   2.0 0
7   .   .   2.0 0
8   .   .   2.0 0
9   .   .   2.0 0
10  .   .   2.0 0
11  .   .   2.0 0
"""

Ls = [1000000, 2000000, 3000000]
#rs = [1e-8, 2e-8, 3e-8]
rate=[1e-7, 1e-8, 1e-7]


ts_chroms = []
pedigree = msprime.parse_pedigree(io.StringIO(ped_txt), sequence_length=1)

for i in range(len(Ls)):
    pedigree.sequence_length = Ls[i]
    rate_map = msprime.RateMap(
        position=[0, round(Ls[i] / 3), 2 * round(Ls[i] / 3), Ls[i]],
        rate=rate)

    ped_ts = msprime.sim_ancestry(
        initial_state=pedigree, model="fixed_pedigree",
        recombination_rate=rate_map, random_seed=i+1)

    ts_chroms.append(
        msprime.sim_ancestry(
            initial_state=ped_ts, population_size=1000,
            recombination_rate=rate_map, model="dtwf", random_seed=i+100))

for i, ts in enumerate(ts_chroms):
    print(f"chromosome {i} has length {ts.sequence_length} and {ts.num_trees} trees")
    tree_seq_mut = msprime.sim_mutations(ts, rate=mut_rate, random_seed=mut_random_seed)
    tree_seq_mut_tree_result = output_path + f"/msprime_chr{i}.trees"
    tree_seq_mut.dump(tree_seq_mut_tree_result)
