import pandas as pd


table = pd.read_csv("../../../../benchmark/benchmark_v4/p1_table/benchmark_table.csv")

inds = list(table.sample(5).index)
inds.extend(list(table[table["name"] == "ABI3_HUMAN"].index))
inds.extend(list(table[table["name"] == "RAPH1_HUMAN"].index))
table.loc[inds].to_csv(
    "../../../../benchmark/benchmark_v4/p1_table/benchmark_table_testset.csv",
    index=False,
)
