# %%
from pathlib import Path

import local_conservation_scores.tools.capra_singh_2007_scores as cs
import pandas as pd
from add_aln_scores_2_table import add_aln_scores_2_table_file
from add_pairwise_embedding_scores_2_table import \
    add_pairwise_embedding_scores_2_df
from add_pairwise_scores_2_table import add_pairwise_scores_2_df

n_cores = 60

pairwise_score_keys = [
    # "fragpair_gapless_lf5_rf5_grantham",
    "fragpair_gapless_lf0_rf0_edssmat50",
    "fragpair_gapless_lf2_rf2_edssmat50",
    "fragpair_gapless_lf5_rf5_edssmat50",
    "fragpair_gapless_lf10_rf10_edssmat50",
    # "fragpair_gapless_lf10_rf10_grantham",
]

# embedding_score_keys = [
#     "fragpair_gapless_embedding_lf5_rf5"
# ]

# def gen_longform_from_df()


table_file = "../p3_conservation/benchmark_table_renamed_ANNOTATED.csv"
df = add_aln_scores_2_table_file(table_file, n_cores=n_cores)
for score_key in pairwise_score_keys:
    df = add_pairwise_scores_2_df(df, score_key=score_key, n_cores=n_cores, score_func=cs.property_entropy, reciprocal_best_match=True)
    df = add_pairwise_scores_2_df(df, score_key=score_key, n_cores=n_cores, score_func=cs.property_entropy, reciprocal_best_match=False)
# for score_key in embedding_score_keys:
#     df = add_pairwise_embedding_scores_2_df(df, score_key=score_key, n_cores=n_cores, score_func=cs.shannon_entropy)
df.to_csv(str(Path(table_file).name).replace(".csv", "_w_scores.csv"), index=False)

# %%








