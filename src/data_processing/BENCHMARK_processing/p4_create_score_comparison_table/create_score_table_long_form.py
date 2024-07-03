# %%
from pathlib import Path
from dataclasses import dataclass

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from typing import Callable
from local_conservation_scores.tools import capra_singh_2007_scores as cs
from attrs import asdict, define, field, validators
from local_config import conservation_pipeline_parameters as conf
# from local_conservation_scores import PairwiseMatrixKmerScoreMethods
from local_conservation_scores import ConservationScoreMethods, PairwiseMatrixMethods
# from local_conservation_scores import ColumnwiseScoreMethods
import pandas as pd
from add_aln_scores_2_table import add_aln_scores_2_df
# from add_pairwise_embedding_scores_2_table import add_pairwise_embedding_scores_2_df
from add_pairwise_scores_2_table import add_pairwise_scores_2_df
import yaml


# %load_ext autoreload
# %autoreload 2
# %%

# PAIRWISEMETHODS = PairwiseMatrixKmerScoreMethods()
# COLSCOREMETHODS = ColumnwiseScoreMethods()
PAIRWISEMATFUNCS = PairwiseMatrixMethods()
ALNSCORES = ConservationScoreMethods()


def load_config(config_file: str) -> conf.PipelineParameters:
    # if config_file is None:
    #     config = conf.PipelineParameters()
    # else:
    # with open(config_file, 'r') as f:
    #     config_dict = yaml.safe_load(f)
    # config = conf.PipelineParameters.from_dict(config_dict)
    with open(config_file, "r") as f:
        config_dict = yaml.safe_load(f)
    config = conf.PipelineParameters.from_dict(config_dict)
    return config


config = load_config("../p3_run_conservation_pipeline/params.yaml")
score_index = 0


aln_score_methods = []
pairwise_score_methods = []
for scoremethod in config.score_methods:
    if hasattr(ALNSCORES, scoremethod.score_function_name):
        aln_score_methods.append(scoremethod)
    elif hasattr(PAIRWISEMATFUNCS, scoremethod.score_function_name):
        pairwise_score_methods.append(scoremethod)
emb_pairwise_score_methods = [
    scoremethod for scoremethod in config.embedding_score_methods
]

# for i in pairwise_score_methods:
# print(i)

# %%
levels = ["Tetrapoda", "Vertebrata", "Metazoa", "Eukaryota"]
df = pd.read_csv(
    "../../../../benchmark/benchmark_v3/p3_conservation/benchmark_table_renamed_ANNOTATED.csv"
)
df = df[df["critical_error"].isnull()]
score_df = pd.DataFrame()

aln_id_vars = [
    "reference_index",
    "score_key",
    "errors",
    "score_function_name",
    "n_bg_scores",
    "bg_STD",
    "bg_mean",
]
for level in levels:
    for scoremethod in aln_score_methods:
        if scoremethod.level is not None:
            if scoremethod.level != level:
                continue
        df_temp = add_aln_scores_2_df(df, scoremethod.score_key, level, n_cores=50)
        df2concat = df_temp.melt(
            id_vars=aln_id_vars,
            value_vars=["hit_scores", "hit_z_scores"],
            var_name="score_type",
            value_name="score_list",
        )
        df2concat["aln_type"] = "MSA - MAFFT"
        df2concat["level"] = level
        df2concat['score_index'] = score_index
        score_df = pd.concat([score_df, df2concat]).reset_index(drop=True)
        score_index += 1

# ==============================================================================
# // pairwise scores
# ==============================================================================
# %%
pairwise_id_vars = [
    "reference_index",
    "score_key",
    "errors",
    "score_function_name",
    "k",
    "similarity_threshold",
    "reciprocal_best_match",
    "columnwise_score_function_name",
    "matrix_to_score_function_name",
    "matrix_name",
    "lflank",
    "rflank",
    "n_bg_scores",
    "bg_STD",
    "bg_mean",
]
kmerscoreobj = config.pairwise_matrix_to_score_params
# variables to modify
# columnwise_score_function_name
# reciprocal_best_match
# similarity_threshold
thresholds = [1.0]
reciprocal_best_match = [True, False]
columnwise_score_function_name = ["shannon_entropy", "property_entropy"]
# columnwise_score_function_name = ["shannon_entropy"]
mat_2_score_configs = []
for threshold in thresholds:
    for r in reciprocal_best_match:
        for c in columnwise_score_function_name:
            mat_2_score_configs.append(
                conf.PairMatrixToScoreConf(
                    matrix_to_score_function_name="pairwise_matrix_to_kmer_scores",
                    columnwise_score_function_name=c,
                    reciprocal_best_match=r,
                    similarity_threshold=threshold,
                )
            )


for i in mat_2_score_configs:
    print(i)
# %%
for level in levels:
    for mat_2_score_config in mat_2_score_configs:
        for scoremethod in pairwise_score_methods:
            if scoremethod.level is not None:
                if scoremethod.level != level:
                    continue
            df_temp = add_pairwise_scores_2_df(df, mat_2_score_config, scoremethod.score_key, level, n_cores=50)
            df2concat = df_temp.melt(
                id_vars=[i for i in pairwise_id_vars if i in df_temp.columns],
                value_vars=["hit_scores", "hit_z_scores"],
                var_name="score_type",
                value_name="score_list",
            )
            df2concat["aln_type"] = "Pairwise"
            df2concat["level"] = level
            df2concat['score_index'] = score_index
            score_df = pd.concat([score_df, df2concat]).reset_index(drop=True)
            score_index += 1


# ==============================================================================
# // embedding pairwise scores
# ==============================================================================
# %%
thresholds = [0.5, 1.0]
reciprocal_best_match = [False]
columnwise_score_function_name = ["shannon_entropy", "property_entropy"]
# columnwise_score_function_name = ["shannon_entropy"]
mat_2_score_configs = []
for threshold in thresholds:
    for r in reciprocal_best_match:
        for c in columnwise_score_function_name:
            mat_2_score_configs.append(
                conf.PairMatrixToScoreConf(
                    matrix_to_score_function_name="pairwise_matrix_to_kmer_scores",
                    columnwise_score_function_name=c,
                    reciprocal_best_match=r,
                    similarity_threshold=threshold,
                )
            )
levels = ["Tetrapoda", "Vertebrata", "Metazoa", "Eukaryota"]
for level in levels:
    for mat_2_score_config in mat_2_score_configs:
        for scoremethod in emb_pairwise_score_methods:
            if scoremethod.level is not None:
                if scoremethod.level != level:
                    continue

            df_temp = add_pairwise_scores_2_df(df, mat_2_score_config, scoremethod.score_key, level, n_cores=50)
            df2concat = df_temp.melt(
                id_vars=[i for i in pairwise_id_vars if i in df_temp.columns],
                value_vars=["hit_scores", "hit_z_scores"],
                var_name="score_type",
                value_name="score_list",
            )
            df2concat["aln_type"] = "Pairwise embedding"
            df2concat['ESM model'] = config.esm_params.model_name
            df2concat["level"] = level
            df2concat['score_index'] = score_index
            score_df = pd.concat([score_df, df2concat]).reset_index(drop=True)
            score_index += 1

score_df.to_csv("../../../../benchmark/benchmark_v3/p3_conservation/score_table_v3.csv", index=False)


# # %%


# # %%
