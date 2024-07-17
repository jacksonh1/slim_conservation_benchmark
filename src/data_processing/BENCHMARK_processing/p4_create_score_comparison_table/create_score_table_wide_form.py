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


# %%
# levels = ["Tetrapoda", "Vertebrata", "Metazoa", "Eukaryota"]
levels = ["Tetrapoda", "Vertebrata", "Metazoa"]
df = pd.read_csv(
    "../../../../benchmark/benchmark_v4/p3_conservation/benchmark_table_ANNOTATED.csv"
)
output_folder = Path(
    "../../../../benchmark/benchmark_v4/p3_conservation/wide_form_tables_with_scores_v3/"
)
output_folder.mkdir(exist_ok=True, parents=True)


df = df[df["critical_error"].isnull()]
score_key_df = pd.DataFrame()

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
        output_file = output_folder / f"{score_index}.csv"
        df_temp = add_aln_scores_2_df(df, scoremethod.score_key, level, n_cores=50)
        df_temp.to_csv(output_file, index=False)
        score_key_df.loc[score_index, "score_index"] = score_index
        score_key_df.loc[score_index, "aln_type"] = "MSA - MAFFT"
        score_key_df.loc[score_index, "level"] = level
        score_key_df.loc[score_index, "score_key"] = scoremethod.score_key
        score_key_df.loc[score_index, "table_file"] = str(output_file.resolve())
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
            output_file = output_folder / f"{score_index}.csv"
            df_temp = add_pairwise_scores_2_df(
                df, mat_2_score_config, scoremethod.score_key, level, n_cores=50
            )
            df_temp.to_csv(output_file, index=False)
            score_key_df.loc[score_index, "table_file"] = str(output_file.resolve())
            score_key_df.loc[score_index, "score_index"] = score_index
            score_key_df.loc[score_index, "aln_type"] = "Pairwise"
            score_key_df.loc[score_index, "level"] = level
            score_key_df.loc[score_index, "score_key"] = scoremethod.score_key
            score_key_df.loc[score_index, "lflank"] = scoremethod.lflank
            score_key_df.loc[score_index, "rflank"] = scoremethod.rflank
            score_key_df.loc[score_index, "reciprocal_best_match"] = (
                mat_2_score_config.reciprocal_best_match
            )
            score_key_df.loc[score_index, "columnwise_score_function_name"] = (
                mat_2_score_config.columnwise_score_function_name
            )
            score_key_df.loc[score_index, "similarity_threshold"] = (
                mat_2_score_config.similarity_threshold
            )
            score_index += 1


# ==============================================================================
# // embedding pairwise scores
# ==============================================================================
# %%
thresholds = [1.0]
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
# levels = ["Tetrapoda", "Vertebrata", "Metazoa", "Eukaryota"]
levels = ["Tetrapoda", "Vertebrata", "Metazoa"]
for level in levels:
    for mat_2_score_config in mat_2_score_configs:
        for scoremethod in emb_pairwise_score_methods:
            if scoremethod.level is not None:
                if scoremethod.level != level:
                    continue
            output_file = output_folder / f"{score_index}.csv"
            df_temp = add_pairwise_scores_2_df(
                df, mat_2_score_config, scoremethod.score_key, level, n_cores=50
            )
            df_temp.to_csv(output_file, index=False)
            score_key_df.loc[score_index, "table_file"] = str(output_file.resolve())
            score_key_df.loc[score_index, "score_index"] = score_index
            score_key_df.loc[score_index, "aln_type"] = "Pairwise embedding"
            score_key_df.loc[score_index, "level"] = level
            score_key_df.loc[score_index, "score_key"] = scoremethod.score_key
            score_key_df.loc[score_index, "lflank"] = scoremethod.lflank
            score_key_df.loc[score_index, "rflank"] = scoremethod.rflank
            score_key_df.loc[score_index, "reciprocal_best_match"] = (
                mat_2_score_config.reciprocal_best_match
            )
            score_key_df.loc[score_index, "columnwise_score_function_name"] = (
                mat_2_score_config.columnwise_score_function_name
            )
            score_key_df.loc[score_index, "similarity_threshold"] = (
                mat_2_score_config.similarity_threshold
            )
            score_index += 1

# score_df.to_csv("../../../../benchmark/benchmark_v4/p3_conservation/score_table_v2.csv", index=False)
score_key_df.to_csv(output_folder / f"score_key.csv", index=False)

# # %%


# # %%
