# %%
import copy
import concurrent.futures
import json
import multiprocessing
import time
from pathlib import Path
from typing import Callable

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_scores.tools.capra_singh_2007_scores as cs
import local_seqtools.general_utils as tools
import numpy as np
import pandas as pd
from local_conservation_scores import PairKmerConservationMethods
from local_conservation_scores import ColumnwiseScoreMethods
from attrs import asdict, define, field, validators
from local_config import conservation_pipeline_parameters as conf
import traceback
from dataclasses import dataclass


PAIRKCONSMETHODS = PairKmerConservationMethods()
COLSCOREMETHODS = ColumnwiseScoreMethods()


# %%

N_CORES = round(multiprocessing.cpu_count())


@define
class Params:
    scoreconfig: conf.PairKmerConservationParams
    score_key: str = "pairk_aln"
    level: str = "Vertebrata"
    keys_for_table: list[str] = field(
        default=["function_name", "matrix_name", "lflank", "rflank"]
    )


@define
class PairwiseScoreResults:
    hit_sequence: str
    hit_scores: list[float]
    hit_z_scores: list[float]
    flank_hit_sequence: str
    flank_hit_scores: list[float]
    flank_hit_z_scores: list[float]
    background_scores: list[float]


def lvlo_2_pairwise_scores(
    lvlo: group_tools.ConserLevel,
    score_key: str,
    # mat2score_func: str = "matrix_json_2_pairwise_scores",
    scoreconfig: conf.PairKmerConservationParams,
):
    pairk_cons_function = PAIRKCONSMETHODS.__getitem__(
        scoreconfig.kmer_conservation_function_name
    )
    col_function = COLSCOREMETHODS.__getitem__(
        scoreconfig.columnwise_score_function_name
    )

    pairdict = lvlo.conservation_scores[score_key]
    flanked_hit_scores = pairk_cons_function(
        pairdict["kmer_aln_file"],
        pairdict["flanked_hit_start_position_in_idr"],
        columnwise_score_func=col_function,
        bg_cutoff=scoreconfig.bg_cutoff,
        bg_kmer_cutoff=scoreconfig.bg_kmer_cutoff,
    )
    hit_slice = slice(
        pairdict["original_hit_st_in_flanked_hit"],
        pairdict["original_hit_end_in_flanked_hit"] + 1,
    )
    scores = PairwiseScoreResults(
        hit_sequence=flanked_hit_scores.hit_sequence[hit_slice],
        hit_scores=flanked_hit_scores.hit_scores[hit_slice],
        hit_z_scores=flanked_hit_scores.hit_z_scores[hit_slice],
        flank_hit_sequence=flanked_hit_scores.hit_sequence,
        flank_hit_scores=flanked_hit_scores.hit_scores,
        flank_hit_z_scores=flanked_hit_scores.hit_z_scores,
        background_scores=flanked_hit_scores.background_scores,
    )
    return scores


def process_row(row, params: Params):
    json_file = row["json_file"]
    og = group_tools.ConserGene(json_file)
    row["score_key"] = params.score_key
    if hasattr(og, "critical_error"):
        row["errors"] = f"critical_error: {og.critical_error}"
        return row
    og.load_levels()
    if og.level_objects is None:
        row["errors"] = f"level {params.level} not in og.level_objects"
        return row
    if params.level not in og.level_objects:
        row["errors"] = f"level {params.level} not in og.level_objects"
        return row
    lvlo = og.level_objects[params.level]
    if params.score_key not in lvlo.conservation_scores:
        row["errors"] = f"score_key {params.score_key} not in lvlo.conservation_scores"
        return row
    score_dict = copy.deepcopy(lvlo.conservation_scores[f"{params.score_key}"])
    row["k"] = len(score_dict["flanked_hit"])
    for k in score_dict:
        if isinstance(score_dict[k], dict):
            for j in score_dict[k]:
                if j in params.keys_for_table:
                    row[f"{j}"] = score_dict[k][j]
            continue
        if k in params.keys_for_table:
            row[f"{k}"] = score_dict[k]
    try:
        scores = lvlo_2_pairwise_scores(
            lvlo=lvlo,
            score_key=params.score_key,
            scoreconfig=params.scoreconfig,
        )
    except ValueError as e:
        row["errors"] = str(e)
        return row
    row["flanked_hit_sequence"] = scores.flank_hit_sequence
    row["flanked_hit_scores"] = scores.flank_hit_scores
    row["flanked_hit_z_scores"] = scores.flank_hit_z_scores
    row["hit_scores"] = scores.hit_scores
    row["hit_z_scores"] = scores.hit_z_scores
    row["n_bg_scores"] = len(scores.background_scores)
    row["bg_STD"] = np.std(scores.background_scores)
    row["bg_mean"] = np.mean(scores.background_scores)
    for k, v in asdict(params.scoreconfig).items():
        row[k] = v
    return row


def process_chunk(chunk, params: Params):
    return chunk.apply(
        process_row,
        params=params,
        axis=1,
    )


def add_pairwise_scores_2_df(
    df: pd.DataFrame,
    scoreconfig: conf.PairKmerConservationParams,
    score_key: str,
    level: str,
    n_cores: int = N_CORES,
):
    params = Params(
        scoreconfig=scoreconfig,
        score_key=score_key,
        level=level,
    )
    chunks = np.array_split(df, n_cores)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        processed_chunks = list(
            executor.map(
                process_chunk,
                chunks,
                [params] * n_cores,
            )
        )
    df_result = pd.concat(processed_chunks)
    return df_result
