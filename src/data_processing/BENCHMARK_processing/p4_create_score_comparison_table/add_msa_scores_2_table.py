# %%
from pathlib import Path

import slim_conservation_scoring.pipeline.group_conservation_objects as group_tools
import pandas as pd
import slim_conservation_scoring.seqtools.general_utils as tools
import numpy as np
import re
import concurrent.futures
import multiprocessing
from attrs import asdict, define, field, validators

# %%


@define
class AlnScoreResults:
    hit_scores: list[float]
    hit_z_scores: list[float]
    background_scores: list[float] | None = None


N_CORES = round(multiprocessing.cpu_count())

# def get_hit_aln_scores(lvlo: group_tools.LevelAlnScore):
#     """
#     returns a list of the scores and a list of the z scores for the hit (non-gap) positions in query sequence
#     """
#     hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
#     hit_scores = lvlo.scores[hit_slice]
#     hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
#     nongap_inds = tools.get_non_gap_indexes(hit_aln_seq)
#     return list(np.array(hit_scores)[nongap_inds])


# def get_hit_aln_z_scores(lvlo: group_tools.LevelAlnScore):
#     hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
#     hit_z_scores = lvlo.z_scores[hit_slice]
#     hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
#     nongap_inds = tools.get_non_gap_indexes(hit_aln_seq)
#     return list(np.array(hit_z_scores)[nongap_inds])


# def json_to_aln_scores(
#     json_file: str | Path,
#     score_key: str = "aln_property_entropy",
#     level: str = "Vertebrata",
# ):
#     og = group_tools.ConserGene(json_file)
#     # print(og.info_dict)
#     output_dict = {}
#     if score_key not in og.info_dict["orthogroups"][level]["conservation_scores"]:
#         output_dict["error"] = f"error - score_key ({score_key}) not in json file"
#         return output_dict
#     score_o = og.get_aln_score_obj(
#         level,
#         score_key,
#     )
#     output_dict = {
#         "hit_scores": score_o.hit_scores,
#         "hit_sequence": og.hit_sequence,
#         "function_name": score_o.function_name,
#     }
#     for k,v in score_o.score_params.items():
#         output_dict[k] = v
#     if score_o.z_score_failure is not None:
#         output_dict['error'] = re.sub(r'not enough background scores to calculate z-score.*', 'not enough background scores to calculate z-score', str(score_o.z_score_failure))
#     else:
#         output_dict['hit_z_scores'] = score_o.hit_z_scores
#     return output_dict


def process_row(row, score_key, level):
    json_file = row["json_file"]
    reference_index = row["reference_index"]
    og = group_tools.ConserGene(json_file)
    row["score_key"] = score_key
    if score_key not in og.info_dict["orthogroups"][level]["conservation_scores"]:
        row["errors"] = f"error - score_key ({score_key}) not in json file"
        return row
    score_o = og.get_aln_score_obj(
        level,
        score_key,
    )
    row["hit_scores"] = score_o.hit_scores
    row["function_name"] = score_o.function_name
    for k, v in score_o.function_params.items():
        row[k] = v
    if score_o.z_score_failure is not None:
        row["errors"] = re.sub(
            r"not enough background scores to calculate z-score.*",
            "not enough background scores to calculate z-score",
            str(score_o.z_score_failure),
        )
    else:
        row["hit_z_scores"] = [float(i) for i in score_o.hit_z_scores]
        row["n_bg_scores"] = len(score_o.bg_scores)
        row["bg_STD"] = np.std(score_o.bg_scores)
        row["bg_mean"] = np.mean(score_o.bg_scores)
    return row


def process_chunk(chunk, score_key, level):
    return chunk.apply(process_row, score_key=score_key, level=level, axis=1)


def add_aln_scores_2_df(
    df: pd.DataFrame,
    score_key: str = "aln_property_entropy",
    level: str = "Vertebrata",
    n_cores: int = N_CORES,
):
    chunks = np.array_split(df, n_cores)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        processed_chunks = list(
            executor.map(
                process_chunk, chunks, [score_key] * n_cores, [level] * n_cores
            )
        )
    df_result = pd.concat(processed_chunks)
    return df_result


def add_aln_scores_2_table_file(
    table_file: str | Path,
    score_key: str = "aln_property_entropy",
    level: str = "Vertebrata",
    n_cores: int = N_CORES,
    output_file: str | Path | None = None,
):
    df = pd.read_csv(table_file)
    df_filtered = df[df["critical_error"].isna()].copy()
    df_result = add_aln_scores_2_df(
        df_filtered, score_key=score_key, level=level, n_cores=n_cores
    )
    # df_filtered = add_aln_scores_2_df(
    #     df_filtered, score_key=score_key, level=level
    # )
    if output_file is None:
        return df_result
    df_result.to_csv(output_file, index=False)


# %%
