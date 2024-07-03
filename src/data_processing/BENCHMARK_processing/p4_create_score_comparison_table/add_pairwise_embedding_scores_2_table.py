# %%
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

# %%

TABLE_FILE = "./benchmark_table_renamed_ANNOTATED.csv"
OUTPUT_FILE = TABLE_FILE.replace(".csv", "_w_pairwise_scores.csv")
SCORE_KEY="fragment_pairwise_gapless"
LEVEL="Vertebrata"
N_CORES = round(multiprocessing.cpu_count())


def import_pairwise_matrices(filepath):
    with open(filepath, "r") as json_file:
        data = json.load(json_file)
        score_df = pd.DataFrame(
            data["score_dataframe"]["data"],
            columns=data["score_dataframe"]["columns"],
            index=data["score_dataframe"]["index"],
        )
        subseq_df = pd.DataFrame(
            data["subseq_dataframe"]["data"],
            columns=data["subseq_dataframe"]["columns"],
            index=data["subseq_dataframe"]["index"],
        )
        position_df = pd.DataFrame(
            data["position_dataframe"]["data"],
            columns=data["position_dataframe"]["columns"],
            index=data["position_dataframe"]["index"],
        )

        matrix_dict = {
            "score_df": score_df,
            "subseq_df": subseq_df,
            "position_df": position_df,
        }
    return matrix_dict


def list_of_strings_to_list_of_columns(seqlist: list[str]):
    seqlist_filtered = [s for s in seqlist if isinstance(s, str)]
    col_strings = []
    for c in range(len(seqlist_filtered[0])):
        col = "".join([seqlist_filtered[i][c] for i in range(len(seqlist_filtered))])
        col_strings.append(col)
    return col_strings


# def list_of_strings_to_list_of_columns(seqlist: list[str]):
#     seqlist_filtered = [s for s in seqlist if isinstance(s, str)]
#     seqarray = np.array([np.array(list(i)) for i in seqlist_filtered])
#     col_strings = []
#     for c in range(len(seqlist_filtered[0])):
#         col = "".join(seqarray[:,c])
#         col_strings.append(col)
#     return col_strings


# def list_of_strings_to_list_of_columns(seqlist):
#     seqlist_filtered = [s for s in seqlist if isinstance(s, str)]
#     seqarray = np.array([list(i) for i in seqlist_filtered])
#     col_strings = [''.join(seqarray[:, c]) for c in range(seqarray.shape[1])]
#     return col_strings


def score_pseudo_aln_columns(
    seqlist: list[str], score_func: Callable = cs.property_entropy
):
    col_strings = list_of_strings_to_list_of_columns(seqlist)
    scores = []
    for c in col_strings:
        scores.append(score_func(c))
    return scores


def subseq_df_2_background_scores(subseq_df: pd.DataFrame, score_func: Callable = cs.property_entropy):
    background_scores = []
    for position in subseq_df.index:
        background_scores.extend(
            score_pseudo_aln_columns(subseq_df.loc[position, :].to_list(), score_func=score_func)
        )
    return background_scores


def json_file_2_pairwise_scores(
    json_file: str | Path,
    score_key: str = "fragment_pairwise_gapless",
    level: str = "Vertebrata",
    score_func: Callable = cs.property_entropy,
):
    og = group_tools.ConserGene(json_file)
    lvlo = og.get_level_obj(level)
    if score_key not in lvlo.conservation_scores:
        raise ValueError(f"score_key {score_key} not in lvlo.conservation_scores")
    pairdict = lvlo.conservation_scores[score_key]
    matrix_dict = import_pairwise_matrices(pairdict["file"])
    background_scores = subseq_df_2_background_scores(matrix_dict["subseq_df"], score_func=score_func)
    # background_scores = matrix_dict["subseq_df"].apply(lambda row: score_pseudo_aln_columns(row.tolist()), axis=1).explode().tolist()
    
    flank_hit_st = pairdict["flanked_hit_start_position_in_idr"]
    assert pairdict["flanked_hit"] == matrix_dict["subseq_df"].loc[flank_hit_st, "reference_kmer"]
    flank_hit_scores = score_pseudo_aln_columns(
        matrix_dict["subseq_df"].loc[flank_hit_st, :].to_list(),
        score_func=score_func,
    )
    hit_slice = slice(
        pairdict["original_hit_st_in_flanked_hit"],
        pairdict["original_hit_end_in_flanked_hit"] + 1,
    )
    flank_hit_z_scores = tools.z_score_comparison(flank_hit_scores, background_scores)
    hit_z_scores = flank_hit_z_scores[hit_slice]
    output_dict = {
        'hit_sequence': og.hit_sequence,
        'hit_z_scores': hit_z_scores,
        'flank_hit_sequence': pairdict['flanked_hit'],
        'flank_hit_z_scores': flank_hit_z_scores,
        'background_scores': background_scores,
        
    }
    return output_dict


def process_row(row, score_key, level, score_func=cs.property_entropy):
    json_file = row["json_file"]
    reference_index = row["reference_index"]
    colname = f"{score_key}-mean_zscore-{level}"
    try:
        score_dict = json_file_2_pairwise_scores(json_file, score_key=score_key, level=level, score_func=score_func)
    except ValueError as e:
        print(f"ValueError: {e}")
        row[f'{colname}-errors'] = str(e)
        return row
    except IndexError as e:
        print(f"IndexError: {e}")
        row[f'{colname}-errors'] = f"{str(e)}: probably no ortholog idrs (probably mostly gaps in the alignment)"
        return row
    except ZeroDivisionError as e:
        print(f"ZeroDivisionError: {e}")
        row[f'{colname}-errors'] = f"{str(e)}: probably no ortholog idrs (probably mostly gaps in the alignment)"
        return row
    row[colname] = np.mean(score_dict["hit_z_scores"])
    return row


def process_chunk(chunk, score_key, level, score_func=cs.property_entropy):
    return chunk.apply(process_row, score_key=score_key, level=level, score_func=score_func, axis=1)


def add_pairwise_embedding_scores_2_df(
    df: pd.DataFrame,
    score_key: str = "fragment_pairwise_gapless",
    level: str = "Vertebrata",
    n_cores: int = N_CORES,
    score_func: Callable = cs.property_entropy,
    # multiprocess=False,
):
    # if not multiprocess:
    #     for i, row in df.iterrows():
    #         print(f"Processing {row['reference_index']}...")
    #         df.loc[i] = process_row(row, score_key=score_key, level=level)
    #         return
    chunks = np.array_split(df, n_cores)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        processed_chunks = list(executor.map(
            process_chunk, chunks, [score_key] * n_cores, [level] * n_cores, [score_func] * n_cores
        ))
    df_result = pd.concat(processed_chunks)
    return df_result


def add_pairwise_embedding_scores_2_table_file(
    table_file: str | Path,
    score_key: str = "fragment_pairwise_gapless",
    level: str = "Vertebrata",
    output_file: str | Path | None = None,
    score_func: Callable = cs.property_entropy,
):
    df = pd.read_csv(table_file)
    df_filtered = df[df["critical_error"].isna()].copy()
    df_filtered = add_pairwise_embedding_scores_2_df(
        df_filtered, score_key=score_key, level=level, score_func=score_func
    )
    if output_file is None:
        return df_filtered
    df_filtered.to_csv(output_file, index=False)


if __name__ == "__main__":
    add_pairwise_embedding_scores_2_table_file(
        table_file=TABLE_FILE,
        score_key=SCORE_KEY,
        level=LEVEL,
        output_file=OUTPUT_FILE,
    )



# %%


# %%

# %%
# def add_pairwise_scores_2_df(
#     df: pd.DataFrame,
#     score_key: str = "fragment_pairwise_gapless",
#     level: str = "Vertebrata",
# ):
#     pairscore_map = {}
#     error_map = {}
#     for i, row in df.iterrows():
#         json_file = row["json_file"]
#         reference_index = row["reference_index"]
#         # if reference_index != 1267:
#         #     continue
#         # a = time.time()
#         print(f"Processing {reference_index}...")
#         try: 
#             pairscore_dict = json_file_2_pairwise_scores(json_file, score_key=score_key, level=level)
#             pairscore_map[reference_index] = np.mean(pairscore_dict['hit_z_scores'])
#         except ValueError as e:
#             print(f"ValueError: {e}")
#             error_map[reference_index] = str(e)
#             continue
#         except IndexError as e:
#             print(f"IndexError: {e}")
#             print(json_file)
#             error_map[reference_index] = f"{str(e)}: probably no ortholog idrs (probably mostly gaps in the alignment)"
#             continue
#         except ZeroDivisionError as e:
#             print(f"ZeroDivisionError: {e}")
#             error_map[reference_index] = f"{str(e)}: probably no ortholog idrs (probably mostly gaps in the alignment)"
#             continue
#         # b = time.time()
#         # print(f"Processed {reference_index} in {b-a} seconds")
#     colname = f"{score_key}-mean_zscore-{level}"
#     df[colname] = df["reference_index"].map(pairscore_map)
#     df[f'{colname}-errors'] = df["reference_index"].map(error_map)
#     return df
