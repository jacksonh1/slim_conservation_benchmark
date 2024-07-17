# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: odb_conservation
#     language: python
#     name: python3
# ---

# %% [markdown]
# # plan
#
# - Get hit scores from the jsons.
#   - need to get the sequences (from alignment file or matrix files) from there as well to build the logo.
# - apply the position weighting to the scores.
#
#
# steps:
# - construct table with score lists
#
#
# - import table
# - add flanking sequence
# - add score lists (by score_key/level)
# - use weights to get ave scores from score lists
# - get difference between aln and pairwise
# - apply a few filters to get likely interesting hits
# - plot logos with 5 flanking residues
#   - aln w/ gaps | aln w/o gaps | pairwise
# -
#
#
#
#
# %%
#
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_scores.tools.pairwise_tools as pairwise_tools
import local_conservation_scores.tools.score_plots as score_plots
import local_seqtools.general_utils as tools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from local_seqtools import pssms

plt.style.use("custom_standard")
# plt.style.use('custom_small')
import seaborn as sns

# pd.options.plotting.backend = "plotly"

# %load_ext autoreload
# %autoreload 2

# %%
def convert_jsonfile_to_relative(json_file):
    return json_file.replace("/home/jch/Documents/07-pt1_pt2_combined/", "../../")


position_weights = {
    "DOC_WW_Pin1_4": np.array([0, 0, 0, 1, 1, 0]),  # ...([ST])P.
    "LIG_AP2alpha_2": np.array([1, 1, 1]),  # DP[FW]
    "LIG_EH_1": np.array([0, 1, 1, 1, 0]),  # .NPF.
    "LIG_SH2_GRB2like": np.array([1, 1, 1, 0]),  # (Y)([EDST]|[MLIVAFYHQW])N.
    "LIG_SH3_CIN85_PxpxPR_1": np.array([1, 0, 1, 0, 1, 1]),  # P.[AP].PR
    "enah_LPPPP_FPPPP": np.array([2, 1, 0, 1, 1]),  # [FWYL]P.[AFILTVYWP]P
    "TRAF6": np.array([0, 0, 0, 1, 0, 1, 0, 0, 1]),  # ...P.E..[FYWDE]
}

table_file = (
    "../../benchmark/benchmark_v3/p3_conservation/benchmark_table_renamed_ANNOTATED.csv"
)
df = pd.read_csv(table_file)
df = df[
    df["ELM_motif_class"] != "LIG_14-3-3_CanoR_1"
]  # this motif has a variable length regex and so it's more difficult to apply any position weighting
df = df[
    [
        "reference_index",
        "ELM_motif_class",
        "Organism",
        "UniprotID",
        "regex",
        "hit_sequence",
        "gene_id",
        "hit start position",
        "hit end position",
        "verified interaction",
        "name",
        "json_file",
        "critical_error",
    ]
]
df = df[df["critical_error"].isna()]
df["json_file"] = df["json_file"].apply(convert_jsonfile_to_relative)
df["weight_array"] = df["ELM_motif_class"].map(position_weights)
# %% [markdown]
# ## add flanking sequence

# %%
import local_env_variables.env_variables as env

data_all_seqrecords_dict = env.load_data_all_odb_seqs()
df["odb_seq"] = df["gene_id"].apply(
    lambda x: (
        str(data_all_seqrecords_dict[x].seq)
        if x in data_all_seqrecords_dict.keys()
        else False
    )
)
df["flanked_hit"] = df.apply(
    lambda x: tools.pad_with_aas_or_gaps(
        x["odb_seq"], x["hit start position"], x["hit end position"] + 1, flank=5
    ),
    axis=1,
    result_type="expand",
)
df = df.drop("odb_seq", axis=1)

# %%
temp = df.loc[0]
og = group_tools.ConserGene(temp["json_file"])
og.load_levels(filepath_converter=convert_jsonfile_to_relative)
# lvlaln=og.get_aln_score_obj('Metazoa', 'aln_shannon_entropy')
# for aa, s in zip(lvlaln.hit_aln_sequence, lvlaln.hit_aln_scores):
#     print(aa, s)
for level, lvlo in og.level_objects.items():
    for scorekey in lvlo.conservation_scores:
        print(level, scorekey)

# %% [markdown]
# ## add score lists for a pairwise score and alignment score at a specific phylogenetic level
# - level - Metazoa
# - score_keys
#   - aln_property_entropy
#   - fragpair_gapless_lf5_rf5_edssmat50

# %%
def json_2_z_score_list(json_file, level, scorekey):
    og = group_tools.ConserGene(
        json_file, filepath_converter=convert_jsonfile_to_relative
    )
    if level not in og.levels_passing_filters:
        return
    lvlo = og.get_level_obj(level, filepath_converter=convert_jsonfile_to_relative)
    if scorekey not in lvlo.conservation_scores:
        return
    if "hit_z_scores" not in lvlo.conservation_scores[scorekey]:
        return
    return lvlo.conservation_scores[scorekey]["hit_z_scores"]

def add_scorelist_2_df(df, level, scorekey):
    colname = f"{level}_{scorekey}_z_scores"
    df[colname] = df["json_file"].apply(
        lambda x: json_2_z_score_list(x, level, scorekey)
    )
    return df

# %%
from attrs import asdict, define, field, validators


@define
class PairwiseScoreResults:
    flanked_hit: str
    flanked_hit_start_position_in_idr: int
    original_hit_st_in_flanked_hit: int
    original_hit_end_in_flanked_hit: int
    score_function_name: str
    score_params: dict
    lflank: int
    rflank: int
    matrix_file: str | Path
    flanked_hit_sequence: str
    flanked_hit_scores: list
    flanked_hit_z_scores: list
    hit_sequence: str
    hit_scores: list
    hit_z_scores: list
    mat2score_params: dict

    def __attrs_post_init__(self):
        self.matrix_file = convert_jsonfile_to_relative(self.matrix_file)


@dataclass
class AlnScoreResults:
    file: str
    score_function_name: str
    score_params: dict
    hit_scores: list
    hit_z_scores: list


def slice_aln_scores(lvlo: group_tools.LevelAlnScore, aln_start, aln_end):
    hit_slice = slice(aln_start, aln_end + 1)
    hit_scores = lvlo.scores[hit_slice]
    hit_z_scores = lvlo.z_scores[hit_slice]
    hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
    return hit_scores, hit_z_scores, hit_aln_seq


def json2logoplot_alnscore(
    jsonfile, score_key, with_gaps=False, axes=None, level="Vertebrata", flank=5
):
    og = group_tools.ConserGene(
        jsonfile, filepath_converter=convert_jsonfile_to_relative
    )
    lvlo = og.get_aln_score_obj(
        level, score_key, filepath_converter=convert_jsonfile_to_relative
    )
    # aln = jch_alignment.jch_alignment(lvlo.aln, og.query_gene_id)
    flst, flend, flhit = tools.pad_hit(
        og.query_idr_sequence,
        og.hit_st_in_idr,
        og.hit_end_in_idr,
        l_flank=flank,
        r_flank=flank,
    )
    query_idr, index = tools.reindex_alignment_str(
        lvlo.query_aln_sequence[lvlo.idr_aln_start : lvlo.idr_aln_end + 1]
    )
    flstaln, flendaln = index[flst], index[flend]
    flanked_hit_scores, flanked_hit_z_scores, flhit_aln_seq = slice_aln_scores(
        lvlo, flstaln + lvlo.idr_aln_start, flendaln + lvlo.idr_aln_start
    )
    idr_aln = lvlo.aln[:, lvlo.idr_aln_start : lvlo.idr_aln_end + 1]
    flhit_aln = idr_aln[:, flstaln : flendaln + 1]

    if not with_gaps:
        seqlist, query_slice, nongapinds = score_plots.strip_gaps_from_slice(
            flhit_aln, flhit_aln_seq
        )
        score_list = list(np.array(flanked_hit_z_scores)[nongapinds])
    else:
        seqlist = [str(i.seq) for i in list(flhit_aln)]
        query_slice = flhit_aln_seq
        score_list = flanked_hit_z_scores
    if axes is None:
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 4))
    score_plots.plot_score_bar_plot(
        ax=axes[0],
        score_list=score_list,
        query_seq=query_slice,
    )
    score_plots.plot_logo(
        ax=axes[1],
        str_list=seqlist,
        tick_label_str=query_slice,
    )
    counts = pssms.alignment_2_counts(seqlist, show_plot=False, heatmap=False)
    return counts


def json2logoplot(
    jsonfile, score_key, rbm: bool = False, axes=None, level="Vertebrata"
):
    og = group_tools.ConserGene(
        jsonfile, filepath_converter=convert_jsonfile_to_relative
    )
    lvlo = og.get_level_obj(level, filepath_converter=convert_jsonfile_to_relative)
    result = PairwiseScoreResults(**lvlo.conservation_scores[score_key])
    mat_dict = pairwise_tools.import_pairwise_matrices(result.matrix_file)
    subseqdf = mat_dict["subseq_dataframe"]
    subseqdf = subseqdf.fillna("-" * len(result.flanked_hit))
    if rbm:
        rbmdf = mat_dict["reciprocal_best_match_dataframe"]
        rbmdf = rbmdf.fillna(False)
        hitdf = pd.concat(
            [
                subseqdf.loc[result.flanked_hit_start_position_in_idr],
                rbmdf.loc[result.flanked_hit_start_position_in_idr],
            ],
            axis=1,
            keys=["subseq", "rbm"],
        )
        hitdf.loc["reference_kmer", "rbm"] = True
        seqlist = hitdf[hitdf["rbm"]]["subseq"].to_list()
    else:
        seqlist = subseqdf.loc[result.flanked_hit_start_position_in_idr].to_list()
    if axes is None:
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 4))
    if not rbm:
        score_plots.plot_score_bar_plot(
            ax=axes[0],
            score_list=result.flanked_hit_z_scores,
            query_seq=result.flanked_hit_sequence,
        )
    score_plots.plot_logo(
        ax=axes[1], str_list=seqlist, tick_label_str=result.flanked_hit_sequence
    )
    counts = pssms.alignment_2_counts(seqlist, show_plot=False, heatmap=False)
    return counts, result


def plots_from_df(df, level, pairkey = "fragpair_gapless_lf5_rf5_edssmat50", output_folder=None, flank=5, rbm=False):
    counter = 0
    for i, row in df.iterrows():
        jsonfile = row["json_file"]
        fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(20, 4))
        _ = json2logoplot_alnscore(
            jsonfile,
            "aln_property_entropy",
            axes=ax[:, 0],
            level=level,
            with_gaps=True,
            flank=flank,
        )
        # automatically uses whatever scoring parameters are in the json files (should be just normal right now)
        _ = json2logoplot_alnscore(
            jsonfile,
            "aln_property_entropy",
            axes=ax[:, 1],
            level=level,
            with_gaps=False,
            flank=flank,
        )
        _ = json2logoplot(
            jsonfile,
            pairkey,
            rbm=rbm,
            axes=ax[:, 2],
            level=level,
        )
        for axi in ax[0, :]:
            axi.set_ylim([-4, 4])
        # fig.suptitle(f'{row["name"]}-{row["UniprotID"]}-{row["reference_index"]}')
        plt.tight_layout()
        if output_folder is not None:
            # output_folder2 = Path(output_folder) / f'{row["reference_index"]}'
            # output_folder2.mkdir(parents=True, exist_ok=True)
            # fig.savefig(
            #     output_folder2
            #     / f'{level}-{row["name"]}-{row["UniprotID"]}-{row["pair-aln"]:.2f}-{row["gene_id"]}.png',
            #     bbox_inches="tight",
            #     dpi=300,
            # )
            fig.savefig(
                output_folder
                / f'{row["reference_index"]}-{row["name"]}-{row["UniprotID"]}-{row["pair"]:.2f}-{row["gene_id"]}.png',
                bbox_inches="tight",
                dpi=300,
            )
            plt.close(fig)
            counter += 1


# %%
def create_score_comparison_df(df, level, score_key1, score_key2):
    dftemp = df.copy()
    dftemp = add_scorelist_2_df(dftemp, level, score_key1)
    dftemp = add_scorelist_2_df(dftemp, level, score_key2)
    df2 = dftemp[
        ~dftemp[
            [
                f"{level}_{score_key1}_z_scores",
                f"{level}_{score_key2}_z_scores",
            ]
        ]
        .isna()
        .any(axis=1)
    ].copy()
    cols = [i for i in df2.columns if "_scores" in i]
    for col in cols:
        df2[col + "_weighted_mean"] = df2.apply(
            lambda x: np.average(x[col], weights=x["weight_array"]), axis=1
        )
    cols = [i for i in df2.columns if "weighted_mean" in i]
    rn = {
        f"{level}_{score_key1}_z_scores_weighted_mean": "aln",
        f"{level}_{score_key2}_z_scores_weighted_mean": "pair",
    }
    df2 = df2.rename(columns=rn)
    df2["pair-aln"] = df2["pair"] - df2["aln"]
    return df2

# %%

# %% [markdown]
# Metazoa plots

# %%
plt.rcParams.update({"font.size": 14})
LEVEL = "Metazoa"
SCORE_KEY1 = "aln_property_entropy"
SCORE_KEY2 = "fragpair_gapless_lf5_rf5_edssmat50"
df2 = create_score_comparison_df(df, LEVEL, SCORE_KEY1, SCORE_KEY2)
df2.to_csv("score_table.csv", index=False)

# %% [markdown]
# We lost a lot of motifs mostly because there weren't enough points in the z-score backgrounds. there are either too many gaps in the alignments or not enough k-mers in the pairwise method.

# %%
MOTIF = "enah_LPPPP_FPPPP"
output_folder = Path(f"./{MOTIF}/") / f"{LEVEL}"
output_folder.mkdir(exist_ok=True, parents=True)

enadf = df2[df2["ELM_motif_class"] == MOTIF].copy()
temp = enadf[(enadf["pair-aln"] > 2) | (enadf["UniprotID"].isin(['A0A5F9YLF9', 'Q8WXX7', 'Q9Y6W5']))].copy()
print(len(temp))
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, flank=5)


# %%
MOTIF = "TRAF6"
output_folder = Path(f"./{MOTIF}/") / f"{LEVEL}"
output_folder.mkdir(exist_ok=True, parents=True)
traf = df2[df2["ELM_motif_class"] == MOTIF].copy()
temp = traf[traf["pair-aln"] > 1].copy()
temp = temp[temp["pair"] > 1].copy()
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, pairkey="fragpair_gapless_lf0_rf0_edssmat50", flank=0)

# %%
temp = traf[traf["verified interaction"]].copy()
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, pairkey="fragpair_gapless_lf0_rf0_edssmat50", flank=0)

# %%
def make_folder_and_get_df(df, motif, level):
    output_folder = Path(f"./{motif}/") / f"{level}"
    output_folder.mkdir(exist_ok=True, parents=True)
    motdf = df[df["ELM_motif_class"] == motif].copy()
    return motdf, output_folder

# %% [markdown]
# Vertebrata plots

# %%
LEVEL = "Vertebrata"
MOTIF = "enah_LPPPP_FPPPP"
df2 = create_score_comparison_df(df, LEVEL, SCORE_KEY1, SCORE_KEY2)
df2.to_csv("score_table_vertebrata.csv", index=False)
# %%
enadf, output_folder = make_folder_and_get_df(df2, MOTIF, LEVEL)
temp = enadf[(enadf["pair-aln"] > 1) | (enadf["UniprotID"].isin(['A0A5F9YLF9', 'Q8WXX7', 'Q9Y6W5'])) | (enadf["reference_index"].isin([2368]))].copy()
temp = temp[(temp["pair"] > 0.5) | (temp["UniprotID"].isin(['A0A5F9YLF9', 'Q8WXX7', 'Q9Y6W5'])) | (temp["reference_index"].isin([2368]))].copy()
print(len(temp))
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, flank=5)
# %%

# %%

MOTIF = "TRAF6"
traf, output_folder = make_folder_and_get_df(df2, MOTIF, LEVEL)
temp = traf[traf["pair-aln"] > 1].copy()
temp = temp[temp["pair"] > 1].copy()
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, pairkey="fragpair_gapless_lf0_rf0_edssmat50", flank=0)

# %%
temp = traf[traf["verified interaction"]].copy()
temp = temp.sort_values("pair", ascending=False)
plots_from_df(temp, LEVEL, output_folder=output_folder, pairkey="fragpair_gapless_lf0_rf0_edssmat50", flank=0)

# %% [markdown]
# Q9Y6W5
# WASF2_HUMAN
# Part of the WAVE complex that regulates lamellipodia formation

# %%
# temp=enadf[enadf['reference_index'] == 2368] # ROBO2
temp=enadf[enadf['UniprotID'] == 'Q70E73'] # raph1
plots_from_df(temp, LEVEL, output_folder=output_folder, flank=5)

# %%

LEVEL = "Metazoa"
MOTIF = "enah_LPPPP_FPPPP"
SCORE_KEY1 = "aln_property_entropy"
SCORE_KEY2 = "fragpair_gapless_lf5_rf5_edssmat50"
df2 = create_score_comparison_df(df, LEVEL, SCORE_KEY1, SCORE_KEY2)
enadf, output_folder = make_folder_and_get_df(df2, MOTIF, LEVEL)
# %%
temp = enadf[enadf["reference_index"].isin([2368])].copy()
plots_from_df(temp, LEVEL, flank=5, rbm=True)

# %%

# %%

# %%

# %%

# %%

# %%

# %%
df2.loc[2887].to_json("./test.json", orient="records")

# %%
np.average(
    [
        -0.2100420216,
        -0.7622876615,
        0.2637430374,
        -0.3815422762,
        -0.9570456609,
        -1.3029297428,
        -1.2251521697,
        -1.4368219208,
        -1.327423687,
    ]
)

# %%
np.average(
    [
        -0.2100420216,
        -0.7622876615,
        0.2637430374,
        -0.3815422762,
        -0.9570456609,
        -1.3029297428,
        -1.2251521697,
        -1.4368219208,
        -1.327423687,
    ],
    weights=[0, 0, 0, 1, 0, 1, 0, 0, 1],
)

# %%
np.average([-0.3815422762, -1.3029297428, -1.327423687])
np.mean([-0.3815422762, -1.3029297428, -1.327423687])

# %%

# %%

# %%

# %%

# %%

# %%
