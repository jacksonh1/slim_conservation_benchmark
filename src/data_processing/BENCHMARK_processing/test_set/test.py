import json
import multiprocessing
import shutil
from functools import partial
from pathlib import Path

import yaml
from attrs import asdict
import copy
from typing import Callable

import local_config.conservation_pipeline_parameters as conf
from local_conservation_analysis_pipeline import (
    s1setup_folder,
    s2define_idrs,
    s3find_hit,
    s4add_lvlinfo,
    s5a_compute_aln_scores,
    s5b_compute_pairwise_matrices,
    s5c_add_hit_scores,
    s6multilevel_plots,
    s7output_aln_slice,
    s8calculate_annotations,
    s9add_annotations2table,
)

# import local_conservation_scores.tools.esm_model as esm_model
from local_seqtools import esm_tools

CONFIG_FILE = "./params.yaml"
N_CORES = multiprocessing.cpu_count()


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


def get_passing_jsons(search_dir):
    search_dir = Path(search_dir)
    json_files = search_dir.rglob("*.json")
    passing_jsons = []
    for json_file in json_files:
        if "clustered" in json_file.name:
            continue
        # print(json_file)
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if "reference_index" in json_dict:
                passing_jsons.append(json_file)
    return passing_jsons


def remove_failed_jsons(json_files):
    passing_jsons = []
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if "reference_index" in json_dict:
                passing_jsons.append(json_file)
    return passing_jsons


def main(config_file):
    config = load_config(config_file)
    config.pairk_conservation_params.kmer_conservation_function_name = (
        "pairwise_matrix_to_kmer_scores_rbm_penalty"
    )
    for scoremethod in config.score_methods:
        s5c_add_hit_scores.compute_hit_conservation_scores(
            json_file="../../../../benchmark/benchmark_v2/test/output/12-9606_0_000b76/12-9606_0_000b76.json",
            scoremethod=scoremethod,
            params=config.pairk_conservation_params,
        )


if __name__ == "__main__":
    main(CONFIG_FILE)
