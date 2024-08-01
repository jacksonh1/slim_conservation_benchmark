import json
from pathlib import Path

json_search_path = (
    "../../../../benchmark/benchmark_v4/p3_conservation/conservation_analysis/"
)


def get_passing_jsons(search_dir, exclude_dir: str | Path | None = None):
    search_dir = Path(search_dir)
    json_files = search_dir.rglob("*.json")
    passing_jsons = []
    for json_file in json_files:
        if "clustered" in json_file.name:
            continue
        if exclude_dir is not None:
            if str(exclude_dir) in str(json_file):
                continue
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if "reference_index" in json_dict:
                passing_jsons.append(json_file)
    return passing_jsons


json_files = get_passing_jsons(json_search_path)


def load_json(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict


def fix_dict(d):
    ogs = d["orthogroups"]
    levels = list(ogs.keys())
    for level in levels:
        for score_key in ogs[level]["conservation_scores"].keys():
            score_dict = ogs[level]["conservation_scores"][score_key]
            if "score_function_name" in score_dict:
                score_dict["function_name"] = score_dict.pop("score_function_name")
            if "score_params" in score_dict:
                score_dict["function_params"] = score_dict.pop("score_params")
            if "matrix_file" in score_dict:
                score_dict["kmer_aln_file"] = score_dict.pop("matrix_file")
            if "mat2score_params" in score_dict:
                score_dict["pairk_conservation_params"] = score_dict.pop(
                    "mat2score_params"
                )
    return d


def save_json(json_file, d):
    with open(json_file, "w") as f:
        json.dump(d, f, indent=4)


for json_file in json_files:
    json_dict = load_json(json_file)
    fixed_dict = fix_dict(json_dict)
    save_json(json_file, fixed_dict)
    print(f"Fixed {json_file}")
