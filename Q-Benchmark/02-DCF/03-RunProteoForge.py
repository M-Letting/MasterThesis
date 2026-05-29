import os
import sys
import time
import shutil
import zipfile
import warnings
import tempfile
import multiprocessing as mp
from urllib.request import urlretrieve

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, "../.."))


def ensure_local_proteoforge(base_dir):
	"""Ensure local ProteoForge package exists; download from GitHub if missing."""
	local_pkg_dir = os.path.join(base_dir, "ProteoForge")
	local_init = os.path.join(local_pkg_dir, "__init__.py")

	if os.path.exists(local_init):
		return

	print("Local ProteoForge not found. Downloading from LangeLab/ProteoForge_Analysis...")
	zip_url = "https://github.com/LangeLab/ProteoForge_Analysis/archive/refs/heads/main.zip"

	with tempfile.TemporaryDirectory(prefix="proteoforge_bootstrap_") as tmp_dir:
		zip_path = os.path.join(tmp_dir, "ProteoForge_Analysis-main.zip")
		urlretrieve(zip_url, zip_path)

		with zipfile.ZipFile(zip_path, "r") as zf:
			zf.extractall(tmp_dir)

		source_pkg_dir = os.path.join(tmp_dir, "ProteoForge_Analysis-main", "ProteoForge")
		if not os.path.isdir(source_pkg_dir):
			raise FileNotFoundError(
				"Could not find 'ProteoForge' folder in downloaded repository archive."
			)

		shutil.copytree(source_pkg_dir, local_pkg_dir, dirs_exist_ok=True)

	if not os.path.exists(local_init):
		raise RuntimeError("ProteoForge bootstrap failed: local package was not created.")


ensure_local_proteoforge(project_root)

os.chdir(project_root)
sys.path.insert(0, project_root)

from ProteoForge import normalize, impute  # noqa: E402
from ProteoForge import disluster  # noqa: E402
from ProteoForge import weight, model  # noqa: E402

input_dir = "./Q-Benchmark/02-DCF/00-Data/02-Prepared"
output_dir = "./Q-Benchmark/02-DCF/00-Data/03-Results/ProteoForge"
os.makedirs(output_dir, exist_ok=True)

proteomaker_datasets = [
	{"name": "ProteoMakerLowNA",   "file": "ProteoMakerLowNA_wide_processed.csv"},
	{"name": "ProteoMakerLowNoNA", "file": "ProteoMakerLowNoNA_wide_processed.csv"},
	{"name": "ProteoMakerMedNA",   "file": "ProteoMakerMedNA_wide_processed.csv"},
	{"name": "ProteoMakerMedNoNA", "file": "ProteoMakerMedNoNA_wide_processed.csv"},
	{"name": "ProteoMakerHighNA",  "file": "ProteoMakerHighNA_wide_processed.csv"},
	{"name": "ProteoMakerHighNoNA","file": "ProteoMakerHighNoNA_wide_processed.csv"},
]

complex_files = [
	{"name": "B-Wich",                                          "file": "B-Wich_WideNoLogNoNorm.csv",                                          "complex_id": "CPX-1099"},
	{"name": "CRD-mediated",                                    "file": "CRD-mediated_WideNoLogNoNorm.csv",                                    "complex_id": "CPX-1080"},
	{"name": "Multiaminoacyl-tRNA",                             "file": "Multiaminoacyl-tRNA_WideNoLogNoNorm.csv",                             "complex_id": "CPX-2469"},
	{"name": "26S-proteasome",                                  "file": "26S-proteasome_WideNoLogNoNorm.csv",                                  "complex_id": "CPX-5993"},
	{"name": "60S-cytosolic-large-ribosomal-subunit",           "file": "60S-cytosolic-large-ribosomal-subunit_WideNoLogNoNorm.csv",           "complex_id": "CPX-5183"},
	{"name": "BSSADCRC",                                        "file": "BSSADCRC_WideNoLogNoNorm.csv",                                        "complex_id": "CPX-1203"},
	{"name": "Dynactin-complex",                                "file": "Dynactin-complex_WideNoLogNoNorm.csv",                                "complex_id": "CPX-26352"},
	{"name": "Dynein-1-complex-variant-1",                      "file": "Dynein-1-complex-variant-1_WideNoLogNoNorm.csv",                      "complex_id": "CPX-5025"},
	{"name": "Eukaryotic-translation-initiation-factor-3-complex", "file": "Eukaryotic-translation-initiation-factor-3-complex_WideNoLogNoNorm.csv", "complex_id": "CPX-6036"},
	{"name": "Intraflagellar-transport-complex-B",              "file": "Intraflagellar-transport-complex-B_WideNoLogNoNorm.csv",              "complex_id": "CPX-5022"},
	{"name": "Laminin-213",                                     "file": "Laminin-213_WideNoLogNoNorm.csv",                                     "complex_id": "CPX-1781"},
	{"name": "Major-Spliceosomal-B",                            "file": "Major-Spliceosomal-B_WideNoLogNoNorm.csv",                            "complex_id": "CPX-26435"},
	{"name": "Nuclear-pore-complex",                            "file": "Nuclear-pore-complex_WideNoLogNoNorm.csv",                            "complex_id": "CPX-873"},
]

model_to_use = "rlm"
correction = {
	"strategy": "two-step",
	"methods": ("bonferroni", "fdr_bh")
}

imputation_scheme = [
	{
		"method": "fill_dense",
		"params": {
			"max_missing": 1,
			"strategy": "mean"
		}
	},
	{
		"method": "downshift",
		"params": {
			"missingness_threshold": 1.0,
			"shift_magnitude": 1.5,
			"low_percentile": 0.10
		}
	},
	{
		"method": "knn",
		"params": {
			"n_neighbors": 5
		}
	}
]

total_cores = mp.cpu_count()
n_jobs = max(1, total_cores // 4)
mp.set_start_method("spawn", force=True)

def format_time_diff(time_diff):
	if time_diff < 60:
		return f"{time_diff:.2f} seconds"
	if time_diff < 3600:
		return f"{time_diff / 60:.2f} minutes"
	return f"{time_diff / 3600:.2f} hours"

def _get_intensity_cols_and_conditions(df, dataset_type):
	if dataset_type == "proteomaker":
		intensity_cols = [c for c in df.columns if pd.Series(c).str.match(r"^C\d+_R\d+$").any()]
		condition_map = {}
		for col in intensity_cols:
			condition = pd.Series(col).str.replace(r"_R\d+$", "", regex=True).iloc[0]
			condition_map.setdefault(condition, []).append(col)
	else:
		intensity_cols = [c for c in df.columns if pd.Series(c).str.match(r"^B\d+_D\d+$").any()]
		condition_map = {}
		for col in intensity_cols:
			condition = pd.Series(col).str.replace(r"^B\d+_", "", regex=True).iloc[0]
			condition_map.setdefault(condition, []).append(col)

	return intensity_cols, condition_map

def build_test_data(df, dataset_type):
	if dataset_type == "proteomaker":
		required = {"Peptidoform", "Accession"}
		is_log2 = True
		protein_col = "Accession"
		peptide_col = "Peptidoform"
	else:
		required = {"Peptide", "complex_id"}
		is_log2 = False
		protein_col = "complex_id"
		peptide_col = "Peptide"

	if not required.issubset(df.columns):
		raise ValueError("Dataset missing required peptide/protein columns")

	df = df.copy()
	intensity_cols, condition_map = _get_intensity_cols_and_conditions(df, dataset_type)
	if not intensity_cols:
		raise ValueError("No intensity columns found for dataset")

	df["protein_id"] = df[protein_col]
	df["peptide_id"] = df[peptide_col]
	df["feature_id"] = df["protein_id"].astype(str) + "||" + df["peptide_id"].astype(str)

	wide_cleaned = (
		df.set_index("feature_id")[intensity_cols]
		.apply(pd.to_numeric, errors="coerce")
	)
	if not is_log2:
		wide_cleaned = wide_cleaned.where(wide_cleaned > 0)

	wide_imputed = impute.run_imputation_pipeline(
		data=wide_cleaned,
		cond_dict=condition_map,
		scheme=imputation_scheme,
		is_log2=is_log2,
		return_log2=False,
		verbose=False
	)

	if is_log2:
		wide_raw = np.power(2.0, wide_cleaned)
	else:
		wide_raw = wide_cleaned.copy()

	melt_imputed = wide_imputed.reset_index().melt(
		id_vars="feature_id",
		var_name="Sample",
		value_name="Intensity"
	)
	melt_cleaned = wide_raw.reset_index().melt(
		id_vars="feature_id",
		var_name="Sample",
		value_name="Intensity(Raw)"
	)

	combined = melt_imputed.merge(
		melt_cleaned,
		on=["feature_id", "Sample"],
		how="left"
	)

	sample_to_condition = {
		sample: cond
		for cond, samples in condition_map.items()
		for sample in samples
	}
	combined["Condition"] = combined["Sample"].map(sample_to_condition)
	combined["isReal"] = combined["Intensity(Raw)"].notna()

	comp_miss = pd.DataFrame(index=wide_cleaned.index)
	for cond, cols in condition_map.items():
		comp_miss[cond] = wide_cleaned[cols].isna().all(axis=1)

	comp_miss = comp_miss.reset_index().melt(
		id_vars="feature_id",
		var_name="Condition",
		value_name="isCompleteMiss"
	)
	combined = combined.merge(
		comp_miss,
		on=["feature_id", "Condition"],
		how="left"
	)

	combined = combined.rename(columns={
		"Intensity": "intensity",
		"Intensity(Raw)": "intensity_raw"
	})
	combined["log10Intensity"] = np.log10(combined["intensity"])
	combined["day"] = combined["Condition"]
	combined["filename"] = combined["Sample"]

	combined = combined.merge(
		df[["feature_id", "protein_id", "peptide_id"]].drop_duplicates(),
		on="feature_id",
		how="left"
	)
	combined = combined.dropna(subset=["intensity"])

	return combined

def run_proteoforge_pipeline(long_df, control_condition):
	cond_run_dict = long_df.groupby("day")["filename"].unique().to_dict()

	test_data = normalize.against_condition(
		long_df,
		cond_run_dict,
		run_col="filename",
		index_cols=["protein_id", "peptide_id"],
		norm_against=control_condition,
		intensity_col="intensity",
		is_log2=False,
		norm_intensity_col="ms1adj"
	)
	print(f"Normalized intensities against {control_condition} (shape {test_data.shape})")

	test_data["log10Intensity"] = np.log10(test_data["intensity"])

	weights_data = weight.generate_weights_data(
		test_data,
		sample_cols=["filename"],
		log_intensity_col="log10Intensity",
		adj_intensity_col="ms1adj",
		control_condition=control_condition,
		condition_col="day",
		protein_col="protein_id",
		peptide_col="peptide_id",
		is_real_col="isReal",
		is_comp_miss_col="isCompleteMiss",
		sparse_imputed_val=10**-5,
		dense_imputed_val=0.5,
		verbose=False,
	)
	print(f"Generated weights (shape {weights_data.shape})")

	test_data["Weight"] = (
		(weights_data["W_Impute"] * 0.90) +
		(weights_data["W_RevTechVar"] * 0.10)
	)

	test_data["peptide_Idx"] = (
		test_data.groupby(["protein_id", "peptide_id"]).cumcount() + 1
	)

	cur_model = model.LinearModel(
		data=test_data,
		protein_col="protein_id",
		peptide_col="peptide_id",
		cond_col="day",
		intensity_col="ms1adj",
		weight_col="Weight",
	)

	test_data = cur_model.run_analysis(
		model_type=model_to_use,
		correction_strategy=correction["strategy"],
		correction_methods=correction["methods"],
		n_jobs=n_jobs
	)
	print(f"Modeling complete (shape {test_data.shape})")

	clusters = disluster.distance_and_cluster(
		data=test_data,
		protein_col="protein_id",
		peptide_col="peptide_id",
		cond_col="day",
		quant_col="ms1adj",
		clustering_params={
			"min_clusters": 1,
			"distance_transform": "corr",
			"clustering_method": "hybrid_outlier_cut",
			"linkage_method": "ward",
			"distance_metric": "euclidean"
		},
		n_jobs=n_jobs,
		verbose=False
	)
	print(f"Clustering complete (shape {clusters.shape})")

	test_data = test_data.merge(
		clusters[["protein_id", "peptide_id", "cluster_label"]],
		on=["protein_id", "peptide_id"],
		how="left"
	).rename(columns={"cluster_label": "ClusterID"})

	return test_data

def main():
	summary_rows = []

	# --- ProteoMaker datasets: processed individually ---
	for ds in proteomaker_datasets:
		input_path = os.path.join(input_dir, ds["file"])
		if not os.path.exists(input_path):
			continue

		output_path = os.path.join(output_dir, f"ProteoForge_{ds['name']}_result.feather")
		if os.path.exists(output_path):
			print(f"Skipping {ds['name']} — result already exists.")
			continue

		print("=" * 80)
		print(f"Processing dataset: {ds['name']}")
		print("=" * 80)
		start_time = time.time()

		df = pd.read_csv(input_path)
		print(f"Loaded {ds['name']} with shape {df.shape}")
		long_df = build_test_data(df, "proteomaker")
		print(f"Built long data with shape {long_df.shape}")

		cond_run_dict = long_df.groupby("day")["filename"].unique().to_dict()
		control_condition = "C1" if "C1" in cond_run_dict else sorted(cond_run_dict.keys())[0]

		test_data = run_proteoforge_pipeline(long_df, control_condition)

		test_data.to_feather(output_path)
		print(f"Saved results to {output_path}")

		discordant_peptides = int((test_data["ClusterID"] > 1).sum())
		proteoforms = test_data.loc[test_data["ClusterID"] > 1, "ClusterID"].nunique()
		summary_rows.append({
			"dataset": ds["name"],
			"method": "ProteoForge",
			"discordant_peptides": discordant_peptides,
			"proteoforms": int(proteoforms)
		})
		print(f"Completed {ds['name']} in {format_time_diff(time.time() - start_time)}")

	# --- Complex datasets: all combined into one run ---
	available_complexes = [
		cf for cf in complex_files
		if os.path.exists(os.path.join(input_dir, cf["file"]))
	]

	complex_output_paths = {
		cf["name"]: os.path.join(output_dir, f"ProteoForge_{cf['name']}_result.feather")
		for cf in available_complexes
	}

	if all(os.path.exists(p) for p in complex_output_paths.values()):
		print("Skipping all complex datasets — results already exist.")
	else:
		print("=" * 80)
		print(f"Processing all {len(available_complexes)} complex datasets combined")
		print("=" * 80)
		start_time = time.time()

		frames = []
		for cf in available_complexes:
			df = pd.read_csv(os.path.join(input_dir, cf["file"]))
			print(f"  Loaded {cf['name']} with shape {df.shape}")
			frames.append(df)

		combined_df = pd.concat(frames, ignore_index=True)
		print(f"Combined complex data shape: {combined_df.shape}")
		print(f"Unique complexes (protein_id): {combined_df['complex_id'].nunique()}")

		long_df = build_test_data(combined_df, "complex")
		print(f"Built long data with shape {long_df.shape}")

		cond_run_dict = long_df.groupby("day")["filename"].unique().to_dict()
		control_condition = "D17" if "D17" in cond_run_dict else sorted(cond_run_dict.keys())[0]

		test_data = run_proteoforge_pipeline(long_df, control_condition)

		for cf in available_complexes:
			complex_result = test_data[test_data["protein_id"] == cf["complex_id"]]
			complex_result.to_feather(complex_output_paths[cf["name"]])
			print(f"Saved results for {cf['name']} to {complex_output_paths[cf['name']]}")

			discordant_peptides = int((complex_result["ClusterID"] > 1).sum())
			proteoforms = complex_result.loc[complex_result["ClusterID"] > 1, "ClusterID"].nunique()
			summary_rows.append({
				"dataset": cf["name"],
				"method": "ProteoForge",
				"discordant_peptides": discordant_peptides,
				"proteoforms": int(proteoforms)
			})

		print(f"Completed all complexes in {format_time_diff(time.time() - start_time)}")

	if summary_rows:
		summary_df = pd.DataFrame(summary_rows)
		summary_df.to_csv(os.path.join(output_dir, "ProteoForge_summary.csv"), index=False)


if __name__ == "__main__":
	main()
