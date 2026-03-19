from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDWriter
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer


class CorrelationFilter(BaseEstimator, TransformerMixin):
    """
    Included so joblib can deserialize models that were trained with this
    custom transformer in the pipeline.
    """

    def __init__(self, threshold: float = 0.8):
        self.threshold = float(threshold)
        self.columns_: list[str] | None = None

    def fit(self, X, y=None):
        Xdf = self._to_df(X).copy()
        Xdf = Xdf.dropna(axis=1, how="all")
        Xdf = Xdf.apply(pd.to_numeric, errors="coerce")
        Xdf = Xdf.replace([np.inf, -np.inf], np.nan).fillna(0.0)

        nunique = Xdf.nunique(dropna=False)
        Xdf = Xdf.loc[:, nunique > 1]

        if Xdf.shape[1] > 1:
            corr = Xdf.corr(method="spearman", numeric_only=True).abs()
            upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
            to_drop = [c for c in upper.columns if (upper[c] > self.threshold).any()]
            Xdf = Xdf.drop(columns=to_drop, errors="ignore")

        self.columns_ = list(Xdf.columns)
        return self

    def transform(self, X):
        Xdf = self._to_df(X).copy()
        Xdf = Xdf.reindex(columns=self.columns_, fill_value=0.0)
        Xdf = Xdf.apply(pd.to_numeric, errors="coerce")
        Xdf = Xdf.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        return Xdf.values

    @staticmethod
    def _to_df(X):
        return X if isinstance(X, pd.DataFrame) else pd.DataFrame(X)


def get_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def get_default_model_path() -> Path:
    return get_repo_root() / "IRIS_model" / "IRIS.pkl"


def load_model(model_path: str | Path):
    import __main__

    model_path = Path(model_path).resolve()

    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")

    if not model_path.is_file():
        raise IsADirectoryError(f"Model path is not a file: {model_path}")

    # Backward compatibility for models pickled when CorrelationFilter
    # lived in a top-level script / __main__ namespace.
    setattr(__main__, "CorrelationFilter", CorrelationFilter)

    return joblib.load(model_path)


def read_scores_file(file_path: str | Path) -> pd.DataFrame:
    file_path = Path(file_path)

    if not file_path.is_file():
        raise FileNotFoundError(f"Descriptor file not found: {file_path}")

    with file_path.open("r") as handle:
        lines = handle.readlines()

    ligand_data_list: list[dict[str, float]] = []
    ligand_data: dict[str, float] = {}

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("Ligand ") and "Scores" in line:
            if ligand_data:
                ligand_data_list.append(ligand_data)
                ligand_data = {}
            continue

        parts = line.split(":", 1)
        if len(parts) != 2:
            continue

        key, value = parts[0].strip(), parts[1].strip()
        try:
            ligand_data[key] = float(value)
        except ValueError:
            pass

    if ligand_data:
        ligand_data_list.append(ligand_data)

    if not ligand_data_list:
        raise ValueError(f"No ligand descriptor blocks parsed from {file_path}")

    return pd.DataFrame(ligand_data_list)


def preprocess_data(df_features: pd.DataFrame) -> pd.DataFrame:
    imputer = SimpleImputer(strategy="most_frequent")
    X = imputer.fit_transform(df_features)
    return pd.DataFrame(X, columns=df_features.columns, index=df_features.index)


def predict_rmsd(model, X: pd.DataFrame) -> np.ndarray:
    return model.predict(X)


def read_sdf_and_attach_ranks(
    sdf_path: str | Path,
    df_complex: pd.DataFrame,
    sanitize: bool = False,
    score_field: str = "SCORE",
):
    sdf_path = Path(sdf_path)

    if not sdf_path.is_file():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    if score_field not in df_complex.columns:
        raise ValueError(
            f"'{score_field}' column not found in descriptors dataframe; "
            "cannot map SDF poses."
        )

    score_to_rank: dict[float, float] = {}
    for _, row in df_complex.iterrows():
        score = row.get(score_field)
        rank = row.get("Predicted_Rank")

        if pd.isna(score) or pd.isna(rank):
            continue

        score = float(score)
        if score not in score_to_rank:
            score_to_rank[score] = float(rank)

    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=sanitize)
    mols = []

    for mol in suppl:
        if mol is None:
            continue

        if not mol.HasProp(score_field):
            continue

        try:
            score = float(mol.GetProp(score_field))
        except Exception:
            continue

        if score not in score_to_rank:
            continue

        mol.SetProp("Predicted_Rank", str(score_to_rank[score]))
        mols.append(mol)

    mols.sort(key=lambda m: float(m.GetProp("Predicted_Rank")))
    return mols


def write_sorted_sdf(sorted_mols, output_sdf_path: str | Path) -> None:
    writer = SDWriter(str(output_sdf_path))
    try:
        for mol in sorted_mols:
            writer.write(mol)
    finally:
        writer.close()


def run_prediction(
    input_dir: str | Path,
    model_path: str | Path | None = None,
) -> Path:
    input_dir = Path(input_dir).resolve()

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    folder_name = input_dir.name
    resolved_model_path = (
        Path(model_path).resolve() if model_path is not None else get_default_model_path()
    )

    descriptors_path = input_dir / "descriptors_all.txt"
    sdf_file_path = input_dir / f"{folder_name}_docking_out_sorted.sdf"
    output_sdf_path = input_dir / f"{folder_name}_IRIS_sorted.sdf"

    if not descriptors_path.is_file():
        raise FileNotFoundError(f"Descriptor file not found: {descriptors_path}")

    if not sdf_file_path.is_file():
        raise FileNotFoundError(f"SDF file not found: {sdf_file_path}")

    model = load_model(resolved_model_path)
    df = read_scores_file(descriptors_path)

    feature_names = getattr(model, "feature_names_in_", None)
    if feature_names is None:
        feature_names = list(df.columns)

    required_features = list(feature_names)

    for feature in required_features:
        if feature not in df.columns:
            df[feature] = 0.0

    df_features = df[required_features].replace([np.inf, -np.inf], 0.0).fillna(0.0)
    X_input = preprocess_data(df_features)
    predicted_rmsd = predict_rmsd(model, X_input)

    df["Predicted_RMSD"] = predicted_rmsd
    df["Predicted_Rank"] = df["Predicted_RMSD"].rank(method="min", ascending=True)

    sorted_mols = read_sdf_and_attach_ranks(
        sdf_file_path,
        df,
        sanitize=False,
        score_field="SCORE",
    )

    if not sorted_mols:
        raise ValueError(
            "No poses could be mapped between SDF and descriptors. "
            "Common causes: SCORE mismatch or missing SCORE field."
        )

    write_sorted_sdf(sorted_mols, output_sdf_path)
    return output_sdf_path