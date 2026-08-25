from pathlib import Path
import pandas as pd

abx_cols = [
    "abx_highcut_id33",
    "abx_basscut_id33",
    "abx_vocalreduce_id33",
    "abx_vocalreduce_id178",
    "abx_basscut_id178",
    "abx_highcut_id178",
    "abx_basscut_id180",
    "abx_highcut_id180",
    "abx_vocalreduce_id180",
]

if __name__ == "__main__":
    data_dir = Path(__file__).resolve().parent
    df = pd.read_csv(data_dir / "data_en.csv")

    participant_abx = df.groupby("participant_id")[abx_cols].first()
    correct_counts = participant_abx.sum(axis=1)
    excluded_ids = correct_counts[correct_counts < 8].index

    screened_df = df[~df["participant_id"].isin(excluded_ids)]
    screened_df.to_csv(data_dir / "data_en_screened.csv", index=False)
