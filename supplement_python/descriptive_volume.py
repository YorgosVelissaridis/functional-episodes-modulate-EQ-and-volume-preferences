from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ------------------------------
# Paths
# ------------------------------
DATA_FILE = Path("../data/data_en_screened.csv")
OUTPUT_DIR = Path("results")
OUT_FIG = OUTPUT_DIR / "descriptive_volume_eq_1db.pdf"
OUT_CSV = OUTPUT_DIR / "descriptive_volume_eq_1db.csv"

# ------------------------------
# Plot Categories and Labels
# ------------------------------
EQ_ORDER_DEFAULT = ["Flat", "BassCut", "HighCut", "VocalReduce", "V-Shape"]


def read_csv_with_fallback(path: Path) -> pd.DataFrame:
    """Read CSV with common encoding fallbacks used on macOS/Windows exports."""
    encodings = ["utf-8-sig", "cp932", "shift_jis"]
    last_error: UnicodeDecodeError | None = None

    for encoding in encodings:
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as exc:
            last_error = exc

    if last_error is not None:
        raise last_error
    raise RuntimeError(f"Failed to read CSV: {path}")


def build_1db_bins(level_series: pd.Series) -> tuple[np.ndarray, list[str]]:
    """Create 1 dB bin edges and human-readable labels."""
    min_edge = int(np.floor(level_series.min()))
    max_edge = int(np.ceil(level_series.max()))

    if min_edge == max_edge:
        max_edge += 1

    edges = np.arange(min_edge, max_edge + 1, 1)
    labels = [f"{edges[i]}-{edges[i + 1]} dB" for i in range(len(edges) - 1)]
    return edges, labels


def plot_eq_by_1db_bin(proportions: pd.DataFrame, eq_order: list[str], out_path: Path) -> None:
    """Plot EQ selection proportions per 1 dB adjusted-level bin."""
    eq_colors = plt.cm.Set2.colors
    eq_color_map = {eq: eq_colors[i % len(eq_colors)] for i, eq in enumerate(eq_order)}
    eq_color_map["Flat"] = "#b0b0b0"

    n_rows = len(proportions.index)
    fig_height = max(4.5, n_rows * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(11, fig_height))

    y_positions = np.arange(n_rows)
    left = np.zeros(n_rows)

    for eq in eq_order:
        vals = proportions[eq].to_numpy()
        bars = ax.barh(
            y_positions,
            vals,
            left=left,
            height=0.78,
            color=eq_color_map[eq],
            label=eq,
        )

        for bar, p in zip(bars, vals):
            if p >= 0.08:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{p * 100:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=9,
                )

        left += vals

    ax.set_yticks(y_positions, labels=proportions.index.tolist())
    ax.set_xlabel("EQ selection proportion (%)", fontsize=13)
    ax.set_ylabel("Adjusted volume level bin", fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
    ax.tick_params(axis="x", labelsize=11)
    ax.tick_params(axis="y", labelsize=10)
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.invert_yaxis()

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        ncol=min(len(eq_order), 5),
        loc="upper right",
        bbox_to_anchor=(0.98, 0.99),
        fontsize=11,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    df = read_csv_with_fallback(DATA_FILE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = df.dropna(subset=["selected_level_db", "selected_eq"]).copy()

    edges, labels = build_1db_bins(df["selected_level_db"])
    df["level_bin"] = pd.cut(
        df["selected_level_db"],
        bins=edges,
        labels=labels,
        right=False,
        include_lowest=True,
    )

    eq_observed = set(df["selected_eq"])
    eq_order = [eq for eq in EQ_ORDER_DEFAULT if eq in eq_observed]
    eq_order += [eq for eq in sorted(eq_observed) if eq not in eq_order]

    counts = (
        df.groupby(["level_bin", "selected_eq"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=eq_order, fill_value=0)
        .reindex(labels, fill_value=0)
    )

    # Show louder bins first (e.g., 5-6 dB at top, then 4-5 dB ...).
    counts = counts.iloc[::-1]

    row_totals = counts.sum(axis=1).replace(0, np.nan)
    proportions = counts.div(row_totals, axis=0).fillna(0.0)

    out_df = proportions.reset_index().rename(columns={"level_bin": "selected_level_db_bin"})
    out_df.to_csv(OUT_CSV, index=False)

    plot_eq_by_1db_bin(
        proportions=proportions,
        eq_order=eq_order,
        out_path=OUT_FIG,
    )


if __name__ == "__main__":
    main()
