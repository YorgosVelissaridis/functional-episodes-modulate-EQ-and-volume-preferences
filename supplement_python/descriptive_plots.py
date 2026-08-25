from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# ------------------------------
# Paths
# ------------------------------
DATA_FILE = Path("../data/data_en_screened.csv")
OUTPUT_DIR = Path("results")

# ------------------------------
# Plot Categories and Labels
# ------------------------------
EQ_ORDER_DEFAULT = ["Flat", "BassCut", "HighCut", "VocalReduce", "V-Shape"]
EPISODE_ORDER_DEFAULT = ["Null", "Enjoyment", "Distraction", "Relaxation", "Focus-Motivation"]
NOISE_ORDER_DEFAULT = ["EV", "DIESEL"]
NOISE_LABEL_MAP = {
    "EV": "EV Engine Noise",
    "DIESEL": "Diesel Engine Noise",
}


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


def draw_eq_stacked(
    ax: plt.Axes,
    counts: pd.DataFrame,
    row_order: list[str],
    eq_order: list[str],
    eq_color_map: dict[str, tuple[float, float, float, float]],
    y_labels: list[str],
    show_xlabel: bool,
    invert_y: bool,
) -> None:
    """Draw one horizontal stacked-bar panel of EQ selection proportions."""
    n_levels = len(row_order)
    bar_height = 0.8
    totals = counts.sum(axis=1)
    y_positions = list(range(n_levels))
    left = pd.Series(0.0, index=row_order)

    for eq in eq_order:
        proportions = counts.loc[row_order, eq] / totals.replace(0, 1)
        bars = ax.barh(
            y_positions,
            proportions.values,
            left=left.values,
            height=bar_height,
            color=eq_color_map[eq],
            label=eq,
        )

        for bar, p in zip(bars, proportions.values):
            if p > 0:
                x_center = bar.get_x() + bar.get_width() / 2
                y_center = bar.get_y() + bar.get_height() / 2
                ax.text(x_center, y_center, f"{p * 100:.1f}%", ha="center", va="center", fontsize=12)

        left += proportions

    ax.set_yticks(y_positions, labels=y_labels)
    ax.set_xlabel("Selection proportion (%)" if show_xlabel else "", fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n_levels - 0.5)
    ax.margins(y=0)
    if show_xlabel:
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
        ax.tick_params(axis="x", labelsize=14)
    else:
        ax.set_xticks([])
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
    ax.tick_params(axis="y", labelsize=16)

    if invert_y:
        ax.invert_yaxis()


def plot_eq_combined(
    counts_base: pd.DataFrame,
    counts_episode: pd.DataFrame,
    counts_noise: pd.DataFrame,
    eq_order: list[str],
    episode_order: list[str],
    noise_order: list[str],
    out_path: Path,
) -> None:
    """Create the three-panel EQ composition figure (baseline / episode / noise)."""
    eq_colors = plt.cm.Set2.colors
    eq_color_map = {eq: eq_colors[i % len(eq_colors)] for i, eq in enumerate(eq_order)}
    eq_color_map["Flat"] = "#b0b0b0"

    # Match the original publication layout geometry.
    bar_height_in = 0.7
    legend_h_in = 0.6
    panel_gap_in = 0.35
    fig_width = 15.0
    left_margin_in = 3.1
    right_margin_in = 0.3
    top_margin_in = 0.15
    bottom_margin_in = 1.15

    panel_h = [
        bar_height_in * 1,
        bar_height_in * len(episode_order),
        bar_height_in * len(noise_order),
    ]
    fig_height = legend_h_in + sum(panel_h) + 2 * panel_gap_in + bottom_margin_in + top_margin_in

    fig = plt.figure(figsize=(fig_width, fig_height))
    axes_left = left_margin_in / fig_width
    axes_width = 1 - (left_margin_in + right_margin_in) / fig_width

    y_top = 1 - top_margin_in / fig_height
    legend_h_norm = legend_h_in / fig_height
    ax_legend = fig.add_axes([axes_left, y_top - legend_h_norm, axes_width, legend_h_norm])
    y_cursor = y_top - legend_h_norm - panel_gap_in / fig_height

    base_h_norm = panel_h[0] / fig_height
    ax_base = fig.add_axes([axes_left, y_cursor - base_h_norm, axes_width, base_h_norm])
    y_cursor -= base_h_norm + panel_gap_in / fig_height

    episode_h_norm = panel_h[1] / fig_height
    ax_episode = fig.add_axes([axes_left, y_cursor - episode_h_norm, axes_width, episode_h_norm])
    y_cursor -= episode_h_norm + panel_gap_in / fig_height

    noise_h_norm = panel_h[2] / fig_height
    ax_noise = fig.add_axes([axes_left, y_cursor - noise_h_norm, axes_width, noise_h_norm])
    ax_legend.axis("off")

    draw_eq_stacked(
        ax=ax_base,
        counts=counts_base,
        row_order=["All"],
        eq_order=eq_order,
        eq_color_map=eq_color_map,
        y_labels=["Perceived Naturalness\nBaseline"],
        show_xlabel=False,
        invert_y=False,
    )

    draw_eq_stacked(
        ax=ax_episode,
        counts=counts_episode,
        row_order=episode_order,
        eq_order=eq_order,
        eq_color_map=eq_color_map,
        y_labels=[f"{e} Episode" for e in episode_order],
        show_xlabel=False,
        invert_y=True,
    )

    draw_eq_stacked(
        ax=ax_noise,
        counts=counts_noise,
        row_order=noise_order,
        eq_order=eq_order,
        eq_color_map=eq_color_map,
        y_labels=[NOISE_LABEL_MAP.get(n, n) for n in noise_order],
        show_xlabel=True,
        invert_y=False,
    )

    handles, labels = ax_base.get_legend_handles_labels()
    ax_legend.legend(
        handles,
        labels,
        frameon=False,
        ncol=len(eq_order),
        loc="center right",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=16,
    )

    fig.savefig(out_path)
    plt.close(fig)


def plot_level_combined(df: pd.DataFrame, episode_order: list[str], noise_order: list[str], out_path: Path) -> None:
    """Create the two-panel violin+box figure for playback level (dB)."""
    fig, (ax_episode, ax_noise) = plt.subplots(
        1,
        2,
        figsize=(12, 5),
        sharey=True,
        gridspec_kw={"width_ratios": [2, 1]},
    )

    # Episode panel
    level_by_episode = [df.loc[df["episode"] == ep, "selected_level_db"] for ep in episode_order]
    episode_violin = ax_episode.violinplot(level_by_episode, showmedians=False, showextrema=False, points=400, bw_method=0.2)
    for body in episode_violin["bodies"]:
        body.set_facecolor("#b0b0b0")
        body.set_alpha(0.8)
    ax_episode.boxplot(level_by_episode, widths=0.2, medianprops={"color": "black", "linewidth": 1.2})
    ax_episode.set_xticks(range(1, len(episode_order) + 1), labels=episode_order)
    ax_episode.set_xlabel("Listening Episode", fontsize=16)
    ax_episode.set_ylabel("Adjusted playback level (dB)", fontsize=16)
    ax_episode.set_yticks(range(-24, 7, 6))
    ax_episode.tick_params(axis="x", labelsize=12)
    ax_episode.tick_params(axis="y", labelsize=12)
    ax_episode.grid(axis="y", linestyle="--", alpha=0.5)

    # Noise panel
    level_by_noise = [df.loc[df["noise"] == noise, "selected_level_db"] for noise in noise_order]
    noise_violin = ax_noise.violinplot(level_by_noise, showmedians=False, showextrema=False, points=400, bw_method=0.2)
    for body in noise_violin["bodies"]:
        body.set_facecolor("#b0b0b0")
        body.set_alpha(0.8)
    ax_noise.boxplot(level_by_noise, widths=0.2, medianprops={"color": "black", "linewidth": 1.2})
    ax_noise.set_xticks(range(1, len(noise_order) + 1), labels=noise_order)
    ax_noise.set_xlabel("Automotive Engine Noise", fontsize=16)
    ax_noise.tick_params(axis="x", labelsize=12)
    ax_noise.tick_params(axis="y", labelsize=12)
    ax_noise.grid(axis="y", linestyle="--", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_eq_noise_episode_grid(
    df: pd.DataFrame,
    eq_order: list[str],
    episode_order: list[str],
    noise_order: list[str],
    out_path: Path,
) -> None:
    """Create side-by-side panels (EV/DIESEL), with episodes on x-axis and stacked bars."""
    eq_colors = plt.cm.Set2.colors
    eq_color_map = {eq: eq_colors[i % len(eq_colors)] for i, eq in enumerate(eq_order)}
    eq_color_map["Flat"] = "#b0b0b0"

    # Proportions per (noise, episode).
    counts = (
        df.groupby(["noise", "episode", "selected_eq"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=eq_order, fill_value=0)
    )

    fig, axes = plt.subplots(
        1,
        len(noise_order),
        figsize=(12, 6.8),
        sharey=True,
    )
    if len(noise_order) == 1:
        axes = [axes]

    x_positions = list(range(len(episode_order)))
    episode_tick_labels = [ep.replace("Focus-Motivation", "Focus-\nMotivation") for ep in episode_order]
    ytick_values = [0, 0.25, 0.5, 0.75, 1.0]
    ytick_labels = ["0%", "25%", "50%", "75%", "100%"]

    for col_idx, noise in enumerate(noise_order):
        ax = axes[col_idx]
        bottoms = pd.Series(0.0, index=episode_order)

        for eq in eq_order:
            heights = []
            for episode in episode_order:
                if (noise, episode) in counts.index:
                    row_counts = counts.loc[(noise, episode)]
                else:
                    row_counts = pd.Series(0, index=eq_order)
                total = float(row_counts.sum())
                p = float(row_counts[eq] / total) if total > 0 else 0.0
                heights.append(p)

            bars = ax.bar(
                x_positions,
                heights,
                width=0.72,
                bottom=bottoms.values,
                color=eq_color_map[eq],
                label=eq,
            )

            for x_idx, (bar, p) in enumerate(zip(bars, heights)):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottoms.iloc[x_idx] + p / 2,
                    f"{p * 100:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=10,
                )

            bottoms += pd.Series(heights, index=episode_order)

        ax.set_title(NOISE_LABEL_MAP.get(noise, noise), fontsize=17, pad=8)
        ax.set_xticks(x_positions, labels=episode_tick_labels)
        ax.tick_params(axis="x", labelsize=13, rotation=0, color="#666666")
        ax.set_ylim(0, 1)
        ax.set_yticks(ytick_values)
        if col_idx == 0:
            ax.set_ylabel("EQ selection proportion (%)", fontsize=16)
            ax.set_yticklabels(ytick_labels)
            ax.tick_params(axis="y", labelsize=14, labelleft=True, left=True, color="#666666")
        else:
            ax.tick_params(axis="y", labelleft=False, left=False, color="#666666")
        # Remove box frame while keeping x/y axes.
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_color("#666666")
        ax.spines["bottom"].set_color("#666666")
        ax.grid(axis="y", linestyle="--", alpha=0.35)

    handles = [plt.Rectangle((0, 0), 1, 1, color=eq_color_map[eq]) for eq in eq_order]
    fig.legend(
        handles,
        eq_order,
        frameon=False,
        ncol=len(eq_order),
        loc="upper right",
        bbox_to_anchor=(0.99, 0.995),
        fontsize=15,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path)
    plt.close(fig)


# ------------------------------
# Main Script
# ------------------------------
def main() -> None:
    df = read_csv_with_fallback(DATA_FILE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Keep only rows used by descriptive figures.
    df = df.dropna(subset=["selected_level_db", "noise", "episode", "selected_eq", "selected_system_eq"]).copy()

    # Preserve manuscript order while keeping only observed levels.
    eq_order = [eq for eq in EQ_ORDER_DEFAULT if eq in set(df["selected_eq"]) or eq in set(df["selected_system_eq"])]
    episode_order = [ep for ep in EPISODE_ORDER_DEFAULT if ep in set(df["episode"])]
    noise_order = [noise for noise in NOISE_ORDER_DEFAULT if noise in set(df["noise"])]

    # EQ counts for baseline / episode / noise panels.
    counts_base = pd.DataFrame(
        [df["selected_system_eq"].value_counts().reindex(eq_order, fill_value=0).values],
        index=["All"],
        columns=eq_order,
    )
    counts_episode = (
        df.groupby(["episode", "selected_eq"]).size().unstack(fill_value=0).reindex(index=episode_order, columns=eq_order, fill_value=0)
    )
    counts_noise = (
        df.groupby(["noise", "selected_eq"]).size().unstack(fill_value=0).reindex(index=noise_order, columns=eq_order, fill_value=0)
    )

    plot_eq_combined(
        counts_base=counts_base,
        counts_episode=counts_episode,
        counts_noise=counts_noise,
        eq_order=eq_order,
        episode_order=episode_order,
        noise_order=noise_order,
        out_path=OUTPUT_DIR / "descriptive_eq.pdf",
    )

    plot_level_combined(
        df=df,
        episode_order=episode_order,
        noise_order=noise_order,
        out_path=OUTPUT_DIR / "descriptive_level.pdf",
    )

    plot_eq_noise_episode_grid(
        df=df,
        eq_order=eq_order,
        episode_order=episode_order,
        noise_order=noise_order,
        out_path=OUTPUT_DIR / "descriptive_eq_grid.pdf",
    )


if __name__ == "__main__":
    main()
