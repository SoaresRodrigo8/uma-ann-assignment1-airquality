import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
import jinja2

path = "C:/Users/Utilizador/Documents/GitHub/uma-ann-assignment1-airquality/data/Raw/AirQualityUCI.csv"

df = pd.read_csv(path, sep=";", decimal=",")

df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
df = df.dropna(how="all")
df = df.replace(-200, np.nan)

ts = pd.to_datetime(df["Date"] + " " + df["Time"], format="%d/%m/%Y %H.%M.%S")

# 1. Duplicate timestamps
print("Duplicates:", ts.duplicated().sum())

# 2. Is it in chronological order?
print("Sorted:", ts.is_monotonic_increasing)

# 3. Compare against a complete hourly range
full_range = pd.date_range(ts.min(), ts.max(), freq="h")
missing = full_range.difference(ts)
print(f"Range: {ts.min()} to {ts.max()}")
print(f"Expected hours: {len(full_range)}, present: {ts.nunique()}, missing: {len(missing)}")

# 4. Step sizes between consecutive rows
steps = ts.diff().value_counts()
print(steps)

num = df.select_dtypes("number").columns

counts = df[num].isna().sum()

table = pd.DataFrame({
    "Missing": counts,
    "% Missing": (counts / len(df) * 100).round(2),
}).sort_values("Missing", ascending=False)

#print(f"\nTotal missing values (after): {table['Missing'].sum()}")
#print(table)

sns.set_theme(style="whitegrid", context="paper", font_scale=3)
color = "#2E6F9E"

target = df["CO(GT)"].dropna()
mean, median = target.mean(), target.median()

fig, (ax_hist, ax_box) = plt.subplots(
    2, 1, figsize=(8, 5), sharex=True,
    gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}
)

# Histogram: bin edges at x.x5 so each bin holds exactly two 0.1-step values
sns.histplot(target, binwidth=0.2, binrange=(0.05, target.max() + 0.2),
             kde=True, color=color, edgecolor="white", linewidth=0.5,
             alpha=0.8, ax=ax_hist)
ax_hist.lines[0].set_color("#1B3F5C")
ax_hist.lines[0].set_linewidth(2)

ax_hist.axvline(mean, color="#C0392B", ls="--", lw=1.5, label=f"Mean = {mean:.2f}")
ax_hist.axvline(median, color="#E67E22", ls=":", lw=2, label=f"Median = {median:.2f}")
ax_hist.legend(frameon=False)
ax_hist.set_ylabel("Count")
ax_hist.set_title("Distribution of CO(GT)", fontweight="bold")

# Boxplot
sns.boxplot(x=target, color=color, width=0.5, linewidth=1.2,
            flierprops=dict(marker="o", markersize=3, alpha=0.4,
                            markerfacecolor=color, markeredgecolor="none"),
            ax=ax_box)
ax_box.set_xlabel("CO concentration (mg/m³)")
ax_box.set_yticks([])

sns.despine(ax=ax_hist)
sns.despine(ax=ax_box, left=True)

plt.savefig("co_distribution.pdf", bbox_inches="tight")
plt.show()

target_col = "CO(GT)"
feature_cols = [c for c in df.select_dtypes("number").columns
                if c not in [target_col, "NMHC(GT)"]]   # NMHC dropped (90% missing)

cols = [target_col] + feature_cols
stats = df[cols].describe().T
stats["skew"] = df[cols].skew()
stats["kurtosis"] = df[cols].kurt()
stats = stats.rename(columns={"25%": "Q1", "50%": "median", "75%": "Q3"})
stats = stats[["mean", "std", "min", "median", "max", "skew"]]

#print(stats.round(2))

# For LaTeX
#print(stats.round(2).to_latex(float_format="%.2f",
     # caption="Descriptive statistics of the target and features.",
      #label="tab:descriptive_stats"))

n = len(cols)
ncols = 4
nrows = int(np.ceil(n / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.2 * nrows))
for ax, col in zip(axes.flat, cols):
    sns.histplot(df[col].dropna(), bins=40, kde=True, color=color,
                 edgecolor="white", linewidth=0.4, ax=ax)
    ax.set_title(col, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
for ax in axes.flat[n:]:
    ax.set_visible(False)

sns.despine()
plt.tight_layout()
plt.savefig("feature_distributions.pdf", bbox_inches="tight")
plt.show()

groups = {
    "Strongly right-skewed": ["NOx(GT)", "C6H6(GT)", "PT08.S3(NOx)"],
    "Moderately right-skewed": ["PT08.S1(CO)", "PT08.S2(NMHC)", "PT08.S5(O3)", "NO2(GT)"],
    "Approximately symmetric": ["PT08.S4(NO2)", "T", "RH", "AH"],
}

palette = ["#2E6F9E", "#E67E22", "#27AE60", "#C0392B"]
linestyles = ["-", "--", "-.", ":"]

z = (df[sum(groups.values(), [])] - df[sum(groups.values(), [])].mean()) \
    / df[sum(groups.values(), [])].std()

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True, sharey=True)

for ax, (title, feats) in zip(axes, groups.items()):
    for col, c, ls in zip(feats, palette, linestyles):
        sns.kdeplot(z[col].dropna(), ax=ax, color=c, linestyle=ls, lw=2.2,
                    fill=True, alpha=0.08, cut=0,
                    label=f"{col} (skew = {df[col].skew():.2f})")
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Standardised value (z-score)")
    ax.legend(frameon=False, fontsize=10, loc="upper right")

axes[0].set_ylabel("Density")
axes[0].set_xlim(-3, 6)
sns.despine()
plt.tight_layout()
plt.savefig("feature_distributions_grouped.pdf", bbox_inches="tight")
plt.show()

corr_target = df[cols].corr(method="spearman")[target_col].drop(target_col) \
                .sort_values(key=abs, ascending=False)

ax = corr_target.plot.barh(figsize=(8, 6), color=color, width=0.7)
ax.axvline(0, color="black", lw=0.8)
ax.invert_yaxis()
ax.set_xlabel("Spearman correlation with CO(GT)")
ax.set_title("Feature correlation with the target", fontweight="bold")
sns.despine()
plt.tight_layout()
plt.savefig("corr_with_target.pdf", bbox_inches="tight")
plt.show()

num = df.select_dtypes("number").columns
q1, q3 = df[num].quantile(0.25), df[num].quantile(0.75)
iqr = q3 - q1
outliers = ((df[num] < q1 - 1.5 * iqr) | (df[num] > q3 + 1.5 * iqr)).sum()
print(pd.DataFrame({"Outliers": outliers,
                    "%": (outliers / df[num].notna().sum() * 100).round(2)}))

upper = (df[num] > q3 + 1.5 * iqr).sum()
lower = (df[num] < q1 - 1.5 * iqr).sum()
print(pd.DataFrame({"upper": upper, "lower": lower}))
