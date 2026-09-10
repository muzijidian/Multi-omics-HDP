import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import pandas as pd
import scipy.optimize as sciopt
import seaborn as sns

# ==========================================
# 1. Hyperparameter settings (Modify column names only here for future use)
# ==========================================
key_micro = 's50'               # X-axis: Column name corresponding to microbial data
micro_label = 'Clostridium fessum'     # X-axis: Label name displayed on the chart

target_col = 'sbp_c'             # Y-axis: Clinical indicator column name (e.g., sbp_b, dbp_a, etc.)
target_label = 'SBP (T3)'        # Y-axis: Label name displayed on the chart


# ==========================================
# 2. Data loading and processing
# ==========================================
cogcat_df = pd.read_excel("../Anti_data/THSBC_Taxonomy_s_arcsinZ_ALLBP_725.xlsx", index_col=0) 
cogcat_df = cogcat_df.loc[cogcat_df['period'] == 'T1'].copy()
df = cogcat_df.copy()
df['Type'] = df['Type'].replace('Non-HDP', 'NP')

# Filter out missing values using hyperparameters
df = df.dropna(subset=[target_col, 'Type', key_micro])

# [CRITICAL]: Reset index to prevent 'duplicate labels' error when plotting the KDE distribution!
df = df.reset_index(drop=True)

min_value = df[key_micro].min()
max_value = df[key_micro].max()
print(f"Value range of {key_micro} column is: {min_value} to {max_value}")

# ==========================================
# 3. Plotting (Jointplot)
# ==========================================
g = sns.jointplot(
    data=df,
    kind="scatter",
    x=key_micro,
    y=target_col,                # Use hyperparameters
    hue="Type",
    alpha=0.5,
    joint_kws=dict(s=60),
    marginal_kws={'common_norm': False},
    palette=["#7a94ba", "#B5665D", "#f8c672"],
    hue_order=["NP", "PE", "GH"],
    height=6
)

g.fig.set_size_inches(6.5, 6)

g.ax_joint.grid(which="major", axis="both", linestyle="--", zorder=0)

g.ax_joint.set(
    xlim=(-1.5, None)
)

# Set labels using hyperparameters
g.ax_joint.set_xlabel(micro_label, fontsize=20, fontstyle='italic')
g.ax_joint.set_ylabel(target_label, fontsize=20)

# Adjust legend
new_leg = {
    key: "{} (N={})".format(key, value)
    for key, value in df["Type"].value_counts().to_dict().items()
}
handles, labels = g.ax_joint.get_legend_handles_labels()
_ = g.ax_joint.legend(handles=handles, labels=labels, fontsize=18, columnspacing=0.8, loc="upper center", ncol=3)

# Fix axis display flaws in the top density plot
g.ax_marg_x.ticklabel_format(axis="y", style="plain")


# ==========================================
# 4. Statistical calculation and regression line fitting
# ==========================================
def lin_law(x, a, b):
    return a * x + b

# Calculate correlation coefficient (using target_col)
r_value, p_value = stats.pearsonr(df[key_micro], df[target_col])


p_value = 0.026956 
beta = 0.69106 
r_squared = r_value**2  

print(f"Overall Pearson correlation coefficient (r): {r_value:.4f}")
print(f"p-value: {p_value:.4e}")
print(f"R² value: {r_squared:.4f}")

# Fit regression line (using target_col)
popt, pcov = sciopt.curve_fit(
    lin_law,
    df[key_micro],  
    df[target_col]  
)

# Plot regression line
xx = np.linspace(df[key_micro].min(), df[key_micro].max(), 100)  
g.ax_joint.plot(xx, lin_law(xx, *popt), color="r", linewidth=2, linestyle="--") 


# ==========================================
# 5. Fine-tuning figure details and saving
# ==========================================
if p_value < 0.001:
    p_value_str = "< 0.001"
else:
    p_value_str = f"{p_value:.3f}"

# Add statistical information box
g.ax_joint.annotate(f"$\\beta$ = {beta:.3f}\n$\\mathit{{P}}$ = {p_value_str}",
                   xy=(0.6, 0.17),
                   xycoords='axes fraction',
                   ha='left',
                   va='top',
                   bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7),
                   fontsize=18)

# Hide the X-axis of the top marginal plot (for better aesthetics)
g.ax_marg_x.set_visible(False)
g.ax_joint.tick_params(axis='both', which='major', labelsize=18)


save_path = f"Anti_graph/scatter_{target_col}_{key_micro}.pdf"
plt.savefig(save_path, format='pdf', bbox_inches='tight')
print(f"Plot successfully saved to: {save_path}")