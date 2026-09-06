
import csv
import itertools
import os
import random
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import mpl_toolkits.axes_grid1.inset_locator
import networkx as nx
import numpy as np
import pandas as pd
import scipy as sci
import scipy.optimize as sciopt
import seaborn as sns

import statsmodels.stats.multitest as smsm
from statannotations.Annotator import Annotator 
from sklearn.metrics import r2_score

# ==========================================
# 1. 超参数设置区 (以后修改列名只改这里)
# ==========================================
key_micro = 's50'               # X轴：微生物数据对应的列名 g130 
micro_label = 'Clostridium fessum'     # X轴：图表显示的标签名字

target_col = 'sbp_c'             # Y轴：临床指标列名 (如 sbp_b, dbp_a 等)
target_label = 'SBP (T3)'        # Y轴：图表显示的标签名字


# ==========================================
# 2. 数据读取与处理
# ==========================================
cogcat_df = pd.read_excel("../Anti_data/THSBC_Taxonomy_s_arcsinZ_ALLBP_725.xlsx", index_col=0) #THSBC_Taxonomy_s_arcsinZ_ALLBP_725 THSBC_ITS_arcsinZ_ALLBP_725
cogcat_df = cogcat_df.loc[cogcat_df['period'] == 'T1'].copy()
df = cogcat_df.copy()
df['Type'] = df['Type'].replace('Non-HDP', 'NP')

# 使用超参数过滤缺失值
df = df.dropna(subset=[target_col, 'Type', key_micro])

# 【非常重要】：重置索引，防止画 kde 分布图时报 duplicate labels 的错！
df = df.reset_index(drop=True)

min_value = df[key_micro].min()
max_value = df[key_micro].max()
print(f"{key_micro} 列的取值范围是: {min_value} 到 {max_value}")

# ==========================================
# 3. 绘图 (Jointplot)
# ==========================================
g = sns.jointplot(
    data = df,
    kind="scatter",
    x=key_micro,
    y=target_col,                # 使用超参数
    hue="Type",
    alpha=0.5,
    joint_kws=dict(s=60),
    marginal_kws={'common_norm': False},
    palette=[ "#7a94ba", "#B5665D", "#f8c672"],
    hue_order=["NP", "PE", "GH"],
    height=6
)

g.fig.set_size_inches(6.5, 6)

g.ax_joint.grid(which="major", axis="both", linestyle="--", zorder=0)

g.ax_joint.set(
    xlim=(-1.5, None)
)

# 使用超参数设置标签
g.ax_joint.set_xlabel(micro_label, fontsize=20, fontstyle='italic')
g.ax_joint.set_ylabel(target_label, fontsize=20)

# 调整图例
new_leg = {
    key: "{} (N={})".format(key, value)
    for key, value in df["Type"].value_counts().to_dict().items()
}
handles, labels = g.ax_joint.get_legend_handles_labels()
_ = g.ax_joint.legend(handles=handles, labels=labels, fontsize=18, columnspacing=0.8, loc="upper center", ncol=3)

# 修复顶部密度图的坐标轴显示瑕疵
g.ax_marg_x.ticklabel_format(axis="y", style="plain")


# ==========================================
# 4. 统计计算与回归线拟合
# ==========================================
def lin_law(x, a, b):
    return a * x + b

# 计算相关系数 (使用超参数 target_col)
r_value, p_value = stats.pearsonr(df[key_micro], df[target_col])


p_value = 0.026956 #0.026956  0.03225
beta = 0.69106  # 0.69106 0.49029
r_squared = r_value**2  

print(f"Overall Pearson correlation coefficient (r): {r_value:.4f}")
print(f"p-value: {p_value:.4e}")
print(f"R² value: {r_squared:.4f}")

# 拟合回归线 (使用超参数 target_col)
popt, pcov = sciopt.curve_fit(
    lin_law,
    df[key_micro],  
    df[target_col]  
)

# 画回归线
xx = np.linspace(df[key_micro].min(), df[key_micro].max(), 100)  
g.ax_joint.plot(xx, lin_law(xx, *popt), color="r", linewidth=2, linestyle="--") 


# ==========================================
# 5. 图形细节微调与保存
# ==========================================
if p_value < 0.001:
    p_value_str = "< 0.001"
else:
    p_value_str = f"{p_value:.3f}"

# 添加统计信息框
g.ax_joint.annotate(f"$\\beta$ = {beta:.3f}\n$\\mathit{{P}}$ = {p_value_str}",
                   xy=(0.6, 0.17),
                   xycoords='axes fraction',
                   ha='left',
                   va='top',
                   bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7),
                   fontsize=18)

# 隐藏顶部边缘图的X轴（让图更美观）
g.ax_marg_x.set_visible(False)
g.ax_joint.tick_params(axis='both', which='major', labelsize=18)


save_path = f"Anti_graph/scatter_{target_col}_{key_micro}.pdf"
plt.savefig(save_path, format='pdf', bbox_inches='tight')
print(f"Plot successfully saved to: {save_path}")
