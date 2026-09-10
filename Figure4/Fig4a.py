from pycirclize import Circos
import numpy as np
import pandas as pd
from collections import OrderedDict
from collections import Counter 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch


def get_vmin_vmax(sheet):
    all_y = sheet.loc[:, 'beta']
    all_y = all_y.sort_values(ascending=False)
    abs_max_y = max(abs(all_y))
    vmin, vmax = -abs_max_y, abs_max_y
    return vmin, vmax


def draw_circos_line(sec_n, deg, text, r=(30, 101), text_r=105, color='k', line=False):
    text_common_kws = dict(ha="left", va="center", size=7, color=color, fontname='Arial', orientation='vertical', adjust_rotation=True)
    if line:
        circos.line(r=r, deg_lim=[deg], ls='-', lw=0.5)
    circos.text(text, r=text_r, deg=deg, fontweight='semibold', **text_common_kws)


class SplitPatch(Patch):
    def __init__(self, *args, **kwargs):
        self.color1 = kwargs.pop('color1', 'white')
        self.color2 = kwargs.pop('color2', 'white')
        self.edgecolor = kwargs.pop('edgecolor', 'black')
        super().__init__(*args, **kwargs)

class SplitPatchHandler(HandlerPatch):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        r1 = Rectangle([x0, y0], width / 2, height, facecolor='#c00000', edgecolor='none')
        r2 = Rectangle([x0 + width / 2, y0], width / 2, height, facecolor='#002060', edgecolor='none')
        return [r1, r2]

        
# Load data
PE_T1 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='PE_T1', index_col='meta_name')
PE_T2 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='PE_T2', index_col='meta_name')
PE_T3 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='PE_T3', index_col='meta_name')
GH_T1 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='GH_T1', index_col='meta_name')
GH_T2 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='GH_T2', index_col='meta_name')
GH_T3 = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='GH_T3', index_col='meta_name')
PE_Pooled = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='Pooled (PE)', index_col='meta_name')
GH_Pooled = pd.read_excel('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet_name='Pooled (GH)', index_col='meta_name')

def remove_row(df):
    return df.drop(index='MEDN1125', errors='ignore')  # Ignore error if MEDN1125 does not exist

# Apply row deletion to each dataframe
PE_T1 = remove_row(PE_T1)
PE_T2 = remove_row(PE_T2)
PE_T3 = remove_row(PE_T3)
GH_T1 = remove_row(GH_T1)
GH_T2 = remove_row(GH_T2)
GH_T3 = remove_row(GH_T3)
PE_Pooled = remove_row(PE_Pooled)
GH_Pooled = remove_row(GH_Pooled)


def format_metabolite_name(name):
    # Handle 'acid/Acid' capitalization
    if ' acid' in name.lower():
        name = name.replace(' Acid', ' acid').replace(' ACID', ' acid')
    
    # Map specific metabolite names
    specific_names = {
        'trans-resveratrol-3-O-sulfate' : 'Trans-resveratrol-3-O-sulfate',
        'Phosphatidylethanolamine lyso alkenyl 16:0': 'LPE(P-16:0)',
        'Glycine deoxycholic acid': 'GDCA',
        'Glycochenodeoxycholic acid': 'GCDCA',
        'Deoxycholic acid': 'DCA',
        'Taurocholic acid': 'TCA',
        'Taurochenodesoxycholic acid': 'TCDCA',
        'Glycoursodeoxycholic acid': 'GUDCA',
        'Glycohyodeoxycholic acid': 'GHDCA'
    }
    
    return specific_names.get(name, name)

# Modify compound names in all dataframes
dfs = [PE_T1, PE_T2, PE_T3, GH_T1, GH_T2, GH_T3, PE_Pooled, GH_Pooled]
for df in dfs:
    if 'sugg_cmpd_name' in df.columns:
        df['sugg_cmpd_name'] = df['sugg_cmpd_name'].apply(format_metabolite_name)


PE_T1 = PE_T1.sort_values(by=['Class.I'])
PE_T2 = PE_T2.loc[PE_T1.index]
PE_T3 = PE_T3.loc[PE_T1.index]
GH_T1 = GH_T1.loc[PE_T1.index]
GH_T2 = GH_T2.loc[PE_T1.index]
GH_T3 = GH_T3.loc[PE_T1.index]
PE_Pooled = PE_Pooled.loc[PE_T1.index]
GH_Pooled = GH_Pooled.loc[PE_T1.index]


colorset = ("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")
order = ["AlAm", "AAD", "BA", "BSD", "CHD", "FA", "GP", "HRC", "NUC", "OAD", "SPL", "Others"]
color_mapping = dict(zip(order, colorset))
sectors = dict(Counter(PE_T1.loc[:, 'Class.I']))
sectors = OrderedDict(sorted(sectors.items(), key=lambda x: order.index(x[0])))

fig = plt.figure(figsize=(12, 12))
circos = Circos(sectors, space=2, start=0, end=270)


for i, sector in enumerate(circos.sectors):
    
    sector_start = sector.start
    sector_end = sector.end
    deg_per_point = sector.deg_size / sector.size
    base_deg = sector.deg_lim[0] + deg_per_point/2
    
    outer_track = sector.add_track((90, 100))
    outer_track.axis(fc=color_mapping[sector.name])
    outer_track.text(sector.name, size=9)

    global_data = PE_T1[PE_T1.loc[:, 'Class.I']==sector.name]
    gobal_y = global_data.loc[:, 'beta']
    gobal_y = gobal_y.sort_values(ascending=False)

    # Ring 2: Proteins distribution
    ProteinsDis_track = sector.add_track((84.5, 89.5), r_pad_ratio=0.8)
    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)
    ProteinsDis_track.axis(fc='#fbffff', ec="thistle")

    # Filter PE_T1 and GH_T1 data for the current sector
    PE_T1_data = PE_T1[PE_T1.loc[:, 'Class.I'] == sector.name].reindex(gobal_y.index)
    GH_T1_data = GH_T1[GH_T1.loc[:, 'Class.I'] == sector.name].reindex(gobal_y.index)

    # Calculate if any significance in T1 for PE or GH (pv_adj_FDR < 0.05)
    PE_T1_significant = PE_T1_data['pv_adj_FDR'] < 0.05
    GH_T1_significant = GH_T1_data['pv_adj_FDR'] < 0.05
    significant_in_T1 = (PE_T1_significant.fillna(False) | GH_T1_significant.fillna(False))

    colors = ['firebrick' if sig else 'orange' for sig in significant_in_T1]

    # Use color list directly when plotting to completely avoid colormap interference
    ProteinsDis_track.scatter(x, y, c=colors, s=13)


    # Ring 5: PE T1 FDR
    PE_T1_fdr_track = sector.add_track((80, 84), r_pad_ratio=0.8)
    PE_T1_fdr_track.axis(fc='aliceblue', ec="aliceblue")
    
    colormap = ['seagreen', 'purple'] 
    markers = ['^', 'x']

    # Get values for symbols
    PE_T1_data = PE_T1[PE_T1.loc[:, 'Class.I']==sector.name]
    colorvalue = PE_T1_data.loc[:, 'pv_adj_FDR']
    colorvalue = colorvalue.reindex(gobal_y.index)

    # Define x and y positions
    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)

    for j in range(len(markers)):
        # Category 1: colorvalue < 0.01
        if j == 0:
            mask = colorvalue < 0.01
        # Category 2: 0.01 <= colorvalue < 0.05
        else:
            mask = (colorvalue >= 0.01) & (colorvalue < 0.05)

        vx_subset = x[mask]
        vy_subset = y[mask]

        if vx_subset.size == 0:
            continue

        # Draw scatter plot
        PE_T1_fdr_track.scatter(vx_subset, vy_subset, c=colormap[j], marker=markers[j], s=18)

    
    # Ring 6: PE T1 beta
    PE_T1_beta_track = sector.add_track((68.5, 79.5))
    x = np.arange(sector.start, sector.end) + 0.5
    
    PE_T1_beta = PE_T1[PE_T1.loc[:, 'Class.I']==sector.name]
    y = PE_T1_beta.loc[:, 'beta']
    y = y.reindex(gobal_y.index)

    abs_max_y = get_vmin_vmax(PE_T1)[1]
    vmin, vmax = -abs_max_y, abs_max_y
    
    PE_T1_beta_track.axis(fc='#f6fffc', ec='black')
    PE_T1_beta_track.bar(x, y, color=np.where(y>0, '#c00000', '#002060'), vmin=vmin, vmax=vmax)


    # Ring 7: GH T1 FDR
    GH_T1_fdr_track = sector.add_track((64, 68), r_pad_ratio=0.8)
    GH_T1_fdr_track.axis(fc='aliceblue', ec="aliceblue")
    
    # Get values for symbols
    GH_T1_data = GH_T1[GH_T1.loc[:, 'Class.I']==sector.name]
    colorvalue = GH_T1_data.loc[:, 'pv_adj_FDR']
    colorvalue = colorvalue.reindex(gobal_y.index)

    # Define x and y positions
    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)

    # Draw symbols based on colorvalue
    for j in range(len(markers)):
        if j == 0:
            mask = colorvalue < 0.01
        else:
            mask = (colorvalue >= 0.01) & (colorvalue < 0.05)

        vx_subset = x[mask]
        vy_subset = y[mask]

        if vx_subset.size == 0:
            continue

        GH_T1_fdr_track.scatter(vx_subset, vy_subset, c=colormap[j], marker=markers[j], s=18)
    

    # Ring 8: GH T1 beta
    GH_T1_beta_track = sector.add_track((52.5, 63.5))
    x = np.arange(sector.start, sector.end) + 0.5
    
    GH_T1_beta = GH_T1[GH_T1.loc[:, 'Class.I']==sector.name]
    y = GH_T1_beta.loc[:, 'beta']
    y = y.reindex(gobal_y.index)

    abs_max_y = get_vmin_vmax(GH_T1)[1]
    vmin, vmax = -abs_max_y, abs_max_y
    
    GH_T1_beta_track.axis(fc='#fffcf7')
    GH_T1_beta_track.bar(x, y, color=np.where(y>0, '#c00000', '#002060'), vmin=vmin, vmax=vmax)
    
    
    datasets = [PE_T2, GH_T2, PE_T3, GH_T3, PE_Pooled, GH_Pooled]

    # Create colormap
    cmap = ListedColormap(['#002060', 'white', '#c00000'])

    # Initialize an empty 2D array to store the entire heatmap data
    full_heatmap_data = np.zeros((len(datasets), sector.size))

    # Set initial position
    start_position = 51.
    track_width = 4

    # Draw heatmap
    for idx, dataset in enumerate(datasets):
        data = dataset.loc[gobal_y.index]  # Ensure order alignment with global_y.index
        data = data[data.loc[:, 'Class.I'] == sector.name]
        beta = data.loc[:, 'beta']
        p_adj_FDR = data.loc[:, 'pv_adj_FDR']
        p_adj_FDR = p_adj_FDR.reindex(gobal_y.index)

        # Create an empty array for the heatmap in the current loop
        heatmap_data = np.full(sector.size, 1)

        # Set heatmap data
        heatmap_data[(p_adj_FDR < 0.05) & (beta > 0)] = 2  # Set beta > 0 to red
        heatmap_data[(p_adj_FDR < 0.05) & (beta < 0)] = 0  # Set beta < 0 to blue (0)

        # Add the heatmap data of the current loop to the full heatmap data
        full_heatmap_data[idx, :] = heatmap_data

    # Create a new track to display the complete heatmap
    inner_track = sector.add_track((start_position - len(datasets) * track_width, start_position))

    # Set axis color
    inner_track.axis(fc='#fffcf7', ec="black")

    # Draw the complete heatmap
    inner_track.heatmap(full_heatmap_data, vmin=0, vmax=2, cmap=cmap, rect_kws=dict(ec="k", lw=0.2, ls="dashed"))


    # Significance data at T1 time point
    temp_pd_raw = pd.concat([PE_T1.loc[:, 'pv_adj_FDR'], GH_T1.loc[:, 'pv_adj_FDR']], axis=1)
    temp_pd_T1_sig = temp_pd_raw < 0.05  # Significance of all groups at T1

    # Significance data at T2 and T3 time points
    significant_in_T2_T3 = pd.concat([
        pd.Series(PE_T2.loc[PE_T2['pv_adj_FDR'] < 0.05].index),
        pd.Series(GH_T2.loc[GH_T2['pv_adj_FDR'] < 0.05].index),
        pd.Series(PE_T3.loc[PE_T3['pv_adj_FDR'] < 0.05].index),
        pd.Series(GH_T3.loc[GH_T3['pv_adj_FDR'] < 0.05].index),
        pd.Series(PE_Pooled.loc[PE_Pooled['pv_adj_FDR'] < 0.05].index),
        pd.Series(GH_Pooled.loc[GH_Pooled['pv_adj_FDR'] < 0.05].index)
    ]).unique()

    # Reindex
    temp_pd = temp_pd_T1_sig.reindex(gobal_y.index)

    for j in range(len(temp_pd.index)):
        deg = base_deg + deg_per_point*j
        current_metabolite = temp_pd.index[j]
        
        # Check if at least one group is significant in T1
        is_significant_T1_any = any(temp_pd_T1_sig.loc[current_metabolite])
        # Check if significant in T2 or T3
        is_in_T2_T3 = current_metabolite in significant_in_T2_T3
        
        if is_significant_T1_any:
            if not is_in_T2_T3:
                # At least one significant in T1, and neither significant in T2/T3 -> black
                draw_circos_line(sec_n=i, deg=deg, 
                            text=PE_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='k')
            else:
                # At least one significant in T1, and also significant in T2 or T3 -> firebrick
                draw_circos_line(sec_n=i, deg=deg, 
                            text=PE_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='firebrick')
        else:
            # None significant in T1 -> gray
            draw_circos_line(sec_n=i, deg=deg, 
                            text=PE_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='gray')


common_r = 100
# Left Legend
text_common_kws = dict(ha="right", va="center", size=8.5)            
circos.text(f"{PE_T1.shape[0]} serum \n  metabolites", r=0, size=12)
circos.text("Metabolite class 1—| ", r=common_r-5, color="black", **text_common_kws)
circos.text("T1 perturbation 2—| ", r=common_r-12, color="black", **text_common_kws)
circos.text("PE (T1): Significance 3—| ", r=common_r-19, color="black", **text_common_kws)
circos.text("PE (T1): Coefficient 4—| ", r=common_r-26, color="black", **text_common_kws)
circos.text("GH (T1): Significance 5—| ", r=common_r-34., color="black", **text_common_kws)
circos.text("GH (T1): Coefficient 6—| ", r=common_r-42, color="black", **text_common_kws)


text_common_kws = dict(ha="right", va="center", size=7)
circos.text("PE (T2) —| ", r=53-4.*1, color="black", **text_common_kws)
circos.text("GH (T2) —| ", r=53-4.*2, color="black", **text_common_kws)
circos.text("PE (T3) —| ", r=53-4.*3, color="black", **text_common_kws)
circos.text("GH (T3) —| ", r=53-4.*4, color="black", **text_common_kws)
circos.text("PE (Pooled) —| ", r=53-4.*5, color="black", **text_common_kws)
circos.text("GH (Pooled) —| ", r=53-4.*6, color="black", **text_common_kws)

fig = circos.plotfig()
 
scatter_legend = circos.ax.legend(
    handles=[
        Line2D([], [], color="firebrick", marker="o", label="Perturbed in T1", ms=5, ls="None"),
        Line2D([], [], color="orange", marker="o", label="Unperturbed in T1", ms=5, ls="None"),

        Line2D([], [], color="seagreen", marker="^", label="FDR<0.01", ms=5, ls="None"),
        Line2D([], [], color="purple", marker="x", label="FDR<0.05", ms=5, ls="None"),

        Patch(color="#c00000", label="Up-regulated"),
        Patch(color="#002060", label="Down-regulated"),
        SplitPatch(label="Significant (FDR < 0.05)", color1='#c00000', color2='#002060', edgecolor="black"),
        
        Line2D([], [], color="black", marker="s", markerfacecolor='None', label="Non-significant", ms=6.5, ls="None"),
        Line2D([], [], marker=r'$\mathbf{A}$', color="none", markerfacecolor="black", markeredgecolor="black", label="T1-specific significant", ms=7, ls="None"),
        Line2D([], [], marker=r'$\mathbf{A}$', color="none", markerfacecolor="firebrick", markeredgecolor="firebrick", label="Persistent significant", ms=7, ls="None"),
        Line2D([], [], marker=r'$\mathbf{A}$', color="none", markerfacecolor="gray", markeredgecolor="gray", label="Non-significant in T1", ms=7, ls="None"),

    ],
    bbox_to_anchor=(0.05, 0.68),
    loc="center",
    ncols=1,
    fontsize=8,
    handler_map={SplitPatch: SplitPatchHandler()},
)
circos.ax.add_artist(scatter_legend)

fig.savefig(
    "new_graph/HDP_Met_FDR0.05_Circos_2026_biaopin.pdf", 
    transparent=True,
    format="pdf",          # Explicitly specify format
    bbox_inches="tight",   # Automatically crop surrounding whitespace to prevent incomplete label display
    pad_inches=0.1,        # Leave a little padding at the cropped edge
    dpi=300                # Set high DPI for rasterized elements
)