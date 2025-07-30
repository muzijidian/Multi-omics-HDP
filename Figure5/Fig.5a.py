from pycirclize import Circos
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from collections import OrderedDict
from collections import Counter 
# np.random.seed(0)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerPatch


def get_vmin_vmax(sheet):
    all_y = sheet.loc[:, 'beta']
    all_y = all_y.sort_values(ascending=False)
    abs_max_y = max(abs(all_y))
    vmin, vmax = -abs_max_y, abs_max_y
    return vmin, vmax


def draw_circos_line(sec_n, deg, text, r=(30, 101), text_r=105, color='k', line=False):
    text_common_kws = dict(ha="left", va="center", size=5, color=color, fontname='Arial', orientation='vertical', adjust_rotation=True)
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

        

HDP_T1 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='HDP_T1', index_col='meta_name')
HDP_T2 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='HDP_T2', index_col='meta_name')
HDP_T3 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='HDP_T3', index_col='meta_name')
PE_T1 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='PE_T1', index_col='meta_name')
PE_T2 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='PE_T2', index_col='meta_name')
PE_T3 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='PE_T3', index_col='meta_name')
GH_T1 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='GH_T1', index_col='meta_name')
GH_T2 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='GH_T2', index_col='meta_name')
GH_T3 = pd.read_excel('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet_name='GH_T3', index_col='meta_name')



def format_metabolite_name(name):

    if ' acid' in name.lower():
        name = name.replace(' Acid', ' acid').replace(' ACID', ' acid')

    specific_names = {
        'Glycerophospho-N-Arachidonoyl Ethanolamine': 'GP-NAE',
        'Phosphatidylethanolamine lyso alkenyl 16:0': 'LPE(P-16:0)'
    }
    
    return specific_names.get(name, name)


dfs = [HDP_T1, HDP_T2, HDP_T3, PE_T1, PE_T2, PE_T3, GH_T1, GH_T2, GH_T3]

for df in dfs:
    if 'sugg_cmpd_name' in df.columns:
        df['sugg_cmpd_name'] = df['sugg_cmpd_name'].apply(format_metabolite_name)



HDP_T1 = HDP_T1.sort_values(by=['Class.I'], )
HDP_T2 = HDP_T2.loc[HDP_T1.index]
HDP_T3 = HDP_T3.loc[HDP_T1.index]
PE_T1 = PE_T1.loc[HDP_T1.index]
PE_T2 = PE_T2.loc[HDP_T1.index]
PE_T3 = PE_T3.loc[HDP_T1.index]
GH_T1 = GH_T1.loc[HDP_T1.index]
GH_T2 = GH_T2.loc[HDP_T1.index]
GH_T3 = GH_T3.loc[HDP_T1.index]

colorset = ("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")#, "#FFD19D"
order = ["AA", "AAM", "BA", "BSD", "CHM", "FA", "GP", "HetC", "HHC", "NTM", "OAD", "SPL", "Others"]
color_mapping = dict(zip(order, colorset))
sectors = dict(Counter(HDP_T1.loc[:, 'Class.I']))
sectors = OrderedDict(sorted(sectors.items(), key=lambda x: order.index(x[0])))
# print(sectors)

sheets = [HDP_T1, HDP_T2, HDP_T3, PE_T1, PE_T2, PE_T3, GH_T1, GH_T2, GH_T3]
vminmax = [get_vmin_vmax(sheet) for sheet in sheets]




ProteinsDis = np.random.randint(0, 3, HDP_T1.shape[0])

fig = plt.figure(figsize=(12, 12))
circos = Circos(sectors, space=2, start=0, end=270)

for i, sector in enumerate(circos.sectors):
  
    sector_start = sector.start
    sector_end = sector.end
    deg_per_point = sector.deg_size / sector.size
    base_deg = sector.deg_lim[0] + deg_per_point/2
    

    outer_track = sector.add_track((97, 100))
    outer_track.axis(fc=color_mapping[sector.name], ec="w")
    outer_track.text(sector.name, size=5)

    global_data = HDP_T1[HDP_T1.loc[:, 'Class.I']==sector.name]
    gobal_y = global_data.loc[:, 'beta']
    gobal_y = gobal_y.sort_values(ascending=False)




    ProteinsDis_track = sector.add_track((93.5, 96.5), r_pad_ratio=0.8)
    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)
    ProteinsDis_track.axis(fc='thistle', ec="thistle")


    HDP_T1_data = HDP_T1[HDP_T1.loc[:, 'Class.I'] == sector.name].reindex(gobal_y.index)
    PE_T1_data = PE_T1[PE_T1.loc[:, 'Class.I'] == sector.name].reindex(gobal_y.index)
    GH_T1_data = GH_T1[GH_T1.loc[:, 'Class.I'] == sector.name].reindex(gobal_y.index)


    HDP_T1_significant = HDP_T1_data['pv_adj_FDR'] < 0.05
    PE_T1_significant = PE_T1_data['pv_adj_FDR'] < 0.05
    GH_T1_significant = GH_T1_data['pv_adj_FDR'] < 0.05


    significant_in_T1 = pd.concat([
        HDP_T1_significant, 
        PE_T1_significant, 
        GH_T1_significant
    ]).groupby(level=0).any()

    significant_in_T1 = significant_in_T1.reindex(gobal_y.index)


    ProteinsDis = np.where(significant_in_T1, 1, 0)
    ProteinsDis = pd.Series(ProteinsDis, index=HDP_T1_data.index).reindex(gobal_y.index).fillna(0).astype(int)


    colormap_P = ['orange', 'firebrick']  
    Pcmap = LinearSegmentedColormap.from_list('Protein', colormap_P)
    ProteinsDis_track.scatter(x, y, c=ProteinsDis, cmap=Pcmap)



    colormap_P = ['orange', 'k', 'firebrick']
    Pcmap = LinearSegmentedColormap.from_list('Protein', colormap_P)
    start = 0
    end = sector.size
    colorvalue = ProteinsDis[start: end]
    start = start + end
    ProteinsDis_track.scatter(x, y, c=colorvalue, cmap=Pcmap)


    HDP_T1_fdr_track = sector.add_track((88, 92), r_pad_ratio=0.8)
    HDP_T1_fdr_track.axis(fc='aliceblue', ec="aliceblue")
    
    colormap = ['gold', 'g', 'purple']
    markers = ['*', '+', 'x']


    HDP_T1_data = HDP_T1[HDP_T1.loc[:, 'Class.I']==sector.name]
    colorvalue = HDP_T1_data.loc[:, 'pv_adj_FDR']
    colorvalue = colorvalue.reindex(gobal_y.index)


    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)

    for j in range(len(markers)):
        if j == 0:
            mask = colorvalue < 0.001
        elif j == 1:
            mask = (colorvalue >= 0.001) & (colorvalue < 0.01)
        else:
            mask = (colorvalue >= 0.01) & (colorvalue < 0.05)

        vx_001 = x[mask]
        vy_001 = y[mask]

        if vx_001.size == 0:
            continue

        HDP_T1_fdr_track.scatter(vx_001, vy_001, c=colormap[j], marker=markers[j])
    

    HDP_T1_beta_track = sector.add_track((77.5, 87.5))
    x = np.arange(sector.start, sector.end) + 0.5
    
    HDP_T1_beta = HDP_T1[HDP_T1.loc[:, 'Class.I']==sector.name]
    y = HDP_T1_data.loc[:, 'beta']
    y = y.reindex(gobal_y.index)

    abs_max_y = get_vmin_vmax(HDP_T1)[1]
    vmin, vmax = -abs_max_y, abs_max_y
    
    HDP_T1_beta_track.axis(fc='lavender', ec='lavender')
    HDP_T1_beta_track.bar(x, y, color=np.where(y>0, '#c00000', '#002060'), vmin=vmin, vmax=vmax)


    PE_T1_fdr_track = sector.add_track((73, 77), r_pad_ratio=0.8)
    PE_T1_fdr_track.axis(fc='aliceblue', ec="aliceblue")
    
    colormap = ['gold', 'g', 'purple']
    markers = ['*', '+', 'x']


    PE_T1_data = PE_T1[PE_T1.loc[:, 'Class.I']==sector.name]
    colorvalue = PE_T1_data.loc[:, 'pv_adj_FDR']
    colorvalue = colorvalue.reindex(gobal_y.index)


    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)


    for j in range(len(markers)):
        if j == 0:
            mask = colorvalue < 0.001
        elif j == 1:
            mask = (colorvalue >= 0.001) & (colorvalue < 0.01)
        else:
            mask = (colorvalue >= 0.01) & (colorvalue < 0.05)

        vx_001 = x[mask]
        vy_001 = y[mask]

        if vx_001.size == 0:
            continue

        PE_T1_fdr_track.scatter(vx_001, vy_001, c=colormap[j], marker=markers[j])
    



    PE_T1_beta_track = sector.add_track((62.5, 72.5))
    x = np.arange(sector.start, sector.end) + 0.5
    
    PE_T1_beta = PE_T1[PE_T1.loc[:, 'Class.I']==sector.name]
    y = PE_T1_data.loc[:, 'beta']
    y = y.reindex(gobal_y.index)

    abs_max_y = get_vmin_vmax(PE_T1)[1]
    vmin, vmax = -abs_max_y, abs_max_y
    
    PE_T1_beta_track.axis(fc='lavender', ec='lavender')
    PE_T1_beta_track.bar(x, y, color=np.where(y>0, '#c00000', '#002060'), vmin=vmin, vmax=vmax)
    
    

    GH_T1_fdr_track = sector.add_track((58, 62), r_pad_ratio=0.8)
    GH_T1_fdr_track.axis(fc='aliceblue', ec="aliceblue")
    
    colormap = ['gold', 'g', 'purple']
    markers = ['*', '+', 'x']


    GH_T1_data = GH_T1[GH_T1.loc[:, 'Class.I']==sector.name]
    colorvalue = GH_T1_data.loc[:, 'pv_adj_FDR']
    colorvalue = colorvalue.reindex(gobal_y.index)



    x = np.arange(sector.start, sector.end) + 0.5
    y = np.zeros_like(x)

    for j in range(len(markers)):
        if j == 0:
            mask = colorvalue < 0.001
        elif j == 1:
            mask = (colorvalue >= 0.001) & (colorvalue < 0.01)
        else:
            mask = (colorvalue >= 0.01) & (colorvalue < 0.05)

        vx_001 = x[mask]
        vy_001 = y[mask]

        if vx_001.size == 0:
            continue

        GH_T1_fdr_track.scatter(vx_001, vy_001, c=colormap[j], marker=markers[j])


    GH_T1_beta_track = sector.add_track((47.5, 57.5))
    x = np.arange(sector.start, sector.end) + 0.5
    
    GH_T1_beta = GH_T1[GH_T1.loc[:, 'Class.I']==sector.name]
    y = GH_T1_data.loc[:, 'beta']
    y = y.reindex(gobal_y.index)

    abs_max_y = get_vmin_vmax(GH_T1)[1]
    vmin, vmax = -abs_max_y, abs_max_y
    
    GH_T1_beta_track.axis(fc='lavender', ec='lavender')
    GH_T1_beta_track.bar(x, y, color=np.where(y>0, '#c00000', '#002060'), vmin=vmin, vmax=vmax)


    datasets = [HDP_T2, PE_T2, GH_T2, HDP_T3, PE_T3, GH_T3]


    cmap = ListedColormap(['#002060', 'white', '#c00000'])



    full_heatmap_data = np.zeros((len(datasets), sector.size))


    start_position = 47
    track_width = 3


    for i, dataset in enumerate(datasets):
        data = dataset.loc[gobal_y.index]  
        data = data[data.loc[:, 'Class.I'] == sector.name]
        beta = data.loc[:, 'beta']
        p_adj_FDR = data.loc[:, 'pv_adj_FDR']
        p_adj_FDR = p_adj_FDR.reindex(gobal_y.index)


        heatmap_data = np.full(sector.size, 1)


        heatmap_data[(p_adj_FDR < 0.05) & (beta > 0)] = 2  
        heatmap_data[(p_adj_FDR < 0.05) & (beta < 0)] = 0  

  
        full_heatmap_data[i, :] = heatmap_data


    inner_track = sector.add_track((start_position - len(datasets) * track_width, start_position))


    inner_track.axis(fc='aliceblue', ec="aliceblue")


    inner_track.heatmap(full_heatmap_data, vmin=0, vmax=2, cmap=cmap, rect_kws=dict(ec="k", lw=0.2, ls="dashed"))



    temp_pd = pd.concat([HDP_T1.loc[:, 'pv_adj_FDR'], PE_T1.loc[:, 'pv_adj_FDR'], GH_T1.loc[:, 'pv_adj_FDR']], axis=1)
    temp_pd1 = temp_pd < 0.05  


    significant_in_T2_T3 = pd.concat([
        pd.Series(HDP_T2.loc[HDP_T2['pv_adj_FDR'] < 0.05].index),
        pd.Series(PE_T2.loc[PE_T2['pv_adj_FDR'] < 0.05].index),
        pd.Series(GH_T2.loc[GH_T2['pv_adj_FDR'] < 0.05].index),
        pd.Series(HDP_T3.loc[HDP_T3['pv_adj_FDR'] < 0.05].index),
        pd.Series(PE_T3.loc[PE_T3['pv_adj_FDR'] < 0.05].index),
        pd.Series(GH_T3.loc[GH_T3['pv_adj_FDR'] < 0.05].index)
    ]).unique()


    temp_pd = temp_pd1.reindex(gobal_y.index)

    for j in range(len(temp_pd.index)):
        deg = base_deg + deg_per_point*j
        current_metabolite = temp_pd.index[j]
    
        is_significant_T1_any = any(temp_pd1.loc[current_metabolite])

        is_in_T2_T3 = current_metabolite in significant_in_T2_T3
        
        if is_significant_T1_any:
            if not is_in_T2_T3:

                draw_circos_line(sec_n=i, deg=deg, 
                            text=HDP_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='firebrick', line=True)
            else:

                draw_circos_line(sec_n=i, deg=deg, 
                            text=HDP_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='k')
        else:

            draw_circos_line(sec_n=i, deg=deg, 
                            text=HDP_T1.loc[current_metabolite, 'sugg_cmpd_name'], 
                            text_r=102, color='gray')


common_r = 100

text_common_kws = dict(ha="right", va="center", size=8)            
circos.text(f"{HDP_T1.shape[0]} serum \n  metabolites", r=0, size=12)
circos.text("Thirteen types of serum metabolites 1—| ", r=common_r-1, color="black", **text_common_kws)
circos.text("Metabolites distribution 2—| ", r=common_r-5.5, color="black", **text_common_kws)
circos.text("HDP: T1 $\mathit{p}$-value 3—|  ", r=common_r-10.5, color="black", **text_common_kws)
circos.text("HDP: T1 coefficients 4—| ", r=common_r-17, color="black", **text_common_kws)
circos.text("PE: T1 $\mathit{p}$-value 5—|  ", r=common_r-25., color="black", **text_common_kws)
circos.text("PE: T1 coefficients 6—| ", r=common_r-32.5, color="black", **text_common_kws)
circos.text("GH: T1 $\mathit{p}$-value 7—|  ", r=common_r-40, color="black", **text_common_kws) 
circos.text("GH: T1 coefficients 8—| ", r=common_r-47, color="black", **text_common_kws)



text_common_kws = dict(ha="right", va="center", size=7)
circos.text("HDP: T2 —| ", r=49-3.*1, color="black", **text_common_kws)
circos.text("PE: T2 —| ", r=49-3.*2, color="black", **text_common_kws)
circos.text("GH: T2 —| ", r=49-3.*3, color="black", **text_common_kws)
circos.text("HDP: T3 —| ", r=49-3.*4, color="black", **text_common_kws)
circos.text("PE: T3 —| ", r=49-3.*5, color="black", **text_common_kws)
circos.text("GH: T3 —| ", r=49-3.*6, color="black", **text_common_kws)


text_common_kws = dict(ha="right", va="center", size=8)



fig = circos.plotfig()
 
scatter_legend = circos.ax.legend(
    handles=[

        Patch(color="#c00000", label="Up-regulated"),
        Patch(color="#002060", label="Down-regulated"),
        
        Line2D([], [], color="gold", marker="*", label="FDR<0.001", ms=5, ls="None"),
        Line2D([], [], color="g", marker="+", label="FDR<0.01", ms=5, ls="None"),
        Line2D([], [], color="purple", marker="x", label="FDR<0.05", ms=5, ls="None"),

  
        SplitPatch(label="Significant metabolites", color1='#c00000', color2='#002060', edgecolor="black"),

        Line2D([], [], color="black", marker="s", markerfacecolor='None', label="Non-significant metabolites", ms=6.5, ls="None"),
    ],
    bbox_to_anchor=(0.18, 0.63),
    loc="center",
    ncols=1,
    fontsize=8,
    handler_map={SplitPatch: SplitPatchHandler()},
)
circos.ax.add_artist(scatter_legend)


