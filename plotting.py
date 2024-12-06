import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests


def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)

def mak_circles(datasets, names, tit, func = np.mean, qs = [10, 20, 30, 40, 50, 60, 70, 80, 90]):
    
    ranges = []
    means = []
    for data in datasets:
        ranges.append(np.percentile(data, q=qs))
        means.append(func(data))
    
    ans = [str(np.round(i,2)) for i in means]
    fontsz = 55 * np.min([ 1, np.min([8/len(list(an)) for an in ans]) ])

    #maxmedian = []
    setoff = []
    xskip = []
    for i in ranges:
        setoff.append(i[0])
        xskip.append(i[-1])
    
    #maxmedian = np.max(maxmedian)
    xskip = np.max(xskip)
    setoff = np.min(setoff)
    
    for j,i in enumerate(ranges):
        ranges[j] = (ranges[j] - setoff) / (xskip - setoff) + 1
        #ranges[j] = ranges[j] / (xskip - setoff) + 1

    lnr = len(ranges)
    half = len(qs) // 2

    plt.figure(figsize=(9*lnr,10))
    ax = plt.gca()
    ax.cla() 
    
    for indx, j in enumerate(np.argsort(means)):#ranges:
        
        for ii in ranges[j][half:]:
            ax.add_patch(plt.Circle((0 + 4.2*indx, 0), ii, alpha=0.2, color="g", fill=True))
            ax.add_patch(plt.Circle((0 + 4.2*indx, 0), ii, alpha=1, lw=1, edgecolor="g", fill=False))
        
        ax.add_patch(plt.Circle((0 + 4.2*indx, 0), ranges[j][0] , alpha=1, color="white", fill=True))

        for ii in ranges[j][:half]:
            ax.add_patch(plt.Circle((0 + 4.2*indx, 0), ii , alpha=0.2, color="white", fill=True))
        
        for ii in ranges[j][:half]:
            ax.add_patch(plt.Circle((0 + 4.2*indx, 0), ii, alpha=1, lw=1, edgecolor="g", fill=False))
        
        ax.annotate(ans[j], xy = (0 + 4.2*indx, 0), size=fontsz, ha='center', va='center')
        ax.annotate(names[j], xy = (0 + 4.2*indx, 2.2), size=45, ha='center', va='center')

    # change default range so that new circles will work
    ax.set_xlim((-2.1, 4.2 * lnr-2.1))
    ax.set_ylim((-2.1, 2.6))
    plt.axis('off')
    plt.title(tit,  fontsize = 60)
    plt.show()
    
    #mak_circles([np.log(st32v), np.log(st44v), np.log(st32s), np.log(st44s)], 
    #        ["st32v", "st44v", "st32s", "st44s"],
    #        "Cell size distributions (log)")


def get_data(dct, datasets, bs, imgs, condition = "morning", comparex = ['evening|0d_ct', 'morning|0d_ct']):
    empty = {i:[] for i in dct }
    #[st33sub.Image, st33vis.Image, df.Image, st44vis.Image]
    for j,i in enumerate(imgs):
        im = np.array( [v.split(" ")[0] for v in i[bs[j]]] )
        for k, v in dct.items():
            empty[k] += list( datasets[j][ np.isin(im, v) ] )


    dd = empty.keys()
    plott = [i for i in dd if condition in i]
    mak_circles([empty[i] for i in plott], 
                plott,
                "Cell size distributions")
    
    plt.pause(0.0001)
    
    b = plt.violinplot([np.log(empty[i]) for i in plott])
    b = plt.xticks(np.arange(len(plott))+1, plott, rotation=45)
    plt.title("Violin plot (log)")

    plt.pause(0.0001)

    me = [np.mean(np.log(empty[i])) for i in plott]
    stdd = [np.std(np.log(empty[i])) for i in plott]
    b = plt.bar(np.arange(len(plott))+1, me)
    plt.errorbar(np.arange(len(plott))+1, me, yerr=stdd, fmt="o", color="k")
    b = plt.xticks(np.arange(len(plott))+1, plott, rotation=45)
    plt.title("Bar plot (log)")

    plt.pause(0.0001)

    me = [np.mean(empty[i]) for i in plott]
    stdd = [np.std(empty[i]) for i in plott]
    b = plt.bar(np.arange(len(plott))+1, me)
    plt.errorbar(np.arange(len(plott))+1, me, yerr=stdd, fmt="o", color="k")
    b = plt.xticks(np.arange(len(plott))+1, plott, rotation=45)
    plt.title("Bar plot")

    plt.pause(0.0001)

    dd = empty.keys()
    b = plt.boxplot([np.log(empty[i]) for i in plott])
    b = plt.xticks(np.arange(len(plott))+1, plott, rotation=45)
    plt.title("boxplot (log)")


    empty2 = {i:{} for i in comparex}#dct }
    empty3 = {i:{} for i in comparex}#dct }
    for i in comparex:#dd:
        for ii in plott:
            empty2[i][ii] = stats.permutation_test((np.log(empty[i]), np.log(empty[ii])), statistic, n_resamples=2000, random_state = 1234).pvalue
            empty3[i][ii] = stats.mannwhitneyu(np.log(empty[i]), np.log(empty[ii])).pvalue



    empty21 = {i:{} for i in comparex}#dct }
    empty31 = {i:{} for i in comparex}#dct }
    for i in comparex:#dd:
        for ii in plott:
            empty21[i][ii] = stats.permutation_test((empty[i], empty[ii]), statistic, n_resamples=2000, random_state = 1234).pvalue
            empty31[i][ii] = stats.mannwhitneyu(empty[i], empty[ii]).pvalue

    return empty2, empty3, empty21, empty31


labs = ["st32s", "st32v", "st44s", "st44v"]

def get_groups(st32visK, dct, dctt, col):
    
    st32visK["groups"] = ""
    st32visK["time"] = ""
    st32visK["lables"] = ""
    im = np.array( [v.split(" ")[0] for v in st32visK.Image] )
    for k, v in dct.items():
        st32visK["time"][np.isin(im, v)] = k
        
    st32visK["groups"] = im
    
    st32visK["lables"] = [dctt[i] if i in dctt else "" for i in list(st32visK[col])]
            

    return st32visK[ st32visK.time != "" ]



if __name__ == "__main__":



    # reading output files from the image segmentation
    st32sub = pd.read_excel("results/Cell_size_sheet_st32sub_191124_fl_fh_final.xlsx")
    st32vis = pd.read_excel("results/Cell_size_sheet_st32vis_191124_fl_fh_final.xlsx")
    st44sub = pd.read_excel("results/Cell_size_sheet_st44sub_201124_fl_hl.xlsx")
    st44vis = pd.read_excel("results/Cell_size_sheet_st44vis_201124_fl_hl.xlsx")
    # reading table with group divisions
    gr = pd.read_excel("results/division.xlsx")

    # noramalizng to 10x magnification
    st32s = st32sub.Cell_size*4
    st32v = st32vis.Cell_size/4
    st44s = st44sub.Cell_size*4
    st44v = st44vis.Cell_size

    # makes boolean arrays where unphysiological segmented ares are removed
    b1 = (st32s >= 200) * (st32s <= 20000)
    b2 = (st32v >= 200) * (st32v <= 20000)
    b3 = (st44s >= 200) * (st44s <= 20000)
    b4 = (st44v >= 200) * (st44v <= 20000)

    # defining new cell size distributions and 
    st32s = st32s[ b1 ] * 0.042
    st32v = st32v[ b2 ] * 0.042
    st44s = st44s[ b3 ] * 0.042
    st44v = st44v[ b4 ] * 0.042

    st32sub = st32sub[b1][['Groups', 'Image', 'Cell_size']]#.to_excel("results/st32sub_10x_cleaned.xlsx", index=False)
    st32vis = st32vis[b2][['Groups', 'Image', 'Cell_size']]#.to_excel("results/st32vis_10x_cleaned.xlsx", index=False)
    st44sub = st44sub[b3][['Groups', 'Image', 'Cell_size']]#.to_excel("results/st44sub_10x_cleaned.xlsx", index=False)
    st44vis = st44vis[b4][['Groups', 'Image', 'Cell_size']]#.to_excel("results/st44vis_10x_cleaned.xlsx", index=False)

    st32sub["Cell_size"] = st32s.to_numpy()
    st32vis["Cell_size"] = st32v.to_numpy()
    st44sub["Cell_size"] = st44s.to_numpy()
    st44vis["Cell_size"] = st44v.to_numpy()
    
    #st32
    dct_to_label = { 'evening-0d_ct': "CT12",
                     'evening-2d': '48',
                     'evening-3d': '72',
                     'evening-4d': '96',
                     'evening-5d': '120',
                     'morning-0d_ct': 'CT0',
                     'morning-2d': '36',
                     'morning-3d': '60',
                     'morning-4d': '84',
                     'morning-5d': '108'}
    
    #st44
    st44sublabel = {'3': 'sedentary mice',
     '44': 'sedentary mice',
     '46': 'sedentary mice',
     '58': 'sedentary mice',
     '60': 'sedentary mice',
     '63': 'sedentary mice',
     '65': 'sedentary mice',
     '6': 'Low-runners',
     '19': 'Low-runners',
     '22': 'Low-runners',
     '45': 'Low-runners',
     '49': 'Low-runners',
     '56': 'Low-runners',
     '66': 'Low-runners',
     '7': 'High-runners',
     '8': 'High-runners',
     '18': 'High-runners',
     '24': 'High-runners',
     '31': 'High-runners',
     '33': 'High-runners',
     '15': '3-day Inactive Low-runners',
     '21': '3-day Inactive Low-runners',
     '23': '3-day Inactive Low-runners',
     '40': '3-day Inactive Low-runners',
     '42': '3-day Inactive Low-runners',
     '50': '3-day Inactive Low-runners',
     '61': '3-day Inactive Low-runners',
     '11': '3-day Inactive High-runners',
     '13': '3-day Inactive High-runners',
     '16': '3-day Inactive High-runners',
     '25': '3-day Inactive High-runners',
     '41': '3-day Inactive High-runners',
     '52': '3-day Inactive High-runners',
     '4': 'sedentary mice*',
     '26': 'sedentary mice*',
     '29': 'sedentary mice*',
     '30': 'sedentary mice*',
     '36': 'sedentary mice*',
     '59': 'sedentary mice*',
     '64': 'sedentary mice*',
     '0': '3-week Low-runners',
     '14': '3-week Low-runners',
     '17': '3-week Low-runners',
     '20': '3-week Low-runners',
     '27': '3-week Low-runners',
     '28': '3-week Low-runners',
     '54': '3-week Low-runners',
     '1': '3-week High-runners',
     '2': '3-week High-runners',
     '9': '3-week High-runners',
     '39': '3-week High-runners',
     '43': '3-week High-runners',
     '48': '3-week High-runners',
     
     '12': 'Group before running',
      '34': 'Group before running',
      '47': 'Group before running',
      '55': 'Group before running',
      '53': 'Group before running'} # <-- inserted to provide right label to 53
    
    
    
    st44vislabel = { '3': 'sedentary mice',
     '44': 'sedentary mice',
     '46': 'sedentary mice',
     '58': 'sedentary mice',
     '60': 'sedentary mice',
     '63': 'sedentary mice',
     '65': 'sedentary mice',
     '6': 'Low-runners',
     '19': 'Low-runners',
     '22': 'Low-runners',
     '45': 'Low-runners',
     '49': 'Low-runners',
     '56': 'Low-runners',
     '7': 'High-runners',
     '8': 'High-runners',
     '18': 'High-runners',
     '31': 'High-runners',
     '33': 'High-runners',
     '51': 'High-runners',
     '15': '3-day Inactive Low-runners',
     '21': '3-day Inactive Low-runners',
     '23': '3-day Inactive Low-runners',
     '40': '3-day Inactive Low-runners',
     '42': '3-day Inactive Low-runners',
     '50': '3-day Inactive Low-runners',
     '61': '3-day Inactive Low-runners',
     '11': '3-day Inactive High-runners',
     '13': '3-day Inactive High-runners',
     '16': '3-day Inactive High-runners',
     '41': '3-day Inactive High-runners',
     '52': '3-day Inactive High-runners',
     '4': 'sedentary mice*',
     '26': 'sedentary mice*',
     '29': 'sedentary mice*',
     '30': 'sedentary mice*',
     '36': 'sedentary mice*',
     '59': 'sedentary mice*',
     '64': 'sedentary mice*',
     '0': '3-week Low-runners',
     '14': '3-week Low-runners',
     '17': '3-week Low-runners',
     '20': '3-week Low-runners',
     '27': '3-week Low-runners',
     '28': '3-week Low-runners',
     '54': '3-week Low-runners',
     '1': '3-week High-runners',
     '2': '3-week High-runners',
     '9': '3-week High-runners',
     '39': '3-week High-runners',
     '43': '3-week High-runners',
     '48': '3-week High-runners',
     
     '12': 'Group before running',
      '34': 'Group before running',
      '47': 'Group before running',
      '55': 'Group before running'}
    
    
    dct = { i: [ str(ii) for ii in gr["ID "][ (gr.group == i.split("-")[1]) * (gr.time == i.split("-")[0]) ] ] for i in np.unique( [u[0]+"-"+u[1] for u in gr[["time", "group"]].to_numpy() ] ) }    
    dct['morning-5d'].append("v55")
    
    st32sub = get_groups(st32sub, dct, dct_to_label, "time")
    st32vis = get_groups(st32vis, dct, dct_to_label, "time")
    
    st44vis = get_groups(st44vis, dct, st44vislabel, "groups")
    st44sub = get_groups(st44sub, dct, st44sublabel, "groups")
    
    
    st32sub.to_excel("results/st32sub_10x_cleaned.xlsx", index=False)
    st32vis.to_excel("results/st32vis_10x_cleaned.xlsx", index=False)
    st44sub.to_excel("results/st44sub_10x_cleaned.xlsx", index=False)
    st44vis.to_excel("results/st44vis_10x_cleaned.xlsx", index=False)

