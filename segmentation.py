# -*- coding: utf-8 -*-

import glob
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from skimage.filters import meijering, rank, threshold_multiotsu
from skimage.morphology import disk
from skimage.restoration import denoise_tv_bregman
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
import pandas as pd
import seaborn as sns
#from scipy import stats
#from statsmodels.stats.multitest import multipletests
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-path", type=str, required=True)
parser.add_argument("-output_name", type=str, required=False, default="output")

args = parser.parse_args()

order = ['group 1', 'group 2', 'group 3', 'group 4', 'group 5', 'group 6', 'group 7', 'group 8', 'group 9', 'group 10']
disk3 = disk(3)

def ul(array):
    q1 = np.percentile(array, 25)
    q3 = np.percentile(array, 75)
    iqr = q3 - q1
    fl = q1 - (1.5*iqr)
    fh = q3 + (1.5*iqr)
    return fl, fh

"""            
def prepare(img_names):

    dct_imgs = {}
    dds = {}
    dct_denoised = {}

    for i in range(len(img_names)):
        print("\t\t", img_names[i])
        img10 = cv2.imread(img_names[i], cv2.IMREAD_GRAYSCALE)
        dct_imgs[img_names[i]] = img10

        img10d = denoise_tv_bregman(img10, weight=10)

        meid = meijering(img10d, sigmas=np.arange(0,4,0.1))

        t = threshold_multiotsu(meid, classes=4)
        tmprr = (meid < t[0])# == False
        dct_denoised[img_names[i]] = tmprr# == False

        dd = ndi.distance_transform_edt(tmprr)
        dds[img_names[i]] = dd

    return dds, dct_denoised, dct_imgs
"""

def prepare(img_names):

    dct_imgs = {}
    dds = {}
    dct_denoised = {}

    for i in range(len(img_names)):
        print("\t\t", img_names[i])
        img10 = cv2.imread(img_names[i], cv2.IMREAD_GRAYSCALE)
        dct_imgs[img_names[i]] = img10

        img10d = denoise_tv_bregman(img10, weight=10)

        meid = meijering(img10d, sigmas=np.arange(0,4,0.1))
        t = threshold_multiotsu(meid, classes=4)
        tmprr = meid < t[0]

        im = np.asarray((1-img10d)*255, dtype=np.uint8)
        im = np.log(im)
        im[np.isinf(im)] = 0
        im[np.isnan(im)] = 0
        im = im / np.max(im)
        im = 1 - im
        im = rank.mean(im, footprint=disk3)
        
        tt = threshold_multiotsu(im, classes=4)
        tmp = im > tt[0]

        im = tmprr * tmp
        
        im = np.asarray(im, dtype=np.uint8)
        im = cv2.morphologyEx(im, cv2.MORPH_OPEN, disk3, iterations=1)

        dct_denoised[img_names[i]] = im

        dd = ndi.distance_transform_edt(im) # before tmprr
        dds[img_names[i]] = dd

    return dds, dct_denoised, dct_imgs


def get_all_distances_above_zero(dct):
    outall = {}
    for k,v in dct.items():
        outall[k] = []
        for kk,vv in v.items():
            outall[k].append(vv[vv > 0])

    return outall#pd.DataFrame(dfs)



def remover_template(ccm, cm, dc, thrs):
    remove = {}
    for k,v in cm.items():
        group = k.split("/")[-2]
        print(k, group)
        areas = dist_clusters[group][k][1]
        conn_comps = dist_clusters[group][k][0]
        zz = np.zeros_like(areas)
        for i in np.arange(1,conn_comps,1)[np.log(v) > thrs["median"][1]]:
            zz[ areas == i ] = 1
            
        remove[k] = zz
        
    return remove

def segment(dist, denoised, rem):
    
    dct_instance_seg = {}
    
    for k,v in dist.items():
        coords10 = peak_local_max(v, labels=denoised[k], footprint=np.ones((15, 15)), min_distance=15)#,
        mask = np.zeros(v.shape, dtype=bool)
        mask[tuple(coords10.T)] = True
        markers, _ = ndi.label(mask)
    
        labels10 = watershed( -v, markers, mask = denoised[k], watershed_line=True)
        
        rema = np.unique(labels10[rem[k] == 1])
        newlab = labels10.copy()
        for xx in rema:
            newlab[newlab == xx] = 0
        dct_instance_seg[k] = newlab
    
    return dct_instance_seg
    

def plotter(dct_imgs, dct_denoised, dct_instance_seg, filename):
    figs = []
    print("Generating images and compiling them into a PDF")
    #text_kwargs = dict(color='r', fontsize=6)
    for k in order:#dct_imgs.items():
        for kk,vv in dct_imgs[k].items():#v.items():
    
            fig, axs = plt.subplots(1,3) #8
            fig.set_size_inches(50, 50) #50,50
            na = " - ".join(kk.split("/")[-2:]) + " \n "
            axs[0].imshow(dct_imgs[k][kk], cmap = "gray")
            axs[0].set_title(na + "Original gray scale")
            axs[1].imshow(dct_denoised[k][kk], cmap = "gray")
            axs[1].set_title(na + "The clean cells")
            axs[2].imshow(dct_instance_seg[k][kk], cmap="gist_earth")
            u = np.unique(dct_instance_seg[k][kk])
            number = len(u)
            axs[2].set_title(na + "Cleaned and segmented cells ({})".format(number))
    
            axs[0].get_xaxis().set_visible(False)
            axs[0].get_yaxis().set_visible(False)
            axs[1].get_xaxis().set_visible(False)
            axs[1].get_yaxis().set_visible(False)
            axs[2].get_xaxis().set_visible(False)
            axs[2].get_yaxis().set_visible(False)
    
            figs.append(fig)
            plt.close()
            plt.pause(0.001)
    
    #filename = "segmentation_results_multi_mei1.pdf"

    pp = PdfPages(filename)
    for fig in figs:
        fig.savefig(pp, format='pdf', bbox_inches='tight')
    pp.close()

def info_gater(dct):
    out = {}
    dfs = {"Groups" : [], "Image": [], "Cell_size": []}
    for k,v in dct.items():
        out[k] = []
        tmp = []
        for kk,vv in v.items():
            u,c = np.unique(vv, return_counts=True)
            out[k].append(c[1:])
            tmp += [kk.split("/")[-1]] * (len(c) - 1)
            
        clsz = np.concatenate(out[k])
        out[k] = clsz
        dfs["Image"] += tmp
        dfs["Cell_size"] += list(clsz)
        dfs["Groups"] += len(tmp) * [k]

    dfs["Cell_size_(log)"] = np.log(dfs["Cell_size"])
    return out, pd.DataFrame(dfs)


def get_cell_count_area(dct, name):
    out_dct = dict(img_names = [], cell_area = [], img_area = [], percentage = [])
    for k,v in dct.items():
        for kk,vv in v.items():
            n = "_".join(kk.split("/")[-2:])
            u,c = np.unique(vv, return_counts=True)
            sc = np.sum(c[1:])
            tot = np.prod(vv.shape)
            out_dct["img_names"].append(n)
            out_dct["cell_area"].append(sc)
            out_dct["img_area"].append(tot)
            out_dct["percentage"].append(sc/tot)
            
    pd.DataFrame(out_dct).to_excel("area_data_"+name+".xlsx")


if __name__ == "__main__":

    main_path = args.path#"drive-download-20230901T081504Z-001/"
    groups = glob.glob1(main_path, "*")
    dct_paths = {}
    print("Finding images using the path '{}'".format(args.path))
    for i in groups:
        dct_paths[i] = []
        mg = main_path + i + "/"
        tmp_img = glob.glob1(mg, "*")
    
        dct_paths[i] = [mg + ii for ii in tmp_img]
        
        
    #dct = {}
    #for i in ['group1', 'group2', 'group3']:
    #    dct[i] = dct_paths[i]
    #dct_paths = dct
    
    
    print("Preparing images for segmentation...")
    dist = {}
    dct_imgs = {}
    dct_denoised = {}
    for k,v in dct_paths.items():
        print("\t",k)
        dist[k] = {}
        dds, denoised, imgs = prepare(v)
        dist[k] = dds
        dct_imgs[k] = imgs
        dct_denoised[k] = denoised
    
    
    print("Collecting image data...")
    oa = get_all_distances_above_zero(dist)
    
    stuff = []
    for k,v in oa.items():
        stuff += v
    
    stuff = np.concatenate(stuff)
    
    print("Estimating 'true background' to find noisy patterns...")
    dist_clusters = {}
    prcnt = np.percentile(stuff, 75)
    for k,v in dist.items():
        dist_clusters[k] = {}
        for kk,vv in v.items():
            tr = vv > prcnt #provide a good true white space
            cr = cv2.connectedComponentsWithStats(np.asarray(tr, dtype=np.uint8))
            dist_clusters[k][kk] = cr
    
    
    print("using all images to remove noise and 'non-cell' areas...")
    size_dist_med = {}
    sizes_dct = {}
    for k,v in dist.items():
        for kk,vv in v.items():
                
            dist_score = vv.copy()
            tmp_size = dist_clusters[k][kk][-2][:,-1][1:].copy()
            con_comp = dist_clusters[k][kk][1].copy()
            mean_dist_score = []
            mean_dist_score_sqr = []
            med_dist_score = []
            med_dist_score_sqr = []
            for i in range(1, len(tmp_size)+1,1):
                concompdist = dist_score[ con_comp == i].copy()
                #concompdist_sqr = concompdist.copy()**2
                med_dist_score.append(np.median(concompdist))
            size_dist_med[kk] = med_dist_score
            sizes_dct[kk] = tmp_size
    
    
    print("combining cell sizes with median distance to zero")
    combined_med = {}
    for k, v in sizes_dct.items():
        combined_med[k] = np.array(size_dist_med[k]) * v
        
    ccombined_med = np.concatenate(list(combined_med.values()))
    
    
    print("finding outliers using upper and lower fences...")
    dct_thrs = {}
    tmp = np.log(ccombined_med)
    fl, fh = ul(tmp)
    dct_thrs["median"] = [fl, fh, np.mean(tmp), np.std(tmp)]
    
    
    rem = remover_template(ccombined_med, combined_med, dist_clusters, dct_thrs)
    
    dct_instance_seg = {}
    for k,v in dist.items():
        dct_instance_seg[k] = segment(v, dct_denoised[k], rem)
    
    plotter(dct_imgs, dct_denoised, dct_instance_seg, "watershed_segmentaion_{}.pdf".format(args.output_name))
    
    
    out, df = info_gater(dct_instance_seg)
    
    get_cell_count_area(dct_instance_seg, args.output_name)
    
    df.to_excel("Cell_size_sheet_{}.xlsx".format(args.output_name))
    
    sns.boxplot(data=df, x="Cell_size", y="Groups")
    plt.title("Cell size distributions")
    plt.savefig("boxplot_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()
    
    sns.boxplot(data=df, x="Cell_size_(log)", y="Groups")
    plt.title("Cell size distributions (log)")
    plt.savefig("boxplot_log_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()
    
    sns.histplot(data=df, x="Cell_size", y="Groups")
    plt.title("Cell size distributions")
    plt.savefig("histplot_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()
    
    sns.histplot(data=df, x="Cell_size_(log)", y="Groups")
    plt.title("Cell size distributions (log)")
    plt.savefig("histplot_log_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()
    
    sns.violinplot(data=df, x="Cell_size", y="Groups")
    plt.title("Cell size distributions")
    plt.savefig("violinplot_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()
    
    sns.violinplot(data=df, x="Cell_size_(log)", y="Groups")
    plt.title("Cell size distributions (log)")
    plt.savefig("violinplot_log_{}.pdf".format(args.output_name))
    plt.pause(0.001)
    plt.close()