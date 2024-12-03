# -*- coding: utf-8 -*-

import zipfile
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from skimage.filters import meijering, rank, threshold_multiotsu, gaussian
from skimage.morphology import disk
from skimage.restoration import denoise_tv_bregman
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-path", type=str, required=True)
parser.add_argument("-output_name", type=str, required=False, default="new_set")

args = parser.parse_args()


#data_path = 'drive-download-20240219T084952Z-001.zip'
dct_imgs = {}
dds = {}
dct_denoised = {}
#kernel = np.ones((5,5),np.uint8)
#kernel3 = np.ones((3,3),np.uint8)
disk5 = disk(5)

def get_zip_files(data_path):

    with zipfile.ZipFile(data_path, 'r') as zfile:
        files = zfile.namelist()

    files_dct = {}
    for i in sorted(files):
        name = i.split("/")[0]
        if name in files_dct:
            files_dct[name].append(i)
        else:
            files_dct[name] = [i]

    return files_dct



def driving_motor(data_path, files_dct):

    with zipfile.ZipFile(data_path, 'r') as zfile:
        # from: https://stackoverflow.com/questions/21357184/using-python-opencv-to-load-image-from-zip
        for k,v in files_dct.items():
            for i in v:
                data = zfile.read(i)
                img = cv2.imdecode(np.frombuffer(data, np.uint8), 1)
                img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
                yield img, i


def driving_motor_test(data_path, files_dct):

    with zipfile.ZipFile(data_path, 'r') as zfile:
        # from: https://stackoverflow.com/questions/21357184/using-python-opencv-to-load-image-from-zip
        #for k,v in files_dct.items():
        for i in files_dct["visFat"][:6]:
            data = zfile.read(i)
            img = cv2.imdecode(np.frombuffer(data, np.uint8), 1)
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            yield img, i


#xtra = {"subfat/sb3_0002_CY5.tiff":7,
#        "subfat/sb10_0002_CY5.tiff": 5        
#        }


    
def prepare(img, name):#, dct_imgs, dds, dct_denoised):

    dct_imgs[name] = 255 - img
    #imo = im.copy()
    im = denoise_tv_bregman(img)#, weight=10)
    mei = meijering(1-im)

    for i in range(3):
        try:
            #meid = gaussian(meid)
            t = threshold_multiotsu(mei, classes=4)
            mei[ mei > t[2] ] = t[2]
            #md = md / np.max(md)
        except:
            break


    mei = np.asarray(mei < t[0], dtype=np.uint8)
    mei = gaussian(mei) > 0
    
    
    im = np.asarray(im*255, dtype=np.uint8)
    im = np.log(im)
    im[np.isinf(im)] = 0#np.min(im)
    im[np.isnan(im)] = 0#np.min(im)
    im = im / np.max(im)
    im = 1 - im
    
    im = rank.mean(im, footprint=disk5)
    im = rank.mean(im, footprint=disk5)
    im = rank.mean(im, footprint=disk5)

    tt = threshold_multiotsu(im, classes=4)
    tmp = im > tt[1]
    im = mei * tmp
    
    
    im = np.asarray(im, dtype=np.uint8)
    
    #im = im * tmp
    #im = cv2.morphologyEx(im, cv2.MORPH_CLOSE, disk5, iterations=1)
    im = cv2.morphologyEx(im, cv2.MORPH_OPEN, disk5, iterations=1)
    dct_denoised[name] = im
    
    meid = ndi.distance_transform_edt(im)

    dds[name] = meid

    
    
def get_all_distances_above_zero(dct):
    outall = {}
    for k,v in dct.items():
        outall[k] = v[v > 0]

    return outall#pd.DataFrame(dfs)


def remover_template(ccm, cm, dc, thrs):
    remove = {}
    for k,v in cm.items():
        #group = k.split("/")[-2]
        print(k)
        areas = dist_clusters[k][1]
        conn_comps = dist_clusters[k][0]
        zz = np.zeros_like(areas)
        for i in np.arange(1,conn_comps,1)[np.log(v) > thrs["median"][1]]:
            zz[ areas == i ] = 1

        remove[k] = zz

    return remove


def ul(array):
    q1 = np.percentile(array, 25)
    q3 = np.percentile(array, 75)
    iqr = q3 - q1
    fl = q1 - (1.5*iqr)
    fh = q3 + (1.5*iqr)
    return fl, fh


def segment(dist, denoised, rem):

    dct_instance_seg = {}

    for k,v in dist.items():
        coords10 = peak_local_max(v, labels=denoised[k], footprint=np.ones((50, 50)), min_distance=15)#,
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


def info_gater(dct, files_dct):
    out = {}
    dfs = {"Groups" : [], "Image": [], "Cell_size": []}
    for k,v in files_dct.items():
        out[k] = []
        tmp = []
        for kk in v:
            u,c = np.unique(dct[kk], return_counts=True)
            out[k].append(c[1:])
            tmp += [kk.split("/")[-1]] * (len(c) - 1)

        clsz = np.concatenate(out[k])
        out[k] = clsz
        dfs["Image"] += tmp
        dfs["Cell_size"] += list(clsz)
        dfs["Groups"] += len(tmp) * [k]

    dfs["Cell_size_(log)"] = np.log(dfs["Cell_size"])
    return out, pd.DataFrame(dfs)


def get_cell_count_area(dct, files_dct, name):
    out_dct = dict(img_names = [], cell_area = [], img_area = [], percentage = [])
    for k,v in files_dct.items():
        for kk in v:
            n = "_".join(kk.split("/")[-2:])
            u,c = np.unique(dct[kk], return_counts=True)
            sc = np.sum(c[1:])
            tot = np.prod(dct[kk].shape)
            out_dct["img_names"].append(n)
            out_dct["cell_area"].append(sc)
            out_dct["img_area"].append(tot)
            out_dct["percentage"].append(sc/tot)

    pd.DataFrame(out_dct).to_excel("area_data_"+name+".xlsx")



def plotter(dct_imgs, dct_denoised, dct_instance_seg, filename):
    figs = []
    print("Generating images and compiling them into a PDF")
    
    for kk,vv in dct_imgs.items():#v.items():

        fig, axs = plt.subplots(1,3) #8
        fig.set_size_inches(50, 50) #50,50
        na = " - ".join(kk.split("/")[-2:]) + " \n "
        axs[0].imshow(dct_imgs[kk], cmap = "gray")
        axs[0].set_title(na + "Original gray scale")
        axs[1].imshow(dct_denoised[kk], cmap = "gray")
        axs[1].set_title(na + "The clean cells")
        axs[2].imshow(dct_instance_seg[kk], cmap="gist_earth")
        u = np.unique(dct_instance_seg[kk])
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




if __name__ == "__main__":

    print("Reading and preparing images...")    
    files_dct = get_zip_files(args.path)

    for img, name in driving_motor(args.path, files_dct):
        prepare(img, name)
    
    
    print("Collecting image data...")
    oa = get_all_distances_above_zero(dds)

    stuff = []
    for k,v in oa.items():
        stuff.append(v)

    stuff = np.concatenate(stuff)
    
    
    dist_clusters = {}
    prcnt = np.percentile(stuff, 75)
    for kk,vv in dds.items():
        #dist_clusters[k] = {}
        #for kk,vv in v.items():
        tr = vv > prcnt#provide a good true white space
        cr = cv2.connectedComponentsWithStats(np.asarray(tr, dtype=np.uint8))
        dist_clusters[kk] = cr
        
        
    print("using all images to remove noise and 'non-cell' areas...")
    size_dist_med = {}
    sizes_dct = {}
    for k,v in dds.items():
        #for kk,vv in v.items():

        dist_score = v.copy()
        tmp_size = dist_clusters[k][-2][:,-1][1:].copy()
        con_comp = dist_clusters[k][1].copy()
        mean_dist_score = []
        mean_dist_score_sqr = []
        med_dist_score = []
        med_dist_score_sqr = []
        for i in range(1, len(tmp_size)+1,1):
            concompdist = dist_score[ con_comp == i].copy()
            #concompdist_sqr = concompdist.copy()**2
            med_dist_score.append(np.median(concompdist))
        size_dist_med[k] = med_dist_score
        sizes_dct[k] = tmp_size


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

    print("Producing tables and figures...")
    rem = remover_template(ccombined_med, combined_med, dist_clusters, dct_thrs)

    dct_instance_seg = segment(dds, dct_denoised, rem)
    
    
    plotter(dct_imgs, dct_denoised, dct_instance_seg, args.output_name+".pdf")
    
    out, df = info_gater(dct_instance_seg, files_dct)
    
    
    get_cell_count_area(dct_instance_seg, files_dct, args.output_name)


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
    

    
