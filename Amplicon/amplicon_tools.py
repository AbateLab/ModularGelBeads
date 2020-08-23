import itertools
import math
import os
import loompy
import h5py
import copy
import umap
import numpy as np
import pandas as pd


from collections import Counter
import seaborn as sns; sns.set(style="white", color_codes=True)
import matplotlib
import matplotlib.colors as mcol
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from IPython.display import display, HTML
from matplotlib.ticker import NullFormatter
from scipy.stats import binned_statistic
from IPython.display import display, HTML
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import scipy.cluster.hierarchy as sch

plt.rcParams['axes.unicode_minus'] = False

class Cluster_cells(object):
    """ this class contains methods to evaluate genotyped amplicons"""

    def __init__(self, allel_traits):
        self.barcodes = []
        self.data_tensor = []
        self.clusters = []
        self.linkage = []
        self.cells = 0
        self.genotypes = 0
        self.allel_traits = allel_traits
        self.m_cell_idt = []

    @classmethod
    def read_hd5(cls, genotype_file, allel_traits=(0, 1, 2, 3), merged=False):
        """r"""
        run = cls(allel_traits)
        run._select_allel_variants(genotype_file.T, allel_traits, merged)
        return run

    def _select_allel_variants(self, genotype, allel_traits, merged):
        """reads the hd5 genotyping outpu file and converts the data in a binary tensor
         of shape m * n * k, with m genetic variants, n cells and k the considered allelic traits 
        (0: homozygote wt, 1: heterozygote, 2: homozygote alternate, 3: unknown/not applicable)"""
        
        self.hd5_data = genotype
        dat = np.array(genotype)
        self.genotypes, self.cells = dat.shape
        try:
            self.data_tensor = np.array([dat == i for i in allel_traits]).astype('int')
        except TypeError:
            self.data_tensor = np.array([dat == allel_traits]).astype('int')
        if merged:
            self.data_tensor = self.data_tensor.sum(axis=0).reshape(1,self.data_tensor.shape[1],self.data_tensor.shape[2])
        return

    def cell_identity(self, sparsity_thresh=0.05, dense_dot=False):
        """returns for each allele trait a n * n cell identity matrix"""
        for i, m in enumerate(self.data_tensor):
            counts = np.sum(m)
            if counts / (self.cells * self.genotypes) < sparsity_thresh:
                A = csr_matrix(m)
                rec = A.T.dot(A).toarray()
            else:
                if dense_dot:
                    print(('matrix {} is not sparse, try dense dot product').format(i))
                    rec = np.dot(m.T, m)
                else:
                    chunks = 300
                    cell_arr = np.arange(self.cells)
                    cell_chunk = [cell_arr[i:i + chunks] for i in range(0, self.cells, chunks)]
                    rec = np.zeros([self.cells, self.cells])
                    for i, k in enumerate(cell_chunk):
                        for j, l in enumerate(cell_chunk):
                            dat1 = np.dot(m[:, k].T, m[:, l])
                            c1 = len(k)
                            c2 = len(l)
                            rec[i * chunks:i * chunks + c1, j * chunks:j * chunks + c2] = dat1

            self.m_cell_idt.append(rec)
        self.m_cell_idt = np.array(self.m_cell_idt)
        return
    
    def cos_similarity(self):
        dat = self.data_tensor.sum(axis=0)
        norm_dat = np.linalg.norm(dat,axis=0)
        self.cos_sim = np.dot(dat.T, dat) / np.dot(norm_dat.reshape(-1,1), norm_dat.reshape(1,-1))

    def angular_similarity(self):
        try:
            arccos = np.arccos(self.cos_sim)
        except AttributeError:
            self.cos_similarity()
            arccos = np.arccos(self.cos_sim)
        arccos[np.isnan(arccos)] = 0
        self.ang_sim = 1 - 2*arccos/np.pi
    
    def jaccard_similarity(self):
        dat = self.data_tensor.sum(axis=0)
        sq_norm_dat = np.linalg.norm(dat,axis=0)**2
        self.jaccard_sim = np.dot(dat.T, dat)/(-np.dot(dat.T, dat)+sq_norm_dat.reshape(-1,1)+sq_norm_dat)
    
    def make_cluster(self, method, data=None, cmap=plt.cm.YlGnBu):
        """rr"""
        try:
            if data == None: pass
            dat = self.m_cell_idt.sum(axis=0)
        except ValueError:
            dat = data
        
        self.linkage = sch.linkage(dat, method=method)
        
        # Compute and plot first dendrogram.
        fig = plt.figure(figsize=(16,16))
        ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
        Z1 = sch.dendrogram(self.linkage, orientation='left')
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        # Compute and plot second dendrogram.
        ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
        Z2 = sch.dendrogram(self.linkage)
        ax2.set_xticks([])
        ax2.set_yticks([])
        
        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        dat = dat[idx1,:]
        dat = dat[:,idx2]
        im = axmatrix.matshow(dat, aspect='auto', origin='lower', cmap=cmap)
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        
        # Plot colorbar.
        axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
        plt.colorbar(im, cax=axcolor)
        plt.savefig('clusters.svg', dpi=600)
        plt.savefig('clusters.png', dpi=600)
        #fig.show()
        self.cell_sort_idx = idx1
        
        return
        
    def retrieve_cluster(self, number):
        """rr"""
        self.clusters = sch.fcluster(self.linkage, number, criterion='maxclust')
        return

def load_genotypes(genotypes_path):
    # load genotyping data from hdf5 compressed file 
    
    with h5py.File(genotypes_path, 'r') as f:
        
        # import hdf5 layers into arrays  
        cell_barcodes = copy.deepcopy([c.decode('utf8').split('.')[0] for c in f['CELL_BARCODES']])
        variants = copy.deepcopy([v.decode('utf8') for v in f['VARIANTS']])
        
        genotypes = pd.DataFrame(np.transpose(f['GT']), index=cell_barcodes, columns=variants).sort_index()
        genotypes.index.name = 'cell_barcode'
        
        quality = pd.DataFrame(np.transpose(f['GQ']), index=cell_barcodes, columns=variants).sort_index()
        quality.index.name = 'cell_barcode'
        
        total_depth = pd.DataFrame(np.transpose(f['DP']), index=cell_barcodes, columns=variants).sort_index()
        total_depth.index.name = 'cell_barcode'
        
        alt_depth = pd.DataFrame(np.transpose(f['AD']), index=cell_barcodes, columns=variants).sort_index()
        alt_depth.index.name = 'cell_barcode'
                      
        # calculate vaf - nan for division by 0
        #vaf = np.divide(alt_depth, total_depth)
        
    return genotypes, quality, total_depth, alt_depth#, vaf

def filter_variants(genotypes, alt_depth, total_depth, quality, min_alt_depth, min_total_depth, min_quality):
    # filters variants from genotyping data based on simple metrics

    genotypes[total_depth < min_total_depth] = 3
    genotypes[((genotypes == 1) | (genotypes == 2)) & (alt_depth < min_alt_depth)] = 3
    genotypes[quality < min_quality] = 3
    genotypes[genotypes.isnull()] = 3

    return genotypes

def load_variants(variants_file_path):
    # load variant annotations tsv file
    variant_info = pd.read_csv(variants_file_path, sep='\t', header=0, index_col=0, low_memory=False)
    variant_info.index.name = 'variant' 
    
    return variant_info


