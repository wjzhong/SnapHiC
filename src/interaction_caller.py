import os
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd
import skimage
import subprocess
from statsmodels.stats.multitest import multipletests
import gc
import sys
import h5py

def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key = lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #chrom_names = list(chrom_lens.keys())
    #chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms

def combine_chrom_interactions(directory):
    headers = "\t".join(["chr1", "x1", "x2","chr2","y1","y2","outlier_count",\
                         "case_avg","pvalue", "tstat", "fdr_dist"])
    output_filename_temp = os.path.join(directory, "combined_significances.bedpe.temp")
    output_filename = os.path.join(directory, "combined_significances.bedpe")
    input_filepattern = directory + '/significances.*.bedpe'
    proc = subprocess.Popen("awk 'FNR>1' " + input_filepattern + ' > ' + output_filename_temp, shell = True)
    proc.communicate()
    with open(output_filename, 'w') as ofile:
        ofile.write(headers + "\n")
    proc = subprocess.Popen(" ".join(["cat",output_filename_temp,">>",output_filename]), shell = True)
    proc.communicate()

def determine_dense_matrix_size(num_cells, dist, binsize, max_mem):
    max_mem_floats = max_mem * 1e9 
    max_mem_floats /= 8
    square_cells = max_mem_floats // num_cells
    mat_size = int(np.floor(np.sqrt(square_cells)) / 4)
    #print(mat_size)
    mat_size = max(int((dist // binsize) + 50), mat_size)
    #print('mat_size:', mat_size)
    #if mat_size < (dist // binsize):
    #    raise "Specified " + str(max_mem) + "GB is not enough for constructing dense matrix with distance " + str(dist) + "."
    return mat_size

def convert_sparse_dataframe_to_dense_matrix(d, mat_size, dist, binsize, upper_limit, num_cells, chrom_size, chrom_filename, tss_promoter_subset, enhancer_subset):
    d['i'] = (d.iloc[:,1] // binsize).astype(int)
    d['j'] = (d.iloc[:,4] // binsize).astype(int)
    #all_rows = set(range(d.shape[0]))
    max_distance_bin = dist // binsize
    chrom_bins = int(chrom_size // binsize)
    for i in range(0, chrom_bins + 1, int(mat_size - max_distance_bin)):
        #matrix_upper_bound = max(0, i - upper_limit)
        matrix_upper_bound = max(0, i)
        #matrix_lower_bound = min(i + mat_size + upper_limit, chrom_bins + 1)
        matrix_lower_bound = min(i + mat_size, chrom_bins + 1)
        keeprows = list(np.where((d['i'] >= matrix_upper_bound) & (d['j'] < matrix_lower_bound))[0])
        d_portion = d.iloc[keeprows, 0:6].reset_index(drop = True)
        d_portion.columns = ['chr1','x1','x2','chr2','y1','y2']
        #print(d_portion, 'd_portions shape')
        #skiprows = all_rows.difference(keeprows)
        hdf_file = h5py.File(chrom_filename + '.cells.hdf', 'r')
        portion = hdf_file[list(hdf_file.keys())[0]]
        #print(type(keeprows))
        #if isinstance(keeprows, list):
        #    print(len(keeprows))
        #    if len(keeprows) > 0:
        #        print(keeprows[0], type(keeprows[0]))
        if len(keeprows) == 0:
            continue
        portion = portion[keeprows, :]
        hdf_file.close()
        if portion.shape[0] == 0:
            continue
        portion = pd.DataFrame(portion)
        #print('portions shape', portion.shape)
        portion = pd.concat([d_portion, portion], axis = 1)
        #print('concatted portion', portion.shape)
        #print(np.where(portion.isnull().sum() > 0))
        
        # filter bin pairs based on the gene's TSS, promoter and enhancer
        d_portion_1 = tss_promoter_subset.merge(d_portion, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr1', 'x1', 'x2'])
        d_portion_2 = tss_promoter_subset.merge(d_portion, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr2', 'y1', 'y2'])
        d_portion_or = pd.concat([d_portion_1, d_portion_2]) 
        d_portion_or = d_portion_or.drop(['chr','z1','z2'], 1).drop_duplicates()    
        if len(d_portion_or) == 0:
            continue
        
        d_portion_AND = tss_promoter_subset.merge(d_portion_1.drop(['chr','z1','z2'], 1), how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr2', 'y1', 'y2'])
        d_portion_AND = d_portion_AND.drop(['chr','z1','z2'], 1).drop_duplicates()    
        
        d_portion_either = pd.concat([d_portion_or,d_portion_AND]).drop_duplicates(keep=False)
        
        d_portion_either_1 = tss_promoter_subset.merge(d_portion_either, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr1', 'x1', 'x2'])
        d_portion_either_1 = d_portion_either_1.drop(['chr','z1','z2'], 1)
        d_portion_either_2 = tss_promoter_subset.merge(d_portion_either, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr2', 'y1', 'y2'])
        d_portion_either_2 = d_portion_either_2.drop(['chr','z1','z2'], 1)
        
        if enhancer_subset:
            d_portion_XOR_1 = enhancer_subset.merge(d_portion_either_1, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr2', 'y1', 'y2'])
            d_portion_XOR_1 = d_portion_XOR_1.drop(['chr','z1','z2'], 1)
            d_portion_XOR_2 = enhancer_subset.merge(d_portion_either_2, how='inner', left_on=['chr', 'z1', 'z2'], right_on=['chr1', 'x1', 'x2'])
            d_portion_XOR_2 = d_portion_XOR_2.drop(['chr','z1','z2'], 1)
        else:
            d_portion_XOR_1 = d_portion_either_1
            d_portion_XOR_2 = d_portion_either_2
        
        d_portion_XOR = pd.concat([d_portion_XOR_1, d_portion_XOR_2])
        
        d_portion_AND_XOR = pd.concat([d_portion_AND, d_portion_XOR])     
        if len(d_portion_AND_XOR) == 0:
            continue
        
        portion = portion.merge(d_portion_AND_XOR, how = 'right', on = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2'])
        
        
        portion['i'] =  (portion.loc[:,'x1'] // binsize).astype(int)
        portion['j'] =  (portion.loc[:,'y1'] // binsize).astype(int)
        #portion_old = d[(d['i'] >= matrix_upper_bound) & \
        #            (d['j'] < matrix_lower_bound)]
        #if portion.shape[0] == 0:
        #    continue
        full_sparse = pd.DataFrame({'i': range(min(portion['i']), max(portion['j'])-1), \
                    'j': range(min(portion['i'])+1, max(portion['j']))})
        portion = portion.merge(full_sparse, on = ['i','j'], how = "outer")
        #print(portion.iloc[:,list(range(7)) + [11, 12]])
        #print("start", matrix_upper_bound, "end", matrix_lower_bound)
        dense_cells = []
        #print('here', portion.columns)
        #print(portion.head())
        #print(portion.dtypes[:20])
        #sys.stdout.flush()
        for cell_index in range(num_cells):
            cell_mat = sp.sparse.csr_matrix((portion.iloc[:, 6 + cell_index], \
                                            ((portion['i'] - matrix_upper_bound), \
                                                (portion['j'] - matrix_upper_bound))), \
                                            shape = (matrix_lower_bound - matrix_upper_bound, \
                                                    matrix_lower_bound - matrix_upper_bound))
            cell_mat = np.array(cell_mat.todense())
            cell_mat[np.tril_indices(cell_mat.shape[0],0)] = np.nan
            dense_cells.append(cell_mat)
        mat_3d = np.stack(dense_cells, axis = -1)
        if matrix_upper_bound == 0:
            #pad_size = abs(i - upper_limit)
            pad_size = abs(i)
            mat_3d = np.pad(mat_3d, ((pad_size, 0), (pad_size, 0), (0, 0)), mode = 'constant', constant_values = np.nan)
        if matrix_lower_bound == chrom_bins + 1:
            #pad_size = upper_limit #- 1
            pad_size = 0 
            mat_3d = np.pad(mat_3d, ((0, pad_size), (0, pad_size), (0, 0)), mode = 'constant', constant_values = np.nan)
        yield mat_3d, i
                        
def get_nth_diag_indices(mat, offset):
    rows, cols_orig = np.diag_indices_from(mat)
    cols = cols_orig.copy()
    if offset > 0:
        cols += offset
        rows = rows[:-offset]
        cols = cols[:-offset]
    return rows, cols

def get_neighbor_counts_matrix(shape, gap_large, gap_small, max_distance):
    a = np.zeros(shape)
    big_width = gap_large*2 + 1
    small_width = gap_small*2 + 1
    area = big_width**2 - small_width**2
    for i in range(0, big_width):
        val = np.sum(list(range(big_width - i)))
        rows, cols = get_nth_diag_indices(a, i + 1)
        a[rows, cols] -= val
        #rows, cols = get_nth_diag_indices(a, a.shape[0]-i)
        rows, cols = get_nth_diag_indices(a, max_distance-i)
        a[rows, cols] -= val
    for  i in range(0, small_width):
        val = np.sum(list(range(small_width - i)))
        rows, cols = get_nth_diag_indices(a, i + 1)
        a[rows, cols] += val
        #rows, cols = get_nth_diag_indices(a, a.shape[0]-i)
        rows, cols = get_nth_diag_indices(a, max_distance-i)
        a[rows, cols] += val
    a += area
    return a

def compute_significances(mat, upper_limit, lower_limit, num_cells, start_index, max_distance_bin):
    gc.collect()
    #sliding window
    #print('in function')
    #sys.stdout.flush()
    
    #remove edge cases that are used only as neighbors
    #mat = mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    #gc.collect()
    
    #apply one-sample t-test across all cells for each point
    zeroMatrix=np.zeros(mat.shape[0:2])
    ttest_results = stats.ttest_1samp(mat, zeroMatrix, axis = 2)
    pvals = ttest_results.pvalue
    tstat = ttest_results.statistic
    ##print(pvals.shape)
    mat = np.mean(mat, axis = -1)
    
    #keep only upper triangle
    mat = np.triu(mat, 1)
    pvals = np.triu(pvals, 1)
    pvals = np.nan_to_num(pvals, nan = 1)
    tstat = np.triu(tstat, 1)
    tstat = np.nan_to_num(tstat, nan = -1000)
    
    #convert matrix to dataframe
    mat = sp.sparse.coo_matrix(mat)
    pvals = sp.sparse.coo_matrix(pvals)
    tstat = sp.sparse.coo_matrix(tstat)
    result_mat = pd.DataFrame({'i': mat.row, 'j': mat.col, 'case_avg': mat.data})
    result_pval = pd.DataFrame({'i': pvals.row, 'j': pvals.col, 'pvalue': pvals.data})
    result_tstat = pd.DataFrame({'i':tstat.row, 'j': tstat.col, 'tstat': tstat.data})
    result = result_mat.merge(result_pval, on = ['i','j'], how = "outer")
    result = result.merge(result_tstat, on = ['i','j'], how = "outer")
    result.loc[:,'pvalue'] = result['pvalue'].fillna(0)
    result.loc[:,'tstat'] = result['tstat'].fillna(-1000)
    result.loc[:,'i'] += (start_index)# + upper_limit)
    result.loc[:,'j'] += (start_index)# + upper_limit)
    result.loc[:,'i'] = result['i'].astype(int)
    result.loc[:,'j'] = result['j'].astype(int)
    result = result[result['j'] - result['i'] <= max_distance_bin]
    return result

def call_interactions(indir, outdir, chrom_lens, binsize, dist, neighborhood_limit_lower = 3, \
                      neighborhood_limit_upper = 5, rank = 0, n_proc = 1, max_mem = 2, logger = None, tss_file = None, promoter_file = None, enhancer_file = None):
    logger.set_rank(rank)
    try:
        os.makedirs(outdir)
    except:
        pass
    
    
    def ceiling_division(n, d):
        return -(n // -d)
    
    def extend(dat):
        temp_list = []
        
        for x in dat.itertuples():
            tempZ1 = x.z1
            tempZ2 = x.z2
            for r in range(tempZ1, tempZ2, int(binsize)):
                temp_list.append([x.chr, r, r + int(binsize)])
        
        dat_extended = pd.DataFrame(temp_list, columns =['chr', 'z1', 'z2'])
        return dat_extended
    
    # load gene's TSS data
    tss = pd.read_csv(tss_file, sep = "\t", header = 0) 
    
    #we use 500 bp to define the promoter region around TSS
    tss_range = 500 
    tss.loc[:,'TSS_low'] = tss.loc[:,'TSS'] - tss_range 
    tss.loc[:,'TSS_up'] = tss.loc[:,'TSS'] + tss_range  
    tss.loc[:,'z1'] = ((tss.loc[:,'TSS_low'] // binsize ) * binsize).astype(int)      
    tss.loc[:,'z2'] = ( ceiling_division(tss.loc[:,'TSS_up'], binsize) * binsize ).astype(int)             
    tss_extended = extend(tss)
    tss_extended = tss_extended.drop_duplicates()
    
    # load promoter data
    if promoter_file:
        promoter = pd.read_csv(promoter_file, sep = "\t", header = None)
        promoter.columns = ['chr', 'start', 'end']     
        promoter.loc[:,'z1'] = ((promoter.loc[:,'start'] // binsize ) * binsize).astype(int)      
        promoter.loc[:,'z2'] = ( ceiling_division(promoter.loc[:,'end'], binsize) * binsize ).astype(int)    
        promoter_extended = extend(promoter)
        promoter_extended = promoter_extended.drop_duplicates()
    else:
        promoter_extended = tss_extended
    
    # load enhancer data
    if enhancer_file:
        enhancer = pd.read_csv(enhancer_file, sep = "\t", header = None)
        enhancer = enhancer.iloc[:, 0:3]
        enhancer.columns = ['chr', 'start', 'end']   
        enhancer.loc[:,'z1'] = ((enhancer.loc[:,'start'] // binsize ) * binsize).astype(int)      
        enhancer.loc[:,'z2'] = ( ceiling_division(enhancer.loc[:,'end'], binsize) * binsize ).astype(int)    
        enhancer_extended = extend(enhancer)
        enhancer_extended = enhancer_extended.drop_duplicates()
    else:
        enhancer_extended = None
    
    tss_promoter = tss_extended.merge(promoter_extended, how='inner', on=['chr', 'z1', 'z2']).drop_duplicates()
    
    
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    #print(rank, proc_chroms)
    #sys.stdout.flush()
    for chrom in proc_chroms:
        logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
        #print(rank, chrom)
        #d = pd.read_csv(chrom_filename, sep = "\t", header = None, usecols = [0,1,2,3,4,5, num_cells + 6])
        ##command = "awk -F '\t' '{print NF; exit}' " + chrom_filename
        ##proc_output = subprocess.check_output(command, shell = True, executable = "/bin/bash")
        ##num_cells = int(proc_output) - 7
        chrom_filename = os.path.join(indir, ".".join([chrom, "normalized", "combined", "bedpe"]))
        with h5py.File(chrom_filename + ".cells.hdf", 'r') as ifile:
            num_cells = ifile[chrom].shape[1]
        logger.write(f'\tprocessor {rank}: detected {num_cells} cells for chromosome {chrom}', \
                             append_time = False, allow_all_ranks = True, verbose_level = 2)
        #print('num_cells', num_cells)
        #sys.stdout.flush()
        d = pd.read_csv(chrom_filename, sep = "\t", header = None)
        #num_cells = d.shape[1] - 7  
        
        matrix_max_size = determine_dense_matrix_size(num_cells, dist, binsize, max_mem)
        #print(rank, matrix_max_size)
        
        tss_promoter_subset = tss_promoter[tss_promoter['chr'] == chrom]
        
        if enhancer_extended:
            enhancer_subset = enhancer_extended[enhancer_extended['chr'] == chrom]
        else:
            enhancer_subset = None
        
        submatrices = convert_sparse_dataframe_to_dense_matrix(d, matrix_max_size, \
                                                               dist, binsize, neighborhood_limit_upper, \
                                                               num_cells, chrom_lens[chrom], chrom_filename, tss_promoter_subset, enhancer_subset)
        max_distance_bin = dist // binsize
        results = []
        #print(matrix_max_size, neighborhood_limit_upper, neighborhood_limit_lower)
        #neighbor_counts_matrix = get_neighbor_counts_matrix((matrix_max_size + neighborhood_limit_upper * 2, \
        #                                                     matrix_max_size + neighborhood_limit_upper * 2), \
        #                                                 neighborhood_limit_upper, \
        #                                                 neighborhood_limit_lower, max_distance_bin)
        #print('num zeros_2d', len(np.where(neighbor_counts_matrix==0)[0]))
        #print('going in for')
        #sys.stdout.flush()
        for i, (submatrix, start_index) in enumerate(submatrices):
            logger.write(f'\tprocessor {rank}: computing background for batch {i} of {chrom}, start index = {start_index}', \
                              verbose_level = 3, allow_all_ranks = True, append_time = False)
            #print('iteration', i)
            #print('start_index', start_index)
            #sys.stdout.flush()
            if i > 0:
                limit =  i * (matrix_max_size - max_distance_bin) #- neighborhood_limit_upper
                #results[-1] = results[-1][results[-1]['i'] < limit]
                results[-1] = results[-1][results[-1]['i'] < start_index]
            #start_index = i * (matrix_max_size - max_distance_bin) - neighborhood_limit_upper
            #print(start_index)
            submat_result = compute_significances(submatrix, neighborhood_limit_upper, \
                                                  neighborhood_limit_lower, num_cells, start_index, \
                                                  max_distance_bin)
            #print('returned')
            results.append(submat_result)
        #print(rank, 'offtheloop')
        #print(rank, len(results))
        results = pd.concat(results, axis = 0)
        #print(rank, chrom, results.shape[0])
        min_index = 0
        max_index = results['j'].max()
        #print(max_index, min_index, results['i'].dtype, results['j'].dtype, neighborhood_limit_upper)
        results = results[(results['i'] >= min_index + neighborhood_limit_upper) & \
                          (results['j'] <= max_index - neighborhood_limit_upper)]
        #print(results.shape[0])
        
        
        def compute_fdr_by_dist(d):
            fdrs = multipletests(list(d['pvalue']), method = 'fdr_bh')[1]
            d.loc[:,'fdr_dist'] = fdrs
            return d
        
        results.loc[:,'i'] = (results['i'] * binsize).astype(int)
        results.loc[:,'j'] = (results['j'] * binsize).astype(int)
        
        #print('finishing', d.shape) 
        d = d.iloc[:, list(range(7))]
        d.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'outlier_count']
        #print(d.head())
        #print(results.head())
        #d = d.merge(results, left_on = ['x1', 'y1'], right_on = ['i', 'j'], how = "outer")
        d = d.merge(results, left_on = ['x1', 'y1'], right_on = ['i', 'j'])
        
        #compute FDR
        d.reset_index(drop = True, inplace = True)  
        d = d.groupby(d['j'] - d['i'], as_index = False).apply(compute_fdr_by_dist)
        
        
        #print(d.shape)
        d.drop(['i', 'j'], axis =1, inplace = True)
        logger.write(f'\tprocessor {rank}: computation for {chrom} completed. writing to file.', \
                             append_time = False, allow_all_ranks = True, verbose_level = 2)
        d.to_csv(os.path.join(outdir, ".".join(["significances", chrom, "bedpe"])), sep = "\t", index = False)   
                        
                        
