import numpy as np
import pandas as pd
import logging
import os
import scipy.sparse as sparse


def load_ld(ld_prefix):
    # load SNPs info
    snps_filename_parquet = ld_prefix + '.parquet'
    snps_filename_gz = ld_prefix + '.gz'
    if os.path.exists(snps_filename_parquet):
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz):
        df_ld_snps = pd.read_table(snps_filename_gz, delim_whitespace=True)
    else:
        raise ValueError('couldn\'t find SNPs file %s or %s' % (snps_filename_parquet, snps_filename_gz))

    # load LD matrix
    R_filename = ld_prefix + '.npz'
    if not os.path.exists(R_filename):
        raise IOError('%s not found' % (R_filename))
    R = sparse.load_npz(R_filename).toarray()
    R = R + R.T
    assert np.allclose(np.diag(R), 1.0)
    logging.info('Loaded an LD matrix for %d SNPs from %s' % (R.shape[0], R_filename))

    # sanity checks
    assert R.shape[0] == R.shape[1]
    if R.shape[0] != df_ld_snps.shape[0]:
        raise ValueError('LD matrix has a different number of SNPs than the SNPs file')

    return R, df_ld_snps
