import numpy as np
import pandas as pd
import logging
import os
import scipy.sparse as sparse
from io import BytesIO
import requests


def find_ld_prefix(chrom, min_pos):
    # min_pos = 40590585; chrom=12
    alkes_url = "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
    bp_starts = list(range(1, 252000001, 1000000))
    bp_ends = [x+3000000 for x in bp_starts]
    i = max([i for i,x in enumerate(bp_starts) if x<=min_pos])
    prefix = "chr%d_%d_%d" % (chrom, bp_starts[i], bp_ends[i])
    ld_prefix = os.path.join(alkes_url, prefix)
    return ld_prefix

def load_ld(ld_prefix, server=True, npz_suffix=""):
    # ld_prefix="https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr12_41000001_44000001"
    # load SNPs info
    snps_filename_parquet = ld_prefix + '.parquet'
    snps_filename_gz = ld_prefix + '.gz'

    if os.path.exists(snps_filename_parquet):
        print(snps_filename_parquet)
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz) | server==True:
        print(snps_filename_gz)
        df_ld_snps = pd.read_csv(snps_filename_gz, delim_whitespace=True)
    else:
        raise ValueError('couldn\'t find SNPs file %s or %s' % (snps_filename_parquet, snps_filename_gz))

    # load LD matrix
    R_filename = ld_prefix + '.npz' + str(npz_suffix)
    print(R_filename)
    if server:
        print("Reading .npz file from alkesgroup server")
        r = requests.get(R_filename, stream=True)
        R = sparse.load_npz(BytesIO(r.raw.read())).toarray()
        # F = HttpFile(R_filename)
        # with np.load(F) as data:
        #     k = data.keys()
        #     i = data.items()
        #     v = data.values()
    else:
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

#
# class HttpFile(object):
#     # https://stackoverflow.com/questions/48355140/load-npz-from-a-http-link
#     import urllib.request as urllib2
#     def __init__(self, url):
#         self.url = url
#         self.offset = 0
#         self._size = -1
#
#     def size(self):
#         if self._size < 0:
#             f = urllib2.urlopen(self.url)
#             self._size = int(f.headers["Content-length"])
#         return self._size
#
#     def read(self, count=-1):
#         req = urllib2.Request(self.url)
#         if count < 0:
#             end = self.size() - 1
#         else:
#             end = self.offset + count - 1
#         req.headers['Range'] = "bytes=%s-%s" % (self.offset, end)
#         f = urllib2.urlopen(req)
#         data = f.read()
#         # FIXME: should check that we got the range expected, etc.
#         chunk = len(data)
#         if count >= 0:
#             assert chunk == count
#         self.offset += chunk
#         return data
#
#     def seek(self, offset, whence=0):
#         if whence == 0:
#             self.offset = offset
#         elif whence == 1:
#             self.offset += offset
#         elif whence == 2:
#             self.offset = self.size() + offset
#         else:
#             raise Exception("Invalid whence")
#
#     def tell(self):
#         return self.offset
#
#
#
