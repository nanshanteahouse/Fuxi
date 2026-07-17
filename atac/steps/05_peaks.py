#!/usr/bin/env python3
"""Step 05: Post-clustering peak calling"""
import sys,os,time,argparse
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..'))
from core.utils import setup_logger,resolve_config,safe_write,validate_adata
import snapatac2 as snap
def main():
 t0=time.time();p=argparse.ArgumentParser();p.add_argument("--config",required=True)
 a=p.parse_args();CFG=resolve_config(a.config)
 log=setup_logger("05_peaks",os.path.join(CFG.log_dir,"05_peaks.log"))
 log.info("Step 05: Post-clustering peak calling")
 data=snap.read(CFG.clustered_h5ad,backed="r")
 log.info("Loaded: %d cells",data.n_obs)
 try:snap.tl.macs3(data,groupby='leiden',qvalue=CFG.atac.peak_qval)
 except Exception as e:log.warning("Per-cluster MACS3 failed: %s",e);snap.tl.macs3(data,qvalue=CFG.atac.peak_qval)
 snap.tl.merge_peaks(data);pd=snap.pp.make_peak_matrix(data,backend='hdf5')
 log.info("Peak matrix: %d x %d",pd.n_obs,pd.n_vars)
 snap.metrics.frip(pd);log.info("FRiP: mean=%.3f median=%.3f",pd.obs['frip'].mean(),pd.obs['frip'].median())
 try:
  os.makedirs(os.path.join(CFG.figure_dir,"05_peaks"),exist_ok=True)
  import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
  fig,ax=plt.subplots(figsize=(6,4));ax.hist(pd.obs['frip'],bins=50);ax.set_xlabel("FRiP")
  plt.savefig(os.path.join(CFG.figure_dir,"05_peaks","frip_distribution.png"),dpi=150);plt.close()
 except Exception as e:log.warning("FRiP histogram: %s",e)
 validate_adata(pd,stage_name="05_peaks",logger=log)
 safe_write(pd,CFG.peak_h5ad,cfg=CFG,compression_override=None);log.info("Done %.1fs",time.time()-t0)
if __name__=='__main__':main()
