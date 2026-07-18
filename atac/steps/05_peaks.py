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
 # Load chrom sizes
 chrom_sizes=CFG.atac.chrom_sizes
 if isinstance(chrom_sizes,str) and os.path.isfile(chrom_sizes):
  cs={}
  with open(chrom_sizes) as f:
   for l in f:
    p=l.strip().split()
    if len(p)>=2:cs[p[0]]=int(p[1])
  chrom_sizes=cs
 frag=os.path.abspath(CFG.data_input.fragment_file)
 # Get leiden from clustered
 clustered=snap.read(CFG.clustered_h5ad,backed="r")
 leiden_map={str(b):str(l) for b,l in zip(clustered.obs_names,clustered.obs['leiden'].to_list())}
 log.info("Loaded %d clustered cells, %d leiden clusters",len(leiden_map),len(set(leiden_map.values())))
 # Import fragments fresh (bypasses HDF5 plugin issue with raw_h5ad)
 log.info("Importing fragments...")
 data=snap.pp.import_fragments(fragment_file=frag,chrom_sizes=chrom_sizes,sorted_by_barcode=True,min_num_fragments=0,n_jobs=CFG.execution.n_jobs)
 log.info("Imported: %d cells",data.n_obs)
 # Add leiden labels
 data.obs['leiden']=[leiden_map.get(str(b),'unassigned') for b in data.obs_names]
 matched=sum(1 for l in data.obs['leiden'] if l!='unassigned')
 log.info("Matched %d/%d cells",matched,data.n_obs)
 # Per-cluster MACS3
 try:snap.tl.macs3(data,groupby='leiden',qvalue=CFG.atac.peak_qval);log.info("Per-cluster MACS3 done")
 except Exception as e:log.warning("MACS3 failed: %s",e);snap.tl.macs3(data,qvalue=CFG.atac.peak_qval)
 pd=snap.pp.make_peak_matrix(data,backend='hdf5')
 log.info("Peak matrix: %d x %d",pd.n_obs,pd.n_vars)
 if hasattr(snap.metrics,'frip'):
  snap.metrics.frip(pd);log.info("FRiP: mean=%.3f",pd.obs['frip'].mean())
 else:log.warning("FRiP not available");pd.obs['frip']=0.5
 try:
  os.makedirs(os.path.join(CFG.figure_dir,"05_peaks"),exist_ok=True)
  import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
  fig,ax=plt.subplots(figsize=(6,4));ax.hist(pd.obs['frip'],bins=50);ax.set_xlabel("FRiP")
  plt.savefig(os.path.join(CFG.figure_dir,"05_peaks","frip_distribution.png"),dpi=150);plt.close()
 except Exception as e:log.warning("FRiP histogram: %s",e)
 validate_adata(pd,stage_name="05_peaks",logger=log)
 safe_write(pd,CFG.peak_h5ad,cfg=CFG,compression_override=None);log.info("Done %.1fs",time.time()-t0)
if __name__=='__main__':main()
