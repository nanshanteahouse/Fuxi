#!/usr/bin/env python3
"""Step 01: Scrublet doublet detection (ATAC)"""
import sys,os,time,argparse
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..'))
from core.utils import setup_logger,resolve_config,safe_write,validate_adata
import snapatac2 as snap
def main():
 t0=time.time();p=argparse.ArgumentParser();p.add_argument("--config",required=True)
 a=p.parse_args();CFG=resolve_config(a.config)
 log=setup_logger("01_doublet",os.path.join(CFG.log_dir,"01_doublet.log"))
 log.info("Step 01: Scrublet doublet detection")
 data=snap.read(CFG.raw_h5ad,backed="r")
 log.info("Loaded: %d cells",data.n_obs)
 try:
  snap.pp.scrublet(data,random_state=CFG.execution.random_seed)
  n=data.obs['predicted_doublet'].sum()
  log.info("Doublets: %d / %d (%.1f%%)",n,data.n_obs,100*n/data.n_obs)
 except Exception as e:
  log.warning("Scrublet failed: %s",e);data.obs['predicted_doublet']=False;data.obs['doublet_scores']=0.0
 validate_adata(data,stage_name="01_doublet",logger=log)
 safe_write(data,CFG.doublet_h5ad,cfg=CFG);log.info("Done %.1fs",time.time()-t0)
if __name__=='__main__':main()
