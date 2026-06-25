import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
"""Quick verification: ortholog conversion across all 12 species."""
import sys, os, re, logging
_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'scripts'))
import scanpy as sc
from rna.ortholog import convert_species_gene_names, detect_gene_id_type
logging.basicConfig(level=logging.WARNING)

ens_pat = re.compile(r'^ENS[A-Z]{0,4}G\d{11}$')
species_list = ['Cow','Ferret','Lizard','Marmoset','Opossum','Peromyscus','Pig','Rhabdomys','Sheep','Squirrel','Tree_shrew','Zebrafish']

print(f'{"Species":>12s}  {"Type":>16s}  {"Genes":>7s}  {"EnsBefore":>9s}  {"EnsAfter":>8s}  {"Unmapped":>9s}  {"Mapped%":>7s}')
print('-' * 85)

for sp in species_list:
    path = os.path.join(_REPO, 'projects', 'rna', '<GSE_ID>', sp, 'results', 'h5ad', '00_raw.h5ad')
    try:
        adata = sc.read_h5ad(path)
        n_genes = adata.n_vars
        gene_type = detect_gene_id_type(adata.var_names)
        n_ens_before = sum(1 for g in adata.var_names if ens_pat.match(str(g)))
        convert_species_gene_names(adata, sp.lower(), cache_dir='data/orthologs')
        n_ens_after = sum(1 for g in adata.var_names if ens_pat.match(str(g)))
        n_unmapped = sum(1 for g in adata.var_names if str(g).startswith('UNMAPPED_'))
        pct = (n_genes - n_unmapped) / n_genes * 100
        print(f'{sp:>12s}  {gene_type:>16s}  {n_genes:>7d}  {n_ens_before:>9d}  {n_ens_after:>8d}  {n_unmapped:>9d}  {pct:6.1f}%')
        del adata
    except Exception as e:
        print(f'{sp:>12s}  ERROR: {e}')
