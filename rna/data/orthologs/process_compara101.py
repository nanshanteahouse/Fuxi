"""Process Compara v101 TSV to extract 1:1 orthologs for 12 species → human.

Streams through the gzipped file line-by-line (never loads into memory).
Maps species Ensembl gene IDs → human Ensembl gene IDs.
A second step maps human gene IDs → gene symbols via BioMart.
"""
import gzip, json, sys, os, time, requests
from pathlib import Path
from collections import defaultdict

# Locate the Compara file relative to the configured data root.
# The file ships inside GSE246169 (multi-omics dataset).  If data_root()
# is not configured, you need to either set FUXI_DATA_ROOT or edit
# COMPARA_PATH below.
try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
    from core.utils import data_root
    COMPARA_PATH = os.path.join(data_root(), "GSE246169", "Compara.101.protein_default.homologies.tsv.gz")
except RuntimeError:
    # data_root() not configured — edit this path for your machine
    COMPARA_PATH = os.environ.get(
        "COMPARA_PATH",
        os.path.join(os.path.dirname(__file__), "Compara.101.protein_default.homologies.tsv.gz"),
    )
CACHE_DIR = Path(__file__).parent  # data/orthologs/

# Scientific names used in Ensembl Compara v101
OUR_SPECIES = {
    "sus_scrofa": "pig",
    "bos_taurus": "cow",
    "ovis_aries": "sheep",
    "mustela_putorius_furo": "ferret",
    "callithrix_jacchus": "marmoset",
    "monodelphis_domestica": "opossum",
    "danio_rerio": "zebrafish",
    "anolis_carolinensis": "lizard",
    "ictidomys_tridecemlineatus": "squirrel",
    "tupaia_belangeri": "tree_shrew",
    "peromyscus_maniculatus": "peromyscus",
    "rhabdomys_pumilio": "rhabdomys",
}

def main():
    t0 = time.time()
    counts = defaultdict(int)     # key → lines read
    mappings = defaultdict(dict)  # key → {species_gene_id: human_gene_id}
    n_total = 0

    print("Scanning Compara TSV for our 12 species...", flush=True)
    with gzip.open(COMPARA_PATH, 'rt') as f:
        f.readline()  # skip header
        for line in f:
            n_total += 1
            if n_total % 5000000 == 0:
                elapsed = time.time() - t0
                print(f"  {n_total//1000000}M rows ({elapsed:.0f}s), "
                      f"kept={sum(len(m) for m in mappings.values())}", flush=True)

            cols = line.split('\t')
            if len(cols) < 8:
                continue

            homology_type = cols[4].strip()
            if homology_type != 'ortholog_one2one':
                continue

            target_species = cols[7].strip()
            if target_species != 'homo_sapiens':
                continue

            source_species = cols[2].strip()
            key = OUR_SPECIES.get(source_species)
            if key is None:
                continue

            src_gene_id = cols[0].strip()
            human_gene_id = cols[5].strip()
            counts[key] += 1

            if src_gene_id and human_gene_id:
                mappings[key][src_gene_id] = human_gene_id

    elapsed = time.time() - t0
    print(f"\nProcessed {n_total:,} rows in {elapsed:.0f}s", flush=True)

    # Report
    for name in sorted(OUR_SPECIES.values()):
        n = len(mappings[name])
        print(f"  {name:>12s}: {n:>8d} orthologs found ({counts[name]:,} 1:1 human rows)")

    # Save intermediate: species_gene_id → human_gene_id
    for name, mapping in mappings.items():
        if mapping:
            cache_file = CACHE_DIR / f"{name}_to_human_orthologs_v101.json"
            cache_file.write_text(json.dumps(mapping, indent=2))
            print(f"  Saved {name}: {len(mapping)} entries → {cache_file}")

    # ── Step 2: human gene ID → gene symbol ──
    print("\nFetching human gene ID → gene symbol mapping from BioMart...")
    human_gene_ids = set()
    for mapping in mappings.values():
        human_gene_ids.update(mapping.values())
    print(f"  {len(human_gene_ids)} unique human gene IDs", flush=True)

    # Batch query: human gene ID → gene name via BioMart
    human_id_to_name = fetch_human_gene_names(human_gene_ids)
    print(f"  Got {len(human_id_to_name)} human gene names", flush=True)

    # ── Step 3: replace human gene ID with gene name ──
    for name, mapping in mappings.items():
        if not mapping:
            continue
        final = {}
        n_renamed = 0
        for src_id, human_id in mapping.items():
            gene_name = human_id_to_name.get(human_id, "")
            if gene_name:
                final[src_id] = gene_name
                n_renamed += 1
            else:
                final[src_id] = human_id  # keep ID if no name found

        cache_file = CACHE_DIR / f"{name}_to_human_orthologs.json"
        cache_file.write_text(json.dumps(final, indent=2))
        print(f"  Final {name}: {len(final)} entries ({n_renamed} with gene names) → {cache_file}")

    # Summary
    total_orth = sum(len(m) for m in mappings.values())
    print(f"\nDone. Total: {total_orth:,} orthologs across {sum(1 for m in mappings.values() if m)} species, {elapsed:.0f}s")


def fetch_human_gene_names(gene_ids, batch_size=5000):
    """Batch query BioMart for human gene ID → gene name.

    Returns {ENSG...: GENE_NAME} dict.
    """
    gene_list = sorted(gene_ids)
    result = {}

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        try:
            names = _query_biomart_batch(batch)
            result.update(names)
            print(f"    batch {i//batch_size+1}: {len(batch)} ids → {len(names)} names", flush=True)
        except Exception as e:
            print(f"    batch {i//batch_size+1}: ERROR {e}", flush=True)
        time.sleep(0.5)

    return result


def _query_biomart_batch(gene_ids):
    """Query BioMart for human gene names from gene IDs."""
    # We need the filter on ensembl_gene_id column. BioMart batch query via POST.
    gene_list_str = ",".join(gene_ids[:5000])  # BioMart has limits

    xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
       count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_gene_id" value="{gene_list_str}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>'''

    resp = requests.post(
        'https://www.ensembl.org/biomart/martservice',
        data={'query': xml},
        timeout=120,
    )
    resp.raise_for_status()

    result = {}
    for line in resp.text.strip().split('\n')[1:]:
        cols = line.split('\t')
        if len(cols) >= 2:
            result[cols[0].strip()] = cols[1].strip()
    return result


if __name__ == '__main__':
    main()
