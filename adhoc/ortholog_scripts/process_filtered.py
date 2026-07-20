"""Process pre-filtered Compara v101 lines into ortholog mappings.

Reads from stdin (pre-filtered to only lines with both ortholog_one2one and homo_sapiens).
Writes {species_name}_to_human_orthologs.json to CACHE_DIR.
"""
import sys, json, time
from pathlib import Path
from collections import defaultdict

CACHE_DIR = Path(__file__).parent  # data/orthologs/

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
    counts = defaultdict(int)
    mappings = defaultdict(dict)  # key → {species_gene_id: human_gene_symbol}
    n_total = 0
    n_human_genes = 0

    # Read from stdin — already filtered to ortholog_one2one + homo_sapiens lines
    for line in sys.stdin:
        n_total += 1
        if n_total % 1000000 == 0:
            elapsed = time.time() - t0
            kept = sum(len(m) for m in mappings.values())
            print(f"  {n_total//1000000}M lines processed ({elapsed:.0f}s), kept={kept}", flush=True)

        cols = line.strip().split('\t')
        if len(cols) < 8:
            continue

        source_species = cols[2].strip()
        key = OUR_SPECIES.get(source_species)
        if key is None:
            continue

        src_gene_id = cols[0].strip()
        human_gene_id = cols[5].strip()
        counts[key] += 1

        # Build mapping: overwrites are fine (latest 1:1 takes priority)
        if src_gene_id and human_gene_id:
            mappings[key][src_gene_id] = human_gene_id
            n_human_genes += 1

    elapsed = time.time() - t0
    print(f"\nProcessed {n_total:,} filtered lines in {elapsed:.0f}s", flush=True)

    for name in sorted(OUR_SPECIES.values()):
        n = len(mappings[name])
        if n > 0 or counts[name] > 0:
            print(f"  {name:>12s}: {n:>8d} unique gene pairs ({counts[name]:,} total lines)", flush=True)

    # ── Step 2: convert human ENSG IDs → gene symbols via BioMart ──
    print("\nFetching human gene ID → gene symbol from BioMart...", flush=True)
    human_gene_ids = set()
    for mapping in mappings.values():
        human_gene_ids.update(mapping.values())
    print(f"  {len(human_gene_ids)} unique human gene IDs", flush=True)

    human_id_to_name = fetch_human_gene_names(human_gene_ids)
    print(f"  Got {len(human_id_to_name)} gene names", flush=True)

    # ── Step 3: write final files ──
    for name, mapping in mappings.items():
        if not mapping:
            continue
        final = {}
        n_with_name = 0
        for src_id, human_id in mapping.items():
            gene_name = human_id_to_name.get(human_id, human_id)
            final[src_id] = gene_name
            if gene_name != human_id:
                n_with_name += 1

        cache_file = CACHE_DIR / f"{name}_to_human_orthologs.json"
        cache_file.write_text(json.dumps(final, indent=2, sort_keys=True))
        print(f"  {name}: {len(final)} entries ({n_with_name} with gene names) → {cache_file.name}", flush=True)

    total_orth = sum(len(m) for m in mappings.values())
    print(f"\nDone. {total_orth:,} total orthologs across {sum(1 for m in mappings.values() if m)} species, {elapsed:.0f}s", flush=True)


def fetch_human_gene_names(gene_ids, batch_size=7000):
    """Batch query BioMart for human gene ID → gene name."""
    import requests
    gene_list = sorted(gene_ids)
    result = {}

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
       count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_gene_id" value="{",".join(batch)}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>'''
        try:
            resp = requests.post('https://www.ensembl.org/biomart/martservice',
                                data={'query': xml}, timeout=120)
            for line in resp.text.strip().split('\n')[1:]:
                cols = line.split('\t')
                if len(cols) >= 2 and cols[1].strip():
                    result[cols[0].strip()] = cols[1].strip()
            print(f"    batch {i//batch_size+1}: {len(batch)} ids → {len(result)-i} new names", flush=True)
        except Exception as e:
            print(f"    batch {i//batch_size+1}: ERROR {e}", flush=True)
        import time as _time
        _time.sleep(0.3)

    return result


if __name__ == '__main__':
    main()
