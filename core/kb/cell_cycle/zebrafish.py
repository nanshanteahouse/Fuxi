"""Cell-cycle gene lists for Danio rerio (loaded from zebrafish.yaml)."""

from core.kb.cell_cycle import _load_yaml_species

S_GENES, G2M_GENES = _load_yaml_species("zebrafish")
