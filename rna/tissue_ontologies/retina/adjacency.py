"""
tissue_ontologies/retina/adjacency.py — Anatomically validated cell-type adjacency map
of the vertebrate retina.

This module defines known direct physical and synaptic connections between retinal
cell types, curated from electron-microscopy connectomics, serial-section
reconstruction, and electrophysiological tracing studies.

References
----------
- Masland, R.H. (2012). The Neuronal Organization of the Retina. *Neuron*, 76(2),
  266–280. DOI: 10.1016/j.neuron.2012.10.002
- Hoon, M., Okawa, H., Della Santina, L., & Wong, R.O.L. (2014). Functional
  architecture of the retina: Development and disease. *Progress in Retinal and
  Eye Research*, 42, 44–84. DOI: 10.1016/j.preteyeres.2014.06.003
- Helmstaedter, M. (2013). Cellular-resolution connectomics: challenges of dense
  neural circuit reconstruction. *Nature Methods*, 10, 501–507.
  DOI: 10.1038/nmeth.2476
- Demb, J.B. & Singer, J.H. (2015). Functional Circuitry of the Retina.
  *Annual Review of Vision Science*, 1, 263–289.
  DOI: 10.1146/annurev-vision-082114-035334

Notes
-----
Each adjacency is a (source, target, adjacency_type) tuple representing a
directed or structurally asymmetric relationship (e.g., synaptic transmission
from presynaptic source to postsynaptic target; glial ensheathment of a neuron).
Bidirectional or symmetric contacts are listed with the most physiologically
relevant direction.

Cell type names follow the canonical keys defined in
``rna/tissue_ontologies/retina/synonyms.py``.
"""

from typing import List, Tuple, Dict

ADJACENCY_TYPES: Dict[str, str] = {
    "synaptic": (
        "Chemical synapse — direct neurotransmitter release from presynaptic "
        "to postsynaptic cell"
    ),
    "gap_junction": (
        "Electrical synapse / gap junction — direct intercellular ion flow "
        "via connexin channels"
    ),
    "physical": (
        "Direct physical contact not mediated by a classical synapse "
        "(e.g., RPE-photoreceptor outer segment apposition)"
    ),
    "ensheathment": (
        "Glial wrapping or support — the source glial cell ensheaths or "
        "surrounds the target neuron (e.g., Mu\"ller glia endfeet)"
    ),
    "modulatory": (
        "Indirect modulation — ephaptic, volume-transmission, or feedback "
        "interaction (e.g., horizontal cell → bipolar cell feedback at the "
        "cone pedicle)"
    ),
}

ADJACENCY: List[Tuple[str, str, str]] = [
    # ═══════════════════════════════════════════════════════════════════════
    # Outer Plexiform Layer (OPL) — photoreceptor → second-order neurons
    # ═══════════════════════════════════════════════════════════════════════
    ("Rod_Photoreceptor", "Bipolar_Cell", "synaptic"),
    ("Cone_Photoreceptor", "Bipolar_Cell", "synaptic"),
    ("Rod_Photoreceptor", "Horizontal_Cell", "synaptic"),
    ("Cone_Photoreceptor", "Horizontal_Cell", "synaptic"),
    # Horizontal cell feedback onto photoreceptors and bipolar cell dendrites
    ("Horizontal_Cell", "Bipolar_Cell", "modulatory"),
    ("Horizontal_Cell", "Rod_Photoreceptor", "modulatory"),
    ("Horizontal_Cell", "Cone_Photoreceptor", "modulatory"),

    # ═══════════════════════════════════════════════════════════════════════
    # Inner Plexiform Layer (IPL) — bipolar cell → amacrine / RGC
    # ═══════════════════════════════════════════════════════════════════════
    ("Bipolar_Cell", "RGC", "synaptic"),
    ("Bipolar_Cell", "Amacrine_Cell", "synaptic"),
    ("Amacrine_Cell", "Bipolar_Cell", "synaptic"),          # feedback inhibition
    ("Amacrine_Cell", "RGC", "synaptic"),                   # feedforward inhibition
    ("Amacrine_Cell", "Amacrine_Cell", "synaptic"),         # lateral inhibition

    # ═══════════════════════════════════════════════════════════════════════
    # Rod Pathway — via AII amacrine cells
    # ═══════════════════════════════════════════════════════════════════════
    ("Bipolar_Cell", "Amacrine_Cell", "synaptic"),          # Rod BC → AII AC
    ("Amacrine_Cell", "Bipolar_Cell", "gap_junction"),      # AII AC → ON Cone BC
    ("Amacrine_Cell", "Bipolar_Cell", "synaptic"),          # AII AC → OFF Cone BC (glycinergic)

    # ═══════════════════════════════════════════════════════════════════════
    # Non-neuronal contacts — RPE, glia, vasculature, immune
    # ═══════════════════════════════════════════════════════════════════════
    # RPE-photoreceptor outer segment physical apposition
    ("RPE", "Rod_Photoreceptor", "physical"),
    ("RPE", "Cone_Photoreceptor", "physical"),

    # Mu"ller glia ensheathment of all major retinal neurons
    ("Muller_Glia", "Rod_Photoreceptor", "ensheathment"),
    ("Muller_Glia", "Cone_Photoreceptor", "ensheathment"),
    ("Muller_Glia", "Bipolar_Cell", "ensheathment"),
    ("Muller_Glia", "RGC", "ensheathment"),
    ("Muller_Glia", "Amacrine_Cell", "ensheathment"),
    ("Muller_Glia", "Horizontal_Cell", "ensheathment"),

    # Blood-retina barrier
    ("Vascular_Endothelial", "Pericyte", "physical"),

    # Microglia — synaptic pruning and immune surveillance
    ("Microglia", "RGC", "modulatory"),
    ("Microglia", "Amacrine_Cell", "modulatory"),

    # Astrocyte — perivascular and nerve-fiber-layer ensheathment
    ("Astrocyte", "RGC", "ensheathment"),
]
