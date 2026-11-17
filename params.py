"""Every parameter of tinyCA3 with its source and recording conditions.

Status values:
  "measured"     - taken from an experimental paper (conditions given)
  "model"        - taken from a published model (which one, and the conditions it was built for)
  "placeholder"  - a stand-in from another cell type or region, to be replaced before any CA3 claim
  "free"         - chosen here, to be calibrated against an acceptance criterion
No value is entered from memory without the status "placeholder" or "free".
"""
from dataclasses import dataclass

CELSIUS = 34.0  # degC. Bezaire et al. 2016 model temperature (status: model). The STP mod's AMPA kinetics are
                # described as "standard CA1 Schaffer ~30 degC" in its header.

# --- Bezaire et al. 2016 (eLife 5:e18566), CA1 model, github.com/mbezaire/ca1 commit a8655fa, CC0 -------------------
# Pyramidal cell (class_poolosyncell.hoc, `mechinit`/`insert_mechs`): somatic channel densities and passive values.
# Status: model (CA1 pyramidal cell, rat, 34 degC); used here as a PLACEHOLDER for CA3 pyramidal and EC/DG/SWT cells.
BEZ_PYR_SOMA = dict(
    gNav=0.032,        # S/cm2, ch_Navp, soma and basal dendrites
    gKdr=0.003,        # S/cm2, ch_Kdrp
    gKvAprox=0.008,    # S/cm2, ch_KvAproxp
    gHCN=0.0006,       # S/cm2, ch_HCNp (soma value; it rises along the apical tree in the template)
    gNavaxon=0.064,    # S/cm2, ch_Navaxonp in the axon
    Rm=28000.0,        # ohm cm2, soma and axon; dendrites Rm/2
    Cm=1.0,            # uF/cm2, soma and axon; dendrites 2*Cm
    Ra_soma=150.0,     # ohm cm, soma and dendrites
    Ra_axon=50.0,      # ohm cm
    ena=55.0, ek=-90.0, e_hcn=-30.0,
    v_rest=-66.0,      # mV, Vrest and e_pas of the template
)

SOURCES = {
    "BEZ_PYR_SOMA": "Bezaire et al. 2016 eLife 5:e18566, class_poolosyncell.hoc (commit a8655fa); model; CA1, rat, 34 C",
    "BEZ_PV": "Bezaire et al. 2016, class_pvbasketcell.hoc (commit a8655fa); model; CA1 PV basket cell, rat, 34 C",
    "STP_CA1_SC": "TMExp2SynNMDACA1_SCExc.mod header: fit to Wierenga & Wadman 2003 Fig. 4A+B; Wistar rat P14-28, CA1",
}


@dataclass(frozen=True)
class Geometry:
    """Reduced geometries chosen for tinyCA3. Status: free (to be calibrated in T1 against rheobase and input
    resistance targets once CA3 data are sourced)."""
    pyr_soma_L: float = 20.0     # um
    pyr_soma_diam: float = 20.0
    pyr_dend_L: float = 300.0    # one cylinder standing for the apical tree (stratum radiatum)
    pyr_dend_diam: float = 3.0
    pyr_axon_L: float = 50.0     # axon initial segment, as Bezaire axon[0] (L 50, diam 1)
    pyr_axon_diam: float = 1.0
    rs_soma_L: float = 20.0      # EC, DG and SWT cells: one compartment
    rs_soma_diam: float = 20.0
