TITLE 5-HT1A-activated G protein-gated inwardly rectifying K+ current (Kir3, GIRK)

COMMENT
Serotonin acting on 5-HT1A receptors hyperpolarizes CA3 pyramidal cells and lowers their firing through Kir3
channels; parvalbumin basket cells are not affected (Johnston, McBain & Fisahn 2014, J Physiol 592:4187,
doi:10.1113/jphysiol.2014.279083). The model inserts it in every pyramidal cell with gmax = 0; run_pilot.py switches
gmax on in the cells 5-HT reaches while 5-HT is present (--ht-girk).

Instantaneous inward rectification of the usual Kir form, relative to E_K:
    i = gmax * (v - ek) / (1 + exp((v - ek - vhalf) / vslope))
vhalf = 0 mV and vslope = 27.7 mV are a least-squares fit to the normalized I-V of the 5-HT-induced current in
Fig. 5B of Johnston et al. 2014 (E_K = -85 mV there; rectification index of the fit -0.41, measured median -0.62).
gmax: run_pilot.py --ht-girk 6e-5 hyperpolarizes a pyramidal cell held near -45 mV by 4.4 mV, against the measured
4.3 +/- 1.6 mV (10 uM 5-HT, active CA3 pyramidal cells near -49 mV, n = 14; their Fig. 4C).
ENDCOMMENT

NEURON {
    SUFFIX ch_GIRK
    USEION k READ ek WRITE ik
    RANGE gmax, vhalf, vslope, g, i
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gmax = 0 (mho/cm2)
    vhalf = 0 (mV)
    vslope = 27.7 (mV)
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    g (mho/cm2)
    i (mA/cm2)
}

BREAKPOINT {
    g = gmax / (1 + exp((v - ek - vhalf) / vslope))
    i = g * (v - ek)
    ik = i
}
