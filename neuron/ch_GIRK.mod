TITLE 5-HT1A-activated G protein-gated inwardly rectifying K+ current (Kir3, GIRK)

COMMENT
Serotonin acting on 5-HT1A receptors hyperpolarizes CA3 pyramidal cells and lowers their firing through Kir3
channels; parvalbumin basket cells are not affected (Johnston, McBain & Fisahn 2014, J Physiol 592:4187,
doi:10.1113/jphysiol.2014.279083). The model inserts it in every pyramidal cell with gmax = 0; run_pilot.py switches
gmax on in the cells 5-HT reaches while 5-HT is present (--ht-girk).

Instantaneous inward rectification of the usual Kir form, relative to E_K:
    i = gmax * (v - ek) / (1 + exp((v - ek - vhalf) / vslope))
vhalf and vslope are free (tinyCA3 status: free); gmax is calibrated to the resting hyperpolarization (cells.py).
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
    vhalf = 15 (mV)
    vslope = 10 (mV)
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
