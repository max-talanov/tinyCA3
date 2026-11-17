# Channel mechanisms from the Bezaire et al. 2016 CA1 model

Copied unchanged from github.com/mbezaire/ca1, commit a8655fa8614687276b82d88599c818dc29cb522f (CC0 licence; also
ModelDB 187604). Reference: Bezaire MJ, Raikov I, Burk K, Vyas D, Soltesz I (2016) eLife 5:e18566.

ch_CavL, ch_CavN, ch_HCNolm, ch_KCaS, ch_Kdrfast, ch_KvA, ch_KvAolm, ch_KvCaB, ch_Nav, ch_Navaxonp, ch_Navbis, ch_leak,
iconc_Ca; for the pyramidal cell also ch_Navp, ch_Kdrp, ch_KvAproxp, ch_KvAdistp, ch_HCNp. Used by `cells.build_bezaire_interneuron` (checked against the original hoc templates in
`tests/test_s3b_bezaire_cells.py`) and `cells.build_bezaire_pyr` (which runs the original template). Temperature dependence is per file (e.g. ch_Navaxonp Q10 2 from 24 C, ch_Kdrfast Q10 3
from 34 C, ch_HCNolm and ch_Nav none); the model runs at 34 C.
