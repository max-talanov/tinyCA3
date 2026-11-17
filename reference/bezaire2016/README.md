Original cell templates of the Bezaire et al. 2016 CA1 model, unchanged, from github.com/mbezaire/ca1 at commit
a8655fa8614687276b82d88599c818dc29cb522f (`cells/`), CC0.
- `class_pvbasketcell.hoc`, `class_bistratifiedcell.hoc`, `class_olmcell.hoc`: the reference for
  `tests/test_s3b_bezaire_cells.py`, which checks the Python port `cells.build_bezaire_interneuron` against them.
- `class_poolosyncell.hoc` (pyramidal cell, reconstructed morphology, 202 sections): run directly by
  `cells.build_bezaire_pyr` (not ported).
