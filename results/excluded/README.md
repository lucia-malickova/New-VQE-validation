# Excluded Cu-S geometry scan

A preliminary five-point Cu-S scan was generated as an exploratory extension. It is **not** part of the validated scientific result.

The campaign validation returned `VERDICT=FAIL`.

Reasons:
- doublet CASSCF did not converge across the five points;
- state-specific quartet CASSCF did not converge across the five points;
- active-orbital continuity dropped below the preset overlap threshold away from the reference geometry.

The validation log is retained for transparency. The scan must not be used as a validated geometry dependence of the doublet-quartet gap or of VQE performance.
