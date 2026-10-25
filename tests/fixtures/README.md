# Legacy regression fixtures

These R sources are copied unchanged from breaktest 0.3.0, the upstream
snapshot used before the multiple-break CHLT change (commit `96e71b7`).
`original-CHLT.R` includes its original `GLS.bt` helper; the other fixtures
contain only the corresponding original function definition.

They deliberately preserve existing single-break behavior, including its
errors and bootstrap quirks. Tests source them into a separate environment
whose parent is the installed breaktest namespace. The original bootstrap
method is bound explicitly so its dispatch cannot resolve to the new method.
The multi-break numerical references do not use these fixtures or the package
regression, detrending, grid or bootstrap helpers.
