# Submission revision v8

The v8 referee audit supersedes the raw null-space-basis mixing interpretation
used in the earlier v6 sensitivity probe.

The deterministic sparse two-dimensional Cu null-space basis is individually
normalized but not mutually orthogonal. Direct index-space O(2) mixing
therefore rescales generators (up to 1.512 in Frobenius norm and 1.732 in
Hartree-Fock action norm) and is not a clean sensitivity test for
single-operator ADAPT.

The v8 audit evaluates 100 independently mixed and explicitly renormalized
bases at every one of the 51 stored states along the final 18q trajectory,
using both coefficient-vector Euclidean normalization and determinant-space
Frobenius normalization. Under either convention the highest-gradient
null-space group changes at only 4/51 prefixes; the maximum changed-basis
fraction is 0.14 and 0.12, respectively. A basis-invariant group-gradient
score agrees with the canonical top group at 50/51 coefficient-metric prefixes
and 49/51 Frobenius-metric prefixes.

The manuscript therefore narrows the basis-dependence claim: the ambiguity is
formal and detectable, but top-group sensitivity is modest along the stored
trajectory after normalization is controlled. The v6 raw-rotation table is
retained only for provenance and is not used for a scientific claim.

Data Availability is also made consistent between the main text and SI.
