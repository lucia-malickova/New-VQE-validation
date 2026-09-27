# Submission v10

Final pre-submission scientific audit.

## New representability test
The canonical 121 parent generators were kept in exactly the same order and
their amplitudes were optimized directly for overlap with the lifted projected
state.

- Four unconstrained-fidelity starts from the best spin-constrained energy
  basin and small perturbations all converge to F_proj = 0.988268496,
  p_D = 0.997155 and approximately 9.875 mEh energy error.
- With the hard constraint p_D >= 0.999, three tested starts converge to
  F_proj = 0.987416189 and approximately 10.817 mEh energy error.
- The analytic fidelity gradient was checked by central finite differences at
  six amplitudes; the maximum absolute discrepancy was 2.8e-11.

These are robust local optima, not certified global maxima. The manuscript
therefore does not claim that the fixed parent sequence is fundamentally
incapable of reproducing the projected state. The supported conclusion is
narrower: same-parameter identification fails catastrophically, and neither
energy nor fidelity reoptimization of the fixed order reproduces the projected
construction in the searches performed here.

## Editorial cleanup
- SLSQP numerical settings are stated explicitly in the SI.
- Remaining "ablation" wording was replaced by "comparison".
- "physical generator" was replaced by "full-space fermionic/parent generator"
  where hardware-native meaning could be inferred.
- Data Availability now matches the public repository.
