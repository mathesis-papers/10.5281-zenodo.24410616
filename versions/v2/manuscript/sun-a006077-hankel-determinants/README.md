# Sun's A006077 Hankel determinants

The note proves Conjecture 3.16 of Zhi-Wei Sun's
*On some determinants with Legendre symbol entries*, arXiv:1308.2900v5,
p. 17. Its title credits Krattenthaler's earlier zero-class identity,
as recorded by Sun in arXiv:1308.2900v6, Remark 3.16.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is included in
`note.tex`; there is no BibTeX database or external figure dependency.

## Mathematical coverage

The note contains an all-index written proof: binomial translation, residue
blocks, integral normalization, binary determinant recurrences, and positive
moments. Its verification section links the six pinned Lean files and maps
their declarations to the precise statements. The additional real-translation
and square-class assertions are written consequences of the displayed identities.
The same residue decomposition determines all eigenvalue signs and gives
nullity exactly one in the zero class. A scale parameter multiplies the
nonzero determinant by `t^(n(n+1)/3)` without changing its inertia for
nonzero real `t`.

Separate complete moment proofs establish Sun v5 Conjecture 3.9(i) for
every positive real harmonic order and both Bernoulli/Euler-square
inequalities of Conjecture 3.17. The order-one harmonic evaluation is
credited to the dated Petrov–Cigler argument and Amdeberhan's independent
proof; the bibliography links the exact contributions. The source coverage
table distinguishes prior integrality from the remaining positivity assertions
for generalized Franel and Apéry determinants. Edgar's positive density
settles the latter's zeta-three sequence, and Adamczewski's Table 1 supplies
the two prime congruences of Conjecture 3.12. The parity ingredients are
credited chronologically to Bacher, Zhu–Sun, Cigler and Guo. It includes an elementary proof
of the order-two Franel assertion.
The prime-square congruence in the next version of Sun's conjecture is separate.

## Exact finite calculations

Python 3 needs no third-party packages:

```sh
python3 verify.py > finite-check.json
```

`finite-check.json` records exact fraction-free determinant calculations from
the literal finite sum for indices 0 through 32, five translation parameters,
coefficient parity through 512, both binary families through size 64, and
integral normalization through index 16. These finite identities supplement
the proof for all indices in the note. Exact rational congruence elimination
also checks inertia through index 20 at three translations and six scales,
including zero. Rational determinant checks cover harmonic orders 1 through
6 at indices 1 through 10 and both squared-number sign assertions at indices
1 through 10. The order-two Franel determinant is checked through index 16.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (`cc-by-4.0`).
