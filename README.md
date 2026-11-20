# Sharp M-SERPT and SERPT constants for two-point job sizes

Code and outputs for the paper "Sharp ratios of SERPT and M-SERPT to the Gittins policy for two-point job sizes" (Y. Shmalo).

The job size takes the value a with probability p and b > a with probability 1 - p, and jobs arrive in an M/G/1 queue at rate lambda.

## Files

- `verify.py` runs the checks for the paper. `python verify.py exact` does the following:
  - checks the symbolic identities;
  - certifies an interval for the constant C_2 = 1.0601734219... and locates its stationary point (lambda a = 0.45061..., p = 0.65713...);
  - evaluates 5046 rational parameter points exactly;
  - reproduces the SOAP coin-flip example with sizes 2 and 14 at lambda = 1/10 (M-SERPT 51/2, Gittins 101/4, ratio 102/101);
  - checks the degenerate cases;
  - follows the ratio along the extremal families for M-SERPT (to C_2) and SERPT (to (7 - 2 sqrt 2)/4 = 1.04289...).

  It takes about ten seconds. `python verify.py simulate` runs the five discrete-event simulations, with seed 7115 and about 98,000 measured jobs per policy and case, and is much slower. `python verify.py` runs both.
- `identities_check.py` checks with SymPy the identities used in the proofs:
  - the derivative of the quintic;
  - the two partial derivatives of the limiting gap and the factorizations of K;
  - the resultant that yields the quintic;
  - the boundary formula and the SERPT denominator.
- `exact_results.txt`, `identities_output.txt` and `simulation_results.jsonl` are the outputs of `verify.py exact`, `identities_check.py` and `verify.py simulate`.

In all five cases the simulated mean response times of M-SERPT, Gittins and SERPT are within 0.7 per cent of the closed forms. The ratios agree more closely, because the three policies are simulated on the same arrivals:
- In the M-SERPT near-extremal case the M-SERPT/Gittins ratio is 1.06017 in theory and 1.06007 in simulation.
- In the SERPT near-extremal case (b = 1001) the SERPT/Gittins ratio is 1.04259 in theory and 1.04225 in simulation.
- At the transition p(b - a) = a both ratios are 1.

## Requirements

Python 3 with SymPy and mpmath. The outputs here were produced with Python 3.14 and SymPy 1.14.

The paper's source and PDF are in `paper/`.
