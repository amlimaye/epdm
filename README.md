# epdm
Implementation of the Enhanced Partial Propensity Direct Method for Gillespie stochastic simulation

An implementation of the EPDM method (algorithm described in this paper: https://arxiv.org/abs/1609.06403) for Gillespie stochastic simulation. Used for a personal project reproducing the results found in this paper: http://www.pnas.org/content/114/36/E7460.abstract.

Simulation is driven by a C++ backend that produces a JSON file with results from a run. This JSON file is parsed Python-side (see: `scripts/parse.py`) for analysis and plotting (see: `plots/flory_distribution.png`).
