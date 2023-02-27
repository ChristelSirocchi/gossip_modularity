# Community-based gossip algorithm for distributed averaging

This repository contains code and data for the manuscript "**Community-based gossip algorithm for distributed averaging**"

It contains two folders:
* simulator
   * convergence_ABM.py: event-driven simulator based on the Python library Simpy;
   * graph_generators.py: functions generating connected graphs for the four chosen graph families (Erdos-Renyi, geometric random, small world, scale-free);
   * FARZ.py: code for the FARZ modular network generator adapted from https://github.com/rabbanyk/FARZ;
   * graph_metrics.py: functions calculating structural graph metrics (global, local and spectral metrics);
   * interaction_methods.py: functions defining neighbour selection criteria (random, community), type of interaction (averaging, averaging with label propagationa), time to next move (Poisson);
   * modularity_util.py: additional functions to calculate modularity metrics, identify boundary nodes, generate null models for modular networks
* tests: networks generated within the study and results
   * synthetic_set: set of 1000 modular structures generated using the FARZ benchmark (Fagnan et al. 2018)
   * results: average convergence rates of the random and community-based gossip algorithms, generating parameters and structural metrics for the 1000 networks in the sythetic set
   * test_run: a test exemplifying how to generate community graphs, execute the algorithms and compute structural features.
 
