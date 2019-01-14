# The Effect of Directional Dominance on Additive Effect Sizes
##### Systems Biology Department, Columbia University
##### W3500 Independent Biological Research

In this repository is the work I did while working as a research assistant with Jeremy Berg in the [Sella Lab](https://sellalab.biology.columbia.edu/) while I was taking the Independent Biology Research class at Columbia. In a nutshell, I derived mathematically an expression for a certain hypothesized bias in a statistical test for polygenic selection, and implemented a bunch of tools in R (and some bash scripts) to investigate the presence of that bias using genome data from the UK Biobank.

Probably the most relevant thing here is the final report for the class, which is ```final.pdf```. Below I also copy the abstract from that paper:

The average effect size $\alpha$ is the estimated effect size for an
allele for a polygenic trait. These effect sizes are estimated using
Genome-Wide Association Studies, and we can estimate signals of natural
selection for polygenic traits by weighting by $\alpha$ in statistical
tests for polygenic adaptation. However, in the presence of
directional dominance, we hypothesized that such statistical tests
could be biased. For example, alleles which are systematically
dominant and have increased in frequency in the recent past will have
deflated average effect sizes, and vice versa. Simulations of the test
statistic $Q_x$ used to detect polygenic adaptation, expanded to
account for dominance effects, show that the effect that directional
dominance has on the statistic is to shift the distribution to the
right. Furthermore, directional dominance causes $\alpha$ to
systematically stretch or shrink genetic value over generations,
further pointing to inflations of deflations of $\alpha$ in accordance
with directional dominance. Preliminary investigations of the UK
Biobank data for height show that height displays both polygenic
adaptation and directional dominance.
