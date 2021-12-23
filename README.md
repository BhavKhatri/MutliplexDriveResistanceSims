# MutliplexDriveResistanceSims
Wright-Fisher combined with Beverton-Holt stochastic simulations of the evolution of multiplex drive resistance. Multinomial sampling has option of fast and accurate hybrid Poisson-Gaussian approximation for simulations at arbitrarily large population sizes. Can simulate up to m=3 gRNAs with NHEJ, de novo mutation and standing variation as options.

Basic model consists of m target sites for cleavage, where each have nomimally the same parameters describing it's evolution (cleavage efficiency, NHEJ rate, mutation rate, and fitness parameters) and the same number of alleles: W — wild type; R — functional resistant; N — non-functional resistant & D — drive. The simulations assume that for m>1, that a single successful cleavage amongst m target sites, followed by homology directed repair is sufficient to copy all m Ds and so for m>1 mosaic alleles WDR is not possible. This means we do not model loss of gRNAs or loss of function of gRNAs.

The simulations are based on the two-sex (diecious) Wright-Fisher model combined with Beverton-Holt population dynamics. In other words, only allele frequencies are tracked (not genotype freqencies) and given allele frequencies in one generation, the next generation allele frequencies are drawn using multinomial sampling with replacement, based on the (relative) fitness, mutation and drive matrices with a fixed population size in a given generation Nt. Beverton-Holt dynamics with absolute fitness of the population and density-dependent competition parameter alpha (tuned to a given carrying capacity/population size before introduction of drive) determines the population size in the next generation. 

These are the simulations used to create results in "Weakly deleterious natural genetic variation amplifies probability of resistance in multiplexed gene drive systems", bioRxiv ..., Bhavin S. Khatri and Austin Burt. (N.B. note that in paper frequency of males is y and females x — here they are reversed)

To run in matlab command prompt:

[resistance,Nend,tau,x,y,t] = TargetSiteResistance(m,s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis)


Inputs:

m — number of gRNAs (m<4)

s — female fitness cost of homozygotes of drive (D..D/D..D) and homozygotes involving at least a single N in the haplotype/allele or heterozygotes at least a single N in the haplotype/allele on both chromosomes 

h — dominance coeff/fitness cost of heterozygotes of W..W/D..D in females

hN — dominance coeff/fitness cost of heterozygotes W..W/N..N, W..W/N..X in females

sigma — female fitness cost of each occurrence of R in a genotype (e.g. w(WR/WR) = (1-σ)^2, w(WW/WR) = 1-σ, w(WR/RR) = (1-σ)^3 )

Rm — absolute or intrinsic growth rate of population assuming population is all Ws at all target sites

K — carrying capacity of population assuming population is all Ws at all target sites

xD — initial frequency of drive in males (D..D)

yD — initial frequency of drive in females (D..D)

epsilon — efficiency of drive cleavage per target site

nu — non-homologous end-joining (NHEJ) rate per generation per individual per target site

mu — mutation rate per cleavage target site per generation per individual per target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9, mu=1e-8)

beta — fraction of functional resistant NHEJ mutations

xi — fraction of functional resistant single nucleotide mutations 

preexisting — preexisting =1 runs simulations for a period 1/sigma as a burn-in period to establish a mutation-selection balance

T — length of simulations in generations (simulations will end earlier if population eliminated)

Gauss — Gauss = 0 uses multinomial distribution; Gauss = 1 uses Poisson-Gaussian hybrid approximation to the multinomial distribution. If Gauss=0 the maximum population size is ~N=1e9 and is very slow — if Gauss=1, the hybrid approx can be used with arbitrarily large N with little speed penalty

Cswitch — Cswitch = 1 if together with Gauss=0 uses GSL C version of multinomial random number generator, which is much quicker than the matlab implementation. Cswitch = 0 uses Matlab's multinomial RNG, unless Gauss=1. This requires compiling supplied c file using mex within Matlab on your platform.

plotfig — plotfig=1 plots allele frequency time series and population dynamics; plotfig=0 doesn't

vis — vis=0 with plotfig=1 plots figures and makes them invisible and saves them (for remote jobs); their visibility state will be 0 and so when opened in Matlab this needs to be changed, e.g. openfig('Figure.fig','visible')
  
Outputs:

resistance — resistance=0 or 1, whether simulation produced resistance (freq(R)>test_threshold); if resistance alleles have not passed threshold by T, then resistance=0

Nend — final population size at the end of simulation

tau — the time when resistance occured (freq(R)>test_threshold) or = NaN if population eliminated

x — matrix whose columns are the frequency of each allele in males at subsequent time points

y — matrix whose columns are the frequency of each allele in females at subsequent time points allele ordering in these vector are as described in each of the functions called for m=1,2,3 below 

t — corresponding vector of times (generations)


