function [resistance,Nend,tau,x,y,t] = TargetSiteResistance(m,s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis)

% MIT License
% 
% Copyright (c) 2021 Bhavin S Khatri
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%%Main function

% This code runs a single instance of Wright-Fisher diecious simulations
% for m(<4) gRNAs, where each cleavage target site has 4 alleles: W - wild type; 
% R - functional resistant; N — non-functional
% resistant; D - drive, but where a mosaic of drive on a single chromosome
% is not possible (e.g. gametes WD is not possible, but only DD)

%Inputs:
% m — number of gRNAs (m<4)
% s — female fitness cost of homozygotes of drive (D..D/D..D) and homozygotes 
%   involving at least a single N in the haplotype/allele or heterozygotes 
%   at least a single N in the haplotype/allele on both chromosomes
% h — dominance coeff/fitness cost of heterozygotes of W..W/D..D in females
% hN — dominance coeff/fitness cost of heterozygotes W..W/N..N, W..W/N..X in females
% sigma — female fitness cost of each occurrence of R in a genotype (e.g.
%   w(WR/WR) = (1-σ)^2, w(WW/WR) = 1-σ, w(WR/RR) = (1-σ)^3 )
% Rm — absolute or intrinsic growth rate of population assuming population
%   is all Ws at all target sites
% K — carrying capacity of population assuming population is all Ws at all target sites
% xD — initial frequency of drive in males (D..D)
% yD — initial frequency of drive in females (D..D)
% epsilon — efficiency of drive cleavage per target site
% nu — non-homologous end-joining (NHEJ) rate per generation per individual
%   per target site
% mu — mutation rate per cleavage target site per generation per individual
%   per target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9, mu=1e-8)
% beta — fraction of functional resistant NHEJ mutations
% xi — fraction of functional resistant single nucleotide mutations
% preexisting — preexisting =1 runs simulations for a period 1/sigma as a
%   burn-in period to establish a mutation-selection balance
% T — length of simulations in generations (simulations will end earlier if population eliminated)
% Gauss — Gauss = 0 uses multinomial distribution; Gauss = 1 uses
%   Poisson-Gaussian hybrid approximation to the multinomial distribution.
%   If Gauss=0 the maximum population size is ~N=1e9 and is very slow — if
%   Gauss=1, the hybrid approx can be used with arbitrarily large N with
%   little speed penalty
% Cswitch — Cswitch = 1 if together with Gauss=0 uses GSL C version of
%   multinomial random number generator, which is much quicker than the
%   matlab implementation. Cswitch = 0 uses Matlab's multinomial RNG,
%   unless Gauss=1
% plotfig — plotfig=1 plots allele frequency time series and population
%   dynamics; plotfig=0 doesn't
% vis — vis=0 with plotfig=1 plots figures and makes them invisible and 
%   saves them (for remote jobs); their visibility state will be 0 and so when
%   opened in Matlab this needs to be changed, e.g. openfig('Figure.fig','visible')

%Outputs:
% resistance — resistance=0 or 1, whether simulation produced resistance 
%   (freq(R)>test_threshold); 
%   if resistance alleles have not passed threshold by T, then resistance=0
% Nend — final population size at the end of simulation
% tau — the time when resistance occured (freq(R)>test_threshold) or = NaN
%   if population eliminated
% x — matrix whose columns are the frequency of each allele in males at subsequent
%   time points
% y — matrix whose columns are the frequency of each allele in females at subsequent
%   time points
%   allele ordering in these vector are as described in each of the functions called for m=1,2,3 below 
% t — corresponding vector of times (generations)

%N.B. note that in paper frequency of males have frequency y and females x


switch m
    case 1
        [resistance,Nend,tau,x,y,t] = TargetSiteResistance_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis);
    case 2
        [resistance,Nend,tau,x,y,t] = TargetSiteResistance_Duplex_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis);
    case 3
        [resistance,Nend,tau,x,y,t] = TargetSiteResistance_Triplex_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis);

end


        
        