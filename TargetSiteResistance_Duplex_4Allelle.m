function [resistance,Nend,tau,x,y,t] = TargetSiteResistance_Duplex_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis)
                                                                            
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



% This code runs a single instance of Wright-Fisher diecious simulations
% for 2 gRNAs, where each cleavage target site has 4 alleles: W - wild type; 
% R - functional resistant; N — non-functional
% resistant; D - drive, but where a mosaic of drive on a single chromosome
% is not possible (e.g. gametes WD is not possible, but only DD)

%Inputs:
% s — female fitness cost of homozygote DD/DD, NX/NX, where X is any
%   allele ~= D and heterozygotes DD/NX
% h — dominance coeff/fitness cost of heterozygote WW/DD in females
% hN — dominance coeff/fitness cost of heterozygotes WW/NN, WW/NX in females
% sigma — female fitness cost of each occurrence of R in a genotype (e.g.
%   w(WR/WR) = (1-σ)^2, w(WW/WR) = 1-σ, w(WR/RR) = (1-σ)^3 )
% Rm — absolute or intrinsic growth rate of population assuming population
%   is all WW
% K — carrying capacity of population assuming population is all WW
% xD — initial frequency of drive in males (DD)
% yD — initial frequency of drive in females (DD)
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
%   For both of these the order of alleles is W R N D — so the 3rd row of x 
%   is the time series of the freq(N)    
% t — corresponding vector of times (generations)

%N.B. note that in paper frequency of males have frequency y and females x

%For m=2 (duplex) need to decided how to order our super-alleles/haplotypes
%The ordering is chosen such that 
%UpperDiagonal of (W,R,N)'*(W,R,N) = [ WW, WR, WN;
%                                      WR, RR, RN;
%                                      WN, RN, NN ]
%and then we take the column order of upper diagonal elements of this rank deficient matrix to give
%the super allele ordering:
%[WW, WR, RR, WN, RN, NN]
% then including drive:
%[WW, WR, RR, WN, RN, NN, DD]


%Fraction of males
psi=0.5;


test_allele = [3,5,6];%RR RN, and NN are resistant
test_threshold = 0.95;

%Density dependent parameter from Beverton-Holt model of population
%dynamics
alpha = K/(Rm-1);


%Dominance coefficients
hD = h;




%Multiplex level
m = 2;

%Number of alleles is n=4 — [W R N D]
n = 4;

%Number of super alleles (haplotypes) ns
%This is the number of combinations (unordered) of 3-alleles  (W,R,N)
%across m sites & is given by n^(m)/m! - where n^(m) is the rising
%factorial
% + 1 drive allele, since drive replaces all allelles at each site
ns =  SumTensorDiagonal(n-1,m) +1;

% %Number of super-alleles including the possibility of drive occurring
% %independently on each site — currently not used?
% Ns = SumTensorDiagonal(n,m);


%Create linear index for upper diagonal of nxn array and (n-1)x(n-1) array

% %nxn
% ii = 1:n^2;
% [i, j] = ind2sub([n,n],ii);
% 
% indn = find(j>=i);

clear ii i j
%(n-1)x(n-1)
ii = 1:(n-1)^2;
[i, j] = ind2sub([n-1,n-1],ii); 
%e.g. for n=3, n-1=2, and i = [1 2 1 2] and j = [1 1 2 2], which represent
%the subscript indices for a 2x2 matrix

%ind_ud (Upper Diagonal) gives an index that filters i and j to give only upper diagonal
%elements — this is useful since for a square matrix A, A(ind_ud) returns
%a 1D row vector of the upper diagonal elements
ind_ud = find(j>=i);

%% Calculate genotype fitness matrices
Wm = ones(ns);
Wf = zeros(ns);

% Heuristic dominance model


%nnR will count the number of R copies there are in each possible haplotype
%WW
nnR = zeros(n-1);
nnW = zeros(n-1);

for i = 1:n-1
    for j=1:n-1
        
        if j>=i
            
            u = [i,j];
            
            % nn is a (n-1)x1 vector of the number of times that allele appears in index u 
            for kk=1:n-1
                nn(kk) = numel(find(u==kk));
            end
            
%             nn
            
            if nn(3) == 0 %i.e. there are no N alleles in haplotype u=[i j]
                nnR(i,j) = nn(2);
                nnW(i,j) = nn(1);
            else
                nnR(i,j) = NaN;
            end
        end
        
    end
end


% Note that the NaNs in nnR are simply to indicate below whether one or both
% chromosomes of diploid genotype has N or D in them, and are treated
% accordingly below

%convert to linear array (column vector)
nR = nnR(ind_ud)';
nW = nnW(ind_ud)';

%Add last element refering to DD haplotype
nR = [nR;NaN];
nW = [nW;NaN];



for i=1:ns
    for j=1:ns
        
        
        
        if i==ns || j==ns
            
            if isnan(nR(i)) && isnan(nR(j))
                Wf(i,j) = 1-s;
            elseif isnan(nR(i))
                if nW(j)==0 %by process of elimination this means resistant RR haplotype and so no somatic costs
                    Wf(i,j) = (1-sigma)^nR(j)*(1-hN*s);
                else
                    Wf(i,j) = (1-sigma)^nR(j)*(1-hD*s);
                end
            elseif isnan(nR(j))
                if nW(i)==0
                    Wf(i,j) = (1-sigma)^nR(i)*(1-hN*s);
                else
                    Wf(i,j) = (1-sigma)^nR(i)*(1-hD*s);
                end
            end
            
        elseif isnan(nR(i)) && isnan(nR(j))
            
            Wf(i,j) = 1-s;
            
        elseif isnan(nR(i))
            Wf(i,j) = (1-sigma)^nR(j)*(1-hN*s);
        elseif isnan(nR(j))
            Wf(i,j) = (1-sigma)^nR(i)*(1-hN*s);
        else
            Wf(i,j) = (1-sigma)^nR(i)*(1-sigma)^nR(j);
        end
        
    end
end



%% Mutation matrix

muvec = [1-mu;xi*mu; (1-xi)*mu];
uR = [0;1;0];
uN = [0;0;1];

MWW = muvec*muvec'+tril(muvec*muvec',-1)';
MuWW = [MWW(ind_ud)';0];

MWR = uR*muvec' + tril(uR*muvec',-1)';
MuWR = [MWR(ind_ud)';0];

MWN = uN*muvec' + tril(uN*muvec',-1)';
MuWN = [MWN(ind_ud)';0];

MRR = uR*uR' + tril(uR*uR',-1)';
MuRR = [MRR(ind_ud)';0];

MRN = uR*uN' + tril(uR*uN',-1)';
MuRN = [MRN(ind_ud)';0];


MNN = uN*uN' + tril(uN*uN',-1)';
MuNN = [MNN(ind_ud)';0];

MuDD = zeros(ns,1);
MuDD(ns)=1;


Mu = [MuWW,MuWR,MuRR,MuWN,MuRN,MuNN,MuDD] - eye(ns);


%% Drive heterozygotes

%Set up drive heterozygotes - Drive allele is now last for allele vector

%For duplex the only genotypes involved in drive are WW\DD, WR\DD, WN\DD -
%the frequency of DD will always be the last element of super-allele
%frequency vector. The frequency of WW, WR, and WN depend on the ordering
%of super-alleles in the super-allele vector:

%The ordering is chosen such that 
%UpperDiagonal of (W,R,N)'*(W,R,N) = [ WW, WR, WN;
%                                      0,  RR, RN;
%                                      0,  0,  NN ]
%and then we take the column order of non-zero elements of this rank deficient matrix to give
%the super allele ordering:
%[WW, WR, RR, WN, RN, NN]
% then including drive:
%[WW, WR, RR, WN, RN, NN, DD]
%So the heterozygotes 
%WW\DD -> 1 & 7
%WR\DD -> 2 & 7
%WN\DD -> 4 & 7
H = [1,2,4;7,7,7];


% Calculate stoichiometry vectors

%This vector kappa is the fraction of gametes expected from W at a single
%site with drive on the other chromosome
kappa = [1-epsilon ; epsilon*nu*beta;epsilon*nu*(1-beta);epsilon*(1-nu)];

%Take outer product - for case of WW\DD 
Kappa  = kappa*kappa';
%Keep only elements on upper diagonal, but weight by their overall
%occurence in matrix - weighting is given by the multinomial coefficient
%(m,[m1,m2,m3...,mn]), where m is the number of indices or level of
%multiplexing and m1 is the number of times allele 1 (W) occurs in the index,
%m2 is the number of times allele 2 occurs in the index (R)...
%in this case m=2 so it is just the binomial coeff

for j = 1:numel(kappa)
    for k = 1:numel(kappa)
        if k>=j
            u = [j,k];
            
            fact = 1;
            
            for kk=1:numel(kappa)
                a = numel(find(u==kk)); %number of times kk appears in the index j,k
                fact = fact*factorial(a);
            end
            
            clear a 
            
            Kappa(j,k) = factorial(m)/fact*Kappa(j,k);
        else
            Kappa(j,k) = 0;
        end
    end
end


%We now need to reduce the size of this matrix by removing the drive rows &
%columns  - and summing over these to give the overall fraction of gametes
%that give DD

% Sum over last row
SDD = sum(Kappa(:,end));

%Reduce size
KKappa = Kappa(1:end-1,1:end-1);

SS = KKappa(ind_ud);

%S is now the stochiometry vector for WW\DD
SS = [SS';SDD];

SS(1) = SS(1)-1; %Stoichiometry vector is the change from Mendelian due to drive

%Stoichiometry for WX\DD - this is effectively kappa, but need to have
%zeros for certain alleles that aren't produced by these reactions 
%WX\DD will produce: WX, RX, NX, & DX = DD
%The ordering of super-alleles is WW,WR, RR, WN, RN, NN, DD
SR = [0;kappa(1)-1;kappa(2);0;kappa(3);0;kappa(4)];
SN = [0;0;0;kappa(1)-1;kappa(2);kappa(3);kappa(4)];

S = [SS,SR,SN]; %Check dimensions...



%% Initial Frequencies
x0 = zeros(ns,1);
y0 = zeros(ns,1);

%For preexisting=1, resistance alleles are potentially in mutation 
%selection balance — for m>1 the distribution of the frequency of WR, RR is 
%not easy to calculate so run each simulation for time 1/σ to implicitly give
%the frequency distribution of resistance alleles when drive is introduced.
%We can calculate approximately the mean frequency and we use this as the 
%starting frequency of this equilbration phase


if preexisting==1
    
    
    %Calculate approximate mean frequency of resistance alleles in mutation
    %selection balance.

    %N.B. That following matrix inversion approach is only valid if all the
    %R containing haplotypes are sufficiently deleterious that they should
    %be at small frequency before drive is introduced, such that non-linear
    %selection terms can be ignored
    %This approach is also only valid as long as the R alleles have no
    %dominance (h=1/2) — although it should be possible to extend this
    %matrix linear approach to include dominance..(?)
    indm = find(Mu(:,1)~=0);
    indm = indm(2:end);
    
    muin = Mu(indm,1);

    
    f = diag((Wf(1,indm)-1)/2);
    
    
    MMu = Mu(indm,indm);
    
    x0(indm) = -inv(MMu+f)*muin;
    y0(indm) = -inv(MMu+f)*muin;
    
    

    %Run simulations for 1/sigma generations *without drive*

    x0(ns) = 0;
    x0(1) = 1-sum(x0); 

    y0(ns) = 0;
    y0(1) = 1-sum(y0);
    
    Tpre = round(1/sigma);
    
    [t,x,y,z,M,F,N,wm,wf,resistance] = TargetSiteResistanceSims_while(ns,K,psi,x0,y0,Wm,Wf,Mu,S,H,Rm,alpha,Tpre,Gauss, test_allele, inf, Cswitch);
    
    
    
    x0 = zeros(ns,1);
    y0 = zeros(ns,1);
    
    x0(indm) = x(indm,end);
    y0(indm) = y(indm,end);
    
    
    
    
    
    x0(7) = xD;
    x0(1) = 1-sum(x0); 

    y0(7) = yD;
    y0(1) = 1-sum(y0);
    
    
    if ~isempty(find(x0<0, 1))
        disp('Error: frequency vector must not have negative entries')
%         input('')
    end
 
    if ~isempty(find(y0<0, 1))
        disp('Error: frequency vector must not have negative entries')
%         input('')
    end
    
    
    

else
    
    x0(7) = xD;
    y0(7) = yD;
    
    x0(1) = 1-xD;
    y0(1) = 1-yD;
    
end






[t,x,y,z,M,F,N,wm,wf,resistance] = TargetSiteResistanceSims_while(ns,K,psi,x0,y0,Wm,Wf,Mu,S,H,Rm,alpha,T,Gauss, test_allele, test_threshold, Cswitch);


Nend = N(end);
tau = t(end);

if plotfig ==1
    
    if vis==0
        f = figure('visible','off')
    else
        figure
    end
    
    col=get(gca,'ColorOrder');
    

    subplot(2,1,1)
    
    for j=1:ns-1
%         jj=j-1;
        plot(t,x(j,:),'Color',col(j,:), 'LineWidth',1,'LineStyle','-');hold on
        plot(t,y(j,:),'Color',col(j,:),'LineWidth',2,'LineStyle','-');hold on
    end
    
    colD = 'k';
    plot(t,x(ns,:),'Color',colD, 'LineWidth',1,'LineStyle','-');hold on
    plot(t,y(ns,:),'Color',colD, 'LineWidth',2,'LineStyle','-');hold on
  
    plot(t(1:end-1),wm(1:end-1),'Color','k','LineWidth',1,'LineStyle','--');
    plot(t(1:end-1),wf(1:end-1),'Color','k','LineWidth',2,'LineStyle','--');
    
    
    hold off
  
    
    title(['$m=2$ gRNAs: $N = ',num2strpow(K),'; s = ',num2str(s),'; h = ',num2str(h),'; h_N = ',num2str(hN),'; \sigma = ',num2str(sigma),...
        '; \nu=',num2strpow(nu),'; \mu = ',num2strpow(mu),...
        '; \beta = ',num2strpow(beta),'; \xi = ',num2strpow(xi),'$'])
%     title(['$N = ',num2strpow(K),'$'])
    
%     set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'FontSize',22)
    
    xlabel('Time $t$ (generations)')
    ylabel('Frequency')
    
    %Legend
    AS = 'WRN';
    
    HS =SuperAlleleString(AS);
    
    
    for k=1:ns
        
        legendstr{2*k-1} = ['Male ',HS{k}];
        legendstr{2*k} = ['Female ',HS{k}];
            
    end
    
    legendstr{2*ns+1} = 'Mean Male Fitness';
    legendstr{2*ns+2} = 'Mean Female Fitness';
    
    hleg=legend(legendstr);
    hleg.FontSize = 12;
    xlim([0 max(t)])
        
    
    



    subplot(2,1,2);
    plot(t,M,'LineWidth',2);hold on
    plot(t,F,'LineWidth',2);hold on
    plot(t,N,'k','LineWidth',2)

    plot(t,K*ones(size(t)),'k--');hold off
    
    xlabel('Time $t$ (generations)')
    ylabel('Population size')
        
 

    legend('Male','Female','Total','Intrinsic carrying capacity')
    
    ylim([0.1 1.2*K])
    xlim([0 max(t)])
    set(gca,'FontSize',22)
    
    if vis==0
        savefig(f,['TargetSiteResistance_Duplex_4Allele_N=',num2strexp(K),'_s=',num2str(s),'_h=',num2str(h),'_hN=',num2str(hN),'_sigma=',num2str(sigma),...
            '_epsilon=',num2str(epsilon),'_nu=',num2str(nu),'_beta=',num2str(beta),'_xi=',num2str(xi),'.fig'])
    end
    
end

end




function S = SumTensorDiagonal(n,d)

%Number of Alleles n
%Number of dimensions d

S = gamma(n+d)/gamma(n)/factorial(d);



end

function SS = SuperAlleleString(a)

S = OuterProductStrings(a);

SS = [];

for k=1:numel(a)
    for j=1:numel(a)
        
        if k>=j
            SS = horzcat(SS,S(j,k));
        end
    end
end


SS = horzcat(SS,'D');

end


function A = num2strexp(A)

A = num2str(A,'%10.1e');

end

function S = num2strpow(s)


if s==0
    
    S = num2str(s);
    
else
    
    pow = floor(log10(s));

    if abs(pow) <5
        S = num2str(s);
    else
        C = round(10^(log10(s)-pow),3);

        if C==1
            S = ['10^{',num2str(pow),'}'];
        else
            S = [num2str(C),'\times',' 10^{',num2str(pow),'}'];
        end

    end


end

end

