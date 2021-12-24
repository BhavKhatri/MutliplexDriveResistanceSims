function [resistance,Nend,tau,x,y,t] = TargetSiteResistance_Triplex_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis)


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
% for 3 gRNAs, where each cleavage target site has 4 alleles: W - wild type; 
% R - functional resistant; N — non-functional
% resistant; D - drive, but where a mosaic of drive on a single chromosome
% is not possible (e.g. gametes WRD is not possible, but only DDD)

%Inputs:
% s — female fitness cost of homozygote DDD/DDD, NNX/NNX,NXY/NXY where X and Y are any
%   allele ~= D and heterozygotes DDD/NNX, DDD/NXY
% h — dominance coeff/fitness cost of heterozygote WWW/DDD in females
% hN — dominance coeff/fitness cost of heterozygotes WWW/NNN, WWW/NNX,WWW/NXY in females
% sigma — female fitness cost of each occurrence of R in a genotype (e.g.
%   w(WWR/WWR) = (1-σ)^2, w(WWW/WWR) = 1-σ, w(WRR/RRR) = (1-σ)^5 )
% Rm — absolute or intrinsic growth rate of population assuming population
%   is all WWW
% K — carrying capacity of population assuming population is all WWW
% xD — initial frequency of drive in males (DDD)
% yD — initial frequency of drive in females (DDD)
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

%For m=3 (triplex) need to decided how to order our super-alleles/haplotypes
%The ordering is chosen such that we take the  the outer/kronecker product 
%(W,R,N)*(W,R,N)*(W,R,N) 
%and then from the "upper diagonal" (i<=j<=k) of the resulting 3D array 
%we take the column order of these elements of this rank deficient matrix to give
%the super allele ordering:
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN]
% then including drive:
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN, DDD]



%Fraction of males
psi=0.5;

% test_allele are RRR RRN, RNN and NNN, which are resistant — their indices
% are calculated below
test_threshold = 0.95;

%Density dependent parameter from Beverton-Holt model of population
%dynamics
alpha = K/(Rm-1);


%Dominance coefficients
hD = h;


%Multiplex level
m = 3;

%Number of alleles is n=4 — [W R N D]
n=4;

%Number of super alleles (haplotypes) ns
%This is the number of combinations (unordered) of 3-alleles  (W,R,N)
%across m sites & is given by n^(m)/m! - where n^(m) is the rising
%factorial
% + 1 drive allele, since drive replaces all allelles at each site
ns =  SumTensorDiagonal(n-1,m) +1;

% %Number of super-alleles including the possibility of drive occurring
% %independently on each site — currently not used?
% Ns = SumTensorDiagonal(n,m);



%Create linear index for "upper diagonal" of (n-1)x(n-1)x(n-1) array

% %nxn
% ii = 1:n^2;
% [i, j] = ind2sub([n,n],ii);
% 
% indn = find(j>=i);

% clear ii i j


%(n-1)x(n-1)x(n-1)
ii = 1:(n-1)^m;
[i, j, k] = ind2sub([n-1,n-1,n-1],ii); 
%e.g. for m=2 and n=3, n-1=2, and i = [1 2 1 2] and j = [1 1 2 2], which represent
%the subscript indices for a 2x2 matrix

%ind_ud (Upper Diagonal) gives an index that filters i and j to give only upper diagonal
%elements — this is useful since for a square matrix A, A(ind_ud) returns
%a 1D row vector of the upper diagonal elements
ind_ud = find(k>=j & j>=i);

%% Calculate genotype fitness matrices
Wm = ones(ns);
Wf = zeros(ns);

% Heuristic dominance model


%nnR will count the number of R copies there are in each possible haplotype
%WW
nnR = zeros(n-1,n-1,n-1);
nnW = zeros(n-1,n-1,n-1);


for i = 1:n-1
    for j=1:n-1
        for k=1:n-1
            
        
            if k>=j && j>=i

                u = [i,j,k];

                % nn is a (n-1)x1 vector of the number of times that allele appears in index u 
                for kk=1:n-1
                    nn(kk) = numel(find(u==kk));
                end

    %             nn

                if nn(3) == 0 %i.e. there are no N alleles in haplotype u=[i j]
                    nnR(i,j,k) = nn(2);
                    nnW(i,j,k) = nn(1);
                else
                    nnR(i,j,k) = NaN;
                end
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

%Add last element refering to DDD haplotype
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

%For triplex I need to loop over all haplotypes and count the number of W R
%and N alleles in each haplotype

muvec = [1-mu;xi*mu; (1-xi)*mu];
uR = [0;1;0];
uN = [0;0;1];

nalleles = zeros(ns-1,n-1);

for jj=1:m
    ndims(jj) = n-1;
end

M = zeros(ns);

for l=1:ns-1
    [i, j, k] = ind2sub(ndims,ind_ud(l));
    
    u = [i,j,k];
    
    % nn is a (n-1)x1 vector of the number of times that allele appears in index u 
    for kk=1:n-1
        nalleles(l,kk) = numel(find(u==kk));
    end
    
    V=1;
    for r=1:3
        
        switch r
        
            case 1
                mvec = muvec;        
            
            case 2
                mvec = uR;
        
            case 3
                mvec = uN;
                        
            
        end
                
            
        
         for rr=1:nalleles(l,r) %If nalleles=0 then this should be by-passed
             V = kron(V,mvec);
         end
         
    end
         
    %Sum over rows of nalleles should always be n-1 and so v should be an 3
    %dimensional for triplex and m-dimensional in general - but flattened to 2D
    
    V = reshape(V,ndims);
    
    %Do analogue of suming transpose of matrix (V+V^T), of non-diagonal
    %part only!
    V = Symmetrise3DTensor(V);

    v = V(ind_ud);
    
    M(:,l) = [v';0];
    
    
end

M(ns,ns) = 1;

Mu = M-eye(ns);




%% Drive heterozygotes — Triplex

%Set up drive heterozygotes - Drive allele is now last for allele vector

%For triplex, the ordering of haplotypes is using column order of the
%"upper diagonal" (i<=j<=k) of a 3D array as described above this results
%in the following order
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN, DDD]

%Now drive only affects gamete production for those heterozygotes with "DDD" 
%that have at least one W in their haplotype, which are 
%WWW/DDD, WWR/DDD, WRR/DDD, WWN/DDD, WRN/DDD, WNN/DDD
%However, approach we'll take is to calculate all combinations using base
%stoichiometry vectors kappa, kR, kN, as below

% Calculate stoichiometry vectors

%This vector kappa is the fraction of gametes expected from W at a single
%site with drive on the other chromosome
kappa = [1-epsilon ; epsilon*nu*beta;epsilon*nu*(1-beta);epsilon*(1-nu)];

%These are the fractions of gametes expected for R and N respectively
kR = [0 1 0 0]';
kN = [0 0 1 0]';
% kD = [0 0 0 1];

%Here for triplex we need ndims4 = [n,n,n], as well as ndims =
%[n-1,n-1,n-1]
for jj=1:m
    ndims4(jj) = n;
end

S = zeros(ns);

%want to loop over all possible haplotypes which could be heterozygote with
%DDD (e.g. XYZ/DDD) where {X,Y,Z} = {W,R,N}, so not including D — since XYD
%is assumed to always be converted to "DDD" and so cannot arise
%So will use linear indexing over ns-1 haplotypes not inc D


for l=1:ns-1
    
    [i, j, k] = ind2sub(ndims,ind_ud(l)); %(e.g. l=3-> u=[1,2,2] =WRR
    
    u = [i,j,k];
    
    % nalleles is a (ns-1)x(n-1) vector:nalleles(l,kk) is number of times that allele kk appears in index u/haplotype l (ell not one) 
    for kk=1:n-1
        nalleles(l,kk) = numel(find(u==kk));
    end
    
    %e.g. for l=3-> u=[1,2,2] =WRR, we have nalleles(l,:) = [1,2,0]
    
    V=1;
    for r=1:n-1
        
        switch r
        
            case 1
                kvec = kappa;   
            
            case 2
                kvec = kR;
        
            case 3
                kvec = kN;
                
                     
            
        end
                
         
        
         for rr=1:nalleles(l,r) %If nalleles=0 then this should be by-passed
             V = kron(V,kvec);
         end
         
    end
         

    %Sum over rows of nalleles should always be n-1 and so v should be an 3
    %dimensional for triplex and m-dimensional in general - but flattened to 2D
    
    V = reshape(V,ndims4);
    
        
    %Do analogue of suming transpose of matrix (V+V^T), of non-diagonal
    %part only!
    V = Symmetrise3DTensor(V); %nxnxn tensor
    
    %Drive entry corresponds to the "last" page of the tensor V, and only
    %the "upper diagonal" (i<=j<=k) 
    
    ii = 1:(n)^m;
    [i, j, k] = ind2sub([n,n,n],ii);

    ind_ud4 = find(k>=j & j>=i & (i==4 | j==4 | k==4));
    
    SD = sum(V(ind_ud4));
    
    %Reduce to (n-1)x(n-1)x(n-1) to focus on haplotypes that have only W,R,N
    V = V(1:n-1,1:n-1,1:n-1);
    
    %Find "upper diagonal"
    v = V(ind_ud);
    
    S(:,l) = [v';SD];
    
    
end


S(ns,ns) = 1;
S = S-eye(ns);


%H are the columns of S which correspond to the heterozygotes we listed
%above that have at least one W in their haplotype
H = find(sum(abs(S))~=0);

%The complement of this also defines which haplotypes are resistant
R = setdiff(1:ns-1,H);
test_allele = R;

S = S(:,H);

H = [H;ns*ones(size(H))];





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
    
    x0(ns) = xD;
    x0(1) = 1-sum(x0); 

    y0(ns) = yD;
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
    
    x0(ns) = xD;
    y0(ns) = yD;
    
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
    c = distinguishable_colors(ns-7);
  
    col = [col;c];

    subplot(2,1,1)
    
    for j=1:ns-1
%         jj=j-1;
        plot(t,x(j,:),'Color',col(j,:), 'LineWidth',1,'LineStyle','-');hold on
        plot(t,y(j,:),'Color',col(j,:),'LineWidth',2,'LineStyle','-');hold on
    end
    
    colD = 'k';
    plot(t,x(ns,:),'Color',colD, 'LineWidth',1,'LineStyle','-');hold on
    plot(t,y(ns,:),'Color',colD, 'LineWidth',2,'LineStyle','-');hold on
  

    plot(t(1:end-1),wm(1:end-1),'Color','k','LineWidth',1,'LineStyle','-.');
    plot(t(1:end-1),wf(1:end-1),'Color','k','LineWidth',2,'LineStyle','-.');
    
    
    hold off
   
    xlabel('Time $t$ (generations)')
    ylabel('Frequency')

    
    title(['$m=3$ gRNAs: $N = ',num2strpow(K),'; s = ',num2str(s),'; h = ',num2str(h),'; h_N = ',num2str(hN),'; \sigma = ',num2str(sigma),...
        '; \nu=',num2strpow(nu),'; \mu = ',num2strpow(mu),...
        '; \beta = ',num2strpow(beta),'; \xi = ',num2strpow(xi),'$'])
%     title(['$N = ',num2strpow(K),'$'])
%     set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'FontSize',22)
    
    %Legend
%     AS = 'WRN';
%     
%     HS =SuperAlleleString(AS);
%     


        
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN, DDD]

HS{1} = 'WWW';
HS{2} = 'WWR';
HS{3} = 'WRR';
HS{4} = 'RRR';
HS{5} = 'WWN';
HS{6} = 'WRN';
HS{7} = 'RRN';
HS{8} = 'WNN';   
HS{6} = 'WRN';
HS{7} = 'RRN';
HS{8} = 'WNN';  
HS{9} = 'RNN';
HS{10} = 'NNN';
HS{11} = 'D';  


    
    for k=1:ns
        
        legendstr{2*k-1} = ['Male ',HS{k}];
        legendstr{2*k} = ['Female ',HS{k}];
            
    end
    
    legendstr{2*ns+1} = 'Mean Male Fitness';
    legendstr{2*ns+2} = 'Mean Female Fitness';
    
    hleg = legend(legendstr, 'NumColumns',2);
    hleg.FontSize = 12;
    xlim([0 max(t)])


    subplot(2,1,2);
    plot(t,M,'LineWidth',2);hold on
    plot(t,F,'LineWidth',2);hold on
    plot(t,N,'k','LineWidth',2)
    plot(t,K*ones(size(t)),'k--');hold off
    
    xlabel('Time $t$ (generations)')
    ylabel('Population size')
    
    legend('Male','Female','Total','Intrinsic carrying capacity','Location','southwest')
    
    ylim([0.1 1.2*K])
    xlim([0 max(t)])
    set(gca,'FontSize',22)
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')

    if vis==0
        savefig(f,['TargetSiteResistance_Triplex_4Allele_N=',num2strexp(K),'_s=',num2str(s),'_h=',num2str(h),'_hN=',num2str(hN),'_sigma=',num2str(sigma),...
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




function S = Symmetrise3DTensor(A)


%This produces a symmetrised tensor, where the sum over the 
%"upper diagonal" i<=j<=k, gives the same as the sum over the whole original tensor 
%This is analogous to S = A'+tril(A,-1)' in 2D, i.e. folding back and
%summing the off-diagonal elements leaving diagonal untouched.


%Size of tensor
[N1,N2,N3] = size(A);

if N1~=N2 | N1~=N3 | N2~=N3
    disp('Error: tensor should be "square" (cube)')
else
   nD=N1;
end

%This needs changing to account for different redundancies of permutations

p = perms([1 2 3]);
[N,~] = size(p);

% N
%Extract diagonal
for n=1:nD
    D(n) = A(n,n,n);
end

S = zeros(nD,nD,nD);

%Sum over all permutations
for n=1:N
    S = S + permute(A,p(n,:));
end


%This overcounts for redundant pairs of indices

for i=1:nD
    for j=1:nD
        for k=1:nD
            u = [i,j,k];
            
            %If any two indices are the same, but 3rd not
            
            if i==j & (k~=i & k~=j) | i==k & (j~=i & j~=i) | j==k & (i~=k & i~=k)
                
                S(i,j,k) = S(i,j,k)/2;
            end
        end
    end
end

                




% S = S/6;


for n=1:nD
    S(n,n,n) = D(n);
end

end





function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end


