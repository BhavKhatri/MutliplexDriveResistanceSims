function [resistance,Nend,tau,x,y,t] = TargetSiteResistance_4Allelle(s,h,hN,sigma,Rm,K,xD,yD,epsilon,nu,mu,beta,xi,preexisting,T,Gauss,Cswitch,plotfig,vis)

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



% This code runs a single instance of Wright-Fisher diecious simulations for 1 gRNA
% with 4 alleles: W - wild type; R - functional resistant; N —
% non-functional resistant; D - drive 

%Inputs:
% s — fitness cost of homozygote D/D and N/N in females
% h — dominance coeff/fitness cost of heterozygote W/D in females
% hN — dominance coeff/fitness cost of heterozygote W/N in females
% sigma — fitness cost of functional resistant heterozygote W/R and half 
%   the fitness cost of R/R, before introduction of drive (dominance
%   coeff=1/2) in females
% Rm — absolute or intrinsic growth rate of population assuming population
%   is all W
% K — carrying capacity of population assuming population is all W
% xD — initial frequency of drive in males
% yD — initial frequency of drive in females
% epsilon — efficiency of drive cleavage
% nu — non-homologous end-joining (NHEJ) rate per generation per individual
% mu — mutation rate per cleavage target site per generation per individual
%   (e.g. if target site contains 10bp and bp mutation rate is 1e-9, mu=1e-8)
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

%Fraction of males to females
psi=0.5;

%Identify which indices correspond resistance alleles
test_allele = [2,3]; %R, N are resistant
test_threshold = 0.95;

%Density dependent parameter from Beverton-Holt model of population
%dynamics
alpha = K/(Rm-1);

%Fitness wij of i/j genotype
w12 = 1-sigma;
w13 = 1-hN*s;
w14 = 1-h*s;
w22 = (1-sigma)^2;
w23 = (1-sigma)*(1-hN*s);
w24 = (1-sigma)*(1-hN*s);
w33 = 1-s;
w34 = 1-s;
w44 = 1-s;


%Fitness matrix of females
Wf = [1 w12 w13 w14;...
    w12 w22 w23 w24;...
    w13 w23 w33 w34;...
    w14 w24 w34 w44];
   
%Fitness matrix of males
Wm = ones(4);


%Mutation matrix
Mu = [-mu 0 0 0;...
    xi*mu 0 0 0;...
(1-xi)*mu 0 0 0;...
        0 0 0 0];




%Set up drive heterozygotes

% Stochiometry vector of *changes* in allele frequency due to drive action
% with NHEJ
S = [-epsilon; epsilon*nu*beta; epsilon*nu*(1-beta); epsilon*(1-nu)];

% Identify which alleles correspond to drive heterozygotes in x and y
H = [1;4];




if preexisting==1

    xR = xi*mu/(sigma/2);
    yR = xi*mu/(sigma/2);
    
    
    %N.B. setting initial frequency of drive to zero is sufficient to turn
    %drive off
    x0 = [1-xR;xR;0;0];
    y0 = [1-yR;yR;0;0];
    
    Tpre = round(1/sigma);
    
    %Run simulations for 1/sigma generations *without drive*
    [t,x,y,z,M,F,N,wm,wf,resistance] = TargetSiteResistanceSims_while(4,K,psi,x0,y0,Wm,Wf,Mu,S,H,Rm,alpha,Tpre,Gauss,test_allele, Inf, Cswitch);
    
    xR = x(2,end);
    yR = y(2,end);
    
    
    x0 = [1-xD-xR;xR;0;xD];
    y0 = [1-yD-yR;yR;0;yD];
    
    
    
    if ~isempty(find(x0<0, 1))
        disp('Error: frequency vector must not have negative entries')
        xR
        x(:,end-10:end)
        yR
        y(:,end-10:end)
        xD
        yD
        resistance
%         input('')
    end
 
    if ~isempty(find(y0<0, 1))
        disp('Error: frequency vector must not have negative entries')
        xR
        x(:,end-10:end)
        yR
        y(:,end-10:end)
        xD
        yD
        resistance
%         input('')
    end
    
    
    
    
    
else
    xR = 0;
    yR = 0;
    
    x0 = [1-xD-xR;xR;0;xD];
    y0 = [1-yD-yR;yR;0;yD];
    
end



[t,x,y,z,M,F,N,wm,wf,resistance] = TargetSiteResistanceSims_while(4,K,psi,x0,y0,Wm,Wf,Mu,S,H,Rm,alpha,T,Gauss, test_allele, test_threshold, Cswitch);


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
    
     
    for j=1:3
        plot(t,x(j,:),'Color',col(j,:), 'LineWidth',1,'LineStyle','-');hold on
        plot(t,y(j,:),'Color',col(j,:),'LineWidth',2,'LineStyle','-');hold on
    end
        plot(t,x(4,:),'Color','k', 'LineWidth',1,'LineStyle','-');hold on
        plot(t,y(4,:),'Color','k', 'LineWidth',2,'LineStyle','-');hold on
  
        plot(t(1:end-1),wm(1:end-1),'Color','k','LineWidth',1,'LineStyle','--');
        plot(t(1:end-1),wf(1:end-1),'Color','k','LineWidth',2,'LineStyle','--');

    
    xlabel('Time $t$ (generations)')
    ylabel('Frequency')
    

    hl1=legend('$x_1$ (Male $\mathsf{W}$)','$y_1$ (Female W)','$x_2$ (Male R)','$y_2$ (Female R)','$x_3$ (Male N)','$y_3$ (Female N)',...
            '$x_4$ (Male D)','$y_4$ (Female D)','$\bar{w}_m$','$\bar{w}_f$',...
            'Location','best');

    hl1.FontSize = 14;
        
   title(['$m=1$ gRNA: $N = ',num2strpow(K),'; s = ',num2str(s),'; h = ',num2str(h),'; h_N = ',num2str(hN),'; \sigma = ',num2str(sigma),...
        '; \nu=',num2strpow(nu),'; \mu = ',num2strpow(mu),...
        '; \beta = ',num2strpow(beta),'; \xi = ',num2strpow(xi),'$'])
%     title(['$N = ',num2strpow(K),'$'])
        
    
    set(gca,'FontSize',22)
%     set(gca,'XScale','log')
    set(gca,'YScale','log')
    


    xlim([0 max(t)])
    
    hold off




    subplot(2,1,2);plot(t,M,t,F);hold on
    plot(t,N,'k')
    plot(t,K*ones(size(t)),'k--');hold off
    
    xlabel('Time $t$ (generations)')
    ylabel('Population size')

    hl2 = legend('Male','Female','Total','Intrinsic carrying capacity');
    hl2.FontSize = 14;

    xlim([0 max(t)])
    ylim([0.1 1.2*K])
    
    set(gca,'FontSize',22)
    
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')

    
    if vis==0
        
        savefig(f,['TargetSiteResistance_4Allele_N=',num2strexp(K),'_s=',num2str(s),'_h=',num2str(h),'_hN=',num2str(hN),'_sigma=',num2str(sigma),...
            '_epsilon=',num2str(epsilon),'_nu=',num2str(nu),'_beta=',num2str(beta),'_xi=',num2str(xi),'.fig'])
    
    end


end



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



