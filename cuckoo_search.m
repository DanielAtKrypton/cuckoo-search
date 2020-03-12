% Cuckoo Search (CS) algorithm by Daniel Kaminski de Souza
% Based on XinShe Yang and Suash Deb algoritm
%
% Programmed by: Daniel Kaminski de Souza at:
%   Universidade Federal do Paraná
%
% Oriented by:
%   Gideon Villar
%   Gustavo Oliveira
%
% Programming dates: Nov 2012
% Last revised: Nov 2012   (simplified version for demo only)
%
% PapersCitation Details:

% 1) X.S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspired
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA,  pp. 210214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf
%
% 2) X.S. Yang, S. Deb, Engineering optimization by cuckoo search,
% Int. J. Mathematical Modelling and Numerical Optimisation,
% Vol. 1, No. 4, 330343 (2010).
% http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
%
% 3) R. N. Mantegna, Fast, accurate algorithm for numerical ...simulation of
% Levy stable stochastic processes. Dipartimento di Energetica ed ...Applica
% zioni di Fisica, Università degli Studi, viale delle Scienze, I90128
% Palermo, Italy (Received 28 October 1993)
%
% This demo program only implements a standard version of
% Cuckoo Search (CS), as the Levy flights and generation of
% new solutions may use slightly different methods.
% The pseudo code was given sequentially (select a cuckoo etc),
% but the implementation here uses Matlab's vector capability,
% which results in neater/better codes and shorter running time.
% This implementation is different and more efficient than the
% the demo code provided in the book by
%    "Yang X. S., NatureInspired Metaheuristic Algoirthms,
%     2nd Edition, Luniver Press, (2010).                 "
%
% ===============================================================
% Notes:
% Different implementations may lead to slightly different
% behavour and/or results, but there is nothing wrong with it,
% as this is the nature of random walks and all metaheuristics.
% Presentation Data Selection:
% 0 = None
% 1 = Numbers
% 2 = Graphic Evolution

function [bestnest,fmin, N_iter]=cuckoo_search(fun, nvars, yd, ...
    varargin)
%CUCKOO_SEARCH    Optimization using Cuckoo Search.

% initialTetaT, lb, ub, IntCon, options
initialTetaT = [];
lb = [];
ub = [];
IntCon =[];
options =[];

minArgs=3;
maxArgs=8;
narginchk(minArgs,maxArgs)

if (length(varargin) >= 1)
    initialTetaT = varargin{1};
end
if (length(varargin) >= 2)
    lb = varargin{2};
end
if (length(varargin) >= 3)
    ub = varargin{3};
end
if (length(varargin) >= 4)
    IntCon = varargin{4};
end
if (length(varargin) >= 5)
    options = varargin{5};
end

% Load default values if needed
if isempty(options)
    %% Estimation Process
    StallGenLimit_Data = 50;
    PopulationSize_Data = 15;
    TolFun_Data = 1e-3;
    DiscoveryRate_Data = 0.25;
    
    % Presentation Data Selection:
    % 0 = None
    % 1 = Numbers
    % 2 = Graphic Evolution
    Presentation_Data = 2;
    
    % LevyFlight Data Selection:
    % 0 = Original
    % 1 = Mantegna
    LevyFlight_Data = 0;
    
    % StepSize Data Selection:
    % 0 = Original randn
    % 1 = Fast randn versor
    StepSize_Data = 0;
    
    options = struct( ...
        'StallGenLimit', StallGenLimit_Data, ...
        'TolFun', TolFun_Data, ...
        'PopulationSize', PopulationSize_Data, ...
        'DiscoveryRate', DiscoveryRate_Data,...
        'PresentationType', Presentation_Data,...
        'LevyFlight', LevyFlight_Data,...
        'StepSize', StepSize_Data...
        );
end

if isempty(lb)
    lb = realmin.*ones(1, nvars);
end
if isempty(ub)
    ub = realmax.*ones(1, nvars);
end

% Discovery rate of alien eggs/solutions
pa=options.DiscoveryRate;

%% Change this if you want to get better results
% Tolerance
StallGenLimit = options.StallGenLimit;

nd = nvars;

% initial solutions
%nest=lb+(ub-lb).*rand(options.PopulationSize, yd, nd);
A = rand(options.PopulationSize, yd, nd);
for i=1:nd
    A(:,:,i) = lb(i)+(ub(i)-lb(i))*A(:,:,i);
end
nest = reshapeNewNest(size(A), lb, ub, A, IntCon);
clear A

clones = nest;
if (isempty(initialTetaT)~=1)
    clones(1,:,:) = initialTetaT;
end

% Get the current best
fitness=Inf*ones(options.PopulationSize,1);
[fmin,bestnest,nest,fitness]=get_best_nest(fun, nest,clones,fitness);
populationSize = size(nest,1);

a = Inf;
b = 0;
i = 0;
switch options.PresentationType
    case 1
        fprintf('Generation\tFunc-count\tBest Penalty\tChange Mean\n');
    case 2
        a = [a fmin];
        b = [b i];
        graphicEvolutionFigure = figure('NumberTitle', 'off', 'Name', 'Cuckoo Search');
        h = semilogy(1,0);
        hold;
        h.Marker = 'o';
        h.MarkerEdgeColor = 'm';
        titleStr = 'Cuckoo estimation evolution';
        title(titleStr)
        xlabel('iteration')
        ylabel('fmin')
        set(h,'XData',b)
        set(h,'YData',a)
end

N_iter=0;
change = Inf*ones(StallGenLimit, 1);
changeMean = Inf;
%% Starting iterations
while (changeMean > options.TolFun)
    % Generate new solutions (but keep the current best)
    new_nest=get_cuckoos(nest,bestnest,lb,ub,IntCon, options);
    [~,~,nest,fitness]=get_best_nest(fun, nest,new_nest,fitness);
    % Update the counter
    N_iter=N_iter+options.PopulationSize;
    % Discovery and randomization
    new_nest=empty_nests(nest,lb,ub,IntCon,pa) ;
    % Evaluate this set of solutions
    [fnew,best,nest,fitness]=get_best_nest(fun, nest,new_nest,fitness);
    % Update the counter again
    N_iter=N_iter+options.PopulationSize;
    % Find the best objective so far
    in = i+1;
    if fnew<fmin   
        a = [a fnew];
        b = [b in];
        currentChange = (a(end-1)-a(end))/(b(end)-b(end-1));
        bestnest=best;
        changeMeanStr = num2str(changeMean);
        switch options.PresentationType
            case 1
%                 disp(changeMeanStr)
                funcCount = in*populationSize;
                bestPenalty = fnew;
                fprintf('\t%d\t\t\t%d\t\t%f\t\t\t%f\n',in, funcCount, bestPenalty, changeMean);
            case 2  % Check if Graphic Evolution selected
                figure(graphicEvolutionFigure);
                title([titleStr ', ' changeMeanStr])
                if length(a) > 60
                    ylim([a(end) a(end-60)]);
                end
                set(h,'XData',b);
                set(h,'YData',a);
        end
        fmin=fnew;
    else
        currentChange = 0;
    end
    if (options.PresentationType == 2)
        xlim([0 in]);
        refreshdata
        drawnow
    end
    index  = (rem(i,StallGenLimit)+1);
    change(index) = currentChange;
    changeMean = mean(change);
    i = i + 1;
end %% End of iterations
%% Postoptimization processing
%% Display all the nests
%disp(strcat('Total number of iterations=',num2str(N_iter)));
%fmin
%bestnest
%change
%mean(change)



%%All subfunctions are list below
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,lb,ub,IntCon, options)
% Levy flights
sizeNest = size(nest);
n=sizeNest(1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, NatureInspired Metaheuristic Algorithms, 2nd Edition, ...Luniver Press, (2010).
reshapeSize = sizeNest(2:end);
levyFlightOption = options.LevyFlight;
stepSizeOption = options.StepSize;
beta=3/2;
parfor j=1:n
    s=reshape(nest(j,:,:), reshapeSize);
    step = [];
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    switch levyFlightOption
        case 0
            sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
            u=randn(size(s))*sigma;
            v=randn(size(s));
            step=u./abs(v).^(1/beta);
        otherwise
            step = mantegna(beta, 1, 1, s);
    end
    
    % In the next equation, the difference factor (s-best) means that
    % when the solution is the best solution, it remains unchanged.
    switch stepSizeOption
        case 0
            stepsize=0.01*step.*(s-best);
            % Here the factor 0.01 comes from the fact that L/100 should ...the typical
            % step size of walks/flights where L is the typical lenghtscale;
            % otherwise, Levy flights may become too aggresive/efficient,
            % which makes new solutions (even) jump out side of the ...design domain
            % (and thus wasting evaluations).
            % Now the actual random walks or flights
            s=s+stepsize.*randn(size(s));
        case 1
            if (s ~= best)
                % stepsize = step;
                stepsize=0.01*step.*(s-best);
                % Now the actual random walks or flights
                % From ...http://stackoverflow.com/questions/9750908/howtogenerateaunitvectorpointinginarandomdirectionwithisotropicdist
                v = randn(size(s));
                randomDistributedVersor = bsxfun(@rdivide,v,sqrt(sum(v.^2)));
                    
            s=s+stepsize.*randomDistributedVersor;
            end
    end
    % Apply simple bounds/limits
    bounded = simplebounds(s,lb,ub);
    nest(j,:,:)=someIntegerParameters(bounded, IntCon);
end

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(fun, nest,newnest,fitness)
    % Evaluating all new solutions
sizeNest = size(nest);
reshapeSize = sizeNest(2:end);
parfor j=1:sizeNest(1)
    currentNest = reshape(newnest(j,:,:), reshapeSize);
    fnew=fun(currentNest);
    if fnew<fitness(j)
        fitness(j)=fnew;
        nest(j,:,:)=currentNest;
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=reshape(nest(K,:,:), reshapeSize);

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,lb,ub,IntCon,pa)
% A fraction of worse nests are discovered with a probability pa
sizeNest = size(nest);
n=sizeNest(1);
% Discovered or nota status vector
K=rand(sizeNest)>pa;

% In the real world, if a cuckoo's egg is very similar to a host's ...eggs, then
% this cuckoo's egg is less likely to be discovered, thus the ...fitness should
% be related to the difference in solutions.  Therefore, it is a ...good idea
% to do a random walk in a biased way with some random step sizes.
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:,:)-nest(randperm(n),:,:));
new_nest=nest+stepsize.*K;
new_nest = reshapeNewNest(sizeNest, lb, ub, new_nest, IntCon);

function new_nest = reshapeNewNest(sizeNest, lb, ub, new_nest, IntCon)
n=sizeNest(1);
signalsNumber=sizeNest(2);
LB = reshape(repmat(lb, n*signalsNumber,1), n, signalsNumber, length(lb));
UB = reshape(repmat(ub, n*signalsNumber,1), n, signalsNumber, length(ub));
new_nest = someIntegerParametersNests(simplebounds(new_nest,LB,UB), IntCon);


% Application of simple constraints
function s=simplebounds(s,Lb_m,Ub_m)
% Apply the lower bound
ns_tmp=s;

try 
    I=ns_tmp<Lb_m;
    ns_tmp(I)=Lb_m(I);

    % Apply the upper bounds
    J=ns_tmp>Ub_m;
    ns_tmp(J)=Ub_m(J);
catch error
    disp(error);
    keyboard;
end

% Update this new move
s=ns_tmp;

function s=someIntegerParameters(s, IntCon)
% Convert to integer variables the ones specified in IntCon. 
for i=1:length(IntCon)
    s(IntCon(i)) = round(s(IntCon(i)));
end

function nests=someIntegerParametersNests(nests, IntCon)
% Convert to integer variables the ones specified in IntCon. 
for i=1:length(IntCon)
    nests(:,:,IntCon(i)) = round(nests(:,:,IntCon(i)));
end