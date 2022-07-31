function [problem] = problem_generator(name, dim, C1, C2)

% generate test multi-objective problems
% input: name  -  'DOD2DK' or 'Lam\"e'
%        dim   - dimension of the search space
%        C1,C2 -  auxiliary parameters
% output: structure problem
%       problem.g   - objective fucntion
%       problem.lb  - lower bound of the search space, for each i=1, ..,dim
%       problem.ub  - upper bound of the search space, for each i=1, ..,dim
%       problem.F   - reference approximation of Pareto front
%       problem.Fx  - approximation  of EP optimal points
%       problem.max - strong maximal point used to compute hypervolume
%           contribution



%auxiliary function
select = @(M,i) M(:,i);

problem.dim = dim;
problem.name = name;

filename = 'ReferenceSolutions';
switch name
        
    case 'Lamé'
        gamma = C1;
        problem.g  = @(x) lame(x, gamma,2) +(pi/gamma).*distH(x);
        problem.lb = 0.*ones(1,dim);
        problem.ub = 1.*ones(1,dim);
        
        % data to compute performance metrics
        problem.Fx = zeros(100,dim);
        problem.Fx(:,1) = linspace(0,1,100)';
        switch gamma
            case 0.25
                load(filename ,'lame0_25')
                problem.F = lame0_25;
            case 1
                load(filename ,'lame1')
                problem.F = lame1;
            case 3
                load(filename,'lame3')
                problem.F = lame3;
            otherwise
                % if not pre-computed reference solution is available, it is
                % created by taking an uniform ditribution over a
                % parametrization of the optimal EP points
                problem.F = problem.g(problem.Fx);
        end
        
        problem.max = [1,1];
        
    case 'DOD2DK'
        K=C1;
        s=C2;
               
        problem.g = @(x) do2dk( x, K, s ) + 10.*distH(x);
        problem.lb = 0.*ones(1,dim);
        problem.ub = 1.*ones(1,dim);
        
        %optimal points are known
        problem.Fx = [linspace(0,1,100)' zeros(100,1)  zeros(100,1)];
        if C1 ==2 && C2 ==1
                load(filename,'DO2DK_2_1')
                problem.F = DO2DK_2_1;
        elseif C1 ==4 && C2 ==2
                load(filename,'DO2DK_4_2')
                problem.F = DO2DK_4_2;
        else
            % if not pre-computed reference solution is available, it is
            % created by taking an uniform ditribution over a
            % parametrization of the optimal EP points
            F = problem.g(problem.Fx);
            [F,dompos,~] = faster_pareto2(F);
            problem.F = F;
            problem.Fx = problem.Fx(dompos);        
        end
        
        problem.max = [10,10];

    otherwise
        message = ['The selected problem does not exists, choose between', ...
             'Lam\"e and DOD2DK. Specifiy one or two additional parameters respectively'];
        error(message)
end


problem.g1 = @(x) select(problem.g(x),1);
problem.g2 = @(x) select(problem.g(x),2);

end





function [ f ] = lame( x, gam, nob )
%LAME(x,gam,nob) A scalable version in the number of objectives of the Lame
%hypersphere problem proposed in "Test problems based on Lam� superspheres"
%by Emmerich and Deutz.
%   x is the matrix of the decision variables. Each row represents one
%   solution, while the number of columns determines the dimensionality of
%   each solution (number of variables) Please note that the following
%   inequality must hold: number of variables >= number of objectives - 1
%   gam determines the curvature of the pareto front
%   nob specifies the number of objectives

% number of decision variables
n = size(x,2);

sinx = sin( pi/2 * x(:,1:(nob-1) ));
cosx = cos( pi/2 * x(:,1:(nob-1) ));

f = ones(size(x,1),nob);

% Main loop for calculating the trigonometric functions. Awesomly
% vectorized code that makes the function completely unreadable, yet highly
% performative
for k=1:(nob-1)
    f(:,k) = f(:,k) .* cosx(:,k);
    f(:,(k+1):nob) = f(:,(k+1):nob) .* repmat(sinx(:,k),1,nob-k);
end

% Exponentiate
f = abs(f); %added by G.Borghi
f = f.^(2/gam);

% Multiply by function g
if nob <= n
    g = sqrt(sum(x(:,nob:n).^2,2));
    f = f .* (1 + repmat(g,1,nob));
end

end

function [ f ] = do2dk( x, K, s )
%DO2DK Calculates DO2DK objective values. Pareto optimal solutions are
%given by x=(x1,0,...,0) for s \in [0,1]. For s > 1 you have to filter
%dominated solutions from the solution set.
%   x is a matrix of decision variables. Each row corresponds to a single
%   solution's decision variable values.
%   K is the number of knees
%   s parameter for skewing

[rows, nObj] = size(x);

r = 5 + 10*( x(:,1) - 0.5).^2 + 1/K * cos(2*K*pi*x(:,1)) * 2^(s/2);

g = 1;

% Prevent division by 0 if nObj = 1;
if nObj > 1
    g = g + 9/(nObj -1) * sum(x(:,2:nObj),2);
end

f = zeros(rows,2);

f(:,1) = g .* r .* (sin(pi * x(:,1)/(2^(s+1)) + (1 + (2^s - 1)/(2^(s+2)))*pi )+1);
f(:,2) = g .* r .* (cos(pi * x(:,1)/2 + pi) + 1);
f = abs(f); %added by G. Borghi

end

function [distance] = distH(x)
%x must be a vector Nxd
projx = x;
projx(x>1) = 1;
projx(x<0) = 0;
distance = vecnorm(projx - x,2,2);

end
