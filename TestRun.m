% Run test

clear all; close all; clc

% select problem parameters ('Lamé' or 'DOD2DK')
name = 'DOD2DK';
dim = 10;           %dimension of the serach space
C1 = 2;              % additional parameter (gamma or k)
C2 = 1;              %additional paramer (s, only for 'DOD2DK')
problem = problem_generator(name, dim, C1, C2);

% set algorithm parameters in the structure ev.-
ev.deltat    = 1e-2;  
ev.alpha     = 1e+6;
ev.lambda    = 1;
ev.sigma     = 3;
ev.tau       = 1e-5;
ev.diffusion = 'aniso';     %options: 'iso' or 'aniso'
ev.kmax      = 5e+3;
ev.potential = 'Morse';     %options: 'Riesz', 'Newtonian', 'Morse' 
ev.potLenght = 20;          %addtional parameter for 'Morse' potential


% intilize particles in the structure in.-
N = 100;
in.x0 = rand(N,problem.dim).*(problem.ub-problem.lb) + problem.lb;
in.w0 = linspace(0,1,N)';


output = AM_CBO(problem, ev, in);

%% plot results
close all;
figure; 
set(gcf,'units','centimeters','position',[15,12,14,7])


%FIGURE 1: front and parameter distribution
%plot computed front
tl1 = tiledlayout(1,2);


nexttile;
scatter(problem.F(:,1),problem.F(:,2),'.r','DisplayName','$F$'); 
hold on
pbaspect([1 1 1])
axis manual

g1 = output.g(:,1); g2 = output.g(:,2);
figf = scatter(g1,g2,'filled','k','MarkerFaceAlpha',.2, ...
    'MarkerEdgeAlpha',.2,'DisplayName','$g(X_i^k)$');

xlabel('$g_1$','interpreter','latex');
ylabel('$g_2$','interpreter','latex');

legend('Interpreter','latex')
title('computed front')

%plot final parameter distribution
nexttile;

nbins = 11;
edges = linspace(0,1,nbins);
trasp = 0.2;

%plot exact distribution
normF = vecnorm(problem.F,1,2);
wexact = problem.F./normF;
wexact1 = wexact(:,1);

histogram(wexact1,edges,'Normalization','probability','FaceColor','r',...
    'FaceAlpha',trasp,'DisplayName', 'optimal');
hold on;
pbaspect([1 1 1])
histogram(1-output.w,edges,'Normalization','probability','FaceColor','b',...
    'FaceAlpha',trasp,'DisplayName', 'computed')
xlabel('$\Omega$','Interpreter','latex')
legend
title('parameters distribution')


% FIGURE 2: evolution of performane metrics
figure;
set(gcf,'units','centimeters','position',[10,2,22,7])


tl2  = tiledlayout(1,3);

nexttile;
semilogy(output.GDtot, 'b');
xlabel('k','interpreter','latex')
title('GD')

nexttile;
plot(output.Utot, 'g');
xlabel('k','interpreter','latex')
title([ev.potential,' potential'])

nexttile;
semilogy(output.IGDtot, 'r');
xlabel('k','interpreter','latex')
title('IGD')




