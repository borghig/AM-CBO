% The AM-CBO algorithm is tested for 5 different problems in the following
% 4 settings:
%   1. No parameters interaction, tau = 0;
%   2. Riesz potential,           tau = 10−5
%   3. Newtonian potential,       tau = 10−3
%   4. Morse potential,           tau = 10−1,  C = 20;

clear all; close all; clc

% statistics are done over several tests 
Mtest = 2; 

% select problems to solve 
pbnames = {'Lamé', 'Lamé', 'Lamé', 'DOD2DK', 'DOD2DK'}; 
C1s    = [0.25       1       3         2        4  ];
C2s    = [   1       1       1         1        2  ];
dim = 10;

%select potentials
potentials = {'NoInt', 'Riesz','Newtonian','Morse'};
potLenghts = [  1      1        1          20     ];
taus       = [  0      1e-5     1e-3       1e-1   ];


% select algorithm parameters
ev.deltat    = 1e-2;  
ev.alpha     = 1e+6;
ev.lambda    = 1;
ev.sigma     = 3;
ev.tau       = [];          %! to be determined later
ev.diffusion = 'aniso';     %options: 'iso' or 'aniso'
ev.kmax      = 1e+4;
ev.potential = [];          %! to be determined later 
ev.potLenght = [];          %! to be determined later
ev.compute_all_stat = true; %eventually compute additional statistics

% intilize particles in the structure in.-
N = 100;
in.x0 = [];                 %! to be determined later             
in.w0 = linspace(0,1,N)';


%preallocation
ZEROS   = zeros(length(pbnames),length(potentials),ev.kmax);
Ustat   = ZEROS ;
IGDstat = ZEROS ;
MORstat = ZEROS ;
NEWstat = ZEROS ;
RIEstat = ZEROS ;
GDstat  = ZEROS ;
hypstat = zeros(length(pbnames),length(potentials)) ;

bar = waitbar(0,'figure1and2');
hline = '-------------------------------------------------------------------------------\n';

tic
for pi = 1: length(pbnames)
    name = cell2mat(pbnames(pi));
    problem = problem_generator(name, dim, C1s(pi), C2s(pi));   
    
    %display
    fprintf(hline);
    fprintf('Problem: %6s, C1 = %3.2f, C2 = %d \n', name,C1s(pi),C2s(pi));
    fprintf(hline);
    fprintf('             iteraction        GD    Riesz    Newt.   Morse         S      IGD\n');
    for ki = 1:length(potentials)

        % set remaining paramters
        ev.tau = taus(ki);
        ev.potential = cell2mat(potentials(ki));
        ev.potLenght = potLenghts(ki);
        in.x0 = rand(N,problem.dim).*(problem.ub-problem.lb) + problem.lb;
        
        % display 1/2
        fprintf('%10s, tau = %1.0e: ', cell2mat(potentials(ki)),taus(ki));
        
        % run the algorithm
        out = AM_CBO_Mtimes(Mtest, problem, ev, in);
        
        % display 2/2
        fprintf('%3.2e %3.2e %3.2e %3.2e %3.2e %3.2e \n',...
            out.GDstat(end), out.IGDstat(end), out.RIEstat(end),...
            out.NEWstat(end), out.MORstat(end), out.hypstat(end))
        
        % store statistics
        GDstat (pi,ki,:) = out.GDstat;
        IGDstat(pi,ki,:) = out.IGDstat; 
        Ustat  (pi,ki,:) = out.Ustat;
        RIEstat(pi,ki,:) = out.RIEstat; 
        NEWstat(pi,ki,:) = out.NEWstat; 
        MORstat(pi,ki,:) = out.MORstat;
        hypstat(pi,ki)   = out.hypstat;
        
        %store best run data
        bestx(pi,ki) = {out.bestx};
        bestw(pi,ki) = {out.bestw};
        bestg(pi,ki) = {out.bestg};
        refsol(pi,ki) = {problem.F};
        
        % update wait bar 
        waitbar(((pi)*(length(potentials)-1) + ki)/(length(pbnames)*length(potentials)),bar);
    end
end
fprintf('\n')
toc

close(bar);

% store data for plotting
filename = sprintf('figure1and2_data.mat');
save(filename);

%% visualization
close all

%decide whether to save the generated figures
savefigures = true;

% load computed data
load('figure1and2_data.mat');

% set some graphics settings
set(0, 'DefaultLineLineWidth', 1.1);
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

colors = [0, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125;	          
                              0.494, 0.184, 0.556;  0.466, 0.674, 0.188];                  
Kmax = ev.kmax;
dt = 2;

%%%%%%%%%%%%%%%%%%%%%%
% PLOT best fronts
figure;
tl1  = tiledlayout(length(pbnames), length(potentials));

% figures iteration
for pi = 1:length(pbnames)
    
    for ki = 1:length(potentials)
        nexttile;
  
        F = cell2mat(refsol(pi,ki));
        g = cell2mat(bestg(pi,ki));
        
        % plot 
        scatter(F(:,1),F(:,2),'.','MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:)); 
        
        axis padded
        axis manual
        
        hold on;
        ll = scatter(g(:,1), g(:,2),15,'filled','k','MarkerFaceAlpha',.3, ...
                'MarkerEdgeAlpha',.3);   
                  
        % add information on problem
        pbaspect([1 1 1]);
        if pi ==1 
            title({cell2mat(potentials(ki)),''},'Interpreter','latex');
        end   
        if ki ==1
            problemname = sprintf('%s, C1=%g', cell2mat(pbnames(pi)),C1s(pi));
            ylabel({problemname,''});
        end    
    end
end

tl1.TileSpacing = 'tight';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0 0.5 1]);
title(tl1,'Computed fronts','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%
% PLOT best fronts' histograms 

figure;
tl2  = tiledlayout(length(pbnames), length(potentials));

nbins = 11;
edges = linspace(0,1,nbins);
trasp = 0.2;

%figures iteration
for pi = 1:length(pbnames)
    
    for ki = 1:length(potentials)
        
        nexttile;
  
        F = cell2mat(refsol(pi,ki));
        w = cell2mat(bestw(pi,ki));
        
        normF = vecnorm(F,1,2);
        wexact = F./normF;
        wexact1 = wexact(:,1); %if w1 =0. point down-right on front
        
        histogram(wexact1,edges,'Normalization','probability','FaceColor','r',...
            'FaceAlpha',trasp, 'DisplayName', 'exact');
        hold on;
        
        axis manual
        histogram(1-w,edges,'Normalization','probability','FaceColor','b',...
            'FaceAlpha',trasp,'DisplayName', 'computed')
        pbaspect([1 1 1])
        set(gca,'XTick',[], 'YTick', [])
        
                          
        % add information on problem
        pbaspect([1 1 1]);
        if pi ==1 
            title({cell2mat(potentials(ki)),''},'Interpreter','latex');
        end   
        if ki ==1
            problemname = sprintf('%s, C1=%g', cell2mat(pbnames(pi)),C1s(pi));
            ylabel({problemname,''});
        end    
    end
end

title(tl2,'Parameters distribution over $\Omega$','Interpreter','latex')
tl2.TileSpacing = 'compact';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0 0.5 1]);

legg = legend(gca,'Orientation', 'Horizontal');
legg.Layout.Tile = 'north';

        
%%%%%%%%%%%%%%%%%%%
% plot the metrics evolution

figure;
tl3 = tiledlayout(length(pbnames),5);

% figures iteration
for pi = 1:length(pbnames)
    for ki = 1: length(potentials)    
        leg = cell2mat(potentials(ki));
        col = colors(ki,:);
        
        %plot 1
        nexttile((pi-1)*5 +1 );
        semilogy(1:dt:Kmax,squeeze(GDstat(pi,ki,1:dt:Kmax)),'Color',col,'DisplayName',leg);
        title('GD')
        xlabel('$k$', 'Interpreter','latex')
        %legend; 
        hold on; %axis padded
        problemname = sprintf('%6s, C1=%g', cell2mat(pbnames(pi)),C1s(pi));
        ylabel(problemname)
        
        %plot 2
        nexttile((pi-1)*5 +2);  
        semilogy(1:dt:Kmax,squeeze(RIEstat(pi,ki,1:dt:Kmax)),'Color',col,'DisplayName',leg);
        title('Riesz en.')
        xlabel('$k$', 'Interpreter','latex')
        %legend('location','southeast'); 
        hold on; %axis padded
        
        %plot 3
        nexttile((pi-1)*5 +3);  
        plot(1:dt:Kmax,squeeze(NEWstat(pi,ki,1:dt:Kmax)),'Color',col,'DisplayName',leg);
        title('Newtonian en.')
        xlabel('$k$', 'Interpreter','latex')
        %legend('location','southeast'); 
        hold on; %axis padded
        
        %plot 3
        nexttile((pi-1)*5 +4);  
        plot(1:dt:Kmax, squeeze(MORstat(pi,ki,1:dt:Kmax)),'Color',col,'DisplayName',leg);
        title('Morse en.')
        xlabel('$k$', 'Interpreter','latex')
        %legend('location','southeast'); 
        hold on; %axis padded
        
        %plot 5
        nexttile((pi-1)*5 +5); 
        semilogy(1:dt:Kmax,squeeze(IGDstat(pi,ki,1:dt:Kmax)),'Color',col,'DisplayName',leg); 
        title('IGD')
        xlabel('$k$', 'Interpreter','latex')
        %legend;
        hold on; %axis padded
  
    end
end

tl3.TileSpacing = 'tight';
legg = legend(gca,'Orientation', 'Horizontal');
legg.Layout.Tile = 'north';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
title(tl3,'Performance metrics evolution','Interpreter','latex')

% save figure
if savefigures
    filename = sprintf('figure1_fronts.pdf');
    exportgraphics(tl1,filename)
    
    filename = sprintf('figure1_hist.pdf');
    exportgraphics(tl2,filename)
    
    filename = sprintf('figure2.pdf');
    exportgraphics(tl3,filename)
end

    
        
       



