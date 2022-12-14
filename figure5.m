% The AM-CBO algorithm is tested for 5 different problems
% for different values of the search space d
% for one proposed energy potential (Morse) 


clear all; close all; clc

% statistics are done over several tests 
Mtest = 5; 

% select problems to solve 
pbnames = {'Lamé', 'Lamé', 'Lamé', 'DOD2DK', 'DOD2DK'}; 
C1s    = [0.25       1       3         2        4  ];
C2s    = [   1       1       1         1        2  ];
dims   = [ 5:5:30 ];

%select potentials
potentials = {'Riesz', 'Newtonian', 'Morse'};
potLenghts = [ 1        1          20     ];
taus       = [  0      1e-5     1e-3       1e-1   ];

% select potential to use
potindex = 3;

% select algorithm parameters
ev.deltat    = 1e-2;  
ev.alpha     = 1e+6;
ev.lambda    = 1;
ev.sigma     = 3;
ev.tau       = taus(potindex);           %! to be determined later
ev.diffusion = 'aniso';      %options: 'iso' or 'aniso'
ev.kmax      = 5e+3;
ev.potential = cell2mat(potentials(potindex));      
ev.potLenght = potLenghts(potindex);           
ev.compute_all_stat = false; %eventually compute additional statistics

% intilize particles in the structure in.-
N  = [];                    %! to be determined later
Nfixed = 100;               % plot left: N fixed to 100
Ncoeff = 10;                % plot right: N = Ncoeff*d
in.x0 = [];                 %! to be determined later             
in.w0 = [];                 %! to be determined later


Nprob = length(pbnames);
Ndim  = length(dims);

%preallocation
ZEROS   = zeros(Nprob,Ndim,2);
Ustat   = ZEROS ;
IGDstat = ZEROS ;
MORstat = ZEROS ;
NEWstat = ZEROS ;
RIEstat = ZEROS ;
GDstat  = ZEROS ;
hypstat = ZEROS ;

bar = waitbar(0,'figure5');
hline = '-------------------------------------------------------------------------------\n';

tic
for pi = 1: Nprob 
    name = cell2mat(pbnames(pi));
    
    %display
    fprintf(hline);
    fprintf('Problem: %6s, C1 = %3.2f, C2 = %d  (%d/%d) \n d = ', name,C1s(pi),C2s(pi),pi,Nprob);
    
    for di = 1:Ndim

        % select dimension and problem
        dim = dims(di);
        problem = problem_generator(name, dim, C1s(pi), C2s(pi));  
        
        % display
        fprintf('%d, ', dim);
        
        % (1/2) compute with fixed N = Nfixed
        N = Nfixed;
        in.w0 = linspace(0,1,N)';
                    
        % run the algorithm
        out = AM_CBO_Mtimes(Mtest, problem, ev, in);

        % store statistics
        GDstat (pi,di,1) = out.GDstat(end);
        IGDstat(pi,di,1) = out.IGDstat(end);
        Ustat  (pi,di,1) = out.Ustat(end);
        RIEstat(pi,di,1) = out.RIEstat(end);
        NEWstat(pi,di,1) = out.NEWstat(end);
        MORstat(pi,di,1) = out.MORstat(end);
        hypstat(pi,di,1) = out.hypstat(end);
        
        % (2/2) compute with increasing N = Ncoeff*dim
        N = Ncoeff*dim;
        in.w0 = linspace(0,1,N)';
                    
        % run the algorithm
        out = AM_CBO_Mtimes(Mtest, problem, ev, in);

        % store statistics
        GDstat (pi,di,2) = out.GDstat(end);
        IGDstat(pi,di,2) = out.IGDstat(end);
        Ustat  (pi,di,2) = out.Ustat(end);
        RIEstat(pi,di,2) = out.RIEstat(end);
        NEWstat(pi,di,2) = out.NEWstat(end);
        MORstat(pi,di,2) = out.MORstat(end);
        hypstat(pi,di,2) = out.hypstat(end);
        

        % update wait bar
        waitbar(((pi-1)*(Ndim-1) + di)/(Nprob*Ndim),bar);
        
    end
    fprintf('\n')
end
fprintf('\n')
toc

close(bar);

% store data for plotting
filename = sprintf('figure5_data.mat');
save(filename);


%% visualization
close all

%decide whether to save the generated figures
savefigures = true;

% load computed data
load('figure5_data.mat');

% set some graphics settings
set(0, 'DefaultLineLineWidth', 1.1);
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

colors = [0, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125;	          
                              0.494, 0.184, 0.556;  0.466, 0.674, 0.188]; 
mark = 'o*xs^';


%%%%%%%%%%%%%%%%%%%%%%
% PLOT GD and IGD
figure;
tl1  = tiledlayout(1, 2);

% figures iteration
for strategy = [1 2]
    nexttile

    switch strategy
        case 1
            stratname = sprintf('$N = %d$',Nfixed);
        case 2
            stratname = sprintf('$N = %d d$',Ncoeff);
    end
    
    for pi = 1:Nprob
        % plot
        lgnd = sprintf('%s,  C1 = %g',cell2mat(pbnames(pi)), C1s(pi));
        semilogy(dims,squeeze(IGDstat(pi,:,strategy)),'-','Marker',mark(pi),...
            'MarkerSize',5,'Color',colors(pi,:),'DisplayName',lgnd);
        hold on
    end

    if strategy ==1
        ylabel('IGD','Interpreter','latex');
    end

    title(stratname,'Interpreter','latex');
    xlabel('d','Interpreter','latex')

    axis padded
    ylim([9e-3,1])
    axis manual
    
end

leg = legend(gca,'Orientation', 'Horizontal','Interpreter','tex');
leg.Layout.Tile = 'north';    
tl1.TileSpacing = 'compact';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.2 0.6 0.6]);
title(tl1,'IGD metric as a function of dimension $d$','Interpreter','latex')

% save figure
if savefigures
    filename = sprintf('figure5.pdf');
    exportgraphics(tl1,filename)
end

    