% The AM-CBO algorithm is tested for 5 different problems
% for different values of the parameter sigma and taus
% for one proposed energy potential (Morse)


clear all; close all; clc

% statistics are done over several tests 
Mtest = 5; 

% select problems to solve 
pbnames = {'Lamé', 'Lamé', 'Lamé', 'DOD2DK', 'DOD2DK'}; 
C1s    = [0.25       1       3         2        4  ];
C2s    = [   1       1       1         1        2  ];
dim = 10;

%select potentials
potentials = {'Riesz', 'Newtonian', 'Morse'};
potLenghts = [ 1        1          20     ];
taus       = [ 0  1e-2  1e-1  1];
sigmas     = [0:5 6:2:12];

% select algorithm parameters
ev.deltat    = 1e-2;  
ev.alpha     = 1e+6;
ev.lambda    = 1;
ev.sigma     = [];           %! to be determined later
ev.tau       = [];           %! to be determined later
ev.diffusion = 'aniso';      %options: 'iso' or 'aniso'
ev.kmax      = 5e+3;
ev.potential = cell2mat(potentials(3));  %restrict to Morse potential 
ev.potLenght = potLenghts(3);           %! to be determined later
ev.compute_all_stat = false; %eventually compute additional statistics

% intilize particles in the structure in.-
N = 100;
in.x0 = [];                 %! to be determined later             
in.w0 = linspace(0,1,N)';


Nprob = length(pbnames);
Ntau  = length(taus(1,:));
Nsigma  = length(sigmas);

%preallocation
ZEROS   = zeros(Nprob,Ntau,Nsigma);
Ustat   = ZEROS ;
IGDstat = ZEROS ;
MORstat = ZEROS ;
NEWstat = ZEROS ;
RIEstat = ZEROS ;
GDstat  = ZEROS ;
hypstat = ZEROS ;

bar = waitbar(0,'figure4');
hline = '-------------------------------------------------------------------------------\n';

tic
for pi = 1: Nprob 
    name = cell2mat(pbnames(pi));
    problem = problem_generator(name, dim, C1s(pi), C2s(pi));   
    
    %display
    fprintf(hline);
    fprintf('Problem: %6s, C1 = %3.2f, C2 = %d  (%d/%d) \n tau = ', name,C1s(pi),C2s(pi),pi,Nprob);
    %fprintf(hline);
    
    for ti = 1:Ntau
        % set tau 
        ev.tau = taus(ti);
        fprintf('%1.0e, ', ev.tau);
        
        for si = 1:Nsigma
            
            %set sigma
            ev.sigma = sigmas(si);
            
            % run the algorithm
            out = AM_CBO_Mtimes(Mtest, problem, ev, in);
            
            % store statistics
            GDstat (pi,ti,si) = out.GDstat(end);
            IGDstat(pi,ti,si) = out.IGDstat(end);
            Ustat  (pi,ti,si) = out.Ustat(end);
            RIEstat(pi,ti,si) = out.RIEstat(end);
            NEWstat(pi,ti,si) = out.NEWstat(end);
            MORstat(pi,ti,si) = out.MORstat(end);
            hypstat(pi,ti,si) = out.hypstat(end);
            
            % update wait bar
            waitbar(((pi-1)*(Nsigma-1)*(Ntau-1) + ti*(Nsigma-1) + si)/(Nprob*Nsigma*Ntau),bar);
        end
    end
    fprintf('\n')
end
fprintf('\n')
toc

close(bar);

% store data for plotting
filename = sprintf('figure4_data.mat');
save(filename);


%% visualization
close all

%decide whether to save the generated figures
savefigures = true;

% load computed data
load('figure4_data.mat');

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
tl1  = tiledlayout(2, Nprob);

% figures iteration
for measure = [1 2]

    switch measure
        case 1
            data = GDstat;
            measname = 'GD';
        case 2
            data = IGDstat;
            measname = 'IGD';
    end
    
    for pi = 1:Nprob
        nexttile;
        for ti = 1:Ntau
            % plot
            legnd = sprintf('tau = %1.0g',taus(ti));
            semilogy(sigmas,squeeze(data(pi,ti,:)),'-','Marker',mark(pi),...
                'MarkerSize',5,'Color',colors(ti,:),'DisplayName',['$\',legnd,'$']);
            hold on
        end
        
        if measure ==1
            titlename = sprintf('%s,  C1 = %g',cell2mat(pbnames(pi)), C1s(pi));
            title({titlename,''},'Interpreter','tex');
        end
        
        if pi ==1
            ylabel(measname);
        end
        xlabel('$\sigma$', 'Interpreter','latex')
        
        axis padded
        axis manual
    end
end

leg = legend(gca,'Orientation', 'Horizontal','Interpreter','latex');
leg.Layout.Tile = 'north';    
tl1.TileSpacing = 'compact';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.2 0.6 0.6]);
title(tl1,{'Performance metrics as a function of $\sigma$ for different $\tau$', ''},'Interpreter','latex')

% save figure
if savefigures
    filename = sprintf('figure4.pdf');
    exportgraphics(tl1,filename)
end

    
       