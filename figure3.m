% The AM-CBO algorithm is tested for 5 different problems
% for different values of the parameter tau
% for all proposed energy potential (Riesz, Newtonian,Morse)


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
taus       = [  10.^(-7:-1);
                10.^(-6:-0);
                10.^(-5: 1) ];

% select algorithm parameters
ev.deltat    = 1e-2;  
ev.alpha     = 1e+6;
ev.lambda    = 1;
ev.sigma     = 3;
ev.tau       = [];           %! to be determined later
ev.diffusion = 'aniso';      %options: 'iso' or 'aniso'
ev.kmax      = 5e+3;
ev.potential = [];           %! to be determined later 
ev.potLenght = [];           %! to be determined later
ev.compute_all_stat = false; %eventually compute additional statistics

% intilize particles in the structure in.-
N = 100;
in.x0 = [];                 %! to be determined later             
in.w0 = linspace(0,1,N)';


Nprob = length(pbnames);
Npot  = length(potentials);
Ntau  = length(taus(1,:));

%preallocation
ZEROS   = zeros(Nprob,Npot,Ntau);
Ustat   = ZEROS ;
IGDstat = ZEROS ;
MORstat = ZEROS ;
NEWstat = ZEROS ;
RIEstat = ZEROS ;
GDstat  = ZEROS ;
hypstat = ZEROS ;

bar = waitbar(0,'figure3');
hline = '-------------------------------------------------------------------------------\n';

tic
for pi = 1: Nprob 
    name = cell2mat(pbnames(pi));
    problem = problem_generator(name, dim, C1s(pi), C2s(pi));   
    
    %display
    fprintf(hline);
    fprintf('Problem: %6s, C1 = %3.2f, C2 = %d  (%d/%d) \n', name,C1s(pi),C2s(pi),pi,Nprob);
    %fprintf(hline);
    
    for ki = 1:Npot

        % set potential paramters
        ev.potential = cell2mat(potentials(ki));
        ev.potLenght = potLenghts(ki);
        
        fprintf('%9s: ', cell2mat(potentials(ki)));
        
        for ti = 1:Ntau
            
            %set tau and initilization
            ev.tau = taus(ki,ti);
            
            fprintf('%1.0e, ', ev.tau);
            
            % run the algorithm
            out = AM_CBO_Mtimes(Mtest, problem, ev, in);
            
            % store statistics
            GDstat (pi,ki,ti) = out.GDstat(end);
            IGDstat(pi,ki,ti) = out.IGDstat(end);
            Ustat  (pi,ki,ti) = out.Ustat(end);
            RIEstat(pi,ki,ti) = out.RIEstat(end);
            NEWstat(pi,ki,ti) = out.NEWstat(end);
            MORstat(pi,ki,ti) = out.MORstat(end);
            hypstat(pi,ki,ti) = out.hypstat(end);
            
            % update wait bar
            waitbar(((pi-1)*(Npot-1)*(Ntau-1) + ki*(Ntau-1) + ti)/(Nprob*Npot*Ntau),bar);
        end
        fprintf('\n')
    end
end
fprintf('\n')
toc

close(bar);

% store data for plotting
filename = sprintf('figure3_data.mat');
save(filename);


%% visualization
close all

%decide whether to save the generated figures
savefigures = true;

% load computed data
load('figure3_data.mat');

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
tl1  = tiledlayout(2, Npot);

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
    
    for ki = 1:Npot
        nexttile;
        for pi = 1:Nprob
            % plot
            lgnd = sprintf('%s,  C1 = %g',cell2mat(pbnames(pi)), C1s(pi));
            loglog(taus(ki,:),squeeze(data(pi,ki,:)),'-','Marker',mark(pi),...
                'MarkerSize',5,'Color',colors(pi,:),'DisplayName',lgnd);
            hold on
        end
        
        if measure ==1
            title({cell2mat(potentials(ki)),''},'Interpreter','latex');
        end
        
        if ki ==1
            ylabel(measname);
        end
        
        axis padded
        axis manual
    end
end

leg = legend(gca,'Orientation', 'Horizontal','Interpreter','tex');
leg.Layout.Tile = 'north';    
tl1.TileSpacing = 'compact';
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.2 0.6 0.6]);
title(tl1,'Performance metrics as a function of $\tau$','Interpreter','latex')

% save figure
if savefigures
    filename = sprintf('figure3.pdf');
    exportgraphics(tl1,filename)
end

    
       


%% old
close all

set(0, 'DefaultLineLineWidth', 1.1);
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
colors = [0, 0.4470, 0.7410;	       
          0.8500, 0.3250, 0.0980;	          
          0.9290, 0.6940, 0.1250;	          
          0.4940, 0.1840, 0.5560;	          
          0.4660, 0.6740, 0.1880];


d =10;
filename = sprintf('ParamTAU_%u.mat',d);
%filename = sprintf('ParamTAU_all_%u.mat',d);

load(filename);
names = {'Riesz en.', 'Newtonian en.', 'Morse en.'};



TT = [ logspace(-6,1,Ntau);
       logspace(-5,2,Ntau);
       logspace(-3,3,Ntau)].*ev.eps; 
   
LAB = {{'0','$10^{-6}$','$10^{-4}$', '$10^{-2}$' };   
    {'0','$10^{-5}$','$10^{-3}$', '$10^{-1}$' };
    {'0','$10^{-4}$','$10^{-2}$', '$10^{0}$' }};
   

PBnames = {'Lam\''e ${0.25}$', 'Lam\''e ${1}$','Lam\''e ${3}$','DO2DK ${2,1}$','DO2DK ${4,2}$'};
mark = 'o*xs^';

yti = [10.^([-2 0 2]);  10.^([-2,-1,0])];


for fig = 1:2 
    figure(fig)
    switch fig
        case 1 %GD metric
            DATA = GDstat;
            titlename = 'GD';
        case 2 %IGD metric
            DATA = IGDstat;
            titlename = 'IGD';
    end
   
    tl = tiledlayout(1,3);

    
    for poti = [1,2,3]
      
        nexttile;
        for pi = 1:length(PBset)
            loglog(TT(poti,:),squeeze(DATA(pi,:,poti)),'--','Marker',mark(pi),...
                'MarkerSize',5,'Color',colors(pi,:),'DisplayName',['PB',mat2str(pi)]);
            hold on;
            rr = loglog(TT(poti,2:end),squeeze(DATA(pi,2:end,poti)),'-','Marker',mark(pi),...
                'MarkerSize',5,'Color',colors(pi,:),'DisplayName',cell2mat(PBnames(pi)));

            l(pi) = rr;
            hold on;
            xticks(TT(poti,1:2:end))
            xlim([-Inf TT(poti,end)])
            title(names(poti))
            xlabel('$\tau$', 'Interpreter','latex');
            %axis padded
            switch fig
                case 1 %GD metric
                
                    ylim([9e-3 5e+2])

                    rang = [-2 0 2];
                    yti = 10.^(rang);
                    yticks(yti);
                    xticklabels(LAB{poti});
                case 2 %IGD metric
                    ylim([9e-3 0.6]);
                    %yticks(10.^([-2 -1 0]))
                    xticklabels(LAB{poti});
            end
        end
       %ylim([9e-3 0.6]);
       %tl.TileSpacing = 'compact';
    end
    
    leg = legend(l,'Orientation', 'Horizontal','Interpreter','latex');
    leg.Layout.Tile = 'north';

     %title and save
    %titlename2 = sprintf('%s, $d = %u$', titlename,d);
    %titlename2 = ['$', titlename,'$'];
    %title(tl,titlename2, 'Interpreter','latex')

    set(gcf,'units','centimeters','position',[5,5,22,6]);
    filename = sprintf('ParOpt_besttau_%s.pdf',titlename);
    %print(filename,'-dpdf','-r0');
    %exportgraphics(tl,filename)
end
