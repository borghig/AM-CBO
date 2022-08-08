function [output] = AM_CBO_Mtimes(Mtest, problem, ev, in)
%{Test the AM-CBO algorithm 
% to the same given problem 
% "Mtest" times to collect statistics 
%}

% tests the same problem Ntest times
N = length(in.w0);

%initialization of quanties of interest 
X = zeros(N,problem.dim, Mtest);  %store positions
W = zeros(N,1,           Mtest);  %store fronts
G = zeros(N,2,           Mtest);  %store parameters

ZEROS = zeros(ev.kmax,1);

Ustat   = ZEROS;  %currently used Energy
GDstat  = ZEROS;  %Generational Distance
IGDstat = ZEROS;  %inverted Generational Distance
RIEstat = ZEROS;  %Riesz energy
NEWstat = ZEROS;  %Newtonian energy
MORstat = ZEROS;  %Riesz energy
hypstat = 0;  %hypervolume contribution

output.igd = zeros(Mtest,1);

parfor m = 1:Mtest
    
    % copy input locally for parallel computation 
    problem_loc = problem;
    in_loc = in;
    ev_loc = ev;
    
    % randomly sample intial data
    in_loc.x0 = rand(N,problem_loc.dim).*(problem_loc.ub-problem_loc.lb) + problem_loc.lb;

    
    % run the algorithm
    test = AM_CBO(problem_loc,ev_loc,in_loc);
    X(:,:,m) = test.x;
    W(:,:,m) = test.w;
    G(:,:,m) = test.g;

    %save statistics
    igd(m)   = test.IGDtot(end);
    Ustat(m) = test.Utot(end);
    
    GDstat  = GDstat  + test.GDtot /Mtest;
    IGDstat = IGDstat + test.IGDtot/Mtest;
    RIEstat = RIEstat + test.enRIE /Mtest;
    NEWstat = NEWstat + test.enNEW /Mtest;
    MORstat = MORstat + test.enMOR /Mtest;
    hypstat = hypstat + test.hypTOT/Mtest;
end


% store best run
[bestIGD, besti] = min(igd);
output.bestIGD = bestIGD;
output.bestx   = squeeze(X(:,:,besti));
output.bestw   = squeeze(W(:,:,besti));
output.bestg   = squeeze(G(:,:,besti));

% return computed statistics 
output.Ustat   = Ustat;
output.GDstat  = GDstat ;
output.IGDstat = IGDstat;
output.Estat   = Ustat;
output.RIEstat = RIEstat;
output.NEWstat = NEWstat;
output.MORstat = MORstat;
output.hypstat = hypstat;

end





