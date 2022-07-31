function [output] = AM_CBO(problem, ev, in)

%{

    %(remark: w represent the first component of the weight vector, the
    % second one is uniquely determined as 1-w)
%}

dim = problem.dim;
N = length(in.w0);
x = in.x0;
w = in.w0;
g = problem.g(x);
g1 = g(:,1); 
g2 = g(:,2);

v.potLenght = ev.potLenght;


% initilize data
ZEROS = zeros(ev.kmax , 1); 
output.Utot   = ZEROS ;
output.IGDtot = ZEROS ;
output.GDtot  = ZEROS ;
output.enMOR  = ZEROS ;
output.enNEW  = ZEROS ;
output.enRIE  = ZEROS ;

yalpha_s = x;

% algorithm iteration
k = 1;
while k <= ev.kmax   

    xnew = x;   
       
    % position update
    for i = 1:N
        % Chebyschev scalarization apporach
        gw = max([w(i)*g(:,1),(1-w(i))*g(:,2)],[],2);       
        yalpha_s(i,:) = compute_yalpha(x,gw, ev.alpha);

        % add deterministic component
        diff =  yalpha_s(i,:) - x(i,:); 
        xnew(i,:) = x(i,:) + ev.lambda*ev.deltat*diff;
        
        % add explorative component 
        switch ev.diffusion
            case 'iso'
                xnew(i,:) = xnew(i,:) + ev.sigma*sqrt(ev.deltat)*norm(diff)*randn(1,dim);
            case 'aniso'
                xnew(i,:) = xnew(i,:) + ev.sigma*sqrt(ev.deltat).*diff.*randn(1,dim);
        end        
    end
    x = xnew;


    % parameters adaptation
    [DU, U] = compute_potential(g,ev);     %compute potential forces
    wnew = w + ev.deltat.*ev.tau.*(DU(:,2) - DU(:,1));
    wnew(wnew<0) = 0;              %projection to \Omega
    wnew(wnew>1) = 1; 
    w = wnew;
     
    g = problem.g(x);
    
    %store evolution of performance metrics
    output.IGDtot(k) = compute_IGD(problem.F,g);
    output.GDtot(k)  = compute_IGD(g,problem.F);
    output.Utot(k)   = U;
    
    v.potential = 'Riesz';
    v.potLenght = 1;  
    [~,E] = compute_potential(g,v);
    output.enRIE(k) = E;
    
    v.potential = 'Newtonian';
    v.potLenght = 1;
    [~,E] = compute_potential(g,v);
    output.enNEW(k) = E;
    
    v.potential = 'Morse';
    v.potLenght = 20;
    [~,E] = compute_potential(g,v);
    output.enMOR(k) = E;

    k=k+1;
end

% to compute hypervolume contribution, eliminate eventually non-dominated
% points
g_ref = faster_pareto2(g');
g_ref(g_ref>problem.max(1)) = problem.max(1); %since problem.max(1) = problem.max(2)
output.hypTOT = lebesgue_measure(g_ref, problem.max);

%return computed solution and final particles 
output.x = x;
output.w = w;
output.g = g;


