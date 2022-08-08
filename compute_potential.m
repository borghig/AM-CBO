function [DU,U,Uvec] = compute_potential(g,ev)

% given g (size Np x m)
% returns the gradient fo the potential DU (size Np x m)
% and potential itself

C1 = 1;
C2 = ev.potLenght;
        

N = length(g);

 
a1 = g(:,1);
A1 = repmat(a1,1,N);
D1 = A1 - A1'; % D1(i,j) = a(i,1) - a(j,1)

a2 = g(:,2);
A2 = repmat(a2,1,N);
D2 = A2 - A2'; % D2(i,j) = a(i,2) - a(j,2)

% D(i,j,:) = a(i,:) - a(j,:) % after squeezing and transposing
D = zeros(N, N, 2);
D(:,:,1) = D1;
D(:,:,2) = D2;

% ND(i,j) = |a(i,:) - a(j,:)|
ND = vecnorm(D,2,3);

switch ev.potential
    case 'Riesz'
        invND = 1./ND;
        invND(invND ==Inf) = 0;

        nom = zeros(size(D)); 
        nom(:,:,1) = invND;
        nom(:,:,2) = invND;
        
        DDU = D.*(nom.^3);
        DU = C1*squeeze(sum(DDU,2)/(N-1));
        %DU = real(DU);
        
        % potential 
        llog = invND;
        llog(llog ==Inf) = 0;
        llog(llog ==-Inf) = 0;
        U = C1*sum(llog,'all')/(N*(N-1));

    case 'Newtonian'
        % AKA Newtonian potential
        
        % here we assume d==2
        % if not, this needs to change
        
        % gradient
        invND = 1./ND;
        invND(invND ==Inf) = 0;

        nom = zeros(size(D)); 
        nom(:,:,1) = invND;
        nom(:,:,2) = invND;
        
        DDU = D.*(nom.^2); %should be .^2 here, but was .^3
        DU = C1*squeeze(sum(DDU,2)/(N-1));
        %DU = real(DU);
        
        % potential 
        llog = -log(ND);
        llog(llog ==Inf) = 0;
        llog(llog ==-Inf) = 0;
        U = C1*sum(llog,'all')/(N*(N-1));
        
    case 'Morse'
        % DU  = -(C1*C2*x/|x|)e^(-|x|*C2)/2

        %calculate e^(-|x|)
        ex = C1.*exp(-ND.*C2);
        expo = zeros(size(D)); 
        expo(:,:,1) = ex;
        expo(:,:,2) = ex;
        
        %calculate -(x/|x|)
        nom  = zeros(size(D)); 
        nom(:,:,1)  = ND;
        nom(:,:,2)  = ND;
        nom(nom==0) = 1; %doesnt matter what we put here, but cannot be 0
        
        %return DU
        DDU = C2.*D./nom.*expo./2;
        DU = squeeze(sum(DDU,2)/(N-1)); %this sums all the entries
        %DU = real(DU);
       
        % potential 
        U = sum(expo,'all')/(N*(N-1));  
        
        Uvec = sum(ex,2);
    case 'NoInt'
        DU = zeros(N,2);
        U = 0;
        Uvec = zeros(N,1);
end
  
end

