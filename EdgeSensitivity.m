%% This script includes code to estimnate how edge addition and edge removal
% can change the stability of game theoretic equilibria 


%Adjacency matrix
W=[0,1,0,0,0,0,1,1;
    1,0,1,0,0,0,1,1;
    0,1,0,0,0,0,1,1;
    0,0,0,0,1,0,1,1;
    0,0,0,1,0,1,1,1;
    0,0,0,0,1,0,1,1;
    1,1,1,1,1,1,0,0;
    1,1,1,1,1,1,0,0];

%number of players
n = length(W(1,:));

%Payoff matrix 
%A = [1,0;0,1]; %................................2-strategy Coordination
%A=[0,-1,2;2,0,-1;-1,2,0]; %.....................Skew RPS
A=[10,1,-10;1,0,1;-10,1,10]; %..................Arbitrary

%number of strategies
m=length(A(:,1));

%strategy profile in question
%must be mXn matrix
u = [1,1,1,0,0,0,0,0;
    0,0,0,0,0,0,1,1;
    0,0,0,1,1,1,0,0];

% Determine discrete sensitivity for each edge and non-edge
DW = EdgeSensitivityDis(W,u,A)

% Determine network which "optimizes" stability of the strategy profile u
W=WOptimize(u,A)
J=Jacobian(W,u,A);
l = max(real(eig(J)))




function Wopt = WOptimize(u,A,tol,iterlimit)
%% WOptimize: basic numerical method to descend the weak gradiant of the
%maximal real part of the eigenvalues in the space of adjacency matrices
%for a particular strategy profile
%
% u:        strategy profile
% A:        Adjacency matrix
% tol:      Tolerence for convergence
%iterlimit: Maximum number of steps
%
% This method takes a very basic mixed forward and backward gradient
% estimation method and to reduce the maximum real part of the eigvalues of
% the Jacobian evaluated at the strategy profile u

    %set defaults
    if(nargin<2)
        disp("insufficient arguments")
    elseif(nargin==3)
        iterlimit = 10;
    else 
        tol = 10^(-5);
        iterlimit = 10000;
    end

    %hardcoded space and time mesh
    dx = 0.001;
    dt = 0.001;

    %number of players
    n=length(u(1,:));

    %number of strategies
    m=length(A(:,1));

    %random noise at the start prevents getting caught in unstable Nash
    %equilibria
    rng("default")

    %Produce a symmetric, simple, random, weighted adjacency matrix.
    r = rand(n,n);
    diag(diag(r));
    r=r-diag(diag(r));
    W0=1/2*(r+r');
    
    %establish flags for finishing process
    iter = 1;
    finished = false;
    Wi = W0;
    toleven=1;
    tolodd=1;
    while not(finished)
        if(mod(iter,2)==0)
            % on even steps use compute the forward difference
            DlDW= EdgeSensitivityContinuous(Wi,dx,u,A);
            iter = iter+1;
            %ensure nothing leaves the bounds of the adjacency matrix space
            Wnew = Wi-dt*DlDW-dt*DlDW';
            Wnew = min(Wnew,ones(n,n));
            Wnew = max(Wnew,zeros(n,n));
            %test for convergence 
            toleven = max(abs(Wnew-Wi),[],[1,2]);
        else
            % on odd steps compute the backward difference
            DlDW = EdgeSensitivityContinuous(Wi,-dx,u,A);
            iter = iter+1;
            %ensure nothing left the bounds of the adjacency matrix space
            Wnew = Wi-dt*DlDW-dt*DlDW';
            Wnew = min(Wnew,ones(n,n));
            Wnew = max(Wnew,zeros(n,n));
            %test for convergence
            tolodd = max(abs(Wnew-Wi),[],[1,2]);
        end
        
        if((toleven<tol && tolodd<tol) || iter >iterlimit)
            finished =true ;
        end
        Wi = Wnew;
    end
    Wopt= Wi;
end


function DW = EdgeSensitivityDis(W,u,A)
%% EdgeSensitivityDis: Determines the change the maximum real part of the
% eigenvalues under edge addition and removal
%
% W:        Adjacency matrix
% u:        Strategy Profile
% A:        Payoff matrix
%
% In the resulting matrix each entry is the change in maximum real part of
% the eigenvalues under edge removal (if that edge was in W) or removal (if
% that edge was not in W).

    
    % number of players
    n = length(W(:,1));

    %preallocate space for result
    DW = zeros(n,n);

    %Compute Initial value
    J= Jacobian(W,u,A)
    lambda0= max(real(eig(J)))

    for v = 1:n
        for w = v+1:n
            %compute new adjacency matrix with edge removed or added
            X = W;
            X(v,w)=not(W(v,w));
            X(w,v)=not(W(w,v));
            
            %compute new jacobian and eigenvalue
            J= Jacobian(X,u,A);
            lambdaprime= max(real(eig(J)));
            
            %add entry to result
            DW(v,w)=lambdaprime-lambda0;
        end
    end
end

function DlDW = EdgeSensitivityContinuous(W,dx,u,A)
%% EdgeSensitivityDis: Estimate the (weak) derivative for each entry in the 
% adjacency matrix through forward (backwards) difference
%
% W:        Adjacency matrix
% u:        Strategy Profile
% A:        Payoff matrix
%
% The convergence of this method is not necessary certain because of the
% discontinuous derivative of the maximum real part of the eigenvalues.
   
    %number of players
    n = length(W(:,1));

    %preallocate gradiant 
    DlDW=zeros(n,n);

    %compute Original Jacobian and eigenvalue
    J=Jacobian(W,u,A);
    l=max(real(eig(J)));
            
    for v = 1:n
        for w = v+1:n
            %generate new adjacency matrix with dx change in the weight of
            %edge (v,w) and ensure that it does not escape the boundes of 
            %the adjacency matrix space
            Xplus = W;
            Xplus(v,w)=min(W(v,w)+dx,1);
            Xplus(v,w)=max(W(v,w)+dx,0);
            Xplus(w,v)=min(W(w,v)+dx,1);
            Xplus(w,v)=max(W(w,v)+dx,0);
           
            %Compute Jacobain and eigenvalue
            Jplus= Jacobian(Xplus,u,A);
            lplus = max(real(eig(Jplus)));
            
            %compute forward difference 
            DlDW(v,w)= (lplus-l)/(dx);
        end
    end
end


function J = Jacobian(W,u,A)
%% Jacobian: computes the jacobian matrix for the ODE system in the associated
% manuscript
%
% W:        Adjacency matrix
% u:        strategy profile
% A:        Payoff matrix
%
%The computation of the Jacobian is detailed in the associated matrix; it
%is an analytical result, not an approximnation of the Jacobian 
    
    %number of players
    n = length(W(:,1));

    %number of strategies
    m = length(A(:,1));
    
    %Compute cumulative strategies in the neighborhood of each player v
    g=zeros(m,n);
    for v=1:n
        for w =1:n
            g(:,v) = g(:,v)+ A*(W(v,w)*u(:,w));
        end
    end

    %Determine that the strategy profile is indeed at equilibrium
    equilibrium = true;
    %ODE determined described in the manuscript
    dudt = u(:,v).*(g(:,v)-(u(:,v)'*g(:,v))*ones(m,1));
    if any(dudt)
        equilibrium = false;
    end

    %If it is an equilibrium compute Jacobian as described in the
    %manuscript
    if equilibrium
        J=zeros(n*m,n*m);
        for v=1:n
            for w = 1:n
                if v==w
                    TMP = diag(g(:,v)-(u(:,v)'*g(:,v))*ones(m,1))-u(:,v)*g(:,v)';
                else
                    TMP=W(v,w)*(A.*(u(:,v)*ones(1,m))-u(:,v)*(A*u(:,v))');
                end
                J((v-1)*m+1:v*m,(w-1)*m+1:w*m)=TMP;
            end
        end
    else
        %If it is not an equilibrium, give up
        disp("not equilibrium")
        J= zeros(n*m,n*m);
    end
end
