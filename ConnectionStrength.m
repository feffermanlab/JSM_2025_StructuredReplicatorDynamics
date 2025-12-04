%% This script finds critical levels of connectedness between two modules
% in a graph depending which change the stability of game theoertic
% equilibria. It is made obsolete by other code in this repository (i.e.
% EdgeSensitivity.m)

% size of the two modules
n=6;
m=4;

% Internal adjacency matrix for each module
Kn= ones(n,n)-eye(n);
Km= ones(m,m)-eye(m);

% Intermodular connectivity adjacency matrix
Knm = [1,1,1,1; 0,1,1,1;0,0,1,1;0,0,0,1;0,0,0,0;0,0,0,0];

% construct strategy profile so that every player in a module plays the same
% strategy 
un = repmat([1;0],1,n);
um = repmat([0;1],1,m);
u = [un,um];

% 2x2 Payoff matrix
A=[1,0;0,1]; %............................Coordination
%A=[0,1;1,0] %............................Anti-Coordination

%Overall Adjacency matrix
W = [Kn,Knm;(Knm'),Km];

%prealocate space for the result
res = zeros(101,1);
for i = 0:100
    a=i/100;
    %construct adjacency matrix with intermediate level of connectedness 
    W=[Kn,a*Knm;(a*Knm'),Km];
    %determine stability 
    lambda = maxeig(W,A,u,n,m);
    res(i+1)=lambda;
end

%Find where the result crosses the x-axis (if that occurs)
[M,I]= max(-abs(res));

%produce figure
figure()
X=linspace(0,1,101);
h=plot(X,res,X,zeros(101,1),0.75,0);
hold on
plot(0.75,0,'o','Color','black')
plot([0.75,0.75],[-3,1],'LineStyle','--')
text(0.275,-2.75,'Critical Connectedness level: 0.75')
set(h(2),'Color','black')
ylabel('Maximum Eigenvalue')
xlabel('intermodule edge weight')
title('Connectivity between K_6 and K_4')



function lambda = maxeig(W,A,u,n,m)
%% Function to determine the maximum eigenvalue of an equilibrium strategy
%  profile 
% W:        Adjacency matrix
% A:        Payoff matrix
% u:        Strategy Profile
% n:        Number of Players in module 1
% m:        Number of Players in module 2

%%This function computes the eigenvalues of the jacobian matrix for the ODE
%%system described in the manuscript and identifies the maximum eigenvalue. 
    % determine cumulative strategies in the nieghborhood of each player
    g=zeros(2,n+m);
    for v=1:(n+m)
        for w =1:(n+m)
            g(:,v) = g(:,v)+ A*(W(v,w)*u(:,w));
        end
    end
    
    % generate the jacobian matrix through method described in manuscript
    J=zeros(2*(n+m),2*(n+m));
        for v=1:(n+m)
            for w = 1:(n+m)
                if v==w
                    TMP = diag(g(:,v)-(u(:,v)'*g(:,v))*ones(2,1))-u(:,v)*g(:,v)';
                else
                    TMP = W(w,v)*(diag(u(:,v))-u(:,v)*u(:,v)');
                end
                J((v-1)*2+1:v*2,(w-1)*2+1:w*2)=TMP;
            end
        end
        %get eigen system for Jacobian matrix and return the maximum real
        %real part of the eigenvalues
        [V,D]= eig(J);
        lambda = max(real(sum(D,1)));
end