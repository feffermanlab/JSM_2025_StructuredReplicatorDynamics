%% This script solves the continuous space structured replicator IVP from 
% manuscript associated to this repository. 
clear all

%PlayerSpace
L=1;
res= 500;
Omega = linspace(0,L,res);

%TimeMesh
tau=0.1;
T=450;
time = linspace(0,T,T*round(1/tau));


%Strategy Space
%A=[0,-1,1;1,0,-1;-1,1,0]; %.........................RPS
%A=[0,1;-1,0]; %.....................................Trivial Domination
%A=[1,0,0;0,1,0;0,0,1]; %............................Coordination
%A=[0,1,2;-1,0,1;-2,-1,0]; %.........................Discoordination
A=[0,1,1;1,0,1;1,1,0]; %............................Anticoordination

%number of strategies
m=length(A(1,:));

%familiarity constant (parameter describing the shape of the familiarity kernel) 
s=0.01;

%initial condition normalized so that u(x) is in the m-1 simplex for all x

%3-strategy
u0 = [sin(Omega).^2+0.01*ones(1,length(Omega));
    cos(Omega).^2-0.01*ones(1,length(Omega));
    ones(1,length(Omega))];

%2-strategy
%u0 = [sin(Omega).^2+0.01*ones(1,length(Omega));
% cos(Omega).^2-0.01*ones(1,length(Omega))]

%normalize
u0 = u0./sum(u0);

%solve the nonlocal IVP
Sol = solve(u0,L,res,T,tau,A,s);

%%This code block generates a figure from two time slices from the solution
% uncomment to run 

% s1 = Sol(:,:,500);
% s2 = Sol(:,:,520);
% 
% h=plot(Omega,s1(1,:),"--",
%       Omega,s1(2,:),"--",
%       Omega,s1(3,:),"--",
%       Omega,s2(1,:),
%       Omega,s2(2,:),
%       Omega,s2(3,:));
% set(h(1),'Color',"#0000a4")
% set(h(2),'Color','#bc272d')
% set(h(3),'Color','#e9c716')
% set(h(4),'Color',"#0000a4")
% set(h(5),'Color','#bc272d')
% set(h(6),'Color','#e9c716')


%%This code block generates an animated figure of the solution running in
%time (This can take several seconds depending on the time mesh because the
%numerical methods are not overly efficient.

h=figure;
s = Sol(:,:,1);
%Start by plotting initial data
plot(Omega,s(1,:),Omega,s(2,:),Omega,s(3,:)) %......3-strategy
%plot(Omega,s(1,:),Omega,s(2,:)) %...................2-strategy

%Must determine the axis before animation 
axis([0 L 0 1])

disp("here")
ax = gca;
ax.NextPlot = 'replaceChildren';
 
loops = length(time);
M(loops) = struct('cdata',[],'colormap',[]);
 
for j = 1:loops
    s= Sol(:,:,j);
    plot(Omega,s(1,:),Omega,s(2,:),Omega,s(3,:)) %.......3-strategy
    %plot(Omega,s(1,:),Omega,s(2,:)) %....................2-Strategy
    drawnow
    M(j) = getframe;
end
 
movie(M);



function Sol=solve(U0,L,res,T,tau,A,s)
%% Solve: This function solves the non-local Initial value problem as described in
% the manuscript associated with this repository
%
% U0:       Initial Data (vector valued function on Omega)
% L:        Length of Omega
% res:      Spatial Resolution of Omega
% T:        Total time 
% tau:      width of each interval in the time mesh
% A:        Payoff Matrix
% s:        Familiarity Parameter
%
% Solve the initial value problem for the parabolic nonlinear equation
% using a standard forward eular method in time and right point quadrature
% for the nonlocality
    
    %Estabilish the space time mesh
    finaltic = round(T/tau);
    finalplayer = res;
    spacemesh = L/(res-1);
    Omega = linspace(0,L,res);

    %number of strategies
    m = length(A(1,:)); 

    %preallocate solution
    Sol = zeros(m,finalplayer, finaltic);
    Sol(:,:,1)=U0;
    
    for t = 2:finaltic
        for i=1:finalplayer
            u = Sol(:,:,t-1);

            %compute nonlocality
            k = ker(Omega(i)-Omega,s);
            gu = sum(u.*repmat(k, [m 1]),2)*spacemesh;
            
            % compute time derivaitve as directed by nonlocal equation in
            % the associated manuscript
            ux=u(:,i);
            du= ux.*(A*gu-ux'*A*gu);
            
            Sol(:,i,t)=Sol(:,i,t-1)+tau*du;
        end
    end
end



%Familiarity Kernel
function K=ker(x,s)
    K = 1/(s*sqrt(2*pi))*exp((x).^2./(-2*s^2));
end


