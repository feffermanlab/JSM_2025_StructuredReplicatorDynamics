%% This script solves the continuous space structured replicator IVP from 
% manuscript associated to this repository. 
clear all

%PlayerSpace
L=1;
res= 1000;
Omega = linspace(0,L,res);

%TimeMesh
tau=1;
T=1000;
time = linspace(0,T,round(T/tau));


%Strategy Space
%A=[1,0;2,-1/2];%....................................Hawk-Dove
%A=[0,-1,1;1,0,-1;-1,1,0]; %.........................RPS
%A=[0,1;-1,0]; %.....................................Trivial Domination
%A=[1,0,0;0,1,0;0,0,1]; %............................Coordination
A=[1,0;0,1]; %......................................Coordination(2)
%A=[0,1,2;-1,0,1;-2,-1,0]; %.........................Discoordination
%A=[0,1,1;1,0,1;1,1,0]; %............................Anticoordination


%number of strategies
m=length(A(1,:));

%familiarity constant (parameter describing the shape of the familiarity kernel) 
s=0.15;

%% initial condition normalized so that u(x) is in the m-1 simplex for all x
%% 3-strategy Initial Conditions

%u0 = [sin(Omega).^2+0.01*ones(1,length(Omega));
%    cos(Omega).^2-0.01*ones(1,length(Omega));
%   ones(1,length(Omega))];

%u0= [Omega.^2+0.01*ones(1,length(Omega));
%    sqrt(Omega)-Omega.^2+0.01*ones(1,length(Omega));
%    ones(1,length(Omega))-sqrt(Omega)+0.01*ones(1,length(Omega))]

%u0 = [0.33+0.2*Omega;0.33-0.2*Omega;0.33*ones(1,length(Omega))];

%u0 = [1/3+0.16*sin(2*pi*Omega);
%      2/3+0.16*sin(2*pi*Omega+2*pi/3)-1/3+0.16*sin(2*pi*Omega);
%       ones(1,length(Omega))-2/3+0.16*sin(2*pi*Omega+2*pi/3)];

%% 2-strategy Initial Conditions
%u0 = [0.5+0.3*sin(2*pi*Omega)+0.1*sin(6*pi*Omega)+0.01*ones(1,length(Omega));
%   0.5-0.3*sin(2*pi*Omega)-0.1*sin(6*pi*Omega)+0.01*ones(1,length(Omega))]

%u0 = [0.5+0.4*Omega.^2.*sin(6*pi*Omega.^2)+0.01;
%      0.5-0.4*Omega.^2.*sin(6*pi*Omega.^2)+0.01];

u0=[0.5+0.2*sin(4*pi*Omega);
    0.5-0.4*sin(4*pi*Omega)];

%u0=[0.5+0.4*Omega.^3;
%    0.5-0.4*Omega.^3];

%u0 = [0.4+0.3*sin(6*pi*Omega.^3)+0.05*Omega+0.01;
%      0.5-0.3*sin(6*pi*Omega.^3)-0.05*Omega+0.01]

%normalize
u0 = u0./sum(u0);

%% solve the nonlocal IVP
Sol = solve(u0,L,res,T,tau,A,s,true);
simplexANI(Sol,2,'FirstTest.gif',["Dove","Hawk"],time,Omega)
%simplexPic(Sol,1,length(time),'HawkDovePatternLong.png',Omega)

%% Functions for the visualization of solutions
function LongSim(u0,Omega,A,s)
%% LongSim this function is used to run the simulation for very long times
% alowing for adjustment of the time mesh as the simulation continues
%
% u0:       Initial Data (vector valued function on Omega)
% Omega     spatial mesh
% A:        Payoff Matrix
%
% This function can be modified to solve the solution is several time
% chunks each with different time resolutions
    res = length(Omega);
    L=Omega(res);
    
    %First Chunk
    Sol = solve(u0,L,res,10000,1,A,s,true);
    u1=Sol(:,:,10000);
    
    %%Second Chunk
    % Sol2= solve(u1,L,res,8000,5,A,s,true);
    % u2=Sol2(:,:,1600);

    %Third Chunk
    % Sol3= solve(u2,L,res,16000,5,A,s,true);
    % u3=Sol3(:,:,3200)

    fig= tiledlayout(1,2);
    fig.TileSpacing = 'compact';
    fig.Padding = 'compact';
    
    nexttile 
    X=[Omega,fliplr(Omega)];
    Y=[zeros(1,res),fliplr(u0(1,:))];
    Z=[u0(1,:),fliplr(ones(1,res))];
    fill(X,Y,[0,0.4470,0.7410]);
    hold on
    fill(X,Z,[0.8500,0.3250,0.0980]);
    hold off
    title('t=0')
    xlabel('\Omega')
    ylabel('\Delta^1')
    
    nexttile
    X=[Omega,fliplr(Omega)];
    Y=[zeros(1,res),fliplr(u1(1,:))];
    Z=[u1(1,:),fliplr(ones(1,res))];
    fill(X,Y,[0,0.4470,0.7410]);
    hold on
    fill(X,Z,[0.8500,0.3250,0.0980]);
    hold off
    xlabel('\Omega')
    legend('Doves','Hawks')
    title('t=1,600')
    
    % nexttile
    % X=[Omega,fliplr(Omega)];
    % Y=[zeros(1,res),fliplr(u2(1,:))];
    % Z=[u2(1,:),fliplr(ones(1,res))];
    % fill(X,Y,[0,0.4470,0.7410]);
    % hold on
    % fill(X,Z,[0.8500,0.3250,0.0980]);
    % hold off
    % xlabel('\Omega')
    % ylabel('\Delta^1')
    % title('t=41,600')
    % 
    % nexttile
    % X=[Omega,fliplr(Omega)];
    % Y=[zeros(1,res),fliplr(u3(1,:))];
    % Z=[u3(1,:),fliplr(ones(1,res))];
    % fill(X,Y,[0,0.4470,0.7410]);
    % hold on
    % fill(X,Z,[0.8500,0.3250,0.0980]);
    % hold off
    % xlabel('\Omega')
    % title('t=81,600')
    
    exportgraphics(gcf,'HawkDoveExtraLong.png','Append',false)
end

function timeSlices(Sol,t1,t2) 
%% timeSlices: This code block generates a figure from two time slices
% from the solution as line graphs superimposed on one another
%
% Sol:      Solution from Solve()
% t1:       Time for frame1
% t2:       Time for frame2
%
    s1 = Sol(:,:,t1);
    s2 = Sol(:,:,t2);
    
    h=plot(Omega,s1(1,:),"--",Omega,s1(2,:),"--",Omega,s1(3,:),"--",Omega,s2(1,:),Omega,s2(2,:),Omega,s2(3,:));
    set(h(1),'Color',"#0000a4")
    set(h(2),'Color','#bc272d')
    set(h(3),'Color','#e9c716')
    set(h(4),'Color',"#0000a4")
    set(h(5),'Color','#bc272d')
    set(h(6),'Color','#e9c716')
end

function simplexANI(Sol,m, filename,legendlabs,time,Omega)
%% simplexANI this function generates an animated gif of the solution through time
%
% Sol           Solution from Solve()
% m             Number of strategies
% filename      Name so save .gif under
% legendlabs    Labels for legend
% time          time mesh used in Solve()
% Omega         Space mesh used in Solve()
    res = length(Omega);
    if(m==2)
        name = filename;
        length(time)
        for i = 1:length(time)
            disp(i)
            si=Sol(1,:,i);
            X=[Omega,fliplr(Omega)];
            Y=[zeros(1,res),fliplr(si)];
            Z=[si,fliplr(ones(1,res))];
            fill(X,Y,[0,0.4470,0.7410]);
            hold on
            fill(X,Z,[0.8500,0.3250,0.0980]);
            hold off
            xlabel('\Omega')
            ylabel('\Delta^1')
            %title(title)
            legend(legendlabs(1),legendlabs(2))
            exportgraphics(gcf,name,'Append',true)
        end
    else
        name = filename;
        for i = 1:length(time)
            si1=Sol(1,:,i);
            si2=Sol(2,:,i);
            X=[Omega,fliplr(Omega)];
            Y1=[zeros(1,res),fliplr(si1)];
            Y2=[si1,fliplr(si1)+fliplr(si2)];
            Y3=[si1+si2,fliplr(ones(1,res))];
            fill(X,Y1,[0,0.4470,0.7410]);
            hold on
            fill(X,Y2,[0.8500,0.3250,0.0980]);
            fill(X,Y3,[0,0.6196,0.4510])
            hold off
            xlabel('\Omega')
            ylabel('\Delta^2')
            legend(legendlabs(1),legendlabs(2),legendlabs(3))
            %title(gca,title)
            exportgraphics(gcf,name,'Append',true)
        end
    end
end

function simplexPic(Sol,t1,t2,filename,Omega)
%% SimplexPic: This function visualizes a solution in such a way that 
% emphasizes each players' strategy is in the m-1 simplex
%
% Sol:      Solution from Solve()
% t1:       Time for frame1
% t2:       Time for frame2
% Filename: Name to same png
% Omega     Spatial Domain
%
% Visualize the solutions of the structured replicator equation in two time
% slices side by side
    res = length(Omega)
    si=Sol(1,:,t1);
    sf =Sol(1,:,t2);
    fig= tiledlayout(1,2);
    fig.TileSpacing = 'compact';
    fig.Padding = 'compact';
    nexttile 
    X=[Omega,fliplr(Omega)];
    Y=[zeros(1,res),fliplr(si)];
    Z=[si,fliplr(ones(1,res))];
    fill(X,Y,[0,0.4470,0.7410]);
    hold on
    fill(X,Z,[0.8500,0.3250,0.0980]);
    hold off
    xlabel('\Omega')
    ylabel('\Delta^1')
    nexttile
    X=[Omega,fliplr(Omega)];
    Y=[zeros(1,res),fliplr(sf)];
    Z=[sf,fliplr(ones(1,res))];
    fill(X,Y,[0,0.4470,0.7410]);
    hold on
    fill(X,Z,[0.8500,0.3250,0.0980]);
    hold off
    xlabel('\Omega')
    exportgraphics(gcf,filename,'Append',false)
end

%% Function to solve non-local equation
function Sol=solve(U0,L,res,T,tau,A,s,circ)
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
% circ:     Boolean variable which determines if the domain is periodic
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
    if circ
        k= circshift(1/(s*sqrt(2*pi))*exp((Omega-(1/2)).^2./(-2*s^2)),round(length(Omega)/2));
        for t = 2:finaltic
            g = zeros(m,res);
            u = Sol(:,:,t-1);
            for j = 1:m
                uj= Sol(j,:,t-1);
                g(j,:)=cconv(k,uj,res)*spacemesh;
                
            end
            resolvent = repmat(diag(u'*A*g)',m,1);
            du = u.*(A*g-resolvent);
            Sol(:,:,t)=Sol(:,:,t-1)+tau*du;
        end
               
    else
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
end

%% Familiarity Kernel
function K=ker(x,s)
    K = (1/(s*sqrt(2*pi))*exp((x).^2./(-2*s^2)))-(1/(s/2*sqrt(2*pi))*exp((x).^2./(-2*(s/2)^2)));
end


