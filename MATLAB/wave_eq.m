%% HOMEWORK 25
% EQUAZIONE DELLE ONDE 1D CON SCHEMA LEAPFROG HYPERBOLIC SYSTEM
%    d^2 u       d^2 u        
%   ------- - c ------- = f
%     dt^2        dx^2      

% Definiamo la geometria del problema con gli intervalli:
% [ X_SX; X_DX ] x [ T_IN, T_FIN ]

% Le condizioni al contorno in termini spaziali-temporali sono
% INITIAL VALUE PROBLEM  --->  u(x,0) = u0
% BOUNDARY VALUE PROBLEM  -->  u(x,t) = ub

% Il dominio è discretizzato in:
% Spazio --> xd = suddivisione in intervalli uguali di [ X_SX; X_DX ]
% Tempo ---> td = suddivisione in intervalli uguali di [ T_IN, T_FIN ]

% Il parametro di governo del problema sarà direttamente descritto da CFL,
% nel caso della equazione d'onda la condizione di stabilità è verificata
% quando CFL = dt*sqrt(c)/h <= 1

% INITIAL VALUE PROBLEM e BOUNDARY VALUE PROBLEM sono descritti da una
% condizione iniziale u0(x)=exp(-10*x^2) e v0(x)=0, considerando il caso
% particolare di f=0, condizioni di Dirichlet omogenee u(a,t)=u(b,t) per
% ogni istante di tempo, e derivata temporale nulla per t=0

% INITIAL VALUE PROBLEM e BOUNDARY VALUE PROBLEM del sistema iperbolico
% sono descritti da condizioni iniziali legate a w1 e w2, infatti posto
% w1=du/dx e w2=du/dt, si pone w2(x,0)=u0'(x) e w2(x,0)=v0(x)


%% INIZIALIZZAZIONE VALORI DEL PROBLEMA, DOMINIO, E DISCRETIZZAZIONE
clear all
close all
clc

% Dominio del problema
t_in=0;
t_fin=7;

x_sx=-2;
x_dx=2;

% Parametri numerici del problema e del sistema iperbolico
f=0;
c=1;

A=[0 -1;-sqrt(c) 0];
Aq=A^2;

h=0.04;
CFL=0.8;
dt = CFL*h/max(eig(A));

% Discretizzazione del dominio
Nt=fix((t_fin-t_in)/dt+1);
td=linspace(t_in,t_fin,Nt);
Nx=fix((x_dx-x_sx)/h+1);
xd=linspace(x_sx,x_dx,Nx);


%% IMPLEMENTAZIONE SCHEMI NUMERICI DIFFERENTI
% Inizializzazione della soluzione: le righe della soluzioni rappresentano
% la soluzione al tempo t-esimo allo scorrere delle colonne in cui è
% memorizzata la soluzione alla posizione x-esima

% Inizializzazione soluzione, comando: u=zeros(Nt,Nx);
% Initial Value Problem, comando:  u(1,:)=u0(xd);
% Boundary Value Problem, comando: u(:,1)=ub(x);


% RISOLUZIONE CON DIVERSI SCHEMI NUMERICI
for scheme=1:3
    % Si inizializza la soluzione
    u=zeros(Nt,Nx);
    
    % INITIAL VALUE PROBLEM: soluzione per t=0 data da u0, si impone quindi
    % all'istante iniziale per ogni x u=u0
    u(1,:)=u0(xd);
    % BOUNDARY VALUE PROBLEM: supponiamo la condizione iniziale data da
    u(:,1)=ub(x_sx);
    u(:,end)=ub(x_dx);
    
    % INITIAL VALUE PROBLEM: si impone derivata temporale nulla per t=0
    t_lim=0;
    n=1;
    u(n+1,:)=u(n,:);

    % SCHEMA LEAP-FROG
    % Fonte: Quarteroni "Calcolo Scientifico"
    if scheme==1
        n=2;
        while t_lim<t_fin
            for j = 2:Nx-1
                u(n+1,j)=2*u(n,j)-u(n-1,j)+c*(dt/h)^2*(u(n,j+1)-2*u(n,j)+u(n,j-1));
            end
            t_lim=t_lim+dt;
            n=n+1;
        end
        u_leapfrog = u;
    end
    
    
    % SCHEMA LAX-WENDROFF
    % Fonte: Quarteroni "Calcolo Scientifico"
    if scheme==2

        % Costruzione del sistema iperbolico
        w1 = zeros(Nt,Nx);
        w2 = zeros(Nt,Nx);

        w1(1,:)=ub1(xd);
        w2(1,:)=zeros(1,Nx);
        
        while t_lim<t_fin
            
                w1(n+1,1)=w1(n,1)+dt/h*(w2(n,2)-w2(n,1));
                w1(n+1,end)=w1(n,end)+dt/h*(w2(n,end)-w2(n,end-1));

            for j = 2:Nx-1
                w1(n+1,j)=w1(n,j)-0.5*(dt/h)*(A(1,1)*(w1(n,j+1)-w1(n,j-1))+A(1,2)*(w2(n,j+1)-w2(n,j-1)))+0.5*(dt/h)^2*(Aq(1,1)*(w1(n,j+1)-2*w1(n,j)+w1(n,j-1))+Aq(1,2)*(w2(n,j+1)-2*w2(n,j)+w2(n,j-1)));
                w2(n+1,j)=w2(n,j)-0.5*(dt/h)*(A(2,1)*(w1(n,j+1)-w1(n,j-1))+A(2,2)*(w2(n,j+1)-w2(n,j-1)))+0.5*(dt/h)^2*(Aq(2,1)*(w1(n,j+1)-2*w1(n,j)+w1(n,j-1))+Aq(2,2)*(w2(n,j+1)-2*w2(n,j)+w2(n,j-1)));
                u(n+1,j)=u(n,j)+dt/2*(w2(n+1,j)+w2(n,j));
            end
            
        t_lim=t_lim+dt;
        n=n+1;    
        end
        u_lw = u;
    end
    
    %METODO LAX-FRIEDRICHS
    % Fonte: Quarteroni "Calcolo Scientifico"
    if scheme==3

        % Costruzione del sistema iperbolico
        w1 = zeros(Nt,Nx);
        w2 = zeros(Nt,Nx);

        w1(1,:)=ub1(xd);
        w2(1,:)=zeros(1,Nx);

        while t_lim<t_fin

            w1(n+1,1)=w1(n,1)+dt/h*(w2(n,2)-w2(n,1));
            w1(n+1,end)=w1(n,end)+dt/h*(w2(n,end)-w2(n,end-1));

            for j = 2:Nx-1
                w1(n+1,j)=0.5*(w1(n,j+1)+w1(n,j-1))+0.5*CFL*(w2(n,j+1)-w2(n,j-1));
                w2(n+1,j)=0.5*(w2(n,j+1)+w2(n,j-1))+0.5*CFL*(w1(n,j+1)-w1(n,j-1));
                u(n+1,j)=u(n,j)+dt*w2(n,j);        
            end
            
            t_lim=t_lim+dt;
            n=n+1;
        end
        u_lf=u;
    end
    
end

%% PLOT DELLE SOLUZIONI
for m=1:length(td)
    figure(1)
    % Plot della soluzione con schema LEAP-FROG
    plot(xd,u_leapfrog(m,:),'b','LineWidth',1),hold on
    axis([-2 2 -1 1])
    % Plot della soluzione con schema LAX-WENDROFF
    plot(xd,u_lw(m,:),'r','LineWidth',1), hold on
    % Plot della soluzione con schema LAX-FRIEDRICHS
    plot(xd,u_lf(m,:),'color',[0.9290 0.6940 0.1250],'LineWidth',1), hold off

    title('WAVE EQUATION: u_{tt}=c^2u_{xx} - comparison between schemes (CFL = 0.8)')
    legend('LF','Lax-W','Lax-F')
end

%% INITIAL VALUE AND BOUNDARY VALUE PROBLEM FUNCTIONS
% Function per le condizioni al contorno: variando i function handle
% presenti all'interno delle function di Matlab è possibile variare il
% problema ai valori iniziali.

function [u_ivp] = u0(x)

   f_ivp=@(v) exp(-10*v.^2);
   u_ivp=f_ivp(x);

end

function [u_bvp] = ub(x)

   f_ivp=@(v) exp(-10*v.^2);
   u_bvp=f_ivp(x);
   
end

function [u_bvp1] = ub1(x)

   f_bvp1=@(t) -20*t.*exp(-10*t.^2);
   u_bvp1=f_bvp1(x);
   
end
