%% HOMEWORK 23
% EQUAZIONE DI AVVEZIONE CON
% CON SCHEMA UPWIND AL SECOND'ORDINE
% 
% Risolvere l'equazione di avvezione
%     du/dt + a du/dx = 0
% con schema upwind al secondo ordine

% Definiamo la geometria del problema con gli intervalli:
% [ X_SX; X_DX ] x [ T_IN, T_FIN ]

% Le condizioni al contorno in termini spaziali-temporali sono
% INITIAL VALUE PROBLEM  --->  u(x,0) = u0
% BOUNDARY VALUE PROBLEM  -->  u(X_SX,t) = ub

% Il dominio è discretizzato in:
% Spazio --> xd = suddivisione in intervalli uguali di [ X_SX; X_DX ]
% Tempo ---> td = suddivisione in intervalli uguali di [ T_IN, T_FIN ]

% Il parametro di governo del problema sarà direttamente descritto da CFL,
% non si richiede in ingresso la velocità "a" prevista nella equazione
% ricordando la condizione CFL = A dt/dx < 1

% INITIAL VALUE PROBLEM e BOUNDARY VALUE PROBLEM descritti a pa. 421 di
% Quarteroni "Calcolo Scientifico", eq.(9.98) e commentati dove necessario
% all'interno di questo codice


clear all
close all
clc


%% SELEZIONE DEL PROBLEMA
% CONDIZIONE INIZIALE: SCEGLIERE 1 PER l=1/2
%                      SCEGLIERE 2 PER l=1
initial_condition = 2;


%% INIZIALIZZAZIONE VALORI DEL PROBLEMA, DOMINIO, E DISCRETIZZAZIONE
% Dominio del problema
t_in=0;
t_fin=10;

x_sx=-1;
x_dx=3;

if initial_condition == 1
    % Parametri numerici del problema: caso 1
    l=1/2;
    a=1;
    CFL=0.8;    
    h=l/8;
elseif initial_condition == 2
    % Parametri numerici del problema: caso 2
    l=1;
    a=1;
    CFL=0.8;    
    h=l/20;
end

dt=1/a*CFL*h;

Nt=ceil((t_fin-t_in)/dt+1);
td=linspace(t_in,t_fin,Nt);
Nx=(x_dx-x_sx)/h+1;
xd=linspace(x_sx,x_dx,Nx);

%% BOUNDARY CONDITIONS E IMPOSIZIONE DEL PROBLEMA
% Inizializzazione della soluzione: le righe della soluzioni rappresentano
% la soluzione al tempo t-esimo allo scorrere delle colonne in cui è
% memorizzata la soluzione alla posizione x-esima

% Comando: u=zeros(Nt,Nx);

% INITIAL VALUE PROBLEM: supponiamo la condizione iniziale data da
%             {  sin(2*pi*x/L) se -L < x < L
%   u0(x,0) = {
%             {      0         se  L < x < 3

% Comando: u(1,:)=u0(xd,l);

% BOUNDARY VALUE PROBLEM: supponiamo la condizione iniziale data da
%             { -sin(2*pi*t/L) se t > 1/a
%   u0(0,t) = {
%             {       0        se t < 1/a

% Comando: u(:,1)=ub(td,a,l);




%% IMPLEMENTAZIONE SCHEMI NUMERICI DIFFERENTI
% Inizializzazione della soluzione: le righe della soluzioni rappresentano
% la soluzione al tempo t-esimo allo scorrere delle colonne in cui è
% memorizzata la soluzione alla posizione x-esima

% Inizializzazione soluzione, comando: u=zeros(Nt,Nx);
% Initial Value Problem, comando:  u(1,:)=u0(xd,l);
% Boundary Value Problem, comando: u(:,1)=ub(td,a,l);

% RISOLUZIONE CON DIVERSI SCHEMI NUMERICI
for scheme=1:4
    
    u=zeros(Nt,Nx);
    u(1,:)=u0(xd,l);
    u(:,1)=ub(td,a,l);
    
    
    %METODO DI LAX-WENDROFF
    % Fonte: Quarteroni, Sacco, Gervasio "Calcolo Scientifico con MATLAB e
    % Octave"
    if scheme==1
        % Algoritmo di Lax-Wendroff
        for n=1:Nt-1
            for j = 2:Nx-1
            u(n+1,j)=u(n,j)-CFL/2*(u(n,j+1)-u(n,j-1))+CFL^2/2*(u(n,j+1)-2*u(n,j)+u(n,j-1));
            end
        j=Nx;
        u(n+1,j)=(1-CFL^2)*u(n,j)+CFL/2*((CFL-1)*(2*u(n,j)-u(n,j-1))+(CFL+1)*u(n,j-1));
        end
        u_lw = u;
    end
    
    
    %METODO DI LAX-FRIEDRICHS
    % Fonte: Quarteroni, Sacco, Gervasio "Calcolo Scientifico con MATLAB e
    % Octave"
    if scheme==2
        % Algoritmo di Lax-Friedrichs
        for n=1:Nt-1
            for j = 2:Nx-1
            u(n+1,j)=0.5*(u(n,j+1)+u(n,j-1))-CFL/2*(u(n,j+1)-u(n,j-1));
            end
            j=Nx;
            u(n+1,j)=0.5*( (CFL-1)*(2*u(n,j)-u(n,j-1)) + (CFL+1)*(u(n,j-1)));
        end
        u_lf = u;
    end
    
      
    %METODO FIRST ORDER UPWIND
    % Fonte: Quarteroni, Sacco, Gervasio "Calcolo Scientifico con MATLAB e
    % Octave"
    if scheme==3
        % Algoritmo Upwind
        for n=1:Nt-1
            for j = 2:Nx-1
            u(n+1,j)=u(n,j)+CFL*(-u(n,j)+u(n,j-1));
            end
            j=Nx;
            u(n+1,j)=u(n+1,j-1);
        end
        u_upw = u;
    end
    
    
    %METODO II UPWIND
    % Algoritmo Walming and Beam (1976), Second Order Upwind
    % Fonte: Hirsch, "Numerical Computation of Internal and External Flows",
    % pag. 325 eq.(7.4.35)
    if scheme==4
        % Upwind I order per i primi due passi spaziali: ad
        % ogni istante di tempo calcoliamo la soluzione per j=1 (x=0) la
        % condizione iniziale, e per j=2 il primo passo spaziale successivo
        j=2;
        for n = 2:Nt-1
            u(n+1,j)=u(n,j)+CFL*(-u(n,j)+u(n,j-1));
        end

        % Algoritmo Upwind II
        for n=1:Nt-1
            for j = 3:Nx-1
            u(n+1,j)=u(n,j)-0.5*CFL*(3*u(n,j)-4*u(n,j-1)+u(n,j-2))+0.5*CFL^2*(u(n,j)-2*u(n,j-1)+u(n,j-2));
            end
        end
        u_upwII = u;
    end
        
end
%% SOLUZIONE ESATTA ALL'ISTANTE T
% Istante di tempo in cui si desidera visualizzare la soluzione: cambiando
% t_ask, si individua la riga della matrice u contenente la soluzione ad
% ogni istante di tempo corrispondente all'istante di tempo desiderato

t_ask=0.4;
ind = find(td==t_ask);


% Calcolo della soluzione esatta
% uex=@(x,t) sin(2*pi*(x-a*t)/l);
U=uex(xd,td(ind),a,l);


%% PLOT DELLE SOLUZIONI
if initial_condition == 1
    figure(1)
    subplot(4,1,1)
    plot(xd,u_lw(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-0.8 t_ask+0.8 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Lax-Wendroff')

    subplot(4,1,2)
    plot(xd,u_lf(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-0.8 t_ask+0.8 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Lax-Friedrichs')

    subplot(4,1,3)
    plot(xd,u_upw(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-0.8 t_ask+0.8 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Upwind')

    subplot(4,1,4)
    plot(xd,u_upwII(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-0.8 t_ask+0.8 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Upwind II')
elseif initial_condition == 2
    figure(1)
    subplot(4,1,1)
    plot(xd,u_lw(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-1.4 t_ask+1.4 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Lax-Wendroff')

    subplot(4,1,2)
    plot(xd,u_lf(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-1.4 t_ask+1.4 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Lax-Friedrichs')

    subplot(4,1,3)
    plot(xd,u_upw(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-1.4 t_ask+1.4 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Upwind')

    subplot(4,1,4)
    plot(xd,u_upwII(ind,:),'LineWidth',1.2), hold on
    axis([t_ask-1.4 t_ask+1.4 -1.2 1.2])
    plot(xd,U,'--','LineWidth',1)
    title('Upwind II')
end



%% INITIAL VALUE AND BOUNDARY VALUE PROBLEM FUNCTIONS
% Function per le condizioni al contorno: variando i function handle
% presenti all'interno delle function di Matlab è possibile variare il
% problema ai valori iniziali. Il BOUNDARY Value Problem viene ricostruito
% tramite la soluzione esatta della equazione di avvezione, sapendo che
%  u_ex(x,t) = u0(x-at), ponendo x=0 si ricava la condizione per il BVP.

function [u_ivp] = u0(x,L)
   u_ivp=[];
   
   f_ivp=@(x) sin(2*pi*x/L);
   
   for k=1:length(x)
       if x(k)>=-L && x(k)<=L
           u_ivp(k)= f_ivp(x(k));
       else
           u_ivp(k)=0;
       end
   end    
end

function [u_bvp] = ub(t,A,L)
   u_ivp=[];
   
   f_bvp=@(t) -sin(2*pi*A*t/L);
   
    for k=1:length(t)
        if t(k) >= 1/A
            u_bvp(k)=f_bvp(t(k));
        else
            u_bvp(k)=0;
        end
    end
end

%% SOLUZIONE ESATTA
% Function per la soluzione esatta: sapendo che u_ex(x,t) = u0(x-at), si
% pone in questa function la soluzione esatta del problema tramite una
% piecewise function, calcolata sull'intero intervallo delle x al tempo t
% indicato

function [u_exact] = uex(x,t,A,L)
   u_exact=[];
   
   f_ex=@(x) sin(2*pi*(x-A*t)/L);
   
    for k=1:length(x)
        if (x(k)-A*t)>=-L && (x(k)-A*t)<=L
            u_exact(k)=f_ex(x(k));
        elseif (x(k)-A*t)>=L && (x(k)-A*t)<=3
            u_exact(k)=0;
        end
    end
end
