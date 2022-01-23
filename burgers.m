%% HOMEWORK 24
% EQUAZIONE DI BURGERS
% CON DIVERSI SCHEMI NUMERICI E DIVERSE CONDIZIONI INIZIALI
% 
% Risolvere l'equazione di Burgers
%     du/dt + u du/dx = 0
% con schema Lax, Lax-Wendroff e Upwind

% Definiamo la geometria del problema con gli intervalli:
% [ X_SX; X_DX ] x [ T_IN, T_FIN ]

% Le condizioni al contorno in termini spaziali-temporali sono
% INITIAL VALUE PROBLEM  --->  u(x,0) = u0
% BOUNDARY VALUE PROBLEM  -->  u(X_SX,t) = ub

% Il dominio è discretizzato in:
% Spazio --> xd = suddivisione in intervalli uguali di [ X_SX; X_DX ]
% Tempo ---> td = suddivisione in intervalli uguali di [ T_IN, T_FIN ]

% Il parametro di governo del problema sarà direttamente descritto da CFL,
% ,la velocità "u" prevista nella equazione varia nel tempo, ma deve valere
% la condizione CFL = u dt/dx < 1

% INITIAL VALUE PROBLEM e BOUNDARY VALUE PROBLEM sono due problemi
% differenti, definiti da una condizione iniziale a gradino e sinusoidale

clear all
close all
clc

%% SELEZIONE DEL PROBLEMA
% CONDIZIONE INIZIALE: SCEGLIERE 1 PER SINUSOIDALE
%                      SCEGLIERE 2 PER GRADINO
initial_condition = 1;


%% INIZIALIZZAZIONE VALORI DEL PROBLEMA, DOMINIO, E DISCRETIZZAZIONE
% Dominio del problema
t_in=0;
t_fin=3;

x_sx=0;
x_dx=4*pi;

% Parametri numerici del problema
a=1; %Velocità variabile
N=1000;
x=linspace(x_sx,x_dx,N+1);
h=(x_dx-x_sx)/N;
CFL=0.6;    
dt=1/a*CFL*h;

% Discretizzazione (iniziale) del dominio
Nt=fix((t_fin-t_in)/dt+1);
td=linspace(t_in,t_fin,Nt);
Nx=fix((x_dx-x_sx)/h+1);
xd=linspace(x_sx,x_dx,Nx);


%% RISOLUZIONE CON DIVERSI SCHEMI NUMERICI
% Per ogni istante di tempo, è come risolvere una serie successiva di
% equazioni di avvezione, con il parametro "a" che varia nel tempo: viene
% quindi ricalcolato ad ogni step il dt adeguato, sulla base del CFL
% imposto a priori, lo step temporale che garantisce la condizione di
% stabilità del sistema

for scheme=1:4
    
    u=zeros(Nt,Nx);
    u(1,:)=u0(xd,initial_condition);
    u(:,1)=ub(td,a,initial_condition);
    
    %METODO DI LAX
    % Fonte: Tannehill "Computational Fluid Mechanics and Heat Transfer", pag. 181 eq.(4.151)
    if scheme==1
        
        t_lim=0;
        n=1;
        while t_lim<t_fin

            for j = 2:Nx-1
                u(n+1,j)=0.5*(u(n,j+1)+u(n,j-1))-dt/(4*h)*((u(n,j+1))^2-(u(n,j-1))^2);
            end
            j=Nx;
            u(n+1,j)=u(n+1,j-1);
           
            t_lim=t_lim+dt;
            n=n+1;
        end
        u_LAX = u;
    end
    
    
    %METODO DI LAX-WENDROFF
    % Fonte: Tannehill "Computational Fluid Mechanics and Heat Transfer", pag. 185 eq.(4.162)
    if scheme==2
        t_lim=0;
        n=1;

        while t_lim<t_fin
            for j = 2:Nx-1
                Fn_j   = 0.5*(u(n,j))^2;
                Fn_jm1 = 0.5*(u(n,j-1))^2;
                Fn_jp1 = 0.5*(u(n,j+1))^2;

                A_jp05 = 0.5*(u(n,j)+u(n,j+1));
                A_jm05 = 0.5*(u(n,j)+u(n,j-1));

                u(n+1,j)=u(n,j)-dt/h*0.5*(Fn_jp1-Fn_jm1)+0.5*(dt/h)^2*(A_jp05*(Fn_jp1-Fn_j)-A_jm05*(Fn_j-Fn_jm1));
            end
        j=Nx;
        u(n+1,j)=u(n+1,j-1);
   
        t_lim=t_lim+dt;
        n=n+1;
        end
        u_LW = u;
    end
    
    
    % METODO PRED/CORR MACCORMACK
    % Fonte: Tannehill "Computational Fluid Mechanics and Heat Transfer", pag. 187 eq.(4.167)
    if scheme==3
        t_lim=0;
        n=1;
        while t_lim<t_fin
    
        for j = 2:Nx-1
            Fn_j   = 0.5*(u(n,j))^2;
            Fn_jp1 = 0.5*(u(n,j+1))^2;
            Fn_jm1 = 0.5*(u(n,j-1))^2;
        
            pred_u=u(n,j)-dt/h*(Fn_jp1-Fn_j);
            pred_Fj= 0.5*(pred_u)^2;
            pred_um1=u(n,j-1)-dt/h*(Fn_j-Fn_jm1);
            pred_Fjm1=.5*(pred_um1)^2;
        
            u(n+1,j)=0.5*(u(n,j)+pred_u-dt/h*(pred_Fj-pred_Fjm1));       
        end
        j=Nx;
        u(n+1,j)=u(n+1,j-1);
        u_MC = u;

        t_lim=t_lim+dt;
        n=n+1;
        end
        
        u_MC = u;

    end
    
    
    % FIRST ORDER UPWIND
    % Fonte: Quarteroni, Sacco, Gervasio "Calcolo Scientifico con MATLAB e Octave"
  
    if scheme==4
        t_lim=0;
        n=1;

        while t_lim<t_fin
           for j = 2:Nx-1
               
               if u(n,j) >= 0
                    u(n+1,j)=u(n,j)-dt/h*(flux(u(n,j))-flux(u(n,j-1)));
               else
                   u(n+1,j)=u(n,j)-dt/h*(flux(u(n,j+1))-flux(u(n,j)));
               end
               
           end

           j=Nx;
           u(n+1,j)=u(n,j)-dt/h*u(n,j)*(u(n,j)-u(n,j-1));
 
           t_lim=t_lim+dt;
           n=n+1;
        end
        u_UPWI = u;
    end
    

end

%% LIVE-PLOT DELLA SOLUZIONE
% Set-up automatico degli assi per la visualizzazione 

if initial_condition==1
    % Plot della soluzione ad un istante qualsiasi per il confronto tra
    % metodi numerici per la risoluzione
    figure(1)
    subplot(4,2,[2 4 6 8])
    plot(xd,u_LAX(200,:),'LineWidth',1.2), hold on
    plot(xd,u_LW(200,:),'LineWidth',1.2), hold on
    plot(xd,u_MC(200,:),'LineWidth',1.2), hold on
    plot(xd,u_UPWI(200,:),'LineWidth',1.2), hold on
    axis([2.8 3.5 -1.5 2])
    title('Comparison between solutions')
    legend('Lax','Lax-Wendroff','Mac-Cormack','Upwind I')
    hold off
    
    % Live-plot della soluzione
    for m=1:length(td)
        figure(1)

        subplot(4,2,1)
        plot(xd,u_LAX(m,:),'LineWidth',1.2)
        axis([0 4*pi -1.2 1.2])
        title('Schema Lax')

        subplot(4,2,3)
        plot(xd,u_LW(m,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.2)
        axis([0 4*pi -1.2 1.2])
        title('Schema Lax-Wendroff')

        subplot(4,2,5)
        plot(xd,u_MC(m,:),'Color',[0.9290,0.6940, 0.1250],'LineWidth',1.2)
        axis([0 4*pi -1.2 1.2])
        title('Schema Mac-Cormack')

        subplot(4,2,7)
        plot(xd,u_UPWI(m,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.2)
        axis([0 4*pi -1.2 1.2])
        title('Schema Upwind I')
    end

end

if initial_condition==2
    % Plot della soluzione ad un istante qualsiasi per il confronto tra
    % metodi numerici per la risoluzione
    figure(1)
    subplot(4,2,[2 4 6 8])
    plot(xd,u_LAX(100,:),'LineWidth',1.2), hold on
    plot(xd,u_LW(100,:),'LineWidth',1.2), hold on
    plot(xd,u_MC(100,:),'LineWidth',1.2), hold on
    plot(xd,u_UPWI(100,:),'LineWidth',1.2), hold on
    axis([0.15 0.6 -0.2 1.2])
    title('Comparison between solutions')
    legend('Lax','Lax-Wendroff','Mac-Cormack','Upwind I')
    hold off
    
    % Live-plot della soluzione
    for m=1:length(td)
    figure(1)
    
    subplot(4,2,1)
    plot(xd,u_LAX(m,:),'LineWidth',1.2)
    axis([0 0.7 -0.2 1.2])
    title('Schema Lax')
    
    subplot(4,2,3)
    plot(xd,u_LW(m,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.2)
    axis([0 0.7 -0.2 1.2])
    title('Schema Lax-Wendroff')
    
    subplot(4,2,5)
    plot(xd,u_MC(m,:),'Color',[0.9290,0.6940, 0.1250],'LineWidth',1.2)
    axis([0 0.7 -0.2 1.2])
    title('Schema Mac-Cormack')
    
    subplot(4,2,7)
    plot(xd,u_UPWI(m,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.2)
    axis([0 0.7 -0.2 1.2])
    title('Schema Upwind I')
    
    end
end

% figure(2)
% for m=1:length(td)
%     plot3(m*ones(1,length(xd)),xd,u_LW(m,:),'k'), hold on
%     plot3((m+1)*ones(1,length(xd)),xd,u_LW((m+1),:),'b'), hold off
%     view([82,20])
%     axis([0 length(td) 0 4*pi -1 1])
%     pause(0.01)
% end


%% INITIAL VALUE AND BOUNDARY VALUE PROBLEM FUNCTIONS
% Function per le condizioni al contorno: variando i function handle
% presenti all'interno delle function di Matlab è possibile variare il
% problema ai valori iniziali. Il BOUNDARY Value Problem viene ricostruito
% tramite la soluzione esatta della equazione di Burgers, sapendo che
%  u_ex(x,t) = u0(x-at), ponendo x=0 si ricava la condizione per il BVP.

function [u_ivp] = u0(x,in)
    % Cond. iniziale sinusoidale == 1
    u_ivp=[];
    for k=1:length(x)
        
    if in == 1
        f_ivp=@(v) sin(v);
        u_ivp=f_ivp(x);
    
    % Cond. iniziale gradino == 2
    elseif in==2
        if x(k)<0
            u_ivp(k)=1;
        else
            u_ivp(k)=0;
        end
    end
    
    end

end

function [u_bvp] = ub(t,A,in)

    % Cond. iniziale sinusoidale == 1
    for k=1:length(t)
        
    if in == 1
        u_bvp(k)=0;
    
    % Cond. iniziale gradino == 2
    elseif in==2
        if -A*t(k)<0
            u_bvp(k)=1;
        else
            u_bvp(k)=0;
        end
    end
    
    end

end

function [flux] = flux(t)
    flux=t^2/2;
end