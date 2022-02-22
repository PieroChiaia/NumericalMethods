%% HOMEWORK 21
% OSCILLATORE DI VAN DER POL
% RISULUZIONE CON DIVERSI SCHEMI NUMERICI

% L'EQUAZIONE DI VAN DER POL E'
%     d^2 x               dx
%   --------- + mu(1-x^2)---- + x = 0
%     dt^2                dt

%% DEFINIZIONE DEL PROBLEMA E DEL SISTEMA DIFFERENZIALE EQUIVALENTE
clear all
clear vars
clc

% Dominio del problema
tspan=[0,100];
x0=[1.5,1];
tstart=cputime;

% Implementazione Sistema di ODE Van Der Pol
mu=3;
vdp=@(t,y) [y(2); mu*y(2)*(1-(y(1))^2)-y(1)];



%% RISOLUZIONE CON DIVERSI SCHEMI NUMERICI
% Risoluzione con ODE23
[t_ode23,y_ode23]=ode23(vdp,tspan,x0);
tend=cputime-tstart;
nstep_ode23=size(t_ode23,1);

% Risoluzione con ODE45
[t_ode45,y_ode45]=ode45(vdp,tspan,x0);
tend=cputime-tstart;
nstep_ode45=size(t_ode45,1);

% Risoluzione con ODE23S
[t_ode23s,y_ode23s]=ode23s(vdp,tspan,x0);
tend=cputime-tstart;
nstep_ode23s=size(t_ode23s,1);

% Risoluzione con ODE15s
[t_ode15s,y_ode15s]=ode15s(vdp,tspan,x0);
tend=cputime-tstart;
nstep_ode15s=size(t_ode15s,1);


%% PLOT DELLO SPAZIO DELLE FASI
figure(1)
subplot(2,2,[1,3]);
plot(y_ode23(:,1),y_ode23(:,2)), hold on
plot(y_ode45(:,1),y_ode45(:,2)), hold on
plot(y_ode23s(:,1),y_ode23s(:,2)), hold on
plot(y_ode15s(:,1),y_ode15s(:,2))
title('Phase space')
xlabel('x') 
ylabel('y') 

%% Plot delle componenti della soluzione
figure(1)
subplot(2,2,2);
plot(t_ode23,y_ode23(:,1)), hold on
plot(t_ode45,y_ode45(:,1)), hold on
plot(t_ode23s,y_ode23s(:,1)), hold on
plot(t_ode15s,y_ode15s(:,1))
title('Time plot: x')
xlabel('t') 
ylabel('x') 

%% Plot delle componenti della soluzione
figure(1)
subplot(2,2,4);
plot(t_ode23,y_ode23(:,2),'r'), hold on
plot(t_ode45,y_ode45(:,2),'b'), hold on
plot(t_ode23s,y_ode23s(:,2),'k'), hold on
plot(t_ode15s,y_ode15s(:,2),'m')
title('Time plot: y')
xlabel('t') 
ylabel('y') 


% Non è possibile individuare evidenti differenze nella soluzione ottenuta
% grazie a diversi schemi numerici: il sistema è ben posto.
