%% ATTRATTORE DI LORENZ
% SIMULARE L'ATTRATTORE DI LORENZ E OSSERVARE LA SENSIBILITA ALLE
% CONDIZIONI INIZIALI DEL SISTEMA DINAMICO

%  { x' = sigma( y - x )
%  { y' = rho x - y -xz
%  { z' = -beta z + xy

%% DEFINIZIONE DEL SISTEMA DI LORENZ: PARAMETRI E DISCRETIZZAZIONE

clear all
close all
clc

% Dominio del problema
t_in=0;
t_fin=100;
h=0.010;
tspan= t_in:h:t_fin;


% Implementazione Sistema di ODE Lorenz Attractor
sigma=10;
rho=28;
beta=8/3;

lorenz=@(t,x) [sigma*(x(2)-x(1)); rho*x(1)-x(2)-x(1)*x(3);-beta*x(3)+x(1)*x(2)];


%% SOLUZIONI E TRAIETTORIE
% Si genera una condizione iniziale random, ogni coordinata può assumere
% valori casuali tra -10 e 10
x0=fix(randi([-10 10],1,3));


% Risoluzione con ODE45, traiettoria non perturbata, si memorizza la
% soluzione nel vettore [ x1, y1, z1 ] per comodità di successive analisi
% (plot, esponente di Lyapunov)
[t,x]=ode45(lorenz,tspan,x0);

x1=x(:,1)';
y1=x(:,2)';
z1=x(:,3)';


% Risoluzione con ODE45, traiettoria perturbata
% Si perturba il sistema perturbando una qualsiasi componente della
% condizione iniziale e si memorizza la soluzione nel vettore [ x2, y2, z2 ]
% per comodità di successive analisi (plot, esponente di Lyapunov)
eps_x=0;
eps_y=0;
eps_z=1e-9;

x0=[x0(1)+eps_x,x0(2)+eps_y,x0(3)+eps_z];
[t_p,x_p]=ode45(lorenz,tspan,x0);

x2=x_p(:,1)';
y2=x_p(:,2)';
z2=x_p(:,3)';



%% DISTANZA TRA LE TRAIETTORIE
% Stima della distanza, istante per istante, delle traiettorie
d = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
figure(1)
subplot(4,2,8);
set(gcf,'color','w');
semilogy(t,d,'k','LineWidth',1)
axis([t_in t_fin 0 1e3])
title('Distance between trajectories during time');
xlabel('t') 
ylabel('d')

%% PLOT 
% Plot della prima coordinata, si/no perturbata
figure(1)
subplot(4,2,2);
plot(t,x1), hold on
plot(t,x2), hold off
title('Coordinate: x(t)')
xlabel('t') 
ylabel('x')

% Plot della seconda coordinata, si/no perturbata
figure(1)
subplot(4,2,4);
plot(t,y1), hold on
plot(t,y2), hold off
title('Coordinate: y(t)')
xlabel('t') 
ylabel('y')

% Plot della terza coordinata, si/noperturbata
figure(1)
subplot(4,2,6);
plot(t,z1), hold on
plot(t,z2), hold off
title('Coordinate: z(t)')
xlabel('t') 
ylabel('z')


% Dai seguenti grafici si osserva che, da un certo istante di tempo in poi,
% non è più possibile predire la dinamica del sistema in quanto, anche per
% una infinitesima perturbazione delle condizioni iniziali, la dinamica e la
% evoluzione del sistema è completamente differente


% Plot dell'attrattore nello spazio
figure(1)
subplot(4,2,[1 3 5 7]);
plot3(x1,y1,z1)
view([50 30]);
title('Lorenz Attractor')
xlabel('x') 
ylabel('y')
zlabel('z')





%% LIVE PLOT DELL'ATTRATTORE
% DECOMMENTARE IL CODICE PER OSSERVARE IL GRAFICO LIVE
% figure(2)
% for i = 1:length(x(:,1))
%     set(gcf,'color','w');
%     plot(x(1:i,1),x(1:i,2),'-r'), hold on
%     pbaspect([1 1 1])
%     plot(x_p(1:i,1),x_p(1:i,2),'-b');
%     pbaspect([1 1 1])
%     xlim([-35 35]);
%     ylim([-35 35]);
%     drawnow
% end
% title('Live-plot in the x-y plane')
% xlabel('x') 
% ylabel('y')
