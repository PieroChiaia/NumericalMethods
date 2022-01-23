% HOMEWORK N.16
% STUDIO DI STABILITA' EQ. DIFFERENZIALE
% TRAMITE CODICE QUARTERONI

close all
clear all
clc

oldpath = path;
path(oldpath,'../Codici Quarteroni')

% Definizione del problema

tspan=[0 100];       %Dominio
y0=1;                %Condizione iniziale
lambda=-1;
f=@(t,y) lambda*y; 
m=10;                % Parametro per regolare le suddivisioni massime
                     % Suddivisioni successive in 2,4,8... sottointervali
                     
y=@(t) exp(lambda*t);%Soluzione esatta


figure(1)

% Risolviamo con Backward Euler implicito
% Condizione di stablità data da h< - Re(lambda)/ abs(lambda)^2
subplot(2,3,1)
for k=1:m
    N(k)=2^k; % Aumento numero intervalli
    hold on
    [tt,u]=beuler(f,tspan,y0,N(k));
    e_be(k)=max(abs(u-y(tt))); % Memorizzo errore massimo
    plot(tt,u);
end
title('Backward Euler Method')
subplot(2,3,4)
semilogy(N,e_be,'r','LineWidth',1);
title('BE: errore massimo')
xlabel('Num. suddivisioni') 
ylabel('Errore') 



% Risolviamo con Forward Euler esplicito
% Condizione di stablità data da h<2/abs(lambda)
subplot(2,3,2)
for k=1:m
    N(k)=2^k; % Aumento numero intervalli
    hold on
    [tt,u]=feuler(f,tspan,y0,N(k)); 
    e_fe(k)=max(abs(u-y(tt))); % Memorizzo errore massimo
    plot(tt,u);
end
title('Forward Euler Method')
subplot(2,3,5)
semilogy(N,e_fe,'r','LineWidth',1);
title('FE: errore massimo')
xlabel('Num. suddivisioni') 
ylabel('Errore') 



% Risolviamo con Crank-Nicholson
% Condizione di stablità data da Re(h*lambda)<0
subplot(2,3,3)
for k=1:m
    N(k)=2^k; % Aumento numero intervalli
    hold on
    [tt,u]=cranknic(f,tspan,y0,N(k));
    e_cn(k)=max(abs(u-y(tt))); % Memorizzo errore massimo
    plot(tt,u);
    
end
title('Crank-Nicholson Method')
subplot(2,3,6)
semilogy(N,e_cn,'r','LineWidth',1);
title('CN: errore massimo')
xlabel('Num. suddivisioni') 
ylabel('Errore') 

% I metodi BE e CN sono stabili, garantiscono una corretta valutazione
% del problema per tutte le suddivisioni scelte, grazie alle loro
% proprietà di stabilità

% Il metodo FE invece fallisce fintanto che il numero di suddivisioni
% dell'intervallo considerato non è sufficiente
