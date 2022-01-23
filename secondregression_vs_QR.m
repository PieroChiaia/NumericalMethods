% HOMEWORK N.12
% FATTORIZZAZIONE QR / REGRESSIONE SECOND'ORDINE

close all
clear all
clc

% Definisco il problema stress-deformazione imponendo 
% un problema ai minimi quadrati ordine 2
% eps(sigma)= a2*sigma^2 + a1*sigma + a0

% Il vettore delle incognite sarà il vettore a=(a2,a1,a0)^T
% La matrice A sarà data dai coefficienti di a2,a1,a0 imponendo
% la condizione di interpolazione per ogni valore di sigma

sigma=[0,0.06,0.14,0.25,0.31,0.47,0.60,0.70];
eps=[0; 0.08; 0.14; 0.2; 0.23; 0.25; 0.28; 0.29];

% Definisco la matrice dei coefficienti: prima colonna, sigma^2, seconda
% colonna sigma^1, terza colonna sigma^0=1, riga per riga
for i=1:length(sigma)
    for j=1:3
        A(i,j)=sigma(i)^(3-j);
    end
end

% Fattorizzazione QR
% Risolviamo il sistema sovradeterminato A a = eps
[Q,R]=qr(A);
Qt=Q(:,1:3); Rt=R(1:3,:); 
a = Rt \ (Qt'*eps);
fprintf('Fattorizzazione QR:\n    a2 = %12.8f\n    a1 = %12.8f\n    a0 = %12.8f\n', a(1), a(2), a(3));

% Calcoliamo i valori nei nodi tramite i coefficienti
% ottenuti grazie  fattorizzazione QR
eps_app=A*a;



% Regressione al secondo ordine
% I nodi sono i valori sperimentali (sigma,eps)
c=polyfit(sigma,eps,2);
fprintf('Regressione al second ordine:\n    a2 = %12.8f\n    a1 = %12.8f\n    a0 = %12.8f\n', c(1), c(2), c(3));
x=linspace(0,sigma(8),100); % discretizzo l'intervallo delle sigma considerato
p=polyval(c,x); % Valori della curva di regressione

% La fattorizzazione QR restituisce gli stessi coefficienti
% interpolanti a2,a1,a0 di quelli ottenibili tramite una 
% interpolazione ai minimi quadrati di grado 2 sugli stessi nodi


% Plot valori sperimentali
plot(sigma,eps,'^','LineWidth',1),
hold on
% Plot QR
plot(sigma,eps_app,'--o','LineWidth',1)
hold on
% Plot regressione II ordine
plot(x,p,'LineWidth',1.2)
title('QRfact vs Second Order Regression')
legend('Valori sperimentali','QR','Lagrange2')

% Nei nodi, i valori di eps ottenuti attraverso stima con QR
% coinciono con quelli ottenuti tramite interpolazione al second'ordine
