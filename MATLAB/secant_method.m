% HOMEWORK N.1
% METODO DELLE SECANTI

% IL CODICE E' STATO VERIFICATO CON ESERCIZI
% DI CALCOLO NUMERICO DI CUI ERA NOTA LA SOLUZIONE

% f=@(x) x-log(2-x+x.^2) IN [0,1], zero=0.561
% f=@(x) x.*exp(3*x)-1-x IN [-2,-1], zero=-1.045
% f=@(x) log(2+sin(x))-x IN [1,2], zero=1.054
% f=@(x) 2*sin(x)+x-1 IN [0,1], zero=0.337
% f=@(x) (sin(x)).^2-x-1 IN [-1,0], zero=-0.641

% E SUCCESSIVAMENTE UTILIZZATO NELLA RISOLUZIONE
% DELLA EQUAZIONE DI STATO ASSEGNATA

close all
clear all
clc

%Dati e costanti del problema
p=3.5e7;
T=300.;
N=1000.;
a1=0.401;
b1=42.7e-6;
kb=1.3806503e-23;

% Inizializzo funzione, dati e intervallo
f=@(x) (p+a1*(N^2./x.^2)).*(x-N*b1)-kb*N*T
a=0.03;
b=0.07;

x=linspace(a,b,300);
fx=f(x);
%plot(x,fx)

tol = 1.e-12;
k=0;
kmax=100;

% Confronto con metodo di bisezione, con
% funzione già implementata
[zero,res,niter]=bisection(f,a,b,tol,100);
% Visualizzo a video il risultato
fprintf('Bisezione: Iter = %5.0f \t Zero = %10.7f \t Errore = %6.3e \n', niter, zero, res);

% Confronto con metodo di Newton, con
% funzione già implementata
df=@(x) p-a1*N^2./x.^2 + 2*a1*N^3*b1./x^(3);
[zero,res,niter]=newton(f,df,0.01,tol,100);
% Visualizzo a video il risultato
fprintf('Newton:    Iter = %5.0f \t Zero = %10.7f \t Errore = %6.3e \n', niter, zero, res);


% Verifico applicabilità del metodo
if fx(1)*fx(end) > 0
    disp('Metodo non applicabile');
elseif fx(1) == 0
    zero = x(1);
    fprintf('Iter = %5.0f \t Zero = %10.7f', k, zero)
elseif fx(end) == 0
    zero = x(end);
    fprintf('Iter = %5.0f \t Zero = %10.7f', k, zero)
else
    % IMPLEMENTO METODO
    fa=f(a); fb= f(b);
    I=b-a;
    while I>=tol && k<kmax
        % Ricalcolo la funzione nei nuovi estremi
        fa=f(a);
        fb=f(b);
        % Calcolo lo zero della retta secante
        c=b-(b-a)/(fb-fa)*fb;
        fc=f(c);
    
        % Costruisco il nuovo intervallo di iterazione
        %aggiornando gli estremi a e b
        if fa*fc < 0
            b=c;
        else
            a=c;
        end
        I=abs(b-a)/2;
        k=k+1;
        res=f(c);
    end
end

% Visualizzo a video il risultato
fprintf('Secanti:   Iter = %5.0f \t Zero = %10.7f \t Errore = %6.3e\n', k, c, res);

% Dalla esecuzione risulta:
%Bisezione: Iter =    35 	 Zero =  0.0427000 	 Errore = -1.009e-04 
%Newton:    Iter =    10 	 Zero =  0.0427000 	 Errore = -4.142e-18 
%Secanti:   Iter =    49 	 Zero =  0.0427000 	 Errore = -4.142e-18
% Il metodo delle secanti risulta accurato ma, in questo caso, il più
% lento

    
