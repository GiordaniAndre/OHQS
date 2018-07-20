hbarra = 7;
m      = .1;
omega  = 1;
pos    = 1;
L      = 9

%Dimensões da função laplaciana

n = 100
x = 0:n-1

volta      = (-2*diag(ones(n,1),0) + diag(ones((n-1),1),1) ...
           + diag(ones((n-1),1),-1));

volta(1,1) = 0; 
volta(2,1) = 0; 
volta(1,2) = 0;

%estacionário, logo f(0) = 0

volta(n,n-1) = 0; 
volta(n-1,n) = 0; 
volta(n,n)   = 0;

%estacionário, logo f(L) = 0

%hamiltoniano

%definir hamiltoniano para o sistema. Oscilador Harmonico Simples

%H|p> = E|p>
H=(-(hbarra.^2)/(2.*m).*(volta.^2) + 1/2*m*(omega^2)*(pos^2));

%encontrar as soluções de H
eig(H); %encontra autosoluções para a matriz hamiltoniana

%autovalores e autovetores
[V,E] = eig(H)
V; %autoestados autovetores
E; %autoenergias autovalores

%E é uma matriz com seus autovalores na diagonal
%diagonalizar para obter as energias
En=diag(E)


FUNCFO = zeros(100,1)
%função densidade de probabilidade para os n autoestados
%estado inicial t=0 a t=2

for n=1:5 %5 primeiros autoestados do oscilador harmonico quântico
    figure
    V(:,n)=V(:,n)/max(abs(V(:,n)))
    for t=0:0.1:2
        onda=cos(t*En(n)/hbarra); %calcula a função de onda
        fo=(onda.*(V(:,n)));
        FUNCFO=FUNCFO + fo
        plot(fo.^2) %densidade de probabilidade
        axis([0,100,0,1]) %eixos variando de 0 a 100 e de -.5 a .5
        pause(0.1)
    end
end


%sem efeitos de entrelaçamento
