%Modelo Monodominio - Método de Volúmenes Finitos (FVM)
%Autor: Osman Villanueva García
%[Licencia GPL versión 3.0](https://github.com/osmanmx/FVM-code/blob/main/LICENSE)
% /***************************************************************************                                        *
%  *   Copyright (C) 2013 by Osman Villanueva Garcia                         *
%  *   Email osman@educart.org                                               *
%  *                                                                         *
%  *   This program is free software: you can redistribute it and/or modify  *
%  *   it under the terms of the GNU General Public License as published by  *
%  *   the Free Software Foundation, either version 3 of the License, or     *
%  *   (at your option) any later version.                                   *
%  *                                                                         *
%  *   This program is distributed in the hope that it will be useful,       *
%  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%  *   GNU General Public License for more details.                          *
%  *                                                                         *
%  *   You should have received a copy of the GNU General Public License     *
%  *   along with this program.  If not, see <https://www.gnu.org/licenses/> *
%  ***************************************************************************/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rutina para generar la malla y establecer condiciones iniciales en t=0
% clear;
xmin=0; %input(’digite el valor minimo de x xmin=’);
xmax=5; %input(’digite el valor maximo de x xmax=’);
ymin=0; %input(’digite el valor minimo de x ymin=’);
ymax=5; %input(’digite el valor maximo de x ymax=’);
n=5; %input(’digite el numero de intervalos n=’);

%vectores de coordenadas de los puntos de la malla
deltax=(xmax-xmin)/n;deltay=(ymax-ymin)/n;
coordenadasx=linspace(xmin,xmax,n+1);coordenadasy=linspace(ymin,ymax,n+1);

%coordenadas de los puntos medios de las celdas
pm(1,:)=linspace(deltax/2,xmax-deltax/2,n);
pm(2,:)=linspace(deltay/2,ymax-deltay/2,n);

%condiciones iniciales
w=0;
t=0;

%constantes
a=0.16875;b=1.0;cm=1.0;B=1.0;
dt=2e-3; %Paso en el tiempo
ka=deltax*deltay; %area de celda
for i=1:n
for j=1:n
x=linspace(coordenadasx(i),coordenadasx(i+1),20);
y=linspace(coordenadasy(i),coordenadasx(i+1),20);
vf=1 - 1./(1+exp(-50*sqrt(x.^2+y.^2)- 0.1)); %Funcion de potencial V_{m} inicial
v=trapz(y,vf);v=v*ones(1,20);
v(i,j,1)=trapz(x,v)/ka;
Hf=a*v(i,j,1)-b*w;Hf=Hf*ones(1,20); %Funcion de H inicial
H=trapz(y,Hf);H=H*ones(1,20);
Hk(i,j)=trapz(x,H)/ka;
wk(i,j)=w*(coordenadasx(i+1)- coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))/ka;
Iion=-gamma*(wk(i,j)-v(i,j)*(1-v(i,j))*(v(i,j)-teta));
Ion(i,j)=Iion*((coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j)))/ka;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculo de Mi de cada celda
int=Mi*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mkl=norm(int);
int=Mi*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mlk=norm(int);
d=Mkl*Mlk*deltax/(deltax/2*Mkl+deltax/2*Mlk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo de la Suma de Datos de celdas vecinas
lamnda=-0.99; %Valor de la conductividad extracelular
Mi=[1,0;0,1]*lamnda; Tensor de conductividad extracelular
%vectores normales
%celdas   1  2  3
%%%%%%%   8  X  4
%%%%%%%   7  6  5

%n1=[-1/sqrt(2) 1/sqrt(2)];%para la celda 1
n2=[0 1];
%n3=[1/sqrt(2) 1/sqrt(2)];%para la celda 3
n4=[1 0];
%n5=[1/sqrt(2) -1/sqrt(2)];%para la celda 5
n6=[0 -1];
%n7=[-1/sqrt(2) -1/sqrt(2)];%para la celda 8
n8=[-1 0];

suma=0;

if i==1 && j==1 % esta en la esquina superior izquierda (suma 2)
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
elseif i==1 && j==n % esta en la esquina superior derecha (suma 2)
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
elseif i==1 % esta en los bordes superiores (suma 3)
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i,j));
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
elseif i==n && j==1 % esta en la esquina inferior izquierda (suma 2)
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i-1,j)-v(i,j));
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
elseif j==1 % esta en los bordes laterales izquierdos (suma 3)
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i-1,j)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
elseif i==n % esta en el borde inferior (suma 3)
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i,j));
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i-1,j)-v(i,j));
elseif i==n && j==n % esta en la esquina inferior derecha (suma 2)
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i-1,j));
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
elseif j==n % esta en el borde lateral derecho (suma 3)
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i-1,j)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i,j));
else % esta en las celdas centrales (suma 4)
normal=n8’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j-1)-v(i,j));
normal=n4’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i,j+1)-v(i,j));
normal=n2’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i-1,j)-v(i,j));
normal=n6’; calculoM;d=(Mkl+Mlk)/(deltax/2*Mkl+deltax/2*Mlk)*deltax;
suma =suma+1/(1+gamma)*d*deltax/deltax+(v(i+1,j)-v(i,j));
\end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programa principal
clear;clc
% 1-datos y calculos iniciales t=0
% 2-calcular potenciales
% 3-iterar t y recalcular todo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcular el w inicial
%Calcular H,inicial integre para calcular H en cada celda
%Con H calcular v inicial
%con v y w calcular Iion e Iapp
%empezar a iterar en tiempo
%calcular los nuevos v, w, H, Iion e Istm
gamma=-100; %dato inicial para ingresar
teta=0.25; %dato inicial para ingresar
x0=2.5; %dato inicial para ingresar
y0=2.5; %dato inicial para ingresar

%producto=gamma/(1+gamma)*Istm
% 1: Programa malla
malla; % calcula v, H, Hk y w
for k=2:20
t=t+dt;
for i=1:n
for j= 1:n
if t>4e-3 && (coordenadasx(i)-x0)^2+(coordenadasy(j)-y0)^2<0.04;
producto=1; %disp(’ya’)
else
producto=0; %es (gamma/(gamma+1)*Istm)
end
%Calculo de V_{m}
sumatoria;
v(i,j,k)=(dt/(B*cm*ka))*(producto*ka-B*ka*Ion(i,j)+(B*cm*ka*v(i,j,k-1))/dt - suma);
%Calculo de w
wk(i,j)=dt*Hk(i,j)+wk(i,j);
%Calculo de Hk
H=a*v(i,j)- b*wk(i,j);
Hk(i,j)=H*((coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j)))/ka;
%Calculo Iion
Iion=-gamma*(wk(i,j)-v(i,j)*(1-v(i,j))*(v(i,j)-teta));
Ion(i,j)=Iion*((coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j)))/ka;
end
end
end

%Condicion de frontera en los bordes externos, es decir el flujo es cero.
for i=1:n
v(1,i,k)=0;
v(n,i,k)=0;
v(i,1,k)=0;
v(i,n,k)=0;
end

\newpage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programa que realiza las graficas
xg=linspace(0,5,n);yg=xg;[xg,yg]=meshgrid(xg,yg);
for z=1:20
subplot(1,2,1)
contourf(xg,yg,v(1:n,1:n,z)), %view(60,30) %Crea las curvas de nivel
subplot(1,2,2)
surf(xg,yg,v(1:n,1:n,z)),view(60,30)
%Crea la superficie V
title(z)
pause(0.3) %velocidad de transicion de las graficas.
end

