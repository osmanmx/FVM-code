%Modelo Bidodominio - Método de Volúmenes Finitos (FVM)
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
%Rutina para generar la malla admisible y establecer condiciones iniciales en t=0
% clear;
xmin=0; %input(’digite el valor minimo de x xmin=’);
xmax=5; %input(’digite el valor maximo de x xmax=’);
ymin=0; %input(’digite el valor minimo de x ymin=’);
ymax=5; %input(’digite el valor maximo de x ymax=’);
n=20; %input(’digite el numero de intervalos n=’);

%vectores de coordenadas de los puntos de la malla
deltax=(xmax-xmin)/n;deltay=(ymax-ymin)/n;
coordenadasx=linspace(xmin,xmax,n+1);coordenadasy=linspace(ymin,ymax,n+1);

%coordenadas de los puntos medios de las celdas
pm(1,:)=linspace(deltax/2,xmax-deltax/2,n);
pm(2,:)=linspace(deltay/2,ymax-deltay/2,n);

%condiciones iniciales
w=0; a1=0.16875;b=1.0;

%constantes
cm=1;B=1/2000;
t=0;
dt=1e-3; %Paso en el tiempo
ka=deltax*deltay; %area de celda
for i=1:n
for j=1:n
x=linspace(coordenadasx(i),coordenadasx(i+1),20);
y=linspace(coordenadasy(i),coordenadasx(i+1),20);
vf=1 - 1./(1+exp(-50*sqrt(x.^2+y.^2)- 0.1)); %Funcion de potencial inicial
vf=trapz(y,vf);vf=vf*ones(1,20);
v(i,j,1)=trapz(x,vf)/ka;
u(i,j,1)=0.0; %Funcion potencial V_e inicial
Hf=a1*v(i,j,1)-b*w;Hf=Hf*ones(1,20);
H=trapz(y,Hf);H=H*ones(1,20);
Hk(i,j)=trapz(x,H)/ka;
wk(i,j)=w*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)- coordenadasy(j))/ka;
Iion=-gamma*(wk(i,j)-v(i,j)*(1-v(i,j))*(v(i,j)-teta));
Ion(i,j)=Iion*((coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j)))/ka;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculo de las Mi de cada celda
Me=[1/6,0;0,10/6]; Tensor Inicial
Mi=[1/24,0;0,1/12];Tensor Inicial
int=Mi*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mkli=norm(int);
int=Mi*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mlki=norm(int);
di=Mkli*Mlki*deltax/(deltax/2*Mkli+deltax/2*Mlki);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculo de Mi y Me en cada celda
int=Me*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mkle=norm(int);
int=Me*(coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j))*normal/ka;
Mlke=norm(int);
de=Mkle*Mlke*deltax/(deltax/2*Mkle+deltax/2*Mlke);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%suma de valores de celdas vecinas para V_e
%Me=[0.1,0;0,0.1];
%Mi=[0.2,0;0,0.2];
%vectores normales
%celdas   1  2  3
%%%%%%%%  8  X  4
%%%%%%%   7  6  5

%n1=[-1/sqrt(2) 1/sqrt(2)];%para la celda 1
n2=[0 1];
%n3=[1/sqrt(2) 1/sqrt(2)];%para la celda 3
n4=[1 0];
%n5=[1/sqrt(2) -1/sqrt(2)];%para la celda 5
n6=[0 -1];
%n7=[-1/sqrt(2) -1/sqrt(2)];%para la celda 8
n8=[-1 0];

%se debe calcular primero Mi y Me, luego se determina di y de, finalmente se hace la suma
%se debe iterar para cada celda y determinar: los valores de Mi y Me
%en cada celda con el programa calculo M2

suma=0;
if i==1 && j==1% esta en la esquina superior izquierda (suma 2)
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));

elseif i==1 && j==n % esta en la esquina superior derecha (suma 2)
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)- u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));

elseif i==1 % esta en los bordes superiores (suma 3)
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)-u(i,j,k-1));
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));

elseif i==n && j==1 %esta en la esquina inferior izquierda (suma 2)
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));

elseif i==n && j==n %esta en la esquina inferior derecha (suma 2)
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)-u(i,j,k-1));
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));

elseif j==1 % esta en los bordes laterales izquierdos (suma 3)
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));

elseif i==n % esta en el borde inferior (suma 3)
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)-u(i,j,k-1));
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));

elseif j==n % esta en el borde lateral derecho (suma 3)
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)-u(i,j,k-1));

else % esta en las celdas centrales (suma 4)
normal=n8’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j-1,k-1)-u(i,j,k-1));
normal=n4’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i,j+1,k-1)-u(i,j,k-1));
normal=n2’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i-1,j,k-1)-u(i,j,k-1));
normal=n6’; calculoM2;de=(Mkle+Mlke)/(deltax/2*Mkle+deltax/2*Mlke)*deltax;
suma =suma+de*deltax/deltax+(u(i+1,j,k-1)-u(i,j,k-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construccion de la matrices de coeficientes para resolver el sistema para potencial extracelular
%vectores normales
%celdas   1  2  3
%%%%%%%%  8  X  4
%%%%%%%   7  6  5

%n1=[-1/sqrt(2) 1/sqrt(2)];%para la celda 1
n2=[0 1];
%n3=[1/sqrt(2) 1/sqrt(2)];%para la celda 3
n4=[1 0];
%n5=[1/sqrt(2) -1/sqrt(2)];%para la celda 5
n6=[0 -1];
%n7=[-1/sqrt(2) -1/sqrt(2)];%para la celda 8
n8=[-1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=zeros(n^2);
%primera diagonal
cont=0;
for jj=n+1:n^2
cont=cont+1;
a(cont,jj)=di+de;
end

%segunda diagonal
cont=0;
for ii=n+1:n^2
cont=cont+1;
a(ii,cont)=di+de;
end

%diagonal principal
%bordes superiores
for jj=2:n
a(jj,jj)=-3*(di+de);a(jj,jj-1)=di+de;a(jj,jj+1)=di+de;
end
%bordes laterales izquierdos
for jj=n+1:n:n^2-1
a(jj,jj)=-3*(di+de);a(jj,jj+1)=di+de;
end
%bordes inferiores
for jj=n^2+1-n:n^2-1
if a(jj,jj)==0
a(jj,jj)=-3*(di+de);a(jj,jj-1)=di+de;a(jj,jj+1)=di+de;
end
end
%bordes laterales derechos
for jj=n:n:n^2-1
a(jj,jj)=-3*(di+de);a(jj,jj-1)=di+de;
end
%esquinas de 2
a(1,1)=-2*(di+de);a(n^2,n^2)=-2*(di+de);
a(n,n)=-2*(di+de);a(n^2-n+1,n^2-n+1)=-2*(di+de);
%restantes que corresponden a los centrales
for jj=n+1:n^2
if a(jj,jj)==0
a(jj,jj)=-4*(di+de);
end
end
%contiguos a la diagonal principal
for jj=n+1:n^2-1
if a(jj,jj)==-4*(di+de); % centrales(antes y despues)
a(jj,jj-1)=di+de;a(jj,jj+1)=di+de;
end
end
cont=0;
for ii=1:n
for jj=1:n
cont=cont+1;
if ii==1 && jj==1 % esta en la esquina superior izquierda (suma 2)

%coeficientes del vector de terminos independientes
b(cont)=ka*Iapp-di*(v(ii,jj+1,k)-v(ii,jj,k-1))-di*(v(ii+1,jj,k-1)-v(ii,jj,k-1));

elseif ii==1 && jj==n % esta en la esquina superior derecha (suma 2)
b(cont)=ka*Iapp-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1))-di*(v(ii+1,jj,k-1)-v(ii,jj,k-1));

elseif ii==1 % esta en los bordes superiores (suma 3)
b(cont)=ka*Iapp-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1))-di*(v(ii,jj+1,k-1)-v(ii,jj,k-1))-
di*(v(ii+1,jj,k-1)-v(ii,jj,k-1));

elseif ii==n && jj==1 % esta en la esquina inferior izquierda (suma 2)
b(cont)=ka*Iapp-di*(v(ii-1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii,jj+1,k-1)-v(ii,jj,k-1));

elseif ii==n && jj==n % esta en la esquina inferior derecha (suma 2)
b(cont)=ka*Iapp-di*(v(ii-1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1));

elseif jj==1 % esta en los bordes laterales izquierdos (suma 3)
b(cont)=ka*Iapp-di*(v(ki-1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii+1,jj,k-1)-v(ii,jj,k-1))-
di*(v(ii,jj+1,k-1)-v(ii,jj,k-1));

elseif ii==n % esta en el borde inferior (suma 3)
b(cont)=ka*Iapp-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1))-di*(v(ii-1,jj,k-1)-v(ii,jj,k-1))-
di*(v(ii,jj+1,k-1)-v(ii,jj,k-1));

elseif jj==n % esta en el borde lateral derecho (suma 3)
b(cont)=ka*Iapp-di*(v(ii-1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1))-
di*(v(ii+1,jj,k-1)-v(ii,jj,k-1));

else % esta en las celdas centrales (suma 4)
b(cont)=ka*Iapp-di*(v(ii-1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii,jj-1,k-1)-v(ii,jj,k-1))-
di*(v(ii+1,jj,k-1)-v(ii,jj,k-1))-di*(v(ii,jj+1,k-1)-v(ii,jj,k-1));
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Se soluciona el sistema de pentadiagonal
bt=b’;
x1=pentsolve(a,bt);
cont=0;
for ii=1:n
for jj=1:n
cont=cont+1;
u(ii,jj,k)=x1(cont);
end
end
%Condicion de flujo en la frontera
for ii=1:n
u(1,ii,k)=0;
u(n,ii,k)=0;
u(ii,1,k)=0;
u(ii,n,k)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programa principal
clear;clc
% 1-datos y calculos iniciales t=0
% 2 calcular potenciales V_{m} y V_{e}
%3-iterar t y recalcular todo
\newpage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcular el u,w,v inicial
%Calcular H,inicial integre para calcular H en cada celda
%Con H calcular
%con v y w calcular Iion e Iapp
%empezar a iterar en tiempo
%calcular los nuevos u,w,v, H, Iion e Istm
gamma=-100; %dato para ingresar 
teta=0.25; %dato para ingresar 
x0=0.5; %dato para ingresar 
y0=0.5; %dato para ingresar 

malla2; % calcula u,v,w, H, Ion en el tiempo 0
for k=2:5
t=t+dt;
if t>1e-3 && (coordenadasx(i)-x0)^2+(coordenadasy(j)-y0)^2<0.04;
Iapp=1;
else
Iapp=0;
end
for i=1:n
for j= 1:n
H=a1*v(i,j,k-1)- b*wk(i,j);
w=dt*Hk+w;sumatoria2;
v(i,j,k)=dt/(B*cm*ka)*(ka*Iapp-B*ka*Ion(i,j)+B*cm*ka/dt*v(i,j,k-1)-suma);
matrices; %Calcula el valor de u
Iion=-gamma*(wk(i,j)-v(i,j)*(1-v(i,j))*(v(i,j)-teta));
Ion(i,j)=Iion*((coordenadasx(i+1)-coordenadasx(i))*(coordenadasy(j+1)-coordenadasy(j)))/ka;
end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programa que realiza las graficas de V_{m}
xg=linspace(0,5,n);yg=xg;[xg,yg]=meshgrid(xg,yg);
for z=1:20
subplot(1,2,1)
contourf(xg,yg,v(1:n,1:n,z)) %Grafica de la superficie para V_{m}
subplot(1,2,2)
surf(xg,yg,v(1:n,1:n,z)),view(60,30) %Curvas de nivel para V_{m}.
title(z)
pause(0.1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Programa que realiza las graficas de V_{e}
xg=linspace(0,5,n);yg=xg;[xg,yg]=meshgrid(xg,yg);
for z=1:20
subplot(1,2,1)
contourf(xg,yg,u(1:n,1:n,z)) %Superficie para V_{e}
subplot(1,2,2)
surf(xg,yg,u(1:n,1:n,z)),view(60,30)%Curvas de nivel para V_{e}.
title(z)
pause(0.1)
end

