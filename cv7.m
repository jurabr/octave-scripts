%
% Linearni stabilita  rovinneho ramu
%

clear; clc

% Dimenze ulohy (2 stupne volnosti:u,v, 2 uzly na prutu):
ndof = 3 ;
puzlu= 2 ;  

% ZADANI - material:
E=210e9 ;
bx=0.01;
hx=0.02;
A=bx*hx;
I=(1/12)*(bx)*(hx^3);

% ZADANI - geometrie:
uzly=[
0 0 ;
0 1 ;
0 2 ;
0 3 ;
1 3 ;
2 3 ;
3 3  ;
4 3 ;
5 3 ; 
5 2 ;
5 1 ;
5 0];

pruty=[
1 2 ;
2 3 ;
3 4 ;
4 5 ;
5 6 ;
6 7 ;
7 8 ;
8 9 ;
9 10;
10 11;
11 12 
];

nuzlu=size(uzly,1);
nprutu=size(pruty,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
12 1 0 ;
12 2 0 ];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
6 2 -1000 ;
9 1 -1000
];
nsil=size(sily,1);

% VYPOCET: --------------------------------

% Vypocet kodovych cisel:
kcis=zeros(nprutu,ndof*puzlu);
for i=1:nprutu
  for j=1:puzlu
    for k=1:ndof
      kcis(i, ((j-1)*ndof+k))=(pruty(i,j)-1)*ndof+k ;
    end
  end
end

% Nulovani matic a vektoru:
velikost = nuzlu*ndof;

K = zeros(velikost);
Kg = zeros(velikost);
u = zeros(velikost,1);
F = zeros(velikost,1);

% Sestaveni a lokalizace matic tuhosti:
for i=1:nprutu
  Keg=zeros(puzlu*ndof);

	% Delka prutu
	dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
	dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
	L = sqrt(dx2 + dy2);

	% Matice tuhosti:
  Keg=[ E*A/L 0 0 -E*A/L 0 0 ;
        0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2 ;
        0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L ;
        -E*A/L 0  0 E*A/L 0 0 ;
        0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2 ;
        0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L ];

	% Transformace:
	s = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2))/L;
	c = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1))/L;
  T = [ c  s  0  0  0  0 ; 
       -s  c  0  0  0  0 ;
        0  0  1  0  0  0 ; 
        0  0  0  c  s  0 ; 
        0  0  0 -s  c  0 ;
        0  0  0  0  0  1 ]  ;
  Ke = T' * Keg * T;
  
	% Samotna lokalizace:
  for j=1:(puzlu*ndof)
    for k=1:(puzlu*ndof)
      K(kcis(i,j),kcis(i,k)) = K(kcis(i,j),kcis(i,k))+Ke(j,k);
    end
  end
end

% Zatizeni:
for i=1:nsil
	iuz=sily(i,1);
	ismer=sily(i,2);
	pos = (ndof*(iuz-1))+ismer;
  F(pos) = F(pos)+sily(i,3);
end

% Podpory:
for i=1:npodpor
	iuz=podpory(i,1);
	ismer=podpory(i,2);
	pos = (ndof*(iuz-1))+ismer;
  u(pos) = u(pos)+podpory(i,3);
	for j=1:velikost
    K(pos,j) = 0.0 ;
    K(j,pos) = 0.0 ;
  end
	K(pos,pos) = 1.0 ;
end

% RESENI soustavy rovnic:
u=K\F;

% Vypocet vysledku na prutech:
for i=1:nprutu
	ue  = zeros(puzlu*ndof,1);
	uel = zeros(puzlu*ndof,1);

	% Ziskani lokalnich vektoru deformaci:
	for j=1:(puzlu*ndof)
		ue(j) = u(kcis(i,j));
	end

	% Transformace:
	dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
	dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
	L = sqrt(dx2 + dy2);

	s = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2))/L;
	c = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1))/L;
  T = [ c  s  0  0  0  0 ; 
       -s  c  0  0  0  0 ;
        0  0  1  0  0  0 ; 
        0  0  0  c  s  0 ; 
        0  0  0 -s  c  0 ;
        0  0  0  0  0  1 ]  ;
  uel = T * ue;

  % Matice tuhosti (jen lokalni):
  Kelb=[ E*A/L 0 0 -E*A/L 0 0 ;
        0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2 ;
        0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L ;
        -E*A/L 0  0 E*A/L 0 0 ;
        0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2 ;
        0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L ];

  % sily v prutech:
  Fe = Kelb * uel;
  %CV9: ziskani normalove sily:
  N=Fe(1);
  %CV9: geometricka matice tuhosti v lok. sour.:
  Kgl=(N/L)*[0 0   0    0 0    0 ;
       0 6/5 L/10 0 -6/5  L/10 ;
       0 L/10 2*L^2/15 0 -L/10 -L^2/30 ;
       0 0    0     0 0   0 ;
       0 -6/5 -L/10 0 6/5 -L/10 ;
       0 L/10 -L^2/30 0 -L/10 2*L^2/15 ];
  %CV9 transformace do globalnich souradnic:
  Kgg = T' * Kgl * T ;
  
  % CV9 Samotna lokalizace:
  for j=1:(puzlu*ndof)
    for k=1:(puzlu*ndof)
      Kg(kcis(i,j),kcis(i,k)) = Kg(kcis(i,j),kcis(i,k))+Kgg(j,k);
    end
  end
end % to je konec smycky pres prvky

%CV9 okrajove podminky Kg:
for i=1:npodpor
	iuz=podpory(i,1);
	ismer=podpory(i,2);
	pos = (ndof*(iuz-1))+ismer;
    %u(pos) = u(pos)+podpory(i,3);
	for j=1:velikost
    Kg(pos,j) = 0.0 ;
    Kg(j,pos) = 0.0 ;
  end
	Kg(pos,pos) = 1.0 ;
end

[v,lambda] = eig(K,Kg);
D=diag(lambda);
for i=1:size(D,1)
if D(i) == 1
    D(i) = 1e999;
end
if D(i) < 0
    D(i) = 1e999;
end
end
disp('MKP:')
[Fmkp,I]=min(D);
Fmkp

% Deformovany tvar (pro primy nosnik):
mult=1;
di=I(1);
a=zeros(nuzlu,1);
b=zeros(nuzlu,1);
for i=1:nuzlu
pos=3*(i-1)+1;
a(i,1)=uzly(i,1)+v(pos,di);
pos=3*(i-1)+2;
b(i,1)=uzly(i,2)+v(pos,di)*mult;
end
plot(a,b);
