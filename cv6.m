% Reseni prihradove konstrukce (predmet NLMECH)
%
% Teorie 2. radu
%

clear; clc

% Dimenze ulohy (2 stupne volnosti:u,v, 2 uzly na prutu):
ndof = 3 ;
puzlu= 2 ;  

% ZADANI - material:
E=210e9 ;
bx=0.001;
hx=0.2;
A=bx*hx;
I=(1/12)*(bx)*(hx^3);

% ZADANI - geometrie:
uzly=[
0.0  0 ;
0.25 0 ;
0.5 0 ;
0.75 0 ;
1.0  0 ;
1.25  0 ;
1.5  0 ;
1.75  0 ;
2.0  0 
];

uzly = 8*uzly;
pruty=[
1 2 ;
2 3 ;
3 4 ;
4 5 ;
5 6 ;
6 7 ;
7 8 ;
8 9 ;
];

nuzlu=size(uzly,1);
nprutu=size(pruty,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
9 2 0 ];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
9 1 -1000 ];
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
M = zeros(velikost);
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
u=K\F ;

% Vypocet vysledku na prutech:
for i=1:nprutu
  Mn=zeros(puzlu*ndof);
  Mb=zeros(puzlu*ndof);
  Mg=zeros(puzlu*ndof);
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
  N = Fe(1);

  % geometricka matice:
  Mn = (N/L)*[ 1  0  0 -1  0  0;
               0  1  0  0 -1  0;
               0  0  0  0  0  0
              -1  0  0  1  0  0;
               0 -1  0  0  1  0
               0  0  0  0  0  0];

  Mb = (N/(30*L))* [0 0   0     0  0  0 ;
                   0 36  3*L   0 -36 3*L ;
                   0 3*L 4*L^2 0 -3*L -(L^2) ;
                   0 0   0     0  0    0 ;
                   0 -36 -3*L  0  36 -3*L ;
                   0  3*L -(L^2) 0 -3*L 4*L^2 ];
  Me = Mn + Mb ;
  Mg = T' * Me * T;

  for j=1:(puzlu*ndof)
    for k=1:(puzlu*ndof)
      M(kcis(i,j),kcis(i,k)) = M(kcis(i,j),kcis(i,k))+Mg(j,k);
    end
  end
end

% CV5 vypocet vlastnich cisel:
% Podpory:
for i=1:npodpor
	iuz=podpory(i,1);
	ismer=podpory(i,2);
	pos = (ndof*(iuz-1))+ismer;
	for j=1:velikost
    M(pos,j) = 0.0 ;
    M(j,pos) = 0.0 ;
  end
	M(pos,pos) = 1.0 ;
end

[v,d] = eig(K, M); % Nefunguje v Octave 2.x :-(

D=diag(d);
disp('Nasobitel zatizeni je:')
dd = D(1);
di = 1;
for iy=1:size(D,1)
  if D(iy) > 1.00001
    if D(iy) == 1.0
      % nic...
    else
      dd = D(iy) ;
      di = iy ;
     break ;
    end
  end
end

for iy=di:size(D,1)
    if D(iy) < dd
        if (D(iy)) > 0
            if D(iy) == 1.0
                % nic...
            else
              dd = D(iy);
              di = iy ;
            end
        end
    end
end
(dd)*1

% Deformovany tvar (pro primy nosnik):
a=zeros(nuzlu,1);
b=zeros(nuzlu,1);
for i=1:nuzlu
a(i,1)=i;
pos=3*(i-1)+2;
b(i,1)=v(pos,di);
end
plot(a,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Euler by spocital [kN]:')
L = uzly(nuzlu,1) -  uzly(1,1);
Frc = (((pi^2)*E*I)/(L^2))/1000
disp('Celkova delka byla:')
L

