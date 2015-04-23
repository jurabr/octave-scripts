% Deska / Slab 
clear; clc

% Dimenze ulohy (3 stupne volnosti:w,fix,fiy; 4 uzly na prvku):
ndof=3   ;
puzlu=4 ;  

% ZADANI - material:
E=20e9 ; % modul pruznosti
nu=0.2 ; % Poisson
h=0.2  ; % tloustka

% ZADANI - geometrie:
uzly=[
-1 -1 ;
1 -1 ;
1 1 ;
-1 1 ;
];
prvky=[
1 2 3 4 ;
];

nuzlu=size(uzly,1);
nprvku=size(prvky,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
1 3 0  ;
2 1 0 ;
2 2 0 ;
2 3 0 ;
];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
3 1 -10000 ;
4 1 -10000 ;
];
nsil=size(sily,1);

% VYPOCET: --------------------------------

% Vypocet kodovych cisel:
kcis=zeros(nprvku,ndof*puzlu);
for i=1:nprvku
  for j=1:puzlu
    for k=1:ndof
      kcis(i, ((j-1)*ndof+k))=(prvky(i,j)-1)*ndof+k ;
    end
  end
end

% Vypis kodovych cisel
kcis;

% Nulovani matic a vektoru:
velikost = nuzlu*ndof;

K = zeros(velikost);
u = zeros(velikost,1);
F = zeros(velikost,1);

% Sestaveni a lokalizace matic tuhosti:
for i=1:nprvku

  % Souradnice uzlu:
  xi=zeros(4,1);
  yi=zeros(4,1);
  for j=1:4
    xi(j)=uzly(prvky(i,j),1); 
    yi(j)=uzly(prvky(i,j),2); 
  end

  % Plocha:
  A=(xi(2)-xi(1))*(yi(3)-yi(2));

  % Matice funkci souradnic:
  S=[
1 xi(1) yi(1) xi(1)^2 xi(1)*yi(1) yi(1)^2 xi(1)^3 xi(1)^2*yi(2) xi(1)*yi(1)^2 yi(1)^3 xi(1)^3*yi(1) xi(1)*yi(1)^3;
0 0 1 0 xi(1) 2*yi(1) 0 xi(1)^2 2*xi(1)*yi(1) 3*yi(1)^2 xi(1)^3 3*xi(1)*yi(1)^2;
0 -1 0 -2*xi(1) -yi(1) 0 -3*xi(1)^2 -2*xi(1)*yi(1)  -yi(1)^2 0 -3*xi(1)^2*yi(1) -yi(1)^3;
1 xi(2) yi(2) xi(2)^2 xi(2)*yi(2) yi(2)^2 xi(2)^3 xi(2)^2*yi(2) xi(2)*yi(2)^2 yi(2)^3 xi(2)^3*yi(2) xi(2)*yi(2)^3;
0 0 1 0 xi(2) 2*yi(2) 0 xi(2)^2 2*xi(2)*yi(2) 3*yi(2)^2 xi(2)^3 3*xi(2)*yi(2)^2;
0 -1 0 -2*xi(2) -yi(2) 0 -3*xi(2)^2 -2*xi(2)*yi(2)  -yi(2)^2 0 -3*xi(2)^2*yi(2) -yi(2)^3;
1 xi(3) yi(3) xi(3)^2 xi(3)*yi(3) yi(3)^2 xi(3)^3 xi(3)^2*yi(3) xi(3)*yi(3)^2 yi(3)^3 xi(3)^3*yi(3) xi(3)*yi(3)^3;
0 0 1 0 xi(3) 2*yi(3) 0 xi(3)^2 2*xi(3)*yi(3) 3*yi(3)^2 xi(3)^3 3*xi(3)*yi(3)^2;
0 -1 0 -2*xi(3) -yi(3) 0 -3*xi(3)^2 -2*xi(3)*yi(3)  -yi(3)^2 0 -3*xi(3)^2*yi(3) -yi(3)^3;
1 xi(4) yi(4) xi(4)^2 xi(4)*yi(4) yi(4)^2 xi(4)^3 xi(4)^2*yi(4) xi(4)*yi(4)^2 yi(4)^3 xi(4)^3*yi(4) xi(4)*yi(4)^3;
0 0 1 0 xi(4) 2*yi(4) 0 xi(4)^2 2*xi(4)*yi(4) 3*yi(4)^2 xi(4)^3 3*xi(4)*yi(4)^2;
0 -1 0 -2*xi(4) -yi(4) 0 -3*xi(4)^2 -2*xi(4)*yi(4)  -yi(4)^2 0 -3*xi(4)^2*yi(4) -yi(4)^3;
];
  Si=inv(S)

  D=((E*h^3)/(12*(1-nu^2)))*[
  1 nu 0 ;
  nu 1 0 ;
  0 0 0.5*(1-nu) ];
  
  % Numericka integrace pro prvek velikosti 1x1 !!!
  Ke=zeros(12,12);
  xj=[-0.5774 ; 0.5775 ; 0.5774; -0.5774] ;
  yj=[-0.5774 ; -0.5775 ; 0.5774; 0.5774] ;

  % Dimenze prvku:
  a = xi(2)-xi(1) ;
  b = yi(3)-yi(2) ; 

  for j=1:4
      B=-[0 0 0 2 0 0 6*xj(j) 2*yj(j) 0 0 6*xj(j)*yj(j) 0 ;
         0 0 0 0 0 2 0 0 2*xj(j) 6*yj(j) 0 6*xj(j)*yj(j) ;
         0 0 0 0 2 0 0 4*xj(j) 4*yj(j) 0 6*xj(j)^2 6*yj(j)^2 ;
        ];

      Ke = Ke +  Si'*B'*D*B*Si*1;
  end;
  Ke = A*Ke
  
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
u=K\F

gmult = 0.025*max(uzly(:))/max(u); % kvuli grafice

% Vypocet vysledku na prvcich:
for i=1:nprvku
  ue  = zeros(puzlu*ndof,1);

  % Ziskani lokalnich vektoru deformaci:
  for j=1:(puzlu*ndof)
    ue(j) = u(kcis(i,j));
  end
end
