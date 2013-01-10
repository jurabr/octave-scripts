% Reseni prostorove prihradove konstrukce 
% 
% je treba jeste overit transforamci Kg, ale neco to pocita
% 

clear; clc

% Dimenze ulohy (2 stupne volnosti:u,v, 2 uzly na prutu):
ndof=3   ;
puzlu=2 ;  

% ZADANI - material:
E=20e9 ;
A=0.1  ;

% ZADANI - geometrie:
uzly=[
0 0 0;
4 0 0;
4 4 0];
uzly1=uzly;
pruty=[
1 2 0.1 0;
2 3 0.1 0;
1 3 0.1 0
];

nuzlu=size(uzly,1);
nprutu=size(pruty,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 2 0 ;
1 3 0 ;
2 1 0 ;
2 2 0 ;
2 3 0 ;
3 1 0 ;
3 3 0 
];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
3 2 -10000 ];

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
Kg= zeros(velikost);
u = zeros(velikost,1);
F = zeros(velikost,1);

% Sestaveni a lokalizace matic tuhosti:
for i=1:nprutu
  % Delka prutu
  dx = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1));
  dy = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2));
  dz = (uzly(pruty(i,2),3)-uzly(pruty(i,1),3));
  L = sqrt(dx^2 + dy^2 + dz^2);
  a = dx / L ;
  b = dy / L ;
  c = dz / L ;
  T = [ a b c 0 0 0 ;
        0 0 0 0 0 0 ;
        0 0 0 0 0 0 ;
        0 0 0 a b c ;
        0 0 0 0 0 0 ;
        0 0 0 0 0 0 ];

  % Matice tuhosti:
  A = pruty(i,3); % Toto je nove pro NLM

  Kelb = ((E*A)/L) * [ 
    1 0 0 -1 0 0 ;
    0 0 0 0 0 0 ;
    0 0 0 0 0 0 ;
    -1 0 0 1 0 0 ;
    0 0 0 0 0 0 ;
    0 0 0 0 0 0 ];

  Ke = T' * Kelb * T ;
    
  for ii=1:3
    for j=1:3
      Ke(ii+3, j) = (-1.0)*Ke(ii,j) ;
      Ke(ii+3, j+3) = Ke(ii,j) ;
      Ke(ii, j+3) = (-1.0)*Ke(ii,j) ;
    end
  end

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

gmult = 0.025*max(uzly(:))/max(abs(u)); % kvuli grafice

% Vypocet vysledku na prutech:
for i=1:nprutu
  ue  = zeros(puzlu*ndof,1);
  uel = zeros(puzlu*ndof,1);

  % Ziskani lokalnich vektoru deformaci:
  for j=1:(puzlu*ndof)
    ue(j) = u(kcis(i,j));
  end

  % Transformace:
  dx = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1));
  dy = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2));
  dz = (uzly(pruty(i,2),3)-uzly(pruty(i,1),3));
  L = sqrt(dx^2 + dy^2 + dz^2);

  % alternativni transformacni matice:
  x = zeros(3,1);
  y = zeros(3,1);
  z = zeros(3,1);

  x(1) = uzly(pruty(i,1),1);
  y(1) = uzly(pruty(i,1),2);
  z(1) = uzly(pruty(i,1),3);
  x(2) = uzly(pruty(i,2),1);
  y(2) = uzly(pruty(i,2),2);
  z(2) = uzly(pruty(i,2),3);
  
  if (y(1)-y(2))*(x(2)-x(1)) == 0
    x(3) = x(2) + 1  ;
    y(3) = y(2)   ;
    z(3) = z(2) + 1  ;
  else
    x(3) = x(2) ;
    y(3) = y(2);
    z(3) = z(2) + 1 ;
  end

  a1 = (x(2)-x(1)) / L ;
  b1 = (y(2)-y(1)) / L ;
  c1 = (z(2)-z(1)) / L ;
  
  A = (y(3)-y(1))*(z(2)-z(1)) - (y(2)-y(1))*(z(3)-z(1)) 
  B = (z(3)-z(1))*(x(2)-x(1)) - (z(2)-z(1))*(x(3)-x(1)) 
  C = (x(3)-x(1))*(y(2)-y(1)) - (x(2)-x(1))*(y(3)-y(1)) 
  
  d = sqrt(A^2 + B^2 + C^2 ) ;
  a2 = A/d ;
  b2 = B/d ;
  c2 = C/d ;

  a3 = b1*c2 - c1*b2 ;
  b3 = a2*c1 - a1*c2 ;
  c3 = a1*b2 - a2*b1 ;
 
  tt = [ a1 b1 c1 ;
         a2 b2 c3 ;
	 a3 b3 c3 ]

  TT = zeros(6,6) ;
  for ii=1:2
    for jj=1:3
      for kk=1:3
      TT( (ii-1)*3+jj,(ii-1)*3+kk ) = tt(jj,kk);
      end
    end
  end

  % Matice tuhosti:
  A = pruty(i,3); % Toto je nove pro NLM
  uel = T * ue ;

  % Matice tuhosti (jen lokalni):
  Kelb = ((E*A)/L) * [ 
    1 0 0 -1 0 0 ;
    0 0 0 0 0 0 ;
    0 0 0 0 0 0 ;
    -1 0 0 1 0 0 ;
    0 0 0 0 0 0 ;
    0 0 0 0 0 0 ];

  % Sila v prutu:
  Fe = Kelb * uel ;
  N = Fe(1)

  % Geometricka matice:
  Kge = (N/L) * [ 0 0 0 0 0 0 ;
                  0 1 0 0 -1 0 ;
	          0 0 1 0 0 -1 ;
	          0 0 0 0 0 0 ;
	          0 -1 0 0 1 0 ;
	          0 0 -1 0 0 1 ];

  Kgg = TT' * Kge * TT ;

  % Lokalizace Kgg do Kg:
  for j=1:(puzlu*ndof)
    for k=1:(puzlu*ndof)
      Kg(kcis(i,j),kcis(i,k)) = Kg(kcis(i,j),kcis(i,k))+Kgg(j,k);
    end
  end

  % Kvuli kresleni:
  souradnice(:,:,i)=[
        uzly(pruty(i,1),1) uzly(pruty(i,1),2) ;
        uzly(pruty(i,2),1) uzly(pruty(i,2),2) ];
    
  souradnicedef(:,:,i)=[
        uzly(pruty(i,1),1)+ue(1)*gmult uzly(pruty(i,1),2)+ue(2)*gmult ;
        uzly(pruty(i,2),1)+ue(3)*gmult uzly(pruty(i,2),2)+ue(4)*gmult ];
end

Kg

% Vypocet vlastnich cisel:
% Podpory:
for i=1:npodpor
	iuz=podpory(i,1);
	ismer=podpory(i,2);
	pos = (ndof*(iuz-1))+ismer;
	for j=1:velikost
    Kg(pos,j) = 0.0 ;
    Kg(j,pos) = 0.0 ;
  end
    Kg(pos,pos) = 1.0 ;
end



[v,d] = eig(K, Kg); % Nefunguje v Octave 2.x :-(

D=diag(d)

break

% % Grafika - vykresleni konstrukce:
 hold on
 for i=1:nprutu
     if pruty(i,3) > 0.05
         % tvar konstrukce:
         plot(souradnice(:,1,i),souradnice(:,2,i),'-r');
         % vykresleni deformaci:
     plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
     else
         if pruty(i,4) < 1
         % tvar konstrukce:
         plot(souradnice(:,1,i),souradnice(:,2,i),'-b');
         % vykresleni deformaci:
         plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
         else
             plot(souradnice(:,1,i),souradnice(:,2,i),'-y');
         end
     end
 end
 axis equal
 hold off
 
%pruty
