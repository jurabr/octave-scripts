% Reseni prihradove konstrukce (predmet NLMECH)
%
% Stabilita prutu (teorie 2. radu)
%
clear; clc

% Dimenze ulohy (2 stupne volnosti:u,v, 2 uzly na prutu):
ndof=2   ;
puzlu=2 ;  

% ZADANI - material:
E=210e9 ;
A=0.1  ;

% ZADANI - geometrie:
uzly=[
0 0 ;
2 0 ;
5 0 ;
8 0 ;
11 0 ;
13 0
2 1.6;
5 1.6;
8 1.6;
11 1.6
];
pruty=[
1 2 0.1 0;
2 3 0.1 0;
3 4 0.1 0;
4 5 0.1 0;
5 6 0.1 0;
1 7 0.1 0;
2 7 0.1 0;
3 8 0.1 0;
4 9 0.1 0;
5 10 0.1 0;
6 10 0.1 0;
7 8 0.1 0;
8 9 0.1 0;
9 10 0.1 0;

2 8 0.01 0;
3 7 0.01 0;
3 9 0.01 0;
4 8 0.01 0;
4 10 0.01 0;
5 9 0.01 0];

nuzlu=size(uzly,1);
nprutu=size(pruty,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
6 1 0 ;
6 2 0 ];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
7 2 -10000 ;
9 2 -10000 ;
10 1 5000 ];
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
  % Delka prutu
  dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
  dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
  L = sqrt(dx2 + dy2);

  % Matice tuhosti:
  S=[1 0 ; 1 L];
  B=[0 1];
  
  Si=inv(S);
  A = pruty(i,3); % Toto je nove pro NLM
  Kel=A*L*Si'*B'*E*B*Si;
  Keg=[Kel(1,1) 0 Kel(1,2) 0 ;
      0 0 0 0 ;
      Kel(2,1) 0 Kel(2,2) 0
      0 0 0 0 ];
  
  % Transformace:
  s = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2))/L;
  c = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1))/L;
  T = [ c s 0 0 ; -s c 0 0 ; 
        0 0 c s ; 0 0 -s c ] ;
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

gmult = 0.025*max(uzly(:))/max(u); % kvuli grafice

% Vypocet vysledku na prutech:
for i=1:nprutu
	ue  = zeros(puzlu*ndof,1);

	% Ziskani lokalnich vektoru deformaci:
	for j=1:(puzlu*ndof)
		ue(j) = u(kcis(i,j));
	end

  % CV3 zapocitani deformaci konstrukce do sourdnic uzlu:
  uzly(pruty(i,1),1) = uzly(pruty(i,1),1) +  ue(1) ;
  uzly(pruty(i,1),2) = uzly(pruty(i,1),2) +  ue(2) ;
  uzly(pruty(i,2),1) = uzly(pruty(i,2),1) +  ue(3) ;
  uzly(pruty(i,2),2) = uzly(pruty(i,2),2) +  ue(4) ;

  % Transformace:
  dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
  dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
  L = sqrt(dx2 + dy2);

  s = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2))/L;
  c = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1))/L;
  T = [ c s 0 0 ; -s c 0 0 ; 
        0 0 c s ; 0 0 -s c ] ;
  uel = T * ue;

  % Matice tuhosti (jen lokalni):
  S=[1 0 ; 1 L];
  B=[0 1];
  Si=inv(S);
  A = pruty(i,3); % Toto je nove pro NLM
  Kel=A*L*Si'*B'*E*B*Si;
  Kelb=[Kel(1,1) 0 Kel(1,2) 0 ;
      0 0 0 0 ;
      Kel(2,1) 0 Kel(2,2) 0
      0 0 0 0 ];
	
  % sily v prutech:
  Fe = Kelb * uel;
  N = Fe(1);

  % CV5 geometricka matice dle prednasky:
  Me = (N/L)*[ 0  0  0  0 ;
               0  1  0 -1 ;
               0  0  0  0 ;
               0 -1  0  1 ];
  Mg = T' * Me * T;
  
  for j=1:(puzlu*ndof)
    for k=1:(puzlu*ndof)
      M(kcis(i,j),kcis(i,k)) = M(kcis(i,j),kcis(i,k))+Mg(j,k);
    end
  end
  
  % Kvuli kresleni:
    souradnice(:,:,i)=[
        uzly(pruty(i,1),1) uzly(pruty(i,1),2) ;
        uzly(pruty(i,2),1) uzly(pruty(i,2),2) ];
end

% CV5 vypocet vlastnich cisel (mocnina nasobitele zatizeni):
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
di = 0;
for iy=1:size(D,1)
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
sqrt(dd)

% CV5 deformovana konstrukce
gmult = 0.05*max(uzly(:))/max(v(:)); % kvuli grafice

for i=1:nprutu
	ue  = zeros(puzlu*ndof,1);
    
    ue(1) = v(2*pruty(i,1)-1,di);
    ue(2) = v(2*pruty(i,1),di);
    ue(3) = v(2*pruty(i,2)-1,di);
    ue(4) = v(2*pruty(i,2),di);
    ue;
    souradnicedef(:,:,i)=[
        uzly(pruty(i,1),1)+ue(1)*gmult uzly(pruty(i,1),2)+ue(2)*gmult ;
        uzly(pruty(i,2),1)+ue(3)*gmult uzly(pruty(i,2),2)+ue(4)*gmult ];
end

 % Grafika - vykresleni konstrukce:
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