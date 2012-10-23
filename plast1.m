% Reseni prihradove konstrukce (predmet NLMECH)
% Nyni umi: 
% * ruzne plochy prutu
% * vykresleni deformaci

% Dimenze ulohy (2 stupne volnosti:u,v, 2 uzly na prutu):
ndof=2   
puzlu=2   

% ZADANI - material:
E=210e9
A=0.1

% ZADANI - geometrie:
uzly=[
0 0 ;
1 0 ;
3 0 ;
4 0 ;
1 1 ;
3 1 ]
pruty=[
1 2 0.01 ;
2 3 0.01 ;
3 4 0.01 ;
1 5 0.01 ;
2 5 0.01 ;
3 6 0.01 ;
4 6 0.01 ;
5 6 0.01 ;
2 6 0.01 ;
3 5 0.01 ]
nuzlu=size(uzly,1)
nprutu=size(pruty,1)

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
4 1 0 ;
4 2 0 ]
npodpor=size(podpory,1)

% ZADANI - sily (uzel,smer,velikost):
sily=[
3 2 -1000000 ;
5 2 -100000 ]
nsil=size(sily,1)

% CV5: ZADANI - mez kluzu:
fy = 235e6 ;
E1 = 210e6 ;

% CV5: pole s moduly pruznosti:
Ecka = zeros(nprutu,1);
for i=1:nprutu
    Ecka(i) = E;
end

NN = zeros(nprutu,1);

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


prirustek = 1/1000 ; %CV5 1/pocet kroku iii
ropucha = 0 ;

% CV5: vypocet - smycka pro zatizeni:
for iii=1:1000

% Nulovani matic a vektoru:
velikost = nuzlu*ndof;
K = zeros(velikost);
u = zeros(velikost,1);
F = zeros(velikost,1);

% CV5: celkove deformace
u_sum = zeros(velikost,1); 

% Sestaveni a lokalizace matic tuhosti:
for i=1:nprutu
	% Delka prutu
	dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
	dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
	L = sqrt(dx2 + dy2);

	% Matice tuhosti:
  S=[1 0 ; 1 L];
  B=[0 1];
  D=[Ecka(i)]; % CV5 uprava na promenne E
  Si=inv(S);
  A = pruty(i,3); % Toto je nove pro NLM
  Kel=A*L*Si'*B'*D*B*Si;
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
    F(pos) = (prirustek*F(pos)) + sily(i,3); %CV5 uprava na prirustek
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

u_sum = u_sum + u ; % CV5: celkova deformace

gmult = 0.025*max(uzly(:))/max(u); % kvuli grafice

% Vypocet vysledku na prutech:
for i=1:nprutu
	ue  = zeros(puzlu*ndof,1);
	uel = zeros(puzlu*ndof,1);
    
    uet  = zeros(puzlu*ndof,1); % CV5
    uelt = zeros(puzlu*ndof,1); % CV5

    
	% Ziskani lokalnich vektoru deformaci:
	for j=1:(puzlu*ndof)
		ue(j) = u(kcis(i,j));
        uet(j) = u_sum(kcis(i,j)); % CV5
	end

	% Transformace:
	dx2 = (uzly(pruty(i,1),1)-uzly(pruty(i,2),1))^2;
	dy2 = (uzly(pruty(i,1),2)-uzly(pruty(i,2),2))^2;
	L = sqrt(dx2 + dy2);
    
    s = (uzly(pruty(i,2),2)-uzly(pruty(i,1),2))/L;
	c = (uzly(pruty(i,2),1)-uzly(pruty(i,1),1))/L;
  T = [ c s 0 0 ; -s c 0 0 ; 
        0 0 c s ; 0 0 -s c ] ;
  uel = T * ue;
  uelt = T * uet; % CV5

  % Matice tuhosti (jen lokalni):
  S=[1 0 ; 1 L];
  B=[0 1];
  D=[E];
  Si=inv(S);
  A = pruty(i,3); % Toto je nove pro NLM
  Kel=A*L*Si'*B'*D*B*Si;
  Keg=[Kel(1,1) 0 Kel(1,2) 0 ;
      0 0 0 0 ;
      Kel(2,1) 0 Kel(2,2) 0
      0 0 0 0 ];
	
    % VYSLEDEK - sily v prutech:
	Fe = Keg * uel;
    
    % CV5: kontrola stavu: pruzny/plasticky:
    sigma = ( NN(i) + Fe(1) ) / A ;
    if abs(sigma) > fy
        disp('Plasticky stav')
        Ecka(i) = E1/(1-(E1/E)) ;
            NN(i) = NN(i) + (E1/E)*Fe(1);
    else
        Ecka(i) = E ;
        NN(i) = NN(i) + Fe(1);
    end
    
    souradnice(:,:,i)=[
        uzly(pruty(i,1),1) uzly(pruty(i,1),2) ;
        uzly(pruty(i,2),1) uzly(pruty(i,2),2) ];
    
    souradnicedef(:,:,i)=[
        uzly(pruty(i,1),1)+ue(1)*gmult uzly(pruty(i,1),2)+ue(2)*gmult ;
        uzly(pruty(i,2),1)+ue(3)*gmult uzly(pruty(i,2),2)+ue(4)*gmult ];

  % CV5: podminka poruseni (epsilon <0.15)
  dxl2 = (uelt(1)-uelt(3))^2;
  dyl2 = (uelt(2)-uelt(4))^2;
  dL =  sqrt(dxl2 + dyl2);
  epsilon = dL / L ;
  if (epsilon > 0.15)
      ropucha = ropucha + 1 ;
  end
end

% CV5: vypise pocet prutu s prekrozenym epsilon
if ropucha > 0
    disp('Porusene pruty')
    ropucha
    iii
    break;
end
end %CV5 smycka pres zatizeni


% Grafika - vykresleni konstrukce:
hold on
for i=1:nprutu
    if Ecka(i,1) > E1/(1-(E1/E))
        % tvar konstrukce:
        plot(souradnice(:,1,i),souradnice(:,2,i),'-r');
        % vykresleni deformaci:
    plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
    else
        % tvar konstrukce:
        plot(souradnice(:,1,i),souradnice(:,2,i),'-b');
        % vykresleni deformaci:
        plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
    end
end
%axis equal
hold off

u_sum
NN
