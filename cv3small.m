% Reseni prihradove konstrukce (predmet NLMECH)
% Nyni umi: 
% * ruzne plochy prutu
% * vykresleni deformaci

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
1 0 ;
3 0 ;
4 0 ;
1 1 ;
3 1 ];
pruty=[
1 2 0.1 0;
2 3 0.1 0;
3 4 0.1 0;
1 5 0.1 0;
2 5 0.1 0;
3 6 0.1 0;
4 6 0.1 0;
5 6 0.1 0;
2 6 0.01 0;
3 5 0.01 0];
nuzlu=size(uzly,1);
nprutu=size(pruty,1);

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 0 ;
1 2 0  ;
4 1 0 ;
4 2 0 ];
npodpor=size(podpory,1);

% ZADANI - sily (uzel,smer,velikost):
sily=[
5 2 -10000 ];
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

for ii=1:10 % CV2: Hlavni cyklus pro iteraci
    pocet_zmen = 0; % CV2: pocet zmen stavu prutu (tazeny/tlaceny)
    disp('ITERACE') % CV2: vypis iterace
    ii
    
K = zeros(velikost);
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
  
  % CV2: test, zda prut byl tlaceny
  if pruty(i,3) < 0.05
    if pruty(i,4) < 1
        D=[E]; % Tazeny
    else
        D=[E/1000]; % Tlaceny
        %disp('Pocitam s tlacenym prutem...')
        %i
    end
  else
      D=[E];
  end
  
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
  F(pos) = F(pos)+sily(i,3);
end

if ii > 1
  F = F-g ;
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

% CV3 vektor pro nevyvazene sily:
g = zeros(velikost,1);

% Vypocet vysledku na prutech:
for i=1:nprutu
	ue  = zeros(puzlu*ndof,1);
	uel = zeros(puzlu*ndof,1);

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
  
  % CV2: test, zda prut byl tlaceny
  if pruty(i,3) < 0.05
    if pruty(i,4) < 1
        D=[E]; % Tazeny
    else
        D=[E/1000]; % Tlaceny
        %disp('Pocitam s tlacenym prutem...')
    end
    else
      D=[E];
  end
  
  Si=inv(S);
  A = pruty(i,3); % Toto je nove pro NLM
  Kel=A*L*Si'*B'*D*B*Si;
  Kelb=[Kel(1,1) 0 Kel(1,2) 0 ;
      0 0 0 0 ;
      Kel(2,1) 0 Kel(2,2) 0
      0 0 0 0 ];
	
  % sily v prutech:
	Fe = Kelb * uel;
    
    % CV2: podminka, zda je prut tazeny:
    if pruty(i,3) < 0.05
        if (Fe(1) > 0.0)
            if pruty(i,4) == 0
                pocet_zmen = pocet_zmen + 1 % TODO ma byt 1 - vypnuto kvuli testovani
                disp ('Nalezen tlaceny prut:')
                i
            end
            pruty(i,4) = 1 ; % TODO ma byt 1 - vypnuto kvuli testovani
        end
        
        % CV2: pripad, ze prut se znovu stane tazenym
        if (Fe(1) < 0.0)
            if pruty(i,4) == 1
                pocet_zmen = pocet_zmen + 1
                pruty(i,4) = 0;
                disp('Zmena stavu na tazeny!')
                i
            end
        end
    end    
        

    % CV3 Vektor nevyvazebych sil:
    Keg = T' * Kelb * T;
    Feg = Keg * ue ;
    for j=1:(puzlu*ndof)
      g(kcis(i,j)) = g(kcis(i,j)) +  Feg(j) ;
    end

    % Kvuli kresleni:
    souradnice(:,:,i)=[
        uzly(pruty(i,1),1) uzly(pruty(i,1),2) ;
        uzly(pruty(i,2),1) uzly(pruty(i,2),2) ];
    
    souradnicedef(:,:,i)=[
        uzly(pruty(i,1),1)+ue(1)*gmult uzly(pruty(i,1),2)+ue(2)*gmult ;
        uzly(pruty(i,2),1)+ue(3)*gmult uzly(pruty(i,2),2)+ue(4)*gmult ];
end

  % CV3 vycisteni g od podporovych reakce:
  for k=1:npodpor
  	iuz=podpory(k,1);
   	ismer=podpory(k,2);
   	pos = (ndof*(iuz-1))+ismer;
   	g(pos) = 0.0 ;
  end
  % CV3 zapocitani zatizeni do g
  for k=1:nsil
	  iuz=sily(k,1);
	  ismer=sily(k,2);
	  pos = (ndof*(iuz-1))+ismer;
    g(pos) = g(pos)-sily(k,3);
  end


  pom_norm = norm(g) / norm(F); 
  if pom_norm > 0.0001
    pocet_zmen = pocet_zmen + 1 
    disp('Nevychazi rovnovaha')
  end

% CV2: ukonceni iterace pokud uz se pruty nemeni:
    if pocet_zmen == 0
        disp('Koncim iteraci: 0 zmen')
        break;
    end
end % CV2: konec hlavniho cyklu pro iteraci (ii)

% % Grafika - vykresleni konstrukce:
% hold on
% for i=1:nprutu
%     if pruty(i,3) > 0.05
%         % tvar konstrukce:
%         plot(souradnice(:,1,i),souradnice(:,2,i),'-r');
%         % vykresleni deformaci:
%     plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
%     else
%         if pruty(i,4) < 1
%         % tvar konstrukce:
%         plot(souradnice(:,1,i),souradnice(:,2,i),'-b');
%         % vykresleni deformaci:
%         plot(souradnicedef(:,1,i),souradnicedef(:,2,i),'-g');
%         else
%             plot(souradnice(:,1,i),souradnice(:,2,i),'-y');
%         end
%     end
% end
% axis equal
% hold off
% 
%pruty
