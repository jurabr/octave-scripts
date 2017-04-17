% Dimenze ulohu (1 stupen volnosti: T, 3 uzly na prvku):
ndof=1   
puzlu=3   

% ZADANI - material:
lambda=5.0
t=1

% ZADANI - geometrie:
uzly=[
0 0 ;
1 0 ;
2 0 ;
0 1 ;
1 1 ;
2 1 ;
0 2 ;
1 2 ;
2 2 ];
prvky=[
1 5 4 ;
1 2 5 ;
2 6 5 ;
2 3 6 ;
4 8 7 ;
4 5 8 ;
5 9 8 ;
5 6 9 ];
nuzlu=size(uzly,1)
nprvku=size(prvky,1)

% ZADANI - podpory (uzel,smer,velikost):
podpory=[
1 1 -5 ;
2 1 -5 ;
3 1 -5 ;
7 1 10 ;
8 1 10 ;
9 1 10 ];
npodpor=size(podpory,1)

% ZADANI - sily (uzel,smer,velikost):
sily=[
3 1 0 ;
4 1 0 ];
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
  % Matice tuhosti:
  x= zeros(3,1);
  y= zeros(3,1);
  
  for j=1:puzlu
      x(j) = uzly(prvky(i,j), 1);
      y(j) = uzly(prvky(i,j), 2);
  end
  
  S=[ x(1) y(1) 1 ;          
      x(2) y(2) 1 ;
      x(3) y(3) 1 ];
  
  B=[1 0 0 ;
     0 1 0 ];
  
  D= -lambda*[
  1 0 ;
  0 1 ];
  
  Si=inv(S);
  
  % Plocha prvku
  A = 0.5*(x(1)*y(2)-x(2)*y(1)+x(2)*y(3)-x(3)*y(2)+x(3)*y(1)-x(1)*y(3));
  
   % kvuli kresleni:
  souradnice(:,:,i)=[x(1) y(1);x(2) y(2); x(3) y(3); x(1) y(1)];
    
  Ke=A*t*Si'*B'*D*B*Si;
     
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
        F(j) = F(j)-podpory(i,3)*K(j,pos);
        K(pos,j) = 0.0 ;
        K(j,pos) = 0.0 ;
    end
	K(pos,pos) = 1.0 ;
    if abs(podpory(i,3)) > 0 
        F(pos) = 1.0*podpory(i,3);
    end
end


% Kresleni:
hold on
for i=1:nprvku
  %plot(souradnice(:,1,i),souradnice(:,2,i),'--r')
end



% RESENI soustavy rovnic (vysledne teploty):
u=K\F

% Malovani barvicek podle teplot:
axis equal
maplen = 12 ; % pocet barev
colormap (jet (maplen)); % barevna mapa
for i=1:nprvku
  c(i) = 0 ;
  for j=1:puzlu % souradnice ve tvaru, co potrebuje "patch"
    xi(i,j) = uzly(prvky(i,j),1) ;
    yi(i,j) = uzly(prvky(i,j),2) ;
    c(i) = c(i) + u(kcis(i,j)) ;
  end
	c(i) = c(i)/puzlu; % prumerna hodnota na prvku
end
cmi = min(c) ; cma=max(c); % max. a min. hodnota
for i=1:nprvku
	c(i) = (((c(i)-cmi) / ((cma-cmi)))*maplen) ;
end
caxis([cmi cma]); % rozsah osy
colorbar; % teplomer
p = patch (xi', yi'); % kresleni ploch
set (p, "cdatamapping", "direct", "facecolor", "flat", "cdata", c ); %vybarveni
