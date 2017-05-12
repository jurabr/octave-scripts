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
7 1 15 ;
8 1 15 ;
9 1 15 ;
];
npodpor=size(podpory,1)

% ZADANI - sily (uzel,smer,velikost):
sily=[
3 1 0 ;
4 1 0 ];
nsil=size(sily,1);

% KRESLENI:
s_x = zeros(puzlu*nprvku,1);
s_y = zeros(puzlu*nprvku,1);

for i=0:(nprvku-1)
  for j=1:puzlu
    s_x(puzlu*i+j) = uzly(prvky(i+1,j),1);
    s_y(puzlu*i+j) = uzly(prvky(i+1,j),2);
  end
end

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
  %souradnice(:,:,i)=[x(1) y(1);x(2) y(2); x(3) y(3); x(1) y(1)];
    
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

% RESENI soustavy rovnic (vysledne teploty):
u=K\F


% Kresleni: ------------------------------------------
hold on
plot(s_x,s_y,'0');

% Pokus o izolinie pro teploty:
t_max = max(u);
t_min = min(u);
t_del = (t_max - t_min) / 4;
t_lim = zeros(5,2);
t_lim(1,2) = 3; % colors
t_lim(2,2) = 5;
t_lim(3,2) = 2;
t_lim(4,2) = 4;
t_lim(5,2) = 1;
t_lim(1,1) = t_min ;
for i=2:5
  t_lim(i,1) = t_min+ (i-1)*t_del ;
end
t_lim

% Isolines 
for i=1:nprvku
  xyt = zeros(puzlu+1,3); % element data
  px  = zeros(2,1);       % line coordinates
  py  = zeros(2,1);       % line coordinates
  for j=1:puzlu
    xyt(j,1) = uzly(prvky(i,j),1) ;
    xyt(j,2) = uzly(prvky(i,j),2) ;
    xyt(j,3) = u(kcis(i,j)) ;
  end
  xyt(puzlu+1,1) = uzly(prvky(i,1),1) ;
  xyt(puzlu+1,2) = uzly(prvky(i,1),2) ;
  xyt(puzlu+1,3) = u(kcis(i,1)) ;

  for k=1:5 % test - isoline
    pos = 0 ;
    for j=1:puzlu
      if (xyt(j,3)<xyt(j+1,3))
      if ( t_lim(k,1) > xyt(j,3) ) 
        if ( t_lim(k,1) < xyt(j+1,3) ) 
          pos = pos+1 ;
          L = xyt(j+1,3)-xyt(j,3);
          if (L == 0)
            mult = 0.5 ;
          else
            mult = 1-abs(t_lim(k,1)-xyt(j,3)) / L ;
          end
          px(pos) = mult * ( xyt(j+1,1)-xyt(j,1)) + xyt(j,1);
          py(pos) = mult * ( xyt(j+1,2)-xyt(j,2)) + xyt(j,2);
          if (pos == 2)
            plot(px,py,num2str(t_lim(k,2)));
            pos = 0 ;
          end
        end
      end
      else
      %----------------------------------------------
        if ( t_lim(k,1) < xyt(j,3) ) 
        if ( t_lim(k,1) > xyt(j+1,3) ) 
          pos = pos+1 ;
          L = abs(xyt(j,3)-xyt(j+1,3));
          if (L < 0.00001)
            mult = 0.5 ;
          else
            mult = abs((t_lim(k,1)-xyt(j+1,3)) / L) ;
          end
          px(pos) = mult * ( xyt(j+1,1)-xyt(j,1)) + xyt(j,1);
          py(pos) = mult * ( xyt(j+1,2)-xyt(j,2)) + xyt(j,2);
          if (pos == 2)
            plot(px,py,num2str(t_lim(k,2)));
            pos = 0 ;
          end
        end
      end
      %----------------------------------------------
      end
    end
  end

end
hold off
