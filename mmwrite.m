function mmwrite(file, a)
% file: filename given like 'stuff.mat'
% a:    matrix a(n,n)

[ilen, jlen] = size(a);
c = 0;

for i=1:ilen
  for j=1:jlen
    if (abs(a(i,j)) > 1e-6)
      c = c + 1 ;
    end
  end
end

b = zeros(c,2);
c1 = 0;

for i=1:ilen
  for j=1:jlen
    if (abs(a(i,j)) > 1e-6)
      c1 = c1 + 1 ;
      b(c1,1) = i ;
      b(c1,2) = j ;
    end
  end
end

file_id = fopen(file,'w');

fprintf(file_id,'%%%%MatrixMarket matrix coordinate real general\n');;
fprintf(file_id,'%% A %ix%i sparse matrix\n', ilen, jlen);
fprintf(file_id,'%i %i %i\n', ilen, jlen, c );

for i=1:c
  fprintf(file_id,'%i %i %e\n', b(i,1), b(i,2), a(b(i,1),b(i,2)));
end

fclose(file_id);
endfunction
