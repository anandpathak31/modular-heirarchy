function hier_ind = hier_ind_calc(Q_mat,N,layers)
m=zeros(N);

for ii = 1:(length(layers)-1)
  k1 = layers(ii).k; k2 = layers(ii+1).k;
  m(k1,k2)=1; m(k2,k1)=1;
end

hier_ind = sum(sum(Q_mat.*m));
end
