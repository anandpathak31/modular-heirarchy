function [layers, heir_ind] = heir_simulated_annealing(adj)

N = length(adj(:,1));
lsize   = floor(N/5);
runtime = 2e7;

T = 10;
lambda = 1/5e5;

L = sum(sum(adj));
indeg = sum(adj,2);
outdeg = sum(adj,1);

Q_mat = (adj-indeg*outdeg/L)/L;

r = randperm(N);
for ii = 1:4
    layers(ii).k = r([lsize*(ii-1)+[1:lsize]]);
end
layers(5).k = r((4*lsize+1):N);

heir_ind = heir_ind_calc(Q_mat,N,layers);

layers_id = zeros(N,1);
update = 1; failed = 0;

for ii = 1:runtime
    k = length(layers);
    tot_nbr = k*(N+k); % total neighbours
    %N*(k-1) + N + k*(k-1)/2  +  k*(k-1)/2 + k;
    
    p(1) = N*(k-1)/tot_nbr;
    p(2) = N/tot_nbr;
    p(3) = (k*(k-1)/2)/tot_nbr;
    p(4) = (k*(k-1)/2)/tot_nbr;
    p(5) = k/tot_nbr;
    prange = cumsum(p);
    
    prob = rand;
    layers_new=layers;
    
    if prob < prange(1)
        % swap a node from one layer to any other
        
        % check if a layer has been an updated in the previous loop
        
        ll=0;
        for jj = 1:k
            pop = length(layers_new(jj).k);
            layers_id((ll+1):(ll+pop))=jj;
            ll=ll+pop;
        end
        
        rand_node = ceil(N*rand);
        old_lyr = layers_id(rand_node);
        kl = length(layers_new(old_lyr).k);
        nn = ceil(kl*rand);
        pkd_node = layers_new(old_lyr).k(nn);
        
        ar_lyr = setdiff(1:k,old_lyr);
        new_lyr_id = ceil((k-1)*rand);
        new_lyr = ar_lyr(new_lyr_id);
        
        layers_new(new_lyr).k = [layers_new(new_lyr).k pkd_node];
        layers_new(old_lyr).k = setdiff(layers_new(old_lyr).k,pkd_node);
        
        kt = length(layers_new);
        len = zeros(1,kt);
        for jj = 1:kt
            len(jj) = length(layers_new(jj).k);
        end
        
        id = find(len==0);
        
        if ~isempty(id)
            id2 = setdiff(1:kt,id);
            layers_new = layers_new(id2);
            update = 1;
        end
        
    elseif (prob>=prange(1)) && (prob<prange(2))
    % pick a node, make a new layer at the end and put it in the new layer

    % check if a layer has been an updated in the previous loop
        
        ll=0;
        for jj = 1:k
            pop = length(layers_new(jj).k);
            layers_id((ll+1):(ll+pop))=jj;
            ll=ll+pop;
        end
        
        rand_node = ceil(N*rand);
        old_lyr = layers_id(rand_node);
        kl=length(layers_new(old_lyr).k);
        nn = ceil(kl*rand);
        pkd_node = layers_new(old_lyr).k(nn);
        
        layers_new(k+1).k = pkd_node;
        layers_new(old_lyr).k = setdiff(layers_new(old_lyr).k,pkd_node);
        
        kt = length(layers_new);
        len = zeros(1,kt);
        for jj =1:kt
            len(jj)=length(layers_new(jj).k);
        end
        
        id = find(len==0);
        
        if ~isempty(id)
            id2 = setdiff(1:kt,id);
            layers_new = layers_new(id2);
            update=1;
        end
        
    elseif (prob>=prange(2)) && (prob<prange(3))
        %pick two layers randomly and merge
        rp = randperm(k);
        l1 = rp(1); l2 = rp(2);
        layers_new(l1).k = [layers_new(l1).k layers_new(l2).k];
        ar = setdiff([1:k],l2);
        layers_new = layers_new(ar);
        
    elseif (prob>=prange(3)) && (prob<prange(4))
        %pick two layers randomly and swap their positions
        rp = randperm(k);
        l1 = rp(1); l2 = rp(2);
        
        swp = layers_new(l1).k;
        layers_new(l1).k = layers_new(l2).k;
        layers_new(l2).k = swp;
        
    elseif (prob>=prange(4))
        %pick a layer randomly and split it using max_cut
        pkd_lyr = ceil(k*rand);
        kk = layers_new(pkd_lyr).k;
        sub_adj = adj(kk,kk);
        sub_adj = sub_adj+sub_adj';
        
        layers_new1=layers_new;
        layers_new2=layers_new;
        
        ws = sum(sum(sub_adj));
        if ws>0 && length(kk)>2
            phase = max_cut(sub_adj);
            l1 = kk(phase(1).k1);
            l2 = kk(phase(1).k2);
        else
            nn = ceil(rand*length(kk));
            l1= kk(1:nn);
            l2= kk((nn+1):length(kk));
        end
        
        layers_new1(k+1).k = [];
        layers_new1(pkd_lyr).k = l1;
        
        layers_new2(k+1).k = [];
        layers_new2(pkd_lyr).k = l2;
        
        swp1_1 = l2; swp1_2 = l1;
        
        for jj = (pkd_lyr+1):(k+1)
            swp2_1 = layers_new1(jj).k;
            layers_new1(jj).k = swp1_1;
            swp1_1 = swp2_1;
            
            swp2_2 = layers_new2(jj).k;
            layers_new2(jj).k = swp1_2;
            swp1_2 = swp2_2;
        end
        
        kt1 = length(layers_new1);  kt2 = length(layers_new2);
        len = zeros(1,kt1);
        for jj = 1:kt1
            len(jj) = length(layers_new1(jj).k);
        end
        
        id = find(len==0);
        
        if ~isempty(id)
            id2 = setdiff(1:kt1,id);
            layers_new1 = layers_new1(id2);
            update = 1;
        end
        
        len = zeros(1,kt2);
        for jj = 1:kt2
            len(jj) = length(layers_new2(jj).k);
        end
        
        id = find(len==0);
        
        if ~isempty(id)
            id2 = setdiff(1:kt2,id);
            layers_new2 = layers_new2(id2);
            update = 1;
        end
        
        h1 = heir_ind_calc(Q_mat,N,layers_new1);
        h2 = heir_ind_calc(Q_mat,N,layers_new2);
        
        if h1>h2
            layers_new = layers_new1;
        else
            layers_new = layers_new2;
        end
    end
    
    heir_ind_new = heir_ind_calc(Q_mat,N,layers_new);
    delta_h = heir_ind_new-heir_ind;
    
    update = 0;
    
    if delta_h > 0
        
        layers = layers_new;
        heir_ind = heir_ind_new;
        update=1; failed=0;
    else
        prob2 = rand;
        
        if prob2 < exp(delta_h/T)
            layers = layers_new;
            heir_ind = heir_ind_new;
            update=1; failed=0;
        end
    end
    
    T = T * exp(-lambda);
    
    if update==0
        failed = failed+1;
    end
    
    if failed > 5*tot_nbr
        break;
    end
end

end