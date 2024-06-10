% This program splits a given symmetric weighted network into two
% sub-networks such that the total connection weight between the two
% sub-networks is the maximum. The algorithm is modified from the Min-cut
% algorithm by Stoer and Wagner, J. ACM, 1997.
% We can call it 'Max-cut'.
% coded by - Anand Pathak 
% ==========================================================================================

function phase = max_cut(adj)

N=length(adj(:,1));

a = ceil(rand*N); %selecting the starting node for expansion of A in every expansion phase


for ii = 1:N
    super_node(ii).k=ii;
end
m=adj;

Ns=length(super_node);
ph = 0;

wmax=0;

while Ns>2 %loop for every expansion phase
    
    ph=ph+1;

for ii = 1:Ns
    if super_node(ii).k==a
        nuc=ii;
    end
end
clust=nuc;

for jj=1:Ns-2   %loop for a single expansion phase 
links = m(:,clust);
col = sum(links,2);
[b,ind]=sort(col);
cc=1;ii=0;

while cc==1    %to ensure the weakly connected node is outside the expanding group A
    ii=ii+1;
    node_next=ind(ii);
    cc=0;
    if ismember(node_next,clust)
        cc=1;
    end
end

clust=[clust;node_next]; %expansion of A

if jj==Ns-2
    s=node_next;
    t=setdiff(1:Ns,clust);  %last two nodes added in the phase
end

end

w = sum(m(t,:));

if w>wmax

  phase(1).k1 = super_node(t).k;
  phase(1).k2 = setdiff(1:N,super_node(t).k);
  phase(1).w = w;                        %cut at the end of the phase
  wmax = w;
end


%combining s-t into single node for next phase
ss=min(s,t);tt=max(s,t);
super_node(ss).k = [super_node(ss).k super_node(tt).k];
ar=setdiff(1:Ns,tt);
super_node = super_node(ar);

m(ss,:)=m(ss,:)+m(tt,:); m(:,ss)=m(:,ss)+m(:,tt);
m(ss,ss)=0;
m = m(ar,ar); 

Ns=length(super_node);
end

end