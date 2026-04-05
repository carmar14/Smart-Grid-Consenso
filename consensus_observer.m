function xhat_i = consensus_observer(i,xhat_all,neighbors,A,B,u,k0,Pi)

sum_term = zeros(size(xhat_all));

for j = neighbors
    sum_term = sum_term + (xhat_all(:,j) - xhat_all(:,i));
end

local_term = Pi*(xhat_all(:,i) - xhat_all(:,i));

xhat_i = A*xhat_all(:,i) + B*u + k0*(sum_term + local_term);

end