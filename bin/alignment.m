function [G,cvx_status] = alignment(movie)

[L,time] = size(movie);

Y = vec(movie) * (vec(movie))';

idx = @(i,k) (i-1)*L+k;

cvx_begin quiet
    cvx_solver sedumi
    variable G(L*time,L*time) semidefinite
    maximize(trace(Y*G))
    subject to
        for i = 1:time
            G(idx(i,1):idx(i,L),idx(i,1):idx(i,L)) == eye(L);
        end
        for i = 2:time
            for j = 1:i-1
                for k = 1:L
                    sum(G(idx(i,k),idx(j,1):idx(j,L))) == 1;
                    for l = 1:(L-1)
                        G(idx(i,mod(k-1+l-1,L)+1),idx(j,l)) - G(idx(i,mod(k-1+l,L)+1),idx(j,l+1)) == 0;
                    end
                    for l = 1:L
                        G(idx(i,k),idx(j,l)) >= 0;
                    end
                end
            end
        end 
cvx_end 

end
