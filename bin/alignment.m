function [X,cvx_status] = alignment(movie,tol)

[L,time] = size(movie);

Copt = movie(:) * (movie(:))';

idx = @(i,k) (i-1)*L+k;

cvx_begin quiet
    cvx_solver sedumi
    % See http://cvxr.com/cvx/doc/solver.html for a detailed description of
    % how to control the tolerance in CVX
    if tol == 1
        cvx_precision low
    elseif tol == 2
        cvx_precision default
    else
        cvx_precision high
    end
    variable X(L*time,L*time) semidefinite
    maximize(trace(Copt*X))
    subject to
        for i = 1:time
            X(idx(i,1):idx(i,L),idx(i,1):idx(i,L)) == eye(L);
        end
        for i = 2:time
            for j = 1:i-1
                for k = 1:L
                    sum(X(idx(i,k),idx(j,1):idx(j,L))) == 1;
                    for l = 1:(L-1)
                        X(idx(i,mod(k-1+l-1,L)+1),idx(j,l)) - X(idx(i,mod(k-1+l,L)+1),idx(j,l+1)) == 0;
                    end
                    for l = 1:L
                        X(idx(i,k),idx(j,l)) >= 0;
                    end
                end
            end
        end 
cvx_end 

end
