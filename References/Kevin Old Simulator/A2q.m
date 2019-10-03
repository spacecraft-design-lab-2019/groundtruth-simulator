function q = A2q(A)
% NOTE: we always put the scalar part at the end!
    % use stanley's method to compute DCM.
    % I.e., first find the largest qi value naively,
    % then use that one to compute the others

    % naive computation:
    tr = trace(A);
    qs = 0.5 * sqrt(1 + tr);
    qv = 0.5 * sqrt(1 + 2*diag(A) - tr);
    % select the maximum value from those:
    [maxval, ind] = max(abs([qv; qs]));

    % Based on which is the max, compute the rest:
    % see Table in section 6 here:
    % http://www.tu-berlin.de/fileadmin/fg169/miscellaneous/Quaternions.pdf
    denom = 4*maxval;
    if ind == 1      % qx
        qs = (A(2,3) - A(3,2)) / denom;
        qx = maxval;
        qy = (A(2,1) + A(1,2)) / denom;
        qz = (A(1,3) + A(3,1)) / denom;
    elseif ind == 2  % qy
        qs = (A(3,1) - A(1,3)) / denom;
        qx = (A(2,1) + A(1,2)) / denom;
        qy = maxval;
        qz = (A(3,2) + A(2,3)) / denom;
    elseif ind == 3  % qz
        qs = (A(1,2) - A(2,1)) / denom;
        qx = (A(1,3) + A(3,1)) / denom;
        qy = (A(3,2) + A(2,3)) / denom;
        qz = maxval;
    elseif ind == 4  % qs
        qs = maxval;
        qx = (A(2,3) - A(3,2)) / denom;
        qy = (A(3,1) - A(1,3)) / denom;
        qz = (A(1,2) - A(2,1)) / denom;
    end

    q = [qx; qy; qz; qs];
end