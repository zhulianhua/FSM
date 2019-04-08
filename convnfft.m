function C = convnfft(A, B)
    % assuming size(A) = size(B) = [n, n];
    n = size(A,1);
    % The circular convolution size
    l = 2^nextpow2(2*n-1);
    % zero padding in each dimension to l before fft
    A = fftn(A,[l,l,l]);
    B = fftn(B,[l,l,l]);
    C = A.*B;  
    C = ifftn(C);
    % truncated the central mode, size(C) will be the same
    %   as size(A)
    t = ceil((n-1)/2)+(1:n);
    % using t = 1:l, the results will be identical to convn
    C = C(t,t,t);
end % convn