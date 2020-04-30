function C = convnfftzhu1D(A, B)
    ABreal = isreal(A) && isreal(B);
    n = size(A,1);

    % The circular convolution size
    l = 2^nextpow2(2*n-1);
    % [l,l,l] means zero padding to [l,l,l]
    A = fft(A,l);
    B = fft(B,l);

    C = A(:).*B(:);

    C = ifft(C);
    
    % Truncate the results
    % valid
    t = ceil((n-1)/2)+(1:n);
    % full
    t = 1:l;
    % same
    t = 1:2*n-1;
    
    if ABreal
        C = real(C(t));
    else
        C = C(t);
    end
end % convnfft