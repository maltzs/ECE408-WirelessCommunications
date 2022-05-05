% Simulates NrvxNtx Rayleigh channels of length N assuming a max Doppler
% shift of fm. Algorithm from "Wireless Communications: Principles and
% Practice," 1st Edition by Rappaport.
function r = rayleigh(fm, N, Nrv, Ntx)
    f = linspace(-fm,fm,N)';

    g = randn(Nrv,Ntx,N/2) + 1j*randn(Nrv,Ntx,N/2);
    g = cat(3,conj(flip(g,3)),g);
    
    % Fading spectrum, uses linearized values at S(1) and S(end) due to
    % infinity values.
    S = 1.5./(pi*fm*sqrt(1-(f(2:end-1)/fm).^2));
    dS = S(2)-S(1);
    S = [S(1)-dS; S; S(end)-dS];
    
    g = reshape(sqrt(S),1,1,[]).*g;
    gi = ifft(real(g),[],3);
    gq = ifft(imag(g),[],3);
    
    % Rayleigh is magnitude of complex Gaussian random variable
    r = sqrt(abs(gi).^2+abs(gq).^2);
end