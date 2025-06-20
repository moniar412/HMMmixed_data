function p = bvnlr( dh, dk, r )
%BVNLR
%  A function for computing bivariate normal probabilities.
%  bvnlr calculates the probability that x < dh and y < dk.
%    parameters
%      dh 1st upper integration limit
%      dk 2nd upper integration limit
%      r   correlation coefficient
%  Example:
%    p = bvnlr( 3, 1, .35 )
%

%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%   This function is based on the method described by
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
%
dh=-dh;
dk=-dk;
if dh ==  inf | dk ==  inf, p = 0;
elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
elseif dk == -inf, p = phid(-dh);
else
    if abs(r) < 0.3
        %       Gauss Legendre points and weights, n =  6
        w = [0.1713244923791705 0.3607615730481384 0.4679139345726904]';
        x = [0.9324695142031522 0.6612093864662647 0.2386191860831970]';
    elseif abs(r) < 0.75
        %       Gauss Legendre points and weights, n = 12
        w = [.04717533638651177 0.1069393259953183 0.1600783285433464 ...
            0.2031674267230659 0.2334925365383547 0.2491470458134029]';
        x = [0.9815606342467191 0.9041172563704750 0.7699026741943050 ...
            0.5873179542866171 0.3678314989981802 0.1252334085114692]';
    else
        %       Gauss Legendre points and weights, n = 20
        w = [.01761400713915212 .04060142980038694 .06267204833410906 ...
             .08327674157670475 0.1019301198172404 0.1181945319615184 ...
            0.1316886384491766 0.1420961093183821 0.1491729864726037 ...
            0.1527533871307259]';
        x = [0.9931285991850949 0.9639719272779138 0.9122344282513259 ...
             0.8391169718222188 0.7463319064601508 0.6360536807265150 ...
             0.5108670019508271 0.3737060887154196 0.2277858511416451 ...
             0.07652652113349733]';
     end
    h = dh; k = dk; hk = h*k; bvn = 0;
    if abs(r) < 0.925 
        hs = ( h*h + k*k )/2; 
        asr = asin(r);
        sn = sin( asr*( 1 - x)/2 );
        bvn = w'*exp( ( sn*hk - hs )./( 1 - sn.*sn ) );
        sn = sin( asr*( 1 + x )/2 );
        bvn = bvn + w'*exp( ( sn*hk - hs )./( 1 - sn.*sn ) );
        bvn = bvn*asr/( 4*pi );
        bvn = bvn + phid(-h)*phid(-k);
    else
        twopi = 2*pi;
        if r < 0
            k = -k; 
            hk = -hk;
        end
        if abs(r) < 1
            as = ( 1 - r )*( 1 + r ); a = sqrt(as); 
            bs = ( h - k )^2;
            c = ( 4 - hk )/8 ; d = ( 12 - hk )/16; 
            asr = -( bs/as + hk )/2;
            if asr > -100
                bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 );
            end
            if hk > -100, b = sqrt(bs); sp = sqrt(twopi)*phid(-b/a);
                bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*( 1 - d*bs/5 )/3 );
            end
            a = a/2;
            xs = ( a - a*x ).^2;
            asr = -( bs./xs + hk )/2;
            ind = asr > -100;
            xs=xs(ind);
            rs = sqrt( 1 - xs );
            sp = ( 1 + c*xs.*( 1 + d*xs ) );
            ep = exp( -hk*( 1 - rs )./( 2*( 1 + rs ) ) )./rs;
            bvn = bvn + a*w(ind)'*(exp(asr(ind)).*( ep - sp ));
            xs = ( a + a*x ).^2;
            asr = -( bs./xs + hk )/2;
            ind = asr > -100;
            xs=xs(ind);
            rs = sqrt( 1 - xs );
            sp = ( 1 + c*xs.*( 1 + d*xs ) );
            ep = exp( -hk*( 1 - rs )./( 2*( 1 + rs ) ) )./rs;
            bvn = bvn + a*w(ind)'*(exp(asr(ind)).*( ep - sp ));
            bvn = -bvn/twopi;
        end
        if r > 0, bvn =  bvn + phid( -max( h, k ) ); end
        if r < 0, bvn = -bvn + max( 0, phid(-h)-phid(-k) ); end
    end
    p = max( 0, min( 1, bvn ) );
end
%
%
function p = phid(z), p = erfc( -z/sqrt(2) )/2; % Normal cdf
%
% end phid