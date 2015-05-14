

roundc(r) = round(real(r)) + im*round(imag(r));
#round{Td}(r::Td) = Td<:Complex?round(real(r)) + im*round(imag(r)):round(r);

#################################################################
function scodes(M,Ctype)
#C = scodes(M,Ctype)
#  Return a constellation of M signals, where <Ctype> is 'PSK', 'QAM',
#  or 'PAM'.
#[C,d] = scodes(M,Ctype)
#  Also return the minimum distance between constellation points.
#[C,d,Cr] = scodes(M,Ctype)
#  Returns the underlying real PAM constellation for square (M=m^2)
#  QAM constellations

# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Wed 13 May 15, 12:55pm
#

#if nargin==0
#  error("Bad arguments. Try again.")
#end
if Ctype=="PSK"
    C=exp(im*2*pi*[1:M]'/M);
elseif Ctype=="QAM"
    if abs(sqrt(M)-round(sqrt(M)))< 1e-3
        # Make QAM constellations if M is a square
        m = round(Int,sqrt(M))
        C = complex(zeros(m,m))
        Cr = complex([-(m-1):2:(m-1);])
        for ix = 1:m
            for jx = 1:m
                C[jx,ix] = Cr[ix] +im*Cr[jx];
            end
        end
    else
        error("Can only handle square QAM constellations.")
        #     [x y] = qaskenco(M);
        #     C = x+im*y;
    end
elseif Ctype=="PAM"
    C = [-(M-1):2:(M-1);];
    Cr = C;
else
  error("Bad args.");
end

C = C[:]; # Turn into column vector
# for QAM, gam = sqrt(2/3*(M-1))
gam = sqrt(C'*C/length(C));

C = C/gam;
#if nargout ==3
  Cr = Cr/gam;
#end

d2 = real(C[2:end]-C[1]);
id = find(x->abs(x)>0,d2);
d = minimum(abs(d2[id]));
return (C, d, Cr, gam)
end

#######################################################################
function mimo_slice(Y,H,C)
M = size(Y,1);
Td = size(Y,2);
Zd = zeros(M,Td);
Mc = length(C);
dy = zeros(M,Mc);
for t=1:Td
    for i = 1:Mc
        dd = Y[:,t]-H[:,t]*C[i];
        dy[:,i] = abs(dd).^2
    end
    for m = 1:M
        (mn,zd) = findmin(dy[m,:],2);
        Zd[m,t] = zd[1]-1.0;
    end
end
return Zd
end

#################################################################
function siso_demod(y,C)
Mc = length(C);
dy = zeros(Mc,1);
for ix = 1:Mc
    dy[ix] = abs(y[1]-C[ix])
end
(m,zd) = findmin(dy);
sd = C[zd];
return sd[1]
end

