function mimoUplink{Ti<:Integer,Tf<:FloatingPoint}(Ns::Ti,Mc::Ti,
     M::Array{Ti,1},K::Array{Ti,1}, Ki::Array{Ti,1},Tt::Array{Ti,1},
     Td::Ti,rho::Array{Tf,1},gamma::Array{Tf,1}=[Inf],delta::Tf=0.0)
#mimoUplink(Ns,Mc,M,K,Tt,Td,rho)
#  Simulate uncoded transmission of QAM signals over a flat-fading Gaussian
#  uplink channel, for a base station with multiple receive
#  antennas. Constelations of size Mc are used with M receive antennas and K
#  txmitters, each with one transmit antenna. A power of "rho" per user is
#  utilized. Tt is the number of training samples to use, Td is the number
#  of data symbols, with Tt+Td being the training period. The resulting
#  symbol error rate is averaged over all users and Ns training periods.
#mimoUplink(Ns,Mc,M,K,Tt,Td,rho,gamma=Inf,delta=0)
#  If gamma is present, it is a CIR. If delta is present, it is
#  the pilot boost in db.
#
# Examples show SER as we deviate away from LTE PUSCH defaults
#  # Performance as a function of SNR
#  mimoUplink(100,4,4,[4],[0],[12],12,[-5.0:5:20;])
#  # Performance as a function of CIR
#  mimoUplink(400,4,4,[2],[1],[12],12,[15],[-15.0:5:25],0)
#  # Performance as a function of the number of antennas
#  mimoUplink(1000,4,[2:8],[2],[0],[4],8,[10])
#  # Performance as a function of Tt
#  mimoUplink(400,4,4,[2],[0],[2:2:16],12,[10])
#  # Performance as a function of the number of users
#  mimoUplink(400,4,4,[1:4],[0],[12],12,[10])
# Add an additional 'gamma' argument to the end of any of these functions to
# simulate with a single interferer of CIR gamma. 

# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Fri 14 Nov 14, 11:57am by cpeel

println("   Mc     M     K    Ki    Tt    Td   rho  gamma")

#algs = getAlgs()
algs = []

# It would be great to print out the alg names here
#for ax = 1:min(length(algs),4)
#    @printf("%10s ",algs[ax][2])
#end
#println

out = cell(0)
if length(rho)>1
  for s = 1:length(rho)
    push!(out,simQAM(Ns,Mc,M[1],K[1],Ki[1],Tt[1],Td,rho[s],gamma[1],delta,algs));
  end
  xval = rho;
#  xlab = L"\rho"* " (SNR in dB)";
  xlab = "rho (SNR in dB)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,CIR=%4.1f",
                   Mc,Ns,M[1],K[1],Ki[1],Tt[1],Td,gamma[1]);
elseif length(gamma)>1
  for s = 1:length(gamma)
    push!(out,simQAM(Ns,Mc,M[1],K[1],Ki[1],Tt[1],Td,rho[1],gamma[s],delta,algs));
  end
  xval = gamma;
  xlab = "\gamma (CIR in dB)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f",
                 Mc,Ns,M[1],K[1],Ki[1],Tt[1],Td,rho[1]);
elseif length(K)>1
  for s = 1:length(K)
    push!(out,simQAM(Ns,Mc,M[1],K[s],Ki[1],Tt[1],Td,rho[1],gamma[1],delta,algs));
  end
  xval = K;
  xlab = "K (# users)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,M[1],Ki[1],Tt[1],Td,rho[1],gamma[1]);
elseif length(Tt)>1
  for s = 1:length(Tt)
    push!(out,simQAM(Ns,Mc,M[1],K[1],Ki[1],Tt[s],Td,rho[1],gamma[1],delta,algs));
  end
  xval = Tt;
  xlab = "Tt (training length)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,M[1],K[1],Ki[1],Td,rho[1],gamma[1]);
elseif length(M)>1
  for s = 1:length(M)
    push!(out,simQAM(Ns,Mc,M[s],K[1],Ki[1],Tt[1],Td,rho[1],gamma[1],delta,algs));
  end
  xval = M;
  xlab = "M (# Rx antennas)";
  tstr = @sprintf("%dQAM,Ns=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,K[1],Ki[1],Tt[1],Td,rho[1],gamma[1]);
elseif length(K)==length(M) && length(K)>1
  for s = 1:length(K)
    push!(out,simQAM(NsMc,M[s],K[s],Ki[1],Tt[1],Td,rho[1],gamma[1],delta,algs));
  end
  xval = K;
  xlab = "M==K (# users, # antennas)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,M[1],Ki[1],Tt[1],Td,rho[1],gamma[1]);
elseif length(Mc)>1
  for s = 1:length(Mc)
    push!(out,simQAM(Ns,Mc[s],M[1],K[1],Ki[1],Tt[1],Td,rho[1],gamma[1],delta,algs));
  end
  xval = log2(Mc);
  xlab = "# constellation bits";
  tstr = "";
else
  out = simQAM(Ns,Mc,M[1],K[1],Ki[1],Tt[1],Td,rho[1],gamma[1],delta,algs);
end

pColor = {"r>-","bo--","kx-.","gd-","c^--","m*-.",
          "rs--","gp-.","bv-","kh--","c+-.","m.-"};
pIdx   = 1;

Ns = size(out,1);

    
#clf()
ser = zeros(Ns);
for a=1:length(out[1][2])
    for k=1:Ns
        ser[k] = out[k][1][a];
#        mi(k) = out(k).mi(a);
    end
#    figure(1)
    semilogy(xval,ser,pColor[pIdx],label=out[1][2][a]);
    hold(true);
##     figure(2)
##     plot(xval,mi,pColor[pIdx]);  hold(true);
    if pIdx == length(pColor)
         pIdx = 1;
    else
         pIdx = pIdx + 1;
    end
end

hold(false);
##figure(1)
xlabel(xlab);
ylabel("SER");
legend(loc=3);
grid(which="both",axis="y")
grid(which="major",axis="x")
title(tstr)

## figure(2)
## xlabel(xlab);
## ylabel("Mut Info");
## legend(out(1).name,4);
## grid on;
## title(tstr)

return
end

#######################################################################
function simQAM(Ns,Mc,M,K,Ki,Tt,Td,rdb,gdb,ddb,algs)

delta   = 10.^(ddb/10); # pilot boost
rho     = 10.^(rdb/10);
gamma   = 10.^(gdb/10);
if K>M
    error("Can""t have K>M!!")
end
T = Tt+Td;
Ka = K+Ki;
#
# Modulation class, modulated bits
#
Mp = convert(Int,round(sqrt(Mc)));
# Information bits, signals, Threshold.
Zreal = rand((0:Mp-1),K+Ki,Td*Ns);
Zimag = rand((0:Mp-1),K+Ki,Td*Ns);
Zi = Zreal*Mp+Zimag;
ZZ = Zi[1:K,:];
(C,d,Cr) = scodes(Mc,"QAM");
gam = sqrt(2/3*(Mc-1));
C = C[:].';

#
# Channel coefs, noise.
#
sigmaSq = 1;
rho     = rho;
# sigmaSq = 1/rho;
# rho     = 1;
HH = (randn(M,K+Ki,Ns) + im*randn(M,K+Ki,Ns))*sqrt(1/2);# Channel
HH[:,1:K,:] = HH[:,1:K,:]*sqrt(rho);                    # Signal power
HH[:,K+1:Ka,:] = HH[:,K+1:Ka,:]*sqrt(rho/gamma);        # Interference power
NN = (randn(M,T*Ns) + im*randn(M,T*Ns))*sqrt(sigmaSq/2);

#
# Training signal
#
St = repmat(eye(K),1,round(Int,ceil(Tt/K)));
St = St[:,1:Tt]*sqrt(delta)
SSi = sign(rand(Ki,Tt,Ns)-.5)/sqrt(Tt);

onesKT = ones(K,Td);
ri = [1:K;];
ii = [(K+1):(2*K);];
ria = [1:K;]';
iia = [(Ka+1):(Ka+K);]';

Nalgs = 5;
ZD = zeros(K,Td*Ns,Nalgs);
algNames = cell(Nalgs)

for ix = 1:Ns
    zix = (ix-1)*Td+1:ix*Td;
    N = NN[:,(ix-1)*T+1:ix*T];
    Z = Zi[:,zix];
    Si = SSi[:,:,ix];
    Sti = [St; Si];
    Sd = C[Z+1];
    S = [Sti Sd];
    H = HH[:,:,ix];
    
    Y = H*S+N; # Channel
    Yt = Y[:,1:Tt];
    Yd = Y[:,Tt+1:Tt+Td];
    Ydc = [real(Yd); imag(Yd)];
    #
    # Genie channel
    #
    Hdes = H[:,1:K];
    Hc = [real(Hdes) -imag(Hdes); imag(Hdes) real(Hdes)];
    HcAll = [real(H) -imag(H); imag(H) real(H)];
    #
    # Training
    #
    Hhat = Yt*St'*inv(St*St');
    Hhatc = [real(Hhat) -imag(Hhat); imag(Hhat) real(Hhat)];
    #
    # additional Training
    #
    Hhat2= Hhat;
    #
    What = (Yt-Hhat2*St);
    Rhat = What*What';
    var2hat = sum(diag(Rhat))/(Tt-1)/M;
    powhat = sum(diag(Hhat2*Hhat2'))/K/M;
    Hhat2 = Yt*St'*inv(St*St'+ I*var2hat/powhat*10);
    

    ax=0;

    ax = ax+1;
    algNames[ax] = "Zero Forcing";
    Hhat = Yt*St'*inv(St*St');
    Wls = Hhat/(Hhat'*Hhat);
    Zd  = mimo_slice(Wls'*Yd,onesKT,C);
    ZD[:,zix,ax] = Zd;

    ax = ax+1;
    algNames[ax]= "MMSE LLL"
    (Bw, Tw, t1, t2) = lll([Hhat; sigmaSq*complex(eye(K))]);
    Ti = inv(Tw);
    Wlll = Bw/(Bw'*Bw);
    Ztmp = Wlll'*[Yd; zeros(K,Td)];
    Stmp = roundc((Ztmp*gam-(1+im)*Ti*onesKT)/2)*2 +(1+im)*Ti*onesKT;
    Zd = mimo_slice(Tw*Stmp/gam,onesKT,C);
    ZD[:,zix,ax] = Zd;


    ax = ax+1;
    algNames[ax] = "MMSE-VBLAST"
    (W,P,B) = vblast(Hhat,1/sigmaSq/K);
    Yp = W*Yd;
    Sp = zeros(Complex{Float64},K,Td);
    for tt=1:Td
        for kx=1:K
            Sp[kx,tt] = siso_demod(Yp[kx,tt] -B[kx,:]*Sp[:,tt],C);
        end
    end
    ZD[:,zix,ax] = mimo_slice(P*Sp,onesKT,C);
    
    ax = ax+1;
    algNames[ax] = "IRC";
    Hhat = (Yt*St')/(St*St');
    Wp = Yt-Hhat*St;     # To be used to estimate the noise variance
    Rhat = Wp*Wp'/(max(Tt-K,1)*M);
    Ri = inv(Rhat + trace(Rhat)/1e14*complex(eye(M)));
    W = Ri*Hhat/(Hhat'*Ri*Hhat);
    Zd = mimo_slice(W'*Yd,W'*Hhat*onesKT,C);
    ZD[:,zix,ax] = Zd;

    ax = ax+1;
    algNames[ax] = "Sphere Decoder";
    z_est = hard_sphere(Ydc,Hhatc/gam,Mp);
    Zd = z_est[ri,:]*Mp + z_est[ii,:]
    ZD[:,zix,ax] = Zd;

end
@printf("%5d %5d %5d %5d %5d %5d %5.1f %5.1f",
          Mc,  M,  K, Ki, Tt, Td,  rdb,  gdb)


ser = zeros(Nalgs);
#Calculate symbol error rate
for ax = 1:Nalgs
    ser[ax] = mean(map(!=,ZZ,ZD[:,:,ax]));
    if ax<6
        @printf("    %7.4f",ser[ax])
    end 
end
@printf("\n")

return ser, algNames
end
#######################################################################
