K=64;
CP = K/4;
P=8;
pilotValue = 3+3j;
allCarriers = 0:1:K-1;
pilotCarriers = allCarriers(1:K/P:end);
pilotCarriers = [pilotCarriers allCarriers(end)];
pilotCarriers = [pilotCarriers;zeros(7,9)];
pilotCarriers = pilotCarriers(:)';
pilotCarriers = pilotCarriers(1:64);
pilotCarriers(end) = 63;
dataCarriers = allCarriers - pilotCarriers;

pilotCarriers = [0 pilotCarriers(find(pilotCarriers(:)~=0))];
dataCarriers = dataCarriers(find(dataCarriers(:)~=0));

mu=16; %16QAM
k = log2(mu);
payloadBits_per_OFDM = length(dataCarriers)*k;

dataIn = randi([0 1],1,payloadBits_per_OFDM);
%Symbol mapping and modulation
dataInMod = qammod(bi2de(reshape(dataIn,payloadBits_per_OFDM/k,k)),mu,'gray');

symbol = zeros(1,K);
symbol(pilotCarriers+1) = pilotValue;
symbol(dataCarriers+1) = dataInMod;

ofdmTD = ifft(symbol);

% Channel Impulse response
channelResponse = [1 0  0.3+0.3j];
_H = fft(channelResponse, K);
snrIndB = 30;

cp = ofdmTD(end-CP:end);
ofdmTDWithCP = [cp ofdmTD];

ofdmTx = ofdmTDWithCP;
temp = conv(ofdmTx,channelResponse);
signalPower = mean(abs(temp).^2);
noiseVar = signalPower*10^(-snrIndB/10);
ofdmRx = temp + sqrt(noiseVar/2).*(randn(size(temp))+1j*randn(size(temp)));

%Remove CP
ofdmRxNoCP = ofdmRx(CP+1:CP+K);
ofdmDemod = fft(ofdmRxNoCP);

%Channel Estimate and interpolation
%pilots = zeros(1,length(allCarriers));
pilots = ofdmDemod(pilotCarriers+1);
_HEstPilots = pilots/pilotValue;

_HEstMag = spline(pilotCarriers+1, abs(_HEstPilots), allCarriers);
_HEstPhase = spline(pilotCarriers+1, angle(_HEstPilots), allCarriers);

_HEst = _HEstMag .* exp(1j * _HEstPhase);

ofdmEqualized = ofdmDemod ./_HEst;

%Euclidean Distance
x = [real(ofdmEqualized);imag(ofdmEqualized)]';
y = [real(symbol);imag(symbol)]';
Xm = sum(x.*x,2);
Ym = sum(y.*y,2)';
d = Xm(:,ones(1,size(y)(1)))+Ym(ones(1,size(x)(1)),:)-2*x*y'
[~, indices] = min(d,[],2);

harddecisions = symbol(indices);
scatterplot(harddecisions);


