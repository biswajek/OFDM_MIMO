function [ofdmSymbolWithCP] = genOFDMSymb(ofdmParams,genPreamble)
  % K - Number of Subcarriers
  % Kdata - Number of Subcarriers fors data
  % CP - Cyclic Prefix 128 fixed
  % mu - modulation 4 fixed
  % genPreamble - flag to indicate if the symbol is preamble


  %mu =4;
  %Generate Payload bits
  k = log2(ofdmParams.mu);
  payloadBits_per_OFDM = ofdmParams.Kdata*k;

  dataIn = randi([0 1],1,payloadBits_per_OFDM);
  %Symbol mapping and modulation
  if genPreamble == 1
    dataInMod = 2*qammod(bi2de(reshape(dataIn,payloadBits_per_OFDM/k,k)),ofdmParams.mu,'gray');
    dataInMod(1:2:end) = 0;
  else
    dataInMod = qammod(bi2de(reshape(dataIn,payloadBits_per_OFDM/k,k)),ofdmParams.mu,'gray');
  endif

  freqDomData = zeros(1,ofdmParams.K);
  offset = (ofdmParams.K - ofdmParams.Kdata)/2;
  freqDomData(offset:offset+length(dataInMod)-1) = dataInMod;

  freqDomData = fftshift(freqDomData);
  symbol = ifft(freqDomData)*sqrt(ofdmParams.K);
  cp = symbol(end-ofdmParams.CP:end);
  ofdmSymbolWithCP = [cp symbol];

