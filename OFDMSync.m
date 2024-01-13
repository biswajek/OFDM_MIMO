function [ofdmFrame] = OFDMSync(ofdmParams)
  % Refer to paper by Schimdl and Cox
  % K - Number of Subcarriers
  % Kdata - Number of Subcarriers fors data
  % CP - Cyclic Prefix
  % mu - modulation
  % ofdmSymbPerFrame - Number of OFDM symbols in the frame
  % numFrames - number of frames of OFDM symvol
  % L - Length of half preamble

  % Time offset to simulate the path delay
  timeOffset = ofdmParams.K/2;

  ofdmSuperFrame = [];
  payloadSuperFrame = [];
  preambleSuperFrame = [];

  for frame = 1:ofdmParams.numFrames
    genPreamble = 1;
    ofdmFramePreamble = genOFDMSymb(ofdmParams,genPreamble);
    ofdmFrame = ofdmFramePreamble;

    genPreamble = 0;
    for i = 1:ofdmParams.ofdmSymbPerFrame - 1
      ofdmFrame = [ofdmFrame genOFDMSymb(ofdmParams,genPreamble)];
    end

    %Add static time offset and noise
    %ofdmFrame = [zeros(1,timeOffset) ofdmFrame];
    ofdmFrame = awgn(ofdmFrame,10);

    % Used for plots
    preamble = zeros(1,length(ofdmFrame));
    payload = zeros(1,length(ofdmFrame));
    preamble(ofdmParams.CP :ofdmParams.CP + ofdmParams.K) = 1;
    for i = 1:ofdmParams.ofdmSymbPerFrame - 1
      payload((i+1)*ofdmParams.CP + i*ofdmParams.K: (i+1)*ofdmParams.CP+(i+1)*ofdmParams.K) = 1;
    end

    ofdmSuperFrame = [ofdmSuperFrame ofdmFrame];
    payloadSuperFrame = [payloadSuperFrame payload];
    preambleSuperFrame = [preambleSuperFrame preamble];

  endfor

  temp = zeros(1,timeOffset);
  temp =awgn(temp,10);
  ofdmSuperFrame = [temp ofdmSuperFrame];

  payloadSuperFrame = [zeros(1,timeOffset) payloadSuperFrame];
  preambleSuperFrame = [zeros(1,timeOffset) preambleSuperFrame];

  %Calculate P(d) as given in the paper
  Pd = zeros(1, length(ofdmSuperFrame)-2*ofdmParams.L);
  for d = 1:length(ofdmSuperFrame)-2*ofdmParams.L
    Pd(d) = sum(conj(ofdmSuperFrame(d:d+ofdmParams.L)).*ofdmSuperFrame((d+ofdmParams.L):(d+2*ofdmParams.L)));
  endfor

  %Calculate R(d) as given in the paper
  Rd = zeros(1, length(ofdmSuperFrame)-2*ofdmParams.L);
  for d = 1:length(ofdmSuperFrame)-2*ofdmParams.L
    Rd(d) = sum(abs(ofdmSuperFrame((d+ofdmParams.L):(d+2*ofdmParams.L))).^2);
  endfor

  M = abs(Pd).^2./Rd.^2;
  plot(abs(ofdmSuperFrame),"c");
  hold on;
  plot(preambleSuperFrame,"k:");
  plot(payloadSuperFrame,"b");
  plot(M,"b--");



