function out = bandpass2(low_cutoff, high_cutoff, samp_rate, inp)
% BANDPASS applies bandpass filter to a signal
%   OUT = BANDPASS(LOW_CUTOFF, HIGH_CUTOFF, SAMP_RATE, INP)
%   LOW_CUTOFF, HIGH_CUTOFF and SAMP_RATE are in Hz. LOW_CUTOFF can be
%   zero. INP is the vector containing the input signal sampled at
%   SAMP_RATE. INP can be line or column vector.

finp = fft(inp);

nlow_cutoff = round(low_cutoff * (length(inp)/samp_rate));
nhigh_cutoff = round(high_cutoff * (length(inp)/samp_rate));

freq_flt = zeros(size(inp));
freq_flt(nlow_cutoff+1:min(nhigh_cutoff+1,end)) = 1;
freq_flt(end-nhigh_cutoff+1:min(end-nlow_cutoff+1,end)) = 1;

out = ifft(finp.*freq_flt);
if (~isreal(out))
    error('Bug in the BANDPASS routine - the output is complex');
end