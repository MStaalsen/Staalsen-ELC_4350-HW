function s = rectPulse(syms, P, t_off);

% s=srrc(syms, beta, P, t_off);
% Generate a Square-Root Raised Cosine Pulse
%      'syms' is 1/2 the length of srrc pulse in symbol durations
%      'beta' is the rolloff factor: beta=0 gives the sinc function
%      'P' is the oversampling factor
%      't_off' is the phase (or timing) offset

if nargin==2, t_off=0; end;                       % if unspecified, offset is 0
k=-syms*P+1e-8+t_off:syms*P+1e-8+t_off;           % sampling indices as a multiple of T/P
%s=rect(k,P/2);
s=rectpuls(k,P/4);