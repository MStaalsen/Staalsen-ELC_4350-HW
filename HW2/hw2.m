%% idsys.m:  idealized transmission system

%TRANSMITTER
% encode text string as T-spaced 4-PAM sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=20;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r=c.*x;                       % modulate message with carrier

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters 
damps=[1 1 0 0 ]; 
b=firpm(fl,fbe,damps);        % create LPF impulse response

x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented 
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % plot soft decisions
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
% decode decision device output to text string
reconstructed_message=pam2letters(mprime)

%% Exercise 9.1. 
% Using idsys.m, examine the effect of using different carrier frequencies. 
% Try fc = 50, 30, 3, 1, 0.5. What are the limiting factors that cause
% some to work and others to fai

fc_s = [50,30,3,1,0.5]

for fc_i = 1:5
    fc = fc_s(fc_i)

    %TRANSMITTER
    % encode text string as T-spaced 4-PAM sequence
    str='01234 I wish I were an Oscar Meyer wiener 56789';
    m=letters2pam(str); N=length(m); % 4-level signal of length N
    % zero pad T-spaced symbol sequence to create upsampled
    % T/M-spaced sequence of scaled T-spaced pulses (T=1)
    M=100;                        % oversampling factor
    mup=zeros(1,N*M);             % Hamming pulse filter with 
    mup(1:M:N*M)=m;               % T/M-spaced impulse response
    p=hamming(M);                 % blip pulse of width M
    x=filter(p,1,mup);            % convolve pulse shape with data
    figure(1+(2*(fc_i-1))), plotspec(x,1/M)    % baseband AM modulation
    t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
    %fc=20;                        % carrier frequency
    c=cos(2*pi*fc*t);             % carrier
    r=c.*x;                       % modulate message with carrier

    %RECEIVER
    % am demodulation of received signal sequence r
    c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
    x2=r.*c2;                      % demod received signal
    fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters 
    damps=[1 1 0 0 ]; 
    b=firpm(fl,fbe,damps);        % create LPF impulse response

    x3=2*filter(b,1,x2);           % LPF and scale signal
    % extract upsampled pulses using correlation implemented 
    % as a convolving filter; filter with pulse and normalize
    y=filter(fliplr(p)/(pow(p)*M),1,x3);
    % set delay to first symbol-sample and increment by M
    z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
    figure(2+(2*(fc_i-1))), plot([1:length(z)],z,'.') % plot soft decisions
    % decision device and symbol matching performance assessment
    mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
    cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
    lmp=length(mprime);
    pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
    % decode decision device output to text string
    reconstructed_message=pam2letters(mprime)

end


%this works for all the tested choices of fc except fc = 0.5. This is
%because a carrier frequency of 0.5 implies that the message is being
%modulated with a frequency which is lower than the content of the signal,
%resulting in aliasing and ultimately the corruption of the signal.

%% Exercise 9.2. Using idsys.m, examine the effect of using different oversampling
% frequencies. Try M = 1000, 25, 10. What are the limiting factors that
% cause some to work and others to fail?

M_s = [1000, 25, 10]

for ndx = 1:3
    M = M_s(ndx)

    %TRANSMITTER
    % encode text string as T-spaced 4-PAM sequence
    str='01234 I wish I were an Oscar Meyer wiener 56789';
    m=letters2pam(str); N=length(m); % 4-level signal of length N
    % zero pad T-spaced symbol sequence to create upsampled
    % T/M-spaced sequence of scaled T-spaced pulses (T=1)
    %M=100;                        % oversampling factor
    mup=zeros(1,N*M);             % Hamming pulse filter with 
    mup(1:M:N*M)=m;               % T/M-spaced impulse response
    p=hamming(M);                 % blip pulse of width M
    x=filter(p,1,mup);            % convolve pulse shape with data
    figure(1+(2*(ndx-1))), plotspec(x,1/M)    % baseband AM modulation
    t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
    fc=20;                        % carrier frequency
    c=cos(2*pi*fc*t);             % carrier
    r=c.*x;                       % modulate message with carrier

    %RECEIVER
    % am demodulation of received signal sequence r
    c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
    x2=r.*c2;                      % demod received signal
    fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters 
    damps=[1 1 0 0 ]; 
    b=firpm(fl,fbe,damps);        % create LPF impulse response

    x3=2*filter(b,1,x2);           % LPF and scale signal
    % extract upsampled pulses using correlation implemented 
    % as a convolving filter; filter with pulse and normalize
    y=filter(fliplr(p)/(pow(p)*M),1,x3);
    % set delay to first symbol-sample and increment by M
    z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
    figure(2+(2*(ndx-1))), plot([1:length(z)],z,'.') % plot soft decisions
    % decision device and symbol matching performance assessment
    mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
    cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
    lmp=length(mprime);
    pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
    % decode decision device output to text string
    reconstructed_message=pam2letters(mprime)

end

% the case of M = 10 does not work. The figure outputs show that signal
% values recovered from modulation are not tightly grouped around the
% values they were transmitted at. In addition, the recovered signal is not
% entirely corrupted: you can still make out the word "Meyer" in the
% output. This indicates that the some recovered values are different
% enough from their original values to be incorrectly quantized. I believe
% this is caused by the lower oversampling frequency effectively producing
% a "low resolution" image of the signal output, which is then upconverted
% without some of the high frequency information in the signal. This would
% result in noise in the signal.

%% Exercise 9.3. What happens if the LPF at the beginning of the receiver is
%removed? What do you think will happen if there are other users present? Try
%adding in “another user” at fc = 30.

%TRANSMITTER 1
% encode text string as T-spaced 4-PAM sequence
str='keldifnt97*jk5p[;ldr21nbaztqfdvclptnvtq3d$njj#5';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=30;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r1=c.*x;                       % modulate message with carrier


%TRANSMITTER 2
% encode text string as T-spaced 4-PAM sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=20;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r2=c.*x;                       % modulate message with carrier

r = 20*r1+r2;

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters 
damps=[1 1 0 0 ]; 
b=firpm(fl,fbe,damps);        % create LPF impulse response

x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented 
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
y2=filter(fliplr(p)/(pow(p)*M),1,2*x2); %y2 has not been through the LPF
% set delay to first symbol-sample and increment by M
z=y2(0.5*fl+M:M:N*M);           % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % plot soft decisions
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
% decode decision device output to text string
reconstructed_message=pam2letters(mprime)


%adding another user at fc = 30 introduces some additional noise when the
%LPF is not in place at the receiver, but seemingly not enough to prevent
%the quantizer from picking up its intended transmission when the two
%signals are of equal power. In this case, the operation of correlating the
%received signal with the hamming pulse is enough to remove the higher
%frequency components introduced by the out-of-band signal. However, if
%this second signal's power is increased by multiplying it with a scalar
%while adding it to the in band signal, the amplitude of the noise which
%persists even after the pulse correlation process is enough to cause
%errors if the LPF is not in place.

%% Exercise 9.4. What are the limits to the LPF design at the beginning of the
%receiver? What is the lowest cutoff frequency that works? The highest?


%TRANSMITTER 1
% encode text string as T-spaced 4-PAM sequence
str='keldifnt97*jk5p[;ldr21nbaztqfdvclptnvtq3d$njj#5';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=30;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r1=c.*x;                       % modulate message with carrier


%TRANSMITTER 2
% encode text string as T-spaced 4-PAM sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=20;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r2=c.*x;                       % modulate message with carrier

r = 20*r1+r2;

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=50; fbe=[0 0.04 0.16 1];      % LPF parameters 
damps=[1 1 0 0 ]; 
b=firpm(fl,fbe,damps);        % create LPF impulse response

x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented 
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
y2=filter(fliplr(p)/(pow(p)*M),1,2*x2); %y2 has not been through the LPF
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % plot soft decisions
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
% decode decision device output to text string
reconstructed_message=pam2letters(mprime)

% This analysis is based on having another out-of-band user at fc=30 with
% sufficient power to cause errors if the LPF is not working correctly. By
% examining the frequency spectrum of x2 before the LPF is applied, we can
% see that these out of band interference peaks are about 4 frequency units
% wide, and centered at +/- 10 frequency units. This means that the LPF 
% needs to filter out frequencies more than 8 frequency units away from 0.
% Since this receiver works over +/- 50 frequency units, and 8/50 = 0.16,
% then designing a LPF to fully cut off any signals above 0.16 should
% entirely eliminate this out of band signal. However, since the
% correlation process seems to filter out any out of band peaks with
% sufficiently low power, so a more lenient LPF is technically possible
% depending on the power of the out of band interference. As far as the
% minimum frequency for the LPF, it must not be low enough to interfere
% with the in band signal. In this case, since this signal has a bandwidth
% of 4 (extending to +/- 2 frequency units away from the origin) the LPF
% cannot cut off much below 2/50 = 0.04.

%% Exercise 9.5. Using the same specifications (fbe=[0 0.1 0.2 1]; 
% damps= [1 1 0 0 ];), how short can you make the LPF? Explain.

%TRANSMITTER 1
% encode text string as T-spaced 4-PAM sequence
str='keldifnt97*jk5p[;ldr21nbaztqfdvclptnvtq3d$njj#5';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=30;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r1=c.*x;                       % modulate message with carrier


%TRANSMITTER 2
% encode text string as T-spaced 4-PAM sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=20;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r2=c.*x;                       % modulate message with carrier

r = 20*r1+r2;

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=4; fbe=[0 0.1 0.2 1];      % LPF parameters 
damps=[1 1 0 0 ]; 
b=firpm(fl,fbe,damps);        % create LPF impulse response

x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented 
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
y2=filter(fliplr(p)/(pow(p)*M),1,2*x2); %y2 has not been through the LPF
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % plot soft decisions
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])'; % quantize alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
% decode decision device output to text string
reconstructed_message=pam2letters(mprime)

% I found that the message could be recovered successfully with an LPF with
% length as short as 4. My explaination for this is that the length of the
% filter is also the order of the filter, and so a shorter filter has fewer
% terms to adjust the filter's shape in a way which approximates an ideal
% filter. Thus a shorter and as a result lower order filter will introduce
% more distortion into the signals it filters, and if this distortion is
% great enough, it will prevent the received message from being correctly
% reconstructed.

%% Better QPSK attempt


%TRANSMITTER
% encode text string as T-spaced 4-PAM sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2qpsk(str); N=length(m); % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled
% T/M-spaced sequence of scaled T-spaced pulses (T=1)
M=100;                        % oversampling factor
mup=zeros(1,N*M);             % Hamming pulse filter with 
mup(1:M:N*M)=m;               % T/M-spaced impulse response
p=hamming(M);                 % blip pulse of width M
x=filter(p,1,mup);            % convolve pulse shape with data
figure(1), plotspec(x,1/M)    % baseband AM modulation
t=1/M:1/M:length(x)/M;        % T/M-spaced time vector
fc=20;                        % carrier frequency
c=cos(2*pi*fc*t);             % carrier
r=c.*x;                       % modulate message with carrier

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);             % synchronized cosine for mixing
x2=r.*c2;                      % demod received signal
fl=50; fbe=[0 0.1 0.2 1];      % LPF parameters 
damps=[1 1 0 0 ]; 
b=firpm(fl,fbe,damps);        % create LPF impulse response

x3=2*filter(b,1,x2);           % LPF and scale signal
% extract upsampled pulses using correlation implemented 
% as a convolving filter; filter with pulse and normalize
y=filter(fliplr(p)/(pow(p)*M),1,x3);
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:N*M);           % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % plot soft decisions
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-1-j,-1+j,1-j,1+j])'; % quantize alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symbol error
% decode decision device output to text string
reconstructed_message=qpsk2letters(mprime)

