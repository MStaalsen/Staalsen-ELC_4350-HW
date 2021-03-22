% f = pam2letters(seq)
% reconstruct string complex numbers of +/- 1 +/- j letters
function f = qpsk2letters(seq_complex)

S = length(seq_complex);
off = mod(S,4);

if off ~= 0
  sprintf('dropping last %i PAM symbols',off)
  seq_complex = seq_complex(1:S-off);
end

%convert complex numbers into +/- 1 and +/-3
seq= zeros(1,length(seq_complex));
for ndx = 1:length(seq_complex)
    val = seq_complex(ndx);
    if val == (-1-1i)
        seq(ndx) = -3;
    elseif val == (-1+1i)
        seq(ndx) = -1;
    elseif val == (1-1i)
        seq(ndx) = 1;
    else 
        seq(ndx) = 3;
    end    
end

N = length(seq_complex)/4;
f=[];
for k = 0:N-1
  f(k+1) = base2dec(char((seq(4*k+1:4*k+4)+99)/2),4);
end

f = char(f);

