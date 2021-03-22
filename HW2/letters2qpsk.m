% f = letters2qpsk(str)
% encode a string of ASCII text into +/-1, +/-3

function f_complex = letters2qpsk(str);           % call as Matlab function
N=length(str);                           % length of string
f=zeros(1,4*N);                          % store 4-PAM coding here
for k=0:N-1                              % change to "base 4"
  f(4*k+1:4*k+4)=2*(dec2base(double(str(k+1)),4,4))-99;
end

f_complex = zeros(1,length(f));
for ndx = 1:length(f)
    val = f(ndx);
    if val == -3
        f_complex(ndx) = -1-1i;
    elseif val == -1
        f_complex(ndx) = -1+1i;
    elseif val == 1
        f_complex(ndx) = 1-1i;
    else 
        f_complex(ndx) = 1+1i;
    end    
end