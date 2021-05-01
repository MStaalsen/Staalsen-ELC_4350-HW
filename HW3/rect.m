function s = rect(x,width);
for i = 1:length(x)
if abs(x(i)) > abs(width)
    s(i) = 0;
else
    s(i) = 1;
end
end
    
    