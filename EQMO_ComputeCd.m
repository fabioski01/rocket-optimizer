function Cd = ComputeCd(M)
if M<=0.6
    Cd = 0.15;
elseif M>0.6 && M<=1.1
    Cd = -4.32*M^3+11.016*M^2-8.5536*M+2.24952;
elseif M>1.1 && M<=1.3
    Cd = -M^2+2.2*M-0.79;
elseif M>1.3 && M<=5
    Cd = 0.16769+0.17636/sqrt(M^2-1);
else
    Cd = 0;
end
end