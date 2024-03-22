function x = mynewt(f, fp, x)

for i = 1:100
    
    chi = f(x)/fp(x);
    
    x = x - chi;
    
    if norm(chi) < 10^-10
        
        break;
        
    end
    
end

end
