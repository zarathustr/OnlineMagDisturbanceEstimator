function K = B2K(B)
        z = [B(2, 3) - B(3, 2); 
             B(3, 1) - B(1, 3); 
             B(1, 2) - B(2, 1)];
        K = [trace(B), z';
             z, B + B' - trace(B) * eye(3)];
end