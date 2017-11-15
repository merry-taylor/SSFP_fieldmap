function I = EllipticalModel2D(I1, I2, I3, I4)

% Iterate through each pixel and calculate M directly; then compare it to
% the maximum magnitude of all four input images.  If it is greater,
% replace it with the complex sum value.
[Mr,Nc]=size(I1);
M = zeros(Mr, Nc);
maximum = max(abs(I1), abs(I2));
maximum = max(abs(I3), maximum);
maximum = max(abs(I4), maximum);
CS = (I1 + I2 + I3 + I4) / 4;
for k = 1:Mr
    for n = 1:Nc
        M(k,n) = ((real(I1(k,n))*imag(I3(k,n)) - real(I3(k,n))*imag(I1(k,n)))*...
            (I2(k,n) - I4(k,n)) - (real(I2(k,n))*imag(I4(k,n)) - real(I4(k,n))*...
            imag(I2(k,n)))*(I1(k,n) - I3(k,n))) / ((real(I1(k,n)) - real(I3(k,n)))*...
            (imag(I2(k,n)) - imag(I4(k,n))) + (real(I2(k,n)) - real(I4(k,n)))*...
            (imag(I3(k,n)) - imag(I1(k,n)))); % Equation (13)
        if (abs(M(k,n)) > maximum(k,n)) || isnan(M(k,n))
            M(k,n) = CS(k,n); % This removes the really big singularities; without this the image is mostly black.
        end
    end
end

% Calculate the weight w for each pixel.
w1 = zeros(Mr, Nc); %bottle is 512 x 256
w2 = zeros(Mr, Nc);
for k = 1:Mr %bottle is 512
    for n = 1:Nc
        numerator1 = 0;
        denominator1 = 0;
        numerator2 = 0;
        denominator2 = 0;
        for x = -2:2
            a = k + x;
            for y = -2:2
                b = n + y;
                if (a < 1) || (b < 1) || (a > Mr) || (b > Nc)
                    
                else
                    numerator1 = numerator1 + conj(I3(a,b) - M(a,b)) * (I3(a,b) - I1(a,b)) + ...
                        conj(I3(a,b) - I1(a,b)) * (I3(a,b) - M(a,b));
                    denominator1 = denominator1 + conj(I1(a,b) - I3(a,b)) * (I1(a,b) - I3(a,b));
                    numerator2 = numerator2 + conj(I4(a,b) - M(a,b)) * (I4(a,b) - I2(a,b)) + ...
                        conj(I4(a,b) - I2(a,b)) * (I4(a,b) - M(a,b));
                    denominator2 = denominator2 + conj(I2(a,b) - I4(a,b)) * (I2(a,b) - I4(a,b));
                end
            end
        end
        w1(k,n) = numerator1 / (2 * denominator1); % Equation (18) - first pair
        w2(k,n) = numerator2 / (2 * denominator2); % Equation (18) - second pair
    end
end

% Calculate the average weighted sum of image pairs.
I = (I1 .* w1 + I3 .* (1 - w1) + I2 .* w2 + I4 .* (1 - w2)) / 2; % Equation (14) - averaged

end