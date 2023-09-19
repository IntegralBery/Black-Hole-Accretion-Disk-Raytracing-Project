format short g
axis equal

pixels = matfile('output.mat').pixels;

normalize_percentile = 95;
[v_res, h_res, ~] = size(pixels);
rgb_pixels = zeros(v_res, h_res, 3);
rgb_array = [];
for i = 1:h_res
    fprintf('\n Now colouring column: %d', i)
    for s = 1:v_res
        redshift = pixels(s, i, 1);
        T = pixels(s, i, 2);
        [r, g, b] = get_rgb(redshift, T);
        rgb_pixels(s, i, :) = [r, g, b];
        
        if r ~= 0.05 && r ~= 0
            rgb_array(end+1:end+3) = [r, g, b];
        end
    end
end
normalize_value = prctile(rgb_array, normalize_percentile);
for i = 1:h_res
    fprintf('\n Now colouring column: %d', i)
    for s = 1:v_res
        rgb = rgb_pixels(s, i, :);
        if rgb(1) ~= 0.05 && rgb(2) ~= 0.05 && rgb(3) ~= 0.05
            rgb_pixels(s, i, :) = rgb/normalize_value;
        end
    end
end

rendered = imshow(rgb_pixels)
saveas(rendered,'img_name.png'); %Remember to change the name of your saved image!



function [R, G, B] = get_rgb(redshift, T)
    if redshift == 0
        if T == 0
            R = 0;
            G = 0;
            B = 0;
        else
            R = 0.05;
            G = 0.05;
            B = 0.05;
        end
    else
        X = integrate(@x_gaussian, redshift, T);
        Y = integrate(@y_gaussian, redshift, T);
        Z = integrate(@z_gaussian, redshift, T);
        
        RGB =  [3.24096994 -1.53738318 -0.49861076;
                      -0.96924364 1.8759675 0.04155506;
                      0.05563008 -0.20397696 1.05697151]*[X; Y; Z];
        RGB = arrayfun(@rgb_transform, RGB);
        R = max(0, RGB(1));
        G = max(0, RGB(2));
        B = max(0, RGB(3));
    end
end

function u = rgb_transform(u0)
    if u0 <= 0.0031308
        u = 323*u0/25;
    else
        u = (211*u0^(5/12) - 11)/200;
    end
end

function y = integrate(f1, redshift, T)
    % integrates over the wavelength spectrum to produce r,g,b values of
    % pixel
    sum = 0;
    intervals = 100;
    lower = 300e-09;
    upper = 800e-09;
    h = (upper-lower)/intervals;
    for i=0:intervals
        x = lower+h*i;
        fi = f1(x*1e+10)*plancks_law(x/redshift, T);
        if (i == 0) || (i == intervals)
            sum = sum + fi;
        elseif mod(i,2)==0
            sum = sum + 4*fi;
        else
            sum = sum + 2*fi;
        end
    end
    y = sum*h/3;
end

function value = gaussian(x, alpha, mu, sigma1, sigma2)
    if x < mu
        sigma = sigma1;
    else 
        sigma = sigma2;
    end
    value = alpha*exp((x-mu)^2 /(-2*sigma^2));
end

function x_bar = x_gaussian(lambda) %These are in units of angstroms, as I realised later
    x_bar = gaussian(lambda,  1.056, 5998, 379, 310) + gaussian(lambda,  0.362, 4420, 160, 267) + gaussian(lambda, -0.065, 5011, 204, 262);
end

function y_bar = y_gaussian(lambda)
    y_bar = gaussian(lambda,  0.821, 5688, 469, 405) + gaussian(lambda,  0.286, 5309, 163, 311);
end

function z_bar = z_gaussian(lambda)
    z_bar = gaussian(lambda,  1.217, 4370, 118, 360) + gaussian(lambda,  0.681, 4590, 260, 138);
end

function M = plancks_law(lambda, T)
    M = (3.741844e-16)./((lambda.^5).*(exp((1.438833e-2)./(lambda * T))-1));
end