function x = rndLaplace(sigma, m, n)
u = rand(m,n) - 0.5;
b = sigma / sqrt(2);
x = -b * sign(u) .* log(1 - 2*abs(u));
end