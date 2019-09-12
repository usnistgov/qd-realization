function x = rndExp(lambda,m,n)
x = -log(1-rand(m,n)) / lambda;
end