function pg = friisPathGain(d,f)
% c/(4*pi) = 23873241.463784
pg = 20 * log10(23873241.463784 / (f*d));
end