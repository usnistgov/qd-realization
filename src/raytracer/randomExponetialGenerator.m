function x = randomExponetialGenerator(lambda)

x = -(1./lambda).*log(1-(rand(length(lambda),1)));

end