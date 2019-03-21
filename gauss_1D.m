% gauss_1D for calculating 1-D gaussians
% call as Y = gauss_1D(parameters, X)
% parameters = [Amp center sigma background]

function gauss = gauss_1D(parameters, domain)

Amp = parameters(1);
center = parameters(2);
sigma = parameters(3);
Background = parameters(4);

gauss = (Amp*exp((-(domain - center).^2)/(2*sigma^2)) + Background);