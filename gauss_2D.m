% gauss_2D for 2-D gaussians 
% call as Y = gauss_2D(parameters, domain) with 
% parameters = [Amp center_x center_y sigma background]

function gauss = gauss_2D(parameters, domain)

Amp = parameters(1);
center_x = parameters(2);
center_y = parameters(3);
sigma_x = parameters(4);
sigma_y = parameters(5);
Background = parameters(6);

if size(domain, 2) == 1
    [D_x, D_y] = meshgrid(domain);
elseif size(domain, 2) == 2
    [D_x, D_y] = meshgrid(domain(:,1), domain(:,2));
end

gauss = Amp*exp(-((((D_x - center_x).^2)./(2*sigma_x.^2)) + (((D_y - center_y).^2)./(2*sigma_y.^2)))) + Background;