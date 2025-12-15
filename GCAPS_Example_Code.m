clc;
clear all;

[~,~,dat] = xlsread('Sample Scattering Data for G-CAPS.xlsx', 'Sheet1');

q= cell2mat(dat(12:47,1)); % q-range for fitting (G-CAPS works only on a certain q-range, find by trial and error and try to keep it consistent across all systems you're comparing
I = cell2mat(dat(12:47,2)); % corresponding I

% Porod Analysis to get the background since G-CAPS is applied to bkgd-subtracted I(q)
q4 = q.^4;
Iq4 = I.*q4;
Porod = polyfit(q4,Iq4,1);
Porod_fit = Porod(1)*q4 + Porod(2);
bkgd = Porod(1);
I = I - bkgd;

% 4-term G-CAPS (I0, sigma, kappa4, kappa6) on the background subtracted sector-averaged intensity
fit_params0 = [6 150 0.01 0.0001];
GCAPS = @(fit_params, q)fit_params(1).*exp(-(q.*fit_params(2)).^2/2).*(1+(fit_params(3).*q.^4.*fit_params(2).^4)/(factorial(4).*fit_params(2).^4)+(fit_params(4).*q.^6.*fit_params(2).^6)/(factorial(6).*fit_params(2).^6));
[fit_params, R, J, CovB, MSE] = nlinfit(q, I, GCAPS, fit_params0);
fit_params_err = sqrt(diag(CovB));
GCAPS_fit = fit_params(1).*exp(-(q.*fit_params(2)).^2/2).*(1+(fit_params(3).*q.^4.*fit_params(2).^4)/(factorial(4).*fit_params(2).^4)+(fit_params(4).*q.^6.*fit_params(2).^6)/(factorial(6).*fit_params(2).^6));
I0 = [fit_params(1) fit_params_err(1)];      % pre-factor multiplying the form-factor, will depend on conc, SLD etc.
sigma = [fit_params(2) fit_params_err(2)];   % width of the segment-density distribution of the polymer
Rg = [(sqrt(3/2)*fit_params(2)) (sqrt(3/2)*fit_params_err(2))]; %Rg=sqrt(3/2)*sigma, Rg is the radius of gyration of the polymer
kappa4 = [fit_params(3) fit_params_err(3)];  % non-normalized kappa4
kappa6 = [fit_params(4) fit_params_err(4)];  % non-normalized kappa6
kappa4_norm = [kappa4(1)/sigma(1)^4 (((1/sigma(1)^4)*kappa4(2)) - ((kappa4(1)/sigma(1)^5)*sigma(2)))];                % normalized kappa4 (kappa4/sigma^4) with error propagated 
kappa6_norm = [kappa6(1)/sigma(1)^6 1/factorial(6)*(((1/sigma(1)^6)*kappa6(2)) - ((kappa6(1)/sigma(1)^7)*sigma(2)))]; % normalized kappa6 (kappa6/sigma^6) with error propagated 
I = I + bkgd;                 % adding back the bkgd
GCAPS_fit = GCAPS_fit + bkgd; % adding back the bkgd to the fit 

% Plotting I(q) vs q - Expt Data with G-CAPS fit
figure(1)
ax = gca();
loglog(q, I, 'bo', 'MarkerSize', 5)
hold on
loglog(q, GCAPS_fit, 'k-', 'LineWidth', 2)
ylim([0 10.0])
xlim([0.004 0.05])
xlabel('q_{avg}','Interpreter','tex','FontSize',12)
ylabel('I(q) (1/cm)','Interpreter','tex','FontSize',12)
legend({'1D sector-averaged intensity (del phi: 10 deg)','G-CAPS fit'})
grid(ax,'on')
set(gca,'Box','on');


% Constructing real-space 1-D segment density distribution from extracted moments
rx = linspace(-1500, 1500, 1000);
realspace_PDF = (1/(sqrt(2)*sigma(1))).*exp(-rx.^2/(2*sigma(1)^2)).*(1 + ((1/4).*kappa4(1)*hermiteH(4,rx/(sqrt(2)*sigma(1)))/(factorial(4)*sigma(1)^4)));

% Plotting simulated real-space 1-D segment density profile of the polymer 
figure(2)
ax = gca();
plot(rx, realspace_PDF, 'ro', 'MarkerSize', 5)
xlabel('r_{x}')
ylabel('Simulated 1-D real-space PDF')

% Displaying the fit parameters from G-CAPS
disp('Rg')
disp(Rg)
disp('kappa4 not normalized  by sigma')
disp(kappa4)
disp('kappa4 normalized by sigma')
disp(kappa4_norm)