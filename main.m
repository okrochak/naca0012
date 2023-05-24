%% Main Post-Processing Script for the NACA0012 Experimental Analysis %%
% A: Alex Krochak
%% 0. Inputs
% Folder names with HWA data
hwa_calib = 'calibration/Calibration_%03d';
hwa_data = 'data/';

fntSz = 15;
lblSz = 15;
Ny = 21; % number of y-measurements
Na = 3; % number of alpha measurements
%% 1. Hotwire Anemometry Analysis
kings = calibration(hwa_calib); % run the calibration routine to get kings coefficients
% Process 1 file to get the size
fileName = 'data/Measurement_-40_+00';
[t, u] = processHWA(fileName);
Nt = length(t);
% Range of alpha / y
yvec = [-40:4:40];
alphavec = [0 5 15];

% Initialize the data structure
hwa.u = zeros(Na,Ny,Nt);
hwa.ubar = zeros(Na,Ny); % average velocity
hwa.sigma = zeros(Na,Ny)
hwa.t = t;
hwa.y = zeros(Na,Ny);
hwa.alpha = zeros(Na,Ny);

for i = 1:Na
    for j = 1:Ny
        a = alphavec(i);
        y = yvec(j);
        % fix the naming lol
        if a > 10
            astr = strcat("+", num2str(a));
        elseif a < -10
            astr = num2str(a);
        elseif (0 < a) && (a < 10)
            astr = strcat("+0", num2str(a));
        elseif (-10 < a) && (a < 0)
            astr = strcat("-0", num2str(abs(a)));
        else
            astr = '+00';
        end
        if y > 10
            ystr = strcat("+", num2str(y));
        elseif y < -10
            ystr = num2str(y);
        elseif (0 < y) && (y < 10)
            ystr = strcat("+0", num2str(y));
        elseif (-10 < y) && (y < 0)
            ystr = strcat("-0", num2str(abs(y)));
        else
            ystr = '+00';
        end
        fileName = strcat(hwa_data,'Measurement_', ystr,'_',astr);
        [~,volt] = processHWA(fileName);
        vel = polyval(kings, volt);
        hwa.u(i,j,:) = vel;
        % filter the outliers
        mean1 = mean(vel);
        vel((vel > mean1 + 3*std1) | (vel < mean1 - 3*std1)) = nan;
        ubar = mean(vel, "omitnan");
        up = vel - ubar; 
        sigma = std(up, "omitnan");
        % assign to main data structure
        hwa.ubar(i,j) = ubar;
        hwa.sigma(i,j) = sigma;
    end
end
% 
% for i = 1:Na
%     for j = 1:Ny
%     plot(hwa.t,squeeze(hwa.u(i,j,:))); hold on;
%     end
% end

%% 3. Plot the results for HWA

figure(10)
plot(yvec,hwa.ubar(1,:)); hold on
plot(yvec,hwa.ubar(2,:));
plot(yvec,hwa.ubar(3,:));
ylabel('$U$ [m/s]','Interpreter','latex','FontSize',fntSz);
xlabel('$y$ [mm]','Interpreter','latex','FontSize',fntSz);
legend('$\alpha = 0$ deg', '$\alpha = 5$ deg', '$\alpha = 15$ deg','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca,'FontSize',fntSz)
