%% Main Post-Processing Script for the NACA0012 Experimental Analysis %%
% A: Alex Krochak
clear all
close all
clc
%% 0. Inputs
% Folder names with HWA data
hwa_calib = '../hotwire/calibration/Calibration_%03d';
hwa_data = '../hotwire/data/';
piv_data = '../piv/data/'; % directory of piv data
piv_img = '../piv/images/'; % directory of piv images
fntSz = 15;
lblSz = 15;
Ny = 21; % number of y-measurements
Na = 3; % number of alpha measurements

wind = 4000; % welch method window size
%% 1. Hotwire Anemometry Analysis
kings = calibration(hwa_calib); % run the calibration routine to get kings coefficients
% Process 1 file to get the size

% Get the high-frequency measurement of the freestream
[ts,vs] = processHWA('../hotwire/CorrelationTest');
us = polyval(kings, vs);
mean1 = mean(us);
std1 = std(us);
us((us > mean1 + 3*std1) | (us < mean1 - 3*std1)) = nan;
us_bar = mean(us, "omitnan");
us_p = us - us_bar; 
sigma = std(us_p, "omitnan");
us_p(isnan(us_p))=0;

% assign to main data structure


fileName = strcat(hwa_data,'Measurement_-40_+00');
[t, ~] = processHWA(fileName);

Nt = length(t);
% Range of alpha / y
yvec = [-40:4:40];
alphavec = [0 5 15];

% Initialize the data structure
hwa.u = zeros(Na,Ny,Nt);
hwa.ubar = zeros(Na,Ny); % average velocity
hwa.rms = zeros(Na,Ny);
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
        std1 = std(vel);
        vel((vel > mean1 + 3*std1) | (vel < mean1 - 3*std1)) = nan;
        ubar = mean(vel, "omitnan");
        up = vel - ubar; 
        up_rms = rms(up, "omitnan");
        % assign to main data structure
        hwa.ubar(i,j) = ubar;
        hwa.rms(i,j) = up_rms;
    end
end

% Analyze the correlation from the correlation test measurement

ucorr = xcorr(us_p);
ucorr = ucorr(100000:end);
ucorr = ucorr/ucorr(1);
plot(ucorr)
val = find(ucorr<0.1,1);
T = ts(val);
facq = 1/(2*T); % acquisition frequency in Hz

% Compute the pre-multiplied energy spectrum
% L = length(us_p);
% uk = fft(us_p)/L;
% uk = uk(1:L/2);
% fs = 10000;  % sampling frequnecy
% frequencies = (1 : L/2) * fs / L; % Frequency vector
% % frequencies = fftshift(frequencies);
% uk_pre = uk .* frequencies';
% spectrum = abs(uk_pre);

L = length(us_p);
uk = pwelch(us_p,wind);
L = length(uk);
fs = 10000;  % sampling frequnecy
frequencies = (1 : L) * fs / L; % Frequency vector
% frequencies = fftshift(frequencies);
uk_pre = uk .* frequencies';
spectrum = abs(uk_pre);



close all
%% 3. Analyze the PIV data 
piv_inst_files = {'AoA_0__SP(32x32_50ov)./B00001.dat', 'AoA_5__SP(32x32_50ov)./B00001.dat', 'AoA_15__SP(32x32_50ov)./B00001.dat'};
piv_figure = {'AoA_0_32x32_50ov','AoA_5_32x32_50ov','AoA_15_32x32_50ov'};
for i = 1:length(piv_inst_files)
    dat_file = piv_inst_files{i};
    mat = importdata(strcat(piv_data,dat_file));
    mat = mat.data;
    Xog = mat(:,1); Yog = mat(:,2);
    mat(mat(:,5) == 0,1:4) = nan; % remove invalid points
%     quiver(mat(:,1),mat(:,2),mat(:,3),mat(:,4))
    X = mat(:,1); Y = mat(:,2); U = mat(:,3); V = mat(:,4);


    nanIndices = isnan(X) | isnan(Y) | isnan(U);

    % Define the grid on which to interpolate
    xi = linspace(min(X), max(X), 100); % X grid points
    yi = linspace(min(Y), max(Y), 100); % Y grid points
    [xi, yi] = meshgrid(xi, yi);
    
    % Perform the interpolation
    Ui = griddata(X(~nanIndices), Y(~nanIndices), U(~nanIndices), xi, yi);
    
    % Create the contour plot
%     contourf(xi, yi, Ui, 'LineWidth', 0.5, 'LineColor', 'k');hold on;

    sp = 3; % quiver spacing
    contLvl = -5:1:13;
    tickVals = linspace(contLvl(1), contLvl(end), length(contLvl));
    %% Define custom colormap
    numColors = length(contLvl);
    % Get the "Hot" colormap
    hotMap = hot(numColors);
    
    % Get the "Cold" colormap
    coldMap = 1 - hotMap;
    
    % Create a custom colormap by interpolating between "Hot" and "Cold"
    hotNcold = [hotMap; coldMap];
    data = Ui; newdata = data; 
    for tv = 2:length(contLvl)
        % find all data between the new tick ranges, one block at a time
        ind = data>contLvl(tv-1) & data<=contLvl(tv);
        % Change the corresponding values in the copied data to scale to the
        % block edges
        newdata(ind) = rescale(data(ind),tickVals(tv-1),tickVals(tv));
    end
        
    set(gcf,'Position',[100 100 1000 600])
    contourf(xi,yi,newdata,tickVals); hold on
    xlim([-20, 160]);ylim([-160, 30]);
    scatter(Xog(nanIndices), Yog(nanIndices), 'Marker', 'square', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none','SizeData', 130);
    quiver(X(1:sp:end),Y(1:sp:end),U(1:sp:end),V(1:sp:end),'k');
    xlim([-20, 160]);ylim([-60, 30]);
    hold off
    C=hotNcold;
    tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
    colormap(flipud(C))
    colorbar('Ticks',tickVals,...
        'TickLabels',[string(contLvl(1:end)) ""],'TickLabelInterpreter','latex');
    % colorbar('Ticks',tickVals,...
    %     'TickLabels',flip({'','5(a)', '4(b)' ,'3(c)', '2(d)', '1(e)' ,'0.5(f)' ,'0(g)' ,'-0.5(h)' ,'-1(i)', '-2(j)', '-3(k)'}) ...
    %     ,'TickLabelInterpreter','latex');
    
    title('$U_{\mathrm{x}}$ [m/s]','interpreter','latex','FontSize',lblSz);
    xlabel('$x$ [mm]','interpreter','latex','FontSize',lblSz);
    ylabel('$y$ [mm]','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
    exportgraphics(gcf,["figures/"+piv_figure{i}+".pdf"], 'Resolution', 300)

end
%% 3. Plot the results for HWA

% Plot the time scale
figure(1)
plot(ts,ucorr); hold on
yline(0);
scatter(T,ucorr(val)); hold off
xlim([0,0.05]); ylim([-0.1,1]);
title('Normalized Auto-correlation Coefficient','Interpreter','latex','FontSize',fntSz);
ylabel('$\rho_{xx}$','Interpreter','latex','FontSize',fntSz);
xlabel('$\tau$ [s]','Interpreter','latex','FontSize',fntSz);
set(gca,'ticklabelinterpreter','latex')
set(gca,'FontSize',fntSz)

% Plot the pre-multiplied spectrum
figure(2)
semilogx(frequencies,10*log10(spectrum))
title('$f|\hat{u}^{\prime}|$','Interpreter','latex','FontSize',fntSz);
xlabel('$f$ [Hz]','Interpreter','latex','FontSize',fntSz);
set(gca,'ticklabelinterpreter','latex')
set(gca,'FontSize',fntSz)

figure(10)
plot(yvec,hwa.ubar(1,:)); hold on
plot(yvec,hwa.ubar(2,:));
plot(yvec,hwa.ubar(3,:)); hold off
ylabel('$U$ [m/s]','Interpreter','latex','FontSize',fntSz);
xlabel('$y$ [mm]','Interpreter','latex','FontSize',fntSz);
title('Mean Velocity - Wake','Interpreter','latex','FontSize',fntSz);
legend('$\alpha = 0$ deg', '$\alpha = 5$ deg', '$\alpha = 15$ deg','interpreter','latex','Location','southwest')
set(gca,'ticklabelinterpreter','latex')
set(gca,'FontSize',fntSz)

figure(11)
plot(yvec,hwa.rms(1,:)); hold on
plot(yvec,hwa.rms(2,:));
plot(yvec,hwa.rms(3,:)); hold off
ylabel('$U_{\mathrm{rms}}$ [m/s]','Interpreter','latex','FontSize',fntSz);
xlabel('$y$ [mm]','Interpreter','latex','FontSize',fntSz);
title('Velocity Fluctutations - Wake','Interpreter','latex','FontSize',fntSz);
legend('$\alpha = 0$ deg', '$\alpha = 5$ deg', '$\alpha = 15$ deg','interpreter','latex','Location','northwest')
set(gca,'ticklabelinterpreter','latex')
set(gca,'FontSize',fntSz)

