function [kings] = calibration(fileName)

%%  Config
urange = 0:2:20; 
fntSz = 15;
lblSz = 15;
% fileName = strcat(hwa_calib, 'Calibration_%03d'); 

for ii = 1:numel(urange)
    fileIn = sprintf(fileName,urange(ii));
    delimiter = ' ';
    startRow = 23;
    formatSpec = '%s';
    try
        fileID = fopen(fileIn,'r');
    catch
        fileID = fopen(fileIn{1},'r');
    end
    tmp = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tmp2 = strrep(tmp{1},',','.');
    tmp2 = str2mat(tmp2);
    vel = tmp2(:,10:17);
    u(ii) = mean(str2num(vel));
    
end

kings = polyfit(u,urange,4);

figure(1)
scatter(urange,u); hold on
plot(polyval(kings,1.1:0.001:1.82),1.1:0.001:1.82)
legend('calibration data', 'calibration polynomial','Interpreter','latex','Location','southeast')
title('Calibration Curve');
xlabel('$U$ [m/s]','Interpreter','latex','FontSize',fntSz)
ylabel('$V$ [volt]','Interpreter','latex','FontSize',fntSz)
% set(gca,)
grid on
hold off
end