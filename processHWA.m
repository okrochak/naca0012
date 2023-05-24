function [t, u] = processHWA(fileIn)

% fileIn = 'CorrelationTest';

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
t = tmp2(:,1:8);
t = str2num(t);
u= tmp2(:,10:17);
u = str2num(u);

%% Use this to calculate proper sampling frequency
% ubar = mean(u);
% up = u-ubar;
% ucorr = xcorr(up);
% ucorr = ucorr(100000:end);
% ucorr = ucorr/ucorr(1);
% plot(ucorr)
% val = find(ucorr<0.1);
% val = val(1);
% T = t(val);
% facq = 1/(2*T);
end
