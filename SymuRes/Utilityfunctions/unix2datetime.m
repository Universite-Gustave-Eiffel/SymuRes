function tdate = unix2datetime(tunix)
% tdate = unix2datetime(tunix)
% Convert unix timestamp to UTC time and date
% https://currentmillis.com/
%
% INPUTS
%---- tunix : integer, unix timestamp [s]
%
% OUTPUTS
%---- tdate : datetime, DD-month-YYYY HH:MM:SS

tref = datetime(1970,1,1,0,0,0); % reference date

tdate = tref + seconds(tunix);

end