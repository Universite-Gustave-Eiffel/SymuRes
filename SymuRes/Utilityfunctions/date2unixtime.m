function tunix = date2unixtime(tdate)
% tunix = date2unixtime(tdate)
% Convert UTC time and date to unix timestamp
% https://currentmillis.com/
%
% INPUTS
%---- tdate : datetime, DD-month-YYYY HH:MM:SS
%
% OUTPUTS
%---- tunix : integer, unix timestamp [s]

tref = datetime(1970,1,1,0,0,0); % reference date

dt = caldiff([tref tdate],'time');
tunix = seconds(time(dt));

end