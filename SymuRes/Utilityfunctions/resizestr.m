function newstr = resizestr(str,N)
% newstr = resizestr(str,N)
% Resize a char or a string to a fixed length (number of characters): add
% spaces at the end if the original string is too short, or cut its end if
% too long.
%
% INPUTS
%---- str : char or string, string to resize
%---- N   : desired new size (number of characters)
%
% OUTPUTS
%---- newstr : char or string, resized string

if length(str) > N
    newstr = str(1:N);
else
    newstr = char(32*ones(1,N)); % string with N spaces (char 32 in ASCII)
    newstr(1:length(str)) = str;
end

end