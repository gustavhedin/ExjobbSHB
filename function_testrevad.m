function output = function_testrevad(x,y)
    output = adr_mul(adr_add(x,y),adr_add(y,ADRev(1))); % times and plus
    %output = adr_div(x,y); % division
    %output = adr_sub(x,y); % sub
end
