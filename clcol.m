function [col] = clcol(inp)

switch inp
    case 'b'
        [x,y,z] = rgbconv('4682B4');
    case 'o' 
        [x,y,z] = rgbconv('FF8C00');
    case 'g'
        [x,y,z] = rgbconv('2E8B57');
    case 'p'
        [x,y,z] = rgbconv('b22222');
end

col = [x y z];