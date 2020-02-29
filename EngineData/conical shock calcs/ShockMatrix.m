% Calculates conditions after shock for a variety of flight conditions
% Saves as file

hangle = 10.3885; % cone half angle (deg) % 10.3885 for third stage

mat = [];
for i = 7:.5:8.5
    for j = 0:.5:10
        delete('input')
        delete('output')
        
        dlmwrite('input',[i hangle j], 'delimiter',' ');
        
        system('Shock.exe')
        pause(5)
        
        temp = dlmread('output');
        
        mat(end+1,:) = [i j temp];
    end
end

dlmwrite('ShockMat',mat, 'delimiter',' ')