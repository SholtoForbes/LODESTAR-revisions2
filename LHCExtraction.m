4% Extract LHC stuff 
cd ./ArchivedResults/LHC-Updated
clear all
%Uncertainties (%)
unc.Isp1 = 1.3;
unc.CL12_subsonic = 17;
unc.CL12_transonic = 28.7;
unc.CL12_supersonic = 12;
unc.CD12_subsonic = 33;
unc.CD12_transonic = 21;
unc.CD12_supersonic = 11;
unc.Cm12_subsonic = 23;
unc.Cm12_transonic = 67.1;
unc.Cm12_supersonic = 22;
unc.Isp2 = 25;
unc.Isp3 = 1.3;

n = 1

cd ./LHC1
load output.mat

for i = 1:length(output)
    if any([isequal(i,3) isequal(i,6) isequal(i,7)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC2
load output.mat

for i = 1:length(output)
if any([isequal(i,2) isequal(i,3) isequal(i,6)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
end
end

cd ..\

cd ./LHC3
load output.mat

for i = 1:length(output)
    if any([isequal(i,2) isequal(i,5)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC4
load output.mat

for i = 1:length(output)
    if any([isequal(i,3)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC5
load output.mat

for i = 1:length(output)
    if any([isequal(i,2)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC6
load output.mat

for i = 1:length(output)
 if any([isequal(i,5) isequal(i,10)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
 end
end

cd ..\

cd ./LHC7
load output.mat

for i = 1:length(output)
    if any([isequal(i,9) isequal(i,10)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC8
load output.mat

for i = 1:length(output)
    if any([isequal(i,3)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC9
load output.mat

for i = 1:length(output)
    if any([isequal(i,7)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC10
load output.mat

for i = 1:length(output)
    if any([isequal(i,6) isequal(i,7) isequal(i,8) isequal(i,10)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
    end
end

cd ..\

cd ./LHC11
load output.mat

for i = 1:length(output)
 if any([isequal(i,1)])
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
 end
end

cd ..\

cd ./LHC12
load output.mat

for i = 1:length(output)

        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1

end

cd ..\

cd ./LHC13
load output.mat

for i = 1:length(output)
 if any([isequal(i,2) isequal(i,4) isequal(i,5) isequal(i,6) isequal(i,7) isequal(i,8) isequal(i,9) isequal(i,10)]) % this makes 100 runs that are analysed
    % Do nothing if a bad run is selected
    else
        objVals(n) = abs(output{i}.result.objective);

        Isp1s(n) = output{i}.result.setup.auxdata.Isp1mod;
        CLs(n) = output{i}.result.setup.auxdata.CL12_subsonicmod;
        CDs(n) = output{i}.result.setup.auxdata.CD12_subsonicmod;
        CMs(n) = output{i}.result.setup.auxdata.Cm12_subsonicmod;
        Isp2s(n) = output{i}.result.setup.auxdata.Isp2mod;
        Isp3s(n) = output{i}.result.setup.auxdata.Isp3mod; 
       
        n = n+1
 end
end

cd ..\



histogram(objVals,7)

sortedVals = sort(objVals)
pd = pdf('Normal',sortedVals, mean(objVals), std(objVals))

% params = gamfit(sortedVals);
% pd = gampdf(sortedVals, params(1), params(2))

plot([sortedVals(1):sortedVals(end)],spline(sortedVals,pd,[sortedVals(1):sortedVals(end)]),'LineWidth',2)
ylabel('Probability Density')
xlabel('Payload-to-Orbit')

% https://au.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2))
CI = CIFcn(objVals,97.5)

% std(objVals)*2

% SEM = std(objVals)/sqrt(length(objVals));               % Standard Error
% ts = tinv([0.025  0.975],length(objVals)-1);      % T-Score
% CI = mean(objVals) + ts*SEM;  

nBoot = 2000; %number of bootstraps
[bci,bmeans] = bootci(nBoot,{@mean,objVals},'alpha',.1,'type','per'); %90 confidence interval
