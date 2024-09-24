% plots result for LAVA twin exp
%
% itimu Ap3 2012

clear all; close all;

% Get the current directory
current_dir = fileparts(mfilename('fullpath'));

experience = "13_exp.ini";

oneabovedir = '../twin_exp/Step_2/Input/';

paramFile = [oneabovedir+experience];

addpath(genpath('m_map'))

fid = fopen(paramFile);
disp(['Reading file ---> ',paramFile]);
%read first line with drifter info 
dum=fgets(fid); %get first line
strdrft=dum(11:end); %just get drifter numbers
drft = cell2mat(textscan(strdrft,'%d'));
%read the rest 
par = textscan(fid, '%s%s');
fclose(fid);

Ndrft = length(drft);
NT    = str2double(cell2mat(par{2}(3)));
Ncols = str2double(cell2mat(par{2}(1)));
Nrows = str2double(cell2mat(par{2}(2)));

velfile = ['../',char(par{2}(9))]; %set velocity input file
gridfile = ['../',char(par{2}(10))]; %set grid input file
drftpfx = ['../',char(par{2}(11))]; %set file prefix for drifters
outdir = ['../',char(par{2}(12))]; %set output dir

%set input files for drifters
for idrft = 1 : Ndrft
    if (idrft<10)
        drftfiles{idrft} = [drftpfx,'_0',num2str(idrft)];
    else
        drftfiles{idrft} = [drftpfx,'_',num2str(idrft)];
    end
end

bckvelfile = [outdir,'/velbck.dat']; %set file for background velocities
estvelfile = [outdir,'/velest.dat']; %set file for velocity estimates
tgeobkfile = [outdir,'/trjbck.dat']; %set file for background trajectories (georef.)
tgeoesfile = [outdir,'/trjest.dat']; %set file for trajectory estimates (georef.)

disp(' ');
disp(['Reading file ---> ',estvelfile]);
velest = load(estvelfile); %load velocity estimates

disp(' ');
disp(['Reading file ---> ',bckvelfile]);
velbck = load(bckvelfile); %load backward velocities

disp(' ');
disp(['Reading file ---> ',tgeoesfile]);
tgeoes = load(tgeoesfile); %load geo trajectory estimates

disp(' ');
disp(['Reading file ---> ',tgeobkfile]);
tgeobk = load(tgeobkfile); %load geo backward trajectories

disp(' ');
disp(['Reading file ---> ',gridfile]);
grid = load(gridfile); %load grid info

if( strcmp(velfile,'../setAllZeros') ) %only drifter mode
   u = zeros(Ncols,Nrows,NT);
   v = zeros(Ncols,Nrows,NT);					  
else
  fid = fopen(velfile,'r'); %load velocities
  disp(' ');
  disp(['Reading file ---> ',velfile]);
  veldata = textscan(fid,'%s%s%f%f');
  fcl = fclose(fid);

  u = veldata{:,3}(:); 
  u = reshape(u,[Ncols*Nrows NT]);
  u = reshape(u,[Ncols Nrows NT]);

  v = veldata{:,4}(:);
  v = reshape(v,[Ncols*Nrows NT]);
  v = reshape(v,[Ncols Nrows NT]);
end  
  
lon = reshape(grid(:,1),[Ncols Nrows]); 
lat = reshape(grid(:,2),[Ncols Nrows]);
msk = reshape(grid(:,3),[Ncols Nrows]);

Ue = zeros(Ncols,Nrows,NT);
Ve = zeros(Ncols,Nrows,NT);
Ub = zeros(Ncols,Nrows,NT);
Vb = zeros(Ncols,Nrows,NT);
trjex = zeros(NT,Ndrft);
trjey = zeros(NT,Ndrft);
trjbx = zeros(NT,Ndrft);
trjby = zeros(NT,Ndrft);

for it = 1 : NT 
    Ue(:,:,it) = reshape(velest((it-1)*Ncols*Nrows+1:it*Ncols*Nrows,1),[Ncols Nrows]);
    Ve(:,:,it) = reshape(velest((it-1)*Ncols*Nrows+1:it*Ncols*Nrows,2),[Ncols Nrows]);
    Ub(:,:,it) = reshape(velbck((it-1)*Ncols*Nrows+1:it*Ncols*Nrows,1),[Ncols Nrows]);
    Vb(:,:,it) = reshape(velbck((it-1)*Ncols*Nrows+1:it*Ncols*Nrows,2),[Ncols Nrows]);    
    %
    trjex  (it,:) = tgeoes(Ndrft*(it-1)+1:Ndrft*it,1); 
    trjey  (it,:) = tgeoes(Ndrft*(it-1)+1:Ndrft*it,2); 
    trjbx  (it,:) = tgeobk(Ndrft*(it-1)+1:Ndrft*it,1); 
    trjby  (it,:) = tgeobk(Ndrft*(it-1)+1:Ndrft*it,2); 
end
    
%swap cols and rows to get rows in the first dimension
lon = permute(lon,[2 1]);
lat = permute(lat,[2 1]);
msk = permute(msk,[2 1]);
Ue = permute(Ue,[2 1 3]);
Ve = permute(Ve,[2 1 3]);
Ub = permute(Ub,[2 1 3]);
Vb = permute(Vb,[2 1 3]);
u = permute(u,[2 1 3]);
v = permute(v,[2 1 3]);

trjbx(trjbx==0) = NaN;
trjby(trjby==0) = NaN;

trjex(trjex==0) = NaN;
trjey(trjey==0) = NaN;

%gets timestep when drifter enters radar domain
itEnter = zeros(Ndrft,1) ;

for idrft = 1 : Ndrft
    itEnter(idrft) = sum(isnan(trjex(:,idrft)))+1;
end

clear velest velbck tgeoes tgeobk

trjobsx = zeros(NT,Ndrft);
trjobsy = zeros(NT,Ndrft);

% Read in the real observed drifter trajectories
% from Input directory this time

for idrft = 1 : Ndrft
   fid = fopen(drftfiles{idrft},'r'); %load drifter #idfrt
   disp(' ');
   disp(['Reading file ---> ',drftfiles{idrft}]);
   trjdata = textscan(fid,'%s%s%f%f');
   fcl1 = fclose(fid);
%
   trjobsx(:,idrft) = trjdata{:,3}(:);
   trjobsy(:,idrft) = trjdata{:,4}(:);
% get times
   timestr{idrft} = [char(trjdata{:,1}(:)),repmat(' ',NT,1),char(trjdata{:,2}(:))];
   time{idrft} = datenum(timestr{idrft},0);
end

% Load the true velocity data (U_truth, V_truth)
truthFile = fullfile(current_dir, 'U_truth.mat');
load(truthFile, 'uRad', 'vRad'); % Ensure this file contains uRad and vRad

% Ajout du chemin vers m_map
addpath('/home/lienen/ismar-lava-5a07b79ebe08_netcdf_TEST/m_map');

% Initialize arrays to store E metrics
itStrAvg = 1; % Start averaging from this index
itEndAvgArray = itStrAvg:NT;
E_metricsArray = zeros(size(itEndAvgArray));

% Calculate the initial reference norm for normalization
initial_avgUtruth = mean(uRad(:, :, itStrAvg:itStrAvg), 3);
initial_avgVtruth = mean(vRad(:, :, itStrAvg:itStrAvg), 3);
% Logical indexing to ignore NaN values
validMask_initial = ~isnan(initial_avgUtruth) & ~isnan(initial_avgVtruth);

% Calculate the initial norm using the validMask
initial_norm = sum(sum((initial_avgUtruth(validMask_initial)).^2 + (initial_avgVtruth(validMask_initial)).^2));

% Loop through each time step to calculate E metrics
for idx = 1:length(itEndAvgArray)
    itEndAvg = itEndAvgArray(idx);
    
    % Calculate the average estimated velocity fields
    avgUe = mean(Ue(:, :, itStrAvg:itEndAvg), 3);
    avgVe = mean(Ve(:, :, itStrAvg:itEndAvg), 3);
    
    % Calculate the average true velocity fields
    avgUtruth = mean(uRad(:, :, itStrAvg:itEndAvg), 3);
    avgVtruth = mean(vRad(:, :, itStrAvg:itEndAvg), 3);
    

    % Ensure the velocity fields have the same dimensions
    assert(all(size(avgUe) == size(avgUtruth)), 'Velocity fields must have the same dimensions');
    assert(all(size(avgVe) == size(avgVtruth)), 'Velocity fields must have the same dimensions');
    
    % Logical indexing to ignore NaN values and apply mask
    validMask = msk & ~isnan(avgUe) & ~isnan(avgVe) & ~isnan(avgUtruth) & ~isnan(avgVtruth);
    
    % Calculate the E metrics
    numerator = sum(sum((avgUe(validMask) - avgUtruth(validMask)).^2 + (avgVe(validMask) - avgVtruth(validMask)).^2));
    denominator = initial_norm; % Use the initial norm for normalization
    
    % Ensure denominator is not zero to avoid division by zero
    if denominator == 0
        E_metricsArray(idx) = NaN; % Or handle as appropriate
    else
        E_metricsArray(idx) = numerator / denominator;
    end
end

% Plot the E metrics
figure(100);
plot(itEndAvgArray, E_metricsArray, 'LineWidth', 2);
xlabel('Time Steps');
ylabel('E metrics');
title('E metrics over Time');
grid on;

% Exporter la matrice E_metricsArray en CSV
csvwrite(fullfile(outdir,'E_metricsArray.csv'), E_metricsArray);
% Save the E metrics figure
E_metrics_file = fullfile(outdir, sprintf('E_metrics.png', experience));
saveas(gcf, E_metrics_file);

disp('E metrics calculated and plotted successfully.');

% Initialize arrays to store L metrics
L_metricsArray = zeros(1, NT); % Array to store L metrics

% Read the true drifter trajectories
truthTrajectories = cell(Ndrft, 1);
for idrft = 1:Ndrft
    fid = fopen(drftfiles{idrft}, 'r'); % load drifter #idfrt
    disp(' ');
    disp(['Reading file ---> ', drftfiles{idrft}]);
    trjdata = textscan(fid, '%s%s%f%f');
    fclose(fid);
    truthTrajectories{idrft} = [trjdata{3}, trjdata{4}]; % Store as [lat, lon] pairs
end

% Calculate L metrics
for t = 1:NT
    sumSqDist = 0;
    numDrifters = 0;
    
    for idrft = 1:Ndrft
        estTraj = [trjex(t, idrft), trjey(t, idrft)];
        if all(~isnan(estTraj))
            truthTraj = truthTrajectories{idrft}(t, :);
            if all(~isnan(truthTraj))
                % Calculate the distance using m_idist
                dist = m_idist(truthTraj(1), truthTraj(2), estTraj(1), estTraj(2));
                sumSqDist = sumSqDist + dist^2;
                numDrifters = numDrifters + 1;
            end
        end
    end
    
    if numDrifters > 0
        L_metricsArray(t) = sqrt(sumSqDist / numDrifters);
    else
        L_metricsArray(t) = NaN;
    end
end

% Integrate L metrics over the period s
L_integrated = zeros(1, NT);
for t = 1:NT
    if t >= 2
        L_integrated(t) = trapz(1:t, L_metricsArray(1:t)) / t;
    else
        L_integrated(t) = L_metricsArray(t);
    end
end

% Export L metrics to CSV
csvwrite(fullfile(outdir, 'L_metricsArray.csv'), L_integrated);

% Plot L metrics
figure(8);
plot(1:NT, L_integrated, 'LineWidth', 2);
xlabel('Time Steps');
ylabel('L metrics (meters)');
title('Integrated L metrics over Time');
grid on;

% Save the L metrics figure
L_metrics_file = fullfile(outdir, sprintf('L_metrics.png', experience));
saveas(gcf, L_metrics_file);

disp('L metrics calculated and plotted successfully.');

% Continue with the rest of the plotting and analysis
scale = 1 ;
reference = 0.25  ;
ref_str = num2str(reference);
skipX = 1 ;
skipY = 1 ;
lonref = -67; % Adjusted to the center of the new study area
latref = 16.5; % Adjusted to the center of the new study area
headlght = 1.5;
shaftwdt = 0.1;

lat_min = 15; % Minimum latitude of the study area
lat_max = 18.5; % Maximum latitude of the study area
lon_min = -63; % Minimum longitude of the study area
lon_max = -52; % Maximum longitude of the study area

figure(101)
m_proj('UTM','long',[lon_min lon_max],'lat',[lat_min lat_max])
hold on
m_grid
m_usercoast('LigCoast','patch',[.7 .7 .7]); 
m_vecNoClip(scale,lon(1:skipX:end,1:skipY:end),lat(1:skipX:end,1:skipY:end),...
          avgUe(1:skipX:end,1:skipY:end).*msk(1:skipX:end,1:skipY:end),...
          avgVe(1:skipX:end,1:skipY:end).*msk(1:skipX:end,1:skipY:end),...
         'b','shaftwidth',shaftwdt, 'headlength',headlght); 
for idrft = 1 : Ndrft
    m_plot(trjobsx(:,idrft),trjobsy(:,idrft),'k-','linewi',2);
    m_plot(trjobsx(1,idrft),trjobsy(1,idrft),'ko');
    m_plot(trjex(:,idrft),trjey(:,idrft),'g-','linewi',2);
    m_plot(trjex(itEnter(idrft),idrft),trjey(itEnter(idrft),idrft),'go'); 
end
hold off
m_vec(scale,lonref,latref,reference,0,'b','key',[ref_str,' m/s'],...
          'shaftwidth',shaftwdt, 'headlength',headlght);  
xlabel('Longitude');
ylabel('Latitude');
wysiwyg
      
% Save the U_est figure
U_est_file = fullfile(outdir, sprintf('U_est_%s.png', experience));
saveas(gcf, U_est_file);

disp('U_est figure plotted and saved successfully.');

% Extract parameters from .ini file
params = struct();
params.drifters = strdrft;
params.jpi = str2double(cell2mat(par{2}(1)));
params.jpj = str2double(cell2mat(par{2}(2)));
params.jpt = str2double(cell2mat(par{2}(3)));
params.rdt = str2double(cell2mat(par{2}(4)));
params.jtlag = str2double(cell2mat(par{2}(5)));
params.lenfil = str2double(cell2mat(par{2}(6)));
params.NiterOut = str2double(cell2mat(par{2}(7)));
params.NiterInn = str2double(cell2mat(par{2}(8)));
params.velFname = char(par{2}(9));

% Calculate statistics for E_metrics and L_metrics
E_metrics_min = min(E_metricsArray);
E_metrics_max = max(E_metricsArray);
E_metrics_mean = mean(E_metricsArray, 'omitnan');
E_metrics_median = median(E_metricsArray, 'omitnan'); % Calculate median

L_metrics_min = min(L_integrated);
L_metrics_max = max(L_integrated);
L_metrics_mean = mean(L_integrated, 'omitnan');
L_metrics_median = median(L_integrated, 'omitnan'); % Calculate median

% Prepare data for CSV
csvData = {
    'Parameter', 'Value';
    'drifters', params.drifters;
    'jpi', params.jpi;
    'jpj', params.jpj;
    'jpt', params.jpt;
    'rdt', params.rdt;
    'jtlag', params.jtlag;
    'lenfil', params.lenfil;
    'NiterOut', params.NiterOut;
    'NiterInn', params.NiterInn;
    'velFname', params.velFname;
    'E_metrics_min', E_metrics_min;
    'E_metrics_max', E_metrics_max;
    'E_metrics_mean', E_metrics_mean;
    'E_metrics_median', E_metrics_median; % Add median
    'L_metrics_min', L_metrics_min;
    'L_metrics_max', L_metrics_max;
    'L_metrics_mean', L_metrics_mean;
    'L_metrics_median', L_metrics_median; % Add median
};

% Define CSV file name
csvFile = fullfile(outdir, sprintf('rapport.csv', experience(1:end-4)));

% Write to CSV file
fid = fopen(csvFile, 'w');
for i = 1:size(csvData, 1)
    if isnumeric(csvData{i, 2})
        fprintf(fid, '%s,%f\n', csvData{i, 1}, csvData{i, 2});
    else
        fprintf(fid, '%s,%s\n', csvData{i, 1}, csvData{i, 2});
    end
end
fclose(fid);

disp('Metrics and parameters saved to CSV successfully.');
