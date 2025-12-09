% plot mat files extracting different meaningful information


%% Plot all together

ua_ink =[1,7,30,100];
color_list = [
   0.85 0.33 0.10;   % orange
   0.13 0.65 0.14;   % green
   0.14 0.44 0.82;   % blue
   0.75 0.15 0.55;   % purple
   0.55 0.52 0.09    % olive
]; % 5 unique colors (rows are RGB)

files = dir('*.mat');
figure;
legendLabels = cell(1, numel(files));
for k = 1:numel(files)
    data = load(files(k).name);
    splitnum = regexp(files(k).name, '\d+', 'match');
    if ~isempty(splitnum)
        ua = str2double(splitnum{1}); % Assumes first number is ua
        color_idx = find(ua_ink == ua, 1);
        if isempty(color_idx)
            plot_color = [0.5 0.5 0.5]; % if not found, use gray
        else
            plot_color = color_list(color_idx, :);
        end
    else
        plot_color = [0.5 0.5 0.5]; % fallback gray
    end
    % Solid for ink, dashed for baseline
    if contains(files(k).name, 'Baseline', 'IgnoreCase', true)
        style = 's--';
    else
        style = 'o-';
    end
    semilogy(data.all_det_r_cm, data.reflectance, style, 'Color', plot_color, 'LineWidth', 1.8);
    hold on;
    legendLabels{k} = files(k).name;
end
hold off;
xlabel('Source–detector distance (cm)');
ylabel('Intensity (a.u.)');
title('Baseline Only Tissue');
legend(legendLabels, 'Location', 'best');
grid on;
setJournalFigure(gcf, 'Interpreter', 'tex');

%% plot signal difference for 0.05% ink
%%%% plot 
b1= load("C:\Users\Arnab Paul.DESKTOP-R4DL8R7\OneDrive - University Of Houston\Documents\MATLAB\MC Simulation\MC_matlab\24Nov_25\AllPlots\Tissue_us 1 Baseline.mat");
b2= load("C:\Users\Arnab Paul.DESKTOP-R4DL8R7\OneDrive - University Of Houston\Documents\MATLAB\MC Simulation\MC_matlab\24Nov_25\AllPlots\DyeGradientTissue_ 1 Baseline.mat"); 
s1= load("C:\Users\Arnab Paul.DESKTOP-R4DL8R7\OneDrive - University Of Houston\Documents\MATLAB\MC Simulation\MC_matlab\24Nov_25\AllPlots\Ink Tissue_ 1 ua ink signal 0.3cm radius.mat");
s2= load("C:\Users\Arnab Paul.DESKTOP-R4DL8R7\OneDrive - University Of Houston\Documents\MATLAB\MC Simulation\MC_matlab\24Nov_25\AllPlots\DyeGradientTissue_ 1 cm-1 ink ua.mat");

figure;
bar(s1.all_det_r_cm,b2.reflectance-s2.reflectance,'BarWidth',0.8, 'FaceColor','y','DisplayName',("Signal contrast after dye")); hold on;
bar(s1.all_det_r_cm,b1.reflectance-s1.reflectance, 'BarWidth',0.5,'FaceColor','b','DisplayName',("Signal contrast No dye"));
xlim([4 6]);
xlabel('SD Distance'); ylabel(' Difference ( corresponding baseline-signal)');
title ('Difference before and after dye for 0.05%ink (\mu_a= 1cm^-1)');
set(gca,'FontSize',16);
legend show;

%% plot baseline -signals

color_list = [
   0.85 0.33 0.10;   % orange
   0.13 0.65 0.14;   % green
   0.14 0.44 0.82;   % blue
   0.75 0.15 0.55;   % purple
   0.55 0.52 0.09    % olive
]; % 5 unique colors (rows are RGB)

figure;
for k=1:numel(files)

    color_idx = mod(k-1, size(color_list,1)) + 1;
    Solid for ink, dashed for baseline
    if contains(files(k).name, 'Baseline', 'IgnoreCase', true)
        style = 's--';
    else
        style = 'o-';
    end
datafile=load(files(k).name);
semilogy(datafile.all_det_r_cm, datafile.reflectance,style,'Color',color_list(color_idx, :),'DisplayName',erase(files(k).name, '.mat'), 'LineWidth',2);
hold on;
end
xlabel('SD distance');
ylabel('Intensity (a.u.)');
legend show;
title('Intensity over distance in diffuse reflectance');

%% Baseline- signal
Dye = load('Dye_alldata.mat');
dye_fields = fieldnames(Dye);
Dye_sorted = struct();  % One struct for all sorted dye data

for f = 1:numel(dye_fields)
    thisstruct = Dye.(dye_fields{f});
    pos_first = arrayfun(@(s) s.Position(1), thisstruct);
    [~, sort_idx] = sort(pos_first);
    Dye_sorted.(dye_fields{f}) = thisstruct(sort_idx);
end

save('Dye_sorted_struct.mat', 'Dye_sorted');
% To use later: S = load('Dye_sorted_struct.mat'); Dye_sorted = S.Dye_sorted;

ND = load('ND_alldata.mat');
nd_fields = fieldnames(ND);
ND_sorted = struct();  % One struct for all sorted ND data

for f = 1:numel(nd_fields)
    thisstruct = ND.(nd_fields{f});
    pos_first = arrayfun(@(s) s.Position(1), thisstruct);
    [~, sort_idx] = sort(pos_first);
    ND_sorted.(nd_fields{f}) = thisstruct(sort_idx);
end

save('ND_sorted_struct.mat', 'ND_sorted');
% To use later: S = load('ND_sorted_struct.mat'); ND_sorted = S.ND_sorted;



%% plot for Dye

baseY = arrayfun(@(s) s.Position(2), Dye_sorted.Dye_base);
baseX = arrayfun(@(s) s.Position(1), Dye_sorted.Dye_base); % SD positions

ink_names    = {'Dye_1cmink','Dye_7cmink','Dye_30cmink','Dye_100cmink'};
ink_conc     = [1, 7, 30, 100]; % Ink concentration values (as x-axis)
n_pos        = numel(baseX);    % Number of SD positions

diff_matrix  = zeros(n_pos, numel(ink_names)); % Each row: SD distance, each col: ink conc

for k = 1:numel(ink_names)
    inkY = arrayfun(@(s) s.Position(2), Dye_sorted.(ink_names{k}));
    diff_matrix(:,k) = inkY - baseY;
end

figure; hold on;
colors = lines(n_pos); % Distinct color for each SD distance

for sd = 1:n_pos
    semilogy(ink_conc, diff_matrix(sd,:), 'o-', ...
        'Color', colors(sd,:), ...
        'LineWidth', 2, ...
        'DisplayName', ['SD Dye = ' num2str(baseX(sd)) ' cm']);
end
xlabel('Ink concentration (cm^{-1})');
ylabel('Reflectance difference (Ink - Baseline)');
title('Ink concentration vs baseline-subtracted signal (by SD distance)');
legend show; grid on; hold off;



%% for No dye data 
% Extract baseline values
baseY_ND = arrayfun(@(s) s.Position(2), ND_sorted.ND_base);
baseX_ND = arrayfun(@(s) s.Position(1), ND_sorted.ND_base); % SD positions

% Ink concentrations and field names
ink_names_ND = {'ND_1cmink','ND_7cmink','ND_30cmink','ND_100cmink'};
ink_conc_ND  = [1, 7, 30, 100]; % x-axis: ink concentration
n_pos_ND     = numel(baseX_ND);

diff_matrix_ND = zeros(n_pos_ND, numel(ink_names_ND));

for k = 1:numel(ink_names_ND)
    inkY_ND = arrayfun(@(s) s.Position(2), ND_sorted.(ink_names_ND{k}));
    diff_matrix_ND(:,k) = inkY_ND - baseY_ND;
end

% Plot on log-y axis
figure; hold on;
colors_ND = lines(n_pos_ND);
for sd = 1:n_pos_ND
    semilogy(ink_conc_ND, diff_matrix_ND(sd,:), 'o-', ...
        'Color', colors_ND(sd,:), ...
        'LineWidth', 2, ...
        'DisplayName', ['SD No dye = ' num2str(baseX_ND(sd)) ' cm']);
end
xlabel('Ink concentration (cm^{-1})');
ylabel('Reflectance difference (Ink - Baseline)');
title('ND: Ink conc. vs baseline-subtracted signal (by SD distance)');
legend show; grid on; hold off;


