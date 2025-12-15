
%% plot for all conditions no dye, dye 1 and dye 2
groups = {
    1:5
    6:10
    11:15
};

groupTitles = {
    'No dye (varying ink)'
    '1-layer dye (varying ink)'
    '2-layer dye (varying ink)'
};

figure('Color','w','Position',[100 100 1600 450]);

for g = 1:3
    subplot(1,3,g); hold on;

    for k = groups{g}
        RR = Diff_refl{k}.map;

        [r, Rr] = radialAverage(RR, xEdges, yEdges);

        plot(r, Rr, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('File %d', k));
    end

    xlabel('Radial distance r (cm)');
    ylabel('log_{10} diffuse reflectance');
    title(groupTitles{g});
    legend('Location','best');
    grid on;
    xlim([0 6]);
    ylim([-2 4]);
end

sgtitle('Azimuthally averaged diffuse reflectance R(r)', 'FontSize', 14);


%% contrast vs baseline per dye condition

groups = {
    1,  2:5
    6,  7:10
    11, 12:15
};

groupTitles = {
    'No dye: ink contrast'
    '1-layer dye: ink contrast'
    '2-layer dye: ink contrast'
};

figure('Color','w','Position',[100 100 1600 450]);

for g = 1:3
    subplot(1,3,g); hold on;

    baseIdx = groups{g,1};
    [r, Rb] = radialAverage(Diff_refl{baseIdx}.map, xEdges, yEdges);

    for k = groups{g,2}
        [~, Rs] = radialAverage(Diff_refl{k}.map, xEdges, yEdges);

        Contrast = Rs - Rb;   % log-scale contrast

        plot(r, Contrast, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('File %d - %d', k, baseIdx));
    end

    xlabel('Radial distance r (cm)');
    ylabel('\Delta log_{10} Reflectance');
    title(groupTitles{g});
    legend('Location','best');
    grid on;
    xlim([0 6]);
    yline(0,'k--');
end

sgtitle('Ink-induced contrast (radial average)', 'FontSize', 14);


%% contrast vs dye depth for a ink strength
figure('Color','w','Position',[100 100 1600 900]);

inkCases = {
    [2 7 12],  [1 6 11],  'Ink \mu_a = 1 cm^{-1}'
    [3 8 13],  [1 6 11],  'Ink \mu_a = 7 cm^{-1}'
    [4 9 14],  [1 6 11],  'Ink \mu_a = 30 cm^{-1}'
    [5 10 15], [1 6 11],  'Ink \mu_a = 100 cm^{-1}'
};

colors = lines(3);
labels = {'No dye','1-layer dye','2-layer dye'};

for s = 1:4
    subplot(2,2,s); hold on;

    sigIdx  = inkCases{s,1};
    baseIdx = inkCases{s,2};

    for i = 1:3
        [r, Rs] = radialAverage(Diff_refl{sigIdx(i)}.map,  xEdges, yEdges);
        [~, Rb] = radialAverage(Diff_refl{baseIdx(i)}.map, xEdges, yEdges);

        Contrast = Rs - Rb;

        plot(r, Contrast, 'LineWidth', 2, ...
            'Color', colors(i,:), ...
            'DisplayName', labels{i});
    end

    xlabel('Radial distance r (cm)');
    ylabel('\Delta log_{10}(Reflectance)');
    title(inkCases{s,3});
    yline(0,'k--');
    grid on;
    xlim([0 4]);
    ylim([-2 0.5]);
    legend('Location','best');
end

sgtitle('Effect of dye depth on ink contrast (radial average)', 'FontSize', 14);


%%
function [rCenters, Ravg] = radialAverage(RR, xEdges, yEdges)
% RR : Nx × Ny map (log10 reflectance)
% Returns azimuthally averaged R(r)

xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

[X,Y] = meshgrid(xCenters, yCenters);
Rmap  = RR';                       % match X,Y orientation
Rrad  = sqrt(X.^2 + Y.^2);

dr = mean(diff(xCenters));
rBins = 0:dr:max(xCenters);
rCenters = (rBins(1:end-1) + rBins(2:end)) / 2;

Ravg = nan(size(rCenters));

for i = 1:numel(rCenters)
    mask = (Rrad >= rBins(i)) & (Rrad < rBins(i+1));
    vals = Rmap(mask);
    if ~isempty(vals)
        Ravg(i) = mean(vals,'omitnan');
    end
end
Ravg = Ravg;
%Ravg = smoothdata(Ravg,'sgolay',1);
end
