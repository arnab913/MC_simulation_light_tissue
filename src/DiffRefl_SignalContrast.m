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

        RR = Diff_refl{k}.map;   % log10 reflectance map

        % Bin centers
        xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
        yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

        % Central y ≈ 0
        [~, iy0] = min(abs(yCenters));
        Rline = RR(:, iy0);

        % Smooth
        Rline = smoothdata(Rline, 'sgolay', 1);

        % x ≥ 0 only
        idxPos = xCenters >= 0;
        xPos   = xCenters(idxPos);
        Rpos   = Rline(idxPos);

        plot(xPos, Rpos, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('File %d', k));
    end

    xlabel('Radial distance r (cm)');
    ylabel('log_{10} diffuse reflectance');
    title(groupTitles{g});
    legend('Location','best');
    grid on;
    xlim([0 6]);
    ylim ([-2 4]);
end

sgtitle('Diffuse reflectance profiles along central axis', 'FontSize', 14);


%% signal contrast for each condition

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

    % ----- Baseline -----
    baseIdx = groups{g,1};
    RR_base = Diff_refl{baseIdx}.map;

    % Bin centers
    xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
    yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
    [~, iy0] = min(abs(yCenters));

    R_base = RR_base(:, iy0);
    R_base = smoothdata(R_base, 'sgolay', 30);

    idxPos = xCenters >= 0;
    xPos   = xCenters(idxPos);
    RbPos  = R_base(idxPos);

    % ----- Signal files -----
    sigFiles = groups{g,2};

    for k = sigFiles
        RR_sig = Diff_refl{k}.map;
        R_sig  = RR_sig(:, iy0);
        R_sig  = smoothdata(R_sig, 'sgolay', 1);

        RsPos = R_sig(idxPos);

        % ---- Contrast (log scale difference) ----
        Contrast = RsPos - RbPos;

        plot(xPos, Contrast, 'LineWidth', 1.8, ...
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

sgtitle('Ink-induced contrast relative to baseline', 'FontSize', 14);


%% signal contrast for each ink strength

figure('Color','w','Position',[100 100 1600 900]);

% ---- Ink conditions ----
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

        % ---- Extract maps ----
        R_sig  = Diff_refl{sigIdx(i)}.map;
        R_base = Diff_refl{baseIdx(i)}.map;

        % ---- Bin centers ----
        xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
        yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

        % ---- Central y ≈ 0 line ----
        [~, iy0] = min(abs(yCenters));
        Rs = smoothdata(R_sig(:,iy0),  'sgolay', 1);
        Rb = smoothdata(R_base(:,iy0), 'sgolay', 1);

        % ---- Keep x ≥ 0 ----
        idx = xCenters >= 0;
        x   = xCenters(idx);
        Contrast = Rs(idx) - Rb(idx);   % log-scale contrast

        plot(x, Contrast, 'LineWidth', 2, ...
            'Color', colors(i,:), ...
            'DisplayName', labels{i});
    end

    xlabel('Radial distance r (cm)');
    ylabel('\Delta log_{10}(Reflectance)');
    title(inkCases{s,3});
    yline(0,'k--');
    grid on;
    xlim([0 4]);
    ylim ([-2 0.5]);
    legend('Location','best');
end

sgtitle('Effect of dye depth on ink contrast at fixed absorption', 'FontSize', 14);
