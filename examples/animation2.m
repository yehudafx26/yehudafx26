function [] = animation2(aniPlot, vid, frameStep, Nf, k, n, N1, N2, u1, u2, h1, h2, xf, L1, ...
    uf, ub, xb, fn, xcoup, Icoup1 )

    if aniPlot == 1
    
    if mod(n,frameStep) == 0

    figure(1); clf;
    hold on;

    % --- String ---
    z_vec1 = (0:N1)*h1;
    z_vec2 = (0:N2)*h2;

    h_s1 = plot([0;u1;0]*10, z_vec1, '-k', 'LineWidth', 1.5);
    h_s2 = plot([0;u2;0]*10, z_vec2, '--k', 'LineWidth', 1.5);
    
    % --- Finger (红点) ---
    ew = 2*30e-4; 
    eh = 0.03 * L1; 
    h_f = rectangle('Position', [uf*3, xf(n)*L1 - eh/2, ew, eh], 'Curvature' ...
        , [1 1], 'FaceColor', 'r', 'EdgeColor', 'r');

    % --- Bridge (黑三角) ---
    xC = xcoup * L1;
    idx = round(xC/h1);  
    h_brd = plot(u1(idx)*10-0.001, xC, 'k>', ...
        'MarkerSize',10, 'MarkerFaceColor','k');

    % --- Bow (蓝线) ---
    h_bow = yline(xb*L1, 'b', 'LineWidth', 2);

    % --- Axes ---
    xlim([-5e-2, 5e-2]); 
    ylim([0, L1]);
    set(gca, 'YDir', 'reverse');

    % ❌ 关闭坐标 & 网格
    axis off
    set(gca, 'visible', 'off')

    h_f_legend = plot(nan, nan, 'ro', 'MarkerFaceColor','r');

    % --- Legend（关键）---
    legend([h_s1, h_s2, h_bow, h_brd, h_f_legend], ...
    {'String 1', 'String 2', 'Bow', 'Bridge', 'Finger'}, ...
    'Location', 'northwest', ...          % 左侧
    'Orientation', 'vertical', ...   % 竖直
    'FontSize', 13, ...              % 字体变大
    'FontWeight', 'normal');         % 不加粗

    drawnow;


    frame = getframe(gcf);
    writeVideo(vid, frame);

    end
    end

end

