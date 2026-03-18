function plotCorr(posterior, saveFigs)
plotParamNames_nounit = {
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{InSAR}}$', ...
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{GPS}}$', ...
  '$V_{\mathrm{HMM}}$', ...
  '$x_{\mathrm{HMM}}$', ...
  '$y_{\mathrm{HMM}}$', ...
  '$d_{\mathrm{HMM}}$', ...
  '$\alpha_{\mathrm{HMM}}$', ...
  '$x_{\mathrm{SC}}$', ...
  '$y_{\mathrm{SC}}$', ...
  '$d_{\mathrm{SC}}$', ...
  '$\alpha_{\mathrm{SC}}$', ...
  '$\phi_{\mathrm{SC}}$', ...
  '$\psi_{\mathrm{SC}}$', ...
  '$\Delta V_{\mathrm{SC}}^{\mathrm{InSAR}}$', ...
  '$\Delta V_{\mathrm{SC}}^{\mathrm{GPS}}$', ...
  '$V_{\mathrm{SC}}$',
};
% assume posterior is already (nParams × nSamples)
nparams = size(posterior,1);
f = figure(12); clf;
tl = tiledlayout(nparams,nparams,'Padding','tight','TileSpacing','tight');
set(gcf,'Color','w');  % white background
set(f, 'Units', 'pixels'); 
set(f, 'Position', [0, 0, 1500, 2000]);

numBins = 30;

for i = 1:nparams
    for j = 1:nparams
        ax = nexttile(tl);
        hold(ax,'on');
        set(ax,'FontSize',8);
        if i > j
            % 2D density in lower triangle
            postx = posterior(j,:);
            posty = posterior(i,:);
            [N, edgesX, edgesY] = histcounts2(postx,posty,numBins);
            centersX = (edgesX(1:end-1)+edgesX(2:end))/2;
            centersY = (edgesY(1:end-1)+edgesY(2:end))/2;
            [X,Y] = meshgrid(centersX,centersY);
            contourf(ax, X, Y, N.', 10, 'LineColor','none');
            colormap(ax,'hot');
            grid(ax,'on');
            ax.YAxis.Visible = 'off';
        elseif i == j
            % 1D histogram on diagonal
            histogram(ax, posterior(i,:), numBins, 'FaceColor',[.2 .6 .5]);
            % title(ax, plotParamNames_nounit{i}, ...
            %       'Interpreter','latex','FontSize',16);
            text(ax, 0.5, 1.45, plotParamNames_nounit{i}, ...
                'Units', 'normalized', ... 
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'latex', ...
                'FontSize', 28, ... 
                'FontWeight', 'bold');
            ax.YAxis.Visible = 'off';
        else
            % blank above diagonal
            axis(ax,'off');
            continue
        end

        % restore tick labels everywhere
        ax.XAxis.Visible = 'off';



        % only bottom row: add xlabel
        if i == nparams
            ax.XLabel.String = plotParamNames_nounit{j};
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.FontSize = 20;
            ax.XAxis.Visible = 'on';
        else
            ax.XLabel.String = '';
        end

        % only first column: add ylabel
        if j == 1
            ax.YLabel.String = plotParamNames_nounit{i};
            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.FontSize = 18;
            ax.YAxis.Visible = 'on';
        else
            ax.YLabel.String = '';
        end
    end
end

sgtitle(tl, "2D Density (lower triangle) & Histograms (diag.) for MCMC", ...
        'FontWeight','normal', 'FontSize', 30);

if saveFigs
    exportgraphics(tl, './PaperFigs/corr_plot.png','Resolution',500);
end



end