function ax = mySubplot(nrow, ncol)
% returns a matrix of axis handles
% to plot in the second subplot, you would use plot(ax(1,2), x, y)

% standard x, y, dx, dy for subplot(111)
x0 = 0.1300;
y0 = 0.1100;
w0 = 0.7750;
h0 = 0.8150;
w = w0 / ncol;
h = h0 / nrow;
figure()

ax = nan(nrow, ncol);
for irow = 1:nrow
    for icol = 1:ncol
        ax(irow, icol) = axes('position', [x0 + (icol - 1) * w, y0 + (nrow - irow) * h, 0.9*w, 0.9*h]);
        
        if irow ~= nrow
            set(gca, 'XTickLabel', '')
        end
    end
end