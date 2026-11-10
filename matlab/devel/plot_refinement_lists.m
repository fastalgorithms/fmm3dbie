
function plot_refinement_lists( ...
    S1,S2,refinelist1,refinelist2,ifplot,iteration)

if ~ifplot
    return
end

plotlist1 = zeros(S1.npatches,1);
plotlist2 = zeros(S2.npatches,1);
plotlist1(refinelist1) = 1;
plotlist2(refinelist2) = 1;

figure(3)
clf
plot(S1,plotlist1)
hold on
plot(S2,plotlist2)
hold off
colorbar
title(sprintf('Patches selected after AMR iteration %d',iteration))

end
