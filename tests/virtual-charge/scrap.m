
figure

semilogx( [100,32e3], [0, 0], 'k');
hold on
plot( [100,4e3], [-1, -1], 'k');
plot( [32e3,180e3], [-30, -30], 'k');
plot( [180e3,1e6], [-55, -55], 'k');

plot( [4e3,4e3], [-1, -80], 'k');
plot( [32e3,32e3], [0, -30], 'k');
plot( [180e3,180e3], [-30, -55], 'k');
xlim([1e2,1e6])
ylim([-80,10])