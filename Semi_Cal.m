%% Calculate the Fidelity and phase shift by changing N and g0/kappa
clear
clc

n = 1000;
N = 10;
rmax = 0.4;
r = linspace(0, rmax, n);
F = zeros(N, n);
phi = zeros(N, n);

for i = 1:N
    for j = 1:n
        [F(i, j), phi(i, j)] = MB_Fidelity(4 * r(j), i);
    end
end

%% Select the pi point from the above plot and
x = []; % g0/kappa
x_idx = []; % index for g0/kappa
x_F = []; % save the Fidelity for certain g0/kappa
y = []; % number N
for i = N:-1:1
    for j = 1:n-1
        if mod(phi(i,j),pi) - mod(phi(i,j+1),pi) > pi / 2
            x = [x r(j)];
            x_idx = [x_idx j];
            x_F = [x_F F(i, j)];
            y = [y i];
            break
        end
    end
end

%% Plots
%close all

figure()
hold all
pcolor(r, 1:N, F); shading flat; hold on;
colormap('hot')
caxis([0 1])
[C,h]=contour(r, 1:N, F,[0.999 0.99 0.95 0.9 0.5], 'Color',[0 0 0],'Linewidth',2); 
%clabel(C,h)
xlabel('g_0/\kappa')
ylabel('N')
set(gca,'YDir','normal')
colorbar()
set(findall(gcf,'-property','FontSize'),'FontSize',18)
plot(x, y, 'k.', 'Markersize', 40)
xlim([0, rmax])
ylim([1, 10])

title('Fidelity')
%%
figure()
hold all
pcolor(r, 1:N, phi/pi); shading flat; hold on;
colormap('hot')
[C,h]=contour(r, 1:N, phi/pi,[0.5 1 1.5 2 2.5 3], 'Color',[1 1 1],'Linewidth',2); 

xlabel('g_0/\kappa')
ylabel('N')
set(gca,'YDir','normal')
colorbar()
set(findall(gcf,'-property','FontSize'),'FontSize',18)
plot(x, y, 'w.', 'Markersize', 40)
xlim([0, rmax])
ylim([1, 10])
title('Phase shift')
caxis([0 7])
%%
figure()
plot(y, x_F, '-o')
grid on
xlabel('N')
ylabel('Fidelity')

%%

%% Together

ff=figure()
hold all
subplot(211)
pcolor(r, 1:N, 1-F); shading flat; hold on;
colormap('hot')
caxis([0 1])
[C,h]=contour(r, 1:N, F,[0.999 0.99 0.95 0.9 0.5], 'Color',[1 1 1],'Linewidth',2); 
%clabel(C,h)
ylabel('N')
set(gca,'YDir','normal')
colorbar()
set(findall(gcf,'-property','FontSize'),'FontSize',18)
plot(x, y, 'w.', 'Markersize', 40)
xlim([0, rmax])
ylim([1, 10])

set(gca,'XTick',[]);
ylim([1 10])


subplot(212)
hold all
pcolor(r, 1:N, phi/pi); shading flat; hold on;
colormap('hot')
[C,h]=contour(r, 1:N, phi/pi,[0.5 1 1.5 2 2.5 3], 'Color',[1 1 1],'Linewidth',2); 

xlabel('g_0/\kappa')
ylabel('N')
set(gca,'YDir','normal')
colorbar()
set(findall(gcf,'-property','FontSize'),'FontSize',18)
plot(x, y, 'w.', 'Markersize', 40)
xlim([0, rmax])
ylim([1, 10])
caxis([0 7])

set(ff, 'Position', [100, 100, 300, 500]);