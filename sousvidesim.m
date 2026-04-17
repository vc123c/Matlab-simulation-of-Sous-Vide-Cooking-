%% ============================================================
%  SOUS VIDE DOUBLE-BAG MODEL  v4
%  EVA vinyl acetate content sweep — tune to 2% surface O2
%
%  Permeability model (log-linear between LDPE and EVA 28%):
%    P0(VA%) = P0_LDPE * exp(VA% * ln(P0_28/P0_LDPE) / 28)
%    Ea(VA%) = 35e3 - 250 * VA%   (linear interp LDPE -> EVA28)
%
%  Evagelos Christofides  |  UCLA CBE  |  2025
% ============================================================

clear; clc; close all;

%% ---- SHARED PARAMETERS -------------------------------------

rho_m   = 1050;   cp_m  = 3500;
k_m     = 0.50;   alpha_m = k_m/(rho_m*cp_m);
D_O2    = 1.6e-9;
L_slab  = 0.0175;       % half-thickness [m]

T_bath  = 60  + 273.15;
h_conv  = 500;

R_gas   = 8.314;
T_ref   = 296.15;       % 23 C
T_cook  = 333.15;       % 60 C

pO2_out = 0.21 * 101325;
C_out   = pO2_out / (R_gas * T_cook);

L_film  = 25e-6;        % 25 um inner bag

t_cook  = 2 * 3600;
N_t     = 400;
tspan   = linspace(0, t_cook, N_t);
t_min   = tspan / 60;
N_x     = 80;
xmesh   = linspace(0, L_slab, N_x);
x_mm    = xmesh * 1000;

%% ---- VA% SWEEP DEFINITION ----------------------------------

P0_LDPE = 3.0e-13;       % mol/m/s/Pa  at 23 C
P0_28   = 1.8e-12;       % EVA 28% VA  at 23 C
Ea_0    = 35e3;          % LDPE  Ea [J/mol]
Ea_28   = 28e3;          % EVA28 Ea [J/mol]

VA_pcts = [5, 9, 14, 18, 21, 24, 28, 33];
n_VA    = numel(VA_pcts);

% Permeability model: log-linear in VA%
P0_fn  = @(va) P0_LDPE * exp(va * log(P0_28/P0_LDPE) / 28);
Ea_fn  = @(va) Ea_0    - (Ea_0 - Ea_28) / 28 * va;
Pk_fn  = @(va) P0_fn(va) * exp(-Ea_fn(va)/R_gas * (1/T_cook - 1/T_ref));
perm_fn = @(va) Pk_fn(va) / L_film;

% Pre-compute permeance for each VA%
perms = arrayfun(perm_fn, VA_pcts);

fprintf('VA%%   P0 (23C)      Perm (60C)     \n');
fprintf('----  -----------  ---------------\n');
for k = 1:n_VA
    fprintf('%3d%%  %.3e   %.3e mol/m2/s/Pa\n', ...
        VA_pcts(k), P0_fn(VA_pcts(k)), perms(k));
end
fprintf('\n');

%% ---- COLORMAP FOR VA% LINES --------------------------------

cmap = cool(n_VA);

%% ---- SOLVE FOR EACH VA% ------------------------------------

sols      = cell(1, n_VA);
eq_O2     = zeros(1, n_VA);    % equilibrium surface O2 [%]

for k = 1:n_VA
    pm = perms(k);
    sols{k} = pdepe(0, @pdefun, @icfun, ...
        @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t, pm, C_out, ...
                                T_bath, h_conv, rho_m, cp_m), ...
        xmesh, tspan);
    C_surf_final = sols{k}(end, end, 2);
    eq_O2(k) = C_surf_final / C_out * 21;
    fprintf('VA%% = %2d%%  |  equil surface O2 = %.3f%%\n', VA_pcts(k), eq_O2(k));
end
fprintf('\n');

%% ---- FIND OPTIMAL VA% (closest to 2%) ----------------------

[~, idx_opt] = min(abs(eq_O2 - 2.0));
VA_opt       = VA_pcts(idx_opt);
fprintf('Closest to 2%%: %d%% VA  ->  %.3f%% surface O2\n\n', VA_opt, eq_O2(idx_opt));

% Fine interpolation to find exact VA% for 2%
if eq_O2(idx_opt) < 2.0 && idx_opt < n_VA
    VA_exact = interp1(eq_O2(idx_opt:idx_opt+1), VA_pcts(idx_opt:idx_opt+1), 2.0);
elseif eq_O2(idx_opt) > 2.0 && idx_opt > 1
    VA_exact = interp1(eq_O2(idx_opt-1:idx_opt), VA_pcts(idx_opt-1:idx_opt), 2.0);
else
    VA_exact = VA_opt;
end
fprintf('Interpolated optimal VA%%: %.1f%% VA\n\n', VA_exact);

%% ---- MESH GRIDS --------------------------------------------

[T_grid, X_grid] = meshgrid(t_min, x_mm);

%% ============================================================
%  FIGURE 1 — Surface O2 sweep + design curve
% ============================================================

fig1 = figure('Position',[60 60 1300 820],'Color','white', ...
    'Name','EVA VA% Sweep');
sgtitle('EVA Vinyl Acetate Content Sweep  —  25 µm Inner Bag  |  Plain Air Outer  |  60 °C', ...
    'FontSize',13,'FontWeight','bold','Color',[0.12 0.12 0.12]);

% ---- Panel 1: Surface O2 vs time, all VA% ------------------
ax1 = subplot(2,2,1);
hold on; grid on; box on;

fill([0 120 120 0],[0 0 2 2],[0.88 0.96 0.88],'EdgeColor','none','FaceAlpha',0.5);
fill([0 120 120 0],[2 2 10 10],[0.96 0.90 0.88],'EdgeColor','none','FaceAlpha',0.3);
text(3, 0.5,  'Too anaerobic','FontSize',8,'Color',[0.25 0.50 0.25]);
text(3, 2.4,  'Rancidity risk','FontSize',8,'Color',[0.65 0.25 0.15]);

for k = 1:n_VA
    pct_surf = sols{k}(:,end,2) / C_out * 21;
    lw = 1.5;
    if VA_pcts(k) == VA_opt, lw = 3.0; end
    plot(t_min, pct_surf, '-', 'Color', cmap(k,:), ...
        'LineWidth', lw, 'DisplayName', sprintf('%d%% VA', VA_pcts(k)));
end
yline(2, '--', 'Color',[0 0.65 0], 'LineWidth', 1.8, 'DisplayName','2% target');

xlabel('Time (min)', 'FontSize',10);
ylabel('O_2 at meat surface (%)', 'FontSize',10);
title('Surface O_2 vs Time — VA% Sweep', 'FontWeight','bold','FontSize',11);
legend('Location','east','FontSize',7,'NumColumns',2);
ylim([0 10]); xlim([0 120]);

% ---- Panel 2: Design curve — equil O2 vs VA% ---------------
ax2 = subplot(2,2,2);
hold on; grid on; box on;

% Smooth interpolated curve
VA_fine    = linspace(VA_pcts(1), VA_pcts(end), 200);
eq_O2_fine = interp1(VA_pcts, eq_O2, VA_fine, 'pchip');

fill([VA_pcts(1) VA_pcts(end) VA_pcts(end) VA_pcts(1)], [0 0 2 2], ...
    [0.88 0.96 0.88],'EdgeColor','none','FaceAlpha',0.5);
fill([VA_pcts(1) VA_pcts(end) VA_pcts(end) VA_pcts(1)], [2 2 max(eq_O2_fine)*1.1 max(eq_O2_fine)*1.1], ...
    [0.96 0.90 0.88],'EdgeColor','none','FaceAlpha',0.3);

plot(VA_fine, eq_O2_fine, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',2.0, ...
    'DisplayName','Interpolated curve');
scatter(VA_pcts, eq_O2, 70, cmap, 'filled', 'DisplayName','Simulated points');

yline(2, '--', 'Color',[0 0.65 0], 'LineWidth',1.8, 'DisplayName','2% target');
xline(VA_exact, '--', 'Color',[0.8 0.4 0.0], 'LineWidth',1.5, ...
    'DisplayName', sprintf('%.1f%% VA (optimal)', VA_exact));

text(VA_exact + 0.5, max(eq_O2_fine)*0.6, ...
    sprintf('%.1f%% VA\n→ 2%% O_2', VA_exact), ...
    'FontSize',9,'Color',[0.75 0.35 0.0],'FontWeight','bold');

xlabel('Vinyl Acetate Content (%)', 'FontSize',10);
ylabel('Equilibrium Surface O_2 (%)', 'FontSize',10);
title('Design Curve — VA% vs Equilibrium O_2', 'FontWeight','bold','FontSize',11);
legend('Location','northwest','FontSize',8);
ylim([0 max(eq_O2_fine)*1.15]);
xlim([VA_pcts(1)-1 VA_pcts(end)+1]);

% ---- Panel 3: Permeance vs VA% (semi-log) ------------------
ax3 = subplot(2,2,3);
hold on; grid on; box on;

VA_plog  = linspace(0, 33, 200);
perm_log = arrayfun(perm_fn, VA_plog);

semilogy(VA_plog, perm_log, '-', 'Color',[0.3 0.3 0.7], 'LineWidth',2.0);
semilogy(VA_pcts, perms, 'o', 'MarkerFaceColor',[0.85 0.40 0.05], ...
    'MarkerEdgeColor','none','MarkerSize',8);
semilogy([0 0],[perm_fn(0) perm_fn(0)],'s','MarkerFaceColor',[0.2 0.4 0.8], ...
    'MarkerEdgeColor','none','MarkerSize',10,'DisplayName','LDPE');
xline(VA_exact,'--','Color',[0.8 0.4 0],'LineWidth',1.5);
xlabel('Vinyl Acetate Content (%)', 'FontSize',10);
ylabel('Permeance  [mol/(m²·s·Pa)]', 'FontSize',10);
title('Permeance vs VA% at 60 °C', 'FontWeight','bold','FontSize',11);
text(0.5, perm_fn(0)*2.5, 'LDPE', 'FontSize',8,'Color',[0.2 0.4 0.8]);
text(VA_exact+0.5, perm_fn(VA_exact)*0.4, ...
    sprintf('Optimal: %.1f%%',VA_exact),'FontSize',8,'Color',[0.75 0.35 0]);
xlim([-1 35]);

% ---- Panel 4: Optimal VA% — O2 profile snapshots -----------
ax4 = subplot(2,2,4);
hold on; grid on; box on;

fill([0 L_slab*1000 L_slab*1000 0],[0 0 2 2], ...
    [0.88 0.96 0.88],'EdgeColor','none','FaceAlpha',0.5);
text(1, 0.5, 'Anaerobic interior','FontSize',8,'Color',[0.25 0.50 0.25]);

times_snap = [1 5 15 30 60 120];
cmap_snap  = copper(numel(times_snap));

pm_opt = perm_fn(VA_opt);
sol_opt = sols{idx_opt};

for i = 1:numel(times_snap)
    [~, tidx] = min(abs(t_min - times_snap(i)));
    pct_prof   = sol_opt(tidx,:,2) / C_out * 21;
    plot(x_mm, pct_prof, '-', 'Color', cmap_snap(i,:), ...
        'LineWidth',2.0, 'DisplayName', sprintf('%d min', times_snap(i)));
end
yline(2,'--','Color',[0 0.65 0],'LineWidth',1.5,'DisplayName','2% target');
xlabel('Depth from centre (mm)','FontSize',10);
ylabel('O_2 (%)','FontSize',10);
title(sprintf('O_2 Profile Through Steak — %d%% VA (%.2f%% equil)', ...
    VA_opt, eq_O2(idx_opt)),'FontWeight','bold','FontSize',11);
legend('Location','northeast','FontSize',7,'NumColumns',2);
ylim([0 5]); xlim([0 L_slab*1000]);

%% ============================================================
%  FIGURE 2 — 3D fields for optimal VA%
% ============================================================

T_opt  = sol_opt(:,:,1);
C_opt  = sol_opt(:,:,2);
pct_3d = C_opt / C_out * 21;

fig2 = figure('Position',[120 80 1300 550],'Color','white', ...
    'Name', sprintf('Optimal EVA %d%% VA — 3D Fields', VA_opt));
sgtitle(sprintf('EVA %d%% VA  (equil: %.2f%% O_2)  —  Optimal Inner Bag Material', ...
    VA_opt, eq_O2(idx_opt)), ...
    'FontSize',13,'FontWeight','bold','Color',[0.85 0.40 0.05]);

% 3D Temperature
ax5 = subplot(1,2,1);
surf(T_grid, X_grid, T_opt'-273.15, 'EdgeColor','none');
shading interp; colormap(ax5, hot);
lighting gouraud; camlight('headlight');
view(-38,28);
xlabel('Time (min)','FontSize',10);
ylabel('Depth (mm)','FontSize',10);
zlabel('Temp (°C)','FontSize',10);
title('Temperature Field  T(x,t)','FontWeight','bold','FontSize',12);
zlim([0 65]); clim([0 65]);
cb = colorbar('Location','eastoutside');
cb.Label.String = 'T (°C)'; cb.FontSize=8;
grid on;

% 3D O2
ax6 = subplot(1,2,2);
surf(T_grid, X_grid, pct_3d', 'EdgeColor','none');
shading interp; colormap(ax6, parula);
lighting gouraud; camlight('headlight');
view(-38,28);
xlabel('Time (min)','FontSize',10);
ylabel('Depth (mm)','FontSize',10);
zlabel('O_2 (%)','FontSize',10);
title(sprintf('O_2 Field  C(x,t)  —  %d%% VA', VA_opt),'FontWeight','bold','FontSize',12);
zlim([0 5]); clim([0 5]);
cb = colorbar('Location','eastoutside');
cb.Label.String = 'O_2 (%)'; cb.FontSize=8;
grid on;
hold on;
% 2% target plane
[tg, xg] = meshgrid([0 t_cook/60],[0 L_slab*1000]);
mesh(tg, xg, 2*ones(size(tg)), ...
    'FaceAlpha',0,'EdgeColor',[0.15 0.75 0.15],'LineStyle','--');

%% ---- PDE FUNCTIONS -----------------------------------------

function [c, f, s] = pdefun(~, ~, u, DuDx)
    alp   = 0.50 / (1050*3500);
    D     = 1.6e-9;
    kox0  = 2e-4; Eaox = 50e3; Rg = 8.314; Tr = 296.15;
    kox   = kox0 * exp(-Eaox/Rg * (1/max(u(1),250) - 1/Tr));
    c = [1;1];
    f = [alp*DuDx(1); D*DuDx(2)];
    s = [0; -kox*u(2)];
end

function u0 = icfun(~)
    u0 = [5+273.15; 0];
end

function [pl,ql,pr,qr] = bcfun(~,~,~,ur,~, perm, C_out, T_bath, h_conv, rho_m, cp_m)
    pl = [0;0]; ql = [1;1];
    pr = zeros(2,1); qr = ones(2,1);
    pr(1) = (h_conv/(rho_m*cp_m)) * (ur(1) - T_bath);
    pr(2) = perm * (ur(2) - C_out);
end