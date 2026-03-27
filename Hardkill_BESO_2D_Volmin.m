%--------------------------------------------------------------------------%
% Hard-kill BESO code for volume minimization with compliance constraint   %
% Constraint: C_optimal = C0 * C_target                                    %
% This code is written by Md Zakirul for educational purpose only          %
%--------------------------------------------------------------------------%

function Hardkill_BESO_2D_Volmin(nelx,nely,C_target,er,rmin)
%% INITIALIZE
nelx = 100; nely = 50;
E0 = 1; nu = 0.3; C_target = 1.6;  % Target compliance = C0 * C_target
er = 0.01; rmin = 5;

%% START ITERATION
x = ones(nely,nelx); vol = 1; iter = 0; change = 1; c = []; index = 0;
C0 = [];  % Initial compliance (will be computed in first iteration)
comp_con = [];  % Absolute compliance constraint
C_ratio_hist = [];  % History of compliance ratio
V_hist = [];  % History of volume fraction

%% INDEXING NODES AND ELEMENT
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);

%% STIFFNESS MATRIX
k = [1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = 1/(1-nu^2)*[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

%% DEFINE LOADS AND SUPPORTS
% Load at center bottom
load_x = round(nelx/2) + 1;  % Center horizontally
load_y = nely + 1;           % Bottom edge (last node in column)
load_node = (load_x-1) * (nely+1) + load_y;
F = sparse(2*load_node, 1, -1, 2*(nely+1)*(nelx+1), 1);  % Downward force
U = zeros(2*(nely+1)*(nelx+1), 1);
% Left support - PINNED (bottom-left corner)
left_node = nely + 1;  % Bottom-left corner node
fixeddofs_left = [2*left_node-1, 2*left_node];  % Fix both x and y
% Right support - ROLLER (bottom-right corner, Node A)
right_node = nelx * (nely+1) + (nely+1);  % Bottom-right corner node
fixeddofs_right = 2*right_node;  % Fix only y direction
% Combine all fixed DOFs
fixeddofs = [fixeddofs_left, fixeddofs_right];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

%% PREPARE for 2D FILTERING
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        [i2,j2] = ndgrid(max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx),max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely));
        e2 = (i2(:)-1)*nely+j2(:);
        iH(index + (1:numel(e2))) = e1;
        jH(index + (1:numel(e2))) = e2;
        sH(index + (1:numel(e2))) = max(0,rmin-sqrt((i1-i2(:)).^2+(j1-j2(:)).^2));
        index = index + numel(e2);
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);

%% MAIN OPTIMIZATION LOOP
while change > 0.001
    iter = iter + 1;
    if iter >1; olddc = dc; end
    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*E0*x(:)',64*nelx*nely,1);
    K = sparse(iK,jK,sK);  K = (K+K')/2;
    U(freedofs) = decomposition(K(freedofs,freedofs),'chol','lower')\F(freedofs);  
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = [c, 0.5.*sum(sum(x*E0.*ce))];
    % COMPUTE INITIAL COMPLIANCE AND SET CONSTRAINT
    if iter == 1
        C0 = c(1);
        comp_con = C0 * C_target;
        disp(['Initial Compliance C0 = ' sprintf('%6.4f',C0)])
        disp(['Target Compliance = C0 * ' sprintf('%4.2f',C_target) ' = ' sprintf('%6.4f',comp_con)])
        disp('-----------------------------------------------------------')
    end
    % RECORD HISTORY
    C_ratio_hist = [C_ratio_hist, c(iter)/C0];
    V_hist = [V_hist, sum(sum(x))/(nelx*nely)];
    
    % SENSITIVITY ANALYSIS FOR VOLUME MINIMIZATION
    dv = ones(nely,nelx);  % Sensitivity of volume w.r.t. density
    dc = -x.*ce;           % Sensitivity of compliance w.r.t. density (negative)
    % Filter sensitivities
    dv(:) = H*dv(:)./Hs;
    dc(:) = H*dc(:)./Hs;
    if iter > 1; dc = (dc+olddc)/2.; end  % STABILIZATION
    % ADAPTIVE VOLUME UPDATE WITH TOLERANCE BAND
    current_comp = c(iter);
    comp_tolerance = 0.01;  % 2% tolerance band around target
    if current_comp > comp_con * (1 + comp_tolerance)
        % Compliance significantly exceeds constraint - add material
        vol = min(vol*(1+er), 1.0);
    elseif current_comp < comp_con * (1 - comp_tolerance)
        % Compliance well below constraint - remove material
        vol = max(vol*(1-er), 0.01);  % Minimum 1% volume
    else
        % Within tolerance band - hold current volume to stabilize
        vol = vol;
    end
    % BESO DESIGN UPDATE - VOLUME MINIMIZATION
    l1 = min(min(-dc./dv)); l2 = max(max(-dc./dv));
    while ((l2-l1)/l2 > 1e-5)
        th = (l1+l2)/2;
        x = max(1e-9,sign(-dc./dv-th));
        if sum(sum(x))-vol*(nelx*nely) > 0
            l1 = th;
        else
            l2 = th;
        end
    end
    if iter > 10
        change = abs(mean(V_hist(iter-4:iter)) - mean(V_hist(iter-9:iter-5))) / mean(V_hist(iter-4:iter));
    end
    % Calculate compliance metrics
    comp_error = abs(c(iter) - comp_con) / comp_con;
    comp_ratio = c(iter) / C0;
    disp([' It.: ' sprintf('%3i',iter) ...
          ' C.: ' sprintf('%6.4f',c(iter)) ...
          ' C/C0: ' sprintf('%5.3f',comp_ratio) ...
          ' Vol.: ' sprintf('%5.3f',sum(sum(x))/(nelx*nely)) ...
          ' Target: ' sprintf('%5.3f',C_target) ...
          ' Err.: ' sprintf('%5.2f%%',comp_error*100) ...
          ' ch.: ' sprintf('%6.4f',change)])
    clf; colormap(summer); imagesc(x); axis equal tight off; pause(1e-6);
end

%% FINAL SUMMARY
disp('===========================================================')
disp(['OPTIMIZATION COMPLETE'])
disp(['Initial Compliance C0    = ' sprintf('%6.4f',C0)])
disp(['Target Ratio C_target    = ' sprintf('%5.3f',C_target)])
disp(['Target Compliance        = ' sprintf('%6.4f',comp_con)])
disp(['Final Compliance         = ' sprintf('%6.4f',c(end))])
disp(['Final Compliance Ratio   = ' sprintf('%5.3f',c(end)/C0)])
disp(['Final Volume Fraction    = ' sprintf('%5.3f',sum(sum(x))/(nelx*nely))])
disp('===========================================================')

%% HISTORICAL PLOT
loop = iter;
iter_vec = 1:loop;
figure; clf;
set(gcf,'Color','w','Units','centimeters','Position',[3 3 18 12]);
yyaxis left
h1 = plot(iter_vec, C_ratio_hist(1:loop), 'o-', 'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74], ...
'MarkerEdgeColor', [0.00 0.45 0.74], 'MarkerSize', 6, 'LineWidth', 2);
hold on;
yline_left = yline(C_target, 'r--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
set(gca,'YScale','log');
ylabel('Compliance Ratio', 'FontSize', 24, 'FontWeight','bold');
set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.00 0.45 0.74], 'XColor', 'k');
yyaxis right
h2 = plot(iter_vec, V_hist(1:loop), 's-', 'Color', [0.49 0.18 0.56], 'MarkerFaceColor', 'none', ...
'MarkerEdgeColor', [0.49 0.18 0.56], 'MarkerSize', 6, 'LineWidth', 2);
ylabel('Volume Fraction', 'FontSize', 24, 'FontWeight','bold');
set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.49 0.18 0.56]);
xlabel('Iteration', 'FontSize', 24, 'FontWeight','bold');
title('Convergence History', 'FontSize', 24, 'FontWeight','bold');
legend([h1, h2, yline_left], ...
    {'Compliance Ratio $C/C_0$','Volume Fraction','Target Compliance'}, 'Location','best', 'FontSize', 24, 'Box','on', 'Interpreter','latex');
cmin = min(C_ratio_hist(1:loop));
cmax = max(C_ratio_hist(1:loop));
vmin = min(V_hist(1:loop));
vmax = max(V_hist(1:loop));
yyaxis left
ylim([cmin*0.9, cmax*1.1]);
yyaxis right
ylim([vmin - 0.05*(vmax-vmin), vmax + 0.05*(vmax-vmin)]);
grid on;
yyaxis left
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
yyaxis right
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
drawnow;