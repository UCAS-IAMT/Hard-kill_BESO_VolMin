%--------------------------------------------------------------------------%
% Hard-kill BESO code for volume minimization with compliance constraint   %
% Constraint: C_optimal = C0 * C_target                                    %
% This code is written by Md Zakirul for educational purpose only          %
%--------------------------------------------------------------------------%

function Hardkill_BESO_3D_Volmin(nelx,nely,nelz,C_target,er,rmin)
%% INITIALIZE
nelx = 80;  % 400mm length (5mm per element)
nely = 8;   % 40mm width (5mm per element)
nelz = 10;  % 50mm height (5mm per element)
E0 = 1; nu = 0.3; C_target = 5.0;  % Target compliance = C0 * C_target
er = 0.01;
rmin = 3;

%% START ITERATION
nElem = nely*nelx*nelz;  % TOTAL ELEMENT NUMBER
x = ones(nElem,1);
dv0 = ones(nElem,1);
vol = 1; iter = 0; change = 1; c = []; maxiter = 1000;
C0 = [];  % Initial compliance (will be computed in first iteration)
comp_con = [];  % Absolute compliance constraint
C_ratio_hist = [];  % History of compliance ratio
V_hist = [];  % History of volume fraction

%% INDEXING NODES AND ELEMENT
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nely,1+nelz,1+nelx);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nElem,1);
edofMat = repmat(edofVec,1,24)+repmat([0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(nely+...
    1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+[-3,-2,-1]],nElem,1);

%% STIFFNESS MATRIX
k = 1/(1+nu)/(2*nu-1)/144 *([-32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;-4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;-6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;-3;-32;-6;-3;-4;-32;6;6;-32;-6;-32] ...
    +nu*[48;0;0;0;-24;-24;-12;0;-12;0;24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;0;0;24;12;-12;12;0;-12;0;-12;-12;0;
    48;24;0;0;12;12;-12;0;24;0;-24;-24;0;0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;-12;-12;-12;-12;0;0;48;0;24;0;-24;0;
    -12;-12;-12;-12;12;0;0;24;12;-12;0;0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;0;-24;24;-12;0;0;-12;12;-12;0;0;-24;
    -12;-12;0;48;0;24;0;0;0;-12;0;-12;-12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;
    24;48;0;12;-12;0;0;-12;0;-12;-12;-12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;
    24;48;-24;0;0;-12;-12;-12;0;-24;0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48]);
KE(tril(ones(24))==1) = k';
KE = reshape(KE,24,24);
KE = KE+KE'-diag(diag(KE));
index_ik = arrayfun(@(i) i:24,1:24,'UniformOutput',false);
index_ik = horzcat(index_ik{:});
index_jk = arrayfun(@(j) repmat(j,1,24-j+1),1:24,'UniformOutput',false);
index_jk = horzcat(index_jk{:});
iK = edofMat(:, index_ik)';
jK = edofMat(:, index_jk)';
index_k = sort([iK(:), jK(:)],2,'descend');  clear iK jK;

%% DEFINE LOADS AND SUPPORTS
% Apply load at center of beam (middle of span in x, middle in y, BOTTOM surface in z)
% Center position: x = nelx/2, y = nely/2, z = 1 (BOTTOM surface)
center_x = round(nelx/2) + 1;  % +1 for node indexing
center_y = round(nely/2) + 1;
bottom_z = 1;  % Bottom surface
% Create force vector: F = 1.12 kN downward (negative z-direction)
center_x = round((nelx+1)/2);
center_y = round((nely+1)/2);
bottom_z = 1;  % bottom surface
loadNode = nodenrs(center_y, bottom_z, center_x);
F = sparse(3*(nely+1)*(nelx+1)*(nelz+1),1);
F(3*loadNode) = -1120;   % z-direction (vertical)
% Fixed supports at both ends (x = 0 and x = nelx)
% Left end (x = 0): fix all nodes on this face
left_nodes = [];
for iy = 1:nely+1
    for iz = 1:nelz+1
        node_num = (iz-1)*(nely+1) + iy;
        left_nodes = [left_nodes, 3*node_num-2, 3*node_num-1, 3*node_num];
    end
end
% Right end (x = nelx): fix all nodes on this face  
right_nodes = [];
for iy = 1:nely+1
    for iz = 1:nelz+1
        node_num = nelx*(nely+1)*(nelz+1) + (iz-1)*(nely+1) + iy;
        right_nodes = [right_nodes, 3*node_num-2, 3*node_num-1, 3*node_num];
    end
end
fixeddofs = [left_nodes, right_nodes];
U = zeros(3*(nely+1)*(nelx+1)*(nelz+1),1);
alldofs = 1:3*(nely+1)*(nelx+1)*(nelz+1);
freedofs = setdiff(alldofs,fixeddofs);
disp('=============================================================')
disp('STRUCTURE CONFIGURATION:')
disp(['Dimensions: ', num2str(nelx*5), 'mm × ', num2str(nely*5), 'mm × ', num2str(nelz*5), 'mm'])
disp(['Elements: ', num2str(nelx), ' × ', num2str(nely), ' × ', num2str(nelz)])
disp('Boundary Conditions: Fixed supports at left and right faces (x=0 and x=400mm)')
disp(['Load: F = 1.12 kN at BOTTOM center (x=', num2str((center_x-1)*5), 'mm, y=', num2str((center_y-1)*5), 'mm, z=0mm)'])
disp('=============================================================')

%% PREPARE for 3D FILTERING
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2+dz.^2));
Hs = imfilter(ones(nely,nelz,nelx),h,'symmetric');

%% MAIN OPTIMIZATION LOOP
while change > 0.001 && iter < maxiter
    iter = iter + 1;
    % FINITE ELEMENT ANALYSIS
    sK = reshape(k(:)*E0*x(:)',length(k)*nElem,1);
    K = sparse(index_k(:,1),index_k(:,2),sK);
    LK = chol(K(freedofs,freedofs), 'lower');
    U(freedofs) = LK'\(LK\F(freedofs));
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nElem,1);
    c = [c, 0.5.*sum(x*E0.*ce)];
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
    V_hist = [V_hist, sum(sum(sum(x)))/(nElem)];
    % SENSITIVITY ANALYSIS FOR VOLUME MINIMIZATION
    dc = -x.*ce*E0;  % Sensitivity of compliance w.r.t. density
    dc = imfilter(reshape(dc,nely,nelz,nelx)./Hs,h,'symmetric');
    dv = imfilter(reshape(dv0,nely,nelz,nelx)./Hs,h,'symmetric');
    if iter > 1; dc = (dc+olddc)/2.; end  % STABILIZATION OF EVOLUTIONARY PROCESS
    % ADAPTIVE VOLUME UPDATE BASED ON COMPLIANCE CONSTRAINT
    current_comp = c(iter);
    comp_tolerance = 0.01;  % 1% tolerance band around target
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
    l1 = min(min(min(-dc./dv))); l2 = max(max(max(-dc./dv)));
    while ((l2-l1)/l2 > 1e-5)
        th = (l1+l2)/2;
        x = max(1e-9,sign(-dc./dv-th));
        if sum(sum(sum(x)))-vol*nElem > 0
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
          ' Vol.: ' sprintf('%5.3f',sum(sum(sum(x)))/(nElem)) ...
          ' Target: ' sprintf('%5.3f',C_target) ...
          ' Err.: ' sprintf('%5.2f%%',comp_error*100) ...
          ' ch.: ' sprintf('%6.4f',change)])
    surf = shiftdim(reshape(x,nely,nelz,nelx),2);
    surf = smooth3(surf,'box',1);
    clf; patch(isosurface(surf,0.5),'FaceColor',[0 127 102]/255,'EdgeColor','none','FaceAlpha', 0.9);
    patch(isocaps(surf,0.5),'FaceColor',[255 255 102]/255,'EdgeColor','none','FaceAlpha', 0.7);
    light('Position', [1 1 1], 'Style', 'infinite');
    light('Position', [-1 -1 1], 'Style', 'infinite');
    lighting phong; 
    material([0.5 0.6 0.4]);
    drawnow; view([110,20]); axis equal tight off; pause(1e-6);
    x = reshape(x,nElem,1);
    olddc = dc;
end

%% FINAL SUMMARY
disp('===========================================================')
disp(['OPTIMIZATION COMPLETE'])
disp(['Initial Compliance C0    = ' sprintf('%6.4f',C0)])
disp(['Target Ratio C_target    = ' sprintf('%5.3f',C_target)])
disp(['Target Compliance        = ' sprintf('%6.4f',comp_con)])
disp(['Final Compliance         = ' sprintf('%6.4f',c(end))])
disp(['Final Compliance Ratio   = ' sprintf('%5.3f',c(end)/C0)])
disp(['Final Volume Fraction    = ' sprintf('%5.3f',sum(sum(sum(x)))/(nElem))])
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
set(gca, 'FontSize', 24, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.00 0.45 0.74], 'XColor', 'k');
yyaxis right
h2 = plot(iter_vec, V_hist(1:loop), 's-', 'Color', [0.49 0.18 0.56], 'MarkerFaceColor', 'none', ...
'MarkerEdgeColor', [0.49 0.18 0.56], 'MarkerSize', 6, 'LineWidth', 2);
ylabel('Volume Fraction', 'FontSize', 24, 'FontWeight','bold');
set(gca, 'FontSize', 24, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.49 0.18 0.56]);
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