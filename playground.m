clear

% GENERAL FLOW CONSTANTS
lx     = 40;      % number of cells in x-direction
ly     = 10;      % number of cells in y-direction
obst_x = lx/5+1;   % position of the cylinder; (exact
obst_y = ly/2+3;   % y-symmetry is avoided)
obst_r = ly/10+1;  % radius of the cylinder
uMax   = 0.1;      % maximum velocity of Poiseuille inflow
Re     = 100;      % Reynolds number
nu     = uMax * 2.*obst_r / Re;  % kinematic viscosity
omega  = 1. / (3*nu+1./2.);      % relaxation parameter
maxT   = 400;  % total number of iterations
tPlot  = 50;      % cycles

% D2Q9 LATTICE CONSTANTS
t  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];
col = [2:(ly-1)];
in  = 1;   % position of inlet
out = lx;  % position of outlet

[y,x] = meshgrid(1:ly,1:lx); % get coordinate of matrix indices
  
obst = ...                   % Location of cylinder
    (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ly]) = 1;    % Location of top/bottom boundary
bbRegion = find(obst); % Boolean mask for bounce-back cells

% INITIAL CONDITION: Poiseuille profile at equilibrium
L = ly-2;
y_phys = y-1.5;
ux = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
uy = zeros(lx,ly);
rho = 1;

for i=1:9
    cu = 3*(cx(i)*ux+cy(i)*uy);
    fIn(i,:,:) = rho .* t(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) );
end

% MAIN LOOP (TIME CYCLES)
for cycle = 1:maxT

    % MACROSCOPIC VARIABLES
    rho = sum(fIn);
    ux  = reshape ( (cx * reshape(fIn,9,lx*ly)), 1,lx,ly) ./rho;
    uy  = reshape ( (cy * reshape(fIn,9,lx*ly)), 1,lx,ly) ./rho;

end

rho