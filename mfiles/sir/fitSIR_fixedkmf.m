function [p,res,J]= fitSIR_fixedkmf(ti,td,M,Sm,p0,lb,ub,kmf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saved some time by nesting the function signalSIR below. Got rid of if
% loops for plotting flags. Using levenberg-marquardt slightly faster, 
% cannot use bounds..
%
% Fits signal for SIR-TSE sequences to a 2 pool model of MT (with fixed kmf);  
% f = free water, m = macromolecular.
% 
% Input:
% ti - array of inversion times (N x 1)
% td - array of delay times (N x 1)
% M  - measured signal (N x 1)
% Sm - numerically estimaged saturation fraction of m pool due to inv pulse
% p0 - initial parameter guess (see p for description)
% lb,ub - lower and upper bounds for p
% magFlag - magnitude data? ('y' or 'n')
% plotFlag - plot data? ('y' or 'n')
% kmf - fixed exchange rate
% 
% Output:
% p - fitted parameters 
%      p(1) = pmf
%      p(2) = R1f = R1m 
%      p(3) = Sf
%      p(4) = M0f
% ci - corresponding 95% CIs
% res - residuals
% J = Jacobian of FUN at p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Make sure ti/td/M are column vectors
if size(ti,1) == 1
    ti = ti';
end
if size(td,1) == 1
    td = td';
end
if size(M,1) == 1
    M = M';
end


% Perform nonlinear least square fit
% options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'UseParallel',1);
options = optimoptions('lsqnonlin',"Display","off","FunctionTolerance",1E-10,...
    "UseParallel",false,"Algorithm","levenberg-marquardt","StepTolerance",1E-10);

% [p,~,res,~,~,~,J] = lsqnonlin(@optfitSIR,p0,lb,ub,options,ti,td,M,Sm,magFlag,plotFlag,kmf);
% [p,~,res,~,~,~,J] = lsqnonlin(@optfitSIR,p0,[],[],options,ti,td,M,Sm,magFlag,plotFlag,kmf);
[p,~,res,~,~,~,J] = lsqnonlin(@optfitSIR,p0,[],[],options,ti,td,M,Sm,kmf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function er = optfitSIR(p,ti,td,M,Sm,magFlag,plotFlag,kmf)
function er = optfitSIR(p,ti,td,M,Sm,kmf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract parameters
M0 = [1 p(1)]*p(4);
R1 = [p(2) p(2)];
S  = [p(3) Sm];

% Calculate signal
Mm = signalSIR(ti,td,kmf,R1,M0,S);
% if magFlag == 'y'
Mm = abs(Mm);
% end

% Misfit
er = Mm - M;


% Plot data % save some time by not evaluating these if loops.
% if plotFlag == 'y'
%     plot(ti,Mm,'b-',ti,M,'rs','LineWidth',2)
    
%     xlabel('t_i (s)'), ylabel('Signal (a.u.)')
%     grid on, drawnow
% end

function [M,J,bfplus,R1plus,bfminus,R1minus,Mfinf] = signalSIR(ti,td,kmf,R1,M0,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = signalSIR(ti,td,kmf,M0,S)
%
% Calculates signal for SIR-TSE sequences assuming 2 pools undergoing MT;
% f = free water, m = macromolecular.
%
% Input:
% ti - array of inversion times (N x 1)
% td - array of delay times (N x 1)
% kmf - exchange rate from m to f
% M0 - M0 for each pool [M0f Mom]'
% R1 - R1 for each pool [R1f R1m]'
% S -  sat/inv effeciency for each pool [Sf Sm]'
%
% Output:
% M - z-mag [Mf Mm]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No longer using my implentation based upon the notation of:
% Dortch et al. NeuroImage 64 (2013) 640?649
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Make sure vectors is the correct size
% if size(M0,1) == 1
%     M0 = M0';
% end
% if size(S,2) == 1
%     S = S';
% end
%
% % Setup exchange matrix
% kfm = kmf*M0(2)/M0(1);
% K = [-kfm kmf ; kfm -kmf];
%
% % Calculate L1
% L1 = -diag(R1) + K;
%
% % Eigen-expansion of L1
% [U1,S1] = eig(L1);
% V1 = eye(size(U1))/U1;
%
% % Diagonalize S
% S = diag(S);
%
% % Set identitiy matrix
% I = eye(2);
%
% % Calculate signal
% Ei = zeros(size(V1));
% Ed = zeros(size(V1));
% M  = zeros(length(ti),2);
% for ii = 1:length(ti)
%     Ei =  U1*diag(exp(diag(S1*ti(ii))))*V1;
%     Ed =  U1*diag(exp(diag(S1*td(ii))))*V1;
%     M(ii,:) = (Ei*S*(I - Ed) + (I - Ei))*M0;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ke's implementation is faster
% Notation from Li et al. MRM 64:491-500 (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (ti,td,kmf,R1,M0,S);
% M0 = [1 p(1)]*p(4);
% R1 = [p(2) p(2)];
% S  = [p(3) Sm];
% 
% % Calculate signal
% Mm = signalSIR(ti,td,kmf,R1,M0,S);

% Get pmf and Mfinf
pmf = M0(2)/M0(1);
Mfinf = M0(1);

% Get R1f/R1m and Sf/Sm
R1f = R1(1);  
R1m = R1(2);
Sf = S(1);    
Sm = S(2);

% Calculate R1+/- in Eq. 4
R1diff  = sqrt((R1f - R1m + (pmf - 1) * kmf)^ 2 + 4 * pmf * kmf^2);
R1plus  = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2;
R1minus = R1plus - R1diff;

% Component amplitude terms for td terms (Eq. 5)
bftdplus = -(R1f - R1minus) / R1diff;
bftdminus = (R1f - R1plus) / R1diff;
bmtdplus  = -(R1m - R1minus) / R1diff;
bmtdminus = (R1m - R1plus) / R1diff;

% Signal recovery during td (Eq. 5)
Mftd = bftdplus * exp(-R1plus * td) + bftdminus * exp(-R1minus * td) + 1 ;
Mmtd = bmtdplus * exp(-R1plus * td) + bmtdminus * exp(-R1minus * td) + 1 ;

% Component amplitude terms for ti terms (Eq. 5)
bfplus= ((Sf * Mftd -1) * (R1f - R1minus)  +  (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff;
bfminus= -((Sf * Mftd -1) * (R1f - R1plus) + (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff;

% Signal equation (Eq. 3)
M = (bfplus .* exp (-R1plus * ti) + bfminus .* exp (-R1minus * ti) + 1) * Mfinf;