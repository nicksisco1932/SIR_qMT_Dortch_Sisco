% High-res SIR analysis for 7T data - Richard Dortch - 4/13/15
clear, close all, clc

cd '/Users/nicholassisco/Documents/Research/Data/SIR7T_NMSS/analysis/mfiles'

% Set path of primary script - the code below assumes you run the script
% from directory
scriptDIR = mfilename('fullpath');
scriptDIR = pwd % NS added 20200206
%
% Add the necessary matlab functions to the search path
indfs = strfind(scriptDIR,filesep);
%analysDIR = scriptDIR(1:indfs(end-1)); % NS changed 20200206
analysDIR = scriptDIR(1:indfs(end));
addpath(genpath(analysDIR));
disp(analysDIR)
%%
% Set path for nifti images used by FSL and .MAT files
niiDIR = [analysDIR 'nifti'];
matDIR = [analysDIR 'sirdata'];

% Load SIR data
disp('Select SIR PAR/REC ... ')
[DATA_SIR,PAR_SIR] = loadPARREC;
DATA_SIR = squeeze(DATA_SIR);
DATA_SIR = DATA_SIR(:,:,:,:,1);
[nx,ny,nz,nti] = size(DATA_SIR);
SIR_PARNAME = PAR_SIR.filename
disp(analysDIR)

% Load B1 data
dataDIR = fileparts(SIR_PARNAME); cd(dataDIR);
disp('Select B1+ PAR/REC ... ')
[DATA_B1,PAR_B1] = loadPARREC;
MAP_B1 = squeeze(DATA_B1(:,:,:,:,:,:,1,3))/1e2;
DATA_B1 = squeeze(DATA_B1(:,:,:,:,:,:,1,1));
B1_PARNAME = PAR_B1.filename
%
% Input subject ID
subjectID = input('Enter subject ID: ','s');

% Set ti/td
ti = [10 10 278 1007]*1e-3;
td = [684 4121 2730 10]*1e-3;

% Inversion pulse dur
pw = 6.0e-3;

% Fixed kmf value
kmf_fixed = 14.5  % mean value in WM from NIMG paper

% MASK (very conservative for now)
MASK = DATA_SIR(:,:,:,end) > threshold3(DATA_SIR(:,:,:,end))/5;

% Denoise - window size from sims - large enough to yield stable results
disp('Denoising SIR data ... ')
[DATA_SIRf,m,p] = denoiseCV(DATA_SIR,[5 5 5],MASK); 

% And remove Gibbs 
disp('Deringing SIR data ... ')
DATA_SIRfg = unring(DATA_SIRf);

% Make SIR data NIFTIs for FSL
reconRes = [PAR_SIR.imgdef.pixel_spacing_x_y.uniq PAR_SIR.imgdef.slice_thickness_in_mm.uniq];
reconOrigin = [0 0 0];
DATA_SIR_NIFTI = permute(DATA_SIRfg,[2 1 3 4]);
for ii = 1:nti
    DATA_SIR_nii = make_nii(DATA_SIR_NIFTI(:,:,:,ii),reconRes,reconOrigin);
    SIR_fname{ii} = [niiDIR filesep subjectID '_TI_' num2str(ii) '.nii'];
    SIR_reg_fname{ii} = [niiDIR filesep subjectID '_TI_' num2str(ii) 'r.nii'];
    save_nii(DATA_SIR_nii,SIR_fname{ii});
end
%%
% Co-register SIR data via FSL
disp('Registering SIR data ... ')
REF = SIR_fname{end}
for ii = 1:nti 
    IN = SIR_fname{ii}
    OUT = SIR_reg_fname{ii}
    unix(['/usr/local/bin/flirt -paddingsize 1 -cost normmi -searchcost normmi -dof 6 -in ' IN ' -ref ' REF ' -out ' OUT ' -nosearch']); %NS 20200206
end

% Load SIR data into MATLAB 
DATA_SIR_REG = zeros(size(DATA_SIR));
for ii = 1:nti
    DATA_SIR_REG_nii = load_nii(SIR_reg_fname{ii});
    DATA_SIR_REG(:,:,:,ii) = DATA_SIR_REG_nii.img;  
end
DATA_SIR_REG = permute(DATA_SIR_REG,[2 1 3 4]);

% Create bias field corrected T1w images for brain extraction and segmentation
% disp('N4 bias field correction ... ')
% for ii = 1:nti
%     N4_fname{ii} = [niiDIR filesep subjectID '_N4_' num2str(ii) '.nii'];
%     unix(['N4BiasFieldCorrection -d 3 -i ' SIR_reg_fname{ii} ' -s 2 -c [100x100x100x100,0.0000000001] -o ' N4_fname{ii}]);
%     N4_nii = load_untouch_nii(N4_fname{ii});
%     DATA_SIR_REG_N4(:,:,:,ii) = double(permute(N4_nii.img,[2 1 3]));
% end
%%
disp('N4 bias field correction ... ')
for ii = 1:nti
    N4_fname{ii} = [niiDIR filesep subjectID '_N4_' num2str(ii) '.nii'];
    unix(['/Users/nicholassisco/Applications/ants/antsbin/ANTS-build/Examples/N4BiasFieldCorrection -d 3 -i ' SIR_reg_fname{ii} ...
        ' -s 2 -c [100x100x100x100,0.0000000001] -o ' N4_fname{ii}]);
    N4_nii = load_untouch_nii(N4_fname{ii});
    DATA_SIR_REG_N4(:,:,:,ii) = double(permute(N4_nii.img,[2 1 3]));
end
% Brain Extraction - second TI gives good results here
disp('Brain extraction ... ')
BET_fname = [niiDIR filesep subjectID '_BET.nii'];
unix(['bet2 ' N4_fname{2} ' ' BET_fname ' -f 0.1 -w 0.9']);  
BET_nii = load_untouch_nii(BET_fname);
BETMASK = double(permute(BET_nii.img,[2 1 3])) > 0;

% WM/GM segmentation - third TI after BET + N4 gives good results here
disp('FAST WM segmentation ... ')
FAST_fname_out = [niiDIR filesep subjectID '_FAST_out.nii'];

FAST_fname_in = [niiDIR filesep subjectID '_FAST_in.nii'];
DATA_FAST_nii = make_nii(permute(DATA_SIR_REG_N4(:,:,:,3).*BETMASK,[2 1 3]),reconRes,reconOrigin);
save_nii(DATA_FAST_nii,FAST_fname_in);

unix(['fast -N -n 3 -o ' FAST_fname_out ' ' FAST_fname_in]);     
SEG_nii = load_untouch_nii([FAST_fname_out(1:end-4) '_seg.nii']);
WMmask = double(permute(SEG_nii.img,[2 1 3])) == 2;

% And erode a few times to be conservative
WMmaske = erodeBWperim(WMmask,3);

% Identify connected region and get statistics for these areas
cc = bwconncomp(WMmaske, 8);
L = labelmatrix(cc);
R = regionprops(cc);

% Remove smaller size regions (these are likely not WM)
area = [R.Area];
indEx = find(area > 500);
WMmaskf = zeros(size(WMmask));
for jj = 1:length(indEx)
    WMmaskf = WMmaskf + (L == indEx(jj));
end

% Make B1 NIFTIs for FSL
disp('Register B1-SIR data ... ')
DATA_B1_NIFTI = permute(DATA_B1,[2 1 3]);
MAP_B1_NIFTI = permute(MAP_B1,[2 1 3]);

DATA_B1_nii = make_nii(DATA_B1_NIFTI,reconRes,reconOrigin);
MAP_B1_nii = make_nii(MAP_B1_NIFTI,reconRes,reconOrigin);

DATA_B1_fname = [niiDIR filesep subjectID '_B1_data.nii'];
DATA_reg_B1_fname = [niiDIR filesep subjectID '_B1_datar.nii'];

MAP_B1_fname = [niiDIR filesep subjectID '_B1_map.nii'];
MAP_reg_B1_fname = [niiDIR filesep subjectID '_B1_mapr.nii'];

save_nii(DATA_B1_nii,DATA_B1_fname)
save_nii(MAP_B1_nii,MAP_B1_fname)

% Find transform between B1 data and SIR
IN = DATA_B1_fname
OUT = DATA_reg_B1_fname
REF
OMAT = [niiDIR filesep subjectID '_B1_2_SIR']
unix(['/usr/local/bin/flirt -paddingsize 1 -cost normmi -searchcost normmi -dof 6 -in ' IN ' -ref ' REF ' -out ' OUT ' -nosearch -omat ' OMAT]); % NS changed 20200206
DATA_B1_nii = load_untouch_nii(OUT);
DATA_B1_REG = double(permute(DATA_B1_nii.img,[2 1 3]));
    
% Apply this transform to the B1 map
IN = MAP_B1_fname
OUT = MAP_reg_B1_fname
unix(['/usr/local/bin/flirt -paddingsize 1 -in ' IN ' -out ' OUT ' -ref ' REF ' -applyxfm -init ' OMAT]); % NS changed 20200206
MAP_B1_nii = load_untouch_nii(OUT);
MAP_B1_REG = double(permute(MAP_B1_nii.img,[2 1 3]));

% Calc absorption lineshape
delta = 0;
R2b = 1./10e-6;         % assume solid pool T2 is 10 us
lineShape = 'gaussian';
g = absorptionLineShape(1./R2b,delta,lineShape);

% Calculate Sm for each voxel
H = ones(10,10,3);
H = H./length(H(:));
MAP_B1s = imfilter(MAP_B1_REG,H,'replicate','conv');    % Smooth B1 to remove ringing/noise
Sm = zeros(nx,ny,nz);
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            if MASK(ii,jj,kk) == 1
                Sm(ii,jj,kk) = calcSm_compPulse(MAP_B1s(ii,jj,kk),g,pw);
            end
        end
    end
end


% Extract vectors of non-masked signals and Sm for parallel implementation 
DATAn = DATA_SIR_REG(:,:,:,end);
MASKfit = MASK.*BETMASK.*(prod(DATA_SIR_REG,4) > 0);
ind = find(MASKfit); 
%%
for ii = 1:nti
    DATAc = DATA_SIR_REG(:,:,:,ii);
    sigv(:,ii) = DATAc(ind) ./ DATAn(ind) ;
end
%%
Smv = Sm(ind);

% Create pool
% tic, warning off
% delete(gcp('nocreate'))
% pool = parpool(6);

% Parallel loop
X0 = [0.1  1   -0.95  1];  % [pmf R1f Sf M0f]
LB = [0    0.3 -1.05  0];
UB = [1    3    0    10];
Xv = zeros(length(ind),4);
tic
parfor ii = 1:length(ind)
%    disp(['Fitting SIR data: Voxel: ' num2str(ii) '/' num2str(length(ind))])
    Xv(ii,:) = fitSIR_fixedkmf(ti,td,sigv(ii,:),Smv(ii),X0,LB,UB,'y','y',kmf_fixed);
end
toc
% % Close pool
% delete(gcp('nocreate'))
% toc, warning on

% Reorder parameters into image space
PSR = zeros(nx,ny,nz); 
PSR(ind) = Xv(:,1)*1e2;
R1f = zeros(nx,ny,nz); 
R1f(ind) = Xv(:,2);
Sf = zeros(nx,ny,nz);  
Sf(ind) = Xv(:,3);

%%
% Display final results
figure, imagesc(myMontage(PSR(:,:,8:4:end)),[0 30]), axis image off, colormap(gray(256)), colorbar
figure, imagesc(myMontage(R1f(:,:,8:4:end)),[0.25 1]), axis image off, colormap(gray(256)), colorbar


% Save results
save(['SIRData_par_' num2str(subjectID)])




