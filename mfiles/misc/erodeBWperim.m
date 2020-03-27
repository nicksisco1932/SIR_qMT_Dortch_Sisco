function MASKe = erodeBWperim(MASK,N)

% Get number of slices
nz = size(MASK,3);

% Loop over number of pixels to remove from edges
for jj = 1:N
    
    % Loop over slices
    for ii = 1:nz
        B(:,:,ii) = bwperim(MASK(:,:,ii));
    end
    MASKe = MASK.*~B;
    
end