function movMontage(IMG,N,t,cmap,xx,yy,cp)

% Play movie
IMG = squeeze(IMG); 
for ii = 1:N
    for jj = 1:size(IMG,3)
        imagesc(IMG(:,:,jj)), axis image off, colormap(cmap)
        title(['Image:' num2str(jj) '; Loop:' num2str(ii)])
        if nargin > 4
            hold on, plot(xx,yy,cp,'LineWidth',2)
        end
        pause(t), hold off
    end
end