function funval = CostgBlur(gB,data,model_stack,flag,std_Gbg,thr)

h = (fspecial('gaussian',[5 5],gB));
for j = 1:size(data,3)
    tmp = imfilter(model_stack(:,:,j),h,'replicate');
    tmp_data =  data(:,:,j);
    
    % thr mask
    %     mask = tmp_data>max(tmp_data(:)).*thr(j);
    mask = tmp_data>thr(j);
    tmp = tmp.*mask;
    tmp_data = tmp_data.*mask;
    
    if flag==1
        % L1 diff
        Diff(j) = sum(abs(tmp(:)-tmp_data(:)));
    elseif flag==2
        %     L2 diff
        Diff(j) = sum((tmp(:)-tmp_data(:)).^2);
    elseif flag==3
        %     P-MLE diff
        Diff(j) = -sum(tmp_data(:).*log(tmp(:)+eps) - tmp(:));
    elseif flag==4
        % Gaussian MLE
        Diff(j) = sum(mask(:).*(log(sqrt(2*pi*(tmp(:)+std_Gbg(j).^2+eps)))+0.5*((tmp_data(:)-tmp(:)).^2)./(tmp(:)+std_Gbg(j).^2+eps)));
    elseif flag==5
        %     corr2 cost
        Diff(j) = -corr2(tmp./norm(tmp(:)),tmp_data./norm(tmp_data(:)));
    elseif flag ==6 % cross-entropy
        imn = @(im) (im - min(im(:)))./(max(im(:)) - min(im(:)));
        tmp = imn(tmp)+eps;
        tmp_data = imn(tmp_data)+eps;
        Diff(j) = - sum(real(tmp_data(:).*log(tmp(:))+(1-tmp_data(:)).*log(1-tmp(:))));
    end
    %     figure(909)
    %     subplot(1,2,1)
    %     imagesc(tmp(round(end/2-50):round(end/2+50),round(end/2-50):round(end/2+50)))
    %     colorbar
    %     subplot(1,2,2)
    %     imagesc(tmp_data(round(end/2-50):round(end/2+50),round(end/2-50):round(end/2+50)))
    %     colorbar
    %     drawnow
end

%% L2
funval = mean(Diff(:));
end

