function data_new = interpolate_data(xx,yy,data_sts,t)
    [nx,ny] = size(xx);
    data_new = zeros(nx,ny,46);
    data_new(:,:,1) = xx;
    data_new(:,:,2) = yy;
    ti = t(1);
    tf = t(2);
    [nel,gll] = size(data_sts(:,:,1));
    xq = reshape(data_sts(:,:,1), nel*gll,1);
    yq = reshape(data_sts(:,:,2), nel*gll,1);
    method = 'nearest';
    for i=3:46
        vq = reshape(data_sts(:,:,i), nel*gll,1);
        data_new(:,:,i) = griddata(xq,yq,vq,xx,yy)*(tf-ti);
    end
    
end
