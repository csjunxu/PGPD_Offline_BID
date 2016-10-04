function      [im_out, par] = PGs2Image(X_hat,W,par)
% Reconstruction
im_out = zeros(par.h,par.w,'double');
im_wei = zeros(par.h,par.w,'double');
r = 1:par.maxr;
c = 1:par.maxc;
k = 0;
for i = 1:par.ps
    for j = 1:par.ps
        k = k+1;
        im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( X_hat(k,:)', [par.maxr par.maxc]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [par.maxr par.maxc]);
    end
end
im_out  =  im_out./im_wei;