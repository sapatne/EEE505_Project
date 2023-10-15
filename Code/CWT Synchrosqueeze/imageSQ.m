function imageSQ(t, ytic, M) 
% IMAGESQ displays time-frequency result of the Synchrosqueezing transform
% and others functions

fz = 20;

S = size(M);
Q = M(:);
q = quantile(Q, 0.998);
M(find(M>q)) = q;

imagesc(t, ytic, M)
axis xy ;
set(gca, 'fontsize', fz);
