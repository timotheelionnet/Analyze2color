function [Dstd,N] = std_nonzeros(D)
Dstd = zeros(size(D,1),1);
N = zeros(size(D,1),1);
for i=1:size(D,1)
    Dcur = D(i,:);
    Dstd(i) = nanstd(Dcur(Dcur~=0));
    N(i) = sum((Dcur~=0).*(~isnan(Dcur)));
end

end