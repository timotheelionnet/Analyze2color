function Dmed = median_nonzeros(D)
Dmed = zeros(size(D,1),1);
N = zeros(size(D,1),1);
for i=1:size(D,1)
    Dcur = D(i,:);
    Dmed(i) = nanmedian(Dcur(Dcur~=0));
    N(i) = sum(Dcur~=0);
end

end