%  This code was written by Tomohiko Mizutani (April 6, 2023)
 
function K = rpp(A, r, ptLst)


[d, n] = size(A);

ordDatTab = [];
disTab = [];

for i=1:n
    
    col = A(:,i);
    [disLst, ordDatLst] = sort(sum(abs(A - repmat(col, 1, n))));
    
    disTab(i, :) = disLst;
    ordDatTab(i, :) = ordDatLst;
    
end


clst = struct('elem', [], 'pt', [], 'diam', []);
[~, L] = sort(ptLst, 'descend');



for j=1:r

    minDiam = inf;
    elemLst = [];

    
    for i=1:n
        
        idx = L(i);
        
	ordDatLst = ordDatTab(idx, :);
	disLst = disTab(idx, :);
    
	sumLst = cumsum(ptLst(ordDatLst));
	
	I = [1:n];
	I(sumLst <= r / (r+1)) = [];

	if isempty(I) == 0
	    
	    headIdx = min(I);
	    diam = disLst(1, headIdx);
	    
	    if diam < minDiam

		elemLst = ordDatLst([1:headIdx]);
		minDiam = diam;
		
	    end	    
	end
    end

    clst(j).elem = elemLst;
    clst(j).pt = ptLst(elemLst);
    clst(j).diam = minDiam;
    
    ptLst(elemLst) = 0;    
    
end


K = getIdxSet_maxPt(clst, r);


end


%%%%%%%%
function [K] = getIdxSet_maxPt(clst, r)

for j=1:r

    
    elemLst = clst(j).elem;
    elemPtLst = clst(j).pt;
    
    if isempty(elemLst) == 0
	    
	[~, idx] = max(clst(j).pt);
	K(j) = clst(j).elem(idx);
	
    else
	
	K(j) = 0;
	
    end
    
end


end







