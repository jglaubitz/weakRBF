function DM = Tools_DistanceMatrix(EvalPoints, RefPoints)
    % ------------------------------------------------------
    %  Computing the distance of the points in EvalPoints to
    %  all the points in RefPoints
    %
    % Input: 
    %        EvalPoints:  (NbPoints x 2)    array,   List of target points x-y coordinates 
    %        RefPoints:   (NbRefPoints x 2) array,   List of reference points x-y coordinates 
    %
    % Output:
    %        DM: (NbPoints x NbRefPoints) array,  Distance matrix between each of the points
    %
    % ------------------------------------------------------
    
    % Initialisation of the matrix upon the target-refence points numbers
    N = size(RefPoints,1);
    K = size(EvalPoints,1);
    DM = zeros(K,N);

    % Determining the distance matrix (version 1)
    %for n = 1:K 
    %    DM(n,:) = vecnorm((EvalPoints(n,:)-RefPoints)', 2);
    %end
    
    % Determining the distance matrix (version 2)
    DM = sqrt((EvalPoints(:,1)-RefPoints(:,1)').^2+(EvalPoints(:,2)-RefPoints(:,2)').^2);
    
end
