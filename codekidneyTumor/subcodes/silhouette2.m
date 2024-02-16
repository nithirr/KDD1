function [s,h] = silhouette2(X, clust, distance, varargin)

if nargin < 2
    error('stats:silhouette:TooFewInputs',...
          'At least two input arguments required.');
end


[idx,cnames] = grp2idx(clust);
nIn = length(idx);
if ~isempty(X) &  (nIn ~= size(X,1))
    error('stats:silhouette:InputSizeMismatch',...
          'The number of rows in X must match the length of CLUST.');
end

if ~isempty(X)
    nans = find(isnan(idx) & any(isnan(X),2));
    if length(nans) > 0
        X(nans,:) = [];
        idx(nans) = [];
    end
else
    nans = find(isnan(idx));
    if length(nans) > 0
        idx(nans) = [];
    end
end
n = length(idx);
p = size(X,2);
k = length(cnames);
count = histc(idx(:)',1:k);
simatrix = 0;  % my adding to accept matrix

if nargin < 3 | isempty(distance)
    distType = 'sqeuclidean';
else
    if ischar(distance)
        distNames = {'euclidean','sqeuclidean','cityblock','cosine',...
                     'correlation','hamming','jaccard'};
        i = strmatch(lower(distance), distNames);
        if length(i) > 1
            error('stats:silhouette:BadDistance',...
                  'Ambiguous ''distance'' argument:  %s.', distance);
        elseif isempty(i)
           
            distType = 'function'; % user-defined distance function name
            distance = str2func(distance);
        else
            distType = distNames{i};
        end

        switch distType
        case 'cosine'
            Xnorm = sqrt(sum(X.^2, 2));
            if any(min(Xnorm) <= eps(max(Xnorm)))
                error('stats:silhouette:InappropriateDistance',...
              ['Some points have small relative magnitudes, making them ', ...
               'effectively zero.\nEither remove those points, or choose ', ...
               'a distance other than cosine.']);
            end
            X = X ./ Xnorm(:,ones(1,p));
        case 'correlation'
            X = X - repmat(mean(X,2),1,p);
            Xnorm = sqrt(sum(X.^2, 2));
            if any(min(Xnorm) <= eps(max(Xnorm)))
                error('stats:silhouette:InappropriateDistance',...
           ['Some points have small relative standard deviations, making ', ...
            'them effectively constant.\nEither remove those points, or ', ...
            'choose a distance other than correlation.']);
            end
            X = X ./ Xnorm(:,ones(1,p));
        end
    elseif isnumeric(distance)
        if (size(distance,1) == 1) & (size(distance,2) == .5*nIn*(nIn-1))
            if any(distance < 0)
                error('stats:silhouette:BadDistanceMatrix',...
                   'A distance matrix must contain only non-negative values.');
            end
            distType = 'matrix'; 
            if length(nans) > 0
                distance = UTMatSub(distance, nans);
            end
        elseif size(distance,1) == size(distance,2)
            simatrix = 1;
            distType = 'matrixfull'; 
        else
            error('stats:silhouette:BadDistanceMatrix',...
    'A distance matrix must be in upper triangular form as created by PDIST.');
        end
    elseif isa(distance,'function_handle') | isa(distance,'inline')
        distType = 'function'; % user-defined distance function
    else
        error('stats:silhouette:BadDistance',...
              ['The ''distance'' argument must be a string, ' ...
               'a numeric matrix, or a function.']);
    end
    distArgs = varargin(1:end); 
end

mbrs = (repmat(1:k,n,1) == repmat(idx,1,k));
myinf = zeros(1,1,class(X));
myinf(1) = Inf;
avgDWithin = repmat(myinf, n, 1);
avgDBetween = repmat(myinf, n, k);
for j = 1:n
    switch distType
    case 'euclidean'
        distj = sqrt(sum((X - X(repmat(j,n,1),:)).^2, 2));
    case 'sqeuclidean'
        distj = sum((X - X(repmat(j,n,1),:)).^2, 2);
    case 'cityblock'
        distj = sum(abs(X - X(repmat(j,n,1),:)), 2);
    case {'cosine','correlation'}
        distj = 1 - (X * X(j,:)');
    case 'hamming'
        distj = sum(X ~= X(repmat(j,n,1),:), 2) / p;
    case 'jaccard'
        Xj = X(repmat(j,n,1),:);
        nz = X ~= 0 | Xj ~= 0;
        ne = X ~= Xj;
        distj = sum(ne & nz, 2) ./ sum(nz, 2);
    case 'matrixfull'
        distj = distance(:, j);
    case 'matrix'
        distj = UTMatCol(distance, j);
    case 'function'
        try
            distj = feval(distance, X(j,:), X, distArgs{:});
        catch
            [errMsg,errID] = lasterr;
            if isa(distance,'inline')
                error('stats:silhouette:DistanceFunctionError',...
                  ['The inline distance function generated the following ', ...
                   'error:\n%s'], lasterr);
            elseif strcmp('MATLAB:UndefinedFunction', errID) ...
                        && ~isempty(strfind(errMsg, func2str(distance)))
                error('stats:silhouette:DistanceFunctionNotFound',...
                      'The distance function ''%s'' was not found.', ...
                      func2str(distance));
            else
                error('stats:silhouette:DistanceFunctionError',...
                  ['The distance function ''%s'' generated the following ', ...
                   'error:\n%s'], func2str(distance),lasterr);
            end
        end;
end
    

        R = isnan(distj);
        distj(R) = 0;
    for i = 1:k
        Q = R(mbrs(:,i));
        L = sum(Q);   
        if i == idx(j)
            avgDWithin(j) = sum(distj(mbrs(:,i))) ./ max(count(i)-1-L, 1);
        else
            avgDBetween(j,i) = sum(distj(mbrs(:,i))) ./ max(count(i)-L, 1);
        end
    end
end


minavgDBetween = min(avgDBetween, [], 2);
silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);

if (nargout == 0) | (nargout > 1)
   
    space = max(floor(.02*n), 2);
    bars = repmat(NaN,space,1);
    for i = 1:k
        bars = [bars; -sort(-silh(idx == i)); repmat(NaN,space,1)];
        tcks(i) = length(bars);
    end
    tcks = tcks - 0.5*(diff([space tcks]) + space - 1);
    
  
    if k > 20
        cnames = '';
    end
    barsh = barh(bars, 1.0);
    axesh = get(barsh(1), 'Parent');
    set(axesh, 'Xlim',[-Inf 1.1], 'Ylim',[1 length(bars)], 'YDir','reverse', 'YTick',tcks, 'YTickLabel',cnames);
    if n > 50
        shading flat
    end
    xlabel('Silhouette Value');
    ylabel('Cluster');
end

if nargout > 0
    s = silh;
end
if nargout > 1
    h = get(axesh, 'Parent');
end

%------------------------------------------------------------------

function Aj = UTMatCol(A, j)

n = (1+sqrt(1+8*length(A)))/2;


Aj = 0;

if j > 1
    ii = 1:(j-1);
    Aj = [A((ii-1).*(n-ii/2)+j-ii)'; Aj];
end

if j < n
    jj = (j+1):n;
    Aj = [Aj; A((j-1).*(n-j/2)+jj-j)'];
end


%------------------------------------------------------------------

function A = UTMatSubmat(A, cut)


n = (1+sqrt(1+8*length(A)))/2;

dels = [];
for j = cut
 
    if j > 1
        ii = 1:(j-1);
        dels = [dels (ii-1).*(n-ii/2)+j-ii];
    end
   
    if j < n
        jj = (j+1):n;
        dels = [dels (j-1).*(n-j/2)+jj-j];
    end
end
A(dels) = [];
