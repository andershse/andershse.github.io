function max_size = findMaxIndependentSet(polys, vars)
    % polys: vektor av polynom i bildens ideal, t.ex. [y3]
    % vars: målvariabler [y1, y2, y3]
    
    n = length(vars);
    varNames = arrayfun(@char, vars, 'UniformOutput', false);
    max_size = 0;
    
    % Support of each leading monomial in the lexicographic Gröbner basis.
    % A leading monomial can contain several variables (for example y1*y2).
    ltVarGroups = cell(1, length(polys));
    for j = 1:length(polys)
        p = expand(polys(j));
        support = false(1,n);
        for v = 1:n
            [c, monomials] = coeffs(p, vars(v));
            degrees = arrayfun(@(m) double(polynomialDegree(m,vars(v))), monomials);
            [degree, idx] = max(degrees);
            support(v) = degree > 0;
            p = c(idx);
        end
        if ~any(support)
            error('The leading-term ideal contains 1; its variety is empty.');
        end
        ltVarGroups{j} = varNames(support);
    end
    
    % 2. Uddins kombinatoriska sökning
    for i = 0:(2^n - 1)
        subsetIdx = bitget(i, 1:n) == 1;
        currentSubset = varNames(subsetIdx);
        
        isIndependent = true;
        for j = 1:length(ltVarGroups)
            varsInLT = ltVarGroups{j};
            if isempty(varsInLT), continue; end
            
            % Om den ledande variabeln finns i vår delmängd är den BEROENDE
            if all(ismember(varsInLT, currentSubset))
                isIndependent = false;
                break;
            end
        end
        
        if isIndependent
            max_size = max(max_size, sum(subsetIdx));
        end
    end
end





