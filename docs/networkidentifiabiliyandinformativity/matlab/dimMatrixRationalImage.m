function d = dimMatrixRationalImage(P_mat, Q_mat, vars)
    [rows, cols] = size(P_mat);
    num_elements = rows * cols;
    
    % Definiera målvariabler y1, y2, y3...
    y_flat = sym('y', [1, num_elements]);
    targetVars = reshape(y_flat, [rows, cols]);
    
    syms w; 
    mappingEqs = Q_mat .* targetVars - P_mat;
    
    % Saturation (nämnare)
    uniqueDenoms = unique(Q_mat);
    uniqueDenoms = uniqueDenoms(~arrayfun(@(x) isempty(symvar(x)), uniqueDenoms));
    if isempty(uniqueDenoms), satEq = []; else satEq = w * prod(uniqueDenoms) - 1; end
    
    % Ideal och elimineringsordning: [vars, w] > [y_flat]
    ideal = [mappingEqs(:).', satEq];
    allVars = [vars, w, y_flat]; 
    
    fprintf('\n--- Beräknar Gröbnerbas ---\n');
    gb = gbasis(ideal, allVars, 'MonomialOrder', 'lex');
    if numel(gb) == 1 && isequal(gb(1), sym(1))
        error('The rational map has an empty admissible domain.');
    end
    
    % Filtrera fram bildens ideal (relationer i y-variablerna)
    imageGB = [];
    for i = 1:length(gb)
        v_in_p = symvar(gb(i));
        % Ett polynom tillhör bilden om det INTE innehåller källvariabler
        if ~any(ismember(v_in_p, [vars, w])) && ~isempty(v_in_p)
            imageGB = [imageGB, gb(i)];
        end
    end
    
    fprintf('Identifierade ekvationer för bilden: '); disp(imageGB);

    if isempty(imageGB)
        d = num_elements;
    else
        % Vi skickar imageGB direkt till sökningen. 
        % findMaxIndependentSet kommer självt att extrahera ledande termer.
        d = findMaxIndependentSet(imageGB, y_flat);
    end
end


