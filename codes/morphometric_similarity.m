function similarity_matrix = morphometric_similarity(data, nbins, use_mixture, similarity_transform)

    if nargin < 3 || isempty(use_mixture)
        use_mixture = 0; 
    end

    if nargin < 4 || isempty(similarity_transform)
        similarity_transform = 'one_minus_jsd';
    end
    assert(any(strcmp(similarity_transform,{'one_minus_jsd','one_minus_sqrt_jsd'})), ...
        'Invalid similarity_transform value');

    all_values = vertcat(data{:});
    all_values = all_values(isfinite(all_values));
    if isempty(all_values)
        error('All input values are NaN/Inf.');
    end


    if use_mixture == 1

        all_pos = all_values(all_values > 0);
        if isempty(all_pos)
            error('All values are zero; mixture KDE on nonzero part is undefined.');
        end

        s  = std(all_pos, 0, 'omitnan');
        iq = iqr(all_pos);
        sigma = min(s, iq/1.34);
        h = 0.9 * sigma * numel(all_pos)^(-1/5);
        if ~(isfinite(h) && h > 0)
            h = 1e-3 + 0.1 * max(eps, s);
        end

        min_pos = min(all_pos);
        max_pos = max(all_pos);
        x_lower = max(0, min_pos - 3*h);
        x_upper = max_pos + 3*h;

        if nargin < 2 || isempty(nbins)
            dx = h/4;
            nbins = ceil((x_upper - x_lower) / dx);
            nbins = max(128, min(2048, nbins));
        end

        x_range = linspace(x_lower, x_upper, nbins);

        R = numel(data);

        PDFs = zeros(R, nbins + 1);

        for r = 1:R
            values = data{r}(:);
            values = values(isfinite(values));
            if isempty(values)
                PDFs(r,:) = 1/(nbins+1); 
                continue;
            end

            nr = numel(values);
            p0 = sum(values == 0) / nr;
            pos = values(values > 0);

            if isempty(pos)

                P = zeros(1, nbins + 1);
                P(1) = 1;
                PDFs(r,:) = P;
                continue;
            end

            f = ksdensity(pos, x_range, 'Bandwidth', h);
            f = max(f, 0);
            area = trapz(x_range, f);
            if ~(isfinite(area) && area > 0)
                f = ones(1, nbins) / nbins;
            else
                f = f / area;
            end

            P = [p0, (1 - p0) * f];
            P = P / sum(P); 
            PDFs(r,:) = P;
        end

        R = size(PDFs, 1);
        JSD_matrix = zeros(R, R);
        eps_val = 1e-12;

        for i = 1:R
            Pi = PDFs(i,:) + eps_val;
            Pi = Pi / sum(Pi);
            for j = i+1:R
                Qj = PDFs(j,:) + eps_val;
                Qj = Qj / sum(Qj);
                M  = 0.5 * (Pi + Qj);
                KL_PM = sum(Pi .* log2(Pi ./ M));
                KL_QM = sum(Qj .* log2(Qj ./ M));
                JSD = 0.5 * (KL_PM + KL_QM);

                JSD = max(0, min(1, JSD));
                JSD_matrix(i,j) = JSD;
                JSD_matrix(j,i) = JSD;
            end
        end

        D = sqrt(JSD_matrix);
        similarity_matrix = 1 - D;
        similarity_matrix = min(1, max(0, similarity_matrix));

        return
    end



    n = numel(all_values);
    s  = std(all_values, 0, 'omitnan');
    iq = iqr(all_values);
    sigma = min(s, iq/1.34);
    h = 0.9 * sigma * n^(-1/5);
    if ~(isfinite(h) && h > 0)
        h = 1e-3 + 0.1 * max(eps, s);
    end

    min_val = min(all_values);
    max_val = max(all_values);
    x_lower = min_val - 3*h;
    x_upper = max_val + 3*h;

    if nargin < 2 || isempty(nbins)
        dx = h/4;
        nbins = ceil((x_upper - x_lower) / dx);
        nbins = max(256, min(2048, nbins));
    end

    x_range = linspace(x_lower, x_upper, nbins);

    R = numel(data);
    PDFs = zeros(R, nbins);
    for r = 1:R
        values = data{r}(:);
        values = values(isfinite(values));
        if isempty(values)
            PDFs(r,:) = 1/nbins;
            continue;
        end
        f = ksdensity(values, x_range, 'Bandwidth', h);
        f = max(f, 0);
        area = trapz(x_range, f);
        if ~(isfinite(area) && area > 0)
            PDFs(r,:) = 1/nbins;
        else
            PDFs(r,:) = f / area;
        end
    end

    JSD_matrix = zeros(R, R);
    eps_val = 1e-12;
    for i = 1:R
        Pi = PDFs(i,:) + eps_val;
        Pi = Pi / sum(Pi);
        for j = i+1:R
            Qj = PDFs(j,:) + eps_val;
            Qj = Qj / sum(Qj);
            M  = 0.5 * (Pi + Qj);
            KL_PM = sum(Pi .* log2(Pi ./ M));
            KL_QM = sum(Qj .* log2(Qj ./ M));
            JSD = 0.5 * (KL_PM + KL_QM);
            JSD_matrix(i,j) = JSD;
            JSD_matrix(j,i) = JSD;
        end
    end

    % ---- Only altered step in this experiment: convert EXACT SAME JSD ----
    % Original: 1 - JSD; sensitivity analysis: 1 - sqrt(JSD).
    % NOTE: do NOT set use_mixture=1: it changes the underlying KDE model!
    JSD_matrix = max(0,min(1,JSD_matrix));
    if strcmp(similarity_transform,'one_minus_sqrt_jsd')
        similarity_matrix = 1 - sqrt(JSD_matrix);
    else
        similarity_matrix = 1 - JSD_matrix;
    end
    similarity_matrix = min(1, max(0, similarity_matrix));
end
