%% Center the displacements that are missing an initial data point 

function centeredu = CenterDisps(medsize, medsizenan, u, GPSNameList, ntime)
    % Make the displacements start at 0 by taking the 20 closest data points to
    % the start
    interpu = zeros(3, length(GPSNameList));
    for i = 1:length(GPSNameList)
        unfilledu = squeeze(u(i, :, :));
        nanu = isnan(unfilledu(:, :));
        nanu = max(nanu, [], 1);
        notnanu = ~isnan(unfilledu);

        if(isnan(unfilledu(1, 1))); unfilledu(:, nanu) = [];
            interpu(:, i) = median(unfilledu(:, 1:medsizenan), 2, 'omitmissing');
        else; unfilledu(:, nanu) = [];
            interpu(:, i) = median(unfilledu(:, 1:medsize), 2, 'omitmissing'); 
        end
        % if(GPSNameList(i) == "V120")
        %     t = 1:ntime;
        %     for j = 1:3
        %         fitobj = fit(t(notnanu(j, :))', squeeze(u(i, j, notnanu(j, :))), 'exp2');
        %         interpu(j, i) = fitobj.a + fitobj.c;%[-0.05; 0.038; 0.038];
        %     end
        % end
    end
    
    for i = 1:ntime
        u(:, :, i) = u(:, :, i) - interpu';
    end
    
    centeredu = u;
end