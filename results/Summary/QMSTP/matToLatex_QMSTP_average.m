function matToLatex_QMSTP_average()
filename = 'TableTex_QMSTP_ave.txt';
d0 = load("QMSTP-ALM.mat");
d3 = load("QMSTP-NAL.mat");
d0 = d0.rrALM;
d3 = d3.rrNAL;

for i = 1:size(d0,1)-1
    if iscell(d0{i+1,1})
        d0{i+1,1} = d0{i+1,1}{1};
    end
end
for i = 1:size(d3,1)-1
    if iscell(d3{i+1,1})
        d3{i+1,1} = d3{i+1,1}{1};
    end
end

fileID = fopen(filename, 'w');
fclose(fileID);
fileID = fopen(filename, 'a');

%% take average

ave_NAL = d3(1,:);
ave_Rie = d0(1,:);
instance_type = {'vsym','sym','esym'};
instance_dim = [435 1225];
count = 0;
for i = 1:3
    for j = 1:2
        count = count + 1;
        name = strcat(['qmstp_',instance_type{i}]);
        dim  = instance_dim(j);
        ave_NAL{count+1,1} = instance_type{i};
        ave_Rie{count+1,1} = instance_type{i};
        %% SDPNAL+
        rowsWithName = cellfun(@(x) ischar(x) && startsWith(x, name), d3(:, 1)) & ...
               cellfun(@(x) isnumeric(x) && x == dim, d3(:, 2));
        matchingRows = d3(rowsWithName, :);
        if size(matchingRows,1) == 10 
            averow = mean(cell2mat(matchingRows(:,2:end)),1);
            for k = 1:size(averow,2)
                ave_NAL{count+1,1+k} = averow(k);
            end
        else
            for k = 1:size(averow,2)
                ave_NAL{count+1,1+k} = NaN;
            end
            ave_Rie{count+1,2} = dim;
        end
        %% RiNNAL+
        rowsWithName = cellfun(@(x) ischar(x) && startsWith(x, name), d0(:, 1)) & ...
               cellfun(@(x) isnumeric(x) && x == dim, d0(:, 2));
        matchingRows = d0(rowsWithName, :);
        if size(matchingRows,1) == 10 
            averow = mean(cell2mat(matchingRows(:,2:end)),1);
            for k = 1:size(averow,2)
                ave_Rie{count+1,1+k} = averow(k);
            end
        else
            for k = 1:size(averow,2)
                ave_Rie{count+1,1+k} = NaN;
            end
            ave_Rie{count+1,2} = dim;
        end
    end
end

d0 = ave_Rie;
d3 = ave_NAL;

count = 0;

for i = [1 3 5 2 4]
        %%
        dd = d0(i+1,:);
        fprintf(fileID,'%s',dd{1});
        fprintf(fileID, '& RiNNAL+ ');
        if ~isnan(dd{14})
            if dd{14}<3600
            fprintf(fileID, '& %.0f, %.0f, %.0f & %.0f ',dd{12},dd{13},dd{17},dd{5});
            fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
            fprintf(fileID, '& %.7e ',dd{7});
            fprintf('\n RiNNAL+: %.2f',dd{7});
            fprintf(fileID, '& %.1f &%.1f \\\\ \n',dd{14},dd{16});
            else
            fprintf(fileID, '& %.0f, %.0f, %.0f & %.0f ',dd{12},dd{13},dd{17},dd{5});
            fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
            fprintf(fileID, '& %.7e ',dd{7});
            fprintf('\n RiNNAL+: %.2f',dd{7});
            fprintf(fileID, '& %.1f &%.1f \\\\ \n',3600.0,dd{16});
            end
        else
            fprintf(fileID, '& - & - & - & - & - & - & - \\\\ \n');
        end
        %%
        Index = i+1;
        fprintf(fileID,'$n=%d$ ',dd{2});
        fprintf(fileID, '& SDPNAL+ '); 
        flag = 0;
        if Index ~= 0
            if d3{Index,15}<3600 
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %.0f, %.0f, %.0f & %.0f',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf('\n SDPNAL+: %.2f',dd{7});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',dd{15});
            else
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %.0f, %.0f, %.0f & %.0f',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf('\n SDPNAL+: %.2f',dd{7});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',min(3600.0,dd{15}));
            end
        end
        if flag == 0 
            fprintf(fileID, '& - & - & - & - & - & - \\\\ [3pt] \n\n');
        end
end
fclose(fileID);
