function matToLatex_ccMSSC()
filename = 'TableTex_ccMSSC.txt';
d0 = load("ccMSSC-ALM.mat");
d3 = load("ccMSSC-NAL.mat");
d0 = d0.rrALM;
d3 = d3.table_NAL;

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

count = 0;
[~,sortid] = sort(cell2mat(d0(2:end,2)));
modfre = 1;%5

for k = 1:length(d0)-1
        i = sortid(k);
        %%
        dd = d0(i+1,:);
        if dd{14}<=3600 && dd{2}>=500 && max([dd{8},dd{9},dd{10}])<=1e-6
        count = count + 1;
        if mod(count,modfre) == 0
        fprintf(fileID, '%d',count);
        fprintf(fileID, '& RiNNAL+ ');
        fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{19},dd{5});
        fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
        fprintf(fileID, '& %.7e ',dd{6});
        fprintf(fileID, '& %.1f &%.1f \\\\ \n',dd{14},dd{18});
        %%
        Index = find(strcmp([d3{:,1}], dd{1}));
        fprintf(fileID,'%d, %d',dd{16},dd{2});
        fprintf(fileID, '& SDPNAL+ '); 
        flag = 0;
        if Index ~= 0
            if d3{Index,15}<3600
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{6});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',dd{15});
                fprintf('\n %d: n = %d, ratio = %.1d',count,dd{2},dd{15}/d0{i+1,14});
            else
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{6});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',3600.0);
                fprintf('\n %d: n = %d, ratio = %.1d',count,dd{2},dd{15}/d0{i+1,14});
            end
        end
        if flag == 0 
            fprintf(fileID, '& - & - & - & - & - & - \\\\ [3pt] \n\n');
        end
        end
        elseif dd{2}>=500
            fprintf('\n cannot solve %s',dd{1});
        end
end
fclose(fileID);
