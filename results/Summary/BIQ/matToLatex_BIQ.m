function matToLatex_BIQ()
filename = 'TableTex_BIQ.txt';
d0 = load("BIQ-ALM.mat");
d3 = load("BIQ-NAL.mat");
d0 = d0.rrALM;
d3 = d3.rrNAL_slack_comp;

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

for i = 1:10:32%1:length(d0)-1
        %%
        dd = d0(i+1,:);
        if mod(i,10) ~= 0
            fprintf(fileID, '%d',mod(i,10));
        else
            fprintf(fileID, '%d',10);
        end
        fprintf(fileID, '& RiNNAL+ ');
        fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{17},dd{5});
        fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
        fprintf(fileID, '& %.7e ',dd{6});
        fprintf(fileID, '& %.1f &%.1f \\\\ \n',dd{14},dd{16});
        %%
        Index = find(strcmp([d3{:,1}], dd{1}));
        fprintf(fileID,'$n=%d$ ',dd{2});
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
            else
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{6});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',3600.0);
            end
        end
        if flag == 0 
            fprintf(fileID, '& - & - & - & - & - & - \\\\ [3pt] \n\n');
        end
end
fclose(fileID);
