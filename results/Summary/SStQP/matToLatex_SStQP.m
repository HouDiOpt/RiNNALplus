function matToLatex_SStQP()
filename = 'TableTex_SStQP.txt';
d0 = load("SStQP-ALM.mat");
d3 = load("SStQP-NAL.mat");

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

count = 0;

for i = [1:8]
        %%
        dd = d0(i+1,:);
        parts = split(dd{1}, '_');
        fprintf(fileID,'%s',parts{1});
        fprintf(fileID, '& RiNNAL+ ');
        if dd{14}<3600
            fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{18},dd{5});
            fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
            fprintf(fileID, '& %.7e ',dd{6});
            fprintf(fileID, '& %.1f &%.1f \\\\ \n',dd{14},dd{17});
        else
            fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{18},dd{5});
            fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
            fprintf(fileID, '& %.7e ',dd{6});
            fprintf(fileID, '& %.1f &%.1f \\\\ \n',3600.0,dd{17});
        end
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
                fprintf('\n %s: n = %d, ratio = %.1d',dd{1},dd{2},dd{15}/d0{i+1,14});
            else
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{6});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',3600.0);
                fprintf('\n %s: n = %d, ratio = %.1d',dd{1},dd{2},dd{15}/d0{i+1,14});
            end
        end
        if flag == 0 
            fprintf('\n %s: n = %d, ratio = %.1d',dd{1},dd{2},3600/d0{i+1,14});
            fprintf(fileID, '& - & - & - & - & - & - \\\\ [3pt] \n\n');
        end
end
fclose(fileID);
