function matToLatex_QMSTP()
filename = 'TableTex_QMSTP_all.txt';
d0 = load("QMSTP-ALM.mat");
d3 = load("QMSTP-NAL.mat");
d0 = d0.rrALM;
d3 = d3.rrNAL;
d3 = d3([1:41 52:61],:);
d3 = d3([1 2:11 32:41 22:31 12:21 42:51],:);

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

for i = 1:size(d3,1)-1 
        %%
        Index = find(strcmp([d0{:,1}], d3{i+1,1}));
        flag = 0;
        parts = split(d3{i+1,1}, '_');
        if parts{4} == '1'
            newnum = '1';
        elseif parts{4} == '10'
            newnum = '2';
        else
            newnum = mat2str(str2double(parts{4})+1);
        end
        fprintf(fileID,'%s',strcat([parts{2},'-',newnum]));
        fprintf(fileID, '& RiNNAL+ ');
        if Index ~= 0 
            if d0{Index,14}<=3600
                flag = 1;
                dd = d0(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{17},dd{5});
                fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf(fileID, '& %.1f &%.1f \\\\ \n',dd{14},dd{16});
            else
                flag = 1;
                dd = d0(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d ',dd{12},dd{13},dd{17},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf(fileID, '& %.1f &%.1f \\\\ \n',3600.0,dd{16});
            end
        end
        if flag == 0
            fprintf(fileID, '& - & - & - & - & - & - \\\\ \n');
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
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',dd{15});
            else
                flag = 1;
                dd = d3(Index,:);
                fprintf(fileID, '& %d, %d, %d & %d',dd{12},dd{13},dd{14},dd{5});
                fprintf(fileID, '& %.2e$^{\\dagger}$ ',max([dd{8},dd{9},dd{10}]));
                fprintf(fileID, '& %.7e ',dd{7});
                fprintf(fileID, '& %.1f & - \\\\ [3pt] \n\n',3600.0);
            end
        end
        if flag == 0 
            fprintf(fileID, '& - & - & - & - & - & - \\\\ [3pt] \n\n');
        end
end
fclose(fileID);
