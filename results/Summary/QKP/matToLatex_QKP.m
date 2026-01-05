function matToLatex_QKP()
filename = 'TableTex_QKP.txt';
d0 = load("QKP-ALM.mat");
d0 = d0.rrALM;

for i = 1:size(d0,1)-1
    if iscell(d0{i+1,1})
        d0{i+1,1} = d0{i+1,1}{1};
    end
end

fileID = fopen(filename, 'w');
fclose(fileID);
fileID = fopen(filename, 'a');

count = 0;

for i = 1:size(d0,1)-1
        dd = d0(i+1,:);
        fprintf(fileID, '%d, %.1f',dd{4},dd{3});
        fprintf(fileID, '& %d, %d, %d & %d ',dd{14},dd{15},dd{18},dd{7});
        fprintf(fileID, '& %.2e ',max([dd{10},dd{11},dd{12}]));
        fprintf(fileID, '& %.7e ',dd{8});
        fprintf(fileID, '& %.1f &%.1f \\\\ [1pt]\n',dd{16},dd{17});
end
fclose(fileID);
