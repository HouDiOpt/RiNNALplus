%%**********************************************************************
%% Overload operator 'ge/ >='
%%
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
%% Corresponding author: Kim-Chuan Toh
%%**********************************************************************

function constr_obj = ge(obj1, obj2)
    if isa(obj1, 'constraint')
        if obj1.status == -2
            if isa(obj2, 'double')
                [dim_m2, dim_n2] = size(obj1.Constant);
                if isscalar(obj2)
                    info.constr_string = strcat(obj1.constr_string, '>=', num2str(obj2));
                else
                    [dim_m1, dim_n1] = size(obj2);
                    if dim_m1 ~= dim_m2 || dim_n1 ~= dim_n2
                        error('Dimensions must agree.');
                    end
                    info.constr_string = strcat(obj1.constr_string, '>=', inputname(2));
                end
                info.symmetric_constr = obj1.symmetric_constr;
                if strcmp(obj1.constr_type, 'affine_constr')
                    info.constr_type = 'chain_constr';
                elseif strcmp(obj1.constr_type, 'upper_bound') || strcmp(obj1.constr_type, 'lower_bound')
                    info.constr_type = 'twodirect_bound';
                end
                info.operator_type = strcat(obj1.operator_type, '/>=');
                info.Operator_Matrix = obj1.Operator_Matrix;
                info.active_block = obj1.active_block;
                if strcmp(obj1.operator_type, '<=')
                    info.Constant.U = [];
                    info.Constant.L = max(obj1.Constant, obj2);
                else
                    info.Constant.U = obj1.Constant;
                    if isscalar(obj2)
                        info.Constant.L = obj2*ones(dim_m2,dim_n2);
                    else
                        info.Constant.L = obj2;
                    end
                end
                info.num_constr = obj1.num_constr;
                info.status = -3;
                constr_obj = constraint(info);
            else
                error('Invalid input: Please refer to chain type affine constraint in User Guide.');
            end
        else
            error('Invalid input: Please refer to chain type affine constraint in User Guide.');
        end
    else
        error('''%s'' is not a constraint object.', inputname(1));
    end
end