%%**********************************************************************
%% class def for var_sdp
%%
%% SDPNAL+: 
%% Copyright (c) 2017 by
%% Yancheng Yuan , Kim-Chuan Toh, Defeng Sun and Xinyuan Zhao
%%**********************************************************************

classdef var_sdp < matlab.mixin.Copyable
    properties
        blk
        blkorg
        L_Mat
        U_Mat
        block_no
        model
    end
    
    methods
        function obj = var_sdp(m,n)
            obj.blkorg{1} = 's';
            obj.blk{1} = 's';
            if (nargin==1)
                n = m; 
            end
            if (nargin <= 2)
                if (m ~= n)
                    error('Invalid Input: ''sdpvar'' must be square.');
                else
                    obj.blkorg{2}= m;
                    obj.blkorg{3}= n;
                    obj.blk{2}=m;                
                end
            end
            obj.L_Mat = -inf*ones(m,n);
            obj.U_Mat = -obj.L_Mat;
            obj.block_no = -1;
            obj.model = [];
        end
    end
end
        