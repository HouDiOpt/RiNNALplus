========================================================
This is the implementation of RiNNAL+ for solving 
SDP-RLT/DNN relaxations of MBQP problems:

Mixed-bianry Quadratic Programming (MBQP):
 min  x'Qx + 2c'x
 s.t. Ax = b                                (equality)
      Bx <= d                             (inequality)
      x_{i} in {0,1}, i in p                  (binary)
      x_i*x_j = 0, (i+1,j+1) in E    (complementarity)
      x >= 0                             (nonnegative)

The SDP-RLT relaxation of (MBQP):
   min  <Q,X> + 2<c,x>               
   s.t. 1. Ax = b                    
        2. AX = bx'                  
        3. diag(X)_B = x_B           
        4. (Bx <= d)               (usually redundant)
        5. BX <= dx'                 
        6. BXB'-Bxd'-dx'B'+dd' >= 0  
        7. Y_{ij} = 0, (i,j) in E    
        8. Y >= 0                    
        9. Y := [1,x';x,X] PSD       

Algorithm    : ALM + PG + Rie BB grad descent
Authors      : Di Hou, Tianyun Tang, Kim-Chuan Toh
Last edition : 2025/04/10
========================================================


===================== COMMENTS =========================
1. Put the path of RiNNALplus in setup_path.m

2. Put the path of SDPNAL+ in setup_path.m

3. Put the path of GUROBI in setup_path.m 
   (only needed for runs_tightness_BIQ.m)

4. MEX file is necessary for SDPNAL+ but not for RiNNAL+
========================================================
