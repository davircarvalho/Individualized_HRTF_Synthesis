clear all;clc
Obj = SOFAload('../Datasets/3D3A/Public-Data/Subject12/Subject12_HRIRs.sofa');
out_pos = [23, -11]; % desired source position (azimuth [deg], elevation [deg], radius [cm])
Obj_out = SOFA_barycentric_interp(Obj, out_pos, 'plot');