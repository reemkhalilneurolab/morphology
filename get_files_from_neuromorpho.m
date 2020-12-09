%==========================================================================
% @brief This fnction calls a python script that queries Nueromorpho.
% 
% Usage:
% python get_SWC.py will output some help information
% python get_SWC.py --region neocortex will download all SWC files from the region neocortex to current dir
% python get_SWC.py --region neocortex --neurons 10 will download the first 10 SWC files of the region neocortex to current dir
% python get_SWC.py --name cnic_001 will download the specified SWC file by name to current dir
% python get_SWC.py --index 1 will download the specified SWC file by index to current dir
% python get_SWC.py --archive Smith will download all SWC files of given archive name Smith to current dir
% python get_SWC.py --species 

% 
% commandStr = 'python C:\Users\rkhalil\Documents\GitHub\PARS-TMD\neuromorpho-master/get_SWC.py --region amygdala --neurons 50';

% commandStr = 'python C:\Users\rkhalil\Documents\GitHub\PARS-TMD\neuromorpho-master/get_SWC.py --archive Luebke --species mouse --region neocortex';
 commandStr = 'python C:\Users\rkhalil\Documents\GitHub\PARS-TMD\neuromorpho-master/get_SWC.py --archive Luebke';
 [status, commandOut] = system(commandStr);
