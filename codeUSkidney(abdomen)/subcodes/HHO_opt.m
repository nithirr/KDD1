function [CNVG]=HHO_opt(Function_name,T,N)

% N=30; % Number of search agents
% 
% Function_name='F1'; % Name of the test function 
% 
% T=500; % Maximum number of iterations
% 
% % Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details1(Function_name);
[Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj);
display(['The best location of HHO is: ', num2str(Rabbit_Location)]);
display(['The best fitness of HHO is: ', num2str(Rabbit_Energy)]);

        



