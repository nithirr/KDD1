
function [Best_pos1,Best_pos2,Best_pos3,Best_pos4,Best_pos5]=Optcomparison(nrun2,SearchAgents_no)

% SearchAgents_no=30; % Number of search agents

Function_name='F8'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

% Max_iteration=500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score1,Best_pos1,AGPSO1_cg_curve]= AGPSO1(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score2,Best_pos2,AGPSO2_cg_curve]= AGPSO2(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score3,Best_pos3,AGPSO3_cg_curve]= AGPSO3(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score4,Best_pos4,PSO_cg_curve]   = PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score5,Best_pos5,IPSO_cg_curve]= IPSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score6,Best_pos6,TACPSO_cg_curve]= TACPSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score7,Best_pos7,MPSO_cg_curve]= MPSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

% figure('Position',[500 500 660 290])
% 
% %Draw search space
% subplot(1,2,1);
% func_plot(Function_name);
% title('Parameter space')
% xlabel('x_1');
% ylabel('x_2');
% zlabel([Function_name,'( x_1 , x_2 )'])
% 
% %Draw convergence curves
% subplot(1,2,2);
% semilogy(AGPSO1_cg_curve,'Color','r')
% hold on
% semilogy(AGPSO2_cg_curve,'Color','b')
% semilogy(AGPSO3_cg_curve,'Color','k')
% semilogy(PSO_cg_curve,'Color','g')
% semilogy(MPSO_cg_curve,'Color','y')
% semilogy(TACPSO_cg_curve,'Color','c')
% semilogy(IPSO_cg_curve,'Color','m')
% 
% title('Objective space')
% xlabel('Iteration');
% ylabel('Best score obtained so far');
% 
% axis tight
% grid on
% box on
% legend('AGPSO1','AGPSO2','AGPSO3', 'PSO', 'MPSO', 'TACPSO', 'IPSO')

display(['The best solution obtained by AGPSO1 is : ', num2str(Best_pos1)]);
display(['The best optimal value obtained by AGPSO1 is : ', num2str(Best_score1)]);

display(['The best solution obtained by AGPSO2 is : ', num2str(Best_pos2)]);
display(['The best optimal value obtained by AGPSO2 is : ', num2str(Best_score2)]);

display(['The best solution obtained by AGPSO3 is : ', num2str(Best_pos3)]);
display(['The best optimal value obtained by AGPSO3 is : ', num2str(Best_score3)]);

display(['The best solution obtained by SPSO is : ', num2str(Best_pos4)]);
display(['The best optimal value obtained by SPSO is : ', num2str(Best_score4)]);

display(['The best solution obtained by MPSO is : ', num2str(Best_pos5)]);
display(['The best optimal value obtained by MPSO is : ', num2str(Best_score5)]);

display(['The best solution obtained by TACPSO is : ', num2str(Best_pos6)]);
display(['The best optimal value obtained by TACPSO is : ', num2str(Best_score6)]);

display(['The best solution obtained by IPSO is : ', num2str(Best_pos1)]);
display(['The best optimal value obtained by IPSO is : ', num2str(Best_score1)]);

        



